"""Thin wrappers around the external programs used by the pipeline.

Every wrapper runs its program with an argument list (no shell), writes the
command and its output to a log file, and raises ``ToolError`` on failure.
"""

from __future__ import annotations

import re
import shlex
import shutil
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from consiter.fasta import low_coverage_intervals

REQUIRED_PROGRAMS = ("bowtie2", "bowtie2-build", "samtools", "gatk", "bedtools")
# Only needed for --caller bcftools, so it is checked for separately.
BCFTOOLS_PROGRAMS = (*REQUIRED_PROGRAMS, "bcftools")

# GATK's suggested hard filters for germline SNPs, minus SOR: SOR also penalises
# overall strand imbalance, which is expected at contig ends and amplicon edges.
# FS still catches allele-specific strand bias. Depth is filtered separately.
VARIANT_FILTERS = {
    "QD_filter": "QD < 2.0",
    "FS_filter": "FS > 60.0",
    "MQ_filter": "MQ < 20.0",
    "ReadPosRankSum_filter": "ReadPosRankSum < -8.0",
}


class ToolError(RuntimeError):
    """An external program exited with a non-zero status."""


def missing_programs(programs: Sequence[str] = REQUIRED_PROGRAMS) -> list[str]:
    """Return the programs that are not on ``PATH``."""
    return [p for p in programs if shutil.which(p) is None]


def _tail(path: Path, n: int = 20) -> str:
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-n:])


def _fail(cmd: Sequence[str], returncode: int, log: Path) -> ToolError:
    return ToolError(
        f"{Path(cmd[0]).name} exited with status {returncode}; full log: {log}\n{_tail(log)}"
    )


def run(cmd: Sequence[str], log: Path, stdout: Path | None = None) -> None:
    """Run ``cmd``, logging stderr (and stdout unless redirected to ``stdout``) to ``log``."""
    cmd = [str(c) for c in cmd]
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w", encoding="utf-8") as lf:
        lf.write(f"$ {shlex.join(cmd)}\n")
        lf.flush()
        if stdout is None:
            proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
        else:
            with open(stdout, "wb") as out:
                proc = subprocess.run(cmd, stdout=out, stderr=lf)
    if proc.returncode != 0:
        raise _fail(cmd, proc.returncode, log)


def _gatk(tool: str, xmx: int, *args: object) -> list[str]:
    return ["gatk", "--java-options", f"-Xmx{xmx}g", tool, *map(str, args)]


def parse_alignment_rate(bowtie2_stderr: str) -> float:
    """Extract the overall alignment rate (percent) from bowtie2's summary."""
    m = re.search(r"([\d.]+)% overall alignment rate", bowtie2_stderr)
    if m is None:
        raise ValueError("bowtie2 output has no 'overall alignment rate' line")
    return float(m.group(1))


def bowtie2_build(ref: Path, log: Path) -> None:
    """Build a bowtie2 index next to ``ref`` (prefix = ``ref``)."""
    run(["bowtie2-build", "--quiet", ref, ref], log)


def bowtie2_map(
    index: Path,
    left: Path,
    right: Path,
    sample: str,
    threads: int,
    out_bam: Path,
    stats: Path,
    log: Path,
    *,
    preset: Sequence[str] = ("--sensitive-local",),
    min_mapq: int = 20,
    proper_pair: bool = True,
) -> float:
    """Map paired reads, keep properly paired MAPQ>=``min_mapq`` reads, sort to ``out_bam``.

    ``preset`` holds the bowtie2 alignment options. bowtie2's summary is written to
    ``stats`` and the overall alignment rate is returned.
    """
    bt2 = [
        "bowtie2", "-p", threads, "--no-unal", *preset,
        "--rg-id", sample, "--rg", f"SM:{sample}", "--rg", "LB:1", "--rg", "PU:1",
        "--rg", "PL:illumina",
        "-x", index, "-1", left, "-2", right,
    ]  # fmt: skip
    view = ["samtools", "view", "-u", *(["-f", "0x2"] if proper_pair else []),
            "-q", str(min_mapq), "-"]  # fmt: skip
    sort = ["samtools", "sort", "-@", threads, "-o", out_bam, "-"]
    cmds = [[str(c) for c in cmd] for cmd in (bt2, view, sort)]

    log.parent.mkdir(parents=True, exist_ok=True)
    with open(stats, "w", encoding="utf-8") as st, open(log, "w", encoding="utf-8") as lf:
        lf.write("$ " + " | ".join(shlex.join(c) for c in cmds) + "\n")
        lf.flush()
        p1 = subprocess.Popen(cmds[0], stdout=subprocess.PIPE, stderr=st)
        p2 = subprocess.Popen(cmds[1], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=lf)
        p1.stdout.close()  # so bowtie2 gets SIGPIPE if samtools dies
        p3 = subprocess.Popen(cmds[2], stdin=p2.stdout, stderr=lf)
        p2.stdout.close()
        p3.wait()
        p2.wait()
        p1.wait()
    # A failure downstream makes upstream stages die of SIGPIPE, so report every
    # failed stage and every relevant log rather than just the first one.
    failed = [
        (cmd, proc, stats if proc is p1 else log)
        for cmd, proc in zip(cmds, (p1, p2, p3), strict=True)
        if proc.returncode != 0
    ]
    if failed:
        msg = "\n".join(
            f"{' '.join(cmd[:2])} exited with status {p.returncode}" for cmd, p, _ in failed
        )
        for path in dict.fromkeys(path for _, _, path in failed):
            msg += f"\n--- {path} ---\n{_tail(path)}"
        raise ToolError(msg)
    return parse_alignment_rate(stats.read_text(encoding="utf-8"))


def index_bam(bam: Path, log: Path) -> None:
    """Index ``bam`` as ``<stem>.bai``, the name MarkDuplicates' ``--CREATE_INDEX`` writes."""
    run(["samtools", "index", "-o", bam.with_suffix(".bai"), bam], log)


def duplicate_fraction(metrics: Path) -> float | None:
    """``PERCENT_DUPLICATION`` from a MarkDuplicates metrics file.

    ``None`` when the file carries no usable value, as it does for an empty BAM.
    """
    with open(metrics, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("LIBRARY\t"):
                continue
            keys = line.rstrip("\n").split("\t")
            values = next(fh, "").rstrip("\n").split("\t")
            try:
                return float(dict(zip(keys, values, strict=False))["PERCENT_DUPLICATION"])
            except (KeyError, ValueError):
                return None
    return None


def mark_duplicates(in_bam: Path, out_bam: Path, metrics: Path, xmx: int, log: Path) -> None:
    """Remove duplicate reads and index the result (``<out_bam stem>.bai``)."""
    run(
        _gatk(
            "MarkDuplicates",
            xmx,
            "-I",
            in_bam,
            "-O",
            out_bam,
            "-M",
            metrics,
            "--REMOVE_DUPLICATES",
            "true",
            "--CREATE_INDEX",
            "true",
            "--VALIDATION_STRINGENCY",
            "LENIENT",
            "--USE_JDK_DEFLATER",
            "true",
            "--USE_JDK_INFLATER",
            "true",
        ),  # fmt: skip
        log,
    )


def prepare_reference(ref: Path, log: Path) -> None:
    """Create the ``.fai`` and ``.dict`` files GATK needs for ``ref``."""
    run(["samtools", "faidx", ref], log.with_suffix(".faidx.log"))
    run(["samtools", "dict", "-o", ref.with_suffix(".dict"), ref], log.with_suffix(".dict.log"))


def haplotype_caller(ref: Path, bam: Path, vcf: Path, ploidy: int, xmx: int, log: Path) -> None:
    run(
        _gatk(
            "HaplotypeCaller",
            xmx,
            "--use-jdk-deflater",
            "--use-jdk-inflater",
            "-R",
            ref,
            "-I",
            bam,
            "-O",
            vcf,
            "--min-base-quality-score",
            20,
            "-ploidy",
            ploidy,
        ),  # fmt: skip
        log,
    )


def variant_filtration(
    ref: Path, vcf: Path, out_vcf: Path, min_depth: int, xmx: int, log: Path
) -> None:
    """Flag variants failing ``VARIANT_FILTERS`` or with depth below ``min_depth``."""
    filters = []
    for name, expr in {"DP_filter": f"DP < {min_depth}", **VARIANT_FILTERS}.items():
        filters += ["--filter-name", name, "--filter-expression", expr]
    run(_gatk("VariantFiltration", xmx, "-R", ref, "-V", vcf, "-O", out_vcf, *filters), log)


def select_snps(ref: Path, vcf: Path, out_vcf: Path, xmx: int, log: Path) -> None:
    run(
        _gatk(
            "SelectVariants",
            xmx,
            "--use-jdk-deflater",
            "--use-jdk-inflater",
            "-R",
            ref,
            "-V",
            vcf,
            "--select-type-to-include",
            "SNP",
            "-O",
            out_vcf,
        ),  # fmt: skip
        log,
    )


def alternate_reference(ref: Path, vcf: Path, out_fasta: Path, xmx: int, log: Path) -> None:
    """Apply the (unfiltered) variants in ``vcf`` to ``ref``."""
    run(_gatk("FastaAlternateReferenceMaker", xmx, "-R", ref, "-V", vcf, "-O", out_fasta), log)


def _bcftools_ploidy(ploidy: int) -> list[str]:
    """``bcftools call``'s ploidy options: it takes presets, not an arbitrary number."""
    if ploidy == 1:
        return ["--ploidy", "1"]
    if ploidy == 2:
        return []  # bcftools calls diploid by default
    raise ValueError(f"the bcftools caller supports --ploidy 1 or 2 (got {ploidy})")


def bcftools_consensus(
    ref: Path,
    bam: Path,
    out_fasta: Path,
    min_depth: int,
    ploidy: int,
    log: Path,
    *,
    snps_only: bool = False,
) -> None:
    """Call variants from the pileup and apply the passing ones to ``ref``.

    A pileup alternative to ``haplotype_caller`` + ``variant_filtration`` +
    ``alternate_reference``; ``min_depth`` is the depth a variant needs to be applied,
    as GATK's ``DP_filter`` is there. Intermediate files are written beside
    ``out_fasta``, so ``--keep-tmp`` exposes the calls.
    """
    ploidy_args = _bcftools_ploidy(ploidy)  # before any work, so a bad ploidy fails fast
    run(["samtools", "faidx", ref], log.with_suffix(".faidx.log"))
    pileup = out_fasta.with_suffix(".pileup.bcf")
    run(
        ["bcftools", "mpileup", "-f", ref, "-a", "AD,DP", "-q", 20, "-Q", 20,
         "-d", 1000000, "-Ob", "-o", pileup, bam],
        log.with_suffix(".mpileup.log"),
    )  # fmt: skip
    calls = out_fasta.with_suffix(".calls.bcf")
    run(["bcftools", "call", *ploidy_args, "-mv", "-Ob", "-o", calls, pileup],
        log.with_suffix(".call.log"))  # fmt: skip
    normed = out_fasta.with_suffix(".norm.bcf")
    run(["bcftools", "norm", "-f", ref, "-Ob", "-o", normed, calls], log.with_suffix(".norm.log"))
    passing = out_fasta.with_suffix(".pass.vcf.gz")
    run(
        ["bcftools", "view", "-i", f"QUAL>=20 && INFO/DP>={min_depth}",
         *(["-V", "indels"] if snps_only else []), "-Oz", "-o", passing, normed],
        log.with_suffix(".view.log"),
    )  # fmt: skip
    run(["bcftools", "index", passing], log.with_suffix(".index.log"))
    run(["bcftools", "consensus", "-f", ref, "-o", out_fasta, passing],
        log.with_suffix(".consensus.log"))  # fmt: skip


@dataclass(frozen=True)
class Coverage:
    """Depth of one BAM over the reference it was mapped to."""

    low_bases: int
    """Bases covered by fewer than ``min_cov`` reads: what masking would replace with N."""
    total_bases: int
    mean_depth: float

    @property
    def breadth(self) -> float:
        """Fraction of the reference covered at least ``min_cov`` deep."""
        return 1 - self.low_bases / self.total_bases if self.total_bases else 0.0

    @property
    def masked_fraction(self) -> float:
        return self.low_bases / self.total_bases if self.total_bases else 0.0


def coverage_summary(bam: Path, min_cov: int, out_bed: Path, log: Path) -> Coverage:
    """Summarise depth, writing the intervals below ``min_cov`` to ``out_bed``.

    The BED is what ``mask_fasta`` masks, so one pass of ``bedtools genomecov`` serves
    both the statistics reported for an iteration and the masking of the final one.
    """
    bedgraph = out_bed.with_suffix(".bedgraph")
    run(
        ["bedtools", "genomecov", "-ibam", bam, "-bga"], log.with_suffix(".genomecov.log"), bedgraph
    )
    low = 0
    with open(bedgraph, encoding="utf-8") as fin, open(out_bed, "w", encoding="utf-8") as fout:
        for line in low_coverage_intervals(fin, min_cov):
            _, start, end = line.split("\t")
            low += int(end) - int(start)
            fout.write(line)
    total, depth_sum = 0, 0.0
    with open(bedgraph, encoding="utf-8") as fin:
        for line in fin:
            _, start, end, depth = line.split("\t")[:4]
            span = int(end) - int(start)
            total += span
            depth_sum += float(depth) * span
    return Coverage(low, total, depth_sum / total if total else 0.0)


def mask_fasta(ref: Path, bed: Path, out_fasta: Path, log: Path) -> None:
    """Write ``ref`` to ``out_fasta`` with the bases in ``bed`` replaced by N."""
    run(
        ["bedtools", "maskfasta", "-fi", ref, "-bed", bed, "-fo", out_fasta],
        log.with_suffix(".maskfasta.log"),
    )
