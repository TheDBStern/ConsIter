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
from pathlib import Path

from consiter.fasta import low_coverage_intervals

REQUIRED_PROGRAMS = ("bowtie2", "bowtie2-build", "samtools", "gatk", "bedtools")

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
) -> float:
    """Map paired reads, keep properly paired MAPQ>=20 reads, sort to ``out_bam``.

    bowtie2's summary is written to ``stats`` and the overall alignment rate is returned.
    """
    bt2 = [
        "bowtie2", "-p", threads, "--no-unal", "--sensitive-local",
        "--rg-id", sample, "--rg", f"SM:{sample}", "--rg", "LB:1", "--rg", "PU:1",
        "--rg", "PL:illumina",
        "-x", index, "-1", left, "-2", right,
    ]  # fmt: skip
    view = ["samtools", "view", "-u", "-f", "0x2", "-q", "20", "-"]
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
    for cmd, proc in zip(cmds, (p1, p2, p3), strict=True):
        if proc.returncode != 0:
            raise _fail(cmd, proc.returncode, stats if proc is p1 else log)
    return parse_alignment_rate(stats.read_text(encoding="utf-8"))


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


def mask_low_coverage(
    ref: Path, bam: Path, min_cov: int, out_fasta: Path, out_bed: Path, log: Path
) -> int:
    """Replace bases covered by fewer than ``min_cov`` reads with N.

    Returns the number of masked bases.
    """
    bedgraph = out_bed.with_suffix(".bedgraph")
    run(
        ["bedtools", "genomecov", "-ibam", bam, "-bga"], log.with_suffix(".genomecov.log"), bedgraph
    )
    masked = 0
    with open(bedgraph, encoding="utf-8") as fin, open(out_bed, "w", encoding="utf-8") as fout:
        for line in low_coverage_intervals(fin, min_cov):
            _, start, end = line.split("\t")
            masked += int(end) - int(start)
            fout.write(line)
    run(
        ["bedtools", "maskfasta", "-fi", ref, "-bed", out_bed, "-fo", out_fasta],
        log.with_suffix(".maskfasta.log"),
    )
    return masked
