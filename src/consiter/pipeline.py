"""The ConsIter iteration loop.

Each iteration maps the reads to the current reference, removes duplicates and
records mapping and coverage statistics. If the alignment rate improved on the
previous iteration (or this is the first one), variants are called and applied to
produce the next reference.

Two convergence tests can end the loop, and ``--stop`` chooses which of them break it:
``rate``, when the overall alignment rate does not improve on the previous iteration,
and ``consensus``, when the caller applies nothing (a fixed point) or returns a
consensus seen before (a cycle). Either way ``max_iter`` consensus rounds is the cap.
The current reference is then masked for low coverage using the reads mapped to it in
this iteration, and written out.

Two depths are used, and they are not the same: a variant needs ``call_min_cov``
reads to be applied to the next reference, while the final consensus is masked at
``min_cov``. Keeping the first low lets the loop repair thin regions, because a
region that is never repaired never recruits the reads that would deepen it; keeping
the second where it is means a thin call only reaches the output through the
reference it shaped.
"""

from __future__ import annotations

import csv
import gzip
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from consiter import tools
from consiter.fasta import read_fasta_names, rename_headers, sequence_hash

if TYPE_CHECKING:
    from consiter.cli import Config

log = logging.getLogger("consiter")

ITERATION_COLUMNS = (
    "iteration", "alignment_rate", "duplicate_fraction", "mean_depth", "breadth",
    "masked_bases", "ref_hash", "changed", "stop_reason",
)  # fmt: skip

LOW_BREADTH = 0.05
"""Breadth at iteration 0 below which iteration is unlikely to achieve anything."""
HIGH_DUPLICATE_FRACTION = 0.8
"""Duplicate fraction above which an amplicon library is the likely explanation."""
HIGH_MASKED_FRACTION = 0.5
"""Masked fraction above which the consensus should not be reported as a success."""


@dataclass(frozen=True)
class Outputs:
    consensus: Path
    consensus_sample_name: Path
    bam: Path
    bai: Path
    iterations_tsv: Path
    iterations: int
    alignment_rates: list[float]
    masked_bases: int
    masked_fraction: float
    breadth: float
    stop_reason: str


def _copy_reference(src: Path, dest: Path) -> None:
    """Copy the reference into the work directory (decompressing gzip), so that
    index files are not written next to the user's input."""
    with open(src, "rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    opener = gzip.open if gzipped else open
    with opener(src, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def _write_iterations(path: Path, rows: list[dict]) -> None:
    """(Re)write the per-iteration statistics, so they survive a later failure."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ITERATION_COLUMNS, delimiter="\t", restval="")
        writer.writeheader()
        writer.writerows(rows)


def _deduplicate(cfg: Config, i: int, work: Path, logs: Path) -> tuple[Path, float | None]:
    """Remove duplicate reads (unless ``--no-dedup``) and index the result."""
    mapped = work / f"iter{i}.bam"
    bam = work / f"iter{i}.rmdup.bam"
    if cfg.no_dedup:
        log.info("Iteration %d: keeping duplicates (--no-dedup)", i)
        shutil.copy(mapped, bam)
        tools.index_bam(bam, logs / f"iter{i}.index.log")
        return bam, None

    log.info("Iteration %d: removing duplicates", i)
    metrics = work / f"iter{i}.rmdup_metrics.txt"
    tools.mark_duplicates(mapped, bam, metrics, cfg.xmx, logs / f"iter{i}.markdup.log")
    duplicates = tools.duplicate_fraction(metrics)
    if duplicates is not None and duplicates > HIGH_DUPLICATE_FRACTION:
        log.warning(
            "Iteration %d: %.1f%% of reads were removed as duplicates. An amplicon "
            "library looks like this, because every read of an amplicon starts at the "
            "same primer boundary; if this is one, rerun with --no-dedup.",
            i, 100 * duplicates,
        )  # fmt: skip
    return bam, duplicates


def run(cfg: Config) -> Outputs:
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    work = cfg.outdir / "tmp"
    try:
        work.mkdir()  # never reuse an existing directory: it is deleted at the end
    except FileExistsError:
        raise ValueError(f"{work} already exists; remove it or choose another --outdir") from None
    logs = work / "logs"
    logs.mkdir()

    ref = work / "ref.iter0.fa"
    _copy_reference(cfg.ref, ref)
    names = read_fasta_names(ref)
    if not names:
        raise ValueError(f"No sequences found in {cfg.ref}")

    tsv = cfg.outdir / f"{cfg.name}.iterations.tsv"
    rates: list[float] = []
    rows: list[dict] = []
    seen: dict[str, int] = {}
    stop_reason = "max_iter"
    for i in range(cfg.max_iter + 1):
        log.info("Iteration %d: mapping reads to %s", i, ref.name)
        tools.bowtie2_build(ref, logs / f"iter{i}.bowtie2-build.log")
        rate = tools.bowtie2_map(
            ref, cfg.left, cfg.right, cfg.name, cfg.threads,
            work / f"iter{i}.bam", work / f"iter{i}.btstats.txt", logs / f"iter{i}.bowtie2.log",
        )  # fmt: skip
        log.info("Iteration %d: alignment rate %.2f%%", i, rate)

        bam, duplicates = _deduplicate(cfg, i, work, logs)
        bed = work / f"iter{i}.lowcov.bed"
        cov = tools.coverage_summary(bam, cfg.min_cov, bed, logs / f"iter{i}.coverage.log")
        log.info(
            "Iteration %d: %.1f%% of the reference covered at >=%dx, mean depth %.1fx",
            i, 100 * cov.breadth, cfg.min_cov, cov.mean_depth,
        )  # fmt: skip
        ref_hash = sequence_hash(ref)
        seen.setdefault(ref_hash, i)
        rows.append({
            "iteration": i,
            "alignment_rate": f"{rate:.2f}",
            "duplicate_fraction": "" if duplicates is None else f"{duplicates:.4f}",
            "mean_depth": f"{cov.mean_depth:.2f}",
            "breadth": f"{cov.breadth:.4f}",
            "masked_bases": cov.low_bases,
            "ref_hash": ref_hash,
        })  # fmt: skip
        _write_iterations(tsv, rows)
        if i == 0 and cov.breadth < LOW_BREADTH:
            log.warning(
                "Only %.1f%% of the reference reaches %dx before any iteration. Either "
                "the reference is too distant for the reads to map, or the library is "
                "too shallow; runs starting this low rarely improve, whatever the "
                "consensus caller.",
                100 * cov.breadth, cfg.min_cov,
            )  # fmt: skip

        improved = not rates or rate > rates[-1]
        rates.append(rate)
        if not improved and cfg.stop in ("both", "rate"):
            stop_reason = "rate"
            log.info("No improvement over iteration %d (%.2f%%); stopping", i - 1, rates[-2])
            break
        if i == cfg.max_iter:
            log.info("Reached the maximum of %d consensus rounds; stopping", cfg.max_iter)
            break

        new_ref = _next_reference(cfg, i, ref, bam, names, work, logs)
        new_hash = sequence_hash(new_ref)
        rows[-1]["changed"] = int(new_hash != ref_hash)
        if cfg.stop in ("both", "consensus"):
            if new_hash == ref_hash:
                stop_reason = "fixed_point"
                log.info("Iteration %d changed nothing; stopping at the fixed point", i)
                break
            if new_hash in seen:
                stop_reason = "cycle"
                log.info(
                    "Iteration %d returned the consensus of iteration %d; stopping",
                    i, seen[new_hash],
                )  # fmt: skip
                break
        ref = new_ref

    rows[-1]["stop_reason"] = stop_reason
    _write_iterations(tsv, rows)
    return _finalize(cfg, ref, bam, bed, cov, names, work, tsv, rates, stop_reason)


def _next_reference(
    cfg: Config, i: int, ref: Path, bam: Path, names: list[str], work: Path, logs: Path
) -> Path:
    """Call variants against ``ref`` and apply the passing ones to make the next reference."""
    log.info("Iteration %d: calling variants with %s", i, cfg.caller)
    raw = work / f"consensus.iter{i}.tmpnames.fa"
    if cfg.caller == "bcftools":
        tools.bcftools_consensus(
            ref, bam, raw, cfg.call_min_cov, cfg.ploidy, logs / f"iter{i}.bcftools.log",
            snps_only=cfg.no_indel,
        )  # fmt: skip
    else:
        _gatk_reference(cfg, i, ref, bam, raw, work, logs)
    log.info("Iteration %d: generating updated reference", i)
    new_ref = work / f"ref.iter{i + 1}.fa"
    rename_headers(raw, new_ref, names)
    return new_ref


def _gatk_reference(
    cfg: Config, i: int, ref: Path, bam: Path, raw: Path, work: Path, logs: Path
) -> None:
    """HaplotypeCaller, hard filters and FastaAlternateReferenceMaker."""
    tools.prepare_reference(ref, logs / f"iter{i}.prepare.log")
    vcf = work / f"iter{i}.vcf"
    tools.haplotype_caller(
        ref, bam, vcf, cfg.ploidy, cfg.xmx, logs / f"iter{i}.haplotypecaller.log"
    )
    filt = work / f"iter{i}.filt.vcf"
    tools.variant_filtration(
        ref, vcf, filt, cfg.call_min_cov, cfg.xmx, logs / f"iter{i}.filter.log"
    )
    if cfg.no_indel:
        snps = work / f"iter{i}.snps.vcf"
        tools.select_snps(ref, filt, snps, cfg.xmx, logs / f"iter{i}.select.log")
        filt = snps
    tools.alternate_reference(ref, filt, raw, cfg.xmx, logs / f"iter{i}.altref.log")


def _finalize(
    cfg: Config,
    ref: Path,
    bam: Path,
    bed: Path,
    cov: tools.Coverage,
    names: list[str],
    work: Path,
    tsv: Path,
    rates: list[float],
    stop_reason: str,
) -> Outputs:
    log.info("Masking %d bases with coverage < %d", cov.low_bases, cfg.min_cov)
    masked_fa = work / "final.masked.fa"
    tools.mask_fasta(ref, bed, masked_fa, work / "logs" / "mask.log")

    out = Outputs(
        consensus=cfg.outdir / f"{cfg.name}.consensus.fa",
        consensus_sample_name=cfg.outdir / f"{cfg.name}.consensus.sample_name.fa",
        bam=cfg.outdir / f"{cfg.name}.bt2.rmdup.bam",
        bai=cfg.outdir / f"{cfg.name}.bt2.rmdup.bai",
        iterations_tsv=tsv,
        iterations=len(rates),
        alignment_rates=rates,
        masked_bases=cov.low_bases,
        masked_fraction=cov.masked_fraction,
        breadth=cov.breadth,
        stop_reason=stop_reason,
    )
    rename_headers(masked_fa, out.consensus, names)
    rename_headers(masked_fa, out.consensus_sample_name, names, prefix=cfg.name)
    shutil.move(bam, out.bam)
    shutil.move(bam.with_suffix(".bai"), out.bai)

    if not cfg.keep_tmp:
        shutil.rmtree(work)

    if out.masked_fraction > HIGH_MASKED_FRACTION:
        log.warning(
            "%.1f%% of the consensus is masked. Neither the alignment rate nor a stable "
            "consensus detects this, so do not read either as success: check the breadth "
            "column of %s, and see --no-dedup for an amplicon library.",
            100 * out.masked_fraction, out.iterations_tsv.name,
        )  # fmt: skip
    log.info("Updated reference (bases < %dx masked): %s", cfg.min_cov, out.consensus)
    log.info("Same, with sample name in headers: %s", out.consensus_sample_name)
    log.info("Duplicate-removed BAM aligned to the updated reference: %s", out.bam)
    log.info("Per-iteration statistics: %s", out.iterations_tsv)
    return out
