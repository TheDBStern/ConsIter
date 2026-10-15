"""The ConsIter iteration loop.

Each iteration maps the reads to the current reference. If the alignment rate
improved on the previous iteration (or this is the first one), variants are
called and applied to produce the next reference. Otherwise, or once
``max_iter`` consensus rounds have been made, the current reference is masked
for low coverage using the reads mapped to it in this iteration, and written out.
"""

from __future__ import annotations

import gzip
import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from consiter import tools
from consiter.fasta import read_fasta_names, rename_headers

if TYPE_CHECKING:
    from consiter.cli import Config

log = logging.getLogger("consiter")


@dataclass(frozen=True)
class Outputs:
    consensus: Path
    consensus_sample_name: Path
    bam: Path
    bai: Path
    iterations: int
    alignment_rates: list[float]
    masked_bases: int


def _copy_reference(src: Path, dest: Path) -> None:
    """Copy the reference into the work directory (decompressing gzip), so that
    index files are not written next to the user's input."""
    with open(src, "rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    opener = gzip.open if gzipped else open
    with opener(src, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)


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

    rates: list[float] = []
    for i in range(cfg.max_iter + 1):
        log.info("Iteration %d: mapping reads to %s", i, ref.name)
        tools.bowtie2_build(ref, logs / f"iter{i}.bowtie2-build.log")
        rate = tools.bowtie2_map(
            ref, cfg.left, cfg.right, cfg.name, cfg.threads,
            work / f"iter{i}.bam", work / f"iter{i}.btstats.txt", logs / f"iter{i}.bowtie2.log",
        )  # fmt: skip
        log.info("Iteration %d: alignment rate %.2f%%", i, rate)

        log.info("Iteration %d: removing duplicates", i)
        bam = work / f"iter{i}.rmdup.bam"
        tools.mark_duplicates(
            work / f"iter{i}.bam", bam, work / f"iter{i}.rmdup_metrics.txt",
            cfg.xmx, logs / f"iter{i}.markdup.log",
        )  # fmt: skip

        improved = not rates or rate > rates[-1]
        rates.append(rate)
        if not improved:
            log.info("No improvement over iteration %d (%.2f%%); stopping", i - 1, rates[-2])
            break
        if i == cfg.max_iter:
            log.info("Reached the maximum of %d consensus rounds; stopping", cfg.max_iter)
            break

        ref = _next_reference(cfg, i, ref, bam, names, work, logs)

    return _finalize(cfg, ref, bam, names, work, len(rates), rates)


def _next_reference(
    cfg: Config, i: int, ref: Path, bam: Path, names: list[str], work: Path, logs: Path
) -> Path:
    """Call and filter variants against ``ref`` and apply them to make the next reference."""
    log.info("Iteration %d: calling variants", i)
    tools.prepare_reference(ref, logs / f"iter{i}.prepare.log")
    vcf = work / f"iter{i}.vcf"
    tools.haplotype_caller(
        ref, bam, vcf, cfg.ploidy, cfg.xmx, logs / f"iter{i}.haplotypecaller.log"
    )
    filt = work / f"iter{i}.filt.vcf"
    tools.variant_filtration(ref, vcf, filt, cfg.min_cov, cfg.xmx, logs / f"iter{i}.filter.log")
    if cfg.no_indel:
        snps = work / f"iter{i}.snps.vcf"
        tools.select_snps(ref, filt, snps, cfg.xmx, logs / f"iter{i}.select.log")
        filt = snps

    log.info("Iteration %d: generating updated reference", i)
    raw = work / f"consensus.iter{i}.tmpnames.fa"
    tools.alternate_reference(ref, filt, raw, cfg.xmx, logs / f"iter{i}.altref.log")
    new_ref = work / f"ref.iter{i + 1}.fa"
    rename_headers(raw, new_ref, names)
    return new_ref


def _finalize(
    cfg: Config,
    ref: Path,
    bam: Path,
    names: list[str],
    work: Path,
    iterations: int,
    rates: list[float],
) -> Outputs:
    log.info("Masking bases with coverage < %d", cfg.min_cov)
    masked_fa = work / "final.masked.fa"
    masked = tools.mask_low_coverage(
        ref,
        bam,
        cfg.min_cov,
        masked_fa,
        work / "final.masked_sites.bed",
        work / "logs" / "mask.log",
    )
    log.info("Masked %d bases", masked)

    out = Outputs(
        consensus=cfg.outdir / f"{cfg.name}.consensus.fa",
        consensus_sample_name=cfg.outdir / f"{cfg.name}.consensus.sample_name.fa",
        bam=cfg.outdir / f"{cfg.name}.bt2.rmdup.bam",
        bai=cfg.outdir / f"{cfg.name}.bt2.rmdup.bai",
        iterations=iterations,
        alignment_rates=rates,
        masked_bases=masked,
    )
    rename_headers(masked_fa, out.consensus, names)
    rename_headers(masked_fa, out.consensus_sample_name, names, prefix=cfg.name)
    shutil.move(bam, out.bam)
    shutil.move(bam.with_suffix(".bai"), out.bai)

    if not cfg.keep_tmp:
        shutil.rmtree(work)

    log.info("Updated reference (bases < %dx masked): %s", cfg.min_cov, out.consensus)
    log.info("Same, with sample name in headers: %s", out.consensus_sample_name)
    log.info("Duplicate-removed BAM aligned to the updated reference: %s", out.bam)
    return out
