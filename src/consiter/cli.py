"""Command-line interface."""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Literal

import tyro

from consiter import pipeline, tools

Arg = tyro.conf.arg


@dataclass(frozen=True)
class Config:
    """ConsIter produces an updated reference genome that more closely matches the
    sequenced viral population by iteratively mapping paired-end reads to the
    reference and generating a new consensus sequence."""

    ref: Annotated[Path, Arg(aliases=("-r",))]
    """Reference genome FASTA (may be gzipped)."""
    name: Annotated[str, Arg(aliases=("-n",))]
    """Sample name, used for read groups and output file names."""
    left: Annotated[Path, Arg(aliases=("-1",))]
    """Left (R1) reads in FASTQ format."""
    right: Annotated[Path, Arg(aliases=("-2",))]
    """Right (R2) reads in FASTQ format."""
    outdir: Annotated[Path, Arg(aliases=("-o",))] = Path(".")
    """Output directory."""
    max_iter: Annotated[int, Arg(aliases=("-i",))] = 10
    """Maximum number of consensus rounds; the loop usually stops earlier, at a fixed point."""
    threads: Annotated[int, Arg(aliases=("-t",))] = 16
    """Number of threads for bowtie2 and samtools."""
    ploidy: Annotated[int, Arg(aliases=("-p",))] = 1
    """Sample ploidy for the variant caller (the bcftools caller supports 1 or 2)."""
    min_cov: Annotated[int, Arg(aliases=("-c",))] = 20
    """Minimum depth for a base not to be masked in the final consensus.
    Coverage breadth is reported at this depth too."""
    call_min_cov: int = 2
    """Minimum depth for a variant to be applied to the next reference.
    Below --min-cov on purpose: the output is still masked at --min-cov, so a low gate lets
    the loop repair thin regions without reporting what it called there."""
    stop: Literal["both", "rate", "consensus"] = "both"
    """Which convergence test ends the loop, besides the --max-iter cap.
    rate: the alignment rate stopped improving.
    consensus: the caller applied nothing (a fixed point) or returned a consensus
    seen before (a cycle). Both are reported either way; this picks what breaks."""
    caller: Literal["gatk", "bcftools"] = "gatk"
    """Variant caller for the consensus rounds.
    bcftools recovers much more of the genome from a distant starting reference;
    GATK is slightly better once the reference is already close."""
    no_dedup: bool = False
    """Keep duplicate reads, as an amplicon library needs.
    Every read of an amplicon starts at the same primer boundary, so it is
    indistinguishable from a PCR duplicate."""
    xmx: int = 50
    """Maximum Java heap size for GATK, in GB."""
    keep_tmp: bool = False
    """Keep the temporary working directory (<outdir>/tmp)."""
    no_indel: bool = False
    """Only apply SNPs, so the updated reference keeps the original coordinates."""


def validate(cfg: Config) -> list[str]:
    """Return a list of problems with ``cfg`` (empty if it is usable)."""
    errors = []
    for label, path in (
        ("reference", cfg.ref),
        ("left reads", cfg.left),
        ("right reads", cfg.right),
    ):
        if not path.is_file():
            errors.append(f"Cannot find {label} file: {path}")
    for label, value in (("--max-iter", cfg.max_iter), ("--threads", cfg.threads),
                         ("--ploidy", cfg.ploidy), ("--xmx", cfg.xmx)):  # fmt: skip
        if value < 1:
            errors.append(f"{label} must be at least 1 (got {value})")
    if (cfg.outdir / "tmp").exists():
        errors.append(
            f"{cfg.outdir / 'tmp'} already exists and would be deleted; "
            "remove it or choose another --outdir"
        )
    for label, value in (("--min-cov", cfg.min_cov), ("--call-min-cov", cfg.call_min_cov)):
        if value < 0:
            errors.append(f"{label} must not be negative (got {value})")
    if cfg.caller == "bcftools" and cfg.ploidy > 2:
        errors.append(f"--caller bcftools supports --ploidy 1 or 2 (got {cfg.ploidy})")
    if any(c in cfg.name for c in "/\\| \t"):
        errors.append(f"--name must not contain '/', '\\', '|' or whitespace (got {cfg.name!r})")
    return errors


def main(args: list[str] | None = None) -> None:
    cfg = tyro.cli(Config, args=args, prog="consiter", config=(tyro.conf.FlagCreatePairsOff,))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s", datefmt="%H:%M:%S")

    errors = validate(cfg)
    required = tools.BCFTOOLS_PROGRAMS if cfg.caller == "bcftools" else tools.REQUIRED_PROGRAMS
    missing = tools.missing_programs(required)
    if missing:
        errors.append(
            "Required programs not found on PATH: "
            + ", ".join(missing)
            + " (activate the conda environment from environment.yml)"
        )
    if errors:
        sys.exit("Error: " + "\n       ".join(errors))

    try:
        pipeline.run(cfg)
    except (tools.ToolError, ValueError) as e:
        sys.exit(f"Error: {e}\nIntermediate files are kept in {cfg.outdir / 'tmp'}")
