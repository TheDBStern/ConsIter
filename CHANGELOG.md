# Changelog

All notable changes to consiter are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.2.0] - 2026-11-17

First packaged release. `ConsIter.py` is replaced by the installable `consiter`
package, and the iteration loop gains new stopping tests, a second variant
caller and per-iteration diagnostics.

### Added
- Installable `src/consiter` package with a `consiter` command (tyro CLI) and
  `python -m consiter`.
- `environment.yml` that installs the tested versions of Bowtie2, Samtools,
  BCFtools, BedTools and GATK4 together with consiter.
- `--caller {gatk,bcftools}`. BCFtools calls from the pileup and recovers more
  of the genome from a distant starting reference; GATK remains the default.
- `--stop {both,rate,consensus}`. Besides the alignment-rate test, the loop can
  stop when the caller applies no new variant or returns a consensus seen in an
  earlier round (a cycle).
- `--call-min-cov` (default 2): the depth a variant needs to be applied to the
  next reference, separate from the masking depth `--min-cov`.
- `--no-dedup` to keep duplicate reads, for amplicon libraries.
- `<name>.iterations.tsv`, rewritten after every iteration: alignment rate,
  duplicate fraction, mean depth, breadth at `--min-cov`, masked bases,
  reference hash, whether the consensus changed, and the stop reason.
- Warnings when coverage is too thin to iterate from, when the duplicate
  fraction suggests an amplicon library, and when more than half of the output
  is masked.
- Gzipped reference input.
- Unit tests, and end-to-end tests on simulated reads (`pytest -m integration`).

### Changed
- `--max-iter` default raised from 5 to 10, since the new stopping tests
  usually end the loop sooner.
- Duplicates are marked with `gatk MarkDuplicates`; a separate Picard install
  is no longer needed.
- External tools run through `subprocess` with checked exit codes. Each step
  logs to `<outdir>/tmp/logs`, and `<outdir>/tmp` is kept when a step fails.
  A failing pipe reports every failed stage and its log.
- Indexes are built in the working directory instead of next to the input
  reference.
- consiter refuses to start if `<outdir>/tmp` already exists, since it is
  deleted on success.
- The SOR hard filter is dropped.

### Fixed
- Steps no longer run in the background under `dash` (`/bin/sh` on
  Debian/Ubuntu); previously Bowtie2 could start before its index existed.
- The final BAM and coverage mask come from reads mapped to the final
  consensus, not to the previous reference.
- FASTA header renaming: `>1` no longer matches `>10`, and sample-name headers
  are no longer duplicated across contigs.

### Removed
- `ConsIter.py`; use the `consiter` command.

## 0.1.0 - 2022-08-30

The original single-script `ConsIter.py` (not tagged): iterative Bowtie2
mapping, GATK HaplotypeCaller with hard filters, optional SNP-only updates
(`--no-indel`), and masking of low-coverage bases in the final consensus.

[Unreleased]: https://github.com/TheDBStern/ConsIter/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/TheDBStern/ConsIter/releases/tag/v0.2.0
