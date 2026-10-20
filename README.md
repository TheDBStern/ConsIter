# ConsIter
Evolving viral populations often diverge substantially from available reference genomes, potentially reducing read mapping accuracy and biasing variant calls.
ConsIter produces an updated reference genome that more closely matches the viral population by iteratively mapping reads to the reference and generating a new consensus sequence.

In each iteration, paired-end reads are mapped to the current reference with Bowtie2, duplicates are removed, variants are called (with GATK4 HaplotypeCaller and hard filters on QD, FS, MQ, ReadPosRankSum and depth, or with `bcftools` from the pileup), and the passing variants are applied to make the next reference.
This continues until a convergence test fires or the maximum number of consensus rounds is reached. There are two tests, and `--stop` chooses which of them end the loop: the **rate** test stops when the overall alignment rate no longer improves on the previous iteration, and the **consensus** test stops when the caller applies nothing (a fixed point) or returns a consensus seen in an earlier round (a cycle). Both are reported in `<name>.iterations.tsv` whichever is selected.
The latest consensus is then written out together with the reads mapped to it, and bases covered by fewer than `--min-cov` reads are masked with `N`.

Two depths are involved and they are deliberately different: a variant needs `--call-min-cov` reads to be applied to the *next* reference, while the output is masked at `--min-cov`.
Keeping the first low lets the loop repair thin regions — a region that is never repaired never recruits the reads that would deepen it — while keeping the second where it is means a thinly covered call reaches the output only through the reference it shaped.

### Installation

Create the conda environment, which installs the external tools and ConsIter itself (as an editable install):
```
conda env create -f environment.yml
conda activate consiter
```

To install into an existing environment instead, make sure Bowtie2, Samtools, BedTools and GATK4 are on your `PATH`, then run `pip install .` from this directory. BCFtools is additionally required for `--caller bcftools`.

#### Software
- [Python](https://www.python.org/) >= 3.10 and [tyro](https://brentyi.github.io/tyro/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) (also used for duplicate removal, so Picard is no longer needed separately)
- [Samtools](http://www.htslib.org/)
- [BCFtools](http://www.htslib.org/) (only for `--caller bcftools`)
- [BedTools](https://bedtools.readthedocs.io/en/latest/)

`environment.yml` pins the versions that have been tested.

### Usage

```
consiter -r reference.fa -n SAMPLE -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o outdir
```

Run `consiter --help` for all options:

| Option | Default | Description |
| --- | --- | --- |
| `-r`, `--ref` | required | Reference genome FASTA (may be gzipped) |
| `-n`, `--name` | required | Sample name, used for read groups and output file names |
| `-1`, `--left` / `-2`, `--right` | required | Paired-end reads in FASTQ format |
| `-o`, `--outdir` | `.` | Output directory |
| `-i`, `--max-iter` | 10 | Maximum number of consensus rounds (the loop usually stops earlier) |
| `-t`, `--threads` | 16 | Threads for Bowtie2 and Samtools |
| `-p`, `--ploidy` | 1 | Ploidy for the variant caller (`bcftools` supports 1 or 2) |
| `-c`, `--min-cov` | 20 | Minimum read depth for a base not to be masked in the final consensus; also the depth at which breadth is reported |
| `--call-min-cov` | 2 | Minimum read depth for a variant to be applied to the next reference |
| `--caller` | `gatk` | Variant caller for the consensus rounds: `gatk` or `bcftools` |
| `--stop` | `both` | Which convergence test ends the loop: `rate`, `consensus` (fixed point or cycle) or `both` |
| `--no-dedup` | off | Keep duplicate reads; the right choice for an amplicon library |
| `--xmx` | 50 | Maximum Java heap size for GATK, in GB |
| `--keep-tmp` | off | Keep the working directory `<outdir>/tmp`, including per-step logs. ConsIter refuses to start if `<outdir>/tmp` already exists |
| `--no-indel` | off | Only apply SNPs, so the updated reference keeps the original coordinates |

By default, insertions and deletions are applied to the updated reference, so its coordinates can differ from the original reference.

This tool only works with **paired-end (PE) reads**.

#### Outputs
- `<name>.consensus.fa`: updated reference with low-coverage bases masked; sequence names match the input reference.
- `<name>.consensus.sample_name.fa`: the same sequences, named `<name>|<reference sequence name>`.
- `<name>.bt2.rmdup.bam` (and `.bai`): the reads aligned to the updated reference, with duplicates removed unless `--no-dedup` was given (the file name does not change).
- `<name>.iterations.tsv`: one row per iteration — alignment rate, duplicate fraction, mean depth, coverage breadth at `--min-cov`, bases that masking would replace with `N`, a hash of the reference, whether the round changed it, and why the loop stopped. It is rewritten after every iteration, so it is there to inspect even if a later step fails.

**Read the breadth column before trusting a consensus.** Neither a high alignment rate nor a consensus that stopped changing can tell a converged run from one whose consensus is mostly `N`: ConsIter warns when coverage is too thin to iterate from, when the duplicate fraction suggests an amplicon library, and when more than half the output is masked, but the numbers are in the TSV either way.

If a step fails, ConsIter stops with the failing program's error, and `<outdir>/tmp` is kept so the logs in `<outdir>/tmp/logs` can be inspected.

#### Choosing settings

The defaults come from the benchmarks in `benchmarks/` (see `benchmarks/FINDINGS.md`), which iterated real reads from five viruses from starting references 3-49% diverged. Four results are worth acting on:

- **`--call-min-cov` is the setting that matters most.** Separating it from `--min-cov` beat the old shared threshold of 20 in all 16 dengue cases tested (median 18% less residual error, up to 48%), and on samples where the defaults already worked it raised breadth and alignment rate without raising disagreement between runs.
- **`--caller bcftools` when the starting reference is distant.** On the hardest ladder it halved the residual error and nearly doubled breadth against HaplotypeCaller, whose accuracy there got *worse* with more reads. Where the reference is already close and the sample is a quasispecies, GATK is slightly better, which is why it stays the default.
- **`--stop consensus` with `--caller bcftools`.** The rate test is reliable with GATK (it chose the best iteration in every simulated run) but cuts a pileup caller short: the same test cost bcftools a mean of 13 extra residual errors, and on the dengue ladder a pileup caller was still gaining a median of 92 bases at iteration 5. The cost is wasted iterations, since a drifting caller can change a few bases per round without ever repeating a consensus.
- **`--no-dedup` for amplicon libraries.** On an ARTIC library, duplicate removal discarded 99.1% of read pairs and the result was an almost entirely `N` consensus that every stopping signal reported as converged.

#### Changes from the 0.1 script
- Installed as the `consiter` command instead of `python3 ConsIter.py`.
- `-ref` is now `-r`/`--ref`, `-xmx` is now `--xmx`, `--keep` is now `--keep-tmp`, and `--noindel` is now `--no-indel`.
- The final BAM and the coverage mask now come from reads mapped to the final consensus (in 0.1, they came from the previous reference).
- The variant depth filter now uses `--call-min-cov` (it was fixed at 20, then `--min-cov`), and the SOR strand-bias filter was removed.
- `--max-iter` defaults to 10 rather than 5, and the loop also stops at a fixed point or a cycle.
- New options: `--call-min-cov`, `--caller`, `--stop` and `--no-dedup`; new output: `<name>.iterations.tsv`.

### Development

```
pytest                  # unit tests
pytest -m integration   # end-to-end runs on simulated data (needs the conda environment)
ruff check . && ruff format --check .
```

### Terms of Use

By using this software, you agree this software is to be used for research purposes only. Any presentation of data analysis using the software will acknowledge the software according to the guidelines below.

Primary author(s): David B. Stern

Organizational contact information: david.stern AT nih.gov

Date of release: XXX

Version: 0.2.0

License details: see LICENSE file

### Disclaimer

A review of this code has been conducted, no critical errors exist, and to the best of the authors knowledge, there are no problematic file paths, no local system configuration details, and no passwords or keys included in this code. This open source software comes as is with absolutely no warranty.
