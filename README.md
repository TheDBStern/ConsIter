# ConsIter
Evolving viral populations often diverge substantially from available reference genomes, potentially reducing read mapping accuracy and biasing variant calls.
ConsIter produces an updated reference genome that more closely matches the viral population by iteratively mapping reads to the reference and generating a new consensus sequence.

In each iteration, paired-end reads are mapped to the current reference with Bowtie2, duplicates are removed, variants are called with GATK4 HaplotypeCaller and hard-filtered, and the passing variants are applied to make the next reference.
This continues until the alignment rate stops improving, or until the maximum number of consensus rounds is reached.
The latest consensus is then written out together with the reads mapped to it, and bases covered by fewer than a given number of reads are masked with `N`.

### Installation

Create the conda environment, which installs the external tools and ConsIter itself (as an editable install):
```
conda env create -f environment.yml
conda activate consiter
```

To install into an existing environment instead, make sure Bowtie2, Samtools, BedTools and GATK4 are on your `PATH`, then run `pip install .` from this directory.

#### Software
- [Python](https://www.python.org/) >= 3.10 and [tyro](https://brentyi.github.io/tyro/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) (also used for duplicate removal, so Picard is no longer needed separately)
- [Samtools](http://www.htslib.org/)
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
| `-i`, `--max-iter` | 5 | Maximum number of consensus rounds |
| `-t`, `--threads` | 16 | Threads for Bowtie2 and Samtools |
| `-p`, `--ploidy` | 1 | Ploidy for GATK HaplotypeCaller |
| `-c`, `--min-cov` | 20 | Mask bases covered by fewer reads than this |
| `--xmx` | 50 | Maximum Java heap size for GATK, in GB |
| `--keep-tmp` | off | Keep the working directory `<outdir>/tmp`, including per-step logs |
| `--no-indel` | off | Only apply SNPs, so the updated reference keeps the original coordinates |

By default, insertions and deletions are applied to the updated reference, so its coordinates can differ from the original reference.

This tool only works with **paired-end (PE) reads**.

#### Outputs
- `<name>.consensus.fa`: updated reference with low-coverage bases masked; sequence names match the input reference.
- `<name>.consensus.sample_name.fa`: the same sequences, named `<name>|<reference sequence name>`.
- `<name>.bt2.rmdup.bam` (and `.bai`): duplicate-removed reads aligned to the updated reference.

If a step fails, ConsIter stops with the failing program's error, and `<outdir>/tmp` is kept so the logs in `<outdir>/tmp/logs` can be inspected.

#### Changes from the 0.1 script
- Installed as the `consiter` command instead of `python3 ConsIter.py`.
- `-ref` is now `-r`/`--ref`, `-xmx` is now `--xmx`, `--keep` is now `--keep-tmp`, and `--noindel` is now `--no-indel`.
- The final BAM and the coverage mask now come from reads mapped to the final consensus (in 0.1, they came from the previous reference).

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
