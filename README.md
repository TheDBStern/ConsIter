# ConsIter
Evolving viral populations often diverge substantially from available reference genomes, potentially reducing read mapping accuracy and biasing variant calls.
This Python script produces an updated reference genome that more closely matches the viral population by iteratively mapping reads to the reference and generating a new consensus sequence.

Reads are mapped to a reference genome with Bowtie2, variants are called with GATK4 HaplotypeCaller, and a new consensus sequence is generated from the variant calls.
This process continues in an iterative manner until there is no improvement in alignment rate, or until the specified maximum number of iterations is reached.

Bases in the final consensus fasta are masked if they are covered by fewer than a specified number of reads.

### Installation

To install the appropriate packages in a conda environment, run
```
conda create -n consiter -y python=3.8 bowtie2=2.4.4-0 gatk4=4.2.2.0-0 picard=2.26.0-0 samtools=1.13-0 bedtools=2.30.0
conda activate consiter
```

#### Software  
- [Python 3](https://www.python.org/download/releases/3.0/)  
- [Bowtie2 v2.4.4](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)  
- [GATK v4.2.2](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)  
- [Picard v2.26](https://broadinstitute.github.io/picard/)  
- [Samtools v1.13](http://www.htslib.org/)  
- [BedTools v2.30](https://bedtools.readthedocs.io/en/latest/)  


### Usage

If not using the above conda environment, please update the python script to specify your system path for bowtie2, picard, GATK4, and samtools.

Help menu is available by typing: `python3 ConsIter.py -h`

Currently, this script is able to add insertions and deletion to the updated reference causing positional mismatches between the original and updated reference. To disable with behavior, use the flag `--noindel`

This script only works with **Paired-End (PE) reads**.

### Terms of Use

By using this software, you agree this software is to be used for research purposes only. Any presentation of data analysis using the software will acknowledge the software according to the guidelines below.

Primary author(s): David B. Stern

Organizational contact information: david.stern AT nih.gov

Date of release: XXX

Version: 0.1

License details: see LICENSE file

### Disclaimer

A review of this code has been conducted, no critical errors exist, and to the best of the authors knowledge, there are no problematic file paths, no local system configuration details, and no passwords or keys included in this code. This open source software comes as is with absolutely no warranty.
