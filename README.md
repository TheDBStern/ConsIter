# ConsIter
Evolving viral populations often diverge substantially from available reference genomes, potentially reducing read mapping accuracy and biasing variant calls. 
This Python script produces an updated reference genome that more closely matches the viral population by iteratively mapping reads to the reference and generating a new consensus sequence.

Reads are mapped to a reference genome with Bowtie2, variants are called with GATK4 HaplotypeCaller, and a new consensus sequence is generated from the variant calls.
This process continues in an iterative manner until there is no improvement in alignment rate, or until the specified maximum number of iterations is reached.

Please update the python script to specify your system path for bowtie2, picard, GATK4, and samtools.

Help menu is available by typing: `ConsIter.py -h`

Currently, this script is able to add insertions and deletion to the updated reference causing positional mismatches between the original and updated reference. To disable with behavior, use the flag `--noindel`

To install the appropriate packages in a conda environment, run 
```
conda create -n ConsIter bowtie2=2.4.4-0 gatk4=4.2.2.0-0 picard=2.26.0-0 samtools=1.13-0
conda activate ConsIter
```
