# assembly_updater
Python script to improve read alignment to a reference genome by updating the reference using read mapping data

Maps reads to a reference genome with Bowtie2, calls variants with GATK4, and generates a new consensus sequence.
This process with continue in an iterative manner until there is no improve in alignment rate, or until the maximum number of iterations is met.

Please update the python script to specify the path for bowtie2, picard, GATK4, and samtools.
