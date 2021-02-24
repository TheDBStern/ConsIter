#!/usr/bin/env python

import argparse
import shutil
import re
import os

# set base commands for software on your system
bowtie2_cmd = "bowtie2"
samtools_cmd = "samtools"
gatk_cmd = "$EBROOTGATK/gatk"
picard_cmd = "${EBROOTPICARD}/picard.jar"

def run_bowtie2_idx(ref_fasta):
    cmd = '%s-build -q %s temp/%s &> /dev/null' % (bowtie2_cmd, ref_fasta, target_name)
    os.system(cmd)

def run_bowtie2(sample_name, ref, iter):
    #run bowtie2
    cmd1 = '%s -p %s \
            --no-unal \
            --rg-id %s \
            --rg SM:%s \
            --rg LB:1 \
            --rg PU:1 \
            --rg PL:illumina \
            --sensitive-local \
            -x %s  \
            -1 %s \
            -2 %s \
            -S tmp/%s.iter%s.sam 2> tmp/btstats.iter%s.txt'%(bowtie2_cmd, args.Threads, sample_name, sample_name, ref, args.Left, args.Right, iter, iter)
    #convert to sorted bam
    cmd2 = '%s view -bh \
            tmp/%s.iter%s.sam | \
            samtools sort -@ %s\
            - > \
            tmp/%s.iter%s.bam'%(samtools_cmd, sample_name,iter,args.Threads,sample_name,iter)
    os.system(cmd1)
    os.system(cmd2)

def align_rate(iter):
    with open('tmp/btstats.iter%s.txt'%iter, 'r') as fh:
        bt2str = fh.read()
        m = re.search('(\d+\.\d+)\% overall alignment rate', bt2str)
        alnrt = m.group(1)
        return(float(alnrt))

def rmdup():
    cmd = 'java -jar \
            %s MarkDuplicates \
            CREATE_INDEX=true \
            USE_JDK_DEFLATER=true \
            USE_JDK_INFLATER=true \
            M=SRR8525886/SRR8525886.aligned_consensus.bt2.rmdup_metrics.txt \
            I=SRR8525886/SRR8525886.aligned_consensus.bt2.bam \
            O=SRR8525886/SRR8525886.aligned_consensus.bt2.rmdup.bam \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT'%(picard_cmd, )

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Python script to improve read alignment to a reference genome by updating the reference using read mapping data')
    parser.add_argument('-ref', dest = 'ref', type = str, required=True,  help = 'Path to reference genome')
    parser.add_argument('-1', dest= 'Left', type = str, required=True, help ='Left reads file in fastq format')
    parser.add_argument('-2', dest= 'Right', type = str, required=True, help ='Right reads file in fastq format')
    parser.add_argument('-o', dest= 'outdir', type = str, required=True, help ='Output directory')
    parser.add_argument('-i', dest= 'maxIter', type = int, default= '5', help ='Maximum number of iterations. Default = 5')
    parser.add_argument('-t', dest= 'Threads', type = int, default= 16, help ='Number of threads to use. Default = 16')
    parser.add_argument('-xmx', dest= 'xmx', type = int, default= 50, help ='Maximum heap size for Java VM, in GB. Default = 50')
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    if not os.path.isdir("%s/tmp"%args.outdir):
        os.mkdir("%s/tmp"%args.outdir)


    #check that all software works
    if shutil.which(bowtie2_cmd) == None:
        print("Could not find bowtie2 command. Adjust python script")
    if shutil.which(samtools_cmd) == None:
        print("Could not find Samtools command. Adjust python script")
    if shutil.which(gatk_cmd) == None:
        print("Could not find GATK command. Adjust python script")
    if shutil.which(picard_cmd) == None:
        print("Could not find Picard command. Adjust python script")
