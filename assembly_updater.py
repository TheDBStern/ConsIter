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
    cmd = '%s-build -q %s %s &> /dev/null' % (bowtie2_cmd, ref_fasta, ref_fasta)
    os.system(cmd)

def run_bowtie2(ref, iter):
    #run bowtie2 and collect alignment stats
    cmd1 = '%s -p %s \
            --no-unal \
            --rg-id tmp\
            --rg SM:tmp \
            --rg LB:1 \
            --rg PU:1 \
            --rg PL:illumina \
            --sensitive-local \
            -x %s  \
            -1 %s \
            -2 %s \
            -S tmp/iter%s.sam 2> tmp/iter%s.btstats.txt'%(bowtie2_cmd, args.Threads, ref, args.Left, args.Right, iter, iter)
    #convert to sorted bam
    cmd2 = '%s view -bh \
            tmp/iter%s.sam | \
            samtools sort -@ %s\
            - > \
            tmp/iter%s.bam'%(samtools_cmd, iter,args.Threads,iter)
    os.system(cmd1)
    os.system(cmd2)

def align_rate(iter):
    with open('tmp/iter%s.btstats.txt'%iter, 'r') as fh:
        bt2str = fh.read()
        m = re.search('(\d+\.\d+)\% overall alignment rate', bt2str)
        alnrt = m.group(1)
        return(float(alnrt))

def rmdup(iter):
    cmd = 'java -jar \
            %s MarkDuplicates \
            CREATE_INDEX=true \
            USE_JDK_DEFLATER=true \
            USE_JDK_INFLATER=true \
            M=tmp/iter%s.rmdup_metrics.txt \
            I=tmp/iter%s.bam \
            O=tmp/iter%s.rmdup.bam \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT'%(picard_cmd, iter,iter,iter)
    os.system(cmd)


def call_variants(ref, iter):
    refname = '.'.join(ref.split('.')[:-1])
    cmd1 = "samtools faidx %s"%ref
    cmd2 = "samtools dict %s > %s.dict"%(ref,refname)
    cmd3 = "%s --java-options '-Xmx%sg' HaplotypeCaller  \
       --use-jdk-deflater --use-jdk-inflater \
       -R %s \
       -I tmp/iter%s.rmdup.bam \
       -O tmp/iter%s.vcf \
       --min-base-quality-score 20 \
       -ploidy 1"%(gatk_cmd, args.xmx, ref, iter, iter)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)

def consensus(ref,iter):
    cmd1 = "%s IndexFeatureFile -I tmp/iter%s.vcf"%(gatk_cmd, iter)
    cmd2 = "%s FastaAlternateReferenceMaker \
       -R %s \
       -O tmp/consensus.iter%s.fa \
       -V tmp/tmp/iter%s.vcf"%(gatk_cmd,ref,iter,iter)
    os.system(cmd1)
    os.system(cmd2)

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

    #run main loop
    iteration = 0

    #reached last iteration
    if iteration == args.maxIter:
        print("Maximum number of iterations reached")
        cmd = ("mv tmp/consensus.iter%s.fa updated_reference.fa"%iteration-1)
        os.system(cmd)
    # in first iteration, map to original reference genome
    elif iteration == 0:
        print("Indexing reference")
        run_bowtie2_idx(args.ref)
        print("Mapping reads")
        run_bowtie2(args.ref, iteration)
        alnrate = align_rate(iteration)
        print("Removing duplicates")
        rmdup(iteration)
        print("Iteration %s alignment rate: %s"%(iteration,alnrate))
        print("Calling variants")
        call_variants(args.ref,iteration)
        print("Generating updated reference")
        consensus(args.ref,iteration)
        iteration +=1
    else:
        #map to updated reference and check if alignment rate is better than last iteration
        print("Indexing updated reference")
        run_bowtie2_idx("tmp/consensus.iter%s.fa"%iteration)
        print("Mapping reads")
        run_bowtie2("tmp/consensus.iter%s.fa"%iteration, iteration)
        alnrate_last = align_rate(iteration-1)
        alnrate = align_rate(iteration)
        print("Iteration %s alignment rate: %s"%(iteration-1,alnrate_last))
        print("Iteration %s alignment rate: %s"%(iteration,alnrate))
        if alnrate > alnrate_last:
            print("Iteration %s alignment rate better than previous iteration\n \
                    Generating new consensus..."%(iteration))
            rmdup(iteration)
            call_variants(args.ref,iteration)
            consensus(args.ref,iteration)
            iteration +=1
        else:
            print("No improvement in alignment rate")
            print("Terminating")
            cmd = ("mv tmp/consensus.iter%s.fa updated_reference.fa"%iteration-1)
            os.system(cmd)
            iteration +=1
