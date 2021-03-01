#!/usr/bin/env python

import argparse
import subprocess
import re
import os
import sys

# set base commands for software on your system
bowtie2_cmd = "bowtie2"
samtools_cmd = "samtools"
gatk_cmd = "$EBROOTGATK/gatk"
picard_cmd = "java -jar ${EBROOTPICARD}/picard.jar"

def run_bowtie2_idx(ref_fasta):
    cmd = '%s-build -q %s %s &> /dev/null' % (bowtie2_cmd, ref_fasta, ref_fasta)
    os.system(cmd)

def run_bowtie2(ref, iter):
    #run bowtie2 and collect alignment stats
    cmd1 = '%s -p %s \
            --no-unal \
            --rg-id %s\
            --rg SM:%s \
            --rg LB:1 \
            --rg PU:1 \
            --rg PL:illumina \
            --sensitive-local \
            -x %s  \
            -1 %s \
            -2 %s \
            -S %s/tmp/iter%s.sam 2> %s/tmp/iter%s.btstats.txt'%(bowtie2_cmd, args.Threads, args.name, args.name, ref, args.Left, args.Right, args.outdir, iter, args.outdir, iter)
    #convert to sorted bam
    cmd2 = '%s view -bh \
            %s/tmp/iter%s.sam | \
            samtools sort -@ %s\
            -o %s/tmp/iter%s.bam \
            - &> /dev/null'%(samtools_cmd, args.outdir, iter, args.Threads, args.outdir, iter)
    os.system(cmd1)
    os.system(cmd2)

def align_rate(iter):
    with open('%s/tmp/iter%s.btstats.txt'%(args.outdir,iter), 'r') as fh:
        bt2str = fh.read()
        m = re.search('(\d+\.\d+)\% overall alignment rate', bt2str)
        alnrt = m.group(1)
        return(float(alnrt))

def rmdup(iter):
    cmd = '%s MarkDuplicates \
            CREATE_INDEX=true \
            USE_JDK_DEFLATER=true \
            USE_JDK_INFLATER=true \
            M=%s/tmp/iter%s.rmdup_metrics.txt \
            I=%s/tmp/iter%s.bam \
            O=%s/tmp/iter%s.rmdup.bam \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT &> /dev/null'%(picard_cmd, args.outdir, iter, args.outdir, iter, args.outdir, iter)
    os.system(cmd)


def call_variants(ref, iter):
    refname = '.'.join(ref.split('.')[:-1])
    cmd1 = "samtools faidx %s"%ref
    cmd2 = "samtools dict %s > %s.dict"%(ref,refname)
    cmd3 = "%s --java-options '-Xmx%sg' HaplotypeCaller  \
       --use-jdk-deflater --use-jdk-inflater \
       -R %s \
       -I %s/tmp/iter%s.rmdup.bam \
       -O %s/tmp/iter%s.vcf \
       --min-base-quality-score 20 \
       -ploidy 1 &> /dev/null"%(gatk_cmd, args.xmx, ref, args.outdir, iter, args.outdir, iter)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)

def consensus(ref,iter):
    cmd1 = "%s IndexFeatureFile -I %s/tmp/iter%s.vcf &> /dev/null"%(gatk_cmd, args.outdir, iter)
    cmd2 = "%s FastaAlternateReferenceMaker \
       -R %s \
       -O %s/tmp/consensus.iter%s.fa \
       -V %s/tmp/iter%s.vcf &> /dev/null"%(gatk_cmd,ref,args.outdir,iter,args.outdir,iter)
    os.system(cmd1)
    os.system(cmd2)

def check_software(command):
    status = subprocess.getstatusoutput(command)
    m = re.search("command not found", str(status[1]))
    if m is None:
        return True
    else:
        return False

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Python script to improve read alignment to a reference genome by updating the reference using read mapping data')
    parser.add_argument('-ref', dest = 'ref', type = str, required=True,  help = 'Path to reference genome')
    parser.add_argument('-n', dest = 'name', type = str, required=True,  help = 'Sample name')
    parser.add_argument('-o', dest = 'outdir', type = str, required=True,  help = 'Output directory')
    parser.add_argument('-1', dest= 'Left', type = str, required=True, help ='Left reads file in fastq format')
    parser.add_argument('-2', dest= 'Right', type = str, required=True, help ='Right reads file in fastq format')
    parser.add_argument('-i', dest= 'maxIter', type=int, default= 5, help ='Maximum number of iterations. Default = 5')
    parser.add_argument('-t', dest= 'Threads', type = int, default= 16, help ='Number of threads to use. Default = 16')
    parser.add_argument('-xmx', dest= 'xmx', type = int, default= 50, help ='Maximum heap size for Java VM, in GB. Default = 50')
    parser.add_argument('--keep',dest= 'keep', action='store_true',help ='Keep temporary directory')
    args = parser.parse_args()

    #make output directory
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    #tmp directory
    if not os.path.isdir(args.outdir+"/tmp"):
        os.mkdir(args.outdir+"/tmp")


    #check that all software works
    if shutil.which(bowtie2_cmd) == None:
        print(sys.exit("Could not find bowtie2 command. Adjust python script"))
    if shutil.which(samtools_cmd) == None:
        print(sys.exit("Could not find Samtools command. Adjust python script"))
    """
    if shutil.which(gatk_cmd) == None:
        print("Could not find GATK command. Adjust python script")
    if shutil.which(picard_cmd) == None:
        print("Could not find Picard command. Adjust python script")
    """

	#Check for valid file format and parameters

    if args.maxIter<=0:
        print(sys.exit("Error: Number of iterations must be greater than 0"))
    if not os.path.isfile(args.ref):
        print(sys.exit("Error: Cannot find reference fasta"))
    if not os.path.isfile(args.Left):
        print(sys.exit("Error: Cannot find left reads file"))
    if not os.path.isfile(args.Right):
        print(sys.exit("Error: Cannot find right reads file"))


    #### run main loop ####
    iteration = 0

    while iteration <= args.maxIter:
        #reached last iteration
        if iteration == args.maxIter:
            print("Iteration %s"%iteration)
            print("Maximum number of iterations reached. Terminating.")
            cmd = ("mv %s/tmp/consensus.iter%s.fa %s/%s.updated_reference.fa"%(args.outdir,(iteration-1),args.outdir,args.name))
            os.system(cmd)
            if not args.keep:
                os.system('rm -rf %s/tmp'%args.outdir)
            iteration +=1
        # in first iteration, map to original reference genome
        elif iteration == 0:
            print("Iteration %s"%iteration)
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
            print("Iteration %s"%iteration)
            print("Indexing updated reference")
            run_bowtie2_idx("%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)))
            print("Mapping reads")
            run_bowtie2("%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration)
            alnrate_last = align_rate(iteration-1)
            alnrate = align_rate(iteration)
            print("Iteration %s alignment rate: %s"%((iteration-1),alnrate_last))
            print("Iteration %s alignment rate: %s"%(iteration,alnrate))
            if alnrate > alnrate_last:
                print("Iteration %s alignment rate better than previous iteration"%(iteration))
                print("Removing duplicates")
                rmdup(iteration)
                print("Calling variants")
                call_variants("%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration)
                print("Generating updated reference")
                consensus("%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration)
                iteration +=1
            else:
                print("No improvement in alignment rate")
                print("Terminating")
                print("Updated reference genome is in %s.updated_reference.fa"%args.name)
                cmd = ("mv %s/tmp/consensus.iter%s.fa %s/%s.updated_reference.fa"%(args.outdir,(iteration-1),args.outdir,args.name))
                os.system(cmd)
                if not args.keep:
                    os.system('rm -rf %s/tmp'%args.outdir)
                break
