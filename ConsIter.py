#!/usr/bin/env python

import argparse
import subprocess
import shutil
import re
import os
import sys

# set base commands for software on your system
######################################################
bowtie2_cmd = "bowtie2"
samtools_cmd = "samtools"
gatk_cmd = "$EBROOTGATK/gatk"
picard_cmd = "java -jar ${EBROOTPICARD}/picard.jar"
######################################################

def run_bowtie2_idx(bt2_cmd,ref_fasta):
    cmd = '%s-build -q %s %s &> /dev/null' % (bt2_cmd, ref_fasta, ref_fasta)
    os.system(cmd)

def run_bowtie2(bt2_cmd, smt_cmd, ref, iter, threads, name, left, right, outdir):
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
            -S %s/tmp/iter%s.sam 2> %s/tmp/iter%s.btstats.txt'%(bt2_cmd, threads, name, name, ref, left, right, outdir, iter, outdir, iter)
    #convert to sorted bam
    cmd2 = '%s view -bh \
            %s/tmp/iter%s.sam | \
            samtools sort -@ %s\
            -o %s/tmp/iter%s.bam \
            - &> /dev/null'%(smt_cmd, outdir, iter, threads, outdir, iter)
    os.system(cmd1)
    os.system(cmd2)

def align_rate(iter, outdir):
    with open('%s/tmp/iter%s.btstats.txt'%(outdir,iter), 'r') as fh:
        bt2str = fh.read()
        m = re.search('(\d+\.\d+)\% overall alignment rate', bt2str)
        alnrt = m.group(1)
        return(float(alnrt))

def rmdup(pic_cmd, iter, outdir):
    cmd = '%s MarkDuplicates \
            CREATE_INDEX=true \
            USE_JDK_DEFLATER=true \
            USE_JDK_INFLATER=true \
            M=%s/tmp/iter%s.rmdup_metrics.txt \
            I=%s/tmp/iter%s.bam \
            O=%s/tmp/iter%s.rmdup.bam \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT &> /dev/null'%(pic_cmd, outdir, iter, outdir, iter, outdir, iter)
    os.system(cmd)


def call_variants(gt_cmd, smt_cmd, ref, iter, xmx, outdir):
    refname = '.'.join(ref.split('.')[:-1])
    cmd1 = "%s faidx %s"%(smt_cmd, ref)
    cmd2 = "%s dict %s > %s.dict"%(smt_cmd,ref,refname)
    cmd3 = "%s --java-options '-Xmx%sg' HaplotypeCaller  \
       --use-jdk-deflater --use-jdk-inflater \
       -R %s \
       -I %s/tmp/iter%s.rmdup.bam \
       -O %s/tmp/iter%s.vcf \
       --min-base-quality-score 20 \
       -ploidy 1 &> /dev/null"%(gt_cmd, xmx, ref, outdir, iter, outdir, iter)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)

def consensus(ref,iter,vcf, outdir):
    cmd1 = "%s IndexFeatureFile -I %s &> /dev/null"%(gatk_cmd, vcf)
    cmd2 = "%s FastaAlternateReferenceMaker \
       -R %s \
       -O %s/tmp/consensus.iter%s.tmpnames.fa \
       -V %s &> /dev/null"%(gatk_cmd,ref,outdir,iter,vcf)
    os.system(cmd1)
    os.system(cmd2)

def select_snps(ref,iter, outdir):
    cmd = "%s SelectVariants \
        --use-jdk-deflater --use-jdk-inflater \
        -R %s \
        -V %s/tmp/iter%s.vcf \
        -select-type SNP \
        -O %s/tmp/iter%s.snps.vcf &> /dev/null"%(gatk_cmd,ref,outdir,iter,outdir,iter)
    os.system(cmd)

def rename(consensus,ref,outfasta):
    fasta = open(consensus, 'r').read()
    new = open(outfasta,'w')
    refdict = open('.'.join(ref.split('.')[:-1])+'.dict','r')
    next(refdict)
    count=1
    for line in refdict:
        name = line.split('\t')[1].replace("SN:","")
        fasta = fasta.replace(">%s"%count, ">%s"%name)
        count+=1
    new.write(fasta)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='ConsIter produces an updated reference genome that more closely matches the sequenced viral population by iteratively mapping reads to the reference and generating a new consensus sequence.')
    parser.add_argument('-ref', dest = 'ref', type = str, required=True,  help = 'Path to reference genome')
    parser.add_argument('-n', dest = 'name', type = str, required=True,  help = 'Sample name')
    parser.add_argument('-o', dest = 'outdir', type = str, required=True,  default=".", help = 'Output directory. Default is the current directory')
    parser.add_argument('-1', dest= 'Left', type = str, required=True, help ='Left reads file in fastq format')
    parser.add_argument('-2', dest= 'Right', type = str, required=True, help ='Right reads file in fastq format')
    parser.add_argument('-i', dest= 'maxIter', type=int, default= 5, help ='Maximum number of iterations. Default = 5')
    parser.add_argument('-t', dest= 'Threads', type = int, default= 16, help ='Number of threads to use. Default = 16')
    parser.add_argument('-xmx', dest= 'xmx', type = int, default= 50, help ='Maximum heap size for Java VM, in GB. Default = 50')
    parser.add_argument('--keep',dest= 'keep', action='store_true',help ='Keep temporary directory')
    parser.add_argument('--noindel',dest= 'noindel', action='store_true',help ='Do not introduce insertions and deletions into new reference')
    args = parser.parse_args()

    #make output directory
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    #tmp directory
    if not os.path.isdir(args.outdir+"/tmp"):
        os.mkdir(args.outdir+"/tmp")


    #check that all software works
    if shutil.which(bowtie2_cmd) == None:
        print(sys.exit("Could not find bowtie2 command. Adjust python script or load module"))
    if shutil.which(samtools_cmd) == None:
        print(sys.exit("Could not find Samtools command. Adjust python script or load module"))
    gatk_stat = subprocess.getstatusoutput(gatk_cmd)
    gatk_m = re.search("Usage", str(gatk_stat[1]))
    if gatk_m == None:
        print(sys.exit("Could not find GATK command. Adjust python script or load module"))
    picard_stat = subprocess.getstatusoutput(picard_cmd)
    picard_m = re.search("PicardCommandLine", str(picard_stat[1]))
    if picard_m == None:
        print(sys.exit("Could not find Picard command. Adjust python script or load module"))

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
            print("Updated reference genome is in %s/%s.updated_reference.fa"%(args.outdir,args.name))
            cmd1 = ("mv %s/tmp/consensus.iter%s.fa %s/%s.updated_reference.fa"%(args.outdir,(iteration-1),args.outdir,args.name))
            os.system(cmd1)
            print("Duplicate-removed bam file is in %s/%s.bt2.rmdup.bam"%(args.outdir,args.name))
            cmd2 = ("mv %s/tmp/iter%s.rmdup.bam %s/%s.bt2.rmdup.bam"%(args.outdir,(iteration-1),args.outdir,args.name))
            os.system(cmd2)
            cmd3 = ("mv %s/tmp/iter%s.rmdup.bai %s/%s.bt2.rmdup.bai"%(args.outdir,(iteration-1),args.outdir,args.name))
            os.system(cmd3)
            if not args.keep:
                os.system('rm -rf %s/tmp'%args.outdir)
            iteration +=1
        # in first iteration, map to original reference genome
        elif iteration == 0:
            print("Iteration %s"%iteration)
            print("Indexing reference")
            run_bowtie2_idx(bowtie2_cmd,args.ref)
            print("Mapping reads")
            run_bowtie2(bowtie2_cmd, samtools_cmd, args.ref, iteration, args.Threads, args.name, args.Left, args.Right, args.outdir)
            alnrate = align_rate(iteration, args.outdir)
            print("Removing duplicates")
            rmdup(picard_cmd,iteration, args.outdir)
            print("Iteration %s alignment rate: %s"%(iteration,alnrate))
            print("Calling variants")
            call_variants(gatk_cmd, samtools_cmd, args.ref,iteration, args.xmx, args.outdir)
            if args.noindel:
                print("Generating updated reference")
                select_snps(args.ref,iteration, args.outdir)
                consensus(args.ref,iteration,"%s/tmp/iter%s.snps.vcf"%(args.outdir,iteration), args.outdir)
                rename('%s/tmp/consensus.iter%s.tmpnames.fa'%(args.outdir,iteration),args.ref,'%s/tmp/consensus.iter%s.fa'%(args.outdir,iteration))
                iteration +=1
            else:
                print("Generating updated reference")
                consensus(args.ref,iteration,"%s/tmp/iter%s.vcf"%(args.outdir,iteration), args.outdir)
                rename('%s/tmp/consensus.iter%s.tmpnames.fa'%(args.outdir,iteration),args.ref,'%s/tmp/consensus.iter%s.fa'%(args.outdir,iteration))
                iteration +=1
        else:
            #map to updated reference and check if alignment rate is better than last iteration
            print("Iteration %s"%iteration)
            print("Indexing updated reference")
            run_bowtie2_idx(bowtie2_cmd,"%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)))
            print("Mapping reads")
            run_bowtie2(bowtie2_cmd, samtools_cmd, "%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration, args.Threads, args.name, args.Left, args.Right, args.outdir)
            alnrate_last = align_rate(iteration-1, args.outdir)
            alnrate = align_rate(iteration, args.outdir)
            print("Iteration %s alignment rate: %s"%((iteration-1),alnrate_last))
            print("Iteration %s alignment rate: %s"%(iteration,alnrate))
            if alnrate > alnrate_last:
                print("Iteration %s alignment rate better than previous iteration"%(iteration))
                print("Removing duplicates")
                rmdup(picard_cmd,iteration, args.outdir)
                print("Calling variants")
                call_variants(gatk_cmd, samtools_cmd, "%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration, args.xmx, args.outdir)
                if args.noindel:
                    print("Generating updated reference")
                    select_snps("%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration, args.outdir)
                    consensus("%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration,"%s/tmp/iter%s.snps.vcf"%(args.outdir,iteration), args.outdir)
                    rename('%s/tmp/consensus.iter%s.tmpnames.fa'%(args.outdir,iteration),args.ref,'%s/tmp/consensus.iter%s.fa'%(args.outdir,iteration))
                    iteration +=1
                else:
                    print("Generating updated reference")
                    consensus("%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration,"%s/tmp/iter%s.vcf"%(args.outdir,iteration), args.outdir)
                    rename('%s/tmp/consensus.iter%s.tmpnames.fa'%(args.outdir,iteration),args.ref,'%s/tmp/consensus.iter%s.fa'%(args.outdir,iteration))
                    iteration +=1
            else:
                print("No improvement in alignment rate")
                print("Generating final consensus sequence and bam file")
                print("Updated reference genome is in %s/%s.updated_reference.fa"%(args.outdir,args.name))
                cmd1 = ("mv %s/tmp/consensus.iter%s.fa %s/%s.updated_reference.fa"%(args.outdir,(iteration-1),args.outdir,args.name))
                os.system(cmd1)
                print("Duplicate-removed bam file is in %s/%s.bt2.rmdup.bam"%(args.outdir,args.name))
                cmd2 = ("mv %s/tmp/iter%s.rmdup.bam %s/%s.bt2.rmdup.bam"%(args.outdir,(iteration-1),args.outdir,args.name))
                os.system(cmd2)
                cmd3 = ("mv %s/tmp/iter%s.rmdup.bai %s/%s.bt2.rmdup.bai"%(args.outdir,(iteration-1),args.outdir,args.name))
                os.system(cmd3)
                if not args.keep:
                    os.system('rm -rf %s/tmp'%args.outdir)
                break
