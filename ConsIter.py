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
gatk_cmd = "gatk"
picard_cmd = "picard"
bedtools_cmd = "bedtools"
######################################################

def run_bowtie2_idx(bt2_cmd,ref_fasta):
    cmd = '%s-build -q %s %s \
            &> /dev/null' % (bt2_cmd, ref_fasta, ref_fasta)
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
            -S %s/tmp/iter%s.sam \
            2> %s/tmp/iter%s.btstats.txt'%(bt2_cmd, threads, name, name, ref, left, right, outdir, iter, outdir, iter)
    #convert to sorted bam
    cmd2 = '%s view -bh \
            -f 0x2 \
            -q 20 \
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
            VALIDATION_STRINGENCY=LENIENT \
            &> /dev/null'%(pic_cmd, outdir, iter, outdir, iter, outdir, iter)
    os.system(cmd)


def call_variants(gt_cmd, smt_cmd, ref, iter, xmx, outdir, ploidy):
    refname = '.'.join(ref.split('.')[:-1])
    cmd1 = "%s faidx %s"%(smt_cmd, ref)
    cmd2 = "%s dict %s > %s.dict"%(smt_cmd,ref,refname)
    cmd3 = "%s --java-options '-Xmx%sg' HaplotypeCaller  \
       --use-jdk-deflater --use-jdk-inflater \
       -R %s \
       -I %s/tmp/iter%s.rmdup.bam \
       -O %s/tmp/iter%s.vcf \
       --min-base-quality-score 20 \
       -ploidy %s \
       &> /dev/null"%(gt_cmd, xmx, ref, outdir, iter, outdir, iter, ploidy)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)

def consensus(gt_cmd,ref,iter,vcf, outdir):
    cmd1 = "%s IndexFeatureFile -I %s &> /dev/null"%(gt_cmd, vcf)
    cmd2 = "%s FastaAlternateReferenceMaker \
       -R %s \
       -O %s/tmp/consensus.iter%s.tmpnames.fa \
       -V %s \
       &> /dev/null"%(gt_cmd,ref,outdir,iter,vcf)
    os.system(cmd1)
    os.system(cmd2)

def select_snps(gt_cmd, ref,iter, outdir):
    cmd = "%s SelectVariants \
        --use-jdk-deflater --use-jdk-inflater \
        -R %s \
        -V %s/tmp/iter%s.filt.vcf \
        -select-type SNP \
        -O %s/tmp/iter%s.snps.vcf &> \
        /dev/null"%(gt_cmd,ref,outdir,iter,outdir,iter)
    os.system(cmd)

def rename_sname(consensus,ref,outfasta,sample):
    fasta = open(consensus, 'r').read()
    new = open(outfasta,'w')
    refdict = open('.'.join(ref.split('.')[:-1])+'.dict','r')
    next(refdict)
    count=1
    for line in refdict:
        name = line.split('\t')[1].replace("SN:","")
        fasta = fasta.replace(">", ">%s|%s"%(sample,name))
        count+=1
    new.write(fasta)

def rename(consensus,ref,outfasta):
    fasta = open(consensus, 'r').read()
    new = open(outfasta,'w')
    refdict = open('.'.join(ref.split('.')[:-1])+'.dict','r')
    next(refdict)
    count=1
    for line in refdict:
        name = line.split('\t')[1].replace("SN:","")
        fasta = fasta.replace(">%s"%count, ">%s"%(name))
        count+=1
    new.write(fasta)

def filter_variants(gt_cmd, ref, iter, outdir):
    cmd = '%s VariantFiltration \
        -R %s \
        -V %s/tmp/iter%s.vcf \
        -O %s/tmp/iter%s.filt.vcf \
        -filter-name "DP_filter" -filter "DP < 20.0" \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
        -filter-name "SOR_filter" -filter "SOR > 3.0" \
        &> /dev/null'%(gt_cmd,ref,outdir,iter,outdir,iter)
    os.system(cmd)

def mask_fasta(bt_cmd, ref, cov, bam, outfasta, outbed):
    cmd1 = '''%s genomecov \
            -ibam %s \
            -bga | \
            awk -v threshold=%s '$4<threshold' | \
            awk '{print $1 "\t" $2 "\t" $3}' \
            > %s'''%(bt_cmd, bam, cov, outbed)
    cmd2 = "%s maskfasta \
            -fi %s \
            -bed %s \
            -fo %s"%(bt_cmd, ref, outbed, outfasta)
    os.system(cmd1)
    os.system(cmd2)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='ConsIter produces an updated reference genome that more closely matches the sequenced viral population by iteratively mapping reads to the reference and generating a new consensus sequence.')
    parser.add_argument('-ref', dest = 'ref', type = str, required=True,  help = 'Path to reference genome')
    parser.add_argument('-n', dest = 'name', type = str, required=True,  help = 'Sample name')
    parser.add_argument('-o', dest = 'outdir', type = str, required=True,  default=".", help = 'Output directory. Default is the current directory')
    parser.add_argument('-1', dest= 'Left', type = str, required=True, help ='Left reads file in fastq format')
    parser.add_argument('-2', dest= 'Right', type = str, required=True, help ='Right reads file in fastq format')
    parser.add_argument('-i', dest= 'maxIter', type=int, default= 5, help ='Maximum number of iterations. Default = 5')
    parser.add_argument('-t', dest= 'Threads', type = int, default= 16, help ='Number of threads to use. Default = 16')
    parser.add_argument('-p', dest= 'ploidy', type = int, default= 1, help ='Sample ploidy used for GATK HaplotypeCaller. Default = 1')
    parser.add_argument('-c', dest= 'cov', type = int, default= 20, help ='Minimum coverage under which to mask bases in the final consensus. Default = 20')
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
    if shutil.which(bedtools_cmd) == None:
        print(sys.exit("Could not find Bedtools command. Adjust python script or load module"))
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
            print("Maximum number of iterations reached. Masking final consensus and terminating.")
            mask_fasta(bedtools_cmd, '%s/tmp/consensus.iter%s.fa'%(args.outdir,iteration-1), args.cov, '%s/tmp/iter%s.rmdup.bam'%(args.outdir,(iteration-1)), '%s/tmp/consensus.iter%s.masked.fa'%(args.outdir,iteration-1), '%s/tmp/consensus.iter%s.masked_sites.bed'%(args.outdir,iteration-1))
            rename_sname('%s/tmp/consensus.iter%s.masked.fa'%(args.outdir,iteration-1),args.ref,'%s/tmp/consensus.iter%s.masked.sample_name.fa'%(args.outdir,iteration-1),args.name)
            cmd0 = ("mv %s/tmp/consensus.iter%s.masked.sample_name.fa %s/%s.consensus.sample_name.fa"%(args.outdir,(iteration-1),args.outdir,args.name))
            os.system(cmd0)
            print("Updated reference genome with bases < %s coverage masked is in %s/%s.consensus.fa"%(args.cov, args.outdir,args.name))
            print("%s/%s.consensus.sample_name.fa includes sample name in fasta line"%(args.outdir,args.name))
            cmd1 = ("mv %s/tmp/consensus.iter%s.masked.fa %s/%s.consensus.fa"%(args.outdir,(iteration-1),args.outdir,args.name))
            os.system(cmd1)
            print("Duplicate-removed bam file is in %s/%s.bt2.rmdup.bam"%(args.outdir,args.name))
            cmd2 = ("mv %s/tmp/iter%s.rmdup.bam %s/%s.bt2.rmdup.bam"%(args.outdir,(iteration-1),args.outdir,args.name))
            os.system(cmd2)
            cmd3 = ("mv %s/tmp/iter%s.rmdup.bai %s/%s.bt2.rmdup.bai"%(args.outdir,(iteration-1),args.outdir,args.name))
            os.system(cmd3)
            print("Filtered VCF file is in %s/%s.filt.vcf"%(args.outdir,args.name))
            cmd4 = ("mv %s/tmp/iter%s.filt.vcf %s/%s.filt.vcf"%(args.outdir,(iteration-1),args.outdir,args.name))
            os.system(cmd4)
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
            call_variants(gatk_cmd, samtools_cmd, args.ref,iteration, args.xmx, args.outdir, args.ploidy)
            print("Filtering variants")
            filter_variants(gatk_cmd,args.ref,iteration,args.outdir)
            if args.noindel:
                print("Generating updated reference")
                select_snps(gatk_cmd, args.ref,iteration, args.outdir)
                consensus(gatk_cmd, args.ref,iteration,"%s/tmp/iter%s.snps.vcf"%(args.outdir,iteration), args.outdir)
                rename('%s/tmp/consensus.iter%s.tmpnames.fa'%(args.outdir,iteration),args.ref,'%s/tmp/consensus.iter%s.fa'%(args.outdir,iteration))
                iteration +=1
            else:
                print("Generating updated reference")
                consensus(gatk_cmd,args.ref,iteration,"%s/tmp/iter%s.filt.vcf"%(args.outdir,iteration), args.outdir)
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
                call_variants(gatk_cmd, samtools_cmd, "%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration, args.xmx, args.outdir, args.ploidy)
                print("Filtering variants")
                filter_variants(gatk_cmd,"%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration,args.outdir)
                if args.noindel:
                    print("Generating updated reference")
                    select_snps(gatk_cmd,"%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration, args.outdir)
                    consensus(gatk_cmd,"%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration,"%s/tmp/iter%s.snps.vcf"%(args.outdir,iteration), args.outdir)
                    rename('%s/tmp/consensus.iter%s.tmpnames.fa'%(args.outdir,iteration),args.ref,'%s/tmp/consensus.iter%s.fa'%(args.outdir,iteration))
                    iteration +=1
                else:
                    print("Generating updated reference")
                    consensus(gatk_cmd,"%s/tmp/consensus.iter%s.fa"%(args.outdir,(iteration-1)),iteration,"%s/tmp/iter%s.filt.vcf"%(args.outdir,iteration), args.outdir)
                    rename('%s/tmp/consensus.iter%s.tmpnames.fa'%(args.outdir,iteration),args.ref,'%s/tmp/consensus.iter%s.fa'%(args.outdir,iteration))
                    iteration +=1
            else:
                print("No improvement in alignment rate")
                print("Generating final consensus sequence and bam file")
                mask_fasta(bedtools_cmd, '%s/tmp/consensus.iter%s.fa'%(args.outdir,iteration-1), args.cov, '%s/tmp/iter%s.rmdup.bam'%(args.outdir,(iteration-1)), '%s/tmp/consensus.iter%s.masked.fa'%(args.outdir,iteration-1), '%s/tmp/consensus.iter%s.masked_sites.bed'%(args.outdir,iteration-1))
                rename_sname('%s/tmp/consensus.iter%s.masked.fa'%(args.outdir,iteration-1),args.ref,'%s/tmp/consensus.iter%s.masked.sample_name.fa'%(args.outdir,iteration-1),args.name)
                cmd0 = ("mv %s/tmp/consensus.iter%s.masked.sample_name.fa %s/%s.consensus.sample_name.fa"%(args.outdir,(iteration-1),args.outdir,args.name))
                os.system(cmd0)
                print("Updated reference genome with bases < %s coverage masked is in %s/%s.consensus.fa"%(args.cov, args.outdir,args.name))
                print("%s/%s.consensus.sample_name.fa includes sample name in fasta line"%(args.outdir,args.name))
                cmd1 = ("mv %s/tmp/consensus.iter%s.masked.fa %s/%s.consensus.fa"%(args.outdir,(iteration-1),args.outdir,args.name))
                os.system(cmd1)
                print("Duplicate-removed bam file is in %s/%s.bt2.rmdup.bam"%(args.outdir,args.name))
                cmd2 = ("mv %s/tmp/iter%s.rmdup.bam %s/%s.bt2.rmdup.bam"%(args.outdir,(iteration-1),args.outdir,args.name))
                os.system(cmd2)
                cmd3 = ("mv %s/tmp/iter%s.rmdup.bai %s/%s.bt2.rmdup.bai"%(args.outdir,(iteration-1),args.outdir,args.name))
                os.system(cmd3)
                print("Filtered VCF file is in %s/%s.filt.vcf"%(args.outdir,args.name))
                cmd4 = ("mv %s/tmp/iter%s.filt.vcf %s/%s.filt.vcf"%(args.outdir,(iteration-1),args.outdir,args.name))
                os.system(cmd4)
                if not args.keep:
                    os.system('rm -rf %s/tmp'%args.outdir)
                break
