#!/usr/bin/env python3

"""ConsIter produces an updated reference genome that more
closely matches the sequenced viral population by iteratively
mapping reads to the reference and generating a new consensus sequence.
"""

import argparse
import subprocess
import shutil
import re
import os
import sys

# set base commands for software on your system
######################################################
BOWTIE2_CMD = "bowtie2"
SAMTOOLS_CMD = "samtools"
GATK_CMD = "gatk"
PICARD_CMD = "picard"
BEDTOOLS_CMD = "bedtools"
######################################################

def run_bowtie2_idx(bt2_cmd,ref_fasta):
    """Index the reference with bowtie2"""
    cmd = (
        f"{bt2_cmd}-build -q {ref_fasta} {ref_fasta} "
        f"&> /dev/null"
        )
    os.system(cmd)

def run_bowtie2(bt2_cmd, smt_cmd, ref, i, threads, name, left, right, outdir):
    """run bowtie2 and collect alignment stats"""
    cmd1 = (
            f"{bt2_cmd} -p {threads} "
            f"--no-unal "
            f"--rg-id {name} "
            f"--rg SM:{name} "
            f"--rg LB:1 "
            f"--rg PU:1 "
            f"--rg PL:illumina "
            f"--sensitive-local "
            f"-x {ref} "
            f"-1 {left} "
            f"-2 {right} "
            f"-S {outdir}/tmp/iter{i}.sam "
            f"2> {outdir}/tmp/iter{i}.btstats.txt"
        )
    #convert to sorted bam
    cmd2 = (
        f"{smt_cmd} view -bh "
        f"-f 0x2 "
        f"-q 20 "
        f"{outdir}/tmp/iter{i}.sam | "
        f"samtools sort -@ {threads} "
        f"-o {outdir}/tmp/iter{i}.bam "
        f"- &> /dev/null"
        )
    os.system(cmd1)
    os.system(cmd2)

def align_rate(i, outdir):
    """Collect alignment rate from bowtie2 output"""
    with open(f"{outdir}/tmp/iter{i}.btstats.txt", "r", encoding="utf-8") as fh:
        bt2str = fh.read()
        m = re.search(r'(\d+\.\d+)\% overall alignment rate', bt2str)
        alnrt = m.group(1)
        return float(alnrt)

def rmdup(pic_cmd, i, outdir):
    """Remove duplicates"""
    cmd = (
        f"{pic_cmd} MarkDuplicates "
        f"CREATE_INDEX=true "
        f"USE_JDK_DEFLATER=true "
        f"USE_JDK_INFLATER=true "
        f"M={outdir}/tmp/iter{i}.rmdup_metrics.txt "
        f"I={outdir}/tmp/iter{i}.bam "
        f"O={outdir}/tmp/iter{i}.rmdup.bam "
        f"REMOVE_DUPLICATES=true "
        f"VALIDATION_STRINGENCY=LENIENT "
        f"&> /dev/null"
        )
    os.system(cmd)


def call_variants(gt_cmd, smt_cmd, ref, i, xmx, outdir, ploidy):
    """Call variants with GATK"""
    refname = '.'.join(ref.split('.')[:-1])
    cmd1 = f"{smt_cmd} faidx {ref}"
    cmd2 = f"{smt_cmd} dict {ref} > {refname}.dict"
    cmd3 = (
            f"{gt_cmd} --java-options '-Xmx{xmx}g' HaplotypeCaller "
            f"--use-jdk-deflater --use-jdk-inflater "
            f"-R {ref} "
            f"-I {outdir}/tmp/iter{i}.rmdup.bam "
            f"-O {outdir}/tmp/iter{i}.vcf "
            f"--min-base-quality-score 20 "
            f"-ploidy {ploidy} "
            f"&> /dev/null"
            )
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)

def gen_consensus(gt_cmd,ref,i,vcf, outdir):
    """Generate the consensus sequence with GATK"""
    cmd1 = f"{gt_cmd} IndexFeatureFile -I {vcf} &> /dev/null"
    cmd2 = (
        f"{gt_cmd} FastaAlternateReferenceMaker "
        f"-R {ref} "
        f"-O {outdir}/tmp/consensus.iter{i}.tmpnames.fa "
        f"-V {vcf} "
        f"&> /dev/null"
       )
    os.system(cmd1)
    os.system(cmd2)

def select_snps(gt_cmd, ref,i, outdir):
    """Select only SNPs from the GATK VCF file"""
    cmd = (
        f"{gt_cmd} SelectVariants "
        f"--use-jdk-deflater --use-jdk-inflater "
        f"-R {ref} "
        f"-V {outdir}/tmp/iter{i}.filt.vcf "
        f"-select-type SNP "
        f"-O {outdir}/tmp/iter{i}.snps.vcf &> /dev/null"
        )
    os.system(cmd)

def rename_sname(consensus,ref,outfasta,sample):
    """Rename the fasta headers with the sample name and reference name"""
    with (
    open(consensus, 'r', encoding="utf-8") as f,
    open(outfasta,'w', encoding="utf-8") as new,
    open('.'.join(ref.split('.')[:-1])+'.dict','r', encoding="utf-8") as refdict
    ):
        fasta = f.read()
        next(refdict)
        count=1
        for line in refdict:
            name = line.split('\t')[1].replace("SN:","")
            fasta = fasta.replace(">", f">{sample}|{name}")
            count+=1
        new.write(fasta)

def rename(consensus,ref,outfasta):
    """Rename the sample readers with only the reference seq name"""
    with (
    open(consensus, 'r', encoding="utf-8") as f,
    open(outfasta,'w', encoding="utf-8") as new,
    open('.'.join(ref.split('.')[:-1])+'.dict','r', encoding="utf-8") as refdict
    ):
        fasta = f.read()
        next(refdict)
        count=1
        for line in refdict:
            name = line.split('\t')[1].replace("SN:","")
            fasta = fasta.replace(f">{count}", f">{name}")
            count+=1
        new.write(fasta)

def filter_variants(gt_cmd, ref, i, outdir):
    """Filter variants based on GATK's suggestions"""
    cmd = (
        f'{gt_cmd} VariantFiltration '
        f'-R {ref} '
        f'-V {outdir}/tmp/iter{i}.vcf '
        f'-O {outdir}/tmp/iter{i}.filt.vcf '
        f'-filter-name "DP_filter" -filter "DP < 20.0" '
        f'-filter-name "QD_filter" -filter "QD < 2.0" '
        f'-filter-name "FS_filter" -filter "FS > 60.0" '
        f'-filter-name "MQ_filter" -filter "MQ < 40.0" '
        f'-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" '
        f'-filter-name "SOR_filter" -filter "SOR > 3.0" '
        f'&> /dev/null'
        )
    os.system(cmd)

def mask_fasta(bt_cmd, ref, cov, bam, outfasta, outbed):
    """Mask sequences with low coverage"""
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
    if shutil.which(BOWTIE2_CMD) is None:
        print(sys.exit("Could not find bowtie2 command. Adjust python script or load module"))
    if shutil.which(SAMTOOLS_CMD) is None:
        print(sys.exit("Could not find Samtools command. Adjust python script or load module"))
    if shutil.which(BEDTOOLS_CMD) is None:
        print(sys.exit("Could not find Bedtools command. Adjust python script or load module"))
    gatk_stat = subprocess.getstatusoutput(GATK_CMD)
    gatk_m = re.search("Usage", str(gatk_stat[1]))
    if gatk_m is None:
        print(sys.exit("Could not find GATK command. Adjust python script or load module"))
    picard_stat = subprocess.getstatusoutput(PICARD_CMD)
    picard_m = re.search("PicardCommandLine", str(picard_stat[1]))
    if picard_m is None:
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
    ITERATION = 0

    while ITERATION <= args.maxIter:
        #reached last ITERATION
        if ITERATION == args.maxIter:
            print(f"Iteration {ITERATION}")
            print("Maximum number of ITERATIONs reached. Masking final consensus and terminating.")
            mask_fasta(BEDTOOLS_CMD,
                        f'{args.outdir}/tmp/consensus.iter{ITERATION-1}.fa',
                        args.cov,
                        f'{args.outdir}/tmp/iter{ITERATION-1}.rmdup.bam',
                        f'{args.outdir}/tmp/consensus.iter{ITERATION-1}.masked.fa',
                        f'{args.outdir}/tmp/consensus.iter{ITERATION-1}.masked_sites.bed')
            rename_sname(f'{args.outdir}/tmp/consensus.iter{ITERATION-1}.masked.fa',
                        args.ref,
                        f'{args.outdir}/tmp/consensus.iter{ITERATION-1}.masked.sample_name.fa',
                        args.name)
            command0 = (
                   f"mv {args.outdir}/tmp/consensus.iter{ITERATION-1}.masked.sample_name.fa "
                   f"{args.outdir}/{args.name}.consensus.sample_name.fa"
                   )
            os.system(command0)
            print(f"Updated reference genome with bases < {args.cov} coverage masked is in {args.outdir}/{args.name}.consensus.fa")
            print(f"{args.outdir}/{args.name}.consensus.sample_name.fa includes sample name in fasta line")
            command1 = f"mv {args.outdir}/tmp/consensus.iter{ITERATION-1}.masked.fa {args.outdir}/{args.name}.consensus.fa"
            os.system(command1)
            print(f"Duplicate-removed bam file is in {args.outdir}/{args.name}.bt2.rmdup.bam")
            command2 = f"mv {args.outdir}/tmp/iter{ITERATION-1}.rmdup.bam {args.outdir}/{args.name}.bt2.rmdup.bam"
            os.system(command2)
            command3 = f"mv {args.outdir}/tmp/iter{ITERATION-1}.rmdup.bai {args.outdir}/{args.name}.bt2.rmdup.bai"
            os.system(command3)
            if not args.keep:
                os.system(f'rm -rf {args.outdir}/tmp')
            ITERATION +=1
        # in first ITERATION, map to original reference genome
        elif ITERATION == 0:
            print(f"Iteration {ITERATION}")
            print("Indexing reference")
            run_bowtie2_idx(BOWTIE2_CMD,args.ref)
            print("Mapping reads")
            run_bowtie2(BOWTIE2_CMD, SAMTOOLS_CMD, args.ref, ITERATION, args.Threads, args.name, args.Left, args.Right, args.outdir)
            alnrate = align_rate(ITERATION, args.outdir)
            print("Removing duplicates")
            rmdup(PICARD_CMD,ITERATION, args.outdir)
            print(f"Iteration {ITERATION} alignment rate: {alnrate}")
            print("Calling variants")
            call_variants(GATK_CMD, SAMTOOLS_CMD, args.ref,ITERATION, args.xmx, args.outdir, args.ploidy)
            print("Filtering variants")
            filter_variants(GATK_CMD,args.ref,ITERATION,args.outdir)
            if args.noindel:
                print("Generating updated reference")
                select_snps(GATK_CMD, args.ref,ITERATION, args.outdir)
                gen_consensus(GATK_CMD, args.ref,ITERATION,f"{args.outdir}/tmp/iter{ITERATION}.snps.vcf", args.outdir)
                rename(f'{args.outdir}/tmp/consensus.iter{ITERATION}.tmpnames.fa',args.ref,f'{args.outdir}/tmp/consensus.iter{ITERATION}.fa')
                ITERATION +=1
            else:
                print("Generating updated reference")
                gen_consensus(GATK_CMD,args.ref,ITERATION,f"{args.outdir}/tmp/iter{ITERATION}.filt.vcf", args.outdir)
                rename(f'{args.outdir}/tmp/consensus.iter{ITERATION}.tmpnames.fa',args.ref,f'{args.outdir}/tmp/consensus.iter{ITERATION}.fa')
                ITERATION +=1
        else:
            #map to updated reference and check if alignment rate is better than last ITERATION
            print(f"Iteration {ITERATION}")
            print("Indexing updated reference")
            run_bowtie2_idx(BOWTIE2_CMD,f"{args.outdir}/tmp/consensus.iter{ITERATION-1}.fa")
            print("Mapping reads")
            run_bowtie2(BOWTIE2_CMD,
                        SAMTOOLS_CMD,
                        f"{args.outdir}/tmp/consensus.iter{ITERATION-1}.fa",
                        ITERATION, args.Threads, args.name, args.Left, args.Right, args.outdir)
            alnrate_last = align_rate(ITERATION-1, args.outdir)
            alnrate = align_rate(ITERATION, args.outdir)
            print(f"Iteration {ITERATION-1} alignment rate: {alnrate_last}")
            print(f"Iteration {ITERATION} alignment rate: {alnrate}")
            if alnrate > alnrate_last:
                #alignment rate is better than last, make updated reference and continue
                print(f"Iteration {ITERATION} alignment rate better than previous iteration")
                print("Removing duplicates")
                rmdup(PICARD_CMD,ITERATION, args.outdir)
                print("Calling variants")
                call_variants(GATK_CMD, SAMTOOLS_CMD,
                            f"{args.outdir}/tmp/consensus.iter{ITERATION-1}.fa",
                            ITERATION, args.xmx, args.outdir, args.ploidy)
                print("Filtering variants")
                filter_variants(GATK_CMD,f"{args.outdir}/tmp/consensus.iter{ITERATION-1}.fa",ITERATION,args.outdir)
                if args.noindel:
                    print("Generating updated reference")
                    select_snps(GATK_CMD,f"{args.outdir}/tmp/consensus.iter{ITERATION-1}.fa",ITERATION, args.outdir)
                    gen_consensus(GATK_CMD,
                              f"{args.outdir}/tmp/consensus.iter{ITERATION-1}.fa",
                              ITERATION,
                              f"{args.outdir}/tmp/iter{ITERATION}.snps.vcf", args.outdir)
                    rename(f'{args.outdir}/tmp/consensus.iter{ITERATION}.tmpnames.fa',
                            args.ref,
                            f'{args.outdir}/tmp/consensus.iter{ITERATION}.fa')
                    ITERATION +=1
                else:
                    print("Generating updated reference")
                    gen_consensus(GATK_CMD,
                              f"{args.outdir}/tmp/consensus.iter{ITERATION-1}.fa",
                              ITERATION,
                              f"{args.outdir}/tmp/iter{ITERATION}.filt.vcf", args.outdir)
                    rename(f'{args.outdir}/tmp/consensus.iter{ITERATION}.tmpnames.fa',
                            args.ref,
                            f'{args.outdir}/tmp/consensus.iter{ITERATION}.fa')
                    ITERATION +=1
            else:
                #alignment rate is not better than last, finish
                print("No improvement in alignment rate")
                print("Generating final consensus sequence and bam file")
                mask_fasta(BEDTOOLS_CMD,
                           f'{args.outdir}/tmp/consensus.iter{ITERATION-1}.fa',
                           args.cov,
                           f'{args.outdir}/tmp/iter{ITERATION-1}.rmdup.bam',
                           f'{args.outdir}/tmp/consensus.iter{ITERATION-1}.masked.fa',
                           f'{args.outdir}/tmp/consensus.iter{ITERATION-1}.masked_sites.bed')
                rename_sname(f'{args.outdir}/tmp/consensus.iter{ITERATION-1}.masked.fa',
                            args.ref,
                            f'{args.outdir}/tmp/consensus.iter{ITERATION-1}.masked.sample_name.fa',args.name)
                command0 = f"mv {args.outdir}/tmp/consensus.iter{ITERATION-1}.masked.sample_name.fa {args.outdir}/{args.name}.consensus.sample_name.fa"
                os.system(command0)
                print(f"Updated reference genome with bases < {args.cov} coverage masked is in {args.outdir}/{args.name}.consensus.fa")
                print(f"{args.outdir}/{args.name}.consensus.sample_name.fa includes sample name in fasta line")
                command1 = f"mv {args.outdir}/tmp/consensus.iter{ITERATION-1}.masked.fa {args.outdir}/{args.name}.consensus.fa"
                os.system(command1)
                print(f"Duplicate-removed bam file is in {args.outdir}/{args.name}.bt2.rmdup.bam")
                command2 = f"mv {args.outdir}/tmp/iter{ITERATION-1}.rmdup.bam {args.outdir}/{args.name}.bt2.rmdup.bam"
                os.system(command2)
                command3 = f"mv {args.outdir}/tmp/iter{ITERATION-1}.rmdup.bai {args.outdir}/{args.name}.bt2.rmdup.bai"
                os.system(command3)
                if not args.keep:
                    os.system(f'rm -rf {args.outdir}/tmp')
                break
