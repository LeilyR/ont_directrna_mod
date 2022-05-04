#!/usr/bin/env python


## Note: Tombo provides tools and a framework to identify consistent shifts in raw nanopore signal.
##  The interpretation of these shifts is left to the user.

import os
import glob
import shutil
import subprocess as sp
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Process Fast5 files to find modification on RNA data using Tombo.')
    parser.add_argument('-f5','--fast5', dest='fast5',
                        type = str, help='path to fast5')
    parser.add_argument('-fq','--fastq', dest='fastq',
                        type = str, help='path to fastq') # if m2s any random path can be used for the fastq
    parser.add_argument('-o','--ouput', dest='output',
                        type = str, help='path to output')
    parser.add_argument('-s','--samples', dest='samples',
                        type = str, nargs = "+",
                        help='sample names')
    parser.add_argument('-rf','--ref_fasta', dest='ref_fasta',
                        type = str, help='referecne fasta')
    parser.add_argument('-rt','--ref_trx', dest='ref_trx',
                        type = str, help='reference transcripts bed file')
    parser.add_argument('-m2s','--multi2single', action = 'store_true',
                        dest='multi2single',
                        help='make single fast5')

    return parser.parse_args()


def main():
    args = parse_args()
    fast5 = os.path.abspath(args.fast5)
    fastq = os.path.abspath(args.fastq)
    # ref_fasta = '/data/repository/organisms/GRCm38_ensembl/genome_fasta/genome.fa'
    # ref_trx = '/data/repository/organisms/GRCm38_ensembl/Ensembl/release-91/genes.bed'
    ref_fasta = args.ref_fasta
    ref_trx = args.ref_trx
    print(args.samples)
    #1. multi_to_single_fast5
    if args.multi2single:
        print("here!")
        os.makedirs(os.path.join(args.output, "singlefastq"), exist_ok = True)
        for sample in args.samples:
            print(sample)
            cmd = "module load ont-fast5-api;"
            cmd += "multi_to_single_fast5 -i "+os.path.join(fast5, sample)
            cmd += " --recursive -t 30 -s "
            cmd += os.path.join(args.output, "singlefast5", sample)
            print(cmd)
            sp.check_output(cmd, shell = True)
            for dir in glob.glob(os.path.join(args.output, "singlefast5", sample, "*")):
                print(dir)
                if os.path.isdir(dir):
                    for i in glob.glob(os.path.join(dir,'*.fast5')):
                        shutil.move(i, os.path.join(args.output, "singlefast5", sample))
                    os.rmdir(dir)
            # 1.2. basecall using guppy with fast5 and fastq output (but with single fast5´s as input)
            cmd = "/data/processing1/leily/guppy-5.0.7/ont-guppy/bin/guppy_basecaller "
            cmd += " -i "+os.path.join(args.output, "singlefast5", sample)
            cmd += " -s "+os.path.join(args.output, "singlefastq", sample)
            cmd += " -r  --flowcell FLO-MIN106 --kit SQK-RNA002 "
            cmd += " --device cuda:0 --num_callers 50 "
            cmd += " --reverse_sequence true --u_substitution true --trim_strategy rna "
            print(cmd)
            sp.check_output(cmd, shell = True)
        fastq = os.path.join(args.output, "singlefastq")
        fast5 = os.path.join(args.output, "singlefast5")
    # 2. tombo preprocess
    for sample in args.samples:
        cmd = "/localenv/rabbani/anaconda/miniconda3/envs/ont/bin/tombo"
        cmd += " preprocess annotate_raw_with_fastqs --fast5-basedir "+os.path.join(fast5, sample)
        cmd += " --fastq-filenames "+os.path.join(fastq, sample, "pass", "*.fastq")
        cmd += " "+os.path.join(fastq, sample, "fail", "*.fastq")
        cmd += " --process 20 --overwrite --basecall-group Basecall_1D_000 --basecall-subgroup BaseCalled_template "
        print(cmd)
        sp.check_output(cmd, shell = True)
    #3. make a transcripts ref
    cmd = "module load bedtools2;"
    cmd += "bedtools getfasta -s -split -name "
    cmd += " -fi "+args.ref_fasta
    cmd += " -bed "+args.ref_trx+">"
    cmd += os.path.join(args.output , "transcripts.fa")
    cmd += ";sed -i \"s/(+)//; s/(-)//\" "
    cmd += os.path.join(args.output , "transcripts.fa")
    print(cmd)
    sp.check_call(cmd, shell=True)
    #4. tombo re-squiggle
    # # Tombo does not perform spliced mapping. Thus a transcriptime reference must be passed to the re-squiggle command for RNA samples.
    for sample in args.samples:
        cmd = "/localenv/rabbani/anaconda/miniconda3/envs/ont/bin/tombo "
        cmd += " resquiggle "
        cmd += os.path.join(fast5, sample)
        cmd += " "+os.path.join(args.output , "transcripts.fa")
        cmd += "  --rna --processes 15 --overwrite "
        print(cmd)
        sp.check_output(cmd, shell = True)
    #5. tombo detect_modifications by smaple comparison
    # default of min test read is 50 but developers suggested a higher value for a more reliable result
    cmd = "/localenv/rabbani/anaconda/miniconda3/envs/ont/bin/tombo "
    cmd +=" detect_modifications level_sample_compare --fast5-basedirs "
    cmd += os.path.join(fast5, "E12-EmxDOT1L-CTR")
    cmd += " --alternate-fast5-basedirs "
    cmd += os.path.join(fast5, "E12-EmxDOT1L-KO")
    cmd += " --processes 10 --minimum-test-reads 100 --store-p-value --statistics-file-basename "
    cmd += os.path.join(args.output, "CTR_vs_KO")
    print(cmd)
    sp.check_output(cmd, shell = True)
    # 6.browser file
    cmd = "/localenv/rabbani/anaconda/miniconda3/envs/ont/bin/tombo text_output browser_files "
    cmd += " --fast5-basedirs "+os.path.join(fast5, "E12-EmxDOT1L-CTR")
    cmd += " --control-fast5-basedirs "+os.path.join(fast5, "E12-EmxDOT1L-KO")
    cmd += " --statistics-filename "+os.path.join(args.output, "CTR_vs_KO.tombo.stats")
    print(cmd)
    # 7. plot
    # box plot fails for some ggplot version conflict
    sp.check_output(cmd, shell = True)
if __name__ == "__main__":
    main()
