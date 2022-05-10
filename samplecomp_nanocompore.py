#!/usr/bin/env python

# This code is used to test nanocompore. The path to data and the file names
# are specific to the test data.
# If you want to use it please take care of making it more generic or changing
# the paths and names to what suits your data.

import os
import glob
import subprocess as sp
import nanocompore as ncp
from nanocompore.SampCompDB import SampCompDB, jhelp



input_path = os.path.join("/scratch/local/test_nanocompore")
ref_gtf = '/data/repository/organisms/GRCm38_ensembl/Ensembl/release-91/genes.gtf'
ref_fasta = '/data/repository/organisms/GRCm38_ensembl/genome_fasta/genome.fa'
ref_trx = '/data/repository/organisms/GRCm38_ensembl/Ensembl/release-91/genes.bed'

## Prepare a fasta file for the reference transcriptome from ref_trx
cmd = "module load bedtools2 samtools; bedtools getfasta -fi "+ ref_fasta
cmd += " -s -split -name -bed "+ref_trx
cmd += " -fo  > "+ os.path.join(input_path, "reference_transcriptome.fa")
cmd += " ;samtools faidx "+os.path.join(input_path, "reference_transcriptome.fa")
print(cmd)
sp.check_output(cmd, shell = True)

# Make a bed file with strands in the names
cmd = "awk 'BEGIN{OFS=FS=\"\t\"}{$4=$4\"(\"$6\")\"; print}' "+ref_trx
cmd += " > "+ os.path.join(input_path,"reference_transcriptome_renamed.bed")
print(cmd)
sp.check_output(cmd, shell = True)

## map to transcriptome
for file in glob.glob(os.path.join(input_path, 'fastq', "*.fastq.gz")):
	name = os.path.basename(file).split(".fastq.gz")[0]
	cmd = "module load minimap2 samtools;"
	cmd += "minimap2 -x map-ont -t 10 "
	cmd += " -a "+os.path.join(input_path, "reference_transcriptome.fa")+" "
	cmd += file +" |  "
	cmd += "samtools view  -bh -t "+ os.path.join(input_path, "reference_transcriptome.fa.fai")
	cmd += " -F 2324 | samtools sort -@ 10 -o "
	cmd +=  os.path.join(input_path, "bam", name+".sorted.bam")+";"
	cmd += "samtools index "+os.path.join(input_path, "bam", name+".sorted.bam")
	print(cmd)
	sp.check_output(cmd, shell = True)

## Realign the raw signal-level data to the kmers of the reference using f5c
for file in glob.glob(os.path.join(input_path, 'fastq', "*.fastq.gz")):
    name = os.path.basename(file).split(".fastq.gz")[0]
    fast5_fail = os.path.join(input_path, "fast5", name, 'fast5_fail')
    fast5_pass = os.path.join(input_path, "fast5", name, 'fast5_pass')
    cmd = "/localenv/rabbani/anaconda/miniconda3/envs/ont/bin/f5c index -t 20 --iop 10 "
    cmd += "-d "+fast5_fail
    cmd += " -d "+fast5_pass
    cmd += " "+file
    print(cmd)
    sp.check_output(cmd, shell = True)
    bam_file = os.path.join(input_path, "bam", name+".sorted.bam")
    cmd = "/localenv/rabbani/anaconda/miniconda3/envs/ont/bin/f5c eventalign -t 50 "
    cmd += " -r "+file
    cmd += " -b "+bam_file+" -g "+os.path.join(input_path, "reference_transcriptome.fa")
    cmd += " --samples --print-read-names --scale-events --rna --disable-cuda=yes --min-mapq 0 | "
    cmd += "/localenv/rabbani/anaconda/miniconda3/envs/ont/bin/nanocompore eventalign_collapse -t 30 -o "
    cmd += os.path.join(input_path, name+'eventalign_collapse')
    print(cmd)
    sp.check_output(cmd, shell = True)

## sample comparison (Very very slow!)
cmd = "/localenv/rabbani/anaconda/miniconda3/envs/ont/bin/nanocompore sampcomp "
cmd += " --file_list1 " + os.path.join(input_path,"E12-EmxDOT1L-CTReventalign_collapse", "out_eventalign_collapse.tsv")
cmd += " --file_list2 " + os.path.join(input_path,"E12-EmxDOT1L-KOeventalign_collapse", "out_eventalign_collapse.tsv")
cmd += " --label1 CTR "
cmd += " --label2 KO "
cmd += " --fasta "+ os.path.join(input_path, "reference_transcriptome.fa")
cmd += " --outpath " + os.path.join(input_path, "nanocompore_KO_condition")
cmd += " --sequence_context 2 "
cmd += " --downsample_high_cov 5000 "
cmd += " --bed "+os.path.join(input_path,"reference_transcriptome_renamed.bed")
cmd += " --allow_warnings "
cmd += " --pvalue_thr 0.05 "
cmd += " --min_coverage 1 "
cmd += " --logit "
cmd += " --nthreads 20 "
cmd += " --log_level 'debug' "
cmd += " --progress "
print(cmd)
sp.check_output(cmd, shell = True)

print("sample comparison is done!")

## Load database
db = SampCompDB(
    db_fn = "/scratch/local/test_nanocompore/nanocompore_KO_condition/outSampComp.db",
    fasta_fn = ref_fasta)

## Print general metadata information
print(db)

## Prit list of references containing valid data
print(db.ref_id_list)

## save all info:
## Reload DB
db = SampCompDB (db_fn = "/scratch/local/test_nanocompore/nanocompore_KO_condition/outSampComp.db",
                 fasta_fn = ref_fasta)

Save report
db.save_report(output_fn = os.path.join(input_path,"report.tsv"))

## save sig positions:
## Reload DB
db = SampCompDB(db_fn = "/scratch/local/test_nanocompore/nanocompore_KO_condition/outSampComp.db",
                 fasta_fn = ref_fasta,
                 bed_fn= os.path.join(input_path,"reference_transcriptome_renamed.bed")
                 )
print(db)

## Print list of references containing valid data
print(db.results)

## Save report
db.save_to_bed(output_fn = os.path.join(input_path,"sig_positions.bed"),
			   pvalue_field= "GMM_logit_pvalue_context_2")
