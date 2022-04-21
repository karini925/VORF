#!/bin/bash

#fastq_qc_prep.sh 

#make sure we are in base conda environment
#need bedtools, samtools and fastp to work

#--------------------------------------------------
#1. split BAM into fastq files
#--------------------------------------------------

output=/home/keren/DATA/RNAseq/postQCrnaseq

#check BAM file
bam_reads=$1

bam_name="${bam_reads%%.*}"

#sort 
samtools sort -n -o $output/$bam_name.sort.bam $bam_reads

cd $output

#check bam file header 
echo check bam file header 
samtools view -H $bam_name.sort.bam

#stats on bam
echo the number of reads in the sample is:
samtools flagstat $bam_name.sort.bam

#split bam 
echo splitting bam into fastq files
bedtools bamtofastq -i $bam_name.sort.bam \
                      -fq $bam_name.end1.fq \
                      -fq2 $bam_name.end2.fq

#read processing 
echo Trimming and performing read-level QC using FASTP
fastp --in1 $bam_name.end1.fq --in2 $bam_name.end2.fq --out1 $bam_name.r1_trimmed.fq --out2 $bam_name.r2_trimmed.fq --report_title $bam_name --json $bam_name.fastp.json --html $bam_name.fastp.html

echo removing intermediate files 
rm $bam_name.sort.bam
rm $bam_name.end1.fq
rm $bam_name.end2.fq

echo now reads ready for mapping!

#to run on all bam files in directory 
#parallel -j12 --verbose "plink --bfile {.} --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump $pgc_summ --clump-snp-field SNP --clump-field P --allow-no-sex --out {.}" ::: *.bam 
