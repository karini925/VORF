#!/bin/bash

#align_human_reads.sh 

#make sure we are in base conda environment
#need bedtools, samtools and fastp to work

#--------------------------------------------------
#1. map fastq files to human reference
#--------------------------------------------------

human_ref=/home/keren/DATA/human_ref/GRCh38_latest_genomic.fna.gz
output=/home/keren/DATA/RNAseq/postQCrnaseq/aligned

cd /home/keren/DATA/RNAseq

#fastq files 

bamfile=$1
id="${bamfile%%.*}"

cd /home/keren/DATA/RNAseq/postQCrnaseq

r1=$id.r1_trimmed.fq
r2=$id.r2_trimmed.fq

echo Mapping reads to the human genome to remove human contaminants

bwa mem -t 8 ${human_ref} ${r1} ${r2} > $output/${id}.sam

#--------------------------------------------------
#2. extract reads that didn't map to human 
#--------------------------------------------------

echo Extracting unmapped reads from bam file

cd $output 

samtools view -b -f12 ${id}.sam > ${id}.unmapped.bam 

#--------------------------------------------------
#3. convert reads back to forward and reverse strand
#--------------------------------------------------

echo Converting BAM file into fastq files of forward and reverse reads

samtools flagstat ${id}.unmapped.bam > ${id}.flagstat
samtools fastq -1 ${id}_1.fastq -2 ${id}_2.fastq ${id}.unmapped.bam 

echo reads ready for assemebly 