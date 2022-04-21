#!/bin/bash

#[1] make human reference index file
cd ~/ANALYSIS/VORF/processing
nohup ./make_human_index.sh &

#[2] process BAM files that we downloaded from https://figshare.com/articles/dataset/VGEA_A_snakemake_pipeline_for_RNA_virus_genome_assembly_from_next_generation_sequencing_data/13009997/3
cd /home/keren/DATA/RNAseq
script=~/ANALYSIS/VORF/processing/fastq_qc_prep.sh  
nohup parallel -j10 --verbose "$script {}" ::: *.bam & 

#[3] map processed reads to the human reference 
cd /home/keren/DATA/RNAseq
script=~/ANALYSIS/VORF/processing/align_human_reads.sh
nohup parallel -j10 --verbose "$script {}" ::: *.bam & 
