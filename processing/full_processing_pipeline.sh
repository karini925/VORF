#!/bin/bash

#[1] make human reference index file
cd ~/ANALYSIS/VORF/processing
nohup ./make_human_index.sh &

#[2] process BAM files that we downloaded from https://figshare.com/articles/dataset/VGEA_A_snakemake_pipeline_for_RNA_virus_genome_assembly_from_next_generation_sequencing_data/13009997/3
cd /home/keren/DATA/RNAseq
script=~/ANALYSIS/VORF/processing/fastq_qc_prep.sh  
nohup parallel -j10 --verbose "$script {}" ::: *.bam & 

#[3] map processed reads to the human reference and extract those that don't map
cd /home/keren/DATA/RNAseq
script=~/ANALYSIS/VORF/processing/align_human_reads.sh
nohup parallel -j2 --verbose "$script {}" ::: *.bam & 

#[4] run IVA assembly on clean reads 
cd /home/keren/DATA/RNAseq
script=~/ANALYSIS/VORF/assemblers/IVA_assembly.sh
nohup parallel -j2 --verbose "$script {}" ::: *.bam & 

#[5] run trinity on clean reads 
cd /home/keren/DATA/RNAseq
script=~/ANALYSIS/VORF/assemblers/trinity_assembly.sh
nohup parallel -j2 --verbose "$script {}" ::: *.bam & 

#[6] run VORF on all bam files
cd /home/keren/DATA/RNAseq
script=~/ANALYSIS/VORF/assemblers/run_vorf.sh
nohup parallel -j8 --verbose "$script {}" ::: *.bam > vorf.out & 

#[7] run blastn on IVA assembly results 
cd /home/keren/DATA/RNAseq
script=~/ANALYSIS/VORF/assemblers/run_blastn_iva.sh
nohup parallel -j8 --verbose "$script {}" ::: *.bam > program.out & 


