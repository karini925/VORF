#!/bin/bash

#fastq_interleave.sh 

#interleave paired-end fastq files to obtain one fastq file which is often what is required 
#by dependency tools 

#Here, we used the code from https://github.com/ekg/interleave-fastq/blob/master/interleave-fastq 
#to achieve this task

#For example 

#in this case, we are looking at HIV sample from https://figshare.com/articles/dataset/VGEA_A_snakemake_pipeline_for_RNA_virus_genome_assembly_from_next_generation_sequencing_data/13009997/3
#bam file was first split into forward/backward paired end individual fsatq files 

python2 interleave_fastq.py ~/DATA/RNAseq/ERR3953893.end1.fq ~/DATA/RNAseq/ERR3953893.end2.fq > ERR3953893.fq
