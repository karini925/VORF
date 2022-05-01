#!/bin/bash

echo Trinity assembly of de-hosted trimmed reads:

#file name (input into the script is the original bam files)
cd /home/keren/DATA/RNAseq
bamfile=$1
id="${bamfile%%.*}"

conda activate trinity-env

PATH_TO_TRINITY='/home/keren/miniconda3/envs/trinity-env/opt/trinity-2.8.5/Trinity'

#paired end reads are found here:
cd /home/keren/DATA/RNAseq/postQCrnaseq/aligned

r1=/home/keren/DATA/RNAseq/postQCrnaseq/aligned/${id}_1.fastq
r2=/home/keren/DATA/RNAseq/postQCrnaseq/aligned/${id}_2.fastq

output=/home/keren/DATA/RNAseq/postQCrnaseq/aligned/trinity
cd $output

mkdir $id 
cd $id

#Since the assembly process is stochastic multiple assemblies are created to be utilized for validation
${PATH_TO_TRINITY} --seqType fa  --left $r1  --right $r2 --max_memory 50G --CPU 8 --full_cleanup  --output ${id}_trinity_combination1
${PATH_TO_TRINITY} --seqType fa  --left $r1  --right $r2 --max_memory 30G --CPU 6  --output ${id}_trinity_combination2
${PATH_TO_TRINITY} --seqType fa  --left $r1  --right $r2 --max_memory 30G --CPU 6 --full_cleanup --output ${id}_trinity_combination3
${PATH_TO_TRINITY} --seqType fa  --left $r1  --right $r2 --max_memory 30G --CPU 6 --full_cleanup --output ${id}_trinity_combination4
${PATH_TO_TRINITY} --seqType fa  --left $r1  --right $r2 --max_memory 30G --CPU 6 --full_cleanup --output ${id}_trinity_combination5
