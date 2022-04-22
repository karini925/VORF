#!/bin/bash

echo IVA assembly of de-hosted trimmed reads:

#file name (input into the script is the original bam files)
cd /home/keren/DATA/RNAseq
bamfile=$1
id="${bamfile%%.*}"

#paired end reads are found here:
cd /home/keren/DATA/RNAseq/postQCrnaseq/aligned

r1=${id}_1.fastq
r2=${id}_2.fastq

iva --reads_fwd $r1 --reads_rev $r2 --threads 6 ${id}_iva
