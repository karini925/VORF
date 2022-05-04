#!/bin/bash

echo blast results from VORF:

#file name (input into the script is the original bam files)
cd /home/keren/DATA/RNAseq
bamfile=$1
id="${bamfile%%.*}"

#results from assembly are found here 
cd /home/keren/DATA/RNAseq/postQCrnaseq/aligned
cd ${id}_iva

#for each fasta file run tblastn
blastn -db nr -query contigs.fasta -out ${id}_output_blast.txt -outfmt 6 -remote 

