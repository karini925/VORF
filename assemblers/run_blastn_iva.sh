#!/bin/bash

echo blast results from IVA:

#file name (input into the script is the original bam files)
cd /home/keren/DATA/RNAseq
bamfile=$1
id="${bamfile%%.*}"

#results from assembly are found here 
cd /home/keren/DATA/RNAseq/postQCrnaseq/aligned
cd ${id}_iva

#for each fasta file run tblastn
blastn -db nr -query test.fasta -out ${id}_output_blast.txt -max_target_seqs 5 -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -remote 

