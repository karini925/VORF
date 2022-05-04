#!/bin/bash

echo blast results from VORF:

#file name (input into the script is the original bam files)
cd /home/keren/DATA/RNAseq
bamfile=$1
id="${bamfile%%.*}"

#results from assembly are found here 
cd /home/keren/DATA/RNAseq/postQCrnaseq/aligned
cd $id

#for each fasta file run tblastn
for file in *.faa ; do
	echo "processing file=$file"
	tblastn -db nr -query $file  -out ${file}_output_blast.txt -max_target_seqs 5 -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -remote 
done

