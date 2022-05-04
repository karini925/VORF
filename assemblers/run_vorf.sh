#!/bin/bash

echo VORF assembly of de-hosted trimmed reads:

#file name (input into the script is the original bam files)
cd /home/keren/DATA/RNAseq
bamfile=$1
id="${bamfile%%.*}"

#paired end reads are found here:
cd /home/keren/DATA/RNAseq/postQCrnaseq/aligned

r1=${id}_1.fastq
r2=${id}_2.fastq

script=/home/keren/ANALYSIS/VORF/main.py

python $script $r1 $r2 

echo Done running VORF, now time to run tblastn to see what we found 

#cd $id

#for file in *.faa ; do
#	echo "processing file=$file"
#	tblastn -db nr -query $file  -out ${file}_output_blast.txt -max_target_seqs 5 -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -remote 
#done
