#!/bin/bash

#humanref_wout_rrna.sh
#mask the rRNA sequence in the human genome 

#where to save masked reference 
output=/home/keren/DATA/RNAseq/human

#If you know the 16S sequences (or those of a close relative), you can align them to the genome to produce a sam file. For example:
mapPacBio=/home/keren/ANALYSIS/bbmap/mapPacBio.sh
bbmask=/home/keren/ANALYSIS/bbmap/bbmask.sh

ribo=/home/keren/DATA/rRNA_genbank_U13369_1.fasta
human=/home/keren/DATA/GRCh38_latest_genomic.fna.gz

$mapPacBio ref=$human in=$ribo out=$output/mapped.sam ambig=all maxindel=20 nodisk usemodulo

#Then you can run BBMask:

$bbmask in=$human out=$output/human_rRNA_masked.fasta sam=$output/mapped.sam masklowentropy=false maskrepeats=false

#remove the mapped.sam file as we don't need it
rm $output/mapped.sam

echo DONE

#torun 
#nohup ./humanref_wout_rrna.sh & #(will keep running in background)
