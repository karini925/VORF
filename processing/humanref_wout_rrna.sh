#!/bin/bash

#humanref_wout_rrna.sh

#If you know the 16S sequences (or those of a close relative), you can align them to the genome to produce a sam file. For example:
mapPacBio=/home/keren/ANALYSIS/bbmap/mapPacBio.sh
bbmask=/home/keren/ANALYSIS/bbmap/bbmask.sh

ribo=/home/keren/DATA/rRNA_genbank_U13369_1.fasta
human=/home/keren/DATA/GRCh38_latest_genomic.fna.gz

$mapPacBio ref=$human in=$ribo out=mapped.sam ambig=all maxindel=20 -Xmx 30G

#Then you can run BBMask:

$bbmask in=$human out=masked.fasta sam=mapped.sam masklowentropy=false maskrepeats=false

#torun nohup ./annotate_vcfs.sh & #(will keep running in background)
