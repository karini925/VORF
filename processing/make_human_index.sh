#!/bin/bash

human=/home/keren/DATA/human_ref/GRCh38_latest_genomic.fna.gz

echo Indexing the human reference genome GRCh38_latest_genomic

bwa index ${human} #just do once 

#nohup ./annotate_vcfs.sh & #(will keep running in background)

