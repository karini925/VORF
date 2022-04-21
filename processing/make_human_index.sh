#!/bin/bash

human=/home/keren/DATA/GRCh38_latest_genomic.fna.gz

echo Indexing the human reference genome GRCh38_latest_genomic

bwa index ${human} #just do once 

