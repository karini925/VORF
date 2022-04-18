# Viral Open Reading Frame assembly (VORF)
Efficient viral genome assembly using mixed sample RNA-seq reads  

## Input data

Currently, we are testing our tool on four types of viruses. 

## Processing of raw reads 

This includes aligning reads to human and keeping those that do not align. Further, Chimera reads are identified by looking at the rRNA reference to increase the number of viral reads available for assembly. 

## Viral ORF assembly 

Processed reads are split into k-mers. Each k-mer has six possible ORF phases to consider. We assemble contigs (transcripts) by initating a De Bruijn graph for each 'ATG' start codon that is found. 

