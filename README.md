# Viral Open Reading Frame assembly (VORF)

Pandemics due to novel pathogens are an imminent threat, yet few methods for their speedy characterization exist. Current assemblies fall short of identifying novel sequences from environmental samples in an organism-agnostic manner.

VORF is an RNA-seq-based open reading frame assembler capable of detecting amino acid sequences from clinical samples without a reference genome. 

## Datasets 

We applied our tool as well as the assembler called [IVA](https://github.com/sanger-pathogens/iva) to nine samples across three viruses: COVID-19, HIV and Lassa. All these are RNA viruses and were downloaded from the [VGEA manuscript](https://figshare.com/articles/dataset/VGEA_A_snakemake_pipeline_for_RNA_virus_genome_assembly_from_next_generation_sequencing_data/13009997/3). An example dataset from COVID-19 is available in our repository [here](data). A summary of all nine samples is available [here](https://docs.google.com/spreadsheets/d/1zvgPzrfHkJR6LYx6D1xI5kHnIWBcRBKIrzlMGILYPk4/edit?usp=sharing).  

## Processing of raw reads 

The data we downloaded from VGEA contained BAM files with unaligned reads. We first split BAM files into forward and reverse fastq files. We then ran some QC on the reads, aligned them to the human reference and kept only those that did not align. All this processing is summarized in the following script: [full_processing_pipeline.sh](processing/full_processing_pipeline.sh). All processing scripts were run one at a time on all BAM files in parallel. 

## Assembly algorithm 

The assembly algorithm was written in python and can be run using the [main.py file](main.py). The methods used for assembly are stored in the [utils.py file](utils.py).



To run VORF, the following packages need to be installed (can easily be done with conda): ```scikit-learn```, ```biopython``` and ```levenshtein``` with a python 3 environment. 

```
#enter directory with fastq input file 
r1=/home/keren/ANALYSIS/VORF/data/CV167_1.fastq
r2=/home/keren/ANALYSIS/VORF/data/CV167_2.fastq

#run VORF 
script=/home/keren/ANALYSIS/VORF/main.py
python $script $r1 $r2 
```
Below, an example of an output is shown for the above samples. 

```
>0 <unknown description>
MVFLVCHGLVSCNNPLAITVLYPHQCLARGTT
>10 <unknown description>
MARGLLQLTNPWQTKNTIKWFLGRQTAGANVRCQEGNNPDRQLRSQMID
>11 <unknown description>
MVWILLLSDMDLSTHALTPWILLNVHSEFVWT
>12 <unknown description>
MASVTSKNTTKARKRSTLLTINVPVSSETNEYISSYSSACAYKGTLVVVVGSS
>13 <unknown description>
MHSTATPRIPPTSSTLKPASINGNLGVKLLDFTADLTGRLRTL
>14 <unknown description>
MATSSFITVLTKKNLAVSHWYFAHMRDKDTKCLTTTTVLNAFEFCYQLNHNMTMRCSSSI
```
