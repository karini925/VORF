# Viral Open Reading Frame assembly (VORF)

Pandemics due to novel pathogens are an imminent threat, yet few methods for their speedy characterization exist. Current assemblies fall short of identifying novel sequences from environmental samples in an organism-agnostic manner.

VORF is an RNA-seq-based open reading frame assembler capable of detecting amino acid sequences from clinical samples without a reference genome. 

## Datasets 

We applied our tool as well the assembler [IVA](https://github.com/sanger-pathogens/iva) to nine samples across three viruses: COVID-19, HIV and Lassa. All these are RNA viruses. An example dataset is available in our repository here. 

## Processing of raw reads 

This includes running QC on reads, aligning them to the human reference and keeping only those that do not align. All this processing is summarized in the following script: [full_processing_pipeline.sh](processing/full_processing_pipeline.sh)

## Assembly algorithm 

Noga fill in a little bit of detail about tool 

To run VORF, the following packages need to be installed (can easily be done with conda): ```scikit-learn```, ```biopython``` and ```levenshtein``` with a python 3 environment. 

```
#enter directory with fastq input file 
r1=CV167_1.fastq
r2=CV167_2.fastq

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
