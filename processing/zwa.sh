#!/bin/bash

#note: this code is based on the tool available here:
#https://github.com/ndovro/ZWA/blob/main/zwa.sh

#Typical run: # ./zwa.sh RAW_READS.fastq reference.fasta
source activate trinity-env 

#raw_reads=$1
#ribos=$2

#testing
raw_reads=/home/keren/DATA/RNAseq/ERR3953893.fq
ribos=/home/keren/DATA/rRNA_genbank_U13369_1.fasta #from https://www.ncbi.nlm.nih.gov/nuccore/U13369.1?report=fasta

#starting directory 
cd /home/keren/DATA/RNAseq/ZWA

#Edit PATHS so that the script can find the appropriate programs
PATH_TO_BBMAP='/home/keren/ANALYSIS/BBMap/sh/reformat.sh'
PATH_TO_fasomerecords='/home/keren/ANALYSIS/faSomeRecords/faSomeRecords.py'
PATH_TO_TRINITY='/home/keren/miniconda3/envs/trinity-env/opt/trinity-2.8.5/Trinity'

#Change number of processes and cores as needed / can be supported by your server
N=100
CPUS=8

raw_cut=`basename ${raw_reads}`
raw_ribos=`basename ${ribos}`

if [ $# -gt 0 ]; then
    echo "Your command line contains $# arguments"
else
    echo "You need to provide raw reads and a ribosomal reference"
    exit 1
fi

# create appropriate directory
dir1=${raw_cut%.fq}_$(date +%F)
mkdir ${dir1}

# copy reads aligned on ribosomal reference to the new dir
cp ${raw_reads} ${dir1}
cp ${ribos} ${dir1}

cd  ${dir1}

#CREATE BLASTN RIBOSOMAL DB
makeblastdb -dbtype 'nucl' -in  ${ribos} -parse_seqids -out ./${raw_ribos}

#ALIGN ON RIBOSOMAL
bwa index ${raw_ribos}
bwa mem -t ${CPUS} ${raw_ribos} ${raw_cut} > aligned.sam
samtools view -b -F 4 aligned.sam >ribo.bam
samtools view -b -f 4 aligned.sam>other.bam

#DETECT HYBRIDS USING SOFTCLIPPING (MAKE SURE YOUR READ NAMES DONT CONTAIN THE LETTER 'S')

samtools view ribo.bam |cut -f 1,6 |grep S | cut -f 1 > dirty.txt

samtools bam2fq ribo.bam > ribo.fastq
samtools bam2fq other.bam > other.fastq

#convert ribo fastq to fasta file
cat ribo.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > ribo.fasta

#THIS KEEPS ONLY THE READS WITH SOFTCLIPPING
${PATH_TO_fasomerecords} --fasta ribo.fasta --list dirty.txt --outfile aligned.fasta

# change fastq to fasta using bbmap's reformat
cat other.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > other.fasta

# create a clean fasta with 1 line per sequence
awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' aligned.fasta > aligned_clean.fasta

#use blastn to locate hybrid reads
bioawk=/home/keren/ANALYSIS/bioawk/bioawk
cat aligned_clean.fasta| parallel --block 50k --recstart '>' --pipe blastn -db ${raw_ribos} -task blastn -outfmt "'6 qseqid sseqid qcovs length mismatch gapopen qstart qend sstart send evalue bitscore'"  -max_target_seqs 1 -max_hsps 1 -query - | $bioawk -t '$4>19 {print $0}' > combined_results.txt

#create directory for clean reads
mkdir clean

#start the cleaning algorithm

echo "STARTING CLEANING"
$bioawk  -t '{print $1,$7,$8}' combined_results.txt | while IFS=$'\t' read dirty_seq hybrid_start hybrid_end
do
((i=i%N)); ((i++==0)) && wait
(
#the sequence we want to clean

target_seq=$($bioawk -c fastx -v r=${dirty_seq} '$name==r {print $seq }' aligned_clean.fasta)
#the target's length
length_seq=${#target_seq}


#see if the hybrid part is at the start, the end or the middle of the read
percentage10=$(perl -E "say ${length_seq}*0.1")
percentage90=$(perl -E "say ${length_seq}*0.9")
percentage5=$(perl -E "say ${length_seq}*0.05")
percentage95=$(perl -E "say ${length_seq}*0.95")

if (( $(echo "${hybrid_end} - ${hybrid_start} > ${length_seq}*0.8" |bc -l) ));
then
	echo ${dirty_seq} >> cleaninglist.txt
	continue
fi

if (( $(echo "${hybrid_start} > ${percentage5} && ${hybrid_end} < ${percentage95}" |bc -l) ));
then
echo ${dirty_seq} >> cleaninglist.txt
	continue
fi

if (( $(echo "${hybrid_end} > ${percentage90}" |bc -l) ));
then
	echo ">"${dirty_seq} >>clean_seqs.fasta
echo ${target_seq:0:${hybrid_start}} >>clean_seqs.fasta
	echo ${dirty_seq} >> cleaninglist.txt
	continue
fi

if (( $(echo "${hybrid_start} < ${percentage10}" |bc -l) ));
then
		echo ">"${dirty_seq} >>clean_seqs.fasta
echo ${target_seq:${hybrid_end}} >>clean_seqs.fasta
echo ${dirty_seq} >> cleaninglist.txt
continue
fi
) &
 
 done
wait

cp clean_seqs.fasta clean/

cd clean
cat ../other.fasta clean_seqs.fasta > for_trinity.fasta

#create theTrinity assembly on clean_hybrids AND bwa's other
#change memory and CPUs to fit your needs
#this section can be modified to use other assemblers
#Since the assembly process is stochastic multiple assemblies are created to be utilized for validation
#${PATH_TO_TRINITY} --seqType fa --single ./for_trinity.fasta --max_memory 50G --CPU ${CPUS} --full_cleanup --output trinity_combination1
#${PATH_TO_TRINITY} --seqType fa --single ./for_trinity.fasta --max_memory 50G --CPU ${CPUS} --full_cleanup --output trinity_combination2
#${PATH_TO_TRINITY} --seqType fa --single ./for_trinity.fasta --max_memory 50G --CPU ${CPUS} --full_cleanup --output trinity_combination3
#${PATH_TO_TRINITY} --seqType fa --single ./for_trinity.fasta --max_memory 50G --CPU ${CPUS} --full_cleanup --output trinity_combination4
#${PATH_TO_TRINITY} --seqType fa --single ./for_trinity.fasta --max_memory 50G --CPU ${CPUS} --full_cleanup --output trinity_combination5

