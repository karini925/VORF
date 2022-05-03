library(data.table)
library(rBLAST)

#setup databse (NR/NT)
setwd("/home/keren/DATA/blast_db")


untar("nr.00.tar.gz", exdir="nr")

untar("nr.01.tar.gz", exdir="nr")

## load a BLAST database (replace db with the location + name of the BLAST DB)
bl <- blast(db="/home/keren/DATA/blast_db/nr")
bl

#load in fasta files from assembly 


## query a sequence using BLAST



#save results 