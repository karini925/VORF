#library(rBlast)
library(data.table)

#setup databse (NR/NT)
setwd("/home/keren/DATA/blast_db")

#download the db
for(i in 0:57){
	print(i)
	if(i < 10){
		i = paste(i,i, sep="")
	}
	file_d=paste("nr.", i, ".tar.gz", sep="")
	## download NR data base from NCBI
	download.file(paste("https://ftp.ncbi.nlm.nih.gov/blast/db/", file_d, sep=""), file_d, mode='wb')
}

## load a BLAST database (replace db with the location + name of the BLAST DB)
#bl <- blast(db="./16S_rRNA_DB/nr")
#bl

#load in fasta files from assembly 


## query a sequence using BLAST



#save results 