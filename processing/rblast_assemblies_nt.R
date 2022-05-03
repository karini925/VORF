#library(rBlast)
library(data.table)
options(timeout=1000000000000000000000)

#setup databse (NR/NT)
setwd("/home/keren/DATA/blast_db")

#download nt
for(i in 0:60){
	print(i)
	if(i < 10){
		i = paste(0,i, sep="")
	}
	file_d=paste("nt.", i, ".tar.gz", sep="")
	## download NR data base from NCBI
	download.file(paste("https://ftp.ncbi.nlm.nih.gov/blast/db/", file_d, sep=""), file_d, mode='wb')
}
