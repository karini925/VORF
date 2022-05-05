#combine all IVA blast results 

library(data.table)
library(phyloR)
library(plyr)
library(dplyr)

setwd("/Users/kerenisaev/Google Drive/My Drive/Columbia/Courses/ComputationalGenomics/Project/Data/IVA/blastx")

blasts=list.files(pattern=".csv")

#add file id to dataset

read_file = function(f){
  ff=read.csv(f, header=F)
  sample_name=unlist(strsplit(f, "_"))[1]
  ff = as.data.table(ff %>% group_by(V1) %>% filter(row_number()==1))
  ff$sample = sample_name
  ff$specie = ""
  for(i in 1:nrow(ff)){
    print(i)
    x=ff$V2[i]
    ff$specie[i] = as.data.table(genbank2uid_tbl(x = x))$name[1]
  }
  ff=as.data.table(ff)
  return(ff)
}

all_res = llply(blasts, read_file, .progress="text")
all_res_full = as.data.table(ldply(all_res))
all_res_full$tool = "IVA"
write.table(all_res_full, file="IVA_blast_results_wspecie_names.csv", quote=F, row.names=F, sep=";")

#total length 
#number of contigs 
#score fourth from the end 

setwd("/Users/kerenisaev/Google Drive/My Drive/Columbia/Courses/ComputationalGenomics/Project/Data/VORF/blastp")
blasts=list.files(pattern=".csv")
all_res = llply(blasts, read_file, .progress="text")
all_res_full_vorf = as.data.table(ldply(all_res))
all_res_full_vorf$tool = "VORF"
write.table(all_res_full_vorf, file="VORF_blast_results_wspecie_names.csv", quote=F, row.names=F, sep=";")

#all_assembs = rbind(all_res_full, all_res_full_vorf)
