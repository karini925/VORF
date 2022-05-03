# %%
import numpy as np 
from sklearn.cluster import DBSCAN
import Levenshtein
import matplotlib.pyplot as plt
import os
import bisect
import pdb
import copy
from utils import *



# %%

datafiles={'COVID':['../Data/HIV/CV167_1.fastq','../Data/HIV/CV167_2.fastq'], 'HIV':['../Data/HIV/CV167_1.fastq','../Data/HIV/CV167_2.fastq'], 'Lassa':['../Data/Lassa/K7_1.fastq', '../Data/Lassa/K7_1.fastq']}
k=8

sequences={}
for name,file_list in datafiles.items():

    t=read_fastq(file_list)
    dict=create_kmer_dictionary(t,k)
    dict_clean=clean_dict(dict)
    proteins=assemble_proteins(copy.copy(dict_clean),k)
    res=write_seq(name, proteins)

    sequences[name]=res



