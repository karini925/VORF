# %%
import numpy as np 
from sklearn.cluster import DBSCAN
import Levenshtein
import matplotlib.pyplot as plt
import os
import bisect
import pdb
import copy
import sys
from utils import *



# %%


k=8

t=read_fastq([str(sys.argv[1]),str(sys.argv[2])])
dict=create_kmer_dictionary(t,k)
dict_clean=clean_dict(dict)
proteins=assemble_proteins(copy.copy(dict_clean),k)
name=str(sys.argv[1]).split('_',1)
res=write_seq(name, proteins)
