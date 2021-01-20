#from process import load_names, merge_datasets
import os
import scanorama as sc
import pandas as pd
import re

from time import time
import numpy as np

NAMESPACE = 'panorama'
VERBOSE = 2

#this function is to extract the cell names as the cells but not the gene expression data is maintained

def index_from_tsv(data_names):
    #from config import data_names


#    print(os.getcwd())
#    for name in data_names:
#        print(name[-14:])
    #curr_path = os.getcwd()
    #os.chdir(path)
    #print(os.listdir())
    dict_pos = {}
# extract well name and use it as a dictionary key
# assign column names (the positions) to the key
    for i in data_names:
        print(i)
        name = re.search("[A-Z]{2}[0-9]+_[C-E][1-2]", i)
        print(name.group(0))
        print("the expdata file ", i, " exists: ", os.path.isfile(i + ".tsv"), "with the current path: ", os.getcwd(), "the attempted path is: ", os.getcwd(), i, ".tsv", sep = '')
        df = pd.read_csv(i + ".tsv", index_col = 0, delimiter="\t")
        print("KEY SHAPE:", df.shape)
        print(df.columns)
        print(i, i + ".tsv", name.group(0))
        dict_pos[name.group(0)] = df.columns
#    os.chdir(curr_path)
    return dict_pos
