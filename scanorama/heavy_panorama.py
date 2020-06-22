from process import load_names, merge_datasets
import os
import scanorama as sc
import pandas as pd
from get_positions import index_from_csv as ext_index
from time import time
import numpy as np
import re
import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering as ac
import umap
import mygene
import argparse
from datetime import datetime

parser =  argparse.ArgumentParser(description="Add inital parameters")
parser.add_argument("--hvg", type=int, help="how many highly variable genes do you want to include?")
parser.add_argument('-f', '--csv_file', type=str, help='csv file list')
args = parser.parse_args()

NAMESPACE = 'panorama'
VERBOSE = 2
PERPLEXITY = 5
PERPLEXITYU = 50
UNION = True
HVG = args.hvg
GENE_MERGE = False
KEEP_DIMENSIONS = True
ADDED_GENES = False
PATH = None
# DIMRED is going to change if the KEEP_DIMENSIONS variable is True
DIMRED = args.hvg
add_genes = set()

if __name__ == '__main__':
    timestamp = datetime.now().strftime("%Y%m%d_%H-%M-%S")
#   from config import data_names
#   datasets, genes_list, n_cells = load_names(data_names)
#   t0 = time()
#   in order to get an idea, what changes if the union is True of False, theUNION loop was created
#   this should create separate folders and create graphs and dataset-dimred files in the specific UNION-True or UNION-False folders
#    for UNION in (True, False):
#
# added genes allow for adding genes from another dataset and using those as a reference for teh dimensional reduction

    # troublesooting the extracting of the gene symbol from the ensembl id
    # my_gene = mygene.MyGeneInfo()
    # var = my_gene.querymany("ENSG00000260669", scopes='ensembl.gene', fields='symbol', species='human')
    # print(str(var))
    # print(type(str(var)))
    #
    # temp = re.search("symbol\W+([A-Za-z0-9.]+)", str(var))
    # print(temp.group(1))
    # print(type(temp))
    # input("happy?")

    if ADDED_GENES == True:
        try:
            print(os.getcwd())
            #add_genes = set()
            error_lines = list()
            with open("data/scanorama_conf/claudio_mg.selected_probes", 'r') as file:
                file.readline()
                for line in file:
                    try:
                        gene = re.search('^[0-9,.-]+([^,]+),[ATCG]{10}', line)
                        add_genes.add(gene.group(1).upper())
                    except:
                        error_lines.append(line)
                file.close()
            add_genes = sorted(add_genes)
            print("ADDED GENES LOADED!")
            print("added genes: ",len(add_genes))
            print("error_lines: ", len(error_lines))
        except:
            print("something went wrong")
        finally:
            print("you're in: ", os.getcwd())
        input("Proceed?")
        print(" Now iterating with UNION = ", UNION)
    else:
        print("Variable ADDED_GENES is ", ADDED_GENES)
# if there is a need to iterate through many different dimensions:
#    for DIMRED in (995, 5000):
#    print(" Now iterating with DIMRED = ", DIMRED, 'type', type(DIMRED)
#    from config import data_names
    data_names = []
    with open(args.csv_file, 'r') as f:
        for line in f:
            data_names.append(line.rstrip())
        f.close
    print('Data names loaded')
    print(data_names)
    datasets, genes_list, n_cells = load_names(data_names)
    print(data_names)
    t0 = time()
#        print('UNION = ', UNION)
# in sc.correct, the added genes argument and path might not be utilized! check that!
    print("ADDED_GENES: ", ADDED_GENES)
    datasets_dimred, datasets, genes, dimred = sc.correct(
            datasets, genes_list, added_gene_list=add_genes, ds_names=data_names, hvg = HVG,
            sigma=15, dimred=DIMRED, return_dimred=True, return_dense=True,
            union=UNION,path=PATH, added_genes=ADDED_GENES, keep_dimensions=KEEP_DIMENSIONS, gene_merge=GENE_MERGE
            )
    genes = genes.tolist()
    print(genes)
# the following lines were added for the UNION loop
    dict_pos = ext_index(data_names)
    path1 = "data/scanorama/results_" + timestamp
    path2 = "Union-" + str(UNION)
    print(path1)
    print(os.path.isdir(path1))
    if os.path.isdir(path1) is False:
        os.mkdir(path1)
    os.chdir(path1)
    if os.path.isdir(path2) is False:
        os.mkdir(path2)
    os.chdir(path2)
    pathdimred = "DIMRED_" + str(dimred)
    if os.path.isdir(pathdimred) is False:
        os.mkdir(pathdimred)
    os.chdir(pathdimred)
    i = 0
    sample_list=[]
    print("writing files: ", dimred, " into: ", pathdimred, sep="")
    for key in dict_pos.keys():
#   print(datasets_dimred[0])
        csv = 'dataset-dimred_' + str(key) + "_" + str(UNION) + str(dimred) + '.csv'
        df = pd.DataFrame(datasets_dimred[i], columns=genes)
        print('csv: ', csv)
        print("key length", len(dict_pos[key]))
        print("dimred shape:", df.shape)
        df.insert(0, "position", dict_pos[key], True)
        df.to_csv(csv, index=0, header=True)
        i += 1
        sample_list.append(csv)
    with open('heavy_panorama_log.txt', 'w') as log:
        log.write("KEEP_DIMENSIONS, " + str(KEEP_DIMENSIONS) + "\n")
        log.write("GENE_MERGE, " + str(GENE_MERGE)+ "\n")
        log.write("files," + str(len(dict_pos.keys()))+ "\n")
        log.write("genes: " + str(len(genes)) + '\n')
        log.close()
    with open('genes_included_log.txt', 'w') as glog:
        for gene in genes:
            glog.write(gene + "\n")
        glog.close()
#   print(df.columns)
    if VERBOSE:
        print('Integrated and batch corrected panoramas in {:.3f}s'
              .format(time() - t0))
    with open("csv_files.txt", 'w') as f:
        for i in sample_list:
            f.write(i + '\n')
        f.close()
    #teh following read lines were done just for verification
    with open("csv_files.txt", 'r') as f:
        print(os.getcwd())
        for line in f:
            print(line)
        f.close()
#        labels = []
#        names = []
#        curr_label = 0
#        for i, a in enumerate(datasets):
#            labels += list(np.zeros(a.shape[0]) + curr_label)
#            names.append(data_names[i])
#            curr_label += 1
#        labels = np.array(labels, dtype=int)
#        for PERPLEXITY in (10, 20):
#            print(" Now iterating DIMRED ", DIMRED, " with PERPLEXITY = ",
#                  PERPLEXITY)
#            for N_ITER in (250, 500):
#                print(" Iterating DIMRED ", DIMRED, " with ITER = ",
#                      N_ITER)
#                sc.visualize(datasets_dimred, labels,
#                             NAMESPACE + '_ds' + "_perplexity_" +
#                             str(PERPLEXITY) + "iterations" + str(N_ITER),
#                             names, multicore_tsne=False, dimred=DIMRED,
#                             perplexity=PERPLEXITY, n_iter=N_ITER)
#         Uncorrected.
#        os.chdir("../../")
#        datasets, genes_list, n_cells = load_names(data_names)
#        datasets, genes = merge_datasets(datasets, genes_list)
#        datasets_dimred = sc.dimensionality_reduce(datasets)
#        print(datasets_dimred)
#        labels = []
#        names = []
#        curr_label = 0
#        for i, a in enumerate(datasets):
#            labels += list(np.zeros(a.shape[0]) + curr_label)
#            names.append(data_names[i])
#            curr_label += 1
#        labels = np.array(labels, dtype=int)
#        os.chdir(path)
#        os.chdir(pathdimred)
#        for PERPLEXITY in (10, 20):
#            print(" Iterating DIMRED ", DIMRED, " with PERPLEXITY = ",
#                  PERPLEXITY)
#            for N_ITER in (250, 500):
#                print(" Iterating DIMRED ", DIMRED, " with N_ITER = ",
#                      N_ITER)
#                sc.visualize(datasets_dimred, labels,
#                              NAMESPACE + '_ds_uncorrected' +
#                              "_perplexity_" + str(PERPLEXITY) +
#                              "iterations" + str(N_ITER), names,
#                              perplexity=PERPLEXITY, n_iter=N_ITER, dimred=DIMRED)
#   added with the UNION loop
        os.chdir("../../../")
