
import numpy as np
import re
import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering as ac
#import hdbscan
import umap
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from UMAP_plot import create_stdata_dictionaries, join_dataframe, cust_plot
from utils import dispersion
#from utils import dispersion


def detect_diversity (clusters,expression_df=None):
    """
    the input is clusters and a dataframe of gene expression from all positions. I calculate a mean expression for each gene in the cluster - every gene has now only one row and I sort them from highest expressed to the lowest. These are stored in a new dataframe (index = cluster, header = genes as they are in the original df). In this new dataframe, I can calculate the dispersion/diversity of the gene expression and get a new list of genes, which are the most diverse througout the clusters
    """
    #expression_df = expression_df.iloc[:,:]
    expression_df.index = clusters
    cluster_dict = dict()
    for i in range(len(expression_df.index)):
        key = str(expression_df.index[i])
        if key not in cluster_dict.keys():
            cluster_dict[key] = expression_df.iloc[i,:]
        elif key in cluster_dict.keys():
            cluster_dict[key] = pd.concat([cluster_dict[key], expression_df.iloc[i,:]], axis=1)
    cols = list()
    for cluster, df in cluster_dict.items():
        cols.append(cluster)
        if isinstance(df, pd.core.frame.DataFrame):
            df = df.mean(axis = 1)
        try:
            div_df = pd.concat([div_df, df], axis=1)
        except:
            div_df = df
        # finally:
        #     input(div_df)
    div_df.columns = cols
    # cols = div_df.columns.tolist()
    # print(cols)
    cols = sorted([int(i) for i in cols])
    cols = [str(i) for i in cols]
    print(div_df)
    div_df = div_df[cols]

    return cluster_dict, div_df

def x_top_genes(div_df, x=5, plot=False):
    assert isinstance(div_df, pd.core.frame.DataFrame)
    #print(type(x))
    assert type(x)==int
    top_genes = list()
    #cols = list(div_df.columns)
    for i in div_df.columns:
        top_genes.extend(div_df.nlargest(x,i).index.tolist())
    top_genes = set(top_genes)
    x_top_df = div_df[div_df.index.isin(top_genes)]
    if plot == True:
        sns.heatmap(x_top_df, annot=False)
        plt.savefig("x_top_heatmap.png", bbox_inches='tight')
        plt.clf()
    #print(x_top_df.shape)
    return x_top_df

def x_std_gene(div_df, x=100, plot=False):
    input(div_df['TH'])
    div_std = div_df.std(axis=1).nlargest(x).index.tolist()
    input(div_df.std(axis=1))
    input(div_df.std(axis=1).nlargest(x))
    input(div_df.std(axis=1).nlargest(x).index)
    print(div_std)
    div_std_df = div_df[div_df.index.isin(div_std)]
    if plot == True:
        sns.heatmap(div_std_df, annot=False)
        plt.savefig("div_std_heatmap.png", bbox_inches='tight')
        plt.clf()
    return div_std_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Arguments for searching diversity in the data")
    parser.add_argument('-d', '--dimensions', type=int, help="in which scanorama folder should the script look for the initial dataframes?")
    parser.add_argument('-g', '--genes', type=int, default=50, help="at how many top-hit genes should it look?")
    parser.add_argument('-t', '--treshold', type=float, default=1, help="what is the treshold for hierarchical clustering")
    parser.add_argument('--timestamp', type=str)
    parser.add_argument("-o", "--output", type=str, help="output folder")
    args = parser.parse_args()

    output_folder = args.output + "/gene_distrubution_clusters_analysis" + str(args.timestamp)
    os.mkdir(output_folder)
    os.chdir(output_folder)

    # create joined dataframe and perform clustering
    data, spots, gen_dim = create_stdata_dictionaries(csv_list="csv_files.txt", path="../scanorama_output_" + args.timestamp)
    CN = join_dataframe(data).iloc[:, 1:]
    ###

    # disp = dispersion(CN)
    # highest_disp_idx = np.argsort(disp[0])[::-1]
    # print(highest_disp_idx)
    # input("cont to dispersion")
    # print(dispersion)
    # input("cont to the rest?")
    # top_genes = set(genes[highest_disp_idx[range(hvg)]])
    # for i in range(len(datasets)):
    #     gene_idx = [ idx for idx, g_i in enumerate(genes)
    #                  if g_i in top_genes ]
    #     datasets[i] = datasets[i][:, gene_idx]
    # genes = np.array(sorted(top_genes))

    ###
    # dendrogram = sch.dendrogram(sch.linkage(CN.iloc[:,:], method='ward'))
    # plt.title("dendrogram")
    # plt.show()
    clusters = ac(n_clusters=None, affinity = 'euclidean', linkage = 'ward', distance_threshold=args.treshold)
    clusters = clusters.fit_predict(CN)
    print(sorted(set(clusters)))
    cluster_tot = len(set(clusters))
    positions = dict()
    cluster_dictionary = dict()
    start = 0
    for key, value in spots.items():
        spot_list = list()
        for i in value:
            pos = re.search("X([0-9.]+)_([0-9.]+)", i)
            end = start + len(value)
            spot_list.append([float(pos.group(1)), float(pos.group(2))])
        end = start + len(value)
        positions[key] = np.array(spot_list)
        cluster_dictionary[key] = clusters[start:end]
        start = end
    CN, div_df = detect_diversity(clusters, CN)
    sns.heatmap(div_df.iloc[:60,:], annot=False)
    plt.savefig("heatmap.png", bbox_inches='tight')
    plt.clf()
    x_top_df = x_top_genes(div_df, plot=True)
    # print(x_std_gene(div_df, plot=True))
