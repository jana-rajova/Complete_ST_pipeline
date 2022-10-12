import scanpy as sc
import numpy as np
import re
import umap
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from matplotlib import cm
from scipy.interpolate import griddata 
# import ST_matrices_to_anndata as Utils
import scanpy as sc
import matplotlib.colors as colors
from matplotlib.pyplot import figure
import copy

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='choose folder with cluster files')
    args = parser.parse_args()


    os.chdir(args.folder)
    file_list = [i for i in os.listdir() if i.endswith("tsv")]
    print(file_list)
    max_clusters = 0

    for combined_df in file_list:
        df = pd.read_csv(combined_df, index_col=False, header=0, sep="\t")
        if max_clusters < len(df['cluster'].unique()):
            max_clusters = len(df['cluster'].unique())
    print(max_clusters, " clusters found among datasets")



    for combined_df in file_list:
        #print(well)
        dataset = re.search("^([A-Za-z0-9-]+)_", combined_df).group(1)
        method = re.search("^[A-Za-z0-9-]+_([A-Za-z]+)_", combined_df).group(1)
        print(method)
        df = pd.read_csv(combined_df, index_col=False, header=0, sep="\t")
        df['well'] = df['well'].str.extract(pat = '([A-Z0-9]+_[C-E][1-2])').astype("str")
        n_clusters = len(df['cluster'].unique())
        if method == 'DESC':
            df.columns = ['', 'well', 'feature', 'cluster', 'umap1', 'umap2']

        cmap = cm.get_cmap('Spectral') # Colour map (there are many others)
        cols = np.linspace(0, 1, max_clusters)
        cols_loc = np.linspace(0, n_clusters/max_clusters, n_clusters)
        handles = [plt.plot([], color=cmap(c), ls="", marker="o")[0] for c in cols]
        labels = list(range(1, max_clusters+1))

        ax = plt.gca()
        cs = cmap(cols_loc[df['cluster']-1])
        print(df['cluster'])
        print(cs)
        sc = ax.scatter(df['umap1'], df['umap2'], c=cs, s=20, edgecolor='None')
        ax.legend(handles[:n_clusters], labels[:n_clusters], bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig(dataset + "_" + method +"_by_cluster_plot.png",  bbox_inches='tight')
        plt.clf()
    

        print(method)
        well_dict = dict()
        i = 0
        for well in df['well'].unique():
            well_dict[well] = i
            i += 1
        df['well_id'] = [well_dict[well] for well in df['well']]
        n_wells = len(df['well'].unique())
        
        cmap = cm.get_cmap('jet') # Colour map (there are many others)
        cols_w = np.linspace(0, 1, n_wells)
        cs_w = cmap(cols_w[df['well_id']])
        plt.scatter(df['umap1'], df['umap2'], c=cs_w, s=20, edgecolor='None')
        
        handles = [plt.plot([], color=cmap(c), ls="", marker="o")[0] for c in cols_w]
        labels = df['well'].unique()
        plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig(dataset + "_" + method +"_by_well_plot.png",  bbox_inches='tight')
        plt.clf()
