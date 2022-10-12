import numpy as np
import sys
import re
import umap
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.pyplot import figure


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='choose folder', default="/media/MNM-NetStorage/OrganizedSeqFiles/ST/Complete_ST_pipeline/data/Seurat_clustered_wells/TX-sel-sub/seurat_TX-sel-sub/")
    args = parser.parse_args()

    os.chdir(args.folder)
    dataset = "TX-sel-sub"
    data = [i for i in os.listdir() if i.endswith("_seurat_top_TH_subclusters_combined.tsv")]
    if len(data) == 1:
        print("Processing ", data[0])
    elif len(data) == 0:
        print("No file found!")
        sys.exit(0)
    elif len(data) >= 2:
        print("More than 1 file found")
        sys.exit(0)

    data_df = pd.read_csv(data[0], header=0, index_col=False, sep="\t")
    data_df['feature'] = data_df['feature'].str.extract(pat='(X[0-9.]+_[0-9.]+)_')
    data_df['X'] = data_df['feature'].str.extract(pat='X([0-9.]+)_').astype('float32')
    data_df['Y'] = data_df['feature'].str.extract(pat='_([0-9.]+)').astype('float32')*-1
    max_cluster = len(data_df['cluster'].unique())
    min_cluster = data_df['cluster'].unique()
    min_cluster = np.sort(min_cluster)[1]

    wells = data_df['well'].unique()
    print(wells)
    print("min:", min_cluster)
    print(data_df['cluster'].unique(), max(data_df['cluster']))

    cmap = cm.get_cmap('rainbow') # Colour map (there are many others)
    cols = np.linspace(0,1,max_cluster-1)
    handles = [plt.plot([], color=cmap(c), ls="", marker="o")[0] for c in cols]
    labels = list(range(min_cluster, max_cluster+min_cluster))
     
    for well in wells:

        print(well)
        plt.figure(figsize=(4,4))
        well_df = data_df.loc[data_df['well'] == well]
        df_TH = well_df[well_df['cluster'] != -1]
        df_TH[['feature', 'cluster']].to_csv(well + "_seurat_clusters_" + dataset + ".tsv",sep="\t")
        n_clusters = len(df_TH['cluster'].unique())
        cs = cmap(cols[df_TH['cluster']-min_cluster])
        df_surround = well_df[well_df['cluster'] == -1]
        plt.scatter(x=df_TH['X'], y=df_TH['Y'], c=cs)
        plt.scatter(x=df_surround['X'], y=df_surround['Y'], color='grey')
        plt.legend(handles[:n_clusters], labels[:n_clusters], bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title(well)
        plt.savefig(well + "_TH_top_cluster_scatter.png", bbox_inches='tight')
        plt.clf()
