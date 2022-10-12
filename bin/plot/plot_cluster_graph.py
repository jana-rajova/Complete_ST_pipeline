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
    parser.add_argument('--human', type=int, default=0)
    args = parser.parse_args()

    banned = ['figures']
    folders = [args.folder] + [args.folder + x + '/' for x in os.listdir(args.folder) if os.path.isdir(args.folder + x)]
    for ban in banned:
        folders = [x for x in folders if ban not in x]

    #print(os.listdir(args.folder))
    print(folders)
    for idx, folder in enumerate(folders):
        print('Processing folder ', idx+1, '/', len(folders), '\n', folder, sep='')

        n_cluster_arr = np.array([])

        combined_well_df = [x for x in os.listdir(folder) if x.endswith('.tsv') and 'combined' in x and not 'region' in x]
        if len(combined_well_df) == 1:
            combined_well_df = combined_well_df[0]
            combined_well_df = pd.read_csv(folder + combined_well_df, header=0, index_col=None, sep='\t')
        elif len(combined_well_df)  > 1:
            print('Conflict of files, more than one present')
            break


        well_list = combined_well_df['well'].unique()
        if len(well_list) != 0:
            
            print(well_list)
            max_clusters = 0
            min_cluster = 1000
            dfs = dict([])

            for well in well_list:
                #print(well)
                well_label = well
                df = combined_well_df[combined_well_df['well'] == well]
                #print(df)
                n_clusters = max(df['cluster'])
                lowest_cluster = min(df['cluster'])
                clust = len(df['cluster'].unique())
                if n_cluster_arr.size == 0:
                    n_cluster_arr = np.array([[well_label, str(clust)]])
                    #print(n_cluster_arr)
                else:
                    n_cluster_arr = np.append(n_cluster_arr, [np.array([well_label, str(clust)])], axis=0)
                #print(n_cluster_arr)
                dfs[well_label] = (df, n_clusters)
                if n_clusters > max_clusters:
                    max_clusters = n_clusters
                if lowest_cluster < min_cluster:
                    min_cluster = lowest_cluster

            #print(n_cluster_arr)
            n_clust_df = pd.DataFrame({'well': n_cluster_arr[:, 0], 'n_clusters': n_cluster_arr[:, 1]})
            n_clust_df.to_csv(folder + "clusters_per_well.tsv", sep='\t')
            cmap = cm.get_cmap('Spectral') # Colour map (there are many others)
            cols = np.linspace(0, 1, max_clusters+1)
            handles = [plt.plot([], color=cmap(c), ls="", marker="o")[0] for c in cols]
            labels = list(range(min_cluster, max_clusters+1))
            # print(max_clusters)
            # print(labels)

            if args.human != 0:
                cluster_human_dict = combined_well_df[['human_content', 'cluster']].groupby('cluster').mean()
                cluster_human_dict = cluster_human_dict['human_content'].to_dict()

                plt.figure(figsize=(5,5*2))
                for well_label, (df, n_clusters) in dfs.items():
                    df['cluster_human_content'] = df['cluster'].map(cluster_human_dict)

                    print(well_label)
                    print(set(df["cluster"].to_list()))
                    print(n_clusters)
                    if min_cluster == 0:
                        n_clusters += 1

                    df["X"] = df['feature'].str.extract(pat = 'X([0-9.]+)_', expand=True).astype("float32")
                    df["Y"] = df['feature'].str.extract(pat = '_([0-9.]+)').astype("float32")*(-1)

                    plt.subplot(3, 1, 1)
                    ax = plt.gca()
                    ax.set_xlim([0,35])
                    ax.set_ylim([-35,0])

                    # Now here's the plot. range(len(df)) just makes the x values 1, 2, 3...
                    # df[0] is then the y values. c sets the colours (same as y values in this
                    # case). s is the marker size.
                    print(labels[:n_clusters])
                    cs = cmap(cols[df['cluster']-min_cluster])
                    sc = ax.scatter(df["X"], df["Y"], c=cs, s=30, edgecolor='None')
                    ax.legend(handles[:n_clusters], labels[:n_clusters], bbox_to_anchor=(1.05, 1), loc='upper left')

                    plt.axis('off')
                    plt.xticks([])
                    plt.yticks([])
                    plt.title(well_label + '\nSeurat Clusters')

                    plt.subplot(3, 1, 2)
                    plt.scatter(df["X"], df["Y"], c=df['cluster_human_content'], cmap='jet', s=40, edgecolor='None')
                    plt.colorbar()
                    plt.axis('off')
                    plt.xticks([])
                    plt.yticks([])
                    plt.title('Average human transcipts content per cluster')

                    plt.savefig(folder + well_label +"_human_cluster_plot.png",  bbox_inches='tight')
                    plt.clf()
            else:
                plt.figure(figsize=(5, 5))
                for well_label, (df, n_clusters) in dfs.items():
                    print(well_label)
                    print(set(df["cluster"].to_list()))
                    print(n_clusters)
                    if min_cluster == 0:
                        n_clusters += 1

                    df["X"] = df['feature'].str.extract(pat='X([0-9.]+)_', expand=True).astype("float32")
                    df["Y"] = df['feature'].str.extract(pat='_([0-9.]+)').astype("float32") * (-1)

                    plt.subplot(1, 1, 1)
                    ax = plt.gca()
                    ax.set_xlim([0, 35])
                    ax.set_ylim([-35, 0])

                    # Now here's the plot. range(len(df)) just makes the x values 1, 2, 3...
                    # df[0] is then the y values. c sets the colours (same as y values in this
                    # case). s is the marker size.
                    print(labels[:n_clusters])
                    cs = cmap(cols[df['cluster'] - min_cluster])
                    sc = ax.scatter(df["X"], df["Y"], c=cs, s=60, edgecolor='None')
                    ax.legend(handles[:n_clusters], labels[:n_clusters], bbox_to_anchor=(1.05, 1), loc='upper left')

                    plt.axis('off')
                    plt.xticks([])
                    plt.yticks([])
                    plt.title(well_label + '\nSeurat Clusters')

                    plt.savefig(folder + well_label + "_cluster_plot.png", bbox_inches='tight')
                    plt.clf()

