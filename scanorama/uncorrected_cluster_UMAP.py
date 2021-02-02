import numpy as np
import re
import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering as ac
import umap
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import ST_matrices_to_anndata as Utils
import scanpy as sc


# genes is again specific to origianl csv files, but maybe we can also do it for scanorama and have the dimensions there instead if I find a way to find out,which genes contributed to which columns??
# contains lines from the csv file
def initiate_anndata_file(csv_list, path):
    samples_full = list()
    sample_list = list()
    naming = ""
    with open(path + csv_list, 'r') as f:
        for line in f:
            samples_full.append(line.rstrip())
        f.close()
    #print(samples_full)
    
    for sample in samples_full:
        try:
            #print(sample)
            name = re.search("([A-Z]{2}[0-9]+_[C-E][0-9])([A-Za-z0-9-_.]*)", sample)
            #print(name.group(1), name.group(2))
            sample_list.append(name.group(1))
            naming = name.group(2)
        except OSError:
            print(sample, " not found!")
    print("adata creation parameters")
    #rint(csv_list)
    #print(path + csv_list + naming)
    adata = Utils.join_df_to_anndata(sample_list=csv_list, path=path, naming="", features=0)

    return adata

def plot_well(adata, df_dir):
    wells = set(adata.obs["well"].to_list())
    for well in wells:
        adata_well = adata[adata.obs["well"]==well, :]
        # clust_max = int(adata.obs['cluster'].max())
        #print("AND  THEEEEEN...?")
        #print(adata_well.obs['cluster'])
        well = re.search('^([A-Z0-9]+_[A-Z0-9]{2})', well).group(1)
        sc.pl.scatter(adata_well, x='X', y='Y', color='cluster', size=250, show=False, title=well + " unified, coordinates", save="_unified_coord_" + well + ".png", frameon=True)
        adata_well.obs[['feature', 'cluster']].to_csv(df_dir + well + "_scanorama_cluster.tsv", sep="\t")
        sc.pl.scatter(adata_well, x='umap1', y='umap2', color='cluster', size=100, show=False, title=well + " unified, UMAP scatter", save="_unified_UMAP_scatter_" + well + ".png", frameon=True)

def plot_separate(adata):
    wells = set(adata.obs["well"].to_list())
    if not os.path.isdir("figures/well-unique/"):
        os.makedirs("figures/well-unique/")
    if not os.path.isdir("figures/unified/"):
        os.makedirs("figures/unified/")
    
    for well in wells:
        adata_well = adata[adata.obs["well"]==well, :]
        sc.pl.scatter(adata_well, x='umap_sep_x', y='umap_sep_y', size=100, color='cluster_sep', show=False, save="_well-unique_UMAP_scatter_" + well + ".png", title=well + " well-unique UMAP scatter", frameon=True)

        sc.pl.scatter(adata_well, x='X', y='Y', color='cluster_sep', size=250, show=False, save="_well-unique_coord_" + well + ".png", title=well + " well-unique, coordinates", frameon=True)

def plot_scatter(adata):
    wells = set(adata.obs["well"].to_list())

    sc.pl.scatter(adata, x='umap1', y='umap2', color='cluster', size=30, title="UMAP with unified clusters", save="_unified_UMAP_cluster_map.png", show=False)

    sc.pl.scatter(adata, x='umap1', y='umap2', color='well', size=30, title="UMAP with coded wells", save="_unified_UMAP_well_map.png", alpha=0.7, show=False)

def umap_separate(adata, thresh_sep):
    wells = adata.obs["well"].unique().to_list()
    hc = list()
    umap_2D = np.array([])

    for well in wells:
        print(well)
        array_well = adata.X[adata.obs["well"]==well, :]
        print(array_well.shape)
        umap_clust = umap.UMAP(n_neighbors=5, min_dist=0.5, n_components=50).fit_transform(array_well)
        hc_c = ac(n_clusters=None, affinity='euclidean', linkage='ward', distance_threshold=thresh_sep)
        hc = hc + list(hc_c.fit_predict(umap_clust).astype(int, copy=False))
        # print(hc)
        # print(type(hc[0]))
        # print(len(hc))
        if not umap_2D.size == 0:
            umap_2D = np.vstack((umap_2D, umap.UMAP(n_neighbors=5, min_dist=0.5, n_components=2).fit_transform(array_well)))
        else:
            umap_2D = umap.UMAP(n_neighbors=5, min_dist=0.5, n_components=2).fit_transform(array_well)
        # print("umap_2D shape")
        # print(umap_2D.shape)
    
    adata.obs['cluster_sep'] = hc
    adata.obs['cluster_sep'] = adata.obs['cluster_sep'].astype('category')
    adata.obs['umap_sep_x'] = umap_2D[:, 0]
    adata.obs['umap_sep_y'] = umap_2D[:, 1]
    print(adata.obs.head())
    print(adata.obs.columns)

    return adata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add threshold and what are the dimensions of the scanorama dataset in quesion')
    parser.add_argument('--threshold_united', '-u', type=int, help='threshold for cluster assignment with all wells clustered together')
    #parser.add_argument('--threshold_separate', '-s', type=int, default=3, help='threshold for cluster assignment in clustering with separate wells')
    # # parser.add_argument('-d', '--dimensions', type=int, help='how many dimensions does the scanorama file have?')
    # parser.add_argument('-c', '--clustering_dimensions', type=int, help='in how many dimensions does the ahc algorithm cluster?')
    parser.add_argument('-t', '--timestamp', type=str, help='scanorama result file timestamp')
    parser.add_argument("-o", "--output", type=str, help="output folder", default="../results/")
    parser.add_argument("-f", "--stfile", type=str, default="/media/MNM-NetStorage/OrganizedSeqFiles/ST/Complete_ST_pipeline/data/ST_files/CN56.txt")
    parser.add_argument("--output_dataframe", type=str, default="CN56_uncorrected_df.tsv")
    args = parser.parse_args()
    sc.set_figure_params(figsize=(5,5))

    output_folder = args.output + "result_" +  args.timestamp + "/uncorrrected_cluster_UMAP_output"
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    os.chdir(output_folder)
    print("Current  directory: ", os.getcwd())

    adata = Utils.join_df_to_anndata(sample_list=args.stfile, 
    path="/media/MNM-NetStorage/OrganizedSeqFiles/ST/Complete_ST_pipeline/data/ST_files/ST_matrix_processed/",
    naming="_stdata.tsv", 
    features=1)

    # adata = umap_separate(adata, thresh_sep=args.threshold_separate)

    adata.obs['well_id'] = adata.obs['well'].cat.codes
    print( adata.obs['well_id'])
    adata.uns['UMAP_cluster_embedding'] = umap.UMAP(n_neighbors=5, min_dist=0.5, n_components=50).fit_transform(adata.X)
    UMAP_2D = umap.UMAP(n_neighbors=5, min_dist=0.5, n_components=2).fit_transform(adata.X)
    adata.obs['umap1'] = UMAP_2D[:, 0]
    adata.obs['umap2'] = UMAP_2D[:, 1]

    adata.obs["X"] = adata.obs["feature"].str.extract("X([0-9.]+)",expand=True).astype('float32')
    adata.obs["Y"] = adata.obs["feature"].str.extract("_([0-9.]+)",expand=True).astype('float32')
    adata.obs["Y"] = adata.obs["Y"]*(-1)

    hc = ac(n_clusters=None, affinity='euclidean', linkage='ward', distance_threshold=args.threshold_united)
    hc = hc.fit_predict(adata.uns['UMAP_cluster_embedding'])
    adata.obs['cluster'] = hc
    adata.obs['cluster'] = adata.obs['cluster'].astype('category')

    # export clusters and features
    df_dir = "uncorrected_cluster/"
    if not os.path.isdir (df_dir):
        os.mkdir(df_dir)
    print("saving final dataframe in", os.getcwd())
    adata.obs.to_csv(args.output_dataframe, sep='\t', header=True)

    print(adata)
    print(adata.obs.head())
    print(adata.obs.columns)

        
    
    # plot_separate(adata)
    plot_well(adata, df_dir=df_dir)
    plot_scatter(adata)
    adata.write("anndata_uncorrected_clustr_combined-threshold-" + str(args.threshold_united) + ".h5ad")





    

