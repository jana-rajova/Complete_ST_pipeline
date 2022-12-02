import argparse
import os
import scanpy as sc
import pandas as pd
from src.classes import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--st_h5ad_folder', type=str,
                        default='../data/st_data_pre_processed/stdata_h5ad/')
    parser.add_argument('-c', '--cluster_folder', type=str, help='file with ST samples',
                        default='../results/Batch_corrections/desc/TX/')
    # parser.add_argument('-o', '--output', type=str, help='folder to be used as an output',
    #                     default='/home/jana/Dropbox (Bjorklund Lab)/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/st_pp/')
    parser.add_argument('--mode', default='client')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', default=37237)
    args = parser.parse_args()

    print(f'Adding cluster information in folder {args.cluster_folder}')
    cluster_file = [args.cluster_folder + x for x in os.listdir(args.cluster_folder) if 'clusters_combined.tsv' in x][0]

    samples = pd.read_csv(cluster_file, sep='\t', header=0)['sample'].unique()
    samples_list = []
    for sample in samples:
        h5_file = f'{args.st_h5ad_folder}{sample}_stdata.h5ad'
        samples_list.append(sc.read_h5ad(h5_file))

    adata = ST_Anndata()
    adata.concat_anndata(samples_list, samples)
    adata.anndata.obs_names = adata.anndata.obs['sample'].astype('str') + '_' + adata.anndata.obs['feature'].astype('str')
    adata.add_cluster_info(cluster_path=cluster_file)
    dataset = args.cluster_folder.split('/')[-2]

    adata.anndata.write(f'{args.cluster_folder}{dataset}_st_adata_cluster.h5ad')

