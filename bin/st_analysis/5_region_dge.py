import pandas as pd
import anndata as ad
import numpy as np
import os
from src.classes import *
import argparse
import scanpy as sc

# import gseapy

sc.set_figure_params(figsize=(5,5))
change_to_original_data = False


def detect_human_regions(adata, cutoff=0.15):
    cluster_human_content = adata.obs.groupby('cluster').mean()
    cluster_human_content = cluster_human_content[cluster_human_content['G_content'] > cutoff]

    max_human = cluster_human_content[cluster_human_content['G_content'] == cluster_human_content['G_content'].max()].index
    partial = [x for x in cluster_human_content.index if x not in max_human]
    adata.obs['partial_human'] = [True if x in partial else False for x in adata.obs.cluster]
    # adata_no_partial = adata.copy()
    # #adata_no_partial = adata_no_partial[adata_no_partial.obs.partial_human == False, :]
    # print(adata_no_partial.X)
    adata.obs['TX'] = [True if x in max_human else False for x in adata.obs.cluster]
    adata.uns['regions'] = ['TX']

    for cluster in max_human:
        adata.uns['regions_dict'][str(cluster)] = 'TX'


def select_regions_clusters(adata, region_marker_dict, regions=['STR', 'CC'], cutoff=1):
    for region in regions:
        print(region)

        adata.obs[region + '_score'] = adata[:, region_marker_dict[region]].X.mean(axis=1)
        adata.obs[region + '_score'] = adata.obs[region + '_score'] - adata.obs[region + '_score'].mean()
        adata.obs[region + '_score'] = adata.obs[region + '_score'] / adata.obs[region + '_score'].std()

        score = adata.obs.groupby('cluster').mean()
        score[region] = score[region + '_score'] > cutoff
        cluster_region_map = dict(score[region])
        adata.obs[region] = adata.obs['cluster'].map(cluster_region_map)
        for cluster, value in cluster_region_map.items():
            if value == True and cluster not in adata.uns['regions_dict'].keys():
                adata.uns['regions_dict'][cluster] = region

    adata.uns['regions'] = adata.uns['regions'] + regions

    return adata


def create_pseudobulk(adata, groupby='cluster', separateby='sample'):
    replicates = adata.obs[separateby].unique()
    obs_cols = adata.obs.columns
    pseudobulk_df = pd.DataFrame()
    for replicate in replicates:
        print('creating pseudobulk for', replicate)
        replicate_adata = adata[adata.obs[separateby] == replicate, :]
        gene_df = pd.DataFrame(replicate_adata.X,
                               columns=replicate_adata.var.index.to_list(),
                               index=replicate_adata.obs[groupby])
        gene_df = gene_df.groupby(groupby).sum()
        obs_df = replicate_adata.obs.groupby(groupby).mean()
        replicate_df = pd.concat([obs_df, gene_df], axis=1)
        replicate_df[separateby] = replicate
        replicate_df = replicate_df.reset_index()
        pseudobulk_df = pseudobulk_df.append(replicate_df)
    adata_pseudobulk = ad.AnnData(X=np.array(pseudobulk_df[replicate_adata.var.index.to_list()]),
                                  obs=pseudobulk_df[[x for x in pseudobulk_df.columns if x in obs_cols]],
                                  var=replicate_adata.var)
    adata_pseudobulk.uns = adata.uns

    return adata_pseudobulk


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--st_h5ad_cluster_file', type=str, nargs='+',
                        default='../results/Batch_corrections/seurat/TX/TX_st_adata_cluster.h5ad')
    parser.add_argument('-t', '--tissue_type', type=str,
                        default='TX')
    parser.add_argument('--mode', default='client')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', default=37237)
    args = parser.parse_args()

    output_folder = f'{"/".join(args.st_h5ad_cluster_file.split("/")[:-1])}/'

    region_marker_dict= {'STR': ['PENK', 'ADORA2A', 'PPP1R1B'],
                         'SN': ['TH', 'PBX1', 'SLC6A3', 'ALDH1A1', 'DDC', 'RET'],
                         'CC': ['MBP']}
    
    adata = sc.read_h5ad(args.st_h5ad_cluster_file)
    adata.obs['cluster'] = adata.obs['cluster'].astype('str') 
    adata.uns['regions_dict'] = dict()

    samples = adata.obs['sample'].unique()
    
    detect_human_regions(adata)

    # flag penumbra of the transplant and the transplant area
    # top human transplant defined as TX, all other ones with above 10% of human content are "partial human"

    # According to Squair et al. 2021, pseudobulk has much better performance when it comes to DEG
    # for this purpose, to prevent highly expressed genes to be flagged as DEG, I combine data from one cluster in each well into a pseudobulk

    pseudobulk_adata = create_pseudobulk(adata, groupby='cluster', separateby='sample')

    # select regions, which you want to compare with (from the region_marker_dict above)
    if args.tissue_type == 'SN':
        reg = ['SN']
    elif args.tissue_type == 'TX':
        reg = ['STR']
    else:
        reg = list(region_marker_dict.keys())
               
    pseudobulk_adata = select_regions_clusters(pseudobulk_adata, region_marker_dict=region_marker_dict, regions=reg)

    pseudobulk_adata.obs['region'] = pseudobulk_adata.obs.cluster.map(pseudobulk_adata.uns['regions_dict']).astype('str')
    pseudobulk_adata.obs.region = pseudobulk_adata.obs['region'].str.replace('nan', 'other')
    pseudobulk_adata.obs.region = pseudobulk_adata.obs.region.astype('category')

    # create a dataframe thast can be used for mapping pseudobulk columns

    if args.tissue_type == 'TX':
        print(pseudobulk_adata.obs.columns)
        pseudobulk_adata.obs['pseubobuk_cluster'] = pseudobulk_adata.obs['sample'] + '_' + pseudobulk_adata.obs['cluster']
        adata.obs['pseubobuk_cluster'] = adata.obs['sample'].astype('str') + '_' + adata.obs['cluster']
        cols_to_add = ['pseubobuk_cluster', 'STR_score', 'G_content', 'RNOG_content', 'region']

        pseudo_np = pseudobulk_adata.obs[cols_to_add].to_numpy()
        # print(pseudo_np)
        for idx, col in enumerate(cols_to_add[1:]):
            pseudo_cluster_mapping = {col: attr for col, attr in pseudo_np[:, [0, idx + 1]]}
            adata.obs[f'{cols_to_add[idx + 1]}_pseudo'] = adata.obs['pseubobuk_cluster'].map(pseudo_cluster_mapping)
        adata.obs.to_csv(f'{output_folder}TX_tissue_arttributes.tsv', sep='\t')

        cols_added = ['G_content_pseudo', 'STR_score_pseudo', 'G_content']
        for col in cols_added:
            for sample in samples:
                plot_ST(adata, sample=sample, color=col, output=f'{output_folder}attribute_plots/', show=False)

    for sample in samples:
        plot_ST(adata, sample=sample, color='cluster', output=f'{output_folder}attribute_plots/', show=False)
        #export adata back
    output_file = f'{args.st_h5ad_cluster_file.split(".h5ad")[0]}_regions.h5ad'
    print(f'Output file: {output_file}')
    pseudobulk_adata.write_h5ad(output_file)

    #create heatmaps
    pseudobulk_adata.raw = pseudobulk_adata.copy()
    sc.pp.normalize_total(pseudobulk_adata)
    sc.pp.log1p(pseudobulk_adata)

    # scale
    pseudobulk_adata.X = (pseudobulk_adata.X - np.mean(pseudobulk_adata.X, axis=0)) / np.std(pseudobulk_adata.X, axis=0)

    sc.pp.highly_variable_genes(pseudobulk_adata, flavor='seurat', n_top_genes=3000, batch_key='sample', subset=True)
    top_hvg = pseudobulk_adata.var.sort_values(by='dispersions_norm', ascending=False)

    sc.tl.rank_genes_groups(pseudobulk_adata, 'region', method='t-test', key_added="t-test", use_raw=True)

    focus_region = False
    if 'SN' in pseudobulk_adata.uns['regions_dict'].values():
        focus_region = 'SN'
        contrast_region = 'other'
    elif 'TX' in pseudobulk_adata.uns['regions_dict'].values():
        focus_region = 'TX'
        contrast_region = 'STR'

    if focus_region:
        os.chdir(output_folder)
        sc.pl.rank_genes_groups(pseudobulk_adata, n_genes=25, sharey=False, key="t-test", show=False,
                                save=f'region_rank_genes_{focus_region}-{contrast_region}_DEG_sum.pdf')
        sc.pl.rank_genes_groups_tracksplot(pseudobulk_adata, groupby='region', key="t-test", groups=[focus_region],
                                           n_genes=100, save=f'tracksplot_{focus_region}-{contrast_region}_DEG_sum.pdf',
                                           show=False)
        #sc.pl.rank_genes_groups_heatmap(adata_pseudo, groupby='region', key="t-test", groups=['SN'], n_genes=100, save=False, standard_scale='var', show_gene_labels=True)
        sc.pl.heatmap(pseudobulk_adata,
                      var_names=list(pseudobulk_adata.uns['t-test']['names'][focus_region][:30]) + list(pseudobulk_adata.uns['t-test']['names'][contrast_region][:10]),
                      groupby='region', standard_scale='var', var_group_rotation=90, show=False,
                      save=f'heatmap_{focus_region}-{contrast_region}_DEG_sum.pdf', show_gene_labels=True,
                      swap_axes=True)