import pandas as pd
import anndata as ad
import numpy as np
import os
import re
import argparse
from sklearn.cluster import AgglomerativeClustering as ac
import scanpy as sc
from matplotlib import cm
import matplotlib.pyplot as plt
import gseapy

sc.set_figure_params(figsize=(5,5))
change_to_original_data = False

def merge_cluster_expression(cluster_paths, expression_path):
    adata = dict()
    max_cluster = 0
    for cluster_path in cluster_paths:
        print(cluster_path)
        cluster_obs = [cluster_path + x for x in os.listdir(cluster_path) if 'combined' in x and '.tsv' in x and not 'region' in x]
        signifier = cluster_path.split('Batch_corrections')[1]
        if len(cluster_obs) == 1:
            cluster_obs = cluster_obs[0]
        else:
            print('conflicting or missing files')
            sys.exit()
        cluster_obs = pd.read_csv(cluster_obs, index_col=None, header=0, sep='\t')
        cluster_obs['cluster_int'] = cluster_obs['cluster']
        cluster_obs['cluster'] += max_cluster
        # print(cluster_obs['cluster_int'])
        wells = cluster_obs['well'].unique()

        for well in wells:
            expr_temp = pd.read_csv(expression_path + well + '_stdata.tsv', index_col=0, header=0, sep='\t')
            expr_temp['gene'] = [x.split('.')[0] for x in expr_temp.index]
            expr_temp = expr_temp.groupby('gene').sum().transpose()
            expr_temp = expr_temp.reindex(cluster_obs[cluster_obs['well'] == well]['feature'])
            adata[well] = ad.AnnData(obs=cluster_obs[cluster_obs['well'] == well],
                                     var=pd.DataFrame(expr_temp.columns, columns=['gene'], index=expr_temp.columns),
                                     X=expr_temp.iloc[:, :].to_numpy())
            adata[well].obs['signifier'] = signifier
            max_cluster = max(max(adata[well].obs.cluster), max_cluster)
    adata = ad.concat(adata, label='well')
    adata.obs.cluster = adata.obs.cluster.astype('str')
    adata.obs.cluster = adata.obs.cluster.astype('category')
    return adata

def detect_human_regions(adata, cutoff=0.1):
    cluster_human_content = adata.obs.groupby('cluster').mean()
    map = dict(cluster_human_content['human_content'] > cutoff)
    cluster_human_content = cluster_human_content[cluster_human_content['human_content'] > cutoff]

    max_human = cluster_human_content[cluster_human_content['human_content'] == cluster_human_content['human_content'].max()].index
    partial = [x for x in cluster_human_content.index if x not in max_human]
    adata.obs['partial_human'] = [True if x in partial else False for x in adata.obs.cluster]
    # adata_no_partial = adata.copy()
    # #adata_no_partial = adata_no_partial[adata_no_partial.obs.partial_human == False, :]
    # print(adata_no_partial.X)
    adata.obs['TX'] = adata.obs['cluster'].map(map)
    adata.uns['regions'] = ['TX']

    for cluster in max_human:
        adata.uns['regions_dict'][cluster] = 'TX'
    print(adata.X)
    return adata

def remove_mt_genes(adata):
    adata.var['MT'] = adata.var.index.str.startswith('MT')
    adata = adata[:, adata.var.MT == False]

    return adata

def select_regions_clusters(adata, region_marker_dict, regions=('STR', 'CC'), cutoff=0.5):
    "expects "
    for region in regions:
        adata.obs[region] = adata[:, region_marker_dict[region]].X.mean(axis=1)
        adata.obs[region] = adata.obs[region + '_score']-adata.obs[region + '_score'].mean()
        adata.obs[region] = adata.obs[region + '_score']/adata.obs[region + '_score'].std()
        score = adata.obs.groupby('cluster').mean()
        score[region] = score[region] > cutoff
        map = dict(score[region])
        adata.obs[region] = adata.obs['cluster'].map(map)
        for cluster, value in map.items():
            if value == True:
                adata.uns['regions_dict'][cluster] = region

    adata.uns['regions'] = adata.uns['regions'] + regions

    return adata

def create_pseudobulk(adata, groupby='cluster', separateby='well'):
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

def export_region_information(adata, output):
    df = adata.obs
    for region in adata.uns['regions']:
        if region != 'TX':
            scores = adata.obs.groupby('cluster')[region + '_score'].mean()
            df[region + '_cluster_score'] = df['cluster'].map(dict(scores))
        else:
            scores = adata.obs.groupby('cluster')['human_content'].mean()
            df['human_content_cluster_score'] = df['cluster'].map(dict(scores))
    # df = df.groupby('cluster').mean()
    df.to_csv(output + 'region_scores_combined.tsv', sep='\t', header=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--clusters', type=str, nargs='+',
                        default=['../results/Batch_corrections/seurat/TX/1/'])
    parser.add_argument('-e', '--expression_file', type=str, default='data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap/stdata/')
    parser.add_argument('--mode', default='client')
    parser.add_argument('--port', default=56541)
    args = parser.parse_args()

    region_marker_dict= {'STR': ['PENK', 'ADORA2A', 'PPP1R1B'],
                         'SN': ['TH', 'PBX1', 'SLC6A3', 'ALDH1A1', 'DDC', 'RET'],
                         'CC': ['MBP']}

    remove_penumbra = True

    # create an anndata object from expression and cluster assignment files, remove mitochondrial genes
    adata = merge_cluster_expression(args.clusters, args.expression_file)
    adata.uns['regions_dict'] = dict()
    adata = remove_mt_genes(adata)

    # flag penumbra of the transplant and the transplant area
    # top human transplant defined as TX, all other ones with above 10% of human content are "partial human"
    adata = detect_human_regions(adata)
    if remove_penumbra:
        adata = adata[adata.obs.partial_human == False, :]
    # nonpartial_adata = adata.copy()

    # According to Squair et al. 2021, pseudobulk has much better performance when it comes to DEG
    # for this purpose, to prevent highly expressed genes to be flagged as DEG, I combine data from one cluster in each well into a pseudobulk

    pseudobulk_adata = create_pseudobulk(adata, groupby='cluster', separateby='well')

    # select regions, which you want to compare with (from the region_marker_dict above)
    pseudobulk_adata = select_regions_clusters(pseudobulk_adata, region_marker_dict=region_marker_dict, regions=['STR'])


    # create a dataframe that can be used in "stereoscope deconvolution.R" script
    adata = select_regions_clusters(adata, region_marker_dict=region_marker_dict, regions=['STR'])
    # export_region_information(adata, output=args.cluster)

    pseudobulk_adata.obs['region'] = pseudobulk_adata.obs.cluster.map(pseudobulk_adata.uns['regions_dict']).astype('str')
    pseudobulk_adata.obs.region = pseudobulk_adata.obs['region'].str.replace('nan', 'other')
    pseudobulk_adata.obs.region = pseudobulk_adata.obs.region.astype('category')


    pseudobulk_adata.raw = pseudobulk_adata.copy()
    sc.pp.normalize_total(pseudobulk_adata)
    sc.pp.log1p(pseudobulk_adata)

    # scale
    pseudobulk_adata.X = (pseudobulk_adata.X - np.mean(pseudobulk_adata.X, axis=0))/np.std(pseudobulk_adata.X, axis=0)

    sc.pp.highly_variable_genes(pseudobulk_adata, flavor='seurat', n_top_genes=3000, batch_key='well', subset=True)
    top_hvg = pseudobulk_adata.var.sort_values(by='dispersions_norm', ascending=False)

    ## uncomment the following line if you want to remove all 'other' areas
    #nonpartial_adata = nonpartial_adata[nonpartial_adata.obs[nonpartial_adata.uns['regions']].max(axis=1) == 1, :]


    sc.tl.rank_genes_groups(pseudobulk_adata, 'region', method='wilcoxon', key_added="wilcoxon", use_raw=True)
    sc.pl.rank_genes_groups(pseudobulk_adata, n_genes=100, sharey=False, key="wilcoxon")
    sc.pl.rank_genes_groups_dotplot(pseudobulk_adata, n_genes=100,
                                    key="wilcoxon",
                                    groupby='region',
                                    groups=['TX'],
                                    save=False,
                                    values_to_plot='scores'
                                    )

    sc.pl.rank_genes_groups_tracksplot(pseudobulk_adata, groupby='region', key="wilcoxon", groups=['TX'], n_genes=100, save=False)
    sc.pl.rank_genes_groups_heatmap(pseudobulk_adata, groupby='region', key="wilcoxon", groups=['TX'], n_genes=100, save=False, standard_scale='var', show_gene_labels=True)
    sc.pl.heatmap(pseudobulk_adata,
                  var_names=list(pseudobulk_adata.uns['wilcoxon']['names']['TX'][:25]) + list(pseudobulk_adata.uns['wilcoxon']['names']['STR'][:10]),
                  groupby='region', standard_scale='var', var_group_rotation=90,
                  save='heatmap_TX_STR_DEG_sum.pdf', show_gene_labels=True,
                  swap_axes=True)

    sc.pl.heatmap(pseudobulk_adata,
                  var_names=list(top_hvg[:25].index) + list(top_hvg[-25:].index),
                  groupby='cluster', standard_scale='var', var_group_rotation=90, save='heatmap_TX_HVG.pdf', swap_axes=True)

    TX_genes = pd.DataFrame(pseudobulk_adata.uns['wilcoxon']['names']['TX'], columns=['gene'])
    TX_genes['pvals_adj'] = pseudobulk_adata.uns['wilcoxon']['pvals_adj']['TX']
    TX_genes = TX_genes[TX_genes['pvals_adj'] < 0.05]['gene'].to_list()

    gene_set_names = gseapy.get_library_name(database='Human')

    enr_res = dict()
    for database in ['KEGG_2021_Human',
                     'PanglaoDB_Augmented_2021',
                     'GO_Biological_Process_2021',
                     'GO_Molecular_Function_2021',
                     'Human_Gene_Atlas']:
        print(database)
        enr_res[database] = gseapy.enrichr(gene_list=TX_genes,
                             organism='Human',
                             gene_sets=database,
                             description='pathway',
                             cutoff=0.5)
        gseapy.plot.dotplot(enr_res[database].res2d, title=database, cmap='viridis_r', top_term=15)
        plt.savefig('figures/' + database + '_dotplot.png', bbox_inches='tight')

    sc.pl.rank_genes_groups_stacked_violin(pseudobulk_adata, groupby='region', key="wilcoxon", groups=['TX'], n_genes=100, save=False)
    sc.pl.rank_genes_groups_matrixplot(pseudobulk_adata, groupby='region', key="wilcoxon", groups=['TX'], save=True, n_genes=100)

    marker_genes = {
        'astrocyte': ['AQP4'],
        'VLMC': ['COL1A1'],
        'Dopaminergic neuron': ['TH', 'ALDH1A1', 'PBX1'],
        'other': ['NORAD', 'PCSK1', 'PEG10', 'BEX1', 'GAPDH'],
        'striatum': ['PPP1R1B', 'ADORA2A', 'MBP']
    }

    sc.pl.matrixplot(pseudobulk_adata, marker_genes, groupby='region', use_raw=False, save=True)

    sc.pl.rank_genes_groups_dotplot(pseudobulk_adata, key="wilcoxon", groups="cluster", var_names=['TH', 'ALDH1A1', 'TUBB', 'COL1A1'])