import os
import pandas as pd
import seaborn as sms
from bin.region_dge import *
import numpy as np
import anndata as ad
import scanorama
from matplotlib.pyplot import figure

plt.rcParams["figure.figsize"] = (5,5)
plt.rcParams["figure.dpi"] = 80

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--clusters', type=str, nargs='+',
                        default=['../results/Batch_corrections/seurat/TX/1/', '../results/Batch_corrections/seurat/str-intact/1/'])
    parser.add_argument('-e', '--expression_file', type=str,
                        default='data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap/stdata/')
    parser.add_argument('--mode', default='client')
    parser.add_argument('--port', default=56541)
    args = parser.parse_args()

    region_marker_dict = {'STR': ['PENK', 'ADORA2A', 'PPP1R1B'],
                          'SN': ['TH', 'PBX1', 'SLC6A3', 'ALDH1A1', 'DDC', 'RET'],
                          'CC': ['MBP', 'TF']}

    """ 
    GET READY
    
    Read in all the dataframes, merge them and put the alg/tiss/res as a 'batch key'
    Give each cluster a serial nuber so that we have unique dictionary (i.e. 19 - seurat/str-intact/1__1)
    create pseudobulks of clusters
    create a label for the cluster and which tissue it comes for (so that it is at least a bit informative)
        generate scores for tissues of interest based on marker expression
        Decide how to select where to put th tissue - theoretically should not be that hard, but you never know with the noise
        merge batch key (str-intact/TX) and region (i.e.STR from str-intact vs STR from TX)
    create spearman correlation between tissues and clusters
    """
    remove_penumbra = False

    # create an anndata object from expression and cluster assignment files, remove mitochondrial genes
    adata = merge_cluster_expression(args.clusters, args.expression_file)
    adata.uns['regions_dict'] = dict()
    adata = remove_mt_genes(adata)

    adata_list = [adata[adata.obs.well == well, : ] for well in adata.obs.well.unique()]

    for a in adata_list:
        sc.pp.normalize_total(a, inplace=True)
        sc.pp.log1p(a)
        a.var['ribo'] = a.var.index.str.startswith(("RPS","RPL"))
        a.var['mito'] = a.var.index.str.startswith("MT-")
        sc.pp.highly_variable_genes(a, flavor="seurat", n_top_genes=3000, inplace=True)

    # is it better to do in on HVG?
    adatas_cor = scanorama.correct_scanpy(list(adata_list), return_dimred=True, dimred=50)

    adata_spatial = ad.concat(
        adatas_cor,
        label='batch',
        uns_merge="unique"
    )
    adata_spatial_copy = adata_spatial.copy()

    sc.pp.neighbors(adata_spatial, use_rep="X_scanorama", knn=20)
    sc.tl.umap(adata_spatial)
    sc.tl.leiden(adata_spatial, key_added="cluster", resolution=1)

    sc.pl.umap(
        adata_spatial, color=["cluster", "well", "signifier", 'PENK', 'PDYN', 'CCK', 'MBP'],
        palette='Spectral', return_fig=True, wspace=0.2
    )


    # plt.savefig('scanorama_UMAP.png', bbox_inches='tight', dpi=400)
    plt.show()
    # flag penumbra of the transplant and the transplant area
    # top human transplant defined as TX, all other ones with above 10% of human content are "partial human"
    adata_spatial = detect_human_regions(adata_spatial)


    # According to Squair et al. 2021, pseudobulk has much better performance when it comes to DEG
    # for this purpose, to prevent highly expressed genes to be flagged as DEG, I combine data from one cluster in each well into a pseudobulk

    pseudobulk_adata = create_pseudobulk(adata_spatial, groupby='cluster', separateby='well')

    # select regions, which you want to compare with (from the region_marker_dict above)
    pseudobulk_adata = select_regions_clusters(pseudobulk_adata, region_marker_dict=region_marker_dict, regions=['STR'])

    # create a dataframe that can be used in "stereoscope deconvolution.R" script
    adata = select_regions_clusters(adata, region_marker_dict=region_marker_dict, regions=['STR'])
    adata.obs['cluster_well'] = adata.obs.cluster.astype('str') + '_' + adata.obs.well.astype('str')
    adata.obs['cluster_well'] = adata.obs['cluster_well'].astype('category')
    sc.tl.dendrogram(adata, groupby='cluster_well', cor_method='spearman', use_rep='X_scanorama')
    sc.pl.dendrogram(adata, groupby='cluster_well')
    df = pd.DataFrame(np.transpose(adata.obs), columns=adata.obs.cluster_well)
    df_corr = df.corr(method='spearman')
    sns.clustermap(df_corr)
    plt.show()

    # export_region_information(adata, output=args.cluster)

    pseudobulk_adata.obs['region'] = pseudobulk_adata.obs.cluster.map(pseudobulk_adata.uns['regions_dict']).astype(
        'str')
    pseudobulk_adata.obs.region = pseudobulk_adata.obs['region'].str.replace('nan', 'other')
    pseudobulk_adata.obs.region = pseudobulk_adata.obs.region.astype('category')

    pseudobulk_adata.raw = pseudobulk_adata.copy()
    sc.pp.normalize_total(pseudobulk_adata)
    sc.pp.log1p(pseudobulk_adata)

    # scale
    pseudobulk_adata.X = (pseudobulk_adata.X - np.mean(pseudobulk_adata.X, axis=0)) / np.std(pseudobulk_adata.X, axis=0)

    sc.pp.highly_variable_genes(pseudobulk_adata, flavor='seurat', n_top_genes=3000, batch_key='well', subset=True)
    top_hvg = pseudobulk_adata.var.sort_values(by='dispersions_norm', ascending=False)

    ## uncomment the following line if you want to remove all 'other' areas
    # nonpartial_adata = nonpartial_adata[nonpartial_adata.obs[nonpartial_adata.uns['regions']].max(axis=1) == 1, :]

    sc.tl.rank_genes_groups(pseudobulk_adata, 'cluster', method='wilcoxon', key_added="wilcoxon", use_raw=True)
    sc.pl.rank_genes_groups(pseudobulk_adata, n_genes=30, sharey=False, key="wilcoxon")
    sc.pl.rank_genes_groups_dotplot(pseudobulk_adata, n_genes=100,
                                    key="wilcoxon",
                                    groupby='cluster',
                                    # groups=['STR'],
                                    save=False,
                                    values_to_plot='scores'
                                    )

    sc.pl.rank_genes_groups_tracksplot(pseudobulk_adata, groupby='region', key="wilcoxon", groups=['TX'], n_genes=100,
                                       save=False)
    sc.pl.rank_genes_groups_heatmap(pseudobulk_adata, groupby='region', key="wilcoxon", groups=['TX'], n_genes=100,
                                    save=False, standard_scale='var', show_gene_labels=True)