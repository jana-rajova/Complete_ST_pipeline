import scanpy as sc
import anndata as ad
import scanorama
import os
import matplotlib.pyplot as plt
from src.classes import *
import argparse
plt.rcParams['axes.grid'] = False
import pandas as pd
# sc.logging.print_versions()
sc.set_figure_params(figsize=(8, 8))
sc.settings.verbosity = 3


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--st_h5ad_folder', type=str,
                        default='../data/st_data_pre_processed/stdata_h5ad/')
    parser.add_argument('-s', '--sample_file', type=str,
                        default='../data/ST_files/TX.txt')
    parser.add_argument('-d', '--dataset', type=str, default='TX')
    parser.add_argument('-o', '--output', type=str,
                        default='../results/Batch_corrections/scanorama/')
    parser.add_argument('-r', '--resolution', type=float, default=1.0)
    parser.add_argument('--mode', default='client')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', default=37237)
    args = parser.parse_args()

    output = args.output + args.dataset + '/'
    os.makedirs(f'{output}plt/', exist_ok=True)

    # load sample from file
    samples = pd.read_csv(args.sample_file, header=None, index_col=None)[0].to_list()
    print(f'Processing sample {" ".join(samples)}')
    sample_list = []
    for sample in samples:
        adata = sc.read_h5ad(f'{args.st_h5ad_folder}{sample}_stdata.h5ad')
        adata.raw = adata.copy()
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000, inplace=True, subset=True)
        sample_list.append(adata)

    adatas_cor = scanorama.correct_scanpy(sample_list, return_dimred=True)

    adata_spatial = ad.concat(
        adatas_cor,
        label='batch',
        uns_merge="first"
    )
    print(adata_spatial.obs)
    adata_spatial_copy = adata_spatial.copy()

    sc.pp.neighbors(adata_spatial, use_rep="X_scanorama", knn=20)
    sc.tl.umap(adata_spatial)
    sc.tl.leiden(adata_spatial, key_added="cluster", resolution=args.resolution)

    sc.pl.umap(
        adata_spatial, color=["cluster", "sample", "slide"],
        palette='Spectral', return_fig=True, wspace=0.2
    )

    plt.savefig(f'{output}plt/{args.dataset}_scanorama_UMAP.pdf', bbox_inches='tight', dpi=400)
    plt.clf()

    for sample in samples:
        print(sample)
        plot_ST(adata_spatial, sample=sample, output=f'{output}plt/', show=False)


    adata_spatial.write_h5ad(f'{output}st_adata.h5ad')

    scanorama_df = pd.DataFrame(adata_spatial.obs[['cluster', 'feature']], columns=['cluster', 'feature'])
    scanorama_df['umap1'] = adata_spatial.obsm['X_umap'][:, 0]
    scanorama_df['umap2'] = adata_spatial.obsm['X_umap'][:, 1]
    scanorama_df['sample'] = adata_spatial.obs['sample']

    scanorama_df.to_csv(f'{output}{args.dataset}_scanorama_clusters_combined.tsv', sep='\t')