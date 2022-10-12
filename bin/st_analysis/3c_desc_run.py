import desc
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import anndata as ad
import argparse
import shutil
import os
from src.classes import *

# if __name__ == "__main__":
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--st_h5ad_folder', type=str,
					default='../data/st_data_pre_processed/stdata_h5ad/')
parser.add_argument('-s', '--sample_file', type=str,
					default='../data/ST_files/STR_normal.txt')
parser.add_argument('-d', '--dataset', type=str, default='STR')
parser.add_argument('-o', '--output', type=str,
					default='../results/Batch_corrections/desc/')
parser.add_argument('-r', '--resolution', type=float, default=0.3)
parser.add_argument('--mode', default='client')
parser.add_argument('--host', default='127.0.0.1')
parser.add_argument('--port', default=37237)
args = parser.parse_args()


output = args.output + args.dataset + '/'
os.makedirs(f'{output}/plt/', exist_ok=True)


dir = f'{output}AE/'
if os.path.exists(dir):
    shutil.rmtree(dir)
os.makedirs(dir)

samples = pd.read_csv(args.sample_file, header=None, index_col=None)[0].to_list()
print(f'Processing sample {" ".join(samples)}')
sample_list = []
for sample in samples:
	adata = sc.read_h5ad(f'{args.st_h5ad_folder}{sample}_stdata.h5ad')
	sample_list.append(adata)

adata_merge = ST_Anndata()
adata_merge.concat_anndata(sample_list, samples)
adata_merge = adata_merge.anndata

adata_merge.raw = adata_merge.copy()
sc.pp.normalize_per_cell(adata_merge, counts_per_cell_after=1e4)
sc.pp.log1p(adata_merge)
sc.pp.highly_variable_genes(adata_merge, n_top_genes=3000, subset=True)
print(adata_merge)

#scale
adata_merge = desc.scale_bygroup(adata_merge, groupby='sample')
# sc.pp.scale(adata_merge, zero_center=True, max_value=3)

adata_merge = desc.train(adata_merge,
				   dims=[adata_merge.shape[1], 32, 16],
				   tol=0.005,
				   n_neighbors=3,
				   batch_size=256,
				   louvain_resolution=[args.resolution],
				   save_dir=dir,
				   do_tsne=False,
				   learning_rate=300,
				   do_umap=True,
				   num_Cores_tsne=4,
				   save_encoder_weights=False)

for sample in samples:
	plot_ST(adata_merge, sample, show=False, output=f'{output}/plt/{sample}_desc_clusters.png', feat_max=[34, 32], color=f'desc_{args.resolution}')

print(adata_merge)
adata_merge.write_h5ad(f'{output}adata_desc.h5ad')
desc_df = adata_merge.obs[['sample', 'feature']]
desc_df['cluster'] = adata_merge.obs[f'desc_{args.resolution}']
desc_df['umap_1'] = adata_merge.obsm[f'X_umap{args.resolution}'][:, 0]
desc_df['umap_2'] = adata_merge.obsm[f'X_umap{args.resolution}'][:, 1]
desc_df.to_csv(f"{output}{args.dataset}_desc_clusters_combined.tsv", sep='\t')

sc.pl.umap(
        adata_merge, color=[f'desc_{args.resolution}', "sample", "slide"],
        palette='Spectral', return_fig=True, wspace=0.2
    )
plt.savefig(f'{output}/plt/umap_desc_clusters.png', bbox_inches='tight')


