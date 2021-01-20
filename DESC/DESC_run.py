import desc
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import anndata as ad
import argparse
import os
import shutil


def wells_to_combined_df(sample_file, input_folder, output):
	samples = list()
	print(os.getcwd())
	concat_well_df = pd.DataFrame()
	with open(sample_file, "r") as well_list:
		for sample in well_list:
			sample = sample.rstrip()
			samples.append(sample)
			print("Processing well: " + sample)
			well_df = pd.read_csv(input_folder + sample + "_stdata.tsv", delimiter="\t", index_col=0, header=0, dtype="object").transpose()
			#ambiguous_entry = list(filter(lambda symbol: (symbol.startswith("__ambiguous")), well_df.columns.to_list()))
			#print(ambiguous_entry)
			#well_df.drop(ambiguous_entry, inplace=True, axis=1)
			well_df = well_df.iloc[:,:].astype("float32", copy=False)
			MultiInd = [np.array([sample]*len(well_df.index)), well_df.index.to_numpy()]
			well_df.reset_index(drop=True, inplace=True)
			well_df.set_index(MultiInd, inplace=True)
			concat_well_df = pd.concat([concat_well_df, well_df])

	concat_well_df.fillna(0, inplace=True)
	print("The final Anndata has the following parameters:")
	obs = concat_well_df.index.to_frame(index=False, name=["well", "feature"])
	print(obs.iloc[:2,:2])
	var = pd.DataFrame(concat_well_df.columns, columns=["gene_symbol"], index=concat_well_df.columns)
	print(concat_well_df.iloc[:,:].shape)
	print(var.iloc[:2,:2])
	adata = ad.AnnData(X=concat_well_df.iloc[:,:].to_numpy(), obs=obs, var=var)
	sc.pp.filter_cells(adata, min_genes=1)
	sc.pp.filter_genes(adata, min_cells=1)
	print(adata)
	adata.write(output + "Concat_wells_stdata.h5ad")

	return samples, adata


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="DESC parser")
	parser.add_argument("-l", "--list_of_wells", type=str, help="list of wells to be processed")
	parser.add_argument("-i", "--input_folder", type=str, help="folder containing matrices for unique wells")
	parser.add_argument("-o", "--output", type=str, help="output folder", default="../results/")
	parser.add_argument("-t", "--timestamp", type=str, help="timestamp")
	parser.add_argument("--top_genes", type=int, help="top highly variable genes argument")
	#parser.add_argument("-r". "--louvain_res", type=str, help="resolution")
	args = parser.parse_args()
	print("start")
	output_folder = args.output + "result_" + args.timestamp + "/DESC/"

	if os.path.isdir(output_folder) == False:
		os.makedirs(output_folder)

	samples, adata = wells_to_combined_df(sample_file=args.list_of_wells, input_folder=args.input_folder, output=output_folder)

		# Look at mitochondrial genes and filter them out
	#
	mito_genes = adata.var["gene_symbol"].str.startswith('MT-')
	#adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1)/ np.sum(adata.X, axis=1)
	# add the total counts per cell as observations-annotation to adata
	#adata.obs['n_counts'] = adata.X.sum(axis=1)
	os.chdir(output_folder)

	#sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, show=False, save="DESC_mito-gene-percentage.png")	
	#No removal of reads under a treshold, can be added
	#adata = adata[adata.obs['n_genes'] < 2500, :]
	#adata = adata[adata.obs['percent_mito'] < 0.05, :]

	# Normalization
	desc.normalize_per_cell(adata, counts_per_cell_after=1e4)

	#log scaling
	desc.log1p(adata)
	adata.raw=adata
	print(adata.obs)
	print(adata.X)

	#Select highly variable genes
	sc.pp.highly_variable_genes(adata, n_top_genes=args.top_genes, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=True)
	print(adata.var)
	adata.var["gene_symbol"].to_csv("adata_var.csv")
	adata = adata[:, adata.var['highly_variable']]
	#adata.var.to_csv("adata_var.csv")
	print(os.getcwd())
	#scale
	desc.scale(adata, zero_center=True, max_value=3)
	print(adata.var)
	print(adata.shape[1])
	adata.var["gene_symbol"].to_csv("adata_var.csv")
	# Train
	if os.path.isdir("../AE"):
		shutil.rmtree("../AE")
	else:
		os.mkdir("../AE")

	adata = desc.train(adata, dims=[adata.shape[1], 32, 16], tol=0.005, n_neighbors=3, batch_size=256, louvain_resolution=[0.8], save_dir="../AE", do_tsne=True, learning_rate=300, do_umap=True, num_Cores_tsne=4, save_encoder_weights=False)
	prob_08=adata.uns["prob_matrix0.8"]
	adata.obs["max.prob0.8"]=np.max(prob_08,axis=1)
	adata.write("anndata_DESC_clustered.h5ad")
	print(adata.obs['feature'])
	# export for plotting
	cluster_dir = "DESC_cluster/"
	if not os.path.isdir(cluster_dir):
		os.mkdir(cluster_dir)
	for well in adata.obs['well'].unique():
		df_well = adata.obs[adata.obs['well']==well]
		df_well['cluster'] = df_well['desc_0.8']
		df_well[['feature', 'cluster']].to_csv(cluster_dir + well + "_DESC_clusters.tsv", sep='\t')
	
	df_all = adata.obs[['well', 'feature']]
	df_all['cluster'] = adata.obs['desc_0.8']
	df_all['UMAP_1'] = adata.obsm['X_umap0.8'][:, 0] 
	df_all['UMAP_2'] = adata.obsm['X_umap0.8'][:, 1] 
	df_all.to_csv("DESC_all_clusters.tsv", sep='\t')

	
	#tSNE plot 
	sc.pl.scatter(adata,basis="tsne0.8",color=['desc_0.8',"max.prob0.8",'TH', 'SLC6A3'], show=False, save="DESC_tSNE.png")
	#Umap plot 
	sc.pl.scatter(adata,basis="umap0.8",color=['desc_0.8',"max.prob0.8",'TH', 'SLC6A3'], show=False, save="DESC_UMAP.png")

	

