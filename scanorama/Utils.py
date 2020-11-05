import pandas as pd
import anndata as ad
import numpy as np
import os
import sys
import scanpy as sc

def join_df_to_anndata(sample_list, path, naming, features):
	"""
	A function that creates a AnnData object from csv or tsv files
	* sample_file = wells list
	* path  = leads to the folder with files
	* name = the naming convention (e.g. _dimred_True600.csv) 
	* features = 0 if features in index, 1 if features in columns
	"""
	concat_df = pd.DataFrame()

	for sample in sample_list:			
		if not os.path.isfile(path + sample + naming):
			print(sample, " file does not exist in path ", path)
			sys.exit()
		elif naming.endswith("csv"):
			df = pd.read_csv(path + sample + naming, index_col=0, header=0)
		elif naming.endswith("tsv"):
			df = pd.read_csv(path + sample + naming, index_col=0, header=0, delimiter="\t")
	
		if features == 0:
			MultiIndex = [np.array([sample]*len(df.index)), df.index.to_list()]
		elif features == 1:
			MultiIndex = [np.array([sample]*len(df.columns)), df.columns]
			df = df.transpose
		df.reset_index(drop=True, inplace=True)
		df.set_index(MultiIndex, inplace=True)
		concat_df = pd.concat([concat_df, df])

	concat_df.fillna(0, inplace=True)
	obs = concat_df.index.to_frame(index=False, name=["well", "feature"])

	var = pd.DataFrame(concat_df.columns, columns=["gene_symbol"], index=concat_df.columns)
	adata = ad.AnnData(X=concat_df.iloc[:,:].to_numpy(), obs=obs, var=var)
	adata.obs['well'] = adata.obs['well'].astype('category')
	sc.pp.filter_cells(adata, min_genes=1)
	sc.pp.filter_genes(adata, min_cells=1)
	#print("Final AnnData object:")
	#print(adata)

	return adata


