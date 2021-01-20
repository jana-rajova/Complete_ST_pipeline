import pandas as pd
import anndata as ad
import numpy as np
import os
import sys
import scanpy as sc
import argparse

def join_df_to_anndata(sample_list, path, naming, features, export=0):
	"""
	A function that creates a AnnData object from csv or tsv files
	* sample_file = wells list
	* path  = leads to the folder with files
	* name = the naming convention (e.g. _dimred_True600.csv) 
	* features = 0 if features in index, 1 if features in columns
	"""
	concat_df = pd.DataFrame()
	# print(os.getcwd())
	# print(path + sample_list)
	with open(sample_list, "r") as sample_list:
		for sample in sample_list:	
			#print(sample_list)
			sample = sample.rstrip()
			sample = sample + naming
			# print(sample)
			print("Processing:", sample)	
			print(path + sample)
			print(os.path.isfile(path + sample))
			if not os.path.isfile(path + sample):
				print(sample, " file does not exist in path ", path)
				sys.exit()
			elif sample.endswith("csv"):
				df = pd.read_csv(path + sample, index_col=0, header=0)
				print("CSV file read")
			elif sample.endswith("tsv"):
				df = pd.read_csv(path + sample, index_col=0, header=0, delimiter="\t")
				print("TSV file read")
			print(df.shape)
			print(type(df))
			if df.size == 0:
				print(sample, "data not loaded")
				sys.exit()
		
			if features == 0:
				MultiIndex = [np.array([sample]*len(df.index)), df.index.to_list()]
			elif features == 1:
				MultiIndex = [np.array([sample]*len(df.columns)), df.columns]
				df = df.transpose()
			#print(type(df))
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
	if export == 1:
		adata.write(path + "adata_simple_export.h5ad")
	#print("Final AnnData object:")
	#print(adata)

	return adata

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Take several tsv or csv ST matrix files and outputs one AnnData instance")
	parser.add_argument("-s", type=str, help="list of wells")
	parser.add_argument("-p", type=str, help="path to the files")
	parser.add_argument("-n", type=str, help="naming of the sample")
	parser.add_argument("-f", type=int, help="0 if features in index, 1 if in columns")
	parser.add_argument("-e", type=int, default=1)
	args = parser.parse_args()

	adata = join_df_to_anndata(sample_list=args.s, path=args.p, naming=args.n, features=args.f, export=args.e)
	print(adata)
