import numpy as np
import pandas as pd
import anndata as ad
import argparse

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="DESC parser")
	parser.add_argument("-l", "--list_of_wells", type=str, help="list of wells to be processed")
	parser.add_argument("-i", "--input_folder", type=str, help="folder containing matrices for unique wells")
	args = parser.parse_args()
	
	concat_well_df = pd.DataFrame()

	with open(args.list_of_wells, "r") as well_list:
		for sample in well_list:
			sample = sample.rstrip()
			print("Processing well: " + sample)
			well_df = pd.read_csv(args.input_folder + sample + "_stdata.tsv", delimiter="\t", index_col=0, header=0, dtype="object").transpose()
			well_df = well_df.iloc[:,:].astype("float32", copy=False)
			MultiInd = [np.array([sample]*len(well_df.index)), well_df.index.to_numpy()]
			well_df.reset_index(drop=True, inplace=True)
			well_df.set_index(MultiInd, inplace=True)
			concat_well_df = pd.concat([concat_well_df, well_df])


	obs = concat_well_df.index.to_frame(index=False, name=["well", "feature"])
	print(obs.iloc[:2,:2])
	var = pd.DataFrame(concat_well_df.columns, columns=["gene_symbol"])
	print(concat_well_df.iloc[:,:].shape)
	print(var.iloc[:2,:2])
	adata = ad.AnnData(X=concat_well_df.iloc[:,:].to_numpy(), obs=obs, var=var)
	print(adata)