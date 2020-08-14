import pandas as pd
import argparse
import os
import re

if __name__ == "__main__":
	# parser = argparse.ArgumentParser(description="Translate Ensembl Ids to Gene Names for the following \\	species... All these species will be merged into a chimeric file") 
	# parser.add_argument("-s", "--species", type=str, help="R=Ratus Norvegicus; H=Homo Sapiens") 
	# args = parser.parse_args()
	os.chdir("data")
	os.chdir("ensembl_files")
	#print(os.getcwd())
	species_files_list = []
	specie = ""
	with open ("ensembl_ID_to_sybol_files.txt", "r") as species:
		for line in species:
			species_files_list.append(line.rstrip())
			res = re.search("^([^_]+)_", line).group(1)
			if specie != "":
				specie += "_" + res
			else:
				specie = res
		species.close()
	print("Creating combined", specie, "file" )
	#print(species_list)

	for i in species_files_list:
		print(i)
		load_df = pd.read_csv(i, header=0, index_col=None)
		load_df.rename(columns={'Gene stable ID':'ID'}, inplace=True)
		print(load_df.columns)
		load_df["Gene name"] = load_df["Gene name"].str.upper()
		print(load_df.head(10))
		try:
			join_ensembl = pd.merge(join_ensembl, load_df, on='Gene name', how = 'outer')

			print(join_ensembl.shape)
		except:
			join_ensembl = load_df
			print(join_ensembl.shape)
	join_ensembl.set_index("Gene name", inplace=True)
	join_ensembl.to_csv(specie + ".csv")
	print(join_ensembl.head(10))


