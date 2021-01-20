import pandas as pd
import os
import numpy as np
import re
import argparse
import seaborn as sns
from matplotlib import pyplot as plt
import math



def remove_mt_rrna(well_list, input_folder, output):
	with open (well_list, "r") as file:
		for well in file:
			well = well.rstrip()
			well_df = pd.read_csv(input_folder + well + "_stdata.tsv", delimiter="\t", index_col=0, header=0)
			# print(well)
			shape_or = well_df.shape
			mito = well_df[well_df.index.str.startswith('MT-')]
			rrna = well_df[well_df.index.str.contains('RRNA')]
			count_drop = mito.to_numpy().sum() + rrna.to_numpy().sum()
			count_tot = well_df.to_numpy().sum()

			# mito_sum = mito.sum(axis=0)/well_df.sum(axis=0)
			# print(mito_sum.sort_values(ascending=False))
			well_df = well_df.drop(well_df[well_df.index.str.startswith('MT-')].index)
			#print(well_df.shape)
			well_df = well_df.drop(well_df[well_df.index.str.contains('RRNA')].index)
			# print(well_df.shape)
			log.append(well)
			log.append(str(int(shape_or[0] - well_df.shape[0])) + " mitochondrial and rRNA gene transcripts dropped in well " +  well)
			log.append("In total dropped " + str(round((count_drop/count_tot)*100, 2))  + " percent of all transcripts in the well! (" + str(int(count_drop)) + " out of " + str(int(count_tot)) +  ")")
			log.append("\n")
			well_df.to_csv(output + well + "_stdata.tsv", sep="\t")


	return log


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Parser to describe genes that have ot be in the sample for it to be let through")
	parser.add_argument("-f", "--file", type=str, help="File containing well names", default="../data/ST_files/ST_files_striatum.txt")
	parser.add_argument("-i", "--input_folder", type=str, default="../data/ST_files/ST_matrix_processed_lower_tresh_27-10/")
	parser.add_argument("-o", "--output", type=str, help="output folder", default="../data/ST_files/ST_matrix_processed_lower_tresh_MT_rRNA_removed/")
	args = parser.parse_args()

	if not os.path.isdir(args.output):
		os.makedirs(args.output)

	log = []
	log = remove_mt_rrna(well_list=args.file, input_folder=args.input_folder, output=args.output)

	with open(args.output + "dropped_genes.txt", 'w+') as l:
		for line in log:
			if line != "\n":
				l.write(line + "\n")
			else:
				l.write(line)