import pandas as pd
import os
import numpy as np
import re
import argparse
import seaborn as sns
from matplotlib import pyplot as plt
import math


def load_st_files(genes, st_file, timestamp, input_folder):
	# st_list = list,
	st_dict = {}
	plot_dict = {}
	#print(os.getcwd())

	with open (st_file, "r") as file:
		for line in file:
			line = line.rstrip()

			try:
				print(os.getcwd())
				print("reading ", input_folder + line + "_stdata.tsv")
				df = pd.read_csv(input_folder + line + "_stdata.tsv", header=0, index_col=0, delimiter="\t")
				#remove_ambiguous(df)
				print("shape of", line, "file:", df.shape)
				df_sum = df.iloc[0:,0:].sum(axis = 1, skipna = True)
				st_dict[line] = df_sum
				df = df.transpose()

				X = []
				Y = []
				for position in df.index.to_list():
					pos = re.search("^X([0-9.]+)_([0-9.]+)", position)
					X.append(float(pos.group(1)))
					Y.append(float(pos.group(2)))
				pos_df = pd.DataFrame({'X': X, 'Y': Y})
				for gene in genes:
					try:
						pos_df[gene] = df[gene].to_list()
					except:
						print(gene, "not found in well", line)
						pos_df[gene] = [0]*len(df.index)
					plot_dict[line] = pos_df
				# print(pos_df.head())
			except:
				print("File", line + "_stdata.tsv not found")

	print("Loaded", str(len(st_dict.keys())) + " wells")

	return st_dict, plot_dict


def remove_non_containing(st_dict, genes, strict):
	wells_passed = []
	wells_not_passed = []
	for well, data in st_dict.items():
		gene_score = []
		print("Processing well", well)
		for gene in genes:
			if not gene in data.index:
				gene_score.append(0)
			elif gene in data.index:
				gene_score.append(data.loc[[gene]].item())
				# print(data.loc[[gene]])
			# print("gene score", gene_score)
		if strict == True:
			if 0 in gene_score:
				wells_not_passed.append(well)
			else:
				wells_passed.append(well)
		else:
			if sum(gene_score) > 0:
				wells_passed.append(well)
			else:
				wells_not_passed.append(well)
	print("Wells passed:", len(wells_passed), "\nWells not passed:", len(wells_not_passed))
	return wells_passed, wells_not_passed

def graphs(must_genes, inf_genes, wells_passed, wells_not_passed, st_dict, st_plot, output, strict, timestamp):
	plt.style.use('seaborn')
	"""
	Making subplots each gene in the wells
	"""
	count = 0
	res_folder = output + "/gene_filter_results"
	all_genes = must_genes + inf_genes
	try:
		os.mkdir(res_folder)
	except:
		pass
	dim2 = int(math.sqrt(len(all_genes)))
	dim1 = int(len(all_genes)/dim2)
	if len(all_genes)%dim1 != 0:
		if dim1 > dim2:
			dim1 += 1
		else:
			dim2 += 1
	for well in st_plot.keys():
		fig, ax = plt.subplots(dim2, dim1, sharex=True, sharey=True, constrained_layout=True)
		fig.set_size_inches(dim1*2.5, dim2*2.5)
		a = 0
		b = 0
		if well in wells_passed:
			sup =  "_passed"
			fig.suptitle(well + '- passed')
			color = 'magma'
			print("Processing passed well", well)
		elif well in wells_not_passed:
			fig.suptitle(well + '- discarded')
			color = 'Blues'
			sup = "_discarded"
			print("Processing discarded well", well)
		for i in range(len(all_genes)):
			title = all_genes[i]
			if all_genes[i] in inf_genes:
				color = 'jet'
			if dim2 > 1:
				plo = ax[a, b].scatter(x=st_plot[well]['X'], y=((-1)*st_plot[well]['Y']), c=st_plot[well][all_genes[i]], cmap=plt.get_cmap(color), s=10)
				ax[a, b].set_title(title)
				if st_plot[well][all_genes[i]].sum() > 0:
					plt.colorbar(plo, ax=ax[a, b])
			else:
				plo = ax[b].scatter(x=st_plot[well]['X'], y=((-1)*st_plot[well]['Y']), c=st_plot[well][all_genes[i]], cmap=plt.get_cmap(color), s=10)
				ax[b].set_title(title)
				if st_plot[all_well][genes[i]].sum() > 0:
					plt.colorbar(plo, ax=ax[b])
			b += 1
			if b >= dim1:
				b = 0
				a += 1
		plt.savefig(res_folder + "/" + well + sup + ".png", bbox_inches='tight')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Parser to describe genes that have ot be in the sample for it to be let through")
	parser.add_argument("-g", "--must_genes", nargs='+', type=str, default='', help="genes that must be present in the dataset for it to be let through")
	parser.add_argument("-i", "--inf_genes", nargs='+', default='', type=str, help="genes that give additional info about the tissue")
	parser.add_argument("--graph", default=True, type=bool, help="Do you want to show graphs? ")
	parser.add_argument("-f", "--file", type=str, help="File with ST pointers")
	parser.add_argument("-t", "--timestamp", type=str, help="From which timestamp results do you want to take the ST files?")
	parser.add_argument("-s", "--strict", type=int, default=1, help="If True all necessary genes must be present, otherwise presence of any necessary gene will suffice")
	parser.add_argument("-o", "--output", type=str, help="output folder", default="../results/")
	parser.add_argument("--input_folder", type=str, default="../data/ST_files/ST_matrix_processed/")
	args = parser.parse_args()

	output_folder = args.output + "result_" + args.timestamp + "/"
	if not os.path.isdir(output_folder):
		os.mkdir(output_folder)

	all_genes = args.must_genes + args.inf_genes
	# print(all_genes)
	if args.strict == 0:
		strict = False
	else:
		strict = True

	st_dict, plot_dict = load_st_files(all_genes, args.file, args.timestamp, input_folder=args.input_folder)
	passed, non_passed = remove_non_containing(st_dict, args.must_genes, strict)
	print("Passed:", passed)
	print("Not Passed:", non_passed)
	graphs(must_genes=args.must_genes, inf_genes=args.inf_genes, wells_passed=passed, wells_not_passed=non_passed, timestamp=args.timestamp, st_dict=st_dict, st_plot=plot_dict, output=output_folder, strict=strict)
	print("Saving file of passed wells")
	pass_out = output_folder + "ST_files_filter_passed_" + args.timestamp + ".txt"
	well_list = output_folder + "ST_wells_filter_passed_" + args.timestamp + ".txt"

	with open (pass_out, 'w+') as output:
		for well in passed:
			output.write(args.input_folder + well + "_stdata\n")

	with open (well_list, 'w+') as output:
		for well in passed:
			output.write(well + "\n")
	print("Filtering samples based on gene presence COMPLETE!")
