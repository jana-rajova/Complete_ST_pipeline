import pandas as pd
import os
import re
import argparse
from datetime import datetime
from matplotlib import pyplot as plt
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable


def load_st_files(genes, input_folder, well_file):
	wells=[]
	with open(well_file, 'r') as f:
		for line in f:
			wells.append(line.rstrip())

	st_dict = {}
	plot_dict = {}

	print("Samples to be loaded:", len(wells))
	for well in wells:
		print(well)
		well_path = input_folder + well + '_stdata.tsv'
		print(well_path)
		if os.path.isfile(well_path):
			df = pd.read_csv(well_path, header=0, index_col=0, delimiter="\t")
			print("Well", well, "with", df.shape[0], "genes and", df.shape[1], "features loaded")
			df_sum = df.iloc[0:,0:].sum(axis = 1, skipna = True)
			st_dict[well] = df_sum
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
					print(gene, "not found in well", well)
					pos_df[gene] = [0]*len(df.index)
				plot_dict[well] = pos_df
			# print(pos_df.head())
		else:
			print("File", well + "_stdata.tsv not found")

	print("Loaded", str(len(st_dict.keys())) + " wells")

	return st_dict, plot_dict


def graphs(well, res_folder, all_genes, dim1, dim2, st_plot):
	print("Processing well", well)

	fig = plt.figure(figsize=(dim1*2, dim2*3.5))

	color = 'PuRd'
	for i in range(len(all_genes)):
		ax = fig.add_subplot(dim1, dim2, i+1)
		im = ax.scatter(x=st_plot[well]['X'], y=((-1)*st_plot[well]['Y']),
							   c=st_plot[well][all_genes[i]],
							   cmap=plt.get_cmap(color),
							   s=7)
		plt.axis('off')
		ax.set_title(all_genes[i])
		ax.set(adjustable='box', aspect='equal')

		divider = make_axes_locatable(ax)
		cax = divider.append_axes('right', size='5%', pad=0.05)
		fig.colorbar(im, cax=cax, orientation='vertical')

	plt.savefig(res_folder + well + ".png")
	# plt.show()
	plt.clf()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--genes", nargs='+', type=str,
						default=['MBP', 'TH', 'SLC6A3', 'PBX1'],
						help="genes whose expression to show")
	parser.add_argument("--graph", default=True, type=bool, help="Do you want to show graphs?")
	parser.add_argument("-o", "--output", type=str, help="output folder", default="../results/gene_visualization/")
	parser.add_argument("--input_folder", type=str, default="../data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap/stdata/")
	parser.add_argument('-f', '--file', type=str, default="../data/ST_files/TX_good_quality.txt")
	parser.add_argument("--mode", default='client')
	parser.add_argument("--port", default=62543)
	args = parser.parse_args()

	plt.style.use('seaborn')
	timestamp = datetime.now().strftime("%d%m%Y_%H%M")

	st_dict, plot_dict = load_st_files(args.genes, input_folder=args.input_folder, well_file=args.file)

	"""
	Making subplots each gene in the wells
	"""
	count = 0
	res_folder = args.output + timestamp + "/"
	os.makedirs(res_folder, exist_ok=True)
	dim2 = int(math.sqrt(len(args.genes)))
	dim1 = math.ceil(len(args.genes)/dim2)

	for well in plot_dict.keys():
		graphs(well=well, res_folder=res_folder,
			   all_genes=args.genes,
			   dim1=dim1, dim2=dim2,
			   st_plot=plot_dict)
