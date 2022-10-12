import pandas as pd
import numpy as np
import os
import re
import argparse
from matplotlib import pyplot as plt
import seaborn as sns
from datetime import datetime
import mygene
from src.utils import *
mg = mygene.MyGeneInfo()

# The script needs three input files. The matrix, Ensembl id to symbol dataframe
# and new positions file
log = list()

def remove_ambiguous(st_file):
	amb_entries = list(filter(lambda symbol:(symbol.startswith("__")), st_file.columns))
	#print(amb_entries)
	st_file.drop(amb_entries, axis=1, inplace=True)
	print(st_file.shape)
	return amb_entries


def df_to_parts(matrix):
	genes = matrix.columns.to_list()
	positions = matrix.index.to_list()
	matrix = matrix.reset_index()
	matrix.columns = [''] * len(matrix.columns)
	matrix.astype('int')
	return genes, positions, matrix



#def ensembl_to_symbol(sample, matrix,ensembl_location, output_folder, log=log):
	# """
	# Translates ensembl ids into a gene name to avoid species bias
	# """
	# print("Translating ensembl IDs to gene names and merging duplicate columns")
	# gene_np = pd.read_csv(ensembl_location,header=0, index_col=False).to_numpy().astype('str')
	# # gene_ref = pd.read_csv(ensembl_location,header=0, index_col=0)

	# genes_ens = matrix.iloc[0,:].values.astype('str')
	# if genes_ens[0].startswith('ENS'):
	# 	print('Translating Ensembl IDs')
	# 	#print("gene_ens", genes_ens.dtype)
	# 	counter_fin=len(genes_ens)
	# 	counter_pos = 0
	# 	genes_symb = []
	#
	# 	for i in genes_ens:
	# 		#print(i)
	# 		counter_pos += 1
	# 		try:
	# 			gene = gene_np[np.where(gene_np == i)[0], 0].item()
	# 			# gene = gene_df[gene_df.eq(i).any(1)]["Gene name"].item()
	# 		except:
	# 			gene = i
	# 		genes_symb.append(gene)
	# 		if counter_pos % 1000 == 0:
	# 			print("processing entry: ", counter_pos, "/", counter_fin, sep='')
	# 	with open(output_folder + sample + "_genes.txt", "w+") as output:
	# 		print("Outputting gene lists into a file")
	# 		for i in range(len(genes_ens)):
	# 			output.write(str(genes_ens[i])+ " == " + str(genes_symb[i]) + '\n')
	# 		output.close()
	#
	# else:
	# 	genes_symb = genes_ens.tolist()
	# duplicates = set([x for x in genes_symb if genes_symb.count(x)>1])
	# matrix.columns = genes_symb
	# matrix = matrix.iloc[1:,:].astype("float64")
	# # print(matrix.shape)
	# matrix = matrix.groupby(matrix.columns, axis=1).sum()
	# print("remove ambiguous")
	# amb_entries = remove_ambiguous(matrix)
	# # print(matrix)
	# df_amb = pd.DataFrame(amb_entries)
	# # print(df_amb)
	# df_amb.to_csv(output_folder + sample + "_ambiguous_entries_removed.csv", index=False, header=False)
	#
	# print("ensembl length: ", len(set(genes_ens)), "; symbol length: ", len(set(genes_symb)), sep='')
	# log.append("Ensembl ID to Gene Name translation:")
	# df_duplicates = pd.DataFrame(duplicates)
	# df_duplicates.to_csv(output_folder + sample + "_duplicate genes.csv", index=False, header=False)
	# log.append(str(len(set(genes_ens))) +  " Ensembl IDs were translated into " + str(len(set(genes_symb))) + " gene names")
	# return matrix, log

def rename_position_index(pos_array, x_pos, y_pos):
	#print(pos_array)
	#print("processing", x_pos)
	if pos_array.size == 0:
		pos_array = np.array([[x_pos, y_pos]])
		print("initializing new posistion array")
			#print(pos_array)
	elif pos_array.size > 0:
		pos_array = np.append(pos_array, np.array([[x_pos, y_pos]]), axis=0)
	return pos_array


def position_correction(matrix, position_df, drop_outside=True, log=log):
	"""
	takes old coordinates from the index and transforms them into the updated ones
	"""
	print("The positions will be updated")
	counter = 0
	list_updated_pos = ["positions"]
	pos_array = np.array([[]], dtype=int)
	for pos in matrix.index[1:]:
		search = re.search("([0-9]+)x([0-9]+)", pos)
		x_pos = int(search.group(1))
		y_pos = int(search.group(2))
		# print("pos")
		# print(x_pos, y_pos)
		# print(pos_array)
		pos_array = rename_position_index(pos_array, x_pos, y_pos)

	for i in range(len(pos_array)):
		new_coordinates = position_df[(position_df['x'] == pos_array[i,0]) & (position_df['y'] == pos_array[i,1])]
		if len(new_coordinates) > 0:
			list_updated_pos.append("X" + str(float(new_coordinates.get(key = "new_x").values.round(2))) + "_" + str(float(new_coordinates.get(key = "new_y").values.round(2))))
		elif drop_outside == True:
			matrix.drop(str(pos_array[i,0]) + "x" + str(pos_array[i,1]), axis=0, inplace=True)
			counter += 1
		elif drop_outside == False:
			list_updated_pos.append("X" + str(pos_array[i,0]) + "_" + str(pos_array[i,1]))
	matrix.reset_index(drop=True, inplace=True)
	matrix["positions"] = list_updated_pos
	matrix = matrix.set_index("positions")
	print(counter, " positions were dropped as they lie outside the tissue", sep="")
	log.append(str(counter) + " positions were dropped as they lie outside the tissue")
	#matrix.to_csv("../data/ST_files/temp/" + sample + "_stdata_temp.csv")
	return matrix, log


def filter_by_expression(matrix, output_folder, min_row_count=300, min_feat_count=2, min_features=4, log=log):
	print("Filtering features and genes by expression. Positions with less than ", min_row_count,
		" transcripts are discarded as well as genes that do not reach at least ", min_feat_count,
		" in at least ", min_features, " features", sep='')
	# matrix_num = matrix.iloc[:, :]

	sums = matrix.sum(axis=1)
	expr_data = sums
	expr_data.sort_values(ascending=False, inplace=True)
	# print(expr_data)
	# print(type(expr_data))
	# print(os.getcwd())
	expr_data.to_csv(output_folder + sample + "_transcripts_per_feature.csv")
	maxim = len(matrix.columns)
	counter = 0
	list_to_drop_row = list()
	for i in sums.index.tolist():
		# print(i)
		# count_arr = np.append()
		if sums[i]<min_row_count:
			counter += 1
			list_to_drop_row.append(i)
	# print(list_to_drop)
	# print(type(list_to_drop[0]))
	

	print(counter, " row(s) droppped as they contained less than the minimum set ", min_row_count, " transcrips", sep="")
	log.append("Filtering by expression parameters:")
	log.append("Minimum row expression count: " + str(min_row_count))
	log.append(str(counter) + " row(s) droppped as they contained less than the minimum set " + str(min_row_count) + " transcrips")
	counter = 0
	dr_counter = 0
	# print(matrix.iloc[:5, :5])
	# print(matrix.mask(matrix >= min_feat_count, np.nan).iloc[:5, :5])
	# print(matrix.mask(matrix >= min_feat_count, np.nan).isnull().iloc[:5, :5])
	# print(matrix.mask(matrix >= min_feat_count, np.nan).isnull().sum(axis=0).iloc[:5])
	sums = matrix.mask(matrix >= min_feat_count, np.nan).isnull().sum(axis=0)
	list_to_drop_col = list()
	expr_data = matrix.sum(axis=0)
	expr_data.sort_values(ascending=False, inplace=True)
	expr_data.to_csv(output_folder + sample + "_transcripts_per_gene.csv")
	# sums.to_csv(sample + "gene_counts.csv")
	for i in sums.index.tolist():
		#print(matrix[i])
		counter += 1
		if counter%3000==0:
			print("processing column ", counter, "/", maxim, sep="")
		if not sums[i] >= min_features:
			dr_counter += 1
			list_to_drop_col.append(i)
	matrix.drop(list_to_drop_col, axis=1, inplace=True)
	matrix.drop(list_to_drop_row, axis=0, inplace=True)

	print(dr_counter, " gene(s) dropped as their expression was not more than ", min_feat_count, " in more than ", min_features, " features", sep="")
	log.append("Minimum count at least " + str(min_feat_count) + " in at least " + str(min_features) + " features")
	log.append(str(counter) + " gene(s) dropped as their expression was not more than " + str(min_feat_count) + " in more than " + str(min_features) + " features")
	#matrix.to_csv("../data/ST_files/temp/" + sample + "_stdata_temp.csv")
	return matrix, log


def quality_control_graph(step, matrix, cutoff, log, output_folder):
	#fig, axs = plt.subplots(ncols=3)
	# try:
	os.chdir(output_folder)
	tg = matrix.sum(axis=0)
	tf = matrix.sum(axis=1)
	gf = pd.Series(np.count_nonzero(matrix, axis=1))

	sns.displot(tf.tolist(),color='green', label="Transcripts per feature")
		#plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0)
	plt.title(step)
	mean = tf.mean()
	median = tf.median()
	mode = tf.mode()[0].item()
	std = [mean + tf.std(), mean - tf.std()]
	log.append("Transcripts per feature; step: " + step)
	log.append("\t Mean: " + str(mean))
	log.append("\t Median: " + str(median))
	log.append("\t Mode: " + str(mode))
	log.append("\t Standard Deviation: " + str(tf.std()))
	log.append("\t Standard Deviation bounds: " + str(std))
	plt.axvline(cutoff, color='black', linestyle=':')
	plt.axvline(mean, color='red', linestyle='--')
	plt.axvline(median, color='grey', linestyle='-')
	plt.axvline(mode, color='blue', linestyle='-')
	plt.axvline(std[0], color='violet', linestyle='-')
	plt.ylabel("Frequency")
	if std[1]>0:
		plt.axvline(std[1], color='violet', linestyle='-')
	plt.legend({'Cutoff':cutoff, 'Mean':mean,'Median':median, 'Mode':mode, 'Standard Deviation':std})
	plt.xlabel("Transcripts per feature")
	plt.savefig(sample +  "_transcripts_per_feature" + "_" + step + ".png", bbox_inches='tight')
	# plt.show()
	plt.clf()

	sns.displot(gf.tolist(), color='pink', label="Genes per feature")
	plt.title(step)
	mean = gf.mean()
	median = gf.median()
	mode = gf.mode()[0].item()
	std = [mean + gf.std(), mean - gf.std()]
	log.append("Genes per feature; step: " + step)
	log.append("\t Mean: " + str(mean))
	log.append("\t Median: " + str(median))
	log.append("\t Mode: " + str(mode))
	log.append("\t Standard Deviation: " + str(gf.std()))
	log.append("\t Standard Deviation bounds: " + str(std))
	plt.axvline(mean, color='red', linestyle='--')
	plt.axvline(median, color='grey', linestyle='-')
	plt.axvline(mode, color='blue', linestyle='-')
	plt.axvline(std[0], color='violet', linestyle='-')
	plt.xlabel("Genes per feature")
	plt.ylabel("Frequency")
	if std[1]>0:
		plt.axvline(std[1], color='violet', linestyle='-')
	plt.legend({'Mean':mean,'Median':median, 'Mode':mode, 'Standard Deviation':std})
	plt.savefig(sample + "_genes_per_feature" + "_" + step + ".png", bbox_inches='tight')
	plt.clf()

	
	os.chdir(home_folder)
	#print("After plotting the results, we're in: ", os.getcwd())
	return log


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description="parse out all the arguments for the ST dataframes")
	parser.add_argument("-s", "--selection", default=True, type=bool, help="Use all spot or only under the tissue (then False)")
	parser.add_argument("-o", "--output", type=str, help="path to the output matrix files folder", default="./data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap_t/")
	parser.add_argument("--original_matrix_folder", type=str, default="./data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap/test/")
	parser.add_argument("--feature_folder", type=str, default="../Images_rev1/")
	parser.add_argument("--mode", default='client')
	parser.add_argument("--port", default=62543)
	args = parser.parse_args()

	#samples = list()
	home_folder = os.getcwd()
	print(home_folder)
	
	#ext_list = ['_noHTN5MAFXX_geneMatrix.tsv']

	samples = list(set(['_'.join(x.split('_', 2)[:2]) for x in os.listdir(args.original_matrix_folder) if x.endswith('.tsv')]))

	print(len(samples), "sample(s) selected")
	print(samples)
	if args.selection == True:
		feat = "selection"
	else:
		feat = "all"
	sample_list = []

	for sample in samples:
		#loading files
		log = [sample]
		start = datetime.now()
		print(start)
		log.append("Started: " + str(start))
		matrix =  [x for x in os.listdir(args.original_matrix_folder) if sample in x]
		print("Files for well found:", matrix)
		if len(matrix) == 1:
			adata = STLoader(args.original_matrix_folder, sample)
			adata.correct_feature_position(args.feature_folder)
			adata = adata.anndata
			sample_list.append(adata)
		elif len(matrix) > 1:
			print("MORE THAN 1 FILE/SAMPLE - dropping ", sample)
		elif len(matrix) < 1:
			print(sample, "Sample loading unsuccessful!")
		# print("execute: ", execute)
		# matrix = matrix.iloc[:200,:5000]

		"""
		Filtering selections:
		"""
		min_row_count = 300
		min_feat_count = 2
		min_features = 2

	adata_all = st_anndata(sample_list, samples)
	adata_all.transalate_ensembl()
	adata_all.merge_gene_symbol_duplicates()
	adata_all = adata_all.anndata
	if not os.path.isdir(args.output):
			print("Creating output folder anew!")
			os.makedirs(args.output)
	for sample in adata_all.obs['sample'].unique():
		adata_s = adata_all[adata_all.obs['sample']==sample, :]
		matrix = pd.DataFrame(adata_s.X, index=adata_s.obs['feature'], columns=[adata_s.var.index.to_list()]).transpose()
		log = quality_control_graph("Before QC", matrix=matrix, cutoff=min_row_count, log=log, output_folder=args.output)
		matrix, log = filter_by_expression(matrix, min_row_count=min_row_count, min_feat_count=min_feat_count, min_features=min_features, log=log, output_folder=args.output)
		# print(matrix.iloc[:5,:3])
		print(sample, " gene expression filtered!", sep='')
		log = quality_control_graph("After QC", matrix=matrix, cutoff=min_row_count, log=log, output_folder=args.output)
		# print(matrix.shape)
		log.append("The final matrix dimensions are: " + str(matrix.shape))
		#print(os.getcwd())
		matrix = matrix.transpose()
		matrix.to_csv(args.output + sample + '_stdata.tsv', sep="\t")
		print(sample, "COMPLETED")
		end = datetime.now()
		print("Finished: ", end, sep='')
		log.append("Finished: " + str(end))
		log.append("Duration: " + str(end-start))
		print("Duration: ", str(end-start), sep='')

		# log_t = log.pop(int(log[0]))
		# log.pop(0)
		# log.append(log_t)
		# log.insert(-1, "Gene names appearing more than once:")
		with open (args.output + sample + "_log.txt", "w+") as l:
			for a in log:
				l.write(str(a) + "\n")
			l.close()
