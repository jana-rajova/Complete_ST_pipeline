import pandas as pd 
import numpy as np
import os
import re
import argparse
from matplotlib import pyplot as plt
import seaborn as sns
from datetime import datetime

# The script needs three input files. The matrix, Ensembl id to symbol dataframe 
# and new positions file
log = list()
def df_to_parts(matrix):
	genes = matrix.columns.to_list()
	positions = matrix.index.to_list()
	matrix = matrix.reset_index()
	matrix.columns = [''] * len(matrix.columns)
	matrix.astype('int')
	return genes, positions, matrix

def ensembl_to_symbol(matrix,ensembl_location, log=log):
	"""
	Translates ensembl ids into a gene name to avoid species bias
	"""
	print("Translating ensembl IDs to gene names and merging duplicate columns")
	gene_df = pd.read_csv(ensembl_location,header=0, index_col=False).to_numpy()
	# print(type(gene_df))
	# print(gene_df)
	# print(gene_df.shape)
	genes_ens = matrix.iloc[0,:].values.tolist()
	counter_fin=len(genes_ens)
	counter_pos = 0
	# genes_symb = genes_ens
	#genes_ens = genes_ens[:10]
	genes_symb = list()
	for i in genes_ens:
		#print(i)
		counter_pos += 1
		try:
			gene = gene_df[np.where(gene_df == i)[0], 0].item()
			#print(gene)
			# gene = gene_df[gene_df.eq(i).any(1)]["Gene name"].item()
		except:
			gene = i
		genes_symb.append(gene)
		if counter_pos % 1000 == 0:
			print("processing entry: ", counter_pos, "/", counter_fin, sep='')
	with open(timestamp_path + sample + "_genes.txt", "w+") as output:
		print("Outputting gene lists into a file")
		for i in range(len(genes_ens)):
			output.write(str(genes_ens[i])+ " == " + str(genes_symb[i]) + '\n')
		output.close()
	duplicates = set([x for x in genes_symb if genes_symb.count(x)>1])
	matrix.columns = genes_symb
	matrix = matrix.iloc[1:,:].astype("float64")
	print(matrix.shape)
	matrix = matrix.groupby(matrix.columns, axis=1).sum()
	print(matrix.shape)
	print("ensembl length: ", len(set(genes_ens)), "; symbol length: ", len(set(genes_symb)), sep='')
	log.append("Ensembl ID to Gene Name translation:")
	log.append(duplicates)
	duplicate_pos = len(log)
	log.insert(0, duplicate_pos)
	# print("Gene names appearing more than once:")
	# print(duplicates)
	log.append(str(len(set(genes_ens))) +  " Ensembl IDs were translated into " + str(len(set(genes_symb))) + " gene names")
	matrix.to_csv("data/ST_files/temp/" + sample + "_stdata_temp.csv")
	return matrix, log

def position_correction(matrix, position_df, drop_outside=True, log=log):
	"""
	takes old coordinates from the index and transforms them into the updated ones
	"""
	print("The positions will be updated")
	counter = 0
	list_updated_pos = ["positions"]
	pos_array = None
	for pos in matrix.index[1:]:
		search = re.search("([0-9]+)x([0-9]+)", pos)
		if not pos_array is None:
			pos_array = np.append(pos_array, [[int(search.group(1)),int(search.group(2))]], axis=0)
			#print(pos_array)
		else:
			pos_array = np.array([[int(search.group(1)),int(search.group(2))]])
			#print(pos_array)
	
	for i in range(len(pos_array)):
		new_coordinates = position_df[(position_df['x']== pos_array[i,0]) & (position_df['y'] == pos_array[i,1])]
		if len(new_coordinates)>0:
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
	matrix.to_csv("data/ST_files/temp/" + sample + "_stdata_temp.csv")
	return matrix, log

def filter_by_expression(matrix, min_row_count=300, min_feat_count=2, min_features=4, log=log):
	print("Filtering features and genes by expression. Positions with less than ", min_row_count, 
		" transcripts are discarded as well as genes that do not reach at least ", min_feat_count, 
		" in at least ", min_features, " features", sep='')
	# matrix_num = matrix.iloc[:, :]

	sums = matrix.sum(axis=1)
	expr_data = sums
	expr_data.sort_values(ascending=False, inplace=True)
	# print(expr_data)
	# print(type(expr_data))
	expr_data.to_csv(timestamp_path + sample + "_transcripts_per_feature.csv")
	maxim = len(matrix.columns)
	counter = 0
	list_to_drop = list()
	for i in sums.index.tolist():
		# print(i)
		# count_arr = np.append() 
		if sums[i]<min_row_count:
			counter += 1
			list_to_drop.append(i)
	# print(list_to_drop)
	# print(type(list_to_drop[0]))
	matrix.drop(list_to_drop, axis=0, inplace=True)
	
	print(counter, " row(s) droppped as they contained less than the minimum set ", min_row_count, " transcrips", sep="")
	log.append("Filtering by expression parameters:")
	log.append("Minimum row expression count: " + str(min_row_count))
	log.append(str(counter) + " row(s) droppped as they contained less than the minimum set " + str(min_row_count) + " transcrips")
	counter = 0
	dr_counter = 0
	
	sums = matrix.mask(matrix >= min_feat_count, np.nan).isnull().sum(axis=0)
	list_to_drop = list()
	expr_data = matrix.sum(axis=0)
	expr_data.sort_values(ascending=False, inplace=True)
	expr_data.to_csv(timestamp_path + sample + "_transcripts_per_gene.csv")
	# sums.to_csv(sample + "gene_counts.csv")
	for i in sums.index.tolist():
		#print(matrix[i])
		counter += 1
		if counter%3000==0:
			print("processing column ", counter, "/", maxim, sep="")
		if not sums[i] >= min_features:
			dr_counter += 1
			list_to_drop.append(i)
	matrix.drop(list_to_drop, axis=1, inplace=True)
	
	print(dr_counter, " gene(s) dropped as their expression was not more than ", min_feat_count, " in more than ", min_features, " features", sep="")
	log.append("Minimum count at least " + str(min_feat_count) + " in at least " + str(min_features) + " features")
	log.append(str(counter) + " gene(s) dropped as their expression was not more than " + str(min_feat_count) + " in more than " + str(min_features) + " features")
	matrix.to_csv("data/ST_files/temp/" + sample + "_stdata_temp.csv")
	return matrix, log	

def quality_control_graph(step, matrix, cutoff, log):
	#fig, axs = plt.subplots(ncols=3)
	try:
		os.chdir(timestamp_path)
		tg = matrix.sum(axis=0)
		tf = matrix.sum(axis=1)
		gf = pd.Series(np.count_nonzero(matrix, axis=1))

		for scale in [True, False]:
			if scale == True:
				sns.distplot(tg.tolist(), kde=True, color='orange', label="Transcripts per gene").set_xscale('log')
			elif scale == False:
				sns.distplot(tg.tolist(), kde=True, color='orange', label="Transcripts per gene")
		#plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0)				
			plt.title(step)
			mean = tg.mean()
			median = tg.median()
			mode = tg.mode()[0].item()
			std = [mean + tg.std(), mean - tg.std()]
			plt.axvline(mean, color='red', linestyle='--')
			plt.axvline(median, color='grey', linestyle='-')
			plt.axvline(mode, color='blue', linestyle='-')
			plt.axvline(std[0], color='violet', linestyle='-')
			plt.ylabel("Frequency")
			if std[1]>0:
				plt.axvline(std[1], color='violet', linestyle='-')
			plt.legend({'Mean':mean,'Median':median, 'Mode':mode, 'Standard Deviation':std})
			if scale == True:
				plt.xlabel("Transcripts per gene - logarithmic")
				plt.savefig(sample + "_transcripts_per_gene_log" + "_" + step + ".png", bbox_inches='tight')
			elif scale == False:
				plt.xlabel("Transcripts per gene")
				plt.savefig(sample + "_transcripts_per_gene" + "_" + step + ".png", bbox_inches='tight')
		log.append("Transcripts per gene; step: " + step)
		log.append("\t Mean: " + str(mean))
		log.append("\t Median: " + str(median))
		log.append("\t Mode: " + str(mode))
		log.append("\t Standard Deviation: " + str(tg.std()))
		log.append("\t Standard Deviation bounds: " + str(std))
		# plt.show()
		plt.clf()


		sns.distplot(tf.tolist(),color='green', label="Transcripts per feature")
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

		sns.distplot(gf.tolist(), color='pink', label="Genes per feature")
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
		# plt.show()
		plt.clf()
	except:
		print(sample, " graphs not created", sep ='')
		log.append("OBS!!! GRAPHS NOT CREATED! step: " + step)	
	finally:
		os.chdir("../../../")
	return log

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description="parse out all the arguments for the ST dataframes")
	parser.add_argument("-f", "--files", type=str, help="text file with samples for processing in csv format")
	parser.add_argument("-e", "--ensembl", default="data/ensembl_files/Homo_Rat.csv", type=str, help="path to the ensembl path")
	# parser.add_argument("-s", "--spots", type=str, help="spot files")
	parser.add_argument("-s", "--selection", default=True, type=bool, help="Use all spot or only under the tissue (then False)")
	parser.add_argument("-t", "--timestamp", type=str, help="timestamp")
	args = parser.parse_args()
	
	#timestamp_path = "results_" + args.timestamp + "/"
	#os.chdir("data/ST_files")
	#os.mkdir(timestamp_path)
	timestamp_path = "data/ST_files/ST_matrix_processed/"
	#os.chdir("../../")

	samples = list()

	#files loaded:
	ensembl_location = args.ensembl
	with open ("data/ST_files/" + args.files, "r") as file:
		for line in file:
			samples.append(line.rstrip())
	print(len(samples), "sample(s) selected")
	print(samples)
	if args.selection == True:
		feat = "selection"
	else:
		feat = "all"
		
	for sample in samples:
		#loading files
		log = [sample]
		start = datetime.now()
		print(start)
		log.append("Started: " + str(start))
		try:
			matrix =  "data/ST_files/original_ST_troubleshoot/" + sample + "_stdata.csv"
			log.append("ST file: " + matrix)
			matrix = pd.read_csv(matrix, header=None, index_col=0,low_memory=False)
			print(sample, "matrix loaded successfully!")
			features = "data/ST_files/original_features/spot_data-" + feat + "-" + sample + ".csv"
			log.append("Spot_data: " + features)
			position_df = pd.read_csv(features, header=0, index_col=False)
			print(sample, feat, " features' positions loaded!")
			log.append("samples loaded successfully")
			log.append("The initial matrix dimensions are: (" + str(len(matrix.index)-1) + ", " + str(len(matrix.columns)) + ")")
			execute = True
		except:
			print(sample, "Sample loading unsuccessful!")
			execute = False
		print("execute: ", execute)
		# matrix = matrix.iloc[:500,:2000]

		"""
		Filtering selections:
		"""
		min_row_count = 300
		min_feat_count = 2
		min_features = 4

		if execute == True:
			print(matrix.shape)	
			matrix, log = position_correction(matrix, position_df, drop_outside=args.selection, log=log)
			print(matrix.iloc[:5,:3])
			print(sample, "'s feature selection and name correction done!", sep='')
			matrix, log = ensembl_to_symbol(matrix, ensembl_location, log=log)
			print(matrix.iloc[:5,:3])
			print(sample, " Ensemb IDs translated done!", sep='')
			log = quality_control_graph("Before QC", matrix=matrix, cutoff=min_row_count, log=log)
			matrix, log = filter_by_expression(matrix, min_row_count=min_row_count, min_feat_count=min_feat_count, min_features=min_features, log=log)
			print(matrix.iloc[:5,:3])
			print(sample, " gene expression filtered!", sep='')
			log = quality_control_graph("After QC", matrix=matrix, cutoff=min_row_count, log=log)
			print(matrix.shape)
			matrix = matrix.transpose()
			log.append("The final matrix dimensions are: " + str(matrix.shape))
			matrix.to_csv(timestamp_path + sample + "_stdata.csv")
			print(sample, "COMPLETED")
			end = datetime.now()
			print("Finished: ", end, sep='')
			log.append("Finished: " + str(end))
			log.append("Duration: " + str(end-start))
			print("Duration: ", str(end-start), sep='')


			#next section is for using catplot
			# gene_counts = pd.DataFrame({'trans_per_gene' : matrix.sum(axis=1).tolist()})
			# feat_counts = pd.DataFrame({'trans_per_feat' : matrix.sum(axis=0).tolist()})
			# sn.catplot(x='trans_per_gene', data=gene_counts, kind='count', color='orange', label="Transcripts per gene")
			# sn.catplot(x="trans_per_feat", data=feat_counts, kind='count', color='red', label="Transcripts per feature")
			

			log_t = log.pop(int(log[0]))
			log.pop(0)
			log.append(log_t)
			log.insert(-1, "Gene names appearing more than once:")
			with open (timestamp_path + sample + "_log.txt", "w+") as l:
				for a in log:
					l.write(str(a) + "\n")
				l.close()