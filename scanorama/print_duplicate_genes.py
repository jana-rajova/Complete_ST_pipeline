import pandas as pd

input_files = dict()
# with open("results/Union-True/DIMRED_500_previous_partial_sets/genes_included_log.txt", "r") as file1:
# 	input_files['old'] = list()
# 	for line in file1:
# 		input_files['old'].append(line.rstrip())
# with open("results/Union-True/DIMRED_500/genes_included_log.txt", "r") as file2:
# 	input_files['new'] = list()
# 	for line in file2:
# 		input_files['new'].append(line.rstrip())	
input_files['after'] = pd.read_csv("data/CN/CN56_D2_stdata.csv", header=0, index_col=0).columns.to_list()
input_files['before'] = pd.read_csv("data/CN/old/CN56_D2_expdata.csv", header=0, index_col=0).index.to_list()
common = [x for x in input_files['after'] if x in input_files['before']]
unique_old = [x for x in input_files['before'] if x not in input_files['after']]
unique_new = [x for x in input_files['after'] if x not in input_files['before']]
print('common')
print(common)
print(len(common))
print('unique old')
print(unique_old)
print(len(unique_old))
print('unique new')
print(unique_new)
print(len(unique_new))