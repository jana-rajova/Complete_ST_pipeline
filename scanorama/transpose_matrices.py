import pandas as pd
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-f', '--file', type=str, default='conf/complete_dataset.txt', help='list of matrices')
	args = parser.parse_args()
	input("have you replaced tabs by commas in the .tsv file?")
	matrix_dict = dict()
	#prefix = 'data/CN/'
	with open(args.file, 'r') as matrix_list:
		for line in matrix_list:
			line = line.strip()
			print(line)
			df_tsv = pd.read_csv(line + '.tsv', header=0, index_col=0).transpose()
			# print(df_tsv.head(10))
			print(df_tsv.head(3))
			df_tsv.to_csv(line + '.tsv')

			df_csv = pd.read_csv(line + '.tsv', header=0, index_col=0).transpose()
			# print(df_csv.head(3))
			df_csv.to_csv(line + '.csv')
	print("don't forget to replace commas in the .tsv file by tabs")