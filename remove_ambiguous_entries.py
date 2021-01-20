import pandas as pd
import argparse

def remove_ambiguous(st_file):
	amb_entries = list(filter(lambda symbol:(symbol.startswith("__ambiguous")), st_file.index.to_list()))
	print("Ambiguous entries")
	print(amb_entries)
	print("A total of ", len(amb_entries), " was dropped!")
	st_file.drop(amb_entries, axis=0, inplace=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-l", "--list", type=str, default="../data/ST_files/ST_list_comp.txt")
	parser.add_argument("-i", "--input_folder", type=str, default="../data/ST_files/ST_matrix_processed_lower_tresh_27-10/")
	args = parser.parse_args()

	with open(args.list, "r") as sample_list:
		for line in sample_list:
			sample = line.rstrip()
			print("Processing sample ", sample)
			try:
				df = pd.read_csv(args.input_folder + sample + "_stdata.tsv", delimiter="\t", index_col=0, header=0)
				or_length = len(df.index)
				remove_ambiguous(df)
				print("Original: ", or_length, " new length: ", len(df.index), " difference: ", or_length-len(df.index))
				df.to_csv(args.input_folder + sample + "_stdata.tsv", sep="\t")
				print("Original dataframe replaced with a dataframe without ambiguous entries!")
			except FileNotFoundError:
				print("File does not exist!")
