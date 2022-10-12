import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
from textwrap import wrap


def remove_ambiguous(df):
	ambg = [col for col in df.columns.to_list() if col.startswith('__ambig')]
	df = df.drop(ambg, axis=1)
	return df

def read_positional_file(pos_file):
	pos_df = pd.read_csv(pos_file, header=0, index_col=None, sep='\t', dtype='str')
	pos_df['xy'] = pos_df['x'] + 'x' + pos_df['y']
	pos_df['acc_xy'] = pos_df['new_x'] + 'x' + pos_df['new_y']
	pos_df = pos_df[['xy', 'acc_xy']]
	return pos_df


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-w', '--well_file', type=str, help='choose list with wells',
						default="../data/ST_files/SN-TX_all.txt")
	parser.add_argument('-f', '--folder', type=str, help='folder with well matrices',
						default="../../STARsolo_reanalysis/ST_matrices_corrected_PetterREFs/")
	parser.add_argument('--feature_folder', type=str,
						default='../../Images_rev1/')
	parser.add_argument("--mode", default='client')
	parser.add_argument("--port", default=62543)
	args = parser.parse_args()

	wells = list()
	with open(args.well_file, 'r') as f:
		for line in f:
			wells.append(line.rstrip())
	print(wells)
	for well in wells:
		well_file = [x for x in os.listdir(args.folder) if x.startswith(well) and x.endswith("stdata.tsv")]
		if len(well_file) != 1:
			print("CONFLICTING FILES FOR WELL", well, "FOUND! Well will not be processed")
			print(well_file)
			execute = 0
			for i in well_file:
				if 'noHTN5MAFXX' in i:
					well_path = args.folder + well_file[0]
					pos_file = args.feature_folder + well + '-spot_data-selection.tsv'
					execute = 1
					print('noHTN5MAFXX file selected')
		elif len(well_file) == 0:
			print(well, "not found")
			execute = 0
		else:	
			well_path = args.folder + well_file[0]
			pos_file = args.feature_folder + well + '-spot_data-selection.tsv'
			execute = 1


#load well dataframe

		if execute == 1:
			print("Reading in well", well)
			well_df = pd.read_csv(well_path, sep="\t", header=0, index_col=None)
			columns = well_df.columns.to_list()
			columns[0] = 'xy'
			well_df.columns = columns
			well_df = remove_ambiguous(well_df)

			# load position dataframe
			pos_df = read_positional_file(pos_file)

			well_df = pos_df.merge(well_df, how='inner', on='xy')
			well_df.set_index('acc_xy', inplace=True)
			try:
				well_df.drop('xy', inplace=True, axis=1)
			except:
				pass
			well_df = well_df.T

			X = [float(pos.split('x')[0]) for pos in well_df.columns]
			Y = [float(pos.split('x')[1]) for pos in well_df.columns]
			X = [x if not (str(x).endswith('0')) else round(x, 1) for x in X]
			Y = [y if not (str(y).endswith('0')) else round(y, 1) for y in Y]
			features_corrected = [f"{x[0]}x{x[1]}" for x in zip(X, Y)]
			print(features_corrected[:50])
			well_df.columns = features_corrected
				
			species = [s.split('0', 1)[0] for s in well_df.index.to_list()]
			well_df['species'] = species
			well_df = well_df.groupby('species').sum().T
			#print('WELL DF')
			#print(well_df)

			well_df = well_df.div(well_df.sum(axis=1), axis=0)
			well_df.to_csv('../results/Species_specific_gene_expression/' + well + '_species_specific_transcript_per_feature.tsv',
						   sep='\t')
			print('OUTPUT SHAPE')
			print(well_df.head())
			print(well_df.shape)
			#print('PERCENTAGE')
			#print(well_df_percentage)
			#print(well_df)
			well_df_clipped = well_df.loc[well_df['ENSG']<max(well_df['ENSG'])/4]
			#print('clipped')
			#print(well_df_clipped)
			X_clipped = [float(pos.split('x')[0]) for pos in well_df_clipped.index.to_list()]
			Y_clipped = [float(pos.split('x')[1]) for pos in well_df_clipped.index.to_list()]
			well_df_clip = well_df.loc[well_df['ENSG']>=max(well_df['ENSG'])/4]
			#print(well_df_clip)
			X_clip = [float(pos.split('x')[0]) for pos in well_df_clip.index.to_list()]
			Y_clip = [float(pos.split('x')[1]) for pos in well_df_clip.index.to_list()]
			#print('species: ', well_df.columns.to_list())
			species = well_df.columns.to_list()

			plt.figure(figsize=(len(species)*5,5*2))
			counter = 1
			for s in species:
				plt.subplot(2,len(species),counter)
				plt.scatter(x=X, y=Y*-1, c=well_df[s].to_list(), cmap='jet', s=50, vmin=0, vmax=1)
				plt.title('\n'.join(wrap(s + " transcript expression", 30)))
				plt.colorbar()
				counter += 1
				plt.subplot(2,len(species),counter)
				plt.scatter(x=X_clipped, y=Y_clipped*-1, c=well_df_clipped[s].to_list(), cmap='jet', s=50)
				plt.colorbar()
				plt.scatter(x=X_clip, y=Y_clip*-1, c=well_df_clip[s], cmap='jet', alpha=0.1, s=50)
				plt.title('\n'.join(wrap(s + " expression map; transcripts in 1st quartile", 30)))
				counter += 1
			plt.savefig('../results/Species_specific_gene_expression/' + well + '_proportion_clipped_species_specific_expression_map.png',
			 			bbox_inches='tight')
			#plt.show()
			# plt.clf()


