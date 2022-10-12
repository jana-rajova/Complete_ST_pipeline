import scanpy as sc
import numpy as np
import re
import umap
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from matplotlib import cm
from scipy.interpolate import griddata 
# import ST_matrices_to_anndata as Utils
import scanpy as sc
import matplotlib.colors as colors
from matplotlib.pyplot import figure
figure(figsize=(6,5))

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=str, help='choose folder with ST files', default='/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap/')
    parser.add_argument('-f', '--file', type=str, default="ST_files.txt")
    parser.add_argument('-o', '--output', type=str, default='STARsolo_reanalysis_Multimapper/')
    parser.add_argument('--orig', type=int, default=1)
    args = parser.parse_args()

    well_list = []
    os.chdir(args.folder)
    path_list = '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/data/ST_files/' + args.file
    output_path = '/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/data/transcript_feature_plots_orig/' + args.output
    if not os.path.isdir(output_path):
        os.makedirs(output_path, exist_ok=True)
    with open(path_list, 'r') as well_file:
        for well in well_file:
            well_list.append(well.rstrip())


    wells_df = pd.DataFrame()
    print('files in ST matrix folder:', os.listdir(args.folder))

    print(well_list)

    for well in well_list:
        df_sum = pd.DataFrame(columns=['sum', 'feature', 'well'])
        #print(well)
        df = pd.read_csv(args.folder + well + "_stdata.tsv", index_col=0, header=0, sep="\t")
        print(well, "read")
        if args.orig == True:
            df = df.transpose()
        print(well, "transposed")
        df_sum['sum'] = df.sum(axis=0)
        df_sum['feature'] = df.columns
        df_sum['well'] = well
        # df_sum['sum_log2'] = round(np.log2(df.sum(axis=0))).astype('int')
        #print(df_sum.head())
        wells_df = wells_df.append(df_sum)
        print(well, "appended")
       
if args.orig == True:
    wells_df['X'] = wells_df['feature'].str.extract('^([0-9]+)x').astype('int32')
    wells_df['Y'] = wells_df['feature'].str.extract('x([0-9]+)').astype('int32')*(-1)
else:   
    wells_df['X'] = wells_df['feature'].str.extract('^X([0-9.]+)_').astype('float32')
    wells_df['Y'] = wells_df['feature'].str.extract('_([0-9.]+)$').astype('float32')*(-1)
print('Dataframes loaded')
wells_df['sum'] = wells_df['sum'].astype('int64')
log = [str(well_list)]
# log.append("Minimum number of transcripts in all wells (normal/log):")
# log.append(str(min(wells_df['sum']))+ '\t' + str(min(wells_df['sum_log2'])))
# log.append("Maximum number of transcripts in all wells (normal/log):")
# log.append(str(max(wells_df['sum'])) + '\t' + str(max(wells_df['sum_log2'])))
log.append(" ")
log.append("well \t min \t max \t min-log \t max-log")

step = max(wells_df['sum']) - min(wells_df['sum'])
# step_log = max(wells_df['sum_log2']) - min(wells_df['sum_log2'])
cmap = cm.get_cmap('jet') # Colour map (there are many others)
cols = np.linspace(0, 1, step+1)
# cols_log = np.linspace(0, 1, step_log+1)





for well in wells_df['well'].unique():
    ax = plt.gca()
    ax.set_xlim([0,35])
    ax.set_ylim([-35,0])
    # log_well = well + "\t" + str(min(wells_df[wells_df['well']==well]['sum'])) + '\t' +  str(max(wells_df[wells_df['well']==well]['sum'])) + "\t" + str(min(wells_df[wells_df['well']==well]['sum_log2'])) + '\t' +  str(max(wells_df[wells_df['well']==well]['sum_log2']))
    # log.append(log_well)

    if not os.path.isdir(output_path + "unadjusted/"):
        os.makedirs(output_path + "unadjusted/")
    cs = cmap(cols[wells_df[wells_df['well']==well]['sum']-min(wells_df['sum'])])
    sc = ax.scatter(wells_df[wells_df['well']==well]["X"], wells_df[wells_df['well']==well]["Y"], c=wells_df[wells_df['well']==well]["sum"], cmap=cmap, s=30, edgecolor='None')
    plt.axis('off')
    plt.colorbar(sc, ax=ax)
    plt.show
    plt.savefig(output_path + "unadjusted/" + well + '_transcript_per_feature.png', bbox_inches='tight')
    plt.clf()


    if not os.path.isdir(output_path + "normal/"):
        os.makedirs(output_path + "normal/")
    cs = cmap(cols[wells_df[wells_df['well']==well]['sum']-min(wells_df['sum'])])
    sc = ax.scatter(wells_df[wells_df['well']==well]["X"], wells_df[wells_df['well']==well]["Y"], c=cs, s=30, edgecolor='None')
    plt.axis('off')
    plt.colorbar(sc, ax=ax)
    plt.show
    plt.savefig(output_path + "normal/" + well + '_transcript_per_feature.png', bbox_inches='tight')
    plt.clf()

    # ax = plt.gca()
    # ax.set_xlim([0,35])
    # ax.set_ylim([-35,0])
    # if not os.path.isdir(output_path + "log2/"):
    #     os.makedirs(output_path + "log2/")
    # cs_log = cmap(cols_log[wells_df[wells_df['well']==well]['sum_log2']-min(wells_df['sum_log2'])])
    # sc_log = ax.scatter(wells_df[wells_df['well']==well]["X"], wells_df[wells_df['well']==well]["Y"], c=cs_log, s=30, edgecolor='None')
    # plt.axis('off')
    # plt.colorbar(sc_log, ax=ax)
    # plt.savefig(output_path + "log2/" + well + '_transcript_per_feature_log.png', bbox_inches='tight')
    
    # plt.clf()

print(wells_df.head())
with open(output_path + 'log.txt', 'w+') as log_file:
    for line in log:
        log_file.write(line + '\n')


    



