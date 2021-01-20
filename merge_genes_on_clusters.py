import re
import pandas as pd
import os
import argparse
import scanpy as sc
from datetime import datetime

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--st_folder', type=str, help='choose folder with ST files', default='/media/MNM-NetStorage/OrganizedSeqFiles/ST/Complete_ST_pipeline/data/ST_files/ST_matrix_processed/')
    parser.add_argument('-c', '--cluster_folder', type=str, default='/media/MNM-NetStorage/OrganizedSeqFiles/ST/Complete_ST_pipeline/data/Seurat_clustered_wells/')
    parser.add_argument('-d', '--dataset', type=str, default="CN56")
    #parser.add_argument('-f', '--file', type=str, default="CN56_lim.txt")
    parser.add_argument('-o', '--output', type=str, default='/media/MNM-NetStorage/OrganizedSeqFiles/ST/Complete_ST_pipeline/results/cluster_gene_expression/')
    args = parser.parse_args()

    methods = ['seurat', 'DESC', 'scanorama']
    
    path_cluster = args.cluster_folder + args.dataset + '/'
    output_path = args.output + args.dataset + '/'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)



    for method in methods:
        wells_df = pd.DataFrame()
        well_list = [i for i in os.listdir(path_cluster + method + '_' + args.dataset) if (i.endswith('tsv') and re.search('^[A-Z0-9]+_[C-E][1-2]_', i))]
        print(well_list)
        print('Processing ', method, ' files')
        i = 0
        for well_file in well_list:
            i+=1
            print('Processing file ', i, '/', len(well_list))
            well = re.search('^([A-Z0-9]+_[C-E][1-2])_', well_file).group(1)
            path_cluster_df = path_cluster + method + '_' + args.dataset + '/' + well_file
            df_cluster = pd.read_csv(path_cluster_df, header=0, index_col='feature', sep='\t')
            df_cluster = df_cluster['cluster']
            df_cluster
            print(df_cluster.head())
            df_gene = pd.read_csv(args.st_folder + well + '_stdata.tsv', header=0, index_col=0, sep='\t').transpose()
            df = pd.concat([df_cluster, df_gene], axis=1)
            print(df.shape, "after concatenation")
            df = df.groupby(['cluster']).sum()
            print(df.shape, "groupby cluster")
            df['well'] = well
            well_col = df.pop('well')
            df.insert(0, 'well', well_col)
            print()
            # print(df.head())
            wells_df = wells_df.append(df)
            wells_df.fillna(0, inplace=True)

        # print(wells_df.head())
        wells_df.to_csv(output_path + args.dataset + '_' + method + '_cluster_expression.tsv', sep='\t')
            
        print('File saved as: ', output_path + args.dataset + '_' + method + '_cluster_expression.tsv')
    
    with open(output_path + "log.txt", 'w+') as log:
            log.write(str(datetime.now()))
            log.write('\n')
            log.write(str(methods))
            log.write('\n\n')
            for well_file in well_list:

                log.write(well_file)




    



