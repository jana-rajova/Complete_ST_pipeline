import re
import pandas as pd
import os
import argparse
import scanpy as sc
from datetime import datetime

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--st_folder', type=str, help='choose folder with ST files', default='../data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap/stdata/')
    parser.add_argument('-c', '--cluster_folder', nargs='+', type=str, default='../results/Seurat_clustered_wells_STARsolo_corrected_MultiMap_PettersREFs/TX/')
    # parser.add_argument('-d', '--dataset', type=str, default="TX")
    parser.add_argument('-o', '--output', type=str, default='../results/cluster_gene_expression/')
    args = parser.parse_args()

    if re.search('seurat', args.cluster_folder, re.IGNORECASE):
        method = 'seurat'
    elif re.search('scanorama', args.cluster_folder, re.IGNORECASE):
        method = 'scanorama'
    elif re.search('DESC', args.cluster_folder, re.IGNORECASE):
        method = 'DESC'
    else:
        method = 'unknown_method'
    
    path_cluster = args.cluster_folder 
    cluster_subfolders = os.listdir(path_cluster)
    if len(cluster_subfolders) == 0:
        folders = [path_cluster]
    else:
        folders = list()
        for subfolder in cluster_subfolders:
            folders.append(path_cluster + subfolder + '/')
    print(folders)
    for folder in folders:
        output_path = args.output + '/DEG_analysis/'
        sample = '_'.join(folder.split('/')[-3:-1])
        output_path = args.output + '/' + method + '_DEG_analysis/' + sample + '/'
        print(sample)
        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        #for method in methods:
        wells_df = pd.DataFrame()
        well_list = [i for i in os.listdir(folder) if (i.endswith('tsv') and re.search('^[A-Z0-9]+_[C-E][1-2]_', i))]
        print(well_list)
        i = 0
        for well_file in well_list:
            i+=1
            print('Processing file ', i, '/', len(well_list))
            well = re.search('^([A-Z0-9]+_[C-E][1-2])_', well_file).group(1)
            path_cluster_df = folder + '/' + well_file
            df_cluster = pd.read_csv(path_cluster_df, header=0, index_col='feature', sep='\t')

            df_cluster = df_cluster['cluster']
            df_cluster
            print(df_cluster.head())
            df_gene = pd.read_csv(args.st_folder + well + '_stdata.tsv', header=0, index_col=0, sep='\t').transpose()
            df = pd.concat([df_cluster, df_gene], axis=1, join='inner')
            print(df.shape, "after concatenation")
            df = df.groupby(['cluster']).sum()
            print(df.shape, "groupby cluster")
            df['well'] = well
            df_count = pd.DataFrame({'cluster': df_cluster.value_counts().index, 'feature_counts': df_cluster.value_counts()})
            print(df_count)
            df = df_count.merge(df, how='inner', on='cluster')
            print(df.head(10))
            well_col = df.pop('well')
            df.insert(0, 'well', well_col)
            print(df.iloc[1:5, 1:5])
            # print(df.head())
            wells_df = wells_df.append(df)
            wells_df.fillna(0, inplace=True)

        # print(wells_df.head())
        wells_df.to_csv(output_path + 'merged_cluster_expression.tsv', sep='\t')
            
        print('File saved as: ', output_path + 'merged_cluster_expression.tsv')
        
        with open(output_path + "log.txt", 'w+') as log:
                log.write(str(datetime.now()))
                log.write('\n')
                #log.write(str(methods))
                #log.write('\n\n')
                for well_file in well_list:

                    log.write(well_file+ "\n")




    



