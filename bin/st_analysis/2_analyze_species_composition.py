import scanpy as sc
from src.classes import *
import pandas as pd
import anndata as ad
import argparse
import os
import re
import pickle
from datetime import date

def temp_orig_col(adata):
    adata.obs['row_orig'] = np.round(adata.obs['array_row'], 0).astype('int').astype('str')
    adata.obs['col_orig'] = np.round(adata.obs['array_col'], 0).astype('int').astype('str')
    adata.obs['orig_feature'] = 'X' + adata.obs['row_orig'] + '_' + adata.obs['col_orig']

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--st_preprocessed_folder', type=str, help='folder containing ST data',
                        default='../data/st_data_pre_processed/stdata_h5ad/')
    parser.add_argument('-s', '--sample_list', type=str, help='txt files with samples',
                        default='../data/ST_files/STR_normal.txt')
    parser.add_argument('-e', '--ensembl_path', type=str, help='folder containing image data',
                        default='../data/ST_files/ST_matrix_STARsolo_PetterREFs_ensembl/')
    parser.add_argument('-o', '--overwrite', type=int, help='overwrite old file',
                        default=1)
    parser.add_argument('--mode', default='client')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', default=37237)
    args = parser.parse_args()

    if args.sample_list != 'None':
        samples = pd.read_csv(args.sample_list, header=None)
        samples = samples[0].to_list()
        print('The following samples will be used')
        print(' '.join(samples))
    else:
        print('All samples in the ST folder will be used')
        samples = [re.search('[STCN]+[0-9]+_[C-E][1-2]', x).group(0) for x in os.listdir(args.st_preprocessed_folder) if 'stdata' in x]
        print(' '.join(samples))

    for sample in samples:
        adata = STLoader()
        adata.anndata = sc.read_h5ad(f'{args.st_preprocessed_folder}{sample}_stdata.h5ad')
        if 'orig_feature'not in adata.anndata.obs.columns:
            # adata.anndata.obs['orig_row'] = [x.split('.')[0] for x in adata.anndata.obs['array_row']]
            # adata.anndata.obs['orig_col'] = [x.split('.')[0] for x in adata.anndata.obs['array_col']]
            adata.anndata.obs['orig_feature'] = 'X'+ adata.anndata.obs['x'].astype('str') + '_' + adata.anndata.obs['y'].astype('str')
        print(adata.anndata.obs['orig_feature'][:5])
        adata.add_species_content(ensembl_path=args.ensembl_path)
        print(adata.anndata.obs)

        # save the datasets
        if args.overwrite == True:
            output = f'{args.st_preprocessed_folder}{sample}_stdata.h5ad'
            print(f'file {output} will be overwritten')


        else:
            output = f'{args.st_preprocessed_folder}{sample}_stdata_ens_added'
            print(f'Anndata file will be saved as {output}.h5ad')

        # pickle.dump(adata.anndata, open(f'{output}.pkl', 'wb'))
        adata.anndata.write_h5ad(f'{output}.h5ad')

