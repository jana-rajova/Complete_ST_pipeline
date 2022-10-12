import pandas as pd
import os 
import numpy as np
import argparse
import re
from scipy.sparse import coo_matrix

def coord_to_sparse(df10x):
    # code barcode and gene entries so that they can be used as mtx coordinates
    df10x[['barcode', 'ensemblID']] = df10x[['barcode', 'ensemblID']].astype('category')
    df10x['barcode_cc'] = df10x['barcode'].cat.codes
    df10x['ensemblID_cc'] = df10x['ensemblID'].cat.codes

    # create reference dataframes for barcodes and genes 
    bc_df = df10x[['barcode_cc','barcode']]
    bc_df = bc_df.groupby(['barcode']).mean()
    bc_df.reset_index(inplace=True)

    gene_df = df10x[['ensemblID_cc','ensemblID']]
    gene_df = gene_df.groupby(['ensemblID']).mean()
    gene_df.reset_index(inplace=True)

    # create coordinate matrix with coordinates 
    print("Creating coordinate matrix")
    col = df10x['ensemblID_cc'].to_numpy()
    row = df10x['barcode_cc'].to_numpy()
    data = df10x['count'].to_numpy()
    n_row = len(df10x['barcode'].unique())
    n_col = len(df10x['ensemblID'].unique())

    coord_mtx = coo_matrix((data, (row, col)), shape=(n_row, n_col))

    # create a sparse matrix from the coordinate matrix and add barcodes and ensembl IDs as index and columns
    arr_10x_sparse = coord_mtx.toarray()
    df_10x_sparse = pd.DataFrame(data=arr_10x_sparse, index=bc_df['barcode'], columns=gene_df['ensemblID']).transpose()

    return df_10x_sparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, nargs='+', default='GSM3891472_rat45_1')
    args = parser.parse_args()
    print(len(args.folder), "samples loaded")

    for dataset in args.folder:
        #output = re.search('(\S*)final_matrix', dataset).group(1)
        print("Processing sample", dataset)
        output  = dataset + dataset[:-1] + "_sparse_mtx.csv"
        dataset = dataset + "/final_matrix.csv"

        df10x = pd.read_csv(dataset, header=None, index_col=False)
        print('Dataset loaded')
        df10x.columns = ['barcode', 'ensemblID', 'symbol', 'count']

        df_10x_sparse = coord_to_sparse(df10x)
        df_10x_sparse.to_csv(output)
        print("Final matrix dimensions:", df_10x_sparse.shape)
        print("Sparse matrix saved as", output)
        




