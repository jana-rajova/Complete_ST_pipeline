import os
import scanpy as sc
import numpy as np
import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str, help='h5ad file to split into observation and matrix dataset')
    parser.add_argument('-o', '--output', type=str, help = 'output folder', default=1)
    args = parser.parse_args()
    print(args)
    
    if args.output == 0:
    	args.output = os.getcwd()
    else:
        args.output = f'{os.getcwd()}/{"/".join(args.file.split("/")[:-1])}'

    print(args.output)

    assert (args.file.split('/')[-1].split('.')[1]) == 'h5ad'
    file_name = args.file.split('/')[-1].split('.')[0]
    output = f'{args.output}/{file_name}_decomposed_reference/'.replace('//', '/')
    print(output)

    os.makedirs(output, exist_ok=True)
    print(args.file)
    print(os.getcwd())

    h5ad_file = sc.read_h5ad(args.file)
    h5ad_file.obs_names_make_unique()
    h5ad_file.obs.index = [f'{x}_{np.random.choice(range(len(h5ad_file.obs.index)), replace=False)}' for x in h5ad_file.obs.index]
    print('Frequency of obs names:')
    print(h5ad_file.obs_names.value_counts().value_counts())


    mat = pd.DataFrame(h5ad_file.X, columns=h5ad_file.var.index, index=h5ad_file.obs.index)
    mat.to_csv(f'{output}/{file_name}_cnt_matrix.csv')
    h5ad_file.obs.to_csv(f'{output}/{file_name}_mta_annotation.csv')

    print(f'Matrix and annotation files are now in {output}')
