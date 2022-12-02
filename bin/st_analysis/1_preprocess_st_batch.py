from src.utils import *
from src.classes import *
import pandas as pd
import argparse
import os
import re
import pickle
from datetime import date

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--folder', type=str, help='folder containing ST dataframes',
                        default='../data/ST_files/ST_matrix_STARsolo_PetterREFs_ensembl/')
    parser.add_argument('-i', '--image_path', type=str, help='folder containing image data',
                        default='../../Images_rev1/')
    parser.add_argument('-s', '--sample_list', type=str, help='folder containing ST samples to process', default='../data/ST_files/ST_files.txt')
    parser.add_argument('-o', '--output', type=str, help='folder to be used as an output',
                        default='../data/st_data_pre_processed/')
    parser.add_argument('--save_separate', default=1, help='save files for each samples separately')
    # parser.add_argument('--mode', default='client')
    # parser.add_argument('--host', default='127.0.0.1')
    # parser.add_argument('--port', default=37237)
    args = parser.parse_args()

    if args.sample_list != 'None':
        samples = pd.read_csv(args.sample_list, header=None)
        samples = samples[0].to_list()
        print('The following samples will be used')
        print(' '.join(samples))
    else:
        print('All samples in the ST folder will be used')
        samples = [re.search('[STCN]+[0-9]+_[C-E][1-2]', x).group(0) for x in os.listdir(args.folder) if 'stdata' in x]
        print(' '.join(samples))

    sample_list = []
    samples_passed = []
    for sample in samples:
        print(f'Loading sample {sample}')
        adata = STLoader()
        adata.load_stdata(args.folder, sample)
        if type(adata.anndata) != int:
            adata.add_image(args.image_path)
            adata.correct_feature_position(args.image_path)
            adata.filter_by_expression()
            # adata.add_species_content()
            adata = adata.anndata
            sample_list.append(adata)
            samples_passed.append(sample)

    adata_st = ST_Anndata()
    adata_st.concat_anndata(sample_list, samples)

    if 'gene_id' not in adata_st.anndata.var.columns:
        adata_st.transalate_ensembl()

    adata_st.anndata = merge_gene_symbol_duplicates(adata_st.anndata, symbol_column='gene_id')
    adata_st.QC(drop_genes=False)
    adata_st = adata_st.anndata

    # create the output_folder and save the datasets
    today = date.today()
    filename = f"{'_'.join(adata_st.obs.slide.unique().tolist())}_joined_ann_stdata_{today.strftime('%d%m%y')}"
    output = f'{args.output}{filename}'

    if args.save_separate == True:
        os.makedirs(f'{args.output}stdata_h5ad/', exist_ok=True)
        os.makedirs(f'{args.output}stdata_tsv/', exist_ok=True)
        for idx, sample in enumerate(samples):
            adata = adata_st[adata_st.obs['sample'] == sample, :]
            adata.write_h5ad(f'{args.output}stdata_h5ad/{sample}_stdata.h5ad')
            adata_df = pd.DataFrame(adata.X, index=adata.obs['feature'], columns=adata.var_names)
            adata_df.to_csv(f'{args.output}stdata_tsv/{sample}_stdata.tsv', sep='\t')



    adata_st.write_h5ad(f'{output}.h5ad')
    pickle.dump(adata_st, open(f'{output}.pkl', 'wb'))

    print('Preprocessed matrices are in '  + args.output.replace(" ", "\ "))