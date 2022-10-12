import os # todo create the environment for complete pipeline or a docker container you can execute in
import sys
import pandas as pd
import re
import argparse
import numpy as np

def get_Wfiles_stereoscope(path):
    subfold = [path + x + '/' for x in os.listdir(path) if os.path.isdir(path + x)]
    Wfiles = [x + os.listdir(x)[0] for x in subfold if len(os.listdir(x))==1 and os.listdir(x)[0].startswith('W.')]

    return Wfiles

def get_Wfiles_cell2loc(path):
    Wfiles = [path + x for x in os.listdir(path) if x.startswith('W.') and x.endswith('.tsv')]

    return Wfiles


def select_human_cluster(combined_cluster_file):
    human_cluster = combined_cluster_file.copy()
    human_cluster = human_cluster.groupby('cluster').mean()
    human_cluster = human_cluster[human_cluster['human_content'] == max(human_cluster['human_content'])].index.to_list()[0]

    return human_cluster


def create_cluster_well_feature_dict(combined_cluster_file, human_cluster):
    well_cluster_feature_dict = dict()
    for well in combined_cluster_file.well.unique():
        print(well)
        temp = combined_cluster_file[combined_cluster_file['well'] == well]
        temp = temp[temp['cluster'] == human_cluster]
        temp = temp['feature'].to_list()
        if 'x' in temp[0]:
            temp = ['X' + x.replace('x', '_') for x in temp]
        well_cluster_feature_dict[well] = temp

    return well_cluster_feature_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--deconvolution', type=str, default='cell2loc', help='options are "cell2loc" and "stereoscope"')
    parser.add_argument('--cell2loc_folder', type=str, default='../cell2loc_test_data/std_out/')
    parser.add_argument('-s', '--stereoscope_folders', nargs='+', type=str, default=['../../../../../Manuscript preparation/stereoscope_relevant_runs/TX_sections/VM_STR_Ctx_reference_final/Lin-selection-VM-STR-Ctx-230821_ST3/75000/',
            '../../../../../Manuscript preparation/stereoscope_relevant_runs/TX_sections/VM_STR_Ctx_reference_final/Lin-selection-VM-STR-Ctx-230821_CN56/75000/'])
    parser.add_argument('-c', '--cluster_folder', type=str, default='../../results/Batch_corrections/seurat/TX/1/')
    parser.add_argument("--mode", default='client')
    parser.add_argument("--port", default=62543)
    args = parser.parse_args()

    combined_cluster_file = [x for x in os.listdir(args.cluster_folder) if 'clusters_combined' in x and '.tsv' in x and not x.startswith('.')]
    print(combined_cluster_file)
    if len(combined_cluster_file) == 1:
        combined_cluster_file = args.cluster_folder + combined_cluster_file[0]
        combined_cluster_file = pd.read_csv(combined_cluster_file, sep='\t', header=0, index_col=None)
    else:
        print('conflicting cluster assignment files, exiting')
        sys.exit()

    human_cluster = select_human_cluster(combined_cluster_file)

    well_cluster_feature_dict = create_cluster_well_feature_dict(combined_cluster_file, human_cluster)
    if args.deconvolution == 'stereoscope':
        for path in args.stereoscope_folders:
           Wfiles = get_Wfiles_stereoscope(path)
    elif args.deconvolution == 'cell2loc':
        Wfiles = get_Wfiles_cell2loc(args.cell2loc_folder)

    if Wfiles:
        output_folder = args.cluster_folder + 'selection_TX_dataframes-' + args.deconvolution + '/'
        os.makedirs(output_folder, exist_ok=True)
       for Wfile in Wfiles:
            well = re.search('([A-Z]{2}\d+_[C-E][1-2])', Wfile).group(1)
            if well in well_cluster_feature_dict.keys():
                print(f'Processing sample {well}')
                Wfile_mod = pd.read_csv(Wfile, header=0, index_col=0, sep='\t')
                non_TX_features = [x for x in Wfile_mod.index if x not in well_cluster_feature_dict[well]]
                Wfile_mod.loc[non_TX_features] = -0.0001
                #Wfile_mod = Wfile_mod.loc[well_cluster_feature_dict[well]]

                Wfile_mod.to_csv(output_folder + well + '_only_selection.tsv',
                                sep='\t')
            else:
                print(f'Skipping sample {well}')




