import pickle
from datetime import date
import os
import loompy
import pandas as pd
import anndata as ad
from src.utils import *
# in the future, import argparse and make this work from the command line
# import argparse


if __name__ == '__main__':
    today = date.today()

    input_folder = '../../data/stereoscope_reference/single_cell_data/original/Linnarsson_mouse_brain_all/'

    output_path = '../../data/stereoscope_reference/single_cell_data/'
    dataset = 'L5_selection_STR_CTX_MID'
    output_folder = f'{output_path}L5_{today.strftime("%d%m%y")}/'
    print(output_folder)
    os.makedirs(output_folder, exist_ok=True)

    L5_all = loompy.connect(f'{input_folder}L5_All.loom')

    # choose columns and row attributes you want to keep
    ca_selection = ['Age', 'AnalysisPool', 'CellID', 'Class', 'ClusterName', 'Description', 'Location_based_on',
                    'Mean Reads per Cell', 'Median UMI Counts per Cell', 'Neurotransmitter', 'Region',
                    'TaxonomyRank1', 'TaxonomyRank2', 'TaxonomyRank3', 'TaxonomyRank4', 'bio_celltype', 'ChipID',
                    'SampleID', 'SampleIndex', 'PCRCycles', 'Strain']
    ra_selection = ['Accession', 'Gene', 'Gene_no_alt']

    # create an anndata object with the attributes to keep
    adata_all = convert_loom_to_anndata(L5_all, ca_list=ca_selection, ra_list=ra_selection)

    # # uncomment this section to keep it as a pickle
    # pickle_loc = '../../data/stereoscope_reference/single_cell_data/original/Linnarsson_mouse_brain_all/L5_all_subset_notation.pkl'
    # with open(pickle_loc, 'wb') as pickle_file:
    #     pickle.dump(adata_all, pickle_file)

    # now rename clusters
    # clusters that we want to focus on
    cluster_to_celltype = {
        'MBDOP2': 'Dopaminergic neurons; mouse',
        'MBDOP1': 'Dopaminergic neurons; mouse',
        'MOL1': 'Oligodendrocytes',
        'COP1': 'Oligodendrocytes',
        'MFOL1': 'Oligodendrocytes',
        'MFOL2': 'Oligodendrocytes',
        'MSN1': 'D1 Medium Spiny Neurons; mouse',
        'MSN2': 'D2 Medium Spiny Neurons; mouse',
        'MSN3': 'D2 Medium Spiny Neurons; mouse',
        'MSN4': 'D1 Medium Spiny Neurons; mouse',
        'MSN5': 'D1/D2 Medium Spiny Neurons, striatum',
        'MSN6': 'D1 Medium Spiny Neurons; mouse',
        'TEGLU1': 'Cortical projection neurons; mouse',
        'TEGLU2': 'Cortical projection neurons; mouse',
        'TEGLU3': 'Cortical projection neurons; mouse',
        'TEGLU4': 'Cortical projection neurons; mouse',
        'TEGLU5': 'Cortical projection neurons; mouse',
        'TEGLU6': 'Cortical projection neurons; mouse',
        'TEGLU7': 'Cortical projection neurons; mouse',
        'TEGLU8': 'Cortical projection neurons; mouse',
        'TEGLU9': 'Cortical projection neurons; mouse',
        'TEGLU10': 'Cortical projection neurons; mouse',
        'TEGLU11': 'Cortical projection neurons; mouse',
        'TEGLU12': 'Cortical projection neurons; mouse',
        'TEGLU13': 'Cortical projection neurons; mouse',
        'TEGLU14': 'Cortical projection neurons; mouse',
        'TEGLU15': 'Cortical projection neurons; mouse',
        'TEGLU16': 'Cortical projection neurons; mouse',
        'TEGLU17': 'Cortical projection neurons; mouse',
        'TEGLU18': 'Cortical projection neurons; mouse',
        'TEGLU19': 'Cortical projection neurons; mouse',
        'TEGLU20': 'Cortical projection neurons; mouse',
        'TECHO': 'Cholinergic interneurons; mouse',
        'DECHO1': 'Cholinergic interneurons; mouse',
        'VLMC1': 'Vascular leptomeningeal cells; mouse',
        'VLMC2': 'Vascular leptomeningeal cells; mouse',
        'ABC': 'Vascular leptomeningeal cells; mouse',
        'ACTE1': 'Telencephalon astrocytes, fibrous; mouse',
        'ACTE2': 'Telencephalon astrocytes, protoplasmic; mouse',
        'ACMB': 'Dorsal midbrain Myoc-expressing astrocyte-like; mouse',
        'ACNT1': 'Non-telencephalon astrocytes, protoplasmic; mouse',
        'ACNT2': 'Non-telencephalon astrocytes, fibrous; mouse',
        'VECA': 'Vascular; mouse',
        'VSMCA': 'Vascular; mouse',
        'PER1': 'Vascular; mouse',
        'PER2': 'Vascular; mouse',
        'PER3': 'Vascular; mouse',
        'VECC': 'Vascular; mouse',
        'VECV': 'Vascular; mouse',
        'PVM1': 'Immune cells; mouse',
        'PVM2': 'Immune cells; mouse',
        'MGL3': 'Immune cells; mouse',
        'MGL2': 'Immune cells; mouse',
        'MGL1': 'Immune cells; mouse',
        'RGDG': 'Dentate gyrus radial glia-like cells',
        'RGSZ': 'Subventricular zone radial glia-like cells'

    }

    # supplemented by Taxonomy Rank 4 (or any other column) anotation
    annotated_column, annotation_column = 'ClusterName', 'TaxonomyRank4'
    description_cluster_map = adata_all.obs[[annotated_column, annotation_column]]
    description_dict = description_cluster_map.groupby(annotated_column).first().to_dict()[annotation_column]

    # select the regions you want to keep in the final dataset
    selected_regions = ['Striatum', 'Midbrain', 'Cortex']


    # First replace the clusters you are interested in (manually annotated) and replace the rest with 'TaxonomyRank4'
    adata_all.obs['Celltype_assigned'] = adata_all.obs['ClusterName'].replace(cluster_to_celltype).replace(description_dict)

    #rename the gene_no_alt to symbol in the index
    adata_all.var.index = adata_all.var['Gene_no_alt']
    adata_all.var.index.name = 'symbol'


    adata_selected = adata_all[adata_all.obs['Region'].str.contains('|'.join(selected_regions)) |
                               adata_all.obs['ClusterName'].str.contains('|'.join(celltypes_to_keep_dict.keys())), :]

    adata_selected = merge_gene_symbol_duplicates(adata_selected, symbol_column='Gene_no_alt')
    adata_selected.var.index = adata_selected.var['Gene_no_alt']
    adata_selected.var.index.name = 'symbol'
    adata_selected.obs['Sample'] = 'Linnarsson'

    export_name = f'{output_folder}{dataset}'

    #pickle.dump(adata_selected, open(f'{export_name}.pkl', 'wb'))
    adata_selected.write_h5ad(f'{export_name}.h5ad')
