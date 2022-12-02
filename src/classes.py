import os
import re
import anndata as ad
import pandas as pd
import numpy as np
import mygene
import seaborn as sns
# from src.utils import *
# import scanpy as sc
import matplotlib.pyplot as plt

from PIL import Image

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFsx
plt.rcParams['axes.grid'] = False


mg = mygene.MyGeneInfo()

class STLoader:
    def __init__(self):
        self.anndata = ad.AnnData()

    def load_stdata(self, st_path, sample):
        sample_file = [st_path + x for x in os.listdir(st_path) if sample in x and '_stdata' in x]
        if len(sample_file) > 1:
            print('Conflicting sample files!')
            print(sample_file)
        elif len(sample_file) < 1:
            print(f'Sample {sample} not found, skipping!')
            self.anndata = 0
        else:
            # read in the dataframe
            st_df = pd.read_csv(sample_file[0], sep='\t', index_col=0, header=0)
            print('dataset loaded')
            st_df = orient_df(st_df)
            print(st_df.shape)
            # print(st_df)
            unite_feature_notation(st_df)
            # print(st_df)

            # merge duplicate genes
            # st_df['genes'] = st_df.index
            # st_df = st_df.groupby('genes').sum()

            # create the var pard of the anndata
            var = pd.DataFrame(st_df.columns)
            if var.iloc[0, 0].startswith('ENS'):
                var_col = 'ensembl_id'
            else:
                var_col = 'gene_id'
            var.columns = [var_col]
            var[var_col] = var[var_col].str.upper()
            var.index = var[var_col]

            obs = pd.DataFrame(st_df.index)
            obs.columns = ['feature']
            # print(obs)
            if '.' not in ' '.join(obs['feature'].to_list()):
                obs['orig_feature'] = obs['feature']

            print(extract_coordinartes(st_df.index.to_list()))
            coords = extract_coordinartes(st_df.index.to_list())


            obs['array_row'] = coords[:, 0].astype('float')
            obs['array_col'] = coords[:, 1].astype('float')
            obs['in_tissue'] = 1
            obs['sample'] = sample
            obs['slide'] = [x.split('_')[0] for x in obs['sample']]
            obs.index = obs['sample'].astype('str') +  obs['feature'].astype('str')

            X = st_df.iloc[:, :].to_numpy()

            self.anndata = ad.AnnData(obs=obs, var=var, X=X, dtype='float32')
            self.anndata.obs['median_gene_feature'] = np.median(np.sum(self.anndata.X > 0, axis=1))
            self.anndata.obs['median_transcript_feature'] = np.median(np.sum(self.anndata.X, axis=1))
            self.anndata.uns['sample'] = sample
            self.anndata.uns['merged'] = False


    def add_image(self, image_path):
        sample = self.anndata.uns['sample']
        image_path = [image_path + file for file in os.listdir(image_path) if
                      sample in file and '_HE.jpg' in file and not file.startswith('.')][0]
        image = Image.open(image_path)

        self.anndata.uns['spatial'] = {sample: {'images': {'hires': np.asarray(image)}}}
        self.anndata.uns['spatial'][sample]['scalefactors'] = {'spot_diameter_fullres': max(image.size) / 40,
                                                               'tissue_hires_scalef': max(image.size) / max(image.size),
                                                               'fiducial_diameter_fullres': 35 / 2.5}

        X = np.array(self.anndata.obs.array_row.astype('float')) / max(self.anndata.obs.array_row.astype('float')) * \
            image.size[0]
        Y = np.array(self.anndata.obs.array_col.astype('float')) / max(self.anndata.obs.array_col.astype('float')) * \
            image.size[1]

        self.anndata.obsm['spatial'] = np.stack((X, Y), axis=1) * 0.95

    def correct_feature_position(self, corr_coordinates_path):
        if not '.' in self.anndata.obs['feature'][0]:
            file_path = [corr_coordinates_path + x for x in os.listdir(corr_coordinates_path) if
                         self.anndata.uns['sample'] in x and 'data-selection.tsv' in x]
            if len(file_path) > 0:
                file_path = file_path[0]
                corr_coordinates_file = pd.read_csv(file_path, sep='\t', index_col=None, header=0)

                corr_coordinates_file['feature'] = 'X' + corr_coordinates_file['x'].astype('str') + '_' + \
                                                   corr_coordinates_file['y'].astype('str')
                obs = self.anndata.obs.merge(corr_coordinates_file, on='feature', how='inner')
                self.anndata = self.anndata[self.anndata.obs.feature.isin(obs['feature']), :]
                self.anndata.obs = self.anndata.obs.merge(corr_coordinates_file, on='feature', how='inner')
                self.anndata.obs['feature'] = 'X' + self.anndata.obs['new_x'].astype('str') + '_' + \
                                       self.anndata.obs['new_y'].astype('str')
                self.anndata.obs.index = self.anndata.obs['sample'].astype('str') + self.anndata.obs['feature'].astype('str')


    def add_species_content(self, ensembl_path='../data/ST_files/ST_matrix_STARsolo_PetterREFs_ensembl/'):
        '''
        Adds species info to the anndata in the class (by default assumes we are in /bin)
        '''
        # find the ensembl file for the sample
        sample = str(self.anndata.obs['sample'].unique()[0])
        # detect correct ensembl file
        ensembl_file = [ensembl_path + x for x in os.listdir(ensembl_path) if sample in x and 'stdata.tsv' in x]
        assert len(ensembl_file) == 1
        ensembl_file = ensembl_file[0]
        print('Ensembl file found')
        ensembl_df = pd.read_csv(ensembl_file, index_col=0, header=0, sep='\t')

        # make the ensenmbl file compatible
        unite_feature_notation(ensembl_df)
        orient_df(ensembl_df)
        # reduce the ENSEMBL names to their species str
        ensembl_df.columns = [(re.search('[A-Z]+', x).group(0)).split('ENS')[1] for x in ensembl_df.columns]
        ensembl_df = ensembl_df.groupby(ensembl_df.columns.values, axis=1).sum()

        # make them into proportion
        ensembl_df = ensembl_df.divide(ensembl_df.sum(axis=1), axis=0)
        ensembl_df.columns = [f'{species}_content' for species in ensembl_df.columns]
        ensembl_df['orig_feature'] = ensembl_df.index

        # join into the observation dataframe
        self.anndata.obs = self.anndata.obs.merge(ensembl_df, how='left', on='orig_feature')
        self.anndata.obs.index = self.anndata.obs['sample'].astype('str') + '_' + self.anndata.obs['feature'].astype('str')

    def filter_by_expression(self, min_row_count=300, min_feat_count=2, min_features=4):
        print("Filtering features and genes by expression. Positions with less than ", min_row_count,
              " transcripts are discarded as well as genes that do not reach at least ", min_feat_count,
              " in at least ", min_features, " features", sep='')

        # dropping rows with less than min_row_count transcripts
        self.anndata.obs['n_transcripts'] = self.anndata.X.sum(axis=1)
        row_dropped = self.anndata[self.anndata.obs['n_transcripts'] < 300, :].shape[0]

        # dropping genes with less than min_feat_count transcripts in min_features
        self.anndata.var["passed_min_feature_presence"] = (self.anndata.X >= min_feat_count).sum(axis=0) >= min_features
        cols_dropped = (self.anndata.var['passed_min_feature_presence'] == False).sum()
        print(f'{row_dropped} features and {cols_dropped} genes dropped')

        self.anndata = self.anndata[self.anndata.obs['n_transcripts'] >= 300,
                                    self.anndata.var['passed_min_feature_presence'] == True]


class ST_Anndata:
    def __init__(self):
        self.anndata = ad.AnnData()

    def concat_anndata(self, sample_list, samples):
        print(type(sample_list[0]))
        assert type(sample_list[0]) == ad._core.anndata.AnnData
        self.anndata = ad.concat(sample_list,
                                 uns_merge='first',
                                 join='outer',
                                 merge='unique',
                                 keys=samples,
                                 fill_value=0)

    def transalate_ensembl(self):
        if self.anndata.var.index[0].startswith('ENS'):
            self.anndata.var['ensembl_id'] = self.anndata.var.index
            translated = mg.getgenes(self.anndata.var['ensembl_id'].to_list(),
                                     scopes='entrezgene', as_dataframe=True)[['name', 'symbol']]
            translated_no_duplicates = translated.reset_index().drop_duplicates(subset='query')
            translated_no_duplicates.set_index('query', inplace=True)
            print(translated_no_duplicates)
            self.anndata.var = self.anndata.var.join(translated_no_duplicates, how='left')
            self.anndata.var['species'] = [x[0] for x in self.anndata.var['ensembl_id'].str.split('0', 1)]
            self.anndata.var.columns = [x.replace('symbol', 'gene_id') for x in self.anndata.var.columns]
            self.anndata.var['gene_id'] = self.anndata.var['gene_id'].str.upper()
            # remove genes for which there was no symbol found
            self.anndata.var.index = self.anndata.var['gene_id']
            genes_to_keep = [gene for gene in self.anndata.var.index if not isNaN(gene)]
            self.anndata = self.anndata[:, self.anndata.var['gene_id'].isin(genes_to_keep)]


    def QC(self, mt_max=0.2, ribo_max=0.1, drop_features=True, drop_genes=True):
        # sc.pp.calculate_qc_metrics(self.anndata, inplace=True)

        if self.anndata.var.index.name == 'gene_id':
            self.anndata.var['MT'] = self.anndata.var['gene_id'].str.startswith('MT-')
            self.anndata.obs['MT_perc'] = np.sum(self.anndata[:, self.anndata.var['MT'] == True].X, axis=1) / np.sum(self.anndata.X, axis=1)

            self.anndata.var['ribo'] = self.anndata.var['gene_id'].str.startswith(("RPS", "RPL"))
            print(self.anndata.var)
            self.anndata.obs['ribo_perc'] = np.sum(self.anndata[:, self.anndata.var['ribo'] == True].X, axis=1) / np.sum(self.anndata.X, axis=1)

            if drop_features == True:
                dropped_mito = self.anndata[self.anndata.obs['MT_perc'] > mt_max, :].obs['feature'].to_list()
                dropped_ribo = self.anndata[self.anndata.obs['ribo_perc'] > ribo_max, :].obs['feature'].to_list()
                dropped_feat = set(dropped_mito + dropped_ribo)
                print(f'{len(dropped_feat)} features will be dropped due to high mitochondrial/ribosomal gene content ({len(dropped_mito)}/{len(dropped_ribo)})')

                self.anndata = self.anndata[self.anndata.obs['MT_perc'] <= mt_max, :]
                self.anndata = self.anndata[self.anndata.obs['ribo_perc'] <= ribo_max, :]
            if drop_genes == True:
                self.anndata = self.anndata[:, self.anndata.var['MT'] == False]
                self.anndata = self.anndata[:, self.anndata.var['ribo'] == False]


    def select_species(self, species):
        if self.anndata.var['ensembl_id'][0].startswith('ENS'):
            self.anndata.var['species'] = [x[0] for x in self.anndata.var['ensembl_id'].str.split('0', 1)]
        self.anndata = self.anndata[:, self.anndata.var['species'] == species]

    def add_cluster_info(self, cluster_path='/home/jana/Dropbox (Bjorklund Lab)/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/Batch_corrections/seurat/CN56_test/CN56_test_seurat_clusters_combined.tsv'):
        '''
        This function adds the desirable cluster information to the merged anndata.
        I put it here because just like the deconvolution, the result is dependent on the set of wells supplied
        '''
        cluster_file = pd.read_csv(cluster_path, sep='\t', header=0, index_col=None)

        shared_columns = [x for x in cluster_file.columns if x in self.anndata.obs.columns]
        self.anndata.obs = self.anndata.obs.merge(cluster_file, how='left', on=shared_columns)

    def split_and_save(self, output='/home/jana/Dropbox (Bjorklund Lab)/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/st_pp/stdata_h5ad/', by='sample'):
        for cat in self.anndata.obs[by].unique():
            adata = self.anndata[self.anndata.obs[by] == cat, :]
            adata.write_h5ad(f'{output}{cat}_stdata.h5ad')
            adata_df = pd.DataFrame(adata.X, )




def isNaN(string):
    return string != string



def merge_gene_symbol_duplicates(adata, symbol_column='gene_id'):
    adata.var["value"] = 0
    var = adata.var
    original_len = len(var)
    X = adata.X
    temp_df = pd.DataFrame(X).transpose()

    # establish gene names without variants
    gene_no_alt = [x.split('.')[0].upper() for x in var[symbol_column]]

    # add the no alternaticve splicing variants into teh var dataframe and make it the index
    var['Gene_no_alt'] = gene_no_alt
    var.index = var['Gene_no_alt']

    temp_df[symbol_column] = gene_no_alt # add the no splicing variant column into the expression dataframe
    temp_df = temp_df.groupby(symbol_column).sum().transpose()


    var = var.drop_duplicates(subset=['Gene_no_alt'])
    var = var.reindex(temp_df.columns) # make sure that the grouping by some miracle did not rearrange the gene positions

    new_len = len(var)

    adata.uns['merged'] = True
    ad_merge = ad.AnnData(X = temp_df.iloc[:, :].to_numpy(),
                                   var = var,
                                   obs = adata.obs,
                                   obsm = adata.obsm,
                                   uns = adata.uns)
    print(f'Scaled from {original_len} genes incl. alternative splicing to {new_len} genes without alternative splicing variants')

    return ad_merge


def extract_coordinartes(features_list):
    if 'X' in features_list[0]:
        coords = np.array([pos[1:].split('_') for pos in features_list])
    elif 'x' in features_list[0]:
        coords = np.array([pos.split('x') for pos in features_list])
    else:
        print('Coordinates in no known format')
    return coords



def unite_feature_notation(st_df):
    passed = False
    for idx, feat in enumerate([st_df.index, st_df.columns]):
        if re.search('X[0-9.-]+_[0-9.-]+', str(feat[0])):
            corrected_index = [x.replace('-', '') for x in feat]
            passed = idx + 1
        elif re.search('[0-9.-]+x[0-9.-]+', str(feat[0])):
            corrected_index = ['X' + x.replace('x', '_') for x in feat]
            corrected_index = [x.replace('-', '') for x in corrected_index]
            passed = idx + 1

        if passed == True:
            if passed == 1:
                st_df.index = corrected_index
            elif idx == 2:
                st_df.columns = corrected_index
    assert passed > 0


def orient_df(st_df):
    if re.search('X[0-9.-]+_[0-9.-]+', str(st_df.columns[0])) or re.search('[0-9.-]+x[0-9.-]+', str(st_df.columns[0])):
        st_df = st_df.transpose()
    return st_df

def plot_ST(adata, sample, show=True, output=False, feat_max=[34, 32], color='cluster', s=70, vmax_global=True):
    sns.set(rc={'figure.figsize': (10, 10)})
    plt_df = pd.DataFrame()
    plt_df['x'] = adata.obs[adata.obs['sample'] == sample]['array_row'].astype(float).to_numpy()
    plt_df['y'] = adata.obs[adata.obs['sample'] == sample]['array_col'].astype(float).to_numpy()

    print()
    if color in adata.obs.columns:
        ncat = len(adata.obs[color].unique())
        plt_df['c'] = adata.obs[adata.obs['sample'] == sample][color].astype(float).to_list()
        vmax = adata.obs[color].astype(float).max()
        stop = False
    elif color in adata.var.index:
        print(f'{color} in var')
        ncat = len(np.unique(adata[:, color].X))
        plt_df['c'] = adata[adata.obs['sample'] == sample, color].X.flatten().tolist()
        vmax = np.max(adata[:, color].X)
        stop = False
    else:
        print(f'{color} not found!')
        stop = True
    if stop == False:
        if vmax_global:
            vmax = vmax
        else:
            vmax =  max(plt_df['c'])
        print(f'vmax = {vmax}')


        if output or show:
            if 'spatial' in adata.uns.keys():
                im = adata.uns['spatial'][sample]['images']['hires']
                plt_df['x'] = plt_df['x'] * (im.shape[1] * 0.95) / feat_max[1]
                plt_df['y'] = plt_df['y'] * (im.shape[0] * 0.95) / feat_max[0]
                plt.imshow(im)
            if ncat > 50:
                plt.scatter(plt_df['x'], plt_df['y'], c=plt_df['c'], cmap='jet', s=s, label=color,
                            vmin=0, vmax=vmax)
                plt.grid(False)
                plt.axis('off')
                plt.colorbar()
            else:
                for colour in sorted(plt_df['c'].unique()):
                    plt_df_t = plt_df[plt_df['c'] == colour]
                    plt.scatter(plt_df_t['x'], plt_df_t['y'], c=plt_df_t['c'].astype(int), cmap='jet', s=s, label=colour,
                                vmin=plt_df['c'].min(), vmax=vmax)
                    plt.grid(False)
                    plt.axis('off')
                    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
            plt.title(sample)
            if output:
                # print(f'Plots will be saved in folder {output}')
                os.makedirs(f'{output}/', exist_ok=True)
                plt.savefig(f'{output}/{sample}_{color}_feature_plot.png', bbox_inches='tight', dpi=500)
                plt.savefig(f'{output}/{sample}_{color}_feature_plot.pdf', bbox_inches='tight', dpi=500)

            if show == True:
                plt.show()
            plt.clf()


def subsample_anndata(anndata, annot_column='Celltype_assigned', counts=[50, 500]):
    print(f'Dataset will be downsampled to contain between {counts[0]} and {counts[1]} cells per celltype')
    anndata_subset = anndata.copy()
    cells_to_keep = []

    for x in anndata_subset.obs[annot_column].unique():
        print(x)
        all_cells = anndata_subset.obs[anndata_subset.obs[annot_column] == x]['CellID'].to_list()
        if len(all_cells) < counts[0]:
            anndata_subset = anndata_subset[anndata_subset.obs[annot_column] != x, :]
            print(f'{x} with {len(all_cells)} cells will be dropped')
        elif len(all_cells) >= counts[0] and len(all_cells) <= counts[1]:
            cells_to_keep += all_cells
            print(f'All {len(all_cells)} cells will be used')
        elif len(all_cells) > counts[1]:
            cells_that_won_the_lottery = np.random.choice(all_cells, size=counts[1], replace=False).tolist()
            print(f'{len(cells_that_won_the_lottery)} cells will be kept out of {len(all_cells)}')
            cells_to_keep += cells_that_won_the_lottery

    anndata_subset = anndata_subset[anndata_subset.obs['CellID'].isin(cells_to_keep), :]
    print(anndata_subset.obs[annot_column].value_counts())

    return anndata_subset


def convert_loom_to_anndata(loom_file, ca_list=[], ra_list=[], ca_index='CellID', ra_index='Accession'):
    attr_lists = [ca_list, ra_list]

    # if attr lists are empy, keep original columns/rows
    for idx, attr_list in enumerate(attr_lists):
        if len(attr_lists[idx]) == 0:
            if idx == 0:
                attr_lists[idx] = loom_file.ca.keys()
            elif idx == 1:
                attr_lists[idx] = loom_file.ra.keys()

    # select index columns for the dataframes
    attr_indexes = [ca_index, ra_index]
    for idx, index in enumerate(attr_indexes):
        if type(index) == int:
            attr_indexes[idx] = attr_lists[idx][index]
        elif type(index) == str:
            assert index in attr_lists[idx]
    print(f'The indeces for var and obs will be assigned to {attr_indexes[0]} and {attr_indexes[1]}')

    # create var and obs dataframes with defined columns and indexes (indices)
    ad_attr = [pd.DataFrame(), pd.DataFrame()]
    for idx, attr_list in enumerate(attr_lists):
        for attr in attr_list:
            if idx == 0:
                ad_attr[idx][attr] = loom_file.ca[attr]
            elif idx == 1:
                ad_attr[idx][attr] = loom_file.ra[attr]
        ad_attr[idx].index = ad_attr[idx][attr_indexes[idx]]

    adata = ad.AnnData(X=loom_file[:, :].T, var=ad_attr[1], obs=ad_attr[0])

    return adata



