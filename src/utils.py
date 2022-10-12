# import os
# import re
# import anndata as ad
# import pandas as pd
# from src.classes import *
# import matplotlib.pyplot as plt
#
# from matplotlib import rcParams
# rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFsx
# rcParams['axes.grid'] = False
#
#
# mg = mygene.MyGeneInfo()
#
# def isNaN(string):
#     return string != string
#
#
#
# def merge_gene_symbol_duplicates (adata, symbol_column='gene_id'):
#     var = adata.var
#     original_len = len(var)
#     X = adata.X
#     temp_df = pd.DataFrame(X).transpose()
#
#     # establish gene names without variants
#     gene_no_alt = [x.split('.')[0].upper() for x in var[symbol_column]]
#
#     # add the no alternaticve splicing variants into teh var dataframe and make it the index
#     var['Gene_no_alt'] = gene_no_alt
#     var.index = gene_no_alt
#     temp_df[symbol_column] = gene_no_alt # add the no splicing variant column into the expression dataframe
#     temp_df = temp_df.groupby(symbol_column).sum().transpose()
#     var = var.drop_duplicates(subset=['Gene_no_alt'])
#     var = var.reindex(temp_df.columns) # make sure that the grouping by some miracle did not rearrange the gene positions
#     new_len = len(var)
#     adata.uns['merged'] = True
#     ad_merge = ad.AnnData(X = temp_df.iloc[:, :].to_numpy(),
#                                    var = var,
#                                    obs = adata.obs,
#                                    obsm = adata.obsm,
#                                    uns = adata.uns)
#     print(f'Scaled from {original_len} genes incl. alternative splicing to {new_len} genes without alternative splicing variants')
#
#     return ad_merge
#
#
# def extract_coordinartes(features_list):
#     if 'X' in features_list[0]:
#         coords = np.array([pos[1:].split('_') for pos in features_list])
#     elif 'x' in st_df.columns[0]:
#         coords = np.array([pos.split('x') for pos in features_list])
#     else:
#         print('Coordinates in no known format')
#     return coords
#
#
#
# def unite_feature_notation(st_df):
#     passed = False
#     for idx, feat in enumerate([st_df.index, st_df.columns]):
#         if re.search('X[0-9.-]+_[0-9.-]+', str(feat[0])):
#             corrected_index = [x.replace('-', '') for x in feat]
#             passed = idx + 1
#         elif re.search('[0-9.-]+x[0-9.-]+', str(feat[0])):
#             corrected_index = ['X' + x.replace('x', '_') for x in feat]
#             corrected_index = [x.replace('-', '') for x in corrected_index]
#             passed = idx + 1
#
#         if passed == True:
#             if passed == 1:
#                 st_df.index = corrected_index
#             elif idx == 2:
#                 st_df.columns = corrected_index
#     assert passed > 0
#
#
# def orient_df(st_df):
#     if re.search('X[0-9.-]+_[0-9.-]+', str(st_df.columns[0])) or re.search('[0-9.-]+x[0-9.-]+', str(st_df.columns[0])):
#         st_df = st_df.transpose()
#
#
# def plot_ST(adata, sample, show=True, output=False, feat_max=[34, 32], color='cluster', s=30):
#     im = adata.uns['spatial'][sample]['images']['hires']
#     x = adata.obs[adata.obs['sample'] == sample]['array_row'].astype(float).to_numpy()
#     x = x * (im.shape[1]*0.95) / feat_max[1]
#     y = adata.obs[adata.obs['sample'] == sample]['array_col'].astype(float).to_numpy()
#     y = y * (im.shape[0]*0.95) / feat_max[0]
#     c = adata.obs[adata.obs['sample'] == sample][color].astype(int).to_list()
#     if output:
#         plt.imshow(im)
#         plt.scatter(x, y * (1), c=c, cmap='jet', s=s)
#         plt.savefig(output, bbox_inches='tight', dpi=500)
#         plt.clf()
#     if show == True:
#         plt.imshow(im)
#         plt.scatter(x, y * (1), c=c, cmap='jet', s=s)
#         plt.show()
#         plt.clf()
#
#
# def subsample_anndata(anndata, annot_column='Celltype_assigned', counts=[50, 500]):
#     print(f'Dataset will be downsampled to contain between {counts[0]} and {counts[1]} cells per celltype')
#     anndata_subset = anndata.copy()
#     cells_to_keep = []
#
#     for x in anndata_subset.obs[annot_column].unique():
#         print(x)
#         all_cells = anndata_subset.obs[anndata_subset.obs[annot_column] == x]['CellID'].to_list()
#         if len(all_cells) < counts[0]:
#             anndata_subset = anndata_subset[anndata_subset.obs[annot_column] != x, :]
#             print(f'{x} with {len(all_cells)} cells will be dropped')
#         elif len(all_cells) >= counts[0] and len(all_cells) <= counts[1]:
#             cells_to_keep += all_cells
#             print(f'All {len(all_cells)} cells will be used')
#         elif len(all_cells) > counts[1]:
#             cells_that_won_the_lottery = np.random.choice(all_cells, size=counts[1], replace=False).tolist()
#             print(f'{len(cells_that_won_the_lottery)} cells will be kept out of {len(all_cells)}')
#             cells_to_keep += cells_that_won_the_lottery
#
#     anndata_subset = anndata_subset[anndata_subset.obs['CellID'].isin(cells_to_keep), :]
#     print(anndata_subset.obs[annot_column].value_counts())
#
#     return anndata_subset
#
#
# def convert_loom_to_anndata(loom_file, ca_list=[], ra_list=[], ca_index='CellID', ra_index='Accession'):
#     attr_lists = [ca_list, ra_list]
#
#     # if attr lists are empy, keep original columns/rows
#     for idx, attr_list in enumerate(attr_lists):
#         if len(attr_lists[idx]) == 0:
#             if idx == 0:
#                 attr_lists[idx] = loom_file.ca.keys()
#             elif idx == 1:
#                 attr_lists[idx] = loom_file.ra.keys()
#
#     # select index columns for the dataframes
#     attr_indexes = [ca_index, ra_index]
#     for idx, index in enumerate(attr_indexes):
#         if type(index) == int:
#             attr_indexes[idx] = attr_lists[idx][index]
#         elif type(index) == str:
#             assert index in attr_lists[idx]
#     print(f'The indeces for var and obs will be assigned to {attr_indexes[0]} and {attr_indexes[1]}')
#
#     # create var and obs dataframes with defined columns and indexes (indices)
#     ad_attr = [pd.DataFrame(), pd.DataFrame()]
#     for idx, attr_list in enumerate(attr_lists):
#         for attr in attr_list:
#             if idx == 0:
#                 ad_attr[idx][attr] = loom_file.ca[attr]
#             elif idx == 1:
#                 ad_attr[idx][attr] = loom_file.ra[attr]
#         ad_attr[idx].index = ad_attr[idx][attr_indexes[idx]]
#
#     adata = ad.AnnData(X=loom_file[:, :].T, var=ad_attr[1], obs=ad_attr[0])
#
#     return adata
#
#
#
