import os
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from textwrap import wrap
import matplotlib as mpl
from matplotlib.pyplot import cm
import argparse
import math as m
from mpl_toolkits.axes_grid1 import make_axes_locatable

sns.set_context("paper")

def plot_singleR(singler_counts, singler_df, dataset, output):
    celltypes = [x for x in singler_counts.index if x != 'Other']
    umap_all = singler_df[['umap1', 'umap2']]

    cols = 3
    rows = m.ceil(len(celltypes)/cols)
    plt.figure(figsize=(10 * cols, 10 * rows))
    s = 20
    for idx, celltype in enumerate(celltypes):
        plt.subplot(rows, cols, idx + 1)
        umap_celltype = singler_df[singler_df['Celltype_SingleR'] == celltype][['umap1', 'umap2']]
        plt.scatter(x=umap_all['umap1'], y=umap_all['umap2'],
                    color='ghostwhite', edgecolor='lavender', linewidths=s/40, s=s*0.5, alpha=0.5)
        plt.scatter(x=umap_celltype['umap1'], y=umap_celltype['umap2'],
                    color='purple', edgecolors='lavender', linewidths=0, s=s, alpha=0.6)

        plt.box(False)
        plt.xticks([])
        plt.yticks([])
        plt.xlabel("\n".join(wrap(celltype, 30)), fontsize=25)

    plt.savefig(output + 'SingleR_UMAP_celltype_distribution_' + dataset + '.pdf', dpi=400, bbox_inches='tight')
    plt.savefig(output + 'SingleR_UMAP_celltype_distribution_' + dataset + '.png', dpi=400, bbox_inches='tight')
    plt.tight_layout()
    plt.clf()
    
    
def plot_celltype_region_violins(celltype_region_df, all_celltypes, singler_counts, output, omit_region=[]):
    kept_regions = [r for r in celltype_region_df['region'].unique() if r not in omit_region]
    celltype_region_df = celltype_region_df[celltype_region_df['region'].isin(kept_regions)]
    singler_celltypes = singler_counts.index.to_list()
    

    celltype_region_df_singler = celltype_region_df.copy()
    non_singler_celltypes = [x for x in all_celltypes if x not in singler_celltypes]
    celltype_region_df_singler['Other'] = celltype_region_df_singler[non_singler_celltypes].sum(axis=1)
    celltype_region_df_singler.drop(non_singler_celltypes, inplace=True, axis=1)

    
    for df in [celltype_region_df_singler, celltype_region_df]:
        celltypes = [x for x in df.columns if x in set((singler_celltypes + all_celltypes))]
        
        # vmax = df[celltypes].max().max()
        cols = 2*2
        rows = m.ceil(len(celltypes)/(cols/2))

        fig, ax = plt.subplots(nrows=rows, ncols=cols,
                               sharey=False,
                               figsize=(3*cols, 4*rows),
                               gridspec_kw={'width_ratios':[4,3]*int((cols/2))})
        r, c = 0, 0

        cmap = cm.get_cmap('Set1')
        handles = [plt.plot([], color=cmap(c), ls="", marker="o")[0] for c in range(len(df['region'].unique()))]
        labels = df['region'].unique()
        sns.despine(trim=False)
        for idx, ct in enumerate(celltypes):
            sns.violinplot(y=ct, x='slide', ax=ax[r, c],
                           hue='region',
                           cut=True,
                           legend=False,
                           palette='Set1',
                           #scale='area',
                           data=df,
                           scale_hue=True,
                           inner=None)
            sns.stripplot(y=ct, x='slide', ax=ax[r, c],
                          hue='region',
                          dodge=True,
                          color='black',
                          alpha=0.05,
                          size=1,
                          data=df)

            sns.violinplot(x='region', y=ct, ax=ax[r, c+1],
                           cut=True,
                           #scale='count',
                           palette='Set2',
                           data=df,
                           inner=None
                         )
            sns.stripplot(x='region', y=ct, ax=ax[r, c+1],
                          jitter=0.1,
                          dodge=True,
                          color='black',
                          alpha=0.05,
                          size=1,
                          data=df)


            ax[r, c].get_legend().remove()
            # ax[r, c].set_ylim(0, vmax)
            # ax[r, c+1].set_ylim(0, vmax)
            ax[r, c+1].get_yaxis().set_visible(False)
            ax[r, c+1].spines['left'].set_visible(False)
            ax[r, c].set_ylabel('')
            ax[r, c+1].set_xlabel('Combined')
            ax_ct = fig.add_subplot(rows, int(cols/2), idx+1, frameon = False)
            ax_ct.set_xticks([])
            ax_ct.set_yticks([])
            if ct in singler_celltypes:
                ax_ct.axhline(y=singler_counts[ct]/singler_counts.sum(),
                              xmin=-0.05,
                              xmax=1.1,
                              c="red",
                              linewidth=2,
                              linestyle='--')
            ax_ct.set_title('\n'.join(wrap(ct, 40)))

            if c+3 >= cols:
                c = 0
                r += 1
            else:
                c = c+2
        # remove unused plots
        if (idx+1)*2 < rows*cols:
            for i in range(1, (rows*cols - (idx+1)*2)+1):
                ax[r, i*-1].set_visible(False)


        ax[0, cols-1].legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left',
                             title='\n'.join(wrap('Feature area/slide', 20)))
        plt.tight_layout(pad=0.4)
        if len(celltypes) == len(singler_celltypes):
            file_path = output + 'Celltype_distribution_violin_SingleR'
        elif len(celltypes) > len(singler_celltypes):
            file_path = output + 'Celltype_distribution_violin_complete'
        for ext in ['.png', '.pdf']:
            plt.savefig(f'{file_path}{ext}', bbox_inches='tight', dpi=500)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--singler_file', type=str, help='SingleR file',
                       default='../results/SingleR/hESC_TX_sc_SingleR_analysis_2000_ast_merged.tsv')
    parser.add_argument('-d', '--deconvolution_file', type=str, 
                        help='single cell deconvolution h5ad file',
                        default='../results/cell2loc_res/L5_CTX_M_STR_description_selection_2000_astro-merge/TX/cell2location_map/sp.h5ad')
    parser.add_argument('-c', '--cluster_h5ad_path', type=str, 
                        help='Expression AnnData object with cluster information added, h5ad file',
                        default='../results/Batch_corrections/seurat/TX/TX_st_adata_cluster.h5ad')
    parser.add_argument('-r', '--cluster_regions_h5ad_path', type=str, 
                        help='h5ad file with clusters connected to regions',
                        default='../results/Batch_corrections/seurat/TX/TX_st_adata_cluster_regions.h5ad')
    args = parser.parse_args()
    region_marker_dict= {'STR': ['PENK', 'ADORA2A', 'PPP1R1B'],
                             'SN': ['TH', 'PBX1', 'SLC6A3', 'ALDH1A1', 'DDC', 'RET']}
                             
                         
output = '/'.join(args.deconvolution_file.split('/')[:-1]) + '/'
os.makedirs(output + '/region_cell_composition_plt/', exist_ok=True)

# load the SingleR annotated hESC dataset
singleR_df = pd.read_csv(args.singler_file, sep='\t', header=0, index_col=0)
singleR_df['Celltype_SingleR'] = singleR_df['Celltype_SingleR'].astype('category')

# create a value count for cell types 
singleR_counts = singleR_df['Celltype_SingleR'].value_counts()
singleR_counts = singleR_counts[singleR_counts>20]

# create folder for the singleR plot
singleR_output_folder = f"{'/'.join(args.singler_file.split('/')[:-1])}/plt/"
os.makedirs(singleR_output_folder, exist_ok=True)
sr_dataset = args.singler_file.split('SingleR_analysis_')[1].replace('.tsv', '')
plot_singleR(singleR_counts, singleR_df, dataset=sr_dataset, output=singleR_output_folder)


# load all other datasets (cluster_region, cluster-gene and deconvolution matrix)
deconv_h5ad = sc.read_h5ad(args.deconvolution_file)
cluster_h5ad = sc.read_h5ad(args.cluster_h5ad_path)
cluster_region_h5ad = sc.read_h5ad(args.cluster_regions_h5ad_path)

# create maping cluster-region dictionary
cluster_region_dict = cluster_region_h5ad.obs[['cluster', 'region']]
cluster_region_dict.set_index('cluster', inplace=True)
cluster_region_dict = cluster_region_dict['region'].to_dict()

# create a 'region' column in the cluster obs dataframe 
cluster_type = type(list(cluster_region_dict.keys())[0])

cluster_h5ad.obs['cluster'] = cluster_h5ad.obs['cluster'].astype(cluster_type)
cluster_h5ad.obs['region'] = cluster_h5ad.obs['cluster'].map(cluster_region_dict)


# create a dataframe for plotting celltype proportions in STR adn TX
strip_df = cluster_h5ad.obs
strip_df.index = strip_df['sample'].astype('str') + '_' + strip_df['feature'].astype('str')
deconv_df = deconv_h5ad.obsm['q05_cell_abundance_w_sf']
deconv_df.columns = [x.replace('q05cell_abundance_w_sf_', '') for x in deconv_df.columns]
deconv_df_perc = deconv_df.div(deconv_df.sum(axis=1), axis=0)
celltypes = deconv_df.columns

strip_df_perc = strip_df.join(deconv_df_perc)

regions = list(strip_df_perc['region'].unique())
# print(strip_df_perc[strip_df_perc['region'].unique())
# export dataframes to bo plotted with R
cols = ['sample', 'feature', 'cluster'] + list(deconv_df_perc.columns)
plot_export_df = strip_df_perc[cols]
plot_export_df['region'] = plot_export_df['cluster'].map(cluster_region_dict)
plot_export_df.to_csv(f'{output}celltype_feature.tsv', sep='\t')

if 'TX' in regions:
    plot_export_TX_df = plot_export_df.copy()
    plot_export_TX_df.loc[(plot_export_TX_df['region'] != 'TX'), list(deconv_df_perc.columns)] = -0.001
    plot_export_TX_df.to_csv(f'{output}celltype_feature_TX_only.tsv', sep='\t')

unique_regions = len(strip_df_perc[strip_df_perc['region'] != 'other']['region'].unique())
if unique_regions > 1:
    strip_df_perc = strip_df_perc[strip_df_perc['region'] != 'other']

genes = ['TH',  'COL1A1', 'RET', 'GFAP']

strip_df_stacked = pd.DataFrame()

for ct in celltypes:
    strip_df_ct = strip_df_perc.copy()
    strip_df_ct['celltype'] = ct
    strip_df_ct['celltype_perc'] = strip_df_perc[ct]
    strip_df_ct[genes] = cluster_h5ad[cluster_h5ad.obs['region'].isin(unique_regions), genes].X
    strip_df_stacked = pd.concat([strip_df_stacked, strip_df_ct])

sns.set(font_scale = 5)
for gene in genes:
    print(f'Plotting distribution of {gene} in regions')
    sns.despine()
    sns.stripplot(data=strip_df_stacked, x='region', y=gene, hue='slide',
                  dodge=True, jitter=0.1, size=1, alpha=0.2, color='black')
    sns.violinplot(data=strip_df_stacked, x='region', y=gene, hue='slide',
                   dodge=True, cut=0, scale='count')
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.savefig(f'{output}{gene}_distribution_regions.png', bbox_inches='tight', dpi=400)
    plt.clf()

# plot with seaborn
sns.set(font_scale = 1)
plot_celltype_region_violins(strip_df_stacked,
                             singler_counts = singleR_counts,
                             all_celltypes=celltypes.to_list(), 
                             output=output + '/region_cell_composition_plt/',
                             omit_region=[])

