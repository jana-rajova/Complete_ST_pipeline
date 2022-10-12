import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
import math as m
import matplotlib as mpl
from matplotlib.pyplot import cm
import re
import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable

sns.set_context("paper", rc={"font.size":14})

# execute from Complete_ST_pipeline/bin
def import_gene__features_mtx(folder):
    feature_dict, gene_feature_dict = dict(), dict()

    wells = [x.split('_stdata')[0] for x in os.listdir(folder) if x.endswith('_stdata.tsv')]
    print(wells)

    for idx, well in enumerate(wells):
        print("Processing well ", idx + 1, "/", len(wells), sep='')
        df = pd.read_csv(folder + well + '_stdata.tsv', index_col=0, header=0, sep='\t')

        # chech if features are in columns
        test_index = df.index[np.random.randint(len(df.index)-1)]
        if re.match('^X[0-9.]+_[0-9.]+]', test_index) or re.match('^[0-9.]+x[0-9.]+]', test_index):
            df = df.transpose()

        # check if feature format in x_position'x'y_position
        if df.columns[np.random.randint(len(df.columns)-1)].startswith('X'):
            df.columns = [x[1:].replace('_', 'x') for x in df.columns]
        if 'x' not in df.columns[np.random.randint(len(df.columns)-1)]:
            print("Unrecognized feature format, check your file")
            break

        # fill in the well slot in the dictionaries
        features = df.columns.to_list()
        feature_dict[well] = [x.split('x') for x in features]
        gene_feature_dict[well] = df.transpose()
        gene_feature_dict[well].index = well + '-' + gene_feature_dict[well].index

    return feature_dict, gene_feature_dict


def plot_celltype_contribution(well, stereoscope_df, folder='./', feature_dict={}, scales=['absolute', 'relative']):
    # load features (always in index for stereoscope result files)
    X = stereoscope_df.index.to_list()

    # check for feature format
    if X[0].startswith('X'):
        X = [x.replace('_', 'x')[1:] for x in X]  # keep as list for filter out later tne unused features
    X = np.array([x.split('x') for x in X]).astype('float')  # keep as list for filter out later the unused features

    # load all features of the well
    Xe = feature_dict[well]
    Xe = np.array([x for x in Xe if x not in list(X)]).astype('float')

    celltypes = stereoscope_df.columns.to_list()

    for scale in scales:
        if scale not in ['absolute', 'relative'] or len(scale) == 0:
            print('The only accepted values for scales are "absolute" or "relative" and at least one is mandatory')
            break

        output = folder + "visuals/" + scale + "/"
        output_pdf = folder + "PDF_figures/" + scale + "/"

        for fold in [output, output_pdf]:
            os.makedirs(fold, exist_ok=True)

        ncols = m.ceil(m.sqrt(len(celltypes)))
        nrows = m.ceil((len(celltypes) / ncols))
        figure, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 3, nrows * 4), constrained_layout=True)
        i, j = 0, 0

        for idx, celltype in enumerate(celltypes):
            alphas = stereoscope_df[celltype].to_numpy()
            rgba_colors = np.zeros((len(X), 4))
            # for red the first column needs to be one
            rgba_colors[:, 0] = 1
            # the fourth column needs to be your alphas
            rgba_colors[:, 3] = alphas
            if len(Xe) > 0:
                axes[i, j].scatter(Xe[:, 0], Xe[:, 1] * -1, color='grey', alpha=0.125, s=30)

            if scale == 'absolute':
                vmax = stereoscope_df.max().max()
                vmin = stereoscope_df.min().min()
                norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
                cmap = plt.get_cmap('jet')
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                plt.colorbar(sm, ax=axes[0, ncols - 1])
                # print("vmin:", vmin, "vmax:", vmax)

                axes[i, j].scatter(X[:, 0], X[:, 1] * -1, c=alphas, cmap='jet', s=30, vmin=vmin, vmax=vmax)

            elif scale == 'relative':
                axes[i, j].scatter(X[:, 0], X[:, 1] * -1, c=alphas, cmap='jet', s=30)
            axes[i, j].set_title("\n".join(wrap(celltype, 20)), fontsize=14)
            axes[i, j].set_xlim([0, 34])
            axes[i, j].set_ylim([0, -34])
            axes[i, j].axis('off')
            # plt.axis('off')
            if j < ncols - 1:
                j += 1
            else:
                j = 0
                i += 1
        for plotx in range(idx, ncols * nrows - 1):
            axes[i, j].axis('off')
            j += 1

        plt.savefig(output + well + "_celltype_distr.png", bbox_inches='tight', dpi=500)
        plt.savefig(output_pdf + well  + "_celltype_distr.pdf", bbox_inches='tight')
        #plt.show()
        plt.clf()


def get_combined_cluster_df(cluster_folder):
    cluster_file = [x for x in os.listdir(cluster_folder) if 'combined' in x and '.tsv' in x ]
    if len(cluster_file) != 1:
        print('conflicting "combined" files found! \n '
              'There can be only one combined .tsv file in the selected seurat cluster folder')
        #sys.exit()
    cluster_file = cluster_file[0]
    cluster_df = pd.read_csv(cluster_folder + cluster_file, sep='\t', index_col=None, header=0)
    if cluster_df['feature'][np.random.randint(len(cluster_df['feature'])-1)].startswith('X'):
        cluster_df['feature'] = [x[1:].replace('_', 'x') for x in cluster_df['feature']]

    cluster_df.index = cluster_df['well'] + '-' + cluster_df['feature']

    return cluster_df


def get_human_clusters(cluster_df_combined, human_cutoff=0.1):
    """
    Finds all clusters with human reads content over the "human_cutoff" threshold and singles out the top one
    """
    human_clusters = cluster_df_combined[['cluster', 'human_content']].groupby('cluster').mean().sort_values(
        by='human_content', ascending=False)
    human_clusters = human_clusters[human_clusters['human_content'] > human_cutoff]
    top_human_cluster = human_clusters.loc[human_clusters['human_content'] == human_clusters['human_content'].max()].index.to_list()
    partial_clusters = human_clusters.index.to_list()
    partial_clusters.remove(top_human_cluster[0])


    return {'TX': top_human_cluster, 'Partial': partial_clusters}


def scale_df(df):
    """
    Normalizes the gene expression data so that is in between 0 and 1 for each gene
    """
    df = df.divide(df.max(axis=0))

    return(df)


def merge_stdata_cluster_df(scaled_dict, cluster_df_combined, marker_dict):
    markers_kept = []
    cluster_wells = cluster_df_combined['well'].unique()
    for markers in marker_dict.values():
        markers_kept += markers
    print(markers_kept)
    combined_scale_df= pd.DataFrame()
    for well in cluster_wells:
        print(well)
        for marker in markers_kept:
            if marker not in scaled_dict[well].columns:
                print(marker)
                print(marker in scaled_dict[well].columns)
                scaled_dict[well][marker] = [0]*len(scaled_dict[well].index)
        combined_scale_df = combined_scale_df.append(scaled_dict[well][markers_kept])
        print(combined_scale_df.head())
        print(combined_scale_df.shape)
    cluster_scaled_gene_df = pd.concat([cluster_df_combined, combined_scale_df], axis=1)

    return cluster_scaled_gene_df


def get_region_clusters_score(scaled_df, region_marker_dict, region, animal_mapping, additional_clusters={}, cutoff=0.5):
    marker_list = region_marker_dict[region]

    # create a grouped dataframe but keep the original one for plotting
    region_df = scaled_df.copy()
    region_df_grouped = region_df.groupby('cluster').mean()
    # scale the grouped dataframe to create a cluster score
    region_df_grouped[marker_list] = region_df_grouped[marker_list].subtract(region_df_grouped[marker_list].mean()).divide(region_df_grouped[marker_list].std())
    region_df_grouped['score'] = region_df_grouped[marker_list].mean(axis=1)
    region_df_grouped.reset_index(inplace=True)
    region_df_grouped['region'] = [region if x > cutoff else 'Other' for x in region_df_grouped['score']]

    # add in the additional clusters
    for name, clusters in additional_clusters.items():
        region_df_grouped.loc[region_df_grouped.cluster.isin(clusters), 'region'] = name

    #create mapping for score and 'in' parameters and create columns in the non-grouped dataframe
    passed_dict=dict(zip(region_df_grouped['cluster'], region_df_grouped['region']))
    score_dict=dict(zip(region_df_grouped['cluster'], region_df_grouped['score']))

    region_df['region'] = region_df['cluster'].map(passed_dict)
    region_df['score'] = region_df['cluster'].map(score_dict)
    region_df['slide'] = [x.split('_')[0] for x in region_df['well']]
    region_df['animal'] = region_df['slide'].map(animal_mapping)
    print('ANIMALS')
    print(region_df['animal'].unique())

    return region_df, region_df_grouped


def plot_region_assignment(region_df, region_df_grouped, region='', output='./'):
    output = output + 'common_analysis/'
    s = 30

    region_test_well = region_df['well'].unique()[np.random.randint(len(region_df['well'].unique())-1)]
    region_test_well_df = region_df[region_df['well'] == region_test_well]
    region_test_well_X = np.array([x.split('x') for x in region_test_well_df['feature']]).astype('float')
    region_test_well_df['x'] = region_test_well_X[:, 0]
    region_test_well_df['y'] = region_test_well_X[:, 1]*-1

    plt.clf()
    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(2,1,1)
    sns.barplot(x='cluster', y='score', hue='region', data=region_df_grouped, palette='viridis', ax=ax1)
    plt.title('Striatal score for non-human clusters')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    ax2 = fig.add_subplot(2,2,3)
    im = ax2.scatter(region_test_well_df['x'], region_test_well_df['y'],
                c=region_test_well_df['score'], cmap='jet',s=s)
    plt.axis('off')
    plt.title('cluster scores')
    ax2.set(adjustable='box', aspect='equal')

    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

    ax3 = fig.add_subplot(2,2,4)
    sns.scatterplot('x', 'y', hue='region', palette='rocket', data=region_test_well_df,
                    s=s, edgecolor='ghostwhite', linewidth=0.3, ax=ax3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title('region assignment')
    ax3.set(adjustable='box', aspect='equal')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output + region + "_clusters_selection_scores_well_visualization.pdf", bbox_inches='tight', dpi=500)
    plt.show()
    plt.clf()


def get_singler_data(singler_ref, singler_umap, singler_cutoff=20):
    singler_counts = singler_ref['labels'].value_counts()

    other = singler_counts[singler_counts <= singler_cutoff].sum()
    singler_counts = singler_counts[singler_counts > singler_cutoff]
    singler_counts['Other'] = other

    singler_ref = singler_ref.join(singler_umap, how='inner')

    return singler_counts, singler_ref


def plot_singleR(singler_counts, singler_df, dataset, output):
    print(singler_counts.index)
    celltypes = [x for x in singler_counts.index if x != 'Other']
    print(celltypes)
    umap_all = singler_df[['UMAP_1', 'UMAP_2']]

    cols = 3
    rows = m.ceil(len(celltypes)/cols)
    plt.figure(figsize=(10 * cols, 10 * rows))
    #plt.subplots_adjust(hspace=.001, wspace=0.001)
    s = 20
    for idx, celltype in enumerate(celltypes):
        plt.subplot(rows, cols, idx + 1)
        umap_celltype = singler_df[singler_df['labels'] == celltype][['UMAP_1', 'UMAP_2']]
        plt.scatter(x=umap_all['UMAP_1'], y=umap_all['UMAP_2'],
                    color='ghostwhite', edgecolor='lavender', linewidths=s/40, s=s*0.5, alpha=0.5)
        plt.scatter(x=umap_celltype['UMAP_1'], y=umap_celltype['UMAP_2'],
                    color='purple', edgecolors='lavender', linewidths=0, s=s, alpha=0.6)

        plt.box(False)
        plt.xticks([])
        plt.yticks([])
        plt.xlabel("\n".join(wrap(celltype, 30)), fontsize=25)

    plt.savefig(output + 'SingleR_UMAP_celltype_distribution_' + dataset + '.pdf', dpi=400, bbox_inches='tight')
    plt.savefig(output + 'SingleR_UMAP_celltype_distribution_' + dataset + '.png', dpi=400, bbox_inches='tight')
    plt.tight_layout()
    plt.show()


def merge_celltypes_region_assignment(folder, cluster_df):
    all_celltypes = []
    # select all W.files found in the folder for all wells of all slides
    # (only runs with the most epochs and the newest files are selected)
    stereo_folders = [folder + x + '/' for x in os.listdir(folder) if x!= 'common_analysis']
    stereo_folders = [x + max(os.listdir(x)) + '/' for x in stereo_folders]
    w_files = []
    for stereo_folder in stereo_folders:
        w_files += [stereo_folder + x + '/'  for x in os.listdir(stereo_folder) if '_stdata' in x]
    w_files = [x + max(os.listdir(x)) for x in w_files]

    W_df = pd.DataFrame()
    for w in w_files:
        temp_df = pd.read_csv(w, sep='\t', index_col=0, header=0)
        if temp_df.index[3].startswith('X'):
            temp_df.index = [x[1:].replace('_', 'x') for x in temp_df.index]

        temp_well = re.search('/([C-T]{2}[0-9]+_[C-E][1-2])_stdata/W.', w).group(1)
        temp_df.index = temp_well + '-' + temp_df.index
        W_df = W_df.append(temp_df)

        all_celltypes.append(temp_df.columns.to_list())
    all_celltypes = list(set.intersection(*map(set, all_celltypes)))

    return W_df.join(cluster_df, how='left'), stereo_folders, all_celltypes


def plot_celltype_region_violins(celltype_region_df, all_celltypes, singler_counts, output, omit_region=[]):
    output = output + 'common_analysis/'
    kept_regions = [r for r in celltype_region_df['region'].unique() if r not in omit_region]
    print(kept_regions)
    celltype_region_df = celltype_region_df[celltype_region_df['region'].isin(kept_regions)]
    singler_celltypes = singler_counts.index.to_list()

    celltype_region_df_singler = celltype_region_df.copy()
    non_singler_celltypes = [x for x in all_celltypes if x not in singler_celltypes]
    celltype_region_df_singler['Other'] = celltype_region_df_singler[non_singler_celltypes].sum(axis=1)
    celltype_region_df_singler.drop(non_singler_celltypes, inplace=True, axis=1)

    for df in [celltype_region_df_singler, celltype_region_df]:
        celltypes = [x for x in df.columns if x in (singler_celltypes + all_celltypes)]
        cols = 3*2
        rows = m.ceil(len(celltypes)/3)

        fig, ax = plt.subplots(nrows=rows, ncols=cols,
                               sharey=False,
                               figsize=(3*cols, 4*rows),
                               gridspec_kw={'width_ratios':[4,3,4,3,4,3]})
        r, c = 0, 0

        cmap= cm.get_cmap('Set1')
        handles = [plt.plot([], color=cmap(c), ls="", marker="o")[0] for c in range(len(df['region'].unique()))]
        labels = df['region'].unique()
        print(labels)
        sns.despine(trim=True)
        for idx, ct in enumerate(celltypes):
            sns.violinplot(y=ct, x='animal', ax=ax[r, c],
                           hue='region',
                           cut=True,
                           legend=False,
                           palette='Set1',
                           #scale='area',
                           data=df,
                           scale_hue=True,
                           inner=None)
            sns.stripplot(y=ct, x='animal', ax=ax[r, c],
                          hue='region',
                          dodge=True,
                          color='black',
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
                          size=1,
                          data=df)


            ax[r, c].get_legend().remove()
            ax[r, c].set_ylim(0, 1)
            ax[r, c+1].set_ylim(0, 1)
            #ax[r, c+1].get_legend().remove()
            ax[r, c+1].get_yaxis().set_visible(False)
            ax[r, c+1].spines['left'].set_visible(False)
            ax[r, c].set_ylabel('')
            ax[r, c+1].set_xlabel('Combined')
            #ax[r, c+1].set_frame_on(False)
            ax_ct = fig.add_subplot(rows, int(cols/2), idx+1, frameon = False)
            ax_ct.set_xticks([])
            ax_ct.set_yticks([])
            if ct in singler_celltypes:
                ax_ct.axhline(y=singler_counts[ct]/singler_counts.sum(),
                              xmin=-0.05,
                              xmax=1.005,
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
                print(r, i*-1)
                ax[r, i*-1].set_visible(False)


        ax[0, cols-1].legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left',
                             title='\n'.join(wrap('Feature area/Animal', 20)))
        plt.tight_layout(pad=0.4)
        if len(celltypes) == len(singler_celltypes):
            file_path = output + 'Celltype_distribution_violin_SingleR.pdf'
        elif len(celltypes) > len(singler_celltypes):
            file_path = output + 'Celltype_distribution_violin_complete.png'
        plt.savefig(file_path, bbox_inches='tight', dpi=500)
        plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sf', '--scatterplot_folders', nargs='+', type=str,
                        default=['../../../../Manuscript preparation/stereoscope_relevant_runs/TX_sections/VM_STR_Ctx_reference_final/Lin-selection-VM-STR-Ctx-230821_ST3/75000/',
            '../../../Manuscript preparation/stereoscope_relevant_runs/TX_sections/VM_STR_Ctx_reference_final/Lin-selection-VM-STR-Ctx-230821_CN56/75000/'])
    parser.add_argument('-snd', '--nigral_deconvolution_folders', type=str,
                        default='../../../../Manuscript preparation/stereoscope_relevant_runs/nigral_sections/VM_STR_Ctx_reference_final/')
    parser.add_argument('-strd', '--striatal_deconvolution_folders', type=str,
                        default='../../../../Manuscript preparation/stereoscope_relevant_runs/TX_sections/VM_STR_Ctx_reference_final/')
    parser.add_argument('-st', '--stdata', type=str,
                        default='../data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap/stdata/')
    parser.add_argument('-str', '--cluster_folder_striatum', type=str,
                        default='../../results/Batch_corrections/seurat/TX/1/')
    parser.add_argument('-sn', '--cluster_folder_nigra', type=str,
                        default='../../results/Batch_corrections/seurat/SN/1/')
    parser.add_argument('--singler_umap', type=str,
                        default='../../../../Manuscript preparation/SingleR-analysis/TX_UMAP_reduction_coordinates.tsv')
    parser.add_argument('--singler_full', type=str,
                        default='../../../../Manuscript preparation/SingleR-analysis/Linnarsson_VM-STR-Ctx_23082021_full/')
    parser.add_argument('--singler_selection', type=str,
                        default='../../../../Manuscript preparation/SingleR-analysis/Linnarsson_VM-STR-Ctx_23082021_selection/')
    parser.add_argument("--mode", default='client')
    parser.add_argument("--port", default=62543)
    args = parser.parse_args()

    region_marker_dict= {'STR': ['PENK', 'ADORA2A', 'PPP1R1B'],
                             'SN': ['TH', 'PBX1', 'SLC6A3', 'ALDH1A1', 'DDC', 'RET']}
    animal_mapping = {'CN56': 'Tolerized', 'ST3': 'Nude', 'CN57': 'WT-3', 'CN53': 'WT-2', 'ST1': 'WT-1'}

    feature_dict, gene_feature_dict = import_gene__features_mtx(folder=args.stdata)

    # create a dictionary with scaled gene-feature matrices
    scaled_gene_feature_dict = {}
    for well, st_df in gene_feature_dict.items():
        scaled_gene_feature_dict[well] = scale_df(st_df)

    ## plot stereoscope results
    # for stereo_folder in args.folders:
    #     print(stereo_folder)
    #
    #     # well folders are all folders, that end with '_stdata'
    #     wells = [x.split('_stdata')[0] for x in os.listdir(stereo_folder)
    #              if (os.path.isdir(os.path.join(stereo_folder, x))) and x.endswith('_stdata')]
    #     print(wells)
    #
    #     for idx, well in enumerate(wells):
    #         print("Processing well ", well, ' (', idx + 1, "/", len(wells), ')', sep='')
    #         W_file = max([x for x in os.listdir(stereo_folder + well + '_stdata') if x.startswith('W')])
    #         W_path = stereo_folder + well + '_stdata/' + W_file
    #         W_df = pd.read_csv(W_path, sep='\t', header=0, index_col=0)
    #
    #         plot_celltype_contribution(well=well,
    #                                    folder=stereo_folder,
    #                                    stereoscope_df=W_df,
    #                                    feature_dict=feature_dict,
    #                                    scales=['absolute', 'relative'])


    # Here begins SingleR output analysis

    singler_umap = pd.read_csv(args.singler_umap, sep='\t', header=0, index_col=0)

    VM_STR_Ctx_final_ref_red = pd.read_csv(args.singler_selection + 'SingleR-result_dataframe.tsv', sep='\t')
    VM_STR_Ctx_final_ref_full = pd.read_csv(args.singler_full + 'SingleR-result_dataframe.tsv', sep='\t')

    # the following block will be only plotted and then the variables are rewritten
    singler_counts, singler_df = get_singler_data(VM_STR_Ctx_final_ref_full,
                                                  singler_umap)
    plot_singleR(singler_counts, singler_df,
                 output=args.singler_full,
                 dataset='full')

    # the variables from the folowing block are used for analysis
    singler_counts, singler_df = get_singler_data(VM_STR_Ctx_final_ref_red, singler_umap)
    plot_singleR(singler_counts, singler_df,
                 output=args.singler_selection,
                 dataset='selection')


    # Here begins the analysis for striatal tissue with TX

    str_cluster_df_combined = get_combined_cluster_df(args.cluster_folder_striatum)
    str_cluster_wells = str_cluster_df_combined['well'].unique()

    human_clusters = get_human_clusters(str_cluster_df_combined)

    str_cluster_scaled_gene_df = merge_stdata_cluster_df(scaled_dict=scaled_gene_feature_dict,
                                                     cluster_df_combined=str_cluster_df_combined,
                                                     marker_dict=region_marker_dict)

    # select striatal clusters based on the selected genes
    str_cluster_stdata_df, str_region_df_grouped = get_region_clusters_score(str_cluster_scaled_gene_df,
                                                      region_marker_dict,
                                                      region='STR',
                                                      additional_clusters=human_clusters,
                                                      animal_mapping=animal_mapping)
    plot_region_assignment(str_cluster_stdata_df,
                           str_region_df_grouped,
                           region='STR',
                           output=args.striatal_deconvolution_folders)

    str_celltype_region_df, str_stereoscope_folders, all_ct  = merge_celltypes_region_assignment(args.striatal_deconvolution_folders,
                                                                                                 str_cluster_stdata_df)
    plot_celltype_region_violins(str_celltype_region_df,
                                 all_ct,
                                 singler_counts,
                                 output=args.striatal_deconvolution_folders,
                                 omit_region=['Other', 'Partial'])

    # plot celltype distribution inside the transplant ('TX')
    for stereoscope_folder in str_stereoscope_folders:
        wells = [x.split('_stdata')[0] for x in os.listdir(stereoscope_folder) if '_stdata' in x]
        for well in wells:
            temp_celltype_region_df = str_celltype_region_df[str_celltype_region_df['region'] == 'TX']
            temp_celltype_region_df = temp_celltype_region_df[temp_celltype_region_df['well'] == well]
            temp_celltype_region_df.index = temp_celltype_region_df['feature']
            temp_celltype_region_df = temp_celltype_region_df[all_ct]

            plot_celltype_contribution(well=well,
                                       folder=stereoscope_folder + "selection_visuals/selected_areas_plots/complete/",
                                       stereoscope_df=temp_celltype_region_df,
                                       feature_dict=feature_dict,
                                       scales=['absolute'])

            temp_celltype_region_df = temp_celltype_region_df[[x for x in temp_celltype_region_df.columns if x in singler_counts.index]]

            plot_celltype_contribution(well=well,
                                       folder=stereoscope_folder + "selection_visuals/selected_areas_plots/singler/",
                                       stereoscope_df=temp_celltype_region_df,
                                       feature_dict=feature_dict,
                                       scales=['absolute'])


    # Here begins analysis of SN containing tissues/wells

    sn_cluster_df_combined = get_combined_cluster_df(args.cluster_folder_nigra)
    sn_cluster_wells = sn_cluster_df_combined['well'].unique()

    sn_cluster_scaled_gene_df = merge_stdata_cluster_df(scaled_dict=scaled_gene_feature_dict,
                                                     cluster_df_combined=sn_cluster_df_combined,
                                                     marker_dict=region_marker_dict)

    # select striatal clusters based on the selected genes
    sn_cluster_stdata_df, sn_region_df_grouped = get_region_clusters_score(sn_cluster_scaled_gene_df,
                                                     region_marker_dict,
                                                     region='SN',
                                                     animal_mapping=animal_mapping)
    plot_region_assignment(sn_cluster_stdata_df,
                           sn_region_df_grouped,
                           output=args.nigral_deconvolution_folders)

    # Create violin plots for cell type distribution in regions

    sn_celltype_region_df, sn_stereoscope_folders, all_ct  = merge_celltypes_region_assignment(args.nigral_deconvolution_folders,
                                                                                               sn_cluster_stdata_df)
    plot_celltype_region_violins(sn_celltype_region_df,
                                 all_ct, singler_counts,
                                 #region='SN',
                                 output=args.nigral_deconvolution_folders)

    # plot celltype distribution inside the transplant ('TX')
    for stereoscope_folder in sn_stereoscope_folders:
        wells = [x.split('_stdata')[0] for x in os.listdir(stereoscope_folder) if '_stdata' in x]
        for well in wells:
            temp_celltype_region_df = sn_celltype_region_df[sn_celltype_region_df['region'] == 'SN']
            temp_celltype_region_df = temp_celltype_region_df[temp_celltype_region_df['well'] == well]
            temp_celltype_region_df.index = temp_celltype_region_df['feature']
            temp_celltype_region_df = temp_celltype_region_df[all_ct]

            plot_celltype_contribution(well=well,
                                       folder=stereoscope_folder + "selection_visuals/selected_areas_plots/complete/",
                                       stereoscope_df=temp_celltype_region_df,
                                       feature_dict=feature_dict,
                                       scales=['absolute'])

            temp_celltype_region_df = temp_celltype_region_df[[x for x in temp_celltype_region_df.columns if x in singler_counts.index]]

            plot_celltype_contribution(well=well,
                                       folder=stereoscope_folder + "selection_visuals/selected_areas_plots/singler/",
                                       stereoscope_df=temp_celltype_region_df,
                                       feature_dict=feature_dict,
                                       scales=['absolute'])
