import math
import scanpy as sc
import anndata as an
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import os

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


def load_wells(wells_path, wells):
    well_dict = {}
    for well in wells:
        well_dict[well] = pd.read_csv(wells_path + well + '_stdata.tsv',
                                      header=0, index_col=0, sep='\t').transpose()
        well_dict[well] = an.AnnData(X=well_dict[well].to_numpy(),
                                     obs=pd.DataFrame(well_dict[well].index, columns=['feature']),
                                     var=pd.DataFrame(well_dict[well].columns, columns=['symbol']))

        well_dict[well].obs.index = well_dict[well].obs.feature
        well_dict[well].obs['in_tissue'] = 1
        well_dict[well].obs['library_id'] = well
        well_dict[well].obs['slide_id'] = well.split('_')[0]
        well_dict[well].obs['array_row'] = np.array(
            [x[1:].split('_')[0] for x in well_dict[well].obs['feature']]).astype('float')
        well_dict[well].obs['array_col'] = np.array(
            [x[1:].split('_')[1] for x in well_dict[well].obs['feature']]).astype(
            'float')
        well_dict[well].obsm['spatial'] = np.array(
            [x[1:].split('_') for x in well_dict[well].obs['feature']]).astype(
            'float')
        well_dict[well].var.index = well_dict[well].var.symbol
        well_dict[well].uns['spatial'] = {well: well}
        well_dict[well].uns['spatial'][well] = {well: well}

    return well_dict


if __name__ == '__main__':
    wells_path = 'data/ST_files/ST_matrix_STARsolo_correct_PetterREFs_Multimap/stdata/'
    well_file = 'data/ST_files/TX_good_quality.txt'
    dataset = 'TX'
    output = '../results/Batch_corrections/non_corrected/'
    resolution = 1

    output = output + dataset + '/' + str(resolution) + '/'
    os.makedirs(output, exist_ok=True)
    wells = []
    with open(well_file, 'r') as f:
        for line in f:
            wells.append(line.rstrip())

    well_dict = load_wells(wells_path, wells)

    for adata in well_dict.values():
        adata.var_names_make_unique()
        sc.pp.calculate_qc_metrics(adata, inplace=True)

    # for name, adata in well_dict.items():
    #     fig, axs = plt.subplots(1, 4, figsize=(12, 3))
    #     fig.suptitle(f"Covariates for filtering: {name}")
    #
    #     sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    #     sns.distplot(
    #         adata.obs["total_counts"][adata.obs["total_counts"] < 20000],
    #         kde=False,
    #         bins=40,
    #         ax=axs[1]
    #     )
    #     sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
    #     sns.distplot(
    #         adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    #         kde=False,
    #         bins=60,
    #         ax=axs[3]
    #     )
    # plt.show()

    for adata in well_dict.values():
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000, inplace=True)

    adata_spatial = an.concat(
        well_dict,
        label='batch',
        uns_merge="unique"
    )

    sc.pp.neighbors(adata_spatial)
    sc.tl.umap(adata_spatial)
    sc.tl.leiden(adata_spatial, key_added="cluster", resolution=1)

    sc.pl.umap(
        adata_spatial, color=["cluster", "library_id", "slide_id"],
        palette='Spectral', return_fig=True, wspace=0.2
    )

    plt.savefig(output + dataset + '_non-corrected_UMAP.pdf', bbox_inches='tight', dpi=400)
    plt.show()
    plt.clf()
    print(adata_spatial.uns.keys())
    cluster_colors = dict(
        zip([str(i) for i in range(18)], adata_spatial.uns["cluster_colors"])
    )

    for well, adata in well_dict.items():
        print(well)
        X = np.array([x[1:].split('_') for x in adata.obs.feature]).astype('float')
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, key_added="cluster", resolution=1, directed=False)

        sc.pl.umap(
            adata, color=["cluster", "library_id", "slide_id"],
            cmap='Spectral', return_fig=True, wspace=0.2
        )

        plt.savefig(output + well + '_non_corrected_UMAP.pdf', bbox_inches='tight', dpi=400)
        plt.show()
        plt.clf()
        cluster_colors = dict(
            zip([str(i) for i in range(18)], adata.uns["cluster_colors"])
        )

        #adata.obs.columns = ['well' if x == 'library_id' else x for x in adata.obs.columns]
        adata.obs.to_csv(output + well + '_non_corrected_clusters.tsv', sep='\t')

        plt.scatter(X[:, 0], -1*X[:, 1], c=adata.obs.cluster.astype('int'), cmap='Spectral', s=100)
        plt.xlim(0,33)
        plt.title(well)
        plt.ylim(0, -35)
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(output + dataset + '_non_corrected_plot.png', bbox_inches='tight', dpi=500)
        plt.savefig(output + dataset + '_non_corrected_plot.pdf', bbox_inches='tight', dpi=500)




