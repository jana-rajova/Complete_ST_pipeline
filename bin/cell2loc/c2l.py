import scanpy as sc
import os
import re
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
import cell2location
from datetime import date
import argparse

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFsx

def select_slide(adata, s, s_col='sample'):
    r""" This function selects the data for one slide from the spatial anndata object.

    :param adata: Anndata object with multiple spatial experiments
    :param s: name of selected experiment
    :param s_col: column in adata.obs listing experiment name for each location
    """

    slide = adata[adata.obs[s_col].isin([s]), :]
    s_keys = list(slide.uns['spatial'].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]

    slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]}

    return slide

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-st', '--st_h5ad_file', type=str,
                        default='/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/results/st_pp/CN56_ST3_joined_ann_stdata_061022.h5ad')
    parser.add_argument('-sc', '--sc_h5ad_file', type=str,
                        default='/media/Dropbox/MNM projects/Spatial transcriptomics project/Data analysis/Complete_ST_pipeline/data/stereoscope_reference/single_cell_data/L5_061022/L5_CTX_M_STR_cluster_selection_subset_50-500.h5ad')
    parser.add_argument('-o', '--output', type=str,
                        default='../results/cell2location/')
    parser.add_argument('--mode', default='client')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', default=37237)
    args = parser.parse_args()
    
    today = date.today()
    # create paths and names to results folders for reference regression and cell2location models
    results_folder = f'{args.output}{today.strftime("%h-%d%m%y")}/'
    ref_run_name = f'{results_folder}reference_signatures/'
    run_name = f'{results_folder}cell2location_map/'
    print(f'Output folder: {results_folder}')
    
    os.makedirs(f'{ref_run_name}plt/', exist_ok=True)
    

    # load h5ad file
    adata_vis = sc.read_h5ad(args.st_h5ad_file)

    samples = adata_vis.obs['sample'].unique().to_list() # extract samples

    # process sc data
    adata_ref = sc.read_h5ad(args.sc_h5ad_file)

    from cell2location.utils.filtering import filter_genes
    selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

    # filter the object
    adata_ref = adata_ref[:, selected].copy()

    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                            batch_key='ChipID',
                            labels_key='Celltype_assigned',
                            #categorical_covariate_keys=['PCRCycles']
                           )
    # create the regression model
    from cell2location.models import RegressionModel
    mod = RegressionModel(adata_ref)
    # view anndata_setup as a sanity check
    # mod.view_anndata_setup()
    
    mod.train(max_epochs=20, use_gpu=True)
    mod.plot_history(20)
    mod.save(f"{ref_run_name}", overwrite=True)

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_ref = mod.export_posterior(
        adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
    )
    # Save model
    mod.save(f"{ref_run_name}", overwrite=True)
    adata_ref.write_h5ad(f"{ref_run_name}sc.h5ad")
    mod.plot_QC()

    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']

    # find shared genes and subset both spatial and reference
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # add the well transcript median to the observations
    well_transcript_average = []
    for sample in samples:
        np_sample = adata_vis[adata_vis.obs['sample'] == sample, :].X
        well_transcript_average += [int(np.median(np.sum(np_sample, axis=1)))]*np_sample.shape[0]
    adata_vis.obs['well_transcript_median'] = well_transcript_average

    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="slide", continuous_covariate_keys=['well_transcript_median'])

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver,
        N_cells_per_location=50,
        detection_alpha=20)
    # mod.view_anndata_setup()
    
    mod.train(max_epochs=1000,
              batch_size=None,
              train_size=1,
              use_gpu=True)

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history(1000)
    plt.legend(labels=['full data training']);

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
    )
    # Save model
    mod.save(f"{run_name}", overwrite=True)

    # mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

    # Save anndata object with results
    adata_file = f"{run_name}/sp.h5ad"
    adata_vis.write(adata_file)
    mod.plot_QC()

    # plotting part

    cell_types = sorted(adata_ref.obs['Celltype_assigned'].unique().to_list())
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    for sample in samples:
        slide = select_slide(adata_vis, sample)
        with mpl.rc_context({'figure.figsize': [4.5, 5]}):
            sc.pl.spatial(slide, cmap='jet',
                          # show first 8 cell types
                          color=cell_types,
                          ncols=4, size=1.3,
                          img_key='hires',
                          vmin=0,
                          return_fig=True
                         )
            plt.savefig(f'{ref_run_name}plt/{sample}_celltypes_counts.png', bbox_inches='tight', dpi=500)

    cell_types_perc = [f'{x}_feat_perc' for x in cell_types]

    perc = adata_vis.obsm['q05_cell_abundance_w_sf'].div(adata_vis.obsm['q05_cell_abundance_w_sf'].sum(axis=1), axis=0)
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = perc

    for sample in samples:
        slide = select_slide(adata_vis, sample)
        vmax = slide.obs[adata_vis.uns['mod']['factor_names']].max().max()
        with mpl.rc_context({'figure.figsize': [4.5, 5]}):
            sc.pl.spatial(slide, cmap='jet',
                          # show first 8 cell types
                          color=cell_types,
                          ncols=4, size=1.3,
                          img_key='hires',
                          # limit color scale at 99.2% quantile of cell abundance                          vmin=0,
                          return_fig=True
                         )
            plt.savefig(f'{ref_run_name}plt/{sample}_celltypes_feature_percentage.png', bbox_inches='tight', dpi=500)