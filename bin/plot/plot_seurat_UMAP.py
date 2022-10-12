import pandas as pd
import os
import sys
import matplotlib.pyplot as plt

if __name__ == '__main__':
    #%%%
    folder = '../../../results/Batch_corrections/seurat/TX/1/'
    comb_file = [x for x in os.listdir(folder) if 'combined' in x and x.endswith('.tsv')]
    if len(comb_file) == 1:
        comb_file = comb_file[0]
    else:
        print('No file or conflict with multiple files!')
        sys.exit()
    
    comb_file_df = pd.read_csv(folder + comb_file, header=0, index_col=None, sep='\t')
    comb_file_df['animal'] = [x.split('_')[0] for x in comb_file_df['well']]

    #%% block 0
    plt.figure(figsize=(5, 5))
    cluster = plt.scatter(comb_file_df['umap1'], comb_file_df['umap2'], c=comb_file_df['cluster'], s=1, cmap='Spectral')
    # produce a legend with the unique colors from the scatter
    plt.legend(*cluster.legend_elements(num=None),
               loc='center right',
               bbox_to_anchor=(1.1, 0.5),
               frameon=False,
               framealpha=0.0,
               handletextpad=0.5,
               labelspacing=0.00,
               borderpad=0.0,
               handlelength=0,
               borderaxespad=0
               )
    plt.ylabel('UMAP2')
    plt.xlabel('UMAP1')
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    # plt.show()
    plt.savefig(folder + 'seurat_cluster_UMAP_combined.pdf', bbox_inches='tight', dpi=400)


    #%% block1

    well_cat = comb_file_df['well'].astype('category').cat.codes
    dict_wells = dict(zip(comb_file_df['well'].astype('category').cat.codes, comb_file_df['well']))
    wells_ordered = [dict_wells[x] for x in range(min(dict_wells.keys()), max(dict_wells.keys())+1)]
    plt.figure(figsize=(5, 5))
    well = plt.scatter(comb_file_df['umap1'], comb_file_df['umap2'], c=well_cat, s=1, cmap='Spectral')
    plt.legend(handles=well.legend_elements(num=None)[0],
               labels=wells_ordered,
               loc='center right',
               bbox_to_anchor=(1.25, 0.5),
               frameon=False,
               framealpha=0.0,
               handletextpad=0.5,
               labelspacing=0.00,
               borderpad=0.0,
               handlelength=0,
               borderaxespad=0
               )
    plt.xticks([])
    plt.yticks([])
    plt.ylabel('UMAP2')
    plt.xlabel('UMAP1')
    plt.tight_layout()
    # plt.show()
    plt.savefig(folder + 'seurat_well_UMAP_combined.pdf', bbox_inches='tight', dpi=400)

    #%% block 2
    animal_cat = comb_file_df['animal'].astype('category').cat.codes
    dict_animal = dict(zip(comb_file_df['animal'].astype('category').cat.codes, comb_file_df['animal']))
    animal_ordered = [dict_animal[x] for x in range(min(dict_animal.keys()), max(dict_animal.keys())+1)]
    plt.figure(figsize=(5, 5))
    animal = plt.scatter(comb_file_df['umap1'], comb_file_df['umap2'], c=animal_cat, s=1, cmap='Spectral')
    plt.legend(handles=animal.legend_elements(num=None)[0],
               labels=animal_ordered,
               loc='center right',
               bbox_to_anchor=(1.15, 0.5),
               frameon=False,
               framealpha=0.0,
               handletextpad=0.5,
               labelspacing=0.00,
               borderpad=0.0,
               handlelength=0,
               borderaxespad=0
               )
    plt.ylabel('UMAP2')
    plt.xlabel('UMAP1')
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    # plt.show()
    plt.savefig(folder + 'seurat_animal_UMAP_combined.pdf', bbox_inches='tight', dpi=400)