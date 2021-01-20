from process import load_names, merge_datasets
import os
import scanorama as sc
import pandas as pd
from get_positions import index_from_csv as index 
from time import time
import numpy as np

NAMESPACE = 'panorama'
VERBOSE = 2
PERPLEXITY = 5
PERPLEXITYU = 50



if __name__ == '__main__':

#    from config import data_names
#    datasets, genes_list, n_cells = load_names(data_names)
#    t0 = time()
    #in order to get an idea, what changes if the union is True of False, theUNION loop was created
    #this should create separate folders and create graphs and dataset-dimred files in teh specific UNION-True or UNION-False
    #folders
    for UNION in (False, True):
        
        from config import data_names
        datasets, genes_list, n_cells = load_names(data_names)
        t0 = time() 
        print('UNION = ', UNION)
        datasets_dimred, datasets, genes = sc.correct(
                datasets, genes_list, ds_names=data_names,
                sigma=15, return_dimred=True, return_dense = True, union=UNION
                )
        #lines 23-26 were added for the UNION loop
        dict_pos = index("data/")
        path = "./Union-" + str(UNION)
        print(path)
        print(os.path.isdir(path))
        if os.path.isdir(path)==False:
            os.mkdir(path)
        os.chdir(path)  
#        
        i = 0
        for key in dict_pos.keys():
            print(datasets_dimred[0])
            csv = 'dataset-dimred_' + str(key) + '.csv'
            df = pd.DataFrame(datasets_dimred[i])
            print('csv: ', csv)
            print("key length", len(dict_pos[key]))
            print("dimred length", len(datasets_dimred[i]))
            df.insert(0, "position", dict_pos[key], True)   
            df.to_csv(csv, index = 0, header = False)
            i += 1
            print(df.columns)
        

        if VERBOSE:
            print('Integrated and batch corrected panoramas in {:.3f}s'
                  .format(time() - t0))
    
        labels = []
        names = []
        curr_label = 0
        for i, a in enumerate(datasets):
            labels += list(np.zeros(a.shape[0]) + curr_label)
            names.append(data_names[i])
            curr_label += 1
        labels = np.array(labels, dtype=int)
        
    
        sc.visualize(datasets_dimred, labels, NAMESPACE + '_ds', names,
                              multicore_tsne=False, perplexity = PERPLEXITY,)
    
#         Uncorrected.
        os.chdir("../")
        datasets, genes_list, n_cells = load_names(data_names)
        datasets, genes = merge_datasets(datasets, genes_list)
        datasets_dimred = sc.dimensionality_reduce(datasets)
        print(datasets_dimred)
        labels = []
        names = []
        curr_label = 0
        for i, a in enumerate(datasets):
            labels += list(np.zeros(a.shape[0]) + curr_label)
            names.append(data_names[i])
            curr_label += 1
        labels = np.array(labels, dtype=int)
        
        os.chdir(path)    
        sc.visualize(datasets_dimred, labels,
                              NAMESPACE + '_ds_uncorrected', names, perplexity = PERPLEXITYU)
        
        #added with the UNION loop
        os.chdir("../")
