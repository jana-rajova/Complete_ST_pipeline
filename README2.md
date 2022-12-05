# ST_hESC_DA_analysis

This repository contain code for analysis of hESC derived dopaminergic transplants with ST.

### 'SC_processing' folder 
This folder should be run first to generate reference datasets for SingleR, cell2location deconvolution and for correlation analysis.

'subset_reference_loom_dataset'and 'agg_loom_subset'are independent of each other and can be run first. 
'split_obs_matrix'scripts should be run afterwards and be given path to the h5ad subset dataset with the flag -f/--file. This generates a reference dataset for SingleR.

### 'ST_processing'

This folder contains scripts for all processing related to the ST dataset besides deconvolution. All parts up to script 5 can be run before deconvolution

### 'hESC_tx_analysis'

Folder containing scripts for SingleR, which is needed for script #6 in 'ST_processing' folder. 
This script also requires a matrix and annotation files in tsv format (created by 'split_obs_matrix.py').

