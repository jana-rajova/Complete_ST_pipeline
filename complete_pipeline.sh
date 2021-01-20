#!/bin/bash
ENSEMBL_TRANS=0
MATRIX_ANN=0
FILTER=0
DIMRED=0
CLUSTER_ANN=1

echo 'Welcome to the Spatial Transcriptomics visualization pipeline!'

# Matrix annotation parameters

MA_FILES="../data/ST_files/ST_files_striatum.txt" 
#TIMESTAMP=$(date +"%Y-%m-%d-%H-%M-%S")
TIMESTAMP="2020-11-20-16-30-27"
MA_OUTPUT="../data/ST_files/ST_ann_CN56_reseq-test/"
MA_INPUT_FOLDER="../data/ST_files/original_ST_all/"

echo The timestamp is: $(date +"%Y-%m-%d-%H:%M:%S.")

# Filtering on gene presence parameters

MUST_GENES="TH PENK"
OPTIONAL_GENES=" CCK MBP PDYN SLC6R3"
STRICT_SELECTION=0
FILT_INPUT_FOLDER=$MA_OUTPUT
GRAPH="True"
OUTPUT="../results/"

# Scanorama parameters

if [ $FILTER -eq 1 ]
then
	SCAN_FILES="../results/result_$TIMESTAMP/ST_files_filter_passed_$TIMESTAMP.txt "
else
	SCAN_FILES="../results/result_$TIMESTAMP/ST_files_filter_passed_$TIMESTAMP.txt "
fi
HVG=600

#Scanorama UMAP clustering

THRESH_UNI=18
THRESH_SEP=5
CLUST_DIM=30

# Cluster analysis parameters

WELL_FILE="../results/result_$TIMESTAMP/ST_wells_filter_passed_$TIMESTAMP.txt "
DESC_H5AD="../results/result_$TIMESTAMP/DESC/anndata_DESC_clustered.h5ad"
SCAN_H5AD="../results/result_$TIMESTAMP/scanorama_cluster_UMAP_output/anndata_scanorama_clustr_combined-threshold-$THRESH_UNI.h5ad"
CLUST_GENE="TH"
MATCH_SD=0


# Pipeline execution


if [ $ENSEMBL_TRANS -eq 1 ]
then
	echo 'Creating chimeric species'
	python bin/ensembl_translate.py
	echo "Combined ENSEMBL to gene name file created"
fi

if [ $MATRIX_ANN -eq 1 ]
then
	echo Starting matrix annotation!
	python matrix_annotation.py -f $MA_FILES -t $TIMESTAMP -o $MA_OUTPUT --original_matrix_folder $MA_INPUT_FOLDER
	echo "Matrix annotation finished"
fi

if [ $FILTER -eq 1 ]
then
	echo Filtering based on gene presence will be performed
	python scanorama/filter_on_gene_presence.py -f $MA_FILES -t $TIMESTAMP -g $MUST_GENES -i $OPTIONAL_GENES --graph $GRAPH -s $STRICT_SELECTION --input_folder $FILT_INPUT_FOLDER -o $OUTPUT
	echo Only samples with predefined genes present used!
fi

if [ $DIMRED -eq 1 ]
then
	echo "Starting Scanorama script"
	python scanorama/process.py $SCAN_FILES
	python scanorama/heavy_panorama.py --hvg $HVG -f $SCAN_FILES -t $TIMESTAMP -o $OUTPUT
	echo "Scanorama processing finished"
	echo "Starting UMAP clustering"
	python scanorama/scanorama_cluster_UMAP.py -u $THRESH_UNI -s $THRESH_SEP -t $TIMESTAMP -o $OUTPUT
	echo "UMAP clustering finished"
	echo Starting DESC
	python DESC/DESC_run.py -l $WELL_FILE -i $MA_OUTPUT --top_genes $HVG -t $TIMESTAMP
	echo DESC finished
fi

if [ $CLUSTER_ANN -eq 1 ]
then
	echo Cluster anotation and differential gene expression analysis started
	python Cluster_variable_genes.py -d $DESC_H5AD -s $SCAN_H5AD -t $TIMESTAMP -g $CLUST_GENE -m $MATCH_SD -o $OUTPUT
	echo DGE performed!
fi