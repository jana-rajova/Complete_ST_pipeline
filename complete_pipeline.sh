#!/bin/bash
ENSEMBL_TRANS='T'
MATRIX_ANN='T'
SCANORAMA='T'
UMAP_ANN='T'

echo 'Welcome to the Spatial Transcriptomics visualization pipeline!'

#Argumentsfor"ensembltranslate"script
#therearenosetArgumentssofar

matrix_annonation_exec=0

#Argumentsformatrixannnotation
#MA_FILES automatically looks for the file in the data/ST_files folder
# It is stupid, but it is too late now
# Also, replace numpy with cupy in the matrix annotation if possible
MA_FILES="ST_files_test.txt" # this file is the list of st original files to be processed
TIMESTAMP=$(date +"%Y-%m-%d-%H-%M-%S")
# TIMESTAMP="2020-08-09-18-41-34"
echo $(date +"%Y-%m-%d-%H:%M:%S.")
# echo $TIMESTAMP

#Argumentsforgenefiltering
FOLDER="../results/run_$TIMESTAMP"
MUST_GENES="TH"
OPTIONAL_GENES=" CCK MBP PENK PDYN SLC6R3"
STRICT_SELECTION=0

GF_FILES="../data/ST_files/ST_files_test.txt" #these are the files to process from this step on
GRAPH="True"

#Arguments for scanorama (runs from heavy panorama script)
HVG=500
SCAN_FILES="$FOLDER/ST_files_filter_passed$TIMESTAMP.txt"
echo $SCAN_FILES

#ArgumentsforUMAP_plot
TRESH_UNI=8
TRESH_SEP=3
#DIM=$HVG
CLUST_DIM=30

#Argumentsforclustergeneexpressionanalysis
GENE_NUM=50
EXP_TRESH=3

#if ["$ENSEMBL_TRANS"=="T"]; then
#echo 'Creating chimeric species'
#python bin/ensembl_translate.py
#echo "Combined ENSEMBL to gene name file created"
#fi

echo 'Starting matrix annotation'
if [ $matrix_annonation_exec -eq 1 ]
then
  echo Matrix annotation will be performed!
  python matrix_annotation.py -f $MA_FILES -t $TIMESTAMP
  echo "Matrix annotation finished"
fi
if [ ! -d "../data/scanorama/input_st_files" ]
then
  mkdir ../data/scanorama/input_st_files
fi
cp ../data/ST_files/ST_matrix_processed/*stdata.csv ../data/scanorama/input_st_files/
rename 's/csv/tsv/' ../data/scanorama/input_st_files/*csv
cp ../data/ST_files/ST_matrix_processed/*stdata.csv ../data/scanorama/input_st_files/
sed -i s'/,/\t/'g ../data/scanorama/input_st_files/*tsv

echo "Filtering wells based on genes present"
mkdir $FOLDER
python scanorama/filter_on_gene_presence.py -t $TIMESTAMP -g $MUST_GENES -i $OPTIONAL_GENES --graph $GRAPH -f $GF_FILES -s $STRICT_SELECTION -o $FOLDER

#read -r -p "Do you want to continue with these results?[y/N]" response
#if [["$response"=~^([yY][eE][sS]|[yY])$]]
#then
# CONT=1
#else
#echo"Exitingscript" | exit1
#fi


echo "Starting Scanorama script"
python scanorama/process.py $SCAN_FILES
python scanorama/heavy_panorama.py --hvg $HVG -f $SCAN_FILES -t $TIMESTAMP -o $FOLDER
echo "Scanorama processing finished"

echo "Starting UMAP clustering"
python scanorama/UMAP_plot.py -u $TRESH_UNI -s $TRESH_SEP -c $CLUST_DIM -t $TIMESTAMP -o $FOLDER
echo "UMAP clustering finished"

# echo "Starting heatmap generation"
# python bin/scanorama/cluster_gene_expression_analysis.py -d $HVG -g $GENE_NUM -t $EXP_TRESH --timestamp $TIMESTAMP -o $FOLDER
# echo "PIPELINE FINISHED"
