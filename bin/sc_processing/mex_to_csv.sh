#!/bin/bash

k=features
for i in  GSM3891473_rat45_1b
do
    echo Processing sample  $i
    awk -F "\t" 'BEGIN { OFS = "," }; {print NR,$1,$2}' $i/$k.tsv | sort -t, -k 1b,1 > $i/numbered_$k.csv
    awk -F "\t" 'BEGIN { OFS = "," }; {print NR,$1}' $i/barcodes.tsv | sort -t, -k 1b,1 > $i/numbered_barcodes.csv
    tail -n +4 $i/matrix.mtx | awk -F " " 'BEGIN { OFS = "," }; {print $1,$2,$3}' | sort -t, -k 1b,1 > $i/gene_sorted_matrix.csv
    tail -n +4 $i/matrix.mtx | awk -F " " 'BEGIN { OFS = "," }; {print $1,$2,$3}' | sort -t, -k 2b,2 > $i/barcode_sorted_matrix.csv
    join -t, -1 1 -2 1 $i/numbered_$k.csv $i/gene_sorted_matrix.csv | cut -d, -f 2,3,4,5 | sort -t, -k 3b,3 | join -t, -1 1 -2 3 $i/numbered_barcodes.csv - | cut -d, -f 2,3,4,5 > $i/final_matrix.csv
    rm -f $i/barcode_sorted_matrix.csv $i/gene_sorted_matrix.csv $i/numbered_barcodes.csv $i/numbered_$k.csv
    echo Sample $i processed!
done
