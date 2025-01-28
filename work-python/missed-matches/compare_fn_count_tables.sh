#!/bin/bash

set -ex 

#data_dir="/group/ag_abi/evelina/DREAM-stellar-benchmark/genome-wise/work/last_0_masked"
#out_dir="default_last_vs_stellar_l50_e3"
data_dir="/group/ag_abi/evelina/DREAM-stellar-benchmark/genome-wise/work/last_m100"
out_dir="m100_last_vs_stellar_l50_e3"
mkdir -p $out_dir

fn_gff="$data_dir/mouse_vs_fly_w1_k1_l1.fn.gff"
all_matches="$data_dir/mouse_vs_fly_w1_k1_l1.tsv "

awk '{print $1 "\t" $9}' $fn_gff  | awk -F';' '{print $1}' | sort | uniq -c | awk '{print $1 "\t" $2 "\t" $3}' > $out_dir/fn_count_table.tsv

grep -v "#" $all_matches | awk '{ if ( $4 < 50 && $3 < 94 ) print }' | awk '{print $2 "\t" $1}' | sort | uniq -c | awk '{print $1 "\t" $2 "\t" $3}' > $out_dir/fp_count_table.tsv
