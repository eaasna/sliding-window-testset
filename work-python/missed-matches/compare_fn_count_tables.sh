#!/bin/bash

set -ex 

data_dir="/group/ag_abi/evelina/DREAM-stellar-benchmark/genome-wise/work/last"
stellar="l50_e2"
min_len=50
min_percid=96

function count_matches() {
	m=$1
	fn_gff="$data_dir/$stellar/mouse_vs_fly_w1_k1_m$m.fn.gff"
	all_matches="$data_dir/mouse_vs_fly_w1_k1_m$m.tsv "
	out_dir="m${m}_last_vs_stellar_$stellar"
	mkdir -p $out_dir

	awk '{print $1 "\t" $9}' $fn_gff  | awk -F';' '{print $1}' | sort | uniq -c | awk '{print $1 "\t" $2 "\t" $3}' > $out_dir/fn_count_table.tsv

	grep -v "#" $all_matches | awk -v l="$min_len" -v p="$min_percid" '{ if ( $4 < l && $3 < p ) print $2 "\t" $1}' | sort | uniq -c | awk '{print $1 "\t" $2 "\t" $3}' > $out_dir/fp_count_table.tsv
}

count_matches 100
count_matches 10
