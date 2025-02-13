#!/bin/bash

set -e

if [[ "$#" -ne 3 ]]; then
	echo "Usage: bash gather_sample_matches.sh <min len> <error rate> <id>"
	exit
fi

dir="/buffer/ag_abi/evelina/1000genomes/phase2/ftp.sra.ebi.ac.uk/vol1/run"
cd $dir

min_len=$1
er=$2
inv_only=0
id=$3

echo $id
grep  $id ../../../igsr_HGSVC2.tsv | awk '{print $1}' | awk -F'/' '{print $6 "/" $7}' > sample_dirs

# compare all local alignments
sample_inv="${id}_l${min_len}_e${er}.gff"
if [ $inv_only -eq 1 ]; then
	# compare preprocessed inversions only
	sample_inv="${id}_potential_inversions_l${min_len}_e${er}.gff"
fi

if [ -f $sample_inv ]; then
	rm $sample_inv
fi
while read dir; do
	echo $dir
	if [ $inv_only -eq 1 ]; then	
		cat $dir/potential_inversions_l${min_len}_e${er}/*.gff > $sample_inv
	else
		cat $dir/l${min_len}_e${er}.gff >> $sample_inv
	fi
done < sample_dirs
rm sample_dirs
