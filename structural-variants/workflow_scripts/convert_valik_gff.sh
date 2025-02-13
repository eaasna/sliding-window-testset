#!/bin/bash

set -e

if [[ "$#" -ne 3 ]]; then
	echo "Usage: bash convert_valik_gff.sh <min len> <error rate> <id>"
	exit
fi

l=$1
e=$2

dir="/buffer/ag_abi/evelina/1000genomes/phase2/ftp.sra.ebi.ac.uk/vol1/run"
cd $dir
work_dir="/group/ag_abi/evelina/DREAM-stellar-benchmark/structural-variants"

id=$3
is_female=0
if [ $id = "HG00732" ]; then
	is_female=1
fi
if [ $id = "NA19238" ]; then
	is_female=1
fi

in="${id}_l${l}_e${e}.gff"
out="${id}_l${l}_e${e}_simple.gff"
	
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 - 100 "\t" $5 + 100 "\t" $6 "\t" $7 "\t" $8 "\tchr1;seq2Range=0,100;eValue=0;cigar=;mutations="}' $in > $out.tmp
if [ -f $out ]; then
	rm $out
fi

while read chr; do
	if [ $is_female -eq 1 ]; then
		if [ $chr != "chrY" ]; then
			awk -v c="$chr" '$1==c' $out.tmp >> $out
		fi
	else
		awk -v c="$chr" '$1==c' $out.tmp >> $out
	fi
done < $work_dir/chrs_in_meta.txt
rm $out.tmp

