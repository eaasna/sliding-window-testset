#!/bin/bash

set -e

if [[ "$#" -ne 5 ]]; then
	echo "Usage: bash find_local_matches.sh <ref> <min len> <er> <reads IN> <matches OUT>"
	exit 1
fi

ref=$1
ref_dir="$(dirname "$ref")"

minlen=$2
er=$3
unmapped_reads=$4
matches=$5

s="1111110110110111111"
numMatches=200
sortThresh=$((numMatches+1))

meta="$ref_dir/human_l${minlen}_e${er}_s${s}.bin"
index="$ref_dir/human_l${minlen}_s${s}.index"

if [ ! -f $meta ]; then
	valik split $ref --fpr 0.005 --out $meta --error-rate $er --pattern $minlen -n 4096 --shape $s
fi

if [ ! -f $index ]; then
	valik build --fast --threads 8 --output $index --ref-meta $meta
fi

read_count=$(grep ">" $unmapped_reads | wc -l | awk '{ print $1 }')
seg_count=$((read_count*5))
sample_dir="$(dirname "$unmapped_reads")"
	
log="$sample_dir/search_valik.time"
( /usr/bin/time -a -o $log -f "%e\t%M\t%x\t%C" \
	valik search --keep-best-repeats --bin-entropy-cutoff 0.25 \
		--split-query --index $index --ref-meta $meta \
		--query $unmapped_reads --error-rate $er --threads 16 \
		--output $matches --cart-max-capacity 100 \
		--numMatches $numMatches --sortThresh $sortThresh \
		--without-parameter-tuning --threshold 31 \
		--seg-count $seg_count --max-queued-carts 1024 \
		--pattern $minlen \
		--verbose &> $matches.search.err )

