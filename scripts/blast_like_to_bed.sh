#!/bin/bash

set -ex

if [[ "$#" -ne 2 ]]; then
	echo "Usage: bash blast_like_to_bed.sh <in> <out>"
	exit 1
fi

in=$1
out=$2


# 1. add plus or minus strand
# 2. flip query start and end if on minus strand
# 3. flip ref start and end if on minus strand
grep -v "#" $1 | \
	awk '{ if($7>$8) $5="minus"; else $5="plus"; print $2 "\t" $9 "\t" $10 "\t" $3 "\t" $5 "\t" $11 "\t" $1 "\t" $7 "\t" $8 ; }' | \
	awk '$8>$9{tmp=$8; $8=$9; $9=tmp} 1' | \
	awk 'BEGIN {OFS = "\t" } $5=="minus"{tmp=$2; $2=$3; $3=tmp} 1' > $2

