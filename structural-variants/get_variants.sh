#!/bin/bash

set -e

data_dir="/buffer/ag_abi/evelina/1000genomes/phase2/"
cd $data_dir
work_dir="/group/ag_abi/evelina/DREAM-stellar-benchmark/structural-variants"

# structural variants
#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/integrated_callset/freeze3.sv.alt.vcf.gz

# convert to simplified BED
grep -v "#" freeze3.sv.alt.vcf | awk -F';' '{print $1 "\t" $3 "\t" $7 }' | sed 's/SVLEN=//g' | sed 's/TIG_STRAND=//g' | awk '{if ($9>0) $10="plus"; else $10="minus"; print $1 "\t" $2 "\t" $2 + $9 "\t" 100 "\t" $10 "\t0.0\tchr1\t0\t100"}' > freeze3.sv.alt.bed

if [ -f freeze3.sv.alt.meta.bed ]; then
	rm freeze3.sv.alt.meta.bed
fi

while read chr; do
	awk -v c="$chr" '$1==c' freeze3.sv.alt.bed >> freeze3.sv.alt.meta.bed
done < $work_dir/chrs_in_meta.txt
rm freeze3.sv.alt.bed

