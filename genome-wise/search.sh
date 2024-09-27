#!/usr/bin/env bash

set -x

valik=/group/ag_abi/evelina/valik/build/bin/valik

ibf_bins=${1}
kmer_size=${2}
query_seg_count=${3}
min_len=${4}
threshold=${5}
timeout=${6}


#ref="/buffer/ag_abi/evelina/mouse/chr1.fa"
#query="/buffer/ag_abi/evelina/fly/dna4.fa"
ref="/buffer/ag_abi/evelina/mouse/ref.fa"
query="/buffer/ag_abi/evelina/fly/query.fa"
ref_meta="meta/mouse_ref_b${ibf_bins}.bin"
index="/dev/shm/genome-wise/mouse_b${ibf_bins}_k${kmer_size}_l${min_len}.index"

prefix="b${ibf_bins}_k${kmer_size}_q${query_seg_count}_l${min_len}_t${threshold}"
out="valik_${prefix}.gff"

er=0.025
threads=16
echo "Splitting reference database"
$valik split $ref --without-parameter-tuning -k $kmer_size --verbose --fpr 0.0005 --out $ref_meta --error-rate $er --pattern $min_len -n $ibf_bins

kmer_cmin=1
kmer_cmax=50
echo "Building IBF"
$valik build --without-parameter-tuning --verbose --kmer $kmer_size --fast --threads $threads --output $index --ref-meta $ref_meta --kmer-count-min $kmer_cmin --kmer-count-max $kmer_cmax

log="valik_manual.time"
rm $out
numMatches=10000
sortThresh=$(($numMatches + 1))
threads=1
echo "Search for local matches"
(timeout $timeout /usr/bin/time -a -o $log -f "%e\t%M\t%x\tvalik-search\tb=${ibf_bins}\tk=#${kmer_size}\tq=${query_seg_count}\tl=${min_len}\ter=${er}\tt=${threshold}" $valik search --without-parameter-tuning --verbose --seg-count $query_seg_count --threshold $threshold --pattern $min_len --split-query --cache-thresholds --numMatches $numMatches --sortThresh $sortThresh --index $index --ref-meta $ref_meta --query $query --error-rate $er --threads $threads --output $out --cart-max-capacity $query_seg_count --max-queued-carts $ibf_bins || touch $out ) > ${timeout}_${prefix}.log

truncate -s -1 $log
echo -ne "\tibf_size=" >> $log
ls -lh $index | awk '{print $5}' >> $log

truncate -s -1 $log
echo -ne "\tmatches=" >> $log
wc -l $out | awk '{print $1}' >> $log

rm $index
#rm *.minimiser
#rm *.header

