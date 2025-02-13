#!/bin/bash

set -e

cd /buffer/ag_abi/evelina/1000genomes/phase2/ftp.sra.ebi.ac.uk/vol1/run
truth="/buffer/ag_abi/evelina/1000genomes/phase2/freeze3.sv.alt.vcf"

min_len=100
er=0.033

for id in "HG00731" "HG00732" "NA19238" "NA19239"; do
	echo $id
		
	true_var_out=$id.true.matches
	true_valik_out=$id.true.valik
	false_valik_out=$id.false.valik.var
	
	# compare preprocessed inversions only
	#sample_inv="${id}_potential_inversions_l${min_len}_e${er}.gff"
	
	# compare all local alignments
	sample_inv="${id}_l${min_len}_e${er}_simple.gff"
	
	if [ -f $true_var_out ]; then
		rm $true_var_out
	fi
	if [ -f $true_valik_out ]; then
		rm $true_valik_out
	fi
	if [ -f $false_valik_out ]; then
		rm $false_valik_out
	fi
	if [ -f $sample_inv ]; then
		rm $sample_inv
	fi

	awk '{print $1}' $sample_inv | sort | uniq > $id.chrs

	true_pos_count=0
	true_count=0
	all_local_count=0
	while read chr; do
		ground_truth_chr_match_count=$(awk -v c=$chr '$1 == c' $truth | wc -l | awk '{print $1}')	
		if [ $ground_truth_chr_match_count -ne 0 ]; then
			echo $chr
			grep -v '#' $truth | awk -v c=$chr '$1 == c {print $2 + 500}' | rev | cut -c4- | rev | sort | uniq > "sample.$chr.ground.truth"
			awk -v c=$chr '$1 == c {print $4 + 500}' $sample_inv | rev | cut -c4- | rev | sort | uniq > "sample.$chr.valik"

			(grep -w -f "sample.$chr.valik" "sample.$chr.ground.truth" > "sample.$chr.true" || \
				touch "sample.$chr.true")

			(grep -w -v -f "sample.$chr.ground.truth" "sample.$chr.valik" >> "$false_valik_out" || \
				touch "$false_valik_out")

			if [ -s sample.$chr.true ]; then
				echo $(cat sample.$chr.true)

				while read pos; do
					(grep "$chr	$pos" $truth >> $true_var_out || true)
					(grep "$chr	Stellar	eps-matches	$pos" $sample_inv >> $true_valik_out || true)
				done < sample.$chr.true

				chr_true_pos_count=$(wc -l sample.$chr.true | awk '{print $1}')	
				true_pos_count=$((true_pos_count+chr_true_pos_count))
			fi
			chr_true_count=$(wc -l sample.$chr.ground.truth | awk '{print $1}')
			true_count=$((true_count+chr_true_count))
		
			chr_all_local_count=$(wc -l sample.$chr.valik | awk '{print $1}')
			all_local_count=$((all_local_count+chr_all_local_count))
			
			rm sample.$chr.true
			rm sample.$chr.valik
			rm sample.$chr.ground.truth
		fi
	done < $id.chrs
	
	echo "true positive count $true_pos_count"
	echo "true count $true_count"
	echo "all local count $all_local_count"
done
