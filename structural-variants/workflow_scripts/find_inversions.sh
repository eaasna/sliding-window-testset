#!/bin/bash

if [[ "$#" -ne 3 ]]; then
	echo "Usage: bash find_inversion.sh <local_matches> <min_len> <er>"
	exit 1
fi

work_dir=$(dirname $1)
matches=$(basename $1)
min_len=$2
er=$3

cd $work_dir 
inv_dir="potential_inversions_l${min_len}_e${er}"
mkdir -p $inv_dir

awk '$7 == "+" {print}' $matches > forward_matches.gff
awk '$7 == "-" {print}' $matches > reverse_matches.gff

awk '{print $9}' $matches | awk -F';' '{print $1}' | sort | uniq > qnames

ind=1
while read id; 
do
	#echo -e "qname\t$id"
	#echo -e "ind\t$ind"
	grep $id $matches > $ind.gff
	awk '{print $1}' $ind.gff | sort | uniq > chr_hits
	while read chr;
	do
		found_inversion=0
		if [ "$chr" != "Concatenated" ]; then
			#echo -e "chr\t$chr"
			grep $id forward_matches.gff | awk -v chr_id="$chr" ' $1==chr_id ' > curr_forward.gff
			grep $id reverse_matches.gff | awk -v chr_id="$chr" ' $1==chr_id ' > curr_reverse.gff
			
			forward_alignments=$(wc -l curr_forward.gff | awk '{print $1}')
			reverse_alignments=$(wc -l curr_reverse.gff | awk '{print $1}')

			# pick out reference regions that have multiple query matches			
			forward_ref_region="$ind.dstart.forward.region"
			awk '{print $4}' curr_forward.gff | cut -c1-3 | sort | uniq -c | awk '$1>1 {print $2}' > $forward_ref_region			
			forward_in_ref_region=$(wc -l $forward_ref_region | awk '{print $1}')

			reverse_ref_region="$ind.dstart.reverse.region"
			awk '{print $4}' curr_reverse.gff | cut -c1-3 | sort | uniq -c | awk '$1>1 {print $2}' > $reverse_ref_region
			reverse_in_ref_region=$(wc -l $reverse_ref_region | awk '{print $1}')
			
			# same_ref_region=$((forward_same_ref_region+reverse_same_ref_region))
			# find read that matches on forward and reverse strand of the same chr
			if [ $forward_in_ref_region -ge 1 ] && [ $reverse_in_ref_region -ge 1 ]; then
				grep $id $matches | awk -v chr_id="$chr" ' $1==chr_id '  > $inv_dir/${chr}_${ind}.gff
		
				forward_query_table="$ind.qstart.forward.tsv"
				reverse_query_table="$ind.qstart.reverse.tsv"	

				grep -f $forward_ref_region $inv_dir/${chr}_$ind.gff | awk '$7 == "+" {print $9}' | awk -F';' '{print $2}' | sed 's/seq2Range=//g' | awk -F',' '{print $1}' | sort -n | uniq -c > $forward_query_table
				grep -f $reverse_ref_region $inv_dir/${chr}_$ind.gff | awk '$7 == "-" {print $9}' | awk -F';' '{print $2}' | sed 's/seq2Range=//g' | awk -F',' '{print $1}' | sort -n | uniq -c > $reverse_query_table
				max_forward=$(tail -n 1 $forward_query_table | awk '{print $2}')
				min_forward=$(head -n 1 $forward_query_table | awk '{print $2}')
				max_reverse=$(tail -n 1 $reverse_query_table | awk '{print $2}')
				min_reverse=$(head -n 1 $reverse_query_table | awk '{print $2}')
				forward_range=$(($max_forward-$min_forward))
				reverse_range=$(($max_reverse-$min_reverse))
				if [ $forward_range -le -1 ]; then	
					echo -e "max forward\t$max_forward"
					echo -e "min forward\t$min_forward"
					echo -e "max reverse\t$max_reverse"
					echo -e "min reverse\t$min_reverse"
					echo -e "Forward range\t$forward_range"
					echo -e "Reverse range\t$reverse_range"
					echo -e "Error range should not be negative"
					exit 1
				fi		
				# require multiple adjacent local matches in the same reference region
				if [ $forward_range -ge 20 ] && [ $reverse_range -ge 20 ]; then
					#echo -e "\tForward\t$forward_alignments"
					#echo -e "\tReverse\t$reverse_alignments"
					#echo "Potential inversions for read $id on $chr"
					found_inversion=1
				fi
				rm $forward_query_table
				rm $reverse_query_table
				
				if [ $found_inversion -eq 0 ]; then
					rm $inv_dir/${chr}_$ind.gff
				fi
			fi
			rm $forward_ref_region
			rm $reverse_ref_region
		fi
	done < chr_hits
	rm $ind.gff
	ind=$(($ind+1))
done < qnames

