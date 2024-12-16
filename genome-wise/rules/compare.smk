rule valik_compare_blast:
	input:
		ref_meta = dream_out + "/meta/b" + str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_len[0]) + "_e" + str(errors[0]) + ".bin",
		test_files = expand(dream_out + "/b{b}_fpr{fpr}_l{min_len}_cmin{cmin}_cmax{cmax}_e{er}_ent{bin_ent}_cap{max_cap}_carts{max_carts}_t{t}_rp{rp}_rl{rl}.gff", b = bin_list, fpr = fpr_list, min_len = min_lens, cmin = cmin_list, cmax = cmax_list, er = errors, bin_ent = bin_entropy_cutoffs, max_cap = cart_max_capacity, max_carts = max_queued_carts, t = search_threads, rp = repeat_periods, rl = repeat_lengths),
		truth_files = expand(blast_out + "/" + run_id + "_e{ev}_k{k}.bed", ev_evalues, k = blast_kmer_lengths)
	output:
		"valik.blast.accuracy"
	threads:
		workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = 10
	shell:
		"""
		echo -e "test-file\tmatches\ttruth-set-matches\ttrue-matches\tmissed\tmin-overlap\ttruth-file" > {output}
		for test in {input.test_files}
		do
			for truth in {input.truth_files}
			do
			
			match_count=`wc -l $test | awk '{{print $1}}'`
			echo -e "$test\t$match_count\t" >> {output}
			
			truncate -s -1 {output}
			../scripts/search_accuracy.sh $truth $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
			tail -n 1 tmp.log >> {output}
			rm tmp.log
	
			truncate -s -1 {output}
			echo -e "\t{params.min_overlap}\t$truth" >> {output}
			done
		done
		"""

rule valik_compare_stellar:
	input:
		ref_meta = dream_out + "/meta/b{b}_fpr{fpr}_l{l}_e{er}.bin",
		test_files = expand(dream_out + "/b{b}_fpr{fpr}_l{{min_len}}_cmin{cmin}_cmax{cmax}_e{{er}}_ent{bin_ent}_cap{max_cap}_carts{max_carts}_t{t}_rp{rp}_rl{rl}.gff", b = bin_list, fpr = fpr_list, cmin = cmin_list, cmax = cmax_list, bin_ent = bin_entropy_cutoffs, max_cap = cart_max_capacity, max_carts = max_queued_carts, t = search_threads, rp = repeat_periods, rl = repeat_lengths)
		truth_file = stellar_out + "/" + run_id + "_l{min_len}_e{er}.bed"
	output:
		temp(dream_out + "/valik.accuracy.l{min_len}.e{er}")
	threads:
		workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = 10,
	shell:
		"""
		echo -e "test-file\tmatches\ttruth-set-matches\ttrue-matches\tmissed\tmin-overlap\ttruth-file" > {output}
		
		for test in {input.test_files}
		do
			match_count=`wc -l $test | awk '{{print $1}}'`
			echo -e "$test\t$match_count\t" >> {output}
			
			truncate -s -1 {output}
			../scripts/search_accuracy.sh {input.truth} $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
			tail -n 1 tmp.log >> {output}
			rm tmp.log
	
			truncate -s -1 {output}
			echo -e "\t{params.min_overlap}\t$truth" >> {output}
			done
		done
		"""

rule valik_gather_stellar_accuracy:
	input:
		expand(dream_out + "/valik.accuracy.l{min_len}.e{er}, min_len = min_lens, er = errors)
	output:
		"valik.stellar.accuracy"
	threads: 1
	shell:
		"cat {input} > {output}"
