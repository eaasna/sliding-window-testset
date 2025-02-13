rule valik_compare_blast:
	input:
		ref_meta = dream_out + "/meta/b" + str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_lens[0]) + "_e" + str(errors[0]) + ".bin",
		test_files = expand(dream_out + "/b{b}_fpr{fpr}_l{min_len}_cmin{cmin}_cmax{cmax}_e{er}_ent{bin_ent}_cap{max_cap}_carts{max_carts}_t{t}_rp{rp}_rl{rl}.gff", b = bin_list, fpr = fpr_list, min_len = min_lens, cmin = cmin_list, cmax = cmax_list, er = errors, bin_ent = bin_entropy_cutoffs, max_cap = cart_max_capacity, max_carts = max_queued_carts, t = search_threads, rp = repeat_periods, rl = repeat_lengths),
		truth_files = expand(blast_out + "/" + run_id + "_e{ev}_k{k}.bed", ev = evalues, k = blast_kmer_lengths)
	output:
		dream_out + "/valik.blast.accuracy"
	threads:
		workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = min_overlap
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
			{shared_script_dir}/search_accuracy.sh $truth $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
			tail -n 1 tmp.log >> {output}
			rm tmp.log
	
			truncate -s -1 {output}
			echo -e "\t{params.min_overlap}\t$truth" >> {output}
			done
		done
		"""

rule valik_kmer_compare_blast:
	input:
		ref_meta = dream_out + "/meta/b" + str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_lens[0]) + "_e" + str(errors[0]) + "_k" + str(valik_kmer_lengths[0]) + ".bin",
		test_files = expand(dream_out + "/b{b}_fpr{fpr}_l{min_len}_cmin{cmin}_cmax{cmax}_e{er}_k{k}_ent{bin_ent}_cap{max_cap}_carts{max_carts}_t{t}_rp{rp}_rl{rl}.gff", b = bin_list, fpr = fpr_list, min_len = min_lens, cmin = cmin_list, cmax = cmax_list, er = errors, k = valik_kmer_lengths, bin_ent = bin_entropy_cutoffs, max_cap = cart_max_capacity, max_carts = max_queued_carts, t = search_threads, rp = repeat_periods, rl = repeat_lengths),
		truth_files = expand(blast_out + "/" + run_id + "_e{ev}_k{k}.bed", ev = evalues, k = blast_kmer_lengths)
	output:
		dream_out + "/valik.blast.accuracy.k"
	threads:
		workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = min_overlap
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
			{shared_script_dir}/search_accuracy.sh $truth $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
			tail -n 1 tmp.log >> {output}
			rm tmp.log
	
			truncate -s -1 {output}
			echo -e "\t{params.min_overlap}\t$truth" >> {output}
			done
		done
		"""

rule valik_compare_stellar:
	input:
		ref_meta = dream_out + "/meta/b" + str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_lens[0]) + "_e" + str(errors[0]) + ".bin",
		test_files = expand(dream_out + "/b{b}_fpr{fpr}_l{{min_len}}_cmin{cmin}_cmax{cmax}_e{{er}}_ent{bin_ent}_cap{max_cap}_carts{max_carts}_t{t}_rp{rp}_rl{rl}.gff", b = bin_list, fpr = fpr_list, cmin = cmin_list, cmax = cmax_list, bin_ent = bin_entropy_cutoffs, max_cap = cart_max_capacity, max_carts = max_queued_carts, t = search_threads, rp = repeat_periods, rl = repeat_lengths),
		truth_file = stellar_out + "/" + run_id + "_l{min_len}_e{er}_rp" + str(repeat_periods[0]) + "_rl" + str(repeat_lengths[0]) + ".gff"
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
			{shared_script_dir}/search_accuracy.sh {input.truth_file} $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
			tail -n 1 tmp.log >> {output}
			rm tmp.log
	
			truncate -s -1 {output}
			echo -e "\t{params.min_overlap}\t{input.truth_file}" >> {output}
		done
		"""

rule valik_gather_stellar_accuracy:
	input:
		expand(dream_out + "/valik.accuracy.l{min_len}.e{er}", min_len = min_lens, er = errors)
	output:
		dream_out + "/valik.stellar.accuracy"
	threads: 1
	shell:
		"cat {input} > {output}"


rule valik_kmer_compare_stellar:
	input:
		ref_meta = dream_out + "/meta/b" + str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_lens[0]) + "_e" + str(errors[0]) + "_k" + str(valik_kmer_lengths[0]) + ".bin",
		test_files = expand(dream_out + "/b{b}_fpr{fpr}_l{{min_len}}_cmin{cmin}_cmax{cmax}_e{{er}}_k{k}_ent{bin_ent}_cap{max_cap}_carts{max_carts}_t{t}_rp{rp}_rl{rl}.gff", b = bin_list, fpr = fpr_list, cmin = cmin_list, cmax = cmax_list, k = valik_kmer_lengths, bin_ent = bin_entropy_cutoffs, max_cap = cart_max_capacity, max_carts = max_queued_carts, t = search_threads, rp = repeat_periods, rl = repeat_lengths),
		truth_file = stellar_out + "/" + run_id + "_l{min_len}_e{er}_rp" + str(repeat_periods[0]) + "_rl" + str(repeat_lengths[0]) + ".gff"
	output:
		temp(dream_out + "/valik.accuracy.l{min_len}.e{er}.k")
	threads:
		workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = min_overlap,
	shell:
		"""
		echo -e "test-file\tmatches\ttruth-set-matches\ttrue-matches\tmissed\tmin-overlap\ttruth-file" > {output}
		
		for test in {input.test_files}
		do
			match_count=`wc -l $test | awk '{{print $1}}'`
			echo -e "$test\t$match_count\t" >> {output}
			
			truncate -s -1 {output}
			{shared_script_dir}/search_accuracy.sh {input.truth_file} $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
			tail -n 1 tmp.log >> {output}
			rm tmp.log
	
			truncate -s -1 {output}
			echo -e "\t{params.min_overlap}\t{input.truth_file}" >> {output}
		done
		"""

rule valik_kmer_gather_stellar_accuracy:
	input:
		expand(dream_out + "/valik.accuracy.l{min_len}.e{er}.k", min_len = min_lens, er = errors)
	output:
		dream_out + "/valik.stellar.accuracy.k"
	threads: 1
	shell:
		"cat {input} > {output}"

rule valik_kmer_thresh_compare_stellar:
	input:
		ref_meta = dream_out + "/meta/b" + str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_lens[0]) + "_e" + str(errors[0]) + "_k" + str(valik_kmer_lengths[0]) + ".bin",
		test_files = expand(dream_out + "/b{b}_fpr{fpr}_l{{min_len}}_cmin{cmin}_cmax{cmax}_e{{er}}_k{k}_thresh{thresh}_ent{bin_ent}_cap{max_cap}_carts{max_carts}_t{t}_rp{rp}_rl{rl}.gff", b = bin_list, fpr = fpr_list, cmin = cmin_list, cmax = cmax_list, k = valik_kmer_lengths, thresh = valik_thresh, bin_ent = bin_entropy_cutoffs, max_cap = cart_max_capacity, max_carts = max_queued_carts, t = search_threads, rp = repeat_periods, rl = repeat_lengths),
		truth_file = stellar_out + "/" + run_id + "_l{min_len}_e{er}_rp" + str(repeat_periods[0]) + "_rl" + str(repeat_lengths[0]) + ".gff"
	output:
		temp(dream_out + "/valik.accuracy.l{min_len}.e{er}.k.t")
	threads:
		workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = min_overlap,
	shell:
		"""
		echo -e "test-file\tmatches\ttruth-set-matches\ttrue-matches\tmissed\tmin-overlap\ttruth-file" > {output}
		
		for test in {input.test_files}
		do
			match_count=`wc -l $test | awk '{{print $1}}'`
			echo -e "$test\t$match_count\t" >> {output}
			
			truncate -s -1 {output}
			{shared_script_dir}/search_accuracy.sh {input.truth_file} $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
			tail -n 1 tmp.log >> {output}
			rm tmp.log
	
			truncate -s -1 {output}
			echo -e "\t{params.min_overlap}\t{input.truth_file}" >> {output}
		done
		"""

rule valik_kmer_thresh_gather_stellar_accuracy:
	input:
		expand(dream_out + "/valik.accuracy.l{min_len}.e{er}.k.t", min_len = min_lens, er = errors)
	output:
		dream_out + "/valik.stellar.accuracy.k.t"
	threads: 1
	shell:
		"cat {input} > {output}"

rule valik_shape_compare_stellar:
	input:
		ref_meta = dream_out + "/meta/b" + str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_lens[0]) + "_e" + str(errors[0]) + "_s" + str(valik_shapes[0]) + ".bin",
		test_files = expand(dream_out + "/b{b}_fpr{fpr}_l{{min_len}}_cmin{cmin}_cmax{cmax}_e{{er}}_s{shape}_thresh{thresh}_ent{bin_ent}_cap{max_cap}_carts{max_carts}_t{t}.gff", b = bin_list, fpr = fpr_list, cmin = cmin_list, cmax = cmax_list, shape = valik_shapes, thresh = valik_thresh, bin_ent = bin_entropy_cutoffs, max_cap = cart_max_capacity, max_carts = max_queued_carts, t = search_threads),
		truth_file = stellar_out + "/" + run_id + "_l{min_len}_e{er}_rp" + str(repeat_periods[0]) + "_rl" + str(repeat_lengths[0]) + ".gff"
	output:
		temp(dream_out + "/valik.accuracy.l{min_len}.e{er}.s")
	threads:
		workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = min_overlap,
	shell:
		"""
		echo -e "test-file\tmatches\ttruth-set-matches\ttrue-matches\tmissed\tmin-overlap\ttruth-file" > {output}
		
		for test in {input.test_files}
		do
			match_count=`wc -l $test | awk '{{print $1}}'`
			echo -e "$test\t$match_count\t" >> {output}
			
			truncate -s -1 {output}
			{shared_script_dir}/search_accuracy.sh {input.truth_file} $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
			tail -n 1 tmp.log >> {output}
			rm tmp.log
	
			truncate -s -1 {output}
			echo -e "\t{params.min_overlap}\t{input.truth_file}" >> {output}
		done
		"""

rule valik_shape_gather_stellar_accuracy:
	input:
		expand(dream_out + "/valik.accuracy.l{min_len}.e{er}.s", min_len = min_lens, er = errors)
	output:
		dream_out + "/valik.stellar.accuracy.s"
	threads: 1
	shell:
		"cat {input} > {output}"

rule blast_compare_stellar:
	input:
		ref_meta = dream_out + "/meta/b"+ str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_lens[0]) + "_e" + str(errors[0]) + ".bin",
		test_files = expand(blast_out + "/" + run_id + "_e{ev}_k{k}.bed", ev = evalues, k = blast_kmer_lengths),
		truth_files = expand(stellar_out + "/" + run_id + "_l{l}_e{er}_rp{rp}_rl{rl}.gff", l = min_lens, er = errors, rp = repeat_periods, rl = repeat_lengths) 
	output:
		blast_out + "/blast.stellar.accuracy"
	threads: workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = min_overlap
		#min_overlap = round(min(min_lens) / 2)
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
			{shared_script_dir}/search_accuracy.sh $truth $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
			tail -n 1 tmp.log >> {output}
			rm tmp.log
	
			truncate -s -1 {output}
			echo -e "\t{params.min_overlap}\t$truth" >> {output}
			done
		done
		"""

rule blast_compare_valik:
	input:
		ref_meta = dream_out + "/meta/b"+ str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_lens[0]) + "_e" + str(errors[0]) + ".bin",
		test_files = expand(blast_out + "/" + run_id + "_e{ev}_k{k}.bed", ev = evalues, k = blast_kmer_lengths),
		truth_files = expand(dream_out + "/b{b}_fpr{fpr}_l{min_len}_cmin{cmin}_cmax{cmax}_e{er}_ent{bin_ent}_cap{max_cap}_carts{max_carts}_t{t}_rp{rp}_rl{rl}.gff", b = bin_list, fpr = fpr_list, min_len = min_lens, cmin = cmin_list, cmax = cmax_list, er = errors, bin_ent = bin_entropy_cutoffs, max_cap = cart_max_capacity, max_carts = max_queued_carts, t = search_threads, rp = repeat_periods, rl = repeat_lengths),
	output:
		blast_out + "/blast.valik.accuracy"
	threads: workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = min_overlap
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
			{shared_script_dir}/search_accuracy.sh $truth $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
			tail -n 1 tmp.log >> {output}
			rm tmp.log
	
			truncate -s -1 {output}
			echo -e "\t{params.min_overlap}\t$truth" >> {output}
			done
		done
		"""

rule lastz_compare_stellar:
	input:
		ref_meta = dream_out + "/meta/b" + str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_lens[0]) + "_e" + str(errors[0]) + ".bin",
		test_files = expand(lastz_out + "/" + run_id + "_s{s}_" + gap_flag + "_" + transition_flag + "_" + str(step_length) + ".bed", s = lastz_seeds),
		truth_files = expand(stellar_out + "/" + run_id + "_l{min_len}_e{er}_rp" + str(repeat_periods[0]) + "_rl" + str(repeat_lengths[0]) + ".gff", min_len = min_lens, er = errors)
	output:
		lastz_out + "/lastz.stellar.accuracy"
	threads:
		workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = min_overlap
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
				{shared_script_dir}/search_accuracy.sh $truth $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
				
				min_len=$(echo $truth | awk -F'_l' '{{print $2}}' | awk -F'_' '{{print $1}}')
				err=$(echo $truth | awk -F'_e' '{{print $2}}' | awk -F'_' '{{print $1}}')
				
				stellar_id="{lastz_out}/l${{min_len}}_e${{err}}"
				mkdir -p $stellar_id

				mv {lastz_out}/{run_id}_s*{gap_flag}*{transition_flag}*.fn.gff $stellar_id/
				mv {lastz_out}/{run_id}_s*{gap_flag}*{transition_flag}*.fp.bed $stellar_id/
				tail -n 1 tmp.log >> {output}
				rm tmp.log
	
				truncate -s -1 {output}
				echo -e "\t{params.min_overlap}\t$truth" >> {output}
			done
		done
		"""

rule last_compare_stellar:
	input:
		ref_meta = dream_out + "/meta/b" + str(bin_list[0]) + "_fpr" + str(fpr_list[0]) + "_l" + str(min_lens[0]) + "_e" + str(errors[0]) + ".bin",
		test_files = expand(last_out + "/" + run_id + "_w{last_w}_k{last_k}_m{last_m}.bed", last_w = last_index_every, last_k = last_query_every, last_m = last_initial_matches),
		truth_files = expand(stellar_out + "/" + run_id + "_l{min_len}_e{er}_rp" + str(repeat_periods[0]) + "_rl" + str(repeat_lengths[0]) + ".gff", min_len = min_lens, er = errors)
	output:
		last_out + "/last.stellar.accuracy"
	threads:
		workflow.cores
	params:
		min_len = min(min_lens),
		min_overlap = min_overlap
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
				{shared_script_dir}/search_accuracy.sh $truth $test {params.min_len} {params.min_overlap} {input.ref_meta} tmp.log
				
				min_len=$(echo $truth | awk -F'_l' '{{print $2}}' | awk -F'_' '{{print $1}}')
				err=$(echo $truth | awk -F'_e' '{{print $2}}' | awk -F'_' '{{print $1}}')
				
				stellar_id="{last_out}/l${{min_len}}_e${{err}}"
				mkdir -p $stellar_id

				mv {last_out}/{run_id}_w*_k*_m*.fn.gff $stellar_id/
				mv {last_out}/{run_id}_w*_k*_m*.fp.bed $stellar_id/
				tail -n 1 tmp.log >> {output}
				rm tmp.log
	
				truncate -s -1 {output}
				echo -e "\t{params.min_overlap}\t$truth" >> {output}
			done
		done
		"""
