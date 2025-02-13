valik_split_log = "valik_split.time"
f = open(valik_split_log, "a")
f.write("time\tmem\texit-code\tcommand\tthreads\tbins\tfpr\terror-rate\tmin-len\tshape\n")
f.close()

rule valik_split_ref:
	input:
		ref = "ref_rep{rep}.fasta"
	output: 
		ref_meta = "meta/ref_rep{rep}_e{er}.bin"
	threads: 1
	shell:
		"""
		( /usr/bin/time -a -o {valik_split_log} -f "%e\t%M\t%x\t%C\t{threads}\t{bins}\t{fpr}\t{wildcards.er}\t{min_len}\t{shape}" valik split {input.ref} --verbose --out {output.ref_meta} --error-rate {wildcards.er} --pattern {min_len} --shape {shape} -n {bins} --fpr {fpr} &> {output}.split.err )
		"""

valik_build_log="valik_build.time"
f = open(valik_build_log, "a")
f.write("time\tmem\texit-code\tcommand\tthreads\tbins\tfpr\terror-rate\tmin-len\tshape\tibf_size\n")
f.close()
rule valik_build:
	input:
		ref = "ref_rep{rep}.fasta",
		ref_meta = "meta/ref_rep{rep}_e{er}.bin"
	output: 
		temp("/dev/shm/rep{rep}_e{er}.index")
	threads: workflow.cores
	benchmark:
		"benchmarks/valik_build_rep{rep}_e{er}.txt"
	shell:
		"""
		( /usr/bin/time -a -o {valik_build_log} -f "%e\t%M\t%x\t%C\t{threads}\t{bins}\t{fpr}\t{wildcards.er}\t{min_len}\t{shape}" valik build --threads {threads} --output {output} --ref-meta {input.ref_meta})

		truncate -s -1 {valik_build_log}
		ls -lh {output} | awk '{{ print "\t" $5 }}' >> {valik_build_log}
		"""

valik_search_log="valik_search.time"
f = open(valik_search_log, "a")
f.write("time\tmem\texit-code\tcommand\tthreads\tbins\tfpr\terror-rate\tmin-len\tshape\tthreshold\tmatches\n")
f.close()
rule valik_search:
	input:
		ibf = "/dev/shm/rep{rep}_e{er}.index",
		query = "query/rep{rep}_e{er}.fasta",
		ref_meta = "meta/ref_rep{rep}_e{er}.bin"
	output:
		"valik/rep{rep}_e{er}.gff"
	threads: workflow.cores
	params:
		t = get_threshold
	benchmark: 
		"benchmarks/valik_rep{rep}_e{er}.txt"
	shell:
		"""
		(/usr/bin/time -a -o {valik_search_log} -f "%e\t%M\t%x\t%C\t{threads}\t{bins}\t{fpr}\t{wildcards.er}\t{min_len}\t{shape}\t{params.t}" \
			valik search --keep-best-repeats --split-query --verbose --cache-thresholds \
				--numMatches {num_matches} --sortThresh {sort_thresh} --time \
				--index {input.ibf} --ref-meta {input.ref_meta} --query {input.query} \
				--error-rate {wildcards.er} --threads {threads} --output {output} \
				--cart-max-capacity {max_capacity} --max-queued-carts {max_carts} \
				--threshold {params.t} &> {output}.err)
		
		truncate -s -1 {valik_search_log}
		wc -l {output} | awk '{{ print $1 }}' >> {valik_search_log}
		"""

