stellar_log = "stellar.time"
f = open(stellar_log, "a")
f.write("time\tmem\texit-code\tcommand\tthreads\terror-rate\tmatches\n")
f.close()

rule distribute_stellar:
	input:
		ref_meta = "meta/ref_rep{rep}_e{er}.bin",
		query = "query/rep{rep}_e{er}.fasta"
	output: 
		"dist_stellar/rep{rep}_e{er}.gff"
	threads: workflow.cores
	benchmark:
		"benchmarks/dist_stellar_rep{rep}_e{er}.txt"
	shell:
		"""	
		(/usr/bin/time -a -o {stellar_log} -f "%e\t%M\t%x\t%C\t{threads}\t{wildcards.er}" \
			valik search  --split-query --verbose \
				--numMatches {num_matches} --sortThresh {sort_thresh} --time \
				--ref-meta {input.ref_meta} --query {input.query} \
				--error-rate {wildcards.er} --threads {threads} --output {output} \
				--stellar-only &> {output}.err)
		
		truncate -s -1 {stellar_log}
		wc -l {output} | awk '{{ print $1 }}' >> {stellar_log}
		"""

rule stellar:
	input:
		ref = "ref_rep{rep}.fasta",
		query = "query/rep{rep}_e{er}.fasta"
	output: 
		"stellar/rep{rep}_e{er}.gff"
	threads: workflow.cores
	benchmark:
		"benchmarks/stellar_rep{rep}_e{er}.txt"
	shell:
		"""
		( timeout 12h /usr/bin/time -a -o {stellar_log} -f "%e\t%M\t%x\tstellar-search\t1\t{wildcards.er}" stellar --time -a dna --numMatches {num_matches}  --sortThresh {sort_thresh} {input.ref} {input.query} -e {wildcards.er} -l {min_len} -o {output} 2> {output}.err || touch {output} )
		
		truncate -s -1 {stellar_log}
		wc -l {output} | awk '{{ print $1 }}' >> {stellar_log}
		"""
		
