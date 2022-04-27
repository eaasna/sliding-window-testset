rule distribute_search:
	input:
		queries = "rep{rep}/queries/e{er}.fastq",
		search_out = "rep{rep}/search/e{er}.out"
	output:
		expand("rep{{rep}}/queries/seg{bin}_e{{er}}.fasta", bin = bin_list)
	benchmark:
		"benchmarks/rep{rep}/dream_stellar/distribute_search_e{er}.txt"
	script:
		"../scripts/distribute_search.py"

rule dream_stellar_search:
	input:
		ref_seg = "rep{rep}/split/seg{bin}.fasta",
		query = "rep{rep}/queries/seg{bin}_e{er}.fasta"
	output:
		"rep{rep}/dream_stellar/seg{bin}_e{er}.gff"
	params:
		e = get_error_rate,
		dummy_list = "benchmarks/rep{rep}/dream_stellar/dummy.txt"
	benchmark:
		"benchmarks/rep{rep}/dream_stellar/seg{bin}_e{er}.txt"
	shell:
		"""
		if [ -s {input.query} ]; then
		        # Search queries for current bin
			stellar {input.ref_seg} {input.query} --forward -e {params.e} -l {pattern} --numMatches {num} --sortThresh {thresh} -a dna -o {output}
		else
			touch {output} # create dummy output
		fi
		"""
	
