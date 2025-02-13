rule simulate_sequences:
	output:
		ref = "ref_rep{rep}.fasta",
		query = temp("random_rep{rep}.fasta")
	params:
		ref_seed = get_seed,
		query_seed = get_seed
	shell:      
		"{script_dir}/simulate_sequences.sh {wildcards.rep} {ref_len} {query_len} {params.ref_seed} {params.query_seed}"

rule simulate_matches:
	input:
		ref = "ref_rep{rep}.fasta",
		query = "random_rep{rep}.fasta"
	output:
		query = "query/rep{rep}_e{er}.fasta",
		matches = "local_matches/rep{rep}_e{er}.fasta",
		truth = "ground_truth/rep{rep}_e{er}.tsv"
	params:
		seed = get_seed
	shell:      
		"{script_dir}/simulate_local_matches.sh {wildcards.rep} {min_len} {max_len} {ref_len} {wildcards.er} {matches} {params.seed}"

