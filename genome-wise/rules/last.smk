last_log = last_out + "/last.time"
f = open(last_log, "a")
f.write("time\tmem\terror-code\tcommand\tthreads\tindex-every\tquery-every\tinitial-matches\tmatches\n")
f.close()

def db_name(wildcards):
	return last_out + "/" + run_id + "_w" + wildcards.last_w

# https://gitlab.com/mcfrith/last/-/blob/main/doc/last-tuning.rst
# -c: prevents lastal using huge amounts of memory and time on many ways of aligning centromeric repeats.
# -P: parallel threads
# -w: index every nth position in the reference 
#	this increases speed the most while reducing sensitivity the least
#	MegaBLAST indexes every 5th position
rule last_index:
	input:
		dir_path(config["ref"]) + "dna4.fasta"
	output:
		last_out + "/" + run_id + "_w{last_w}.suf"
	threads: workflow.cores
	params:
		db = db_name
	shell:
		"""
		( /usr/bin/time -a -o {last_log} -f "%e\t%M\t%x\t%C\t{threads}\t{wildcards.last_w}" \
			lastdb {params.db} -c -P{threads} -U {max_tandem_repeats} {input} 2> last_index.err )
		"""

# https://gitlab.com/mcfrith/last/-/blob/main/doc/last-seeds.rst
# -k: match every nth position in the query
#	this decreases memory and disk use the most while reducing sensitivity the least.
# -l: minimum length of initial exact seeds
# -u: seeding scheme?? not recongnized
#	YASS: DNA default suited for long and weak similarities
#	RY16: reduce run time and memory use by only seeking seeds at 1/16 of positions
rule last_search:
	input:
		ref = dir_path(config["ref"]) + "dna4.fasta",
		db = last_out + "/" + run_id + "_w{last_w}.suf",
		query = dir_path(config["query"]) + "dna4.fasta"
	output:
		last_out + "/" + run_id + "_w{last_w}_k{last_k}_m{last_m}.tsv"
	threads: workflow.cores
	params:
		db = db_name
	shell:
		"""
		(timeout 24h /usr/bin/time -a -o {last_log} -f "%e\t%M\t%x\t%C\t{threads}\t{wildcards.last_w}\t{wildcards.last_k}\t{wildcards.last_m}" \
			lastal -P{threads} {params.db} {input.query} \
				-k {wildcards.last_k} \
				-f BlastTab -m {wildcards.last_m} > {output})

		truncate -s -1 {last_log}
		grep -v [#] {output} | wc -l | awk '{{print "\t" $1}}' >> {last_log}
		"""

# BLAST reverse strand matches dstart > dend
# LAST reverse strand matches qstart > qend
# Stellar reverse strand matches (dstart < dend) & (qstart < qend)
rule last_convert_to_blast:
	input:
		last_out + "/" + run_id + "_w{last_w}_k{last_k}_m{last_m}.tsv"
	output:
		last_out + "/" + run_id + "_w{last_w}_k{last_k}_m{last_m}.bed"
	threads: 1
	shell:
		"{shared_script_dir}/blast_like_to_bed.sh {input} {output}"

