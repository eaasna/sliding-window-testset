lastz_log = lastz_out + "/lastz.time"
f = open(lastz_log, "a")
f.write("time\tmem\terror-code\tcommand\tthreads\tseed\tgap-flag\ttransition-flag\tstep\tmatches\n")
f.close()

rule lastz_search:
	input:
		ref = ancient(dir_path(config["ref"]) + "dna4.fasta"),
		query = ancient(dir_path(config["query"]) + "dna4.fasta")
	output:
		lastz_out + "/" + run_id + "_s{s}_" + gap_flag + "_" + transition_flag + "_" + str(step_length) + ".tsv"
	threads: workflow.cores
	params:
		flags = "--" + gap_flag + " --" + transition_flag 
	shell:
		"""
		(timeout 24h /usr/bin/time -a -o {lastz_log} -f "%e\t%M\t%x\t%C\t1\t{wildcards.s}\t{gap_flag}\t{transition_flag}\t{step_length}" lastz_32 {input.ref}[multiple] {input.query} {params.flags} --hspthresh={hsp} \
				--seed={wildcards.s} --step={step_length} --progress=1 --format=blastn > {output})

		truncate -s -1 {lastz_log}
		wc -l {output} | awk '{{print "\t" $1}}' >> {lastz_log}
		"""

rule lastz_convert_to_blast:
	input:
		lastz_out + "/" + run_id + "_s{s}_" + gap_flag + "_" + transition_flag + "_" + str(step_length) + ".tsv"
	output:
		lastz_out + "/" + run_id + "_s{s}_" + gap_flag + "_" + transition_flag + "_" + str(step_length) + ".bed"
	threads: 1
	shell:
		"{shared_script_dir}/blast_like_to_bed.sh {input} {output}"
