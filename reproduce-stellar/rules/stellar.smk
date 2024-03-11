rule stellar:
        input:
                ref = "ref_rep{rep}.fasta",
                query = "query/rep{rep}_e{er}.fasta",
		mutex = "blast_table1.tsv"
        output: 
                "stellar/rep{rep}_e{er}.gff"
        benchmark:
                "benchmarks/stellar_rep{rep}_e{er}.txt"
        shell:
                """( timeout 1h /usr/bin/time -a -o stellar.time -f "%e\t%M\t%x\tstellar-search\t{wildcards.er}" ../../../stellar3/build/bin/stellar -a dna --numMatches {num_matches}  --sortThresh {sort_thresh} {input.ref} {input.query} -e {wildcards.er} -l {min_len} -o {output} || touch {output} )"""
