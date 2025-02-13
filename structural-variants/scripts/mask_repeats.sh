#!/bin/bash

set -e

execs=(dust st_dna5todna4)
for exec in "${execs[@]}"; do
    if ! which ${exec} &>/dev/null; then
        echo "${exec} is not available"
        echo ""
        echo "make sure \"${execs[@]}\" are reachable via the \${PATH} variable"
        echo ""

        echo "try "
	echo 'export PATH="/group/ag_abi/evelina/seqan3_tools/build/bin:/group/ag_abi/evelina/meme/bin:/group/ag_abi/evelina/meme/libexec/meme-5.5.6:$PATH"'

        exit 127
    fi
done

work_dir=$1

for sample in $work_dir/*/*/unmapped.fa;do
	sample_dir="$(dirname "$sample")"
	
	masked="$sample_dir/masked.fa"
	echo "Processing $sample"
	if [ ! -f $masked ]; then
		dust $sample > $masked
	fi

	seq_len=$(grep -v ">" $masked | wc -c | awk '{print $1}')
	repeat_len=$(grep -v ">" $masked | grep -o N | wc -l | awk '{print $1}')

	echo "Total sequence length $seq_len"
	echo "Repeat length $repeat_len"
	frac=$(bc <<< "scale=3; $repeat_len/$seq_len")
	echo "Fraction of repeats $frac"

	dna4="$sample_dir/dna4.fa"
	if [ ! -f $dna4 ]; then
		st_dna5todna4 $masked > $dna4	
	fi
done

