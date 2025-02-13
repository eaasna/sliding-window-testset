#!/bin/bash

set -ex

execs=(st_dna5todna4)
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

# this does not seem necessary
work_dir="/buffer/ag_abi/evelina/1000genomes/phase2/ftp.sra.ebi.ac.uk/vol1/run"
for sample in $work_dir/*/*/unmapped.fa;do
	sample_dir="$(dirname "$sample")"
	
	dna4="$sample_dir/unmasked_dna4.fa"
	if [ ! -f $dna4 ]; then
		st_dna5todna4 $sample > $dna4	
	fi
done

