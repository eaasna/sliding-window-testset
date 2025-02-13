#!/bin/bash


work_dir="/buffer/ag_abi/evelina/1000genomes/phase2/ftp.sra.ebi.ac.uk/vol1/run/ERR386/ERR3861389"
cd $work_dir
bam_in="chr16.sam"
gff_out="chr16_from_sam.gff"

samtools view $bam_in | awk '{print $3 "\tpbm2\tmatches\t" $4 "\t" $4+ length($10) "\t100\t+\t.\t" $1 ";"}' > $gff_out
