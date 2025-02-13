There are two ways to check variant overlap:

Using the evaluate c++ executable matches each local match to possible variants
+ precise variant pos

cd /buffer/ag_abi/evelina/1000genomes/phase2/ftp.sra.ebi.ac.uk/vol1/run

evaluate --truth ../../../freeze3.sv.alt.meta.bed --test HG00731_l100_e0.033_simple.gff --ref-meta /group/ag_abi/evelina/DREAM-stellar-benchmark/genome-wise/human/dream/meta/b2048_fpr0.005_l100_e1.bin --overlap 0 --min-len 10

awk '{print $4 + 50}' HG00731_l100_e0.033_simple.fp.gff | rev | cut -c3- | rev | sort | uniq | wc -l

inversion_calling_accuracy.sh 
+ collapse alignments that differ by only few basepairs into a single variant
