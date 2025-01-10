#!/bin/bash

#creating index file

grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d "" -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat gencode.v45.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode
echo "indexing is done"



for r1 in ./Zahnreich_FASTQ_20240701/*R1_001.fastq.gz
do
r2=${r1/R1_001.fastq.gz/R2_001.fastq.gz}
sample_name=$(basename "$r1" "_R1_001.fastq.gz")
echo "$sample_name"
echo "working on $r1"
salmon quant -i salmon_index -l A -1 "$r1" -2 "$r2" --validateMappings --gcBias -o salmon_quant/"$sample_name"/
done
echo "quantification is done"

