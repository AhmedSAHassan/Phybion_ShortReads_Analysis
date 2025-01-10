#!/bin/bash
		
for sampleID in Zahnreich_FASTQ_20240701/*.fastq.gz; do
    sample_name=$(basename "$sampleID" ".fastq.gz")
    echo "working on $sample_name"
    mkdir -p ./fastqc/"$sample_name"
    fastqc "$sampleID" -o ./fastqc/"$sample_name"/
done
echo "fastqc finished running!"


for sampleID in trimmed/t_*/*paired.fq.gz; do
    sample_name=$(basename "$sampleID" "paired.fq.gz")
    echo "working on $sample_name"
    mkdir -p ./fastqc_trimmed/"$sample_name"
    fastqc "$sampleID" -o ./fastqc_trimmed/"$sample_name"/
done
echo "fastqc finished running!"



for sampleID in fastp/t_*/*.fq.gz; do
    sample_name=$(basename "$sampleID" ".fq.gz")
    echo "working on $sample_name"
    mkdir -p ./fastqc_fastp/"$sample_name"
    fastqc "$sampleID" -o ./fastqc_fastp/"$sample_name"/
done
echo "fastqc finished running!"
