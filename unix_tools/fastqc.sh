#!/bin/bash
		
for sampleID in Zahnreich_FASTQ_20240701/*.fastq.gz; do
    sample_name=$(basename "$sampleID" ".fastq.gz")
    echo "working on $sample_name"
    mkdir -p ./fastqc/"$sample_name"
    fastqc "$sampleID" -o ./fastqc/"$sample_name"/
done
echo "fastqc finished running!"


