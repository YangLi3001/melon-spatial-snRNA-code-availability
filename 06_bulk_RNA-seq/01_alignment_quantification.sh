#!/usr/bin/env bash
set -euo pipefail

# Inputs: paired-end clean reads and the genome annotation.
reference=${REFERENCE:-ref/genome.fa}
annotation=${ANNOTATION:-ref/genes.gff}
fastq_dir=${FASTQ_DIR:-data/fastq}
out_dir=${OUT_DIR:-results/map}

samples=(KO1 KO2 WT1 WT2 WT3)
mkdir -p "$out_dir/quant"
hisat2-build "$reference" "$out_dir/genome" > "$out_dir/hisat2-build.log" 2>&1

for sample in "${samples[@]}"; do
    hisat2 -t -p 2 --dta -x "$out_dir/genome" \
        -1 "$fastq_dir/${sample}_R1.fq.gz" \
        -2 "$fastq_dir/${sample}_R2.fq.gz" \
        -S "$out_dir/$sample.sam" > "$out_dir/$sample.log" 2>&1
    samtools sort -@ 2 "$out_dir/$sample.sam" -o "$out_dir/$sample.bam"
    quant_dir="$out_dir/quant/$sample"
    mkdir -p "$quant_dir"
    stringtie -p 5 -B -G "$annotation" -e "$out_dir/$sample.bam" \
        -o "$quant_dir/$sample.gtf" \
        -A "$quant_dir/${sample}_abundance.txt" -l "$sample" \
        > "$quant_dir/run.log" 2>&1
done
