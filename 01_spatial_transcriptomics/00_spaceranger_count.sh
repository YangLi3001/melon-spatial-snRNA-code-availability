#!/usr/bin/env bash
genome=ref/DHL92_V4
thread=20

run_count() {
  spaceranger count --id="$id" \
    --transcriptome="$genome" \
    --fastqs="$fastq_path" \
    --image="$image_file" \
    --slide=V11F03-104 \
    --slidefile=ref/slide.txt \
    --area="$area" \
    --localcores="$thread" \
    --localmem=150
}

id=ST1
fastq_path=data/ST1
image_file=data/ST1.jpg
area=A1
run_count

id=ST2
fastq_path=data/ST2
image_file=data/ST2.jpg
area=B1
run_count
