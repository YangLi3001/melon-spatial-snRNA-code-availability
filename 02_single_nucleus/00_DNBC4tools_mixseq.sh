#!/usr/bin/env bash
. activate DNBC4tools

i=P2
DNBC4tools run --mixseq \
  --cDNAfastq1 "data/$i/cDNA_R1.fq.gz" \
  --cDNAfastq2 "data/$i/cDNA_R2.fq.gz" \
  --oligofastq1 "data/$i/oligo_R1.fq.gz" \
  --oligofastq2 "data/$i/oligo_R2.fq.gz" \
  --starIndexDir ref/DHL92_V4 \
  --gtf ref/genes.gtf \
  --name "$i" \
  --thread 48 &>"$i.log"
