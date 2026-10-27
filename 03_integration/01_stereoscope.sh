#!/usr/bin/env bash

for i in "MO" "C" "P" "U"
do
  stereoscope run \
    --sc_cnt "$i-sc.tsv" \
    --sc_labels "$i-sc_label.tsv" \
    -o "$i-out_marker" \
    --st_cnt "$i-spatial.tsv" \
    -gl marker.list \
    &>"$i-stereoscope.log"
done
