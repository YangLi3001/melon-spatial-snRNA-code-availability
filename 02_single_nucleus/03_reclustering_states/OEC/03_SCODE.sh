#!/usr/bin/env bash
# Run from the SCODE R implementation directory with the study input matrices.
# SCODE: https://github.com/hmatsu1226/SCODE
# Arguments: genes, latent dimensions, cells, iterations, independent runs.

ruby run_R.rb expression_matrix.txt pseudotime.txt SCODE_out \
  50 6 3071 100 100 >run.log 2>run.err
