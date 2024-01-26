#!/bin/bash
set -x
mkdir -p CABCHSV fitresults pseudodata fitlogs
module load singularity
export LC_ALL=C
export JULIA_PKGDIR=$(pwd)/J
export JULIA_DEPOT_PATH=$(pwd)/J
export SINGULARITY_TMPDIR=$(pwd)/tmp
export SINGULARITY_CACHEDIR=$(pwd)/tmp
mkdir -p tmp cache
JULIA='singularity exec  -B '$(pwd)':'$(pwd)' docker://ghcr.io/andriish/high-x-parton-pdf-inference:latest julia'
$JULIA bin/install.jl
