#!/bin/bash
set -x
mkdir -p CABCHSV fitresults pseudodata

export LC_ALL=C
export JULIA_PKGDIR=$(pwd)/J
export JULIA_DEPOT_PATH=$(pwd)/J
export SINGULARITY_TMPDIR=$(pwd)/tmp
export SINGULARITY_CACHEDIR=$(pwd)/tmp

JULIA='singularity exec -B '$(pwd)':'$(pwd)' docker://ghcr.io/andriish/high-x-parton-pdf-inference:latest julia'
#JULIA='singularity exec -B '$(pwd)':'$(pwd)' MPP-julia-fedora39-x86_64-v1.sif julia'
#$JULIA bin/install.jl

SCRIPTPATH=$($JULIA -e 'using PartonDensity; print(string(dirname(pathof(PartonDensity)),"/../utils/"))')
echo $SCRIPTPATH

mkdir -p figures

       $JULIA bin/fig567_2.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  &


       $JULIA bin/fig34.jl -w u -s 47 -p  Dirichlet -d simulation-Dirichlet-44 -f fit-Dirichlet-1-45-simulation-Dirichlet-44 --priorshift=1 &
       $JULIA bin/fig34.jl -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-44 -f fit-Dirichlet-1-45-simulation-Dirichlet-44 --priorshift=1  &



   ##    $JULIA bin/fig34.jl -w u -s 47 -p  Dirichlet -d simulation-Dirichlet-44 -f fit-Dirichlet-2-45-simulation-Dirichlet-44   --priorshift=2&
   ##    $JULIA bin/fig34.jl -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-44 -f fit-Dirichlet-2-45-simulation-Dirichlet-44   --priorshift=2&
