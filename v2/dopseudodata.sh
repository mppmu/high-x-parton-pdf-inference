#!/bin/bash
set -x
module load singularity
mkdir -p CABCHSV fitresults pseudodata

export LC_ALL=C
export JULIA_PKGDIR=$(pwd)/J
export JULIA_DEPOT_PATH=$(pwd)/J
export SINGULARITY_TMPDIR=$(pwd)/tmp
export SINGULARITY_CACHEDIR=$(pwd)/tmp

JULIA='singularity exec -B '$(pwd)':'$(pwd)' docker://ghcr.io/andriish/high-x-parton-pdf-inference:latest julia'
$JULIA bin/install.jl

SCRIPTPATH=$($JULIA -e 'using PartonDensity; print(string(dirname(pathof(PartonDensity)),"/../utils/"))')
echo $SCRIPTPATH


  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 42001 -p Dirichlet -f 1.0 &> logs/r1.log  &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 42005 -p Dirichlet -f 5.0 &> logs/r5.log  & 
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 42010 -p Dirichlet -f 10.0 &> logs/r10.log &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 42020 -p Dirichlet -f 20.0 &> logs/r20.log  &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 42050 -p Dirichlet -f 50.0 &> logs/r50.log  &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 42100 -p Dirichlet -f 100.0 &> logs/r100.log  &



  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 42 -p Dirichlet  &> logs/s3.log &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 42 -p Valence  &> logs/s4.log &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 42 -p Bernstein  &> logs/s5.log &

  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 43 -p Dirichlet &> logs/s6.log  &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 43 -p Valence &> logs/s7.log  &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 43 -p Bernstein &> logs/s8.log  &

  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 44 -p Dirichlet &> logs/s9.log  &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 44 -p Valence &> logs/s10.log  &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 44 -p Bernstein &> logs/s11.log  &

