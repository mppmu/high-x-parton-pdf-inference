#!/bin/bash
set -x
mkdir -p CABCHSV fitresults pseudodata

export LC_ALL=C
export JULIA_PKGDIR=$(pwd)/J
export JULIA_DEPOT_PATH=$(pwd)/J

JULIA='singularity exec -B '$(pwd)':'$(pwd)' docker://ghcr.io/andriish/high-x-parton-pdf-inference:latest julia'
$JULIA bin/install.jl

SCRIPTPATH=$($JULIA -e 'using PartonDensity; print(string(dirname(pathof(PartonDensity)),"/../utils/"))')
echo $SCRIPTPATH

mkdir -p figures
       $JULIA bin/fig8.jl     -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  &
       $JULIA bin/fig7.jl     -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  &
#       $JULIA bin/fig567.jl   -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  &
       $JULIA bin/fig567_2.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  &
       $JULIA bin/fig2.jl     -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  &
#       $JULIA bin/fig10.jl    -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  & 
       $JULIA bin/fig10_2.jl  -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  & 



       $JULIA bin/fig34.jl -w x -s 47 -p  Dirichlet -d simulation-Dirichlet-43 -f fit-Dirichlet-0-45-simulation-Dirichlet-43  &
       $JULIA bin/fig34.jl -w u -s 47 -p  Dirichlet -d simulation-Dirichlet-43 -f fit-Dirichlet-0-45-simulation-Dirichlet-43  &
       $JULIA bin/fig34.jl -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-43 -f fit-Dirichlet-0-45-simulation-Dirichlet-43  &





       $JULIA bin/fig34.jl -w x -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  &
       $JULIA bin/fig34.jl -w u -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  &
       $JULIA bin/fig34.jl -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-42001 -f fit-Dirichlet-0-45-simulation-Dirichlet-42001  &



       $JULIA bin/fig34.jl -w x -s 47 -p  Dirichlet -d simulation-Dirichlet-42005 -f fit-Dirichlet-0-45-simulation-Dirichlet-42005  &
       $JULIA bin/fig34.jl -w u -s 47 -p  Dirichlet -d simulation-Dirichlet-42005 -f fit-Dirichlet-0-45-simulation-Dirichlet-42005  &
       $JULIA bin/fig34.jl -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-42005 -f fit-Dirichlet-0-45-simulation-Dirichlet-42005  &


       $JULIA bin/fig34.jl -w x -s 47 -p  Dirichlet -d simulation-Dirichlet-42010 -f fit-Dirichlet-0-45-simulation-Dirichlet-42010  &
       $JULIA bin/fig34.jl -w u -s 47 -p  Dirichlet -d simulation-Dirichlet-42010 -f fit-Dirichlet-0-45-simulation-Dirichlet-42010  &
       $JULIA bin/fig34.jl -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-42010 -f fit-Dirichlet-0-45-simulation-Dirichlet-42010  &




       $JULIA bin/fig34.jl -w x -s 47 -p  Dirichlet -d simulation-Dirichlet-42020 -f fit-Dirichlet-0-45-simulation-Dirichlet-42020  &
       $JULIA bin/fig34.jl -w u -s 47 -p  Dirichlet -d simulation-Dirichlet-42020 -f fit-Dirichlet-0-45-simulation-Dirichlet-42020  &
       $JULIA bin/fig34.jl -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-42020 -f fit-Dirichlet-0-45-simulation-Dirichlet-42020  &

