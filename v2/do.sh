#!/bin/bash
set -x
#bin/mysysimage.jl
JULIA='julia --sysimage sys_forplots.so'

 $JULIA bin/generatepseudodata.jl -s 43 -p Dirichlet &
 $JULIA bin/generatepseudodata.jl -s 43 -p Valence &

 $JULIA bin/generatepseudodata.jl -s 42 -p Dirichlet &
 $JULIA bin/generatepseudodata.jl -s 42 -p Valence &


 $JULIA bin/generatepseudodata.jl -s 42 -p Bernstein &


#wait $(jobs -p)


#$JULIA bin/PDFfit.jl -s 45 -p  Bernstein -d simulation-Bernstein-42 -n 1000 #Does not work so far

$JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 100000 &
$JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 100000 &
$JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 100000  --priorshift=1 &
$JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 100000  --priorshift=2 &
$JULIA bin/PDFfit.jl -s 45 -p  Valence -d simulation-Valence-42 -n 100000 &

#wait $(jobs -p)

 $JULIA bin/valenceplots.jl -s 47 -p  Valence -d simulation-Valence-42 -f fit-Valence-0-45-simulation-Valence-42 
 $JULIA bin/allplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42 
 $JULIA bin/allplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-1-45-simulation-Dirichlet-42 
 $JULIA bin/allplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-2-45-simulation-Dirichlet-42 
 $JULIA bin/gofplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42 
 
 
 #~/.julia/packages/PartonDensity/VMoA4/src/parametrisations/bernstein.jl