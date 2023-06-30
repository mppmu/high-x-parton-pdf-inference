#!/bin/bash
set -x
#bin/mysysimage.jl
JULIA='julia --sysimage sys_forplots2.so'
#exit

  $JULIA bin/generatepseudodata.jl -s 42 -p Dirichlet &
  $JULIA bin/generatepseudodata.jl -s 42 -p Valence &
  $JULIA bin/generatepseudodata.jl -s 42 -p Bernstein &

  $JULIA bin/generatepseudodata.jl -s 43 -p Dirichlet &
  $JULIA bin/generatepseudodata.jl -s 43 -p Valence &
  $JULIA bin/generatepseudodata.jl -s 43 -p Bernstein &

  $JULIA bin/generatepseudodata.jl -s 44 -p Dirichlet &
  $JULIA bin/generatepseudodata.jl -s 44 -p Valence &
  $JULIA bin/generatepseudodata.jl -s 44 -p Bernstein &


wait $(jobs -p)


 $JULIA bin/PDFfit.jl -s 45 -p  Bernstein -d simulation-Bernstein-42 -n 100000
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 100000 &
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 100000 &
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 100000  --priorshift=1 &
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 100000  --priorshift=2 &
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 100000  --priorshift=1 &
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 100000  --priorshift=2 &
 $JULIA bin/PDFfit.jl -s 45 -p  Valence -d simulation-Valence-42 -n 100000 &

wait $(jobs -p)

 $JULIA bin/bernsteinplots.jl -s 47 -p  Bernstein -d simulation-Bernstein-42 -f fit-Bernstein-0-45-simulation-Bernstein-42 &
   $JULIA bin/valenceplots.jl -s 47 -p  Valence -d simulation-Valence-42 -f fit-Valence-0-45-simulation-Valence-42  &
     $JULIA bin/cornerplot.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  &
       $JULIA bin/gofplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  &

       $JULIA bin/allplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  &
       $JULIA bin/allplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-1-45-simulation-Dirichlet-42 --priorshift=1  &
       $JULIA bin/allplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-2-45-simulation-Dirichlet-42 --priorshift=2 &

       $JULIA bin/allplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-44  &
       $JULIA bin/allplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-1-45-simulation-Dirichlet-44 --priorshift=1  &
       $JULIA bin/allplots.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-2-45-simulation-Dirichlet-44 --priorshift=2 &
 
 wait $(jobs -p)
 
 echo "OK, DONE!!!"
 
 #~/.julia/packages/PartonDensity/VMoA4/src/parametrisations/bernstein.jl