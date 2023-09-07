#!/bin/bash
set -x
#bin/mysysimage.jl
#JULIA='julia --sysimage sys_forplots2.so'

mkdir -p CABCHSV

JULIA='julia '
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
#exit

 $JULIA bin/PDFfit.jl -s 45 -p  Bernstein -d simulation-Bernstein-42 -n 250000 -c 4                 &> logs/1.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 250000 -c 4                 &> logs/2.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 250000 -c 4  --priorshift=1 &> logs/3.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 250000 -c 4  --priorshift=2 &> logs/4.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 250000 -c 4                 &> logs/5.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 250000 -c 4  --priorshift=1 &> logs/6.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 250000 -c 4  --priorshift=2 &> logs/7.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 250000 -c 4                 &> logs/8.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 250000 -c 4  --priorshift=1 &> logs/9.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 250000 -c 4  --priorshift=2 &> logs/10.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Valence   -d simulation-Valence-42   -n 250000 -c 4                 &> logs/11.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Valence   -d simulation-Dirichlet-42 -n 250000 -c 4                 &> logs/12.log&
 $JULIA bin/PDFfit.jl -s 45 -p  Bernstein -d simulation-Dirichlet-42 -n 250000 -c 4                 &> logs/13.log&
wait $(jobs -p)

        #$JULIA bin/fig8.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  &
        #$JULIA bin/fig567.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  &
        #$JULIA bin/fig2.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  & 

 
       $JULIA bin/fig34.jl      -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-43  &
       $JULIA bin/fig34.jl -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-43  &

       $JULIA bin/fig34.jl       -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  &
       $JULIA bin/fig34.jl  -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  &
       
       $JULIA bin/fig34.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-1-45-simulation-Dirichlet-42 --priorshift=1  &
       $JULIA bin/fig34.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-2-45-simulation-Dirichlet-42 --priorshift=2 &

       $JULIA bin/fig34.jl -s 47 -p  Bernstein -d simulation-Bernstein-42 -f fit-Bernstein-0-45-simulation-Bernstein-42 &
       $JULIA bin/fig34.jl -s 47 -p  Valence -d simulation-Valence-42 -f fit-Valence-0-45-simulation-Valence-42  &


       $JULIA bin/fig34.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-1-45-simulation-Dirichlet-43 --priorshift=1  &
       $JULIA bin/fig34.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-2-45-simulation-Dirichlet-43 --priorshift=2 &





 
 wait $(jobs -p)


mkdir -p figures2


cp figures/fig2-corner-fit-Dirichlet-0-45-simulation-Dirichlet-42_v3.pdf  figures2/fig2.pdf

cp figures/fig34d-fit-Dirichlet-0-45-simulation-Dirichlet-42_v2.pdf figures2/fig3a.pdf
cp figures/fig34d-fit-Dirichlet-0-45-simulation-Dirichlet-43_v2.pdf figures2/fig3b.pdf
cp figures/fig34-fit-Dirichlet-0-45-simulation-Dirichlet-42_v2.pdf  figures2/fig3c.pdf
cp figures/fig34-fit-Dirichlet-0-45-simulation-Dirichlet-43_v2.pdf  figures2/fig3d.pdf


cp figures/fig34-fit-Dirichlet-1-45-simulation-Dirichlet-42_v2.pdf   figures2/fig4a.pdf
cp figures/fig34-fit-Dirichlet-2-45-simulation-Dirichlet-42_v2.pdf   figures2/fig4b.pdf
cp figures/fig34-fit-Bernstein-0-45-simulation-Bernstein-42_v2.pdf  figures2/fig4c.pdf
cp figures/fig34-fit-Valence-0-45-simulation-Valence-42_v2.pdf      figures2/fig4d.pdf

cp figures/fig5-momentum-corr-fit-Dirichlet-0-45-simulation-Dirichlet-42_v2.pdf figures2/fig5.pdf
cp figures/'fig6-parton-xf(x)-fit-Dirichlet-0-45-simulation-Dirichlet-42_v2.pdf' figures2/fig6.pdf
cp figures/fig7-fit-Dirichlet-0-45-simulation-Dirichlet-42_v2.pdf figures2/fig7.pdf
cp figures/fig8-chisq-pvalue-fit-Dirichlet-0-45-simulation-Dirichlet-42_v2.pdf figures2/fig8.pdf


  echo "OK, DONE!!!"
 
 
 #~/.julia/packages/PartonDensity/VMoA4/src/parametrisations/bernstein.jl
