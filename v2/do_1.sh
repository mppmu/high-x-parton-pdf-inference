#!/bin/bash
set -x
mkdir -p CABCHSV fitresults pseudodata

export LC_ALL=C
export JULIA_PKGDIR=$(pwd)/J
export JULIA_DEPOT_PATH=$(pwd)/J
JULIA='julia '
#  $JULIA bin/install.jl
JULIA='julia '
SCRIPTPATH=$($JULIA -e 'using PartonDensity; print(string(dirname(pathof(PartonDensity)),"/../utils/"))')
echo $SCRIPTPATH
#exit




  #$JULIA $SCRIPTPATH/generatepseudodata.jl -s 42 -p Dirichlet > logs/s1.log &
  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 1042 -p Dirichlet -f 1.0  > logs/s2.log&

  $JULIA $SCRIPTPATH/generatepseudodata.jl -s 10042 -p Dirichlet -f 100.0 > logs/s3.log &

  #$JULIA $SCRIPTPATH/generatepseudodata.jl -s 42 -p Valence  > logs/s4.log&
  #$JULIA $SCRIPTPATH/generatepseudodata.jl -s 42 -p Bernstein  > logs/s5.log&

  #$JULIA $SCRIPTPATH/generatepseudodata.jl -s 43 -p Dirichlet > logs/s6.log &
  #$JULIA $SCRIPTPATH/generatepseudodata.jl -s 43 -p Valence > logs/s7.log &
  #$JULIA $SCRIPTPATH/generatepseudodata.jl -s 43 -p Bernstein > logs/s8.log &

  #$JULIA $SCRIPTPATH/generatepseudodata.jl -s 44 -p Dirichlet > logs/s9.log &
  #$JULIA $SCRIPTPATH/generatepseudodata.jl -s 44 -p Valence > logs/s10.log &
  #$JULIA $SCRIPTPATH/generatepseudodata.jl -s 44 -p Bernstein > logs/s11.log &


wait $(jobs -p)


 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Bernstein -d simulation-Bernstein-42 -n 250000 -c 4                  &> logs/1.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 405 -p  Dirichlet -d simulation-Dirichlet-42 -n 2500 -c 4   --max_ncycles=1  --nsteps_per_cycle=10 --nsteps_final=10 --strict=false               &> logs/02.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 250000 -c 4                  &> logs/2.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 1045 -p  Dirichlet -d simulation-Dirichlet-1042 -n 250000 -c 4                  &> logs/102.log&


 $JULIA $SCRIPTPATH/PDFfit.jl -s 10045 -p  Dirichlet -d simulation-Dirichlet-10042 -n 250000 -c 4                  > logs/10002.log&

 $JULIA $SCRIPTPATH/PDFfit.jl -s 10045 -p  Dirichlet -d simulation-Dirichlet-1042 -n 250000 -c 4                  > logs/1002.log&

 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 250000 -c 4  --priorshift=1 &> logs/3.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 250000 -c 4  --priorshift=2 &> logs/4.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 250000 -c 4                 &> logs/5.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 250000 -c 4  --priorshift=1 &> logs/6.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 250000 -c 4  --priorshift=2 &> logs/7.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 250000 -c 4                 &> logs/8.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 250000 -c 4  --priorshift=1 &> logs/9.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 250000 -c 4  --priorshift=2 &> logs/10.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Valence   -d simulation-Valence-42   -n 250000 -c 4                 &> logs/11.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Valence   -d simulation-Dirichlet-42 -n 250000 -c 4                 &> logs/12.log&
 #$JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Bernstein -d simulation-Dirichlet-42 -n 250000 -c 4                 &> logs/13.log&
wait $(jobs -p)
exit

#$JULIA bin/fig34.jl  -w d  -s 47 -p  Dirichlet -d simulation-Dirichlet-10042 -f fit-Dirichlet-0-10045-simulation-Dirichlet-10042

#exit
mkdir -p figures


       $JULIA bin/fig8.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42   

#exit

       $JULIA bin/fig567.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42 
#       
     
#exit
        $JULIA bin/fig2.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  
       $JULIA bin/fig10.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42 

#exit

#exit
#wait $(jobs -p)
#exit




       $JULIA bin/fig34.jl      -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-43  > logs/p34f1.log &
       $JULIA bin/fig34.jl -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-43 > logs/p34f2.log &

       $JULIA bin/fig34.jl       -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  > logs/p34f3.log&
       $JULIA bin/fig34.jl  -w d -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-0-45-simulation-Dirichlet-42  > logs/p34f4.log&
       
       $JULIA bin/fig34.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-1-45-simulation-Dirichlet-42 --priorshift=1  > logs/p34f5.log&
       $JULIA bin/fig34.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-2-45-simulation-Dirichlet-42 --priorshift=2 > logs/p34f6.log&

       $JULIA bin/fig34.jl -s 47 -p  Bernstein -d simulation-Bernstein-42 -f fit-Bernstein-0-45-simulation-Bernstein-42 > logs/p34f7.log&
       $JULIA bin/fig34.jl -s 47 -p  Valence -d simulation-Valence-42 -f fit-Valence-0-45-simulation-Valence-42  > logs/p34f8.log&


       $JULIA bin/fig34.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-1-45-simulation-Dirichlet-43 --priorshift=1  > logs/p34f9.log&
       $JULIA bin/fig34.jl -s 47 -p  Dirichlet -d simulation-Dirichlet-42 -f fit-Dirichlet-2-45-simulation-Dirichlet-43 --priorshift=2 > logs/p34f10.log&


       $JULIA bin/fig34.jl        -s 47 -p  Dirichlet -d simulation-Dirichlet-1042 -f fit-Dirichlet-0-1045-simulation-Dirichlet-1042  > logs/p34f11.log& 
       $JULIA bin/fig34.jl  -w d  -s 47 -p  Dirichlet -d simulation-Dirichlet-1042 -f fit-Dirichlet-0-1045-simulation-Dirichlet-1042  > logs/p34f12.log& 


       $JULIA bin/fig34.jl        -s 47 -p  Dirichlet -d simulation-Dirichlet-10042 -f fit-Dirichlet-0-10045-simulation-Dirichlet-10042  > logs/p34f13.log& 
       $JULIA bin/fig34.jl  -w d  -s 47 -p  Dirichlet -d simulation-Dirichlet-10042 -f fit-Dirichlet-0-10045-simulation-Dirichlet-10042  > logs/p34f14.log& 


wait $(jobs -p)
#exit

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
