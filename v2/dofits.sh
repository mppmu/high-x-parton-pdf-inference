#!/bin/bash
set -x
mkdir -p CABCHSV fitresults pseudodata fitlogs
module load singularity
export LC_ALL=C
export JULIA_PKGDIR=$(pwd)/J
export JULIA_DEPOT_PATH=$(pwd)/J
export SINGULARITY_TMPDIR=$(pwd)/tmp
export SINGULARITY_CACHEDIR=$(pwd)/tmp

JULIA='singularity exec  -B '$(pwd)':'$(pwd)' docker://ghcr.io/mppmu/high-x-parton-pdf-inference:latest julia'
$JULIA bin/install.jl

SCRIPTPATH=$($JULIA -e 'using PartonDensity; print(string(dirname(pathof(PartonDensity)),"/../utils/"))')
echo $SCRIPTPATH

JULIA="sbatch $(pwd)/doslurm.sh julia --compiled-modules=no "

#SCRIPTPATH=/raven/$SCRIPTPATH

 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42001 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42005 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42010 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42020 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42050 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42100 -n 200000 -c 4

 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichletnosyst -d simulation-Dirichlet-42001 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichletnosyst -d simulation-Dirichlet-42005 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichletnosyst -d simulation-Dirichlet-42010 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichletnosyst -d simulation-Dirichlet-42020 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichletnosyst -d simulation-Dirichlet-42050 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichletnosyst -d simulation-Dirichlet-42100 -n 200000 -c 4


 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Bernstein -d simulation-Bernstein-42 -n 200000 -c 4
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 200000 -c 4

 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 200000 -c 4  --priorshift=1
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-42 -n 200000 -c 4  --priorshift=2
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 200000 -c 4                
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 200000 -c 4  --priorshift=1 
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-43 -n 200000 -c 4  --priorshift=2 
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 200000 -c 4                 
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 200000 -c 4  --priorshift=1 
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Dirichlet -d simulation-Dirichlet-44 -n 200000 -c 4  --priorshift=2 
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Valence   -d simulation-Valence-42   -n 200000 -c 4                 
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Valence   -d simulation-Dirichlet-42 -n 200000 -c 4                
 $JULIA $SCRIPTPATH/PDFfit.jl -s 45 -p  Bernstein -d simulation-Dirichlet-42 -n 200000 -c 4                
