#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J test_slurm
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#
#SBATCH --mail-type=none
#SBATCH --mail-user=userid@example.mpg.de
#
# Wall clock limit (max. is 24 hours):
#SBATCH --time=12:00:00
#SBATCH --mem=8000MB

module purge
module load apptainer

set -x
mkdir -p CABCHSV fitresults pseudodata

srun  singularity exec -B $(pwd):$(pwd) --env JULIA_DEPOT_PATH=$(pwd):/opt/julia docker://ghcr.io/andriish/high-x-parton-pdf-inference:latest $JULIA  "$@"