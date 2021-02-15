#!/bin/bash
#SBATCH --job-name=scriven3D
#SBATCH --output=log.multiRegionPhaseChangeFlow
#SBATCH --ntasks=32
#SBATCH --partition=orbit

# run solver
mpirun multiRegionPhaseChangeFlow -parallel

#------------------------------------------------------------------------------
