#!/bin/bash
#SBATCH --job-name=vel_static
#SBATCH --output=log.interFlow
#SBATCH --ntasks=8
#SBATCH --partition=orbit

# run solver
of1812_hpcxv23

. $WM_PROJECT_DIR/bin/tools/RunFunctions

decomposePar
mpirun interFlow -parallel 

#------------------------------------------------------------------------------
