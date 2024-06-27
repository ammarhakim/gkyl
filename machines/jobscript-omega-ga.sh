#!/bin/bash
#SBATCH -p medium
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 20
#SBATCH -n 20

#a Example cpu jobscript for GA Omega Cluster. 20 useable cores/node.

module load defaults
module unload python
module unload mdsplus
module load conda
conda activate pgkyl

export gComDir="$HOME/gkylsoft/gkyl/bin"
srun -n 20 $gComDir/gkyl *.lua


