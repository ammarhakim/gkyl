#!/bin/bash
#SBATCH -A [Project]
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -t 0:30:00
#SBATCH -N 4
#a SBATCH -n 512
#SBATCH --ntasks-per-node 128


#a Example cpu jobscript for perlmutter. 128 useable cores/node (256 logical cores but it doesn't seem perlmutter has openmp).
#a The lua input file for this had decompcuts={64,2,4} (512 cores).
#a (change the queue to regular instead of debug of course for a longer job)

module load cray-mpich
module load PrgEnv-gnu/8.3.3
module unload adios
module unload zlib
module unload darshan

export gComDir="$HOME/gkylsoft/gkyl/bin"
srun -n 512 --cpu-bind=cores -c 1 $gComDir/gkyl -e 'finalTime=5.0e-8;numFrames=2' test_lowdim.lua


