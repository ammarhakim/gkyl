#!/bin/bash -l

#.Declare a name for this job (preferably under 16 characters).
#SBATCH -J gkyl

#.Request the queue
#SBATCH -q regular

#.Number of nodes to request (Perlmutter has 64 cores and 4 GPUs per node)
#.Beware that there's also a shared queue for sharing a node (when you need
#.2 GPUs or less), which instead of -N uses -c or -n.
#SBATCH -N 1

#.Specify GPU needs:
#SBATCH --constraint gpu
#SBATCH --gpus 2

#.Request wall time
#SBATCH -t 00:30:00

#.Specify the account to charge. Some accounts need the _g suffix for GPUs.
#SBATCH --account=m6666

#.Mail is sent to you when the job starts and when it terminates or aborts.
#SBATCH --mail-user=jdoe@msn.com
#SBATCH --mail-type=END,FAIL,REQUEUE

module load python/3.9-anaconda-2021.11
module load openmpi/5.0.0rc12
module load cudatoolkit/12.0
module load nccl/2.18.3-cu12
module unload darshan

# For some reason we need to specify the full path to gkyl command in jobscript.
export gComDir="/global/homes/m/jdoe/gkylsoft/gkyl/bin"

#.Two ways to launch the job. Using mpirun may be faster.
echo 'srun --mpi=pmix -n 2 --gpus 2 '$gComDir'/gkyl -g input_file.lua'
srun --mpi=pmix -n 2 --gpus 2 $gComDir/gkyl -g input_file.lua
#echo 'mpirun -np 2 '$gComDir'/gkyl -g input_file.lua'
#mpirun -np 2 $gComDir/gkyl -g input_file.lua

exit 0
