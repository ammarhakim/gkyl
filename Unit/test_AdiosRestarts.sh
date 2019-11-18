#!/bin/bash
#.
#.Test simulation restart Adios infrastructure for parallel jobs.
#.
#.In this script the user must specify:
#.  endTime: final simulation timetime.
#.  frames:  number of frames to output
#.  decomp:  MPI decomposition.
#.for each of 3 runs of the input file chosen. The name of the input file
#.must be the same for all 3 runs, and must match the string passed to
#.the file comparison (test_AdiosRestartsCompare) through `simulationName`.
#.Also pass the frame number through 'compareFrame' to test_AdiosRestartsCompare,
#.it must be the same number as the frame of the files copied after run #1.
#.
#.Lines commented are to switch between serial tests on a laptop,
#.and tests of parallel restarts on a cluster.
#.

#module load intel/2018-01

#.For some reason we need to specify the full path to gkyl command.
export gComDir="$HOME/gkylsoft/gkyl/bin"
#export mpiComDir="$HOME/gkylsoft/openmpi-3.1.2/bin/"

#.Run #1: do a first simulation to produce reference data.
#$mpiComDir/mpirun -n 2 $gComDir/gkyl -e "endTime=10.0; frames=2; decomp={1}" two-stream-p2.lua
$gComDir/gkyl -e "endTime=10.0; frames=2; decomp={1}" two-stream-p2.lua

#.Copy the distribution function and its moments into reference files.
cp two-stream-p2_elc_2.bp two-stream-p2_elc_2ref.bp
cp two-stream-p2_elc_M0_2.bp two-stream-p2_elc_M0_2ref.bp
cp two-stream-p2_elc_M1i_2.bp two-stream-p2_elc_M1i_2ref.bp
cp two-stream-p2_elc_M2_2.bp two-stream-p2_elc_M2_2ref.bp
cp two-stream-p2_field_2.bp two-stream-p2_field_2ref.bp

#.Run #2: first half of the simulation with same number of processes as reference run.
#$mpiComDir/mpirun -n 2 $gComDir/gkyl -e "endTime=5.0; frames=1; decomp={2}" two-stream-p2.lua
$gComDir/gkyl -e "endTime=5.0; frames=1; decomp={1}" two-stream-p2.lua

#.Run #3: restart simulation with a different number of processes.
#$mpiComDir/mpirun -n 4 $gComDir/gkyl -e "endTime=10.0; frames=2; decomp={4}" two-stream-p2.lua restart
$gComDir/gkyl -e "endTime=10.0; frames=2; decomp={1}" two-stream-p2.lua restart

#.Run a Lua script that compares files from runs #1 and #3.
$gComDir/gkyl -e "simulationName=\"two-stream-p2\"; compareFrame=2" test_AdiosRestartsCompare.lua
