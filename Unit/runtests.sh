GKYL=$HOME/gkylsoft/gkyl/bin/xgkyl
MPIEXEC=$HOME/gkylsoft/openmpi-2.0.1/bin/mpiexec

# Adios
cmd="$GKYL test_Adios.lua"
echo $cmd
$cmd

# AdiosCartFieldIo
cmd="$GKYL test_AdiosCartFieldIo.lua"
echo $cmd
$cmd

# Advection
cmd="$GKYL test_Advection.lua"
echo $cmd
$cmd

# Alloc
cmd="$GKYL test_Alloc.lua"
echo $cmd
$cmd

# AllocShared
cmd="$MPIEXEC -n 4 $GKYL test_AllocShared.lua"
echo $cmd
$cmd

# Basis
cmd="$GKYL test_Basis.lua"
echo $cmd
$cmd

# BoundaryCondition
cmd="$GKYL test_BoundaryCondition.lua"
echo $cmd
$cmd

# Burgers
cmd="$GKYL test_Burgers.lua"
echo $cmd
$cmd

# CartDecomp
cmd="$GKYL test_CartDecomp.lua"
echo $cmd
$cmd

# CartDecompNeigh
cmd="$GKYL test_CartDecompNeigh.lua"
echo $cmd
$cmd

# CartDecompNeighShared
cmd="$MPIEXEC -n 4 $GKYL test_CartDecompShared.lua"
echo $cmd
$cmd

# CartField
cmd="$GKYL test_CartField.lua"
echo $cmd
$cmd

# CartGrid
cmd="$GKYL test_CartGrid.lua"
echo $cmd
$cmd

# CartGridPar
cmd="$MPIEXEC -n 4 $GKYL test_CartGridPar.lua"
echo $cmd
$cmd

# DistFuncMomentCalc
cmd="$GKYL test_DistFuncMomentCalc.lua"
echo $cmd
$cmd

# DynVector
cmd="$GKYL test_DynVector.lua"
echo $cmd
$cmd

# Euler
cmd="$GKYL test_Euler.lua"
echo $cmd
$cmd

# FiveMomentSrc
cmd="$GKYL test_FiveMomentSrc.lua"
echo $cmd
$cmd

# GaussQuad
cmd="$GKYL test_GaussQuad.lua"
echo $cmd
$cmd

# Lin
cmd="$GKYL test_Lin.lua"
echo $cmd
$cmd

# LinearDecomp
cmd="$GKYL test_LinearDecomp.lua"
echo $cmd
$cmd

# Mpi (n=1)
cmd="$GKYL test_Mpi.lua"
echo $cmd
$cmd

# Mpi (n=4)
cmd="$MPIEXEC -n 4 $GKYL test_Mpi.lua"
echo $cmd
$cmd

# Mpi (n=2)
cmd="$MPIEXEC -n 2 $GKYL test_Mpi.lua"
echo $cmd
$cmd

# ParCartField (n=2)
cmd="$MPIEXEC -n 2 $GKYL test_ParCartField.lua"
echo $cmd
$cmd

# ParCartField (n=3)
cmd="$MPIEXEC -n 3 $GKYL test_ParCartField.lua"
echo $cmd
$cmd

# ParCartField (n=4)
cmd="$MPIEXEC -n 4 $GKYL test_ParCartField.lua"
echo $cmd
$cmd

# ParCartGrid (n=2)
cmd="$MPIEXEC -n 2 $GKYL test_ParCartGrid.lua"
echo $cmd
$cmd

# ParCartGrid (n=4)
cmd="$MPIEXEC -n 4 $GKYL test_ParCartGrid.lua"
echo $cmd
$cmd

# PerfMaxwell
cmd="$GKYL test_PerfMaxwell.lua"
echo $cmd
$cmd

# ProjectOnBasis
cmd="$GKYL test_ProjectOnBasis.lua"
echo $cmd
$cmd

# Range
cmd="$GKYL test_Range.lua"
echo $cmd
$cmd

# SparseTriples
cmd="$GKYL test_SparseTriples.lua"
echo $cmd
$cmd

# Template
cmd="$GKYL test_Template.lua"
echo $cmd
$cmd

# TenMoment
cmd="$GKYL test_TenMoment.lua"
echo $cmd
$cmd

# cfuncs
cmd="$GKYL test_cfuncs.lua"
echo $cmd
$cmd
