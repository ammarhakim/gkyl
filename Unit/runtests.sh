GKYL=$HOME/gkylsoft/gkyl/bin/xgkyl
MPIEXEC=$HOME/gkylsoft/openmpi-3.0.0/bin/mpiexec

ctext=`tput setaf 5`
red=`tput setaf 1`
reset=`tput sgr0`

# Adios
cmd="$GKYL test_Adios.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Adios (parallel)
cmd="$MPIEXEC -n 4 $GKYL test_Adios.lua"
echo "${ctext}${cmd}${reset}"
echo "${red}$cmd NOT RUN!!!!${reset}"

# AdiosCartFieldIo
cmd="$GKYL test_AdiosCartFieldIo.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Advection
cmd="$GKYL test_Advection.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Alloc
cmd="$GKYL test_Alloc.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# AllocShared
cmd="$MPIEXEC -n 4 $GKYL test_AllocShared.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Basis
cmd="$GKYL test_Basis.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# BoundaryCondition
cmd="$GKYL test_BoundaryCondition.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Burgers
cmd="$GKYL test_Burgers.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# CartDecomp
cmd="$GKYL test_CartDecomp.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# CartDecompNeigh
cmd="$GKYL test_CartDecompNeigh.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# CartDecompNeighShared
cmd="$MPIEXEC -n 4 $GKYL test_CartDecompShared.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# CartField
cmd="$GKYL test_CartField.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# CartGrid
cmd="$GKYL test_CartGrid.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# CartGridPar
cmd="$MPIEXEC -n 4 $GKYL test_CartGridPar.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# DistFuncMomentCalc
cmd="$GKYL test_DistFuncMomentCalc.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# DynVector
cmd="$GKYL test_DynVector.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Euler
cmd="$GKYL test_Euler.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# FiveMomentSrc
cmd="$GKYL test_FiveMomentSrc.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# GaussQuad
cmd="$GKYL test_GaussQuad.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Lin
cmd="$GKYL test_Lin.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# LinearDecomp
cmd="$GKYL test_LinearDecomp.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Mpi (n=1)
cmd="$GKYL test_Mpi.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Mpi (n=4)
cmd="$MPIEXEC -n 4 $GKYL test_Mpi.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Mpi (n=2)
cmd="$MPIEXEC -n 2 $GKYL test_Mpi.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# ParCartField (n=2)
cmd="$MPIEXEC -n 2 $GKYL test_ParCartField.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# ParCartField (n=3)
cmd="$MPIEXEC -n 3 $GKYL test_ParCartField.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# ParCartField (n=4)
cmd="$MPIEXEC -n 4 $GKYL test_ParCartField.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# ParCartGrid (n=2)
cmd="$MPIEXEC -n 2 $GKYL test_ParCartGrid.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# ParCartGrid (n=4)
cmd="$MPIEXEC -n 4 $GKYL test_ParCartGrid.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# PerfMaxwell
cmd="$GKYL test_PerfMaxwell.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# ProjectOnBasis
cmd="$GKYL test_ProjectOnBasis.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Range
cmd="$GKYL test_Range.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# SparseTriples
cmd="$GKYL test_SparseTriples.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# Template
cmd="$GKYL test_Template.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# TenMoment
cmd="$GKYL test_TenMoment.lua"
echo "${ctext}${cmd}${reset}"
$cmd

# cfuncs
cmd="$GKYL test_cfuncs.lua"
echo "${ctext}${cmd}${reset}"
$cmd
