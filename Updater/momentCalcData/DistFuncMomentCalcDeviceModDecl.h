// Gkyl ------------------------------------------------------------------------
//
// Header to C functions calling CUDA kernels that compute moments of the
// distribution function.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <GkylRectCart.h>
#include <GkylRange.h>
#include <GkylCudaFuncs.h>
#include <GkylCartField.h>
#include <GkylBasisTypes.h>

#ifndef DIST_FUNC_MOMENT_CALC_DEVICE_MOD_DELC_H 
#define DIST_FUNC_MOMENT_CALC_DEVICE_MOD_DELC_H 

extern "C" { 

void cuda_MomentCalc1x1vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x1vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x1vSer_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc1x2vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x2vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x2vSer_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vSer_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc1x3vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x3vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x3vSer_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vSer_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc2x2vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x2vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x2vSer_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vSer_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc2x3vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x3vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x3vSer_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vSer_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc3x3vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc3x3vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc3x3vSer_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vSer_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc1x1vMax_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x1vMax_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x1vMax_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vMax_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc1x2vMax_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x2vMax_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x2vMax_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vMax_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc1x3vMax_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x3vMax_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x3vMax_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vMax_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc2x2vMax_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x2vMax_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x2vMax_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vMax_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc2x3vMax_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x3vMax_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x3vMax_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vMax_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc3x3vMax_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc3x3vMax_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc3x3vMax_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vMax_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc1x1vTensor_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x1vTensor_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x1vTensor_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vTensor_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc1x2vTensor_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x2vTensor_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x2vTensor_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x2vTensor_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc1x3vTensor_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x3vTensor_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc1x3vTensor_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x3vTensor_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc2x2vTensor_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x2vTensor_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x2vTensor_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x2vTensor_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc2x3vTensor_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x3vTensor_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc2x3vTensor_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc2x3vTensor_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 


void cuda_MomentCalc3x3vTensor_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc3x3vTensor_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

void cuda_MomentCalc3x3vTensor_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc3x3vTensor_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 



 
} 
#endif 
