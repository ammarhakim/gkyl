// Gkyl ------------------------------------------------------------------------
//
// Header to C functions called by CUDA kernels which compute moments.
// Ideally this wouldn't exist and we would just modify the existing header file.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------
//
#ifndef DIST_FUNC_MOMENT_CALC_MOD_DECL_DEVICE_H
#define DIST_FUNC_MOMENT_CALC_MOD_DECL_DEVICE_H 

extern "C" { 

__host__ __device__ void MomentCalc1x1vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 

} 
#endif 


