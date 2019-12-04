// Gkyl ------------------------------------------------------------------------
//
// Header to C functions calling CUDA kernels which compute moments. 
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------
//
#ifndef DIST_FUNC_MOMENT_CALC_DEVICE_H 
#define DIST_FUNC_MOMENT_CALC_DEVICE_H 

extern "C" { 

void calcMom1x1vSer_M0_P1(int numBlocks, int numThreads, const double *w, const double *dxv, const double *f, double *out); 

} 
#endif 


