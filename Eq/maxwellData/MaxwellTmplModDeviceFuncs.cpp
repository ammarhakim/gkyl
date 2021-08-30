// Gkyl ------------------------------------------------------------------------
//
// Code to copy various kernels from device to host
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <MaxwellTmplModDecl.h>

using namespace Gkyl;



__device__ Maxwell_volumeTerm_t p_Maxwell_volumeTerm_1xSerP1 =
  &Gkyl::MaxwellModDecl<1,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Maxwell_surfTerm_t p_Maxwell_surfTerm_1xSerP1 =
  &Gkyl::MaxwellModDecl<1,1,G_SERENDIPITY_C>::surfTerm;


__device__ Maxwell_volumeTerm_t p_Maxwell_volumeTerm_1xSerP2 =
  &Gkyl::MaxwellModDecl<1,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Maxwell_surfTerm_t p_Maxwell_surfTerm_1xSerP2 =
  &Gkyl::MaxwellModDecl<1,2,G_SERENDIPITY_C>::surfTerm;


__device__ Maxwell_volumeTerm_t p_Maxwell_volumeTerm_2xSerP1 =
  &Gkyl::MaxwellModDecl<2,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Maxwell_surfTerm_t p_Maxwell_surfTerm_2xSerP1 =
  &Gkyl::MaxwellModDecl<2,1,G_SERENDIPITY_C>::surfTerm;


__device__ Maxwell_volumeTerm_t p_Maxwell_volumeTerm_2xSerP2 =
  &Gkyl::MaxwellModDecl<2,2,G_SERENDIPITY_C>::volumeTerm;

__device__ Maxwell_surfTerm_t p_Maxwell_surfTerm_2xSerP2 =
  &Gkyl::MaxwellModDecl<2,2,G_SERENDIPITY_C>::surfTerm;




__device__ Maxwell_volumeTerm_t p_Maxwell_volumeTerm_3xSerP1 =
  &Gkyl::MaxwellModDecl<3,1,G_SERENDIPITY_C>::volumeTerm;

__device__ Maxwell_surfTerm_t p_Maxwell_surfTerm_3xSerP1 =
  &Gkyl::MaxwellModDecl<3,1,G_SERENDIPITY_C>::surfTerm;


Maxwell_volumeTerm_t
Maxwell_getVolumeTerm(int cdim, int polyOrder, int basisType) {
  Maxwell_volumeTerm_t func;
  
  if (cdim == 1) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_volumeTerm_1xSerP1, sizeof(Maxwell_volumeTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_volumeTerm_1xSerP2, sizeof(Maxwell_volumeTerm_t));
      }
  }
  else if (cdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_volumeTerm_2xSerP1, sizeof(Maxwell_volumeTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_volumeTerm_2xSerP2, sizeof(Maxwell_volumeTerm_t));
      }      
  }
  else if (cdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_volumeTerm_3xSerP1, sizeof(Maxwell_volumeTerm_t));
      }
  }
  return func;
}

Maxwell_surfTerm_t
Maxwell_getSurfTerm(int cdim, int polyOrder, int basisType) {
  Maxwell_surfTerm_t func;
  
  if (cdim == 1) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_surfTerm_1xSerP1, sizeof(Maxwell_surfTerm_t));
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_surfTerm_1xSerP2, sizeof(Maxwell_surfTerm_t));
      }
  }
  else if (cdim == 2) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_surfTerm_2xSerP1, sizeof(Maxwell_surfTerm_t));        
      }
      else if (polyOrder == 2) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_surfTerm_2xSerP2, sizeof(Maxwell_surfTerm_t));
      }      
  }
  else if (cdim == 3) {
      if (polyOrder == 1) {
        cudaMemcpyFromSymbol(&func, p_Maxwell_surfTerm_3xSerP1, sizeof(Maxwell_surfTerm_t));
      }
  }
  return func;
}
