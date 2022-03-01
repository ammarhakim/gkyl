// Gkyl ------------------------------------------------------------------------
//
// C header file for the Vlasov equation kernels.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef ALPHAGENGEO_MOD_DECL_H 
#define ALPHAGENGEO_MOD_DECL_H 

#include <cmath> 

extern "C" {

    void AlphaGenGeoSer1x3vP1(const double *w, const double *dxv, const double *tvComp, double const *gxx, double const *gxy, double const *gxz, double const *gyy, double const *gyz, double const *gzz, const double *jacobGeo, double *alphaGeo);

  void AlphaGenGeoSer3x3vP1(const double *w, const double *dxv, const double *tvComp, double const *gxx, double const *gxy, double const *gxz, double const *gyy, double const *gyz, double const *gzz, const double *jacobGeo, double *alphaGeo);

 
} 
#endif 
