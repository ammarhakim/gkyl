#include "AlphaGenGeoModDecl.h" 
#include <math.h> 
void AlphaGenGeoSer1x3vP1(const double *w, const double *dxv, const double *tvComp, double const *gxx, double const *gxy, double const *gyy, double const *gxz, double const *gyz, double const *gzz, double *alphaGeo) 
{ 
  // w[NDIM]:    Cell-center coordinates.
  // dxv[NDIM]:  Cell spacing.
  // tvComp[18]: Components for tangent basis vectors.
  // gij[2]:    Contravariant components of metric tensor.
  // alpha:      Output alpha field.


  alphaGeo[0] = 4.0*w[1]; 
  alphaGeo[2] = 1.154700538379252*dxv[1]; 

} 
