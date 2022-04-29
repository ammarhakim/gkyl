#include "AlphaGenGeoModDecl.h" 
#include <math.h> 
void AlphaGenGeoSer1x3vP1(const double *w, const double *dxv, const double *tvComp, double const *gxx, double const *gxy, double const *gxz, double const *gyy, double const *gyz, double const *gzz, double *alphaGeo) 
{ 
  // w[NDIM]:    Cell-center coordinates.
  // dxv[NDIM]:  Cell spacing.
  // tvComp[18]: Components for tangent basis vectors.
  // gij[2]:    Contravariant components of metric tensor.
  // alpha:      Output alpha field.


  alphaGeo[0] = 2.0*(w[3]*(gxz[1]*tvComp[9]+gxz[0]*tvComp[8])+(gxy[1]*w[3]+gxz[1]*w[2])*tvComp[7]+(gxy[0]*w[3]+gxz[0]*w[2])*tvComp[6]+(gxx[1]*w[3]+gxy[1]*w[2]+gxz[1]*w[1])*tvComp[5]+(gxx[0]*w[3]+gxy[0]*w[2]+gxz[0]*w[1])*tvComp[4]+(gxx[1]*w[2]+gxy[1]*w[1])*tvComp[3]+gxx[0]*tvComp[2]*w[2]+w[1]*(gxy[0]*tvComp[2]+gxx[1]*tvComp[1]+gxx[0]*tvComp[0])); 
  alphaGeo[1] = 2.0*(w[3]*(gxz[0]*tvComp[9]+gxz[1]*tvComp[8])+(gxy[0]*w[3]+gxz[0]*w[2])*tvComp[7]+(gxy[1]*w[3]+gxz[1]*w[2])*tvComp[6]+(gxx[0]*w[3]+gxy[0]*w[2]+gxz[0]*w[1])*tvComp[5]+(gxx[1]*w[3]+gxy[1]*w[2]+gxz[1]*w[1])*tvComp[4]+(gxx[0]*w[2]+gxy[0]*w[1])*tvComp[3]+gxx[1]*tvComp[2]*w[2]+w[1]*(gxy[1]*tvComp[2]+gxx[0]*tvComp[1]+tvComp[0]*gxx[1])); 
  alphaGeo[2] = 0.5773502691896258*dxv[1]*(gxz[1]*tvComp[5]+gxz[0]*tvComp[4]+gxy[1]*tvComp[3]+gxy[0]*tvComp[2]+gxx[1]*tvComp[1]+gxx[0]*tvComp[0]); 
  alphaGeo[3] = 0.5773502691896258*dxv[2]*(gxz[1]*tvComp[7]+gxz[0]*tvComp[6]+gxy[1]*tvComp[5]+gxy[0]*tvComp[4]+gxx[1]*tvComp[3]+gxx[0]*tvComp[2]); 
  alphaGeo[4] = 0.5773502691896258*dxv[3]*(gxz[1]*tvComp[9]+gxz[0]*tvComp[8]+gxy[1]*tvComp[7]+gxy[0]*tvComp[6]+gxx[1]*tvComp[5]+gxx[0]*tvComp[4]); 
  alphaGeo[5] = 0.5773502691896258*dxv[1]*(gxz[0]*tvComp[5]+gxz[1]*tvComp[4]+gxy[0]*tvComp[3]+gxy[1]*tvComp[2]+gxx[0]*tvComp[1]+tvComp[0]*gxx[1]); 
  alphaGeo[6] = 0.5773502691896258*dxv[2]*(gxz[0]*tvComp[7]+gxz[1]*tvComp[6]+gxy[0]*tvComp[5]+gxy[1]*tvComp[4]+gxx[0]*tvComp[3]+gxx[1]*tvComp[2]); 
  alphaGeo[8] = 0.5773502691896258*dxv[3]*(gxz[0]*tvComp[9]+gxz[1]*tvComp[8]+gxy[0]*tvComp[7]+gxy[1]*tvComp[6]+gxx[0]*tvComp[5]+gxx[1]*tvComp[4]); 

} 
