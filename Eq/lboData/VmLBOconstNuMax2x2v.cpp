#include <VmLBOModDecl.h> 
double VmLBOconstNuVol2x2vMaxP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0] = 2.0/dxv[2]; 
  rdvSq4[0] = 4.0/(dxv[2]*dxv[2]); 
  rdv2[1] = 2.0/dxv[3]; 
  rdvSq4[1] = 4.0/(dxv[3]*dxv[3]); 

  double drBar[6]; 
  drBar[0] = nuU[0]-2.0*w[2]*nu; 
  drBar[1] = nuU[1]; 
  drBar[2] = nuU[2]; 
  drBar[3] = nuU[3]-2.0*w[3]*nu; 
  drBar[4] = nuU[4]; 
  drBar[5] = nuU[5]; 

  out[3] += rdv2[0]*(0.8660254037844386*f[2]*drBar[2]+0.8660254037844386*f[1]*drBar[1]+0.8660254037844386*f[0]*drBar[0])-1.0*f[3]*nu; 
  out[4] += rdv2[1]*(0.8660254037844386*f[2]*drBar[5]+0.8660254037844386*f[1]*drBar[4]+0.8660254037844386*f[0]*drBar[3])-1.0*f[4]*nu; 

  double nuVtSqP[0]; 
  nuVtSqP[0] = 2.0*nuVtSq[0]; 
  nuVtSqP[1] = 2.0*nuVtSq[1]; 
  nuVtSqP[2] = 2.0*nuVtSq[2]; 
  const double nuVtSqMid = 0.25*nuVtSqP[0]; 
return nuVtSqMid*(rdvSq4[0]+rdvSq4[1]); 

} 
