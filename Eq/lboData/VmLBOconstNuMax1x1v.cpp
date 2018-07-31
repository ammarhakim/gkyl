#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x1vMaxP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = 4.0/(dxv[1]*dxv[1]); 

  double drBar[2]; 
  drBar[0] = nuU[0]-1.414213562373095*w[1]*nu; 
  drBar[1] = nuU[1]; 

  out[2] += rdv2[0]*(1.224744871391589*f[1]*drBar[1]+1.224744871391589*f[0]*drBar[0])-1.0*f[2]*nu; 

  double nuVtSqP[0]; 
  nuVtSqP[0] = 1.414213562373095*nuVtSq[0]; 
  nuVtSqP[1] = 1.414213562373095*nuVtSq[1]; 
  const double nuVtSqMid = 0.5*nuVtSqP[0]; 
return nuVtSqMid*rdvSq4[0]; 

} 
