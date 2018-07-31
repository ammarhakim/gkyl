#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x2vSerP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = 4.0/(dxv[1]*dxv[1]); 
  rdv2[1] = 2.0/dxv[2]; 
  rdvSq4[1] = 4.0/(dxv[2]*dxv[2]); 

  double drBar[4]; 
  drBar[0] = nuU[0]-1.414213562373095*w[1]*nu; 
  drBar[1] = nuU[1]; 
  drBar[2] = nuU[2]-1.414213562373095*w[2]*nu; 
  drBar[3] = nuU[3]; 

  out[2] += rdv2[0]*(1.224744871391589*f[1]*drBar[1]+1.224744871391589*f[0]*drBar[0])-1.0*f[2]*nu; 
  out[3] += rdv2[1]*(1.224744871391589*f[1]*drBar[3]+1.224744871391589*f[0]*drBar[2])-1.0*f[3]*nu; 
  out[4] += rdv2[0]*(1.224744871391589*f[0]*drBar[1]+1.224744871391589*drBar[0]*f[1])-1.0*f[4]*nu; 
  out[5] += rdv2[1]*(1.224744871391589*f[0]*drBar[3]+1.224744871391589*f[1]*drBar[2])-1.0*f[5]*nu; 
  out[6] += (-2.0*f[6]*nu)+rdv2[0]*(1.224744871391589*drBar[1]*f[5]+1.224744871391589*drBar[0]*f[3])+rdv2[1]*(1.224744871391589*drBar[3]*f[4]+1.224744871391589*f[2]*drBar[2]); 
  out[7] += (-2.0*f[7]*nu)+rdv2[0]*(1.224744871391589*drBar[0]*f[5]+1.224744871391589*drBar[1]*f[3])+rdv2[1]*(1.224744871391589*drBar[2]*f[4]+1.224744871391589*f[2]*drBar[3]); 

  double nuVtSqP[0]; 
  nuVtSqP[0] = 2.0*nuVtSq[0]; 
  nuVtSqP[1] = 2.0*nuVtSq[1]; 
  const double nuVtSqMid = 0.3535533905932737*nuVtSqP[0]; 
return nuVtSqMid*(rdvSq4[0]+rdvSq4[1]); 

} 
