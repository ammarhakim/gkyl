#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x1vSerP1(const double m_, const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 

  double drBar[2]; 
  drBar[0] = nuU[0]-1.414213562373095*w[1]*nu; 
  drBar[1] = nuU[1]; 

  out[2] += rdv2[0]*(1.224744871391589*f[1]*drBar[1]+1.224744871391589*f[0]*drBar[0])-1.0*f[2]*nu; 
  out[3] += rdv2[0]*(1.224744871391589*f[0]*drBar[1]+1.224744871391589*drBar[0]*f[1])-1.0*f[3]*nu; 

return nu*rdvSq4[0]*0.5; 

} 
double GkLBOconstNuVol1x1vSerP2(const double m_, const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 

  double drBar[3]; 
  drBar[0] = nuU[0]-1.414213562373095*w[1]*nu; 
  drBar[1] = nuU[1]; 
  drBar[2] = nuU[2]; 

  out[2] += rdv2[0]*(1.224744871391589*drBar[2]*f[4]+1.224744871391589*f[1]*drBar[1]+1.224744871391589*f[0]*drBar[0])-1.0*f[2]*nu; 
  out[3] += rdv2[0]*(1.095445115010332*drBar[1]*f[4]+1.095445115010332*f[1]*drBar[2]+1.224744871391589*f[0]*drBar[1]+1.224744871391589*drBar[0]*f[1])-1.0*f[3]*nu; 
  out[5] += ((-2.0*f[5])-2.23606797749979*f[0])*nu+rdv2[0]*(2.738612787525831*drBar[2]*f[6]+2.738612787525831*drBar[1]*f[3]+2.738612787525831*drBar[0]*f[2])+rdvSq4[0]*(4.743416490252569*nuVtSq[2]*f[4]+4.743416490252569*f[1]*nuVtSq[1]+4.743416490252569*f[0]*nuVtSq[0]); 
  out[6] += rdv2[0]*(0.7824607964359517*drBar[2]*f[4]+1.224744871391589*drBar[0]*f[4]+1.224744871391589*f[0]*drBar[2]+1.095445115010332*f[1]*drBar[1])-1.0*f[6]*nu; 
  out[7] += ((-2.0*f[7])-2.23606797749979*f[1])*nu+rdv2[0]*(2.449489742783178*drBar[1]*f[6]+2.449489742783178*drBar[2]*f[3]+2.738612787525831*drBar[0]*f[3]+2.738612787525831*drBar[1]*f[2])+rdvSq4[0]*(4.242640687119286*nuVtSq[1]*f[4]+4.242640687119286*f[1]*nuVtSq[2]+4.743416490252569*f[0]*nuVtSq[1]+4.743416490252569*nuVtSq[0]*f[1]); 

return nu*rdvSq4[0]*0.5; 

} 
