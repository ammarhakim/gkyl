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

  out[2] += rdv2[0]*(1.224744871391589*drBar[1]*f[1]+1.224744871391589*drBar[0]*f[0])-1.0*f[2]*nu; 

return nu*rdvSq4[0]*0.5; 

} 
double VmLBOconstNuVol1x1vMaxP2(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = 4.0/(dxv[1]*dxv[1]); 

  double drBar[3]; 
  drBar[0] = nuU[0]-1.414213562373095*w[1]*nu; 
  drBar[1] = nuU[1]; 
  drBar[2] = nuU[2]; 

  out[2] += rdv2[0]*(1.224744871391589*drBar[2]*f[4]+1.224744871391589*drBar[1]*f[1]+1.224744871391589*drBar[0]*f[0])-1.0*f[2]*nu; 
  out[3] += rdv2[0]*(1.095445115010332*drBar[1]*f[4]+1.095445115010332*f[1]*drBar[2]+1.224744871391589*drBar[0]*f[1]+1.224744871391589*f[0]*drBar[1])-1.0*f[3]*nu; 
  out[5] += (-2.0*f[5]*nu)-2.23606797749979*f[0]*nu+rdvSq4[0]*(4.743416490252569*nuVtSq[2]*f[4]+4.743416490252569*f[1]*nuVtSq[1]+4.743416490252569*f[0]*nuVtSq[0])+rdv2[0]*(2.738612787525831*drBar[1]*f[3]+2.738612787525831*drBar[0]*f[2]); 

return nu*rdvSq4[0]*0.5; 

} 
double VmLBOconstNuVol1x1vMaxP3(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = 4.0/(dxv[1]*dxv[1]); 

  double drBar[4]; 
  drBar[0] = nuU[0]-1.414213562373095*w[1]*nu; 
  drBar[1] = nuU[1]; 
  drBar[2] = nuU[2]; 
  drBar[3] = nuU[3]; 

  out[2] += rdv2[0]*(1.224744871391589*drBar[3]*f[8]+1.224744871391589*drBar[2]*f[4]+1.224744871391589*drBar[1]*f[1]+1.224744871391589*drBar[0]*f[0])-1.0*f[2]*nu; 
  out[3] += rdv2[0]*(1.075705748400954*drBar[2]*f[8]+1.075705748400954*drBar[3]*f[4]+1.095445115010332*drBar[1]*f[4]+1.095445115010332*f[1]*drBar[2]+1.224744871391589*drBar[0]*f[1]+1.224744871391589*f[0]*drBar[1])-1.0*f[3]*nu; 
  out[5] += (-2.0*f[5]*nu)-2.23606797749979*f[0]*nu+rdvSq4[0]*(4.743416490252569*nuVtSq[3]*f[8]+4.743416490252569*nuVtSq[2]*f[4]+4.743416490252569*f[1]*nuVtSq[1]+4.743416490252569*f[0]*nuVtSq[0])+rdv2[0]*(2.738612787525831*drBar[2]*f[6]+2.738612787525831*drBar[1]*f[3]+2.738612787525831*drBar[0]*f[2]); 
  out[6] += rdv2[0]*(0.7302967433402215*drBar[3]*f[8]+1.075705748400954*drBar[1]*f[8]+0.7824607964359517*drBar[2]*f[4]+1.224744871391589*drBar[0]*f[4]+1.075705748400954*f[1]*drBar[3]+1.224744871391589*f[0]*drBar[2]+1.095445115010332*drBar[1]*f[1])-1.0*f[6]*nu; 
  out[7] += (-2.0*f[7]*nu)-2.23606797749979*f[1]*nu+rdvSq4[0]*(4.16619044897648*nuVtSq[2]*f[8]+4.16619044897648*nuVtSq[3]*f[4]+4.242640687119286*nuVtSq[1]*f[4]+4.242640687119286*f[1]*nuVtSq[2]+4.743416490252569*f[0]*nuVtSq[1]+4.743416490252569*nuVtSq[0]*f[1])+rdv2[0]*(2.405351177211819*drBar[3]*f[6]+2.449489742783178*drBar[1]*f[6]+2.449489742783178*drBar[2]*f[3]+2.738612787525831*drBar[0]*f[3]+2.738612787525831*drBar[1]*f[2]); 
  out[9] += (-3.0*f[9]*nu)-4.58257569495584*f[2]*nu+rdv2[0]*(1.870828693386971*drBar[3]*f[8]+4.183300132670379*drBar[1]*f[7]+4.183300132670378*drBar[0]*f[5]+1.870828693386971*drBar[2]*f[4]+1.870828693386971*drBar[1]*f[1]+1.870828693386971*drBar[0]*f[0])+rdvSq4[0]*(16.20185174601965*nuVtSq[2]*f[6]+16.20185174601965*nuVtSq[1]*f[3]+16.20185174601965*nuVtSq[0]*f[2]); 

return nu*rdvSq4[0]*0.5; 

} 
