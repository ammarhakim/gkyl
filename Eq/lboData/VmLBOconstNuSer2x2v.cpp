#include <VmLBOModDecl.h> 
double VmLBOconstNuVol2x2vSerP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0] = 2.0/dxv[2]; 
  rdvSq4[0] = 4.0/(dxv[2]*dxv[2]); 
  rdv2[1] = 2.0/dxv[3]; 
  rdvSq4[1] = 4.0/(dxv[3]*dxv[3]); 

  double drBar[8]; 
  drBar[0] = nuU[0]-2.0*w[2]*nu; 
  drBar[1] = nuU[1]; 
  drBar[2] = nuU[2]; 
  drBar[3] = nuU[3]; 
  drBar[4] = nuU[4]-2.0*w[3]*nu; 
  drBar[5] = nuU[5]; 
  drBar[6] = nuU[6]; 
  drBar[7] = nuU[7]; 

  out[3] += rdv2[0]*(0.8660254037844386*drBar[3]*f[5]+0.8660254037844386*f[2]*drBar[2]+0.8660254037844386*f[1]*drBar[1]+0.8660254037844386*f[0]*drBar[0])-1.0*f[3]*nu; 
  out[4] += rdv2[1]*(0.8660254037844386*f[5]*drBar[7]+0.8660254037844386*f[2]*drBar[6]+0.8660254037844386*f[1]*drBar[5]+0.8660254037844386*f[0]*drBar[4])-1.0*f[4]*nu; 
  out[6] += rdv2[0]*(0.8660254037844386*drBar[2]*f[5]+0.8660254037844386*f[2]*drBar[3]+0.8660254037844386*f[0]*drBar[1]+0.8660254037844386*drBar[0]*f[1])-1.0*f[6]*nu; 
  out[7] += rdv2[0]*(0.8660254037844386*drBar[1]*f[5]+0.8660254037844386*f[1]*drBar[3]+0.8660254037844386*f[0]*drBar[2]+0.8660254037844386*drBar[0]*f[2])-1.0*f[7]*nu; 
  out[8] += rdv2[1]*(0.8660254037844386*f[2]*drBar[7]+0.8660254037844386*f[5]*drBar[6]+0.8660254037844386*f[0]*drBar[5]+0.8660254037844386*f[1]*drBar[4])-1.0*f[8]*nu; 
  out[9] += rdv2[1]*(0.8660254037844386*f[1]*drBar[7]+0.8660254037844386*f[0]*drBar[6]+0.8660254037844386*f[5]*drBar[5]+0.8660254037844386*f[2]*drBar[4])-1.0*f[9]*nu; 
  out[10] += (-2.0*f[10]*nu)+rdv2[0]*(0.8660254037844386*drBar[3]*f[12]+0.8660254037844386*drBar[2]*f[9]+0.8660254037844386*drBar[1]*f[8]+0.8660254037844386*drBar[0]*f[4])+rdv2[1]*(0.8660254037844386*drBar[7]*f[11]+0.8660254037844386*drBar[6]*f[7]+0.8660254037844386*drBar[5]*f[6]+0.8660254037844386*f[3]*drBar[4]); 
  out[11] += rdv2[0]*(0.8660254037844386*drBar[0]*f[5]+0.8660254037844386*f[0]*drBar[3]+0.8660254037844386*f[1]*drBar[2]+0.8660254037844386*drBar[1]*f[2])-1.0*f[11]*nu; 
  out[12] += rdv2[1]*(0.8660254037844386*f[0]*drBar[7]+0.8660254037844386*f[1]*drBar[6]+0.8660254037844386*f[2]*drBar[5]+0.8660254037844386*drBar[4]*f[5])-1.0*f[12]*nu; 
  out[13] += (-2.0*f[13]*nu)+rdv2[0]*(0.8660254037844386*drBar[2]*f[12]+0.8660254037844386*drBar[3]*f[9]+0.8660254037844386*drBar[0]*f[8]+0.8660254037844386*drBar[1]*f[4])+rdv2[1]*(0.8660254037844386*drBar[6]*f[11]+0.8660254037844386*f[7]*drBar[7]+0.8660254037844386*drBar[4]*f[6]+0.8660254037844386*f[3]*drBar[5]); 
  out[14] += (-2.0*f[14]*nu)+rdv2[0]*(0.8660254037844386*drBar[1]*f[12]+0.8660254037844386*drBar[0]*f[9]+0.8660254037844386*drBar[3]*f[8]+0.8660254037844386*drBar[2]*f[4])+rdv2[1]*(0.8660254037844386*drBar[5]*f[11]+0.8660254037844386*f[6]*drBar[7]+0.8660254037844386*drBar[4]*f[7]+0.8660254037844386*f[3]*drBar[6]); 
  out[15] += (-2.0*f[15]*nu)+rdv2[0]*(0.8660254037844386*drBar[0]*f[12]+0.8660254037844386*drBar[1]*f[9]+0.8660254037844386*drBar[2]*f[8]+0.8660254037844386*drBar[3]*f[4])+rdv2[1]*(0.8660254037844386*drBar[4]*f[11]+0.8660254037844386*f[3]*drBar[7]+0.8660254037844386*drBar[5]*f[7]+0.8660254037844386*f[6]*drBar[6]); 

  double nuVtSqP[0]; 
  nuVtSqP[0] = 2.0*nuVtSq[0]; 
  nuVtSqP[1] = 2.0*nuVtSq[1]; 
  nuVtSqP[2] = 2.0*nuVtSq[2]; 
  nuVtSqP[5] = 2.0*nuVtSq[3]; 
  const double nuVtSqMid = 0.25*nuVtSqP[0]; 
return nuVtSqMid*(rdvSq4[0]+rdvSq4[1]); 

} 
