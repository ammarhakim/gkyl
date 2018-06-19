#include <VmLBOModDecl.h> 
double VmLBOconstNuVol2x3vMaxP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[3]; 
  double rdvSq4[3]; 
  rdv2[0] = 2.0/dxv[2]; 
  rdvSq4[0] = 4.0/(dxv[2]*dxv[2]); 
  rdv2[1] = 2.0/dxv[3]; 
  rdvSq4[1] = 4.0/(dxv[3]*dxv[3]); 
  rdv2[2] = 2.0/dxv[4]; 
  rdvSq4[2] = 4.0/(dxv[4]*dxv[4]); 

  double drBar[9]; 
  drBar[0] = nuU[0]-2.0*w[2]*nu; 
  drBar[1] = nuU[1]; 
  drBar[2] = nuU[2]; 
  drBar[3] = nuU[3]-2.0*w[3]*nu; 
  drBar[4] = nuU[4]; 
  drBar[5] = nuU[5]; 
  drBar[6] = nuU[6]-2.0*w[4]*nu; 
  drBar[7] = nuU[7]; 
  drBar[8] = nuU[8]; 

  out[3] += rdv2[0]*(0.8660254037844386*drBar[2]*f[2]+0.8660254037844386*drBar[1]*f[1]+0.8660254037844386*drBar[0]*f[0])-1.0*f[3]*nu; 
  out[4] += rdv2[1]*(0.8660254037844386*f[2]*drBar[5]+0.8660254037844386*f[1]*drBar[4]+0.8660254037844386*f[0]*drBar[3])-1.0*f[4]*nu; 
  out[5] += rdv2[2]*(0.8660254037844386*f[2]*drBar[8]+0.8660254037844386*f[1]*drBar[7]+0.8660254037844386*f[0]*drBar[6])-1.0*f[5]*nu; 

return nu*(rdvSq4[0]+rdvSq4[1]+rdvSq4[2])*0.5; 

} 
double VmLBOconstNuVol2x3vMaxP2(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[3]; 
  double rdvSq4[3]; 
  rdv2[0] = 2.0/dxv[2]; 
  rdvSq4[0] = 4.0/(dxv[2]*dxv[2]); 
  rdv2[1] = 2.0/dxv[3]; 
  rdvSq4[1] = 4.0/(dxv[3]*dxv[3]); 
  rdv2[2] = 2.0/dxv[4]; 
  rdvSq4[2] = 4.0/(dxv[4]*dxv[4]); 

  double drBar[18]; 
  drBar[0] = nuU[0]-2.0*w[2]*nu; 
  drBar[1] = nuU[1]; 
  drBar[2] = nuU[2]; 
  drBar[3] = nuU[3]; 
  drBar[4] = nuU[4]; 
  drBar[5] = nuU[5]; 
  drBar[6] = nuU[6]-2.0*w[3]*nu; 
  drBar[7] = nuU[7]; 
  drBar[8] = nuU[8]; 
  drBar[9] = nuU[9]; 
  drBar[10] = nuU[10]; 
  drBar[11] = nuU[11]; 
  drBar[12] = nuU[12]-2.0*w[4]*nu; 
  drBar[13] = nuU[13]; 
  drBar[14] = nuU[14]; 
  drBar[15] = nuU[15]; 
  drBar[16] = nuU[16]; 
  drBar[17] = nuU[17]; 

  out[3] += rdv2[0]*(0.8660254037844386*drBar[5]*f[17]+0.8660254037844386*drBar[4]*f[16]+0.8660254037844386*drBar[3]*f[6]+0.8660254037844386*drBar[2]*f[2]+0.8660254037844386*drBar[1]*f[1]+0.8660254037844386*drBar[0]*f[0])-1.0*f[3]*nu; 
  out[4] += rdv2[1]*(0.8660254037844386*drBar[11]*f[17]+0.8660254037844386*drBar[10]*f[16]+0.8660254037844386*f[6]*drBar[9]+0.8660254037844386*f[2]*drBar[8]+0.8660254037844386*f[1]*drBar[7]+0.8660254037844386*f[0]*drBar[6])-1.0*f[4]*nu; 
  out[5] += rdv2[2]*(0.8660254037844386*drBar[17]*f[17]+0.8660254037844386*drBar[16]*f[16]+0.8660254037844386*f[6]*drBar[15]+0.8660254037844386*f[2]*drBar[14]+0.8660254037844386*f[1]*drBar[13]+0.8660254037844386*f[0]*drBar[12])-1.0*f[5]*nu; 
  out[7] += rdv2[0]*(0.7745966692414833*drBar[1]*f[16]+0.8660254037844386*drBar[2]*f[6]+0.7745966692414833*f[1]*drBar[4]+0.8660254037844386*f[2]*drBar[3]+0.8660254037844386*drBar[0]*f[1]+0.8660254037844386*f[0]*drBar[1])-1.0*f[7]*nu; 
  out[8] += rdv2[0]*(0.7745966692414833*drBar[2]*f[17]+0.8660254037844386*drBar[1]*f[6]+0.7745966692414833*f[2]*drBar[5]+0.8660254037844386*f[1]*drBar[3]+0.8660254037844386*drBar[0]*f[2]+0.8660254037844386*f[0]*drBar[2])-1.0*f[8]*nu; 
  out[9] += rdv2[1]*(0.7745966692414833*drBar[7]*f[16]+0.7745966692414833*f[1]*drBar[10]+0.8660254037844386*f[2]*drBar[9]+0.8660254037844386*f[6]*drBar[8]+0.8660254037844386*f[0]*drBar[7]+0.8660254037844386*f[1]*drBar[6])-1.0*f[9]*nu; 
  out[10] += rdv2[1]*(0.7745966692414833*drBar[8]*f[17]+0.7745966692414833*f[2]*drBar[11]+0.8660254037844386*f[1]*drBar[9]+0.8660254037844386*f[0]*drBar[8]+0.8660254037844386*f[6]*drBar[7]+0.8660254037844386*f[2]*drBar[6])-1.0*f[10]*nu; 
  out[11] += (-2.0*f[11]*nu)+rdv2[0]*(0.8660254037844386*drBar[2]*f[10]+0.8660254037844386*drBar[1]*f[9]+0.8660254037844386*drBar[0]*f[4])+rdv2[1]*(0.8660254037844386*drBar[8]*f[8]+0.8660254037844386*drBar[7]*f[7]+0.8660254037844386*f[3]*drBar[6]); 
  out[12] += rdv2[2]*(0.7745966692414833*drBar[13]*f[16]+0.7745966692414833*f[1]*drBar[16]+0.8660254037844386*f[2]*drBar[15]+0.8660254037844386*f[6]*drBar[14]+0.8660254037844386*f[0]*drBar[13]+0.8660254037844386*f[1]*drBar[12])-1.0*f[12]*nu; 
  out[13] += rdv2[2]*(0.7745966692414833*drBar[14]*f[17]+0.7745966692414833*f[2]*drBar[17]+0.8660254037844386*f[1]*drBar[15]+0.8660254037844386*f[0]*drBar[14]+0.8660254037844386*f[6]*drBar[13]+0.8660254037844386*f[2]*drBar[12])-1.0*f[13]*nu; 
  out[14] += (-2.0*f[14]*nu)+rdv2[2]*(0.8660254037844386*f[8]*drBar[14]+0.8660254037844386*f[7]*drBar[13]+0.8660254037844386*f[3]*drBar[12])+rdv2[0]*(0.8660254037844386*drBar[2]*f[13]+0.8660254037844386*drBar[1]*f[12]+0.8660254037844386*drBar[0]*f[5]); 
  out[15] += (-2.0*f[15]*nu)+rdv2[2]*(0.8660254037844386*f[10]*drBar[14]+0.8660254037844386*f[9]*drBar[13]+0.8660254037844386*f[4]*drBar[12])+rdv2[1]*(0.8660254037844386*drBar[8]*f[13]+0.8660254037844386*drBar[7]*f[12]+0.8660254037844386*f[5]*drBar[6]); 
  out[18] += (-2.0*f[18]*nu)-2.23606797749979*f[0]*nu+rdvSq4[0]*(3.354101966249685*nuVtSq[5]*f[17]+3.354101966249685*nuVtSq[4]*f[16]+3.354101966249685*nuVtSq[3]*f[6]+3.354101966249685*f[2]*nuVtSq[2]+3.354101966249685*f[1]*nuVtSq[1]+3.354101966249685*f[0]*nuVtSq[0])+rdv2[0]*(1.936491673103709*drBar[2]*f[8]+1.936491673103709*drBar[1]*f[7]+1.936491673103709*drBar[0]*f[3]); 
  out[19] += (-2.0*f[19]*nu)-2.23606797749979*f[0]*nu+rdvSq4[1]*(3.354101966249685*nuVtSq[5]*f[17]+3.354101966249685*nuVtSq[4]*f[16]+3.354101966249685*nuVtSq[3]*f[6]+3.354101966249685*f[2]*nuVtSq[2]+3.354101966249685*f[1]*nuVtSq[1]+3.354101966249685*f[0]*nuVtSq[0])+rdv2[1]*(1.936491673103709*drBar[8]*f[10]+1.936491673103709*drBar[7]*f[9]+1.936491673103709*f[4]*drBar[6]); 
  out[20] += (-2.0*f[20]*nu)-2.23606797749979*f[0]*nu+rdvSq4[2]*(3.354101966249685*nuVtSq[5]*f[17]+3.354101966249685*nuVtSq[4]*f[16]+3.354101966249685*nuVtSq[3]*f[6]+3.354101966249685*f[2]*nuVtSq[2]+3.354101966249685*f[1]*nuVtSq[1]+3.354101966249685*f[0]*nuVtSq[0])+rdv2[2]*(1.936491673103709*f[13]*drBar[14]+1.936491673103709*f[12]*drBar[13]+1.936491673103709*f[5]*drBar[12]); 

return nu*(rdvSq4[0]+rdvSq4[1]+rdvSq4[2])*0.5; 

} 
