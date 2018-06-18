#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x2vMaxP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
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

  out[2] += rdv2[0]*(1.224744871391589*drBar[1]*f[1]+1.224744871391589*drBar[0]*f[0])-1.0*f[2]*nu; 
  out[3] += rdv2[1]*(1.224744871391589*f[1]*drBar[3]+1.224744871391589*f[0]*drBar[2])-1.0*f[3]*nu; 

return nu*(rdvSq4[0]+rdvSq4[1])*0.5; 

} 
double VmLBOconstNuVol1x2vMaxP2(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = 4.0/(dxv[1]*dxv[1]); 
  rdv2[1] = 2.0/dxv[2]; 
  rdvSq4[1] = 4.0/(dxv[2]*dxv[2]); 

  double drBar[6]; 
  drBar[0] = nuU[0]-1.414213562373095*w[1]*nu; 
  drBar[1] = nuU[1]; 
  drBar[2] = nuU[2]; 
  drBar[3] = nuU[3]-1.414213562373095*w[2]*nu; 
  drBar[4] = nuU[4]; 
  drBar[5] = nuU[5]; 

  out[2] += rdv2[0]*(1.224744871391589*drBar[2]*f[7]+1.224744871391589*drBar[1]*f[1]+1.224744871391589*drBar[0]*f[0])-1.0*f[2]*nu; 
  out[3] += rdv2[1]*(1.224744871391589*drBar[5]*f[7]+1.224744871391589*f[1]*drBar[4]+1.224744871391589*f[0]*drBar[3])-1.0*f[3]*nu; 
  out[4] += rdv2[0]*(1.095445115010332*drBar[1]*f[7]+1.095445115010332*f[1]*drBar[2]+1.224744871391589*drBar[0]*f[1]+1.224744871391589*f[0]*drBar[1])-1.0*f[4]*nu; 
  out[5] += rdv2[1]*(1.095445115010332*drBar[4]*f[7]+1.095445115010332*f[1]*drBar[5]+1.224744871391589*f[0]*drBar[4]+1.224744871391589*f[1]*drBar[3])-1.0*f[5]*nu; 
  out[6] += (-2.0*f[6]*nu)+rdv2[0]*(1.224744871391589*drBar[1]*f[5]+1.224744871391589*drBar[0]*f[3])+rdv2[1]*(1.224744871391589*drBar[4]*f[4]+1.224744871391589*f[2]*drBar[3]); 
  out[8] += (-2.0*f[8]*nu)-2.23606797749979*f[0]*nu+rdvSq4[0]*(4.743416490252569*nuVtSq[2]*f[7]+4.743416490252569*f[1]*nuVtSq[1]+4.743416490252569*f[0]*nuVtSq[0])+rdv2[0]*(2.738612787525831*drBar[1]*f[4]+2.738612787525831*drBar[0]*f[2]); 
  out[9] += (-2.0*f[9]*nu)-2.23606797749979*f[0]*nu+rdvSq4[1]*(4.743416490252569*nuVtSq[2]*f[7]+4.743416490252569*f[1]*nuVtSq[1]+4.743416490252569*f[0]*nuVtSq[0])+rdv2[1]*(2.738612787525831*drBar[4]*f[5]+2.738612787525831*drBar[3]*f[3]); 

return nu*(rdvSq4[0]+rdvSq4[1])*0.5; 

} 
double VmLBOconstNuVol1x2vMaxP3(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = 4.0/(dxv[1]*dxv[1]); 
  rdv2[1] = 2.0/dxv[2]; 
  rdvSq4[1] = 4.0/(dxv[2]*dxv[2]); 

  double drBar[8]; 
  drBar[0] = nuU[0]-1.414213562373095*w[1]*nu; 
  drBar[1] = nuU[1]; 
  drBar[2] = nuU[2]; 
  drBar[3] = nuU[3]; 
  drBar[4] = nuU[4]-1.414213562373095*w[2]*nu; 
  drBar[5] = nuU[5]; 
  drBar[6] = nuU[6]; 
  drBar[7] = nuU[7]; 

  out[2] += rdv2[0]*(1.224744871391589*drBar[3]*f[17]+1.224744871391589*drBar[2]*f[7]+1.224744871391589*drBar[1]*f[1]+1.224744871391589*drBar[0]*f[0])-1.0*f[2]*nu; 
  out[3] += rdv2[1]*(1.224744871391589*drBar[7]*f[17]+1.224744871391589*drBar[6]*f[7]+1.224744871391589*f[1]*drBar[5]+1.224744871391589*f[0]*drBar[4])-1.0*f[3]*nu; 
  out[4] += rdv2[0]*(1.075705748400954*drBar[2]*f[17]+1.075705748400954*drBar[3]*f[7]+1.095445115010332*drBar[1]*f[7]+1.095445115010332*f[1]*drBar[2]+1.224744871391589*drBar[0]*f[1]+1.224744871391589*f[0]*drBar[1])-1.0*f[4]*nu; 
  out[5] += rdv2[1]*(1.075705748400954*drBar[6]*f[17]+1.075705748400954*drBar[7]*f[7]+1.095445115010332*drBar[5]*f[7]+1.095445115010332*f[1]*drBar[6]+1.224744871391589*f[0]*drBar[5]+1.224744871391589*f[1]*drBar[4])-1.0*f[5]*nu; 
  out[6] += (-2.0*f[6]*nu)+rdv2[0]*(1.224744871391589*drBar[2]*f[13]+1.224744871391589*drBar[1]*f[5]+1.224744871391589*drBar[0]*f[3])+rdv2[1]*(1.224744871391589*drBar[6]*f[11]+1.224744871391589*f[4]*drBar[5]+1.224744871391589*f[2]*drBar[4]); 
  out[8] += (-2.0*f[8]*nu)-2.23606797749979*f[0]*nu+rdvSq4[0]*(4.743416490252569*nuVtSq[3]*f[17]+4.743416490252569*nuVtSq[2]*f[7]+4.743416490252569*f[1]*nuVtSq[1]+4.743416490252569*f[0]*nuVtSq[0])+rdv2[0]*(2.738612787525831*drBar[2]*f[11]+2.738612787525831*drBar[1]*f[4]+2.738612787525831*drBar[0]*f[2]); 
  out[9] += (-2.0*f[9]*nu)-2.23606797749979*f[0]*nu+rdvSq4[1]*(4.743416490252569*nuVtSq[3]*f[17]+4.743416490252569*nuVtSq[2]*f[7]+4.743416490252569*f[1]*nuVtSq[1]+4.743416490252569*f[0]*nuVtSq[0])+rdv2[1]*(2.738612787525831*drBar[6]*f[13]+2.738612787525831*drBar[5]*f[5]+2.738612787525831*f[3]*drBar[4]); 
  out[10] += (-2.0*f[10]*nu)+rdv2[0]*(1.075705748400954*drBar[3]*f[13]+1.095445115010332*drBar[1]*f[13]+1.095445115010332*drBar[2]*f[5]+1.224744871391589*drBar[0]*f[5]+1.224744871391589*drBar[1]*f[3])+rdv2[1]*(1.075705748400954*drBar[7]*f[11]+1.095445115010332*drBar[5]*f[11]+1.095445115010332*f[4]*drBar[6]+1.224744871391589*f[2]*drBar[5]+1.224744871391589*drBar[4]*f[4]); 
  out[11] += rdv2[0]*(0.7302967433402215*drBar[3]*f[17]+1.075705748400954*drBar[1]*f[17]+0.7824607964359517*drBar[2]*f[7]+1.224744871391589*drBar[0]*f[7]+1.075705748400954*f[1]*drBar[3]+1.224744871391589*f[0]*drBar[2]+1.095445115010332*drBar[1]*f[1])-1.0*f[11]*nu; 
  out[12] += (-2.0*f[12]*nu)-2.23606797749979*f[1]*nu+rdvSq4[0]*(4.16619044897648*nuVtSq[2]*f[17]+4.16619044897648*nuVtSq[3]*f[7]+4.242640687119286*nuVtSq[1]*f[7]+4.242640687119286*f[1]*nuVtSq[2]+4.743416490252569*f[0]*nuVtSq[1]+4.743416490252569*nuVtSq[0]*f[1])+rdv2[0]*(2.405351177211819*drBar[3]*f[11]+2.449489742783178*drBar[1]*f[11]+2.449489742783178*drBar[2]*f[4]+2.738612787525831*drBar[0]*f[4]+2.738612787525831*drBar[1]*f[2]); 
  out[13] += rdv2[1]*(0.7302967433402215*drBar[7]*f[17]+1.075705748400954*drBar[5]*f[17]+0.7824607964359517*drBar[6]*f[7]+1.224744871391589*drBar[4]*f[7]+1.075705748400954*f[1]*drBar[7]+1.224744871391589*f[0]*drBar[6]+1.095445115010332*f[1]*drBar[5])-1.0*f[13]*nu; 
  out[14] += (-3.0*f[14]*nu)-2.23606797749979*f[3]*nu+rdvSq4[0]*(4.743416490252569*nuVtSq[2]*f[13]+4.743416490252569*nuVtSq[1]*f[5]+4.743416490252569*nuVtSq[0]*f[3])+rdv2[1]*(1.224744871391589*drBar[5]*f[12]+1.224744871391589*drBar[4]*f[8])+rdv2[0]*(2.738612787525831*drBar[1]*f[10]+2.738612787525831*drBar[0]*f[6]); 
  out[15] += (-2.0*f[15]*nu)-2.23606797749979*f[1]*nu+rdvSq4[1]*(4.16619044897648*nuVtSq[2]*f[17]+4.16619044897648*nuVtSq[3]*f[7]+4.242640687119286*nuVtSq[1]*f[7]+4.242640687119286*f[1]*nuVtSq[2]+4.743416490252569*f[0]*nuVtSq[1]+4.743416490252569*nuVtSq[0]*f[1])+rdv2[1]*(2.405351177211819*drBar[7]*f[13]+2.449489742783178*drBar[5]*f[13]+2.449489742783178*f[5]*drBar[6]+2.738612787525831*drBar[4]*f[5]+2.738612787525831*f[3]*drBar[5]); 
  out[16] += (-3.0*f[16]*nu)-2.23606797749979*f[2]*nu+rdv2[0]*(1.224744871391589*drBar[1]*f[15]+1.224744871391589*drBar[0]*f[9])+rdvSq4[1]*(4.743416490252569*nuVtSq[2]*f[11]+4.743416490252569*nuVtSq[1]*f[4]+4.743416490252569*nuVtSq[0]*f[2])+rdv2[1]*(2.738612787525831*drBar[5]*f[10]+2.738612787525831*drBar[4]*f[6]); 
  out[18] += (-3.0*f[18]*nu)-4.58257569495584*f[2]*nu+rdv2[0]*(1.870828693386971*drBar[3]*f[17]+4.183300132670379*drBar[1]*f[12]+4.183300132670378*drBar[0]*f[8]+1.870828693386971*drBar[2]*f[7]+1.870828693386971*drBar[1]*f[1]+1.870828693386971*drBar[0]*f[0])+rdvSq4[0]*(16.20185174601965*nuVtSq[2]*f[11]+16.20185174601965*nuVtSq[1]*f[4]+16.20185174601965*nuVtSq[0]*f[2]); 
  out[19] += (-3.0*f[19]*nu)-4.58257569495584*f[3]*nu+rdv2[1]*(1.870828693386971*drBar[7]*f[17]+4.183300132670379*drBar[5]*f[15]+4.183300132670378*drBar[4]*f[9]+1.870828693386971*drBar[6]*f[7]+1.870828693386971*f[1]*drBar[5]+1.870828693386971*f[0]*drBar[4])+rdvSq4[1]*(16.20185174601965*nuVtSq[2]*f[13]+16.20185174601965*nuVtSq[1]*f[5]+16.20185174601965*nuVtSq[0]*f[3]); 

return nu*(rdvSq4[0]+rdvSq4[1])*0.5; 

} 
