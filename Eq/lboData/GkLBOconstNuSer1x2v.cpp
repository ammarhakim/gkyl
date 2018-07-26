#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x2vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[1] = 2.0/dxv[2]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 

  double alpha0[8]; 
  alpha0[0] = 2.0*rdv2[0]*nuU[0]-2.828427124746191*rdv2[0]*w[1]*nu; 
  alpha0[1] = 2.0*rdv2[0]*nuU[1]; 
  alpha0[2] = -1.632993161855453*nu; 

  double alpha1[8]; 
  alpha1[0] = 2.0*nuVtSq[0]*rdvSq4[0]; 
  alpha1[1] = 2.0*rdvSq4[0]*nuVtSq[1]; 

  double mufac0[8]; 
  mufac0[0] = (-5.656854249492382*rdv2[1]*w[2]*nu)+2.828427124746191*BmagInv[1]*rdv2[1]*nuVtSq[1]*m_+2.828427124746191*BmagInv[0]*nuVtSq[0]*rdv2[1]*m_; 
  mufac0[1] = 2.828427124746191*BmagInv[0]*rdv2[1]*nuVtSq[1]*m_+2.828427124746191*nuVtSq[0]*BmagInv[1]*rdv2[1]*m_; 
  mufac0[3] = -3.265986323710906*nu; 

  double mufac1[8]; 
  mufac1[0] = 2.828427124746191*BmagInv[1]*nuVtSq[1]*rdvSq4[1]*w[2]*m_+2.828427124746191*BmagInv[0]*nuVtSq[0]*rdvSq4[1]*w[2]*m_; 
  mufac1[1] = 2.828427124746191*BmagInv[0]*nuVtSq[1]*rdvSq4[1]*w[2]*m_+2.828427124746191*nuVtSq[0]*BmagInv[1]*rdvSq4[1]*w[2]*m_; 
  mufac1[3] = 1.632993161855453*BmagInv[1]*rdv2[1]*nuVtSq[1]*m_+1.632993161855453*BmagInv[0]*nuVtSq[0]*rdv2[1]*m_; 
  mufac1[5] = 1.632993161855453*BmagInv[0]*rdv2[1]*nuVtSq[1]*m_+1.632993161855453*nuVtSq[0]*BmagInv[1]*rdv2[1]*m_; 

  out[2] += 0.6123724356957944*(alpha0[2]*f[2]+alpha0[1]*f[1]+alpha0[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[3]*mufac0[3]+f[1]*mufac0[1]+f[0]*mufac0[0]); 
  out[4] += 0.6123724356957944*(alpha0[2]*f[4]+alpha0[0]*f[1]+f[0]*alpha0[1]); 
  out[5] += 0.6123724356957944*(mufac0[3]*f[5]+f[0]*mufac0[1]+mufac0[0]*f[1]); 
  out[6] += 0.6123724356957944*((mufac0[3]+alpha0[2])*f[6]+alpha0[1]*f[5]+mufac0[1]*f[4]+alpha0[0]*f[3]+mufac0[0]*f[2]); 
  out[7] += 0.6123724356957944*((mufac0[3]+alpha0[2])*f[7]+alpha0[0]*f[5]+mufac0[0]*f[4]+alpha0[1]*f[3]+mufac0[1]*f[2]); 

  const double alpha1Mid = 0.3535533905932737*alpha1[0]; 
  const double mufac1Mid = 0.3535533905932737*mufac1[0]; 
  return alpha1Mid + mufac1Mid; 

} 
double GkLBOconstNuVol1x2vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[1] = 2.0/dxv[2]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 

  double alpha0[20]; 
  alpha0[0] = 2.0*rdv2[0]*nuU[0]-2.828427124746191*rdv2[0]*w[1]*nu; 
  alpha0[1] = 2.0*rdv2[0]*nuU[1]; 
  alpha0[2] = -1.632993161855453*nu; 
  alpha0[7] = 2.0*rdv2[0]*nuU[2]; 

  double alpha1[20]; 
  alpha1[0] = 2.0*nuVtSq[0]*rdvSq4[0]; 
  alpha1[1] = 2.0*rdvSq4[0]*nuVtSq[1]; 
  alpha1[7] = 2.0*rdvSq4[0]*nuVtSq[2]; 

  double mufac0[20]; 
  mufac0[0] = (-5.656854249492382*rdv2[1]*w[2]*nu)+2.828427124746191*rdv2[1]*BmagInv[2]*nuVtSq[2]*m_+2.828427124746191*BmagInv[1]*rdv2[1]*nuVtSq[1]*m_+2.828427124746191*BmagInv[0]*nuVtSq[0]*rdv2[1]*m_; 
  mufac0[1] = 2.529822128134704*BmagInv[1]*rdv2[1]*nuVtSq[2]*m_+2.529822128134704*rdv2[1]*nuVtSq[1]*BmagInv[2]*m_+2.828427124746191*BmagInv[0]*rdv2[1]*nuVtSq[1]*m_+2.828427124746191*nuVtSq[0]*BmagInv[1]*rdv2[1]*m_; 
  mufac0[3] = -3.265986323710906*nu; 
  mufac0[7] = 1.807015805810503*rdv2[1]*BmagInv[2]*nuVtSq[2]*m_+2.828427124746191*BmagInv[0]*rdv2[1]*nuVtSq[2]*m_+2.828427124746191*nuVtSq[0]*rdv2[1]*BmagInv[2]*m_+2.529822128134704*BmagInv[1]*rdv2[1]*nuVtSq[1]*m_; 

  double mufac1[20]; 
  mufac1[0] = 2.828427124746191*rdvSq4[1]*BmagInv[2]*w[2]*nuVtSq[2]*m_+2.828427124746191*BmagInv[1]*nuVtSq[1]*rdvSq4[1]*w[2]*m_+2.828427124746191*BmagInv[0]*nuVtSq[0]*rdvSq4[1]*w[2]*m_; 
  mufac1[1] = 2.529822128134704*BmagInv[1]*rdvSq4[1]*w[2]*nuVtSq[2]*m_+2.529822128134704*nuVtSq[1]*rdvSq4[1]*BmagInv[2]*w[2]*m_+2.828427124746191*BmagInv[0]*nuVtSq[1]*rdvSq4[1]*w[2]*m_+2.828427124746191*nuVtSq[0]*BmagInv[1]*rdvSq4[1]*w[2]*m_; 
  mufac1[3] = 1.632993161855453*rdv2[1]*BmagInv[2]*nuVtSq[2]*m_+1.632993161855453*BmagInv[1]*rdv2[1]*nuVtSq[1]*m_+1.632993161855453*BmagInv[0]*nuVtSq[0]*rdv2[1]*m_; 
  mufac1[5] = 1.460593486680444*BmagInv[1]*rdv2[1]*nuVtSq[2]*m_+1.460593486680444*rdv2[1]*nuVtSq[1]*BmagInv[2]*m_+1.632993161855453*BmagInv[0]*rdv2[1]*nuVtSq[1]*m_+1.632993161855453*nuVtSq[0]*BmagInv[1]*rdv2[1]*m_; 
  mufac1[7] = 1.807015805810503*rdvSq4[1]*BmagInv[2]*w[2]*nuVtSq[2]*m_+2.828427124746191*BmagInv[0]*rdvSq4[1]*w[2]*nuVtSq[2]*m_+2.828427124746191*nuVtSq[0]*rdvSq4[1]*BmagInv[2]*w[2]*m_+2.529822128134704*BmagInv[1]*nuVtSq[1]*rdvSq4[1]*w[2]*m_; 
  mufac1[13] = 1.043281061914602*rdv2[1]*BmagInv[2]*nuVtSq[2]*m_+1.632993161855453*BmagInv[0]*rdv2[1]*nuVtSq[2]*m_+1.632993161855453*nuVtSq[0]*rdv2[1]*BmagInv[2]*m_+1.460593486680443*BmagInv[1]*rdv2[1]*nuVtSq[1]*m_; 

  out[2] += 0.6123724356957944*(alpha0[7]*f[7]+alpha0[2]*f[2]+alpha0[1]*f[1]+alpha0[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*mufac0[7]+f[3]*mufac0[3]+f[1]*mufac0[1]+f[0]*mufac0[0]); 
  out[4] += 0.5477225575051661*(alpha0[1]*f[7]+f[1]*alpha0[7])+0.6123724356957944*(alpha0[2]*f[4]+alpha0[0]*f[1]+f[0]*alpha0[1]); 
  out[5] += 0.5477225575051661*(f[1]*mufac0[7]+mufac0[1]*f[7])+0.6123724356957944*(mufac0[3]*f[5]+f[0]*mufac0[1]+mufac0[0]*f[1]); 
  out[6] += 0.6123724356957944*(alpha0[7]*f[13]+mufac0[7]*f[11]+(mufac0[3]+alpha0[2])*f[6]+alpha0[1]*f[5]+mufac0[1]*f[4]+alpha0[0]*f[3]+mufac0[0]*f[2]); 
  out[8] += 1.369306393762915*alpha0[7]*f[11]+1.224744871391589*alpha0[2]*f[8]+2.371708245126284*alpha1[7]*f[7]+1.369306393762915*(alpha0[1]*f[4]+alpha0[0]*f[2]+f[0]*alpha0[2])+2.371708245126284*(alpha1[1]*f[1]+alpha1[0]*f[0]); 
  out[9] += f[13]*(2.371708245126284*mufac1[13]+1.369306393762915*mufac0[7])+1.224744871391589*mufac0[3]*f[9]+2.371708245126284*f[7]*mufac1[7]+f[5]*(2.371708245126284*mufac1[5]+1.369306393762915*mufac0[1])+2.371708245126284*f[3]*mufac1[3]+1.369306393762915*(f[0]*mufac0[3]+mufac0[0]*f[3])+2.371708245126284*(f[1]*mufac1[1]+f[0]*mufac1[0]); 
  out[10] += 0.5477225575051661*(alpha0[1]*f[13]+mufac0[1]*f[11])+0.6123724356957944*(mufac0[3]+alpha0[2])*f[10]+0.5477225575051661*(f[4]*mufac0[7]+f[5]*alpha0[7])+0.6123724356957944*(alpha0[0]*f[5]+mufac0[0]*f[4]+alpha0[1]*f[3]+mufac0[1]*f[2]); 
  out[11] += 0.6123724356957944*alpha0[2]*f[11]+0.3912303982179757*alpha0[7]*f[7]+0.6123724356957944*(alpha0[0]*f[7]+f[0]*alpha0[7])+0.5477225575051661*alpha0[1]*f[1]; 
  out[12] += 1.224744871391589*(alpha0[2]*f[12]+alpha0[1]*f[11])+2.121320343559642*(alpha1[1]*f[7]+f[1]*alpha1[7])+1.224744871391589*f[4]*alpha0[7]+1.369306393762915*(alpha0[0]*f[4]+alpha0[1]*f[2]+f[1]*alpha0[2])+2.371708245126284*(alpha1[0]*f[1]+f[0]*alpha1[1]); 
  out[13] += 0.6123724356957944*mufac0[3]*f[13]+0.3912303982179757*f[7]*mufac0[7]+0.6123724356957944*(f[0]*mufac0[7]+mufac0[0]*f[7])+0.5477225575051661*f[1]*mufac0[1]; 
  out[14] += 1.369306393762915*alpha0[7]*f[17]+(0.6123724356957944*mufac0[3]+1.224744871391589*alpha0[2])*f[14]+2.371708245126284*alpha1[7]*f[13]+0.6123724356957944*mufac0[1]*f[12]+1.369306393762915*alpha0[1]*f[10]+0.6123724356957944*mufac0[0]*f[8]+1.369306393762915*alpha0[0]*f[6]+2.371708245126284*alpha1[1]*f[5]+(1.369306393762915*alpha0[2]+2.371708245126284*alpha1[0])*f[3]; 
  out[15] += 1.224744871391589*mufac0[3]*f[15]+2.121320343559642*f[5]*mufac1[13]+(2.121320343559642*mufac1[5]+1.224744871391589*mufac0[1])*f[13]+2.121320343559642*f[1]*mufac1[7]+1.224744871391589*f[5]*mufac0[7]+2.121320343559642*mufac1[1]*f[7]+2.371708245126284*(f[3]*mufac1[5]+mufac1[3]*f[5])+1.369306393762915*(mufac0[0]*f[5]+f[1]*mufac0[3]+mufac0[1]*f[3])+2.371708245126284*(f[0]*mufac1[1]+mufac1[0]*f[1]); 
  out[16] += (2.371708245126284*mufac1[13]+1.369306393762915*mufac0[7])*f[17]+1.224744871391589*mufac0[3]*f[16]+0.6123724356957944*(alpha0[2]*f[16]+alpha0[1]*f[15])+2.371708245126284*mufac1[7]*f[11]+(2.371708245126284*mufac1[5]+1.369306393762915*mufac0[1])*f[10]+0.6123724356957944*alpha0[0]*f[9]+(2.371708245126284*mufac1[3]+1.369306393762915*mufac0[0])*f[6]+2.371708245126284*mufac1[1]*f[4]+f[2]*(1.369306393762915*mufac0[3]+2.371708245126284*mufac1[0]); 
  out[17] += 0.6123724356957944*(mufac0[3]+alpha0[2])*f[17]+(0.3912303982179757*alpha0[7]+0.6123724356957944*alpha0[0])*f[13]+0.3912303982179757*mufac0[7]*f[11]+0.6123724356957944*(mufac0[0]*f[11]+f[2]*mufac0[7]+f[3]*alpha0[7])+0.5477225575051661*(alpha0[1]*f[5]+mufac0[1]*f[4]); 
  out[18] += 0.6123724356957944*mufac0[3]*f[18]+1.224744871391589*(alpha0[2]*f[18]+alpha0[1]*f[17])+2.121320343559642*alpha1[1]*f[13]+(0.5477225575051661*mufac0[7]+0.6123724356957944*mufac0[0])*f[12]+(1.224744871391589*alpha0[7]+1.369306393762915*alpha0[0])*f[10]+0.6123724356957944*mufac0[1]*f[8]+2.121320343559642*f[5]*alpha1[7]+1.369306393762915*(alpha0[1]*f[6]+alpha0[2]*f[5])+2.371708245126284*(alpha1[0]*f[5]+alpha1[1]*f[3]); 
  out[19] += (1.224744871391589*mufac0[3]+0.6123724356957944*alpha0[2])*f[19]+(2.121320343559642*mufac1[5]+1.224744871391589*mufac0[1])*f[17]+(0.5477225575051661*alpha0[7]+0.6123724356957944*alpha0[0])*f[15]+2.121320343559642*(f[10]*mufac1[13]+mufac1[1]*f[11])+(1.224744871391589*mufac0[7]+2.371708245126284*mufac1[3]+1.369306393762915*mufac0[0])*f[10]+0.6123724356957944*alpha0[1]*f[9]+2.121320343559642*f[4]*mufac1[7]+2.371708245126284*mufac1[5]*f[6]+1.369306393762915*(mufac0[1]*f[6]+mufac0[3]*f[4])+2.371708245126284*(mufac1[0]*f[4]+mufac1[1]*f[2]); 

  const double alpha1Mid = 0.3535533905932737*alpha1[0]-0.3952847075210473*alpha1[7]; 
  const double mufac1Mid = 0.3535533905932737*mufac1[0]-0.3952847075210473*mufac1[7]; 
  return alpha1Mid + mufac1Mid; 

} 
