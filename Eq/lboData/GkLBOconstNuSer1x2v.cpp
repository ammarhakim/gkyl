#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x2vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f[9]:      Input distribution function. 
  // out[9]:    Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0]   = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[0]   = rdv2[0]; 
  rdv2[1]   = 2.0/dxv[2]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 
  rdv2[1]   = rdv2[1]; 

  double alphaVpar[9]; 
  alphaVpar[0] = rdv2[0]*(2.0*nuUSum[0]-2.828427124746191*w[1]*nuSum); 
  alphaVpar[1] = 2.0*rdv2[0]*nuUSum[1]; 
  alphaVpar[2] = -1.632993161855453*nuSum; 

  double facVpar[2]; 
  facVpar[0] = nuVtSqSum[0]*rdvSq4[0]; 
  facVpar[1] = rdvSq4[0]*nuVtSqSum[1]; 

  double alphaMu[9]; 
  alphaMu[0] = (-5.656854249492382*rdv2[1]*w[2]*nuSum)+2.828427124746191*BmagInv[1]*nuVtSqSum[1]*rdv2[1]*m_+2.828427124746191*BmagInv[0]*nuVtSqSum[0]*rdv2[1]*m_; 
  alphaMu[1] = 2.828427124746191*BmagInv[0]*nuVtSqSum[1]*rdv2[1]*m_+2.828427124746191*nuVtSqSum[0]*BmagInv[1]*rdv2[1]*m_; 
  alphaMu[3] = -3.265986323710906*nuSum; 

  out[2] += 0.6123724356957944*(alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.6123724356957944*(alphaMu[3]*f[3]+alphaMu[1]*f[1]+alphaMu[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphaVpar[2]*f[4]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[5] += 0.6123724356957944*(alphaMu[3]*f[5]+alphaMu[0]*f[1]+f[0]*alphaMu[1]); 
  out[6] += 0.6123724356957944*((alphaMu[3]+alphaVpar[2])*f[6]+alphaVpar[1]*f[5]+alphaMu[1]*f[4]+alphaVpar[0]*f[3]+alphaMu[0]*f[2]); 
  out[7] += 0.6123724356957944*((alphaMu[3]+alphaVpar[2])*f[7]+alphaVpar[0]*f[5]+alphaMu[0]*f[4]+alphaVpar[1]*f[3]+alphaMu[1]*f[2]); 
  out[8] += 1.224744871391589*alphaVpar[2]*f[8]+1.369306393762915*(alphaVpar[1]*f[4]+alphaVpar[0]*f[2]+f[0]*alphaVpar[2])+4.743416490252569*(f[1]*facVpar[1]+f[0]*facVpar[0]); 

  return std::abs(0.3535533905932737*alphaVpar[0]) + 2.0*rdv2[1]*w[2]*nuSum+1.333333333333333*(BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*rdvSq4[1]*w[2]*m_+0.9428090415820625*nuVtSqSum[0]*rdvSq4[0]; 

} 
double GkLBOconstNuVol1x2vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f[20]:      Input distribution function. 
  // out[20]:    Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0]   = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[0]   = rdv2[0]; 
  rdv2[1]   = 2.0/dxv[2]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 
  rdv2[1]   = rdv2[1]; 

  double alphaVpar[20]; 
  alphaVpar[0] = rdv2[0]*(2.0*nuUSum[0]-2.828427124746191*w[1]*nuSum); 
  alphaVpar[1] = 2.0*rdv2[0]*nuUSum[1]; 
  alphaVpar[2] = -1.632993161855453*nuSum; 
  alphaVpar[7] = 2.0*rdv2[0]*nuUSum[2]; 

  double facVpar[3]; 
  facVpar[0] = nuVtSqSum[0]*rdvSq4[0]; 
  facVpar[1] = rdvSq4[0]*nuVtSqSum[1]; 
  facVpar[2] = rdvSq4[0]*nuVtSqSum[2]; 

  double alphaMu[20]; 
  alphaMu[0] = (-5.656854249492382*rdv2[1]*w[2]*nuSum)+2.828427124746191*rdv2[1]*BmagInv[2]*nuVtSqSum[2]*m_+2.828427124746191*BmagInv[1]*nuVtSqSum[1]*rdv2[1]*m_+2.828427124746191*BmagInv[0]*nuVtSqSum[0]*rdv2[1]*m_; 
  alphaMu[1] = 2.529822128134704*BmagInv[1]*rdv2[1]*nuVtSqSum[2]*m_+2.529822128134704*nuVtSqSum[1]*rdv2[1]*BmagInv[2]*m_+2.828427124746191*BmagInv[0]*nuVtSqSum[1]*rdv2[1]*m_+2.828427124746191*nuVtSqSum[0]*BmagInv[1]*rdv2[1]*m_; 
  alphaMu[3] = -3.265986323710906*nuSum; 
  alphaMu[7] = 1.807015805810503*rdv2[1]*BmagInv[2]*nuVtSqSum[2]*m_+2.828427124746191*BmagInv[0]*rdv2[1]*nuVtSqSum[2]*m_+2.828427124746191*nuVtSqSum[0]*rdv2[1]*BmagInv[2]*m_+2.529822128134704*BmagInv[1]*nuVtSqSum[1]*rdv2[1]*m_; 

  double facMu[20]; 
  facMu[0] = 2.828427124746191*rdvSq4[1]*BmagInv[2]*nuVtSqSum[2]*w[2]*m_+2.828427124746191*BmagInv[1]*nuVtSqSum[1]*rdvSq4[1]*w[2]*m_+2.828427124746191*BmagInv[0]*nuVtSqSum[0]*rdvSq4[1]*w[2]*m_; 
  facMu[1] = 2.529822128134704*BmagInv[1]*rdvSq4[1]*nuVtSqSum[2]*w[2]*m_+2.529822128134704*nuVtSqSum[1]*rdvSq4[1]*BmagInv[2]*w[2]*m_+2.828427124746191*BmagInv[0]*nuVtSqSum[1]*rdvSq4[1]*w[2]*m_+2.828427124746191*nuVtSqSum[0]*BmagInv[1]*rdvSq4[1]*w[2]*m_; 
  facMu[3] = 1.632993161855453*rdv2[1]*BmagInv[2]*nuVtSqSum[2]*m_+1.632993161855453*BmagInv[1]*nuVtSqSum[1]*rdv2[1]*m_+1.632993161855453*BmagInv[0]*nuVtSqSum[0]*rdv2[1]*m_; 
  facMu[5] = 1.460593486680444*BmagInv[1]*rdv2[1]*nuVtSqSum[2]*m_+1.460593486680444*nuVtSqSum[1]*rdv2[1]*BmagInv[2]*m_+1.632993161855453*BmagInv[0]*nuVtSqSum[1]*rdv2[1]*m_+1.632993161855453*nuVtSqSum[0]*BmagInv[1]*rdv2[1]*m_; 
  facMu[7] = 1.807015805810503*rdvSq4[1]*BmagInv[2]*nuVtSqSum[2]*w[2]*m_+2.828427124746191*BmagInv[0]*rdvSq4[1]*nuVtSqSum[2]*w[2]*m_+2.828427124746191*nuVtSqSum[0]*rdvSq4[1]*BmagInv[2]*w[2]*m_+2.529822128134704*BmagInv[1]*nuVtSqSum[1]*rdvSq4[1]*w[2]*m_; 
  facMu[13] = 1.043281061914602*rdv2[1]*BmagInv[2]*nuVtSqSum[2]*m_+1.632993161855453*BmagInv[0]*rdv2[1]*nuVtSqSum[2]*m_+1.632993161855453*nuVtSqSum[0]*rdv2[1]*BmagInv[2]*m_+1.460593486680443*BmagInv[1]*nuVtSqSum[1]*rdv2[1]*m_; 

  out[2] += 0.6123724356957944*(alphaVpar[7]*f[7]+alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.6123724356957944*(alphaMu[7]*f[7]+alphaMu[3]*f[3]+alphaMu[1]*f[1]+alphaMu[0]*f[0]); 
  out[4] += 0.5477225575051661*(alphaVpar[1]*f[7]+f[1]*alphaVpar[7])+0.6123724356957944*(alphaVpar[2]*f[4]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[5] += 0.5477225575051661*(alphaMu[1]*f[7]+f[1]*alphaMu[7])+0.6123724356957944*(alphaMu[3]*f[5]+alphaMu[0]*f[1]+f[0]*alphaMu[1]); 
  out[6] += 0.6123724356957944*(alphaVpar[7]*f[13]+alphaMu[7]*f[11]+(alphaMu[3]+alphaVpar[2])*f[6]+alphaVpar[1]*f[5]+alphaMu[1]*f[4]+alphaVpar[0]*f[3]+alphaMu[0]*f[2]); 
  out[8] += 1.369306393762915*alphaVpar[7]*f[11]+1.224744871391589*alphaVpar[2]*f[8]+4.743416490252569*facVpar[2]*f[7]+1.369306393762915*(alphaVpar[1]*f[4]+alphaVpar[0]*f[2]+f[0]*alphaVpar[2])+4.743416490252569*(f[1]*facVpar[1]+f[0]*facVpar[0]); 
  out[9] += f[13]*(2.371708245126284*facMu[13]+1.369306393762915*alphaMu[7])+1.224744871391589*alphaMu[3]*f[9]+2.371708245126284*f[7]*facMu[7]+f[5]*(2.371708245126284*facMu[5]+1.369306393762915*alphaMu[1])+2.371708245126284*f[3]*facMu[3]+1.369306393762915*(alphaMu[0]*f[3]+f[0]*alphaMu[3])+2.371708245126284*(f[1]*facMu[1]+f[0]*facMu[0]); 
  out[10] += 0.5477225575051661*(alphaVpar[1]*f[13]+alphaMu[1]*f[11])+0.6123724356957944*(alphaMu[3]+alphaVpar[2])*f[10]+0.5477225575051661*(f[5]*alphaVpar[7]+f[4]*alphaMu[7])+0.6123724356957944*(alphaVpar[0]*f[5]+alphaMu[0]*f[4]+alphaVpar[1]*f[3]+alphaMu[1]*f[2]); 
  out[11] += 0.6123724356957944*alphaVpar[2]*f[11]+0.3912303982179757*alphaVpar[7]*f[7]+0.6123724356957944*(alphaVpar[0]*f[7]+f[0]*alphaVpar[7])+0.5477225575051661*alphaVpar[1]*f[1]; 
  out[12] += 1.224744871391589*(alphaVpar[2]*f[12]+alphaVpar[1]*f[11])+4.242640687119286*facVpar[1]*f[7]+f[4]*(1.224744871391589*alphaVpar[7]+1.369306393762915*alphaVpar[0])+4.242640687119286*f[1]*facVpar[2]+1.369306393762915*(alphaVpar[1]*f[2]+f[1]*alphaVpar[2])+4.743416490252569*(f[0]*facVpar[1]+facVpar[0]*f[1]); 
  out[13] += 0.6123724356957944*alphaMu[3]*f[13]+0.3912303982179757*alphaMu[7]*f[7]+0.6123724356957944*(alphaMu[0]*f[7]+f[0]*alphaMu[7])+0.5477225575051661*alphaMu[1]*f[1]; 
  out[14] += 1.369306393762915*alphaVpar[7]*f[17]+(0.6123724356957944*alphaMu[3]+1.224744871391589*alphaVpar[2])*f[14]+4.743416490252569*facVpar[2]*f[13]+0.6123724356957944*alphaMu[1]*f[12]+1.369306393762915*alphaVpar[1]*f[10]+0.6123724356957944*alphaMu[0]*f[8]+1.369306393762915*alphaVpar[0]*f[6]+4.743416490252569*facVpar[1]*f[5]+(1.369306393762915*alphaVpar[2]+4.743416490252569*facVpar[0])*f[3]; 
  out[15] += 1.224744871391589*alphaMu[3]*f[15]+2.121320343559642*f[5]*facMu[13]+(2.121320343559642*facMu[5]+1.224744871391589*alphaMu[1])*f[13]+2.121320343559642*(f[1]*facMu[7]+facMu[1]*f[7])+1.224744871391589*f[5]*alphaMu[7]+2.371708245126284*(f[3]*facMu[5]+facMu[3]*f[5])+1.369306393762915*(alphaMu[0]*f[5]+alphaMu[1]*f[3]+f[1]*alphaMu[3])+2.371708245126284*(f[0]*facMu[1]+facMu[0]*f[1]); 
  out[16] += (2.371708245126284*facMu[13]+1.369306393762915*alphaMu[7])*f[17]+1.224744871391589*alphaMu[3]*f[16]+0.6123724356957944*(alphaVpar[2]*f[16]+alphaVpar[1]*f[15])+2.371708245126284*facMu[7]*f[11]+(2.371708245126284*facMu[5]+1.369306393762915*alphaMu[1])*f[10]+0.6123724356957944*alphaVpar[0]*f[9]+(2.371708245126284*facMu[3]+1.369306393762915*alphaMu[0])*f[6]+2.371708245126284*facMu[1]*f[4]+f[2]*(1.369306393762915*alphaMu[3]+2.371708245126284*facMu[0]); 
  out[17] += 0.6123724356957944*(alphaMu[3]+alphaVpar[2])*f[17]+(0.3912303982179757*alphaVpar[7]+0.6123724356957944*alphaVpar[0])*f[13]+0.3912303982179757*alphaMu[7]*f[11]+0.6123724356957944*(alphaMu[0]*f[11]+f[3]*alphaVpar[7]+f[2]*alphaMu[7])+0.5477225575051661*(alphaVpar[1]*f[5]+alphaMu[1]*f[4]); 
  out[18] += 0.6123724356957944*alphaMu[3]*f[18]+1.224744871391589*(alphaVpar[2]*f[18]+alphaVpar[1]*f[17])+4.242640687119286*facVpar[1]*f[13]+(0.5477225575051661*alphaMu[7]+0.6123724356957944*alphaMu[0])*f[12]+(1.224744871391589*alphaVpar[7]+1.369306393762915*alphaVpar[0])*f[10]+0.6123724356957944*alphaMu[1]*f[8]+1.369306393762915*alphaVpar[1]*f[6]+(4.242640687119286*facVpar[2]+1.369306393762915*alphaVpar[2])*f[5]+4.743416490252569*(facVpar[0]*f[5]+facVpar[1]*f[3]); 
  out[19] += (1.224744871391589*alphaMu[3]+0.6123724356957944*alphaVpar[2])*f[19]+(2.121320343559642*facMu[5]+1.224744871391589*alphaMu[1])*f[17]+(0.5477225575051661*alphaVpar[7]+0.6123724356957944*alphaVpar[0])*f[15]+2.121320343559642*(f[10]*facMu[13]+facMu[1]*f[11])+(1.224744871391589*alphaMu[7]+2.371708245126284*facMu[3]+1.369306393762915*alphaMu[0])*f[10]+0.6123724356957944*alphaVpar[1]*f[9]+2.121320343559642*f[4]*facMu[7]+2.371708245126284*facMu[5]*f[6]+1.369306393762915*(alphaMu[1]*f[6]+alphaMu[3]*f[4])+2.371708245126284*(facMu[0]*f[4]+facMu[1]*f[2]); 

  return std::abs(0.3535533905932737*alphaVpar[0]-0.3952847075210473*alphaVpar[7]) + 2.0*rdv2[1]*w[2]*nuSum-0.7115124735378848*facMu[7]-1.42302494707577*facVpar[2]+1.272792206135784*facVpar[0]+0.6363961030678922*facMu[0]; 

} 
