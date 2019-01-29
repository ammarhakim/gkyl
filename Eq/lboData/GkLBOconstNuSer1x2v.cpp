#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x2vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[3]:    Cell-center coordinates. 
  // dxv[3]:  Cell spacing. 
  // nu:      diffusion coefficient (collisionality). 
  // u[2]:    bulk velocity. 
  // vtSq[2]: thermal speed squared. 
  // f[8]:    Input distribution function. 
  // out[8]:  Incremented output 
  double rdv2nu[2]; 
  double rdvSq4nu[2]; 
  rdv2nu[0] = 2.0/dxv[1]; 
  rdvSq4nu[0] = nu*rdv2nu[0]*rdv2nu[0]; 
  rdv2nu[0] = nu*rdv2nu[0]; 
  rdv2nu[1] = 2.0/dxv[2]; 
  rdvSq4nu[1] = nu*rdv2nu[1]*rdv2nu[1]; 
  rdv2nu[1] = nu*rdv2nu[1]; 

  double alphaVpar[8]; 
  alphaVpar[0] = rdv2nu[0]*(2.0*u[0]-2.828427124746191*w[1]); 
  alphaVpar[1] = 2.0*rdv2nu[0]*u[1]; 
  alphaVpar[2] = -1.632993161855453*nu; 

  double alphaMu[8]; 
  alphaMu[0] = 2.828427124746191*BmagInv[1]*rdv2nu[1]*vtSq[1]*m_+2.828427124746191*BmagInv[0]*vtSq[0]*rdv2nu[1]*m_-5.656854249492382*rdv2nu[1]*w[2]; 
  alphaMu[1] = 2.828427124746191*BmagInv[0]*rdv2nu[1]*vtSq[1]*m_+2.828427124746191*vtSq[0]*BmagInv[1]*rdv2nu[1]*m_; 
  alphaMu[3] = -3.265986323710906*nu; 

  out[2] += 0.6123724356957944*(alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.6123724356957944*(alphaMu[3]*f[3]+alphaMu[1]*f[1]+alphaMu[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphaVpar[2]*f[4]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[5] += 0.6123724356957944*(alphaMu[3]*f[5]+alphaMu[0]*f[1]+f[0]*alphaMu[1]); 
  out[6] += 0.6123724356957944*((alphaMu[3]+alphaVpar[2])*f[6]+alphaVpar[1]*f[5]+alphaMu[1]*f[4]+alphaVpar[0]*f[3]+alphaMu[0]*f[2]); 
  out[7] += 0.6123724356957944*((alphaMu[3]+alphaVpar[2])*f[7]+alphaVpar[0]*f[5]+alphaMu[0]*f[4]+alphaVpar[1]*f[3]+alphaMu[1]*f[2]); 

  return std::abs(0.3535533905932737*alphaVpar[0]) + 1.333333333333333*rdvSq4nu[1]*(BmagInv[1]*vtSq[1]+BmagInv[0]*vtSq[0])*w[2]*m_+2.0*rdv2nu[1]*w[2]+0.9428090415820625*rdvSq4nu[0]*vtSq[0]; 

} 
double GkLBOconstNuVol1x2vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[3]:    Cell-center coordinates. 
  // dxv[3]:  Cell spacing. 
  // nu:      diffusion coefficient (collisionality). 
  // u[3]:    bulk velocity. 
  // vtSq[3]: thermal speed squared. 
  // f[20]:    Input distribution function. 
  // out[20]:  Incremented output 
  double rdv2nu[2]; 
  double rdvSq4nu[2]; 
  rdv2nu[0] = 2.0/dxv[1]; 
  rdvSq4nu[0] = nu*rdv2nu[0]*rdv2nu[0]; 
  rdv2nu[0] = nu*rdv2nu[0]; 
  rdv2nu[1] = 2.0/dxv[2]; 
  rdvSq4nu[1] = nu*rdv2nu[1]*rdv2nu[1]; 
  rdv2nu[1] = nu*rdv2nu[1]; 

  double alphaVpar[20]; 
  alphaVpar[0] = rdv2nu[0]*(2.0*u[0]-2.828427124746191*w[1]); 
  alphaVpar[1] = 2.0*rdv2nu[0]*u[1]; 
  alphaVpar[2] = -1.632993161855453*nu; 
  alphaVpar[7] = 2.0*rdv2nu[0]*u[2]; 

  double facVpar[3]; 
  facVpar[0] = rdvSq4nu[0]*vtSq[0]; 
  facVpar[1] = rdvSq4nu[0]*vtSq[1]; 
  facVpar[2] = rdvSq4nu[0]*vtSq[2]; 

  double alphaMu[20]; 
  alphaMu[0] = 2.828427124746191*rdv2nu[1]*BmagInv[2]*vtSq[2]*m_+2.828427124746191*BmagInv[1]*rdv2nu[1]*vtSq[1]*m_+2.828427124746191*BmagInv[0]*vtSq[0]*rdv2nu[1]*m_-5.656854249492382*rdv2nu[1]*w[2]; 
  alphaMu[1] = 2.529822128134704*BmagInv[1]*rdv2nu[1]*vtSq[2]*m_+2.529822128134704*rdv2nu[1]*vtSq[1]*BmagInv[2]*m_+2.828427124746191*BmagInv[0]*rdv2nu[1]*vtSq[1]*m_+2.828427124746191*vtSq[0]*BmagInv[1]*rdv2nu[1]*m_; 
  alphaMu[3] = -3.265986323710906*nu; 
  alphaMu[7] = 1.807015805810503*rdv2nu[1]*BmagInv[2]*vtSq[2]*m_+2.828427124746191*BmagInv[0]*rdv2nu[1]*vtSq[2]*m_+2.828427124746191*vtSq[0]*rdv2nu[1]*BmagInv[2]*m_+2.529822128134704*BmagInv[1]*rdv2nu[1]*vtSq[1]*m_; 

  double facMu[20]; 
  facMu[0] = 2.828427124746191*rdvSq4nu[1]*BmagInv[2]*vtSq[2]*w[2]*m_+2.828427124746191*BmagInv[1]*rdvSq4nu[1]*vtSq[1]*w[2]*m_+2.828427124746191*BmagInv[0]*vtSq[0]*rdvSq4nu[1]*w[2]*m_; 
  facMu[1] = 2.529822128134704*BmagInv[1]*rdvSq4nu[1]*vtSq[2]*w[2]*m_+2.529822128134704*rdvSq4nu[1]*vtSq[1]*BmagInv[2]*w[2]*m_+2.828427124746191*BmagInv[0]*rdvSq4nu[1]*vtSq[1]*w[2]*m_+2.828427124746191*vtSq[0]*BmagInv[1]*rdvSq4nu[1]*w[2]*m_; 
  facMu[3] = 1.632993161855453*rdv2nu[1]*BmagInv[2]*vtSq[2]*m_+1.632993161855453*BmagInv[1]*rdv2nu[1]*vtSq[1]*m_+1.632993161855453*BmagInv[0]*vtSq[0]*rdv2nu[1]*m_; 
  facMu[5] = 1.460593486680444*BmagInv[1]*rdv2nu[1]*vtSq[2]*m_+1.460593486680444*rdv2nu[1]*vtSq[1]*BmagInv[2]*m_+1.632993161855453*BmagInv[0]*rdv2nu[1]*vtSq[1]*m_+1.632993161855453*vtSq[0]*BmagInv[1]*rdv2nu[1]*m_; 
  facMu[7] = 1.807015805810503*rdvSq4nu[1]*BmagInv[2]*vtSq[2]*w[2]*m_+2.828427124746191*BmagInv[0]*rdvSq4nu[1]*vtSq[2]*w[2]*m_+2.828427124746191*vtSq[0]*rdvSq4nu[1]*BmagInv[2]*w[2]*m_+2.529822128134704*BmagInv[1]*rdvSq4nu[1]*vtSq[1]*w[2]*m_; 
  facMu[13] = 1.043281061914602*rdv2nu[1]*BmagInv[2]*vtSq[2]*m_+1.632993161855453*BmagInv[0]*rdv2nu[1]*vtSq[2]*m_+1.632993161855453*vtSq[0]*rdv2nu[1]*BmagInv[2]*m_+1.460593486680443*BmagInv[1]*rdv2nu[1]*vtSq[1]*m_; 

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

  return std::abs(0.3535533905932737*alphaVpar[0]-0.3952847075210473*alphaVpar[7]) + (-0.7115124735378848*facMu[7])+2.0*rdv2nu[1]*w[2]-1.42302494707577*facVpar[2]+1.272792206135784*facVpar[0]+0.6363961030678922*facMu[0]; 

} 
double GkLBOconstNuVol1x2vSerP3(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[3]:    Cell-center coordinates. 
  // dxv[3]:  Cell spacing. 
  // nu:      diffusion coefficient (collisionality). 
  // u[4]:    bulk velocity. 
  // vtSq[4]: thermal speed squared. 
  // f[32]:    Input distribution function. 
  // out[32]:  Incremented output 
  double rdv2nu[2]; 
  double rdvSq4nu[2]; 
  rdv2nu[0] = 2.0/dxv[1]; 
  rdvSq4nu[0] = nu*rdv2nu[0]*rdv2nu[0]; 
  rdv2nu[0] = nu*rdv2nu[0]; 
  rdv2nu[1] = 2.0/dxv[2]; 
  rdvSq4nu[1] = nu*rdv2nu[1]*rdv2nu[1]; 
  rdv2nu[1] = nu*rdv2nu[1]; 

  double alphaVpar[32]; 
  alphaVpar[0] = rdv2nu[0]*(2.0*u[0]-2.828427124746191*w[1]); 
  alphaVpar[1] = 2.0*rdv2nu[0]*u[1]; 
  alphaVpar[2] = -1.632993161855453*nu; 
  alphaVpar[7] = 2.0*rdv2nu[0]*u[2]; 
  alphaVpar[17] = 2.0*rdv2nu[0]*u[3]; 

  double facVpar[4]; 
  facVpar[0] = rdvSq4nu[0]*vtSq[0]; 
  facVpar[1] = rdvSq4nu[0]*vtSq[1]; 
  facVpar[2] = rdvSq4nu[0]*vtSq[2]; 
  facVpar[3] = rdvSq4nu[0]*vtSq[3]; 

  double alphaMu[32]; 
  alphaMu[0] = 2.828427124746191*rdv2nu[1]*BmagInv[3]*vtSq[3]*m_+2.828427124746191*rdv2nu[1]*BmagInv[2]*vtSq[2]*m_+2.828427124746191*BmagInv[1]*rdv2nu[1]*vtSq[1]*m_+2.828427124746191*BmagInv[0]*vtSq[0]*rdv2nu[1]*m_-5.656854249492382*rdv2nu[1]*w[2]; 
  alphaMu[1] = 2.484236013632475*rdv2nu[1]*BmagInv[2]*vtSq[3]*m_+2.484236013632475*rdv2nu[1]*vtSq[2]*BmagInv[3]*m_+2.529822128134704*BmagInv[1]*rdv2nu[1]*vtSq[2]*m_+2.529822128134704*rdv2nu[1]*vtSq[1]*BmagInv[2]*m_+2.828427124746191*BmagInv[0]*rdv2nu[1]*vtSq[1]*m_+2.828427124746191*vtSq[0]*BmagInv[1]*rdv2nu[1]*m_; 
  alphaMu[3] = -3.265986323710906*nu; 
  alphaMu[7] = 1.686548085423136*rdv2nu[1]*BmagInv[3]*vtSq[3]*m_+2.484236013632475*BmagInv[1]*rdv2nu[1]*vtSq[3]*m_+2.484236013632475*rdv2nu[1]*vtSq[1]*BmagInv[3]*m_+1.807015805810503*rdv2nu[1]*BmagInv[2]*vtSq[2]*m_+2.828427124746191*BmagInv[0]*rdv2nu[1]*vtSq[2]*m_+2.828427124746191*vtSq[0]*rdv2nu[1]*BmagInv[2]*m_+2.529822128134704*BmagInv[1]*rdv2nu[1]*vtSq[1]*m_; 
  alphaMu[17] = 1.686548085423136*rdv2nu[1]*BmagInv[2]*vtSq[3]*m_+2.828427124746191*BmagInv[0]*rdv2nu[1]*vtSq[3]*m_+1.686548085423136*rdv2nu[1]*vtSq[2]*BmagInv[3]*m_+2.828427124746191*vtSq[0]*rdv2nu[1]*BmagInv[3]*m_+2.484236013632475*BmagInv[1]*rdv2nu[1]*vtSq[2]*m_+2.484236013632475*rdv2nu[1]*vtSq[1]*BmagInv[2]*m_; 

  double facMu[32]; 
  facMu[0] = 2.828427124746191*rdvSq4nu[1]*w[2]*BmagInv[3]*vtSq[3]*m_+2.828427124746191*rdvSq4nu[1]*BmagInv[2]*vtSq[2]*w[2]*m_+2.828427124746191*BmagInv[1]*rdvSq4nu[1]*vtSq[1]*w[2]*m_+2.828427124746191*BmagInv[0]*vtSq[0]*rdvSq4nu[1]*w[2]*m_; 
  facMu[1] = 2.484236013632475*rdvSq4nu[1]*BmagInv[2]*w[2]*vtSq[3]*m_+2.484236013632475*rdvSq4nu[1]*vtSq[2]*w[2]*BmagInv[3]*m_+2.529822128134704*BmagInv[1]*rdvSq4nu[1]*vtSq[2]*w[2]*m_+2.529822128134704*rdvSq4nu[1]*vtSq[1]*BmagInv[2]*w[2]*m_+2.828427124746191*BmagInv[0]*rdvSq4nu[1]*vtSq[1]*w[2]*m_+2.828427124746191*vtSq[0]*BmagInv[1]*rdvSq4nu[1]*w[2]*m_; 
  facMu[3] = 1.632993161855453*rdv2nu[1]*BmagInv[3]*vtSq[3]*m_+1.632993161855453*rdv2nu[1]*BmagInv[2]*vtSq[2]*m_+1.632993161855453*BmagInv[1]*rdv2nu[1]*vtSq[1]*m_+1.632993161855453*BmagInv[0]*vtSq[0]*rdv2nu[1]*m_; 
  facMu[5] = 1.434274331201272*rdv2nu[1]*BmagInv[2]*vtSq[3]*m_+1.434274331201272*rdv2nu[1]*vtSq[2]*BmagInv[3]*m_+1.460593486680444*BmagInv[1]*rdv2nu[1]*vtSq[2]*m_+1.460593486680444*rdv2nu[1]*vtSq[1]*BmagInv[2]*m_+1.632993161855453*BmagInv[0]*rdv2nu[1]*vtSq[1]*m_+1.632993161855453*vtSq[0]*BmagInv[1]*rdv2nu[1]*m_; 
  facMu[7] = 1.686548085423136*rdvSq4nu[1]*w[2]*BmagInv[3]*vtSq[3]*m_+2.484236013632475*BmagInv[1]*rdvSq4nu[1]*w[2]*vtSq[3]*m_+2.484236013632475*rdvSq4nu[1]*vtSq[1]*w[2]*BmagInv[3]*m_+1.807015805810503*rdvSq4nu[1]*BmagInv[2]*vtSq[2]*w[2]*m_+2.828427124746191*BmagInv[0]*rdvSq4nu[1]*vtSq[2]*w[2]*m_+2.828427124746191*vtSq[0]*rdvSq4nu[1]*BmagInv[2]*w[2]*m_+2.529822128134704*BmagInv[1]*rdvSq4nu[1]*vtSq[1]*w[2]*m_; 
  facMu[13] = 0.9737289911202957*rdv2nu[1]*BmagInv[3]*vtSq[3]*m_+1.434274331201273*BmagInv[1]*rdv2nu[1]*vtSq[3]*m_+1.434274331201273*rdv2nu[1]*vtSq[1]*BmagInv[3]*m_+1.043281061914602*rdv2nu[1]*BmagInv[2]*vtSq[2]*m_+1.632993161855453*BmagInv[0]*rdv2nu[1]*vtSq[2]*m_+1.632993161855453*vtSq[0]*rdv2nu[1]*BmagInv[2]*m_+1.460593486680443*BmagInv[1]*rdv2nu[1]*vtSq[1]*m_; 
  facMu[17] = 1.686548085423136*rdvSq4nu[1]*BmagInv[2]*w[2]*vtSq[3]*m_+2.828427124746191*BmagInv[0]*rdvSq4nu[1]*w[2]*vtSq[3]*m_+1.686548085423136*rdvSq4nu[1]*vtSq[2]*w[2]*BmagInv[3]*m_+2.828427124746191*vtSq[0]*rdvSq4nu[1]*w[2]*BmagInv[3]*m_+2.484236013632475*BmagInv[1]*rdvSq4nu[1]*vtSq[2]*w[2]*m_+2.484236013632475*rdvSq4nu[1]*vtSq[1]*BmagInv[2]*w[2]*m_; 
  facMu[25] = 0.9737289911202957*rdv2nu[1]*BmagInv[2]*vtSq[3]*m_+1.632993161855452*BmagInv[0]*rdv2nu[1]*vtSq[3]*m_+0.9737289911202957*rdv2nu[1]*vtSq[2]*BmagInv[3]*m_+1.632993161855452*vtSq[0]*rdv2nu[1]*BmagInv[3]*m_+1.434274331201273*BmagInv[1]*rdv2nu[1]*vtSq[2]*m_+1.434274331201273*rdv2nu[1]*vtSq[1]*BmagInv[2]*m_; 

  out[2] += 0.6123724356957944*(alphaVpar[17]*f[17]+alphaVpar[7]*f[7]+alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.6123724356957944*(alphaMu[17]*f[17]+alphaMu[7]*f[7]+alphaMu[3]*f[3]+alphaMu[1]*f[1]+alphaMu[0]*f[0]); 
  out[4] += 0.537852874200477*(alphaVpar[7]*f[17]+f[7]*alphaVpar[17])+0.5477225575051661*(alphaVpar[1]*f[7]+f[1]*alphaVpar[7])+0.6123724356957944*(alphaVpar[2]*f[4]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[5] += 0.537852874200477*(alphaMu[7]*f[17]+f[7]*alphaMu[17])+0.5477225575051661*(alphaMu[1]*f[7]+f[1]*alphaMu[7])+0.6123724356957944*(alphaMu[3]*f[5]+alphaMu[0]*f[1]+f[0]*alphaMu[1]); 
  out[6] += 0.6123724356957944*(alphaVpar[17]*f[25]+alphaMu[17]*f[23]+alphaVpar[7]*f[13]+alphaMu[7]*f[11]+(alphaMu[3]+alphaVpar[2])*f[6]+alphaVpar[1]*f[5]+alphaMu[1]*f[4]+alphaVpar[0]*f[3]+alphaMu[0]*f[2]); 
  out[8] += 1.369306393762915*alphaVpar[17]*f[23]+4.743416490252569*facVpar[3]*f[17]+1.369306393762915*alphaVpar[7]*f[11]+1.224744871391589*alphaVpar[2]*f[8]+4.743416490252569*facVpar[2]*f[7]+1.369306393762915*(alphaVpar[1]*f[4]+alphaVpar[0]*f[2]+f[0]*alphaVpar[2])+4.743416490252569*(f[1]*facVpar[1]+f[0]*facVpar[0]); 
  out[9] += f[25]*(2.371708245126284*facMu[25]+1.369306393762915*alphaMu[17])+2.371708245126284*f[17]*facMu[17]+f[13]*(2.371708245126284*facMu[13]+1.369306393762915*alphaMu[7])+1.224744871391589*alphaMu[3]*f[9]+2.371708245126284*f[7]*facMu[7]+f[5]*(2.371708245126284*facMu[5]+1.369306393762915*alphaMu[1])+2.371708245126284*f[3]*facMu[3]+1.369306393762915*(alphaMu[0]*f[3]+f[0]*alphaMu[3])+2.371708245126284*(f[1]*facMu[1]+f[0]*facMu[0]); 
  out[10] += 0.537852874200477*(alphaVpar[7]*f[25]+alphaMu[7]*f[23]+f[13]*alphaVpar[17]+f[11]*alphaMu[17])+0.5477225575051661*(alphaVpar[1]*f[13]+alphaMu[1]*f[11])+0.6123724356957944*(alphaMu[3]+alphaVpar[2])*f[10]+0.5477225575051661*(f[5]*alphaVpar[7]+f[4]*alphaMu[7])+0.6123724356957944*(alphaVpar[0]*f[5]+alphaMu[0]*f[4]+alphaVpar[1]*f[3]+alphaMu[1]*f[2]); 
  out[11] += 0.3651483716701108*alphaVpar[17]*f[17]+0.537852874200477*(alphaVpar[1]*f[17]+f[1]*alphaVpar[17])+0.6123724356957944*alphaVpar[2]*f[11]+0.3912303982179757*alphaVpar[7]*f[7]+0.6123724356957944*(alphaVpar[0]*f[7]+f[0]*alphaVpar[7])+0.5477225575051661*alphaVpar[1]*f[1]; 
  out[12] += 1.202675588605909*alphaVpar[7]*f[23]+4.166190448976479*facVpar[2]*f[17]+1.202675588605909*f[11]*alphaVpar[17]+1.224744871391589*(alphaVpar[2]*f[12]+alphaVpar[1]*f[11])+(4.166190448976479*facVpar[3]+4.242640687119286*facVpar[1])*f[7]+f[4]*(1.224744871391589*alphaVpar[7]+1.369306393762915*alphaVpar[0])+4.242640687119286*f[1]*facVpar[2]+1.369306393762915*(alphaVpar[1]*f[2]+f[1]*alphaVpar[2])+4.743416490252569*(f[0]*facVpar[1]+facVpar[0]*f[1]); 
  out[13] += 0.3651483716701108*alphaMu[17]*f[17]+0.537852874200477*(alphaMu[1]*f[17]+f[1]*alphaMu[17])+0.6123724356957944*alphaMu[3]*f[13]+0.3912303982179757*alphaMu[7]*f[7]+0.6123724356957944*(alphaMu[0]*f[7]+f[0]*alphaMu[7])+0.5477225575051661*alphaMu[1]*f[1]; 
  out[14] += 1.369306393762915*alphaVpar[17]*f[29]+4.743416490252569*facVpar[3]*f[25]+1.369306393762915*alphaVpar[7]*f[20]+(0.6123724356957944*alphaMu[3]+1.224744871391589*alphaVpar[2])*f[14]+4.743416490252569*facVpar[2]*f[13]+0.6123724356957944*alphaMu[1]*f[12]+1.369306393762915*alphaVpar[1]*f[10]+0.6123724356957944*alphaMu[0]*f[8]+1.369306393762915*alphaVpar[0]*f[6]+4.743416490252569*facVpar[1]*f[5]+(1.369306393762915*alphaVpar[2]+4.743416490252569*facVpar[0])*f[3]; 
  out[15] += 2.08309522448824*f[13]*facMu[25]+(2.08309522448824*facMu[13]+1.202675588605909*alphaMu[7])*f[25]+2.08309522448824*(f[7]*facMu[17]+facMu[7]*f[17])+1.202675588605909*f[13]*alphaMu[17]+1.224744871391589*alphaMu[3]*f[15]+2.121320343559642*f[5]*facMu[13]+(2.121320343559642*facMu[5]+1.224744871391589*alphaMu[1])*f[13]+2.121320343559642*(f[1]*facMu[7]+facMu[1]*f[7])+1.224744871391589*f[5]*alphaMu[7]+2.371708245126284*(f[3]*facMu[5]+facMu[3]*f[5])+1.369306393762915*(alphaMu[0]*f[5]+alphaMu[1]*f[3]+f[1]*alphaMu[3])+2.371708245126284*(f[0]*facMu[1]+facMu[0]*f[1]); 
  out[16] += (2.371708245126284*facMu[25]+1.369306393762915*alphaMu[17])*f[29]+2.371708245126284*facMu[17]*f[23]+(2.371708245126284*facMu[13]+1.369306393762915*alphaMu[7])*f[20]+1.224744871391589*alphaMu[3]*f[16]+0.6123724356957944*(alphaVpar[2]*f[16]+alphaVpar[1]*f[15])+2.371708245126284*facMu[7]*f[11]+(2.371708245126284*facMu[5]+1.369306393762915*alphaMu[1])*f[10]+0.6123724356957944*alphaVpar[0]*f[9]+(2.371708245126284*facMu[3]+1.369306393762915*alphaMu[0])*f[6]+2.371708245126284*facMu[1]*f[4]+f[2]*(1.369306393762915*alphaMu[3]+2.371708245126284*facMu[0]); 
  out[18] += 16.20185174601965*facVpar[3]*f[23]+1.837117307087383*alphaVpar[2]*f[18]+0.9354143466934851*alphaVpar[17]*f[17]+2.091650066335188*alphaVpar[1]*f[12]+16.20185174601965*facVpar[2]*f[11]+2.091650066335188*alphaVpar[0]*f[8]+0.9354143466934851*alphaVpar[7]*f[7]+16.20185174601965*facVpar[1]*f[4]+(2.806243040080455*alphaVpar[2]+16.20185174601965*facVpar[0])*f[2]+0.9354143466934851*(alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[19] += 8.100925873009823*(f[17]*facMu[25]+facMu[17]*f[25])+1.837117307087383*alphaMu[3]*f[19]+0.9354143466934851*alphaMu[17]*f[17]+(7.24568837309472*facMu[5]+2.091650066335188*alphaMu[1])*f[15]+8.100925873009823*(f[7]*facMu[13]+facMu[7]*f[13])+(7.24568837309472*facMu[3]+2.091650066335188*alphaMu[0])*f[9]+0.9354143466934851*alphaMu[7]*f[7]+8.100925873009823*(f[1]*facMu[5]+facMu[1]*f[5]+f[0]*facMu[3])+(2.806243040080455*alphaMu[3]+8.100925873009823*facMu[0])*f[3]+0.9354143466934851*(alphaMu[1]*f[1]+alphaMu[0]*f[0]); 
  out[20] += (0.3651483716701108*alphaVpar[17]+0.537852874200477*alphaVpar[1])*f[25]+(0.3651483716701108*alphaMu[17]+0.537852874200477*alphaMu[1])*f[23]+0.6123724356957944*(alphaMu[3]+alphaVpar[2])*f[20]+0.537852874200477*(f[5]*alphaVpar[17]+f[4]*alphaMu[17])+(0.3912303982179757*alphaVpar[7]+0.6123724356957944*alphaVpar[0])*f[13]+0.3912303982179757*alphaMu[7]*f[11]+0.6123724356957944*(alphaMu[0]*f[11]+f[3]*alphaVpar[7]+f[2]*alphaMu[7])+0.5477225575051661*(alphaVpar[1]*f[5]+alphaMu[1]*f[4]); 
  out[21] += 1.202675588605909*alphaVpar[7]*f[29]+4.166190448976479*facVpar[2]*f[25]+(0.6123724356957944*alphaMu[3]+1.224744871391589*alphaVpar[2])*f[21]+(1.202675588605909*alphaVpar[17]+1.224744871391589*alphaVpar[1])*f[20]+(4.166190448976479*facVpar[3]+4.242640687119286*facVpar[1])*f[13]+(0.5477225575051661*alphaMu[7]+0.6123724356957944*alphaMu[0])*f[12]+(1.224744871391589*alphaVpar[7]+1.369306393762915*alphaVpar[0])*f[10]+0.6123724356957944*alphaMu[1]*f[8]+1.369306393762915*alphaVpar[1]*f[6]+(4.242640687119286*facVpar[2]+1.369306393762915*alphaVpar[2])*f[5]+4.743416490252569*(facVpar[0]*f[5]+facVpar[1]*f[3]); 
  out[22] += (2.08309522448824*facMu[13]+1.202675588605909*alphaMu[7])*f[29]+2.08309522448824*(f[20]*facMu[25]+facMu[7]*f[23])+(1.224744871391589*alphaMu[3]+0.6123724356957944*alphaVpar[2])*f[22]+(1.202675588605909*alphaMu[17]+2.121320343559642*facMu[5]+1.224744871391589*alphaMu[1])*f[20]+2.08309522448824*f[11]*facMu[17]+(0.5477225575051661*alphaVpar[7]+0.6123724356957944*alphaVpar[0])*f[15]+2.121320343559642*(f[10]*facMu[13]+facMu[1]*f[11])+(1.224744871391589*alphaMu[7]+2.371708245126284*facMu[3]+1.369306393762915*alphaMu[0])*f[10]+0.6123724356957944*alphaVpar[1]*f[9]+2.121320343559642*f[4]*facMu[7]+2.371708245126284*facMu[5]*f[6]+1.369306393762915*(alphaMu[1]*f[6]+alphaMu[3]*f[4])+2.371708245126284*(facMu[0]*f[4]+facMu[1]*f[2]); 
  out[23] += 0.6123724356957944*alphaVpar[2]*f[23]+(0.3651483716701108*alphaVpar[7]+0.6123724356957944*alphaVpar[0])*f[17]+(0.3651483716701108*f[7]+0.6123724356957944*f[0])*alphaVpar[17]+0.537852874200477*(alphaVpar[1]*f[7]+f[1]*alphaVpar[7]); 
  out[24] += 1.837117307087383*alphaVpar[2]*f[24]+14.23024947075771*facVpar[2]*f[23]+0.8215838362577489*(alphaVpar[7]*f[17]+f[7]*alphaVpar[17])+(1.870828693386971*alphaVpar[7]+2.091650066335188*alphaVpar[0])*f[12]+(14.23024947075771*facVpar[3]+14.49137674618944*facVpar[1])*f[11]+2.091650066335188*alphaVpar[1]*f[8]+0.8366600265340755*(alphaVpar[1]*f[7]+f[1]*alphaVpar[7])+(14.49137674618944*facVpar[2]+2.806243040080455*alphaVpar[2])*f[4]+16.20185174601965*(facVpar[0]*f[4]+facVpar[1]*f[2])+0.9354143466934851*(alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[25] += 0.6123724356957944*alphaMu[3]*f[25]+(0.3651483716701108*alphaMu[7]+0.6123724356957944*alphaMu[0])*f[17]+(0.3651483716701108*f[7]+0.6123724356957944*f[0])*alphaMu[17]+0.537852874200477*(alphaMu[1]*f[7]+f[1]*alphaMu[7]); 
  out[26] += 16.20185174601965*facVpar[3]*f[29]+(0.6123724356957944*alphaMu[3]+1.837117307087383*alphaVpar[2])*f[26]+0.9354143466934851*alphaVpar[17]*f[25]+0.6123724356957944*alphaMu[1]*f[24]+2.091650066335188*alphaVpar[1]*f[21]+16.20185174601965*facVpar[2]*f[20]+0.6123724356957944*alphaMu[0]*f[18]+2.091650066335188*alphaVpar[0]*f[14]+0.9354143466934851*alphaVpar[7]*f[13]+16.20185174601965*facVpar[1]*f[10]+(2.806243040080455*alphaVpar[2]+16.20185174601965*facVpar[0])*f[6]+0.9354143466934851*(alphaVpar[1]*f[5]+alphaVpar[0]*f[3]); 
  out[27] += 1.837117307087383*alphaMu[3]*f[27]+7.115124735378852*(f[7]*facMu[25]+facMu[7]*f[25]+f[13]*facMu[17]+facMu[13]*f[17])+0.8215838362577489*(alphaMu[7]*f[17]+f[7]*alphaMu[17])+(6.480740698407861*facMu[13]+1.870828693386971*alphaMu[7]+7.24568837309472*facMu[3]+2.091650066335188*alphaMu[0])*f[15]+7.24568837309472*(f[1]*facMu[13]+facMu[1]*f[13])+(7.24568837309472*facMu[5]+2.091650066335188*alphaMu[1])*f[9]+7.24568837309472*(f[5]*facMu[7]+facMu[5]*f[7])+0.8366600265340755*(alphaMu[1]*f[7]+f[1]*alphaMu[7])+8.100925873009823*f[0]*facMu[5]+2.806243040080455*alphaMu[3]*f[5]+8.100925873009823*(facMu[0]*f[5]+f[1]*facMu[3]+facMu[1]*f[3])+0.9354143466934851*(alphaMu[0]*f[1]+f[0]*alphaMu[1]); 
  out[28] += 8.100925873009823*facMu[17]*f[29]+1.837117307087383*alphaMu[3]*f[28]+0.6123724356957944*(alphaVpar[2]*f[28]+alphaVpar[1]*f[27])+f[23]*(8.100925873009823*facMu[25]+0.9354143466934851*alphaMu[17])+(7.24568837309472*facMu[5]+2.091650066335188*alphaMu[1])*f[22]+8.100925873009823*facMu[7]*f[20]+0.6123724356957944*alphaVpar[0]*f[19]+(7.24568837309472*facMu[3]+2.091650066335188*alphaMu[0])*f[16]+f[11]*(8.100925873009823*facMu[13]+0.9354143466934851*alphaMu[7])+8.100925873009823*facMu[1]*f[10]+(2.806243040080455*alphaMu[3]+8.100925873009823*facMu[0])*f[6]+f[4]*(8.100925873009823*facMu[5]+0.9354143466934851*alphaMu[1])+f[2]*(8.100925873009823*facMu[3]+0.9354143466934851*alphaMu[0]); 
  out[29] += 0.6123724356957944*(alphaMu[3]+alphaVpar[2])*f[29]+(0.3651483716701108*alphaVpar[7]+0.6123724356957944*alphaVpar[0])*f[25]+(0.3651483716701108*alphaMu[7]+0.6123724356957944*alphaMu[0])*f[23]+(0.3651483716701108*f[13]+0.6123724356957944*f[3])*alphaVpar[17]+(0.3651483716701108*f[11]+0.6123724356957944*f[2])*alphaMu[17]+0.537852874200477*(alphaVpar[1]*f[13]+alphaMu[1]*f[11]+f[5]*alphaVpar[7]+f[4]*alphaMu[7]); 
  out[30] += (0.6123724356957944*alphaMu[3]+1.837117307087383*alphaVpar[2])*f[30]+14.23024947075771*facVpar[2]*f[29]+0.8215838362577489*alphaVpar[7]*f[25]+(0.5477225575051661*alphaMu[7]+0.6123724356957944*alphaMu[0])*f[24]+(1.870828693386971*alphaVpar[7]+2.091650066335188*alphaVpar[0])*f[21]+(14.23024947075771*facVpar[3]+14.49137674618944*facVpar[1])*f[20]+0.6123724356957944*alphaMu[1]*f[18]+0.8215838362577489*f[13]*alphaVpar[17]+alphaVpar[1]*(2.091650066335188*f[14]+0.8366600265340755*f[13])+(14.49137674618944*facVpar[2]+2.806243040080455*alphaVpar[2]+16.20185174601965*facVpar[0])*f[10]+0.8366600265340755*f[5]*alphaVpar[7]+16.20185174601965*facVpar[1]*f[6]+0.9354143466934851*(alphaVpar[0]*f[5]+alphaVpar[1]*f[3]); 
  out[31] += (1.837117307087383*alphaMu[3]+0.6123724356957944*alphaVpar[2])*f[31]+7.115124735378852*facMu[7]*f[29]+(0.5477225575051661*alphaVpar[7]+0.6123724356957944*alphaVpar[0])*f[27]+7.115124735378852*f[11]*facMu[25]+(7.115124735378852*facMu[13]+0.8215838362577489*alphaMu[7])*f[23]+(6.480740698407861*facMu[13]+1.870828693386971*alphaMu[7]+7.24568837309472*facMu[3]+2.091650066335188*alphaMu[0])*f[22]+(7.115124735378852*facMu[17]+7.24568837309472*facMu[1])*f[20]+0.6123724356957944*alphaVpar[1]*f[19]+0.8215838362577489*f[11]*alphaMu[17]+(7.24568837309472*facMu[5]+2.091650066335188*alphaMu[1])*f[16]+7.24568837309472*f[4]*facMu[13]+(7.24568837309472*facMu[5]+0.8366600265340755*alphaMu[1])*f[11]+(7.24568837309472*facMu[7]+2.806243040080455*alphaMu[3]+8.100925873009823*facMu[0])*f[10]+0.8366600265340755*f[4]*alphaMu[7]+8.100925873009823*(facMu[1]*f[6]+f[2]*facMu[5]+facMu[3]*f[4])+0.9354143466934851*(alphaMu[0]*f[4]+alphaMu[1]*f[2]); 

  return std::abs(0.3535533905932737*alphaVpar[0]-0.3952847075210473*alphaVpar[7]) + (-0.9035079029052505*facMu[7])+2.0*rdv2nu[1]*w[2]-1.807015805810502*facVpar[2]+1.616244071283536*facVpar[0]+0.8081220356417678*facMu[0]; 

} 
