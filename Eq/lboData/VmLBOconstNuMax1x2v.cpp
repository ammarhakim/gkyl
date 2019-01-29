#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x2vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[1]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[1]*dxv[1]); 
  const double rdvy2nu = 2.0*nu/dxv[2]; 
  const double rdvySq4nu = 4.0*nu/(dxv[2]*dxv[2]); 

  double alphaDrag[8]; 
  // Expand rdv2nu*(vx-ux) in phase basis.
  alphaDrag[0] = (2.0*u[0]-2.828427124746191*w[1])*rdvx2nu; 
  alphaDrag[1] = 2.0*u[1]*rdvx2nu; 
  alphaDrag[2] = -0.8164965809277261*dxv[1]*rdvx2nu; 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alphaDrag[4] = (2.0*u[2]-2.828427124746191*w[2])*rdvy2nu; 
  alphaDrag[5] = 2.0*u[3]*rdvy2nu; 
  alphaDrag[7] = -0.8164965809277261*dxv[2]*rdvy2nu; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[3]*alphaDrag[7]+f[1]*alphaDrag[5]+f[0]*alphaDrag[4]); 

  return std::abs(0.1767766952966368*alphaDrag[0])+std::abs(0.1767766952966368*alphaDrag[4])+std::abs(0.9428090415820636*vtSq[0]*rdvxSq4nu)+std::abs(0.9428090415820636*vtSq[0]*rdvySq4nu); 

} 
double VmLBOconstNuVol1x2vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[1]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[1]*dxv[1]); 
  const double rdvy2nu = 2.0*nu/dxv[2]; 
  const double rdvySq4nu = 4.0*nu/(dxv[2]*dxv[2]); 

  double alphaDrag[20]; 
  // Expand rdv2nu*(vx-ux) in phase basis.
  alphaDrag[0] = (2.0*u[0]-2.828427124746191*w[1])*rdvx2nu; 
  alphaDrag[1] = 2.0*u[1]*rdvx2nu; 
  alphaDrag[2] = -0.8164965809277261*dxv[1]*rdvx2nu; 
  alphaDrag[7] = 2.0*u[2]*rdvx2nu; 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alphaDrag[10] = (2.0*u[3]-2.828427124746191*w[2])*rdvy2nu; 
  alphaDrag[11] = 2.0*u[4]*rdvy2nu; 
  alphaDrag[13] = -0.8164965809277261*dxv[2]*rdvy2nu; 
  alphaDrag[17] = 2.0*u[5]*rdvy2nu; 

  double facDiff[3]; 
  // Expand nu*vthSq in phase basis.
  facDiff[0] = vtSq[0]; 
  facDiff[1] = vtSq[1]; 
  facDiff[2] = vtSq[2]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alphaDrag[7]*f[7]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*alphaDrag[17]+f[3]*alphaDrag[13]+f[1]*alphaDrag[11]+f[0]*alphaDrag[10]); 
  out[4] += 0.5477225575051661*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.6123724356957944*(alphaDrag[2]*f[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.5477225575051661*f[1]*alphaDrag[17]+0.6123724356957944*f[5]*alphaDrag[13]+0.5477225575051661*f[7]*alphaDrag[11]+0.6123724356957944*(f[0]*alphaDrag[11]+f[1]*alphaDrag[10]); 
  out[6] += 0.6123724356957944*(f[6]*alphaDrag[13]+f[4]*alphaDrag[11]+f[2]*alphaDrag[10]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[8] += 4.743416490252569*(facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4nu+1.224744871391589*alphaDrag[2]*f[8]+1.369306393762915*(alphaDrag[1]*f[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 4.743416490252569*(facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4nu+1.224744871391589*f[9]*alphaDrag[13]+1.369306393762915*(f[0]*alphaDrag[13]+f[5]*alphaDrag[11]+f[3]*alphaDrag[10]); 

  return std::abs(0.1767766952966368*alphaDrag[0]-0.1976423537605236*alphaDrag[7])+std::abs(0.1767766952966368*alphaDrag[10]-0.1976423537605236*alphaDrag[17])+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvxSq4nu)+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvySq4nu); 

} 
double VmLBOconstNuVol1x2vMaxP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[1]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[1]*dxv[1]); 
  const double rdvy2nu = 2.0*nu/dxv[2]; 
  const double rdvySq4nu = 4.0*nu/(dxv[2]*dxv[2]); 

  double alphaDrag[40]; 
  // Expand rdv2nu*(vx-ux) in phase basis.
  alphaDrag[0] = (2.0*u[0]-2.828427124746191*w[1])*rdvx2nu; 
  alphaDrag[1] = 2.0*u[1]*rdvx2nu; 
  alphaDrag[2] = -0.8164965809277261*dxv[1]*rdvx2nu; 
  alphaDrag[7] = 2.0*u[2]*rdvx2nu; 
  alphaDrag[17] = 2.0*u[3]*rdvx2nu; 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alphaDrag[20] = (2.0*u[4]-2.828427124746191*w[2])*rdvy2nu; 
  alphaDrag[21] = 2.0*u[5]*rdvy2nu; 
  alphaDrag[23] = -0.8164965809277261*dxv[2]*rdvy2nu; 
  alphaDrag[27] = 2.0*u[6]*rdvy2nu; 
  alphaDrag[37] = 2.0*u[7]*rdvy2nu; 

  double facDiff[4]; 
  // Expand nu*vthSq in phase basis.
  facDiff[0] = vtSq[0]; 
  facDiff[1] = vtSq[1]; 
  facDiff[2] = vtSq[2]; 
  facDiff[3] = vtSq[3]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alphaDrag[17]*f[17]+alphaDrag[7]*f[7]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[17]*alphaDrag[37]+f[7]*alphaDrag[27]+f[3]*alphaDrag[23]+f[1]*alphaDrag[21]+f[0]*alphaDrag[20]); 
  out[4] += 0.537852874200477*(alphaDrag[7]*f[17]+f[7]*alphaDrag[17])+0.5477225575051661*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.6123724356957944*(alphaDrag[2]*f[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.537852874200477*f[7]*alphaDrag[37]+(0.537852874200477*f[17]+0.5477225575051661*f[1])*alphaDrag[27]+0.6123724356957944*f[5]*alphaDrag[23]+0.5477225575051661*f[7]*alphaDrag[21]+0.6123724356957944*(f[0]*alphaDrag[21]+f[1]*alphaDrag[20]); 
  out[6] += 0.6123724356957944*(f[11]*alphaDrag[27]+f[6]*alphaDrag[23]+f[4]*alphaDrag[21]+f[2]*alphaDrag[20]+alphaDrag[7]*f[13]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[8] += 4.743416490252569*(facDiff[3]*f[17]+facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4nu+1.369306393762915*alphaDrag[7]*f[11]+1.224744871391589*alphaDrag[2]*f[8]+1.369306393762915*(alphaDrag[1]*f[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 4.743416490252569*(facDiff[3]*f[17]+facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4nu+1.369306393762915*f[13]*alphaDrag[27]+1.224744871391589*f[9]*alphaDrag[23]+1.369306393762915*(f[0]*alphaDrag[23]+f[5]*alphaDrag[21]+f[3]*alphaDrag[20]); 
  out[10] += 0.537852874200477*f[11]*alphaDrag[37]+0.5477225575051661*f[4]*alphaDrag[27]+0.6123724356957944*f[10]*alphaDrag[23]+0.5477225575051661*f[11]*alphaDrag[21]+0.6123724356957944*(f[2]*alphaDrag[21]+f[4]*alphaDrag[20])+f[13]*(0.537852874200477*alphaDrag[17]+0.5477225575051661*alphaDrag[1])+0.6123724356957944*alphaDrag[2]*f[10]+0.5477225575051661*f[5]*alphaDrag[7]+0.6123724356957944*(alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[11] += 0.3651483716701108*alphaDrag[17]*f[17]+0.537852874200477*(alphaDrag[1]*f[17]+f[1]*alphaDrag[17])+0.6123724356957944*alphaDrag[2]*f[11]+0.3912303982179757*alphaDrag[7]*f[7]+0.6123724356957944*(alphaDrag[0]*f[7]+f[0]*alphaDrag[7])+0.5477225575051661*alphaDrag[1]*f[1]; 
  out[12] += (4.166190448976479*facDiff[2]*f[17]+4.242640687119286*(facDiff[1]*f[7]+f[1]*facDiff[2])+4.166190448976479*facDiff[3]*f[7]+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4nu+1.202675588605909*f[11]*alphaDrag[17]+1.224744871391589*(alphaDrag[2]*f[12]+alphaDrag[1]*f[11]+f[4]*alphaDrag[7])+1.369306393762915*(alphaDrag[0]*f[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[13] += (0.3651483716701108*f[17]+0.537852874200477*f[1])*alphaDrag[37]+0.3912303982179757*f[7]*alphaDrag[27]+0.6123724356957944*(f[0]*alphaDrag[27]+f[13]*alphaDrag[23])+(0.537852874200477*f[17]+0.5477225575051661*f[1])*alphaDrag[21]+0.6123724356957944*f[7]*alphaDrag[20]; 
  out[14] += 4.743416490252569*(facDiff[2]*f[13]+facDiff[1]*f[5]+facDiff[0]*f[3])*rdvxSq4nu+0.6123724356957944*(f[14]*alphaDrag[23]+f[12]*alphaDrag[21]+f[8]*alphaDrag[20])+1.224744871391589*alphaDrag[2]*f[14]+1.369306393762915*(alphaDrag[1]*f[10]+alphaDrag[0]*f[6]+alphaDrag[2]*f[3]); 
  out[15] += (4.166190448976479*facDiff[2]*f[17]+4.242640687119286*(facDiff[1]*f[7]+f[1]*facDiff[2])+4.166190448976479*facDiff[3]*f[7]+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvySq4nu+1.202675588605909*f[13]*alphaDrag[37]+1.224744871391589*f[5]*alphaDrag[27]+(1.224744871391589*f[15]+1.369306393762915*f[1])*alphaDrag[23]+1.224744871391589*f[13]*alphaDrag[21]+1.369306393762915*(f[3]*alphaDrag[21]+f[5]*alphaDrag[20]); 
  out[16] += 4.743416490252569*(facDiff[2]*f[11]+facDiff[1]*f[4]+facDiff[0]*f[2])*rdvySq4nu+1.224744871391589*f[16]*alphaDrag[23]+1.369306393762915*(f[2]*alphaDrag[23]+f[10]*alphaDrag[21]+f[6]*alphaDrag[20])+0.6123724356957944*(alphaDrag[2]*f[16]+alphaDrag[1]*f[15]+alphaDrag[0]*f[9]); 
  out[18] += 16.20185174601965*(facDiff[2]*f[11]+facDiff[1]*f[4]+facDiff[0]*f[2])*rdvxSq4nu+1.837117307087383*alphaDrag[2]*f[18]+0.9354143466934851*alphaDrag[17]*f[17]+2.091650066335188*(alphaDrag[1]*f[12]+alphaDrag[0]*f[8])+0.9354143466934851*alphaDrag[7]*f[7]+2.806243040080455*alphaDrag[2]*f[2]+0.9354143466934851*(alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[19] += 16.20185174601965*(facDiff[2]*f[13]+facDiff[1]*f[5]+facDiff[0]*f[3])*rdvySq4nu+0.9354143466934851*(f[17]*alphaDrag[37]+f[7]*alphaDrag[27])+(1.837117307087383*f[19]+2.806243040080455*f[3])*alphaDrag[23]+(2.091650066335188*f[15]+0.9354143466934851*f[1])*alphaDrag[21]+(2.091650066335188*f[9]+0.9354143466934851*f[0])*alphaDrag[20]; 

  return std::abs(0.1767766952966368*alphaDrag[0]-0.1976423537605236*alphaDrag[7])+std::abs(0.1767766952966368*alphaDrag[20]-0.1976423537605236*alphaDrag[27])+std::abs((1.616244071283538*facDiff[0]-1.807015805810503*facDiff[2])*rdvxSq4nu)+std::abs((1.616244071283538*facDiff[0]-1.807015805810503*facDiff[2])*rdvySq4nu); 

} 
