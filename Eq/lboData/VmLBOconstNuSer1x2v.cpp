#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x2vSerP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdv20 = 2/dxv[1]; 
  const double rdvSq40 = 4/(dxv[1]*dxv[1]); 
  const double rdv21 = 2/dxv[2]; 
  const double rdvSq41 = 4/(dxv[2]*dxv[2]); 

  double alpha_mid = 0.0; 
  double alpha_drag[16]; 
  double alpha_diffusion[8]; 

  alpha_drag[0] = (2.0*nuU[0]-2.828427124746191*w[1]*nu)*rdv20; 
  alpha_drag[1] = 2.0*nuU[1]*rdv20; 
  alpha_drag[2] = -0.8164965809277261*dxv[1]*nu*rdv20; 
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[0]); 
  alpha_drag[8] = (2.0*nuU[8]-2.828427124746191*w[2]*nu)*rdv21; 
  alpha_drag[9] = 2.0*nuU[9]*rdv21; 
  alpha_drag[11] = -0.8164965809277261*dxv[2]*nu*rdv21; 
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[8]); 
  alpha_diffusion[0] = 2.0*nuVtSq[0]; 
  alpha_diffusion[1] = 2.0*nuVtSq[1]; 
  alpha_mid += std::abs(0.08838834764831842*alpha_diffusion[0]*rdvSq40); 
  alpha_mid += std::abs(0.08838834764831842*alpha_diffusion[0]*rdvSq41); 
  out[2] += 0.6123724356957944*(alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[3]*alpha_drag[11]+f[1]*alpha_drag[9]+f[0]*alpha_drag[8]); 
  out[4] += 0.6123724356957944*(alpha_drag[2]*f[4]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[5] += 0.6123724356957944*(f[5]*alpha_drag[11]+f[0]*alpha_drag[9]+f[1]*alpha_drag[8]); 
  out[6] += 0.6123724356957944*(f[6]*alpha_drag[11]+f[4]*alpha_drag[9]+f[2]*alpha_drag[8]+alpha_drag[2]*f[6]+alpha_drag[1]*f[5]+alpha_drag[0]*f[3]); 
  out[7] += 0.6123724356957944*(f[7]*alpha_drag[11]+f[2]*alpha_drag[9]+f[4]*alpha_drag[8]+alpha_drag[2]*f[7]+alpha_drag[0]*f[5]+alpha_drag[1]*f[3]); 
return alpha_mid; 

} 
double VmLBOconstNuVol1x2vSerP2(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdv20 = 2/dxv[1]; 
  const double rdvSq40 = 4/(dxv[1]*dxv[1]); 
  const double rdv21 = 2/dxv[2]; 
  const double rdvSq41 = 4/(dxv[2]*dxv[2]); 

  double alpha_mid = 0.0; 
  double alpha_drag[40]; 
  double alpha_diffusion[20]; 

  alpha_drag[0] = (2.0*nuU[0]-2.828427124746191*w[1]*nu)*rdv20; 
  alpha_drag[1] = 2.0*nuU[1]*rdv20; 
  alpha_drag[2] = -0.8164965809277261*dxv[1]*nu*rdv20; 
  alpha_drag[7] = 2.0*nuU[2]*rdv20; 
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[0]-0.1976423537605236*alpha_drag[7]); 
  alpha_drag[20] = (2.0*nuU[20]-2.828427124746191*w[2]*nu)*rdv21; 
  alpha_drag[21] = 2.0*nuU[21]*rdv21; 
  alpha_drag[23] = -0.8164965809277261*dxv[2]*nu*rdv21; 
  alpha_drag[27] = 2.0*nuU[22]*rdv21; 
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[20]-0.1976423537605236*alpha_drag[27]); 
  alpha_diffusion[0] = 2.0*nuVtSq[0]; 
  alpha_diffusion[1] = 2.0*nuVtSq[1]; 
  alpha_diffusion[7] = 2.0*nuVtSq[2]; 
  alpha_mid += std::abs((0.08838834764831842*alpha_diffusion[0]-0.09882117688026182*alpha_diffusion[7])*rdvSq40); 
  alpha_mid += std::abs((0.08838834764831842*alpha_diffusion[0]-0.09882117688026182*alpha_diffusion[7])*rdvSq41); 
  out[2] += 0.6123724356957944*(alpha_drag[7]*f[7]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*alpha_drag[27]+f[3]*alpha_drag[23]+f[1]*alpha_drag[21]+f[0]*alpha_drag[20]); 
  out[4] += 0.5477225575051661*(alpha_drag[1]*f[7]+f[1]*alpha_drag[7])+0.6123724356957944*(alpha_drag[2]*f[4]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[5] += 0.5477225575051661*f[1]*alpha_drag[27]+0.6123724356957944*f[5]*alpha_drag[23]+0.5477225575051661*f[7]*alpha_drag[21]+0.6123724356957944*(f[0]*alpha_drag[21]+f[1]*alpha_drag[20]); 
  out[6] += 0.6123724356957944*(f[11]*alpha_drag[27]+f[6]*alpha_drag[23]+f[4]*alpha_drag[21]+f[2]*alpha_drag[20]+alpha_drag[7]*f[13]+alpha_drag[2]*f[6]+alpha_drag[1]*f[5]+alpha_drag[0]*f[3]); 
  out[8] += 2.371708245126284*(alpha_diffusion[7]*f[7]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq40+1.369306393762915*alpha_drag[7]*f[11]+1.224744871391589*alpha_drag[2]*f[8]+1.369306393762915*(alpha_drag[1]*f[4]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[9] += 2.371708245126284*(alpha_diffusion[7]*f[7]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq41+1.369306393762915*f[13]*alpha_drag[27]+1.224744871391589*f[9]*alpha_drag[23]+1.369306393762915*(f[0]*alpha_drag[23]+f[5]*alpha_drag[21]+f[3]*alpha_drag[20]); 
  out[10] += 0.5477225575051661*f[4]*alpha_drag[27]+0.6123724356957944*f[10]*alpha_drag[23]+0.5477225575051661*f[11]*alpha_drag[21]+0.6123724356957944*(f[2]*alpha_drag[21]+f[4]*alpha_drag[20])+0.5477225575051661*alpha_drag[1]*f[13]+0.6123724356957944*alpha_drag[2]*f[10]+0.5477225575051661*f[5]*alpha_drag[7]+0.6123724356957944*(alpha_drag[0]*f[5]+alpha_drag[1]*f[3]); 
  out[11] += 0.6123724356957944*alpha_drag[2]*f[11]+0.3912303982179757*alpha_drag[7]*f[7]+0.6123724356957944*(alpha_drag[0]*f[7]+f[0]*alpha_drag[7])+0.5477225575051661*alpha_drag[1]*f[1]; 
  out[12] += (2.121320343559642*(alpha_diffusion[1]*f[7]+f[1]*alpha_diffusion[7])+2.371708245126284*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvSq40+1.224744871391589*(alpha_drag[2]*f[12]+alpha_drag[1]*f[11]+f[4]*alpha_drag[7])+1.369306393762915*(alpha_drag[0]*f[4]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 
  out[13] += 0.3912303982179757*f[7]*alpha_drag[27]+0.6123724356957944*(f[0]*alpha_drag[27]+f[13]*alpha_drag[23])+0.5477225575051661*f[1]*alpha_drag[21]+0.6123724356957944*f[7]*alpha_drag[20]; 
  out[14] += 2.371708245126284*(alpha_diffusion[7]*f[13]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[3])*rdvSq40+0.6123724356957944*(f[14]*alpha_drag[23]+f[12]*alpha_drag[21]+f[8]*alpha_drag[20])+1.369306393762915*alpha_drag[7]*f[17]+1.224744871391589*alpha_drag[2]*f[14]+1.369306393762915*(alpha_drag[1]*f[10]+alpha_drag[0]*f[6]+alpha_drag[2]*f[3]); 
  out[15] += (2.121320343559642*(alpha_diffusion[1]*f[7]+f[1]*alpha_diffusion[7])+2.371708245126284*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvSq41+1.224744871391589*f[5]*alpha_drag[27]+(1.224744871391589*f[15]+1.369306393762915*f[1])*alpha_drag[23]+1.224744871391589*f[13]*alpha_drag[21]+1.369306393762915*(f[3]*alpha_drag[21]+f[5]*alpha_drag[20]); 
  out[16] += 2.371708245126284*(alpha_diffusion[7]*f[11]+alpha_diffusion[1]*f[4]+alpha_diffusion[0]*f[2])*rdvSq41+1.369306393762915*f[17]*alpha_drag[27]+1.224744871391589*f[16]*alpha_drag[23]+1.369306393762915*(f[2]*alpha_drag[23]+f[10]*alpha_drag[21]+f[6]*alpha_drag[20])+0.6123724356957944*(alpha_drag[2]*f[16]+alpha_drag[1]*f[15]+alpha_drag[0]*f[9]); 
  out[17] += 0.3912303982179757*f[11]*alpha_drag[27]+0.6123724356957944*(f[2]*alpha_drag[27]+f[17]*alpha_drag[23])+0.5477225575051661*f[4]*alpha_drag[21]+0.6123724356957944*(f[11]*alpha_drag[20]+alpha_drag[2]*f[17])+0.3912303982179757*alpha_drag[7]*f[13]+0.6123724356957944*(alpha_drag[0]*f[13]+f[3]*alpha_drag[7])+0.5477225575051661*alpha_drag[1]*f[5]; 
  out[18] += (2.121320343559642*(alpha_diffusion[1]*f[13]+f[5]*alpha_diffusion[7])+2.371708245126284*(alpha_diffusion[0]*f[5]+alpha_diffusion[1]*f[3]))*rdvSq40+0.5477225575051661*f[12]*alpha_drag[27]+0.6123724356957944*(f[18]*alpha_drag[23]+f[8]*alpha_drag[21]+f[12]*alpha_drag[20])+1.224744871391589*(alpha_drag[2]*f[18]+alpha_drag[1]*f[17]+alpha_drag[7]*f[10])+1.369306393762915*(alpha_drag[0]*f[10]+alpha_drag[1]*f[6]+alpha_drag[2]*f[5]); 
  out[19] += (2.121320343559642*(alpha_diffusion[1]*f[11]+f[4]*alpha_diffusion[7])+2.371708245126284*(alpha_diffusion[0]*f[4]+alpha_diffusion[1]*f[2]))*rdvSq41+1.224744871391589*f[10]*alpha_drag[27]+(1.224744871391589*f[19]+1.369306393762915*f[4])*alpha_drag[23]+1.224744871391589*f[17]*alpha_drag[21]+1.369306393762915*(f[6]*alpha_drag[21]+f[10]*alpha_drag[20])+0.6123724356957944*alpha_drag[2]*f[19]+0.5477225575051661*alpha_drag[7]*f[15]+0.6123724356957944*(alpha_drag[0]*f[15]+alpha_drag[1]*f[9]); 
return alpha_mid; 

} 
double VmLBOconstNuVol1x2vSerP3(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdv20 = 2/dxv[1]; 
  const double rdvSq40 = 4/(dxv[1]*dxv[1]); 
  const double rdv21 = 2/dxv[2]; 
  const double rdvSq41 = 4/(dxv[2]*dxv[2]); 

  double alpha_mid = 0.0; 
  double alpha_drag[64]; 
  double alpha_diffusion[32]; 

  alpha_drag[0] = (2.0*nuU[0]-2.828427124746191*w[1]*nu)*rdv20; 
  alpha_drag[1] = 2.0*nuU[1]*rdv20; 
  alpha_drag[2] = -0.8164965809277261*dxv[1]*nu*rdv20; 
  alpha_drag[7] = 2.0*nuU[2]*rdv20; 
  alpha_drag[17] = 2.0*nuU[3]*rdv20; 
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[0]-0.1976423537605236*alpha_drag[7]); 
  alpha_drag[32] = (2.0*nuU[32]-2.828427124746191*w[2]*nu)*rdv21; 
  alpha_drag[33] = 2.0*nuU[33]*rdv21; 
  alpha_drag[35] = -0.8164965809277261*dxv[2]*nu*rdv21; 
  alpha_drag[39] = 2.0*nuU[34]*rdv21; 
  alpha_drag[49] = 2.0*nuU[35]*rdv21; 
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[32]-0.1976423537605236*alpha_drag[39]); 
  alpha_diffusion[0] = 2.0*nuVtSq[0]; 
  alpha_diffusion[1] = 2.0*nuVtSq[1]; 
  alpha_diffusion[7] = 2.0*nuVtSq[2]; 
  alpha_diffusion[17] = 2.0*nuVtSq[3]; 
  alpha_mid += std::abs((0.08838834764831842*alpha_diffusion[0]-0.09882117688026182*alpha_diffusion[7])*rdvSq40); 
  alpha_mid += std::abs((0.08838834764831842*alpha_diffusion[0]-0.09882117688026182*alpha_diffusion[7])*rdvSq41); 
  out[2] += 0.6123724356957944*(alpha_drag[17]*f[17]+alpha_drag[7]*f[7]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[17]*alpha_drag[49]+f[7]*alpha_drag[39]+f[3]*alpha_drag[35]+f[1]*alpha_drag[33]+f[0]*alpha_drag[32]); 
  out[4] += 0.537852874200477*(alpha_drag[7]*f[17]+f[7]*alpha_drag[17])+0.5477225575051661*(alpha_drag[1]*f[7]+f[1]*alpha_drag[7])+0.6123724356957944*(alpha_drag[2]*f[4]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[5] += 0.537852874200477*f[7]*alpha_drag[49]+(0.537852874200477*f[17]+0.5477225575051661*f[1])*alpha_drag[39]+0.6123724356957944*f[5]*alpha_drag[35]+0.5477225575051661*f[7]*alpha_drag[33]+0.6123724356957944*(f[0]*alpha_drag[33]+f[1]*alpha_drag[32]); 
  out[6] += 0.6123724356957944*(f[23]*alpha_drag[49]+f[11]*alpha_drag[39]+f[6]*alpha_drag[35]+f[4]*alpha_drag[33]+f[2]*alpha_drag[32]+alpha_drag[17]*f[25]+alpha_drag[7]*f[13]+alpha_drag[2]*f[6]+alpha_drag[1]*f[5]+alpha_drag[0]*f[3]); 
  out[8] += 2.371708245126284*(alpha_diffusion[17]*f[17]+alpha_diffusion[7]*f[7]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq40+1.369306393762915*(alpha_drag[17]*f[23]+alpha_drag[7]*f[11])+1.224744871391589*alpha_drag[2]*f[8]+1.369306393762915*(alpha_drag[1]*f[4]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[9] += 2.371708245126284*(alpha_diffusion[17]*f[17]+alpha_diffusion[7]*f[7]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq41+1.369306393762915*(f[25]*alpha_drag[49]+f[13]*alpha_drag[39])+1.224744871391589*f[9]*alpha_drag[35]+1.369306393762915*(f[0]*alpha_drag[35]+f[5]*alpha_drag[33]+f[3]*alpha_drag[32]); 
  out[10] += 0.537852874200477*f[11]*alpha_drag[49]+(0.537852874200477*f[23]+0.5477225575051661*f[4])*alpha_drag[39]+0.6123724356957944*f[10]*alpha_drag[35]+0.5477225575051661*f[11]*alpha_drag[33]+0.6123724356957944*(f[2]*alpha_drag[33]+f[4]*alpha_drag[32])+0.537852874200477*alpha_drag[7]*f[25]+f[13]*(0.537852874200477*alpha_drag[17]+0.5477225575051661*alpha_drag[1])+0.6123724356957944*alpha_drag[2]*f[10]+0.5477225575051661*f[5]*alpha_drag[7]+0.6123724356957944*(alpha_drag[0]*f[5]+alpha_drag[1]*f[3]); 
  out[11] += 0.3651483716701108*alpha_drag[17]*f[17]+0.537852874200477*(alpha_drag[1]*f[17]+f[1]*alpha_drag[17])+0.6123724356957944*alpha_drag[2]*f[11]+0.3912303982179757*alpha_drag[7]*f[7]+0.6123724356957944*(alpha_drag[0]*f[7]+f[0]*alpha_drag[7])+0.5477225575051661*alpha_drag[1]*f[1]; 
  out[12] += (2.08309522448824*(alpha_diffusion[7]*f[17]+f[7]*alpha_diffusion[17])+2.121320343559642*(alpha_diffusion[1]*f[7]+f[1]*alpha_diffusion[7])+2.371708245126284*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvSq40+1.202675588605909*(alpha_drag[7]*f[23]+f[11]*alpha_drag[17])+1.224744871391589*(alpha_drag[2]*f[12]+alpha_drag[1]*f[11]+f[4]*alpha_drag[7])+1.369306393762915*(alpha_drag[0]*f[4]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 
  out[13] += (0.3651483716701108*f[17]+0.537852874200477*f[1])*alpha_drag[49]+0.3912303982179757*f[7]*alpha_drag[39]+0.6123724356957944*(f[0]*alpha_drag[39]+f[13]*alpha_drag[35])+(0.537852874200477*f[17]+0.5477225575051661*f[1])*alpha_drag[33]+0.6123724356957944*f[7]*alpha_drag[32]; 
  out[14] += 2.371708245126284*(alpha_diffusion[17]*f[25]+alpha_diffusion[7]*f[13]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[3])*rdvSq40+0.6123724356957944*(f[14]*alpha_drag[35]+f[12]*alpha_drag[33]+f[8]*alpha_drag[32])+1.369306393762915*(alpha_drag[17]*f[29]+alpha_drag[7]*f[20])+1.224744871391589*alpha_drag[2]*f[14]+1.369306393762915*(alpha_drag[1]*f[10]+alpha_drag[0]*f[6]+alpha_drag[2]*f[3]); 
  out[15] += (2.08309522448824*(alpha_diffusion[7]*f[17]+f[7]*alpha_diffusion[17])+2.121320343559642*(alpha_diffusion[1]*f[7]+f[1]*alpha_diffusion[7])+2.371708245126284*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvSq41+1.202675588605909*f[13]*alpha_drag[49]+(1.202675588605909*f[25]+1.224744871391589*f[5])*alpha_drag[39]+(1.224744871391589*f[15]+1.369306393762915*f[1])*alpha_drag[35]+1.224744871391589*f[13]*alpha_drag[33]+1.369306393762915*(f[3]*alpha_drag[33]+f[5]*alpha_drag[32]); 
  out[16] += 2.371708245126284*(alpha_diffusion[17]*f[23]+alpha_diffusion[7]*f[11]+alpha_diffusion[1]*f[4]+alpha_diffusion[0]*f[2])*rdvSq41+1.369306393762915*(f[29]*alpha_drag[49]+f[20]*alpha_drag[39])+1.224744871391589*f[16]*alpha_drag[35]+1.369306393762915*(f[2]*alpha_drag[35]+f[10]*alpha_drag[33]+f[6]*alpha_drag[32])+0.6123724356957944*(alpha_drag[2]*f[16]+alpha_drag[1]*f[15]+alpha_drag[0]*f[9]); 
  out[18] += 8.100925873009823*(alpha_diffusion[17]*f[23]+alpha_diffusion[7]*f[11]+alpha_diffusion[1]*f[4]+alpha_diffusion[0]*f[2])*rdvSq40+1.837117307087383*alpha_drag[2]*f[18]+0.9354143466934851*alpha_drag[17]*f[17]+2.091650066335188*(alpha_drag[1]*f[12]+alpha_drag[0]*f[8])+0.9354143466934851*alpha_drag[7]*f[7]+2.806243040080455*alpha_drag[2]*f[2]+0.9354143466934851*(alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[19] += 8.100925873009823*(alpha_diffusion[17]*f[25]+alpha_diffusion[7]*f[13]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[3])*rdvSq41+0.9354143466934851*(f[17]*alpha_drag[49]+f[7]*alpha_drag[39])+(1.837117307087383*f[19]+2.806243040080455*f[3])*alpha_drag[35]+(2.091650066335188*f[15]+0.9354143466934851*f[1])*alpha_drag[33]+(2.091650066335188*f[9]+0.9354143466934851*f[0])*alpha_drag[32]; 
  out[20] += (0.3651483716701108*f[23]+0.537852874200477*f[4])*alpha_drag[49]+0.3912303982179757*f[11]*alpha_drag[39]+0.6123724356957944*(f[2]*alpha_drag[39]+f[20]*alpha_drag[35])+(0.537852874200477*f[23]+0.5477225575051661*f[4])*alpha_drag[33]+0.6123724356957944*f[11]*alpha_drag[32]+(0.3651483716701108*alpha_drag[17]+0.537852874200477*alpha_drag[1])*f[25]+0.6123724356957944*alpha_drag[2]*f[20]+0.537852874200477*f[5]*alpha_drag[17]+0.3912303982179757*alpha_drag[7]*f[13]+0.6123724356957944*(alpha_drag[0]*f[13]+f[3]*alpha_drag[7])+0.5477225575051661*alpha_drag[1]*f[5]; 
  out[21] += (2.08309522448824*(alpha_diffusion[7]*f[25]+f[13]*alpha_diffusion[17])+2.121320343559642*(alpha_diffusion[1]*f[13]+f[5]*alpha_diffusion[7])+2.371708245126284*(alpha_diffusion[0]*f[5]+alpha_diffusion[1]*f[3]))*rdvSq40+0.5477225575051661*f[12]*alpha_drag[39]+0.6123724356957944*(f[21]*alpha_drag[35]+f[8]*alpha_drag[33]+f[12]*alpha_drag[32])+1.202675588605909*alpha_drag[7]*f[29]+1.224744871391589*alpha_drag[2]*f[21]+1.202675588605909*alpha_drag[17]*f[20]+1.224744871391589*(alpha_drag[1]*f[20]+alpha_drag[7]*f[10])+1.369306393762915*(alpha_drag[0]*f[10]+alpha_drag[1]*f[6]+alpha_drag[2]*f[5]); 
  out[22] += (2.08309522448824*(alpha_diffusion[7]*f[23]+f[11]*alpha_diffusion[17])+2.121320343559642*(alpha_diffusion[1]*f[11]+f[4]*alpha_diffusion[7])+2.371708245126284*(alpha_diffusion[0]*f[4]+alpha_diffusion[1]*f[2]))*rdvSq41+1.202675588605909*f[20]*alpha_drag[49]+(1.202675588605909*f[29]+1.224744871391589*f[10])*alpha_drag[39]+(1.224744871391589*f[22]+1.369306393762915*f[4])*alpha_drag[35]+1.224744871391589*f[20]*alpha_drag[33]+1.369306393762915*(f[6]*alpha_drag[33]+f[10]*alpha_drag[32])+0.6123724356957944*alpha_drag[2]*f[22]+0.5477225575051661*alpha_drag[7]*f[15]+0.6123724356957944*(alpha_drag[0]*f[15]+alpha_drag[1]*f[9]); 
  out[23] += 0.6123724356957944*alpha_drag[2]*f[23]+(0.3651483716701108*alpha_drag[7]+0.6123724356957944*alpha_drag[0])*f[17]+(0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha_drag[17]+0.537852874200477*(alpha_drag[1]*f[7]+f[1]*alpha_drag[7]); 
  out[24] += (7.115124735378852*(alpha_diffusion[7]*f[23]+f[11]*alpha_diffusion[17])+7.24568837309472*(alpha_diffusion[1]*f[11]+f[4]*alpha_diffusion[7])+8.100925873009823*(alpha_diffusion[0]*f[4]+alpha_diffusion[1]*f[2]))*rdvSq40+1.837117307087383*alpha_drag[2]*f[24]+0.8215838362577489*(alpha_drag[7]*f[17]+f[7]*alpha_drag[17])+1.870828693386971*alpha_drag[7]*f[12]+2.091650066335188*(alpha_drag[0]*f[12]+alpha_drag[1]*f[8])+0.8366600265340755*(alpha_drag[1]*f[7]+f[1]*alpha_drag[7])+2.806243040080455*alpha_drag[2]*f[4]+0.9354143466934851*(alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[25] += (0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha_drag[49]+(0.3651483716701108*f[17]+0.537852874200477*f[1])*alpha_drag[39]+0.6123724356957944*f[25]*alpha_drag[35]+0.537852874200477*f[7]*alpha_drag[33]+0.6123724356957944*f[17]*alpha_drag[32]; 
  out[26] += 8.100925873009823*(alpha_diffusion[17]*f[29]+alpha_diffusion[7]*f[20]+alpha_diffusion[1]*f[10]+alpha_diffusion[0]*f[6])*rdvSq40+0.6123724356957944*(f[26]*alpha_drag[35]+f[24]*alpha_drag[33]+f[18]*alpha_drag[32])+1.837117307087383*alpha_drag[2]*f[26]+0.9354143466934851*alpha_drag[17]*f[25]+2.091650066335188*(alpha_drag[1]*f[21]+alpha_drag[0]*f[14])+0.9354143466934851*alpha_drag[7]*f[13]+2.806243040080455*alpha_drag[2]*f[6]+0.9354143466934851*(alpha_drag[1]*f[5]+alpha_drag[0]*f[3]); 
  out[27] += (7.115124735378852*(alpha_diffusion[7]*f[25]+f[13]*alpha_diffusion[17])+7.24568837309472*(alpha_diffusion[1]*f[13]+f[5]*alpha_diffusion[7])+8.100925873009823*(alpha_diffusion[0]*f[5]+alpha_diffusion[1]*f[3]))*rdvSq41+0.8215838362577489*f[7]*alpha_drag[49]+(0.8215838362577489*f[17]+1.870828693386971*f[15]+0.8366600265340755*f[1])*alpha_drag[39]+(1.837117307087383*f[27]+2.806243040080455*f[5])*alpha_drag[35]+(2.091650066335188*f[9]+0.8366600265340755*f[7]+0.9354143466934851*f[0])*alpha_drag[33]+(2.091650066335188*f[15]+0.9354143466934851*f[1])*alpha_drag[32]; 
  out[28] += 8.100925873009823*(alpha_diffusion[17]*f[29]+alpha_diffusion[7]*f[20]+alpha_diffusion[1]*f[10]+alpha_diffusion[0]*f[6])*rdvSq41+0.9354143466934851*(f[23]*alpha_drag[49]+f[11]*alpha_drag[39])+(1.837117307087383*f[28]+2.806243040080455*f[6])*alpha_drag[35]+(2.091650066335188*f[22]+0.9354143466934851*f[4])*alpha_drag[33]+(2.091650066335188*f[16]+0.9354143466934851*f[2])*alpha_drag[32]+0.6123724356957944*(alpha_drag[2]*f[28]+alpha_drag[1]*f[27]+alpha_drag[0]*f[19]); 
  out[29] += (0.3651483716701108*f[11]+0.6123724356957944*f[2])*alpha_drag[49]+(0.3651483716701108*f[23]+0.537852874200477*f[4])*alpha_drag[39]+0.6123724356957944*f[29]*alpha_drag[35]+0.537852874200477*f[11]*alpha_drag[33]+0.6123724356957944*(f[23]*alpha_drag[32]+alpha_drag[2]*f[29])+(0.3651483716701108*alpha_drag[7]+0.6123724356957944*alpha_drag[0])*f[25]+(0.3651483716701108*f[13]+0.6123724356957944*f[3])*alpha_drag[17]+0.537852874200477*(alpha_drag[1]*f[13]+f[5]*alpha_drag[7]); 
  out[30] += (7.115124735378852*alpha_diffusion[7]*f[29]+(7.115124735378852*alpha_diffusion[17]+7.24568837309472*alpha_diffusion[1])*f[20]+8.100925873009823*(alpha_diffusion[0]*f[10]+alpha_diffusion[1]*f[6])+7.24568837309472*alpha_diffusion[7]*f[10])*rdvSq40+0.5477225575051661*f[24]*alpha_drag[39]+0.6123724356957944*(f[30]*alpha_drag[35]+f[18]*alpha_drag[33]+f[24]*alpha_drag[32])+1.837117307087383*alpha_drag[2]*f[30]+0.8215838362577489*alpha_drag[7]*f[25]+(1.870828693386971*alpha_drag[7]+2.091650066335188*alpha_drag[0])*f[21]+0.8215838362577489*f[13]*alpha_drag[17]+alpha_drag[1]*(2.091650066335188*f[14]+0.8366600265340755*f[13])+2.806243040080455*alpha_drag[2]*f[10]+0.8366600265340755*f[5]*alpha_drag[7]+0.9354143466934851*(alpha_drag[0]*f[5]+alpha_drag[1]*f[3]); 
  out[31] += (7.115124735378852*alpha_diffusion[7]*f[29]+(7.115124735378852*alpha_diffusion[17]+7.24568837309472*alpha_diffusion[1])*f[20]+8.100925873009823*(alpha_diffusion[0]*f[10]+alpha_diffusion[1]*f[6])+7.24568837309472*alpha_diffusion[7]*f[10])*rdvSq41+0.8215838362577489*f[11]*alpha_drag[49]+(0.8215838362577489*f[23]+1.870828693386971*f[22]+0.8366600265340755*f[4])*alpha_drag[39]+(1.837117307087383*f[31]+2.806243040080455*f[10])*alpha_drag[35]+(2.091650066335188*f[16]+0.8366600265340755*f[11]+0.9354143466934851*f[2])*alpha_drag[33]+(2.091650066335188*f[22]+0.9354143466934851*f[4])*alpha_drag[32]+0.6123724356957944*alpha_drag[2]*f[31]+0.5477225575051661*alpha_drag[7]*f[27]+0.6123724356957944*(alpha_drag[0]*f[27]+alpha_drag[1]*f[19]); 
return alpha_mid; 

} 
