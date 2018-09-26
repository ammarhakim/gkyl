#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x2vMaxP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdvx2 = 2/dxv[1]; 
  const double rdvxSq4 = 4/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2/dxv[2]; 
  const double rdvySq4 = 4/(dxv[2]*dxv[2]); 

  double alpha_mid = 0.0; 
  double alpha_drag[8]; 
  double alpha_diffusion[4]; 

  // Expand rdv2*nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.0*nuU[0]-2.828427124746191*w[1]*nu)*rdvx2; 
  alpha_drag[1] = 2.0*nuU[1]*rdvx2; 
  alpha_drag[2] = -0.8164965809277261*dxv[1]*nu*rdvx2; 

  // cvard[1]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[0]); 

  // Expand rdv2*nu*(vy-uy) in phase basis.
  alpha_drag[4] = (2.0*nuU[2]-2.828427124746191*w[2]*nu)*rdvy2; 
  alpha_drag[5] = 2.0*nuU[3]*rdvy2; 
  alpha_drag[7] = -0.8164965809277261*dxv[2]*nu*rdvy2; 

  // cvard[2]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[4]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.0*nuVtSq[0]; 
  alpha_diffusion[1] = 2.0*nuVtSq[1]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.4714045207910317*alpha_diffusion[0]*rdvxSq4); 
  alpha_mid += std::abs(0.4714045207910317*alpha_diffusion[0]*rdvySq4); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[3]*alpha_drag[7]+f[1]*alpha_drag[5]+f[0]*alpha_drag[4]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol1x2vMaxP2(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdvx2 = 2/dxv[1]; 
  const double rdvxSq4 = 4/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2/dxv[2]; 
  const double rdvySq4 = 4/(dxv[2]*dxv[2]); 

  double alpha_mid = 0.0; 
  double alpha_drag[20]; 
  double alpha_diffusion[10]; 

  // Expand rdv2*nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.0*nuU[0]-2.828427124746191*w[1]*nu)*rdvx2; 
  alpha_drag[1] = 2.0*nuU[1]*rdvx2; 
  alpha_drag[2] = -0.8164965809277261*dxv[1]*nu*rdvx2; 
  alpha_drag[7] = 2.0*nuU[2]*rdvx2; 

  // cvard[1]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[0]-0.1976423537605236*alpha_drag[7]); 

  // Expand rdv2*nu*(vy-uy) in phase basis.
  alpha_drag[10] = (2.0*nuU[3]-2.828427124746191*w[2]*nu)*rdvy2; 
  alpha_drag[11] = 2.0*nuU[4]*rdvy2; 
  alpha_drag[13] = -0.8164965809277261*dxv[2]*nu*rdvy2; 
  alpha_drag[17] = 2.0*nuU[5]*rdvy2; 

  // cvard[2]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[10]-0.1976423537605236*alpha_drag[17]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.0*nuVtSq[0]; 
  alpha_diffusion[1] = 2.0*nuVtSq[1]; 
  alpha_diffusion[7] = 2.0*nuVtSq[2]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((0.6363961030678926*alpha_diffusion[0]-0.711512473537885*alpha_diffusion[7])*rdvxSq4); 
  alpha_mid += std::abs((0.6363961030678926*alpha_diffusion[0]-0.711512473537885*alpha_diffusion[7])*rdvySq4); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alpha_drag[7]*f[7]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*alpha_drag[17]+f[3]*alpha_drag[13]+f[1]*alpha_drag[11]+f[0]*alpha_drag[10]); 
  out[4] += 0.5477225575051661*(alpha_drag[1]*f[7]+f[1]*alpha_drag[7])+0.6123724356957944*(alpha_drag[2]*f[4]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[5] += 0.5477225575051661*f[1]*alpha_drag[17]+0.6123724356957944*f[5]*alpha_drag[13]+0.5477225575051661*f[7]*alpha_drag[11]+0.6123724356957944*(f[0]*alpha_drag[11]+f[1]*alpha_drag[10]); 
  out[6] += 0.6123724356957944*(f[6]*alpha_drag[13]+f[4]*alpha_drag[11]+f[2]*alpha_drag[10]+alpha_drag[2]*f[6]+alpha_drag[1]*f[5]+alpha_drag[0]*f[3]); 
  out[8] += 2.371708245126284*(alpha_diffusion[7]*f[7]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4+1.224744871391589*alpha_drag[2]*f[8]+1.369306393762915*(alpha_drag[1]*f[4]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[9] += 2.371708245126284*(alpha_diffusion[7]*f[7]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvySq4+1.224744871391589*f[9]*alpha_drag[13]+1.369306393762915*(f[0]*alpha_drag[13]+f[5]*alpha_drag[11]+f[3]*alpha_drag[10]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol1x2vMaxP3(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdvx2 = 2/dxv[1]; 
  const double rdvxSq4 = 4/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2/dxv[2]; 
  const double rdvySq4 = 4/(dxv[2]*dxv[2]); 

  double alpha_mid = 0.0; 
  double alpha_drag[40]; 
  double alpha_diffusion[20]; 

  // Expand rdv2*nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.0*nuU[0]-2.828427124746191*w[1]*nu)*rdvx2; 
  alpha_drag[1] = 2.0*nuU[1]*rdvx2; 
  alpha_drag[2] = -0.8164965809277261*dxv[1]*nu*rdvx2; 
  alpha_drag[7] = 2.0*nuU[2]*rdvx2; 
  alpha_drag[17] = 2.0*nuU[3]*rdvx2; 

  // cvard[1]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[0]-0.1976423537605236*alpha_drag[7]); 

  // Expand rdv2*nu*(vy-uy) in phase basis.
  alpha_drag[20] = (2.0*nuU[4]-2.828427124746191*w[2]*nu)*rdvy2; 
  alpha_drag[21] = 2.0*nuU[5]*rdvy2; 
  alpha_drag[23] = -0.8164965809277261*dxv[2]*nu*rdvy2; 
  alpha_drag[27] = 2.0*nuU[6]*rdvy2; 
  alpha_drag[37] = 2.0*nuU[7]*rdvy2; 

  // cvard[2]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.1767766952966368*alpha_drag[20]-0.1976423537605236*alpha_drag[27]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.0*nuVtSq[0]; 
  alpha_diffusion[1] = 2.0*nuVtSq[1]; 
  alpha_diffusion[7] = 2.0*nuVtSq[2]; 
  alpha_diffusion[17] = 2.0*nuVtSq[3]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((0.8081220356417689*alpha_diffusion[0]-0.9035079029052515*alpha_diffusion[7])*rdvxSq4); 
  alpha_mid += std::abs((0.8081220356417689*alpha_diffusion[0]-0.9035079029052515*alpha_diffusion[7])*rdvySq4); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alpha_drag[17]*f[17]+alpha_drag[7]*f[7]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[17]*alpha_drag[37]+f[7]*alpha_drag[27]+f[3]*alpha_drag[23]+f[1]*alpha_drag[21]+f[0]*alpha_drag[20]); 
  out[4] += 0.537852874200477*(alpha_drag[7]*f[17]+f[7]*alpha_drag[17])+0.5477225575051661*(alpha_drag[1]*f[7]+f[1]*alpha_drag[7])+0.6123724356957944*(alpha_drag[2]*f[4]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[5] += 0.537852874200477*f[7]*alpha_drag[37]+(0.537852874200477*f[17]+0.5477225575051661*f[1])*alpha_drag[27]+0.6123724356957944*f[5]*alpha_drag[23]+0.5477225575051661*f[7]*alpha_drag[21]+0.6123724356957944*(f[0]*alpha_drag[21]+f[1]*alpha_drag[20]); 
  out[6] += 0.6123724356957944*(f[11]*alpha_drag[27]+f[6]*alpha_drag[23]+f[4]*alpha_drag[21]+f[2]*alpha_drag[20]+alpha_drag[7]*f[13]+alpha_drag[2]*f[6]+alpha_drag[1]*f[5]+alpha_drag[0]*f[3]); 
  out[8] += 2.371708245126284*(alpha_diffusion[17]*f[17]+alpha_diffusion[7]*f[7]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4+1.369306393762915*alpha_drag[7]*f[11]+1.224744871391589*alpha_drag[2]*f[8]+1.369306393762915*(alpha_drag[1]*f[4]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[9] += 2.371708245126284*(alpha_diffusion[17]*f[17]+alpha_diffusion[7]*f[7]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvySq4+1.369306393762915*f[13]*alpha_drag[27]+1.224744871391589*f[9]*alpha_drag[23]+1.369306393762915*(f[0]*alpha_drag[23]+f[5]*alpha_drag[21]+f[3]*alpha_drag[20]); 
  out[10] += 0.537852874200477*f[11]*alpha_drag[37]+0.5477225575051661*f[4]*alpha_drag[27]+0.6123724356957944*f[10]*alpha_drag[23]+0.5477225575051661*f[11]*alpha_drag[21]+0.6123724356957944*(f[2]*alpha_drag[21]+f[4]*alpha_drag[20])+f[13]*(0.537852874200477*alpha_drag[17]+0.5477225575051661*alpha_drag[1])+0.6123724356957944*alpha_drag[2]*f[10]+0.5477225575051661*f[5]*alpha_drag[7]+0.6123724356957944*(alpha_drag[0]*f[5]+alpha_drag[1]*f[3]); 
  out[11] += 0.3651483716701108*alpha_drag[17]*f[17]+0.537852874200477*(alpha_drag[1]*f[17]+f[1]*alpha_drag[17])+0.6123724356957944*alpha_drag[2]*f[11]+0.3912303982179757*alpha_drag[7]*f[7]+0.6123724356957944*(alpha_drag[0]*f[7]+f[0]*alpha_drag[7])+0.5477225575051661*alpha_drag[1]*f[1]; 
  out[12] += (2.08309522448824*(alpha_diffusion[7]*f[17]+f[7]*alpha_diffusion[17])+2.121320343559642*(alpha_diffusion[1]*f[7]+f[1]*alpha_diffusion[7])+2.371708245126284*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvxSq4+1.202675588605909*f[11]*alpha_drag[17]+1.224744871391589*(alpha_drag[2]*f[12]+alpha_drag[1]*f[11]+f[4]*alpha_drag[7])+1.369306393762915*(alpha_drag[0]*f[4]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 
  out[13] += (0.3651483716701108*f[17]+0.537852874200477*f[1])*alpha_drag[37]+0.3912303982179757*f[7]*alpha_drag[27]+0.6123724356957944*(f[0]*alpha_drag[27]+f[13]*alpha_drag[23])+(0.537852874200477*f[17]+0.5477225575051661*f[1])*alpha_drag[21]+0.6123724356957944*f[7]*alpha_drag[20]; 
  out[14] += 2.371708245126284*(alpha_diffusion[7]*f[13]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[3])*rdvxSq4+0.6123724356957944*(f[14]*alpha_drag[23]+f[12]*alpha_drag[21]+f[8]*alpha_drag[20])+1.224744871391589*alpha_drag[2]*f[14]+1.369306393762915*(alpha_drag[1]*f[10]+alpha_drag[0]*f[6]+alpha_drag[2]*f[3]); 
  out[15] += (2.08309522448824*(alpha_diffusion[7]*f[17]+f[7]*alpha_diffusion[17])+2.121320343559642*(alpha_diffusion[1]*f[7]+f[1]*alpha_diffusion[7])+2.371708245126284*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvySq4+1.202675588605909*f[13]*alpha_drag[37]+1.224744871391589*f[5]*alpha_drag[27]+(1.224744871391589*f[15]+1.369306393762915*f[1])*alpha_drag[23]+1.224744871391589*f[13]*alpha_drag[21]+1.369306393762915*(f[3]*alpha_drag[21]+f[5]*alpha_drag[20]); 
  out[16] += 2.371708245126284*(alpha_diffusion[7]*f[11]+alpha_diffusion[1]*f[4]+alpha_diffusion[0]*f[2])*rdvySq4+1.224744871391589*f[16]*alpha_drag[23]+1.369306393762915*(f[2]*alpha_drag[23]+f[10]*alpha_drag[21]+f[6]*alpha_drag[20])+0.6123724356957944*(alpha_drag[2]*f[16]+alpha_drag[1]*f[15]+alpha_drag[0]*f[9]); 
  out[18] += 8.100925873009823*(alpha_diffusion[7]*f[11]+alpha_diffusion[1]*f[4]+alpha_diffusion[0]*f[2])*rdvxSq4+1.837117307087383*alpha_drag[2]*f[18]+0.9354143466934851*alpha_drag[17]*f[17]+2.091650066335188*(alpha_drag[1]*f[12]+alpha_drag[0]*f[8])+0.9354143466934851*alpha_drag[7]*f[7]+2.806243040080455*alpha_drag[2]*f[2]+0.9354143466934851*(alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[19] += 8.100925873009823*(alpha_diffusion[7]*f[13]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[3])*rdvySq4+0.9354143466934851*(f[17]*alpha_drag[37]+f[7]*alpha_drag[27])+(1.837117307087383*f[19]+2.806243040080455*f[3])*alpha_drag[23]+(2.091650066335188*f[15]+0.9354143466934851*f[1])*alpha_drag[21]+(2.091650066335188*f[9]+0.9354143466934851*f[0])*alpha_drag[20]; 

  return alpha_mid; 

} 
