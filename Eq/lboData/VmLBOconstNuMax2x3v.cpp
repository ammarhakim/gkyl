#include <VmLBOModDecl.h> 
double VmLBOconstNuVol2x3vMaxP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdvx2 = 2/dxv[2]; 
  const double rdvxSq4 = 4/(dxv[2]*dxv[2]); 
  const double rdvy2 = 2/dxv[3]; 
  const double rdvySq4 = 4/(dxv[3]*dxv[3]); 
  const double rdvz2 = 2/dxv[4]; 
  const double rdvzSq4 = 4/(dxv[4]*dxv[4]); 

  double alpha_mid = 0.0; 
  double alpha_drag[18]; 
  double alpha_diffusion[6]; 

  // Expand rdv2*nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.828427124746191*nuU[0]-5.656854249492382*w[2]*nu)*rdvx2; 
  alpha_drag[1] = 2.828427124746191*nuU[1]*rdvx2; 
  alpha_drag[2] = 2.828427124746191*nuU[2]*rdvx2; 
  alpha_drag[3] = -1.632993161855453*dxv[2]*nu*rdvx2; 

  // cvard[1]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[0]); 

  // Expand rdv2*nu*(vy-uy) in phase basis.
  alpha_drag[6] = (2.828427124746191*nuU[3]-5.656854249492382*w[3]*nu)*rdvy2; 
  alpha_drag[7] = 2.828427124746191*nuU[4]*rdvy2; 
  alpha_drag[8] = 2.828427124746191*nuU[5]*rdvy2; 
  alpha_drag[10] = -1.632993161855453*dxv[3]*nu*rdvy2; 

  // cvard[2]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[6]); 

  // Expand rdv2*nu*(vz-uz) in phase basis.
  alpha_drag[12] = (2.828427124746191*nuU[6]-5.656854249492382*w[4]*nu)*rdvz2; 
  alpha_drag[13] = 2.828427124746191*nuU[7]*rdvz2; 
  alpha_drag[14] = 2.828427124746191*nuU[8]*rdvz2; 
  alpha_drag[17] = -1.632993161855453*dxv[4]*nu*rdvz2; 

  // cvard[3]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[12]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.828427124746191*nuVtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*nuVtSq[1]; 
  alpha_diffusion[2] = 2.828427124746191*nuVtSq[2]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.2357022603955158*alpha_diffusion[0]*rdvxSq4); 
  alpha_mid += std::abs(0.2357022603955158*alpha_diffusion[0]*rdvySq4); 
  alpha_mid += std::abs(0.2357022603955158*alpha_diffusion[0]*rdvzSq4); 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.3061862178478971*(alpha_drag[3]*f[3]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[4]*alpha_drag[10]+f[2]*alpha_drag[8]+f[1]*alpha_drag[7]+f[0]*alpha_drag[6]); 
  out[5] += 0.3061862178478971*(f[5]*alpha_drag[17]+f[2]*alpha_drag[14]+f[1]*alpha_drag[13]+f[0]*alpha_drag[12]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol2x3vMaxP2(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdvx2 = 2/dxv[2]; 
  const double rdvxSq4 = 4/(dxv[2]*dxv[2]); 
  const double rdvy2 = 2/dxv[3]; 
  const double rdvySq4 = 4/(dxv[3]*dxv[3]); 
  const double rdvz2 = 2/dxv[4]; 
  const double rdvzSq4 = 4/(dxv[4]*dxv[4]); 

  double alpha_mid = 0.0; 
  double alpha_drag[63]; 
  double alpha_diffusion[21]; 

  // Expand rdv2*nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.828427124746191*nuU[0]-5.656854249492382*w[2]*nu)*rdvx2; 
  alpha_drag[1] = 2.828427124746191*nuU[1]*rdvx2; 
  alpha_drag[2] = 2.828427124746191*nuU[2]*rdvx2; 
  alpha_drag[3] = -1.632993161855453*dxv[2]*nu*rdvx2; 
  alpha_drag[6] = 2.828427124746191*nuU[3]*rdvx2; 
  alpha_drag[16] = 2.828427124746191*nuU[4]*rdvx2; 
  alpha_drag[17] = 2.828427124746191*nuU[5]*rdvx2; 

  // cvard[1]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[0]-0.09882117688026182*(alpha_drag[17]+alpha_drag[16])); 

  // Expand rdv2*nu*(vy-uy) in phase basis.
  alpha_drag[21] = (2.828427124746191*nuU[6]-5.656854249492382*w[3]*nu)*rdvy2; 
  alpha_drag[22] = 2.828427124746191*nuU[7]*rdvy2; 
  alpha_drag[23] = 2.828427124746191*nuU[8]*rdvy2; 
  alpha_drag[25] = -1.632993161855453*dxv[3]*nu*rdvy2; 
  alpha_drag[27] = 2.828427124746191*nuU[9]*rdvy2; 
  alpha_drag[37] = 2.828427124746191*nuU[10]*rdvy2; 
  alpha_drag[38] = 2.828427124746191*nuU[11]*rdvy2; 

  // cvard[2]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[21]-0.09882117688026182*(alpha_drag[38]+alpha_drag[37])); 

  // Expand rdv2*nu*(vz-uz) in phase basis.
  alpha_drag[42] = (2.828427124746191*nuU[12]-5.656854249492382*w[4]*nu)*rdvz2; 
  alpha_drag[43] = 2.828427124746191*nuU[13]*rdvz2; 
  alpha_drag[44] = 2.828427124746191*nuU[14]*rdvz2; 
  alpha_drag[47] = -1.632993161855453*dxv[4]*nu*rdvz2; 
  alpha_drag[48] = 2.828427124746191*nuU[15]*rdvz2; 
  alpha_drag[58] = 2.828427124746191*nuU[16]*rdvz2; 
  alpha_drag[59] = 2.828427124746191*nuU[17]*rdvz2; 

  // cvard[3]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[42]-0.09882117688026182*(alpha_drag[59]+alpha_drag[58])); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.828427124746191*nuVtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*nuVtSq[1]; 
  alpha_diffusion[2] = 2.828427124746191*nuVtSq[2]; 
  alpha_diffusion[6] = 2.828427124746191*nuVtSq[3]; 
  alpha_diffusion[16] = 2.828427124746191*nuVtSq[4]; 
  alpha_diffusion[17] = 2.828427124746191*nuVtSq[5]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((0.3181980515339463*alpha_diffusion[0]-0.3557562367689425*(alpha_diffusion[17]+alpha_diffusion[16]))*rdvxSq4); 
  alpha_mid += std::abs((0.3181980515339463*alpha_diffusion[0]-0.3557562367689425*(alpha_diffusion[17]+alpha_diffusion[16]))*rdvySq4); 
  alpha_mid += std::abs((0.3181980515339463*alpha_diffusion[0]-0.3557562367689425*(alpha_diffusion[17]+alpha_diffusion[16]))*rdvzSq4); 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.3061862178478971*(alpha_drag[17]*f[17]+alpha_drag[16]*f[16]+alpha_drag[6]*f[6]+alpha_drag[3]*f[3]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[17]*alpha_drag[38]+f[16]*alpha_drag[37]+f[6]*alpha_drag[27]+f[4]*alpha_drag[25]+f[2]*alpha_drag[23]+f[1]*alpha_drag[22]+f[0]*alpha_drag[21]); 
  out[5] += 0.3061862178478971*(f[17]*alpha_drag[59]+f[16]*alpha_drag[58]+f[6]*alpha_drag[48]+f[5]*alpha_drag[47]+f[2]*alpha_drag[44]+f[1]*alpha_drag[43]+f[0]*alpha_drag[42]); 
  out[7] += 0.273861278752583*(alpha_drag[1]*f[16]+f[1]*alpha_drag[16])+0.3061862178478971*(alpha_drag[3]*f[7]+alpha_drag[2]*f[6]+f[2]*alpha_drag[6]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[8] += 0.273861278752583*(alpha_drag[2]*f[17]+f[2]*alpha_drag[17])+0.3061862178478971*(alpha_drag[3]*f[8]+alpha_drag[1]*f[6]+f[1]*alpha_drag[6]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[9] += 0.273861278752583*f[1]*alpha_drag[37]+0.3061862178478971*(f[2]*alpha_drag[27]+f[9]*alpha_drag[25]+f[6]*alpha_drag[23])+0.273861278752583*f[16]*alpha_drag[22]+0.3061862178478971*(f[0]*alpha_drag[22]+f[1]*alpha_drag[21]); 
  out[10] += 0.273861278752583*f[2]*alpha_drag[38]+0.3061862178478971*(f[1]*alpha_drag[27]+f[10]*alpha_drag[25])+0.273861278752583*f[17]*alpha_drag[23]+0.3061862178478971*(f[0]*alpha_drag[23]+f[6]*alpha_drag[22]+f[2]*alpha_drag[21]); 
  out[11] += 0.3061862178478971*(f[11]*alpha_drag[25]+f[8]*alpha_drag[23]+f[7]*alpha_drag[22]+f[3]*alpha_drag[21]+alpha_drag[3]*f[11]+alpha_drag[2]*f[10]+alpha_drag[1]*f[9]+alpha_drag[0]*f[4]); 
  out[12] += 0.273861278752583*f[1]*alpha_drag[58]+0.3061862178478971*(f[2]*alpha_drag[48]+f[12]*alpha_drag[47]+f[6]*alpha_drag[44])+0.273861278752583*f[16]*alpha_drag[43]+0.3061862178478971*(f[0]*alpha_drag[43]+f[1]*alpha_drag[42]); 
  out[13] += 0.273861278752583*f[2]*alpha_drag[59]+0.3061862178478971*(f[1]*alpha_drag[48]+f[13]*alpha_drag[47])+0.273861278752583*f[17]*alpha_drag[44]+0.3061862178478971*(f[0]*alpha_drag[44]+f[6]*alpha_drag[43]+f[2]*alpha_drag[42]); 
  out[14] += 0.3061862178478971*(f[14]*alpha_drag[47]+f[8]*alpha_drag[44]+f[7]*alpha_drag[43]+f[3]*alpha_drag[42]+alpha_drag[3]*f[14]+alpha_drag[2]*f[13]+alpha_drag[1]*f[12]+alpha_drag[0]*f[5]); 
  out[15] += 0.3061862178478971*(f[15]*alpha_drag[47]+f[10]*alpha_drag[44]+f[9]*alpha_drag[43]+f[4]*alpha_drag[42]+f[15]*alpha_drag[25]+f[13]*alpha_drag[23]+f[12]*alpha_drag[22]+f[5]*alpha_drag[21]); 
  out[18] += 1.185854122563142*(alpha_diffusion[17]*f[17]+alpha_diffusion[16]*f[16]+alpha_diffusion[6]*f[6]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4+0.6123724356957944*alpha_drag[3]*f[18]+0.6846531968814573*(alpha_drag[2]*f[8]+alpha_drag[1]*f[7]+alpha_drag[0]*f[3]+f[0]*alpha_drag[3]); 
  out[19] += 1.185854122563142*(alpha_diffusion[17]*f[17]+alpha_diffusion[16]*f[16]+alpha_diffusion[6]*f[6]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvySq4+0.6123724356957944*f[19]*alpha_drag[25]+0.6846531968814573*(f[0]*alpha_drag[25]+f[10]*alpha_drag[23]+f[9]*alpha_drag[22]+f[4]*alpha_drag[21]); 
  out[20] += 1.185854122563142*(alpha_diffusion[17]*f[17]+alpha_diffusion[16]*f[16]+alpha_diffusion[6]*f[6]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvzSq4+0.6123724356957944*f[20]*alpha_drag[47]+0.6846531968814573*(f[0]*alpha_drag[47]+f[13]*alpha_drag[44]+f[12]*alpha_drag[43]+f[5]*alpha_drag[42]); 

  return alpha_mid; 

} 
