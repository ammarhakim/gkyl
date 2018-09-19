#include <VmLBOModDecl.h> 
double VmLBOconstNuVol3x3vMaxP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdvx2 = 2/dxv[3]; 
  const double rdvxSq4 = 4/(dxv[3]*dxv[3]); 
  const double rdvy2 = 2/dxv[4]; 
  const double rdvySq4 = 4/(dxv[4]*dxv[4]); 
  const double rdvz2 = 2/dxv[5]; 
  const double rdvzSq4 = 4/(dxv[5]*dxv[5]); 

  double alpha_mid = 0.0; 
  double alpha_drag[21]; 
  double alpha_diffusion[7]; 

  // Expand rdv2*nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.828427124746191*nuU[0]-8.0*w[3]*nu)*rdvx2; 
  alpha_drag[1] = 2.828427124746191*nuU[1]*rdvx2; 
  alpha_drag[2] = 2.828427124746191*nuU[2]*rdvx2; 
  alpha_drag[3] = 2.828427124746191*nuU[3]*rdvx2; 
  alpha_drag[4] = -2.309401076758503*dxv[3]*nu*rdvx2; 

  // cvard[1]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0625*alpha_drag[0]); 

  // Expand rdv2*nu*(vy-uy) in phase basis.
  alpha_drag[7] = (2.828427124746191*nuU[4]-8.0*w[4]*nu)*rdvy2; 
  alpha_drag[8] = 2.828427124746191*nuU[5]*rdvy2; 
  alpha_drag[9] = 2.828427124746191*nuU[6]*rdvy2; 
  alpha_drag[10] = 2.828427124746191*nuU[7]*rdvy2; 
  alpha_drag[12] = -2.309401076758503*dxv[4]*nu*rdvy2; 

  // cvard[2]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0625*alpha_drag[7]); 

  // Expand rdv2*nu*(vz-uz) in phase basis.
  alpha_drag[14] = (2.828427124746191*nuU[8]-8.0*w[5]*nu)*rdvz2; 
  alpha_drag[15] = 2.828427124746191*nuU[9]*rdvz2; 
  alpha_drag[16] = 2.828427124746191*nuU[10]*rdvz2; 
  alpha_drag[17] = 2.828427124746191*nuU[11]*rdvz2; 
  alpha_drag[20] = -2.309401076758503*dxv[5]*nu*rdvz2; 

  // cvard[3]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0625*alpha_drag[14]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.828427124746191*nuVtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*nuVtSq[1]; 
  alpha_diffusion[2] = 2.828427124746191*nuVtSq[2]; 
  alpha_diffusion[3] = 2.828427124746191*nuVtSq[3]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.1666666666666667*alpha_diffusion[0]*rdvxSq4); 
  alpha_mid += std::abs(0.1666666666666667*alpha_diffusion[0]*rdvySq4); 
  alpha_mid += std::abs(0.1666666666666667*alpha_diffusion[0]*rdvzSq4); 

  // Put together updates due to drag and diffusion terms.
  out[4] += 0.2165063509461096*(alpha_drag[4]*f[4]+alpha_drag[3]*f[3]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[5] += 0.2165063509461096*(f[5]*alpha_drag[12]+f[3]*alpha_drag[10]+f[2]*alpha_drag[9]+f[1]*alpha_drag[8]+f[0]*alpha_drag[7]); 
  out[6] += 0.2165063509461096*(f[6]*alpha_drag[20]+f[3]*alpha_drag[17]+f[2]*alpha_drag[16]+f[1]*alpha_drag[15]+f[0]*alpha_drag[14]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol3x3vMaxP2(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdvx2 = 2/dxv[3]; 
  const double rdvxSq4 = 4/(dxv[3]*dxv[3]); 
  const double rdvy2 = 2/dxv[4]; 
  const double rdvySq4 = 4/(dxv[4]*dxv[4]); 
  const double rdvz2 = 2/dxv[5]; 
  const double rdvzSq4 = 4/(dxv[5]*dxv[5]); 

  double alpha_mid = 0.0; 
  double alpha_drag[84]; 
  double alpha_diffusion[28]; 

  // Expand rdv2*nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.828427124746191*nuU[0]-8.0*w[3]*nu)*rdvx2; 
  alpha_drag[1] = 2.828427124746191*nuU[1]*rdvx2; 
  alpha_drag[2] = 2.828427124746191*nuU[2]*rdvx2; 
  alpha_drag[3] = 2.828427124746191*nuU[3]*rdvx2; 
  alpha_drag[4] = -2.309401076758503*dxv[3]*nu*rdvx2; 
  alpha_drag[7] = 2.828427124746191*nuU[4]*rdvx2; 
  alpha_drag[8] = 2.828427124746191*nuU[5]*rdvx2; 
  alpha_drag[9] = 2.828427124746191*nuU[6]*rdvx2; 
  alpha_drag[22] = 2.828427124746191*nuU[7]*rdvx2; 
  alpha_drag[23] = 2.828427124746191*nuU[8]*rdvx2; 
  alpha_drag[24] = 2.828427124746191*nuU[9]*rdvx2; 

  // cvard[1]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0625*alpha_drag[0]-0.06987712429686843*(alpha_drag[24]+alpha_drag[23]+alpha_drag[22])); 

  // Expand rdv2*nu*(vy-uy) in phase basis.
  alpha_drag[28] = (2.828427124746191*nuU[10]-8.0*w[4]*nu)*rdvy2; 
  alpha_drag[29] = 2.828427124746191*nuU[11]*rdvy2; 
  alpha_drag[30] = 2.828427124746191*nuU[12]*rdvy2; 
  alpha_drag[31] = 2.828427124746191*nuU[13]*rdvy2; 
  alpha_drag[33] = -2.309401076758503*dxv[4]*nu*rdvy2; 
  alpha_drag[35] = 2.828427124746191*nuU[14]*rdvy2; 
  alpha_drag[36] = 2.828427124746191*nuU[15]*rdvy2; 
  alpha_drag[37] = 2.828427124746191*nuU[16]*rdvy2; 
  alpha_drag[50] = 2.828427124746191*nuU[17]*rdvy2; 
  alpha_drag[51] = 2.828427124746191*nuU[18]*rdvy2; 
  alpha_drag[52] = 2.828427124746191*nuU[19]*rdvy2; 

  // cvard[2]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0625*alpha_drag[28]-0.06987712429686843*(alpha_drag[52]+alpha_drag[51]+alpha_drag[50])); 

  // Expand rdv2*nu*(vz-uz) in phase basis.
  alpha_drag[56] = (2.828427124746191*nuU[20]-8.0*w[5]*nu)*rdvz2; 
  alpha_drag[57] = 2.828427124746191*nuU[21]*rdvz2; 
  alpha_drag[58] = 2.828427124746191*nuU[22]*rdvz2; 
  alpha_drag[59] = 2.828427124746191*nuU[23]*rdvz2; 
  alpha_drag[62] = -2.309401076758503*dxv[5]*nu*rdvz2; 
  alpha_drag[63] = 2.828427124746191*nuU[24]*rdvz2; 
  alpha_drag[64] = 2.828427124746191*nuU[25]*rdvz2; 
  alpha_drag[65] = 2.828427124746191*nuU[26]*rdvz2; 
  alpha_drag[78] = 2.828427124746191*nuU[27]*rdvz2; 
  alpha_drag[79] = 2.828427124746191*nuU[28]*rdvz2; 
  alpha_drag[80] = 2.828427124746191*nuU[29]*rdvz2; 

  // cvard[3]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.0625*alpha_drag[56]-0.06987712429686843*(alpha_drag[80]+alpha_drag[79]+alpha_drag[78])); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.828427124746191*nuVtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*nuVtSq[1]; 
  alpha_diffusion[2] = 2.828427124746191*nuVtSq[2]; 
  alpha_diffusion[3] = 2.828427124746191*nuVtSq[3]; 
  alpha_diffusion[7] = 2.828427124746191*nuVtSq[4]; 
  alpha_diffusion[8] = 2.828427124746191*nuVtSq[5]; 
  alpha_diffusion[9] = 2.828427124746191*nuVtSq[6]; 
  alpha_diffusion[22] = 2.828427124746191*nuVtSq[7]; 
  alpha_diffusion[23] = 2.828427124746191*nuVtSq[8]; 
  alpha_diffusion[24] = 2.828427124746191*nuVtSq[9]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((0.225*alpha_diffusion[0]-0.2515576474687263*(alpha_diffusion[24]+alpha_diffusion[23]+alpha_diffusion[22]))*rdvxSq4); 
  alpha_mid += std::abs((0.225*alpha_diffusion[0]-0.2515576474687263*(alpha_diffusion[24]+alpha_diffusion[23]+alpha_diffusion[22]))*rdvySq4); 
  alpha_mid += std::abs((0.225*alpha_diffusion[0]-0.2515576474687263*(alpha_diffusion[24]+alpha_diffusion[23]+alpha_diffusion[22]))*rdvzSq4); 

  // Put together updates due to drag and diffusion terms.
  out[4] += 0.2165063509461096*(alpha_drag[24]*f[24]+alpha_drag[23]*f[23]+alpha_drag[22]*f[22]+alpha_drag[9]*f[9]+alpha_drag[8]*f[8]+alpha_drag[7]*f[7]+alpha_drag[4]*f[4]+alpha_drag[3]*f[3]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[5] += 0.2165063509461096*(f[24]*alpha_drag[52]+f[23]*alpha_drag[51]+f[22]*alpha_drag[50]+f[9]*alpha_drag[37]+f[8]*alpha_drag[36]+f[7]*alpha_drag[35]+f[5]*alpha_drag[33]+f[3]*alpha_drag[31]+f[2]*alpha_drag[30]+f[1]*alpha_drag[29]+f[0]*alpha_drag[28]); 
  out[6] += 0.2165063509461096*(f[24]*alpha_drag[80]+f[23]*alpha_drag[79]+f[22]*alpha_drag[78]+f[9]*alpha_drag[65]+f[8]*alpha_drag[64]+f[7]*alpha_drag[63]+f[6]*alpha_drag[62]+f[3]*alpha_drag[59]+f[2]*alpha_drag[58]+f[1]*alpha_drag[57]+f[0]*alpha_drag[56]); 
  out[10] += 0.1936491673103708*(alpha_drag[1]*f[22]+f[1]*alpha_drag[22])+0.2165063509461096*(alpha_drag[4]*f[10]+alpha_drag[3]*f[8]+f[3]*alpha_drag[8]+alpha_drag[2]*f[7]+f[2]*alpha_drag[7]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[11] += 0.1936491673103708*(alpha_drag[2]*f[23]+f[2]*alpha_drag[23])+0.2165063509461096*(alpha_drag[4]*f[11]+alpha_drag[3]*f[9]+f[3]*alpha_drag[9]+alpha_drag[1]*f[7]+f[1]*alpha_drag[7]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[12] += 0.1936491673103708*(alpha_drag[3]*f[24]+f[3]*alpha_drag[24])+0.2165063509461096*(alpha_drag[4]*f[12]+alpha_drag[2]*f[9]+f[2]*alpha_drag[9]+alpha_drag[1]*f[8]+f[1]*alpha_drag[8]+alpha_drag[0]*f[3]+f[0]*alpha_drag[3]); 
  out[13] += 0.1936491673103708*f[1]*alpha_drag[50]+0.2165063509461096*(f[3]*alpha_drag[36]+f[2]*alpha_drag[35]+f[13]*alpha_drag[33]+f[8]*alpha_drag[31]+f[7]*alpha_drag[30])+0.1936491673103708*f[22]*alpha_drag[29]+0.2165063509461096*(f[0]*alpha_drag[29]+f[1]*alpha_drag[28]); 
  out[14] += 0.1936491673103708*f[2]*alpha_drag[51]+0.2165063509461096*(f[3]*alpha_drag[37]+f[1]*alpha_drag[35]+f[14]*alpha_drag[33]+f[9]*alpha_drag[31])+0.1936491673103708*f[23]*alpha_drag[30]+0.2165063509461096*(f[0]*alpha_drag[30]+f[7]*alpha_drag[29]+f[2]*alpha_drag[28]); 
  out[15] += 0.1936491673103708*f[3]*alpha_drag[52]+0.2165063509461096*(f[2]*alpha_drag[37]+f[1]*alpha_drag[36]+f[15]*alpha_drag[33])+0.1936491673103708*f[24]*alpha_drag[31]+0.2165063509461096*(f[0]*alpha_drag[31]+f[9]*alpha_drag[30]+f[8]*alpha_drag[29]+f[3]*alpha_drag[28]); 
  out[16] += 0.2165063509461096*(f[16]*alpha_drag[33]+f[12]*alpha_drag[31]+f[11]*alpha_drag[30]+f[10]*alpha_drag[29]+f[4]*alpha_drag[28]+alpha_drag[4]*f[16]+alpha_drag[3]*f[15]+alpha_drag[2]*f[14]+alpha_drag[1]*f[13]+alpha_drag[0]*f[5]); 
  out[17] += 0.1936491673103708*f[1]*alpha_drag[78]+0.2165063509461096*(f[3]*alpha_drag[64]+f[2]*alpha_drag[63]+f[17]*alpha_drag[62]+f[8]*alpha_drag[59]+f[7]*alpha_drag[58])+0.1936491673103708*f[22]*alpha_drag[57]+0.2165063509461096*(f[0]*alpha_drag[57]+f[1]*alpha_drag[56]); 
  out[18] += 0.1936491673103708*f[2]*alpha_drag[79]+0.2165063509461096*(f[3]*alpha_drag[65]+f[1]*alpha_drag[63]+f[18]*alpha_drag[62]+f[9]*alpha_drag[59])+0.1936491673103708*f[23]*alpha_drag[58]+0.2165063509461096*(f[0]*alpha_drag[58]+f[7]*alpha_drag[57]+f[2]*alpha_drag[56]); 
  out[19] += 0.1936491673103708*f[3]*alpha_drag[80]+0.2165063509461096*(f[2]*alpha_drag[65]+f[1]*alpha_drag[64]+f[19]*alpha_drag[62])+0.1936491673103708*f[24]*alpha_drag[59]+0.2165063509461096*(f[0]*alpha_drag[59]+f[9]*alpha_drag[58]+f[8]*alpha_drag[57]+f[3]*alpha_drag[56]); 
  out[20] += 0.2165063509461096*(f[20]*alpha_drag[62]+f[12]*alpha_drag[59]+f[11]*alpha_drag[58]+f[10]*alpha_drag[57]+f[4]*alpha_drag[56]+alpha_drag[4]*f[20]+alpha_drag[3]*f[19]+alpha_drag[2]*f[18]+alpha_drag[1]*f[17]+alpha_drag[0]*f[6]); 
  out[21] += 0.2165063509461096*(f[21]*alpha_drag[62]+f[15]*alpha_drag[59]+f[14]*alpha_drag[58]+f[13]*alpha_drag[57]+f[5]*alpha_drag[56]+f[21]*alpha_drag[33]+f[19]*alpha_drag[31]+f[18]*alpha_drag[30]+f[17]*alpha_drag[29]+f[6]*alpha_drag[28]); 
  out[25] += 0.8385254915624212*(alpha_diffusion[24]*f[24]+alpha_diffusion[23]*f[23]+alpha_diffusion[22]*f[22]+alpha_diffusion[9]*f[9]+alpha_diffusion[8]*f[8]+alpha_diffusion[7]*f[7]+alpha_diffusion[3]*f[3]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4+0.4330127018922193*alpha_drag[4]*f[25]+0.4841229182759271*(alpha_drag[3]*f[12]+alpha_drag[2]*f[11]+alpha_drag[1]*f[10]+alpha_drag[0]*f[4]+f[0]*alpha_drag[4]); 
  out[26] += 0.8385254915624212*(alpha_diffusion[24]*f[24]+alpha_diffusion[23]*f[23]+alpha_diffusion[22]*f[22]+alpha_diffusion[9]*f[9]+alpha_diffusion[8]*f[8]+alpha_diffusion[7]*f[7]+alpha_diffusion[3]*f[3]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvySq4+0.4330127018922193*f[26]*alpha_drag[33]+0.4841229182759271*(f[0]*alpha_drag[33]+f[15]*alpha_drag[31]+f[14]*alpha_drag[30]+f[13]*alpha_drag[29]+f[5]*alpha_drag[28]); 
  out[27] += 0.8385254915624212*(alpha_diffusion[24]*f[24]+alpha_diffusion[23]*f[23]+alpha_diffusion[22]*f[22]+alpha_diffusion[9]*f[9]+alpha_diffusion[8]*f[8]+alpha_diffusion[7]*f[7]+alpha_diffusion[3]*f[3]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvzSq4+0.4330127018922193*f[27]*alpha_drag[62]+0.4841229182759271*(f[0]*alpha_drag[62]+f[19]*alpha_drag[59]+f[18]*alpha_drag[58]+f[17]*alpha_drag[57]+f[6]*alpha_drag[56]); 

  return alpha_mid; 

} 
