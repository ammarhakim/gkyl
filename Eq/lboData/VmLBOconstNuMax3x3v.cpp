#include <VmLBOModDecl.h> 
double VmLBOconstNuVol3x3vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[6]:   Cell-center coordinates. 
  // dxv[6]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[3]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[3]*dxv[3]); 
  const double rdvy2nu = 2.0*nu/dxv[4]; 
  const double rdvySq4nu = 4.0*nu/(dxv[4]*dxv[4]); 
  const double rdvz2nu = 2.0*nu/dxv[5]; 
  const double rdvzSq4nu = 4.0*nu/(dxv[5]*dxv[5]); 

  double alphaDrag[21]; 
  // Expand rdv2nu*(vx-ux) in phase basis.
  alphaDrag[0] = (2.828427124746191*u[0]-8.0*w[3])*rdvx2nu; 
  alphaDrag[1] = 2.828427124746191*u[1]*rdvx2nu; 
  alphaDrag[2] = 2.828427124746191*u[2]*rdvx2nu; 
  alphaDrag[3] = 2.828427124746191*u[3]*rdvx2nu; 
  alphaDrag[4] = -2.309401076758503*dxv[3]*rdvx2nu; 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alphaDrag[7] = (2.828427124746191*u[4]-8.0*w[4])*rdvy2nu; 
  alphaDrag[8] = 2.828427124746191*u[5]*rdvy2nu; 
  alphaDrag[9] = 2.828427124746191*u[6]*rdvy2nu; 
  alphaDrag[10] = 2.828427124746191*u[7]*rdvy2nu; 
  alphaDrag[12] = -2.309401076758503*dxv[4]*rdvy2nu; 

  // Expand rdv2nu*(vz-uz) in phase basis.
  alphaDrag[14] = (2.828427124746191*u[8]-8.0*w[5])*rdvz2nu; 
  alphaDrag[15] = 2.828427124746191*u[9]*rdvz2nu; 
  alphaDrag[16] = 2.828427124746191*u[10]*rdvz2nu; 
  alphaDrag[17] = 2.828427124746191*u[11]*rdvz2nu; 
  alphaDrag[20] = -2.309401076758503*dxv[5]*rdvz2nu; 

  // Put together updates due to drag and diffusion terms.
  out[4] += 0.2165063509461096*(alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[5] += 0.2165063509461096*(f[5]*alphaDrag[12]+f[3]*alphaDrag[10]+f[2]*alphaDrag[9]+f[1]*alphaDrag[8]+f[0]*alphaDrag[7]); 
  out[6] += 0.2165063509461096*(f[6]*alphaDrag[20]+f[3]*alphaDrag[17]+f[2]*alphaDrag[16]+f[1]*alphaDrag[15]+f[0]*alphaDrag[14]); 

  return std::abs(0.0625*alphaDrag[0])+std::abs(0.0625*alphaDrag[7])+std::abs(0.0625*alphaDrag[14])+std::abs(0.4714045207910317*vtSq[0]*rdvxSq4nu)+std::abs(0.4714045207910317*vtSq[0]*rdvySq4nu)+std::abs(0.4714045207910317*vtSq[0]*rdvzSq4nu); 

} 
double VmLBOconstNuVol3x3vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[6]:   Cell-center coordinates. 
  // dxv[6]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[3]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[3]*dxv[3]); 
  const double rdvy2nu = 2.0*nu/dxv[4]; 
  const double rdvySq4nu = 4.0*nu/(dxv[4]*dxv[4]); 
  const double rdvz2nu = 2.0*nu/dxv[5]; 
  const double rdvzSq4nu = 4.0*nu/(dxv[5]*dxv[5]); 

  double alphaDrag[84]; 
  // Expand rdv2nu*(vx-ux) in phase basis.
  alphaDrag[0] = (2.828427124746191*u[0]-8.0*w[3])*rdvx2nu; 
  alphaDrag[1] = 2.828427124746191*u[1]*rdvx2nu; 
  alphaDrag[2] = 2.828427124746191*u[2]*rdvx2nu; 
  alphaDrag[3] = 2.828427124746191*u[3]*rdvx2nu; 
  alphaDrag[4] = -2.309401076758503*dxv[3]*rdvx2nu; 
  alphaDrag[7] = 2.828427124746191*u[4]*rdvx2nu; 
  alphaDrag[8] = 2.828427124746191*u[5]*rdvx2nu; 
  alphaDrag[9] = 2.828427124746191*u[6]*rdvx2nu; 
  alphaDrag[22] = 2.828427124746191*u[7]*rdvx2nu; 
  alphaDrag[23] = 2.828427124746191*u[8]*rdvx2nu; 
  alphaDrag[24] = 2.828427124746191*u[9]*rdvx2nu; 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alphaDrag[28] = (2.828427124746191*u[10]-8.0*w[4])*rdvy2nu; 
  alphaDrag[29] = 2.828427124746191*u[11]*rdvy2nu; 
  alphaDrag[30] = 2.828427124746191*u[12]*rdvy2nu; 
  alphaDrag[31] = 2.828427124746191*u[13]*rdvy2nu; 
  alphaDrag[33] = -2.309401076758503*dxv[4]*rdvy2nu; 
  alphaDrag[35] = 2.828427124746191*u[14]*rdvy2nu; 
  alphaDrag[36] = 2.828427124746191*u[15]*rdvy2nu; 
  alphaDrag[37] = 2.828427124746191*u[16]*rdvy2nu; 
  alphaDrag[50] = 2.828427124746191*u[17]*rdvy2nu; 
  alphaDrag[51] = 2.828427124746191*u[18]*rdvy2nu; 
  alphaDrag[52] = 2.828427124746191*u[19]*rdvy2nu; 

  // Expand rdv2nu*(vz-uz) in phase basis.
  alphaDrag[56] = (2.828427124746191*u[20]-8.0*w[5])*rdvz2nu; 
  alphaDrag[57] = 2.828427124746191*u[21]*rdvz2nu; 
  alphaDrag[58] = 2.828427124746191*u[22]*rdvz2nu; 
  alphaDrag[59] = 2.828427124746191*u[23]*rdvz2nu; 
  alphaDrag[62] = -2.309401076758503*dxv[5]*rdvz2nu; 
  alphaDrag[63] = 2.828427124746191*u[24]*rdvz2nu; 
  alphaDrag[64] = 2.828427124746191*u[25]*rdvz2nu; 
  alphaDrag[65] = 2.828427124746191*u[26]*rdvz2nu; 
  alphaDrag[78] = 2.828427124746191*u[27]*rdvz2nu; 
  alphaDrag[79] = 2.828427124746191*u[28]*rdvz2nu; 
  alphaDrag[80] = 2.828427124746191*u[29]*rdvz2nu; 

  double facDiff[10]; 
  // Expand nu*vthSq in phase basis.
  facDiff[0] = vtSq[0]; 
  facDiff[1] = vtSq[1]; 
  facDiff[2] = vtSq[2]; 
  facDiff[3] = vtSq[3]; 
  facDiff[4] = vtSq[4]; 
  facDiff[5] = vtSq[5]; 
  facDiff[6] = vtSq[6]; 
  facDiff[7] = vtSq[7]; 
  facDiff[8] = vtSq[8]; 
  facDiff[9] = vtSq[9]; 

  // Put together updates due to drag and diffusion terms.
  out[4] += 0.2165063509461096*(alphaDrag[24]*f[24]+alphaDrag[23]*f[23]+alphaDrag[22]*f[22]+alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[5] += 0.2165063509461096*(f[24]*alphaDrag[52]+f[23]*alphaDrag[51]+f[22]*alphaDrag[50]+f[9]*alphaDrag[37]+f[8]*alphaDrag[36]+f[7]*alphaDrag[35]+f[5]*alphaDrag[33]+f[3]*alphaDrag[31]+f[2]*alphaDrag[30]+f[1]*alphaDrag[29]+f[0]*alphaDrag[28]); 
  out[6] += 0.2165063509461096*(f[24]*alphaDrag[80]+f[23]*alphaDrag[79]+f[22]*alphaDrag[78]+f[9]*alphaDrag[65]+f[8]*alphaDrag[64]+f[7]*alphaDrag[63]+f[6]*alphaDrag[62]+f[3]*alphaDrag[59]+f[2]*alphaDrag[58]+f[1]*alphaDrag[57]+f[0]*alphaDrag[56]); 
  out[10] += 0.1936491673103708*(alphaDrag[1]*f[22]+f[1]*alphaDrag[22])+0.2165063509461096*(alphaDrag[4]*f[10]+alphaDrag[3]*f[8]+f[3]*alphaDrag[8]+alphaDrag[2]*f[7]+f[2]*alphaDrag[7]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[11] += 0.1936491673103708*(alphaDrag[2]*f[23]+f[2]*alphaDrag[23])+0.2165063509461096*(alphaDrag[4]*f[11]+alphaDrag[3]*f[9]+f[3]*alphaDrag[9]+alphaDrag[1]*f[7]+f[1]*alphaDrag[7]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[12] += 0.1936491673103708*(alphaDrag[3]*f[24]+f[3]*alphaDrag[24])+0.2165063509461096*(alphaDrag[4]*f[12]+alphaDrag[2]*f[9]+f[2]*alphaDrag[9]+alphaDrag[1]*f[8]+f[1]*alphaDrag[8]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[13] += 0.1936491673103708*f[1]*alphaDrag[50]+0.2165063509461096*(f[3]*alphaDrag[36]+f[2]*alphaDrag[35]+f[13]*alphaDrag[33]+f[8]*alphaDrag[31]+f[7]*alphaDrag[30])+0.1936491673103708*f[22]*alphaDrag[29]+0.2165063509461096*(f[0]*alphaDrag[29]+f[1]*alphaDrag[28]); 
  out[14] += 0.1936491673103708*f[2]*alphaDrag[51]+0.2165063509461096*(f[3]*alphaDrag[37]+f[1]*alphaDrag[35]+f[14]*alphaDrag[33]+f[9]*alphaDrag[31])+0.1936491673103708*f[23]*alphaDrag[30]+0.2165063509461096*(f[0]*alphaDrag[30]+f[7]*alphaDrag[29]+f[2]*alphaDrag[28]); 
  out[15] += 0.1936491673103708*f[3]*alphaDrag[52]+0.2165063509461096*(f[2]*alphaDrag[37]+f[1]*alphaDrag[36]+f[15]*alphaDrag[33])+0.1936491673103708*f[24]*alphaDrag[31]+0.2165063509461096*(f[0]*alphaDrag[31]+f[9]*alphaDrag[30]+f[8]*alphaDrag[29]+f[3]*alphaDrag[28]); 
  out[16] += 0.2165063509461096*(f[16]*alphaDrag[33]+f[12]*alphaDrag[31]+f[11]*alphaDrag[30]+f[10]*alphaDrag[29]+f[4]*alphaDrag[28]+alphaDrag[4]*f[16]+alphaDrag[3]*f[15]+alphaDrag[2]*f[14]+alphaDrag[1]*f[13]+alphaDrag[0]*f[5]); 
  out[17] += 0.1936491673103708*f[1]*alphaDrag[78]+0.2165063509461096*(f[3]*alphaDrag[64]+f[2]*alphaDrag[63]+f[17]*alphaDrag[62]+f[8]*alphaDrag[59]+f[7]*alphaDrag[58])+0.1936491673103708*f[22]*alphaDrag[57]+0.2165063509461096*(f[0]*alphaDrag[57]+f[1]*alphaDrag[56]); 
  out[18] += 0.1936491673103708*f[2]*alphaDrag[79]+0.2165063509461096*(f[3]*alphaDrag[65]+f[1]*alphaDrag[63]+f[18]*alphaDrag[62]+f[9]*alphaDrag[59])+0.1936491673103708*f[23]*alphaDrag[58]+0.2165063509461096*(f[0]*alphaDrag[58]+f[7]*alphaDrag[57]+f[2]*alphaDrag[56]); 
  out[19] += 0.1936491673103708*f[3]*alphaDrag[80]+0.2165063509461096*(f[2]*alphaDrag[65]+f[1]*alphaDrag[64]+f[19]*alphaDrag[62])+0.1936491673103708*f[24]*alphaDrag[59]+0.2165063509461096*(f[0]*alphaDrag[59]+f[9]*alphaDrag[58]+f[8]*alphaDrag[57]+f[3]*alphaDrag[56]); 
  out[20] += 0.2165063509461096*(f[20]*alphaDrag[62]+f[12]*alphaDrag[59]+f[11]*alphaDrag[58]+f[10]*alphaDrag[57]+f[4]*alphaDrag[56]+alphaDrag[4]*f[20]+alphaDrag[3]*f[19]+alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[0]*f[6]); 
  out[21] += 0.2165063509461096*(f[21]*alphaDrag[62]+f[15]*alphaDrag[59]+f[14]*alphaDrag[58]+f[13]*alphaDrag[57]+f[5]*alphaDrag[56]+f[21]*alphaDrag[33]+f[19]*alphaDrag[31]+f[18]*alphaDrag[30]+f[17]*alphaDrag[29]+f[6]*alphaDrag[28]); 
  out[25] += 2.371708245126284*(facDiff[9]*f[24]+facDiff[8]*f[23]+facDiff[7]*f[22]+facDiff[6]*f[9]+facDiff[5]*f[8]+facDiff[4]*f[7]+f[3]*facDiff[3]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4nu+0.4330127018922193*alphaDrag[4]*f[25]+0.4841229182759271*(alphaDrag[3]*f[12]+alphaDrag[2]*f[11]+alphaDrag[1]*f[10]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]); 
  out[26] += 2.371708245126284*(facDiff[9]*f[24]+facDiff[8]*f[23]+facDiff[7]*f[22]+facDiff[6]*f[9]+facDiff[5]*f[8]+facDiff[4]*f[7]+f[3]*facDiff[3]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4nu+0.4330127018922193*f[26]*alphaDrag[33]+0.4841229182759271*(f[0]*alphaDrag[33]+f[15]*alphaDrag[31]+f[14]*alphaDrag[30]+f[13]*alphaDrag[29]+f[5]*alphaDrag[28]); 
  out[27] += 2.371708245126284*(facDiff[9]*f[24]+facDiff[8]*f[23]+facDiff[7]*f[22]+facDiff[6]*f[9]+facDiff[5]*f[8]+facDiff[4]*f[7]+f[3]*facDiff[3]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvzSq4nu+0.4330127018922193*f[27]*alphaDrag[62]+0.4841229182759271*(f[0]*alphaDrag[62]+f[19]*alphaDrag[59]+f[18]*alphaDrag[58]+f[17]*alphaDrag[57]+f[6]*alphaDrag[56]); 

  return std::abs(0.0625*alphaDrag[0]-0.06987712429686843*(alphaDrag[24]+alphaDrag[23]+alphaDrag[22]))+std::abs(0.0625*alphaDrag[28]-0.06987712429686843*(alphaDrag[52]+alphaDrag[51]+alphaDrag[50]))+std::abs(0.0625*alphaDrag[56]-0.06987712429686843*(alphaDrag[80]+alphaDrag[79]+alphaDrag[78]))+std::abs((0.6363961030678926*facDiff[0]-0.711512473537885*(facDiff[9]+facDiff[8]+facDiff[7]))*rdvxSq4nu)+std::abs((0.6363961030678926*facDiff[0]-0.711512473537885*(facDiff[9]+facDiff[8]+facDiff[7]))*rdvySq4nu)+std::abs((0.6363961030678926*facDiff[0]-0.711512473537885*(facDiff[9]+facDiff[8]+facDiff[7]))*rdvzSq4nu); 

} 
