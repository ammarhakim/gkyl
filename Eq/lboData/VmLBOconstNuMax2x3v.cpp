#include <VmLBOModDecl.h> 
double VmLBOconstNuVol2x3vMaxP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdv20 = 2/dxv[2]; 
  const double rdvSq40 = 4/(dxv[2]*dxv[2]); 
  const double rdv21 = 2/dxv[3]; 
  const double rdvSq41 = 4/(dxv[3]*dxv[3]); 
  const double rdv22 = 2/dxv[4]; 
  const double rdvSq42 = 4/(dxv[4]*dxv[4]); 

  double alpha_mid = 0.0; 
  double alpha_drag[18]; 
  double alpha_diffusion[6]; 

  alpha_drag[0] = (2.828427124746191*nuU[0]-5.656854249492382*w[2]*nu)*rdv20; 
  alpha_drag[1] = 2.828427124746191*nuU[1]*rdv20; 
  alpha_drag[2] = 2.828427124746191*nuU[2]*rdv20; 
  alpha_drag[3] = -1.632993161855453*dxv[2]*nu*rdv20; 
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[0]); 
  alpha_drag[6] = (2.828427124746191*nuU[6]-5.656854249492382*w[3]*nu)*rdv21; 
  alpha_drag[7] = 2.828427124746191*nuU[7]*rdv21; 
  alpha_drag[8] = 2.828427124746191*nuU[8]*rdv21; 
  alpha_drag[10] = -1.632993161855453*dxv[3]*nu*rdv21; 
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[6]); 
  alpha_drag[12] = (2.828427124746191*nuU[12]-5.656854249492382*w[4]*nu)*rdv22; 
  alpha_drag[13] = 2.828427124746191*nuU[13]*rdv22; 
  alpha_drag[14] = 2.828427124746191*nuU[14]*rdv22; 
  alpha_drag[17] = -1.632993161855453*dxv[4]*nu*rdv22; 
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[12]); 
  alpha_diffusion[0] = 2.828427124746191*nuVtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*nuVtSq[1]; 
  alpha_diffusion[2] = 2.828427124746191*nuVtSq[2]; 
  alpha_mid += std::abs(0.0441941738241592*alpha_diffusion[0]*rdvSq40); 
  alpha_mid += std::abs(0.0441941738241592*alpha_diffusion[0]*rdvSq41); 
  alpha_mid += std::abs(0.0441941738241592*alpha_diffusion[0]*rdvSq42); 
  out[3] += 0.3061862178478971*(alpha_drag[3]*f[3]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[4]*alpha_drag[10]+f[2]*alpha_drag[8]+f[1]*alpha_drag[7]+f[0]*alpha_drag[6]); 
  out[5] += 0.3061862178478971*(f[5]*alpha_drag[17]+f[2]*alpha_drag[14]+f[1]*alpha_drag[13]+f[0]*alpha_drag[12]); 
return alpha_mid; 

} 
double VmLBOconstNuVol2x3vMaxP2(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdv20 = 2/dxv[2]; 
  const double rdvSq40 = 4/(dxv[2]*dxv[2]); 
  const double rdv21 = 2/dxv[3]; 
  const double rdvSq41 = 4/(dxv[3]*dxv[3]); 
  const double rdv22 = 2/dxv[4]; 
  const double rdvSq42 = 4/(dxv[4]*dxv[4]); 

  double alpha_mid = 0.0; 
  double alpha_drag[63]; 
  double alpha_diffusion[21]; 

  alpha_drag[0] = (2.828427124746191*nuU[0]-5.656854249492382*w[2]*nu)*rdv20; 
  alpha_drag[1] = 2.828427124746191*nuU[1]*rdv20; 
  alpha_drag[2] = 2.828427124746191*nuU[2]*rdv20; 
  alpha_drag[3] = -1.632993161855453*dxv[2]*nu*rdv20; 
  alpha_drag[6] = 2.828427124746191*nuU[3]*rdv20; 
  alpha_drag[16] = 2.828427124746191*nuU[4]*rdv20; 
  alpha_drag[17] = 2.828427124746191*nuU[5]*rdv20; 
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[0]-0.09882117688026182*(alpha_drag[17]+alpha_drag[16])); 
  alpha_drag[21] = (2.828427124746191*nuU[21]-5.656854249492382*w[3]*nu)*rdv21; 
  alpha_drag[22] = 2.828427124746191*nuU[22]*rdv21; 
  alpha_drag[23] = 2.828427124746191*nuU[23]*rdv21; 
  alpha_drag[25] = -1.632993161855453*dxv[3]*nu*rdv21; 
  alpha_drag[27] = 2.828427124746191*nuU[24]*rdv21; 
  alpha_drag[37] = 2.828427124746191*nuU[25]*rdv21; 
  alpha_drag[38] = 2.828427124746191*nuU[26]*rdv21; 
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[21]-0.09882117688026182*(alpha_drag[38]+alpha_drag[37])); 
  alpha_drag[42] = (2.828427124746191*nuU[42]-5.656854249492382*w[4]*nu)*rdv22; 
  alpha_drag[43] = 2.828427124746191*nuU[43]*rdv22; 
  alpha_drag[44] = 2.828427124746191*nuU[44]*rdv22; 
  alpha_drag[47] = -1.632993161855453*dxv[4]*nu*rdv22; 
  alpha_drag[48] = 2.828427124746191*nuU[45]*rdv22; 
  alpha_drag[58] = 2.828427124746191*nuU[46]*rdv22; 
  alpha_drag[59] = 2.828427124746191*nuU[47]*rdv22; 
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[42]-0.09882117688026182*(alpha_drag[59]+alpha_drag[58])); 
  alpha_diffusion[0] = 2.828427124746191*nuVtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*nuVtSq[1]; 
  alpha_diffusion[2] = 2.828427124746191*nuVtSq[2]; 
  alpha_diffusion[6] = 2.828427124746191*nuVtSq[3]; 
  alpha_diffusion[16] = 2.828427124746191*nuVtSq[4]; 
  alpha_diffusion[17] = 2.828427124746191*nuVtSq[5]; 
  alpha_mid += std::abs((0.0441941738241592*alpha_diffusion[0]-0.04941058844013091*(alpha_diffusion[17]+alpha_diffusion[16]))*rdvSq40); 
  alpha_mid += std::abs((0.0441941738241592*alpha_diffusion[0]-0.04941058844013091*(alpha_diffusion[17]+alpha_diffusion[16]))*rdvSq41); 
  alpha_mid += std::abs((0.0441941738241592*alpha_diffusion[0]-0.04941058844013091*(alpha_diffusion[17]+alpha_diffusion[16]))*rdvSq42); 
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
  out[18] += 1.185854122563142*(alpha_diffusion[17]*f[17]+alpha_diffusion[16]*f[16]+alpha_diffusion[6]*f[6]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq40+0.6123724356957944*alpha_drag[3]*f[18]+0.6846531968814573*(alpha_drag[2]*f[8]+alpha_drag[1]*f[7]+alpha_drag[0]*f[3]+f[0]*alpha_drag[3]); 
  out[19] += 1.185854122563142*(alpha_diffusion[17]*f[17]+alpha_diffusion[16]*f[16]+alpha_diffusion[6]*f[6]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq41+0.6123724356957944*f[19]*alpha_drag[25]+0.6846531968814573*(f[0]*alpha_drag[25]+f[10]*alpha_drag[23]+f[9]*alpha_drag[22]+f[4]*alpha_drag[21]); 
  out[20] += 1.185854122563142*(alpha_diffusion[17]*f[17]+alpha_diffusion[16]*f[16]+alpha_diffusion[6]*f[6]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq42+0.6123724356957944*f[20]*alpha_drag[47]+0.6846531968814573*(f[0]*alpha_drag[47]+f[13]*alpha_drag[44]+f[12]*alpha_drag[43]+f[5]*alpha_drag[42]); 
return alpha_mid; 

} 
double VmLBOconstNuVol2x3vMaxP3(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdv20 = 2/dxv[2]; 
  const double rdvSq40 = 4/(dxv[2]*dxv[2]); 
  const double rdv21 = 2/dxv[3]; 
  const double rdvSq41 = 4/(dxv[3]*dxv[3]); 
  const double rdv22 = 2/dxv[4]; 
  const double rdvSq42 = 4/(dxv[4]*dxv[4]); 

  double alpha_mid = 0.0; 
  double alpha_drag[168]; 
  double alpha_diffusion[56]; 

  alpha_drag[0] = (2.828427124746191*nuU[0]-5.656854249492382*w[2]*nu)*rdv20; 
  alpha_drag[1] = 2.828427124746191*nuU[1]*rdv20; 
  alpha_drag[2] = 2.828427124746191*nuU[2]*rdv20; 
  alpha_drag[3] = -1.632993161855453*dxv[2]*nu*rdv20; 
  alpha_drag[6] = 2.828427124746191*nuU[3]*rdv20; 
  alpha_drag[16] = 2.828427124746191*nuU[4]*rdv20; 
  alpha_drag[17] = 2.828427124746191*nuU[5]*rdv20; 
  alpha_drag[31] = 2.828427124746191*nuU[6]*rdv20; 
  alpha_drag[32] = 2.828427124746191*nuU[7]*rdv20; 
  alpha_drag[51] = 2.828427124746191*nuU[8]*rdv20; 
  alpha_drag[52] = 2.828427124746191*nuU[9]*rdv20; 
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[0]-0.09882117688026182*(alpha_drag[17]+alpha_drag[16])); 
  alpha_drag[56] = (2.828427124746191*nuU[56]-5.656854249492382*w[3]*nu)*rdv21; 
  alpha_drag[57] = 2.828427124746191*nuU[57]*rdv21; 
  alpha_drag[58] = 2.828427124746191*nuU[58]*rdv21; 
  alpha_drag[60] = -1.632993161855453*dxv[3]*nu*rdv21; 
  alpha_drag[62] = 2.828427124746191*nuU[59]*rdv21; 
  alpha_drag[72] = 2.828427124746191*nuU[60]*rdv21; 
  alpha_drag[73] = 2.828427124746191*nuU[61]*rdv21; 
  alpha_drag[87] = 2.828427124746191*nuU[62]*rdv21; 
  alpha_drag[88] = 2.828427124746191*nuU[63]*rdv21; 
  alpha_drag[107] = 2.828427124746191*nuU[64]*rdv21; 
  alpha_drag[108] = 2.828427124746191*nuU[65]*rdv21; 
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[56]-0.09882117688026182*(alpha_drag[73]+alpha_drag[72])); 
  alpha_drag[112] = (2.828427124746191*nuU[112]-5.656854249492382*w[4]*nu)*rdv22; 
  alpha_drag[113] = 2.828427124746191*nuU[113]*rdv22; 
  alpha_drag[114] = 2.828427124746191*nuU[114]*rdv22; 
  alpha_drag[117] = -1.632993161855453*dxv[4]*nu*rdv22; 
  alpha_drag[118] = 2.828427124746191*nuU[115]*rdv22; 
  alpha_drag[128] = 2.828427124746191*nuU[116]*rdv22; 
  alpha_drag[129] = 2.828427124746191*nuU[117]*rdv22; 
  alpha_drag[143] = 2.828427124746191*nuU[118]*rdv22; 
  alpha_drag[144] = 2.828427124746191*nuU[119]*rdv22; 
  alpha_drag[163] = 2.828427124746191*nuU[120]*rdv22; 
  alpha_drag[164] = 2.828427124746191*nuU[121]*rdv22; 
  alpha_mid += std::abs(0.0883883476483184*alpha_drag[112]-0.09882117688026182*(alpha_drag[129]+alpha_drag[128])); 
  alpha_diffusion[0] = 2.828427124746191*nuVtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*nuVtSq[1]; 
  alpha_diffusion[2] = 2.828427124746191*nuVtSq[2]; 
  alpha_diffusion[6] = 2.828427124746191*nuVtSq[3]; 
  alpha_diffusion[16] = 2.828427124746191*nuVtSq[4]; 
  alpha_diffusion[17] = 2.828427124746191*nuVtSq[5]; 
  alpha_diffusion[31] = 2.828427124746191*nuVtSq[6]; 
  alpha_diffusion[32] = 2.828427124746191*nuVtSq[7]; 
  alpha_diffusion[51] = 2.828427124746191*nuVtSq[8]; 
  alpha_diffusion[52] = 2.828427124746191*nuVtSq[9]; 
  alpha_mid += std::abs((0.0441941738241592*alpha_diffusion[0]-0.04941058844013091*(alpha_diffusion[17]+alpha_diffusion[16]))*rdvSq40); 
  alpha_mid += std::abs((0.0441941738241592*alpha_diffusion[0]-0.04941058844013091*(alpha_diffusion[17]+alpha_diffusion[16]))*rdvSq41); 
  alpha_mid += std::abs((0.0441941738241592*alpha_diffusion[0]-0.04941058844013091*(alpha_diffusion[17]+alpha_diffusion[16]))*rdvSq42); 
  out[3] += 0.3061862178478971*(alpha_drag[52]*f[52]+alpha_drag[51]*f[51]+alpha_drag[32]*f[32]+alpha_drag[31]*f[31]+alpha_drag[17]*f[17]+alpha_drag[16]*f[16]+alpha_drag[6]*f[6]+alpha_drag[3]*f[3]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[52]*alpha_drag[108]+f[51]*alpha_drag[107]+f[32]*alpha_drag[88]+f[31]*alpha_drag[87]+f[17]*alpha_drag[73]+f[16]*alpha_drag[72]+f[6]*alpha_drag[62]+f[4]*alpha_drag[60]+f[2]*alpha_drag[58]+f[1]*alpha_drag[57]+f[0]*alpha_drag[56]); 
  out[5] += 0.3061862178478971*(f[52]*alpha_drag[164]+f[51]*alpha_drag[163]+f[32]*alpha_drag[144]+f[31]*alpha_drag[143]+f[17]*alpha_drag[129]+f[16]*alpha_drag[128]+f[6]*alpha_drag[118]+f[5]*alpha_drag[117]+f[2]*alpha_drag[114]+f[1]*alpha_drag[113]+f[0]*alpha_drag[112]); 
  out[7] += 0.2689264371002384*(alpha_drag[16]*f[51]+f[16]*alpha_drag[51])+0.3061862178478971*(alpha_drag[17]*f[32]+f[17]*alpha_drag[32])+0.273861278752583*(alpha_drag[6]*f[31]+f[6]*alpha_drag[31]+alpha_drag[1]*f[16]+f[1]*alpha_drag[16])+0.3061862178478971*(alpha_drag[3]*f[7]+alpha_drag[2]*f[6]+f[2]*alpha_drag[6]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[8] += 0.2689264371002384*(alpha_drag[17]*f[52]+f[17]*alpha_drag[52])+0.273861278752583*(alpha_drag[6]*f[32]+f[6]*alpha_drag[32])+0.3061862178478971*(alpha_drag[16]*f[31]+f[16]*alpha_drag[31])+0.273861278752583*(alpha_drag[2]*f[17]+f[2]*alpha_drag[17])+0.3061862178478971*(alpha_drag[3]*f[8]+alpha_drag[1]*f[6]+f[1]*alpha_drag[6]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[9] += 0.2689264371002384*f[16]*alpha_drag[107]+0.3061862178478971*f[17]*alpha_drag[88]+0.273861278752583*f[6]*alpha_drag[87]+0.3061862178478971*f[32]*alpha_drag[73]+0.2689264371002384*f[51]*alpha_drag[72]+0.273861278752583*(f[1]*alpha_drag[72]+f[31]*alpha_drag[62])+0.3061862178478971*(f[2]*alpha_drag[62]+f[9]*alpha_drag[60]+f[6]*alpha_drag[58])+0.273861278752583*f[16]*alpha_drag[57]+0.3061862178478971*(f[0]*alpha_drag[57]+f[1]*alpha_drag[56]); 
  out[10] += 0.2689264371002384*f[17]*alpha_drag[108]+0.273861278752583*f[6]*alpha_drag[88]+0.3061862178478971*f[16]*alpha_drag[87]+(0.2689264371002384*f[52]+0.273861278752583*f[2])*alpha_drag[73]+0.3061862178478971*f[31]*alpha_drag[72]+0.273861278752583*f[32]*alpha_drag[62]+0.3061862178478971*(f[1]*alpha_drag[62]+f[10]*alpha_drag[60])+0.273861278752583*f[17]*alpha_drag[58]+0.3061862178478971*(f[0]*alpha_drag[58]+f[6]*alpha_drag[57]+f[2]*alpha_drag[56]); 
  out[11] += 0.3061862178478971*(f[34]*alpha_drag[73]+f[33]*alpha_drag[72]+f[21]*alpha_drag[62]+f[11]*alpha_drag[60]+f[8]*alpha_drag[58]+f[7]*alpha_drag[57]+f[3]*alpha_drag[56]+alpha_drag[17]*f[38]+alpha_drag[16]*f[37]+alpha_drag[6]*f[22]+alpha_drag[3]*f[11]+alpha_drag[2]*f[10]+alpha_drag[1]*f[9]+alpha_drag[0]*f[4]); 
  out[12] += 0.2689264371002384*f[16]*alpha_drag[163]+0.3061862178478971*f[17]*alpha_drag[144]+0.273861278752583*f[6]*alpha_drag[143]+0.3061862178478971*f[32]*alpha_drag[129]+0.2689264371002384*f[51]*alpha_drag[128]+0.273861278752583*(f[1]*alpha_drag[128]+f[31]*alpha_drag[118])+0.3061862178478971*(f[2]*alpha_drag[118]+f[12]*alpha_drag[117]+f[6]*alpha_drag[114])+0.273861278752583*f[16]*alpha_drag[113]+0.3061862178478971*(f[0]*alpha_drag[113]+f[1]*alpha_drag[112]); 
  out[13] += 0.2689264371002384*f[17]*alpha_drag[164]+0.273861278752583*f[6]*alpha_drag[144]+0.3061862178478971*f[16]*alpha_drag[143]+(0.2689264371002384*f[52]+0.273861278752583*f[2])*alpha_drag[129]+0.3061862178478971*f[31]*alpha_drag[128]+0.273861278752583*f[32]*alpha_drag[118]+0.3061862178478971*(f[1]*alpha_drag[118]+f[13]*alpha_drag[117])+0.273861278752583*f[17]*alpha_drag[114]+0.3061862178478971*(f[0]*alpha_drag[114]+f[6]*alpha_drag[113]+f[2]*alpha_drag[112]); 
  out[14] += 0.3061862178478971*(f[34]*alpha_drag[129]+f[33]*alpha_drag[128]+f[21]*alpha_drag[118]+f[14]*alpha_drag[117]+f[8]*alpha_drag[114]+f[7]*alpha_drag[113]+f[3]*alpha_drag[112]+alpha_drag[17]*f[44]+alpha_drag[16]*f[43]+alpha_drag[6]*f[25]+alpha_drag[3]*f[14]+alpha_drag[2]*f[13]+alpha_drag[1]*f[12]+alpha_drag[0]*f[5]); 
  out[15] += 0.3061862178478971*(f[38]*alpha_drag[129]+f[37]*alpha_drag[128]+f[22]*alpha_drag[118]+f[15]*alpha_drag[117]+f[10]*alpha_drag[114]+f[9]*alpha_drag[113]+f[4]*alpha_drag[112]+f[44]*alpha_drag[73]+f[43]*alpha_drag[72]+f[25]*alpha_drag[62]+f[15]*alpha_drag[60]+f[13]*alpha_drag[58]+f[12]*alpha_drag[57]+f[5]*alpha_drag[56]); 
  out[18] += 1.185854122563142*(alpha_diffusion[52]*f[52]+alpha_diffusion[51]*f[51]+alpha_diffusion[32]*f[32]+alpha_diffusion[31]*f[31]+alpha_diffusion[17]*f[17]+alpha_diffusion[16]*f[16]+alpha_diffusion[6]*f[6]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq40+0.6846531968814573*(alpha_drag[17]*f[34]+alpha_drag[16]*f[33]+alpha_drag[6]*f[21])+0.6123724356957944*alpha_drag[3]*f[18]+0.6846531968814573*(alpha_drag[2]*f[8]+alpha_drag[1]*f[7]+alpha_drag[0]*f[3]+f[0]*alpha_drag[3]); 
  out[19] += 1.185854122563142*(alpha_diffusion[52]*f[52]+alpha_diffusion[51]*f[51]+alpha_diffusion[32]*f[32]+alpha_diffusion[31]*f[31]+alpha_diffusion[17]*f[17]+alpha_diffusion[16]*f[16]+alpha_diffusion[6]*f[6]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq41+0.6846531968814573*(f[38]*alpha_drag[73]+f[37]*alpha_drag[72]+f[22]*alpha_drag[62])+0.6123724356957944*f[19]*alpha_drag[60]+0.6846531968814573*(f[0]*alpha_drag[60]+f[10]*alpha_drag[58]+f[9]*alpha_drag[57]+f[4]*alpha_drag[56]); 
  out[20] += 1.185854122563142*(alpha_diffusion[52]*f[52]+alpha_diffusion[51]*f[51]+alpha_diffusion[32]*f[32]+alpha_diffusion[31]*f[31]+alpha_diffusion[17]*f[17]+alpha_diffusion[16]*f[16]+alpha_diffusion[6]*f[6]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq42+0.6846531968814573*(f[44]*alpha_drag[129]+f[43]*alpha_drag[128]+f[25]*alpha_drag[118])+0.6123724356957944*f[20]*alpha_drag[117]+0.6846531968814573*(f[0]*alpha_drag[117]+f[13]*alpha_drag[114]+f[12]*alpha_drag[113]+f[5]*alpha_drag[112]); 
  out[21] += 0.2689264371002384*(alpha_drag[32]*f[52]+f[32]*alpha_drag[52]+alpha_drag[31]*f[51]+f[31]*alpha_drag[51])+(0.2449489742783178*alpha_drag[31]+0.273861278752583*alpha_drag[2])*f[32]+0.2449489742783178*f[31]*alpha_drag[32]+0.273861278752583*(f[2]*alpha_drag[32]+alpha_drag[1]*f[31]+f[1]*alpha_drag[31])+0.3061862178478971*alpha_drag[3]*f[21]+0.273861278752583*(alpha_drag[6]*f[17]+f[6]*alpha_drag[17]+alpha_drag[6]*f[16]+f[6]*alpha_drag[16])+0.3061862178478971*(alpha_drag[0]*f[6]+f[0]*alpha_drag[6]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 
  out[22] += 0.2689264371002384*(f[32]*alpha_drag[108]+f[31]*alpha_drag[107])+(0.2689264371002384*f[52]+0.2449489742783178*f[31]+0.273861278752583*f[2])*alpha_drag[88]+(0.2689264371002384*f[51]+0.2449489742783178*f[32])*alpha_drag[87]+0.273861278752583*(f[1]*alpha_drag[87]+f[6]*(alpha_drag[73]+alpha_drag[72])+(f[17]+f[16])*alpha_drag[62])+0.3061862178478971*(f[0]*alpha_drag[62]+f[22]*alpha_drag[60])+(0.273861278752583*f[32]+0.3061862178478971*f[1])*alpha_drag[58]+0.273861278752583*f[31]*alpha_drag[57]+0.3061862178478971*(f[2]*alpha_drag[57]+f[6]*alpha_drag[56]); 
  out[23] += 0.2689264371002384*f[33]*alpha_drag[107]+0.3061862178478971*f[34]*alpha_drag[88]+0.273861278752583*(f[21]*alpha_drag[87]+f[7]*alpha_drag[72])+0.3061862178478971*(f[8]*alpha_drag[62]+f[23]*alpha_drag[60]+f[21]*alpha_drag[58])+0.273861278752583*f[33]*alpha_drag[57]+0.3061862178478971*(f[3]*alpha_drag[57]+f[7]*alpha_drag[56])+0.2689264371002384*f[37]*alpha_drag[51]+0.3061862178478971*alpha_drag[32]*f[38]+0.273861278752583*(alpha_drag[1]*f[37]+f[22]*alpha_drag[31])+0.3061862178478971*(alpha_drag[3]*f[23]+alpha_drag[2]*f[22])+0.273861278752583*f[9]*alpha_drag[16]+0.3061862178478971*(alpha_drag[6]*f[10]+alpha_drag[0]*f[9]+alpha_drag[1]*f[4]); 
  out[24] += 0.2689264371002384*f[34]*alpha_drag[108]+0.273861278752583*f[21]*alpha_drag[88]+0.3061862178478971*f[33]*alpha_drag[87]+0.273861278752583*f[8]*alpha_drag[73]+0.3061862178478971*(f[7]*alpha_drag[62]+f[24]*alpha_drag[60])+0.273861278752583*f[34]*alpha_drag[58]+0.3061862178478971*(f[3]*alpha_drag[58]+f[21]*alpha_drag[57]+f[8]*alpha_drag[56])+f[38]*(0.2689264371002384*alpha_drag[52]+0.273861278752583*alpha_drag[2])+0.3061862178478971*alpha_drag[31]*f[37]+0.273861278752583*f[22]*alpha_drag[32]+0.3061862178478971*(alpha_drag[3]*f[24]+alpha_drag[1]*f[22])+0.273861278752583*f[10]*alpha_drag[17]+0.3061862178478971*(alpha_drag[0]*f[10]+alpha_drag[6]*f[9]+alpha_drag[2]*f[4]); 
  out[25] += 0.2689264371002384*(f[32]*alpha_drag[164]+f[31]*alpha_drag[163])+(0.2689264371002384*f[52]+0.2449489742783178*f[31]+0.273861278752583*f[2])*alpha_drag[144]+(0.2689264371002384*f[51]+0.2449489742783178*f[32])*alpha_drag[143]+0.273861278752583*(f[1]*alpha_drag[143]+f[6]*(alpha_drag[129]+alpha_drag[128])+(f[17]+f[16])*alpha_drag[118])+0.3061862178478971*(f[0]*alpha_drag[118]+f[25]*alpha_drag[117])+(0.273861278752583*f[32]+0.3061862178478971*f[1])*alpha_drag[114]+0.273861278752583*f[31]*alpha_drag[113]+0.3061862178478971*(f[2]*alpha_drag[113]+f[6]*alpha_drag[112]); 
  out[26] += 0.2689264371002384*f[33]*alpha_drag[163]+0.3061862178478971*f[34]*alpha_drag[144]+0.273861278752583*(f[21]*alpha_drag[143]+f[7]*alpha_drag[128])+0.3061862178478971*(f[8]*alpha_drag[118]+f[26]*alpha_drag[117]+f[21]*alpha_drag[114])+0.273861278752583*f[33]*alpha_drag[113]+0.3061862178478971*(f[3]*alpha_drag[113]+f[7]*alpha_drag[112])+0.2689264371002384*f[43]*alpha_drag[51]+0.3061862178478971*alpha_drag[32]*f[44]+0.273861278752583*(alpha_drag[1]*f[43]+f[25]*alpha_drag[31])+0.3061862178478971*(alpha_drag[3]*f[26]+alpha_drag[2]*f[25])+0.273861278752583*f[12]*alpha_drag[16]+0.3061862178478971*(alpha_drag[6]*f[13]+alpha_drag[0]*f[12]+alpha_drag[1]*f[5]); 
  out[27] += 0.2689264371002384*f[34]*alpha_drag[164]+0.273861278752583*f[21]*alpha_drag[144]+0.3061862178478971*f[33]*alpha_drag[143]+0.273861278752583*f[8]*alpha_drag[129]+0.3061862178478971*(f[7]*alpha_drag[118]+f[27]*alpha_drag[117])+0.273861278752583*f[34]*alpha_drag[114]+0.3061862178478971*(f[3]*alpha_drag[114]+f[21]*alpha_drag[113]+f[8]*alpha_drag[112])+f[44]*(0.2689264371002384*alpha_drag[52]+0.273861278752583*alpha_drag[2])+0.3061862178478971*alpha_drag[31]*f[43]+0.273861278752583*f[25]*alpha_drag[32]+0.3061862178478971*(alpha_drag[3]*f[27]+alpha_drag[1]*f[25])+0.273861278752583*f[13]*alpha_drag[17]+0.3061862178478971*(alpha_drag[0]*f[13]+alpha_drag[6]*f[12]+alpha_drag[2]*f[5]); 
  out[28] += 0.2689264371002384*f[37]*alpha_drag[163]+0.3061862178478971*f[38]*alpha_drag[144]+0.273861278752583*(f[22]*alpha_drag[143]+f[9]*alpha_drag[128])+0.3061862178478971*(f[10]*alpha_drag[118]+f[28]*alpha_drag[117]+f[22]*alpha_drag[114])+0.273861278752583*f[37]*alpha_drag[113]+0.3061862178478971*(f[4]*alpha_drag[113]+f[9]*alpha_drag[112])+0.2689264371002384*f[43]*alpha_drag[107]+0.3061862178478971*f[44]*alpha_drag[88]+0.273861278752583*(f[25]*alpha_drag[87]+f[12]*alpha_drag[72])+0.3061862178478971*(f[13]*alpha_drag[62]+f[28]*alpha_drag[60]+f[25]*alpha_drag[58])+0.273861278752583*f[43]*alpha_drag[57]+0.3061862178478971*(f[5]*alpha_drag[57]+f[12]*alpha_drag[56]); 
  out[29] += 0.2689264371002384*f[38]*alpha_drag[164]+0.273861278752583*f[22]*alpha_drag[144]+0.3061862178478971*f[37]*alpha_drag[143]+0.273861278752583*f[10]*alpha_drag[129]+0.3061862178478971*(f[9]*alpha_drag[118]+f[29]*alpha_drag[117])+0.273861278752583*f[38]*alpha_drag[114]+0.3061862178478971*(f[4]*alpha_drag[114]+f[22]*alpha_drag[113]+f[10]*alpha_drag[112])+0.2689264371002384*f[44]*alpha_drag[108]+0.273861278752583*f[25]*alpha_drag[88]+0.3061862178478971*f[43]*alpha_drag[87]+0.273861278752583*f[13]*alpha_drag[73]+0.3061862178478971*(f[12]*alpha_drag[62]+f[29]*alpha_drag[60])+0.273861278752583*f[44]*alpha_drag[58]+0.3061862178478971*(f[5]*alpha_drag[58]+f[25]*alpha_drag[57]+f[13]*alpha_drag[56]); 
  out[30] += 0.3061862178478971*(f[30]*alpha_drag[117]+f[24]*alpha_drag[114]+f[23]*alpha_drag[113]+f[11]*alpha_drag[112]+f[30]*alpha_drag[60]+f[27]*alpha_drag[58]+f[26]*alpha_drag[57]+f[14]*alpha_drag[56]+alpha_drag[3]*f[30]+alpha_drag[2]*f[29]+alpha_drag[1]*f[28]+alpha_drag[0]*f[15]); 
  out[33] += 0.1825741858350554*alpha_drag[51]*f[51]+0.2689264371002384*(alpha_drag[1]*f[51]+f[1]*alpha_drag[51])+0.3061862178478971*alpha_drag[3]*f[33]+0.273861278752583*alpha_drag[32]*f[32]+0.1956151991089878*alpha_drag[31]*f[31]+0.3061862178478971*(alpha_drag[2]*f[31]+f[2]*alpha_drag[31])+0.1956151991089878*alpha_drag[16]*f[16]+0.3061862178478971*(alpha_drag[0]*f[16]+f[0]*alpha_drag[16])+0.273861278752583*(alpha_drag[6]*f[6]+alpha_drag[1]*f[1]); 
  out[34] += 0.1825741858350554*alpha_drag[52]*f[52]+0.2689264371002384*(alpha_drag[2]*f[52]+f[2]*alpha_drag[52])+0.3061862178478971*alpha_drag[3]*f[34]+0.1956151991089878*alpha_drag[32]*f[32]+0.3061862178478971*(alpha_drag[1]*f[32]+f[1]*alpha_drag[32])+0.273861278752583*alpha_drag[31]*f[31]+0.1956151991089878*alpha_drag[17]*f[17]+0.3061862178478971*(alpha_drag[0]*f[17]+f[0]*alpha_drag[17])+0.273861278752583*(alpha_drag[6]*f[6]+alpha_drag[2]*f[2]); 
  out[35] += (1.04154761224412*(alpha_diffusion[16]*f[51]+f[16]*alpha_diffusion[51])+1.185854122563142*(alpha_diffusion[17]*f[32]+f[17]*alpha_diffusion[32])+1.060660171779821*(alpha_diffusion[6]*f[31]+f[6]*alpha_diffusion[31]+alpha_diffusion[1]*f[16]+f[1]*alpha_diffusion[16])+1.185854122563142*(alpha_diffusion[2]*f[6]+f[2]*alpha_diffusion[6]+alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvSq40+0.6013377943029545*f[33]*alpha_drag[51]+0.6123724356957944*alpha_drag[3]*f[35]+0.6846531968814573*alpha_drag[32]*f[34]+0.6123724356957944*alpha_drag[1]*f[33]+f[21]*(0.6123724356957944*alpha_drag[31]+0.6846531968814573*alpha_drag[2])+0.6123724356957944*f[7]*alpha_drag[16]+0.6846531968814573*(alpha_drag[6]*f[8]+alpha_drag[0]*f[7]+alpha_drag[1]*f[3]+f[1]*alpha_drag[3]); 
  out[36] += (1.04154761224412*(alpha_diffusion[17]*f[52]+f[17]*alpha_diffusion[52])+1.060660171779821*(alpha_diffusion[6]*f[32]+f[6]*alpha_diffusion[32])+1.185854122563142*(alpha_diffusion[16]*f[31]+f[16]*alpha_diffusion[31])+1.060660171779821*(alpha_diffusion[2]*f[17]+f[2]*alpha_diffusion[17])+1.185854122563142*(alpha_diffusion[1]*f[6]+f[1]*alpha_diffusion[6]+alpha_diffusion[0]*f[2]+f[0]*alpha_diffusion[2]))*rdvSq40+0.6013377943029545*f[34]*alpha_drag[52]+0.6123724356957944*(alpha_drag[3]*f[36]+alpha_drag[2]*f[34])+0.6846531968814573*alpha_drag[31]*f[33]+f[21]*(0.6123724356957944*alpha_drag[32]+0.6846531968814573*alpha_drag[1])+0.6123724356957944*f[8]*alpha_drag[17]+0.6846531968814573*(alpha_drag[0]*f[8]+alpha_drag[6]*f[7]+alpha_drag[2]*f[3]+f[2]*alpha_drag[3]); 
  out[37] += (0.1825741858350554*f[51]+0.2689264371002384*f[1])*alpha_drag[107]+0.273861278752583*f[32]*alpha_drag[88]+(0.1956151991089878*f[31]+0.3061862178478971*f[2])*alpha_drag[87]+(0.1956151991089878*f[16]+0.3061862178478971*f[0])*alpha_drag[72]+0.273861278752583*f[6]*alpha_drag[62]+0.3061862178478971*(f[37]*alpha_drag[60]+f[31]*alpha_drag[58])+(0.2689264371002384*f[51]+0.273861278752583*f[1])*alpha_drag[57]+0.3061862178478971*f[16]*alpha_drag[56]; 
  out[38] += (0.1825741858350554*f[52]+0.2689264371002384*f[2])*alpha_drag[108]+(0.1956151991089878*f[32]+0.3061862178478971*f[1])*alpha_drag[88]+0.273861278752583*f[31]*alpha_drag[87]+(0.1956151991089878*f[17]+0.3061862178478971*f[0])*alpha_drag[73]+0.273861278752583*f[6]*alpha_drag[62]+0.3061862178478971*f[38]*alpha_drag[60]+(0.2689264371002384*f[52]+0.273861278752583*f[2])*alpha_drag[58]+0.3061862178478971*(f[32]*alpha_drag[57]+f[17]*alpha_drag[56]); 
  out[39] += 1.185854122563142*(alpha_diffusion[17]*f[38]+alpha_diffusion[16]*f[37]+alpha_diffusion[6]*f[22]+alpha_diffusion[2]*f[10]+alpha_diffusion[1]*f[9]+alpha_diffusion[0]*f[4])*rdvSq40+0.3061862178478971*(f[39]*alpha_drag[60]+f[36]*alpha_drag[58]+f[35]*alpha_drag[57]+f[18]*alpha_drag[56])+0.6123724356957944*alpha_drag[3]*f[39]+0.6846531968814573*(alpha_drag[2]*f[24]+alpha_drag[1]*f[23]+alpha_drag[0]*f[11]+alpha_drag[3]*f[4]); 
  out[40] += (1.04154761224412*(alpha_diffusion[16]*f[51]+f[16]*alpha_diffusion[51])+1.185854122563142*(alpha_diffusion[17]*f[32]+f[17]*alpha_diffusion[32])+1.060660171779821*(alpha_diffusion[6]*f[31]+f[6]*alpha_diffusion[31]+alpha_diffusion[1]*f[16]+f[1]*alpha_diffusion[16])+1.185854122563142*(alpha_diffusion[2]*f[6]+f[2]*alpha_diffusion[6]+alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvSq41+0.6013377943029545*f[37]*alpha_drag[107]+0.6846531968814573*f[38]*alpha_drag[88]+0.6123724356957944*(f[22]*alpha_drag[87]+f[9]*alpha_drag[72])+0.6846531968814573*f[10]*alpha_drag[62]+0.6123724356957944*f[40]*alpha_drag[60]+0.6846531968814573*(f[1]*alpha_drag[60]+f[22]*alpha_drag[58])+0.6123724356957944*f[37]*alpha_drag[57]+0.6846531968814573*(f[4]*alpha_drag[57]+f[9]*alpha_drag[56]); 
  out[41] += (1.04154761224412*(alpha_diffusion[17]*f[52]+f[17]*alpha_diffusion[52])+1.060660171779821*(alpha_diffusion[6]*f[32]+f[6]*alpha_diffusion[32])+1.185854122563142*(alpha_diffusion[16]*f[31]+f[16]*alpha_diffusion[31])+1.060660171779821*(alpha_diffusion[2]*f[17]+f[2]*alpha_diffusion[17])+1.185854122563142*(alpha_diffusion[1]*f[6]+f[1]*alpha_diffusion[6]+alpha_diffusion[0]*f[2]+f[0]*alpha_diffusion[2]))*rdvSq41+0.6013377943029545*f[38]*alpha_drag[108]+0.6123724356957944*f[22]*alpha_drag[88]+0.6846531968814573*f[37]*alpha_drag[87]+0.6123724356957944*f[10]*alpha_drag[73]+0.6846531968814573*f[9]*alpha_drag[62]+(0.6123724356957944*f[41]+0.6846531968814573*f[2])*alpha_drag[60]+0.6123724356957944*f[38]*alpha_drag[58]+0.6846531968814573*(f[4]*alpha_drag[58]+f[22]*alpha_drag[57]+f[10]*alpha_drag[56]); 
  out[42] += 1.185854122563142*(alpha_diffusion[17]*f[34]+alpha_diffusion[16]*f[33]+alpha_diffusion[6]*f[21]+alpha_diffusion[2]*f[8]+alpha_diffusion[1]*f[7]+alpha_diffusion[0]*f[3])*rdvSq41+0.6123724356957944*f[42]*alpha_drag[60]+0.6846531968814573*(f[3]*alpha_drag[60]+f[24]*alpha_drag[58]+f[23]*alpha_drag[57]+f[11]*alpha_drag[56])+0.3061862178478971*(alpha_drag[3]*f[42]+alpha_drag[2]*f[41]+alpha_drag[1]*f[40]+alpha_drag[0]*f[19]); 
  out[43] += (0.1825741858350554*f[51]+0.2689264371002384*f[1])*alpha_drag[163]+0.273861278752583*f[32]*alpha_drag[144]+(0.1956151991089878*f[31]+0.3061862178478971*f[2])*alpha_drag[143]+(0.1956151991089878*f[16]+0.3061862178478971*f[0])*alpha_drag[128]+0.273861278752583*f[6]*alpha_drag[118]+0.3061862178478971*(f[43]*alpha_drag[117]+f[31]*alpha_drag[114])+(0.2689264371002384*f[51]+0.273861278752583*f[1])*alpha_drag[113]+0.3061862178478971*f[16]*alpha_drag[112]; 
  out[44] += (0.1825741858350554*f[52]+0.2689264371002384*f[2])*alpha_drag[164]+(0.1956151991089878*f[32]+0.3061862178478971*f[1])*alpha_drag[144]+0.273861278752583*f[31]*alpha_drag[143]+(0.1956151991089878*f[17]+0.3061862178478971*f[0])*alpha_drag[129]+0.273861278752583*f[6]*alpha_drag[118]+0.3061862178478971*f[44]*alpha_drag[117]+(0.2689264371002384*f[52]+0.273861278752583*f[2])*alpha_drag[114]+0.3061862178478971*(f[32]*alpha_drag[113]+f[17]*alpha_drag[112]); 
  out[45] += 1.185854122563142*(alpha_diffusion[17]*f[44]+alpha_diffusion[16]*f[43]+alpha_diffusion[6]*f[25]+alpha_diffusion[2]*f[13]+alpha_diffusion[1]*f[12]+alpha_diffusion[0]*f[5])*rdvSq40+0.3061862178478971*(f[45]*alpha_drag[117]+f[36]*alpha_drag[114]+f[35]*alpha_drag[113]+f[18]*alpha_drag[112])+0.6123724356957944*alpha_drag[3]*f[45]+0.6846531968814573*(alpha_drag[2]*f[27]+alpha_drag[1]*f[26]+alpha_drag[0]*f[14]+alpha_drag[3]*f[5]); 
  out[46] += 1.185854122563142*(alpha_diffusion[17]*f[44]+alpha_diffusion[16]*f[43]+alpha_diffusion[6]*f[25]+alpha_diffusion[2]*f[13]+alpha_diffusion[1]*f[12]+alpha_diffusion[0]*f[5])*rdvSq41+0.3061862178478971*(f[46]*alpha_drag[117]+f[41]*alpha_drag[114]+f[40]*alpha_drag[113]+f[19]*alpha_drag[112])+0.6123724356957944*f[46]*alpha_drag[60]+0.6846531968814573*(f[5]*alpha_drag[60]+f[29]*alpha_drag[58]+f[28]*alpha_drag[57]+f[15]*alpha_drag[56]); 
  out[47] += (1.04154761224412*(alpha_diffusion[16]*f[51]+f[16]*alpha_diffusion[51])+1.185854122563142*(alpha_diffusion[17]*f[32]+f[17]*alpha_diffusion[32])+1.060660171779821*(alpha_diffusion[6]*f[31]+f[6]*alpha_diffusion[31]+alpha_diffusion[1]*f[16]+f[1]*alpha_diffusion[16])+1.185854122563142*(alpha_diffusion[2]*f[6]+f[2]*alpha_diffusion[6]+alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvSq42+0.6013377943029545*f[43]*alpha_drag[163]+0.6846531968814573*f[44]*alpha_drag[144]+0.6123724356957944*(f[25]*alpha_drag[143]+f[12]*alpha_drag[128])+0.6846531968814573*f[13]*alpha_drag[118]+0.6123724356957944*f[47]*alpha_drag[117]+0.6846531968814573*(f[1]*alpha_drag[117]+f[25]*alpha_drag[114])+0.6123724356957944*f[43]*alpha_drag[113]+0.6846531968814573*(f[5]*alpha_drag[113]+f[12]*alpha_drag[112]); 
  out[48] += (1.04154761224412*(alpha_diffusion[17]*f[52]+f[17]*alpha_diffusion[52])+1.060660171779821*(alpha_diffusion[6]*f[32]+f[6]*alpha_diffusion[32])+1.185854122563142*(alpha_diffusion[16]*f[31]+f[16]*alpha_diffusion[31])+1.060660171779821*(alpha_diffusion[2]*f[17]+f[2]*alpha_diffusion[17])+1.185854122563142*(alpha_diffusion[1]*f[6]+f[1]*alpha_diffusion[6]+alpha_diffusion[0]*f[2]+f[0]*alpha_diffusion[2]))*rdvSq42+0.6013377943029545*f[44]*alpha_drag[164]+0.6123724356957944*f[25]*alpha_drag[144]+0.6846531968814573*f[43]*alpha_drag[143]+0.6123724356957944*f[13]*alpha_drag[129]+0.6846531968814573*f[12]*alpha_drag[118]+(0.6123724356957944*f[48]+0.6846531968814573*f[2])*alpha_drag[117]+0.6123724356957944*f[44]*alpha_drag[114]+0.6846531968814573*(f[5]*alpha_drag[114]+f[25]*alpha_drag[113]+f[13]*alpha_drag[112]); 
  out[49] += 1.185854122563142*(alpha_diffusion[17]*f[34]+alpha_diffusion[16]*f[33]+alpha_diffusion[6]*f[21]+alpha_diffusion[2]*f[8]+alpha_diffusion[1]*f[7]+alpha_diffusion[0]*f[3])*rdvSq42+0.6123724356957944*f[49]*alpha_drag[117]+0.6846531968814573*(f[3]*alpha_drag[117]+f[27]*alpha_drag[114]+f[26]*alpha_drag[113]+f[14]*alpha_drag[112])+0.3061862178478971*(alpha_drag[3]*f[49]+alpha_drag[2]*f[48]+alpha_drag[1]*f[47]+alpha_drag[0]*f[20]); 
  out[50] += 1.185854122563142*(alpha_diffusion[17]*f[38]+alpha_diffusion[16]*f[37]+alpha_diffusion[6]*f[22]+alpha_diffusion[2]*f[10]+alpha_diffusion[1]*f[9]+alpha_diffusion[0]*f[4])*rdvSq42+0.6123724356957944*f[50]*alpha_drag[117]+0.6846531968814573*(f[4]*alpha_drag[117]+f[29]*alpha_drag[114]+f[28]*alpha_drag[113]+f[15]*alpha_drag[112])+0.3061862178478971*(f[50]*alpha_drag[60]+f[48]*alpha_drag[58]+f[47]*alpha_drag[57]+f[20]*alpha_drag[56]); 
  out[53] += 4.050462936504911*(alpha_diffusion[17]*f[34]+alpha_diffusion[16]*f[33]+alpha_diffusion[6]*f[21]+alpha_diffusion[2]*f[8]+alpha_diffusion[1]*f[7]+alpha_diffusion[0]*f[3])*rdvSq40+0.9185586535436913*alpha_drag[3]*f[53]+0.4677071733467425*(alpha_drag[52]*f[52]+alpha_drag[51]*f[51])+1.045825033167594*(alpha_drag[2]*f[36]+alpha_drag[1]*f[35])+0.4677071733467425*(alpha_drag[32]*f[32]+alpha_drag[31]*f[31])+1.045825033167594*alpha_drag[0]*f[18]+0.4677071733467425*(alpha_drag[17]*f[17]+alpha_drag[16]*f[16]+alpha_drag[6]*f[6])+1.403121520040228*alpha_drag[3]*f[3]+0.4677071733467425*(alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[54] += 4.050462936504911*(alpha_diffusion[17]*f[38]+alpha_diffusion[16]*f[37]+alpha_diffusion[6]*f[22]+alpha_diffusion[2]*f[10]+alpha_diffusion[1]*f[9]+alpha_diffusion[0]*f[4])*rdvSq41+0.4677071733467425*(f[52]*alpha_drag[108]+f[51]*alpha_drag[107]+f[32]*alpha_drag[88]+f[31]*alpha_drag[87]+f[17]*alpha_drag[73]+f[16]*alpha_drag[72]+f[6]*alpha_drag[62])+(0.9185586535436913*f[54]+1.403121520040228*f[4])*alpha_drag[60]+(1.045825033167594*f[41]+0.4677071733467425*f[2])*alpha_drag[58]+(1.045825033167594*f[40]+0.4677071733467425*f[1])*alpha_drag[57]+(1.045825033167594*f[19]+0.4677071733467425*f[0])*alpha_drag[56]; 
  out[55] += 4.050462936504911*(alpha_diffusion[17]*f[44]+alpha_diffusion[16]*f[43]+alpha_diffusion[6]*f[25]+alpha_diffusion[2]*f[13]+alpha_diffusion[1]*f[12]+alpha_diffusion[0]*f[5])*rdvSq42+0.4677071733467425*(f[52]*alpha_drag[164]+f[51]*alpha_drag[163]+f[32]*alpha_drag[144]+f[31]*alpha_drag[143]+f[17]*alpha_drag[129]+f[16]*alpha_drag[128]+f[6]*alpha_drag[118])+(0.9185586535436913*f[55]+1.403121520040228*f[5])*alpha_drag[117]+(1.045825033167594*f[48]+0.4677071733467425*f[2])*alpha_drag[114]+(1.045825033167594*f[47]+0.4677071733467425*f[1])*alpha_drag[113]+(1.045825033167594*f[20]+0.4677071733467425*f[0])*alpha_drag[112]; 
return alpha_mid; 

} 
