#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x3vMaxP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdv20 = 2/dxv[1]; 
  const double rdvSq40 = 4/(dxv[1]*dxv[1]); 
  const double rdv21 = 2/dxv[2]; 
  const double rdvSq41 = 4/(dxv[2]*dxv[2]); 
  const double rdv22 = 2/dxv[3]; 
  const double rdvSq42 = 4/(dxv[3]*dxv[3]); 

  double alpha_mid = 0.0; 
  double alpha_drag[15]; 
  double alpha_diffusion[5]; 

  alpha_drag[0] = (2.828427124746191*nuU[0]-4.0*w[1]*nu)*rdv20; 
  alpha_drag[1] = 2.828427124746191*nuU[1]*rdv20; 
  alpha_drag[2] = -1.154700538379252*dxv[1]*nu*rdv20; 
  alpha_mid += std::abs(0.125*alpha_drag[0]); 
  alpha_drag[5] = (2.828427124746191*nuU[5]-4.0*w[2]*nu)*rdv21; 
  alpha_drag[6] = 2.828427124746191*nuU[6]*rdv21; 
  alpha_drag[8] = -1.154700538379252*dxv[2]*nu*rdv21; 
  alpha_mid += std::abs(0.125*alpha_drag[5]); 
  alpha_drag[10] = (2.828427124746191*nuU[10]-4.0*w[3]*nu)*rdv22; 
  alpha_drag[11] = 2.828427124746191*nuU[11]*rdv22; 
  alpha_drag[14] = -1.154700538379252*dxv[3]*nu*rdv22; 
  alpha_mid += std::abs(0.125*alpha_drag[10]); 
  alpha_diffusion[0] = 2.828427124746191*nuVtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*nuVtSq[1]; 
  alpha_mid += std::abs(0.0625*alpha_diffusion[0]*rdvSq40); 
  alpha_mid += std::abs(0.0625*alpha_diffusion[0]*rdvSq41); 
  alpha_mid += std::abs(0.0625*alpha_diffusion[0]*rdvSq42); 
  out[2] += 0.4330127018922193*(alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[3]*alpha_drag[8]+f[1]*alpha_drag[6]+f[0]*alpha_drag[5]); 
  out[4] += 0.4330127018922193*(f[4]*alpha_drag[14]+f[1]*alpha_drag[11]+f[0]*alpha_drag[10]); 
return alpha_mid; 

} 
double VmLBOconstNuVol1x3vMaxP2(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdv20 = 2/dxv[1]; 
  const double rdvSq40 = 4/(dxv[1]*dxv[1]); 
  const double rdv21 = 2/dxv[2]; 
  const double rdvSq41 = 4/(dxv[2]*dxv[2]); 
  const double rdv22 = 2/dxv[3]; 
  const double rdvSq42 = 4/(dxv[3]*dxv[3]); 

  double alpha_mid = 0.0; 
  double alpha_drag[45]; 
  double alpha_diffusion[15]; 

  alpha_drag[0] = (2.828427124746191*nuU[0]-4.0*w[1]*nu)*rdv20; 
  alpha_drag[1] = 2.828427124746191*nuU[1]*rdv20; 
  alpha_drag[2] = -1.154700538379252*dxv[1]*nu*rdv20; 
  alpha_drag[11] = 2.828427124746191*nuU[2]*rdv20; 
  alpha_mid += std::abs(0.125*alpha_drag[0]-0.1397542485937369*alpha_drag[11]); 
  alpha_drag[15] = (2.828427124746191*nuU[15]-4.0*w[2]*nu)*rdv21; 
  alpha_drag[16] = 2.828427124746191*nuU[16]*rdv21; 
  alpha_drag[18] = -1.154700538379252*dxv[2]*nu*rdv21; 
  alpha_drag[26] = 2.828427124746191*nuU[17]*rdv21; 
  alpha_mid += std::abs(0.125*alpha_drag[15]-0.1397542485937369*alpha_drag[26]); 
  alpha_drag[30] = (2.828427124746191*nuU[30]-4.0*w[3]*nu)*rdv22; 
  alpha_drag[31] = 2.828427124746191*nuU[31]*rdv22; 
  alpha_drag[34] = -1.154700538379252*dxv[3]*nu*rdv22; 
  alpha_drag[41] = 2.828427124746191*nuU[32]*rdv22; 
  alpha_mid += std::abs(0.125*alpha_drag[30]-0.1397542485937369*alpha_drag[41]); 
  alpha_diffusion[0] = 2.828427124746191*nuVtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*nuVtSq[1]; 
  alpha_diffusion[11] = 2.828427124746191*nuVtSq[2]; 
  alpha_mid += std::abs((0.0625*alpha_diffusion[0]-0.06987712429686843*alpha_diffusion[11])*rdvSq40); 
  alpha_mid += std::abs((0.0625*alpha_diffusion[0]-0.06987712429686843*alpha_diffusion[11])*rdvSq41); 
  alpha_mid += std::abs((0.0625*alpha_diffusion[0]-0.06987712429686843*alpha_diffusion[11])*rdvSq42); 
  out[2] += 0.4330127018922193*(alpha_drag[11]*f[11]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[11]*alpha_drag[26]+f[3]*alpha_drag[18]+f[1]*alpha_drag[16]+f[0]*alpha_drag[15]); 
  out[4] += 0.4330127018922193*(f[11]*alpha_drag[41]+f[4]*alpha_drag[34]+f[1]*alpha_drag[31]+f[0]*alpha_drag[30]); 
  out[5] += 0.3872983346207416*(alpha_drag[1]*f[11]+f[1]*alpha_drag[11])+0.4330127018922193*(alpha_drag[2]*f[5]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[6] += 0.3872983346207416*f[1]*alpha_drag[26]+0.4330127018922193*f[6]*alpha_drag[18]+0.3872983346207416*f[11]*alpha_drag[16]+0.4330127018922193*(f[0]*alpha_drag[16]+f[1]*alpha_drag[15]); 
  out[7] += 0.4330127018922193*(f[7]*alpha_drag[18]+f[5]*alpha_drag[16]+f[2]*alpha_drag[15]+alpha_drag[2]*f[7]+alpha_drag[1]*f[6]+alpha_drag[0]*f[3]); 
  out[8] += 0.3872983346207416*f[1]*alpha_drag[41]+0.4330127018922193*f[8]*alpha_drag[34]+0.3872983346207416*f[11]*alpha_drag[31]+0.4330127018922193*(f[0]*alpha_drag[31]+f[1]*alpha_drag[30]); 
  out[9] += 0.4330127018922193*(f[9]*alpha_drag[34]+f[5]*alpha_drag[31]+f[2]*alpha_drag[30]+alpha_drag[2]*f[9]+alpha_drag[1]*f[8]+alpha_drag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[10]*alpha_drag[34]+f[6]*alpha_drag[31]+f[3]*alpha_drag[30]+f[10]*alpha_drag[18]+f[8]*alpha_drag[16]+f[4]*alpha_drag[15]); 
  out[12] += 1.677050983124842*(alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq40+0.8660254037844386*alpha_drag[2]*f[12]+0.9682458365518543*(alpha_drag[1]*f[5]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[13] += 1.677050983124842*(alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq41+0.8660254037844386*f[13]*alpha_drag[18]+0.9682458365518543*(f[0]*alpha_drag[18]+f[6]*alpha_drag[16]+f[3]*alpha_drag[15]); 
  out[14] += 1.677050983124842*(alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq42+0.8660254037844386*f[14]*alpha_drag[34]+0.9682458365518543*(f[0]*alpha_drag[34]+f[8]*alpha_drag[31]+f[4]*alpha_drag[30]); 
return alpha_mid; 

} 
double VmLBOconstNuVol1x3vMaxP3(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdv20 = 2/dxv[1]; 
  const double rdvSq40 = 4/(dxv[1]*dxv[1]); 
  const double rdv21 = 2/dxv[2]; 
  const double rdvSq41 = 4/(dxv[2]*dxv[2]); 
  const double rdv22 = 2/dxv[3]; 
  const double rdvSq42 = 4/(dxv[3]*dxv[3]); 

  double alpha_mid = 0.0; 
  double alpha_drag[105]; 
  double alpha_diffusion[35]; 

  alpha_drag[0] = (2.828427124746191*nuU[0]-4.0*w[1]*nu)*rdv20; 
  alpha_drag[1] = 2.828427124746191*nuU[1]*rdv20; 
  alpha_drag[2] = -1.154700538379252*dxv[1]*nu*rdv20; 
  alpha_drag[11] = 2.828427124746191*nuU[2]*rdv20; 
  alpha_drag[31] = 2.828427124746191*nuU[3]*rdv20; 
  alpha_mid += std::abs(0.125*alpha_drag[0]-0.1397542485937369*alpha_drag[11]); 
  alpha_drag[35] = (2.828427124746191*nuU[35]-4.0*w[2]*nu)*rdv21; 
  alpha_drag[36] = 2.828427124746191*nuU[36]*rdv21; 
  alpha_drag[38] = -1.154700538379252*dxv[2]*nu*rdv21; 
  alpha_drag[46] = 2.828427124746191*nuU[37]*rdv21; 
  alpha_drag[66] = 2.828427124746191*nuU[38]*rdv21; 
  alpha_mid += std::abs(0.125*alpha_drag[35]-0.1397542485937369*alpha_drag[46]); 
  alpha_drag[70] = (2.828427124746191*nuU[70]-4.0*w[3]*nu)*rdv22; 
  alpha_drag[71] = 2.828427124746191*nuU[71]*rdv22; 
  alpha_drag[74] = -1.154700538379252*dxv[3]*nu*rdv22; 
  alpha_drag[81] = 2.828427124746191*nuU[72]*rdv22; 
  alpha_drag[101] = 2.828427124746191*nuU[73]*rdv22; 
  alpha_mid += std::abs(0.125*alpha_drag[70]-0.1397542485937369*alpha_drag[81]); 
  alpha_diffusion[0] = 2.828427124746191*nuVtSq[0]; 
  alpha_diffusion[1] = 2.828427124746191*nuVtSq[1]; 
  alpha_diffusion[11] = 2.828427124746191*nuVtSq[2]; 
  alpha_diffusion[31] = 2.828427124746191*nuVtSq[3]; 
  alpha_mid += std::abs((0.0625*alpha_diffusion[0]-0.06987712429686843*alpha_diffusion[11])*rdvSq40); 
  alpha_mid += std::abs((0.0625*alpha_diffusion[0]-0.06987712429686843*alpha_diffusion[11])*rdvSq41); 
  alpha_mid += std::abs((0.0625*alpha_diffusion[0]-0.06987712429686843*alpha_diffusion[11])*rdvSq42); 
  out[2] += 0.4330127018922193*(alpha_drag[31]*f[31]+alpha_drag[11]*f[11]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[31]*alpha_drag[66]+f[11]*alpha_drag[46]+f[3]*alpha_drag[38]+f[1]*alpha_drag[36]+f[0]*alpha_drag[35]); 
  out[4] += 0.4330127018922193*(f[31]*alpha_drag[101]+f[11]*alpha_drag[81]+f[4]*alpha_drag[74]+f[1]*alpha_drag[71]+f[0]*alpha_drag[70]); 
  out[5] += 0.3803194146278324*(alpha_drag[11]*f[31]+f[11]*alpha_drag[31])+0.3872983346207416*(alpha_drag[1]*f[11]+f[1]*alpha_drag[11])+0.4330127018922193*(alpha_drag[2]*f[5]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[6] += 0.3803194146278324*f[11]*alpha_drag[66]+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_drag[46]+0.4330127018922193*f[6]*alpha_drag[38]+0.3872983346207416*f[11]*alpha_drag[36]+0.4330127018922193*(f[0]*alpha_drag[36]+f[1]*alpha_drag[35]); 
  out[7] += 0.4330127018922193*(f[19]*alpha_drag[46]+f[7]*alpha_drag[38]+f[5]*alpha_drag[36]+f[2]*alpha_drag[35]+alpha_drag[11]*f[21]+alpha_drag[2]*f[7]+alpha_drag[1]*f[6]+alpha_drag[0]*f[3]); 
  out[8] += 0.3803194146278324*f[11]*alpha_drag[101]+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_drag[81]+0.4330127018922193*f[8]*alpha_drag[74]+0.3872983346207416*f[11]*alpha_drag[71]+0.4330127018922193*(f[0]*alpha_drag[71]+f[1]*alpha_drag[70]); 
  out[9] += 0.4330127018922193*(f[19]*alpha_drag[81]+f[9]*alpha_drag[74]+f[5]*alpha_drag[71]+f[2]*alpha_drag[70]+alpha_drag[11]*f[25]+alpha_drag[2]*f[9]+alpha_drag[1]*f[8]+alpha_drag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[21]*alpha_drag[81]+f[10]*alpha_drag[74]+f[6]*alpha_drag[71]+f[3]*alpha_drag[70]+f[25]*alpha_drag[46]+f[10]*alpha_drag[38]+f[8]*alpha_drag[36]+f[4]*alpha_drag[35]); 
  out[12] += 1.677050983124842*(alpha_diffusion[31]*f[31]+alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq40+0.9682458365518543*alpha_drag[11]*f[19]+0.8660254037844386*alpha_drag[2]*f[12]+0.9682458365518543*(alpha_drag[1]*f[5]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[13] += 1.677050983124842*(alpha_diffusion[31]*f[31]+alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq41+0.9682458365518543*f[21]*alpha_drag[46]+0.8660254037844386*f[13]*alpha_drag[38]+0.9682458365518543*(f[0]*alpha_drag[38]+f[6]*alpha_drag[36]+f[3]*alpha_drag[35]); 
  out[14] += 1.677050983124842*(alpha_diffusion[31]*f[31]+alpha_diffusion[11]*f[11]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvSq42+0.9682458365518543*f[25]*alpha_drag[81]+0.8660254037844386*f[14]*alpha_drag[74]+0.9682458365518543*(f[0]*alpha_drag[74]+f[8]*alpha_drag[71]+f[4]*alpha_drag[70]); 
  out[15] += 0.3803194146278324*f[19]*alpha_drag[66]+0.3872983346207416*f[5]*alpha_drag[46]+0.4330127018922193*f[15]*alpha_drag[38]+0.3872983346207416*f[19]*alpha_drag[36]+0.4330127018922193*(f[2]*alpha_drag[36]+f[5]*alpha_drag[35])+f[21]*(0.3803194146278324*alpha_drag[31]+0.3872983346207416*alpha_drag[1])+0.4330127018922193*alpha_drag[2]*f[15]+0.3872983346207416*f[6]*alpha_drag[11]+0.4330127018922193*(alpha_drag[0]*f[6]+alpha_drag[1]*f[3]); 
  out[16] += 0.3803194146278324*f[19]*alpha_drag[101]+0.3872983346207416*f[5]*alpha_drag[81]+0.4330127018922193*f[16]*alpha_drag[74]+0.3872983346207416*f[19]*alpha_drag[71]+0.4330127018922193*(f[2]*alpha_drag[71]+f[5]*alpha_drag[70])+f[25]*(0.3803194146278324*alpha_drag[31]+0.3872983346207416*alpha_drag[1])+0.4330127018922193*alpha_drag[2]*f[16]+0.3872983346207416*f[8]*alpha_drag[11]+0.4330127018922193*(alpha_drag[0]*f[8]+alpha_drag[1]*f[4]); 
  out[17] += 0.3803194146278324*f[21]*alpha_drag[101]+0.3872983346207416*f[6]*alpha_drag[81]+0.4330127018922193*f[17]*alpha_drag[74]+0.3872983346207416*f[21]*alpha_drag[71]+0.4330127018922193*(f[3]*alpha_drag[71]+f[6]*alpha_drag[70])+0.3803194146278324*f[25]*alpha_drag[66]+0.3872983346207416*f[8]*alpha_drag[46]+0.4330127018922193*f[17]*alpha_drag[38]+0.3872983346207416*f[25]*alpha_drag[36]+0.4330127018922193*(f[4]*alpha_drag[36]+f[8]*alpha_drag[35]); 
  out[18] += 0.4330127018922193*(f[18]*alpha_drag[74]+f[15]*alpha_drag[71]+f[7]*alpha_drag[70]+f[18]*alpha_drag[38]+f[16]*alpha_drag[36]+f[9]*alpha_drag[35]+alpha_drag[2]*f[18]+alpha_drag[1]*f[17]+alpha_drag[0]*f[10]); 
  out[19] += 0.2581988897471612*alpha_drag[31]*f[31]+0.3803194146278324*(alpha_drag[1]*f[31]+f[1]*alpha_drag[31])+0.4330127018922193*alpha_drag[2]*f[19]+0.276641667586244*alpha_drag[11]*f[11]+0.4330127018922193*(alpha_drag[0]*f[11]+f[0]*alpha_drag[11])+0.3872983346207416*alpha_drag[1]*f[1]; 
  out[20] += (1.472970759092948*(alpha_diffusion[11]*f[31]+f[11]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvSq40+0.8504200642707612*f[19]*alpha_drag[31]+0.8660254037844386*(alpha_drag[2]*f[20]+alpha_drag[1]*f[19]+f[5]*alpha_drag[11])+0.9682458365518543*(alpha_drag[0]*f[5]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 
  out[21] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_drag[66]+0.276641667586244*f[11]*alpha_drag[46]+0.4330127018922193*(f[0]*alpha_drag[46]+f[21]*alpha_drag[38])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_drag[36]+0.4330127018922193*f[11]*alpha_drag[35]; 
  out[22] += 1.677050983124842*(alpha_diffusion[11]*f[21]+alpha_diffusion[1]*f[6]+alpha_diffusion[0]*f[3])*rdvSq40+0.4330127018922193*(f[22]*alpha_drag[38]+f[20]*alpha_drag[36]+f[12]*alpha_drag[35])+0.8660254037844386*alpha_drag[2]*f[22]+0.9682458365518543*(alpha_drag[1]*f[15]+alpha_drag[0]*f[7]+alpha_drag[2]*f[3]); 
  out[23] += (1.472970759092948*(alpha_diffusion[11]*f[31]+f[11]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvSq41+0.8504200642707612*f[21]*alpha_drag[66]+0.8660254037844386*f[6]*alpha_drag[46]+(0.8660254037844386*f[23]+0.9682458365518543*f[1])*alpha_drag[38]+0.8660254037844386*f[21]*alpha_drag[36]+0.9682458365518543*(f[3]*alpha_drag[36]+f[6]*alpha_drag[35]); 
  out[24] += 1.677050983124842*(alpha_diffusion[11]*f[19]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[2])*rdvSq41+0.8660254037844386*f[24]*alpha_drag[38]+0.9682458365518543*(f[2]*alpha_drag[38]+f[15]*alpha_drag[36]+f[7]*alpha_drag[35])+0.4330127018922193*(alpha_drag[2]*f[24]+alpha_drag[1]*f[23]+alpha_drag[0]*f[13]); 
  out[25] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_drag[101]+0.276641667586244*f[11]*alpha_drag[81]+0.4330127018922193*(f[0]*alpha_drag[81]+f[25]*alpha_drag[74])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_drag[71]+0.4330127018922193*f[11]*alpha_drag[70]; 
  out[26] += 1.677050983124842*(alpha_diffusion[11]*f[25]+alpha_diffusion[1]*f[8]+alpha_diffusion[0]*f[4])*rdvSq40+0.4330127018922193*(f[26]*alpha_drag[74]+f[20]*alpha_drag[71]+f[12]*alpha_drag[70])+0.8660254037844386*alpha_drag[2]*f[26]+0.9682458365518543*(alpha_drag[1]*f[16]+alpha_drag[0]*f[9]+alpha_drag[2]*f[4]); 
  out[27] += 1.677050983124842*(alpha_diffusion[11]*f[25]+alpha_diffusion[1]*f[8]+alpha_diffusion[0]*f[4])*rdvSq41+0.4330127018922193*(f[27]*alpha_drag[74]+f[23]*alpha_drag[71]+f[13]*alpha_drag[70])+0.8660254037844386*f[27]*alpha_drag[38]+0.9682458365518543*(f[4]*alpha_drag[38]+f[17]*alpha_drag[36]+f[10]*alpha_drag[35]); 
  out[28] += (1.472970759092948*(alpha_diffusion[11]*f[31]+f[11]*alpha_diffusion[31])+1.5*(alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvSq42+0.8504200642707612*f[25]*alpha_drag[101]+0.8660254037844386*f[8]*alpha_drag[81]+(0.8660254037844386*f[28]+0.9682458365518543*f[1])*alpha_drag[74]+0.8660254037844386*f[25]*alpha_drag[71]+0.9682458365518543*(f[4]*alpha_drag[71]+f[8]*alpha_drag[70]); 
  out[29] += 1.677050983124842*(alpha_diffusion[11]*f[19]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[2])*rdvSq42+0.8660254037844386*f[29]*alpha_drag[74]+0.9682458365518543*(f[2]*alpha_drag[74]+f[16]*alpha_drag[71]+f[9]*alpha_drag[70])+0.4330127018922193*(alpha_drag[2]*f[29]+alpha_drag[1]*f[28]+alpha_drag[0]*f[14]); 
  out[30] += 1.677050983124842*(alpha_diffusion[11]*f[21]+alpha_diffusion[1]*f[6]+alpha_diffusion[0]*f[3])*rdvSq42+0.8660254037844386*f[30]*alpha_drag[74]+0.9682458365518543*(f[3]*alpha_drag[74]+f[17]*alpha_drag[71]+f[10]*alpha_drag[70])+0.4330127018922193*(f[30]*alpha_drag[38]+f[28]*alpha_drag[36]+f[14]*alpha_drag[35]); 
  out[32] += 5.7282196186948*(alpha_diffusion[11]*f[19]+alpha_diffusion[1]*f[5]+alpha_diffusion[0]*f[2])*rdvSq40+1.299038105676658*alpha_drag[2]*f[32]+0.6614378277661477*alpha_drag[31]*f[31]+1.479019945774904*(alpha_drag[1]*f[20]+alpha_drag[0]*f[12])+0.6614378277661477*alpha_drag[11]*f[11]+1.984313483298443*alpha_drag[2]*f[2]+0.6614378277661477*(alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[33] += 5.7282196186948*(alpha_diffusion[11]*f[21]+alpha_diffusion[1]*f[6]+alpha_diffusion[0]*f[3])*rdvSq41+0.6614378277661477*(f[31]*alpha_drag[66]+f[11]*alpha_drag[46])+(1.299038105676658*f[33]+1.984313483298443*f[3])*alpha_drag[38]+(1.479019945774904*f[23]+0.6614378277661477*f[1])*alpha_drag[36]+(1.479019945774904*f[13]+0.6614378277661477*f[0])*alpha_drag[35]; 
  out[34] += 5.7282196186948*(alpha_diffusion[11]*f[25]+alpha_diffusion[1]*f[8]+alpha_diffusion[0]*f[4])*rdvSq42+0.6614378277661477*(f[31]*alpha_drag[101]+f[11]*alpha_drag[81])+(1.299038105676658*f[34]+1.984313483298443*f[4])*alpha_drag[74]+(1.479019945774904*f[28]+0.6614378277661477*f[1])*alpha_drag[71]+(1.479019945774904*f[14]+0.6614378277661477*f[0])*alpha_drag[70]; 
return alpha_mid; 

} 
