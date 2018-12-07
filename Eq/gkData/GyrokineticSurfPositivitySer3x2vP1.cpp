#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity3x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_z = 2.0/dxv[2]; 
  double dfac_v = 2.0/dxv[3]; 
  double dfac_m = 2.0/dxv[4]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wz = w[2]; 
  double wv = w[3]; 
  double wm = w[4]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[32]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.0625*(2.828427124746191*BdriftX[0]*m_*wv2+BmagInv[0]*(3.0*Phi[4]-1.732050807568877*Phi[2])*dfac_y*q_))/q_; 

  double alpha[16]; 
  alpha[0] = (0.5*(2.828427124746191*BdriftX[0]*m_*wv2+BmagInv[0]*(3.0*Phi[4]-1.732050807568877*Phi[2])*dfac_y*q_))/q_; 
  alpha[2] = 0.5*BmagInv[0]*(3.0*Phi[7]-1.732050807568877*Phi[6])*dfac_y; 
  alpha[3] = (0.8164965809277261*BdriftX[0]*m_*wv)/(dfac_v*q_); 
  double f0Quad[16]; 
  double f1Quad[16]; 
  double limQuad[16]; 
  // determine upwinding at each surface quadrature node 
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[0] = 0.25*(fl[30]-1.0*(fl[25]+fl[24]+fl[22]+fl[19])+fl[15]+fl[14]+fl[13]+fl[11]+fl[10]+fl[8]-1.0*(fl[5]+fl[4]+fl[3]+fl[2])+fl[0]); 
    f1Quad[0] = -0.25*(fl[31]-1.0*(fl[29]+fl[28]+fl[27]+fl[26])+fl[23]+fl[21]+fl[20]+fl[18]+fl[17]+fl[16]-1.0*(fl[12]+fl[9]+fl[7]+fl[6])+fl[1]); 
    limQuad[0] = fl[0]/cflL; 
  } else {
    f0Quad[0] = 0.25*(fr[30]-1.0*(fr[25]+fr[24]+fr[22]+fr[19])+fr[15]+fr[14]+fr[13]+fr[11]+fr[10]+fr[8]-1.0*(fr[5]+fr[4]+fr[3]+fr[2])+fr[0]); 
    f1Quad[0] = 0.25*(fr[31]-1.0*(fr[29]+fr[28]+fr[27]+fr[26])+fr[23]+fr[21]+fr[20]+fr[18]+fr[17]+fr[16]-1.0*(fr[12]+fr[9]+fr[7]+fr[6])+fr[1]); 
    limQuad[0] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[1] = -0.25*(fl[30]+fl[25]-1.0*(fl[24]+fl[22]+fl[19]+fl[15]+fl[14])+fl[13]-1.0*fl[11]+fl[10]+fl[8]+fl[5]+fl[4]+fl[3]-1.0*(fl[2]+fl[0])); 
    f1Quad[1] = 0.25*(fl[31]+fl[29]-1.0*(fl[28]+fl[27]+fl[26]+fl[23]+fl[21])+fl[20]-1.0*fl[18]+fl[17]+fl[16]+fl[12]+fl[9]+fl[7]-1.0*(fl[6]+fl[1])); 
    limQuad[1] = fl[0]/cflL; 
  } else {
    f0Quad[1] = -0.25*(fr[30]+fr[25]-1.0*(fr[24]+fr[22]+fr[19]+fr[15]+fr[14])+fr[13]-1.0*fr[11]+fr[10]+fr[8]+fr[5]+fr[4]+fr[3]-1.0*(fr[2]+fr[0])); 
    f1Quad[1] = -0.25*(fr[31]+fr[29]-1.0*(fr[28]+fr[27]+fr[26]+fr[23]+fr[21])+fr[20]-1.0*fr[18]+fr[17]+fr[16]+fr[12]+fr[9]+fr[7]-1.0*(fr[6]+fr[1])); 
    limQuad[1] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[2] = -0.25*(fl[30]-1.0*fl[25]+fl[24]-1.0*(fl[22]+fl[19]+fl[15])+fl[14]-1.0*fl[13]+fl[11]-1.0*fl[10]+fl[8]+fl[5]+fl[4]-1.0*fl[3]+fl[2]-1.0*fl[0]); 
    f1Quad[2] = 0.25*(fl[31]-1.0*fl[29]+fl[28]-1.0*(fl[27]+fl[26]+fl[23])+fl[21]-1.0*fl[20]+fl[18]-1.0*fl[17]+fl[16]+fl[12]+fl[9]-1.0*fl[7]+fl[6]-1.0*fl[1]); 
    limQuad[2] = fl[0]/cflL; 
  } else {
    f0Quad[2] = -0.25*(fr[30]-1.0*fr[25]+fr[24]-1.0*(fr[22]+fr[19]+fr[15])+fr[14]-1.0*fr[13]+fr[11]-1.0*fr[10]+fr[8]+fr[5]+fr[4]-1.0*fr[3]+fr[2]-1.0*fr[0]); 
    f1Quad[2] = -0.25*(fr[31]-1.0*fr[29]+fr[28]-1.0*(fr[27]+fr[26]+fr[23])+fr[21]-1.0*fr[20]+fr[18]-1.0*fr[17]+fr[16]+fr[12]+fr[9]-1.0*fr[7]+fr[6]-1.0*fr[1]); 
    limQuad[2] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[3] = 0.25*(fl[30]+fl[25]+fl[24]-1.0*(fl[22]+fl[19])+fl[15]-1.0*(fl[14]+fl[13]+fl[11]+fl[10])+fl[8]-1.0*(fl[5]+fl[4])+fl[3]+fl[2]+fl[0]); 
    f1Quad[3] = -0.25*(fl[31]+fl[29]+fl[28]-1.0*(fl[27]+fl[26])+fl[23]-1.0*(fl[21]+fl[20]+fl[18]+fl[17])+fl[16]-1.0*(fl[12]+fl[9])+fl[7]+fl[6]+fl[1]); 
    limQuad[3] = fl[0]/cflL; 
  } else {
    f0Quad[3] = 0.25*(fr[30]+fr[25]+fr[24]-1.0*(fr[22]+fr[19])+fr[15]-1.0*(fr[14]+fr[13]+fr[11]+fr[10])+fr[8]-1.0*(fr[5]+fr[4])+fr[3]+fr[2]+fr[0]); 
    f1Quad[3] = 0.25*(fr[31]+fr[29]+fr[28]-1.0*(fr[27]+fr[26])+fr[23]-1.0*(fr[21]+fr[20]+fr[18]+fr[17])+fr[16]-1.0*(fr[12]+fr[9])+fr[7]+fr[6]+fr[1]); 
    limQuad[3] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[4] = -0.25*(fl[30]-1.0*(fl[25]+fl[24])+fl[22]-1.0*fl[19]+fl[15]-1.0*(fl[14]+fl[13])+fl[11]+fl[10]-1.0*fl[8]+fl[5]-1.0*fl[4]+fl[3]+fl[2]-1.0*fl[0]); 
    f1Quad[4] = 0.25*(fl[31]-1.0*(fl[29]+fl[28])+fl[27]-1.0*fl[26]+fl[23]-1.0*(fl[21]+fl[20])+fl[18]+fl[17]-1.0*fl[16]+fl[12]-1.0*fl[9]+fl[7]+fl[6]-1.0*fl[1]); 
    limQuad[4] = fl[0]/cflL; 
  } else {
    f0Quad[4] = -0.25*(fr[30]-1.0*(fr[25]+fr[24])+fr[22]-1.0*fr[19]+fr[15]-1.0*(fr[14]+fr[13])+fr[11]+fr[10]-1.0*fr[8]+fr[5]-1.0*fr[4]+fr[3]+fr[2]-1.0*fr[0]); 
    f1Quad[4] = -0.25*(fr[31]-1.0*(fr[29]+fr[28])+fr[27]-1.0*fr[26]+fr[23]-1.0*(fr[21]+fr[20])+fr[18]+fr[17]-1.0*fr[16]+fr[12]-1.0*fr[9]+fr[7]+fr[6]-1.0*fr[1]); 
    limQuad[4] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[5] = 0.25*(fl[30]+fl[25]-1.0*fl[24]+fl[22]-1.0*(fl[19]+fl[15])+fl[14]-1.0*(fl[13]+fl[11])+fl[10]-1.0*(fl[8]+fl[5])+fl[4]-1.0*fl[3]+fl[2]+fl[0]); 
    f1Quad[5] = -0.25*(fl[31]+fl[29]-1.0*fl[28]+fl[27]-1.0*(fl[26]+fl[23])+fl[21]-1.0*(fl[20]+fl[18])+fl[17]-1.0*(fl[16]+fl[12])+fl[9]-1.0*fl[7]+fl[6]+fl[1]); 
    limQuad[5] = fl[0]/cflL; 
  } else {
    f0Quad[5] = 0.25*(fr[30]+fr[25]-1.0*fr[24]+fr[22]-1.0*(fr[19]+fr[15])+fr[14]-1.0*(fr[13]+fr[11])+fr[10]-1.0*(fr[8]+fr[5])+fr[4]-1.0*fr[3]+fr[2]+fr[0]); 
    f1Quad[5] = 0.25*(fr[31]+fr[29]-1.0*fr[28]+fr[27]-1.0*(fr[26]+fr[23])+fr[21]-1.0*(fr[20]+fr[18])+fr[17]-1.0*(fr[16]+fr[12])+fr[9]-1.0*fr[7]+fr[6]+fr[1]); 
    limQuad[5] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[6] = 0.25*(fl[30]-1.0*fl[25]+fl[24]+fl[22]-1.0*(fl[19]+fl[15]+fl[14])+fl[13]+fl[11]-1.0*(fl[10]+fl[8]+fl[5])+fl[4]+fl[3]-1.0*fl[2]+fl[0]); 
    f1Quad[6] = -0.25*(fl[31]-1.0*fl[29]+fl[28]+fl[27]-1.0*(fl[26]+fl[23]+fl[21])+fl[20]+fl[18]-1.0*(fl[17]+fl[16]+fl[12])+fl[9]+fl[7]-1.0*fl[6]+fl[1]); 
    limQuad[6] = fl[0]/cflL; 
  } else {
    f0Quad[6] = 0.25*(fr[30]-1.0*fr[25]+fr[24]+fr[22]-1.0*(fr[19]+fr[15]+fr[14])+fr[13]+fr[11]-1.0*(fr[10]+fr[8]+fr[5])+fr[4]+fr[3]-1.0*fr[2]+fr[0]); 
    f1Quad[6] = 0.25*(fr[31]-1.0*fr[29]+fr[28]+fr[27]-1.0*(fr[26]+fr[23]+fr[21])+fr[20]+fr[18]-1.0*(fr[17]+fr[16]+fr[12])+fr[9]+fr[7]-1.0*fr[6]+fr[1]); 
    limQuad[6] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[7] = -0.25*(fl[30]+fl[25]+fl[24]+fl[22]-1.0*fl[19]+fl[15]+fl[14]+fl[13]-1.0*(fl[11]+fl[10]+fl[8])+fl[5]-1.0*(fl[4]+fl[3]+fl[2]+fl[0])); 
    f1Quad[7] = 0.25*(fl[31]+fl[29]+fl[28]+fl[27]-1.0*fl[26]+fl[23]+fl[21]+fl[20]-1.0*(fl[18]+fl[17]+fl[16])+fl[12]-1.0*(fl[9]+fl[7]+fl[6]+fl[1])); 
    limQuad[7] = fl[0]/cflL; 
  } else {
    f0Quad[7] = -0.25*(fr[30]+fr[25]+fr[24]+fr[22]-1.0*fr[19]+fr[15]+fr[14]+fr[13]-1.0*(fr[11]+fr[10]+fr[8])+fr[5]-1.0*(fr[4]+fr[3]+fr[2]+fr[0])); 
    f1Quad[7] = -0.25*(fr[31]+fr[29]+fr[28]+fr[27]-1.0*fr[26]+fr[23]+fr[21]+fr[20]-1.0*(fr[18]+fr[17]+fr[16])+fr[12]-1.0*(fr[9]+fr[7]+fr[6]+fr[1])); 
    limQuad[7] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[8] = -0.25*(fl[30]-1.0*(fl[25]+fl[24]+fl[22])+fl[19]+fl[15]+fl[14]+fl[13]-1.0*(fl[11]+fl[10]+fl[8]+fl[5])+fl[4]+fl[3]+fl[2]-1.0*fl[0]); 
    f1Quad[8] = 0.25*(fl[31]-1.0*(fl[29]+fl[28]+fl[27])+fl[26]+fl[23]+fl[21]+fl[20]-1.0*(fl[18]+fl[17]+fl[16]+fl[12])+fl[9]+fl[7]+fl[6]-1.0*fl[1]); 
    limQuad[8] = fl[0]/cflL; 
  } else {
    f0Quad[8] = -0.25*(fr[30]-1.0*(fr[25]+fr[24]+fr[22])+fr[19]+fr[15]+fr[14]+fr[13]-1.0*(fr[11]+fr[10]+fr[8]+fr[5])+fr[4]+fr[3]+fr[2]-1.0*fr[0]); 
    f1Quad[8] = -0.25*(fr[31]-1.0*(fr[29]+fr[28]+fr[27])+fr[26]+fr[23]+fr[21]+fr[20]-1.0*(fr[18]+fr[17]+fr[16]+fr[12])+fr[9]+fr[7]+fr[6]-1.0*fr[1]); 
    limQuad[8] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[9] = 0.25*(fl[30]+fl[25]-1.0*(fl[24]+fl[22])+fl[19]-1.0*(fl[15]+fl[14])+fl[13]+fl[11]-1.0*(fl[10]+fl[8])+fl[5]-1.0*(fl[4]+fl[3])+fl[2]+fl[0]); 
    f1Quad[9] = -0.25*(fl[31]+fl[29]-1.0*(fl[28]+fl[27])+fl[26]-1.0*(fl[23]+fl[21])+fl[20]+fl[18]-1.0*(fl[17]+fl[16])+fl[12]-1.0*(fl[9]+fl[7])+fl[6]+fl[1]); 
    limQuad[9] = fl[0]/cflL; 
  } else {
    f0Quad[9] = 0.25*(fr[30]+fr[25]-1.0*(fr[24]+fr[22])+fr[19]-1.0*(fr[15]+fr[14])+fr[13]+fr[11]-1.0*(fr[10]+fr[8])+fr[5]-1.0*(fr[4]+fr[3])+fr[2]+fr[0]); 
    f1Quad[9] = 0.25*(fr[31]+fr[29]-1.0*(fr[28]+fr[27])+fr[26]-1.0*(fr[23]+fr[21])+fr[20]+fr[18]-1.0*(fr[17]+fr[16])+fr[12]-1.0*(fr[9]+fr[7])+fr[6]+fr[1]); 
    limQuad[9] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[10] = 0.25*(fl[30]-1.0*fl[25]+fl[24]-1.0*fl[22]+fl[19]-1.0*fl[15]+fl[14]-1.0*(fl[13]+fl[11])+fl[10]-1.0*fl[8]+fl[5]-1.0*fl[4]+fl[3]-1.0*fl[2]+fl[0]); 
    f1Quad[10] = -0.25*(fl[31]-1.0*fl[29]+fl[28]-1.0*fl[27]+fl[26]-1.0*fl[23]+fl[21]-1.0*(fl[20]+fl[18])+fl[17]-1.0*fl[16]+fl[12]-1.0*fl[9]+fl[7]-1.0*fl[6]+fl[1]); 
    limQuad[10] = fl[0]/cflL; 
  } else {
    f0Quad[10] = 0.25*(fr[30]-1.0*fr[25]+fr[24]-1.0*fr[22]+fr[19]-1.0*fr[15]+fr[14]-1.0*(fr[13]+fr[11])+fr[10]-1.0*fr[8]+fr[5]-1.0*fr[4]+fr[3]-1.0*fr[2]+fr[0]); 
    f1Quad[10] = 0.25*(fr[31]-1.0*fr[29]+fr[28]-1.0*fr[27]+fr[26]-1.0*fr[23]+fr[21]-1.0*(fr[20]+fr[18])+fr[17]-1.0*fr[16]+fr[12]-1.0*fr[9]+fr[7]-1.0*fr[6]+fr[1]); 
    limQuad[10] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[11] = -0.25*(fl[30]+fl[25]+fl[24]-1.0*fl[22]+fl[19]+fl[15]-1.0*(fl[14]+fl[13])+fl[11]+fl[10]-1.0*(fl[8]+fl[5])+fl[4]-1.0*(fl[3]+fl[2]+fl[0])); 
    f1Quad[11] = 0.25*(fl[31]+fl[29]+fl[28]-1.0*fl[27]+fl[26]+fl[23]-1.0*(fl[21]+fl[20])+fl[18]+fl[17]-1.0*(fl[16]+fl[12])+fl[9]-1.0*(fl[7]+fl[6]+fl[1])); 
    limQuad[11] = fl[0]/cflL; 
  } else {
    f0Quad[11] = -0.25*(fr[30]+fr[25]+fr[24]-1.0*fr[22]+fr[19]+fr[15]-1.0*(fr[14]+fr[13])+fr[11]+fr[10]-1.0*(fr[8]+fr[5])+fr[4]-1.0*(fr[3]+fr[2]+fr[0])); 
    f1Quad[11] = -0.25*(fr[31]+fr[29]+fr[28]-1.0*fr[27]+fr[26]+fr[23]-1.0*(fr[21]+fr[20])+fr[18]+fr[17]-1.0*(fr[16]+fr[12])+fr[9]-1.0*(fr[7]+fr[6]+fr[1])); 
    limQuad[11] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[12] = 0.25*(fl[30]-1.0*(fl[25]+fl[24])+fl[22]+fl[19]+fl[15]-1.0*(fl[14]+fl[13]+fl[11]+fl[10])+fl[8]+fl[5]+fl[4]-1.0*(fl[3]+fl[2])+fl[0]); 
    f1Quad[12] = -0.25*(fl[31]-1.0*(fl[29]+fl[28])+fl[27]+fl[26]+fl[23]-1.0*(fl[21]+fl[20]+fl[18]+fl[17])+fl[16]+fl[12]+fl[9]-1.0*(fl[7]+fl[6])+fl[1]); 
    limQuad[12] = fl[0]/cflL; 
  } else {
    f0Quad[12] = 0.25*(fr[30]-1.0*(fr[25]+fr[24])+fr[22]+fr[19]+fr[15]-1.0*(fr[14]+fr[13]+fr[11]+fr[10])+fr[8]+fr[5]+fr[4]-1.0*(fr[3]+fr[2])+fr[0]); 
    f1Quad[12] = 0.25*(fr[31]-1.0*(fr[29]+fr[28])+fr[27]+fr[26]+fr[23]-1.0*(fr[21]+fr[20]+fr[18]+fr[17])+fr[16]+fr[12]+fr[9]-1.0*(fr[7]+fr[6])+fr[1]); 
    limQuad[12] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[13] = -0.25*(fl[30]+fl[25]-1.0*fl[24]+fl[22]+fl[19]-1.0*fl[15]+fl[14]-1.0*fl[13]+fl[11]-1.0*fl[10]+fl[8]-1.0*(fl[5]+fl[4])+fl[3]-1.0*(fl[2]+fl[0])); 
    f1Quad[13] = 0.25*(fl[31]+fl[29]-1.0*fl[28]+fl[27]+fl[26]-1.0*fl[23]+fl[21]-1.0*fl[20]+fl[18]-1.0*fl[17]+fl[16]-1.0*(fl[12]+fl[9])+fl[7]-1.0*(fl[6]+fl[1])); 
    limQuad[13] = fl[0]/cflL; 
  } else {
    f0Quad[13] = -0.25*(fr[30]+fr[25]-1.0*fr[24]+fr[22]+fr[19]-1.0*fr[15]+fr[14]-1.0*fr[13]+fr[11]-1.0*fr[10]+fr[8]-1.0*(fr[5]+fr[4])+fr[3]-1.0*(fr[2]+fr[0])); 
    f1Quad[13] = -0.25*(fr[31]+fr[29]-1.0*fr[28]+fr[27]+fr[26]-1.0*fr[23]+fr[21]-1.0*fr[20]+fr[18]-1.0*fr[17]+fr[16]-1.0*(fr[12]+fr[9])+fr[7]-1.0*(fr[6]+fr[1])); 
    limQuad[13] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[14] = -0.25*(fl[30]-1.0*fl[25]+fl[24]+fl[22]+fl[19]-1.0*(fl[15]+fl[14])+fl[13]-1.0*fl[11]+fl[10]+fl[8]-1.0*(fl[5]+fl[4]+fl[3])+fl[2]-1.0*fl[0]); 
    f1Quad[14] = 0.25*(fl[31]-1.0*fl[29]+fl[28]+fl[27]+fl[26]-1.0*(fl[23]+fl[21])+fl[20]-1.0*fl[18]+fl[17]+fl[16]-1.0*(fl[12]+fl[9]+fl[7])+fl[6]-1.0*fl[1]); 
    limQuad[14] = fl[0]/cflL; 
  } else {
    f0Quad[14] = -0.25*(fr[30]-1.0*fr[25]+fr[24]+fr[22]+fr[19]-1.0*(fr[15]+fr[14])+fr[13]-1.0*fr[11]+fr[10]+fr[8]-1.0*(fr[5]+fr[4]+fr[3])+fr[2]-1.0*fr[0]); 
    f1Quad[14] = -0.25*(fr[31]-1.0*fr[29]+fr[28]+fr[27]+fr[26]-1.0*(fr[23]+fr[21])+fr[20]-1.0*fr[18]+fr[17]+fr[16]-1.0*(fr[12]+fr[9]+fr[7])+fr[6]-1.0*fr[1]); 
    limQuad[14] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[15] = 0.25*(fl[30]+fl[25]+fl[24]+fl[22]+fl[19]+fl[15]+fl[14]+fl[13]+fl[11]+fl[10]+fl[8]+fl[5]+fl[4]+fl[3]+fl[2]+fl[0]); 
    f1Quad[15] = -0.25*(fl[31]+fl[29]+fl[28]+fl[27]+fl[26]+fl[23]+fl[21]+fl[20]+fl[18]+fl[17]+fl[16]+fl[12]+fl[9]+fl[7]+fl[6]+fl[1]); 
    limQuad[15] = fl[0]/cflL; 
  } else {
    f0Quad[15] = 0.25*(fr[30]+fr[25]+fr[24]+fr[22]+fr[19]+fr[15]+fr[14]+fr[13]+fr[11]+fr[10]+fr[8]+fr[5]+fr[4]+fr[3]+fr[2]+fr[0]); 
    f1Quad[15] = 0.25*(fr[31]+fr[29]+fr[28]+fr[27]+fr[26]+fr[23]+fr[21]+fr[20]+fr[18]+fr[17]+fr[16]+fr[12]+fr[9]+fr[7]+fr[6]+fr[1]); 
    limQuad[15] = fr[0]/cflR; 
  }
  double fhat[32]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[2] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[3] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[6] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[7] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[8] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[9] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[10] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[11] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[12] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[13] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[14] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[15] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[16] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[17] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[18] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[19] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[20] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[21] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[22] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[23] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[24] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[25] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[26] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[27] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[28] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[29] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[30] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[31] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[16];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrl[0] = (1.414213562373095*(3.0*fhat[31]-5.196152422706631*(fhat[29]+fhat[28]+fhat[27]+fhat[26])+9.0*(fhat[23]+fhat[21]+fhat[20]+fhat[18]+fhat[17]+fhat[16])-15.58845726811989*(fhat[12]+fhat[9]+fhat[7]+fhat[6])+27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*(1.732050807568877*fhat[30]-3.0*(fhat[25]+fhat[24]+fhat[22]+fhat[19])+5.196152422706631*(fhat[15]+fhat[14]+fhat[13]+fhat[11]+fhat[10]+fhat[8])-9.0*(fhat[5]+fhat[4]+fhat[3]+fhat[2])+15.58845726811989*fhat[0])); 
  rCtrl[1] = -(1.414213562373095*(3.0*fhat[31]+5.196152422706631*(fhat[29]-1.0*(fhat[28]+fhat[27]+fhat[26]))+9.0*((-1.0*(fhat[23]+fhat[21]))+fhat[20]-1.0*fhat[18]+fhat[17]+fhat[16])+15.58845726811989*(fhat[12]+fhat[9]+fhat[7])-1.0*(15.58845726811989*fhat[6]+27.0*fhat[1])))/(124.7076581449591*EPSILON+1.414213562373095*((-1.732050807568877*fhat[30])+3.0*((-1.0*fhat[25])+fhat[24]+fhat[22]+fhat[19])+5.196152422706631*(fhat[15]+fhat[14]-1.0*fhat[13]+fhat[11])-1.0*(5.196152422706631*(fhat[10]+fhat[8])+9.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[2])-15.58845726811989*fhat[0]))); 
  rCtrl[2] = -(1.414213562373095*(3.0*fhat[31]+5.196152422706631*((-1.0*fhat[29])+fhat[28]-1.0*(fhat[27]+fhat[26]))+9.0*((-1.0*fhat[23])+fhat[21]-1.0*fhat[20]+fhat[18]-1.0*fhat[17]+fhat[16])+15.58845726811989*(fhat[12]+fhat[9]-1.0*fhat[7]+fhat[6])-27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*((-1.732050807568877*fhat[30])+3.0*(fhat[25]-1.0*fhat[24]+fhat[22]+fhat[19])+5.196152422706631*(fhat[15]-1.0*fhat[14]+fhat[13]-1.0*fhat[11]+fhat[10]-1.0*fhat[8])+9.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[2])+15.58845726811989*fhat[0])); 
  rCtrl[3] = (1.414213562373095*(3.0*fhat[31]+5.196152422706631*(fhat[29]+fhat[28]-1.0*(fhat[27]+fhat[26]))+9.0*(fhat[23]-1.0*(fhat[21]+fhat[20]+fhat[18]+fhat[17]-1.0*fhat[16]))+15.58845726811989*((-1.0*(fhat[12]+fhat[9]))+fhat[7]+fhat[6])+27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*(1.732050807568877*fhat[30]+3.0*(fhat[25]+fhat[24]-1.0*(fhat[22]+fhat[19]))+5.196152422706631*(fhat[15]-1.0*(fhat[14]+fhat[13]+fhat[11]+fhat[10]-1.0*fhat[8]))+9.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[2])+15.58845726811989*fhat[0])); 
  rCtrl[4] = -(1.414213562373095*(3.0*fhat[31]+5.196152422706631*((-1.0*(fhat[29]+fhat[28]))+fhat[27]-1.0*fhat[26])+9.0*(fhat[23]-1.0*(fhat[21]+fhat[20]-1.0*fhat[18])+fhat[17]-1.0*fhat[16])+15.58845726811989*(fhat[12]-1.0*fhat[9]+fhat[7]+fhat[6])-27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*((-1.732050807568877*fhat[30])+3.0*(fhat[25]+fhat[24]-1.0*fhat[22]+fhat[19])+5.196152422706631*((-1.0*fhat[15])+fhat[14]+fhat[13]-1.0*(fhat[11]+fhat[10]-1.0*fhat[8]))+9.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[2]))+15.58845726811989*fhat[0])); 
  rCtrl[5] = (1.414213562373095*(3.0*fhat[31]+5.196152422706631*(fhat[29]-1.0*fhat[28]+fhat[27]-1.0*fhat[26])+9.0*((-1.0*fhat[23])+fhat[21]-1.0*(fhat[20]+fhat[18]-1.0*fhat[17]+fhat[16]))+15.58845726811989*((-1.0*fhat[12])+fhat[9]-1.0*fhat[7]+fhat[6])+27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*(1.732050807568877*fhat[30]+3.0*(fhat[25]-1.0*fhat[24]+fhat[22]-1.0*fhat[19])+5.196152422706631*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8]))+9.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[2])+15.58845726811989*fhat[0])); 
  rCtrl[6] = (1.414213562373095*(3.0*fhat[31]+5.196152422706631*((-1.0*fhat[29])+fhat[28]+fhat[27]-1.0*fhat[26])+9.0*((-1.0*(fhat[23]+fhat[21]))+fhat[20]+fhat[18]-1.0*(fhat[17]+fhat[16]))+15.58845726811989*((-1.0*fhat[12])+fhat[9]+fhat[7]-1.0*fhat[6])+27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*(1.732050807568877*fhat[30]+3.0*((-1.0*fhat[25])+fhat[24]+fhat[22]-1.0*fhat[19])+5.196152422706631*((-1.0*(fhat[15]+fhat[14]))+fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]))+9.0*((-1.0*fhat[5])+fhat[4]+fhat[3]-1.0*fhat[2])+15.58845726811989*fhat[0])); 
  rCtrl[7] = -(1.414213562373095*(3.0*fhat[31]+5.196152422706631*(fhat[29]+fhat[28]+fhat[27]-1.0*fhat[26])+9.0*(fhat[23]+fhat[21]+fhat[20]-1.0*(fhat[18]+fhat[17]+fhat[16]))+15.58845726811989*fhat[12]-1.0*(15.58845726811989*(fhat[9]+fhat[7]+fhat[6])+27.0*fhat[1])))/(124.7076581449591*EPSILON+1.414213562373095*((-1.0*(1.732050807568877*fhat[30]+3.0*(fhat[25]+fhat[24]+fhat[22]-1.0*fhat[19])))+5.196152422706631*((-1.0*(fhat[15]+fhat[14]+fhat[13]-1.0*fhat[11]))+fhat[10]+fhat[8])+9.0*((-1.0*fhat[5])+fhat[4]+fhat[3]+fhat[2])+15.58845726811989*fhat[0])); 
  rCtrl[8] = -(1.414213562373095*(3.0*fhat[31]-5.196152422706631*(fhat[29]+fhat[28]+fhat[27]-1.0*fhat[26])+9.0*(fhat[23]+fhat[21]+fhat[20]-1.0*(fhat[18]+fhat[17]+fhat[16]))+15.58845726811989*((-1.0*fhat[12])+fhat[9]+fhat[7]+fhat[6])-27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*((-1.732050807568877*fhat[30])+3.0*(fhat[25]+fhat[24]+fhat[22]-1.0*fhat[19])+5.196152422706631*((-1.0*(fhat[15]+fhat[14]+fhat[13]-1.0*fhat[11]))+fhat[10]+fhat[8])+9.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[2]))+15.58845726811989*fhat[0])); 
  rCtrl[9] = (1.414213562373095*(3.0*fhat[31]+5.196152422706631*(fhat[29]-1.0*(fhat[28]+fhat[27]-1.0*fhat[26]))+9.0*((-1.0*(fhat[23]+fhat[21]))+fhat[20]+fhat[18]-1.0*(fhat[17]+fhat[16]))+15.58845726811989*(fhat[12]-1.0*(fhat[9]+fhat[7]-1.0*fhat[6]))+27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*(1.732050807568877*fhat[30]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[22]-1.0*fhat[19]))+5.196152422706631*((-1.0*(fhat[15]+fhat[14]))+fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]))+9.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[2]))+15.58845726811989*fhat[0])); 
  rCtrl[10] = (1.414213562373095*(3.0*fhat[31]+5.196152422706631*((-1.0*fhat[29])+fhat[28]-1.0*fhat[27]+fhat[26])+9.0*((-1.0*fhat[23])+fhat[21]-1.0*(fhat[20]+fhat[18]-1.0*fhat[17]+fhat[16]))+15.58845726811989*(fhat[12]-1.0*fhat[9]+fhat[7]-1.0*fhat[6])+27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*(1.732050807568877*fhat[30]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*fhat[22]+fhat[19])+5.196152422706631*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8]))+9.0*(fhat[5]-1.0*fhat[4]+fhat[3]-1.0*fhat[2])+15.58845726811989*fhat[0])); 
  rCtrl[11] = -(1.414213562373095*(3.0*fhat[31]+5.196152422706631*(fhat[29]+fhat[28]-1.0*fhat[27]+fhat[26])+9.0*(fhat[23]-1.0*(fhat[21]+fhat[20]-1.0*fhat[18])+fhat[17]-1.0*fhat[16])+15.58845726811989*(fhat[9]-1.0*fhat[12])-1.0*(15.58845726811989*(fhat[7]+fhat[6])+27.0*fhat[1])))/(124.7076581449591*EPSILON+1.414213562373095*((-1.732050807568877*fhat[30])+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[22]-1.0*fhat[19])+5.196152422706631*((-1.0*fhat[15])+fhat[14]+fhat[13]-1.0*(fhat[11]+fhat[10]-1.0*fhat[8]))+9.0*(fhat[5]-1.0*fhat[4]+fhat[3]+fhat[2])+15.58845726811989*fhat[0])); 
  rCtrl[12] = (1.414213562373095*(3.0*fhat[31]+5.196152422706631*((-1.0*(fhat[29]+fhat[28]))+fhat[27]+fhat[26])+9.0*(fhat[23]-1.0*(fhat[21]+fhat[20]+fhat[18]+fhat[17]-1.0*fhat[16]))+15.58845726811989*(fhat[12]+fhat[9]-1.0*(fhat[7]+fhat[6]))+27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*(1.732050807568877*fhat[30]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[22]+fhat[19])+5.196152422706631*(fhat[15]-1.0*(fhat[14]+fhat[13]+fhat[11]+fhat[10]-1.0*fhat[8]))+9.0*(fhat[5]+fhat[4]-1.0*(fhat[3]+fhat[2]))+15.58845726811989*fhat[0])); 
  rCtrl[13] = -(1.414213562373095*(3.0*fhat[31]+5.196152422706631*(fhat[29]-1.0*fhat[28]+fhat[27]+fhat[26])+9.0*((-1.0*fhat[23])+fhat[21]-1.0*fhat[20]+fhat[18]-1.0*fhat[17]+fhat[16])+15.58845726811989*(fhat[7]-1.0*(fhat[12]+fhat[9]))-1.0*(15.58845726811989*fhat[6]+27.0*fhat[1])))/(124.7076581449591*EPSILON+1.414213562373095*((-1.732050807568877*fhat[30])+3.0*((-1.0*fhat[25])+fhat[24]-1.0*(fhat[22]+fhat[19]))+5.196152422706631*(fhat[15]-1.0*fhat[14]+fhat[13]-1.0*fhat[11]+fhat[10]-1.0*fhat[8])+9.0*(fhat[5]+fhat[4]-1.0*fhat[3]+fhat[2])+15.58845726811989*fhat[0])); 
  rCtrl[14] = -(1.414213562373095*(3.0*fhat[31]+5.196152422706631*((-1.0*fhat[29])+fhat[28]+fhat[27]+fhat[26])+9.0*((-1.0*(fhat[23]+fhat[21]))+fhat[20]-1.0*fhat[18]+fhat[17]+fhat[16])-15.58845726811989*(fhat[12]+fhat[9]+fhat[7]-1.0*fhat[6])-27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*((-1.732050807568877*fhat[30])+3.0*(fhat[25]-1.0*(fhat[24]+fhat[22]+fhat[19]))+5.196152422706631*(fhat[15]+fhat[14]-1.0*fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]))+9.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[2])+15.58845726811989*fhat[0])); 
  rCtrl[15] = (1.414213562373095*(3.0*fhat[31]+5.196152422706631*(fhat[29]+fhat[28]+fhat[27]+fhat[26])+9.0*(fhat[23]+fhat[21]+fhat[20]+fhat[18]+fhat[17]+fhat[16])+15.58845726811989*(fhat[12]+fhat[9]+fhat[7]+fhat[6])+27.0*fhat[1]))/(124.7076581449591*EPSILON+1.414213562373095*(1.732050807568877*fhat[30]+3.0*(fhat[25]+fhat[24]+fhat[22]+fhat[19])+5.196152422706631*(fhat[15]+fhat[14]+fhat[13]+fhat[11]+fhat[10]+fhat[8])+9.0*(fhat[5]+fhat[4]+fhat[3]+fhat[2])+15.58845726811989*fhat[0])); 
  double fhatCtrl[16];  // fhatCtrl = anti-limited fhat evaluated at each control node on x surface 
  fhatCtrl[0] = 0.01134023029066286*(1.732050807568877*fhat[30]-3.0*(fhat[25]+fhat[24]+fhat[22]+fhat[19])+5.196152422706631*(fhat[15]+fhat[14]+fhat[13]+fhat[11]+fhat[10]+fhat[8])-9.0*(fhat[5]+fhat[4]+fhat[3]+fhat[2])+15.58845726811989*fhat[0])*limTheta(rCtrl[0],-1.0,EPSILON); 
  fhatCtrl[1] = -0.01134023029066286*(1.732050807568877*fhat[30]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[22]+fhat[19]))+5.196152422706631*((-1.0*(fhat[15]+fhat[14]))+fhat[13]-1.0*fhat[11]+fhat[10]+fhat[8])+9.0*(fhat[5]+fhat[4]+fhat[3])-1.0*(9.0*fhat[2]+15.58845726811989*fhat[0]))*limTheta(rCtrl[1],-1.0,EPSILON); 
  fhatCtrl[2] = -0.01134023029066286*(1.732050807568877*fhat[30]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*(fhat[22]+fhat[19]))+5.196152422706631*((-1.0*fhat[15])+fhat[14]-1.0*fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8])+9.0*(fhat[5]+fhat[4]-1.0*fhat[3]+fhat[2])-15.58845726811989*fhat[0])*limTheta(rCtrl[2],-1.0,EPSILON); 
  fhatCtrl[3] = 0.01134023029066286*(1.732050807568877*fhat[30]+3.0*(fhat[25]+fhat[24]-1.0*(fhat[22]+fhat[19]))+5.196152422706631*(fhat[15]-1.0*(fhat[14]+fhat[13]+fhat[11]+fhat[10]-1.0*fhat[8]))+9.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[2])+15.58845726811989*fhat[0])*limTheta(rCtrl[3],-1.0,EPSILON); 
  fhatCtrl[4] = -0.01134023029066286*(1.732050807568877*fhat[30]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[22]-1.0*fhat[19])+5.196152422706631*(fhat[15]-1.0*(fhat[14]+fhat[13]-1.0*fhat[11])+fhat[10]-1.0*fhat[8])+9.0*(fhat[5]-1.0*fhat[4]+fhat[3]+fhat[2])-15.58845726811989*fhat[0])*limTheta(rCtrl[4],-1.0,EPSILON); 
  fhatCtrl[5] = 0.01134023029066286*(1.732050807568877*fhat[30]+3.0*(fhat[25]-1.0*fhat[24]+fhat[22]-1.0*fhat[19])+5.196152422706631*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8]))+9.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[2])+15.58845726811989*fhat[0])*limTheta(rCtrl[5],-1.0,EPSILON); 
  fhatCtrl[6] = 0.01134023029066286*(1.732050807568877*fhat[30]+3.0*((-1.0*fhat[25])+fhat[24]+fhat[22]-1.0*fhat[19])+5.196152422706631*((-1.0*(fhat[15]+fhat[14]))+fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]))+9.0*((-1.0*fhat[5])+fhat[4]+fhat[3]-1.0*fhat[2])+15.58845726811989*fhat[0])*limTheta(rCtrl[6],-1.0,EPSILON); 
  fhatCtrl[7] = -0.01134023029066286*(1.732050807568877*fhat[30]+3.0*(fhat[25]+fhat[24]+fhat[22]-1.0*fhat[19])+5.196152422706631*(fhat[15]+fhat[14]+fhat[13]-1.0*(fhat[11]+fhat[10]+fhat[8]))+9.0*fhat[5]-1.0*(9.0*(fhat[4]+fhat[3]+fhat[2])+15.58845726811989*fhat[0]))*limTheta(rCtrl[7],-1.0,EPSILON); 
  fhatCtrl[8] = -0.01134023029066286*(1.732050807568877*fhat[30]-3.0*(fhat[25]+fhat[24]+fhat[22]-1.0*fhat[19])+5.196152422706631*(fhat[15]+fhat[14]+fhat[13]-1.0*(fhat[11]+fhat[10]+fhat[8]))+9.0*((-1.0*fhat[5])+fhat[4]+fhat[3]+fhat[2])-15.58845726811989*fhat[0])*limTheta(rCtrl[8],-1.0,EPSILON); 
  fhatCtrl[9] = 0.01134023029066286*(1.732050807568877*fhat[30]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[22]-1.0*fhat[19]))+5.196152422706631*((-1.0*(fhat[15]+fhat[14]))+fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]))+9.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[2]))+15.58845726811989*fhat[0])*limTheta(rCtrl[9],-1.0,EPSILON); 
  fhatCtrl[10] = 0.01134023029066286*(1.732050807568877*fhat[30]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*fhat[22]+fhat[19])+5.196152422706631*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8]))+9.0*(fhat[5]-1.0*fhat[4]+fhat[3]-1.0*fhat[2])+15.58845726811989*fhat[0])*limTheta(rCtrl[10],-1.0,EPSILON); 
  fhatCtrl[11] = -0.01134023029066286*(1.732050807568877*fhat[30]+3.0*(fhat[25]+fhat[24]-1.0*fhat[22]+fhat[19])+5.196152422706631*(fhat[15]-1.0*(fhat[14]+fhat[13]-1.0*fhat[11])+fhat[10]-1.0*fhat[8])+9.0*(fhat[4]-1.0*fhat[5])-1.0*(9.0*(fhat[3]+fhat[2])+15.58845726811989*fhat[0]))*limTheta(rCtrl[11],-1.0,EPSILON); 
  fhatCtrl[12] = 0.01134023029066286*(1.732050807568877*fhat[30]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[22]+fhat[19])+5.196152422706631*(fhat[15]-1.0*(fhat[14]+fhat[13]+fhat[11]+fhat[10]-1.0*fhat[8]))+9.0*(fhat[5]+fhat[4]-1.0*(fhat[3]+fhat[2]))+15.58845726811989*fhat[0])*limTheta(rCtrl[12],-1.0,EPSILON); 
  fhatCtrl[13] = -0.01134023029066286*(1.732050807568877*fhat[30]+3.0*(fhat[25]-1.0*fhat[24]+fhat[22]+fhat[19])+5.196152422706631*((-1.0*fhat[15])+fhat[14]-1.0*fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8])+9.0*(fhat[3]-1.0*(fhat[5]+fhat[4]))-1.0*(9.0*fhat[2]+15.58845726811989*fhat[0]))*limTheta(rCtrl[13],-1.0,EPSILON); 
  fhatCtrl[14] = -0.01134023029066286*(1.732050807568877*fhat[30]+3.0*((-1.0*fhat[25])+fhat[24]+fhat[22]+fhat[19])+5.196152422706631*((-1.0*(fhat[15]+fhat[14]))+fhat[13]-1.0*fhat[11]+fhat[10]+fhat[8])-9.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[2])-15.58845726811989*fhat[0])*limTheta(rCtrl[14],-1.0,EPSILON); 
  fhatCtrl[15] = 0.01134023029066286*(1.732050807568877*fhat[30]+3.0*(fhat[25]+fhat[24]+fhat[22]+fhat[19])+5.196152422706631*(fhat[15]+fhat[14]+fhat[13]+fhat[11]+fhat[10]+fhat[8])+9.0*(fhat[5]+fhat[4]+fhat[3]+fhat[2])+15.58845726811989*fhat[0])*limTheta(rCtrl[15],-1.0,EPSILON); 
  double fhatAL[16];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.25*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.4330127018922193*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[5] = 0.75*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[6] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[8] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[9] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[10] = 0.75*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[11] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[12] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[13] = 1.299038105676658*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[14] = 1.299038105676658*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[15] = 2.25*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[16]; 
  fhatALQuad[0] = fmin(0.25*((-0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7])))-1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[0]); 
  fhatALQuad[1] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*((-0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+0.5773502691896258*(fhatAL[13]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[1]); 
  fhatALQuad[2] = fmin(0.25*(0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[2]); 
  fhatALQuad[3] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[3]); 
  fhatALQuad[4] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*(0.5773502691896258*fhatAL[3]-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[4]+fhatAL[2]+fhatAL[1]))-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[4]); 
  fhatALQuad[5] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]))+0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[5]); 
  fhatALQuad[6] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2])-0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[6]); 
  fhatALQuad[7] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[7]); 
  fhatALQuad[8] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]+fhatAL[1]))+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[8]); 
  fhatALQuad[9] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]))+0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[9]); 
  fhatALQuad[10] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[10]); 
  fhatALQuad[11] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[4]+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[11]); 
  fhatALQuad[12] = fmin(0.25*((-0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5])))+1.732050807568877*(0.5773502691896258*(fhatAL[4]+fhatAL[3])-0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[12]); 
  fhatALQuad[13] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*(0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1])-0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[13]); 
  fhatALQuad[14] = fmin(0.25*(0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[1])+0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[14]); 
  fhatALQuad[15] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[15]); 
  fhatAL[0] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[5] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[6] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[8] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[9] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[10] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[11] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[12] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[13] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[14] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[15] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.1767766952966368*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_x; 
  incr[1] = -0.3061862178478971*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_x; 
  incr[2] = 0.1767766952966368*(alpha[3]*fhatAL[6]+alpha[2]*fhatAL[5]+alpha[0]*fhatAL[1])*dfac_x; 
  incr[3] = 0.1767766952966368*(alpha[3]*fhatAL[7]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_x; 
  incr[4] = 0.1767766952966368*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_x; 
  incr[5] = 0.1767766952966368*(alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[0]*fhatAL[4])*dfac_x; 
  incr[6] = -0.3061862178478971*(alpha[3]*fhatAL[6]+alpha[2]*fhatAL[5]+alpha[0]*fhatAL[1])*dfac_x; 
  incr[7] = -0.3061862178478971*(alpha[3]*fhatAL[7]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_x; 
  incr[8] = 0.1767766952966368*(alpha[3]*fhatAL[11]+alpha[0]*fhatAL[5]+fhatAL[1]*alpha[2])*dfac_x; 
  incr[9] = -0.3061862178478971*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_x; 
  incr[10] = 0.1767766952966368*(alpha[2]*fhatAL[11]+alpha[0]*fhatAL[6]+fhatAL[1]*alpha[3])*dfac_x; 
  incr[11] = 0.1767766952966368*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_x; 
  incr[12] = -0.3061862178478971*(alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[0]*fhatAL[4])*dfac_x; 
  incr[13] = 0.1767766952966368*(alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[0]*fhatAL[8])*dfac_x; 
  incr[14] = 0.1767766952966368*(alpha[3]*fhatAL[14]+alpha[0]*fhatAL[9]+alpha[2]*fhatAL[4])*dfac_x; 
  incr[15] = 0.1767766952966368*(alpha[2]*fhatAL[14]+alpha[0]*fhatAL[10]+alpha[3]*fhatAL[4])*dfac_x; 
  incr[16] = -0.3061862178478971*(alpha[3]*fhatAL[11]+alpha[0]*fhatAL[5]+fhatAL[1]*alpha[2])*dfac_x; 
  incr[17] = -0.3061862178478971*(alpha[2]*fhatAL[11]+alpha[0]*fhatAL[6]+fhatAL[1]*alpha[3])*dfac_x; 
  incr[18] = -0.3061862178478971*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_x; 
  incr[19] = 0.1767766952966368*(alpha[0]*fhatAL[11]+alpha[2]*fhatAL[6]+alpha[3]*fhatAL[5])*dfac_x; 
  incr[20] = -0.3061862178478971*(alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[0]*fhatAL[8])*dfac_x; 
  incr[21] = -0.3061862178478971*(alpha[3]*fhatAL[14]+alpha[0]*fhatAL[9]+alpha[2]*fhatAL[4])*dfac_x; 
  incr[22] = 0.1767766952966368*(alpha[3]*fhatAL[15]+alpha[0]*fhatAL[12]+alpha[2]*fhatAL[8])*dfac_x; 
  incr[23] = -0.3061862178478971*(alpha[2]*fhatAL[14]+alpha[0]*fhatAL[10]+alpha[3]*fhatAL[4])*dfac_x; 
  incr[24] = 0.1767766952966368*(alpha[2]*fhatAL[15]+alpha[0]*fhatAL[13]+alpha[3]*fhatAL[8])*dfac_x; 
  incr[25] = 0.1767766952966368*(alpha[0]*fhatAL[14]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9])*dfac_x; 
  incr[26] = -0.3061862178478971*(alpha[0]*fhatAL[11]+alpha[2]*fhatAL[6]+alpha[3]*fhatAL[5])*dfac_x; 
  incr[27] = -0.3061862178478971*(alpha[3]*fhatAL[15]+alpha[0]*fhatAL[12]+alpha[2]*fhatAL[8])*dfac_x; 
  incr[28] = -0.3061862178478971*(alpha[2]*fhatAL[15]+alpha[0]*fhatAL[13]+alpha[3]*fhatAL[8])*dfac_x; 
  incr[29] = -0.3061862178478971*(alpha[0]*fhatAL[14]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9])*dfac_x; 
  incr[30] = 0.1767766952966368*(alpha[0]*fhatAL[15]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12])*dfac_x; 
  incr[31] = -0.3061862178478971*(alpha[0]*fhatAL[15]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 
  outr[28] += incr[28]; 
  outr[29] += incr[29]; 
  outr[30] += incr[30]; 
  outr[31] += incr[31]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += -1.0*incr[19]; 
  outl[20] += incr[20]; 
  outl[21] += incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += incr[26]; 
  outl[27] += incr[27]; 
  outl[28] += incr[28]; 
  outl[29] += incr[29]; 
  outl[30] += -1.0*incr[30]; 
  outl[31] += incr[31]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity3x2vSer_Y_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_z = 2.0/dxv[2]; 
  double dfac_v = 2.0/dxv[3]; 
  double dfac_m = 2.0/dxv[4]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wz = w[2]; 
  double wv = w[3]; 
  double wm = w[4]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[32]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.0625*(2.828427124746191*BdriftY[0]*m_*wv2+BmagInv[0]*(1.732050807568877*Phi[1]-3.0*Phi[4])*dfac_x*q_))/q_; 

  double alpha[16]; 
  alpha[0] = (0.5*(2.828427124746191*BdriftY[0]*m_*wv2+BmagInv[0]*(1.732050807568877*Phi[1]-3.0*Phi[4])*dfac_x*q_))/q_; 
  alpha[2] = 0.5*BmagInv[0]*(1.732050807568877*Phi[5]-3.0*Phi[7])*dfac_x; 
  alpha[3] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  double f0Quad[16]; 
  double f1Quad[16]; 
  double limQuad[16]; 
  // determine upwinding at each surface quadrature node 
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[0] = 0.25*(fl[29]-1.0*(fl[25]+fl[23]+fl[21]+fl[18])+fl[15]+fl[14]+fl[12]+fl[11]+fl[9]+fl[7]-1.0*(fl[5]+fl[4]+fl[3]+fl[1])+fl[0]); 
    f1Quad[0] = -0.25*(fl[31]-1.0*(fl[30]+fl[28]+fl[27]+fl[26])+fl[24]+fl[22]+fl[20]+fl[19]+fl[17]+fl[16]-1.0*(fl[13]+fl[10]+fl[8]+fl[6])+fl[2]); 
    limQuad[0] = fl[0]/cflL; 
  } else {
    f0Quad[0] = 0.25*(fr[29]-1.0*(fr[25]+fr[23]+fr[21]+fr[18])+fr[15]+fr[14]+fr[12]+fr[11]+fr[9]+fr[7]-1.0*(fr[5]+fr[4]+fr[3]+fr[1])+fr[0]); 
    f1Quad[0] = 0.25*(fr[31]-1.0*(fr[30]+fr[28]+fr[27]+fr[26])+fr[24]+fr[22]+fr[20]+fr[19]+fr[17]+fr[16]-1.0*(fr[13]+fr[10]+fr[8]+fr[6])+fr[2]); 
    limQuad[0] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[1] = -0.25*(fl[29]+fl[25]-1.0*(fl[23]+fl[21]+fl[18]+fl[15]+fl[14])+fl[12]-1.0*fl[11]+fl[9]+fl[7]+fl[5]+fl[4]+fl[3]-1.0*(fl[1]+fl[0])); 
    f1Quad[1] = 0.25*(fl[31]+fl[30]-1.0*(fl[28]+fl[27]+fl[26]+fl[24]+fl[22])+fl[20]-1.0*fl[19]+fl[17]+fl[16]+fl[13]+fl[10]+fl[8]-1.0*(fl[6]+fl[2])); 
    limQuad[1] = fl[0]/cflL; 
  } else {
    f0Quad[1] = -0.25*(fr[29]+fr[25]-1.0*(fr[23]+fr[21]+fr[18]+fr[15]+fr[14])+fr[12]-1.0*fr[11]+fr[9]+fr[7]+fr[5]+fr[4]+fr[3]-1.0*(fr[1]+fr[0])); 
    f1Quad[1] = -0.25*(fr[31]+fr[30]-1.0*(fr[28]+fr[27]+fr[26]+fr[24]+fr[22])+fr[20]-1.0*fr[19]+fr[17]+fr[16]+fr[13]+fr[10]+fr[8]-1.0*(fr[6]+fr[2])); 
    limQuad[1] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[2] = -0.25*(fl[29]-1.0*fl[25]+fl[23]-1.0*(fl[21]+fl[18]+fl[15])+fl[14]-1.0*fl[12]+fl[11]-1.0*fl[9]+fl[7]+fl[5]+fl[4]-1.0*fl[3]+fl[1]-1.0*fl[0]); 
    f1Quad[2] = 0.25*(fl[31]-1.0*fl[30]+fl[28]-1.0*(fl[27]+fl[26]+fl[24])+fl[22]-1.0*fl[20]+fl[19]-1.0*fl[17]+fl[16]+fl[13]+fl[10]-1.0*fl[8]+fl[6]-1.0*fl[2]); 
    limQuad[2] = fl[0]/cflL; 
  } else {
    f0Quad[2] = -0.25*(fr[29]-1.0*fr[25]+fr[23]-1.0*(fr[21]+fr[18]+fr[15])+fr[14]-1.0*fr[12]+fr[11]-1.0*fr[9]+fr[7]+fr[5]+fr[4]-1.0*fr[3]+fr[1]-1.0*fr[0]); 
    f1Quad[2] = -0.25*(fr[31]-1.0*fr[30]+fr[28]-1.0*(fr[27]+fr[26]+fr[24])+fr[22]-1.0*fr[20]+fr[19]-1.0*fr[17]+fr[16]+fr[13]+fr[10]-1.0*fr[8]+fr[6]-1.0*fr[2]); 
    limQuad[2] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[3] = 0.25*(fl[29]+fl[25]+fl[23]-1.0*(fl[21]+fl[18])+fl[15]-1.0*(fl[14]+fl[12]+fl[11]+fl[9])+fl[7]-1.0*(fl[5]+fl[4])+fl[3]+fl[1]+fl[0]); 
    f1Quad[3] = -0.25*(fl[31]+fl[30]+fl[28]-1.0*(fl[27]+fl[26])+fl[24]-1.0*(fl[22]+fl[20]+fl[19]+fl[17])+fl[16]-1.0*(fl[13]+fl[10])+fl[8]+fl[6]+fl[2]); 
    limQuad[3] = fl[0]/cflL; 
  } else {
    f0Quad[3] = 0.25*(fr[29]+fr[25]+fr[23]-1.0*(fr[21]+fr[18])+fr[15]-1.0*(fr[14]+fr[12]+fr[11]+fr[9])+fr[7]-1.0*(fr[5]+fr[4])+fr[3]+fr[1]+fr[0]); 
    f1Quad[3] = 0.25*(fr[31]+fr[30]+fr[28]-1.0*(fr[27]+fr[26])+fr[24]-1.0*(fr[22]+fr[20]+fr[19]+fr[17])+fr[16]-1.0*(fr[13]+fr[10])+fr[8]+fr[6]+fr[2]); 
    limQuad[3] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[4] = -0.25*(fl[29]-1.0*(fl[25]+fl[23])+fl[21]-1.0*fl[18]+fl[15]-1.0*(fl[14]+fl[12])+fl[11]+fl[9]-1.0*fl[7]+fl[5]-1.0*fl[4]+fl[3]+fl[1]-1.0*fl[0]); 
    f1Quad[4] = 0.25*(fl[31]-1.0*(fl[30]+fl[28])+fl[27]-1.0*fl[26]+fl[24]-1.0*(fl[22]+fl[20])+fl[19]+fl[17]-1.0*fl[16]+fl[13]-1.0*fl[10]+fl[8]+fl[6]-1.0*fl[2]); 
    limQuad[4] = fl[0]/cflL; 
  } else {
    f0Quad[4] = -0.25*(fr[29]-1.0*(fr[25]+fr[23])+fr[21]-1.0*fr[18]+fr[15]-1.0*(fr[14]+fr[12])+fr[11]+fr[9]-1.0*fr[7]+fr[5]-1.0*fr[4]+fr[3]+fr[1]-1.0*fr[0]); 
    f1Quad[4] = -0.25*(fr[31]-1.0*(fr[30]+fr[28])+fr[27]-1.0*fr[26]+fr[24]-1.0*(fr[22]+fr[20])+fr[19]+fr[17]-1.0*fr[16]+fr[13]-1.0*fr[10]+fr[8]+fr[6]-1.0*fr[2]); 
    limQuad[4] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[5] = 0.25*(fl[29]+fl[25]-1.0*fl[23]+fl[21]-1.0*(fl[18]+fl[15])+fl[14]-1.0*(fl[12]+fl[11])+fl[9]-1.0*(fl[7]+fl[5])+fl[4]-1.0*fl[3]+fl[1]+fl[0]); 
    f1Quad[5] = -0.25*(fl[31]+fl[30]-1.0*fl[28]+fl[27]-1.0*(fl[26]+fl[24])+fl[22]-1.0*(fl[20]+fl[19])+fl[17]-1.0*(fl[16]+fl[13])+fl[10]-1.0*fl[8]+fl[6]+fl[2]); 
    limQuad[5] = fl[0]/cflL; 
  } else {
    f0Quad[5] = 0.25*(fr[29]+fr[25]-1.0*fr[23]+fr[21]-1.0*(fr[18]+fr[15])+fr[14]-1.0*(fr[12]+fr[11])+fr[9]-1.0*(fr[7]+fr[5])+fr[4]-1.0*fr[3]+fr[1]+fr[0]); 
    f1Quad[5] = 0.25*(fr[31]+fr[30]-1.0*fr[28]+fr[27]-1.0*(fr[26]+fr[24])+fr[22]-1.0*(fr[20]+fr[19])+fr[17]-1.0*(fr[16]+fr[13])+fr[10]-1.0*fr[8]+fr[6]+fr[2]); 
    limQuad[5] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[6] = 0.25*(fl[29]-1.0*fl[25]+fl[23]+fl[21]-1.0*(fl[18]+fl[15]+fl[14])+fl[12]+fl[11]-1.0*(fl[9]+fl[7]+fl[5])+fl[4]+fl[3]-1.0*fl[1]+fl[0]); 
    f1Quad[6] = -0.25*(fl[31]-1.0*fl[30]+fl[28]+fl[27]-1.0*(fl[26]+fl[24]+fl[22])+fl[20]+fl[19]-1.0*(fl[17]+fl[16]+fl[13])+fl[10]+fl[8]-1.0*fl[6]+fl[2]); 
    limQuad[6] = fl[0]/cflL; 
  } else {
    f0Quad[6] = 0.25*(fr[29]-1.0*fr[25]+fr[23]+fr[21]-1.0*(fr[18]+fr[15]+fr[14])+fr[12]+fr[11]-1.0*(fr[9]+fr[7]+fr[5])+fr[4]+fr[3]-1.0*fr[1]+fr[0]); 
    f1Quad[6] = 0.25*(fr[31]-1.0*fr[30]+fr[28]+fr[27]-1.0*(fr[26]+fr[24]+fr[22])+fr[20]+fr[19]-1.0*(fr[17]+fr[16]+fr[13])+fr[10]+fr[8]-1.0*fr[6]+fr[2]); 
    limQuad[6] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[7] = -0.25*(fl[29]+fl[25]+fl[23]+fl[21]-1.0*fl[18]+fl[15]+fl[14]+fl[12]-1.0*(fl[11]+fl[9]+fl[7])+fl[5]-1.0*(fl[4]+fl[3]+fl[1]+fl[0])); 
    f1Quad[7] = 0.25*(fl[31]+fl[30]+fl[28]+fl[27]-1.0*fl[26]+fl[24]+fl[22]+fl[20]-1.0*(fl[19]+fl[17]+fl[16])+fl[13]-1.0*(fl[10]+fl[8]+fl[6]+fl[2])); 
    limQuad[7] = fl[0]/cflL; 
  } else {
    f0Quad[7] = -0.25*(fr[29]+fr[25]+fr[23]+fr[21]-1.0*fr[18]+fr[15]+fr[14]+fr[12]-1.0*(fr[11]+fr[9]+fr[7])+fr[5]-1.0*(fr[4]+fr[3]+fr[1]+fr[0])); 
    f1Quad[7] = -0.25*(fr[31]+fr[30]+fr[28]+fr[27]-1.0*fr[26]+fr[24]+fr[22]+fr[20]-1.0*(fr[19]+fr[17]+fr[16])+fr[13]-1.0*(fr[10]+fr[8]+fr[6]+fr[2])); 
    limQuad[7] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[8] = -0.25*(fl[29]-1.0*(fl[25]+fl[23]+fl[21])+fl[18]+fl[15]+fl[14]+fl[12]-1.0*(fl[11]+fl[9]+fl[7]+fl[5])+fl[4]+fl[3]+fl[1]-1.0*fl[0]); 
    f1Quad[8] = 0.25*(fl[31]-1.0*(fl[30]+fl[28]+fl[27])+fl[26]+fl[24]+fl[22]+fl[20]-1.0*(fl[19]+fl[17]+fl[16]+fl[13])+fl[10]+fl[8]+fl[6]-1.0*fl[2]); 
    limQuad[8] = fl[0]/cflL; 
  } else {
    f0Quad[8] = -0.25*(fr[29]-1.0*(fr[25]+fr[23]+fr[21])+fr[18]+fr[15]+fr[14]+fr[12]-1.0*(fr[11]+fr[9]+fr[7]+fr[5])+fr[4]+fr[3]+fr[1]-1.0*fr[0]); 
    f1Quad[8] = -0.25*(fr[31]-1.0*(fr[30]+fr[28]+fr[27])+fr[26]+fr[24]+fr[22]+fr[20]-1.0*(fr[19]+fr[17]+fr[16]+fr[13])+fr[10]+fr[8]+fr[6]-1.0*fr[2]); 
    limQuad[8] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[9] = 0.25*(fl[29]+fl[25]-1.0*(fl[23]+fl[21])+fl[18]-1.0*(fl[15]+fl[14])+fl[12]+fl[11]-1.0*(fl[9]+fl[7])+fl[5]-1.0*(fl[4]+fl[3])+fl[1]+fl[0]); 
    f1Quad[9] = -0.25*(fl[31]+fl[30]-1.0*(fl[28]+fl[27])+fl[26]-1.0*(fl[24]+fl[22])+fl[20]+fl[19]-1.0*(fl[17]+fl[16])+fl[13]-1.0*(fl[10]+fl[8])+fl[6]+fl[2]); 
    limQuad[9] = fl[0]/cflL; 
  } else {
    f0Quad[9] = 0.25*(fr[29]+fr[25]-1.0*(fr[23]+fr[21])+fr[18]-1.0*(fr[15]+fr[14])+fr[12]+fr[11]-1.0*(fr[9]+fr[7])+fr[5]-1.0*(fr[4]+fr[3])+fr[1]+fr[0]); 
    f1Quad[9] = 0.25*(fr[31]+fr[30]-1.0*(fr[28]+fr[27])+fr[26]-1.0*(fr[24]+fr[22])+fr[20]+fr[19]-1.0*(fr[17]+fr[16])+fr[13]-1.0*(fr[10]+fr[8])+fr[6]+fr[2]); 
    limQuad[9] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[10] = 0.25*(fl[29]-1.0*fl[25]+fl[23]-1.0*fl[21]+fl[18]-1.0*fl[15]+fl[14]-1.0*(fl[12]+fl[11])+fl[9]-1.0*fl[7]+fl[5]-1.0*fl[4]+fl[3]-1.0*fl[1]+fl[0]); 
    f1Quad[10] = -0.25*(fl[31]-1.0*fl[30]+fl[28]-1.0*fl[27]+fl[26]-1.0*fl[24]+fl[22]-1.0*(fl[20]+fl[19])+fl[17]-1.0*fl[16]+fl[13]-1.0*fl[10]+fl[8]-1.0*fl[6]+fl[2]); 
    limQuad[10] = fl[0]/cflL; 
  } else {
    f0Quad[10] = 0.25*(fr[29]-1.0*fr[25]+fr[23]-1.0*fr[21]+fr[18]-1.0*fr[15]+fr[14]-1.0*(fr[12]+fr[11])+fr[9]-1.0*fr[7]+fr[5]-1.0*fr[4]+fr[3]-1.0*fr[1]+fr[0]); 
    f1Quad[10] = 0.25*(fr[31]-1.0*fr[30]+fr[28]-1.0*fr[27]+fr[26]-1.0*fr[24]+fr[22]-1.0*(fr[20]+fr[19])+fr[17]-1.0*fr[16]+fr[13]-1.0*fr[10]+fr[8]-1.0*fr[6]+fr[2]); 
    limQuad[10] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[11] = -0.25*(fl[29]+fl[25]+fl[23]-1.0*fl[21]+fl[18]+fl[15]-1.0*(fl[14]+fl[12])+fl[11]+fl[9]-1.0*(fl[7]+fl[5])+fl[4]-1.0*(fl[3]+fl[1]+fl[0])); 
    f1Quad[11] = 0.25*(fl[31]+fl[30]+fl[28]-1.0*fl[27]+fl[26]+fl[24]-1.0*(fl[22]+fl[20])+fl[19]+fl[17]-1.0*(fl[16]+fl[13])+fl[10]-1.0*(fl[8]+fl[6]+fl[2])); 
    limQuad[11] = fl[0]/cflL; 
  } else {
    f0Quad[11] = -0.25*(fr[29]+fr[25]+fr[23]-1.0*fr[21]+fr[18]+fr[15]-1.0*(fr[14]+fr[12])+fr[11]+fr[9]-1.0*(fr[7]+fr[5])+fr[4]-1.0*(fr[3]+fr[1]+fr[0])); 
    f1Quad[11] = -0.25*(fr[31]+fr[30]+fr[28]-1.0*fr[27]+fr[26]+fr[24]-1.0*(fr[22]+fr[20])+fr[19]+fr[17]-1.0*(fr[16]+fr[13])+fr[10]-1.0*(fr[8]+fr[6]+fr[2])); 
    limQuad[11] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[12] = 0.25*(fl[29]-1.0*(fl[25]+fl[23])+fl[21]+fl[18]+fl[15]-1.0*(fl[14]+fl[12]+fl[11]+fl[9])+fl[7]+fl[5]+fl[4]-1.0*(fl[3]+fl[1])+fl[0]); 
    f1Quad[12] = -0.25*(fl[31]-1.0*(fl[30]+fl[28])+fl[27]+fl[26]+fl[24]-1.0*(fl[22]+fl[20]+fl[19]+fl[17])+fl[16]+fl[13]+fl[10]-1.0*(fl[8]+fl[6])+fl[2]); 
    limQuad[12] = fl[0]/cflL; 
  } else {
    f0Quad[12] = 0.25*(fr[29]-1.0*(fr[25]+fr[23])+fr[21]+fr[18]+fr[15]-1.0*(fr[14]+fr[12]+fr[11]+fr[9])+fr[7]+fr[5]+fr[4]-1.0*(fr[3]+fr[1])+fr[0]); 
    f1Quad[12] = 0.25*(fr[31]-1.0*(fr[30]+fr[28])+fr[27]+fr[26]+fr[24]-1.0*(fr[22]+fr[20]+fr[19]+fr[17])+fr[16]+fr[13]+fr[10]-1.0*(fr[8]+fr[6])+fr[2]); 
    limQuad[12] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[13] = -0.25*(fl[29]+fl[25]-1.0*fl[23]+fl[21]+fl[18]-1.0*fl[15]+fl[14]-1.0*fl[12]+fl[11]-1.0*fl[9]+fl[7]-1.0*(fl[5]+fl[4])+fl[3]-1.0*(fl[1]+fl[0])); 
    f1Quad[13] = 0.25*(fl[31]+fl[30]-1.0*fl[28]+fl[27]+fl[26]-1.0*fl[24]+fl[22]-1.0*fl[20]+fl[19]-1.0*fl[17]+fl[16]-1.0*(fl[13]+fl[10])+fl[8]-1.0*(fl[6]+fl[2])); 
    limQuad[13] = fl[0]/cflL; 
  } else {
    f0Quad[13] = -0.25*(fr[29]+fr[25]-1.0*fr[23]+fr[21]+fr[18]-1.0*fr[15]+fr[14]-1.0*fr[12]+fr[11]-1.0*fr[9]+fr[7]-1.0*(fr[5]+fr[4])+fr[3]-1.0*(fr[1]+fr[0])); 
    f1Quad[13] = -0.25*(fr[31]+fr[30]-1.0*fr[28]+fr[27]+fr[26]-1.0*fr[24]+fr[22]-1.0*fr[20]+fr[19]-1.0*fr[17]+fr[16]-1.0*(fr[13]+fr[10])+fr[8]-1.0*(fr[6]+fr[2])); 
    limQuad[13] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[14] = -0.25*(fl[29]-1.0*fl[25]+fl[23]+fl[21]+fl[18]-1.0*(fl[15]+fl[14])+fl[12]-1.0*fl[11]+fl[9]+fl[7]-1.0*(fl[5]+fl[4]+fl[3])+fl[1]-1.0*fl[0]); 
    f1Quad[14] = 0.25*(fl[31]-1.0*fl[30]+fl[28]+fl[27]+fl[26]-1.0*(fl[24]+fl[22])+fl[20]-1.0*fl[19]+fl[17]+fl[16]-1.0*(fl[13]+fl[10]+fl[8])+fl[6]-1.0*fl[2]); 
    limQuad[14] = fl[0]/cflL; 
  } else {
    f0Quad[14] = -0.25*(fr[29]-1.0*fr[25]+fr[23]+fr[21]+fr[18]-1.0*(fr[15]+fr[14])+fr[12]-1.0*fr[11]+fr[9]+fr[7]-1.0*(fr[5]+fr[4]+fr[3])+fr[1]-1.0*fr[0]); 
    f1Quad[14] = -0.25*(fr[31]-1.0*fr[30]+fr[28]+fr[27]+fr[26]-1.0*(fr[24]+fr[22])+fr[20]-1.0*fr[19]+fr[17]+fr[16]-1.0*(fr[13]+fr[10]+fr[8])+fr[6]-1.0*fr[2]); 
    limQuad[14] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[15] = 0.25*(fl[29]+fl[25]+fl[23]+fl[21]+fl[18]+fl[15]+fl[14]+fl[12]+fl[11]+fl[9]+fl[7]+fl[5]+fl[4]+fl[3]+fl[1]+fl[0]); 
    f1Quad[15] = -0.25*(fl[31]+fl[30]+fl[28]+fl[27]+fl[26]+fl[24]+fl[22]+fl[20]+fl[19]+fl[17]+fl[16]+fl[13]+fl[10]+fl[8]+fl[6]+fl[2]); 
    limQuad[15] = fl[0]/cflL; 
  } else {
    f0Quad[15] = 0.25*(fr[29]+fr[25]+fr[23]+fr[21]+fr[18]+fr[15]+fr[14]+fr[12]+fr[11]+fr[9]+fr[7]+fr[5]+fr[4]+fr[3]+fr[1]+fr[0]); 
    f1Quad[15] = 0.25*(fr[31]+fr[30]+fr[28]+fr[27]+fr[26]+fr[24]+fr[22]+fr[20]+fr[19]+fr[17]+fr[16]+fr[13]+fr[10]+fr[8]+fr[6]+fr[2]); 
    limQuad[15] = fr[0]/cflR; 
  }
  double fhat[32]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[3] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[6] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[7] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[8] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[9] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[10] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[11] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[12] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[13] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[14] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[15] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[16] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[17] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[18] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[19] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[20] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[21] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[22] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[23] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[24] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[25] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[26] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[27] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[28] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[29] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[30] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[31] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[16];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]+fhat[22]+fhat[20]+fhat[19]+fhat[17]+fhat[16]+3.0*fhat[2]))-3.0*(fhat[30]+fhat[28]+fhat[27]+fhat[26]+3.0*(fhat[13]+fhat[10]+fhat[8]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[29]-1.732050807568877*(fhat[25]+fhat[23]+fhat[21]+fhat[18]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[15]+fhat[14]+fhat[12]+fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0]))); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[24]+fhat[22]))+fhat[20]-1.0*fhat[19]+fhat[17]+fhat[16]-3.0*fhat[2]))+3.0*(fhat[30]-1.0*(fhat[28]+fhat[27]+fhat[26]-3.0*(fhat[13]+fhat[10]+fhat[8]-1.0*fhat[6])))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[15]+fhat[14]-1.0*fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0]))-1.0*(fhat[29]+1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21]+fhat[18])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[1]))))); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[24])+fhat[22]-1.0*fhat[20]+fhat[19]-1.0*fhat[17]+fhat[16]-3.0*fhat[2]))+3.0*((-1.0*fhat[30])+fhat[28]-1.0*(fhat[27]+fhat[26]-3.0*(fhat[13]+fhat[10]-1.0*fhat[8]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[29])+1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]+fhat[18]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[1]))+3.0*(fhat[15]-1.0*fhat[14]+fhat[12]-1.0*fhat[11]+fhat[9]-1.0*fhat[7]+3.0*fhat[0]))); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]-1.0*(fhat[22]+fhat[20]+fhat[19]+fhat[17]-1.0*(fhat[16]+3.0*fhat[2]))))+3.0*(fhat[30]+fhat[28]-1.0*(fhat[27]+fhat[26]-3.0*((-1.0*(fhat[13]+fhat[10]))+fhat[8]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]-1.0*(fhat[21]+fhat[18])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]+fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0]))))); 
  rCtrl[4] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]-1.0*(fhat[22]+fhat[20]-1.0*fhat[19])+fhat[17]-1.0*(fhat[16]+3.0*fhat[2])))+3.0*((-1.0*(fhat[30]+fhat[28]))+fhat[27]-1.0*fhat[26]+3.0*(fhat[13]-1.0*fhat[10]+fhat[8]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[29])+1.732050807568877*(fhat[25]+fhat[23]-1.0*fhat[21]+fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[1])))+3.0*((-1.0*fhat[15])+fhat[14]+fhat[12]-1.0*(fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0]))))); 
  rCtrl[5] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[24])+fhat[22]-1.0*(fhat[20]+fhat[19]-1.0*fhat[17]+fhat[16])+3.0*fhat[2]))+3.0*(fhat[30]-1.0*fhat[28]+fhat[27]-1.0*fhat[26]+3.0*((-1.0*fhat[13])+fhat[10]-1.0*fhat[8]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[29]+1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]-1.0*fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7])+3.0*fhat[0]))); 
  rCtrl[6] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[24]+fhat[22]))+fhat[20]+fhat[19]-1.0*(fhat[17]+fhat[16]-3.0*fhat[2])))+3.0*((-1.0*fhat[30])+fhat[28]+fhat[27]-1.0*fhat[26]+3.0*((-1.0*fhat[13])+fhat[10]+fhat[8]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[29]-1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21])+fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0])))); 
  rCtrl[7] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]+fhat[22]+fhat[20]-1.0*(fhat[19]+fhat[17]+fhat[16]+3.0*fhat[2])))+3.0*(fhat[30]+fhat[28]+fhat[27]-1.0*fhat[26]+3.0*(fhat[13]-1.0*(fhat[10]+fhat[8]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*(fhat[15]+fhat[14]+fhat[12]-1.0*fhat[11]))+fhat[9]+fhat[7]+3.0*fhat[0])-1.0*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]+fhat[21]-1.0*fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[1])))))); 
  rCtrl[8] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]+fhat[22]+fhat[20]-1.0*(fhat[19]+fhat[17]+fhat[16]+3.0*fhat[2])))+3.0*(3.0*((-1.0*fhat[13])+fhat[10]+fhat[8]+fhat[6])-1.0*(fhat[30]+fhat[28]+fhat[27]-1.0*fhat[26]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[29])+1.732050807568877*(fhat[25]+fhat[23]+fhat[21]-1.0*fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[14]+fhat[12]-1.0*fhat[11]))+fhat[9]+fhat[7]+3.0*fhat[0]))); 
  rCtrl[9] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[24]+fhat[22]))+fhat[20]+fhat[19]-1.0*(fhat[17]+fhat[16]-3.0*fhat[2])))+3.0*(fhat[30]-1.0*(fhat[28]+fhat[27]-1.0*fhat[26])+3.0*(fhat[13]-1.0*(fhat[10]+fhat[8]-1.0*fhat[6])))))/(72.0*EPSILON+1.414213562373095*(fhat[29]+1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21])+fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0])))); 
  rCtrl[10] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[24])+fhat[22]-1.0*(fhat[20]+fhat[19]-1.0*fhat[17]+fhat[16])+3.0*fhat[2]))+3.0*((-1.0*fhat[30])+fhat[28]-1.0*fhat[27]+fhat[26]+3.0*(fhat[13]-1.0*fhat[10]+fhat[8]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[29]-1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]-1.0*fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7])+3.0*fhat[0]))); 
  rCtrl[11] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]-1.0*(fhat[22]+fhat[20]-1.0*fhat[19])+fhat[17]-1.0*(fhat[16]+3.0*fhat[2])))+3.0*(fhat[30]+fhat[28]-1.0*fhat[27]+fhat[26]+3.0*((-1.0*fhat[13])+fhat[10]-1.0*(fhat[8]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*fhat[15])+fhat[14]+fhat[12]-1.0*(fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0])))-1.0*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]-1.0*fhat[21]+fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[1])))))); 
  rCtrl[12] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]-1.0*(fhat[22]+fhat[20]+fhat[19]+fhat[17]-1.0*(fhat[16]+3.0*fhat[2]))))+3.0*((-1.0*(fhat[30]+fhat[28]))+fhat[27]+fhat[26]+3.0*(fhat[13]+fhat[10]-1.0*(fhat[8]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(fhat[29]-1.732050807568877*(fhat[25]+fhat[23]-1.0*(fhat[21]+fhat[18])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]+fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0]))))); 
  rCtrl[13] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[24])+fhat[22]-1.0*fhat[20]+fhat[19]-1.0*fhat[17]+fhat[16]-3.0*fhat[2]))+3.0*(fhat[30]-1.0*fhat[28]+fhat[27]+fhat[26]+3.0*((-1.0*(fhat[13]+fhat[10]))+fhat[8]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[15]-1.0*fhat[14]+fhat[12]-1.0*fhat[11]+fhat[9]-1.0*fhat[7]+3.0*fhat[0])-1.0*(fhat[29]+1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]+fhat[18]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[1]))))); 
  rCtrl[14] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[24]+fhat[22]))+fhat[20]-1.0*fhat[19]+fhat[17]+fhat[16]-3.0*fhat[2]))+3.0*((-1.0*fhat[30])+fhat[28]+fhat[27]+fhat[26]-3.0*(fhat[13]+fhat[10]+fhat[8]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[29])+1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21]+fhat[18])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[1]))+3.0*(fhat[15]+fhat[14]-1.0*fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0])))); 
  rCtrl[15] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]+fhat[22]+fhat[20]+fhat[19]+fhat[17]+fhat[16]+3.0*fhat[2]))+3.0*(fhat[30]+fhat[28]+fhat[27]+fhat[26]+3.0*(fhat[13]+fhat[10]+fhat[8]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]+fhat[21]+fhat[18]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[15]+fhat[14]+fhat[12]+fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0]))); 
  double fhatCtrl[16];  // fhatCtrl = anti-limited fhat evaluated at each control node on y surface 
  fhatCtrl[0] = -0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[23]+fhat[21]+fhat[18]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[1]))-1.0*(fhat[29]+3.0*(fhat[15]+fhat[14]+fhat[12]+fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0])))*limTheta(rCtrl[0],-1.0,EPSILON); 
  fhatCtrl[1] = -0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21]+fhat[18])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[1]))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]-1.0*fhat[11]+fhat[9]+fhat[7]-3.0*fhat[0]))*limTheta(rCtrl[1],-1.0,EPSILON); 
  fhatCtrl[2] = 0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]+fhat[18]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[1]))-1.0*(fhat[29]+3.0*((-1.0*fhat[15])+fhat[14]-1.0*fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7]-3.0*fhat[0])))*limTheta(rCtrl[2],-1.0,EPSILON); 
  fhatCtrl[3] = 0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]-1.0*(fhat[21]+fhat[18])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]+fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0]))))*limTheta(rCtrl[3],-1.0,EPSILON); 
  fhatCtrl[4] = 0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[23]-1.0*fhat[21]+fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[1])))-1.0*(fhat[29]+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]-1.0*fhat[11])+fhat[9]-1.0*(fhat[7]+3.0*fhat[0]))))*limTheta(rCtrl[4],-1.0,EPSILON); 
  fhatCtrl[5] = 0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]-1.0*fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7])+3.0*fhat[0]))*limTheta(rCtrl[5],-1.0,EPSILON); 
  fhatCtrl[6] = -0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21])+fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[1])))-1.0*(fhat[29]+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0]))))*limTheta(rCtrl[6],-1.0,EPSILON); 
  fhatCtrl[7] = -0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]+fhat[21]-1.0*fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[1])))+3.0*(fhat[15]+fhat[14]+fhat[12]-1.0*(fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0])))*limTheta(rCtrl[7],-1.0,EPSILON); 
  fhatCtrl[8] = 0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[23]+fhat[21]-1.0*fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[1])))-1.0*(fhat[29]+3.0*(fhat[15]+fhat[14]+fhat[12]-1.0*(fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0]))))*limTheta(rCtrl[8],-1.0,EPSILON); 
  fhatCtrl[9] = 0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21])+fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0])))*limTheta(rCtrl[9],-1.0,EPSILON); 
  fhatCtrl[10] = -0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]-1.0*fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[1]))-1.0*(fhat[29]+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7])+3.0*fhat[0])))*limTheta(rCtrl[10],-1.0,EPSILON); 
  fhatCtrl[11] = -0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]-1.0*fhat[21]+fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[1])))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]-1.0*fhat[11])+fhat[9]-1.0*(fhat[7]+3.0*fhat[0])))*limTheta(rCtrl[11],-1.0,EPSILON); 
  fhatCtrl[12] = -0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[23]-1.0*(fhat[21]+fhat[18])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[1]))-1.0*(fhat[29]+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]+fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0])))))*limTheta(rCtrl[12],-1.0,EPSILON); 
  fhatCtrl[13] = -0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]+fhat[18]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[1]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7]-3.0*fhat[0]))*limTheta(rCtrl[13],-1.0,EPSILON); 
  fhatCtrl[14] = 0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21]+fhat[18])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[1]))-1.0*(fhat[29]+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]-1.0*fhat[11]+fhat[9]+fhat[7]-3.0*fhat[0])))*limTheta(rCtrl[14],-1.0,EPSILON); 
  fhatCtrl[15] = 0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]+fhat[21]+fhat[18]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[15]+fhat[14]+fhat[12]+fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0]))*limTheta(rCtrl[15],-1.0,EPSILON); 
  double fhatAL[16];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.25*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.4330127018922193*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[5] = 0.75*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[6] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[8] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[9] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[10] = 0.75*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[11] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[12] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[13] = 1.299038105676658*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[14] = 1.299038105676658*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[15] = 2.25*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[16]; 
  fhatALQuad[0] = fmin(0.25*((-0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7])))-1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[0]); 
  fhatALQuad[1] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*((-0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+0.5773502691896258*(fhatAL[13]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[1]); 
  fhatALQuad[2] = fmin(0.25*(0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[2]); 
  fhatALQuad[3] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[3]); 
  fhatALQuad[4] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*(0.5773502691896258*fhatAL[3]-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[4]+fhatAL[2]+fhatAL[1]))-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[4]); 
  fhatALQuad[5] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]))+0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[5]); 
  fhatALQuad[6] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2])-0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[6]); 
  fhatALQuad[7] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[7]); 
  fhatALQuad[8] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]+fhatAL[1]))+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[8]); 
  fhatALQuad[9] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]))+0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[9]); 
  fhatALQuad[10] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[10]); 
  fhatALQuad[11] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[4]+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[11]); 
  fhatALQuad[12] = fmin(0.25*((-0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5])))+1.732050807568877*(0.5773502691896258*(fhatAL[4]+fhatAL[3])-0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[12]); 
  fhatALQuad[13] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*(0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1])-0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[13]); 
  fhatALQuad[14] = fmin(0.25*(0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[1])+0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[14]); 
  fhatALQuad[15] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[15]); 
  fhatAL[0] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[5] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[6] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[8] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[9] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[10] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[11] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[12] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[13] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[14] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[15] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.1767766952966368*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_y; 
  incr[1] = 0.1767766952966368*(alpha[3]*fhatAL[6]+alpha[2]*fhatAL[5]+alpha[0]*fhatAL[1])*dfac_y; 
  incr[2] = -0.3061862178478971*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_y; 
  incr[3] = 0.1767766952966368*(alpha[3]*fhatAL[7]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_y; 
  incr[4] = 0.1767766952966368*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_y; 
  incr[5] = 0.1767766952966368*(alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[0]*fhatAL[4])*dfac_y; 
  incr[6] = -0.3061862178478971*(alpha[3]*fhatAL[6]+alpha[2]*fhatAL[5]+alpha[0]*fhatAL[1])*dfac_y; 
  incr[7] = 0.1767766952966368*(alpha[3]*fhatAL[11]+alpha[0]*fhatAL[5]+fhatAL[1]*alpha[2])*dfac_y; 
  incr[8] = -0.3061862178478971*(alpha[3]*fhatAL[7]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_y; 
  incr[9] = 0.1767766952966368*(alpha[2]*fhatAL[11]+alpha[0]*fhatAL[6]+fhatAL[1]*alpha[3])*dfac_y; 
  incr[10] = -0.3061862178478971*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_y; 
  incr[11] = 0.1767766952966368*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_y; 
  incr[12] = 0.1767766952966368*(alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[0]*fhatAL[8])*dfac_y; 
  incr[13] = -0.3061862178478971*(alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[0]*fhatAL[4])*dfac_y; 
  incr[14] = 0.1767766952966368*(alpha[3]*fhatAL[14]+alpha[0]*fhatAL[9]+alpha[2]*fhatAL[4])*dfac_y; 
  incr[15] = 0.1767766952966368*(alpha[2]*fhatAL[14]+alpha[0]*fhatAL[10]+alpha[3]*fhatAL[4])*dfac_y; 
  incr[16] = -0.3061862178478971*(alpha[3]*fhatAL[11]+alpha[0]*fhatAL[5]+fhatAL[1]*alpha[2])*dfac_y; 
  incr[17] = -0.3061862178478971*(alpha[2]*fhatAL[11]+alpha[0]*fhatAL[6]+fhatAL[1]*alpha[3])*dfac_y; 
  incr[18] = 0.1767766952966368*(alpha[0]*fhatAL[11]+alpha[2]*fhatAL[6]+alpha[3]*fhatAL[5])*dfac_y; 
  incr[19] = -0.3061862178478971*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_y; 
  incr[20] = -0.3061862178478971*(alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[0]*fhatAL[8])*dfac_y; 
  incr[21] = 0.1767766952966368*(alpha[3]*fhatAL[15]+alpha[0]*fhatAL[12]+alpha[2]*fhatAL[8])*dfac_y; 
  incr[22] = -0.3061862178478971*(alpha[3]*fhatAL[14]+alpha[0]*fhatAL[9]+alpha[2]*fhatAL[4])*dfac_y; 
  incr[23] = 0.1767766952966368*(alpha[2]*fhatAL[15]+alpha[0]*fhatAL[13]+alpha[3]*fhatAL[8])*dfac_y; 
  incr[24] = -0.3061862178478971*(alpha[2]*fhatAL[14]+alpha[0]*fhatAL[10]+alpha[3]*fhatAL[4])*dfac_y; 
  incr[25] = 0.1767766952966368*(alpha[0]*fhatAL[14]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9])*dfac_y; 
  incr[26] = -0.3061862178478971*(alpha[0]*fhatAL[11]+alpha[2]*fhatAL[6]+alpha[3]*fhatAL[5])*dfac_y; 
  incr[27] = -0.3061862178478971*(alpha[3]*fhatAL[15]+alpha[0]*fhatAL[12]+alpha[2]*fhatAL[8])*dfac_y; 
  incr[28] = -0.3061862178478971*(alpha[2]*fhatAL[15]+alpha[0]*fhatAL[13]+alpha[3]*fhatAL[8])*dfac_y; 
  incr[29] = 0.1767766952966368*(alpha[0]*fhatAL[15]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12])*dfac_y; 
  incr[30] = -0.3061862178478971*(alpha[0]*fhatAL[14]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9])*dfac_y; 
  incr[31] = -0.3061862178478971*(alpha[0]*fhatAL[15]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 
  outr[28] += incr[28]; 
  outr[29] += incr[29]; 
  outr[30] += incr[30]; 
  outr[31] += incr[31]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += incr[19]; 
  outl[20] += incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += incr[22]; 
  outl[23] += -1.0*incr[23]; 
  outl[24] += incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += incr[26]; 
  outl[27] += incr[27]; 
  outl[28] += incr[28]; 
  outl[29] += -1.0*incr[29]; 
  outl[30] += incr[30]; 
  outl[31] += incr[31]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity3x2vSer_Z_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_z = 2.0/dxv[2]; 
  double dfac_v = 2.0/dxv[3]; 
  double dfac_m = 2.0/dxv[4]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wz = w[2]; 
  double wv = w[3]; 
  double wm = w[4]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[32]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.1767766952966368*Gradpar[0]*wv; 

  double alpha[16]; 
  alpha[0] = 1.414213562373095*Gradpar[0]*wv; 
  double f0Quad[16]; 
  double f1Quad[16]; 
  double limQuad[16]; 
  // determine upwinding at each surface quadrature node 
  if(0.25*alpha[0] > 0) {
    f0Quad[0] = 0.25*(fl[28]-1.0*(fl[24]+fl[23]+fl[20]+fl[17])+fl[15]+fl[13]+fl[12]+fl[10]+fl[9]+fl[6]-1.0*(fl[5]+fl[4]+fl[2]+fl[1])+fl[0]); 
    f1Quad[0] = -0.25*(fl[31]-1.0*(fl[30]+fl[29]+fl[27]+fl[26])+fl[25]+fl[22]+fl[21]+fl[19]+fl[18]+fl[16]-1.0*(fl[14]+fl[11]+fl[8]+fl[7])+fl[3]); 
    limQuad[0] = fl[0]/cflL; 
  } else {
    f0Quad[0] = 0.25*(fr[28]-1.0*(fr[24]+fr[23]+fr[20]+fr[17])+fr[15]+fr[13]+fr[12]+fr[10]+fr[9]+fr[6]-1.0*(fr[5]+fr[4]+fr[2]+fr[1])+fr[0]); 
    f1Quad[0] = 0.25*(fr[31]-1.0*(fr[30]+fr[29]+fr[27]+fr[26])+fr[25]+fr[22]+fr[21]+fr[19]+fr[18]+fr[16]-1.0*(fr[14]+fr[11]+fr[8]+fr[7])+fr[3]); 
    limQuad[0] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[1] = -0.25*(fl[28]+fl[24]-1.0*(fl[23]+fl[20]+fl[17]+fl[15]+fl[13])+fl[12]-1.0*fl[10]+fl[9]+fl[6]+fl[5]+fl[4]+fl[2]-1.0*(fl[1]+fl[0])); 
    f1Quad[1] = 0.25*(fl[31]+fl[30]-1.0*(fl[29]+fl[27]+fl[26]+fl[25]+fl[22])+fl[21]-1.0*fl[19]+fl[18]+fl[16]+fl[14]+fl[11]+fl[8]-1.0*(fl[7]+fl[3])); 
    limQuad[1] = fl[0]/cflL; 
  } else {
    f0Quad[1] = -0.25*(fr[28]+fr[24]-1.0*(fr[23]+fr[20]+fr[17]+fr[15]+fr[13])+fr[12]-1.0*fr[10]+fr[9]+fr[6]+fr[5]+fr[4]+fr[2]-1.0*(fr[1]+fr[0])); 
    f1Quad[1] = -0.25*(fr[31]+fr[30]-1.0*(fr[29]+fr[27]+fr[26]+fr[25]+fr[22])+fr[21]-1.0*fr[19]+fr[18]+fr[16]+fr[14]+fr[11]+fr[8]-1.0*(fr[7]+fr[3])); 
    limQuad[1] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[2] = -0.25*(fl[28]-1.0*fl[24]+fl[23]-1.0*(fl[20]+fl[17]+fl[15])+fl[13]-1.0*fl[12]+fl[10]-1.0*fl[9]+fl[6]+fl[5]+fl[4]-1.0*fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[2] = 0.25*(fl[31]-1.0*fl[30]+fl[29]-1.0*(fl[27]+fl[26]+fl[25])+fl[22]-1.0*fl[21]+fl[19]-1.0*fl[18]+fl[16]+fl[14]+fl[11]-1.0*fl[8]+fl[7]-1.0*fl[3]); 
    limQuad[2] = fl[0]/cflL; 
  } else {
    f0Quad[2] = -0.25*(fr[28]-1.0*fr[24]+fr[23]-1.0*(fr[20]+fr[17]+fr[15])+fr[13]-1.0*fr[12]+fr[10]-1.0*fr[9]+fr[6]+fr[5]+fr[4]-1.0*fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[2] = -0.25*(fr[31]-1.0*fr[30]+fr[29]-1.0*(fr[27]+fr[26]+fr[25])+fr[22]-1.0*fr[21]+fr[19]-1.0*fr[18]+fr[16]+fr[14]+fr[11]-1.0*fr[8]+fr[7]-1.0*fr[3]); 
    limQuad[2] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[3] = 0.25*(fl[28]+fl[24]+fl[23]-1.0*(fl[20]+fl[17])+fl[15]-1.0*(fl[13]+fl[12]+fl[10]+fl[9])+fl[6]-1.0*(fl[5]+fl[4])+fl[2]+fl[1]+fl[0]); 
    f1Quad[3] = -0.25*(fl[31]+fl[30]+fl[29]-1.0*(fl[27]+fl[26])+fl[25]-1.0*(fl[22]+fl[21]+fl[19]+fl[18])+fl[16]-1.0*(fl[14]+fl[11])+fl[8]+fl[7]+fl[3]); 
    limQuad[3] = fl[0]/cflL; 
  } else {
    f0Quad[3] = 0.25*(fr[28]+fr[24]+fr[23]-1.0*(fr[20]+fr[17])+fr[15]-1.0*(fr[13]+fr[12]+fr[10]+fr[9])+fr[6]-1.0*(fr[5]+fr[4])+fr[2]+fr[1]+fr[0]); 
    f1Quad[3] = 0.25*(fr[31]+fr[30]+fr[29]-1.0*(fr[27]+fr[26])+fr[25]-1.0*(fr[22]+fr[21]+fr[19]+fr[18])+fr[16]-1.0*(fr[14]+fr[11])+fr[8]+fr[7]+fr[3]); 
    limQuad[3] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[4] = -0.25*(fl[28]-1.0*(fl[24]+fl[23])+fl[20]-1.0*fl[17]+fl[15]-1.0*(fl[13]+fl[12])+fl[10]+fl[9]-1.0*fl[6]+fl[5]-1.0*fl[4]+fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[4] = 0.25*(fl[31]-1.0*(fl[30]+fl[29])+fl[27]-1.0*fl[26]+fl[25]-1.0*(fl[22]+fl[21])+fl[19]+fl[18]-1.0*fl[16]+fl[14]-1.0*fl[11]+fl[8]+fl[7]-1.0*fl[3]); 
    limQuad[4] = fl[0]/cflL; 
  } else {
    f0Quad[4] = -0.25*(fr[28]-1.0*(fr[24]+fr[23])+fr[20]-1.0*fr[17]+fr[15]-1.0*(fr[13]+fr[12])+fr[10]+fr[9]-1.0*fr[6]+fr[5]-1.0*fr[4]+fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[4] = -0.25*(fr[31]-1.0*(fr[30]+fr[29])+fr[27]-1.0*fr[26]+fr[25]-1.0*(fr[22]+fr[21])+fr[19]+fr[18]-1.0*fr[16]+fr[14]-1.0*fr[11]+fr[8]+fr[7]-1.0*fr[3]); 
    limQuad[4] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[5] = 0.25*(fl[28]+fl[24]-1.0*fl[23]+fl[20]-1.0*(fl[17]+fl[15])+fl[13]-1.0*(fl[12]+fl[10])+fl[9]-1.0*(fl[6]+fl[5])+fl[4]-1.0*fl[2]+fl[1]+fl[0]); 
    f1Quad[5] = -0.25*(fl[31]+fl[30]-1.0*fl[29]+fl[27]-1.0*(fl[26]+fl[25])+fl[22]-1.0*(fl[21]+fl[19])+fl[18]-1.0*(fl[16]+fl[14])+fl[11]-1.0*fl[8]+fl[7]+fl[3]); 
    limQuad[5] = fl[0]/cflL; 
  } else {
    f0Quad[5] = 0.25*(fr[28]+fr[24]-1.0*fr[23]+fr[20]-1.0*(fr[17]+fr[15])+fr[13]-1.0*(fr[12]+fr[10])+fr[9]-1.0*(fr[6]+fr[5])+fr[4]-1.0*fr[2]+fr[1]+fr[0]); 
    f1Quad[5] = 0.25*(fr[31]+fr[30]-1.0*fr[29]+fr[27]-1.0*(fr[26]+fr[25])+fr[22]-1.0*(fr[21]+fr[19])+fr[18]-1.0*(fr[16]+fr[14])+fr[11]-1.0*fr[8]+fr[7]+fr[3]); 
    limQuad[5] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[6] = 0.25*(fl[28]-1.0*fl[24]+fl[23]+fl[20]-1.0*(fl[17]+fl[15]+fl[13])+fl[12]+fl[10]-1.0*(fl[9]+fl[6]+fl[5])+fl[4]+fl[2]-1.0*fl[1]+fl[0]); 
    f1Quad[6] = -0.25*(fl[31]-1.0*fl[30]+fl[29]+fl[27]-1.0*(fl[26]+fl[25]+fl[22])+fl[21]+fl[19]-1.0*(fl[18]+fl[16]+fl[14])+fl[11]+fl[8]-1.0*fl[7]+fl[3]); 
    limQuad[6] = fl[0]/cflL; 
  } else {
    f0Quad[6] = 0.25*(fr[28]-1.0*fr[24]+fr[23]+fr[20]-1.0*(fr[17]+fr[15]+fr[13])+fr[12]+fr[10]-1.0*(fr[9]+fr[6]+fr[5])+fr[4]+fr[2]-1.0*fr[1]+fr[0]); 
    f1Quad[6] = 0.25*(fr[31]-1.0*fr[30]+fr[29]+fr[27]-1.0*(fr[26]+fr[25]+fr[22])+fr[21]+fr[19]-1.0*(fr[18]+fr[16]+fr[14])+fr[11]+fr[8]-1.0*fr[7]+fr[3]); 
    limQuad[6] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[7] = -0.25*(fl[28]+fl[24]+fl[23]+fl[20]-1.0*fl[17]+fl[15]+fl[13]+fl[12]-1.0*(fl[10]+fl[9]+fl[6])+fl[5]-1.0*(fl[4]+fl[2]+fl[1]+fl[0])); 
    f1Quad[7] = 0.25*(fl[31]+fl[30]+fl[29]+fl[27]-1.0*fl[26]+fl[25]+fl[22]+fl[21]-1.0*(fl[19]+fl[18]+fl[16])+fl[14]-1.0*(fl[11]+fl[8]+fl[7]+fl[3])); 
    limQuad[7] = fl[0]/cflL; 
  } else {
    f0Quad[7] = -0.25*(fr[28]+fr[24]+fr[23]+fr[20]-1.0*fr[17]+fr[15]+fr[13]+fr[12]-1.0*(fr[10]+fr[9]+fr[6])+fr[5]-1.0*(fr[4]+fr[2]+fr[1]+fr[0])); 
    f1Quad[7] = -0.25*(fr[31]+fr[30]+fr[29]+fr[27]-1.0*fr[26]+fr[25]+fr[22]+fr[21]-1.0*(fr[19]+fr[18]+fr[16])+fr[14]-1.0*(fr[11]+fr[8]+fr[7]+fr[3])); 
    limQuad[7] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[8] = -0.25*(fl[28]-1.0*(fl[24]+fl[23]+fl[20])+fl[17]+fl[15]+fl[13]+fl[12]-1.0*(fl[10]+fl[9]+fl[6]+fl[5])+fl[4]+fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[8] = 0.25*(fl[31]-1.0*(fl[30]+fl[29]+fl[27])+fl[26]+fl[25]+fl[22]+fl[21]-1.0*(fl[19]+fl[18]+fl[16]+fl[14])+fl[11]+fl[8]+fl[7]-1.0*fl[3]); 
    limQuad[8] = fl[0]/cflL; 
  } else {
    f0Quad[8] = -0.25*(fr[28]-1.0*(fr[24]+fr[23]+fr[20])+fr[17]+fr[15]+fr[13]+fr[12]-1.0*(fr[10]+fr[9]+fr[6]+fr[5])+fr[4]+fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[8] = -0.25*(fr[31]-1.0*(fr[30]+fr[29]+fr[27])+fr[26]+fr[25]+fr[22]+fr[21]-1.0*(fr[19]+fr[18]+fr[16]+fr[14])+fr[11]+fr[8]+fr[7]-1.0*fr[3]); 
    limQuad[8] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[9] = 0.25*(fl[28]+fl[24]-1.0*(fl[23]+fl[20])+fl[17]-1.0*(fl[15]+fl[13])+fl[12]+fl[10]-1.0*(fl[9]+fl[6])+fl[5]-1.0*(fl[4]+fl[2])+fl[1]+fl[0]); 
    f1Quad[9] = -0.25*(fl[31]+fl[30]-1.0*(fl[29]+fl[27])+fl[26]-1.0*(fl[25]+fl[22])+fl[21]+fl[19]-1.0*(fl[18]+fl[16])+fl[14]-1.0*(fl[11]+fl[8])+fl[7]+fl[3]); 
    limQuad[9] = fl[0]/cflL; 
  } else {
    f0Quad[9] = 0.25*(fr[28]+fr[24]-1.0*(fr[23]+fr[20])+fr[17]-1.0*(fr[15]+fr[13])+fr[12]+fr[10]-1.0*(fr[9]+fr[6])+fr[5]-1.0*(fr[4]+fr[2])+fr[1]+fr[0]); 
    f1Quad[9] = 0.25*(fr[31]+fr[30]-1.0*(fr[29]+fr[27])+fr[26]-1.0*(fr[25]+fr[22])+fr[21]+fr[19]-1.0*(fr[18]+fr[16])+fr[14]-1.0*(fr[11]+fr[8])+fr[7]+fr[3]); 
    limQuad[9] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[10] = 0.25*(fl[28]-1.0*fl[24]+fl[23]-1.0*fl[20]+fl[17]-1.0*fl[15]+fl[13]-1.0*(fl[12]+fl[10])+fl[9]-1.0*fl[6]+fl[5]-1.0*fl[4]+fl[2]-1.0*fl[1]+fl[0]); 
    f1Quad[10] = -0.25*(fl[31]-1.0*fl[30]+fl[29]-1.0*fl[27]+fl[26]-1.0*fl[25]+fl[22]-1.0*(fl[21]+fl[19])+fl[18]-1.0*fl[16]+fl[14]-1.0*fl[11]+fl[8]-1.0*fl[7]+fl[3]); 
    limQuad[10] = fl[0]/cflL; 
  } else {
    f0Quad[10] = 0.25*(fr[28]-1.0*fr[24]+fr[23]-1.0*fr[20]+fr[17]-1.0*fr[15]+fr[13]-1.0*(fr[12]+fr[10])+fr[9]-1.0*fr[6]+fr[5]-1.0*fr[4]+fr[2]-1.0*fr[1]+fr[0]); 
    f1Quad[10] = 0.25*(fr[31]-1.0*fr[30]+fr[29]-1.0*fr[27]+fr[26]-1.0*fr[25]+fr[22]-1.0*(fr[21]+fr[19])+fr[18]-1.0*fr[16]+fr[14]-1.0*fr[11]+fr[8]-1.0*fr[7]+fr[3]); 
    limQuad[10] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[11] = -0.25*(fl[28]+fl[24]+fl[23]-1.0*fl[20]+fl[17]+fl[15]-1.0*(fl[13]+fl[12])+fl[10]+fl[9]-1.0*(fl[6]+fl[5])+fl[4]-1.0*(fl[2]+fl[1]+fl[0])); 
    f1Quad[11] = 0.25*(fl[31]+fl[30]+fl[29]-1.0*fl[27]+fl[26]+fl[25]-1.0*(fl[22]+fl[21])+fl[19]+fl[18]-1.0*(fl[16]+fl[14])+fl[11]-1.0*(fl[8]+fl[7]+fl[3])); 
    limQuad[11] = fl[0]/cflL; 
  } else {
    f0Quad[11] = -0.25*(fr[28]+fr[24]+fr[23]-1.0*fr[20]+fr[17]+fr[15]-1.0*(fr[13]+fr[12])+fr[10]+fr[9]-1.0*(fr[6]+fr[5])+fr[4]-1.0*(fr[2]+fr[1]+fr[0])); 
    f1Quad[11] = -0.25*(fr[31]+fr[30]+fr[29]-1.0*fr[27]+fr[26]+fr[25]-1.0*(fr[22]+fr[21])+fr[19]+fr[18]-1.0*(fr[16]+fr[14])+fr[11]-1.0*(fr[8]+fr[7]+fr[3])); 
    limQuad[11] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[12] = 0.25*(fl[28]-1.0*(fl[24]+fl[23])+fl[20]+fl[17]+fl[15]-1.0*(fl[13]+fl[12]+fl[10]+fl[9])+fl[6]+fl[5]+fl[4]-1.0*(fl[2]+fl[1])+fl[0]); 
    f1Quad[12] = -0.25*(fl[31]-1.0*(fl[30]+fl[29])+fl[27]+fl[26]+fl[25]-1.0*(fl[22]+fl[21]+fl[19]+fl[18])+fl[16]+fl[14]+fl[11]-1.0*(fl[8]+fl[7])+fl[3]); 
    limQuad[12] = fl[0]/cflL; 
  } else {
    f0Quad[12] = 0.25*(fr[28]-1.0*(fr[24]+fr[23])+fr[20]+fr[17]+fr[15]-1.0*(fr[13]+fr[12]+fr[10]+fr[9])+fr[6]+fr[5]+fr[4]-1.0*(fr[2]+fr[1])+fr[0]); 
    f1Quad[12] = 0.25*(fr[31]-1.0*(fr[30]+fr[29])+fr[27]+fr[26]+fr[25]-1.0*(fr[22]+fr[21]+fr[19]+fr[18])+fr[16]+fr[14]+fr[11]-1.0*(fr[8]+fr[7])+fr[3]); 
    limQuad[12] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[13] = -0.25*(fl[28]+fl[24]-1.0*fl[23]+fl[20]+fl[17]-1.0*fl[15]+fl[13]-1.0*fl[12]+fl[10]-1.0*fl[9]+fl[6]-1.0*(fl[5]+fl[4])+fl[2]-1.0*(fl[1]+fl[0])); 
    f1Quad[13] = 0.25*(fl[31]+fl[30]-1.0*fl[29]+fl[27]+fl[26]-1.0*fl[25]+fl[22]-1.0*fl[21]+fl[19]-1.0*fl[18]+fl[16]-1.0*(fl[14]+fl[11])+fl[8]-1.0*(fl[7]+fl[3])); 
    limQuad[13] = fl[0]/cflL; 
  } else {
    f0Quad[13] = -0.25*(fr[28]+fr[24]-1.0*fr[23]+fr[20]+fr[17]-1.0*fr[15]+fr[13]-1.0*fr[12]+fr[10]-1.0*fr[9]+fr[6]-1.0*(fr[5]+fr[4])+fr[2]-1.0*(fr[1]+fr[0])); 
    f1Quad[13] = -0.25*(fr[31]+fr[30]-1.0*fr[29]+fr[27]+fr[26]-1.0*fr[25]+fr[22]-1.0*fr[21]+fr[19]-1.0*fr[18]+fr[16]-1.0*(fr[14]+fr[11])+fr[8]-1.0*(fr[7]+fr[3])); 
    limQuad[13] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[14] = -0.25*(fl[28]-1.0*fl[24]+fl[23]+fl[20]+fl[17]-1.0*(fl[15]+fl[13])+fl[12]-1.0*fl[10]+fl[9]+fl[6]-1.0*(fl[5]+fl[4]+fl[2])+fl[1]-1.0*fl[0]); 
    f1Quad[14] = 0.25*(fl[31]-1.0*fl[30]+fl[29]+fl[27]+fl[26]-1.0*(fl[25]+fl[22])+fl[21]-1.0*fl[19]+fl[18]+fl[16]-1.0*(fl[14]+fl[11]+fl[8])+fl[7]-1.0*fl[3]); 
    limQuad[14] = fl[0]/cflL; 
  } else {
    f0Quad[14] = -0.25*(fr[28]-1.0*fr[24]+fr[23]+fr[20]+fr[17]-1.0*(fr[15]+fr[13])+fr[12]-1.0*fr[10]+fr[9]+fr[6]-1.0*(fr[5]+fr[4]+fr[2])+fr[1]-1.0*fr[0]); 
    f1Quad[14] = -0.25*(fr[31]-1.0*fr[30]+fr[29]+fr[27]+fr[26]-1.0*(fr[25]+fr[22])+fr[21]-1.0*fr[19]+fr[18]+fr[16]-1.0*(fr[14]+fr[11]+fr[8])+fr[7]-1.0*fr[3]); 
    limQuad[14] = fr[0]/cflR; 
  }
  if(0.25*alpha[0] > 0) {
    f0Quad[15] = 0.25*(fl[28]+fl[24]+fl[23]+fl[20]+fl[17]+fl[15]+fl[13]+fl[12]+fl[10]+fl[9]+fl[6]+fl[5]+fl[4]+fl[2]+fl[1]+fl[0]); 
    f1Quad[15] = -0.25*(fl[31]+fl[30]+fl[29]+fl[27]+fl[26]+fl[25]+fl[22]+fl[21]+fl[19]+fl[18]+fl[16]+fl[14]+fl[11]+fl[8]+fl[7]+fl[3]); 
    limQuad[15] = fl[0]/cflL; 
  } else {
    f0Quad[15] = 0.25*(fr[28]+fr[24]+fr[23]+fr[20]+fr[17]+fr[15]+fr[13]+fr[12]+fr[10]+fr[9]+fr[6]+fr[5]+fr[4]+fr[2]+fr[1]+fr[0]); 
    f1Quad[15] = 0.25*(fr[31]+fr[30]+fr[29]+fr[27]+fr[26]+fr[25]+fr[22]+fr[21]+fr[19]+fr[18]+fr[16]+fr[14]+fr[11]+fr[8]+fr[7]+fr[3]); 
    limQuad[15] = fr[0]/cflR; 
  }
  double fhat[32]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[3] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[4] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[6] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[7] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[8] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[9] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[10] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[11] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[12] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[13] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[14] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[15] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[16] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[17] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[18] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[19] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[20] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[21] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[22] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[23] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[24] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[25] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[26] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[27] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[28] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[29] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[30] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[31] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[16];  // rCtrl=f1/f0 at each control node in dimensions other than z 
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[22]+fhat[21]+fhat[19]+fhat[18]+fhat[16]+3.0*fhat[3]))-3.0*(fhat[30]+fhat[29]+fhat[27]+fhat[26]+3.0*(fhat[14]+fhat[11]+fhat[8]+fhat[7]))))/(72.0*EPSILON+1.414213562373095*(fhat[28]-1.732050807568877*(fhat[24]+fhat[23]+fhat[20]+fhat[17]+3.0*(fhat[5]+fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[15]+fhat[13]+fhat[12]+fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0]))); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[22]))+fhat[21]-1.0*fhat[19]+fhat[18]+fhat[16]-3.0*fhat[3]))+3.0*(fhat[30]-1.0*(fhat[29]+fhat[27]+fhat[26]-3.0*(fhat[14]+fhat[11]+fhat[8]-1.0*fhat[7])))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[15]+fhat[13]-1.0*fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0]))-1.0*(fhat[28]+1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20]+fhat[17])+3.0*(fhat[5]+fhat[4]+fhat[2]-1.0*fhat[1]))))); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[22]-1.0*fhat[21]+fhat[19]-1.0*fhat[18]+fhat[16]-3.0*fhat[3]))+3.0*((-1.0*fhat[30])+fhat[29]-1.0*(fhat[27]+fhat[26]-3.0*(fhat[14]+fhat[11]-1.0*fhat[8]+fhat[7])))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[28])+1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]+fhat[17]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]-1.0*fhat[1]))+3.0*(fhat[15]-1.0*fhat[13]+fhat[12]-1.0*fhat[10]+fhat[9]-1.0*fhat[6]+3.0*fhat[0]))); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[22]+fhat[21]+fhat[19]+fhat[18]-1.0*(fhat[16]+3.0*fhat[3]))))+3.0*(fhat[30]+fhat[29]-1.0*(fhat[27]+fhat[26]-3.0*((-1.0*(fhat[14]+fhat[11]))+fhat[8]+fhat[7])))))/(72.0*EPSILON+1.414213562373095*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]-1.0*(fhat[20]+fhat[17])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]+fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[4] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[22]+fhat[21]-1.0*fhat[19])+fhat[18]-1.0*(fhat[16]+3.0*fhat[3])))+3.0*((-1.0*(fhat[30]+fhat[29]))+fhat[27]-1.0*fhat[26]+3.0*(fhat[14]-1.0*fhat[11]+fhat[8]+fhat[7]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[28])+1.732050807568877*(fhat[24]+fhat[23]-1.0*fhat[20]+fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[2]+fhat[1])))+3.0*((-1.0*fhat[15])+fhat[13]+fhat[12]-1.0*(fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[5] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[22]-1.0*(fhat[21]+fhat[19]-1.0*fhat[18]+fhat[16])+3.0*fhat[3]))+3.0*(fhat[30]-1.0*fhat[29]+fhat[27]-1.0*fhat[26]+3.0*((-1.0*fhat[14])+fhat[11]-1.0*fhat[8]+fhat[7]))))/(72.0*EPSILON+1.414213562373095*(fhat[28]+1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]-1.0*fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[13]-1.0*(fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6])+3.0*fhat[0]))); 
  rCtrl[6] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[22]))+fhat[21]+fhat[19]-1.0*(fhat[18]+fhat[16]-3.0*fhat[3])))+3.0*((-1.0*fhat[30])+fhat[29]+fhat[27]-1.0*fhat[26]+3.0*((-1.0*fhat[14])+fhat[11]+fhat[8]-1.0*fhat[7]))))/(72.0*EPSILON+1.414213562373095*(fhat[28]-1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20])+fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[7] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[22]+fhat[21]-1.0*(fhat[19]+fhat[18]+fhat[16]+3.0*fhat[3])))+3.0*(fhat[30]+fhat[29]+fhat[27]-1.0*fhat[26]+3.0*(fhat[14]-1.0*(fhat[11]+fhat[8]+fhat[7])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*(fhat[15]+fhat[13]+fhat[12]-1.0*fhat[10]))+fhat[9]+fhat[6]+3.0*fhat[0])-1.0*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]+fhat[20]-1.0*fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]+fhat[1])))))); 
  rCtrl[8] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[22]+fhat[21]-1.0*(fhat[19]+fhat[18]+fhat[16]+3.0*fhat[3])))+3.0*(3.0*((-1.0*fhat[14])+fhat[11]+fhat[8]+fhat[7])-1.0*(fhat[30]+fhat[29]+fhat[27]-1.0*fhat[26]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[28])+1.732050807568877*(fhat[24]+fhat[23]+fhat[20]-1.0*fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]+fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[13]+fhat[12]-1.0*fhat[10]))+fhat[9]+fhat[6]+3.0*fhat[0]))); 
  rCtrl[9] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[22]))+fhat[21]+fhat[19]-1.0*(fhat[18]+fhat[16]-3.0*fhat[3])))+3.0*(fhat[30]-1.0*(fhat[29]+fhat[27]-1.0*fhat[26])+3.0*(fhat[14]-1.0*(fhat[11]+fhat[8]-1.0*fhat[7])))))/(72.0*EPSILON+1.414213562373095*(fhat[28]+1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20])+fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[10] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[22]-1.0*(fhat[21]+fhat[19]-1.0*fhat[18]+fhat[16])+3.0*fhat[3]))+3.0*((-1.0*fhat[30])+fhat[29]-1.0*fhat[27]+fhat[26]+3.0*(fhat[14]-1.0*fhat[11]+fhat[8]-1.0*fhat[7]))))/(72.0*EPSILON+1.414213562373095*(fhat[28]-1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]-1.0*fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[13]-1.0*(fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6])+3.0*fhat[0]))); 
  rCtrl[11] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[22]+fhat[21]-1.0*fhat[19])+fhat[18]-1.0*(fhat[16]+3.0*fhat[3])))+3.0*(fhat[30]+fhat[29]-1.0*fhat[27]+fhat[26]+3.0*((-1.0*fhat[14])+fhat[11]-1.0*(fhat[8]+fhat[7])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*fhat[15])+fhat[13]+fhat[12]-1.0*(fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0])))-1.0*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]-1.0*fhat[20]+fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[2]+fhat[1])))))); 
  rCtrl[12] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[22]+fhat[21]+fhat[19]+fhat[18]-1.0*(fhat[16]+3.0*fhat[3]))))+3.0*((-1.0*(fhat[30]+fhat[29]))+fhat[27]+fhat[26]+3.0*(fhat[14]+fhat[11]-1.0*(fhat[8]+fhat[7])))))/(72.0*EPSILON+1.414213562373095*(fhat[28]-1.732050807568877*(fhat[24]+fhat[23]-1.0*(fhat[20]+fhat[17])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]+fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[13] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[22]-1.0*fhat[21]+fhat[19]-1.0*fhat[18]+fhat[16]-3.0*fhat[3]))+3.0*(fhat[30]-1.0*fhat[29]+fhat[27]+fhat[26]+3.0*((-1.0*(fhat[14]+fhat[11]))+fhat[8]-1.0*fhat[7]))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[15]-1.0*fhat[13]+fhat[12]-1.0*fhat[10]+fhat[9]-1.0*fhat[6]+3.0*fhat[0])-1.0*(fhat[28]+1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]+fhat[17]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]-1.0*fhat[1]))))); 
  rCtrl[14] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[22]))+fhat[21]-1.0*fhat[19]+fhat[18]+fhat[16]-3.0*fhat[3]))+3.0*((-1.0*fhat[30])+fhat[29]+fhat[27]+fhat[26]-3.0*(fhat[14]+fhat[11]+fhat[8]-1.0*fhat[7]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[28])+1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20]+fhat[17])+3.0*(fhat[5]+fhat[4]+fhat[2]-1.0*fhat[1]))+3.0*(fhat[15]+fhat[13]-1.0*fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[15] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[22]+fhat[21]+fhat[19]+fhat[18]+fhat[16]+3.0*fhat[3]))+3.0*(fhat[30]+fhat[29]+fhat[27]+fhat[26]+3.0*(fhat[14]+fhat[11]+fhat[8]+fhat[7]))))/(72.0*EPSILON+1.414213562373095*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]+fhat[20]+fhat[17]+3.0*(fhat[5]+fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[15]+fhat[13]+fhat[12]+fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0]))); 
  double fhatCtrl[16];  // fhatCtrl = anti-limited fhat evaluated at each control node on z surface 
  fhatCtrl[0] = -0.01964185503295965*(1.732050807568877*(fhat[24]+fhat[23]+fhat[20]+fhat[17]+3.0*(fhat[5]+fhat[4]+fhat[2]+fhat[1]))-1.0*(fhat[28]+3.0*(fhat[15]+fhat[13]+fhat[12]+fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[0],-1.0,EPSILON); 
  fhatCtrl[1] = -0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20]+fhat[17])+3.0*(fhat[5]+fhat[4]+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]-1.0*fhat[10]+fhat[9]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[1],-1.0,EPSILON); 
  fhatCtrl[2] = 0.01964185503295965*(1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]+fhat[17]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]-1.0*fhat[1]))-1.0*(fhat[28]+3.0*((-1.0*fhat[15])+fhat[13]-1.0*fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[2],-1.0,EPSILON); 
  fhatCtrl[3] = 0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]-1.0*(fhat[20]+fhat[17])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]+fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[3],-1.0,EPSILON); 
  fhatCtrl[4] = 0.01964185503295965*(1.732050807568877*(fhat[24]+fhat[23]-1.0*fhat[20]+fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[2]+fhat[1])))-1.0*(fhat[28]+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]-1.0*fhat[10])+fhat[9]-1.0*(fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[4],-1.0,EPSILON); 
  fhatCtrl[5] = 0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]-1.0*fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[13]-1.0*(fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6])+3.0*fhat[0]))*limTheta(rCtrl[5],-1.0,EPSILON); 
  fhatCtrl[6] = -0.01964185503295965*(1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20])+fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]-1.0*fhat[1])))-1.0*(fhat[28]+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0]))))*limTheta(rCtrl[6],-1.0,EPSILON); 
  fhatCtrl[7] = -0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]+fhat[20]-1.0*fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]+fhat[1])))+3.0*(fhat[15]+fhat[13]+fhat[12]-1.0*(fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[7],-1.0,EPSILON); 
  fhatCtrl[8] = 0.01964185503295965*(1.732050807568877*(fhat[24]+fhat[23]+fhat[20]-1.0*fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]+fhat[1])))-1.0*(fhat[28]+3.0*(fhat[15]+fhat[13]+fhat[12]-1.0*(fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[8],-1.0,EPSILON); 
  fhatCtrl[9] = 0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20])+fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[9],-1.0,EPSILON); 
  fhatCtrl[10] = -0.01964185503295965*(1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]-1.0*fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[2]+fhat[1]))-1.0*(fhat[28]+3.0*((-1.0*fhat[15])+fhat[13]-1.0*(fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6])+3.0*fhat[0])))*limTheta(rCtrl[10],-1.0,EPSILON); 
  fhatCtrl[11] = -0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]-1.0*fhat[20]+fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[2]+fhat[1])))+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]-1.0*fhat[10])+fhat[9]-1.0*(fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[11],-1.0,EPSILON); 
  fhatCtrl[12] = -0.01964185503295965*(1.732050807568877*(fhat[24]+fhat[23]-1.0*(fhat[20]+fhat[17])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]+fhat[1]))-1.0*(fhat[28]+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]+fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0])))))*limTheta(rCtrl[12],-1.0,EPSILON); 
  fhatCtrl[13] = -0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]+fhat[17]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*fhat[15])+fhat[13]-1.0*fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[13],-1.0,EPSILON); 
  fhatCtrl[14] = 0.01964185503295965*(1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20]+fhat[17])+3.0*(fhat[5]+fhat[4]+fhat[2]-1.0*fhat[1]))-1.0*(fhat[28]+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]-1.0*fhat[10]+fhat[9]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[14],-1.0,EPSILON); 
  fhatCtrl[15] = 0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]+fhat[20]+fhat[17]+3.0*(fhat[5]+fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[15]+fhat[13]+fhat[12]+fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0]))*limTheta(rCtrl[15],-1.0,EPSILON); 
  double fhatAL[16];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.25*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.4330127018922193*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[5] = 0.75*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[6] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[8] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[9] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[10] = 0.75*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[11] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[12] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[13] = 1.299038105676658*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[14] = 1.299038105676658*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[15] = 2.25*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[16]; 
  fhatALQuad[0] = fmin(0.25*((-0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7])))-1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[0]); 
  fhatALQuad[1] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*((-0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+0.5773502691896258*(fhatAL[13]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[1]); 
  fhatALQuad[2] = fmin(0.25*(0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[2]); 
  fhatALQuad[3] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[3]); 
  fhatALQuad[4] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*(0.5773502691896258*fhatAL[3]-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[4]+fhatAL[2]+fhatAL[1]))-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[4]); 
  fhatALQuad[5] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]))+0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[5]); 
  fhatALQuad[6] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2])-0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[6]); 
  fhatALQuad[7] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[7]); 
  fhatALQuad[8] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]+fhatAL[1]))+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[8]); 
  fhatALQuad[9] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]))+0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[9]); 
  fhatALQuad[10] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[10]); 
  fhatALQuad[11] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[4]+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[11]); 
  fhatALQuad[12] = fmin(0.25*((-0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5])))+1.732050807568877*(0.5773502691896258*(fhatAL[4]+fhatAL[3])-0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[12]); 
  fhatALQuad[13] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*(0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1])-0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[13]); 
  fhatALQuad[14] = fmin(0.25*(0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[1])+0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[14]); 
  fhatALQuad[15] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[15]); 
  fhatAL[0] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[5] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[6] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[8] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[9] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[10] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[11] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[12] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[13] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[14] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[15] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.1767766952966368*alpha[0]*fhatAL[0]*dfac_z; 
  incr[1] = 0.1767766952966368*alpha[0]*fhatAL[1]*dfac_z; 
  incr[2] = 0.1767766952966368*alpha[0]*fhatAL[2]*dfac_z; 
  incr[3] = -0.3061862178478971*alpha[0]*fhatAL[0]*dfac_z; 
  incr[4] = 0.1767766952966368*alpha[0]*fhatAL[3]*dfac_z; 
  incr[5] = 0.1767766952966368*alpha[0]*fhatAL[4]*dfac_z; 
  incr[6] = 0.1767766952966368*alpha[0]*fhatAL[5]*dfac_z; 
  incr[7] = -0.3061862178478971*alpha[0]*fhatAL[1]*dfac_z; 
  incr[8] = -0.3061862178478971*alpha[0]*fhatAL[2]*dfac_z; 
  incr[9] = 0.1767766952966368*alpha[0]*fhatAL[6]*dfac_z; 
  incr[10] = 0.1767766952966368*alpha[0]*fhatAL[7]*dfac_z; 
  incr[11] = -0.3061862178478971*alpha[0]*fhatAL[3]*dfac_z; 
  incr[12] = 0.1767766952966368*alpha[0]*fhatAL[8]*dfac_z; 
  incr[13] = 0.1767766952966368*alpha[0]*fhatAL[9]*dfac_z; 
  incr[14] = -0.3061862178478971*alpha[0]*fhatAL[4]*dfac_z; 
  incr[15] = 0.1767766952966368*alpha[0]*fhatAL[10]*dfac_z; 
  incr[16] = -0.3061862178478971*alpha[0]*fhatAL[5]*dfac_z; 
  incr[17] = 0.1767766952966368*alpha[0]*fhatAL[11]*dfac_z; 
  incr[18] = -0.3061862178478971*alpha[0]*fhatAL[6]*dfac_z; 
  incr[19] = -0.3061862178478971*alpha[0]*fhatAL[7]*dfac_z; 
  incr[20] = 0.1767766952966368*alpha[0]*fhatAL[12]*dfac_z; 
  incr[21] = -0.3061862178478971*alpha[0]*fhatAL[8]*dfac_z; 
  incr[22] = -0.3061862178478971*alpha[0]*fhatAL[9]*dfac_z; 
  incr[23] = 0.1767766952966368*alpha[0]*fhatAL[13]*dfac_z; 
  incr[24] = 0.1767766952966368*alpha[0]*fhatAL[14]*dfac_z; 
  incr[25] = -0.3061862178478971*alpha[0]*fhatAL[10]*dfac_z; 
  incr[26] = -0.3061862178478971*alpha[0]*fhatAL[11]*dfac_z; 
  incr[27] = -0.3061862178478971*alpha[0]*fhatAL[12]*dfac_z; 
  incr[28] = 0.1767766952966368*alpha[0]*fhatAL[15]*dfac_z; 
  incr[29] = -0.3061862178478971*alpha[0]*fhatAL[13]*dfac_z; 
  incr[30] = -0.3061862178478971*alpha[0]*fhatAL[14]*dfac_z; 
  incr[31] = -0.3061862178478971*alpha[0]*fhatAL[15]*dfac_z; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 
  outr[28] += incr[28]; 
  outr[29] += incr[29]; 
  outr[30] += incr[30]; 
  outr[31] += incr[31]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += incr[19]; 
  outl[20] += -1.0*incr[20]; 
  outl[21] += incr[21]; 
  outl[22] += incr[22]; 
  outl[23] += -1.0*incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += incr[25]; 
  outl[26] += incr[26]; 
  outl[27] += incr[27]; 
  outl[28] += -1.0*incr[28]; 
  outl[29] += incr[29]; 
  outl[30] += incr[30]; 
  outl[31] += incr[31]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity3x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_z = 2.0/dxv[2]; 
  double dfac_v = 2.0/dxv[3]; 
  double dfac_m = 2.0/dxv[4]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wz = w[2]; 
  double wv = w[3]; 
  double wm = w[4]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[32]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.1082531754730548*(dfac_v*((BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*wv+Gradpar[0]*Phi[3]*dfac_z*q_)-1.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_))/(dfac_v*m_); 

  double alpha[16]; 
  alpha[0] = -(0.8660254037844386*(dfac_v*((BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_*wv+Gradpar[0]*Phi[3]*dfac_z*q_)-1.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_))/(dfac_v*m_); 
  alpha[1] = -(0.8660254037844386*(dfac_v*(BdriftY[0]*Phi[4]*dfac_y*m_*wv+Gradpar[0]*Phi[5]*dfac_z*q_)-1.0*BdriftY[0]*Phi[4]*dfac_y*m_))/(dfac_v*m_); 
  alpha[2] = -(0.8660254037844386*(dfac_v*(BdriftX[0]*Phi[4]*dfac_x*m_*wv+Gradpar[0]*Phi[6]*dfac_z*q_)-1.0*BdriftX[0]*Phi[4]*dfac_x*m_))/(dfac_v*m_); 
  alpha[3] = -(0.8660254037844386*(BdriftY[0]*Phi[6]*dfac_y+BdriftX[0]*Phi[5]*dfac_x)*(dfac_v*wv-1.0))/dfac_v; 
  alpha[5] = -(0.8660254037844386*Gradpar[0]*Phi[7]*dfac_z*q_)/m_; 
  alpha[6] = -(0.8660254037844386*BdriftY[0]*Phi[7]*dfac_y*(dfac_v*wv-1.0))/dfac_v; 
  alpha[7] = -(0.8660254037844386*BdriftX[0]*Phi[7]*dfac_x*(dfac_v*wv-1.0))/dfac_v; 
  double f0Quad[16]; 
  double f1Quad[16]; 
  double limQuad[16]; 
  // determine upwinding at each surface quadrature node 
  if(0.25*(alpha[7]+alpha[6]+alpha[5])-0.25*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[0] = 0.25*(fl[27]-1.0*(fl[22]+fl[21]+fl[20]+fl[16])+fl[14]+fl[13]+fl[12]+fl[8]+fl[7]+fl[6]-1.0*(fl[5]+fl[3]+fl[2]+fl[1])+fl[0]); 
    f1Quad[0] = -0.25*(fl[31]-1.0*(fl[30]+fl[29]+fl[28]+fl[26])+fl[25]+fl[24]+fl[23]+fl[19]+fl[18]+fl[17]-1.0*(fl[15]+fl[11]+fl[10]+fl[9])+fl[4]); 
    limQuad[0] = fl[0]/cflL; 
  } else {
    f0Quad[0] = 0.25*(fr[27]-1.0*(fr[22]+fr[21]+fr[20]+fr[16])+fr[14]+fr[13]+fr[12]+fr[8]+fr[7]+fr[6]-1.0*(fr[5]+fr[3]+fr[2]+fr[1])+fr[0]); 
    f1Quad[0] = 0.25*(fr[31]-1.0*(fr[30]+fr[29]+fr[28]+fr[26])+fr[25]+fr[24]+fr[23]+fr[19]+fr[18]+fr[17]-1.0*(fr[15]+fr[11]+fr[10]+fr[9])+fr[4]); 
    limQuad[0] = fr[0]/cflR; 
  }
  if(0.25*alpha[7]-0.25*(alpha[6]+alpha[5]+alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[1] = -0.25*(fl[27]+fl[22]-1.0*(fl[21]+fl[20]+fl[16]+fl[14]+fl[13])+fl[12]-1.0*fl[8]+fl[7]+fl[6]+fl[5]+fl[3]+fl[2]-1.0*(fl[1]+fl[0])); 
    f1Quad[1] = 0.25*(fl[31]+fl[30]-1.0*(fl[29]+fl[28]+fl[26]+fl[25]+fl[24])+fl[23]-1.0*fl[19]+fl[18]+fl[17]+fl[15]+fl[11]+fl[10]-1.0*(fl[9]+fl[4])); 
    limQuad[1] = fl[0]/cflL; 
  } else {
    f0Quad[1] = -0.25*(fr[27]+fr[22]-1.0*(fr[21]+fr[20]+fr[16]+fr[14]+fr[13])+fr[12]-1.0*fr[8]+fr[7]+fr[6]+fr[5]+fr[3]+fr[2]-1.0*(fr[1]+fr[0])); 
    f1Quad[1] = -0.25*(fr[31]+fr[30]-1.0*(fr[29]+fr[28]+fr[26]+fr[25]+fr[24])+fr[23]-1.0*fr[19]+fr[18]+fr[17]+fr[15]+fr[11]+fr[10]-1.0*(fr[9]+fr[4])); 
    limQuad[1] = fr[0]/cflR; 
  }
  if((-0.25*alpha[7])+0.25*alpha[6]-0.25*(alpha[5]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[2] = -0.25*(fl[27]-1.0*fl[22]+fl[21]-1.0*(fl[20]+fl[16]+fl[14])+fl[13]-1.0*fl[12]+fl[8]-1.0*fl[7]+fl[6]+fl[5]+fl[3]-1.0*fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[2] = 0.25*(fl[31]-1.0*fl[30]+fl[29]-1.0*(fl[28]+fl[26]+fl[25])+fl[24]-1.0*fl[23]+fl[19]-1.0*fl[18]+fl[17]+fl[15]+fl[11]-1.0*fl[10]+fl[9]-1.0*fl[4]); 
    limQuad[2] = fl[0]/cflL; 
  } else {
    f0Quad[2] = -0.25*(fr[27]-1.0*fr[22]+fr[21]-1.0*(fr[20]+fr[16]+fr[14])+fr[13]-1.0*fr[12]+fr[8]-1.0*fr[7]+fr[6]+fr[5]+fr[3]-1.0*fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[2] = -0.25*(fr[31]-1.0*fr[30]+fr[29]-1.0*(fr[28]+fr[26]+fr[25])+fr[24]-1.0*fr[23]+fr[19]-1.0*fr[18]+fr[17]+fr[15]+fr[11]-1.0*fr[10]+fr[9]-1.0*fr[4]); 
    limQuad[2] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[7]+alpha[6]))+0.25*alpha[5]-0.25*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[3] = 0.25*(fl[27]+fl[22]+fl[21]-1.0*(fl[20]+fl[16])+fl[14]-1.0*(fl[13]+fl[12]+fl[8]+fl[7])+fl[6]-1.0*(fl[5]+fl[3])+fl[2]+fl[1]+fl[0]); 
    f1Quad[3] = -0.25*(fl[31]+fl[30]+fl[29]-1.0*(fl[28]+fl[26])+fl[25]-1.0*(fl[24]+fl[23]+fl[19]+fl[18])+fl[17]-1.0*(fl[15]+fl[11])+fl[10]+fl[9]+fl[4]); 
    limQuad[3] = fl[0]/cflL; 
  } else {
    f0Quad[3] = 0.25*(fr[27]+fr[22]+fr[21]-1.0*(fr[20]+fr[16])+fr[14]-1.0*(fr[13]+fr[12]+fr[8]+fr[7])+fr[6]-1.0*(fr[5]+fr[3])+fr[2]+fr[1]+fr[0]); 
    f1Quad[3] = 0.25*(fr[31]+fr[30]+fr[29]-1.0*(fr[28]+fr[26])+fr[25]-1.0*(fr[24]+fr[23]+fr[19]+fr[18])+fr[17]-1.0*(fr[15]+fr[11])+fr[10]+fr[9]+fr[4]); 
    limQuad[3] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[7]+alpha[6]))+0.25*(alpha[5]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[4] = -0.25*(fl[27]-1.0*(fl[22]+fl[21])+fl[20]-1.0*fl[16]+fl[14]-1.0*(fl[13]+fl[12])+fl[8]+fl[7]-1.0*fl[6]+fl[5]-1.0*fl[3]+fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[4] = 0.25*(fl[31]-1.0*(fl[30]+fl[29])+fl[28]-1.0*fl[26]+fl[25]-1.0*(fl[24]+fl[23])+fl[19]+fl[18]-1.0*fl[17]+fl[15]-1.0*fl[11]+fl[10]+fl[9]-1.0*fl[4]); 
    limQuad[4] = fl[0]/cflL; 
  } else {
    f0Quad[4] = -0.25*(fr[27]-1.0*(fr[22]+fr[21])+fr[20]-1.0*fr[16]+fr[14]-1.0*(fr[13]+fr[12])+fr[8]+fr[7]-1.0*fr[6]+fr[5]-1.0*fr[3]+fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[4] = -0.25*(fr[31]-1.0*(fr[30]+fr[29])+fr[28]-1.0*fr[26]+fr[25]-1.0*(fr[24]+fr[23])+fr[19]+fr[18]-1.0*fr[17]+fr[15]-1.0*fr[11]+fr[10]+fr[9]-1.0*fr[4]); 
    limQuad[4] = fr[0]/cflR; 
  }
  if((-0.25*alpha[7])+0.25*alpha[6]-0.25*alpha[5]+0.25*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[5] = 0.25*(fl[27]+fl[22]-1.0*fl[21]+fl[20]-1.0*(fl[16]+fl[14])+fl[13]-1.0*(fl[12]+fl[8])+fl[7]-1.0*(fl[6]+fl[5])+fl[3]-1.0*fl[2]+fl[1]+fl[0]); 
    f1Quad[5] = -0.25*(fl[31]+fl[30]-1.0*fl[29]+fl[28]-1.0*(fl[26]+fl[25])+fl[24]-1.0*(fl[23]+fl[19])+fl[18]-1.0*(fl[17]+fl[15])+fl[11]-1.0*fl[10]+fl[9]+fl[4]); 
    limQuad[5] = fl[0]/cflL; 
  } else {
    f0Quad[5] = 0.25*(fr[27]+fr[22]-1.0*fr[21]+fr[20]-1.0*(fr[16]+fr[14])+fr[13]-1.0*(fr[12]+fr[8])+fr[7]-1.0*(fr[6]+fr[5])+fr[3]-1.0*fr[2]+fr[1]+fr[0]); 
    f1Quad[5] = 0.25*(fr[31]+fr[30]-1.0*fr[29]+fr[28]-1.0*(fr[26]+fr[25])+fr[24]-1.0*(fr[23]+fr[19])+fr[18]-1.0*(fr[17]+fr[15])+fr[11]-1.0*fr[10]+fr[9]+fr[4]); 
    limQuad[5] = fr[0]/cflR; 
  }
  if(0.25*alpha[7]-0.25*(alpha[6]+alpha[5])+0.25*(alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[6] = 0.25*(fl[27]-1.0*fl[22]+fl[21]+fl[20]-1.0*(fl[16]+fl[14]+fl[13])+fl[12]+fl[8]-1.0*(fl[7]+fl[6]+fl[5])+fl[3]+fl[2]-1.0*fl[1]+fl[0]); 
    f1Quad[6] = -0.25*(fl[31]-1.0*fl[30]+fl[29]+fl[28]-1.0*(fl[26]+fl[25]+fl[24])+fl[23]+fl[19]-1.0*(fl[18]+fl[17]+fl[15])+fl[11]+fl[10]-1.0*fl[9]+fl[4]); 
    limQuad[6] = fl[0]/cflL; 
  } else {
    f0Quad[6] = 0.25*(fr[27]-1.0*fr[22]+fr[21]+fr[20]-1.0*(fr[16]+fr[14]+fr[13])+fr[12]+fr[8]-1.0*(fr[7]+fr[6]+fr[5])+fr[3]+fr[2]-1.0*fr[1]+fr[0]); 
    f1Quad[6] = 0.25*(fr[31]-1.0*fr[30]+fr[29]+fr[28]-1.0*(fr[26]+fr[25]+fr[24])+fr[23]+fr[19]-1.0*(fr[18]+fr[17]+fr[15])+fr[11]+fr[10]-1.0*fr[9]+fr[4]); 
    limQuad[6] = fr[0]/cflR; 
  }
  if(0.25*(alpha[7]+alpha[6]+alpha[5]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[7] = -0.25*(fl[27]+fl[22]+fl[21]+fl[20]-1.0*fl[16]+fl[14]+fl[13]+fl[12]-1.0*(fl[8]+fl[7]+fl[6])+fl[5]-1.0*(fl[3]+fl[2]+fl[1]+fl[0])); 
    f1Quad[7] = 0.25*(fl[31]+fl[30]+fl[29]+fl[28]-1.0*fl[26]+fl[25]+fl[24]+fl[23]-1.0*(fl[19]+fl[18]+fl[17])+fl[15]-1.0*(fl[11]+fl[10]+fl[9]+fl[4])); 
    limQuad[7] = fl[0]/cflL; 
  } else {
    f0Quad[7] = -0.25*(fr[27]+fr[22]+fr[21]+fr[20]-1.0*fr[16]+fr[14]+fr[13]+fr[12]-1.0*(fr[8]+fr[7]+fr[6])+fr[5]-1.0*(fr[3]+fr[2]+fr[1]+fr[0])); 
    f1Quad[7] = -0.25*(fr[31]+fr[30]+fr[29]+fr[28]-1.0*fr[26]+fr[25]+fr[24]+fr[23]-1.0*(fr[19]+fr[18]+fr[17])+fr[15]-1.0*(fr[11]+fr[10]+fr[9]+fr[4])); 
    limQuad[7] = fr[0]/cflR; 
  }
  if(0.25*(alpha[7]+alpha[6]+alpha[5])-0.25*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[8] = -0.25*(fl[27]-1.0*(fl[22]+fl[21]+fl[20])+fl[16]+fl[14]+fl[13]+fl[12]-1.0*(fl[8]+fl[7]+fl[6]+fl[5])+fl[3]+fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[8] = 0.25*(fl[31]-1.0*(fl[30]+fl[29]+fl[28])+fl[26]+fl[25]+fl[24]+fl[23]-1.0*(fl[19]+fl[18]+fl[17]+fl[15])+fl[11]+fl[10]+fl[9]-1.0*fl[4]); 
    limQuad[8] = fl[0]/cflL; 
  } else {
    f0Quad[8] = -0.25*(fr[27]-1.0*(fr[22]+fr[21]+fr[20])+fr[16]+fr[14]+fr[13]+fr[12]-1.0*(fr[8]+fr[7]+fr[6]+fr[5])+fr[3]+fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[8] = -0.25*(fr[31]-1.0*(fr[30]+fr[29]+fr[28])+fr[26]+fr[25]+fr[24]+fr[23]-1.0*(fr[19]+fr[18]+fr[17]+fr[15])+fr[11]+fr[10]+fr[9]-1.0*fr[4]); 
    limQuad[8] = fr[0]/cflR; 
  }
  if(0.25*alpha[7]-0.25*(alpha[6]+alpha[5]+alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[9] = 0.25*(fl[27]+fl[22]-1.0*(fl[21]+fl[20])+fl[16]-1.0*(fl[14]+fl[13])+fl[12]+fl[8]-1.0*(fl[7]+fl[6])+fl[5]-1.0*(fl[3]+fl[2])+fl[1]+fl[0]); 
    f1Quad[9] = -0.25*(fl[31]+fl[30]-1.0*(fl[29]+fl[28])+fl[26]-1.0*(fl[25]+fl[24])+fl[23]+fl[19]-1.0*(fl[18]+fl[17])+fl[15]-1.0*(fl[11]+fl[10])+fl[9]+fl[4]); 
    limQuad[9] = fl[0]/cflL; 
  } else {
    f0Quad[9] = 0.25*(fr[27]+fr[22]-1.0*(fr[21]+fr[20])+fr[16]-1.0*(fr[14]+fr[13])+fr[12]+fr[8]-1.0*(fr[7]+fr[6])+fr[5]-1.0*(fr[3]+fr[2])+fr[1]+fr[0]); 
    f1Quad[9] = 0.25*(fr[31]+fr[30]-1.0*(fr[29]+fr[28])+fr[26]-1.0*(fr[25]+fr[24])+fr[23]+fr[19]-1.0*(fr[18]+fr[17])+fr[15]-1.0*(fr[11]+fr[10])+fr[9]+fr[4]); 
    limQuad[9] = fr[0]/cflR; 
  }
  if((-0.25*alpha[7])+0.25*alpha[6]-0.25*(alpha[5]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[10] = 0.25*(fl[27]-1.0*fl[22]+fl[21]-1.0*fl[20]+fl[16]-1.0*fl[14]+fl[13]-1.0*(fl[12]+fl[8])+fl[7]-1.0*fl[6]+fl[5]-1.0*fl[3]+fl[2]-1.0*fl[1]+fl[0]); 
    f1Quad[10] = -0.25*(fl[31]-1.0*fl[30]+fl[29]-1.0*fl[28]+fl[26]-1.0*fl[25]+fl[24]-1.0*(fl[23]+fl[19])+fl[18]-1.0*fl[17]+fl[15]-1.0*fl[11]+fl[10]-1.0*fl[9]+fl[4]); 
    limQuad[10] = fl[0]/cflL; 
  } else {
    f0Quad[10] = 0.25*(fr[27]-1.0*fr[22]+fr[21]-1.0*fr[20]+fr[16]-1.0*fr[14]+fr[13]-1.0*(fr[12]+fr[8])+fr[7]-1.0*fr[6]+fr[5]-1.0*fr[3]+fr[2]-1.0*fr[1]+fr[0]); 
    f1Quad[10] = 0.25*(fr[31]-1.0*fr[30]+fr[29]-1.0*fr[28]+fr[26]-1.0*fr[25]+fr[24]-1.0*(fr[23]+fr[19])+fr[18]-1.0*fr[17]+fr[15]-1.0*fr[11]+fr[10]-1.0*fr[9]+fr[4]); 
    limQuad[10] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[7]+alpha[6]))+0.25*alpha[5]-0.25*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[11] = -0.25*(fl[27]+fl[22]+fl[21]-1.0*fl[20]+fl[16]+fl[14]-1.0*(fl[13]+fl[12])+fl[8]+fl[7]-1.0*(fl[6]+fl[5])+fl[3]-1.0*(fl[2]+fl[1]+fl[0])); 
    f1Quad[11] = 0.25*(fl[31]+fl[30]+fl[29]-1.0*fl[28]+fl[26]+fl[25]-1.0*(fl[24]+fl[23])+fl[19]+fl[18]-1.0*(fl[17]+fl[15])+fl[11]-1.0*(fl[10]+fl[9]+fl[4])); 
    limQuad[11] = fl[0]/cflL; 
  } else {
    f0Quad[11] = -0.25*(fr[27]+fr[22]+fr[21]-1.0*fr[20]+fr[16]+fr[14]-1.0*(fr[13]+fr[12])+fr[8]+fr[7]-1.0*(fr[6]+fr[5])+fr[3]-1.0*(fr[2]+fr[1]+fr[0])); 
    f1Quad[11] = -0.25*(fr[31]+fr[30]+fr[29]-1.0*fr[28]+fr[26]+fr[25]-1.0*(fr[24]+fr[23])+fr[19]+fr[18]-1.0*(fr[17]+fr[15])+fr[11]-1.0*(fr[10]+fr[9]+fr[4])); 
    limQuad[11] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[7]+alpha[6]))+0.25*(alpha[5]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[12] = 0.25*(fl[27]-1.0*(fl[22]+fl[21])+fl[20]+fl[16]+fl[14]-1.0*(fl[13]+fl[12]+fl[8]+fl[7])+fl[6]+fl[5]+fl[3]-1.0*(fl[2]+fl[1])+fl[0]); 
    f1Quad[12] = -0.25*(fl[31]-1.0*(fl[30]+fl[29])+fl[28]+fl[26]+fl[25]-1.0*(fl[24]+fl[23]+fl[19]+fl[18])+fl[17]+fl[15]+fl[11]-1.0*(fl[10]+fl[9])+fl[4]); 
    limQuad[12] = fl[0]/cflL; 
  } else {
    f0Quad[12] = 0.25*(fr[27]-1.0*(fr[22]+fr[21])+fr[20]+fr[16]+fr[14]-1.0*(fr[13]+fr[12]+fr[8]+fr[7])+fr[6]+fr[5]+fr[3]-1.0*(fr[2]+fr[1])+fr[0]); 
    f1Quad[12] = 0.25*(fr[31]-1.0*(fr[30]+fr[29])+fr[28]+fr[26]+fr[25]-1.0*(fr[24]+fr[23]+fr[19]+fr[18])+fr[17]+fr[15]+fr[11]-1.0*(fr[10]+fr[9])+fr[4]); 
    limQuad[12] = fr[0]/cflR; 
  }
  if((-0.25*alpha[7])+0.25*alpha[6]-0.25*alpha[5]+0.25*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[13] = -0.25*(fl[27]+fl[22]-1.0*fl[21]+fl[20]+fl[16]-1.0*fl[14]+fl[13]-1.0*fl[12]+fl[8]-1.0*fl[7]+fl[6]-1.0*(fl[5]+fl[3])+fl[2]-1.0*(fl[1]+fl[0])); 
    f1Quad[13] = 0.25*(fl[31]+fl[30]-1.0*fl[29]+fl[28]+fl[26]-1.0*fl[25]+fl[24]-1.0*fl[23]+fl[19]-1.0*fl[18]+fl[17]-1.0*(fl[15]+fl[11])+fl[10]-1.0*(fl[9]+fl[4])); 
    limQuad[13] = fl[0]/cflL; 
  } else {
    f0Quad[13] = -0.25*(fr[27]+fr[22]-1.0*fr[21]+fr[20]+fr[16]-1.0*fr[14]+fr[13]-1.0*fr[12]+fr[8]-1.0*fr[7]+fr[6]-1.0*(fr[5]+fr[3])+fr[2]-1.0*(fr[1]+fr[0])); 
    f1Quad[13] = -0.25*(fr[31]+fr[30]-1.0*fr[29]+fr[28]+fr[26]-1.0*fr[25]+fr[24]-1.0*fr[23]+fr[19]-1.0*fr[18]+fr[17]-1.0*(fr[15]+fr[11])+fr[10]-1.0*(fr[9]+fr[4])); 
    limQuad[13] = fr[0]/cflR; 
  }
  if(0.25*alpha[7]-0.25*(alpha[6]+alpha[5])+0.25*(alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[14] = -0.25*(fl[27]-1.0*fl[22]+fl[21]+fl[20]+fl[16]-1.0*(fl[14]+fl[13])+fl[12]-1.0*fl[8]+fl[7]+fl[6]-1.0*(fl[5]+fl[3]+fl[2])+fl[1]-1.0*fl[0]); 
    f1Quad[14] = 0.25*(fl[31]-1.0*fl[30]+fl[29]+fl[28]+fl[26]-1.0*(fl[25]+fl[24])+fl[23]-1.0*fl[19]+fl[18]+fl[17]-1.0*(fl[15]+fl[11]+fl[10])+fl[9]-1.0*fl[4]); 
    limQuad[14] = fl[0]/cflL; 
  } else {
    f0Quad[14] = -0.25*(fr[27]-1.0*fr[22]+fr[21]+fr[20]+fr[16]-1.0*(fr[14]+fr[13])+fr[12]-1.0*fr[8]+fr[7]+fr[6]-1.0*(fr[5]+fr[3]+fr[2])+fr[1]-1.0*fr[0]); 
    f1Quad[14] = -0.25*(fr[31]-1.0*fr[30]+fr[29]+fr[28]+fr[26]-1.0*(fr[25]+fr[24])+fr[23]-1.0*fr[19]+fr[18]+fr[17]-1.0*(fr[15]+fr[11]+fr[10])+fr[9]-1.0*fr[4]); 
    limQuad[14] = fr[0]/cflR; 
  }
  if(0.25*(alpha[7]+alpha[6]+alpha[5]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[15] = 0.25*(fl[27]+fl[22]+fl[21]+fl[20]+fl[16]+fl[14]+fl[13]+fl[12]+fl[8]+fl[7]+fl[6]+fl[5]+fl[3]+fl[2]+fl[1]+fl[0]); 
    f1Quad[15] = -0.25*(fl[31]+fl[30]+fl[29]+fl[28]+fl[26]+fl[25]+fl[24]+fl[23]+fl[19]+fl[18]+fl[17]+fl[15]+fl[11]+fl[10]+fl[9]+fl[4]); 
    limQuad[15] = fl[0]/cflL; 
  } else {
    f0Quad[15] = 0.25*(fr[27]+fr[22]+fr[21]+fr[20]+fr[16]+fr[14]+fr[13]+fr[12]+fr[8]+fr[7]+fr[6]+fr[5]+fr[3]+fr[2]+fr[1]+fr[0]); 
    f1Quad[15] = 0.25*(fr[31]+fr[30]+fr[29]+fr[28]+fr[26]+fr[25]+fr[24]+fr[23]+fr[19]+fr[18]+fr[17]+fr[15]+fr[11]+fr[10]+fr[9]+fr[4]); 
    limQuad[15] = fr[0]/cflR; 
  }
  double fhat[32]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[3] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[5] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[6] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[7] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[8] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[9] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[10] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[11] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[12] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[13] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[14] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[15] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[16] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[17] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[18] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[19] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[20] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[21] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[22] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[23] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[24] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[25] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[26] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[27] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[28] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[29] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[30] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[31] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[16];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[24]+fhat[23]+fhat[19]+fhat[18]+fhat[17]+3.0*fhat[4]))-3.0*(fhat[30]+fhat[29]+fhat[28]+fhat[26]+3.0*(fhat[15]+fhat[11]+fhat[10]+fhat[9]))))/(72.0*EPSILON+1.414213562373095*(fhat[27]-1.732050807568877*(fhat[22]+fhat[21]+fhat[20]+fhat[16]+3.0*(fhat[5]+fhat[3]+fhat[2]+fhat[1]))+3.0*(fhat[14]+fhat[13]+fhat[12]+fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0]))); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[23]-1.0*fhat[19]+fhat[18]+fhat[17]-3.0*fhat[4]))+3.0*(fhat[30]-1.0*(fhat[29]+fhat[28]+fhat[26]-3.0*(fhat[15]+fhat[11]+fhat[10]-1.0*fhat[9])))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[14]+fhat[13]-1.0*fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0]))-1.0*(fhat[27]+1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20]+fhat[16])+3.0*(fhat[5]+fhat[3]+fhat[2]-1.0*fhat[1]))))); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*fhat[23]+fhat[19]-1.0*fhat[18]+fhat[17]-3.0*fhat[4]))+3.0*((-1.0*fhat[30])+fhat[29]-1.0*(fhat[28]+fhat[26]-3.0*(fhat[15]+fhat[11]-1.0*fhat[10]+fhat[9])))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[27])+1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]+fhat[16]+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]-1.0*fhat[1]))+3.0*(fhat[14]-1.0*fhat[13]+fhat[12]-1.0*fhat[8]+fhat[7]-1.0*fhat[6]+3.0*fhat[0]))); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[23]+fhat[19]+fhat[18]-1.0*(fhat[17]+3.0*fhat[4]))))+3.0*(fhat[30]+fhat[29]-1.0*(fhat[28]+fhat[26]-3.0*((-1.0*(fhat[15]+fhat[11]))+fhat[10]+fhat[9])))))/(72.0*EPSILON+1.414213562373095*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]-1.0*(fhat[20]+fhat[16])+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]+fhat[1]))+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]+fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[4] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[23]-1.0*fhat[19])+fhat[18]-1.0*(fhat[17]+3.0*fhat[4])))+3.0*((-1.0*(fhat[30]+fhat[29]))+fhat[28]-1.0*fhat[26]+3.0*(fhat[15]-1.0*fhat[11]+fhat[10]+fhat[9]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[27])+1.732050807568877*(fhat[22]+fhat[21]-1.0*fhat[20]+fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*(fhat[2]+fhat[1])))+3.0*((-1.0*fhat[14])+fhat[13]+fhat[12]-1.0*(fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[5] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*(fhat[23]+fhat[19]-1.0*fhat[18]+fhat[17])+3.0*fhat[4]))+3.0*(fhat[30]-1.0*fhat[29]+fhat[28]-1.0*fhat[26]+3.0*((-1.0*fhat[15])+fhat[11]-1.0*fhat[10]+fhat[9]))))/(72.0*EPSILON+1.414213562373095*(fhat[27]+1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]-1.0*fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[14])+fhat[13]-1.0*(fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6])+3.0*fhat[0]))); 
  rCtrl[6] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[23]+fhat[19]-1.0*(fhat[18]+fhat[17]-3.0*fhat[4])))+3.0*((-1.0*fhat[30])+fhat[29]+fhat[28]-1.0*fhat[26]+3.0*((-1.0*fhat[15])+fhat[11]+fhat[10]-1.0*fhat[9]))))/(72.0*EPSILON+1.414213562373095*(fhat[27]-1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20])+fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[7] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[24]+fhat[23]-1.0*(fhat[19]+fhat[18]+fhat[17]+3.0*fhat[4])))+3.0*(fhat[30]+fhat[29]+fhat[28]-1.0*fhat[26]+3.0*(fhat[15]-1.0*(fhat[11]+fhat[10]+fhat[9])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*(fhat[14]+fhat[13]+fhat[12]-1.0*fhat[8]))+fhat[7]+fhat[6]+3.0*fhat[0])-1.0*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]+fhat[20]-1.0*fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]+fhat[1])))))); 
  rCtrl[8] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[24]+fhat[23]-1.0*(fhat[19]+fhat[18]+fhat[17]+3.0*fhat[4])))+3.0*(3.0*((-1.0*fhat[15])+fhat[11]+fhat[10]+fhat[9])-1.0*(fhat[30]+fhat[29]+fhat[28]-1.0*fhat[26]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[27])+1.732050807568877*(fhat[22]+fhat[21]+fhat[20]-1.0*fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]+fhat[1])))+3.0*((-1.0*(fhat[14]+fhat[13]+fhat[12]-1.0*fhat[8]))+fhat[7]+fhat[6]+3.0*fhat[0]))); 
  rCtrl[9] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[23]+fhat[19]-1.0*(fhat[18]+fhat[17]-3.0*fhat[4])))+3.0*(fhat[30]-1.0*(fhat[29]+fhat[28]-1.0*fhat[26])+3.0*(fhat[15]-1.0*(fhat[11]+fhat[10]-1.0*fhat[9])))))/(72.0*EPSILON+1.414213562373095*(fhat[27]+1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20])+fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[10] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*(fhat[23]+fhat[19]-1.0*fhat[18]+fhat[17])+3.0*fhat[4]))+3.0*((-1.0*fhat[30])+fhat[29]-1.0*fhat[28]+fhat[26]+3.0*(fhat[15]-1.0*fhat[11]+fhat[10]-1.0*fhat[9]))))/(72.0*EPSILON+1.414213562373095*(fhat[27]-1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]-1.0*fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[14])+fhat[13]-1.0*(fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6])+3.0*fhat[0]))); 
  rCtrl[11] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[23]-1.0*fhat[19])+fhat[18]-1.0*(fhat[17]+3.0*fhat[4])))+3.0*(fhat[30]+fhat[29]-1.0*fhat[28]+fhat[26]+3.0*((-1.0*fhat[15])+fhat[11]-1.0*(fhat[10]+fhat[9])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*fhat[14])+fhat[13]+fhat[12]-1.0*(fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0])))-1.0*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]-1.0*fhat[20]+fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*(fhat[2]+fhat[1])))))); 
  rCtrl[12] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[23]+fhat[19]+fhat[18]-1.0*(fhat[17]+3.0*fhat[4]))))+3.0*((-1.0*(fhat[30]+fhat[29]))+fhat[28]+fhat[26]+3.0*(fhat[15]+fhat[11]-1.0*(fhat[10]+fhat[9])))))/(72.0*EPSILON+1.414213562373095*(fhat[27]-1.732050807568877*(fhat[22]+fhat[21]-1.0*(fhat[20]+fhat[16])+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]+fhat[1]))+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]+fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[13] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*fhat[23]+fhat[19]-1.0*fhat[18]+fhat[17]-3.0*fhat[4]))+3.0*(fhat[30]-1.0*fhat[29]+fhat[28]+fhat[26]+3.0*((-1.0*(fhat[15]+fhat[11]))+fhat[10]-1.0*fhat[9]))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[14]-1.0*fhat[13]+fhat[12]-1.0*fhat[8]+fhat[7]-1.0*fhat[6]+3.0*fhat[0])-1.0*(fhat[27]+1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]+fhat[16]+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]-1.0*fhat[1]))))); 
  rCtrl[14] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[23]-1.0*fhat[19]+fhat[18]+fhat[17]-3.0*fhat[4]))+3.0*((-1.0*fhat[30])+fhat[29]+fhat[28]+fhat[26]-3.0*(fhat[15]+fhat[11]+fhat[10]-1.0*fhat[9]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[27])+1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20]+fhat[16])+3.0*(fhat[5]+fhat[3]+fhat[2]-1.0*fhat[1]))+3.0*(fhat[14]+fhat[13]-1.0*fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[15] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[24]+fhat[23]+fhat[19]+fhat[18]+fhat[17]+3.0*fhat[4]))+3.0*(fhat[30]+fhat[29]+fhat[28]+fhat[26]+3.0*(fhat[15]+fhat[11]+fhat[10]+fhat[9]))))/(72.0*EPSILON+1.414213562373095*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]+fhat[20]+fhat[16]+3.0*(fhat[5]+fhat[3]+fhat[2]+fhat[1]))+3.0*(fhat[14]+fhat[13]+fhat[12]+fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0]))); 
  double fhatCtrl[16];  // fhatCtrl = anti-limited fhat evaluated at each control node on vx surface 
  fhatCtrl[0] = -0.01964185503295965*(1.732050807568877*(fhat[22]+fhat[21]+fhat[20]+fhat[16]+3.0*(fhat[5]+fhat[3]+fhat[2]+fhat[1]))-1.0*(fhat[27]+3.0*(fhat[14]+fhat[13]+fhat[12]+fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[0],-1.0,EPSILON); 
  fhatCtrl[1] = -0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20]+fhat[16])+3.0*(fhat[5]+fhat[3]+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]-1.0*fhat[8]+fhat[7]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[1],-1.0,EPSILON); 
  fhatCtrl[2] = 0.01964185503295965*(1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]+fhat[16]+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]-1.0*fhat[1]))-1.0*(fhat[27]+3.0*((-1.0*fhat[14])+fhat[13]-1.0*fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[2],-1.0,EPSILON); 
  fhatCtrl[3] = 0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]-1.0*(fhat[20]+fhat[16])+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]+fhat[1]))+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]+fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[3],-1.0,EPSILON); 
  fhatCtrl[4] = 0.01964185503295965*(1.732050807568877*(fhat[22]+fhat[21]-1.0*fhat[20]+fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*(fhat[2]+fhat[1])))-1.0*(fhat[27]+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]-1.0*fhat[8])+fhat[7]-1.0*(fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[4],-1.0,EPSILON); 
  fhatCtrl[5] = 0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]-1.0*fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[14])+fhat[13]-1.0*(fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6])+3.0*fhat[0]))*limTheta(rCtrl[5],-1.0,EPSILON); 
  fhatCtrl[6] = -0.01964185503295965*(1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20])+fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]-1.0*fhat[1])))-1.0*(fhat[27]+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0]))))*limTheta(rCtrl[6],-1.0,EPSILON); 
  fhatCtrl[7] = -0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]+fhat[20]-1.0*fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]+fhat[1])))+3.0*(fhat[14]+fhat[13]+fhat[12]-1.0*(fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[7],-1.0,EPSILON); 
  fhatCtrl[8] = 0.01964185503295965*(1.732050807568877*(fhat[22]+fhat[21]+fhat[20]-1.0*fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]+fhat[1])))-1.0*(fhat[27]+3.0*(fhat[14]+fhat[13]+fhat[12]-1.0*(fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[8],-1.0,EPSILON); 
  fhatCtrl[9] = 0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20])+fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[9],-1.0,EPSILON); 
  fhatCtrl[10] = -0.01964185503295965*(1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]-1.0*fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*fhat[2]+fhat[1]))-1.0*(fhat[27]+3.0*((-1.0*fhat[14])+fhat[13]-1.0*(fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6])+3.0*fhat[0])))*limTheta(rCtrl[10],-1.0,EPSILON); 
  fhatCtrl[11] = -0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]-1.0*fhat[20]+fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*(fhat[2]+fhat[1])))+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]-1.0*fhat[8])+fhat[7]-1.0*(fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[11],-1.0,EPSILON); 
  fhatCtrl[12] = -0.01964185503295965*(1.732050807568877*(fhat[22]+fhat[21]-1.0*(fhat[20]+fhat[16])+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]+fhat[1]))-1.0*(fhat[27]+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]+fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0])))))*limTheta(rCtrl[12],-1.0,EPSILON); 
  fhatCtrl[13] = -0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]+fhat[16]+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*fhat[14])+fhat[13]-1.0*fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[13],-1.0,EPSILON); 
  fhatCtrl[14] = 0.01964185503295965*(1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20]+fhat[16])+3.0*(fhat[5]+fhat[3]+fhat[2]-1.0*fhat[1]))-1.0*(fhat[27]+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]-1.0*fhat[8]+fhat[7]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[14],-1.0,EPSILON); 
  fhatCtrl[15] = 0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]+fhat[20]+fhat[16]+3.0*(fhat[5]+fhat[3]+fhat[2]+fhat[1]))+3.0*(fhat[14]+fhat[13]+fhat[12]+fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0]))*limTheta(rCtrl[15],-1.0,EPSILON); 
  double fhatAL[16];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.25*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.4330127018922193*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[5] = 0.75*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[6] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[8] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[9] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[10] = 0.75*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[11] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[12] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[13] = 1.299038105676658*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[14] = 1.299038105676658*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[15] = 2.25*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[16]; 
  fhatALQuad[0] = fmin(0.25*((-0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[7]))-1.732050807568877*(fhatAL[10]+fhatAL[6])))-1.0*(1.0*(fhatAL[14]+fhatAL[11])+fhatAL[13]+fhatAL[12]+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[0]); 
  fhatALQuad[1] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[10]+1.732050807568877*fhatAL[6]))+1.732050807568877*(0.5773502691896258*fhatAL[1]-0.5773502691896258*((-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[14]))-1.0*(fhatAL[13]+fhatAL[12])+fhatAL[4]+fhatAL[3]+fhatAL[2]))-0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[9])-1.0*fhatAL[8]+fhatAL[0]), limQuad[1]); 
  fhatALQuad[2] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[7])-1.732050807568877*(fhatAL[10]+fhatAL[6])))+1.732050807568877*((-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[11]))+fhatAL[13]+fhatAL[3]))+0.5773502691896258*(fhatAL[12]+fhatAL[2])-0.5773502691896258*(fhatAL[4]+fhatAL[1]))-1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[2]); 
  fhatALQuad[3] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15])-1.732050807568877*fhatAL[10]+1.732050807568877*fhatAL[6]))+1.732050807568877*((-0.5773502691896258*(0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[14])-1.0*fhatAL[13]+fhatAL[3]))+0.5773502691896258*(fhatAL[2]-1.0*fhatAL[12])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[1])+0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[9])-1.0*fhatAL[8]+fhatAL[0]), limQuad[3]); 
  fhatALQuad[4] = fmin(0.25*(0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[7]))-1.732050807568877*(fhatAL[10]+fhatAL[6]))+1.732050807568877*(0.5773502691896258*(1.0*(fhatAL[14]+fhatAL[11])+fhatAL[13]+fhatAL[3])-0.5773502691896258*(fhatAL[12]+fhatAL[4]+fhatAL[2]+fhatAL[1]))+1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[4]); 
  fhatALQuad[5] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[10]+1.732050807568877*fhatAL[6])+1.732050807568877*(0.5773502691896258*((-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[14]))-1.0*fhatAL[13]+fhatAL[3])-0.5773502691896258*((-1.0*fhatAL[12])+fhatAL[4]+fhatAL[2])+0.5773502691896258*fhatAL[1])-0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[9])-1.0*fhatAL[8]+fhatAL[0]), limQuad[5]); 
  fhatALQuad[6] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[7])-1.732050807568877*(fhatAL[10]+fhatAL[6]))+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[11]))+fhatAL[13]+fhatAL[12]+fhatAL[3]+fhatAL[2])-0.5773502691896258*(fhatAL[4]+fhatAL[1]))-1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[6]); 
  fhatALQuad[7] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15])-1.732050807568877*fhatAL[10]+1.732050807568877*fhatAL[6])+1.732050807568877*(0.5773502691896258*(0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[14])-1.0*(fhatAL[13]+fhatAL[12])+fhatAL[3]+fhatAL[2])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[1])+0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[9])-1.0*fhatAL[8]+fhatAL[0]), limQuad[7]); 
  fhatALQuad[8] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[10]-1.732050807568877*fhatAL[6]))+1.732050807568877*((-0.5773502691896258*((-0.5773502691896258*(1.732050807568877*fhatAL[14]-1.732050807568877*fhatAL[11]))-1.0*(fhatAL[13]+fhatAL[12])+fhatAL[3]+fhatAL[2]))+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[1])-0.5773502691896258*(1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[5])-1.0*fhatAL[8]+fhatAL[0]), limQuad[8]); 
  fhatALQuad[9] = fmin(0.25*((-0.5773502691896258*(1.732050807568877*(fhatAL[10]+fhatAL[6])-1.732050807568878*(fhatAL[15]+fhatAL[7])))+1.732050807568877*(0.5773502691896258*(fhatAL[4]+fhatAL[1])-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[11]))+fhatAL[13]+fhatAL[12]+fhatAL[3]+fhatAL[2]))-1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[9]); 
  fhatALQuad[10] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15])+1.732050807568877*fhatAL[10]-1.732050807568877*fhatAL[6]))+1.732050807568877*((-0.5773502691896258*(0.5773502691896258*(1.732050807568877*fhatAL[14]-1.732050807568877*fhatAL[11])-1.0*fhatAL[13]+fhatAL[3]))+0.5773502691896258*((-1.0*fhatAL[12])+fhatAL[4]+fhatAL[2])-0.5773502691896258*fhatAL[1])+0.5773502691896258*(1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[5])-1.0*fhatAL[8]+fhatAL[0]), limQuad[10]); 
  fhatALQuad[11] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[7])+1.732050807568877*(fhatAL[10]+fhatAL[6])))+1.732050807568877*(0.5773502691896258*(fhatAL[12]+fhatAL[4]+fhatAL[2]+fhatAL[1])-0.5773502691896258*(1.0*(fhatAL[14]+fhatAL[11])+fhatAL[13]+fhatAL[3]))+1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[11]); 
  fhatALQuad[12] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[10]-1.732050807568877*fhatAL[6])+1.732050807568877*(0.5773502691896258*((-0.5773502691896258*(1.732050807568877*fhatAL[14]-1.732050807568877*fhatAL[11]))-1.0*fhatAL[13]+fhatAL[3])-0.5773502691896258*(fhatAL[2]-1.0*fhatAL[12])+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[1])-0.5773502691896258*(1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[5])-1.0*fhatAL[8]+fhatAL[0]), limQuad[12]); 
  fhatALQuad[13] = fmin(0.25*(0.5773502691896258*(1.732050807568877*(fhatAL[10]+fhatAL[6])-1.732050807568878*(fhatAL[15]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[11]))+fhatAL[13]+fhatAL[3])-0.5773502691896258*(fhatAL[12]+fhatAL[2])+0.5773502691896258*(fhatAL[4]+fhatAL[1]))-1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[13]); 
  fhatALQuad[14] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15])+1.732050807568877*fhatAL[10]-1.732050807568877*fhatAL[6])+1.732050807568877*(0.5773502691896258*(0.5773502691896258*(1.732050807568877*fhatAL[14]-1.732050807568877*fhatAL[11])-1.0*(fhatAL[13]+fhatAL[12])+fhatAL[4]+fhatAL[3]+fhatAL[2])-0.5773502691896258*fhatAL[1])+0.5773502691896258*(1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[5])-1.0*fhatAL[8]+fhatAL[0]), limQuad[14]); 
  fhatALQuad[15] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[7])+1.732050807568877*(fhatAL[10]+fhatAL[6]))+1.0*(1.0*(fhatAL[14]+fhatAL[11])+fhatAL[13]+fhatAL[12]+fhatAL[9]+fhatAL[5]+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[8]+fhatAL[0]), limQuad[15]); 
  fhatAL[0] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[5] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[6] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[8] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[9] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[10] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[11] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[12] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[13] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[14] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[15] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.1767766952966368*(alpha[7]*fhatAL[7]+alpha[6]*fhatAL[6]+alpha[5]*fhatAL[5]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(alpha[7]*fhatAL[11]+alpha[3]*fhatAL[6]+fhatAL[3]*alpha[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(alpha[6]*fhatAL[11]+alpha[3]*fhatAL[7]+fhatAL[3]*alpha[7]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[3] = 0.1767766952966368*(alpha[5]*fhatAL[11]+alpha[2]*fhatAL[7]+fhatAL[2]*alpha[7]+alpha[1]*fhatAL[6]+fhatAL[1]*alpha[6]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_v; 
  incr[4] = -0.3061862178478971*(alpha[7]*fhatAL[7]+alpha[6]*fhatAL[6]+alpha[5]*fhatAL[5]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[5] = 0.1767766952966368*(alpha[7]*fhatAL[14]+alpha[6]*fhatAL[13]+alpha[5]*fhatAL[12]+alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[1]*fhatAL[8]+alpha[0]*fhatAL[4])*dfac_v; 
  incr[6] = 0.1767766952966368*(alpha[3]*fhatAL[11]+alpha[6]*fhatAL[7]+fhatAL[6]*alpha[7]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 
  incr[7] = 0.1767766952966368*(alpha[2]*fhatAL[11]+alpha[5]*fhatAL[7]+fhatAL[5]*alpha[7]+alpha[0]*fhatAL[6]+fhatAL[0]*alpha[6]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_v; 
  incr[8] = 0.1767766952966368*(alpha[1]*fhatAL[11]+alpha[0]*fhatAL[7]+fhatAL[0]*alpha[7]+alpha[5]*fhatAL[6]+fhatAL[5]*alpha[6]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_v; 
  incr[9] = -0.3061862178478971*(alpha[7]*fhatAL[11]+alpha[3]*fhatAL[6]+fhatAL[3]*alpha[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[10] = -0.3061862178478971*(alpha[6]*fhatAL[11]+alpha[3]*fhatAL[7]+fhatAL[3]*alpha[7]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[11] = -0.3061862178478971*(alpha[5]*fhatAL[11]+alpha[2]*fhatAL[7]+fhatAL[2]*alpha[7]+alpha[1]*fhatAL[6]+fhatAL[1]*alpha[6]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_v; 
  incr[12] = 0.1767766952966368*(alpha[7]*fhatAL[15]+alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[6]*fhatAL[10]+alpha[5]*fhatAL[9]+alpha[0]*fhatAL[8]+alpha[1]*fhatAL[4])*dfac_v; 
  incr[13] = 0.1767766952966368*(alpha[6]*fhatAL[15]+alpha[3]*fhatAL[14]+alpha[1]*fhatAL[12]+alpha[7]*fhatAL[10]+alpha[0]*fhatAL[9]+alpha[5]*fhatAL[8]+alpha[2]*fhatAL[4])*dfac_v; 
  incr[14] = 0.1767766952966368*(alpha[5]*fhatAL[15]+alpha[2]*fhatAL[14]+alpha[1]*fhatAL[13]+alpha[0]*fhatAL[10]+alpha[7]*fhatAL[9]+alpha[6]*fhatAL[8]+alpha[3]*fhatAL[4])*dfac_v; 
  incr[15] = -0.3061862178478971*(alpha[7]*fhatAL[14]+alpha[6]*fhatAL[13]+alpha[5]*fhatAL[12]+alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[1]*fhatAL[8]+alpha[0]*fhatAL[4])*dfac_v; 
  incr[16] = 0.1767766952966368*(alpha[0]*fhatAL[11]+alpha[1]*fhatAL[7]+fhatAL[1]*alpha[7]+alpha[2]*fhatAL[6]+fhatAL[2]*alpha[6]+alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5])*dfac_v; 
  incr[17] = -0.3061862178478971*(alpha[3]*fhatAL[11]+alpha[6]*fhatAL[7]+fhatAL[6]*alpha[7]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 
  incr[18] = -0.3061862178478971*(alpha[2]*fhatAL[11]+alpha[5]*fhatAL[7]+fhatAL[5]*alpha[7]+alpha[0]*fhatAL[6]+fhatAL[0]*alpha[6]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_v; 
  incr[19] = -0.3061862178478971*(alpha[1]*fhatAL[11]+alpha[0]*fhatAL[7]+fhatAL[0]*alpha[7]+alpha[5]*fhatAL[6]+fhatAL[5]*alpha[6]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_v; 
  incr[20] = 0.1767766952966368*(alpha[3]*fhatAL[15]+alpha[6]*fhatAL[14]+alpha[7]*fhatAL[13]+alpha[0]*fhatAL[12]+alpha[1]*fhatAL[9]+alpha[2]*fhatAL[8]+fhatAL[4]*alpha[5])*dfac_v; 
  incr[21] = 0.1767766952966368*(alpha[2]*fhatAL[15]+alpha[5]*fhatAL[14]+alpha[0]*fhatAL[13]+alpha[7]*fhatAL[12]+alpha[1]*fhatAL[10]+alpha[3]*fhatAL[8]+fhatAL[4]*alpha[6])*dfac_v; 
  incr[22] = 0.1767766952966368*(alpha[1]*fhatAL[15]+alpha[0]*fhatAL[14]+alpha[5]*fhatAL[13]+alpha[6]*fhatAL[12]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9]+fhatAL[4]*alpha[7])*dfac_v; 
  incr[23] = -0.3061862178478971*(alpha[7]*fhatAL[15]+alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[6]*fhatAL[10]+alpha[5]*fhatAL[9]+alpha[0]*fhatAL[8]+alpha[1]*fhatAL[4])*dfac_v; 
  incr[24] = -0.3061862178478971*(alpha[6]*fhatAL[15]+alpha[3]*fhatAL[14]+alpha[1]*fhatAL[12]+alpha[7]*fhatAL[10]+alpha[0]*fhatAL[9]+alpha[5]*fhatAL[8]+alpha[2]*fhatAL[4])*dfac_v; 
  incr[25] = -0.3061862178478971*(alpha[5]*fhatAL[15]+alpha[2]*fhatAL[14]+alpha[1]*fhatAL[13]+alpha[0]*fhatAL[10]+alpha[7]*fhatAL[9]+alpha[6]*fhatAL[8]+alpha[3]*fhatAL[4])*dfac_v; 
  incr[26] = -0.3061862178478971*(alpha[0]*fhatAL[11]+alpha[1]*fhatAL[7]+fhatAL[1]*alpha[7]+alpha[2]*fhatAL[6]+fhatAL[2]*alpha[6]+alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5])*dfac_v; 
  incr[27] = 0.1767766952966368*(alpha[0]*fhatAL[15]+alpha[1]*fhatAL[14]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12]+alpha[5]*fhatAL[10]+alpha[6]*fhatAL[9]+alpha[7]*fhatAL[8])*dfac_v; 
  incr[28] = -0.3061862178478971*(alpha[3]*fhatAL[15]+alpha[6]*fhatAL[14]+alpha[7]*fhatAL[13]+alpha[0]*fhatAL[12]+alpha[1]*fhatAL[9]+alpha[2]*fhatAL[8]+fhatAL[4]*alpha[5])*dfac_v; 
  incr[29] = -0.3061862178478971*(alpha[2]*fhatAL[15]+alpha[5]*fhatAL[14]+alpha[0]*fhatAL[13]+alpha[7]*fhatAL[12]+alpha[1]*fhatAL[10]+alpha[3]*fhatAL[8]+fhatAL[4]*alpha[6])*dfac_v; 
  incr[30] = -0.3061862178478971*(alpha[1]*fhatAL[15]+alpha[0]*fhatAL[14]+alpha[5]*fhatAL[13]+alpha[6]*fhatAL[12]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9]+fhatAL[4]*alpha[7])*dfac_v; 
  incr[31] = -0.3061862178478971*(alpha[0]*fhatAL[15]+alpha[1]*fhatAL[14]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12]+alpha[5]*fhatAL[10]+alpha[6]*fhatAL[9]+alpha[7]*fhatAL[8])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 
  outr[28] += incr[28]; 
  outr[29] += incr[29]; 
  outr[30] += incr[30]; 
  outr[31] += incr[31]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += incr[19]; 
  outl[20] += -1.0*incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += incr[23]; 
  outl[24] += incr[24]; 
  outl[25] += incr[25]; 
  outl[26] += incr[26]; 
  outl[27] += -1.0*incr[27]; 
  outl[28] += incr[28]; 
  outl[29] += incr[29]; 
  outl[30] += incr[30]; 
  outl[31] += incr[31]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity3x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_z = 2.0/dxv[2]; 
  double dfac_v = 2.0/dxv[3]; 
  double dfac_m = 2.0/dxv[4]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wz = w[2]; 
  double wv = w[3]; 
  double wm = w[4]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[32]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.0625*(1.732050807568877*(2.828427124746191*BdriftX[1]*m_*wv2+(3.0*BmagInv[1]*Phi[4]+BmagInv[0]*Phi[2])*dfac_y*q_)-2.828427124746191*BdriftX[0]*m_*wv2-3.0*(BmagInv[0]*Phi[4]+BmagInv[1]*Phi[2])*dfac_y*q_))/q_; 

  double alpha[16]; 
  alpha[0] = -(0.5*(1.732050807568877*(2.828427124746191*BdriftX[1]*m_*wv2+(3.0*BmagInv[1]*Phi[4]+BmagInv[0]*Phi[2])*dfac_y*q_)-2.828427124746191*BdriftX[0]*m_*wv2-3.0*(BmagInv[0]*Phi[4]+BmagInv[1]*Phi[2])*dfac_y*q_))/q_; 
  alpha[2] = -0.5*(1.732050807568877*(3.0*BmagInv[1]*Phi[7]+BmagInv[0]*Phi[6])-3.0*(BmagInv[0]*Phi[7]+BmagInv[1]*Phi[6]))*dfac_y; 
  alpha[3] = (0.3333333333333333*(2.449489742783178*BdriftX[0]-4.242640687119286*BdriftX[1])*m_*wv)/(dfac_v*q_); 
  double f0Quad[16]; 
  double f1Quad[16]; 
  double limQuad[16]; 
  // determine upwinding at each surface quadrature node 
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[0] = 0.25*(fl[30]-1.0*(fl[25]+fl[24]+fl[22]+fl[19])+fl[15]+fl[14]+fl[13]+fl[11]+fl[10]+fl[8]-1.0*(fl[5]+fl[4]+fl[3]+fl[2])+fl[0]); 
    f1Quad[0] = -0.25*(fl[31]-1.0*(fl[29]+fl[28]+fl[27]+fl[26])+fl[23]+fl[21]+fl[20]+fl[18]+fl[17]+fl[16]-1.0*(fl[12]+fl[9]+fl[7]+fl[6])+fl[1]); 
    limQuad[0] = fl[0]/cflL; 
  } else {
    f0Quad[0] = 0.25*(fr[30]-1.0*(fr[25]+fr[24]+fr[22]+fr[19])+fr[15]+fr[14]+fr[13]+fr[11]+fr[10]+fr[8]-1.0*(fr[5]+fr[4]+fr[3]+fr[2])+fr[0]); 
    f1Quad[0] = 0.25*(fr[31]-1.0*(fr[29]+fr[28]+fr[27]+fr[26])+fr[23]+fr[21]+fr[20]+fr[18]+fr[17]+fr[16]-1.0*(fr[12]+fr[9]+fr[7]+fr[6])+fr[1]); 
    limQuad[0] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[1] = -0.25*(fl[30]+fl[25]-1.0*(fl[24]+fl[22]+fl[19]+fl[15]+fl[14])+fl[13]-1.0*fl[11]+fl[10]+fl[8]+fl[5]+fl[4]+fl[3]-1.0*(fl[2]+fl[0])); 
    f1Quad[1] = 0.25*(fl[31]+fl[29]-1.0*(fl[28]+fl[27]+fl[26]+fl[23]+fl[21])+fl[20]-1.0*fl[18]+fl[17]+fl[16]+fl[12]+fl[9]+fl[7]-1.0*(fl[6]+fl[1])); 
    limQuad[1] = fl[0]/cflL; 
  } else {
    f0Quad[1] = -0.25*(fr[30]+fr[25]-1.0*(fr[24]+fr[22]+fr[19]+fr[15]+fr[14])+fr[13]-1.0*fr[11]+fr[10]+fr[8]+fr[5]+fr[4]+fr[3]-1.0*(fr[2]+fr[0])); 
    f1Quad[1] = -0.25*(fr[31]+fr[29]-1.0*(fr[28]+fr[27]+fr[26]+fr[23]+fr[21])+fr[20]-1.0*fr[18]+fr[17]+fr[16]+fr[12]+fr[9]+fr[7]-1.0*(fr[6]+fr[1])); 
    limQuad[1] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[2] = -0.25*(fl[30]-1.0*fl[25]+fl[24]-1.0*(fl[22]+fl[19]+fl[15])+fl[14]-1.0*fl[13]+fl[11]-1.0*fl[10]+fl[8]+fl[5]+fl[4]-1.0*fl[3]+fl[2]-1.0*fl[0]); 
    f1Quad[2] = 0.25*(fl[31]-1.0*fl[29]+fl[28]-1.0*(fl[27]+fl[26]+fl[23])+fl[21]-1.0*fl[20]+fl[18]-1.0*fl[17]+fl[16]+fl[12]+fl[9]-1.0*fl[7]+fl[6]-1.0*fl[1]); 
    limQuad[2] = fl[0]/cflL; 
  } else {
    f0Quad[2] = -0.25*(fr[30]-1.0*fr[25]+fr[24]-1.0*(fr[22]+fr[19]+fr[15])+fr[14]-1.0*fr[13]+fr[11]-1.0*fr[10]+fr[8]+fr[5]+fr[4]-1.0*fr[3]+fr[2]-1.0*fr[0]); 
    f1Quad[2] = -0.25*(fr[31]-1.0*fr[29]+fr[28]-1.0*(fr[27]+fr[26]+fr[23])+fr[21]-1.0*fr[20]+fr[18]-1.0*fr[17]+fr[16]+fr[12]+fr[9]-1.0*fr[7]+fr[6]-1.0*fr[1]); 
    limQuad[2] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[3] = 0.25*(fl[30]+fl[25]+fl[24]-1.0*(fl[22]+fl[19])+fl[15]-1.0*(fl[14]+fl[13]+fl[11]+fl[10])+fl[8]-1.0*(fl[5]+fl[4])+fl[3]+fl[2]+fl[0]); 
    f1Quad[3] = -0.25*(fl[31]+fl[29]+fl[28]-1.0*(fl[27]+fl[26])+fl[23]-1.0*(fl[21]+fl[20]+fl[18]+fl[17])+fl[16]-1.0*(fl[12]+fl[9])+fl[7]+fl[6]+fl[1]); 
    limQuad[3] = fl[0]/cflL; 
  } else {
    f0Quad[3] = 0.25*(fr[30]+fr[25]+fr[24]-1.0*(fr[22]+fr[19])+fr[15]-1.0*(fr[14]+fr[13]+fr[11]+fr[10])+fr[8]-1.0*(fr[5]+fr[4])+fr[3]+fr[2]+fr[0]); 
    f1Quad[3] = 0.25*(fr[31]+fr[29]+fr[28]-1.0*(fr[27]+fr[26])+fr[23]-1.0*(fr[21]+fr[20]+fr[18]+fr[17])+fr[16]-1.0*(fr[12]+fr[9])+fr[7]+fr[6]+fr[1]); 
    limQuad[3] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[4] = -0.25*(fl[30]-1.0*(fl[25]+fl[24])+fl[22]-1.0*fl[19]+fl[15]-1.0*(fl[14]+fl[13])+fl[11]+fl[10]-1.0*fl[8]+fl[5]-1.0*fl[4]+fl[3]+fl[2]-1.0*fl[0]); 
    f1Quad[4] = 0.25*(fl[31]-1.0*(fl[29]+fl[28])+fl[27]-1.0*fl[26]+fl[23]-1.0*(fl[21]+fl[20])+fl[18]+fl[17]-1.0*fl[16]+fl[12]-1.0*fl[9]+fl[7]+fl[6]-1.0*fl[1]); 
    limQuad[4] = fl[0]/cflL; 
  } else {
    f0Quad[4] = -0.25*(fr[30]-1.0*(fr[25]+fr[24])+fr[22]-1.0*fr[19]+fr[15]-1.0*(fr[14]+fr[13])+fr[11]+fr[10]-1.0*fr[8]+fr[5]-1.0*fr[4]+fr[3]+fr[2]-1.0*fr[0]); 
    f1Quad[4] = -0.25*(fr[31]-1.0*(fr[29]+fr[28])+fr[27]-1.0*fr[26]+fr[23]-1.0*(fr[21]+fr[20])+fr[18]+fr[17]-1.0*fr[16]+fr[12]-1.0*fr[9]+fr[7]+fr[6]-1.0*fr[1]); 
    limQuad[4] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[5] = 0.25*(fl[30]+fl[25]-1.0*fl[24]+fl[22]-1.0*(fl[19]+fl[15])+fl[14]-1.0*(fl[13]+fl[11])+fl[10]-1.0*(fl[8]+fl[5])+fl[4]-1.0*fl[3]+fl[2]+fl[0]); 
    f1Quad[5] = -0.25*(fl[31]+fl[29]-1.0*fl[28]+fl[27]-1.0*(fl[26]+fl[23])+fl[21]-1.0*(fl[20]+fl[18])+fl[17]-1.0*(fl[16]+fl[12])+fl[9]-1.0*fl[7]+fl[6]+fl[1]); 
    limQuad[5] = fl[0]/cflL; 
  } else {
    f0Quad[5] = 0.25*(fr[30]+fr[25]-1.0*fr[24]+fr[22]-1.0*(fr[19]+fr[15])+fr[14]-1.0*(fr[13]+fr[11])+fr[10]-1.0*(fr[8]+fr[5])+fr[4]-1.0*fr[3]+fr[2]+fr[0]); 
    f1Quad[5] = 0.25*(fr[31]+fr[29]-1.0*fr[28]+fr[27]-1.0*(fr[26]+fr[23])+fr[21]-1.0*(fr[20]+fr[18])+fr[17]-1.0*(fr[16]+fr[12])+fr[9]-1.0*fr[7]+fr[6]+fr[1]); 
    limQuad[5] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[6] = 0.25*(fl[30]-1.0*fl[25]+fl[24]+fl[22]-1.0*(fl[19]+fl[15]+fl[14])+fl[13]+fl[11]-1.0*(fl[10]+fl[8]+fl[5])+fl[4]+fl[3]-1.0*fl[2]+fl[0]); 
    f1Quad[6] = -0.25*(fl[31]-1.0*fl[29]+fl[28]+fl[27]-1.0*(fl[26]+fl[23]+fl[21])+fl[20]+fl[18]-1.0*(fl[17]+fl[16]+fl[12])+fl[9]+fl[7]-1.0*fl[6]+fl[1]); 
    limQuad[6] = fl[0]/cflL; 
  } else {
    f0Quad[6] = 0.25*(fr[30]-1.0*fr[25]+fr[24]+fr[22]-1.0*(fr[19]+fr[15]+fr[14])+fr[13]+fr[11]-1.0*(fr[10]+fr[8]+fr[5])+fr[4]+fr[3]-1.0*fr[2]+fr[0]); 
    f1Quad[6] = 0.25*(fr[31]-1.0*fr[29]+fr[28]+fr[27]-1.0*(fr[26]+fr[23]+fr[21])+fr[20]+fr[18]-1.0*(fr[17]+fr[16]+fr[12])+fr[9]+fr[7]-1.0*fr[6]+fr[1]); 
    limQuad[6] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[7] = -0.25*(fl[30]+fl[25]+fl[24]+fl[22]-1.0*fl[19]+fl[15]+fl[14]+fl[13]-1.0*(fl[11]+fl[10]+fl[8])+fl[5]-1.0*(fl[4]+fl[3]+fl[2]+fl[0])); 
    f1Quad[7] = 0.25*(fl[31]+fl[29]+fl[28]+fl[27]-1.0*fl[26]+fl[23]+fl[21]+fl[20]-1.0*(fl[18]+fl[17]+fl[16])+fl[12]-1.0*(fl[9]+fl[7]+fl[6]+fl[1])); 
    limQuad[7] = fl[0]/cflL; 
  } else {
    f0Quad[7] = -0.25*(fr[30]+fr[25]+fr[24]+fr[22]-1.0*fr[19]+fr[15]+fr[14]+fr[13]-1.0*(fr[11]+fr[10]+fr[8])+fr[5]-1.0*(fr[4]+fr[3]+fr[2]+fr[0])); 
    f1Quad[7] = -0.25*(fr[31]+fr[29]+fr[28]+fr[27]-1.0*fr[26]+fr[23]+fr[21]+fr[20]-1.0*(fr[18]+fr[17]+fr[16])+fr[12]-1.0*(fr[9]+fr[7]+fr[6]+fr[1])); 
    limQuad[7] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[8] = -0.25*(fl[30]-1.0*(fl[25]+fl[24]+fl[22])+fl[19]+fl[15]+fl[14]+fl[13]-1.0*(fl[11]+fl[10]+fl[8]+fl[5])+fl[4]+fl[3]+fl[2]-1.0*fl[0]); 
    f1Quad[8] = 0.25*(fl[31]-1.0*(fl[29]+fl[28]+fl[27])+fl[26]+fl[23]+fl[21]+fl[20]-1.0*(fl[18]+fl[17]+fl[16]+fl[12])+fl[9]+fl[7]+fl[6]-1.0*fl[1]); 
    limQuad[8] = fl[0]/cflL; 
  } else {
    f0Quad[8] = -0.25*(fr[30]-1.0*(fr[25]+fr[24]+fr[22])+fr[19]+fr[15]+fr[14]+fr[13]-1.0*(fr[11]+fr[10]+fr[8]+fr[5])+fr[4]+fr[3]+fr[2]-1.0*fr[0]); 
    f1Quad[8] = -0.25*(fr[31]-1.0*(fr[29]+fr[28]+fr[27])+fr[26]+fr[23]+fr[21]+fr[20]-1.0*(fr[18]+fr[17]+fr[16]+fr[12])+fr[9]+fr[7]+fr[6]-1.0*fr[1]); 
    limQuad[8] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*(alpha[3]+alpha[2]) > 0) {
    f0Quad[9] = 0.25*(fl[30]+fl[25]-1.0*(fl[24]+fl[22])+fl[19]-1.0*(fl[15]+fl[14])+fl[13]+fl[11]-1.0*(fl[10]+fl[8])+fl[5]-1.0*(fl[4]+fl[3])+fl[2]+fl[0]); 
    f1Quad[9] = -0.25*(fl[31]+fl[29]-1.0*(fl[28]+fl[27])+fl[26]-1.0*(fl[23]+fl[21])+fl[20]+fl[18]-1.0*(fl[17]+fl[16])+fl[12]-1.0*(fl[9]+fl[7])+fl[6]+fl[1]); 
    limQuad[9] = fl[0]/cflL; 
  } else {
    f0Quad[9] = 0.25*(fr[30]+fr[25]-1.0*(fr[24]+fr[22])+fr[19]-1.0*(fr[15]+fr[14])+fr[13]+fr[11]-1.0*(fr[10]+fr[8])+fr[5]-1.0*(fr[4]+fr[3])+fr[2]+fr[0]); 
    f1Quad[9] = 0.25*(fr[31]+fr[29]-1.0*(fr[28]+fr[27])+fr[26]-1.0*(fr[23]+fr[21])+fr[20]+fr[18]-1.0*(fr[17]+fr[16])+fr[12]-1.0*(fr[9]+fr[7])+fr[6]+fr[1]); 
    limQuad[9] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[10] = 0.25*(fl[30]-1.0*fl[25]+fl[24]-1.0*fl[22]+fl[19]-1.0*fl[15]+fl[14]-1.0*(fl[13]+fl[11])+fl[10]-1.0*fl[8]+fl[5]-1.0*fl[4]+fl[3]-1.0*fl[2]+fl[0]); 
    f1Quad[10] = -0.25*(fl[31]-1.0*fl[29]+fl[28]-1.0*fl[27]+fl[26]-1.0*fl[23]+fl[21]-1.0*(fl[20]+fl[18])+fl[17]-1.0*fl[16]+fl[12]-1.0*fl[9]+fl[7]-1.0*fl[6]+fl[1]); 
    limQuad[10] = fl[0]/cflL; 
  } else {
    f0Quad[10] = 0.25*(fr[30]-1.0*fr[25]+fr[24]-1.0*fr[22]+fr[19]-1.0*fr[15]+fr[14]-1.0*(fr[13]+fr[11])+fr[10]-1.0*fr[8]+fr[5]-1.0*fr[4]+fr[3]-1.0*fr[2]+fr[0]); 
    f1Quad[10] = 0.25*(fr[31]-1.0*fr[29]+fr[28]-1.0*fr[27]+fr[26]-1.0*fr[23]+fr[21]-1.0*(fr[20]+fr[18])+fr[17]-1.0*fr[16]+fr[12]-1.0*fr[9]+fr[7]-1.0*fr[6]+fr[1]); 
    limQuad[10] = fr[0]/cflR; 
  }
  if(0.25*(alpha[2]+alpha[0])-0.25*alpha[3] > 0) {
    f0Quad[11] = -0.25*(fl[30]+fl[25]+fl[24]-1.0*fl[22]+fl[19]+fl[15]-1.0*(fl[14]+fl[13])+fl[11]+fl[10]-1.0*(fl[8]+fl[5])+fl[4]-1.0*(fl[3]+fl[2]+fl[0])); 
    f1Quad[11] = 0.25*(fl[31]+fl[29]+fl[28]-1.0*fl[27]+fl[26]+fl[23]-1.0*(fl[21]+fl[20])+fl[18]+fl[17]-1.0*(fl[16]+fl[12])+fl[9]-1.0*(fl[7]+fl[6]+fl[1])); 
    limQuad[11] = fl[0]/cflL; 
  } else {
    f0Quad[11] = -0.25*(fr[30]+fr[25]+fr[24]-1.0*fr[22]+fr[19]+fr[15]-1.0*(fr[14]+fr[13])+fr[11]+fr[10]-1.0*(fr[8]+fr[5])+fr[4]-1.0*(fr[3]+fr[2]+fr[0])); 
    f1Quad[11] = -0.25*(fr[31]+fr[29]+fr[28]-1.0*fr[27]+fr[26]+fr[23]-1.0*(fr[21]+fr[20])+fr[18]+fr[17]-1.0*(fr[16]+fr[12])+fr[9]-1.0*(fr[7]+fr[6]+fr[1])); 
    limQuad[11] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[12] = 0.25*(fl[30]-1.0*(fl[25]+fl[24])+fl[22]+fl[19]+fl[15]-1.0*(fl[14]+fl[13]+fl[11]+fl[10])+fl[8]+fl[5]+fl[4]-1.0*(fl[3]+fl[2])+fl[0]); 
    f1Quad[12] = -0.25*(fl[31]-1.0*(fl[29]+fl[28])+fl[27]+fl[26]+fl[23]-1.0*(fl[21]+fl[20]+fl[18]+fl[17])+fl[16]+fl[12]+fl[9]-1.0*(fl[7]+fl[6])+fl[1]); 
    limQuad[12] = fl[0]/cflL; 
  } else {
    f0Quad[12] = 0.25*(fr[30]-1.0*(fr[25]+fr[24])+fr[22]+fr[19]+fr[15]-1.0*(fr[14]+fr[13]+fr[11]+fr[10])+fr[8]+fr[5]+fr[4]-1.0*(fr[3]+fr[2])+fr[0]); 
    f1Quad[12] = 0.25*(fr[31]-1.0*(fr[29]+fr[28])+fr[27]+fr[26]+fr[23]-1.0*(fr[21]+fr[20]+fr[18]+fr[17])+fr[16]+fr[12]+fr[9]-1.0*(fr[7]+fr[6])+fr[1]); 
    limQuad[12] = fr[0]/cflR; 
  }
  if(0.25*alpha[3]-0.25*alpha[2]+0.25*alpha[0] > 0) {
    f0Quad[13] = -0.25*(fl[30]+fl[25]-1.0*fl[24]+fl[22]+fl[19]-1.0*fl[15]+fl[14]-1.0*fl[13]+fl[11]-1.0*fl[10]+fl[8]-1.0*(fl[5]+fl[4])+fl[3]-1.0*(fl[2]+fl[0])); 
    f1Quad[13] = 0.25*(fl[31]+fl[29]-1.0*fl[28]+fl[27]+fl[26]-1.0*fl[23]+fl[21]-1.0*fl[20]+fl[18]-1.0*fl[17]+fl[16]-1.0*(fl[12]+fl[9])+fl[7]-1.0*(fl[6]+fl[1])); 
    limQuad[13] = fl[0]/cflL; 
  } else {
    f0Quad[13] = -0.25*(fr[30]+fr[25]-1.0*fr[24]+fr[22]+fr[19]-1.0*fr[15]+fr[14]-1.0*fr[13]+fr[11]-1.0*fr[10]+fr[8]-1.0*(fr[5]+fr[4])+fr[3]-1.0*(fr[2]+fr[0])); 
    f1Quad[13] = -0.25*(fr[31]+fr[29]-1.0*fr[28]+fr[27]+fr[26]-1.0*fr[23]+fr[21]-1.0*fr[20]+fr[18]-1.0*fr[17]+fr[16]-1.0*(fr[12]+fr[9])+fr[7]-1.0*(fr[6]+fr[1])); 
    limQuad[13] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[14] = -0.25*(fl[30]-1.0*fl[25]+fl[24]+fl[22]+fl[19]-1.0*(fl[15]+fl[14])+fl[13]-1.0*fl[11]+fl[10]+fl[8]-1.0*(fl[5]+fl[4]+fl[3])+fl[2]-1.0*fl[0]); 
    f1Quad[14] = 0.25*(fl[31]-1.0*fl[29]+fl[28]+fl[27]+fl[26]-1.0*(fl[23]+fl[21])+fl[20]-1.0*fl[18]+fl[17]+fl[16]-1.0*(fl[12]+fl[9]+fl[7])+fl[6]-1.0*fl[1]); 
    limQuad[14] = fl[0]/cflL; 
  } else {
    f0Quad[14] = -0.25*(fr[30]-1.0*fr[25]+fr[24]+fr[22]+fr[19]-1.0*(fr[15]+fr[14])+fr[13]-1.0*fr[11]+fr[10]+fr[8]-1.0*(fr[5]+fr[4]+fr[3])+fr[2]-1.0*fr[0]); 
    f1Quad[14] = -0.25*(fr[31]-1.0*fr[29]+fr[28]+fr[27]+fr[26]-1.0*(fr[23]+fr[21])+fr[20]-1.0*fr[18]+fr[17]+fr[16]-1.0*(fr[12]+fr[9]+fr[7])+fr[6]-1.0*fr[1]); 
    limQuad[14] = fr[0]/cflR; 
  }
  if(0.25*(alpha[3]+alpha[2]+alpha[0]) > 0) {
    f0Quad[15] = 0.25*(fl[30]+fl[25]+fl[24]+fl[22]+fl[19]+fl[15]+fl[14]+fl[13]+fl[11]+fl[10]+fl[8]+fl[5]+fl[4]+fl[3]+fl[2]+fl[0]); 
    f1Quad[15] = -0.25*(fl[31]+fl[29]+fl[28]+fl[27]+fl[26]+fl[23]+fl[21]+fl[20]+fl[18]+fl[17]+fl[16]+fl[12]+fl[9]+fl[7]+fl[6]+fl[1]); 
    limQuad[15] = fl[0]/cflL; 
  } else {
    f0Quad[15] = 0.25*(fr[30]+fr[25]+fr[24]+fr[22]+fr[19]+fr[15]+fr[14]+fr[13]+fr[11]+fr[10]+fr[8]+fr[5]+fr[4]+fr[3]+fr[2]+fr[0]); 
    f1Quad[15] = 0.25*(fr[31]+fr[29]+fr[28]+fr[27]+fr[26]+fr[23]+fr[21]+fr[20]+fr[18]+fr[17]+fr[16]+fr[12]+fr[9]+fr[7]+fr[6]+fr[1]); 
    limQuad[15] = fr[0]/cflR; 
  }
  double fhat[32]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[2] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[3] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[6] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[7] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[8] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[9] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[10] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[11] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[12] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[13] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[14] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[15] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[16] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[17] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[18] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[19] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[20] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[21] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[22] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[23] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[24] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[25] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[26] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[27] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[28] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[29] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[30] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[31] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[16];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[23]+fhat[21]+fhat[20]+fhat[18]+fhat[17]+fhat[16]+3.0*fhat[1]))-3.0*(fhat[29]+fhat[28]+fhat[27]+fhat[26]+3.0*(fhat[12]+fhat[9]+fhat[7]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[30]-1.732050807568877*(fhat[25]+fhat[24]+fhat[22]+fhat[19]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[2]))+3.0*(fhat[15]+fhat[14]+fhat[13]+fhat[11]+fhat[10]+fhat[8]+3.0*fhat[0]))); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[23]+fhat[21]))+fhat[20]-1.0*fhat[18]+fhat[17]+fhat[16]-3.0*fhat[1]))+3.0*(fhat[29]-1.0*(fhat[28]+fhat[27]+fhat[26]-3.0*(fhat[12]+fhat[9]+fhat[7]-1.0*fhat[6])))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[15]+fhat[14]-1.0*fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]-3.0*fhat[0]))-1.0*(fhat[30]+1.732050807568877*(fhat[25]-1.0*(fhat[24]+fhat[22]+fhat[19])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[2]))))); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[23])+fhat[21]-1.0*fhat[20]+fhat[18]-1.0*fhat[17]+fhat[16]-3.0*fhat[1]))+3.0*((-1.0*fhat[29])+fhat[28]-1.0*(fhat[27]+fhat[26]-3.0*(fhat[12]+fhat[9]-1.0*fhat[7]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[30])+1.732050807568877*(fhat[25]-1.0*fhat[24]+fhat[22]+fhat[19]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[2]))+3.0*(fhat[15]-1.0*fhat[14]+fhat[13]-1.0*fhat[11]+fhat[10]-1.0*fhat[8]+3.0*fhat[0]))); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[23]-1.0*(fhat[21]+fhat[20]+fhat[18]+fhat[17]-1.0*(fhat[16]+3.0*fhat[1]))))+3.0*(fhat[29]+fhat[28]-1.0*(fhat[27]+fhat[26]-3.0*((-1.0*(fhat[12]+fhat[9]))+fhat[7]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(fhat[30]+1.732050807568877*(fhat[25]+fhat[24]-1.0*(fhat[22]+fhat[19])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[2]))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[13]+fhat[11]+fhat[10]-1.0*(fhat[8]+3.0*fhat[0]))))); 
  rCtrl[4] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[23]-1.0*(fhat[21]+fhat[20]-1.0*fhat[18])+fhat[17]-1.0*(fhat[16]+3.0*fhat[1])))+3.0*((-1.0*(fhat[29]+fhat[28]))+fhat[27]-1.0*fhat[26]+3.0*(fhat[12]-1.0*fhat[9]+fhat[7]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[30])+1.732050807568877*(fhat[25]+fhat[24]-1.0*fhat[22]+fhat[19]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[2])))+3.0*((-1.0*fhat[15])+fhat[14]+fhat[13]-1.0*(fhat[11]+fhat[10]-1.0*(fhat[8]+3.0*fhat[0]))))); 
  rCtrl[5] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[23])+fhat[21]-1.0*(fhat[20]+fhat[18]-1.0*fhat[17]+fhat[16])+3.0*fhat[1]))+3.0*(fhat[29]-1.0*fhat[28]+fhat[27]-1.0*fhat[26]+3.0*((-1.0*fhat[12])+fhat[9]-1.0*fhat[7]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[30]+1.732050807568877*(fhat[25]-1.0*fhat[24]+fhat[22]-1.0*fhat[19]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[2]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8])+3.0*fhat[0]))); 
  rCtrl[6] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[23]+fhat[21]))+fhat[20]+fhat[18]-1.0*(fhat[17]+fhat[16]-3.0*fhat[1])))+3.0*((-1.0*fhat[29])+fhat[28]+fhat[27]-1.0*fhat[26]+3.0*((-1.0*fhat[12])+fhat[9]+fhat[7]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[30]-1.732050807568877*(fhat[25]-1.0*(fhat[24]+fhat[22])+fhat[19]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[2])))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]-3.0*fhat[0])))); 
  rCtrl[7] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[23]+fhat[21]+fhat[20]-1.0*(fhat[18]+fhat[17]+fhat[16]+3.0*fhat[1])))+3.0*(fhat[29]+fhat[28]+fhat[27]-1.0*fhat[26]+3.0*(fhat[12]-1.0*(fhat[9]+fhat[7]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*(fhat[15]+fhat[14]+fhat[13]-1.0*fhat[11]))+fhat[10]+fhat[8]+3.0*fhat[0])-1.0*(fhat[30]+1.732050807568877*(fhat[25]+fhat[24]+fhat[22]-1.0*fhat[19]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[2])))))); 
  rCtrl[8] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[23]+fhat[21]+fhat[20]-1.0*(fhat[18]+fhat[17]+fhat[16]+3.0*fhat[1])))+3.0*(3.0*((-1.0*fhat[12])+fhat[9]+fhat[7]+fhat[6])-1.0*(fhat[29]+fhat[28]+fhat[27]-1.0*fhat[26]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[30])+1.732050807568877*(fhat[25]+fhat[24]+fhat[22]-1.0*fhat[19]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[2])))+3.0*((-1.0*(fhat[15]+fhat[14]+fhat[13]-1.0*fhat[11]))+fhat[10]+fhat[8]+3.0*fhat[0]))); 
  rCtrl[9] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[23]+fhat[21]))+fhat[20]+fhat[18]-1.0*(fhat[17]+fhat[16]-3.0*fhat[1])))+3.0*(fhat[29]-1.0*(fhat[28]+fhat[27]-1.0*fhat[26])+3.0*(fhat[12]-1.0*(fhat[9]+fhat[7]-1.0*fhat[6])))))/(72.0*EPSILON+1.414213562373095*(fhat[30]+1.732050807568877*(fhat[25]-1.0*(fhat[24]+fhat[22])+fhat[19]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[2])))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]-3.0*fhat[0])))); 
  rCtrl[10] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[23])+fhat[21]-1.0*(fhat[20]+fhat[18]-1.0*fhat[17]+fhat[16])+3.0*fhat[1]))+3.0*((-1.0*fhat[29])+fhat[28]-1.0*fhat[27]+fhat[26]+3.0*(fhat[12]-1.0*fhat[9]+fhat[7]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[30]-1.732050807568877*(fhat[25]-1.0*fhat[24]+fhat[22]-1.0*fhat[19]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[2]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8])+3.0*fhat[0]))); 
  rCtrl[11] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[23]-1.0*(fhat[21]+fhat[20]-1.0*fhat[18])+fhat[17]-1.0*(fhat[16]+3.0*fhat[1])))+3.0*(fhat[29]+fhat[28]-1.0*fhat[27]+fhat[26]+3.0*((-1.0*fhat[12])+fhat[9]-1.0*(fhat[7]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*fhat[15])+fhat[14]+fhat[13]-1.0*(fhat[11]+fhat[10]-1.0*(fhat[8]+3.0*fhat[0])))-1.0*(fhat[30]+1.732050807568877*(fhat[25]+fhat[24]-1.0*fhat[22]+fhat[19]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[2])))))); 
  rCtrl[12] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[23]-1.0*(fhat[21]+fhat[20]+fhat[18]+fhat[17]-1.0*(fhat[16]+3.0*fhat[1]))))+3.0*((-1.0*(fhat[29]+fhat[28]))+fhat[27]+fhat[26]+3.0*(fhat[12]+fhat[9]-1.0*(fhat[7]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(fhat[30]-1.732050807568877*(fhat[25]+fhat[24]-1.0*(fhat[22]+fhat[19])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[2]))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[13]+fhat[11]+fhat[10]-1.0*(fhat[8]+3.0*fhat[0]))))); 
  rCtrl[13] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[23])+fhat[21]-1.0*fhat[20]+fhat[18]-1.0*fhat[17]+fhat[16]-3.0*fhat[1]))+3.0*(fhat[29]-1.0*fhat[28]+fhat[27]+fhat[26]+3.0*((-1.0*(fhat[12]+fhat[9]))+fhat[7]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[15]-1.0*fhat[14]+fhat[13]-1.0*fhat[11]+fhat[10]-1.0*fhat[8]+3.0*fhat[0])-1.0*(fhat[30]+1.732050807568877*(fhat[25]-1.0*fhat[24]+fhat[22]+fhat[19]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[2]))))); 
  rCtrl[14] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[23]+fhat[21]))+fhat[20]-1.0*fhat[18]+fhat[17]+fhat[16]-3.0*fhat[1]))+3.0*((-1.0*fhat[29])+fhat[28]+fhat[27]+fhat[26]-3.0*(fhat[12]+fhat[9]+fhat[7]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[30])+1.732050807568877*(fhat[25]-1.0*(fhat[24]+fhat[22]+fhat[19])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[2]))+3.0*(fhat[15]+fhat[14]-1.0*fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]-3.0*fhat[0])))); 
  rCtrl[15] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[23]+fhat[21]+fhat[20]+fhat[18]+fhat[17]+fhat[16]+3.0*fhat[1]))+3.0*(fhat[29]+fhat[28]+fhat[27]+fhat[26]+3.0*(fhat[12]+fhat[9]+fhat[7]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[30]+1.732050807568877*(fhat[25]+fhat[24]+fhat[22]+fhat[19]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[2]))+3.0*(fhat[15]+fhat[14]+fhat[13]+fhat[11]+fhat[10]+fhat[8]+3.0*fhat[0]))); 
  double fhatCtrl[16];  // fhatCtrl = anti-limited fhat evaluated at each control node on x surface 
  fhatCtrl[0] = -0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[24]+fhat[22]+fhat[19]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[2]))-1.0*(fhat[30]+3.0*(fhat[15]+fhat[14]+fhat[13]+fhat[11]+fhat[10]+fhat[8]+3.0*fhat[0])))*limTheta(rCtrl[0],-1.0,EPSILON); 
  fhatCtrl[1] = -0.01964185503295965*(fhat[30]+1.732050807568877*(fhat[25]-1.0*(fhat[24]+fhat[22]+fhat[19])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[2]))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[13]-1.0*fhat[11]+fhat[10]+fhat[8]-3.0*fhat[0]))*limTheta(rCtrl[1],-1.0,EPSILON); 
  fhatCtrl[2] = 0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*fhat[24]+fhat[22]+fhat[19]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[2]))-1.0*(fhat[30]+3.0*((-1.0*fhat[15])+fhat[14]-1.0*fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8]-3.0*fhat[0])))*limTheta(rCtrl[2],-1.0,EPSILON); 
  fhatCtrl[3] = 0.01964185503295965*(fhat[30]+1.732050807568877*(fhat[25]+fhat[24]-1.0*(fhat[22]+fhat[19])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[2]))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[13]+fhat[11]+fhat[10]-1.0*(fhat[8]+3.0*fhat[0]))))*limTheta(rCtrl[3],-1.0,EPSILON); 
  fhatCtrl[4] = 0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[24]-1.0*fhat[22]+fhat[19]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[2])))-1.0*(fhat[30]+3.0*(fhat[15]-1.0*(fhat[14]+fhat[13]-1.0*fhat[11])+fhat[10]-1.0*(fhat[8]+3.0*fhat[0]))))*limTheta(rCtrl[4],-1.0,EPSILON); 
  fhatCtrl[5] = 0.01964185503295965*(fhat[30]+1.732050807568877*(fhat[25]-1.0*fhat[24]+fhat[22]-1.0*fhat[19]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[2]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8])+3.0*fhat[0]))*limTheta(rCtrl[5],-1.0,EPSILON); 
  fhatCtrl[6] = -0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*(fhat[24]+fhat[22])+fhat[19]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[2])))-1.0*(fhat[30]+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]-3.0*fhat[0]))))*limTheta(rCtrl[6],-1.0,EPSILON); 
  fhatCtrl[7] = -0.01964185503295965*(fhat[30]+1.732050807568877*(fhat[25]+fhat[24]+fhat[22]-1.0*fhat[19]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[2])))+3.0*(fhat[15]+fhat[14]+fhat[13]-1.0*(fhat[11]+fhat[10]+fhat[8]+3.0*fhat[0])))*limTheta(rCtrl[7],-1.0,EPSILON); 
  fhatCtrl[8] = 0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[24]+fhat[22]-1.0*fhat[19]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[2])))-1.0*(fhat[30]+3.0*(fhat[15]+fhat[14]+fhat[13]-1.0*(fhat[11]+fhat[10]+fhat[8]+3.0*fhat[0]))))*limTheta(rCtrl[8],-1.0,EPSILON); 
  fhatCtrl[9] = 0.01964185503295965*(fhat[30]+1.732050807568877*(fhat[25]-1.0*(fhat[24]+fhat[22])+fhat[19]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[2])))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[13]+fhat[11]-1.0*(fhat[10]+fhat[8]-3.0*fhat[0])))*limTheta(rCtrl[9],-1.0,EPSILON); 
  fhatCtrl[10] = -0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*fhat[24]+fhat[22]-1.0*fhat[19]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[2]))-1.0*(fhat[30]+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8])+3.0*fhat[0])))*limTheta(rCtrl[10],-1.0,EPSILON); 
  fhatCtrl[11] = -0.01964185503295965*(fhat[30]+1.732050807568877*(fhat[25]+fhat[24]-1.0*fhat[22]+fhat[19]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[2])))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[13]-1.0*fhat[11])+fhat[10]-1.0*(fhat[8]+3.0*fhat[0])))*limTheta(rCtrl[11],-1.0,EPSILON); 
  fhatCtrl[12] = -0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[24]-1.0*(fhat[22]+fhat[19])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[2]))-1.0*(fhat[30]+3.0*(fhat[15]-1.0*(fhat[14]+fhat[13]+fhat[11]+fhat[10]-1.0*(fhat[8]+3.0*fhat[0])))))*limTheta(rCtrl[12],-1.0,EPSILON); 
  fhatCtrl[13] = -0.01964185503295965*(fhat[30]+1.732050807568877*(fhat[25]-1.0*fhat[24]+fhat[22]+fhat[19]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[2]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*fhat[13]+fhat[11]-1.0*fhat[10]+fhat[8]-3.0*fhat[0]))*limTheta(rCtrl[13],-1.0,EPSILON); 
  fhatCtrl[14] = 0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*(fhat[24]+fhat[22]+fhat[19])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[2]))-1.0*(fhat[30]+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[13]-1.0*fhat[11]+fhat[10]+fhat[8]-3.0*fhat[0])))*limTheta(rCtrl[14],-1.0,EPSILON); 
  fhatCtrl[15] = 0.01964185503295965*(fhat[30]+1.732050807568877*(fhat[25]+fhat[24]+fhat[22]+fhat[19]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[2]))+3.0*(fhat[15]+fhat[14]+fhat[13]+fhat[11]+fhat[10]+fhat[8]+3.0*fhat[0]))*limTheta(rCtrl[15],-1.0,EPSILON); 
  double fhatAL[16];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.25*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.4330127018922193*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[5] = 0.75*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[6] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[8] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[9] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[10] = 0.75*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[11] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[12] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[13] = 1.299038105676658*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[14] = 1.299038105676658*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[15] = 2.25*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[16]; 
  fhatALQuad[0] = fmin(0.25*((-0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7])))-1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[0]); 
  fhatALQuad[1] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*((-0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+0.5773502691896258*(fhatAL[13]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[1]); 
  fhatALQuad[2] = fmin(0.25*(0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[2]); 
  fhatALQuad[3] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[3]); 
  fhatALQuad[4] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*(0.5773502691896258*fhatAL[3]-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[4]+fhatAL[2]+fhatAL[1]))-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[4]); 
  fhatALQuad[5] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]))+0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[5]); 
  fhatALQuad[6] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2])-0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[6]); 
  fhatALQuad[7] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[7]); 
  fhatALQuad[8] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]+fhatAL[1]))+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[8]); 
  fhatALQuad[9] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]))+0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[9]); 
  fhatALQuad[10] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[10]); 
  fhatALQuad[11] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[4]+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[11]); 
  fhatALQuad[12] = fmin(0.25*((-0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5])))+1.732050807568877*(0.5773502691896258*(fhatAL[4]+fhatAL[3])-0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[12]); 
  fhatALQuad[13] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*(0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1])-0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[13]); 
  fhatALQuad[14] = fmin(0.25*(0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[1])+0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[14]); 
  fhatALQuad[15] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[15]); 
  fhatAL[0] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[5] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[6] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[8] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[9] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[10] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[11] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[12] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[13] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[14] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[15] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.1767766952966368*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_x; 
  incr[1] = -0.3061862178478971*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[0]*fhatAL[0])*dfac_x; 
  incr[2] = 0.1767766952966368*(alpha[3]*fhatAL[6]+alpha[2]*fhatAL[5]+alpha[0]*fhatAL[1])*dfac_x; 
  incr[3] = 0.1767766952966368*(alpha[3]*fhatAL[7]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_x; 
  incr[4] = 0.1767766952966368*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_x; 
  incr[5] = 0.1767766952966368*(alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[0]*fhatAL[4])*dfac_x; 
  incr[6] = -0.3061862178478971*(alpha[3]*fhatAL[6]+alpha[2]*fhatAL[5]+alpha[0]*fhatAL[1])*dfac_x; 
  incr[7] = -0.3061862178478971*(alpha[3]*fhatAL[7]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_x; 
  incr[8] = 0.1767766952966368*(alpha[3]*fhatAL[11]+alpha[0]*fhatAL[5]+fhatAL[1]*alpha[2])*dfac_x; 
  incr[9] = -0.3061862178478971*(alpha[2]*fhatAL[7]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_x; 
  incr[10] = 0.1767766952966368*(alpha[2]*fhatAL[11]+alpha[0]*fhatAL[6]+fhatAL[1]*alpha[3])*dfac_x; 
  incr[11] = 0.1767766952966368*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_x; 
  incr[12] = -0.3061862178478971*(alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[0]*fhatAL[4])*dfac_x; 
  incr[13] = 0.1767766952966368*(alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[0]*fhatAL[8])*dfac_x; 
  incr[14] = 0.1767766952966368*(alpha[3]*fhatAL[14]+alpha[0]*fhatAL[9]+alpha[2]*fhatAL[4])*dfac_x; 
  incr[15] = 0.1767766952966368*(alpha[2]*fhatAL[14]+alpha[0]*fhatAL[10]+alpha[3]*fhatAL[4])*dfac_x; 
  incr[16] = -0.3061862178478971*(alpha[3]*fhatAL[11]+alpha[0]*fhatAL[5]+fhatAL[1]*alpha[2])*dfac_x; 
  incr[17] = -0.3061862178478971*(alpha[2]*fhatAL[11]+alpha[0]*fhatAL[6]+fhatAL[1]*alpha[3])*dfac_x; 
  incr[18] = -0.3061862178478971*(alpha[0]*fhatAL[7]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_x; 
  incr[19] = 0.1767766952966368*(alpha[0]*fhatAL[11]+alpha[2]*fhatAL[6]+alpha[3]*fhatAL[5])*dfac_x; 
  incr[20] = -0.3061862178478971*(alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[0]*fhatAL[8])*dfac_x; 
  incr[21] = -0.3061862178478971*(alpha[3]*fhatAL[14]+alpha[0]*fhatAL[9]+alpha[2]*fhatAL[4])*dfac_x; 
  incr[22] = 0.1767766952966368*(alpha[3]*fhatAL[15]+alpha[0]*fhatAL[12]+alpha[2]*fhatAL[8])*dfac_x; 
  incr[23] = -0.3061862178478971*(alpha[2]*fhatAL[14]+alpha[0]*fhatAL[10]+alpha[3]*fhatAL[4])*dfac_x; 
  incr[24] = 0.1767766952966368*(alpha[2]*fhatAL[15]+alpha[0]*fhatAL[13]+alpha[3]*fhatAL[8])*dfac_x; 
  incr[25] = 0.1767766952966368*(alpha[0]*fhatAL[14]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9])*dfac_x; 
  incr[26] = -0.3061862178478971*(alpha[0]*fhatAL[11]+alpha[2]*fhatAL[6]+alpha[3]*fhatAL[5])*dfac_x; 
  incr[27] = -0.3061862178478971*(alpha[3]*fhatAL[15]+alpha[0]*fhatAL[12]+alpha[2]*fhatAL[8])*dfac_x; 
  incr[28] = -0.3061862178478971*(alpha[2]*fhatAL[15]+alpha[0]*fhatAL[13]+alpha[3]*fhatAL[8])*dfac_x; 
  incr[29] = -0.3061862178478971*(alpha[0]*fhatAL[14]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9])*dfac_x; 
  incr[30] = 0.1767766952966368*(alpha[0]*fhatAL[15]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12])*dfac_x; 
  incr[31] = -0.3061862178478971*(alpha[0]*fhatAL[15]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 
  outr[28] += incr[28]; 
  outr[29] += incr[29]; 
  outr[30] += incr[30]; 
  outr[31] += incr[31]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += -1.0*incr[19]; 
  outl[20] += incr[20]; 
  outl[21] += incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += incr[26]; 
  outl[27] += incr[27]; 
  outl[28] += incr[28]; 
  outl[29] += incr[29]; 
  outl[30] += -1.0*incr[30]; 
  outl[31] += incr[31]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity3x2vSer_Y_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_z = 2.0/dxv[2]; 
  double dfac_v = 2.0/dxv[3]; 
  double dfac_m = 2.0/dxv[4]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wz = w[2]; 
  double wv = w[3]; 
  double wm = w[4]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[32]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.0625*(2.828427124746191*BdriftY[0]*m_*wv2+BmagInv[0]*dfac_x*(1.732050807568877*(Bmag[1]*wm+Phi[1]*q_)-3.0*Phi[4]*q_)))/q_; 

  double alpha[16]; 
  alpha[0] = (0.5*(2.828427124746191*BdriftY[0]*m_*wv2+BmagInv[0]*dfac_x*(1.732050807568877*(Bmag[1]*wm+Phi[1]*q_)-3.0*Phi[4]*q_)))/q_; 
  alpha[1] = (0.5*(2.828427124746191*BdriftY[1]*m_*wv2+BmagInv[1]*dfac_x*(1.732050807568877*(Bmag[1]*wm+Phi[1]*q_)-3.0*Phi[4]*q_)))/q_; 
  alpha[2] = 0.5*BmagInv[0]*(1.732050807568877*Phi[5]-3.0*Phi[7])*dfac_x; 
  alpha[3] = (0.8164965809277261*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alpha[4] = (0.5*BmagInv[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alpha[5] = 0.5*BmagInv[1]*(1.732050807568877*Phi[5]-3.0*Phi[7])*dfac_x; 
  alpha[6] = (0.8164965809277261*BdriftY[1]*m_*wv)/(dfac_v*q_); 
  alpha[8] = (0.5*Bmag[1]*BmagInv[1]*dfac_x)/(dfac_m*q_); 
  double f0Quad[16]; 
  double f1Quad[16]; 
  double limQuad[16]; 
  // determine upwinding at each surface quadrature node 
  if(0.25*(alpha[8]+alpha[6]+alpha[5])-0.25*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[0] = 0.25*(fl[29]-1.0*(fl[25]+fl[23]+fl[21]+fl[18])+fl[15]+fl[14]+fl[12]+fl[11]+fl[9]+fl[7]-1.0*(fl[5]+fl[4]+fl[3]+fl[1])+fl[0]); 
    f1Quad[0] = -0.25*(fl[31]-1.0*(fl[30]+fl[28]+fl[27]+fl[26])+fl[24]+fl[22]+fl[20]+fl[19]+fl[17]+fl[16]-1.0*(fl[13]+fl[10]+fl[8]+fl[6])+fl[2]); 
    limQuad[0] = fl[0]/cflL; 
  } else {
    f0Quad[0] = 0.25*(fr[29]-1.0*(fr[25]+fr[23]+fr[21]+fr[18])+fr[15]+fr[14]+fr[12]+fr[11]+fr[9]+fr[7]-1.0*(fr[5]+fr[4]+fr[3]+fr[1])+fr[0]); 
    f1Quad[0] = 0.25*(fr[31]-1.0*(fr[30]+fr[28]+fr[27]+fr[26])+fr[24]+fr[22]+fr[20]+fr[19]+fr[17]+fr[16]-1.0*(fr[13]+fr[10]+fr[8]+fr[6])+fr[2]); 
    limQuad[0] = fr[0]/cflR; 
  }
  if(0.25*(alpha[1]+alpha[0])-0.25*(alpha[8]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]) > 0) {
    f0Quad[1] = -0.25*(fl[29]+fl[25]-1.0*(fl[23]+fl[21]+fl[18]+fl[15]+fl[14])+fl[12]-1.0*fl[11]+fl[9]+fl[7]+fl[5]+fl[4]+fl[3]-1.0*(fl[1]+fl[0])); 
    f1Quad[1] = 0.25*(fl[31]+fl[30]-1.0*(fl[28]+fl[27]+fl[26]+fl[24]+fl[22])+fl[20]-1.0*fl[19]+fl[17]+fl[16]+fl[13]+fl[10]+fl[8]-1.0*(fl[6]+fl[2])); 
    limQuad[1] = fl[0]/cflL; 
  } else {
    f0Quad[1] = -0.25*(fr[29]+fr[25]-1.0*(fr[23]+fr[21]+fr[18]+fr[15]+fr[14])+fr[12]-1.0*fr[11]+fr[9]+fr[7]+fr[5]+fr[4]+fr[3]-1.0*(fr[1]+fr[0])); 
    f1Quad[1] = -0.25*(fr[31]+fr[30]-1.0*(fr[28]+fr[27]+fr[26]+fr[24]+fr[22])+fr[20]-1.0*fr[19]+fr[17]+fr[16]+fr[13]+fr[10]+fr[8]-1.0*(fr[6]+fr[2])); 
    limQuad[1] = fr[0]/cflR; 
  }
  if(0.25*(alpha[8]+alpha[6])-0.25*(alpha[5]+alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[2] = -0.25*(fl[29]-1.0*fl[25]+fl[23]-1.0*(fl[21]+fl[18]+fl[15])+fl[14]-1.0*fl[12]+fl[11]-1.0*fl[9]+fl[7]+fl[5]+fl[4]-1.0*fl[3]+fl[1]-1.0*fl[0]); 
    f1Quad[2] = 0.25*(fl[31]-1.0*fl[30]+fl[28]-1.0*(fl[27]+fl[26]+fl[24])+fl[22]-1.0*fl[20]+fl[19]-1.0*fl[17]+fl[16]+fl[13]+fl[10]-1.0*fl[8]+fl[6]-1.0*fl[2]); 
    limQuad[2] = fl[0]/cflL; 
  } else {
    f0Quad[2] = -0.25*(fr[29]-1.0*fr[25]+fr[23]-1.0*(fr[21]+fr[18]+fr[15])+fr[14]-1.0*fr[12]+fr[11]-1.0*fr[9]+fr[7]+fr[5]+fr[4]-1.0*fr[3]+fr[1]-1.0*fr[0]); 
    f1Quad[2] = -0.25*(fr[31]-1.0*fr[30]+fr[28]-1.0*(fr[27]+fr[26]+fr[24])+fr[22]-1.0*fr[20]+fr[19]-1.0*fr[17]+fr[16]+fr[13]+fr[10]-1.0*fr[8]+fr[6]-1.0*fr[2]); 
    limQuad[2] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[8]+alpha[6]))+0.25*alpha[5]-0.25*(alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[3] = 0.25*(fl[29]+fl[25]+fl[23]-1.0*(fl[21]+fl[18])+fl[15]-1.0*(fl[14]+fl[12]+fl[11]+fl[9])+fl[7]-1.0*(fl[5]+fl[4])+fl[3]+fl[1]+fl[0]); 
    f1Quad[3] = -0.25*(fl[31]+fl[30]+fl[28]-1.0*(fl[27]+fl[26])+fl[24]-1.0*(fl[22]+fl[20]+fl[19]+fl[17])+fl[16]-1.0*(fl[13]+fl[10])+fl[8]+fl[6]+fl[2]); 
    limQuad[3] = fl[0]/cflL; 
  } else {
    f0Quad[3] = 0.25*(fr[29]+fr[25]+fr[23]-1.0*(fr[21]+fr[18])+fr[15]-1.0*(fr[14]+fr[12]+fr[11]+fr[9])+fr[7]-1.0*(fr[5]+fr[4])+fr[3]+fr[1]+fr[0]); 
    f1Quad[3] = 0.25*(fr[31]+fr[30]+fr[28]-1.0*(fr[27]+fr[26])+fr[24]-1.0*(fr[22]+fr[20]+fr[19]+fr[17])+fr[16]-1.0*(fr[13]+fr[10])+fr[8]+fr[6]+fr[2]); 
    limQuad[3] = fr[0]/cflR; 
  }
  if(0.25*alpha[8]-0.25*alpha[6]+0.25*alpha[5]-0.25*alpha[4]+0.25*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[4] = -0.25*(fl[29]-1.0*(fl[25]+fl[23])+fl[21]-1.0*fl[18]+fl[15]-1.0*(fl[14]+fl[12])+fl[11]+fl[9]-1.0*fl[7]+fl[5]-1.0*fl[4]+fl[3]+fl[1]-1.0*fl[0]); 
    f1Quad[4] = 0.25*(fl[31]-1.0*(fl[30]+fl[28])+fl[27]-1.0*fl[26]+fl[24]-1.0*(fl[22]+fl[20])+fl[19]+fl[17]-1.0*fl[16]+fl[13]-1.0*fl[10]+fl[8]+fl[6]-1.0*fl[2]); 
    limQuad[4] = fl[0]/cflL; 
  } else {
    f0Quad[4] = -0.25*(fr[29]-1.0*(fr[25]+fr[23])+fr[21]-1.0*fr[18]+fr[15]-1.0*(fr[14]+fr[12])+fr[11]+fr[9]-1.0*fr[7]+fr[5]-1.0*fr[4]+fr[3]+fr[1]-1.0*fr[0]); 
    f1Quad[4] = -0.25*(fr[31]-1.0*(fr[30]+fr[28])+fr[27]-1.0*fr[26]+fr[24]-1.0*(fr[22]+fr[20])+fr[19]+fr[17]-1.0*fr[16]+fr[13]-1.0*fr[10]+fr[8]+fr[6]-1.0*fr[2]); 
    limQuad[4] = fr[0]/cflR; 
  }
  if((-0.25*alpha[8])+0.25*alpha[6]-0.25*(alpha[5]+alpha[4])+0.25*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[5] = 0.25*(fl[29]+fl[25]-1.0*fl[23]+fl[21]-1.0*(fl[18]+fl[15])+fl[14]-1.0*(fl[12]+fl[11])+fl[9]-1.0*(fl[7]+fl[5])+fl[4]-1.0*fl[3]+fl[1]+fl[0]); 
    f1Quad[5] = -0.25*(fl[31]+fl[30]-1.0*fl[28]+fl[27]-1.0*(fl[26]+fl[24])+fl[22]-1.0*(fl[20]+fl[19])+fl[17]-1.0*(fl[16]+fl[13])+fl[10]-1.0*fl[8]+fl[6]+fl[2]); 
    limQuad[5] = fl[0]/cflL; 
  } else {
    f0Quad[5] = 0.25*(fr[29]+fr[25]-1.0*fr[23]+fr[21]-1.0*(fr[18]+fr[15])+fr[14]-1.0*(fr[12]+fr[11])+fr[9]-1.0*(fr[7]+fr[5])+fr[4]-1.0*fr[3]+fr[1]+fr[0]); 
    f1Quad[5] = 0.25*(fr[31]+fr[30]-1.0*fr[28]+fr[27]-1.0*(fr[26]+fr[24])+fr[22]-1.0*(fr[20]+fr[19])+fr[17]-1.0*(fr[16]+fr[13])+fr[10]-1.0*fr[8]+fr[6]+fr[2]); 
    limQuad[5] = fr[0]/cflR; 
  }
  if(0.25*alpha[8]-0.25*(alpha[6]+alpha[5]+alpha[4])+0.25*(alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[6] = 0.25*(fl[29]-1.0*fl[25]+fl[23]+fl[21]-1.0*(fl[18]+fl[15]+fl[14])+fl[12]+fl[11]-1.0*(fl[9]+fl[7]+fl[5])+fl[4]+fl[3]-1.0*fl[1]+fl[0]); 
    f1Quad[6] = -0.25*(fl[31]-1.0*fl[30]+fl[28]+fl[27]-1.0*(fl[26]+fl[24]+fl[22])+fl[20]+fl[19]-1.0*(fl[17]+fl[16]+fl[13])+fl[10]+fl[8]-1.0*fl[6]+fl[2]); 
    limQuad[6] = fl[0]/cflL; 
  } else {
    f0Quad[6] = 0.25*(fr[29]-1.0*fr[25]+fr[23]+fr[21]-1.0*(fr[18]+fr[15]+fr[14])+fr[12]+fr[11]-1.0*(fr[9]+fr[7]+fr[5])+fr[4]+fr[3]-1.0*fr[1]+fr[0]); 
    f1Quad[6] = 0.25*(fr[31]-1.0*fr[30]+fr[28]+fr[27]-1.0*(fr[26]+fr[24]+fr[22])+fr[20]+fr[19]-1.0*(fr[17]+fr[16]+fr[13])+fr[10]+fr[8]-1.0*fr[6]+fr[2]); 
    limQuad[6] = fr[0]/cflR; 
  }
  if((-0.25*alpha[8])+0.25*(alpha[6]+alpha[5])-0.25*alpha[4]+0.25*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[7] = -0.25*(fl[29]+fl[25]+fl[23]+fl[21]-1.0*fl[18]+fl[15]+fl[14]+fl[12]-1.0*(fl[11]+fl[9]+fl[7])+fl[5]-1.0*(fl[4]+fl[3]+fl[1]+fl[0])); 
    f1Quad[7] = 0.25*(fl[31]+fl[30]+fl[28]+fl[27]-1.0*fl[26]+fl[24]+fl[22]+fl[20]-1.0*(fl[19]+fl[17]+fl[16])+fl[13]-1.0*(fl[10]+fl[8]+fl[6]+fl[2])); 
    limQuad[7] = fl[0]/cflL; 
  } else {
    f0Quad[7] = -0.25*(fr[29]+fr[25]+fr[23]+fr[21]-1.0*fr[18]+fr[15]+fr[14]+fr[12]-1.0*(fr[11]+fr[9]+fr[7])+fr[5]-1.0*(fr[4]+fr[3]+fr[1]+fr[0])); 
    f1Quad[7] = -0.25*(fr[31]+fr[30]+fr[28]+fr[27]-1.0*fr[26]+fr[24]+fr[22]+fr[20]-1.0*(fr[19]+fr[17]+fr[16])+fr[13]-1.0*(fr[10]+fr[8]+fr[6]+fr[2])); 
    limQuad[7] = fr[0]/cflR; 
  }
  if((-0.25*alpha[8])+0.25*(alpha[6]+alpha[5]+alpha[4])-0.25*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[8] = -0.25*(fl[29]-1.0*(fl[25]+fl[23]+fl[21])+fl[18]+fl[15]+fl[14]+fl[12]-1.0*(fl[11]+fl[9]+fl[7]+fl[5])+fl[4]+fl[3]+fl[1]-1.0*fl[0]); 
    f1Quad[8] = 0.25*(fl[31]-1.0*(fl[30]+fl[28]+fl[27])+fl[26]+fl[24]+fl[22]+fl[20]-1.0*(fl[19]+fl[17]+fl[16]+fl[13])+fl[10]+fl[8]+fl[6]-1.0*fl[2]); 
    limQuad[8] = fl[0]/cflL; 
  } else {
    f0Quad[8] = -0.25*(fr[29]-1.0*(fr[25]+fr[23]+fr[21])+fr[18]+fr[15]+fr[14]+fr[12]-1.0*(fr[11]+fr[9]+fr[7]+fr[5])+fr[4]+fr[3]+fr[1]-1.0*fr[0]); 
    f1Quad[8] = -0.25*(fr[31]-1.0*(fr[30]+fr[28]+fr[27])+fr[26]+fr[24]+fr[22]+fr[20]-1.0*(fr[19]+fr[17]+fr[16]+fr[13])+fr[10]+fr[8]+fr[6]-1.0*fr[2]); 
    limQuad[8] = fr[0]/cflR; 
  }
  if(0.25*alpha[8]-0.25*(alpha[6]+alpha[5])+0.25*alpha[4]-0.25*(alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[9] = 0.25*(fl[29]+fl[25]-1.0*(fl[23]+fl[21])+fl[18]-1.0*(fl[15]+fl[14])+fl[12]+fl[11]-1.0*(fl[9]+fl[7])+fl[5]-1.0*(fl[4]+fl[3])+fl[1]+fl[0]); 
    f1Quad[9] = -0.25*(fl[31]+fl[30]-1.0*(fl[28]+fl[27])+fl[26]-1.0*(fl[24]+fl[22])+fl[20]+fl[19]-1.0*(fl[17]+fl[16])+fl[13]-1.0*(fl[10]+fl[8])+fl[6]+fl[2]); 
    limQuad[9] = fl[0]/cflL; 
  } else {
    f0Quad[9] = 0.25*(fr[29]+fr[25]-1.0*(fr[23]+fr[21])+fr[18]-1.0*(fr[15]+fr[14])+fr[12]+fr[11]-1.0*(fr[9]+fr[7])+fr[5]-1.0*(fr[4]+fr[3])+fr[1]+fr[0]); 
    f1Quad[9] = 0.25*(fr[31]+fr[30]-1.0*(fr[28]+fr[27])+fr[26]-1.0*(fr[24]+fr[22])+fr[20]+fr[19]-1.0*(fr[17]+fr[16])+fr[13]-1.0*(fr[10]+fr[8])+fr[6]+fr[2]); 
    limQuad[9] = fr[0]/cflR; 
  }
  if((-0.25*alpha[8])+0.25*alpha[6]-0.25*alpha[5]+0.25*alpha[4]-0.25*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[10] = 0.25*(fl[29]-1.0*fl[25]+fl[23]-1.0*fl[21]+fl[18]-1.0*fl[15]+fl[14]-1.0*(fl[12]+fl[11])+fl[9]-1.0*fl[7]+fl[5]-1.0*fl[4]+fl[3]-1.0*fl[1]+fl[0]); 
    f1Quad[10] = -0.25*(fl[31]-1.0*fl[30]+fl[28]-1.0*fl[27]+fl[26]-1.0*fl[24]+fl[22]-1.0*(fl[20]+fl[19])+fl[17]-1.0*fl[16]+fl[13]-1.0*fl[10]+fl[8]-1.0*fl[6]+fl[2]); 
    limQuad[10] = fl[0]/cflL; 
  } else {
    f0Quad[10] = 0.25*(fr[29]-1.0*fr[25]+fr[23]-1.0*fr[21]+fr[18]-1.0*fr[15]+fr[14]-1.0*(fr[12]+fr[11])+fr[9]-1.0*fr[7]+fr[5]-1.0*fr[4]+fr[3]-1.0*fr[1]+fr[0]); 
    f1Quad[10] = 0.25*(fr[31]-1.0*fr[30]+fr[28]-1.0*fr[27]+fr[26]-1.0*fr[24]+fr[22]-1.0*(fr[20]+fr[19])+fr[17]-1.0*fr[16]+fr[13]-1.0*fr[10]+fr[8]-1.0*fr[6]+fr[2]); 
    limQuad[10] = fr[0]/cflR; 
  }
  if(0.25*alpha[8]-0.25*alpha[6]+0.25*(alpha[5]+alpha[4])-0.25*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[11] = -0.25*(fl[29]+fl[25]+fl[23]-1.0*fl[21]+fl[18]+fl[15]-1.0*(fl[14]+fl[12])+fl[11]+fl[9]-1.0*(fl[7]+fl[5])+fl[4]-1.0*(fl[3]+fl[1]+fl[0])); 
    f1Quad[11] = 0.25*(fl[31]+fl[30]+fl[28]-1.0*fl[27]+fl[26]+fl[24]-1.0*(fl[22]+fl[20])+fl[19]+fl[17]-1.0*(fl[16]+fl[13])+fl[10]-1.0*(fl[8]+fl[6]+fl[2])); 
    limQuad[11] = fl[0]/cflL; 
  } else {
    f0Quad[11] = -0.25*(fr[29]+fr[25]+fr[23]-1.0*fr[21]+fr[18]+fr[15]-1.0*(fr[14]+fr[12])+fr[11]+fr[9]-1.0*(fr[7]+fr[5])+fr[4]-1.0*(fr[3]+fr[1]+fr[0])); 
    f1Quad[11] = -0.25*(fr[31]+fr[30]+fr[28]-1.0*fr[27]+fr[26]+fr[24]-1.0*(fr[22]+fr[20])+fr[19]+fr[17]-1.0*(fr[16]+fr[13])+fr[10]-1.0*(fr[8]+fr[6]+fr[2])); 
    limQuad[11] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[8]+alpha[6]))+0.25*(alpha[5]+alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[12] = 0.25*(fl[29]-1.0*(fl[25]+fl[23])+fl[21]+fl[18]+fl[15]-1.0*(fl[14]+fl[12]+fl[11]+fl[9])+fl[7]+fl[5]+fl[4]-1.0*(fl[3]+fl[1])+fl[0]); 
    f1Quad[12] = -0.25*(fl[31]-1.0*(fl[30]+fl[28])+fl[27]+fl[26]+fl[24]-1.0*(fl[22]+fl[20]+fl[19]+fl[17])+fl[16]+fl[13]+fl[10]-1.0*(fl[8]+fl[6])+fl[2]); 
    limQuad[12] = fl[0]/cflL; 
  } else {
    f0Quad[12] = 0.25*(fr[29]-1.0*(fr[25]+fr[23])+fr[21]+fr[18]+fr[15]-1.0*(fr[14]+fr[12]+fr[11]+fr[9])+fr[7]+fr[5]+fr[4]-1.0*(fr[3]+fr[1])+fr[0]); 
    f1Quad[12] = 0.25*(fr[31]-1.0*(fr[30]+fr[28])+fr[27]+fr[26]+fr[24]-1.0*(fr[22]+fr[20]+fr[19]+fr[17])+fr[16]+fr[13]+fr[10]-1.0*(fr[8]+fr[6])+fr[2]); 
    limQuad[12] = fr[0]/cflR; 
  }
  if(0.25*(alpha[8]+alpha[6])-0.25*alpha[5]+0.25*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[13] = -0.25*(fl[29]+fl[25]-1.0*fl[23]+fl[21]+fl[18]-1.0*fl[15]+fl[14]-1.0*fl[12]+fl[11]-1.0*fl[9]+fl[7]-1.0*(fl[5]+fl[4])+fl[3]-1.0*(fl[1]+fl[0])); 
    f1Quad[13] = 0.25*(fl[31]+fl[30]-1.0*fl[28]+fl[27]+fl[26]-1.0*fl[24]+fl[22]-1.0*fl[20]+fl[19]-1.0*fl[17]+fl[16]-1.0*(fl[13]+fl[10])+fl[8]-1.0*(fl[6]+fl[2])); 
    limQuad[13] = fl[0]/cflL; 
  } else {
    f0Quad[13] = -0.25*(fr[29]+fr[25]-1.0*fr[23]+fr[21]+fr[18]-1.0*fr[15]+fr[14]-1.0*fr[12]+fr[11]-1.0*fr[9]+fr[7]-1.0*(fr[5]+fr[4])+fr[3]-1.0*(fr[1]+fr[0])); 
    f1Quad[13] = -0.25*(fr[31]+fr[30]-1.0*fr[28]+fr[27]+fr[26]-1.0*fr[24]+fr[22]-1.0*fr[20]+fr[19]-1.0*fr[17]+fr[16]-1.0*(fr[13]+fr[10])+fr[8]-1.0*(fr[6]+fr[2])); 
    limQuad[13] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[8]+alpha[6]+alpha[5]))+0.25*(alpha[4]+alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[14] = -0.25*(fl[29]-1.0*fl[25]+fl[23]+fl[21]+fl[18]-1.0*(fl[15]+fl[14])+fl[12]-1.0*fl[11]+fl[9]+fl[7]-1.0*(fl[5]+fl[4]+fl[3])+fl[1]-1.0*fl[0]); 
    f1Quad[14] = 0.25*(fl[31]-1.0*fl[30]+fl[28]+fl[27]+fl[26]-1.0*(fl[24]+fl[22])+fl[20]-1.0*fl[19]+fl[17]+fl[16]-1.0*(fl[13]+fl[10]+fl[8])+fl[6]-1.0*fl[2]); 
    limQuad[14] = fl[0]/cflL; 
  } else {
    f0Quad[14] = -0.25*(fr[29]-1.0*fr[25]+fr[23]+fr[21]+fr[18]-1.0*(fr[15]+fr[14])+fr[12]-1.0*fr[11]+fr[9]+fr[7]-1.0*(fr[5]+fr[4]+fr[3])+fr[1]-1.0*fr[0]); 
    f1Quad[14] = -0.25*(fr[31]-1.0*fr[30]+fr[28]+fr[27]+fr[26]-1.0*(fr[24]+fr[22])+fr[20]-1.0*fr[19]+fr[17]+fr[16]-1.0*(fr[13]+fr[10]+fr[8])+fr[6]-1.0*fr[2]); 
    limQuad[14] = fr[0]/cflR; 
  }
  if(0.25*(alpha[8]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[15] = 0.25*(fl[29]+fl[25]+fl[23]+fl[21]+fl[18]+fl[15]+fl[14]+fl[12]+fl[11]+fl[9]+fl[7]+fl[5]+fl[4]+fl[3]+fl[1]+fl[0]); 
    f1Quad[15] = -0.25*(fl[31]+fl[30]+fl[28]+fl[27]+fl[26]+fl[24]+fl[22]+fl[20]+fl[19]+fl[17]+fl[16]+fl[13]+fl[10]+fl[8]+fl[6]+fl[2]); 
    limQuad[15] = fl[0]/cflL; 
  } else {
    f0Quad[15] = 0.25*(fr[29]+fr[25]+fr[23]+fr[21]+fr[18]+fr[15]+fr[14]+fr[12]+fr[11]+fr[9]+fr[7]+fr[5]+fr[4]+fr[3]+fr[1]+fr[0]); 
    f1Quad[15] = 0.25*(fr[31]+fr[30]+fr[28]+fr[27]+fr[26]+fr[24]+fr[22]+fr[20]+fr[19]+fr[17]+fr[16]+fr[13]+fr[10]+fr[8]+fr[6]+fr[2]); 
    limQuad[15] = fr[0]/cflR; 
  }
  double fhat[32]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[3] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[6] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[7] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[8] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[9] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[10] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[11] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[12] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[13] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[14] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[15] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[16] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[17] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[18] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[19] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[20] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[21] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[22] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[23] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[24] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[25] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[26] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[27] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[28] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[29] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[30] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[31] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[16];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]+fhat[22]+fhat[20]+fhat[19]+fhat[17]+fhat[16]+3.0*fhat[2]))-3.0*(fhat[30]+fhat[28]+fhat[27]+fhat[26]+3.0*(fhat[13]+fhat[10]+fhat[8]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[29]-1.732050807568877*(fhat[25]+fhat[23]+fhat[21]+fhat[18]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[15]+fhat[14]+fhat[12]+fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0]))); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[24]+fhat[22]))+fhat[20]-1.0*fhat[19]+fhat[17]+fhat[16]-3.0*fhat[2]))+3.0*(fhat[30]-1.0*(fhat[28]+fhat[27]+fhat[26]-3.0*(fhat[13]+fhat[10]+fhat[8]-1.0*fhat[6])))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[15]+fhat[14]-1.0*fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0]))-1.0*(fhat[29]+1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21]+fhat[18])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[1]))))); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[24])+fhat[22]-1.0*fhat[20]+fhat[19]-1.0*fhat[17]+fhat[16]-3.0*fhat[2]))+3.0*((-1.0*fhat[30])+fhat[28]-1.0*(fhat[27]+fhat[26]-3.0*(fhat[13]+fhat[10]-1.0*fhat[8]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[29])+1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]+fhat[18]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[1]))+3.0*(fhat[15]-1.0*fhat[14]+fhat[12]-1.0*fhat[11]+fhat[9]-1.0*fhat[7]+3.0*fhat[0]))); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]-1.0*(fhat[22]+fhat[20]+fhat[19]+fhat[17]-1.0*(fhat[16]+3.0*fhat[2]))))+3.0*(fhat[30]+fhat[28]-1.0*(fhat[27]+fhat[26]-3.0*((-1.0*(fhat[13]+fhat[10]))+fhat[8]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]-1.0*(fhat[21]+fhat[18])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]+fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0]))))); 
  rCtrl[4] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]-1.0*(fhat[22]+fhat[20]-1.0*fhat[19])+fhat[17]-1.0*(fhat[16]+3.0*fhat[2])))+3.0*((-1.0*(fhat[30]+fhat[28]))+fhat[27]-1.0*fhat[26]+3.0*(fhat[13]-1.0*fhat[10]+fhat[8]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[29])+1.732050807568877*(fhat[25]+fhat[23]-1.0*fhat[21]+fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[1])))+3.0*((-1.0*fhat[15])+fhat[14]+fhat[12]-1.0*(fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0]))))); 
  rCtrl[5] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[24])+fhat[22]-1.0*(fhat[20]+fhat[19]-1.0*fhat[17]+fhat[16])+3.0*fhat[2]))+3.0*(fhat[30]-1.0*fhat[28]+fhat[27]-1.0*fhat[26]+3.0*((-1.0*fhat[13])+fhat[10]-1.0*fhat[8]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[29]+1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]-1.0*fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7])+3.0*fhat[0]))); 
  rCtrl[6] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[24]+fhat[22]))+fhat[20]+fhat[19]-1.0*(fhat[17]+fhat[16]-3.0*fhat[2])))+3.0*((-1.0*fhat[30])+fhat[28]+fhat[27]-1.0*fhat[26]+3.0*((-1.0*fhat[13])+fhat[10]+fhat[8]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[29]-1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21])+fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0])))); 
  rCtrl[7] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]+fhat[22]+fhat[20]-1.0*(fhat[19]+fhat[17]+fhat[16]+3.0*fhat[2])))+3.0*(fhat[30]+fhat[28]+fhat[27]-1.0*fhat[26]+3.0*(fhat[13]-1.0*(fhat[10]+fhat[8]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*(fhat[15]+fhat[14]+fhat[12]-1.0*fhat[11]))+fhat[9]+fhat[7]+3.0*fhat[0])-1.0*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]+fhat[21]-1.0*fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[1])))))); 
  rCtrl[8] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]+fhat[22]+fhat[20]-1.0*(fhat[19]+fhat[17]+fhat[16]+3.0*fhat[2])))+3.0*(3.0*((-1.0*fhat[13])+fhat[10]+fhat[8]+fhat[6])-1.0*(fhat[30]+fhat[28]+fhat[27]-1.0*fhat[26]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[29])+1.732050807568877*(fhat[25]+fhat[23]+fhat[21]-1.0*fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[14]+fhat[12]-1.0*fhat[11]))+fhat[9]+fhat[7]+3.0*fhat[0]))); 
  rCtrl[9] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[24]+fhat[22]))+fhat[20]+fhat[19]-1.0*(fhat[17]+fhat[16]-3.0*fhat[2])))+3.0*(fhat[30]-1.0*(fhat[28]+fhat[27]-1.0*fhat[26])+3.0*(fhat[13]-1.0*(fhat[10]+fhat[8]-1.0*fhat[6])))))/(72.0*EPSILON+1.414213562373095*(fhat[29]+1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21])+fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0])))); 
  rCtrl[10] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[24])+fhat[22]-1.0*(fhat[20]+fhat[19]-1.0*fhat[17]+fhat[16])+3.0*fhat[2]))+3.0*((-1.0*fhat[30])+fhat[28]-1.0*fhat[27]+fhat[26]+3.0*(fhat[13]-1.0*fhat[10]+fhat[8]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[29]-1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]-1.0*fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7])+3.0*fhat[0]))); 
  rCtrl[11] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]-1.0*(fhat[22]+fhat[20]-1.0*fhat[19])+fhat[17]-1.0*(fhat[16]+3.0*fhat[2])))+3.0*(fhat[30]+fhat[28]-1.0*fhat[27]+fhat[26]+3.0*((-1.0*fhat[13])+fhat[10]-1.0*(fhat[8]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*fhat[15])+fhat[14]+fhat[12]-1.0*(fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0])))-1.0*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]-1.0*fhat[21]+fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[1])))))); 
  rCtrl[12] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]-1.0*(fhat[22]+fhat[20]+fhat[19]+fhat[17]-1.0*(fhat[16]+3.0*fhat[2]))))+3.0*((-1.0*(fhat[30]+fhat[28]))+fhat[27]+fhat[26]+3.0*(fhat[13]+fhat[10]-1.0*(fhat[8]+fhat[6])))))/(72.0*EPSILON+1.414213562373095*(fhat[29]-1.732050807568877*(fhat[25]+fhat[23]-1.0*(fhat[21]+fhat[18])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]+fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0]))))); 
  rCtrl[13] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[24])+fhat[22]-1.0*fhat[20]+fhat[19]-1.0*fhat[17]+fhat[16]-3.0*fhat[2]))+3.0*(fhat[30]-1.0*fhat[28]+fhat[27]+fhat[26]+3.0*((-1.0*(fhat[13]+fhat[10]))+fhat[8]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[15]-1.0*fhat[14]+fhat[12]-1.0*fhat[11]+fhat[9]-1.0*fhat[7]+3.0*fhat[0])-1.0*(fhat[29]+1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]+fhat[18]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[1]))))); 
  rCtrl[14] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[24]+fhat[22]))+fhat[20]-1.0*fhat[19]+fhat[17]+fhat[16]-3.0*fhat[2]))+3.0*((-1.0*fhat[30])+fhat[28]+fhat[27]+fhat[26]-3.0*(fhat[13]+fhat[10]+fhat[8]-1.0*fhat[6]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[29])+1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21]+fhat[18])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[1]))+3.0*(fhat[15]+fhat[14]-1.0*fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0])))); 
  rCtrl[15] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[24]+fhat[22]+fhat[20]+fhat[19]+fhat[17]+fhat[16]+3.0*fhat[2]))+3.0*(fhat[30]+fhat[28]+fhat[27]+fhat[26]+3.0*(fhat[13]+fhat[10]+fhat[8]+fhat[6]))))/(72.0*EPSILON+1.414213562373095*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]+fhat[21]+fhat[18]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[15]+fhat[14]+fhat[12]+fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0]))); 
  double fhatCtrl[16];  // fhatCtrl = anti-limited fhat evaluated at each control node on y surface 
  fhatCtrl[0] = -0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[23]+fhat[21]+fhat[18]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[1]))-1.0*(fhat[29]+3.0*(fhat[15]+fhat[14]+fhat[12]+fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0])))*limTheta(rCtrl[0],-1.0,EPSILON); 
  fhatCtrl[1] = -0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21]+fhat[18])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[1]))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]-1.0*fhat[11]+fhat[9]+fhat[7]-3.0*fhat[0]))*limTheta(rCtrl[1],-1.0,EPSILON); 
  fhatCtrl[2] = 0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]+fhat[18]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[1]))-1.0*(fhat[29]+3.0*((-1.0*fhat[15])+fhat[14]-1.0*fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7]-3.0*fhat[0])))*limTheta(rCtrl[2],-1.0,EPSILON); 
  fhatCtrl[3] = 0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]-1.0*(fhat[21]+fhat[18])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]+fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0]))))*limTheta(rCtrl[3],-1.0,EPSILON); 
  fhatCtrl[4] = 0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[23]-1.0*fhat[21]+fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[1])))-1.0*(fhat[29]+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]-1.0*fhat[11])+fhat[9]-1.0*(fhat[7]+3.0*fhat[0]))))*limTheta(rCtrl[4],-1.0,EPSILON); 
  fhatCtrl[5] = 0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]-1.0*fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7])+3.0*fhat[0]))*limTheta(rCtrl[5],-1.0,EPSILON); 
  fhatCtrl[6] = -0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21])+fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[1])))-1.0*(fhat[29]+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0]))))*limTheta(rCtrl[6],-1.0,EPSILON); 
  fhatCtrl[7] = -0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]+fhat[21]-1.0*fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[1])))+3.0*(fhat[15]+fhat[14]+fhat[12]-1.0*(fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0])))*limTheta(rCtrl[7],-1.0,EPSILON); 
  fhatCtrl[8] = 0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[23]+fhat[21]-1.0*fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]+fhat[1])))-1.0*(fhat[29]+3.0*(fhat[15]+fhat[14]+fhat[12]-1.0*(fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0]))))*limTheta(rCtrl[8],-1.0,EPSILON); 
  fhatCtrl[9] = 0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21])+fhat[18]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[3]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]+fhat[11]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0])))*limTheta(rCtrl[9],-1.0,EPSILON); 
  fhatCtrl[10] = -0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]-1.0*fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[3]+fhat[1]))-1.0*(fhat[29]+3.0*((-1.0*fhat[15])+fhat[14]-1.0*(fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7])+3.0*fhat[0])))*limTheta(rCtrl[10],-1.0,EPSILON); 
  fhatCtrl[11] = -0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]-1.0*fhat[21]+fhat[18]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[3]+fhat[1])))+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]-1.0*fhat[11])+fhat[9]-1.0*(fhat[7]+3.0*fhat[0])))*limTheta(rCtrl[11],-1.0,EPSILON); 
  fhatCtrl[12] = -0.01964185503295965*(1.732050807568877*(fhat[25]+fhat[23]-1.0*(fhat[21]+fhat[18])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]+fhat[1]))-1.0*(fhat[29]+3.0*(fhat[15]-1.0*(fhat[14]+fhat[12]+fhat[11]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0])))))*limTheta(rCtrl[12],-1.0,EPSILON); 
  fhatCtrl[13] = -0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]-1.0*fhat[23]+fhat[21]+fhat[18]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[3]-1.0*fhat[1]))+3.0*((-1.0*fhat[15])+fhat[14]-1.0*fhat[12]+fhat[11]-1.0*fhat[9]+fhat[7]-3.0*fhat[0]))*limTheta(rCtrl[13],-1.0,EPSILON); 
  fhatCtrl[14] = 0.01964185503295965*(1.732050807568877*(fhat[25]-1.0*(fhat[23]+fhat[21]+fhat[18])+3.0*(fhat[5]+fhat[4]+fhat[3]-1.0*fhat[1]))-1.0*(fhat[29]+3.0*((-1.0*(fhat[15]+fhat[14]))+fhat[12]-1.0*fhat[11]+fhat[9]+fhat[7]-3.0*fhat[0])))*limTheta(rCtrl[14],-1.0,EPSILON); 
  fhatCtrl[15] = 0.01964185503295965*(fhat[29]+1.732050807568877*(fhat[25]+fhat[23]+fhat[21]+fhat[18]+3.0*(fhat[5]+fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[15]+fhat[14]+fhat[12]+fhat[11]+fhat[9]+fhat[7]+3.0*fhat[0]))*limTheta(rCtrl[15],-1.0,EPSILON); 
  double fhatAL[16];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.25*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.4330127018922193*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[5] = 0.75*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[6] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[8] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[9] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[10] = 0.75*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[11] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[12] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[13] = 1.299038105676658*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[14] = 1.299038105676658*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[15] = 2.25*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[16]; 
  fhatALQuad[0] = fmin(0.25*((-0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7])))-1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[0]); 
  fhatALQuad[1] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*((-0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+0.5773502691896258*(fhatAL[13]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[1]); 
  fhatALQuad[2] = fmin(0.25*(0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[2]); 
  fhatALQuad[3] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[3]); 
  fhatALQuad[4] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*(0.5773502691896258*fhatAL[3]-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[4]+fhatAL[2]+fhatAL[1]))-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[4]); 
  fhatALQuad[5] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]))+0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[5]); 
  fhatALQuad[6] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2])-0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[6]); 
  fhatALQuad[7] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[7]); 
  fhatALQuad[8] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]+fhatAL[1]))+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[8]); 
  fhatALQuad[9] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]))+0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[9]); 
  fhatALQuad[10] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[10]); 
  fhatALQuad[11] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[4]+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[11]); 
  fhatALQuad[12] = fmin(0.25*((-0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5])))+1.732050807568877*(0.5773502691896258*(fhatAL[4]+fhatAL[3])-0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[12]); 
  fhatALQuad[13] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*(0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1])-0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[13]); 
  fhatALQuad[14] = fmin(0.25*(0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[1])+0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[14]); 
  fhatALQuad[15] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[15]); 
  fhatAL[0] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[5] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[6] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[8] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[9] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[10] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[11] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[12] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[13] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[14] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[15] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.1767766952966368*(alpha[8]*fhatAL[8]+alpha[6]*fhatAL[6]+alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_y; 
  incr[1] = 0.1767766952966368*(alpha[4]*fhatAL[8]+fhatAL[4]*alpha[8]+alpha[3]*fhatAL[6]+fhatAL[3]*alpha[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_y; 
  incr[2] = -0.3061862178478971*(alpha[8]*fhatAL[8]+alpha[6]*fhatAL[6]+alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_y; 
  incr[3] = 0.1767766952966368*(alpha[8]*fhatAL[12]+alpha[6]*fhatAL[11]+alpha[4]*fhatAL[9]+alpha[3]*fhatAL[7]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_y; 
  incr[4] = 0.1767766952966368*(alpha[8]*fhatAL[13]+alpha[5]*fhatAL[11]+alpha[4]*fhatAL[10]+alpha[2]*fhatAL[7]+alpha[1]*fhatAL[6]+fhatAL[1]*alpha[6]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_y; 
  incr[5] = 0.1767766952966368*(alpha[6]*fhatAL[13]+alpha[5]*fhatAL[12]+alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[1]*fhatAL[8]+fhatAL[1]*alpha[8]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4])*dfac_y; 
  incr[6] = -0.3061862178478971*(alpha[4]*fhatAL[8]+fhatAL[4]*alpha[8]+alpha[3]*fhatAL[6]+fhatAL[3]*alpha[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_y; 
  incr[7] = 0.1767766952966368*(alpha[4]*fhatAL[12]+alpha[3]*fhatAL[11]+alpha[8]*fhatAL[9]+alpha[6]*fhatAL[7]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_y; 
  incr[8] = -0.3061862178478971*(alpha[8]*fhatAL[12]+alpha[6]*fhatAL[11]+alpha[4]*fhatAL[9]+alpha[3]*fhatAL[7]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_y; 
  incr[9] = 0.1767766952966368*(alpha[4]*fhatAL[13]+alpha[2]*fhatAL[11]+alpha[8]*fhatAL[10]+alpha[5]*fhatAL[7]+alpha[0]*fhatAL[6]+fhatAL[0]*alpha[6]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_y; 
  incr[10] = -0.3061862178478971*(alpha[8]*fhatAL[13]+alpha[5]*fhatAL[11]+alpha[4]*fhatAL[10]+alpha[2]*fhatAL[7]+alpha[1]*fhatAL[6]+fhatAL[1]*alpha[6]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_y; 
  incr[11] = 0.1767766952966368*(alpha[8]*fhatAL[15]+alpha[4]*fhatAL[14]+alpha[1]*fhatAL[11]+alpha[0]*fhatAL[7]+alpha[5]*fhatAL[6]+fhatAL[5]*alpha[6]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_y; 
  incr[12] = 0.1767766952966368*(alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[6]*fhatAL[10]+alpha[5]*fhatAL[9]+alpha[0]*fhatAL[8]+fhatAL[0]*alpha[8]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4])*dfac_y; 
  incr[13] = -0.3061862178478971*(alpha[6]*fhatAL[13]+alpha[5]*fhatAL[12]+alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[1]*fhatAL[8]+fhatAL[1]*alpha[8]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4])*dfac_y; 
  incr[14] = 0.1767766952966368*(alpha[6]*fhatAL[15]+alpha[3]*fhatAL[14]+alpha[1]*fhatAL[12]+alpha[0]*fhatAL[9]+alpha[5]*fhatAL[8]+fhatAL[5]*alpha[8]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4])*dfac_y; 
  incr[15] = 0.1767766952966368*(alpha[5]*fhatAL[15]+alpha[2]*fhatAL[14]+alpha[1]*fhatAL[13]+alpha[0]*fhatAL[10]+alpha[6]*fhatAL[8]+fhatAL[6]*alpha[8]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])*dfac_y; 
  incr[16] = -0.3061862178478971*(alpha[4]*fhatAL[12]+alpha[3]*fhatAL[11]+alpha[8]*fhatAL[9]+alpha[6]*fhatAL[7]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_y; 
  incr[17] = -0.3061862178478971*(alpha[4]*fhatAL[13]+alpha[2]*fhatAL[11]+alpha[8]*fhatAL[10]+alpha[5]*fhatAL[7]+alpha[0]*fhatAL[6]+fhatAL[0]*alpha[6]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_y; 
  incr[18] = 0.1767766952966368*(alpha[4]*fhatAL[15]+alpha[8]*fhatAL[14]+alpha[0]*fhatAL[11]+alpha[1]*fhatAL[7]+alpha[2]*fhatAL[6]+fhatAL[2]*alpha[6]+alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5])*dfac_y; 
  incr[19] = -0.3061862178478971*(alpha[8]*fhatAL[15]+alpha[4]*fhatAL[14]+alpha[1]*fhatAL[11]+alpha[0]*fhatAL[7]+alpha[5]*fhatAL[6]+fhatAL[5]*alpha[6]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_y; 
  incr[20] = -0.3061862178478971*(alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[6]*fhatAL[10]+alpha[5]*fhatAL[9]+alpha[0]*fhatAL[8]+fhatAL[0]*alpha[8]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4])*dfac_y; 
  incr[21] = 0.1767766952966368*(alpha[3]*fhatAL[15]+alpha[6]*fhatAL[14]+alpha[0]*fhatAL[12]+alpha[1]*fhatAL[9]+alpha[2]*fhatAL[8]+fhatAL[2]*alpha[8]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5])*dfac_y; 
  incr[22] = -0.3061862178478971*(alpha[6]*fhatAL[15]+alpha[3]*fhatAL[14]+alpha[1]*fhatAL[12]+alpha[0]*fhatAL[9]+alpha[5]*fhatAL[8]+fhatAL[5]*alpha[8]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4])*dfac_y; 
  incr[23] = 0.1767766952966368*(alpha[2]*fhatAL[15]+alpha[5]*fhatAL[14]+alpha[0]*fhatAL[13]+alpha[1]*fhatAL[10]+alpha[3]*fhatAL[8]+fhatAL[3]*alpha[8]+alpha[4]*fhatAL[6]+fhatAL[4]*alpha[6])*dfac_y; 
  incr[24] = -0.3061862178478971*(alpha[5]*fhatAL[15]+alpha[2]*fhatAL[14]+alpha[1]*fhatAL[13]+alpha[0]*fhatAL[10]+alpha[6]*fhatAL[8]+fhatAL[6]*alpha[8]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])*dfac_y; 
  incr[25] = 0.1767766952966368*(alpha[1]*fhatAL[15]+alpha[0]*fhatAL[14]+alpha[5]*fhatAL[13]+alpha[6]*fhatAL[12]+alpha[8]*fhatAL[11]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9]+alpha[4]*fhatAL[7])*dfac_y; 
  incr[26] = -0.3061862178478971*(alpha[4]*fhatAL[15]+alpha[8]*fhatAL[14]+alpha[0]*fhatAL[11]+alpha[1]*fhatAL[7]+alpha[2]*fhatAL[6]+fhatAL[2]*alpha[6]+alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5])*dfac_y; 
  incr[27] = -0.3061862178478971*(alpha[3]*fhatAL[15]+alpha[6]*fhatAL[14]+alpha[0]*fhatAL[12]+alpha[1]*fhatAL[9]+alpha[2]*fhatAL[8]+fhatAL[2]*alpha[8]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5])*dfac_y; 
  incr[28] = -0.3061862178478971*(alpha[2]*fhatAL[15]+alpha[5]*fhatAL[14]+alpha[0]*fhatAL[13]+alpha[1]*fhatAL[10]+alpha[3]*fhatAL[8]+fhatAL[3]*alpha[8]+alpha[4]*fhatAL[6]+fhatAL[4]*alpha[6])*dfac_y; 
  incr[29] = 0.1767766952966368*(alpha[0]*fhatAL[15]+alpha[1]*fhatAL[14]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12]+alpha[4]*fhatAL[11]+alpha[5]*fhatAL[10]+alpha[6]*fhatAL[9]+fhatAL[7]*alpha[8])*dfac_y; 
  incr[30] = -0.3061862178478971*(alpha[1]*fhatAL[15]+alpha[0]*fhatAL[14]+alpha[5]*fhatAL[13]+alpha[6]*fhatAL[12]+alpha[8]*fhatAL[11]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9]+alpha[4]*fhatAL[7])*dfac_y; 
  incr[31] = -0.3061862178478971*(alpha[0]*fhatAL[15]+alpha[1]*fhatAL[14]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12]+alpha[4]*fhatAL[11]+alpha[5]*fhatAL[10]+alpha[6]*fhatAL[9]+fhatAL[7]*alpha[8])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 
  outr[28] += incr[28]; 
  outr[29] += incr[29]; 
  outr[30] += incr[30]; 
  outr[31] += incr[31]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += incr[19]; 
  outl[20] += incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += incr[22]; 
  outl[23] += -1.0*incr[23]; 
  outl[24] += incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += incr[26]; 
  outl[27] += incr[27]; 
  outl[28] += incr[28]; 
  outl[29] += -1.0*incr[29]; 
  outl[30] += incr[30]; 
  outl[31] += incr[31]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity3x2vSer_Z_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_z = 2.0/dxv[2]; 
  double dfac_v = 2.0/dxv[3]; 
  double dfac_m = 2.0/dxv[4]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wz = w[2]; 
  double wv = w[3]; 
  double wm = w[4]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[32]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.1767766952966368*Gradpar[0]*wv; 

  double alpha[16]; 
  alpha[0] = 1.414213562373095*Gradpar[0]*wv; 
  alpha[1] = 1.414213562373095*Gradpar[1]*wv; 
  double f0Quad[16]; 
  double f1Quad[16]; 
  double limQuad[16]; 
  // determine upwinding at each surface quadrature node 
  if(0.25*alpha[0]-0.25*alpha[1] > 0) {
    f0Quad[0] = 0.25*(fl[28]-1.0*(fl[24]+fl[23]+fl[20]+fl[17])+fl[15]+fl[13]+fl[12]+fl[10]+fl[9]+fl[6]-1.0*(fl[5]+fl[4]+fl[2]+fl[1])+fl[0]); 
    f1Quad[0] = -0.25*(fl[31]-1.0*(fl[30]+fl[29]+fl[27]+fl[26])+fl[25]+fl[22]+fl[21]+fl[19]+fl[18]+fl[16]-1.0*(fl[14]+fl[11]+fl[8]+fl[7])+fl[3]); 
    limQuad[0] = fl[0]/cflL; 
  } else {
    f0Quad[0] = 0.25*(fr[28]-1.0*(fr[24]+fr[23]+fr[20]+fr[17])+fr[15]+fr[13]+fr[12]+fr[10]+fr[9]+fr[6]-1.0*(fr[5]+fr[4]+fr[2]+fr[1])+fr[0]); 
    f1Quad[0] = 0.25*(fr[31]-1.0*(fr[30]+fr[29]+fr[27]+fr[26])+fr[25]+fr[22]+fr[21]+fr[19]+fr[18]+fr[16]-1.0*(fr[14]+fr[11]+fr[8]+fr[7])+fr[3]); 
    limQuad[0] = fr[0]/cflR; 
  }
  if(0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[1] = -0.25*(fl[28]+fl[24]-1.0*(fl[23]+fl[20]+fl[17]+fl[15]+fl[13])+fl[12]-1.0*fl[10]+fl[9]+fl[6]+fl[5]+fl[4]+fl[2]-1.0*(fl[1]+fl[0])); 
    f1Quad[1] = 0.25*(fl[31]+fl[30]-1.0*(fl[29]+fl[27]+fl[26]+fl[25]+fl[22])+fl[21]-1.0*fl[19]+fl[18]+fl[16]+fl[14]+fl[11]+fl[8]-1.0*(fl[7]+fl[3])); 
    limQuad[1] = fl[0]/cflL; 
  } else {
    f0Quad[1] = -0.25*(fr[28]+fr[24]-1.0*(fr[23]+fr[20]+fr[17]+fr[15]+fr[13])+fr[12]-1.0*fr[10]+fr[9]+fr[6]+fr[5]+fr[4]+fr[2]-1.0*(fr[1]+fr[0])); 
    f1Quad[1] = -0.25*(fr[31]+fr[30]-1.0*(fr[29]+fr[27]+fr[26]+fr[25]+fr[22])+fr[21]-1.0*fr[19]+fr[18]+fr[16]+fr[14]+fr[11]+fr[8]-1.0*(fr[7]+fr[3])); 
    limQuad[1] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*alpha[1] > 0) {
    f0Quad[2] = -0.25*(fl[28]-1.0*fl[24]+fl[23]-1.0*(fl[20]+fl[17]+fl[15])+fl[13]-1.0*fl[12]+fl[10]-1.0*fl[9]+fl[6]+fl[5]+fl[4]-1.0*fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[2] = 0.25*(fl[31]-1.0*fl[30]+fl[29]-1.0*(fl[27]+fl[26]+fl[25])+fl[22]-1.0*fl[21]+fl[19]-1.0*fl[18]+fl[16]+fl[14]+fl[11]-1.0*fl[8]+fl[7]-1.0*fl[3]); 
    limQuad[2] = fl[0]/cflL; 
  } else {
    f0Quad[2] = -0.25*(fr[28]-1.0*fr[24]+fr[23]-1.0*(fr[20]+fr[17]+fr[15])+fr[13]-1.0*fr[12]+fr[10]-1.0*fr[9]+fr[6]+fr[5]+fr[4]-1.0*fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[2] = -0.25*(fr[31]-1.0*fr[30]+fr[29]-1.0*(fr[27]+fr[26]+fr[25])+fr[22]-1.0*fr[21]+fr[19]-1.0*fr[18]+fr[16]+fr[14]+fr[11]-1.0*fr[8]+fr[7]-1.0*fr[3]); 
    limQuad[2] = fr[0]/cflR; 
  }
  if(0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[3] = 0.25*(fl[28]+fl[24]+fl[23]-1.0*(fl[20]+fl[17])+fl[15]-1.0*(fl[13]+fl[12]+fl[10]+fl[9])+fl[6]-1.0*(fl[5]+fl[4])+fl[2]+fl[1]+fl[0]); 
    f1Quad[3] = -0.25*(fl[31]+fl[30]+fl[29]-1.0*(fl[27]+fl[26])+fl[25]-1.0*(fl[22]+fl[21]+fl[19]+fl[18])+fl[16]-1.0*(fl[14]+fl[11])+fl[8]+fl[7]+fl[3]); 
    limQuad[3] = fl[0]/cflL; 
  } else {
    f0Quad[3] = 0.25*(fr[28]+fr[24]+fr[23]-1.0*(fr[20]+fr[17])+fr[15]-1.0*(fr[13]+fr[12]+fr[10]+fr[9])+fr[6]-1.0*(fr[5]+fr[4])+fr[2]+fr[1]+fr[0]); 
    f1Quad[3] = 0.25*(fr[31]+fr[30]+fr[29]-1.0*(fr[27]+fr[26])+fr[25]-1.0*(fr[22]+fr[21]+fr[19]+fr[18])+fr[16]-1.0*(fr[14]+fr[11])+fr[8]+fr[7]+fr[3]); 
    limQuad[3] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*alpha[1] > 0) {
    f0Quad[4] = -0.25*(fl[28]-1.0*(fl[24]+fl[23])+fl[20]-1.0*fl[17]+fl[15]-1.0*(fl[13]+fl[12])+fl[10]+fl[9]-1.0*fl[6]+fl[5]-1.0*fl[4]+fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[4] = 0.25*(fl[31]-1.0*(fl[30]+fl[29])+fl[27]-1.0*fl[26]+fl[25]-1.0*(fl[22]+fl[21])+fl[19]+fl[18]-1.0*fl[16]+fl[14]-1.0*fl[11]+fl[8]+fl[7]-1.0*fl[3]); 
    limQuad[4] = fl[0]/cflL; 
  } else {
    f0Quad[4] = -0.25*(fr[28]-1.0*(fr[24]+fr[23])+fr[20]-1.0*fr[17]+fr[15]-1.0*(fr[13]+fr[12])+fr[10]+fr[9]-1.0*fr[6]+fr[5]-1.0*fr[4]+fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[4] = -0.25*(fr[31]-1.0*(fr[30]+fr[29])+fr[27]-1.0*fr[26]+fr[25]-1.0*(fr[22]+fr[21])+fr[19]+fr[18]-1.0*fr[16]+fr[14]-1.0*fr[11]+fr[8]+fr[7]-1.0*fr[3]); 
    limQuad[4] = fr[0]/cflR; 
  }
  if(0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[5] = 0.25*(fl[28]+fl[24]-1.0*fl[23]+fl[20]-1.0*(fl[17]+fl[15])+fl[13]-1.0*(fl[12]+fl[10])+fl[9]-1.0*(fl[6]+fl[5])+fl[4]-1.0*fl[2]+fl[1]+fl[0]); 
    f1Quad[5] = -0.25*(fl[31]+fl[30]-1.0*fl[29]+fl[27]-1.0*(fl[26]+fl[25])+fl[22]-1.0*(fl[21]+fl[19])+fl[18]-1.0*(fl[16]+fl[14])+fl[11]-1.0*fl[8]+fl[7]+fl[3]); 
    limQuad[5] = fl[0]/cflL; 
  } else {
    f0Quad[5] = 0.25*(fr[28]+fr[24]-1.0*fr[23]+fr[20]-1.0*(fr[17]+fr[15])+fr[13]-1.0*(fr[12]+fr[10])+fr[9]-1.0*(fr[6]+fr[5])+fr[4]-1.0*fr[2]+fr[1]+fr[0]); 
    f1Quad[5] = 0.25*(fr[31]+fr[30]-1.0*fr[29]+fr[27]-1.0*(fr[26]+fr[25])+fr[22]-1.0*(fr[21]+fr[19])+fr[18]-1.0*(fr[16]+fr[14])+fr[11]-1.0*fr[8]+fr[7]+fr[3]); 
    limQuad[5] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*alpha[1] > 0) {
    f0Quad[6] = 0.25*(fl[28]-1.0*fl[24]+fl[23]+fl[20]-1.0*(fl[17]+fl[15]+fl[13])+fl[12]+fl[10]-1.0*(fl[9]+fl[6]+fl[5])+fl[4]+fl[2]-1.0*fl[1]+fl[0]); 
    f1Quad[6] = -0.25*(fl[31]-1.0*fl[30]+fl[29]+fl[27]-1.0*(fl[26]+fl[25]+fl[22])+fl[21]+fl[19]-1.0*(fl[18]+fl[16]+fl[14])+fl[11]+fl[8]-1.0*fl[7]+fl[3]); 
    limQuad[6] = fl[0]/cflL; 
  } else {
    f0Quad[6] = 0.25*(fr[28]-1.0*fr[24]+fr[23]+fr[20]-1.0*(fr[17]+fr[15]+fr[13])+fr[12]+fr[10]-1.0*(fr[9]+fr[6]+fr[5])+fr[4]+fr[2]-1.0*fr[1]+fr[0]); 
    f1Quad[6] = 0.25*(fr[31]-1.0*fr[30]+fr[29]+fr[27]-1.0*(fr[26]+fr[25]+fr[22])+fr[21]+fr[19]-1.0*(fr[18]+fr[16]+fr[14])+fr[11]+fr[8]-1.0*fr[7]+fr[3]); 
    limQuad[6] = fr[0]/cflR; 
  }
  if(0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[7] = -0.25*(fl[28]+fl[24]+fl[23]+fl[20]-1.0*fl[17]+fl[15]+fl[13]+fl[12]-1.0*(fl[10]+fl[9]+fl[6])+fl[5]-1.0*(fl[4]+fl[2]+fl[1]+fl[0])); 
    f1Quad[7] = 0.25*(fl[31]+fl[30]+fl[29]+fl[27]-1.0*fl[26]+fl[25]+fl[22]+fl[21]-1.0*(fl[19]+fl[18]+fl[16])+fl[14]-1.0*(fl[11]+fl[8]+fl[7]+fl[3])); 
    limQuad[7] = fl[0]/cflL; 
  } else {
    f0Quad[7] = -0.25*(fr[28]+fr[24]+fr[23]+fr[20]-1.0*fr[17]+fr[15]+fr[13]+fr[12]-1.0*(fr[10]+fr[9]+fr[6])+fr[5]-1.0*(fr[4]+fr[2]+fr[1]+fr[0])); 
    f1Quad[7] = -0.25*(fr[31]+fr[30]+fr[29]+fr[27]-1.0*fr[26]+fr[25]+fr[22]+fr[21]-1.0*(fr[19]+fr[18]+fr[16])+fr[14]-1.0*(fr[11]+fr[8]+fr[7]+fr[3])); 
    limQuad[7] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*alpha[1] > 0) {
    f0Quad[8] = -0.25*(fl[28]-1.0*(fl[24]+fl[23]+fl[20])+fl[17]+fl[15]+fl[13]+fl[12]-1.0*(fl[10]+fl[9]+fl[6]+fl[5])+fl[4]+fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[8] = 0.25*(fl[31]-1.0*(fl[30]+fl[29]+fl[27])+fl[26]+fl[25]+fl[22]+fl[21]-1.0*(fl[19]+fl[18]+fl[16]+fl[14])+fl[11]+fl[8]+fl[7]-1.0*fl[3]); 
    limQuad[8] = fl[0]/cflL; 
  } else {
    f0Quad[8] = -0.25*(fr[28]-1.0*(fr[24]+fr[23]+fr[20])+fr[17]+fr[15]+fr[13]+fr[12]-1.0*(fr[10]+fr[9]+fr[6]+fr[5])+fr[4]+fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[8] = -0.25*(fr[31]-1.0*(fr[30]+fr[29]+fr[27])+fr[26]+fr[25]+fr[22]+fr[21]-1.0*(fr[19]+fr[18]+fr[16]+fr[14])+fr[11]+fr[8]+fr[7]-1.0*fr[3]); 
    limQuad[8] = fr[0]/cflR; 
  }
  if(0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[9] = 0.25*(fl[28]+fl[24]-1.0*(fl[23]+fl[20])+fl[17]-1.0*(fl[15]+fl[13])+fl[12]+fl[10]-1.0*(fl[9]+fl[6])+fl[5]-1.0*(fl[4]+fl[2])+fl[1]+fl[0]); 
    f1Quad[9] = -0.25*(fl[31]+fl[30]-1.0*(fl[29]+fl[27])+fl[26]-1.0*(fl[25]+fl[22])+fl[21]+fl[19]-1.0*(fl[18]+fl[16])+fl[14]-1.0*(fl[11]+fl[8])+fl[7]+fl[3]); 
    limQuad[9] = fl[0]/cflL; 
  } else {
    f0Quad[9] = 0.25*(fr[28]+fr[24]-1.0*(fr[23]+fr[20])+fr[17]-1.0*(fr[15]+fr[13])+fr[12]+fr[10]-1.0*(fr[9]+fr[6])+fr[5]-1.0*(fr[4]+fr[2])+fr[1]+fr[0]); 
    f1Quad[9] = 0.25*(fr[31]+fr[30]-1.0*(fr[29]+fr[27])+fr[26]-1.0*(fr[25]+fr[22])+fr[21]+fr[19]-1.0*(fr[18]+fr[16])+fr[14]-1.0*(fr[11]+fr[8])+fr[7]+fr[3]); 
    limQuad[9] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*alpha[1] > 0) {
    f0Quad[10] = 0.25*(fl[28]-1.0*fl[24]+fl[23]-1.0*fl[20]+fl[17]-1.0*fl[15]+fl[13]-1.0*(fl[12]+fl[10])+fl[9]-1.0*fl[6]+fl[5]-1.0*fl[4]+fl[2]-1.0*fl[1]+fl[0]); 
    f1Quad[10] = -0.25*(fl[31]-1.0*fl[30]+fl[29]-1.0*fl[27]+fl[26]-1.0*fl[25]+fl[22]-1.0*(fl[21]+fl[19])+fl[18]-1.0*fl[16]+fl[14]-1.0*fl[11]+fl[8]-1.0*fl[7]+fl[3]); 
    limQuad[10] = fl[0]/cflL; 
  } else {
    f0Quad[10] = 0.25*(fr[28]-1.0*fr[24]+fr[23]-1.0*fr[20]+fr[17]-1.0*fr[15]+fr[13]-1.0*(fr[12]+fr[10])+fr[9]-1.0*fr[6]+fr[5]-1.0*fr[4]+fr[2]-1.0*fr[1]+fr[0]); 
    f1Quad[10] = 0.25*(fr[31]-1.0*fr[30]+fr[29]-1.0*fr[27]+fr[26]-1.0*fr[25]+fr[22]-1.0*(fr[21]+fr[19])+fr[18]-1.0*fr[16]+fr[14]-1.0*fr[11]+fr[8]-1.0*fr[7]+fr[3]); 
    limQuad[10] = fr[0]/cflR; 
  }
  if(0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[11] = -0.25*(fl[28]+fl[24]+fl[23]-1.0*fl[20]+fl[17]+fl[15]-1.0*(fl[13]+fl[12])+fl[10]+fl[9]-1.0*(fl[6]+fl[5])+fl[4]-1.0*(fl[2]+fl[1]+fl[0])); 
    f1Quad[11] = 0.25*(fl[31]+fl[30]+fl[29]-1.0*fl[27]+fl[26]+fl[25]-1.0*(fl[22]+fl[21])+fl[19]+fl[18]-1.0*(fl[16]+fl[14])+fl[11]-1.0*(fl[8]+fl[7]+fl[3])); 
    limQuad[11] = fl[0]/cflL; 
  } else {
    f0Quad[11] = -0.25*(fr[28]+fr[24]+fr[23]-1.0*fr[20]+fr[17]+fr[15]-1.0*(fr[13]+fr[12])+fr[10]+fr[9]-1.0*(fr[6]+fr[5])+fr[4]-1.0*(fr[2]+fr[1]+fr[0])); 
    f1Quad[11] = -0.25*(fr[31]+fr[30]+fr[29]-1.0*fr[27]+fr[26]+fr[25]-1.0*(fr[22]+fr[21])+fr[19]+fr[18]-1.0*(fr[16]+fr[14])+fr[11]-1.0*(fr[8]+fr[7]+fr[3])); 
    limQuad[11] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*alpha[1] > 0) {
    f0Quad[12] = 0.25*(fl[28]-1.0*(fl[24]+fl[23])+fl[20]+fl[17]+fl[15]-1.0*(fl[13]+fl[12]+fl[10]+fl[9])+fl[6]+fl[5]+fl[4]-1.0*(fl[2]+fl[1])+fl[0]); 
    f1Quad[12] = -0.25*(fl[31]-1.0*(fl[30]+fl[29])+fl[27]+fl[26]+fl[25]-1.0*(fl[22]+fl[21]+fl[19]+fl[18])+fl[16]+fl[14]+fl[11]-1.0*(fl[8]+fl[7])+fl[3]); 
    limQuad[12] = fl[0]/cflL; 
  } else {
    f0Quad[12] = 0.25*(fr[28]-1.0*(fr[24]+fr[23])+fr[20]+fr[17]+fr[15]-1.0*(fr[13]+fr[12]+fr[10]+fr[9])+fr[6]+fr[5]+fr[4]-1.0*(fr[2]+fr[1])+fr[0]); 
    f1Quad[12] = 0.25*(fr[31]-1.0*(fr[30]+fr[29])+fr[27]+fr[26]+fr[25]-1.0*(fr[22]+fr[21]+fr[19]+fr[18])+fr[16]+fr[14]+fr[11]-1.0*(fr[8]+fr[7])+fr[3]); 
    limQuad[12] = fr[0]/cflR; 
  }
  if(0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[13] = -0.25*(fl[28]+fl[24]-1.0*fl[23]+fl[20]+fl[17]-1.0*fl[15]+fl[13]-1.0*fl[12]+fl[10]-1.0*fl[9]+fl[6]-1.0*(fl[5]+fl[4])+fl[2]-1.0*(fl[1]+fl[0])); 
    f1Quad[13] = 0.25*(fl[31]+fl[30]-1.0*fl[29]+fl[27]+fl[26]-1.0*fl[25]+fl[22]-1.0*fl[21]+fl[19]-1.0*fl[18]+fl[16]-1.0*(fl[14]+fl[11])+fl[8]-1.0*(fl[7]+fl[3])); 
    limQuad[13] = fl[0]/cflL; 
  } else {
    f0Quad[13] = -0.25*(fr[28]+fr[24]-1.0*fr[23]+fr[20]+fr[17]-1.0*fr[15]+fr[13]-1.0*fr[12]+fr[10]-1.0*fr[9]+fr[6]-1.0*(fr[5]+fr[4])+fr[2]-1.0*(fr[1]+fr[0])); 
    f1Quad[13] = -0.25*(fr[31]+fr[30]-1.0*fr[29]+fr[27]+fr[26]-1.0*fr[25]+fr[22]-1.0*fr[21]+fr[19]-1.0*fr[18]+fr[16]-1.0*(fr[14]+fr[11])+fr[8]-1.0*(fr[7]+fr[3])); 
    limQuad[13] = fr[0]/cflR; 
  }
  if(0.25*alpha[0]-0.25*alpha[1] > 0) {
    f0Quad[14] = -0.25*(fl[28]-1.0*fl[24]+fl[23]+fl[20]+fl[17]-1.0*(fl[15]+fl[13])+fl[12]-1.0*fl[10]+fl[9]+fl[6]-1.0*(fl[5]+fl[4]+fl[2])+fl[1]-1.0*fl[0]); 
    f1Quad[14] = 0.25*(fl[31]-1.0*fl[30]+fl[29]+fl[27]+fl[26]-1.0*(fl[25]+fl[22])+fl[21]-1.0*fl[19]+fl[18]+fl[16]-1.0*(fl[14]+fl[11]+fl[8])+fl[7]-1.0*fl[3]); 
    limQuad[14] = fl[0]/cflL; 
  } else {
    f0Quad[14] = -0.25*(fr[28]-1.0*fr[24]+fr[23]+fr[20]+fr[17]-1.0*(fr[15]+fr[13])+fr[12]-1.0*fr[10]+fr[9]+fr[6]-1.0*(fr[5]+fr[4]+fr[2])+fr[1]-1.0*fr[0]); 
    f1Quad[14] = -0.25*(fr[31]-1.0*fr[30]+fr[29]+fr[27]+fr[26]-1.0*(fr[25]+fr[22])+fr[21]-1.0*fr[19]+fr[18]+fr[16]-1.0*(fr[14]+fr[11]+fr[8])+fr[7]-1.0*fr[3]); 
    limQuad[14] = fr[0]/cflR; 
  }
  if(0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[15] = 0.25*(fl[28]+fl[24]+fl[23]+fl[20]+fl[17]+fl[15]+fl[13]+fl[12]+fl[10]+fl[9]+fl[6]+fl[5]+fl[4]+fl[2]+fl[1]+fl[0]); 
    f1Quad[15] = -0.25*(fl[31]+fl[30]+fl[29]+fl[27]+fl[26]+fl[25]+fl[22]+fl[21]+fl[19]+fl[18]+fl[16]+fl[14]+fl[11]+fl[8]+fl[7]+fl[3]); 
    limQuad[15] = fl[0]/cflL; 
  } else {
    f0Quad[15] = 0.25*(fr[28]+fr[24]+fr[23]+fr[20]+fr[17]+fr[15]+fr[13]+fr[12]+fr[10]+fr[9]+fr[6]+fr[5]+fr[4]+fr[2]+fr[1]+fr[0]); 
    f1Quad[15] = 0.25*(fr[31]+fr[30]+fr[29]+fr[27]+fr[26]+fr[25]+fr[22]+fr[21]+fr[19]+fr[18]+fr[16]+fr[14]+fr[11]+fr[8]+fr[7]+fr[3]); 
    limQuad[15] = fr[0]/cflR; 
  }
  double fhat[32]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[3] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[4] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[6] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[7] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[8] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[9] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[10] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[11] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[12] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[13] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[14] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[15] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[16] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[17] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[18] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[19] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[20] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[21] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[22] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[23] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[24] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[25] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[26] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[27] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[28] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[29] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[30] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[31] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[16];  // rCtrl=f1/f0 at each control node in dimensions other than z 
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[22]+fhat[21]+fhat[19]+fhat[18]+fhat[16]+3.0*fhat[3]))-3.0*(fhat[30]+fhat[29]+fhat[27]+fhat[26]+3.0*(fhat[14]+fhat[11]+fhat[8]+fhat[7]))))/(72.0*EPSILON+1.414213562373095*(fhat[28]-1.732050807568877*(fhat[24]+fhat[23]+fhat[20]+fhat[17]+3.0*(fhat[5]+fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[15]+fhat[13]+fhat[12]+fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0]))); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[22]))+fhat[21]-1.0*fhat[19]+fhat[18]+fhat[16]-3.0*fhat[3]))+3.0*(fhat[30]-1.0*(fhat[29]+fhat[27]+fhat[26]-3.0*(fhat[14]+fhat[11]+fhat[8]-1.0*fhat[7])))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[15]+fhat[13]-1.0*fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0]))-1.0*(fhat[28]+1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20]+fhat[17])+3.0*(fhat[5]+fhat[4]+fhat[2]-1.0*fhat[1]))))); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[22]-1.0*fhat[21]+fhat[19]-1.0*fhat[18]+fhat[16]-3.0*fhat[3]))+3.0*((-1.0*fhat[30])+fhat[29]-1.0*(fhat[27]+fhat[26]-3.0*(fhat[14]+fhat[11]-1.0*fhat[8]+fhat[7])))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[28])+1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]+fhat[17]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]-1.0*fhat[1]))+3.0*(fhat[15]-1.0*fhat[13]+fhat[12]-1.0*fhat[10]+fhat[9]-1.0*fhat[6]+3.0*fhat[0]))); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[22]+fhat[21]+fhat[19]+fhat[18]-1.0*(fhat[16]+3.0*fhat[3]))))+3.0*(fhat[30]+fhat[29]-1.0*(fhat[27]+fhat[26]-3.0*((-1.0*(fhat[14]+fhat[11]))+fhat[8]+fhat[7])))))/(72.0*EPSILON+1.414213562373095*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]-1.0*(fhat[20]+fhat[17])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]+fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[4] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[22]+fhat[21]-1.0*fhat[19])+fhat[18]-1.0*(fhat[16]+3.0*fhat[3])))+3.0*((-1.0*(fhat[30]+fhat[29]))+fhat[27]-1.0*fhat[26]+3.0*(fhat[14]-1.0*fhat[11]+fhat[8]+fhat[7]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[28])+1.732050807568877*(fhat[24]+fhat[23]-1.0*fhat[20]+fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[2]+fhat[1])))+3.0*((-1.0*fhat[15])+fhat[13]+fhat[12]-1.0*(fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[5] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[22]-1.0*(fhat[21]+fhat[19]-1.0*fhat[18]+fhat[16])+3.0*fhat[3]))+3.0*(fhat[30]-1.0*fhat[29]+fhat[27]-1.0*fhat[26]+3.0*((-1.0*fhat[14])+fhat[11]-1.0*fhat[8]+fhat[7]))))/(72.0*EPSILON+1.414213562373095*(fhat[28]+1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]-1.0*fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[13]-1.0*(fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6])+3.0*fhat[0]))); 
  rCtrl[6] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[22]))+fhat[21]+fhat[19]-1.0*(fhat[18]+fhat[16]-3.0*fhat[3])))+3.0*((-1.0*fhat[30])+fhat[29]+fhat[27]-1.0*fhat[26]+3.0*((-1.0*fhat[14])+fhat[11]+fhat[8]-1.0*fhat[7]))))/(72.0*EPSILON+1.414213562373095*(fhat[28]-1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20])+fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[7] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[22]+fhat[21]-1.0*(fhat[19]+fhat[18]+fhat[16]+3.0*fhat[3])))+3.0*(fhat[30]+fhat[29]+fhat[27]-1.0*fhat[26]+3.0*(fhat[14]-1.0*(fhat[11]+fhat[8]+fhat[7])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*(fhat[15]+fhat[13]+fhat[12]-1.0*fhat[10]))+fhat[9]+fhat[6]+3.0*fhat[0])-1.0*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]+fhat[20]-1.0*fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]+fhat[1])))))); 
  rCtrl[8] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[22]+fhat[21]-1.0*(fhat[19]+fhat[18]+fhat[16]+3.0*fhat[3])))+3.0*(3.0*((-1.0*fhat[14])+fhat[11]+fhat[8]+fhat[7])-1.0*(fhat[30]+fhat[29]+fhat[27]-1.0*fhat[26]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[28])+1.732050807568877*(fhat[24]+fhat[23]+fhat[20]-1.0*fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]+fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[13]+fhat[12]-1.0*fhat[10]))+fhat[9]+fhat[6]+3.0*fhat[0]))); 
  rCtrl[9] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[22]))+fhat[21]+fhat[19]-1.0*(fhat[18]+fhat[16]-3.0*fhat[3])))+3.0*(fhat[30]-1.0*(fhat[29]+fhat[27]-1.0*fhat[26])+3.0*(fhat[14]-1.0*(fhat[11]+fhat[8]-1.0*fhat[7])))))/(72.0*EPSILON+1.414213562373095*(fhat[28]+1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20])+fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[10] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[22]-1.0*(fhat[21]+fhat[19]-1.0*fhat[18]+fhat[16])+3.0*fhat[3]))+3.0*((-1.0*fhat[30])+fhat[29]-1.0*fhat[27]+fhat[26]+3.0*(fhat[14]-1.0*fhat[11]+fhat[8]-1.0*fhat[7]))))/(72.0*EPSILON+1.414213562373095*(fhat[28]-1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]-1.0*fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[13]-1.0*(fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6])+3.0*fhat[0]))); 
  rCtrl[11] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[22]+fhat[21]-1.0*fhat[19])+fhat[18]-1.0*(fhat[16]+3.0*fhat[3])))+3.0*(fhat[30]+fhat[29]-1.0*fhat[27]+fhat[26]+3.0*((-1.0*fhat[14])+fhat[11]-1.0*(fhat[8]+fhat[7])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*fhat[15])+fhat[13]+fhat[12]-1.0*(fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0])))-1.0*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]-1.0*fhat[20]+fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[2]+fhat[1])))))); 
  rCtrl[12] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[22]+fhat[21]+fhat[19]+fhat[18]-1.0*(fhat[16]+3.0*fhat[3]))))+3.0*((-1.0*(fhat[30]+fhat[29]))+fhat[27]+fhat[26]+3.0*(fhat[14]+fhat[11]-1.0*(fhat[8]+fhat[7])))))/(72.0*EPSILON+1.414213562373095*(fhat[28]-1.732050807568877*(fhat[24]+fhat[23]-1.0*(fhat[20]+fhat[17])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]+fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[13] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[22]-1.0*fhat[21]+fhat[19]-1.0*fhat[18]+fhat[16]-3.0*fhat[3]))+3.0*(fhat[30]-1.0*fhat[29]+fhat[27]+fhat[26]+3.0*((-1.0*(fhat[14]+fhat[11]))+fhat[8]-1.0*fhat[7]))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[15]-1.0*fhat[13]+fhat[12]-1.0*fhat[10]+fhat[9]-1.0*fhat[6]+3.0*fhat[0])-1.0*(fhat[28]+1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]+fhat[17]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]-1.0*fhat[1]))))); 
  rCtrl[14] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[22]))+fhat[21]-1.0*fhat[19]+fhat[18]+fhat[16]-3.0*fhat[3]))+3.0*((-1.0*fhat[30])+fhat[29]+fhat[27]+fhat[26]-3.0*(fhat[14]+fhat[11]+fhat[8]-1.0*fhat[7]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[28])+1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20]+fhat[17])+3.0*(fhat[5]+fhat[4]+fhat[2]-1.0*fhat[1]))+3.0*(fhat[15]+fhat[13]-1.0*fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[15] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[22]+fhat[21]+fhat[19]+fhat[18]+fhat[16]+3.0*fhat[3]))+3.0*(fhat[30]+fhat[29]+fhat[27]+fhat[26]+3.0*(fhat[14]+fhat[11]+fhat[8]+fhat[7]))))/(72.0*EPSILON+1.414213562373095*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]+fhat[20]+fhat[17]+3.0*(fhat[5]+fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[15]+fhat[13]+fhat[12]+fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0]))); 
  double fhatCtrl[16];  // fhatCtrl = anti-limited fhat evaluated at each control node on z surface 
  fhatCtrl[0] = -0.01964185503295965*(1.732050807568877*(fhat[24]+fhat[23]+fhat[20]+fhat[17]+3.0*(fhat[5]+fhat[4]+fhat[2]+fhat[1]))-1.0*(fhat[28]+3.0*(fhat[15]+fhat[13]+fhat[12]+fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[0],-1.0,EPSILON); 
  fhatCtrl[1] = -0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20]+fhat[17])+3.0*(fhat[5]+fhat[4]+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]-1.0*fhat[10]+fhat[9]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[1],-1.0,EPSILON); 
  fhatCtrl[2] = 0.01964185503295965*(1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]+fhat[17]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]-1.0*fhat[1]))-1.0*(fhat[28]+3.0*((-1.0*fhat[15])+fhat[13]-1.0*fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[2],-1.0,EPSILON); 
  fhatCtrl[3] = 0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]-1.0*(fhat[20]+fhat[17])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]+fhat[1]))+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]+fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[3],-1.0,EPSILON); 
  fhatCtrl[4] = 0.01964185503295965*(1.732050807568877*(fhat[24]+fhat[23]-1.0*fhat[20]+fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[2]+fhat[1])))-1.0*(fhat[28]+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]-1.0*fhat[10])+fhat[9]-1.0*(fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[4],-1.0,EPSILON); 
  fhatCtrl[5] = 0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]-1.0*fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[15])+fhat[13]-1.0*(fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6])+3.0*fhat[0]))*limTheta(rCtrl[5],-1.0,EPSILON); 
  fhatCtrl[6] = -0.01964185503295965*(1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20])+fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]-1.0*fhat[1])))-1.0*(fhat[28]+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0]))))*limTheta(rCtrl[6],-1.0,EPSILON); 
  fhatCtrl[7] = -0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]+fhat[20]-1.0*fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]+fhat[1])))+3.0*(fhat[15]+fhat[13]+fhat[12]-1.0*(fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[7],-1.0,EPSILON); 
  fhatCtrl[8] = 0.01964185503295965*(1.732050807568877*(fhat[24]+fhat[23]+fhat[20]-1.0*fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]+fhat[1])))-1.0*(fhat[28]+3.0*(fhat[15]+fhat[13]+fhat[12]-1.0*(fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[8],-1.0,EPSILON); 
  fhatCtrl[9] = 0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20])+fhat[17]+3.0*(fhat[5]-1.0*(fhat[4]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]+fhat[10]-1.0*(fhat[9]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[9],-1.0,EPSILON); 
  fhatCtrl[10] = -0.01964185503295965*(1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]-1.0*fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*fhat[2]+fhat[1]))-1.0*(fhat[28]+3.0*((-1.0*fhat[15])+fhat[13]-1.0*(fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6])+3.0*fhat[0])))*limTheta(rCtrl[10],-1.0,EPSILON); 
  fhatCtrl[11] = -0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]-1.0*fhat[20]+fhat[17]+3.0*((-1.0*fhat[5])+fhat[4]-1.0*(fhat[2]+fhat[1])))+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]-1.0*fhat[10])+fhat[9]-1.0*(fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[11],-1.0,EPSILON); 
  fhatCtrl[12] = -0.01964185503295965*(1.732050807568877*(fhat[24]+fhat[23]-1.0*(fhat[20]+fhat[17])+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]+fhat[1]))-1.0*(fhat[28]+3.0*(fhat[15]-1.0*(fhat[13]+fhat[12]+fhat[10]+fhat[9]-1.0*(fhat[6]+3.0*fhat[0])))))*limTheta(rCtrl[12],-1.0,EPSILON); 
  fhatCtrl[13] = -0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]-1.0*fhat[23]+fhat[20]+fhat[17]+3.0*((-1.0*(fhat[5]+fhat[4]))+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*fhat[15])+fhat[13]-1.0*fhat[12]+fhat[10]-1.0*fhat[9]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[13],-1.0,EPSILON); 
  fhatCtrl[14] = 0.01964185503295965*(1.732050807568877*(fhat[24]-1.0*(fhat[23]+fhat[20]+fhat[17])+3.0*(fhat[5]+fhat[4]+fhat[2]-1.0*fhat[1]))-1.0*(fhat[28]+3.0*((-1.0*(fhat[15]+fhat[13]))+fhat[12]-1.0*fhat[10]+fhat[9]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[14],-1.0,EPSILON); 
  fhatCtrl[15] = 0.01964185503295965*(fhat[28]+1.732050807568877*(fhat[24]+fhat[23]+fhat[20]+fhat[17]+3.0*(fhat[5]+fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[15]+fhat[13]+fhat[12]+fhat[10]+fhat[9]+fhat[6]+3.0*fhat[0]))*limTheta(rCtrl[15],-1.0,EPSILON); 
  double fhatAL[16];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.25*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.4330127018922193*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[5] = 0.75*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[6] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[8] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[9] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[10] = 0.75*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[11] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[12] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[13] = 1.299038105676658*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[14] = 1.299038105676658*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[15] = 2.25*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[16]; 
  fhatALQuad[0] = fmin(0.25*((-0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7])))-1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[0]); 
  fhatALQuad[1] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*((-0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+0.5773502691896258*(fhatAL[13]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[1]); 
  fhatALQuad[2] = fmin(0.25*(0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[5]))-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[2]); 
  fhatALQuad[3] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])-1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1])-0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[3]); 
  fhatALQuad[4] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*(0.5773502691896258*fhatAL[3]-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[4]+fhatAL[2]+fhatAL[1]))-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[4]); 
  fhatALQuad[5] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]))+0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[5]); 
  fhatALQuad[6] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2])-0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[6]); 
  fhatALQuad[7] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])-1.732050807568877*fhatAL[9]+1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[12])+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[8])+fhatAL[0]), limQuad[7]); 
  fhatALQuad[8] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]+fhatAL[1]))+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[8]); 
  fhatALQuad[9] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7]))+1.732050807568877*((-0.5773502691896258*((-1.0*fhatAL[14])+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2]))+0.5773502691896258*((-1.0*fhatAL[13])+fhatAL[4]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[9]); 
  fhatALQuad[10] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*fhatAL[14])-0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[1]-1.0*fhatAL[13])+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]-0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[10]); 
  fhatALQuad[11] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[5]-3.0*fhatAL[15])+1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[7])+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[13]))+0.5773502691896258*(1.732050807568877*fhatAL[12]-1.732050807568877*fhatAL[11])+fhatAL[4]+fhatAL[2]+fhatAL[1])-0.5773502691896258*fhatAL[3])-1.0*fhatAL[10]+0.5773502691896258*(1.732050807568877*fhatAL[8]-1.732050807568877*fhatAL[6])+fhatAL[0]), limQuad[11]); 
  fhatALQuad[12] = fmin(0.25*((-0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5])))+1.732050807568877*(0.5773502691896258*(fhatAL[4]+fhatAL[3])-0.5773502691896258*(fhatAL[14]+fhatAL[13]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]+fhatAL[1]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[12]); 
  fhatALQuad[13] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7])))+1.732050807568877*(0.5773502691896258*(fhatAL[13]+fhatAL[4]+fhatAL[3]+fhatAL[1])-0.5773502691896258*(fhatAL[14]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2]))+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[13]); 
  fhatALQuad[14] = fmin(0.25*(0.5773502691896258*(1.732050807568877*(fhatAL[9]+fhatAL[7])-1.732050807568878*(fhatAL[15]+fhatAL[5]))+1.732050807568877*(0.5773502691896258*(fhatAL[14]-1.0*(fhatAL[12]+fhatAL[11])+fhatAL[2])-0.5773502691896258*(fhatAL[13]+fhatAL[1])+0.5773502691896258*(fhatAL[4]+fhatAL[3]))+fhatAL[10]-1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[14]); 
  fhatALQuad[15] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[5])+1.732050807568877*(fhatAL[9]+fhatAL[7]))+1.0*(fhatAL[14]+fhatAL[13]+1.0*(fhatAL[12]+fhatAL[11])+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[10]+1.0*(fhatAL[8]+fhatAL[6])+fhatAL[0]), limQuad[15]); 
  fhatAL[0] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[5] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[6] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[8] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[9] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[10] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[11] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[12] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[13] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[14] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[15] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.1767766952966368*(alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_z; 
  incr[1] = 0.1767766952966368*(alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_z; 
  incr[2] = 0.1767766952966368*(alpha[1]*fhatAL[5]+alpha[0]*fhatAL[2])*dfac_z; 
  incr[3] = -0.3061862178478971*(alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_z; 
  incr[4] = 0.1767766952966368*(alpha[1]*fhatAL[6]+alpha[0]*fhatAL[3])*dfac_z; 
  incr[5] = 0.1767766952966368*(alpha[1]*fhatAL[8]+alpha[0]*fhatAL[4])*dfac_z; 
  incr[6] = 0.1767766952966368*(alpha[0]*fhatAL[5]+alpha[1]*fhatAL[2])*dfac_z; 
  incr[7] = -0.3061862178478971*(alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_z; 
  incr[8] = -0.3061862178478971*(alpha[1]*fhatAL[5]+alpha[0]*fhatAL[2])*dfac_z; 
  incr[9] = 0.1767766952966368*(alpha[0]*fhatAL[6]+alpha[1]*fhatAL[3])*dfac_z; 
  incr[10] = 0.1767766952966368*(alpha[1]*fhatAL[11]+alpha[0]*fhatAL[7])*dfac_z; 
  incr[11] = -0.3061862178478971*(alpha[1]*fhatAL[6]+alpha[0]*fhatAL[3])*dfac_z; 
  incr[12] = 0.1767766952966368*(alpha[0]*fhatAL[8]+alpha[1]*fhatAL[4])*dfac_z; 
  incr[13] = 0.1767766952966368*(alpha[1]*fhatAL[12]+alpha[0]*fhatAL[9])*dfac_z; 
  incr[14] = -0.3061862178478971*(alpha[1]*fhatAL[8]+alpha[0]*fhatAL[4])*dfac_z; 
  incr[15] = 0.1767766952966368*(alpha[1]*fhatAL[13]+alpha[0]*fhatAL[10])*dfac_z; 
  incr[16] = -0.3061862178478971*(alpha[0]*fhatAL[5]+alpha[1]*fhatAL[2])*dfac_z; 
  incr[17] = 0.1767766952966368*(alpha[0]*fhatAL[11]+alpha[1]*fhatAL[7])*dfac_z; 
  incr[18] = -0.3061862178478971*(alpha[0]*fhatAL[6]+alpha[1]*fhatAL[3])*dfac_z; 
  incr[19] = -0.3061862178478971*(alpha[1]*fhatAL[11]+alpha[0]*fhatAL[7])*dfac_z; 
  incr[20] = 0.1767766952966368*(alpha[0]*fhatAL[12]+alpha[1]*fhatAL[9])*dfac_z; 
  incr[21] = -0.3061862178478971*(alpha[0]*fhatAL[8]+alpha[1]*fhatAL[4])*dfac_z; 
  incr[22] = -0.3061862178478971*(alpha[1]*fhatAL[12]+alpha[0]*fhatAL[9])*dfac_z; 
  incr[23] = 0.1767766952966368*(alpha[0]*fhatAL[13]+alpha[1]*fhatAL[10])*dfac_z; 
  incr[24] = 0.1767766952966368*(alpha[1]*fhatAL[15]+alpha[0]*fhatAL[14])*dfac_z; 
  incr[25] = -0.3061862178478971*(alpha[1]*fhatAL[13]+alpha[0]*fhatAL[10])*dfac_z; 
  incr[26] = -0.3061862178478971*(alpha[0]*fhatAL[11]+alpha[1]*fhatAL[7])*dfac_z; 
  incr[27] = -0.3061862178478971*(alpha[0]*fhatAL[12]+alpha[1]*fhatAL[9])*dfac_z; 
  incr[28] = 0.1767766952966368*(alpha[0]*fhatAL[15]+alpha[1]*fhatAL[14])*dfac_z; 
  incr[29] = -0.3061862178478971*(alpha[0]*fhatAL[13]+alpha[1]*fhatAL[10])*dfac_z; 
  incr[30] = -0.3061862178478971*(alpha[1]*fhatAL[15]+alpha[0]*fhatAL[14])*dfac_z; 
  incr[31] = -0.3061862178478971*(alpha[0]*fhatAL[15]+alpha[1]*fhatAL[14])*dfac_z; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 
  outr[28] += incr[28]; 
  outr[29] += incr[29]; 
  outr[30] += incr[30]; 
  outr[31] += incr[31]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += incr[19]; 
  outl[20] += -1.0*incr[20]; 
  outl[21] += incr[21]; 
  outl[22] += incr[22]; 
  outl[23] += -1.0*incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += incr[25]; 
  outl[26] += incr[26]; 
  outl[27] += incr[27]; 
  outl[28] += -1.0*incr[28]; 
  outl[29] += incr[29]; 
  outl[30] += incr[30]; 
  outl[31] += incr[31]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity3x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_z = 2.0/dxv[2]; 
  double dfac_v = 2.0/dxv[3]; 
  double dfac_m = 2.0/dxv[4]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wz = w[2]; 
  double wv = w[3]; 
  double wm = w[4]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[32]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.1082531754730548*(m_*(dfac_v*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[4]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*wv-1.0*BdriftX[0]*Bmag[1]*dfac_x*wm)+q_*((Gradpar[1]*Phi[5]+Gradpar[0]*Phi[3])*dfac_v*dfac_z*q_-1.0*((BdriftY[1]*Phi[4]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*m_)))/(dfac_v*m_*q_); 

  double alpha[16]; 
  alpha[0] = -(0.8660254037844386*(m_*(BdriftX[0]*Bmag[1]*dfac_x*wm+((BdriftY[1]*Phi[4]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*q_)*(dfac_v*wv-1.0)+(Gradpar[1]*Phi[5]+Gradpar[0]*Phi[3])*dfac_v*dfac_z*q2))/(dfac_v*m_*q_); 
  alpha[1] = -(0.8660254037844386*(m_*(BdriftX[1]*Bmag[1]*dfac_x*wm+((BdriftY[0]*Phi[4]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*q_)*(dfac_v*wv-1.0)+(Gradpar[0]*Phi[5]+Gradpar[1]*Phi[3])*dfac_v*dfac_z*q2))/(dfac_v*m_*q_); 
  alpha[2] = -(0.8660254037844386*(dfac_v*(BdriftX[0]*Phi[4]*dfac_x*m_*wv+(Gradpar[1]*Phi[7]+Gradpar[0]*Phi[6])*dfac_z*q_)-1.0*BdriftX[0]*Phi[4]*dfac_x*m_))/(dfac_v*m_); 
  alpha[3] = -(0.8660254037844386*((BdriftY[1]*Phi[7]+BdriftY[0]*Phi[6])*dfac_y+BdriftX[0]*Phi[5]*dfac_x)*(dfac_v*wv-1.0))/dfac_v; 
  alpha[4] = -(0.5*BdriftX[0]*Bmag[1]*dfac_x*(dfac_v*wv-1.0))/(dfac_m*dfac_v*q_); 
  alpha[5] = -(0.8660254037844386*(dfac_v*(BdriftX[1]*Phi[4]*dfac_x*m_*wv+(Gradpar[0]*Phi[7]+Gradpar[1]*Phi[6])*dfac_z*q_)-1.0*BdriftX[1]*Phi[4]*dfac_x*m_))/(dfac_v*m_); 
  alpha[6] = -(0.8660254037844386*((BdriftY[0]*Phi[7]+BdriftY[1]*Phi[6])*dfac_y+BdriftX[1]*Phi[5]*dfac_x)*(dfac_v*wv-1.0))/dfac_v; 
  alpha[7] = -(0.8660254037844386*BdriftX[0]*Phi[7]*dfac_x*(dfac_v*wv-1.0))/dfac_v; 
  alpha[8] = -(0.5*BdriftX[1]*Bmag[1]*dfac_x*(dfac_v*wv-1.0))/(dfac_m*dfac_v*q_); 
  alpha[11] = -(0.8660254037844386*BdriftX[1]*Phi[7]*dfac_x*(dfac_v*wv-1.0))/dfac_v; 
  double f0Quad[16]; 
  double f1Quad[16]; 
  double limQuad[16]; 
  // determine upwinding at each surface quadrature node 
  if((-0.25*alpha[11])+0.25*(alpha[8]+alpha[7]+alpha[6]+alpha[5])-0.25*(alpha[4]+alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[0] = 0.25*(fl[27]-1.0*(fl[22]+fl[21]+fl[20]+fl[16])+fl[14]+fl[13]+fl[12]+fl[8]+fl[7]+fl[6]-1.0*(fl[5]+fl[3]+fl[2]+fl[1])+fl[0]); 
    f1Quad[0] = -0.25*(fl[31]-1.0*(fl[30]+fl[29]+fl[28]+fl[26])+fl[25]+fl[24]+fl[23]+fl[19]+fl[18]+fl[17]-1.0*(fl[15]+fl[11]+fl[10]+fl[9])+fl[4]); 
    limQuad[0] = fl[0]/cflL; 
  } else {
    f0Quad[0] = 0.25*(fr[27]-1.0*(fr[22]+fr[21]+fr[20]+fr[16])+fr[14]+fr[13]+fr[12]+fr[8]+fr[7]+fr[6]-1.0*(fr[5]+fr[3]+fr[2]+fr[1])+fr[0]); 
    f1Quad[0] = 0.25*(fr[31]-1.0*(fr[30]+fr[29]+fr[28]+fr[26])+fr[25]+fr[24]+fr[23]+fr[19]+fr[18]+fr[17]-1.0*(fr[15]+fr[11]+fr[10]+fr[9])+fr[4]); 
    limQuad[0] = fr[0]/cflR; 
  }
  if(0.25*alpha[11]-0.25*alpha[8]+0.25*alpha[7]-0.25*(alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[1] = -0.25*(fl[27]+fl[22]-1.0*(fl[21]+fl[20]+fl[16]+fl[14]+fl[13])+fl[12]-1.0*fl[8]+fl[7]+fl[6]+fl[5]+fl[3]+fl[2]-1.0*(fl[1]+fl[0])); 
    f1Quad[1] = 0.25*(fl[31]+fl[30]-1.0*(fl[29]+fl[28]+fl[26]+fl[25]+fl[24])+fl[23]-1.0*fl[19]+fl[18]+fl[17]+fl[15]+fl[11]+fl[10]-1.0*(fl[9]+fl[4])); 
    limQuad[1] = fl[0]/cflL; 
  } else {
    f0Quad[1] = -0.25*(fr[27]+fr[22]-1.0*(fr[21]+fr[20]+fr[16]+fr[14]+fr[13])+fr[12]-1.0*fr[8]+fr[7]+fr[6]+fr[5]+fr[3]+fr[2]-1.0*(fr[1]+fr[0])); 
    f1Quad[1] = -0.25*(fr[31]+fr[30]-1.0*(fr[29]+fr[28]+fr[26]+fr[25]+fr[24])+fr[23]-1.0*fr[19]+fr[18]+fr[17]+fr[15]+fr[11]+fr[10]-1.0*(fr[9]+fr[4])); 
    limQuad[1] = fr[0]/cflR; 
  }
  if(0.25*(alpha[11]+alpha[8])-0.25*alpha[7]+0.25*alpha[6]-0.25*(alpha[5]+alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[2] = -0.25*(fl[27]-1.0*fl[22]+fl[21]-1.0*(fl[20]+fl[16]+fl[14])+fl[13]-1.0*fl[12]+fl[8]-1.0*fl[7]+fl[6]+fl[5]+fl[3]-1.0*fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[2] = 0.25*(fl[31]-1.0*fl[30]+fl[29]-1.0*(fl[28]+fl[26]+fl[25])+fl[24]-1.0*fl[23]+fl[19]-1.0*fl[18]+fl[17]+fl[15]+fl[11]-1.0*fl[10]+fl[9]-1.0*fl[4]); 
    limQuad[2] = fl[0]/cflL; 
  } else {
    f0Quad[2] = -0.25*(fr[27]-1.0*fr[22]+fr[21]-1.0*(fr[20]+fr[16]+fr[14])+fr[13]-1.0*fr[12]+fr[8]-1.0*fr[7]+fr[6]+fr[5]+fr[3]-1.0*fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[2] = -0.25*(fr[31]-1.0*fr[30]+fr[29]-1.0*(fr[28]+fr[26]+fr[25])+fr[24]-1.0*fr[23]+fr[19]-1.0*fr[18]+fr[17]+fr[15]+fr[11]-1.0*fr[10]+fr[9]-1.0*fr[4]); 
    limQuad[2] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[11]+alpha[8]+alpha[7]+alpha[6]))+0.25*alpha[5]-0.25*(alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[3] = 0.25*(fl[27]+fl[22]+fl[21]-1.0*(fl[20]+fl[16])+fl[14]-1.0*(fl[13]+fl[12]+fl[8]+fl[7])+fl[6]-1.0*(fl[5]+fl[3])+fl[2]+fl[1]+fl[0]); 
    f1Quad[3] = -0.25*(fl[31]+fl[30]+fl[29]-1.0*(fl[28]+fl[26])+fl[25]-1.0*(fl[24]+fl[23]+fl[19]+fl[18])+fl[17]-1.0*(fl[15]+fl[11])+fl[10]+fl[9]+fl[4]); 
    limQuad[3] = fl[0]/cflL; 
  } else {
    f0Quad[3] = 0.25*(fr[27]+fr[22]+fr[21]-1.0*(fr[20]+fr[16])+fr[14]-1.0*(fr[13]+fr[12]+fr[8]+fr[7])+fr[6]-1.0*(fr[5]+fr[3])+fr[2]+fr[1]+fr[0]); 
    f1Quad[3] = 0.25*(fr[31]+fr[30]+fr[29]-1.0*(fr[28]+fr[26])+fr[25]-1.0*(fr[24]+fr[23]+fr[19]+fr[18])+fr[17]-1.0*(fr[15]+fr[11])+fr[10]+fr[9]+fr[4]); 
    limQuad[3] = fr[0]/cflR; 
  }
  if(0.25*(alpha[11]+alpha[8])-0.25*(alpha[7]+alpha[6])+0.25*alpha[5]-0.25*alpha[4]+0.25*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[4] = -0.25*(fl[27]-1.0*(fl[22]+fl[21])+fl[20]-1.0*fl[16]+fl[14]-1.0*(fl[13]+fl[12])+fl[8]+fl[7]-1.0*fl[6]+fl[5]-1.0*fl[3]+fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[4] = 0.25*(fl[31]-1.0*(fl[30]+fl[29])+fl[28]-1.0*fl[26]+fl[25]-1.0*(fl[24]+fl[23])+fl[19]+fl[18]-1.0*fl[17]+fl[15]-1.0*fl[11]+fl[10]+fl[9]-1.0*fl[4]); 
    limQuad[4] = fl[0]/cflL; 
  } else {
    f0Quad[4] = -0.25*(fr[27]-1.0*(fr[22]+fr[21])+fr[20]-1.0*fr[16]+fr[14]-1.0*(fr[13]+fr[12])+fr[8]+fr[7]-1.0*fr[6]+fr[5]-1.0*fr[3]+fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[4] = -0.25*(fr[31]-1.0*(fr[30]+fr[29])+fr[28]-1.0*fr[26]+fr[25]-1.0*(fr[24]+fr[23])+fr[19]+fr[18]-1.0*fr[17]+fr[15]-1.0*fr[11]+fr[10]+fr[9]-1.0*fr[4]); 
    limQuad[4] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[11]+alpha[8]+alpha[7]))+0.25*alpha[6]-0.25*(alpha[5]+alpha[4])+0.25*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[5] = 0.25*(fl[27]+fl[22]-1.0*fl[21]+fl[20]-1.0*(fl[16]+fl[14])+fl[13]-1.0*(fl[12]+fl[8])+fl[7]-1.0*(fl[6]+fl[5])+fl[3]-1.0*fl[2]+fl[1]+fl[0]); 
    f1Quad[5] = -0.25*(fl[31]+fl[30]-1.0*fl[29]+fl[28]-1.0*(fl[26]+fl[25])+fl[24]-1.0*(fl[23]+fl[19])+fl[18]-1.0*(fl[17]+fl[15])+fl[11]-1.0*fl[10]+fl[9]+fl[4]); 
    limQuad[5] = fl[0]/cflL; 
  } else {
    f0Quad[5] = 0.25*(fr[27]+fr[22]-1.0*fr[21]+fr[20]-1.0*(fr[16]+fr[14])+fr[13]-1.0*(fr[12]+fr[8])+fr[7]-1.0*(fr[6]+fr[5])+fr[3]-1.0*fr[2]+fr[1]+fr[0]); 
    f1Quad[5] = 0.25*(fr[31]+fr[30]-1.0*fr[29]+fr[28]-1.0*(fr[26]+fr[25])+fr[24]-1.0*(fr[23]+fr[19])+fr[18]-1.0*(fr[17]+fr[15])+fr[11]-1.0*fr[10]+fr[9]+fr[4]); 
    limQuad[5] = fr[0]/cflR; 
  }
  if((-0.25*alpha[11])+0.25*(alpha[8]+alpha[7])-0.25*(alpha[6]+alpha[5]+alpha[4])+0.25*(alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[6] = 0.25*(fl[27]-1.0*fl[22]+fl[21]+fl[20]-1.0*(fl[16]+fl[14]+fl[13])+fl[12]+fl[8]-1.0*(fl[7]+fl[6]+fl[5])+fl[3]+fl[2]-1.0*fl[1]+fl[0]); 
    f1Quad[6] = -0.25*(fl[31]-1.0*fl[30]+fl[29]+fl[28]-1.0*(fl[26]+fl[25]+fl[24])+fl[23]+fl[19]-1.0*(fl[18]+fl[17]+fl[15])+fl[11]+fl[10]-1.0*fl[9]+fl[4]); 
    limQuad[6] = fl[0]/cflL; 
  } else {
    f0Quad[6] = 0.25*(fr[27]-1.0*fr[22]+fr[21]+fr[20]-1.0*(fr[16]+fr[14]+fr[13])+fr[12]+fr[8]-1.0*(fr[7]+fr[6]+fr[5])+fr[3]+fr[2]-1.0*fr[1]+fr[0]); 
    f1Quad[6] = 0.25*(fr[31]-1.0*fr[30]+fr[29]+fr[28]-1.0*(fr[26]+fr[25]+fr[24])+fr[23]+fr[19]-1.0*(fr[18]+fr[17]+fr[15])+fr[11]+fr[10]-1.0*fr[9]+fr[4]); 
    limQuad[6] = fr[0]/cflR; 
  }
  if(0.25*alpha[11]-0.25*alpha[8]+0.25*(alpha[7]+alpha[6]+alpha[5])-0.25*alpha[4]+0.25*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[7] = -0.25*(fl[27]+fl[22]+fl[21]+fl[20]-1.0*fl[16]+fl[14]+fl[13]+fl[12]-1.0*(fl[8]+fl[7]+fl[6])+fl[5]-1.0*(fl[3]+fl[2]+fl[1]+fl[0])); 
    f1Quad[7] = 0.25*(fl[31]+fl[30]+fl[29]+fl[28]-1.0*fl[26]+fl[25]+fl[24]+fl[23]-1.0*(fl[19]+fl[18]+fl[17])+fl[15]-1.0*(fl[11]+fl[10]+fl[9]+fl[4])); 
    limQuad[7] = fl[0]/cflL; 
  } else {
    f0Quad[7] = -0.25*(fr[27]+fr[22]+fr[21]+fr[20]-1.0*fr[16]+fr[14]+fr[13]+fr[12]-1.0*(fr[8]+fr[7]+fr[6])+fr[5]-1.0*(fr[3]+fr[2]+fr[1]+fr[0])); 
    f1Quad[7] = -0.25*(fr[31]+fr[30]+fr[29]+fr[28]-1.0*fr[26]+fr[25]+fr[24]+fr[23]-1.0*(fr[19]+fr[18]+fr[17])+fr[15]-1.0*(fr[11]+fr[10]+fr[9]+fr[4])); 
    limQuad[7] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[11]+alpha[8]))+0.25*(alpha[7]+alpha[6]+alpha[5]+alpha[4])-0.25*(alpha[3]+alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[8] = -0.25*(fl[27]-1.0*(fl[22]+fl[21]+fl[20])+fl[16]+fl[14]+fl[13]+fl[12]-1.0*(fl[8]+fl[7]+fl[6]+fl[5])+fl[3]+fl[2]+fl[1]-1.0*fl[0]); 
    f1Quad[8] = 0.25*(fl[31]-1.0*(fl[30]+fl[29]+fl[28])+fl[26]+fl[25]+fl[24]+fl[23]-1.0*(fl[19]+fl[18]+fl[17]+fl[15])+fl[11]+fl[10]+fl[9]-1.0*fl[4]); 
    limQuad[8] = fl[0]/cflL; 
  } else {
    f0Quad[8] = -0.25*(fr[27]-1.0*(fr[22]+fr[21]+fr[20])+fr[16]+fr[14]+fr[13]+fr[12]-1.0*(fr[8]+fr[7]+fr[6]+fr[5])+fr[3]+fr[2]+fr[1]-1.0*fr[0]); 
    f1Quad[8] = -0.25*(fr[31]-1.0*(fr[30]+fr[29]+fr[28])+fr[26]+fr[25]+fr[24]+fr[23]-1.0*(fr[19]+fr[18]+fr[17]+fr[15])+fr[11]+fr[10]+fr[9]-1.0*fr[4]); 
    limQuad[8] = fr[0]/cflR; 
  }
  if(0.25*(alpha[11]+alpha[8]+alpha[7])-0.25*(alpha[6]+alpha[5])+0.25*alpha[4]-0.25*(alpha[3]+alpha[2])+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[9] = 0.25*(fl[27]+fl[22]-1.0*(fl[21]+fl[20])+fl[16]-1.0*(fl[14]+fl[13])+fl[12]+fl[8]-1.0*(fl[7]+fl[6])+fl[5]-1.0*(fl[3]+fl[2])+fl[1]+fl[0]); 
    f1Quad[9] = -0.25*(fl[31]+fl[30]-1.0*(fl[29]+fl[28])+fl[26]-1.0*(fl[25]+fl[24])+fl[23]+fl[19]-1.0*(fl[18]+fl[17])+fl[15]-1.0*(fl[11]+fl[10])+fl[9]+fl[4]); 
    limQuad[9] = fl[0]/cflL; 
  } else {
    f0Quad[9] = 0.25*(fr[27]+fr[22]-1.0*(fr[21]+fr[20])+fr[16]-1.0*(fr[14]+fr[13])+fr[12]+fr[8]-1.0*(fr[7]+fr[6])+fr[5]-1.0*(fr[3]+fr[2])+fr[1]+fr[0]); 
    f1Quad[9] = 0.25*(fr[31]+fr[30]-1.0*(fr[29]+fr[28])+fr[26]-1.0*(fr[25]+fr[24])+fr[23]+fr[19]-1.0*(fr[18]+fr[17])+fr[15]-1.0*(fr[11]+fr[10])+fr[9]+fr[4]); 
    limQuad[9] = fr[0]/cflR; 
  }
  if(0.25*alpha[11]-0.25*(alpha[8]+alpha[7])+0.25*alpha[6]-0.25*alpha[5]+0.25*alpha[4]-0.25*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[10] = 0.25*(fl[27]-1.0*fl[22]+fl[21]-1.0*fl[20]+fl[16]-1.0*fl[14]+fl[13]-1.0*(fl[12]+fl[8])+fl[7]-1.0*fl[6]+fl[5]-1.0*fl[3]+fl[2]-1.0*fl[1]+fl[0]); 
    f1Quad[10] = -0.25*(fl[31]-1.0*fl[30]+fl[29]-1.0*fl[28]+fl[26]-1.0*fl[25]+fl[24]-1.0*(fl[23]+fl[19])+fl[18]-1.0*fl[17]+fl[15]-1.0*fl[11]+fl[10]-1.0*fl[9]+fl[4]); 
    limQuad[10] = fl[0]/cflL; 
  } else {
    f0Quad[10] = 0.25*(fr[27]-1.0*fr[22]+fr[21]-1.0*fr[20]+fr[16]-1.0*fr[14]+fr[13]-1.0*(fr[12]+fr[8])+fr[7]-1.0*fr[6]+fr[5]-1.0*fr[3]+fr[2]-1.0*fr[1]+fr[0]); 
    f1Quad[10] = 0.25*(fr[31]-1.0*fr[30]+fr[29]-1.0*fr[28]+fr[26]-1.0*fr[25]+fr[24]-1.0*(fr[23]+fr[19])+fr[18]-1.0*fr[17]+fr[15]-1.0*fr[11]+fr[10]-1.0*fr[9]+fr[4]); 
    limQuad[10] = fr[0]/cflR; 
  }
  if((-0.25*alpha[11])+0.25*alpha[8]-0.25*(alpha[7]+alpha[6])+0.25*(alpha[5]+alpha[4])-0.25*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[11] = -0.25*(fl[27]+fl[22]+fl[21]-1.0*fl[20]+fl[16]+fl[14]-1.0*(fl[13]+fl[12])+fl[8]+fl[7]-1.0*(fl[6]+fl[5])+fl[3]-1.0*(fl[2]+fl[1]+fl[0])); 
    f1Quad[11] = 0.25*(fl[31]+fl[30]+fl[29]-1.0*fl[28]+fl[26]+fl[25]-1.0*(fl[24]+fl[23])+fl[19]+fl[18]-1.0*(fl[17]+fl[15])+fl[11]-1.0*(fl[10]+fl[9]+fl[4])); 
    limQuad[11] = fl[0]/cflL; 
  } else {
    f0Quad[11] = -0.25*(fr[27]+fr[22]+fr[21]-1.0*fr[20]+fr[16]+fr[14]-1.0*(fr[13]+fr[12])+fr[8]+fr[7]-1.0*(fr[6]+fr[5])+fr[3]-1.0*(fr[2]+fr[1]+fr[0])); 
    f1Quad[11] = -0.25*(fr[31]+fr[30]+fr[29]-1.0*fr[28]+fr[26]+fr[25]-1.0*(fr[24]+fr[23])+fr[19]+fr[18]-1.0*(fr[17]+fr[15])+fr[11]-1.0*(fr[10]+fr[9]+fr[4])); 
    limQuad[11] = fr[0]/cflR; 
  }
  if(0.25*alpha[11]-0.25*(alpha[8]+alpha[7]+alpha[6])+0.25*(alpha[5]+alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) {
    f0Quad[12] = 0.25*(fl[27]-1.0*(fl[22]+fl[21])+fl[20]+fl[16]+fl[14]-1.0*(fl[13]+fl[12]+fl[8]+fl[7])+fl[6]+fl[5]+fl[3]-1.0*(fl[2]+fl[1])+fl[0]); 
    f1Quad[12] = -0.25*(fl[31]-1.0*(fl[30]+fl[29])+fl[28]+fl[26]+fl[25]-1.0*(fl[24]+fl[23]+fl[19]+fl[18])+fl[17]+fl[15]+fl[11]-1.0*(fl[10]+fl[9])+fl[4]); 
    limQuad[12] = fl[0]/cflL; 
  } else {
    f0Quad[12] = 0.25*(fr[27]-1.0*(fr[22]+fr[21])+fr[20]+fr[16]+fr[14]-1.0*(fr[13]+fr[12]+fr[8]+fr[7])+fr[6]+fr[5]+fr[3]-1.0*(fr[2]+fr[1])+fr[0]); 
    f1Quad[12] = 0.25*(fr[31]-1.0*(fr[30]+fr[29])+fr[28]+fr[26]+fr[25]-1.0*(fr[24]+fr[23]+fr[19]+fr[18])+fr[17]+fr[15]+fr[11]-1.0*(fr[10]+fr[9])+fr[4]); 
    limQuad[12] = fr[0]/cflR; 
  }
  if((-0.25*alpha[11])+0.25*alpha[8]-0.25*alpha[7]+0.25*alpha[6]-0.25*alpha[5]+0.25*(alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) {
    f0Quad[13] = -0.25*(fl[27]+fl[22]-1.0*fl[21]+fl[20]+fl[16]-1.0*fl[14]+fl[13]-1.0*fl[12]+fl[8]-1.0*fl[7]+fl[6]-1.0*(fl[5]+fl[3])+fl[2]-1.0*(fl[1]+fl[0])); 
    f1Quad[13] = 0.25*(fl[31]+fl[30]-1.0*fl[29]+fl[28]+fl[26]-1.0*fl[25]+fl[24]-1.0*fl[23]+fl[19]-1.0*fl[18]+fl[17]-1.0*(fl[15]+fl[11])+fl[10]-1.0*(fl[9]+fl[4])); 
    limQuad[13] = fl[0]/cflL; 
  } else {
    f0Quad[13] = -0.25*(fr[27]+fr[22]-1.0*fr[21]+fr[20]+fr[16]-1.0*fr[14]+fr[13]-1.0*fr[12]+fr[8]-1.0*fr[7]+fr[6]-1.0*(fr[5]+fr[3])+fr[2]-1.0*(fr[1]+fr[0])); 
    f1Quad[13] = -0.25*(fr[31]+fr[30]-1.0*fr[29]+fr[28]+fr[26]-1.0*fr[25]+fr[24]-1.0*fr[23]+fr[19]-1.0*fr[18]+fr[17]-1.0*(fr[15]+fr[11])+fr[10]-1.0*(fr[9]+fr[4])); 
    limQuad[13] = fr[0]/cflR; 
  }
  if((-0.25*(alpha[11]+alpha[8]))+0.25*alpha[7]-0.25*(alpha[6]+alpha[5])+0.25*(alpha[4]+alpha[3]+alpha[2])-0.25*alpha[1]+0.25*alpha[0] > 0) {
    f0Quad[14] = -0.25*(fl[27]-1.0*fl[22]+fl[21]+fl[20]+fl[16]-1.0*(fl[14]+fl[13])+fl[12]-1.0*fl[8]+fl[7]+fl[6]-1.0*(fl[5]+fl[3]+fl[2])+fl[1]-1.0*fl[0]); 
    f1Quad[14] = 0.25*(fl[31]-1.0*fl[30]+fl[29]+fl[28]+fl[26]-1.0*(fl[25]+fl[24])+fl[23]-1.0*fl[19]+fl[18]+fl[17]-1.0*(fl[15]+fl[11]+fl[10])+fl[9]-1.0*fl[4]); 
    limQuad[14] = fl[0]/cflL; 
  } else {
    f0Quad[14] = -0.25*(fr[27]-1.0*fr[22]+fr[21]+fr[20]+fr[16]-1.0*(fr[14]+fr[13])+fr[12]-1.0*fr[8]+fr[7]+fr[6]-1.0*(fr[5]+fr[3]+fr[2])+fr[1]-1.0*fr[0]); 
    f1Quad[14] = -0.25*(fr[31]-1.0*fr[30]+fr[29]+fr[28]+fr[26]-1.0*(fr[25]+fr[24])+fr[23]-1.0*fr[19]+fr[18]+fr[17]-1.0*(fr[15]+fr[11]+fr[10])+fr[9]-1.0*fr[4]); 
    limQuad[14] = fr[0]/cflR; 
  }
  if(0.25*(alpha[11]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) {
    f0Quad[15] = 0.25*(fl[27]+fl[22]+fl[21]+fl[20]+fl[16]+fl[14]+fl[13]+fl[12]+fl[8]+fl[7]+fl[6]+fl[5]+fl[3]+fl[2]+fl[1]+fl[0]); 
    f1Quad[15] = -0.25*(fl[31]+fl[30]+fl[29]+fl[28]+fl[26]+fl[25]+fl[24]+fl[23]+fl[19]+fl[18]+fl[17]+fl[15]+fl[11]+fl[10]+fl[9]+fl[4]); 
    limQuad[15] = fl[0]/cflL; 
  } else {
    f0Quad[15] = 0.25*(fr[27]+fr[22]+fr[21]+fr[20]+fr[16]+fr[14]+fr[13]+fr[12]+fr[8]+fr[7]+fr[6]+fr[5]+fr[3]+fr[2]+fr[1]+fr[0]); 
    f1Quad[15] = 0.25*(fr[31]+fr[30]+fr[29]+fr[28]+fr[26]+fr[25]+fr[24]+fr[23]+fr[19]+fr[18]+fr[17]+fr[15]+fr[11]+fr[10]+fr[9]+fr[4]); 
    limQuad[15] = fr[0]/cflR; 
  }
  double fhat[32]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[3] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8])+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[5] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[6] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[7] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]+f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[8] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[9] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[10] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[11] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8])+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[12] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*f0Quad[12]+f0Quad[11]-1.0*f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[13] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12])+f0Quad[11]+f0Quad[10]-1.0*(f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[14] = 0.25*(f0Quad[15]+f0Quad[14]+f0Quad[13]+f0Quad[12]-1.0*(f0Quad[11]+f0Quad[10]+f0Quad[9]+f0Quad[8]+f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[15] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[16] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*f0Quad[8]+f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[17] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[18] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]+f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[19] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[20] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]+f0Quad[11]-1.0*(f0Quad[10]+f0Quad[9])+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[21] = 0.25*(f0Quad[15]-1.0*f0Quad[14]+f0Quad[13]-1.0*(f0Quad[12]+f0Quad[11])+f0Quad[10]-1.0*f0Quad[9]+f0Quad[8]-1.0*f0Quad[7]+f0Quad[6]-1.0*f0Quad[5]+f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[22] = 0.25*(f0Quad[15]+f0Quad[14]-1.0*(f0Quad[13]+f0Quad[12]+f0Quad[11]+f0Quad[10])+f0Quad[9]+f0Quad[8]-1.0*(f0Quad[7]+f0Quad[6])+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[23] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*f1Quad[12]+f1Quad[11]-1.0*f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[24] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12])+f1Quad[11]+f1Quad[10]-1.0*(f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[25] = 0.25*(f1Quad[15]+f1Quad[14]+f1Quad[13]+f1Quad[12]-1.0*(f1Quad[11]+f1Quad[10]+f1Quad[9]+f1Quad[8]+f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[26] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*f1Quad[8]+f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[27] = 0.25*(f0Quad[15]-1.0*(f0Quad[14]+f0Quad[13])+f0Quad[12]-1.0*f0Quad[11]+f0Quad[10]+f0Quad[9]-1.0*(f0Quad[8]+f0Quad[7])+f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[28] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]+f1Quad[11]-1.0*(f1Quad[10]+f1Quad[9])+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[29] = 0.25*(f1Quad[15]-1.0*f1Quad[14]+f1Quad[13]-1.0*(f1Quad[12]+f1Quad[11])+f1Quad[10]-1.0*f1Quad[9]+f1Quad[8]-1.0*f1Quad[7]+f1Quad[6]-1.0*f1Quad[5]+f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[30] = 0.25*(f1Quad[15]+f1Quad[14]-1.0*(f1Quad[13]+f1Quad[12]+f1Quad[11]+f1Quad[10])+f1Quad[9]+f1Quad[8]-1.0*(f1Quad[7]+f1Quad[6])+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[31] = 0.25*(f1Quad[15]-1.0*(f1Quad[14]+f1Quad[13])+f1Quad[12]-1.0*f1Quad[11]+f1Quad[10]+f1Quad[9]-1.0*(f1Quad[8]+f1Quad[7])+f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[16];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[24]+fhat[23]+fhat[19]+fhat[18]+fhat[17]+3.0*fhat[4]))-3.0*(fhat[30]+fhat[29]+fhat[28]+fhat[26]+3.0*(fhat[15]+fhat[11]+fhat[10]+fhat[9]))))/(72.0*EPSILON+1.414213562373095*(fhat[27]-1.732050807568877*(fhat[22]+fhat[21]+fhat[20]+fhat[16]+3.0*(fhat[5]+fhat[3]+fhat[2]+fhat[1]))+3.0*(fhat[14]+fhat[13]+fhat[12]+fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0]))); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[23]-1.0*fhat[19]+fhat[18]+fhat[17]-3.0*fhat[4]))+3.0*(fhat[30]-1.0*(fhat[29]+fhat[28]+fhat[26]-3.0*(fhat[15]+fhat[11]+fhat[10]-1.0*fhat[9])))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[14]+fhat[13]-1.0*fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0]))-1.0*(fhat[27]+1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20]+fhat[16])+3.0*(fhat[5]+fhat[3]+fhat[2]-1.0*fhat[1]))))); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*fhat[23]+fhat[19]-1.0*fhat[18]+fhat[17]-3.0*fhat[4]))+3.0*((-1.0*fhat[30])+fhat[29]-1.0*(fhat[28]+fhat[26]-3.0*(fhat[15]+fhat[11]-1.0*fhat[10]+fhat[9])))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[27])+1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]+fhat[16]+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]-1.0*fhat[1]))+3.0*(fhat[14]-1.0*fhat[13]+fhat[12]-1.0*fhat[8]+fhat[7]-1.0*fhat[6]+3.0*fhat[0]))); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[23]+fhat[19]+fhat[18]-1.0*(fhat[17]+3.0*fhat[4]))))+3.0*(fhat[30]+fhat[29]-1.0*(fhat[28]+fhat[26]-3.0*((-1.0*(fhat[15]+fhat[11]))+fhat[10]+fhat[9])))))/(72.0*EPSILON+1.414213562373095*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]-1.0*(fhat[20]+fhat[16])+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]+fhat[1]))+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]+fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[4] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[23]-1.0*fhat[19])+fhat[18]-1.0*(fhat[17]+3.0*fhat[4])))+3.0*((-1.0*(fhat[30]+fhat[29]))+fhat[28]-1.0*fhat[26]+3.0*(fhat[15]-1.0*fhat[11]+fhat[10]+fhat[9]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[27])+1.732050807568877*(fhat[22]+fhat[21]-1.0*fhat[20]+fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*(fhat[2]+fhat[1])))+3.0*((-1.0*fhat[14])+fhat[13]+fhat[12]-1.0*(fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[5] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*(fhat[23]+fhat[19]-1.0*fhat[18]+fhat[17])+3.0*fhat[4]))+3.0*(fhat[30]-1.0*fhat[29]+fhat[28]-1.0*fhat[26]+3.0*((-1.0*fhat[15])+fhat[11]-1.0*fhat[10]+fhat[9]))))/(72.0*EPSILON+1.414213562373095*(fhat[27]+1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]-1.0*fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[14])+fhat[13]-1.0*(fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6])+3.0*fhat[0]))); 
  rCtrl[6] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[23]+fhat[19]-1.0*(fhat[18]+fhat[17]-3.0*fhat[4])))+3.0*((-1.0*fhat[30])+fhat[29]+fhat[28]-1.0*fhat[26]+3.0*((-1.0*fhat[15])+fhat[11]+fhat[10]-1.0*fhat[9]))))/(72.0*EPSILON+1.414213562373095*(fhat[27]-1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20])+fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[7] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[24]+fhat[23]-1.0*(fhat[19]+fhat[18]+fhat[17]+3.0*fhat[4])))+3.0*(fhat[30]+fhat[29]+fhat[28]-1.0*fhat[26]+3.0*(fhat[15]-1.0*(fhat[11]+fhat[10]+fhat[9])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*(fhat[14]+fhat[13]+fhat[12]-1.0*fhat[8]))+fhat[7]+fhat[6]+3.0*fhat[0])-1.0*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]+fhat[20]-1.0*fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]+fhat[1])))))); 
  rCtrl[8] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[24]+fhat[23]-1.0*(fhat[19]+fhat[18]+fhat[17]+3.0*fhat[4])))+3.0*(3.0*((-1.0*fhat[15])+fhat[11]+fhat[10]+fhat[9])-1.0*(fhat[30]+fhat[29]+fhat[28]-1.0*fhat[26]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[27])+1.732050807568877*(fhat[22]+fhat[21]+fhat[20]-1.0*fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]+fhat[1])))+3.0*((-1.0*(fhat[14]+fhat[13]+fhat[12]-1.0*fhat[8]))+fhat[7]+fhat[6]+3.0*fhat[0]))); 
  rCtrl[9] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[23]+fhat[19]-1.0*(fhat[18]+fhat[17]-3.0*fhat[4])))+3.0*(fhat[30]-1.0*(fhat[29]+fhat[28]-1.0*fhat[26])+3.0*(fhat[15]-1.0*(fhat[11]+fhat[10]-1.0*fhat[9])))))/(72.0*EPSILON+1.414213562373095*(fhat[27]+1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20])+fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[10] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*(fhat[23]+fhat[19]-1.0*fhat[18]+fhat[17])+3.0*fhat[4]))+3.0*((-1.0*fhat[30])+fhat[29]-1.0*fhat[28]+fhat[26]+3.0*(fhat[15]-1.0*fhat[11]+fhat[10]-1.0*fhat[9]))))/(72.0*EPSILON+1.414213562373095*(fhat[27]-1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]-1.0*fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[14])+fhat[13]-1.0*(fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6])+3.0*fhat[0]))); 
  rCtrl[11] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[23]-1.0*fhat[19])+fhat[18]-1.0*(fhat[17]+3.0*fhat[4])))+3.0*(fhat[30]+fhat[29]-1.0*fhat[28]+fhat[26]+3.0*((-1.0*fhat[15])+fhat[11]-1.0*(fhat[10]+fhat[9])))))/(72.0*EPSILON+1.414213562373095*(3.0*((-1.0*fhat[14])+fhat[13]+fhat[12]-1.0*(fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0])))-1.0*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]-1.0*fhat[20]+fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*(fhat[2]+fhat[1])))))); 
  rCtrl[12] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]-1.0*(fhat[24]+fhat[23]+fhat[19]+fhat[18]-1.0*(fhat[17]+3.0*fhat[4]))))+3.0*((-1.0*(fhat[30]+fhat[29]))+fhat[28]+fhat[26]+3.0*(fhat[15]+fhat[11]-1.0*(fhat[10]+fhat[9])))))/(72.0*EPSILON+1.414213562373095*(fhat[27]-1.732050807568877*(fhat[22]+fhat[21]-1.0*(fhat[20]+fhat[16])+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]+fhat[1]))+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]+fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0]))))); 
  rCtrl[13] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*fhat[25])+fhat[24]-1.0*fhat[23]+fhat[19]-1.0*fhat[18]+fhat[17]-3.0*fhat[4]))+3.0*(fhat[30]-1.0*fhat[29]+fhat[28]+fhat[26]+3.0*((-1.0*(fhat[15]+fhat[11]))+fhat[10]-1.0*fhat[9]))))/(72.0*EPSILON+1.414213562373095*(3.0*(fhat[14]-1.0*fhat[13]+fhat[12]-1.0*fhat[8]+fhat[7]-1.0*fhat[6]+3.0*fhat[0])-1.0*(fhat[27]+1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]+fhat[16]+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]-1.0*fhat[1]))))); 
  rCtrl[14] = -(1.414213562373095*(1.732050807568877*(fhat[31]+3.0*((-1.0*(fhat[25]+fhat[24]))+fhat[23]-1.0*fhat[19]+fhat[18]+fhat[17]-3.0*fhat[4]))+3.0*((-1.0*fhat[30])+fhat[29]+fhat[28]+fhat[26]-3.0*(fhat[15]+fhat[11]+fhat[10]-1.0*fhat[9]))))/(72.0*EPSILON+1.414213562373095*((-1.0*fhat[27])+1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20]+fhat[16])+3.0*(fhat[5]+fhat[3]+fhat[2]-1.0*fhat[1]))+3.0*(fhat[14]+fhat[13]-1.0*fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0])))); 
  rCtrl[15] = (1.414213562373095*(1.732050807568877*(fhat[31]+3.0*(fhat[25]+fhat[24]+fhat[23]+fhat[19]+fhat[18]+fhat[17]+3.0*fhat[4]))+3.0*(fhat[30]+fhat[29]+fhat[28]+fhat[26]+3.0*(fhat[15]+fhat[11]+fhat[10]+fhat[9]))))/(72.0*EPSILON+1.414213562373095*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]+fhat[20]+fhat[16]+3.0*(fhat[5]+fhat[3]+fhat[2]+fhat[1]))+3.0*(fhat[14]+fhat[13]+fhat[12]+fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0]))); 
  double fhatCtrl[16];  // fhatCtrl = anti-limited fhat evaluated at each control node on vx surface 
  fhatCtrl[0] = -0.01964185503295965*(1.732050807568877*(fhat[22]+fhat[21]+fhat[20]+fhat[16]+3.0*(fhat[5]+fhat[3]+fhat[2]+fhat[1]))-1.0*(fhat[27]+3.0*(fhat[14]+fhat[13]+fhat[12]+fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[0],-1.0,EPSILON); 
  fhatCtrl[1] = -0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20]+fhat[16])+3.0*(fhat[5]+fhat[3]+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]-1.0*fhat[8]+fhat[7]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[1],-1.0,EPSILON); 
  fhatCtrl[2] = 0.01964185503295965*(1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]+fhat[16]+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]-1.0*fhat[1]))-1.0*(fhat[27]+3.0*((-1.0*fhat[14])+fhat[13]-1.0*fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[2],-1.0,EPSILON); 
  fhatCtrl[3] = 0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]-1.0*(fhat[20]+fhat[16])+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]+fhat[1]))+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]+fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[3],-1.0,EPSILON); 
  fhatCtrl[4] = 0.01964185503295965*(1.732050807568877*(fhat[22]+fhat[21]-1.0*fhat[20]+fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*(fhat[2]+fhat[1])))-1.0*(fhat[27]+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]-1.0*fhat[8])+fhat[7]-1.0*(fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[4],-1.0,EPSILON); 
  fhatCtrl[5] = 0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]-1.0*fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*fhat[2]+fhat[1]))+3.0*((-1.0*fhat[14])+fhat[13]-1.0*(fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6])+3.0*fhat[0]))*limTheta(rCtrl[5],-1.0,EPSILON); 
  fhatCtrl[6] = -0.01964185503295965*(1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20])+fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]-1.0*fhat[1])))-1.0*(fhat[27]+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0]))))*limTheta(rCtrl[6],-1.0,EPSILON); 
  fhatCtrl[7] = -0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]+fhat[20]-1.0*fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]+fhat[1])))+3.0*(fhat[14]+fhat[13]+fhat[12]-1.0*(fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[7],-1.0,EPSILON); 
  fhatCtrl[8] = 0.01964185503295965*(1.732050807568877*(fhat[22]+fhat[21]+fhat[20]-1.0*fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]+fhat[1])))-1.0*(fhat[27]+3.0*(fhat[14]+fhat[13]+fhat[12]-1.0*(fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0]))))*limTheta(rCtrl[8],-1.0,EPSILON); 
  fhatCtrl[9] = 0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20])+fhat[16]+3.0*(fhat[5]-1.0*(fhat[3]+fhat[2]-1.0*fhat[1])))+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]+fhat[8]-1.0*(fhat[7]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[9],-1.0,EPSILON); 
  fhatCtrl[10] = -0.01964185503295965*(1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]-1.0*fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*fhat[2]+fhat[1]))-1.0*(fhat[27]+3.0*((-1.0*fhat[14])+fhat[13]-1.0*(fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6])+3.0*fhat[0])))*limTheta(rCtrl[10],-1.0,EPSILON); 
  fhatCtrl[11] = -0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]-1.0*fhat[20]+fhat[16]+3.0*((-1.0*fhat[5])+fhat[3]-1.0*(fhat[2]+fhat[1])))+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]-1.0*fhat[8])+fhat[7]-1.0*(fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[11],-1.0,EPSILON); 
  fhatCtrl[12] = -0.01964185503295965*(1.732050807568877*(fhat[22]+fhat[21]-1.0*(fhat[20]+fhat[16])+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]+fhat[1]))-1.0*(fhat[27]+3.0*(fhat[14]-1.0*(fhat[13]+fhat[12]+fhat[8]+fhat[7]-1.0*(fhat[6]+3.0*fhat[0])))))*limTheta(rCtrl[12],-1.0,EPSILON); 
  fhatCtrl[13] = -0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]-1.0*fhat[21]+fhat[20]+fhat[16]+3.0*((-1.0*(fhat[5]+fhat[3]))+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*fhat[14])+fhat[13]-1.0*fhat[12]+fhat[8]-1.0*fhat[7]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[13],-1.0,EPSILON); 
  fhatCtrl[14] = 0.01964185503295965*(1.732050807568877*(fhat[22]-1.0*(fhat[21]+fhat[20]+fhat[16])+3.0*(fhat[5]+fhat[3]+fhat[2]-1.0*fhat[1]))-1.0*(fhat[27]+3.0*((-1.0*(fhat[14]+fhat[13]))+fhat[12]-1.0*fhat[8]+fhat[7]+fhat[6]-3.0*fhat[0])))*limTheta(rCtrl[14],-1.0,EPSILON); 
  fhatCtrl[15] = 0.01964185503295965*(fhat[27]+1.732050807568877*(fhat[22]+fhat[21]+fhat[20]+fhat[16]+3.0*(fhat[5]+fhat[3]+fhat[2]+fhat[1]))+3.0*(fhat[14]+fhat[13]+fhat[12]+fhat[8]+fhat[7]+fhat[6]+3.0*fhat[0]))*limTheta(rCtrl[15],-1.0,EPSILON); 
  double fhatAL[16];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.25*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.4330127018922193*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8])+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 0.4330127018922193*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[5] = 0.75*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[6] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[8] = 0.75*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*fhatCtrl[12]+fhatCtrl[11]-1.0*fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[9] = 0.75*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12])+fhatCtrl[11]+fhatCtrl[10]-1.0*(fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[10] = 0.75*(fhatCtrl[15]+fhatCtrl[14]+fhatCtrl[13]+fhatCtrl[12]-1.0*(fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]+fhatCtrl[8]+fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[11] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*fhatCtrl[8]+fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[12] = 1.299038105676658*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]+fhatCtrl[11]-1.0*(fhatCtrl[10]+fhatCtrl[9])+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[13] = 1.299038105676658*(fhatCtrl[15]-1.0*fhatCtrl[14]+fhatCtrl[13]-1.0*(fhatCtrl[12]+fhatCtrl[11])+fhatCtrl[10]-1.0*fhatCtrl[9]+fhatCtrl[8]-1.0*fhatCtrl[7]+fhatCtrl[6]-1.0*fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[14] = 1.299038105676658*(fhatCtrl[15]+fhatCtrl[14]-1.0*(fhatCtrl[13]+fhatCtrl[12]+fhatCtrl[11]+fhatCtrl[10])+fhatCtrl[9]+fhatCtrl[8]-1.0*(fhatCtrl[7]+fhatCtrl[6])+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[15] = 2.25*(fhatCtrl[15]-1.0*(fhatCtrl[14]+fhatCtrl[13])+fhatCtrl[12]-1.0*fhatCtrl[11]+fhatCtrl[10]+fhatCtrl[9]-1.0*(fhatCtrl[8]+fhatCtrl[7])+fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[16]; 
  fhatALQuad[0] = fmin(0.25*((-0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[7]))-1.732050807568877*(fhatAL[10]+fhatAL[6])))-1.0*(1.0*(fhatAL[14]+fhatAL[11])+fhatAL[13]+fhatAL[12]+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[0]); 
  fhatALQuad[1] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[10]+1.732050807568877*fhatAL[6]))+1.732050807568877*(0.5773502691896258*fhatAL[1]-0.5773502691896258*((-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[14]))-1.0*(fhatAL[13]+fhatAL[12])+fhatAL[4]+fhatAL[3]+fhatAL[2]))-0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[9])-1.0*fhatAL[8]+fhatAL[0]), limQuad[1]); 
  fhatALQuad[2] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[7])-1.732050807568877*(fhatAL[10]+fhatAL[6])))+1.732050807568877*((-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[11]))+fhatAL[13]+fhatAL[3]))+0.5773502691896258*(fhatAL[12]+fhatAL[2])-0.5773502691896258*(fhatAL[4]+fhatAL[1]))-1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[2]); 
  fhatALQuad[3] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15])-1.732050807568877*fhatAL[10]+1.732050807568877*fhatAL[6]))+1.732050807568877*((-0.5773502691896258*(0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[14])-1.0*fhatAL[13]+fhatAL[3]))+0.5773502691896258*(fhatAL[2]-1.0*fhatAL[12])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[1])+0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[9])-1.0*fhatAL[8]+fhatAL[0]), limQuad[3]); 
  fhatALQuad[4] = fmin(0.25*(0.5773502691896258*((-1.732050807568878*(fhatAL[15]+fhatAL[7]))-1.732050807568877*(fhatAL[10]+fhatAL[6]))+1.732050807568877*(0.5773502691896258*(1.0*(fhatAL[14]+fhatAL[11])+fhatAL[13]+fhatAL[3])-0.5773502691896258*(fhatAL[12]+fhatAL[4]+fhatAL[2]+fhatAL[1]))+1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[4]); 
  fhatALQuad[5] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15]))-1.732050807568877*fhatAL[10]+1.732050807568877*fhatAL[6])+1.732050807568877*(0.5773502691896258*((-0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[14]))-1.0*fhatAL[13]+fhatAL[3])-0.5773502691896258*((-1.0*fhatAL[12])+fhatAL[4]+fhatAL[2])+0.5773502691896258*fhatAL[1])-0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[9])-1.0*fhatAL[8]+fhatAL[0]), limQuad[5]); 
  fhatALQuad[6] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[7])-1.732050807568877*(fhatAL[10]+fhatAL[6]))+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[11]))+fhatAL[13]+fhatAL[12]+fhatAL[3]+fhatAL[2])-0.5773502691896258*(fhatAL[4]+fhatAL[1]))-1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[6]); 
  fhatALQuad[7] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15])-1.732050807568877*fhatAL[10]+1.732050807568877*fhatAL[6])+1.732050807568877*(0.5773502691896258*(0.5773502691896258*(1.732050807568877*fhatAL[11]-1.732050807568877*fhatAL[14])-1.0*(fhatAL[13]+fhatAL[12])+fhatAL[3]+fhatAL[2])-0.5773502691896258*fhatAL[4]+0.5773502691896258*fhatAL[1])+0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[9])-1.0*fhatAL[8]+fhatAL[0]), limQuad[7]); 
  fhatALQuad[8] = fmin(0.25*((-0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[10]-1.732050807568877*fhatAL[6]))+1.732050807568877*((-0.5773502691896258*((-0.5773502691896258*(1.732050807568877*fhatAL[14]-1.732050807568877*fhatAL[11]))-1.0*(fhatAL[13]+fhatAL[12])+fhatAL[3]+fhatAL[2]))+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[1])-0.5773502691896258*(1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[5])-1.0*fhatAL[8]+fhatAL[0]), limQuad[8]); 
  fhatALQuad[9] = fmin(0.25*((-0.5773502691896258*(1.732050807568877*(fhatAL[10]+fhatAL[6])-1.732050807568878*(fhatAL[15]+fhatAL[7])))+1.732050807568877*(0.5773502691896258*(fhatAL[4]+fhatAL[1])-0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[11]))+fhatAL[13]+fhatAL[12]+fhatAL[3]+fhatAL[2]))-1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[9]); 
  fhatALQuad[10] = fmin(0.25*((-0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15])+1.732050807568877*fhatAL[10]-1.732050807568877*fhatAL[6]))+1.732050807568877*((-0.5773502691896258*(0.5773502691896258*(1.732050807568877*fhatAL[14]-1.732050807568877*fhatAL[11])-1.0*fhatAL[13]+fhatAL[3]))+0.5773502691896258*((-1.0*fhatAL[12])+fhatAL[4]+fhatAL[2])-0.5773502691896258*fhatAL[1])+0.5773502691896258*(1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[5])-1.0*fhatAL[8]+fhatAL[0]), limQuad[10]); 
  fhatALQuad[11] = fmin(0.25*((-0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[7])+1.732050807568877*(fhatAL[10]+fhatAL[6])))+1.732050807568877*(0.5773502691896258*(fhatAL[12]+fhatAL[4]+fhatAL[2]+fhatAL[1])-0.5773502691896258*(1.0*(fhatAL[14]+fhatAL[11])+fhatAL[13]+fhatAL[3]))+1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[11]); 
  fhatALQuad[12] = fmin(0.25*(0.5773502691896258*((-0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15]))+1.732050807568877*fhatAL[10]-1.732050807568877*fhatAL[6])+1.732050807568877*(0.5773502691896258*((-0.5773502691896258*(1.732050807568877*fhatAL[14]-1.732050807568877*fhatAL[11]))-1.0*fhatAL[13]+fhatAL[3])-0.5773502691896258*(fhatAL[2]-1.0*fhatAL[12])+0.5773502691896258*fhatAL[4]-0.5773502691896258*fhatAL[1])-0.5773502691896258*(1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[5])-1.0*fhatAL[8]+fhatAL[0]), limQuad[12]); 
  fhatALQuad[13] = fmin(0.25*(0.5773502691896258*(1.732050807568877*(fhatAL[10]+fhatAL[6])-1.732050807568878*(fhatAL[15]+fhatAL[7]))+1.732050807568877*(0.5773502691896258*((-1.0*(fhatAL[14]+fhatAL[11]))+fhatAL[13]+fhatAL[3])-0.5773502691896258*(fhatAL[12]+fhatAL[2])+0.5773502691896258*(fhatAL[4]+fhatAL[1]))-1.0*(fhatAL[9]+fhatAL[5])+fhatAL[8]+fhatAL[0]), limQuad[13]); 
  fhatALQuad[14] = fmin(0.25*(0.5773502691896258*(0.5773502691896258*(3.0*fhatAL[7]-3.0*fhatAL[15])+1.732050807568877*fhatAL[10]-1.732050807568877*fhatAL[6])+1.732050807568877*(0.5773502691896258*(0.5773502691896258*(1.732050807568877*fhatAL[14]-1.732050807568877*fhatAL[11])-1.0*(fhatAL[13]+fhatAL[12])+fhatAL[4]+fhatAL[3]+fhatAL[2])-0.5773502691896258*fhatAL[1])+0.5773502691896258*(1.732050807568877*fhatAL[9]-1.732050807568877*fhatAL[5])-1.0*fhatAL[8]+fhatAL[0]), limQuad[14]); 
  fhatALQuad[15] = fmin(0.25*(0.5773502691896258*(1.732050807568878*(fhatAL[15]+fhatAL[7])+1.732050807568877*(fhatAL[10]+fhatAL[6]))+1.0*(1.0*(fhatAL[14]+fhatAL[11])+fhatAL[13]+fhatAL[12]+fhatAL[9]+fhatAL[5]+fhatAL[4]+fhatAL[3]+fhatAL[2]+fhatAL[1])+fhatAL[8]+fhatAL[0]), limQuad[15]); 
  fhatAL[0] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8])+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[5] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[6] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[8] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*fhatALQuad[12]+fhatALQuad[11]-1.0*fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[9] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12])+fhatALQuad[11]+fhatALQuad[10]-1.0*(fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[10] = 0.25*(fhatALQuad[15]+fhatALQuad[14]+fhatALQuad[13]+fhatALQuad[12]-1.0*(fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]+fhatALQuad[8]+fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[11] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*fhatALQuad[8]+fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[12] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]+fhatALQuad[11]-1.0*(fhatALQuad[10]+fhatALQuad[9])+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[13] = 0.25*(fhatALQuad[15]-1.0*fhatALQuad[14]+fhatALQuad[13]-1.0*(fhatALQuad[12]+fhatALQuad[11])+fhatALQuad[10]-1.0*fhatALQuad[9]+fhatALQuad[8]-1.0*fhatALQuad[7]+fhatALQuad[6]-1.0*fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[14] = 0.25*(fhatALQuad[15]+fhatALQuad[14]-1.0*(fhatALQuad[13]+fhatALQuad[12]+fhatALQuad[11]+fhatALQuad[10])+fhatALQuad[9]+fhatALQuad[8]-1.0*(fhatALQuad[7]+fhatALQuad[6])+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[15] = 0.25*(fhatALQuad[15]-1.0*(fhatALQuad[14]+fhatALQuad[13])+fhatALQuad[12]-1.0*fhatALQuad[11]+fhatALQuad[10]+fhatALQuad[9]-1.0*(fhatALQuad[8]+fhatALQuad[7])+fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.1767766952966368*(alpha[11]*fhatAL[11]+alpha[8]*fhatAL[8]+alpha[7]*fhatAL[7]+alpha[6]*fhatAL[6]+alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(alpha[7]*fhatAL[11]+fhatAL[7]*alpha[11]+alpha[4]*fhatAL[8]+fhatAL[4]*alpha[8]+alpha[3]*fhatAL[6]+fhatAL[3]*alpha[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(alpha[8]*fhatAL[12]+alpha[6]*fhatAL[11]+fhatAL[6]*alpha[11]+alpha[4]*fhatAL[9]+alpha[3]*fhatAL[7]+fhatAL[3]*alpha[7]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[3] = 0.1767766952966368*(alpha[8]*fhatAL[13]+alpha[5]*fhatAL[11]+fhatAL[5]*alpha[11]+alpha[4]*fhatAL[10]+alpha[2]*fhatAL[7]+fhatAL[2]*alpha[7]+alpha[1]*fhatAL[6]+fhatAL[1]*alpha[6]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_v; 
  incr[4] = -0.3061862178478971*(alpha[11]*fhatAL[11]+alpha[8]*fhatAL[8]+alpha[7]*fhatAL[7]+alpha[6]*fhatAL[6]+alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[5] = 0.1767766952966368*(alpha[11]*fhatAL[15]+alpha[7]*fhatAL[14]+alpha[6]*fhatAL[13]+alpha[5]*fhatAL[12]+alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[1]*fhatAL[8]+fhatAL[1]*alpha[8]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4])*dfac_v; 
  incr[6] = 0.1767766952966368*(alpha[4]*fhatAL[12]+alpha[3]*fhatAL[11]+fhatAL[3]*alpha[11]+alpha[8]*fhatAL[9]+alpha[6]*fhatAL[7]+fhatAL[6]*alpha[7]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 
  incr[7] = 0.1767766952966368*(alpha[4]*fhatAL[13]+alpha[2]*fhatAL[11]+fhatAL[2]*alpha[11]+alpha[8]*fhatAL[10]+alpha[5]*fhatAL[7]+fhatAL[5]*alpha[7]+alpha[0]*fhatAL[6]+fhatAL[0]*alpha[6]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_v; 
  incr[8] = 0.1767766952966368*(alpha[8]*fhatAL[15]+alpha[4]*fhatAL[14]+alpha[1]*fhatAL[11]+fhatAL[1]*alpha[11]+alpha[0]*fhatAL[7]+fhatAL[0]*alpha[7]+alpha[5]*fhatAL[6]+fhatAL[5]*alpha[6]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_v; 
  incr[9] = -0.3061862178478971*(alpha[7]*fhatAL[11]+fhatAL[7]*alpha[11]+alpha[4]*fhatAL[8]+fhatAL[4]*alpha[8]+alpha[3]*fhatAL[6]+fhatAL[3]*alpha[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[10] = -0.3061862178478971*(alpha[8]*fhatAL[12]+alpha[6]*fhatAL[11]+fhatAL[6]*alpha[11]+alpha[4]*fhatAL[9]+alpha[3]*fhatAL[7]+fhatAL[3]*alpha[7]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[11] = -0.3061862178478971*(alpha[8]*fhatAL[13]+alpha[5]*fhatAL[11]+fhatAL[5]*alpha[11]+alpha[4]*fhatAL[10]+alpha[2]*fhatAL[7]+fhatAL[2]*alpha[7]+alpha[1]*fhatAL[6]+fhatAL[1]*alpha[6]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_v; 
  incr[12] = 0.1767766952966368*(alpha[7]*fhatAL[15]+alpha[11]*fhatAL[14]+alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[6]*fhatAL[10]+alpha[5]*fhatAL[9]+alpha[0]*fhatAL[8]+fhatAL[0]*alpha[8]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4])*dfac_v; 
  incr[13] = 0.1767766952966368*(alpha[6]*fhatAL[15]+alpha[3]*fhatAL[14]+alpha[11]*fhatAL[13]+alpha[1]*fhatAL[12]+alpha[7]*fhatAL[10]+alpha[0]*fhatAL[9]+alpha[5]*fhatAL[8]+fhatAL[5]*alpha[8]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4])*dfac_v; 
  incr[14] = 0.1767766952966368*(alpha[5]*fhatAL[15]+alpha[2]*fhatAL[14]+alpha[1]*fhatAL[13]+alpha[11]*fhatAL[12]+alpha[0]*fhatAL[10]+alpha[7]*fhatAL[9]+alpha[6]*fhatAL[8]+fhatAL[6]*alpha[8]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])*dfac_v; 
  incr[15] = -0.3061862178478971*(alpha[11]*fhatAL[15]+alpha[7]*fhatAL[14]+alpha[6]*fhatAL[13]+alpha[5]*fhatAL[12]+alpha[3]*fhatAL[10]+alpha[2]*fhatAL[9]+alpha[1]*fhatAL[8]+fhatAL[1]*alpha[8]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4])*dfac_v; 
  incr[16] = 0.1767766952966368*(alpha[4]*fhatAL[15]+alpha[8]*fhatAL[14]+alpha[0]*fhatAL[11]+fhatAL[0]*alpha[11]+alpha[1]*fhatAL[7]+fhatAL[1]*alpha[7]+alpha[2]*fhatAL[6]+fhatAL[2]*alpha[6]+alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5])*dfac_v; 
  incr[17] = -0.3061862178478971*(alpha[4]*fhatAL[12]+alpha[3]*fhatAL[11]+fhatAL[3]*alpha[11]+alpha[8]*fhatAL[9]+alpha[6]*fhatAL[7]+fhatAL[6]*alpha[7]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 
  incr[18] = -0.3061862178478971*(alpha[4]*fhatAL[13]+alpha[2]*fhatAL[11]+fhatAL[2]*alpha[11]+alpha[8]*fhatAL[10]+alpha[5]*fhatAL[7]+fhatAL[5]*alpha[7]+alpha[0]*fhatAL[6]+fhatAL[0]*alpha[6]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_v; 
  incr[19] = -0.3061862178478971*(alpha[8]*fhatAL[15]+alpha[4]*fhatAL[14]+alpha[1]*fhatAL[11]+fhatAL[1]*alpha[11]+alpha[0]*fhatAL[7]+fhatAL[0]*alpha[7]+alpha[5]*fhatAL[6]+fhatAL[5]*alpha[6]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_v; 
  incr[20] = 0.1767766952966368*(alpha[3]*fhatAL[15]+alpha[6]*fhatAL[14]+alpha[7]*fhatAL[13]+alpha[0]*fhatAL[12]+fhatAL[10]*alpha[11]+alpha[1]*fhatAL[9]+alpha[2]*fhatAL[8]+fhatAL[2]*alpha[8]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5])*dfac_v; 
  incr[21] = 0.1767766952966368*(alpha[2]*fhatAL[15]+alpha[5]*fhatAL[14]+alpha[0]*fhatAL[13]+alpha[7]*fhatAL[12]+fhatAL[9]*alpha[11]+alpha[1]*fhatAL[10]+alpha[3]*fhatAL[8]+fhatAL[3]*alpha[8]+alpha[4]*fhatAL[6]+fhatAL[4]*alpha[6])*dfac_v; 
  incr[22] = 0.1767766952966368*(alpha[1]*fhatAL[15]+alpha[0]*fhatAL[14]+alpha[5]*fhatAL[13]+alpha[6]*fhatAL[12]+alpha[8]*fhatAL[11]+fhatAL[8]*alpha[11]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9]+alpha[4]*fhatAL[7]+fhatAL[4]*alpha[7])*dfac_v; 
  incr[23] = -0.3061862178478971*(alpha[7]*fhatAL[15]+alpha[11]*fhatAL[14]+alpha[3]*fhatAL[13]+alpha[2]*fhatAL[12]+alpha[6]*fhatAL[10]+alpha[5]*fhatAL[9]+alpha[0]*fhatAL[8]+fhatAL[0]*alpha[8]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4])*dfac_v; 
  incr[24] = -0.3061862178478971*(alpha[6]*fhatAL[15]+alpha[3]*fhatAL[14]+alpha[11]*fhatAL[13]+alpha[1]*fhatAL[12]+alpha[7]*fhatAL[10]+alpha[0]*fhatAL[9]+alpha[5]*fhatAL[8]+fhatAL[5]*alpha[8]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4])*dfac_v; 
  incr[25] = -0.3061862178478971*(alpha[5]*fhatAL[15]+alpha[2]*fhatAL[14]+alpha[1]*fhatAL[13]+alpha[11]*fhatAL[12]+alpha[0]*fhatAL[10]+alpha[7]*fhatAL[9]+alpha[6]*fhatAL[8]+fhatAL[6]*alpha[8]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])*dfac_v; 
  incr[26] = -0.3061862178478971*(alpha[4]*fhatAL[15]+alpha[8]*fhatAL[14]+alpha[0]*fhatAL[11]+fhatAL[0]*alpha[11]+alpha[1]*fhatAL[7]+fhatAL[1]*alpha[7]+alpha[2]*fhatAL[6]+fhatAL[2]*alpha[6]+alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5])*dfac_v; 
  incr[27] = 0.1767766952966368*(alpha[0]*fhatAL[15]+alpha[1]*fhatAL[14]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12]+alpha[4]*fhatAL[11]+fhatAL[4]*alpha[11]+alpha[5]*fhatAL[10]+alpha[6]*fhatAL[9]+alpha[7]*fhatAL[8]+fhatAL[7]*alpha[8])*dfac_v; 
  incr[28] = -0.3061862178478971*(alpha[3]*fhatAL[15]+alpha[6]*fhatAL[14]+alpha[7]*fhatAL[13]+alpha[0]*fhatAL[12]+fhatAL[10]*alpha[11]+alpha[1]*fhatAL[9]+alpha[2]*fhatAL[8]+fhatAL[2]*alpha[8]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5])*dfac_v; 
  incr[29] = -0.3061862178478971*(alpha[2]*fhatAL[15]+alpha[5]*fhatAL[14]+alpha[0]*fhatAL[13]+alpha[7]*fhatAL[12]+fhatAL[9]*alpha[11]+alpha[1]*fhatAL[10]+alpha[3]*fhatAL[8]+fhatAL[3]*alpha[8]+alpha[4]*fhatAL[6]+fhatAL[4]*alpha[6])*dfac_v; 
  incr[30] = -0.3061862178478971*(alpha[1]*fhatAL[15]+alpha[0]*fhatAL[14]+alpha[5]*fhatAL[13]+alpha[6]*fhatAL[12]+alpha[8]*fhatAL[11]+fhatAL[8]*alpha[11]+alpha[2]*fhatAL[10]+alpha[3]*fhatAL[9]+alpha[4]*fhatAL[7]+fhatAL[4]*alpha[7])*dfac_v; 
  incr[31] = -0.3061862178478971*(alpha[0]*fhatAL[15]+alpha[1]*fhatAL[14]+alpha[2]*fhatAL[13]+alpha[3]*fhatAL[12]+alpha[4]*fhatAL[11]+fhatAL[4]*alpha[11]+alpha[5]*fhatAL[10]+alpha[6]*fhatAL[9]+alpha[7]*fhatAL[8]+fhatAL[7]*alpha[8])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 
  outr[28] += incr[28]; 
  outr[29] += incr[29]; 
  outr[30] += incr[30]; 
  outr[31] += incr[31]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += incr[19]; 
  outl[20] += -1.0*incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += incr[23]; 
  outl[24] += incr[24]; 
  outl[25] += incr[25]; 
  outl[26] += incr[26]; 
  outl[27] += -1.0*incr[27]; 
  outl[28] += incr[28]; 
  outl[29] += incr[29]; 
  outl[30] += incr[30]; 
  outl[31] += incr[31]; 
  return std::abs(alpha0); 
} 
