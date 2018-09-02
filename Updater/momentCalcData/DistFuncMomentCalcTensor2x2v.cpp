#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc2x2vTensor_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
} 
void MomentCalc2x2vTensor_M0_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
  out[6] += 2.0*f[19]*volFact; 
  out[7] += 2.0*f[20]*volFact; 
  out[8] += 2.0*f[44]*volFact; 
} 
void MomentCalc2x2vTensor_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1)*volFact; 
  out[1] += (2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1)*volFact; 
  out[2] += (2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1)*volFact; 
  out[3] += (2.0*f[5]*wx1+0.5773502691896258*f[11]*dv1)*volFact; 
  out[4] += (2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2)*volFact; 
  out[5] += (2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2)*volFact; 
  out[6] += (2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2)*volFact; 
  out[7] += (2.0*f[5]*wx2+0.5773502691896258*f[12]*dv2)*volFact; 
} 
void MomentCalc2x2vTensor_M1i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1)*volFact; 
  out[1] += (2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1)*volFact; 
  out[2] += (2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1)*volFact; 
  out[3] += (2.0*f[5]*wx1+0.5773502691896258*f[15]*dv1)*volFact; 
  out[4] += (2.0*f[11]*wx1+0.5773502691896257*f[21]*dv1)*volFact; 
  out[5] += (2.0*f[12]*wx1+0.5773502691896257*f[22]*dv1)*volFact; 
  out[6] += (2.0*f[19]*wx1+0.5773502691896257*f[32]*dv1)*volFact; 
  out[7] += (2.0*f[20]*wx1+0.5773502691896257*f[33]*dv1)*volFact; 
  out[8] += (2.0*f[44]*wx1+0.5773502691896258*f[54]*dv1)*volFact; 
  out[9] += (2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2)*volFact; 
  out[10] += (2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2)*volFact; 
  out[11] += (2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2)*volFact; 
  out[12] += (2.0*f[5]*wx2+0.5773502691896258*f[16]*dv2)*volFact; 
  out[13] += (2.0*f[11]*wx2+0.5773502691896257*f[25]*dv2)*volFact; 
  out[14] += (2.0*f[12]*wx2+0.5773502691896257*f[26]*dv2)*volFact; 
  out[15] += (2.0*f[19]*wx2+0.5773502691896257*f[35]*dv2)*volFact; 
  out[16] += (2.0*f[20]*wx2+0.5773502691896257*f[36]*dv2)*volFact; 
  out[17] += (2.0*f[44]*wx2+0.5773502691896258*f[57]*dv2)*volFact; 
} 
void MomentCalc2x2vTensor_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  out[3] += (2.0*f[5]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[5]*dv1_sq)*volFact; 
  out[4] += (2.0*f[0]*wx1*wx2+0.5773502691896258*f[3]*dv1*wx2+0.5773502691896258*f[4]*dv2*wx1+0.1666666666666667*f[10]*dv1*dv2)*volFact; 
  out[5] += (2.0*f[1]*wx1*wx2+0.5773502691896258*f[6]*dv1*wx2+0.5773502691896258*f[8]*dv2*wx1+0.1666666666666667*f[13]*dv1*dv2)*volFact; 
  out[6] += (2.0*f[2]*wx1*wx2+0.5773502691896258*f[7]*dv1*wx2+0.5773502691896258*f[9]*dv2*wx1+0.1666666666666667*f[14]*dv1*dv2)*volFact; 
  out[7] += (2.0*f[5]*wx1*wx2+0.5773502691896258*f[11]*dv1*wx2+0.5773502691896258*f[12]*dv2*wx1+0.1666666666666667*f[15]*dv1*dv2)*volFact; 
  out[8] += (2.0*f[0]*wx2_sq+1.154700538379252*f[4]*dv2*wx2+0.1666666666666667*f[0]*dv2_sq)*volFact; 
  out[9] += (2.0*f[1]*wx2_sq+1.154700538379252*f[8]*dv2*wx2+0.1666666666666667*f[1]*dv2_sq)*volFact; 
  out[10] += (2.0*f[2]*wx2_sq+1.154700538379252*f[9]*dv2*wx2+0.1666666666666667*f[2]*dv2_sq)*volFact; 
  out[11] += (2.0*f[5]*wx2_sq+1.154700538379252*f[12]*dv2*wx2+0.1666666666666667*f[5]*dv2_sq)*volFact; 
} 
void MomentCalc2x2vTensor_M2ij_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[23]*dv1_sq+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[24]*dv1_sq+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  out[3] += (2.0*f[5]*wx1_sq+1.154700538379252*f[15]*dv1*wx1+0.149071198499986*f[34]*dv1_sq+0.1666666666666667*f[5]*dv1_sq)*volFact; 
  out[4] += (2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.149071198499986*f[45]*dv1_sq+0.1666666666666667*f[11]*dv1_sq)*volFact; 
  out[5] += (2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.149071198499986*f[46]*dv1_sq+0.1666666666666667*f[12]*dv1_sq)*volFact; 
  out[6] += (2.0*f[19]*wx1_sq+1.154700538379251*f[32]*dv1*wx1+0.149071198499986*f[55]*dv1_sq+0.1666666666666667*f[19]*dv1_sq)*volFact; 
  out[7] += (2.0*f[20]*wx1_sq+1.154700538379251*f[33]*dv1*wx1+0.149071198499986*f[56]*dv1_sq+0.1666666666666667*f[20]*dv1_sq)*volFact; 
  out[8] += (2.0*f[44]*wx1_sq+1.154700538379252*f[54]*dv1*wx1+0.149071198499986*f[72]*dv1_sq+0.1666666666666667*f[44]*dv1_sq)*volFact; 
  out[9] += (2.0*f[0]*wx1*wx2+0.5773502691896258*f[3]*dv1*wx2+0.5773502691896258*f[4]*dv2*wx1+0.1666666666666667*f[10]*dv1*dv2)*volFact; 
  out[10] += (2.0*f[1]*wx1*wx2+0.5773502691896258*f[6]*dv1*wx2+0.5773502691896258*f[8]*dv2*wx1+0.1666666666666667*f[17]*dv1*dv2)*volFact; 
  out[11] += (2.0*f[2]*wx1*wx2+0.5773502691896258*f[7]*dv1*wx2+0.5773502691896258*f[9]*dv2*wx1+0.1666666666666667*f[18]*dv1*dv2)*volFact; 
  out[12] += (2.0*f[5]*wx1*wx2+0.5773502691896258*f[15]*dv1*wx2+0.5773502691896258*f[16]*dv2*wx1+0.1666666666666667*f[31]*dv1*dv2)*volFact; 
  out[13] += (2.0*f[11]*wx1*wx2+0.5773502691896257*f[21]*dv1*wx2+0.5773502691896257*f[25]*dv2*wx1+0.1666666666666667*f[37]*dv1*dv2)*volFact; 
  out[14] += (2.0*f[12]*wx1*wx2+0.5773502691896257*f[22]*dv1*wx2+0.5773502691896257*f[26]*dv2*wx1+0.1666666666666667*f[38]*dv1*dv2)*volFact; 
  out[15] += (2.0*f[19]*wx1*wx2+0.5773502691896257*f[32]*dv1*wx2+0.5773502691896257*f[35]*dv2*wx1+0.1666666666666667*f[50]*dv1*dv2)*volFact; 
  out[16] += (2.0*f[20]*wx1*wx2+0.5773502691896257*f[33]*dv1*wx2+0.5773502691896257*f[36]*dv2*wx1+0.1666666666666667*f[51]*dv1*dv2)*volFact; 
  out[17] += (2.0*f[44]*wx1*wx2+0.5773502691896258*f[54]*dv1*wx2+0.5773502691896258*f[57]*dv2*wx1+0.1666666666666667*f[66]*dv1*dv2)*volFact; 
  out[18] += (2.0*f[0]*wx2_sq+1.154700538379252*f[4]*dv2*wx2+0.149071198499986*f[14]*dv2_sq+0.1666666666666667*f[0]*dv2_sq)*volFact; 
  out[19] += (2.0*f[1]*wx2_sq+1.154700538379252*f[8]*dv2*wx2+0.149071198499986*f[28]*dv2_sq+0.1666666666666667*f[1]*dv2_sq)*volFact; 
  out[20] += (2.0*f[2]*wx2_sq+1.154700538379252*f[9]*dv2*wx2+0.149071198499986*f[29]*dv2_sq+0.1666666666666667*f[2]*dv2_sq)*volFact; 
  out[21] += (2.0*f[5]*wx2_sq+1.154700538379252*f[16]*dv2*wx2+0.149071198499986*f[41]*dv2_sq+0.1666666666666667*f[5]*dv2_sq)*volFact; 
  out[22] += (2.0*f[11]*wx2_sq+1.154700538379251*f[25]*dv2*wx2+0.149071198499986*f[47]*dv2_sq+0.1666666666666667*f[11]*dv2_sq)*volFact; 
  out[23] += (2.0*f[12]*wx2_sq+1.154700538379251*f[26]*dv2*wx2+0.149071198499986*f[48]*dv2_sq+0.1666666666666667*f[12]*dv2_sq)*volFact; 
  out[24] += (2.0*f[19]*wx2_sq+1.154700538379251*f[35]*dv2*wx2+0.149071198499986*f[60]*dv2_sq+0.1666666666666667*f[19]*dv2_sq)*volFact; 
  out[25] += (2.0*f[20]*wx2_sq+1.154700538379251*f[36]*dv2*wx2+0.149071198499986*f[61]*dv2_sq+0.1666666666666667*f[20]*dv2_sq)*volFact; 
  out[26] += (2.0*f[44]*wx2_sq+1.154700538379252*f[57]*dv2*wx2+0.149071198499986*f[73]*dv2_sq+0.1666666666666667*f[44]*dv2_sq)*volFact; 
} 
void MomentCalc2x2vTensor_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx2_sq+1.154700538379252*f[4]*dv2*wx2+2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.1666666666666667*f[0]*dv2_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx2_sq+1.154700538379252*f[8]*dv2*wx2+2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.1666666666666667*f[1]*dv2_sq+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx2_sq+1.154700538379252*f[9]*dv2*wx2+2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.1666666666666667*f[2]*dv2_sq+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  out[3] += (2.0*f[5]*wx2_sq+1.154700538379252*f[12]*dv2*wx2+2.0*f[5]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[5]*dv2_sq+0.1666666666666667*f[5]*dv1_sq)*volFact; 
} 
void MomentCalc2x2vTensor_M2_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx2_sq+1.154700538379252*f[4]*dv2*wx2+2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[14]*dv2_sq+0.1666666666666667*f[0]*dv2_sq+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx2_sq+1.154700538379252*f[8]*dv2*wx2+2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[28]*dv2_sq+0.1666666666666667*f[1]*dv2_sq+0.149071198499986*f[23]*dv1_sq+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx2_sq+1.154700538379252*f[9]*dv2*wx2+2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[29]*dv2_sq+0.1666666666666667*f[2]*dv2_sq+0.149071198499986*f[24]*dv1_sq+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  out[3] += (2.0*f[5]*wx2_sq+1.154700538379252*f[16]*dv2*wx2+2.0*f[5]*wx1_sq+1.154700538379252*f[15]*dv1*wx1+0.149071198499986*f[41]*dv2_sq+0.1666666666666667*f[5]*dv2_sq+0.149071198499986*f[34]*dv1_sq+0.1666666666666667*f[5]*dv1_sq)*volFact; 
  out[4] += (2.0*f[11]*wx2_sq+1.154700538379251*f[25]*dv2*wx2+2.0*f[11]*wx1_sq+1.154700538379251*f[21]*dv1*wx1+0.149071198499986*f[47]*dv2_sq+0.1666666666666667*f[11]*dv2_sq+0.149071198499986*f[45]*dv1_sq+0.1666666666666667*f[11]*dv1_sq)*volFact; 
  out[5] += (2.0*f[12]*wx2_sq+1.154700538379251*f[26]*dv2*wx2+2.0*f[12]*wx1_sq+1.154700538379251*f[22]*dv1*wx1+0.149071198499986*f[48]*dv2_sq+0.1666666666666667*f[12]*dv2_sq+0.149071198499986*f[46]*dv1_sq+0.1666666666666667*f[12]*dv1_sq)*volFact; 
  out[6] += (2.0*f[19]*wx2_sq+1.154700538379251*f[35]*dv2*wx2+2.0*f[19]*wx1_sq+1.154700538379251*f[32]*dv1*wx1+0.149071198499986*f[60]*dv2_sq+0.1666666666666667*f[19]*dv2_sq+0.149071198499986*f[55]*dv1_sq+0.1666666666666667*f[19]*dv1_sq)*volFact; 
  out[7] += (2.0*f[20]*wx2_sq+1.154700538379251*f[36]*dv2*wx2+2.0*f[20]*wx1_sq+1.154700538379251*f[33]*dv1*wx1+0.149071198499986*f[61]*dv2_sq+0.1666666666666667*f[20]*dv2_sq+0.149071198499986*f[56]*dv1_sq+0.1666666666666667*f[20]*dv1_sq)*volFact; 
  out[8] += (2.0*f[44]*wx2_sq+1.154700538379252*f[57]*dv2*wx2+2.0*f[44]*wx1_sq+1.154700538379252*f[54]*dv1*wx1+0.149071198499986*f[73]*dv2_sq+0.1666666666666667*f[44]*dv2_sq+0.149071198499986*f[72]*dv1_sq+0.1666666666666667*f[44]*dv1_sq)*volFact; 
} 
void MomentCalc2x2vTensor_M3i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += (2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[3]*dv1*wx2_sq+1.154700538379252*f[4]*dv2*wx1*wx2+0.3333333333333333*f[10]*dv1*dv2*wx2+2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[3]*dv1*wx1_sq+0.1666666666666667*f[0]*dv2_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.04811252243246882*f[3]*dv1*dv2_sq+0.08660254037844387*f[3]*dv1*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[6]*dv1*wx2_sq+1.154700538379252*f[8]*dv2*wx1*wx2+0.3333333333333333*f[13]*dv1*dv2*wx2+2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[6]*dv1*wx1_sq+0.1666666666666667*f[1]*dv2_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.04811252243246882*f[6]*dv1*dv2_sq+0.08660254037844387*f[6]*dv1*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1*wx2_sq+0.5773502691896258*f[7]*dv1*wx2_sq+1.154700538379252*f[9]*dv2*wx1*wx2+0.3333333333333333*f[14]*dv1*dv2*wx2+2.0*f[2]*wx1*wx1_sq+1.732050807568877*f[7]*dv1*wx1_sq+0.1666666666666667*f[2]*dv2_sq*wx1+0.5*f[2]*dv1_sq*wx1+0.04811252243246882*f[7]*dv1*dv2_sq+0.08660254037844387*f[7]*dv1*dv1_sq)*volFact; 
  out[3] += (2.0*f[5]*wx1*wx2_sq+0.5773502691896258*f[11]*dv1*wx2_sq+1.154700538379252*f[12]*dv2*wx1*wx2+0.3333333333333333*f[15]*dv1*dv2*wx2+2.0*f[5]*wx1*wx1_sq+1.732050807568877*f[11]*dv1*wx1_sq+0.1666666666666667*f[5]*dv2_sq*wx1+0.5*f[5]*dv1_sq*wx1+0.04811252243246882*f[11]*dv1*dv2_sq+0.08660254037844387*f[11]*dv1*dv1_sq)*volFact; 
  out[4] += (2.0*f[0]*wx2*wx2_sq+1.732050807568877*f[4]*dv2*wx2_sq+2.0*f[0]*wx1_sq*wx2+1.154700538379252*f[3]*dv1*wx1*wx2+0.5*f[0]*dv2_sq*wx2+0.1666666666666667*f[0]*dv1_sq*wx2+0.5773502691896258*f[4]*dv2*wx1_sq+0.3333333333333333*f[10]*dv1*dv2*wx1+0.08660254037844387*f[4]*dv2*dv2_sq+0.04811252243246882*f[4]*dv1_sq*dv2)*volFact; 
  out[5] += (2.0*f[1]*wx2*wx2_sq+1.732050807568877*f[8]*dv2*wx2_sq+2.0*f[1]*wx1_sq*wx2+1.154700538379252*f[6]*dv1*wx1*wx2+0.5*f[1]*dv2_sq*wx2+0.1666666666666667*f[1]*dv1_sq*wx2+0.5773502691896258*f[8]*dv2*wx1_sq+0.3333333333333333*f[13]*dv1*dv2*wx1+0.08660254037844387*f[8]*dv2*dv2_sq+0.04811252243246882*f[8]*dv1_sq*dv2)*volFact; 
  out[6] += (2.0*f[2]*wx2*wx2_sq+1.732050807568877*f[9]*dv2*wx2_sq+2.0*f[2]*wx1_sq*wx2+1.154700538379252*f[7]*dv1*wx1*wx2+0.5*f[2]*dv2_sq*wx2+0.1666666666666667*f[2]*dv1_sq*wx2+0.5773502691896258*f[9]*dv2*wx1_sq+0.3333333333333333*f[14]*dv1*dv2*wx1+0.08660254037844387*f[9]*dv2*dv2_sq+0.04811252243246882*f[9]*dv1_sq*dv2)*volFact; 
  out[7] += (2.0*f[5]*wx2*wx2_sq+1.732050807568877*f[12]*dv2*wx2_sq+2.0*f[5]*wx1_sq*wx2+1.154700538379252*f[11]*dv1*wx1*wx2+0.5*f[5]*dv2_sq*wx2+0.1666666666666667*f[5]*dv1_sq*wx2+0.5773502691896258*f[12]*dv2*wx1_sq+0.3333333333333333*f[15]*dv1*dv2*wx1+0.08660254037844387*f[12]*dv2*dv2_sq+0.04811252243246882*f[12]*dv1_sq*dv2)*volFact; 
} 
void MomentCalc2x2vTensor_M3i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += (2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[3]*dv1*wx2_sq+1.154700538379252*f[4]*dv2*wx1*wx2+0.3333333333333333*f[10]*dv1*dv2*wx2+2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[3]*dv1*wx1_sq+0.149071198499986*f[14]*dv2_sq*wx1+0.1666666666666667*f[0]*dv2_sq*wx1+0.4472135954999579*f[13]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.04303314829119351*f[30]*dv1*dv2_sq+0.04811252243246882*f[3]*dv1*dv2_sq+0.08660254037844387*f[3]*dv1*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[6]*dv1*wx2_sq+1.154700538379252*f[8]*dv2*wx1*wx2+0.3333333333333333*f[17]*dv1*dv2*wx2+2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[6]*dv1*wx1_sq+0.149071198499986*f[28]*dv2_sq*wx1+0.1666666666666667*f[1]*dv2_sq*wx1+0.447213595499958*f[23]*dv1_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.04303314829119353*f[42]*dv1*dv2_sq+0.04811252243246882*f[6]*dv1*dv2_sq+0.08660254037844387*f[6]*dv1*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1*wx2_sq+0.5773502691896258*f[7]*dv1*wx2_sq+1.154700538379252*f[9]*dv2*wx1*wx2+0.3333333333333333*f[18]*dv1*dv2*wx2+2.0*f[2]*wx1*wx1_sq+1.732050807568877*f[7]*dv1*wx1_sq+0.149071198499986*f[29]*dv2_sq*wx1+0.1666666666666667*f[2]*dv2_sq*wx1+0.447213595499958*f[24]*dv1_sq*wx1+0.5*f[2]*dv1_sq*wx1+0.04303314829119353*f[43]*dv1*dv2_sq+0.04811252243246882*f[7]*dv1*dv2_sq+0.08660254037844387*f[7]*dv1*dv1_sq)*volFact; 
  out[3] += (2.0*f[5]*wx1*wx2_sq+0.5773502691896258*f[15]*dv1*wx2_sq+1.154700538379252*f[16]*dv2*wx1*wx2+0.3333333333333333*f[31]*dv1*dv2*wx2+2.0*f[5]*wx1*wx1_sq+1.732050807568877*f[15]*dv1*wx1_sq+0.149071198499986*f[41]*dv2_sq*wx1+0.1666666666666667*f[5]*dv2_sq*wx1+0.4472135954999579*f[34]*dv1_sq*wx1+0.5*f[5]*dv1_sq*wx1+0.04303314829119351*f[53]*dv1*dv2_sq+0.04811252243246882*f[15]*dv1*dv2_sq+0.08660254037844387*f[15]*dv1*dv1_sq)*volFact; 
  out[4] += (2.0*f[11]*wx1*wx2_sq+0.5773502691896257*f[21]*dv1*wx2_sq+1.154700538379251*f[25]*dv2*wx1*wx2+0.3333333333333333*f[37]*dv1*dv2*wx2+2.0*f[11]*wx1*wx1_sq+1.732050807568877*f[21]*dv1*wx1_sq+0.149071198499986*f[47]*dv2_sq*wx1+0.1666666666666667*f[11]*dv2_sq*wx1+0.4472135954999579*f[45]*dv1_sq*wx1+0.5*f[11]*dv1_sq*wx1+0.04303314829119353*f[62]*dv1*dv2_sq+0.04811252243246881*f[21]*dv1*dv2_sq+0.08660254037844385*f[21]*dv1*dv1_sq)*volFact; 
  out[5] += (2.0*f[12]*wx1*wx2_sq+0.5773502691896257*f[22]*dv1*wx2_sq+1.154700538379251*f[26]*dv2*wx1*wx2+0.3333333333333333*f[38]*dv1*dv2*wx2+2.0*f[12]*wx1*wx1_sq+1.732050807568877*f[22]*dv1*wx1_sq+0.149071198499986*f[48]*dv2_sq*wx1+0.1666666666666667*f[12]*dv2_sq*wx1+0.4472135954999579*f[46]*dv1_sq*wx1+0.5*f[12]*dv1_sq*wx1+0.04303314829119353*f[63]*dv1*dv2_sq+0.04811252243246881*f[22]*dv1*dv2_sq+0.08660254037844385*f[22]*dv1*dv1_sq)*volFact; 
  out[6] += (2.0*f[19]*wx1*wx2_sq+0.5773502691896257*f[32]*dv1*wx2_sq+1.154700538379251*f[35]*dv2*wx1*wx2+0.3333333333333333*f[50]*dv1*dv2*wx2+2.0*f[19]*wx1*wx1_sq+1.732050807568877*f[32]*dv1*wx1_sq+0.149071198499986*f[60]*dv2_sq*wx1+0.1666666666666667*f[19]*dv2_sq*wx1+0.447213595499958*f[55]*dv1_sq*wx1+0.5*f[19]*dv1_sq*wx1+0.04303314829119351*f[69]*dv1*dv2_sq+0.04811252243246881*f[32]*dv1*dv2_sq+0.08660254037844385*f[32]*dv1*dv1_sq)*volFact; 
  out[7] += (2.0*f[20]*wx1*wx2_sq+0.5773502691896257*f[33]*dv1*wx2_sq+1.154700538379251*f[36]*dv2*wx1*wx2+0.3333333333333333*f[51]*dv1*dv2*wx2+2.0*f[20]*wx1*wx1_sq+1.732050807568877*f[33]*dv1*wx1_sq+0.149071198499986*f[61]*dv2_sq*wx1+0.1666666666666667*f[20]*dv2_sq*wx1+0.447213595499958*f[56]*dv1_sq*wx1+0.5*f[20]*dv1_sq*wx1+0.04303314829119351*f[70]*dv1*dv2_sq+0.04811252243246881*f[33]*dv1*dv2_sq+0.08660254037844385*f[33]*dv1*dv1_sq)*volFact; 
  out[8] += (2.0*f[44]*wx1*wx2_sq+0.5773502691896258*f[54]*dv1*wx2_sq+1.154700538379252*f[57]*dv2*wx1*wx2+0.3333333333333333*f[66]*dv1*dv2*wx2+2.0*f[44]*wx1*wx1_sq+1.732050807568877*f[54]*dv1*wx1_sq+0.149071198499986*f[73]*dv2_sq*wx1+0.1666666666666667*f[44]*dv2_sq*wx1+0.4472135954999579*f[72]*dv1_sq*wx1+0.5*f[44]*dv1_sq*wx1+0.04303314829119351*f[77]*dv1*dv2_sq+0.04811252243246882*f[54]*dv1*dv2_sq+0.08660254037844387*f[54]*dv1*dv1_sq)*volFact; 
  out[9] += (2.0*f[0]*wx2*wx2_sq+1.732050807568877*f[4]*dv2*wx2_sq+2.0*f[0]*wx1_sq*wx2+1.154700538379252*f[3]*dv1*wx1*wx2+0.4472135954999579*f[14]*dv2_sq*wx2+0.5*f[0]*dv2_sq*wx2+0.149071198499986*f[13]*dv1_sq*wx2+0.1666666666666667*f[0]*dv1_sq*wx2+0.5773502691896258*f[4]*dv2*wx1_sq+0.3333333333333333*f[10]*dv1*dv2*wx1+0.08660254037844387*f[4]*dv2*dv2_sq+0.04303314829119351*f[27]*dv1_sq*dv2+0.04811252243246882*f[4]*dv1_sq*dv2)*volFact; 
  out[10] += (2.0*f[1]*wx2*wx2_sq+1.732050807568877*f[8]*dv2*wx2_sq+2.0*f[1]*wx1_sq*wx2+1.154700538379252*f[6]*dv1*wx1*wx2+0.447213595499958*f[28]*dv2_sq*wx2+0.5*f[1]*dv2_sq*wx2+0.149071198499986*f[23]*dv1_sq*wx2+0.1666666666666667*f[1]*dv1_sq*wx2+0.5773502691896258*f[8]*dv2*wx1_sq+0.3333333333333333*f[17]*dv1*dv2*wx1+0.08660254037844387*f[8]*dv2*dv2_sq+0.04303314829119353*f[39]*dv1_sq*dv2+0.04811252243246882*f[8]*dv1_sq*dv2)*volFact; 
  out[11] += (2.0*f[2]*wx2*wx2_sq+1.732050807568877*f[9]*dv2*wx2_sq+2.0*f[2]*wx1_sq*wx2+1.154700538379252*f[7]*dv1*wx1*wx2+0.447213595499958*f[29]*dv2_sq*wx2+0.5*f[2]*dv2_sq*wx2+0.149071198499986*f[24]*dv1_sq*wx2+0.1666666666666667*f[2]*dv1_sq*wx2+0.5773502691896258*f[9]*dv2*wx1_sq+0.3333333333333333*f[18]*dv1*dv2*wx1+0.08660254037844387*f[9]*dv2*dv2_sq+0.04303314829119353*f[40]*dv1_sq*dv2+0.04811252243246882*f[9]*dv1_sq*dv2)*volFact; 
  out[12] += (2.0*f[5]*wx2*wx2_sq+1.732050807568877*f[16]*dv2*wx2_sq+2.0*f[5]*wx1_sq*wx2+1.154700538379252*f[15]*dv1*wx1*wx2+0.4472135954999579*f[41]*dv2_sq*wx2+0.5*f[5]*dv2_sq*wx2+0.149071198499986*f[34]*dv1_sq*wx2+0.1666666666666667*f[5]*dv1_sq*wx2+0.5773502691896258*f[16]*dv2*wx1_sq+0.3333333333333333*f[31]*dv1*dv2*wx1+0.08660254037844387*f[16]*dv2*dv2_sq+0.04303314829119351*f[52]*dv1_sq*dv2+0.04811252243246882*f[16]*dv1_sq*dv2)*volFact; 
  out[13] += (2.0*f[11]*wx2*wx2_sq+1.732050807568877*f[25]*dv2*wx2_sq+2.0*f[11]*wx1_sq*wx2+1.154700538379251*f[21]*dv1*wx1*wx2+0.4472135954999579*f[47]*dv2_sq*wx2+0.5*f[11]*dv2_sq*wx2+0.149071198499986*f[45]*dv1_sq*wx2+0.1666666666666667*f[11]*dv1_sq*wx2+0.5773502691896257*f[25]*dv2*wx1_sq+0.3333333333333333*f[37]*dv1*dv2*wx1+0.08660254037844385*f[25]*dv2*dv2_sq+0.04303314829119353*f[58]*dv1_sq*dv2+0.04811252243246881*f[25]*dv1_sq*dv2)*volFact; 
  out[14] += (2.0*f[12]*wx2*wx2_sq+1.732050807568877*f[26]*dv2*wx2_sq+2.0*f[12]*wx1_sq*wx2+1.154700538379251*f[22]*dv1*wx1*wx2+0.4472135954999579*f[48]*dv2_sq*wx2+0.5*f[12]*dv2_sq*wx2+0.149071198499986*f[46]*dv1_sq*wx2+0.1666666666666667*f[12]*dv1_sq*wx2+0.5773502691896257*f[26]*dv2*wx1_sq+0.3333333333333333*f[38]*dv1*dv2*wx1+0.08660254037844385*f[26]*dv2*dv2_sq+0.04303314829119353*f[59]*dv1_sq*dv2+0.04811252243246881*f[26]*dv1_sq*dv2)*volFact; 
  out[15] += (2.0*f[19]*wx2*wx2_sq+1.732050807568877*f[35]*dv2*wx2_sq+2.0*f[19]*wx1_sq*wx2+1.154700538379251*f[32]*dv1*wx1*wx2+0.447213595499958*f[60]*dv2_sq*wx2+0.5*f[19]*dv2_sq*wx2+0.149071198499986*f[55]*dv1_sq*wx2+0.1666666666666667*f[19]*dv1_sq*wx2+0.5773502691896257*f[35]*dv2*wx1_sq+0.3333333333333333*f[50]*dv1*dv2*wx1+0.08660254037844385*f[35]*dv2*dv2_sq+0.04303314829119351*f[67]*dv1_sq*dv2+0.04811252243246881*f[35]*dv1_sq*dv2)*volFact; 
  out[16] += (2.0*f[20]*wx2*wx2_sq+1.732050807568877*f[36]*dv2*wx2_sq+2.0*f[20]*wx1_sq*wx2+1.154700538379251*f[33]*dv1*wx1*wx2+0.447213595499958*f[61]*dv2_sq*wx2+0.5*f[20]*dv2_sq*wx2+0.149071198499986*f[56]*dv1_sq*wx2+0.1666666666666667*f[20]*dv1_sq*wx2+0.5773502691896257*f[36]*dv2*wx1_sq+0.3333333333333333*f[51]*dv1*dv2*wx1+0.08660254037844385*f[36]*dv2*dv2_sq+0.04303314829119351*f[68]*dv1_sq*dv2+0.04811252243246881*f[36]*dv1_sq*dv2)*volFact; 
  out[17] += (2.0*f[44]*wx2*wx2_sq+1.732050807568877*f[57]*dv2*wx2_sq+2.0*f[44]*wx1_sq*wx2+1.154700538379252*f[54]*dv1*wx1*wx2+0.4472135954999579*f[73]*dv2_sq*wx2+0.5*f[44]*dv2_sq*wx2+0.149071198499986*f[72]*dv1_sq*wx2+0.1666666666666667*f[44]*dv1_sq*wx2+0.5773502691896258*f[57]*dv2*wx1_sq+0.3333333333333333*f[66]*dv1*dv2*wx1+0.08660254037844387*f[57]*dv2*dv2_sq+0.04303314829119351*f[76]*dv1_sq*dv2+0.04811252243246882*f[57]*dv1_sq*dv2)*volFact; 
} 
void MomentCalc2x2vTensor_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double tempM0[4], tempM1i[8]; 

  tempM0[0] = 2.0*f[0]*volFact; 
  tempM0[1] = 2.0*f[1]*volFact; 
  tempM0[2] = 2.0*f[2]*volFact; 
  tempM0[3] = 2.0*f[5]*volFact; 

  tempM1i[0] = 0.5773502691896258*f[3]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.5773502691896258*f[6]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.5773502691896258*f[7]*dv1*volFact+tempM0[2]*wx1; 
  tempM1i[3] = 0.5773502691896258*f[11]*dv1*volFact+tempM0[3]*wx1; 
  tempM1i[4] = 0.5773502691896258*f[4]*dv2*volFact+tempM0[0]*wx2; 
  tempM1i[5] = 0.5773502691896258*f[8]*dv2*volFact+tempM0[1]*wx2; 
  tempM1i[6] = 0.5773502691896258*f[9]*dv2*volFact+tempM0[2]*wx2; 
  tempM1i[7] = 0.5773502691896258*f[12]*dv2*volFact+tempM0[3]*wx2; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM0[2] += tempM0[2]; 
  outM0[3] += tempM0[3]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM1i[2] += tempM1i[2]; 
  outM1i[3] += tempM1i[3]; 
  outM1i[4] += tempM1i[4]; 
  outM1i[5] += tempM1i[5]; 
  outM1i[6] += tempM1i[6]; 
  outM1i[7] += tempM1i[7]; 
  outM2[0] += (0.1666666666666667*f[0]*dv2_sq+0.1666666666666667*f[0]*dv1_sq)*volFact+tempM0[0]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[4]*wx2+2.0*tempM1i[0]*wx1; 
  outM2[1] += (0.1666666666666667*f[1]*dv2_sq+0.1666666666666667*f[1]*dv1_sq)*volFact+tempM0[1]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[5]*wx2+2.0*tempM1i[1]*wx1; 
  outM2[2] += (0.1666666666666667*f[2]*dv2_sq+0.1666666666666667*f[2]*dv1_sq)*volFact+tempM0[2]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[6]*wx2+2.0*tempM1i[2]*wx1; 
  outM2[3] += (0.1666666666666667*f[5]*dv2_sq+0.1666666666666667*f[5]*dv1_sq)*volFact+tempM0[3]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[7]*wx2+2.0*tempM1i[3]*wx1; 
} 
void MomentCalc2x2vTensor_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double tempM0[9], tempM1i[18]; 

  tempM0[0] = 2.0*f[0]*volFact; 
  tempM0[1] = 2.0*f[1]*volFact; 
  tempM0[2] = 2.0*f[2]*volFact; 
  tempM0[3] = 2.0*f[5]*volFact; 
  tempM0[4] = 2.0*f[11]*volFact; 
  tempM0[5] = 2.0*f[12]*volFact; 
  tempM0[6] = 2.0*f[19]*volFact; 
  tempM0[7] = 2.0*f[20]*volFact; 
  tempM0[8] = 2.0*f[44]*volFact; 

  tempM1i[0] = 0.5773502691896258*f[3]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.5773502691896258*f[6]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.5773502691896258*f[7]*dv1*volFact+tempM0[2]*wx1; 
  tempM1i[3] = 0.5773502691896258*f[15]*dv1*volFact+tempM0[3]*wx1; 
  tempM1i[4] = 0.5773502691896257*f[21]*dv1*volFact+tempM0[4]*wx1; 
  tempM1i[5] = 0.5773502691896257*f[22]*dv1*volFact+tempM0[5]*wx1; 
  tempM1i[6] = 0.5773502691896257*f[32]*dv1*volFact+tempM0[6]*wx1; 
  tempM1i[7] = 0.5773502691896257*f[33]*dv1*volFact+tempM0[7]*wx1; 
  tempM1i[8] = 0.5773502691896258*f[54]*dv1*volFact+tempM0[8]*wx1; 
  tempM1i[9] = 0.5773502691896258*f[4]*dv2*volFact+tempM0[0]*wx2; 
  tempM1i[10] = 0.5773502691896258*f[8]*dv2*volFact+tempM0[1]*wx2; 
  tempM1i[11] = 0.5773502691896258*f[9]*dv2*volFact+tempM0[2]*wx2; 
  tempM1i[12] = 0.5773502691896258*f[16]*dv2*volFact+tempM0[3]*wx2; 
  tempM1i[13] = 0.5773502691896257*f[25]*dv2*volFact+tempM0[4]*wx2; 
  tempM1i[14] = 0.5773502691896257*f[26]*dv2*volFact+tempM0[5]*wx2; 
  tempM1i[15] = 0.5773502691896257*f[35]*dv2*volFact+tempM0[6]*wx2; 
  tempM1i[16] = 0.5773502691896257*f[36]*dv2*volFact+tempM0[7]*wx2; 
  tempM1i[17] = 0.5773502691896258*f[57]*dv2*volFact+tempM0[8]*wx2; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM0[2] += tempM0[2]; 
  outM0[3] += tempM0[3]; 
  outM0[4] += tempM0[4]; 
  outM0[5] += tempM0[5]; 
  outM0[6] += tempM0[6]; 
  outM0[7] += tempM0[7]; 
  outM0[8] += tempM0[8]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM1i[2] += tempM1i[2]; 
  outM1i[3] += tempM1i[3]; 
  outM1i[4] += tempM1i[4]; 
  outM1i[5] += tempM1i[5]; 
  outM1i[6] += tempM1i[6]; 
  outM1i[7] += tempM1i[7]; 
  outM1i[8] += tempM1i[8]; 
  outM1i[9] += tempM1i[9]; 
  outM1i[10] += tempM1i[10]; 
  outM1i[11] += tempM1i[11]; 
  outM1i[12] += tempM1i[12]; 
  outM1i[13] += tempM1i[13]; 
  outM1i[14] += tempM1i[14]; 
  outM1i[15] += tempM1i[15]; 
  outM1i[16] += tempM1i[16]; 
  outM1i[17] += tempM1i[17]; 
  outM2[0] += (0.149071198499986*f[14]*dv2_sq+0.1666666666666667*f[0]*dv2_sq+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact+tempM0[0]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[9]*wx2+2.0*tempM1i[0]*wx1; 
  outM2[1] += (0.149071198499986*f[28]*dv2_sq+0.1666666666666667*f[1]*dv2_sq+0.149071198499986*f[23]*dv1_sq+0.1666666666666667*f[1]*dv1_sq)*volFact+tempM0[1]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[10]*wx2+2.0*tempM1i[1]*wx1; 
  outM2[2] += (0.149071198499986*f[29]*dv2_sq+0.1666666666666667*f[2]*dv2_sq+0.149071198499986*f[24]*dv1_sq+0.1666666666666667*f[2]*dv1_sq)*volFact+tempM0[2]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[11]*wx2+2.0*tempM1i[2]*wx1; 
  outM2[3] += (0.149071198499986*f[41]*dv2_sq+0.1666666666666667*f[5]*dv2_sq+0.149071198499986*f[34]*dv1_sq+0.1666666666666667*f[5]*dv1_sq)*volFact+tempM0[3]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[12]*wx2+2.0*tempM1i[3]*wx1; 
  outM2[4] += (0.149071198499986*f[47]*dv2_sq+0.1666666666666667*f[11]*dv2_sq+0.149071198499986*f[45]*dv1_sq+0.1666666666666667*f[11]*dv1_sq)*volFact+tempM0[4]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[13]*wx2+2.0*tempM1i[4]*wx1; 
  outM2[5] += (0.149071198499986*f[48]*dv2_sq+0.1666666666666667*f[12]*dv2_sq+0.149071198499986*f[46]*dv1_sq+0.1666666666666667*f[12]*dv1_sq)*volFact+tempM0[5]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[14]*wx2+2.0*tempM1i[5]*wx1; 
  outM2[6] += (0.149071198499986*f[60]*dv2_sq+0.1666666666666667*f[19]*dv2_sq+0.149071198499986*f[55]*dv1_sq+0.1666666666666667*f[19]*dv1_sq)*volFact+tempM0[6]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[15]*wx2+2.0*tempM1i[6]*wx1; 
  outM2[7] += (0.149071198499986*f[61]*dv2_sq+0.1666666666666667*f[20]*dv2_sq+0.149071198499986*f[56]*dv1_sq+0.1666666666666667*f[20]*dv1_sq)*volFact+tempM0[7]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[16]*wx2+2.0*tempM1i[7]*wx1; 
  outM2[8] += (0.149071198499986*f[73]*dv2_sq+0.1666666666666667*f[44]*dv2_sq+0.149071198499986*f[72]*dv1_sq+0.1666666666666667*f[44]*dv1_sq)*volFact+tempM0[8]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[17]*wx2+2.0*tempM1i[8]*wx1; 
} 
void MomentCalc2x2vTensor_StarMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double tempM0[4], tempM1i[8]; 

  tempM0[0] = 2.0*f[0]*volFact; 
  tempM0[1] = 2.0*f[1]*volFact; 
  tempM0[2] = 2.0*f[2]*volFact; 
  tempM0[3] = 2.0*f[5]*volFact; 

  tempM1i[0] = 0.5773502691896258*f[3]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.5773502691896258*f[6]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.5773502691896258*f[7]*dv1*volFact+tempM0[2]*wx1; 
  tempM1i[3] = 0.5773502691896258*f[11]*dv1*volFact+tempM0[3]*wx1; 
  tempM1i[4] = 0.5773502691896258*f[4]*dv2*volFact+tempM0[0]*wx2; 
  tempM1i[5] = 0.5773502691896258*f[8]*dv2*volFact+tempM0[1]*wx2; 
  tempM1i[6] = 0.5773502691896258*f[9]*dv2*volFact+tempM0[2]*wx2; 
  tempM1i[7] = 0.5773502691896258*f[12]*dv2*volFact+tempM0[3]*wx2; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM0[2] += tempM0[2]; 
  outM0[3] += tempM0[3]; 
  outM1i[0] += tempM0[0]*wx1; 
  outM1i[1] += tempM0[1]*wx1; 
  outM1i[2] += tempM0[2]*wx1; 
  outM1i[3] += tempM0[3]*wx1; 
  outM1i[4] += tempM0[0]*wx2; 
  outM1i[5] += tempM0[1]*wx2; 
  outM1i[6] += tempM0[2]*wx2; 
  outM1i[7] += tempM0[3]*wx2; 
  outM2[0] += tempM0[0]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[4]*wx2+2.0*tempM1i[0]*wx1; 
  outM2[1] += tempM0[1]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[5]*wx2+2.0*tempM1i[1]*wx1; 
  outM2[2] += tempM0[2]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[6]*wx2+2.0*tempM1i[2]*wx1; 
  outM2[3] += tempM0[3]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[7]*wx2+2.0*tempM1i[3]*wx1; 
} 
