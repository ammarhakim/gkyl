#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc2x2vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
} 
void MomentCalc2x2vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
} 
void MomentCalc2x2vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out) 
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
  out[8] += 2.0*f[31]*volFact; 
  out[9] += 2.0*f[32]*volFact; 
} 
void MomentCalc2x2vMax_M0_P4(const double *w, const double *dxv, const double *f, double *out) 
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
  out[8] += 2.0*f[31]*volFact; 
  out[9] += 2.0*f[32]*volFact; 
  out[10] += 2.0*f[48]*volFact; 
  out[11] += 2.0*f[54]*volFact; 
  out[12] += 2.0*f[55]*volFact; 
  out[13] += 2.0*f[66]*volFact; 
  out[14] += 2.0*f[67]*volFact; 
} 
void MomentCalc2x2vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += 2.0*f[0]*wx1*volFact+0.5773502691896258*f[3]*dv1*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact; 
  out[2] += 2.0*f[2]*wx1*volFact; 
  out[3] += 2.0*f[0]*wx2*volFact+0.5773502691896258*f[4]*dv2*volFact; 
  out[4] += 2.0*f[1]*wx2*volFact; 
  out[5] += 2.0*f[2]*wx2*volFact; 
} 
void MomentCalc2x2vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += 2.0*f[0]*wx1*volFact+0.5773502691896258*f[3]*dv1*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact+0.5773502691896258*f[6]*dv1*volFact; 
  out[2] += 2.0*f[2]*wx1*volFact+0.5773502691896258*f[7]*dv1*volFact; 
  out[3] += 2.0*f[5]*wx1*volFact; 
  out[4] += 2.0*f[11]*wx1*volFact; 
  out[5] += 2.0*f[12]*wx1*volFact; 
  out[6] += 2.0*f[0]*wx2*volFact+0.5773502691896258*f[4]*dv2*volFact; 
  out[7] += 2.0*f[1]*wx2*volFact+0.5773502691896258*f[8]*dv2*volFact; 
  out[8] += 2.0*f[2]*wx2*volFact+0.5773502691896258*f[9]*dv2*volFact; 
  out[9] += 2.0*f[5]*wx2*volFact; 
  out[10] += 2.0*f[11]*wx2*volFact; 
  out[11] += 2.0*f[12]*wx2*volFact; 
} 
void MomentCalc2x2vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += 2.0*f[0]*wx1*volFact+0.5773502691896258*f[3]*dv1*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact+0.5773502691896258*f[6]*dv1*volFact; 
  out[2] += 2.0*f[2]*wx1*volFact+0.5773502691896258*f[7]*dv1*volFact; 
  out[3] += 2.0*f[5]*wx1*volFact+0.5773502691896258*f[15]*dv1*volFact; 
  out[4] += 2.0*f[11]*wx1*volFact+0.5773502691896257*f[21]*dv1*volFact; 
  out[5] += 2.0*f[12]*wx1*volFact+0.5773502691896257*f[22]*dv1*volFact; 
  out[6] += 2.0*f[19]*wx1*volFact; 
  out[7] += 2.0*f[20]*wx1*volFact; 
  out[8] += 2.0*f[31]*wx1*volFact; 
  out[9] += 2.0*f[32]*wx1*volFact; 
  out[10] += 2.0*f[0]*wx2*volFact+0.5773502691896258*f[4]*dv2*volFact; 
  out[11] += 2.0*f[1]*wx2*volFact+0.5773502691896258*f[8]*dv2*volFact; 
  out[12] += 2.0*f[2]*wx2*volFact+0.5773502691896258*f[9]*dv2*volFact; 
  out[13] += 2.0*f[5]*wx2*volFact+0.5773502691896258*f[16]*dv2*volFact; 
  out[14] += 2.0*f[11]*wx2*volFact+0.5773502691896257*f[25]*dv2*volFact; 
  out[15] += 2.0*f[12]*wx2*volFact+0.5773502691896257*f[26]*dv2*volFact; 
  out[16] += 2.0*f[19]*wx2*volFact; 
  out[17] += 2.0*f[20]*wx2*volFact; 
  out[18] += 2.0*f[31]*wx2*volFact; 
  out[19] += 2.0*f[32]*wx2*volFact; 
} 
void MomentCalc2x2vMax_M1i_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += 2.0*f[0]*wx1*volFact+0.5773502691896258*f[3]*dv1*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact+0.5773502691896258*f[6]*dv1*volFact; 
  out[2] += 2.0*f[2]*wx1*volFact+0.5773502691896258*f[7]*dv1*volFact; 
  out[3] += 2.0*f[5]*wx1*volFact+0.5773502691896258*f[15]*dv1*volFact; 
  out[4] += 2.0*f[11]*wx1*volFact+0.5773502691896257*f[21]*dv1*volFact; 
  out[5] += 2.0*f[12]*wx1*volFact+0.5773502691896257*f[22]*dv1*volFact; 
  out[6] += 2.0*f[19]*wx1*volFact+0.5773502691896257*f[36]*dv1*volFact; 
  out[7] += 2.0*f[20]*wx1*volFact+0.5773502691896257*f[37]*dv1*volFact; 
  out[8] += 2.0*f[31]*wx1*volFact+0.5773502691896256*f[56]*dv1*volFact; 
  out[9] += 2.0*f[32]*wx1*volFact+0.5773502691896256*f[57]*dv1*volFact; 
  out[10] += 2.0*f[48]*wx1*volFact; 
  out[11] += 2.0*f[54]*wx1*volFact; 
  out[12] += 2.0*f[55]*wx1*volFact; 
  out[13] += 2.0*f[66]*wx1*volFact; 
  out[14] += 2.0*f[67]*wx1*volFact; 
  out[15] += 2.0*f[0]*wx2*volFact+0.5773502691896258*f[4]*dv2*volFact; 
  out[16] += 2.0*f[1]*wx2*volFact+0.5773502691896258*f[8]*dv2*volFact; 
  out[17] += 2.0*f[2]*wx2*volFact+0.5773502691896258*f[9]*dv2*volFact; 
  out[18] += 2.0*f[5]*wx2*volFact+0.5773502691896258*f[16]*dv2*volFact; 
  out[19] += 2.0*f[11]*wx2*volFact+0.5773502691896257*f[25]*dv2*volFact; 
  out[20] += 2.0*f[12]*wx2*volFact+0.5773502691896257*f[26]*dv2*volFact; 
  out[21] += 2.0*f[19]*wx2*volFact+0.5773502691896257*f[39]*dv2*volFact; 
  out[22] += 2.0*f[20]*wx2*volFact+0.5773502691896257*f[40]*dv2*volFact; 
  out[23] += 2.0*f[31]*wx2*volFact+0.5773502691896256*f[60]*dv2*volFact; 
  out[24] += 2.0*f[32]*wx2*volFact+0.5773502691896256*f[61]*dv2*volFact; 
  out[25] += 2.0*f[48]*wx2*volFact; 
  out[26] += 2.0*f[54]*wx2*volFact; 
  out[27] += 2.0*f[55]*wx2*volFact; 
  out[28] += 2.0*f[66]*wx2*volFact; 
  out[29] += 2.0*f[67]*wx2*volFact; 
} 
void MomentCalc2x2vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[3]*dv1*wx1*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx1_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[0]*wx1*wx2*volFact+0.5773502691896258*f[3]*dv1*wx2*volFact+0.5773502691896258*f[4]*dv2*wx1*volFact; 
  out[4] += 2.0*f[1]*wx1*wx2*volFact; 
  out[5] += 2.0*f[2]*wx1*wx2*volFact; 
  out[6] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx2*volFact+0.1666666666666667*f[0]*dv2_sq*volFact; 
  out[7] += 2.0*f[1]*wx2_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact; 
  out[8] += 2.0*f[2]*wx2_sq*volFact+0.1666666666666667*f[2]*dv2_sq*volFact; 
} 
void MomentCalc2x2vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[3]*dv1*wx1*volFact+0.149071198499986*f[13]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[6]*dv1*wx1*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx1_sq*volFact+1.154700538379252*f[7]*dv1*wx1*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*wx1_sq*volFact+0.1666666666666667*f[5]*dv1_sq*volFact; 
  out[4] += 2.0*f[11]*wx1_sq*volFact+0.1666666666666667*f[11]*dv1_sq*volFact; 
  out[5] += 2.0*f[12]*wx1_sq*volFact+0.1666666666666667*f[12]*dv1_sq*volFact; 
  out[6] += 2.0*f[0]*wx1*wx2*volFact+0.5773502691896258*f[3]*dv1*wx2*volFact+0.5773502691896258*f[4]*dv2*wx1*volFact+0.1666666666666667*f[10]*dv1*dv2*volFact; 
  out[7] += 2.0*f[1]*wx1*wx2*volFact+0.5773502691896258*f[6]*dv1*wx2*volFact+0.5773502691896258*f[8]*dv2*wx1*volFact; 
  out[8] += 2.0*f[2]*wx1*wx2*volFact+0.5773502691896258*f[7]*dv1*wx2*volFact+0.5773502691896258*f[9]*dv2*wx1*volFact; 
  out[9] += 2.0*f[5]*wx1*wx2*volFact; 
  out[10] += 2.0*f[11]*wx1*wx2*volFact; 
  out[11] += 2.0*f[12]*wx1*wx2*volFact; 
  out[12] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx2*volFact+0.149071198499986*f[14]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact; 
  out[13] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[8]*dv2*wx2*volFact+0.1666666666666667*f[1]*dv2_sq*volFact; 
  out[14] += 2.0*f[2]*wx2_sq*volFact+1.154700538379252*f[9]*dv2*wx2*volFact+0.1666666666666667*f[2]*dv2_sq*volFact; 
  out[15] += 2.0*f[5]*wx2_sq*volFact+0.1666666666666667*f[5]*dv2_sq*volFact; 
  out[16] += 2.0*f[11]*wx2_sq*volFact+0.1666666666666667*f[11]*dv2_sq*volFact; 
  out[17] += 2.0*f[12]*wx2_sq*volFact+0.1666666666666667*f[12]*dv2_sq*volFact; 
} 
void MomentCalc2x2vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[3]*dv1*wx1*volFact+0.149071198499986*f[13]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[6]*dv1*wx1*volFact+0.149071198499986*f[23]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx1_sq*volFact+1.154700538379252*f[7]*dv1*wx1*volFact+0.149071198499986*f[24]*dv1_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*wx1_sq*volFact+1.154700538379252*f[15]*dv1*wx1*volFact+0.1666666666666667*f[5]*dv1_sq*volFact; 
  out[4] += 2.0*f[11]*wx1_sq*volFact+1.154700538379251*f[21]*dv1*wx1*volFact+0.1666666666666667*f[11]*dv1_sq*volFact; 
  out[5] += 2.0*f[12]*wx1_sq*volFact+1.154700538379251*f[22]*dv1*wx1*volFact+0.1666666666666667*f[12]*dv1_sq*volFact; 
  out[6] += 2.0*f[19]*wx1_sq*volFact+0.1666666666666667*f[19]*dv1_sq*volFact; 
  out[7] += 2.0*f[20]*wx1_sq*volFact+0.1666666666666667*f[20]*dv1_sq*volFact; 
  out[8] += 2.0*f[31]*wx1_sq*volFact+0.1666666666666667*f[31]*dv1_sq*volFact; 
  out[9] += 2.0*f[32]*wx1_sq*volFact+0.1666666666666667*f[32]*dv1_sq*volFact; 
  out[10] += 2.0*f[0]*wx1*wx2*volFact+0.5773502691896258*f[3]*dv1*wx2*volFact+0.5773502691896258*f[4]*dv2*wx1*volFact+0.1666666666666667*f[10]*dv1*dv2*volFact; 
  out[11] += 2.0*f[1]*wx1*wx2*volFact+0.5773502691896258*f[6]*dv1*wx2*volFact+0.5773502691896258*f[8]*dv2*wx1*volFact+0.1666666666666667*f[17]*dv1*dv2*volFact; 
  out[12] += 2.0*f[2]*wx1*wx2*volFact+0.5773502691896258*f[7]*dv1*wx2*volFact+0.5773502691896258*f[9]*dv2*wx1*volFact+0.1666666666666667*f[18]*dv1*dv2*volFact; 
  out[13] += 2.0*f[5]*wx1*wx2*volFact+0.5773502691896258*f[15]*dv1*wx2*volFact+0.5773502691896258*f[16]*dv2*wx1*volFact; 
  out[14] += 2.0*f[11]*wx1*wx2*volFact+0.5773502691896257*f[21]*dv1*wx2*volFact+0.5773502691896257*f[25]*dv2*wx1*volFact; 
  out[15] += 2.0*f[12]*wx1*wx2*volFact+0.5773502691896257*f[22]*dv1*wx2*volFact+0.5773502691896257*f[26]*dv2*wx1*volFact; 
  out[16] += 2.0*f[19]*wx1*wx2*volFact; 
  out[17] += 2.0*f[20]*wx1*wx2*volFact; 
  out[18] += 2.0*f[31]*wx1*wx2*volFact; 
  out[19] += 2.0*f[32]*wx1*wx2*volFact; 
  out[20] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx2*volFact+0.149071198499986*f[14]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact; 
  out[21] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[8]*dv2*wx2*volFact+0.149071198499986*f[28]*dv2_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact; 
  out[22] += 2.0*f[2]*wx2_sq*volFact+1.154700538379252*f[9]*dv2*wx2*volFact+0.149071198499986*f[29]*dv2_sq*volFact+0.1666666666666667*f[2]*dv2_sq*volFact; 
  out[23] += 2.0*f[5]*wx2_sq*volFact+1.154700538379252*f[16]*dv2*wx2*volFact+0.1666666666666667*f[5]*dv2_sq*volFact; 
  out[24] += 2.0*f[11]*wx2_sq*volFact+1.154700538379251*f[25]*dv2*wx2*volFact+0.1666666666666667*f[11]*dv2_sq*volFact; 
  out[25] += 2.0*f[12]*wx2_sq*volFact+1.154700538379251*f[26]*dv2*wx2*volFact+0.1666666666666667*f[12]*dv2_sq*volFact; 
  out[26] += 2.0*f[19]*wx2_sq*volFact+0.1666666666666667*f[19]*dv2_sq*volFact; 
  out[27] += 2.0*f[20]*wx2_sq*volFact+0.1666666666666667*f[20]*dv2_sq*volFact; 
  out[28] += 2.0*f[31]*wx2_sq*volFact+0.1666666666666667*f[31]*dv2_sq*volFact; 
  out[29] += 2.0*f[32]*wx2_sq*volFact+0.1666666666666667*f[32]*dv2_sq*volFact; 
} 
void MomentCalc2x2vMax_M2ij_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[3]*dv1*wx1*volFact+0.149071198499986*f[13]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[6]*dv1*wx1*volFact+0.149071198499986*f[23]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx1_sq*volFact+1.154700538379252*f[7]*dv1*wx1*volFact+0.149071198499986*f[24]*dv1_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*wx1_sq*volFact+1.154700538379252*f[15]*dv1*wx1*volFact+0.149071198499986*f[38]*dv1_sq*volFact+0.1666666666666667*f[5]*dv1_sq*volFact; 
  out[4] += 2.0*f[11]*wx1_sq*volFact+1.154700538379251*f[21]*dv1*wx1*volFact+0.149071198499986*f[49]*dv1_sq*volFact+0.1666666666666667*f[11]*dv1_sq*volFact; 
  out[5] += 2.0*f[12]*wx1_sq*volFact+1.154700538379251*f[22]*dv1*wx1*volFact+0.149071198499986*f[50]*dv1_sq*volFact+0.1666666666666667*f[12]*dv1_sq*volFact; 
  out[6] += 2.0*f[19]*wx1_sq*volFact+1.154700538379251*f[36]*dv1*wx1*volFact+0.1666666666666667*f[19]*dv1_sq*volFact; 
  out[7] += 2.0*f[20]*wx1_sq*volFact+1.154700538379251*f[37]*dv1*wx1*volFact+0.1666666666666667*f[20]*dv1_sq*volFact; 
  out[8] += 2.0*f[31]*wx1_sq*volFact+1.154700538379251*f[56]*dv1*wx1*volFact+0.1666666666666667*f[31]*dv1_sq*volFact; 
  out[9] += 2.0*f[32]*wx1_sq*volFact+1.154700538379251*f[57]*dv1*wx1*volFact+0.1666666666666667*f[32]*dv1_sq*volFact; 
  out[10] += 2.0*f[48]*wx1_sq*volFact+0.1666666666666667*f[48]*dv1_sq*volFact; 
  out[11] += 2.0*f[54]*wx1_sq*volFact+0.1666666666666667*f[54]*dv1_sq*volFact; 
  out[12] += 2.0*f[55]*wx1_sq*volFact+0.1666666666666667*f[55]*dv1_sq*volFact; 
  out[13] += 2.0*f[66]*wx1_sq*volFact+0.1666666666666667*f[66]*dv1_sq*volFact; 
  out[14] += 2.0*f[67]*wx1_sq*volFact+0.1666666666666667*f[67]*dv1_sq*volFact; 
  out[15] += 2.0*f[0]*wx1*wx2*volFact+0.5773502691896258*f[3]*dv1*wx2*volFact+0.5773502691896258*f[4]*dv2*wx1*volFact+0.1666666666666667*f[10]*dv1*dv2*volFact; 
  out[16] += 2.0*f[1]*wx1*wx2*volFact+0.5773502691896258*f[6]*dv1*wx2*volFact+0.5773502691896258*f[8]*dv2*wx1*volFact+0.1666666666666667*f[17]*dv1*dv2*volFact; 
  out[17] += 2.0*f[2]*wx1*wx2*volFact+0.5773502691896258*f[7]*dv1*wx2*volFact+0.5773502691896258*f[9]*dv2*wx1*volFact+0.1666666666666667*f[18]*dv1*dv2*volFact; 
  out[18] += 2.0*f[5]*wx1*wx2*volFact+0.5773502691896258*f[15]*dv1*wx2*volFact+0.5773502691896258*f[16]*dv2*wx1*volFact+0.1666666666666667*f[35]*dv1*dv2*volFact; 
  out[19] += 2.0*f[11]*wx1*wx2*volFact+0.5773502691896257*f[21]*dv1*wx2*volFact+0.5773502691896257*f[25]*dv2*wx1*volFact+0.1666666666666667*f[41]*dv1*dv2*volFact; 
  out[20] += 2.0*f[12]*wx1*wx2*volFact+0.5773502691896257*f[22]*dv1*wx2*volFact+0.5773502691896257*f[26]*dv2*wx1*volFact+0.1666666666666667*f[42]*dv1*dv2*volFact; 
  out[21] += 2.0*f[19]*wx1*wx2*volFact+0.5773502691896257*f[36]*dv1*wx2*volFact+0.5773502691896257*f[39]*dv2*wx1*volFact; 
  out[22] += 2.0*f[20]*wx1*wx2*volFact+0.5773502691896257*f[37]*dv1*wx2*volFact+0.5773502691896257*f[40]*dv2*wx1*volFact; 
  out[23] += 2.0*f[31]*wx1*wx2*volFact+0.5773502691896256*f[56]*dv1*wx2*volFact+0.5773502691896256*f[60]*dv2*wx1*volFact; 
  out[24] += 2.0*f[32]*wx1*wx2*volFact+0.5773502691896256*f[57]*dv1*wx2*volFact+0.5773502691896256*f[61]*dv2*wx1*volFact; 
  out[25] += 2.0*f[48]*wx1*wx2*volFact; 
  out[26] += 2.0*f[54]*wx1*wx2*volFact; 
  out[27] += 2.0*f[55]*wx1*wx2*volFact; 
  out[28] += 2.0*f[66]*wx1*wx2*volFact; 
  out[29] += 2.0*f[67]*wx1*wx2*volFact; 
  out[30] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx2*volFact+0.149071198499986*f[14]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact; 
  out[31] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[8]*dv2*wx2*volFact+0.149071198499986*f[28]*dv2_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact; 
  out[32] += 2.0*f[2]*wx2_sq*volFact+1.154700538379252*f[9]*dv2*wx2*volFact+0.149071198499986*f[29]*dv2_sq*volFact+0.1666666666666667*f[2]*dv2_sq*volFact; 
  out[33] += 2.0*f[5]*wx2_sq*volFact+1.154700538379252*f[16]*dv2*wx2*volFact+0.149071198499986*f[45]*dv2_sq*volFact+0.1666666666666667*f[5]*dv2_sq*volFact; 
  out[34] += 2.0*f[11]*wx2_sq*volFact+1.154700538379251*f[25]*dv2*wx2*volFact+0.149071198499986*f[51]*dv2_sq*volFact+0.1666666666666667*f[11]*dv2_sq*volFact; 
  out[35] += 2.0*f[12]*wx2_sq*volFact+1.154700538379251*f[26]*dv2*wx2*volFact+0.149071198499986*f[52]*dv2_sq*volFact+0.1666666666666667*f[12]*dv2_sq*volFact; 
  out[36] += 2.0*f[19]*wx2_sq*volFact+1.154700538379251*f[39]*dv2*wx2*volFact+0.1666666666666667*f[19]*dv2_sq*volFact; 
  out[37] += 2.0*f[20]*wx2_sq*volFact+1.154700538379251*f[40]*dv2*wx2*volFact+0.1666666666666667*f[20]*dv2_sq*volFact; 
  out[38] += 2.0*f[31]*wx2_sq*volFact+1.154700538379251*f[60]*dv2*wx2*volFact+0.1666666666666667*f[31]*dv2_sq*volFact; 
  out[39] += 2.0*f[32]*wx2_sq*volFact+1.154700538379251*f[61]*dv2*wx2*volFact+0.1666666666666667*f[32]*dv2_sq*volFact; 
  out[40] += 2.0*f[48]*wx2_sq*volFact+0.1666666666666667*f[48]*dv2_sq*volFact; 
  out[41] += 2.0*f[54]*wx2_sq*volFact+0.1666666666666667*f[54]*dv2_sq*volFact; 
  out[42] += 2.0*f[55]*wx2_sq*volFact+0.1666666666666667*f[55]*dv2_sq*volFact; 
  out[43] += 2.0*f[66]*wx2_sq*volFact+0.1666666666666667*f[66]*dv2_sq*volFact; 
  out[44] += 2.0*f[67]*wx2_sq*volFact+0.1666666666666667*f[67]*dv2_sq*volFact; 
} 
void MomentCalc2x2vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx2*volFact+2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[3]*dv1*wx1*volFact+0.1666666666666667*f[0]*dv2_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx2_sq*volFact+2.0*f[1]*wx1_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx2_sq*volFact+2.0*f[2]*wx1_sq*volFact+0.1666666666666667*f[2]*dv2_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
} 
void MomentCalc2x2vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx2*volFact+2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[3]*dv1*wx1*volFact+0.149071198499986*f[14]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact+0.149071198499986*f[13]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[8]*dv2*wx2*volFact+2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[6]*dv1*wx1*volFact+0.1666666666666667*f[1]*dv2_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx2_sq*volFact+1.154700538379252*f[9]*dv2*wx2*volFact+2.0*f[2]*wx1_sq*volFact+1.154700538379252*f[7]*dv1*wx1*volFact+0.1666666666666667*f[2]*dv2_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*wx2_sq*volFact+2.0*f[5]*wx1_sq*volFact+0.1666666666666667*f[5]*dv2_sq*volFact+0.1666666666666667*f[5]*dv1_sq*volFact; 
  out[4] += 2.0*f[11]*wx2_sq*volFact+2.0*f[11]*wx1_sq*volFact+0.1666666666666667*f[11]*dv2_sq*volFact+0.1666666666666667*f[11]*dv1_sq*volFact; 
  out[5] += 2.0*f[12]*wx2_sq*volFact+2.0*f[12]*wx1_sq*volFact+0.1666666666666667*f[12]*dv2_sq*volFact+0.1666666666666667*f[12]*dv1_sq*volFact; 
} 
void MomentCalc2x2vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx2*volFact+2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[3]*dv1*wx1*volFact+0.149071198499986*f[14]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact+0.149071198499986*f[13]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[8]*dv2*wx2*volFact+2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[6]*dv1*wx1*volFact+0.149071198499986*f[28]*dv2_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact+0.149071198499986*f[23]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx2_sq*volFact+1.154700538379252*f[9]*dv2*wx2*volFact+2.0*f[2]*wx1_sq*volFact+1.154700538379252*f[7]*dv1*wx1*volFact+0.149071198499986*f[29]*dv2_sq*volFact+0.1666666666666667*f[2]*dv2_sq*volFact+0.149071198499986*f[24]*dv1_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*wx2_sq*volFact+1.154700538379252*f[16]*dv2*wx2*volFact+2.0*f[5]*wx1_sq*volFact+1.154700538379252*f[15]*dv1*wx1*volFact+0.1666666666666667*f[5]*dv2_sq*volFact+0.1666666666666667*f[5]*dv1_sq*volFact; 
  out[4] += 2.0*f[11]*wx2_sq*volFact+1.154700538379251*f[25]*dv2*wx2*volFact+2.0*f[11]*wx1_sq*volFact+1.154700538379251*f[21]*dv1*wx1*volFact+0.1666666666666667*f[11]*dv2_sq*volFact+0.1666666666666667*f[11]*dv1_sq*volFact; 
  out[5] += 2.0*f[12]*wx2_sq*volFact+1.154700538379251*f[26]*dv2*wx2*volFact+2.0*f[12]*wx1_sq*volFact+1.154700538379251*f[22]*dv1*wx1*volFact+0.1666666666666667*f[12]*dv2_sq*volFact+0.1666666666666667*f[12]*dv1_sq*volFact; 
  out[6] += 2.0*f[19]*wx2_sq*volFact+2.0*f[19]*wx1_sq*volFact+0.1666666666666667*f[19]*dv2_sq*volFact+0.1666666666666667*f[19]*dv1_sq*volFact; 
  out[7] += 2.0*f[20]*wx2_sq*volFact+2.0*f[20]*wx1_sq*volFact+0.1666666666666667*f[20]*dv2_sq*volFact+0.1666666666666667*f[20]*dv1_sq*volFact; 
  out[8] += 2.0*f[31]*wx2_sq*volFact+2.0*f[31]*wx1_sq*volFact+0.1666666666666667*f[31]*dv2_sq*volFact+0.1666666666666667*f[31]*dv1_sq*volFact; 
  out[9] += 2.0*f[32]*wx2_sq*volFact+2.0*f[32]*wx1_sq*volFact+0.1666666666666667*f[32]*dv2_sq*volFact+0.1666666666666667*f[32]*dv1_sq*volFact; 
} 
void MomentCalc2x2vMax_M2_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx2*volFact+2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[3]*dv1*wx1*volFact+0.149071198499986*f[14]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact+0.149071198499986*f[13]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[8]*dv2*wx2*volFact+2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[6]*dv1*wx1*volFact+0.149071198499986*f[28]*dv2_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact+0.149071198499986*f[23]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx2_sq*volFact+1.154700538379252*f[9]*dv2*wx2*volFact+2.0*f[2]*wx1_sq*volFact+1.154700538379252*f[7]*dv1*wx1*volFact+0.149071198499986*f[29]*dv2_sq*volFact+0.1666666666666667*f[2]*dv2_sq*volFact+0.149071198499986*f[24]*dv1_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*wx2_sq*volFact+1.154700538379252*f[16]*dv2*wx2*volFact+2.0*f[5]*wx1_sq*volFact+1.154700538379252*f[15]*dv1*wx1*volFact+0.149071198499986*f[45]*dv2_sq*volFact+0.1666666666666667*f[5]*dv2_sq*volFact+0.149071198499986*f[38]*dv1_sq*volFact+0.1666666666666667*f[5]*dv1_sq*volFact; 
  out[4] += 2.0*f[11]*wx2_sq*volFact+1.154700538379251*f[25]*dv2*wx2*volFact+2.0*f[11]*wx1_sq*volFact+1.154700538379251*f[21]*dv1*wx1*volFact+0.149071198499986*f[51]*dv2_sq*volFact+0.1666666666666667*f[11]*dv2_sq*volFact+0.149071198499986*f[49]*dv1_sq*volFact+0.1666666666666667*f[11]*dv1_sq*volFact; 
  out[5] += 2.0*f[12]*wx2_sq*volFact+1.154700538379251*f[26]*dv2*wx2*volFact+2.0*f[12]*wx1_sq*volFact+1.154700538379251*f[22]*dv1*wx1*volFact+0.149071198499986*f[52]*dv2_sq*volFact+0.1666666666666667*f[12]*dv2_sq*volFact+0.149071198499986*f[50]*dv1_sq*volFact+0.1666666666666667*f[12]*dv1_sq*volFact; 
  out[6] += 2.0*f[19]*wx2_sq*volFact+1.154700538379251*f[39]*dv2*wx2*volFact+2.0*f[19]*wx1_sq*volFact+1.154700538379251*f[36]*dv1*wx1*volFact+0.1666666666666667*f[19]*dv2_sq*volFact+0.1666666666666667*f[19]*dv1_sq*volFact; 
  out[7] += 2.0*f[20]*wx2_sq*volFact+1.154700538379251*f[40]*dv2*wx2*volFact+2.0*f[20]*wx1_sq*volFact+1.154700538379251*f[37]*dv1*wx1*volFact+0.1666666666666667*f[20]*dv2_sq*volFact+0.1666666666666667*f[20]*dv1_sq*volFact; 
  out[8] += 2.0*f[31]*wx2_sq*volFact+1.154700538379251*f[60]*dv2*wx2*volFact+2.0*f[31]*wx1_sq*volFact+1.154700538379251*f[56]*dv1*wx1*volFact+0.1666666666666667*f[31]*dv2_sq*volFact+0.1666666666666667*f[31]*dv1_sq*volFact; 
  out[9] += 2.0*f[32]*wx2_sq*volFact+1.154700538379251*f[61]*dv2*wx2*volFact+2.0*f[32]*wx1_sq*volFact+1.154700538379251*f[57]*dv1*wx1*volFact+0.1666666666666667*f[32]*dv2_sq*volFact+0.1666666666666667*f[32]*dv1_sq*volFact; 
  out[10] += 2.0*f[48]*wx2_sq*volFact+2.0*f[48]*wx1_sq*volFact+0.1666666666666667*f[48]*dv2_sq*volFact+0.1666666666666667*f[48]*dv1_sq*volFact; 
  out[11] += 2.0*f[54]*wx2_sq*volFact+2.0*f[54]*wx1_sq*volFact+0.1666666666666667*f[54]*dv2_sq*volFact+0.1666666666666667*f[54]*dv1_sq*volFact; 
  out[12] += 2.0*f[55]*wx2_sq*volFact+2.0*f[55]*wx1_sq*volFact+0.1666666666666667*f[55]*dv2_sq*volFact+0.1666666666666667*f[55]*dv1_sq*volFact; 
  out[13] += 2.0*f[66]*wx2_sq*volFact+2.0*f[66]*wx1_sq*volFact+0.1666666666666667*f[66]*dv2_sq*volFact+0.1666666666666667*f[66]*dv1_sq*volFact; 
  out[14] += 2.0*f[67]*wx2_sq*volFact+2.0*f[67]*wx1_sq*volFact+0.1666666666666667*f[67]*dv2_sq*volFact+0.1666666666666667*f[67]*dv1_sq*volFact; 
} 
void MomentCalc2x2vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*wx1*wx2_sq*volFact+0.5773502691896258*f[3]*dv1*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx1*wx2*volFact+2.0*f[0]*wx1*wx1_sq*volFact+1.732050807568877*f[3]*dv1*wx1_sq*volFact+0.1666666666666667*f[0]*dv2_sq*wx1*volFact+0.5*f[0]*dv1_sq*wx1*volFact+0.04811252243246882*f[3]*dv1*dv2_sq*volFact+0.08660254037844387*f[3]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1*wx2_sq*volFact+2.0*f[1]*wx1*wx1_sq*volFact+0.1666666666666667*f[1]*dv2_sq*wx1*volFact+0.5*f[1]*dv1_sq*wx1*volFact; 
  out[2] += 2.0*f[2]*wx1*wx2_sq*volFact+2.0*f[2]*wx1*wx1_sq*volFact+0.1666666666666667*f[2]*dv2_sq*wx1*volFact+0.5*f[2]*dv1_sq*wx1*volFact; 
  out[3] += 2.0*f[0]*wx2*wx2_sq*volFact+1.732050807568877*f[4]*dv2*wx2_sq*volFact+2.0*f[0]*wx1_sq*wx2*volFact+1.154700538379252*f[3]*dv1*wx1*wx2*volFact+0.5*f[0]*dv2_sq*wx2*volFact+0.1666666666666667*f[0]*dv1_sq*wx2*volFact+0.5773502691896258*f[4]*dv2*wx1_sq*volFact+0.08660254037844387*f[4]*dv2*dv2_sq*volFact+0.04811252243246882*f[4]*dv1_sq*dv2*volFact; 
  out[4] += 2.0*f[1]*wx2*wx2_sq*volFact+2.0*f[1]*wx1_sq*wx2*volFact+0.5*f[1]*dv2_sq*wx2*volFact+0.1666666666666667*f[1]*dv1_sq*wx2*volFact; 
  out[5] += 2.0*f[2]*wx2*wx2_sq*volFact+2.0*f[2]*wx1_sq*wx2*volFact+0.5*f[2]*dv2_sq*wx2*volFact+0.1666666666666667*f[2]*dv1_sq*wx2*volFact; 
} 
void MomentCalc2x2vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*wx1*wx2_sq*volFact+0.5773502691896258*f[3]*dv1*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx1*wx2*volFact+0.3333333333333333*f[10]*dv1*dv2*wx2*volFact+2.0*f[0]*wx1*wx1_sq*volFact+1.732050807568877*f[3]*dv1*wx1_sq*volFact+0.149071198499986*f[14]*dv2_sq*wx1*volFact+0.1666666666666667*f[0]*dv2_sq*wx1*volFact+0.4472135954999579*f[13]*dv1_sq*wx1*volFact+0.5*f[0]*dv1_sq*wx1*volFact+0.04811252243246882*f[3]*dv1*dv2_sq*volFact+0.08660254037844387*f[3]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1*wx2_sq*volFact+0.5773502691896258*f[6]*dv1*wx2_sq*volFact+1.154700538379252*f[8]*dv2*wx1*wx2*volFact+2.0*f[1]*wx1*wx1_sq*volFact+1.732050807568877*f[6]*dv1*wx1_sq*volFact+0.1666666666666667*f[1]*dv2_sq*wx1*volFact+0.5*f[1]*dv1_sq*wx1*volFact+0.04811252243246882*f[6]*dv1*dv2_sq*volFact+0.08660254037844387*f[6]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx1*wx2_sq*volFact+0.5773502691896258*f[7]*dv1*wx2_sq*volFact+1.154700538379252*f[9]*dv2*wx1*wx2*volFact+2.0*f[2]*wx1*wx1_sq*volFact+1.732050807568877*f[7]*dv1*wx1_sq*volFact+0.1666666666666667*f[2]*dv2_sq*wx1*volFact+0.5*f[2]*dv1_sq*wx1*volFact+0.04811252243246882*f[7]*dv1*dv2_sq*volFact+0.08660254037844387*f[7]*dv1*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*wx1*wx2_sq*volFact+2.0*f[5]*wx1*wx1_sq*volFact+0.1666666666666667*f[5]*dv2_sq*wx1*volFact+0.5*f[5]*dv1_sq*wx1*volFact; 
  out[4] += 2.0*f[11]*wx1*wx2_sq*volFact+2.0*f[11]*wx1*wx1_sq*volFact+0.1666666666666667*f[11]*dv2_sq*wx1*volFact+0.5*f[11]*dv1_sq*wx1*volFact; 
  out[5] += 2.0*f[12]*wx1*wx2_sq*volFact+2.0*f[12]*wx1*wx1_sq*volFact+0.1666666666666667*f[12]*dv2_sq*wx1*volFact+0.5*f[12]*dv1_sq*wx1*volFact; 
  out[6] += 2.0*f[0]*wx2*wx2_sq*volFact+1.732050807568877*f[4]*dv2*wx2_sq*volFact+2.0*f[0]*wx1_sq*wx2*volFact+1.154700538379252*f[3]*dv1*wx1*wx2*volFact+0.4472135954999579*f[14]*dv2_sq*wx2*volFact+0.5*f[0]*dv2_sq*wx2*volFact+0.149071198499986*f[13]*dv1_sq*wx2*volFact+0.1666666666666667*f[0]*dv1_sq*wx2*volFact+0.5773502691896258*f[4]*dv2*wx1_sq*volFact+0.3333333333333333*f[10]*dv1*dv2*wx1*volFact+0.08660254037844387*f[4]*dv2*dv2_sq*volFact+0.04811252243246882*f[4]*dv1_sq*dv2*volFact; 
  out[7] += 2.0*f[1]*wx2*wx2_sq*volFact+1.732050807568877*f[8]*dv2*wx2_sq*volFact+2.0*f[1]*wx1_sq*wx2*volFact+1.154700538379252*f[6]*dv1*wx1*wx2*volFact+0.5*f[1]*dv2_sq*wx2*volFact+0.1666666666666667*f[1]*dv1_sq*wx2*volFact+0.5773502691896258*f[8]*dv2*wx1_sq*volFact+0.08660254037844387*f[8]*dv2*dv2_sq*volFact+0.04811252243246882*f[8]*dv1_sq*dv2*volFact; 
  out[8] += 2.0*f[2]*wx2*wx2_sq*volFact+1.732050807568877*f[9]*dv2*wx2_sq*volFact+2.0*f[2]*wx1_sq*wx2*volFact+1.154700538379252*f[7]*dv1*wx1*wx2*volFact+0.5*f[2]*dv2_sq*wx2*volFact+0.1666666666666667*f[2]*dv1_sq*wx2*volFact+0.5773502691896258*f[9]*dv2*wx1_sq*volFact+0.08660254037844387*f[9]*dv2*dv2_sq*volFact+0.04811252243246882*f[9]*dv1_sq*dv2*volFact; 
  out[9] += 2.0*f[5]*wx2*wx2_sq*volFact+2.0*f[5]*wx1_sq*wx2*volFact+0.5*f[5]*dv2_sq*wx2*volFact+0.1666666666666667*f[5]*dv1_sq*wx2*volFact; 
  out[10] += 2.0*f[11]*wx2*wx2_sq*volFact+2.0*f[11]*wx1_sq*wx2*volFact+0.5*f[11]*dv2_sq*wx2*volFact+0.1666666666666667*f[11]*dv1_sq*wx2*volFact; 
  out[11] += 2.0*f[12]*wx2*wx2_sq*volFact+2.0*f[12]*wx1_sq*wx2*volFact+0.5*f[12]*dv2_sq*wx2*volFact+0.1666666666666667*f[12]*dv1_sq*wx2*volFact; 
} 
void MomentCalc2x2vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*wx1*wx2_sq*volFact+0.5773502691896258*f[3]*dv1*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx1*wx2*volFact+0.3333333333333333*f[10]*dv1*dv2*wx2*volFact+2.0*f[0]*wx1*wx1_sq*volFact+1.732050807568877*f[3]*dv1*wx1_sq*volFact+0.149071198499986*f[14]*dv2_sq*wx1*volFact+0.1666666666666667*f[0]*dv2_sq*wx1*volFact+0.4472135954999579*f[13]*dv1_sq*wx1*volFact+0.5*f[0]*dv1_sq*wx1*volFact+0.04303314829119351*f[30]*dv1*dv2_sq*volFact+0.04811252243246882*f[3]*dv1*dv2_sq*volFact+0.03779644730092272*f[33]*dv1*dv1_sq*volFact+0.08660254037844387*f[3]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1*wx2_sq*volFact+0.5773502691896258*f[6]*dv1*wx2_sq*volFact+1.154700538379252*f[8]*dv2*wx1*wx2*volFact+0.3333333333333333*f[17]*dv1*dv2*wx2*volFact+2.0*f[1]*wx1*wx1_sq*volFact+1.732050807568877*f[6]*dv1*wx1_sq*volFact+0.149071198499986*f[28]*dv2_sq*wx1*volFact+0.1666666666666667*f[1]*dv2_sq*wx1*volFact+0.447213595499958*f[23]*dv1_sq*wx1*volFact+0.5*f[1]*dv1_sq*wx1*volFact+0.04811252243246882*f[6]*dv1*dv2_sq*volFact+0.08660254037844387*f[6]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx1*wx2_sq*volFact+0.5773502691896258*f[7]*dv1*wx2_sq*volFact+1.154700538379252*f[9]*dv2*wx1*wx2*volFact+0.3333333333333333*f[18]*dv1*dv2*wx2*volFact+2.0*f[2]*wx1*wx1_sq*volFact+1.732050807568877*f[7]*dv1*wx1_sq*volFact+0.149071198499986*f[29]*dv2_sq*wx1*volFact+0.1666666666666667*f[2]*dv2_sq*wx1*volFact+0.447213595499958*f[24]*dv1_sq*wx1*volFact+0.5*f[2]*dv1_sq*wx1*volFact+0.04811252243246882*f[7]*dv1*dv2_sq*volFact+0.08660254037844387*f[7]*dv1*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*wx1*wx2_sq*volFact+0.5773502691896258*f[15]*dv1*wx2_sq*volFact+1.154700538379252*f[16]*dv2*wx1*wx2*volFact+2.0*f[5]*wx1*wx1_sq*volFact+1.732050807568877*f[15]*dv1*wx1_sq*volFact+0.1666666666666667*f[5]*dv2_sq*wx1*volFact+0.5*f[5]*dv1_sq*wx1*volFact+0.04811252243246882*f[15]*dv1*dv2_sq*volFact+0.08660254037844387*f[15]*dv1*dv1_sq*volFact; 
  out[4] += 2.0*f[11]*wx1*wx2_sq*volFact+0.5773502691896257*f[21]*dv1*wx2_sq*volFact+1.154700538379251*f[25]*dv2*wx1*wx2*volFact+2.0*f[11]*wx1*wx1_sq*volFact+1.732050807568877*f[21]*dv1*wx1_sq*volFact+0.1666666666666667*f[11]*dv2_sq*wx1*volFact+0.5*f[11]*dv1_sq*wx1*volFact+0.04811252243246881*f[21]*dv1*dv2_sq*volFact+0.08660254037844385*f[21]*dv1*dv1_sq*volFact; 
  out[5] += 2.0*f[12]*wx1*wx2_sq*volFact+0.5773502691896257*f[22]*dv1*wx2_sq*volFact+1.154700538379251*f[26]*dv2*wx1*wx2*volFact+2.0*f[12]*wx1*wx1_sq*volFact+1.732050807568877*f[22]*dv1*wx1_sq*volFact+0.1666666666666667*f[12]*dv2_sq*wx1*volFact+0.5*f[12]*dv1_sq*wx1*volFact+0.04811252243246881*f[22]*dv1*dv2_sq*volFact+0.08660254037844385*f[22]*dv1*dv1_sq*volFact; 
  out[6] += 2.0*f[19]*wx1*wx2_sq*volFact+2.0*f[19]*wx1*wx1_sq*volFact+0.1666666666666667*f[19]*dv2_sq*wx1*volFact+0.5*f[19]*dv1_sq*wx1*volFact; 
  out[7] += 2.0*f[20]*wx1*wx2_sq*volFact+2.0*f[20]*wx1*wx1_sq*volFact+0.1666666666666667*f[20]*dv2_sq*wx1*volFact+0.5*f[20]*dv1_sq*wx1*volFact; 
  out[8] += 2.0*f[31]*wx1*wx2_sq*volFact+2.0*f[31]*wx1*wx1_sq*volFact+0.1666666666666667*f[31]*dv2_sq*wx1*volFact+0.5*f[31]*dv1_sq*wx1*volFact; 
  out[9] += 2.0*f[32]*wx1*wx2_sq*volFact+2.0*f[32]*wx1*wx1_sq*volFact+0.1666666666666667*f[32]*dv2_sq*wx1*volFact+0.5*f[32]*dv1_sq*wx1*volFact; 
  out[10] += 2.0*f[0]*wx2*wx2_sq*volFact+1.732050807568877*f[4]*dv2*wx2_sq*volFact+2.0*f[0]*wx1_sq*wx2*volFact+1.154700538379252*f[3]*dv1*wx1*wx2*volFact+0.4472135954999579*f[14]*dv2_sq*wx2*volFact+0.5*f[0]*dv2_sq*wx2*volFact+0.149071198499986*f[13]*dv1_sq*wx2*volFact+0.1666666666666667*f[0]*dv1_sq*wx2*volFact+0.5773502691896258*f[4]*dv2*wx1_sq*volFact+0.3333333333333333*f[10]*dv1*dv2*wx1*volFact+0.03779644730092272*f[34]*dv2*dv2_sq*volFact+0.08660254037844387*f[4]*dv2*dv2_sq*volFact+0.04303314829119351*f[27]*dv1_sq*dv2*volFact+0.04811252243246882*f[4]*dv1_sq*dv2*volFact; 
  out[11] += 2.0*f[1]*wx2*wx2_sq*volFact+1.732050807568877*f[8]*dv2*wx2_sq*volFact+2.0*f[1]*wx1_sq*wx2*volFact+1.154700538379252*f[6]*dv1*wx1*wx2*volFact+0.447213595499958*f[28]*dv2_sq*wx2*volFact+0.5*f[1]*dv2_sq*wx2*volFact+0.149071198499986*f[23]*dv1_sq*wx2*volFact+0.1666666666666667*f[1]*dv1_sq*wx2*volFact+0.5773502691896258*f[8]*dv2*wx1_sq*volFact+0.3333333333333333*f[17]*dv1*dv2*wx1*volFact+0.08660254037844387*f[8]*dv2*dv2_sq*volFact+0.04811252243246882*f[8]*dv1_sq*dv2*volFact; 
  out[12] += 2.0*f[2]*wx2*wx2_sq*volFact+1.732050807568877*f[9]*dv2*wx2_sq*volFact+2.0*f[2]*wx1_sq*wx2*volFact+1.154700538379252*f[7]*dv1*wx1*wx2*volFact+0.447213595499958*f[29]*dv2_sq*wx2*volFact+0.5*f[2]*dv2_sq*wx2*volFact+0.149071198499986*f[24]*dv1_sq*wx2*volFact+0.1666666666666667*f[2]*dv1_sq*wx2*volFact+0.5773502691896258*f[9]*dv2*wx1_sq*volFact+0.3333333333333333*f[18]*dv1*dv2*wx1*volFact+0.08660254037844387*f[9]*dv2*dv2_sq*volFact+0.04811252243246882*f[9]*dv1_sq*dv2*volFact; 
  out[13] += 2.0*f[5]*wx2*wx2_sq*volFact+1.732050807568877*f[16]*dv2*wx2_sq*volFact+2.0*f[5]*wx1_sq*wx2*volFact+1.154700538379252*f[15]*dv1*wx1*wx2*volFact+0.5*f[5]*dv2_sq*wx2*volFact+0.1666666666666667*f[5]*dv1_sq*wx2*volFact+0.5773502691896258*f[16]*dv2*wx1_sq*volFact+0.08660254037844387*f[16]*dv2*dv2_sq*volFact+0.04811252243246882*f[16]*dv1_sq*dv2*volFact; 
  out[14] += 2.0*f[11]*wx2*wx2_sq*volFact+1.732050807568877*f[25]*dv2*wx2_sq*volFact+2.0*f[11]*wx1_sq*wx2*volFact+1.154700538379251*f[21]*dv1*wx1*wx2*volFact+0.5*f[11]*dv2_sq*wx2*volFact+0.1666666666666667*f[11]*dv1_sq*wx2*volFact+0.5773502691896257*f[25]*dv2*wx1_sq*volFact+0.08660254037844385*f[25]*dv2*dv2_sq*volFact+0.04811252243246881*f[25]*dv1_sq*dv2*volFact; 
  out[15] += 2.0*f[12]*wx2*wx2_sq*volFact+1.732050807568877*f[26]*dv2*wx2_sq*volFact+2.0*f[12]*wx1_sq*wx2*volFact+1.154700538379251*f[22]*dv1*wx1*wx2*volFact+0.5*f[12]*dv2_sq*wx2*volFact+0.1666666666666667*f[12]*dv1_sq*wx2*volFact+0.5773502691896257*f[26]*dv2*wx1_sq*volFact+0.08660254037844385*f[26]*dv2*dv2_sq*volFact+0.04811252243246881*f[26]*dv1_sq*dv2*volFact; 
  out[16] += 2.0*f[19]*wx2*wx2_sq*volFact+2.0*f[19]*wx1_sq*wx2*volFact+0.5*f[19]*dv2_sq*wx2*volFact+0.1666666666666667*f[19]*dv1_sq*wx2*volFact; 
  out[17] += 2.0*f[20]*wx2*wx2_sq*volFact+2.0*f[20]*wx1_sq*wx2*volFact+0.5*f[20]*dv2_sq*wx2*volFact+0.1666666666666667*f[20]*dv1_sq*wx2*volFact; 
  out[18] += 2.0*f[31]*wx2*wx2_sq*volFact+2.0*f[31]*wx1_sq*wx2*volFact+0.5*f[31]*dv2_sq*wx2*volFact+0.1666666666666667*f[31]*dv1_sq*wx2*volFact; 
  out[19] += 2.0*f[32]*wx2*wx2_sq*volFact+2.0*f[32]*wx1_sq*wx2*volFact+0.5*f[32]*dv2_sq*wx2*volFact+0.1666666666666667*f[32]*dv1_sq*wx2*volFact; 
} 
void MomentCalc2x2vMax_M3i_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*wx1*wx2_sq*volFact+0.5773502691896258*f[3]*dv1*wx2_sq*volFact+1.154700538379252*f[4]*dv2*wx1*wx2*volFact+0.3333333333333333*f[10]*dv1*dv2*wx2*volFact+2.0*f[0]*wx1*wx1_sq*volFact+1.732050807568877*f[3]*dv1*wx1_sq*volFact+0.149071198499986*f[14]*dv2_sq*wx1*volFact+0.1666666666666667*f[0]*dv2_sq*wx1*volFact+0.4472135954999579*f[13]*dv1_sq*wx1*volFact+0.5*f[0]*dv1_sq*wx1*volFact+0.04303314829119351*f[30]*dv1*dv2_sq*volFact+0.04811252243246882*f[3]*dv1*dv2_sq*volFact+0.03779644730092272*f[33]*dv1*dv1_sq*volFact+0.08660254037844387*f[3]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1*wx2_sq*volFact+0.5773502691896258*f[6]*dv1*wx2_sq*volFact+1.154700538379252*f[8]*dv2*wx1*wx2*volFact+0.3333333333333333*f[17]*dv1*dv2*wx2*volFact+2.0*f[1]*wx1*wx1_sq*volFact+1.732050807568877*f[6]*dv1*wx1_sq*volFact+0.149071198499986*f[28]*dv2_sq*wx1*volFact+0.1666666666666667*f[1]*dv2_sq*wx1*volFact+0.447213595499958*f[23]*dv1_sq*wx1*volFact+0.5*f[1]*dv1_sq*wx1*volFact+0.04303314829119353*f[46]*dv1*dv2_sq*volFact+0.04811252243246882*f[6]*dv1*dv2_sq*volFact+0.03779644730092273*f[58]*dv1*dv1_sq*volFact+0.08660254037844387*f[6]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*wx1*wx2_sq*volFact+0.5773502691896258*f[7]*dv1*wx2_sq*volFact+1.154700538379252*f[9]*dv2*wx1*wx2*volFact+0.3333333333333333*f[18]*dv1*dv2*wx2*volFact+2.0*f[2]*wx1*wx1_sq*volFact+1.732050807568877*f[7]*dv1*wx1_sq*volFact+0.149071198499986*f[29]*dv2_sq*wx1*volFact+0.1666666666666667*f[2]*dv2_sq*wx1*volFact+0.447213595499958*f[24]*dv1_sq*wx1*volFact+0.5*f[2]*dv1_sq*wx1*volFact+0.04303314829119353*f[47]*dv1*dv2_sq*volFact+0.04811252243246882*f[7]*dv1*dv2_sq*volFact+0.03779644730092273*f[59]*dv1*dv1_sq*volFact+0.08660254037844387*f[7]*dv1*dv1_sq*volFact; 
  out[3] += 2.0*f[5]*wx1*wx2_sq*volFact+0.5773502691896258*f[15]*dv1*wx2_sq*volFact+1.154700538379252*f[16]*dv2*wx1*wx2*volFact+0.3333333333333333*f[35]*dv1*dv2*wx2*volFact+2.0*f[5]*wx1*wx1_sq*volFact+1.732050807568877*f[15]*dv1*wx1_sq*volFact+0.149071198499986*f[45]*dv2_sq*wx1*volFact+0.1666666666666667*f[5]*dv2_sq*wx1*volFact+0.4472135954999579*f[38]*dv1_sq*wx1*volFact+0.5*f[5]*dv1_sq*wx1*volFact+0.04811252243246882*f[15]*dv1*dv2_sq*volFact+0.08660254037844387*f[15]*dv1*dv1_sq*volFact; 
  out[4] += 2.0*f[11]*wx1*wx2_sq*volFact+0.5773502691896257*f[21]*dv1*wx2_sq*volFact+1.154700538379251*f[25]*dv2*wx1*wx2*volFact+0.3333333333333333*f[41]*dv1*dv2*wx2*volFact+2.0*f[11]*wx1*wx1_sq*volFact+1.732050807568877*f[21]*dv1*wx1_sq*volFact+0.149071198499986*f[51]*dv2_sq*wx1*volFact+0.1666666666666667*f[11]*dv2_sq*wx1*volFact+0.4472135954999579*f[49]*dv1_sq*wx1*volFact+0.5*f[11]*dv1_sq*wx1*volFact+0.04811252243246881*f[21]*dv1*dv2_sq*volFact+0.08660254037844385*f[21]*dv1*dv1_sq*volFact; 
  out[5] += 2.0*f[12]*wx1*wx2_sq*volFact+0.5773502691896257*f[22]*dv1*wx2_sq*volFact+1.154700538379251*f[26]*dv2*wx1*wx2*volFact+0.3333333333333333*f[42]*dv1*dv2*wx2*volFact+2.0*f[12]*wx1*wx1_sq*volFact+1.732050807568877*f[22]*dv1*wx1_sq*volFact+0.149071198499986*f[52]*dv2_sq*wx1*volFact+0.1666666666666667*f[12]*dv2_sq*wx1*volFact+0.4472135954999579*f[50]*dv1_sq*wx1*volFact+0.5*f[12]*dv1_sq*wx1*volFact+0.04811252243246881*f[22]*dv1*dv2_sq*volFact+0.08660254037844385*f[22]*dv1*dv1_sq*volFact; 
  out[6] += 2.0*f[19]*wx1*wx2_sq*volFact+0.5773502691896257*f[36]*dv1*wx2_sq*volFact+1.154700538379251*f[39]*dv2*wx1*wx2*volFact+2.0*f[19]*wx1*wx1_sq*volFact+1.732050807568877*f[36]*dv1*wx1_sq*volFact+0.1666666666666667*f[19]*dv2_sq*wx1*volFact+0.5*f[19]*dv1_sq*wx1*volFact+0.04811252243246881*f[36]*dv1*dv2_sq*volFact+0.08660254037844385*f[36]*dv1*dv1_sq*volFact; 
  out[7] += 2.0*f[20]*wx1*wx2_sq*volFact+0.5773502691896257*f[37]*dv1*wx2_sq*volFact+1.154700538379251*f[40]*dv2*wx1*wx2*volFact+2.0*f[20]*wx1*wx1_sq*volFact+1.732050807568877*f[37]*dv1*wx1_sq*volFact+0.1666666666666667*f[20]*dv2_sq*wx1*volFact+0.5*f[20]*dv1_sq*wx1*volFact+0.04811252243246881*f[37]*dv1*dv2_sq*volFact+0.08660254037844385*f[37]*dv1*dv1_sq*volFact; 
  out[8] += 2.0*f[31]*wx1*wx2_sq*volFact+0.5773502691896256*f[56]*dv1*wx2_sq*volFact+1.154700538379251*f[60]*dv2*wx1*wx2*volFact+2.0*f[31]*wx1*wx1_sq*volFact+1.732050807568877*f[56]*dv1*wx1_sq*volFact+0.1666666666666667*f[31]*dv2_sq*wx1*volFact+0.5*f[31]*dv1_sq*wx1*volFact+0.0481125224324688*f[56]*dv1*dv2_sq*volFact+0.08660254037844385*f[56]*dv1*dv1_sq*volFact; 
  out[9] += 2.0*f[32]*wx1*wx2_sq*volFact+0.5773502691896256*f[57]*dv1*wx2_sq*volFact+1.154700538379251*f[61]*dv2*wx1*wx2*volFact+2.0*f[32]*wx1*wx1_sq*volFact+1.732050807568877*f[57]*dv1*wx1_sq*volFact+0.1666666666666667*f[32]*dv2_sq*wx1*volFact+0.5*f[32]*dv1_sq*wx1*volFact+0.0481125224324688*f[57]*dv1*dv2_sq*volFact+0.08660254037844385*f[57]*dv1*dv1_sq*volFact; 
  out[10] += 2.0*f[48]*wx1*wx2_sq*volFact+2.0*f[48]*wx1*wx1_sq*volFact+0.1666666666666667*f[48]*dv2_sq*wx1*volFact+0.5*f[48]*dv1_sq*wx1*volFact; 
  out[11] += 2.0*f[54]*wx1*wx2_sq*volFact+2.0*f[54]*wx1*wx1_sq*volFact+0.1666666666666667*f[54]*dv2_sq*wx1*volFact+0.5*f[54]*dv1_sq*wx1*volFact; 
  out[12] += 2.0*f[55]*wx1*wx2_sq*volFact+2.0*f[55]*wx1*wx1_sq*volFact+0.1666666666666667*f[55]*dv2_sq*wx1*volFact+0.5*f[55]*dv1_sq*wx1*volFact; 
  out[13] += 2.0*f[66]*wx1*wx2_sq*volFact+2.0*f[66]*wx1*wx1_sq*volFact+0.1666666666666667*f[66]*dv2_sq*wx1*volFact+0.5*f[66]*dv1_sq*wx1*volFact; 
  out[14] += 2.0*f[67]*wx1*wx2_sq*volFact+2.0*f[67]*wx1*wx1_sq*volFact+0.1666666666666667*f[67]*dv2_sq*wx1*volFact+0.5*f[67]*dv1_sq*wx1*volFact; 
  out[15] += 2.0*f[0]*wx2*wx2_sq*volFact+1.732050807568877*f[4]*dv2*wx2_sq*volFact+2.0*f[0]*wx1_sq*wx2*volFact+1.154700538379252*f[3]*dv1*wx1*wx2*volFact+0.4472135954999579*f[14]*dv2_sq*wx2*volFact+0.5*f[0]*dv2_sq*wx2*volFact+0.149071198499986*f[13]*dv1_sq*wx2*volFact+0.1666666666666667*f[0]*dv1_sq*wx2*volFact+0.5773502691896258*f[4]*dv2*wx1_sq*volFact+0.3333333333333333*f[10]*dv1*dv2*wx1*volFact+0.03779644730092272*f[34]*dv2*dv2_sq*volFact+0.08660254037844387*f[4]*dv2*dv2_sq*volFact+0.04303314829119351*f[27]*dv1_sq*dv2*volFact+0.04811252243246882*f[4]*dv1_sq*dv2*volFact; 
  out[16] += 2.0*f[1]*wx2*wx2_sq*volFact+1.732050807568877*f[8]*dv2*wx2_sq*volFact+2.0*f[1]*wx1_sq*wx2*volFact+1.154700538379252*f[6]*dv1*wx1*wx2*volFact+0.447213595499958*f[28]*dv2_sq*wx2*volFact+0.5*f[1]*dv2_sq*wx2*volFact+0.149071198499986*f[23]*dv1_sq*wx2*volFact+0.1666666666666667*f[1]*dv1_sq*wx2*volFact+0.5773502691896258*f[8]*dv2*wx1_sq*volFact+0.3333333333333333*f[17]*dv1*dv2*wx1*volFact+0.03779644730092273*f[63]*dv2*dv2_sq*volFact+0.08660254037844387*f[8]*dv2*dv2_sq*volFact+0.04303314829119353*f[43]*dv1_sq*dv2*volFact+0.04811252243246882*f[8]*dv1_sq*dv2*volFact; 
  out[17] += 2.0*f[2]*wx2*wx2_sq*volFact+1.732050807568877*f[9]*dv2*wx2_sq*volFact+2.0*f[2]*wx1_sq*wx2*volFact+1.154700538379252*f[7]*dv1*wx1*wx2*volFact+0.447213595499958*f[29]*dv2_sq*wx2*volFact+0.5*f[2]*dv2_sq*wx2*volFact+0.149071198499986*f[24]*dv1_sq*wx2*volFact+0.1666666666666667*f[2]*dv1_sq*wx2*volFact+0.5773502691896258*f[9]*dv2*wx1_sq*volFact+0.3333333333333333*f[18]*dv1*dv2*wx1*volFact+0.03779644730092273*f[64]*dv2*dv2_sq*volFact+0.08660254037844387*f[9]*dv2*dv2_sq*volFact+0.04303314829119353*f[44]*dv1_sq*dv2*volFact+0.04811252243246882*f[9]*dv1_sq*dv2*volFact; 
  out[18] += 2.0*f[5]*wx2*wx2_sq*volFact+1.732050807568877*f[16]*dv2*wx2_sq*volFact+2.0*f[5]*wx1_sq*wx2*volFact+1.154700538379252*f[15]*dv1*wx1*wx2*volFact+0.4472135954999579*f[45]*dv2_sq*wx2*volFact+0.5*f[5]*dv2_sq*wx2*volFact+0.149071198499986*f[38]*dv1_sq*wx2*volFact+0.1666666666666667*f[5]*dv1_sq*wx2*volFact+0.5773502691896258*f[16]*dv2*wx1_sq*volFact+0.3333333333333333*f[35]*dv1*dv2*wx1*volFact+0.08660254037844387*f[16]*dv2*dv2_sq*volFact+0.04811252243246882*f[16]*dv1_sq*dv2*volFact; 
  out[19] += 2.0*f[11]*wx2*wx2_sq*volFact+1.732050807568877*f[25]*dv2*wx2_sq*volFact+2.0*f[11]*wx1_sq*wx2*volFact+1.154700538379251*f[21]*dv1*wx1*wx2*volFact+0.4472135954999579*f[51]*dv2_sq*wx2*volFact+0.5*f[11]*dv2_sq*wx2*volFact+0.149071198499986*f[49]*dv1_sq*wx2*volFact+0.1666666666666667*f[11]*dv1_sq*wx2*volFact+0.5773502691896257*f[25]*dv2*wx1_sq*volFact+0.3333333333333333*f[41]*dv1*dv2*wx1*volFact+0.08660254037844385*f[25]*dv2*dv2_sq*volFact+0.04811252243246881*f[25]*dv1_sq*dv2*volFact; 
  out[20] += 2.0*f[12]*wx2*wx2_sq*volFact+1.732050807568877*f[26]*dv2*wx2_sq*volFact+2.0*f[12]*wx1_sq*wx2*volFact+1.154700538379251*f[22]*dv1*wx1*wx2*volFact+0.4472135954999579*f[52]*dv2_sq*wx2*volFact+0.5*f[12]*dv2_sq*wx2*volFact+0.149071198499986*f[50]*dv1_sq*wx2*volFact+0.1666666666666667*f[12]*dv1_sq*wx2*volFact+0.5773502691896257*f[26]*dv2*wx1_sq*volFact+0.3333333333333333*f[42]*dv1*dv2*wx1*volFact+0.08660254037844385*f[26]*dv2*dv2_sq*volFact+0.04811252243246881*f[26]*dv1_sq*dv2*volFact; 
  out[21] += 2.0*f[19]*wx2*wx2_sq*volFact+1.732050807568877*f[39]*dv2*wx2_sq*volFact+2.0*f[19]*wx1_sq*wx2*volFact+1.154700538379251*f[36]*dv1*wx1*wx2*volFact+0.5*f[19]*dv2_sq*wx2*volFact+0.1666666666666667*f[19]*dv1_sq*wx2*volFact+0.5773502691896257*f[39]*dv2*wx1_sq*volFact+0.08660254037844385*f[39]*dv2*dv2_sq*volFact+0.04811252243246881*f[39]*dv1_sq*dv2*volFact; 
  out[22] += 2.0*f[20]*wx2*wx2_sq*volFact+1.732050807568877*f[40]*dv2*wx2_sq*volFact+2.0*f[20]*wx1_sq*wx2*volFact+1.154700538379251*f[37]*dv1*wx1*wx2*volFact+0.5*f[20]*dv2_sq*wx2*volFact+0.1666666666666667*f[20]*dv1_sq*wx2*volFact+0.5773502691896257*f[40]*dv2*wx1_sq*volFact+0.08660254037844385*f[40]*dv2*dv2_sq*volFact+0.04811252243246881*f[40]*dv1_sq*dv2*volFact; 
  out[23] += 2.0*f[31]*wx2*wx2_sq*volFact+1.732050807568877*f[60]*dv2*wx2_sq*volFact+2.0*f[31]*wx1_sq*wx2*volFact+1.154700538379251*f[56]*dv1*wx1*wx2*volFact+0.5*f[31]*dv2_sq*wx2*volFact+0.1666666666666667*f[31]*dv1_sq*wx2*volFact+0.5773502691896256*f[60]*dv2*wx1_sq*volFact+0.08660254037844385*f[60]*dv2*dv2_sq*volFact+0.0481125224324688*f[60]*dv1_sq*dv2*volFact; 
  out[24] += 2.0*f[32]*wx2*wx2_sq*volFact+1.732050807568877*f[61]*dv2*wx2_sq*volFact+2.0*f[32]*wx1_sq*wx2*volFact+1.154700538379251*f[57]*dv1*wx1*wx2*volFact+0.5*f[32]*dv2_sq*wx2*volFact+0.1666666666666667*f[32]*dv1_sq*wx2*volFact+0.5773502691896256*f[61]*dv2*wx1_sq*volFact+0.08660254037844385*f[61]*dv2*dv2_sq*volFact+0.0481125224324688*f[61]*dv1_sq*dv2*volFact; 
  out[25] += 2.0*f[48]*wx2*wx2_sq*volFact+2.0*f[48]*wx1_sq*wx2*volFact+0.5*f[48]*dv2_sq*wx2*volFact+0.1666666666666667*f[48]*dv1_sq*wx2*volFact; 
  out[26] += 2.0*f[54]*wx2*wx2_sq*volFact+2.0*f[54]*wx1_sq*wx2*volFact+0.5*f[54]*dv2_sq*wx2*volFact+0.1666666666666667*f[54]*dv1_sq*wx2*volFact; 
  out[27] += 2.0*f[55]*wx2*wx2_sq*volFact+2.0*f[55]*wx1_sq*wx2*volFact+0.5*f[55]*dv2_sq*wx2*volFact+0.1666666666666667*f[55]*dv1_sq*wx2*volFact; 
  out[28] += 2.0*f[66]*wx2*wx2_sq*volFact+2.0*f[66]*wx1_sq*wx2*volFact+0.5*f[66]*dv2_sq*wx2*volFact+0.1666666666666667*f[66]*dv1_sq*wx2*volFact; 
  out[29] += 2.0*f[67]*wx2*wx2_sq*volFact+2.0*f[67]*wx1_sq*wx2*volFact+0.5*f[67]*dv2_sq*wx2*volFact+0.1666666666666667*f[67]*dv1_sq*wx2*volFact; 
} 
