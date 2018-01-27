#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc1x2vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
void MomentCalc1x2vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
} 
void MomentCalc1x2vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
  out[3] += 2.0*f[17]*volFact; 
} 
void MomentCalc1x2vMax_M0_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
  out[3] += 2.0*f[17]*volFact; 
  out[4] += 2.0*f[32]*volFact; 
} 
void MomentCalc1x2vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*wx1*volFact+0.5773502691896258*f[2]*dv1*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact; 
  out[2] += 2.0*f[0]*wx2*volFact+0.5773502691896258*f[3]*dv2*volFact; 
  out[3] += 2.0*f[1]*wx2*volFact; 
} 
void MomentCalc1x2vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*wx1*volFact+0.5773502691896258*f[2]*dv1*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact+0.5773502691896258*f[4]*dv1*volFact; 
  out[2] += 2.0*f[7]*wx1*volFact; 
  out[3] += 2.0*f[0]*wx2*volFact+0.5773502691896258*f[3]*dv2*volFact; 
  out[4] += 2.0*f[1]*wx2*volFact+0.5773502691896258*f[5]*dv2*volFact; 
  out[5] += 2.0*f[7]*wx2*volFact; 
} 
void MomentCalc1x2vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*wx1*volFact+0.5773502691896258*f[2]*dv1*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact+0.5773502691896258*f[4]*dv1*volFact; 
  out[2] += 2.0*f[7]*wx1*volFact+0.5773502691896257*f[11]*dv1*volFact; 
  out[3] += 2.0*f[17]*wx1*volFact; 
  out[4] += 2.0*f[0]*wx2*volFact+0.5773502691896258*f[3]*dv2*volFact; 
  out[5] += 2.0*f[1]*wx2*volFact+0.5773502691896258*f[5]*dv2*volFact; 
  out[6] += 2.0*f[7]*wx2*volFact+0.5773502691896257*f[13]*dv2*volFact; 
  out[7] += 2.0*f[17]*wx2*volFact; 
} 
void MomentCalc1x2vMax_M1i_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*wx1*volFact+0.5773502691896258*f[2]*dv1*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact+0.5773502691896258*f[4]*dv1*volFact; 
  out[2] += 2.0*f[7]*wx1*volFact+0.5773502691896257*f[11]*dv1*volFact; 
  out[3] += 2.0*f[17]*wx1*volFact+0.5773502691896256*f[26]*dv1*volFact; 
  out[4] += 2.0*f[32]*wx1*volFact; 
  out[5] += 2.0*f[0]*wx2*volFact+0.5773502691896258*f[3]*dv2*volFact; 
  out[6] += 2.0*f[1]*wx2*volFact+0.5773502691896258*f[5]*dv2*volFact; 
  out[7] += 2.0*f[7]*wx2*volFact+0.5773502691896257*f[13]*dv2*volFact; 
  out[8] += 2.0*f[17]*wx2*volFact+0.5773502691896256*f[28]*dv2*volFact; 
  out[9] += 2.0*f[32]*wx2*volFact; 
} 
void MomentCalc1x2vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[0]*wx1*wx2*volFact+0.5773502691896258*f[2]*dv1*wx2*volFact+0.5773502691896258*f[3]*dv2*wx1*volFact; 
  out[3] += 2.0*f[1]*wx1*wx2*volFact; 
  out[4] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx2*volFact+0.1666666666666667*f[0]*dv2_sq*volFact; 
  out[5] += 2.0*f[1]*wx2_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact; 
} 
void MomentCalc1x2vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.149071198499986*f[8]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[4]*dv1*wx1*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*wx1_sq*volFact+0.1666666666666667*f[7]*dv1_sq*volFact; 
  out[3] += 2.0*f[0]*wx1*wx2*volFact+0.5773502691896258*f[2]*dv1*wx2*volFact+0.5773502691896258*f[3]*dv2*wx1*volFact+0.1666666666666667*f[6]*dv1*dv2*volFact; 
  out[4] += 2.0*f[1]*wx1*wx2*volFact+0.5773502691896258*f[4]*dv1*wx2*volFact+0.5773502691896258*f[5]*dv2*wx1*volFact; 
  out[5] += 2.0*f[7]*wx1*wx2*volFact; 
  out[6] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx2*volFact+0.149071198499986*f[9]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact; 
  out[7] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[5]*dv2*wx2*volFact+0.1666666666666667*f[1]*dv2_sq*volFact; 
  out[8] += 2.0*f[7]*wx2_sq*volFact+0.1666666666666667*f[7]*dv2_sq*volFact; 
} 
void MomentCalc1x2vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.149071198499986*f[8]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[4]*dv1*wx1*volFact+0.149071198499986*f[12]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*wx1_sq*volFact+1.154700538379251*f[11]*dv1*wx1*volFact+0.1666666666666667*f[7]*dv1_sq*volFact; 
  out[3] += 2.0*f[17]*wx1_sq*volFact+0.1666666666666667*f[17]*dv1_sq*volFact; 
  out[4] += 2.0*f[0]*wx1*wx2*volFact+0.5773502691896258*f[2]*dv1*wx2*volFact+0.5773502691896258*f[3]*dv2*wx1*volFact+0.1666666666666667*f[6]*dv1*dv2*volFact; 
  out[5] += 2.0*f[1]*wx1*wx2*volFact+0.5773502691896258*f[4]*dv1*wx2*volFact+0.5773502691896258*f[5]*dv2*wx1*volFact+0.1666666666666667*f[10]*dv1*dv2*volFact; 
  out[6] += 2.0*f[7]*wx1*wx2*volFact+0.5773502691896257*f[11]*dv1*wx2*volFact+0.5773502691896257*f[13]*dv2*wx1*volFact; 
  out[7] += 2.0*f[17]*wx1*wx2*volFact; 
  out[8] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx2*volFact+0.149071198499986*f[9]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact; 
  out[9] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[5]*dv2*wx2*volFact+0.149071198499986*f[15]*dv2_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact; 
  out[10] += 2.0*f[7]*wx2_sq*volFact+1.154700538379251*f[13]*dv2*wx2*volFact+0.1666666666666667*f[7]*dv2_sq*volFact; 
  out[11] += 2.0*f[17]*wx2_sq*volFact+0.1666666666666667*f[17]*dv2_sq*volFact; 
} 
void MomentCalc1x2vMax_M2ij_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.149071198499986*f[8]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[4]*dv1*wx1*volFact+0.149071198499986*f[12]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*wx1_sq*volFact+1.154700538379251*f[11]*dv1*wx1*volFact+0.149071198499986*f[23]*dv1_sq*volFact+0.1666666666666667*f[7]*dv1_sq*volFact; 
  out[3] += 2.0*f[17]*wx1_sq*volFact+1.154700538379251*f[26]*dv1*wx1*volFact+0.1666666666666667*f[17]*dv1_sq*volFact; 
  out[4] += 2.0*f[32]*wx1_sq*volFact+0.1666666666666667*f[32]*dv1_sq*volFact; 
  out[5] += 2.0*f[0]*wx1*wx2*volFact+0.5773502691896258*f[2]*dv1*wx2*volFact+0.5773502691896258*f[3]*dv2*wx1*volFact+0.1666666666666667*f[6]*dv1*dv2*volFact; 
  out[6] += 2.0*f[1]*wx1*wx2*volFact+0.5773502691896258*f[4]*dv1*wx2*volFact+0.5773502691896258*f[5]*dv2*wx1*volFact+0.1666666666666667*f[10]*dv1*dv2*volFact; 
  out[7] += 2.0*f[7]*wx1*wx2*volFact+0.5773502691896257*f[11]*dv1*wx2*volFact+0.5773502691896257*f[13]*dv2*wx1*volFact+0.1666666666666667*f[20]*dv1*dv2*volFact; 
  out[8] += 2.0*f[17]*wx1*wx2*volFact+0.5773502691896256*f[26]*dv1*wx2*volFact+0.5773502691896256*f[28]*dv2*wx1*volFact; 
  out[9] += 2.0*f[32]*wx1*wx2*volFact; 
  out[10] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx2*volFact+0.149071198499986*f[9]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact; 
  out[11] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[5]*dv2*wx2*volFact+0.149071198499986*f[15]*dv2_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact; 
  out[12] += 2.0*f[7]*wx2_sq*volFact+1.154700538379251*f[13]*dv2*wx2*volFact+0.149071198499986*f[24]*dv2_sq*volFact+0.1666666666666667*f[7]*dv2_sq*volFact; 
  out[13] += 2.0*f[17]*wx2_sq*volFact+1.154700538379251*f[28]*dv2*wx2*volFact+0.1666666666666667*f[17]*dv2_sq*volFact; 
  out[14] += 2.0*f[32]*wx2_sq*volFact+0.1666666666666667*f[32]*dv2_sq*volFact; 
} 
void MomentCalc1x2vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx2*volFact+2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.1666666666666667*f[0]*dv2_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx2_sq*volFact+2.0*f[1]*wx1_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
} 
void MomentCalc1x2vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx2*volFact+2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.149071198499986*f[9]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact+0.149071198499986*f[8]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[5]*dv2*wx2*volFact+2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[4]*dv1*wx1*volFact+0.1666666666666667*f[1]*dv2_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*wx2_sq*volFact+2.0*f[7]*wx1_sq*volFact+0.1666666666666667*f[7]*dv2_sq*volFact+0.1666666666666667*f[7]*dv1_sq*volFact; 
} 
void MomentCalc1x2vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx2*volFact+2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.149071198499986*f[9]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact+0.149071198499986*f[8]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[5]*dv2*wx2*volFact+2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[4]*dv1*wx1*volFact+0.149071198499986*f[15]*dv2_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact+0.149071198499986*f[12]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*wx2_sq*volFact+1.154700538379251*f[13]*dv2*wx2*volFact+2.0*f[7]*wx1_sq*volFact+1.154700538379251*f[11]*dv1*wx1*volFact+0.1666666666666667*f[7]*dv2_sq*volFact+0.1666666666666667*f[7]*dv1_sq*volFact; 
  out[3] += 2.0*f[17]*wx2_sq*volFact+2.0*f[17]*wx1_sq*volFact+0.1666666666666667*f[17]*dv2_sq*volFact+0.1666666666666667*f[17]*dv1_sq*volFact; 
} 
void MomentCalc1x2vMax_M2_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx2*volFact+2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.149071198499986*f[9]*dv2_sq*volFact+0.1666666666666667*f[0]*dv2_sq*volFact+0.149071198499986*f[8]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx2_sq*volFact+1.154700538379252*f[5]*dv2*wx2*volFact+2.0*f[1]*wx1_sq*volFact+1.154700538379252*f[4]*dv1*wx1*volFact+0.149071198499986*f[15]*dv2_sq*volFact+0.1666666666666667*f[1]*dv2_sq*volFact+0.149071198499986*f[12]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*wx2_sq*volFact+1.154700538379251*f[13]*dv2*wx2*volFact+2.0*f[7]*wx1_sq*volFact+1.154700538379251*f[11]*dv1*wx1*volFact+0.149071198499986*f[24]*dv2_sq*volFact+0.1666666666666667*f[7]*dv2_sq*volFact+0.149071198499986*f[23]*dv1_sq*volFact+0.1666666666666667*f[7]*dv1_sq*volFact; 
  out[3] += 2.0*f[17]*wx2_sq*volFact+1.154700538379251*f[28]*dv2*wx2*volFact+2.0*f[17]*wx1_sq*volFact+1.154700538379251*f[26]*dv1*wx1*volFact+0.1666666666666667*f[17]*dv2_sq*volFact+0.1666666666666667*f[17]*dv1_sq*volFact; 
  out[4] += 2.0*f[32]*wx2_sq*volFact+2.0*f[32]*wx1_sq*volFact+0.1666666666666667*f[32]*dv2_sq*volFact+0.1666666666666667*f[32]*dv1_sq*volFact; 
} 
void MomentCalc1x2vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*wx1*wx2_sq*volFact+0.5773502691896258*f[2]*dv1*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx1*wx2*volFact+2.0*f[0]*wx1*wx1_sq*volFact+1.732050807568877*f[2]*dv1*wx1_sq*volFact+0.1666666666666667*f[0]*dv2_sq*wx1*volFact+0.5*f[0]*dv1_sq*wx1*volFact+0.04811252243246882*f[2]*dv1*dv2_sq*volFact+0.08660254037844387*f[2]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1*wx2_sq*volFact+2.0*f[1]*wx1*wx1_sq*volFact+0.1666666666666667*f[1]*dv2_sq*wx1*volFact+0.5*f[1]*dv1_sq*wx1*volFact; 
  out[2] += 2.0*f[0]*wx2*wx2_sq*volFact+1.732050807568877*f[3]*dv2*wx2_sq*volFact+2.0*f[0]*wx1_sq*wx2*volFact+1.154700538379252*f[2]*dv1*wx1*wx2*volFact+0.5*f[0]*dv2_sq*wx2*volFact+0.1666666666666667*f[0]*dv1_sq*wx2*volFact+0.5773502691896258*f[3]*dv2*wx1_sq*volFact+0.08660254037844387*f[3]*dv2*dv2_sq*volFact+0.04811252243246882*f[3]*dv1_sq*dv2*volFact; 
  out[3] += 2.0*f[1]*wx2*wx2_sq*volFact+2.0*f[1]*wx1_sq*wx2*volFact+0.5*f[1]*dv2_sq*wx2*volFact+0.1666666666666667*f[1]*dv1_sq*wx2*volFact; 
} 
void MomentCalc1x2vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*wx1*wx2_sq*volFact+0.5773502691896258*f[2]*dv1*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx1*wx2*volFact+0.3333333333333333*f[6]*dv1*dv2*wx2*volFact+2.0*f[0]*wx1*wx1_sq*volFact+1.732050807568877*f[2]*dv1*wx1_sq*volFact+0.149071198499986*f[9]*dv2_sq*wx1*volFact+0.1666666666666667*f[0]*dv2_sq*wx1*volFact+0.4472135954999579*f[8]*dv1_sq*wx1*volFact+0.5*f[0]*dv1_sq*wx1*volFact+0.04811252243246882*f[2]*dv1*dv2_sq*volFact+0.08660254037844387*f[2]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1*wx2_sq*volFact+0.5773502691896258*f[4]*dv1*wx2_sq*volFact+1.154700538379252*f[5]*dv2*wx1*wx2*volFact+2.0*f[1]*wx1*wx1_sq*volFact+1.732050807568877*f[4]*dv1*wx1_sq*volFact+0.1666666666666667*f[1]*dv2_sq*wx1*volFact+0.5*f[1]*dv1_sq*wx1*volFact+0.04811252243246882*f[4]*dv1*dv2_sq*volFact+0.08660254037844387*f[4]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*wx1*wx2_sq*volFact+2.0*f[7]*wx1*wx1_sq*volFact+0.1666666666666667*f[7]*dv2_sq*wx1*volFact+0.5*f[7]*dv1_sq*wx1*volFact; 
  out[3] += 2.0*f[0]*wx2*wx2_sq*volFact+1.732050807568877*f[3]*dv2*wx2_sq*volFact+2.0*f[0]*wx1_sq*wx2*volFact+1.154700538379252*f[2]*dv1*wx1*wx2*volFact+0.4472135954999579*f[9]*dv2_sq*wx2*volFact+0.5*f[0]*dv2_sq*wx2*volFact+0.149071198499986*f[8]*dv1_sq*wx2*volFact+0.1666666666666667*f[0]*dv1_sq*wx2*volFact+0.5773502691896258*f[3]*dv2*wx1_sq*volFact+0.3333333333333333*f[6]*dv1*dv2*wx1*volFact+0.08660254037844387*f[3]*dv2*dv2_sq*volFact+0.04811252243246882*f[3]*dv1_sq*dv2*volFact; 
  out[4] += 2.0*f[1]*wx2*wx2_sq*volFact+1.732050807568877*f[5]*dv2*wx2_sq*volFact+2.0*f[1]*wx1_sq*wx2*volFact+1.154700538379252*f[4]*dv1*wx1*wx2*volFact+0.5*f[1]*dv2_sq*wx2*volFact+0.1666666666666667*f[1]*dv1_sq*wx2*volFact+0.5773502691896258*f[5]*dv2*wx1_sq*volFact+0.08660254037844387*f[5]*dv2*dv2_sq*volFact+0.04811252243246882*f[5]*dv1_sq*dv2*volFact; 
  out[5] += 2.0*f[7]*wx2*wx2_sq*volFact+2.0*f[7]*wx1_sq*wx2*volFact+0.5*f[7]*dv2_sq*wx2*volFact+0.1666666666666667*f[7]*dv1_sq*wx2*volFact; 
} 
void MomentCalc1x2vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*wx1*wx2_sq*volFact+0.5773502691896258*f[2]*dv1*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx1*wx2*volFact+0.3333333333333333*f[6]*dv1*dv2*wx2*volFact+2.0*f[0]*wx1*wx1_sq*volFact+1.732050807568877*f[2]*dv1*wx1_sq*volFact+0.149071198499986*f[9]*dv2_sq*wx1*volFact+0.1666666666666667*f[0]*dv2_sq*wx1*volFact+0.4472135954999579*f[8]*dv1_sq*wx1*volFact+0.5*f[0]*dv1_sq*wx1*volFact+0.04303314829119351*f[16]*dv1*dv2_sq*volFact+0.04811252243246882*f[2]*dv1*dv2_sq*volFact+0.03779644730092272*f[18]*dv1*dv1_sq*volFact+0.08660254037844387*f[2]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1*wx2_sq*volFact+0.5773502691896258*f[4]*dv1*wx2_sq*volFact+1.154700538379252*f[5]*dv2*wx1*wx2*volFact+0.3333333333333333*f[10]*dv1*dv2*wx2*volFact+2.0*f[1]*wx1*wx1_sq*volFact+1.732050807568877*f[4]*dv1*wx1_sq*volFact+0.149071198499986*f[15]*dv2_sq*wx1*volFact+0.1666666666666667*f[1]*dv2_sq*wx1*volFact+0.447213595499958*f[12]*dv1_sq*wx1*volFact+0.5*f[1]*dv1_sq*wx1*volFact+0.04811252243246882*f[4]*dv1*dv2_sq*volFact+0.08660254037844387*f[4]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*wx1*wx2_sq*volFact+0.5773502691896257*f[11]*dv1*wx2_sq*volFact+1.154700538379251*f[13]*dv2*wx1*wx2*volFact+2.0*f[7]*wx1*wx1_sq*volFact+1.732050807568877*f[11]*dv1*wx1_sq*volFact+0.1666666666666667*f[7]*dv2_sq*wx1*volFact+0.5*f[7]*dv1_sq*wx1*volFact+0.04811252243246881*f[11]*dv1*dv2_sq*volFact+0.08660254037844385*f[11]*dv1*dv1_sq*volFact; 
  out[3] += 2.0*f[17]*wx1*wx2_sq*volFact+2.0*f[17]*wx1*wx1_sq*volFact+0.1666666666666667*f[17]*dv2_sq*wx1*volFact+0.5*f[17]*dv1_sq*wx1*volFact; 
  out[4] += 2.0*f[0]*wx2*wx2_sq*volFact+1.732050807568877*f[3]*dv2*wx2_sq*volFact+2.0*f[0]*wx1_sq*wx2*volFact+1.154700538379252*f[2]*dv1*wx1*wx2*volFact+0.4472135954999579*f[9]*dv2_sq*wx2*volFact+0.5*f[0]*dv2_sq*wx2*volFact+0.149071198499986*f[8]*dv1_sq*wx2*volFact+0.1666666666666667*f[0]*dv1_sq*wx2*volFact+0.5773502691896258*f[3]*dv2*wx1_sq*volFact+0.3333333333333333*f[6]*dv1*dv2*wx1*volFact+0.03779644730092272*f[19]*dv2*dv2_sq*volFact+0.08660254037844387*f[3]*dv2*dv2_sq*volFact+0.04303314829119351*f[14]*dv1_sq*dv2*volFact+0.04811252243246882*f[3]*dv1_sq*dv2*volFact; 
  out[5] += 2.0*f[1]*wx2*wx2_sq*volFact+1.732050807568877*f[5]*dv2*wx2_sq*volFact+2.0*f[1]*wx1_sq*wx2*volFact+1.154700538379252*f[4]*dv1*wx1*wx2*volFact+0.447213595499958*f[15]*dv2_sq*wx2*volFact+0.5*f[1]*dv2_sq*wx2*volFact+0.149071198499986*f[12]*dv1_sq*wx2*volFact+0.1666666666666667*f[1]*dv1_sq*wx2*volFact+0.5773502691896258*f[5]*dv2*wx1_sq*volFact+0.3333333333333333*f[10]*dv1*dv2*wx1*volFact+0.08660254037844387*f[5]*dv2*dv2_sq*volFact+0.04811252243246882*f[5]*dv1_sq*dv2*volFact; 
  out[6] += 2.0*f[7]*wx2*wx2_sq*volFact+1.732050807568877*f[13]*dv2*wx2_sq*volFact+2.0*f[7]*wx1_sq*wx2*volFact+1.154700538379251*f[11]*dv1*wx1*wx2*volFact+0.5*f[7]*dv2_sq*wx2*volFact+0.1666666666666667*f[7]*dv1_sq*wx2*volFact+0.5773502691896257*f[13]*dv2*wx1_sq*volFact+0.08660254037844385*f[13]*dv2*dv2_sq*volFact+0.04811252243246881*f[13]*dv1_sq*dv2*volFact; 
  out[7] += 2.0*f[17]*wx2*wx2_sq*volFact+2.0*f[17]*wx1_sq*wx2*volFact+0.5*f[17]*dv2_sq*wx2*volFact+0.1666666666666667*f[17]*dv1_sq*wx2*volFact; 
} 
void MomentCalc1x2vMax_M3i_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*wx1*wx2_sq*volFact+0.5773502691896258*f[2]*dv1*wx2_sq*volFact+1.154700538379252*f[3]*dv2*wx1*wx2*volFact+0.3333333333333333*f[6]*dv1*dv2*wx2*volFact+2.0*f[0]*wx1*wx1_sq*volFact+1.732050807568877*f[2]*dv1*wx1_sq*volFact+0.149071198499986*f[9]*dv2_sq*wx1*volFact+0.1666666666666667*f[0]*dv2_sq*wx1*volFact+0.4472135954999579*f[8]*dv1_sq*wx1*volFact+0.5*f[0]*dv1_sq*wx1*volFact+0.04303314829119351*f[16]*dv1*dv2_sq*volFact+0.04811252243246882*f[2]*dv1*dv2_sq*volFact+0.03779644730092272*f[18]*dv1*dv1_sq*volFact+0.08660254037844387*f[2]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*wx1*wx2_sq*volFact+0.5773502691896258*f[4]*dv1*wx2_sq*volFact+1.154700538379252*f[5]*dv2*wx1*wx2*volFact+0.3333333333333333*f[10]*dv1*dv2*wx2*volFact+2.0*f[1]*wx1*wx1_sq*volFact+1.732050807568877*f[4]*dv1*wx1_sq*volFact+0.149071198499986*f[15]*dv2_sq*wx1*volFact+0.1666666666666667*f[1]*dv2_sq*wx1*volFact+0.447213595499958*f[12]*dv1_sq*wx1*volFact+0.5*f[1]*dv1_sq*wx1*volFact+0.04303314829119353*f[22]*dv1*dv2_sq*volFact+0.04811252243246882*f[4]*dv1*dv2_sq*volFact+0.03779644730092273*f[27]*dv1*dv1_sq*volFact+0.08660254037844387*f[4]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*wx1*wx2_sq*volFact+0.5773502691896257*f[11]*dv1*wx2_sq*volFact+1.154700538379251*f[13]*dv2*wx1*wx2*volFact+0.3333333333333333*f[20]*dv1*dv2*wx2*volFact+2.0*f[7]*wx1*wx1_sq*volFact+1.732050807568877*f[11]*dv1*wx1_sq*volFact+0.149071198499986*f[24]*dv2_sq*wx1*volFact+0.1666666666666667*f[7]*dv2_sq*wx1*volFact+0.4472135954999579*f[23]*dv1_sq*wx1*volFact+0.5*f[7]*dv1_sq*wx1*volFact+0.04811252243246881*f[11]*dv1*dv2_sq*volFact+0.08660254037844385*f[11]*dv1*dv1_sq*volFact; 
  out[3] += 2.0*f[17]*wx1*wx2_sq*volFact+0.5773502691896256*f[26]*dv1*wx2_sq*volFact+1.154700538379251*f[28]*dv2*wx1*wx2*volFact+2.0*f[17]*wx1*wx1_sq*volFact+1.732050807568877*f[26]*dv1*wx1_sq*volFact+0.1666666666666667*f[17]*dv2_sq*wx1*volFact+0.5*f[17]*dv1_sq*wx1*volFact+0.0481125224324688*f[26]*dv1*dv2_sq*volFact+0.08660254037844385*f[26]*dv1*dv1_sq*volFact; 
  out[4] += 2.0*f[32]*wx1*wx2_sq*volFact+2.0*f[32]*wx1*wx1_sq*volFact+0.1666666666666667*f[32]*dv2_sq*wx1*volFact+0.5*f[32]*dv1_sq*wx1*volFact; 
  out[5] += 2.0*f[0]*wx2*wx2_sq*volFact+1.732050807568877*f[3]*dv2*wx2_sq*volFact+2.0*f[0]*wx1_sq*wx2*volFact+1.154700538379252*f[2]*dv1*wx1*wx2*volFact+0.4472135954999579*f[9]*dv2_sq*wx2*volFact+0.5*f[0]*dv2_sq*wx2*volFact+0.149071198499986*f[8]*dv1_sq*wx2*volFact+0.1666666666666667*f[0]*dv1_sq*wx2*volFact+0.5773502691896258*f[3]*dv2*wx1_sq*volFact+0.3333333333333333*f[6]*dv1*dv2*wx1*volFact+0.03779644730092272*f[19]*dv2*dv2_sq*volFact+0.08660254037844387*f[3]*dv2*dv2_sq*volFact+0.04303314829119351*f[14]*dv1_sq*dv2*volFact+0.04811252243246882*f[3]*dv1_sq*dv2*volFact; 
  out[6] += 2.0*f[1]*wx2*wx2_sq*volFact+1.732050807568877*f[5]*dv2*wx2_sq*volFact+2.0*f[1]*wx1_sq*wx2*volFact+1.154700538379252*f[4]*dv1*wx1*wx2*volFact+0.447213595499958*f[15]*dv2_sq*wx2*volFact+0.5*f[1]*dv2_sq*wx2*volFact+0.149071198499986*f[12]*dv1_sq*wx2*volFact+0.1666666666666667*f[1]*dv1_sq*wx2*volFact+0.5773502691896258*f[5]*dv2*wx1_sq*volFact+0.3333333333333333*f[10]*dv1*dv2*wx1*volFact+0.03779644730092273*f[30]*dv2*dv2_sq*volFact+0.08660254037844387*f[5]*dv2*dv2_sq*volFact+0.04303314829119353*f[21]*dv1_sq*dv2*volFact+0.04811252243246882*f[5]*dv1_sq*dv2*volFact; 
  out[7] += 2.0*f[7]*wx2*wx2_sq*volFact+1.732050807568877*f[13]*dv2*wx2_sq*volFact+2.0*f[7]*wx1_sq*wx2*volFact+1.154700538379251*f[11]*dv1*wx1*wx2*volFact+0.4472135954999579*f[24]*dv2_sq*wx2*volFact+0.5*f[7]*dv2_sq*wx2*volFact+0.149071198499986*f[23]*dv1_sq*wx2*volFact+0.1666666666666667*f[7]*dv1_sq*wx2*volFact+0.5773502691896257*f[13]*dv2*wx1_sq*volFact+0.3333333333333333*f[20]*dv1*dv2*wx1*volFact+0.08660254037844385*f[13]*dv2*dv2_sq*volFact+0.04811252243246881*f[13]*dv1_sq*dv2*volFact; 
  out[8] += 2.0*f[17]*wx2*wx2_sq*volFact+1.732050807568877*f[28]*dv2*wx2_sq*volFact+2.0*f[17]*wx1_sq*wx2*volFact+1.154700538379251*f[26]*dv1*wx1*wx2*volFact+0.5*f[17]*dv2_sq*wx2*volFact+0.1666666666666667*f[17]*dv1_sq*wx2*volFact+0.5773502691896256*f[28]*dv2*wx1_sq*volFact+0.08660254037844385*f[28]*dv2*dv2_sq*volFact+0.0481125224324688*f[28]*dv1_sq*dv2*volFact; 
  out[9] += 2.0*f[32]*wx2*wx2_sq*volFact+2.0*f[32]*wx1_sq*wx2*volFact+0.5*f[32]*dv2_sq*wx2*volFact+0.1666666666666667*f[32]*dv1_sq*wx2*volFact; 
} 
