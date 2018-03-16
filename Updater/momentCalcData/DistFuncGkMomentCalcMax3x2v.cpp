#include <math.h> 
#include <DistFuncMomentCalcModDecl.h> 
void GkMomentCalc3x2vMax_Dens_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[3]*volFact; 
} 
void GkMomentCalc3x2vMax_Dens_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[3]*volFact; 
  out[4] += 2.0*f[6]*volFact; 
  out[5] += 2.0*f[7]*volFact; 
  out[6] += 2.0*f[8]*volFact; 
  out[7] += 2.0*f[16]*volFact; 
  out[8] += 2.0*f[17]*volFact; 
  out[9] += 2.0*f[18]*volFact; 
} 
void GkMomentCalc3x2vMax_Dens_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[3]*volFact; 
  out[4] += 2.0*f[6]*volFact; 
  out[5] += 2.0*f[7]*volFact; 
  out[6] += 2.0*f[8]*volFact; 
  out[7] += 2.0*f[16]*volFact; 
  out[8] += 2.0*f[17]*volFact; 
  out[9] += 2.0*f[18]*volFact; 
  out[10] += 2.0*f[21]*volFact; 
  out[11] += 2.0*f[31]*volFact; 
  out[12] += 2.0*f[32]*volFact; 
  out[13] += 2.0*f[33]*volFact; 
  out[14] += 2.0*f[34]*volFact; 
  out[15] += 2.0*f[35]*volFact; 
  out[16] += 2.0*f[36]*volFact; 
  out[17] += 2.0*f[51]*volFact; 
  out[18] += 2.0*f[52]*volFact; 
  out[19] += 2.0*f[53]*volFact; 
} 
void GkMomentCalc3x2vMax_Dens_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[3]*volFact; 
  out[4] += 2.0*f[6]*volFact; 
  out[5] += 2.0*f[7]*volFact; 
  out[6] += 2.0*f[8]*volFact; 
  out[7] += 2.0*f[16]*volFact; 
  out[8] += 2.0*f[17]*volFact; 
  out[9] += 2.0*f[18]*volFact; 
  out[10] += 2.0*f[21]*volFact; 
  out[11] += 2.0*f[31]*volFact; 
  out[12] += 2.0*f[32]*volFact; 
  out[13] += 2.0*f[33]*volFact; 
  out[14] += 2.0*f[34]*volFact; 
  out[15] += 2.0*f[35]*volFact; 
  out[16] += 2.0*f[36]*volFact; 
  out[17] += 2.0*f[51]*volFact; 
  out[18] += 2.0*f[52]*volFact; 
  out[19] += 2.0*f[53]*volFact; 
  out[20] += 2.0*f[61]*volFact; 
  out[21] += 2.0*f[62]*volFact; 
  out[22] += 2.0*f[63]*volFact; 
  out[23] += 2.0*f[91]*volFact; 
  out[24] += 2.0*f[92]*volFact; 
  out[25] += 2.0*f[93]*volFact; 
  out[26] += 2.0*f[101]*volFact; 
  out[27] += 2.0*f[102]*volFact; 
  out[28] += 2.0*f[103]*volFact; 
  out[29] += 2.0*f[104]*volFact; 
  out[30] += 2.0*f[105]*volFact; 
  out[31] += 2.0*f[106]*volFact; 
  out[32] += 2.0*f[121]*volFact; 
  out[33] += 2.0*f[122]*volFact; 
  out[34] += 2.0*f[123]*volFact; 
} 
void GkMomentCalc3x2vMax_Upar_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx1+0.5773502691896258*f[4]*dv1*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1; 
  out[2] += 2.0*f[2]*volFact*wx1; 
  out[3] += 2.0*f[3]*volFact*wx1; 
} 
void GkMomentCalc3x2vMax_Upar_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx1+0.5773502691896258*f[4]*dv1*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1+0.5773502691896258*f[9]*dv1*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1+0.5773502691896258*f[10]*dv1*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1+0.5773502691896258*f[11]*dv1*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1; 
  out[5] += 2.0*f[7]*volFact*wx1; 
  out[6] += 2.0*f[8]*volFact*wx1; 
  out[7] += 2.0*f[16]*volFact*wx1; 
  out[8] += 2.0*f[17]*volFact*wx1; 
  out[9] += 2.0*f[18]*volFact*wx1; 
} 
void GkMomentCalc3x2vMax_Upar_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx1+0.5773502691896258*f[4]*dv1*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1+0.5773502691896258*f[9]*dv1*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1+0.5773502691896258*f[10]*dv1*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1+0.5773502691896258*f[11]*dv1*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1+0.5773502691896258*f[22]*dv1*volFact; 
  out[5] += 2.0*f[7]*volFact*wx1+0.5773502691896258*f[23]*dv1*volFact; 
  out[6] += 2.0*f[8]*volFact*wx1+0.5773502691896258*f[24]*dv1*volFact; 
  out[7] += 2.0*f[16]*volFact*wx1+0.5773502691896257*f[37]*dv1*volFact; 
  out[8] += 2.0*f[17]*volFact*wx1+0.5773502691896257*f[38]*dv1*volFact; 
  out[9] += 2.0*f[18]*volFact*wx1+0.5773502691896257*f[39]*dv1*volFact; 
  out[10] += 2.0*f[21]*volFact*wx1; 
  out[11] += 2.0*f[31]*volFact*wx1; 
  out[12] += 2.0*f[32]*volFact*wx1; 
  out[13] += 2.0*f[33]*volFact*wx1; 
  out[14] += 2.0*f[34]*volFact*wx1; 
  out[15] += 2.0*f[35]*volFact*wx1; 
  out[16] += 2.0*f[36]*volFact*wx1; 
  out[17] += 2.0*f[51]*volFact*wx1; 
  out[18] += 2.0*f[52]*volFact*wx1; 
  out[19] += 2.0*f[53]*volFact*wx1; 
} 
void GkMomentCalc3x2vMax_Upar_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx1+0.5773502691896258*f[4]*dv1*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1+0.5773502691896258*f[9]*dv1*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1+0.5773502691896258*f[10]*dv1*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1+0.5773502691896258*f[11]*dv1*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1+0.5773502691896258*f[22]*dv1*volFact; 
  out[5] += 2.0*f[7]*volFact*wx1+0.5773502691896258*f[23]*dv1*volFact; 
  out[6] += 2.0*f[8]*volFact*wx1+0.5773502691896258*f[24]*dv1*volFact; 
  out[7] += 2.0*f[16]*volFact*wx1+0.5773502691896257*f[37]*dv1*volFact; 
  out[8] += 2.0*f[17]*volFact*wx1+0.5773502691896257*f[38]*dv1*volFact; 
  out[9] += 2.0*f[18]*volFact*wx1+0.5773502691896257*f[39]*dv1*volFact; 
  out[10] += 2.0*f[21]*volFact*wx1+0.5773502691896258*f[56]*dv1*volFact; 
  out[11] += 2.0*f[31]*volFact*wx1+0.5773502691896257*f[64]*dv1*volFact; 
  out[12] += 2.0*f[32]*volFact*wx1+0.5773502691896257*f[65]*dv1*volFact; 
  out[13] += 2.0*f[33]*volFact*wx1+0.5773502691896257*f[66]*dv1*volFact; 
  out[14] += 2.0*f[34]*volFact*wx1+0.5773502691896257*f[67]*dv1*volFact; 
  out[15] += 2.0*f[35]*volFact*wx1+0.5773502691896257*f[68]*dv1*volFact; 
  out[16] += 2.0*f[36]*volFact*wx1+0.5773502691896257*f[69]*dv1*volFact; 
  out[17] += 2.0*f[51]*volFact*wx1+0.5773502691896256*f[107]*dv1*volFact; 
  out[18] += 2.0*f[52]*volFact*wx1+0.5773502691896256*f[108]*dv1*volFact; 
  out[19] += 2.0*f[53]*volFact*wx1+0.5773502691896256*f[109]*dv1*volFact; 
  out[20] += 2.0*f[61]*volFact*wx1; 
  out[21] += 2.0*f[62]*volFact*wx1; 
  out[22] += 2.0*f[63]*volFact*wx1; 
  out[23] += 2.0*f[91]*volFact*wx1; 
  out[24] += 2.0*f[92]*volFact*wx1; 
  out[25] += 2.0*f[93]*volFact*wx1; 
  out[26] += 2.0*f[101]*volFact*wx1; 
  out[27] += 2.0*f[102]*volFact*wx1; 
  out[28] += 2.0*f[103]*volFact*wx1; 
  out[29] += 2.0*f[104]*volFact*wx1; 
  out[30] += 2.0*f[105]*volFact*wx1; 
  out[31] += 2.0*f[106]*volFact*wx1; 
  out[32] += 2.0*f[121]*volFact*wx1; 
  out[33] += 2.0*f[122]*volFact*wx1; 
  out[34] += 2.0*f[123]*volFact*wx1; 
} 
void GkMomentCalc3x2vMax_Ppar_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[4]*dv1*volFact*wx1+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1_sq+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1_sq+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1_sq+0.1666666666666667*f[3]*dv1_sq*volFact; 
} 
void GkMomentCalc3x2vMax_Ppar_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[4]*dv1*volFact*wx1+0.149071198499986*f[19]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1_sq+1.154700538379252*f[9]*dv1*volFact*wx1+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1_sq+1.154700538379252*f[10]*dv1*volFact*wx1+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1_sq+1.154700538379252*f[11]*dv1*volFact*wx1+0.1666666666666667*f[3]*dv1_sq*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1_sq+0.1666666666666667*f[6]*dv1_sq*volFact; 
  out[5] += 2.0*f[7]*volFact*wx1_sq+0.1666666666666667*f[7]*dv1_sq*volFact; 
  out[6] += 2.0*f[8]*volFact*wx1_sq+0.1666666666666667*f[8]*dv1_sq*volFact; 
  out[7] += 2.0*f[16]*volFact*wx1_sq+0.1666666666666667*f[16]*dv1_sq*volFact; 
  out[8] += 2.0*f[17]*volFact*wx1_sq+0.1666666666666667*f[17]*dv1_sq*volFact; 
  out[9] += 2.0*f[18]*volFact*wx1_sq+0.1666666666666667*f[18]*dv1_sq*volFact; 
} 
void GkMomentCalc3x2vMax_Ppar_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[4]*dv1*volFact*wx1+0.149071198499986*f[19]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1_sq+1.154700538379252*f[9]*dv1*volFact*wx1+0.149071198499986*f[40]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1_sq+1.154700538379252*f[10]*dv1*volFact*wx1+0.149071198499986*f[41]*dv1_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1_sq+1.154700538379252*f[11]*dv1*volFact*wx1+0.149071198499986*f[42]*dv1_sq*volFact+0.1666666666666667*f[3]*dv1_sq*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1_sq+1.154700538379252*f[22]*dv1*volFact*wx1+0.1666666666666667*f[6]*dv1_sq*volFact; 
  out[5] += 2.0*f[7]*volFact*wx1_sq+1.154700538379252*f[23]*dv1*volFact*wx1+0.1666666666666667*f[7]*dv1_sq*volFact; 
  out[6] += 2.0*f[8]*volFact*wx1_sq+1.154700538379252*f[24]*dv1*volFact*wx1+0.1666666666666667*f[8]*dv1_sq*volFact; 
  out[7] += 2.0*f[16]*volFact*wx1_sq+1.154700538379251*f[37]*dv1*volFact*wx1+0.1666666666666667*f[16]*dv1_sq*volFact; 
  out[8] += 2.0*f[17]*volFact*wx1_sq+1.154700538379251*f[38]*dv1*volFact*wx1+0.1666666666666667*f[17]*dv1_sq*volFact; 
  out[9] += 2.0*f[18]*volFact*wx1_sq+1.154700538379251*f[39]*dv1*volFact*wx1+0.1666666666666667*f[18]*dv1_sq*volFact; 
  out[10] += 2.0*f[21]*volFact*wx1_sq+0.1666666666666667*f[21]*dv1_sq*volFact; 
  out[11] += 2.0*f[31]*volFact*wx1_sq+0.1666666666666667*f[31]*dv1_sq*volFact; 
  out[12] += 2.0*f[32]*volFact*wx1_sq+0.1666666666666667*f[32]*dv1_sq*volFact; 
  out[13] += 2.0*f[33]*volFact*wx1_sq+0.1666666666666667*f[33]*dv1_sq*volFact; 
  out[14] += 2.0*f[34]*volFact*wx1_sq+0.1666666666666667*f[34]*dv1_sq*volFact; 
  out[15] += 2.0*f[35]*volFact*wx1_sq+0.1666666666666667*f[35]*dv1_sq*volFact; 
  out[16] += 2.0*f[36]*volFact*wx1_sq+0.1666666666666667*f[36]*dv1_sq*volFact; 
  out[17] += 2.0*f[51]*volFact*wx1_sq+0.1666666666666667*f[51]*dv1_sq*volFact; 
  out[18] += 2.0*f[52]*volFact*wx1_sq+0.1666666666666667*f[52]*dv1_sq*volFact; 
  out[19] += 2.0*f[53]*volFact*wx1_sq+0.1666666666666667*f[53]*dv1_sq*volFact; 
} 
void GkMomentCalc3x2vMax_Ppar_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[4]*dv1*volFact*wx1+0.149071198499986*f[19]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1_sq+1.154700538379252*f[9]*dv1*volFact*wx1+0.149071198499986*f[40]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1_sq+1.154700538379252*f[10]*dv1*volFact*wx1+0.149071198499986*f[41]*dv1_sq*volFact+0.1666666666666667*f[2]*dv1_sq*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1_sq+1.154700538379252*f[11]*dv1*volFact*wx1+0.149071198499986*f[42]*dv1_sq*volFact+0.1666666666666667*f[3]*dv1_sq*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1_sq+1.154700538379252*f[22]*dv1*volFact*wx1+0.149071198499986*f[70]*dv1_sq*volFact+0.1666666666666667*f[6]*dv1_sq*volFact; 
  out[5] += 2.0*f[7]*volFact*wx1_sq+1.154700538379252*f[23]*dv1*volFact*wx1+0.149071198499986*f[71]*dv1_sq*volFact+0.1666666666666667*f[7]*dv1_sq*volFact; 
  out[6] += 2.0*f[8]*volFact*wx1_sq+1.154700538379252*f[24]*dv1*volFact*wx1+0.149071198499986*f[72]*dv1_sq*volFact+0.1666666666666667*f[8]*dv1_sq*volFact; 
  out[7] += 2.0*f[16]*volFact*wx1_sq+1.154700538379251*f[37]*dv1*volFact*wx1+0.149071198499986*f[94]*dv1_sq*volFact+0.1666666666666667*f[16]*dv1_sq*volFact; 
  out[8] += 2.0*f[17]*volFact*wx1_sq+1.154700538379251*f[38]*dv1*volFact*wx1+0.149071198499986*f[95]*dv1_sq*volFact+0.1666666666666667*f[17]*dv1_sq*volFact; 
  out[9] += 2.0*f[18]*volFact*wx1_sq+1.154700538379251*f[39]*dv1*volFact*wx1+0.149071198499986*f[96]*dv1_sq*volFact+0.1666666666666667*f[18]*dv1_sq*volFact; 
  out[10] += 2.0*f[21]*volFact*wx1_sq+1.154700538379252*f[56]*dv1*volFact*wx1+0.1666666666666667*f[21]*dv1_sq*volFact; 
  out[11] += 2.0*f[31]*volFact*wx1_sq+1.154700538379251*f[64]*dv1*volFact*wx1+0.1666666666666667*f[31]*dv1_sq*volFact; 
  out[12] += 2.0*f[32]*volFact*wx1_sq+1.154700538379251*f[65]*dv1*volFact*wx1+0.1666666666666667*f[32]*dv1_sq*volFact; 
  out[13] += 2.0*f[33]*volFact*wx1_sq+1.154700538379251*f[66]*dv1*volFact*wx1+0.1666666666666667*f[33]*dv1_sq*volFact; 
  out[14] += 2.0*f[34]*volFact*wx1_sq+1.154700538379251*f[67]*dv1*volFact*wx1+0.1666666666666667*f[34]*dv1_sq*volFact; 
  out[15] += 2.0*f[35]*volFact*wx1_sq+1.154700538379251*f[68]*dv1*volFact*wx1+0.1666666666666667*f[35]*dv1_sq*volFact; 
  out[16] += 2.0*f[36]*volFact*wx1_sq+1.154700538379251*f[69]*dv1*volFact*wx1+0.1666666666666667*f[36]*dv1_sq*volFact; 
  out[17] += 2.0*f[51]*volFact*wx1_sq+1.154700538379251*f[107]*dv1*volFact*wx1+0.1666666666666667*f[51]*dv1_sq*volFact; 
  out[18] += 2.0*f[52]*volFact*wx1_sq+1.154700538379251*f[108]*dv1*volFact*wx1+0.1666666666666667*f[52]*dv1_sq*volFact; 
  out[19] += 2.0*f[53]*volFact*wx1_sq+1.154700538379251*f[109]*dv1*volFact*wx1+0.1666666666666667*f[53]*dv1_sq*volFact; 
  out[20] += 2.0*f[61]*volFact*wx1_sq+0.1666666666666667*f[61]*dv1_sq*volFact; 
  out[21] += 2.0*f[62]*volFact*wx1_sq+0.1666666666666667*f[62]*dv1_sq*volFact; 
  out[22] += 2.0*f[63]*volFact*wx1_sq+0.1666666666666667*f[63]*dv1_sq*volFact; 
  out[23] += 2.0*f[91]*volFact*wx1_sq+0.1666666666666667*f[91]*dv1_sq*volFact; 
  out[24] += 2.0*f[92]*volFact*wx1_sq+0.1666666666666667*f[92]*dv1_sq*volFact; 
  out[25] += 2.0*f[93]*volFact*wx1_sq+0.1666666666666667*f[93]*dv1_sq*volFact; 
  out[26] += 2.0*f[101]*volFact*wx1_sq+0.1666666666666667*f[101]*dv1_sq*volFact; 
  out[27] += 2.0*f[102]*volFact*wx1_sq+0.1666666666666667*f[102]*dv1_sq*volFact; 
  out[28] += 2.0*f[103]*volFact*wx1_sq+0.1666666666666667*f[103]*dv1_sq*volFact; 
  out[29] += 2.0*f[104]*volFact*wx1_sq+0.1666666666666667*f[104]*dv1_sq*volFact; 
  out[30] += 2.0*f[105]*volFact*wx1_sq+0.1666666666666667*f[105]*dv1_sq*volFact; 
  out[31] += 2.0*f[106]*volFact*wx1_sq+0.1666666666666667*f[106]*dv1_sq*volFact; 
  out[32] += 2.0*f[121]*volFact*wx1_sq+0.1666666666666667*f[121]*dv1_sq*volFact; 
  out[33] += 2.0*f[122]*volFact*wx1_sq+0.1666666666666667*f[122]*dv1_sq*volFact; 
  out[34] += 2.0*f[123]*volFact*wx1_sq+0.1666666666666667*f[123]*dv1_sq*volFact; 
} 
void GkMomentCalc3x2vMax_Pperp_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx2; 
  out[2] += 2.0*f[2]*volFact*wx2; 
  out[3] += 2.0*f[3]*volFact*wx2; 
} 
void GkMomentCalc3x2vMax_Pperp_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx2+0.5773502691896258*f[12]*dv2*volFact; 
  out[2] += 2.0*f[2]*volFact*wx2+0.5773502691896258*f[13]*dv2*volFact; 
  out[3] += 2.0*f[3]*volFact*wx2+0.5773502691896258*f[14]*dv2*volFact; 
  out[4] += 2.0*f[6]*volFact*wx2; 
  out[5] += 2.0*f[7]*volFact*wx2; 
  out[6] += 2.0*f[8]*volFact*wx2; 
  out[7] += 2.0*f[16]*volFact*wx2; 
  out[8] += 2.0*f[17]*volFact*wx2; 
  out[9] += 2.0*f[18]*volFact*wx2; 
} 
void GkMomentCalc3x2vMax_Pperp_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx2+0.5773502691896258*f[12]*dv2*volFact; 
  out[2] += 2.0*f[2]*volFact*wx2+0.5773502691896258*f[13]*dv2*volFact; 
  out[3] += 2.0*f[3]*volFact*wx2+0.5773502691896258*f[14]*dv2*volFact; 
  out[4] += 2.0*f[6]*volFact*wx2+0.5773502691896258*f[25]*dv2*volFact; 
  out[5] += 2.0*f[7]*volFact*wx2+0.5773502691896258*f[26]*dv2*volFact; 
  out[6] += 2.0*f[8]*volFact*wx2+0.5773502691896258*f[27]*dv2*volFact; 
  out[7] += 2.0*f[16]*volFact*wx2+0.5773502691896257*f[43]*dv2*volFact; 
  out[8] += 2.0*f[17]*volFact*wx2+0.5773502691896257*f[44]*dv2*volFact; 
  out[9] += 2.0*f[18]*volFact*wx2+0.5773502691896257*f[45]*dv2*volFact; 
  out[10] += 2.0*f[21]*volFact*wx2; 
  out[11] += 2.0*f[31]*volFact*wx2; 
  out[12] += 2.0*f[32]*volFact*wx2; 
  out[13] += 2.0*f[33]*volFact*wx2; 
  out[14] += 2.0*f[34]*volFact*wx2; 
  out[15] += 2.0*f[35]*volFact*wx2; 
  out[16] += 2.0*f[36]*volFact*wx2; 
  out[17] += 2.0*f[51]*volFact*wx2; 
  out[18] += 2.0*f[52]*volFact*wx2; 
  out[19] += 2.0*f[53]*volFact*wx2; 
} 
void GkMomentCalc3x2vMax_Pperp_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx2+0.5773502691896258*f[12]*dv2*volFact; 
  out[2] += 2.0*f[2]*volFact*wx2+0.5773502691896258*f[13]*dv2*volFact; 
  out[3] += 2.0*f[3]*volFact*wx2+0.5773502691896258*f[14]*dv2*volFact; 
  out[4] += 2.0*f[6]*volFact*wx2+0.5773502691896258*f[25]*dv2*volFact; 
  out[5] += 2.0*f[7]*volFact*wx2+0.5773502691896258*f[26]*dv2*volFact; 
  out[6] += 2.0*f[8]*volFact*wx2+0.5773502691896258*f[27]*dv2*volFact; 
  out[7] += 2.0*f[16]*volFact*wx2+0.5773502691896257*f[43]*dv2*volFact; 
  out[8] += 2.0*f[17]*volFact*wx2+0.5773502691896257*f[44]*dv2*volFact; 
  out[9] += 2.0*f[18]*volFact*wx2+0.5773502691896257*f[45]*dv2*volFact; 
  out[10] += 2.0*f[21]*volFact*wx2+0.5773502691896258*f[57]*dv2*volFact; 
  out[11] += 2.0*f[31]*volFact*wx2+0.5773502691896257*f[73]*dv2*volFact; 
  out[12] += 2.0*f[32]*volFact*wx2+0.5773502691896257*f[74]*dv2*volFact; 
  out[13] += 2.0*f[33]*volFact*wx2+0.5773502691896257*f[75]*dv2*volFact; 
  out[14] += 2.0*f[34]*volFact*wx2+0.5773502691896257*f[76]*dv2*volFact; 
  out[15] += 2.0*f[35]*volFact*wx2+0.5773502691896257*f[77]*dv2*volFact; 
  out[16] += 2.0*f[36]*volFact*wx2+0.5773502691896257*f[78]*dv2*volFact; 
  out[17] += 2.0*f[51]*volFact*wx2+0.5773502691896256*f[113]*dv2*volFact; 
  out[18] += 2.0*f[52]*volFact*wx2+0.5773502691896256*f[114]*dv2*volFact; 
  out[19] += 2.0*f[53]*volFact*wx2+0.5773502691896256*f[115]*dv2*volFact; 
  out[20] += 2.0*f[61]*volFact*wx2; 
  out[21] += 2.0*f[62]*volFact*wx2; 
  out[22] += 2.0*f[63]*volFact*wx2; 
  out[23] += 2.0*f[91]*volFact*wx2; 
  out[24] += 2.0*f[92]*volFact*wx2; 
  out[25] += 2.0*f[93]*volFact*wx2; 
  out[26] += 2.0*f[101]*volFact*wx2; 
  out[27] += 2.0*f[102]*volFact*wx2; 
  out[28] += 2.0*f[103]*volFact*wx2; 
  out[29] += 2.0*f[104]*volFact*wx2; 
  out[30] += 2.0*f[105]*volFact*wx2; 
  out[31] += 2.0*f[106]*volFact*wx2; 
  out[32] += 2.0*f[121]*volFact*wx2; 
  out[33] += 2.0*f[122]*volFact*wx2; 
  out[34] += 2.0*f[123]*volFact*wx2; 
} 
void GkMomentCalc3x2vMax_Qpar_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1*wx1_sq+1.732050807568877*f[4]*dv1*volFact*wx1_sq+0.5*f[0]*dv1_sq*volFact*wx1+0.08660254037844387*f[4]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx1_sq+0.5*f[1]*dv1_sq*volFact*wx1; 
  out[2] += 2.0*f[2]*volFact*wx1*wx1_sq+0.5*f[2]*dv1_sq*volFact*wx1; 
  out[3] += 2.0*f[3]*volFact*wx1*wx1_sq+0.5*f[3]*dv1_sq*volFact*wx1; 
} 
void GkMomentCalc3x2vMax_Qpar_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1*wx1_sq+1.732050807568877*f[4]*dv1*volFact*wx1_sq+0.4472135954999579*f[19]*dv1_sq*volFact*wx1+0.5*f[0]*dv1_sq*volFact*wx1+0.08660254037844387*f[4]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx1_sq+1.732050807568877*f[9]*dv1*volFact*wx1_sq+0.5*f[1]*dv1_sq*volFact*wx1+0.08660254037844387*f[9]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1*wx1_sq+1.732050807568877*f[10]*dv1*volFact*wx1_sq+0.5*f[2]*dv1_sq*volFact*wx1+0.08660254037844387*f[10]*dv1*dv1_sq*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1*wx1_sq+1.732050807568877*f[11]*dv1*volFact*wx1_sq+0.5*f[3]*dv1_sq*volFact*wx1+0.08660254037844387*f[11]*dv1*dv1_sq*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1*wx1_sq+0.5*f[6]*dv1_sq*volFact*wx1; 
  out[5] += 2.0*f[7]*volFact*wx1*wx1_sq+0.5*f[7]*dv1_sq*volFact*wx1; 
  out[6] += 2.0*f[8]*volFact*wx1*wx1_sq+0.5*f[8]*dv1_sq*volFact*wx1; 
  out[7] += 2.0*f[16]*volFact*wx1*wx1_sq+0.5*f[16]*dv1_sq*volFact*wx1; 
  out[8] += 2.0*f[17]*volFact*wx1*wx1_sq+0.5*f[17]*dv1_sq*volFact*wx1; 
  out[9] += 2.0*f[18]*volFact*wx1*wx1_sq+0.5*f[18]*dv1_sq*volFact*wx1; 
} 
void GkMomentCalc3x2vMax_Qpar_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1*wx1_sq+1.732050807568877*f[4]*dv1*volFact*wx1_sq+0.4472135954999579*f[19]*dv1_sq*volFact*wx1+0.5*f[0]*dv1_sq*volFact*wx1+0.03779644730092272*f[54]*dv1*dv1_sq*volFact+0.08660254037844387*f[4]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx1_sq+1.732050807568877*f[9]*dv1*volFact*wx1_sq+0.447213595499958*f[40]*dv1_sq*volFact*wx1+0.5*f[1]*dv1_sq*volFact*wx1+0.08660254037844387*f[9]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1*wx1_sq+1.732050807568877*f[10]*dv1*volFact*wx1_sq+0.447213595499958*f[41]*dv1_sq*volFact*wx1+0.5*f[2]*dv1_sq*volFact*wx1+0.08660254037844387*f[10]*dv1*dv1_sq*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1*wx1_sq+1.732050807568877*f[11]*dv1*volFact*wx1_sq+0.447213595499958*f[42]*dv1_sq*volFact*wx1+0.5*f[3]*dv1_sq*volFact*wx1+0.08660254037844387*f[11]*dv1*dv1_sq*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1*wx1_sq+1.732050807568877*f[22]*dv1*volFact*wx1_sq+0.5*f[6]*dv1_sq*volFact*wx1+0.08660254037844387*f[22]*dv1*dv1_sq*volFact; 
  out[5] += 2.0*f[7]*volFact*wx1*wx1_sq+1.732050807568877*f[23]*dv1*volFact*wx1_sq+0.5*f[7]*dv1_sq*volFact*wx1+0.08660254037844387*f[23]*dv1*dv1_sq*volFact; 
  out[6] += 2.0*f[8]*volFact*wx1*wx1_sq+1.732050807568877*f[24]*dv1*volFact*wx1_sq+0.5*f[8]*dv1_sq*volFact*wx1+0.08660254037844387*f[24]*dv1*dv1_sq*volFact; 
  out[7] += 2.0*f[16]*volFact*wx1*wx1_sq+1.732050807568877*f[37]*dv1*volFact*wx1_sq+0.5*f[16]*dv1_sq*volFact*wx1+0.08660254037844385*f[37]*dv1*dv1_sq*volFact; 
  out[8] += 2.0*f[17]*volFact*wx1*wx1_sq+1.732050807568877*f[38]*dv1*volFact*wx1_sq+0.5*f[17]*dv1_sq*volFact*wx1+0.08660254037844385*f[38]*dv1*dv1_sq*volFact; 
  out[9] += 2.0*f[18]*volFact*wx1*wx1_sq+1.732050807568877*f[39]*dv1*volFact*wx1_sq+0.5*f[18]*dv1_sq*volFact*wx1+0.08660254037844385*f[39]*dv1*dv1_sq*volFact; 
  out[10] += 2.0*f[21]*volFact*wx1*wx1_sq+0.5*f[21]*dv1_sq*volFact*wx1; 
  out[11] += 2.0*f[31]*volFact*wx1*wx1_sq+0.5*f[31]*dv1_sq*volFact*wx1; 
  out[12] += 2.0*f[32]*volFact*wx1*wx1_sq+0.5*f[32]*dv1_sq*volFact*wx1; 
  out[13] += 2.0*f[33]*volFact*wx1*wx1_sq+0.5*f[33]*dv1_sq*volFact*wx1; 
  out[14] += 2.0*f[34]*volFact*wx1*wx1_sq+0.5*f[34]*dv1_sq*volFact*wx1; 
  out[15] += 2.0*f[35]*volFact*wx1*wx1_sq+0.5*f[35]*dv1_sq*volFact*wx1; 
  out[16] += 2.0*f[36]*volFact*wx1*wx1_sq+0.5*f[36]*dv1_sq*volFact*wx1; 
  out[17] += 2.0*f[51]*volFact*wx1*wx1_sq+0.5*f[51]*dv1_sq*volFact*wx1; 
  out[18] += 2.0*f[52]*volFact*wx1*wx1_sq+0.5*f[52]*dv1_sq*volFact*wx1; 
  out[19] += 2.0*f[53]*volFact*wx1*wx1_sq+0.5*f[53]*dv1_sq*volFact*wx1; 
} 
void GkMomentCalc3x2vMax_Qpar_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1*wx1_sq+1.732050807568877*f[4]*dv1*volFact*wx1_sq+0.4472135954999579*f[19]*dv1_sq*volFact*wx1+0.5*f[0]*dv1_sq*volFact*wx1+0.03779644730092272*f[54]*dv1*dv1_sq*volFact+0.08660254037844387*f[4]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx1_sq+1.732050807568877*f[9]*dv1*volFact*wx1_sq+0.447213595499958*f[40]*dv1_sq*volFact*wx1+0.5*f[1]*dv1_sq*volFact*wx1+0.03779644730092273*f[110]*dv1*dv1_sq*volFact+0.08660254037844387*f[9]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1*wx1_sq+1.732050807568877*f[10]*dv1*volFact*wx1_sq+0.447213595499958*f[41]*dv1_sq*volFact*wx1+0.5*f[2]*dv1_sq*volFact*wx1+0.03779644730092273*f[111]*dv1*dv1_sq*volFact+0.08660254037844387*f[10]*dv1*dv1_sq*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1*wx1_sq+1.732050807568877*f[11]*dv1*volFact*wx1_sq+0.447213595499958*f[42]*dv1_sq*volFact*wx1+0.5*f[3]*dv1_sq*volFact*wx1+0.03779644730092273*f[112]*dv1*dv1_sq*volFact+0.08660254037844387*f[11]*dv1*dv1_sq*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1*wx1_sq+1.732050807568877*f[22]*dv1*volFact*wx1_sq+0.4472135954999579*f[70]*dv1_sq*volFact*wx1+0.5*f[6]*dv1_sq*volFact*wx1+0.08660254037844387*f[22]*dv1*dv1_sq*volFact; 
  out[5] += 2.0*f[7]*volFact*wx1*wx1_sq+1.732050807568877*f[23]*dv1*volFact*wx1_sq+0.4472135954999579*f[71]*dv1_sq*volFact*wx1+0.5*f[7]*dv1_sq*volFact*wx1+0.08660254037844387*f[23]*dv1*dv1_sq*volFact; 
  out[6] += 2.0*f[8]*volFact*wx1*wx1_sq+1.732050807568877*f[24]*dv1*volFact*wx1_sq+0.4472135954999579*f[72]*dv1_sq*volFact*wx1+0.5*f[8]*dv1_sq*volFact*wx1+0.08660254037844387*f[24]*dv1*dv1_sq*volFact; 
  out[7] += 2.0*f[16]*volFact*wx1*wx1_sq+1.732050807568877*f[37]*dv1*volFact*wx1_sq+0.4472135954999579*f[94]*dv1_sq*volFact*wx1+0.5*f[16]*dv1_sq*volFact*wx1+0.08660254037844385*f[37]*dv1*dv1_sq*volFact; 
  out[8] += 2.0*f[17]*volFact*wx1*wx1_sq+1.732050807568877*f[38]*dv1*volFact*wx1_sq+0.4472135954999579*f[95]*dv1_sq*volFact*wx1+0.5*f[17]*dv1_sq*volFact*wx1+0.08660254037844385*f[38]*dv1*dv1_sq*volFact; 
  out[9] += 2.0*f[18]*volFact*wx1*wx1_sq+1.732050807568877*f[39]*dv1*volFact*wx1_sq+0.4472135954999579*f[96]*dv1_sq*volFact*wx1+0.5*f[18]*dv1_sq*volFact*wx1+0.08660254037844385*f[39]*dv1*dv1_sq*volFact; 
  out[10] += 2.0*f[21]*volFact*wx1*wx1_sq+1.732050807568877*f[56]*dv1*volFact*wx1_sq+0.5*f[21]*dv1_sq*volFact*wx1+0.08660254037844387*f[56]*dv1*dv1_sq*volFact; 
  out[11] += 2.0*f[31]*volFact*wx1*wx1_sq+1.732050807568877*f[64]*dv1*volFact*wx1_sq+0.5*f[31]*dv1_sq*volFact*wx1+0.08660254037844385*f[64]*dv1*dv1_sq*volFact; 
  out[12] += 2.0*f[32]*volFact*wx1*wx1_sq+1.732050807568877*f[65]*dv1*volFact*wx1_sq+0.5*f[32]*dv1_sq*volFact*wx1+0.08660254037844385*f[65]*dv1*dv1_sq*volFact; 
  out[13] += 2.0*f[33]*volFact*wx1*wx1_sq+1.732050807568877*f[66]*dv1*volFact*wx1_sq+0.5*f[33]*dv1_sq*volFact*wx1+0.08660254037844385*f[66]*dv1*dv1_sq*volFact; 
  out[14] += 2.0*f[34]*volFact*wx1*wx1_sq+1.732050807568877*f[67]*dv1*volFact*wx1_sq+0.5*f[34]*dv1_sq*volFact*wx1+0.08660254037844385*f[67]*dv1*dv1_sq*volFact; 
  out[15] += 2.0*f[35]*volFact*wx1*wx1_sq+1.732050807568877*f[68]*dv1*volFact*wx1_sq+0.5*f[35]*dv1_sq*volFact*wx1+0.08660254037844385*f[68]*dv1*dv1_sq*volFact; 
  out[16] += 2.0*f[36]*volFact*wx1*wx1_sq+1.732050807568877*f[69]*dv1*volFact*wx1_sq+0.5*f[36]*dv1_sq*volFact*wx1+0.08660254037844385*f[69]*dv1*dv1_sq*volFact; 
  out[17] += 2.0*f[51]*volFact*wx1*wx1_sq+1.732050807568877*f[107]*dv1*volFact*wx1_sq+0.5*f[51]*dv1_sq*volFact*wx1+0.08660254037844385*f[107]*dv1*dv1_sq*volFact; 
  out[18] += 2.0*f[52]*volFact*wx1*wx1_sq+1.732050807568877*f[108]*dv1*volFact*wx1_sq+0.5*f[52]*dv1_sq*volFact*wx1+0.08660254037844385*f[108]*dv1*dv1_sq*volFact; 
  out[19] += 2.0*f[53]*volFact*wx1*wx1_sq+1.732050807568877*f[109]*dv1*volFact*wx1_sq+0.5*f[53]*dv1_sq*volFact*wx1+0.08660254037844385*f[109]*dv1*dv1_sq*volFact; 
  out[20] += 2.0*f[61]*volFact*wx1*wx1_sq+0.5*f[61]*dv1_sq*volFact*wx1; 
  out[21] += 2.0*f[62]*volFact*wx1*wx1_sq+0.5*f[62]*dv1_sq*volFact*wx1; 
  out[22] += 2.0*f[63]*volFact*wx1*wx1_sq+0.5*f[63]*dv1_sq*volFact*wx1; 
  out[23] += 2.0*f[91]*volFact*wx1*wx1_sq+0.5*f[91]*dv1_sq*volFact*wx1; 
  out[24] += 2.0*f[92]*volFact*wx1*wx1_sq+0.5*f[92]*dv1_sq*volFact*wx1; 
  out[25] += 2.0*f[93]*volFact*wx1*wx1_sq+0.5*f[93]*dv1_sq*volFact*wx1; 
  out[26] += 2.0*f[101]*volFact*wx1*wx1_sq+0.5*f[101]*dv1_sq*volFact*wx1; 
  out[27] += 2.0*f[102]*volFact*wx1*wx1_sq+0.5*f[102]*dv1_sq*volFact*wx1; 
  out[28] += 2.0*f[103]*volFact*wx1*wx1_sq+0.5*f[103]*dv1_sq*volFact*wx1; 
  out[29] += 2.0*f[104]*volFact*wx1*wx1_sq+0.5*f[104]*dv1_sq*volFact*wx1; 
  out[30] += 2.0*f[105]*volFact*wx1*wx1_sq+0.5*f[105]*dv1_sq*volFact*wx1; 
  out[31] += 2.0*f[106]*volFact*wx1*wx1_sq+0.5*f[106]*dv1_sq*volFact*wx1; 
  out[32] += 2.0*f[121]*volFact*wx1*wx1_sq+0.5*f[121]*dv1_sq*volFact*wx1; 
  out[33] += 2.0*f[122]*volFact*wx1*wx1_sq+0.5*f[122]*dv1_sq*volFact*wx1; 
  out[34] += 2.0*f[123]*volFact*wx1*wx1_sq+0.5*f[123]*dv1_sq*volFact*wx1; 
} 
void GkMomentCalc3x2vMax_Qperp_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx1*wx2+0.5773502691896258*f[4]*dv1*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact*wx1; 
  out[1] += 2.0*f[1]*volFact*wx1*wx2; 
  out[2] += 2.0*f[2]*volFact*wx1*wx2; 
  out[3] += 2.0*f[3]*volFact*wx1*wx2; 
} 
void GkMomentCalc3x2vMax_Qperp_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx1*wx2+0.5773502691896258*f[4]*dv1*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact*wx1+0.1666666666666667*f[15]*dv1*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx2+0.5773502691896258*f[9]*dv1*volFact*wx2+0.5773502691896258*f[12]*dv2*volFact*wx1; 
  out[2] += 2.0*f[2]*volFact*wx1*wx2+0.5773502691896258*f[10]*dv1*volFact*wx2+0.5773502691896258*f[13]*dv2*volFact*wx1; 
  out[3] += 2.0*f[3]*volFact*wx1*wx2+0.5773502691896258*f[11]*dv1*volFact*wx2+0.5773502691896258*f[14]*dv2*volFact*wx1; 
  out[4] += 2.0*f[6]*volFact*wx1*wx2; 
  out[5] += 2.0*f[7]*volFact*wx1*wx2; 
  out[6] += 2.0*f[8]*volFact*wx1*wx2; 
  out[7] += 2.0*f[16]*volFact*wx1*wx2; 
  out[8] += 2.0*f[17]*volFact*wx1*wx2; 
  out[9] += 2.0*f[18]*volFact*wx1*wx2; 
} 
void GkMomentCalc3x2vMax_Qperp_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx1*wx2+0.5773502691896258*f[4]*dv1*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact*wx1+0.1666666666666667*f[15]*dv1*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx2+0.5773502691896258*f[9]*dv1*volFact*wx2+0.5773502691896258*f[12]*dv2*volFact*wx1+0.1666666666666667*f[28]*dv1*dv2*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1*wx2+0.5773502691896258*f[10]*dv1*volFact*wx2+0.5773502691896258*f[13]*dv2*volFact*wx1+0.1666666666666667*f[29]*dv1*dv2*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1*wx2+0.5773502691896258*f[11]*dv1*volFact*wx2+0.5773502691896258*f[14]*dv2*volFact*wx1+0.1666666666666667*f[30]*dv1*dv2*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1*wx2+0.5773502691896258*f[22]*dv1*volFact*wx2+0.5773502691896258*f[25]*dv2*volFact*wx1; 
  out[5] += 2.0*f[7]*volFact*wx1*wx2+0.5773502691896258*f[23]*dv1*volFact*wx2+0.5773502691896258*f[26]*dv2*volFact*wx1; 
  out[6] += 2.0*f[8]*volFact*wx1*wx2+0.5773502691896258*f[24]*dv1*volFact*wx2+0.5773502691896258*f[27]*dv2*volFact*wx1; 
  out[7] += 2.0*f[16]*volFact*wx1*wx2+0.5773502691896257*f[37]*dv1*volFact*wx2+0.5773502691896257*f[43]*dv2*volFact*wx1; 
  out[8] += 2.0*f[17]*volFact*wx1*wx2+0.5773502691896257*f[38]*dv1*volFact*wx2+0.5773502691896257*f[44]*dv2*volFact*wx1; 
  out[9] += 2.0*f[18]*volFact*wx1*wx2+0.5773502691896257*f[39]*dv1*volFact*wx2+0.5773502691896257*f[45]*dv2*volFact*wx1; 
  out[10] += 2.0*f[21]*volFact*wx1*wx2; 
  out[11] += 2.0*f[31]*volFact*wx1*wx2; 
  out[12] += 2.0*f[32]*volFact*wx1*wx2; 
  out[13] += 2.0*f[33]*volFact*wx1*wx2; 
  out[14] += 2.0*f[34]*volFact*wx1*wx2; 
  out[15] += 2.0*f[35]*volFact*wx1*wx2; 
  out[16] += 2.0*f[36]*volFact*wx1*wx2; 
  out[17] += 2.0*f[51]*volFact*wx1*wx2; 
  out[18] += 2.0*f[52]*volFact*wx1*wx2; 
  out[19] += 2.0*f[53]*volFact*wx1*wx2; 
} 
void GkMomentCalc3x2vMax_Qperp_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += 2.0*f[0]*volFact*wx1*wx2+0.5773502691896258*f[4]*dv1*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact*wx1+0.1666666666666667*f[15]*dv1*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx2+0.5773502691896258*f[9]*dv1*volFact*wx2+0.5773502691896258*f[12]*dv2*volFact*wx1+0.1666666666666667*f[28]*dv1*dv2*volFact; 
  out[2] += 2.0*f[2]*volFact*wx1*wx2+0.5773502691896258*f[10]*dv1*volFact*wx2+0.5773502691896258*f[13]*dv2*volFact*wx1+0.1666666666666667*f[29]*dv1*dv2*volFact; 
  out[3] += 2.0*f[3]*volFact*wx1*wx2+0.5773502691896258*f[11]*dv1*volFact*wx2+0.5773502691896258*f[14]*dv2*volFact*wx1+0.1666666666666667*f[30]*dv1*dv2*volFact; 
  out[4] += 2.0*f[6]*volFact*wx1*wx2+0.5773502691896258*f[22]*dv1*volFact*wx2+0.5773502691896258*f[25]*dv2*volFact*wx1+0.1666666666666667*f[58]*dv1*dv2*volFact; 
  out[5] += 2.0*f[7]*volFact*wx1*wx2+0.5773502691896258*f[23]*dv1*volFact*wx2+0.5773502691896258*f[26]*dv2*volFact*wx1+0.1666666666666667*f[59]*dv1*dv2*volFact; 
  out[6] += 2.0*f[8]*volFact*wx1*wx2+0.5773502691896258*f[24]*dv1*volFact*wx2+0.5773502691896258*f[27]*dv2*volFact*wx1+0.1666666666666667*f[60]*dv1*dv2*volFact; 
  out[7] += 2.0*f[16]*volFact*wx1*wx2+0.5773502691896257*f[37]*dv1*volFact*wx2+0.5773502691896257*f[43]*dv2*volFact*wx1+0.1666666666666667*f[79]*dv1*dv2*volFact; 
  out[8] += 2.0*f[17]*volFact*wx1*wx2+0.5773502691896257*f[38]*dv1*volFact*wx2+0.5773502691896257*f[44]*dv2*volFact*wx1+0.1666666666666667*f[80]*dv1*dv2*volFact; 
  out[9] += 2.0*f[18]*volFact*wx1*wx2+0.5773502691896257*f[39]*dv1*volFact*wx2+0.5773502691896257*f[45]*dv2*volFact*wx1+0.1666666666666667*f[81]*dv1*dv2*volFact; 
  out[10] += 2.0*f[21]*volFact*wx1*wx2+0.5773502691896258*f[56]*dv1*volFact*wx2+0.5773502691896258*f[57]*dv2*volFact*wx1; 
  out[11] += 2.0*f[31]*volFact*wx1*wx2+0.5773502691896257*f[64]*dv1*volFact*wx2+0.5773502691896257*f[73]*dv2*volFact*wx1; 
  out[12] += 2.0*f[32]*volFact*wx1*wx2+0.5773502691896257*f[65]*dv1*volFact*wx2+0.5773502691896257*f[74]*dv2*volFact*wx1; 
  out[13] += 2.0*f[33]*volFact*wx1*wx2+0.5773502691896257*f[66]*dv1*volFact*wx2+0.5773502691896257*f[75]*dv2*volFact*wx1; 
  out[14] += 2.0*f[34]*volFact*wx1*wx2+0.5773502691896257*f[67]*dv1*volFact*wx2+0.5773502691896257*f[76]*dv2*volFact*wx1; 
  out[15] += 2.0*f[35]*volFact*wx1*wx2+0.5773502691896257*f[68]*dv1*volFact*wx2+0.5773502691896257*f[77]*dv2*volFact*wx1; 
  out[16] += 2.0*f[36]*volFact*wx1*wx2+0.5773502691896257*f[69]*dv1*volFact*wx2+0.5773502691896257*f[78]*dv2*volFact*wx1; 
  out[17] += 2.0*f[51]*volFact*wx1*wx2+0.5773502691896256*f[107]*dv1*volFact*wx2+0.5773502691896256*f[113]*dv2*volFact*wx1; 
  out[18] += 2.0*f[52]*volFact*wx1*wx2+0.5773502691896256*f[108]*dv1*volFact*wx2+0.5773502691896256*f[114]*dv2*volFact*wx1; 
  out[19] += 2.0*f[53]*volFact*wx1*wx2+0.5773502691896256*f[109]*dv1*volFact*wx2+0.5773502691896256*f[115]*dv2*volFact*wx1; 
  out[20] += 2.0*f[61]*volFact*wx1*wx2; 
  out[21] += 2.0*f[62]*volFact*wx1*wx2; 
  out[22] += 2.0*f[63]*volFact*wx1*wx2; 
  out[23] += 2.0*f[91]*volFact*wx1*wx2; 
  out[24] += 2.0*f[92]*volFact*wx1*wx2; 
  out[25] += 2.0*f[93]*volFact*wx1*wx2; 
  out[26] += 2.0*f[101]*volFact*wx1*wx2; 
  out[27] += 2.0*f[102]*volFact*wx1*wx2; 
  out[28] += 2.0*f[103]*volFact*wx1*wx2; 
  out[29] += 2.0*f[104]*volFact*wx1*wx2; 
  out[30] += 2.0*f[105]*volFact*wx1*wx2; 
  out[31] += 2.0*f[106]*volFact*wx1*wx2; 
  out[32] += 2.0*f[121]*volFact*wx1*wx2; 
  out[33] += 2.0*f[122]*volFact*wx1*wx2; 
  out[34] += 2.0*f[123]*volFact*wx1*wx2; 
} 
