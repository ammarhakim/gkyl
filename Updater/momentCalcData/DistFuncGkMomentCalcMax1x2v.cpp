#include <math.h> 
#include <DistFuncMomentCalcModDecl.h> 
void GkMomentCalc1x2vMax_Dens_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
void GkMomentCalc1x2vMax_Dens_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
} 
void GkMomentCalc1x2vMax_Dens_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
  out[3] += 2.0*f[17]*volFact; 
} 
void GkMomentCalc1x2vMax_Dens_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
  out[3] += 2.0*f[17]*volFact; 
  out[4] += 2.0*f[32]*volFact; 
} 
void GkMomentCalc1x2vMax_Upar_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx1+0.5773502691896258*f[2]*dv1*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1; 
} 
void GkMomentCalc1x2vMax_Upar_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx1+0.5773502691896258*f[2]*dv1*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1+0.5773502691896258*f[4]*dv1*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1; 
} 
void GkMomentCalc1x2vMax_Upar_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx1+0.5773502691896258*f[2]*dv1*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1+0.5773502691896258*f[4]*dv1*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1+0.5773502691896257*f[11]*dv1*volFact; 
  out[3] += 2.0*f[17]*volFact*wx1; 
} 
void GkMomentCalc1x2vMax_Upar_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx1+0.5773502691896258*f[2]*dv1*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1+0.5773502691896258*f[4]*dv1*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1+0.5773502691896257*f[11]*dv1*volFact; 
  out[3] += 2.0*f[17]*volFact*wx1+0.5773502691896256*f[26]*dv1*volFact; 
  out[4] += 2.0*f[32]*volFact*wx1; 
} 
void GkMomentCalc1x2vMax_Ppar_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[2]*dv1*volFact*wx1+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1_sq+0.1666666666666667*f[1]*dv1_sq*volFact; 
} 
void GkMomentCalc1x2vMax_Ppar_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[2]*dv1*volFact*wx1+0.149071198499986*f[8]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1_sq+1.154700538379252*f[4]*dv1*volFact*wx1+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1_sq+0.1666666666666667*f[7]*dv1_sq*volFact; 
} 
void GkMomentCalc1x2vMax_Ppar_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[2]*dv1*volFact*wx1+0.149071198499986*f[8]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1_sq+1.154700538379252*f[4]*dv1*volFact*wx1+0.149071198499986*f[12]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1_sq+1.154700538379251*f[11]*dv1*volFact*wx1+0.1666666666666667*f[7]*dv1_sq*volFact; 
  out[3] += 2.0*f[17]*volFact*wx1_sq+0.1666666666666667*f[17]*dv1_sq*volFact; 
} 
void GkMomentCalc1x2vMax_Ppar_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1_sq+1.154700538379252*f[2]*dv1*volFact*wx1+0.149071198499986*f[8]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1_sq+1.154700538379252*f[4]*dv1*volFact*wx1+0.149071198499986*f[12]*dv1_sq*volFact+0.1666666666666667*f[1]*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1_sq+1.154700538379251*f[11]*dv1*volFact*wx1+0.149071198499986*f[23]*dv1_sq*volFact+0.1666666666666667*f[7]*dv1_sq*volFact; 
  out[3] += 2.0*f[17]*volFact*wx1_sq+1.154700538379251*f[26]*dv1*volFact*wx1+0.1666666666666667*f[17]*dv1_sq*volFact; 
  out[4] += 2.0*f[32]*volFact*wx1_sq+0.1666666666666667*f[32]*dv1_sq*volFact; 
} 
void GkMomentCalc1x2vMax_Pperp_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx2+0.5773502691896258*f[3]*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx2; 
} 
void GkMomentCalc1x2vMax_Pperp_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx2+0.5773502691896258*f[3]*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact; 
  out[2] += 2.0*f[7]*volFact*wx2; 
} 
void GkMomentCalc1x2vMax_Pperp_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx2+0.5773502691896258*f[3]*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact; 
  out[2] += 2.0*f[7]*volFact*wx2+0.5773502691896257*f[13]*dv2*volFact; 
  out[3] += 2.0*f[17]*volFact*wx2; 
} 
void GkMomentCalc1x2vMax_Pperp_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx2+0.5773502691896258*f[3]*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact; 
  out[2] += 2.0*f[7]*volFact*wx2+0.5773502691896257*f[13]*dv2*volFact; 
  out[3] += 2.0*f[17]*volFact*wx2+0.5773502691896256*f[28]*dv2*volFact; 
  out[4] += 2.0*f[32]*volFact*wx2; 
} 
void GkMomentCalc1x2vMax_Qpar_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1*wx1_sq+1.732050807568877*f[2]*dv1*volFact*wx1_sq+0.5*f[0]*dv1_sq*volFact*wx1+0.08660254037844387*f[2]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx1_sq+0.5*f[1]*dv1_sq*volFact*wx1; 
} 
void GkMomentCalc1x2vMax_Qpar_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1*wx1_sq+1.732050807568877*f[2]*dv1*volFact*wx1_sq+0.4472135954999579*f[8]*dv1_sq*volFact*wx1+0.5*f[0]*dv1_sq*volFact*wx1+0.08660254037844387*f[2]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx1_sq+1.732050807568877*f[4]*dv1*volFact*wx1_sq+0.5*f[1]*dv1_sq*volFact*wx1+0.08660254037844387*f[4]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1*wx1_sq+0.5*f[7]*dv1_sq*volFact*wx1; 
} 
void GkMomentCalc1x2vMax_Qpar_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1*wx1_sq+1.732050807568877*f[2]*dv1*volFact*wx1_sq+0.4472135954999579*f[8]*dv1_sq*volFact*wx1+0.5*f[0]*dv1_sq*volFact*wx1+0.03779644730092272*f[18]*dv1*dv1_sq*volFact+0.08660254037844387*f[2]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx1_sq+1.732050807568877*f[4]*dv1*volFact*wx1_sq+0.447213595499958*f[12]*dv1_sq*volFact*wx1+0.5*f[1]*dv1_sq*volFact*wx1+0.08660254037844387*f[4]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1*wx1_sq+1.732050807568877*f[11]*dv1*volFact*wx1_sq+0.5*f[7]*dv1_sq*volFact*wx1+0.08660254037844385*f[11]*dv1*dv1_sq*volFact; 
  out[3] += 2.0*f[17]*volFact*wx1*wx1_sq+0.5*f[17]*dv1_sq*volFact*wx1; 
} 
void GkMomentCalc1x2vMax_Qpar_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += 2.0*f[0]*volFact*wx1*wx1_sq+1.732050807568877*f[2]*dv1*volFact*wx1_sq+0.4472135954999579*f[8]*dv1_sq*volFact*wx1+0.5*f[0]*dv1_sq*volFact*wx1+0.03779644730092272*f[18]*dv1*dv1_sq*volFact+0.08660254037844387*f[2]*dv1*dv1_sq*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx1_sq+1.732050807568877*f[4]*dv1*volFact*wx1_sq+0.447213595499958*f[12]*dv1_sq*volFact*wx1+0.5*f[1]*dv1_sq*volFact*wx1+0.03779644730092273*f[27]*dv1*dv1_sq*volFact+0.08660254037844387*f[4]*dv1*dv1_sq*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1*wx1_sq+1.732050807568877*f[11]*dv1*volFact*wx1_sq+0.4472135954999579*f[23]*dv1_sq*volFact*wx1+0.5*f[7]*dv1_sq*volFact*wx1+0.08660254037844385*f[11]*dv1*dv1_sq*volFact; 
  out[3] += 2.0*f[17]*volFact*wx1*wx1_sq+1.732050807568877*f[26]*dv1*volFact*wx1_sq+0.5*f[17]*dv1_sq*volFact*wx1+0.08660254037844385*f[26]*dv1*dv1_sq*volFact; 
  out[4] += 2.0*f[32]*volFact*wx1*wx1_sq+0.5*f[32]*dv1_sq*volFact*wx1; 
} 
void GkMomentCalc1x2vMax_Qperp_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx1*wx2+0.5773502691896258*f[2]*dv1*volFact*wx2+0.5773502691896258*f[3]*dv2*volFact*wx1; 
  out[1] += 2.0*f[1]*volFact*wx1*wx2; 
} 
void GkMomentCalc1x2vMax_Qperp_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx1*wx2+0.5773502691896258*f[2]*dv1*volFact*wx2+0.5773502691896258*f[3]*dv2*volFact*wx1+0.1666666666666667*f[6]*dv1*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx2+0.5773502691896258*f[4]*dv1*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact*wx1; 
  out[2] += 2.0*f[7]*volFact*wx1*wx2; 
} 
void GkMomentCalc1x2vMax_Qperp_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx1*wx2+0.5773502691896258*f[2]*dv1*volFact*wx2+0.5773502691896258*f[3]*dv2*volFact*wx1+0.1666666666666667*f[6]*dv1*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx2+0.5773502691896258*f[4]*dv1*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact*wx1+0.1666666666666667*f[10]*dv1*dv2*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1*wx2+0.5773502691896257*f[11]*dv1*volFact*wx2+0.5773502691896257*f[13]*dv2*volFact*wx1; 
  out[3] += 2.0*f[17]*volFact*wx1*wx2; 
} 
void GkMomentCalc1x2vMax_Qperp_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += 2.0*f[0]*volFact*wx1*wx2+0.5773502691896258*f[2]*dv1*volFact*wx2+0.5773502691896258*f[3]*dv2*volFact*wx1+0.1666666666666667*f[6]*dv1*dv2*volFact; 
  out[1] += 2.0*f[1]*volFact*wx1*wx2+0.5773502691896258*f[4]*dv1*volFact*wx2+0.5773502691896258*f[5]*dv2*volFact*wx1+0.1666666666666667*f[10]*dv1*dv2*volFact; 
  out[2] += 2.0*f[7]*volFact*wx1*wx2+0.5773502691896257*f[11]*dv1*volFact*wx2+0.5773502691896257*f[13]*dv2*volFact*wx1+0.1666666666666667*f[20]*dv1*dv2*volFact; 
  out[3] += 2.0*f[17]*volFact*wx1*wx2+0.5773502691896256*f[26]*dv1*volFact*wx2+0.5773502691896256*f[28]*dv2*volFact*wx1; 
  out[4] += 2.0*f[32]*volFact*wx1*wx2; 
} 
