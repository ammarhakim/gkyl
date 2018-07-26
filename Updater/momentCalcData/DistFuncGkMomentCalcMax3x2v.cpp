#include <math.h> 
#include <DistFuncMomentCalcModDecl.h> 
void GkMomentCalc3x2vMax_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[3]*volFact; 
} 
void GkMomentCalc3x2vMax_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
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
void GkMomentCalc3x2vMax_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += (2.0*f[0]*wx1+0.5773502691896258*f[4]*dv1)*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact; 
  out[2] += 2.0*f[2]*wx1*volFact; 
  out[3] += 2.0*f[3]*wx1*volFact; 
} 
void GkMomentCalc3x2vMax_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += (2.0*f[0]*wx1+0.5773502691896258*f[4]*dv1)*volFact; 
  out[1] += (2.0*f[1]*wx1+0.5773502691896258*f[9]*dv1)*volFact; 
  out[2] += (2.0*f[2]*wx1+0.5773502691896258*f[10]*dv1)*volFact; 
  out[3] += (2.0*f[3]*wx1+0.5773502691896258*f[11]*dv1)*volFact; 
  out[4] += 2.0*f[6]*wx1*volFact; 
  out[5] += 2.0*f[7]*wx1*volFact; 
  out[6] += 2.0*f[8]*wx1*volFact; 
  out[7] += 2.0*f[16]*wx1*volFact; 
  out[8] += 2.0*f[17]*wx1*volFact; 
  out[9] += 2.0*f[18]*wx1*volFact; 
} 
void GkMomentCalc3x2vMax_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double zeta[6]; 

  zeta[0] = 5.656854249492382*wx1_sq+0.4714045207910317*dv1_sq; 
  zeta[4] = 3.265986323710906*dv1*wx1; 
  out[0] += (0.3535533905932737*f[4]*zeta[4]+0.3535533905932737*f[0]*zeta[0])*volFact; 
  out[1] += 0.3535533905932737*zeta[0]*f[1]*volFact; 
  out[2] += 0.3535533905932737*zeta[0]*f[2]*volFact; 
  out[3] += 0.3535533905932737*zeta[0]*f[3]*volFact; 
  double tmp[4]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[1] = 2.0*f[1]*wx2; 
  tmp[2] = 2.0*f[2]*wx2; 
  tmp[3] = 2.0*f[3]*wx2; 
  out[0] += (2.0*(0.3535533905932737*Bmag[3]*tmp[3]+0.3535533905932737*Bmag[2]*tmp[2]+0.3535533905932737*Bmag[1]*tmp[1]+0.3535533905932737*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.3535533905932737*Bmag[0]*tmp[1]+0.3535533905932737*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.3535533905932737*Bmag[0]*tmp[2]+0.3535533905932737*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.3535533905932737*Bmag[0]*tmp[3]+0.3535533905932737*tmp[0]*Bmag[3])*volFact)/m_; 
} 
void GkMomentCalc3x2vMax_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.149071198499986*f[19]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[9]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1_sq+1.154700538379252*f[10]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  out[3] += (2.0*f[3]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[3]*dv1_sq)*volFact; 
  out[4] += (2.0*f[6]*wx1_sq+0.1666666666666667*f[6]*dv1_sq)*volFact; 
  out[5] += (2.0*f[7]*wx1_sq+0.1666666666666667*f[7]*dv1_sq)*volFact; 
  out[6] += (2.0*f[8]*wx1_sq+0.1666666666666667*f[8]*dv1_sq)*volFact; 
  out[7] += (2.0*f[16]*wx1_sq+0.1666666666666667*f[16]*dv1_sq)*volFact; 
  out[8] += (2.0*f[17]*wx1_sq+0.1666666666666667*f[17]*dv1_sq)*volFact; 
  out[9] += (2.0*f[18]*wx1_sq+0.1666666666666667*f[18]*dv1_sq)*volFact; 
  double tmp[10]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[12]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[13]*dv2; 
  tmp[3] = 2.0*f[3]*wx2+0.5773502691896258*f[14]*dv2; 
  tmp[4] = 2.0*f[6]*wx2; 
  tmp[5] = 2.0*f[7]*wx2; 
  tmp[6] = 2.0*f[8]*wx2; 
  tmp[7] = 2.0*f[16]*wx2; 
  tmp[8] = 2.0*f[17]*wx2; 
  tmp[9] = 2.0*f[18]*wx2; 
  out[0] += (2.0*(0.3535533905932737*Bmag[9]*tmp[9]+0.3535533905932737*Bmag[8]*tmp[8]+0.3535533905932737*Bmag[7]*tmp[7]+0.3535533905932737*Bmag[6]*tmp[6]+0.3535533905932737*Bmag[5]*tmp[5]+0.3535533905932737*Bmag[4]*tmp[4]+0.3535533905932737*Bmag[3]*tmp[3]+0.3535533905932737*Bmag[2]*tmp[2]+0.3535533905932737*Bmag[1]*tmp[1]+0.3535533905932737*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.3162277660168379*Bmag[1]*tmp[7]+0.3162277660168379*tmp[1]*Bmag[7]+0.3535533905932737*Bmag[3]*tmp[5]+0.3535533905932737*tmp[3]*Bmag[5]+0.3535533905932737*Bmag[2]*tmp[4]+0.3535533905932737*tmp[2]*Bmag[4]+0.3535533905932737*Bmag[0]*tmp[1]+0.3535533905932737*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.3162277660168379*Bmag[2]*tmp[8]+0.3162277660168379*tmp[2]*Bmag[8]+0.3535533905932737*Bmag[3]*tmp[6]+0.3535533905932737*tmp[3]*Bmag[6]+0.3535533905932737*Bmag[1]*tmp[4]+0.3535533905932737*tmp[1]*Bmag[4]+0.3535533905932737*Bmag[0]*tmp[2]+0.3535533905932737*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.3162277660168379*Bmag[3]*tmp[9]+0.3162277660168379*tmp[3]*Bmag[9]+0.3535533905932737*Bmag[2]*tmp[6]+0.3535533905932737*tmp[2]*Bmag[6]+0.3535533905932737*Bmag[1]*tmp[5]+0.3535533905932737*tmp[1]*Bmag[5]+0.3535533905932737*Bmag[0]*tmp[3]+0.3535533905932737*tmp[0]*Bmag[3])*volFact)/m_; 
  out[4] += (2.0*(0.3162277660168379*Bmag[4]*tmp[8]+0.3162277660168379*tmp[4]*Bmag[8]+0.3162277660168379*Bmag[4]*tmp[7]+0.3162277660168379*tmp[4]*Bmag[7]+0.3535533905932737*Bmag[5]*tmp[6]+0.3535533905932737*tmp[5]*Bmag[6]+0.3535533905932737*Bmag[0]*tmp[4]+0.3535533905932737*tmp[0]*Bmag[4]+0.3535533905932737*Bmag[1]*tmp[2]+0.3535533905932737*tmp[1]*Bmag[2])*volFact)/m_; 
  out[5] += (2.0*(0.3162277660168379*Bmag[5]*tmp[9]+0.3162277660168379*tmp[5]*Bmag[9]+0.3162277660168379*Bmag[5]*tmp[7]+0.3162277660168379*tmp[5]*Bmag[7]+0.3535533905932737*Bmag[4]*tmp[6]+0.3535533905932737*tmp[4]*Bmag[6]+0.3535533905932737*Bmag[0]*tmp[5]+0.3535533905932737*tmp[0]*Bmag[5]+0.3535533905932737*Bmag[1]*tmp[3]+0.3535533905932737*tmp[1]*Bmag[3])*volFact)/m_; 
  out[6] += (2.0*(0.3162277660168379*Bmag[6]*tmp[9]+0.3162277660168379*tmp[6]*Bmag[9]+0.3162277660168379*Bmag[6]*tmp[8]+0.3162277660168379*tmp[6]*Bmag[8]+0.3535533905932737*Bmag[0]*tmp[6]+0.3535533905932737*tmp[0]*Bmag[6]+0.3535533905932737*Bmag[4]*tmp[5]+0.3535533905932737*tmp[4]*Bmag[5]+0.3535533905932737*Bmag[2]*tmp[3]+0.3535533905932737*tmp[2]*Bmag[3])*volFact)/m_; 
  out[7] += (2.0*(0.2258769757263128*Bmag[7]*tmp[7]+0.3535533905932737*Bmag[0]*tmp[7]+0.3535533905932737*tmp[0]*Bmag[7]+0.3162277660168379*Bmag[5]*tmp[5]+0.3162277660168379*Bmag[4]*tmp[4]+0.3162277660168379*Bmag[1]*tmp[1])*volFact)/m_; 
  out[8] += (2.0*(0.2258769757263128*Bmag[8]*tmp[8]+0.3535533905932737*Bmag[0]*tmp[8]+0.3535533905932737*tmp[0]*Bmag[8]+0.3162277660168379*Bmag[6]*tmp[6]+0.3162277660168379*Bmag[4]*tmp[4]+0.3162277660168379*Bmag[2]*tmp[2])*volFact)/m_; 
  out[9] += (2.0*(0.2258769757263128*Bmag[9]*tmp[9]+0.3535533905932737*Bmag[0]*tmp[9]+0.3535533905932737*tmp[0]*Bmag[9]+0.3162277660168379*Bmag[6]*tmp[6]+0.3162277660168379*Bmag[5]*tmp[5]+0.3162277660168379*Bmag[3]*tmp[3])*volFact)/m_; 
} 
void GkMomentCalc3x2vMax_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double zeta[6]; 

  zeta[0] = 5.656854249492382*wx1_sq+0.4714045207910317*dv1_sq; 
  zeta[4] = 3.265986323710906*dv1*wx1; 
  out[0] += (0.3535533905932737*f[4]*zeta[4]+0.3535533905932737*f[0]*zeta[0])*volFact; 
  out[1] += 0.3535533905932737*zeta[0]*f[1]*volFact; 
  out[2] += 0.3535533905932737*zeta[0]*f[2]*volFact; 
  out[3] += 0.3535533905932737*zeta[0]*f[3]*volFact; 
} 
void GkMomentCalc3x2vMax_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.149071198499986*f[19]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[9]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1_sq+1.154700538379252*f[10]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  out[3] += (2.0*f[3]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[3]*dv1_sq)*volFact; 
  out[4] += (2.0*f[6]*wx1_sq+0.1666666666666667*f[6]*dv1_sq)*volFact; 
  out[5] += (2.0*f[7]*wx1_sq+0.1666666666666667*f[7]*dv1_sq)*volFact; 
  out[6] += (2.0*f[8]*wx1_sq+0.1666666666666667*f[8]*dv1_sq)*volFact; 
  out[7] += (2.0*f[16]*wx1_sq+0.1666666666666667*f[16]*dv1_sq)*volFact; 
  out[8] += (2.0*f[17]*wx1_sq+0.1666666666666667*f[17]*dv1_sq)*volFact; 
  out[9] += (2.0*f[18]*wx1_sq+0.1666666666666667*f[18]*dv1_sq)*volFact; 
} 
void GkMomentCalc3x2vMax_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  double tmp[4]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[1] = 2.0*f[1]*wx2; 
  tmp[2] = 2.0*f[2]*wx2; 
  tmp[3] = 2.0*f[3]*wx2; 
  out[0] += ((0.3535533905932737*Bmag[3]*tmp[3]+0.3535533905932737*Bmag[2]*tmp[2]+0.3535533905932737*Bmag[1]*tmp[1]+0.3535533905932737*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.3535533905932737*Bmag[0]*tmp[1]+0.3535533905932737*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += ((0.3535533905932737*Bmag[0]*tmp[2]+0.3535533905932737*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += ((0.3535533905932737*Bmag[0]*tmp[3]+0.3535533905932737*tmp[0]*Bmag[3])*volFact)/m_; 
} 
void GkMomentCalc3x2vMax_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  double tmp[10]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[12]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[13]*dv2; 
  tmp[3] = 2.0*f[3]*wx2+0.5773502691896258*f[14]*dv2; 
  tmp[4] = 2.0*f[6]*wx2; 
  tmp[5] = 2.0*f[7]*wx2; 
  tmp[6] = 2.0*f[8]*wx2; 
  tmp[7] = 2.0*f[16]*wx2; 
  tmp[8] = 2.0*f[17]*wx2; 
  tmp[9] = 2.0*f[18]*wx2; 
  out[0] += ((0.3535533905932737*Bmag[9]*tmp[9]+0.3535533905932737*Bmag[8]*tmp[8]+0.3535533905932737*Bmag[7]*tmp[7]+0.3535533905932737*Bmag[6]*tmp[6]+0.3535533905932737*Bmag[5]*tmp[5]+0.3535533905932737*Bmag[4]*tmp[4]+0.3535533905932737*Bmag[3]*tmp[3]+0.3535533905932737*Bmag[2]*tmp[2]+0.3535533905932737*Bmag[1]*tmp[1]+0.3535533905932737*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.3162277660168379*Bmag[1]*tmp[7]+0.3162277660168379*tmp[1]*Bmag[7]+0.3535533905932737*Bmag[3]*tmp[5]+0.3535533905932737*tmp[3]*Bmag[5]+0.3535533905932737*Bmag[2]*tmp[4]+0.3535533905932737*tmp[2]*Bmag[4]+0.3535533905932737*Bmag[0]*tmp[1]+0.3535533905932737*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += ((0.3162277660168379*Bmag[2]*tmp[8]+0.3162277660168379*tmp[2]*Bmag[8]+0.3535533905932737*Bmag[3]*tmp[6]+0.3535533905932737*tmp[3]*Bmag[6]+0.3535533905932737*Bmag[1]*tmp[4]+0.3535533905932737*tmp[1]*Bmag[4]+0.3535533905932737*Bmag[0]*tmp[2]+0.3535533905932737*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += ((0.3162277660168379*Bmag[3]*tmp[9]+0.3162277660168379*tmp[3]*Bmag[9]+0.3535533905932737*Bmag[2]*tmp[6]+0.3535533905932737*tmp[2]*Bmag[6]+0.3535533905932737*Bmag[1]*tmp[5]+0.3535533905932737*tmp[1]*Bmag[5]+0.3535533905932737*Bmag[0]*tmp[3]+0.3535533905932737*tmp[0]*Bmag[3])*volFact)/m_; 
  out[4] += ((0.3162277660168379*Bmag[4]*tmp[8]+0.3162277660168379*tmp[4]*Bmag[8]+0.3162277660168379*Bmag[4]*tmp[7]+0.3162277660168379*tmp[4]*Bmag[7]+0.3535533905932737*Bmag[5]*tmp[6]+0.3535533905932737*tmp[5]*Bmag[6]+0.3535533905932737*Bmag[0]*tmp[4]+0.3535533905932737*tmp[0]*Bmag[4]+0.3535533905932737*Bmag[1]*tmp[2]+0.3535533905932737*tmp[1]*Bmag[2])*volFact)/m_; 
  out[5] += ((0.3162277660168379*Bmag[5]*tmp[9]+0.3162277660168379*tmp[5]*Bmag[9]+0.3162277660168379*Bmag[5]*tmp[7]+0.3162277660168379*tmp[5]*Bmag[7]+0.3535533905932737*Bmag[4]*tmp[6]+0.3535533905932737*tmp[4]*Bmag[6]+0.3535533905932737*Bmag[0]*tmp[5]+0.3535533905932737*tmp[0]*Bmag[5]+0.3535533905932737*Bmag[1]*tmp[3]+0.3535533905932737*tmp[1]*Bmag[3])*volFact)/m_; 
  out[6] += ((0.3162277660168379*Bmag[6]*tmp[9]+0.3162277660168379*tmp[6]*Bmag[9]+0.3162277660168379*Bmag[6]*tmp[8]+0.3162277660168379*tmp[6]*Bmag[8]+0.3535533905932737*Bmag[0]*tmp[6]+0.3535533905932737*tmp[0]*Bmag[6]+0.3535533905932737*Bmag[4]*tmp[5]+0.3535533905932737*tmp[4]*Bmag[5]+0.3535533905932737*Bmag[2]*tmp[3]+0.3535533905932737*tmp[2]*Bmag[3])*volFact)/m_; 
  out[7] += ((0.2258769757263128*Bmag[7]*tmp[7]+0.3535533905932737*Bmag[0]*tmp[7]+0.3535533905932737*tmp[0]*Bmag[7]+0.3162277660168379*Bmag[5]*tmp[5]+0.3162277660168379*Bmag[4]*tmp[4]+0.3162277660168379*Bmag[1]*tmp[1])*volFact)/m_; 
  out[8] += ((0.2258769757263128*Bmag[8]*tmp[8]+0.3535533905932737*Bmag[0]*tmp[8]+0.3535533905932737*tmp[0]*Bmag[8]+0.3162277660168379*Bmag[6]*tmp[6]+0.3162277660168379*Bmag[4]*tmp[4]+0.3162277660168379*Bmag[2]*tmp[2])*volFact)/m_; 
  out[9] += ((0.2258769757263128*Bmag[9]*tmp[9]+0.3535533905932737*Bmag[0]*tmp[9]+0.3535533905932737*tmp[0]*Bmag[9]+0.3162277660168379*Bmag[6]*tmp[6]+0.3162277660168379*Bmag[5]*tmp[5]+0.3162277660168379*Bmag[3]*tmp[3])*volFact)/m_; 
} 
void GkMomentCalc3x2vMax_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  double zeta[6]; 

  zeta[0] = 5.656854249492382*wx1_cu+1.414213562373095*dv1_sq*wx1; 
  zeta[4] = 4.898979485566357*dv1*wx1_sq+0.2449489742783178*dv1_cu; 
  out[0] += (0.3535533905932737*f[4]*zeta[4]+0.3535533905932737*f[0]*zeta[0])*volFact; 
  out[1] += 0.3535533905932737*zeta[0]*f[1]*volFact; 
  out[2] += 0.3535533905932737*zeta[0]*f[2]*volFact; 
  out[3] += 0.3535533905932737*zeta[0]*f[3]*volFact; 
} 
void GkMomentCalc3x2vMax_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  double zeta[21]; 

  zeta[0] = 5.656854249492382*wx1_cu+1.414213562373095*dv1_sq*wx1; 
  zeta[4] = 4.898979485566357*dv1*wx1_sq+0.2449489742783178*dv1_cu; 
  zeta[19] = 1.264911064067352*dv1_sq*wx1; 
  out[0] += (0.3535533905932737*f[19]*zeta[19]+0.3535533905932737*f[4]*zeta[4]+0.3535533905932737*f[0]*zeta[0])*volFact; 
  out[1] += (0.3535533905932737*zeta[4]*f[9]+0.3535533905932737*zeta[0]*f[1])*volFact; 
  out[2] += (0.3535533905932737*zeta[4]*f[10]+0.3535533905932737*zeta[0]*f[2])*volFact; 
  out[3] += (0.3535533905932737*zeta[4]*f[11]+0.3535533905932737*zeta[0]*f[3])*volFact; 
  out[4] += 0.3535533905932737*zeta[0]*f[6]*volFact; 
  out[5] += 0.3535533905932737*zeta[0]*f[7]*volFact; 
  out[6] += 0.3535533905932737*zeta[0]*f[8]*volFact; 
  out[7] += 0.3535533905932737*zeta[0]*f[16]*volFact; 
  out[8] += 0.3535533905932737*zeta[0]*f[17]*volFact; 
  out[9] += 0.3535533905932737*zeta[0]*f[18]*volFact; 
} 
void GkMomentCalc3x2vMax_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += ((0.7071067811865475*Bmag[3]*f[3]*wx1*wx2+0.7071067811865475*Bmag[2]*f[2]*wx1*wx2+0.7071067811865475*Bmag[1]*f[1]*wx1*wx2+0.7071067811865475*Bmag[0]*f[0]*wx1*wx2+0.2041241452319315*Bmag[0]*f[4]*dv1*wx2+0.2041241452319315*Bmag[0]*f[5]*dv2*wx1)*volFact)/m_; 
  out[1] += ((0.7071067811865475*Bmag[0]*f[1]*wx1*wx2+0.7071067811865475*f[0]*Bmag[1]*wx1*wx2+0.2041241452319315*Bmag[1]*f[4]*dv1*wx2+0.2041241452319315*Bmag[1]*f[5]*dv2*wx1)*volFact)/m_; 
  out[2] += ((0.7071067811865475*Bmag[0]*f[2]*wx1*wx2+0.7071067811865475*f[0]*Bmag[2]*wx1*wx2+0.2041241452319315*Bmag[2]*f[4]*dv1*wx2+0.2041241452319315*Bmag[2]*f[5]*dv2*wx1)*volFact)/m_; 
  out[3] += ((0.7071067811865475*Bmag[0]*f[3]*wx1*wx2+0.7071067811865475*f[0]*Bmag[3]*wx1*wx2+0.2041241452319315*Bmag[3]*f[4]*dv1*wx2+0.2041241452319315*Bmag[3]*f[5]*dv2*wx1)*volFact)/m_; 
} 
void GkMomentCalc3x2vMax_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  out[0] += ((0.7071067811865475*Bmag[9]*f[18]*wx1*wx2+0.7071067811865475*Bmag[8]*f[17]*wx1*wx2+0.7071067811865475*Bmag[7]*f[16]*wx1*wx2+0.7071067811865475*Bmag[6]*f[8]*wx1*wx2+0.7071067811865475*Bmag[5]*f[7]*wx1*wx2+0.7071067811865475*Bmag[4]*f[6]*wx1*wx2+0.7071067811865475*Bmag[3]*f[3]*wx1*wx2+0.7071067811865475*Bmag[2]*f[2]*wx1*wx2+0.7071067811865475*Bmag[1]*f[1]*wx1*wx2+0.7071067811865475*Bmag[0]*f[0]*wx1*wx2+0.2041241452319315*Bmag[3]*f[11]*dv1*wx2+0.2041241452319315*Bmag[2]*f[10]*dv1*wx2+0.2041241452319315*Bmag[1]*f[9]*dv1*wx2+0.2041241452319315*Bmag[0]*f[4]*dv1*wx2+0.2041241452319315*Bmag[3]*f[14]*dv2*wx1+0.2041241452319315*Bmag[2]*f[13]*dv2*wx1+0.2041241452319315*Bmag[1]*f[12]*dv2*wx1+0.2041241452319315*Bmag[0]*f[5]*dv2*wx1+0.05892556509887893*Bmag[0]*f[15]*dv1*dv2)*volFact)/m_; 
  out[1] += ((0.6324555320336759*Bmag[1]*f[16]*wx1*wx2+0.7071067811865475*Bmag[3]*f[7]*wx1*wx2+0.6324555320336759*f[1]*Bmag[7]*wx1*wx2+0.7071067811865475*Bmag[2]*f[6]*wx1*wx2+0.7071067811865475*f[3]*Bmag[5]*wx1*wx2+0.7071067811865475*f[2]*Bmag[4]*wx1*wx2+0.7071067811865475*Bmag[0]*f[1]*wx1*wx2+0.7071067811865475*f[0]*Bmag[1]*wx1*wx2+0.2041241452319315*Bmag[5]*f[11]*dv1*wx2+0.2041241452319315*Bmag[4]*f[10]*dv1*wx2+0.1825741858350554*Bmag[7]*f[9]*dv1*wx2+0.2041241452319315*Bmag[0]*f[9]*dv1*wx2+0.2041241452319315*Bmag[1]*f[4]*dv1*wx2+0.2041241452319315*Bmag[5]*f[14]*dv2*wx1+0.2041241452319315*Bmag[4]*f[13]*dv2*wx1+0.1825741858350554*Bmag[7]*f[12]*dv2*wx1+0.2041241452319315*Bmag[0]*f[12]*dv2*wx1+0.2041241452319315*Bmag[1]*f[5]*dv2*wx1+0.05892556509887893*Bmag[1]*f[15]*dv1*dv2)*volFact)/m_; 
  out[2] += ((0.6324555320336759*Bmag[2]*f[17]*wx1*wx2+0.7071067811865475*Bmag[3]*f[8]*wx1*wx2+0.6324555320336759*f[2]*Bmag[8]*wx1*wx2+0.7071067811865475*Bmag[1]*f[6]*wx1*wx2+0.7071067811865475*f[3]*Bmag[6]*wx1*wx2+0.7071067811865475*f[1]*Bmag[4]*wx1*wx2+0.7071067811865475*Bmag[0]*f[2]*wx1*wx2+0.7071067811865475*f[0]*Bmag[2]*wx1*wx2+0.2041241452319315*Bmag[6]*f[11]*dv1*wx2+0.1825741858350554*Bmag[8]*f[10]*dv1*wx2+0.2041241452319315*Bmag[0]*f[10]*dv1*wx2+0.2041241452319315*Bmag[4]*f[9]*dv1*wx2+0.2041241452319315*Bmag[2]*f[4]*dv1*wx2+0.2041241452319315*Bmag[6]*f[14]*dv2*wx1+0.1825741858350554*Bmag[8]*f[13]*dv2*wx1+0.2041241452319315*Bmag[0]*f[13]*dv2*wx1+0.2041241452319315*Bmag[4]*f[12]*dv2*wx1+0.2041241452319315*Bmag[2]*f[5]*dv2*wx1+0.05892556509887893*Bmag[2]*f[15]*dv1*dv2)*volFact)/m_; 
  out[3] += ((0.6324555320336759*Bmag[3]*f[18]*wx1*wx2+0.6324555320336759*f[3]*Bmag[9]*wx1*wx2+0.7071067811865475*Bmag[2]*f[8]*wx1*wx2+0.7071067811865475*Bmag[1]*f[7]*wx1*wx2+0.7071067811865475*f[2]*Bmag[6]*wx1*wx2+0.7071067811865475*f[1]*Bmag[5]*wx1*wx2+0.7071067811865475*Bmag[0]*f[3]*wx1*wx2+0.7071067811865475*f[0]*Bmag[3]*wx1*wx2+0.1825741858350554*Bmag[9]*f[11]*dv1*wx2+0.2041241452319315*Bmag[0]*f[11]*dv1*wx2+0.2041241452319315*Bmag[6]*f[10]*dv1*wx2+0.2041241452319315*Bmag[5]*f[9]*dv1*wx2+0.2041241452319315*Bmag[3]*f[4]*dv1*wx2+0.1825741858350554*Bmag[9]*f[14]*dv2*wx1+0.2041241452319315*Bmag[0]*f[14]*dv2*wx1+0.2041241452319315*Bmag[6]*f[13]*dv2*wx1+0.2041241452319315*Bmag[5]*f[12]*dv2*wx1+0.2041241452319315*Bmag[3]*f[5]*dv2*wx1+0.05892556509887893*Bmag[3]*f[15]*dv1*dv2)*volFact)/m_; 
  out[4] += ((0.6324555320336759*Bmag[4]*f[17]*wx1*wx2+0.6324555320336759*Bmag[4]*f[16]*wx1*wx2+0.7071067811865475*Bmag[5]*f[8]*wx1*wx2+0.6324555320336759*f[6]*Bmag[8]*wx1*wx2+0.7071067811865475*Bmag[6]*f[7]*wx1*wx2+0.6324555320336759*f[6]*Bmag[7]*wx1*wx2+0.7071067811865475*Bmag[0]*f[6]*wx1*wx2+0.7071067811865475*f[0]*Bmag[4]*wx1*wx2+0.7071067811865475*Bmag[1]*f[2]*wx1*wx2+0.7071067811865475*f[1]*Bmag[2]*wx1*wx2+0.2041241452319315*Bmag[1]*f[10]*dv1*wx2+0.2041241452319315*Bmag[2]*f[9]*dv1*wx2+0.2041241452319315*Bmag[4]*f[4]*dv1*wx2+0.2041241452319315*Bmag[1]*f[13]*dv2*wx1+0.2041241452319315*Bmag[2]*f[12]*dv2*wx1+0.2041241452319315*Bmag[4]*f[5]*dv2*wx1+0.05892556509887893*Bmag[4]*f[15]*dv1*dv2)*volFact)/m_; 
  out[5] += ((0.6324555320336759*Bmag[5]*f[18]*wx1*wx2+0.6324555320336759*Bmag[5]*f[16]*wx1*wx2+0.6324555320336759*f[7]*Bmag[9]*wx1*wx2+0.7071067811865475*Bmag[4]*f[8]*wx1*wx2+0.6324555320336759*Bmag[7]*f[7]*wx1*wx2+0.7071067811865475*Bmag[0]*f[7]*wx1*wx2+0.7071067811865475*Bmag[6]*f[6]*wx1*wx2+0.7071067811865475*f[0]*Bmag[5]*wx1*wx2+0.7071067811865475*Bmag[1]*f[3]*wx1*wx2+0.7071067811865475*f[1]*Bmag[3]*wx1*wx2+0.2041241452319315*Bmag[1]*f[11]*dv1*wx2+0.2041241452319315*Bmag[3]*f[9]*dv1*wx2+0.2041241452319315*f[4]*Bmag[5]*dv1*wx2+0.2041241452319315*Bmag[1]*f[14]*dv2*wx1+0.2041241452319315*Bmag[3]*f[12]*dv2*wx1+0.2041241452319315*Bmag[5]*f[5]*dv2*wx1+0.05892556509887893*Bmag[5]*f[15]*dv1*dv2)*volFact)/m_; 
  out[6] += ((0.6324555320336759*Bmag[6]*f[18]*wx1*wx2+0.6324555320336759*Bmag[6]*f[17]*wx1*wx2+0.6324555320336759*f[8]*Bmag[9]*wx1*wx2+0.6324555320336759*Bmag[8]*f[8]*wx1*wx2+0.7071067811865475*Bmag[0]*f[8]*wx1*wx2+0.7071067811865475*Bmag[4]*f[7]*wx1*wx2+0.7071067811865475*Bmag[5]*f[6]*wx1*wx2+0.7071067811865475*f[0]*Bmag[6]*wx1*wx2+0.7071067811865475*Bmag[2]*f[3]*wx1*wx2+0.7071067811865475*f[2]*Bmag[3]*wx1*wx2+0.2041241452319315*Bmag[2]*f[11]*dv1*wx2+0.2041241452319315*Bmag[3]*f[10]*dv1*wx2+0.2041241452319315*f[4]*Bmag[6]*dv1*wx2+0.2041241452319315*Bmag[2]*f[14]*dv2*wx1+0.2041241452319315*Bmag[3]*f[13]*dv2*wx1+0.2041241452319315*f[5]*Bmag[6]*dv2*wx1+0.05892556509887893*Bmag[6]*f[15]*dv1*dv2)*volFact)/m_; 
  out[7] += ((0.4517539514526256*Bmag[7]*f[16]*wx1*wx2+0.7071067811865475*Bmag[0]*f[16]*wx1*wx2+0.6324555320336759*Bmag[5]*f[7]*wx1*wx2+0.7071067811865475*f[0]*Bmag[7]*wx1*wx2+0.6324555320336759*Bmag[4]*f[6]*wx1*wx2+0.6324555320336759*Bmag[1]*f[1]*wx1*wx2+0.1825741858350554*Bmag[1]*f[9]*dv1*wx2+0.2041241452319315*f[4]*Bmag[7]*dv1*wx2+0.1825741858350554*Bmag[1]*f[12]*dv2*wx1+0.2041241452319315*f[5]*Bmag[7]*dv2*wx1+0.05892556509887893*Bmag[7]*f[15]*dv1*dv2)*volFact)/m_; 
  out[8] += ((0.4517539514526256*Bmag[8]*f[17]*wx1*wx2+0.7071067811865475*Bmag[0]*f[17]*wx1*wx2+0.6324555320336759*Bmag[6]*f[8]*wx1*wx2+0.7071067811865475*f[0]*Bmag[8]*wx1*wx2+0.6324555320336759*Bmag[4]*f[6]*wx1*wx2+0.6324555320336759*Bmag[2]*f[2]*wx1*wx2+0.1825741858350554*Bmag[2]*f[10]*dv1*wx2+0.2041241452319315*f[4]*Bmag[8]*dv1*wx2+0.1825741858350554*Bmag[2]*f[13]*dv2*wx1+0.2041241452319315*f[5]*Bmag[8]*dv2*wx1+0.05892556509887893*Bmag[8]*f[15]*dv1*dv2)*volFact)/m_; 
  out[9] += ((0.4517539514526256*Bmag[9]*f[18]*wx1*wx2+0.7071067811865475*Bmag[0]*f[18]*wx1*wx2+0.7071067811865475*f[0]*Bmag[9]*wx1*wx2+0.6324555320336759*Bmag[6]*f[8]*wx1*wx2+0.6324555320336759*Bmag[5]*f[7]*wx1*wx2+0.6324555320336759*Bmag[3]*f[3]*wx1*wx2+0.1825741858350554*Bmag[3]*f[11]*dv1*wx2+0.2041241452319315*f[4]*Bmag[9]*dv1*wx2+0.1825741858350554*Bmag[3]*f[14]*dv2*wx1+0.2041241452319315*f[5]*Bmag[9]*dv2*wx1+0.05892556509887893*Bmag[9]*f[15]*dv1*dv2)*volFact)/m_; 
} 
void GkMomentCalc3x2vMax_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM0[2] += 2.0*f[2]*volFact; 
  outM0[3] += 2.0*f[3]*volFact; 
  outM1[0] += 0.3333333333333333*(6.0*f[0]*wx1+1.732050807568877*f[4]*dv1)*volFact; 
  outM1[1] += 2.0*f[1]*wx1*volFact; 
  outM1[2] += 2.0*f[2]*wx1*volFact; 
  outM1[3] += 2.0*f[3]*wx1*volFact; 
  double zeta[6]; 

  zeta[0] = 5.656854249492382*wx1_sq+0.4714045207910317*dv1_sq; 
  zeta[4] = 3.265986323710906*dv1*wx1; 
  outM2[0] += (0.3535533905932737*f[4]*zeta[4]+0.3535533905932737*f[0]*zeta[0])*volFact; 
  outM2[1] += 0.3535533905932737*zeta[0]*f[1]*volFact; 
  outM2[2] += 0.3535533905932737*zeta[0]*f[2]*volFact; 
  outM2[3] += 0.3535533905932737*zeta[0]*f[3]*volFact; 
  double tmp[4]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[1] = 2.0*f[1]*wx2; 
  tmp[2] = 2.0*f[2]*wx2; 
  tmp[3] = 2.0*f[3]*wx2; 
  outM2[0] += (2.0*(0.3535533905932737*Bmag[3]*tmp[3]+0.3535533905932737*Bmag[2]*tmp[2]+0.3535533905932737*Bmag[1]*tmp[1]+0.3535533905932737*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.3535533905932737*Bmag[0]*tmp[1]+0.3535533905932737*tmp[0]*Bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.3535533905932737*Bmag[0]*tmp[2]+0.3535533905932737*tmp[0]*Bmag[2])*volFact)/m_; 
  outM2[3] += (2.0*(0.3535533905932737*Bmag[0]*tmp[3]+0.3535533905932737*tmp[0]*Bmag[3])*volFact)/m_; 
} 
void GkMomentCalc3x2vMax_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM0[2] += 2.0*f[2]*volFact; 
  outM0[3] += 2.0*f[3]*volFact; 
  outM0[4] += 2.0*f[6]*volFact; 
  outM0[5] += 2.0*f[7]*volFact; 
  outM0[6] += 2.0*f[8]*volFact; 
  outM0[7] += 2.0*f[16]*volFact; 
  outM0[8] += 2.0*f[17]*volFact; 
  outM0[9] += 2.0*f[18]*volFact; 
  outM1[0] += 0.3333333333333333*(6.0*f[0]*wx1+1.732050807568877*f[4]*dv1)*volFact; 
  outM1[1] += 0.3333333333333333*(6.0*f[1]*wx1+1.732050807568877*f[9]*dv1)*volFact; 
  outM1[2] += 0.3333333333333333*(6.0*f[2]*wx1+1.732050807568877*f[10]*dv1)*volFact; 
  outM1[3] += 0.3333333333333333*(6.0*f[3]*wx1+1.732050807568877*f[11]*dv1)*volFact; 
  outM1[4] += 2.0*f[6]*wx1*volFact; 
  outM1[5] += 2.0*f[7]*wx1*volFact; 
  outM1[6] += 2.0*f[8]*wx1*volFact; 
  outM1[7] += 2.0*f[16]*wx1*volFact; 
  outM1[8] += 2.0*f[17]*wx1*volFact; 
  outM1[9] += 2.0*f[18]*wx1*volFact; 
  outM2[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.149071198499986*f[19]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  outM2[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[9]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  outM2[2] += (2.0*f[2]*wx1_sq+1.154700538379252*f[10]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  outM2[3] += (2.0*f[3]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.1666666666666667*f[3]*dv1_sq)*volFact; 
  outM2[4] += (2.0*f[6]*wx1_sq+0.1666666666666667*f[6]*dv1_sq)*volFact; 
  outM2[5] += (2.0*f[7]*wx1_sq+0.1666666666666667*f[7]*dv1_sq)*volFact; 
  outM2[6] += (2.0*f[8]*wx1_sq+0.1666666666666667*f[8]*dv1_sq)*volFact; 
  outM2[7] += (2.0*f[16]*wx1_sq+0.1666666666666667*f[16]*dv1_sq)*volFact; 
  outM2[8] += (2.0*f[17]*wx1_sq+0.1666666666666667*f[17]*dv1_sq)*volFact; 
  outM2[9] += (2.0*f[18]*wx1_sq+0.1666666666666667*f[18]*dv1_sq)*volFact; 
  double tmp[10]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[12]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[13]*dv2; 
  tmp[3] = 2.0*f[3]*wx2+0.5773502691896258*f[14]*dv2; 
  tmp[4] = 2.0*f[6]*wx2; 
  tmp[5] = 2.0*f[7]*wx2; 
  tmp[6] = 2.0*f[8]*wx2; 
  tmp[7] = 2.0*f[16]*wx2; 
  tmp[8] = 2.0*f[17]*wx2; 
  tmp[9] = 2.0*f[18]*wx2; 
  outM2[0] += (2.0*(0.3535533905932737*Bmag[9]*tmp[9]+0.3535533905932737*Bmag[8]*tmp[8]+0.3535533905932737*Bmag[7]*tmp[7]+0.3535533905932737*Bmag[6]*tmp[6]+0.3535533905932737*Bmag[5]*tmp[5]+0.3535533905932737*Bmag[4]*tmp[4]+0.3535533905932737*Bmag[3]*tmp[3]+0.3535533905932737*Bmag[2]*tmp[2]+0.3535533905932737*Bmag[1]*tmp[1]+0.3535533905932737*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.3162277660168379*Bmag[1]*tmp[7]+0.3162277660168379*tmp[1]*Bmag[7]+0.3535533905932737*Bmag[3]*tmp[5]+0.3535533905932737*tmp[3]*Bmag[5]+0.3535533905932737*Bmag[2]*tmp[4]+0.3535533905932737*tmp[2]*Bmag[4]+0.3535533905932737*Bmag[0]*tmp[1]+0.3535533905932737*tmp[0]*Bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.3162277660168379*Bmag[2]*tmp[8]+0.3162277660168379*tmp[2]*Bmag[8]+0.3535533905932737*Bmag[3]*tmp[6]+0.3535533905932737*tmp[3]*Bmag[6]+0.3535533905932737*Bmag[1]*tmp[4]+0.3535533905932737*tmp[1]*Bmag[4]+0.3535533905932737*Bmag[0]*tmp[2]+0.3535533905932737*tmp[0]*Bmag[2])*volFact)/m_; 
  outM2[3] += (2.0*(0.3162277660168379*Bmag[3]*tmp[9]+0.3162277660168379*tmp[3]*Bmag[9]+0.3535533905932737*Bmag[2]*tmp[6]+0.3535533905932737*tmp[2]*Bmag[6]+0.3535533905932737*Bmag[1]*tmp[5]+0.3535533905932737*tmp[1]*Bmag[5]+0.3535533905932737*Bmag[0]*tmp[3]+0.3535533905932737*tmp[0]*Bmag[3])*volFact)/m_; 
  outM2[4] += (2.0*(0.3162277660168379*Bmag[4]*tmp[8]+0.3162277660168379*tmp[4]*Bmag[8]+0.3162277660168379*Bmag[4]*tmp[7]+0.3162277660168379*tmp[4]*Bmag[7]+0.3535533905932737*Bmag[5]*tmp[6]+0.3535533905932737*tmp[5]*Bmag[6]+0.3535533905932737*Bmag[0]*tmp[4]+0.3535533905932737*tmp[0]*Bmag[4]+0.3535533905932737*Bmag[1]*tmp[2]+0.3535533905932737*tmp[1]*Bmag[2])*volFact)/m_; 
  outM2[5] += (2.0*(0.3162277660168379*Bmag[5]*tmp[9]+0.3162277660168379*tmp[5]*Bmag[9]+0.3162277660168379*Bmag[5]*tmp[7]+0.3162277660168379*tmp[5]*Bmag[7]+0.3535533905932737*Bmag[4]*tmp[6]+0.3535533905932737*tmp[4]*Bmag[6]+0.3535533905932737*Bmag[0]*tmp[5]+0.3535533905932737*tmp[0]*Bmag[5]+0.3535533905932737*Bmag[1]*tmp[3]+0.3535533905932737*tmp[1]*Bmag[3])*volFact)/m_; 
  outM2[6] += (2.0*(0.3162277660168379*Bmag[6]*tmp[9]+0.3162277660168379*tmp[6]*Bmag[9]+0.3162277660168379*Bmag[6]*tmp[8]+0.3162277660168379*tmp[6]*Bmag[8]+0.3535533905932737*Bmag[0]*tmp[6]+0.3535533905932737*tmp[0]*Bmag[6]+0.3535533905932737*Bmag[4]*tmp[5]+0.3535533905932737*tmp[4]*Bmag[5]+0.3535533905932737*Bmag[2]*tmp[3]+0.3535533905932737*tmp[2]*Bmag[3])*volFact)/m_; 
  outM2[7] += (2.0*(0.2258769757263128*Bmag[7]*tmp[7]+0.3535533905932737*Bmag[0]*tmp[7]+0.3535533905932737*tmp[0]*Bmag[7]+0.3162277660168379*Bmag[5]*tmp[5]+0.3162277660168379*Bmag[4]*tmp[4]+0.3162277660168379*Bmag[1]*tmp[1])*volFact)/m_; 
  outM2[8] += (2.0*(0.2258769757263128*Bmag[8]*tmp[8]+0.3535533905932737*Bmag[0]*tmp[8]+0.3535533905932737*tmp[0]*Bmag[8]+0.3162277660168379*Bmag[6]*tmp[6]+0.3162277660168379*Bmag[4]*tmp[4]+0.3162277660168379*Bmag[2]*tmp[2])*volFact)/m_; 
  outM2[9] += (2.0*(0.2258769757263128*Bmag[9]*tmp[9]+0.3535533905932737*Bmag[0]*tmp[9]+0.3535533905932737*tmp[0]*Bmag[9]+0.3162277660168379*Bmag[6]*tmp[6]+0.3162277660168379*Bmag[5]*tmp[5]+0.3162277660168379*Bmag[3]*tmp[3])*volFact)/m_; 
} 
