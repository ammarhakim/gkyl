#include <math.h> 
#include <DistFuncMomentCalcModDecl.h> 
void GkMomentCalc2x2vMax_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
} 
void GkMomentCalc2x2vMax_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += 2.0*f[11]*volFact; 
  out[5] += 2.0*f[12]*volFact; 
} 
void GkMomentCalc2x2vMax_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1)*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact; 
  out[2] += 2.0*f[2]*wx1*volFact; 
} 
void GkMomentCalc2x2vMax_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1)*volFact; 
  out[1] += (2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1)*volFact; 
  out[2] += (2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1)*volFact; 
  out[3] += 2.0*f[5]*wx1*volFact; 
  out[4] += 2.0*f[11]*wx1*volFact; 
  out[5] += 2.0*f[12]*wx1*volFact; 
} 
void GkMomentCalc2x2vMax_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1_sq+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1_sq+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  double tmp[3]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2; 
  tmp[2] = 2.0*f[2]*wx2; 
  out[0] += (2.0*(0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
} 
void GkMomentCalc2x2vMax_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  out[3] += (2.0*f[5]*wx1_sq+0.1666666666666667*f[5]*dv1_sq)*volFact; 
  out[4] += (2.0*f[11]*wx1_sq+0.1666666666666667*f[11]*dv1_sq)*volFact; 
  out[5] += (2.0*f[12]*wx1_sq+0.1666666666666667*f[12]*dv1_sq)*volFact; 
  double tmp[6]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2; 
  tmp[4] = 2.0*f[11]*wx2; 
  tmp[5] = 2.0*f[12]*wx2; 
  out[0] += (2.0*(0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  out[4] += (2.0*(0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += (2.0*(0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
} 
void GkMomentCalc2x2vMax_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1_sq+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1_sq+0.1666666666666667*f[2]*dv1_sq)*volFact; 
} 
void GkMomentCalc2x2vMax_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  out[3] += (2.0*f[5]*wx1_sq+0.1666666666666667*f[5]*dv1_sq)*volFact; 
  out[4] += (2.0*f[11]*wx1_sq+0.1666666666666667*f[11]*dv1_sq)*volFact; 
  out[5] += (2.0*f[12]*wx1_sq+0.1666666666666667*f[12]*dv1_sq)*volFact; 
} 
void GkMomentCalc2x2vMax_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  double tmp[3]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2; 
  tmp[2] = 2.0*f[2]*wx2; 
  out[0] += ((0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += ((0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
} 
void GkMomentCalc2x2vMax_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  double tmp[6]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2; 
  tmp[4] = 2.0*f[11]*wx2; 
  tmp[5] = 2.0*f[12]*wx2; 
  out[0] += ((0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += ((0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  out[3] += ((0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  out[4] += ((0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  out[5] += ((0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
} 
void GkMomentCalc2x2vMax_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += (2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[3]*dv1*wx1_sq+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[3]*dv1*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1*wx1_sq+0.5*f[1]*dv1_sq*wx1)*volFact; 
  out[2] += (2.0*f[2]*wx1*wx1_sq+0.5*f[2]*dv1_sq*wx1)*volFact; 
} 
void GkMomentCalc2x2vMax_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += (2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[3]*dv1*wx1_sq+0.4472135954999579*f[13]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[3]*dv1*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[6]*dv1*wx1_sq+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[6]*dv1*dv1_sq)*volFact; 
  out[2] += (2.0*f[2]*wx1*wx1_sq+1.732050807568877*f[7]*dv1*wx1_sq+0.5*f[2]*dv1_sq*wx1+0.08660254037844387*f[7]*dv1*dv1_sq)*volFact; 
  out[3] += (2.0*f[5]*wx1*wx1_sq+0.5*f[5]*dv1_sq*wx1)*volFact; 
  out[4] += (2.0*f[11]*wx1*wx1_sq+0.5*f[11]*dv1_sq*wx1)*volFact; 
  out[5] += (2.0*f[12]*wx1*wx1_sq+0.5*f[12]*dv1_sq*wx1)*volFact; 
} 
void GkMomentCalc2x2vMax_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += ((Bmag[2]*f[2]*wx1*wx2+Bmag[1]*f[1]*wx1*wx2+Bmag[0]*f[0]*wx1*wx2+0.2886751345948129*Bmag[0]*f[3]*dv1*wx2+0.2886751345948129*Bmag[0]*f[4]*dv2*wx1)*volFact)/m_; 
  out[1] += ((Bmag[0]*f[1]*wx1*wx2+f[0]*Bmag[1]*wx1*wx2+0.2886751345948129*Bmag[1]*f[3]*dv1*wx2+0.2886751345948129*Bmag[1]*f[4]*dv2*wx1)*volFact)/m_; 
  out[2] += ((Bmag[0]*f[2]*wx1*wx2+f[0]*Bmag[2]*wx1*wx2+0.2886751345948129*Bmag[2]*f[3]*dv1*wx2+0.2886751345948129*Bmag[2]*f[4]*dv2*wx1)*volFact)/m_; 
} 
void GkMomentCalc2x2vMax_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += ((Bmag[5]*f[12]*wx1*wx2+Bmag[4]*f[11]*wx1*wx2+Bmag[3]*f[5]*wx1*wx2+Bmag[2]*f[2]*wx1*wx2+Bmag[1]*f[1]*wx1*wx2+Bmag[0]*f[0]*wx1*wx2+0.2886751345948129*Bmag[2]*f[7]*dv1*wx2+0.2886751345948129*Bmag[1]*f[6]*dv1*wx2+0.2886751345948129*Bmag[0]*f[3]*dv1*wx2+0.2886751345948129*Bmag[2]*f[9]*dv2*wx1+0.2886751345948129*Bmag[1]*f[8]*dv2*wx1+0.2886751345948129*Bmag[0]*f[4]*dv2*wx1+0.08333333333333333*Bmag[0]*f[10]*dv1*dv2)*volFact)/m_; 
  out[1] += ((0.8944271909999159*Bmag[1]*f[11]*wx1*wx2+Bmag[2]*f[5]*wx1*wx2+0.8944271909999159*f[1]*Bmag[4]*wx1*wx2+f[2]*Bmag[3]*wx1*wx2+Bmag[0]*f[1]*wx1*wx2+f[0]*Bmag[1]*wx1*wx2+0.2886751345948129*Bmag[3]*f[7]*dv1*wx2+0.2581988897471612*Bmag[4]*f[6]*dv1*wx2+0.2886751345948129*Bmag[0]*f[6]*dv1*wx2+0.2886751345948129*Bmag[1]*f[3]*dv1*wx2+0.2886751345948129*Bmag[3]*f[9]*dv2*wx1+0.2581988897471612*Bmag[4]*f[8]*dv2*wx1+0.2886751345948129*Bmag[0]*f[8]*dv2*wx1+0.2886751345948129*Bmag[1]*f[4]*dv2*wx1+0.08333333333333333*Bmag[1]*f[10]*dv1*dv2)*volFact)/m_; 
  out[2] += ((0.8944271909999159*Bmag[2]*f[12]*wx1*wx2+Bmag[1]*f[5]*wx1*wx2+0.8944271909999159*f[2]*Bmag[5]*wx1*wx2+f[1]*Bmag[3]*wx1*wx2+Bmag[0]*f[2]*wx1*wx2+f[0]*Bmag[2]*wx1*wx2+0.2581988897471612*Bmag[5]*f[7]*dv1*wx2+0.2886751345948129*Bmag[0]*f[7]*dv1*wx2+0.2886751345948129*Bmag[3]*f[6]*dv1*wx2+0.2886751345948129*Bmag[2]*f[3]*dv1*wx2+0.2581988897471612*Bmag[5]*f[9]*dv2*wx1+0.2886751345948129*Bmag[0]*f[9]*dv2*wx1+0.2886751345948129*Bmag[3]*f[8]*dv2*wx1+0.2886751345948129*Bmag[2]*f[4]*dv2*wx1+0.08333333333333333*Bmag[2]*f[10]*dv1*dv2)*volFact)/m_; 
  out[3] += ((0.8944271909999159*Bmag[3]*f[12]*wx1*wx2+0.8944271909999159*Bmag[3]*f[11]*wx1*wx2+0.8944271909999159*Bmag[5]*f[5]*wx1*wx2+0.8944271909999159*Bmag[4]*f[5]*wx1*wx2+Bmag[0]*f[5]*wx1*wx2+f[0]*Bmag[3]*wx1*wx2+Bmag[1]*f[2]*wx1*wx2+f[1]*Bmag[2]*wx1*wx2+0.2886751345948129*Bmag[1]*f[7]*dv1*wx2+0.2886751345948129*Bmag[2]*f[6]*dv1*wx2+0.2886751345948129*Bmag[3]*f[3]*dv1*wx2+0.2886751345948129*Bmag[1]*f[9]*dv2*wx1+0.2886751345948129*Bmag[2]*f[8]*dv2*wx1+0.2886751345948129*Bmag[3]*f[4]*dv2*wx1+0.08333333333333333*Bmag[3]*f[10]*dv1*dv2)*volFact)/m_; 
  out[4] += ((0.6388765649999399*Bmag[4]*f[11]*wx1*wx2+Bmag[0]*f[11]*wx1*wx2+0.8944271909999159*Bmag[3]*f[5]*wx1*wx2+f[0]*Bmag[4]*wx1*wx2+0.8944271909999159*Bmag[1]*f[1]*wx1*wx2+0.2581988897471612*Bmag[1]*f[6]*dv1*wx2+0.2886751345948129*f[3]*Bmag[4]*dv1*wx2+0.2581988897471612*Bmag[1]*f[8]*dv2*wx1+0.2886751345948129*Bmag[4]*f[4]*dv2*wx1+0.08333333333333333*Bmag[4]*f[10]*dv1*dv2)*volFact)/m_; 
  out[5] += ((0.6388765649999399*Bmag[5]*f[12]*wx1*wx2+Bmag[0]*f[12]*wx1*wx2+0.8944271909999159*Bmag[3]*f[5]*wx1*wx2+f[0]*Bmag[5]*wx1*wx2+0.8944271909999159*Bmag[2]*f[2]*wx1*wx2+0.2581988897471612*Bmag[2]*f[7]*dv1*wx2+0.2886751345948129*f[3]*Bmag[5]*dv1*wx2+0.2581988897471612*Bmag[2]*f[9]*dv2*wx1+0.2886751345948129*f[4]*Bmag[5]*dv2*wx1+0.08333333333333333*Bmag[5]*f[10]*dv1*dv2)*volFact)/m_; 
} 
void GkMomentCalc2x2vMax_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM0[2] += 2.0*f[2]*volFact; 
  outM1[0] += 0.3333333333333333*(6.0*f[0]*wx1+1.732050807568877*f[3]*dv1)*volFact; 
  outM1[1] += 2.0*f[1]*wx1*volFact; 
  outM1[2] += 2.0*f[2]*wx1*volFact; 
  outM2[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  outM2[1] += (2.0*f[1]*wx1_sq+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  outM2[2] += (2.0*f[2]*wx1_sq+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  double tmp[3]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2; 
  tmp[2] = 2.0*f[2]*wx2; 
  outM2[0] += (2.0*(0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
} 
void GkMomentCalc2x2vMax_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM0[2] += 2.0*f[2]*volFact; 
  outM0[3] += 2.0*f[5]*volFact; 
  outM0[4] += 2.0*f[11]*volFact; 
  outM0[5] += 2.0*f[12]*volFact; 
  outM1[0] += 0.3333333333333333*(6.0*f[0]*wx1+1.732050807568877*f[3]*dv1)*volFact; 
  outM1[1] += 0.3333333333333333*(6.0*f[1]*wx1+1.732050807568877*f[6]*dv1)*volFact; 
  outM1[2] += 0.3333333333333333*(6.0*f[2]*wx1+1.732050807568877*f[7]*dv1)*volFact; 
  outM1[3] += 2.0*f[5]*wx1*volFact; 
  outM1[4] += 2.0*f[11]*wx1*volFact; 
  outM1[5] += 2.0*f[12]*wx1*volFact; 
  outM2[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[13]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  outM2[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  outM2[2] += (2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.1666666666666667*f[2]*dv1_sq)*volFact; 
  outM2[3] += (2.0*f[5]*wx1_sq+0.1666666666666667*f[5]*dv1_sq)*volFact; 
  outM2[4] += (2.0*f[11]*wx1_sq+0.1666666666666667*f[11]*dv1_sq)*volFact; 
  outM2[5] += (2.0*f[12]*wx1_sq+0.1666666666666667*f[12]*dv1_sq)*volFact; 
  double tmp[6]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2; 
  tmp[4] = 2.0*f[11]*wx2; 
  tmp[5] = 2.0*f[12]*wx2; 
  outM2[0] += (2.0*(0.5*Bmag[5]*tmp[5]+0.5*Bmag[4]*tmp[4]+0.5*Bmag[3]*tmp[3]+0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.4472135954999579*Bmag[1]*tmp[4]+0.4472135954999579*tmp[1]*Bmag[4]+0.5*Bmag[2]*tmp[3]+0.5*tmp[2]*Bmag[3]+0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.4472135954999579*Bmag[2]*tmp[5]+0.4472135954999579*tmp[2]*Bmag[5]+0.5*Bmag[1]*tmp[3]+0.5*tmp[1]*Bmag[3]+0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
  outM2[3] += (2.0*(0.4472135954999579*Bmag[3]*tmp[5]+0.4472135954999579*tmp[3]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[4]+0.4472135954999579*tmp[3]*Bmag[4]+0.5*Bmag[0]*tmp[3]+0.5*tmp[0]*Bmag[3]+0.5*Bmag[1]*tmp[2]+0.5*tmp[1]*Bmag[2])*volFact)/m_; 
  outM2[4] += (2.0*(0.31943828249997*Bmag[4]*tmp[4]+0.5*Bmag[0]*tmp[4]+0.5*tmp[0]*Bmag[4]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[1]*tmp[1])*volFact)/m_; 
  outM2[5] += (2.0*(0.31943828249997*Bmag[5]*tmp[5]+0.5*Bmag[0]*tmp[5]+0.5*tmp[0]*Bmag[5]+0.4472135954999579*Bmag[3]*tmp[3]+0.4472135954999579*Bmag[2]*tmp[2])*volFact)/m_; 
} 
void GkMomentCalc2x2vMax_StarMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM0[2] += 2.0*f[2]*volFact; 
  outM1[0] += 2.0*f[0]*wx1*volFact; 
  outM1[1] += 2.0*f[1]*wx1*volFact; 
  outM1[2] += 2.0*f[2]*wx1*volFact; 
  outM2[0] += (2.0*f[0]*wx1_sq+0.5773502691896258*f[3]*dv1*wx1)*volFact; 
  outM2[1] += 2.0*f[1]*wx1_sq*volFact; 
  outM2[2] += 2.0*f[2]*wx1_sq*volFact; 
  double tmp[3]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2; 
  tmp[2] = 2.0*f[2]*wx2; 
  outM2[0] += (2.0*(0.5*Bmag[2]*tmp[2]+0.5*Bmag[1]*tmp[1]+0.5*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.5*Bmag[0]*tmp[1]+0.5*tmp[0]*Bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.5*Bmag[0]*tmp[2]+0.5*tmp[0]*Bmag[2])*volFact)/m_; 
} 
