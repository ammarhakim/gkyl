#include <math.h> 
#include <DistFuncMomentCalcModDecl.h> 
void GkMomentCalc1x2vMax_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
void GkMomentCalc1x2vMax_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
} 
void GkMomentCalc1x2vMax_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += (2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1)*volFact; 
  out[1] += 2.0*f[1]*wx1*volFact; 
} 
void GkMomentCalc1x2vMax_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += (2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1)*volFact; 
  out[1] += (2.0*f[1]*wx1+0.5773502691896258*f[4]*dv1)*volFact; 
  out[2] += 2.0*f[7]*wx1*volFact; 
} 
void GkMomentCalc1x2vMax_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double zeta[4]; 

  zeta[0] = 2.828427124746191*wx1_sq+0.2357022603955158*dv1_sq; 
  zeta[2] = 1.632993161855453*dv1*wx1; 
  out[0] += (0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  out[1] += 0.7071067811865475*zeta[0]*f[1]*volFact; 
  double tmp[2]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2; 
  tmp[1] = 2.0*f[1]*wx2; 
  out[0] += (2.0*(0.7071067811865475*Bmag[1]*tmp[1]+0.7071067811865475*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.7071067811865475*Bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*Bmag[1])*volFact)/m_; 
} 
void GkMomentCalc1x2vMax_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.149071198499986*f[8]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[7]*wx1_sq+0.1666666666666667*f[7]*dv1_sq)*volFact; 
  double tmp[3]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[2] = 2.0*f[7]*wx2; 
  out[0] += (2.0*(0.7071067811865475*Bmag[2]*tmp[2]+0.7071067811865475*Bmag[1]*tmp[1]+0.7071067811865475*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.6324555320336759*Bmag[1]*tmp[2]+0.6324555320336759*tmp[1]*Bmag[2]+0.7071067811865475*Bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.4517539514526256*Bmag[2]*tmp[2]+0.7071067811865475*Bmag[0]*tmp[2]+0.7071067811865475*tmp[0]*Bmag[2]+0.6324555320336759*Bmag[1]*tmp[1])*volFact)/m_; 
} 
void GkMomentCalc1x2vMax_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double zeta[4]; 

  zeta[0] = 2.828427124746191*wx1_sq+0.2357022603955158*dv1_sq; 
  zeta[2] = 1.632993161855453*dv1*wx1; 
  out[0] += (0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  out[1] += 0.7071067811865475*zeta[0]*f[1]*volFact; 
} 
void GkMomentCalc1x2vMax_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.149071198499986*f[8]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  out[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  out[2] += (2.0*f[7]*wx1_sq+0.1666666666666667*f[7]*dv1_sq)*volFact; 
} 
void GkMomentCalc1x2vMax_M2perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  double tmp[2]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2; 
  tmp[1] = 2.0*f[1]*wx2; 
  out[0] += ((0.7071067811865475*Bmag[1]*tmp[1]+0.7071067811865475*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.7071067811865475*Bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*Bmag[1])*volFact)/m_; 
} 
void GkMomentCalc1x2vMax_M2perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  double tmp[3]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[2] = 2.0*f[7]*wx2; 
  out[0] += ((0.7071067811865475*Bmag[2]*tmp[2]+0.7071067811865475*Bmag[1]*tmp[1]+0.7071067811865475*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.6324555320336759*Bmag[1]*tmp[2]+0.6324555320336759*tmp[1]*Bmag[2]+0.7071067811865475*Bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*Bmag[1])*volFact)/m_; 
  out[2] += ((0.4517539514526256*Bmag[2]*tmp[2]+0.7071067811865475*Bmag[0]*tmp[2]+0.7071067811865475*tmp[0]*Bmag[2]+0.6324555320336759*Bmag[1]*tmp[1])*volFact)/m_; 
} 
void GkMomentCalc1x2vMax_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  double zeta[4]; 

  zeta[0] = 2.828427124746191*wx1_cu+0.7071067811865475*dv1_sq*wx1; 
  zeta[2] = 2.449489742783178*dv1*wx1_sq+0.1224744871391589*dv1_cu; 
  out[0] += (0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  out[1] += 0.7071067811865475*zeta[0]*f[1]*volFact; 
} 
void GkMomentCalc1x2vMax_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  double zeta[10]; 

  zeta[0] = 2.828427124746191*wx1_cu+0.7071067811865475*dv1_sq*wx1; 
  zeta[2] = 2.449489742783178*dv1*wx1_sq+0.1224744871391589*dv1_cu; 
  zeta[8] = 0.6324555320336759*dv1_sq*wx1; 
  out[0] += (0.7071067811865475*f[8]*zeta[8]+0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  out[1] += (0.7071067811865475*zeta[2]*f[4]+0.7071067811865475*zeta[0]*f[1])*volFact; 
  out[2] += 0.7071067811865475*zeta[0]*f[7]*volFact; 
} 
void GkMomentCalc1x2vMax_M3perp_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += ((1.414213562373095*Bmag[1]*f[1]*wx1*wx2+1.414213562373095*Bmag[0]*f[0]*wx1*wx2+0.408248290463863*Bmag[0]*f[2]*dv1*wx2+0.408248290463863*Bmag[0]*f[3]*dv2*wx1)*volFact)/m_; 
  out[1] += ((1.414213562373095*Bmag[0]*f[1]*wx1*wx2+1.414213562373095*f[0]*Bmag[1]*wx1*wx2+0.408248290463863*Bmag[1]*f[2]*dv1*wx2+0.408248290463863*Bmag[1]*f[3]*dv2*wx1)*volFact)/m_; 
} 
void GkMomentCalc1x2vMax_M3perp_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += ((1.414213562373095*Bmag[2]*f[7]*wx1*wx2+1.414213562373095*Bmag[1]*f[1]*wx1*wx2+1.414213562373095*Bmag[0]*f[0]*wx1*wx2+0.408248290463863*Bmag[1]*f[4]*dv1*wx2+0.408248290463863*Bmag[0]*f[2]*dv1*wx2+0.408248290463863*Bmag[1]*f[5]*dv2*wx1+0.408248290463863*Bmag[0]*f[3]*dv2*wx1+0.1178511301977579*Bmag[0]*f[6]*dv1*dv2)*volFact)/m_; 
  out[1] += ((1.264911064067352*Bmag[1]*f[7]*wx1*wx2+1.264911064067352*f[1]*Bmag[2]*wx1*wx2+1.414213562373095*Bmag[0]*f[1]*wx1*wx2+1.414213562373095*f[0]*Bmag[1]*wx1*wx2+0.3651483716701108*Bmag[2]*f[4]*dv1*wx2+0.408248290463863*Bmag[0]*f[4]*dv1*wx2+0.408248290463863*Bmag[1]*f[2]*dv1*wx2+0.3651483716701108*Bmag[2]*f[5]*dv2*wx1+0.408248290463863*Bmag[0]*f[5]*dv2*wx1+0.408248290463863*Bmag[1]*f[3]*dv2*wx1+0.1178511301977579*Bmag[1]*f[6]*dv1*dv2)*volFact)/m_; 
  out[2] += ((0.9035079029052515*Bmag[2]*f[7]*wx1*wx2+1.414213562373095*Bmag[0]*f[7]*wx1*wx2+1.414213562373095*f[0]*Bmag[2]*wx1*wx2+1.264911064067352*Bmag[1]*f[1]*wx1*wx2+0.3651483716701108*Bmag[1]*f[4]*dv1*wx2+0.408248290463863*Bmag[2]*f[2]*dv1*wx2+0.3651483716701108*Bmag[1]*f[5]*dv2*wx1+0.408248290463863*Bmag[2]*f[3]*dv2*wx1+0.1178511301977579*Bmag[2]*f[6]*dv1*dv2)*volFact)/m_; 
} 
void GkMomentCalc1x2vMax_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM1[0] += 0.3333333333333333*(6.0*f[0]*wx1+1.732050807568877*f[2]*dv1)*volFact; 
  outM1[1] += 2.0*f[1]*wx1*volFact; 
  double zeta[4]; 

  zeta[0] = 2.828427124746191*wx1_sq+0.2357022603955158*dv1_sq; 
  zeta[2] = 1.632993161855453*dv1*wx1; 
  outM2[0] += (0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  outM2[1] += 0.7071067811865475*zeta[0]*f[1]*volFact; 
  double tmp[2]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2; 
  tmp[1] = 2.0*f[1]*wx2; 
  outM2[0] += (2.0*(0.7071067811865475*Bmag[1]*tmp[1]+0.7071067811865475*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.7071067811865475*Bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*Bmag[1])*volFact)/m_; 
} 
void GkMomentCalc1x2vMax_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM0[2] += 2.0*f[7]*volFact; 
  outM1[0] += 0.3333333333333333*(6.0*f[0]*wx1+1.732050807568877*f[2]*dv1)*volFact; 
  outM1[1] += 0.3333333333333333*(6.0*f[1]*wx1+1.732050807568877*f[4]*dv1)*volFact; 
  outM1[2] += 2.0*f[7]*wx1*volFact; 
  outM2[0] += (2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.149071198499986*f[8]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  outM2[1] += (2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  outM2[2] += (2.0*f[7]*wx1_sq+0.1666666666666667*f[7]*dv1_sq)*volFact; 
  double tmp[3]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2; 
  tmp[2] = 2.0*f[7]*wx2; 
  outM2[0] += (2.0*(0.7071067811865475*Bmag[2]*tmp[2]+0.7071067811865475*Bmag[1]*tmp[1]+0.7071067811865475*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.6324555320336759*Bmag[1]*tmp[2]+0.6324555320336759*tmp[1]*Bmag[2]+0.7071067811865475*Bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*Bmag[1])*volFact)/m_; 
  outM2[2] += (2.0*(0.4517539514526256*Bmag[2]*tmp[2]+0.7071067811865475*Bmag[0]*tmp[2]+0.7071067811865475*tmp[0]*Bmag[2]+0.6324555320336759*Bmag[1]*tmp[1])*volFact)/m_; 
} 
