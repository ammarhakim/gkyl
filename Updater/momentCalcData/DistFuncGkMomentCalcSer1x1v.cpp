#include <math.h> 
#include <DistFuncMomentCalcModDecl.h> 
void GkMomentCalc1x1vSer_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
} 
void GkMomentCalc1x1vSer_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact; 
} 
void GkMomentCalc1x1vSer_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += (1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1)*volFact; 
} 
void GkMomentCalc1x1vSer_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += (1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1+0.408248290463863*f[6]*dv1)*volFact; 
} 
void GkMomentCalc1x1vSer_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double zeta[4]; 

  zeta[0] = 2.0*wx1_sq+0.1666666666666667*dv1_sq; 
  zeta[2] = 1.154700538379252*dv1*wx1; 
  out[0] += (0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  out[1] += (0.7071067811865475*zeta[2]*f[3]+0.7071067811865475*zeta[0]*f[1])*volFact; 
} 
void GkMomentCalc1x1vSer_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq)*volFact; 
} 
void GkMomentCalc1x1vSer_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double zeta[4]; 

  zeta[0] = 2.0*wx1_sq+0.1666666666666667*dv1_sq; 
  zeta[2] = 1.154700538379252*dv1*wx1; 
  out[0] += (0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  out[1] += (0.7071067811865475*zeta[2]*f[3]+0.7071067811865475*zeta[0]*f[1])*volFact; 
} 
void GkMomentCalc1x1vSer_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq)*volFact; 
} 
void GkMomentCalc1x1vSer_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  double zeta[4]; 

  zeta[0] = 2.0*wx1_cu+0.5*dv1_sq*wx1; 
  zeta[2] = 1.732050807568877*dv1*wx1_sq+0.08660254037844387*dv1_cu; 
  out[0] += (0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  out[1] += (0.7071067811865475*zeta[2]*f[3]+0.7071067811865475*zeta[0]*f[1])*volFact; 
} 
void GkMomentCalc1x1vSer_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  double zeta[8]; 

  zeta[0] = 2.0*wx1_cu+0.5*dv1_sq*wx1; 
  zeta[2] = 1.732050807568877*dv1*wx1_sq+0.08660254037844387*dv1_cu; 
  zeta[5] = 0.4472135954999579*dv1_sq*wx1; 
  out[0] += (0.7071067811865475*f[5]*zeta[5]+0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  out[1] += (0.7071067811865475*zeta[5]*f[7]+0.7071067811865475*zeta[2]*f[3]+0.7071067811865475*zeta[0]*f[1])*volFact; 
  out[2] += (0.7071067811865475*zeta[2]*f[6]+0.7071067811865475*zeta[0]*f[4])*volFact; 
} 
void GkMomentCalc1x1vSer_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  outM0[0] += 1.414213562373095*f[0]*volFact; 
  outM0[1] += 1.414213562373095*f[1]*volFact; 
  outM1[0] += 0.2357022603955158*(6.0*f[0]*wx1+1.732050807568877*f[2]*dv1)*volFact; 
  outM1[1] += 0.2357022603955158*(6.0*f[1]*wx1+1.732050807568877*f[3]*dv1)*volFact; 
  double zeta[4]; 

  zeta[0] = 2.0*wx1_sq+0.1666666666666667*dv1_sq; 
  zeta[2] = 1.154700538379252*dv1*wx1; 
  outM2[0] += (0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  outM2[1] += (0.7071067811865475*zeta[2]*f[3]+0.7071067811865475*zeta[0]*f[1])*volFact; 
} 
void GkMomentCalc1x1vSer_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  outM0[0] += 1.414213562373095*f[0]*volFact; 
  outM0[1] += 1.414213562373095*f[1]*volFact; 
  outM0[2] += 1.414213562373095*f[4]*volFact; 
  outM1[0] += 0.2357022603955158*(6.0*f[0]*wx1+1.732050807568877*f[2]*dv1)*volFact; 
  outM1[1] += 0.2357022603955158*(6.0*f[1]*wx1+1.732050807568877*f[3]*dv1)*volFact; 
  outM1[2] += 0.04714045207910316*(30.0*f[4]*wx1+8.660254037844387*f[6]*dv1)*volFact; 
  outM2[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  outM2[1] += (1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact; 
  outM2[2] += (1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq)*volFact; 
} 
