#include <math.h> 
#include <DistFuncMomentCalcModDecl.h> 
void GkMomentCalc1x1vMax_M0_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
} 
void GkMomentCalc1x1vMax_M0_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact; 
} 
void GkMomentCalc1x1vMax_M1_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += (1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1)*volFact; 
  out[1] += 1.414213562373095*f[1]*wx1*volFact; 
} 
void GkMomentCalc1x1vMax_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += (1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1)*volFact; 
  out[2] += 1.414213562373095*f[4]*wx1*volFact; 
} 
void GkMomentCalc1x1vMax_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact; 
} 
void GkMomentCalc1x1vMax_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.1178511301977579*f[1]*dv1_sq)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1_sq+0.1178511301977579*f[4]*dv1_sq)*volFact; 
} 
void GkMomentCalc1x1vMax_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact; 
} 
void GkMomentCalc1x1vMax_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.1178511301977579*f[1]*dv1_sq)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1_sq+0.1178511301977579*f[4]*dv1_sq)*volFact; 
} 
void GkMomentCalc1x1vMax_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1*wx1_sq+1.224744871391589*f[2]*dv1*wx1_sq+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[2]*dv1*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1*wx1_sq+0.3535533905932737*f[1]*dv1_sq*wx1)*volFact; 
} 
void GkMomentCalc1x1vMax_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1*wx1_sq+1.224744871391589*f[2]*dv1*wx1_sq+0.3162277660168379*f[5]*dv1_sq*wx1+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[2]*dv1*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1*wx1_sq+1.224744871391589*f[3]*dv1*wx1_sq+0.3535533905932737*f[1]*dv1_sq*wx1+0.06123724356957942*f[3]*dv1*dv1_sq)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1*wx1_sq+0.3535533905932737*f[4]*dv1_sq*wx1)*volFact; 
} 
void GkMomentCalc1x1vMax_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  outM0[0] += 1.414213562373095*f[0]*volFact; 
  outM0[1] += 1.414213562373095*f[1]*volFact; 
  outM1[0] += 0.2357022603955158*(6.0*f[0]*wx1+1.732050807568877*f[2]*dv1)*volFact; 
  outM1[1] += 1.414213562373095*f[1]*wx1*volFact; 
  outM2[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  outM2[1] += (1.414213562373095*f[1]*wx1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact; 
} 
void GkMomentCalc1x1vMax_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  outM0[0] += 1.414213562373095*f[0]*volFact; 
  outM0[1] += 1.414213562373095*f[1]*volFact; 
  outM0[2] += 1.414213562373095*f[4]*volFact; 
  outM1[0] += 0.2357022603955158*(6.0*f[0]*wx1+1.732050807568877*f[2]*dv1)*volFact; 
  outM1[1] += 0.2357022603955158*(6.0*f[1]*wx1+1.732050807568877*f[3]*dv1)*volFact; 
  outM1[2] += 1.414213562373095*f[4]*wx1*volFact; 
  outM2[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  outM2[1] += (1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.1178511301977579*f[1]*dv1_sq)*volFact; 
  outM2[2] += (1.414213562373095*f[4]*wx1_sq+0.1178511301977579*f[4]*dv1_sq)*volFact; 
} 
void GkMomentCalc1x1vMax_StarMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  outM0[0] += 1.414213562373095*f[0]*volFact; 
  outM0[1] += 1.414213562373095*f[1]*volFact; 
  outM1[0] += 1.414213562373095*f[0]*wx1*volFact; 
  outM1[1] += 1.414213562373095*f[1]*wx1*volFact; 
  outM2[0] += (1.414213562373095*f[0]*wx1_sq+0.408248290463863*f[2]*dv1*wx1)*volFact; 
  outM2[1] += 1.414213562373095*f[1]*wx1_sq*volFact; 
} 
