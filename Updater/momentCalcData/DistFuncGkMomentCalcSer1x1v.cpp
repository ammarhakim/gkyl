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
  out[0] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1); 
} 
void GkMomentCalc1x1vSer_M1_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1); 
  out[2] += volFact*(1.414213562373095*f[4]*wx1+0.408248290463863*f[6]*dv1); 
} 
void GkMomentCalc1x1vSer_M1proj_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += 1.414213562373095*f[0]*volFact*wx1; 
  out[1] += 1.414213562373095*f[1]*volFact*wx1; 
} 
void GkMomentCalc1x1vSer_M1proj_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1); 
  out[2] += volFact*(1.414213562373095*f[4]*wx1+0.408248290463863*f[6]*dv1); 
} 
void GkMomentCalc1x1vSer_M2_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.1178511301977579*f[1]*dv1_sq); 
} 
void GkMomentCalc1x1vSer_M2_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq); 
  out[2] += volFact*(1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq); 
} 
void GkMomentCalc1x1vSer_M2par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.1178511301977579*f[1]*dv1_sq); 
} 
void GkMomentCalc1x1vSer_M2par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq); 
  out[2] += volFact*(1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq); 
} 
void GkMomentCalc1x1vSer_M3par_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1*wx1_sq+1.224744871391589*f[2]*dv1*wx1_sq+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1*wx1_sq+1.224744871391589*f[3]*dv1*wx1_sq+0.3535533905932737*f[1]*dv1_sq*wx1+0.06123724356957942*f[3]*dv1*dv1_sq); 
} 
void GkMomentCalc1x1vSer_M3par_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += volFact*(1.414213562373095*f[0]*wx1*wx1_sq+1.224744871391589*f[2]*dv1*wx1_sq+0.3162277660168379*f[5]*dv1_sq*wx1+0.3535533905932737*f[0]*dv1_sq*wx1+0.06123724356957942*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(1.414213562373095*f[1]*wx1*wx1_sq+1.224744871391589*f[3]*dv1*wx1_sq+0.3162277660168379*f[7]*dv1_sq*wx1+0.3535533905932737*f[1]*dv1_sq*wx1+0.06123724356957942*f[3]*dv1*dv1_sq); 
  out[2] += volFact*(1.414213562373095*f[4]*wx1*wx1_sq+1.224744871391589*f[6]*dv1*wx1_sq+0.3535533905932737*f[4]*dv1_sq*wx1+0.06123724356957942*f[6]*dv1*dv1_sq); 
} 
void GkMomentCalc1x1vSer_ThreeMoments_P1(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  outM0[0] += 1.414213562373095*f[0]*volFact; 
  outM0[1] += 1.414213562373095*f[1]*volFact; 
  outM1[0] += 0.2357022603955158*volFact*(6.0*f[0]*wx1+1.732050807568877*f[2]*dv1); 
  outM1[1] += 0.2357022603955158*volFact*(6.0*f[1]*wx1+1.732050807568877*f[3]*dv1); 
  outM2[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.1178511301977579*f[0]*dv1_sq); 
  outM2[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.1178511301977579*f[1]*dv1_sq); 
} 
void GkMomentCalc1x1vSer_ThreeMoments_P2(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *outM0, double *outM1, double *outM2) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  outM0[0] += 1.414213562373095*f[0]*volFact; 
  outM0[1] += 1.414213562373095*f[1]*volFact; 
  outM0[2] += 1.414213562373095*f[4]*volFact; 
  outM1[0] += 0.2357022603955158*volFact*(6.0*f[0]*wx1+1.732050807568877*f[2]*dv1); 
  outM1[1] += 0.2357022603955158*volFact*(6.0*f[1]*wx1+1.732050807568877*f[3]*dv1); 
  outM1[2] += 0.04714045207910316*volFact*(30.0*f[4]*wx1+8.660254037844387*f[6]*dv1); 
  outM2[0] += volFact*(1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq); 
  outM2[1] += volFact*(1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq); 
  outM2[2] += volFact*(1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq); 
} 
