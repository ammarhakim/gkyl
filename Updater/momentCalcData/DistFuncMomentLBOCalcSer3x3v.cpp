#include <DistFuncMomentCalcModDecl.h> 
void VmM0Star3x3vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[6]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[4]*dxvl[5]*(wr[3]-wl[3]); 
 
  out[0] += ((-0.8164965809277261*fr[4])+0.8164965809277261*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[10])+0.8164965809277261*fl[10]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += ((-0.8164965809277261*fr[11])+0.8164965809277261*fl[11]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
  out[3] += ((-0.8164965809277261*fr[12])+0.8164965809277261*fl[12]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3])*dS; 
  out[4] += ((-0.8164965809277261*fr[23])+0.8164965809277261*fl[23]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7])*dS; 
  out[5] += ((-0.8164965809277261*fr[24])+0.8164965809277261*fl[24]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8])*dS; 
  out[6] += ((-0.8164965809277261*fr[25])+0.8164965809277261*fl[25]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9])*dS; 
  out[7] += ((-0.8164965809277261*fr[42])+0.8164965809277261*fl[42]+0.7071067811865475*fr[22]+0.7071067811865475*fl[22])*dS; 
 
} 
 
void VmM0Star3x3vSer_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[6]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[3]*dxvl[5]*(wr[4]-wl[4]); 
 
  out[0] += ((-0.8164965809277261*fr[5])+0.8164965809277261*fl[5]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[13])+0.8164965809277261*fl[13]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += ((-0.8164965809277261*fr[14])+0.8164965809277261*fl[14]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
  out[3] += ((-0.8164965809277261*fr[15])+0.8164965809277261*fl[15]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3])*dS; 
  out[4] += ((-0.8164965809277261*fr[26])+0.8164965809277261*fl[26]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7])*dS; 
  out[5] += ((-0.8164965809277261*fr[27])+0.8164965809277261*fl[27]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8])*dS; 
  out[6] += ((-0.8164965809277261*fr[28])+0.8164965809277261*fl[28]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9])*dS; 
  out[7] += ((-0.8164965809277261*fr[43])+0.8164965809277261*fl[43]+0.7071067811865475*fr[22]+0.7071067811865475*fl[22])*dS; 
 
} 
 
void VmM0Star3x3vSer_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[6]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[3]*dxvl[4]*(wr[5]-wl[5]); 
 
  out[0] += ((-0.8164965809277261*fr[6])+0.8164965809277261*fl[6]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[17])+0.8164965809277261*fl[17]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += ((-0.8164965809277261*fr[18])+0.8164965809277261*fl[18]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
  out[3] += ((-0.8164965809277261*fr[19])+0.8164965809277261*fl[19]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3])*dS; 
  out[4] += ((-0.8164965809277261*fr[32])+0.8164965809277261*fl[32]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7])*dS; 
  out[5] += ((-0.8164965809277261*fr[33])+0.8164965809277261*fl[33]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8])*dS; 
  out[6] += ((-0.8164965809277261*fr[34])+0.8164965809277261*fl[34]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9])*dS; 
  out[7] += ((-0.8164965809277261*fr[47])+0.8164965809277261*fl[47]+0.7071067811865475*fr[22]+0.7071067811865475*fl[22])*dS; 
 
} 
 
void VmM1iM2Star3x3vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2) 
{ 
  // w[6]:    Cell-center coordinates. 
  // dxv[6]:  Cell length in each direciton. 
  // intFac:  for gyrokinetics (not used in Vm). 
  // m_:      mass (not used in Vm). 
  // Bmag[8]: Magnetic field magnitude (not used in Vm). 
  // f:       Distribution function. 
  // outM1i:  Contribution to M_1^star from this cell. 
  // outM2:   Contribution to M_2^star from this cell. 
 
  const double volFact = 0.125*dxv[3]*dxv[4]*dxv[5]; 
  double wvSq[3]; 
  wvSq[0]  = w[3]*w[3]; 
  wvSq[1]  = w[4]*w[4]; 
  wvSq[2]  = w[5]*w[5]; 
  double dvSq[3]; 
  dvSq[0] = dxv[3]*dxv[3]; 
  dvSq[1] = dxv[4]*dxv[4]; 
  dvSq[2] = dxv[5]*dxv[5]; 
  double tempM0[8]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 
  tempM0[2] = 2.828427124746191*f[2]*volFact; 
  tempM0[3] = 2.828427124746191*f[3]*volFact; 
  tempM0[4] = 2.828427124746191*f[7]*volFact; 
  tempM0[5] = 2.828427124746191*f[8]*volFact; 
  tempM0[6] = 2.828427124746191*f[9]*volFact; 
  tempM0[7] = 2.828427124746191*f[22]*volFact; 

  outM1i[0] += tempM0[0]*w[3]; 
  outM1i[1] += tempM0[1]*w[3]; 
  outM1i[2] += tempM0[2]*w[3]; 
  outM1i[3] += tempM0[3]*w[3]; 
  outM1i[4] += w[3]*tempM0[4]; 
  outM1i[5] += w[3]*tempM0[5]; 
  outM1i[6] += w[3]*tempM0[6]; 
  outM1i[7] += w[3]*tempM0[7]; 
  outM1i[8] += tempM0[0]*w[4]; 
  outM1i[9] += tempM0[1]*w[4]; 
  outM1i[10] += tempM0[2]*w[4]; 
  outM1i[11] += tempM0[3]*w[4]; 
  outM1i[12] += tempM0[4]*w[4]; 
  outM1i[13] += w[4]*tempM0[5]; 
  outM1i[14] += w[4]*tempM0[6]; 
  outM1i[15] += w[4]*tempM0[7]; 
  outM1i[16] += tempM0[0]*w[5]; 
  outM1i[17] += tempM0[1]*w[5]; 
  outM1i[18] += tempM0[2]*w[5]; 
  outM1i[19] += tempM0[3]*w[5]; 
  outM1i[20] += tempM0[4]*w[5]; 
  outM1i[21] += tempM0[5]*w[5]; 
  outM1i[22] += w[5]*tempM0[6]; 
  outM1i[23] += w[5]*tempM0[7]; 

  outM2[0] += (0.8164965809277261*dxv[5]*w[5]*f[6]+0.8164965809277261*dxv[4]*w[4]*f[5]+0.8164965809277261*dxv[3]*w[3]*f[4])*volFact+tempM0[0]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[1] += (0.8164965809277261*dxv[5]*w[5]*f[17]+0.8164965809277261*dxv[4]*w[4]*f[13]+0.8164965809277261*dxv[3]*w[3]*f[10])*volFact+tempM0[1]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[2] += (0.8164965809277261*dxv[5]*w[5]*f[18]+0.8164965809277261*dxv[4]*w[4]*f[14]+0.8164965809277261*dxv[3]*w[3]*f[11])*volFact+tempM0[2]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[3] += (0.8164965809277261*dxv[5]*w[5]*f[19]+0.8164965809277261*dxv[4]*w[4]*f[15]+0.8164965809277261*dxv[3]*w[3]*f[12])*volFact+tempM0[3]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[4] += (0.8164965809277261*dxv[5]*w[5]*f[32]+0.8164965809277261*dxv[4]*w[4]*f[26]+0.8164965809277261*dxv[3]*w[3]*f[23])*volFact+tempM0[4]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[5] += (0.8164965809277261*dxv[5]*w[5]*f[33]+0.8164965809277261*dxv[4]*w[4]*f[27]+0.8164965809277261*dxv[3]*w[3]*f[24])*volFact+tempM0[5]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[6] += (0.8164965809277261*dxv[5]*w[5]*f[34]+0.8164965809277261*dxv[4]*w[4]*f[28]+0.8164965809277261*dxv[3]*w[3]*f[25])*volFact+(wvSq[2]+wvSq[1]+wvSq[0])*tempM0[6]; 
  outM2[7] += (0.8164965809277261*dxv[5]*w[5]*f[47]+0.8164965809277261*dxv[4]*w[4]*f[43]+0.8164965809277261*dxv[3]*w[3]*f[42])*volFact+(wvSq[2]+wvSq[1]+wvSq[0])*tempM0[7]; 
 
} 
void VmBoundaryIntegral3x3vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[64]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[4]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[11]-1.414213562373095*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[3])*dS; 
    out[4] += (2.449489742783178*fIn[23]-1.414213562373095*fIn[7])*dS; 
    out[5] += (2.449489742783178*fIn[24]-1.414213562373095*fIn[8])*dS; 
    out[6] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[9])*dS; 
    out[7] += (2.449489742783178*fIn[42]-1.414213562373095*fIn[22])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[11]+1.414213562373095*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[3])*dS; 
    out[4] += (2.449489742783178*fIn[23]+1.414213562373095*fIn[7])*dS; 
    out[5] += (2.449489742783178*fIn[24]+1.414213562373095*fIn[8])*dS; 
    out[6] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[9])*dS; 
    out[7] += (2.449489742783178*fIn[42]+1.414213562373095*fIn[22])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[256]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[4]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[25])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[1] += ((-3.16227766016838*fIn[57])+2.449489742783178*fIn[10]-1.414213562373095*fIn[1])*dS; 
    out[2] += ((-3.16227766016838*fIn[58])+2.449489742783178*fIn[11]-1.414213562373095*fIn[2])*dS; 
    out[3] += ((-3.16227766016838*fIn[59])+2.449489742783178*fIn[12]-1.414213562373095*fIn[3])*dS; 
    out[4] += ((-3.16227766016838*fIn[102])+2.449489742783178*fIn[29]-1.414213562373095*fIn[7])*dS; 
    out[5] += ((-3.16227766016838*fIn[103])+2.449489742783178*fIn[30]-1.414213562373095*fIn[8])*dS; 
    out[6] += ((-3.16227766016838*fIn[104])+2.449489742783178*fIn[31]-1.414213562373095*fIn[9])*dS; 
    out[7] += (2.449489742783178*fIn[54]-1.414213562373095*fIn[22])*dS; 
    out[8] += (2.449489742783178*fIn[55]-1.414213562373095*fIn[23])*dS; 
    out[9] += (2.449489742783178*fIn[56]-1.414213562373095*fIn[24])*dS; 
    out[10] += ((-3.16227766016838*fIn[162])+2.449489742783178*fIn[78]-1.414213562373095*fIn[28])*dS; 
    out[11] += (2.449489742783178*fIn[96]-1.414213562373095*fIn[48])*dS; 
    out[12] += (2.449489742783178*fIn[97]-1.414213562373095*fIn[49])*dS; 
    out[13] += (2.449489742783178*fIn[98]-1.414213562373095*fIn[50])*dS; 
    out[14] += (2.449489742783178*fIn[99]-1.414213562373095*fIn[51])*dS; 
    out[15] += (2.449489742783178*fIn[100]-1.414213562373095*fIn[52])*dS; 
    out[16] += (2.449489742783178*fIn[101]-1.414213562373095*fIn[53])*dS; 
    out[17] += (2.449489742783178*fIn[159]-1.414213562373095*fIn[93])*dS; 
    out[18] += (2.449489742783178*fIn[160]-1.414213562373095*fIn[94])*dS; 
    out[19] += (2.449489742783178*fIn[161]-1.414213562373095*fIn[95])*dS; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[25]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[1] += (3.16227766016838*fIn[57]+2.449489742783178*fIn[10]+1.414213562373095*fIn[1])*dS; 
    out[2] += (3.16227766016838*fIn[58]+2.449489742783178*fIn[11]+1.414213562373095*fIn[2])*dS; 
    out[3] += (3.16227766016838*fIn[59]+2.449489742783178*fIn[12]+1.414213562373095*fIn[3])*dS; 
    out[4] += (3.16227766016838*fIn[102]+2.449489742783178*fIn[29]+1.414213562373095*fIn[7])*dS; 
    out[5] += (3.16227766016838*fIn[103]+2.449489742783178*fIn[30]+1.414213562373095*fIn[8])*dS; 
    out[6] += (3.16227766016838*fIn[104]+2.449489742783178*fIn[31]+1.414213562373095*fIn[9])*dS; 
    out[7] += (2.449489742783178*fIn[54]+1.414213562373095*fIn[22])*dS; 
    out[8] += (2.449489742783178*fIn[55]+1.414213562373095*fIn[23])*dS; 
    out[9] += (2.449489742783178*fIn[56]+1.414213562373095*fIn[24])*dS; 
    out[10] += (3.16227766016838*fIn[162]+2.449489742783178*fIn[78]+1.414213562373095*fIn[28])*dS; 
    out[11] += (2.449489742783178*fIn[96]+1.414213562373095*fIn[48])*dS; 
    out[12] += (2.449489742783178*fIn[97]+1.414213562373095*fIn[49])*dS; 
    out[13] += (2.449489742783178*fIn[98]+1.414213562373095*fIn[50])*dS; 
    out[14] += (2.449489742783178*fIn[99]+1.414213562373095*fIn[51])*dS; 
    out[15] += (2.449489742783178*fIn[100]+1.414213562373095*fIn[52])*dS; 
    out[16] += (2.449489742783178*fIn[101]+1.414213562373095*fIn[53])*dS; 
    out[17] += (2.449489742783178*fIn[159]+1.414213562373095*fIn[93])*dS; 
    out[18] += (2.449489742783178*fIn[160]+1.414213562373095*fIn[94])*dS; 
    out[19] += (2.449489742783178*fIn[161]+1.414213562373095*fIn[95])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[64]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[4]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[4]-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[10]-0.7071067811865475*fIn[1]*dxv[3])*dS; 
    out[2] += (2.449489742783178*fIn[11]-1.414213562373095*fIn[2])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[11]-0.7071067811865475*fIn[2]*dxv[3])*dS; 
    out[3] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[3])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[12]-0.7071067811865475*dxv[3]*fIn[3])*dS; 
    out[4] += (2.449489742783178*fIn[23]-1.414213562373095*fIn[7])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[23]-0.7071067811865475*dxv[3]*fIn[7])*dS; 
    out[5] += (2.449489742783178*fIn[24]-1.414213562373095*fIn[8])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[24]-0.7071067811865475*dxv[3]*fIn[8])*dS; 
    out[6] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[9])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[25]-0.7071067811865475*dxv[3]*fIn[9])*dS; 
    out[7] += (2.449489742783178*fIn[42]-1.414213562373095*fIn[22])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[42]-0.7071067811865475*dxv[3]*fIn[22])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[4])-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[10])-0.7071067811865475*fIn[1]*dxv[3])*dS; 
    out[2] += (2.449489742783178*fIn[11]+1.414213562373095*fIn[2])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[11])-0.7071067811865475*fIn[2]*dxv[3])*dS; 
    out[3] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[3])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[12])-0.7071067811865475*dxv[3]*fIn[3])*dS; 
    out[4] += (2.449489742783178*fIn[23]+1.414213562373095*fIn[7])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[23])-0.7071067811865475*dxv[3]*fIn[7])*dS; 
    out[5] += (2.449489742783178*fIn[24]+1.414213562373095*fIn[8])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[24])-0.7071067811865475*dxv[3]*fIn[8])*dS; 
    out[6] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[9])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[25])-0.7071067811865475*dxv[3]*fIn[9])*dS; 
    out[7] += (2.449489742783178*fIn[42]+1.414213562373095*fIn[22])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[42])-0.7071067811865475*dxv[3]*fIn[22])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[256]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[4]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[25])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[57])+2.449489742783178*fIn[10]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += ((-3.16227766016838*fIn[58])+2.449489742783178*fIn[11]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += ((-3.16227766016838*fIn[59])+2.449489742783178*fIn[12]-1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += ((-3.16227766016838*fIn[102])+2.449489742783178*fIn[29]-1.414213562373095*fIn[7])*dS*vBoundary; 
    out[5] += ((-3.16227766016838*fIn[103])+2.449489742783178*fIn[30]-1.414213562373095*fIn[8])*dS*vBoundary; 
    out[6] += ((-3.16227766016838*fIn[104])+2.449489742783178*fIn[31]-1.414213562373095*fIn[9])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[54]-1.414213562373095*fIn[22])*dS*vBoundary; 
    out[8] += (2.449489742783178*fIn[55]-1.414213562373095*fIn[23])*dS*vBoundary; 
    out[9] += (2.449489742783178*fIn[56]-1.414213562373095*fIn[24])*dS*vBoundary; 
    out[10] += ((-3.16227766016838*fIn[162])+2.449489742783178*fIn[78]-1.414213562373095*fIn[28])*dS*vBoundary; 
    out[11] += (2.449489742783178*fIn[96]-1.414213562373095*fIn[48])*dS*vBoundary; 
    out[12] += (2.449489742783178*fIn[97]-1.414213562373095*fIn[49])*dS*vBoundary; 
    out[13] += (2.449489742783178*fIn[98]-1.414213562373095*fIn[50])*dS*vBoundary; 
    out[14] += (2.449489742783178*fIn[99]-1.414213562373095*fIn[51])*dS*vBoundary; 
    out[15] += (2.449489742783178*fIn[100]-1.414213562373095*fIn[52])*dS*vBoundary; 
    out[16] += (2.449489742783178*fIn[101]-1.414213562373095*fIn[53])*dS*vBoundary; 
    out[17] += (2.449489742783178*fIn[159]-1.414213562373095*fIn[93])*dS*vBoundary; 
    out[18] += (2.449489742783178*fIn[160]-1.414213562373095*fIn[94])*dS*vBoundary; 
    out[19] += (2.449489742783178*fIn[161]-1.414213562373095*fIn[95])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[25]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[57]+2.449489742783178*fIn[10]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (3.16227766016838*fIn[58]+2.449489742783178*fIn[11]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (3.16227766016838*fIn[59]+2.449489742783178*fIn[12]+1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += (3.16227766016838*fIn[102]+2.449489742783178*fIn[29]+1.414213562373095*fIn[7])*dS*vBoundary; 
    out[5] += (3.16227766016838*fIn[103]+2.449489742783178*fIn[30]+1.414213562373095*fIn[8])*dS*vBoundary; 
    out[6] += (3.16227766016838*fIn[104]+2.449489742783178*fIn[31]+1.414213562373095*fIn[9])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[54]+1.414213562373095*fIn[22])*dS*vBoundary; 
    out[8] += (2.449489742783178*fIn[55]+1.414213562373095*fIn[23])*dS*vBoundary; 
    out[9] += (2.449489742783178*fIn[56]+1.414213562373095*fIn[24])*dS*vBoundary; 
    out[10] += (3.16227766016838*fIn[162]+2.449489742783178*fIn[78]+1.414213562373095*fIn[28])*dS*vBoundary; 
    out[11] += (2.449489742783178*fIn[96]+1.414213562373095*fIn[48])*dS*vBoundary; 
    out[12] += (2.449489742783178*fIn[97]+1.414213562373095*fIn[49])*dS*vBoundary; 
    out[13] += (2.449489742783178*fIn[98]+1.414213562373095*fIn[50])*dS*vBoundary; 
    out[14] += (2.449489742783178*fIn[99]+1.414213562373095*fIn[51])*dS*vBoundary; 
    out[15] += (2.449489742783178*fIn[100]+1.414213562373095*fIn[52])*dS*vBoundary; 
    out[16] += (2.449489742783178*fIn[101]+1.414213562373095*fIn[53])*dS*vBoundary; 
    out[17] += (2.449489742783178*fIn[159]+1.414213562373095*fIn[93])*dS*vBoundary; 
    out[18] += (2.449489742783178*fIn[160]+1.414213562373095*fIn[94])*dS*vBoundary; 
    out[19] += (2.449489742783178*fIn[161]+1.414213562373095*fIn[95])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[64]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[5]; 
 
  if (atLower) {
 
    out[8] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
    out[9] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[14]-1.414213562373095*fIn[2])*dS; 
    out[11] += (2.449489742783178*fIn[15]-1.414213562373095*fIn[3])*dS; 
    out[12] += (2.449489742783178*fIn[26]-1.414213562373095*fIn[7])*dS; 
    out[13] += (2.449489742783178*fIn[27]-1.414213562373095*fIn[8])*dS; 
    out[14] += (2.449489742783178*fIn[28]-1.414213562373095*fIn[9])*dS; 
    out[15] += (2.449489742783178*fIn[43]-1.414213562373095*fIn[22])*dS; 
 
  } else {
 
    out[8] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS; 
    out[9] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[14]+1.414213562373095*fIn[2])*dS; 
    out[11] += (2.449489742783178*fIn[15]+1.414213562373095*fIn[3])*dS; 
    out[12] += (2.449489742783178*fIn[26]+1.414213562373095*fIn[7])*dS; 
    out[13] += (2.449489742783178*fIn[27]+1.414213562373095*fIn[8])*dS; 
    out[14] += (2.449489742783178*fIn[28]+1.414213562373095*fIn[9])*dS; 
    out[15] += (2.449489742783178*fIn[43]+1.414213562373095*fIn[22])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[256]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[5]; 
 
  if (atLower) {
 
    out[20] += ((-3.16227766016838*fIn[26])+2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
    out[21] += ((-3.16227766016838*fIn[64])+2.449489742783178*fIn[13]-1.414213562373095*fIn[1])*dS; 
    out[22] += ((-3.16227766016838*fIn[65])+2.449489742783178*fIn[14]-1.414213562373095*fIn[2])*dS; 
    out[23] += ((-3.16227766016838*fIn[66])+2.449489742783178*fIn[15]-1.414213562373095*fIn[3])*dS; 
    out[24] += ((-3.16227766016838*fIn[117])+2.449489742783178*fIn[32]-1.414213562373095*fIn[7])*dS; 
    out[25] += ((-3.16227766016838*fIn[118])+2.449489742783178*fIn[33]-1.414213562373095*fIn[8])*dS; 
    out[26] += ((-3.16227766016838*fIn[119])+2.449489742783178*fIn[34]-1.414213562373095*fIn[9])*dS; 
    out[27] += (2.449489742783178*fIn[60]-1.414213562373095*fIn[22])*dS; 
    out[28] += (2.449489742783178*fIn[61]-1.414213562373095*fIn[23])*dS; 
    out[29] += (2.449489742783178*fIn[62]-1.414213562373095*fIn[24])*dS; 
    out[30] += ((-3.16227766016838*fIn[175])+2.449489742783178*fIn[79]-1.414213562373095*fIn[28])*dS; 
    out[31] += (2.449489742783178*fIn[105]-1.414213562373095*fIn[48])*dS; 
    out[32] += (2.449489742783178*fIn[106]-1.414213562373095*fIn[49])*dS; 
    out[33] += (2.449489742783178*fIn[107]-1.414213562373095*fIn[50])*dS; 
    out[34] += (2.449489742783178*fIn[108]-1.414213562373095*fIn[51])*dS; 
    out[35] += (2.449489742783178*fIn[109]-1.414213562373095*fIn[52])*dS; 
    out[36] += (2.449489742783178*fIn[110]-1.414213562373095*fIn[53])*dS; 
    out[37] += (2.449489742783178*fIn[163]-1.414213562373095*fIn[93])*dS; 
    out[38] += (2.449489742783178*fIn[164]-1.414213562373095*fIn[94])*dS; 
    out[39] += (2.449489742783178*fIn[165]-1.414213562373095*fIn[95])*dS; 
 
  } else {
 
    out[20] += (3.16227766016838*fIn[26]+2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS; 
    out[21] += (3.16227766016838*fIn[64]+2.449489742783178*fIn[13]+1.414213562373095*fIn[1])*dS; 
    out[22] += (3.16227766016838*fIn[65]+2.449489742783178*fIn[14]+1.414213562373095*fIn[2])*dS; 
    out[23] += (3.16227766016838*fIn[66]+2.449489742783178*fIn[15]+1.414213562373095*fIn[3])*dS; 
    out[24] += (3.16227766016838*fIn[117]+2.449489742783178*fIn[32]+1.414213562373095*fIn[7])*dS; 
    out[25] += (3.16227766016838*fIn[118]+2.449489742783178*fIn[33]+1.414213562373095*fIn[8])*dS; 
    out[26] += (3.16227766016838*fIn[119]+2.449489742783178*fIn[34]+1.414213562373095*fIn[9])*dS; 
    out[27] += (2.449489742783178*fIn[60]+1.414213562373095*fIn[22])*dS; 
    out[28] += (2.449489742783178*fIn[61]+1.414213562373095*fIn[23])*dS; 
    out[29] += (2.449489742783178*fIn[62]+1.414213562373095*fIn[24])*dS; 
    out[30] += (3.16227766016838*fIn[175]+2.449489742783178*fIn[79]+1.414213562373095*fIn[28])*dS; 
    out[31] += (2.449489742783178*fIn[105]+1.414213562373095*fIn[48])*dS; 
    out[32] += (2.449489742783178*fIn[106]+1.414213562373095*fIn[49])*dS; 
    out[33] += (2.449489742783178*fIn[107]+1.414213562373095*fIn[50])*dS; 
    out[34] += (2.449489742783178*fIn[108]+1.414213562373095*fIn[51])*dS; 
    out[35] += (2.449489742783178*fIn[109]+1.414213562373095*fIn[52])*dS; 
    out[36] += (2.449489742783178*fIn[110]+1.414213562373095*fIn[53])*dS; 
    out[37] += (2.449489742783178*fIn[163]+1.414213562373095*fIn[93])*dS; 
    out[38] += (2.449489742783178*fIn[164]+1.414213562373095*fIn[94])*dS; 
    out[39] += (2.449489742783178*fIn[165]+1.414213562373095*fIn[95])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[64]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[5]-0.7071067811865475*fIn[0]*dxv[4])*dS; 
    out[1] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[13]-0.7071067811865475*fIn[1]*dxv[4])*dS; 
    out[2] += (2.449489742783178*fIn[14]-1.414213562373095*fIn[2])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[14]-0.7071067811865475*fIn[2]*dxv[4])*dS; 
    out[3] += (2.449489742783178*fIn[15]-1.414213562373095*fIn[3])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[15]-0.7071067811865475*fIn[3]*dxv[4])*dS; 
    out[4] += (2.449489742783178*fIn[26]-1.414213562373095*fIn[7])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[26]-0.7071067811865475*dxv[4]*fIn[7])*dS; 
    out[5] += (2.449489742783178*fIn[27]-1.414213562373095*fIn[8])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[27]-0.7071067811865475*dxv[4]*fIn[8])*dS; 
    out[6] += (2.449489742783178*fIn[28]-1.414213562373095*fIn[9])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[28]-0.7071067811865475*dxv[4]*fIn[9])*dS; 
    out[7] += (2.449489742783178*fIn[43]-1.414213562373095*fIn[22])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[43]-0.7071067811865475*dxv[4]*fIn[22])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[5])-0.7071067811865475*fIn[0]*dxv[4])*dS; 
    out[1] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[13])-0.7071067811865475*fIn[1]*dxv[4])*dS; 
    out[2] += (2.449489742783178*fIn[14]+1.414213562373095*fIn[2])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[14])-0.7071067811865475*fIn[2]*dxv[4])*dS; 
    out[3] += (2.449489742783178*fIn[15]+1.414213562373095*fIn[3])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[15])-0.7071067811865475*fIn[3]*dxv[4])*dS; 
    out[4] += (2.449489742783178*fIn[26]+1.414213562373095*fIn[7])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[26])-0.7071067811865475*dxv[4]*fIn[7])*dS; 
    out[5] += (2.449489742783178*fIn[27]+1.414213562373095*fIn[8])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[27])-0.7071067811865475*dxv[4]*fIn[8])*dS; 
    out[6] += (2.449489742783178*fIn[28]+1.414213562373095*fIn[9])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[28])-0.7071067811865475*dxv[4]*fIn[9])*dS; 
    out[7] += (2.449489742783178*fIn[43]+1.414213562373095*fIn[22])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[43])-0.7071067811865475*dxv[4]*fIn[22])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[256]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[26])+2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[64])+2.449489742783178*fIn[13]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += ((-3.16227766016838*fIn[65])+2.449489742783178*fIn[14]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += ((-3.16227766016838*fIn[66])+2.449489742783178*fIn[15]-1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += ((-3.16227766016838*fIn[117])+2.449489742783178*fIn[32]-1.414213562373095*fIn[7])*dS*vBoundary; 
    out[5] += ((-3.16227766016838*fIn[118])+2.449489742783178*fIn[33]-1.414213562373095*fIn[8])*dS*vBoundary; 
    out[6] += ((-3.16227766016838*fIn[119])+2.449489742783178*fIn[34]-1.414213562373095*fIn[9])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[60]-1.414213562373095*fIn[22])*dS*vBoundary; 
    out[8] += (2.449489742783178*fIn[61]-1.414213562373095*fIn[23])*dS*vBoundary; 
    out[9] += (2.449489742783178*fIn[62]-1.414213562373095*fIn[24])*dS*vBoundary; 
    out[10] += ((-3.16227766016838*fIn[175])+2.449489742783178*fIn[79]-1.414213562373095*fIn[28])*dS*vBoundary; 
    out[11] += (2.449489742783178*fIn[105]-1.414213562373095*fIn[48])*dS*vBoundary; 
    out[12] += (2.449489742783178*fIn[106]-1.414213562373095*fIn[49])*dS*vBoundary; 
    out[13] += (2.449489742783178*fIn[107]-1.414213562373095*fIn[50])*dS*vBoundary; 
    out[14] += (2.449489742783178*fIn[108]-1.414213562373095*fIn[51])*dS*vBoundary; 
    out[15] += (2.449489742783178*fIn[109]-1.414213562373095*fIn[52])*dS*vBoundary; 
    out[16] += (2.449489742783178*fIn[110]-1.414213562373095*fIn[53])*dS*vBoundary; 
    out[17] += (2.449489742783178*fIn[163]-1.414213562373095*fIn[93])*dS*vBoundary; 
    out[18] += (2.449489742783178*fIn[164]-1.414213562373095*fIn[94])*dS*vBoundary; 
    out[19] += (2.449489742783178*fIn[165]-1.414213562373095*fIn[95])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[26]+2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[64]+2.449489742783178*fIn[13]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (3.16227766016838*fIn[65]+2.449489742783178*fIn[14]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (3.16227766016838*fIn[66]+2.449489742783178*fIn[15]+1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += (3.16227766016838*fIn[117]+2.449489742783178*fIn[32]+1.414213562373095*fIn[7])*dS*vBoundary; 
    out[5] += (3.16227766016838*fIn[118]+2.449489742783178*fIn[33]+1.414213562373095*fIn[8])*dS*vBoundary; 
    out[6] += (3.16227766016838*fIn[119]+2.449489742783178*fIn[34]+1.414213562373095*fIn[9])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[60]+1.414213562373095*fIn[22])*dS*vBoundary; 
    out[8] += (2.449489742783178*fIn[61]+1.414213562373095*fIn[23])*dS*vBoundary; 
    out[9] += (2.449489742783178*fIn[62]+1.414213562373095*fIn[24])*dS*vBoundary; 
    out[10] += (3.16227766016838*fIn[175]+2.449489742783178*fIn[79]+1.414213562373095*fIn[28])*dS*vBoundary; 
    out[11] += (2.449489742783178*fIn[105]+1.414213562373095*fIn[48])*dS*vBoundary; 
    out[12] += (2.449489742783178*fIn[106]+1.414213562373095*fIn[49])*dS*vBoundary; 
    out[13] += (2.449489742783178*fIn[107]+1.414213562373095*fIn[50])*dS*vBoundary; 
    out[14] += (2.449489742783178*fIn[108]+1.414213562373095*fIn[51])*dS*vBoundary; 
    out[15] += (2.449489742783178*fIn[109]+1.414213562373095*fIn[52])*dS*vBoundary; 
    out[16] += (2.449489742783178*fIn[110]+1.414213562373095*fIn[53])*dS*vBoundary; 
    out[17] += (2.449489742783178*fIn[163]+1.414213562373095*fIn[93])*dS*vBoundary; 
    out[18] += (2.449489742783178*fIn[164]+1.414213562373095*fIn[94])*dS*vBoundary; 
    out[19] += (2.449489742783178*fIn[165]+1.414213562373095*fIn[95])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[64]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[16] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[0])*dS; 
    out[17] += (2.449489742783178*fIn[17]-1.414213562373095*fIn[1])*dS; 
    out[18] += (2.449489742783178*fIn[18]-1.414213562373095*fIn[2])*dS; 
    out[19] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[3])*dS; 
    out[20] += (2.449489742783178*fIn[32]-1.414213562373095*fIn[7])*dS; 
    out[21] += (2.449489742783178*fIn[33]-1.414213562373095*fIn[8])*dS; 
    out[22] += (2.449489742783178*fIn[34]-1.414213562373095*fIn[9])*dS; 
    out[23] += (2.449489742783178*fIn[47]-1.414213562373095*fIn[22])*dS; 
 
  } else {
 
    out[16] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[0])*dS; 
    out[17] += (2.449489742783178*fIn[17]+1.414213562373095*fIn[1])*dS; 
    out[18] += (2.449489742783178*fIn[18]+1.414213562373095*fIn[2])*dS; 
    out[19] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[3])*dS; 
    out[20] += (2.449489742783178*fIn[32]+1.414213562373095*fIn[7])*dS; 
    out[21] += (2.449489742783178*fIn[33]+1.414213562373095*fIn[8])*dS; 
    out[22] += (2.449489742783178*fIn[34]+1.414213562373095*fIn[9])*dS; 
    out[23] += (2.449489742783178*fIn[47]+1.414213562373095*fIn[22])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[256]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[40] += ((-3.16227766016838*fIn[27])+2.449489742783178*fIn[6]-1.414213562373095*fIn[0])*dS; 
    out[41] += ((-3.16227766016838*fIn[73])+2.449489742783178*fIn[17]-1.414213562373095*fIn[1])*dS; 
    out[42] += ((-3.16227766016838*fIn[74])+2.449489742783178*fIn[18]-1.414213562373095*fIn[2])*dS; 
    out[43] += ((-3.16227766016838*fIn[75])+2.449489742783178*fIn[19]-1.414213562373095*fIn[3])*dS; 
    out[44] += ((-3.16227766016838*fIn[143])+2.449489742783178*fIn[38]-1.414213562373095*fIn[7])*dS; 
    out[45] += ((-3.16227766016838*fIn[144])+2.449489742783178*fIn[39]-1.414213562373095*fIn[8])*dS; 
    out[46] += ((-3.16227766016838*fIn[145])+2.449489742783178*fIn[40]-1.414213562373095*fIn[9])*dS; 
    out[47] += (2.449489742783178*fIn[68]-1.414213562373095*fIn[22])*dS; 
    out[48] += (2.449489742783178*fIn[69]-1.414213562373095*fIn[23])*dS; 
    out[49] += (2.449489742783178*fIn[70]-1.414213562373095*fIn[24])*dS; 
    out[50] += ((-3.16227766016838*fIn[209])+2.449489742783178*fIn[83]-1.414213562373095*fIn[28])*dS; 
    out[51] += (2.449489742783178*fIn[123]-1.414213562373095*fIn[48])*dS; 
    out[52] += (2.449489742783178*fIn[124]-1.414213562373095*fIn[49])*dS; 
    out[53] += (2.449489742783178*fIn[125]-1.414213562373095*fIn[50])*dS; 
    out[54] += (2.449489742783178*fIn[126]-1.414213562373095*fIn[51])*dS; 
    out[55] += (2.449489742783178*fIn[127]-1.414213562373095*fIn[52])*dS; 
    out[56] += (2.449489742783178*fIn[128]-1.414213562373095*fIn[53])*dS; 
    out[57] += (2.449489742783178*fIn[179]-1.414213562373095*fIn[93])*dS; 
    out[58] += (2.449489742783178*fIn[180]-1.414213562373095*fIn[94])*dS; 
    out[59] += (2.449489742783178*fIn[181]-1.414213562373095*fIn[95])*dS; 
 
  } else {
 
    out[40] += (3.16227766016838*fIn[27]+2.449489742783178*fIn[6]+1.414213562373095*fIn[0])*dS; 
    out[41] += (3.16227766016838*fIn[73]+2.449489742783178*fIn[17]+1.414213562373095*fIn[1])*dS; 
    out[42] += (3.16227766016838*fIn[74]+2.449489742783178*fIn[18]+1.414213562373095*fIn[2])*dS; 
    out[43] += (3.16227766016838*fIn[75]+2.449489742783178*fIn[19]+1.414213562373095*fIn[3])*dS; 
    out[44] += (3.16227766016838*fIn[143]+2.449489742783178*fIn[38]+1.414213562373095*fIn[7])*dS; 
    out[45] += (3.16227766016838*fIn[144]+2.449489742783178*fIn[39]+1.414213562373095*fIn[8])*dS; 
    out[46] += (3.16227766016838*fIn[145]+2.449489742783178*fIn[40]+1.414213562373095*fIn[9])*dS; 
    out[47] += (2.449489742783178*fIn[68]+1.414213562373095*fIn[22])*dS; 
    out[48] += (2.449489742783178*fIn[69]+1.414213562373095*fIn[23])*dS; 
    out[49] += (2.449489742783178*fIn[70]+1.414213562373095*fIn[24])*dS; 
    out[50] += (3.16227766016838*fIn[209]+2.449489742783178*fIn[83]+1.414213562373095*fIn[28])*dS; 
    out[51] += (2.449489742783178*fIn[123]+1.414213562373095*fIn[48])*dS; 
    out[52] += (2.449489742783178*fIn[124]+1.414213562373095*fIn[49])*dS; 
    out[53] += (2.449489742783178*fIn[125]+1.414213562373095*fIn[50])*dS; 
    out[54] += (2.449489742783178*fIn[126]+1.414213562373095*fIn[51])*dS; 
    out[55] += (2.449489742783178*fIn[127]+1.414213562373095*fIn[52])*dS; 
    out[56] += (2.449489742783178*fIn[128]+1.414213562373095*fIn[53])*dS; 
    out[57] += (2.449489742783178*fIn[179]+1.414213562373095*fIn[93])*dS; 
    out[58] += (2.449489742783178*fIn[180]+1.414213562373095*fIn[94])*dS; 
    out[59] += (2.449489742783178*fIn[181]+1.414213562373095*fIn[95])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[64]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[5]*fIn[6]-0.7071067811865475*fIn[0]*dxv[5])*dS; 
    out[1] += (2.449489742783178*fIn[17]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[5]*fIn[17]-0.7071067811865475*fIn[1]*dxv[5])*dS; 
    out[2] += (2.449489742783178*fIn[18]-1.414213562373095*fIn[2])*dS*vBoundary+(1.224744871391589*dxv[5]*fIn[18]-0.7071067811865475*fIn[2]*dxv[5])*dS; 
    out[3] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[3])*dS*vBoundary+(1.224744871391589*dxv[5]*fIn[19]-0.7071067811865475*fIn[3]*dxv[5])*dS; 
    out[4] += (2.449489742783178*fIn[32]-1.414213562373095*fIn[7])*dS*vBoundary+(1.224744871391589*dxv[5]*fIn[32]-0.7071067811865475*dxv[5]*fIn[7])*dS; 
    out[5] += (2.449489742783178*fIn[33]-1.414213562373095*fIn[8])*dS*vBoundary+(1.224744871391589*dxv[5]*fIn[33]-0.7071067811865475*dxv[5]*fIn[8])*dS; 
    out[6] += (2.449489742783178*fIn[34]-1.414213562373095*fIn[9])*dS*vBoundary+(1.224744871391589*dxv[5]*fIn[34]-0.7071067811865475*dxv[5]*fIn[9])*dS; 
    out[7] += (2.449489742783178*fIn[47]-1.414213562373095*fIn[22])*dS*vBoundary+(1.224744871391589*dxv[5]*fIn[47]-0.7071067811865475*dxv[5]*fIn[22])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[5]*fIn[6])-0.7071067811865475*fIn[0]*dxv[5])*dS; 
    out[1] += (2.449489742783178*fIn[17]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[5]*fIn[17])-0.7071067811865475*fIn[1]*dxv[5])*dS; 
    out[2] += (2.449489742783178*fIn[18]+1.414213562373095*fIn[2])*dS*vBoundary+((-1.224744871391589*dxv[5]*fIn[18])-0.7071067811865475*fIn[2]*dxv[5])*dS; 
    out[3] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[3])*dS*vBoundary+((-1.224744871391589*dxv[5]*fIn[19])-0.7071067811865475*fIn[3]*dxv[5])*dS; 
    out[4] += (2.449489742783178*fIn[32]+1.414213562373095*fIn[7])*dS*vBoundary+((-1.224744871391589*dxv[5]*fIn[32])-0.7071067811865475*dxv[5]*fIn[7])*dS; 
    out[5] += (2.449489742783178*fIn[33]+1.414213562373095*fIn[8])*dS*vBoundary+((-1.224744871391589*dxv[5]*fIn[33])-0.7071067811865475*dxv[5]*fIn[8])*dS; 
    out[6] += (2.449489742783178*fIn[34]+1.414213562373095*fIn[9])*dS*vBoundary+((-1.224744871391589*dxv[5]*fIn[34])-0.7071067811865475*dxv[5]*fIn[9])*dS; 
    out[7] += (2.449489742783178*fIn[47]+1.414213562373095*fIn[22])*dS*vBoundary+((-1.224744871391589*dxv[5]*fIn[47])-0.7071067811865475*dxv[5]*fIn[22])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vSer_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[256]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[27])+2.449489742783178*fIn[6]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[73])+2.449489742783178*fIn[17]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += ((-3.16227766016838*fIn[74])+2.449489742783178*fIn[18]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += ((-3.16227766016838*fIn[75])+2.449489742783178*fIn[19]-1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += ((-3.16227766016838*fIn[143])+2.449489742783178*fIn[38]-1.414213562373095*fIn[7])*dS*vBoundary; 
    out[5] += ((-3.16227766016838*fIn[144])+2.449489742783178*fIn[39]-1.414213562373095*fIn[8])*dS*vBoundary; 
    out[6] += ((-3.16227766016838*fIn[145])+2.449489742783178*fIn[40]-1.414213562373095*fIn[9])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[68]-1.414213562373095*fIn[22])*dS*vBoundary; 
    out[8] += (2.449489742783178*fIn[69]-1.414213562373095*fIn[23])*dS*vBoundary; 
    out[9] += (2.449489742783178*fIn[70]-1.414213562373095*fIn[24])*dS*vBoundary; 
    out[10] += ((-3.16227766016838*fIn[209])+2.449489742783178*fIn[83]-1.414213562373095*fIn[28])*dS*vBoundary; 
    out[11] += (2.449489742783178*fIn[123]-1.414213562373095*fIn[48])*dS*vBoundary; 
    out[12] += (2.449489742783178*fIn[124]-1.414213562373095*fIn[49])*dS*vBoundary; 
    out[13] += (2.449489742783178*fIn[125]-1.414213562373095*fIn[50])*dS*vBoundary; 
    out[14] += (2.449489742783178*fIn[126]-1.414213562373095*fIn[51])*dS*vBoundary; 
    out[15] += (2.449489742783178*fIn[127]-1.414213562373095*fIn[52])*dS*vBoundary; 
    out[16] += (2.449489742783178*fIn[128]-1.414213562373095*fIn[53])*dS*vBoundary; 
    out[17] += (2.449489742783178*fIn[179]-1.414213562373095*fIn[93])*dS*vBoundary; 
    out[18] += (2.449489742783178*fIn[180]-1.414213562373095*fIn[94])*dS*vBoundary; 
    out[19] += (2.449489742783178*fIn[181]-1.414213562373095*fIn[95])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[27]+2.449489742783178*fIn[6]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[73]+2.449489742783178*fIn[17]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (3.16227766016838*fIn[74]+2.449489742783178*fIn[18]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (3.16227766016838*fIn[75]+2.449489742783178*fIn[19]+1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += (3.16227766016838*fIn[143]+2.449489742783178*fIn[38]+1.414213562373095*fIn[7])*dS*vBoundary; 
    out[5] += (3.16227766016838*fIn[144]+2.449489742783178*fIn[39]+1.414213562373095*fIn[8])*dS*vBoundary; 
    out[6] += (3.16227766016838*fIn[145]+2.449489742783178*fIn[40]+1.414213562373095*fIn[9])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[68]+1.414213562373095*fIn[22])*dS*vBoundary; 
    out[8] += (2.449489742783178*fIn[69]+1.414213562373095*fIn[23])*dS*vBoundary; 
    out[9] += (2.449489742783178*fIn[70]+1.414213562373095*fIn[24])*dS*vBoundary; 
    out[10] += (3.16227766016838*fIn[209]+2.449489742783178*fIn[83]+1.414213562373095*fIn[28])*dS*vBoundary; 
    out[11] += (2.449489742783178*fIn[123]+1.414213562373095*fIn[48])*dS*vBoundary; 
    out[12] += (2.449489742783178*fIn[124]+1.414213562373095*fIn[49])*dS*vBoundary; 
    out[13] += (2.449489742783178*fIn[125]+1.414213562373095*fIn[50])*dS*vBoundary; 
    out[14] += (2.449489742783178*fIn[126]+1.414213562373095*fIn[51])*dS*vBoundary; 
    out[15] += (2.449489742783178*fIn[127]+1.414213562373095*fIn[52])*dS*vBoundary; 
    out[16] += (2.449489742783178*fIn[128]+1.414213562373095*fIn[53])*dS*vBoundary; 
    out[17] += (2.449489742783178*fIn[179]+1.414213562373095*fIn[93])*dS*vBoundary; 
    out[18] += (2.449489742783178*fIn[180]+1.414213562373095*fIn[94])*dS*vBoundary; 
    out[19] += (2.449489742783178*fIn[181]+1.414213562373095*fIn[95])*dS*vBoundary; 
 
  }
 
} 
 
