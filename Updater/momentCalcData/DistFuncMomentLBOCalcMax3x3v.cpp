#include <DistFuncMomentCalcModDecl.h> 
void VmM0Star3x3vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[6]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[4]*dxvl[5]*(wr[3]-wl[3]); 
 
  out[0] += ((-0.8164965809277261*fr[4])+0.8164965809277261*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += (0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += (0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
  out[3] += (0.7071067811865475*fr[3]+0.7071067811865475*fl[3])*dS; 
 
} 
 
void VmM0Star3x3vMax_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[6]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[3]*dxvl[5]*(wr[4]-wl[4]); 
 
  out[0] += ((-0.8164965809277261*fr[5])+0.8164965809277261*fl[5]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += (0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += (0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
  out[3] += (0.7071067811865475*fr[3]+0.7071067811865475*fl[3])*dS; 
 
} 
 
void VmM0Star3x3vMax_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[6]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[3]*dxvl[4]*(wr[5]-wl[5]); 
 
  out[0] += ((-0.8164965809277261*fr[6])+0.8164965809277261*fl[6]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += (0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += (0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
  out[3] += (0.7071067811865475*fr[3]+0.7071067811865475*fl[3])*dS; 
 
} 
 
void VmM1iM2Star3x3vMax(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2) 
{ 
  // w[6]:    Cell-center coordinates. 
  // dxv[6]:  Cell length in each direciton. 
  // intFac:  for gyrokinetics (not used in Vm). 
  // m_:      mass (not used in Vm). 
  // Bmag[4]: Magnetic field magnitude (not used in Vm). 
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
  double tempM0[4]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 
  tempM0[2] = 2.828427124746191*f[2]*volFact; 
  tempM0[3] = 2.828427124746191*f[3]*volFact; 

  outM1i[0] += tempM0[0]*w[3]; 
  outM1i[1] += tempM0[1]*w[3]; 
  outM1i[2] += tempM0[2]*w[3]; 
  outM1i[3] += w[3]*tempM0[3]; 
  outM1i[4] += tempM0[0]*w[4]; 
  outM1i[5] += tempM0[1]*w[4]; 
  outM1i[6] += tempM0[2]*w[4]; 
  outM1i[7] += tempM0[3]*w[4]; 
  outM1i[8] += tempM0[0]*w[5]; 
  outM1i[9] += tempM0[1]*w[5]; 
  outM1i[10] += tempM0[2]*w[5]; 
  outM1i[11] += tempM0[3]*w[5]; 

  outM2[0] += (0.8164965809277261*dxv[5]*w[5]*f[6]+0.8164965809277261*dxv[4]*w[4]*f[5]+0.8164965809277261*dxv[3]*w[3]*f[4])*volFact+tempM0[0]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[1] += tempM0[1]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[2] += tempM0[2]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[3] += tempM0[3]*(wvSq[2]+wvSq[1]+wvSq[0]); 
 
} 
void VmBoundaryIntegral3x3vMax_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[7]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[4]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[1] += -1.414213562373095*fIn[1]*dS; 
    out[2] += -1.414213562373095*fIn[2]*dS; 
    out[3] += -1.414213562373095*fIn[3]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS; 
    out[2] += 1.414213562373095*fIn[2]*dS; 
    out[3] += 1.414213562373095*fIn[3]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[28]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[4]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[25])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[11]-1.414213562373095*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[3])*dS; 
    out[4] += -1.414213562373095*fIn[7]*dS; 
    out[5] += -1.414213562373095*fIn[8]*dS; 
    out[6] += -1.414213562373095*fIn[9]*dS; 
    out[7] += -1.414213562373095*fIn[22]*dS; 
    out[8] += -1.414213562373095*fIn[23]*dS; 
    out[9] += -1.414213562373095*fIn[24]*dS; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[25]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[11]+1.414213562373095*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[3])*dS; 
    out[4] += 1.414213562373095*fIn[7]*dS; 
    out[5] += 1.414213562373095*fIn[8]*dS; 
    out[6] += 1.414213562373095*fIn[9]*dS; 
    out[7] += 1.414213562373095*fIn[22]*dS; 
    out[8] += 1.414213562373095*fIn[23]*dS; 
    out[9] += 1.414213562373095*fIn[24]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[7]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[4]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[4]-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (-1.414213562373095*fIn[1]*dS*vBoundary)-0.7071067811865475*fIn[1]*dxv[3]*dS; 
    out[2] += (-1.414213562373095*fIn[2]*dS*vBoundary)-0.7071067811865475*fIn[2]*dxv[3]*dS; 
    out[3] += (-1.414213562373095*fIn[3]*dS*vBoundary)-0.7071067811865475*dxv[3]*fIn[3]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[4])-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS*vBoundary-0.7071067811865475*fIn[1]*dxv[3]*dS; 
    out[2] += 1.414213562373095*fIn[2]*dS*vBoundary-0.7071067811865475*fIn[2]*dxv[3]*dS; 
    out[3] += 1.414213562373095*fIn[3]*dS*vBoundary-0.7071067811865475*dxv[3]*fIn[3]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[28]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[4]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[25])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[11]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += -1.414213562373095*fIn[7]*dS*vBoundary; 
    out[5] += -1.414213562373095*fIn[8]*dS*vBoundary; 
    out[6] += -1.414213562373095*fIn[9]*dS*vBoundary; 
    out[7] += -1.414213562373095*fIn[22]*dS*vBoundary; 
    out[8] += -1.414213562373095*fIn[23]*dS*vBoundary; 
    out[9] += -1.414213562373095*fIn[24]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[25]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[11]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += 1.414213562373095*fIn[7]*dS*vBoundary; 
    out[5] += 1.414213562373095*fIn[8]*dS*vBoundary; 
    out[6] += 1.414213562373095*fIn[9]*dS*vBoundary; 
    out[7] += 1.414213562373095*fIn[22]*dS*vBoundary; 
    out[8] += 1.414213562373095*fIn[23]*dS*vBoundary; 
    out[9] += 1.414213562373095*fIn[24]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[7]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[5]; 
 
  if (atLower) {
 
    out[4] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
    out[5] += -1.414213562373095*fIn[1]*dS; 
    out[6] += -1.414213562373095*fIn[2]*dS; 
    out[7] += -1.414213562373095*fIn[3]*dS; 
 
  } else {
 
    out[4] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS; 
    out[5] += 1.414213562373095*fIn[1]*dS; 
    out[6] += 1.414213562373095*fIn[2]*dS; 
    out[7] += 1.414213562373095*fIn[3]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[28]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[5]; 
 
  if (atLower) {
 
    out[10] += ((-3.16227766016838*fIn[26])+2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
    out[11] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[1])*dS; 
    out[12] += (2.449489742783178*fIn[14]-1.414213562373095*fIn[2])*dS; 
    out[13] += (2.449489742783178*fIn[15]-1.414213562373095*fIn[3])*dS; 
    out[14] += -1.414213562373095*fIn[7]*dS; 
    out[15] += -1.414213562373095*fIn[8]*dS; 
    out[16] += -1.414213562373095*fIn[9]*dS; 
    out[17] += -1.414213562373095*fIn[22]*dS; 
    out[18] += -1.414213562373095*fIn[23]*dS; 
    out[19] += -1.414213562373095*fIn[24]*dS; 
 
  } else {
 
    out[10] += (3.16227766016838*fIn[26]+2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS; 
    out[11] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[1])*dS; 
    out[12] += (2.449489742783178*fIn[14]+1.414213562373095*fIn[2])*dS; 
    out[13] += (2.449489742783178*fIn[15]+1.414213562373095*fIn[3])*dS; 
    out[14] += 1.414213562373095*fIn[7]*dS; 
    out[15] += 1.414213562373095*fIn[8]*dS; 
    out[16] += 1.414213562373095*fIn[9]*dS; 
    out[17] += 1.414213562373095*fIn[22]*dS; 
    out[18] += 1.414213562373095*fIn[23]*dS; 
    out[19] += 1.414213562373095*fIn[24]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[7]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[5]-0.7071067811865475*fIn[0]*dxv[4])*dS; 
    out[1] += (-1.414213562373095*fIn[1]*dS*vBoundary)-0.7071067811865475*fIn[1]*dxv[4]*dS; 
    out[2] += (-1.414213562373095*fIn[2]*dS*vBoundary)-0.7071067811865475*fIn[2]*dxv[4]*dS; 
    out[3] += (-1.414213562373095*fIn[3]*dS*vBoundary)-0.7071067811865475*fIn[3]*dxv[4]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[5])-0.7071067811865475*fIn[0]*dxv[4])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS*vBoundary-0.7071067811865475*fIn[1]*dxv[4]*dS; 
    out[2] += 1.414213562373095*fIn[2]*dS*vBoundary-0.7071067811865475*fIn[2]*dxv[4]*dS; 
    out[3] += 1.414213562373095*fIn[3]*dS*vBoundary-0.7071067811865475*fIn[3]*dxv[4]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[28]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[5]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[26])+2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[14]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[15]-1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += -1.414213562373095*fIn[7]*dS*vBoundary; 
    out[5] += -1.414213562373095*fIn[8]*dS*vBoundary; 
    out[6] += -1.414213562373095*fIn[9]*dS*vBoundary; 
    out[7] += -1.414213562373095*fIn[22]*dS*vBoundary; 
    out[8] += -1.414213562373095*fIn[23]*dS*vBoundary; 
    out[9] += -1.414213562373095*fIn[24]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[26]+2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[14]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[15]+1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += 1.414213562373095*fIn[7]*dS*vBoundary; 
    out[5] += 1.414213562373095*fIn[8]*dS*vBoundary; 
    out[6] += 1.414213562373095*fIn[9]*dS*vBoundary; 
    out[7] += 1.414213562373095*fIn[22]*dS*vBoundary; 
    out[8] += 1.414213562373095*fIn[23]*dS*vBoundary; 
    out[9] += 1.414213562373095*fIn[24]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[7]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[8] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[0])*dS; 
    out[9] += -1.414213562373095*fIn[1]*dS; 
    out[10] += -1.414213562373095*fIn[2]*dS; 
    out[11] += -1.414213562373095*fIn[3]*dS; 
 
  } else {
 
    out[8] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[0])*dS; 
    out[9] += 1.414213562373095*fIn[1]*dS; 
    out[10] += 1.414213562373095*fIn[2]*dS; 
    out[11] += 1.414213562373095*fIn[3]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[28]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[20] += ((-3.16227766016838*fIn[27])+2.449489742783178*fIn[6]-1.414213562373095*fIn[0])*dS; 
    out[21] += (2.449489742783178*fIn[17]-1.414213562373095*fIn[1])*dS; 
    out[22] += (2.449489742783178*fIn[18]-1.414213562373095*fIn[2])*dS; 
    out[23] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[3])*dS; 
    out[24] += -1.414213562373095*fIn[7]*dS; 
    out[25] += -1.414213562373095*fIn[8]*dS; 
    out[26] += -1.414213562373095*fIn[9]*dS; 
    out[27] += -1.414213562373095*fIn[22]*dS; 
    out[28] += -1.414213562373095*fIn[23]*dS; 
    out[29] += -1.414213562373095*fIn[24]*dS; 
 
  } else {
 
    out[20] += (3.16227766016838*fIn[27]+2.449489742783178*fIn[6]+1.414213562373095*fIn[0])*dS; 
    out[21] += (2.449489742783178*fIn[17]+1.414213562373095*fIn[1])*dS; 
    out[22] += (2.449489742783178*fIn[18]+1.414213562373095*fIn[2])*dS; 
    out[23] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[3])*dS; 
    out[24] += 1.414213562373095*fIn[7]*dS; 
    out[25] += 1.414213562373095*fIn[8]*dS; 
    out[26] += 1.414213562373095*fIn[9]*dS; 
    out[27] += 1.414213562373095*fIn[22]*dS; 
    out[28] += 1.414213562373095*fIn[23]*dS; 
    out[29] += 1.414213562373095*fIn[24]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[7]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[5]*fIn[6]-0.7071067811865475*fIn[0]*dxv[5])*dS; 
    out[1] += (-1.414213562373095*fIn[1]*dS*vBoundary)-0.7071067811865475*fIn[1]*dxv[5]*dS; 
    out[2] += (-1.414213562373095*fIn[2]*dS*vBoundary)-0.7071067811865475*fIn[2]*dxv[5]*dS; 
    out[3] += (-1.414213562373095*fIn[3]*dS*vBoundary)-0.7071067811865475*fIn[3]*dxv[5]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[5]*fIn[6])-0.7071067811865475*fIn[0]*dxv[5])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS*vBoundary-0.7071067811865475*fIn[1]*dxv[5]*dS; 
    out[2] += 1.414213562373095*fIn[2]*dS*vBoundary-0.7071067811865475*fIn[2]*dxv[5]*dS; 
    out[3] += 1.414213562373095*fIn[3]*dS*vBoundary-0.7071067811865475*fIn[3]*dxv[5]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral3x3vMax_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[6]:    cell length in each direciton. 
  // fIn[28]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[27])+2.449489742783178*fIn[6]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[17]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[18]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += -1.414213562373095*fIn[7]*dS*vBoundary; 
    out[5] += -1.414213562373095*fIn[8]*dS*vBoundary; 
    out[6] += -1.414213562373095*fIn[9]*dS*vBoundary; 
    out[7] += -1.414213562373095*fIn[22]*dS*vBoundary; 
    out[8] += -1.414213562373095*fIn[23]*dS*vBoundary; 
    out[9] += -1.414213562373095*fIn[24]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[27]+2.449489742783178*fIn[6]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[17]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[18]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[3])*dS*vBoundary; 
    out[4] += 1.414213562373095*fIn[7]*dS*vBoundary; 
    out[5] += 1.414213562373095*fIn[8]*dS*vBoundary; 
    out[6] += 1.414213562373095*fIn[9]*dS*vBoundary; 
    out[7] += 1.414213562373095*fIn[22]*dS*vBoundary; 
    out[8] += 1.414213562373095*fIn[23]*dS*vBoundary; 
    out[9] += 1.414213562373095*fIn[24]*dS*vBoundary; 
 
  }
 
} 
 
