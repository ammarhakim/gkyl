
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
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
  outM1i[3] += tempM0[3]*w[3]; 
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
 
void VmSelfPrimMoments3x3vMax_P1(binOpData_t* data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1:       moments of the distribution function. 
  // m0S,m1S,m1S: star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[4]; 
  double m1r[12]; 
  double m0Sr[4]; 
  double m1Sr[12]; 
  double m2Sr[4]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = m1[4]; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    m1r[8] = m1[8]; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m0Sr[2] = 0.0; 
    m0Sr[3] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = 0.0; 
    m1Sr[3] = 0.0; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = 0.0; 
    m1Sr[6] = 0.0; 
    m1Sr[7] = 0.0; 
    m1Sr[8] = m1S[8]; 
    m1Sr[9] = 0.0; 
    m1Sr[10] = 0.0; 
    m1Sr[11] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
    m2Sr[2] = 0.0; 
    m2Sr[3] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m1r[4] = m1[4]; 
    m1r[5] = m1[5]; 
    m1r[6] = m1[6]; 
    m1r[7] = m1[7]; 
    m1r[8] = m1[8]; 
    m1r[9] = m1[9]; 
    m1r[10] = m1[10]; 
    m1r[11] = m1[11]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m0Sr[2] = m0S[2]; 
    m0Sr[3] = m0S[3]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = m1S[5]; 
    m1Sr[6] = m1S[6]; 
    m1Sr[7] = m1S[7]; 
    m1Sr[8] = m1S[8]; 
    m1Sr[9] = m1S[9]; 
    m1Sr[10] = m1S[10]; 
    m1Sr[11] = m1S[11]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
    m2Sr[2] = m2S[2]; 
    m2Sr[3] = m2S[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(16,16); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(0,1) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(0,2) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(0,3) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(1,0) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(1,1) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(2,0) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(2,2) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(3,0) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(3,3) = 0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,12) = -0.3535533905932737*cM[0]; 
  data->AEM_S(0,13) = -0.3535533905932737*cM[1]; 
  data->AEM_S(0,14) = -0.3535533905932737*cM[2]; 
  data->AEM_S(0,15) = -0.3535533905932737*cM[3]; 
  data->AEM_S(1,12) = -0.3535533905932737*cM[1]; 
  data->AEM_S(1,13) = -0.3535533905932737*cM[0]; 
  data->AEM_S(2,12) = -0.3535533905932737*cM[2]; 
  data->AEM_S(2,14) = -0.3535533905932737*cM[0]; 
  data->AEM_S(3,12) = -0.3535533905932737*cM[3]; 
  data->AEM_S(3,15) = -0.3535533905932737*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(12,0) = 0.3535533905932737*m1Sr[0]; 
  data->AEM_S(12,1) = 0.3535533905932737*m1Sr[1]; 
  data->AEM_S(12,2) = 0.3535533905932737*m1Sr[2]; 
  data->AEM_S(12,3) = 0.3535533905932737*m1Sr[3]; 
  data->AEM_S(13,0) = 0.3535533905932737*m1Sr[1]; 
  data->AEM_S(13,1) = 0.3535533905932737*m1Sr[0]; 
  data->AEM_S(14,0) = 0.3535533905932737*m1Sr[2]; 
  data->AEM_S(14,2) = 0.3535533905932737*m1Sr[0]; 
  data->AEM_S(15,0) = 0.3535533905932737*m1Sr[3]; 
  data->AEM_S(15,3) = 0.3535533905932737*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(4,4) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(4,5) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(4,6) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(4,7) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(5,4) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(5,5) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(6,4) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(6,6) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(7,4) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(7,7) = 0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(4,12) = -0.3535533905932737*cM[4]; 
  data->AEM_S(4,13) = -0.3535533905932737*cM[5]; 
  data->AEM_S(4,14) = -0.3535533905932737*cM[6]; 
  data->AEM_S(4,15) = -0.3535533905932737*cM[7]; 
  data->AEM_S(5,12) = -0.3535533905932737*cM[5]; 
  data->AEM_S(5,13) = -0.3535533905932737*cM[4]; 
  data->AEM_S(6,12) = -0.3535533905932737*cM[6]; 
  data->AEM_S(6,14) = -0.3535533905932737*cM[4]; 
  data->AEM_S(7,12) = -0.3535533905932737*cM[7]; 
  data->AEM_S(7,15) = -0.3535533905932737*cM[4]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(12,4) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(12,5) = 0.3535533905932737*m1Sr[5]; 
  data->AEM_S(12,6) = 0.3535533905932737*m1Sr[6]; 
  data->AEM_S(12,7) = 0.3535533905932737*m1Sr[7]; 
  data->AEM_S(13,4) = 0.3535533905932737*m1Sr[5]; 
  data->AEM_S(13,5) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(14,4) = 0.3535533905932737*m1Sr[6]; 
  data->AEM_S(14,6) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(15,4) = 0.3535533905932737*m1Sr[7]; 
  data->AEM_S(15,7) = 0.3535533905932737*m1Sr[4]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(8,8) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(8,9) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(8,10) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(8,11) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(9,8) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(9,9) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(10,8) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(10,10) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(11,8) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(11,11) = 0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(8,12) = -0.3535533905932737*cM[8]; 
  data->AEM_S(8,13) = -0.3535533905932737*cM[9]; 
  data->AEM_S(8,14) = -0.3535533905932737*cM[10]; 
  data->AEM_S(8,15) = -0.3535533905932737*cM[11]; 
  data->AEM_S(9,12) = -0.3535533905932737*cM[9]; 
  data->AEM_S(9,13) = -0.3535533905932737*cM[8]; 
  data->AEM_S(10,12) = -0.3535533905932737*cM[10]; 
  data->AEM_S(10,14) = -0.3535533905932737*cM[8]; 
  data->AEM_S(11,12) = -0.3535533905932737*cM[11]; 
  data->AEM_S(11,15) = -0.3535533905932737*cM[8]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(12,8) = 0.3535533905932737*m1Sr[8]; 
  data->AEM_S(12,9) = 0.3535533905932737*m1Sr[9]; 
  data->AEM_S(12,10) = 0.3535533905932737*m1Sr[10]; 
  data->AEM_S(12,11) = 0.3535533905932737*m1Sr[11]; 
  data->AEM_S(13,8) = 0.3535533905932737*m1Sr[9]; 
  data->AEM_S(13,9) = 0.3535533905932737*m1Sr[8]; 
  data->AEM_S(14,8) = 0.3535533905932737*m1Sr[10]; 
  data->AEM_S(14,10) = 0.3535533905932737*m1Sr[8]; 
  data->AEM_S(15,8) = 0.3535533905932737*m1Sr[11]; 
  data->AEM_S(15,11) = 0.3535533905932737*m1Sr[8]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(12,12) = 0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(12,13) = 0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(12,14) = 0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(12,15) = 0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(13,12) = 0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(13,13) = 0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(14,12) = 0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(14,14) = 0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(15,12) = 0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(15,15) = 0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<4,8>(0,4).setZero(); 
  data->AEM_S.block<8,4>(4,0).setZero(); 
  data->AEM_S.block<4,4>(4,8).setZero(); 
  data->AEM_S.block<4,4>(8,4).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m2Sr[0],m2Sr[1],m2Sr[2],m2Sr[3]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,12,1) = data->u_S.segment<12>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = data->u_S.segment<4>(12); 
 
} 
 
void VmSelfPrimMoments3x3vMax_P2(binOpData_t* data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[10]; 
  double m1r[30]; 
  double m2r[10]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m0r[4] = 0.0; 
    m0r[5] = 0.0; 
    m0r[6] = 0.0; 
    m0r[7] = 0.0; 
    m0r[8] = 0.0; 
    m0r[9] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    m1r[9] = 0.0; 
    m1r[10] = m1[10]; 
    m1r[11] = 0.0; 
    m1r[12] = 0.0; 
    m1r[13] = 0.0; 
    m1r[14] = 0.0; 
    m1r[15] = 0.0; 
    m1r[16] = 0.0; 
    m1r[17] = 0.0; 
    m1r[18] = 0.0; 
    m1r[19] = 0.0; 
    m1r[20] = m1[20]; 
    m1r[21] = 0.0; 
    m1r[22] = 0.0; 
    m1r[23] = 0.0; 
    m1r[24] = 0.0; 
    m1r[25] = 0.0; 
    m1r[26] = 0.0; 
    m1r[27] = 0.0; 
    m1r[28] = 0.0; 
    m1r[29] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    m2r[4] = 0.0; 
    m2r[5] = 0.0; 
    m2r[6] = 0.0; 
    m2r[7] = 0.0; 
    m2r[8] = 0.0; 
    m2r[9] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
    m0r[4] = m0[4]; 
    m0r[5] = m0[5]; 
    m0r[6] = m0[6]; 
    m0r[7] = m0[7]; 
    m0r[8] = m0[8]; 
    m0r[9] = m0[9]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m1r[4] = m1[4]; 
    m1r[5] = m1[5]; 
    m1r[6] = m1[6]; 
    m1r[7] = m1[7]; 
    m1r[8] = m1[8]; 
    m1r[9] = m1[9]; 
    m1r[10] = m1[10]; 
    m1r[11] = m1[11]; 
    m1r[12] = m1[12]; 
    m1r[13] = m1[13]; 
    m1r[14] = m1[14]; 
    m1r[15] = m1[15]; 
    m1r[16] = m1[16]; 
    m1r[17] = m1[17]; 
    m1r[18] = m1[18]; 
    m1r[19] = m1[19]; 
    m1r[20] = m1[20]; 
    m1r[21] = m1[21]; 
    m1r[22] = m1[22]; 
    m1r[23] = m1[23]; 
    m1r[24] = m1[24]; 
    m1r[25] = m1[25]; 
    m1r[26] = m1[26]; 
    m1r[27] = m1[27]; 
    m1r[28] = m1[28]; 
    m1r[29] = m1[29]; 
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
    m2r[4] = m2[4]; 
    m2r[5] = m2[5]; 
    m2r[6] = m2[6]; 
    m2r[7] = m2[7]; 
    m2r[8] = m2[8]; 
    m2r[9] = m2[9]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(40,40); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(0,1) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(0,2) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(0,3) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(0,4) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(0,5) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(0,6) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(0,7) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(0,8) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(0,9) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(1,0) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(1,1) = 0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(1,2) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(1,3) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(1,4) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(1,5) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(1,7) = 0.3162277660168379*m0r[1]; 
  data->AEM_S(2,0) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(2,1) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(2,2) = 0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(2,3) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(2,4) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(2,6) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(2,8) = 0.3162277660168379*m0r[2]; 
  data->AEM_S(3,0) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(3,1) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(3,2) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(3,3) = 0.3162277660168379*m0r[9]+0.3535533905932737*m0r[0]; 
  data->AEM_S(3,5) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(3,6) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(3,9) = 0.3162277660168379*m0r[3]; 
  data->AEM_S(4,0) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(4,1) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(4,2) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(4,4) = 0.3162277660168379*m0r[8]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(4,5) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(4,6) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(4,7) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(4,8) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(5,0) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(5,1) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(5,3) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(5,4) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(5,5) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(5,6) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(5,7) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(5,9) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(6,0) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(6,2) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(6,3) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(6,4) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(6,5) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(6,6) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(6,8) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(6,9) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(7,0) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(7,1) = 0.3162277660168379*m0r[1]; 
  data->AEM_S(7,4) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(7,5) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(7,7) = 0.2258769757263128*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(8,0) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(8,2) = 0.3162277660168379*m0r[2]; 
  data->AEM_S(8,4) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(8,6) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(8,8) = 0.2258769757263128*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(9,0) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(9,3) = 0.3162277660168379*m0r[3]; 
  data->AEM_S(9,5) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(9,6) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(9,9) = 0.2258769757263128*m0r[9]+0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,30) = -0.3535533905932737*cM[0]; 
  data->AEM_S(0,31) = -0.3535533905932737*cM[1]; 
  data->AEM_S(0,32) = -0.3535533905932737*cM[2]; 
  data->AEM_S(0,33) = -0.3535533905932737*cM[3]; 
  data->AEM_S(0,34) = -0.3535533905932737*cM[4]; 
  data->AEM_S(0,35) = -0.3535533905932737*cM[5]; 
  data->AEM_S(0,36) = -0.3535533905932737*cM[6]; 
  data->AEM_S(0,37) = -0.3535533905932737*cM[7]; 
  data->AEM_S(0,38) = -0.3535533905932737*cM[8]; 
  data->AEM_S(0,39) = -0.3535533905932737*cM[9]; 
  data->AEM_S(1,30) = -0.3535533905932737*cM[1]; 
  data->AEM_S(1,31) = (-0.3162277660168379*cM[7])-0.3535533905932737*cM[0]; 
  data->AEM_S(1,32) = -0.3535533905932737*cM[4]; 
  data->AEM_S(1,33) = -0.3535533905932737*cM[5]; 
  data->AEM_S(1,34) = -0.3535533905932737*cM[2]; 
  data->AEM_S(1,35) = -0.3535533905932737*cM[3]; 
  data->AEM_S(1,37) = -0.3162277660168379*cM[1]; 
  data->AEM_S(2,30) = -0.3535533905932737*cM[2]; 
  data->AEM_S(2,31) = -0.3535533905932737*cM[4]; 
  data->AEM_S(2,32) = (-0.3162277660168379*cM[8])-0.3535533905932737*cM[0]; 
  data->AEM_S(2,33) = -0.3535533905932737*cM[6]; 
  data->AEM_S(2,34) = -0.3535533905932737*cM[1]; 
  data->AEM_S(2,36) = -0.3535533905932737*cM[3]; 
  data->AEM_S(2,38) = -0.3162277660168379*cM[2]; 
  data->AEM_S(3,30) = -0.3535533905932737*cM[3]; 
  data->AEM_S(3,31) = -0.3535533905932737*cM[5]; 
  data->AEM_S(3,32) = -0.3535533905932737*cM[6]; 
  data->AEM_S(3,33) = (-0.3162277660168379*cM[9])-0.3535533905932737*cM[0]; 
  data->AEM_S(3,35) = -0.3535533905932737*cM[1]; 
  data->AEM_S(3,36) = -0.3535533905932737*cM[2]; 
  data->AEM_S(3,39) = -0.3162277660168379*cM[3]; 
  data->AEM_S(4,30) = -0.3535533905932737*cM[4]; 
  data->AEM_S(4,31) = -0.3535533905932737*cM[2]; 
  data->AEM_S(4,32) = -0.3535533905932737*cM[1]; 
  data->AEM_S(4,34) = (-0.3162277660168379*cM[8])-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(4,35) = -0.3535533905932737*cM[6]; 
  data->AEM_S(4,36) = -0.3535533905932737*cM[5]; 
  data->AEM_S(4,37) = -0.3162277660168379*cM[4]; 
  data->AEM_S(4,38) = -0.3162277660168379*cM[4]; 
  data->AEM_S(5,30) = -0.3535533905932737*cM[5]; 
  data->AEM_S(5,31) = -0.3535533905932737*cM[3]; 
  data->AEM_S(5,33) = -0.3535533905932737*cM[1]; 
  data->AEM_S(5,34) = -0.3535533905932737*cM[6]; 
  data->AEM_S(5,35) = (-0.3162277660168379*cM[9])-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(5,36) = -0.3535533905932737*cM[4]; 
  data->AEM_S(5,37) = -0.3162277660168379*cM[5]; 
  data->AEM_S(5,39) = -0.3162277660168379*cM[5]; 
  data->AEM_S(6,30) = -0.3535533905932737*cM[6]; 
  data->AEM_S(6,32) = -0.3535533905932737*cM[3]; 
  data->AEM_S(6,33) = -0.3535533905932737*cM[2]; 
  data->AEM_S(6,34) = -0.3535533905932737*cM[5]; 
  data->AEM_S(6,35) = -0.3535533905932737*cM[4]; 
  data->AEM_S(6,36) = (-0.3162277660168379*cM[9])-0.3162277660168379*cM[8]-0.3535533905932737*cM[0]; 
  data->AEM_S(6,38) = -0.3162277660168379*cM[6]; 
  data->AEM_S(6,39) = -0.3162277660168379*cM[6]; 
  data->AEM_S(7,30) = -0.3535533905932737*cM[7]; 
  data->AEM_S(7,31) = -0.3162277660168379*cM[1]; 
  data->AEM_S(7,34) = -0.3162277660168379*cM[4]; 
  data->AEM_S(7,35) = -0.3162277660168379*cM[5]; 
  data->AEM_S(7,37) = (-0.2258769757263128*cM[7])-0.3535533905932737*cM[0]; 
  data->AEM_S(8,30) = -0.3535533905932737*cM[8]; 
  data->AEM_S(8,32) = -0.3162277660168379*cM[2]; 
  data->AEM_S(8,34) = -0.3162277660168379*cM[4]; 
  data->AEM_S(8,36) = -0.3162277660168379*cM[6]; 
  data->AEM_S(8,38) = (-0.2258769757263128*cM[8])-0.3535533905932737*cM[0]; 
  data->AEM_S(9,30) = -0.3535533905932737*cM[9]; 
  data->AEM_S(9,33) = -0.3162277660168379*cM[3]; 
  data->AEM_S(9,35) = -0.3162277660168379*cM[5]; 
  data->AEM_S(9,36) = -0.3162277660168379*cM[6]; 
  data->AEM_S(9,39) = (-0.2258769757263128*cM[9])-0.3535533905932737*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(30,0) = 0.3535533905932737*m1r[0]; 
  data->AEM_S(30,1) = 0.3535533905932737*m1r[1]; 
  data->AEM_S(30,2) = 0.3535533905932737*m1r[2]; 
  data->AEM_S(30,3) = 0.3535533905932737*m1r[3]; 
  data->AEM_S(30,4) = 0.3535533905932737*m1r[4]; 
  data->AEM_S(30,5) = 0.3535533905932737*m1r[5]; 
  data->AEM_S(30,6) = 0.3535533905932737*m1r[6]; 
  data->AEM_S(30,7) = 0.3535533905932737*m1r[7]; 
  data->AEM_S(30,8) = 0.3535533905932737*m1r[8]; 
  data->AEM_S(30,9) = 0.3535533905932737*m1r[9]; 
  data->AEM_S(31,0) = 0.3535533905932737*m1r[1]; 
  data->AEM_S(31,1) = 0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(31,2) = 0.3535533905932737*m1r[4]; 
  data->AEM_S(31,3) = 0.3535533905932737*m1r[5]; 
  data->AEM_S(31,4) = 0.3535533905932737*m1r[2]; 
  data->AEM_S(31,5) = 0.3535533905932737*m1r[3]; 
  data->AEM_S(31,7) = 0.3162277660168379*m1r[1]; 
  data->AEM_S(32,0) = 0.3535533905932737*m1r[2]; 
  data->AEM_S(32,1) = 0.3535533905932737*m1r[4]; 
  data->AEM_S(32,2) = 0.3162277660168379*m1r[8]+0.3535533905932737*m1r[0]; 
  data->AEM_S(32,3) = 0.3535533905932737*m1r[6]; 
  data->AEM_S(32,4) = 0.3535533905932737*m1r[1]; 
  data->AEM_S(32,6) = 0.3535533905932737*m1r[3]; 
  data->AEM_S(32,8) = 0.3162277660168379*m1r[2]; 
  data->AEM_S(33,0) = 0.3535533905932737*m1r[3]; 
  data->AEM_S(33,1) = 0.3535533905932737*m1r[5]; 
  data->AEM_S(33,2) = 0.3535533905932737*m1r[6]; 
  data->AEM_S(33,3) = 0.3162277660168379*m1r[9]+0.3535533905932737*m1r[0]; 
  data->AEM_S(33,5) = 0.3535533905932737*m1r[1]; 
  data->AEM_S(33,6) = 0.3535533905932737*m1r[2]; 
  data->AEM_S(33,9) = 0.3162277660168379*m1r[3]; 
  data->AEM_S(34,0) = 0.3535533905932737*m1r[4]; 
  data->AEM_S(34,1) = 0.3535533905932737*m1r[2]; 
  data->AEM_S(34,2) = 0.3535533905932737*m1r[1]; 
  data->AEM_S(34,4) = 0.3162277660168379*m1r[8]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(34,5) = 0.3535533905932737*m1r[6]; 
  data->AEM_S(34,6) = 0.3535533905932737*m1r[5]; 
  data->AEM_S(34,7) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(34,8) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(35,0) = 0.3535533905932737*m1r[5]; 
  data->AEM_S(35,1) = 0.3535533905932737*m1r[3]; 
  data->AEM_S(35,3) = 0.3535533905932737*m1r[1]; 
  data->AEM_S(35,4) = 0.3535533905932737*m1r[6]; 
  data->AEM_S(35,5) = 0.3162277660168379*m1r[9]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(35,6) = 0.3535533905932737*m1r[4]; 
  data->AEM_S(35,7) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(35,9) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(36,0) = 0.3535533905932737*m1r[6]; 
  data->AEM_S(36,2) = 0.3535533905932737*m1r[3]; 
  data->AEM_S(36,3) = 0.3535533905932737*m1r[2]; 
  data->AEM_S(36,4) = 0.3535533905932737*m1r[5]; 
  data->AEM_S(36,5) = 0.3535533905932737*m1r[4]; 
  data->AEM_S(36,6) = 0.3162277660168379*m1r[9]+0.3162277660168379*m1r[8]+0.3535533905932737*m1r[0]; 
  data->AEM_S(36,8) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(36,9) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(37,0) = 0.3535533905932737*m1r[7]; 
  data->AEM_S(37,1) = 0.3162277660168379*m1r[1]; 
  data->AEM_S(37,4) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(37,5) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(37,7) = 0.2258769757263128*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(38,0) = 0.3535533905932737*m1r[8]; 
  data->AEM_S(38,2) = 0.3162277660168379*m1r[2]; 
  data->AEM_S(38,4) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(38,6) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(38,8) = 0.2258769757263128*m1r[8]+0.3535533905932737*m1r[0]; 
  data->AEM_S(39,0) = 0.3535533905932737*m1r[9]; 
  data->AEM_S(39,3) = 0.3162277660168379*m1r[3]; 
  data->AEM_S(39,5) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(39,6) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(39,9) = 0.2258769757263128*m1r[9]+0.3535533905932737*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(10,10) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(10,11) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(10,12) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(10,13) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(10,14) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(10,15) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(10,16) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(10,17) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(10,18) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(10,19) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(11,10) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(11,11) = 0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(11,12) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(11,13) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(11,14) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(11,15) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(11,17) = 0.3162277660168379*m0r[1]; 
  data->AEM_S(12,10) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(12,11) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(12,12) = 0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(12,13) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(12,14) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(12,16) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(12,18) = 0.3162277660168379*m0r[2]; 
  data->AEM_S(13,10) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(13,11) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(13,12) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(13,13) = 0.3162277660168379*m0r[9]+0.3535533905932737*m0r[0]; 
  data->AEM_S(13,15) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(13,16) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(13,19) = 0.3162277660168379*m0r[3]; 
  data->AEM_S(14,10) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(14,11) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(14,12) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(14,14) = 0.3162277660168379*m0r[8]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(14,15) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(14,16) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(14,17) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(14,18) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(15,10) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(15,11) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(15,13) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(15,14) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(15,15) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(15,16) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(15,17) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(15,19) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(16,10) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(16,12) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(16,13) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(16,14) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(16,15) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(16,16) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(16,18) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(16,19) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(17,10) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(17,11) = 0.3162277660168379*m0r[1]; 
  data->AEM_S(17,14) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(17,15) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(17,17) = 0.2258769757263128*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(18,10) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(18,12) = 0.3162277660168379*m0r[2]; 
  data->AEM_S(18,14) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(18,16) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(18,18) = 0.2258769757263128*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(19,10) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(19,13) = 0.3162277660168379*m0r[3]; 
  data->AEM_S(19,15) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(19,16) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(19,19) = 0.2258769757263128*m0r[9]+0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(10,30) = -0.3535533905932737*cM[10]; 
  data->AEM_S(10,31) = -0.3535533905932737*cM[11]; 
  data->AEM_S(10,32) = -0.3535533905932737*cM[12]; 
  data->AEM_S(10,33) = -0.3535533905932737*cM[13]; 
  data->AEM_S(10,34) = -0.3535533905932737*cM[14]; 
  data->AEM_S(10,35) = -0.3535533905932737*cM[15]; 
  data->AEM_S(10,36) = -0.3535533905932737*cM[16]; 
  data->AEM_S(10,37) = -0.3535533905932737*cM[17]; 
  data->AEM_S(10,38) = -0.3535533905932737*cM[18]; 
  data->AEM_S(10,39) = -0.3535533905932737*cM[19]; 
  data->AEM_S(11,30) = -0.3535533905932737*cM[11]; 
  data->AEM_S(11,31) = (-0.3162277660168379*cM[17])-0.3535533905932737*cM[10]; 
  data->AEM_S(11,32) = -0.3535533905932737*cM[14]; 
  data->AEM_S(11,33) = -0.3535533905932737*cM[15]; 
  data->AEM_S(11,34) = -0.3535533905932737*cM[12]; 
  data->AEM_S(11,35) = -0.3535533905932737*cM[13]; 
  data->AEM_S(11,37) = -0.3162277660168379*cM[11]; 
  data->AEM_S(12,30) = -0.3535533905932737*cM[12]; 
  data->AEM_S(12,31) = -0.3535533905932737*cM[14]; 
  data->AEM_S(12,32) = (-0.3162277660168379*cM[18])-0.3535533905932737*cM[10]; 
  data->AEM_S(12,33) = -0.3535533905932737*cM[16]; 
  data->AEM_S(12,34) = -0.3535533905932737*cM[11]; 
  data->AEM_S(12,36) = -0.3535533905932737*cM[13]; 
  data->AEM_S(12,38) = -0.3162277660168379*cM[12]; 
  data->AEM_S(13,30) = -0.3535533905932737*cM[13]; 
  data->AEM_S(13,31) = -0.3535533905932737*cM[15]; 
  data->AEM_S(13,32) = -0.3535533905932737*cM[16]; 
  data->AEM_S(13,33) = (-0.3162277660168379*cM[19])-0.3535533905932737*cM[10]; 
  data->AEM_S(13,35) = -0.3535533905932737*cM[11]; 
  data->AEM_S(13,36) = -0.3535533905932737*cM[12]; 
  data->AEM_S(13,39) = -0.3162277660168379*cM[13]; 
  data->AEM_S(14,30) = -0.3535533905932737*cM[14]; 
  data->AEM_S(14,31) = -0.3535533905932737*cM[12]; 
  data->AEM_S(14,32) = -0.3535533905932737*cM[11]; 
  data->AEM_S(14,34) = (-0.3162277660168379*cM[18])-0.3162277660168379*cM[17]-0.3535533905932737*cM[10]; 
  data->AEM_S(14,35) = -0.3535533905932737*cM[16]; 
  data->AEM_S(14,36) = -0.3535533905932737*cM[15]; 
  data->AEM_S(14,37) = -0.3162277660168379*cM[14]; 
  data->AEM_S(14,38) = -0.3162277660168379*cM[14]; 
  data->AEM_S(15,30) = -0.3535533905932737*cM[15]; 
  data->AEM_S(15,31) = -0.3535533905932737*cM[13]; 
  data->AEM_S(15,33) = -0.3535533905932737*cM[11]; 
  data->AEM_S(15,34) = -0.3535533905932737*cM[16]; 
  data->AEM_S(15,35) = (-0.3162277660168379*cM[19])-0.3162277660168379*cM[17]-0.3535533905932737*cM[10]; 
  data->AEM_S(15,36) = -0.3535533905932737*cM[14]; 
  data->AEM_S(15,37) = -0.3162277660168379*cM[15]; 
  data->AEM_S(15,39) = -0.3162277660168379*cM[15]; 
  data->AEM_S(16,30) = -0.3535533905932737*cM[16]; 
  data->AEM_S(16,32) = -0.3535533905932737*cM[13]; 
  data->AEM_S(16,33) = -0.3535533905932737*cM[12]; 
  data->AEM_S(16,34) = -0.3535533905932737*cM[15]; 
  data->AEM_S(16,35) = -0.3535533905932737*cM[14]; 
  data->AEM_S(16,36) = (-0.3162277660168379*cM[19])-0.3162277660168379*cM[18]-0.3535533905932737*cM[10]; 
  data->AEM_S(16,38) = -0.3162277660168379*cM[16]; 
  data->AEM_S(16,39) = -0.3162277660168379*cM[16]; 
  data->AEM_S(17,30) = -0.3535533905932737*cM[17]; 
  data->AEM_S(17,31) = -0.3162277660168379*cM[11]; 
  data->AEM_S(17,34) = -0.3162277660168379*cM[14]; 
  data->AEM_S(17,35) = -0.3162277660168379*cM[15]; 
  data->AEM_S(17,37) = (-0.2258769757263128*cM[17])-0.3535533905932737*cM[10]; 
  data->AEM_S(18,30) = -0.3535533905932737*cM[18]; 
  data->AEM_S(18,32) = -0.3162277660168379*cM[12]; 
  data->AEM_S(18,34) = -0.3162277660168379*cM[14]; 
  data->AEM_S(18,36) = -0.3162277660168379*cM[16]; 
  data->AEM_S(18,38) = (-0.2258769757263128*cM[18])-0.3535533905932737*cM[10]; 
  data->AEM_S(19,30) = -0.3535533905932737*cM[19]; 
  data->AEM_S(19,33) = -0.3162277660168379*cM[13]; 
  data->AEM_S(19,35) = -0.3162277660168379*cM[15]; 
  data->AEM_S(19,36) = -0.3162277660168379*cM[16]; 
  data->AEM_S(19,39) = (-0.2258769757263128*cM[19])-0.3535533905932737*cM[10]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(30,10) = 0.3535533905932737*m1r[10]; 
  data->AEM_S(30,11) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(30,12) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(30,13) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(30,14) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(30,15) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(30,16) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(30,17) = 0.3535533905932737*m1r[17]; 
  data->AEM_S(30,18) = 0.3535533905932737*m1r[18]; 
  data->AEM_S(30,19) = 0.3535533905932737*m1r[19]; 
  data->AEM_S(31,10) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(31,11) = 0.3162277660168379*m1r[17]+0.3535533905932737*m1r[10]; 
  data->AEM_S(31,12) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(31,13) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(31,14) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(31,15) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(31,17) = 0.3162277660168379*m1r[11]; 
  data->AEM_S(32,10) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(32,11) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(32,12) = 0.3162277660168379*m1r[18]+0.3535533905932737*m1r[10]; 
  data->AEM_S(32,13) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(32,14) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(32,16) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(32,18) = 0.3162277660168379*m1r[12]; 
  data->AEM_S(33,10) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(33,11) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(33,12) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(33,13) = 0.3162277660168379*m1r[19]+0.3535533905932737*m1r[10]; 
  data->AEM_S(33,15) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(33,16) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(33,19) = 0.3162277660168379*m1r[13]; 
  data->AEM_S(34,10) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(34,11) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(34,12) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(34,14) = 0.3162277660168379*m1r[18]+0.3162277660168379*m1r[17]+0.3535533905932737*m1r[10]; 
  data->AEM_S(34,15) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(34,16) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(34,17) = 0.3162277660168379*m1r[14]; 
  data->AEM_S(34,18) = 0.3162277660168379*m1r[14]; 
  data->AEM_S(35,10) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(35,11) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(35,13) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(35,14) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(35,15) = 0.3162277660168379*m1r[19]+0.3162277660168379*m1r[17]+0.3535533905932737*m1r[10]; 
  data->AEM_S(35,16) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(35,17) = 0.3162277660168379*m1r[15]; 
  data->AEM_S(35,19) = 0.3162277660168379*m1r[15]; 
  data->AEM_S(36,10) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(36,12) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(36,13) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(36,14) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(36,15) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(36,16) = 0.3162277660168379*m1r[19]+0.3162277660168379*m1r[18]+0.3535533905932737*m1r[10]; 
  data->AEM_S(36,18) = 0.3162277660168379*m1r[16]; 
  data->AEM_S(36,19) = 0.3162277660168379*m1r[16]; 
  data->AEM_S(37,10) = 0.3535533905932737*m1r[17]; 
  data->AEM_S(37,11) = 0.3162277660168379*m1r[11]; 
  data->AEM_S(37,14) = 0.3162277660168379*m1r[14]; 
  data->AEM_S(37,15) = 0.3162277660168379*m1r[15]; 
  data->AEM_S(37,17) = 0.2258769757263128*m1r[17]+0.3535533905932737*m1r[10]; 
  data->AEM_S(38,10) = 0.3535533905932737*m1r[18]; 
  data->AEM_S(38,12) = 0.3162277660168379*m1r[12]; 
  data->AEM_S(38,14) = 0.3162277660168379*m1r[14]; 
  data->AEM_S(38,16) = 0.3162277660168379*m1r[16]; 
  data->AEM_S(38,18) = 0.2258769757263128*m1r[18]+0.3535533905932737*m1r[10]; 
  data->AEM_S(39,10) = 0.3535533905932737*m1r[19]; 
  data->AEM_S(39,13) = 0.3162277660168379*m1r[13]; 
  data->AEM_S(39,15) = 0.3162277660168379*m1r[15]; 
  data->AEM_S(39,16) = 0.3162277660168379*m1r[16]; 
  data->AEM_S(39,19) = 0.2258769757263128*m1r[19]+0.3535533905932737*m1r[10]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(20,20) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(20,21) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(20,22) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(20,23) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(20,24) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(20,25) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(20,26) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(20,27) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(20,28) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(20,29) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(21,20) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(21,21) = 0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(21,22) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(21,23) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(21,24) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(21,25) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(21,27) = 0.3162277660168379*m0r[1]; 
  data->AEM_S(22,20) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(22,21) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(22,22) = 0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(22,23) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(22,24) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(22,26) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(22,28) = 0.3162277660168379*m0r[2]; 
  data->AEM_S(23,20) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(23,21) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(23,22) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(23,23) = 0.3162277660168379*m0r[9]+0.3535533905932737*m0r[0]; 
  data->AEM_S(23,25) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(23,26) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(23,29) = 0.3162277660168379*m0r[3]; 
  data->AEM_S(24,20) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(24,21) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(24,22) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(24,24) = 0.3162277660168379*m0r[8]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(24,25) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(24,26) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(24,27) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(24,28) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(25,20) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(25,21) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(25,23) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(25,24) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(25,25) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(25,26) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(25,27) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(25,29) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(26,20) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(26,22) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(26,23) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(26,24) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(26,25) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(26,26) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(26,28) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(26,29) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(27,20) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(27,21) = 0.3162277660168379*m0r[1]; 
  data->AEM_S(27,24) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(27,25) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(27,27) = 0.2258769757263128*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(28,20) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(28,22) = 0.3162277660168379*m0r[2]; 
  data->AEM_S(28,24) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(28,26) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(28,28) = 0.2258769757263128*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(29,20) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(29,23) = 0.3162277660168379*m0r[3]; 
  data->AEM_S(29,25) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(29,26) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(29,29) = 0.2258769757263128*m0r[9]+0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(20,30) = -0.3535533905932737*cM[20]; 
  data->AEM_S(20,31) = -0.3535533905932737*cM[21]; 
  data->AEM_S(20,32) = -0.3535533905932737*cM[22]; 
  data->AEM_S(20,33) = -0.3535533905932737*cM[23]; 
  data->AEM_S(20,34) = -0.3535533905932737*cM[24]; 
  data->AEM_S(20,35) = -0.3535533905932737*cM[25]; 
  data->AEM_S(20,36) = -0.3535533905932737*cM[26]; 
  data->AEM_S(20,37) = -0.3535533905932737*cM[27]; 
  data->AEM_S(20,38) = -0.3535533905932737*cM[28]; 
  data->AEM_S(20,39) = -0.3535533905932737*cM[29]; 
  data->AEM_S(21,30) = -0.3535533905932737*cM[21]; 
  data->AEM_S(21,31) = (-0.3162277660168379*cM[27])-0.3535533905932737*cM[20]; 
  data->AEM_S(21,32) = -0.3535533905932737*cM[24]; 
  data->AEM_S(21,33) = -0.3535533905932737*cM[25]; 
  data->AEM_S(21,34) = -0.3535533905932737*cM[22]; 
  data->AEM_S(21,35) = -0.3535533905932737*cM[23]; 
  data->AEM_S(21,37) = -0.3162277660168379*cM[21]; 
  data->AEM_S(22,30) = -0.3535533905932737*cM[22]; 
  data->AEM_S(22,31) = -0.3535533905932737*cM[24]; 
  data->AEM_S(22,32) = (-0.3162277660168379*cM[28])-0.3535533905932737*cM[20]; 
  data->AEM_S(22,33) = -0.3535533905932737*cM[26]; 
  data->AEM_S(22,34) = -0.3535533905932737*cM[21]; 
  data->AEM_S(22,36) = -0.3535533905932737*cM[23]; 
  data->AEM_S(22,38) = -0.3162277660168379*cM[22]; 
  data->AEM_S(23,30) = -0.3535533905932737*cM[23]; 
  data->AEM_S(23,31) = -0.3535533905932737*cM[25]; 
  data->AEM_S(23,32) = -0.3535533905932737*cM[26]; 
  data->AEM_S(23,33) = (-0.3162277660168379*cM[29])-0.3535533905932737*cM[20]; 
  data->AEM_S(23,35) = -0.3535533905932737*cM[21]; 
  data->AEM_S(23,36) = -0.3535533905932737*cM[22]; 
  data->AEM_S(23,39) = -0.3162277660168379*cM[23]; 
  data->AEM_S(24,30) = -0.3535533905932737*cM[24]; 
  data->AEM_S(24,31) = -0.3535533905932737*cM[22]; 
  data->AEM_S(24,32) = -0.3535533905932737*cM[21]; 
  data->AEM_S(24,34) = (-0.3162277660168379*cM[28])-0.3162277660168379*cM[27]-0.3535533905932737*cM[20]; 
  data->AEM_S(24,35) = -0.3535533905932737*cM[26]; 
  data->AEM_S(24,36) = -0.3535533905932737*cM[25]; 
  data->AEM_S(24,37) = -0.3162277660168379*cM[24]; 
  data->AEM_S(24,38) = -0.3162277660168379*cM[24]; 
  data->AEM_S(25,30) = -0.3535533905932737*cM[25]; 
  data->AEM_S(25,31) = -0.3535533905932737*cM[23]; 
  data->AEM_S(25,33) = -0.3535533905932737*cM[21]; 
  data->AEM_S(25,34) = -0.3535533905932737*cM[26]; 
  data->AEM_S(25,35) = (-0.3162277660168379*cM[29])-0.3162277660168379*cM[27]-0.3535533905932737*cM[20]; 
  data->AEM_S(25,36) = -0.3535533905932737*cM[24]; 
  data->AEM_S(25,37) = -0.3162277660168379*cM[25]; 
  data->AEM_S(25,39) = -0.3162277660168379*cM[25]; 
  data->AEM_S(26,30) = -0.3535533905932737*cM[26]; 
  data->AEM_S(26,32) = -0.3535533905932737*cM[23]; 
  data->AEM_S(26,33) = -0.3535533905932737*cM[22]; 
  data->AEM_S(26,34) = -0.3535533905932737*cM[25]; 
  data->AEM_S(26,35) = -0.3535533905932737*cM[24]; 
  data->AEM_S(26,36) = (-0.3162277660168379*cM[29])-0.3162277660168379*cM[28]-0.3535533905932737*cM[20]; 
  data->AEM_S(26,38) = -0.3162277660168379*cM[26]; 
  data->AEM_S(26,39) = -0.3162277660168379*cM[26]; 
  data->AEM_S(27,30) = -0.3535533905932737*cM[27]; 
  data->AEM_S(27,31) = -0.3162277660168379*cM[21]; 
  data->AEM_S(27,34) = -0.3162277660168379*cM[24]; 
  data->AEM_S(27,35) = -0.3162277660168379*cM[25]; 
  data->AEM_S(27,37) = (-0.2258769757263128*cM[27])-0.3535533905932737*cM[20]; 
  data->AEM_S(28,30) = -0.3535533905932737*cM[28]; 
  data->AEM_S(28,32) = -0.3162277660168379*cM[22]; 
  data->AEM_S(28,34) = -0.3162277660168379*cM[24]; 
  data->AEM_S(28,36) = -0.3162277660168379*cM[26]; 
  data->AEM_S(28,38) = (-0.2258769757263128*cM[28])-0.3535533905932737*cM[20]; 
  data->AEM_S(29,30) = -0.3535533905932737*cM[29]; 
  data->AEM_S(29,33) = -0.3162277660168379*cM[23]; 
  data->AEM_S(29,35) = -0.3162277660168379*cM[25]; 
  data->AEM_S(29,36) = -0.3162277660168379*cM[26]; 
  data->AEM_S(29,39) = (-0.2258769757263128*cM[29])-0.3535533905932737*cM[20]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(30,20) = 0.3535533905932737*m1r[20]; 
  data->AEM_S(30,21) = 0.3535533905932737*m1r[21]; 
  data->AEM_S(30,22) = 0.3535533905932737*m1r[22]; 
  data->AEM_S(30,23) = 0.3535533905932737*m1r[23]; 
  data->AEM_S(30,24) = 0.3535533905932737*m1r[24]; 
  data->AEM_S(30,25) = 0.3535533905932737*m1r[25]; 
  data->AEM_S(30,26) = 0.3535533905932737*m1r[26]; 
  data->AEM_S(30,27) = 0.3535533905932737*m1r[27]; 
  data->AEM_S(30,28) = 0.3535533905932737*m1r[28]; 
  data->AEM_S(30,29) = 0.3535533905932737*m1r[29]; 
  data->AEM_S(31,20) = 0.3535533905932737*m1r[21]; 
  data->AEM_S(31,21) = 0.3162277660168379*m1r[27]+0.3535533905932737*m1r[20]; 
  data->AEM_S(31,22) = 0.3535533905932737*m1r[24]; 
  data->AEM_S(31,23) = 0.3535533905932737*m1r[25]; 
  data->AEM_S(31,24) = 0.3535533905932737*m1r[22]; 
  data->AEM_S(31,25) = 0.3535533905932737*m1r[23]; 
  data->AEM_S(31,27) = 0.3162277660168379*m1r[21]; 
  data->AEM_S(32,20) = 0.3535533905932737*m1r[22]; 
  data->AEM_S(32,21) = 0.3535533905932737*m1r[24]; 
  data->AEM_S(32,22) = 0.3162277660168379*m1r[28]+0.3535533905932737*m1r[20]; 
  data->AEM_S(32,23) = 0.3535533905932737*m1r[26]; 
  data->AEM_S(32,24) = 0.3535533905932737*m1r[21]; 
  data->AEM_S(32,26) = 0.3535533905932737*m1r[23]; 
  data->AEM_S(32,28) = 0.3162277660168379*m1r[22]; 
  data->AEM_S(33,20) = 0.3535533905932737*m1r[23]; 
  data->AEM_S(33,21) = 0.3535533905932737*m1r[25]; 
  data->AEM_S(33,22) = 0.3535533905932737*m1r[26]; 
  data->AEM_S(33,23) = 0.3162277660168379*m1r[29]+0.3535533905932737*m1r[20]; 
  data->AEM_S(33,25) = 0.3535533905932737*m1r[21]; 
  data->AEM_S(33,26) = 0.3535533905932737*m1r[22]; 
  data->AEM_S(33,29) = 0.3162277660168379*m1r[23]; 
  data->AEM_S(34,20) = 0.3535533905932737*m1r[24]; 
  data->AEM_S(34,21) = 0.3535533905932737*m1r[22]; 
  data->AEM_S(34,22) = 0.3535533905932737*m1r[21]; 
  data->AEM_S(34,24) = 0.3162277660168379*m1r[28]+0.3162277660168379*m1r[27]+0.3535533905932737*m1r[20]; 
  data->AEM_S(34,25) = 0.3535533905932737*m1r[26]; 
  data->AEM_S(34,26) = 0.3535533905932737*m1r[25]; 
  data->AEM_S(34,27) = 0.3162277660168379*m1r[24]; 
  data->AEM_S(34,28) = 0.3162277660168379*m1r[24]; 
  data->AEM_S(35,20) = 0.3535533905932737*m1r[25]; 
  data->AEM_S(35,21) = 0.3535533905932737*m1r[23]; 
  data->AEM_S(35,23) = 0.3535533905932737*m1r[21]; 
  data->AEM_S(35,24) = 0.3535533905932737*m1r[26]; 
  data->AEM_S(35,25) = 0.3162277660168379*m1r[29]+0.3162277660168379*m1r[27]+0.3535533905932737*m1r[20]; 
  data->AEM_S(35,26) = 0.3535533905932737*m1r[24]; 
  data->AEM_S(35,27) = 0.3162277660168379*m1r[25]; 
  data->AEM_S(35,29) = 0.3162277660168379*m1r[25]; 
  data->AEM_S(36,20) = 0.3535533905932737*m1r[26]; 
  data->AEM_S(36,22) = 0.3535533905932737*m1r[23]; 
  data->AEM_S(36,23) = 0.3535533905932737*m1r[22]; 
  data->AEM_S(36,24) = 0.3535533905932737*m1r[25]; 
  data->AEM_S(36,25) = 0.3535533905932737*m1r[24]; 
  data->AEM_S(36,26) = 0.3162277660168379*m1r[29]+0.3162277660168379*m1r[28]+0.3535533905932737*m1r[20]; 
  data->AEM_S(36,28) = 0.3162277660168379*m1r[26]; 
  data->AEM_S(36,29) = 0.3162277660168379*m1r[26]; 
  data->AEM_S(37,20) = 0.3535533905932737*m1r[27]; 
  data->AEM_S(37,21) = 0.3162277660168379*m1r[21]; 
  data->AEM_S(37,24) = 0.3162277660168379*m1r[24]; 
  data->AEM_S(37,25) = 0.3162277660168379*m1r[25]; 
  data->AEM_S(37,27) = 0.2258769757263128*m1r[27]+0.3535533905932737*m1r[20]; 
  data->AEM_S(38,20) = 0.3535533905932737*m1r[28]; 
  data->AEM_S(38,22) = 0.3162277660168379*m1r[22]; 
  data->AEM_S(38,24) = 0.3162277660168379*m1r[24]; 
  data->AEM_S(38,26) = 0.3162277660168379*m1r[26]; 
  data->AEM_S(38,28) = 0.2258769757263128*m1r[28]+0.3535533905932737*m1r[20]; 
  data->AEM_S(39,20) = 0.3535533905932737*m1r[29]; 
  data->AEM_S(39,23) = 0.3162277660168379*m1r[23]; 
  data->AEM_S(39,25) = 0.3162277660168379*m1r[25]; 
  data->AEM_S(39,26) = 0.3162277660168379*m1r[26]; 
  data->AEM_S(39,29) = 0.2258769757263128*m1r[29]+0.3535533905932737*m1r[20]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(30,30) = 1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(30,31) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(30,32) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(30,33) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(30,34) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(30,35) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(30,36) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(30,37) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(30,38) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  data->AEM_S(30,39) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  data->AEM_S(31,30) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(31,31) = 0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(31,32) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(31,33) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(31,34) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(31,35) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(31,37) = 0.9486832980505137*m0r[1]-0.3162277660168379*cE[1]; 
  data->AEM_S(32,30) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(32,31) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(32,32) = 0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(32,33) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(32,34) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(32,36) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(32,38) = 0.9486832980505137*m0r[2]-0.3162277660168379*cE[2]; 
  data->AEM_S(33,30) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(33,31) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(33,32) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(33,33) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(33,35) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(33,36) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(33,39) = 0.9486832980505137*m0r[3]-0.3162277660168379*cE[3]; 
  data->AEM_S(34,30) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(34,31) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(34,32) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(34,34) = 0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(34,35) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(34,36) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(34,37) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(34,38) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(35,30) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(35,31) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(35,33) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(35,34) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(35,35) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(35,36) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(35,37) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(35,39) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(36,30) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(36,32) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(36,33) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(36,34) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(36,35) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(36,36) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(36,38) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(36,39) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(37,30) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(37,31) = 0.9486832980505137*m0r[1]-0.3162277660168379*cE[1]; 
  data->AEM_S(37,34) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(37,35) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(37,37) = 0.6776309271789384*m0r[7]-0.2258769757263128*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(38,30) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  data->AEM_S(38,32) = 0.9486832980505137*m0r[2]-0.3162277660168379*cE[2]; 
  data->AEM_S(38,34) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(38,36) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(38,38) = 0.6776309271789384*m0r[8]-0.2258769757263128*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(39,30) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  data->AEM_S(39,33) = 0.9486832980505137*m0r[3]-0.3162277660168379*cE[3]; 
  data->AEM_S(39,35) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(39,36) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(39,39) = 0.6776309271789384*m0r[9]-0.2258769757263128*cE[9]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<10,20>(0,10).setZero(); 
  data->AEM_S.block<20,10>(10,0).setZero(); 
  data->AEM_S.block<10,10>(10,20).setZero(); 
  data->AEM_S.block<10,10>(20,10).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m1r[12],m1r[13],m1r[14],m1r[15],m1r[16],m1r[17],m1r[18],m1r[19],m1r[20],m1r[21],m1r[22],m1r[23],m1r[24],m1r[25],m1r[26],m1r[27],m1r[28],m1r[29],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5],m2r[6],m2r[7],m2r[8],m2r[9]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,30,1) = data->u_S.segment<30>(0); 
 
  Eigen::Map<VectorXd>(vtSq,10,1) = data->u_S.segment<10>(30); 
 
} 
 
