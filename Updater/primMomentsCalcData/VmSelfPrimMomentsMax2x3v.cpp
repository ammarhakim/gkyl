#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmM0Star2x3vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[5]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[3]*dxvl[4]*(wr[2]-wl[2]); 
 
  out[0] += ((-0.8164965809277261*fr[3])+0.8164965809277261*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += (0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += (0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
 
} 
 
void VmM0Star2x3vMax_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[5]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[2]*dxvl[4]*(wr[3]-wl[3]); 
 
  out[0] += ((-0.8164965809277261*fr[4])+0.8164965809277261*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += (0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += (0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
 
} 
 
void VmM0Star2x3vMax_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[5]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[2]*dxvl[3]*(wr[4]-wl[4]); 
 
  out[0] += ((-0.8164965809277261*fr[5])+0.8164965809277261*fl[5]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += (0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += (0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
 
} 
 
void VmM1iM2Star2x3vMax(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2) 
{ 
  // w[5]:    Cell-center coordinates. 
  // dxv[5]:  Cell length in each direciton. 
  // intFac:  for gyrokinetics (not used in Vm). 
  // m_:      mass (not used in Vm). 
  // Bmag[3]: Magnetic field magnitude (not used in Vm). 
  // f:       Distribution function. 
  // outM1i:  Contribution to M_1^star from this cell. 
  // outM2:   Contribution to M_2^star from this cell. 
 
  const double volFact = 0.125*dxv[2]*dxv[3]*dxv[4]; 
  double wvSq[3]; 
  wvSq[0]  = w[2]*w[2]; 
  wvSq[1]  = w[3]*w[3]; 
  wvSq[2]  = w[4]*w[4]; 
  double dvSq[3]; 
  dvSq[0] = dxv[2]*dxv[2]; 
  dvSq[1] = dxv[3]*dxv[3]; 
  dvSq[2] = dxv[4]*dxv[4]; 
  double tempM0[3]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 
  tempM0[2] = 2.828427124746191*f[2]*volFact; 

  outM1i[0] += tempM0[0]*w[2]; 
  outM1i[1] += tempM0[1]*w[2]; 
  outM1i[2] += tempM0[2]*w[2]; 
  outM1i[3] += tempM0[0]*w[3]; 
  outM1i[4] += tempM0[1]*w[3]; 
  outM1i[5] += tempM0[2]*w[3]; 
  outM1i[6] += tempM0[0]*w[4]; 
  outM1i[7] += tempM0[1]*w[4]; 
  outM1i[8] += tempM0[2]*w[4]; 

  outM2[0] += (0.8164965809277261*dxv[4]*w[4]*f[5]+0.8164965809277261*dxv[3]*w[3]*f[4]+0.8164965809277261*dxv[2]*w[2]*f[3])*volFact+tempM0[0]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[1] += tempM0[1]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[2] += tempM0[2]*(wvSq[2]+wvSq[1]+wvSq[0]); 
 
} 
void VmBoundaryIntegral2x3vMax_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[6]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[1] += -1.414213562373095*fIn[1]*dS; 
    out[2] += -1.414213562373095*fIn[2]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS; 
    out[2] += 1.414213562373095*fIn[2]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[21]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[18])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[7]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[2])*dS; 
    out[3] += -1.414213562373095*fIn[6]*dS; 
    out[4] += -1.414213562373095*fIn[16]*dS; 
    out[5] += -1.414213562373095*fIn[17]*dS; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[18]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[7]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[2])*dS; 
    out[3] += 1.414213562373095*fIn[6]*dS; 
    out[4] += 1.414213562373095*fIn[16]*dS; 
    out[5] += 1.414213562373095*fIn[17]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[6]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[3]-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += (-1.414213562373095*fIn[1]*dS*vBoundary)-0.7071067811865475*fIn[1]*dxv[2]*dS; 
    out[2] += (-1.414213562373095*fIn[2]*dS*vBoundary)-0.7071067811865475*dxv[2]*fIn[2]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[3])-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS*vBoundary-0.7071067811865475*fIn[1]*dxv[2]*dS; 
    out[2] += 1.414213562373095*fIn[2]*dS*vBoundary-0.7071067811865475*dxv[2]*fIn[2]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[21]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[18])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[7]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += -1.414213562373095*fIn[6]*dS*vBoundary; 
    out[4] += -1.414213562373095*fIn[16]*dS*vBoundary; 
    out[5] += -1.414213562373095*fIn[17]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[18]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[7]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += 1.414213562373095*fIn[6]*dS*vBoundary; 
    out[4] += 1.414213562373095*fIn[16]*dS*vBoundary; 
    out[5] += 1.414213562373095*fIn[17]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[6]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  if (atLower) {
 
    out[3] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[4] += -1.414213562373095*fIn[1]*dS; 
    out[5] += -1.414213562373095*fIn[2]*dS; 
 
  } else {
 
    out[3] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[4] += 1.414213562373095*fIn[1]*dS; 
    out[5] += 1.414213562373095*fIn[2]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[21]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  if (atLower) {
 
    out[6] += ((-3.16227766016838*fIn[19])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[7] += (2.449489742783178*fIn[9]-1.414213562373095*fIn[1])*dS; 
    out[8] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[2])*dS; 
    out[9] += -1.414213562373095*fIn[6]*dS; 
    out[10] += -1.414213562373095*fIn[16]*dS; 
    out[11] += -1.414213562373095*fIn[17]*dS; 
 
  } else {
 
    out[6] += (3.16227766016838*fIn[19]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[7] += (2.449489742783178*fIn[9]+1.414213562373095*fIn[1])*dS; 
    out[8] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[2])*dS; 
    out[9] += 1.414213562373095*fIn[6]*dS; 
    out[10] += 1.414213562373095*fIn[16]*dS; 
    out[11] += 1.414213562373095*fIn[17]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[6]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[4]-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (-1.414213562373095*fIn[1]*dS*vBoundary)-0.7071067811865475*fIn[1]*dxv[3]*dS; 
    out[2] += (-1.414213562373095*fIn[2]*dS*vBoundary)-0.7071067811865475*fIn[2]*dxv[3]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[4])-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS*vBoundary-0.7071067811865475*fIn[1]*dxv[3]*dS; 
    out[2] += 1.414213562373095*fIn[2]*dS*vBoundary-0.7071067811865475*fIn[2]*dxv[3]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[21]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[19])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[9]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += -1.414213562373095*fIn[6]*dS*vBoundary; 
    out[4] += -1.414213562373095*fIn[16]*dS*vBoundary; 
    out[5] += -1.414213562373095*fIn[17]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[19]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[9]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += 1.414213562373095*fIn[6]*dS*vBoundary; 
    out[4] += 1.414213562373095*fIn[16]*dS*vBoundary; 
    out[5] += 1.414213562373095*fIn[17]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[6]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[6] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
    out[7] += -1.414213562373095*fIn[1]*dS; 
    out[8] += -1.414213562373095*fIn[2]*dS; 
 
  } else {
 
    out[6] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS; 
    out[7] += 1.414213562373095*fIn[1]*dS; 
    out[8] += 1.414213562373095*fIn[2]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[21]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[12] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
    out[13] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[1])*dS; 
    out[14] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[2])*dS; 
    out[15] += -1.414213562373095*fIn[6]*dS; 
    out[16] += -1.414213562373095*fIn[16]*dS; 
    out[17] += -1.414213562373095*fIn[17]*dS; 
 
  } else {
 
    out[12] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS; 
    out[13] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[1])*dS; 
    out[14] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[2])*dS; 
    out[15] += 1.414213562373095*fIn[6]*dS; 
    out[16] += 1.414213562373095*fIn[16]*dS; 
    out[17] += 1.414213562373095*fIn[17]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[6]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[5]-0.7071067811865475*fIn[0]*dxv[4])*dS; 
    out[1] += (-1.414213562373095*fIn[1]*dS*vBoundary)-0.7071067811865475*fIn[1]*dxv[4]*dS; 
    out[2] += (-1.414213562373095*fIn[2]*dS*vBoundary)-0.7071067811865475*fIn[2]*dxv[4]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[5])-0.7071067811865475*fIn[0]*dxv[4])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS*vBoundary-0.7071067811865475*fIn[1]*dxv[4]*dS; 
    out[2] += 1.414213562373095*fIn[2]*dS*vBoundary-0.7071067811865475*fIn[2]*dxv[4]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vMax_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[21]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += -1.414213562373095*fIn[6]*dS*vBoundary; 
    out[4] += -1.414213562373095*fIn[16]*dS*vBoundary; 
    out[5] += -1.414213562373095*fIn[17]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += 1.414213562373095*fIn[6]*dS*vBoundary; 
    out[4] += 1.414213562373095*fIn[16]*dS*vBoundary; 
    out[5] += 1.414213562373095*fIn[17]*dS*vBoundary; 
 
  }
 
} 
 
void VmSelfPrimMoments2x3vMax_P1(const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1:       moments of the distribution function. 
  // m0S,m1S,m1S: star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-0.8660254037844386*m0[2])-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0[2])-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0[2])+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0[2])+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[3]; 
  double m1r[9]; 
  double m0Sr[3]; 
  double m1Sr[9]; 
  double m2Sr[3]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = m1[3]; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    m1r[6] = m1[6]; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m0Sr[2] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = 0.0; 
    m1Sr[3] = m1S[3]; 
    m1Sr[4] = 0.0; 
    m1Sr[5] = 0.0; 
    m1Sr[6] = m1S[6]; 
    m1Sr[7] = 0.0; 
    m1Sr[8] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
    m2Sr[2] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m1r[4] = m1[4]; 
    m1r[5] = m1[5]; 
    m1r[6] = m1[6]; 
    m1r[7] = m1[7]; 
    m1r[8] = m1[8]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m0Sr[2] = m0S[2]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = m1S[5]; 
    m1Sr[6] = m1S[6]; 
    m1Sr[7] = m1S[7]; 
    m1Sr[8] = m1S[8]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
    m2Sr[2] = m2S[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(12,12); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(12);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(12);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0r[0]; 
  BigAEM(0,1) = 0.5*m0r[1]; 
  BigAEM(0,2) = 0.5*m0r[2]; 
  BigAEM(1,0) = 0.5*m0r[1]; 
  BigAEM(1,1) = 0.5*m0r[0]; 
  BigAEM(2,0) = 0.5*m0r[2]; 
  BigAEM(2,2) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,9) = -0.5*cM[0]; 
  BigAEM(0,10) = -0.5*cM[1]; 
  BigAEM(0,11) = -0.5*cM[2]; 
  BigAEM(1,9) = -0.5*cM[1]; 
  BigAEM(1,10) = -0.5*cM[0]; 
  BigAEM(2,9) = -0.5*cM[2]; 
  BigAEM(2,11) = -0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(9,0) = 0.5*m1Sr[0]; 
  BigAEM(9,1) = 0.5*m1Sr[1]; 
  BigAEM(9,2) = 0.5*m1Sr[2]; 
  BigAEM(10,0) = 0.5*m1Sr[1]; 
  BigAEM(10,1) = 0.5*m1Sr[0]; 
  BigAEM(11,0) = 0.5*m1Sr[2]; 
  BigAEM(11,2) = 0.5*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(3,3) = 0.5*m0r[0]; 
  BigAEM(3,4) = 0.5*m0r[1]; 
  BigAEM(3,5) = 0.5*m0r[2]; 
  BigAEM(4,3) = 0.5*m0r[1]; 
  BigAEM(4,4) = 0.5*m0r[0]; 
  BigAEM(5,3) = 0.5*m0r[2]; 
  BigAEM(5,5) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(3,9) = -0.5*cM[3]; 
  BigAEM(3,10) = -0.5*cM[4]; 
  BigAEM(3,11) = -0.5*cM[5]; 
  BigAEM(4,9) = -0.5*cM[4]; 
  BigAEM(4,10) = -0.5*cM[3]; 
  BigAEM(5,9) = -0.5*cM[5]; 
  BigAEM(5,11) = -0.5*cM[3]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(9,3) = 0.5*m1Sr[3]; 
  BigAEM(9,4) = 0.5*m1Sr[4]; 
  BigAEM(9,5) = 0.5*m1Sr[5]; 
  BigAEM(10,3) = 0.5*m1Sr[4]; 
  BigAEM(10,4) = 0.5*m1Sr[3]; 
  BigAEM(11,3) = 0.5*m1Sr[5]; 
  BigAEM(11,5) = 0.5*m1Sr[3]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  BigAEM(6,6) = 0.5*m0r[0]; 
  BigAEM(6,7) = 0.5*m0r[1]; 
  BigAEM(6,8) = 0.5*m0r[2]; 
  BigAEM(7,6) = 0.5*m0r[1]; 
  BigAEM(7,7) = 0.5*m0r[0]; 
  BigAEM(8,6) = 0.5*m0r[2]; 
  BigAEM(8,8) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  BigAEM(6,9) = -0.5*cM[6]; 
  BigAEM(6,10) = -0.5*cM[7]; 
  BigAEM(6,11) = -0.5*cM[8]; 
  BigAEM(7,9) = -0.5*cM[7]; 
  BigAEM(7,10) = -0.5*cM[6]; 
  BigAEM(8,9) = -0.5*cM[8]; 
  BigAEM(8,11) = -0.5*cM[6]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  BigAEM(9,6) = 0.5*m1Sr[6]; 
  BigAEM(9,7) = 0.5*m1Sr[7]; 
  BigAEM(9,8) = 0.5*m1Sr[8]; 
  BigAEM(10,6) = 0.5*m1Sr[7]; 
  BigAEM(10,7) = 0.5*m1Sr[6]; 
  BigAEM(11,6) = 0.5*m1Sr[8]; 
  BigAEM(11,8) = 0.5*m1Sr[6]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(9,9) = 0.5*m0Sr[0]-0.5*cE[0]; 
  BigAEM(9,10) = 0.5*m0Sr[1]-0.5*cE[1]; 
  BigAEM(9,11) = 0.5*m0Sr[2]-0.5*cE[2]; 
  BigAEM(10,9) = 0.5*m0Sr[1]-0.5*cE[1]; 
  BigAEM(10,10) = 0.5*m0Sr[0]-0.5*cE[0]; 
  BigAEM(11,9) = 0.5*m0Sr[2]-0.5*cE[2]; 
  BigAEM(11,11) = 0.5*m0Sr[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<3,6>(0,3).setZero(); 
  BigAEM.block<6,3>(3,0).setZero(); 
  BigAEM.block<3,3>(3,6).setZero(); 
  BigAEM.block<3,3>(6,3).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m2Sr[0],m2Sr[1],m2Sr[2]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,9,1) = xEV.segment<9>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = xEV.segment<3>(9); 
 
} 
 
void VmSelfPrimMoments2x3vMax_P2(const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[6]; 
  double m1r[18]; 
  double m2r[6]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m0r[4] = 0.0; 
    m0r[5] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    m1r[6] = m1[6]; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    m1r[12] = m1[12]; 
    m1r[13] = 0.0; 
    m1r[14] = 0.0; 
    m1r[15] = 0.0; 
    m1r[16] = 0.0; 
    m1r[17] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    m2r[4] = 0.0; 
    m2r[5] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
    m0r[4] = m0[4]; 
    m0r[5] = m0[5]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
    m2r[4] = m2[4]; 
    m2r[5] = m2[5]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(24,24); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(24);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(24);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0r[0]; 
  BigAEM(0,1) = 0.5*m0r[1]; 
  BigAEM(0,2) = 0.5*m0r[2]; 
  BigAEM(0,3) = 0.5*m0r[3]; 
  BigAEM(0,4) = 0.5*m0r[4]; 
  BigAEM(0,5) = 0.5*m0r[5]; 
  BigAEM(1,0) = 0.5*m0r[1]; 
  BigAEM(1,1) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(1,2) = 0.5*m0r[3]; 
  BigAEM(1,3) = 0.5*m0r[2]; 
  BigAEM(1,4) = 0.4472135954999579*m0r[1]; 
  BigAEM(2,0) = 0.5*m0r[2]; 
  BigAEM(2,1) = 0.5*m0r[3]; 
  BigAEM(2,2) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(2,3) = 0.5*m0r[1]; 
  BigAEM(2,5) = 0.4472135954999579*m0r[2]; 
  BigAEM(3,0) = 0.5*m0r[3]; 
  BigAEM(3,1) = 0.5*m0r[2]; 
  BigAEM(3,2) = 0.5*m0r[1]; 
  BigAEM(3,3) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(3,4) = 0.4472135954999579*m0r[3]; 
  BigAEM(3,5) = 0.4472135954999579*m0r[3]; 
  BigAEM(4,0) = 0.5*m0r[4]; 
  BigAEM(4,1) = 0.4472135954999579*m0r[1]; 
  BigAEM(4,3) = 0.4472135954999579*m0r[3]; 
  BigAEM(4,4) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(5,0) = 0.5*m0r[5]; 
  BigAEM(5,2) = 0.4472135954999579*m0r[2]; 
  BigAEM(5,3) = 0.4472135954999579*m0r[3]; 
  BigAEM(5,5) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,18) = -0.5*cM[0]; 
  BigAEM(0,19) = -0.5*cM[1]; 
  BigAEM(0,20) = -0.5*cM[2]; 
  BigAEM(0,21) = -0.5*cM[3]; 
  BigAEM(0,22) = -0.5*cM[4]; 
  BigAEM(0,23) = -0.5*cM[5]; 
  BigAEM(1,18) = -0.5*cM[1]; 
  BigAEM(1,19) = (-0.4472135954999579*cM[4])-0.5*cM[0]; 
  BigAEM(1,20) = -0.5*cM[3]; 
  BigAEM(1,21) = -0.5*cM[2]; 
  BigAEM(1,22) = -0.4472135954999579*cM[1]; 
  BigAEM(2,18) = -0.5*cM[2]; 
  BigAEM(2,19) = -0.5*cM[3]; 
  BigAEM(2,20) = (-0.4472135954999579*cM[5])-0.5*cM[0]; 
  BigAEM(2,21) = -0.5*cM[1]; 
  BigAEM(2,23) = -0.4472135954999579*cM[2]; 
  BigAEM(3,18) = -0.5*cM[3]; 
  BigAEM(3,19) = -0.5*cM[2]; 
  BigAEM(3,20) = -0.5*cM[1]; 
  BigAEM(3,21) = (-0.4472135954999579*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(3,22) = -0.4472135954999579*cM[3]; 
  BigAEM(3,23) = -0.4472135954999579*cM[3]; 
  BigAEM(4,18) = -0.5*cM[4]; 
  BigAEM(4,19) = -0.4472135954999579*cM[1]; 
  BigAEM(4,21) = -0.4472135954999579*cM[3]; 
  BigAEM(4,22) = (-0.31943828249997*cM[4])-0.5*cM[0]; 
  BigAEM(5,18) = -0.5*cM[5]; 
  BigAEM(5,20) = -0.4472135954999579*cM[2]; 
  BigAEM(5,21) = -0.4472135954999579*cM[3]; 
  BigAEM(5,23) = (-0.31943828249997*cM[5])-0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(18,0) = 0.5*m1r[0]; 
  BigAEM(18,1) = 0.5*m1r[1]; 
  BigAEM(18,2) = 0.5*m1r[2]; 
  BigAEM(18,3) = 0.5*m1r[3]; 
  BigAEM(18,4) = 0.5*m1r[4]; 
  BigAEM(18,5) = 0.5*m1r[5]; 
  BigAEM(19,0) = 0.5*m1r[1]; 
  BigAEM(19,1) = 0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(19,2) = 0.5*m1r[3]; 
  BigAEM(19,3) = 0.5*m1r[2]; 
  BigAEM(19,4) = 0.4472135954999579*m1r[1]; 
  BigAEM(20,0) = 0.5*m1r[2]; 
  BigAEM(20,1) = 0.5*m1r[3]; 
  BigAEM(20,2) = 0.4472135954999579*m1r[5]+0.5*m1r[0]; 
  BigAEM(20,3) = 0.5*m1r[1]; 
  BigAEM(20,5) = 0.4472135954999579*m1r[2]; 
  BigAEM(21,0) = 0.5*m1r[3]; 
  BigAEM(21,1) = 0.5*m1r[2]; 
  BigAEM(21,2) = 0.5*m1r[1]; 
  BigAEM(21,3) = 0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(21,4) = 0.4472135954999579*m1r[3]; 
  BigAEM(21,5) = 0.4472135954999579*m1r[3]; 
  BigAEM(22,0) = 0.5*m1r[4]; 
  BigAEM(22,1) = 0.4472135954999579*m1r[1]; 
  BigAEM(22,3) = 0.4472135954999579*m1r[3]; 
  BigAEM(22,4) = 0.31943828249997*m1r[4]+0.5*m1r[0]; 
  BigAEM(23,0) = 0.5*m1r[5]; 
  BigAEM(23,2) = 0.4472135954999579*m1r[2]; 
  BigAEM(23,3) = 0.4472135954999579*m1r[3]; 
  BigAEM(23,5) = 0.31943828249997*m1r[5]+0.5*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(6,6) = 0.5*m0r[0]; 
  BigAEM(6,7) = 0.5*m0r[1]; 
  BigAEM(6,8) = 0.5*m0r[2]; 
  BigAEM(6,9) = 0.5*m0r[3]; 
  BigAEM(6,10) = 0.5*m0r[4]; 
  BigAEM(6,11) = 0.5*m0r[5]; 
  BigAEM(7,6) = 0.5*m0r[1]; 
  BigAEM(7,7) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(7,8) = 0.5*m0r[3]; 
  BigAEM(7,9) = 0.5*m0r[2]; 
  BigAEM(7,10) = 0.4472135954999579*m0r[1]; 
  BigAEM(8,6) = 0.5*m0r[2]; 
  BigAEM(8,7) = 0.5*m0r[3]; 
  BigAEM(8,8) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(8,9) = 0.5*m0r[1]; 
  BigAEM(8,11) = 0.4472135954999579*m0r[2]; 
  BigAEM(9,6) = 0.5*m0r[3]; 
  BigAEM(9,7) = 0.5*m0r[2]; 
  BigAEM(9,8) = 0.5*m0r[1]; 
  BigAEM(9,9) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(9,10) = 0.4472135954999579*m0r[3]; 
  BigAEM(9,11) = 0.4472135954999579*m0r[3]; 
  BigAEM(10,6) = 0.5*m0r[4]; 
  BigAEM(10,7) = 0.4472135954999579*m0r[1]; 
  BigAEM(10,9) = 0.4472135954999579*m0r[3]; 
  BigAEM(10,10) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(11,6) = 0.5*m0r[5]; 
  BigAEM(11,8) = 0.4472135954999579*m0r[2]; 
  BigAEM(11,9) = 0.4472135954999579*m0r[3]; 
  BigAEM(11,11) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(6,18) = -0.5*cM[6]; 
  BigAEM(6,19) = -0.5*cM[7]; 
  BigAEM(6,20) = -0.5*cM[8]; 
  BigAEM(6,21) = -0.5*cM[9]; 
  BigAEM(6,22) = -0.5*cM[10]; 
  BigAEM(6,23) = -0.5*cM[11]; 
  BigAEM(7,18) = -0.5*cM[7]; 
  BigAEM(7,19) = (-0.4472135954999579*cM[10])-0.5*cM[6]; 
  BigAEM(7,20) = -0.5*cM[9]; 
  BigAEM(7,21) = -0.5*cM[8]; 
  BigAEM(7,22) = -0.4472135954999579*cM[7]; 
  BigAEM(8,18) = -0.5*cM[8]; 
  BigAEM(8,19) = -0.5*cM[9]; 
  BigAEM(8,20) = (-0.4472135954999579*cM[11])-0.5*cM[6]; 
  BigAEM(8,21) = -0.5*cM[7]; 
  BigAEM(8,23) = -0.4472135954999579*cM[8]; 
  BigAEM(9,18) = -0.5*cM[9]; 
  BigAEM(9,19) = -0.5*cM[8]; 
  BigAEM(9,20) = -0.5*cM[7]; 
  BigAEM(9,21) = (-0.4472135954999579*cM[11])-0.4472135954999579*cM[10]-0.5*cM[6]; 
  BigAEM(9,22) = -0.4472135954999579*cM[9]; 
  BigAEM(9,23) = -0.4472135954999579*cM[9]; 
  BigAEM(10,18) = -0.5*cM[10]; 
  BigAEM(10,19) = -0.4472135954999579*cM[7]; 
  BigAEM(10,21) = -0.4472135954999579*cM[9]; 
  BigAEM(10,22) = (-0.31943828249997*cM[10])-0.5*cM[6]; 
  BigAEM(11,18) = -0.5*cM[11]; 
  BigAEM(11,20) = -0.4472135954999579*cM[8]; 
  BigAEM(11,21) = -0.4472135954999579*cM[9]; 
  BigAEM(11,23) = (-0.31943828249997*cM[11])-0.5*cM[6]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(18,6) = 0.5*m1r[6]; 
  BigAEM(18,7) = 0.5*m1r[7]; 
  BigAEM(18,8) = 0.5*m1r[8]; 
  BigAEM(18,9) = 0.5*m1r[9]; 
  BigAEM(18,10) = 0.5*m1r[10]; 
  BigAEM(18,11) = 0.5*m1r[11]; 
  BigAEM(19,6) = 0.5*m1r[7]; 
  BigAEM(19,7) = 0.4472135954999579*m1r[10]+0.5*m1r[6]; 
  BigAEM(19,8) = 0.5*m1r[9]; 
  BigAEM(19,9) = 0.5*m1r[8]; 
  BigAEM(19,10) = 0.4472135954999579*m1r[7]; 
  BigAEM(20,6) = 0.5*m1r[8]; 
  BigAEM(20,7) = 0.5*m1r[9]; 
  BigAEM(20,8) = 0.4472135954999579*m1r[11]+0.5*m1r[6]; 
  BigAEM(20,9) = 0.5*m1r[7]; 
  BigAEM(20,11) = 0.4472135954999579*m1r[8]; 
  BigAEM(21,6) = 0.5*m1r[9]; 
  BigAEM(21,7) = 0.5*m1r[8]; 
  BigAEM(21,8) = 0.5*m1r[7]; 
  BigAEM(21,9) = 0.4472135954999579*m1r[11]+0.4472135954999579*m1r[10]+0.5*m1r[6]; 
  BigAEM(21,10) = 0.4472135954999579*m1r[9]; 
  BigAEM(21,11) = 0.4472135954999579*m1r[9]; 
  BigAEM(22,6) = 0.5*m1r[10]; 
  BigAEM(22,7) = 0.4472135954999579*m1r[7]; 
  BigAEM(22,9) = 0.4472135954999579*m1r[9]; 
  BigAEM(22,10) = 0.31943828249997*m1r[10]+0.5*m1r[6]; 
  BigAEM(23,6) = 0.5*m1r[11]; 
  BigAEM(23,8) = 0.4472135954999579*m1r[8]; 
  BigAEM(23,9) = 0.4472135954999579*m1r[9]; 
  BigAEM(23,11) = 0.31943828249997*m1r[11]+0.5*m1r[6]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  BigAEM(12,12) = 0.5*m0r[0]; 
  BigAEM(12,13) = 0.5*m0r[1]; 
  BigAEM(12,14) = 0.5*m0r[2]; 
  BigAEM(12,15) = 0.5*m0r[3]; 
  BigAEM(12,16) = 0.5*m0r[4]; 
  BigAEM(12,17) = 0.5*m0r[5]; 
  BigAEM(13,12) = 0.5*m0r[1]; 
  BigAEM(13,13) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(13,14) = 0.5*m0r[3]; 
  BigAEM(13,15) = 0.5*m0r[2]; 
  BigAEM(13,16) = 0.4472135954999579*m0r[1]; 
  BigAEM(14,12) = 0.5*m0r[2]; 
  BigAEM(14,13) = 0.5*m0r[3]; 
  BigAEM(14,14) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(14,15) = 0.5*m0r[1]; 
  BigAEM(14,17) = 0.4472135954999579*m0r[2]; 
  BigAEM(15,12) = 0.5*m0r[3]; 
  BigAEM(15,13) = 0.5*m0r[2]; 
  BigAEM(15,14) = 0.5*m0r[1]; 
  BigAEM(15,15) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(15,16) = 0.4472135954999579*m0r[3]; 
  BigAEM(15,17) = 0.4472135954999579*m0r[3]; 
  BigAEM(16,12) = 0.5*m0r[4]; 
  BigAEM(16,13) = 0.4472135954999579*m0r[1]; 
  BigAEM(16,15) = 0.4472135954999579*m0r[3]; 
  BigAEM(16,16) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(17,12) = 0.5*m0r[5]; 
  BigAEM(17,14) = 0.4472135954999579*m0r[2]; 
  BigAEM(17,15) = 0.4472135954999579*m0r[3]; 
  BigAEM(17,17) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  BigAEM(12,18) = -0.5*cM[12]; 
  BigAEM(12,19) = -0.5*cM[13]; 
  BigAEM(12,20) = -0.5*cM[14]; 
  BigAEM(12,21) = -0.5*cM[15]; 
  BigAEM(12,22) = -0.5*cM[16]; 
  BigAEM(12,23) = -0.5*cM[17]; 
  BigAEM(13,18) = -0.5*cM[13]; 
  BigAEM(13,19) = (-0.4472135954999579*cM[16])-0.5*cM[12]; 
  BigAEM(13,20) = -0.5*cM[15]; 
  BigAEM(13,21) = -0.5*cM[14]; 
  BigAEM(13,22) = -0.4472135954999579*cM[13]; 
  BigAEM(14,18) = -0.5*cM[14]; 
  BigAEM(14,19) = -0.5*cM[15]; 
  BigAEM(14,20) = (-0.4472135954999579*cM[17])-0.5*cM[12]; 
  BigAEM(14,21) = -0.5*cM[13]; 
  BigAEM(14,23) = -0.4472135954999579*cM[14]; 
  BigAEM(15,18) = -0.5*cM[15]; 
  BigAEM(15,19) = -0.5*cM[14]; 
  BigAEM(15,20) = -0.5*cM[13]; 
  BigAEM(15,21) = (-0.4472135954999579*cM[17])-0.4472135954999579*cM[16]-0.5*cM[12]; 
  BigAEM(15,22) = -0.4472135954999579*cM[15]; 
  BigAEM(15,23) = -0.4472135954999579*cM[15]; 
  BigAEM(16,18) = -0.5*cM[16]; 
  BigAEM(16,19) = -0.4472135954999579*cM[13]; 
  BigAEM(16,21) = -0.4472135954999579*cM[15]; 
  BigAEM(16,22) = (-0.31943828249997*cM[16])-0.5*cM[12]; 
  BigAEM(17,18) = -0.5*cM[17]; 
  BigAEM(17,20) = -0.4472135954999579*cM[14]; 
  BigAEM(17,21) = -0.4472135954999579*cM[15]; 
  BigAEM(17,23) = (-0.31943828249997*cM[17])-0.5*cM[12]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  BigAEM(18,12) = 0.5*m1r[12]; 
  BigAEM(18,13) = 0.5*m1r[13]; 
  BigAEM(18,14) = 0.5*m1r[14]; 
  BigAEM(18,15) = 0.5*m1r[15]; 
  BigAEM(18,16) = 0.5*m1r[16]; 
  BigAEM(18,17) = 0.5*m1r[17]; 
  BigAEM(19,12) = 0.5*m1r[13]; 
  BigAEM(19,13) = 0.4472135954999579*m1r[16]+0.5*m1r[12]; 
  BigAEM(19,14) = 0.5*m1r[15]; 
  BigAEM(19,15) = 0.5*m1r[14]; 
  BigAEM(19,16) = 0.4472135954999579*m1r[13]; 
  BigAEM(20,12) = 0.5*m1r[14]; 
  BigAEM(20,13) = 0.5*m1r[15]; 
  BigAEM(20,14) = 0.4472135954999579*m1r[17]+0.5*m1r[12]; 
  BigAEM(20,15) = 0.5*m1r[13]; 
  BigAEM(20,17) = 0.4472135954999579*m1r[14]; 
  BigAEM(21,12) = 0.5*m1r[15]; 
  BigAEM(21,13) = 0.5*m1r[14]; 
  BigAEM(21,14) = 0.5*m1r[13]; 
  BigAEM(21,15) = 0.4472135954999579*m1r[17]+0.4472135954999579*m1r[16]+0.5*m1r[12]; 
  BigAEM(21,16) = 0.4472135954999579*m1r[15]; 
  BigAEM(21,17) = 0.4472135954999579*m1r[15]; 
  BigAEM(22,12) = 0.5*m1r[16]; 
  BigAEM(22,13) = 0.4472135954999579*m1r[13]; 
  BigAEM(22,15) = 0.4472135954999579*m1r[15]; 
  BigAEM(22,16) = 0.31943828249997*m1r[16]+0.5*m1r[12]; 
  BigAEM(23,12) = 0.5*m1r[17]; 
  BigAEM(23,14) = 0.4472135954999579*m1r[14]; 
  BigAEM(23,15) = 0.4472135954999579*m1r[15]; 
  BigAEM(23,17) = 0.31943828249997*m1r[17]+0.5*m1r[12]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(18,18) = 1.5*m0r[0]-0.5*cE[0]; 
  BigAEM(18,19) = 1.5*m0r[1]-0.5*cE[1]; 
  BigAEM(18,20) = 1.5*m0r[2]-0.5*cE[2]; 
  BigAEM(18,21) = 1.5*m0r[3]-0.5*cE[3]; 
  BigAEM(18,22) = 1.5*m0r[4]-0.5*cE[4]; 
  BigAEM(18,23) = 1.5*m0r[5]-0.5*cE[5]; 
  BigAEM(19,18) = 1.5*m0r[1]-0.5*cE[1]; 
  BigAEM(19,19) = 1.341640786499874*m0r[4]-0.4472135954999579*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
  BigAEM(19,20) = 1.5*m0r[3]-0.5*cE[3]; 
  BigAEM(19,21) = 1.5*m0r[2]-0.5*cE[2]; 
  BigAEM(19,22) = 1.341640786499874*m0r[1]-0.4472135954999579*cE[1]; 
  BigAEM(20,18) = 1.5*m0r[2]-0.5*cE[2]; 
  BigAEM(20,19) = 1.5*m0r[3]-0.5*cE[3]; 
  BigAEM(20,20) = 1.341640786499874*m0r[5]-0.4472135954999579*cE[5]+1.5*m0r[0]-0.5*cE[0]; 
  BigAEM(20,21) = 1.5*m0r[1]-0.5*cE[1]; 
  BigAEM(20,23) = 1.341640786499874*m0r[2]-0.4472135954999579*cE[2]; 
  BigAEM(21,18) = 1.5*m0r[3]-0.5*cE[3]; 
  BigAEM(21,19) = 1.5*m0r[2]-0.5*cE[2]; 
  BigAEM(21,20) = 1.5*m0r[1]-0.5*cE[1]; 
  BigAEM(21,21) = 1.341640786499874*m0r[5]-0.4472135954999579*cE[5]+1.341640786499874*m0r[4]-0.4472135954999579*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
  BigAEM(21,22) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(21,23) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(22,18) = 1.5*m0r[4]-0.5*cE[4]; 
  BigAEM(22,19) = 1.341640786499874*m0r[1]-0.4472135954999579*cE[1]; 
  BigAEM(22,21) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(22,22) = 0.9583148474999099*m0r[4]-0.31943828249997*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
  BigAEM(23,18) = 1.5*m0r[5]-0.5*cE[5]; 
  BigAEM(23,20) = 1.341640786499874*m0r[2]-0.4472135954999579*cE[2]; 
  BigAEM(23,21) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(23,23) = 0.9583148474999099*m0r[5]-0.31943828249997*cE[5]+1.5*m0r[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<6,12>(0,6).setZero(); 
  BigAEM.block<12,6>(6,0).setZero(); 
  BigAEM.block<6,6>(6,12).setZero(); 
  BigAEM.block<6,6>(12,6).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m1r[12],m1r[13],m1r[14],m1r[15],m1r[16],m1r[17],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,18,1) = xEV.segment<18>(0); 
 
  Eigen::Map<VectorXd>(vtSq,6,1) = xEV.segment<6>(18); 
 
} 
 
