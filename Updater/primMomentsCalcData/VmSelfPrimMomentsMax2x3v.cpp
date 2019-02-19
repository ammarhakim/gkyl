
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
 
void VmSelfPrimMoments2x3vMax_P1(binOpData_t* data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
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
  data->AEM_S = Eigen::MatrixXd::Zero(12,12); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.5*m0r[0]; 
  data->AEM_S(0,1) = 0.5*m0r[1]; 
  data->AEM_S(0,2) = 0.5*m0r[2]; 
  data->AEM_S(1,0) = 0.5*m0r[1]; 
  data->AEM_S(1,1) = 0.5*m0r[0]; 
  data->AEM_S(2,0) = 0.5*m0r[2]; 
  data->AEM_S(2,2) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,9) = -0.5*cM[0]; 
  data->AEM_S(0,10) = -0.5*cM[1]; 
  data->AEM_S(0,11) = -0.5*cM[2]; 
  data->AEM_S(1,9) = -0.5*cM[1]; 
  data->AEM_S(1,10) = -0.5*cM[0]; 
  data->AEM_S(2,9) = -0.5*cM[2]; 
  data->AEM_S(2,11) = -0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(9,0) = 0.5*m1Sr[0]; 
  data->AEM_S(9,1) = 0.5*m1Sr[1]; 
  data->AEM_S(9,2) = 0.5*m1Sr[2]; 
  data->AEM_S(10,0) = 0.5*m1Sr[1]; 
  data->AEM_S(10,1) = 0.5*m1Sr[0]; 
  data->AEM_S(11,0) = 0.5*m1Sr[2]; 
  data->AEM_S(11,2) = 0.5*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(3,3) = 0.5*m0r[0]; 
  data->AEM_S(3,4) = 0.5*m0r[1]; 
  data->AEM_S(3,5) = 0.5*m0r[2]; 
  data->AEM_S(4,3) = 0.5*m0r[1]; 
  data->AEM_S(4,4) = 0.5*m0r[0]; 
  data->AEM_S(5,3) = 0.5*m0r[2]; 
  data->AEM_S(5,5) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(3,9) = -0.5*cM[3]; 
  data->AEM_S(3,10) = -0.5*cM[4]; 
  data->AEM_S(3,11) = -0.5*cM[5]; 
  data->AEM_S(4,9) = -0.5*cM[4]; 
  data->AEM_S(4,10) = -0.5*cM[3]; 
  data->AEM_S(5,9) = -0.5*cM[5]; 
  data->AEM_S(5,11) = -0.5*cM[3]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(9,3) = 0.5*m1Sr[3]; 
  data->AEM_S(9,4) = 0.5*m1Sr[4]; 
  data->AEM_S(9,5) = 0.5*m1Sr[5]; 
  data->AEM_S(10,3) = 0.5*m1Sr[4]; 
  data->AEM_S(10,4) = 0.5*m1Sr[3]; 
  data->AEM_S(11,3) = 0.5*m1Sr[5]; 
  data->AEM_S(11,5) = 0.5*m1Sr[3]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(6,6) = 0.5*m0r[0]; 
  data->AEM_S(6,7) = 0.5*m0r[1]; 
  data->AEM_S(6,8) = 0.5*m0r[2]; 
  data->AEM_S(7,6) = 0.5*m0r[1]; 
  data->AEM_S(7,7) = 0.5*m0r[0]; 
  data->AEM_S(8,6) = 0.5*m0r[2]; 
  data->AEM_S(8,8) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(6,9) = -0.5*cM[6]; 
  data->AEM_S(6,10) = -0.5*cM[7]; 
  data->AEM_S(6,11) = -0.5*cM[8]; 
  data->AEM_S(7,9) = -0.5*cM[7]; 
  data->AEM_S(7,10) = -0.5*cM[6]; 
  data->AEM_S(8,9) = -0.5*cM[8]; 
  data->AEM_S(8,11) = -0.5*cM[6]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(9,6) = 0.5*m1Sr[6]; 
  data->AEM_S(9,7) = 0.5*m1Sr[7]; 
  data->AEM_S(9,8) = 0.5*m1Sr[8]; 
  data->AEM_S(10,6) = 0.5*m1Sr[7]; 
  data->AEM_S(10,7) = 0.5*m1Sr[6]; 
  data->AEM_S(11,6) = 0.5*m1Sr[8]; 
  data->AEM_S(11,8) = 0.5*m1Sr[6]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(9,9) = 0.5*m0Sr[0]-0.5*cE[0]; 
  data->AEM_S(9,10) = 0.5*m0Sr[1]-0.5*cE[1]; 
  data->AEM_S(9,11) = 0.5*m0Sr[2]-0.5*cE[2]; 
  data->AEM_S(10,9) = 0.5*m0Sr[1]-0.5*cE[1]; 
  data->AEM_S(10,10) = 0.5*m0Sr[0]-0.5*cE[0]; 
  data->AEM_S(11,9) = 0.5*m0Sr[2]-0.5*cE[2]; 
  data->AEM_S(11,11) = 0.5*m0Sr[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<3,6>(0,3).setZero(); 
  data->AEM_S.block<6,3>(3,0).setZero(); 
  data->AEM_S.block<3,3>(3,6).setZero(); 
  data->AEM_S.block<3,3>(6,3).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m2Sr[0],m2Sr[1],m2Sr[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,9,1) = data->u_S.segment<9>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = data->u_S.segment<3>(9); 
 
} 
 
void VmSelfPrimMoments2x3vMax_P2(binOpData_t* data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
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
  data->AEM_S = Eigen::MatrixXd::Zero(24,24); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.5*m0r[0]; 
  data->AEM_S(0,1) = 0.5*m0r[1]; 
  data->AEM_S(0,2) = 0.5*m0r[2]; 
  data->AEM_S(0,3) = 0.5*m0r[3]; 
  data->AEM_S(0,4) = 0.5*m0r[4]; 
  data->AEM_S(0,5) = 0.5*m0r[5]; 
  data->AEM_S(1,0) = 0.5*m0r[1]; 
  data->AEM_S(1,1) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(1,2) = 0.5*m0r[3]; 
  data->AEM_S(1,3) = 0.5*m0r[2]; 
  data->AEM_S(1,4) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(2,0) = 0.5*m0r[2]; 
  data->AEM_S(2,1) = 0.5*m0r[3]; 
  data->AEM_S(2,2) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(2,3) = 0.5*m0r[1]; 
  data->AEM_S(2,5) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(3,0) = 0.5*m0r[3]; 
  data->AEM_S(3,1) = 0.5*m0r[2]; 
  data->AEM_S(3,2) = 0.5*m0r[1]; 
  data->AEM_S(3,3) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(3,4) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(3,5) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(4,0) = 0.5*m0r[4]; 
  data->AEM_S(4,1) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(4,3) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(4,4) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(5,0) = 0.5*m0r[5]; 
  data->AEM_S(5,2) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(5,3) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(5,5) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,18) = -0.5*cM[0]; 
  data->AEM_S(0,19) = -0.5*cM[1]; 
  data->AEM_S(0,20) = -0.5*cM[2]; 
  data->AEM_S(0,21) = -0.5*cM[3]; 
  data->AEM_S(0,22) = -0.5*cM[4]; 
  data->AEM_S(0,23) = -0.5*cM[5]; 
  data->AEM_S(1,18) = -0.5*cM[1]; 
  data->AEM_S(1,19) = (-0.4472135954999579*cM[4])-0.5*cM[0]; 
  data->AEM_S(1,20) = -0.5*cM[3]; 
  data->AEM_S(1,21) = -0.5*cM[2]; 
  data->AEM_S(1,22) = -0.4472135954999579*cM[1]; 
  data->AEM_S(2,18) = -0.5*cM[2]; 
  data->AEM_S(2,19) = -0.5*cM[3]; 
  data->AEM_S(2,20) = (-0.4472135954999579*cM[5])-0.5*cM[0]; 
  data->AEM_S(2,21) = -0.5*cM[1]; 
  data->AEM_S(2,23) = -0.4472135954999579*cM[2]; 
  data->AEM_S(3,18) = -0.5*cM[3]; 
  data->AEM_S(3,19) = -0.5*cM[2]; 
  data->AEM_S(3,20) = -0.5*cM[1]; 
  data->AEM_S(3,21) = (-0.4472135954999579*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  data->AEM_S(3,22) = -0.4472135954999579*cM[3]; 
  data->AEM_S(3,23) = -0.4472135954999579*cM[3]; 
  data->AEM_S(4,18) = -0.5*cM[4]; 
  data->AEM_S(4,19) = -0.4472135954999579*cM[1]; 
  data->AEM_S(4,21) = -0.4472135954999579*cM[3]; 
  data->AEM_S(4,22) = (-0.31943828249997*cM[4])-0.5*cM[0]; 
  data->AEM_S(5,18) = -0.5*cM[5]; 
  data->AEM_S(5,20) = -0.4472135954999579*cM[2]; 
  data->AEM_S(5,21) = -0.4472135954999579*cM[3]; 
  data->AEM_S(5,23) = (-0.31943828249997*cM[5])-0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(18,0) = 0.5*m1r[0]; 
  data->AEM_S(18,1) = 0.5*m1r[1]; 
  data->AEM_S(18,2) = 0.5*m1r[2]; 
  data->AEM_S(18,3) = 0.5*m1r[3]; 
  data->AEM_S(18,4) = 0.5*m1r[4]; 
  data->AEM_S(18,5) = 0.5*m1r[5]; 
  data->AEM_S(19,0) = 0.5*m1r[1]; 
  data->AEM_S(19,1) = 0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(19,2) = 0.5*m1r[3]; 
  data->AEM_S(19,3) = 0.5*m1r[2]; 
  data->AEM_S(19,4) = 0.4472135954999579*m1r[1]; 
  data->AEM_S(20,0) = 0.5*m1r[2]; 
  data->AEM_S(20,1) = 0.5*m1r[3]; 
  data->AEM_S(20,2) = 0.4472135954999579*m1r[5]+0.5*m1r[0]; 
  data->AEM_S(20,3) = 0.5*m1r[1]; 
  data->AEM_S(20,5) = 0.4472135954999579*m1r[2]; 
  data->AEM_S(21,0) = 0.5*m1r[3]; 
  data->AEM_S(21,1) = 0.5*m1r[2]; 
  data->AEM_S(21,2) = 0.5*m1r[1]; 
  data->AEM_S(21,3) = 0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(21,4) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(21,5) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(22,0) = 0.5*m1r[4]; 
  data->AEM_S(22,1) = 0.4472135954999579*m1r[1]; 
  data->AEM_S(22,3) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(22,4) = 0.31943828249997*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(23,0) = 0.5*m1r[5]; 
  data->AEM_S(23,2) = 0.4472135954999579*m1r[2]; 
  data->AEM_S(23,3) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(23,5) = 0.31943828249997*m1r[5]+0.5*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(6,6) = 0.5*m0r[0]; 
  data->AEM_S(6,7) = 0.5*m0r[1]; 
  data->AEM_S(6,8) = 0.5*m0r[2]; 
  data->AEM_S(6,9) = 0.5*m0r[3]; 
  data->AEM_S(6,10) = 0.5*m0r[4]; 
  data->AEM_S(6,11) = 0.5*m0r[5]; 
  data->AEM_S(7,6) = 0.5*m0r[1]; 
  data->AEM_S(7,7) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(7,8) = 0.5*m0r[3]; 
  data->AEM_S(7,9) = 0.5*m0r[2]; 
  data->AEM_S(7,10) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(8,6) = 0.5*m0r[2]; 
  data->AEM_S(8,7) = 0.5*m0r[3]; 
  data->AEM_S(8,8) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(8,9) = 0.5*m0r[1]; 
  data->AEM_S(8,11) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(9,6) = 0.5*m0r[3]; 
  data->AEM_S(9,7) = 0.5*m0r[2]; 
  data->AEM_S(9,8) = 0.5*m0r[1]; 
  data->AEM_S(9,9) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(9,10) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(9,11) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(10,6) = 0.5*m0r[4]; 
  data->AEM_S(10,7) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(10,9) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(10,10) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(11,6) = 0.5*m0r[5]; 
  data->AEM_S(11,8) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(11,9) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(11,11) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(6,18) = -0.5*cM[6]; 
  data->AEM_S(6,19) = -0.5*cM[7]; 
  data->AEM_S(6,20) = -0.5*cM[8]; 
  data->AEM_S(6,21) = -0.5*cM[9]; 
  data->AEM_S(6,22) = -0.5*cM[10]; 
  data->AEM_S(6,23) = -0.5*cM[11]; 
  data->AEM_S(7,18) = -0.5*cM[7]; 
  data->AEM_S(7,19) = (-0.4472135954999579*cM[10])-0.5*cM[6]; 
  data->AEM_S(7,20) = -0.5*cM[9]; 
  data->AEM_S(7,21) = -0.5*cM[8]; 
  data->AEM_S(7,22) = -0.4472135954999579*cM[7]; 
  data->AEM_S(8,18) = -0.5*cM[8]; 
  data->AEM_S(8,19) = -0.5*cM[9]; 
  data->AEM_S(8,20) = (-0.4472135954999579*cM[11])-0.5*cM[6]; 
  data->AEM_S(8,21) = -0.5*cM[7]; 
  data->AEM_S(8,23) = -0.4472135954999579*cM[8]; 
  data->AEM_S(9,18) = -0.5*cM[9]; 
  data->AEM_S(9,19) = -0.5*cM[8]; 
  data->AEM_S(9,20) = -0.5*cM[7]; 
  data->AEM_S(9,21) = (-0.4472135954999579*cM[11])-0.4472135954999579*cM[10]-0.5*cM[6]; 
  data->AEM_S(9,22) = -0.4472135954999579*cM[9]; 
  data->AEM_S(9,23) = -0.4472135954999579*cM[9]; 
  data->AEM_S(10,18) = -0.5*cM[10]; 
  data->AEM_S(10,19) = -0.4472135954999579*cM[7]; 
  data->AEM_S(10,21) = -0.4472135954999579*cM[9]; 
  data->AEM_S(10,22) = (-0.31943828249997*cM[10])-0.5*cM[6]; 
  data->AEM_S(11,18) = -0.5*cM[11]; 
  data->AEM_S(11,20) = -0.4472135954999579*cM[8]; 
  data->AEM_S(11,21) = -0.4472135954999579*cM[9]; 
  data->AEM_S(11,23) = (-0.31943828249997*cM[11])-0.5*cM[6]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(18,6) = 0.5*m1r[6]; 
  data->AEM_S(18,7) = 0.5*m1r[7]; 
  data->AEM_S(18,8) = 0.5*m1r[8]; 
  data->AEM_S(18,9) = 0.5*m1r[9]; 
  data->AEM_S(18,10) = 0.5*m1r[10]; 
  data->AEM_S(18,11) = 0.5*m1r[11]; 
  data->AEM_S(19,6) = 0.5*m1r[7]; 
  data->AEM_S(19,7) = 0.4472135954999579*m1r[10]+0.5*m1r[6]; 
  data->AEM_S(19,8) = 0.5*m1r[9]; 
  data->AEM_S(19,9) = 0.5*m1r[8]; 
  data->AEM_S(19,10) = 0.4472135954999579*m1r[7]; 
  data->AEM_S(20,6) = 0.5*m1r[8]; 
  data->AEM_S(20,7) = 0.5*m1r[9]; 
  data->AEM_S(20,8) = 0.4472135954999579*m1r[11]+0.5*m1r[6]; 
  data->AEM_S(20,9) = 0.5*m1r[7]; 
  data->AEM_S(20,11) = 0.4472135954999579*m1r[8]; 
  data->AEM_S(21,6) = 0.5*m1r[9]; 
  data->AEM_S(21,7) = 0.5*m1r[8]; 
  data->AEM_S(21,8) = 0.5*m1r[7]; 
  data->AEM_S(21,9) = 0.4472135954999579*m1r[11]+0.4472135954999579*m1r[10]+0.5*m1r[6]; 
  data->AEM_S(21,10) = 0.4472135954999579*m1r[9]; 
  data->AEM_S(21,11) = 0.4472135954999579*m1r[9]; 
  data->AEM_S(22,6) = 0.5*m1r[10]; 
  data->AEM_S(22,7) = 0.4472135954999579*m1r[7]; 
  data->AEM_S(22,9) = 0.4472135954999579*m1r[9]; 
  data->AEM_S(22,10) = 0.31943828249997*m1r[10]+0.5*m1r[6]; 
  data->AEM_S(23,6) = 0.5*m1r[11]; 
  data->AEM_S(23,8) = 0.4472135954999579*m1r[8]; 
  data->AEM_S(23,9) = 0.4472135954999579*m1r[9]; 
  data->AEM_S(23,11) = 0.31943828249997*m1r[11]+0.5*m1r[6]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(12,12) = 0.5*m0r[0]; 
  data->AEM_S(12,13) = 0.5*m0r[1]; 
  data->AEM_S(12,14) = 0.5*m0r[2]; 
  data->AEM_S(12,15) = 0.5*m0r[3]; 
  data->AEM_S(12,16) = 0.5*m0r[4]; 
  data->AEM_S(12,17) = 0.5*m0r[5]; 
  data->AEM_S(13,12) = 0.5*m0r[1]; 
  data->AEM_S(13,13) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(13,14) = 0.5*m0r[3]; 
  data->AEM_S(13,15) = 0.5*m0r[2]; 
  data->AEM_S(13,16) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(14,12) = 0.5*m0r[2]; 
  data->AEM_S(14,13) = 0.5*m0r[3]; 
  data->AEM_S(14,14) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(14,15) = 0.5*m0r[1]; 
  data->AEM_S(14,17) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(15,12) = 0.5*m0r[3]; 
  data->AEM_S(15,13) = 0.5*m0r[2]; 
  data->AEM_S(15,14) = 0.5*m0r[1]; 
  data->AEM_S(15,15) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(15,16) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(15,17) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(16,12) = 0.5*m0r[4]; 
  data->AEM_S(16,13) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(16,15) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(16,16) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(17,12) = 0.5*m0r[5]; 
  data->AEM_S(17,14) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(17,15) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(17,17) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(12,18) = -0.5*cM[12]; 
  data->AEM_S(12,19) = -0.5*cM[13]; 
  data->AEM_S(12,20) = -0.5*cM[14]; 
  data->AEM_S(12,21) = -0.5*cM[15]; 
  data->AEM_S(12,22) = -0.5*cM[16]; 
  data->AEM_S(12,23) = -0.5*cM[17]; 
  data->AEM_S(13,18) = -0.5*cM[13]; 
  data->AEM_S(13,19) = (-0.4472135954999579*cM[16])-0.5*cM[12]; 
  data->AEM_S(13,20) = -0.5*cM[15]; 
  data->AEM_S(13,21) = -0.5*cM[14]; 
  data->AEM_S(13,22) = -0.4472135954999579*cM[13]; 
  data->AEM_S(14,18) = -0.5*cM[14]; 
  data->AEM_S(14,19) = -0.5*cM[15]; 
  data->AEM_S(14,20) = (-0.4472135954999579*cM[17])-0.5*cM[12]; 
  data->AEM_S(14,21) = -0.5*cM[13]; 
  data->AEM_S(14,23) = -0.4472135954999579*cM[14]; 
  data->AEM_S(15,18) = -0.5*cM[15]; 
  data->AEM_S(15,19) = -0.5*cM[14]; 
  data->AEM_S(15,20) = -0.5*cM[13]; 
  data->AEM_S(15,21) = (-0.4472135954999579*cM[17])-0.4472135954999579*cM[16]-0.5*cM[12]; 
  data->AEM_S(15,22) = -0.4472135954999579*cM[15]; 
  data->AEM_S(15,23) = -0.4472135954999579*cM[15]; 
  data->AEM_S(16,18) = -0.5*cM[16]; 
  data->AEM_S(16,19) = -0.4472135954999579*cM[13]; 
  data->AEM_S(16,21) = -0.4472135954999579*cM[15]; 
  data->AEM_S(16,22) = (-0.31943828249997*cM[16])-0.5*cM[12]; 
  data->AEM_S(17,18) = -0.5*cM[17]; 
  data->AEM_S(17,20) = -0.4472135954999579*cM[14]; 
  data->AEM_S(17,21) = -0.4472135954999579*cM[15]; 
  data->AEM_S(17,23) = (-0.31943828249997*cM[17])-0.5*cM[12]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(18,12) = 0.5*m1r[12]; 
  data->AEM_S(18,13) = 0.5*m1r[13]; 
  data->AEM_S(18,14) = 0.5*m1r[14]; 
  data->AEM_S(18,15) = 0.5*m1r[15]; 
  data->AEM_S(18,16) = 0.5*m1r[16]; 
  data->AEM_S(18,17) = 0.5*m1r[17]; 
  data->AEM_S(19,12) = 0.5*m1r[13]; 
  data->AEM_S(19,13) = 0.4472135954999579*m1r[16]+0.5*m1r[12]; 
  data->AEM_S(19,14) = 0.5*m1r[15]; 
  data->AEM_S(19,15) = 0.5*m1r[14]; 
  data->AEM_S(19,16) = 0.4472135954999579*m1r[13]; 
  data->AEM_S(20,12) = 0.5*m1r[14]; 
  data->AEM_S(20,13) = 0.5*m1r[15]; 
  data->AEM_S(20,14) = 0.4472135954999579*m1r[17]+0.5*m1r[12]; 
  data->AEM_S(20,15) = 0.5*m1r[13]; 
  data->AEM_S(20,17) = 0.4472135954999579*m1r[14]; 
  data->AEM_S(21,12) = 0.5*m1r[15]; 
  data->AEM_S(21,13) = 0.5*m1r[14]; 
  data->AEM_S(21,14) = 0.5*m1r[13]; 
  data->AEM_S(21,15) = 0.4472135954999579*m1r[17]+0.4472135954999579*m1r[16]+0.5*m1r[12]; 
  data->AEM_S(21,16) = 0.4472135954999579*m1r[15]; 
  data->AEM_S(21,17) = 0.4472135954999579*m1r[15]; 
  data->AEM_S(22,12) = 0.5*m1r[16]; 
  data->AEM_S(22,13) = 0.4472135954999579*m1r[13]; 
  data->AEM_S(22,15) = 0.4472135954999579*m1r[15]; 
  data->AEM_S(22,16) = 0.31943828249997*m1r[16]+0.5*m1r[12]; 
  data->AEM_S(23,12) = 0.5*m1r[17]; 
  data->AEM_S(23,14) = 0.4472135954999579*m1r[14]; 
  data->AEM_S(23,15) = 0.4472135954999579*m1r[15]; 
  data->AEM_S(23,17) = 0.31943828249997*m1r[17]+0.5*m1r[12]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(18,18) = 1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(18,19) = 1.5*m0r[1]-0.5*cE[1]; 
  data->AEM_S(18,20) = 1.5*m0r[2]-0.5*cE[2]; 
  data->AEM_S(18,21) = 1.5*m0r[3]-0.5*cE[3]; 
  data->AEM_S(18,22) = 1.5*m0r[4]-0.5*cE[4]; 
  data->AEM_S(18,23) = 1.5*m0r[5]-0.5*cE[5]; 
  data->AEM_S(19,18) = 1.5*m0r[1]-0.5*cE[1]; 
  data->AEM_S(19,19) = 1.341640786499874*m0r[4]-0.4472135954999579*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(19,20) = 1.5*m0r[3]-0.5*cE[3]; 
  data->AEM_S(19,21) = 1.5*m0r[2]-0.5*cE[2]; 
  data->AEM_S(19,22) = 1.341640786499874*m0r[1]-0.4472135954999579*cE[1]; 
  data->AEM_S(20,18) = 1.5*m0r[2]-0.5*cE[2]; 
  data->AEM_S(20,19) = 1.5*m0r[3]-0.5*cE[3]; 
  data->AEM_S(20,20) = 1.341640786499874*m0r[5]-0.4472135954999579*cE[5]+1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(20,21) = 1.5*m0r[1]-0.5*cE[1]; 
  data->AEM_S(20,23) = 1.341640786499874*m0r[2]-0.4472135954999579*cE[2]; 
  data->AEM_S(21,18) = 1.5*m0r[3]-0.5*cE[3]; 
  data->AEM_S(21,19) = 1.5*m0r[2]-0.5*cE[2]; 
  data->AEM_S(21,20) = 1.5*m0r[1]-0.5*cE[1]; 
  data->AEM_S(21,21) = 1.341640786499874*m0r[5]-0.4472135954999579*cE[5]+1.341640786499874*m0r[4]-0.4472135954999579*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(21,22) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  data->AEM_S(21,23) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  data->AEM_S(22,18) = 1.5*m0r[4]-0.5*cE[4]; 
  data->AEM_S(22,19) = 1.341640786499874*m0r[1]-0.4472135954999579*cE[1]; 
  data->AEM_S(22,21) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  data->AEM_S(22,22) = 0.9583148474999099*m0r[4]-0.31943828249997*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(23,18) = 1.5*m0r[5]-0.5*cE[5]; 
  data->AEM_S(23,20) = 1.341640786499874*m0r[2]-0.4472135954999579*cE[2]; 
  data->AEM_S(23,21) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  data->AEM_S(23,23) = 0.9583148474999099*m0r[5]-0.31943828249997*cE[5]+1.5*m0r[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<6,12>(0,6).setZero(); 
  data->AEM_S.block<12,6>(6,0).setZero(); 
  data->AEM_S.block<6,6>(6,12).setZero(); 
  data->AEM_S.block<6,6>(12,6).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m1r[12],m1r[13],m1r[14],m1r[15],m1r[16],m1r[17],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,18,1) = data->u_S.segment<18>(0); 
 
  Eigen::Map<VectorXd>(vtSq,6,1) = data->u_S.segment<6>(18); 
 
} 
 
