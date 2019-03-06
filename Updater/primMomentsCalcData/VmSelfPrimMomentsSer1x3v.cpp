
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmM0Star1x3vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[2]*dxvl[3]*(wr[1]-wl[1]); 
 
  out[0] += ((-0.8164965809277261*fr[2])+0.8164965809277261*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[5])+0.8164965809277261*fl[5]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
 
} 
 
void VmM0Star1x3vSer_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[1]*dxvl[3]*(wr[2]-wl[2]); 
 
  out[0] += ((-0.8164965809277261*fr[3])+0.8164965809277261*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[6])+0.8164965809277261*fl[6]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
 
} 
 
void VmM0Star1x3vSer_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[1]*dxvl[2]*(wr[3]-wl[3]); 
 
  out[0] += ((-0.8164965809277261*fr[4])+0.8164965809277261*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[8])+0.8164965809277261*fl[8]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
 
} 
 
void VmM1iM2Star1x3vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2) 
{ 
  // w[4]:    Cell-center coordinates. 
  // dxv[4]:  Cell length in each direciton. 
  // intFac:  for gyrokinetics (not used in Vm). 
  // m_:      mass (not used in Vm). 
  // Bmag[2]: Magnetic field magnitude (not used in Vm). 
  // f:       Distribution function. 
  // outM1i:  Contribution to M_1^star from this cell. 
  // outM2:   Contribution to M_2^star from this cell. 
 
  const double volFact = 0.125*dxv[1]*dxv[2]*dxv[3]; 
  double wvSq[3]; 
  wvSq[0]  = w[1]*w[1]; 
  wvSq[1]  = w[2]*w[2]; 
  wvSq[2]  = w[3]*w[3]; 
  double dvSq[3]; 
  dvSq[0] = dxv[1]*dxv[1]; 
  dvSq[1] = dxv[2]*dxv[2]; 
  dvSq[2] = dxv[3]*dxv[3]; 
  double tempM0[2]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 

  outM1i[0] += tempM0[0]*w[1]; 
  outM1i[1] += tempM0[1]*w[1]; 
  outM1i[2] += tempM0[0]*w[2]; 
  outM1i[3] += tempM0[1]*w[2]; 
  outM1i[4] += tempM0[0]*w[3]; 
  outM1i[5] += tempM0[1]*w[3]; 

  outM2[0] += (0.8164965809277261*dxv[3]*w[3]*f[4]+0.8164965809277261*dxv[2]*w[2]*f[3]+0.8164965809277261*dxv[1]*w[1]*f[2])*volFact+tempM0[0]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[1] += (0.8164965809277261*dxv[3]*w[3]*f[8]+0.8164965809277261*dxv[2]*w[2]*f[6]+0.8164965809277261*dxv[1]*w[1]*f[5])*volFact+tempM0[1]*(wvSq[2]+wvSq[1]+wvSq[0]); 
 
} 
void VmBoundaryIntegral1x3vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[16]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[48]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[12])+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[80]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (3.741657386773942*fIn[32]-3.16227766016838*fIn[12]+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += (3.741657386773942*fIn[49]-3.16227766016838*fIn[20]+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS; 
    out[3] += (2.449489742783178*fIn[48]-1.414213562373095*fIn[31])*dS; 
 
  } else {
 
    out[0] += (3.741657386773942*fIn[32]+3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += (3.741657386773942*fIn[49]+3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS; 
    out[3] += (2.449489742783178*fIn[48]+1.414213562373095*fIn[31])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[16]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[1]*fIn[2]-0.7071067811865475*fIn[0]*dxv[1])*dS; 
    out[1] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[1]*fIn[5]-0.7071067811865475*dxv[1]*fIn[1])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[1]*fIn[2])-0.7071067811865475*fIn[0]*dxv[1])*dS; 
    out[1] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[1]*fIn[5])-0.7071067811865475*dxv[1]*fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[48]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[12])+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[80]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (3.741657386773942*fIn[32]-3.16227766016838*fIn[12]+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.741657386773942*fIn[49]-3.16227766016838*fIn[20]+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[48]-1.414213562373095*fIn[31])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.741657386773942*fIn[32]+3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.741657386773942*fIn[49]+3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[48]+1.414213562373095*fIn[31])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[16]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[2] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[3] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
 
  } else {
 
    out[2] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[3] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[48]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[3] += ((-3.16227766016838*fIn[13])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[4] += ((-3.16227766016838*fIn[23])+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
    out[5] += (2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS; 
 
  } else {
 
    out[3] += (3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[4] += (3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
    out[5] += (2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[80]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[4] += (3.741657386773942*fIn[33]-3.16227766016838*fIn[13]+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[5] += (3.741657386773942*fIn[52]-3.16227766016838*fIn[23]+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
    out[6] += (2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS; 
    out[7] += (2.449489742783178*fIn[50]-1.414213562373095*fIn[31])*dS; 
 
  } else {
 
    out[4] += (3.741657386773942*fIn[33]+3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[5] += (3.741657386773942*fIn[52]+3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
    out[6] += (2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS; 
    out[7] += (2.449489742783178*fIn[50]+1.414213562373095*fIn[31])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[16]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[3]-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[6]-0.7071067811865475*fIn[1]*dxv[2])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[3])-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[6])-0.7071067811865475*fIn[1]*dxv[2])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[48]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[13])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[23])+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[80]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (3.741657386773942*fIn[33]-3.16227766016838*fIn[13]+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.741657386773942*fIn[52]-3.16227766016838*fIn[23]+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[50]-1.414213562373095*fIn[31])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.741657386773942*fIn[33]+3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.741657386773942*fIn[52]+3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[50]+1.414213562373095*fIn[31])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[16]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[4] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
 
  } else {
 
    out[4] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[48]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[6] += ((-3.16227766016838*fIn[14])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[7] += ((-3.16227766016838*fIn[28])+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
    out[8] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS; 
 
  } else {
 
    out[6] += (3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[7] += (3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
    out[8] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_F_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[80]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[8] += (3.741657386773942*fIn[34]-3.16227766016838*fIn[14]+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[9] += (3.741657386773942*fIn[57]-3.16227766016838*fIn[28]+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS; 
    out[11] += (2.449489742783178*fIn[54]-1.414213562373095*fIn[31])*dS; 
 
  } else {
 
    out[8] += (3.741657386773942*fIn[34]+3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[9] += (3.741657386773942*fIn[57]+3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS; 
    out[11] += (2.449489742783178*fIn[54]+1.414213562373095*fIn[31])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[16]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[4]-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[8]-0.7071067811865475*fIn[1]*dxv[3])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[4])-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[8])-0.7071067811865475*fIn[1]*dxv[3])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[48]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[14])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[28])+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vSer_vF_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[80]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (3.741657386773942*fIn[34]-3.16227766016838*fIn[14]+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.741657386773942*fIn[57]-3.16227766016838*fIn[28]+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[54]-1.414213562373095*fIn[31])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.741657386773942*fIn[34]+3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.741657386773942*fIn[57]+3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += (2.449489742783178*fIn[54]+1.414213562373095*fIn[31])*dS*vBoundary; 
 
  }
 
} 
 
void VmSelfPrimMoments1x3vSer_P1(binOpData_t* data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1:       moments of the distribution function. 
  // m0S,m1S,m1S: star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7071067811865475*m0[0]-1.224744871391589*m0[1] < 0) { 
    cellAvg = true;
  }
  if (1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[2]; 
  double m1r[6]; 
  double m0Sr[2]; 
  double m1Sr[6]; 
  double m2Sr[2]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = m1[2]; 
    m1r[3] = 0.0; 
    m1r[4] = m1[4]; 
    m1r[5] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = 0.0; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m1r[4] = m1[4]; 
    m1r[5] = m1[5]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = m1S[5]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(8,8); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,0) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,1) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,6) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,7) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,6) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,7) = -0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(6,0) = 0.7071067811865475*m1Sr[0]; 
  data->AEM_S(6,1) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(7,0) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(7,1) = 0.7071067811865475*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(2,2) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(2,3) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(3,2) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(3,3) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(2,6) = -0.7071067811865475*cM[2]; 
  data->AEM_S(2,7) = -0.7071067811865475*cM[3]; 
  data->AEM_S(3,6) = -0.7071067811865475*cM[3]; 
  data->AEM_S(3,7) = -0.7071067811865475*cM[2]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(6,2) = 0.7071067811865475*m1Sr[2]; 
  data->AEM_S(6,3) = 0.7071067811865475*m1Sr[3]; 
  data->AEM_S(7,2) = 0.7071067811865475*m1Sr[3]; 
  data->AEM_S(7,3) = 0.7071067811865475*m1Sr[2]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(4,4) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(4,5) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(5,4) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(5,5) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(4,6) = -0.7071067811865475*cM[4]; 
  data->AEM_S(4,7) = -0.7071067811865475*cM[5]; 
  data->AEM_S(5,6) = -0.7071067811865475*cM[5]; 
  data->AEM_S(5,7) = -0.7071067811865475*cM[4]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(6,4) = 0.7071067811865475*m1Sr[4]; 
  data->AEM_S(6,5) = 0.7071067811865475*m1Sr[5]; 
  data->AEM_S(7,4) = 0.7071067811865475*m1Sr[5]; 
  data->AEM_S(7,5) = 0.7071067811865475*m1Sr[4]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(6,6) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(6,7) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(7,6) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(7,7) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<2,4>(0,2).setZero(); 
  data->AEM_S.block<4,2>(2,0).setZero(); 
  data->AEM_S.block<2,2>(2,4).setZero(); 
  data->AEM_S.block<2,2>(4,2).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m2Sr[0],m2Sr[1]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,6,1) = data->u_S.segment<6>(0); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = data->u_S.segment<2>(6); 
 
} 
 
void VmSelfPrimMoments1x3vSer_P2(binOpData_t* data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (1.58113883008419*m0[2]-1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.58113883008419*m0[2]+1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[3]; 
  double m1r[9]; 
  double m2r[3]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(12,12); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(0,2) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(1,0) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,1) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(1,2) = 0.6324555320336759*m0r[1]; 
  data->AEM_S(2,0) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(2,1) = 0.6324555320336759*m0r[1]; 
  data->AEM_S(2,2) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,9) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,10) = -0.7071067811865475*cM[1]; 
  data->AEM_S(0,11) = -0.7071067811865475*cM[2]; 
  data->AEM_S(1,9) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,10) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(1,11) = -0.6324555320336759*cM[1]; 
  data->AEM_S(2,9) = -0.7071067811865475*cM[2]; 
  data->AEM_S(2,10) = -0.6324555320336759*cM[1]; 
  data->AEM_S(2,11) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(9,0) = 0.7071067811865475*m1r[0]; 
  data->AEM_S(9,1) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(9,2) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(10,0) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(10,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(10,2) = 0.6324555320336759*m1r[1]; 
  data->AEM_S(11,0) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(11,1) = 0.6324555320336759*m1r[1]; 
  data->AEM_S(11,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(3,3) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(3,4) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(3,5) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(4,3) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(4,4) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(4,5) = 0.6324555320336759*m0r[1]; 
  data->AEM_S(5,3) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(5,4) = 0.6324555320336759*m0r[1]; 
  data->AEM_S(5,5) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(3,9) = -0.7071067811865475*cM[3]; 
  data->AEM_S(3,10) = -0.7071067811865475*cM[4]; 
  data->AEM_S(3,11) = -0.7071067811865475*cM[5]; 
  data->AEM_S(4,9) = -0.7071067811865475*cM[4]; 
  data->AEM_S(4,10) = (-0.6324555320336759*cM[5])-0.7071067811865475*cM[3]; 
  data->AEM_S(4,11) = -0.6324555320336759*cM[4]; 
  data->AEM_S(5,9) = -0.7071067811865475*cM[5]; 
  data->AEM_S(5,10) = -0.6324555320336759*cM[4]; 
  data->AEM_S(5,11) = (-0.4517539514526256*cM[5])-0.7071067811865475*cM[3]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(9,3) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(9,4) = 0.7071067811865475*m1r[4]; 
  data->AEM_S(9,5) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(10,3) = 0.7071067811865475*m1r[4]; 
  data->AEM_S(10,4) = 0.6324555320336759*m1r[5]+0.7071067811865475*m1r[3]; 
  data->AEM_S(10,5) = 0.6324555320336759*m1r[4]; 
  data->AEM_S(11,3) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(11,4) = 0.6324555320336759*m1r[4]; 
  data->AEM_S(11,5) = 0.4517539514526256*m1r[5]+0.7071067811865475*m1r[3]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(6,6) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(6,7) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(6,8) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(7,6) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(7,7) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(7,8) = 0.6324555320336759*m0r[1]; 
  data->AEM_S(8,6) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(8,7) = 0.6324555320336759*m0r[1]; 
  data->AEM_S(8,8) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(6,9) = -0.7071067811865475*cM[6]; 
  data->AEM_S(6,10) = -0.7071067811865475*cM[7]; 
  data->AEM_S(6,11) = -0.7071067811865475*cM[8]; 
  data->AEM_S(7,9) = -0.7071067811865475*cM[7]; 
  data->AEM_S(7,10) = (-0.6324555320336759*cM[8])-0.7071067811865475*cM[6]; 
  data->AEM_S(7,11) = -0.6324555320336759*cM[7]; 
  data->AEM_S(8,9) = -0.7071067811865475*cM[8]; 
  data->AEM_S(8,10) = -0.6324555320336759*cM[7]; 
  data->AEM_S(8,11) = (-0.4517539514526256*cM[8])-0.7071067811865475*cM[6]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(9,6) = 0.7071067811865475*m1r[6]; 
  data->AEM_S(9,7) = 0.7071067811865475*m1r[7]; 
  data->AEM_S(9,8) = 0.7071067811865475*m1r[8]; 
  data->AEM_S(10,6) = 0.7071067811865475*m1r[7]; 
  data->AEM_S(10,7) = 0.6324555320336759*m1r[8]+0.7071067811865475*m1r[6]; 
  data->AEM_S(10,8) = 0.6324555320336759*m1r[7]; 
  data->AEM_S(11,6) = 0.7071067811865475*m1r[8]; 
  data->AEM_S(11,7) = 0.6324555320336759*m1r[7]; 
  data->AEM_S(11,8) = 0.4517539514526256*m1r[8]+0.7071067811865475*m1r[6]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(9,9) = 2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(9,10) = 2.121320343559642*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(9,11) = 2.121320343559642*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(10,9) = 2.121320343559642*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(10,10) = 1.897366596101028*m0r[2]-0.6324555320336759*cE[2]+2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(10,11) = 1.897366596101028*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(11,9) = 2.121320343559642*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(11,10) = 1.897366596101028*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(11,11) = 1.355261854357877*m0r[2]-0.4517539514526256*cE[2]+2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<3,6>(0,3).setZero(); 
  data->AEM_S.block<6,3>(3,0).setZero(); 
  data->AEM_S.block<3,3>(3,6).setZero(); 
  data->AEM_S.block<3,3>(6,3).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m2r[0],m2r[1],m2r[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,9,1) = data->u_S.segment<9>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = data->u_S.segment<3>(9); 
 
} 
 
void VmSelfPrimMoments1x3vSer_P3(binOpData_t* data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-1.870828693386971*m0[3])+1.58113883008419*m0[2]-1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.870828693386971*m0[3]+1.58113883008419*m0[2]+1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[4]; 
  double m1r[12]; 
  double m2r[4]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(16,16); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(0,2) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(0,3) = 0.7071067811865475*m0r[3]; 
  data->AEM_S(1,0) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,1) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(1,2) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  data->AEM_S(1,3) = 0.6210590034081186*m0r[2]; 
  data->AEM_S(2,0) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(2,1) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  data->AEM_S(2,2) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(2,3) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  data->AEM_S(3,0) = 0.7071067811865475*m0r[3]; 
  data->AEM_S(3,1) = 0.6210590034081186*m0r[2]; 
  data->AEM_S(3,2) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  data->AEM_S(3,3) = 0.421637021355784*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,12) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,13) = -0.7071067811865475*cM[1]; 
  data->AEM_S(0,14) = -0.7071067811865475*cM[2]; 
  data->AEM_S(0,15) = -0.7071067811865475*cM[3]; 
  data->AEM_S(1,12) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,13) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(1,14) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  data->AEM_S(1,15) = -0.6210590034081186*cM[2]; 
  data->AEM_S(2,12) = -0.7071067811865475*cM[2]; 
  data->AEM_S(2,13) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  data->AEM_S(2,14) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(2,15) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  data->AEM_S(3,12) = -0.7071067811865475*cM[3]; 
  data->AEM_S(3,13) = -0.6210590034081186*cM[2]; 
  data->AEM_S(3,14) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  data->AEM_S(3,15) = (-0.421637021355784*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(12,0) = 0.7071067811865475*m1r[0]; 
  data->AEM_S(12,1) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(12,2) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(12,3) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(13,0) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(13,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(13,2) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  data->AEM_S(13,3) = 0.6210590034081186*m1r[2]; 
  data->AEM_S(14,0) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(14,1) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  data->AEM_S(14,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(14,3) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  data->AEM_S(15,0) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(15,1) = 0.6210590034081186*m1r[2]; 
  data->AEM_S(15,2) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  data->AEM_S(15,3) = 0.421637021355784*m1r[2]+0.7071067811865475*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(4,4) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(4,5) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(4,6) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(4,7) = 0.7071067811865475*m0r[3]; 
  data->AEM_S(5,4) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(5,5) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(5,6) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  data->AEM_S(5,7) = 0.6210590034081186*m0r[2]; 
  data->AEM_S(6,4) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(6,5) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  data->AEM_S(6,6) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(6,7) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  data->AEM_S(7,4) = 0.7071067811865475*m0r[3]; 
  data->AEM_S(7,5) = 0.6210590034081186*m0r[2]; 
  data->AEM_S(7,6) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  data->AEM_S(7,7) = 0.421637021355784*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(4,12) = -0.7071067811865475*cM[4]; 
  data->AEM_S(4,13) = -0.7071067811865475*cM[5]; 
  data->AEM_S(4,14) = -0.7071067811865475*cM[6]; 
  data->AEM_S(4,15) = -0.7071067811865475*cM[7]; 
  data->AEM_S(5,12) = -0.7071067811865475*cM[5]; 
  data->AEM_S(5,13) = (-0.6324555320336759*cM[6])-0.7071067811865475*cM[4]; 
  data->AEM_S(5,14) = (-0.6210590034081186*cM[7])-0.6324555320336759*cM[5]; 
  data->AEM_S(5,15) = -0.6210590034081186*cM[6]; 
  data->AEM_S(6,12) = -0.7071067811865475*cM[6]; 
  data->AEM_S(6,13) = (-0.6210590034081186*cM[7])-0.6324555320336759*cM[5]; 
  data->AEM_S(6,14) = (-0.4517539514526256*cM[6])-0.7071067811865475*cM[4]; 
  data->AEM_S(6,15) = (-0.421637021355784*cM[7])-0.6210590034081186*cM[5]; 
  data->AEM_S(7,12) = -0.7071067811865475*cM[7]; 
  data->AEM_S(7,13) = -0.6210590034081186*cM[6]; 
  data->AEM_S(7,14) = (-0.421637021355784*cM[7])-0.6210590034081186*cM[5]; 
  data->AEM_S(7,15) = (-0.421637021355784*cM[6])-0.7071067811865475*cM[4]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(12,4) = 0.7071067811865475*m1r[4]; 
  data->AEM_S(12,5) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(12,6) = 0.7071067811865475*m1r[6]; 
  data->AEM_S(12,7) = 0.7071067811865475*m1r[7]; 
  data->AEM_S(13,4) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(13,5) = 0.6324555320336759*m1r[6]+0.7071067811865475*m1r[4]; 
  data->AEM_S(13,6) = 0.6210590034081186*m1r[7]+0.6324555320336759*m1r[5]; 
  data->AEM_S(13,7) = 0.6210590034081186*m1r[6]; 
  data->AEM_S(14,4) = 0.7071067811865475*m1r[6]; 
  data->AEM_S(14,5) = 0.6210590034081186*m1r[7]+0.6324555320336759*m1r[5]; 
  data->AEM_S(14,6) = 0.4517539514526256*m1r[6]+0.7071067811865475*m1r[4]; 
  data->AEM_S(14,7) = 0.421637021355784*m1r[7]+0.6210590034081186*m1r[5]; 
  data->AEM_S(15,4) = 0.7071067811865475*m1r[7]; 
  data->AEM_S(15,5) = 0.6210590034081186*m1r[6]; 
  data->AEM_S(15,6) = 0.421637021355784*m1r[7]+0.6210590034081186*m1r[5]; 
  data->AEM_S(15,7) = 0.421637021355784*m1r[6]+0.7071067811865475*m1r[4]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(8,8) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(8,9) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(8,10) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(8,11) = 0.7071067811865475*m0r[3]; 
  data->AEM_S(9,8) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(9,9) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(9,10) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  data->AEM_S(9,11) = 0.6210590034081186*m0r[2]; 
  data->AEM_S(10,8) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(10,9) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  data->AEM_S(10,10) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(10,11) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  data->AEM_S(11,8) = 0.7071067811865475*m0r[3]; 
  data->AEM_S(11,9) = 0.6210590034081186*m0r[2]; 
  data->AEM_S(11,10) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  data->AEM_S(11,11) = 0.421637021355784*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(8,12) = -0.7071067811865475*cM[8]; 
  data->AEM_S(8,13) = -0.7071067811865475*cM[9]; 
  data->AEM_S(8,14) = -0.7071067811865475*cM[10]; 
  data->AEM_S(8,15) = -0.7071067811865475*cM[11]; 
  data->AEM_S(9,12) = -0.7071067811865475*cM[9]; 
  data->AEM_S(9,13) = (-0.6324555320336759*cM[10])-0.7071067811865475*cM[8]; 
  data->AEM_S(9,14) = (-0.6210590034081186*cM[11])-0.6324555320336759*cM[9]; 
  data->AEM_S(9,15) = -0.6210590034081186*cM[10]; 
  data->AEM_S(10,12) = -0.7071067811865475*cM[10]; 
  data->AEM_S(10,13) = (-0.6210590034081186*cM[11])-0.6324555320336759*cM[9]; 
  data->AEM_S(10,14) = (-0.4517539514526256*cM[10])-0.7071067811865475*cM[8]; 
  data->AEM_S(10,15) = (-0.421637021355784*cM[11])-0.6210590034081186*cM[9]; 
  data->AEM_S(11,12) = -0.7071067811865475*cM[11]; 
  data->AEM_S(11,13) = -0.6210590034081186*cM[10]; 
  data->AEM_S(11,14) = (-0.421637021355784*cM[11])-0.6210590034081186*cM[9]; 
  data->AEM_S(11,15) = (-0.421637021355784*cM[10])-0.7071067811865475*cM[8]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(12,8) = 0.7071067811865475*m1r[8]; 
  data->AEM_S(12,9) = 0.7071067811865475*m1r[9]; 
  data->AEM_S(12,10) = 0.7071067811865475*m1r[10]; 
  data->AEM_S(12,11) = 0.7071067811865475*m1r[11]; 
  data->AEM_S(13,8) = 0.7071067811865475*m1r[9]; 
  data->AEM_S(13,9) = 0.6324555320336759*m1r[10]+0.7071067811865475*m1r[8]; 
  data->AEM_S(13,10) = 0.6210590034081186*m1r[11]+0.6324555320336759*m1r[9]; 
  data->AEM_S(13,11) = 0.6210590034081186*m1r[10]; 
  data->AEM_S(14,8) = 0.7071067811865475*m1r[10]; 
  data->AEM_S(14,9) = 0.6210590034081186*m1r[11]+0.6324555320336759*m1r[9]; 
  data->AEM_S(14,10) = 0.4517539514526256*m1r[10]+0.7071067811865475*m1r[8]; 
  data->AEM_S(14,11) = 0.421637021355784*m1r[11]+0.6210590034081186*m1r[9]; 
  data->AEM_S(15,8) = 0.7071067811865475*m1r[11]; 
  data->AEM_S(15,9) = 0.6210590034081186*m1r[10]; 
  data->AEM_S(15,10) = 0.421637021355784*m1r[11]+0.6210590034081186*m1r[9]; 
  data->AEM_S(15,11) = 0.421637021355784*m1r[10]+0.7071067811865475*m1r[8]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(12,12) = 2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(12,13) = 2.121320343559642*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(12,14) = 2.121320343559642*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(12,15) = 2.121320343559642*m0r[3]-0.7071067811865475*cE[3]; 
  data->AEM_S(13,12) = 2.121320343559642*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(13,13) = 1.897366596101028*m0r[2]-0.6324555320336759*cE[2]+2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(13,14) = 1.863177010224355*m0r[3]-0.6210590034081186*cE[3]+1.897366596101028*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(13,15) = 1.863177010224355*m0r[2]-0.6210590034081186*cE[2]; 
  data->AEM_S(14,12) = 2.121320343559642*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(14,13) = 1.863177010224355*m0r[3]-0.6210590034081186*cE[3]+1.897366596101028*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(14,14) = 1.355261854357877*m0r[2]-0.4517539514526256*cE[2]+2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(14,15) = 1.264911064067352*m0r[3]-0.421637021355784*cE[3]+1.863177010224355*m0r[1]-0.6210590034081186*cE[1]; 
  data->AEM_S(15,12) = 2.121320343559642*m0r[3]-0.7071067811865475*cE[3]; 
  data->AEM_S(15,13) = 1.863177010224355*m0r[2]-0.6210590034081186*cE[2]; 
  data->AEM_S(15,14) = 1.264911064067352*m0r[3]-0.421637021355784*cE[3]+1.863177010224355*m0r[1]-0.6210590034081186*cE[1]; 
  data->AEM_S(15,15) = 1.264911064067352*m0r[2]-0.421637021355784*cE[2]+2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<4,8>(0,4).setZero(); 
  data->AEM_S.block<8,4>(4,0).setZero(); 
  data->AEM_S.block<4,4>(4,8).setZero(); 
  data->AEM_S.block<4,4>(8,4).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m2r[0],m2r[1],m2r[2],m2r[3]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,12,1) = data->u_S.segment<12>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = data->u_S.segment<4>(12); 
 
} 
 
