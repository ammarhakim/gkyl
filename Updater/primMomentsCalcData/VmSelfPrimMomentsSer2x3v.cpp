
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmM0Star2x3vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[5]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[3]*dxvl[4]*(wr[2]-wl[2]); 
 
  out[0] += ((-0.8164965809277261*fr[3])+0.8164965809277261*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[7])+0.8164965809277261*fl[7]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += ((-0.8164965809277261*fr[8])+0.8164965809277261*fl[8]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
  out[3] += ((-0.8164965809277261*fr[16])+0.8164965809277261*fl[16]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6])*dS; 
 
} 
 
void VmM0Star2x3vSer_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[5]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[2]*dxvl[4]*(wr[3]-wl[3]); 
 
  out[0] += ((-0.8164965809277261*fr[4])+0.8164965809277261*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[9])+0.8164965809277261*fl[9]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += ((-0.8164965809277261*fr[10])+0.8164965809277261*fl[10]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
  out[3] += ((-0.8164965809277261*fr[17])+0.8164965809277261*fl[17]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6])*dS; 
 
} 
 
void VmM0Star2x3vSer_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[5]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[2]*dxvl[3]*(wr[4]-wl[4]); 
 
  out[0] += ((-0.8164965809277261*fr[5])+0.8164965809277261*fl[5]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[12])+0.8164965809277261*fl[12]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
  out[2] += ((-0.8164965809277261*fr[13])+0.8164965809277261*fl[13]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2])*dS; 
  out[3] += ((-0.8164965809277261*fr[20])+0.8164965809277261*fl[20]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6])*dS; 
 
} 
 
void VmM1iM2Star2x3vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2) 
{ 
  // w[5]:    Cell-center coordinates. 
  // dxv[5]:  Cell length in each direciton. 
  // intFac:  for gyrokinetics (not used in Vm). 
  // m_:      mass (not used in Vm). 
  // Bmag[4]: Magnetic field magnitude (not used in Vm). 
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
  double tempM0[4]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 
  tempM0[2] = 2.828427124746191*f[2]*volFact; 
  tempM0[3] = 2.828427124746191*f[6]*volFact; 

  outM1i[0] += tempM0[0]*w[2]; 
  outM1i[1] += tempM0[1]*w[2]; 
  outM1i[2] += tempM0[2]*w[2]; 
  outM1i[3] += w[2]*tempM0[3]; 
  outM1i[4] += tempM0[0]*w[3]; 
  outM1i[5] += tempM0[1]*w[3]; 
  outM1i[6] += tempM0[2]*w[3]; 
  outM1i[7] += tempM0[3]*w[3]; 
  outM1i[8] += tempM0[0]*w[4]; 
  outM1i[9] += tempM0[1]*w[4]; 
  outM1i[10] += tempM0[2]*w[4]; 
  outM1i[11] += tempM0[3]*w[4]; 

  outM2[0] += (0.8164965809277261*dxv[4]*w[4]*f[5]+0.8164965809277261*dxv[3]*w[3]*f[4]+0.8164965809277261*dxv[2]*w[2]*f[3])*volFact+tempM0[0]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[1] += (0.8164965809277261*dxv[4]*w[4]*f[12]+0.8164965809277261*dxv[3]*w[3]*f[9]+0.8164965809277261*dxv[2]*w[2]*f[7])*volFact+tempM0[1]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[2] += (0.8164965809277261*dxv[4]*w[4]*f[13]+0.8164965809277261*dxv[3]*w[3]*f[10]+0.8164965809277261*dxv[2]*w[2]*f[8])*volFact+tempM0[2]*(wvSq[2]+wvSq[1]+wvSq[0]); 
  outM2[3] += (0.8164965809277261*dxv[4]*w[4]*f[20]+0.8164965809277261*dxv[3]*w[3]*f[17]+0.8164965809277261*dxv[2]*w[2]*f[16])*volFact+tempM0[3]*(wvSq[2]+wvSq[1]+wvSq[0]); 
 
} 
void VmBoundaryIntegral2x3vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[7]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[16]-1.414213562373095*fIn[6])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[7]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[16]+1.414213562373095*fIn[6])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[112]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[18])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[1] += ((-3.16227766016838*fIn[35])+2.449489742783178*fIn[7]-1.414213562373095*fIn[1])*dS; 
    out[2] += ((-3.16227766016838*fIn[36])+2.449489742783178*fIn[8]-1.414213562373095*fIn[2])*dS; 
    out[3] += ((-3.16227766016838*fIn[58])+2.449489742783178*fIn[21]-1.414213562373095*fIn[6])*dS; 
    out[4] += (2.449489742783178*fIn[33]-1.414213562373095*fIn[16])*dS; 
    out[5] += (2.449489742783178*fIn[34]-1.414213562373095*fIn[17])*dS; 
    out[6] += (2.449489742783178*fIn[56]-1.414213562373095*fIn[31])*dS; 
    out[7] += (2.449489742783178*fIn[57]-1.414213562373095*fIn[32])*dS; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[18]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[1] += (3.16227766016838*fIn[35]+2.449489742783178*fIn[7]+1.414213562373095*fIn[1])*dS; 
    out[2] += (3.16227766016838*fIn[36]+2.449489742783178*fIn[8]+1.414213562373095*fIn[2])*dS; 
    out[3] += (3.16227766016838*fIn[58]+2.449489742783178*fIn[21]+1.414213562373095*fIn[6])*dS; 
    out[4] += (2.449489742783178*fIn[33]+1.414213562373095*fIn[16])*dS; 
    out[5] += (2.449489742783178*fIn[34]+1.414213562373095*fIn[17])*dS; 
    out[6] += (2.449489742783178*fIn[56]+1.414213562373095*fIn[31])*dS; 
    out[7] += (2.449489742783178*fIn[57]+1.414213562373095*fIn[32])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[3]-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += (2.449489742783178*fIn[7]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[7]-0.7071067811865475*fIn[1]*dxv[2])*dS; 
    out[2] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[2])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[8]-0.7071067811865475*dxv[2]*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[16]-1.414213562373095*fIn[6])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[16]-0.7071067811865475*dxv[2]*fIn[6])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[3])-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += (2.449489742783178*fIn[7]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[7])-0.7071067811865475*fIn[1]*dxv[2])*dS; 
    out[2] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[2])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[8])-0.7071067811865475*dxv[2]*fIn[2])*dS; 
    out[3] += (2.449489742783178*fIn[16]+1.414213562373095*fIn[6])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[16])-0.7071067811865475*dxv[2]*fIn[6])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[112]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[18])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[35])+2.449489742783178*fIn[7]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += ((-3.16227766016838*fIn[36])+2.449489742783178*fIn[8]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += ((-3.16227766016838*fIn[58])+2.449489742783178*fIn[21]-1.414213562373095*fIn[6])*dS*vBoundary; 
    out[4] += (2.449489742783178*fIn[33]-1.414213562373095*fIn[16])*dS*vBoundary; 
    out[5] += (2.449489742783178*fIn[34]-1.414213562373095*fIn[17])*dS*vBoundary; 
    out[6] += (2.449489742783178*fIn[56]-1.414213562373095*fIn[31])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[57]-1.414213562373095*fIn[32])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[18]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[35]+2.449489742783178*fIn[7]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (3.16227766016838*fIn[36]+2.449489742783178*fIn[8]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (3.16227766016838*fIn[58]+2.449489742783178*fIn[21]+1.414213562373095*fIn[6])*dS*vBoundary; 
    out[4] += (2.449489742783178*fIn[33]+1.414213562373095*fIn[16])*dS*vBoundary; 
    out[5] += (2.449489742783178*fIn[34]+1.414213562373095*fIn[17])*dS*vBoundary; 
    out[6] += (2.449489742783178*fIn[56]+1.414213562373095*fIn[31])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[57]+1.414213562373095*fIn[32])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  if (atLower) {
 
    out[4] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[9]-1.414213562373095*fIn[1])*dS; 
    out[6] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[2])*dS; 
    out[7] += (2.449489742783178*fIn[17]-1.414213562373095*fIn[6])*dS; 
 
  } else {
 
    out[4] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[5] += (2.449489742783178*fIn[9]+1.414213562373095*fIn[1])*dS; 
    out[6] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[2])*dS; 
    out[7] += (2.449489742783178*fIn[17]+1.414213562373095*fIn[6])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[112]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  if (atLower) {
 
    out[8] += ((-3.16227766016838*fIn[19])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[9] += ((-3.16227766016838*fIn[40])+2.449489742783178*fIn[9]-1.414213562373095*fIn[1])*dS; 
    out[10] += ((-3.16227766016838*fIn[41])+2.449489742783178*fIn[10]-1.414213562373095*fIn[2])*dS; 
    out[11] += ((-3.16227766016838*fIn[65])+2.449489742783178*fIn[22]-1.414213562373095*fIn[6])*dS; 
    out[12] += (2.449489742783178*fIn[37]-1.414213562373095*fIn[16])*dS; 
    out[13] += (2.449489742783178*fIn[38]-1.414213562373095*fIn[17])*dS; 
    out[14] += (2.449489742783178*fIn[59]-1.414213562373095*fIn[31])*dS; 
    out[15] += (2.449489742783178*fIn[60]-1.414213562373095*fIn[32])*dS; 
 
  } else {
 
    out[8] += (3.16227766016838*fIn[19]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[9] += (3.16227766016838*fIn[40]+2.449489742783178*fIn[9]+1.414213562373095*fIn[1])*dS; 
    out[10] += (3.16227766016838*fIn[41]+2.449489742783178*fIn[10]+1.414213562373095*fIn[2])*dS; 
    out[11] += (3.16227766016838*fIn[65]+2.449489742783178*fIn[22]+1.414213562373095*fIn[6])*dS; 
    out[12] += (2.449489742783178*fIn[37]+1.414213562373095*fIn[16])*dS; 
    out[13] += (2.449489742783178*fIn[38]+1.414213562373095*fIn[17])*dS; 
    out[14] += (2.449489742783178*fIn[59]+1.414213562373095*fIn[31])*dS; 
    out[15] += (2.449489742783178*fIn[60]+1.414213562373095*fIn[32])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[4]-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (2.449489742783178*fIn[9]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[9]-0.7071067811865475*fIn[1]*dxv[3])*dS; 
    out[2] += (2.449489742783178*fIn[10]-1.414213562373095*fIn[2])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[10]-0.7071067811865475*fIn[2]*dxv[3])*dS; 
    out[3] += (2.449489742783178*fIn[17]-1.414213562373095*fIn[6])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[17]-0.7071067811865475*dxv[3]*fIn[6])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[4])-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (2.449489742783178*fIn[9]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[9])-0.7071067811865475*fIn[1]*dxv[3])*dS; 
    out[2] += (2.449489742783178*fIn[10]+1.414213562373095*fIn[2])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[10])-0.7071067811865475*fIn[2]*dxv[3])*dS; 
    out[3] += (2.449489742783178*fIn[17]+1.414213562373095*fIn[6])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[17])-0.7071067811865475*dxv[3]*fIn[6])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[112]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[19])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[40])+2.449489742783178*fIn[9]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += ((-3.16227766016838*fIn[41])+2.449489742783178*fIn[10]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += ((-3.16227766016838*fIn[65])+2.449489742783178*fIn[22]-1.414213562373095*fIn[6])*dS*vBoundary; 
    out[4] += (2.449489742783178*fIn[37]-1.414213562373095*fIn[16])*dS*vBoundary; 
    out[5] += (2.449489742783178*fIn[38]-1.414213562373095*fIn[17])*dS*vBoundary; 
    out[6] += (2.449489742783178*fIn[59]-1.414213562373095*fIn[31])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[60]-1.414213562373095*fIn[32])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[19]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[40]+2.449489742783178*fIn[9]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (3.16227766016838*fIn[41]+2.449489742783178*fIn[10]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (3.16227766016838*fIn[65]+2.449489742783178*fIn[22]+1.414213562373095*fIn[6])*dS*vBoundary; 
    out[4] += (2.449489742783178*fIn[37]+1.414213562373095*fIn[16])*dS*vBoundary; 
    out[5] += (2.449489742783178*fIn[38]+1.414213562373095*fIn[17])*dS*vBoundary; 
    out[6] += (2.449489742783178*fIn[59]+1.414213562373095*fIn[31])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[60]+1.414213562373095*fIn[32])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[8] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
    out[9] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[2])*dS; 
    out[11] += (2.449489742783178*fIn[20]-1.414213562373095*fIn[6])*dS; 
 
  } else {
 
    out[8] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS; 
    out[9] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[2])*dS; 
    out[11] += (2.449489742783178*fIn[20]+1.414213562373095*fIn[6])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[112]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[16] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS; 
    out[17] += ((-3.16227766016838*fIn[47])+2.449489742783178*fIn[12]-1.414213562373095*fIn[1])*dS; 
    out[18] += ((-3.16227766016838*fIn[48])+2.449489742783178*fIn[13]-1.414213562373095*fIn[2])*dS; 
    out[19] += ((-3.16227766016838*fIn[80])+2.449489742783178*fIn[25]-1.414213562373095*fIn[6])*dS; 
    out[20] += (2.449489742783178*fIn[43]-1.414213562373095*fIn[16])*dS; 
    out[21] += (2.449489742783178*fIn[44]-1.414213562373095*fIn[17])*dS; 
    out[22] += (2.449489742783178*fIn[68]-1.414213562373095*fIn[31])*dS; 
    out[23] += (2.449489742783178*fIn[69]-1.414213562373095*fIn[32])*dS; 
 
  } else {
 
    out[16] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS; 
    out[17] += (3.16227766016838*fIn[47]+2.449489742783178*fIn[12]+1.414213562373095*fIn[1])*dS; 
    out[18] += (3.16227766016838*fIn[48]+2.449489742783178*fIn[13]+1.414213562373095*fIn[2])*dS; 
    out[19] += (3.16227766016838*fIn[80]+2.449489742783178*fIn[25]+1.414213562373095*fIn[6])*dS; 
    out[20] += (2.449489742783178*fIn[43]+1.414213562373095*fIn[16])*dS; 
    out[21] += (2.449489742783178*fIn[44]+1.414213562373095*fIn[17])*dS; 
    out[22] += (2.449489742783178*fIn[68]+1.414213562373095*fIn[31])*dS; 
    out[23] += (2.449489742783178*fIn[69]+1.414213562373095*fIn[32])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[5]-0.7071067811865475*fIn[0]*dxv[4])*dS; 
    out[1] += (2.449489742783178*fIn[12]-1.414213562373095*fIn[1])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[12]-0.7071067811865475*fIn[1]*dxv[4])*dS; 
    out[2] += (2.449489742783178*fIn[13]-1.414213562373095*fIn[2])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[13]-0.7071067811865475*fIn[2]*dxv[4])*dS; 
    out[3] += (2.449489742783178*fIn[20]-1.414213562373095*fIn[6])*dS*vBoundary+(1.224744871391589*dxv[4]*fIn[20]-0.7071067811865475*dxv[4]*fIn[6])*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[5])-0.7071067811865475*fIn[0]*dxv[4])*dS; 
    out[1] += (2.449489742783178*fIn[12]+1.414213562373095*fIn[1])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[12])-0.7071067811865475*fIn[1]*dxv[4])*dS; 
    out[2] += (2.449489742783178*fIn[13]+1.414213562373095*fIn[2])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[13])-0.7071067811865475*fIn[2]*dxv[4])*dS; 
    out[3] += (2.449489742783178*fIn[20]+1.414213562373095*fIn[6])*dS*vBoundary+((-1.224744871391589*dxv[4]*fIn[20])-0.7071067811865475*dxv[4]*fIn[6])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x3vSer_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[112]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[47])+2.449489742783178*fIn[12]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += ((-3.16227766016838*fIn[48])+2.449489742783178*fIn[13]-1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += ((-3.16227766016838*fIn[80])+2.449489742783178*fIn[25]-1.414213562373095*fIn[6])*dS*vBoundary; 
    out[4] += (2.449489742783178*fIn[43]-1.414213562373095*fIn[16])*dS*vBoundary; 
    out[5] += (2.449489742783178*fIn[44]-1.414213562373095*fIn[17])*dS*vBoundary; 
    out[6] += (2.449489742783178*fIn[68]-1.414213562373095*fIn[31])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[69]-1.414213562373095*fIn[32])*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[47]+2.449489742783178*fIn[12]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (3.16227766016838*fIn[48]+2.449489742783178*fIn[13]+1.414213562373095*fIn[2])*dS*vBoundary; 
    out[3] += (3.16227766016838*fIn[80]+2.449489742783178*fIn[25]+1.414213562373095*fIn[6])*dS*vBoundary; 
    out[4] += (2.449489742783178*fIn[43]+1.414213562373095*fIn[16])*dS*vBoundary; 
    out[5] += (2.449489742783178*fIn[44]+1.414213562373095*fIn[17])*dS*vBoundary; 
    out[6] += (2.449489742783178*fIn[68]+1.414213562373095*fIn[31])*dS*vBoundary; 
    out[7] += (2.449489742783178*fIn[69]+1.414213562373095*fIn[32])*dS*vBoundary; 
 
  }
 
} 
 
void VmSelfPrimMoments2x3vSer_P1(binOpData_t* data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1:       moments of the distribution function. 
  // m0S,m1S,m1S: star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.5*m0[3])-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.5*m0[3])-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
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
  data->AEM_S(0,0) = 0.5*m0r[0]; 
  data->AEM_S(0,1) = 0.5*m0r[1]; 
  data->AEM_S(0,2) = 0.5*m0r[2]; 
  data->AEM_S(0,3) = 0.5*m0r[3]; 
  data->AEM_S(1,0) = 0.5*m0r[1]; 
  data->AEM_S(1,1) = 0.5*m0r[0]; 
  data->AEM_S(1,2) = 0.5*m0r[3]; 
  data->AEM_S(1,3) = 0.5*m0r[2]; 
  data->AEM_S(2,0) = 0.5*m0r[2]; 
  data->AEM_S(2,1) = 0.5*m0r[3]; 
  data->AEM_S(2,2) = 0.5*m0r[0]; 
  data->AEM_S(2,3) = 0.5*m0r[1]; 
  data->AEM_S(3,0) = 0.5*m0r[3]; 
  data->AEM_S(3,1) = 0.5*m0r[2]; 
  data->AEM_S(3,2) = 0.5*m0r[1]; 
  data->AEM_S(3,3) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,12) = -0.5*cM[0]; 
  data->AEM_S(0,13) = -0.5*cM[1]; 
  data->AEM_S(0,14) = -0.5*cM[2]; 
  data->AEM_S(0,15) = -0.5*cM[3]; 
  data->AEM_S(1,12) = -0.5*cM[1]; 
  data->AEM_S(1,13) = -0.5*cM[0]; 
  data->AEM_S(1,14) = -0.5*cM[3]; 
  data->AEM_S(1,15) = -0.5*cM[2]; 
  data->AEM_S(2,12) = -0.5*cM[2]; 
  data->AEM_S(2,13) = -0.5*cM[3]; 
  data->AEM_S(2,14) = -0.5*cM[0]; 
  data->AEM_S(2,15) = -0.5*cM[1]; 
  data->AEM_S(3,12) = -0.5*cM[3]; 
  data->AEM_S(3,13) = -0.5*cM[2]; 
  data->AEM_S(3,14) = -0.5*cM[1]; 
  data->AEM_S(3,15) = -0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(12,0) = 0.5*m1Sr[0]; 
  data->AEM_S(12,1) = 0.5*m1Sr[1]; 
  data->AEM_S(12,2) = 0.5*m1Sr[2]; 
  data->AEM_S(12,3) = 0.5*m1Sr[3]; 
  data->AEM_S(13,0) = 0.5*m1Sr[1]; 
  data->AEM_S(13,1) = 0.5*m1Sr[0]; 
  data->AEM_S(13,2) = 0.5*m1Sr[3]; 
  data->AEM_S(13,3) = 0.5*m1Sr[2]; 
  data->AEM_S(14,0) = 0.5*m1Sr[2]; 
  data->AEM_S(14,1) = 0.5*m1Sr[3]; 
  data->AEM_S(14,2) = 0.5*m1Sr[0]; 
  data->AEM_S(14,3) = 0.5*m1Sr[1]; 
  data->AEM_S(15,0) = 0.5*m1Sr[3]; 
  data->AEM_S(15,1) = 0.5*m1Sr[2]; 
  data->AEM_S(15,2) = 0.5*m1Sr[1]; 
  data->AEM_S(15,3) = 0.5*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(4,4) = 0.5*m0r[0]; 
  data->AEM_S(4,5) = 0.5*m0r[1]; 
  data->AEM_S(4,6) = 0.5*m0r[2]; 
  data->AEM_S(4,7) = 0.5*m0r[3]; 
  data->AEM_S(5,4) = 0.5*m0r[1]; 
  data->AEM_S(5,5) = 0.5*m0r[0]; 
  data->AEM_S(5,6) = 0.5*m0r[3]; 
  data->AEM_S(5,7) = 0.5*m0r[2]; 
  data->AEM_S(6,4) = 0.5*m0r[2]; 
  data->AEM_S(6,5) = 0.5*m0r[3]; 
  data->AEM_S(6,6) = 0.5*m0r[0]; 
  data->AEM_S(6,7) = 0.5*m0r[1]; 
  data->AEM_S(7,4) = 0.5*m0r[3]; 
  data->AEM_S(7,5) = 0.5*m0r[2]; 
  data->AEM_S(7,6) = 0.5*m0r[1]; 
  data->AEM_S(7,7) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(4,12) = -0.5*cM[4]; 
  data->AEM_S(4,13) = -0.5*cM[5]; 
  data->AEM_S(4,14) = -0.5*cM[6]; 
  data->AEM_S(4,15) = -0.5*cM[7]; 
  data->AEM_S(5,12) = -0.5*cM[5]; 
  data->AEM_S(5,13) = -0.5*cM[4]; 
  data->AEM_S(5,14) = -0.5*cM[7]; 
  data->AEM_S(5,15) = -0.5*cM[6]; 
  data->AEM_S(6,12) = -0.5*cM[6]; 
  data->AEM_S(6,13) = -0.5*cM[7]; 
  data->AEM_S(6,14) = -0.5*cM[4]; 
  data->AEM_S(6,15) = -0.5*cM[5]; 
  data->AEM_S(7,12) = -0.5*cM[7]; 
  data->AEM_S(7,13) = -0.5*cM[6]; 
  data->AEM_S(7,14) = -0.5*cM[5]; 
  data->AEM_S(7,15) = -0.5*cM[4]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(12,4) = 0.5*m1Sr[4]; 
  data->AEM_S(12,5) = 0.5*m1Sr[5]; 
  data->AEM_S(12,6) = 0.5*m1Sr[6]; 
  data->AEM_S(12,7) = 0.5*m1Sr[7]; 
  data->AEM_S(13,4) = 0.5*m1Sr[5]; 
  data->AEM_S(13,5) = 0.5*m1Sr[4]; 
  data->AEM_S(13,6) = 0.5*m1Sr[7]; 
  data->AEM_S(13,7) = 0.5*m1Sr[6]; 
  data->AEM_S(14,4) = 0.5*m1Sr[6]; 
  data->AEM_S(14,5) = 0.5*m1Sr[7]; 
  data->AEM_S(14,6) = 0.5*m1Sr[4]; 
  data->AEM_S(14,7) = 0.5*m1Sr[5]; 
  data->AEM_S(15,4) = 0.5*m1Sr[7]; 
  data->AEM_S(15,5) = 0.5*m1Sr[6]; 
  data->AEM_S(15,6) = 0.5*m1Sr[5]; 
  data->AEM_S(15,7) = 0.5*m1Sr[4]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(8,8) = 0.5*m0r[0]; 
  data->AEM_S(8,9) = 0.5*m0r[1]; 
  data->AEM_S(8,10) = 0.5*m0r[2]; 
  data->AEM_S(8,11) = 0.5*m0r[3]; 
  data->AEM_S(9,8) = 0.5*m0r[1]; 
  data->AEM_S(9,9) = 0.5*m0r[0]; 
  data->AEM_S(9,10) = 0.5*m0r[3]; 
  data->AEM_S(9,11) = 0.5*m0r[2]; 
  data->AEM_S(10,8) = 0.5*m0r[2]; 
  data->AEM_S(10,9) = 0.5*m0r[3]; 
  data->AEM_S(10,10) = 0.5*m0r[0]; 
  data->AEM_S(10,11) = 0.5*m0r[1]; 
  data->AEM_S(11,8) = 0.5*m0r[3]; 
  data->AEM_S(11,9) = 0.5*m0r[2]; 
  data->AEM_S(11,10) = 0.5*m0r[1]; 
  data->AEM_S(11,11) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(8,12) = -0.5*cM[8]; 
  data->AEM_S(8,13) = -0.5*cM[9]; 
  data->AEM_S(8,14) = -0.5*cM[10]; 
  data->AEM_S(8,15) = -0.5*cM[11]; 
  data->AEM_S(9,12) = -0.5*cM[9]; 
  data->AEM_S(9,13) = -0.5*cM[8]; 
  data->AEM_S(9,14) = -0.5*cM[11]; 
  data->AEM_S(9,15) = -0.5*cM[10]; 
  data->AEM_S(10,12) = -0.5*cM[10]; 
  data->AEM_S(10,13) = -0.5*cM[11]; 
  data->AEM_S(10,14) = -0.5*cM[8]; 
  data->AEM_S(10,15) = -0.5*cM[9]; 
  data->AEM_S(11,12) = -0.5*cM[11]; 
  data->AEM_S(11,13) = -0.5*cM[10]; 
  data->AEM_S(11,14) = -0.5*cM[9]; 
  data->AEM_S(11,15) = -0.5*cM[8]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(12,8) = 0.5*m1Sr[8]; 
  data->AEM_S(12,9) = 0.5*m1Sr[9]; 
  data->AEM_S(12,10) = 0.5*m1Sr[10]; 
  data->AEM_S(12,11) = 0.5*m1Sr[11]; 
  data->AEM_S(13,8) = 0.5*m1Sr[9]; 
  data->AEM_S(13,9) = 0.5*m1Sr[8]; 
  data->AEM_S(13,10) = 0.5*m1Sr[11]; 
  data->AEM_S(13,11) = 0.5*m1Sr[10]; 
  data->AEM_S(14,8) = 0.5*m1Sr[10]; 
  data->AEM_S(14,9) = 0.5*m1Sr[11]; 
  data->AEM_S(14,10) = 0.5*m1Sr[8]; 
  data->AEM_S(14,11) = 0.5*m1Sr[9]; 
  data->AEM_S(15,8) = 0.5*m1Sr[11]; 
  data->AEM_S(15,9) = 0.5*m1Sr[10]; 
  data->AEM_S(15,10) = 0.5*m1Sr[9]; 
  data->AEM_S(15,11) = 0.5*m1Sr[8]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(12,12) = 0.5*m0Sr[0]-0.5*cE[0]; 
  data->AEM_S(12,13) = 0.5*m0Sr[1]-0.5*cE[1]; 
  data->AEM_S(12,14) = 0.5*m0Sr[2]-0.5*cE[2]; 
  data->AEM_S(12,15) = 0.5*m0Sr[3]-0.5*cE[3]; 
  data->AEM_S(13,12) = 0.5*m0Sr[1]-0.5*cE[1]; 
  data->AEM_S(13,13) = 0.5*m0Sr[0]-0.5*cE[0]; 
  data->AEM_S(13,14) = 0.5*m0Sr[3]-0.5*cE[3]; 
  data->AEM_S(13,15) = 0.5*m0Sr[2]-0.5*cE[2]; 
  data->AEM_S(14,12) = 0.5*m0Sr[2]-0.5*cE[2]; 
  data->AEM_S(14,13) = 0.5*m0Sr[3]-0.5*cE[3]; 
  data->AEM_S(14,14) = 0.5*m0Sr[0]-0.5*cE[0]; 
  data->AEM_S(14,15) = 0.5*m0Sr[1]-0.5*cE[1]; 
  data->AEM_S(15,12) = 0.5*m0Sr[3]-0.5*cE[3]; 
  data->AEM_S(15,13) = 0.5*m0Sr[2]-0.5*cE[2]; 
  data->AEM_S(15,14) = 0.5*m0Sr[1]-0.5*cE[1]; 
  data->AEM_S(15,15) = 0.5*m0Sr[0]-0.5*cE[0]; 
 
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
 
void VmSelfPrimMoments2x3vSer_P2(binOpData_t* data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-1.936491673103709*m0[7])-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.936491673103709*m0[7])-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[8]; 
  double m1r[24]; 
  double m2r[8]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m0r[4] = 0.0; 
    m0r[5] = 0.0; 
    m0r[6] = 0.0; 
    m0r[7] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    m1r[8] = m1[8]; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    m1r[12] = 0.0; 
    m1r[13] = 0.0; 
    m1r[14] = 0.0; 
    m1r[15] = 0.0; 
    m1r[16] = m1[16]; 
    m1r[17] = 0.0; 
    m1r[18] = 0.0; 
    m1r[19] = 0.0; 
    m1r[20] = 0.0; 
    m1r[21] = 0.0; 
    m1r[22] = 0.0; 
    m1r[23] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    m2r[4] = 0.0; 
    m2r[5] = 0.0; 
    m2r[6] = 0.0; 
    m2r[7] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
    m0r[4] = m0[4]; 
    m0r[5] = m0[5]; 
    m0r[6] = m0[6]; 
    m0r[7] = m0[7]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
    m2r[4] = m2[4]; 
    m2r[5] = m2[5]; 
    m2r[6] = m2[6]; 
    m2r[7] = m2[7]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(32,32); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.5*m0r[0]; 
  data->AEM_S(0,1) = 0.5*m0r[1]; 
  data->AEM_S(0,2) = 0.5*m0r[2]; 
  data->AEM_S(0,3) = 0.5*m0r[3]; 
  data->AEM_S(0,4) = 0.5*m0r[4]; 
  data->AEM_S(0,5) = 0.5*m0r[5]; 
  data->AEM_S(0,6) = 0.5*m0r[6]; 
  data->AEM_S(0,7) = 0.5*m0r[7]; 
  data->AEM_S(1,0) = 0.5*m0r[1]; 
  data->AEM_S(1,1) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(1,2) = 0.5*m0r[3]; 
  data->AEM_S(1,3) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  data->AEM_S(1,4) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(1,5) = 0.5000000000000001*m0r[7]; 
  data->AEM_S(1,6) = 0.447213595499958*m0r[3]; 
  data->AEM_S(1,7) = 0.5000000000000001*m0r[5]; 
  data->AEM_S(2,0) = 0.5*m0r[2]; 
  data->AEM_S(2,1) = 0.5*m0r[3]; 
  data->AEM_S(2,2) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(2,3) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  data->AEM_S(2,4) = 0.5000000000000001*m0r[6]; 
  data->AEM_S(2,5) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(2,6) = 0.5000000000000001*m0r[4]; 
  data->AEM_S(2,7) = 0.447213595499958*m0r[3]; 
  data->AEM_S(3,0) = 0.5*m0r[3]; 
  data->AEM_S(3,1) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  data->AEM_S(3,2) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  data->AEM_S(3,3) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(3,4) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(3,5) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(3,6) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  data->AEM_S(3,7) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  data->AEM_S(4,0) = 0.5*m0r[4]; 
  data->AEM_S(4,1) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(4,2) = 0.5000000000000001*m0r[6]; 
  data->AEM_S(4,3) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(4,4) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(4,6) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  data->AEM_S(4,7) = 0.4472135954999579*m0r[7]; 
  data->AEM_S(5,0) = 0.5*m0r[5]; 
  data->AEM_S(5,1) = 0.5000000000000001*m0r[7]; 
  data->AEM_S(5,2) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(5,3) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(5,5) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(5,6) = 0.4472135954999579*m0r[6]; 
  data->AEM_S(5,7) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  data->AEM_S(6,0) = 0.5*m0r[6]; 
  data->AEM_S(6,1) = 0.447213595499958*m0r[3]; 
  data->AEM_S(6,2) = 0.5000000000000001*m0r[4]; 
  data->AEM_S(6,3) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  data->AEM_S(6,4) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  data->AEM_S(6,5) = 0.4472135954999579*m0r[6]; 
  data->AEM_S(6,6) = 0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(6,7) = 0.4*m0r[3]; 
  data->AEM_S(7,0) = 0.5*m0r[7]; 
  data->AEM_S(7,1) = 0.5000000000000001*m0r[5]; 
  data->AEM_S(7,2) = 0.447213595499958*m0r[3]; 
  data->AEM_S(7,3) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  data->AEM_S(7,4) = 0.4472135954999579*m0r[7]; 
  data->AEM_S(7,5) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  data->AEM_S(7,6) = 0.4*m0r[3]; 
  data->AEM_S(7,7) = 0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,24) = -0.5*cM[0]; 
  data->AEM_S(0,25) = -0.5*cM[1]; 
  data->AEM_S(0,26) = -0.5*cM[2]; 
  data->AEM_S(0,27) = -0.5*cM[3]; 
  data->AEM_S(0,28) = -0.5*cM[4]; 
  data->AEM_S(0,29) = -0.5*cM[5]; 
  data->AEM_S(0,30) = -0.5*cM[6]; 
  data->AEM_S(0,31) = -0.5*cM[7]; 
  data->AEM_S(1,24) = -0.5*cM[1]; 
  data->AEM_S(1,25) = (-0.4472135954999579*cM[4])-0.5*cM[0]; 
  data->AEM_S(1,26) = -0.5*cM[3]; 
  data->AEM_S(1,27) = (-0.447213595499958*cM[6])-0.5*cM[2]; 
  data->AEM_S(1,28) = -0.4472135954999579*cM[1]; 
  data->AEM_S(1,29) = -0.5000000000000001*cM[7]; 
  data->AEM_S(1,30) = -0.447213595499958*cM[3]; 
  data->AEM_S(1,31) = -0.5000000000000001*cM[5]; 
  data->AEM_S(2,24) = -0.5*cM[2]; 
  data->AEM_S(2,25) = -0.5*cM[3]; 
  data->AEM_S(2,26) = (-0.4472135954999579*cM[5])-0.5*cM[0]; 
  data->AEM_S(2,27) = (-0.447213595499958*cM[7])-0.5*cM[1]; 
  data->AEM_S(2,28) = -0.5000000000000001*cM[6]; 
  data->AEM_S(2,29) = -0.4472135954999579*cM[2]; 
  data->AEM_S(2,30) = -0.5000000000000001*cM[4]; 
  data->AEM_S(2,31) = -0.447213595499958*cM[3]; 
  data->AEM_S(3,24) = -0.5*cM[3]; 
  data->AEM_S(3,25) = (-0.447213595499958*cM[6])-0.5*cM[2]; 
  data->AEM_S(3,26) = (-0.447213595499958*cM[7])-0.5*cM[1]; 
  data->AEM_S(3,27) = (-0.4472135954999579*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  data->AEM_S(3,28) = -0.4472135954999579*cM[3]; 
  data->AEM_S(3,29) = -0.4472135954999579*cM[3]; 
  data->AEM_S(3,30) = (-0.4*cM[7])-0.447213595499958*cM[1]; 
  data->AEM_S(3,31) = (-0.4*cM[6])-0.447213595499958*cM[2]; 
  data->AEM_S(4,24) = -0.5*cM[4]; 
  data->AEM_S(4,25) = -0.4472135954999579*cM[1]; 
  data->AEM_S(4,26) = -0.5000000000000001*cM[6]; 
  data->AEM_S(4,27) = -0.4472135954999579*cM[3]; 
  data->AEM_S(4,28) = (-0.31943828249997*cM[4])-0.5*cM[0]; 
  data->AEM_S(4,30) = (-0.31943828249997*cM[6])-0.5000000000000001*cM[2]; 
  data->AEM_S(4,31) = -0.4472135954999579*cM[7]; 
  data->AEM_S(5,24) = -0.5*cM[5]; 
  data->AEM_S(5,25) = -0.5000000000000001*cM[7]; 
  data->AEM_S(5,26) = -0.4472135954999579*cM[2]; 
  data->AEM_S(5,27) = -0.4472135954999579*cM[3]; 
  data->AEM_S(5,29) = (-0.31943828249997*cM[5])-0.5*cM[0]; 
  data->AEM_S(5,30) = -0.4472135954999579*cM[6]; 
  data->AEM_S(5,31) = (-0.31943828249997*cM[7])-0.5000000000000001*cM[1]; 
  data->AEM_S(6,24) = -0.5*cM[6]; 
  data->AEM_S(6,25) = -0.447213595499958*cM[3]; 
  data->AEM_S(6,26) = -0.5000000000000001*cM[4]; 
  data->AEM_S(6,27) = (-0.4*cM[7])-0.447213595499958*cM[1]; 
  data->AEM_S(6,28) = (-0.31943828249997*cM[6])-0.5000000000000001*cM[2]; 
  data->AEM_S(6,29) = -0.4472135954999579*cM[6]; 
  data->AEM_S(6,30) = (-0.4472135954999579*cM[5])-0.31943828249997*cM[4]-0.5*cM[0]; 
  data->AEM_S(6,31) = -0.4*cM[3]; 
  data->AEM_S(7,24) = -0.5*cM[7]; 
  data->AEM_S(7,25) = -0.5000000000000001*cM[5]; 
  data->AEM_S(7,26) = -0.447213595499958*cM[3]; 
  data->AEM_S(7,27) = (-0.4*cM[6])-0.447213595499958*cM[2]; 
  data->AEM_S(7,28) = -0.4472135954999579*cM[7]; 
  data->AEM_S(7,29) = (-0.31943828249997*cM[7])-0.5000000000000001*cM[1]; 
  data->AEM_S(7,30) = -0.4*cM[3]; 
  data->AEM_S(7,31) = (-0.31943828249997*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(24,0) = 0.5*m1r[0]; 
  data->AEM_S(24,1) = 0.5*m1r[1]; 
  data->AEM_S(24,2) = 0.5*m1r[2]; 
  data->AEM_S(24,3) = 0.5*m1r[3]; 
  data->AEM_S(24,4) = 0.5*m1r[4]; 
  data->AEM_S(24,5) = 0.5*m1r[5]; 
  data->AEM_S(24,6) = 0.5*m1r[6]; 
  data->AEM_S(24,7) = 0.5*m1r[7]; 
  data->AEM_S(25,0) = 0.5*m1r[1]; 
  data->AEM_S(25,1) = 0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(25,2) = 0.5*m1r[3]; 
  data->AEM_S(25,3) = 0.447213595499958*m1r[6]+0.5*m1r[2]; 
  data->AEM_S(25,4) = 0.4472135954999579*m1r[1]; 
  data->AEM_S(25,5) = 0.5000000000000001*m1r[7]; 
  data->AEM_S(25,6) = 0.447213595499958*m1r[3]; 
  data->AEM_S(25,7) = 0.5000000000000001*m1r[5]; 
  data->AEM_S(26,0) = 0.5*m1r[2]; 
  data->AEM_S(26,1) = 0.5*m1r[3]; 
  data->AEM_S(26,2) = 0.4472135954999579*m1r[5]+0.5*m1r[0]; 
  data->AEM_S(26,3) = 0.447213595499958*m1r[7]+0.5*m1r[1]; 
  data->AEM_S(26,4) = 0.5000000000000001*m1r[6]; 
  data->AEM_S(26,5) = 0.4472135954999579*m1r[2]; 
  data->AEM_S(26,6) = 0.5000000000000001*m1r[4]; 
  data->AEM_S(26,7) = 0.447213595499958*m1r[3]; 
  data->AEM_S(27,0) = 0.5*m1r[3]; 
  data->AEM_S(27,1) = 0.447213595499958*m1r[6]+0.5*m1r[2]; 
  data->AEM_S(27,2) = 0.447213595499958*m1r[7]+0.5*m1r[1]; 
  data->AEM_S(27,3) = 0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(27,4) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(27,5) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(27,6) = 0.4*m1r[7]+0.447213595499958*m1r[1]; 
  data->AEM_S(27,7) = 0.4*m1r[6]+0.447213595499958*m1r[2]; 
  data->AEM_S(28,0) = 0.5*m1r[4]; 
  data->AEM_S(28,1) = 0.4472135954999579*m1r[1]; 
  data->AEM_S(28,2) = 0.5000000000000001*m1r[6]; 
  data->AEM_S(28,3) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(28,4) = 0.31943828249997*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(28,6) = 0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]; 
  data->AEM_S(28,7) = 0.4472135954999579*m1r[7]; 
  data->AEM_S(29,0) = 0.5*m1r[5]; 
  data->AEM_S(29,1) = 0.5000000000000001*m1r[7]; 
  data->AEM_S(29,2) = 0.4472135954999579*m1r[2]; 
  data->AEM_S(29,3) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(29,5) = 0.31943828249997*m1r[5]+0.5*m1r[0]; 
  data->AEM_S(29,6) = 0.4472135954999579*m1r[6]; 
  data->AEM_S(29,7) = 0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]; 
  data->AEM_S(30,0) = 0.5*m1r[6]; 
  data->AEM_S(30,1) = 0.447213595499958*m1r[3]; 
  data->AEM_S(30,2) = 0.5000000000000001*m1r[4]; 
  data->AEM_S(30,3) = 0.4*m1r[7]+0.447213595499958*m1r[1]; 
  data->AEM_S(30,4) = 0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]; 
  data->AEM_S(30,5) = 0.4472135954999579*m1r[6]; 
  data->AEM_S(30,6) = 0.4472135954999579*m1r[5]+0.31943828249997*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(30,7) = 0.4*m1r[3]; 
  data->AEM_S(31,0) = 0.5*m1r[7]; 
  data->AEM_S(31,1) = 0.5000000000000001*m1r[5]; 
  data->AEM_S(31,2) = 0.447213595499958*m1r[3]; 
  data->AEM_S(31,3) = 0.4*m1r[6]+0.447213595499958*m1r[2]; 
  data->AEM_S(31,4) = 0.4472135954999579*m1r[7]; 
  data->AEM_S(31,5) = 0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]; 
  data->AEM_S(31,6) = 0.4*m1r[3]; 
  data->AEM_S(31,7) = 0.31943828249997*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(8,8) = 0.5*m0r[0]; 
  data->AEM_S(8,9) = 0.5*m0r[1]; 
  data->AEM_S(8,10) = 0.5*m0r[2]; 
  data->AEM_S(8,11) = 0.5*m0r[3]; 
  data->AEM_S(8,12) = 0.5*m0r[4]; 
  data->AEM_S(8,13) = 0.5*m0r[5]; 
  data->AEM_S(8,14) = 0.5*m0r[6]; 
  data->AEM_S(8,15) = 0.5*m0r[7]; 
  data->AEM_S(9,8) = 0.5*m0r[1]; 
  data->AEM_S(9,9) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(9,10) = 0.5*m0r[3]; 
  data->AEM_S(9,11) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  data->AEM_S(9,12) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(9,13) = 0.5000000000000001*m0r[7]; 
  data->AEM_S(9,14) = 0.447213595499958*m0r[3]; 
  data->AEM_S(9,15) = 0.5000000000000001*m0r[5]; 
  data->AEM_S(10,8) = 0.5*m0r[2]; 
  data->AEM_S(10,9) = 0.5*m0r[3]; 
  data->AEM_S(10,10) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(10,11) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  data->AEM_S(10,12) = 0.5000000000000001*m0r[6]; 
  data->AEM_S(10,13) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(10,14) = 0.5000000000000001*m0r[4]; 
  data->AEM_S(10,15) = 0.447213595499958*m0r[3]; 
  data->AEM_S(11,8) = 0.5*m0r[3]; 
  data->AEM_S(11,9) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  data->AEM_S(11,10) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  data->AEM_S(11,11) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(11,12) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(11,13) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(11,14) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  data->AEM_S(11,15) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  data->AEM_S(12,8) = 0.5*m0r[4]; 
  data->AEM_S(12,9) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(12,10) = 0.5000000000000001*m0r[6]; 
  data->AEM_S(12,11) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(12,12) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(12,14) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  data->AEM_S(12,15) = 0.4472135954999579*m0r[7]; 
  data->AEM_S(13,8) = 0.5*m0r[5]; 
  data->AEM_S(13,9) = 0.5000000000000001*m0r[7]; 
  data->AEM_S(13,10) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(13,11) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(13,13) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(13,14) = 0.4472135954999579*m0r[6]; 
  data->AEM_S(13,15) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  data->AEM_S(14,8) = 0.5*m0r[6]; 
  data->AEM_S(14,9) = 0.447213595499958*m0r[3]; 
  data->AEM_S(14,10) = 0.5000000000000001*m0r[4]; 
  data->AEM_S(14,11) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  data->AEM_S(14,12) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  data->AEM_S(14,13) = 0.4472135954999579*m0r[6]; 
  data->AEM_S(14,14) = 0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(14,15) = 0.4*m0r[3]; 
  data->AEM_S(15,8) = 0.5*m0r[7]; 
  data->AEM_S(15,9) = 0.5000000000000001*m0r[5]; 
  data->AEM_S(15,10) = 0.447213595499958*m0r[3]; 
  data->AEM_S(15,11) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  data->AEM_S(15,12) = 0.4472135954999579*m0r[7]; 
  data->AEM_S(15,13) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  data->AEM_S(15,14) = 0.4*m0r[3]; 
  data->AEM_S(15,15) = 0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(8,24) = -0.5*cM[8]; 
  data->AEM_S(8,25) = -0.5*cM[9]; 
  data->AEM_S(8,26) = -0.5*cM[10]; 
  data->AEM_S(8,27) = -0.5*cM[11]; 
  data->AEM_S(8,28) = -0.5*cM[12]; 
  data->AEM_S(8,29) = -0.5*cM[13]; 
  data->AEM_S(8,30) = -0.5*cM[14]; 
  data->AEM_S(8,31) = -0.5*cM[15]; 
  data->AEM_S(9,24) = -0.5*cM[9]; 
  data->AEM_S(9,25) = (-0.4472135954999579*cM[12])-0.5*cM[8]; 
  data->AEM_S(9,26) = -0.5*cM[11]; 
  data->AEM_S(9,27) = (-0.447213595499958*cM[14])-0.5*cM[10]; 
  data->AEM_S(9,28) = -0.4472135954999579*cM[9]; 
  data->AEM_S(9,29) = -0.5000000000000001*cM[15]; 
  data->AEM_S(9,30) = -0.447213595499958*cM[11]; 
  data->AEM_S(9,31) = -0.5000000000000001*cM[13]; 
  data->AEM_S(10,24) = -0.5*cM[10]; 
  data->AEM_S(10,25) = -0.5*cM[11]; 
  data->AEM_S(10,26) = (-0.4472135954999579*cM[13])-0.5*cM[8]; 
  data->AEM_S(10,27) = (-0.447213595499958*cM[15])-0.5*cM[9]; 
  data->AEM_S(10,28) = -0.5000000000000001*cM[14]; 
  data->AEM_S(10,29) = -0.4472135954999579*cM[10]; 
  data->AEM_S(10,30) = -0.5000000000000001*cM[12]; 
  data->AEM_S(10,31) = -0.447213595499958*cM[11]; 
  data->AEM_S(11,24) = -0.5*cM[11]; 
  data->AEM_S(11,25) = (-0.447213595499958*cM[14])-0.5*cM[10]; 
  data->AEM_S(11,26) = (-0.447213595499958*cM[15])-0.5*cM[9]; 
  data->AEM_S(11,27) = (-0.4472135954999579*cM[13])-0.4472135954999579*cM[12]-0.5*cM[8]; 
  data->AEM_S(11,28) = -0.4472135954999579*cM[11]; 
  data->AEM_S(11,29) = -0.4472135954999579*cM[11]; 
  data->AEM_S(11,30) = (-0.4*cM[15])-0.447213595499958*cM[9]; 
  data->AEM_S(11,31) = (-0.4*cM[14])-0.447213595499958*cM[10]; 
  data->AEM_S(12,24) = -0.5*cM[12]; 
  data->AEM_S(12,25) = -0.4472135954999579*cM[9]; 
  data->AEM_S(12,26) = -0.5000000000000001*cM[14]; 
  data->AEM_S(12,27) = -0.4472135954999579*cM[11]; 
  data->AEM_S(12,28) = (-0.31943828249997*cM[12])-0.5*cM[8]; 
  data->AEM_S(12,30) = (-0.31943828249997*cM[14])-0.5000000000000001*cM[10]; 
  data->AEM_S(12,31) = -0.4472135954999579*cM[15]; 
  data->AEM_S(13,24) = -0.5*cM[13]; 
  data->AEM_S(13,25) = -0.5000000000000001*cM[15]; 
  data->AEM_S(13,26) = -0.4472135954999579*cM[10]; 
  data->AEM_S(13,27) = -0.4472135954999579*cM[11]; 
  data->AEM_S(13,29) = (-0.31943828249997*cM[13])-0.5*cM[8]; 
  data->AEM_S(13,30) = -0.4472135954999579*cM[14]; 
  data->AEM_S(13,31) = (-0.31943828249997*cM[15])-0.5000000000000001*cM[9]; 
  data->AEM_S(14,24) = -0.5*cM[14]; 
  data->AEM_S(14,25) = -0.447213595499958*cM[11]; 
  data->AEM_S(14,26) = -0.5000000000000001*cM[12]; 
  data->AEM_S(14,27) = (-0.4*cM[15])-0.447213595499958*cM[9]; 
  data->AEM_S(14,28) = (-0.31943828249997*cM[14])-0.5000000000000001*cM[10]; 
  data->AEM_S(14,29) = -0.4472135954999579*cM[14]; 
  data->AEM_S(14,30) = (-0.4472135954999579*cM[13])-0.31943828249997*cM[12]-0.5*cM[8]; 
  data->AEM_S(14,31) = -0.4*cM[11]; 
  data->AEM_S(15,24) = -0.5*cM[15]; 
  data->AEM_S(15,25) = -0.5000000000000001*cM[13]; 
  data->AEM_S(15,26) = -0.447213595499958*cM[11]; 
  data->AEM_S(15,27) = (-0.4*cM[14])-0.447213595499958*cM[10]; 
  data->AEM_S(15,28) = -0.4472135954999579*cM[15]; 
  data->AEM_S(15,29) = (-0.31943828249997*cM[15])-0.5000000000000001*cM[9]; 
  data->AEM_S(15,30) = -0.4*cM[11]; 
  data->AEM_S(15,31) = (-0.31943828249997*cM[13])-0.4472135954999579*cM[12]-0.5*cM[8]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(24,8) = 0.5*m1r[8]; 
  data->AEM_S(24,9) = 0.5*m1r[9]; 
  data->AEM_S(24,10) = 0.5*m1r[10]; 
  data->AEM_S(24,11) = 0.5*m1r[11]; 
  data->AEM_S(24,12) = 0.5*m1r[12]; 
  data->AEM_S(24,13) = 0.5*m1r[13]; 
  data->AEM_S(24,14) = 0.5*m1r[14]; 
  data->AEM_S(24,15) = 0.5*m1r[15]; 
  data->AEM_S(25,8) = 0.5*m1r[9]; 
  data->AEM_S(25,9) = 0.4472135954999579*m1r[12]+0.5*m1r[8]; 
  data->AEM_S(25,10) = 0.5*m1r[11]; 
  data->AEM_S(25,11) = 0.447213595499958*m1r[14]+0.5*m1r[10]; 
  data->AEM_S(25,12) = 0.4472135954999579*m1r[9]; 
  data->AEM_S(25,13) = 0.5000000000000001*m1r[15]; 
  data->AEM_S(25,14) = 0.447213595499958*m1r[11]; 
  data->AEM_S(25,15) = 0.5000000000000001*m1r[13]; 
  data->AEM_S(26,8) = 0.5*m1r[10]; 
  data->AEM_S(26,9) = 0.5*m1r[11]; 
  data->AEM_S(26,10) = 0.4472135954999579*m1r[13]+0.5*m1r[8]; 
  data->AEM_S(26,11) = 0.447213595499958*m1r[15]+0.5*m1r[9]; 
  data->AEM_S(26,12) = 0.5000000000000001*m1r[14]; 
  data->AEM_S(26,13) = 0.4472135954999579*m1r[10]; 
  data->AEM_S(26,14) = 0.5000000000000001*m1r[12]; 
  data->AEM_S(26,15) = 0.447213595499958*m1r[11]; 
  data->AEM_S(27,8) = 0.5*m1r[11]; 
  data->AEM_S(27,9) = 0.447213595499958*m1r[14]+0.5*m1r[10]; 
  data->AEM_S(27,10) = 0.447213595499958*m1r[15]+0.5*m1r[9]; 
  data->AEM_S(27,11) = 0.4472135954999579*m1r[13]+0.4472135954999579*m1r[12]+0.5*m1r[8]; 
  data->AEM_S(27,12) = 0.4472135954999579*m1r[11]; 
  data->AEM_S(27,13) = 0.4472135954999579*m1r[11]; 
  data->AEM_S(27,14) = 0.4*m1r[15]+0.447213595499958*m1r[9]; 
  data->AEM_S(27,15) = 0.4*m1r[14]+0.447213595499958*m1r[10]; 
  data->AEM_S(28,8) = 0.5*m1r[12]; 
  data->AEM_S(28,9) = 0.4472135954999579*m1r[9]; 
  data->AEM_S(28,10) = 0.5000000000000001*m1r[14]; 
  data->AEM_S(28,11) = 0.4472135954999579*m1r[11]; 
  data->AEM_S(28,12) = 0.31943828249997*m1r[12]+0.5*m1r[8]; 
  data->AEM_S(28,14) = 0.31943828249997*m1r[14]+0.5000000000000001*m1r[10]; 
  data->AEM_S(28,15) = 0.4472135954999579*m1r[15]; 
  data->AEM_S(29,8) = 0.5*m1r[13]; 
  data->AEM_S(29,9) = 0.5000000000000001*m1r[15]; 
  data->AEM_S(29,10) = 0.4472135954999579*m1r[10]; 
  data->AEM_S(29,11) = 0.4472135954999579*m1r[11]; 
  data->AEM_S(29,13) = 0.31943828249997*m1r[13]+0.5*m1r[8]; 
  data->AEM_S(29,14) = 0.4472135954999579*m1r[14]; 
  data->AEM_S(29,15) = 0.31943828249997*m1r[15]+0.5000000000000001*m1r[9]; 
  data->AEM_S(30,8) = 0.5*m1r[14]; 
  data->AEM_S(30,9) = 0.447213595499958*m1r[11]; 
  data->AEM_S(30,10) = 0.5000000000000001*m1r[12]; 
  data->AEM_S(30,11) = 0.4*m1r[15]+0.447213595499958*m1r[9]; 
  data->AEM_S(30,12) = 0.31943828249997*m1r[14]+0.5000000000000001*m1r[10]; 
  data->AEM_S(30,13) = 0.4472135954999579*m1r[14]; 
  data->AEM_S(30,14) = 0.4472135954999579*m1r[13]+0.31943828249997*m1r[12]+0.5*m1r[8]; 
  data->AEM_S(30,15) = 0.4*m1r[11]; 
  data->AEM_S(31,8) = 0.5*m1r[15]; 
  data->AEM_S(31,9) = 0.5000000000000001*m1r[13]; 
  data->AEM_S(31,10) = 0.447213595499958*m1r[11]; 
  data->AEM_S(31,11) = 0.4*m1r[14]+0.447213595499958*m1r[10]; 
  data->AEM_S(31,12) = 0.4472135954999579*m1r[15]; 
  data->AEM_S(31,13) = 0.31943828249997*m1r[15]+0.5000000000000001*m1r[9]; 
  data->AEM_S(31,14) = 0.4*m1r[11]; 
  data->AEM_S(31,15) = 0.31943828249997*m1r[13]+0.4472135954999579*m1r[12]+0.5*m1r[8]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(16,16) = 0.5*m0r[0]; 
  data->AEM_S(16,17) = 0.5*m0r[1]; 
  data->AEM_S(16,18) = 0.5*m0r[2]; 
  data->AEM_S(16,19) = 0.5*m0r[3]; 
  data->AEM_S(16,20) = 0.5*m0r[4]; 
  data->AEM_S(16,21) = 0.5*m0r[5]; 
  data->AEM_S(16,22) = 0.5*m0r[6]; 
  data->AEM_S(16,23) = 0.5*m0r[7]; 
  data->AEM_S(17,16) = 0.5*m0r[1]; 
  data->AEM_S(17,17) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(17,18) = 0.5*m0r[3]; 
  data->AEM_S(17,19) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  data->AEM_S(17,20) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(17,21) = 0.5000000000000001*m0r[7]; 
  data->AEM_S(17,22) = 0.447213595499958*m0r[3]; 
  data->AEM_S(17,23) = 0.5000000000000001*m0r[5]; 
  data->AEM_S(18,16) = 0.5*m0r[2]; 
  data->AEM_S(18,17) = 0.5*m0r[3]; 
  data->AEM_S(18,18) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(18,19) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  data->AEM_S(18,20) = 0.5000000000000001*m0r[6]; 
  data->AEM_S(18,21) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(18,22) = 0.5000000000000001*m0r[4]; 
  data->AEM_S(18,23) = 0.447213595499958*m0r[3]; 
  data->AEM_S(19,16) = 0.5*m0r[3]; 
  data->AEM_S(19,17) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  data->AEM_S(19,18) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  data->AEM_S(19,19) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(19,20) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(19,21) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(19,22) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  data->AEM_S(19,23) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  data->AEM_S(20,16) = 0.5*m0r[4]; 
  data->AEM_S(20,17) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(20,18) = 0.5000000000000001*m0r[6]; 
  data->AEM_S(20,19) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(20,20) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(20,22) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  data->AEM_S(20,23) = 0.4472135954999579*m0r[7]; 
  data->AEM_S(21,16) = 0.5*m0r[5]; 
  data->AEM_S(21,17) = 0.5000000000000001*m0r[7]; 
  data->AEM_S(21,18) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(21,19) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(21,21) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(21,22) = 0.4472135954999579*m0r[6]; 
  data->AEM_S(21,23) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  data->AEM_S(22,16) = 0.5*m0r[6]; 
  data->AEM_S(22,17) = 0.447213595499958*m0r[3]; 
  data->AEM_S(22,18) = 0.5000000000000001*m0r[4]; 
  data->AEM_S(22,19) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  data->AEM_S(22,20) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  data->AEM_S(22,21) = 0.4472135954999579*m0r[6]; 
  data->AEM_S(22,22) = 0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(22,23) = 0.4*m0r[3]; 
  data->AEM_S(23,16) = 0.5*m0r[7]; 
  data->AEM_S(23,17) = 0.5000000000000001*m0r[5]; 
  data->AEM_S(23,18) = 0.447213595499958*m0r[3]; 
  data->AEM_S(23,19) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  data->AEM_S(23,20) = 0.4472135954999579*m0r[7]; 
  data->AEM_S(23,21) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  data->AEM_S(23,22) = 0.4*m0r[3]; 
  data->AEM_S(23,23) = 0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(16,24) = -0.5*cM[16]; 
  data->AEM_S(16,25) = -0.5*cM[17]; 
  data->AEM_S(16,26) = -0.5*cM[18]; 
  data->AEM_S(16,27) = -0.5*cM[19]; 
  data->AEM_S(16,28) = -0.5*cM[20]; 
  data->AEM_S(16,29) = -0.5*cM[21]; 
  data->AEM_S(16,30) = -0.5*cM[22]; 
  data->AEM_S(16,31) = -0.5*cM[23]; 
  data->AEM_S(17,24) = -0.5*cM[17]; 
  data->AEM_S(17,25) = (-0.4472135954999579*cM[20])-0.5*cM[16]; 
  data->AEM_S(17,26) = -0.5*cM[19]; 
  data->AEM_S(17,27) = (-0.447213595499958*cM[22])-0.5*cM[18]; 
  data->AEM_S(17,28) = -0.4472135954999579*cM[17]; 
  data->AEM_S(17,29) = -0.5000000000000001*cM[23]; 
  data->AEM_S(17,30) = -0.447213595499958*cM[19]; 
  data->AEM_S(17,31) = -0.5000000000000001*cM[21]; 
  data->AEM_S(18,24) = -0.5*cM[18]; 
  data->AEM_S(18,25) = -0.5*cM[19]; 
  data->AEM_S(18,26) = (-0.4472135954999579*cM[21])-0.5*cM[16]; 
  data->AEM_S(18,27) = (-0.447213595499958*cM[23])-0.5*cM[17]; 
  data->AEM_S(18,28) = -0.5000000000000001*cM[22]; 
  data->AEM_S(18,29) = -0.4472135954999579*cM[18]; 
  data->AEM_S(18,30) = -0.5000000000000001*cM[20]; 
  data->AEM_S(18,31) = -0.447213595499958*cM[19]; 
  data->AEM_S(19,24) = -0.5*cM[19]; 
  data->AEM_S(19,25) = (-0.447213595499958*cM[22])-0.5*cM[18]; 
  data->AEM_S(19,26) = (-0.447213595499958*cM[23])-0.5*cM[17]; 
  data->AEM_S(19,27) = (-0.4472135954999579*cM[21])-0.4472135954999579*cM[20]-0.5*cM[16]; 
  data->AEM_S(19,28) = -0.4472135954999579*cM[19]; 
  data->AEM_S(19,29) = -0.4472135954999579*cM[19]; 
  data->AEM_S(19,30) = (-0.4*cM[23])-0.447213595499958*cM[17]; 
  data->AEM_S(19,31) = (-0.4*cM[22])-0.447213595499958*cM[18]; 
  data->AEM_S(20,24) = -0.5*cM[20]; 
  data->AEM_S(20,25) = -0.4472135954999579*cM[17]; 
  data->AEM_S(20,26) = -0.5000000000000001*cM[22]; 
  data->AEM_S(20,27) = -0.4472135954999579*cM[19]; 
  data->AEM_S(20,28) = (-0.31943828249997*cM[20])-0.5*cM[16]; 
  data->AEM_S(20,30) = (-0.31943828249997*cM[22])-0.5000000000000001*cM[18]; 
  data->AEM_S(20,31) = -0.4472135954999579*cM[23]; 
  data->AEM_S(21,24) = -0.5*cM[21]; 
  data->AEM_S(21,25) = -0.5000000000000001*cM[23]; 
  data->AEM_S(21,26) = -0.4472135954999579*cM[18]; 
  data->AEM_S(21,27) = -0.4472135954999579*cM[19]; 
  data->AEM_S(21,29) = (-0.31943828249997*cM[21])-0.5*cM[16]; 
  data->AEM_S(21,30) = -0.4472135954999579*cM[22]; 
  data->AEM_S(21,31) = (-0.31943828249997*cM[23])-0.5000000000000001*cM[17]; 
  data->AEM_S(22,24) = -0.5*cM[22]; 
  data->AEM_S(22,25) = -0.447213595499958*cM[19]; 
  data->AEM_S(22,26) = -0.5000000000000001*cM[20]; 
  data->AEM_S(22,27) = (-0.4*cM[23])-0.447213595499958*cM[17]; 
  data->AEM_S(22,28) = (-0.31943828249997*cM[22])-0.5000000000000001*cM[18]; 
  data->AEM_S(22,29) = -0.4472135954999579*cM[22]; 
  data->AEM_S(22,30) = (-0.4472135954999579*cM[21])-0.31943828249997*cM[20]-0.5*cM[16]; 
  data->AEM_S(22,31) = -0.4*cM[19]; 
  data->AEM_S(23,24) = -0.5*cM[23]; 
  data->AEM_S(23,25) = -0.5000000000000001*cM[21]; 
  data->AEM_S(23,26) = -0.447213595499958*cM[19]; 
  data->AEM_S(23,27) = (-0.4*cM[22])-0.447213595499958*cM[18]; 
  data->AEM_S(23,28) = -0.4472135954999579*cM[23]; 
  data->AEM_S(23,29) = (-0.31943828249997*cM[23])-0.5000000000000001*cM[17]; 
  data->AEM_S(23,30) = -0.4*cM[19]; 
  data->AEM_S(23,31) = (-0.31943828249997*cM[21])-0.4472135954999579*cM[20]-0.5*cM[16]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(24,16) = 0.5*m1r[16]; 
  data->AEM_S(24,17) = 0.5*m1r[17]; 
  data->AEM_S(24,18) = 0.5*m1r[18]; 
  data->AEM_S(24,19) = 0.5*m1r[19]; 
  data->AEM_S(24,20) = 0.5*m1r[20]; 
  data->AEM_S(24,21) = 0.5*m1r[21]; 
  data->AEM_S(24,22) = 0.5*m1r[22]; 
  data->AEM_S(24,23) = 0.5*m1r[23]; 
  data->AEM_S(25,16) = 0.5*m1r[17]; 
  data->AEM_S(25,17) = 0.4472135954999579*m1r[20]+0.5*m1r[16]; 
  data->AEM_S(25,18) = 0.5*m1r[19]; 
  data->AEM_S(25,19) = 0.447213595499958*m1r[22]+0.5*m1r[18]; 
  data->AEM_S(25,20) = 0.4472135954999579*m1r[17]; 
  data->AEM_S(25,21) = 0.5000000000000001*m1r[23]; 
  data->AEM_S(25,22) = 0.447213595499958*m1r[19]; 
  data->AEM_S(25,23) = 0.5000000000000001*m1r[21]; 
  data->AEM_S(26,16) = 0.5*m1r[18]; 
  data->AEM_S(26,17) = 0.5*m1r[19]; 
  data->AEM_S(26,18) = 0.4472135954999579*m1r[21]+0.5*m1r[16]; 
  data->AEM_S(26,19) = 0.447213595499958*m1r[23]+0.5*m1r[17]; 
  data->AEM_S(26,20) = 0.5000000000000001*m1r[22]; 
  data->AEM_S(26,21) = 0.4472135954999579*m1r[18]; 
  data->AEM_S(26,22) = 0.5000000000000001*m1r[20]; 
  data->AEM_S(26,23) = 0.447213595499958*m1r[19]; 
  data->AEM_S(27,16) = 0.5*m1r[19]; 
  data->AEM_S(27,17) = 0.447213595499958*m1r[22]+0.5*m1r[18]; 
  data->AEM_S(27,18) = 0.447213595499958*m1r[23]+0.5*m1r[17]; 
  data->AEM_S(27,19) = 0.4472135954999579*m1r[21]+0.4472135954999579*m1r[20]+0.5*m1r[16]; 
  data->AEM_S(27,20) = 0.4472135954999579*m1r[19]; 
  data->AEM_S(27,21) = 0.4472135954999579*m1r[19]; 
  data->AEM_S(27,22) = 0.4*m1r[23]+0.447213595499958*m1r[17]; 
  data->AEM_S(27,23) = 0.4*m1r[22]+0.447213595499958*m1r[18]; 
  data->AEM_S(28,16) = 0.5*m1r[20]; 
  data->AEM_S(28,17) = 0.4472135954999579*m1r[17]; 
  data->AEM_S(28,18) = 0.5000000000000001*m1r[22]; 
  data->AEM_S(28,19) = 0.4472135954999579*m1r[19]; 
  data->AEM_S(28,20) = 0.31943828249997*m1r[20]+0.5*m1r[16]; 
  data->AEM_S(28,22) = 0.31943828249997*m1r[22]+0.5000000000000001*m1r[18]; 
  data->AEM_S(28,23) = 0.4472135954999579*m1r[23]; 
  data->AEM_S(29,16) = 0.5*m1r[21]; 
  data->AEM_S(29,17) = 0.5000000000000001*m1r[23]; 
  data->AEM_S(29,18) = 0.4472135954999579*m1r[18]; 
  data->AEM_S(29,19) = 0.4472135954999579*m1r[19]; 
  data->AEM_S(29,21) = 0.31943828249997*m1r[21]+0.5*m1r[16]; 
  data->AEM_S(29,22) = 0.4472135954999579*m1r[22]; 
  data->AEM_S(29,23) = 0.31943828249997*m1r[23]+0.5000000000000001*m1r[17]; 
  data->AEM_S(30,16) = 0.5*m1r[22]; 
  data->AEM_S(30,17) = 0.447213595499958*m1r[19]; 
  data->AEM_S(30,18) = 0.5000000000000001*m1r[20]; 
  data->AEM_S(30,19) = 0.4*m1r[23]+0.447213595499958*m1r[17]; 
  data->AEM_S(30,20) = 0.31943828249997*m1r[22]+0.5000000000000001*m1r[18]; 
  data->AEM_S(30,21) = 0.4472135954999579*m1r[22]; 
  data->AEM_S(30,22) = 0.4472135954999579*m1r[21]+0.31943828249997*m1r[20]+0.5*m1r[16]; 
  data->AEM_S(30,23) = 0.4*m1r[19]; 
  data->AEM_S(31,16) = 0.5*m1r[23]; 
  data->AEM_S(31,17) = 0.5000000000000001*m1r[21]; 
  data->AEM_S(31,18) = 0.447213595499958*m1r[19]; 
  data->AEM_S(31,19) = 0.4*m1r[22]+0.447213595499958*m1r[18]; 
  data->AEM_S(31,20) = 0.4472135954999579*m1r[23]; 
  data->AEM_S(31,21) = 0.31943828249997*m1r[23]+0.5000000000000001*m1r[17]; 
  data->AEM_S(31,22) = 0.4*m1r[19]; 
  data->AEM_S(31,23) = 0.31943828249997*m1r[21]+0.4472135954999579*m1r[20]+0.5*m1r[16]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(24,24) = 1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(24,25) = 1.5*m0r[1]-0.5*cE[1]; 
  data->AEM_S(24,26) = 1.5*m0r[2]-0.5*cE[2]; 
  data->AEM_S(24,27) = 1.5*m0r[3]-0.5*cE[3]; 
  data->AEM_S(24,28) = 1.5*m0r[4]-0.5*cE[4]; 
  data->AEM_S(24,29) = 1.5*m0r[5]-0.5*cE[5]; 
  data->AEM_S(24,30) = 1.5*m0r[6]-0.5*cE[6]; 
  data->AEM_S(24,31) = 1.5*m0r[7]-0.5*cE[7]; 
  data->AEM_S(25,24) = 1.5*m0r[1]-0.5*cE[1]; 
  data->AEM_S(25,25) = 1.341640786499874*m0r[4]-0.4472135954999579*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(25,26) = 1.5*m0r[3]-0.5*cE[3]; 
  data->AEM_S(25,27) = 1.341640786499874*m0r[6]-0.447213595499958*cE[6]+1.5*m0r[2]-0.5*cE[2]; 
  data->AEM_S(25,28) = 1.341640786499874*m0r[1]-0.4472135954999579*cE[1]; 
  data->AEM_S(25,29) = 1.5*m0r[7]-0.5000000000000001*cE[7]; 
  data->AEM_S(25,30) = 1.341640786499874*m0r[3]-0.447213595499958*cE[3]; 
  data->AEM_S(25,31) = 1.5*m0r[5]-0.5000000000000001*cE[5]; 
  data->AEM_S(26,24) = 1.5*m0r[2]-0.5*cE[2]; 
  data->AEM_S(26,25) = 1.5*m0r[3]-0.5*cE[3]; 
  data->AEM_S(26,26) = 1.341640786499874*m0r[5]-0.4472135954999579*cE[5]+1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(26,27) = 1.341640786499874*m0r[7]-0.447213595499958*cE[7]+1.5*m0r[1]-0.5*cE[1]; 
  data->AEM_S(26,28) = 1.5*m0r[6]-0.5000000000000001*cE[6]; 
  data->AEM_S(26,29) = 1.341640786499874*m0r[2]-0.4472135954999579*cE[2]; 
  data->AEM_S(26,30) = 1.5*m0r[4]-0.5000000000000001*cE[4]; 
  data->AEM_S(26,31) = 1.341640786499874*m0r[3]-0.447213595499958*cE[3]; 
  data->AEM_S(27,24) = 1.5*m0r[3]-0.5*cE[3]; 
  data->AEM_S(27,25) = 1.341640786499874*m0r[6]-0.447213595499958*cE[6]+1.5*m0r[2]-0.5*cE[2]; 
  data->AEM_S(27,26) = 1.341640786499874*m0r[7]-0.447213595499958*cE[7]+1.5*m0r[1]-0.5*cE[1]; 
  data->AEM_S(27,27) = 1.341640786499874*m0r[5]-0.4472135954999579*cE[5]+1.341640786499874*m0r[4]-0.4472135954999579*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(27,28) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  data->AEM_S(27,29) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  data->AEM_S(27,30) = 1.2*m0r[7]-0.4*cE[7]+1.341640786499874*m0r[1]-0.447213595499958*cE[1]; 
  data->AEM_S(27,31) = 1.2*m0r[6]-0.4*cE[6]+1.341640786499874*m0r[2]-0.447213595499958*cE[2]; 
  data->AEM_S(28,24) = 1.5*m0r[4]-0.5*cE[4]; 
  data->AEM_S(28,25) = 1.341640786499874*m0r[1]-0.4472135954999579*cE[1]; 
  data->AEM_S(28,26) = 1.5*m0r[6]-0.5000000000000001*cE[6]; 
  data->AEM_S(28,27) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  data->AEM_S(28,28) = 0.9583148474999099*m0r[4]-0.31943828249997*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(28,30) = 0.9583148474999099*m0r[6]-0.31943828249997*cE[6]+1.5*m0r[2]-0.5000000000000001*cE[2]; 
  data->AEM_S(28,31) = 1.341640786499874*m0r[7]-0.4472135954999579*cE[7]; 
  data->AEM_S(29,24) = 1.5*m0r[5]-0.5*cE[5]; 
  data->AEM_S(29,25) = 1.5*m0r[7]-0.5000000000000001*cE[7]; 
  data->AEM_S(29,26) = 1.341640786499874*m0r[2]-0.4472135954999579*cE[2]; 
  data->AEM_S(29,27) = 1.341640786499874*m0r[3]-0.4472135954999579*cE[3]; 
  data->AEM_S(29,29) = 0.9583148474999099*m0r[5]-0.31943828249997*cE[5]+1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(29,30) = 1.341640786499874*m0r[6]-0.4472135954999579*cE[6]; 
  data->AEM_S(29,31) = 0.9583148474999099*m0r[7]-0.31943828249997*cE[7]+1.5*m0r[1]-0.5000000000000001*cE[1]; 
  data->AEM_S(30,24) = 1.5*m0r[6]-0.5*cE[6]; 
  data->AEM_S(30,25) = 1.341640786499874*m0r[3]-0.447213595499958*cE[3]; 
  data->AEM_S(30,26) = 1.5*m0r[4]-0.5000000000000001*cE[4]; 
  data->AEM_S(30,27) = 1.2*m0r[7]-0.4*cE[7]+1.341640786499874*m0r[1]-0.447213595499958*cE[1]; 
  data->AEM_S(30,28) = 0.9583148474999099*m0r[6]-0.31943828249997*cE[6]+1.5*m0r[2]-0.5000000000000001*cE[2]; 
  data->AEM_S(30,29) = 1.341640786499874*m0r[6]-0.4472135954999579*cE[6]; 
  data->AEM_S(30,30) = 1.341640786499874*m0r[5]-0.4472135954999579*cE[5]+0.9583148474999099*m0r[4]-0.31943828249997*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
  data->AEM_S(30,31) = 1.2*m0r[3]-0.4*cE[3]; 
  data->AEM_S(31,24) = 1.5*m0r[7]-0.5*cE[7]; 
  data->AEM_S(31,25) = 1.5*m0r[5]-0.5000000000000001*cE[5]; 
  data->AEM_S(31,26) = 1.341640786499874*m0r[3]-0.447213595499958*cE[3]; 
  data->AEM_S(31,27) = 1.2*m0r[6]-0.4*cE[6]+1.341640786499874*m0r[2]-0.447213595499958*cE[2]; 
  data->AEM_S(31,28) = 1.341640786499874*m0r[7]-0.4472135954999579*cE[7]; 
  data->AEM_S(31,29) = 0.9583148474999099*m0r[7]-0.31943828249997*cE[7]+1.5*m0r[1]-0.5000000000000001*cE[1]; 
  data->AEM_S(31,30) = 1.2*m0r[3]-0.4*cE[3]; 
  data->AEM_S(31,31) = 0.9583148474999099*m0r[5]-0.31943828249997*cE[5]+1.341640786499874*m0r[4]-0.4472135954999579*cE[4]+1.5*m0r[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<8,16>(0,8).setZero(); 
  data->AEM_S.block<16,8>(8,0).setZero(); 
  data->AEM_S.block<8,8>(8,16).setZero(); 
  data->AEM_S.block<8,8>(16,8).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m1r[12],m1r[13],m1r[14],m1r[15],m1r[16],m1r[17],m1r[18],m1r[19],m1r[20],m1r[21],m1r[22],m1r[23],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5],m2r[6],m2r[7]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,24,1) = data->u_S.segment<24>(0); 
 
  Eigen::Map<VectorXd>(vtSq,8,1) = data->u_S.segment<8>(24); 
 
} 
 
