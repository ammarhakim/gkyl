#include <DistFuncMomentCalcModDecl.h> 
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
  outM1i[2] += w[2]*tempM0[2]; 
  outM1i[3] += w[2]*tempM0[3]; 
  outM1i[4] += tempM0[0]*w[3]; 
  outM1i[5] += tempM0[1]*w[3]; 
  outM1i[6] += tempM0[2]*w[3]; 
  outM1i[7] += w[3]*tempM0[3]; 
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
 
