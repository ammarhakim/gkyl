#include <DistFuncMomentCalcModDecl.h> 
void VmM0Star1x3vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direction. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[2]*dxvl[3]*(wr[1]-wl[1]); 
 

  out[0] += ((-0.8164965809277261*fr[2])+0.8164965809277261*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += (0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
 
} 
 
void VmM0Star1x3vMax_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direction. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[1]*dxvl[3]*(wr[2]-wl[2]); 
 

  out[0] += ((-0.8164965809277261*fr[3])+0.8164965809277261*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += (0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
 
} 
 
void VmM0Star1x3vMax_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direction. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[1]*dxvl[2]*(wr[3]-wl[3]); 
 

  out[0] += ((-0.8164965809277261*fr[4])+0.8164965809277261*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += (0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
 
} 
 
void VmM0StarNonuniform1x3vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direction. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[2]*dxvl[3]*(wr[1]-wl[1]); 
 
  const double dxvl1R2 = pow(dxvl[1],2);
  const double dxvl1R3 = pow(dxvl[1],3);
  const double dxvr1R2 = pow(dxvr[1],2);
  const double dxvr1R3 = pow(dxvr[1],3);

  out[0] += ((-(17.32050807568877*dxvl1R2*dxvr[1]*fr[2])/(4.242640687119286*dxvr1R3+12.72792206135786*dxvl[1]*dxvr1R2+12.72792206135786*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3))-(10.39230484541326*dxvl1R3*fr[2])/(4.242640687119286*dxvr1R3+12.72792206135786*dxvl[1]*dxvr1R2+12.72792206135786*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)+(10.39230484541326*dxvr1R3*fl[2])/(4.242640687119286*dxvr1R3+12.72792206135786*dxvl[1]*dxvr1R2+12.72792206135786*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)+(17.32050807568877*dxvl[1]*dxvr1R2*fl[2])/(4.242640687119286*dxvr1R3+12.72792206135786*dxvl[1]*dxvr1R2+12.72792206135786*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)+(6.0*fl[0]*dxvr1R3)/(4.242640687119286*dxvr1R3+12.72792206135786*dxvl[1]*dxvr1R2+12.72792206135786*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)+(18.0*fl[0]*dxvl[1]*dxvr1R2)/(4.242640687119286*dxvr1R3+12.72792206135786*dxvl[1]*dxvr1R2+12.72792206135786*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)+(18.0*fr[0]*dxvl1R2*dxvr[1])/(4.242640687119286*dxvr1R3+12.72792206135786*dxvl[1]*dxvr1R2+12.72792206135786*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3)+(6.0*fr[0]*dxvl1R3)/(4.242640687119286*dxvr1R3+12.72792206135786*dxvl[1]*dxvr1R2+12.72792206135786*dxvl1R2*dxvr[1]+4.242640687119286*dxvl1R3))*dS; 
  out[1] += ((4.242640687119286*dxvl1R2*dxvr[1]*fr[1])/(dxvr1R3+3.0*dxvl[1]*dxvr1R2+3.0*dxvl1R2*dxvr[1]+dxvl1R3)+(1.414213562373095*dxvl1R3*fr[1])/(dxvr1R3+3.0*dxvl[1]*dxvr1R2+3.0*dxvl1R2*dxvr[1]+dxvl1R3)+(1.414213562373095*dxvr1R3*fl[1])/(dxvr1R3+3.0*dxvl[1]*dxvr1R2+3.0*dxvl1R2*dxvr[1]+dxvl1R3)+(4.242640687119286*dxvl[1]*dxvr1R2*fl[1])/(dxvr1R3+3.0*dxvl[1]*dxvr1R2+3.0*dxvl1R2*dxvr[1]+dxvl1R3))*dS; 
 
} 
 
void VmM0StarNonuniform1x3vMax_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direction. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[1]*dxvl[3]*(wr[2]-wl[2]); 
 
  const double dxvl2R2 = pow(dxvl[2],2);
  const double dxvl2R3 = pow(dxvl[2],3);
  const double dxvr2R2 = pow(dxvr[2],2);
  const double dxvr2R3 = pow(dxvr[2],3);

  out[0] += ((-(17.32050807568877*dxvl2R2*dxvr[2]*fr[3])/(4.242640687119286*dxvr2R3+12.72792206135786*dxvl[2]*dxvr2R2+12.72792206135786*dxvl2R2*dxvr[2]+4.242640687119286*dxvl2R3))-(10.39230484541326*dxvl2R3*fr[3])/(4.242640687119286*dxvr2R3+12.72792206135786*dxvl[2]*dxvr2R2+12.72792206135786*dxvl2R2*dxvr[2]+4.242640687119286*dxvl2R3)+(10.39230484541326*dxvr2R3*fl[3])/(4.242640687119286*dxvr2R3+12.72792206135786*dxvl[2]*dxvr2R2+12.72792206135786*dxvl2R2*dxvr[2]+4.242640687119286*dxvl2R3)+(17.32050807568877*dxvl[2]*dxvr2R2*fl[3])/(4.242640687119286*dxvr2R3+12.72792206135786*dxvl[2]*dxvr2R2+12.72792206135786*dxvl2R2*dxvr[2]+4.242640687119286*dxvl2R3)+(6.0*fl[0]*dxvr2R3)/(4.242640687119286*dxvr2R3+12.72792206135786*dxvl[2]*dxvr2R2+12.72792206135786*dxvl2R2*dxvr[2]+4.242640687119286*dxvl2R3)+(18.0*fl[0]*dxvl[2]*dxvr2R2)/(4.242640687119286*dxvr2R3+12.72792206135786*dxvl[2]*dxvr2R2+12.72792206135786*dxvl2R2*dxvr[2]+4.242640687119286*dxvl2R3)+(18.0*fr[0]*dxvl2R2*dxvr[2])/(4.242640687119286*dxvr2R3+12.72792206135786*dxvl[2]*dxvr2R2+12.72792206135786*dxvl2R2*dxvr[2]+4.242640687119286*dxvl2R3)+(6.0*fr[0]*dxvl2R3)/(4.242640687119286*dxvr2R3+12.72792206135786*dxvl[2]*dxvr2R2+12.72792206135786*dxvl2R2*dxvr[2]+4.242640687119286*dxvl2R3))*dS; 
  out[1] += ((1.414213562373095*fl[1]*dxvr2R3)/(dxvr2R3+3.0*dxvl[2]*dxvr2R2+3.0*dxvl2R2*dxvr[2]+dxvl2R3)+(4.242640687119286*fl[1]*dxvl[2]*dxvr2R2)/(dxvr2R3+3.0*dxvl[2]*dxvr2R2+3.0*dxvl2R2*dxvr[2]+dxvl2R3)+(4.242640687119286*fr[1]*dxvl2R2*dxvr[2])/(dxvr2R3+3.0*dxvl[2]*dxvr2R2+3.0*dxvl2R2*dxvr[2]+dxvl2R3)+(1.414213562373095*fr[1]*dxvl2R3)/(dxvr2R3+3.0*dxvl[2]*dxvr2R2+3.0*dxvl2R2*dxvr[2]+dxvl2R3))*dS; 
 
} 
 
void VmM0StarNonuniform1x3vMax_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direction. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[1]*dxvl[2]*(wr[3]-wl[3]); 
 
  const double dxvl3R2 = pow(dxvl[3],2);
  const double dxvl3R3 = pow(dxvl[3],3);
  const double dxvr3R2 = pow(dxvr[3],2);
  const double dxvr3R3 = pow(dxvr[3],3);

  out[0] += ((-(8.660254037844386*dxvl3R2*dxvr3R2*fr[4])/(4.242640687119286*dxvr3R3+12.72792206135786*dxvl[3]*dxvr3R2+12.72792206135786*dxvl3R2*dxvr[3]+4.242640687119286*dxvl3R3))-(5.196152422706631*dxvl3R3*dxvr[3]*fr[4])/(4.242640687119286*dxvr3R3+12.72792206135786*dxvl[3]*dxvr3R2+12.72792206135786*dxvl3R2*dxvr[3]+4.242640687119286*dxvl3R3)+(5.196152422706631*dxvl[3]*dxvr3R3*fl[4])/(4.242640687119286*dxvr3R3+12.72792206135786*dxvl[3]*dxvr3R2+12.72792206135786*dxvl3R2*dxvr[3]+4.242640687119286*dxvl3R3)+(8.660254037844386*dxvl3R2*dxvr3R2*fl[4])/(4.242640687119286*dxvr3R3+12.72792206135786*dxvl[3]*dxvr3R2+12.72792206135786*dxvl3R2*dxvr[3]+4.242640687119286*dxvl3R3)+(6.0*fl[0]*dxvr3R3)/(4.242640687119286*dxvr3R3+12.72792206135786*dxvl[3]*dxvr3R2+12.72792206135786*dxvl3R2*dxvr[3]+4.242640687119286*dxvl3R3)+(18.0*fl[0]*dxvl[3]*dxvr3R2)/(4.242640687119286*dxvr3R3+12.72792206135786*dxvl[3]*dxvr3R2+12.72792206135786*dxvl3R2*dxvr[3]+4.242640687119286*dxvl3R3)+(18.0*fr[0]*dxvl3R2*dxvr[3])/(4.242640687119286*dxvr3R3+12.72792206135786*dxvl[3]*dxvr3R2+12.72792206135786*dxvl3R2*dxvr[3]+4.242640687119286*dxvl3R3)+(6.0*fr[0]*dxvl3R3)/(4.242640687119286*dxvr3R3+12.72792206135786*dxvl[3]*dxvr3R2+12.72792206135786*dxvl3R2*dxvr[3]+4.242640687119286*dxvl3R3))*dS; 
  out[1] += ((1.414213562373095*fl[1]*dxvr3R3)/(dxvr3R3+3.0*dxvl[3]*dxvr3R2+3.0*dxvl3R2*dxvr[3]+dxvl3R3)+(4.242640687119286*fl[1]*dxvl[3]*dxvr3R2)/(dxvr3R3+3.0*dxvl[3]*dxvr3R2+3.0*dxvl3R2*dxvr[3]+dxvl3R3)+(4.242640687119286*fr[1]*dxvl3R2*dxvr[3])/(dxvr3R3+3.0*dxvl[3]*dxvr3R2+3.0*dxvl3R2*dxvr[3]+dxvl3R3)+(1.414213562373095*fr[1]*dxvl3R3)/(dxvr3R3+3.0*dxvl[3]*dxvr3R2+3.0*dxvl3R2*dxvr[3]+dxvl3R3))*dS; 
 
} 
 
void VmM1iM2Star1x3vMax(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2) 
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
  outM2[1] += tempM0[1]*(wvSq[2]+wvSq[1]+wvSq[0]); 
 
} 
void VmBoundaryIntegral1x3vMax_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += -1.414213562373095*fIn[1]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[15]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[12])+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
    out[2] += -1.414213562373095*fIn[11]*dS; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
    out[2] += 1.414213562373095*fIn[11]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[35]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (3.741657386773942*fIn[32]-3.16227766016838*fIn[12]+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS; 
    out[1] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS; 
    out[3] += -1.414213562373095*fIn[31]*dS; 
 
  } else {
 
    out[0] += (3.741657386773942*fIn[32]+3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS; 
    out[1] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS; 
    out[2] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS; 
    out[3] += 1.414213562373095*fIn[31]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[1]*fIn[2]-0.7071067811865475*fIn[0]*dxv[1])*dS; 
    out[1] += (-1.414213562373095*fIn[1]*dS*vBoundary)-0.7071067811865475*dxv[1]*fIn[1]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[1]*fIn[2])-0.7071067811865475*fIn[0]*dxv[1])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS*vBoundary-0.7071067811865475*dxv[1]*fIn[1]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[15]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[12])+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += -1.414213562373095*fIn[11]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += 1.414213562373095*fIn[11]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[35]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (3.741657386773942*fIn[32]-3.16227766016838*fIn[12]+2.449489742783178*fIn[2]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[20])+2.449489742783178*fIn[5]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[19]-1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += -1.414213562373095*fIn[31]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.741657386773942*fIn[32]+3.16227766016838*fIn[12]+2.449489742783178*fIn[2]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[20]+2.449489742783178*fIn[5]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[19]+1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += 1.414213562373095*fIn[31]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[2] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[3] += -1.414213562373095*fIn[1]*dS; 
 
  } else {
 
    out[2] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[3] += 1.414213562373095*fIn[1]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[15]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[3] += ((-3.16227766016838*fIn[13])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[4] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
    out[5] += -1.414213562373095*fIn[11]*dS; 
 
  } else {
 
    out[3] += (3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[4] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
    out[5] += 1.414213562373095*fIn[11]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[35]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[4] += (3.741657386773942*fIn[33]-3.16227766016838*fIn[13]+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS; 
    out[5] += ((-3.16227766016838*fIn[23])+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS; 
    out[6] += (2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS; 
    out[7] += -1.414213562373095*fIn[31]*dS; 
 
  } else {
 
    out[4] += (3.741657386773942*fIn[33]+3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS; 
    out[5] += (3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS; 
    out[6] += (2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS; 
    out[7] += 1.414213562373095*fIn[31]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[2]*fIn[3]-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += (-1.414213562373095*fIn[1]*dS*vBoundary)-0.7071067811865475*fIn[1]*dxv[2]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[2]*fIn[3])-0.7071067811865475*fIn[0]*dxv[2])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS*vBoundary-0.7071067811865475*fIn[1]*dxv[2]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[15]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[13])+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += -1.414213562373095*fIn[11]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += 1.414213562373095*fIn[11]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[35]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (3.741657386773942*fIn[33]-3.16227766016838*fIn[13]+2.449489742783178*fIn[3]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[23])+2.449489742783178*fIn[6]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[21]-1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += -1.414213562373095*fIn[31]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.741657386773942*fIn[33]+3.16227766016838*fIn[13]+2.449489742783178*fIn[3]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[23]+2.449489742783178*fIn[6]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[21]+1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += 1.414213562373095*fIn[31]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[4] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[5] += -1.414213562373095*fIn[1]*dS; 
 
  } else {
 
    out[4] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[5] += 1.414213562373095*fIn[1]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[15]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[6] += ((-3.16227766016838*fIn[14])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[7] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
    out[8] += -1.414213562373095*fIn[11]*dS; 
 
  } else {
 
    out[6] += (3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[7] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
    out[8] += 1.414213562373095*fIn[11]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_F_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[35]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[8] += (3.741657386773942*fIn[34]-3.16227766016838*fIn[14]+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS; 
    out[9] += ((-3.16227766016838*fIn[28])+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS; 
    out[11] += -1.414213562373095*fIn[31]*dS; 
 
  } else {
 
    out[8] += (3.741657386773942*fIn[34]+3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS; 
    out[9] += (3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS; 
    out[10] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS; 
    out[11] += 1.414213562373095*fIn[31]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary+(1.224744871391589*dxv[3]*fIn[4]-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += (-1.414213562373095*fIn[1]*dS*vBoundary)-0.7071067811865475*fIn[1]*dxv[3]*dS; 
 
  } else {
 
    out[0] += (2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary+((-1.224744871391589*dxv[3]*fIn[4])-0.7071067811865475*fIn[0]*dxv[3])*dS; 
    out[1] += 1.414213562373095*fIn[1]*dS*vBoundary-0.7071067811865475*fIn[1]*dxv[3]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[15]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[0] += ((-3.16227766016838*fIn[14])+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += -1.414213562373095*fIn[11]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += 1.414213562373095*fIn[11]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x3vMax_vF_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:   cell length in each direction. 
  // fIn[35]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (3.741657386773942*fIn[34]-3.16227766016838*fIn[14]+2.449489742783178*fIn[4]-1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += ((-3.16227766016838*fIn[28])+2.449489742783178*fIn[8]-1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[25]-1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += -1.414213562373095*fIn[31]*dS*vBoundary; 
 
  } else {
 
    out[0] += (3.741657386773942*fIn[34]+3.16227766016838*fIn[14]+2.449489742783178*fIn[4]+1.414213562373095*fIn[0])*dS*vBoundary; 
    out[1] += (3.16227766016838*fIn[28]+2.449489742783178*fIn[8]+1.414213562373095*fIn[1])*dS*vBoundary; 
    out[2] += (2.449489742783178*fIn[25]+1.414213562373095*fIn[11])*dS*vBoundary; 
    out[3] += 1.414213562373095*fIn[31]*dS*vBoundary; 
 
  }
 
} 
 
