#include <DistFuncMomentCalcModDecl.h> 
void VmM0Star2x2vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[3]*(wr[2]-wl[2]); 
 
  out[0] += ((-0.5773502691896258*fr[3])+0.5773502691896258*fl[3]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += ((-0.5773502691896258*fr[6])+0.5773502691896258*fl[6]+0.5*fr[1]+0.5*fl[1])*dS; 
  out[2] += ((-0.5773502691896258*fr[7])+0.5773502691896258*fl[7]+0.5*fr[2]+0.5*fl[2])*dS; 
  out[3] += ((-0.5773502691896258*fr[11])+0.5773502691896258*fl[11]+0.5*fr[5]+0.5*fl[5])*dS; 
 
} 
 
void VmM0Star2x2vSer_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[2]*(wr[3]-wl[3]); 
 
  out[0] += ((-0.5773502691896258*fr[4])+0.5773502691896258*fl[4]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += ((-0.5773502691896258*fr[8])+0.5773502691896258*fl[8]+0.5*fr[1]+0.5*fl[1])*dS; 
  out[2] += ((-0.5773502691896258*fr[9])+0.5773502691896258*fl[9]+0.5*fr[2]+0.5*fl[2])*dS; 
  out[3] += ((-0.5773502691896258*fr[12])+0.5773502691896258*fl[12]+0.5*fr[5]+0.5*fl[5])*dS; 
 
} 
 
void VmM1iM2Star2x2vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2) 
{ 
  // w[4]:    Cell-center coordinates. 
  // dxv[4]:  Cell length in each direciton. 
  // intFac:  for gyrokinetics (not used in Vm). 
  // m_:      mass (not used in Vm). 
  // Bmag[4]: Magnetic field magnitude (not used in Vm). 
  // f:       Distribution function. 
  // outM1i:  Contribution to M_1^star from this cell. 
  // outM2:   Contribution to M_2^star from this cell. 
 
  const double volFact = 0.25*dxv[2]*dxv[3]; 
  double wvSq[2]; 
  wvSq[0]  = w[2]*w[2]; 
  wvSq[1]  = w[3]*w[3]; 
  double dvSq[2]; 
  dvSq[0] = dxv[2]*dxv[2]; 
  dvSq[1] = dxv[3]*dxv[3]; 
  double tempM0[4]; 

  tempM0[0] = 2.0*f[0]*volFact; 
  tempM0[1] = 2.0*f[1]*volFact; 
  tempM0[2] = 2.0*f[2]*volFact; 
  tempM0[3] = 2.0*f[5]*volFact; 

  outM1i[0] += tempM0[0]*w[2]; 
  outM1i[1] += tempM0[1]*w[2]; 
  outM1i[2] += tempM0[2]*w[2]; 
  outM1i[3] += w[2]*tempM0[3]; 
  outM1i[4] += tempM0[0]*w[3]; 
  outM1i[5] += tempM0[1]*w[3]; 
  outM1i[6] += tempM0[2]*w[3]; 
  outM1i[7] += tempM0[3]*w[3]; 

  outM2[0] += (0.5773502691896258*dxv[3]*w[3]*f[4]+0.5773502691896258*dxv[2]*w[2]*f[3])*volFact+tempM0[0]*(wvSq[1]+wvSq[0]); 
  outM2[1] += (0.5773502691896258*dxv[3]*w[3]*f[8]+0.5773502691896258*dxv[2]*w[2]*f[6])*volFact+tempM0[1]*(wvSq[1]+wvSq[0]); 
  outM2[2] += (0.5773502691896258*dxv[3]*w[3]*f[9]+0.5773502691896258*dxv[2]*w[2]*f[7])*volFact+tempM0[2]*(wvSq[1]+wvSq[0]); 
  outM2[3] += (0.5773502691896258*dxv[3]*w[3]*f[12]+0.5773502691896258*dxv[2]*w[2]*f[11])*volFact+tempM0[3]*(wvSq[1]+wvSq[0]); 
 
} 
void VmBoundaryIntegral2x2vSer_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[17]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[16])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]-1.0*fIn[5])*dS; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[16]+1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[6]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[7]+fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]+fIn[5])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[17]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[16])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[6]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[7]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[11]-1.0*fIn[5])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[16]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[6]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[7]+fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[11]+fIn[5])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_F_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[17]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[4] += (1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[5] += (1.732050807568877*fIn[8]-1.0*fIn[1])*dS; 
    out[6] += (1.732050807568877*fIn[9]-1.0*fIn[2])*dS; 
    out[7] += (1.732050807568877*fIn[12]-1.0*fIn[5])*dS; 
 
  } else {
 
    out[4] += (1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[5] += (1.732050807568877*fIn[8]+fIn[1])*dS; 
    out[6] += (1.732050807568877*fIn[9]+fIn[2])*dS; 
    out[7] += (1.732050807568877*fIn[12]+fIn[5])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[17]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[8]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[9]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[12]-1.0*fIn[5])*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[8]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[9]+fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[12]+fIn[5])*dS*vBoundary; 
 
  }
 
} 
 
