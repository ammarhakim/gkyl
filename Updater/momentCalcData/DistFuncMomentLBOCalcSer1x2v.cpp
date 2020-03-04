#include <DistFuncMomentCalcModDecl.h> 
void VmM0Star1x2vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[3]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[2]*(wr[1]-wl[1]); 
 
  out[0] += ((-0.5773502691896258*fr[2])+0.5773502691896258*fl[2]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += ((-0.5773502691896258*fr[4])+0.5773502691896258*fl[4]+0.5*fr[1]+0.5*fl[1])*dS; 
 
} 
 
void VmM0Star1x2vSer_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[3]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[1]*(wr[2]-wl[2]); 
 
  out[0] += ((-0.5773502691896258*fr[3])+0.5773502691896258*fl[3]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += ((-0.5773502691896258*fr[5])+0.5773502691896258*fl[5]+0.5*fr[1]+0.5*fl[1])*dS; 
 
} 
 
void VmM1iM2Star1x2vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2) 
{ 
  // w[3]:    Cell-center coordinates. 
  // dxv[3]:  Cell length in each direciton. 
  // intFac:  for gyrokinetics (not used in Vm). 
  // m_:      mass (not used in Vm). 
  // Bmag[2]: Magnetic field magnitude (not used in Vm). 
  // f:       Distribution function. 
  // outM1i:  Contribution to M_1^star from this cell. 
  // outM2:   Contribution to M_2^star from this cell. 
 
  const double volFact = 0.25*dxv[1]*dxv[2]; 
  double wvSq[2]; 
  wvSq[0]  = w[1]*w[1]; 
  wvSq[1]  = w[2]*w[2]; 
  double dvSq[2]; 
  dvSq[0] = dxv[1]*dxv[1]; 
  dvSq[1] = dxv[2]*dxv[2]; 
  double tempM0[2]; 

  tempM0[0] = 2.0*f[0]*volFact; 
  tempM0[1] = 2.0*f[1]*volFact; 

  outM1i[0] += tempM0[0]*w[1]; 
  outM1i[1] += tempM0[1]*w[1]; 
  outM1i[2] += tempM0[0]*w[2]; 
  outM1i[3] += tempM0[1]*w[2]; 

  outM2[0] += (0.5773502691896258*dxv[2]*w[2]*f[3]+0.5773502691896258*dxv[1]*w[1]*f[2])*volFact+tempM0[0]*(wvSq[1]+wvSq[0]); 
  outM2[1] += (0.5773502691896258*dxv[2]*w[2]*f[5]+0.5773502691896258*dxv[1]*w[1]*f[4])*volFact+tempM0[1]*(wvSq[1]+wvSq[0]); 
 
} 
void VmBoundaryIntegral1x2vSer_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[9]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[8])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[4]+fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[9]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[8])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[4]-1.0*fIn[1])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[4]+fIn[1])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_F_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[9]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]*intFac; 
 
  if (atLower) {
 
    out[2] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[3] += (1.732050807568877*fIn[5]-1.0*fIn[1])*dS; 
 
  } else {
 
    out[2] += (1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[3] += (1.732050807568877*fIn[5]+fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[9]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[5]-1.0*fIn[1])*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[5]+fIn[1])*dS*vBoundary; 
 
  }
 
} 
 
