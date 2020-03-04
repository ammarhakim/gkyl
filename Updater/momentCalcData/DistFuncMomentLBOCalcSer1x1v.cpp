#include <DistFuncMomentCalcModDecl.h> 
void VmM0Star1x1vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[2]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 1.0*(wr[1]-wl[1]); 
 
  out[0] += ((-0.408248290463863*fr[2])+0.408248290463863*fl[2]+0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*dS; 
  out[1] += ((-0.408248290463863*fr[3])+0.408248290463863*fl[3]+0.3535533905932737*fr[1]+0.3535533905932737*fl[1])*dS; 
 
} 
 
void VmM1iM2Star1x1vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2) 
{ 
  // w[2]:    Cell-center coordinates. 
  // dxv[2]:  Cell length in each direciton. 
  // intFac:  for gyrokinetics (not used in Vm). 
  // m_:      mass (not used in Vm). 
  // Bmag[2]: Magnetic field magnitude (not used in Vm). 
  // f:       Distribution function. 
  // outM1i:  Contribution to M_1^star from this cell. 
  // outM2:   Contribution to M_2^star from this cell. 
 
  const double volFact = 0.5*dxv[1]; 
  double wvSq[1]; 
  wvSq[0]  = w[1]*w[1]; 
  double dvSq[1]; 
  dvSq[0] = dxv[1]*dxv[1]; 
  double tempM0[2]; 

  tempM0[0] = 1.414213562373095*f[0]*volFact; 
  tempM0[1] = 1.414213562373095*f[1]*volFact; 

  outM1i[0] += tempM0[0]*w[1]; 
  outM1i[1] += tempM0[1]*w[1]; 

  outM2[0] += 0.408248290463863*dxv[1]*w[1]*f[2]*volFact+tempM0[0]*wvSq[0]; 
  outM2[1] += 0.408248290463863*dxv[1]*w[1]*f[3]*volFact+tempM0[1]*wvSq[0]; 
 
} 
void VmBoundaryIntegral1x1vSer_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*intFac; 
 
  if (atLower) {
 
    out[0] += ((-1.58113883008419*fIn[4])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.58113883008419*fIn[4]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x1vSer_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*intFac; 
 
  if (atLower) {
 
    out[0] += ((-1.58113883008419*fIn[4])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += (1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.58113883008419*fIn[4]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += (1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS*vBoundary; 
 
  }
 
} 
 
