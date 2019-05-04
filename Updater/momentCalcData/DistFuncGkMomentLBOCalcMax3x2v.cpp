#include <DistFuncMomentCalcModDecl.h> 
void GkM0Star3x2vMax_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[5]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[4]*intFac*(wr[3]-wl[3]); 
 
  out[0] += ((-0.5773502691896258*fr[4])+0.5773502691896258*fl[4]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += (0.5*fr[1]+0.5*fl[1])*dS; 
  out[2] += (0.5*fr[2]+0.5*fl[2])*dS; 
  out[3] += (0.5*fr[3]+0.5*fl[3])*dS; 
 
} 
 
void GkM1iM2Star3x2vMax(const double *w, const double *dxv, const double intFac, const double m_, const double *Bmag, const double *f, double *outM1i, double *outM2) 
{ 
  // w[5]:    Cell-center coordinates. 
  // dxv[5]:  Cell length in each direciton. 
  // intFac:  =2pi/m for gyrokinetics. 
  // m_:      mass. 
  // Bmag[4]: Magnetic field magnitude. 
  // f:       Distribution function. 
  // outM1i:  Contribution to M_1^star from this cell. 
  // outM2:   Contribution to M_2^star from this cell. 
 
  const double volFact = intFac*0.25*dxv[3]*dxv[4]; 
  double wvSq[2]; 
  wvSq[0]  = w[3]*w[3]; 
  wvSq[1]  = w[4]*w[4]; 
  double dvSq[2]; 
  dvSq[0] = dxv[3]*dxv[3]; 
  dvSq[1] = dxv[4]*dxv[4]; 
 
  outM1i[0] += 2.0*f[0]*w[3]*volFact; 
  outM1i[1] += 2.0*f[1]*w[3]*volFact; 
  outM1i[2] += 2.0*f[2]*w[3]*volFact; 
  outM1i[3] += 2.0*f[3]*w[3]*volFact; 
 
  double tmp[4]; 
  tmp[0] = 0.5773502691896258*dxv[4]*f[5]+2.0*f[0]*w[4]; 
  tmp[1] = 2.0*f[1]*w[4]; 
  tmp[2] = 2.0*f[2]*w[4]; 
  tmp[3] = 2.0*f[3]*w[4]; 
 
  outM2[0] += ((0.7071067811865474*Bmag[3]*tmp[3]+0.7071067811865474*Bmag[2]*tmp[2]+0.7071067811865474*Bmag[1]*tmp[1]+0.7071067811865474*Bmag[0]*tmp[0])/m_+0.5773502691896258*dxv[3]*w[3]*f[4]+2.0*f[0]*wvSq[0])*volFact; 
  outM2[1] += ((0.7071067811865474*Bmag[0]*tmp[1]+0.7071067811865474*tmp[0]*Bmag[1])/m_+2.0*f[1]*wvSq[0])*volFact; 
  outM2[2] += ((0.7071067811865474*Bmag[0]*tmp[2]+0.7071067811865474*tmp[0]*Bmag[2])/m_+2.0*f[2]*wvSq[0])*volFact; 
  outM2[3] += ((0.7071067811865474*Bmag[0]*tmp[3]+0.7071067811865474*tmp[0]*Bmag[3])/m_+2.0*f[3]*wvSq[0])*volFact; 
 
} 
void GkBoundaryIntegral3x2vMax_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[6]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[1] += -1.0*fIn[1]*dS; 
    out[2] += -1.0*fIn[2]*dS; 
    out[3] += -1.0*fIn[3]*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[1] += fIn[1]*dS; 
    out[2] += fIn[2]*dS; 
    out[3] += fIn[3]*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral3x2vMax_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[21]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[19])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[9]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[10]-1.0*fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]-1.0*fIn[3])*dS; 
    out[4] += -1.0*fIn[6]*dS; 
    out[5] += -1.0*fIn[7]*dS; 
    out[6] += -1.0*fIn[8]*dS; 
    out[7] += -1.0*fIn[16]*dS; 
    out[8] += -1.0*fIn[17]*dS; 
    out[9] += -1.0*fIn[18]*dS; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[19]+1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[9]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[10]+fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]+fIn[3])*dS; 
    out[4] += fIn[6]*dS; 
    out[5] += fIn[7]*dS; 
    out[6] += fIn[8]*dS; 
    out[7] += fIn[16]*dS; 
    out[8] += fIn[17]*dS; 
    out[9] += fIn[18]*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral3x2vMax_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[6]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[4]-0.5*fIn[0]*dxv[3])*dS; 
    out[1] += (-1.0*fIn[1]*dS*vBoundary)-0.5*fIn[1]*dxv[3]*dS; 
    out[2] += (-1.0*fIn[2]*dS*vBoundary)-0.5*fIn[2]*dxv[3]*dS; 
    out[3] += (-1.0*fIn[3]*dS*vBoundary)-0.5*dxv[3]*fIn[3]*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[4])-0.5*fIn[0]*dxv[3])*dS; 
    out[1] += fIn[1]*dS*vBoundary-0.5*fIn[1]*dxv[3]*dS; 
    out[2] += fIn[2]*dS*vBoundary-0.5*fIn[2]*dxv[3]*dS; 
    out[3] += fIn[3]*dS*vBoundary-0.5*dxv[3]*fIn[3]*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral3x2vMax_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[21]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[19])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[9]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[10]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[11]-1.0*fIn[3])*dS*vBoundary; 
    out[4] += -1.0*fIn[6]*dS*vBoundary; 
    out[5] += -1.0*fIn[7]*dS*vBoundary; 
    out[6] += -1.0*fIn[8]*dS*vBoundary; 
    out[7] += -1.0*fIn[16]*dS*vBoundary; 
    out[8] += -1.0*fIn[17]*dS*vBoundary; 
    out[9] += -1.0*fIn[18]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[19]+1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[9]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[10]+fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[11]+fIn[3])*dS*vBoundary; 
    out[4] += fIn[6]*dS*vBoundary; 
    out[5] += fIn[7]*dS*vBoundary; 
    out[6] += fIn[8]*dS*vBoundary; 
    out[7] += fIn[16]*dS*vBoundary; 
    out[8] += fIn[17]*dS*vBoundary; 
    out[9] += fIn[18]*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral3x2vMax_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[6]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[5]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += -1.0*fIn[1]*dS*vBoundary; 
    out[2] += -1.0*fIn[2]*dS*vBoundary; 
    out[3] += -1.0*fIn[3]*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[5]+fIn[0])*dS*vBoundary; 
    out[1] += fIn[1]*dS*vBoundary; 
    out[2] += fIn[2]*dS*vBoundary; 
    out[3] += fIn[3]*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral3x2vMax_vF_VY_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[21]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[20])+1.732050807568877*fIn[5]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[12]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[14]-1.0*fIn[3])*dS*vBoundary; 
    out[4] += -1.0*fIn[6]*dS*vBoundary; 
    out[5] += -1.0*fIn[7]*dS*vBoundary; 
    out[6] += -1.0*fIn[8]*dS*vBoundary; 
    out[7] += -1.0*fIn[16]*dS*vBoundary; 
    out[8] += -1.0*fIn[17]*dS*vBoundary; 
    out[9] += -1.0*fIn[18]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[20]+1.732050807568877*fIn[5]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[12]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]+fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[14]+fIn[3])*dS*vBoundary; 
    out[4] += fIn[6]*dS*vBoundary; 
    out[5] += fIn[7]*dS*vBoundary; 
    out[6] += fIn[8]*dS*vBoundary; 
    out[7] += fIn[16]*dS*vBoundary; 
    out[8] += fIn[17]*dS*vBoundary; 
    out[9] += fIn[18]*dS*vBoundary; 
 
  }
 
} 
 
