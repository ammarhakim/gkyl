#include <DistFuncMomentCalcModDecl.h> 
#include <cmath> 
void GkM0Star1x2vMax_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[3]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[2]*intFac*(wr[1]-wl[1]); 
 
  out[0] += ((-0.5773502691896258*fr[2])+0.5773502691896258*fl[2]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += (0.5*fr[1]+0.5*fl[1])*dS; 
 
} 
 
void GkM1iM2Star1x2vMax(const double *w, const double *dxv, const double intFac, const double m_, const double *Bmag, const double *f, double *outM1i, double *outM2) 
{ 
  // w[3]:    Cell-center coordinates. 
  // dxv[3]:  Cell length in each direciton. 
  // intFac:  =2pi/m for gyrokinetics. 
  // m_:      mass. 
  // Bmag[2]: Magnetic field magnitude. 
  // f:       Distribution function. 
  // outM1i:  Contribution to M_1^star from this cell. 
  // outM2:   Contribution to M_2^star from this cell. 
 
  const double volFact = intFac*0.25*dxv[1]*dxv[2]; 
  double wvSq[2]; 
  wvSq[0]  = w[1]*w[1]; 
  wvSq[1]  = w[2]*w[2]; 
  double dvSq[2]; 
  dvSq[0] = dxv[1]*dxv[1]; 
  dvSq[1] = dxv[2]*dxv[2]; 
 
  outM1i[0] += 2.0*f[0]*w[1]*volFact; 
  outM1i[1] += 2.0*f[1]*w[1]*volFact; 
 
  double tmp[2]; 
  tmp[0] = 0.5773502691896258*dxv[2]*f[3]+2.0*f[0]*w[2]; 
  tmp[1] = 2.0*f[1]*w[2]; 
 
  outM2[0] += ((1.414213562373095*Bmag[1]*tmp[1]+1.414213562373095*Bmag[0]*tmp[0])/m_+0.5773502691896258*dxv[1]*w[1]*f[2]+2.0*f[0]*wvSq[0])*volFact; 
  outM2[1] += ((1.414213562373095*Bmag[0]*tmp[1]+1.414213562373095*tmp[0]*Bmag[1])/m_+2.0*f[1]*wvSq[0])*volFact; 
 
} 
void GkBoundaryIntegral1x2vMax_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[4])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += -1.0*fIn[1]*dS; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[4]+1.732050807568877*fIn[2]+fIn[0])*dS; 
    out[1] += fIn[1]*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vMax_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[10]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[8])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
    out[2] += -1.0*fIn[7]*dS; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[4]+fIn[1])*dS; 
    out[2] += fIn[7]*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vMax_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[20]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[18]-2.23606797749979*fIn[8]+1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += ((-2.23606797749979*fIn[12])+1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS; 
    out[3] += -1.0*fIn[17]*dS; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[18]+2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS; 
    out[1] += (2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS; 
    out[3] += fIn[17]*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vMax_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[4])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += -1.0*fIn[1]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[4]+1.732050807568877*fIn[2]+fIn[0])*dS*vBoundary; 
    out[1] += fIn[1]*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vMax_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[10]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[8])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[4]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += -1.0*fIn[7]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[4]+fIn[1])*dS*vBoundary; 
    out[2] += fIn[7]*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vMax_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[20]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]*intFac; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[18]-2.23606797749979*fIn[8]+1.732050807568877*fIn[2]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[12])+1.732050807568877*fIn[4]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS*vBoundary; 
    out[3] += -1.0*fIn[17]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[18]+2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS*vBoundary; 
    out[3] += fIn[17]*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vMax_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += -1.0*fIn[1]*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += fIn[1]*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vMax_vF_VY_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[10]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[9])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[5]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += -1.0*fIn[7]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[9]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[5]+fIn[1])*dS*vBoundary; 
    out[2] += fIn[7]*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral1x2vMax_vF_VY_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[20]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]*intFac; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[19]-2.23606797749979*fIn[9]+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[15])+1.732050807568877*fIn[5]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]-1.0*fIn[7])*dS*vBoundary; 
    out[3] += -1.0*fIn[17]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[19]+2.23606797749979*fIn[9]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[15]+1.732050807568877*fIn[5]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]+fIn[7])*dS*vBoundary; 
    out[3] += fIn[17]*dS*vBoundary; 
 
  }
 
} 
 
