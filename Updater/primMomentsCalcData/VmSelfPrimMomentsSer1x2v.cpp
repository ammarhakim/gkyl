
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
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
void VmBoundaryIntegral1x2vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[2]+fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[4]+fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[20]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[8])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += ((-2.23606797749979*fIn[12])+1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS; 
    out[1] += (2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[18]-2.23606797749979*fIn[8]+1.732050807568877*fIn[2]-1.0*fIn[0])*dS; 
    out[1] += (2.645751311064591*fIn[24]-2.23606797749979*fIn[12]+1.732050807568877*fIn[4]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS; 
    out[3] += (1.732050807568877*fIn[23]-1.0*fIn[17])*dS; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[18]+2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS; 
    out[1] += (2.645751311064591*fIn[24]+2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS; 
    out[3] += (1.732050807568877*fIn[23]+fIn[17])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[2]-1.0*fIn[0])*dS*vBoundary+(0.8660254037844386*dxv[1]*fIn[2]-0.5*fIn[0]*dxv[1])*dS; 
    out[1] += (1.732050807568877*fIn[4]-1.0*fIn[1])*dS*vBoundary+(0.8660254037844386*dxv[1]*fIn[4]-0.5*dxv[1]*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[2]+fIn[0])*dS*vBoundary+((-0.8660254037844386*dxv[1]*fIn[2])-0.5*fIn[0]*dxv[1])*dS; 
    out[1] += (1.732050807568877*fIn[4]+fIn[1])*dS*vBoundary+((-0.8660254037844386*dxv[1]*fIn[4])-0.5*dxv[1]*fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[20]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[8])+1.732050807568877*fIn[2]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[12])+1.732050807568877*fIn[4]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[18]-2.23606797749979*fIn[8]+1.732050807568877*fIn[2]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[24]-2.23606797749979*fIn[12]+1.732050807568877*fIn[4]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[11]-1.0*fIn[7])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[23]-1.0*fIn[17])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[18]+2.23606797749979*fIn[8]+1.732050807568877*fIn[2]+fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[24]+2.23606797749979*fIn[12]+1.732050807568877*fIn[4]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[11]+fIn[7])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[23]+fIn[17])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  if (atLower) {
 
    out[2] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[3] += (1.732050807568877*fIn[5]-1.0*fIn[1])*dS; 
 
  } else {
 
    out[2] += (1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[3] += (1.732050807568877*fIn[5]+fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[20]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  if (atLower) {
 
    out[3] += ((-2.23606797749979*fIn[9])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[4] += ((-2.23606797749979*fIn[15])+1.732050807568877*fIn[5]-1.0*fIn[1])*dS; 
    out[5] += (1.732050807568877*fIn[13]-1.0*fIn[7])*dS; 
 
  } else {
 
    out[3] += (2.23606797749979*fIn[9]+1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[4] += (2.23606797749979*fIn[15]+1.732050807568877*fIn[5]+fIn[1])*dS; 
    out[5] += (1.732050807568877*fIn[13]+fIn[7])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  if (atLower) {
 
    out[4] += (2.645751311064591*fIn[19]-2.23606797749979*fIn[9]+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[5] += (2.645751311064591*fIn[27]-2.23606797749979*fIn[15]+1.732050807568877*fIn[5]-1.0*fIn[1])*dS; 
    out[6] += (1.732050807568877*fIn[13]-1.0*fIn[7])*dS; 
    out[7] += (1.732050807568877*fIn[25]-1.0*fIn[17])*dS; 
 
  } else {
 
    out[4] += (2.645751311064591*fIn[19]+2.23606797749979*fIn[9]+1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[5] += (2.645751311064591*fIn[27]+2.23606797749979*fIn[15]+1.732050807568877*fIn[5]+fIn[1])*dS; 
    out[6] += (1.732050807568877*fIn[13]+fIn[7])*dS; 
    out[7] += (1.732050807568877*fIn[25]+fIn[17])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary+(0.8660254037844386*dxv[2]*fIn[3]-0.5*fIn[0]*dxv[2])*dS; 
    out[1] += (1.732050807568877*fIn[5]-1.0*fIn[1])*dS*vBoundary+(0.8660254037844386*dxv[2]*fIn[5]-0.5*fIn[1]*dxv[2])*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary+((-0.8660254037844386*dxv[2]*fIn[3])-0.5*fIn[0]*dxv[2])*dS; 
    out[1] += (1.732050807568877*fIn[5]+fIn[1])*dS*vBoundary+((-0.8660254037844386*dxv[2]*fIn[5])-0.5*fIn[1]*dxv[2])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[20]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[9])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[15])+1.732050807568877*fIn[5]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]-1.0*fIn[7])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[9]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[15]+1.732050807568877*fIn[5]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]+fIn[7])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x2vSer_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[3]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[19]-2.23606797749979*fIn[9]+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[27]-2.23606797749979*fIn[15]+1.732050807568877*fIn[5]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]-1.0*fIn[7])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[25]-1.0*fIn[17])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[19]+2.23606797749979*fIn[9]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[27]+2.23606797749979*fIn[15]+1.732050807568877*fIn[5]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]+fIn[7])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[25]+fIn[17])*dS*vBoundary; 
 
  }
 
} 
 
void VmSelfPrimMoments1x2vSer_P1(binOpData_t* data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
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
  double m1r[4]; 
  double m0Sr[2]; 
  double m1Sr[4]; 
  double m2Sr[2]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = m1[2]; 
    m1r[3] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(6,6); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,0) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,1) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,4) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,5) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,4) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,5) = -0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(4,0) = 0.7071067811865475*m1Sr[0]; 
  data->AEM_S(4,1) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(5,0) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(5,1) = 0.7071067811865475*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(2,2) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(2,3) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(3,2) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(3,3) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(2,4) = -0.7071067811865475*cM[2]; 
  data->AEM_S(2,5) = -0.7071067811865475*cM[3]; 
  data->AEM_S(3,4) = -0.7071067811865475*cM[3]; 
  data->AEM_S(3,5) = -0.7071067811865475*cM[2]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(4,2) = 0.7071067811865475*m1Sr[2]; 
  data->AEM_S(4,3) = 0.7071067811865475*m1Sr[3]; 
  data->AEM_S(5,2) = 0.7071067811865475*m1Sr[3]; 
  data->AEM_S(5,3) = 0.7071067811865475*m1Sr[2]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(4,4) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(4,5) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(5,4) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(5,5) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<2,2>(0,2).setZero(); 
  data->AEM_S.block<2,2>(2,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m2Sr[0],m2Sr[1]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,4,1) = data->u_S.segment<4>(0); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = data->u_S.segment<2>(4); 
 
} 
 
void VmSelfPrimMoments1x2vSer_P2(binOpData_t* data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
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
  double m1r[6]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(9,9); 
  
  
 
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
  data->AEM_S(0,6) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,7) = -0.7071067811865475*cM[1]; 
  data->AEM_S(0,8) = -0.7071067811865475*cM[2]; 
  data->AEM_S(1,6) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,7) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(1,8) = -0.6324555320336759*cM[1]; 
  data->AEM_S(2,6) = -0.7071067811865475*cM[2]; 
  data->AEM_S(2,7) = -0.6324555320336759*cM[1]; 
  data->AEM_S(2,8) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(6,0) = 0.7071067811865475*m1r[0]; 
  data->AEM_S(6,1) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(6,2) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(7,0) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(7,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(7,2) = 0.6324555320336759*m1r[1]; 
  data->AEM_S(8,0) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(8,1) = 0.6324555320336759*m1r[1]; 
  data->AEM_S(8,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
 
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
  data->AEM_S(3,6) = -0.7071067811865475*cM[3]; 
  data->AEM_S(3,7) = -0.7071067811865475*cM[4]; 
  data->AEM_S(3,8) = -0.7071067811865475*cM[5]; 
  data->AEM_S(4,6) = -0.7071067811865475*cM[4]; 
  data->AEM_S(4,7) = (-0.6324555320336759*cM[5])-0.7071067811865475*cM[3]; 
  data->AEM_S(4,8) = -0.6324555320336759*cM[4]; 
  data->AEM_S(5,6) = -0.7071067811865475*cM[5]; 
  data->AEM_S(5,7) = -0.6324555320336759*cM[4]; 
  data->AEM_S(5,8) = (-0.4517539514526256*cM[5])-0.7071067811865475*cM[3]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(6,3) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(6,4) = 0.7071067811865475*m1r[4]; 
  data->AEM_S(6,5) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(7,3) = 0.7071067811865475*m1r[4]; 
  data->AEM_S(7,4) = 0.6324555320336759*m1r[5]+0.7071067811865475*m1r[3]; 
  data->AEM_S(7,5) = 0.6324555320336759*m1r[4]; 
  data->AEM_S(8,3) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(8,4) = 0.6324555320336759*m1r[4]; 
  data->AEM_S(8,5) = 0.4517539514526256*m1r[5]+0.7071067811865475*m1r[3]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(6,6) = 1.414213562373095*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(6,7) = 1.414213562373095*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(6,8) = 1.414213562373095*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(7,6) = 1.414213562373095*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(7,7) = 1.264911064067352*m0r[2]-0.6324555320336759*cE[2]+1.414213562373095*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(7,8) = 1.264911064067352*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(8,6) = 1.414213562373095*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(8,7) = 1.264911064067352*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(8,8) = 0.9035079029052515*m0r[2]-0.4517539514526256*cE[2]+1.414213562373095*m0r[0]-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<3,3>(0,3).setZero(); 
  data->AEM_S.block<3,3>(3,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m2r[0],m2r[1],m2r[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,6,1) = data->u_S.segment<6>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = data->u_S.segment<3>(6); 
 
} 
 
void VmSelfPrimMoments1x2vSer_P3(binOpData_t* data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
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
  double m1r[8]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(12,12); 
  
  
 
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
  data->AEM_S(0,8) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,9) = -0.7071067811865475*cM[1]; 
  data->AEM_S(0,10) = -0.7071067811865475*cM[2]; 
  data->AEM_S(0,11) = -0.7071067811865475*cM[3]; 
  data->AEM_S(1,8) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,9) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(1,10) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  data->AEM_S(1,11) = -0.6210590034081186*cM[2]; 
  data->AEM_S(2,8) = -0.7071067811865475*cM[2]; 
  data->AEM_S(2,9) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  data->AEM_S(2,10) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(2,11) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  data->AEM_S(3,8) = -0.7071067811865475*cM[3]; 
  data->AEM_S(3,9) = -0.6210590034081186*cM[2]; 
  data->AEM_S(3,10) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  data->AEM_S(3,11) = (-0.421637021355784*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(8,0) = 0.7071067811865475*m1r[0]; 
  data->AEM_S(8,1) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(8,2) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(8,3) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(9,0) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(9,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(9,2) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  data->AEM_S(9,3) = 0.6210590034081186*m1r[2]; 
  data->AEM_S(10,0) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(10,1) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  data->AEM_S(10,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(10,3) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  data->AEM_S(11,0) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(11,1) = 0.6210590034081186*m1r[2]; 
  data->AEM_S(11,2) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  data->AEM_S(11,3) = 0.421637021355784*m1r[2]+0.7071067811865475*m1r[0]; 
 
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
  data->AEM_S(4,8) = -0.7071067811865475*cM[4]; 
  data->AEM_S(4,9) = -0.7071067811865475*cM[5]; 
  data->AEM_S(4,10) = -0.7071067811865475*cM[6]; 
  data->AEM_S(4,11) = -0.7071067811865475*cM[7]; 
  data->AEM_S(5,8) = -0.7071067811865475*cM[5]; 
  data->AEM_S(5,9) = (-0.6324555320336759*cM[6])-0.7071067811865475*cM[4]; 
  data->AEM_S(5,10) = (-0.6210590034081186*cM[7])-0.6324555320336759*cM[5]; 
  data->AEM_S(5,11) = -0.6210590034081186*cM[6]; 
  data->AEM_S(6,8) = -0.7071067811865475*cM[6]; 
  data->AEM_S(6,9) = (-0.6210590034081186*cM[7])-0.6324555320336759*cM[5]; 
  data->AEM_S(6,10) = (-0.4517539514526256*cM[6])-0.7071067811865475*cM[4]; 
  data->AEM_S(6,11) = (-0.421637021355784*cM[7])-0.6210590034081186*cM[5]; 
  data->AEM_S(7,8) = -0.7071067811865475*cM[7]; 
  data->AEM_S(7,9) = -0.6210590034081186*cM[6]; 
  data->AEM_S(7,10) = (-0.421637021355784*cM[7])-0.6210590034081186*cM[5]; 
  data->AEM_S(7,11) = (-0.421637021355784*cM[6])-0.7071067811865475*cM[4]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(8,4) = 0.7071067811865475*m1r[4]; 
  data->AEM_S(8,5) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(8,6) = 0.7071067811865475*m1r[6]; 
  data->AEM_S(8,7) = 0.7071067811865475*m1r[7]; 
  data->AEM_S(9,4) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(9,5) = 0.6324555320336759*m1r[6]+0.7071067811865475*m1r[4]; 
  data->AEM_S(9,6) = 0.6210590034081186*m1r[7]+0.6324555320336759*m1r[5]; 
  data->AEM_S(9,7) = 0.6210590034081186*m1r[6]; 
  data->AEM_S(10,4) = 0.7071067811865475*m1r[6]; 
  data->AEM_S(10,5) = 0.6210590034081186*m1r[7]+0.6324555320336759*m1r[5]; 
  data->AEM_S(10,6) = 0.4517539514526256*m1r[6]+0.7071067811865475*m1r[4]; 
  data->AEM_S(10,7) = 0.421637021355784*m1r[7]+0.6210590034081186*m1r[5]; 
  data->AEM_S(11,4) = 0.7071067811865475*m1r[7]; 
  data->AEM_S(11,5) = 0.6210590034081186*m1r[6]; 
  data->AEM_S(11,6) = 0.421637021355784*m1r[7]+0.6210590034081186*m1r[5]; 
  data->AEM_S(11,7) = 0.421637021355784*m1r[6]+0.7071067811865475*m1r[4]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(8,8) = 1.414213562373095*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(8,9) = 1.414213562373095*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(8,10) = 1.414213562373095*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(8,11) = 1.414213562373095*m0r[3]-0.7071067811865475*cE[3]; 
  data->AEM_S(9,8) = 1.414213562373095*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(9,9) = 1.264911064067352*m0r[2]-0.6324555320336759*cE[2]+1.414213562373095*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(9,10) = 1.242118006816237*m0r[3]-0.6210590034081186*cE[3]+1.264911064067352*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(9,11) = 1.242118006816237*m0r[2]-0.6210590034081186*cE[2]; 
  data->AEM_S(10,8) = 1.414213562373095*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(10,9) = 1.242118006816237*m0r[3]-0.6210590034081186*cE[3]+1.264911064067352*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(10,10) = 0.9035079029052515*m0r[2]-0.4517539514526256*cE[2]+1.414213562373095*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(10,11) = 0.8432740427115681*m0r[3]-0.421637021355784*cE[3]+1.242118006816237*m0r[1]-0.6210590034081186*cE[1]; 
  data->AEM_S(11,8) = 1.414213562373095*m0r[3]-0.7071067811865475*cE[3]; 
  data->AEM_S(11,9) = 1.242118006816237*m0r[2]-0.6210590034081186*cE[2]; 
  data->AEM_S(11,10) = 0.8432740427115681*m0r[3]-0.421637021355784*cE[3]+1.242118006816237*m0r[1]-0.6210590034081186*cE[1]; 
  data->AEM_S(11,11) = 0.8432740427115681*m0r[2]-0.421637021355784*cE[2]+1.414213562373095*m0r[0]-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<4,4>(0,4).setZero(); 
  data->AEM_S.block<4,4>(4,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m2r[0],m2r[1],m2r[2],m2r[3]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,8,1) = data->u_S.segment<8>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = data->u_S.segment<4>(8); 
 
} 
 
