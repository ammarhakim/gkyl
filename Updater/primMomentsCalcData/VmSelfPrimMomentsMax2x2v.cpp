#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmM0Star2x2vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[3]*(wr[2]-wl[2]); 
 
  out[0] += ((-0.5773502691896258*fr[3])+0.5773502691896258*fl[3]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += (0.5*fr[1]+0.5*fl[1])*dS; 
  out[2] += (0.5*fr[2]+0.5*fl[2])*dS; 
 
} 
 
void VmM0Star2x2vMax_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[2]*(wr[3]-wl[3]); 
 
  out[0] += ((-0.5773502691896258*fr[4])+0.5773502691896258*fl[4]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += (0.5*fr[1]+0.5*fl[1])*dS; 
  out[2] += (0.5*fr[2]+0.5*fl[2])*dS; 
 
} 
 
void VmM1iM2Star2x2vMax(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2) 
{ 
  // w[4]:    Cell-center coordinates. 
  // dxv[4]:  Cell length in each direciton. 
  // intFac:  for gyrokinetics (not used in Vm). 
  // m_:      mass (not used in Vm). 
  // Bmag[3]: Magnetic field magnitude (not used in Vm). 
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
  double tempM0[3]; 

  tempM0[0] = 2.0*f[0]*volFact; 
  tempM0[1] = 2.0*f[1]*volFact; 
  tempM0[2] = 2.0*f[2]*volFact; 

  outM1i[0] += tempM0[0]*w[2]; 
  outM1i[1] += tempM0[1]*w[2]; 
  outM1i[2] += tempM0[2]*w[2]; 
  outM1i[3] += tempM0[0]*w[3]; 
  outM1i[4] += tempM0[1]*w[3]; 
  outM1i[5] += tempM0[2]*w[3]; 

  outM2[0] += (0.5773502691896258*dxv[3]*w[3]*f[4]+0.5773502691896258*dxv[2]*w[2]*f[3])*volFact+tempM0[0]*(wvSq[1]+wvSq[0]); 
  outM2[1] += tempM0[1]*(wvSq[1]+wvSq[0]); 
  outM2[2] += tempM0[2]*(wvSq[1]+wvSq[0]); 
 
} 
void VmBoundaryIntegral2x2vMax_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += -1.0*fIn[1]*dS; 
    out[2] += -1.0*fIn[2]*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[1] += fIn[1]*dS; 
    out[2] += fIn[2]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[15]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[13])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
    out[3] += -1.0*fIn[5]*dS; 
    out[4] += -1.0*fIn[11]*dS; 
    out[5] += -1.0*fIn[12]*dS; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[6]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[7]+fIn[2])*dS; 
    out[3] += fIn[5]*dS; 
    out[4] += fIn[11]*dS; 
    out[5] += fIn[12]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[35]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[33]-2.23606797749979*fIn[13]+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += ((-2.23606797749979*fIn[23])+1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
    out[2] += ((-2.23606797749979*fIn[24])+1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[15]-1.0*fIn[5])*dS; 
    out[4] += (1.732050807568877*fIn[21]-1.0*fIn[11])*dS; 
    out[5] += (1.732050807568877*fIn[22]-1.0*fIn[12])*dS; 
    out[6] += -1.0*fIn[19]*dS; 
    out[7] += -1.0*fIn[20]*dS; 
    out[8] += -1.0*fIn[31]*dS; 
    out[9] += -1.0*fIn[32]*dS; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[33]+2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[1] += (2.23606797749979*fIn[23]+1.732050807568877*fIn[6]+fIn[1])*dS; 
    out[2] += (2.23606797749979*fIn[24]+1.732050807568877*fIn[7]+fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[15]+fIn[5])*dS; 
    out[4] += (1.732050807568877*fIn[21]+fIn[11])*dS; 
    out[5] += (1.732050807568877*fIn[22]+fIn[12])*dS; 
    out[6] += fIn[19]*dS; 
    out[7] += fIn[20]*dS; 
    out[8] += fIn[31]*dS; 
    out[9] += fIn[32]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary+(0.8660254037844386*dxv[2]*fIn[3]-0.5*fIn[0]*dxv[2])*dS; 
    out[1] += (-1.0*fIn[1]*dS*vBoundary)-0.5*fIn[1]*dxv[2]*dS; 
    out[2] += (-1.0*fIn[2]*dS*vBoundary)-0.5*dxv[2]*fIn[2]*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary+((-0.8660254037844386*dxv[2]*fIn[3])-0.5*fIn[0]*dxv[2])*dS; 
    out[1] += fIn[1]*dS*vBoundary-0.5*fIn[1]*dxv[2]*dS; 
    out[2] += fIn[2]*dS*vBoundary-0.5*dxv[2]*fIn[2]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[15]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[13])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[6]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[7]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += -1.0*fIn[5]*dS*vBoundary; 
    out[4] += -1.0*fIn[11]*dS*vBoundary; 
    out[5] += -1.0*fIn[12]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[6]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[7]+fIn[2])*dS*vBoundary; 
    out[3] += fIn[5]*dS*vBoundary; 
    out[4] += fIn[11]*dS*vBoundary; 
    out[5] += fIn[12]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[35]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[33]-2.23606797749979*fIn[13]+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[23])+1.732050807568877*fIn[6]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += ((-2.23606797749979*fIn[24])+1.732050807568877*fIn[7]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[15]-1.0*fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[21]-1.0*fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[22]-1.0*fIn[12])*dS*vBoundary; 
    out[6] += -1.0*fIn[19]*dS*vBoundary; 
    out[7] += -1.0*fIn[20]*dS*vBoundary; 
    out[8] += -1.0*fIn[31]*dS*vBoundary; 
    out[9] += -1.0*fIn[32]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[33]+2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[23]+1.732050807568877*fIn[6]+fIn[1])*dS*vBoundary; 
    out[2] += (2.23606797749979*fIn[24]+1.732050807568877*fIn[7]+fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[15]+fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[21]+fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[22]+fIn[12])*dS*vBoundary; 
    out[6] += fIn[19]*dS*vBoundary; 
    out[7] += fIn[20]*dS*vBoundary; 
    out[8] += fIn[31]*dS*vBoundary; 
    out[9] += fIn[32]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[3] += (1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[4] += -1.0*fIn[1]*dS; 
    out[5] += -1.0*fIn[2]*dS; 
 
  } else {
 
    out[3] += (1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[4] += fIn[1]*dS; 
    out[5] += fIn[2]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[15]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[6] += ((-2.23606797749979*fIn[14])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[7] += (1.732050807568877*fIn[8]-1.0*fIn[1])*dS; 
    out[8] += (1.732050807568877*fIn[9]-1.0*fIn[2])*dS; 
    out[9] += -1.0*fIn[5]*dS; 
    out[10] += -1.0*fIn[11]*dS; 
    out[11] += -1.0*fIn[12]*dS; 
 
  } else {
 
    out[6] += (2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[7] += (1.732050807568877*fIn[8]+fIn[1])*dS; 
    out[8] += (1.732050807568877*fIn[9]+fIn[2])*dS; 
    out[9] += fIn[5]*dS; 
    out[10] += fIn[11]*dS; 
    out[11] += fIn[12]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[35]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[10] += (2.645751311064591*fIn[34]-2.23606797749979*fIn[14]+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[11] += ((-2.23606797749979*fIn[28])+1.732050807568877*fIn[8]-1.0*fIn[1])*dS; 
    out[12] += ((-2.23606797749979*fIn[29])+1.732050807568877*fIn[9]-1.0*fIn[2])*dS; 
    out[13] += (1.732050807568877*fIn[16]-1.0*fIn[5])*dS; 
    out[14] += (1.732050807568877*fIn[25]-1.0*fIn[11])*dS; 
    out[15] += (1.732050807568877*fIn[26]-1.0*fIn[12])*dS; 
    out[16] += -1.0*fIn[19]*dS; 
    out[17] += -1.0*fIn[20]*dS; 
    out[18] += -1.0*fIn[31]*dS; 
    out[19] += -1.0*fIn[32]*dS; 
 
  } else {
 
    out[10] += (2.645751311064591*fIn[34]+2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[11] += (2.23606797749979*fIn[28]+1.732050807568877*fIn[8]+fIn[1])*dS; 
    out[12] += (2.23606797749979*fIn[29]+1.732050807568877*fIn[9]+fIn[2])*dS; 
    out[13] += (1.732050807568877*fIn[16]+fIn[5])*dS; 
    out[14] += (1.732050807568877*fIn[25]+fIn[11])*dS; 
    out[15] += (1.732050807568877*fIn[26]+fIn[12])*dS; 
    out[16] += fIn[19]*dS; 
    out[17] += fIn[20]*dS; 
    out[18] += fIn[31]*dS; 
    out[19] += fIn[32]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[5]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[4]-0.5*fIn[0]*dxv[3])*dS; 
    out[1] += (-1.0*fIn[1]*dS*vBoundary)-0.5*fIn[1]*dxv[3]*dS; 
    out[2] += (-1.0*fIn[2]*dS*vBoundary)-0.5*fIn[2]*dxv[3]*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[4])-0.5*fIn[0]*dxv[3])*dS; 
    out[1] += fIn[1]*dS*vBoundary-0.5*fIn[1]*dxv[3]*dS; 
    out[2] += fIn[2]*dS*vBoundary-0.5*fIn[2]*dxv[3]*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[15]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[14])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[8]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[9]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += -1.0*fIn[5]*dS*vBoundary; 
    out[4] += -1.0*fIn[11]*dS*vBoundary; 
    out[5] += -1.0*fIn[12]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[8]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[9]+fIn[2])*dS*vBoundary; 
    out[3] += fIn[5]*dS*vBoundary; 
    out[4] += fIn[11]*dS*vBoundary; 
    out[5] += fIn[12]*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vMax_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[35]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[34]-2.23606797749979*fIn[14]+1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[28])+1.732050807568877*fIn[8]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += ((-2.23606797749979*fIn[29])+1.732050807568877*fIn[9]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[16]-1.0*fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[25]-1.0*fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[26]-1.0*fIn[12])*dS*vBoundary; 
    out[6] += -1.0*fIn[19]*dS*vBoundary; 
    out[7] += -1.0*fIn[20]*dS*vBoundary; 
    out[8] += -1.0*fIn[31]*dS*vBoundary; 
    out[9] += -1.0*fIn[32]*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[34]+2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[28]+1.732050807568877*fIn[8]+fIn[1])*dS*vBoundary; 
    out[2] += (2.23606797749979*fIn[29]+1.732050807568877*fIn[9]+fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[16]+fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[25]+fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[26]+fIn[12])*dS*vBoundary; 
    out[6] += fIn[19]*dS*vBoundary; 
    out[7] += fIn[20]*dS*vBoundary; 
    out[8] += fIn[31]*dS*vBoundary; 
    out[9] += fIn[32]*dS*vBoundary; 
 
  }
 
} 
 
void VmSelfPrimMoments2x2vMax_P1(const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1:       moments of the distribution function. 
  // m0S,m1S,m1S: star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-0.8660254037844386*m0[2])-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0[2])-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0[2])+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0[2])+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[3]; 
  double m1r[6]; 
  double m0Sr[3]; 
  double m1Sr[6]; 
  double m2Sr[3]; 
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
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m0Sr[2] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = 0.0; 
    m1Sr[3] = m1S[3]; 
    m1Sr[4] = 0.0; 
    m1Sr[5] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
    m2Sr[2] = 0.0; 
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
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m0Sr[2] = m0S[2]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = m1S[5]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
    m2Sr[2] = m2S[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(9,9); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(9);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(9);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0r[0]; 
  BigAEM(0,1) = 0.5*m0r[1]; 
  BigAEM(0,2) = 0.5*m0r[2]; 
  BigAEM(1,0) = 0.5*m0r[1]; 
  BigAEM(1,1) = 0.5*m0r[0]; 
  BigAEM(2,0) = 0.5*m0r[2]; 
  BigAEM(2,2) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,6) = -0.5*cM[0]; 
  BigAEM(0,7) = -0.5*cM[1]; 
  BigAEM(0,8) = -0.5*cM[2]; 
  BigAEM(1,6) = -0.5*cM[1]; 
  BigAEM(1,7) = -0.5*cM[0]; 
  BigAEM(2,6) = -0.5*cM[2]; 
  BigAEM(2,8) = -0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(6,0) = 0.5*m1Sr[0]; 
  BigAEM(6,1) = 0.5*m1Sr[1]; 
  BigAEM(6,2) = 0.5*m1Sr[2]; 
  BigAEM(7,0) = 0.5*m1Sr[1]; 
  BigAEM(7,1) = 0.5*m1Sr[0]; 
  BigAEM(8,0) = 0.5*m1Sr[2]; 
  BigAEM(8,2) = 0.5*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(3,3) = 0.5*m0r[0]; 
  BigAEM(3,4) = 0.5*m0r[1]; 
  BigAEM(3,5) = 0.5*m0r[2]; 
  BigAEM(4,3) = 0.5*m0r[1]; 
  BigAEM(4,4) = 0.5*m0r[0]; 
  BigAEM(5,3) = 0.5*m0r[2]; 
  BigAEM(5,5) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(3,6) = -0.5*cM[3]; 
  BigAEM(3,7) = -0.5*cM[4]; 
  BigAEM(3,8) = -0.5*cM[5]; 
  BigAEM(4,6) = -0.5*cM[4]; 
  BigAEM(4,7) = -0.5*cM[3]; 
  BigAEM(5,6) = -0.5*cM[5]; 
  BigAEM(5,8) = -0.5*cM[3]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(6,3) = 0.5*m1Sr[3]; 
  BigAEM(6,4) = 0.5*m1Sr[4]; 
  BigAEM(6,5) = 0.5*m1Sr[5]; 
  BigAEM(7,3) = 0.5*m1Sr[4]; 
  BigAEM(7,4) = 0.5*m1Sr[3]; 
  BigAEM(8,3) = 0.5*m1Sr[5]; 
  BigAEM(8,5) = 0.5*m1Sr[3]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(6,6) = 0.5*m0Sr[0]-0.5*cE[0]; 
  BigAEM(6,7) = 0.5*m0Sr[1]-0.5*cE[1]; 
  BigAEM(6,8) = 0.5*m0Sr[2]-0.5*cE[2]; 
  BigAEM(7,6) = 0.5*m0Sr[1]-0.5*cE[1]; 
  BigAEM(7,7) = 0.5*m0Sr[0]-0.5*cE[0]; 
  BigAEM(8,6) = 0.5*m0Sr[2]-0.5*cE[2]; 
  BigAEM(8,8) = 0.5*m0Sr[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<3,3>(0,3).setZero(); 
  BigAEM.block<3,3>(3,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m2Sr[0],m2Sr[1],m2Sr[2]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,6,1) = xEV.segment<6>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = xEV.segment<3>(6); 
 
} 
 
void VmSelfPrimMoments2x2vMax_P2(const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[6]; 
  double m1r[12]; 
  double m2r[6]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m0r[4] = 0.0; 
    m0r[5] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    m1r[6] = m1[6]; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    m2r[4] = 0.0; 
    m2r[5] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
    m0r[4] = m0[4]; 
    m0r[5] = m0[5]; 
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
    m2r[4] = m2[4]; 
    m2r[5] = m2[5]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(18,18); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(18);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(18);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0r[0]; 
  BigAEM(0,1) = 0.5*m0r[1]; 
  BigAEM(0,2) = 0.5*m0r[2]; 
  BigAEM(0,3) = 0.5*m0r[3]; 
  BigAEM(0,4) = 0.5*m0r[4]; 
  BigAEM(0,5) = 0.5*m0r[5]; 
  BigAEM(1,0) = 0.5*m0r[1]; 
  BigAEM(1,1) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(1,2) = 0.5*m0r[3]; 
  BigAEM(1,3) = 0.5*m0r[2]; 
  BigAEM(1,4) = 0.4472135954999579*m0r[1]; 
  BigAEM(2,0) = 0.5*m0r[2]; 
  BigAEM(2,1) = 0.5*m0r[3]; 
  BigAEM(2,2) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(2,3) = 0.5*m0r[1]; 
  BigAEM(2,5) = 0.4472135954999579*m0r[2]; 
  BigAEM(3,0) = 0.5*m0r[3]; 
  BigAEM(3,1) = 0.5*m0r[2]; 
  BigAEM(3,2) = 0.5*m0r[1]; 
  BigAEM(3,3) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(3,4) = 0.4472135954999579*m0r[3]; 
  BigAEM(3,5) = 0.4472135954999579*m0r[3]; 
  BigAEM(4,0) = 0.5*m0r[4]; 
  BigAEM(4,1) = 0.4472135954999579*m0r[1]; 
  BigAEM(4,3) = 0.4472135954999579*m0r[3]; 
  BigAEM(4,4) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(5,0) = 0.5*m0r[5]; 
  BigAEM(5,2) = 0.4472135954999579*m0r[2]; 
  BigAEM(5,3) = 0.4472135954999579*m0r[3]; 
  BigAEM(5,5) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,12) = -0.5*cM[0]; 
  BigAEM(0,13) = -0.5*cM[1]; 
  BigAEM(0,14) = -0.5*cM[2]; 
  BigAEM(0,15) = -0.5*cM[3]; 
  BigAEM(0,16) = -0.5*cM[4]; 
  BigAEM(0,17) = -0.5*cM[5]; 
  BigAEM(1,12) = -0.5*cM[1]; 
  BigAEM(1,13) = (-0.4472135954999579*cM[4])-0.5*cM[0]; 
  BigAEM(1,14) = -0.5*cM[3]; 
  BigAEM(1,15) = -0.5*cM[2]; 
  BigAEM(1,16) = -0.4472135954999579*cM[1]; 
  BigAEM(2,12) = -0.5*cM[2]; 
  BigAEM(2,13) = -0.5*cM[3]; 
  BigAEM(2,14) = (-0.4472135954999579*cM[5])-0.5*cM[0]; 
  BigAEM(2,15) = -0.5*cM[1]; 
  BigAEM(2,17) = -0.4472135954999579*cM[2]; 
  BigAEM(3,12) = -0.5*cM[3]; 
  BigAEM(3,13) = -0.5*cM[2]; 
  BigAEM(3,14) = -0.5*cM[1]; 
  BigAEM(3,15) = (-0.4472135954999579*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(3,16) = -0.4472135954999579*cM[3]; 
  BigAEM(3,17) = -0.4472135954999579*cM[3]; 
  BigAEM(4,12) = -0.5*cM[4]; 
  BigAEM(4,13) = -0.4472135954999579*cM[1]; 
  BigAEM(4,15) = -0.4472135954999579*cM[3]; 
  BigAEM(4,16) = (-0.31943828249997*cM[4])-0.5*cM[0]; 
  BigAEM(5,12) = -0.5*cM[5]; 
  BigAEM(5,14) = -0.4472135954999579*cM[2]; 
  BigAEM(5,15) = -0.4472135954999579*cM[3]; 
  BigAEM(5,17) = (-0.31943828249997*cM[5])-0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(12,0) = 0.5*m1r[0]; 
  BigAEM(12,1) = 0.5*m1r[1]; 
  BigAEM(12,2) = 0.5*m1r[2]; 
  BigAEM(12,3) = 0.5*m1r[3]; 
  BigAEM(12,4) = 0.5*m1r[4]; 
  BigAEM(12,5) = 0.5*m1r[5]; 
  BigAEM(13,0) = 0.5*m1r[1]; 
  BigAEM(13,1) = 0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(13,2) = 0.5*m1r[3]; 
  BigAEM(13,3) = 0.5*m1r[2]; 
  BigAEM(13,4) = 0.4472135954999579*m1r[1]; 
  BigAEM(14,0) = 0.5*m1r[2]; 
  BigAEM(14,1) = 0.5*m1r[3]; 
  BigAEM(14,2) = 0.4472135954999579*m1r[5]+0.5*m1r[0]; 
  BigAEM(14,3) = 0.5*m1r[1]; 
  BigAEM(14,5) = 0.4472135954999579*m1r[2]; 
  BigAEM(15,0) = 0.5*m1r[3]; 
  BigAEM(15,1) = 0.5*m1r[2]; 
  BigAEM(15,2) = 0.5*m1r[1]; 
  BigAEM(15,3) = 0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(15,4) = 0.4472135954999579*m1r[3]; 
  BigAEM(15,5) = 0.4472135954999579*m1r[3]; 
  BigAEM(16,0) = 0.5*m1r[4]; 
  BigAEM(16,1) = 0.4472135954999579*m1r[1]; 
  BigAEM(16,3) = 0.4472135954999579*m1r[3]; 
  BigAEM(16,4) = 0.31943828249997*m1r[4]+0.5*m1r[0]; 
  BigAEM(17,0) = 0.5*m1r[5]; 
  BigAEM(17,2) = 0.4472135954999579*m1r[2]; 
  BigAEM(17,3) = 0.4472135954999579*m1r[3]; 
  BigAEM(17,5) = 0.31943828249997*m1r[5]+0.5*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(6,6) = 0.5*m0r[0]; 
  BigAEM(6,7) = 0.5*m0r[1]; 
  BigAEM(6,8) = 0.5*m0r[2]; 
  BigAEM(6,9) = 0.5*m0r[3]; 
  BigAEM(6,10) = 0.5*m0r[4]; 
  BigAEM(6,11) = 0.5*m0r[5]; 
  BigAEM(7,6) = 0.5*m0r[1]; 
  BigAEM(7,7) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(7,8) = 0.5*m0r[3]; 
  BigAEM(7,9) = 0.5*m0r[2]; 
  BigAEM(7,10) = 0.4472135954999579*m0r[1]; 
  BigAEM(8,6) = 0.5*m0r[2]; 
  BigAEM(8,7) = 0.5*m0r[3]; 
  BigAEM(8,8) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(8,9) = 0.5*m0r[1]; 
  BigAEM(8,11) = 0.4472135954999579*m0r[2]; 
  BigAEM(9,6) = 0.5*m0r[3]; 
  BigAEM(9,7) = 0.5*m0r[2]; 
  BigAEM(9,8) = 0.5*m0r[1]; 
  BigAEM(9,9) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(9,10) = 0.4472135954999579*m0r[3]; 
  BigAEM(9,11) = 0.4472135954999579*m0r[3]; 
  BigAEM(10,6) = 0.5*m0r[4]; 
  BigAEM(10,7) = 0.4472135954999579*m0r[1]; 
  BigAEM(10,9) = 0.4472135954999579*m0r[3]; 
  BigAEM(10,10) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(11,6) = 0.5*m0r[5]; 
  BigAEM(11,8) = 0.4472135954999579*m0r[2]; 
  BigAEM(11,9) = 0.4472135954999579*m0r[3]; 
  BigAEM(11,11) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(6,12) = -0.5*cM[6]; 
  BigAEM(6,13) = -0.5*cM[7]; 
  BigAEM(6,14) = -0.5*cM[8]; 
  BigAEM(6,15) = -0.5*cM[9]; 
  BigAEM(6,16) = -0.5*cM[10]; 
  BigAEM(6,17) = -0.5*cM[11]; 
  BigAEM(7,12) = -0.5*cM[7]; 
  BigAEM(7,13) = (-0.4472135954999579*cM[10])-0.5*cM[6]; 
  BigAEM(7,14) = -0.5*cM[9]; 
  BigAEM(7,15) = -0.5*cM[8]; 
  BigAEM(7,16) = -0.4472135954999579*cM[7]; 
  BigAEM(8,12) = -0.5*cM[8]; 
  BigAEM(8,13) = -0.5*cM[9]; 
  BigAEM(8,14) = (-0.4472135954999579*cM[11])-0.5*cM[6]; 
  BigAEM(8,15) = -0.5*cM[7]; 
  BigAEM(8,17) = -0.4472135954999579*cM[8]; 
  BigAEM(9,12) = -0.5*cM[9]; 
  BigAEM(9,13) = -0.5*cM[8]; 
  BigAEM(9,14) = -0.5*cM[7]; 
  BigAEM(9,15) = (-0.4472135954999579*cM[11])-0.4472135954999579*cM[10]-0.5*cM[6]; 
  BigAEM(9,16) = -0.4472135954999579*cM[9]; 
  BigAEM(9,17) = -0.4472135954999579*cM[9]; 
  BigAEM(10,12) = -0.5*cM[10]; 
  BigAEM(10,13) = -0.4472135954999579*cM[7]; 
  BigAEM(10,15) = -0.4472135954999579*cM[9]; 
  BigAEM(10,16) = (-0.31943828249997*cM[10])-0.5*cM[6]; 
  BigAEM(11,12) = -0.5*cM[11]; 
  BigAEM(11,14) = -0.4472135954999579*cM[8]; 
  BigAEM(11,15) = -0.4472135954999579*cM[9]; 
  BigAEM(11,17) = (-0.31943828249997*cM[11])-0.5*cM[6]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(12,6) = 0.5*m1r[6]; 
  BigAEM(12,7) = 0.5*m1r[7]; 
  BigAEM(12,8) = 0.5*m1r[8]; 
  BigAEM(12,9) = 0.5*m1r[9]; 
  BigAEM(12,10) = 0.5*m1r[10]; 
  BigAEM(12,11) = 0.5*m1r[11]; 
  BigAEM(13,6) = 0.5*m1r[7]; 
  BigAEM(13,7) = 0.4472135954999579*m1r[10]+0.5*m1r[6]; 
  BigAEM(13,8) = 0.5*m1r[9]; 
  BigAEM(13,9) = 0.5*m1r[8]; 
  BigAEM(13,10) = 0.4472135954999579*m1r[7]; 
  BigAEM(14,6) = 0.5*m1r[8]; 
  BigAEM(14,7) = 0.5*m1r[9]; 
  BigAEM(14,8) = 0.4472135954999579*m1r[11]+0.5*m1r[6]; 
  BigAEM(14,9) = 0.5*m1r[7]; 
  BigAEM(14,11) = 0.4472135954999579*m1r[8]; 
  BigAEM(15,6) = 0.5*m1r[9]; 
  BigAEM(15,7) = 0.5*m1r[8]; 
  BigAEM(15,8) = 0.5*m1r[7]; 
  BigAEM(15,9) = 0.4472135954999579*m1r[11]+0.4472135954999579*m1r[10]+0.5*m1r[6]; 
  BigAEM(15,10) = 0.4472135954999579*m1r[9]; 
  BigAEM(15,11) = 0.4472135954999579*m1r[9]; 
  BigAEM(16,6) = 0.5*m1r[10]; 
  BigAEM(16,7) = 0.4472135954999579*m1r[7]; 
  BigAEM(16,9) = 0.4472135954999579*m1r[9]; 
  BigAEM(16,10) = 0.31943828249997*m1r[10]+0.5*m1r[6]; 
  BigAEM(17,6) = 0.5*m1r[11]; 
  BigAEM(17,8) = 0.4472135954999579*m1r[8]; 
  BigAEM(17,9) = 0.4472135954999579*m1r[9]; 
  BigAEM(17,11) = 0.31943828249997*m1r[11]+0.5*m1r[6]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(12,12) = m0r[0]-0.5*cE[0]; 
  BigAEM(12,13) = m0r[1]-0.5*cE[1]; 
  BigAEM(12,14) = m0r[2]-0.5*cE[2]; 
  BigAEM(12,15) = m0r[3]-0.5*cE[3]; 
  BigAEM(12,16) = m0r[4]-0.5*cE[4]; 
  BigAEM(12,17) = m0r[5]-0.5*cE[5]; 
  BigAEM(13,12) = m0r[1]-0.5*cE[1]; 
  BigAEM(13,13) = 0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(13,14) = m0r[3]-0.5*cE[3]; 
  BigAEM(13,15) = m0r[2]-0.5*cE[2]; 
  BigAEM(13,16) = 0.8944271909999159*m0r[1]-0.4472135954999579*cE[1]; 
  BigAEM(14,12) = m0r[2]-0.5*cE[2]; 
  BigAEM(14,13) = m0r[3]-0.5*cE[3]; 
  BigAEM(14,14) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+m0r[0]-0.5*cE[0]; 
  BigAEM(14,15) = m0r[1]-0.5*cE[1]; 
  BigAEM(14,17) = 0.8944271909999159*m0r[2]-0.4472135954999579*cE[2]; 
  BigAEM(15,12) = m0r[3]-0.5*cE[3]; 
  BigAEM(15,13) = m0r[2]-0.5*cE[2]; 
  BigAEM(15,14) = m0r[1]-0.5*cE[1]; 
  BigAEM(15,15) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(15,16) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(15,17) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(16,12) = m0r[4]-0.5*cE[4]; 
  BigAEM(16,13) = 0.8944271909999159*m0r[1]-0.4472135954999579*cE[1]; 
  BigAEM(16,15) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(16,16) = 0.6388765649999399*m0r[4]-0.31943828249997*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(17,12) = m0r[5]-0.5*cE[5]; 
  BigAEM(17,14) = 0.8944271909999159*m0r[2]-0.4472135954999579*cE[2]; 
  BigAEM(17,15) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(17,17) = 0.6388765649999399*m0r[5]-0.31943828249997*cE[5]+m0r[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<6,6>(0,6).setZero(); 
  BigAEM.block<6,6>(6,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,12,1) = xEV.segment<12>(0); 
 
  Eigen::Map<VectorXd>(vtSq,6,1) = xEV.segment<6>(12); 
 
} 
 
void VmSelfPrimMoments2x2vMax_P3(const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-1.322875655532295*m0[9])-1.322875655532295*m0[8]-1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.322875655532295*m0[9])-1.322875655532295*m0[8]-1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.322875655532295*m0[9])+1.322875655532295*m0[8]+1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.322875655532295*m0[9])+1.322875655532295*m0[8]+1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[10]; 
  double m1r[20]; 
  double m2r[10]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m0r[4] = 0.0; 
    m0r[5] = 0.0; 
    m0r[6] = 0.0; 
    m0r[7] = 0.0; 
    m0r[8] = 0.0; 
    m0r[9] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    m1r[9] = 0.0; 
    m1r[10] = m1[10]; 
    m1r[11] = 0.0; 
    m1r[12] = 0.0; 
    m1r[13] = 0.0; 
    m1r[14] = 0.0; 
    m1r[15] = 0.0; 
    m1r[16] = 0.0; 
    m1r[17] = 0.0; 
    m1r[18] = 0.0; 
    m1r[19] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    m2r[4] = 0.0; 
    m2r[5] = 0.0; 
    m2r[6] = 0.0; 
    m2r[7] = 0.0; 
    m2r[8] = 0.0; 
    m2r[9] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
    m0r[4] = m0[4]; 
    m0r[5] = m0[5]; 
    m0r[6] = m0[6]; 
    m0r[7] = m0[7]; 
    m0r[8] = m0[8]; 
    m0r[9] = m0[9]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
    m2r[4] = m2[4]; 
    m2r[5] = m2[5]; 
    m2r[6] = m2[6]; 
    m2r[7] = m2[7]; 
    m2r[8] = m2[8]; 
    m2r[9] = m2[9]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(30,30); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(30);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(30);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0r[0]; 
  BigAEM(0,1) = 0.5*m0r[1]; 
  BigAEM(0,2) = 0.5*m0r[2]; 
  BigAEM(0,3) = 0.5*m0r[3]; 
  BigAEM(0,4) = 0.5*m0r[4]; 
  BigAEM(0,5) = 0.5*m0r[5]; 
  BigAEM(0,6) = 0.5*m0r[6]; 
  BigAEM(0,7) = 0.5*m0r[7]; 
  BigAEM(0,8) = 0.5*m0r[8]; 
  BigAEM(0,9) = 0.5*m0r[9]; 
  BigAEM(1,0) = 0.5*m0r[1]; 
  BigAEM(1,1) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(1,2) = 0.5*m0r[3]; 
  BigAEM(1,3) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(1,4) = 0.4391550328268398*m0r[8]+0.4472135954999579*m0r[1]; 
  BigAEM(1,5) = 0.5000000000000001*m0r[7]; 
  BigAEM(1,6) = 0.447213595499958*m0r[3]; 
  BigAEM(1,7) = 0.5000000000000001*m0r[5]; 
  BigAEM(1,8) = 0.4391550328268398*m0r[4]; 
  BigAEM(2,0) = 0.5*m0r[2]; 
  BigAEM(2,1) = 0.5*m0r[3]; 
  BigAEM(2,2) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(2,3) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(2,4) = 0.5000000000000001*m0r[6]; 
  BigAEM(2,5) = 0.4391550328268398*m0r[9]+0.4472135954999579*m0r[2]; 
  BigAEM(2,6) = 0.5000000000000001*m0r[4]; 
  BigAEM(2,7) = 0.447213595499958*m0r[3]; 
  BigAEM(2,9) = 0.4391550328268398*m0r[5]; 
  BigAEM(3,0) = 0.5*m0r[3]; 
  BigAEM(3,1) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(3,2) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(3,3) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(3,4) = 0.4472135954999579*m0r[3]; 
  BigAEM(3,5) = 0.4472135954999579*m0r[3]; 
  BigAEM(3,6) = 0.4391550328268399*m0r[8]+0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(3,7) = 0.4391550328268399*m0r[9]+0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(3,8) = 0.4391550328268399*m0r[6]; 
  BigAEM(3,9) = 0.4391550328268399*m0r[7]; 
  BigAEM(4,0) = 0.5*m0r[4]; 
  BigAEM(4,1) = 0.4391550328268398*m0r[8]+0.4472135954999579*m0r[1]; 
  BigAEM(4,2) = 0.5000000000000001*m0r[6]; 
  BigAEM(4,3) = 0.4472135954999579*m0r[3]; 
  BigAEM(4,4) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(4,6) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(4,7) = 0.4472135954999579*m0r[7]; 
  BigAEM(4,8) = 0.2981423969999719*m0r[8]+0.4391550328268398*m0r[1]; 
  BigAEM(5,0) = 0.5*m0r[5]; 
  BigAEM(5,1) = 0.5000000000000001*m0r[7]; 
  BigAEM(5,2) = 0.4391550328268398*m0r[9]+0.4472135954999579*m0r[2]; 
  BigAEM(5,3) = 0.4472135954999579*m0r[3]; 
  BigAEM(5,5) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
  BigAEM(5,6) = 0.4472135954999579*m0r[6]; 
  BigAEM(5,7) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(5,9) = 0.2981423969999719*m0r[9]+0.4391550328268398*m0r[2]; 
  BigAEM(6,0) = 0.5*m0r[6]; 
  BigAEM(6,1) = 0.447213595499958*m0r[3]; 
  BigAEM(6,2) = 0.5000000000000001*m0r[4]; 
  BigAEM(6,3) = 0.4391550328268399*m0r[8]+0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(6,4) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(6,5) = 0.4472135954999579*m0r[6]; 
  BigAEM(6,6) = 0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(6,7) = 0.4*m0r[3]; 
  BigAEM(6,8) = 0.4391550328268399*m0r[3]; 
  BigAEM(7,0) = 0.5*m0r[7]; 
  BigAEM(7,1) = 0.5000000000000001*m0r[5]; 
  BigAEM(7,2) = 0.447213595499958*m0r[3]; 
  BigAEM(7,3) = 0.4391550328268399*m0r[9]+0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(7,4) = 0.4472135954999579*m0r[7]; 
  BigAEM(7,5) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(7,6) = 0.4*m0r[3]; 
  BigAEM(7,7) = 0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(7,9) = 0.4391550328268399*m0r[3]; 
  BigAEM(8,0) = 0.5*m0r[8]; 
  BigAEM(8,1) = 0.4391550328268398*m0r[4]; 
  BigAEM(8,3) = 0.4391550328268399*m0r[6]; 
  BigAEM(8,4) = 0.2981423969999719*m0r[8]+0.4391550328268398*m0r[1]; 
  BigAEM(8,6) = 0.4391550328268399*m0r[3]; 
  BigAEM(8,8) = 0.2981423969999719*m0r[4]+0.5*m0r[0]; 
  BigAEM(9,0) = 0.5*m0r[9]; 
  BigAEM(9,2) = 0.4391550328268398*m0r[5]; 
  BigAEM(9,3) = 0.4391550328268399*m0r[7]; 
  BigAEM(9,5) = 0.2981423969999719*m0r[9]+0.4391550328268398*m0r[2]; 
  BigAEM(9,7) = 0.4391550328268399*m0r[3]; 
  BigAEM(9,9) = 0.2981423969999719*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,20) = -0.5*cM[0]; 
  BigAEM(0,21) = -0.5*cM[1]; 
  BigAEM(0,22) = -0.5*cM[2]; 
  BigAEM(0,23) = -0.5*cM[3]; 
  BigAEM(0,24) = -0.5*cM[4]; 
  BigAEM(0,25) = -0.5*cM[5]; 
  BigAEM(0,26) = -0.5*cM[6]; 
  BigAEM(0,27) = -0.5*cM[7]; 
  BigAEM(0,28) = -0.5*cM[8]; 
  BigAEM(0,29) = -0.5*cM[9]; 
  BigAEM(1,20) = -0.5*cM[1]; 
  BigAEM(1,21) = (-0.4472135954999579*cM[4])-0.5*cM[0]; 
  BigAEM(1,22) = -0.5*cM[3]; 
  BigAEM(1,23) = (-0.447213595499958*cM[6])-0.5*cM[2]; 
  BigAEM(1,24) = (-0.4391550328268398*cM[8])-0.4472135954999579*cM[1]; 
  BigAEM(1,25) = -0.5000000000000001*cM[7]; 
  BigAEM(1,26) = -0.447213595499958*cM[3]; 
  BigAEM(1,27) = -0.5000000000000001*cM[5]; 
  BigAEM(1,28) = -0.4391550328268398*cM[4]; 
  BigAEM(2,20) = -0.5*cM[2]; 
  BigAEM(2,21) = -0.5*cM[3]; 
  BigAEM(2,22) = (-0.4472135954999579*cM[5])-0.5*cM[0]; 
  BigAEM(2,23) = (-0.447213595499958*cM[7])-0.5*cM[1]; 
  BigAEM(2,24) = -0.5000000000000001*cM[6]; 
  BigAEM(2,25) = (-0.4391550328268398*cM[9])-0.4472135954999579*cM[2]; 
  BigAEM(2,26) = -0.5000000000000001*cM[4]; 
  BigAEM(2,27) = -0.447213595499958*cM[3]; 
  BigAEM(2,29) = -0.4391550328268398*cM[5]; 
  BigAEM(3,20) = -0.5*cM[3]; 
  BigAEM(3,21) = (-0.447213595499958*cM[6])-0.5*cM[2]; 
  BigAEM(3,22) = (-0.447213595499958*cM[7])-0.5*cM[1]; 
  BigAEM(3,23) = (-0.4472135954999579*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(3,24) = -0.4472135954999579*cM[3]; 
  BigAEM(3,25) = -0.4472135954999579*cM[3]; 
  BigAEM(3,26) = (-0.4391550328268399*cM[8])-0.4*cM[7]-0.447213595499958*cM[1]; 
  BigAEM(3,27) = (-0.4391550328268399*cM[9])-0.4*cM[6]-0.447213595499958*cM[2]; 
  BigAEM(3,28) = -0.4391550328268399*cM[6]; 
  BigAEM(3,29) = -0.4391550328268399*cM[7]; 
  BigAEM(4,20) = -0.5*cM[4]; 
  BigAEM(4,21) = (-0.4391550328268398*cM[8])-0.4472135954999579*cM[1]; 
  BigAEM(4,22) = -0.5000000000000001*cM[6]; 
  BigAEM(4,23) = -0.4472135954999579*cM[3]; 
  BigAEM(4,24) = (-0.31943828249997*cM[4])-0.5*cM[0]; 
  BigAEM(4,26) = (-0.31943828249997*cM[6])-0.5000000000000001*cM[2]; 
  BigAEM(4,27) = -0.4472135954999579*cM[7]; 
  BigAEM(4,28) = (-0.2981423969999719*cM[8])-0.4391550328268398*cM[1]; 
  BigAEM(5,20) = -0.5*cM[5]; 
  BigAEM(5,21) = -0.5000000000000001*cM[7]; 
  BigAEM(5,22) = (-0.4391550328268398*cM[9])-0.4472135954999579*cM[2]; 
  BigAEM(5,23) = -0.4472135954999579*cM[3]; 
  BigAEM(5,25) = (-0.31943828249997*cM[5])-0.5*cM[0]; 
  BigAEM(5,26) = -0.4472135954999579*cM[6]; 
  BigAEM(5,27) = (-0.31943828249997*cM[7])-0.5000000000000001*cM[1]; 
  BigAEM(5,29) = (-0.2981423969999719*cM[9])-0.4391550328268398*cM[2]; 
  BigAEM(6,20) = -0.5*cM[6]; 
  BigAEM(6,21) = -0.447213595499958*cM[3]; 
  BigAEM(6,22) = -0.5000000000000001*cM[4]; 
  BigAEM(6,23) = (-0.4391550328268399*cM[8])-0.4*cM[7]-0.447213595499958*cM[1]; 
  BigAEM(6,24) = (-0.31943828249997*cM[6])-0.5000000000000001*cM[2]; 
  BigAEM(6,25) = -0.4472135954999579*cM[6]; 
  BigAEM(6,26) = (-0.4472135954999579*cM[5])-0.31943828249997*cM[4]-0.5*cM[0]; 
  BigAEM(6,27) = -0.4*cM[3]; 
  BigAEM(6,28) = -0.4391550328268399*cM[3]; 
  BigAEM(7,20) = -0.5*cM[7]; 
  BigAEM(7,21) = -0.5000000000000001*cM[5]; 
  BigAEM(7,22) = -0.447213595499958*cM[3]; 
  BigAEM(7,23) = (-0.4391550328268399*cM[9])-0.4*cM[6]-0.447213595499958*cM[2]; 
  BigAEM(7,24) = -0.4472135954999579*cM[7]; 
  BigAEM(7,25) = (-0.31943828249997*cM[7])-0.5000000000000001*cM[1]; 
  BigAEM(7,26) = -0.4*cM[3]; 
  BigAEM(7,27) = (-0.31943828249997*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(7,29) = -0.4391550328268399*cM[3]; 
  BigAEM(8,20) = -0.5*cM[8]; 
  BigAEM(8,21) = -0.4391550328268398*cM[4]; 
  BigAEM(8,23) = -0.4391550328268399*cM[6]; 
  BigAEM(8,24) = (-0.2981423969999719*cM[8])-0.4391550328268398*cM[1]; 
  BigAEM(8,26) = -0.4391550328268399*cM[3]; 
  BigAEM(8,28) = (-0.2981423969999719*cM[4])-0.5*cM[0]; 
  BigAEM(9,20) = -0.5*cM[9]; 
  BigAEM(9,22) = -0.4391550328268398*cM[5]; 
  BigAEM(9,23) = -0.4391550328268399*cM[7]; 
  BigAEM(9,25) = (-0.2981423969999719*cM[9])-0.4391550328268398*cM[2]; 
  BigAEM(9,27) = -0.4391550328268399*cM[3]; 
  BigAEM(9,29) = (-0.2981423969999719*cM[5])-0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(20,0) = 0.5*m1r[0]; 
  BigAEM(20,1) = 0.5*m1r[1]; 
  BigAEM(20,2) = 0.5*m1r[2]; 
  BigAEM(20,3) = 0.5*m1r[3]; 
  BigAEM(20,4) = 0.5*m1r[4]; 
  BigAEM(20,5) = 0.5*m1r[5]; 
  BigAEM(20,6) = 0.5*m1r[6]; 
  BigAEM(20,7) = 0.5*m1r[7]; 
  BigAEM(20,8) = 0.5*m1r[8]; 
  BigAEM(20,9) = 0.5*m1r[9]; 
  BigAEM(21,0) = 0.5*m1r[1]; 
  BigAEM(21,1) = 0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(21,2) = 0.5*m1r[3]; 
  BigAEM(21,3) = 0.447213595499958*m1r[6]+0.5*m1r[2]; 
  BigAEM(21,4) = 0.4391550328268398*m1r[8]+0.4472135954999579*m1r[1]; 
  BigAEM(21,5) = 0.5000000000000001*m1r[7]; 
  BigAEM(21,6) = 0.447213595499958*m1r[3]; 
  BigAEM(21,7) = 0.5000000000000001*m1r[5]; 
  BigAEM(21,8) = 0.4391550328268398*m1r[4]; 
  BigAEM(22,0) = 0.5*m1r[2]; 
  BigAEM(22,1) = 0.5*m1r[3]; 
  BigAEM(22,2) = 0.4472135954999579*m1r[5]+0.5*m1r[0]; 
  BigAEM(22,3) = 0.447213595499958*m1r[7]+0.5*m1r[1]; 
  BigAEM(22,4) = 0.5000000000000001*m1r[6]; 
  BigAEM(22,5) = 0.4391550328268398*m1r[9]+0.4472135954999579*m1r[2]; 
  BigAEM(22,6) = 0.5000000000000001*m1r[4]; 
  BigAEM(22,7) = 0.447213595499958*m1r[3]; 
  BigAEM(22,9) = 0.4391550328268398*m1r[5]; 
  BigAEM(23,0) = 0.5*m1r[3]; 
  BigAEM(23,1) = 0.447213595499958*m1r[6]+0.5*m1r[2]; 
  BigAEM(23,2) = 0.447213595499958*m1r[7]+0.5*m1r[1]; 
  BigAEM(23,3) = 0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(23,4) = 0.4472135954999579*m1r[3]; 
  BigAEM(23,5) = 0.4472135954999579*m1r[3]; 
  BigAEM(23,6) = 0.4391550328268399*m1r[8]+0.4*m1r[7]+0.447213595499958*m1r[1]; 
  BigAEM(23,7) = 0.4391550328268399*m1r[9]+0.4*m1r[6]+0.447213595499958*m1r[2]; 
  BigAEM(23,8) = 0.4391550328268399*m1r[6]; 
  BigAEM(23,9) = 0.4391550328268399*m1r[7]; 
  BigAEM(24,0) = 0.5*m1r[4]; 
  BigAEM(24,1) = 0.4391550328268398*m1r[8]+0.4472135954999579*m1r[1]; 
  BigAEM(24,2) = 0.5000000000000001*m1r[6]; 
  BigAEM(24,3) = 0.4472135954999579*m1r[3]; 
  BigAEM(24,4) = 0.31943828249997*m1r[4]+0.5*m1r[0]; 
  BigAEM(24,6) = 0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]; 
  BigAEM(24,7) = 0.4472135954999579*m1r[7]; 
  BigAEM(24,8) = 0.2981423969999719*m1r[8]+0.4391550328268398*m1r[1]; 
  BigAEM(25,0) = 0.5*m1r[5]; 
  BigAEM(25,1) = 0.5000000000000001*m1r[7]; 
  BigAEM(25,2) = 0.4391550328268398*m1r[9]+0.4472135954999579*m1r[2]; 
  BigAEM(25,3) = 0.4472135954999579*m1r[3]; 
  BigAEM(25,5) = 0.31943828249997*m1r[5]+0.5*m1r[0]; 
  BigAEM(25,6) = 0.4472135954999579*m1r[6]; 
  BigAEM(25,7) = 0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]; 
  BigAEM(25,9) = 0.2981423969999719*m1r[9]+0.4391550328268398*m1r[2]; 
  BigAEM(26,0) = 0.5*m1r[6]; 
  BigAEM(26,1) = 0.447213595499958*m1r[3]; 
  BigAEM(26,2) = 0.5000000000000001*m1r[4]; 
  BigAEM(26,3) = 0.4391550328268399*m1r[8]+0.4*m1r[7]+0.447213595499958*m1r[1]; 
  BigAEM(26,4) = 0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]; 
  BigAEM(26,5) = 0.4472135954999579*m1r[6]; 
  BigAEM(26,6) = 0.4472135954999579*m1r[5]+0.31943828249997*m1r[4]+0.5*m1r[0]; 
  BigAEM(26,7) = 0.4*m1r[3]; 
  BigAEM(26,8) = 0.4391550328268399*m1r[3]; 
  BigAEM(27,0) = 0.5*m1r[7]; 
  BigAEM(27,1) = 0.5000000000000001*m1r[5]; 
  BigAEM(27,2) = 0.447213595499958*m1r[3]; 
  BigAEM(27,3) = 0.4391550328268399*m1r[9]+0.4*m1r[6]+0.447213595499958*m1r[2]; 
  BigAEM(27,4) = 0.4472135954999579*m1r[7]; 
  BigAEM(27,5) = 0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]; 
  BigAEM(27,6) = 0.4*m1r[3]; 
  BigAEM(27,7) = 0.31943828249997*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(27,9) = 0.4391550328268399*m1r[3]; 
  BigAEM(28,0) = 0.5*m1r[8]; 
  BigAEM(28,1) = 0.4391550328268398*m1r[4]; 
  BigAEM(28,3) = 0.4391550328268399*m1r[6]; 
  BigAEM(28,4) = 0.2981423969999719*m1r[8]+0.4391550328268398*m1r[1]; 
  BigAEM(28,6) = 0.4391550328268399*m1r[3]; 
  BigAEM(28,8) = 0.2981423969999719*m1r[4]+0.5*m1r[0]; 
  BigAEM(29,0) = 0.5*m1r[9]; 
  BigAEM(29,2) = 0.4391550328268398*m1r[5]; 
  BigAEM(29,3) = 0.4391550328268399*m1r[7]; 
  BigAEM(29,5) = 0.2981423969999719*m1r[9]+0.4391550328268398*m1r[2]; 
  BigAEM(29,7) = 0.4391550328268399*m1r[3]; 
  BigAEM(29,9) = 0.2981423969999719*m1r[5]+0.5*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(10,10) = 0.5*m0r[0]; 
  BigAEM(10,11) = 0.5*m0r[1]; 
  BigAEM(10,12) = 0.5*m0r[2]; 
  BigAEM(10,13) = 0.5*m0r[3]; 
  BigAEM(10,14) = 0.5*m0r[4]; 
  BigAEM(10,15) = 0.5*m0r[5]; 
  BigAEM(10,16) = 0.5*m0r[6]; 
  BigAEM(10,17) = 0.5*m0r[7]; 
  BigAEM(10,18) = 0.5*m0r[8]; 
  BigAEM(10,19) = 0.5*m0r[9]; 
  BigAEM(11,10) = 0.5*m0r[1]; 
  BigAEM(11,11) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(11,12) = 0.5*m0r[3]; 
  BigAEM(11,13) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(11,14) = 0.4391550328268398*m0r[8]+0.4472135954999579*m0r[1]; 
  BigAEM(11,15) = 0.5000000000000001*m0r[7]; 
  BigAEM(11,16) = 0.447213595499958*m0r[3]; 
  BigAEM(11,17) = 0.5000000000000001*m0r[5]; 
  BigAEM(11,18) = 0.4391550328268398*m0r[4]; 
  BigAEM(12,10) = 0.5*m0r[2]; 
  BigAEM(12,11) = 0.5*m0r[3]; 
  BigAEM(12,12) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(12,13) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(12,14) = 0.5000000000000001*m0r[6]; 
  BigAEM(12,15) = 0.4391550328268398*m0r[9]+0.4472135954999579*m0r[2]; 
  BigAEM(12,16) = 0.5000000000000001*m0r[4]; 
  BigAEM(12,17) = 0.447213595499958*m0r[3]; 
  BigAEM(12,19) = 0.4391550328268398*m0r[5]; 
  BigAEM(13,10) = 0.5*m0r[3]; 
  BigAEM(13,11) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(13,12) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(13,13) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(13,14) = 0.4472135954999579*m0r[3]; 
  BigAEM(13,15) = 0.4472135954999579*m0r[3]; 
  BigAEM(13,16) = 0.4391550328268399*m0r[8]+0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(13,17) = 0.4391550328268399*m0r[9]+0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(13,18) = 0.4391550328268399*m0r[6]; 
  BigAEM(13,19) = 0.4391550328268399*m0r[7]; 
  BigAEM(14,10) = 0.5*m0r[4]; 
  BigAEM(14,11) = 0.4391550328268398*m0r[8]+0.4472135954999579*m0r[1]; 
  BigAEM(14,12) = 0.5000000000000001*m0r[6]; 
  BigAEM(14,13) = 0.4472135954999579*m0r[3]; 
  BigAEM(14,14) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(14,16) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(14,17) = 0.4472135954999579*m0r[7]; 
  BigAEM(14,18) = 0.2981423969999719*m0r[8]+0.4391550328268398*m0r[1]; 
  BigAEM(15,10) = 0.5*m0r[5]; 
  BigAEM(15,11) = 0.5000000000000001*m0r[7]; 
  BigAEM(15,12) = 0.4391550328268398*m0r[9]+0.4472135954999579*m0r[2]; 
  BigAEM(15,13) = 0.4472135954999579*m0r[3]; 
  BigAEM(15,15) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
  BigAEM(15,16) = 0.4472135954999579*m0r[6]; 
  BigAEM(15,17) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(15,19) = 0.2981423969999719*m0r[9]+0.4391550328268398*m0r[2]; 
  BigAEM(16,10) = 0.5*m0r[6]; 
  BigAEM(16,11) = 0.447213595499958*m0r[3]; 
  BigAEM(16,12) = 0.5000000000000001*m0r[4]; 
  BigAEM(16,13) = 0.4391550328268399*m0r[8]+0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(16,14) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(16,15) = 0.4472135954999579*m0r[6]; 
  BigAEM(16,16) = 0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(16,17) = 0.4*m0r[3]; 
  BigAEM(16,18) = 0.4391550328268399*m0r[3]; 
  BigAEM(17,10) = 0.5*m0r[7]; 
  BigAEM(17,11) = 0.5000000000000001*m0r[5]; 
  BigAEM(17,12) = 0.447213595499958*m0r[3]; 
  BigAEM(17,13) = 0.4391550328268399*m0r[9]+0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(17,14) = 0.4472135954999579*m0r[7]; 
  BigAEM(17,15) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(17,16) = 0.4*m0r[3]; 
  BigAEM(17,17) = 0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(17,19) = 0.4391550328268399*m0r[3]; 
  BigAEM(18,10) = 0.5*m0r[8]; 
  BigAEM(18,11) = 0.4391550328268398*m0r[4]; 
  BigAEM(18,13) = 0.4391550328268399*m0r[6]; 
  BigAEM(18,14) = 0.2981423969999719*m0r[8]+0.4391550328268398*m0r[1]; 
  BigAEM(18,16) = 0.4391550328268399*m0r[3]; 
  BigAEM(18,18) = 0.2981423969999719*m0r[4]+0.5*m0r[0]; 
  BigAEM(19,10) = 0.5*m0r[9]; 
  BigAEM(19,12) = 0.4391550328268398*m0r[5]; 
  BigAEM(19,13) = 0.4391550328268399*m0r[7]; 
  BigAEM(19,15) = 0.2981423969999719*m0r[9]+0.4391550328268398*m0r[2]; 
  BigAEM(19,17) = 0.4391550328268399*m0r[3]; 
  BigAEM(19,19) = 0.2981423969999719*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(10,20) = -0.5*cM[10]; 
  BigAEM(10,21) = -0.5*cM[11]; 
  BigAEM(10,22) = -0.5*cM[12]; 
  BigAEM(10,23) = -0.5*cM[13]; 
  BigAEM(10,24) = -0.5*cM[14]; 
  BigAEM(10,25) = -0.5*cM[15]; 
  BigAEM(10,26) = -0.5*cM[16]; 
  BigAEM(10,27) = -0.5*cM[17]; 
  BigAEM(10,28) = -0.5*cM[18]; 
  BigAEM(10,29) = -0.5*cM[19]; 
  BigAEM(11,20) = -0.5*cM[11]; 
  BigAEM(11,21) = (-0.4472135954999579*cM[14])-0.5*cM[10]; 
  BigAEM(11,22) = -0.5*cM[13]; 
  BigAEM(11,23) = (-0.447213595499958*cM[16])-0.5*cM[12]; 
  BigAEM(11,24) = (-0.4391550328268398*cM[18])-0.4472135954999579*cM[11]; 
  BigAEM(11,25) = -0.5000000000000001*cM[17]; 
  BigAEM(11,26) = -0.447213595499958*cM[13]; 
  BigAEM(11,27) = -0.5000000000000001*cM[15]; 
  BigAEM(11,28) = -0.4391550328268398*cM[14]; 
  BigAEM(12,20) = -0.5*cM[12]; 
  BigAEM(12,21) = -0.5*cM[13]; 
  BigAEM(12,22) = (-0.4472135954999579*cM[15])-0.5*cM[10]; 
  BigAEM(12,23) = (-0.447213595499958*cM[17])-0.5*cM[11]; 
  BigAEM(12,24) = -0.5000000000000001*cM[16]; 
  BigAEM(12,25) = (-0.4391550328268398*cM[19])-0.4472135954999579*cM[12]; 
  BigAEM(12,26) = -0.5000000000000001*cM[14]; 
  BigAEM(12,27) = -0.447213595499958*cM[13]; 
  BigAEM(12,29) = -0.4391550328268398*cM[15]; 
  BigAEM(13,20) = -0.5*cM[13]; 
  BigAEM(13,21) = (-0.447213595499958*cM[16])-0.5*cM[12]; 
  BigAEM(13,22) = (-0.447213595499958*cM[17])-0.5*cM[11]; 
  BigAEM(13,23) = (-0.4472135954999579*cM[15])-0.4472135954999579*cM[14]-0.5*cM[10]; 
  BigAEM(13,24) = -0.4472135954999579*cM[13]; 
  BigAEM(13,25) = -0.4472135954999579*cM[13]; 
  BigAEM(13,26) = (-0.4391550328268399*cM[18])-0.4*cM[17]-0.447213595499958*cM[11]; 
  BigAEM(13,27) = (-0.4391550328268399*cM[19])-0.4*cM[16]-0.447213595499958*cM[12]; 
  BigAEM(13,28) = -0.4391550328268399*cM[16]; 
  BigAEM(13,29) = -0.4391550328268399*cM[17]; 
  BigAEM(14,20) = -0.5*cM[14]; 
  BigAEM(14,21) = (-0.4391550328268398*cM[18])-0.4472135954999579*cM[11]; 
  BigAEM(14,22) = -0.5000000000000001*cM[16]; 
  BigAEM(14,23) = -0.4472135954999579*cM[13]; 
  BigAEM(14,24) = (-0.31943828249997*cM[14])-0.5*cM[10]; 
  BigAEM(14,26) = (-0.31943828249997*cM[16])-0.5000000000000001*cM[12]; 
  BigAEM(14,27) = -0.4472135954999579*cM[17]; 
  BigAEM(14,28) = (-0.2981423969999719*cM[18])-0.4391550328268398*cM[11]; 
  BigAEM(15,20) = -0.5*cM[15]; 
  BigAEM(15,21) = -0.5000000000000001*cM[17]; 
  BigAEM(15,22) = (-0.4391550328268398*cM[19])-0.4472135954999579*cM[12]; 
  BigAEM(15,23) = -0.4472135954999579*cM[13]; 
  BigAEM(15,25) = (-0.31943828249997*cM[15])-0.5*cM[10]; 
  BigAEM(15,26) = -0.4472135954999579*cM[16]; 
  BigAEM(15,27) = (-0.31943828249997*cM[17])-0.5000000000000001*cM[11]; 
  BigAEM(15,29) = (-0.2981423969999719*cM[19])-0.4391550328268398*cM[12]; 
  BigAEM(16,20) = -0.5*cM[16]; 
  BigAEM(16,21) = -0.447213595499958*cM[13]; 
  BigAEM(16,22) = -0.5000000000000001*cM[14]; 
  BigAEM(16,23) = (-0.4391550328268399*cM[18])-0.4*cM[17]-0.447213595499958*cM[11]; 
  BigAEM(16,24) = (-0.31943828249997*cM[16])-0.5000000000000001*cM[12]; 
  BigAEM(16,25) = -0.4472135954999579*cM[16]; 
  BigAEM(16,26) = (-0.4472135954999579*cM[15])-0.31943828249997*cM[14]-0.5*cM[10]; 
  BigAEM(16,27) = -0.4*cM[13]; 
  BigAEM(16,28) = -0.4391550328268399*cM[13]; 
  BigAEM(17,20) = -0.5*cM[17]; 
  BigAEM(17,21) = -0.5000000000000001*cM[15]; 
  BigAEM(17,22) = -0.447213595499958*cM[13]; 
  BigAEM(17,23) = (-0.4391550328268399*cM[19])-0.4*cM[16]-0.447213595499958*cM[12]; 
  BigAEM(17,24) = -0.4472135954999579*cM[17]; 
  BigAEM(17,25) = (-0.31943828249997*cM[17])-0.5000000000000001*cM[11]; 
  BigAEM(17,26) = -0.4*cM[13]; 
  BigAEM(17,27) = (-0.31943828249997*cM[15])-0.4472135954999579*cM[14]-0.5*cM[10]; 
  BigAEM(17,29) = -0.4391550328268399*cM[13]; 
  BigAEM(18,20) = -0.5*cM[18]; 
  BigAEM(18,21) = -0.4391550328268398*cM[14]; 
  BigAEM(18,23) = -0.4391550328268399*cM[16]; 
  BigAEM(18,24) = (-0.2981423969999719*cM[18])-0.4391550328268398*cM[11]; 
  BigAEM(18,26) = -0.4391550328268399*cM[13]; 
  BigAEM(18,28) = (-0.2981423969999719*cM[14])-0.5*cM[10]; 
  BigAEM(19,20) = -0.5*cM[19]; 
  BigAEM(19,22) = -0.4391550328268398*cM[15]; 
  BigAEM(19,23) = -0.4391550328268399*cM[17]; 
  BigAEM(19,25) = (-0.2981423969999719*cM[19])-0.4391550328268398*cM[12]; 
  BigAEM(19,27) = -0.4391550328268399*cM[13]; 
  BigAEM(19,29) = (-0.2981423969999719*cM[15])-0.5*cM[10]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(20,10) = 0.5*m1r[10]; 
  BigAEM(20,11) = 0.5*m1r[11]; 
  BigAEM(20,12) = 0.5*m1r[12]; 
  BigAEM(20,13) = 0.5*m1r[13]; 
  BigAEM(20,14) = 0.5*m1r[14]; 
  BigAEM(20,15) = 0.5*m1r[15]; 
  BigAEM(20,16) = 0.5*m1r[16]; 
  BigAEM(20,17) = 0.5*m1r[17]; 
  BigAEM(20,18) = 0.5*m1r[18]; 
  BigAEM(20,19) = 0.5*m1r[19]; 
  BigAEM(21,10) = 0.5*m1r[11]; 
  BigAEM(21,11) = 0.4472135954999579*m1r[14]+0.5*m1r[10]; 
  BigAEM(21,12) = 0.5*m1r[13]; 
  BigAEM(21,13) = 0.447213595499958*m1r[16]+0.5*m1r[12]; 
  BigAEM(21,14) = 0.4391550328268398*m1r[18]+0.4472135954999579*m1r[11]; 
  BigAEM(21,15) = 0.5000000000000001*m1r[17]; 
  BigAEM(21,16) = 0.447213595499958*m1r[13]; 
  BigAEM(21,17) = 0.5000000000000001*m1r[15]; 
  BigAEM(21,18) = 0.4391550328268398*m1r[14]; 
  BigAEM(22,10) = 0.5*m1r[12]; 
  BigAEM(22,11) = 0.5*m1r[13]; 
  BigAEM(22,12) = 0.4472135954999579*m1r[15]+0.5*m1r[10]; 
  BigAEM(22,13) = 0.447213595499958*m1r[17]+0.5*m1r[11]; 
  BigAEM(22,14) = 0.5000000000000001*m1r[16]; 
  BigAEM(22,15) = 0.4391550328268398*m1r[19]+0.4472135954999579*m1r[12]; 
  BigAEM(22,16) = 0.5000000000000001*m1r[14]; 
  BigAEM(22,17) = 0.447213595499958*m1r[13]; 
  BigAEM(22,19) = 0.4391550328268398*m1r[15]; 
  BigAEM(23,10) = 0.5*m1r[13]; 
  BigAEM(23,11) = 0.447213595499958*m1r[16]+0.5*m1r[12]; 
  BigAEM(23,12) = 0.447213595499958*m1r[17]+0.5*m1r[11]; 
  BigAEM(23,13) = 0.4472135954999579*m1r[15]+0.4472135954999579*m1r[14]+0.5*m1r[10]; 
  BigAEM(23,14) = 0.4472135954999579*m1r[13]; 
  BigAEM(23,15) = 0.4472135954999579*m1r[13]; 
  BigAEM(23,16) = 0.4391550328268399*m1r[18]+0.4*m1r[17]+0.447213595499958*m1r[11]; 
  BigAEM(23,17) = 0.4391550328268399*m1r[19]+0.4*m1r[16]+0.447213595499958*m1r[12]; 
  BigAEM(23,18) = 0.4391550328268399*m1r[16]; 
  BigAEM(23,19) = 0.4391550328268399*m1r[17]; 
  BigAEM(24,10) = 0.5*m1r[14]; 
  BigAEM(24,11) = 0.4391550328268398*m1r[18]+0.4472135954999579*m1r[11]; 
  BigAEM(24,12) = 0.5000000000000001*m1r[16]; 
  BigAEM(24,13) = 0.4472135954999579*m1r[13]; 
  BigAEM(24,14) = 0.31943828249997*m1r[14]+0.5*m1r[10]; 
  BigAEM(24,16) = 0.31943828249997*m1r[16]+0.5000000000000001*m1r[12]; 
  BigAEM(24,17) = 0.4472135954999579*m1r[17]; 
  BigAEM(24,18) = 0.2981423969999719*m1r[18]+0.4391550328268398*m1r[11]; 
  BigAEM(25,10) = 0.5*m1r[15]; 
  BigAEM(25,11) = 0.5000000000000001*m1r[17]; 
  BigAEM(25,12) = 0.4391550328268398*m1r[19]+0.4472135954999579*m1r[12]; 
  BigAEM(25,13) = 0.4472135954999579*m1r[13]; 
  BigAEM(25,15) = 0.31943828249997*m1r[15]+0.5*m1r[10]; 
  BigAEM(25,16) = 0.4472135954999579*m1r[16]; 
  BigAEM(25,17) = 0.31943828249997*m1r[17]+0.5000000000000001*m1r[11]; 
  BigAEM(25,19) = 0.2981423969999719*m1r[19]+0.4391550328268398*m1r[12]; 
  BigAEM(26,10) = 0.5*m1r[16]; 
  BigAEM(26,11) = 0.447213595499958*m1r[13]; 
  BigAEM(26,12) = 0.5000000000000001*m1r[14]; 
  BigAEM(26,13) = 0.4391550328268399*m1r[18]+0.4*m1r[17]+0.447213595499958*m1r[11]; 
  BigAEM(26,14) = 0.31943828249997*m1r[16]+0.5000000000000001*m1r[12]; 
  BigAEM(26,15) = 0.4472135954999579*m1r[16]; 
  BigAEM(26,16) = 0.4472135954999579*m1r[15]+0.31943828249997*m1r[14]+0.5*m1r[10]; 
  BigAEM(26,17) = 0.4*m1r[13]; 
  BigAEM(26,18) = 0.4391550328268399*m1r[13]; 
  BigAEM(27,10) = 0.5*m1r[17]; 
  BigAEM(27,11) = 0.5000000000000001*m1r[15]; 
  BigAEM(27,12) = 0.447213595499958*m1r[13]; 
  BigAEM(27,13) = 0.4391550328268399*m1r[19]+0.4*m1r[16]+0.447213595499958*m1r[12]; 
  BigAEM(27,14) = 0.4472135954999579*m1r[17]; 
  BigAEM(27,15) = 0.31943828249997*m1r[17]+0.5000000000000001*m1r[11]; 
  BigAEM(27,16) = 0.4*m1r[13]; 
  BigAEM(27,17) = 0.31943828249997*m1r[15]+0.4472135954999579*m1r[14]+0.5*m1r[10]; 
  BigAEM(27,19) = 0.4391550328268399*m1r[13]; 
  BigAEM(28,10) = 0.5*m1r[18]; 
  BigAEM(28,11) = 0.4391550328268398*m1r[14]; 
  BigAEM(28,13) = 0.4391550328268399*m1r[16]; 
  BigAEM(28,14) = 0.2981423969999719*m1r[18]+0.4391550328268398*m1r[11]; 
  BigAEM(28,16) = 0.4391550328268399*m1r[13]; 
  BigAEM(28,18) = 0.2981423969999719*m1r[14]+0.5*m1r[10]; 
  BigAEM(29,10) = 0.5*m1r[19]; 
  BigAEM(29,12) = 0.4391550328268398*m1r[15]; 
  BigAEM(29,13) = 0.4391550328268399*m1r[17]; 
  BigAEM(29,15) = 0.2981423969999719*m1r[19]+0.4391550328268398*m1r[12]; 
  BigAEM(29,17) = 0.4391550328268399*m1r[13]; 
  BigAEM(29,19) = 0.2981423969999719*m1r[15]+0.5*m1r[10]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(20,20) = m0r[0]-0.5*cE[0]; 
  BigAEM(20,21) = m0r[1]-0.5*cE[1]; 
  BigAEM(20,22) = m0r[2]-0.5*cE[2]; 
  BigAEM(20,23) = m0r[3]-0.5*cE[3]; 
  BigAEM(20,24) = m0r[4]-0.5*cE[4]; 
  BigAEM(20,25) = m0r[5]-0.5*cE[5]; 
  BigAEM(20,26) = m0r[6]-0.5*cE[6]; 
  BigAEM(20,27) = m0r[7]-0.5*cE[7]; 
  BigAEM(20,28) = m0r[8]-0.5*cE[8]; 
  BigAEM(20,29) = m0r[9]-0.5*cE[9]; 
  BigAEM(21,20) = m0r[1]-0.5*cE[1]; 
  BigAEM(21,21) = 0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(21,22) = m0r[3]-0.5*cE[3]; 
  BigAEM(21,23) = 0.8944271909999161*m0r[6]-0.447213595499958*cE[6]+m0r[2]-0.5*cE[2]; 
  BigAEM(21,24) = 0.8783100656536796*m0r[8]-0.4391550328268398*cE[8]+0.8944271909999159*m0r[1]-0.4472135954999579*cE[1]; 
  BigAEM(21,25) = 1.0*m0r[7]-0.5000000000000001*cE[7]; 
  BigAEM(21,26) = 0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(21,27) = 1.0*m0r[5]-0.5000000000000001*cE[5]; 
  BigAEM(21,28) = 0.8783100656536796*m0r[4]-0.4391550328268398*cE[4]; 
  BigAEM(22,20) = m0r[2]-0.5*cE[2]; 
  BigAEM(22,21) = m0r[3]-0.5*cE[3]; 
  BigAEM(22,22) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+m0r[0]-0.5*cE[0]; 
  BigAEM(22,23) = 0.8944271909999161*m0r[7]-0.447213595499958*cE[7]+m0r[1]-0.5*cE[1]; 
  BigAEM(22,24) = 1.0*m0r[6]-0.5000000000000001*cE[6]; 
  BigAEM(22,25) = 0.8783100656536796*m0r[9]-0.4391550328268398*cE[9]+0.8944271909999159*m0r[2]-0.4472135954999579*cE[2]; 
  BigAEM(22,26) = 1.0*m0r[4]-0.5000000000000001*cE[4]; 
  BigAEM(22,27) = 0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(22,29) = 0.8783100656536796*m0r[5]-0.4391550328268398*cE[5]; 
  BigAEM(23,20) = m0r[3]-0.5*cE[3]; 
  BigAEM(23,21) = 0.8944271909999161*m0r[6]-0.447213595499958*cE[6]+m0r[2]-0.5*cE[2]; 
  BigAEM(23,22) = 0.8944271909999161*m0r[7]-0.447213595499958*cE[7]+m0r[1]-0.5*cE[1]; 
  BigAEM(23,23) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(23,24) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(23,25) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(23,26) = 0.8783100656536798*m0r[8]-0.4391550328268399*cE[8]+0.8*m0r[7]-0.4*cE[7]+0.8944271909999161*m0r[1]-0.447213595499958*cE[1]; 
  BigAEM(23,27) = 0.8783100656536798*m0r[9]-0.4391550328268399*cE[9]+0.8*m0r[6]-0.4*cE[6]+0.8944271909999161*m0r[2]-0.447213595499958*cE[2]; 
  BigAEM(23,28) = 0.8783100656536798*m0r[6]-0.4391550328268399*cE[6]; 
  BigAEM(23,29) = 0.8783100656536798*m0r[7]-0.4391550328268399*cE[7]; 
  BigAEM(24,20) = m0r[4]-0.5*cE[4]; 
  BigAEM(24,21) = 0.8783100656536796*m0r[8]-0.4391550328268398*cE[8]+0.8944271909999159*m0r[1]-0.4472135954999579*cE[1]; 
  BigAEM(24,22) = 1.0*m0r[6]-0.5000000000000001*cE[6]; 
  BigAEM(24,23) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(24,24) = 0.6388765649999399*m0r[4]-0.31943828249997*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(24,26) = 0.6388765649999399*m0r[6]-0.31943828249997*cE[6]+1.0*m0r[2]-0.5000000000000001*cE[2]; 
  BigAEM(24,27) = 0.8944271909999159*m0r[7]-0.4472135954999579*cE[7]; 
  BigAEM(24,28) = 0.5962847939999438*m0r[8]-0.2981423969999719*cE[8]+0.8783100656536796*m0r[1]-0.4391550328268398*cE[1]; 
  BigAEM(25,20) = m0r[5]-0.5*cE[5]; 
  BigAEM(25,21) = 1.0*m0r[7]-0.5000000000000001*cE[7]; 
  BigAEM(25,22) = 0.8783100656536796*m0r[9]-0.4391550328268398*cE[9]+0.8944271909999159*m0r[2]-0.4472135954999579*cE[2]; 
  BigAEM(25,23) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(25,25) = 0.6388765649999399*m0r[5]-0.31943828249997*cE[5]+m0r[0]-0.5*cE[0]; 
  BigAEM(25,26) = 0.8944271909999159*m0r[6]-0.4472135954999579*cE[6]; 
  BigAEM(25,27) = 0.6388765649999399*m0r[7]-0.31943828249997*cE[7]+1.0*m0r[1]-0.5000000000000001*cE[1]; 
  BigAEM(25,29) = 0.5962847939999438*m0r[9]-0.2981423969999719*cE[9]+0.8783100656536796*m0r[2]-0.4391550328268398*cE[2]; 
  BigAEM(26,20) = m0r[6]-0.5*cE[6]; 
  BigAEM(26,21) = 0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(26,22) = 1.0*m0r[4]-0.5000000000000001*cE[4]; 
  BigAEM(26,23) = 0.8783100656536798*m0r[8]-0.4391550328268399*cE[8]+0.8*m0r[7]-0.4*cE[7]+0.8944271909999161*m0r[1]-0.447213595499958*cE[1]; 
  BigAEM(26,24) = 0.6388765649999399*m0r[6]-0.31943828249997*cE[6]+1.0*m0r[2]-0.5000000000000001*cE[2]; 
  BigAEM(26,25) = 0.8944271909999159*m0r[6]-0.4472135954999579*cE[6]; 
  BigAEM(26,26) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+0.6388765649999399*m0r[4]-0.31943828249997*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(26,27) = 0.8*m0r[3]-0.4*cE[3]; 
  BigAEM(26,28) = 0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(27,20) = m0r[7]-0.5*cE[7]; 
  BigAEM(27,21) = 1.0*m0r[5]-0.5000000000000001*cE[5]; 
  BigAEM(27,22) = 0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(27,23) = 0.8783100656536798*m0r[9]-0.4391550328268399*cE[9]+0.8*m0r[6]-0.4*cE[6]+0.8944271909999161*m0r[2]-0.447213595499958*cE[2]; 
  BigAEM(27,24) = 0.8944271909999159*m0r[7]-0.4472135954999579*cE[7]; 
  BigAEM(27,25) = 0.6388765649999399*m0r[7]-0.31943828249997*cE[7]+1.0*m0r[1]-0.5000000000000001*cE[1]; 
  BigAEM(27,26) = 0.8*m0r[3]-0.4*cE[3]; 
  BigAEM(27,27) = 0.6388765649999399*m0r[5]-0.31943828249997*cE[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(27,29) = 0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(28,20) = m0r[8]-0.5*cE[8]; 
  BigAEM(28,21) = 0.8783100656536796*m0r[4]-0.4391550328268398*cE[4]; 
  BigAEM(28,23) = 0.8783100656536798*m0r[6]-0.4391550328268399*cE[6]; 
  BigAEM(28,24) = 0.5962847939999438*m0r[8]-0.2981423969999719*cE[8]+0.8783100656536796*m0r[1]-0.4391550328268398*cE[1]; 
  BigAEM(28,26) = 0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(28,28) = 0.5962847939999438*m0r[4]-0.2981423969999719*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(29,20) = m0r[9]-0.5*cE[9]; 
  BigAEM(29,22) = 0.8783100656536796*m0r[5]-0.4391550328268398*cE[5]; 
  BigAEM(29,23) = 0.8783100656536798*m0r[7]-0.4391550328268399*cE[7]; 
  BigAEM(29,25) = 0.5962847939999438*m0r[9]-0.2981423969999719*cE[9]+0.8783100656536796*m0r[2]-0.4391550328268398*cE[2]; 
  BigAEM(29,27) = 0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(29,29) = 0.5962847939999438*m0r[5]-0.2981423969999719*cE[5]+m0r[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<10,10>(0,10).setZero(); 
  BigAEM.block<10,10>(10,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m1r[12],m1r[13],m1r[14],m1r[15],m1r[16],m1r[17],m1r[18],m1r[19],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5],m2r[6],m2r[7],m2r[8],m2r[9]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,20,1) = xEV.segment<20>(0); 
 
  Eigen::Map<VectorXd>(vtSq,10,1) = xEV.segment<10>(20); 
 
} 
 
