#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
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
void VmBoundaryIntegral2x2vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[16]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]-1.0*fIn[5])*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[6]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[7]+fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]+fIn[5])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[48]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[13])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += ((-2.23606797749979*fIn[23])+1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
    out[2] += ((-2.23606797749979*fIn[24])+1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
    out[3] += ((-2.23606797749979*fIn[34])+1.732050807568877*fIn[15]-1.0*fIn[5])*dS; 
    out[4] += (1.732050807568877*fIn[21]-1.0*fIn[11])*dS; 
    out[5] += (1.732050807568877*fIn[22]-1.0*fIn[12])*dS; 
    out[6] += (1.732050807568877*fIn[32]-1.0*fIn[19])*dS; 
    out[7] += (1.732050807568877*fIn[33]-1.0*fIn[20])*dS; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[1] += (2.23606797749979*fIn[23]+1.732050807568877*fIn[6]+fIn[1])*dS; 
    out[2] += (2.23606797749979*fIn[24]+1.732050807568877*fIn[7]+fIn[2])*dS; 
    out[3] += (2.23606797749979*fIn[34]+1.732050807568877*fIn[15]+fIn[5])*dS; 
    out[4] += (1.732050807568877*fIn[21]+fIn[11])*dS; 
    out[5] += (1.732050807568877*fIn[22]+fIn[12])*dS; 
    out[6] += (1.732050807568877*fIn[32]+fIn[19])*dS; 
    out[7] += (1.732050807568877*fIn[33]+fIn[20])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[80]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[33]-2.23606797749979*fIn[13]+1.732050807568877*fIn[3]-1.0*fIn[0])*dS; 
    out[1] += (2.645751311064591*fIn[52]-2.23606797749979*fIn[23]+1.732050807568877*fIn[6]-1.0*fIn[1])*dS; 
    out[2] += (2.645751311064591*fIn[53]-2.23606797749979*fIn[24]+1.732050807568877*fIn[7]-1.0*fIn[2])*dS; 
    out[3] += (2.645751311064591*fIn[66]-2.23606797749979*fIn[38]+1.732050807568877*fIn[15]-1.0*fIn[5])*dS; 
    out[4] += (1.732050807568877*fIn[21]-1.0*fIn[11])*dS; 
    out[5] += (1.732050807568877*fIn[22]-1.0*fIn[12])*dS; 
    out[6] += (1.732050807568877*fIn[36]-1.0*fIn[19])*dS; 
    out[7] += (1.732050807568877*fIn[37]-1.0*fIn[20])*dS; 
    out[8] += (1.732050807568877*fIn[50]-1.0*fIn[31])*dS; 
    out[9] += (1.732050807568877*fIn[51]-1.0*fIn[32])*dS; 
    out[10] += (1.732050807568877*fIn[64]-1.0*fIn[48])*dS; 
    out[11] += (1.732050807568877*fIn[65]-1.0*fIn[49])*dS; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[33]+2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS; 
    out[1] += (2.645751311064591*fIn[52]+2.23606797749979*fIn[23]+1.732050807568877*fIn[6]+fIn[1])*dS; 
    out[2] += (2.645751311064591*fIn[53]+2.23606797749979*fIn[24]+1.732050807568877*fIn[7]+fIn[2])*dS; 
    out[3] += (2.645751311064591*fIn[66]+2.23606797749979*fIn[38]+1.732050807568877*fIn[15]+fIn[5])*dS; 
    out[4] += (1.732050807568877*fIn[21]+fIn[11])*dS; 
    out[5] += (1.732050807568877*fIn[22]+fIn[12])*dS; 
    out[6] += (1.732050807568877*fIn[36]+fIn[19])*dS; 
    out[7] += (1.732050807568877*fIn[37]+fIn[20])*dS; 
    out[8] += (1.732050807568877*fIn[50]+fIn[31])*dS; 
    out[9] += (1.732050807568877*fIn[51]+fIn[32])*dS; 
    out[10] += (1.732050807568877*fIn[64]+fIn[48])*dS; 
    out[11] += (1.732050807568877*fIn[65]+fIn[49])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[16]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary+(0.8660254037844386*dxv[2]*fIn[3]-0.5*fIn[0]*dxv[2])*dS; 
    out[1] += (1.732050807568877*fIn[6]-1.0*fIn[1])*dS*vBoundary+(0.8660254037844386*dxv[2]*fIn[6]-0.5*fIn[1]*dxv[2])*dS; 
    out[2] += (1.732050807568877*fIn[7]-1.0*fIn[2])*dS*vBoundary+(0.8660254037844386*dxv[2]*fIn[7]-0.5*dxv[2]*fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]-1.0*fIn[5])*dS*vBoundary+(0.8660254037844386*dxv[2]*fIn[11]-0.5*dxv[2]*fIn[5])*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary+((-0.8660254037844386*dxv[2]*fIn[3])-0.5*fIn[0]*dxv[2])*dS; 
    out[1] += (1.732050807568877*fIn[6]+fIn[1])*dS*vBoundary+((-0.8660254037844386*dxv[2]*fIn[6])-0.5*fIn[1]*dxv[2])*dS; 
    out[2] += (1.732050807568877*fIn[7]+fIn[2])*dS*vBoundary+((-0.8660254037844386*dxv[2]*fIn[7])-0.5*dxv[2]*fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]+fIn[5])*dS*vBoundary+((-0.8660254037844386*dxv[2]*fIn[11])-0.5*dxv[2]*fIn[5])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[48]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[13])+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[23])+1.732050807568877*fIn[6]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += ((-2.23606797749979*fIn[24])+1.732050807568877*fIn[7]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += ((-2.23606797749979*fIn[34])+1.732050807568877*fIn[15]-1.0*fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[21]-1.0*fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[22]-1.0*fIn[12])*dS*vBoundary; 
    out[6] += (1.732050807568877*fIn[32]-1.0*fIn[19])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[33]-1.0*fIn[20])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[23]+1.732050807568877*fIn[6]+fIn[1])*dS*vBoundary; 
    out[2] += (2.23606797749979*fIn[24]+1.732050807568877*fIn[7]+fIn[2])*dS*vBoundary; 
    out[3] += (2.23606797749979*fIn[34]+1.732050807568877*fIn[15]+fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[21]+fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[22]+fIn[12])*dS*vBoundary; 
    out[6] += (1.732050807568877*fIn[32]+fIn[19])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[33]+fIn[20])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[80]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[33]-2.23606797749979*fIn[13]+1.732050807568877*fIn[3]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[52]-2.23606797749979*fIn[23]+1.732050807568877*fIn[6]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (2.645751311064591*fIn[53]-2.23606797749979*fIn[24]+1.732050807568877*fIn[7]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += (2.645751311064591*fIn[66]-2.23606797749979*fIn[38]+1.732050807568877*fIn[15]-1.0*fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[21]-1.0*fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[22]-1.0*fIn[12])*dS*vBoundary; 
    out[6] += (1.732050807568877*fIn[36]-1.0*fIn[19])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[37]-1.0*fIn[20])*dS*vBoundary; 
    out[8] += (1.732050807568877*fIn[50]-1.0*fIn[31])*dS*vBoundary; 
    out[9] += (1.732050807568877*fIn[51]-1.0*fIn[32])*dS*vBoundary; 
    out[10] += (1.732050807568877*fIn[64]-1.0*fIn[48])*dS*vBoundary; 
    out[11] += (1.732050807568877*fIn[65]-1.0*fIn[49])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[33]+2.23606797749979*fIn[13]+1.732050807568877*fIn[3]+fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[52]+2.23606797749979*fIn[23]+1.732050807568877*fIn[6]+fIn[1])*dS*vBoundary; 
    out[2] += (2.645751311064591*fIn[53]+2.23606797749979*fIn[24]+1.732050807568877*fIn[7]+fIn[2])*dS*vBoundary; 
    out[3] += (2.645751311064591*fIn[66]+2.23606797749979*fIn[38]+1.732050807568877*fIn[15]+fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[21]+fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[22]+fIn[12])*dS*vBoundary; 
    out[6] += (1.732050807568877*fIn[36]+fIn[19])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[37]+fIn[20])*dS*vBoundary; 
    out[8] += (1.732050807568877*fIn[50]+fIn[31])*dS*vBoundary; 
    out[9] += (1.732050807568877*fIn[51]+fIn[32])*dS*vBoundary; 
    out[10] += (1.732050807568877*fIn[64]+fIn[48])*dS*vBoundary; 
    out[11] += (1.732050807568877*fIn[65]+fIn[49])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[16]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
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
 
void VmBoundaryIntegral2x2vSer_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[48]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[8] += ((-2.23606797749979*fIn[14])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[9] += ((-2.23606797749979*fIn[28])+1.732050807568877*fIn[8]-1.0*fIn[1])*dS; 
    out[10] += ((-2.23606797749979*fIn[29])+1.732050807568877*fIn[9]-1.0*fIn[2])*dS; 
    out[11] += ((-2.23606797749979*fIn[41])+1.732050807568877*fIn[16]-1.0*fIn[5])*dS; 
    out[12] += (1.732050807568877*fIn[25]-1.0*fIn[11])*dS; 
    out[13] += (1.732050807568877*fIn[26]-1.0*fIn[12])*dS; 
    out[14] += (1.732050807568877*fIn[35]-1.0*fIn[19])*dS; 
    out[15] += (1.732050807568877*fIn[36]-1.0*fIn[20])*dS; 
 
  } else {
 
    out[8] += (2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[9] += (2.23606797749979*fIn[28]+1.732050807568877*fIn[8]+fIn[1])*dS; 
    out[10] += (2.23606797749979*fIn[29]+1.732050807568877*fIn[9]+fIn[2])*dS; 
    out[11] += (2.23606797749979*fIn[41]+1.732050807568877*fIn[16]+fIn[5])*dS; 
    out[12] += (1.732050807568877*fIn[25]+fIn[11])*dS; 
    out[13] += (1.732050807568877*fIn[26]+fIn[12])*dS; 
    out[14] += (1.732050807568877*fIn[35]+fIn[19])*dS; 
    out[15] += (1.732050807568877*fIn[36]+fIn[20])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[80]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[12] += (2.645751311064591*fIn[34]-2.23606797749979*fIn[14]+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[13] += (2.645751311064591*fIn[57]-2.23606797749979*fIn[28]+1.732050807568877*fIn[8]-1.0*fIn[1])*dS; 
    out[14] += (2.645751311064591*fIn[58]-2.23606797749979*fIn[29]+1.732050807568877*fIn[9]-1.0*fIn[2])*dS; 
    out[15] += (2.645751311064591*fIn[73]-2.23606797749979*fIn[45]+1.732050807568877*fIn[16]-1.0*fIn[5])*dS; 
    out[16] += (1.732050807568877*fIn[25]-1.0*fIn[11])*dS; 
    out[17] += (1.732050807568877*fIn[26]-1.0*fIn[12])*dS; 
    out[18] += (1.732050807568877*fIn[39]-1.0*fIn[19])*dS; 
    out[19] += (1.732050807568877*fIn[40]-1.0*fIn[20])*dS; 
    out[20] += (1.732050807568877*fIn[54]-1.0*fIn[31])*dS; 
    out[21] += (1.732050807568877*fIn[55]-1.0*fIn[32])*dS; 
    out[22] += (1.732050807568877*fIn[67]-1.0*fIn[48])*dS; 
    out[23] += (1.732050807568877*fIn[68]-1.0*fIn[49])*dS; 
 
  } else {
 
    out[12] += (2.645751311064591*fIn[34]+2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[13] += (2.645751311064591*fIn[57]+2.23606797749979*fIn[28]+1.732050807568877*fIn[8]+fIn[1])*dS; 
    out[14] += (2.645751311064591*fIn[58]+2.23606797749979*fIn[29]+1.732050807568877*fIn[9]+fIn[2])*dS; 
    out[15] += (2.645751311064591*fIn[73]+2.23606797749979*fIn[45]+1.732050807568877*fIn[16]+fIn[5])*dS; 
    out[16] += (1.732050807568877*fIn[25]+fIn[11])*dS; 
    out[17] += (1.732050807568877*fIn[26]+fIn[12])*dS; 
    out[18] += (1.732050807568877*fIn[39]+fIn[19])*dS; 
    out[19] += (1.732050807568877*fIn[40]+fIn[20])*dS; 
    out[20] += (1.732050807568877*fIn[54]+fIn[31])*dS; 
    out[21] += (1.732050807568877*fIn[55]+fIn[32])*dS; 
    out[22] += (1.732050807568877*fIn[67]+fIn[48])*dS; 
    out[23] += (1.732050807568877*fIn[68]+fIn[49])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[16]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[4]-0.5*fIn[0]*dxv[3])*dS; 
    out[1] += (1.732050807568877*fIn[8]-1.0*fIn[1])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[8]-0.5*fIn[1]*dxv[3])*dS; 
    out[2] += (1.732050807568877*fIn[9]-1.0*fIn[2])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[9]-0.5*fIn[2]*dxv[3])*dS; 
    out[3] += (1.732050807568877*fIn[12]-1.0*fIn[5])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[12]-0.5*dxv[3]*fIn[5])*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[4])-0.5*fIn[0]*dxv[3])*dS; 
    out[1] += (1.732050807568877*fIn[8]+fIn[1])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[8])-0.5*fIn[1]*dxv[3])*dS; 
    out[2] += (1.732050807568877*fIn[9]+fIn[2])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[9])-0.5*fIn[2]*dxv[3])*dS; 
    out[3] += (1.732050807568877*fIn[12]+fIn[5])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[12])-0.5*dxv[3]*fIn[5])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[48]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[14])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[28])+1.732050807568877*fIn[8]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += ((-2.23606797749979*fIn[29])+1.732050807568877*fIn[9]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += ((-2.23606797749979*fIn[41])+1.732050807568877*fIn[16]-1.0*fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[25]-1.0*fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[26]-1.0*fIn[12])*dS*vBoundary; 
    out[6] += (1.732050807568877*fIn[35]-1.0*fIn[19])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[36]-1.0*fIn[20])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[28]+1.732050807568877*fIn[8]+fIn[1])*dS*vBoundary; 
    out[2] += (2.23606797749979*fIn[29]+1.732050807568877*fIn[9]+fIn[2])*dS*vBoundary; 
    out[3] += (2.23606797749979*fIn[41]+1.732050807568877*fIn[16]+fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[25]+fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[26]+fIn[12])*dS*vBoundary; 
    out[6] += (1.732050807568877*fIn[35]+fIn[19])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[36]+fIn[20])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral2x2vSer_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[4]:    cell length in each direciton. 
  // fIn[80]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  if (atLower) {
 
    out[0] += (2.645751311064591*fIn[34]-2.23606797749979*fIn[14]+1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[57]-2.23606797749979*fIn[28]+1.732050807568877*fIn[8]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (2.645751311064591*fIn[58]-2.23606797749979*fIn[29]+1.732050807568877*fIn[9]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += (2.645751311064591*fIn[73]-2.23606797749979*fIn[45]+1.732050807568877*fIn[16]-1.0*fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[25]-1.0*fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[26]-1.0*fIn[12])*dS*vBoundary; 
    out[6] += (1.732050807568877*fIn[39]-1.0*fIn[19])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[40]-1.0*fIn[20])*dS*vBoundary; 
    out[8] += (1.732050807568877*fIn[54]-1.0*fIn[31])*dS*vBoundary; 
    out[9] += (1.732050807568877*fIn[55]-1.0*fIn[32])*dS*vBoundary; 
    out[10] += (1.732050807568877*fIn[67]-1.0*fIn[48])*dS*vBoundary; 
    out[11] += (1.732050807568877*fIn[68]-1.0*fIn[49])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.645751311064591*fIn[34]+2.23606797749979*fIn[14]+1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary; 
    out[1] += (2.645751311064591*fIn[57]+2.23606797749979*fIn[28]+1.732050807568877*fIn[8]+fIn[1])*dS*vBoundary; 
    out[2] += (2.645751311064591*fIn[58]+2.23606797749979*fIn[29]+1.732050807568877*fIn[9]+fIn[2])*dS*vBoundary; 
    out[3] += (2.645751311064591*fIn[73]+2.23606797749979*fIn[45]+1.732050807568877*fIn[16]+fIn[5])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[25]+fIn[11])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[26]+fIn[12])*dS*vBoundary; 
    out[6] += (1.732050807568877*fIn[39]+fIn[19])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[40]+fIn[20])*dS*vBoundary; 
    out[8] += (1.732050807568877*fIn[54]+fIn[31])*dS*vBoundary; 
    out[9] += (1.732050807568877*fIn[55]+fIn[32])*dS*vBoundary; 
    out[10] += (1.732050807568877*fIn[67]+fIn[48])*dS*vBoundary; 
    out[11] += (1.732050807568877*fIn[68]+fIn[49])*dS*vBoundary; 
 
  }
 
} 
 
void VmSelfPrimMoments2x2vSer_P1(const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
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
  double m1r[8]; 
  double m0Sr[4]; 
  double m1Sr[8]; 
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
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
    m2Sr[2] = m2S[2]; 
    m2Sr[3] = m2S[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(12,12); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(12);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(12);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0r[0]; 
  BigAEM(0,1) = 0.5*m0r[1]; 
  BigAEM(0,2) = 0.5*m0r[2]; 
  BigAEM(0,3) = 0.5*m0r[3]; 
  BigAEM(1,0) = 0.5*m0r[1]; 
  BigAEM(1,1) = 0.5*m0r[0]; 
  BigAEM(1,2) = 0.5*m0r[3]; 
  BigAEM(1,3) = 0.5*m0r[2]; 
  BigAEM(2,0) = 0.5*m0r[2]; 
  BigAEM(2,1) = 0.5*m0r[3]; 
  BigAEM(2,2) = 0.5*m0r[0]; 
  BigAEM(2,3) = 0.5*m0r[1]; 
  BigAEM(3,0) = 0.5*m0r[3]; 
  BigAEM(3,1) = 0.5*m0r[2]; 
  BigAEM(3,2) = 0.5*m0r[1]; 
  BigAEM(3,3) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,8) = -0.5*cM[0]; 
  BigAEM(0,9) = -0.5*cM[1]; 
  BigAEM(0,10) = -0.5*cM[2]; 
  BigAEM(0,11) = -0.5*cM[3]; 
  BigAEM(1,8) = -0.5*cM[1]; 
  BigAEM(1,9) = -0.5*cM[0]; 
  BigAEM(1,10) = -0.5*cM[3]; 
  BigAEM(1,11) = -0.5*cM[2]; 
  BigAEM(2,8) = -0.5*cM[2]; 
  BigAEM(2,9) = -0.5*cM[3]; 
  BigAEM(2,10) = -0.5*cM[0]; 
  BigAEM(2,11) = -0.5*cM[1]; 
  BigAEM(3,8) = -0.5*cM[3]; 
  BigAEM(3,9) = -0.5*cM[2]; 
  BigAEM(3,10) = -0.5*cM[1]; 
  BigAEM(3,11) = -0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(8,0) = 0.5*m1Sr[0]; 
  BigAEM(8,1) = 0.5*m1Sr[1]; 
  BigAEM(8,2) = 0.5*m1Sr[2]; 
  BigAEM(8,3) = 0.5*m1Sr[3]; 
  BigAEM(9,0) = 0.5*m1Sr[1]; 
  BigAEM(9,1) = 0.5*m1Sr[0]; 
  BigAEM(9,2) = 0.5*m1Sr[3]; 
  BigAEM(9,3) = 0.5*m1Sr[2]; 
  BigAEM(10,0) = 0.5*m1Sr[2]; 
  BigAEM(10,1) = 0.5*m1Sr[3]; 
  BigAEM(10,2) = 0.5*m1Sr[0]; 
  BigAEM(10,3) = 0.5*m1Sr[1]; 
  BigAEM(11,0) = 0.5*m1Sr[3]; 
  BigAEM(11,1) = 0.5*m1Sr[2]; 
  BigAEM(11,2) = 0.5*m1Sr[1]; 
  BigAEM(11,3) = 0.5*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(4,4) = 0.5*m0r[0]; 
  BigAEM(4,5) = 0.5*m0r[1]; 
  BigAEM(4,6) = 0.5*m0r[2]; 
  BigAEM(4,7) = 0.5*m0r[3]; 
  BigAEM(5,4) = 0.5*m0r[1]; 
  BigAEM(5,5) = 0.5*m0r[0]; 
  BigAEM(5,6) = 0.5*m0r[3]; 
  BigAEM(5,7) = 0.5*m0r[2]; 
  BigAEM(6,4) = 0.5*m0r[2]; 
  BigAEM(6,5) = 0.5*m0r[3]; 
  BigAEM(6,6) = 0.5*m0r[0]; 
  BigAEM(6,7) = 0.5*m0r[1]; 
  BigAEM(7,4) = 0.5*m0r[3]; 
  BigAEM(7,5) = 0.5*m0r[2]; 
  BigAEM(7,6) = 0.5*m0r[1]; 
  BigAEM(7,7) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(4,8) = -0.5*cM[4]; 
  BigAEM(4,9) = -0.5*cM[5]; 
  BigAEM(4,10) = -0.5*cM[6]; 
  BigAEM(4,11) = -0.5*cM[7]; 
  BigAEM(5,8) = -0.5*cM[5]; 
  BigAEM(5,9) = -0.5*cM[4]; 
  BigAEM(5,10) = -0.5*cM[7]; 
  BigAEM(5,11) = -0.5*cM[6]; 
  BigAEM(6,8) = -0.5*cM[6]; 
  BigAEM(6,9) = -0.5*cM[7]; 
  BigAEM(6,10) = -0.5*cM[4]; 
  BigAEM(6,11) = -0.5*cM[5]; 
  BigAEM(7,8) = -0.5*cM[7]; 
  BigAEM(7,9) = -0.5*cM[6]; 
  BigAEM(7,10) = -0.5*cM[5]; 
  BigAEM(7,11) = -0.5*cM[4]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(8,4) = 0.5*m1Sr[4]; 
  BigAEM(8,5) = 0.5*m1Sr[5]; 
  BigAEM(8,6) = 0.5*m1Sr[6]; 
  BigAEM(8,7) = 0.5*m1Sr[7]; 
  BigAEM(9,4) = 0.5*m1Sr[5]; 
  BigAEM(9,5) = 0.5*m1Sr[4]; 
  BigAEM(9,6) = 0.5*m1Sr[7]; 
  BigAEM(9,7) = 0.5*m1Sr[6]; 
  BigAEM(10,4) = 0.5*m1Sr[6]; 
  BigAEM(10,5) = 0.5*m1Sr[7]; 
  BigAEM(10,6) = 0.5*m1Sr[4]; 
  BigAEM(10,7) = 0.5*m1Sr[5]; 
  BigAEM(11,4) = 0.5*m1Sr[7]; 
  BigAEM(11,5) = 0.5*m1Sr[6]; 
  BigAEM(11,6) = 0.5*m1Sr[5]; 
  BigAEM(11,7) = 0.5*m1Sr[4]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(8,8) = 0.5*m0Sr[0]-0.5*cE[0]; 
  BigAEM(8,9) = 0.5*m0Sr[1]-0.5*cE[1]; 
  BigAEM(8,10) = 0.5*m0Sr[2]-0.5*cE[2]; 
  BigAEM(8,11) = 0.5*m0Sr[3]-0.5*cE[3]; 
  BigAEM(9,8) = 0.5*m0Sr[1]-0.5*cE[1]; 
  BigAEM(9,9) = 0.5*m0Sr[0]-0.5*cE[0]; 
  BigAEM(9,10) = 0.5*m0Sr[3]-0.5*cE[3]; 
  BigAEM(9,11) = 0.5*m0Sr[2]-0.5*cE[2]; 
  BigAEM(10,8) = 0.5*m0Sr[2]-0.5*cE[2]; 
  BigAEM(10,9) = 0.5*m0Sr[3]-0.5*cE[3]; 
  BigAEM(10,10) = 0.5*m0Sr[0]-0.5*cE[0]; 
  BigAEM(10,11) = 0.5*m0Sr[1]-0.5*cE[1]; 
  BigAEM(11,8) = 0.5*m0Sr[3]-0.5*cE[3]; 
  BigAEM(11,9) = 0.5*m0Sr[2]-0.5*cE[2]; 
  BigAEM(11,10) = 0.5*m0Sr[1]-0.5*cE[1]; 
  BigAEM(11,11) = 0.5*m0Sr[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<4,4>(0,4).setZero(); 
  BigAEM.block<4,4>(4,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m2Sr[0],m2Sr[1],m2Sr[2],m2Sr[3]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,8,1) = xEV.segment<8>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = xEV.segment<4>(8); 
 
} 
 
void VmSelfPrimMoments2x2vSer_P2(const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
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
  double m1r[16]; 
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
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(24,24); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(24);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(24);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0r[0]; 
  BigAEM(0,1) = 0.5*m0r[1]; 
  BigAEM(0,2) = 0.5*m0r[2]; 
  BigAEM(0,3) = 0.5*m0r[3]; 
  BigAEM(0,4) = 0.5*m0r[4]; 
  BigAEM(0,5) = 0.5*m0r[5]; 
  BigAEM(0,6) = 0.5*m0r[6]; 
  BigAEM(0,7) = 0.5*m0r[7]; 
  BigAEM(1,0) = 0.5*m0r[1]; 
  BigAEM(1,1) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(1,2) = 0.5*m0r[3]; 
  BigAEM(1,3) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(1,4) = 0.4472135954999579*m0r[1]; 
  BigAEM(1,5) = 0.5000000000000001*m0r[7]; 
  BigAEM(1,6) = 0.447213595499958*m0r[3]; 
  BigAEM(1,7) = 0.5000000000000001*m0r[5]; 
  BigAEM(2,0) = 0.5*m0r[2]; 
  BigAEM(2,1) = 0.5*m0r[3]; 
  BigAEM(2,2) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(2,3) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(2,4) = 0.5000000000000001*m0r[6]; 
  BigAEM(2,5) = 0.4472135954999579*m0r[2]; 
  BigAEM(2,6) = 0.5000000000000001*m0r[4]; 
  BigAEM(2,7) = 0.447213595499958*m0r[3]; 
  BigAEM(3,0) = 0.5*m0r[3]; 
  BigAEM(3,1) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(3,2) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(3,3) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(3,4) = 0.4472135954999579*m0r[3]; 
  BigAEM(3,5) = 0.4472135954999579*m0r[3]; 
  BigAEM(3,6) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(3,7) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(4,0) = 0.5*m0r[4]; 
  BigAEM(4,1) = 0.4472135954999579*m0r[1]; 
  BigAEM(4,2) = 0.5000000000000001*m0r[6]; 
  BigAEM(4,3) = 0.4472135954999579*m0r[3]; 
  BigAEM(4,4) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(4,6) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(4,7) = 0.4472135954999579*m0r[7]; 
  BigAEM(5,0) = 0.5*m0r[5]; 
  BigAEM(5,1) = 0.5000000000000001*m0r[7]; 
  BigAEM(5,2) = 0.4472135954999579*m0r[2]; 
  BigAEM(5,3) = 0.4472135954999579*m0r[3]; 
  BigAEM(5,5) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
  BigAEM(5,6) = 0.4472135954999579*m0r[6]; 
  BigAEM(5,7) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(6,0) = 0.5*m0r[6]; 
  BigAEM(6,1) = 0.447213595499958*m0r[3]; 
  BigAEM(6,2) = 0.5000000000000001*m0r[4]; 
  BigAEM(6,3) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(6,4) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(6,5) = 0.4472135954999579*m0r[6]; 
  BigAEM(6,6) = 0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(6,7) = 0.4*m0r[3]; 
  BigAEM(7,0) = 0.5*m0r[7]; 
  BigAEM(7,1) = 0.5000000000000001*m0r[5]; 
  BigAEM(7,2) = 0.447213595499958*m0r[3]; 
  BigAEM(7,3) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(7,4) = 0.4472135954999579*m0r[7]; 
  BigAEM(7,5) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(7,6) = 0.4*m0r[3]; 
  BigAEM(7,7) = 0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,16) = -0.5*cM[0]; 
  BigAEM(0,17) = -0.5*cM[1]; 
  BigAEM(0,18) = -0.5*cM[2]; 
  BigAEM(0,19) = -0.5*cM[3]; 
  BigAEM(0,20) = -0.5*cM[4]; 
  BigAEM(0,21) = -0.5*cM[5]; 
  BigAEM(0,22) = -0.5*cM[6]; 
  BigAEM(0,23) = -0.5*cM[7]; 
  BigAEM(1,16) = -0.5*cM[1]; 
  BigAEM(1,17) = (-0.4472135954999579*cM[4])-0.5*cM[0]; 
  BigAEM(1,18) = -0.5*cM[3]; 
  BigAEM(1,19) = (-0.447213595499958*cM[6])-0.5*cM[2]; 
  BigAEM(1,20) = -0.4472135954999579*cM[1]; 
  BigAEM(1,21) = -0.5000000000000001*cM[7]; 
  BigAEM(1,22) = -0.447213595499958*cM[3]; 
  BigAEM(1,23) = -0.5000000000000001*cM[5]; 
  BigAEM(2,16) = -0.5*cM[2]; 
  BigAEM(2,17) = -0.5*cM[3]; 
  BigAEM(2,18) = (-0.4472135954999579*cM[5])-0.5*cM[0]; 
  BigAEM(2,19) = (-0.447213595499958*cM[7])-0.5*cM[1]; 
  BigAEM(2,20) = -0.5000000000000001*cM[6]; 
  BigAEM(2,21) = -0.4472135954999579*cM[2]; 
  BigAEM(2,22) = -0.5000000000000001*cM[4]; 
  BigAEM(2,23) = -0.447213595499958*cM[3]; 
  BigAEM(3,16) = -0.5*cM[3]; 
  BigAEM(3,17) = (-0.447213595499958*cM[6])-0.5*cM[2]; 
  BigAEM(3,18) = (-0.447213595499958*cM[7])-0.5*cM[1]; 
  BigAEM(3,19) = (-0.4472135954999579*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(3,20) = -0.4472135954999579*cM[3]; 
  BigAEM(3,21) = -0.4472135954999579*cM[3]; 
  BigAEM(3,22) = (-0.4*cM[7])-0.447213595499958*cM[1]; 
  BigAEM(3,23) = (-0.4*cM[6])-0.447213595499958*cM[2]; 
  BigAEM(4,16) = -0.5*cM[4]; 
  BigAEM(4,17) = -0.4472135954999579*cM[1]; 
  BigAEM(4,18) = -0.5000000000000001*cM[6]; 
  BigAEM(4,19) = -0.4472135954999579*cM[3]; 
  BigAEM(4,20) = (-0.31943828249997*cM[4])-0.5*cM[0]; 
  BigAEM(4,22) = (-0.31943828249997*cM[6])-0.5000000000000001*cM[2]; 
  BigAEM(4,23) = -0.4472135954999579*cM[7]; 
  BigAEM(5,16) = -0.5*cM[5]; 
  BigAEM(5,17) = -0.5000000000000001*cM[7]; 
  BigAEM(5,18) = -0.4472135954999579*cM[2]; 
  BigAEM(5,19) = -0.4472135954999579*cM[3]; 
  BigAEM(5,21) = (-0.31943828249997*cM[5])-0.5*cM[0]; 
  BigAEM(5,22) = -0.4472135954999579*cM[6]; 
  BigAEM(5,23) = (-0.31943828249997*cM[7])-0.5000000000000001*cM[1]; 
  BigAEM(6,16) = -0.5*cM[6]; 
  BigAEM(6,17) = -0.447213595499958*cM[3]; 
  BigAEM(6,18) = -0.5000000000000001*cM[4]; 
  BigAEM(6,19) = (-0.4*cM[7])-0.447213595499958*cM[1]; 
  BigAEM(6,20) = (-0.31943828249997*cM[6])-0.5000000000000001*cM[2]; 
  BigAEM(6,21) = -0.4472135954999579*cM[6]; 
  BigAEM(6,22) = (-0.4472135954999579*cM[5])-0.31943828249997*cM[4]-0.5*cM[0]; 
  BigAEM(6,23) = -0.4*cM[3]; 
  BigAEM(7,16) = -0.5*cM[7]; 
  BigAEM(7,17) = -0.5000000000000001*cM[5]; 
  BigAEM(7,18) = -0.447213595499958*cM[3]; 
  BigAEM(7,19) = (-0.4*cM[6])-0.447213595499958*cM[2]; 
  BigAEM(7,20) = -0.4472135954999579*cM[7]; 
  BigAEM(7,21) = (-0.31943828249997*cM[7])-0.5000000000000001*cM[1]; 
  BigAEM(7,22) = -0.4*cM[3]; 
  BigAEM(7,23) = (-0.31943828249997*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(16,0) = 0.5*m1r[0]; 
  BigAEM(16,1) = 0.5*m1r[1]; 
  BigAEM(16,2) = 0.5*m1r[2]; 
  BigAEM(16,3) = 0.5*m1r[3]; 
  BigAEM(16,4) = 0.5*m1r[4]; 
  BigAEM(16,5) = 0.5*m1r[5]; 
  BigAEM(16,6) = 0.5*m1r[6]; 
  BigAEM(16,7) = 0.5*m1r[7]; 
  BigAEM(17,0) = 0.5*m1r[1]; 
  BigAEM(17,1) = 0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(17,2) = 0.5*m1r[3]; 
  BigAEM(17,3) = 0.447213595499958*m1r[6]+0.5*m1r[2]; 
  BigAEM(17,4) = 0.4472135954999579*m1r[1]; 
  BigAEM(17,5) = 0.5000000000000001*m1r[7]; 
  BigAEM(17,6) = 0.447213595499958*m1r[3]; 
  BigAEM(17,7) = 0.5000000000000001*m1r[5]; 
  BigAEM(18,0) = 0.5*m1r[2]; 
  BigAEM(18,1) = 0.5*m1r[3]; 
  BigAEM(18,2) = 0.4472135954999579*m1r[5]+0.5*m1r[0]; 
  BigAEM(18,3) = 0.447213595499958*m1r[7]+0.5*m1r[1]; 
  BigAEM(18,4) = 0.5000000000000001*m1r[6]; 
  BigAEM(18,5) = 0.4472135954999579*m1r[2]; 
  BigAEM(18,6) = 0.5000000000000001*m1r[4]; 
  BigAEM(18,7) = 0.447213595499958*m1r[3]; 
  BigAEM(19,0) = 0.5*m1r[3]; 
  BigAEM(19,1) = 0.447213595499958*m1r[6]+0.5*m1r[2]; 
  BigAEM(19,2) = 0.447213595499958*m1r[7]+0.5*m1r[1]; 
  BigAEM(19,3) = 0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(19,4) = 0.4472135954999579*m1r[3]; 
  BigAEM(19,5) = 0.4472135954999579*m1r[3]; 
  BigAEM(19,6) = 0.4*m1r[7]+0.447213595499958*m1r[1]; 
  BigAEM(19,7) = 0.4*m1r[6]+0.447213595499958*m1r[2]; 
  BigAEM(20,0) = 0.5*m1r[4]; 
  BigAEM(20,1) = 0.4472135954999579*m1r[1]; 
  BigAEM(20,2) = 0.5000000000000001*m1r[6]; 
  BigAEM(20,3) = 0.4472135954999579*m1r[3]; 
  BigAEM(20,4) = 0.31943828249997*m1r[4]+0.5*m1r[0]; 
  BigAEM(20,6) = 0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]; 
  BigAEM(20,7) = 0.4472135954999579*m1r[7]; 
  BigAEM(21,0) = 0.5*m1r[5]; 
  BigAEM(21,1) = 0.5000000000000001*m1r[7]; 
  BigAEM(21,2) = 0.4472135954999579*m1r[2]; 
  BigAEM(21,3) = 0.4472135954999579*m1r[3]; 
  BigAEM(21,5) = 0.31943828249997*m1r[5]+0.5*m1r[0]; 
  BigAEM(21,6) = 0.4472135954999579*m1r[6]; 
  BigAEM(21,7) = 0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]; 
  BigAEM(22,0) = 0.5*m1r[6]; 
  BigAEM(22,1) = 0.447213595499958*m1r[3]; 
  BigAEM(22,2) = 0.5000000000000001*m1r[4]; 
  BigAEM(22,3) = 0.4*m1r[7]+0.447213595499958*m1r[1]; 
  BigAEM(22,4) = 0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]; 
  BigAEM(22,5) = 0.4472135954999579*m1r[6]; 
  BigAEM(22,6) = 0.4472135954999579*m1r[5]+0.31943828249997*m1r[4]+0.5*m1r[0]; 
  BigAEM(22,7) = 0.4*m1r[3]; 
  BigAEM(23,0) = 0.5*m1r[7]; 
  BigAEM(23,1) = 0.5000000000000001*m1r[5]; 
  BigAEM(23,2) = 0.447213595499958*m1r[3]; 
  BigAEM(23,3) = 0.4*m1r[6]+0.447213595499958*m1r[2]; 
  BigAEM(23,4) = 0.4472135954999579*m1r[7]; 
  BigAEM(23,5) = 0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]; 
  BigAEM(23,6) = 0.4*m1r[3]; 
  BigAEM(23,7) = 0.31943828249997*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(8,8) = 0.5*m0r[0]; 
  BigAEM(8,9) = 0.5*m0r[1]; 
  BigAEM(8,10) = 0.5*m0r[2]; 
  BigAEM(8,11) = 0.5*m0r[3]; 
  BigAEM(8,12) = 0.5*m0r[4]; 
  BigAEM(8,13) = 0.5*m0r[5]; 
  BigAEM(8,14) = 0.5*m0r[6]; 
  BigAEM(8,15) = 0.5*m0r[7]; 
  BigAEM(9,8) = 0.5*m0r[1]; 
  BigAEM(9,9) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(9,10) = 0.5*m0r[3]; 
  BigAEM(9,11) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(9,12) = 0.4472135954999579*m0r[1]; 
  BigAEM(9,13) = 0.5000000000000001*m0r[7]; 
  BigAEM(9,14) = 0.447213595499958*m0r[3]; 
  BigAEM(9,15) = 0.5000000000000001*m0r[5]; 
  BigAEM(10,8) = 0.5*m0r[2]; 
  BigAEM(10,9) = 0.5*m0r[3]; 
  BigAEM(10,10) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(10,11) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(10,12) = 0.5000000000000001*m0r[6]; 
  BigAEM(10,13) = 0.4472135954999579*m0r[2]; 
  BigAEM(10,14) = 0.5000000000000001*m0r[4]; 
  BigAEM(10,15) = 0.447213595499958*m0r[3]; 
  BigAEM(11,8) = 0.5*m0r[3]; 
  BigAEM(11,9) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(11,10) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(11,11) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(11,12) = 0.4472135954999579*m0r[3]; 
  BigAEM(11,13) = 0.4472135954999579*m0r[3]; 
  BigAEM(11,14) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(11,15) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(12,8) = 0.5*m0r[4]; 
  BigAEM(12,9) = 0.4472135954999579*m0r[1]; 
  BigAEM(12,10) = 0.5000000000000001*m0r[6]; 
  BigAEM(12,11) = 0.4472135954999579*m0r[3]; 
  BigAEM(12,12) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(12,14) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(12,15) = 0.4472135954999579*m0r[7]; 
  BigAEM(13,8) = 0.5*m0r[5]; 
  BigAEM(13,9) = 0.5000000000000001*m0r[7]; 
  BigAEM(13,10) = 0.4472135954999579*m0r[2]; 
  BigAEM(13,11) = 0.4472135954999579*m0r[3]; 
  BigAEM(13,13) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
  BigAEM(13,14) = 0.4472135954999579*m0r[6]; 
  BigAEM(13,15) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(14,8) = 0.5*m0r[6]; 
  BigAEM(14,9) = 0.447213595499958*m0r[3]; 
  BigAEM(14,10) = 0.5000000000000001*m0r[4]; 
  BigAEM(14,11) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(14,12) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(14,13) = 0.4472135954999579*m0r[6]; 
  BigAEM(14,14) = 0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(14,15) = 0.4*m0r[3]; 
  BigAEM(15,8) = 0.5*m0r[7]; 
  BigAEM(15,9) = 0.5000000000000001*m0r[5]; 
  BigAEM(15,10) = 0.447213595499958*m0r[3]; 
  BigAEM(15,11) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(15,12) = 0.4472135954999579*m0r[7]; 
  BigAEM(15,13) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(15,14) = 0.4*m0r[3]; 
  BigAEM(15,15) = 0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(8,16) = -0.5*cM[8]; 
  BigAEM(8,17) = -0.5*cM[9]; 
  BigAEM(8,18) = -0.5*cM[10]; 
  BigAEM(8,19) = -0.5*cM[11]; 
  BigAEM(8,20) = -0.5*cM[12]; 
  BigAEM(8,21) = -0.5*cM[13]; 
  BigAEM(8,22) = -0.5*cM[14]; 
  BigAEM(8,23) = -0.5*cM[15]; 
  BigAEM(9,16) = -0.5*cM[9]; 
  BigAEM(9,17) = (-0.4472135954999579*cM[12])-0.5*cM[8]; 
  BigAEM(9,18) = -0.5*cM[11]; 
  BigAEM(9,19) = (-0.447213595499958*cM[14])-0.5*cM[10]; 
  BigAEM(9,20) = -0.4472135954999579*cM[9]; 
  BigAEM(9,21) = -0.5000000000000001*cM[15]; 
  BigAEM(9,22) = -0.447213595499958*cM[11]; 
  BigAEM(9,23) = -0.5000000000000001*cM[13]; 
  BigAEM(10,16) = -0.5*cM[10]; 
  BigAEM(10,17) = -0.5*cM[11]; 
  BigAEM(10,18) = (-0.4472135954999579*cM[13])-0.5*cM[8]; 
  BigAEM(10,19) = (-0.447213595499958*cM[15])-0.5*cM[9]; 
  BigAEM(10,20) = -0.5000000000000001*cM[14]; 
  BigAEM(10,21) = -0.4472135954999579*cM[10]; 
  BigAEM(10,22) = -0.5000000000000001*cM[12]; 
  BigAEM(10,23) = -0.447213595499958*cM[11]; 
  BigAEM(11,16) = -0.5*cM[11]; 
  BigAEM(11,17) = (-0.447213595499958*cM[14])-0.5*cM[10]; 
  BigAEM(11,18) = (-0.447213595499958*cM[15])-0.5*cM[9]; 
  BigAEM(11,19) = (-0.4472135954999579*cM[13])-0.4472135954999579*cM[12]-0.5*cM[8]; 
  BigAEM(11,20) = -0.4472135954999579*cM[11]; 
  BigAEM(11,21) = -0.4472135954999579*cM[11]; 
  BigAEM(11,22) = (-0.4*cM[15])-0.447213595499958*cM[9]; 
  BigAEM(11,23) = (-0.4*cM[14])-0.447213595499958*cM[10]; 
  BigAEM(12,16) = -0.5*cM[12]; 
  BigAEM(12,17) = -0.4472135954999579*cM[9]; 
  BigAEM(12,18) = -0.5000000000000001*cM[14]; 
  BigAEM(12,19) = -0.4472135954999579*cM[11]; 
  BigAEM(12,20) = (-0.31943828249997*cM[12])-0.5*cM[8]; 
  BigAEM(12,22) = (-0.31943828249997*cM[14])-0.5000000000000001*cM[10]; 
  BigAEM(12,23) = -0.4472135954999579*cM[15]; 
  BigAEM(13,16) = -0.5*cM[13]; 
  BigAEM(13,17) = -0.5000000000000001*cM[15]; 
  BigAEM(13,18) = -0.4472135954999579*cM[10]; 
  BigAEM(13,19) = -0.4472135954999579*cM[11]; 
  BigAEM(13,21) = (-0.31943828249997*cM[13])-0.5*cM[8]; 
  BigAEM(13,22) = -0.4472135954999579*cM[14]; 
  BigAEM(13,23) = (-0.31943828249997*cM[15])-0.5000000000000001*cM[9]; 
  BigAEM(14,16) = -0.5*cM[14]; 
  BigAEM(14,17) = -0.447213595499958*cM[11]; 
  BigAEM(14,18) = -0.5000000000000001*cM[12]; 
  BigAEM(14,19) = (-0.4*cM[15])-0.447213595499958*cM[9]; 
  BigAEM(14,20) = (-0.31943828249997*cM[14])-0.5000000000000001*cM[10]; 
  BigAEM(14,21) = -0.4472135954999579*cM[14]; 
  BigAEM(14,22) = (-0.4472135954999579*cM[13])-0.31943828249997*cM[12]-0.5*cM[8]; 
  BigAEM(14,23) = -0.4*cM[11]; 
  BigAEM(15,16) = -0.5*cM[15]; 
  BigAEM(15,17) = -0.5000000000000001*cM[13]; 
  BigAEM(15,18) = -0.447213595499958*cM[11]; 
  BigAEM(15,19) = (-0.4*cM[14])-0.447213595499958*cM[10]; 
  BigAEM(15,20) = -0.4472135954999579*cM[15]; 
  BigAEM(15,21) = (-0.31943828249997*cM[15])-0.5000000000000001*cM[9]; 
  BigAEM(15,22) = -0.4*cM[11]; 
  BigAEM(15,23) = (-0.31943828249997*cM[13])-0.4472135954999579*cM[12]-0.5*cM[8]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(16,8) = 0.5*m1r[8]; 
  BigAEM(16,9) = 0.5*m1r[9]; 
  BigAEM(16,10) = 0.5*m1r[10]; 
  BigAEM(16,11) = 0.5*m1r[11]; 
  BigAEM(16,12) = 0.5*m1r[12]; 
  BigAEM(16,13) = 0.5*m1r[13]; 
  BigAEM(16,14) = 0.5*m1r[14]; 
  BigAEM(16,15) = 0.5*m1r[15]; 
  BigAEM(17,8) = 0.5*m1r[9]; 
  BigAEM(17,9) = 0.4472135954999579*m1r[12]+0.5*m1r[8]; 
  BigAEM(17,10) = 0.5*m1r[11]; 
  BigAEM(17,11) = 0.447213595499958*m1r[14]+0.5*m1r[10]; 
  BigAEM(17,12) = 0.4472135954999579*m1r[9]; 
  BigAEM(17,13) = 0.5000000000000001*m1r[15]; 
  BigAEM(17,14) = 0.447213595499958*m1r[11]; 
  BigAEM(17,15) = 0.5000000000000001*m1r[13]; 
  BigAEM(18,8) = 0.5*m1r[10]; 
  BigAEM(18,9) = 0.5*m1r[11]; 
  BigAEM(18,10) = 0.4472135954999579*m1r[13]+0.5*m1r[8]; 
  BigAEM(18,11) = 0.447213595499958*m1r[15]+0.5*m1r[9]; 
  BigAEM(18,12) = 0.5000000000000001*m1r[14]; 
  BigAEM(18,13) = 0.4472135954999579*m1r[10]; 
  BigAEM(18,14) = 0.5000000000000001*m1r[12]; 
  BigAEM(18,15) = 0.447213595499958*m1r[11]; 
  BigAEM(19,8) = 0.5*m1r[11]; 
  BigAEM(19,9) = 0.447213595499958*m1r[14]+0.5*m1r[10]; 
  BigAEM(19,10) = 0.447213595499958*m1r[15]+0.5*m1r[9]; 
  BigAEM(19,11) = 0.4472135954999579*m1r[13]+0.4472135954999579*m1r[12]+0.5*m1r[8]; 
  BigAEM(19,12) = 0.4472135954999579*m1r[11]; 
  BigAEM(19,13) = 0.4472135954999579*m1r[11]; 
  BigAEM(19,14) = 0.4*m1r[15]+0.447213595499958*m1r[9]; 
  BigAEM(19,15) = 0.4*m1r[14]+0.447213595499958*m1r[10]; 
  BigAEM(20,8) = 0.5*m1r[12]; 
  BigAEM(20,9) = 0.4472135954999579*m1r[9]; 
  BigAEM(20,10) = 0.5000000000000001*m1r[14]; 
  BigAEM(20,11) = 0.4472135954999579*m1r[11]; 
  BigAEM(20,12) = 0.31943828249997*m1r[12]+0.5*m1r[8]; 
  BigAEM(20,14) = 0.31943828249997*m1r[14]+0.5000000000000001*m1r[10]; 
  BigAEM(20,15) = 0.4472135954999579*m1r[15]; 
  BigAEM(21,8) = 0.5*m1r[13]; 
  BigAEM(21,9) = 0.5000000000000001*m1r[15]; 
  BigAEM(21,10) = 0.4472135954999579*m1r[10]; 
  BigAEM(21,11) = 0.4472135954999579*m1r[11]; 
  BigAEM(21,13) = 0.31943828249997*m1r[13]+0.5*m1r[8]; 
  BigAEM(21,14) = 0.4472135954999579*m1r[14]; 
  BigAEM(21,15) = 0.31943828249997*m1r[15]+0.5000000000000001*m1r[9]; 
  BigAEM(22,8) = 0.5*m1r[14]; 
  BigAEM(22,9) = 0.447213595499958*m1r[11]; 
  BigAEM(22,10) = 0.5000000000000001*m1r[12]; 
  BigAEM(22,11) = 0.4*m1r[15]+0.447213595499958*m1r[9]; 
  BigAEM(22,12) = 0.31943828249997*m1r[14]+0.5000000000000001*m1r[10]; 
  BigAEM(22,13) = 0.4472135954999579*m1r[14]; 
  BigAEM(22,14) = 0.4472135954999579*m1r[13]+0.31943828249997*m1r[12]+0.5*m1r[8]; 
  BigAEM(22,15) = 0.4*m1r[11]; 
  BigAEM(23,8) = 0.5*m1r[15]; 
  BigAEM(23,9) = 0.5000000000000001*m1r[13]; 
  BigAEM(23,10) = 0.447213595499958*m1r[11]; 
  BigAEM(23,11) = 0.4*m1r[14]+0.447213595499958*m1r[10]; 
  BigAEM(23,12) = 0.4472135954999579*m1r[15]; 
  BigAEM(23,13) = 0.31943828249997*m1r[15]+0.5000000000000001*m1r[9]; 
  BigAEM(23,14) = 0.4*m1r[11]; 
  BigAEM(23,15) = 0.31943828249997*m1r[13]+0.4472135954999579*m1r[12]+0.5*m1r[8]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(16,16) = m0r[0]-0.5*cE[0]; 
  BigAEM(16,17) = m0r[1]-0.5*cE[1]; 
  BigAEM(16,18) = m0r[2]-0.5*cE[2]; 
  BigAEM(16,19) = m0r[3]-0.5*cE[3]; 
  BigAEM(16,20) = m0r[4]-0.5*cE[4]; 
  BigAEM(16,21) = m0r[5]-0.5*cE[5]; 
  BigAEM(16,22) = m0r[6]-0.5*cE[6]; 
  BigAEM(16,23) = m0r[7]-0.5*cE[7]; 
  BigAEM(17,16) = m0r[1]-0.5*cE[1]; 
  BigAEM(17,17) = 0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(17,18) = m0r[3]-0.5*cE[3]; 
  BigAEM(17,19) = 0.8944271909999161*m0r[6]-0.447213595499958*cE[6]+m0r[2]-0.5*cE[2]; 
  BigAEM(17,20) = 0.8944271909999159*m0r[1]-0.4472135954999579*cE[1]; 
  BigAEM(17,21) = 1.0*m0r[7]-0.5000000000000001*cE[7]; 
  BigAEM(17,22) = 0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(17,23) = 1.0*m0r[5]-0.5000000000000001*cE[5]; 
  BigAEM(18,16) = m0r[2]-0.5*cE[2]; 
  BigAEM(18,17) = m0r[3]-0.5*cE[3]; 
  BigAEM(18,18) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+m0r[0]-0.5*cE[0]; 
  BigAEM(18,19) = 0.8944271909999161*m0r[7]-0.447213595499958*cE[7]+m0r[1]-0.5*cE[1]; 
  BigAEM(18,20) = 1.0*m0r[6]-0.5000000000000001*cE[6]; 
  BigAEM(18,21) = 0.8944271909999159*m0r[2]-0.4472135954999579*cE[2]; 
  BigAEM(18,22) = 1.0*m0r[4]-0.5000000000000001*cE[4]; 
  BigAEM(18,23) = 0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(19,16) = m0r[3]-0.5*cE[3]; 
  BigAEM(19,17) = 0.8944271909999161*m0r[6]-0.447213595499958*cE[6]+m0r[2]-0.5*cE[2]; 
  BigAEM(19,18) = 0.8944271909999161*m0r[7]-0.447213595499958*cE[7]+m0r[1]-0.5*cE[1]; 
  BigAEM(19,19) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(19,20) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(19,21) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(19,22) = 0.8*m0r[7]-0.4*cE[7]+0.8944271909999161*m0r[1]-0.447213595499958*cE[1]; 
  BigAEM(19,23) = 0.8*m0r[6]-0.4*cE[6]+0.8944271909999161*m0r[2]-0.447213595499958*cE[2]; 
  BigAEM(20,16) = m0r[4]-0.5*cE[4]; 
  BigAEM(20,17) = 0.8944271909999159*m0r[1]-0.4472135954999579*cE[1]; 
  BigAEM(20,18) = 1.0*m0r[6]-0.5000000000000001*cE[6]; 
  BigAEM(20,19) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(20,20) = 0.6388765649999399*m0r[4]-0.31943828249997*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(20,22) = 0.6388765649999399*m0r[6]-0.31943828249997*cE[6]+1.0*m0r[2]-0.5000000000000001*cE[2]; 
  BigAEM(20,23) = 0.8944271909999159*m0r[7]-0.4472135954999579*cE[7]; 
  BigAEM(21,16) = m0r[5]-0.5*cE[5]; 
  BigAEM(21,17) = 1.0*m0r[7]-0.5000000000000001*cE[7]; 
  BigAEM(21,18) = 0.8944271909999159*m0r[2]-0.4472135954999579*cE[2]; 
  BigAEM(21,19) = 0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(21,21) = 0.6388765649999399*m0r[5]-0.31943828249997*cE[5]+m0r[0]-0.5*cE[0]; 
  BigAEM(21,22) = 0.8944271909999159*m0r[6]-0.4472135954999579*cE[6]; 
  BigAEM(21,23) = 0.6388765649999399*m0r[7]-0.31943828249997*cE[7]+1.0*m0r[1]-0.5000000000000001*cE[1]; 
  BigAEM(22,16) = m0r[6]-0.5*cE[6]; 
  BigAEM(22,17) = 0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(22,18) = 1.0*m0r[4]-0.5000000000000001*cE[4]; 
  BigAEM(22,19) = 0.8*m0r[7]-0.4*cE[7]+0.8944271909999161*m0r[1]-0.447213595499958*cE[1]; 
  BigAEM(22,20) = 0.6388765649999399*m0r[6]-0.31943828249997*cE[6]+1.0*m0r[2]-0.5000000000000001*cE[2]; 
  BigAEM(22,21) = 0.8944271909999159*m0r[6]-0.4472135954999579*cE[6]; 
  BigAEM(22,22) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+0.6388765649999399*m0r[4]-0.31943828249997*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(22,23) = 0.8*m0r[3]-0.4*cE[3]; 
  BigAEM(23,16) = m0r[7]-0.5*cE[7]; 
  BigAEM(23,17) = 1.0*m0r[5]-0.5000000000000001*cE[5]; 
  BigAEM(23,18) = 0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(23,19) = 0.8*m0r[6]-0.4*cE[6]+0.8944271909999161*m0r[2]-0.447213595499958*cE[2]; 
  BigAEM(23,20) = 0.8944271909999159*m0r[7]-0.4472135954999579*cE[7]; 
  BigAEM(23,21) = 0.6388765649999399*m0r[7]-0.31943828249997*cE[7]+1.0*m0r[1]-0.5000000000000001*cE[1]; 
  BigAEM(23,22) = 0.8*m0r[3]-0.4*cE[3]; 
  BigAEM(23,23) = 0.6388765649999399*m0r[5]-0.31943828249997*cE[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<8,8>(0,8).setZero(); 
  BigAEM.block<8,8>(8,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m1r[12],m1r[13],m1r[14],m1r[15],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5],m2r[6],m2r[7]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,16,1) = xEV.segment<16>(0); 
 
  Eigen::Map<VectorXd>(vtSq,8,1) = xEV.segment<8>(16); 
 
} 
 
void VmSelfPrimMoments2x2vSer_P3(const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (2.29128784747792*m0[11]+2.29128784747792*m0[10]-1.322875655532295*m0[9]-1.322875655532295*m0[8]-1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if (2.29128784747792*m0[11]+2.29128784747792*m0[10]-1.322875655532295*m0[9]-1.322875655532295*m0[8]-1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-2.29128784747792*m0[11])-2.29128784747792*m0[10]-1.322875655532295*m0[9]+1.322875655532295*m0[8]+1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-2.29128784747792*m0[11])-2.29128784747792*m0[10]-1.322875655532295*m0[9]+1.322875655532295*m0[8]+1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[12]; 
  double m1r[24]; 
  double m2r[12]; 
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
    m0r[10] = 0.0; 
    m0r[11] = 0.0; 
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
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    m1r[12] = m1[12]; 
    m1r[13] = 0.0; 
    m1r[14] = 0.0; 
    m1r[15] = 0.0; 
    m1r[16] = 0.0; 
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
    m2r[8] = 0.0; 
    m2r[9] = 0.0; 
    m2r[10] = 0.0; 
    m2r[11] = 0.0; 
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
    m0r[10] = m0[10]; 
    m0r[11] = m0[11]; 
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
    m2r[8] = m2[8]; 
    m2r[9] = m2[9]; 
    m2r[10] = m2[10]; 
    m2r[11] = m2[11]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(36,36); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(36);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(36);  
 
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
  BigAEM(0,10) = 0.5*m0r[10]; 
  BigAEM(0,11) = 0.5*m0r[11]; 
  BigAEM(1,0) = 0.5*m0r[1]; 
  BigAEM(1,1) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(1,2) = 0.5*m0r[3]; 
  BigAEM(1,3) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(1,4) = 0.4391550328268398*m0r[8]+0.4472135954999579*m0r[1]; 
  BigAEM(1,5) = 0.5000000000000001*m0r[7]; 
  BigAEM(1,6) = 0.4391550328268399*m0r[10]+0.447213595499958*m0r[3]; 
  BigAEM(1,7) = 0.5000000000000001*m0r[5]; 
  BigAEM(1,8) = 0.4391550328268398*m0r[4]; 
  BigAEM(1,9) = 0.5*m0r[11]; 
  BigAEM(1,10) = 0.4391550328268399*m0r[6]; 
  BigAEM(1,11) = 0.5*m0r[9]; 
  BigAEM(2,0) = 0.5*m0r[2]; 
  BigAEM(2,1) = 0.5*m0r[3]; 
  BigAEM(2,2) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(2,3) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(2,4) = 0.5000000000000001*m0r[6]; 
  BigAEM(2,5) = 0.4391550328268398*m0r[9]+0.4472135954999579*m0r[2]; 
  BigAEM(2,6) = 0.5000000000000001*m0r[4]; 
  BigAEM(2,7) = 0.4391550328268399*m0r[11]+0.447213595499958*m0r[3]; 
  BigAEM(2,8) = 0.5*m0r[10]; 
  BigAEM(2,9) = 0.4391550328268398*m0r[5]; 
  BigAEM(2,10) = 0.5*m0r[8]; 
  BigAEM(2,11) = 0.4391550328268399*m0r[7]; 
  BigAEM(3,0) = 0.5*m0r[3]; 
  BigAEM(3,1) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(3,2) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(3,3) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(3,4) = 0.4391550328268399*m0r[10]+0.4472135954999579*m0r[3]; 
  BigAEM(3,5) = 0.4391550328268399*m0r[11]+0.4472135954999579*m0r[3]; 
  BigAEM(3,6) = 0.4391550328268399*m0r[8]+0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(3,7) = 0.4391550328268399*m0r[9]+0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(3,8) = 0.4391550328268399*m0r[6]; 
  BigAEM(3,9) = 0.4391550328268399*m0r[7]; 
  BigAEM(3,10) = 0.4391550328268399*m0r[4]; 
  BigAEM(3,11) = 0.4391550328268399*m0r[5]; 
  BigAEM(4,0) = 0.5*m0r[4]; 
  BigAEM(4,1) = 0.4391550328268398*m0r[8]+0.4472135954999579*m0r[1]; 
  BigAEM(4,2) = 0.5000000000000001*m0r[6]; 
  BigAEM(4,3) = 0.4391550328268399*m0r[10]+0.4472135954999579*m0r[3]; 
  BigAEM(4,4) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(4,6) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(4,7) = 0.4472135954999579*m0r[7]; 
  BigAEM(4,8) = 0.2981423969999719*m0r[8]+0.4391550328268398*m0r[1]; 
  BigAEM(4,10) = 0.2981423969999719*m0r[10]+0.4391550328268399*m0r[3]; 
  BigAEM(4,11) = 0.4472135954999579*m0r[11]; 
  BigAEM(5,0) = 0.5*m0r[5]; 
  BigAEM(5,1) = 0.5000000000000001*m0r[7]; 
  BigAEM(5,2) = 0.4391550328268398*m0r[9]+0.4472135954999579*m0r[2]; 
  BigAEM(5,3) = 0.4391550328268399*m0r[11]+0.4472135954999579*m0r[3]; 
  BigAEM(5,5) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
  BigAEM(5,6) = 0.4472135954999579*m0r[6]; 
  BigAEM(5,7) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(5,9) = 0.2981423969999719*m0r[9]+0.4391550328268398*m0r[2]; 
  BigAEM(5,10) = 0.4472135954999579*m0r[10]; 
  BigAEM(5,11) = 0.2981423969999719*m0r[11]+0.4391550328268399*m0r[3]; 
  BigAEM(6,0) = 0.5*m0r[6]; 
  BigAEM(6,1) = 0.4391550328268399*m0r[10]+0.447213595499958*m0r[3]; 
  BigAEM(6,2) = 0.5000000000000001*m0r[4]; 
  BigAEM(6,3) = 0.4391550328268399*m0r[8]+0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(6,4) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(6,5) = 0.4472135954999579*m0r[6]; 
  BigAEM(6,6) = 0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(6,7) = 0.3927922024247863*m0r[11]+0.3927922024247863*m0r[10]+0.4*m0r[3]; 
  BigAEM(6,8) = 0.2981423969999719*m0r[10]+0.4391550328268399*m0r[3]; 
  BigAEM(6,10) = 0.2981423969999719*m0r[8]+0.3927922024247863*m0r[7]+0.4391550328268399*m0r[1]; 
  BigAEM(6,11) = 0.3927922024247863*m0r[7]; 
  BigAEM(7,0) = 0.5*m0r[7]; 
  BigAEM(7,1) = 0.5000000000000001*m0r[5]; 
  BigAEM(7,2) = 0.4391550328268399*m0r[11]+0.447213595499958*m0r[3]; 
  BigAEM(7,3) = 0.4391550328268399*m0r[9]+0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(7,4) = 0.4472135954999579*m0r[7]; 
  BigAEM(7,5) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(7,6) = 0.3927922024247863*m0r[11]+0.3927922024247863*m0r[10]+0.4*m0r[3]; 
  BigAEM(7,7) = 0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(7,9) = 0.2981423969999719*m0r[11]+0.4391550328268399*m0r[3]; 
  BigAEM(7,10) = 0.3927922024247863*m0r[6]; 
  BigAEM(7,11) = 0.2981423969999719*m0r[9]+0.3927922024247863*m0r[6]+0.4391550328268399*m0r[2]; 
  BigAEM(8,0) = 0.5*m0r[8]; 
  BigAEM(8,1) = 0.4391550328268398*m0r[4]; 
  BigAEM(8,2) = 0.5*m0r[10]; 
  BigAEM(8,3) = 0.4391550328268399*m0r[6]; 
  BigAEM(8,4) = 0.2981423969999719*m0r[8]+0.4391550328268398*m0r[1]; 
  BigAEM(8,6) = 0.2981423969999719*m0r[10]+0.4391550328268399*m0r[3]; 
  BigAEM(8,8) = 0.2981423969999719*m0r[4]+0.5*m0r[0]; 
  BigAEM(8,10) = 0.2981423969999719*m0r[6]+0.5*m0r[2]; 
  BigAEM(9,0) = 0.5*m0r[9]; 
  BigAEM(9,1) = 0.5*m0r[11]; 
  BigAEM(9,2) = 0.4391550328268398*m0r[5]; 
  BigAEM(9,3) = 0.4391550328268399*m0r[7]; 
  BigAEM(9,5) = 0.2981423969999719*m0r[9]+0.4391550328268398*m0r[2]; 
  BigAEM(9,7) = 0.2981423969999719*m0r[11]+0.4391550328268399*m0r[3]; 
  BigAEM(9,9) = 0.2981423969999719*m0r[5]+0.5*m0r[0]; 
  BigAEM(9,11) = 0.2981423969999719*m0r[7]+0.5*m0r[1]; 
  BigAEM(10,0) = 0.5*m0r[10]; 
  BigAEM(10,1) = 0.4391550328268399*m0r[6]; 
  BigAEM(10,2) = 0.5*m0r[8]; 
  BigAEM(10,3) = 0.4391550328268399*m0r[4]; 
  BigAEM(10,4) = 0.2981423969999719*m0r[10]+0.4391550328268399*m0r[3]; 
  BigAEM(10,5) = 0.4472135954999579*m0r[10]; 
  BigAEM(10,6) = 0.2981423969999719*m0r[8]+0.3927922024247863*m0r[7]+0.4391550328268399*m0r[1]; 
  BigAEM(10,7) = 0.3927922024247863*m0r[6]; 
  BigAEM(10,8) = 0.2981423969999719*m0r[6]+0.5*m0r[2]; 
  BigAEM(10,10) = 0.4472135954999579*m0r[5]+0.2981423969999719*m0r[4]+0.5*m0r[0]; 
  BigAEM(11,0) = 0.5*m0r[11]; 
  BigAEM(11,1) = 0.5*m0r[9]; 
  BigAEM(11,2) = 0.4391550328268399*m0r[7]; 
  BigAEM(11,3) = 0.4391550328268399*m0r[5]; 
  BigAEM(11,4) = 0.4472135954999579*m0r[11]; 
  BigAEM(11,5) = 0.2981423969999719*m0r[11]+0.4391550328268399*m0r[3]; 
  BigAEM(11,6) = 0.3927922024247863*m0r[7]; 
  BigAEM(11,7) = 0.2981423969999719*m0r[9]+0.3927922024247863*m0r[6]+0.4391550328268399*m0r[2]; 
  BigAEM(11,9) = 0.2981423969999719*m0r[7]+0.5*m0r[1]; 
  BigAEM(11,11) = 0.2981423969999719*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,24) = -0.5*cM[0]; 
  BigAEM(0,25) = -0.5*cM[1]; 
  BigAEM(0,26) = -0.5*cM[2]; 
  BigAEM(0,27) = -0.5*cM[3]; 
  BigAEM(0,28) = -0.5*cM[4]; 
  BigAEM(0,29) = -0.5*cM[5]; 
  BigAEM(0,30) = -0.5*cM[6]; 
  BigAEM(0,31) = -0.5*cM[7]; 
  BigAEM(0,32) = -0.5*cM[8]; 
  BigAEM(0,33) = -0.5*cM[9]; 
  BigAEM(0,34) = -0.5*cM[10]; 
  BigAEM(0,35) = -0.5*cM[11]; 
  BigAEM(1,24) = -0.5*cM[1]; 
  BigAEM(1,25) = (-0.4472135954999579*cM[4])-0.5*cM[0]; 
  BigAEM(1,26) = -0.5*cM[3]; 
  BigAEM(1,27) = (-0.447213595499958*cM[6])-0.5*cM[2]; 
  BigAEM(1,28) = (-0.4391550328268398*cM[8])-0.4472135954999579*cM[1]; 
  BigAEM(1,29) = -0.5000000000000001*cM[7]; 
  BigAEM(1,30) = (-0.4391550328268399*cM[10])-0.447213595499958*cM[3]; 
  BigAEM(1,31) = -0.5000000000000001*cM[5]; 
  BigAEM(1,32) = -0.4391550328268398*cM[4]; 
  BigAEM(1,33) = -0.5*cM[11]; 
  BigAEM(1,34) = -0.4391550328268399*cM[6]; 
  BigAEM(1,35) = -0.5*cM[9]; 
  BigAEM(2,24) = -0.5*cM[2]; 
  BigAEM(2,25) = -0.5*cM[3]; 
  BigAEM(2,26) = (-0.4472135954999579*cM[5])-0.5*cM[0]; 
  BigAEM(2,27) = (-0.447213595499958*cM[7])-0.5*cM[1]; 
  BigAEM(2,28) = -0.5000000000000001*cM[6]; 
  BigAEM(2,29) = (-0.4391550328268398*cM[9])-0.4472135954999579*cM[2]; 
  BigAEM(2,30) = -0.5000000000000001*cM[4]; 
  BigAEM(2,31) = (-0.4391550328268399*cM[11])-0.447213595499958*cM[3]; 
  BigAEM(2,32) = -0.5*cM[10]; 
  BigAEM(2,33) = -0.4391550328268398*cM[5]; 
  BigAEM(2,34) = -0.5*cM[8]; 
  BigAEM(2,35) = -0.4391550328268399*cM[7]; 
  BigAEM(3,24) = -0.5*cM[3]; 
  BigAEM(3,25) = (-0.447213595499958*cM[6])-0.5*cM[2]; 
  BigAEM(3,26) = (-0.447213595499958*cM[7])-0.5*cM[1]; 
  BigAEM(3,27) = (-0.4472135954999579*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(3,28) = (-0.4391550328268399*cM[10])-0.4472135954999579*cM[3]; 
  BigAEM(3,29) = (-0.4391550328268399*cM[11])-0.4472135954999579*cM[3]; 
  BigAEM(3,30) = (-0.4391550328268399*cM[8])-0.4*cM[7]-0.447213595499958*cM[1]; 
  BigAEM(3,31) = (-0.4391550328268399*cM[9])-0.4*cM[6]-0.447213595499958*cM[2]; 
  BigAEM(3,32) = -0.4391550328268399*cM[6]; 
  BigAEM(3,33) = -0.4391550328268399*cM[7]; 
  BigAEM(3,34) = -0.4391550328268399*cM[4]; 
  BigAEM(3,35) = -0.4391550328268399*cM[5]; 
  BigAEM(4,24) = -0.5*cM[4]; 
  BigAEM(4,25) = (-0.4391550328268398*cM[8])-0.4472135954999579*cM[1]; 
  BigAEM(4,26) = -0.5000000000000001*cM[6]; 
  BigAEM(4,27) = (-0.4391550328268399*cM[10])-0.4472135954999579*cM[3]; 
  BigAEM(4,28) = (-0.31943828249997*cM[4])-0.5*cM[0]; 
  BigAEM(4,30) = (-0.31943828249997*cM[6])-0.5000000000000001*cM[2]; 
  BigAEM(4,31) = -0.4472135954999579*cM[7]; 
  BigAEM(4,32) = (-0.2981423969999719*cM[8])-0.4391550328268398*cM[1]; 
  BigAEM(4,34) = (-0.2981423969999719*cM[10])-0.4391550328268399*cM[3]; 
  BigAEM(4,35) = -0.4472135954999579*cM[11]; 
  BigAEM(5,24) = -0.5*cM[5]; 
  BigAEM(5,25) = -0.5000000000000001*cM[7]; 
  BigAEM(5,26) = (-0.4391550328268398*cM[9])-0.4472135954999579*cM[2]; 
  BigAEM(5,27) = (-0.4391550328268399*cM[11])-0.4472135954999579*cM[3]; 
  BigAEM(5,29) = (-0.31943828249997*cM[5])-0.5*cM[0]; 
  BigAEM(5,30) = -0.4472135954999579*cM[6]; 
  BigAEM(5,31) = (-0.31943828249997*cM[7])-0.5000000000000001*cM[1]; 
  BigAEM(5,33) = (-0.2981423969999719*cM[9])-0.4391550328268398*cM[2]; 
  BigAEM(5,34) = -0.4472135954999579*cM[10]; 
  BigAEM(5,35) = (-0.2981423969999719*cM[11])-0.4391550328268399*cM[3]; 
  BigAEM(6,24) = -0.5*cM[6]; 
  BigAEM(6,25) = (-0.4391550328268399*cM[10])-0.447213595499958*cM[3]; 
  BigAEM(6,26) = -0.5000000000000001*cM[4]; 
  BigAEM(6,27) = (-0.4391550328268399*cM[8])-0.4*cM[7]-0.447213595499958*cM[1]; 
  BigAEM(6,28) = (-0.31943828249997*cM[6])-0.5000000000000001*cM[2]; 
  BigAEM(6,29) = -0.4472135954999579*cM[6]; 
  BigAEM(6,30) = (-0.4472135954999579*cM[5])-0.31943828249997*cM[4]-0.5*cM[0]; 
  BigAEM(6,31) = (-0.3927922024247863*cM[11])-0.3927922024247863*cM[10]-0.4*cM[3]; 
  BigAEM(6,32) = (-0.2981423969999719*cM[10])-0.4391550328268399*cM[3]; 
  BigAEM(6,34) = (-0.2981423969999719*cM[8])-0.3927922024247863*cM[7]-0.4391550328268399*cM[1]; 
  BigAEM(6,35) = -0.3927922024247863*cM[7]; 
  BigAEM(7,24) = -0.5*cM[7]; 
  BigAEM(7,25) = -0.5000000000000001*cM[5]; 
  BigAEM(7,26) = (-0.4391550328268399*cM[11])-0.447213595499958*cM[3]; 
  BigAEM(7,27) = (-0.4391550328268399*cM[9])-0.4*cM[6]-0.447213595499958*cM[2]; 
  BigAEM(7,28) = -0.4472135954999579*cM[7]; 
  BigAEM(7,29) = (-0.31943828249997*cM[7])-0.5000000000000001*cM[1]; 
  BigAEM(7,30) = (-0.3927922024247863*cM[11])-0.3927922024247863*cM[10]-0.4*cM[3]; 
  BigAEM(7,31) = (-0.31943828249997*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(7,33) = (-0.2981423969999719*cM[11])-0.4391550328268399*cM[3]; 
  BigAEM(7,34) = -0.3927922024247863*cM[6]; 
  BigAEM(7,35) = (-0.2981423969999719*cM[9])-0.3927922024247863*cM[6]-0.4391550328268399*cM[2]; 
  BigAEM(8,24) = -0.5*cM[8]; 
  BigAEM(8,25) = -0.4391550328268398*cM[4]; 
  BigAEM(8,26) = -0.5*cM[10]; 
  BigAEM(8,27) = -0.4391550328268399*cM[6]; 
  BigAEM(8,28) = (-0.2981423969999719*cM[8])-0.4391550328268398*cM[1]; 
  BigAEM(8,30) = (-0.2981423969999719*cM[10])-0.4391550328268399*cM[3]; 
  BigAEM(8,32) = (-0.2981423969999719*cM[4])-0.5*cM[0]; 
  BigAEM(8,34) = (-0.2981423969999719*cM[6])-0.5*cM[2]; 
  BigAEM(9,24) = -0.5*cM[9]; 
  BigAEM(9,25) = -0.5*cM[11]; 
  BigAEM(9,26) = -0.4391550328268398*cM[5]; 
  BigAEM(9,27) = -0.4391550328268399*cM[7]; 
  BigAEM(9,29) = (-0.2981423969999719*cM[9])-0.4391550328268398*cM[2]; 
  BigAEM(9,31) = (-0.2981423969999719*cM[11])-0.4391550328268399*cM[3]; 
  BigAEM(9,33) = (-0.2981423969999719*cM[5])-0.5*cM[0]; 
  BigAEM(9,35) = (-0.2981423969999719*cM[7])-0.5*cM[1]; 
  BigAEM(10,24) = -0.5*cM[10]; 
  BigAEM(10,25) = -0.4391550328268399*cM[6]; 
  BigAEM(10,26) = -0.5*cM[8]; 
  BigAEM(10,27) = -0.4391550328268399*cM[4]; 
  BigAEM(10,28) = (-0.2981423969999719*cM[10])-0.4391550328268399*cM[3]; 
  BigAEM(10,29) = -0.4472135954999579*cM[10]; 
  BigAEM(10,30) = (-0.2981423969999719*cM[8])-0.3927922024247863*cM[7]-0.4391550328268399*cM[1]; 
  BigAEM(10,31) = -0.3927922024247863*cM[6]; 
  BigAEM(10,32) = (-0.2981423969999719*cM[6])-0.5*cM[2]; 
  BigAEM(10,34) = (-0.4472135954999579*cM[5])-0.2981423969999719*cM[4]-0.5*cM[0]; 
  BigAEM(11,24) = -0.5*cM[11]; 
  BigAEM(11,25) = -0.5*cM[9]; 
  BigAEM(11,26) = -0.4391550328268399*cM[7]; 
  BigAEM(11,27) = -0.4391550328268399*cM[5]; 
  BigAEM(11,28) = -0.4472135954999579*cM[11]; 
  BigAEM(11,29) = (-0.2981423969999719*cM[11])-0.4391550328268399*cM[3]; 
  BigAEM(11,30) = -0.3927922024247863*cM[7]; 
  BigAEM(11,31) = (-0.2981423969999719*cM[9])-0.3927922024247863*cM[6]-0.4391550328268399*cM[2]; 
  BigAEM(11,33) = (-0.2981423969999719*cM[7])-0.5*cM[1]; 
  BigAEM(11,35) = (-0.2981423969999719*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(24,0) = 0.5*m1r[0]; 
  BigAEM(24,1) = 0.5*m1r[1]; 
  BigAEM(24,2) = 0.5*m1r[2]; 
  BigAEM(24,3) = 0.5*m1r[3]; 
  BigAEM(24,4) = 0.5*m1r[4]; 
  BigAEM(24,5) = 0.5*m1r[5]; 
  BigAEM(24,6) = 0.5*m1r[6]; 
  BigAEM(24,7) = 0.5*m1r[7]; 
  BigAEM(24,8) = 0.5*m1r[8]; 
  BigAEM(24,9) = 0.5*m1r[9]; 
  BigAEM(24,10) = 0.5*m1r[10]; 
  BigAEM(24,11) = 0.5*m1r[11]; 
  BigAEM(25,0) = 0.5*m1r[1]; 
  BigAEM(25,1) = 0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(25,2) = 0.5*m1r[3]; 
  BigAEM(25,3) = 0.447213595499958*m1r[6]+0.5*m1r[2]; 
  BigAEM(25,4) = 0.4391550328268398*m1r[8]+0.4472135954999579*m1r[1]; 
  BigAEM(25,5) = 0.5000000000000001*m1r[7]; 
  BigAEM(25,6) = 0.4391550328268399*m1r[10]+0.447213595499958*m1r[3]; 
  BigAEM(25,7) = 0.5000000000000001*m1r[5]; 
  BigAEM(25,8) = 0.4391550328268398*m1r[4]; 
  BigAEM(25,9) = 0.5*m1r[11]; 
  BigAEM(25,10) = 0.4391550328268399*m1r[6]; 
  BigAEM(25,11) = 0.5*m1r[9]; 
  BigAEM(26,0) = 0.5*m1r[2]; 
  BigAEM(26,1) = 0.5*m1r[3]; 
  BigAEM(26,2) = 0.4472135954999579*m1r[5]+0.5*m1r[0]; 
  BigAEM(26,3) = 0.447213595499958*m1r[7]+0.5*m1r[1]; 
  BigAEM(26,4) = 0.5000000000000001*m1r[6]; 
  BigAEM(26,5) = 0.4391550328268398*m1r[9]+0.4472135954999579*m1r[2]; 
  BigAEM(26,6) = 0.5000000000000001*m1r[4]; 
  BigAEM(26,7) = 0.4391550328268399*m1r[11]+0.447213595499958*m1r[3]; 
  BigAEM(26,8) = 0.5*m1r[10]; 
  BigAEM(26,9) = 0.4391550328268398*m1r[5]; 
  BigAEM(26,10) = 0.5*m1r[8]; 
  BigAEM(26,11) = 0.4391550328268399*m1r[7]; 
  BigAEM(27,0) = 0.5*m1r[3]; 
  BigAEM(27,1) = 0.447213595499958*m1r[6]+0.5*m1r[2]; 
  BigAEM(27,2) = 0.447213595499958*m1r[7]+0.5*m1r[1]; 
  BigAEM(27,3) = 0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(27,4) = 0.4391550328268399*m1r[10]+0.4472135954999579*m1r[3]; 
  BigAEM(27,5) = 0.4391550328268399*m1r[11]+0.4472135954999579*m1r[3]; 
  BigAEM(27,6) = 0.4391550328268399*m1r[8]+0.4*m1r[7]+0.447213595499958*m1r[1]; 
  BigAEM(27,7) = 0.4391550328268399*m1r[9]+0.4*m1r[6]+0.447213595499958*m1r[2]; 
  BigAEM(27,8) = 0.4391550328268399*m1r[6]; 
  BigAEM(27,9) = 0.4391550328268399*m1r[7]; 
  BigAEM(27,10) = 0.4391550328268399*m1r[4]; 
  BigAEM(27,11) = 0.4391550328268399*m1r[5]; 
  BigAEM(28,0) = 0.5*m1r[4]; 
  BigAEM(28,1) = 0.4391550328268398*m1r[8]+0.4472135954999579*m1r[1]; 
  BigAEM(28,2) = 0.5000000000000001*m1r[6]; 
  BigAEM(28,3) = 0.4391550328268399*m1r[10]+0.4472135954999579*m1r[3]; 
  BigAEM(28,4) = 0.31943828249997*m1r[4]+0.5*m1r[0]; 
  BigAEM(28,6) = 0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]; 
  BigAEM(28,7) = 0.4472135954999579*m1r[7]; 
  BigAEM(28,8) = 0.2981423969999719*m1r[8]+0.4391550328268398*m1r[1]; 
  BigAEM(28,10) = 0.2981423969999719*m1r[10]+0.4391550328268399*m1r[3]; 
  BigAEM(28,11) = 0.4472135954999579*m1r[11]; 
  BigAEM(29,0) = 0.5*m1r[5]; 
  BigAEM(29,1) = 0.5000000000000001*m1r[7]; 
  BigAEM(29,2) = 0.4391550328268398*m1r[9]+0.4472135954999579*m1r[2]; 
  BigAEM(29,3) = 0.4391550328268399*m1r[11]+0.4472135954999579*m1r[3]; 
  BigAEM(29,5) = 0.31943828249997*m1r[5]+0.5*m1r[0]; 
  BigAEM(29,6) = 0.4472135954999579*m1r[6]; 
  BigAEM(29,7) = 0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]; 
  BigAEM(29,9) = 0.2981423969999719*m1r[9]+0.4391550328268398*m1r[2]; 
  BigAEM(29,10) = 0.4472135954999579*m1r[10]; 
  BigAEM(29,11) = 0.2981423969999719*m1r[11]+0.4391550328268399*m1r[3]; 
  BigAEM(30,0) = 0.5*m1r[6]; 
  BigAEM(30,1) = 0.4391550328268399*m1r[10]+0.447213595499958*m1r[3]; 
  BigAEM(30,2) = 0.5000000000000001*m1r[4]; 
  BigAEM(30,3) = 0.4391550328268399*m1r[8]+0.4*m1r[7]+0.447213595499958*m1r[1]; 
  BigAEM(30,4) = 0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]; 
  BigAEM(30,5) = 0.4472135954999579*m1r[6]; 
  BigAEM(30,6) = 0.4472135954999579*m1r[5]+0.31943828249997*m1r[4]+0.5*m1r[0]; 
  BigAEM(30,7) = 0.3927922024247863*m1r[11]+0.3927922024247863*m1r[10]+0.4*m1r[3]; 
  BigAEM(30,8) = 0.2981423969999719*m1r[10]+0.4391550328268399*m1r[3]; 
  BigAEM(30,10) = 0.2981423969999719*m1r[8]+0.3927922024247863*m1r[7]+0.4391550328268399*m1r[1]; 
  BigAEM(30,11) = 0.3927922024247863*m1r[7]; 
  BigAEM(31,0) = 0.5*m1r[7]; 
  BigAEM(31,1) = 0.5000000000000001*m1r[5]; 
  BigAEM(31,2) = 0.4391550328268399*m1r[11]+0.447213595499958*m1r[3]; 
  BigAEM(31,3) = 0.4391550328268399*m1r[9]+0.4*m1r[6]+0.447213595499958*m1r[2]; 
  BigAEM(31,4) = 0.4472135954999579*m1r[7]; 
  BigAEM(31,5) = 0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]; 
  BigAEM(31,6) = 0.3927922024247863*m1r[11]+0.3927922024247863*m1r[10]+0.4*m1r[3]; 
  BigAEM(31,7) = 0.31943828249997*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  BigAEM(31,9) = 0.2981423969999719*m1r[11]+0.4391550328268399*m1r[3]; 
  BigAEM(31,10) = 0.3927922024247863*m1r[6]; 
  BigAEM(31,11) = 0.2981423969999719*m1r[9]+0.3927922024247863*m1r[6]+0.4391550328268399*m1r[2]; 
  BigAEM(32,0) = 0.5*m1r[8]; 
  BigAEM(32,1) = 0.4391550328268398*m1r[4]; 
  BigAEM(32,2) = 0.5*m1r[10]; 
  BigAEM(32,3) = 0.4391550328268399*m1r[6]; 
  BigAEM(32,4) = 0.2981423969999719*m1r[8]+0.4391550328268398*m1r[1]; 
  BigAEM(32,6) = 0.2981423969999719*m1r[10]+0.4391550328268399*m1r[3]; 
  BigAEM(32,8) = 0.2981423969999719*m1r[4]+0.5*m1r[0]; 
  BigAEM(32,10) = 0.2981423969999719*m1r[6]+0.5*m1r[2]; 
  BigAEM(33,0) = 0.5*m1r[9]; 
  BigAEM(33,1) = 0.5*m1r[11]; 
  BigAEM(33,2) = 0.4391550328268398*m1r[5]; 
  BigAEM(33,3) = 0.4391550328268399*m1r[7]; 
  BigAEM(33,5) = 0.2981423969999719*m1r[9]+0.4391550328268398*m1r[2]; 
  BigAEM(33,7) = 0.2981423969999719*m1r[11]+0.4391550328268399*m1r[3]; 
  BigAEM(33,9) = 0.2981423969999719*m1r[5]+0.5*m1r[0]; 
  BigAEM(33,11) = 0.2981423969999719*m1r[7]+0.5*m1r[1]; 
  BigAEM(34,0) = 0.5*m1r[10]; 
  BigAEM(34,1) = 0.4391550328268399*m1r[6]; 
  BigAEM(34,2) = 0.5*m1r[8]; 
  BigAEM(34,3) = 0.4391550328268399*m1r[4]; 
  BigAEM(34,4) = 0.2981423969999719*m1r[10]+0.4391550328268399*m1r[3]; 
  BigAEM(34,5) = 0.4472135954999579*m1r[10]; 
  BigAEM(34,6) = 0.2981423969999719*m1r[8]+0.3927922024247863*m1r[7]+0.4391550328268399*m1r[1]; 
  BigAEM(34,7) = 0.3927922024247863*m1r[6]; 
  BigAEM(34,8) = 0.2981423969999719*m1r[6]+0.5*m1r[2]; 
  BigAEM(34,10) = 0.4472135954999579*m1r[5]+0.2981423969999719*m1r[4]+0.5*m1r[0]; 
  BigAEM(35,0) = 0.5*m1r[11]; 
  BigAEM(35,1) = 0.5*m1r[9]; 
  BigAEM(35,2) = 0.4391550328268399*m1r[7]; 
  BigAEM(35,3) = 0.4391550328268399*m1r[5]; 
  BigAEM(35,4) = 0.4472135954999579*m1r[11]; 
  BigAEM(35,5) = 0.2981423969999719*m1r[11]+0.4391550328268399*m1r[3]; 
  BigAEM(35,6) = 0.3927922024247863*m1r[7]; 
  BigAEM(35,7) = 0.2981423969999719*m1r[9]+0.3927922024247863*m1r[6]+0.4391550328268399*m1r[2]; 
  BigAEM(35,9) = 0.2981423969999719*m1r[7]+0.5*m1r[1]; 
  BigAEM(35,11) = 0.2981423969999719*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(12,12) = 0.5*m0r[0]; 
  BigAEM(12,13) = 0.5*m0r[1]; 
  BigAEM(12,14) = 0.5*m0r[2]; 
  BigAEM(12,15) = 0.5*m0r[3]; 
  BigAEM(12,16) = 0.5*m0r[4]; 
  BigAEM(12,17) = 0.5*m0r[5]; 
  BigAEM(12,18) = 0.5*m0r[6]; 
  BigAEM(12,19) = 0.5*m0r[7]; 
  BigAEM(12,20) = 0.5*m0r[8]; 
  BigAEM(12,21) = 0.5*m0r[9]; 
  BigAEM(12,22) = 0.5*m0r[10]; 
  BigAEM(12,23) = 0.5*m0r[11]; 
  BigAEM(13,12) = 0.5*m0r[1]; 
  BigAEM(13,13) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(13,14) = 0.5*m0r[3]; 
  BigAEM(13,15) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(13,16) = 0.4391550328268398*m0r[8]+0.4472135954999579*m0r[1]; 
  BigAEM(13,17) = 0.5000000000000001*m0r[7]; 
  BigAEM(13,18) = 0.4391550328268399*m0r[10]+0.447213595499958*m0r[3]; 
  BigAEM(13,19) = 0.5000000000000001*m0r[5]; 
  BigAEM(13,20) = 0.4391550328268398*m0r[4]; 
  BigAEM(13,21) = 0.5*m0r[11]; 
  BigAEM(13,22) = 0.4391550328268399*m0r[6]; 
  BigAEM(13,23) = 0.5*m0r[9]; 
  BigAEM(14,12) = 0.5*m0r[2]; 
  BigAEM(14,13) = 0.5*m0r[3]; 
  BigAEM(14,14) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  BigAEM(14,15) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(14,16) = 0.5000000000000001*m0r[6]; 
  BigAEM(14,17) = 0.4391550328268398*m0r[9]+0.4472135954999579*m0r[2]; 
  BigAEM(14,18) = 0.5000000000000001*m0r[4]; 
  BigAEM(14,19) = 0.4391550328268399*m0r[11]+0.447213595499958*m0r[3]; 
  BigAEM(14,20) = 0.5*m0r[10]; 
  BigAEM(14,21) = 0.4391550328268398*m0r[5]; 
  BigAEM(14,22) = 0.5*m0r[8]; 
  BigAEM(14,23) = 0.4391550328268399*m0r[7]; 
  BigAEM(15,12) = 0.5*m0r[3]; 
  BigAEM(15,13) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  BigAEM(15,14) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  BigAEM(15,15) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(15,16) = 0.4391550328268399*m0r[10]+0.4472135954999579*m0r[3]; 
  BigAEM(15,17) = 0.4391550328268399*m0r[11]+0.4472135954999579*m0r[3]; 
  BigAEM(15,18) = 0.4391550328268399*m0r[8]+0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(15,19) = 0.4391550328268399*m0r[9]+0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(15,20) = 0.4391550328268399*m0r[6]; 
  BigAEM(15,21) = 0.4391550328268399*m0r[7]; 
  BigAEM(15,22) = 0.4391550328268399*m0r[4]; 
  BigAEM(15,23) = 0.4391550328268399*m0r[5]; 
  BigAEM(16,12) = 0.5*m0r[4]; 
  BigAEM(16,13) = 0.4391550328268398*m0r[8]+0.4472135954999579*m0r[1]; 
  BigAEM(16,14) = 0.5000000000000001*m0r[6]; 
  BigAEM(16,15) = 0.4391550328268399*m0r[10]+0.4472135954999579*m0r[3]; 
  BigAEM(16,16) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(16,18) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(16,19) = 0.4472135954999579*m0r[7]; 
  BigAEM(16,20) = 0.2981423969999719*m0r[8]+0.4391550328268398*m0r[1]; 
  BigAEM(16,22) = 0.2981423969999719*m0r[10]+0.4391550328268399*m0r[3]; 
  BigAEM(16,23) = 0.4472135954999579*m0r[11]; 
  BigAEM(17,12) = 0.5*m0r[5]; 
  BigAEM(17,13) = 0.5000000000000001*m0r[7]; 
  BigAEM(17,14) = 0.4391550328268398*m0r[9]+0.4472135954999579*m0r[2]; 
  BigAEM(17,15) = 0.4391550328268399*m0r[11]+0.4472135954999579*m0r[3]; 
  BigAEM(17,17) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
  BigAEM(17,18) = 0.4472135954999579*m0r[6]; 
  BigAEM(17,19) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(17,21) = 0.2981423969999719*m0r[9]+0.4391550328268398*m0r[2]; 
  BigAEM(17,22) = 0.4472135954999579*m0r[10]; 
  BigAEM(17,23) = 0.2981423969999719*m0r[11]+0.4391550328268399*m0r[3]; 
  BigAEM(18,12) = 0.5*m0r[6]; 
  BigAEM(18,13) = 0.4391550328268399*m0r[10]+0.447213595499958*m0r[3]; 
  BigAEM(18,14) = 0.5000000000000001*m0r[4]; 
  BigAEM(18,15) = 0.4391550328268399*m0r[8]+0.4*m0r[7]+0.447213595499958*m0r[1]; 
  BigAEM(18,16) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  BigAEM(18,17) = 0.4472135954999579*m0r[6]; 
  BigAEM(18,18) = 0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]; 
  BigAEM(18,19) = 0.3927922024247863*m0r[11]+0.3927922024247863*m0r[10]+0.4*m0r[3]; 
  BigAEM(18,20) = 0.2981423969999719*m0r[10]+0.4391550328268399*m0r[3]; 
  BigAEM(18,22) = 0.2981423969999719*m0r[8]+0.3927922024247863*m0r[7]+0.4391550328268399*m0r[1]; 
  BigAEM(18,23) = 0.3927922024247863*m0r[7]; 
  BigAEM(19,12) = 0.5*m0r[7]; 
  BigAEM(19,13) = 0.5000000000000001*m0r[5]; 
  BigAEM(19,14) = 0.4391550328268399*m0r[11]+0.447213595499958*m0r[3]; 
  BigAEM(19,15) = 0.4391550328268399*m0r[9]+0.4*m0r[6]+0.447213595499958*m0r[2]; 
  BigAEM(19,16) = 0.4472135954999579*m0r[7]; 
  BigAEM(19,17) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  BigAEM(19,18) = 0.3927922024247863*m0r[11]+0.3927922024247863*m0r[10]+0.4*m0r[3]; 
  BigAEM(19,19) = 0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  BigAEM(19,21) = 0.2981423969999719*m0r[11]+0.4391550328268399*m0r[3]; 
  BigAEM(19,22) = 0.3927922024247863*m0r[6]; 
  BigAEM(19,23) = 0.2981423969999719*m0r[9]+0.3927922024247863*m0r[6]+0.4391550328268399*m0r[2]; 
  BigAEM(20,12) = 0.5*m0r[8]; 
  BigAEM(20,13) = 0.4391550328268398*m0r[4]; 
  BigAEM(20,14) = 0.5*m0r[10]; 
  BigAEM(20,15) = 0.4391550328268399*m0r[6]; 
  BigAEM(20,16) = 0.2981423969999719*m0r[8]+0.4391550328268398*m0r[1]; 
  BigAEM(20,18) = 0.2981423969999719*m0r[10]+0.4391550328268399*m0r[3]; 
  BigAEM(20,20) = 0.2981423969999719*m0r[4]+0.5*m0r[0]; 
  BigAEM(20,22) = 0.2981423969999719*m0r[6]+0.5*m0r[2]; 
  BigAEM(21,12) = 0.5*m0r[9]; 
  BigAEM(21,13) = 0.5*m0r[11]; 
  BigAEM(21,14) = 0.4391550328268398*m0r[5]; 
  BigAEM(21,15) = 0.4391550328268399*m0r[7]; 
  BigAEM(21,17) = 0.2981423969999719*m0r[9]+0.4391550328268398*m0r[2]; 
  BigAEM(21,19) = 0.2981423969999719*m0r[11]+0.4391550328268399*m0r[3]; 
  BigAEM(21,21) = 0.2981423969999719*m0r[5]+0.5*m0r[0]; 
  BigAEM(21,23) = 0.2981423969999719*m0r[7]+0.5*m0r[1]; 
  BigAEM(22,12) = 0.5*m0r[10]; 
  BigAEM(22,13) = 0.4391550328268399*m0r[6]; 
  BigAEM(22,14) = 0.5*m0r[8]; 
  BigAEM(22,15) = 0.4391550328268399*m0r[4]; 
  BigAEM(22,16) = 0.2981423969999719*m0r[10]+0.4391550328268399*m0r[3]; 
  BigAEM(22,17) = 0.4472135954999579*m0r[10]; 
  BigAEM(22,18) = 0.2981423969999719*m0r[8]+0.3927922024247863*m0r[7]+0.4391550328268399*m0r[1]; 
  BigAEM(22,19) = 0.3927922024247863*m0r[6]; 
  BigAEM(22,20) = 0.2981423969999719*m0r[6]+0.5*m0r[2]; 
  BigAEM(22,22) = 0.4472135954999579*m0r[5]+0.2981423969999719*m0r[4]+0.5*m0r[0]; 
  BigAEM(23,12) = 0.5*m0r[11]; 
  BigAEM(23,13) = 0.5*m0r[9]; 
  BigAEM(23,14) = 0.4391550328268399*m0r[7]; 
  BigAEM(23,15) = 0.4391550328268399*m0r[5]; 
  BigAEM(23,16) = 0.4472135954999579*m0r[11]; 
  BigAEM(23,17) = 0.2981423969999719*m0r[11]+0.4391550328268399*m0r[3]; 
  BigAEM(23,18) = 0.3927922024247863*m0r[7]; 
  BigAEM(23,19) = 0.2981423969999719*m0r[9]+0.3927922024247863*m0r[6]+0.4391550328268399*m0r[2]; 
  BigAEM(23,21) = 0.2981423969999719*m0r[7]+0.5*m0r[1]; 
  BigAEM(23,23) = 0.2981423969999719*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(12,24) = -0.5*cM[12]; 
  BigAEM(12,25) = -0.5*cM[13]; 
  BigAEM(12,26) = -0.5*cM[14]; 
  BigAEM(12,27) = -0.5*cM[15]; 
  BigAEM(12,28) = -0.5*cM[16]; 
  BigAEM(12,29) = -0.5*cM[17]; 
  BigAEM(12,30) = -0.5*cM[18]; 
  BigAEM(12,31) = -0.5*cM[19]; 
  BigAEM(12,32) = -0.5*cM[20]; 
  BigAEM(12,33) = -0.5*cM[21]; 
  BigAEM(12,34) = -0.5*cM[22]; 
  BigAEM(12,35) = -0.5*cM[23]; 
  BigAEM(13,24) = -0.5*cM[13]; 
  BigAEM(13,25) = (-0.4472135954999579*cM[16])-0.5*cM[12]; 
  BigAEM(13,26) = -0.5*cM[15]; 
  BigAEM(13,27) = (-0.447213595499958*cM[18])-0.5*cM[14]; 
  BigAEM(13,28) = (-0.4391550328268398*cM[20])-0.4472135954999579*cM[13]; 
  BigAEM(13,29) = -0.5000000000000001*cM[19]; 
  BigAEM(13,30) = (-0.4391550328268399*cM[22])-0.447213595499958*cM[15]; 
  BigAEM(13,31) = -0.5000000000000001*cM[17]; 
  BigAEM(13,32) = -0.4391550328268398*cM[16]; 
  BigAEM(13,33) = -0.5*cM[23]; 
  BigAEM(13,34) = -0.4391550328268399*cM[18]; 
  BigAEM(13,35) = -0.5*cM[21]; 
  BigAEM(14,24) = -0.5*cM[14]; 
  BigAEM(14,25) = -0.5*cM[15]; 
  BigAEM(14,26) = (-0.4472135954999579*cM[17])-0.5*cM[12]; 
  BigAEM(14,27) = (-0.447213595499958*cM[19])-0.5*cM[13]; 
  BigAEM(14,28) = -0.5000000000000001*cM[18]; 
  BigAEM(14,29) = (-0.4391550328268398*cM[21])-0.4472135954999579*cM[14]; 
  BigAEM(14,30) = -0.5000000000000001*cM[16]; 
  BigAEM(14,31) = (-0.4391550328268399*cM[23])-0.447213595499958*cM[15]; 
  BigAEM(14,32) = -0.5*cM[22]; 
  BigAEM(14,33) = -0.4391550328268398*cM[17]; 
  BigAEM(14,34) = -0.5*cM[20]; 
  BigAEM(14,35) = -0.4391550328268399*cM[19]; 
  BigAEM(15,24) = -0.5*cM[15]; 
  BigAEM(15,25) = (-0.447213595499958*cM[18])-0.5*cM[14]; 
  BigAEM(15,26) = (-0.447213595499958*cM[19])-0.5*cM[13]; 
  BigAEM(15,27) = (-0.4472135954999579*cM[17])-0.4472135954999579*cM[16]-0.5*cM[12]; 
  BigAEM(15,28) = (-0.4391550328268399*cM[22])-0.4472135954999579*cM[15]; 
  BigAEM(15,29) = (-0.4391550328268399*cM[23])-0.4472135954999579*cM[15]; 
  BigAEM(15,30) = (-0.4391550328268399*cM[20])-0.4*cM[19]-0.447213595499958*cM[13]; 
  BigAEM(15,31) = (-0.4391550328268399*cM[21])-0.4*cM[18]-0.447213595499958*cM[14]; 
  BigAEM(15,32) = -0.4391550328268399*cM[18]; 
  BigAEM(15,33) = -0.4391550328268399*cM[19]; 
  BigAEM(15,34) = -0.4391550328268399*cM[16]; 
  BigAEM(15,35) = -0.4391550328268399*cM[17]; 
  BigAEM(16,24) = -0.5*cM[16]; 
  BigAEM(16,25) = (-0.4391550328268398*cM[20])-0.4472135954999579*cM[13]; 
  BigAEM(16,26) = -0.5000000000000001*cM[18]; 
  BigAEM(16,27) = (-0.4391550328268399*cM[22])-0.4472135954999579*cM[15]; 
  BigAEM(16,28) = (-0.31943828249997*cM[16])-0.5*cM[12]; 
  BigAEM(16,30) = (-0.31943828249997*cM[18])-0.5000000000000001*cM[14]; 
  BigAEM(16,31) = -0.4472135954999579*cM[19]; 
  BigAEM(16,32) = (-0.2981423969999719*cM[20])-0.4391550328268398*cM[13]; 
  BigAEM(16,34) = (-0.2981423969999719*cM[22])-0.4391550328268399*cM[15]; 
  BigAEM(16,35) = -0.4472135954999579*cM[23]; 
  BigAEM(17,24) = -0.5*cM[17]; 
  BigAEM(17,25) = -0.5000000000000001*cM[19]; 
  BigAEM(17,26) = (-0.4391550328268398*cM[21])-0.4472135954999579*cM[14]; 
  BigAEM(17,27) = (-0.4391550328268399*cM[23])-0.4472135954999579*cM[15]; 
  BigAEM(17,29) = (-0.31943828249997*cM[17])-0.5*cM[12]; 
  BigAEM(17,30) = -0.4472135954999579*cM[18]; 
  BigAEM(17,31) = (-0.31943828249997*cM[19])-0.5000000000000001*cM[13]; 
  BigAEM(17,33) = (-0.2981423969999719*cM[21])-0.4391550328268398*cM[14]; 
  BigAEM(17,34) = -0.4472135954999579*cM[22]; 
  BigAEM(17,35) = (-0.2981423969999719*cM[23])-0.4391550328268399*cM[15]; 
  BigAEM(18,24) = -0.5*cM[18]; 
  BigAEM(18,25) = (-0.4391550328268399*cM[22])-0.447213595499958*cM[15]; 
  BigAEM(18,26) = -0.5000000000000001*cM[16]; 
  BigAEM(18,27) = (-0.4391550328268399*cM[20])-0.4*cM[19]-0.447213595499958*cM[13]; 
  BigAEM(18,28) = (-0.31943828249997*cM[18])-0.5000000000000001*cM[14]; 
  BigAEM(18,29) = -0.4472135954999579*cM[18]; 
  BigAEM(18,30) = (-0.4472135954999579*cM[17])-0.31943828249997*cM[16]-0.5*cM[12]; 
  BigAEM(18,31) = (-0.3927922024247863*cM[23])-0.3927922024247863*cM[22]-0.4*cM[15]; 
  BigAEM(18,32) = (-0.2981423969999719*cM[22])-0.4391550328268399*cM[15]; 
  BigAEM(18,34) = (-0.2981423969999719*cM[20])-0.3927922024247863*cM[19]-0.4391550328268399*cM[13]; 
  BigAEM(18,35) = -0.3927922024247863*cM[19]; 
  BigAEM(19,24) = -0.5*cM[19]; 
  BigAEM(19,25) = -0.5000000000000001*cM[17]; 
  BigAEM(19,26) = (-0.4391550328268399*cM[23])-0.447213595499958*cM[15]; 
  BigAEM(19,27) = (-0.4391550328268399*cM[21])-0.4*cM[18]-0.447213595499958*cM[14]; 
  BigAEM(19,28) = -0.4472135954999579*cM[19]; 
  BigAEM(19,29) = (-0.31943828249997*cM[19])-0.5000000000000001*cM[13]; 
  BigAEM(19,30) = (-0.3927922024247863*cM[23])-0.3927922024247863*cM[22]-0.4*cM[15]; 
  BigAEM(19,31) = (-0.31943828249997*cM[17])-0.4472135954999579*cM[16]-0.5*cM[12]; 
  BigAEM(19,33) = (-0.2981423969999719*cM[23])-0.4391550328268399*cM[15]; 
  BigAEM(19,34) = -0.3927922024247863*cM[18]; 
  BigAEM(19,35) = (-0.2981423969999719*cM[21])-0.3927922024247863*cM[18]-0.4391550328268399*cM[14]; 
  BigAEM(20,24) = -0.5*cM[20]; 
  BigAEM(20,25) = -0.4391550328268398*cM[16]; 
  BigAEM(20,26) = -0.5*cM[22]; 
  BigAEM(20,27) = -0.4391550328268399*cM[18]; 
  BigAEM(20,28) = (-0.2981423969999719*cM[20])-0.4391550328268398*cM[13]; 
  BigAEM(20,30) = (-0.2981423969999719*cM[22])-0.4391550328268399*cM[15]; 
  BigAEM(20,32) = (-0.2981423969999719*cM[16])-0.5*cM[12]; 
  BigAEM(20,34) = (-0.2981423969999719*cM[18])-0.5*cM[14]; 
  BigAEM(21,24) = -0.5*cM[21]; 
  BigAEM(21,25) = -0.5*cM[23]; 
  BigAEM(21,26) = -0.4391550328268398*cM[17]; 
  BigAEM(21,27) = -0.4391550328268399*cM[19]; 
  BigAEM(21,29) = (-0.2981423969999719*cM[21])-0.4391550328268398*cM[14]; 
  BigAEM(21,31) = (-0.2981423969999719*cM[23])-0.4391550328268399*cM[15]; 
  BigAEM(21,33) = (-0.2981423969999719*cM[17])-0.5*cM[12]; 
  BigAEM(21,35) = (-0.2981423969999719*cM[19])-0.5*cM[13]; 
  BigAEM(22,24) = -0.5*cM[22]; 
  BigAEM(22,25) = -0.4391550328268399*cM[18]; 
  BigAEM(22,26) = -0.5*cM[20]; 
  BigAEM(22,27) = -0.4391550328268399*cM[16]; 
  BigAEM(22,28) = (-0.2981423969999719*cM[22])-0.4391550328268399*cM[15]; 
  BigAEM(22,29) = -0.4472135954999579*cM[22]; 
  BigAEM(22,30) = (-0.2981423969999719*cM[20])-0.3927922024247863*cM[19]-0.4391550328268399*cM[13]; 
  BigAEM(22,31) = -0.3927922024247863*cM[18]; 
  BigAEM(22,32) = (-0.2981423969999719*cM[18])-0.5*cM[14]; 
  BigAEM(22,34) = (-0.4472135954999579*cM[17])-0.2981423969999719*cM[16]-0.5*cM[12]; 
  BigAEM(23,24) = -0.5*cM[23]; 
  BigAEM(23,25) = -0.5*cM[21]; 
  BigAEM(23,26) = -0.4391550328268399*cM[19]; 
  BigAEM(23,27) = -0.4391550328268399*cM[17]; 
  BigAEM(23,28) = -0.4472135954999579*cM[23]; 
  BigAEM(23,29) = (-0.2981423969999719*cM[23])-0.4391550328268399*cM[15]; 
  BigAEM(23,30) = -0.3927922024247863*cM[19]; 
  BigAEM(23,31) = (-0.2981423969999719*cM[21])-0.3927922024247863*cM[18]-0.4391550328268399*cM[14]; 
  BigAEM(23,33) = (-0.2981423969999719*cM[19])-0.5*cM[13]; 
  BigAEM(23,35) = (-0.2981423969999719*cM[17])-0.4472135954999579*cM[16]-0.5*cM[12]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(24,12) = 0.5*m1r[12]; 
  BigAEM(24,13) = 0.5*m1r[13]; 
  BigAEM(24,14) = 0.5*m1r[14]; 
  BigAEM(24,15) = 0.5*m1r[15]; 
  BigAEM(24,16) = 0.5*m1r[16]; 
  BigAEM(24,17) = 0.5*m1r[17]; 
  BigAEM(24,18) = 0.5*m1r[18]; 
  BigAEM(24,19) = 0.5*m1r[19]; 
  BigAEM(24,20) = 0.5*m1r[20]; 
  BigAEM(24,21) = 0.5*m1r[21]; 
  BigAEM(24,22) = 0.5*m1r[22]; 
  BigAEM(24,23) = 0.5*m1r[23]; 
  BigAEM(25,12) = 0.5*m1r[13]; 
  BigAEM(25,13) = 0.4472135954999579*m1r[16]+0.5*m1r[12]; 
  BigAEM(25,14) = 0.5*m1r[15]; 
  BigAEM(25,15) = 0.447213595499958*m1r[18]+0.5*m1r[14]; 
  BigAEM(25,16) = 0.4391550328268398*m1r[20]+0.4472135954999579*m1r[13]; 
  BigAEM(25,17) = 0.5000000000000001*m1r[19]; 
  BigAEM(25,18) = 0.4391550328268399*m1r[22]+0.447213595499958*m1r[15]; 
  BigAEM(25,19) = 0.5000000000000001*m1r[17]; 
  BigAEM(25,20) = 0.4391550328268398*m1r[16]; 
  BigAEM(25,21) = 0.5*m1r[23]; 
  BigAEM(25,22) = 0.4391550328268399*m1r[18]; 
  BigAEM(25,23) = 0.5*m1r[21]; 
  BigAEM(26,12) = 0.5*m1r[14]; 
  BigAEM(26,13) = 0.5*m1r[15]; 
  BigAEM(26,14) = 0.4472135954999579*m1r[17]+0.5*m1r[12]; 
  BigAEM(26,15) = 0.447213595499958*m1r[19]+0.5*m1r[13]; 
  BigAEM(26,16) = 0.5000000000000001*m1r[18]; 
  BigAEM(26,17) = 0.4391550328268398*m1r[21]+0.4472135954999579*m1r[14]; 
  BigAEM(26,18) = 0.5000000000000001*m1r[16]; 
  BigAEM(26,19) = 0.4391550328268399*m1r[23]+0.447213595499958*m1r[15]; 
  BigAEM(26,20) = 0.5*m1r[22]; 
  BigAEM(26,21) = 0.4391550328268398*m1r[17]; 
  BigAEM(26,22) = 0.5*m1r[20]; 
  BigAEM(26,23) = 0.4391550328268399*m1r[19]; 
  BigAEM(27,12) = 0.5*m1r[15]; 
  BigAEM(27,13) = 0.447213595499958*m1r[18]+0.5*m1r[14]; 
  BigAEM(27,14) = 0.447213595499958*m1r[19]+0.5*m1r[13]; 
  BigAEM(27,15) = 0.4472135954999579*m1r[17]+0.4472135954999579*m1r[16]+0.5*m1r[12]; 
  BigAEM(27,16) = 0.4391550328268399*m1r[22]+0.4472135954999579*m1r[15]; 
  BigAEM(27,17) = 0.4391550328268399*m1r[23]+0.4472135954999579*m1r[15]; 
  BigAEM(27,18) = 0.4391550328268399*m1r[20]+0.4*m1r[19]+0.447213595499958*m1r[13]; 
  BigAEM(27,19) = 0.4391550328268399*m1r[21]+0.4*m1r[18]+0.447213595499958*m1r[14]; 
  BigAEM(27,20) = 0.4391550328268399*m1r[18]; 
  BigAEM(27,21) = 0.4391550328268399*m1r[19]; 
  BigAEM(27,22) = 0.4391550328268399*m1r[16]; 
  BigAEM(27,23) = 0.4391550328268399*m1r[17]; 
  BigAEM(28,12) = 0.5*m1r[16]; 
  BigAEM(28,13) = 0.4391550328268398*m1r[20]+0.4472135954999579*m1r[13]; 
  BigAEM(28,14) = 0.5000000000000001*m1r[18]; 
  BigAEM(28,15) = 0.4391550328268399*m1r[22]+0.4472135954999579*m1r[15]; 
  BigAEM(28,16) = 0.31943828249997*m1r[16]+0.5*m1r[12]; 
  BigAEM(28,18) = 0.31943828249997*m1r[18]+0.5000000000000001*m1r[14]; 
  BigAEM(28,19) = 0.4472135954999579*m1r[19]; 
  BigAEM(28,20) = 0.2981423969999719*m1r[20]+0.4391550328268398*m1r[13]; 
  BigAEM(28,22) = 0.2981423969999719*m1r[22]+0.4391550328268399*m1r[15]; 
  BigAEM(28,23) = 0.4472135954999579*m1r[23]; 
  BigAEM(29,12) = 0.5*m1r[17]; 
  BigAEM(29,13) = 0.5000000000000001*m1r[19]; 
  BigAEM(29,14) = 0.4391550328268398*m1r[21]+0.4472135954999579*m1r[14]; 
  BigAEM(29,15) = 0.4391550328268399*m1r[23]+0.4472135954999579*m1r[15]; 
  BigAEM(29,17) = 0.31943828249997*m1r[17]+0.5*m1r[12]; 
  BigAEM(29,18) = 0.4472135954999579*m1r[18]; 
  BigAEM(29,19) = 0.31943828249997*m1r[19]+0.5000000000000001*m1r[13]; 
  BigAEM(29,21) = 0.2981423969999719*m1r[21]+0.4391550328268398*m1r[14]; 
  BigAEM(29,22) = 0.4472135954999579*m1r[22]; 
  BigAEM(29,23) = 0.2981423969999719*m1r[23]+0.4391550328268399*m1r[15]; 
  BigAEM(30,12) = 0.5*m1r[18]; 
  BigAEM(30,13) = 0.4391550328268399*m1r[22]+0.447213595499958*m1r[15]; 
  BigAEM(30,14) = 0.5000000000000001*m1r[16]; 
  BigAEM(30,15) = 0.4391550328268399*m1r[20]+0.4*m1r[19]+0.447213595499958*m1r[13]; 
  BigAEM(30,16) = 0.31943828249997*m1r[18]+0.5000000000000001*m1r[14]; 
  BigAEM(30,17) = 0.4472135954999579*m1r[18]; 
  BigAEM(30,18) = 0.4472135954999579*m1r[17]+0.31943828249997*m1r[16]+0.5*m1r[12]; 
  BigAEM(30,19) = 0.3927922024247863*m1r[23]+0.3927922024247863*m1r[22]+0.4*m1r[15]; 
  BigAEM(30,20) = 0.2981423969999719*m1r[22]+0.4391550328268399*m1r[15]; 
  BigAEM(30,22) = 0.2981423969999719*m1r[20]+0.3927922024247863*m1r[19]+0.4391550328268399*m1r[13]; 
  BigAEM(30,23) = 0.3927922024247863*m1r[19]; 
  BigAEM(31,12) = 0.5*m1r[19]; 
  BigAEM(31,13) = 0.5000000000000001*m1r[17]; 
  BigAEM(31,14) = 0.4391550328268399*m1r[23]+0.447213595499958*m1r[15]; 
  BigAEM(31,15) = 0.4391550328268399*m1r[21]+0.4*m1r[18]+0.447213595499958*m1r[14]; 
  BigAEM(31,16) = 0.4472135954999579*m1r[19]; 
  BigAEM(31,17) = 0.31943828249997*m1r[19]+0.5000000000000001*m1r[13]; 
  BigAEM(31,18) = 0.3927922024247863*m1r[23]+0.3927922024247863*m1r[22]+0.4*m1r[15]; 
  BigAEM(31,19) = 0.31943828249997*m1r[17]+0.4472135954999579*m1r[16]+0.5*m1r[12]; 
  BigAEM(31,21) = 0.2981423969999719*m1r[23]+0.4391550328268399*m1r[15]; 
  BigAEM(31,22) = 0.3927922024247863*m1r[18]; 
  BigAEM(31,23) = 0.2981423969999719*m1r[21]+0.3927922024247863*m1r[18]+0.4391550328268399*m1r[14]; 
  BigAEM(32,12) = 0.5*m1r[20]; 
  BigAEM(32,13) = 0.4391550328268398*m1r[16]; 
  BigAEM(32,14) = 0.5*m1r[22]; 
  BigAEM(32,15) = 0.4391550328268399*m1r[18]; 
  BigAEM(32,16) = 0.2981423969999719*m1r[20]+0.4391550328268398*m1r[13]; 
  BigAEM(32,18) = 0.2981423969999719*m1r[22]+0.4391550328268399*m1r[15]; 
  BigAEM(32,20) = 0.2981423969999719*m1r[16]+0.5*m1r[12]; 
  BigAEM(32,22) = 0.2981423969999719*m1r[18]+0.5*m1r[14]; 
  BigAEM(33,12) = 0.5*m1r[21]; 
  BigAEM(33,13) = 0.5*m1r[23]; 
  BigAEM(33,14) = 0.4391550328268398*m1r[17]; 
  BigAEM(33,15) = 0.4391550328268399*m1r[19]; 
  BigAEM(33,17) = 0.2981423969999719*m1r[21]+0.4391550328268398*m1r[14]; 
  BigAEM(33,19) = 0.2981423969999719*m1r[23]+0.4391550328268399*m1r[15]; 
  BigAEM(33,21) = 0.2981423969999719*m1r[17]+0.5*m1r[12]; 
  BigAEM(33,23) = 0.2981423969999719*m1r[19]+0.5*m1r[13]; 
  BigAEM(34,12) = 0.5*m1r[22]; 
  BigAEM(34,13) = 0.4391550328268399*m1r[18]; 
  BigAEM(34,14) = 0.5*m1r[20]; 
  BigAEM(34,15) = 0.4391550328268399*m1r[16]; 
  BigAEM(34,16) = 0.2981423969999719*m1r[22]+0.4391550328268399*m1r[15]; 
  BigAEM(34,17) = 0.4472135954999579*m1r[22]; 
  BigAEM(34,18) = 0.2981423969999719*m1r[20]+0.3927922024247863*m1r[19]+0.4391550328268399*m1r[13]; 
  BigAEM(34,19) = 0.3927922024247863*m1r[18]; 
  BigAEM(34,20) = 0.2981423969999719*m1r[18]+0.5*m1r[14]; 
  BigAEM(34,22) = 0.4472135954999579*m1r[17]+0.2981423969999719*m1r[16]+0.5*m1r[12]; 
  BigAEM(35,12) = 0.5*m1r[23]; 
  BigAEM(35,13) = 0.5*m1r[21]; 
  BigAEM(35,14) = 0.4391550328268399*m1r[19]; 
  BigAEM(35,15) = 0.4391550328268399*m1r[17]; 
  BigAEM(35,16) = 0.4472135954999579*m1r[23]; 
  BigAEM(35,17) = 0.2981423969999719*m1r[23]+0.4391550328268399*m1r[15]; 
  BigAEM(35,18) = 0.3927922024247863*m1r[19]; 
  BigAEM(35,19) = 0.2981423969999719*m1r[21]+0.3927922024247863*m1r[18]+0.4391550328268399*m1r[14]; 
  BigAEM(35,21) = 0.2981423969999719*m1r[19]+0.5*m1r[13]; 
  BigAEM(35,23) = 0.2981423969999719*m1r[17]+0.4472135954999579*m1r[16]+0.5*m1r[12]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(24,24) = m0r[0]-0.5*cE[0]; 
  BigAEM(24,25) = m0r[1]-0.5*cE[1]; 
  BigAEM(24,26) = m0r[2]-0.5*cE[2]; 
  BigAEM(24,27) = m0r[3]-0.5*cE[3]; 
  BigAEM(24,28) = m0r[4]-0.5*cE[4]; 
  BigAEM(24,29) = m0r[5]-0.5*cE[5]; 
  BigAEM(24,30) = m0r[6]-0.5*cE[6]; 
  BigAEM(24,31) = m0r[7]-0.5*cE[7]; 
  BigAEM(24,32) = m0r[8]-0.5*cE[8]; 
  BigAEM(24,33) = m0r[9]-0.5*cE[9]; 
  BigAEM(24,34) = m0r[10]-0.5*cE[10]; 
  BigAEM(24,35) = m0r[11]-0.5*cE[11]; 
  BigAEM(25,24) = m0r[1]-0.5*cE[1]; 
  BigAEM(25,25) = 0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(25,26) = m0r[3]-0.5*cE[3]; 
  BigAEM(25,27) = 0.8944271909999161*m0r[6]-0.447213595499958*cE[6]+m0r[2]-0.5*cE[2]; 
  BigAEM(25,28) = 0.8783100656536796*m0r[8]-0.4391550328268398*cE[8]+0.8944271909999159*m0r[1]-0.4472135954999579*cE[1]; 
  BigAEM(25,29) = 1.0*m0r[7]-0.5000000000000001*cE[7]; 
  BigAEM(25,30) = 0.8783100656536798*m0r[10]-0.4391550328268399*cE[10]+0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(25,31) = 1.0*m0r[5]-0.5000000000000001*cE[5]; 
  BigAEM(25,32) = 0.8783100656536796*m0r[4]-0.4391550328268398*cE[4]; 
  BigAEM(25,33) = 1.0*m0r[11]-0.5*cE[11]; 
  BigAEM(25,34) = 0.8783100656536798*m0r[6]-0.4391550328268399*cE[6]; 
  BigAEM(25,35) = 1.0*m0r[9]-0.5*cE[9]; 
  BigAEM(26,24) = m0r[2]-0.5*cE[2]; 
  BigAEM(26,25) = m0r[3]-0.5*cE[3]; 
  BigAEM(26,26) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+m0r[0]-0.5*cE[0]; 
  BigAEM(26,27) = 0.8944271909999161*m0r[7]-0.447213595499958*cE[7]+m0r[1]-0.5*cE[1]; 
  BigAEM(26,28) = 1.0*m0r[6]-0.5000000000000001*cE[6]; 
  BigAEM(26,29) = 0.8783100656536796*m0r[9]-0.4391550328268398*cE[9]+0.8944271909999159*m0r[2]-0.4472135954999579*cE[2]; 
  BigAEM(26,30) = 1.0*m0r[4]-0.5000000000000001*cE[4]; 
  BigAEM(26,31) = 0.8783100656536798*m0r[11]-0.4391550328268399*cE[11]+0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(26,32) = 1.0*m0r[10]-0.5*cE[10]; 
  BigAEM(26,33) = 0.8783100656536796*m0r[5]-0.4391550328268398*cE[5]; 
  BigAEM(26,34) = 1.0*m0r[8]-0.5*cE[8]; 
  BigAEM(26,35) = 0.8783100656536798*m0r[7]-0.4391550328268399*cE[7]; 
  BigAEM(27,24) = m0r[3]-0.5*cE[3]; 
  BigAEM(27,25) = 0.8944271909999161*m0r[6]-0.447213595499958*cE[6]+m0r[2]-0.5*cE[2]; 
  BigAEM(27,26) = 0.8944271909999161*m0r[7]-0.447213595499958*cE[7]+m0r[1]-0.5*cE[1]; 
  BigAEM(27,27) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(27,28) = 0.8783100656536798*m0r[10]-0.4391550328268399*cE[10]+0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(27,29) = 0.8783100656536798*m0r[11]-0.4391550328268399*cE[11]+0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(27,30) = 0.8783100656536798*m0r[8]-0.4391550328268399*cE[8]+0.8*m0r[7]-0.4*cE[7]+0.8944271909999161*m0r[1]-0.447213595499958*cE[1]; 
  BigAEM(27,31) = 0.8783100656536798*m0r[9]-0.4391550328268399*cE[9]+0.8*m0r[6]-0.4*cE[6]+0.8944271909999161*m0r[2]-0.447213595499958*cE[2]; 
  BigAEM(27,32) = 0.8783100656536798*m0r[6]-0.4391550328268399*cE[6]; 
  BigAEM(27,33) = 0.8783100656536798*m0r[7]-0.4391550328268399*cE[7]; 
  BigAEM(27,34) = 0.8783100656536798*m0r[4]-0.4391550328268399*cE[4]; 
  BigAEM(27,35) = 0.8783100656536798*m0r[5]-0.4391550328268399*cE[5]; 
  BigAEM(28,24) = m0r[4]-0.5*cE[4]; 
  BigAEM(28,25) = 0.8783100656536796*m0r[8]-0.4391550328268398*cE[8]+0.8944271909999159*m0r[1]-0.4472135954999579*cE[1]; 
  BigAEM(28,26) = 1.0*m0r[6]-0.5000000000000001*cE[6]; 
  BigAEM(28,27) = 0.8783100656536798*m0r[10]-0.4391550328268399*cE[10]+0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(28,28) = 0.6388765649999399*m0r[4]-0.31943828249997*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(28,30) = 0.6388765649999399*m0r[6]-0.31943828249997*cE[6]+1.0*m0r[2]-0.5000000000000001*cE[2]; 
  BigAEM(28,31) = 0.8944271909999159*m0r[7]-0.4472135954999579*cE[7]; 
  BigAEM(28,32) = 0.5962847939999438*m0r[8]-0.2981423969999719*cE[8]+0.8783100656536796*m0r[1]-0.4391550328268398*cE[1]; 
  BigAEM(28,34) = 0.5962847939999438*m0r[10]-0.2981423969999719*cE[10]+0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(28,35) = 0.8944271909999159*m0r[11]-0.4472135954999579*cE[11]; 
  BigAEM(29,24) = m0r[5]-0.5*cE[5]; 
  BigAEM(29,25) = 1.0*m0r[7]-0.5000000000000001*cE[7]; 
  BigAEM(29,26) = 0.8783100656536796*m0r[9]-0.4391550328268398*cE[9]+0.8944271909999159*m0r[2]-0.4472135954999579*cE[2]; 
  BigAEM(29,27) = 0.8783100656536798*m0r[11]-0.4391550328268399*cE[11]+0.8944271909999159*m0r[3]-0.4472135954999579*cE[3]; 
  BigAEM(29,29) = 0.6388765649999399*m0r[5]-0.31943828249997*cE[5]+m0r[0]-0.5*cE[0]; 
  BigAEM(29,30) = 0.8944271909999159*m0r[6]-0.4472135954999579*cE[6]; 
  BigAEM(29,31) = 0.6388765649999399*m0r[7]-0.31943828249997*cE[7]+1.0*m0r[1]-0.5000000000000001*cE[1]; 
  BigAEM(29,33) = 0.5962847939999438*m0r[9]-0.2981423969999719*cE[9]+0.8783100656536796*m0r[2]-0.4391550328268398*cE[2]; 
  BigAEM(29,34) = 0.8944271909999159*m0r[10]-0.4472135954999579*cE[10]; 
  BigAEM(29,35) = 0.5962847939999438*m0r[11]-0.2981423969999719*cE[11]+0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(30,24) = m0r[6]-0.5*cE[6]; 
  BigAEM(30,25) = 0.8783100656536798*m0r[10]-0.4391550328268399*cE[10]+0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(30,26) = 1.0*m0r[4]-0.5000000000000001*cE[4]; 
  BigAEM(30,27) = 0.8783100656536798*m0r[8]-0.4391550328268399*cE[8]+0.8*m0r[7]-0.4*cE[7]+0.8944271909999161*m0r[1]-0.447213595499958*cE[1]; 
  BigAEM(30,28) = 0.6388765649999399*m0r[6]-0.31943828249997*cE[6]+1.0*m0r[2]-0.5000000000000001*cE[2]; 
  BigAEM(30,29) = 0.8944271909999159*m0r[6]-0.4472135954999579*cE[6]; 
  BigAEM(30,30) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+0.6388765649999399*m0r[4]-0.31943828249997*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(30,31) = 0.7855844048495726*m0r[11]-0.3927922024247863*cE[11]+0.7855844048495726*m0r[10]-0.3927922024247863*cE[10]+0.8*m0r[3]-0.4*cE[3]; 
  BigAEM(30,32) = 0.5962847939999437*m0r[10]-0.2981423969999719*cE[10]+0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(30,34) = 0.5962847939999437*m0r[8]-0.2981423969999719*cE[8]+0.7855844048495726*m0r[7]-0.3927922024247863*cE[7]+0.8783100656536798*m0r[1]-0.4391550328268399*cE[1]; 
  BigAEM(30,35) = 0.7855844048495726*m0r[7]-0.3927922024247863*cE[7]; 
  BigAEM(31,24) = m0r[7]-0.5*cE[7]; 
  BigAEM(31,25) = 1.0*m0r[5]-0.5000000000000001*cE[5]; 
  BigAEM(31,26) = 0.8783100656536798*m0r[11]-0.4391550328268399*cE[11]+0.8944271909999161*m0r[3]-0.447213595499958*cE[3]; 
  BigAEM(31,27) = 0.8783100656536798*m0r[9]-0.4391550328268399*cE[9]+0.8*m0r[6]-0.4*cE[6]+0.8944271909999161*m0r[2]-0.447213595499958*cE[2]; 
  BigAEM(31,28) = 0.8944271909999159*m0r[7]-0.4472135954999579*cE[7]; 
  BigAEM(31,29) = 0.6388765649999399*m0r[7]-0.31943828249997*cE[7]+1.0*m0r[1]-0.5000000000000001*cE[1]; 
  BigAEM(31,30) = 0.7855844048495726*m0r[11]-0.3927922024247863*cE[11]+0.7855844048495726*m0r[10]-0.3927922024247863*cE[10]+0.8*m0r[3]-0.4*cE[3]; 
  BigAEM(31,31) = 0.6388765649999399*m0r[5]-0.31943828249997*cE[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(31,33) = 0.5962847939999437*m0r[11]-0.2981423969999719*cE[11]+0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(31,34) = 0.7855844048495726*m0r[6]-0.3927922024247863*cE[6]; 
  BigAEM(31,35) = 0.5962847939999437*m0r[9]-0.2981423969999719*cE[9]+0.7855844048495726*m0r[6]-0.3927922024247863*cE[6]+0.8783100656536798*m0r[2]-0.4391550328268399*cE[2]; 
  BigAEM(32,24) = m0r[8]-0.5*cE[8]; 
  BigAEM(32,25) = 0.8783100656536796*m0r[4]-0.4391550328268398*cE[4]; 
  BigAEM(32,26) = 1.0*m0r[10]-0.5*cE[10]; 
  BigAEM(32,27) = 0.8783100656536798*m0r[6]-0.4391550328268399*cE[6]; 
  BigAEM(32,28) = 0.5962847939999438*m0r[8]-0.2981423969999719*cE[8]+0.8783100656536796*m0r[1]-0.4391550328268398*cE[1]; 
  BigAEM(32,30) = 0.5962847939999437*m0r[10]-0.2981423969999719*cE[10]+0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(32,32) = 0.5962847939999438*m0r[4]-0.2981423969999719*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(32,34) = 0.5962847939999437*m0r[6]-0.2981423969999719*cE[6]+1.0*m0r[2]-0.5*cE[2]; 
  BigAEM(33,24) = m0r[9]-0.5*cE[9]; 
  BigAEM(33,25) = 1.0*m0r[11]-0.5*cE[11]; 
  BigAEM(33,26) = 0.8783100656536796*m0r[5]-0.4391550328268398*cE[5]; 
  BigAEM(33,27) = 0.8783100656536798*m0r[7]-0.4391550328268399*cE[7]; 
  BigAEM(33,29) = 0.5962847939999438*m0r[9]-0.2981423969999719*cE[9]+0.8783100656536796*m0r[2]-0.4391550328268398*cE[2]; 
  BigAEM(33,31) = 0.5962847939999437*m0r[11]-0.2981423969999719*cE[11]+0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(33,33) = 0.5962847939999438*m0r[5]-0.2981423969999719*cE[5]+m0r[0]-0.5*cE[0]; 
  BigAEM(33,35) = 0.5962847939999437*m0r[7]-0.2981423969999719*cE[7]+1.0*m0r[1]-0.5*cE[1]; 
  BigAEM(34,24) = m0r[10]-0.5*cE[10]; 
  BigAEM(34,25) = 0.8783100656536798*m0r[6]-0.4391550328268399*cE[6]; 
  BigAEM(34,26) = 1.0*m0r[8]-0.5*cE[8]; 
  BigAEM(34,27) = 0.8783100656536798*m0r[4]-0.4391550328268399*cE[4]; 
  BigAEM(34,28) = 0.5962847939999438*m0r[10]-0.2981423969999719*cE[10]+0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(34,29) = 0.8944271909999159*m0r[10]-0.4472135954999579*cE[10]; 
  BigAEM(34,30) = 0.5962847939999437*m0r[8]-0.2981423969999719*cE[8]+0.7855844048495726*m0r[7]-0.3927922024247863*cE[7]+0.8783100656536798*m0r[1]-0.4391550328268399*cE[1]; 
  BigAEM(34,31) = 0.7855844048495726*m0r[6]-0.3927922024247863*cE[6]; 
  BigAEM(34,32) = 0.5962847939999437*m0r[6]-0.2981423969999719*cE[6]+1.0*m0r[2]-0.5*cE[2]; 
  BigAEM(34,34) = 0.8944271909999159*m0r[5]-0.4472135954999579*cE[5]+0.5962847939999438*m0r[4]-0.2981423969999719*cE[4]+m0r[0]-0.5*cE[0]; 
  BigAEM(35,24) = m0r[11]-0.5*cE[11]; 
  BigAEM(35,25) = 1.0*m0r[9]-0.5*cE[9]; 
  BigAEM(35,26) = 0.8783100656536798*m0r[7]-0.4391550328268399*cE[7]; 
  BigAEM(35,27) = 0.8783100656536798*m0r[5]-0.4391550328268399*cE[5]; 
  BigAEM(35,28) = 0.8944271909999159*m0r[11]-0.4472135954999579*cE[11]; 
  BigAEM(35,29) = 0.5962847939999438*m0r[11]-0.2981423969999719*cE[11]+0.8783100656536798*m0r[3]-0.4391550328268399*cE[3]; 
  BigAEM(35,30) = 0.7855844048495726*m0r[7]-0.3927922024247863*cE[7]; 
  BigAEM(35,31) = 0.5962847939999437*m0r[9]-0.2981423969999719*cE[9]+0.7855844048495726*m0r[6]-0.3927922024247863*cE[6]+0.8783100656536798*m0r[2]-0.4391550328268399*cE[2]; 
  BigAEM(35,33) = 0.5962847939999437*m0r[7]-0.2981423969999719*cE[7]+1.0*m0r[1]-0.5*cE[1]; 
  BigAEM(35,35) = 0.5962847939999438*m0r[5]-0.2981423969999719*cE[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cE[4]+m0r[0]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<12,12>(0,12).setZero(); 
  BigAEM.block<12,12>(12,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m1r[12],m1r[13],m1r[14],m1r[15],m1r[16],m1r[17],m1r[18],m1r[19],m1r[20],m1r[21],m1r[22],m1r[23],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5],m2r[6],m2r[7],m2r[8],m2r[9],m2r[10],m2r[11]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,24,1) = xEV.segment<24>(0); 
 
  Eigen::Map<VectorXd>(vtSq,12,1) = xEV.segment<12>(24); 
 
} 
 
