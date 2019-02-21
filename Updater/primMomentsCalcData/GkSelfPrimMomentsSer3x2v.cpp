
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void GkM0Star3x2vSer_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[5]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[4]*intFac*(wr[3]-wl[3]); 
 
  out[0] += ((-0.5773502691896258*fr[4])+0.5773502691896258*fl[4]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += ((-0.5773502691896258*fr[9])+0.5773502691896258*fl[9]+0.5*fr[1]+0.5*fl[1])*dS; 
  out[2] += ((-0.5773502691896258*fr[10])+0.5773502691896258*fl[10]+0.5*fr[2]+0.5*fl[2])*dS; 
  out[3] += ((-0.5773502691896258*fr[11])+0.5773502691896258*fl[11]+0.5*fr[3]+0.5*fl[3])*dS; 
  out[4] += ((-0.5773502691896258*fr[17])+0.5773502691896258*fl[17]+0.5*fr[6]+0.5*fl[6])*dS; 
  out[5] += ((-0.5773502691896258*fr[18])+0.5773502691896258*fl[18]+0.5*fr[7]+0.5*fl[7])*dS; 
  out[6] += ((-0.5773502691896258*fr[19])+0.5773502691896258*fl[19]+0.5*fr[8]+0.5*fl[8])*dS; 
  out[7] += ((-0.5773502691896258*fr[26])+0.5773502691896258*fl[26]+0.5*fr[16]+0.5*fl[16])*dS; 
 
} 
 
void GkM1iM2Star3x2vSer(const double *w, const double *dxv, const double intFac, const double m_, const double *Bmag, const double *f, double *outM1i, double *outM2) 
{ 
  // w[5]:    Cell-center coordinates. 
  // dxv[5]:  Cell length in each direciton. 
  // intFac:  =2pi/m for gyrokinetics. 
  // m_:      mass. 
  // Bmag[8]: Magnetic field magnitude. 
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
  outM1i[4] += 2.0*w[3]*f[6]*volFact; 
  outM1i[5] += 2.0*w[3]*f[7]*volFact; 
  outM1i[6] += 2.0*w[3]*f[8]*volFact; 
  outM1i[7] += 2.0*w[3]*f[16]*volFact; 
 
  double tmp[8]; 
  tmp[0] = 0.5773502691896258*dxv[4]*f[5]+2.0*f[0]*w[4]; 
  tmp[1] = 0.5773502691896258*dxv[4]*f[12]+2.0*f[1]*w[4]; 
  tmp[2] = 0.5773502691896258*dxv[4]*f[13]+2.0*f[2]*w[4]; 
  tmp[3] = 0.5773502691896258*dxv[4]*f[14]+2.0*f[3]*w[4]; 
  tmp[4] = 0.5773502691896258*dxv[4]*f[20]+2.0*w[4]*f[6]; 
  tmp[5] = 0.5773502691896258*dxv[4]*f[21]+2.0*w[4]*f[7]; 
  tmp[6] = 0.5773502691896258*dxv[4]*f[22]+2.0*w[4]*f[8]; 
  tmp[7] = 0.5773502691896258*dxv[4]*f[27]+2.0*w[4]*f[16]; 
 
  outM2[0] += ((0.7071067811865474*Bmag[7]*tmp[7]+0.7071067811865474*Bmag[6]*tmp[6]+0.7071067811865474*Bmag[5]*tmp[5]+0.7071067811865474*Bmag[4]*tmp[4]+0.7071067811865474*Bmag[3]*tmp[3]+0.7071067811865474*Bmag[2]*tmp[2]+0.7071067811865474*Bmag[1]*tmp[1]+0.7071067811865474*Bmag[0]*tmp[0])/m_+0.5773502691896258*dxv[3]*w[3]*f[4]+2.0*f[0]*wvSq[0])*volFact; 
  outM2[1] += ((0.7071067811865474*Bmag[6]*tmp[7]+0.7071067811865474*tmp[6]*Bmag[7]+0.7071067811865474*Bmag[3]*tmp[5]+0.7071067811865474*tmp[3]*Bmag[5]+0.7071067811865474*Bmag[2]*tmp[4]+0.7071067811865474*tmp[2]*Bmag[4]+0.7071067811865474*Bmag[0]*tmp[1]+0.7071067811865474*tmp[0]*Bmag[1])/m_+0.5773502691896258*dxv[3]*w[3]*f[9]+2.0*f[1]*wvSq[0])*volFact; 
  outM2[2] += ((0.7071067811865474*Bmag[5]*tmp[7]+0.7071067811865474*tmp[5]*Bmag[7]+0.7071067811865474*Bmag[3]*tmp[6]+0.7071067811865474*tmp[3]*Bmag[6]+0.7071067811865474*Bmag[1]*tmp[4]+0.7071067811865474*tmp[1]*Bmag[4]+0.7071067811865474*Bmag[0]*tmp[2]+0.7071067811865474*tmp[0]*Bmag[2])/m_+0.5773502691896258*dxv[3]*w[3]*f[10]+2.0*f[2]*wvSq[0])*volFact; 
  outM2[3] += ((0.7071067811865474*Bmag[4]*tmp[7]+0.7071067811865474*tmp[4]*Bmag[7]+0.7071067811865474*Bmag[2]*tmp[6]+0.7071067811865474*tmp[2]*Bmag[6]+0.7071067811865474*Bmag[1]*tmp[5]+0.7071067811865474*tmp[1]*Bmag[5]+0.7071067811865474*Bmag[0]*tmp[3]+0.7071067811865474*tmp[0]*Bmag[3])/m_+0.5773502691896258*dxv[3]*w[3]*f[11]+2.0*f[3]*wvSq[0])*volFact; 
  outM2[4] += ((0.7071067811865474*Bmag[3]*tmp[7]+0.7071067811865474*tmp[3]*Bmag[7]+0.7071067811865474*Bmag[5]*tmp[6]+0.7071067811865474*tmp[5]*Bmag[6]+0.7071067811865474*Bmag[0]*tmp[4]+0.7071067811865474*tmp[0]*Bmag[4]+0.7071067811865474*Bmag[1]*tmp[2]+0.7071067811865474*tmp[1]*Bmag[2])/m_+0.5773502691896258*dxv[3]*w[3]*f[17]+2.0*wvSq[0]*f[6])*volFact; 
  outM2[5] += ((0.7071067811865474*Bmag[2]*tmp[7]+0.7071067811865474*tmp[2]*Bmag[7]+0.7071067811865474*Bmag[4]*tmp[6]+0.7071067811865474*tmp[4]*Bmag[6]+0.7071067811865474*Bmag[0]*tmp[5]+0.7071067811865474*tmp[0]*Bmag[5]+0.7071067811865474*Bmag[1]*tmp[3]+0.7071067811865474*tmp[1]*Bmag[3])/m_+0.5773502691896258*dxv[3]*w[3]*f[18]+2.0*wvSq[0]*f[7])*volFact; 
  outM2[6] += ((0.7071067811865474*Bmag[1]*tmp[7]+0.7071067811865474*tmp[1]*Bmag[7]+0.7071067811865474*Bmag[0]*tmp[6]+0.7071067811865474*tmp[0]*Bmag[6]+0.7071067811865474*Bmag[4]*tmp[5]+0.7071067811865474*tmp[4]*Bmag[5]+0.7071067811865474*Bmag[2]*tmp[3]+0.7071067811865474*tmp[2]*Bmag[3])/m_+0.5773502691896258*dxv[3]*w[3]*f[19]+2.0*wvSq[0]*f[8])*volFact; 
  outM2[7] += ((0.7071067811865474*Bmag[0]*tmp[7]+0.7071067811865474*tmp[0]*Bmag[7]+0.7071067811865474*Bmag[1]*tmp[6]+0.7071067811865474*tmp[1]*Bmag[6]+0.7071067811865474*Bmag[2]*tmp[5]+0.7071067811865474*tmp[2]*Bmag[5]+0.7071067811865474*Bmag[3]*tmp[4]+0.7071067811865474*tmp[3]*Bmag[4])/m_+0.5773502691896258*dxv[3]*w[3]*f[26]+2.0*wvSq[0]*f[16])*volFact; 
 
} 
void GkBoundaryIntegral3x2vSer_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[9]-1.0*fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[10]-1.0*fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]-1.0*fIn[3])*dS; 
    out[4] += (1.732050807568877*fIn[17]-1.0*fIn[6])*dS; 
    out[5] += (1.732050807568877*fIn[18]-1.0*fIn[7])*dS; 
    out[6] += (1.732050807568877*fIn[19]-1.0*fIn[8])*dS; 
    out[7] += (1.732050807568877*fIn[26]-1.0*fIn[16])*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[1] += (1.732050807568877*fIn[9]+fIn[1])*dS; 
    out[2] += (1.732050807568877*fIn[10]+fIn[2])*dS; 
    out[3] += (1.732050807568877*fIn[11]+fIn[3])*dS; 
    out[4] += (1.732050807568877*fIn[17]+fIn[6])*dS; 
    out[5] += (1.732050807568877*fIn[18]+fIn[7])*dS; 
    out[6] += (1.732050807568877*fIn[19]+fIn[8])*dS; 
    out[7] += (1.732050807568877*fIn[26]+fIn[16])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral3x2vSer_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[112]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[19])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS; 
    out[1] += ((-2.23606797749979*fIn[40])+1.732050807568877*fIn[9]-1.0*fIn[1])*dS; 
    out[2] += ((-2.23606797749979*fIn[41])+1.732050807568877*fIn[10]-1.0*fIn[2])*dS; 
    out[3] += ((-2.23606797749979*fIn[42])+1.732050807568877*fIn[11]-1.0*fIn[3])*dS; 
    out[4] += ((-2.23606797749979*fIn[65])+1.732050807568877*fIn[22]-1.0*fIn[6])*dS; 
    out[5] += ((-2.23606797749979*fIn[66])+1.732050807568877*fIn[23]-1.0*fIn[7])*dS; 
    out[6] += ((-2.23606797749979*fIn[67])+1.732050807568877*fIn[24]-1.0*fIn[8])*dS; 
    out[7] += (1.732050807568877*fIn[37]-1.0*fIn[16])*dS; 
    out[8] += (1.732050807568877*fIn[38]-1.0*fIn[17])*dS; 
    out[9] += (1.732050807568877*fIn[39]-1.0*fIn[18])*dS; 
    out[10] += ((-2.23606797749979*fIn[90])+1.732050807568877*fIn[51]-1.0*fIn[21])*dS; 
    out[11] += (1.732050807568877*fIn[59]-1.0*fIn[31])*dS; 
    out[12] += (1.732050807568877*fIn[60]-1.0*fIn[32])*dS; 
    out[13] += (1.732050807568877*fIn[61]-1.0*fIn[33])*dS; 
    out[14] += (1.732050807568877*fIn[62]-1.0*fIn[34])*dS; 
    out[15] += (1.732050807568877*fIn[63]-1.0*fIn[35])*dS; 
    out[16] += (1.732050807568877*fIn[64]-1.0*fIn[36])*dS; 
    out[17] += (1.732050807568877*fIn[87]-1.0*fIn[56])*dS; 
    out[18] += (1.732050807568877*fIn[88]-1.0*fIn[57])*dS; 
    out[19] += (1.732050807568877*fIn[89]-1.0*fIn[58])*dS; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[19]+1.732050807568877*fIn[4]+fIn[0])*dS; 
    out[1] += (2.23606797749979*fIn[40]+1.732050807568877*fIn[9]+fIn[1])*dS; 
    out[2] += (2.23606797749979*fIn[41]+1.732050807568877*fIn[10]+fIn[2])*dS; 
    out[3] += (2.23606797749979*fIn[42]+1.732050807568877*fIn[11]+fIn[3])*dS; 
    out[4] += (2.23606797749979*fIn[65]+1.732050807568877*fIn[22]+fIn[6])*dS; 
    out[5] += (2.23606797749979*fIn[66]+1.732050807568877*fIn[23]+fIn[7])*dS; 
    out[6] += (2.23606797749979*fIn[67]+1.732050807568877*fIn[24]+fIn[8])*dS; 
    out[7] += (1.732050807568877*fIn[37]+fIn[16])*dS; 
    out[8] += (1.732050807568877*fIn[38]+fIn[17])*dS; 
    out[9] += (1.732050807568877*fIn[39]+fIn[18])*dS; 
    out[10] += (2.23606797749979*fIn[90]+1.732050807568877*fIn[51]+fIn[21])*dS; 
    out[11] += (1.732050807568877*fIn[59]+fIn[31])*dS; 
    out[12] += (1.732050807568877*fIn[60]+fIn[32])*dS; 
    out[13] += (1.732050807568877*fIn[61]+fIn[33])*dS; 
    out[14] += (1.732050807568877*fIn[62]+fIn[34])*dS; 
    out[15] += (1.732050807568877*fIn[63]+fIn[35])*dS; 
    out[16] += (1.732050807568877*fIn[64]+fIn[36])*dS; 
    out[17] += (1.732050807568877*fIn[87]+fIn[56])*dS; 
    out[18] += (1.732050807568877*fIn[88]+fIn[57])*dS; 
    out[19] += (1.732050807568877*fIn[89]+fIn[58])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral3x2vSer_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[4]-0.5*fIn[0]*dxv[3])*dS; 
    out[1] += (1.732050807568877*fIn[9]-1.0*fIn[1])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[9]-0.5*fIn[1]*dxv[3])*dS; 
    out[2] += (1.732050807568877*fIn[10]-1.0*fIn[2])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[10]-0.5*fIn[2]*dxv[3])*dS; 
    out[3] += (1.732050807568877*fIn[11]-1.0*fIn[3])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[11]-0.5*dxv[3]*fIn[3])*dS; 
    out[4] += (1.732050807568877*fIn[17]-1.0*fIn[6])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[17]-0.5*dxv[3]*fIn[6])*dS; 
    out[5] += (1.732050807568877*fIn[18]-1.0*fIn[7])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[18]-0.5*dxv[3]*fIn[7])*dS; 
    out[6] += (1.732050807568877*fIn[19]-1.0*fIn[8])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[19]-0.5*dxv[3]*fIn[8])*dS; 
    out[7] += (1.732050807568877*fIn[26]-1.0*fIn[16])*dS*vBoundary+(0.8660254037844386*dxv[3]*fIn[26]-0.5*dxv[3]*fIn[16])*dS; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[4])-0.5*fIn[0]*dxv[3])*dS; 
    out[1] += (1.732050807568877*fIn[9]+fIn[1])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[9])-0.5*fIn[1]*dxv[3])*dS; 
    out[2] += (1.732050807568877*fIn[10]+fIn[2])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[10])-0.5*fIn[2]*dxv[3])*dS; 
    out[3] += (1.732050807568877*fIn[11]+fIn[3])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[11])-0.5*dxv[3]*fIn[3])*dS; 
    out[4] += (1.732050807568877*fIn[17]+fIn[6])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[17])-0.5*dxv[3]*fIn[6])*dS; 
    out[5] += (1.732050807568877*fIn[18]+fIn[7])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[18])-0.5*dxv[3]*fIn[7])*dS; 
    out[6] += (1.732050807568877*fIn[19]+fIn[8])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[19])-0.5*dxv[3]*fIn[8])*dS; 
    out[7] += (1.732050807568877*fIn[26]+fIn[16])*dS*vBoundary+((-0.8660254037844386*dxv[3]*fIn[26])-0.5*dxv[3]*fIn[16])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral3x2vSer_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[112]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[19])+1.732050807568877*fIn[4]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[40])+1.732050807568877*fIn[9]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += ((-2.23606797749979*fIn[41])+1.732050807568877*fIn[10]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += ((-2.23606797749979*fIn[42])+1.732050807568877*fIn[11]-1.0*fIn[3])*dS*vBoundary; 
    out[4] += ((-2.23606797749979*fIn[65])+1.732050807568877*fIn[22]-1.0*fIn[6])*dS*vBoundary; 
    out[5] += ((-2.23606797749979*fIn[66])+1.732050807568877*fIn[23]-1.0*fIn[7])*dS*vBoundary; 
    out[6] += ((-2.23606797749979*fIn[67])+1.732050807568877*fIn[24]-1.0*fIn[8])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[37]-1.0*fIn[16])*dS*vBoundary; 
    out[8] += (1.732050807568877*fIn[38]-1.0*fIn[17])*dS*vBoundary; 
    out[9] += (1.732050807568877*fIn[39]-1.0*fIn[18])*dS*vBoundary; 
    out[10] += ((-2.23606797749979*fIn[90])+1.732050807568877*fIn[51]-1.0*fIn[21])*dS*vBoundary; 
    out[11] += (1.732050807568877*fIn[59]-1.0*fIn[31])*dS*vBoundary; 
    out[12] += (1.732050807568877*fIn[60]-1.0*fIn[32])*dS*vBoundary; 
    out[13] += (1.732050807568877*fIn[61]-1.0*fIn[33])*dS*vBoundary; 
    out[14] += (1.732050807568877*fIn[62]-1.0*fIn[34])*dS*vBoundary; 
    out[15] += (1.732050807568877*fIn[63]-1.0*fIn[35])*dS*vBoundary; 
    out[16] += (1.732050807568877*fIn[64]-1.0*fIn[36])*dS*vBoundary; 
    out[17] += (1.732050807568877*fIn[87]-1.0*fIn[56])*dS*vBoundary; 
    out[18] += (1.732050807568877*fIn[88]-1.0*fIn[57])*dS*vBoundary; 
    out[19] += (1.732050807568877*fIn[89]-1.0*fIn[58])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[19]+1.732050807568877*fIn[4]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[40]+1.732050807568877*fIn[9]+fIn[1])*dS*vBoundary; 
    out[2] += (2.23606797749979*fIn[41]+1.732050807568877*fIn[10]+fIn[2])*dS*vBoundary; 
    out[3] += (2.23606797749979*fIn[42]+1.732050807568877*fIn[11]+fIn[3])*dS*vBoundary; 
    out[4] += (2.23606797749979*fIn[65]+1.732050807568877*fIn[22]+fIn[6])*dS*vBoundary; 
    out[5] += (2.23606797749979*fIn[66]+1.732050807568877*fIn[23]+fIn[7])*dS*vBoundary; 
    out[6] += (2.23606797749979*fIn[67]+1.732050807568877*fIn[24]+fIn[8])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[37]+fIn[16])*dS*vBoundary; 
    out[8] += (1.732050807568877*fIn[38]+fIn[17])*dS*vBoundary; 
    out[9] += (1.732050807568877*fIn[39]+fIn[18])*dS*vBoundary; 
    out[10] += (2.23606797749979*fIn[90]+1.732050807568877*fIn[51]+fIn[21])*dS*vBoundary; 
    out[11] += (1.732050807568877*fIn[59]+fIn[31])*dS*vBoundary; 
    out[12] += (1.732050807568877*fIn[60]+fIn[32])*dS*vBoundary; 
    out[13] += (1.732050807568877*fIn[61]+fIn[33])*dS*vBoundary; 
    out[14] += (1.732050807568877*fIn[62]+fIn[34])*dS*vBoundary; 
    out[15] += (1.732050807568877*fIn[63]+fIn[35])*dS*vBoundary; 
    out[16] += (1.732050807568877*fIn[64]+fIn[36])*dS*vBoundary; 
    out[17] += (1.732050807568877*fIn[87]+fIn[56])*dS*vBoundary; 
    out[18] += (1.732050807568877*fIn[88]+fIn[57])*dS*vBoundary; 
    out[19] += (1.732050807568877*fIn[89]+fIn[58])*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral3x2vSer_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[32]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]*intFac; 
 
  if (atLower) {
 
    out[0] += (1.732050807568877*fIn[5]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[12]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[14]-1.0*fIn[3])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[20]-1.0*fIn[6])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[21]-1.0*fIn[7])*dS*vBoundary; 
    out[6] += (1.732050807568877*fIn[22]-1.0*fIn[8])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[27]-1.0*fIn[16])*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.732050807568877*fIn[5]+fIn[0])*dS*vBoundary; 
    out[1] += (1.732050807568877*fIn[12]+fIn[1])*dS*vBoundary; 
    out[2] += (1.732050807568877*fIn[13]+fIn[2])*dS*vBoundary; 
    out[3] += (1.732050807568877*fIn[14]+fIn[3])*dS*vBoundary; 
    out[4] += (1.732050807568877*fIn[20]+fIn[6])*dS*vBoundary; 
    out[5] += (1.732050807568877*fIn[21]+fIn[7])*dS*vBoundary; 
    out[6] += (1.732050807568877*fIn[22]+fIn[8])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[27]+fIn[16])*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral3x2vSer_vF_VY_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[5]:    cell length in each direciton. 
  // fIn[112]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]*intFac; 
 
  if (atLower) {
 
    out[0] += ((-2.23606797749979*fIn[20])+1.732050807568877*fIn[5]-1.0*fIn[0])*dS*vBoundary; 
    out[1] += ((-2.23606797749979*fIn[47])+1.732050807568877*fIn[12]-1.0*fIn[1])*dS*vBoundary; 
    out[2] += ((-2.23606797749979*fIn[48])+1.732050807568877*fIn[13]-1.0*fIn[2])*dS*vBoundary; 
    out[3] += ((-2.23606797749979*fIn[49])+1.732050807568877*fIn[14]-1.0*fIn[3])*dS*vBoundary; 
    out[4] += ((-2.23606797749979*fIn[80])+1.732050807568877*fIn[25]-1.0*fIn[6])*dS*vBoundary; 
    out[5] += ((-2.23606797749979*fIn[81])+1.732050807568877*fIn[26]-1.0*fIn[7])*dS*vBoundary; 
    out[6] += ((-2.23606797749979*fIn[82])+1.732050807568877*fIn[27]-1.0*fIn[8])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[43]-1.0*fIn[16])*dS*vBoundary; 
    out[8] += (1.732050807568877*fIn[44]-1.0*fIn[17])*dS*vBoundary; 
    out[9] += (1.732050807568877*fIn[45]-1.0*fIn[18])*dS*vBoundary; 
    out[10] += ((-2.23606797749979*fIn[103])+1.732050807568877*fIn[52]-1.0*fIn[21])*dS*vBoundary; 
    out[11] += (1.732050807568877*fIn[68]-1.0*fIn[31])*dS*vBoundary; 
    out[12] += (1.732050807568877*fIn[69]-1.0*fIn[32])*dS*vBoundary; 
    out[13] += (1.732050807568877*fIn[70]-1.0*fIn[33])*dS*vBoundary; 
    out[14] += (1.732050807568877*fIn[71]-1.0*fIn[34])*dS*vBoundary; 
    out[15] += (1.732050807568877*fIn[72]-1.0*fIn[35])*dS*vBoundary; 
    out[16] += (1.732050807568877*fIn[73]-1.0*fIn[36])*dS*vBoundary; 
    out[17] += (1.732050807568877*fIn[91]-1.0*fIn[56])*dS*vBoundary; 
    out[18] += (1.732050807568877*fIn[92]-1.0*fIn[57])*dS*vBoundary; 
    out[19] += (1.732050807568877*fIn[93]-1.0*fIn[58])*dS*vBoundary; 
 
  } else {
 
    out[0] += (2.23606797749979*fIn[20]+1.732050807568877*fIn[5]+fIn[0])*dS*vBoundary; 
    out[1] += (2.23606797749979*fIn[47]+1.732050807568877*fIn[12]+fIn[1])*dS*vBoundary; 
    out[2] += (2.23606797749979*fIn[48]+1.732050807568877*fIn[13]+fIn[2])*dS*vBoundary; 
    out[3] += (2.23606797749979*fIn[49]+1.732050807568877*fIn[14]+fIn[3])*dS*vBoundary; 
    out[4] += (2.23606797749979*fIn[80]+1.732050807568877*fIn[25]+fIn[6])*dS*vBoundary; 
    out[5] += (2.23606797749979*fIn[81]+1.732050807568877*fIn[26]+fIn[7])*dS*vBoundary; 
    out[6] += (2.23606797749979*fIn[82]+1.732050807568877*fIn[27]+fIn[8])*dS*vBoundary; 
    out[7] += (1.732050807568877*fIn[43]+fIn[16])*dS*vBoundary; 
    out[8] += (1.732050807568877*fIn[44]+fIn[17])*dS*vBoundary; 
    out[9] += (1.732050807568877*fIn[45]+fIn[18])*dS*vBoundary; 
    out[10] += (2.23606797749979*fIn[103]+1.732050807568877*fIn[52]+fIn[21])*dS*vBoundary; 
    out[11] += (1.732050807568877*fIn[68]+fIn[31])*dS*vBoundary; 
    out[12] += (1.732050807568877*fIn[69]+fIn[32])*dS*vBoundary; 
    out[13] += (1.732050807568877*fIn[70]+fIn[33])*dS*vBoundary; 
    out[14] += (1.732050807568877*fIn[71]+fIn[34])*dS*vBoundary; 
    out[15] += (1.732050807568877*fIn[72]+fIn[35])*dS*vBoundary; 
    out[16] += (1.732050807568877*fIn[73]+fIn[36])*dS*vBoundary; 
    out[17] += (1.732050807568877*fIn[91]+fIn[56])*dS*vBoundary; 
    out[18] += (1.732050807568877*fIn[92]+fIn[57])*dS*vBoundary; 
    out[19] += (1.732050807568877*fIn[93]+fIn[58])*dS*vBoundary; 
 
  }
 
} 
 
void GkSelfPrimMoments3x2vSer_P1(binOpData_t* data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1:       moments of the distribution function. 
  // m0S,m1S,m1S: star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-1.837117307087383*m0[7])+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.837117307087383*m0[7])+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.837117307087383*m0[7])+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.837117307087383*m0[7])+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.837117307087383*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.837117307087383*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.837117307087383*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.837117307087383*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[8]; 
  double m1r[8]; 
  double m0Sr[8]; 
  double m1Sr[8]; 
  double m2Sr[8]; 
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
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m0Sr[2] = 0.0; 
    m0Sr[3] = 0.0; 
    m0Sr[4] = 0.0; 
    m0Sr[5] = 0.0; 
    m0Sr[6] = 0.0; 
    m0Sr[7] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = 0.0; 
    m1Sr[3] = 0.0; 
    m1Sr[4] = 0.0; 
    m1Sr[5] = 0.0; 
    m1Sr[6] = 0.0; 
    m1Sr[7] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
    m2Sr[2] = 0.0; 
    m2Sr[3] = 0.0; 
    m2Sr[4] = 0.0; 
    m2Sr[5] = 0.0; 
    m2Sr[6] = 0.0; 
    m2Sr[7] = 0.0; 
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
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m0Sr[2] = m0S[2]; 
    m0Sr[3] = m0S[3]; 
    m0Sr[4] = m0S[4]; 
    m0Sr[5] = m0S[5]; 
    m0Sr[6] = m0S[6]; 
    m0Sr[7] = m0S[7]; 
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
    m2Sr[4] = m2S[4]; 
    m2Sr[5] = m2S[5]; 
    m2Sr[6] = m2S[6]; 
    m2Sr[7] = m2S[7]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(16,16); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(0,1) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(0,2) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(0,3) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(0,4) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(0,5) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(0,6) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(0,7) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(1,0) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(1,1) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(1,2) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(1,3) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(1,4) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(1,5) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(1,6) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(1,7) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(2,0) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(2,1) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(2,2) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(2,3) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(2,4) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(2,5) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(2,6) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(2,7) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(3,0) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(3,1) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(3,2) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(3,3) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(3,4) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(3,5) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(3,6) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(3,7) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(4,0) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(4,1) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(4,2) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(4,3) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(4,4) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(4,5) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(4,6) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(4,7) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(5,0) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(5,1) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(5,2) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(5,3) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(5,4) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(5,5) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(5,6) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(5,7) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(6,0) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(6,1) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(6,2) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(6,3) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(6,4) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(6,5) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(6,6) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(6,7) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(7,0) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(7,1) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(7,2) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(7,3) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(7,4) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(7,5) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(7,6) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(7,7) = 0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,8) = -0.3535533905932737*cM[0]; 
  data->AEM_S(0,9) = -0.3535533905932737*cM[1]; 
  data->AEM_S(0,10) = -0.3535533905932737*cM[2]; 
  data->AEM_S(0,11) = -0.3535533905932737*cM[3]; 
  data->AEM_S(0,12) = -0.3535533905932737*cM[4]; 
  data->AEM_S(0,13) = -0.3535533905932737*cM[5]; 
  data->AEM_S(0,14) = -0.3535533905932737*cM[6]; 
  data->AEM_S(0,15) = -0.3535533905932737*cM[7]; 
  data->AEM_S(1,8) = -0.3535533905932737*cM[1]; 
  data->AEM_S(1,9) = -0.3535533905932737*cM[0]; 
  data->AEM_S(1,10) = -0.3535533905932737*cM[4]; 
  data->AEM_S(1,11) = -0.3535533905932737*cM[5]; 
  data->AEM_S(1,12) = -0.3535533905932737*cM[2]; 
  data->AEM_S(1,13) = -0.3535533905932737*cM[3]; 
  data->AEM_S(1,14) = -0.3535533905932737*cM[7]; 
  data->AEM_S(1,15) = -0.3535533905932737*cM[6]; 
  data->AEM_S(2,8) = -0.3535533905932737*cM[2]; 
  data->AEM_S(2,9) = -0.3535533905932737*cM[4]; 
  data->AEM_S(2,10) = -0.3535533905932737*cM[0]; 
  data->AEM_S(2,11) = -0.3535533905932737*cM[6]; 
  data->AEM_S(2,12) = -0.3535533905932737*cM[1]; 
  data->AEM_S(2,13) = -0.3535533905932737*cM[7]; 
  data->AEM_S(2,14) = -0.3535533905932737*cM[3]; 
  data->AEM_S(2,15) = -0.3535533905932737*cM[5]; 
  data->AEM_S(3,8) = -0.3535533905932737*cM[3]; 
  data->AEM_S(3,9) = -0.3535533905932737*cM[5]; 
  data->AEM_S(3,10) = -0.3535533905932737*cM[6]; 
  data->AEM_S(3,11) = -0.3535533905932737*cM[0]; 
  data->AEM_S(3,12) = -0.3535533905932737*cM[7]; 
  data->AEM_S(3,13) = -0.3535533905932737*cM[1]; 
  data->AEM_S(3,14) = -0.3535533905932737*cM[2]; 
  data->AEM_S(3,15) = -0.3535533905932737*cM[4]; 
  data->AEM_S(4,8) = -0.3535533905932737*cM[4]; 
  data->AEM_S(4,9) = -0.3535533905932737*cM[2]; 
  data->AEM_S(4,10) = -0.3535533905932737*cM[1]; 
  data->AEM_S(4,11) = -0.3535533905932737*cM[7]; 
  data->AEM_S(4,12) = -0.3535533905932737*cM[0]; 
  data->AEM_S(4,13) = -0.3535533905932737*cM[6]; 
  data->AEM_S(4,14) = -0.3535533905932737*cM[5]; 
  data->AEM_S(4,15) = -0.3535533905932737*cM[3]; 
  data->AEM_S(5,8) = -0.3535533905932737*cM[5]; 
  data->AEM_S(5,9) = -0.3535533905932737*cM[3]; 
  data->AEM_S(5,10) = -0.3535533905932737*cM[7]; 
  data->AEM_S(5,11) = -0.3535533905932737*cM[1]; 
  data->AEM_S(5,12) = -0.3535533905932737*cM[6]; 
  data->AEM_S(5,13) = -0.3535533905932737*cM[0]; 
  data->AEM_S(5,14) = -0.3535533905932737*cM[4]; 
  data->AEM_S(5,15) = -0.3535533905932737*cM[2]; 
  data->AEM_S(6,8) = -0.3535533905932737*cM[6]; 
  data->AEM_S(6,9) = -0.3535533905932737*cM[7]; 
  data->AEM_S(6,10) = -0.3535533905932737*cM[3]; 
  data->AEM_S(6,11) = -0.3535533905932737*cM[2]; 
  data->AEM_S(6,12) = -0.3535533905932737*cM[5]; 
  data->AEM_S(6,13) = -0.3535533905932737*cM[4]; 
  data->AEM_S(6,14) = -0.3535533905932737*cM[0]; 
  data->AEM_S(6,15) = -0.3535533905932737*cM[1]; 
  data->AEM_S(7,8) = -0.3535533905932737*cM[7]; 
  data->AEM_S(7,9) = -0.3535533905932737*cM[6]; 
  data->AEM_S(7,10) = -0.3535533905932737*cM[5]; 
  data->AEM_S(7,11) = -0.3535533905932737*cM[4]; 
  data->AEM_S(7,12) = -0.3535533905932737*cM[3]; 
  data->AEM_S(7,13) = -0.3535533905932737*cM[2]; 
  data->AEM_S(7,14) = -0.3535533905932737*cM[1]; 
  data->AEM_S(7,15) = -0.3535533905932737*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(8,0) = 0.3535533905932737*m1Sr[0]; 
  data->AEM_S(8,1) = 0.3535533905932737*m1Sr[1]; 
  data->AEM_S(8,2) = 0.3535533905932737*m1Sr[2]; 
  data->AEM_S(8,3) = 0.3535533905932737*m1Sr[3]; 
  data->AEM_S(8,4) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(8,5) = 0.3535533905932737*m1Sr[5]; 
  data->AEM_S(8,6) = 0.3535533905932737*m1Sr[6]; 
  data->AEM_S(8,7) = 0.3535533905932737*m1Sr[7]; 
  data->AEM_S(9,0) = 0.3535533905932737*m1Sr[1]; 
  data->AEM_S(9,1) = 0.3535533905932737*m1Sr[0]; 
  data->AEM_S(9,2) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(9,3) = 0.3535533905932737*m1Sr[5]; 
  data->AEM_S(9,4) = 0.3535533905932737*m1Sr[2]; 
  data->AEM_S(9,5) = 0.3535533905932737*m1Sr[3]; 
  data->AEM_S(9,6) = 0.3535533905932737*m1Sr[7]; 
  data->AEM_S(9,7) = 0.3535533905932737*m1Sr[6]; 
  data->AEM_S(10,0) = 0.3535533905932737*m1Sr[2]; 
  data->AEM_S(10,1) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(10,2) = 0.3535533905932737*m1Sr[0]; 
  data->AEM_S(10,3) = 0.3535533905932737*m1Sr[6]; 
  data->AEM_S(10,4) = 0.3535533905932737*m1Sr[1]; 
  data->AEM_S(10,5) = 0.3535533905932737*m1Sr[7]; 
  data->AEM_S(10,6) = 0.3535533905932737*m1Sr[3]; 
  data->AEM_S(10,7) = 0.3535533905932737*m1Sr[5]; 
  data->AEM_S(11,0) = 0.3535533905932737*m1Sr[3]; 
  data->AEM_S(11,1) = 0.3535533905932737*m1Sr[5]; 
  data->AEM_S(11,2) = 0.3535533905932737*m1Sr[6]; 
  data->AEM_S(11,3) = 0.3535533905932737*m1Sr[0]; 
  data->AEM_S(11,4) = 0.3535533905932737*m1Sr[7]; 
  data->AEM_S(11,5) = 0.3535533905932737*m1Sr[1]; 
  data->AEM_S(11,6) = 0.3535533905932737*m1Sr[2]; 
  data->AEM_S(11,7) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(12,0) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(12,1) = 0.3535533905932737*m1Sr[2]; 
  data->AEM_S(12,2) = 0.3535533905932737*m1Sr[1]; 
  data->AEM_S(12,3) = 0.3535533905932737*m1Sr[7]; 
  data->AEM_S(12,4) = 0.3535533905932737*m1Sr[0]; 
  data->AEM_S(12,5) = 0.3535533905932737*m1Sr[6]; 
  data->AEM_S(12,6) = 0.3535533905932737*m1Sr[5]; 
  data->AEM_S(12,7) = 0.3535533905932737*m1Sr[3]; 
  data->AEM_S(13,0) = 0.3535533905932737*m1Sr[5]; 
  data->AEM_S(13,1) = 0.3535533905932737*m1Sr[3]; 
  data->AEM_S(13,2) = 0.3535533905932737*m1Sr[7]; 
  data->AEM_S(13,3) = 0.3535533905932737*m1Sr[1]; 
  data->AEM_S(13,4) = 0.3535533905932737*m1Sr[6]; 
  data->AEM_S(13,5) = 0.3535533905932737*m1Sr[0]; 
  data->AEM_S(13,6) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(13,7) = 0.3535533905932737*m1Sr[2]; 
  data->AEM_S(14,0) = 0.3535533905932737*m1Sr[6]; 
  data->AEM_S(14,1) = 0.3535533905932737*m1Sr[7]; 
  data->AEM_S(14,2) = 0.3535533905932737*m1Sr[3]; 
  data->AEM_S(14,3) = 0.3535533905932737*m1Sr[2]; 
  data->AEM_S(14,4) = 0.3535533905932737*m1Sr[5]; 
  data->AEM_S(14,5) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(14,6) = 0.3535533905932737*m1Sr[0]; 
  data->AEM_S(14,7) = 0.3535533905932737*m1Sr[1]; 
  data->AEM_S(15,0) = 0.3535533905932737*m1Sr[7]; 
  data->AEM_S(15,1) = 0.3535533905932737*m1Sr[6]; 
  data->AEM_S(15,2) = 0.3535533905932737*m1Sr[5]; 
  data->AEM_S(15,3) = 0.3535533905932737*m1Sr[4]; 
  data->AEM_S(15,4) = 0.3535533905932737*m1Sr[3]; 
  data->AEM_S(15,5) = 0.3535533905932737*m1Sr[2]; 
  data->AEM_S(15,6) = 0.3535533905932737*m1Sr[1]; 
  data->AEM_S(15,7) = 0.3535533905932737*m1Sr[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(8,8) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(8,9) = 0.7071067811865475*m0r[1]+0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(8,10) = 0.7071067811865475*m0r[2]+0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(8,11) = 0.7071067811865475*m0r[3]+0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(8,12) = 0.7071067811865475*m0r[4]+0.3535533905932737*m0Sr[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(8,13) = 0.7071067811865475*m0r[5]+0.3535533905932737*m0Sr[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(8,14) = 0.7071067811865475*m0r[6]+0.3535533905932737*m0Sr[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(8,15) = 0.7071067811865475*m0r[7]+0.3535533905932737*m0Sr[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(9,8) = 0.7071067811865475*m0r[1]+0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(9,9) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(9,10) = 0.7071067811865475*m0r[4]+0.3535533905932737*m0Sr[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(9,11) = 0.7071067811865475*m0r[5]+0.3535533905932737*m0Sr[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(9,12) = 0.7071067811865475*m0r[2]+0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(9,13) = 0.7071067811865475*m0r[3]+0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(9,14) = 0.7071067811865475*m0r[7]+0.3535533905932737*m0Sr[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(9,15) = 0.7071067811865475*m0r[6]+0.3535533905932737*m0Sr[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(10,8) = 0.7071067811865475*m0r[2]+0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(10,9) = 0.7071067811865475*m0r[4]+0.3535533905932737*m0Sr[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(10,10) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(10,11) = 0.7071067811865475*m0r[6]+0.3535533905932737*m0Sr[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(10,12) = 0.7071067811865475*m0r[1]+0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(10,13) = 0.7071067811865475*m0r[7]+0.3535533905932737*m0Sr[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(10,14) = 0.7071067811865475*m0r[3]+0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(10,15) = 0.7071067811865475*m0r[5]+0.3535533905932737*m0Sr[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(11,8) = 0.7071067811865475*m0r[3]+0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(11,9) = 0.7071067811865475*m0r[5]+0.3535533905932737*m0Sr[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(11,10) = 0.7071067811865475*m0r[6]+0.3535533905932737*m0Sr[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(11,11) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(11,12) = 0.7071067811865475*m0r[7]+0.3535533905932737*m0Sr[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(11,13) = 0.7071067811865475*m0r[1]+0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(11,14) = 0.7071067811865475*m0r[2]+0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(11,15) = 0.7071067811865475*m0r[4]+0.3535533905932737*m0Sr[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(12,8) = 0.7071067811865475*m0r[4]+0.3535533905932737*m0Sr[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(12,9) = 0.7071067811865475*m0r[2]+0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(12,10) = 0.7071067811865475*m0r[1]+0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(12,11) = 0.7071067811865475*m0r[7]+0.3535533905932737*m0Sr[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(12,12) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(12,13) = 0.7071067811865475*m0r[6]+0.3535533905932737*m0Sr[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(12,14) = 0.7071067811865475*m0r[5]+0.3535533905932737*m0Sr[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(12,15) = 0.7071067811865475*m0r[3]+0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(13,8) = 0.7071067811865475*m0r[5]+0.3535533905932737*m0Sr[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(13,9) = 0.7071067811865475*m0r[3]+0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(13,10) = 0.7071067811865475*m0r[7]+0.3535533905932737*m0Sr[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(13,11) = 0.7071067811865475*m0r[1]+0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(13,12) = 0.7071067811865475*m0r[6]+0.3535533905932737*m0Sr[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(13,13) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(13,14) = 0.7071067811865475*m0r[4]+0.3535533905932737*m0Sr[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(13,15) = 0.7071067811865475*m0r[2]+0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(14,8) = 0.7071067811865475*m0r[6]+0.3535533905932737*m0Sr[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(14,9) = 0.7071067811865475*m0r[7]+0.3535533905932737*m0Sr[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(14,10) = 0.7071067811865475*m0r[3]+0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(14,11) = 0.7071067811865475*m0r[2]+0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(14,12) = 0.7071067811865475*m0r[5]+0.3535533905932737*m0Sr[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(14,13) = 0.7071067811865475*m0r[4]+0.3535533905932737*m0Sr[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(14,14) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(14,15) = 0.7071067811865475*m0r[1]+0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(15,8) = 0.7071067811865475*m0r[7]+0.3535533905932737*m0Sr[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(15,9) = 0.7071067811865475*m0r[6]+0.3535533905932737*m0Sr[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(15,10) = 0.7071067811865475*m0r[5]+0.3535533905932737*m0Sr[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(15,11) = 0.7071067811865475*m0r[4]+0.3535533905932737*m0Sr[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(15,12) = 0.7071067811865475*m0r[3]+0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(15,13) = 0.7071067811865475*m0r[2]+0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(15,14) = 0.7071067811865475*m0r[1]+0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(15,15) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m2Sr[0],m2Sr[1],m2Sr[2],m2Sr[3],m2Sr[4],m2Sr[5],m2Sr[6],m2Sr[7]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,8,1) = data->u_S.segment<8>(0); 
 
  Eigen::Map<VectorXd>(vtSq,8,1) = data->u_S.segment<8>(8); 
 
} 
 
void GkSelfPrimMoments3x2vSer_P2(binOpData_t* data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (2.371708245126284*m0[19]+2.371708245126284*m0[18]+2.371708245126284*m0[17]-1.369306393762915*m0[16]-1.369306393762915*m0[15]-1.369306393762915*m0[14]-1.369306393762915*m0[13]-1.369306393762915*m0[12]-1.369306393762915*m0[11]-1.837117307087383*m0[10]+0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (2.371708245126284*m0[19]+2.371708245126284*m0[18]+2.371708245126284*m0[17]-1.369306393762915*m0[16]-1.369306393762915*m0[15]-1.369306393762915*m0[14]-1.369306393762915*m0[13]-1.369306393762915*m0[12]-1.369306393762915*m0[11]-1.837117307087383*m0[10]+0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (2.371708245126284*m0[19]+2.371708245126284*m0[18]+2.371708245126284*m0[17]-1.369306393762915*m0[16]-1.369306393762915*m0[15]-1.369306393762915*m0[14]-1.369306393762915*m0[13]-1.369306393762915*m0[12]-1.369306393762915*m0[11]-1.837117307087383*m0[10]+0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (2.371708245126284*m0[19]+2.371708245126284*m0[18]+2.371708245126284*m0[17]-1.369306393762915*m0[16]-1.369306393762915*m0[15]-1.369306393762915*m0[14]-1.369306393762915*m0[13]-1.369306393762915*m0[12]-1.369306393762915*m0[11]-1.837117307087383*m0[10]+0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-2.371708245126284*m0[19])-2.371708245126284*m0[18]+2.371708245126284*m0[17]-1.369306393762915*m0[16]+1.369306393762915*m0[15]-1.369306393762915*m0[14]-1.369306393762915*m0[13]+1.369306393762915*m0[12]-1.369306393762915*m0[11]+1.837117307087383*m0[10]+0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-2.371708245126284*m0[19])-2.371708245126284*m0[18]+2.371708245126284*m0[17]-1.369306393762915*m0[16]+1.369306393762915*m0[15]-1.369306393762915*m0[14]-1.369306393762915*m0[13]+1.369306393762915*m0[12]-1.369306393762915*m0[11]+1.837117307087383*m0[10]+0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-2.371708245126284*m0[19])-2.371708245126284*m0[18]+2.371708245126284*m0[17]-1.369306393762915*m0[16]+1.369306393762915*m0[15]-1.369306393762915*m0[14]-1.369306393762915*m0[13]+1.369306393762915*m0[12]-1.369306393762915*m0[11]+1.837117307087383*m0[10]+0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-2.371708245126284*m0[19])-2.371708245126284*m0[18]+2.371708245126284*m0[17]-1.369306393762915*m0[16]+1.369306393762915*m0[15]-1.369306393762915*m0[14]-1.369306393762915*m0[13]+1.369306393762915*m0[12]-1.369306393762915*m0[11]+1.837117307087383*m0[10]+0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[20]; 
  double m1r[20]; 
  double m2r[20]; 
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
    m0r[12] = 0.0; 
    m0r[13] = 0.0; 
    m0r[14] = 0.0; 
    m0r[15] = 0.0; 
    m0r[16] = 0.0; 
    m0r[17] = 0.0; 
    m0r[18] = 0.0; 
    m0r[19] = 0.0; 
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
    m2r[10] = 0.0; 
    m2r[11] = 0.0; 
    m2r[12] = 0.0; 
    m2r[13] = 0.0; 
    m2r[14] = 0.0; 
    m2r[15] = 0.0; 
    m2r[16] = 0.0; 
    m2r[17] = 0.0; 
    m2r[18] = 0.0; 
    m2r[19] = 0.0; 
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
    m0r[12] = m0[12]; 
    m0r[13] = m0[13]; 
    m0r[14] = m0[14]; 
    m0r[15] = m0[15]; 
    m0r[16] = m0[16]; 
    m0r[17] = m0[17]; 
    m0r[18] = m0[18]; 
    m0r[19] = m0[19]; 
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
    m2r[10] = m2[10]; 
    m2r[11] = m2[11]; 
    m2r[12] = m2[12]; 
    m2r[13] = m2[13]; 
    m2r[14] = m2[14]; 
    m2r[15] = m2[15]; 
    m2r[16] = m2[16]; 
    m2r[17] = m2[17]; 
    m2r[18] = m2[18]; 
    m2r[19] = m2[19]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(40,40); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.3535533905932737*m0r[0]; 
  data->AEM_S(0,1) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(0,2) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(0,3) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(0,4) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(0,5) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(0,6) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(0,7) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(0,8) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(0,9) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(0,10) = 0.3535533905932737*m0r[10]; 
  data->AEM_S(0,11) = 0.3535533905932737*m0r[11]; 
  data->AEM_S(0,12) = 0.3535533905932737*m0r[12]; 
  data->AEM_S(0,13) = 0.3535533905932737*m0r[13]; 
  data->AEM_S(0,14) = 0.3535533905932737*m0r[14]; 
  data->AEM_S(0,15) = 0.3535533905932737*m0r[15]; 
  data->AEM_S(0,16) = 0.3535533905932737*m0r[16]; 
  data->AEM_S(0,17) = 0.3535533905932737*m0r[17]; 
  data->AEM_S(0,18) = 0.3535533905932737*m0r[18]; 
  data->AEM_S(0,19) = 0.3535533905932737*m0r[19]; 
  data->AEM_S(1,0) = 0.3535533905932737*m0r[1]; 
  data->AEM_S(1,1) = 0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(1,2) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(1,3) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(1,4) = 0.3162277660168379*m0r[11]+0.3535533905932737*m0r[2]; 
  data->AEM_S(1,5) = 0.3162277660168379*m0r[13]+0.3535533905932737*m0r[3]; 
  data->AEM_S(1,6) = 0.3535533905932737*m0r[10]; 
  data->AEM_S(1,7) = 0.3162277660168379*m0r[1]; 
  data->AEM_S(1,8) = 0.3535533905932737*m0r[12]; 
  data->AEM_S(1,9) = 0.3535533905932737*m0r[15]; 
  data->AEM_S(1,10) = 0.3162277660168379*m0r[17]+0.3535533905932737*m0r[6]; 
  data->AEM_S(1,11) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(1,12) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(1,13) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(1,14) = 0.3535533905932737*m0r[18]; 
  data->AEM_S(1,15) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(1,16) = 0.3535533905932737*m0r[19]; 
  data->AEM_S(1,17) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(1,18) = 0.3535533905932737*m0r[14]; 
  data->AEM_S(1,19) = 0.3535533905932737*m0r[16]; 
  data->AEM_S(2,0) = 0.3535533905932737*m0r[2]; 
  data->AEM_S(2,1) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(2,2) = 0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(2,3) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(2,4) = 0.3162277660168379*m0r[12]+0.3535533905932737*m0r[1]; 
  data->AEM_S(2,5) = 0.3535533905932737*m0r[10]; 
  data->AEM_S(2,6) = 0.3162277660168379*m0r[14]+0.3535533905932737*m0r[3]; 
  data->AEM_S(2,7) = 0.3535533905932737*m0r[11]; 
  data->AEM_S(2,8) = 0.3162277660168379*m0r[2]; 
  data->AEM_S(2,9) = 0.3535533905932737*m0r[16]; 
  data->AEM_S(2,10) = 0.3162277660168379*m0r[18]+0.3535533905932737*m0r[5]; 
  data->AEM_S(2,11) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(2,12) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(2,13) = 0.3535533905932737*m0r[17]; 
  data->AEM_S(2,14) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(2,15) = 0.3535533905932737*m0r[19]; 
  data->AEM_S(2,16) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(2,17) = 0.3535533905932737*m0r[13]; 
  data->AEM_S(2,18) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(2,19) = 0.3535533905932737*m0r[15]; 
  data->AEM_S(3,0) = 0.3535533905932737*m0r[3]; 
  data->AEM_S(3,1) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(3,2) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(3,3) = 0.3162277660168379*m0r[9]+0.3535533905932737*m0r[0]; 
  data->AEM_S(3,4) = 0.3535533905932737*m0r[10]; 
  data->AEM_S(3,5) = 0.3162277660168379*m0r[15]+0.3535533905932737*m0r[1]; 
  data->AEM_S(3,6) = 0.3162277660168379*m0r[16]+0.3535533905932737*m0r[2]; 
  data->AEM_S(3,7) = 0.3535533905932737*m0r[13]; 
  data->AEM_S(3,8) = 0.3535533905932737*m0r[14]; 
  data->AEM_S(3,9) = 0.3162277660168379*m0r[3]; 
  data->AEM_S(3,10) = 0.3162277660168379*m0r[19]+0.3535533905932737*m0r[4]; 
  data->AEM_S(3,11) = 0.3535533905932737*m0r[17]; 
  data->AEM_S(3,12) = 0.3535533905932737*m0r[18]; 
  data->AEM_S(3,13) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(3,14) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(3,15) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(3,16) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(3,17) = 0.3535533905932737*m0r[11]; 
  data->AEM_S(3,18) = 0.3535533905932737*m0r[12]; 
  data->AEM_S(3,19) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(4,0) = 0.3535533905932737*m0r[4]; 
  data->AEM_S(4,1) = 0.3162277660168379*m0r[11]+0.3535533905932737*m0r[2]; 
  data->AEM_S(4,2) = 0.3162277660168379*m0r[12]+0.3535533905932737*m0r[1]; 
  data->AEM_S(4,3) = 0.3535533905932737*m0r[10]; 
  data->AEM_S(4,4) = 0.3162277660168379*m0r[8]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(4,5) = 0.3162277660168379*m0r[17]+0.3535533905932737*m0r[6]; 
  data->AEM_S(4,6) = 0.3162277660168379*m0r[18]+0.3535533905932737*m0r[5]; 
  data->AEM_S(4,7) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(4,8) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(4,9) = 0.3535533905932737*m0r[19]; 
  data->AEM_S(4,10) = 0.3162277660168379*m0r[14]+0.3162277660168379*m0r[13]+0.3535533905932737*m0r[3]; 
  data->AEM_S(4,11) = 0.2828427124746191*m0r[12]+0.3162277660168379*m0r[1]; 
  data->AEM_S(4,12) = 0.2828427124746191*m0r[11]+0.3162277660168379*m0r[2]; 
  data->AEM_S(4,13) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(4,14) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(4,15) = 0.3535533905932737*m0r[16]; 
  data->AEM_S(4,16) = 0.3535533905932737*m0r[15]; 
  data->AEM_S(4,17) = 0.2828427124746191*m0r[18]+0.3162277660168379*m0r[5]; 
  data->AEM_S(4,18) = 0.2828427124746191*m0r[17]+0.3162277660168379*m0r[6]; 
  data->AEM_S(4,19) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(5,0) = 0.3535533905932737*m0r[5]; 
  data->AEM_S(5,1) = 0.3162277660168379*m0r[13]+0.3535533905932737*m0r[3]; 
  data->AEM_S(5,2) = 0.3535533905932737*m0r[10]; 
  data->AEM_S(5,3) = 0.3162277660168379*m0r[15]+0.3535533905932737*m0r[1]; 
  data->AEM_S(5,4) = 0.3162277660168379*m0r[17]+0.3535533905932737*m0r[6]; 
  data->AEM_S(5,5) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(5,6) = 0.3162277660168379*m0r[19]+0.3535533905932737*m0r[4]; 
  data->AEM_S(5,7) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(5,8) = 0.3535533905932737*m0r[18]; 
  data->AEM_S(5,9) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(5,10) = 0.3162277660168379*m0r[16]+0.3162277660168379*m0r[11]+0.3535533905932737*m0r[2]; 
  data->AEM_S(5,11) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(5,12) = 0.3535533905932737*m0r[14]; 
  data->AEM_S(5,13) = 0.2828427124746191*m0r[15]+0.3162277660168379*m0r[1]; 
  data->AEM_S(5,14) = 0.3535533905932737*m0r[12]; 
  data->AEM_S(5,15) = 0.2828427124746191*m0r[13]+0.3162277660168379*m0r[3]; 
  data->AEM_S(5,16) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(5,17) = 0.2828427124746191*m0r[19]+0.3162277660168379*m0r[4]; 
  data->AEM_S(5,18) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(5,19) = 0.2828427124746191*m0r[17]+0.3162277660168379*m0r[6]; 
  data->AEM_S(6,0) = 0.3535533905932737*m0r[6]; 
  data->AEM_S(6,1) = 0.3535533905932737*m0r[10]; 
  data->AEM_S(6,2) = 0.3162277660168379*m0r[14]+0.3535533905932737*m0r[3]; 
  data->AEM_S(6,3) = 0.3162277660168379*m0r[16]+0.3535533905932737*m0r[2]; 
  data->AEM_S(6,4) = 0.3162277660168379*m0r[18]+0.3535533905932737*m0r[5]; 
  data->AEM_S(6,5) = 0.3162277660168379*m0r[19]+0.3535533905932737*m0r[4]; 
  data->AEM_S(6,6) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(6,7) = 0.3535533905932737*m0r[17]; 
  data->AEM_S(6,8) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(6,9) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(6,10) = 0.3162277660168379*m0r[15]+0.3162277660168379*m0r[12]+0.3535533905932737*m0r[1]; 
  data->AEM_S(6,11) = 0.3535533905932737*m0r[13]; 
  data->AEM_S(6,12) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(6,13) = 0.3535533905932737*m0r[11]; 
  data->AEM_S(6,14) = 0.2828427124746191*m0r[16]+0.3162277660168379*m0r[2]; 
  data->AEM_S(6,15) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(6,16) = 0.2828427124746191*m0r[14]+0.3162277660168379*m0r[3]; 
  data->AEM_S(6,17) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(6,18) = 0.2828427124746191*m0r[19]+0.3162277660168379*m0r[4]; 
  data->AEM_S(6,19) = 0.2828427124746191*m0r[18]+0.3162277660168379*m0r[5]; 
  data->AEM_S(7,0) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(7,1) = 0.3162277660168379*m0r[1]; 
  data->AEM_S(7,2) = 0.3535533905932737*m0r[11]; 
  data->AEM_S(7,3) = 0.3535533905932737*m0r[13]; 
  data->AEM_S(7,4) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(7,5) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(7,6) = 0.3535533905932737*m0r[17]; 
  data->AEM_S(7,7) = 0.2258769757263128*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(7,10) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(7,11) = 0.2258769757263128*m0r[11]+0.3535533905932737*m0r[2]; 
  data->AEM_S(7,12) = 0.3162277660168379*m0r[12]; 
  data->AEM_S(7,13) = 0.2258769757263128*m0r[13]+0.3535533905932737*m0r[3]; 
  data->AEM_S(7,15) = 0.3162277660168379*m0r[15]; 
  data->AEM_S(7,17) = 0.2258769757263128*m0r[17]+0.3535533905932737*m0r[6]; 
  data->AEM_S(7,18) = 0.3162277660168379*m0r[18]; 
  data->AEM_S(7,19) = 0.3162277660168379*m0r[19]; 
  data->AEM_S(8,0) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(8,1) = 0.3535533905932737*m0r[12]; 
  data->AEM_S(8,2) = 0.3162277660168379*m0r[2]; 
  data->AEM_S(8,3) = 0.3535533905932737*m0r[14]; 
  data->AEM_S(8,4) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(8,5) = 0.3535533905932737*m0r[18]; 
  data->AEM_S(8,6) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(8,8) = 0.2258769757263128*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(8,10) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(8,11) = 0.3162277660168379*m0r[11]; 
  data->AEM_S(8,12) = 0.2258769757263128*m0r[12]+0.3535533905932737*m0r[1]; 
  data->AEM_S(8,14) = 0.2258769757263128*m0r[14]+0.3535533905932737*m0r[3]; 
  data->AEM_S(8,16) = 0.3162277660168379*m0r[16]; 
  data->AEM_S(8,17) = 0.3162277660168379*m0r[17]; 
  data->AEM_S(8,18) = 0.2258769757263128*m0r[18]+0.3535533905932737*m0r[5]; 
  data->AEM_S(8,19) = 0.3162277660168379*m0r[19]; 
  data->AEM_S(9,0) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(9,1) = 0.3535533905932737*m0r[15]; 
  data->AEM_S(9,2) = 0.3535533905932737*m0r[16]; 
  data->AEM_S(9,3) = 0.3162277660168379*m0r[3]; 
  data->AEM_S(9,4) = 0.3535533905932737*m0r[19]; 
  data->AEM_S(9,5) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(9,6) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(9,9) = 0.2258769757263128*m0r[9]+0.3535533905932737*m0r[0]; 
  data->AEM_S(9,10) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(9,13) = 0.3162277660168379*m0r[13]; 
  data->AEM_S(9,14) = 0.3162277660168379*m0r[14]; 
  data->AEM_S(9,15) = 0.2258769757263128*m0r[15]+0.3535533905932737*m0r[1]; 
  data->AEM_S(9,16) = 0.2258769757263128*m0r[16]+0.3535533905932737*m0r[2]; 
  data->AEM_S(9,17) = 0.3162277660168379*m0r[17]; 
  data->AEM_S(9,18) = 0.3162277660168379*m0r[18]; 
  data->AEM_S(9,19) = 0.2258769757263128*m0r[19]+0.3535533905932737*m0r[4]; 
  data->AEM_S(10,0) = 0.3535533905932737*m0r[10]; 
  data->AEM_S(10,1) = 0.3162277660168379*m0r[17]+0.3535533905932737*m0r[6]; 
  data->AEM_S(10,2) = 0.3162277660168379*m0r[18]+0.3535533905932737*m0r[5]; 
  data->AEM_S(10,3) = 0.3162277660168379*m0r[19]+0.3535533905932737*m0r[4]; 
  data->AEM_S(10,4) = 0.3162277660168379*m0r[14]+0.3162277660168379*m0r[13]+0.3535533905932737*m0r[3]; 
  data->AEM_S(10,5) = 0.3162277660168379*m0r[16]+0.3162277660168379*m0r[11]+0.3535533905932737*m0r[2]; 
  data->AEM_S(10,6) = 0.3162277660168379*m0r[15]+0.3162277660168379*m0r[12]+0.3535533905932737*m0r[1]; 
  data->AEM_S(10,7) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(10,8) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(10,9) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(10,10) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[8]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(10,11) = 0.282842712474619*m0r[18]+0.3162277660168379*m0r[5]; 
  data->AEM_S(10,12) = 0.282842712474619*m0r[17]+0.3162277660168379*m0r[6]; 
  data->AEM_S(10,13) = 0.282842712474619*m0r[19]+0.3162277660168379*m0r[4]; 
  data->AEM_S(10,14) = 0.282842712474619*m0r[19]+0.3162277660168379*m0r[4]; 
  data->AEM_S(10,15) = 0.282842712474619*m0r[17]+0.3162277660168379*m0r[6]; 
  data->AEM_S(10,16) = 0.282842712474619*m0r[18]+0.3162277660168379*m0r[5]; 
  data->AEM_S(10,17) = 0.282842712474619*m0r[15]+0.282842712474619*m0r[12]+0.3162277660168379*m0r[1]; 
  data->AEM_S(10,18) = 0.282842712474619*m0r[16]+0.282842712474619*m0r[11]+0.3162277660168379*m0r[2]; 
  data->AEM_S(10,19) = 0.282842712474619*m0r[14]+0.282842712474619*m0r[13]+0.3162277660168379*m0r[3]; 
  data->AEM_S(11,0) = 0.3535533905932737*m0r[11]; 
  data->AEM_S(11,1) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(11,2) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(11,3) = 0.3535533905932737*m0r[17]; 
  data->AEM_S(11,4) = 0.2828427124746191*m0r[12]+0.3162277660168379*m0r[1]; 
  data->AEM_S(11,5) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(11,6) = 0.3535533905932737*m0r[13]; 
  data->AEM_S(11,7) = 0.2258769757263128*m0r[11]+0.3535533905932737*m0r[2]; 
  data->AEM_S(11,8) = 0.3162277660168379*m0r[11]; 
  data->AEM_S(11,10) = 0.282842712474619*m0r[18]+0.3162277660168379*m0r[5]; 
  data->AEM_S(11,11) = 0.3162277660168379*m0r[8]+0.2258769757263128*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(11,12) = 0.2828427124746191*m0r[4]; 
  data->AEM_S(11,13) = 0.2258769757263128*m0r[17]+0.3535533905932737*m0r[6]; 
  data->AEM_S(11,14) = 0.3162277660168379*m0r[17]; 
  data->AEM_S(11,15) = 0.3162277660168379*m0r[19]; 
  data->AEM_S(11,17) = 0.3162277660168379*m0r[14]+0.2258769757263128*m0r[13]+0.3535533905932737*m0r[3]; 
  data->AEM_S(11,18) = 0.282842712474619*m0r[10]; 
  data->AEM_S(11,19) = 0.3162277660168379*m0r[15]; 
  data->AEM_S(12,0) = 0.3535533905932737*m0r[12]; 
  data->AEM_S(12,1) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(12,2) = 0.3162277660168379*m0r[4]; 
  data->AEM_S(12,3) = 0.3535533905932737*m0r[18]; 
  data->AEM_S(12,4) = 0.2828427124746191*m0r[11]+0.3162277660168379*m0r[2]; 
  data->AEM_S(12,5) = 0.3535533905932737*m0r[14]; 
  data->AEM_S(12,6) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(12,7) = 0.3162277660168379*m0r[12]; 
  data->AEM_S(12,8) = 0.2258769757263128*m0r[12]+0.3535533905932737*m0r[1]; 
  data->AEM_S(12,10) = 0.282842712474619*m0r[17]+0.3162277660168379*m0r[6]; 
  data->AEM_S(12,11) = 0.2828427124746191*m0r[4]; 
  data->AEM_S(12,12) = 0.2258769757263128*m0r[8]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(12,13) = 0.3162277660168379*m0r[18]; 
  data->AEM_S(12,14) = 0.2258769757263128*m0r[18]+0.3535533905932737*m0r[5]; 
  data->AEM_S(12,16) = 0.3162277660168379*m0r[19]; 
  data->AEM_S(12,17) = 0.282842712474619*m0r[10]; 
  data->AEM_S(12,18) = 0.2258769757263128*m0r[14]+0.3162277660168379*m0r[13]+0.3535533905932737*m0r[3]; 
  data->AEM_S(12,19) = 0.3162277660168379*m0r[16]; 
  data->AEM_S(13,0) = 0.3535533905932737*m0r[13]; 
  data->AEM_S(13,1) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(13,2) = 0.3535533905932737*m0r[17]; 
  data->AEM_S(13,3) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(13,4) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(13,5) = 0.2828427124746191*m0r[15]+0.3162277660168379*m0r[1]; 
  data->AEM_S(13,6) = 0.3535533905932737*m0r[11]; 
  data->AEM_S(13,7) = 0.2258769757263128*m0r[13]+0.3535533905932737*m0r[3]; 
  data->AEM_S(13,9) = 0.3162277660168379*m0r[13]; 
  data->AEM_S(13,10) = 0.282842712474619*m0r[19]+0.3162277660168379*m0r[4]; 
  data->AEM_S(13,11) = 0.2258769757263128*m0r[17]+0.3535533905932737*m0r[6]; 
  data->AEM_S(13,12) = 0.3162277660168379*m0r[18]; 
  data->AEM_S(13,13) = 0.3162277660168379*m0r[9]+0.2258769757263128*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(13,15) = 0.2828427124746191*m0r[5]; 
  data->AEM_S(13,16) = 0.3162277660168379*m0r[17]; 
  data->AEM_S(13,17) = 0.3162277660168379*m0r[16]+0.2258769757263128*m0r[11]+0.3535533905932737*m0r[2]; 
  data->AEM_S(13,18) = 0.3162277660168379*m0r[12]; 
  data->AEM_S(13,19) = 0.282842712474619*m0r[10]; 
  data->AEM_S(14,0) = 0.3535533905932737*m0r[14]; 
  data->AEM_S(14,1) = 0.3535533905932737*m0r[18]; 
  data->AEM_S(14,2) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(14,3) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(14,4) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(14,5) = 0.3535533905932737*m0r[12]; 
  data->AEM_S(14,6) = 0.2828427124746191*m0r[16]+0.3162277660168379*m0r[2]; 
  data->AEM_S(14,8) = 0.2258769757263128*m0r[14]+0.3535533905932737*m0r[3]; 
  data->AEM_S(14,9) = 0.3162277660168379*m0r[14]; 
  data->AEM_S(14,10) = 0.282842712474619*m0r[19]+0.3162277660168379*m0r[4]; 
  data->AEM_S(14,11) = 0.3162277660168379*m0r[17]; 
  data->AEM_S(14,12) = 0.2258769757263128*m0r[18]+0.3535533905932737*m0r[5]; 
  data->AEM_S(14,14) = 0.3162277660168379*m0r[9]+0.2258769757263128*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(14,15) = 0.3162277660168379*m0r[18]; 
  data->AEM_S(14,16) = 0.2828427124746191*m0r[6]; 
  data->AEM_S(14,17) = 0.3162277660168379*m0r[11]; 
  data->AEM_S(14,18) = 0.3162277660168379*m0r[15]+0.2258769757263128*m0r[12]+0.3535533905932737*m0r[1]; 
  data->AEM_S(14,19) = 0.282842712474619*m0r[10]; 
  data->AEM_S(15,0) = 0.3535533905932737*m0r[15]; 
  data->AEM_S(15,1) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(15,2) = 0.3535533905932737*m0r[19]; 
  data->AEM_S(15,3) = 0.3162277660168379*m0r[5]; 
  data->AEM_S(15,4) = 0.3535533905932737*m0r[16]; 
  data->AEM_S(15,5) = 0.2828427124746191*m0r[13]+0.3162277660168379*m0r[3]; 
  data->AEM_S(15,6) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(15,7) = 0.3162277660168379*m0r[15]; 
  data->AEM_S(15,9) = 0.2258769757263128*m0r[15]+0.3535533905932737*m0r[1]; 
  data->AEM_S(15,10) = 0.282842712474619*m0r[17]+0.3162277660168379*m0r[6]; 
  data->AEM_S(15,11) = 0.3162277660168379*m0r[19]; 
  data->AEM_S(15,13) = 0.2828427124746191*m0r[5]; 
  data->AEM_S(15,14) = 0.3162277660168379*m0r[18]; 
  data->AEM_S(15,15) = 0.2258769757263128*m0r[9]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(15,16) = 0.2258769757263128*m0r[19]+0.3535533905932737*m0r[4]; 
  data->AEM_S(15,17) = 0.282842712474619*m0r[10]; 
  data->AEM_S(15,18) = 0.3162277660168379*m0r[14]; 
  data->AEM_S(15,19) = 0.2258769757263128*m0r[16]+0.3162277660168379*m0r[11]+0.3535533905932737*m0r[2]; 
  data->AEM_S(16,0) = 0.3535533905932737*m0r[16]; 
  data->AEM_S(16,1) = 0.3535533905932737*m0r[19]; 
  data->AEM_S(16,2) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(16,3) = 0.3162277660168379*m0r[6]; 
  data->AEM_S(16,4) = 0.3535533905932737*m0r[15]; 
  data->AEM_S(16,5) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(16,6) = 0.2828427124746191*m0r[14]+0.3162277660168379*m0r[3]; 
  data->AEM_S(16,8) = 0.3162277660168379*m0r[16]; 
  data->AEM_S(16,9) = 0.2258769757263128*m0r[16]+0.3535533905932737*m0r[2]; 
  data->AEM_S(16,10) = 0.282842712474619*m0r[18]+0.3162277660168379*m0r[5]; 
  data->AEM_S(16,12) = 0.3162277660168379*m0r[19]; 
  data->AEM_S(16,13) = 0.3162277660168379*m0r[17]; 
  data->AEM_S(16,14) = 0.2828427124746191*m0r[6]; 
  data->AEM_S(16,15) = 0.2258769757263128*m0r[19]+0.3535533905932737*m0r[4]; 
  data->AEM_S(16,16) = 0.2258769757263128*m0r[9]+0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  data->AEM_S(16,17) = 0.3162277660168379*m0r[13]; 
  data->AEM_S(16,18) = 0.282842712474619*m0r[10]; 
  data->AEM_S(16,19) = 0.2258769757263128*m0r[15]+0.3162277660168379*m0r[12]+0.3535533905932737*m0r[1]; 
  data->AEM_S(17,0) = 0.3535533905932737*m0r[17]; 
  data->AEM_S(17,1) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(17,2) = 0.3535533905932737*m0r[13]; 
  data->AEM_S(17,3) = 0.3535533905932737*m0r[11]; 
  data->AEM_S(17,4) = 0.2828427124746191*m0r[18]+0.3162277660168379*m0r[5]; 
  data->AEM_S(17,5) = 0.2828427124746191*m0r[19]+0.3162277660168379*m0r[4]; 
  data->AEM_S(17,6) = 0.3535533905932737*m0r[7]; 
  data->AEM_S(17,7) = 0.2258769757263128*m0r[17]+0.3535533905932737*m0r[6]; 
  data->AEM_S(17,8) = 0.3162277660168379*m0r[17]; 
  data->AEM_S(17,9) = 0.3162277660168379*m0r[17]; 
  data->AEM_S(17,10) = 0.282842712474619*m0r[15]+0.282842712474619*m0r[12]+0.3162277660168379*m0r[1]; 
  data->AEM_S(17,11) = 0.3162277660168379*m0r[14]+0.2258769757263128*m0r[13]+0.3535533905932737*m0r[3]; 
  data->AEM_S(17,12) = 0.282842712474619*m0r[10]; 
  data->AEM_S(17,13) = 0.3162277660168379*m0r[16]+0.2258769757263128*m0r[11]+0.3535533905932737*m0r[2]; 
  data->AEM_S(17,14) = 0.3162277660168379*m0r[11]; 
  data->AEM_S(17,15) = 0.282842712474619*m0r[10]; 
  data->AEM_S(17,16) = 0.3162277660168379*m0r[13]; 
  data->AEM_S(17,17) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[8]+0.2258769757263128*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(17,18) = 0.2529822128134704*m0r[19]+0.2828427124746191*m0r[4]; 
  data->AEM_S(17,19) = 0.2529822128134704*m0r[18]+0.2828427124746191*m0r[5]; 
  data->AEM_S(18,0) = 0.3535533905932737*m0r[18]; 
  data->AEM_S(18,1) = 0.3535533905932737*m0r[14]; 
  data->AEM_S(18,2) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(18,3) = 0.3535533905932737*m0r[12]; 
  data->AEM_S(18,4) = 0.2828427124746191*m0r[17]+0.3162277660168379*m0r[6]; 
  data->AEM_S(18,5) = 0.3535533905932737*m0r[8]; 
  data->AEM_S(18,6) = 0.2828427124746191*m0r[19]+0.3162277660168379*m0r[4]; 
  data->AEM_S(18,7) = 0.3162277660168379*m0r[18]; 
  data->AEM_S(18,8) = 0.2258769757263128*m0r[18]+0.3535533905932737*m0r[5]; 
  data->AEM_S(18,9) = 0.3162277660168379*m0r[18]; 
  data->AEM_S(18,10) = 0.282842712474619*m0r[16]+0.282842712474619*m0r[11]+0.3162277660168379*m0r[2]; 
  data->AEM_S(18,11) = 0.282842712474619*m0r[10]; 
  data->AEM_S(18,12) = 0.2258769757263128*m0r[14]+0.3162277660168379*m0r[13]+0.3535533905932737*m0r[3]; 
  data->AEM_S(18,13) = 0.3162277660168379*m0r[12]; 
  data->AEM_S(18,14) = 0.3162277660168379*m0r[15]+0.2258769757263128*m0r[12]+0.3535533905932737*m0r[1]; 
  data->AEM_S(18,15) = 0.3162277660168379*m0r[14]; 
  data->AEM_S(18,16) = 0.282842712474619*m0r[10]; 
  data->AEM_S(18,17) = 0.2529822128134704*m0r[19]+0.2828427124746191*m0r[4]; 
  data->AEM_S(18,18) = 0.3162277660168379*m0r[9]+0.2258769757263128*m0r[8]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  data->AEM_S(18,19) = 0.2529822128134704*m0r[17]+0.2828427124746191*m0r[6]; 
  data->AEM_S(19,0) = 0.3535533905932737*m0r[19]; 
  data->AEM_S(19,1) = 0.3535533905932737*m0r[16]; 
  data->AEM_S(19,2) = 0.3535533905932737*m0r[15]; 
  data->AEM_S(19,3) = 0.3162277660168379*m0r[10]; 
  data->AEM_S(19,4) = 0.3535533905932737*m0r[9]; 
  data->AEM_S(19,5) = 0.2828427124746191*m0r[17]+0.3162277660168379*m0r[6]; 
  data->AEM_S(19,6) = 0.2828427124746191*m0r[18]+0.3162277660168379*m0r[5]; 
  data->AEM_S(19,7) = 0.3162277660168379*m0r[19]; 
  data->AEM_S(19,8) = 0.3162277660168379*m0r[19]; 
  data->AEM_S(19,9) = 0.2258769757263128*m0r[19]+0.3535533905932737*m0r[4]; 
  data->AEM_S(19,10) = 0.282842712474619*m0r[14]+0.282842712474619*m0r[13]+0.3162277660168379*m0r[3]; 
  data->AEM_S(19,11) = 0.3162277660168379*m0r[15]; 
  data->AEM_S(19,12) = 0.3162277660168379*m0r[16]; 
  data->AEM_S(19,13) = 0.282842712474619*m0r[10]; 
  data->AEM_S(19,14) = 0.282842712474619*m0r[10]; 
  data->AEM_S(19,15) = 0.2258769757263128*m0r[16]+0.3162277660168379*m0r[11]+0.3535533905932737*m0r[2]; 
  data->AEM_S(19,16) = 0.2258769757263128*m0r[15]+0.3162277660168379*m0r[12]+0.3535533905932737*m0r[1]; 
  data->AEM_S(19,17) = 0.2529822128134704*m0r[18]+0.2828427124746191*m0r[5]; 
  data->AEM_S(19,18) = 0.2529822128134704*m0r[17]+0.2828427124746191*m0r[6]; 
  data->AEM_S(19,19) = 0.2258769757263128*m0r[9]+0.3162277660168379*m0r[8]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,20) = -0.3535533905932737*cM[0]; 
  data->AEM_S(0,21) = -0.3535533905932737*cM[1]; 
  data->AEM_S(0,22) = -0.3535533905932737*cM[2]; 
  data->AEM_S(0,23) = -0.3535533905932737*cM[3]; 
  data->AEM_S(0,24) = -0.3535533905932737*cM[4]; 
  data->AEM_S(0,25) = -0.3535533905932737*cM[5]; 
  data->AEM_S(0,26) = -0.3535533905932737*cM[6]; 
  data->AEM_S(0,27) = -0.3535533905932737*cM[7]; 
  data->AEM_S(0,28) = -0.3535533905932737*cM[8]; 
  data->AEM_S(0,29) = -0.3535533905932737*cM[9]; 
  data->AEM_S(0,30) = -0.3535533905932737*cM[10]; 
  data->AEM_S(0,31) = -0.3535533905932737*cM[11]; 
  data->AEM_S(0,32) = -0.3535533905932737*cM[12]; 
  data->AEM_S(0,33) = -0.3535533905932737*cM[13]; 
  data->AEM_S(0,34) = -0.3535533905932737*cM[14]; 
  data->AEM_S(0,35) = -0.3535533905932737*cM[15]; 
  data->AEM_S(0,36) = -0.3535533905932737*cM[16]; 
  data->AEM_S(0,37) = -0.3535533905932737*cM[17]; 
  data->AEM_S(0,38) = -0.3535533905932737*cM[18]; 
  data->AEM_S(0,39) = -0.3535533905932737*cM[19]; 
  data->AEM_S(1,20) = -0.3535533905932737*cM[1]; 
  data->AEM_S(1,21) = (-0.3162277660168379*cM[7])-0.3535533905932737*cM[0]; 
  data->AEM_S(1,22) = -0.3535533905932737*cM[4]; 
  data->AEM_S(1,23) = -0.3535533905932737*cM[5]; 
  data->AEM_S(1,24) = (-0.3162277660168379*cM[11])-0.3535533905932737*cM[2]; 
  data->AEM_S(1,25) = (-0.3162277660168379*cM[13])-0.3535533905932737*cM[3]; 
  data->AEM_S(1,26) = -0.3535533905932737*cM[10]; 
  data->AEM_S(1,27) = -0.3162277660168379*cM[1]; 
  data->AEM_S(1,28) = -0.3535533905932737*cM[12]; 
  data->AEM_S(1,29) = -0.3535533905932737*cM[15]; 
  data->AEM_S(1,30) = (-0.3162277660168379*cM[17])-0.3535533905932737*cM[6]; 
  data->AEM_S(1,31) = -0.3162277660168379*cM[4]; 
  data->AEM_S(1,32) = -0.3535533905932737*cM[8]; 
  data->AEM_S(1,33) = -0.3162277660168379*cM[5]; 
  data->AEM_S(1,34) = -0.3535533905932737*cM[18]; 
  data->AEM_S(1,35) = -0.3535533905932737*cM[9]; 
  data->AEM_S(1,36) = -0.3535533905932737*cM[19]; 
  data->AEM_S(1,37) = -0.3162277660168379*cM[10]; 
  data->AEM_S(1,38) = -0.3535533905932737*cM[14]; 
  data->AEM_S(1,39) = -0.3535533905932737*cM[16]; 
  data->AEM_S(2,20) = -0.3535533905932737*cM[2]; 
  data->AEM_S(2,21) = -0.3535533905932737*cM[4]; 
  data->AEM_S(2,22) = (-0.3162277660168379*cM[8])-0.3535533905932737*cM[0]; 
  data->AEM_S(2,23) = -0.3535533905932737*cM[6]; 
  data->AEM_S(2,24) = (-0.3162277660168379*cM[12])-0.3535533905932737*cM[1]; 
  data->AEM_S(2,25) = -0.3535533905932737*cM[10]; 
  data->AEM_S(2,26) = (-0.3162277660168379*cM[14])-0.3535533905932737*cM[3]; 
  data->AEM_S(2,27) = -0.3535533905932737*cM[11]; 
  data->AEM_S(2,28) = -0.3162277660168379*cM[2]; 
  data->AEM_S(2,29) = -0.3535533905932737*cM[16]; 
  data->AEM_S(2,30) = (-0.3162277660168379*cM[18])-0.3535533905932737*cM[5]; 
  data->AEM_S(2,31) = -0.3535533905932737*cM[7]; 
  data->AEM_S(2,32) = -0.3162277660168379*cM[4]; 
  data->AEM_S(2,33) = -0.3535533905932737*cM[17]; 
  data->AEM_S(2,34) = -0.3162277660168379*cM[6]; 
  data->AEM_S(2,35) = -0.3535533905932737*cM[19]; 
  data->AEM_S(2,36) = -0.3535533905932737*cM[9]; 
  data->AEM_S(2,37) = -0.3535533905932737*cM[13]; 
  data->AEM_S(2,38) = -0.3162277660168379*cM[10]; 
  data->AEM_S(2,39) = -0.3535533905932737*cM[15]; 
  data->AEM_S(3,20) = -0.3535533905932737*cM[3]; 
  data->AEM_S(3,21) = -0.3535533905932737*cM[5]; 
  data->AEM_S(3,22) = -0.3535533905932737*cM[6]; 
  data->AEM_S(3,23) = (-0.3162277660168379*cM[9])-0.3535533905932737*cM[0]; 
  data->AEM_S(3,24) = -0.3535533905932737*cM[10]; 
  data->AEM_S(3,25) = (-0.3162277660168379*cM[15])-0.3535533905932737*cM[1]; 
  data->AEM_S(3,26) = (-0.3162277660168379*cM[16])-0.3535533905932737*cM[2]; 
  data->AEM_S(3,27) = -0.3535533905932737*cM[13]; 
  data->AEM_S(3,28) = -0.3535533905932737*cM[14]; 
  data->AEM_S(3,29) = -0.3162277660168379*cM[3]; 
  data->AEM_S(3,30) = (-0.3162277660168379*cM[19])-0.3535533905932737*cM[4]; 
  data->AEM_S(3,31) = -0.3535533905932737*cM[17]; 
  data->AEM_S(3,32) = -0.3535533905932737*cM[18]; 
  data->AEM_S(3,33) = -0.3535533905932737*cM[7]; 
  data->AEM_S(3,34) = -0.3535533905932737*cM[8]; 
  data->AEM_S(3,35) = -0.3162277660168379*cM[5]; 
  data->AEM_S(3,36) = -0.3162277660168379*cM[6]; 
  data->AEM_S(3,37) = -0.3535533905932737*cM[11]; 
  data->AEM_S(3,38) = -0.3535533905932737*cM[12]; 
  data->AEM_S(3,39) = -0.3162277660168379*cM[10]; 
  data->AEM_S(4,20) = -0.3535533905932737*cM[4]; 
  data->AEM_S(4,21) = (-0.3162277660168379*cM[11])-0.3535533905932737*cM[2]; 
  data->AEM_S(4,22) = (-0.3162277660168379*cM[12])-0.3535533905932737*cM[1]; 
  data->AEM_S(4,23) = -0.3535533905932737*cM[10]; 
  data->AEM_S(4,24) = (-0.3162277660168379*cM[8])-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(4,25) = (-0.3162277660168379*cM[17])-0.3535533905932737*cM[6]; 
  data->AEM_S(4,26) = (-0.3162277660168379*cM[18])-0.3535533905932737*cM[5]; 
  data->AEM_S(4,27) = -0.3162277660168379*cM[4]; 
  data->AEM_S(4,28) = -0.3162277660168379*cM[4]; 
  data->AEM_S(4,29) = -0.3535533905932737*cM[19]; 
  data->AEM_S(4,30) = (-0.3162277660168379*cM[14])-0.3162277660168379*cM[13]-0.3535533905932737*cM[3]; 
  data->AEM_S(4,31) = (-0.2828427124746191*cM[12])-0.3162277660168379*cM[1]; 
  data->AEM_S(4,32) = (-0.2828427124746191*cM[11])-0.3162277660168379*cM[2]; 
  data->AEM_S(4,33) = -0.3162277660168379*cM[10]; 
  data->AEM_S(4,34) = -0.3162277660168379*cM[10]; 
  data->AEM_S(4,35) = -0.3535533905932737*cM[16]; 
  data->AEM_S(4,36) = -0.3535533905932737*cM[15]; 
  data->AEM_S(4,37) = (-0.2828427124746191*cM[18])-0.3162277660168379*cM[5]; 
  data->AEM_S(4,38) = (-0.2828427124746191*cM[17])-0.3162277660168379*cM[6]; 
  data->AEM_S(4,39) = -0.3535533905932737*cM[9]; 
  data->AEM_S(5,20) = -0.3535533905932737*cM[5]; 
  data->AEM_S(5,21) = (-0.3162277660168379*cM[13])-0.3535533905932737*cM[3]; 
  data->AEM_S(5,22) = -0.3535533905932737*cM[10]; 
  data->AEM_S(5,23) = (-0.3162277660168379*cM[15])-0.3535533905932737*cM[1]; 
  data->AEM_S(5,24) = (-0.3162277660168379*cM[17])-0.3535533905932737*cM[6]; 
  data->AEM_S(5,25) = (-0.3162277660168379*cM[9])-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(5,26) = (-0.3162277660168379*cM[19])-0.3535533905932737*cM[4]; 
  data->AEM_S(5,27) = -0.3162277660168379*cM[5]; 
  data->AEM_S(5,28) = -0.3535533905932737*cM[18]; 
  data->AEM_S(5,29) = -0.3162277660168379*cM[5]; 
  data->AEM_S(5,30) = (-0.3162277660168379*cM[16])-0.3162277660168379*cM[11]-0.3535533905932737*cM[2]; 
  data->AEM_S(5,31) = -0.3162277660168379*cM[10]; 
  data->AEM_S(5,32) = -0.3535533905932737*cM[14]; 
  data->AEM_S(5,33) = (-0.2828427124746191*cM[15])-0.3162277660168379*cM[1]; 
  data->AEM_S(5,34) = -0.3535533905932737*cM[12]; 
  data->AEM_S(5,35) = (-0.2828427124746191*cM[13])-0.3162277660168379*cM[3]; 
  data->AEM_S(5,36) = -0.3162277660168379*cM[10]; 
  data->AEM_S(5,37) = (-0.2828427124746191*cM[19])-0.3162277660168379*cM[4]; 
  data->AEM_S(5,38) = -0.3535533905932737*cM[8]; 
  data->AEM_S(5,39) = (-0.2828427124746191*cM[17])-0.3162277660168379*cM[6]; 
  data->AEM_S(6,20) = -0.3535533905932737*cM[6]; 
  data->AEM_S(6,21) = -0.3535533905932737*cM[10]; 
  data->AEM_S(6,22) = (-0.3162277660168379*cM[14])-0.3535533905932737*cM[3]; 
  data->AEM_S(6,23) = (-0.3162277660168379*cM[16])-0.3535533905932737*cM[2]; 
  data->AEM_S(6,24) = (-0.3162277660168379*cM[18])-0.3535533905932737*cM[5]; 
  data->AEM_S(6,25) = (-0.3162277660168379*cM[19])-0.3535533905932737*cM[4]; 
  data->AEM_S(6,26) = (-0.3162277660168379*cM[9])-0.3162277660168379*cM[8]-0.3535533905932737*cM[0]; 
  data->AEM_S(6,27) = -0.3535533905932737*cM[17]; 
  data->AEM_S(6,28) = -0.3162277660168379*cM[6]; 
  data->AEM_S(6,29) = -0.3162277660168379*cM[6]; 
  data->AEM_S(6,30) = (-0.3162277660168379*cM[15])-0.3162277660168379*cM[12]-0.3535533905932737*cM[1]; 
  data->AEM_S(6,31) = -0.3535533905932737*cM[13]; 
  data->AEM_S(6,32) = -0.3162277660168379*cM[10]; 
  data->AEM_S(6,33) = -0.3535533905932737*cM[11]; 
  data->AEM_S(6,34) = (-0.2828427124746191*cM[16])-0.3162277660168379*cM[2]; 
  data->AEM_S(6,35) = -0.3162277660168379*cM[10]; 
  data->AEM_S(6,36) = (-0.2828427124746191*cM[14])-0.3162277660168379*cM[3]; 
  data->AEM_S(6,37) = -0.3535533905932737*cM[7]; 
  data->AEM_S(6,38) = (-0.2828427124746191*cM[19])-0.3162277660168379*cM[4]; 
  data->AEM_S(6,39) = (-0.2828427124746191*cM[18])-0.3162277660168379*cM[5]; 
  data->AEM_S(7,20) = -0.3535533905932737*cM[7]; 
  data->AEM_S(7,21) = -0.3162277660168379*cM[1]; 
  data->AEM_S(7,22) = -0.3535533905932737*cM[11]; 
  data->AEM_S(7,23) = -0.3535533905932737*cM[13]; 
  data->AEM_S(7,24) = -0.3162277660168379*cM[4]; 
  data->AEM_S(7,25) = -0.3162277660168379*cM[5]; 
  data->AEM_S(7,26) = -0.3535533905932737*cM[17]; 
  data->AEM_S(7,27) = (-0.2258769757263128*cM[7])-0.3535533905932737*cM[0]; 
  data->AEM_S(7,30) = -0.3162277660168379*cM[10]; 
  data->AEM_S(7,31) = (-0.2258769757263128*cM[11])-0.3535533905932737*cM[2]; 
  data->AEM_S(7,32) = -0.3162277660168379*cM[12]; 
  data->AEM_S(7,33) = (-0.2258769757263128*cM[13])-0.3535533905932737*cM[3]; 
  data->AEM_S(7,35) = -0.3162277660168379*cM[15]; 
  data->AEM_S(7,37) = (-0.2258769757263128*cM[17])-0.3535533905932737*cM[6]; 
  data->AEM_S(7,38) = -0.3162277660168379*cM[18]; 
  data->AEM_S(7,39) = -0.3162277660168379*cM[19]; 
  data->AEM_S(8,20) = -0.3535533905932737*cM[8]; 
  data->AEM_S(8,21) = -0.3535533905932737*cM[12]; 
  data->AEM_S(8,22) = -0.3162277660168379*cM[2]; 
  data->AEM_S(8,23) = -0.3535533905932737*cM[14]; 
  data->AEM_S(8,24) = -0.3162277660168379*cM[4]; 
  data->AEM_S(8,25) = -0.3535533905932737*cM[18]; 
  data->AEM_S(8,26) = -0.3162277660168379*cM[6]; 
  data->AEM_S(8,28) = (-0.2258769757263128*cM[8])-0.3535533905932737*cM[0]; 
  data->AEM_S(8,30) = -0.3162277660168379*cM[10]; 
  data->AEM_S(8,31) = -0.3162277660168379*cM[11]; 
  data->AEM_S(8,32) = (-0.2258769757263128*cM[12])-0.3535533905932737*cM[1]; 
  data->AEM_S(8,34) = (-0.2258769757263128*cM[14])-0.3535533905932737*cM[3]; 
  data->AEM_S(8,36) = -0.3162277660168379*cM[16]; 
  data->AEM_S(8,37) = -0.3162277660168379*cM[17]; 
  data->AEM_S(8,38) = (-0.2258769757263128*cM[18])-0.3535533905932737*cM[5]; 
  data->AEM_S(8,39) = -0.3162277660168379*cM[19]; 
  data->AEM_S(9,20) = -0.3535533905932737*cM[9]; 
  data->AEM_S(9,21) = -0.3535533905932737*cM[15]; 
  data->AEM_S(9,22) = -0.3535533905932737*cM[16]; 
  data->AEM_S(9,23) = -0.3162277660168379*cM[3]; 
  data->AEM_S(9,24) = -0.3535533905932737*cM[19]; 
  data->AEM_S(9,25) = -0.3162277660168379*cM[5]; 
  data->AEM_S(9,26) = -0.3162277660168379*cM[6]; 
  data->AEM_S(9,29) = (-0.2258769757263128*cM[9])-0.3535533905932737*cM[0]; 
  data->AEM_S(9,30) = -0.3162277660168379*cM[10]; 
  data->AEM_S(9,33) = -0.3162277660168379*cM[13]; 
  data->AEM_S(9,34) = -0.3162277660168379*cM[14]; 
  data->AEM_S(9,35) = (-0.2258769757263128*cM[15])-0.3535533905932737*cM[1]; 
  data->AEM_S(9,36) = (-0.2258769757263128*cM[16])-0.3535533905932737*cM[2]; 
  data->AEM_S(9,37) = -0.3162277660168379*cM[17]; 
  data->AEM_S(9,38) = -0.3162277660168379*cM[18]; 
  data->AEM_S(9,39) = (-0.2258769757263128*cM[19])-0.3535533905932737*cM[4]; 
  data->AEM_S(10,20) = -0.3535533905932737*cM[10]; 
  data->AEM_S(10,21) = (-0.3162277660168379*cM[17])-0.3535533905932737*cM[6]; 
  data->AEM_S(10,22) = (-0.3162277660168379*cM[18])-0.3535533905932737*cM[5]; 
  data->AEM_S(10,23) = (-0.3162277660168379*cM[19])-0.3535533905932737*cM[4]; 
  data->AEM_S(10,24) = (-0.3162277660168379*cM[14])-0.3162277660168379*cM[13]-0.3535533905932737*cM[3]; 
  data->AEM_S(10,25) = (-0.3162277660168379*cM[16])-0.3162277660168379*cM[11]-0.3535533905932737*cM[2]; 
  data->AEM_S(10,26) = (-0.3162277660168379*cM[15])-0.3162277660168379*cM[12]-0.3535533905932737*cM[1]; 
  data->AEM_S(10,27) = -0.3162277660168379*cM[10]; 
  data->AEM_S(10,28) = -0.3162277660168379*cM[10]; 
  data->AEM_S(10,29) = -0.3162277660168379*cM[10]; 
  data->AEM_S(10,30) = (-0.3162277660168379*cM[9])-0.3162277660168379*cM[8]-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(10,31) = (-0.282842712474619*cM[18])-0.3162277660168379*cM[5]; 
  data->AEM_S(10,32) = (-0.282842712474619*cM[17])-0.3162277660168379*cM[6]; 
  data->AEM_S(10,33) = (-0.282842712474619*cM[19])-0.3162277660168379*cM[4]; 
  data->AEM_S(10,34) = (-0.282842712474619*cM[19])-0.3162277660168379*cM[4]; 
  data->AEM_S(10,35) = (-0.282842712474619*cM[17])-0.3162277660168379*cM[6]; 
  data->AEM_S(10,36) = (-0.282842712474619*cM[18])-0.3162277660168379*cM[5]; 
  data->AEM_S(10,37) = (-0.282842712474619*cM[15])-0.282842712474619*cM[12]-0.3162277660168379*cM[1]; 
  data->AEM_S(10,38) = (-0.282842712474619*cM[16])-0.282842712474619*cM[11]-0.3162277660168379*cM[2]; 
  data->AEM_S(10,39) = (-0.282842712474619*cM[14])-0.282842712474619*cM[13]-0.3162277660168379*cM[3]; 
  data->AEM_S(11,20) = -0.3535533905932737*cM[11]; 
  data->AEM_S(11,21) = -0.3162277660168379*cM[4]; 
  data->AEM_S(11,22) = -0.3535533905932737*cM[7]; 
  data->AEM_S(11,23) = -0.3535533905932737*cM[17]; 
  data->AEM_S(11,24) = (-0.2828427124746191*cM[12])-0.3162277660168379*cM[1]; 
  data->AEM_S(11,25) = -0.3162277660168379*cM[10]; 
  data->AEM_S(11,26) = -0.3535533905932737*cM[13]; 
  data->AEM_S(11,27) = (-0.2258769757263128*cM[11])-0.3535533905932737*cM[2]; 
  data->AEM_S(11,28) = -0.3162277660168379*cM[11]; 
  data->AEM_S(11,30) = (-0.282842712474619*cM[18])-0.3162277660168379*cM[5]; 
  data->AEM_S(11,31) = (-0.3162277660168379*cM[8])-0.2258769757263128*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(11,32) = -0.2828427124746191*cM[4]; 
  data->AEM_S(11,33) = (-0.2258769757263128*cM[17])-0.3535533905932737*cM[6]; 
  data->AEM_S(11,34) = -0.3162277660168379*cM[17]; 
  data->AEM_S(11,35) = -0.3162277660168379*cM[19]; 
  data->AEM_S(11,37) = (-0.3162277660168379*cM[14])-0.2258769757263128*cM[13]-0.3535533905932737*cM[3]; 
  data->AEM_S(11,38) = -0.282842712474619*cM[10]; 
  data->AEM_S(11,39) = -0.3162277660168379*cM[15]; 
  data->AEM_S(12,20) = -0.3535533905932737*cM[12]; 
  data->AEM_S(12,21) = -0.3535533905932737*cM[8]; 
  data->AEM_S(12,22) = -0.3162277660168379*cM[4]; 
  data->AEM_S(12,23) = -0.3535533905932737*cM[18]; 
  data->AEM_S(12,24) = (-0.2828427124746191*cM[11])-0.3162277660168379*cM[2]; 
  data->AEM_S(12,25) = -0.3535533905932737*cM[14]; 
  data->AEM_S(12,26) = -0.3162277660168379*cM[10]; 
  data->AEM_S(12,27) = -0.3162277660168379*cM[12]; 
  data->AEM_S(12,28) = (-0.2258769757263128*cM[12])-0.3535533905932737*cM[1]; 
  data->AEM_S(12,30) = (-0.282842712474619*cM[17])-0.3162277660168379*cM[6]; 
  data->AEM_S(12,31) = -0.2828427124746191*cM[4]; 
  data->AEM_S(12,32) = (-0.2258769757263128*cM[8])-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(12,33) = -0.3162277660168379*cM[18]; 
  data->AEM_S(12,34) = (-0.2258769757263128*cM[18])-0.3535533905932737*cM[5]; 
  data->AEM_S(12,36) = -0.3162277660168379*cM[19]; 
  data->AEM_S(12,37) = -0.282842712474619*cM[10]; 
  data->AEM_S(12,38) = (-0.2258769757263128*cM[14])-0.3162277660168379*cM[13]-0.3535533905932737*cM[3]; 
  data->AEM_S(12,39) = -0.3162277660168379*cM[16]; 
  data->AEM_S(13,20) = -0.3535533905932737*cM[13]; 
  data->AEM_S(13,21) = -0.3162277660168379*cM[5]; 
  data->AEM_S(13,22) = -0.3535533905932737*cM[17]; 
  data->AEM_S(13,23) = -0.3535533905932737*cM[7]; 
  data->AEM_S(13,24) = -0.3162277660168379*cM[10]; 
  data->AEM_S(13,25) = (-0.2828427124746191*cM[15])-0.3162277660168379*cM[1]; 
  data->AEM_S(13,26) = -0.3535533905932737*cM[11]; 
  data->AEM_S(13,27) = (-0.2258769757263128*cM[13])-0.3535533905932737*cM[3]; 
  data->AEM_S(13,29) = -0.3162277660168379*cM[13]; 
  data->AEM_S(13,30) = (-0.282842712474619*cM[19])-0.3162277660168379*cM[4]; 
  data->AEM_S(13,31) = (-0.2258769757263128*cM[17])-0.3535533905932737*cM[6]; 
  data->AEM_S(13,32) = -0.3162277660168379*cM[18]; 
  data->AEM_S(13,33) = (-0.3162277660168379*cM[9])-0.2258769757263128*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(13,35) = -0.2828427124746191*cM[5]; 
  data->AEM_S(13,36) = -0.3162277660168379*cM[17]; 
  data->AEM_S(13,37) = (-0.3162277660168379*cM[16])-0.2258769757263128*cM[11]-0.3535533905932737*cM[2]; 
  data->AEM_S(13,38) = -0.3162277660168379*cM[12]; 
  data->AEM_S(13,39) = -0.282842712474619*cM[10]; 
  data->AEM_S(14,20) = -0.3535533905932737*cM[14]; 
  data->AEM_S(14,21) = -0.3535533905932737*cM[18]; 
  data->AEM_S(14,22) = -0.3162277660168379*cM[6]; 
  data->AEM_S(14,23) = -0.3535533905932737*cM[8]; 
  data->AEM_S(14,24) = -0.3162277660168379*cM[10]; 
  data->AEM_S(14,25) = -0.3535533905932737*cM[12]; 
  data->AEM_S(14,26) = (-0.2828427124746191*cM[16])-0.3162277660168379*cM[2]; 
  data->AEM_S(14,28) = (-0.2258769757263128*cM[14])-0.3535533905932737*cM[3]; 
  data->AEM_S(14,29) = -0.3162277660168379*cM[14]; 
  data->AEM_S(14,30) = (-0.282842712474619*cM[19])-0.3162277660168379*cM[4]; 
  data->AEM_S(14,31) = -0.3162277660168379*cM[17]; 
  data->AEM_S(14,32) = (-0.2258769757263128*cM[18])-0.3535533905932737*cM[5]; 
  data->AEM_S(14,34) = (-0.3162277660168379*cM[9])-0.2258769757263128*cM[8]-0.3535533905932737*cM[0]; 
  data->AEM_S(14,35) = -0.3162277660168379*cM[18]; 
  data->AEM_S(14,36) = -0.2828427124746191*cM[6]; 
  data->AEM_S(14,37) = -0.3162277660168379*cM[11]; 
  data->AEM_S(14,38) = (-0.3162277660168379*cM[15])-0.2258769757263128*cM[12]-0.3535533905932737*cM[1]; 
  data->AEM_S(14,39) = -0.282842712474619*cM[10]; 
  data->AEM_S(15,20) = -0.3535533905932737*cM[15]; 
  data->AEM_S(15,21) = -0.3535533905932737*cM[9]; 
  data->AEM_S(15,22) = -0.3535533905932737*cM[19]; 
  data->AEM_S(15,23) = -0.3162277660168379*cM[5]; 
  data->AEM_S(15,24) = -0.3535533905932737*cM[16]; 
  data->AEM_S(15,25) = (-0.2828427124746191*cM[13])-0.3162277660168379*cM[3]; 
  data->AEM_S(15,26) = -0.3162277660168379*cM[10]; 
  data->AEM_S(15,27) = -0.3162277660168379*cM[15]; 
  data->AEM_S(15,29) = (-0.2258769757263128*cM[15])-0.3535533905932737*cM[1]; 
  data->AEM_S(15,30) = (-0.282842712474619*cM[17])-0.3162277660168379*cM[6]; 
  data->AEM_S(15,31) = -0.3162277660168379*cM[19]; 
  data->AEM_S(15,33) = -0.2828427124746191*cM[5]; 
  data->AEM_S(15,34) = -0.3162277660168379*cM[18]; 
  data->AEM_S(15,35) = (-0.2258769757263128*cM[9])-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(15,36) = (-0.2258769757263128*cM[19])-0.3535533905932737*cM[4]; 
  data->AEM_S(15,37) = -0.282842712474619*cM[10]; 
  data->AEM_S(15,38) = -0.3162277660168379*cM[14]; 
  data->AEM_S(15,39) = (-0.2258769757263128*cM[16])-0.3162277660168379*cM[11]-0.3535533905932737*cM[2]; 
  data->AEM_S(16,20) = -0.3535533905932737*cM[16]; 
  data->AEM_S(16,21) = -0.3535533905932737*cM[19]; 
  data->AEM_S(16,22) = -0.3535533905932737*cM[9]; 
  data->AEM_S(16,23) = -0.3162277660168379*cM[6]; 
  data->AEM_S(16,24) = -0.3535533905932737*cM[15]; 
  data->AEM_S(16,25) = -0.3162277660168379*cM[10]; 
  data->AEM_S(16,26) = (-0.2828427124746191*cM[14])-0.3162277660168379*cM[3]; 
  data->AEM_S(16,28) = -0.3162277660168379*cM[16]; 
  data->AEM_S(16,29) = (-0.2258769757263128*cM[16])-0.3535533905932737*cM[2]; 
  data->AEM_S(16,30) = (-0.282842712474619*cM[18])-0.3162277660168379*cM[5]; 
  data->AEM_S(16,32) = -0.3162277660168379*cM[19]; 
  data->AEM_S(16,33) = -0.3162277660168379*cM[17]; 
  data->AEM_S(16,34) = -0.2828427124746191*cM[6]; 
  data->AEM_S(16,35) = (-0.2258769757263128*cM[19])-0.3535533905932737*cM[4]; 
  data->AEM_S(16,36) = (-0.2258769757263128*cM[9])-0.3162277660168379*cM[8]-0.3535533905932737*cM[0]; 
  data->AEM_S(16,37) = -0.3162277660168379*cM[13]; 
  data->AEM_S(16,38) = -0.282842712474619*cM[10]; 
  data->AEM_S(16,39) = (-0.2258769757263128*cM[15])-0.3162277660168379*cM[12]-0.3535533905932737*cM[1]; 
  data->AEM_S(17,20) = -0.3535533905932737*cM[17]; 
  data->AEM_S(17,21) = -0.3162277660168379*cM[10]; 
  data->AEM_S(17,22) = -0.3535533905932737*cM[13]; 
  data->AEM_S(17,23) = -0.3535533905932737*cM[11]; 
  data->AEM_S(17,24) = (-0.2828427124746191*cM[18])-0.3162277660168379*cM[5]; 
  data->AEM_S(17,25) = (-0.2828427124746191*cM[19])-0.3162277660168379*cM[4]; 
  data->AEM_S(17,26) = -0.3535533905932737*cM[7]; 
  data->AEM_S(17,27) = (-0.2258769757263128*cM[17])-0.3535533905932737*cM[6]; 
  data->AEM_S(17,28) = -0.3162277660168379*cM[17]; 
  data->AEM_S(17,29) = -0.3162277660168379*cM[17]; 
  data->AEM_S(17,30) = (-0.282842712474619*cM[15])-0.282842712474619*cM[12]-0.3162277660168379*cM[1]; 
  data->AEM_S(17,31) = (-0.3162277660168379*cM[14])-0.2258769757263128*cM[13]-0.3535533905932737*cM[3]; 
  data->AEM_S(17,32) = -0.282842712474619*cM[10]; 
  data->AEM_S(17,33) = (-0.3162277660168379*cM[16])-0.2258769757263128*cM[11]-0.3535533905932737*cM[2]; 
  data->AEM_S(17,34) = -0.3162277660168379*cM[11]; 
  data->AEM_S(17,35) = -0.282842712474619*cM[10]; 
  data->AEM_S(17,36) = -0.3162277660168379*cM[13]; 
  data->AEM_S(17,37) = (-0.3162277660168379*cM[9])-0.3162277660168379*cM[8]-0.2258769757263128*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(17,38) = (-0.2529822128134704*cM[19])-0.2828427124746191*cM[4]; 
  data->AEM_S(17,39) = (-0.2529822128134704*cM[18])-0.2828427124746191*cM[5]; 
  data->AEM_S(18,20) = -0.3535533905932737*cM[18]; 
  data->AEM_S(18,21) = -0.3535533905932737*cM[14]; 
  data->AEM_S(18,22) = -0.3162277660168379*cM[10]; 
  data->AEM_S(18,23) = -0.3535533905932737*cM[12]; 
  data->AEM_S(18,24) = (-0.2828427124746191*cM[17])-0.3162277660168379*cM[6]; 
  data->AEM_S(18,25) = -0.3535533905932737*cM[8]; 
  data->AEM_S(18,26) = (-0.2828427124746191*cM[19])-0.3162277660168379*cM[4]; 
  data->AEM_S(18,27) = -0.3162277660168379*cM[18]; 
  data->AEM_S(18,28) = (-0.2258769757263128*cM[18])-0.3535533905932737*cM[5]; 
  data->AEM_S(18,29) = -0.3162277660168379*cM[18]; 
  data->AEM_S(18,30) = (-0.282842712474619*cM[16])-0.282842712474619*cM[11]-0.3162277660168379*cM[2]; 
  data->AEM_S(18,31) = -0.282842712474619*cM[10]; 
  data->AEM_S(18,32) = (-0.2258769757263128*cM[14])-0.3162277660168379*cM[13]-0.3535533905932737*cM[3]; 
  data->AEM_S(18,33) = -0.3162277660168379*cM[12]; 
  data->AEM_S(18,34) = (-0.3162277660168379*cM[15])-0.2258769757263128*cM[12]-0.3535533905932737*cM[1]; 
  data->AEM_S(18,35) = -0.3162277660168379*cM[14]; 
  data->AEM_S(18,36) = -0.282842712474619*cM[10]; 
  data->AEM_S(18,37) = (-0.2529822128134704*cM[19])-0.2828427124746191*cM[4]; 
  data->AEM_S(18,38) = (-0.3162277660168379*cM[9])-0.2258769757263128*cM[8]-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  data->AEM_S(18,39) = (-0.2529822128134704*cM[17])-0.2828427124746191*cM[6]; 
  data->AEM_S(19,20) = -0.3535533905932737*cM[19]; 
  data->AEM_S(19,21) = -0.3535533905932737*cM[16]; 
  data->AEM_S(19,22) = -0.3535533905932737*cM[15]; 
  data->AEM_S(19,23) = -0.3162277660168379*cM[10]; 
  data->AEM_S(19,24) = -0.3535533905932737*cM[9]; 
  data->AEM_S(19,25) = (-0.2828427124746191*cM[17])-0.3162277660168379*cM[6]; 
  data->AEM_S(19,26) = (-0.2828427124746191*cM[18])-0.3162277660168379*cM[5]; 
  data->AEM_S(19,27) = -0.3162277660168379*cM[19]; 
  data->AEM_S(19,28) = -0.3162277660168379*cM[19]; 
  data->AEM_S(19,29) = (-0.2258769757263128*cM[19])-0.3535533905932737*cM[4]; 
  data->AEM_S(19,30) = (-0.282842712474619*cM[14])-0.282842712474619*cM[13]-0.3162277660168379*cM[3]; 
  data->AEM_S(19,31) = -0.3162277660168379*cM[15]; 
  data->AEM_S(19,32) = -0.3162277660168379*cM[16]; 
  data->AEM_S(19,33) = -0.282842712474619*cM[10]; 
  data->AEM_S(19,34) = -0.282842712474619*cM[10]; 
  data->AEM_S(19,35) = (-0.2258769757263128*cM[16])-0.3162277660168379*cM[11]-0.3535533905932737*cM[2]; 
  data->AEM_S(19,36) = (-0.2258769757263128*cM[15])-0.3162277660168379*cM[12]-0.3535533905932737*cM[1]; 
  data->AEM_S(19,37) = (-0.2529822128134704*cM[18])-0.2828427124746191*cM[5]; 
  data->AEM_S(19,38) = (-0.2529822128134704*cM[17])-0.2828427124746191*cM[6]; 
  data->AEM_S(19,39) = (-0.2258769757263128*cM[9])-0.3162277660168379*cM[8]-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(20,0) = 0.3535533905932737*m1r[0]; 
  data->AEM_S(20,1) = 0.3535533905932737*m1r[1]; 
  data->AEM_S(20,2) = 0.3535533905932737*m1r[2]; 
  data->AEM_S(20,3) = 0.3535533905932737*m1r[3]; 
  data->AEM_S(20,4) = 0.3535533905932737*m1r[4]; 
  data->AEM_S(20,5) = 0.3535533905932737*m1r[5]; 
  data->AEM_S(20,6) = 0.3535533905932737*m1r[6]; 
  data->AEM_S(20,7) = 0.3535533905932737*m1r[7]; 
  data->AEM_S(20,8) = 0.3535533905932737*m1r[8]; 
  data->AEM_S(20,9) = 0.3535533905932737*m1r[9]; 
  data->AEM_S(20,10) = 0.3535533905932737*m1r[10]; 
  data->AEM_S(20,11) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(20,12) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(20,13) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(20,14) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(20,15) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(20,16) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(20,17) = 0.3535533905932737*m1r[17]; 
  data->AEM_S(20,18) = 0.3535533905932737*m1r[18]; 
  data->AEM_S(20,19) = 0.3535533905932737*m1r[19]; 
  data->AEM_S(21,0) = 0.3535533905932737*m1r[1]; 
  data->AEM_S(21,1) = 0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(21,2) = 0.3535533905932737*m1r[4]; 
  data->AEM_S(21,3) = 0.3535533905932737*m1r[5]; 
  data->AEM_S(21,4) = 0.3162277660168379*m1r[11]+0.3535533905932737*m1r[2]; 
  data->AEM_S(21,5) = 0.3162277660168379*m1r[13]+0.3535533905932737*m1r[3]; 
  data->AEM_S(21,6) = 0.3535533905932737*m1r[10]; 
  data->AEM_S(21,7) = 0.3162277660168379*m1r[1]; 
  data->AEM_S(21,8) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(21,9) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(21,10) = 0.3162277660168379*m1r[17]+0.3535533905932737*m1r[6]; 
  data->AEM_S(21,11) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(21,12) = 0.3535533905932737*m1r[8]; 
  data->AEM_S(21,13) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(21,14) = 0.3535533905932737*m1r[18]; 
  data->AEM_S(21,15) = 0.3535533905932737*m1r[9]; 
  data->AEM_S(21,16) = 0.3535533905932737*m1r[19]; 
  data->AEM_S(21,17) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(21,18) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(21,19) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(22,0) = 0.3535533905932737*m1r[2]; 
  data->AEM_S(22,1) = 0.3535533905932737*m1r[4]; 
  data->AEM_S(22,2) = 0.3162277660168379*m1r[8]+0.3535533905932737*m1r[0]; 
  data->AEM_S(22,3) = 0.3535533905932737*m1r[6]; 
  data->AEM_S(22,4) = 0.3162277660168379*m1r[12]+0.3535533905932737*m1r[1]; 
  data->AEM_S(22,5) = 0.3535533905932737*m1r[10]; 
  data->AEM_S(22,6) = 0.3162277660168379*m1r[14]+0.3535533905932737*m1r[3]; 
  data->AEM_S(22,7) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(22,8) = 0.3162277660168379*m1r[2]; 
  data->AEM_S(22,9) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(22,10) = 0.3162277660168379*m1r[18]+0.3535533905932737*m1r[5]; 
  data->AEM_S(22,11) = 0.3535533905932737*m1r[7]; 
  data->AEM_S(22,12) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(22,13) = 0.3535533905932737*m1r[17]; 
  data->AEM_S(22,14) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(22,15) = 0.3535533905932737*m1r[19]; 
  data->AEM_S(22,16) = 0.3535533905932737*m1r[9]; 
  data->AEM_S(22,17) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(22,18) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(22,19) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(23,0) = 0.3535533905932737*m1r[3]; 
  data->AEM_S(23,1) = 0.3535533905932737*m1r[5]; 
  data->AEM_S(23,2) = 0.3535533905932737*m1r[6]; 
  data->AEM_S(23,3) = 0.3162277660168379*m1r[9]+0.3535533905932737*m1r[0]; 
  data->AEM_S(23,4) = 0.3535533905932737*m1r[10]; 
  data->AEM_S(23,5) = 0.3162277660168379*m1r[15]+0.3535533905932737*m1r[1]; 
  data->AEM_S(23,6) = 0.3162277660168379*m1r[16]+0.3535533905932737*m1r[2]; 
  data->AEM_S(23,7) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(23,8) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(23,9) = 0.3162277660168379*m1r[3]; 
  data->AEM_S(23,10) = 0.3162277660168379*m1r[19]+0.3535533905932737*m1r[4]; 
  data->AEM_S(23,11) = 0.3535533905932737*m1r[17]; 
  data->AEM_S(23,12) = 0.3535533905932737*m1r[18]; 
  data->AEM_S(23,13) = 0.3535533905932737*m1r[7]; 
  data->AEM_S(23,14) = 0.3535533905932737*m1r[8]; 
  data->AEM_S(23,15) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(23,16) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(23,17) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(23,18) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(23,19) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(24,0) = 0.3535533905932737*m1r[4]; 
  data->AEM_S(24,1) = 0.3162277660168379*m1r[11]+0.3535533905932737*m1r[2]; 
  data->AEM_S(24,2) = 0.3162277660168379*m1r[12]+0.3535533905932737*m1r[1]; 
  data->AEM_S(24,3) = 0.3535533905932737*m1r[10]; 
  data->AEM_S(24,4) = 0.3162277660168379*m1r[8]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(24,5) = 0.3162277660168379*m1r[17]+0.3535533905932737*m1r[6]; 
  data->AEM_S(24,6) = 0.3162277660168379*m1r[18]+0.3535533905932737*m1r[5]; 
  data->AEM_S(24,7) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(24,8) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(24,9) = 0.3535533905932737*m1r[19]; 
  data->AEM_S(24,10) = 0.3162277660168379*m1r[14]+0.3162277660168379*m1r[13]+0.3535533905932737*m1r[3]; 
  data->AEM_S(24,11) = 0.2828427124746191*m1r[12]+0.3162277660168379*m1r[1]; 
  data->AEM_S(24,12) = 0.2828427124746191*m1r[11]+0.3162277660168379*m1r[2]; 
  data->AEM_S(24,13) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(24,14) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(24,15) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(24,16) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(24,17) = 0.2828427124746191*m1r[18]+0.3162277660168379*m1r[5]; 
  data->AEM_S(24,18) = 0.2828427124746191*m1r[17]+0.3162277660168379*m1r[6]; 
  data->AEM_S(24,19) = 0.3535533905932737*m1r[9]; 
  data->AEM_S(25,0) = 0.3535533905932737*m1r[5]; 
  data->AEM_S(25,1) = 0.3162277660168379*m1r[13]+0.3535533905932737*m1r[3]; 
  data->AEM_S(25,2) = 0.3535533905932737*m1r[10]; 
  data->AEM_S(25,3) = 0.3162277660168379*m1r[15]+0.3535533905932737*m1r[1]; 
  data->AEM_S(25,4) = 0.3162277660168379*m1r[17]+0.3535533905932737*m1r[6]; 
  data->AEM_S(25,5) = 0.3162277660168379*m1r[9]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(25,6) = 0.3162277660168379*m1r[19]+0.3535533905932737*m1r[4]; 
  data->AEM_S(25,7) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(25,8) = 0.3535533905932737*m1r[18]; 
  data->AEM_S(25,9) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(25,10) = 0.3162277660168379*m1r[16]+0.3162277660168379*m1r[11]+0.3535533905932737*m1r[2]; 
  data->AEM_S(25,11) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(25,12) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(25,13) = 0.2828427124746191*m1r[15]+0.3162277660168379*m1r[1]; 
  data->AEM_S(25,14) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(25,15) = 0.2828427124746191*m1r[13]+0.3162277660168379*m1r[3]; 
  data->AEM_S(25,16) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(25,17) = 0.2828427124746191*m1r[19]+0.3162277660168379*m1r[4]; 
  data->AEM_S(25,18) = 0.3535533905932737*m1r[8]; 
  data->AEM_S(25,19) = 0.2828427124746191*m1r[17]+0.3162277660168379*m1r[6]; 
  data->AEM_S(26,0) = 0.3535533905932737*m1r[6]; 
  data->AEM_S(26,1) = 0.3535533905932737*m1r[10]; 
  data->AEM_S(26,2) = 0.3162277660168379*m1r[14]+0.3535533905932737*m1r[3]; 
  data->AEM_S(26,3) = 0.3162277660168379*m1r[16]+0.3535533905932737*m1r[2]; 
  data->AEM_S(26,4) = 0.3162277660168379*m1r[18]+0.3535533905932737*m1r[5]; 
  data->AEM_S(26,5) = 0.3162277660168379*m1r[19]+0.3535533905932737*m1r[4]; 
  data->AEM_S(26,6) = 0.3162277660168379*m1r[9]+0.3162277660168379*m1r[8]+0.3535533905932737*m1r[0]; 
  data->AEM_S(26,7) = 0.3535533905932737*m1r[17]; 
  data->AEM_S(26,8) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(26,9) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(26,10) = 0.3162277660168379*m1r[15]+0.3162277660168379*m1r[12]+0.3535533905932737*m1r[1]; 
  data->AEM_S(26,11) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(26,12) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(26,13) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(26,14) = 0.2828427124746191*m1r[16]+0.3162277660168379*m1r[2]; 
  data->AEM_S(26,15) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(26,16) = 0.2828427124746191*m1r[14]+0.3162277660168379*m1r[3]; 
  data->AEM_S(26,17) = 0.3535533905932737*m1r[7]; 
  data->AEM_S(26,18) = 0.2828427124746191*m1r[19]+0.3162277660168379*m1r[4]; 
  data->AEM_S(26,19) = 0.2828427124746191*m1r[18]+0.3162277660168379*m1r[5]; 
  data->AEM_S(27,0) = 0.3535533905932737*m1r[7]; 
  data->AEM_S(27,1) = 0.3162277660168379*m1r[1]; 
  data->AEM_S(27,2) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(27,3) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(27,4) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(27,5) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(27,6) = 0.3535533905932737*m1r[17]; 
  data->AEM_S(27,7) = 0.2258769757263128*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(27,10) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(27,11) = 0.2258769757263128*m1r[11]+0.3535533905932737*m1r[2]; 
  data->AEM_S(27,12) = 0.3162277660168379*m1r[12]; 
  data->AEM_S(27,13) = 0.2258769757263128*m1r[13]+0.3535533905932737*m1r[3]; 
  data->AEM_S(27,15) = 0.3162277660168379*m1r[15]; 
  data->AEM_S(27,17) = 0.2258769757263128*m1r[17]+0.3535533905932737*m1r[6]; 
  data->AEM_S(27,18) = 0.3162277660168379*m1r[18]; 
  data->AEM_S(27,19) = 0.3162277660168379*m1r[19]; 
  data->AEM_S(28,0) = 0.3535533905932737*m1r[8]; 
  data->AEM_S(28,1) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(28,2) = 0.3162277660168379*m1r[2]; 
  data->AEM_S(28,3) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(28,4) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(28,5) = 0.3535533905932737*m1r[18]; 
  data->AEM_S(28,6) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(28,8) = 0.2258769757263128*m1r[8]+0.3535533905932737*m1r[0]; 
  data->AEM_S(28,10) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(28,11) = 0.3162277660168379*m1r[11]; 
  data->AEM_S(28,12) = 0.2258769757263128*m1r[12]+0.3535533905932737*m1r[1]; 
  data->AEM_S(28,14) = 0.2258769757263128*m1r[14]+0.3535533905932737*m1r[3]; 
  data->AEM_S(28,16) = 0.3162277660168379*m1r[16]; 
  data->AEM_S(28,17) = 0.3162277660168379*m1r[17]; 
  data->AEM_S(28,18) = 0.2258769757263128*m1r[18]+0.3535533905932737*m1r[5]; 
  data->AEM_S(28,19) = 0.3162277660168379*m1r[19]; 
  data->AEM_S(29,0) = 0.3535533905932737*m1r[9]; 
  data->AEM_S(29,1) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(29,2) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(29,3) = 0.3162277660168379*m1r[3]; 
  data->AEM_S(29,4) = 0.3535533905932737*m1r[19]; 
  data->AEM_S(29,5) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(29,6) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(29,9) = 0.2258769757263128*m1r[9]+0.3535533905932737*m1r[0]; 
  data->AEM_S(29,10) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(29,13) = 0.3162277660168379*m1r[13]; 
  data->AEM_S(29,14) = 0.3162277660168379*m1r[14]; 
  data->AEM_S(29,15) = 0.2258769757263128*m1r[15]+0.3535533905932737*m1r[1]; 
  data->AEM_S(29,16) = 0.2258769757263128*m1r[16]+0.3535533905932737*m1r[2]; 
  data->AEM_S(29,17) = 0.3162277660168379*m1r[17]; 
  data->AEM_S(29,18) = 0.3162277660168379*m1r[18]; 
  data->AEM_S(29,19) = 0.2258769757263128*m1r[19]+0.3535533905932737*m1r[4]; 
  data->AEM_S(30,0) = 0.3535533905932737*m1r[10]; 
  data->AEM_S(30,1) = 0.3162277660168379*m1r[17]+0.3535533905932737*m1r[6]; 
  data->AEM_S(30,2) = 0.3162277660168379*m1r[18]+0.3535533905932737*m1r[5]; 
  data->AEM_S(30,3) = 0.3162277660168379*m1r[19]+0.3535533905932737*m1r[4]; 
  data->AEM_S(30,4) = 0.3162277660168379*m1r[14]+0.3162277660168379*m1r[13]+0.3535533905932737*m1r[3]; 
  data->AEM_S(30,5) = 0.3162277660168379*m1r[16]+0.3162277660168379*m1r[11]+0.3535533905932737*m1r[2]; 
  data->AEM_S(30,6) = 0.3162277660168379*m1r[15]+0.3162277660168379*m1r[12]+0.3535533905932737*m1r[1]; 
  data->AEM_S(30,7) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(30,8) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(30,9) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(30,10) = 0.3162277660168379*m1r[9]+0.3162277660168379*m1r[8]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(30,11) = 0.282842712474619*m1r[18]+0.3162277660168379*m1r[5]; 
  data->AEM_S(30,12) = 0.282842712474619*m1r[17]+0.3162277660168379*m1r[6]; 
  data->AEM_S(30,13) = 0.282842712474619*m1r[19]+0.3162277660168379*m1r[4]; 
  data->AEM_S(30,14) = 0.282842712474619*m1r[19]+0.3162277660168379*m1r[4]; 
  data->AEM_S(30,15) = 0.282842712474619*m1r[17]+0.3162277660168379*m1r[6]; 
  data->AEM_S(30,16) = 0.282842712474619*m1r[18]+0.3162277660168379*m1r[5]; 
  data->AEM_S(30,17) = 0.282842712474619*m1r[15]+0.282842712474619*m1r[12]+0.3162277660168379*m1r[1]; 
  data->AEM_S(30,18) = 0.282842712474619*m1r[16]+0.282842712474619*m1r[11]+0.3162277660168379*m1r[2]; 
  data->AEM_S(30,19) = 0.282842712474619*m1r[14]+0.282842712474619*m1r[13]+0.3162277660168379*m1r[3]; 
  data->AEM_S(31,0) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(31,1) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(31,2) = 0.3535533905932737*m1r[7]; 
  data->AEM_S(31,3) = 0.3535533905932737*m1r[17]; 
  data->AEM_S(31,4) = 0.2828427124746191*m1r[12]+0.3162277660168379*m1r[1]; 
  data->AEM_S(31,5) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(31,6) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(31,7) = 0.2258769757263128*m1r[11]+0.3535533905932737*m1r[2]; 
  data->AEM_S(31,8) = 0.3162277660168379*m1r[11]; 
  data->AEM_S(31,10) = 0.282842712474619*m1r[18]+0.3162277660168379*m1r[5]; 
  data->AEM_S(31,11) = 0.3162277660168379*m1r[8]+0.2258769757263128*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(31,12) = 0.2828427124746191*m1r[4]; 
  data->AEM_S(31,13) = 0.2258769757263128*m1r[17]+0.3535533905932737*m1r[6]; 
  data->AEM_S(31,14) = 0.3162277660168379*m1r[17]; 
  data->AEM_S(31,15) = 0.3162277660168379*m1r[19]; 
  data->AEM_S(31,17) = 0.3162277660168379*m1r[14]+0.2258769757263128*m1r[13]+0.3535533905932737*m1r[3]; 
  data->AEM_S(31,18) = 0.282842712474619*m1r[10]; 
  data->AEM_S(31,19) = 0.3162277660168379*m1r[15]; 
  data->AEM_S(32,0) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(32,1) = 0.3535533905932737*m1r[8]; 
  data->AEM_S(32,2) = 0.3162277660168379*m1r[4]; 
  data->AEM_S(32,3) = 0.3535533905932737*m1r[18]; 
  data->AEM_S(32,4) = 0.2828427124746191*m1r[11]+0.3162277660168379*m1r[2]; 
  data->AEM_S(32,5) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(32,6) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(32,7) = 0.3162277660168379*m1r[12]; 
  data->AEM_S(32,8) = 0.2258769757263128*m1r[12]+0.3535533905932737*m1r[1]; 
  data->AEM_S(32,10) = 0.282842712474619*m1r[17]+0.3162277660168379*m1r[6]; 
  data->AEM_S(32,11) = 0.2828427124746191*m1r[4]; 
  data->AEM_S(32,12) = 0.2258769757263128*m1r[8]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(32,13) = 0.3162277660168379*m1r[18]; 
  data->AEM_S(32,14) = 0.2258769757263128*m1r[18]+0.3535533905932737*m1r[5]; 
  data->AEM_S(32,16) = 0.3162277660168379*m1r[19]; 
  data->AEM_S(32,17) = 0.282842712474619*m1r[10]; 
  data->AEM_S(32,18) = 0.2258769757263128*m1r[14]+0.3162277660168379*m1r[13]+0.3535533905932737*m1r[3]; 
  data->AEM_S(32,19) = 0.3162277660168379*m1r[16]; 
  data->AEM_S(33,0) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(33,1) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(33,2) = 0.3535533905932737*m1r[17]; 
  data->AEM_S(33,3) = 0.3535533905932737*m1r[7]; 
  data->AEM_S(33,4) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(33,5) = 0.2828427124746191*m1r[15]+0.3162277660168379*m1r[1]; 
  data->AEM_S(33,6) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(33,7) = 0.2258769757263128*m1r[13]+0.3535533905932737*m1r[3]; 
  data->AEM_S(33,9) = 0.3162277660168379*m1r[13]; 
  data->AEM_S(33,10) = 0.282842712474619*m1r[19]+0.3162277660168379*m1r[4]; 
  data->AEM_S(33,11) = 0.2258769757263128*m1r[17]+0.3535533905932737*m1r[6]; 
  data->AEM_S(33,12) = 0.3162277660168379*m1r[18]; 
  data->AEM_S(33,13) = 0.3162277660168379*m1r[9]+0.2258769757263128*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(33,15) = 0.2828427124746191*m1r[5]; 
  data->AEM_S(33,16) = 0.3162277660168379*m1r[17]; 
  data->AEM_S(33,17) = 0.3162277660168379*m1r[16]+0.2258769757263128*m1r[11]+0.3535533905932737*m1r[2]; 
  data->AEM_S(33,18) = 0.3162277660168379*m1r[12]; 
  data->AEM_S(33,19) = 0.282842712474619*m1r[10]; 
  data->AEM_S(34,0) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(34,1) = 0.3535533905932737*m1r[18]; 
  data->AEM_S(34,2) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(34,3) = 0.3535533905932737*m1r[8]; 
  data->AEM_S(34,4) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(34,5) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(34,6) = 0.2828427124746191*m1r[16]+0.3162277660168379*m1r[2]; 
  data->AEM_S(34,8) = 0.2258769757263128*m1r[14]+0.3535533905932737*m1r[3]; 
  data->AEM_S(34,9) = 0.3162277660168379*m1r[14]; 
  data->AEM_S(34,10) = 0.282842712474619*m1r[19]+0.3162277660168379*m1r[4]; 
  data->AEM_S(34,11) = 0.3162277660168379*m1r[17]; 
  data->AEM_S(34,12) = 0.2258769757263128*m1r[18]+0.3535533905932737*m1r[5]; 
  data->AEM_S(34,14) = 0.3162277660168379*m1r[9]+0.2258769757263128*m1r[8]+0.3535533905932737*m1r[0]; 
  data->AEM_S(34,15) = 0.3162277660168379*m1r[18]; 
  data->AEM_S(34,16) = 0.2828427124746191*m1r[6]; 
  data->AEM_S(34,17) = 0.3162277660168379*m1r[11]; 
  data->AEM_S(34,18) = 0.3162277660168379*m1r[15]+0.2258769757263128*m1r[12]+0.3535533905932737*m1r[1]; 
  data->AEM_S(34,19) = 0.282842712474619*m1r[10]; 
  data->AEM_S(35,0) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(35,1) = 0.3535533905932737*m1r[9]; 
  data->AEM_S(35,2) = 0.3535533905932737*m1r[19]; 
  data->AEM_S(35,3) = 0.3162277660168379*m1r[5]; 
  data->AEM_S(35,4) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(35,5) = 0.2828427124746191*m1r[13]+0.3162277660168379*m1r[3]; 
  data->AEM_S(35,6) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(35,7) = 0.3162277660168379*m1r[15]; 
  data->AEM_S(35,9) = 0.2258769757263128*m1r[15]+0.3535533905932737*m1r[1]; 
  data->AEM_S(35,10) = 0.282842712474619*m1r[17]+0.3162277660168379*m1r[6]; 
  data->AEM_S(35,11) = 0.3162277660168379*m1r[19]; 
  data->AEM_S(35,13) = 0.2828427124746191*m1r[5]; 
  data->AEM_S(35,14) = 0.3162277660168379*m1r[18]; 
  data->AEM_S(35,15) = 0.2258769757263128*m1r[9]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(35,16) = 0.2258769757263128*m1r[19]+0.3535533905932737*m1r[4]; 
  data->AEM_S(35,17) = 0.282842712474619*m1r[10]; 
  data->AEM_S(35,18) = 0.3162277660168379*m1r[14]; 
  data->AEM_S(35,19) = 0.2258769757263128*m1r[16]+0.3162277660168379*m1r[11]+0.3535533905932737*m1r[2]; 
  data->AEM_S(36,0) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(36,1) = 0.3535533905932737*m1r[19]; 
  data->AEM_S(36,2) = 0.3535533905932737*m1r[9]; 
  data->AEM_S(36,3) = 0.3162277660168379*m1r[6]; 
  data->AEM_S(36,4) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(36,5) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(36,6) = 0.2828427124746191*m1r[14]+0.3162277660168379*m1r[3]; 
  data->AEM_S(36,8) = 0.3162277660168379*m1r[16]; 
  data->AEM_S(36,9) = 0.2258769757263128*m1r[16]+0.3535533905932737*m1r[2]; 
  data->AEM_S(36,10) = 0.282842712474619*m1r[18]+0.3162277660168379*m1r[5]; 
  data->AEM_S(36,12) = 0.3162277660168379*m1r[19]; 
  data->AEM_S(36,13) = 0.3162277660168379*m1r[17]; 
  data->AEM_S(36,14) = 0.2828427124746191*m1r[6]; 
  data->AEM_S(36,15) = 0.2258769757263128*m1r[19]+0.3535533905932737*m1r[4]; 
  data->AEM_S(36,16) = 0.2258769757263128*m1r[9]+0.3162277660168379*m1r[8]+0.3535533905932737*m1r[0]; 
  data->AEM_S(36,17) = 0.3162277660168379*m1r[13]; 
  data->AEM_S(36,18) = 0.282842712474619*m1r[10]; 
  data->AEM_S(36,19) = 0.2258769757263128*m1r[15]+0.3162277660168379*m1r[12]+0.3535533905932737*m1r[1]; 
  data->AEM_S(37,0) = 0.3535533905932737*m1r[17]; 
  data->AEM_S(37,1) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(37,2) = 0.3535533905932737*m1r[13]; 
  data->AEM_S(37,3) = 0.3535533905932737*m1r[11]; 
  data->AEM_S(37,4) = 0.2828427124746191*m1r[18]+0.3162277660168379*m1r[5]; 
  data->AEM_S(37,5) = 0.2828427124746191*m1r[19]+0.3162277660168379*m1r[4]; 
  data->AEM_S(37,6) = 0.3535533905932737*m1r[7]; 
  data->AEM_S(37,7) = 0.2258769757263128*m1r[17]+0.3535533905932737*m1r[6]; 
  data->AEM_S(37,8) = 0.3162277660168379*m1r[17]; 
  data->AEM_S(37,9) = 0.3162277660168379*m1r[17]; 
  data->AEM_S(37,10) = 0.282842712474619*m1r[15]+0.282842712474619*m1r[12]+0.3162277660168379*m1r[1]; 
  data->AEM_S(37,11) = 0.3162277660168379*m1r[14]+0.2258769757263128*m1r[13]+0.3535533905932737*m1r[3]; 
  data->AEM_S(37,12) = 0.282842712474619*m1r[10]; 
  data->AEM_S(37,13) = 0.3162277660168379*m1r[16]+0.2258769757263128*m1r[11]+0.3535533905932737*m1r[2]; 
  data->AEM_S(37,14) = 0.3162277660168379*m1r[11]; 
  data->AEM_S(37,15) = 0.282842712474619*m1r[10]; 
  data->AEM_S(37,16) = 0.3162277660168379*m1r[13]; 
  data->AEM_S(37,17) = 0.3162277660168379*m1r[9]+0.3162277660168379*m1r[8]+0.2258769757263128*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(37,18) = 0.2529822128134704*m1r[19]+0.2828427124746191*m1r[4]; 
  data->AEM_S(37,19) = 0.2529822128134704*m1r[18]+0.2828427124746191*m1r[5]; 
  data->AEM_S(38,0) = 0.3535533905932737*m1r[18]; 
  data->AEM_S(38,1) = 0.3535533905932737*m1r[14]; 
  data->AEM_S(38,2) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(38,3) = 0.3535533905932737*m1r[12]; 
  data->AEM_S(38,4) = 0.2828427124746191*m1r[17]+0.3162277660168379*m1r[6]; 
  data->AEM_S(38,5) = 0.3535533905932737*m1r[8]; 
  data->AEM_S(38,6) = 0.2828427124746191*m1r[19]+0.3162277660168379*m1r[4]; 
  data->AEM_S(38,7) = 0.3162277660168379*m1r[18]; 
  data->AEM_S(38,8) = 0.2258769757263128*m1r[18]+0.3535533905932737*m1r[5]; 
  data->AEM_S(38,9) = 0.3162277660168379*m1r[18]; 
  data->AEM_S(38,10) = 0.282842712474619*m1r[16]+0.282842712474619*m1r[11]+0.3162277660168379*m1r[2]; 
  data->AEM_S(38,11) = 0.282842712474619*m1r[10]; 
  data->AEM_S(38,12) = 0.2258769757263128*m1r[14]+0.3162277660168379*m1r[13]+0.3535533905932737*m1r[3]; 
  data->AEM_S(38,13) = 0.3162277660168379*m1r[12]; 
  data->AEM_S(38,14) = 0.3162277660168379*m1r[15]+0.2258769757263128*m1r[12]+0.3535533905932737*m1r[1]; 
  data->AEM_S(38,15) = 0.3162277660168379*m1r[14]; 
  data->AEM_S(38,16) = 0.282842712474619*m1r[10]; 
  data->AEM_S(38,17) = 0.2529822128134704*m1r[19]+0.2828427124746191*m1r[4]; 
  data->AEM_S(38,18) = 0.3162277660168379*m1r[9]+0.2258769757263128*m1r[8]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  data->AEM_S(38,19) = 0.2529822128134704*m1r[17]+0.2828427124746191*m1r[6]; 
  data->AEM_S(39,0) = 0.3535533905932737*m1r[19]; 
  data->AEM_S(39,1) = 0.3535533905932737*m1r[16]; 
  data->AEM_S(39,2) = 0.3535533905932737*m1r[15]; 
  data->AEM_S(39,3) = 0.3162277660168379*m1r[10]; 
  data->AEM_S(39,4) = 0.3535533905932737*m1r[9]; 
  data->AEM_S(39,5) = 0.2828427124746191*m1r[17]+0.3162277660168379*m1r[6]; 
  data->AEM_S(39,6) = 0.2828427124746191*m1r[18]+0.3162277660168379*m1r[5]; 
  data->AEM_S(39,7) = 0.3162277660168379*m1r[19]; 
  data->AEM_S(39,8) = 0.3162277660168379*m1r[19]; 
  data->AEM_S(39,9) = 0.2258769757263128*m1r[19]+0.3535533905932737*m1r[4]; 
  data->AEM_S(39,10) = 0.282842712474619*m1r[14]+0.282842712474619*m1r[13]+0.3162277660168379*m1r[3]; 
  data->AEM_S(39,11) = 0.3162277660168379*m1r[15]; 
  data->AEM_S(39,12) = 0.3162277660168379*m1r[16]; 
  data->AEM_S(39,13) = 0.282842712474619*m1r[10]; 
  data->AEM_S(39,14) = 0.282842712474619*m1r[10]; 
  data->AEM_S(39,15) = 0.2258769757263128*m1r[16]+0.3162277660168379*m1r[11]+0.3535533905932737*m1r[2]; 
  data->AEM_S(39,16) = 0.2258769757263128*m1r[15]+0.3162277660168379*m1r[12]+0.3535533905932737*m1r[1]; 
  data->AEM_S(39,17) = 0.2529822128134704*m1r[18]+0.2828427124746191*m1r[5]; 
  data->AEM_S(39,18) = 0.2529822128134704*m1r[17]+0.2828427124746191*m1r[6]; 
  data->AEM_S(39,19) = 0.2258769757263128*m1r[9]+0.3162277660168379*m1r[8]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(20,20) = 1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(20,21) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(20,22) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(20,23) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(20,24) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(20,25) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(20,26) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(20,27) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(20,28) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  data->AEM_S(20,29) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  data->AEM_S(20,30) = 1.060660171779821*m0r[10]-0.3535533905932737*cE[10]; 
  data->AEM_S(20,31) = 1.060660171779821*m0r[11]-0.3535533905932737*cE[11]; 
  data->AEM_S(20,32) = 1.060660171779821*m0r[12]-0.3535533905932737*cE[12]; 
  data->AEM_S(20,33) = 1.060660171779821*m0r[13]-0.3535533905932737*cE[13]; 
  data->AEM_S(20,34) = 1.060660171779821*m0r[14]-0.3535533905932737*cE[14]; 
  data->AEM_S(20,35) = 1.060660171779821*m0r[15]-0.3535533905932737*cE[15]; 
  data->AEM_S(20,36) = 1.060660171779821*m0r[16]-0.3535533905932737*cE[16]; 
  data->AEM_S(20,37) = 1.060660171779821*m0r[17]-0.3535533905932737*cE[17]; 
  data->AEM_S(20,38) = 1.060660171779821*m0r[18]-0.3535533905932737*cE[18]; 
  data->AEM_S(20,39) = 1.060660171779821*m0r[19]-0.3535533905932737*cE[19]; 
  data->AEM_S(21,20) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(21,21) = 0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(21,22) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(21,23) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(21,24) = 0.9486832980505138*m0r[11]-0.3162277660168379*cE[11]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(21,25) = 0.9486832980505138*m0r[13]-0.3162277660168379*cE[13]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(21,26) = 1.060660171779821*m0r[10]-0.3535533905932737*cE[10]; 
  data->AEM_S(21,27) = 0.9486832980505137*m0r[1]-0.3162277660168379*cE[1]; 
  data->AEM_S(21,28) = 1.060660171779821*m0r[12]-0.3535533905932737*cE[12]; 
  data->AEM_S(21,29) = 1.060660171779821*m0r[15]-0.3535533905932737*cE[15]; 
  data->AEM_S(21,30) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]+1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(21,31) = 0.9486832980505138*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(21,32) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  data->AEM_S(21,33) = 0.9486832980505138*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(21,34) = 1.060660171779821*m0r[18]-0.3535533905932737*cE[18]; 
  data->AEM_S(21,35) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  data->AEM_S(21,36) = 1.060660171779821*m0r[19]-0.3535533905932737*cE[19]; 
  data->AEM_S(21,37) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(21,38) = 1.060660171779821*m0r[14]-0.3535533905932737*cE[14]; 
  data->AEM_S(21,39) = 1.060660171779821*m0r[16]-0.3535533905932737*cE[16]; 
  data->AEM_S(22,20) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(22,21) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(22,22) = 0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(22,23) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(22,24) = 0.9486832980505138*m0r[12]-0.3162277660168379*cE[12]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(22,25) = 1.060660171779821*m0r[10]-0.3535533905932737*cE[10]; 
  data->AEM_S(22,26) = 0.9486832980505138*m0r[14]-0.3162277660168379*cE[14]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(22,27) = 1.060660171779821*m0r[11]-0.3535533905932737*cE[11]; 
  data->AEM_S(22,28) = 0.9486832980505137*m0r[2]-0.3162277660168379*cE[2]; 
  data->AEM_S(22,29) = 1.060660171779821*m0r[16]-0.3535533905932737*cE[16]; 
  data->AEM_S(22,30) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]+1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(22,31) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(22,32) = 0.9486832980505138*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(22,33) = 1.060660171779821*m0r[17]-0.3535533905932737*cE[17]; 
  data->AEM_S(22,34) = 0.9486832980505138*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(22,35) = 1.060660171779821*m0r[19]-0.3535533905932737*cE[19]; 
  data->AEM_S(22,36) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  data->AEM_S(22,37) = 1.060660171779821*m0r[13]-0.3535533905932737*cE[13]; 
  data->AEM_S(22,38) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(22,39) = 1.060660171779821*m0r[15]-0.3535533905932737*cE[15]; 
  data->AEM_S(23,20) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(23,21) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(23,22) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(23,23) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(23,24) = 1.060660171779821*m0r[10]-0.3535533905932737*cE[10]; 
  data->AEM_S(23,25) = 0.9486832980505138*m0r[15]-0.3162277660168379*cE[15]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(23,26) = 0.9486832980505138*m0r[16]-0.3162277660168379*cE[16]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(23,27) = 1.060660171779821*m0r[13]-0.3535533905932737*cE[13]; 
  data->AEM_S(23,28) = 1.060660171779821*m0r[14]-0.3535533905932737*cE[14]; 
  data->AEM_S(23,29) = 0.9486832980505137*m0r[3]-0.3162277660168379*cE[3]; 
  data->AEM_S(23,30) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]+1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(23,31) = 1.060660171779821*m0r[17]-0.3535533905932737*cE[17]; 
  data->AEM_S(23,32) = 1.060660171779821*m0r[18]-0.3535533905932737*cE[18]; 
  data->AEM_S(23,33) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(23,34) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  data->AEM_S(23,35) = 0.9486832980505138*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(23,36) = 0.9486832980505138*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(23,37) = 1.060660171779821*m0r[11]-0.3535533905932737*cE[11]; 
  data->AEM_S(23,38) = 1.060660171779821*m0r[12]-0.3535533905932737*cE[12]; 
  data->AEM_S(23,39) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(24,20) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(24,21) = 0.9486832980505138*m0r[11]-0.3162277660168379*cE[11]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(24,22) = 0.9486832980505138*m0r[12]-0.3162277660168379*cE[12]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(24,23) = 1.060660171779821*m0r[10]-0.3535533905932737*cE[10]; 
  data->AEM_S(24,24) = 0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(24,25) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]+1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(24,26) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]+1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(24,27) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(24,28) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(24,29) = 1.060660171779821*m0r[19]-0.3535533905932737*cE[19]; 
  data->AEM_S(24,30) = 0.9486832980505138*m0r[14]-0.3162277660168379*cE[14]+0.9486832980505138*m0r[13]-0.3162277660168379*cE[13]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(24,31) = 0.848528137423857*m0r[12]-0.2828427124746191*cE[12]+0.9486832980505138*m0r[1]-0.3162277660168379*cE[1]; 
  data->AEM_S(24,32) = 0.848528137423857*m0r[11]-0.2828427124746191*cE[11]+0.9486832980505138*m0r[2]-0.3162277660168379*cE[2]; 
  data->AEM_S(24,33) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(24,34) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(24,35) = 1.060660171779821*m0r[16]-0.3535533905932737*cE[16]; 
  data->AEM_S(24,36) = 1.060660171779821*m0r[15]-0.3535533905932737*cE[15]; 
  data->AEM_S(24,37) = 0.848528137423857*m0r[18]-0.2828427124746191*cE[18]+0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(24,38) = 0.848528137423857*m0r[17]-0.2828427124746191*cE[17]+0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(24,39) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  data->AEM_S(25,20) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(25,21) = 0.9486832980505138*m0r[13]-0.3162277660168379*cE[13]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(25,22) = 1.060660171779821*m0r[10]-0.3535533905932737*cE[10]; 
  data->AEM_S(25,23) = 0.9486832980505138*m0r[15]-0.3162277660168379*cE[15]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(25,24) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]+1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(25,25) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(25,26) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]+1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(25,27) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(25,28) = 1.060660171779821*m0r[18]-0.3535533905932737*cE[18]; 
  data->AEM_S(25,29) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(25,30) = 0.9486832980505138*m0r[16]-0.3162277660168379*cE[16]+0.9486832980505138*m0r[11]-0.3162277660168379*cE[11]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(25,31) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(25,32) = 1.060660171779821*m0r[14]-0.3535533905932737*cE[14]; 
  data->AEM_S(25,33) = 0.848528137423857*m0r[15]-0.2828427124746191*cE[15]+0.9486832980505138*m0r[1]-0.3162277660168379*cE[1]; 
  data->AEM_S(25,34) = 1.060660171779821*m0r[12]-0.3535533905932737*cE[12]; 
  data->AEM_S(25,35) = 0.848528137423857*m0r[13]-0.2828427124746191*cE[13]+0.9486832980505138*m0r[3]-0.3162277660168379*cE[3]; 
  data->AEM_S(25,36) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(25,37) = 0.848528137423857*m0r[19]-0.2828427124746191*cE[19]+0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(25,38) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  data->AEM_S(25,39) = 0.848528137423857*m0r[17]-0.2828427124746191*cE[17]+0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(26,20) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(26,21) = 1.060660171779821*m0r[10]-0.3535533905932737*cE[10]; 
  data->AEM_S(26,22) = 0.9486832980505138*m0r[14]-0.3162277660168379*cE[14]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(26,23) = 0.9486832980505138*m0r[16]-0.3162277660168379*cE[16]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(26,24) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]+1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(26,25) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]+1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(26,26) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(26,27) = 1.060660171779821*m0r[17]-0.3535533905932737*cE[17]; 
  data->AEM_S(26,28) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(26,29) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(26,30) = 0.9486832980505138*m0r[15]-0.3162277660168379*cE[15]+0.9486832980505138*m0r[12]-0.3162277660168379*cE[12]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(26,31) = 1.060660171779821*m0r[13]-0.3535533905932737*cE[13]; 
  data->AEM_S(26,32) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(26,33) = 1.060660171779821*m0r[11]-0.3535533905932737*cE[11]; 
  data->AEM_S(26,34) = 0.848528137423857*m0r[16]-0.2828427124746191*cE[16]+0.9486832980505138*m0r[2]-0.3162277660168379*cE[2]; 
  data->AEM_S(26,35) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(26,36) = 0.848528137423857*m0r[14]-0.2828427124746191*cE[14]+0.9486832980505138*m0r[3]-0.3162277660168379*cE[3]; 
  data->AEM_S(26,37) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(26,38) = 0.848528137423857*m0r[19]-0.2828427124746191*cE[19]+0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(26,39) = 0.848528137423857*m0r[18]-0.2828427124746191*cE[18]+0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(27,20) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(27,21) = 0.9486832980505137*m0r[1]-0.3162277660168379*cE[1]; 
  data->AEM_S(27,22) = 1.060660171779821*m0r[11]-0.3535533905932737*cE[11]; 
  data->AEM_S(27,23) = 1.060660171779821*m0r[13]-0.3535533905932737*cE[13]; 
  data->AEM_S(27,24) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(27,25) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(27,26) = 1.060660171779821*m0r[17]-0.3535533905932737*cE[17]; 
  data->AEM_S(27,27) = 0.6776309271789384*m0r[7]-0.2258769757263128*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(27,30) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(27,31) = 0.6776309271789384*m0r[11]-0.2258769757263128*cE[11]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(27,32) = 0.9486832980505137*m0r[12]-0.3162277660168379*cE[12]; 
  data->AEM_S(27,33) = 0.6776309271789384*m0r[13]-0.2258769757263128*cE[13]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(27,35) = 0.9486832980505137*m0r[15]-0.3162277660168379*cE[15]; 
  data->AEM_S(27,37) = 0.6776309271789384*m0r[17]-0.2258769757263128*cE[17]+1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(27,38) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]; 
  data->AEM_S(27,39) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]; 
  data->AEM_S(28,20) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  data->AEM_S(28,21) = 1.060660171779821*m0r[12]-0.3535533905932737*cE[12]; 
  data->AEM_S(28,22) = 0.9486832980505137*m0r[2]-0.3162277660168379*cE[2]; 
  data->AEM_S(28,23) = 1.060660171779821*m0r[14]-0.3535533905932737*cE[14]; 
  data->AEM_S(28,24) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(28,25) = 1.060660171779821*m0r[18]-0.3535533905932737*cE[18]; 
  data->AEM_S(28,26) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(28,28) = 0.6776309271789384*m0r[8]-0.2258769757263128*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(28,30) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(28,31) = 0.9486832980505137*m0r[11]-0.3162277660168379*cE[11]; 
  data->AEM_S(28,32) = 0.6776309271789384*m0r[12]-0.2258769757263128*cE[12]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(28,34) = 0.6776309271789384*m0r[14]-0.2258769757263128*cE[14]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(28,36) = 0.9486832980505137*m0r[16]-0.3162277660168379*cE[16]; 
  data->AEM_S(28,37) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]; 
  data->AEM_S(28,38) = 0.6776309271789384*m0r[18]-0.2258769757263128*cE[18]+1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(28,39) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]; 
  data->AEM_S(29,20) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  data->AEM_S(29,21) = 1.060660171779821*m0r[15]-0.3535533905932737*cE[15]; 
  data->AEM_S(29,22) = 1.060660171779821*m0r[16]-0.3535533905932737*cE[16]; 
  data->AEM_S(29,23) = 0.9486832980505137*m0r[3]-0.3162277660168379*cE[3]; 
  data->AEM_S(29,24) = 1.060660171779821*m0r[19]-0.3535533905932737*cE[19]; 
  data->AEM_S(29,25) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(29,26) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(29,29) = 0.6776309271789384*m0r[9]-0.2258769757263128*cE[9]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(29,30) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(29,33) = 0.9486832980505137*m0r[13]-0.3162277660168379*cE[13]; 
  data->AEM_S(29,34) = 0.9486832980505137*m0r[14]-0.3162277660168379*cE[14]; 
  data->AEM_S(29,35) = 0.6776309271789384*m0r[15]-0.2258769757263128*cE[15]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(29,36) = 0.6776309271789384*m0r[16]-0.2258769757263128*cE[16]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(29,37) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]; 
  data->AEM_S(29,38) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]; 
  data->AEM_S(29,39) = 0.6776309271789384*m0r[19]-0.2258769757263128*cE[19]+1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(30,20) = 1.060660171779821*m0r[10]-0.3535533905932737*cE[10]; 
  data->AEM_S(30,21) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]+1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(30,22) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]+1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(30,23) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]+1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(30,24) = 0.9486832980505138*m0r[14]-0.3162277660168379*cE[14]+0.9486832980505138*m0r[13]-0.3162277660168379*cE[13]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(30,25) = 0.9486832980505138*m0r[16]-0.3162277660168379*cE[16]+0.9486832980505138*m0r[11]-0.3162277660168379*cE[11]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(30,26) = 0.9486832980505138*m0r[15]-0.3162277660168379*cE[15]+0.9486832980505138*m0r[12]-0.3162277660168379*cE[12]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(30,27) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(30,28) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(30,29) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(30,30) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(30,31) = 0.848528137423857*m0r[18]-0.282842712474619*cE[18]+0.9486832980505138*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(30,32) = 0.848528137423857*m0r[17]-0.282842712474619*cE[17]+0.9486832980505138*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(30,33) = 0.848528137423857*m0r[19]-0.282842712474619*cE[19]+0.9486832980505138*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(30,34) = 0.848528137423857*m0r[19]-0.282842712474619*cE[19]+0.9486832980505138*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(30,35) = 0.848528137423857*m0r[17]-0.282842712474619*cE[17]+0.9486832980505138*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(30,36) = 0.848528137423857*m0r[18]-0.282842712474619*cE[18]+0.9486832980505138*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(30,37) = 0.848528137423857*m0r[15]-0.282842712474619*cE[15]+0.848528137423857*m0r[12]-0.282842712474619*cE[12]+0.9486832980505137*m0r[1]-0.3162277660168379*cE[1]; 
  data->AEM_S(30,38) = 0.848528137423857*m0r[16]-0.282842712474619*cE[16]+0.848528137423857*m0r[11]-0.282842712474619*cE[11]+0.9486832980505137*m0r[2]-0.3162277660168379*cE[2]; 
  data->AEM_S(30,39) = 0.848528137423857*m0r[14]-0.282842712474619*cE[14]+0.848528137423857*m0r[13]-0.282842712474619*cE[13]+0.9486832980505137*m0r[3]-0.3162277660168379*cE[3]; 
  data->AEM_S(31,20) = 1.060660171779821*m0r[11]-0.3535533905932737*cE[11]; 
  data->AEM_S(31,21) = 0.9486832980505138*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(31,22) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(31,23) = 1.060660171779821*m0r[17]-0.3535533905932737*cE[17]; 
  data->AEM_S(31,24) = 0.848528137423857*m0r[12]-0.2828427124746191*cE[12]+0.9486832980505138*m0r[1]-0.3162277660168379*cE[1]; 
  data->AEM_S(31,25) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(31,26) = 1.060660171779821*m0r[13]-0.3535533905932737*cE[13]; 
  data->AEM_S(31,27) = 0.6776309271789384*m0r[11]-0.2258769757263128*cE[11]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(31,28) = 0.9486832980505137*m0r[11]-0.3162277660168379*cE[11]; 
  data->AEM_S(31,30) = 0.848528137423857*m0r[18]-0.282842712474619*cE[18]+0.9486832980505138*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(31,31) = 0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+0.6776309271789384*m0r[7]-0.2258769757263128*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(31,32) = 0.848528137423857*m0r[4]-0.2828427124746191*cE[4]; 
  data->AEM_S(31,33) = 0.6776309271789384*m0r[17]-0.2258769757263128*cE[17]+1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(31,34) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]; 
  data->AEM_S(31,35) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]; 
  data->AEM_S(31,37) = 0.9486832980505137*m0r[14]-0.3162277660168379*cE[14]+0.6776309271789384*m0r[13]-0.2258769757263128*cE[13]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(31,38) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(31,39) = 0.9486832980505137*m0r[15]-0.3162277660168379*cE[15]; 
  data->AEM_S(32,20) = 1.060660171779821*m0r[12]-0.3535533905932737*cE[12]; 
  data->AEM_S(32,21) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  data->AEM_S(32,22) = 0.9486832980505138*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(32,23) = 1.060660171779821*m0r[18]-0.3535533905932737*cE[18]; 
  data->AEM_S(32,24) = 0.848528137423857*m0r[11]-0.2828427124746191*cE[11]+0.9486832980505138*m0r[2]-0.3162277660168379*cE[2]; 
  data->AEM_S(32,25) = 1.060660171779821*m0r[14]-0.3535533905932737*cE[14]; 
  data->AEM_S(32,26) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(32,27) = 0.9486832980505137*m0r[12]-0.3162277660168379*cE[12]; 
  data->AEM_S(32,28) = 0.6776309271789384*m0r[12]-0.2258769757263128*cE[12]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(32,30) = 0.848528137423857*m0r[17]-0.282842712474619*cE[17]+0.9486832980505138*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(32,31) = 0.848528137423857*m0r[4]-0.2828427124746191*cE[4]; 
  data->AEM_S(32,32) = 0.6776309271789384*m0r[8]-0.2258769757263128*cE[8]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(32,33) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]; 
  data->AEM_S(32,34) = 0.6776309271789384*m0r[18]-0.2258769757263128*cE[18]+1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(32,36) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]; 
  data->AEM_S(32,37) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(32,38) = 0.6776309271789384*m0r[14]-0.2258769757263128*cE[14]+0.9486832980505137*m0r[13]-0.3162277660168379*cE[13]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(32,39) = 0.9486832980505137*m0r[16]-0.3162277660168379*cE[16]; 
  data->AEM_S(33,20) = 1.060660171779821*m0r[13]-0.3535533905932737*cE[13]; 
  data->AEM_S(33,21) = 0.9486832980505138*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(33,22) = 1.060660171779821*m0r[17]-0.3535533905932737*cE[17]; 
  data->AEM_S(33,23) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(33,24) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(33,25) = 0.848528137423857*m0r[15]-0.2828427124746191*cE[15]+0.9486832980505138*m0r[1]-0.3162277660168379*cE[1]; 
  data->AEM_S(33,26) = 1.060660171779821*m0r[11]-0.3535533905932737*cE[11]; 
  data->AEM_S(33,27) = 0.6776309271789384*m0r[13]-0.2258769757263128*cE[13]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(33,29) = 0.9486832980505137*m0r[13]-0.3162277660168379*cE[13]; 
  data->AEM_S(33,30) = 0.848528137423857*m0r[19]-0.282842712474619*cE[19]+0.9486832980505138*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(33,31) = 0.6776309271789384*m0r[17]-0.2258769757263128*cE[17]+1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(33,32) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]; 
  data->AEM_S(33,33) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.6776309271789384*m0r[7]-0.2258769757263128*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(33,35) = 0.848528137423857*m0r[5]-0.2828427124746191*cE[5]; 
  data->AEM_S(33,36) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]; 
  data->AEM_S(33,37) = 0.9486832980505137*m0r[16]-0.3162277660168379*cE[16]+0.6776309271789384*m0r[11]-0.2258769757263128*cE[11]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(33,38) = 0.9486832980505137*m0r[12]-0.3162277660168379*cE[12]; 
  data->AEM_S(33,39) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(34,20) = 1.060660171779821*m0r[14]-0.3535533905932737*cE[14]; 
  data->AEM_S(34,21) = 1.060660171779821*m0r[18]-0.3535533905932737*cE[18]; 
  data->AEM_S(34,22) = 0.9486832980505138*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(34,23) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  data->AEM_S(34,24) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(34,25) = 1.060660171779821*m0r[12]-0.3535533905932737*cE[12]; 
  data->AEM_S(34,26) = 0.848528137423857*m0r[16]-0.2828427124746191*cE[16]+0.9486832980505138*m0r[2]-0.3162277660168379*cE[2]; 
  data->AEM_S(34,28) = 0.6776309271789384*m0r[14]-0.2258769757263128*cE[14]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(34,29) = 0.9486832980505137*m0r[14]-0.3162277660168379*cE[14]; 
  data->AEM_S(34,30) = 0.848528137423857*m0r[19]-0.282842712474619*cE[19]+0.9486832980505138*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(34,31) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]; 
  data->AEM_S(34,32) = 0.6776309271789384*m0r[18]-0.2258769757263128*cE[18]+1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(34,34) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.6776309271789384*m0r[8]-0.2258769757263128*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(34,35) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]; 
  data->AEM_S(34,36) = 0.848528137423857*m0r[6]-0.2828427124746191*cE[6]; 
  data->AEM_S(34,37) = 0.9486832980505137*m0r[11]-0.3162277660168379*cE[11]; 
  data->AEM_S(34,38) = 0.9486832980505137*m0r[15]-0.3162277660168379*cE[15]+0.6776309271789384*m0r[12]-0.2258769757263128*cE[12]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(34,39) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(35,20) = 1.060660171779821*m0r[15]-0.3535533905932737*cE[15]; 
  data->AEM_S(35,21) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  data->AEM_S(35,22) = 1.060660171779821*m0r[19]-0.3535533905932737*cE[19]; 
  data->AEM_S(35,23) = 0.9486832980505138*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(35,24) = 1.060660171779821*m0r[16]-0.3535533905932737*cE[16]; 
  data->AEM_S(35,25) = 0.848528137423857*m0r[13]-0.2828427124746191*cE[13]+0.9486832980505138*m0r[3]-0.3162277660168379*cE[3]; 
  data->AEM_S(35,26) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(35,27) = 0.9486832980505137*m0r[15]-0.3162277660168379*cE[15]; 
  data->AEM_S(35,29) = 0.6776309271789384*m0r[15]-0.2258769757263128*cE[15]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(35,30) = 0.848528137423857*m0r[17]-0.282842712474619*cE[17]+0.9486832980505138*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(35,31) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]; 
  data->AEM_S(35,33) = 0.848528137423857*m0r[5]-0.2828427124746191*cE[5]; 
  data->AEM_S(35,34) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]; 
  data->AEM_S(35,35) = 0.6776309271789384*m0r[9]-0.2258769757263128*cE[9]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(35,36) = 0.6776309271789384*m0r[19]-0.2258769757263128*cE[19]+1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(35,37) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(35,38) = 0.9486832980505137*m0r[14]-0.3162277660168379*cE[14]; 
  data->AEM_S(35,39) = 0.6776309271789384*m0r[16]-0.2258769757263128*cE[16]+0.9486832980505137*m0r[11]-0.3162277660168379*cE[11]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(36,20) = 1.060660171779821*m0r[16]-0.3535533905932737*cE[16]; 
  data->AEM_S(36,21) = 1.060660171779821*m0r[19]-0.3535533905932737*cE[19]; 
  data->AEM_S(36,22) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  data->AEM_S(36,23) = 0.9486832980505138*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(36,24) = 1.060660171779821*m0r[15]-0.3535533905932737*cE[15]; 
  data->AEM_S(36,25) = 0.9486832980505138*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(36,26) = 0.848528137423857*m0r[14]-0.2828427124746191*cE[14]+0.9486832980505138*m0r[3]-0.3162277660168379*cE[3]; 
  data->AEM_S(36,28) = 0.9486832980505137*m0r[16]-0.3162277660168379*cE[16]; 
  data->AEM_S(36,29) = 0.6776309271789384*m0r[16]-0.2258769757263128*cE[16]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(36,30) = 0.848528137423857*m0r[18]-0.282842712474619*cE[18]+0.9486832980505138*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(36,32) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]; 
  data->AEM_S(36,33) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]; 
  data->AEM_S(36,34) = 0.848528137423857*m0r[6]-0.2828427124746191*cE[6]; 
  data->AEM_S(36,35) = 0.6776309271789384*m0r[19]-0.2258769757263128*cE[19]+1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(36,36) = 0.6776309271789384*m0r[9]-0.2258769757263128*cE[9]+0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(36,37) = 0.9486832980505137*m0r[13]-0.3162277660168379*cE[13]; 
  data->AEM_S(36,38) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(36,39) = 0.6776309271789384*m0r[15]-0.2258769757263128*cE[15]+0.9486832980505137*m0r[12]-0.3162277660168379*cE[12]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(37,20) = 1.060660171779821*m0r[17]-0.3535533905932737*cE[17]; 
  data->AEM_S(37,21) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(37,22) = 1.060660171779821*m0r[13]-0.3535533905932737*cE[13]; 
  data->AEM_S(37,23) = 1.060660171779821*m0r[11]-0.3535533905932737*cE[11]; 
  data->AEM_S(37,24) = 0.848528137423857*m0r[18]-0.2828427124746191*cE[18]+0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(37,25) = 0.848528137423857*m0r[19]-0.2828427124746191*cE[19]+0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(37,26) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  data->AEM_S(37,27) = 0.6776309271789384*m0r[17]-0.2258769757263128*cE[17]+1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  data->AEM_S(37,28) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]; 
  data->AEM_S(37,29) = 0.9486832980505137*m0r[17]-0.3162277660168379*cE[17]; 
  data->AEM_S(37,30) = 0.848528137423857*m0r[15]-0.282842712474619*cE[15]+0.848528137423857*m0r[12]-0.282842712474619*cE[12]+0.9486832980505137*m0r[1]-0.3162277660168379*cE[1]; 
  data->AEM_S(37,31) = 0.9486832980505137*m0r[14]-0.3162277660168379*cE[14]+0.6776309271789384*m0r[13]-0.2258769757263128*cE[13]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(37,32) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(37,33) = 0.9486832980505137*m0r[16]-0.3162277660168379*cE[16]+0.6776309271789384*m0r[11]-0.2258769757263128*cE[11]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(37,34) = 0.9486832980505137*m0r[11]-0.3162277660168379*cE[11]; 
  data->AEM_S(37,35) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(37,36) = 0.9486832980505137*m0r[13]-0.3162277660168379*cE[13]; 
  data->AEM_S(37,37) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+0.6776309271789384*m0r[7]-0.2258769757263128*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(37,38) = 0.758946638440411*m0r[19]-0.2529822128134704*cE[19]+0.848528137423857*m0r[4]-0.2828427124746191*cE[4]; 
  data->AEM_S(37,39) = 0.758946638440411*m0r[18]-0.2529822128134704*cE[18]+0.848528137423857*m0r[5]-0.2828427124746191*cE[5]; 
  data->AEM_S(38,20) = 1.060660171779821*m0r[18]-0.3535533905932737*cE[18]; 
  data->AEM_S(38,21) = 1.060660171779821*m0r[14]-0.3535533905932737*cE[14]; 
  data->AEM_S(38,22) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(38,23) = 1.060660171779821*m0r[12]-0.3535533905932737*cE[12]; 
  data->AEM_S(38,24) = 0.848528137423857*m0r[17]-0.2828427124746191*cE[17]+0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(38,25) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  data->AEM_S(38,26) = 0.848528137423857*m0r[19]-0.2828427124746191*cE[19]+0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  data->AEM_S(38,27) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]; 
  data->AEM_S(38,28) = 0.6776309271789384*m0r[18]-0.2258769757263128*cE[18]+1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  data->AEM_S(38,29) = 0.9486832980505137*m0r[18]-0.3162277660168379*cE[18]; 
  data->AEM_S(38,30) = 0.848528137423857*m0r[16]-0.282842712474619*cE[16]+0.848528137423857*m0r[11]-0.282842712474619*cE[11]+0.9486832980505137*m0r[2]-0.3162277660168379*cE[2]; 
  data->AEM_S(38,31) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(38,32) = 0.6776309271789384*m0r[14]-0.2258769757263128*cE[14]+0.9486832980505137*m0r[13]-0.3162277660168379*cE[13]+1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  data->AEM_S(38,33) = 0.9486832980505137*m0r[12]-0.3162277660168379*cE[12]; 
  data->AEM_S(38,34) = 0.9486832980505137*m0r[15]-0.3162277660168379*cE[15]+0.6776309271789384*m0r[12]-0.2258769757263128*cE[12]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(38,35) = 0.9486832980505137*m0r[14]-0.3162277660168379*cE[14]; 
  data->AEM_S(38,36) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(38,37) = 0.758946638440411*m0r[19]-0.2529822128134704*cE[19]+0.848528137423857*m0r[4]-0.2828427124746191*cE[4]; 
  data->AEM_S(38,38) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.6776309271789384*m0r[8]-0.2258769757263128*cE[8]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  data->AEM_S(38,39) = 0.758946638440411*m0r[17]-0.2529822128134704*cE[17]+0.848528137423857*m0r[6]-0.2828427124746191*cE[6]; 
  data->AEM_S(39,20) = 1.060660171779821*m0r[19]-0.3535533905932737*cE[19]; 
  data->AEM_S(39,21) = 1.060660171779821*m0r[16]-0.3535533905932737*cE[16]; 
  data->AEM_S(39,22) = 1.060660171779821*m0r[15]-0.3535533905932737*cE[15]; 
  data->AEM_S(39,23) = 0.9486832980505137*m0r[10]-0.3162277660168379*cE[10]; 
  data->AEM_S(39,24) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  data->AEM_S(39,25) = 0.848528137423857*m0r[17]-0.2828427124746191*cE[17]+0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  data->AEM_S(39,26) = 0.848528137423857*m0r[18]-0.2828427124746191*cE[18]+0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  data->AEM_S(39,27) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]; 
  data->AEM_S(39,28) = 0.9486832980505137*m0r[19]-0.3162277660168379*cE[19]; 
  data->AEM_S(39,29) = 0.6776309271789384*m0r[19]-0.2258769757263128*cE[19]+1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  data->AEM_S(39,30) = 0.848528137423857*m0r[14]-0.282842712474619*cE[14]+0.848528137423857*m0r[13]-0.282842712474619*cE[13]+0.9486832980505137*m0r[3]-0.3162277660168379*cE[3]; 
  data->AEM_S(39,31) = 0.9486832980505137*m0r[15]-0.3162277660168379*cE[15]; 
  data->AEM_S(39,32) = 0.9486832980505137*m0r[16]-0.3162277660168379*cE[16]; 
  data->AEM_S(39,33) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(39,34) = 0.848528137423857*m0r[10]-0.282842712474619*cE[10]; 
  data->AEM_S(39,35) = 0.6776309271789384*m0r[16]-0.2258769757263128*cE[16]+0.9486832980505137*m0r[11]-0.3162277660168379*cE[11]+1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  data->AEM_S(39,36) = 0.6776309271789384*m0r[15]-0.2258769757263128*cE[15]+0.9486832980505137*m0r[12]-0.3162277660168379*cE[12]+1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  data->AEM_S(39,37) = 0.758946638440411*m0r[18]-0.2529822128134704*cE[18]+0.848528137423857*m0r[5]-0.2828427124746191*cE[5]; 
  data->AEM_S(39,38) = 0.758946638440411*m0r[17]-0.2529822128134704*cE[17]+0.848528137423857*m0r[6]-0.2828427124746191*cE[6]; 
  data->AEM_S(39,39) = 0.6776309271789384*m0r[9]-0.2258769757263128*cE[9]+0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m1r[12],m1r[13],m1r[14],m1r[15],m1r[16],m1r[17],m1r[18],m1r[19],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5],m2r[6],m2r[7],m2r[8],m2r[9],m2r[10],m2r[11],m2r[12],m2r[13],m2r[14],m2r[15],m2r[16],m2r[17],m2r[18],m2r[19]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,20,1) = data->u_S.segment<20>(0); 
 
  Eigen::Map<VectorXd>(vtSq,20,1) = data->u_S.segment<20>(20); 
 
} 
 
