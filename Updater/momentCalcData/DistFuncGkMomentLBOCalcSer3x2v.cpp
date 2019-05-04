#include <DistFuncMomentCalcModDecl.h> 
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
 
