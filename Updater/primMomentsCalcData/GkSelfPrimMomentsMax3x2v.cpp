#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
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
 
void GkSelfPrimMoments3x2vMax_P1(const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1:       moments of the distribution function. 
  // m0S,m1S,m1S: star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[4]; 
  double m1r[4]; 
  double m0Sr[4]; 
  double m1Sr[4]; 
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
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m0Sr[2] = 0.0; 
    m0Sr[3] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = 0.0; 
    m1Sr[3] = 0.0; 
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
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m0Sr[2] = m0S[2]; 
    m0Sr[3] = m0S[3]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
    m2Sr[2] = m2S[2]; 
    m2Sr[3] = m2S[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(8,8); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(8);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(8);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.3535533905932737*m0r[0]; 
  BigAEM(0,1) = 0.3535533905932737*m0r[1]; 
  BigAEM(0,2) = 0.3535533905932737*m0r[2]; 
  BigAEM(0,3) = 0.3535533905932737*m0r[3]; 
  BigAEM(1,0) = 0.3535533905932737*m0r[1]; 
  BigAEM(1,1) = 0.3535533905932737*m0r[0]; 
  BigAEM(2,0) = 0.3535533905932737*m0r[2]; 
  BigAEM(2,2) = 0.3535533905932737*m0r[0]; 
  BigAEM(3,0) = 0.3535533905932737*m0r[3]; 
  BigAEM(3,3) = 0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,4) = -0.3535533905932737*cM[0]; 
  BigAEM(0,5) = -0.3535533905932737*cM[1]; 
  BigAEM(0,6) = -0.3535533905932737*cM[2]; 
  BigAEM(0,7) = -0.3535533905932737*cM[3]; 
  BigAEM(1,4) = -0.3535533905932737*cM[1]; 
  BigAEM(1,5) = -0.3535533905932737*cM[0]; 
  BigAEM(2,4) = -0.3535533905932737*cM[2]; 
  BigAEM(2,6) = -0.3535533905932737*cM[0]; 
  BigAEM(3,4) = -0.3535533905932737*cM[3]; 
  BigAEM(3,7) = -0.3535533905932737*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(4,0) = 0.3535533905932737*m1Sr[0]; 
  BigAEM(4,1) = 0.3535533905932737*m1Sr[1]; 
  BigAEM(4,2) = 0.3535533905932737*m1Sr[2]; 
  BigAEM(4,3) = 0.3535533905932737*m1Sr[3]; 
  BigAEM(5,0) = 0.3535533905932737*m1Sr[1]; 
  BigAEM(5,1) = 0.3535533905932737*m1Sr[0]; 
  BigAEM(6,0) = 0.3535533905932737*m1Sr[2]; 
  BigAEM(6,2) = 0.3535533905932737*m1Sr[0]; 
  BigAEM(7,0) = 0.3535533905932737*m1Sr[3]; 
  BigAEM(7,3) = 0.3535533905932737*m1Sr[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(4,4) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  BigAEM(4,5) = 0.7071067811865475*m0r[1]+0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  BigAEM(4,6) = 0.7071067811865475*m0r[2]+0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  BigAEM(4,7) = 0.7071067811865475*m0r[3]+0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  BigAEM(5,4) = 0.7071067811865475*m0r[1]+0.3535533905932737*m0Sr[1]-0.3535533905932737*cE[1]; 
  BigAEM(5,5) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  BigAEM(6,4) = 0.7071067811865475*m0r[2]+0.3535533905932737*m0Sr[2]-0.3535533905932737*cE[2]; 
  BigAEM(6,6) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
  BigAEM(7,4) = 0.7071067811865475*m0r[3]+0.3535533905932737*m0Sr[3]-0.3535533905932737*cE[3]; 
  BigAEM(7,7) = 0.7071067811865475*m0r[0]+0.3535533905932737*m0Sr[0]-0.3535533905932737*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m2Sr[0],m2Sr[1],m2Sr[2],m2Sr[3]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,4,1) = xEV.segment<4>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = xEV.segment<4>(4); 
 
} 
 
void GkSelfPrimMoments3x2vMax_P2(const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[10]; 
  double m1r[10]; 
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
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(20,20); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(20);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(20);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.3535533905932737*m0r[0]; 
  BigAEM(0,1) = 0.3535533905932737*m0r[1]; 
  BigAEM(0,2) = 0.3535533905932737*m0r[2]; 
  BigAEM(0,3) = 0.3535533905932737*m0r[3]; 
  BigAEM(0,4) = 0.3535533905932737*m0r[4]; 
  BigAEM(0,5) = 0.3535533905932737*m0r[5]; 
  BigAEM(0,6) = 0.3535533905932737*m0r[6]; 
  BigAEM(0,7) = 0.3535533905932737*m0r[7]; 
  BigAEM(0,8) = 0.3535533905932737*m0r[8]; 
  BigAEM(0,9) = 0.3535533905932737*m0r[9]; 
  BigAEM(1,0) = 0.3535533905932737*m0r[1]; 
  BigAEM(1,1) = 0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  BigAEM(1,2) = 0.3535533905932737*m0r[4]; 
  BigAEM(1,3) = 0.3535533905932737*m0r[5]; 
  BigAEM(1,4) = 0.3535533905932737*m0r[2]; 
  BigAEM(1,5) = 0.3535533905932737*m0r[3]; 
  BigAEM(1,7) = 0.3162277660168379*m0r[1]; 
  BigAEM(2,0) = 0.3535533905932737*m0r[2]; 
  BigAEM(2,1) = 0.3535533905932737*m0r[4]; 
  BigAEM(2,2) = 0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  BigAEM(2,3) = 0.3535533905932737*m0r[6]; 
  BigAEM(2,4) = 0.3535533905932737*m0r[1]; 
  BigAEM(2,6) = 0.3535533905932737*m0r[3]; 
  BigAEM(2,8) = 0.3162277660168379*m0r[2]; 
  BigAEM(3,0) = 0.3535533905932737*m0r[3]; 
  BigAEM(3,1) = 0.3535533905932737*m0r[5]; 
  BigAEM(3,2) = 0.3535533905932737*m0r[6]; 
  BigAEM(3,3) = 0.3162277660168379*m0r[9]+0.3535533905932737*m0r[0]; 
  BigAEM(3,5) = 0.3535533905932737*m0r[1]; 
  BigAEM(3,6) = 0.3535533905932737*m0r[2]; 
  BigAEM(3,9) = 0.3162277660168379*m0r[3]; 
  BigAEM(4,0) = 0.3535533905932737*m0r[4]; 
  BigAEM(4,1) = 0.3535533905932737*m0r[2]; 
  BigAEM(4,2) = 0.3535533905932737*m0r[1]; 
  BigAEM(4,4) = 0.3162277660168379*m0r[8]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  BigAEM(4,5) = 0.3535533905932737*m0r[6]; 
  BigAEM(4,6) = 0.3535533905932737*m0r[5]; 
  BigAEM(4,7) = 0.3162277660168379*m0r[4]; 
  BigAEM(4,8) = 0.3162277660168379*m0r[4]; 
  BigAEM(5,0) = 0.3535533905932737*m0r[5]; 
  BigAEM(5,1) = 0.3535533905932737*m0r[3]; 
  BigAEM(5,3) = 0.3535533905932737*m0r[1]; 
  BigAEM(5,4) = 0.3535533905932737*m0r[6]; 
  BigAEM(5,5) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  BigAEM(5,6) = 0.3535533905932737*m0r[4]; 
  BigAEM(5,7) = 0.3162277660168379*m0r[5]; 
  BigAEM(5,9) = 0.3162277660168379*m0r[5]; 
  BigAEM(6,0) = 0.3535533905932737*m0r[6]; 
  BigAEM(6,2) = 0.3535533905932737*m0r[3]; 
  BigAEM(6,3) = 0.3535533905932737*m0r[2]; 
  BigAEM(6,4) = 0.3535533905932737*m0r[5]; 
  BigAEM(6,5) = 0.3535533905932737*m0r[4]; 
  BigAEM(6,6) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  BigAEM(6,8) = 0.3162277660168379*m0r[6]; 
  BigAEM(6,9) = 0.3162277660168379*m0r[6]; 
  BigAEM(7,0) = 0.3535533905932737*m0r[7]; 
  BigAEM(7,1) = 0.3162277660168379*m0r[1]; 
  BigAEM(7,4) = 0.3162277660168379*m0r[4]; 
  BigAEM(7,5) = 0.3162277660168379*m0r[5]; 
  BigAEM(7,7) = 0.2258769757263128*m0r[7]+0.3535533905932737*m0r[0]; 
  BigAEM(8,0) = 0.3535533905932737*m0r[8]; 
  BigAEM(8,2) = 0.3162277660168379*m0r[2]; 
  BigAEM(8,4) = 0.3162277660168379*m0r[4]; 
  BigAEM(8,6) = 0.3162277660168379*m0r[6]; 
  BigAEM(8,8) = 0.2258769757263128*m0r[8]+0.3535533905932737*m0r[0]; 
  BigAEM(9,0) = 0.3535533905932737*m0r[9]; 
  BigAEM(9,3) = 0.3162277660168379*m0r[3]; 
  BigAEM(9,5) = 0.3162277660168379*m0r[5]; 
  BigAEM(9,6) = 0.3162277660168379*m0r[6]; 
  BigAEM(9,9) = 0.2258769757263128*m0r[9]+0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,10) = -0.3535533905932737*cM[0]; 
  BigAEM(0,11) = -0.3535533905932737*cM[1]; 
  BigAEM(0,12) = -0.3535533905932737*cM[2]; 
  BigAEM(0,13) = -0.3535533905932737*cM[3]; 
  BigAEM(0,14) = -0.3535533905932737*cM[4]; 
  BigAEM(0,15) = -0.3535533905932737*cM[5]; 
  BigAEM(0,16) = -0.3535533905932737*cM[6]; 
  BigAEM(0,17) = -0.3535533905932737*cM[7]; 
  BigAEM(0,18) = -0.3535533905932737*cM[8]; 
  BigAEM(0,19) = -0.3535533905932737*cM[9]; 
  BigAEM(1,10) = -0.3535533905932737*cM[1]; 
  BigAEM(1,11) = (-0.3162277660168379*cM[7])-0.3535533905932737*cM[0]; 
  BigAEM(1,12) = -0.3535533905932737*cM[4]; 
  BigAEM(1,13) = -0.3535533905932737*cM[5]; 
  BigAEM(1,14) = -0.3535533905932737*cM[2]; 
  BigAEM(1,15) = -0.3535533905932737*cM[3]; 
  BigAEM(1,17) = -0.3162277660168379*cM[1]; 
  BigAEM(2,10) = -0.3535533905932737*cM[2]; 
  BigAEM(2,11) = -0.3535533905932737*cM[4]; 
  BigAEM(2,12) = (-0.3162277660168379*cM[8])-0.3535533905932737*cM[0]; 
  BigAEM(2,13) = -0.3535533905932737*cM[6]; 
  BigAEM(2,14) = -0.3535533905932737*cM[1]; 
  BigAEM(2,16) = -0.3535533905932737*cM[3]; 
  BigAEM(2,18) = -0.3162277660168379*cM[2]; 
  BigAEM(3,10) = -0.3535533905932737*cM[3]; 
  BigAEM(3,11) = -0.3535533905932737*cM[5]; 
  BigAEM(3,12) = -0.3535533905932737*cM[6]; 
  BigAEM(3,13) = (-0.3162277660168379*cM[9])-0.3535533905932737*cM[0]; 
  BigAEM(3,15) = -0.3535533905932737*cM[1]; 
  BigAEM(3,16) = -0.3535533905932737*cM[2]; 
  BigAEM(3,19) = -0.3162277660168379*cM[3]; 
  BigAEM(4,10) = -0.3535533905932737*cM[4]; 
  BigAEM(4,11) = -0.3535533905932737*cM[2]; 
  BigAEM(4,12) = -0.3535533905932737*cM[1]; 
  BigAEM(4,14) = (-0.3162277660168379*cM[8])-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  BigAEM(4,15) = -0.3535533905932737*cM[6]; 
  BigAEM(4,16) = -0.3535533905932737*cM[5]; 
  BigAEM(4,17) = -0.3162277660168379*cM[4]; 
  BigAEM(4,18) = -0.3162277660168379*cM[4]; 
  BigAEM(5,10) = -0.3535533905932737*cM[5]; 
  BigAEM(5,11) = -0.3535533905932737*cM[3]; 
  BigAEM(5,13) = -0.3535533905932737*cM[1]; 
  BigAEM(5,14) = -0.3535533905932737*cM[6]; 
  BigAEM(5,15) = (-0.3162277660168379*cM[9])-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  BigAEM(5,16) = -0.3535533905932737*cM[4]; 
  BigAEM(5,17) = -0.3162277660168379*cM[5]; 
  BigAEM(5,19) = -0.3162277660168379*cM[5]; 
  BigAEM(6,10) = -0.3535533905932737*cM[6]; 
  BigAEM(6,12) = -0.3535533905932737*cM[3]; 
  BigAEM(6,13) = -0.3535533905932737*cM[2]; 
  BigAEM(6,14) = -0.3535533905932737*cM[5]; 
  BigAEM(6,15) = -0.3535533905932737*cM[4]; 
  BigAEM(6,16) = (-0.3162277660168379*cM[9])-0.3162277660168379*cM[8]-0.3535533905932737*cM[0]; 
  BigAEM(6,18) = -0.3162277660168379*cM[6]; 
  BigAEM(6,19) = -0.3162277660168379*cM[6]; 
  BigAEM(7,10) = -0.3535533905932737*cM[7]; 
  BigAEM(7,11) = -0.3162277660168379*cM[1]; 
  BigAEM(7,14) = -0.3162277660168379*cM[4]; 
  BigAEM(7,15) = -0.3162277660168379*cM[5]; 
  BigAEM(7,17) = (-0.2258769757263128*cM[7])-0.3535533905932737*cM[0]; 
  BigAEM(8,10) = -0.3535533905932737*cM[8]; 
  BigAEM(8,12) = -0.3162277660168379*cM[2]; 
  BigAEM(8,14) = -0.3162277660168379*cM[4]; 
  BigAEM(8,16) = -0.3162277660168379*cM[6]; 
  BigAEM(8,18) = (-0.2258769757263128*cM[8])-0.3535533905932737*cM[0]; 
  BigAEM(9,10) = -0.3535533905932737*cM[9]; 
  BigAEM(9,13) = -0.3162277660168379*cM[3]; 
  BigAEM(9,15) = -0.3162277660168379*cM[5]; 
  BigAEM(9,16) = -0.3162277660168379*cM[6]; 
  BigAEM(9,19) = (-0.2258769757263128*cM[9])-0.3535533905932737*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(10,0) = 0.3535533905932737*m1r[0]; 
  BigAEM(10,1) = 0.3535533905932737*m1r[1]; 
  BigAEM(10,2) = 0.3535533905932737*m1r[2]; 
  BigAEM(10,3) = 0.3535533905932737*m1r[3]; 
  BigAEM(10,4) = 0.3535533905932737*m1r[4]; 
  BigAEM(10,5) = 0.3535533905932737*m1r[5]; 
  BigAEM(10,6) = 0.3535533905932737*m1r[6]; 
  BigAEM(10,7) = 0.3535533905932737*m1r[7]; 
  BigAEM(10,8) = 0.3535533905932737*m1r[8]; 
  BigAEM(10,9) = 0.3535533905932737*m1r[9]; 
  BigAEM(11,0) = 0.3535533905932737*m1r[1]; 
  BigAEM(11,1) = 0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  BigAEM(11,2) = 0.3535533905932737*m1r[4]; 
  BigAEM(11,3) = 0.3535533905932737*m1r[5]; 
  BigAEM(11,4) = 0.3535533905932737*m1r[2]; 
  BigAEM(11,5) = 0.3535533905932737*m1r[3]; 
  BigAEM(11,7) = 0.3162277660168379*m1r[1]; 
  BigAEM(12,0) = 0.3535533905932737*m1r[2]; 
  BigAEM(12,1) = 0.3535533905932737*m1r[4]; 
  BigAEM(12,2) = 0.3162277660168379*m1r[8]+0.3535533905932737*m1r[0]; 
  BigAEM(12,3) = 0.3535533905932737*m1r[6]; 
  BigAEM(12,4) = 0.3535533905932737*m1r[1]; 
  BigAEM(12,6) = 0.3535533905932737*m1r[3]; 
  BigAEM(12,8) = 0.3162277660168379*m1r[2]; 
  BigAEM(13,0) = 0.3535533905932737*m1r[3]; 
  BigAEM(13,1) = 0.3535533905932737*m1r[5]; 
  BigAEM(13,2) = 0.3535533905932737*m1r[6]; 
  BigAEM(13,3) = 0.3162277660168379*m1r[9]+0.3535533905932737*m1r[0]; 
  BigAEM(13,5) = 0.3535533905932737*m1r[1]; 
  BigAEM(13,6) = 0.3535533905932737*m1r[2]; 
  BigAEM(13,9) = 0.3162277660168379*m1r[3]; 
  BigAEM(14,0) = 0.3535533905932737*m1r[4]; 
  BigAEM(14,1) = 0.3535533905932737*m1r[2]; 
  BigAEM(14,2) = 0.3535533905932737*m1r[1]; 
  BigAEM(14,4) = 0.3162277660168379*m1r[8]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  BigAEM(14,5) = 0.3535533905932737*m1r[6]; 
  BigAEM(14,6) = 0.3535533905932737*m1r[5]; 
  BigAEM(14,7) = 0.3162277660168379*m1r[4]; 
  BigAEM(14,8) = 0.3162277660168379*m1r[4]; 
  BigAEM(15,0) = 0.3535533905932737*m1r[5]; 
  BigAEM(15,1) = 0.3535533905932737*m1r[3]; 
  BigAEM(15,3) = 0.3535533905932737*m1r[1]; 
  BigAEM(15,4) = 0.3535533905932737*m1r[6]; 
  BigAEM(15,5) = 0.3162277660168379*m1r[9]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  BigAEM(15,6) = 0.3535533905932737*m1r[4]; 
  BigAEM(15,7) = 0.3162277660168379*m1r[5]; 
  BigAEM(15,9) = 0.3162277660168379*m1r[5]; 
  BigAEM(16,0) = 0.3535533905932737*m1r[6]; 
  BigAEM(16,2) = 0.3535533905932737*m1r[3]; 
  BigAEM(16,3) = 0.3535533905932737*m1r[2]; 
  BigAEM(16,4) = 0.3535533905932737*m1r[5]; 
  BigAEM(16,5) = 0.3535533905932737*m1r[4]; 
  BigAEM(16,6) = 0.3162277660168379*m1r[9]+0.3162277660168379*m1r[8]+0.3535533905932737*m1r[0]; 
  BigAEM(16,8) = 0.3162277660168379*m1r[6]; 
  BigAEM(16,9) = 0.3162277660168379*m1r[6]; 
  BigAEM(17,0) = 0.3535533905932737*m1r[7]; 
  BigAEM(17,1) = 0.3162277660168379*m1r[1]; 
  BigAEM(17,4) = 0.3162277660168379*m1r[4]; 
  BigAEM(17,5) = 0.3162277660168379*m1r[5]; 
  BigAEM(17,7) = 0.2258769757263128*m1r[7]+0.3535533905932737*m1r[0]; 
  BigAEM(18,0) = 0.3535533905932737*m1r[8]; 
  BigAEM(18,2) = 0.3162277660168379*m1r[2]; 
  BigAEM(18,4) = 0.3162277660168379*m1r[4]; 
  BigAEM(18,6) = 0.3162277660168379*m1r[6]; 
  BigAEM(18,8) = 0.2258769757263128*m1r[8]+0.3535533905932737*m1r[0]; 
  BigAEM(19,0) = 0.3535533905932737*m1r[9]; 
  BigAEM(19,3) = 0.3162277660168379*m1r[3]; 
  BigAEM(19,5) = 0.3162277660168379*m1r[5]; 
  BigAEM(19,6) = 0.3162277660168379*m1r[6]; 
  BigAEM(19,9) = 0.2258769757263128*m1r[9]+0.3535533905932737*m1r[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(10,10) = 1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  BigAEM(10,11) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  BigAEM(10,12) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  BigAEM(10,13) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  BigAEM(10,14) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  BigAEM(10,15) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  BigAEM(10,16) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  BigAEM(10,17) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  BigAEM(10,18) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  BigAEM(10,19) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  BigAEM(11,10) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  BigAEM(11,11) = 0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  BigAEM(11,12) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  BigAEM(11,13) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  BigAEM(11,14) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  BigAEM(11,15) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  BigAEM(11,17) = 0.9486832980505137*m0r[1]-0.3162277660168379*cE[1]; 
  BigAEM(12,10) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  BigAEM(12,11) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  BigAEM(12,12) = 0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  BigAEM(12,13) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  BigAEM(12,14) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  BigAEM(12,16) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  BigAEM(12,18) = 0.9486832980505137*m0r[2]-0.3162277660168379*cE[2]; 
  BigAEM(13,10) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  BigAEM(13,11) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  BigAEM(13,12) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  BigAEM(13,13) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  BigAEM(13,15) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  BigAEM(13,16) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  BigAEM(13,19) = 0.9486832980505137*m0r[3]-0.3162277660168379*cE[3]; 
  BigAEM(14,10) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  BigAEM(14,11) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  BigAEM(14,12) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  BigAEM(14,14) = 0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  BigAEM(14,15) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  BigAEM(14,16) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  BigAEM(14,17) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  BigAEM(14,18) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  BigAEM(15,10) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  BigAEM(15,11) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  BigAEM(15,13) = 1.060660171779821*m0r[1]-0.3535533905932737*cE[1]; 
  BigAEM(15,14) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  BigAEM(15,15) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.9486832980505137*m0r[7]-0.3162277660168379*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  BigAEM(15,16) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  BigAEM(15,17) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  BigAEM(15,19) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  BigAEM(16,10) = 1.060660171779821*m0r[6]-0.3535533905932737*cE[6]; 
  BigAEM(16,12) = 1.060660171779821*m0r[3]-0.3535533905932737*cE[3]; 
  BigAEM(16,13) = 1.060660171779821*m0r[2]-0.3535533905932737*cE[2]; 
  BigAEM(16,14) = 1.060660171779821*m0r[5]-0.3535533905932737*cE[5]; 
  BigAEM(16,15) = 1.060660171779821*m0r[4]-0.3535533905932737*cE[4]; 
  BigAEM(16,16) = 0.9486832980505137*m0r[9]-0.3162277660168379*cE[9]+0.9486832980505137*m0r[8]-0.3162277660168379*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  BigAEM(16,18) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  BigAEM(16,19) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  BigAEM(17,10) = 1.060660171779821*m0r[7]-0.3535533905932737*cE[7]; 
  BigAEM(17,11) = 0.9486832980505137*m0r[1]-0.3162277660168379*cE[1]; 
  BigAEM(17,14) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  BigAEM(17,15) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  BigAEM(17,17) = 0.6776309271789384*m0r[7]-0.2258769757263128*cE[7]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  BigAEM(18,10) = 1.060660171779821*m0r[8]-0.3535533905932737*cE[8]; 
  BigAEM(18,12) = 0.9486832980505137*m0r[2]-0.3162277660168379*cE[2]; 
  BigAEM(18,14) = 0.9486832980505137*m0r[4]-0.3162277660168379*cE[4]; 
  BigAEM(18,16) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  BigAEM(18,18) = 0.6776309271789384*m0r[8]-0.2258769757263128*cE[8]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
  BigAEM(19,10) = 1.060660171779821*m0r[9]-0.3535533905932737*cE[9]; 
  BigAEM(19,13) = 0.9486832980505137*m0r[3]-0.3162277660168379*cE[3]; 
  BigAEM(19,15) = 0.9486832980505137*m0r[5]-0.3162277660168379*cE[5]; 
  BigAEM(19,16) = 0.9486832980505137*m0r[6]-0.3162277660168379*cE[6]; 
  BigAEM(19,19) = 0.6776309271789384*m0r[9]-0.2258769757263128*cE[9]+1.060660171779821*m0r[0]-0.3535533905932737*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5],m2r[6],m2r[7],m2r[8],m2r[9]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,10,1) = xEV.segment<10>(0); 
 
  Eigen::Map<VectorXd>(vtSq,10,1) = xEV.segment<10>(10); 
 
} 
 
