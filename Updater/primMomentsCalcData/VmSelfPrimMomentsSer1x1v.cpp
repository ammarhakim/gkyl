
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
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
void VmBoundaryIntegral1x1vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[4]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0; 
 
  if (atLower) {
 
    out[0] += (1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x1vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0; 
 
  if (atLower) {
 
    out[0] += ((-1.58113883008419*fIn[5])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[1] += ((-1.58113883008419*fIn[7])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS; 
 
  } else {
 
    out[0] += (1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x1vSer_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[12]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0; 
 
  if (atLower) {
 
    out[0] += (1.870828693386971*fIn[9]-1.58113883008419*fIn[5]+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.870828693386971*fIn[11]-1.58113883008419*fIn[7]+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS; 
    out[3] += (1.224744871391589*fIn[10]-0.7071067811865475*fIn[8])*dS; 
 
  } else {
 
    out[0] += (1.870828693386971*fIn[9]+1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.870828693386971*fIn[11]+1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS; 
    out[3] += (1.224744871391589*fIn[10]+0.7071067811865475*fIn[8])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x1vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[4]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0; 
 
  if (atLower) {
 
    out[0] += (1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS*vBoundary+(0.6123724356957944*dxv[1]*fIn[2]-0.3535533905932737*fIn[0]*dxv[1])*dS; 
    out[1] += (1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS*vBoundary+(0.6123724356957944*dxv[1]*fIn[3]-0.3535533905932737*dxv[1]*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS*vBoundary+((-0.6123724356957944*dxv[1]*fIn[2])-0.3535533905932737*fIn[0]*dxv[1])*dS; 
    out[1] += (1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS*vBoundary+((-0.6123724356957944*dxv[1]*fIn[3])-0.3535533905932737*dxv[1]*fIn[1])*dS; 
 
  }
 
} 
 
void VmBoundaryIntegral1x1vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0; 
 
  if (atLower) {
 
    out[0] += ((-1.58113883008419*fIn[5])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += ((-1.58113883008419*fIn[7])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS*vBoundary; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += (1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS*vBoundary; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS*vBoundary; 
 
  }
 
} 
 
void VmBoundaryIntegral1x1vSer_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[12]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0; 
 
  if (atLower) {
 
    out[0] += (1.870828693386971*fIn[9]-1.58113883008419*fIn[5]+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += (1.870828693386971*fIn[11]-1.58113883008419*fIn[7]+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS*vBoundary; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS*vBoundary; 
    out[3] += (1.224744871391589*fIn[10]-0.7071067811865475*fIn[8])*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.870828693386971*fIn[9]+1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += (1.870828693386971*fIn[11]+1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS*vBoundary; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS*vBoundary; 
    out[3] += (1.224744871391589*fIn[10]+0.7071067811865475*fIn[8])*dS*vBoundary; 
 
  }
 
} 
 
void VmSelfPrimMoments1x1vSer_P1(binOpData_t* data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
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
  double m1r[2]; 
  double m0Sr[2]; 
  double m1Sr[2]; 
  double m2Sr[2]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(4,4); 
  
  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,0) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,1) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,2) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,3) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,2) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,3) = -0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(2,0) = 0.7071067811865475*m1Sr[0]; 
  data->AEM_S(2,1) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(3,0) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(3,1) = 0.7071067811865475*m1Sr[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(2,2) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(2,3) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(3,2) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(3,3) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m2Sr[0],m2Sr[1]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,2,1) = data->u_S.segment<2>(0); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = data->u_S.segment<2>(2); 
 
} 
 
void VmSelfPrimMoments1x1vSer_P2(binOpData_t* data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
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
  double m1r[3]; 
  double m2r[3]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(6,6); 
  
  
 
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
  data->AEM_S(0,3) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,4) = -0.7071067811865475*cM[1]; 
  data->AEM_S(0,5) = -0.7071067811865475*cM[2]; 
  data->AEM_S(1,3) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,4) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(1,5) = -0.6324555320336759*cM[1]; 
  data->AEM_S(2,3) = -0.7071067811865475*cM[2]; 
  data->AEM_S(2,4) = -0.6324555320336759*cM[1]; 
  data->AEM_S(2,5) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(3,0) = 0.7071067811865475*m1r[0]; 
  data->AEM_S(3,1) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(3,2) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(4,0) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(4,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(4,2) = 0.6324555320336759*m1r[1]; 
  data->AEM_S(5,0) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(5,1) = 0.6324555320336759*m1r[1]; 
  data->AEM_S(5,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(3,3) = 0.7071067811865475*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(3,4) = 0.7071067811865475*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(3,5) = 0.7071067811865475*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(4,3) = 0.7071067811865475*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(4,4) = 0.6324555320336759*m0r[2]-0.6324555320336759*cE[2]+0.7071067811865475*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(4,5) = 0.6324555320336759*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(5,3) = 0.7071067811865475*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(5,4) = 0.6324555320336759*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(5,5) = 0.4517539514526256*m0r[2]-0.4517539514526256*cE[2]+0.7071067811865475*m0r[0]-0.7071067811865475*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m2r[0],m2r[1],m2r[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,3,1) = data->u_S.segment<3>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = data->u_S.segment<3>(3); 
 
} 
 
void VmSelfPrimMoments1x1vSer_P3(binOpData_t* data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
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
  double m1r[4]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(8,8); 
  
  
 
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
  data->AEM_S(0,4) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,5) = -0.7071067811865475*cM[1]; 
  data->AEM_S(0,6) = -0.7071067811865475*cM[2]; 
  data->AEM_S(0,7) = -0.7071067811865475*cM[3]; 
  data->AEM_S(1,4) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,5) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(1,6) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  data->AEM_S(1,7) = -0.6210590034081186*cM[2]; 
  data->AEM_S(2,4) = -0.7071067811865475*cM[2]; 
  data->AEM_S(2,5) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  data->AEM_S(2,6) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(2,7) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  data->AEM_S(3,4) = -0.7071067811865475*cM[3]; 
  data->AEM_S(3,5) = -0.6210590034081186*cM[2]; 
  data->AEM_S(3,6) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  data->AEM_S(3,7) = (-0.421637021355784*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(4,0) = 0.7071067811865475*m1r[0]; 
  data->AEM_S(4,1) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(4,2) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(4,3) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(5,0) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(5,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(5,2) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  data->AEM_S(5,3) = 0.6210590034081186*m1r[2]; 
  data->AEM_S(6,0) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(6,1) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  data->AEM_S(6,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(6,3) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  data->AEM_S(7,0) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(7,1) = 0.6210590034081186*m1r[2]; 
  data->AEM_S(7,2) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  data->AEM_S(7,3) = 0.421637021355784*m1r[2]+0.7071067811865475*m1r[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(4,4) = 0.7071067811865475*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(4,5) = 0.7071067811865475*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(4,6) = 0.7071067811865475*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(4,7) = 0.7071067811865475*m0r[3]-0.7071067811865475*cE[3]; 
  data->AEM_S(5,4) = 0.7071067811865475*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(5,5) = 0.6324555320336759*m0r[2]-0.6324555320336759*cE[2]+0.7071067811865475*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(5,6) = 0.6210590034081186*m0r[3]-0.6210590034081186*cE[3]+0.6324555320336759*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(5,7) = 0.6210590034081186*m0r[2]-0.6210590034081186*cE[2]; 
  data->AEM_S(6,4) = 0.7071067811865475*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(6,5) = 0.6210590034081186*m0r[3]-0.6210590034081186*cE[3]+0.6324555320336759*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(6,6) = 0.4517539514526256*m0r[2]-0.4517539514526256*cE[2]+0.7071067811865475*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(6,7) = 0.421637021355784*m0r[3]-0.421637021355784*cE[3]+0.6210590034081186*m0r[1]-0.6210590034081186*cE[1]; 
  data->AEM_S(7,4) = 0.7071067811865475*m0r[3]-0.7071067811865475*cE[3]; 
  data->AEM_S(7,5) = 0.6210590034081186*m0r[2]-0.6210590034081186*cE[2]; 
  data->AEM_S(7,6) = 0.421637021355784*m0r[3]-0.421637021355784*cE[3]+0.6210590034081186*m0r[1]-0.6210590034081186*cE[1]; 
  data->AEM_S(7,7) = 0.421637021355784*m0r[2]-0.421637021355784*cE[2]+0.7071067811865475*m0r[0]-0.7071067811865475*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m2r[0],m2r[1],m2r[2],m2r[3]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,4,1) = data->u_S.segment<4>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = data->u_S.segment<4>(4); 
 
} 
 
