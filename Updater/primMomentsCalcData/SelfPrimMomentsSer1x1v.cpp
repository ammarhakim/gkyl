#include <math.h> 
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void SelfPrimMoments1x1vSer_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE: vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(4,4); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(4);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(4);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.7071067811865475*m0[0]; 
  BigAEM(0,1) = 0.7071067811865475*m0[1]; 
  BigAEM(1,0) = 0.7071067811865475*m0[1]; 
  BigAEM(1,1) = 0.7071067811865475*m0[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,2) = -0.7071067811865475*cM[0]; 
  BigAEM(0,3) = -0.7071067811865475*cM[1]; 
  BigAEM(1,2) = -0.7071067811865475*cM[1]; 
  BigAEM(1,3) = -0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(2,0) = 0.7071067811865475*m1[0]; 
  BigAEM(2,1) = 0.7071067811865475*m1[1]; 
  BigAEM(3,0) = 0.7071067811865475*m1[1]; 
  BigAEM(3,1) = 0.7071067811865475*m1[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(2,2) = 0.7071067811865475*m0[0]*pVdim-0.7071067811865475*cE[0]; 
  BigAEM(2,3) = 0.7071067811865475*m0[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(3,2) = 0.7071067811865475*m0[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(3,3) = 0.7071067811865475*m0[0]*pVdim-0.7071067811865475*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1[0],m1[1],m2[0],m2[1]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,2,1) = xEV.segment<2>(0); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = xEV.segment<2>(2); 
 
} 
 
void SelfPrimMoments1x1vSer_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE: vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(6,6); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(6);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(6);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.7071067811865475*m0[0]; 
  BigAEM(0,1) = 0.7071067811865475*m0[1]; 
  BigAEM(0,2) = 0.7071067811865475*m0[2]; 
  BigAEM(1,0) = 0.7071067811865475*m0[1]; 
  BigAEM(1,1) = 0.6324555320336759*m0[2]+0.7071067811865475*m0[0]; 
  BigAEM(1,2) = 0.6324555320336759*m0[1]; 
  BigAEM(2,0) = 0.7071067811865475*m0[2]; 
  BigAEM(2,1) = 0.6324555320336759*m0[1]; 
  BigAEM(2,2) = 0.4517539514526256*m0[2]+0.7071067811865475*m0[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,3) = -0.7071067811865475*cM[0]; 
  BigAEM(0,4) = -0.7071067811865475*cM[1]; 
  BigAEM(0,5) = -0.7071067811865475*cM[2]; 
  BigAEM(1,3) = -0.7071067811865475*cM[1]; 
  BigAEM(1,4) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  BigAEM(1,5) = -0.6324555320336759*cM[1]; 
  BigAEM(2,3) = -0.7071067811865475*cM[2]; 
  BigAEM(2,4) = -0.6324555320336759*cM[1]; 
  BigAEM(2,5) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(3,0) = 0.7071067811865475*m1[0]; 
  BigAEM(3,1) = 0.7071067811865475*m1[1]; 
  BigAEM(3,2) = 0.7071067811865475*m1[2]; 
  BigAEM(4,0) = 0.7071067811865475*m1[1]; 
  BigAEM(4,1) = 0.6324555320336759*m1[2]+0.7071067811865475*m1[0]; 
  BigAEM(4,2) = 0.6324555320336759*m1[1]; 
  BigAEM(5,0) = 0.7071067811865475*m1[2]; 
  BigAEM(5,1) = 0.6324555320336759*m1[1]; 
  BigAEM(5,2) = 0.4517539514526256*m1[2]+0.7071067811865475*m1[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(3,3) = 0.7071067811865475*m0[0]*pVdim-0.7071067811865475*cE[0]; 
  BigAEM(3,4) = 0.7071067811865475*m0[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(3,5) = 0.7071067811865475*m0[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(4,3) = 0.7071067811865475*m0[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(4,4) = 0.6324555320336759*m0[2]*pVdim+0.7071067811865475*m0[0]*pVdim-0.6324555320336759*cE[2]-0.7071067811865475*cE[0]; 
  BigAEM(4,5) = 0.6324555320336759*m0[1]*pVdim-0.6324555320336759*cE[1]; 
  BigAEM(5,3) = 0.7071067811865475*m0[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(5,4) = 0.6324555320336759*m0[1]*pVdim-0.6324555320336759*cE[1]; 
  BigAEM(5,5) = 0.4517539514526256*m0[2]*pVdim+0.7071067811865475*m0[0]*pVdim-0.4517539514526256*cE[2]-0.7071067811865475*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1[0],m1[1],m1[2],m2[0],m2[1],m2[2]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,3,1) = xEV.segment<3>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = xEV.segment<3>(3); 
 
} 
 
void SelfPrimMoments1x1vSer_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE: vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(8,8); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(8);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(8);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.7071067811865475*m0[0]; 
  BigAEM(0,1) = 0.7071067811865475*m0[1]; 
  BigAEM(0,2) = 0.7071067811865475*m0[2]; 
  BigAEM(0,3) = 0.7071067811865475*m0[3]; 
  BigAEM(1,0) = 0.7071067811865475*m0[1]; 
  BigAEM(1,1) = 0.6324555320336759*m0[2]+0.7071067811865475*m0[0]; 
  BigAEM(1,2) = 0.6210590034081186*m0[3]+0.6324555320336759*m0[1]; 
  BigAEM(1,3) = 0.6210590034081186*m0[2]; 
  BigAEM(2,0) = 0.7071067811865475*m0[2]; 
  BigAEM(2,1) = 0.6210590034081186*m0[3]+0.6324555320336759*m0[1]; 
  BigAEM(2,2) = 0.4517539514526256*m0[2]+0.7071067811865475*m0[0]; 
  BigAEM(2,3) = 0.421637021355784*m0[3]+0.6210590034081186*m0[1]; 
  BigAEM(3,0) = 0.7071067811865475*m0[3]; 
  BigAEM(3,1) = 0.6210590034081186*m0[2]; 
  BigAEM(3,2) = 0.421637021355784*m0[3]+0.6210590034081186*m0[1]; 
  BigAEM(3,3) = 0.421637021355784*m0[2]+0.7071067811865475*m0[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,4) = -0.7071067811865475*cM[0]; 
  BigAEM(0,5) = -0.7071067811865475*cM[1]; 
  BigAEM(0,6) = -0.7071067811865475*cM[2]; 
  BigAEM(0,7) = -0.7071067811865475*cM[3]; 
  BigAEM(1,4) = -0.7071067811865475*cM[1]; 
  BigAEM(1,5) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  BigAEM(1,6) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  BigAEM(1,7) = -0.6210590034081186*cM[2]; 
  BigAEM(2,4) = -0.7071067811865475*cM[2]; 
  BigAEM(2,5) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  BigAEM(2,6) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
  BigAEM(2,7) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  BigAEM(3,4) = -0.7071067811865475*cM[3]; 
  BigAEM(3,5) = -0.6210590034081186*cM[2]; 
  BigAEM(3,6) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  BigAEM(3,7) = (-0.421637021355784*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(4,0) = 0.7071067811865475*m1[0]; 
  BigAEM(4,1) = 0.7071067811865475*m1[1]; 
  BigAEM(4,2) = 0.7071067811865475*m1[2]; 
  BigAEM(4,3) = 0.7071067811865475*m1[3]; 
  BigAEM(5,0) = 0.7071067811865475*m1[1]; 
  BigAEM(5,1) = 0.6324555320336759*m1[2]+0.7071067811865475*m1[0]; 
  BigAEM(5,2) = 0.6210590034081186*m1[3]+0.6324555320336759*m1[1]; 
  BigAEM(5,3) = 0.6210590034081186*m1[2]; 
  BigAEM(6,0) = 0.7071067811865475*m1[2]; 
  BigAEM(6,1) = 0.6210590034081186*m1[3]+0.6324555320336759*m1[1]; 
  BigAEM(6,2) = 0.4517539514526256*m1[2]+0.7071067811865475*m1[0]; 
  BigAEM(6,3) = 0.421637021355784*m1[3]+0.6210590034081186*m1[1]; 
  BigAEM(7,0) = 0.7071067811865475*m1[3]; 
  BigAEM(7,1) = 0.6210590034081186*m1[2]; 
  BigAEM(7,2) = 0.421637021355784*m1[3]+0.6210590034081186*m1[1]; 
  BigAEM(7,3) = 0.421637021355784*m1[2]+0.7071067811865475*m1[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(4,4) = 0.7071067811865475*m0[0]*pVdim-0.7071067811865475*cE[0]; 
  BigAEM(4,5) = 0.7071067811865475*m0[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(4,6) = 0.7071067811865475*m0[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(4,7) = 0.7071067811865475*m0[3]*pVdim-0.7071067811865475*cE[3]; 
  BigAEM(5,4) = 0.7071067811865475*m0[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(5,5) = 0.6324555320336759*m0[2]*pVdim+0.7071067811865475*m0[0]*pVdim-0.6324555320336759*cE[2]-0.7071067811865475*cE[0]; 
  BigAEM(5,6) = 0.6210590034081186*m0[3]*pVdim+0.6324555320336759*m0[1]*pVdim-0.6210590034081186*cE[3]-0.6324555320336759*cE[1]; 
  BigAEM(5,7) = 0.6210590034081186*m0[2]*pVdim-0.6210590034081186*cE[2]; 
  BigAEM(6,4) = 0.7071067811865475*m0[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(6,5) = 0.6210590034081186*m0[3]*pVdim+0.6324555320336759*m0[1]*pVdim-0.6210590034081186*cE[3]-0.6324555320336759*cE[1]; 
  BigAEM(6,6) = 0.4517539514526256*m0[2]*pVdim+0.7071067811865475*m0[0]*pVdim-0.4517539514526256*cE[2]-0.7071067811865475*cE[0]; 
  BigAEM(6,7) = 0.421637021355784*m0[3]*pVdim+0.6210590034081186*m0[1]*pVdim-0.421637021355784*cE[3]-0.6210590034081186*cE[1]; 
  BigAEM(7,4) = 0.7071067811865475*m0[3]*pVdim-0.7071067811865475*cE[3]; 
  BigAEM(7,5) = 0.6210590034081186*m0[2]*pVdim-0.6210590034081186*cE[2]; 
  BigAEM(7,6) = 0.421637021355784*m0[3]*pVdim+0.6210590034081186*m0[1]*pVdim-0.421637021355784*cE[3]-0.6210590034081186*cE[1]; 
  BigAEM(7,7) = 0.421637021355784*m0[2]*pVdim+0.7071067811865475*m0[0]*pVdim-0.421637021355784*cE[2]-0.7071067811865475*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1[0],m1[1],m1[2],m1[3],m2[0],m2[1],m2[2],m2[3]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,4,1) = xEV.segment<4>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = xEV.segment<4>(4); 
 
} 
 
void BoundaryIntegral1x1vSer_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[2]:             cell length in each direciton. 
  // fvmax[4], fvmin[4]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*1.0; 
 
  out[0] += (1.224744871391589*fvmin[2]*dS+1.224744871391589*fvmax[2]*dS-0.7071067811865475*fvmin[0]*dS+0.7071067811865475*fvmax[0]*dS)*intFac; 
  out[1] += (1.224744871391589*fvmin[3]*dS+1.224744871391589*fvmax[3]*dS-0.7071067811865475*fvmin[1]*dS+0.7071067811865475*fvmax[1]*dS)*intFac; 
 
} 
 
void BoundaryIntegral1x1vSer_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[2]:             cell length in each direciton. 
  // fvmax[8], fvmin[8]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*1.0; 
 
  out[0] += ((-1.58113883008419*fvmin[5]*dS)+1.58113883008419*fvmax[5]*dS+1.224744871391589*fvmin[2]*dS+1.224744871391589*fvmax[2]*dS-0.7071067811865475*fvmin[0]*dS+0.7071067811865475*fvmax[0]*dS)*intFac; 
  out[1] += ((-1.58113883008419*fvmin[7]*dS)+1.58113883008419*fvmax[7]*dS+1.224744871391589*fvmin[3]*dS+1.224744871391589*fvmax[3]*dS-0.7071067811865475*fvmin[1]*dS+0.7071067811865475*fvmax[1]*dS)*intFac; 
  out[2] += (1.224744871391589*fvmin[6]*dS+1.224744871391589*fvmax[6]*dS-0.7071067811865475*fvmin[4]*dS+0.7071067811865475*fvmax[4]*dS)*intFac; 
 
} 
 
void BoundaryIntegral1x1vSer_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[2]:             cell length in each direciton. 
  // fvmax[12], fvmin[12]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*1.0; 
 
  out[0] += (1.870828693386971*fvmin[9]*dS+1.870828693386971*fvmax[9]*dS-1.58113883008419*fvmin[5]*dS+1.58113883008419*fvmax[5]*dS+1.224744871391589*fvmin[2]*dS+1.224744871391589*fvmax[2]*dS-0.7071067811865475*fvmin[0]*dS+0.7071067811865475*fvmax[0]*dS)*intFac; 
  out[1] += (1.870828693386971*fvmin[11]*dS+1.870828693386971*fvmax[11]*dS-1.58113883008419*fvmin[7]*dS+1.58113883008419*fvmax[7]*dS+1.224744871391589*fvmin[3]*dS+1.224744871391589*fvmax[3]*dS-0.7071067811865475*fvmin[1]*dS+0.7071067811865475*fvmax[1]*dS)*intFac; 
  out[2] += (1.224744871391589*fvmin[6]*dS+1.224744871391589*fvmax[6]*dS-0.7071067811865475*fvmin[4]*dS+0.7071067811865475*fvmax[4]*dS)*intFac; 
  out[3] += (1.224744871391589*fvmin[10]*dS+1.224744871391589*fvmax[10]*dS-0.7071067811865475*fvmin[8]*dS+0.7071067811865475*fvmax[8]*dS)*intFac; 
 
} 
 
void BoundaryIntegral1x1vSer_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[2]:             cell length in each direciton. 
  // fvmax[4], fvmin[4]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*1.0; 
 
  out[0] += intFac*(1.224744871391589*fvmin[2]*dS*vmin-0.7071067811865475*fvmin[0]*dS*vmin+1.224744871391589*fvmax[2]*dS*vmax+0.7071067811865475*fvmax[0]*dS*vmax); 
  out[1] += intFac*(1.224744871391589*fvmin[3]*dS*vmin-0.7071067811865475*fvmin[1]*dS*vmin+1.224744871391589*fvmax[3]*dS*vmax+0.7071067811865475*fvmax[1]*dS*vmax); 
 
} 
 
void BoundaryIntegral1x1vSer_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[2]:             cell length in each direciton. 
  // fvmax[8], fvmin[8]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*1.0; 
 
  out[0] += intFac*((-1.58113883008419*fvmin[5]*dS*vmin)+1.224744871391589*fvmin[2]*dS*vmin-0.7071067811865475*fvmin[0]*dS*vmin+1.58113883008419*fvmax[5]*dS*vmax+1.224744871391589*fvmax[2]*dS*vmax+0.7071067811865475*fvmax[0]*dS*vmax); 
  out[1] += intFac*((-1.58113883008419*fvmin[7]*dS*vmin)+1.224744871391589*fvmin[3]*dS*vmin-0.7071067811865475*fvmin[1]*dS*vmin+1.58113883008419*fvmax[7]*dS*vmax+1.224744871391589*fvmax[3]*dS*vmax+0.7071067811865475*fvmax[1]*dS*vmax); 
  out[2] += intFac*(1.224744871391589*fvmin[6]*dS*vmin-0.7071067811865475*fvmin[4]*dS*vmin+1.224744871391589*fvmax[6]*dS*vmax+0.7071067811865475*fvmax[4]*dS*vmax); 
 
} 
 
void BoundaryIntegral1x1vSer_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[2]:             cell length in each direciton. 
  // fvmax[12], fvmin[12]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*1.0; 
 
  out[0] += intFac*(1.870828693386971*fvmin[9]*dS*vmin-1.58113883008419*fvmin[5]*dS*vmin+1.224744871391589*fvmin[2]*dS*vmin-0.7071067811865475*fvmin[0]*dS*vmin+1.870828693386971*fvmax[9]*dS*vmax+1.58113883008419*fvmax[5]*dS*vmax+1.224744871391589*fvmax[2]*dS*vmax+0.7071067811865475*fvmax[0]*dS*vmax); 
  out[1] += intFac*(1.870828693386971*fvmin[11]*dS*vmin-1.58113883008419*fvmin[7]*dS*vmin+1.224744871391589*fvmin[3]*dS*vmin-0.7071067811865475*fvmin[1]*dS*vmin+1.870828693386971*fvmax[11]*dS*vmax+1.58113883008419*fvmax[7]*dS*vmax+1.224744871391589*fvmax[3]*dS*vmax+0.7071067811865475*fvmax[1]*dS*vmax); 
  out[2] += intFac*(1.224744871391589*fvmin[6]*dS*vmin-0.7071067811865475*fvmin[4]*dS*vmin+1.224744871391589*fvmax[6]*dS*vmax+0.7071067811865475*fvmax[4]*dS*vmax); 
  out[3] += intFac*(1.224744871391589*fvmin[10]*dS*vmin-0.7071067811865475*fvmin[8]*dS*vmin+1.224744871391589*fvmax[10]*dS*vmax+0.7071067811865475*fvmax[8]*dS*vmax); 
 
} 
 
