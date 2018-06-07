#include <math.h> 
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void CrossPrimMoments_VmeiLBO_1x2vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uei, double *vtSqei) 
{ 
  // mRat:           mass ratio = m_i/m_e. 
  // ne, ue, vtSqe:  electron number density, bulk flow velocity and T_e/m_e. 
  // ni, ui, vtSqi:  ion number density, bulk flow velocity and T_e/m_e. 
  // u_ei:           bulk flow velocity for electron-ion collision term in electron equation. 
  // vtSq_ei:        squared thermal speed, T_ei/m_e, for electron-ion collision term in electron equation. 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(2,2); 
  Eigen::VectorXd bEV(2); 
  Eigen::VectorXd xEV(2); 
 
  // ....... Compute density ratio through weak division ni/ne .......... // 
  double nRat[2]; 
 
  AEM(0,0) = 0.7071067811865475*ne[0]; 
  AEM(0,1) = 0.7071067811865475*ne[1]; 
  AEM(1,0) = 0.7071067811865475*ne[1]; 
  AEM(1,1) = 0.7071067811865475*ne[0]; 
 
  bEV << ni[0],ni[1]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(nRat,2,1) = xEV; 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    uei[0] = 0.7071067811865475*nRat[1]*ui[a0+1]-0.7071067811865475*nRat[1]*ue[a0+1]+0.7071067811865475*nRat[0]*ui[a0]-0.7071067811865475*nRat[0]*ue[a0]+ue[a0]; 
    uei[1] = 0.7071067811865475*nRat[0]*ui[a0+1]-0.7071067811865475*nRat[0]*ue[a0+1]+ue[a0+1]+0.7071067811865475*nRat[1]*ui[a0]-0.7071067811865475*nRat[1]*ue[a0]; 
 
  } 
 
  // ..... Get the relative speed squared (ue-ui)^2 ..... // 
  double uRelSq[2]; 
  for (unsigned short int k=0; k<2; k++) 
  { 
    uRelSq[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    uRelSq[0] += 0.7071067811865475*ui[a0+1]*ui[a0+1]-1.414213562373095*ue[a0+1]*ui[a0+1]+0.7071067811865475*ue[a0+1]*ue[a0+1]+0.7071067811865475*ui[a0]*ui[a0]-1.414213562373095*ue[a0]*ui[a0]+0.7071067811865475*ue[a0]*ue[a0]; 
    uRelSq[1] += 1.414213562373095*ui[a0]*ui[a0+1]-1.414213562373095*ue[a0]*ui[a0+1]-1.414213562373095*ui[a0]*ue[a0+1]+1.414213562373095*ue[a0]*ue[a0+1]; 
  } 
 
  double rmRat = 1.0/mRat; 
  vtSqei[0] = ((-1.414213562373095*nRat[1]*vtSqe[1])-1.414213562373095*nRat[0]*vtSqe[0])*rmRat+1.414213562373095*nRat[1]*vtSqi[1]+0.3535533905932737*nRat[1]*uRelSq[1]+1.414213562373095*nRat[0]*vtSqi[0]+vtSqe[0]+0.3535533905932737*nRat[0]*uRelSq[0]; 
  vtSqei[1] = ((-1.414213562373095*nRat[0]*vtSqe[1])-1.414213562373095*vtSqe[0]*nRat[1])*rmRat+1.414213562373095*nRat[0]*vtSqi[1]+vtSqe[1]+0.3535533905932737*nRat[0]*uRelSq[1]+1.414213562373095*vtSqi[0]*nRat[1]+0.3535533905932737*uRelSq[0]*nRat[1]; 
 
} 
 
void CrossPrimMoments_VmeiLBO_1x2vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uei, double *vtSqei) 
{ 
  // mRat:           mass ratio = m_i/m_e. 
  // ne, ue, vtSqe:  electron number density, bulk flow velocity and T_e/m_e. 
  // ni, ui, vtSqi:  ion number density, bulk flow velocity and T_e/m_e. 
  // u_ei:           bulk flow velocity for electron-ion collision term in electron equation. 
  // vtSq_ei:        squared thermal speed, T_ei/m_e, for electron-ion collision term in electron equation. 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(3,3); 
  Eigen::VectorXd bEV(3); 
  Eigen::VectorXd xEV(3); 
 
  // ....... Compute density ratio through weak division ni/ne .......... // 
  double nRat[3]; 
 
  AEM(0,0) = 0.7071067811865475*ne[0]; 
  AEM(0,1) = 0.7071067811865475*ne[1]; 
  AEM(0,2) = 0.7071067811865475*ne[2]; 
  AEM(1,0) = 0.7071067811865475*ne[1]; 
  AEM(1,1) = 0.6324555320336759*ne[2]+0.7071067811865475*ne[0]; 
  AEM(1,2) = 0.6324555320336759*ne[1]; 
  AEM(2,0) = 0.7071067811865475*ne[2]; 
  AEM(2,1) = 0.6324555320336759*ne[1]; 
  AEM(2,2) = 0.4517539514526256*ne[2]+0.7071067811865475*ne[0]; 
 
  bEV << ni[0],ni[1],ni[2]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(nRat,3,1) = xEV; 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    uei[0] = 0.7071067811865475*nRat[2]*ui[a0+2]-0.7071067811865475*nRat[2]*ue[a0+2]+0.7071067811865475*nRat[1]*ui[a0+1]-0.7071067811865475*nRat[1]*ue[a0+1]+0.7071067811865475*nRat[0]*ui[a0]-0.7071067811865475*nRat[0]*ue[a0]+ue[a0]; 
    uei[1] = 0.6324555320336759*nRat[1]*ui[a0+2]-0.6324555320336759*nRat[1]*ue[a0+2]+0.6324555320336759*nRat[2]*ui[a0+1]+0.7071067811865475*nRat[0]*ui[a0+1]-0.6324555320336759*nRat[2]*ue[a0+1]-0.7071067811865475*nRat[0]*ue[a0+1]+ue[a0+1]+0.7071067811865475*nRat[1]*ui[a0]-0.7071067811865475*nRat[1]*ue[a0]; 
    uei[2] = 0.4517539514526256*nRat[2]*ui[a0+2]+0.7071067811865475*nRat[0]*ui[a0+2]-0.4517539514526256*nRat[2]*ue[a0+2]-0.7071067811865475*nRat[0]*ue[a0+2]+ue[a0+2]+0.6324555320336759*nRat[1]*ui[a0+1]-0.6324555320336759*nRat[1]*ue[a0+1]+0.7071067811865475*nRat[2]*ui[a0]-0.7071067811865475*nRat[2]*ue[a0]; 
 
  } 
 
  // ..... Get the relative speed squared (ue-ui)^2 ..... // 
  double uRelSq[3]; 
  for (unsigned short int k=0; k<3; k++) 
  { 
    uRelSq[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    uRelSq[0] += 0.7071067811865475*ui[a0+2]*ui[a0+2]-1.414213562373095*ue[a0+2]*ui[a0+2]+0.7071067811865475*ue[a0+2]*ue[a0+2]+0.7071067811865475*ui[a0+1]*ui[a0+1]-1.414213562373095*ue[a0+1]*ui[a0+1]+0.7071067811865475*ue[a0+1]*ue[a0+1]+0.7071067811865475*ui[a0]*ui[a0]-1.414213562373095*ue[a0]*ui[a0]+0.7071067811865475*ue[a0]*ue[a0]; 
    uRelSq[1] += 1.264911064067352*ui[a0+1]*ui[a0+2]-1.264911064067352*ue[a0+1]*ui[a0+2]-1.264911064067352*ui[a0+1]*ue[a0+2]+1.264911064067352*ue[a0+1]*ue[a0+2]+1.414213562373095*ui[a0]*ui[a0+1]-1.414213562373095*ue[a0]*ui[a0+1]-1.414213562373095*ui[a0]*ue[a0+1]+1.414213562373095*ue[a0]*ue[a0+1]; 
    uRelSq[2] += 0.4517539514526256*ui[a0+2]*ui[a0+2]-0.9035079029052515*ue[a0+2]*ui[a0+2]+1.414213562373095*ui[a0]*ui[a0+2]-1.414213562373095*ue[a0]*ui[a0+2]+0.4517539514526256*ue[a0+2]*ue[a0+2]-1.414213562373095*ui[a0]*ue[a0+2]+1.414213562373095*ue[a0]*ue[a0+2]+0.6324555320336759*ui[a0+1]*ui[a0+1]-1.264911064067352*ue[a0+1]*ui[a0+1]+0.6324555320336759*ue[a0+1]*ue[a0+1]; 
  } 
 
  double rmRat = 1.0/mRat; 
  vtSqei[0] = ((-1.414213562373095*nRat[2]*vtSqe[2])-1.414213562373095*nRat[1]*vtSqe[1]-1.414213562373095*nRat[0]*vtSqe[0])*rmRat+1.414213562373095*nRat[2]*vtSqi[2]+0.3535533905932737*nRat[2]*uRelSq[2]+1.414213562373095*nRat[1]*vtSqi[1]+0.3535533905932737*nRat[1]*uRelSq[1]+1.414213562373095*nRat[0]*vtSqi[0]+vtSqe[0]+0.3535533905932737*nRat[0]*uRelSq[0]; 
  vtSqei[1] = ((-1.264911064067352*nRat[1]*vtSqe[2])-1.264911064067352*vtSqe[1]*nRat[2]-1.414213562373095*nRat[0]*vtSqe[1]-1.414213562373095*vtSqe[0]*nRat[1])*rmRat+1.264911064067352*nRat[1]*vtSqi[2]+0.3162277660168379*nRat[1]*uRelSq[2]+1.264911064067352*vtSqi[1]*nRat[2]+0.3162277660168379*uRelSq[1]*nRat[2]+1.414213562373095*nRat[0]*vtSqi[1]+vtSqe[1]+0.3535533905932737*nRat[0]*uRelSq[1]+1.414213562373095*vtSqi[0]*nRat[1]+0.3535533905932737*uRelSq[0]*nRat[1]; 
  vtSqei[2] = ((-0.9035079029052515*nRat[2]*vtSqe[2])-1.414213562373095*nRat[0]*vtSqe[2]-1.414213562373095*vtSqe[0]*nRat[2]-1.264911064067352*nRat[1]*vtSqe[1])*rmRat+0.9035079029052515*nRat[2]*vtSqi[2]+1.414213562373095*nRat[0]*vtSqi[2]+vtSqe[2]+0.2258769757263128*nRat[2]*uRelSq[2]+0.3535533905932737*nRat[0]*uRelSq[2]+1.414213562373095*vtSqi[0]*nRat[2]+0.3535533905932737*uRelSq[0]*nRat[2]+1.264911064067352*nRat[1]*vtSqi[1]+0.3162277660168379*nRat[1]*uRelSq[1]; 
 
} 
 
void CrossPrimMoments_VmeiLBO_1x2vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uei, double *vtSqei) 
{ 
  // mRat:           mass ratio = m_i/m_e. 
  // ne, ue, vtSqe:  electron number density, bulk flow velocity and T_e/m_e. 
  // ni, ui, vtSqi:  ion number density, bulk flow velocity and T_e/m_e. 
  // u_ei:           bulk flow velocity for electron-ion collision term in electron equation. 
  // vtSq_ei:        squared thermal speed, T_ei/m_e, for electron-ion collision term in electron equation. 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(4,4); 
  Eigen::VectorXd bEV(4); 
  Eigen::VectorXd xEV(4); 
 
  // ....... Compute density ratio through weak division ni/ne .......... // 
  double nRat[4]; 
 
  AEM(0,0) = 0.7071067811865475*ne[0]; 
  AEM(0,1) = 0.7071067811865475*ne[1]; 
  AEM(0,2) = 0.7071067811865475*ne[2]; 
  AEM(0,3) = 0.7071067811865475*ne[3]; 
  AEM(1,0) = 0.7071067811865475*ne[1]; 
  AEM(1,1) = 0.6324555320336759*ne[2]+0.7071067811865475*ne[0]; 
  AEM(1,2) = 0.6210590034081186*ne[3]+0.6324555320336759*ne[1]; 
  AEM(1,3) = 0.6210590034081186*ne[2]; 
  AEM(2,0) = 0.7071067811865475*ne[2]; 
  AEM(2,1) = 0.6210590034081186*ne[3]+0.6324555320336759*ne[1]; 
  AEM(2,2) = 0.4517539514526256*ne[2]+0.7071067811865475*ne[0]; 
  AEM(2,3) = 0.421637021355784*ne[3]+0.6210590034081186*ne[1]; 
  AEM(3,0) = 0.7071067811865475*ne[3]; 
  AEM(3,1) = 0.6210590034081186*ne[2]; 
  AEM(3,2) = 0.421637021355784*ne[3]+0.6210590034081186*ne[1]; 
  AEM(3,3) = 0.421637021355784*ne[2]+0.7071067811865475*ne[0]; 
 
  bEV << ni[0],ni[1],ni[2],ni[3]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(nRat,4,1) = xEV; 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uei[0] = 0.7071067811865475*nRat[3]*ui[a0+3]-0.7071067811865475*nRat[3]*ue[a0+3]+0.7071067811865475*nRat[2]*ui[a0+2]-0.7071067811865475*nRat[2]*ue[a0+2]+0.7071067811865475*nRat[1]*ui[a0+1]-0.7071067811865475*nRat[1]*ue[a0+1]+0.7071067811865475*nRat[0]*ui[a0]-0.7071067811865475*nRat[0]*ue[a0]+ue[a0]; 
    uei[1] = 0.6210590034081186*nRat[2]*ui[a0+3]-0.6210590034081186*nRat[2]*ue[a0+3]+0.6210590034081186*nRat[3]*ui[a0+2]+0.6324555320336759*nRat[1]*ui[a0+2]-0.6210590034081186*nRat[3]*ue[a0+2]-0.6324555320336759*nRat[1]*ue[a0+2]+0.6324555320336759*nRat[2]*ui[a0+1]+0.7071067811865475*nRat[0]*ui[a0+1]-0.6324555320336759*nRat[2]*ue[a0+1]-0.7071067811865475*nRat[0]*ue[a0+1]+ue[a0+1]+0.7071067811865475*nRat[1]*ui[a0]-0.7071067811865475*nRat[1]*ue[a0]; 
    uei[2] = 0.421637021355784*nRat[3]*ui[a0+3]+0.6210590034081186*nRat[1]*ui[a0+3]-0.421637021355784*nRat[3]*ue[a0+3]-0.6210590034081186*nRat[1]*ue[a0+3]+0.4517539514526256*nRat[2]*ui[a0+2]+0.7071067811865475*nRat[0]*ui[a0+2]-0.4517539514526256*nRat[2]*ue[a0+2]-0.7071067811865475*nRat[0]*ue[a0+2]+ue[a0+2]+0.6210590034081186*nRat[3]*ui[a0+1]+0.6324555320336759*nRat[1]*ui[a0+1]-0.6210590034081186*nRat[3]*ue[a0+1]-0.6324555320336759*nRat[1]*ue[a0+1]+0.7071067811865475*nRat[2]*ui[a0]-0.7071067811865475*nRat[2]*ue[a0]; 
    uei[3] = 0.421637021355784*nRat[2]*ui[a0+3]+0.7071067811865475*nRat[0]*ui[a0+3]-0.421637021355784*nRat[2]*ue[a0+3]-0.7071067811865475*nRat[0]*ue[a0+3]+ue[a0+3]+0.421637021355784*nRat[3]*ui[a0+2]+0.6210590034081186*nRat[1]*ui[a0+2]-0.421637021355784*nRat[3]*ue[a0+2]-0.6210590034081186*nRat[1]*ue[a0+2]+0.6210590034081186*nRat[2]*ui[a0+1]-0.6210590034081186*nRat[2]*ue[a0+1]+0.7071067811865475*nRat[3]*ui[a0]-0.7071067811865475*nRat[3]*ue[a0]; 
 
  } 
 
  // ..... Get the relative speed squared (ue-ui)^2 ..... // 
  double uRelSq[4]; 
  for (unsigned short int k=0; k<4; k++) 
  { 
    uRelSq[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uRelSq[0] += 0.7071067811865475*ui[a0+3]*ui[a0+3]-1.414213562373095*ue[a0+3]*ui[a0+3]+0.7071067811865475*ue[a0+3]*ue[a0+3]+0.7071067811865475*ui[a0+2]*ui[a0+2]-1.414213562373095*ue[a0+2]*ui[a0+2]+0.7071067811865475*ue[a0+2]*ue[a0+2]+0.7071067811865475*ui[a0+1]*ui[a0+1]-1.414213562373095*ue[a0+1]*ui[a0+1]+0.7071067811865475*ue[a0+1]*ue[a0+1]+0.7071067811865475*ui[a0]*ui[a0]-1.414213562373095*ue[a0]*ui[a0]+0.7071067811865475*ue[a0]*ue[a0]; 
    uRelSq[1] += 1.242118006816237*ui[a0+2]*ui[a0+3]-1.242118006816237*ue[a0+2]*ui[a0+3]-1.242118006816237*ui[a0+2]*ue[a0+3]+1.242118006816237*ue[a0+2]*ue[a0+3]+1.264911064067352*ui[a0+1]*ui[a0+2]-1.264911064067352*ue[a0+1]*ui[a0+2]-1.264911064067352*ui[a0+1]*ue[a0+2]+1.264911064067352*ue[a0+1]*ue[a0+2]+1.414213562373095*ui[a0]*ui[a0+1]-1.414213562373095*ue[a0]*ui[a0+1]-1.414213562373095*ui[a0]*ue[a0+1]+1.414213562373095*ue[a0]*ue[a0+1]; 
    uRelSq[2] += 0.421637021355784*ui[a0+3]*ui[a0+3]-0.8432740427115681*ue[a0+3]*ui[a0+3]+1.242118006816237*ui[a0+1]*ui[a0+3]-1.242118006816237*ue[a0+1]*ui[a0+3]+0.421637021355784*ue[a0+3]*ue[a0+3]-1.242118006816237*ui[a0+1]*ue[a0+3]+1.242118006816237*ue[a0+1]*ue[a0+3]+0.4517539514526256*ui[a0+2]*ui[a0+2]-0.9035079029052515*ue[a0+2]*ui[a0+2]+1.414213562373095*ui[a0]*ui[a0+2]-1.414213562373095*ue[a0]*ui[a0+2]+0.4517539514526256*ue[a0+2]*ue[a0+2]-1.414213562373095*ui[a0]*ue[a0+2]+1.414213562373095*ue[a0]*ue[a0+2]+0.6324555320336759*ui[a0+1]*ui[a0+1]-1.264911064067352*ue[a0+1]*ui[a0+1]+0.6324555320336759*ue[a0+1]*ue[a0+1]; 
    uRelSq[3] += 0.8432740427115681*ui[a0+2]*ui[a0+3]-0.8432740427115681*ue[a0+2]*ui[a0+3]+1.414213562373095*ui[a0]*ui[a0+3]-1.414213562373095*ue[a0]*ui[a0+3]-0.8432740427115681*ui[a0+2]*ue[a0+3]+0.8432740427115681*ue[a0+2]*ue[a0+3]-1.414213562373095*ui[a0]*ue[a0+3]+1.414213562373095*ue[a0]*ue[a0+3]+1.242118006816237*ui[a0+1]*ui[a0+2]-1.242118006816237*ue[a0+1]*ui[a0+2]-1.242118006816237*ui[a0+1]*ue[a0+2]+1.242118006816237*ue[a0+1]*ue[a0+2]; 
  } 
 
  double rmRat = 1.0/mRat; 
  vtSqei[0] = ((-1.414213562373095*nRat[3]*vtSqe[3])-1.414213562373095*nRat[2]*vtSqe[2]-1.414213562373095*nRat[1]*vtSqe[1]-1.414213562373095*nRat[0]*vtSqe[0])*rmRat+1.414213562373095*nRat[3]*vtSqi[3]+0.3535533905932737*nRat[3]*uRelSq[3]+1.414213562373095*nRat[2]*vtSqi[2]+0.3535533905932737*nRat[2]*uRelSq[2]+1.414213562373095*nRat[1]*vtSqi[1]+0.3535533905932737*nRat[1]*uRelSq[1]+1.414213562373095*nRat[0]*vtSqi[0]+vtSqe[0]+0.3535533905932737*nRat[0]*uRelSq[0]; 
  vtSqei[1] = ((-1.242118006816237*nRat[2]*vtSqe[3])-1.242118006816237*vtSqe[2]*nRat[3]-1.264911064067352*nRat[1]*vtSqe[2]-1.264911064067352*vtSqe[1]*nRat[2]-1.414213562373095*nRat[0]*vtSqe[1]-1.414213562373095*vtSqe[0]*nRat[1])*rmRat+1.242118006816237*nRat[2]*vtSqi[3]+0.3105295017040592*nRat[2]*uRelSq[3]+1.242118006816237*vtSqi[2]*nRat[3]+0.3105295017040592*uRelSq[2]*nRat[3]+1.264911064067352*nRat[1]*vtSqi[2]+0.3162277660168379*nRat[1]*uRelSq[2]+1.264911064067352*vtSqi[1]*nRat[2]+0.3162277660168379*uRelSq[1]*nRat[2]+1.414213562373095*nRat[0]*vtSqi[1]+vtSqe[1]+0.3535533905932737*nRat[0]*uRelSq[1]+1.414213562373095*vtSqi[0]*nRat[1]+0.3535533905932737*uRelSq[0]*nRat[1]; 
  vtSqei[2] = ((-0.8432740427115681*nRat[3]*vtSqe[3])-1.242118006816237*nRat[1]*vtSqe[3]-1.242118006816237*vtSqe[1]*nRat[3]-0.9035079029052515*nRat[2]*vtSqe[2]-1.414213562373095*nRat[0]*vtSqe[2]-1.414213562373095*vtSqe[0]*nRat[2]-1.264911064067352*nRat[1]*vtSqe[1])*rmRat+0.8432740427115681*nRat[3]*vtSqi[3]+1.242118006816237*nRat[1]*vtSqi[3]+0.210818510677892*nRat[3]*uRelSq[3]+0.3105295017040592*nRat[1]*uRelSq[3]+1.242118006816237*vtSqi[1]*nRat[3]+0.3105295017040592*uRelSq[1]*nRat[3]+0.9035079029052515*nRat[2]*vtSqi[2]+1.414213562373095*nRat[0]*vtSqi[2]+vtSqe[2]+0.2258769757263128*nRat[2]*uRelSq[2]+0.3535533905932737*nRat[0]*uRelSq[2]+1.414213562373095*vtSqi[0]*nRat[2]+0.3535533905932737*uRelSq[0]*nRat[2]+1.264911064067352*nRat[1]*vtSqi[1]+0.3162277660168379*nRat[1]*uRelSq[1]; 
  vtSqei[3] = ((-0.8432740427115681*nRat[2]*vtSqe[3])-1.414213562373095*nRat[0]*vtSqe[3]-0.8432740427115681*vtSqe[2]*nRat[3]-1.414213562373095*vtSqe[0]*nRat[3]-1.242118006816237*nRat[1]*vtSqe[2]-1.242118006816237*vtSqe[1]*nRat[2])*rmRat+0.8432740427115681*nRat[2]*vtSqi[3]+1.414213562373095*nRat[0]*vtSqi[3]+vtSqe[3]+0.210818510677892*nRat[2]*uRelSq[3]+0.3535533905932737*nRat[0]*uRelSq[3]+0.8432740427115681*vtSqi[2]*nRat[3]+0.210818510677892*uRelSq[2]*nRat[3]+1.414213562373095*vtSqi[0]*nRat[3]+0.3535533905932737*uRelSq[0]*nRat[3]+1.242118006816237*nRat[1]*vtSqi[2]+0.3105295017040592*nRat[1]*uRelSq[2]+1.242118006816237*vtSqi[1]*nRat[2]+0.3105295017040592*uRelSq[1]*nRat[2]; 
 
} 
 
void CrossPrimMoments_VmeiLBO_1x2vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uei, double *vtSqei) 
{ 
  // mRat:           mass ratio = m_i/m_e. 
  // ne, ue, vtSqe:  electron number density, bulk flow velocity and T_e/m_e. 
  // ni, ui, vtSqi:  ion number density, bulk flow velocity and T_e/m_e. 
  // u_ei:           bulk flow velocity for electron-ion collision term in electron equation. 
  // vtSq_ei:        squared thermal speed, T_ei/m_e, for electron-ion collision term in electron equation. 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(5,5); 
  Eigen::VectorXd bEV(5); 
  Eigen::VectorXd xEV(5); 
 
  // ....... Compute density ratio through weak division ni/ne .......... // 
  double nRat[5]; 
 
  AEM(0,0) = 0.7071067811865475*ne[0]; 
  AEM(0,1) = 0.7071067811865475*ne[1]; 
  AEM(0,2) = 0.7071067811865475*ne[2]; 
  AEM(0,3) = 0.7071067811865475*ne[3]; 
  AEM(0,4) = 0.7071067811865475*ne[4]; 
  AEM(1,0) = 0.7071067811865475*ne[1]; 
  AEM(1,1) = 0.6324555320336759*ne[2]+0.7071067811865475*ne[0]; 
  AEM(1,2) = 0.6210590034081186*ne[3]+0.6324555320336759*ne[1]; 
  AEM(1,3) = 0.6172133998483679*ne[4]+0.6210590034081186*ne[2]; 
  AEM(1,4) = 0.6172133998483679*ne[3]; 
  AEM(2,0) = 0.7071067811865475*ne[2]; 
  AEM(2,1) = 0.6210590034081186*ne[3]+0.6324555320336759*ne[1]; 
  AEM(2,2) = 0.6060915267313265*ne[4]+0.4517539514526256*ne[2]+0.7071067811865475*ne[0]; 
  AEM(2,3) = 0.421637021355784*ne[3]+0.6210590034081186*ne[1]; 
  AEM(2,4) = 0.410685410411478*ne[4]+0.6060915267313265*ne[2]; 
  AEM(3,0) = 0.7071067811865475*ne[3]; 
  AEM(3,1) = 0.6172133998483679*ne[4]+0.6210590034081186*ne[2]; 
  AEM(3,2) = 0.421637021355784*ne[3]+0.6210590034081186*ne[1]; 
  AEM(3,3) = 0.385694607919935*ne[4]+0.421637021355784*ne[2]+0.7071067811865475*ne[0]; 
  AEM(3,4) = 0.385694607919935*ne[3]+0.6172133998483679*ne[1]; 
  AEM(4,0) = 0.7071067811865475*ne[4]; 
  AEM(4,1) = 0.6172133998483679*ne[3]; 
  AEM(4,2) = 0.410685410411478*ne[4]+0.6060915267313265*ne[2]; 
  AEM(4,3) = 0.385694607919935*ne[3]+0.6172133998483679*ne[1]; 
  AEM(4,4) = 0.3433105850715905*ne[4]+0.410685410411478*ne[2]+0.7071067811865475*ne[0]; 
 
  bEV << ni[0],ni[1],ni[2],ni[3],ni[4]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(nRat,5,1) = xEV; 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 5*vd; 
    uei[0] = 0.7071067811865475*nRat[4]*ui[a0+4]-0.7071067811865475*nRat[4]*ue[a0+4]+0.7071067811865475*nRat[3]*ui[a0+3]-0.7071067811865475*nRat[3]*ue[a0+3]+0.7071067811865475*nRat[2]*ui[a0+2]-0.7071067811865475*nRat[2]*ue[a0+2]+0.7071067811865475*nRat[1]*ui[a0+1]-0.7071067811865475*nRat[1]*ue[a0+1]+0.7071067811865475*nRat[0]*ui[a0]-0.7071067811865475*nRat[0]*ue[a0]+ue[a0]; 
    uei[1] = 0.6172133998483679*nRat[3]*ui[a0+4]-0.6172133998483679*nRat[3]*ue[a0+4]+0.6172133998483679*nRat[4]*ui[a0+3]+0.6210590034081186*nRat[2]*ui[a0+3]-0.6172133998483679*nRat[4]*ue[a0+3]-0.6210590034081186*nRat[2]*ue[a0+3]+0.6210590034081186*nRat[3]*ui[a0+2]+0.6324555320336759*nRat[1]*ui[a0+2]-0.6210590034081186*nRat[3]*ue[a0+2]-0.6324555320336759*nRat[1]*ue[a0+2]+0.6324555320336759*nRat[2]*ui[a0+1]+0.7071067811865475*nRat[0]*ui[a0+1]-0.6324555320336759*nRat[2]*ue[a0+1]-0.7071067811865475*nRat[0]*ue[a0+1]+ue[a0+1]+0.7071067811865475*nRat[1]*ui[a0]-0.7071067811865475*nRat[1]*ue[a0]; 
    uei[2] = 0.410685410411478*nRat[4]*ui[a0+4]+0.6060915267313265*nRat[2]*ui[a0+4]-0.410685410411478*nRat[4]*ue[a0+4]-0.6060915267313265*nRat[2]*ue[a0+4]+0.421637021355784*nRat[3]*ui[a0+3]+0.6210590034081186*nRat[1]*ui[a0+3]-0.421637021355784*nRat[3]*ue[a0+3]-0.6210590034081186*nRat[1]*ue[a0+3]+0.6060915267313265*nRat[4]*ui[a0+2]+0.4517539514526256*nRat[2]*ui[a0+2]+0.7071067811865475*nRat[0]*ui[a0+2]-0.6060915267313265*nRat[4]*ue[a0+2]-0.4517539514526256*nRat[2]*ue[a0+2]-0.7071067811865475*nRat[0]*ue[a0+2]+ue[a0+2]+0.6210590034081186*nRat[3]*ui[a0+1]+0.6324555320336759*nRat[1]*ui[a0+1]-0.6210590034081186*nRat[3]*ue[a0+1]-0.6324555320336759*nRat[1]*ue[a0+1]+0.7071067811865475*nRat[2]*ui[a0]-0.7071067811865475*nRat[2]*ue[a0]; 
    uei[3] = 0.385694607919935*nRat[3]*ui[a0+4]+0.6172133998483679*nRat[1]*ui[a0+4]-0.385694607919935*nRat[3]*ue[a0+4]-0.6172133998483679*nRat[1]*ue[a0+4]+0.385694607919935*nRat[4]*ui[a0+3]+0.421637021355784*nRat[2]*ui[a0+3]+0.7071067811865475*nRat[0]*ui[a0+3]-0.385694607919935*nRat[4]*ue[a0+3]-0.421637021355784*nRat[2]*ue[a0+3]-0.7071067811865475*nRat[0]*ue[a0+3]+ue[a0+3]+0.421637021355784*nRat[3]*ui[a0+2]+0.6210590034081186*nRat[1]*ui[a0+2]-0.421637021355784*nRat[3]*ue[a0+2]-0.6210590034081186*nRat[1]*ue[a0+2]+0.6172133998483679*nRat[4]*ui[a0+1]+0.6210590034081186*nRat[2]*ui[a0+1]-0.6172133998483679*nRat[4]*ue[a0+1]-0.6210590034081186*nRat[2]*ue[a0+1]+0.7071067811865475*nRat[3]*ui[a0]-0.7071067811865475*nRat[3]*ue[a0]; 
    uei[4] = 0.3433105850715905*nRat[4]*ui[a0+4]+0.410685410411478*nRat[2]*ui[a0+4]+0.7071067811865475*nRat[0]*ui[a0+4]-0.3433105850715905*nRat[4]*ue[a0+4]-0.410685410411478*nRat[2]*ue[a0+4]-0.7071067811865475*nRat[0]*ue[a0+4]+ue[a0+4]+0.385694607919935*nRat[3]*ui[a0+3]+0.6172133998483679*nRat[1]*ui[a0+3]-0.385694607919935*nRat[3]*ue[a0+3]-0.6172133998483679*nRat[1]*ue[a0+3]+0.410685410411478*nRat[4]*ui[a0+2]+0.6060915267313265*nRat[2]*ui[a0+2]-0.410685410411478*nRat[4]*ue[a0+2]-0.6060915267313265*nRat[2]*ue[a0+2]+0.6172133998483679*nRat[3]*ui[a0+1]-0.6172133998483679*nRat[3]*ue[a0+1]+0.7071067811865475*nRat[4]*ui[a0]-0.7071067811865475*nRat[4]*ue[a0]; 
 
  } 
 
  // ..... Get the relative speed squared (ue-ui)^2 ..... // 
  double uRelSq[5]; 
  for (unsigned short int k=0; k<5; k++) 
  { 
    uRelSq[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 5*vd; 
    uRelSq[0] += 0.7071067811865475*ui[a0+4]*ui[a0+4]-1.414213562373095*ue[a0+4]*ui[a0+4]+0.7071067811865475*ue[a0+4]*ue[a0+4]+0.7071067811865475*ui[a0+3]*ui[a0+3]-1.414213562373095*ue[a0+3]*ui[a0+3]+0.7071067811865475*ue[a0+3]*ue[a0+3]+0.7071067811865475*ui[a0+2]*ui[a0+2]-1.414213562373095*ue[a0+2]*ui[a0+2]+0.7071067811865475*ue[a0+2]*ue[a0+2]+0.7071067811865475*ui[a0+1]*ui[a0+1]-1.414213562373095*ue[a0+1]*ui[a0+1]+0.7071067811865475*ue[a0+1]*ue[a0+1]+0.7071067811865475*ui[a0]*ui[a0]-1.414213562373095*ue[a0]*ui[a0]+0.7071067811865475*ue[a0]*ue[a0]; 
    uRelSq[1] += 1.234426799696736*ui[a0+3]*ui[a0+4]-1.234426799696736*ue[a0+3]*ui[a0+4]-1.234426799696736*ui[a0+3]*ue[a0+4]+1.234426799696736*ue[a0+3]*ue[a0+4]+1.242118006816237*ui[a0+2]*ui[a0+3]-1.242118006816237*ue[a0+2]*ui[a0+3]-1.242118006816237*ui[a0+2]*ue[a0+3]+1.242118006816237*ue[a0+2]*ue[a0+3]+1.264911064067352*ui[a0+1]*ui[a0+2]-1.264911064067352*ue[a0+1]*ui[a0+2]-1.264911064067352*ui[a0+1]*ue[a0+2]+1.264911064067352*ue[a0+1]*ue[a0+2]+1.414213562373095*ui[a0]*ui[a0+1]-1.414213562373095*ue[a0]*ui[a0+1]-1.414213562373095*ui[a0]*ue[a0+1]+1.414213562373095*ue[a0]*ue[a0+1]; 
    uRelSq[2] += 0.410685410411478*ui[a0+4]*ui[a0+4]-0.8213708208229562*ue[a0+4]*ui[a0+4]+1.212183053462653*ui[a0+2]*ui[a0+4]-1.212183053462653*ue[a0+2]*ui[a0+4]+0.410685410411478*ue[a0+4]*ue[a0+4]-1.212183053462653*ui[a0+2]*ue[a0+4]+1.212183053462653*ue[a0+2]*ue[a0+4]+0.421637021355784*ui[a0+3]*ui[a0+3]-0.8432740427115681*ue[a0+3]*ui[a0+3]+1.242118006816237*ui[a0+1]*ui[a0+3]-1.242118006816237*ue[a0+1]*ui[a0+3]+0.421637021355784*ue[a0+3]*ue[a0+3]-1.242118006816237*ui[a0+1]*ue[a0+3]+1.242118006816237*ue[a0+1]*ue[a0+3]+0.4517539514526256*ui[a0+2]*ui[a0+2]-0.9035079029052515*ue[a0+2]*ui[a0+2]+1.414213562373095*ui[a0]*ui[a0+2]-1.414213562373095*ue[a0]*ui[a0+2]+0.4517539514526256*ue[a0+2]*ue[a0+2]-1.414213562373095*ui[a0]*ue[a0+2]+1.414213562373095*ue[a0]*ue[a0+2]+0.6324555320336759*ui[a0+1]*ui[a0+1]-1.264911064067352*ue[a0+1]*ui[a0+1]+0.6324555320336759*ue[a0+1]*ue[a0+1]; 
    uRelSq[3] += 0.7713892158398702*ui[a0+3]*ui[a0+4]-0.7713892158398702*ue[a0+3]*ui[a0+4]+1.234426799696736*ui[a0+1]*ui[a0+4]-1.234426799696736*ue[a0+1]*ui[a0+4]-0.7713892158398702*ui[a0+3]*ue[a0+4]+0.7713892158398702*ue[a0+3]*ue[a0+4]-1.234426799696736*ui[a0+1]*ue[a0+4]+1.234426799696736*ue[a0+1]*ue[a0+4]+0.8432740427115681*ui[a0+2]*ui[a0+3]-0.8432740427115681*ue[a0+2]*ui[a0+3]+1.414213562373095*ui[a0]*ui[a0+3]-1.414213562373095*ue[a0]*ui[a0+3]-0.8432740427115681*ui[a0+2]*ue[a0+3]+0.8432740427115681*ue[a0+2]*ue[a0+3]-1.414213562373095*ui[a0]*ue[a0+3]+1.414213562373095*ue[a0]*ue[a0+3]+1.242118006816237*ui[a0+1]*ui[a0+2]-1.242118006816237*ue[a0+1]*ui[a0+2]-1.242118006816237*ui[a0+1]*ue[a0+2]+1.242118006816237*ue[a0+1]*ue[a0+2]; 
    uRelSq[4] += 0.3433105850715905*ui[a0+4]*ui[a0+4]-0.6866211701431811*ue[a0+4]*ui[a0+4]+0.8213708208229562*ui[a0+2]*ui[a0+4]-0.8213708208229562*ue[a0+2]*ui[a0+4]+1.414213562373095*ui[a0]*ui[a0+4]-1.414213562373095*ue[a0]*ui[a0+4]+0.3433105850715905*ue[a0+4]*ue[a0+4]-0.8213708208229562*ui[a0+2]*ue[a0+4]+0.8213708208229562*ue[a0+2]*ue[a0+4]-1.414213562373095*ui[a0]*ue[a0+4]+1.414213562373095*ue[a0]*ue[a0+4]+0.385694607919935*ui[a0+3]*ui[a0+3]-0.7713892158398702*ue[a0+3]*ui[a0+3]+1.234426799696736*ui[a0+1]*ui[a0+3]-1.234426799696736*ue[a0+1]*ui[a0+3]+0.385694607919935*ue[a0+3]*ue[a0+3]-1.234426799696736*ui[a0+1]*ue[a0+3]+1.234426799696736*ue[a0+1]*ue[a0+3]+0.6060915267313265*ui[a0+2]*ui[a0+2]-1.212183053462653*ue[a0+2]*ui[a0+2]+0.6060915267313265*ue[a0+2]*ue[a0+2]; 
  } 
 
  double rmRat = 1.0/mRat; 
  vtSqei[0] = ((-1.414213562373095*nRat[4]*vtSqe[4])-1.414213562373095*nRat[3]*vtSqe[3]-1.414213562373095*nRat[2]*vtSqe[2]-1.414213562373095*nRat[1]*vtSqe[1]-1.414213562373095*nRat[0]*vtSqe[0])*rmRat+1.414213562373095*nRat[4]*vtSqi[4]+0.3535533905932737*nRat[4]*uRelSq[4]+1.414213562373095*nRat[3]*vtSqi[3]+0.3535533905932737*nRat[3]*uRelSq[3]+1.414213562373095*nRat[2]*vtSqi[2]+0.3535533905932737*nRat[2]*uRelSq[2]+1.414213562373095*nRat[1]*vtSqi[1]+0.3535533905932737*nRat[1]*uRelSq[1]+1.414213562373095*nRat[0]*vtSqi[0]+vtSqe[0]+0.3535533905932737*nRat[0]*uRelSq[0]; 
  vtSqei[1] = ((-1.234426799696736*nRat[3]*vtSqe[4])-1.234426799696736*vtSqe[3]*nRat[4]-1.242118006816237*nRat[2]*vtSqe[3]-1.242118006816237*vtSqe[2]*nRat[3]-1.264911064067352*nRat[1]*vtSqe[2]-1.264911064067352*vtSqe[1]*nRat[2]-1.414213562373095*nRat[0]*vtSqe[1]-1.414213562373095*vtSqe[0]*nRat[1])*rmRat+1.234426799696736*nRat[3]*vtSqi[4]+0.3086066999241838*nRat[3]*uRelSq[4]+1.234426799696736*vtSqi[3]*nRat[4]+0.3086066999241838*uRelSq[3]*nRat[4]+1.242118006816237*nRat[2]*vtSqi[3]+0.3105295017040592*nRat[2]*uRelSq[3]+1.242118006816237*vtSqi[2]*nRat[3]+0.3105295017040592*uRelSq[2]*nRat[3]+1.264911064067352*nRat[1]*vtSqi[2]+0.3162277660168379*nRat[1]*uRelSq[2]+1.264911064067352*vtSqi[1]*nRat[2]+0.3162277660168379*uRelSq[1]*nRat[2]+1.414213562373095*nRat[0]*vtSqi[1]+vtSqe[1]+0.3535533905932737*nRat[0]*uRelSq[1]+1.414213562373095*vtSqi[0]*nRat[1]+0.3535533905932737*uRelSq[0]*nRat[1]; 
  vtSqei[2] = ((-0.8213708208229562*nRat[4]*vtSqe[4])-1.212183053462653*nRat[2]*vtSqe[4]-1.212183053462653*vtSqe[2]*nRat[4]-0.8432740427115681*nRat[3]*vtSqe[3]-1.242118006816237*nRat[1]*vtSqe[3]-1.242118006816237*vtSqe[1]*nRat[3]-0.9035079029052515*nRat[2]*vtSqe[2]-1.414213562373095*nRat[0]*vtSqe[2]-1.414213562373095*vtSqe[0]*nRat[2]-1.264911064067352*nRat[1]*vtSqe[1])*rmRat+0.8213708208229562*nRat[4]*vtSqi[4]+1.212183053462653*nRat[2]*vtSqi[4]+0.205342705205739*nRat[4]*uRelSq[4]+0.3030457633656632*nRat[2]*uRelSq[4]+1.212183053462653*vtSqi[2]*nRat[4]+0.3030457633656632*uRelSq[2]*nRat[4]+0.8432740427115681*nRat[3]*vtSqi[3]+1.242118006816237*nRat[1]*vtSqi[3]+0.210818510677892*nRat[3]*uRelSq[3]+0.3105295017040592*nRat[1]*uRelSq[3]+1.242118006816237*vtSqi[1]*nRat[3]+0.3105295017040592*uRelSq[1]*nRat[3]+0.9035079029052515*nRat[2]*vtSqi[2]+1.414213562373095*nRat[0]*vtSqi[2]+vtSqe[2]+0.2258769757263128*nRat[2]*uRelSq[2]+0.3535533905932737*nRat[0]*uRelSq[2]+1.414213562373095*vtSqi[0]*nRat[2]+0.3535533905932737*uRelSq[0]*nRat[2]+1.264911064067352*nRat[1]*vtSqi[1]+0.3162277660168379*nRat[1]*uRelSq[1]; 
  vtSqei[3] = ((-0.7713892158398702*nRat[3]*vtSqe[4])-1.234426799696736*nRat[1]*vtSqe[4]-0.7713892158398702*vtSqe[3]*nRat[4]-1.234426799696736*vtSqe[1]*nRat[4]-0.8432740427115681*nRat[2]*vtSqe[3]-1.414213562373095*nRat[0]*vtSqe[3]-0.8432740427115681*vtSqe[2]*nRat[3]-1.414213562373095*vtSqe[0]*nRat[3]-1.242118006816237*nRat[1]*vtSqe[2]-1.242118006816237*vtSqe[1]*nRat[2])*rmRat+0.7713892158398702*nRat[3]*vtSqi[4]+1.234426799696736*nRat[1]*vtSqi[4]+0.1928473039599675*nRat[3]*uRelSq[4]+0.3086066999241838*nRat[1]*uRelSq[4]+0.7713892158398702*vtSqi[3]*nRat[4]+0.1928473039599675*uRelSq[3]*nRat[4]+1.234426799696736*vtSqi[1]*nRat[4]+0.3086066999241838*uRelSq[1]*nRat[4]+0.8432740427115681*nRat[2]*vtSqi[3]+1.414213562373095*nRat[0]*vtSqi[3]+vtSqe[3]+0.210818510677892*nRat[2]*uRelSq[3]+0.3535533905932737*nRat[0]*uRelSq[3]+0.8432740427115681*vtSqi[2]*nRat[3]+0.210818510677892*uRelSq[2]*nRat[3]+1.414213562373095*vtSqi[0]*nRat[3]+0.3535533905932737*uRelSq[0]*nRat[3]+1.242118006816237*nRat[1]*vtSqi[2]+0.3105295017040592*nRat[1]*uRelSq[2]+1.242118006816237*vtSqi[1]*nRat[2]+0.3105295017040592*uRelSq[1]*nRat[2]; 
  vtSqei[4] = ((-0.6866211701431811*nRat[4]*vtSqe[4])-0.8213708208229562*nRat[2]*vtSqe[4]-1.414213562373095*nRat[0]*vtSqe[4]-0.8213708208229562*vtSqe[2]*nRat[4]-1.414213562373095*vtSqe[0]*nRat[4]-0.7713892158398702*nRat[3]*vtSqe[3]-1.234426799696736*nRat[1]*vtSqe[3]-1.234426799696736*vtSqe[1]*nRat[3]-1.212183053462653*nRat[2]*vtSqe[2])*rmRat+0.6866211701431811*nRat[4]*vtSqi[4]+0.8213708208229562*nRat[2]*vtSqi[4]+1.414213562373095*nRat[0]*vtSqi[4]+vtSqe[4]+0.1716552925357952*nRat[4]*uRelSq[4]+0.205342705205739*nRat[2]*uRelSq[4]+0.3535533905932737*nRat[0]*uRelSq[4]+0.8213708208229562*vtSqi[2]*nRat[4]+0.205342705205739*uRelSq[2]*nRat[4]+1.414213562373095*vtSqi[0]*nRat[4]+0.3535533905932737*uRelSq[0]*nRat[4]+0.7713892158398702*nRat[3]*vtSqi[3]+1.234426799696736*nRat[1]*vtSqi[3]+0.1928473039599675*nRat[3]*uRelSq[3]+0.3086066999241838*nRat[1]*uRelSq[3]+1.234426799696736*vtSqi[1]*nRat[3]+0.3086066999241838*uRelSq[1]*nRat[3]+1.212183053462653*nRat[2]*vtSqi[2]+0.3030457633656632*nRat[2]*uRelSq[2]; 
 
} 
 
void CrossPrimMoments_VmieLBO_1x2vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uei, double *vtSqie) 
{ 
  // mRat:           mass ratio = m_i/m_e. 
  // ne, ue, vtSqe:  electron number density, bulk flow velocity and T_e/m_e. 
  // ni, ui, vtSqi:  ion number density, bulk flow velocity and T_e/m_e. 
  // u_ie:           bulk flow velocity for ion-electron collision term in ion equation. 
  // vtSq_ie:        squared thermal speed, T_ie/m_i, for ion-electron collision term in ion equation. 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(2,2); 
  Eigen::VectorXd bEV(2); 
  Eigen::VectorXd xEV(2); 
 
  // ....... Compute density ratio through weak division ni/ne .......... // 
  double nRat[2]; 
 
  AEM(0,0) = 0.7071067811865475*ne[0]; 
  AEM(0,1) = 0.7071067811865475*ne[1]; 
  AEM(1,0) = 0.7071067811865475*ne[1]; 
  AEM(1,1) = 0.7071067811865475*ne[0]; 
 
  bEV << ni[0],ni[1]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(nRat,2,1) = xEV; 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    uei[0] = (-0.7071067811865475*nRat[1]*ui[a0+1])+0.7071067811865475*nRat[1]*ue[a0+1]-0.7071067811865475*nRat[0]*ui[a0]+ui[a0]+0.7071067811865475*nRat[0]*ue[a0]; 
    uei[1] = (-0.7071067811865475*nRat[0]*ui[a0+1])+ui[a0+1]+0.7071067811865475*nRat[0]*ue[a0+1]-0.7071067811865475*nRat[1]*ui[a0]+0.7071067811865475*nRat[1]*ue[a0]; 
 
  } 
 
  // ..... Get the relative speed squared (ue-ui)^2 ..... // 
  double uRelSq[2]; 
  for (unsigned short int k=0; k<2; k++) 
  { 
    uRelSq[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    uRelSq[0] += 0.7071067811865475*ui[a0+1]*ui[a0+1]-1.414213562373095*ue[a0+1]*ui[a0+1]+0.7071067811865475*ue[a0+1]*ue[a0+1]+0.7071067811865475*ui[a0]*ui[a0]-1.414213562373095*ue[a0]*ui[a0]+0.7071067811865475*ue[a0]*ue[a0]; 
    uRelSq[1] += 1.414213562373095*ui[a0]*ui[a0+1]-1.414213562373095*ue[a0]*ui[a0+1]-1.414213562373095*ui[a0]*ue[a0+1]+1.414213562373095*ue[a0]*ue[a0+1]; 
  } 
 
  double rmRat = 1.0/mRat; 
  vtSqie[0] = (1.414213562373095*nRat[1]*vtSqe[1]+1.414213562373095*nRat[0]*vtSqe[0])*rmRat-1.414213562373095*nRat[1]*vtSqi[1]-0.1178511301977579*nRat[1]*uRelSq[1]-1.414213562373095*nRat[0]*vtSqi[0]+vtSqi[0]-0.1178511301977579*nRat[0]*uRelSq[0]; 
  vtSqie[1] = (1.414213562373095*nRat[0]*vtSqe[1]+1.414213562373095*vtSqe[0]*nRat[1])*rmRat-1.414213562373095*nRat[0]*vtSqi[1]+vtSqi[1]-0.1178511301977579*nRat[0]*uRelSq[1]-1.414213562373095*vtSqi[0]*nRat[1]-0.1178511301977579*uRelSq[0]*nRat[1]; 
 
} 
 
void CrossPrimMoments_VmieLBO_1x2vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uei, double *vtSqie) 
{ 
  // mRat:           mass ratio = m_i/m_e. 
  // ne, ue, vtSqe:  electron number density, bulk flow velocity and T_e/m_e. 
  // ni, ui, vtSqi:  ion number density, bulk flow velocity and T_e/m_e. 
  // u_ie:           bulk flow velocity for ion-electron collision term in ion equation. 
  // vtSq_ie:        squared thermal speed, T_ie/m_i, for ion-electron collision term in ion equation. 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(3,3); 
  Eigen::VectorXd bEV(3); 
  Eigen::VectorXd xEV(3); 
 
  // ....... Compute density ratio through weak division ni/ne .......... // 
  double nRat[3]; 
 
  AEM(0,0) = 0.7071067811865475*ne[0]; 
  AEM(0,1) = 0.7071067811865475*ne[1]; 
  AEM(0,2) = 0.7071067811865475*ne[2]; 
  AEM(1,0) = 0.7071067811865475*ne[1]; 
  AEM(1,1) = 0.6324555320336759*ne[2]+0.7071067811865475*ne[0]; 
  AEM(1,2) = 0.6324555320336759*ne[1]; 
  AEM(2,0) = 0.7071067811865475*ne[2]; 
  AEM(2,1) = 0.6324555320336759*ne[1]; 
  AEM(2,2) = 0.4517539514526256*ne[2]+0.7071067811865475*ne[0]; 
 
  bEV << ni[0],ni[1],ni[2]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(nRat,3,1) = xEV; 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    uei[0] = (-0.7071067811865475*nRat[2]*ui[a0+2])+0.7071067811865475*nRat[2]*ue[a0+2]-0.7071067811865475*nRat[1]*ui[a0+1]+0.7071067811865475*nRat[1]*ue[a0+1]-0.7071067811865475*nRat[0]*ui[a0]+ui[a0]+0.7071067811865475*nRat[0]*ue[a0]; 
    uei[1] = (-0.6324555320336759*nRat[1]*ui[a0+2])+0.6324555320336759*nRat[1]*ue[a0+2]-0.6324555320336759*nRat[2]*ui[a0+1]-0.7071067811865475*nRat[0]*ui[a0+1]+ui[a0+1]+0.6324555320336759*nRat[2]*ue[a0+1]+0.7071067811865475*nRat[0]*ue[a0+1]-0.7071067811865475*nRat[1]*ui[a0]+0.7071067811865475*nRat[1]*ue[a0]; 
    uei[2] = (-0.4517539514526256*nRat[2]*ui[a0+2])-0.7071067811865475*nRat[0]*ui[a0+2]+ui[a0+2]+0.4517539514526256*nRat[2]*ue[a0+2]+0.7071067811865475*nRat[0]*ue[a0+2]-0.6324555320336759*nRat[1]*ui[a0+1]+0.6324555320336759*nRat[1]*ue[a0+1]-0.7071067811865475*nRat[2]*ui[a0]+0.7071067811865475*nRat[2]*ue[a0]; 
 
  } 
 
  // ..... Get the relative speed squared (ue-ui)^2 ..... // 
  double uRelSq[3]; 
  for (unsigned short int k=0; k<3; k++) 
  { 
    uRelSq[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    uRelSq[0] += 0.7071067811865475*ui[a0+2]*ui[a0+2]-1.414213562373095*ue[a0+2]*ui[a0+2]+0.7071067811865475*ue[a0+2]*ue[a0+2]+0.7071067811865475*ui[a0+1]*ui[a0+1]-1.414213562373095*ue[a0+1]*ui[a0+1]+0.7071067811865475*ue[a0+1]*ue[a0+1]+0.7071067811865475*ui[a0]*ui[a0]-1.414213562373095*ue[a0]*ui[a0]+0.7071067811865475*ue[a0]*ue[a0]; 
    uRelSq[1] += 1.264911064067352*ui[a0+1]*ui[a0+2]-1.264911064067352*ue[a0+1]*ui[a0+2]-1.264911064067352*ui[a0+1]*ue[a0+2]+1.264911064067352*ue[a0+1]*ue[a0+2]+1.414213562373095*ui[a0]*ui[a0+1]-1.414213562373095*ue[a0]*ui[a0+1]-1.414213562373095*ui[a0]*ue[a0+1]+1.414213562373095*ue[a0]*ue[a0+1]; 
    uRelSq[2] += 0.4517539514526256*ui[a0+2]*ui[a0+2]-0.9035079029052515*ue[a0+2]*ui[a0+2]+1.414213562373095*ui[a0]*ui[a0+2]-1.414213562373095*ue[a0]*ui[a0+2]+0.4517539514526256*ue[a0+2]*ue[a0+2]-1.414213562373095*ui[a0]*ue[a0+2]+1.414213562373095*ue[a0]*ue[a0+2]+0.6324555320336759*ui[a0+1]*ui[a0+1]-1.264911064067352*ue[a0+1]*ui[a0+1]+0.6324555320336759*ue[a0+1]*ue[a0+1]; 
  } 
 
  double rmRat = 1.0/mRat; 
  vtSqie[0] = (1.414213562373095*nRat[2]*vtSqe[2]+1.414213562373095*nRat[1]*vtSqe[1]+1.414213562373095*nRat[0]*vtSqe[0])*rmRat-1.414213562373095*nRat[2]*vtSqi[2]-0.1178511301977579*nRat[2]*uRelSq[2]-1.414213562373095*nRat[1]*vtSqi[1]-0.1178511301977579*nRat[1]*uRelSq[1]-1.414213562373095*nRat[0]*vtSqi[0]+vtSqi[0]-0.1178511301977579*nRat[0]*uRelSq[0]; 
  vtSqie[1] = (1.264911064067352*nRat[1]*vtSqe[2]+1.264911064067352*vtSqe[1]*nRat[2]+1.414213562373095*nRat[0]*vtSqe[1]+1.414213562373095*vtSqe[0]*nRat[1])*rmRat-1.264911064067352*nRat[1]*vtSqi[2]-0.105409255338946*nRat[1]*uRelSq[2]-1.264911064067352*vtSqi[1]*nRat[2]-0.105409255338946*uRelSq[1]*nRat[2]-1.414213562373095*nRat[0]*vtSqi[1]+vtSqi[1]-0.1178511301977579*nRat[0]*uRelSq[1]-1.414213562373095*vtSqi[0]*nRat[1]-0.1178511301977579*uRelSq[0]*nRat[1]; 
  vtSqie[2] = (0.9035079029052515*nRat[2]*vtSqe[2]+1.414213562373095*nRat[0]*vtSqe[2]+1.414213562373095*vtSqe[0]*nRat[2]+1.264911064067352*nRat[1]*vtSqe[1])*rmRat-0.9035079029052515*nRat[2]*vtSqi[2]-1.414213562373095*nRat[0]*vtSqi[2]+vtSqi[2]-0.07529232524210427*nRat[2]*uRelSq[2]-0.1178511301977579*nRat[0]*uRelSq[2]-1.414213562373095*vtSqi[0]*nRat[2]-0.1178511301977579*uRelSq[0]*nRat[2]-1.264911064067352*nRat[1]*vtSqi[1]-0.105409255338946*nRat[1]*uRelSq[1]; 
 
} 
 
void CrossPrimMoments_VmieLBO_1x2vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uei, double *vtSqie) 
{ 
  // mRat:           mass ratio = m_i/m_e. 
  // ne, ue, vtSqe:  electron number density, bulk flow velocity and T_e/m_e. 
  // ni, ui, vtSqi:  ion number density, bulk flow velocity and T_e/m_e. 
  // u_ie:           bulk flow velocity for ion-electron collision term in ion equation. 
  // vtSq_ie:        squared thermal speed, T_ie/m_i, for ion-electron collision term in ion equation. 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(4,4); 
  Eigen::VectorXd bEV(4); 
  Eigen::VectorXd xEV(4); 
 
  // ....... Compute density ratio through weak division ni/ne .......... // 
  double nRat[4]; 
 
  AEM(0,0) = 0.7071067811865475*ne[0]; 
  AEM(0,1) = 0.7071067811865475*ne[1]; 
  AEM(0,2) = 0.7071067811865475*ne[2]; 
  AEM(0,3) = 0.7071067811865475*ne[3]; 
  AEM(1,0) = 0.7071067811865475*ne[1]; 
  AEM(1,1) = 0.6324555320336759*ne[2]+0.7071067811865475*ne[0]; 
  AEM(1,2) = 0.6210590034081186*ne[3]+0.6324555320336759*ne[1]; 
  AEM(1,3) = 0.6210590034081186*ne[2]; 
  AEM(2,0) = 0.7071067811865475*ne[2]; 
  AEM(2,1) = 0.6210590034081186*ne[3]+0.6324555320336759*ne[1]; 
  AEM(2,2) = 0.4517539514526256*ne[2]+0.7071067811865475*ne[0]; 
  AEM(2,3) = 0.421637021355784*ne[3]+0.6210590034081186*ne[1]; 
  AEM(3,0) = 0.7071067811865475*ne[3]; 
  AEM(3,1) = 0.6210590034081186*ne[2]; 
  AEM(3,2) = 0.421637021355784*ne[3]+0.6210590034081186*ne[1]; 
  AEM(3,3) = 0.421637021355784*ne[2]+0.7071067811865475*ne[0]; 
 
  bEV << ni[0],ni[1],ni[2],ni[3]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(nRat,4,1) = xEV; 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uei[0] = (-0.7071067811865475*nRat[3]*ui[a0+3])+0.7071067811865475*nRat[3]*ue[a0+3]-0.7071067811865475*nRat[2]*ui[a0+2]+0.7071067811865475*nRat[2]*ue[a0+2]-0.7071067811865475*nRat[1]*ui[a0+1]+0.7071067811865475*nRat[1]*ue[a0+1]-0.7071067811865475*nRat[0]*ui[a0]+ui[a0]+0.7071067811865475*nRat[0]*ue[a0]; 
    uei[1] = (-0.6210590034081186*nRat[2]*ui[a0+3])+0.6210590034081186*nRat[2]*ue[a0+3]-0.6210590034081186*nRat[3]*ui[a0+2]-0.6324555320336759*nRat[1]*ui[a0+2]+0.6210590034081186*nRat[3]*ue[a0+2]+0.6324555320336759*nRat[1]*ue[a0+2]-0.6324555320336759*nRat[2]*ui[a0+1]-0.7071067811865475*nRat[0]*ui[a0+1]+ui[a0+1]+0.6324555320336759*nRat[2]*ue[a0+1]+0.7071067811865475*nRat[0]*ue[a0+1]-0.7071067811865475*nRat[1]*ui[a0]+0.7071067811865475*nRat[1]*ue[a0]; 
    uei[2] = (-0.421637021355784*nRat[3]*ui[a0+3])-0.6210590034081186*nRat[1]*ui[a0+3]+0.421637021355784*nRat[3]*ue[a0+3]+0.6210590034081186*nRat[1]*ue[a0+3]-0.4517539514526256*nRat[2]*ui[a0+2]-0.7071067811865475*nRat[0]*ui[a0+2]+ui[a0+2]+0.4517539514526256*nRat[2]*ue[a0+2]+0.7071067811865475*nRat[0]*ue[a0+2]-0.6210590034081186*nRat[3]*ui[a0+1]-0.6324555320336759*nRat[1]*ui[a0+1]+0.6210590034081186*nRat[3]*ue[a0+1]+0.6324555320336759*nRat[1]*ue[a0+1]-0.7071067811865475*nRat[2]*ui[a0]+0.7071067811865475*nRat[2]*ue[a0]; 
    uei[3] = (-0.421637021355784*nRat[2]*ui[a0+3])-0.7071067811865475*nRat[0]*ui[a0+3]+ui[a0+3]+0.421637021355784*nRat[2]*ue[a0+3]+0.7071067811865475*nRat[0]*ue[a0+3]-0.421637021355784*nRat[3]*ui[a0+2]-0.6210590034081186*nRat[1]*ui[a0+2]+0.421637021355784*nRat[3]*ue[a0+2]+0.6210590034081186*nRat[1]*ue[a0+2]-0.6210590034081186*nRat[2]*ui[a0+1]+0.6210590034081186*nRat[2]*ue[a0+1]-0.7071067811865475*nRat[3]*ui[a0]+0.7071067811865475*nRat[3]*ue[a0]; 
 
  } 
 
  // ..... Get the relative speed squared (ue-ui)^2 ..... // 
  double uRelSq[4]; 
  for (unsigned short int k=0; k<4; k++) 
  { 
    uRelSq[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uRelSq[0] += 0.7071067811865475*ui[a0+3]*ui[a0+3]-1.414213562373095*ue[a0+3]*ui[a0+3]+0.7071067811865475*ue[a0+3]*ue[a0+3]+0.7071067811865475*ui[a0+2]*ui[a0+2]-1.414213562373095*ue[a0+2]*ui[a0+2]+0.7071067811865475*ue[a0+2]*ue[a0+2]+0.7071067811865475*ui[a0+1]*ui[a0+1]-1.414213562373095*ue[a0+1]*ui[a0+1]+0.7071067811865475*ue[a0+1]*ue[a0+1]+0.7071067811865475*ui[a0]*ui[a0]-1.414213562373095*ue[a0]*ui[a0]+0.7071067811865475*ue[a0]*ue[a0]; 
    uRelSq[1] += 1.242118006816237*ui[a0+2]*ui[a0+3]-1.242118006816237*ue[a0+2]*ui[a0+3]-1.242118006816237*ui[a0+2]*ue[a0+3]+1.242118006816237*ue[a0+2]*ue[a0+3]+1.264911064067352*ui[a0+1]*ui[a0+2]-1.264911064067352*ue[a0+1]*ui[a0+2]-1.264911064067352*ui[a0+1]*ue[a0+2]+1.264911064067352*ue[a0+1]*ue[a0+2]+1.414213562373095*ui[a0]*ui[a0+1]-1.414213562373095*ue[a0]*ui[a0+1]-1.414213562373095*ui[a0]*ue[a0+1]+1.414213562373095*ue[a0]*ue[a0+1]; 
    uRelSq[2] += 0.421637021355784*ui[a0+3]*ui[a0+3]-0.8432740427115681*ue[a0+3]*ui[a0+3]+1.242118006816237*ui[a0+1]*ui[a0+3]-1.242118006816237*ue[a0+1]*ui[a0+3]+0.421637021355784*ue[a0+3]*ue[a0+3]-1.242118006816237*ui[a0+1]*ue[a0+3]+1.242118006816237*ue[a0+1]*ue[a0+3]+0.4517539514526256*ui[a0+2]*ui[a0+2]-0.9035079029052515*ue[a0+2]*ui[a0+2]+1.414213562373095*ui[a0]*ui[a0+2]-1.414213562373095*ue[a0]*ui[a0+2]+0.4517539514526256*ue[a0+2]*ue[a0+2]-1.414213562373095*ui[a0]*ue[a0+2]+1.414213562373095*ue[a0]*ue[a0+2]+0.6324555320336759*ui[a0+1]*ui[a0+1]-1.264911064067352*ue[a0+1]*ui[a0+1]+0.6324555320336759*ue[a0+1]*ue[a0+1]; 
    uRelSq[3] += 0.8432740427115681*ui[a0+2]*ui[a0+3]-0.8432740427115681*ue[a0+2]*ui[a0+3]+1.414213562373095*ui[a0]*ui[a0+3]-1.414213562373095*ue[a0]*ui[a0+3]-0.8432740427115681*ui[a0+2]*ue[a0+3]+0.8432740427115681*ue[a0+2]*ue[a0+3]-1.414213562373095*ui[a0]*ue[a0+3]+1.414213562373095*ue[a0]*ue[a0+3]+1.242118006816237*ui[a0+1]*ui[a0+2]-1.242118006816237*ue[a0+1]*ui[a0+2]-1.242118006816237*ui[a0+1]*ue[a0+2]+1.242118006816237*ue[a0+1]*ue[a0+2]; 
  } 
 
  double rmRat = 1.0/mRat; 
  vtSqie[0] = (1.414213562373095*nRat[3]*vtSqe[3]+1.414213562373095*nRat[2]*vtSqe[2]+1.414213562373095*nRat[1]*vtSqe[1]+1.414213562373095*nRat[0]*vtSqe[0])*rmRat-1.414213562373095*nRat[3]*vtSqi[3]-0.1178511301977579*nRat[3]*uRelSq[3]-1.414213562373095*nRat[2]*vtSqi[2]-0.1178511301977579*nRat[2]*uRelSq[2]-1.414213562373095*nRat[1]*vtSqi[1]-0.1178511301977579*nRat[1]*uRelSq[1]-1.414213562373095*nRat[0]*vtSqi[0]+vtSqi[0]-0.1178511301977579*nRat[0]*uRelSq[0]; 
  vtSqie[1] = (1.242118006816237*nRat[2]*vtSqe[3]+1.242118006816237*vtSqe[2]*nRat[3]+1.264911064067352*nRat[1]*vtSqe[2]+1.264911064067352*vtSqe[1]*nRat[2]+1.414213562373095*nRat[0]*vtSqe[1]+1.414213562373095*vtSqe[0]*nRat[1])*rmRat-1.242118006816237*nRat[2]*vtSqi[3]-0.1035098339013531*nRat[2]*uRelSq[3]-1.242118006816237*vtSqi[2]*nRat[3]-0.1035098339013531*uRelSq[2]*nRat[3]-1.264911064067352*nRat[1]*vtSqi[2]-0.105409255338946*nRat[1]*uRelSq[2]-1.264911064067352*vtSqi[1]*nRat[2]-0.105409255338946*uRelSq[1]*nRat[2]-1.414213562373095*nRat[0]*vtSqi[1]+vtSqi[1]-0.1178511301977579*nRat[0]*uRelSq[1]-1.414213562373095*vtSqi[0]*nRat[1]-0.1178511301977579*uRelSq[0]*nRat[1]; 
  vtSqie[2] = (0.8432740427115681*nRat[3]*vtSqe[3]+1.242118006816237*nRat[1]*vtSqe[3]+1.242118006816237*vtSqe[1]*nRat[3]+0.9035079029052515*nRat[2]*vtSqe[2]+1.414213562373095*nRat[0]*vtSqe[2]+1.414213562373095*vtSqe[0]*nRat[2]+1.264911064067352*nRat[1]*vtSqe[1])*rmRat-0.8432740427115681*nRat[3]*vtSqi[3]-1.242118006816237*nRat[1]*vtSqi[3]-0.07027283689263064*nRat[3]*uRelSq[3]-0.1035098339013531*nRat[1]*uRelSq[3]-1.242118006816237*vtSqi[1]*nRat[3]-0.1035098339013531*uRelSq[1]*nRat[3]-0.9035079029052515*nRat[2]*vtSqi[2]-1.414213562373095*nRat[0]*vtSqi[2]+vtSqi[2]-0.07529232524210427*nRat[2]*uRelSq[2]-0.1178511301977579*nRat[0]*uRelSq[2]-1.414213562373095*vtSqi[0]*nRat[2]-0.1178511301977579*uRelSq[0]*nRat[2]-1.264911064067352*nRat[1]*vtSqi[1]-0.105409255338946*nRat[1]*uRelSq[1]; 
  vtSqie[3] = (0.8432740427115681*nRat[2]*vtSqe[3]+1.414213562373095*nRat[0]*vtSqe[3]+0.8432740427115681*vtSqe[2]*nRat[3]+1.414213562373095*vtSqe[0]*nRat[3]+1.242118006816237*nRat[1]*vtSqe[2]+1.242118006816237*vtSqe[1]*nRat[2])*rmRat-0.8432740427115681*nRat[2]*vtSqi[3]-1.414213562373095*nRat[0]*vtSqi[3]+vtSqi[3]-0.07027283689263064*nRat[2]*uRelSq[3]-0.1178511301977579*nRat[0]*uRelSq[3]-0.8432740427115681*vtSqi[2]*nRat[3]-0.07027283689263064*uRelSq[2]*nRat[3]-1.414213562373095*vtSqi[0]*nRat[3]-0.1178511301977579*uRelSq[0]*nRat[3]-1.242118006816237*nRat[1]*vtSqi[2]-0.1035098339013531*nRat[1]*uRelSq[2]-1.242118006816237*vtSqi[1]*nRat[2]-0.1035098339013531*uRelSq[1]*nRat[2]; 
 
} 
 
void CrossPrimMoments_VmieLBO_1x2vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uei, double *vtSqie) 
{ 
  // mRat:           mass ratio = m_i/m_e. 
  // ne, ue, vtSqe:  electron number density, bulk flow velocity and T_e/m_e. 
  // ni, ui, vtSqi:  ion number density, bulk flow velocity and T_e/m_e. 
  // u_ie:           bulk flow velocity for ion-electron collision term in ion equation. 
  // vtSq_ie:        squared thermal speed, T_ie/m_i, for ion-electron collision term in ion equation. 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(5,5); 
  Eigen::VectorXd bEV(5); 
  Eigen::VectorXd xEV(5); 
 
  // ....... Compute density ratio through weak division ni/ne .......... // 
  double nRat[5]; 
 
  AEM(0,0) = 0.7071067811865475*ne[0]; 
  AEM(0,1) = 0.7071067811865475*ne[1]; 
  AEM(0,2) = 0.7071067811865475*ne[2]; 
  AEM(0,3) = 0.7071067811865475*ne[3]; 
  AEM(0,4) = 0.7071067811865475*ne[4]; 
  AEM(1,0) = 0.7071067811865475*ne[1]; 
  AEM(1,1) = 0.6324555320336759*ne[2]+0.7071067811865475*ne[0]; 
  AEM(1,2) = 0.6210590034081186*ne[3]+0.6324555320336759*ne[1]; 
  AEM(1,3) = 0.6172133998483679*ne[4]+0.6210590034081186*ne[2]; 
  AEM(1,4) = 0.6172133998483679*ne[3]; 
  AEM(2,0) = 0.7071067811865475*ne[2]; 
  AEM(2,1) = 0.6210590034081186*ne[3]+0.6324555320336759*ne[1]; 
  AEM(2,2) = 0.6060915267313265*ne[4]+0.4517539514526256*ne[2]+0.7071067811865475*ne[0]; 
  AEM(2,3) = 0.421637021355784*ne[3]+0.6210590034081186*ne[1]; 
  AEM(2,4) = 0.410685410411478*ne[4]+0.6060915267313265*ne[2]; 
  AEM(3,0) = 0.7071067811865475*ne[3]; 
  AEM(3,1) = 0.6172133998483679*ne[4]+0.6210590034081186*ne[2]; 
  AEM(3,2) = 0.421637021355784*ne[3]+0.6210590034081186*ne[1]; 
  AEM(3,3) = 0.385694607919935*ne[4]+0.421637021355784*ne[2]+0.7071067811865475*ne[0]; 
  AEM(3,4) = 0.385694607919935*ne[3]+0.6172133998483679*ne[1]; 
  AEM(4,0) = 0.7071067811865475*ne[4]; 
  AEM(4,1) = 0.6172133998483679*ne[3]; 
  AEM(4,2) = 0.410685410411478*ne[4]+0.6060915267313265*ne[2]; 
  AEM(4,3) = 0.385694607919935*ne[3]+0.6172133998483679*ne[1]; 
  AEM(4,4) = 0.3433105850715905*ne[4]+0.410685410411478*ne[2]+0.7071067811865475*ne[0]; 
 
  bEV << ni[0],ni[1],ni[2],ni[3],ni[4]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(nRat,5,1) = xEV; 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 5*vd; 
    uei[0] = (-0.7071067811865475*nRat[4]*ui[a0+4])+0.7071067811865475*nRat[4]*ue[a0+4]-0.7071067811865475*nRat[3]*ui[a0+3]+0.7071067811865475*nRat[3]*ue[a0+3]-0.7071067811865475*nRat[2]*ui[a0+2]+0.7071067811865475*nRat[2]*ue[a0+2]-0.7071067811865475*nRat[1]*ui[a0+1]+0.7071067811865475*nRat[1]*ue[a0+1]-0.7071067811865475*nRat[0]*ui[a0]+ui[a0]+0.7071067811865475*nRat[0]*ue[a0]; 
    uei[1] = (-0.6172133998483679*nRat[3]*ui[a0+4])+0.6172133998483679*nRat[3]*ue[a0+4]-0.6172133998483679*nRat[4]*ui[a0+3]-0.6210590034081186*nRat[2]*ui[a0+3]+0.6172133998483679*nRat[4]*ue[a0+3]+0.6210590034081186*nRat[2]*ue[a0+3]-0.6210590034081186*nRat[3]*ui[a0+2]-0.6324555320336759*nRat[1]*ui[a0+2]+0.6210590034081186*nRat[3]*ue[a0+2]+0.6324555320336759*nRat[1]*ue[a0+2]-0.6324555320336759*nRat[2]*ui[a0+1]-0.7071067811865475*nRat[0]*ui[a0+1]+ui[a0+1]+0.6324555320336759*nRat[2]*ue[a0+1]+0.7071067811865475*nRat[0]*ue[a0+1]-0.7071067811865475*nRat[1]*ui[a0]+0.7071067811865475*nRat[1]*ue[a0]; 
    uei[2] = (-0.410685410411478*nRat[4]*ui[a0+4])-0.6060915267313265*nRat[2]*ui[a0+4]+0.410685410411478*nRat[4]*ue[a0+4]+0.6060915267313265*nRat[2]*ue[a0+4]-0.421637021355784*nRat[3]*ui[a0+3]-0.6210590034081186*nRat[1]*ui[a0+3]+0.421637021355784*nRat[3]*ue[a0+3]+0.6210590034081186*nRat[1]*ue[a0+3]-0.6060915267313265*nRat[4]*ui[a0+2]-0.4517539514526256*nRat[2]*ui[a0+2]-0.7071067811865475*nRat[0]*ui[a0+2]+ui[a0+2]+0.6060915267313265*nRat[4]*ue[a0+2]+0.4517539514526256*nRat[2]*ue[a0+2]+0.7071067811865475*nRat[0]*ue[a0+2]-0.6210590034081186*nRat[3]*ui[a0+1]-0.6324555320336759*nRat[1]*ui[a0+1]+0.6210590034081186*nRat[3]*ue[a0+1]+0.6324555320336759*nRat[1]*ue[a0+1]-0.7071067811865475*nRat[2]*ui[a0]+0.7071067811865475*nRat[2]*ue[a0]; 
    uei[3] = (-0.385694607919935*nRat[3]*ui[a0+4])-0.6172133998483679*nRat[1]*ui[a0+4]+0.385694607919935*nRat[3]*ue[a0+4]+0.6172133998483679*nRat[1]*ue[a0+4]-0.385694607919935*nRat[4]*ui[a0+3]-0.421637021355784*nRat[2]*ui[a0+3]-0.7071067811865475*nRat[0]*ui[a0+3]+ui[a0+3]+0.385694607919935*nRat[4]*ue[a0+3]+0.421637021355784*nRat[2]*ue[a0+3]+0.7071067811865475*nRat[0]*ue[a0+3]-0.421637021355784*nRat[3]*ui[a0+2]-0.6210590034081186*nRat[1]*ui[a0+2]+0.421637021355784*nRat[3]*ue[a0+2]+0.6210590034081186*nRat[1]*ue[a0+2]-0.6172133998483679*nRat[4]*ui[a0+1]-0.6210590034081186*nRat[2]*ui[a0+1]+0.6172133998483679*nRat[4]*ue[a0+1]+0.6210590034081186*nRat[2]*ue[a0+1]-0.7071067811865475*nRat[3]*ui[a0]+0.7071067811865475*nRat[3]*ue[a0]; 
    uei[4] = (-0.3433105850715905*nRat[4]*ui[a0+4])-0.410685410411478*nRat[2]*ui[a0+4]-0.7071067811865475*nRat[0]*ui[a0+4]+ui[a0+4]+0.3433105850715905*nRat[4]*ue[a0+4]+0.410685410411478*nRat[2]*ue[a0+4]+0.7071067811865475*nRat[0]*ue[a0+4]-0.385694607919935*nRat[3]*ui[a0+3]-0.6172133998483679*nRat[1]*ui[a0+3]+0.385694607919935*nRat[3]*ue[a0+3]+0.6172133998483679*nRat[1]*ue[a0+3]-0.410685410411478*nRat[4]*ui[a0+2]-0.6060915267313265*nRat[2]*ui[a0+2]+0.410685410411478*nRat[4]*ue[a0+2]+0.6060915267313265*nRat[2]*ue[a0+2]-0.6172133998483679*nRat[3]*ui[a0+1]+0.6172133998483679*nRat[3]*ue[a0+1]-0.7071067811865475*nRat[4]*ui[a0]+0.7071067811865475*nRat[4]*ue[a0]; 
 
  } 
 
  // ..... Get the relative speed squared (ue-ui)^2 ..... // 
  double uRelSq[5]; 
  for (unsigned short int k=0; k<5; k++) 
  { 
    uRelSq[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 5*vd; 
    uRelSq[0] += 0.7071067811865475*ui[a0+4]*ui[a0+4]-1.414213562373095*ue[a0+4]*ui[a0+4]+0.7071067811865475*ue[a0+4]*ue[a0+4]+0.7071067811865475*ui[a0+3]*ui[a0+3]-1.414213562373095*ue[a0+3]*ui[a0+3]+0.7071067811865475*ue[a0+3]*ue[a0+3]+0.7071067811865475*ui[a0+2]*ui[a0+2]-1.414213562373095*ue[a0+2]*ui[a0+2]+0.7071067811865475*ue[a0+2]*ue[a0+2]+0.7071067811865475*ui[a0+1]*ui[a0+1]-1.414213562373095*ue[a0+1]*ui[a0+1]+0.7071067811865475*ue[a0+1]*ue[a0+1]+0.7071067811865475*ui[a0]*ui[a0]-1.414213562373095*ue[a0]*ui[a0]+0.7071067811865475*ue[a0]*ue[a0]; 
    uRelSq[1] += 1.234426799696736*ui[a0+3]*ui[a0+4]-1.234426799696736*ue[a0+3]*ui[a0+4]-1.234426799696736*ui[a0+3]*ue[a0+4]+1.234426799696736*ue[a0+3]*ue[a0+4]+1.242118006816237*ui[a0+2]*ui[a0+3]-1.242118006816237*ue[a0+2]*ui[a0+3]-1.242118006816237*ui[a0+2]*ue[a0+3]+1.242118006816237*ue[a0+2]*ue[a0+3]+1.264911064067352*ui[a0+1]*ui[a0+2]-1.264911064067352*ue[a0+1]*ui[a0+2]-1.264911064067352*ui[a0+1]*ue[a0+2]+1.264911064067352*ue[a0+1]*ue[a0+2]+1.414213562373095*ui[a0]*ui[a0+1]-1.414213562373095*ue[a0]*ui[a0+1]-1.414213562373095*ui[a0]*ue[a0+1]+1.414213562373095*ue[a0]*ue[a0+1]; 
    uRelSq[2] += 0.410685410411478*ui[a0+4]*ui[a0+4]-0.8213708208229562*ue[a0+4]*ui[a0+4]+1.212183053462653*ui[a0+2]*ui[a0+4]-1.212183053462653*ue[a0+2]*ui[a0+4]+0.410685410411478*ue[a0+4]*ue[a0+4]-1.212183053462653*ui[a0+2]*ue[a0+4]+1.212183053462653*ue[a0+2]*ue[a0+4]+0.421637021355784*ui[a0+3]*ui[a0+3]-0.8432740427115681*ue[a0+3]*ui[a0+3]+1.242118006816237*ui[a0+1]*ui[a0+3]-1.242118006816237*ue[a0+1]*ui[a0+3]+0.421637021355784*ue[a0+3]*ue[a0+3]-1.242118006816237*ui[a0+1]*ue[a0+3]+1.242118006816237*ue[a0+1]*ue[a0+3]+0.4517539514526256*ui[a0+2]*ui[a0+2]-0.9035079029052515*ue[a0+2]*ui[a0+2]+1.414213562373095*ui[a0]*ui[a0+2]-1.414213562373095*ue[a0]*ui[a0+2]+0.4517539514526256*ue[a0+2]*ue[a0+2]-1.414213562373095*ui[a0]*ue[a0+2]+1.414213562373095*ue[a0]*ue[a0+2]+0.6324555320336759*ui[a0+1]*ui[a0+1]-1.264911064067352*ue[a0+1]*ui[a0+1]+0.6324555320336759*ue[a0+1]*ue[a0+1]; 
    uRelSq[3] += 0.7713892158398702*ui[a0+3]*ui[a0+4]-0.7713892158398702*ue[a0+3]*ui[a0+4]+1.234426799696736*ui[a0+1]*ui[a0+4]-1.234426799696736*ue[a0+1]*ui[a0+4]-0.7713892158398702*ui[a0+3]*ue[a0+4]+0.7713892158398702*ue[a0+3]*ue[a0+4]-1.234426799696736*ui[a0+1]*ue[a0+4]+1.234426799696736*ue[a0+1]*ue[a0+4]+0.8432740427115681*ui[a0+2]*ui[a0+3]-0.8432740427115681*ue[a0+2]*ui[a0+3]+1.414213562373095*ui[a0]*ui[a0+3]-1.414213562373095*ue[a0]*ui[a0+3]-0.8432740427115681*ui[a0+2]*ue[a0+3]+0.8432740427115681*ue[a0+2]*ue[a0+3]-1.414213562373095*ui[a0]*ue[a0+3]+1.414213562373095*ue[a0]*ue[a0+3]+1.242118006816237*ui[a0+1]*ui[a0+2]-1.242118006816237*ue[a0+1]*ui[a0+2]-1.242118006816237*ui[a0+1]*ue[a0+2]+1.242118006816237*ue[a0+1]*ue[a0+2]; 
    uRelSq[4] += 0.3433105850715905*ui[a0+4]*ui[a0+4]-0.6866211701431811*ue[a0+4]*ui[a0+4]+0.8213708208229562*ui[a0+2]*ui[a0+4]-0.8213708208229562*ue[a0+2]*ui[a0+4]+1.414213562373095*ui[a0]*ui[a0+4]-1.414213562373095*ue[a0]*ui[a0+4]+0.3433105850715905*ue[a0+4]*ue[a0+4]-0.8213708208229562*ui[a0+2]*ue[a0+4]+0.8213708208229562*ue[a0+2]*ue[a0+4]-1.414213562373095*ui[a0]*ue[a0+4]+1.414213562373095*ue[a0]*ue[a0+4]+0.385694607919935*ui[a0+3]*ui[a0+3]-0.7713892158398702*ue[a0+3]*ui[a0+3]+1.234426799696736*ui[a0+1]*ui[a0+3]-1.234426799696736*ue[a0+1]*ui[a0+3]+0.385694607919935*ue[a0+3]*ue[a0+3]-1.234426799696736*ui[a0+1]*ue[a0+3]+1.234426799696736*ue[a0+1]*ue[a0+3]+0.6060915267313265*ui[a0+2]*ui[a0+2]-1.212183053462653*ue[a0+2]*ui[a0+2]+0.6060915267313265*ue[a0+2]*ue[a0+2]; 
  } 
 
  double rmRat = 1.0/mRat; 
  vtSqie[0] = (1.414213562373095*nRat[4]*vtSqe[4]+1.414213562373095*nRat[3]*vtSqe[3]+1.414213562373095*nRat[2]*vtSqe[2]+1.414213562373095*nRat[1]*vtSqe[1]+1.414213562373095*nRat[0]*vtSqe[0])*rmRat-1.414213562373095*nRat[4]*vtSqi[4]-0.1178511301977579*nRat[4]*uRelSq[4]-1.414213562373095*nRat[3]*vtSqi[3]-0.1178511301977579*nRat[3]*uRelSq[3]-1.414213562373095*nRat[2]*vtSqi[2]-0.1178511301977579*nRat[2]*uRelSq[2]-1.414213562373095*nRat[1]*vtSqi[1]-0.1178511301977579*nRat[1]*uRelSq[1]-1.414213562373095*nRat[0]*vtSqi[0]+vtSqi[0]-0.1178511301977579*nRat[0]*uRelSq[0]; 
  vtSqie[1] = (1.234426799696736*nRat[3]*vtSqe[4]+1.234426799696736*vtSqe[3]*nRat[4]+1.242118006816237*nRat[2]*vtSqe[3]+1.242118006816237*vtSqe[2]*nRat[3]+1.264911064067352*nRat[1]*vtSqe[2]+1.264911064067352*vtSqe[1]*nRat[2]+1.414213562373095*nRat[0]*vtSqe[1]+1.414213562373095*vtSqe[0]*nRat[1])*rmRat-1.234426799696736*nRat[3]*vtSqi[4]-0.102868899974728*nRat[3]*uRelSq[4]-1.234426799696736*vtSqi[3]*nRat[4]-0.102868899974728*uRelSq[3]*nRat[4]-1.242118006816237*nRat[2]*vtSqi[3]-0.1035098339013531*nRat[2]*uRelSq[3]-1.242118006816237*vtSqi[2]*nRat[3]-0.1035098339013531*uRelSq[2]*nRat[3]-1.264911064067352*nRat[1]*vtSqi[2]-0.105409255338946*nRat[1]*uRelSq[2]-1.264911064067352*vtSqi[1]*nRat[2]-0.105409255338946*uRelSq[1]*nRat[2]-1.414213562373095*nRat[0]*vtSqi[1]+vtSqi[1]-0.1178511301977579*nRat[0]*uRelSq[1]-1.414213562373095*vtSqi[0]*nRat[1]-0.1178511301977579*uRelSq[0]*nRat[1]; 
  vtSqie[2] = (0.8213708208229562*nRat[4]*vtSqe[4]+1.212183053462653*nRat[2]*vtSqe[4]+1.212183053462653*vtSqe[2]*nRat[4]+0.8432740427115681*nRat[3]*vtSqe[3]+1.242118006816237*nRat[1]*vtSqe[3]+1.242118006816237*vtSqe[1]*nRat[3]+0.9035079029052515*nRat[2]*vtSqe[2]+1.414213562373095*nRat[0]*vtSqe[2]+1.414213562373095*vtSqe[0]*nRat[2]+1.264911064067352*nRat[1]*vtSqe[1])*rmRat-0.8213708208229562*nRat[4]*vtSqi[4]-1.212183053462653*nRat[2]*vtSqi[4]-0.06844756840191299*nRat[4]*uRelSq[4]-0.1010152544552211*nRat[2]*uRelSq[4]-1.212183053462653*vtSqi[2]*nRat[4]-0.1010152544552211*uRelSq[2]*nRat[4]-0.8432740427115681*nRat[3]*vtSqi[3]-1.242118006816237*nRat[1]*vtSqi[3]-0.07027283689263064*nRat[3]*uRelSq[3]-0.1035098339013531*nRat[1]*uRelSq[3]-1.242118006816237*vtSqi[1]*nRat[3]-0.1035098339013531*uRelSq[1]*nRat[3]-0.9035079029052515*nRat[2]*vtSqi[2]-1.414213562373095*nRat[0]*vtSqi[2]+vtSqi[2]-0.07529232524210427*nRat[2]*uRelSq[2]-0.1178511301977579*nRat[0]*uRelSq[2]-1.414213562373095*vtSqi[0]*nRat[2]-0.1178511301977579*uRelSq[0]*nRat[2]-1.264911064067352*nRat[1]*vtSqi[1]-0.105409255338946*nRat[1]*uRelSq[1]; 
  vtSqie[3] = (0.7713892158398702*nRat[3]*vtSqe[4]+1.234426799696736*nRat[1]*vtSqe[4]+0.7713892158398702*vtSqe[3]*nRat[4]+1.234426799696736*vtSqe[1]*nRat[4]+0.8432740427115681*nRat[2]*vtSqe[3]+1.414213562373095*nRat[0]*vtSqe[3]+0.8432740427115681*vtSqe[2]*nRat[3]+1.414213562373095*vtSqe[0]*nRat[3]+1.242118006816237*nRat[1]*vtSqe[2]+1.242118006816237*vtSqe[1]*nRat[2])*rmRat-0.7713892158398702*nRat[3]*vtSqi[4]-1.234426799696736*nRat[1]*vtSqi[4]-0.06428243465332249*nRat[3]*uRelSq[4]-0.102868899974728*nRat[1]*uRelSq[4]-0.7713892158398702*vtSqi[3]*nRat[4]-0.06428243465332249*uRelSq[3]*nRat[4]-1.234426799696736*vtSqi[1]*nRat[4]-0.102868899974728*uRelSq[1]*nRat[4]-0.8432740427115681*nRat[2]*vtSqi[3]-1.414213562373095*nRat[0]*vtSqi[3]+vtSqi[3]-0.07027283689263064*nRat[2]*uRelSq[3]-0.1178511301977579*nRat[0]*uRelSq[3]-0.8432740427115681*vtSqi[2]*nRat[3]-0.07027283689263064*uRelSq[2]*nRat[3]-1.414213562373095*vtSqi[0]*nRat[3]-0.1178511301977579*uRelSq[0]*nRat[3]-1.242118006816237*nRat[1]*vtSqi[2]-0.1035098339013531*nRat[1]*uRelSq[2]-1.242118006816237*vtSqi[1]*nRat[2]-0.1035098339013531*uRelSq[1]*nRat[2]; 
  vtSqie[4] = (0.6866211701431811*nRat[4]*vtSqe[4]+0.8213708208229562*nRat[2]*vtSqe[4]+1.414213562373095*nRat[0]*vtSqe[4]+0.8213708208229562*vtSqe[2]*nRat[4]+1.414213562373095*vtSqe[0]*nRat[4]+0.7713892158398702*nRat[3]*vtSqe[3]+1.234426799696736*nRat[1]*vtSqe[3]+1.234426799696736*vtSqe[1]*nRat[3]+1.212183053462653*nRat[2]*vtSqe[2])*rmRat-0.6866211701431811*nRat[4]*vtSqi[4]-0.8213708208229562*nRat[2]*vtSqi[4]-1.414213562373095*nRat[0]*vtSqi[4]+vtSqi[4]-0.05721843084526507*nRat[4]*uRelSq[4]-0.06844756840191299*nRat[2]*uRelSq[4]-0.1178511301977579*nRat[0]*uRelSq[4]-0.8213708208229562*vtSqi[2]*nRat[4]-0.06844756840191299*uRelSq[2]*nRat[4]-1.414213562373095*vtSqi[0]*nRat[4]-0.1178511301977579*uRelSq[0]*nRat[4]-0.7713892158398702*nRat[3]*vtSqi[3]-1.234426799696736*nRat[1]*vtSqi[3]-0.06428243465332249*nRat[3]*uRelSq[3]-0.102868899974728*nRat[1]*uRelSq[3]-1.234426799696736*vtSqi[1]*nRat[3]-0.102868899974728*uRelSq[1]*nRat[3]-1.212183053462653*nRat[2]*vtSqi[2]-0.1010152544552211*nRat[2]*uRelSq[2]; 
 
} 
 
