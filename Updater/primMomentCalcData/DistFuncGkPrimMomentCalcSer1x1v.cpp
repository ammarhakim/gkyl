#include <math.h> 
#include <DistFuncPrimMomentCalcModDecl.h> 
 
using namespace Eigen; 
 
void GkPrimMomentCalc1x1vSer_Upar_P1(const double *m0, const double *m1, const double *m2, double *out) 
{ 
  // Define Eigen Matrix with triple basis tensor dotted with m0 vector. 
  Eigen::MatrixXd m0EM(2,2); 
  m0EM(0,0) = 0.7071067811865475*m0[0]; 
  m0EM(0,1) = 0.7071067811865475*m0[1]; 
  m0EM(1,0) = 0.7071067811865475*m0[1]; 
  m0EM(1,1) = 0.7071067811865475*m0[0]; 
 
  // Define Eigen Vector with coefficients of m1. 
  Eigen::Map<const VectorXd> m1EV(m1,2); 
  ::new (&m1EV) Eigen::Map<const VectorXd>(m1,2); 
 
  // Solve the system of equations. 
  Eigen::VectorXd parU(2); 
  parU = m0EM.colPivHouseholderQr().solve(m1EV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,2,1) = parU; 
 
} 
