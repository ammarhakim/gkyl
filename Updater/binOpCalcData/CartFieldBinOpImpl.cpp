#include <CartFieldBinOpModDecl.h> 

binOpData_t::binOpData_t(int nbasis_S, int nbasis_D) { 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  AEM_S = Eigen::MatrixXd::Zero(nbasis_S,nbasis_S); 
  // Declare Eigen Vector with coefficients of B. 
  BEV_S = Eigen::VectorXd::Zero(nbasis_S);  
  // Declare vector with solution to system of equations. 
  u_S = Eigen::VectorXd::Zero(nbasis_S);  

  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  AEM_D = Eigen::MatrixXd::Zero(nbasis_D,nbasis_D); 
  // Declare Eigen Vector with coefficients of B. 
  BEV_D = Eigen::VectorXd::Zero(nbasis_D);  
  // Declare vector with solution to system of equations. 
  u_D = Eigen::VectorXd::Zero(nbasis_D);  
} 

extern "C" void* new_binOpData_t(int nbasis_S, int nbasis_D) {
  binOpData_t* b = new binOpData_t(nbasis_S, nbasis_D);
  return reinterpret_cast<void*>(b);
}
