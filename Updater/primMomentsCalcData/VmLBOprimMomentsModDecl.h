#ifndef VM_LBO_PRIMMOMENTS_MOD_DECL_H 
#define VM_LBO_PRIMMOMENTS_MOD_DECL_H 
 
// Eigen include statements. 
#include <Eigen/Dense> 
 
extern "C" { 
void VmLBOconstNuPrimMoments1x1vSer_P1(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, double *u, double *vtSq); 


} 
#endif 
