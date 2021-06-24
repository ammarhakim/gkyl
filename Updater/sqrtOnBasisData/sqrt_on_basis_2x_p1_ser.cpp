#include <sqrt_on_basis_mod_decl.h>

void sqrt_on_basis_gauss_2x_p1_ser(const double qExp, const double *fIn, double *out) 
{ 
  // qExp: exponent in sqrt(f)^q.
  // fIn:  input field.
  // out:  output field.
double sqrtfRq[4];
  sqrtfRq[0] = pow(sqrt(0.4999999999999999*fIn[3]-0.4999999999999999*fIn[2]-0.4999999999999999*fIn[1]+0.5*fIn[0]),qExp); 
  sqrtfRq[1] = pow(sqrt((-0.4999999999999999*fIn[3])+0.4999999999999999*fIn[2]-0.4999999999999999*fIn[1]+0.5*fIn[0]),qExp); 
  sqrtfRq[2] = pow(sqrt((-0.4999999999999999*fIn[3])-0.4999999999999999*fIn[2]+0.4999999999999999*fIn[1]+0.5*fIn[0]),qExp); 
  sqrtfRq[3] = pow(sqrt(0.4999999999999999*fIn[3]+0.4999999999999999*fIn[2]+0.4999999999999999*fIn[1]+0.5*fIn[0]),qExp); 

  out[0] = 0.5*(sqrtfRq[3]+sqrtfRq[2]+sqrtfRq[1]+sqrtfRq[0]); 
  out[1] = 0.4999999999999999*(sqrtfRq[3]+sqrtfRq[2])-0.4999999999999999*(sqrtfRq[1]+sqrtfRq[0]); 
  out[2] = 0.4999999999999999*sqrtfRq[3]-0.4999999999999999*sqrtfRq[2]+0.4999999999999999*sqrtfRq[1]-0.4999999999999999*sqrtfRq[0]; 
  out[3] = 0.4999999999999999*sqrtfRq[3]-0.4999999999999999*(sqrtfRq[2]+sqrtfRq[1])+0.4999999999999999*sqrtfRq[0]; 
}
