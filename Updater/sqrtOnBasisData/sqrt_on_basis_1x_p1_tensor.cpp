#include <sqrt_on_basis_mod_decl.h>

void sqrt_on_basis_gauss_1x_p1_tensor(const double qExp, const double *fIn, double *out) 
{ 
  // qExp: exponent in sqrt(f)^q.
  // fIn:  input field.
  // out:  output field.
double sqrtfRq[2];
  sqrtfRq[0] = pow(sqrt(0.7071067811865475*fIn[0]-0.7071067811865474*fIn[1]),qExp); 
  sqrtfRq[1] = pow(sqrt(0.7071067811865474*fIn[1]+0.7071067811865475*fIn[0]),qExp); 

  out[0] = 0.7071067811865475*(sqrtfRq[1]+sqrtfRq[0]); 
  out[1] = 0.7071067811865474*sqrtfRq[1]-0.7071067811865474*sqrtfRq[0]; 
}
