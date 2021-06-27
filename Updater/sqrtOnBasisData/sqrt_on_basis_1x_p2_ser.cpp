#include <sqrt_on_basis_mod_decl.h>

void sqrt_on_basis_gauss_1x_p2_ser(const double qExp, const double *fIn, double *out) 
{ 
  // qExp: exponent in sqrt(f)^q.
  // fIn:  input field.
  // out:  output field.
double sqrtfRq[3];
  sqrtfRq[0] = pow(sqrt(0.7071067811865475*fIn[0]-0.7905694150420947*fIn[2]),qExp); 
  sqrtfRq[1] = pow(sqrt(0.6324555320336759*fIn[2]-0.9486832980505137*fIn[1]+0.7071067811865475*fIn[0]),qExp); 
  sqrtfRq[2] = pow(sqrt(0.6324555320336759*fIn[2]+0.9486832980505137*fIn[1]+0.7071067811865475*fIn[0]),qExp); 

  out[0] = 0.392837100659193*(sqrtfRq[2]+sqrtfRq[1])+0.6285393610547089*sqrtfRq[0]; 
  out[1] = 0.5270462766947299*sqrtfRq[2]-0.5270462766947299*sqrtfRq[1]; 
  out[2] = 0.3513641844631533*(sqrtfRq[2]+sqrtfRq[1])-0.7027283689263064*sqrtfRq[0]; 
}
