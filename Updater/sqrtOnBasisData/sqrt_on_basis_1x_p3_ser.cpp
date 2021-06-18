#include <sqrt_on_basis_mod_decl.h>

void sqrt_on_basis_gauss_1x_p3_ser(const double qExp, const double *fIn, double *out) 
{ 
  // qExp: exponent in sqrt(f)^q.
  // fIn:  input field.
  // out:  output field.
double sqrtfRq[4];
  sqrtfRq[0] = pow(sqrt(0.7702725556588816*fIn[3]-0.5164305132317774*fIn[2]-0.416390039500913*fIn[1]+0.7071067811865475*fIn[0]),qExp); 
  sqrtfRq[1] = pow(sqrt((-0.7702725556588816*fIn[3])-0.5164305132317774*fIn[2]+0.416390039500913*fIn[1]+0.7071067811865475*fIn[0]),qExp); 
  sqrtfRq[2] = pow(sqrt(0.5701294036773671*fIn[3]+0.9681844646844028*fIn[2]+1.054672281193885*fIn[1]+0.7071067811865475*fIn[0]),qExp); 
  sqrtfRq[3] = pow(sqrt((-0.5701294036773671*fIn[3])+0.9681844646844028*fIn[2]-1.054672281193885*fIn[1]+0.7071067811865475*fIn[0]),qExp); 

  out[0] = 0.2459705198652899*(sqrtfRq[3]+sqrtfRq[2])+0.4611362613212575*(sqrtfRq[1]+sqrtfRq[0]); 
  out[1] = (-0.3668728630454641*sqrtfRq[3])+0.3668728630454641*sqrtfRq[2]+0.2715467467935446*sqrtfRq[1]-0.2715467467935446*sqrtfRq[0]; 
  out[2] = 0.3367876570272816*(sqrtfRq[3]+sqrtfRq[2])-0.3367876570272816*(sqrtfRq[1]+sqrtfRq[0]); 
  out[3] = (-0.1983222754244995*sqrtfRq[3])+0.1983222754244995*sqrtfRq[2]-0.5023295150965305*sqrtfRq[1]+0.5023295150965305*sqrtfRq[0]; 
}
