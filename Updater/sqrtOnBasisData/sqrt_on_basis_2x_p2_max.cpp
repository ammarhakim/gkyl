#include <sqrt_on_basis_mod_decl.h>

void sqrt_on_basis_gauss_2x_p2_max(const double qExp, const double *fIn, double *out) 
{ 
  // qExp: exponent in sqrt(f)^q.
  // fIn:  input field.
  // out:  output field.
double sqrtfRq[9];
  sqrtfRq[0] = pow(sqrt((-0.5590169943749475*fIn[5])-0.5590169943749475*fIn[4]+0.5*fIn[0]),qExp); 
  sqrtfRq[1] = pow(sqrt(0.4472135954999581*fIn[5]-0.5590169943749475*fIn[4]-0.6708203932499369*fIn[2]+0.5*fIn[0]),qExp); 
  sqrtfRq[2] = pow(sqrt(0.4472135954999581*fIn[5]-0.5590169943749475*fIn[4]+0.6708203932499369*fIn[2]+0.5*fIn[0]),qExp); 
  sqrtfRq[3] = pow(sqrt((-0.5590169943749475*fIn[5])+0.4472135954999581*fIn[4]-0.6708203932499369*fIn[1]+0.5*fIn[0]),qExp); 
  sqrtfRq[4] = pow(sqrt(0.4472135954999581*fIn[5]+0.4472135954999581*fIn[4]+0.9*fIn[3]-0.6708203932499369*fIn[2]-0.6708203932499369*fIn[1]+0.5*fIn[0]),qExp); 
  sqrtfRq[5] = pow(sqrt(0.4472135954999581*fIn[5]+0.4472135954999581*fIn[4]-0.9*fIn[3]+0.6708203932499369*fIn[2]-0.6708203932499369*fIn[1]+0.5*fIn[0]),qExp); 
  sqrtfRq[6] = pow(sqrt((-0.5590169943749475*fIn[5])+0.4472135954999581*fIn[4]+0.6708203932499369*fIn[1]+0.5*fIn[0]),qExp); 
  sqrtfRq[7] = pow(sqrt(0.4472135954999581*fIn[5]+0.4472135954999581*fIn[4]-0.9*fIn[3]-0.6708203932499369*fIn[2]+0.6708203932499369*fIn[1]+0.5*fIn[0]),qExp); 
  sqrtfRq[8] = pow(sqrt(0.4472135954999581*fIn[5]+0.4472135954999581*fIn[4]+0.9*fIn[3]+0.6708203932499369*fIn[2]+0.6708203932499369*fIn[1]+0.5*fIn[0]),qExp); 

  out[0] = 0.154320987654321*(sqrtfRq[8]+sqrtfRq[7])+0.2469135802469136*sqrtfRq[6]+0.154320987654321*(sqrtfRq[5]+sqrtfRq[4])+0.2469135802469136*(sqrtfRq[3]+sqrtfRq[2]+sqrtfRq[1])+0.3950617283950617*sqrtfRq[0]; 
  out[1] = 0.2070433312499806*(sqrtfRq[8]+sqrtfRq[7])+0.3312693299999688*sqrtfRq[6]-0.2070433312499806*(sqrtfRq[5]+sqrtfRq[4])-0.3312693299999688*sqrtfRq[3]; 
  out[2] = 0.2070433312499806*sqrtfRq[8]-0.2070433312499806*sqrtfRq[7]+0.2070433312499806*sqrtfRq[5]-0.2070433312499806*sqrtfRq[4]+0.3312693299999688*sqrtfRq[2]-0.3312693299999688*sqrtfRq[1]; 
  out[3] = 0.2777777777777778*sqrtfRq[8]-0.2777777777777778*(sqrtfRq[7]+sqrtfRq[5])+0.2777777777777778*sqrtfRq[4]; 
  out[4] = 0.1380288874999871*(sqrtfRq[8]+sqrtfRq[7])+0.2208462199999793*sqrtfRq[6]+0.1380288874999871*(sqrtfRq[5]+sqrtfRq[4])+0.2208462199999793*sqrtfRq[3]-0.276057774999974*(sqrtfRq[2]+sqrtfRq[1])-0.4416924399999584*sqrtfRq[0]; 
  out[5] = 0.1380288874999871*(sqrtfRq[8]+sqrtfRq[7])-0.276057774999974*sqrtfRq[6]+0.1380288874999871*(sqrtfRq[5]+sqrtfRq[4])-0.276057774999974*sqrtfRq[3]+0.2208462199999793*(sqrtfRq[2]+sqrtfRq[1])-0.4416924399999584*sqrtfRq[0]; 
}
