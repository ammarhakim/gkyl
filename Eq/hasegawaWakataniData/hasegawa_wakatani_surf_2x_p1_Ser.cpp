#include <hasegawa_wakatani_mod_decl.h>

double hasegawa_wakatani_surf_2x_p1_Ser_x(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR) 
{ 
  // C_: adiabaticity parameter (T_e*kpar^2/(e^2*n_0*eta*omega_ci)).
  // kappa_: normalized density gradient (rho_s/L_n).
  // xc[2]: cell-center coordinates.
  // dx[2]: cell spacing.
  // phi[4]: electrostatic potential.
  // fIn[8]: input fields (vorticity and density).
  // out[8]: output increment (dy/dt).

  double rdx2 = 2.0/dx[0]; 
  double rdy2 = 2.0/dx[1]; 

  // Surface-averaged phase velocity in this direction.
  double alpha0 = -0.5*(3.0*phi[3]-1.732050807568877*phi[2])*rdy2; 

  double alpha[2]; 
  alpha[0] = -0.5*(4.242640687119286*phi[3]-2.449489742783178*phi[2])*rdy2; 

  double alphaQuad;
  double fUpwindQuad[4];
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fUpwindQuad[0] = 0.5*(((-0.8660254037844386*(fInR[3]+fInL[3]))+0.4999999999999999*fInR[2]-0.4999999999999999*fInL[2]+0.8660254037844386*(fInR[1]+fInL[1])-0.5*fInR[0]+0.5*fInL[0])*sgn(alphaQuad)+0.8660254037844386*fInR[3]-0.8660254037844386*fInL[3]-0.4999999999999999*(fInR[2]+fInL[2])-0.8660254037844386*fInR[1]+0.8660254037844386*fInL[1]+0.5*(fInR[0]+fInL[0])); 
  fUpwindQuad[2] = 0.5*((-2.0*(xc[0]-0.5*dx[0])*kappa_)+((-0.8660254037844386*(fInR[7]+fInL[7]))+0.4999999999999999*fInR[6]-0.4999999999999999*fInL[6]+0.8660254037844386*(fInR[5]+fInL[5])-0.5*fInR[4]+0.5*fInL[4])*sgn(alphaQuad)+0.8660254037844386*fInR[7]-0.8660254037844386*fInL[7]-0.4999999999999999*(fInR[6]+fInL[6])-0.8660254037844386*fInR[5]+0.8660254037844386*fInL[5]+0.5*(fInR[4]+fInL[4])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fUpwindQuad[1] = 0.5*((0.8660254037844386*(fInR[3]+fInL[3])-0.4999999999999999*fInR[2]+0.4999999999999999*fInL[2]+0.8660254037844386*(fInR[1]+fInL[1])-0.5*fInR[0]+0.5*fInL[0])*sgn(alphaQuad)-0.8660254037844386*fInR[3]+0.8660254037844386*fInL[3]+0.4999999999999999*(fInR[2]+fInL[2])-0.8660254037844386*fInR[1]+0.8660254037844386*fInL[1]+0.5*(fInR[0]+fInL[0])); 
  fUpwindQuad[3] = 0.5*((-2.0*(xc[0]-0.5*dx[0])*kappa_)+(0.8660254037844386*(fInR[7]+fInL[7])-0.4999999999999999*fInR[6]+0.4999999999999999*fInL[6]+0.8660254037844386*(fInR[5]+fInL[5])-0.5*fInR[4]+0.5*fInL[4])*sgn(alphaQuad)-0.8660254037844386*fInR[7]+0.8660254037844386*fInL[7]+0.4999999999999999*(fInR[6]+fInL[6])-0.8660254037844386*fInR[5]+0.8660254037844386*fInL[5]+0.5*(fInR[4]+fInL[4])); 

  double fUpwind[4];
  fUpwind[0] = 0.7071067811865475*(fUpwindQuad[1]+fUpwindQuad[0]);
  fUpwind[1] = 0.7071067811865475*(fUpwindQuad[1]-1.0*fUpwindQuad[0]);
  fUpwind[2] = 0.7071067811865475*(fUpwindQuad[3]+fUpwindQuad[2]);
  fUpwind[3] = 0.7071067811865475*(fUpwindQuad[3]-1.0*fUpwindQuad[2]);

  double incr[8]; 
  incr[0] = 0.5*alpha[0]*fUpwind[0]*rdx2; 
  incr[1] = -0.8660254037844386*alpha[0]*fUpwind[0]*rdx2; 
  incr[2] = 0.5*alpha[0]*fUpwind[1]*rdx2; 
  incr[3] = -0.8660254037844386*alpha[0]*fUpwind[1]*rdx2; 
  incr[4] = 0.5*alpha[0]*fUpwind[2]*rdx2; 
  incr[5] = -0.8660254037844386*alpha[0]*fUpwind[2]*rdx2; 
  incr[6] = 0.5*alpha[0]*fUpwind[3]*rdx2; 
  incr[7] = -0.8660254037844386*alpha[0]*fUpwind[3]*rdx2; 

  outR[0] += incr[0]; 
  outR[1] += incr[1]; 
  outR[2] += incr[2]; 
  outR[3] += incr[3]; 
  outR[4] += incr[4]; 
  outR[5] += incr[5]; 
  outR[6] += incr[6]; 
  outR[7] += incr[7]; 

  outL[0] += -1.0*incr[0]; 
  outL[1] += incr[1]; 
  outL[2] += -1.0*incr[2]; 
  outL[3] += incr[3]; 
  outL[4] += -1.0*incr[4]; 
  outL[5] += incr[5]; 
  outL[6] += -1.0*incr[6]; 
  outL[7] += incr[7]; 

  return std::abs(alpha0); 
} 
double hasegawa_wakatani_surf_2x_p1_Ser_y(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR) 
{ 
  // C_: adiabaticity parameter (T_e*kpar^2/(e^2*n_0*eta*omega_ci)).
  // kappa_: normalized density gradient (rho_s/L_n).
  // xc[2]: cell-center coordinates.
  // dx[2]: cell spacing.
  // phi[4]: electrostatic potential.
  // fIn[8]: input fields (vorticity and density).
  // out[8]: output increment (dy/dt).

  double rdx2 = 2.0/dx[0]; 
  double rdy2 = 2.0/dx[1]; 

  // Surface-averaged phase velocity in this direction.
  double alpha0 = 0.5*(3.0*phi[3]-1.732050807568877*phi[1])*rdx2; 

  double alpha[2]; 
  alpha[0] = 0.5*(4.242640687119286*phi[3]-2.449489742783178*phi[1])*rdx2; 

  double alphaQuad;
  double fUpwindQuad[4];
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fUpwindQuad[0] = 0.5*(((-0.8660254037844386*(fInR[3]+fInL[3]))+0.8660254037844386*(fInR[2]+fInL[2])+0.4999999999999999*fInR[1]-0.4999999999999999*fInL[1]-0.5*fInR[0]+0.5*fInL[0])*sgn(alphaQuad)+0.8660254037844386*fInR[3]-0.8660254037844386*(fInL[3]+fInR[2])+0.8660254037844386*fInL[2]-0.4999999999999999*(fInR[1]+fInL[1])+0.5*(fInR[0]+fInL[0])); 
  fUpwindQuad[2] = 0.5*((-2.0*(xc[0]-0.5*dx[0])*kappa_)+((-0.8660254037844386*(fInR[7]+fInL[7]))+0.8660254037844386*(fInR[6]+fInL[6])+0.4999999999999999*fInR[5]-0.4999999999999999*fInL[5]-0.5*fInR[4]+0.5*fInL[4])*sgn(alphaQuad)+0.8660254037844386*fInR[7]-0.8660254037844386*(fInL[7]+fInR[6])+0.8660254037844386*fInL[6]-0.4999999999999999*(fInR[5]+fInL[5])+0.5*(fInR[4]+fInL[4])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fUpwindQuad[1] = 0.5*((0.8660254037844386*(fInR[3]+fInL[3]+fInR[2]+fInL[2])-0.4999999999999999*fInR[1]+0.4999999999999999*fInL[1]-0.5*fInR[0]+0.5*fInL[0])*sgn(alphaQuad)-0.8660254037844386*(fInR[3]+fInR[2])+0.8660254037844386*(fInL[3]+fInL[2])+0.4999999999999999*(fInR[1]+fInL[1])+0.5*(fInR[0]+fInL[0])); 
  fUpwindQuad[3] = 0.5*((-2.0*(xc[0]-0.5*dx[0])*kappa_)+(0.8660254037844386*(fInR[7]+fInL[7]+fInR[6]+fInL[6])-0.4999999999999999*fInR[5]+0.4999999999999999*fInL[5]-0.5*fInR[4]+0.5*fInL[4])*sgn(alphaQuad)-0.8660254037844386*(fInR[7]+fInR[6])+0.8660254037844386*(fInL[7]+fInL[6])+0.4999999999999999*(fInR[5]+fInL[5])+0.5*(fInR[4]+fInL[4])); 

  double fUpwind[4];
  fUpwind[0] = 0.7071067811865475*(fUpwindQuad[1]+fUpwindQuad[0]);
  fUpwind[1] = 0.7071067811865475*(fUpwindQuad[1]-1.0*fUpwindQuad[0]);
  fUpwind[2] = 0.7071067811865475*(fUpwindQuad[3]+fUpwindQuad[2]);
  fUpwind[3] = 0.7071067811865475*(fUpwindQuad[3]-1.0*fUpwindQuad[2]);

  double incr[8]; 
  incr[0] = 0.5*alpha[0]*fUpwind[0]*rdy2; 
  incr[1] = 0.5*alpha[0]*fUpwind[1]*rdy2; 
  incr[2] = -0.8660254037844386*alpha[0]*fUpwind[0]*rdy2; 
  incr[3] = -0.8660254037844386*alpha[0]*fUpwind[1]*rdy2; 
  incr[4] = 0.5*alpha[0]*fUpwind[2]*rdy2; 
  incr[5] = 0.5*alpha[0]*fUpwind[3]*rdy2; 
  incr[6] = -0.8660254037844386*alpha[0]*fUpwind[2]*rdy2; 
  incr[7] = -0.8660254037844386*alpha[0]*fUpwind[3]*rdy2; 

  outR[0] += incr[0]; 
  outR[1] += incr[1]; 
  outR[2] += incr[2]; 
  outR[3] += incr[3]; 
  outR[4] += incr[4]; 
  outR[5] += incr[5]; 
  outR[6] += incr[6]; 
  outR[7] += incr[7]; 

  outL[0] += -1.0*incr[0]; 
  outL[1] += -1.0*incr[1]; 
  outL[2] += incr[2]; 
  outL[3] += incr[3]; 
  outL[4] += -1.0*incr[4]; 
  outL[5] += -1.0*incr[5]; 
  outL[6] += incr[6]; 
  outL[7] += incr[7]; 

  return std::abs(alpha0); 
} 
