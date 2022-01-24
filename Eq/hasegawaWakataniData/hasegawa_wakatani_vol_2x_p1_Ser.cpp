#include <hasegawa_wakatani_mod_decl.h>

double hasegawa_wakatani_vol_2x_p1_Ser(const double C_, const double kappa_, const double *xc, const double *dx, const double *phi, const double *fIn, double *out) 
{ 
  // C_: adiabaticity parameter (T_e*kpar^2/(e^2*n_0*eta*omega_ci)).
  // kappa_: normalized density gradient (rho_s/L_n).
  // xc[2]: cell-center coordinates.
  // dx[2]: cell spacing.
  // phi[4]: electrostatic potential.
  // fIn[8]: input fields (vorticity and density).
  // out[8]: output increment (dy/dt).

  double rDx2 = 2.0/dx[0]; 
  double rDy2 = 2.0/dx[1]; 
  double cflFreq = 0.0;
  double alphaL  = 0.0;
  double alphaR  = 0.0;
  double alpha_x[4]; 
  alpha_x[0] = 1.732050807568877*phi[2]*rDx2*rDy2; 
  alpha_x[1] = 1.732050807568877*phi[3]*rDx2*rDy2; 
  // Evaluate alpha at left surface quadrature points. 
  alphaL = 0.5*alpha_x[0]-0.8660254037844386*alpha_x[1]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.5*alpha_x[0]-0.8660254037844386*alpha_x[1]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points. 
  alphaR = 0.8660254037844386*alpha_x[1]+0.5*alpha_x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.8660254037844386*alpha_x[1]+0.5*alpha_x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha_y[4]; 
  alpha_y[0] = -1.732050807568877*phi[1]*rDx2*rDy2; 
  alpha_y[2] = -1.732050807568877*phi[3]*rDx2*rDy2; 
  // Evaluate alpha at left surface quadrature points. 
  alphaL = 0.5*alpha_y[0]-0.8660254037844386*alpha_y[2]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.5*alpha_y[0]-0.8660254037844386*alpha_y[2]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points. 
  alphaR = 0.8660254037844386*alpha_y[2]+0.5*alpha_y[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.8660254037844386*alpha_y[2]+0.5*alpha_y[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[0] += -1.0*(fIn[4]-1.0*phi[0])*C_; 
  out[1] += -0.5*((2.0*fIn[5]-2.0*phi[1])*C_-1.732050807568877*(alpha_x[1]*fIn[1]+alpha_x[0]*fIn[0])); 
  out[2] += -0.5*((2.0*fIn[6]-2.0*phi[2])*C_-1.732050807568877*(alpha_y[2]*fIn[2]+alpha_y[0]*fIn[0])); 
  out[3] += -0.5*((2.0*fIn[7]-2.0*phi[3])*C_-1.732050807568877*((alpha_y[2]+alpha_x[1])*fIn[3]+alpha_x[0]*fIn[2]+alpha_y[0]*fIn[1])); 

  out[4] += -1.0*(fIn[4]-1.0*phi[0])*C_; 
  out[5] += -0.5*(2.0*alpha_x[1]*kappa_+(2.0*fIn[5]-2.0*phi[1])*C_-1.732050807568877*(alpha_x[1]*fIn[5]+alpha_x[0]*fIn[4])); 
  out[6] += -0.5*((2.0*fIn[6]-2.0*phi[2])*C_-1.732050807568877*(alpha_y[2]*fIn[6]+alpha_y[0]*fIn[4])); 
  out[7] += -0.5*(2.0*alpha_y[0]*kappa_+(2.0*fIn[7]-2.0*phi[3])*C_-1.732050807568877*((alpha_y[2]+alpha_x[1])*fIn[7]+alpha_x[0]*fIn[6]+alpha_y[0]*fIn[5])); 

  return cflFreq; 
} 
