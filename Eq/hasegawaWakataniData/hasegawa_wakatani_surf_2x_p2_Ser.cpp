#include <hasegawa_wakatani_mod_decl.h>

double hasegawa_wakatani_surf_2x_p2_Ser_x(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR) 
{ 
  // C_: adiabaticity parameter (T_e*kpar^2/(e^2*n_0*eta*omega_ci)).
  // kappa_: normalized density gradient (rho_s/L_n).
  // xc[2]: cell-center coordinates.
  // dx[2]: cell spacing.
  // phi[8]: electrostatic potential.
  // fIn[16]: input fields (vorticity and density).
  // out[16]: output increment (dy/dt).

  double rdx2 = 2.0/dx[0]; 
  double rdy2 = 2.0/dx[1]; 

  // Surface-averaged phase velocity in this direction.
  double alpha0 = 0.5*(3.872983346207417*phi[6]-3.0*phi[3]+1.732050807568877*phi[2])*rdy2; 

  double alpha[3]; 
  alpha[0] = 0.7071067811865475*(3.872983346207417*phi[6]-3.0*phi[3]+1.732050807568877*phi[2])*rdy2; 
  alpha[1] = -0.7071067811865475*(6.708203932499369*phi[7]-3.872983346207417*phi[5])*rdy2; 

  double alphaQuad;
  double fUpwindQuad[6];
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fUpwindQuad[0] = 0.5*(((-0.9682458365518543*(fInR[7]+fInL[7]))+0.5590169943749475*fInR[5]-0.5590169943749475*fInL[5]-1.118033988749895*fInR[4]+1.118033988749895*fInL[4]+0.8660254037844386*(fInR[1]+fInL[1])-0.5*fInR[0]+0.5*fInL[0])*sgn(alphaQuad)+0.9682458365518543*fInR[7]-0.9682458365518543*fInL[7]-0.5590169943749475*(fInR[5]+fInL[5])+1.118033988749895*(fInR[4]+fInL[4])-0.8660254037844386*fInR[1]+0.8660254037844386*fInL[1]+0.5*(fInR[0]+fInL[0])); 
  fUpwindQuad[3] = 0.5*((-2.0*(xc[0]-0.5*dx[0])*kappa_)+((-0.9682458365518543*(fInR[15]+fInL[15]))+0.5590169943749475*fInR[13]-0.5590169943749475*fInL[13]-1.118033988749895*fInR[12]+1.118033988749895*fInL[12]+0.8660254037844386*(fInR[9]+fInL[9])-0.5*fInR[8]+0.5*fInL[8])*sgn(alphaQuad)+0.9682458365518543*fInR[15]-0.9682458365518543*fInL[15]-0.5590169943749475*(fInR[13]+fInL[13])+1.118033988749895*(fInR[12]+fInL[12])-0.8660254037844386*fInR[9]+0.8660254037844386*fInL[9]+0.5*(fInR[8]+fInL[8])); 
  alphaQuad = 0.7071067811865475*alpha[0]-0.9486832980505137*alpha[1]; 
  fUpwindQuad[1] = 0.5*((0.7745966692414834*(fInR[7]+fInL[7])+1.5*fInR[6]-1.5*fInL[6]-0.447213595499958*fInR[5]+0.447213595499958*fInL[5]-1.118033988749895*fInR[4]+1.118033988749895*fInL[4]-1.161895003862225*(fInR[3]+fInL[3])+0.6708203932499369*fInR[2]-0.6708203932499369*fInL[2]+0.8660254037844386*(fInR[1]+fInL[1])-0.5*fInR[0]+0.5*fInL[0])*sgn(alphaQuad)-0.7745966692414834*fInR[7]+0.7745966692414834*fInL[7]-1.5*(fInR[6]+fInL[6])+0.447213595499958*(fInR[5]+fInL[5])+1.118033988749895*(fInR[4]+fInL[4])+1.161895003862225*fInR[3]-1.161895003862225*fInL[3]-0.6708203932499369*(fInR[2]+fInL[2])-0.8660254037844386*fInR[1]+0.8660254037844386*fInL[1]+0.5*(fInR[0]+fInL[0])); 
  fUpwindQuad[4] = 0.5*((-2.0*(xc[0]-0.5*dx[0])*kappa_)+(0.7745966692414834*(fInR[15]+fInL[15])+1.5*fInR[14]-1.5*fInL[14]-0.447213595499958*fInR[13]+0.447213595499958*fInL[13]-1.118033988749895*fInR[12]+1.118033988749895*fInL[12]-1.161895003862225*(fInR[11]+fInL[11])+0.6708203932499369*fInR[10]-0.6708203932499369*fInL[10]+0.8660254037844386*(fInR[9]+fInL[9])-0.5*fInR[8]+0.5*fInL[8])*sgn(alphaQuad)-0.7745966692414834*fInR[15]+0.7745966692414834*fInL[15]-1.5*(fInR[14]+fInL[14])+0.447213595499958*(fInR[13]+fInL[13])+1.118033988749895*(fInR[12]+fInL[12])+1.161895003862225*fInR[11]-1.161895003862225*fInL[11]-0.6708203932499369*(fInR[10]+fInL[10])-0.8660254037844386*fInR[9]+0.8660254037844386*fInL[9]+0.5*(fInR[8]+fInL[8])); 
  alphaQuad = 0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0]; 
  fUpwindQuad[2] = 0.5*((0.7745966692414834*(fInR[7]+fInL[7])-1.5*fInR[6]+1.5*fInL[6]-0.447213595499958*fInR[5]+0.447213595499958*fInL[5]-1.118033988749895*fInR[4]+1.118033988749895*fInL[4]+1.161895003862225*(fInR[3]+fInL[3])-0.6708203932499369*fInR[2]+0.6708203932499369*fInL[2]+0.8660254037844386*(fInR[1]+fInL[1])-0.5*fInR[0]+0.5*fInL[0])*sgn(alphaQuad)-0.7745966692414834*fInR[7]+0.7745966692414834*fInL[7]+1.5*(fInR[6]+fInL[6])+0.447213595499958*(fInR[5]+fInL[5])+1.118033988749895*(fInR[4]+fInL[4])-1.161895003862225*fInR[3]+1.161895003862225*fInL[3]+0.6708203932499369*(fInR[2]+fInL[2])-0.8660254037844386*fInR[1]+0.8660254037844386*fInL[1]+0.5*(fInR[0]+fInL[0])); 
  fUpwindQuad[5] = 0.5*((-2.0*(xc[0]-0.5*dx[0])*kappa_)+(0.7745966692414834*(fInR[15]+fInL[15])-1.5*fInR[14]+1.5*fInL[14]-0.447213595499958*fInR[13]+0.447213595499958*fInL[13]-1.118033988749895*fInR[12]+1.118033988749895*fInL[12]+1.161895003862225*(fInR[11]+fInL[11])-0.6708203932499369*fInR[10]+0.6708203932499369*fInL[10]+0.8660254037844386*(fInR[9]+fInL[9])-0.5*fInR[8]+0.5*fInL[8])*sgn(alphaQuad)-0.7745966692414834*fInR[15]+0.7745966692414834*fInL[15]+1.5*(fInR[14]+fInL[14])+0.447213595499958*(fInR[13]+fInL[13])+1.118033988749895*(fInR[12]+fInL[12])-1.161895003862225*fInR[11]+1.161895003862225*fInL[11]+0.6708203932499369*(fInR[10]+fInL[10])-0.8660254037844386*fInR[9]+0.8660254037844386*fInL[9]+0.5*(fInR[8]+fInL[8])); 

  double fUpwind[6];
  fUpwind[0] = 0.07856742013183861*(5.0*fUpwindQuad[2]+5.0*fUpwindQuad[1]+8.0*fUpwindQuad[0]);
  fUpwind[1] = 0.5270462766947309*(fUpwindQuad[2]-1.0*fUpwindQuad[1]);
  fUpwind[2] = 0.3513641844631533*(fUpwindQuad[2]+fUpwindQuad[1]-2.0*fUpwindQuad[0]);
  fUpwind[3] = 0.07856742013183861*(5.0*fUpwindQuad[5]+5.0*fUpwindQuad[4]+8.0*fUpwindQuad[3]);
  fUpwind[4] = 0.5270462766947309*(fUpwindQuad[5]-1.0*fUpwindQuad[4]);
  fUpwind[5] = 0.3513641844631533*(fUpwindQuad[5]+fUpwindQuad[4]-2.0*fUpwindQuad[3]);

  double incr[16]; 
  incr[0] = 0.5*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0])*rdx2; 
  incr[1] = -0.8660254037844386*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0])*rdx2; 
  incr[2] = 0.1*(4.47213595499958*alpha[1]*fUpwind[2]+5.0*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]))*rdx2; 
  incr[3] = -0.1732050807568877*(4.47213595499958*alpha[1]*fUpwind[2]+5.0*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]))*rdx2; 
  incr[4] = 1.118033988749895*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0])*rdx2; 
  incr[5] = 0.1*(5.0*alpha[0]*fUpwind[2]+4.47213595499958*alpha[1]*fUpwind[1])*rdx2; 
  incr[6] = 0.223606797749979*(4.47213595499958*alpha[1]*fUpwind[2]+5.0*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]))*rdx2; 
  incr[7] = -0.3872983346207417*(2.23606797749979*alpha[0]*fUpwind[2]+2.0*alpha[1]*fUpwind[1])*rdx2; 
  incr[8] = 0.5*(alpha[1]*fUpwind[4]+alpha[0]*fUpwind[3])*rdx2; 
  incr[9] = -0.8660254037844386*(alpha[1]*fUpwind[4]+alpha[0]*fUpwind[3])*rdx2; 
  incr[10] = 0.1*(4.47213595499958*alpha[1]*fUpwind[5]+5.0*(alpha[0]*fUpwind[4]+alpha[1]*fUpwind[3]))*rdx2; 
  incr[11] = -0.1732050807568877*(4.47213595499958*alpha[1]*fUpwind[5]+5.0*(alpha[0]*fUpwind[4]+alpha[1]*fUpwind[3]))*rdx2; 
  incr[12] = 1.118033988749895*(alpha[1]*fUpwind[4]+alpha[0]*fUpwind[3])*rdx2; 
  incr[13] = 0.1*(5.0*alpha[0]*fUpwind[5]+4.47213595499958*alpha[1]*fUpwind[4])*rdx2; 
  incr[14] = 0.223606797749979*(4.47213595499958*alpha[1]*fUpwind[5]+5.0*(alpha[0]*fUpwind[4]+alpha[1]*fUpwind[3]))*rdx2; 
  incr[15] = -0.3872983346207417*(2.23606797749979*alpha[0]*fUpwind[5]+2.0*alpha[1]*fUpwind[4])*rdx2; 

  outR[0] += incr[0]; 
  outR[1] += incr[1]; 
  outR[2] += incr[2]; 
  outR[3] += incr[3]; 
  outR[4] += incr[4]; 
  outR[5] += incr[5]; 
  outR[6] += incr[6]; 
  outR[7] += incr[7]; 
  outR[8] += incr[8]; 
  outR[9] += incr[9]; 
  outR[10] += incr[10]; 
  outR[11] += incr[11]; 
  outR[12] += incr[12]; 
  outR[13] += incr[13]; 
  outR[14] += incr[14]; 
  outR[15] += incr[15]; 

  outL[0] += -1.0*incr[0]; 
  outL[1] += incr[1]; 
  outL[2] += -1.0*incr[2]; 
  outL[3] += incr[3]; 
  outL[4] += -1.0*incr[4]; 
  outL[5] += -1.0*incr[5]; 
  outL[6] += -1.0*incr[6]; 
  outL[7] += incr[7]; 
  outL[8] += -1.0*incr[8]; 
  outL[9] += incr[9]; 
  outL[10] += -1.0*incr[10]; 
  outL[11] += incr[11]; 
  outL[12] += -1.0*incr[12]; 
  outL[13] += -1.0*incr[13]; 
  outL[14] += -1.0*incr[14]; 
  outL[15] += incr[15]; 

  return std::abs(alpha0); 
} 
double hasegawa_wakatani_surf_2x_p2_Ser_y(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR) 
{ 
  // C_: adiabaticity parameter (T_e*kpar^2/(e^2*n_0*eta*omega_ci)).
  // kappa_: normalized density gradient (rho_s/L_n).
  // xc[2]: cell-center coordinates.
  // dx[2]: cell spacing.
  // phi[8]: electrostatic potential.
  // fIn[16]: input fields (vorticity and density).
  // out[16]: output increment (dy/dt).

  double rdx2 = 2.0/dx[0]; 
  double rdy2 = 2.0/dx[1]; 

  // Surface-averaged phase velocity in this direction.
  double alpha0 = -0.5*(3.872983346207417*phi[7]-3.0*phi[3]+1.732050807568877*phi[1])*rdx2; 

  double alpha[3]; 
  alpha[0] = -0.7071067811865475*(3.872983346207417*phi[7]-3.0*phi[3]+1.732050807568877*phi[1])*rdx2; 
  alpha[1] = 0.7071067811865475*(6.708203932499369*phi[6]-3.872983346207417*phi[4])*rdx2; 

  double alphaQuad;
  double fUpwindQuad[6];
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fUpwindQuad[0] = 0.5*(((-0.9682458365518543*(fInR[6]+fInL[6]))-1.118033988749895*fInR[5]+1.118033988749895*fInL[5]+0.5590169943749475*fInR[4]-0.5590169943749475*fInL[4]+0.8660254037844386*(fInR[2]+fInL[2])-0.5*fInR[0]+0.5*fInL[0])*sgn(alphaQuad)+0.9682458365518543*fInR[6]-0.9682458365518543*fInL[6]+1.118033988749895*(fInR[5]+fInL[5])-0.5590169943749475*(fInR[4]+fInL[4])-0.8660254037844386*fInR[2]+0.8660254037844386*fInL[2]+0.5*(fInR[0]+fInL[0])); 
  fUpwindQuad[3] = 0.5*((-2.0*(xc[0]-0.5*dx[0])*kappa_)+((-0.9682458365518543*(fInR[14]+fInL[14]))-1.118033988749895*fInR[13]+1.118033988749895*fInL[13]+0.5590169943749475*fInR[12]-0.5590169943749475*fInL[12]+0.8660254037844386*(fInR[10]+fInL[10])-0.5*fInR[8]+0.5*fInL[8])*sgn(alphaQuad)+0.9682458365518543*fInR[14]-0.9682458365518543*fInL[14]+1.118033988749895*(fInR[13]+fInL[13])-0.5590169943749475*(fInR[12]+fInL[12])-0.8660254037844386*fInR[10]+0.8660254037844386*fInL[10]+0.5*(fInR[8]+fInL[8])); 
  alphaQuad = 0.7071067811865475*alpha[0]-0.9486832980505137*alpha[1]; 
  fUpwindQuad[1] = 0.5*((1.5*fInR[7]-1.5*fInL[7]+0.7745966692414834*(fInR[6]+fInL[6])-1.118033988749895*fInR[5]+1.118033988749895*fInL[5]-0.447213595499958*fInR[4]+0.447213595499958*fInL[4]-1.161895003862225*(fInR[3]+fInL[3])+0.8660254037844386*(fInR[2]+fInL[2])+0.6708203932499369*fInR[1]-0.6708203932499369*fInL[1]-0.5*fInR[0]+0.5*fInL[0])*sgn(alphaQuad)-1.5*(fInR[7]+fInL[7])-0.7745966692414834*fInR[6]+0.7745966692414834*fInL[6]+1.118033988749895*(fInR[5]+fInL[5])+0.447213595499958*(fInR[4]+fInL[4])+1.161895003862225*fInR[3]-1.161895003862225*fInL[3]-0.8660254037844386*fInR[2]+0.8660254037844386*fInL[2]-0.6708203932499369*(fInR[1]+fInL[1])+0.5*(fInR[0]+fInL[0])); 
  fUpwindQuad[4] = 0.5*((-2.0*(xc[0]-0.5*dx[0])*kappa_)+(1.5*fInR[15]-1.5*fInL[15]+0.7745966692414834*(fInR[14]+fInL[14])-1.118033988749895*fInR[13]+1.118033988749895*fInL[13]-0.447213595499958*fInR[12]+0.447213595499958*fInL[12]-1.161895003862225*(fInR[11]+fInL[11])+0.8660254037844386*(fInR[10]+fInL[10])+0.6708203932499369*fInR[9]-0.6708203932499369*fInL[9]-0.5*fInR[8]+0.5*fInL[8])*sgn(alphaQuad)-1.5*(fInR[15]+fInL[15])-0.7745966692414834*fInR[14]+0.7745966692414834*fInL[14]+1.118033988749895*(fInR[13]+fInL[13])+0.447213595499958*(fInR[12]+fInL[12])+1.161895003862225*fInR[11]-1.161895003862225*fInL[11]-0.8660254037844386*fInR[10]+0.8660254037844386*fInL[10]-0.6708203932499369*(fInR[9]+fInL[9])+0.5*(fInR[8]+fInL[8])); 
  alphaQuad = 0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0]; 
  fUpwindQuad[2] = 0.5*(((-1.5*fInR[7])+1.5*fInL[7]+0.7745966692414834*(fInR[6]+fInL[6])-1.118033988749895*fInR[5]+1.118033988749895*fInL[5]-0.447213595499958*fInR[4]+0.447213595499958*fInL[4]+1.161895003862225*(fInR[3]+fInL[3])+0.8660254037844386*(fInR[2]+fInL[2])-0.6708203932499369*fInR[1]+0.6708203932499369*fInL[1]-0.5*fInR[0]+0.5*fInL[0])*sgn(alphaQuad)+1.5*(fInR[7]+fInL[7])-0.7745966692414834*fInR[6]+0.7745966692414834*fInL[6]+1.118033988749895*(fInR[5]+fInL[5])+0.447213595499958*(fInR[4]+fInL[4])-1.161895003862225*fInR[3]+1.161895003862225*fInL[3]-0.8660254037844386*fInR[2]+0.8660254037844386*fInL[2]+0.6708203932499369*(fInR[1]+fInL[1])+0.5*(fInR[0]+fInL[0])); 
  fUpwindQuad[5] = 0.5*((-2.0*(xc[0]-0.5*dx[0])*kappa_)+((-1.5*fInR[15])+1.5*fInL[15]+0.7745966692414834*(fInR[14]+fInL[14])-1.118033988749895*fInR[13]+1.118033988749895*fInL[13]-0.447213595499958*fInR[12]+0.447213595499958*fInL[12]+1.161895003862225*(fInR[11]+fInL[11])+0.8660254037844386*(fInR[10]+fInL[10])-0.6708203932499369*fInR[9]+0.6708203932499369*fInL[9]-0.5*fInR[8]+0.5*fInL[8])*sgn(alphaQuad)+1.5*(fInR[15]+fInL[15])-0.7745966692414834*fInR[14]+0.7745966692414834*fInL[14]+1.118033988749895*(fInR[13]+fInL[13])+0.447213595499958*(fInR[12]+fInL[12])-1.161895003862225*fInR[11]+1.161895003862225*fInL[11]-0.8660254037844386*fInR[10]+0.8660254037844386*fInL[10]+0.6708203932499369*(fInR[9]+fInL[9])+0.5*(fInR[8]+fInL[8])); 

  double fUpwind[6];
  fUpwind[0] = 0.07856742013183861*(5.0*fUpwindQuad[2]+5.0*fUpwindQuad[1]+8.0*fUpwindQuad[0]);
  fUpwind[1] = 0.5270462766947309*(fUpwindQuad[2]-1.0*fUpwindQuad[1]);
  fUpwind[2] = 0.3513641844631533*(fUpwindQuad[2]+fUpwindQuad[1]-2.0*fUpwindQuad[0]);
  fUpwind[3] = 0.07856742013183861*(5.0*fUpwindQuad[5]+5.0*fUpwindQuad[4]+8.0*fUpwindQuad[3]);
  fUpwind[4] = 0.5270462766947309*(fUpwindQuad[5]-1.0*fUpwindQuad[4]);
  fUpwind[5] = 0.3513641844631533*(fUpwindQuad[5]+fUpwindQuad[4]-2.0*fUpwindQuad[3]);

  double incr[16]; 
  incr[0] = 0.5*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0])*rdy2; 
  incr[1] = 0.1*(4.47213595499958*alpha[1]*fUpwind[2]+5.0*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]))*rdy2; 
  incr[2] = -0.8660254037844386*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0])*rdy2; 
  incr[3] = -0.1732050807568877*(4.47213595499958*alpha[1]*fUpwind[2]+5.0*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]))*rdy2; 
  incr[4] = 0.1*(5.0*alpha[0]*fUpwind[2]+4.47213595499958*alpha[1]*fUpwind[1])*rdy2; 
  incr[5] = 1.118033988749895*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0])*rdy2; 
  incr[6] = -0.3872983346207417*(2.23606797749979*alpha[0]*fUpwind[2]+2.0*alpha[1]*fUpwind[1])*rdy2; 
  incr[7] = 0.223606797749979*(4.47213595499958*alpha[1]*fUpwind[2]+5.0*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]))*rdy2; 
  incr[8] = 0.5*(alpha[1]*fUpwind[4]+alpha[0]*fUpwind[3])*rdy2; 
  incr[9] = 0.1*(4.47213595499958*alpha[1]*fUpwind[5]+5.0*(alpha[0]*fUpwind[4]+alpha[1]*fUpwind[3]))*rdy2; 
  incr[10] = -0.8660254037844386*(alpha[1]*fUpwind[4]+alpha[0]*fUpwind[3])*rdy2; 
  incr[11] = -0.1732050807568877*(4.47213595499958*alpha[1]*fUpwind[5]+5.0*(alpha[0]*fUpwind[4]+alpha[1]*fUpwind[3]))*rdy2; 
  incr[12] = 0.1*(5.0*alpha[0]*fUpwind[5]+4.47213595499958*alpha[1]*fUpwind[4])*rdy2; 
  incr[13] = 1.118033988749895*(alpha[1]*fUpwind[4]+alpha[0]*fUpwind[3])*rdy2; 
  incr[14] = -0.3872983346207417*(2.23606797749979*alpha[0]*fUpwind[5]+2.0*alpha[1]*fUpwind[4])*rdy2; 
  incr[15] = 0.223606797749979*(4.47213595499958*alpha[1]*fUpwind[5]+5.0*(alpha[0]*fUpwind[4]+alpha[1]*fUpwind[3]))*rdy2; 

  outR[0] += incr[0]; 
  outR[1] += incr[1]; 
  outR[2] += incr[2]; 
  outR[3] += incr[3]; 
  outR[4] += incr[4]; 
  outR[5] += incr[5]; 
  outR[6] += incr[6]; 
  outR[7] += incr[7]; 
  outR[8] += incr[8]; 
  outR[9] += incr[9]; 
  outR[10] += incr[10]; 
  outR[11] += incr[11]; 
  outR[12] += incr[12]; 
  outR[13] += incr[13]; 
  outR[14] += incr[14]; 
  outR[15] += incr[15]; 

  outL[0] += -1.0*incr[0]; 
  outL[1] += -1.0*incr[1]; 
  outL[2] += incr[2]; 
  outL[3] += incr[3]; 
  outL[4] += -1.0*incr[4]; 
  outL[5] += -1.0*incr[5]; 
  outL[6] += incr[6]; 
  outL[7] += -1.0*incr[7]; 
  outL[8] += -1.0*incr[8]; 
  outL[9] += -1.0*incr[9]; 
  outL[10] += incr[10]; 
  outL[11] += incr[11]; 
  outL[12] += -1.0*incr[12]; 
  outL[13] += -1.0*incr[13]; 
  outL[14] += incr[14]; 
  outL[15] += -1.0*incr[15]; 

  return std::abs(alpha0); 
} 
