#include <hasegawa_wakatani_mod_decl.h>

double hasegawa_wakatani_vol_2x_p2_Ser(const double C_, const double kappa_, const double *xc, const double *dx, const double *phi, const double *fIn, double *out) 
{ 
  // C_: adiabaticity parameter (T_e*kpar^2/(e^2*n_0*eta*omega_ci)).
  // kappa_: normalized density gradient (rho_s/L_n).
  // xc[2]: cell-center coordinates.
  // dx[2]: cell spacing.
  // phi[8]: electrostatic potential.
  // fIn[16]: input fields (vorticity and density).
  // out[16]: output increment (dy/dt).

  double rDx2 = 2.0/dx[0]; 
  double rDy2 = 2.0/dx[1]; 
  double cflFreq = 0.0;
  double alphaL  = 0.0;
  double alphaR  = 0.0;
  double alpha_x[8]; 
  alpha_x[0] = 1.732050807568877*phi[2]*rDx2*rDy2; 
  alpha_x[1] = 1.732050807568877*phi[3]*rDx2*rDy2; 
  alpha_x[2] = 3.872983346207417*phi[5]*rDx2*rDy2; 
  alpha_x[3] = 3.872983346207417*phi[7]*rDx2*rDy2; 
  alpha_x[4] = 1.732050807568877*phi[6]*rDx2*rDy2; 
  // Evaluate alpha at left surface quadrature points. 
  alphaL = 1.118033988749895*alpha_x[4]+1.161895003862225*alpha_x[3]-0.6708203932499369*alpha_x[2]-0.8660254037844386*alpha_x[1]+0.5*alpha_x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 1.118033988749895*alpha_x[4]-0.8660254037844386*alpha_x[1]+0.5*alpha_x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 1.118033988749895*alpha_x[4]-1.161895003862225*alpha_x[3]+0.6708203932499369*alpha_x[2]-0.8660254037844386*alpha_x[1]+0.5*alpha_x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points. 
  alphaR = 1.118033988749895*alpha_x[4]-1.161895003862225*alpha_x[3]-0.6708203932499369*alpha_x[2]+0.8660254037844386*alpha_x[1]+0.5*alpha_x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 1.118033988749895*alpha_x[4]+0.8660254037844386*alpha_x[1]+0.5*alpha_x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 1.118033988749895*alpha_x[4]+1.161895003862225*alpha_x[3]+0.6708203932499369*alpha_x[2]+0.8660254037844386*alpha_x[1]+0.5*alpha_x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha_y[8]; 
  alpha_y[0] = -1.732050807568877*phi[1]*rDx2*rDy2; 
  alpha_y[1] = -3.872983346207417*phi[4]*rDx2*rDy2; 
  alpha_y[2] = -1.732050807568877*phi[3]*rDx2*rDy2; 
  alpha_y[3] = -3.872983346207417*phi[6]*rDx2*rDy2; 
  alpha_y[5] = -1.732050807568877*phi[7]*rDx2*rDy2; 
  // Evaluate alpha at left surface quadrature points. 
  alphaL = 1.118033988749895*alpha_y[5]+1.161895003862225*alpha_y[3]-0.8660254037844386*alpha_y[2]-0.6708203932499369*alpha_y[1]+0.5*alpha_y[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 1.118033988749895*alpha_y[5]-0.8660254037844386*alpha_y[2]+0.5*alpha_y[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 1.118033988749895*alpha_y[5]-1.161895003862225*alpha_y[3]-0.8660254037844386*alpha_y[2]+0.6708203932499369*alpha_y[1]+0.5*alpha_y[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points. 
  alphaR = 1.118033988749895*alpha_y[5]-1.161895003862225*alpha_y[3]+0.8660254037844386*alpha_y[2]-0.6708203932499369*alpha_y[1]+0.5*alpha_y[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 1.118033988749895*alpha_y[5]+0.8660254037844386*alpha_y[2]+0.5*alpha_y[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 1.118033988749895*alpha_y[5]+1.161895003862225*alpha_y[3]+0.8660254037844386*alpha_y[2]+0.6708203932499369*alpha_y[1]+0.5*alpha_y[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[0] += -1.0*(fIn[8]-1.0*phi[0])*C_; 
  out[1] += -0.5*((2.0*fIn[9]-2.0*phi[1])*C_-1.732050807568877*(alpha_x[4]*fIn[4]+alpha_x[3]*fIn[3]+alpha_x[2]*fIn[2]+alpha_x[1]*fIn[1]+alpha_x[0]*fIn[0])); 
  out[2] += -0.5*((2.0*fIn[10]-2.0*phi[2])*C_-1.732050807568877*(alpha_y[5]*fIn[5]+alpha_y[3]*fIn[3]+alpha_y[2]*fIn[2]+alpha_y[1]*fIn[1]+alpha_y[0]*fIn[0])); 
  out[3] += -0.1*((10.0*fIn[11]-10.0*phi[3])*C_+((-8.660254037844387*alpha_y[5])-7.745966692414834*alpha_x[3])*fIn[7]-8.660254037844387*alpha_x[4]*fIn[6]-7.745966692414834*(alpha_y[3]*fIn[6]+alpha_x[2]*fIn[5]+alpha_y[1]*fIn[4])-8.660254037844386*((alpha_y[2]+alpha_x[1])*fIn[3]+fIn[2]*alpha_y[3]+fIn[1]*alpha_x[3]+alpha_x[0]*fIn[2]+fIn[0]*alpha_x[2]+alpha_y[0]*fIn[1]+fIn[0]*alpha_y[1])); 
  out[4] += -0.1*((10.0*fIn[12]-10.0*phi[4])*C_-17.32050807568877*alpha_x[3]*fIn[6]-17.32050807568877*(alpha_x[1]*fIn[4]+fIn[1]*alpha_x[4])-19.36491673103709*(alpha_x[2]*fIn[3]+fIn[2]*alpha_x[3]+alpha_x[0]*fIn[1]+fIn[0]*alpha_x[1])); 
  out[5] += -0.1*((10.0*fIn[13]-10.0*phi[5])*C_-17.32050807568877*alpha_y[3]*fIn[7]-17.32050807568877*(alpha_y[2]*fIn[5]+fIn[2]*alpha_y[5])-19.36491673103709*(alpha_y[1]*fIn[3]+fIn[1]*alpha_y[3]+alpha_y[0]*fIn[2]+fIn[0]*alpha_y[2])); 
  out[6] += -0.1*((10.0*fIn[14]-10.0*phi[6])*C_-17.32050807568877*alpha_x[2]*fIn[7]+((-8.660254037844386*alpha_y[2])-17.32050807568877*alpha_x[1])*fIn[6]-17.32050807568877*alpha_x[3]*fIn[5]+((-17.32050807568877*alpha_x[3])-8.660254037844387*alpha_y[0])*fIn[4]+fIn[3]*((-17.32050807568877*alpha_x[4])-7.745966692414834*alpha_y[3])-19.36491673103708*(alpha_x[0]*fIn[3]+fIn[0]*alpha_x[3]+alpha_x[1]*fIn[2])+fIn[1]*((-19.36491673103708*alpha_x[2])-7.745966692414834*alpha_y[1])); 
  out[7] += -0.1*((10.0*fIn[15]-10.0*phi[7])*C_+((-17.32050807568877*alpha_y[2])-8.660254037844386*alpha_x[1])*fIn[7]-17.32050807568877*alpha_y[1]*fIn[6]+((-17.32050807568877*alpha_y[3])-8.660254037844387*alpha_x[0])*fIn[5]-17.32050807568877*(fIn[3]*alpha_y[5]+alpha_y[3]*fIn[4])-7.745966692414834*alpha_x[3]*fIn[3]-19.36491673103708*(alpha_y[0]*fIn[3]+fIn[0]*alpha_y[3])-7.745966692414834*alpha_x[2]*fIn[2]-19.36491673103708*(alpha_y[1]*fIn[2]+fIn[1]*alpha_y[2])); 

  out[8] += -1.0*(fIn[8]-1.0*phi[0])*C_; 
  out[9] += -0.5*((dx[0]*alpha_x[1]+3.464101615137754*alpha_x[0]*xc[0])*kappa_+(2.0*fIn[9]-2.0*phi[1])*C_-1.732050807568877*(alpha_x[4]*fIn[12]+alpha_x[3]*fIn[11]+alpha_x[2]*fIn[10]+alpha_x[1]*fIn[9]+alpha_x[0]*fIn[8])); 
  out[10] += -0.5*((dx[0]*alpha_y[1]+3.464101615137754*alpha_y[0]*xc[0])*kappa_+(2.0*fIn[10]-2.0*phi[2])*C_-1.732050807568877*(alpha_y[5]*fIn[13]+alpha_y[3]*fIn[11]+alpha_y[2]*fIn[10]+alpha_y[1]*fIn[9]+alpha_y[0]*fIn[8])); 
  out[11] += -0.1*((5.0*dx[0]*alpha_x[3]+17.32050807568877*xc[0]*(alpha_x[2]+alpha_y[1])+5.0*alpha_y[0]*dx[0])*kappa_+(10.0*fIn[11]-10.0*phi[3])*C_+((-8.660254037844387*alpha_y[5])-7.745966692414834*alpha_x[3])*fIn[15]-8.660254037844387*alpha_x[4]*fIn[14]-7.745966692414834*(alpha_y[3]*fIn[14]+alpha_x[2]*fIn[13]+alpha_y[1]*fIn[12])-8.660254037844386*((alpha_y[2]+alpha_x[1])*fIn[11]+(alpha_y[3]+alpha_x[0])*fIn[10]+(alpha_x[3]+alpha_y[0])*fIn[9]+(alpha_x[2]+alpha_y[1])*fIn[8])); 
  out[12] += -0.1*((10.0*dx[0]*alpha_x[4]+38.72983346207418*xc[0]*alpha_x[1]+11.18033988749895*alpha_x[0]*dx[0])*kappa_+(10.0*fIn[12]-10.0*phi[4])*C_-17.32050807568877*alpha_x[3]*fIn[14]-17.32050807568877*alpha_x[1]*fIn[12]-19.36491673103709*(alpha_x[2]*fIn[11]+alpha_x[3]*fIn[10])-17.32050807568877*alpha_x[4]*fIn[9]-19.36491673103709*(alpha_x[0]*fIn[9]+alpha_x[1]*fIn[8])); 
  out[13] += -0.1*((11.18033988749895*dx[0]*alpha_y[3]+38.72983346207418*xc[0]*alpha_y[2])*kappa_+(10.0*fIn[13]-10.0*phi[5])*C_-17.32050807568877*alpha_y[3]*fIn[15]-17.32050807568877*alpha_y[2]*fIn[13]-19.36491673103709*alpha_y[1]*fIn[11]-17.32050807568877*alpha_y[5]*fIn[10]-19.36491673103709*(alpha_y[0]*fIn[10]+alpha_y[3]*fIn[9]+alpha_y[2]*fIn[8])); 
  out[14] += -0.03333333333333333*((116.1895003862225*xc[0]*alpha_x[3]+dx[0]*(33.54101966249684*alpha_x[2]+13.41640786499874*alpha_y[1]))*kappa_+(30.0*fIn[14]-30.0*phi[6])*C_-51.96152422706631*alpha_x[2]*fIn[15]+((-25.98076211353316*alpha_y[2])-51.96152422706631*alpha_x[1])*fIn[14]-51.96152422706632*alpha_x[3]*fIn[13]+((-51.96152422706632*alpha_x[3])-25.98076211353316*alpha_y[0])*fIn[12]+((-51.96152422706632*alpha_x[4])-23.2379000772445*alpha_y[3])*fIn[11]-58.09475019311126*(alpha_x[0]*fIn[11]+alpha_x[1]*fIn[10])+((-58.09475019311126*alpha_x[2])-23.2379000772445*alpha_y[1])*fIn[9]-58.09475019311126*alpha_x[3]*fIn[8]); 
  out[15] += -0.03333333333333333*((116.1895003862225*xc[0]*alpha_y[3]+33.54101966249684*dx[0]*alpha_y[2])*kappa_+(30.0*fIn[15]-30.0*phi[7])*C_+((-51.96152422706631*alpha_y[2])-25.98076211353316*alpha_x[1])*fIn[15]-51.96152422706631*alpha_y[1]*fIn[14]+((-51.96152422706632*alpha_y[3])-25.98076211353316*alpha_x[0])*fIn[13]-51.96152422706632*alpha_y[3]*fIn[12]+((-51.96152422706632*alpha_y[5])-23.2379000772445*alpha_x[3]-58.09475019311126*alpha_y[0])*fIn[11]-23.2379000772445*alpha_x[2]*fIn[10]-58.09475019311126*(alpha_y[1]*fIn[10]+alpha_y[2]*fIn[9]+alpha_y[3]*fIn[8])); 

  return cflFreq; 
} 
