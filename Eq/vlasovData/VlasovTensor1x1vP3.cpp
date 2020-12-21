#include <VlasovModDecl.h> 
__host__ __device__ double VlasovVol1x1vTensorP3(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // EM:        Input EM-field.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[16]; 
  double alpha_vdim[16]; 

  alpha_cdim[0] = 4.0*w0dx0; 
  alpha_cdim[2] = 1.154700538379252*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = 1.414213562373095*E0[0]*dv10; 
  alpha_vdim[1] = 1.414213562373095*E0[1]*dv10; 
  alpha_vdim[4] = 1.414213562373095*E0[2]*dv10; 
  alpha_vdim[8] = 1.414213562373095*E0[3]*dv10; 
  alpha_mid += std::abs(0.25*alpha_vdim[0]-0.2795084971874737*alpha_vdim[4]); 

  out[1] += 0.8660254037844386*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.8660254037844386*(alpha_vdim[8]*f[8]+alpha_vdim[4]*f[4]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.7606388292556648*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8])+0.7745966692414833*(alpha_cdim[2]*f[5]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+0.8660254037844386*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[4] += 1.936491673103709*(alpha_cdim[2]*f[3]+alpha_cdim[0]*f[1]); 
  out[5] += 1.936491673103709*(alpha_vdim[8]*f[11]+alpha_vdim[4]*f[6]+alpha_vdim[1]*f[3]+alpha_vdim[0]*f[2]); 
  out[6] += 0.5163977794943223*alpha_vdim[8]*f[8]+0.7606388292556648*(alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8])+1.732050807568877*alpha_cdim[2]*f[7]+0.5532833351724881*alpha_vdim[4]*f[4]+0.8660254037844386*(alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4])+1.936491673103709*alpha_cdim[0]*f[3]+f[1]*(1.936491673103709*alpha_cdim[2]+0.7745966692414833*alpha_vdim[1]); 
  out[7] += 1.700840128541522*alpha_vdim[4]*f[11]+0.7606388292556648*alpha_cdim[2]*f[9]+f[6]*(1.700840128541522*alpha_vdim[8]+1.732050807568877*alpha_vdim[1])+0.8660254037844386*alpha_cdim[0]*f[5]+f[3]*(1.732050807568877*alpha_vdim[4]+1.936491673103709*alpha_vdim[0])+(0.7745966692414833*alpha_cdim[2]+1.936491673103709*alpha_vdim[1])*f[2]; 
  out[8] += 2.958039891549809*(alpha_cdim[2]*f[6]+alpha_cdim[0]*f[4])+1.322875655532295*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[9] += 2.958039891549809*(alpha_vdim[8]*f[13]+alpha_vdim[4]*f[10])+1.322875655532295*alpha_vdim[8]*f[8]+2.958039891549809*(alpha_vdim[1]*f[7]+alpha_vdim[0]*f[5])+1.322875655532295*(alpha_vdim[4]*f[4]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[10] += 1.700840128541522*alpha_cdim[2]*f[12]+1.154700538379252*alpha_vdim[8]*f[11]+1.700840128541522*(alpha_vdim[1]*f[11]+f[3]*alpha_vdim[8])+1.936491673103709*alpha_cdim[0]*f[7]+1.237179148263484*alpha_vdim[4]*f[6]+1.936491673103709*(alpha_vdim[0]*f[6]+f[2]*alpha_vdim[4])+1.732050807568877*(alpha_cdim[2]+alpha_vdim[1])*f[3]; 
  out[11] += 2.645751311064591*alpha_cdim[2]*f[10]+(0.5163977794943223*alpha_vdim[4]+0.8660254037844386*alpha_vdim[0])*f[8]+(0.5163977794943223*f[4]+0.8660254037844386*f[0])*alpha_vdim[8]+2.958039891549809*alpha_cdim[0]*f[6]+alpha_cdim[2]*(1.183215956619923*f[5]+2.958039891549809*f[4])+0.7606388292556648*(alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+1.322875655532295*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[12] += 2.598076211353316*alpha_vdim[4]*f[13]+(2.598076211353316*alpha_vdim[8]+2.645751311064591*alpha_vdim[1])*f[10]+0.8660254037844386*alpha_cdim[0]*f[9]+1.161895003862225*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8])+(2.645751311064591*alpha_vdim[4]+2.958039891549809*alpha_vdim[0])*f[7]+(0.7606388292556648*alpha_cdim[2]+2.958039891549809*alpha_vdim[1])*f[5]+1.183215956619923*(alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+1.322875655532295*(alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[13] += 2.598076211353316*alpha_cdim[2]*f[14]+(1.154700538379252*alpha_vdim[4]+1.936491673103709*alpha_vdim[0])*f[11]+2.958039891549809*alpha_cdim[0]*f[10]+1.161895003862225*alpha_cdim[2]*f[9]+(1.154700538379252*f[6]+1.936491673103709*f[2])*alpha_vdim[8]+(2.645751311064591*alpha_cdim[2]+1.700840128541522*alpha_vdim[1])*f[6]+1.322875655532295*alpha_cdim[0]*f[5]+1.700840128541522*f[3]*alpha_vdim[4]+1.183215956619923*alpha_cdim[2]*f[2]; 
  out[14] += (1.763834207376394*alpha_vdim[8]+2.598076211353316*alpha_vdim[1])*f[13]+1.936491673103709*alpha_cdim[0]*f[12]+(1.889822365046136*alpha_vdim[4]+2.958039891549809*alpha_vdim[0])*f[10]+(0.7888106377466154*alpha_vdim[8]+1.161895003862225*alpha_vdim[1])*f[8]+(2.598076211353316*f[7]+1.161895003862225*f[1])*alpha_vdim[8]+(1.700840128541522*alpha_cdim[2]+2.645751311064591*alpha_vdim[1])*f[7]+alpha_vdim[4]*(2.958039891549809*f[5]+0.8451542547285166*f[4])+1.322875655532295*(alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4])+1.183215956619923*alpha_vdim[1]*f[1]; 
  out[15] += 2.958039891549809*alpha_cdim[0]*f[14]+(1.763834207376394*alpha_vdim[4]+2.958039891549809*alpha_vdim[0])*f[13]+(1.763834207376394*alpha_vdim[8]+2.598076211353316*(alpha_cdim[2]+alpha_vdim[1]))*f[10]+1.322875655532295*alpha_cdim[0]*f[9]+(0.7888106377466154*alpha_vdim[4]+1.322875655532295*alpha_vdim[0])*f[8]+(2.958039891549809*f[5]+0.7888106377466154*f[4]+1.322875655532295*f[0])*alpha_vdim[8]+2.598076211353316*alpha_vdim[4]*f[7]+1.161895003862225*(alpha_cdim[2]*f[5]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4]); 

  return alpha_mid; 
} 
