#include <VlasovModDecl.h> 
double VlasovVolElcMag1x1vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double alpha1[3]; 

  alpha1[0] = 1.414213562373095*E0[0]*dv10; 
  alpha1[1] = 1.414213562373095*E0[1]*dv10; 
  const double amid1 = 0.5*alpha1[0]; 
  out[2] += 0.8660254037844386*(alpha1[1]*f[1]+alpha1[0]*f[0]); 
return 0.5*std::abs(amid1); 
} 
double VlasovVolElcMag1x1vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double alpha1[6]; 

  alpha1[0] = 1.414213562373095*E0[0]*dv10; 
  alpha1[1] = 1.414213562373095*E0[1]*dv10; 
  alpha1[4] = 1.414213562373095*E0[2]*dv10; 
  const double amid1 = 0.5*alpha1[0]-0.5590169943749475*alpha1[4]; 
  out[2] += 0.8660254037844386*(alpha1[4]*f[4]+alpha1[1]*f[1]+alpha1[0]*f[0]); 
  out[3] += 0.7745966692414833*(alpha1[1]*f[4]+f[1]*alpha1[4])+0.8660254037844386*(alpha1[0]*f[1]+f[0]*alpha1[1]); 
  out[5] += 1.936491673103709*(alpha1[1]*f[3]+alpha1[0]*f[2]); 
return 0.5*std::abs(amid1); 
} 
double VlasovVolElcMag1x1vMaxP3(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double alpha1[10]; 

  alpha1[0] = 1.414213562373095*E0[0]*dv10; 
  alpha1[1] = 1.414213562373095*E0[1]*dv10; 
  alpha1[4] = 1.414213562373095*E0[2]*dv10; 
  alpha1[8] = 1.414213562373095*E0[3]*dv10; 
  const double amid1 = 0.5*alpha1[0]-0.5590169943749475*alpha1[4]; 
  out[2] += 0.8660254037844386*(alpha1[8]*f[8]+alpha1[4]*f[4]+alpha1[1]*f[1]+alpha1[0]*f[0]); 
  out[3] += 0.7606388292556648*(alpha1[4]*f[8]+f[4]*alpha1[8])+0.7745966692414833*(alpha1[1]*f[4]+f[1]*alpha1[4])+0.8660254037844386*(alpha1[0]*f[1]+f[0]*alpha1[1]); 
  out[5] += 1.936491673103709*(alpha1[4]*f[6]+alpha1[1]*f[3]+alpha1[0]*f[2]); 
  out[6] += 0.5163977794943223*alpha1[8]*f[8]+0.7606388292556648*(alpha1[1]*f[8]+f[1]*alpha1[8])+0.5532833351724881*alpha1[4]*f[4]+0.8660254037844386*(alpha1[0]*f[4]+f[0]*alpha1[4])+0.7745966692414833*alpha1[1]*f[1]; 
  out[7] += 1.700840128541522*f[6]*alpha1[8]+1.732050807568877*(alpha1[1]*f[6]+f[3]*alpha1[4])+1.936491673103709*(alpha1[0]*f[3]+alpha1[1]*f[2]); 
  out[9] += 1.322875655532295*alpha1[8]*f[8]+2.958039891549809*(alpha1[1]*f[7]+alpha1[0]*f[5])+1.322875655532295*(alpha1[4]*f[4]+alpha1[1]*f[1]+alpha1[0]*f[0]); 
return 0.5*std::abs(amid1); 
} 
