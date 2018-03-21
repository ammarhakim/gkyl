#include <cmath> 
#include <CanonicalModDecl.h> 
double CanonicalVol2x2vMaxP1(const double *w, const double *dxv, const double *H, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double dxv2 = dxv[2]; 
  double dxv3 = dxv[3]; 
  double dxInv0 = 1.0/dxv[0]; 
  double dxInv1 = 1.0/dxv[1]; 
  double dvInv0 = 1.0/dxv[2]; 
  double dvInv1 = 1.0/dxv[3]; 
  out[1] += 0.1875*f[0]*H[3]*dxv1*dxv3; 
  out[2] += 0.1875*f[0]*H[4]*dxv0*dxv2; 
  out[3] += -0.1875*f[0]*H[1]*dxv1*dxv3; 
  out[4] += -0.1875*f[0]*H[2]*dxv0*dxv2; 
  double cflFreq = 0.0; 
  // calculate phase space speed alpha in each direction, and add to cfl 
  // here alpha is cell-averaged 
  double alpha = 0.0; 
  alpha = fabs(1.732050807568877*H[3]*dxv1*dxv3);  // alpha_x0 
  cflFreq += alpha*dxInv0; 
  alpha = fabs(1.732050807568877*H[4]*dxv0*dxv2);  // alpha_x1 
  cflFreq += alpha*dxInv1; 
  alpha = fabs(-1.732050807568877*H[1]*dxv1*dxv3);  // alpha_v0 
  cflFreq += alpha*dvInv0; 
  alpha = fabs(-1.732050807568877*H[2]*dxv0*dxv2);  // alpha_v1 
  cflFreq += alpha*dvInv1; 
  return cflFreq; 
} 
double CanonicalVol2x2vMaxP2(const double *w, const double *dxv, const double *H, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double dxv2 = dxv[2]; 
  double dxv3 = dxv[3]; 
  double dxInv0 = 1.0/dxv[0]; 
  double dxInv1 = 1.0/dxv[1]; 
  double dvInv0 = 1.0/dxv[2]; 
  double dvInv1 = 1.0/dxv[3]; 
  out[1] += 0.1875*(2.23606797749979*f[3]*H[13]+f[4]*H[10]+f[2]*H[7]+f[1]*H[6]+f[0]*H[3])*dxv1*dxv3; 
  out[2] += 0.1875*(2.23606797749979*f[4]*H[14]+f[3]*H[10]+f[2]*H[9]+f[1]*H[8]+f[0]*H[4])*dxv0*dxv2; 
  out[3] += -0.1875*(2.23606797749979*f[1]*H[11]+f[4]*H[8]+f[3]*H[6]+f[2]*H[5]+f[0]*H[1])*dxv1*dxv3; 
  out[4] += -0.1875*(2.23606797749979*f[2]*H[12]+f[4]*H[9]+f[3]*H[7]+f[1]*H[5]+f[0]*H[2])*dxv0*dxv2; 
  out[5] += 0.0125*(3.0*(11.18033988749895*f[7]*H[13]+4.47213595499958*H[7]*f[12]+5.0*f[9]*H[10]+5.0*f[0]*H[7]+5.0*f[5]*H[6]+5.0*f[2]*H[3])*dxv1*dxv3+3.0*(11.18033988749895*f[8]*H[14]+4.47213595499958*H[8]*f[11]+5.0*f[6]*H[10]+5.0*f[5]*H[9]+5.0*f[0]*H[8]+5.0*f[1]*H[4])*dxv0*dxv2); 
  out[6] += 0.1875*(2.0*H[13]*f[13]+2.23606797749979*f[0]*H[13]-2.0*H[11]*f[11]-2.23606797749979*f[0]*H[11]+H[10]*f[10]-1.0*H[8]*f[8]+H[7]*f[7]-1.0*H[5]*f[5]+H[3]*f[3]-1.0*H[1]*f[1])*dxv1*dxv3; 
  out[7] += -0.0125*(3.0*(4.47213595499958*H[5]*f[12]+11.18033988749895*f[5]*H[11]+5.0*H[8]*f[9]+5.0*H[6]*f[7]+5.0*f[0]*H[5]+5.0*H[1]*f[2])*dxv1*dxv3-3.0*(11.18033988749895*f[10]*H[14]+4.47213595499958*H[10]*f[13]+5.0*f[0]*H[10]+5.0*f[7]*H[9]+5.0*f[6]*H[8]+5.0*f[3]*H[4])*dxv0*dxv2); 
  out[8] += 0.0125*(3.0*(4.47213595499958*H[10]*f[14]+11.18033988749895*f[10]*H[13]+5.0*f[0]*H[10]+5.0*H[7]*f[9]+5.0*H[6]*f[8]+5.0*H[3]*f[4])*dxv1*dxv3-3.0*(11.18033988749895*f[5]*H[12]+4.47213595499958*H[5]*f[11]+5.0*f[8]*H[9]+5.0*f[6]*H[7]+5.0*f[0]*H[5]+5.0*f[1]*H[2])*dxv0*dxv2); 
  out[9] += 0.1875*(2.0*H[14]*f[14]+2.23606797749979*f[0]*H[14]-2.0*H[12]*f[12]-2.23606797749979*f[0]*H[12]+H[10]*f[10]+H[8]*f[8]-1.0*H[7]*f[7]-1.0*H[5]*f[5]+H[4]*f[4]-1.0*H[2]*f[2])*dxv0*dxv2; 
  out[10] += -0.0125*(3.0*(4.47213595499958*H[8]*f[14]+11.18033988749895*f[8]*H[11]+5.0*H[6]*f[10]+5.0*H[5]*f[9]+5.0*f[0]*H[8]+5.0*H[1]*f[4])*dxv1*dxv3+3.0*(4.47213595499958*H[7]*f[13]+11.18033988749895*f[7]*H[12]+5.0*H[9]*f[10]+5.0*f[0]*H[7]+5.0*H[5]*f[6]+5.0*H[2]*f[3])*dxv0*dxv2); 
  out[11] += 0.1875*(5.0*f[6]*H[13]+2.0*H[6]*f[11]+2.23606797749979*f[8]*H[10]+2.23606797749979*f[5]*H[7]+2.23606797749979*f[0]*H[6]+2.23606797749979*f[1]*H[3])*dxv1*dxv3; 
  out[12] += 0.1875*(5.0*f[9]*H[14]+2.0*H[9]*f[12]+2.23606797749979*f[7]*H[10]+2.23606797749979*f[0]*H[9]+2.23606797749979*f[5]*H[8]+2.23606797749979*f[2]*H[4])*dxv0*dxv2; 
  out[13] += -0.1875*(2.0*H[6]*f[13]+5.0*f[6]*H[11]+2.23606797749979*H[8]*f[10]+2.23606797749979*H[5]*f[7]+2.23606797749979*f[0]*H[6]+2.23606797749979*H[1]*f[3])*dxv1*dxv3; 
  out[14] += -0.1875*(2.0*H[9]*f[14]+5.0*f[9]*H[12]+2.23606797749979*H[7]*f[10]+2.23606797749979*f[0]*H[9]+2.23606797749979*H[5]*f[8]+2.23606797749979*H[2]*f[4])*dxv0*dxv2; 
  double cflFreq = 0.0; 
  // calculate phase space speed alpha in each direction, and add to cfl 
  // here alpha is cell-averaged 
  double alpha = 0.0; 
  alpha = fabs(1.732050807568877*H[3]*dxv1*dxv3);  // alpha_x0 
  cflFreq += alpha*dxInv0; 
  alpha = fabs(1.732050807568877*H[4]*dxv0*dxv2);  // alpha_x1 
  cflFreq += alpha*dxInv1; 
  alpha = fabs(-1.732050807568877*H[1]*dxv1*dxv3);  // alpha_v0 
  cflFreq += alpha*dvInv0; 
  alpha = fabs(-1.732050807568877*H[2]*dxv0*dxv2);  // alpha_v1 
  cflFreq += alpha*dvInv1; 
  return cflFreq; 
} 
