#include <VlasovModDecl.h> 
__host__ __device__ double VlasovUpwindSurfElcMag2x2vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[2]; 
  double dv10r = 2/dxvr[2]; 

  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double *E0 = &EM[0]; 
  const double *B2 = &EM[20]; 

  double alpha[8]; 
  double incr[16]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 1.414213562373095*(B2[2]*wv2+E0[2]); 
  alpha[3] = 0.408248290463863*B2[0]*dv2; 
  alpha[4] = 1.414213562373095*(B2[3]*wv2+E0[3]); 
  alpha[5] = 0.408248290463863*B2[1]*dv2; 
  alpha[6] = 0.408248290463863*B2[2]*dv2; 
  alpha[7] = 0.408248290463863*B2[3]*dv2; 

  const double amid = 0.3535533905932737*alpha[0]; 

  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = (-0.3535533905932737*alpha[7])+0.3535533905932737*(alpha[6]+alpha[5]+alpha[4])-0.3535533905932737*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fr[15]+fl[15]))+0.4330127018922193*(fr[14]+fl[14]+fr[13]+fl[13])+0.25*fr[12]-0.25*fl[12]+0.4330127018922193*(fr[11]+fl[11])-0.4330127018922193*(fr[10]+fl[10])-0.25*(fr[9]+fr[8])+0.25*(fl[9]+fl[8])-0.4330127018922193*(fr[7]+fl[7]+fr[6]+fl[6])-0.25*fr[5]+0.25*(fl[5]+fr[4])-0.25*fl[4]+0.4330127018922193*(fr[3]+fl[3])+0.25*(fr[2]+fr[1])-0.25*(fl[2]+fl[1]+fr[0])+0.25*fl[0])*sgn(alphaOrdR)+0.4330127018922193*fr[15]-0.4330127018922193*(fl[15]+fr[14]+fr[13])+0.4330127018922193*(fl[14]+fl[13])-0.25*(fr[12]+fl[12])-0.4330127018922193*fr[11]+0.4330127018922193*(fl[11]+fr[10])-0.4330127018922193*fl[10]+0.25*(fr[9]+fl[9]+fr[8]+fl[8])+0.4330127018922193*(fr[7]+fr[6])-0.4330127018922193*(fl[7]+fl[6])+0.25*(fr[5]+fl[5])-0.25*(fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]-0.25*(fr[2]+fl[2]+fr[1]+fl[1])+0.25*(fr[0]+fl[0])); 
  alphaOrdR = 0.3535533905932737*(alpha[7]+alpha[6])-0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14])-0.4330127018922193*(fr[13]+fl[13])-0.25*fr[12]+0.25*fl[12]-0.4330127018922193*(fr[11]+fl[11]+fr[10]+fl[10])-0.25*fr[9]+0.25*(fl[9]+fr[8])-0.25*fl[8]-0.4330127018922193*(fr[7]+fl[7])+0.4330127018922193*(fr[6]+fl[6])+0.25*(fr[5]+fr[4])-0.25*(fl[5]+fl[4])+0.4330127018922193*(fr[3]+fl[3])+0.25*fr[2]-0.25*(fl[2]+fr[1]+fr[0])+0.25*(fl[1]+fl[0]))*sgn(alphaOrdR)-0.4330127018922193*(fr[15]+fr[14])+0.4330127018922193*(fl[15]+fl[14]+fr[13])-0.4330127018922193*fl[13]+0.25*(fr[12]+fl[12])+0.4330127018922193*(fr[11]+fr[10])-0.4330127018922193*(fl[11]+fl[10])+0.25*(fr[9]+fl[9])-0.25*(fr[8]+fl[8])+0.4330127018922193*fr[7]-0.4330127018922193*(fl[7]+fr[6])+0.4330127018922193*fl[6]-0.25*(fr[5]+fl[5]+fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]-0.25*(fr[2]+fl[2])+0.25*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaOrdR = 0.3535533905932737*alpha[7]-0.3535533905932737*alpha[6]+0.3535533905932737*alpha[5]-0.3535533905932737*(alpha[4]+alpha[3])+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fr[15]+fl[15])-0.4330127018922193*(fr[14]+fl[14])+0.4330127018922193*(fr[13]+fl[13])-0.25*fr[12]+0.25*fl[12]-0.4330127018922193*(fr[11]+fl[11]+fr[10]+fl[10])+0.25*fr[9]-0.25*(fl[9]+fr[8])+0.25*fl[8]+0.4330127018922193*(fr[7]+fl[7])-0.4330127018922193*(fr[6]+fl[6])+0.25*(fr[5]+fr[4])-0.25*(fl[5]+fl[4])+0.4330127018922193*(fr[3]+fl[3])-0.25*fr[2]+0.25*(fl[2]+fr[1])-0.25*(fl[1]+fr[0])+0.25*fl[0])*sgn(alphaOrdR)-0.4330127018922193*fr[15]+0.4330127018922193*(fl[15]+fr[14])-0.4330127018922193*(fl[14]+fr[13])+0.4330127018922193*fl[13]+0.25*(fr[12]+fl[12])+0.4330127018922193*(fr[11]+fr[10])-0.4330127018922193*(fl[11]+fl[10])-0.25*(fr[9]+fl[9])+0.25*(fr[8]+fl[8])-0.4330127018922193*fr[7]+0.4330127018922193*(fl[7]+fr[6])-0.4330127018922193*fl[6]-0.25*(fr[5]+fl[5]+fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]+0.25*(fr[2]+fl[2])-0.25*(fr[1]+fl[1])+0.25*(fr[0]+fl[0])); 
  alphaOrdR = (-0.3535533905932737*(alpha[7]+alpha[6]+alpha[5]))+0.3535533905932737*alpha[4]-0.3535533905932737*alpha[3]+0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]))+0.25*fr[12]-0.25*fl[12]+0.4330127018922193*(fr[11]+fl[11])-0.4330127018922193*(fr[10]+fl[10])+0.25*(fr[9]+fr[8])-0.25*(fl[9]+fl[8])+0.4330127018922193*(fr[7]+fl[7]+fr[6]+fl[6])-0.25*fr[5]+0.25*(fl[5]+fr[4])-0.25*fl[4]+0.4330127018922193*(fr[3]+fl[3])-0.25*(fr[2]+fr[1]+fr[0])+0.25*(fl[2]+fl[1]+fl[0]))*sgn(alphaOrdR)+0.4330127018922193*(fr[15]+fr[14]+fr[13])-0.4330127018922193*(fl[15]+fl[14]+fl[13])-0.25*(fr[12]+fl[12])-0.4330127018922193*fr[11]+0.4330127018922193*(fl[11]+fr[10])-0.4330127018922193*fl[10]-0.25*(fr[9]+fl[9]+fr[8]+fl[8])-0.4330127018922193*(fr[7]+fr[6])+0.4330127018922193*(fl[7]+fl[6])+0.25*(fr[5]+fl[5])-0.25*(fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]+0.25*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaOrdR = 0.3535533905932737*alpha[7]-0.3535533905932737*(alpha[6]+alpha[5])+0.3535533905932737*(alpha[4]+alpha[3])-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fr[15]+fl[15])-0.4330127018922193*(fr[14]+fl[14]+fr[13]+fl[13])-0.25*fr[12]+0.25*fl[12]+0.4330127018922193*(fr[11]+fl[11]+fr[10]+fl[10])+0.25*(fr[9]+fr[8])-0.25*(fl[9]+fl[8])-0.4330127018922193*(fr[7]+fl[7]+fr[6]+fl[6])-0.25*(fr[5]+fr[4])+0.25*(fl[5]+fl[4])+0.4330127018922193*(fr[3]+fl[3])+0.25*(fr[2]+fr[1])-0.25*(fl[2]+fl[1]+fr[0])+0.25*fl[0])*sgn(alphaOrdR)-0.4330127018922193*fr[15]+0.4330127018922193*(fl[15]+fr[14]+fr[13])-0.4330127018922193*(fl[14]+fl[13])+0.25*(fr[12]+fl[12])-0.4330127018922193*(fr[11]+fr[10])+0.4330127018922193*(fl[11]+fl[10])-0.25*(fr[9]+fl[9]+fr[8]+fl[8])+0.4330127018922193*(fr[7]+fr[6])-0.4330127018922193*(fl[7]+fl[6])+0.25*(fr[5]+fl[5]+fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]-0.25*(fr[2]+fl[2]+fr[1]+fl[1])+0.25*(fr[0]+fl[0])); 
  alphaOrdR = (-0.3535533905932737*(alpha[7]+alpha[6]))+0.3535533905932737*alpha[5]-0.3535533905932737*alpha[4]+0.3535533905932737*alpha[3]-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]))+0.4330127018922193*(fr[13]+fl[13])+0.25*fr[12]-0.25*fl[12]-0.4330127018922193*(fr[11]+fl[11])+0.4330127018922193*(fr[10]+fl[10])+0.25*fr[9]-0.25*(fl[9]+fr[8])+0.25*fl[8]-0.4330127018922193*(fr[7]+fl[7])+0.4330127018922193*(fr[6]+fl[6])+0.25*fr[5]-0.25*(fl[5]+fr[4])+0.25*fl[4]+0.4330127018922193*(fr[3]+fl[3])+0.25*fr[2]-0.25*(fl[2]+fr[1]+fr[0])+0.25*(fl[1]+fl[0]))*sgn(alphaOrdR)+0.4330127018922193*(fr[15]+fr[14])-0.4330127018922193*(fl[15]+fl[14]+fr[13])+0.4330127018922193*fl[13]-0.25*(fr[12]+fl[12])+0.4330127018922193*fr[11]-0.4330127018922193*(fl[11]+fr[10])+0.4330127018922193*fl[10]-0.25*(fr[9]+fl[9])+0.25*(fr[8]+fl[8])+0.4330127018922193*fr[7]-0.4330127018922193*(fl[7]+fr[6])+0.4330127018922193*fl[6]-0.25*(fr[5]+fl[5])+0.25*(fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]-0.25*(fr[2]+fl[2])+0.25*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaOrdR = (-0.3535533905932737*alpha[7])+0.3535533905932737*alpha[6]-0.3535533905932737*(alpha[5]+alpha[4])+0.3535533905932737*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fr[15]+fl[15]))+0.4330127018922193*(fr[14]+fl[14])-0.4330127018922193*(fr[13]+fl[13])+0.25*fr[12]-0.25*fl[12]-0.4330127018922193*(fr[11]+fl[11])+0.4330127018922193*(fr[10]+fl[10])-0.25*fr[9]+0.25*(fl[9]+fr[8])-0.25*fl[8]+0.4330127018922193*(fr[7]+fl[7])-0.4330127018922193*(fr[6]+fl[6])+0.25*fr[5]-0.25*(fl[5]+fr[4])+0.25*fl[4]+0.4330127018922193*(fr[3]+fl[3])-0.25*fr[2]+0.25*(fl[2]+fr[1])-0.25*(fl[1]+fr[0])+0.25*fl[0])*sgn(alphaOrdR)+0.4330127018922193*fr[15]-0.4330127018922193*(fl[15]+fr[14])+0.4330127018922193*(fl[14]+fr[13])-0.4330127018922193*fl[13]-0.25*(fr[12]+fl[12])+0.4330127018922193*fr[11]-0.4330127018922193*(fl[11]+fr[10])+0.4330127018922193*fl[10]+0.25*(fr[9]+fl[9])-0.25*(fr[8]+fl[8])-0.4330127018922193*fr[7]+0.4330127018922193*(fl[7]+fr[6])-0.4330127018922193*fl[6]-0.25*(fr[5]+fl[5])+0.25*(fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]+0.25*(fr[2]+fl[2])-0.25*(fr[1]+fl[1])+0.25*(fr[0]+fl[0])); 
  alphaOrdR = 0.3535533905932737*(alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13])-0.25*fr[12]+0.25*fl[12]+0.4330127018922193*(fr[11]+fl[11]+fr[10]+fl[10])-0.25*(fr[9]+fr[8])+0.25*(fl[9]+fl[8])+0.4330127018922193*(fr[7]+fl[7]+fr[6]+fl[6])-0.25*(fr[5]+fr[4])+0.25*(fl[5]+fl[4])+0.4330127018922193*(fr[3]+fl[3])-0.25*(fr[2]+fr[1]+fr[0])+0.25*(fl[2]+fl[1]+fl[0]))*sgn(alphaOrdR)-0.4330127018922193*(fr[15]+fr[14]+fr[13])+0.4330127018922193*(fl[15]+fl[14]+fl[13])+0.25*(fr[12]+fl[12])-0.4330127018922193*(fr[11]+fr[10])+0.4330127018922193*(fl[11]+fl[10])+0.25*(fr[9]+fl[9]+fr[8]+fl[8])-0.4330127018922193*(fr[7]+fr[6])+0.4330127018922193*(fl[7]+fl[6])+0.25*(fr[5]+fl[5]+fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]+0.25*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alpha[7]*fUp[7]+alpha[6]*fUp[6]+alpha[5]*fUp[5]+alpha[4]*fUp[4]+alpha[3]*fUp[3]+alpha[2]*fUp[2]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[1] = 0.25*(alpha[6]*fUp[7]+fUp[6]*alpha[7]+alpha[3]*fUp[5]+fUp[3]*alpha[5]+alpha[2]*fUp[4]+fUp[2]*alpha[4]+alpha[0]*fUp[1]+fUp[0]*alpha[1]); 
  incr[2] = 0.25*(alpha[5]*fUp[7]+fUp[5]*alpha[7]+alpha[3]*fUp[6]+fUp[3]*alpha[6]+alpha[1]*fUp[4]+fUp[1]*alpha[4]+alpha[0]*fUp[2]+fUp[0]*alpha[2]); 
  incr[3] = -0.4330127018922193*(alpha[7]*fUp[7]+alpha[6]*fUp[6]+alpha[5]*fUp[5]+alpha[4]*fUp[4]+alpha[3]*fUp[3]+alpha[2]*fUp[2]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[4] = 0.25*(alpha[4]*fUp[7]+fUp[4]*alpha[7]+alpha[2]*fUp[6]+fUp[2]*alpha[6]+alpha[1]*fUp[5]+fUp[1]*alpha[5]+alpha[0]*fUp[3]+fUp[0]*alpha[3]); 
  incr[5] = 0.25*(alpha[3]*fUp[7]+fUp[3]*alpha[7]+alpha[5]*fUp[6]+fUp[5]*alpha[6]+alpha[0]*fUp[4]+fUp[0]*alpha[4]+alpha[1]*fUp[2]+fUp[1]*alpha[2]); 
  incr[6] = -0.4330127018922193*(alpha[6]*fUp[7]+fUp[6]*alpha[7]+alpha[3]*fUp[5]+fUp[3]*alpha[5]+alpha[2]*fUp[4]+fUp[2]*alpha[4]+alpha[0]*fUp[1]+fUp[0]*alpha[1]); 
  incr[7] = -0.4330127018922193*(alpha[5]*fUp[7]+fUp[5]*alpha[7]+alpha[3]*fUp[6]+fUp[3]*alpha[6]+alpha[1]*fUp[4]+fUp[1]*alpha[4]+alpha[0]*fUp[2]+fUp[0]*alpha[2]); 
  incr[8] = 0.25*(alpha[2]*fUp[7]+fUp[2]*alpha[7]+alpha[4]*fUp[6]+fUp[4]*alpha[6]+alpha[0]*fUp[5]+fUp[0]*alpha[5]+alpha[1]*fUp[3]+fUp[1]*alpha[3]); 
  incr[9] = 0.25*(alpha[1]*fUp[7]+fUp[1]*alpha[7]+alpha[0]*fUp[6]+fUp[0]*alpha[6]+alpha[4]*fUp[5]+fUp[4]*alpha[5]+alpha[2]*fUp[3]+fUp[2]*alpha[3]); 
  incr[10] = -0.4330127018922193*(alpha[4]*fUp[7]+fUp[4]*alpha[7]+alpha[2]*fUp[6]+fUp[2]*alpha[6]+alpha[1]*fUp[5]+fUp[1]*alpha[5]+alpha[0]*fUp[3]+fUp[0]*alpha[3]); 
  incr[11] = -0.4330127018922193*(alpha[3]*fUp[7]+fUp[3]*alpha[7]+alpha[5]*fUp[6]+fUp[5]*alpha[6]+alpha[0]*fUp[4]+fUp[0]*alpha[4]+alpha[1]*fUp[2]+fUp[1]*alpha[2]); 
  incr[12] = 0.25*(alpha[0]*fUp[7]+fUp[0]*alpha[7]+alpha[1]*fUp[6]+fUp[1]*alpha[6]+alpha[2]*fUp[5]+fUp[2]*alpha[5]+alpha[3]*fUp[4]+fUp[3]*alpha[4]); 
  incr[13] = -0.4330127018922193*(alpha[2]*fUp[7]+fUp[2]*alpha[7]+alpha[4]*fUp[6]+fUp[4]*alpha[6]+alpha[0]*fUp[5]+fUp[0]*alpha[5]+alpha[1]*fUp[3]+fUp[1]*alpha[3]); 
  incr[14] = -0.4330127018922193*(alpha[1]*fUp[7]+fUp[1]*alpha[7]+alpha[0]*fUp[6]+fUp[0]*alpha[6]+alpha[4]*fUp[5]+fUp[4]*alpha[5]+alpha[2]*fUp[3]+fUp[2]*alpha[3]); 
  incr[15] = -0.4330127018922193*(alpha[0]*fUp[7]+fUp[0]*alpha[7]+alpha[1]*fUp[6]+fUp[1]*alpha[6]+alpha[2]*fUp[5]+fUp[2]*alpha[5]+alpha[3]*fUp[4]+fUp[3]*alpha[4]); 

  outr[0] += incr[0]*dv10r; 
  outr[1] += incr[1]*dv10r; 
  outr[2] += incr[2]*dv10r; 
  outr[3] += incr[3]*dv10r; 
  outr[4] += incr[4]*dv10r; 
  outr[5] += incr[5]*dv10r; 
  outr[6] += incr[6]*dv10r; 
  outr[7] += incr[7]*dv10r; 
  outr[8] += incr[8]*dv10r; 
  outr[9] += incr[9]*dv10r; 
  outr[10] += incr[10]*dv10r; 
  outr[11] += incr[11]*dv10r; 
  outr[12] += incr[12]*dv10r; 
  outr[13] += incr[13]*dv10r; 
  outr[14] += incr[14]*dv10r; 
  outr[15] += incr[15]*dv10r; 

  outl[0] += -1.0*incr[0]*dv10l; 
  outl[1] += -1.0*incr[1]*dv10l; 
  outl[2] += -1.0*incr[2]*dv10l; 
  outl[3] += incr[3]*dv10l; 
  outl[4] += -1.0*incr[4]*dv10l; 
  outl[5] += -1.0*incr[5]*dv10l; 
  outl[6] += incr[6]*dv10l; 
  outl[7] += incr[7]*dv10l; 
  outl[8] += -1.0*incr[8]*dv10l; 
  outl[9] += -1.0*incr[9]*dv10l; 
  outl[10] += incr[10]*dv10l; 
  outl[11] += incr[11]*dv10l; 
  outl[12] += -1.0*incr[12]*dv10l; 
  outl[13] += incr[13]*dv10l; 
  outl[14] += incr[14]*dv10l; 
  outl[15] += incr[15]*dv10l; 


  return std::abs(amid); 
} 
__host__ __device__ double VlasovUpwindSurfElcMag2x2vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[3]; 
  double dv11r = 2/dxvr[3]; 

  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double *E1 = &EM[4]; 
  const double *B2 = &EM[20]; 

  double alpha[8]; 
  double incr[16]; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = 1.414213562373095*E1[2]-1.414213562373095*B2[2]*wv1; 
  alpha[3] = -0.408248290463863*B2[0]*dv1; 
  alpha[4] = 1.414213562373095*E1[3]-1.414213562373095*B2[3]*wv1; 
  alpha[5] = -0.408248290463863*B2[1]*dv1; 
  alpha[6] = -0.408248290463863*B2[2]*dv1; 
  alpha[7] = -0.408248290463863*B2[3]*dv1; 

  const double amid = 0.3535533905932737*alpha[0]; 

  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = (-0.3535533905932737*alpha[7])+0.3535533905932737*(alpha[6]+alpha[5]+alpha[4])-0.3535533905932737*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fr[15]+fl[15]))+0.4330127018922193*(fr[14]+fl[14]+fr[13]+fl[13]+fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9]+fr[8]+fl[8])-0.25*(fr[7]+fr[6]+fr[5])+0.25*(fl[7]+fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])+0.25*(fr[3]+fr[2]+fr[1])-0.25*(fl[3]+fl[2]+fl[1]+fr[0])+0.25*fl[0])*sgn(alphaOrdR)+0.4330127018922193*fr[15]-0.4330127018922193*(fl[15]+fr[14]+fr[13]+fr[12])+0.4330127018922193*(fl[14]+fl[13]+fl[12])-0.25*(fr[11]+fl[11])+0.4330127018922193*(fr[10]+fr[9]+fr[8])-0.4330127018922193*(fl[10]+fl[9]+fl[8])+0.25*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3]+fr[2]+fl[2]+fr[1]+fl[1])+0.25*(fr[0]+fl[0])); 
  alphaOrdR = 0.3535533905932737*(alpha[7]+alpha[6])-0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2])+0.3535533905932737*(alpha[1]+alpha[0]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14])-0.4330127018922193*(fr[13]+fl[13]+fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9])+0.4330127018922193*(fr[8]+fl[8])-0.25*fr[7]+0.25*(fl[7]+fr[6]+fr[5])-0.25*(fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])+0.25*(fr[3]+fr[2])-0.25*(fl[3]+fl[2]+fr[1]+fr[0])+0.25*(fl[1]+fl[0]))*sgn(alphaOrdR)-0.4330127018922193*(fr[15]+fr[14])+0.4330127018922193*(fl[15]+fl[14]+fr[13]+fr[12])-0.4330127018922193*(fl[13]+fl[12])+0.25*(fr[11]+fl[11])+0.4330127018922193*(fr[10]+fr[9])-0.4330127018922193*(fl[10]+fl[9]+fr[8])+0.4330127018922193*fl[8]+0.25*(fr[7]+fl[7])-0.25*(fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3]+fr[2]+fl[2])+0.25*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaOrdR = 0.3535533905932737*alpha[7]-0.3535533905932737*alpha[6]+0.3535533905932737*alpha[5]-0.3535533905932737*(alpha[4]+alpha[3])+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fr[15]+fl[15])-0.4330127018922193*(fr[14]+fl[14])+0.4330127018922193*(fr[13]+fl[13])-0.4330127018922193*(fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10])+0.4330127018922193*(fr[9]+fl[9])-0.4330127018922193*(fr[8]+fl[8])+0.25*fr[7]-0.25*(fl[7]+fr[6])+0.25*(fl[6]+fr[5])-0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])+0.25*fr[3]-0.25*(fl[3]+fr[2])+0.25*(fl[2]+fr[1])-0.25*(fl[1]+fr[0])+0.25*fl[0])*sgn(alphaOrdR)-0.4330127018922193*fr[15]+0.4330127018922193*(fl[15]+fr[14])-0.4330127018922193*(fl[14]+fr[13])+0.4330127018922193*(fl[13]+fr[12])-0.4330127018922193*fl[12]+0.25*(fr[11]+fl[11])+0.4330127018922193*fr[10]-0.4330127018922193*(fl[10]+fr[9])+0.4330127018922193*(fl[9]+fr[8])-0.4330127018922193*fl[8]-0.25*(fr[7]+fl[7])+0.25*(fr[6]+fl[6])-0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3])+0.25*(fr[2]+fl[2])-0.25*(fr[1]+fl[1])+0.25*(fr[0]+fl[0])); 
  alphaOrdR = (-0.3535533905932737*(alpha[7]+alpha[6]+alpha[5]))+0.3535533905932737*alpha[4]-0.3535533905932737*alpha[3]+0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]))+0.4330127018922193*(fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10])+0.4330127018922193*(fr[9]+fl[9]+fr[8]+fl[8])+0.25*(fr[7]+fr[6])-0.25*(fl[7]+fl[6]+fr[5])+0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])+0.25*fr[3]-0.25*(fl[3]+fr[2]+fr[1]+fr[0])+0.25*(fl[2]+fl[1]+fl[0]))*sgn(alphaOrdR)+0.4330127018922193*(fr[15]+fr[14]+fr[13])-0.4330127018922193*(fl[15]+fl[14]+fl[13]+fr[12])+0.4330127018922193*fl[12]-0.25*(fr[11]+fl[11])+0.4330127018922193*fr[10]-0.4330127018922193*(fl[10]+fr[9]+fr[8])+0.4330127018922193*(fl[9]+fl[8])-0.25*(fr[7]+fl[7]+fr[6]+fl[6])+0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3])+0.25*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaOrdR = 0.3535533905932737*alpha[7]-0.3535533905932737*(alpha[6]+alpha[5])+0.3535533905932737*(alpha[4]+alpha[3])-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fr[15]+fl[15])-0.4330127018922193*(fr[14]+fl[14]+fr[13]+fl[13])+0.4330127018922193*(fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10])-0.4330127018922193*(fr[9]+fl[9]+fr[8]+fl[8])+0.25*(fr[7]+fr[6])-0.25*(fl[7]+fl[6]+fr[5])+0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])-0.25*fr[3]+0.25*(fl[3]+fr[2]+fr[1])-0.25*(fl[2]+fl[1]+fr[0])+0.25*fl[0])*sgn(alphaOrdR)-0.4330127018922193*fr[15]+0.4330127018922193*(fl[15]+fr[14]+fr[13])-0.4330127018922193*(fl[14]+fl[13]+fr[12])+0.4330127018922193*fl[12]+0.25*(fr[11]+fl[11])-0.4330127018922193*fr[10]+0.4330127018922193*(fl[10]+fr[9]+fr[8])-0.4330127018922193*(fl[9]+fl[8])-0.25*(fr[7]+fl[7]+fr[6]+fl[6])+0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3])-0.25*(fr[2]+fl[2]+fr[1]+fl[1])+0.25*(fr[0]+fl[0])); 
  alphaOrdR = (-0.3535533905932737*(alpha[7]+alpha[6]))+0.3535533905932737*alpha[5]-0.3535533905932737*alpha[4]+0.3535533905932737*alpha[3]-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]))+0.4330127018922193*(fr[13]+fl[13])-0.4330127018922193*(fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10])-0.4330127018922193*(fr[9]+fl[9])+0.4330127018922193*(fr[8]+fl[8])+0.25*fr[7]-0.25*(fl[7]+fr[6])+0.25*(fl[6]+fr[5])-0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])-0.25*fr[3]+0.25*(fl[3]+fr[2])-0.25*(fl[2]+fr[1]+fr[0])+0.25*(fl[1]+fl[0]))*sgn(alphaOrdR)+0.4330127018922193*(fr[15]+fr[14])-0.4330127018922193*(fl[15]+fl[14]+fr[13])+0.4330127018922193*(fl[13]+fr[12])-0.4330127018922193*fl[12]-0.25*(fr[11]+fl[11])-0.4330127018922193*fr[10]+0.4330127018922193*(fl[10]+fr[9])-0.4330127018922193*(fl[9]+fr[8])+0.4330127018922193*fl[8]-0.25*(fr[7]+fl[7])+0.25*(fr[6]+fl[6])-0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3])-0.25*(fr[2]+fl[2])+0.25*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaOrdR = (-0.3535533905932737*alpha[7])+0.3535533905932737*alpha[6]-0.3535533905932737*(alpha[5]+alpha[4])+0.3535533905932737*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fr[15]+fl[15]))+0.4330127018922193*(fr[14]+fl[14])-0.4330127018922193*(fr[13]+fl[13]+fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9])-0.4330127018922193*(fr[8]+fl[8])-0.25*fr[7]+0.25*(fl[7]+fr[6]+fr[5])-0.25*(fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])-0.25*(fr[3]+fr[2])+0.25*(fl[3]+fl[2]+fr[1])-0.25*(fl[1]+fr[0])+0.25*fl[0])*sgn(alphaOrdR)+0.4330127018922193*fr[15]-0.4330127018922193*(fl[15]+fr[14])+0.4330127018922193*(fl[14]+fr[13]+fr[12])-0.4330127018922193*(fl[13]+fl[12])-0.25*(fr[11]+fl[11])-0.4330127018922193*(fr[10]+fr[9])+0.4330127018922193*(fl[10]+fl[9]+fr[8])-0.4330127018922193*fl[8]+0.25*(fr[7]+fl[7])-0.25*(fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3]+fr[2]+fl[2])-0.25*(fr[1]+fl[1])+0.25*(fr[0]+fl[0])); 
  alphaOrdR = 0.3535533905932737*(alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]+fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9]+fr[8]+fl[8])-0.25*(fr[7]+fr[6]+fr[5])+0.25*(fl[7]+fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])-0.25*(fr[3]+fr[2]+fr[1]+fr[0])+0.25*(fl[3]+fl[2]+fl[1]+fl[0]))*sgn(alphaOrdR)-0.4330127018922193*(fr[15]+fr[14]+fr[13]+fr[12])+0.4330127018922193*(fl[15]+fl[14]+fl[13]+fl[12])+0.25*(fr[11]+fl[11])-0.4330127018922193*(fr[10]+fr[9]+fr[8])+0.4330127018922193*(fl[10]+fl[9]+fl[8])+0.25*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3]+fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alpha[7]*fUp[7]+alpha[6]*fUp[6]+alpha[5]*fUp[5]+alpha[4]*fUp[4]+alpha[3]*fUp[3]+alpha[2]*fUp[2]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[1] = 0.25*(alpha[6]*fUp[7]+fUp[6]*alpha[7]+alpha[3]*fUp[5]+fUp[3]*alpha[5]+alpha[2]*fUp[4]+fUp[2]*alpha[4]+alpha[0]*fUp[1]+fUp[0]*alpha[1]); 
  incr[2] = 0.25*(alpha[5]*fUp[7]+fUp[5]*alpha[7]+alpha[3]*fUp[6]+fUp[3]*alpha[6]+alpha[1]*fUp[4]+fUp[1]*alpha[4]+alpha[0]*fUp[2]+fUp[0]*alpha[2]); 
  incr[3] = 0.25*(alpha[4]*fUp[7]+fUp[4]*alpha[7]+alpha[2]*fUp[6]+fUp[2]*alpha[6]+alpha[1]*fUp[5]+fUp[1]*alpha[5]+alpha[0]*fUp[3]+fUp[0]*alpha[3]); 
  incr[4] = -0.4330127018922193*(alpha[7]*fUp[7]+alpha[6]*fUp[6]+alpha[5]*fUp[5]+alpha[4]*fUp[4]+alpha[3]*fUp[3]+alpha[2]*fUp[2]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[5] = 0.25*(alpha[3]*fUp[7]+fUp[3]*alpha[7]+alpha[5]*fUp[6]+fUp[5]*alpha[6]+alpha[0]*fUp[4]+fUp[0]*alpha[4]+alpha[1]*fUp[2]+fUp[1]*alpha[2]); 
  incr[6] = 0.25*(alpha[2]*fUp[7]+fUp[2]*alpha[7]+alpha[4]*fUp[6]+fUp[4]*alpha[6]+alpha[0]*fUp[5]+fUp[0]*alpha[5]+alpha[1]*fUp[3]+fUp[1]*alpha[3]); 
  incr[7] = 0.25*(alpha[1]*fUp[7]+fUp[1]*alpha[7]+alpha[0]*fUp[6]+fUp[0]*alpha[6]+alpha[4]*fUp[5]+fUp[4]*alpha[5]+alpha[2]*fUp[3]+fUp[2]*alpha[3]); 
  incr[8] = -0.4330127018922193*(alpha[6]*fUp[7]+fUp[6]*alpha[7]+alpha[3]*fUp[5]+fUp[3]*alpha[5]+alpha[2]*fUp[4]+fUp[2]*alpha[4]+alpha[0]*fUp[1]+fUp[0]*alpha[1]); 
  incr[9] = -0.4330127018922193*(alpha[5]*fUp[7]+fUp[5]*alpha[7]+alpha[3]*fUp[6]+fUp[3]*alpha[6]+alpha[1]*fUp[4]+fUp[1]*alpha[4]+alpha[0]*fUp[2]+fUp[0]*alpha[2]); 
  incr[10] = -0.4330127018922193*(alpha[4]*fUp[7]+fUp[4]*alpha[7]+alpha[2]*fUp[6]+fUp[2]*alpha[6]+alpha[1]*fUp[5]+fUp[1]*alpha[5]+alpha[0]*fUp[3]+fUp[0]*alpha[3]); 
  incr[11] = 0.25*(alpha[0]*fUp[7]+fUp[0]*alpha[7]+alpha[1]*fUp[6]+fUp[1]*alpha[6]+alpha[2]*fUp[5]+fUp[2]*alpha[5]+alpha[3]*fUp[4]+fUp[3]*alpha[4]); 
  incr[12] = -0.4330127018922193*(alpha[3]*fUp[7]+fUp[3]*alpha[7]+alpha[5]*fUp[6]+fUp[5]*alpha[6]+alpha[0]*fUp[4]+fUp[0]*alpha[4]+alpha[1]*fUp[2]+fUp[1]*alpha[2]); 
  incr[13] = -0.4330127018922193*(alpha[2]*fUp[7]+fUp[2]*alpha[7]+alpha[4]*fUp[6]+fUp[4]*alpha[6]+alpha[0]*fUp[5]+fUp[0]*alpha[5]+alpha[1]*fUp[3]+fUp[1]*alpha[3]); 
  incr[14] = -0.4330127018922193*(alpha[1]*fUp[7]+fUp[1]*alpha[7]+alpha[0]*fUp[6]+fUp[0]*alpha[6]+alpha[4]*fUp[5]+fUp[4]*alpha[5]+alpha[2]*fUp[3]+fUp[2]*alpha[3]); 
  incr[15] = -0.4330127018922193*(alpha[0]*fUp[7]+fUp[0]*alpha[7]+alpha[1]*fUp[6]+fUp[1]*alpha[6]+alpha[2]*fUp[5]+fUp[2]*alpha[5]+alpha[3]*fUp[4]+fUp[3]*alpha[4]); 

  outr[0] += incr[0]*dv11r; 
  outr[1] += incr[1]*dv11r; 
  outr[2] += incr[2]*dv11r; 
  outr[3] += incr[3]*dv11r; 
  outr[4] += incr[4]*dv11r; 
  outr[5] += incr[5]*dv11r; 
  outr[6] += incr[6]*dv11r; 
  outr[7] += incr[7]*dv11r; 
  outr[8] += incr[8]*dv11r; 
  outr[9] += incr[9]*dv11r; 
  outr[10] += incr[10]*dv11r; 
  outr[11] += incr[11]*dv11r; 
  outr[12] += incr[12]*dv11r; 
  outr[13] += incr[13]*dv11r; 
  outr[14] += incr[14]*dv11r; 
  outr[15] += incr[15]*dv11r; 

  outl[0] += -1.0*incr[0]*dv11l; 
  outl[1] += -1.0*incr[1]*dv11l; 
  outl[2] += -1.0*incr[2]*dv11l; 
  outl[3] += -1.0*incr[3]*dv11l; 
  outl[4] += incr[4]*dv11l; 
  outl[5] += -1.0*incr[5]*dv11l; 
  outl[6] += -1.0*incr[6]*dv11l; 
  outl[7] += -1.0*incr[7]*dv11l; 
  outl[8] += incr[8]*dv11l; 
  outl[9] += incr[9]*dv11l; 
  outl[10] += incr[10]*dv11l; 
  outl[11] += -1.0*incr[11]*dv11l; 
  outl[12] += incr[12]*dv11l; 
  outl[13] += incr[13]*dv11l; 
  outl[14] += incr[14]*dv11l; 
  outl[15] += incr[15]*dv11l; 


  return std::abs(amid); 
} 
