#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream2x3vSer_X_P1(const int stride_f, const double* __restrict__ wl, const double* __restrict__ wr, const double* __restrict__ dxvl, const double* __restrict__ dxvr, const double* __restrict__ fl, const double* __restrict__ fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[2]; 
  double wxl = wl[2]; 

  double dvxr = dxvr[2]; 
  double wxr = wr[2]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[32]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[stride_f*1]+fl[stride_f*0])*wxl+(0.5*fl[stride_f*7]+0.2886751345948129*fl[stride_f*3])*dvxl; 
  incr[1] = ((-3.0*fl[stride_f*1])-1.732050807568877*fl[stride_f*0])*wxl+((-0.8660254037844386*fl[stride_f*7])-0.5*fl[stride_f*3])*dvxl; 
  incr[2] = (1.732050807568877*fl[stride_f*6]+fl[stride_f*2])*wxl+(0.5*fl[stride_f*16]+0.2886751345948129*fl[stride_f*8])*dvxl; 
  incr[3] = (1.732050807568877*fl[stride_f*7]+fl[stride_f*3])*wxl+(0.5*fl[stride_f*1]+0.2886751345948129*fl[stride_f*0])*dvxl; 
  incr[4] = (1.732050807568877*fl[stride_f*9]+fl[stride_f*4])*wxl+(0.5*fl[stride_f*18]+0.2886751345948129*fl[stride_f*11])*dvxl; 
  incr[5] = (1.732050807568877*fl[stride_f*12]+fl[stride_f*5])*wxl+(0.5*fl[stride_f*21]+0.2886751345948129*fl[stride_f*14])*dvxl; 
  incr[6] = ((-3.0*fl[stride_f*6])-1.732050807568877*fl[stride_f*2])*wxl+((-0.8660254037844386*fl[stride_f*16])-0.5*fl[stride_f*8])*dvxl; 
  incr[7] = ((-3.0*fl[stride_f*7])-1.732050807568877*fl[stride_f*3])*wxl+((-0.8660254037844386*fl[stride_f*1])-0.5*fl[stride_f*0])*dvxl; 
  incr[8] = (1.732050807568877*fl[stride_f*16]+fl[stride_f*8])*wxl+(0.5*fl[stride_f*6]+0.2886751345948129*fl[stride_f*2])*dvxl; 
  incr[9] = ((-3.0*fl[stride_f*9])-1.732050807568877*fl[stride_f*4])*wxl+((-0.8660254037844386*fl[stride_f*18])-0.5*fl[stride_f*11])*dvxl; 
  incr[10] = (1.732050807568877*fl[stride_f*17]+fl[stride_f*10])*wxl+(0.5*fl[stride_f*26]+0.2886751345948129*fl[stride_f*19])*dvxl; 
  incr[11] = (1.732050807568877*fl[stride_f*18]+fl[stride_f*11])*wxl+(0.5*fl[stride_f*9]+0.2886751345948129*fl[stride_f*4])*dvxl; 
  incr[12] = ((-3.0*fl[stride_f*12])-1.732050807568877*fl[stride_f*5])*wxl+((-0.8660254037844386*fl[stride_f*21])-0.5*fl[stride_f*14])*dvxl; 
  incr[13] = (1.732050807568877*fl[stride_f*20]+fl[stride_f*13])*wxl+(0.5*fl[stride_f*27]+0.2886751345948129*fl[stride_f*22])*dvxl; 
  incr[14] = (1.732050807568877*fl[stride_f*21]+fl[stride_f*14])*wxl+(0.5*fl[stride_f*12]+0.2886751345948129*fl[stride_f*5])*dvxl; 
  incr[15] = (1.732050807568877*fl[stride_f*23]+fl[stride_f*15])*wxl+(0.5*fl[stride_f*29]+0.2886751345948129*fl[stride_f*25])*dvxl; 
  incr[16] = ((-3.0*fl[stride_f*16])-1.732050807568877*fl[stride_f*8])*wxl+((-0.8660254037844386*fl[stride_f*6])-0.5*fl[stride_f*2])*dvxl; 
  incr[17] = ((-3.0*fl[stride_f*17])-1.732050807568877*fl[stride_f*10])*wxl+((-0.8660254037844386*fl[stride_f*26])-0.5*fl[stride_f*19])*dvxl; 
  incr[18] = ((-3.0*fl[stride_f*18])-1.732050807568877*fl[stride_f*11])*wxl+((-0.8660254037844386*fl[stride_f*9])-0.5*fl[stride_f*4])*dvxl; 
  incr[19] = (1.732050807568877*fl[stride_f*26]+fl[stride_f*19])*wxl+(0.5*fl[stride_f*17]+0.2886751345948129*fl[stride_f*10])*dvxl; 
  incr[20] = ((-3.0*fl[stride_f*20])-1.732050807568877*fl[stride_f*13])*wxl+((-0.8660254037844386*fl[stride_f*27])-0.5*fl[stride_f*22])*dvxl; 
  incr[21] = ((-3.0*fl[stride_f*21])-1.732050807568877*fl[stride_f*14])*wxl+((-0.8660254037844386*fl[stride_f*12])-0.5*fl[stride_f*5])*dvxl; 
  incr[22] = (1.732050807568877*fl[stride_f*27]+fl[stride_f*22])*wxl+(0.5*fl[stride_f*20]+0.2886751345948129*fl[stride_f*13])*dvxl; 
  incr[23] = ((-3.0*fl[stride_f*23])-1.732050807568877*fl[stride_f*15])*wxl+((-0.8660254037844386*fl[stride_f*29])-0.5*fl[stride_f*25])*dvxl; 
  incr[24] = (1.732050807568877*fl[stride_f*28]+fl[stride_f*24])*wxl+(0.5*fl[stride_f*31]+0.2886751345948129*fl[stride_f*30])*dvxl; 
  incr[25] = (1.732050807568877*fl[stride_f*29]+fl[stride_f*25])*wxl+(0.5*fl[stride_f*23]+0.2886751345948129*fl[stride_f*15])*dvxl; 
  incr[26] = ((-3.0*fl[stride_f*26])-1.732050807568877*fl[stride_f*19])*wxl+((-0.8660254037844386*fl[stride_f*17])-0.5*fl[stride_f*10])*dvxl; 
  incr[27] = ((-3.0*fl[stride_f*27])-1.732050807568877*fl[stride_f*22])*wxl+((-0.8660254037844386*fl[stride_f*20])-0.5*fl[stride_f*13])*dvxl; 
  incr[28] = ((-3.0*fl[stride_f*28])-1.732050807568877*fl[stride_f*24])*wxl+((-0.8660254037844386*fl[stride_f*31])-0.5*fl[stride_f*30])*dvxl; 
  incr[29] = ((-3.0*fl[stride_f*29])-1.732050807568877*fl[stride_f*25])*wxl+((-0.8660254037844386*fl[stride_f*23])-0.5*fl[stride_f*15])*dvxl; 
  incr[30] = (1.732050807568877*fl[stride_f*31]+fl[stride_f*30])*wxl+(0.5*fl[stride_f*28]+0.2886751345948129*fl[stride_f*24])*dvxl; 
  incr[31] = ((-3.0*fl[stride_f*31])-1.732050807568877*fl[stride_f*30])*wxl+((-0.8660254037844386*fl[stride_f*28])-0.5*fl[stride_f*24])*dvxl; 

  if(outr) {
    outr[stride_f*0] += incr[0]*dxr; 
    outr[stride_f*1] += incr[1]*dxr; 
    outr[stride_f*2] += incr[2]*dxr; 
    outr[stride_f*3] += incr[3]*dxr; 
    outr[stride_f*4] += incr[4]*dxr; 
    outr[stride_f*5] += incr[5]*dxr; 
    outr[stride_f*6] += incr[6]*dxr; 
    outr[stride_f*7] += incr[7]*dxr; 
    outr[stride_f*8] += incr[8]*dxr; 
    outr[stride_f*9] += incr[9]*dxr; 
    outr[stride_f*10] += incr[10]*dxr; 
    outr[stride_f*11] += incr[11]*dxr; 
    outr[stride_f*12] += incr[12]*dxr; 
    outr[stride_f*13] += incr[13]*dxr; 
    outr[stride_f*14] += incr[14]*dxr; 
    outr[stride_f*15] += incr[15]*dxr; 
    outr[stride_f*16] += incr[16]*dxr; 
    outr[stride_f*17] += incr[17]*dxr; 
    outr[stride_f*18] += incr[18]*dxr; 
    outr[stride_f*19] += incr[19]*dxr; 
    outr[stride_f*20] += incr[20]*dxr; 
    outr[stride_f*21] += incr[21]*dxr; 
    outr[stride_f*22] += incr[22]*dxr; 
    outr[stride_f*23] += incr[23]*dxr; 
    outr[stride_f*24] += incr[24]*dxr; 
    outr[stride_f*25] += incr[25]*dxr; 
    outr[stride_f*26] += incr[26]*dxr; 
    outr[stride_f*27] += incr[27]*dxr; 
    outr[stride_f*28] += incr[28]*dxr; 
    outr[stride_f*29] += incr[29]*dxr; 
    outr[stride_f*30] += incr[30]*dxr; 
    outr[stride_f*31] += incr[31]*dxr; 
  }
  
  if(outl) {
    outl[stride_f*0] += -1.0*incr[0]*dxl; 
    outl[stride_f*1] += incr[1]*dxl; 
    outl[stride_f*2] += -1.0*incr[2]*dxl; 
    outl[stride_f*3] += -1.0*incr[3]*dxl; 
    outl[stride_f*4] += -1.0*incr[4]*dxl; 
    outl[stride_f*5] += -1.0*incr[5]*dxl; 
    outl[stride_f*6] += incr[6]*dxl; 
    outl[stride_f*7] += incr[7]*dxl; 
    outl[stride_f*8] += -1.0*incr[8]*dxl; 
    outl[stride_f*9] += incr[9]*dxl; 
    outl[stride_f*10] += -1.0*incr[10]*dxl; 
    outl[stride_f*11] += -1.0*incr[11]*dxl; 
    outl[stride_f*12] += incr[12]*dxl; 
    outl[stride_f*13] += -1.0*incr[13]*dxl; 
    outl[stride_f*14] += -1.0*incr[14]*dxl; 
    outl[stride_f*15] += -1.0*incr[15]*dxl; 
    outl[stride_f*16] += incr[16]*dxl; 
    outl[stride_f*17] += incr[17]*dxl; 
    outl[stride_f*18] += incr[18]*dxl; 
    outl[stride_f*19] += -1.0*incr[19]*dxl; 
    outl[stride_f*20] += incr[20]*dxl; 
    outl[stride_f*21] += incr[21]*dxl; 
    outl[stride_f*22] += -1.0*incr[22]*dxl; 
    outl[stride_f*23] += incr[23]*dxl; 
    outl[stride_f*24] += -1.0*incr[24]*dxl; 
    outl[stride_f*25] += -1.0*incr[25]*dxl; 
    outl[stride_f*26] += incr[26]*dxl; 
    outl[stride_f*27] += incr[27]*dxl; 
    outl[stride_f*28] += incr[28]*dxl; 
    outl[stride_f*29] += incr[29]*dxl; 
    outl[stride_f*30] += -1.0*incr[30]*dxl; 
    outl[stride_f*31] += incr[31]*dxl; 
  }
  } else { 
  incr[0] = (fr[stride_f*0]-1.732050807568877*fr[stride_f*1])*wxr+(0.2886751345948129*fr[stride_f*3]-0.5*fr[stride_f*7])*dvxr; 
  incr[1] = (3.0*fr[stride_f*1]-1.732050807568877*fr[stride_f*0])*wxr+(0.8660254037844386*fr[stride_f*7]-0.5*fr[stride_f*3])*dvxr; 
  incr[2] = (fr[stride_f*2]-1.732050807568877*fr[stride_f*6])*wxr+(0.2886751345948129*fr[stride_f*8]-0.5*fr[stride_f*16])*dvxr; 
  incr[3] = (fr[stride_f*3]-1.732050807568877*fr[stride_f*7])*wxr+(0.2886751345948129*fr[stride_f*0]-0.5*fr[stride_f*1])*dvxr; 
  incr[4] = (fr[stride_f*4]-1.732050807568877*fr[stride_f*9])*wxr+(0.2886751345948129*fr[stride_f*11]-0.5*fr[stride_f*18])*dvxr; 
  incr[5] = (fr[stride_f*5]-1.732050807568877*fr[stride_f*12])*wxr+(0.2886751345948129*fr[stride_f*14]-0.5*fr[stride_f*21])*dvxr; 
  incr[6] = (3.0*fr[stride_f*6]-1.732050807568877*fr[stride_f*2])*wxr+(0.8660254037844386*fr[stride_f*16]-0.5*fr[stride_f*8])*dvxr; 
  incr[7] = (3.0*fr[stride_f*7]-1.732050807568877*fr[stride_f*3])*wxr+(0.8660254037844386*fr[stride_f*1]-0.5*fr[stride_f*0])*dvxr; 
  incr[8] = (fr[stride_f*8]-1.732050807568877*fr[stride_f*16])*wxr+(0.2886751345948129*fr[stride_f*2]-0.5*fr[stride_f*6])*dvxr; 
  incr[9] = (3.0*fr[stride_f*9]-1.732050807568877*fr[stride_f*4])*wxr+(0.8660254037844386*fr[stride_f*18]-0.5*fr[stride_f*11])*dvxr; 
  incr[10] = (fr[stride_f*10]-1.732050807568877*fr[stride_f*17])*wxr+(0.2886751345948129*fr[stride_f*19]-0.5*fr[stride_f*26])*dvxr; 
  incr[11] = (fr[stride_f*11]-1.732050807568877*fr[stride_f*18])*wxr+(0.2886751345948129*fr[stride_f*4]-0.5*fr[stride_f*9])*dvxr; 
  incr[12] = (3.0*fr[stride_f*12]-1.732050807568877*fr[stride_f*5])*wxr+(0.8660254037844386*fr[stride_f*21]-0.5*fr[stride_f*14])*dvxr; 
  incr[13] = (fr[stride_f*13]-1.732050807568877*fr[stride_f*20])*wxr+(0.2886751345948129*fr[stride_f*22]-0.5*fr[stride_f*27])*dvxr; 
  incr[14] = (fr[stride_f*14]-1.732050807568877*fr[stride_f*21])*wxr+(0.2886751345948129*fr[stride_f*5]-0.5*fr[stride_f*12])*dvxr; 
  incr[15] = (fr[stride_f*15]-1.732050807568877*fr[stride_f*23])*wxr+(0.2886751345948129*fr[stride_f*25]-0.5*fr[stride_f*29])*dvxr; 
  incr[16] = (3.0*fr[stride_f*16]-1.732050807568877*fr[stride_f*8])*wxr+(0.8660254037844386*fr[stride_f*6]-0.5*fr[stride_f*2])*dvxr; 
  incr[17] = (3.0*fr[stride_f*17]-1.732050807568877*fr[stride_f*10])*wxr+(0.8660254037844386*fr[stride_f*26]-0.5*fr[stride_f*19])*dvxr; 
  incr[18] = (3.0*fr[stride_f*18]-1.732050807568877*fr[stride_f*11])*wxr+(0.8660254037844386*fr[stride_f*9]-0.5*fr[stride_f*4])*dvxr; 
  incr[19] = (fr[stride_f*19]-1.732050807568877*fr[stride_f*26])*wxr+(0.2886751345948129*fr[stride_f*10]-0.5*fr[stride_f*17])*dvxr; 
  incr[20] = (3.0*fr[stride_f*20]-1.732050807568877*fr[stride_f*13])*wxr+(0.8660254037844386*fr[stride_f*27]-0.5*fr[stride_f*22])*dvxr; 
  incr[21] = (3.0*fr[stride_f*21]-1.732050807568877*fr[stride_f*14])*wxr+(0.8660254037844386*fr[stride_f*12]-0.5*fr[stride_f*5])*dvxr; 
  incr[22] = (fr[stride_f*22]-1.732050807568877*fr[stride_f*27])*wxr+(0.2886751345948129*fr[stride_f*13]-0.5*fr[stride_f*20])*dvxr; 
  incr[23] = (3.0*fr[stride_f*23]-1.732050807568877*fr[stride_f*15])*wxr+(0.8660254037844386*fr[stride_f*29]-0.5*fr[stride_f*25])*dvxr; 
  incr[24] = (fr[stride_f*24]-1.732050807568877*fr[stride_f*28])*wxr+(0.2886751345948129*fr[stride_f*30]-0.5*fr[stride_f*31])*dvxr; 
  incr[25] = (fr[stride_f*25]-1.732050807568877*fr[stride_f*29])*wxr+(0.2886751345948129*fr[stride_f*15]-0.5*fr[stride_f*23])*dvxr; 
  incr[26] = (3.0*fr[stride_f*26]-1.732050807568877*fr[stride_f*19])*wxr+(0.8660254037844386*fr[stride_f*17]-0.5*fr[stride_f*10])*dvxr; 
  incr[27] = (3.0*fr[stride_f*27]-1.732050807568877*fr[stride_f*22])*wxr+(0.8660254037844386*fr[stride_f*20]-0.5*fr[stride_f*13])*dvxr; 
  incr[28] = (3.0*fr[stride_f*28]-1.732050807568877*fr[stride_f*24])*wxr+(0.8660254037844386*fr[stride_f*31]-0.5*fr[stride_f*30])*dvxr; 
  incr[29] = (3.0*fr[stride_f*29]-1.732050807568877*fr[stride_f*25])*wxr+(0.8660254037844386*fr[stride_f*23]-0.5*fr[stride_f*15])*dvxr; 
  incr[30] = (fr[stride_f*30]-1.732050807568877*fr[stride_f*31])*wxr+(0.2886751345948129*fr[stride_f*24]-0.5*fr[stride_f*28])*dvxr; 
  incr[31] = (3.0*fr[stride_f*31]-1.732050807568877*fr[stride_f*30])*wxr+(0.8660254037844386*fr[stride_f*28]-0.5*fr[stride_f*24])*dvxr; 

  if(outr) {
    outr[stride_f*0] += incr[0]*dxr; 
    outr[stride_f*1] += incr[1]*dxr; 
    outr[stride_f*2] += incr[2]*dxr; 
    outr[stride_f*3] += incr[3]*dxr; 
    outr[stride_f*4] += incr[4]*dxr; 
    outr[stride_f*5] += incr[5]*dxr; 
    outr[stride_f*6] += incr[6]*dxr; 
    outr[stride_f*7] += incr[7]*dxr; 
    outr[stride_f*8] += incr[8]*dxr; 
    outr[stride_f*9] += incr[9]*dxr; 
    outr[stride_f*10] += incr[10]*dxr; 
    outr[stride_f*11] += incr[11]*dxr; 
    outr[stride_f*12] += incr[12]*dxr; 
    outr[stride_f*13] += incr[13]*dxr; 
    outr[stride_f*14] += incr[14]*dxr; 
    outr[stride_f*15] += incr[15]*dxr; 
    outr[stride_f*16] += incr[16]*dxr; 
    outr[stride_f*17] += incr[17]*dxr; 
    outr[stride_f*18] += incr[18]*dxr; 
    outr[stride_f*19] += incr[19]*dxr; 
    outr[stride_f*20] += incr[20]*dxr; 
    outr[stride_f*21] += incr[21]*dxr; 
    outr[stride_f*22] += incr[22]*dxr; 
    outr[stride_f*23] += incr[23]*dxr; 
    outr[stride_f*24] += incr[24]*dxr; 
    outr[stride_f*25] += incr[25]*dxr; 
    outr[stride_f*26] += incr[26]*dxr; 
    outr[stride_f*27] += incr[27]*dxr; 
    outr[stride_f*28] += incr[28]*dxr; 
    outr[stride_f*29] += incr[29]*dxr; 
    outr[stride_f*30] += incr[30]*dxr; 
    outr[stride_f*31] += incr[31]*dxr; 
  }
  
  if(outl) {
    outl[stride_f*0] += -1.0*incr[0]*dxl; 
    outl[stride_f*1] += incr[1]*dxl; 
    outl[stride_f*2] += -1.0*incr[2]*dxl; 
    outl[stride_f*3] += -1.0*incr[3]*dxl; 
    outl[stride_f*4] += -1.0*incr[4]*dxl; 
    outl[stride_f*5] += -1.0*incr[5]*dxl; 
    outl[stride_f*6] += incr[6]*dxl; 
    outl[stride_f*7] += incr[7]*dxl; 
    outl[stride_f*8] += -1.0*incr[8]*dxl; 
    outl[stride_f*9] += incr[9]*dxl; 
    outl[stride_f*10] += -1.0*incr[10]*dxl; 
    outl[stride_f*11] += -1.0*incr[11]*dxl; 
    outl[stride_f*12] += incr[12]*dxl; 
    outl[stride_f*13] += -1.0*incr[13]*dxl; 
    outl[stride_f*14] += -1.0*incr[14]*dxl; 
    outl[stride_f*15] += -1.0*incr[15]*dxl; 
    outl[stride_f*16] += incr[16]*dxl; 
    outl[stride_f*17] += incr[17]*dxl; 
    outl[stride_f*18] += incr[18]*dxl; 
    outl[stride_f*19] += -1.0*incr[19]*dxl; 
    outl[stride_f*20] += incr[20]*dxl; 
    outl[stride_f*21] += incr[21]*dxl; 
    outl[stride_f*22] += -1.0*incr[22]*dxl; 
    outl[stride_f*23] += incr[23]*dxl; 
    outl[stride_f*24] += -1.0*incr[24]*dxl; 
    outl[stride_f*25] += -1.0*incr[25]*dxl; 
    outl[stride_f*26] += incr[26]*dxl; 
    outl[stride_f*27] += incr[27]*dxl; 
    outl[stride_f*28] += incr[28]*dxl; 
    outl[stride_f*29] += incr[29]*dxl; 
    outl[stride_f*30] += -1.0*incr[30]*dxl; 
    outl[stride_f*31] += incr[31]*dxl; 
  }
  } 
} 
__host__ __device__ void VlasovSurfStream2x3vSer_Y_P1(const int stride_f, const double* __restrict__ wl, const double* __restrict__ wr, const double* __restrict__ dxvl, const double* __restrict__ dxvr, const double* __restrict__ fl, const double* __restrict__ fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[3]; 
  double wxl = wl[3]; 

  double dvxr = dxvr[3]; 
  double wxr = wr[3]; 

  double dxl = 1.0/dxvl[1]; 
  double dxr = 1.0/dxvr[1]; 

  double incr[32]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[stride_f*2]+fl[stride_f*0])*wxl+(0.5*fl[stride_f*10]+0.2886751345948129*fl[stride_f*4])*dvxl; 
  incr[1] = (1.732050807568877*fl[stride_f*6]+fl[stride_f*1])*wxl+(0.5*fl[stride_f*17]+0.2886751345948129*fl[stride_f*9])*dvxl; 
  incr[2] = ((-3.0*fl[stride_f*2])-1.732050807568877*fl[stride_f*0])*wxl+((-0.8660254037844386*fl[stride_f*10])-0.5*fl[stride_f*4])*dvxl; 
  incr[3] = (1.732050807568877*fl[stride_f*8]+fl[stride_f*3])*wxl+(0.5*fl[stride_f*19]+0.2886751345948129*fl[stride_f*11])*dvxl; 
  incr[4] = (1.732050807568877*fl[stride_f*10]+fl[stride_f*4])*wxl+(0.5*fl[stride_f*2]+0.2886751345948129*fl[stride_f*0])*dvxl; 
  incr[5] = (1.732050807568877*fl[stride_f*13]+fl[stride_f*5])*wxl+(0.5*fl[stride_f*24]+0.2886751345948129*fl[stride_f*15])*dvxl; 
  incr[6] = ((-3.0*fl[stride_f*6])-1.732050807568877*fl[stride_f*1])*wxl+((-0.8660254037844386*fl[stride_f*17])-0.5*fl[stride_f*9])*dvxl; 
  incr[7] = (1.732050807568877*fl[stride_f*16]+fl[stride_f*7])*wxl+(0.5*fl[stride_f*26]+0.2886751345948129*fl[stride_f*18])*dvxl; 
  incr[8] = ((-3.0*fl[stride_f*8])-1.732050807568877*fl[stride_f*3])*wxl+((-0.8660254037844386*fl[stride_f*19])-0.5*fl[stride_f*11])*dvxl; 
  incr[9] = (1.732050807568877*fl[stride_f*17]+fl[stride_f*9])*wxl+(0.5*fl[stride_f*6]+0.2886751345948129*fl[stride_f*1])*dvxl; 
  incr[10] = ((-3.0*fl[stride_f*10])-1.732050807568877*fl[stride_f*4])*wxl+((-0.8660254037844386*fl[stride_f*2])-0.5*fl[stride_f*0])*dvxl; 
  incr[11] = (1.732050807568877*fl[stride_f*19]+fl[stride_f*11])*wxl+(0.5*fl[stride_f*8]+0.2886751345948129*fl[stride_f*3])*dvxl; 
  incr[12] = (1.732050807568877*fl[stride_f*20]+fl[stride_f*12])*wxl+(0.5*fl[stride_f*28]+0.2886751345948129*fl[stride_f*23])*dvxl; 
  incr[13] = ((-3.0*fl[stride_f*13])-1.732050807568877*fl[stride_f*5])*wxl+((-0.8660254037844386*fl[stride_f*24])-0.5*fl[stride_f*15])*dvxl; 
  incr[14] = (1.732050807568877*fl[stride_f*22]+fl[stride_f*14])*wxl+(0.5*fl[stride_f*30]+0.2886751345948129*fl[stride_f*25])*dvxl; 
  incr[15] = (1.732050807568877*fl[stride_f*24]+fl[stride_f*15])*wxl+(0.5*fl[stride_f*13]+0.2886751345948129*fl[stride_f*5])*dvxl; 
  incr[16] = ((-3.0*fl[stride_f*16])-1.732050807568877*fl[stride_f*7])*wxl+((-0.8660254037844386*fl[stride_f*26])-0.5*fl[stride_f*18])*dvxl; 
  incr[17] = ((-3.0*fl[stride_f*17])-1.732050807568877*fl[stride_f*9])*wxl+((-0.8660254037844386*fl[stride_f*6])-0.5*fl[stride_f*1])*dvxl; 
  incr[18] = (1.732050807568877*fl[stride_f*26]+fl[stride_f*18])*wxl+(0.5*fl[stride_f*16]+0.2886751345948129*fl[stride_f*7])*dvxl; 
  incr[19] = ((-3.0*fl[stride_f*19])-1.732050807568877*fl[stride_f*11])*wxl+((-0.8660254037844386*fl[stride_f*8])-0.5*fl[stride_f*3])*dvxl; 
  incr[20] = ((-3.0*fl[stride_f*20])-1.732050807568877*fl[stride_f*12])*wxl+((-0.8660254037844386*fl[stride_f*28])-0.5*fl[stride_f*23])*dvxl; 
  incr[21] = (1.732050807568877*fl[stride_f*27]+fl[stride_f*21])*wxl+(0.5*fl[stride_f*31]+0.2886751345948129*fl[stride_f*29])*dvxl; 
  incr[22] = ((-3.0*fl[stride_f*22])-1.732050807568877*fl[stride_f*14])*wxl+((-0.8660254037844386*fl[stride_f*30])-0.5*fl[stride_f*25])*dvxl; 
  incr[23] = (1.732050807568877*fl[stride_f*28]+fl[stride_f*23])*wxl+(0.5*fl[stride_f*20]+0.2886751345948129*fl[stride_f*12])*dvxl; 
  incr[24] = ((-3.0*fl[stride_f*24])-1.732050807568877*fl[stride_f*15])*wxl+((-0.8660254037844386*fl[stride_f*13])-0.5*fl[stride_f*5])*dvxl; 
  incr[25] = (1.732050807568877*fl[stride_f*30]+fl[stride_f*25])*wxl+(0.5*fl[stride_f*22]+0.2886751345948129*fl[stride_f*14])*dvxl; 
  incr[26] = ((-3.0*fl[stride_f*26])-1.732050807568877*fl[stride_f*18])*wxl+((-0.8660254037844386*fl[stride_f*16])-0.5*fl[stride_f*7])*dvxl; 
  incr[27] = ((-3.0*fl[stride_f*27])-1.732050807568877*fl[stride_f*21])*wxl+((-0.8660254037844386*fl[stride_f*31])-0.5*fl[stride_f*29])*dvxl; 
  incr[28] = ((-3.0*fl[stride_f*28])-1.732050807568877*fl[stride_f*23])*wxl+((-0.8660254037844386*fl[stride_f*20])-0.5*fl[stride_f*12])*dvxl; 
  incr[29] = (1.732050807568877*fl[stride_f*31]+fl[stride_f*29])*wxl+(0.5*fl[stride_f*27]+0.2886751345948129*fl[stride_f*21])*dvxl; 
  incr[30] = ((-3.0*fl[stride_f*30])-1.732050807568877*fl[stride_f*25])*wxl+((-0.8660254037844386*fl[stride_f*22])-0.5*fl[stride_f*14])*dvxl; 
  incr[31] = ((-3.0*fl[stride_f*31])-1.732050807568877*fl[stride_f*29])*wxl+((-0.8660254037844386*fl[stride_f*27])-0.5*fl[stride_f*21])*dvxl; 

if(outr) {
  outr[stride_f*0] += incr[0]*dxr; 
  outr[stride_f*1] += incr[1]*dxr; 
  outr[stride_f*2] += incr[2]*dxr; 
  outr[stride_f*3] += incr[3]*dxr; 
  outr[stride_f*4] += incr[4]*dxr; 
  outr[stride_f*5] += incr[5]*dxr; 
  outr[stride_f*6] += incr[6]*dxr; 
  outr[stride_f*7] += incr[7]*dxr; 
  outr[stride_f*8] += incr[8]*dxr; 
  outr[stride_f*9] += incr[9]*dxr; 
  outr[stride_f*10] += incr[10]*dxr; 
  outr[stride_f*11] += incr[11]*dxr; 
  outr[stride_f*12] += incr[12]*dxr; 
  outr[stride_f*13] += incr[13]*dxr; 
  outr[stride_f*14] += incr[14]*dxr; 
  outr[stride_f*15] += incr[15]*dxr; 
  outr[stride_f*16] += incr[16]*dxr; 
  outr[stride_f*17] += incr[17]*dxr; 
  outr[stride_f*18] += incr[18]*dxr; 
  outr[stride_f*19] += incr[19]*dxr; 
  outr[stride_f*20] += incr[20]*dxr; 
  outr[stride_f*21] += incr[21]*dxr; 
  outr[stride_f*22] += incr[22]*dxr; 
  outr[stride_f*23] += incr[23]*dxr; 
  outr[stride_f*24] += incr[24]*dxr; 
  outr[stride_f*25] += incr[25]*dxr; 
  outr[stride_f*26] += incr[26]*dxr; 
  outr[stride_f*27] += incr[27]*dxr; 
  outr[stride_f*28] += incr[28]*dxr; 
  outr[stride_f*29] += incr[29]*dxr; 
  outr[stride_f*30] += incr[30]*dxr; 
  outr[stride_f*31] += incr[31]*dxr; 
}

if(outl) {
  outl[stride_f*0] += -1.0*incr[0]*dxl; 
  outl[stride_f*1] += -1.0*incr[1]*dxl; 
  outl[stride_f*2] += incr[2]*dxl; 
  outl[stride_f*3] += -1.0*incr[3]*dxl; 
  outl[stride_f*4] += -1.0*incr[4]*dxl; 
  outl[stride_f*5] += -1.0*incr[5]*dxl; 
  outl[stride_f*6] += incr[6]*dxl; 
  outl[stride_f*7] += -1.0*incr[7]*dxl; 
  outl[stride_f*8] += incr[8]*dxl; 
  outl[stride_f*9] += -1.0*incr[9]*dxl; 
  outl[stride_f*10] += incr[10]*dxl; 
  outl[stride_f*11] += -1.0*incr[11]*dxl; 
  outl[stride_f*12] += -1.0*incr[12]*dxl; 
  outl[stride_f*13] += incr[13]*dxl; 
  outl[stride_f*14] += -1.0*incr[14]*dxl; 
  outl[stride_f*15] += -1.0*incr[15]*dxl; 
  outl[stride_f*16] += incr[16]*dxl; 
  outl[stride_f*17] += incr[17]*dxl; 
  outl[stride_f*18] += -1.0*incr[18]*dxl; 
  outl[stride_f*19] += incr[19]*dxl; 
  outl[stride_f*20] += incr[20]*dxl; 
  outl[stride_f*21] += -1.0*incr[21]*dxl; 
  outl[stride_f*22] += incr[22]*dxl; 
  outl[stride_f*23] += -1.0*incr[23]*dxl; 
  outl[stride_f*24] += incr[24]*dxl; 
  outl[stride_f*25] += -1.0*incr[25]*dxl; 
  outl[stride_f*26] += incr[26]*dxl; 
  outl[stride_f*27] += incr[27]*dxl; 
  outl[stride_f*28] += incr[28]*dxl; 
  outl[stride_f*29] += -1.0*incr[29]*dxl; 
  outl[stride_f*30] += incr[30]*dxl; 
  outl[stride_f*31] += incr[31]*dxl; 
}
  } else { 
  incr[0] = (fr[stride_f*0]-1.732050807568877*fr[stride_f*2])*wxr+(0.2886751345948129*fr[stride_f*4]-0.5*fr[stride_f*10])*dvxr; 
  incr[1] = (fr[stride_f*1]-1.732050807568877*fr[stride_f*6])*wxr+(0.2886751345948129*fr[stride_f*9]-0.5*fr[stride_f*17])*dvxr; 
  incr[2] = (3.0*fr[stride_f*2]-1.732050807568877*fr[stride_f*0])*wxr+(0.8660254037844386*fr[stride_f*10]-0.5*fr[stride_f*4])*dvxr; 
  incr[3] = (fr[stride_f*3]-1.732050807568877*fr[stride_f*8])*wxr+(0.2886751345948129*fr[stride_f*11]-0.5*fr[stride_f*19])*dvxr; 
  incr[4] = (fr[stride_f*4]-1.732050807568877*fr[stride_f*10])*wxr+(0.2886751345948129*fr[stride_f*0]-0.5*fr[stride_f*2])*dvxr; 
  incr[5] = (fr[stride_f*5]-1.732050807568877*fr[stride_f*13])*wxr+(0.2886751345948129*fr[stride_f*15]-0.5*fr[stride_f*24])*dvxr; 
  incr[6] = (3.0*fr[stride_f*6]-1.732050807568877*fr[stride_f*1])*wxr+(0.8660254037844386*fr[stride_f*17]-0.5*fr[stride_f*9])*dvxr; 
  incr[7] = (fr[stride_f*7]-1.732050807568877*fr[stride_f*16])*wxr+(0.2886751345948129*fr[stride_f*18]-0.5*fr[stride_f*26])*dvxr; 
  incr[8] = (3.0*fr[stride_f*8]-1.732050807568877*fr[stride_f*3])*wxr+(0.8660254037844386*fr[stride_f*19]-0.5*fr[stride_f*11])*dvxr; 
  incr[9] = (fr[stride_f*9]-1.732050807568877*fr[stride_f*17])*wxr+(0.2886751345948129*fr[stride_f*1]-0.5*fr[stride_f*6])*dvxr; 
  incr[10] = (3.0*fr[stride_f*10]-1.732050807568877*fr[stride_f*4])*wxr+(0.8660254037844386*fr[stride_f*2]-0.5*fr[stride_f*0])*dvxr; 
  incr[11] = (fr[stride_f*11]-1.732050807568877*fr[stride_f*19])*wxr+(0.2886751345948129*fr[stride_f*3]-0.5*fr[stride_f*8])*dvxr; 
  incr[12] = (fr[stride_f*12]-1.732050807568877*fr[stride_f*20])*wxr+(0.2886751345948129*fr[stride_f*23]-0.5*fr[stride_f*28])*dvxr; 
  incr[13] = (3.0*fr[stride_f*13]-1.732050807568877*fr[stride_f*5])*wxr+(0.8660254037844386*fr[stride_f*24]-0.5*fr[stride_f*15])*dvxr; 
  incr[14] = (fr[stride_f*14]-1.732050807568877*fr[stride_f*22])*wxr+(0.2886751345948129*fr[stride_f*25]-0.5*fr[stride_f*30])*dvxr; 
  incr[15] = (fr[stride_f*15]-1.732050807568877*fr[stride_f*24])*wxr+(0.2886751345948129*fr[stride_f*5]-0.5*fr[stride_f*13])*dvxr; 
  incr[16] = (3.0*fr[stride_f*16]-1.732050807568877*fr[stride_f*7])*wxr+(0.8660254037844386*fr[stride_f*26]-0.5*fr[stride_f*18])*dvxr; 
  incr[17] = (3.0*fr[stride_f*17]-1.732050807568877*fr[stride_f*9])*wxr+(0.8660254037844386*fr[stride_f*6]-0.5*fr[stride_f*1])*dvxr; 
  incr[18] = (fr[stride_f*18]-1.732050807568877*fr[stride_f*26])*wxr+(0.2886751345948129*fr[stride_f*7]-0.5*fr[stride_f*16])*dvxr; 
  incr[19] = (3.0*fr[stride_f*19]-1.732050807568877*fr[stride_f*11])*wxr+(0.8660254037844386*fr[stride_f*8]-0.5*fr[stride_f*3])*dvxr; 
  incr[20] = (3.0*fr[stride_f*20]-1.732050807568877*fr[stride_f*12])*wxr+(0.8660254037844386*fr[stride_f*28]-0.5*fr[stride_f*23])*dvxr; 
  incr[21] = (fr[stride_f*21]-1.732050807568877*fr[stride_f*27])*wxr+(0.2886751345948129*fr[stride_f*29]-0.5*fr[stride_f*31])*dvxr; 
  incr[22] = (3.0*fr[stride_f*22]-1.732050807568877*fr[stride_f*14])*wxr+(0.8660254037844386*fr[stride_f*30]-0.5*fr[stride_f*25])*dvxr; 
  incr[23] = (fr[stride_f*23]-1.732050807568877*fr[stride_f*28])*wxr+(0.2886751345948129*fr[stride_f*12]-0.5*fr[stride_f*20])*dvxr; 
  incr[24] = (3.0*fr[stride_f*24]-1.732050807568877*fr[stride_f*15])*wxr+(0.8660254037844386*fr[stride_f*13]-0.5*fr[stride_f*5])*dvxr; 
  incr[25] = (fr[stride_f*25]-1.732050807568877*fr[stride_f*30])*wxr+(0.2886751345948129*fr[stride_f*14]-0.5*fr[stride_f*22])*dvxr; 
  incr[26] = (3.0*fr[stride_f*26]-1.732050807568877*fr[stride_f*18])*wxr+(0.8660254037844386*fr[stride_f*16]-0.5*fr[stride_f*7])*dvxr; 
  incr[27] = (3.0*fr[stride_f*27]-1.732050807568877*fr[stride_f*21])*wxr+(0.8660254037844386*fr[stride_f*31]-0.5*fr[stride_f*29])*dvxr; 
  incr[28] = (3.0*fr[stride_f*28]-1.732050807568877*fr[stride_f*23])*wxr+(0.8660254037844386*fr[stride_f*20]-0.5*fr[stride_f*12])*dvxr; 
  incr[29] = (fr[stride_f*29]-1.732050807568877*fr[stride_f*31])*wxr+(0.2886751345948129*fr[stride_f*21]-0.5*fr[stride_f*27])*dvxr; 
  incr[30] = (3.0*fr[stride_f*30]-1.732050807568877*fr[stride_f*25])*wxr+(0.8660254037844386*fr[stride_f*22]-0.5*fr[stride_f*14])*dvxr; 
  incr[31] = (3.0*fr[stride_f*31]-1.732050807568877*fr[stride_f*29])*wxr+(0.8660254037844386*fr[stride_f*27]-0.5*fr[stride_f*21])*dvxr; 

if(outr) {
  outr[stride_f*0] += incr[0]*dxr; 
  outr[stride_f*1] += incr[1]*dxr; 
  outr[stride_f*2] += incr[2]*dxr; 
  outr[stride_f*3] += incr[3]*dxr; 
  outr[stride_f*4] += incr[4]*dxr; 
  outr[stride_f*5] += incr[5]*dxr; 
  outr[stride_f*6] += incr[6]*dxr; 
  outr[stride_f*7] += incr[7]*dxr; 
  outr[stride_f*8] += incr[8]*dxr; 
  outr[stride_f*9] += incr[9]*dxr; 
  outr[stride_f*10] += incr[10]*dxr; 
  outr[stride_f*11] += incr[11]*dxr; 
  outr[stride_f*12] += incr[12]*dxr; 
  outr[stride_f*13] += incr[13]*dxr; 
  outr[stride_f*14] += incr[14]*dxr; 
  outr[stride_f*15] += incr[15]*dxr; 
  outr[stride_f*16] += incr[16]*dxr; 
  outr[stride_f*17] += incr[17]*dxr; 
  outr[stride_f*18] += incr[18]*dxr; 
  outr[stride_f*19] += incr[19]*dxr; 
  outr[stride_f*20] += incr[20]*dxr; 
  outr[stride_f*21] += incr[21]*dxr; 
  outr[stride_f*22] += incr[22]*dxr; 
  outr[stride_f*23] += incr[23]*dxr; 
  outr[stride_f*24] += incr[24]*dxr; 
  outr[stride_f*25] += incr[25]*dxr; 
  outr[stride_f*26] += incr[26]*dxr; 
  outr[stride_f*27] += incr[27]*dxr; 
  outr[stride_f*28] += incr[28]*dxr; 
  outr[stride_f*29] += incr[29]*dxr; 
  outr[stride_f*30] += incr[30]*dxr; 
  outr[stride_f*31] += incr[31]*dxr; 
}

if(outl) {
  outl[stride_f*0] += -1.0*incr[0]*dxl; 
  outl[stride_f*1] += -1.0*incr[1]*dxl; 
  outl[stride_f*2] += incr[2]*dxl; 
  outl[stride_f*3] += -1.0*incr[3]*dxl; 
  outl[stride_f*4] += -1.0*incr[4]*dxl; 
  outl[stride_f*5] += -1.0*incr[5]*dxl; 
  outl[stride_f*6] += incr[6]*dxl; 
  outl[stride_f*7] += -1.0*incr[7]*dxl; 
  outl[stride_f*8] += incr[8]*dxl; 
  outl[stride_f*9] += -1.0*incr[9]*dxl; 
  outl[stride_f*10] += incr[10]*dxl; 
  outl[stride_f*11] += -1.0*incr[11]*dxl; 
  outl[stride_f*12] += -1.0*incr[12]*dxl; 
  outl[stride_f*13] += incr[13]*dxl; 
  outl[stride_f*14] += -1.0*incr[14]*dxl; 
  outl[stride_f*15] += -1.0*incr[15]*dxl; 
  outl[stride_f*16] += incr[16]*dxl; 
  outl[stride_f*17] += incr[17]*dxl; 
  outl[stride_f*18] += -1.0*incr[18]*dxl; 
  outl[stride_f*19] += incr[19]*dxl; 
  outl[stride_f*20] += incr[20]*dxl; 
  outl[stride_f*21] += -1.0*incr[21]*dxl; 
  outl[stride_f*22] += incr[22]*dxl; 
  outl[stride_f*23] += -1.0*incr[23]*dxl; 
  outl[stride_f*24] += incr[24]*dxl; 
  outl[stride_f*25] += -1.0*incr[25]*dxl; 
  outl[stride_f*26] += incr[26]*dxl; 
  outl[stride_f*27] += incr[27]*dxl; 
  outl[stride_f*28] += incr[28]*dxl; 
  outl[stride_f*29] += -1.0*incr[29]*dxl; 
  outl[stride_f*30] += incr[30]*dxl; 
  outl[stride_f*31] += incr[31]*dxl; 
}
  } 
} 
