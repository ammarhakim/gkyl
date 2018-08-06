#include <VlasovModDecl.h> 
double VlasovSurfElcMag1x2vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double Ghat[8]; 

  for(unsigned int i=0; i<8; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[8]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  double fjump[8]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(-1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(-1*fr[7]-fl[7]); 
  double alpha[8]; 

  alpha[0] = 2.0*B2[0]*wv2+2.0*E0[0]; 
  alpha[1] = 2.0*B2[1]*wv2+2.0*E0[1]; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[5] = 0.5773502691896258*B2[1]*dv2; 
  const double amid = 0.3535533905932737*alpha[0]; 
  Ghat[0] += 0.3061862178478971*(alpha[5]*favg[7]+alpha[3]*favg[6])+0.1767766952966368*alpha[5]*favg[5]+0.3061862178478971*alpha[1]*favg[4]+0.1767766952966368*alpha[3]*favg[3]-0.8660254037844386*fjump[2]+0.3061862178478971*alpha[0]*favg[2]+0.1767766952966368*alpha[1]*favg[1]-0.5*fjump[0]+0.1767766952966368*alpha[0]*favg[0]; 
  Ghat[1] += 0.3061862178478971*(alpha[3]*favg[7]+alpha[5]*favg[6])+0.1767766952966368*(alpha[3]*favg[5]+favg[3]*alpha[5])-0.8660254037844386*fjump[4]+0.3061862178478971*(alpha[0]*favg[4]+alpha[1]*favg[2])-0.5*fjump[1]+0.1767766952966368*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[3] += 0.3061862178478971*alpha[1]*favg[7]-0.8660254037844386*fjump[6]+0.3061862178478971*alpha[0]*favg[6]+0.1767766952966368*alpha[1]*favg[5]+(0.3061862178478971*favg[4]+0.1767766952966368*favg[1])*alpha[5]-0.5*fjump[3]+0.1767766952966368*alpha[0]*favg[3]+(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[3]; 
  Ghat[5] += (-0.8660254037844386*fjump[7])+0.3061862178478971*(alpha[0]*favg[7]+alpha[1]*favg[6])-0.5*fjump[5]+0.1767766952966368*alpha[0]*favg[5]+(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[5]+0.3061862178478971*alpha[3]*favg[4]+0.1767766952966368*(alpha[1]*favg[3]+favg[1]*alpha[3]); 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[3] += 0.5*Ghat[3]*dv10r; 
  outr[4] += -0.8660254037844386*Ghat[1]*dv10r; 
  outr[5] += 0.5*Ghat[5]*dv10r; 
  outr[6] += -0.8660254037844386*Ghat[3]*dv10r; 
  outr[7] += -0.8660254037844386*Ghat[5]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[3] += -0.5*Ghat[3]*dv10l; 
  outl[4] += -0.8660254037844386*Ghat[1]*dv10l; 
  outl[5] += -0.5*Ghat[5]*dv10l; 
  outl[6] += -0.8660254037844386*Ghat[3]*dv10l; 
  outl[7] += -0.8660254037844386*Ghat[5]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[20]; 

  for(unsigned int i=0; i<20; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[20]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = -1*fr[16]+fl[16]; 
  favg[17] = -1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = -1*fr[19]+fl[19]; 
  double fjump[20]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(-1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(1*fr[9]-fl[9]); 
  fjump[10] = amax*(-1*fr[10]-fl[10]); 
  fjump[11] = amax*(-1*fr[11]-fl[11]); 
  fjump[12] = amax*(1*fr[12]-fl[12]); 
  fjump[13] = amax*(1*fr[13]-fl[13]); 
  fjump[14] = amax*(1*fr[14]-fl[14]); 
  fjump[15] = amax*(1*fr[15]-fl[15]); 
  fjump[16] = amax*(-1*fr[16]-fl[16]); 
  fjump[17] = amax*(-1*fr[17]-fl[17]); 
  fjump[18] = amax*(1*fr[18]-fl[18]); 
  fjump[19] = amax*(-1*fr[19]-fl[19]); 
  double alpha[20]; 

  alpha[0] = 2.0*B2[0]*wv2+2.0*E0[0]; 
  alpha[1] = 2.0*B2[1]*wv2+2.0*E0[1]; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[5] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 2.0*B2[2]*wv2+2.0*E0[2]; 
  alpha[13] = 0.5773502691896257*B2[2]*dv2; 
  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 
  Ghat[0] += 0.3952847075210473*alpha[5]*favg[18]+0.3061862178478971*alpha[13]*favg[17]+0.3952847075210473*alpha[3]*favg[14]+0.1767766952966368*alpha[13]*favg[13]+0.3952847075210473*alpha[1]*favg[12]+0.3061862178478971*(alpha[7]*favg[11]+alpha[5]*favg[10])-1.118033988749895*fjump[8]+0.3952847075210473*alpha[0]*favg[8]+0.1767766952966368*alpha[7]*favg[7]+0.3061862178478971*alpha[3]*favg[6]+0.1767766952966368*alpha[5]*favg[5]+0.3061862178478971*alpha[1]*favg[4]+0.1767766952966368*alpha[3]*favg[3]-0.8660254037844386*fjump[2]+0.3061862178478971*alpha[0]*favg[2]+0.1767766952966368*alpha[1]*favg[1]-0.5*fjump[0]+0.1767766952966368*alpha[0]*favg[0]; 
  Ghat[1] += (0.3535533905932737*alpha[13]+0.3952847075210473*alpha[3])*favg[18]+alpha[5]*(0.273861278752583*favg[17]+0.3952847075210473*favg[14]+0.1581138830084189*favg[13])+(0.273861278752583*favg[10]+0.1581138830084189*favg[5])*alpha[13]-1.118033988749895*fjump[12]+(0.3535533905932737*alpha[7]+0.3952847075210473*alpha[0])*favg[12]+0.273861278752583*alpha[1]*favg[11]+0.3061862178478971*alpha[3]*favg[10]+alpha[1]*(0.3952847075210473*favg[8]+0.1581138830084189*favg[7])+(0.273861278752583*favg[4]+0.1581138830084189*favg[1])*alpha[7]+0.3061862178478971*alpha[5]*favg[6]+0.1767766952966368*(alpha[3]*favg[5]+favg[3]*alpha[5])-0.8660254037844386*fjump[4]+0.3061862178478971*(alpha[0]*favg[4]+alpha[1]*favg[2])-0.5*fjump[1]+0.1767766952966368*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[3] += 0.273861278752583*alpha[5]*favg[19]+0.3952847075210473*alpha[1]*favg[18]+0.3061862178478971*alpha[7]*favg[17]+0.273861278752583*alpha[3]*favg[16]+0.1581138830084189*alpha[5]*favg[15]-1.118033988749895*fjump[14]+0.3952847075210473*alpha[0]*favg[14]+0.1767766952966368*alpha[7]*favg[13]+(0.3061862178478971*favg[11]+0.1767766952966368*favg[7])*alpha[13]+0.3952847075210473*alpha[5]*favg[12]+0.3061862178478971*alpha[1]*favg[10]+alpha[3]*(0.1581138830084189*favg[9]+0.3952847075210473*favg[8])-0.8660254037844386*fjump[6]+0.3061862178478971*alpha[0]*favg[6]+0.1767766952966368*alpha[1]*favg[5]+(0.3061862178478971*favg[4]+0.1767766952966368*favg[1])*alpha[5]-0.5*fjump[3]+0.1767766952966368*alpha[0]*favg[3]+(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[3]; 
  Ghat[5] += (0.2449489742783178*alpha[13]+0.273861278752583*alpha[3])*favg[19]-1.118033988749895*fjump[18]+(0.3535533905932737*alpha[7]+0.3952847075210473*alpha[0])*favg[18]+0.273861278752583*(alpha[1]*favg[17]+alpha[5]*favg[16])+(0.1414213562373095*alpha[13]+0.1581138830084189*alpha[3])*favg[15]+alpha[1]*(0.3952847075210473*favg[14]+0.1581138830084189*favg[13])+(0.3535533905932737*favg[12]+0.273861278752583*favg[4]+0.1581138830084189*favg[1])*alpha[13]+0.3952847075210473*alpha[3]*favg[12]+0.273861278752583*alpha[5]*favg[11]-0.8660254037844386*fjump[10]+(0.273861278752583*alpha[7]+0.3061862178478971*alpha[0])*favg[10]+alpha[5]*(0.1581138830084189*favg[9]+0.3952847075210473*favg[8])+0.1581138830084189*(alpha[5]*favg[7]+favg[5]*alpha[7])+0.3061862178478971*alpha[1]*favg[6]-0.5*fjump[5]+0.1767766952966368*alpha[0]*favg[5]+(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[5]+0.3061862178478971*alpha[3]*favg[4]+0.1767766952966368*(alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] += 0.3535533905932737*alpha[5]*favg[18]+(0.1956151991089878*alpha[13]+0.3061862178478971*alpha[3])*favg[17]+0.3952847075210473*alpha[13]*favg[14]+(0.1129384878631564*alpha[13]+0.1767766952966368*alpha[3])*favg[13]+(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])*alpha[13]+0.3535533905932737*alpha[1]*favg[12]-0.8660254037844386*fjump[11]+(0.1956151991089878*alpha[7]+0.3061862178478971*alpha[0])*favg[11]+0.273861278752583*alpha[5]*favg[10]+0.3952847075210473*alpha[7]*favg[8]-0.5*fjump[7]+(0.1129384878631564*alpha[7]+0.1767766952966368*alpha[0])*favg[7]+(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[7]+0.1581138830084189*alpha[5]*favg[5]+alpha[1]*(0.273861278752583*favg[4]+0.1581138830084189*favg[1]); 
  Ghat[9] += 0.3061862178478971*alpha[1]*favg[19]+0.3535533905932737*alpha[5]*favg[18]+0.273861278752583*alpha[13]*favg[17]-0.8660254037844386*fjump[16]+0.3061862178478971*alpha[0]*favg[16]+0.1767766952966368*alpha[1]*favg[15]+0.3535533905932737*alpha[3]*favg[14]+0.1581138830084189*alpha[13]*favg[13]+0.273861278752583*alpha[5]*favg[10]-0.5*fjump[9]+0.1767766952966368*alpha[0]*favg[9]+0.273861278752583*alpha[3]*favg[6]+0.1581138830084189*(alpha[5]*favg[5]+alpha[3]*favg[3]); 
  Ghat[13] += 0.2449489742783178*alpha[5]*favg[19]+0.3535533905932737*alpha[1]*favg[18]-0.8660254037844386*fjump[17]+(0.1956151991089878*alpha[7]+0.3061862178478971*alpha[0])*favg[17]+0.273861278752583*alpha[13]*favg[16]+0.1414213562373095*alpha[5]*favg[15]+0.3952847075210473*alpha[7]*favg[14]-0.5*fjump[13]+(0.1129384878631564*alpha[7]+0.1767766952966368*alpha[0])*favg[13]+(0.1956151991089878*favg[11]+0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.1129384878631564*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[13]+0.3535533905932737*alpha[5]*favg[12]+0.3061862178478971*alpha[3]*favg[11]+0.273861278752583*alpha[1]*favg[10]+0.1767766952966368*alpha[3]*favg[7]+(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])*alpha[7]+0.1581138830084189*alpha[1]*favg[5]+(0.273861278752583*favg[4]+0.1581138830084189*favg[1])*alpha[5]; 
  Ghat[15] += (-0.8660254037844386*fjump[19])+(0.273861278752583*alpha[7]+0.3061862178478971*alpha[0])*favg[19]+(0.3162277660168379*alpha[13]+0.3535533905932737*alpha[3])*favg[18]+0.2449489742783178*alpha[5]*favg[17]+0.3061862178478971*alpha[1]*favg[16]-0.5*fjump[15]+(0.1581138830084189*alpha[7]+0.1767766952966368*alpha[0])*favg[15]+alpha[5]*(0.3535533905932737*favg[14]+0.1414213562373095*favg[13])+(0.2449489742783178*favg[10]+0.1414213562373095*favg[5])*alpha[13]+0.273861278752583*alpha[3]*favg[10]+0.1767766952966368*alpha[1]*favg[9]+0.273861278752583*alpha[5]*favg[6]+0.1581138830084189*(alpha[3]*favg[5]+favg[3]*alpha[5]); 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[3] += 0.5*Ghat[3]*dv10r; 
  outr[4] += -0.8660254037844386*Ghat[1]*dv10r; 
  outr[5] += 0.5*Ghat[5]*dv10r; 
  outr[6] += -0.8660254037844386*Ghat[3]*dv10r; 
  outr[7] += 0.5*Ghat[7]*dv10r; 
  outr[8] += 1.118033988749895*Ghat[0]*dv10r; 
  outr[9] += 0.5*Ghat[9]*dv10r; 
  outr[10] += -0.8660254037844386*Ghat[5]*dv10r; 
  outr[11] += -0.8660254037844387*Ghat[7]*dv10r; 
  outr[12] += 1.118033988749895*Ghat[1]*dv10r; 
  outr[13] += 0.5*Ghat[13]*dv10r; 
  outr[14] += 1.118033988749895*Ghat[3]*dv10r; 
  outr[15] += 0.5*Ghat[15]*dv10r; 
  outr[16] += -0.8660254037844387*Ghat[9]*dv10r; 
  outr[17] += -0.8660254037844387*Ghat[13]*dv10r; 
  outr[18] += 1.118033988749895*Ghat[5]*dv10r; 
  outr[19] += -0.8660254037844387*Ghat[15]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[3] += -0.5*Ghat[3]*dv10l; 
  outl[4] += -0.8660254037844386*Ghat[1]*dv10l; 
  outl[5] += -0.5*Ghat[5]*dv10l; 
  outl[6] += -0.8660254037844386*Ghat[3]*dv10l; 
  outl[7] += -0.5*Ghat[7]*dv10l; 
  outl[8] += -1.118033988749895*Ghat[0]*dv10l; 
  outl[9] += -0.5*Ghat[9]*dv10l; 
  outl[10] += -0.8660254037844386*Ghat[5]*dv10l; 
  outl[11] += -0.8660254037844387*Ghat[7]*dv10l; 
  outl[12] += -1.118033988749895*Ghat[1]*dv10l; 
  outl[13] += -0.5*Ghat[13]*dv10l; 
  outl[14] += -1.118033988749895*Ghat[3]*dv10l; 
  outl[15] += -0.5*Ghat[15]*dv10l; 
  outl[16] += -0.8660254037844387*Ghat[9]*dv10l; 
  outl[17] += -0.8660254037844387*Ghat[13]*dv10l; 
  outl[18] += -1.118033988749895*Ghat[5]*dv10l; 
  outl[19] += -0.8660254037844387*Ghat[15]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double Ghat[32]; 

  for(unsigned int i=0; i<32; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[32]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = -1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = -1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = -1*fr[20]+fl[20]; 
  favg[21] = 1*fr[21]+fl[21]; 
  favg[22] = -1*fr[22]+fl[22]; 
  favg[23] = -1*fr[23]+fl[23]; 
  favg[24] = -1*fr[24]+fl[24]; 
  favg[25] = 1*fr[25]+fl[25]; 
  favg[26] = -1*fr[26]+fl[26]; 
  favg[27] = 1*fr[27]+fl[27]; 
  favg[28] = -1*fr[28]+fl[28]; 
  favg[29] = -1*fr[29]+fl[29]; 
  favg[30] = -1*fr[30]+fl[30]; 
  favg[31] = -1*fr[31]+fl[31]; 
  double fjump[32]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(-1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(1*fr[9]-fl[9]); 
  fjump[10] = amax*(-1*fr[10]-fl[10]); 
  fjump[11] = amax*(-1*fr[11]-fl[11]); 
  fjump[12] = amax*(1*fr[12]-fl[12]); 
  fjump[13] = amax*(1*fr[13]-fl[13]); 
  fjump[14] = amax*(1*fr[14]-fl[14]); 
  fjump[15] = amax*(1*fr[15]-fl[15]); 
  fjump[16] = amax*(-1*fr[16]-fl[16]); 
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(-1*fr[18]-fl[18]); 
  fjump[19] = amax*(1*fr[19]-fl[19]); 
  fjump[20] = amax*(-1*fr[20]-fl[20]); 
  fjump[21] = amax*(1*fr[21]-fl[21]); 
  fjump[22] = amax*(-1*fr[22]-fl[22]); 
  fjump[23] = amax*(-1*fr[23]-fl[23]); 
  fjump[24] = amax*(-1*fr[24]-fl[24]); 
  fjump[25] = amax*(1*fr[25]-fl[25]); 
  fjump[26] = amax*(-1*fr[26]-fl[26]); 
  fjump[27] = amax*(1*fr[27]-fl[27]); 
  fjump[28] = amax*(-1*fr[28]-fl[28]); 
  fjump[29] = amax*(-1*fr[29]-fl[29]); 
  fjump[30] = amax*(-1*fr[30]-fl[30]); 
  fjump[31] = amax*(-1*fr[31]-fl[31]); 
  double alpha[32]; 

  alpha[0] = 2.0*B2[0]*wv2+2.0*E0[0]; 
  alpha[1] = 2.0*B2[1]*wv2+2.0*E0[1]; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[5] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 2.0*B2[2]*wv2+2.0*E0[2]; 
  alpha[13] = 0.5773502691896257*B2[2]*dv2; 
  alpha[17] = 2.0*B2[3]*wv2+2.0*E0[3]; 
  alpha[25] = 0.5773502691896256*B2[3]*dv2; 
  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 
  Ghat[0] += 0.4677071733467425*alpha[5]*favg[30]+0.3061862178478971*alpha[25]*favg[29]+0.4677071733467425*alpha[3]*favg[26]+0.1767766952966368*alpha[25]*favg[25]+0.4677071733467425*alpha[1]*favg[24]+0.3061862178478971*alpha[17]*favg[23]+0.3952847075210473*alpha[5]*favg[21]+0.3061862178478971*alpha[13]*favg[20]-1.322875655532295*fjump[18]+0.4677071733467425*alpha[0]*favg[18]+0.1767766952966368*alpha[17]*favg[17]+0.3952847075210473*alpha[3]*favg[14]+0.1767766952966368*alpha[13]*favg[13]+0.3952847075210473*alpha[1]*favg[12]+0.3061862178478971*(alpha[7]*favg[11]+alpha[5]*favg[10])-1.118033988749895*fjump[8]+0.3952847075210473*alpha[0]*favg[8]+0.1767766952966368*alpha[7]*favg[7]+0.3061862178478971*alpha[3]*favg[6]+0.1767766952966368*alpha[5]*favg[5]+0.3061862178478971*alpha[1]*favg[4]+0.1767766952966368*alpha[3]*favg[3]-0.8660254037844386*fjump[2]+0.3061862178478971*alpha[0]*favg[2]+0.1767766952966368*alpha[1]*favg[1]-0.5*fjump[0]+0.1767766952966368*alpha[0]*favg[0]; 
  Ghat[1] += (0.4183300132670377*alpha[13]+0.4677071733467425*alpha[3])*favg[30]+0.2689264371002384*alpha[13]*favg[29]+0.4677071733467425*alpha[5]*favg[26]+0.1552647508520296*alpha[13]*favg[25]+(0.2689264371002384*favg[20]+0.1552647508520296*favg[13])*alpha[25]-1.322875655532295*fjump[24]+(0.4183300132670377*alpha[7]+0.4677071733467425*alpha[0])*favg[24]+0.2689264371002384*alpha[7]*favg[23]+(0.3535533905932737*alpha[13]+0.3952847075210473*alpha[3])*favg[21]+0.273861278752583*alpha[5]*favg[20]+0.4677071733467425*alpha[1]*favg[18]+0.1552647508520296*alpha[7]*favg[17]+(0.2689264371002384*favg[11]+0.1552647508520296*favg[7])*alpha[17]+alpha[5]*(0.3952847075210473*favg[14]+0.1581138830084189*favg[13])+(0.273861278752583*favg[10]+0.1581138830084189*favg[5])*alpha[13]-1.118033988749895*fjump[12]+(0.3535533905932737*alpha[7]+0.3952847075210473*alpha[0])*favg[12]+0.273861278752583*alpha[1]*favg[11]+0.3061862178478971*alpha[3]*favg[10]+alpha[1]*(0.3952847075210473*favg[8]+0.1581138830084189*favg[7])+(0.273861278752583*favg[4]+0.1581138830084189*favg[1])*alpha[7]+0.3061862178478971*alpha[5]*favg[6]+0.1767766952966368*(alpha[3]*favg[5]+favg[3]*alpha[5])-0.8660254037844386*fjump[4]+0.3061862178478971*(alpha[0]*favg[4]+alpha[1]*favg[2])-0.5*fjump[1]+0.1767766952966368*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[3] += 0.4677071733467425*alpha[1]*favg[30]+0.3061862178478971*alpha[17]*favg[29]-1.322875655532295*fjump[26]+0.4677071733467425*alpha[0]*favg[26]+0.1767766952966368*alpha[17]*favg[25]+(0.3061862178478971*favg[23]+0.1767766952966368*favg[17])*alpha[25]+alpha[5]*(0.4677071733467425*favg[24]+0.273861278752583*favg[22])+0.3952847075210473*alpha[1]*favg[21]+0.3061862178478971*alpha[7]*favg[20]+alpha[3]*(0.4677071733467425*favg[18]+0.273861278752583*favg[16])+0.1581138830084189*alpha[5]*favg[15]-1.118033988749895*fjump[14]+0.3952847075210473*alpha[0]*favg[14]+0.1767766952966368*alpha[7]*favg[13]+(0.3061862178478971*favg[11]+0.1767766952966368*favg[7])*alpha[13]+0.3952847075210473*alpha[5]*favg[12]+0.3061862178478971*alpha[1]*favg[10]+alpha[3]*(0.1581138830084189*favg[9]+0.3952847075210473*favg[8])-0.8660254037844386*fjump[6]+0.3061862178478971*alpha[0]*favg[6]+0.1767766952966368*alpha[1]*favg[5]+(0.3061862178478971*favg[4]+0.1767766952966368*favg[1])*alpha[5]-0.5*fjump[3]+0.1767766952966368*alpha[0]*favg[3]+(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[3]; 
  Ghat[5] += (-1.322875655532295*fjump[30])+(0.4183300132670377*alpha[7]+0.4677071733467425*alpha[0])*favg[30]+0.2689264371002384*alpha[7]*favg[29]+0.4677071733467425*alpha[1]*favg[26]+0.1552647508520296*alpha[7]*favg[25]+(0.2689264371002384*favg[11]+0.1552647508520296*favg[7])*alpha[25]+(0.4183300132670377*alpha[13]+0.4677071733467425*alpha[3])*favg[24]+0.2689264371002384*alpha[13]*favg[23]+(0.2449489742783178*alpha[13]+0.273861278752583*alpha[3])*favg[22]-1.118033988749895*fjump[21]+(0.3535533905932737*alpha[7]+0.3952847075210473*alpha[0])*favg[21]+(0.2689264371002384*alpha[17]+0.273861278752583*alpha[1])*favg[20]+0.4677071733467425*alpha[5]*favg[18]+0.1552647508520296*(alpha[13]*favg[17]+favg[13]*alpha[17])+0.273861278752583*alpha[5]*favg[16]+(0.1414213562373095*alpha[13]+0.1581138830084189*alpha[3])*favg[15]+alpha[1]*(0.3952847075210473*favg[14]+0.1581138830084189*favg[13])+(0.3535533905932737*favg[12]+0.273861278752583*favg[4]+0.1581138830084189*favg[1])*alpha[13]+0.3952847075210473*alpha[3]*favg[12]+0.273861278752583*alpha[5]*favg[11]-0.8660254037844386*fjump[10]+(0.273861278752583*alpha[7]+0.3061862178478971*alpha[0])*favg[10]+alpha[5]*(0.1581138830084189*favg[9]+0.3952847075210473*favg[8])+0.1581138830084189*(alpha[5]*favg[7]+favg[5]*alpha[7])+0.3061862178478971*alpha[1]*favg[6]-0.5*fjump[5]+0.1767766952966368*alpha[0]*favg[5]+(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[5]+0.3061862178478971*alpha[3]*favg[4]+0.1767766952966368*(alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] += (0.4107919181288743*alpha[25]+0.4183300132670377*alpha[5])*favg[30]+(0.1825741858350554*alpha[25]+0.2689264371002384*alpha[5])*favg[29]+0.4677071733467425*alpha[13]*favg[26]+(0.105409255338946*alpha[25]+0.1552647508520296*alpha[5])*favg[25]+(0.3471825374147066*favg[21]+0.2689264371002384*favg[10]+0.1552647508520296*favg[5])*alpha[25]+(0.4107919181288743*alpha[17]+0.4183300132670377*alpha[1])*favg[24]+(0.1825741858350554*alpha[17]+0.2689264371002384*alpha[1])*favg[23]+0.3535533905932737*alpha[5]*favg[21]+(0.1956151991089878*alpha[13]+0.3061862178478971*alpha[3])*favg[20]+0.4677071733467425*alpha[7]*favg[18]+(0.105409255338946*alpha[17]+0.1552647508520296*alpha[1])*favg[17]+(0.3471825374147066*favg[12]+0.2689264371002384*favg[4]+0.1552647508520296*favg[1])*alpha[17]+0.3952847075210473*alpha[13]*favg[14]+(0.1129384878631564*alpha[13]+0.1767766952966368*alpha[3])*favg[13]+(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])*alpha[13]+0.3535533905932737*alpha[1]*favg[12]-0.8660254037844386*fjump[11]+(0.1956151991089878*alpha[7]+0.3061862178478971*alpha[0])*favg[11]+0.273861278752583*alpha[5]*favg[10]+0.3952847075210473*alpha[7]*favg[8]-0.5*fjump[7]+(0.1129384878631564*alpha[7]+0.1767766952966368*alpha[0])*favg[7]+(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[7]+0.1581138830084189*alpha[5]*favg[5]+alpha[1]*(0.273861278752583*favg[4]+0.1581138830084189*favg[1]); 
  Ghat[9] += alpha[5]*(0.2689264371002384*favg[31]+0.4183300132670377*favg[30])+0.273861278752583*alpha[25]*favg[29]+0.2689264371002384*alpha[3]*favg[28]+0.1552647508520296*alpha[5]*favg[27]+0.4183300132670377*alpha[3]*favg[26]+0.1581138830084189*alpha[25]*favg[25]+0.3061862178478971*alpha[1]*favg[22]+0.3535533905932737*alpha[5]*favg[21]+0.273861278752583*alpha[13]*favg[20]+0.1552647508520296*alpha[3]*favg[19]-0.8660254037844386*fjump[16]+0.3061862178478971*alpha[0]*favg[16]+0.1767766952966368*alpha[1]*favg[15]+0.3535533905932737*alpha[3]*favg[14]+0.1581138830084189*alpha[13]*favg[13]+0.273861278752583*alpha[5]*favg[10]-0.5*fjump[9]+0.1767766952966368*alpha[0]*favg[9]+0.273861278752583*alpha[3]*favg[6]+0.1581138830084189*(alpha[5]*favg[5]+alpha[3]*favg[3]); 
  Ghat[13] += (0.4107919181288743*alpha[17]+0.4183300132670377*alpha[1])*favg[30]+(0.1825741858350554*alpha[17]+0.2689264371002384*alpha[1])*favg[29]+0.4677071733467425*alpha[7]*favg[26]+(0.105409255338946*alpha[17]+0.1552647508520296*alpha[1])*favg[25]+(0.4107919181288743*favg[24]+0.1825741858350554*favg[23]+0.2405351177211819*favg[22]+0.105409255338946*favg[17]+0.1388730149658826*favg[15]+0.3471825374147066*favg[12]+0.2689264371002384*favg[4]+0.1552647508520296*favg[1])*alpha[25]+alpha[5]*(0.4183300132670377*favg[24]+0.2689264371002384*favg[23]+0.2449489742783178*favg[22])+(0.3471825374147066*alpha[17]+0.3535533905932737*alpha[1])*favg[21]-0.8660254037844386*fjump[20]+(0.1956151991089878*alpha[7]+0.3061862178478971*alpha[0])*favg[20]+0.4677071733467425*alpha[13]*favg[18]+0.1552647508520296*alpha[5]*favg[17]+(0.2689264371002384*favg[10]+0.1552647508520296*favg[5])*alpha[17]+0.273861278752583*alpha[13]*favg[16]+0.1414213562373095*alpha[5]*favg[15]+0.3952847075210473*alpha[7]*favg[14]-0.5*fjump[13]+(0.1129384878631564*alpha[7]+0.1767766952966368*alpha[0])*favg[13]+(0.1956151991089878*favg[11]+0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.1129384878631564*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[13]+0.3535533905932737*alpha[5]*favg[12]+0.3061862178478971*alpha[3]*favg[11]+0.273861278752583*alpha[1]*favg[10]+0.1767766952966368*alpha[3]*favg[7]+(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])*alpha[7]+0.1581138830084189*alpha[1]*favg[5]+(0.273861278752583*favg[4]+0.1581138830084189*favg[1])*alpha[5]; 
  Ghat[15] += (0.2405351177211819*alpha[13]+0.2689264371002384*alpha[3])*favg[31]+(0.3741657386773942*alpha[13]+0.4183300132670377*alpha[3])*favg[30]+0.2405351177211819*alpha[13]*favg[29]+0.2689264371002384*alpha[5]*favg[28]+(0.1388730149658826*alpha[13]+0.1552647508520296*alpha[3])*favg[27]+0.4183300132670377*alpha[5]*favg[26]+0.1388730149658826*alpha[13]*favg[25]+(0.2405351177211819*favg[20]+0.1388730149658826*favg[13])*alpha[25]-0.8660254037844386*fjump[22]+(0.273861278752583*alpha[7]+0.3061862178478971*alpha[0])*favg[22]+(0.3162277660168379*alpha[13]+0.3535533905932737*alpha[3])*favg[21]+alpha[5]*(0.2449489742783178*favg[20]+0.1552647508520296*favg[19])+0.3061862178478971*alpha[1]*favg[16]-0.5*fjump[15]+(0.1581138830084189*alpha[7]+0.1767766952966368*alpha[0])*favg[15]+alpha[5]*(0.3535533905932737*favg[14]+0.1414213562373095*favg[13])+(0.2449489742783178*favg[10]+0.1414213562373095*favg[5])*alpha[13]+0.273861278752583*alpha[3]*favg[10]+0.1767766952966368*alpha[1]*favg[9]+0.273861278752583*alpha[5]*favg[6]+0.1581138830084189*(alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[17] += 0.4107919181288743*alpha[13]*favg[30]+(0.1825741858350554*alpha[13]+0.3061862178478971*alpha[3])*favg[29]+0.4677071733467425*alpha[25]*favg[26]+(0.105409255338946*alpha[13]+0.1767766952966368*alpha[3])*favg[25]+(0.1825741858350554*favg[20]+0.3952847075210473*favg[14]+0.105409255338946*favg[13]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])*alpha[25]+0.4107919181288743*alpha[7]*favg[24]-0.8660254037844386*fjump[23]+(0.1825741858350554*alpha[7]+0.3061862178478971*alpha[0])*favg[23]+0.3471825374147066*alpha[13]*favg[21]+0.2689264371002384*alpha[5]*favg[20]+0.4677071733467425*alpha[17]*favg[18]-0.5*fjump[17]+(0.105409255338946*alpha[7]+0.1767766952966368*alpha[0])*favg[17]+(0.1825741858350554*favg[11]+0.3952847075210473*favg[8]+0.105409255338946*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[17]+0.1552647508520296*alpha[5]*favg[13]+(0.2689264371002384*favg[10]+0.1552647508520296*favg[5])*alpha[13]+0.3471825374147066*alpha[7]*favg[12]+alpha[1]*(0.2689264371002384*favg[11]+0.1552647508520296*favg[7])+(0.2689264371002384*favg[4]+0.1552647508520296*favg[1])*alpha[7]; 
  Ghat[19] += 0.3061862178478971*alpha[1]*favg[31]-0.8660254037844386*fjump[28]+0.3061862178478971*alpha[0]*favg[28]+0.1767766952966368*alpha[1]*favg[27]+0.2689264371002384*alpha[5]*favg[22]-0.5*fjump[19]+0.1767766952966368*alpha[0]*favg[19]+0.2689264371002384*alpha[3]*favg[16]+0.1552647508520296*(alpha[5]*favg[15]+alpha[3]*favg[9]); 
  Ghat[25] += 0.4107919181288743*alpha[7]*favg[30]-0.8660254037844386*fjump[29]+(0.1825741858350554*alpha[7]+0.3061862178478971*alpha[0])*favg[29]+0.4677071733467425*alpha[17]*favg[26]-0.5*fjump[25]+(0.105409255338946*alpha[7]+0.1767766952966368*alpha[0])*favg[25]+(0.4677071733467425*favg[18]+0.273861278752583*favg[16]+0.1825741858350554*favg[11]+0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.105409255338946*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[25]+0.4107919181288743*alpha[13]*favg[24]+(0.1825741858350554*alpha[13]+0.3061862178478971*alpha[3])*favg[23]+0.2405351177211819*alpha[13]*favg[22]+0.3471825374147066*alpha[7]*favg[21]+(0.1825741858350554*alpha[17]+0.2689264371002384*alpha[1])*favg[20]+(0.105409255338946*alpha[13]+0.1767766952966368*alpha[3])*favg[17]+(0.3952847075210473*favg[14]+0.105409255338946*favg[13]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])*alpha[17]+0.1388730149658826*alpha[13]*favg[15]+0.1552647508520296*alpha[1]*favg[13]+(0.3471825374147066*favg[12]+0.2689264371002384*favg[4]+0.1552647508520296*favg[1])*alpha[13]+0.2689264371002384*(alpha[5]*favg[11]+alpha[7]*favg[10])+0.1552647508520296*(alpha[5]*favg[7]+favg[5]*alpha[7]); 
  Ghat[27] += (-0.8660254037844386*fjump[31])+0.273861278752583*alpha[7]*favg[31]+0.3061862178478971*(alpha[0]*favg[31]+alpha[1]*favg[28])-0.5*fjump[27]+(0.1581138830084189*alpha[7]+0.1767766952966368*alpha[0])*favg[27]+(0.2405351177211819*alpha[13]+0.2689264371002384*alpha[3])*favg[22]+0.1767766952966368*alpha[1]*favg[19]+0.2689264371002384*alpha[5]*favg[16]+0.1388730149658826*alpha[13]*favg[15]+0.1552647508520296*(alpha[3]*favg[15]+alpha[5]*favg[9]); 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[3] += 0.5*Ghat[3]*dv10r; 
  outr[4] += -0.8660254037844386*Ghat[1]*dv10r; 
  outr[5] += 0.5*Ghat[5]*dv10r; 
  outr[6] += -0.8660254037844386*Ghat[3]*dv10r; 
  outr[7] += 0.5*Ghat[7]*dv10r; 
  outr[8] += 1.118033988749895*Ghat[0]*dv10r; 
  outr[9] += 0.5*Ghat[9]*dv10r; 
  outr[10] += -0.8660254037844386*Ghat[5]*dv10r; 
  outr[11] += -0.8660254037844387*Ghat[7]*dv10r; 
  outr[12] += 1.118033988749895*Ghat[1]*dv10r; 
  outr[13] += 0.5*Ghat[13]*dv10r; 
  outr[14] += 1.118033988749895*Ghat[3]*dv10r; 
  outr[15] += 0.5*Ghat[15]*dv10r; 
  outr[16] += -0.8660254037844387*Ghat[9]*dv10r; 
  outr[17] += 0.5*Ghat[17]*dv10r; 
  outr[18] += -1.322875655532295*Ghat[0]*dv10r; 
  outr[19] += 0.5*Ghat[19]*dv10r; 
  outr[20] += -0.8660254037844387*Ghat[13]*dv10r; 
  outr[21] += 1.118033988749895*Ghat[5]*dv10r; 
  outr[22] += -0.8660254037844387*Ghat[15]*dv10r; 
  outr[23] += -0.8660254037844386*Ghat[17]*dv10r; 
  outr[24] += -1.322875655532295*Ghat[1]*dv10r; 
  outr[25] += 0.5*Ghat[25]*dv10r; 
  outr[26] += -1.322875655532295*Ghat[3]*dv10r; 
  outr[27] += 0.5*Ghat[27]*dv10r; 
  outr[28] += -0.8660254037844386*Ghat[19]*dv10r; 
  outr[29] += -0.8660254037844386*Ghat[25]*dv10r; 
  outr[30] += -1.322875655532295*Ghat[5]*dv10r; 
  outr[31] += -0.8660254037844386*Ghat[27]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[3] += -0.5*Ghat[3]*dv10l; 
  outl[4] += -0.8660254037844386*Ghat[1]*dv10l; 
  outl[5] += -0.5*Ghat[5]*dv10l; 
  outl[6] += -0.8660254037844386*Ghat[3]*dv10l; 
  outl[7] += -0.5*Ghat[7]*dv10l; 
  outl[8] += -1.118033988749895*Ghat[0]*dv10l; 
  outl[9] += -0.5*Ghat[9]*dv10l; 
  outl[10] += -0.8660254037844386*Ghat[5]*dv10l; 
  outl[11] += -0.8660254037844387*Ghat[7]*dv10l; 
  outl[12] += -1.118033988749895*Ghat[1]*dv10l; 
  outl[13] += -0.5*Ghat[13]*dv10l; 
  outl[14] += -1.118033988749895*Ghat[3]*dv10l; 
  outl[15] += -0.5*Ghat[15]*dv10l; 
  outl[16] += -0.8660254037844387*Ghat[9]*dv10l; 
  outl[17] += -0.5*Ghat[17]*dv10l; 
  outl[18] += -1.322875655532295*Ghat[0]*dv10l; 
  outl[19] += -0.5*Ghat[19]*dv10l; 
  outl[20] += -0.8660254037844387*Ghat[13]*dv10l; 
  outl[21] += -1.118033988749895*Ghat[5]*dv10l; 
  outl[22] += -0.8660254037844387*Ghat[15]*dv10l; 
  outl[23] += -0.8660254037844386*Ghat[17]*dv10l; 
  outl[24] += -1.322875655532295*Ghat[1]*dv10l; 
  outl[25] += -0.5*Ghat[25]*dv10l; 
  outl[26] += -1.322875655532295*Ghat[3]*dv10l; 
  outl[27] += -0.5*Ghat[27]*dv10l; 
  outl[28] += -0.8660254037844386*Ghat[19]*dv10l; 
  outl[29] += -0.8660254037844386*Ghat[25]*dv10l; 
  outl[30] += -1.322875655532295*Ghat[5]*dv10l; 
  outl[31] += -0.8660254037844386*Ghat[27]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[2]; 
  double dv11r = 2/dxvr[2]; 
  const double *E1 = &EM[2]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double Ghat[8]; 

  for(unsigned int i=0; i<8; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[8]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  double fjump[8]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(-1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(-1*fr[7]-fl[7]); 
  double alpha[8]; 

  alpha[0] = 2.0*E1[0]-2.0*B2[0]*wv1; 
  alpha[1] = 2.0*E1[1]-2.0*B2[1]*wv1; 
  alpha[2] = -0.5773502691896258*B2[0]*dv1; 
  alpha[4] = -0.5773502691896258*B2[1]*dv1; 
  const double amid = 0.3535533905932737*alpha[0]; 
  Ghat[0] += 0.3061862178478971*(alpha[4]*favg[7]+alpha[2]*favg[6]+alpha[1]*favg[5])+0.1767766952966368*alpha[4]*favg[4]-0.8660254037844386*fjump[3]+0.3061862178478971*alpha[0]*favg[3]+0.1767766952966368*(alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.1767766952966368*alpha[0]*favg[0]; 
  Ghat[1] += 0.3061862178478971*(alpha[2]*favg[7]+alpha[4]*favg[6])-0.8660254037844386*fjump[5]+0.3061862178478971*alpha[0]*favg[5]+0.1767766952966368*(alpha[2]*favg[4]+favg[2]*alpha[4])+0.3061862178478971*alpha[1]*favg[3]-0.5*fjump[1]+0.1767766952966368*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += 0.3061862178478971*alpha[1]*favg[7]-0.8660254037844386*fjump[6]+0.3061862178478971*(alpha[0]*favg[6]+alpha[4]*favg[5])+0.1767766952966368*(alpha[1]*favg[4]+favg[1]*alpha[4])+0.3061862178478971*alpha[2]*favg[3]-0.5*fjump[2]+0.1767766952966368*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[4] += (-0.8660254037844386*fjump[7])+0.3061862178478971*(alpha[0]*favg[7]+alpha[1]*favg[6]+alpha[2]*favg[5])-0.5*fjump[4]+0.1767766952966368*alpha[0]*favg[4]+0.3061862178478971*favg[3]*alpha[4]+0.1767766952966368*(favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]); 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += -0.8660254037844386*Ghat[0]*dv11r; 
  outr[4] += 0.5*Ghat[4]*dv11r; 
  outr[5] += -0.8660254037844386*Ghat[1]*dv11r; 
  outr[6] += -0.8660254037844386*Ghat[2]*dv11r; 
  outr[7] += -0.8660254037844386*Ghat[4]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.8660254037844386*Ghat[0]*dv11l; 
  outl[4] += -0.5*Ghat[4]*dv11l; 
  outl[5] += -0.8660254037844386*Ghat[1]*dv11l; 
  outl[6] += -0.8660254037844386*Ghat[2]*dv11l; 
  outl[7] += -0.8660254037844386*Ghat[4]*dv11l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[2]; 
  double dv11r = 2/dxvr[2]; 
  const double *E1 = &EM[3]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[20]; 

  for(unsigned int i=0; i<20; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[20]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = -1*fr[17]+fl[17]; 
  favg[18] = -1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  double fjump[20]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(-1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(1*fr[9]-fl[9]); 
  fjump[10] = amax*(-1*fr[10]-fl[10]); 
  fjump[11] = amax*(1*fr[11]-fl[11]); 
  fjump[12] = amax*(1*fr[12]-fl[12]); 
  fjump[13] = amax*(-1*fr[13]-fl[13]); 
  fjump[14] = amax*(-1*fr[14]-fl[14]); 
  fjump[15] = amax*(1*fr[15]-fl[15]); 
  fjump[16] = amax*(1*fr[16]-fl[16]); 
  fjump[17] = amax*(-1*fr[17]-fl[17]); 
  fjump[18] = amax*(-1*fr[18]-fl[18]); 
  fjump[19] = amax*(1*fr[19]-fl[19]); 
  double alpha[20]; 

  alpha[0] = 2.0*E1[0]-2.0*B2[0]*wv1; 
  alpha[1] = 2.0*E1[1]-2.0*B2[1]*wv1; 
  alpha[2] = -0.5773502691896258*B2[0]*dv1; 
  alpha[4] = -0.5773502691896258*B2[1]*dv1; 
  alpha[7] = 2.0*E1[2]-2.0*B2[2]*wv1; 
  alpha[11] = -0.5773502691896257*B2[2]*dv1; 
  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 
  Ghat[0] += 0.3952847075210473*alpha[4]*favg[19]+0.3061862178478971*alpha[11]*favg[17]+0.3952847075210473*(alpha[2]*favg[16]+alpha[1]*favg[15])+0.3061862178478971*alpha[7]*favg[13]+0.1767766952966368*alpha[11]*favg[11]+0.3061862178478971*alpha[4]*favg[10]-1.118033988749895*fjump[9]+0.3952847075210473*alpha[0]*favg[9]+0.1767766952966368*alpha[7]*favg[7]+0.3061862178478971*(alpha[2]*favg[6]+alpha[1]*favg[5])+0.1767766952966368*alpha[4]*favg[4]-0.8660254037844386*fjump[3]+0.3061862178478971*alpha[0]*favg[3]+0.1767766952966368*(alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.1767766952966368*alpha[0]*favg[0]; 
  Ghat[1] += (0.3535533905932737*alpha[11]+0.3952847075210473*alpha[2])*favg[19]+alpha[4]*(0.273861278752583*favg[17]+0.3952847075210473*favg[16])-1.118033988749895*fjump[15]+(0.3535533905932737*alpha[7]+0.3952847075210473*alpha[0])*favg[15]+0.273861278752583*alpha[1]*favg[13]+0.1581138830084189*alpha[4]*favg[11]+(0.273861278752583*favg[10]+0.1581138830084189*favg[4])*alpha[11]+0.3061862178478971*alpha[2]*favg[10]+alpha[1]*(0.3952847075210473*favg[9]+0.1581138830084189*favg[7])+(0.273861278752583*favg[5]+0.1581138830084189*favg[1])*alpha[7]+0.3061862178478971*alpha[4]*favg[6]-0.8660254037844386*fjump[5]+0.3061862178478971*alpha[0]*favg[5]+0.1767766952966368*(alpha[2]*favg[4]+favg[2]*alpha[4])+0.3061862178478971*alpha[1]*favg[3]-0.5*fjump[1]+0.1767766952966368*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += 0.3952847075210473*alpha[1]*favg[19]+0.273861278752583*alpha[4]*favg[18]+0.3061862178478971*alpha[7]*favg[17]-1.118033988749895*fjump[16]+0.3952847075210473*(alpha[0]*favg[16]+alpha[4]*favg[15])+0.273861278752583*alpha[2]*favg[14]+0.3061862178478971*alpha[11]*favg[13]+0.1581138830084189*alpha[4]*favg[12]+0.1767766952966368*(alpha[7]*favg[11]+favg[7]*alpha[11])+0.3061862178478971*alpha[1]*favg[10]+alpha[2]*(0.3952847075210473*favg[9]+0.1581138830084189*favg[8])-0.8660254037844386*fjump[6]+0.3061862178478971*(alpha[0]*favg[6]+alpha[4]*favg[5])+0.1767766952966368*(alpha[1]*favg[4]+favg[1]*alpha[4])+0.3061862178478971*alpha[2]*favg[3]-0.5*fjump[2]+0.1767766952966368*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[4] += (-1.118033988749895*fjump[19])+(0.3535533905932737*alpha[7]+0.3952847075210473*alpha[0])*favg[19]+(0.2449489742783178*alpha[11]+0.273861278752583*alpha[2])*favg[18]+alpha[1]*(0.273861278752583*favg[17]+0.3952847075210473*favg[16])+(0.3535533905932737*alpha[11]+0.3952847075210473*alpha[2])*favg[15]+0.273861278752583*alpha[4]*(favg[14]+favg[13])+0.1414213562373095*alpha[11]*favg[12]+0.1581138830084189*(alpha[2]*favg[12]+alpha[1]*favg[11])+(0.273861278752583*favg[5]+0.1581138830084189*favg[1])*alpha[11]-0.8660254037844386*fjump[10]+(0.273861278752583*alpha[7]+0.3061862178478971*alpha[0])*favg[10]+0.3952847075210473*alpha[4]*favg[9]+0.1581138830084189*(alpha[4]*(favg[8]+favg[7])+favg[4]*alpha[7])+0.3061862178478971*(alpha[1]*favg[6]+alpha[2]*favg[5])-0.5*fjump[4]+0.1767766952966368*alpha[0]*favg[4]+0.3061862178478971*favg[3]*alpha[4]+0.1767766952966368*(favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[7] += 0.3535533905932737*alpha[4]*favg[19]+(0.1956151991089878*alpha[11]+0.3061862178478971*alpha[2])*favg[17]+0.3952847075210473*alpha[11]*favg[16]+0.3535533905932737*alpha[1]*favg[15]-0.8660254037844386*fjump[13]+(0.1956151991089878*alpha[7]+0.3061862178478971*alpha[0])*favg[13]+(0.1129384878631564*alpha[11]+0.1767766952966368*alpha[2])*favg[11]+(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])*alpha[11]+0.273861278752583*alpha[4]*favg[10]+0.3952847075210473*alpha[7]*favg[9]-0.5*fjump[7]+(0.1129384878631564*alpha[7]+0.1767766952966368*alpha[0])*favg[7]+(0.3061862178478971*favg[3]+0.1767766952966368*favg[0])*alpha[7]+0.273861278752583*alpha[1]*favg[5]+0.1581138830084189*(alpha[4]*favg[4]+alpha[1]*favg[1]); 
  Ghat[8] += 0.3535533905932737*alpha[4]*favg[19]+0.3061862178478971*alpha[1]*favg[18]+0.273861278752583*alpha[11]*favg[17]+0.3535533905932737*alpha[2]*favg[16]-0.8660254037844386*fjump[14]+0.3061862178478971*alpha[0]*favg[14]+0.1767766952966368*alpha[1]*favg[12]+0.1581138830084189*alpha[11]*favg[11]+0.273861278752583*alpha[4]*favg[10]-0.5*fjump[8]+0.1767766952966368*alpha[0]*favg[8]+0.273861278752583*alpha[2]*favg[6]+0.1581138830084189*(alpha[4]*favg[4]+alpha[2]*favg[2]); 
  Ghat[11] += 0.3535533905932737*alpha[1]*favg[19]+0.2449489742783178*alpha[4]*favg[18]-0.8660254037844386*fjump[17]+(0.1956151991089878*alpha[7]+0.3061862178478971*alpha[0])*favg[17]+0.3952847075210473*alpha[7]*favg[16]+0.3535533905932737*alpha[4]*favg[15]+0.273861278752583*alpha[11]*favg[14]+(0.1956151991089878*alpha[11]+0.3061862178478971*alpha[2])*favg[13]+0.1414213562373095*alpha[4]*favg[12]-0.5*fjump[11]+(0.1129384878631564*alpha[7]+0.1767766952966368*alpha[0])*favg[11]+(0.3952847075210473*favg[9]+0.1581138830084189*favg[8]+0.1129384878631564*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])*alpha[11]+0.273861278752583*alpha[1]*favg[10]+0.1767766952966368*alpha[2]*favg[7]+(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])*alpha[7]+0.273861278752583*alpha[4]*favg[5]+0.1581138830084189*(alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[12] += (0.3162277660168379*alpha[11]+0.3535533905932737*alpha[2])*favg[19]-0.8660254037844386*fjump[18]+(0.273861278752583*alpha[7]+0.3061862178478971*alpha[0])*favg[18]+alpha[4]*(0.2449489742783178*favg[17]+0.3535533905932737*favg[16])+0.3061862178478971*alpha[1]*favg[14]-0.5*fjump[12]+(0.1581138830084189*alpha[7]+0.1767766952966368*alpha[0])*favg[12]+0.1414213562373095*alpha[4]*favg[11]+(0.2449489742783178*favg[10]+0.1414213562373095*favg[4])*alpha[11]+0.273861278752583*alpha[2]*favg[10]+0.1767766952966368*alpha[1]*favg[8]+0.273861278752583*alpha[4]*favg[6]+0.1581138830084189*(alpha[2]*favg[4]+favg[2]*alpha[4]); 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += -0.8660254037844386*Ghat[0]*dv11r; 
  outr[4] += 0.5*Ghat[4]*dv11r; 
  outr[5] += -0.8660254037844386*Ghat[1]*dv11r; 
  outr[6] += -0.8660254037844386*Ghat[2]*dv11r; 
  outr[7] += 0.5*Ghat[7]*dv11r; 
  outr[8] += 0.5*Ghat[8]*dv11r; 
  outr[9] += 1.118033988749895*Ghat[0]*dv11r; 
  outr[10] += -0.8660254037844386*Ghat[4]*dv11r; 
  outr[11] += 0.5*Ghat[11]*dv11r; 
  outr[12] += 0.5*Ghat[12]*dv11r; 
  outr[13] += -0.8660254037844387*Ghat[7]*dv11r; 
  outr[14] += -0.8660254037844387*Ghat[8]*dv11r; 
  outr[15] += 1.118033988749895*Ghat[1]*dv11r; 
  outr[16] += 1.118033988749895*Ghat[2]*dv11r; 
  outr[17] += -0.8660254037844387*Ghat[11]*dv11r; 
  outr[18] += -0.8660254037844387*Ghat[12]*dv11r; 
  outr[19] += 1.118033988749895*Ghat[4]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.8660254037844386*Ghat[0]*dv11l; 
  outl[4] += -0.5*Ghat[4]*dv11l; 
  outl[5] += -0.8660254037844386*Ghat[1]*dv11l; 
  outl[6] += -0.8660254037844386*Ghat[2]*dv11l; 
  outl[7] += -0.5*Ghat[7]*dv11l; 
  outl[8] += -0.5*Ghat[8]*dv11l; 
  outl[9] += -1.118033988749895*Ghat[0]*dv11l; 
  outl[10] += -0.8660254037844386*Ghat[4]*dv11l; 
  outl[11] += -0.5*Ghat[11]*dv11l; 
  outl[12] += -0.5*Ghat[12]*dv11l; 
  outl[13] += -0.8660254037844387*Ghat[7]*dv11l; 
  outl[14] += -0.8660254037844387*Ghat[8]*dv11l; 
  outl[15] += -1.118033988749895*Ghat[1]*dv11l; 
  outl[16] += -1.118033988749895*Ghat[2]*dv11l; 
  outl[17] += -0.8660254037844387*Ghat[11]*dv11l; 
  outl[18] += -0.8660254037844387*Ghat[12]*dv11l; 
  outl[19] += -1.118033988749895*Ghat[4]*dv11l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vSer_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[2]; 
  double dv11r = 2/dxvr[2]; 
  const double *E1 = &EM[4]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double Ghat[32]; 

  for(unsigned int i=0; i<32; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[32]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = -1*fr[19]+fl[19]; 
  favg[20] = -1*fr[20]+fl[20]; 
  favg[21] = -1*fr[21]+fl[21]; 
  favg[22] = 1*fr[22]+fl[22]; 
  favg[23] = 1*fr[23]+fl[23]; 
  favg[24] = 1*fr[24]+fl[24]; 
  favg[25] = -1*fr[25]+fl[25]; 
  favg[26] = -1*fr[26]+fl[26]; 
  favg[27] = -1*fr[27]+fl[27]; 
  favg[28] = -1*fr[28]+fl[28]; 
  favg[29] = -1*fr[29]+fl[29]; 
  favg[30] = -1*fr[30]+fl[30]; 
  favg[31] = -1*fr[31]+fl[31]; 
  double fjump[32]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(-1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(1*fr[9]-fl[9]); 
  fjump[10] = amax*(-1*fr[10]-fl[10]); 
  fjump[11] = amax*(1*fr[11]-fl[11]); 
  fjump[12] = amax*(1*fr[12]-fl[12]); 
  fjump[13] = amax*(-1*fr[13]-fl[13]); 
  fjump[14] = amax*(-1*fr[14]-fl[14]); 
  fjump[15] = amax*(1*fr[15]-fl[15]); 
  fjump[16] = amax*(1*fr[16]-fl[16]); 
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(1*fr[18]-fl[18]); 
  fjump[19] = amax*(-1*fr[19]-fl[19]); 
  fjump[20] = amax*(-1*fr[20]-fl[20]); 
  fjump[21] = amax*(-1*fr[21]-fl[21]); 
  fjump[22] = amax*(1*fr[22]-fl[22]); 
  fjump[23] = amax*(1*fr[23]-fl[23]); 
  fjump[24] = amax*(1*fr[24]-fl[24]); 
  fjump[25] = amax*(-1*fr[25]-fl[25]); 
  fjump[26] = amax*(-1*fr[26]-fl[26]); 
  fjump[27] = amax*(-1*fr[27]-fl[27]); 
  fjump[28] = amax*(-1*fr[28]-fl[28]); 
  fjump[29] = amax*(-1*fr[29]-fl[29]); 
  fjump[30] = amax*(-1*fr[30]-fl[30]); 
  fjump[31] = amax*(-1*fr[31]-fl[31]); 
  double alpha[32]; 

  alpha[0] = 2.0*E1[0]-2.0*B2[0]*wv1; 
  alpha[1] = 2.0*E1[1]-2.0*B2[1]*wv1; 
  alpha[2] = -0.5773502691896258*B2[0]*dv1; 
  alpha[4] = -0.5773502691896258*B2[1]*dv1; 
  alpha[7] = 2.0*E1[2]-2.0*B2[2]*wv1; 
  alpha[11] = -0.5773502691896257*B2[2]*dv1; 
  alpha[17] = 2.0*E1[3]-2.0*B2[3]*wv1; 
  alpha[23] = -0.5773502691896256*B2[3]*dv1; 
  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 
  Ghat[0] += 0.4677071733467425*alpha[4]*favg[31]+0.3061862178478971*alpha[23]*favg[29]+0.4677071733467425*(alpha[2]*favg[28]+alpha[1]*favg[27])+0.3061862178478971*alpha[17]*favg[25]+0.1767766952966368*alpha[23]*favg[23]+0.3952847075210473*alpha[4]*favg[22]+0.3061862178478971*alpha[11]*favg[20]-1.322875655532295*fjump[19]+0.4677071733467425*alpha[0]*favg[19]+0.1767766952966368*alpha[17]*favg[17]+0.3952847075210473*(alpha[2]*favg[16]+alpha[1]*favg[15])+0.3061862178478971*alpha[7]*favg[13]+0.1767766952966368*alpha[11]*favg[11]+0.3061862178478971*alpha[4]*favg[10]-1.118033988749895*fjump[9]+0.3952847075210473*alpha[0]*favg[9]+0.1767766952966368*alpha[7]*favg[7]+0.3061862178478971*(alpha[2]*favg[6]+alpha[1]*favg[5])+0.1767766952966368*alpha[4]*favg[4]-0.8660254037844386*fjump[3]+0.3061862178478971*alpha[0]*favg[3]+0.1767766952966368*(alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.1767766952966368*alpha[0]*favg[0]; 
  Ghat[1] += (0.4183300132670377*alpha[11]+0.4677071733467425*alpha[2])*favg[31]+0.2689264371002384*alpha[11]*favg[29]+0.4677071733467425*alpha[4]*favg[28]-1.322875655532295*fjump[27]+(0.4183300132670377*alpha[7]+0.4677071733467425*alpha[0])*favg[27]+0.2689264371002384*alpha[7]*favg[25]+0.1552647508520296*alpha[11]*favg[23]+(0.2689264371002384*favg[20]+0.1552647508520296*favg[11])*alpha[23]+(0.3535533905932737*alpha[11]+0.3952847075210473*alpha[2])*favg[22]+0.273861278752583*alpha[4]*favg[20]+0.4677071733467425*alpha[1]*favg[19]+0.1552647508520296*alpha[7]*favg[17]+(0.2689264371002384*favg[13]+0.1552647508520296*favg[7])*alpha[17]+0.3952847075210473*alpha[4]*favg[16]-1.118033988749895*fjump[15]+(0.3535533905932737*alpha[7]+0.3952847075210473*alpha[0])*favg[15]+0.273861278752583*alpha[1]*favg[13]+0.1581138830084189*alpha[4]*favg[11]+(0.273861278752583*favg[10]+0.1581138830084189*favg[4])*alpha[11]+0.3061862178478971*alpha[2]*favg[10]+alpha[1]*(0.3952847075210473*favg[9]+0.1581138830084189*favg[7])+(0.273861278752583*favg[5]+0.1581138830084189*favg[1])*alpha[7]+0.3061862178478971*alpha[4]*favg[6]-0.8660254037844386*fjump[5]+0.3061862178478971*alpha[0]*favg[5]+0.1767766952966368*(alpha[2]*favg[4]+favg[2]*alpha[4])+0.3061862178478971*alpha[1]*favg[3]-0.5*fjump[1]+0.1767766952966368*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += 0.4677071733467425*alpha[1]*favg[31]+0.3061862178478971*alpha[17]*favg[29]-1.322875655532295*fjump[28]+0.4677071733467425*(alpha[0]*favg[28]+alpha[4]*favg[27])+0.3061862178478971*alpha[23]*favg[25]+0.1767766952966368*(alpha[17]*favg[23]+favg[17]*alpha[23])+0.3952847075210473*alpha[1]*favg[22]+0.273861278752583*alpha[4]*favg[21]+0.3061862178478971*alpha[7]*favg[20]+0.4677071733467425*alpha[2]*favg[19]-1.118033988749895*fjump[16]+0.3952847075210473*(alpha[0]*favg[16]+alpha[4]*favg[15])+0.273861278752583*alpha[2]*favg[14]+0.3061862178478971*alpha[11]*favg[13]+0.1581138830084189*alpha[4]*favg[12]+0.1767766952966368*(alpha[7]*favg[11]+favg[7]*alpha[11])+0.3061862178478971*alpha[1]*favg[10]+alpha[2]*(0.3952847075210473*favg[9]+0.1581138830084189*favg[8])-0.8660254037844386*fjump[6]+0.3061862178478971*(alpha[0]*favg[6]+alpha[4]*favg[5])+0.1767766952966368*(alpha[1]*favg[4]+favg[1]*alpha[4])+0.3061862178478971*alpha[2]*favg[3]-0.5*fjump[2]+0.1767766952966368*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[4] += (-1.322875655532295*fjump[31])+(0.4183300132670377*alpha[7]+0.4677071733467425*alpha[0])*favg[31]+0.2689264371002384*alpha[7]*favg[29]+0.4677071733467425*alpha[1]*favg[28]+(0.4183300132670377*alpha[11]+0.4677071733467425*alpha[2])*favg[27]+0.2689264371002384*alpha[11]*favg[25]+0.1552647508520296*alpha[7]*favg[23]+(0.2689264371002384*favg[13]+0.1552647508520296*favg[7])*alpha[23]-1.118033988749895*fjump[22]+(0.3535533905932737*alpha[7]+0.3952847075210473*alpha[0])*favg[22]+(0.2449489742783178*alpha[11]+0.273861278752583*alpha[2])*favg[21]+(0.2689264371002384*alpha[17]+0.273861278752583*alpha[1])*favg[20]+0.4677071733467425*alpha[4]*favg[19]+0.1552647508520296*(alpha[11]*favg[17]+favg[11]*alpha[17])+0.3952847075210473*alpha[1]*favg[16]+(0.3535533905932737*alpha[11]+0.3952847075210473*alpha[2])*favg[15]+0.273861278752583*alpha[4]*(favg[14]+favg[13])+0.1414213562373095*alpha[11]*favg[12]+0.1581138830084189*(alpha[2]*favg[12]+alpha[1]*favg[11])+(0.273861278752583*favg[5]+0.1581138830084189*favg[1])*alpha[11]-0.8660254037844386*fjump[10]+(0.273861278752583*alpha[7]+0.3061862178478971*alpha[0])*favg[10]+0.3952847075210473*alpha[4]*favg[9]+0.1581138830084189*(alpha[4]*(favg[8]+favg[7])+favg[4]*alpha[7])+0.3061862178478971*(alpha[1]*favg[6]+alpha[2]*favg[5])-0.5*fjump[4]+0.1767766952966368*alpha[0]*favg[4]+0.3061862178478971*favg[3]*alpha[4]+0.1767766952966368*(favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[7] += (0.4107919181288743*alpha[23]+0.4183300132670377*alpha[4])*favg[31]+(0.1825741858350554*alpha[23]+0.2689264371002384*alpha[4])*favg[29]+0.4677071733467425*alpha[11]*favg[28]+(0.4107919181288743*alpha[17]+0.4183300132670377*alpha[1])*favg[27]+(0.1825741858350554*alpha[17]+0.2689264371002384*alpha[1])*favg[25]+(0.105409255338946*alpha[23]+0.1552647508520296*alpha[4])*favg[23]+(0.3471825374147066*favg[22]+0.2689264371002384*favg[10]+0.1552647508520296*favg[4])*alpha[23]+0.3535533905932737*alpha[4]*favg[22]+(0.1956151991089878*alpha[11]+0.3061862178478971*alpha[2])*favg[20]+0.4677071733467425*alpha[7]*favg[19]+(0.105409255338946*alpha[17]+0.1552647508520296*alpha[1])*favg[17]+(0.3471825374147066*favg[15]+0.2689264371002384*favg[5]+0.1552647508520296*favg[1])*alpha[17]+0.3952847075210473*alpha[11]*favg[16]+0.3535533905932737*alpha[1]*favg[15]-0.8660254037844386*fjump[13]+(0.1956151991089878*alpha[7]+0.3061862178478971*alpha[0])*favg[13]+(0.1129384878631564*alpha[11]+0.1767766952966368*alpha[2])*favg[11]+(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])*alpha[11]+0.273861278752583*alpha[4]*favg[10]+0.3952847075210473*alpha[7]*favg[9]-0.5*fjump[7]+(0.1129384878631564*alpha[7]+0.1767766952966368*alpha[0])*favg[7]+(0.3061862178478971*favg[3]+0.1767766952966368*favg[0])*alpha[7]+0.273861278752583*alpha[1]*favg[5]+0.1581138830084189*(alpha[4]*favg[4]+alpha[1]*favg[1]); 
  Ghat[8] += alpha[4]*(0.4183300132670377*favg[31]+0.2689264371002384*favg[30])+0.273861278752583*alpha[23]*favg[29]+alpha[2]*(0.4183300132670377*favg[28]+0.2689264371002384*favg[26])+0.1552647508520296*alpha[4]*favg[24]+0.1581138830084189*alpha[23]*favg[23]+0.3535533905932737*alpha[4]*favg[22]+0.3061862178478971*alpha[1]*favg[21]+0.273861278752583*alpha[11]*favg[20]+alpha[2]*(0.1552647508520296*favg[18]+0.3535533905932737*favg[16])-0.8660254037844386*fjump[14]+0.3061862178478971*alpha[0]*favg[14]+0.1767766952966368*alpha[1]*favg[12]+0.1581138830084189*alpha[11]*favg[11]+0.273861278752583*alpha[4]*favg[10]-0.5*fjump[8]+0.1767766952966368*alpha[0]*favg[8]+0.273861278752583*alpha[2]*favg[6]+0.1581138830084189*(alpha[4]*favg[4]+alpha[2]*favg[2]); 
  Ghat[11] += (0.4107919181288743*alpha[17]+0.4183300132670377*alpha[1])*favg[31]+(0.1825741858350554*alpha[17]+0.2689264371002384*alpha[1])*favg[29]+0.4677071733467425*alpha[7]*favg[28]+(0.4107919181288743*alpha[23]+0.4183300132670377*alpha[4])*favg[27]+(0.1825741858350554*alpha[23]+0.2689264371002384*alpha[4])*favg[25]+(0.105409255338946*alpha[17]+0.1552647508520296*alpha[1])*favg[23]+(0.2405351177211819*favg[21]+0.105409255338946*favg[17]+0.3471825374147066*favg[15]+0.1388730149658826*favg[12]+0.2689264371002384*favg[5]+0.1552647508520296*favg[1])*alpha[23]+(0.3471825374147066*alpha[17]+0.3535533905932737*alpha[1])*favg[22]+0.2449489742783178*alpha[4]*favg[21]-0.8660254037844386*fjump[20]+(0.1956151991089878*alpha[7]+0.3061862178478971*alpha[0])*favg[20]+0.4677071733467425*alpha[11]*favg[19]+0.1552647508520296*alpha[4]*favg[17]+(0.2689264371002384*favg[10]+0.1552647508520296*favg[4])*alpha[17]+0.3952847075210473*alpha[7]*favg[16]+0.3535533905932737*alpha[4]*favg[15]+0.273861278752583*alpha[11]*favg[14]+(0.1956151991089878*alpha[11]+0.3061862178478971*alpha[2])*favg[13]+0.1414213562373095*alpha[4]*favg[12]-0.5*fjump[11]+(0.1129384878631564*alpha[7]+0.1767766952966368*alpha[0])*favg[11]+(0.3952847075210473*favg[9]+0.1581138830084189*favg[8]+0.1129384878631564*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])*alpha[11]+0.273861278752583*alpha[1]*favg[10]+0.1767766952966368*alpha[2]*favg[7]+(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])*alpha[7]+0.273861278752583*alpha[4]*favg[5]+0.1581138830084189*(alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[12] += (0.3741657386773942*alpha[11]+0.4183300132670377*alpha[2])*favg[31]+(0.2405351177211819*alpha[11]+0.2689264371002384*alpha[2])*favg[30]+0.2405351177211819*alpha[11]*favg[29]+alpha[4]*(0.4183300132670377*favg[28]+0.2689264371002384*favg[26])+(0.1388730149658826*alpha[11]+0.1552647508520296*alpha[2])*favg[24]+0.1388730149658826*alpha[11]*favg[23]+(0.2405351177211819*favg[20]+0.1388730149658826*favg[11])*alpha[23]+(0.3162277660168379*alpha[11]+0.3535533905932737*alpha[2])*favg[22]-0.8660254037844386*fjump[21]+(0.273861278752583*alpha[7]+0.3061862178478971*alpha[0])*favg[21]+alpha[4]*(0.2449489742783178*favg[20]+0.1552647508520296*favg[18]+0.3535533905932737*favg[16])+0.3061862178478971*alpha[1]*favg[14]-0.5*fjump[12]+(0.1581138830084189*alpha[7]+0.1767766952966368*alpha[0])*favg[12]+0.1414213562373095*alpha[4]*favg[11]+(0.2449489742783178*favg[10]+0.1414213562373095*favg[4])*alpha[11]+0.273861278752583*alpha[2]*favg[10]+0.1767766952966368*alpha[1]*favg[8]+0.273861278752583*alpha[4]*favg[6]+0.1581138830084189*(alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[17] += 0.4107919181288743*alpha[11]*favg[31]+(0.1825741858350554*alpha[11]+0.3061862178478971*alpha[2])*favg[29]+0.4677071733467425*alpha[23]*favg[28]+0.4107919181288743*alpha[7]*favg[27]-0.8660254037844386*fjump[25]+(0.1825741858350554*alpha[7]+0.3061862178478971*alpha[0])*favg[25]+(0.105409255338946*alpha[11]+0.1767766952966368*alpha[2])*favg[23]+(0.1825741858350554*favg[20]+0.3952847075210473*favg[16]+0.105409255338946*favg[11]+0.3061862178478971*favg[6]+0.1767766952966368*favg[2])*alpha[23]+0.3471825374147066*alpha[11]*favg[22]+0.2689264371002384*alpha[4]*favg[20]+0.4677071733467425*alpha[17]*favg[19]-0.5*fjump[17]+(0.105409255338946*alpha[7]+0.1767766952966368*alpha[0])*favg[17]+(0.1825741858350554*favg[13]+0.3952847075210473*favg[9]+0.105409255338946*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])*alpha[17]+0.3471825374147066*alpha[7]*favg[15]+0.2689264371002384*alpha[1]*favg[13]+0.1552647508520296*alpha[4]*favg[11]+0.2689264371002384*favg[10]*alpha[11]+0.1552647508520296*(favg[4]*alpha[11]+alpha[1]*favg[7])+(0.2689264371002384*favg[5]+0.1552647508520296*favg[1])*alpha[7]; 
  Ghat[18] += 0.3061862178478971*alpha[1]*favg[30]-0.8660254037844386*fjump[26]+0.3061862178478971*alpha[0]*favg[26]+0.1767766952966368*alpha[1]*favg[24]+0.2689264371002384*alpha[4]*favg[21]-0.5*fjump[18]+0.1767766952966368*alpha[0]*favg[18]+0.2689264371002384*alpha[2]*favg[14]+0.1552647508520296*(alpha[4]*favg[12]+alpha[2]*favg[8]); 
  Ghat[23] += 0.4107919181288743*alpha[7]*favg[31]-0.8660254037844386*fjump[29]+(0.1825741858350554*alpha[7]+0.3061862178478971*alpha[0])*favg[29]+0.4677071733467425*alpha[17]*favg[28]+0.4107919181288743*alpha[11]*favg[27]+(0.1825741858350554*alpha[11]+0.3061862178478971*alpha[2])*favg[25]-0.5*fjump[23]+(0.105409255338946*alpha[7]+0.1767766952966368*alpha[0])*favg[23]+(0.4677071733467425*favg[19]+0.273861278752583*favg[14]+0.1825741858350554*favg[13]+0.3952847075210473*favg[9]+0.1581138830084189*favg[8]+0.105409255338946*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])*alpha[23]+0.3471825374147066*alpha[7]*favg[22]+0.2405351177211819*alpha[11]*favg[21]+(0.1825741858350554*alpha[17]+0.2689264371002384*alpha[1])*favg[20]+(0.105409255338946*alpha[11]+0.1767766952966368*alpha[2])*favg[17]+(0.3952847075210473*favg[16]+0.105409255338946*favg[11]+0.3061862178478971*favg[6]+0.1767766952966368*favg[2])*alpha[17]+0.3471825374147066*alpha[11]*favg[15]+0.2689264371002384*alpha[4]*favg[13]+0.1388730149658826*alpha[11]*favg[12]+0.1552647508520296*alpha[1]*favg[11]+(0.2689264371002384*favg[5]+0.1552647508520296*favg[1])*alpha[11]+0.2689264371002384*alpha[7]*favg[10]+0.1552647508520296*(alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[24] += (-0.8660254037844386*fjump[30])+0.273861278752583*alpha[7]*favg[30]+0.3061862178478971*(alpha[0]*favg[30]+alpha[1]*favg[26])-0.5*fjump[24]+(0.1581138830084189*alpha[7]+0.1767766952966368*alpha[0])*favg[24]+(0.2405351177211819*alpha[11]+0.2689264371002384*alpha[2])*favg[21]+0.1767766952966368*alpha[1]*favg[18]+0.2689264371002384*alpha[4]*favg[14]+0.1388730149658826*alpha[11]*favg[12]+0.1552647508520296*(alpha[2]*favg[12]+alpha[4]*favg[8]); 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += -0.8660254037844386*Ghat[0]*dv11r; 
  outr[4] += 0.5*Ghat[4]*dv11r; 
  outr[5] += -0.8660254037844386*Ghat[1]*dv11r; 
  outr[6] += -0.8660254037844386*Ghat[2]*dv11r; 
  outr[7] += 0.5*Ghat[7]*dv11r; 
  outr[8] += 0.5*Ghat[8]*dv11r; 
  outr[9] += 1.118033988749895*Ghat[0]*dv11r; 
  outr[10] += -0.8660254037844386*Ghat[4]*dv11r; 
  outr[11] += 0.5*Ghat[11]*dv11r; 
  outr[12] += 0.5*Ghat[12]*dv11r; 
  outr[13] += -0.8660254037844387*Ghat[7]*dv11r; 
  outr[14] += -0.8660254037844387*Ghat[8]*dv11r; 
  outr[15] += 1.118033988749895*Ghat[1]*dv11r; 
  outr[16] += 1.118033988749895*Ghat[2]*dv11r; 
  outr[17] += 0.5*Ghat[17]*dv11r; 
  outr[18] += 0.5*Ghat[18]*dv11r; 
  outr[19] += -1.322875655532295*Ghat[0]*dv11r; 
  outr[20] += -0.8660254037844387*Ghat[11]*dv11r; 
  outr[21] += -0.8660254037844387*Ghat[12]*dv11r; 
  outr[22] += 1.118033988749895*Ghat[4]*dv11r; 
  outr[23] += 0.5*Ghat[23]*dv11r; 
  outr[24] += 0.5*Ghat[24]*dv11r; 
  outr[25] += -0.8660254037844386*Ghat[17]*dv11r; 
  outr[26] += -0.8660254037844386*Ghat[18]*dv11r; 
  outr[27] += -1.322875655532295*Ghat[1]*dv11r; 
  outr[28] += -1.322875655532295*Ghat[2]*dv11r; 
  outr[29] += -0.8660254037844386*Ghat[23]*dv11r; 
  outr[30] += -0.8660254037844386*Ghat[24]*dv11r; 
  outr[31] += -1.322875655532295*Ghat[4]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.8660254037844386*Ghat[0]*dv11l; 
  outl[4] += -0.5*Ghat[4]*dv11l; 
  outl[5] += -0.8660254037844386*Ghat[1]*dv11l; 
  outl[6] += -0.8660254037844386*Ghat[2]*dv11l; 
  outl[7] += -0.5*Ghat[7]*dv11l; 
  outl[8] += -0.5*Ghat[8]*dv11l; 
  outl[9] += -1.118033988749895*Ghat[0]*dv11l; 
  outl[10] += -0.8660254037844386*Ghat[4]*dv11l; 
  outl[11] += -0.5*Ghat[11]*dv11l; 
  outl[12] += -0.5*Ghat[12]*dv11l; 
  outl[13] += -0.8660254037844387*Ghat[7]*dv11l; 
  outl[14] += -0.8660254037844387*Ghat[8]*dv11l; 
  outl[15] += -1.118033988749895*Ghat[1]*dv11l; 
  outl[16] += -1.118033988749895*Ghat[2]*dv11l; 
  outl[17] += -0.5*Ghat[17]*dv11l; 
  outl[18] += -0.5*Ghat[18]*dv11l; 
  outl[19] += -1.322875655532295*Ghat[0]*dv11l; 
  outl[20] += -0.8660254037844387*Ghat[11]*dv11l; 
  outl[21] += -0.8660254037844387*Ghat[12]*dv11l; 
  outl[22] += -1.118033988749895*Ghat[4]*dv11l; 
  outl[23] += -0.5*Ghat[23]*dv11l; 
  outl[24] += -0.5*Ghat[24]*dv11l; 
  outl[25] += -0.8660254037844386*Ghat[17]*dv11l; 
  outl[26] += -0.8660254037844386*Ghat[18]*dv11l; 
  outl[27] += -1.322875655532295*Ghat[1]*dv11l; 
  outl[28] += -1.322875655532295*Ghat[2]*dv11l; 
  outl[29] += -0.8660254037844386*Ghat[23]*dv11l; 
  outl[30] += -0.8660254037844386*Ghat[24]*dv11l; 
  outl[31] += -1.322875655532295*Ghat[4]*dv11l; 
return std::abs(amid); 
} 
