#include <VlasovModDecl.h> 
double VlasovSurfElcMag1x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  double Ghat[4]; 

  double alpha[4]; 

  for(unsigned int i=0; i<4; ++i){ 

    Ghat[i]=0.0; 

    alpha[i]=0.0; 

  }; 

  double favg[4]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  double fjump[4]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  alpha[0] += 2.0*(B2[0]*wv2+E0[0]); 
  alpha[1] += 2.0*(B2[1]*wv2+E0[1]); 
  alpha[3] += 0.5773502691896258*B2[0]*dv2; 
  const double amid = 0.3535533905932737*alpha[0]; 
  Ghat[0] = 0.1767766952966368*alpha[3]*favg[3]-0.8660254037844386*fjump[2]+alpha[0]*(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+0.1767766952966368*alpha[1]*favg[1]-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.5*fjump[1]+0.1767766952966368*alpha[0]*favg[1]; 
  Ghat[3] = (-0.5*fjump[3])+0.1767766952966368*alpha[0]*favg[3]+(0.3061862178478971*favg[2]+0.1767766952966368*favg[0])*alpha[3]; 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[3] += 0.5*Ghat[3]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[3] += -0.5*Ghat[3]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  double Ghat[10]; 

  double alpha[10]; 

  for(unsigned int i=0; i<10; ++i){ 

    Ghat[i]=0.0; 

    alpha[i]=0.0; 

  }; 

  double favg[10]; 

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
  double fjump[10]; 

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
  alpha[0] += 2.0*(B2[0]*wv2+E0[0]); 
  alpha[1] += 2.0*(B2[1]*wv2+E0[1]); 
  alpha[3] += 0.5773502691896258*B2[0]*dv2; 
  alpha[5] += 0.5773502691896258*B2[1]*dv2; 
  alpha[7] += 2.0*(B2[2]*wv2+E0[2]); 
  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 
  Ghat[0] = (-1.118033988749895*fjump[8])+alpha[0]*(0.3952847075210473*favg[8]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+0.1767766952966368*alpha[7]*favg[7]+alpha[3]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+0.1767766952966368*alpha[5]*favg[5]+alpha[1]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[1])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.3952847075210473*favg[8]+0.1581138830084189*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+(0.273861278752583*favg[4]+0.1581138830084189*favg[1])*alpha[7]+alpha[5]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+0.1767766952966368*alpha[3]*favg[5]-0.8660254037844386*fjump[4]+alpha[0]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[1])-0.5*fjump[1]; 
  Ghat[3] = alpha[3]*(0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[6]+alpha[0]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+0.1767766952966368*alpha[1]*favg[5]+(0.3061862178478971*favg[4]+0.1767766952966368*favg[1])*alpha[5]-0.5*fjump[3]; 
  Ghat[5] = alpha[5]*(0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.1581138830084189*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+0.1581138830084189*favg[5]*alpha[7]+alpha[1]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[3])-0.5*fjump[5]+0.1767766952966368*alpha[0]*favg[5]+alpha[3]*(0.3061862178478971*favg[4]+0.1767766952966368*favg[1]); 
  Ghat[7] = alpha[7]*(0.3952847075210473*favg[8]+0.1129384878631564*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.5*fjump[7]+0.1767766952966368*alpha[0]*favg[7]+0.1581138830084189*alpha[5]*favg[5]+alpha[1]*(0.273861278752583*favg[4]+0.1581138830084189*favg[1]); 
  Ghat[9] = (-0.5*fjump[9])+0.1767766952966368*alpha[0]*favg[9]+alpha[3]*(0.273861278752583*favg[6]+0.1581138830084189*favg[3])+0.1581138830084189*alpha[5]*favg[5]; 

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
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  double Ghat[20]; 

  double alpha[20]; 

  for(unsigned int i=0; i<20; ++i){ 

    Ghat[i]=0.0; 

    alpha[i]=0.0; 

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
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = -1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
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
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(-1*fr[18]-fl[18]); 
  fjump[19] = amax*(1*fr[19]-fl[19]); 
  alpha[0] += 2.0*(B2[0]*wv2+E0[0]); 
  alpha[1] += 2.0*(B2[1]*wv2+E0[1]); 
  alpha[3] += 0.5773502691896258*B2[0]*dv2; 
  alpha[5] += 0.5773502691896258*B2[1]*dv2; 
  alpha[7] += 2.0*(B2[2]*wv2+E0[2]); 
  alpha[13] += 0.5773502691896258*B2[2]*dv2; 
  alpha[17] += 2.0*(B2[3]*wv2+E0[3]); 
  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 
  Ghat[0] = (-1.322875655532295*fjump[18])+alpha[0]*(0.4677071733467425*favg[18]+0.3952847075210473*favg[8]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+0.1767766952966368*alpha[17]*favg[17]+alpha[3]*(0.3952847075210473*favg[14]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+0.1767766952966368*alpha[13]*favg[13]+alpha[1]*(0.3952847075210473*favg[12]+0.3061862178478971*favg[4]+0.1767766952966368*favg[1])+alpha[7]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[7])+alpha[5]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[5])-1.118033988749895*fjump[8]-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.4677071733467425*favg[18]+0.273861278752583*favg[11]+0.3952847075210473*favg[8]+0.1581138830084189*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+alpha[7]*(0.1552647508520296*favg[17]+0.3535533905932737*favg[12]+0.273861278752583*favg[4]+0.1581138830084189*favg[1])+(0.2689264371002384*favg[11]+0.1552647508520296*favg[7])*alpha[17]+alpha[5]*(0.3952847075210473*favg[14]+0.1581138830084189*favg[13]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+(0.273861278752583*favg[10]+0.1581138830084189*favg[5])*alpha[13]-1.118033988749895*fjump[12]+alpha[0]*(0.3952847075210473*favg[12]+0.3061862178478971*favg[4]+0.1767766952966368*favg[1])+alpha[3]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[5])-0.8660254037844386*fjump[4]-0.5*fjump[1]; 
  Ghat[3] = alpha[3]*(0.4677071733467425*favg[18]+0.273861278752583*favg[16]+0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+alpha[5]*(0.1581138830084189*favg[15]+0.3952847075210473*favg[12]+0.3061862178478971*favg[4]+0.1767766952966368*favg[1])-1.118033988749895*fjump[14]+alpha[0]*(0.3952847075210473*favg[14]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+0.1767766952966368*alpha[7]*favg[13]+(0.3061862178478971*favg[11]+0.1767766952966368*favg[7])*alpha[13]+alpha[1]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[5])-0.8660254037844386*fjump[6]-0.5*fjump[3]; 
  Ghat[5] = alpha[5]*(0.4677071733467425*favg[18]+0.273861278752583*(favg[16]+favg[11])+0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.1581138830084189*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+alpha[13]*(0.1552647508520296*favg[17]+0.1414213562373095*favg[15]+0.3535533905932737*favg[12]+0.273861278752583*favg[4]+0.1581138830084189*favg[1])+0.1552647508520296*favg[13]*alpha[17]+alpha[3]*(0.1581138830084189*favg[15]+0.3952847075210473*favg[12]+0.3061862178478971*favg[4]+0.1767766952966368*favg[1])+alpha[1]*(0.3952847075210473*favg[14]+0.1581138830084189*favg[13]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])-0.8660254037844386*fjump[10]+alpha[0]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[5])+alpha[7]*(0.273861278752583*favg[10]+0.1581138830084189*favg[5])-0.5*fjump[5]; 
  Ghat[7] = alpha[7]*(0.4677071733467425*favg[18]+0.1956151991089878*favg[11]+0.3952847075210473*favg[8]+0.1129384878631564*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+alpha[1]*(0.1552647508520296*favg[17]+0.3535533905932737*favg[12]+0.273861278752583*favg[4]+0.1581138830084189*favg[1])+alpha[17]*(0.105409255338946*favg[17]+0.3471825374147066*favg[12]+0.2689264371002384*favg[4]+0.1552647508520296*favg[1])+alpha[13]*(0.3952847075210473*favg[14]+0.1129384878631564*favg[13]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+0.1767766952966368*alpha[3]*favg[13]-0.8660254037844386*fjump[11]+alpha[0]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[7])+alpha[5]*(0.273861278752583*favg[10]+0.1581138830084189*favg[5])-0.5*fjump[7]; 
  Ghat[9] = alpha[3]*(0.1552647508520296*favg[19]+0.3535533905932737*favg[14]+0.273861278752583*favg[6]+0.1581138830084189*favg[3])-0.8660254037844386*fjump[16]+alpha[0]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[9])+0.1767766952966368*alpha[1]*favg[15]+0.1581138830084189*alpha[13]*favg[13]+alpha[5]*(0.273861278752583*favg[10]+0.1581138830084189*favg[5])-0.5*fjump[9]; 
  Ghat[13] = alpha[13]*(0.4677071733467425*favg[18]+0.273861278752583*favg[16]+0.1956151991089878*favg[11]+0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.1129384878631564*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])+alpha[5]*(0.1552647508520296*favg[17]+0.1414213562373095*favg[15]+0.3535533905932737*favg[12]+0.273861278752583*favg[4]+0.1581138830084189*favg[1])+(0.2689264371002384*favg[10]+0.1552647508520296*favg[5])*alpha[17]+alpha[7]*(0.3952847075210473*favg[14]+0.1129384878631564*favg[13]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])-0.5*fjump[13]+0.1767766952966368*alpha[0]*favg[13]+alpha[3]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[7])+alpha[1]*(0.273861278752583*favg[10]+0.1581138830084189*favg[5]); 
  Ghat[15] = alpha[5]*(0.1552647508520296*favg[19]+0.3535533905932737*favg[14]+0.1414213562373095*favg[13]+0.273861278752583*favg[6]+0.1581138830084189*favg[3])+alpha[1]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[9])-0.5*fjump[15]+(0.1581138830084189*alpha[7]+0.1767766952966368*alpha[0])*favg[15]+(0.2449489742783178*favg[10]+0.1414213562373095*favg[5])*alpha[13]+alpha[3]*(0.273861278752583*favg[10]+0.1581138830084189*favg[5]); 
  Ghat[17] = alpha[17]*(0.4677071733467425*favg[18]+0.1825741858350554*favg[11]+0.3952847075210473*favg[8]+0.105409255338946*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.5*fjump[17]+0.1767766952966368*alpha[0]*favg[17]+alpha[7]*(0.105409255338946*favg[17]+0.3471825374147066*favg[12]+0.2689264371002384*favg[4]+0.1552647508520296*favg[1])+0.1552647508520296*alpha[5]*favg[13]+(0.2689264371002384*favg[10]+0.1552647508520296*favg[5])*alpha[13]+alpha[1]*(0.2689264371002384*favg[11]+0.1552647508520296*favg[7]); 
  Ghat[19] = (-0.5*fjump[19])+0.1767766952966368*alpha[0]*favg[19]+alpha[3]*(0.2689264371002384*favg[16]+0.1552647508520296*favg[9])+0.1552647508520296*alpha[5]*favg[15]; 

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
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  double Ghat[4]; 

  double alpha[4]; 

  for(unsigned int i=0; i<4; ++i){ 

    Ghat[i]=0.0; 

    alpha[i]=0.0; 

  }; 

  double favg[4]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  double fjump[4]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  alpha[0] += 2.0*E1[0]-2.0*B2[0]*wv1; 
  alpha[1] += 2.0*E1[1]-2.0*B2[1]*wv1; 
  alpha[2] += -0.5773502691896258*B2[0]*dv1; 
  const double amid = 0.3535533905932737*alpha[0]; 
  Ghat[0] = (-0.8660254037844386*fjump[3])+alpha[0]*(0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+0.1767766952966368*(alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.3061862178478971*favg[3]+0.1767766952966368*favg[0])-0.5*fjump[1]+0.1767766952966368*alpha[0]*favg[1]; 
  Ghat[2] = alpha[2]*(0.3061862178478971*favg[3]+0.1767766952966368*favg[0])-0.5*fjump[2]+0.1767766952966368*alpha[0]*favg[2]; 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += -0.8660254037844386*Ghat[0]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.8660254037844386*Ghat[0]*dv11l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  double Ghat[10]; 

  double alpha[10]; 

  for(unsigned int i=0; i<10; ++i){ 

    Ghat[i]=0.0; 

    alpha[i]=0.0; 

  }; 

  double favg[10]; 

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
  double fjump[10]; 

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
  alpha[0] += 2.0*E1[0]-2.0*B2[0]*wv1; 
  alpha[1] += 2.0*E1[1]-2.0*B2[1]*wv1; 
  alpha[2] += -0.5773502691896258*B2[0]*dv1; 
  alpha[4] += -0.5773502691896258*B2[1]*dv1; 
  alpha[7] += 2.0*E1[2]-2.0*B2[2]*wv1; 
  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 
  Ghat[0] = (-1.118033988749895*fjump[9])+alpha[0]*(0.3952847075210473*favg[9]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+0.1767766952966368*alpha[7]*favg[7]+alpha[2]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])+alpha[1]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[1])+0.1767766952966368*alpha[4]*favg[4]-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.3952847075210473*favg[9]+0.1581138830084189*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+(0.273861278752583*favg[5]+0.1581138830084189*favg[1])*alpha[7]+alpha[4]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])-0.8660254037844386*fjump[5]+alpha[0]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[1])+0.1767766952966368*alpha[2]*favg[4]-0.5*fjump[1]; 
  Ghat[2] = alpha[2]*(0.3952847075210473*favg[9]+0.1581138830084189*favg[8]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[6]+alpha[0]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])+alpha[4]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[1])+0.1767766952966368*alpha[1]*favg[4]-0.5*fjump[2]; 
  Ghat[4] = alpha[4]*(0.3952847075210473*favg[9]+0.1581138830084189*(favg[8]+favg[7])+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+0.1581138830084189*favg[4]*alpha[7]+alpha[1]*(0.3061862178478971*favg[6]+0.1767766952966368*favg[2])+alpha[2]*(0.3061862178478971*favg[5]+0.1767766952966368*favg[1])-0.5*fjump[4]+0.1767766952966368*alpha[0]*favg[4]; 
  Ghat[7] = alpha[7]*(0.3952847075210473*favg[9]+0.1129384878631564*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])-0.5*fjump[7]+0.1767766952966368*alpha[0]*favg[7]+alpha[1]*(0.273861278752583*favg[5]+0.1581138830084189*favg[1])+0.1581138830084189*alpha[4]*favg[4]; 
  Ghat[8] = (-0.5*fjump[8])+0.1767766952966368*alpha[0]*favg[8]+alpha[2]*(0.273861278752583*favg[6]+0.1581138830084189*favg[2])+0.1581138830084189*alpha[4]*favg[4]; 

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
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  double Ghat[20]; 

  double alpha[20]; 

  for(unsigned int i=0; i<20; ++i){ 

    Ghat[i]=0.0; 

    alpha[i]=0.0; 

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
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = -1*fr[19]+fl[19]; 
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
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(1*fr[18]-fl[18]); 
  fjump[19] = amax*(-1*fr[19]-fl[19]); 
  alpha[0] += 2.0*E1[0]-2.0*B2[0]*wv1; 
  alpha[1] += 2.0*E1[1]-2.0*B2[1]*wv1; 
  alpha[2] += -0.5773502691896258*B2[0]*dv1; 
  alpha[4] += -0.5773502691896258*B2[1]*dv1; 
  alpha[7] += 2.0*E1[2]-2.0*B2[2]*wv1; 
  alpha[11] += -0.5773502691896258*B2[2]*dv1; 
  alpha[17] += 2.0*E1[3]-2.0*B2[3]*wv1; 
  const double amid = 0.3535533905932737*alpha[0]-0.3952847075210473*alpha[7]; 
  Ghat[0] = (-1.322875655532295*fjump[19])+alpha[0]*(0.4677071733467425*favg[19]+0.3952847075210473*favg[9]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+0.1767766952966368*alpha[17]*favg[17]+alpha[2]*(0.3952847075210473*favg[16]+0.3061862178478971*favg[6]+0.1767766952966368*favg[2])+alpha[1]*(0.3952847075210473*favg[15]+0.3061862178478971*favg[5]+0.1767766952966368*favg[1])+alpha[7]*(0.3061862178478971*favg[13]+0.1767766952966368*favg[7])+0.1767766952966368*alpha[11]*favg[11]+alpha[4]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[4])-1.118033988749895*fjump[9]-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = alpha[1]*(0.4677071733467425*favg[19]+0.273861278752583*favg[13]+0.3952847075210473*favg[9]+0.1581138830084189*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+alpha[7]*(0.1552647508520296*favg[17]+0.3535533905932737*favg[15]+0.273861278752583*favg[5]+0.1581138830084189*favg[1])+(0.2689264371002384*favg[13]+0.1552647508520296*favg[7])*alpha[17]+alpha[4]*(0.3952847075210473*favg[16]+0.1581138830084189*favg[11]+0.3061862178478971*favg[6]+0.1767766952966368*favg[2])-1.118033988749895*fjump[15]+alpha[0]*(0.3952847075210473*favg[15]+0.3061862178478971*favg[5]+0.1767766952966368*favg[1])+(0.273861278752583*favg[10]+0.1581138830084189*favg[4])*alpha[11]+alpha[2]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[4])-0.8660254037844386*fjump[5]-0.5*fjump[1]; 
  Ghat[2] = alpha[2]*(0.4677071733467425*favg[19]+0.273861278752583*favg[14]+0.3952847075210473*favg[9]+0.1581138830084189*favg[8]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])-1.118033988749895*fjump[16]+alpha[0]*(0.3952847075210473*favg[16]+0.3061862178478971*favg[6]+0.1767766952966368*favg[2])+alpha[4]*(0.3952847075210473*favg[15]+0.1581138830084189*favg[12]+0.3061862178478971*favg[5]+0.1767766952966368*favg[1])+alpha[11]*(0.3061862178478971*favg[13]+0.1767766952966368*favg[7])+0.1767766952966368*alpha[7]*favg[11]+alpha[1]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[4])-0.8660254037844386*fjump[6]-0.5*fjump[2]; 
  Ghat[4] = alpha[4]*(0.4677071733467425*favg[19]+0.273861278752583*(favg[14]+favg[13])+0.3952847075210473*favg[9]+0.1581138830084189*(favg[8]+favg[7])+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+alpha[11]*(0.1552647508520296*favg[17]+0.3535533905932737*favg[15]+0.1414213562373095*favg[12]+0.273861278752583*favg[5]+0.1581138830084189*favg[1])+0.1552647508520296*favg[11]*alpha[17]+alpha[1]*(0.3952847075210473*favg[16]+0.1581138830084189*favg[11]+0.3061862178478971*favg[6]+0.1767766952966368*favg[2])+alpha[2]*(0.3952847075210473*favg[15]+0.1581138830084189*favg[12]+0.3061862178478971*favg[5]+0.1767766952966368*favg[1])-0.8660254037844386*fjump[10]+alpha[0]*(0.3061862178478971*favg[10]+0.1767766952966368*favg[4])+alpha[7]*(0.273861278752583*favg[10]+0.1581138830084189*favg[4])-0.5*fjump[4]; 
  Ghat[7] = alpha[7]*(0.4677071733467425*favg[19]+0.1956151991089878*favg[13]+0.3952847075210473*favg[9]+0.1129384878631564*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+alpha[1]*(0.1552647508520296*favg[17]+0.3535533905932737*favg[15]+0.273861278752583*favg[5]+0.1581138830084189*favg[1])+alpha[17]*(0.105409255338946*favg[17]+0.3471825374147066*favg[15]+0.2689264371002384*favg[5]+0.1552647508520296*favg[1])+alpha[11]*(0.3952847075210473*favg[16]+0.1129384878631564*favg[11]+0.3061862178478971*favg[6]+0.1767766952966368*favg[2])-0.8660254037844386*fjump[13]+alpha[0]*(0.3061862178478971*favg[13]+0.1767766952966368*favg[7])+0.1767766952966368*alpha[2]*favg[11]+alpha[4]*(0.273861278752583*favg[10]+0.1581138830084189*favg[4])-0.5*fjump[7]; 
  Ghat[8] = alpha[2]*(0.1552647508520296*favg[18]+0.3535533905932737*favg[16]+0.273861278752583*favg[6]+0.1581138830084189*favg[2])-0.8660254037844386*fjump[14]+alpha[0]*(0.3061862178478971*favg[14]+0.1767766952966368*favg[8])+0.1767766952966368*alpha[1]*favg[12]+0.1581138830084189*alpha[11]*favg[11]+alpha[4]*(0.273861278752583*favg[10]+0.1581138830084189*favg[4])-0.5*fjump[8]; 
  Ghat[11] = alpha[11]*(0.4677071733467425*favg[19]+0.273861278752583*favg[14]+0.1956151991089878*favg[13]+0.3952847075210473*favg[9]+0.1581138830084189*favg[8]+0.1129384878631564*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])+alpha[4]*(0.1552647508520296*favg[17]+0.3535533905932737*favg[15]+0.1414213562373095*favg[12]+0.273861278752583*favg[5]+0.1581138830084189*favg[1])+(0.2689264371002384*favg[10]+0.1552647508520296*favg[4])*alpha[17]+alpha[7]*(0.3952847075210473*favg[16]+0.1129384878631564*favg[11]+0.3061862178478971*favg[6]+0.1767766952966368*favg[2])+alpha[2]*(0.3061862178478971*favg[13]+0.1767766952966368*favg[7])-0.5*fjump[11]+0.1767766952966368*alpha[0]*favg[11]+alpha[1]*(0.273861278752583*favg[10]+0.1581138830084189*favg[4]); 
  Ghat[12] = alpha[4]*(0.1552647508520296*favg[18]+0.3535533905932737*favg[16]+0.1414213562373095*favg[11]+0.273861278752583*favg[6]+0.1581138830084189*favg[2])+alpha[1]*(0.3061862178478971*favg[14]+0.1767766952966368*favg[8])-0.5*fjump[12]+(0.1581138830084189*alpha[7]+0.1767766952966368*alpha[0])*favg[12]+(0.2449489742783178*favg[10]+0.1414213562373095*favg[4])*alpha[11]+alpha[2]*(0.273861278752583*favg[10]+0.1581138830084189*favg[4]); 
  Ghat[17] = alpha[17]*(0.4677071733467425*favg[19]+0.1825741858350554*favg[13]+0.3952847075210473*favg[9]+0.105409255338946*favg[7]+0.3061862178478971*favg[3]+0.1767766952966368*favg[0])-0.5*fjump[17]+0.1767766952966368*alpha[0]*favg[17]+alpha[7]*(0.105409255338946*favg[17]+0.3471825374147066*favg[15]+0.2689264371002384*favg[5]+0.1552647508520296*favg[1])+alpha[1]*(0.2689264371002384*favg[13]+0.1552647508520296*favg[7])+0.1552647508520296*alpha[4]*favg[11]+(0.2689264371002384*favg[10]+0.1552647508520296*favg[4])*alpha[11]; 
  Ghat[18] = (-0.5*fjump[18])+0.1767766952966368*alpha[0]*favg[18]+alpha[2]*(0.2689264371002384*favg[14]+0.1552647508520296*favg[8])+0.1552647508520296*alpha[4]*favg[12]; 

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
return std::abs(amid); 
} 
