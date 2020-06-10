//#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag2x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[2]; 
  double dv10r = 2/dxvr[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 
  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double Ghat[32]; 
  double alpha[32]; 

  double favg[32]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = -1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = -1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = -1*fr[18]+fl[18]; 
  favg[19] = -1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  favg[21] = -1*fr[21]+fl[21]; 
  favg[22] = -1*fr[22]+fl[22]; 
  favg[23] = 1*fr[23]+fl[23]; 
  favg[24] = 1*fr[24]+fl[24]; 
  favg[25] = -1*fr[25]+fl[25]; 
  favg[26] = -1*fr[26]+fl[26]; 
  favg[27] = -1*fr[27]+fl[27]; 
  favg[28] = 1*fr[28]+fl[28]; 
  favg[29] = -1*fr[29]+fl[29]; 
  favg[30] = -1*fr[30]+fl[30]; 
  favg[31] = -1*fr[31]+fl[31]; 

  double fjump[32]; 
  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(1*fr[6]-fl[6]); 
  fjump[7] = amax*(-1*fr[7]-fl[7]); 
  fjump[8] = amax*(-1*fr[8]-fl[8]); 
  fjump[9] = amax*(1*fr[9]-fl[9]); 
  fjump[10] = amax*(1*fr[10]-fl[10]); 
  fjump[11] = amax*(-1*fr[11]-fl[11]); 
  fjump[12] = amax*(1*fr[12]-fl[12]); 
  fjump[13] = amax*(1*fr[13]-fl[13]); 
  fjump[14] = amax*(-1*fr[14]-fl[14]); 
  fjump[15] = amax*(1*fr[15]-fl[15]); 
  fjump[16] = amax*(-1*fr[16]-fl[16]); 
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(-1*fr[18]-fl[18]); 
  fjump[19] = amax*(-1*fr[19]-fl[19]); 
  fjump[20] = amax*(1*fr[20]-fl[20]); 
  fjump[21] = amax*(-1*fr[21]-fl[21]); 
  fjump[22] = amax*(-1*fr[22]-fl[22]); 
  fjump[23] = amax*(1*fr[23]-fl[23]); 
  fjump[24] = amax*(1*fr[24]-fl[24]); 
  fjump[25] = amax*(-1*fr[25]-fl[25]); 
  fjump[26] = amax*(-1*fr[26]-fl[26]); 
  fjump[27] = amax*(-1*fr[27]-fl[27]); 
  fjump[28] = amax*(1*fr[28]-fl[28]); 
  fjump[29] = amax*(-1*fr[29]-fl[29]); 
  fjump[30] = amax*(-1*fr[30]-fl[30]); 
  fjump[31] = amax*(-1*fr[31]-fl[31]); 

  alpha[0] = 2.828427124746191*(B2[0]*wv2+E0[0])-2.828427124746191*B1[0]*wv3; 
  alpha[1] = 2.828427124746191*(B2[1]*wv2+E0[1])-2.828427124746191*B1[1]*wv3; 
  alpha[2] = 2.828427124746191*(B2[2]*wv2+E0[2])-2.828427124746191*B1[2]*wv3; 
  alpha[4] = 0.8164965809277261*B2[0]*dv2; 
  alpha[5] = -0.8164965809277261*B1[0]*dv3; 
  alpha[6] = 2.828427124746191*(B2[3]*wv2+E0[3])-2.828427124746191*B1[3]*wv3; 
  alpha[9] = 0.8164965809277261*B2[1]*dv2; 
  alpha[10] = 0.8164965809277261*B2[2]*dv2; 
  alpha[12] = -0.8164965809277261*B1[1]*dv3; 
  alpha[13] = -0.8164965809277261*B1[2]*dv3; 
  alpha[17] = 0.8164965809277261*B2[3]*dv2; 
  alpha[20] = -0.8164965809277261*B1[3]*dv3; 
  const double amid = 0.1767766952966368*alpha[0]; 

  Ghat[0] = alpha[20]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[17]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[13]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[12]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+alpha[10]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])+alpha[9]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])+alpha[6]*(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])+alpha[5]*(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])+alpha[4]*(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])+alpha[2]*(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])+alpha[1]*(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])-0.8660254037844386*fjump[3]+alpha[0]*(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])-0.5*fjump[0]; 
  Ghat[1] = alpha[13]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[10]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[20]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[5]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+alpha[17]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])+alpha[4]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])+alpha[2]*(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])+alpha[12]*(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])+alpha[9]*(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])+alpha[6]*(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])-0.8660254037844386*fjump[7]+alpha[0]*(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])+alpha[1]*(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])-0.5*fjump[1]; 
  Ghat[2] = alpha[12]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[9]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[5]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[20]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+alpha[4]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])+alpha[17]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])+alpha[1]*(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])+alpha[13]*(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])+alpha[10]*(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])-0.8660254037844386*fjump[8]+alpha[0]*(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])+alpha[6]*(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])+alpha[2]*(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])-0.5*fjump[2]; 
  Ghat[4] = alpha[20]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[13]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[12]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])+alpha[6]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[5]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])+alpha[2]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])+alpha[1]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])+(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])*alpha[17]-0.8660254037844386*fjump[11]+alpha[0]*(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])+(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])*alpha[10]+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[9]-0.5*fjump[4]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[4]; 
  Ghat[5] = alpha[17]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[10]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[9]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])+alpha[6]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[4]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])+alpha[2]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[1]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])*alpha[20]-0.8660254037844386*fjump[14]+alpha[0]*(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])+(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])*alpha[13]+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[12]-0.5*fjump[5]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[5]; 
  Ghat[6] = alpha[5]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[4]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[12]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[13]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])*alpha[20]+alpha[9]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])+alpha[10]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])+(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])*alpha[17]-0.8660254037844386*fjump[16]+alpha[0]*(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])+alpha[1]*(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])+alpha[2]*(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])-0.5*fjump[6]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[6]; 
  Ghat[9] = alpha[13]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[20]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[5]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])+alpha[2]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[12]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])+alpha[6]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])-0.8660254037844386*fjump[18]+alpha[0]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])+(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])*alpha[17]+alpha[10]*(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])+alpha[1]*(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])-0.5*fjump[9]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[9]+alpha[4]*(0.1530931089239486*favg[7]+0.0883883476483184*favg[1]); 
  Ghat[10] = alpha[12]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[5]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[20]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])+alpha[1]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[13]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])-0.8660254037844386*fjump[19]+alpha[0]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])+alpha[6]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[17]+alpha[9]*(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])+alpha[2]*(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])-0.5*fjump[10]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[10]+alpha[4]*(0.1530931089239486*favg[8]+0.0883883476483184*favg[2]); 
  Ghat[12] = alpha[10]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[17]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[4]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])+alpha[2]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[9]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])+alpha[6]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])-0.8660254037844386*fjump[21]+alpha[0]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])*alpha[20]+alpha[13]*(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])+alpha[1]*(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])-0.5*fjump[12]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[12]+alpha[5]*(0.1530931089239486*favg[7]+0.0883883476483184*favg[1]); 
  Ghat[13] = alpha[9]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[4]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[17]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])+alpha[1]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[10]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])-0.8660254037844386*fjump[22]+alpha[0]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[6]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[20]+alpha[12]*(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])+alpha[2]*(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])-0.5*fjump[13]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[13]+alpha[5]*(0.1530931089239486*favg[8]+0.0883883476483184*favg[2]); 
  Ghat[15] = alpha[6]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[2]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[1]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])+alpha[17]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[20]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])-0.8660254037844386*fjump[25]+alpha[0]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])+alpha[10]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[9]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+alpha[13]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])+alpha[12]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])-0.5*fjump[15]+alpha[4]*(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])+alpha[5]*(0.1530931089239486*favg[11]+0.0883883476483184*favg[4]); 
  Ghat[17] = alpha[5]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[12]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[13]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])-0.8660254037844386*fjump[26]+alpha[0]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[20]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])+alpha[1]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])+alpha[2]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])-0.5*fjump[17]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[17]+alpha[4]*(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])+alpha[6]*(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[10]+(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])*alpha[9]; 
  Ghat[20] = alpha[4]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[9]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[10]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])-0.8660254037844386*fjump[27]+alpha[0]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[17]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])+alpha[1]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[2]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])-0.5*fjump[20]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[20]+alpha[5]*(0.1530931089239486*favg[16]+0.0883883476483184*favg[6])+alpha[6]*(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[13]+(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])*alpha[12]; 
  Ghat[23] = alpha[2]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[6]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])-0.8660254037844386*fjump[29]+alpha[0]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])+alpha[10]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[13]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[1]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])-0.5*fjump[23]+alpha[17]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[4]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])*alpha[20]+alpha[5]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])+alpha[9]*(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])+(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])*alpha[12]; 
  Ghat[24] = alpha[1]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])-0.8660254037844386*fjump[30]+alpha[0]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[6]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])+alpha[9]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[12]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[2]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])-0.5*fjump[24]+alpha[4]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[17]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])*alpha[20]+alpha[5]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])+alpha[10]*(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])+(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])*alpha[13]; 
  Ghat[28] = (-0.8660254037844386*fjump[31])+alpha[0]*(0.1530931089239486*favg[31]+0.0883883476483184*favg[28])+alpha[1]*(0.1530931089239486*favg[30]+0.0883883476483184*favg[24])+alpha[2]*(0.1530931089239486*favg[29]+0.0883883476483184*favg[23])-0.5*fjump[28]+alpha[4]*(0.1530931089239486*favg[27]+0.0883883476483184*favg[20])+alpha[5]*(0.1530931089239486*favg[26]+0.0883883476483184*favg[17])+alpha[6]*(0.1530931089239486*favg[25]+0.0883883476483184*favg[15])+alpha[9]*(0.1530931089239486*favg[22]+0.0883883476483184*favg[13])+alpha[10]*(0.1530931089239486*favg[21]+0.0883883476483184*favg[12])+(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])*alpha[20]+alpha[12]*(0.1530931089239486*favg[19]+0.0883883476483184*favg[10])+alpha[13]*(0.1530931089239486*favg[18]+0.0883883476483184*favg[9])+(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])*alpha[17]; 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += 0.5*Ghat[2]*dv10r; 
  outr[3] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[4] += 0.5*Ghat[4]*dv10r; 
  outr[5] += 0.5*Ghat[5]*dv10r; 
  outr[6] += 0.5*Ghat[6]*dv10r; 
  outr[7] += -0.8660254037844386*Ghat[1]*dv10r; 
  outr[8] += -0.8660254037844386*Ghat[2]*dv10r; 
  outr[9] += 0.5*Ghat[9]*dv10r; 
  outr[10] += 0.5*Ghat[10]*dv10r; 
  outr[11] += -0.8660254037844386*Ghat[4]*dv10r; 
  outr[12] += 0.5*Ghat[12]*dv10r; 
  outr[13] += 0.5*Ghat[13]*dv10r; 
  outr[14] += -0.8660254037844386*Ghat[5]*dv10r; 
  outr[15] += 0.5*Ghat[15]*dv10r; 
  outr[16] += -0.8660254037844386*Ghat[6]*dv10r; 
  outr[17] += 0.5*Ghat[17]*dv10r; 
  outr[18] += -0.8660254037844386*Ghat[9]*dv10r; 
  outr[19] += -0.8660254037844386*Ghat[10]*dv10r; 
  outr[20] += 0.5*Ghat[20]*dv10r; 
  outr[21] += -0.8660254037844386*Ghat[12]*dv10r; 
  outr[22] += -0.8660254037844386*Ghat[13]*dv10r; 
  outr[23] += 0.5*Ghat[23]*dv10r; 
  outr[24] += 0.5*Ghat[24]*dv10r; 
  outr[25] += -0.8660254037844386*Ghat[15]*dv10r; 
  outr[26] += -0.8660254037844386*Ghat[17]*dv10r; 
  outr[27] += -0.8660254037844386*Ghat[20]*dv10r; 
  outr[28] += 0.5*Ghat[28]*dv10r; 
  outr[29] += -0.8660254037844386*Ghat[23]*dv10r; 
  outr[30] += -0.8660254037844386*Ghat[24]*dv10r; 
  outr[31] += -0.8660254037844386*Ghat[28]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.5*Ghat[2]*dv10l; 
  outl[3] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[4] += -0.5*Ghat[4]*dv10l; 
  outl[5] += -0.5*Ghat[5]*dv10l; 
  outl[6] += -0.5*Ghat[6]*dv10l; 
  outl[7] += -0.8660254037844386*Ghat[1]*dv10l; 
  outl[8] += -0.8660254037844386*Ghat[2]*dv10l; 
  outl[9] += -0.5*Ghat[9]*dv10l; 
  outl[10] += -0.5*Ghat[10]*dv10l; 
  outl[11] += -0.8660254037844386*Ghat[4]*dv10l; 
  outl[12] += -0.5*Ghat[12]*dv10l; 
  outl[13] += -0.5*Ghat[13]*dv10l; 
  outl[14] += -0.8660254037844386*Ghat[5]*dv10l; 
  outl[15] += -0.5*Ghat[15]*dv10l; 
  outl[16] += -0.8660254037844386*Ghat[6]*dv10l; 
  outl[17] += -0.5*Ghat[17]*dv10l; 
  outl[18] += -0.8660254037844386*Ghat[9]*dv10l; 
  outl[19] += -0.8660254037844386*Ghat[10]*dv10l; 
  outl[20] += -0.5*Ghat[20]*dv10l; 
  outl[21] += -0.8660254037844386*Ghat[12]*dv10l; 
  outl[22] += -0.8660254037844386*Ghat[13]*dv10l; 
  outl[23] += -0.5*Ghat[23]*dv10l; 
  outl[24] += -0.5*Ghat[24]*dv10l; 
  outl[25] += -0.8660254037844386*Ghat[15]*dv10l; 
  outl[26] += -0.8660254037844386*Ghat[17]*dv10l; 
  outl[27] += -0.8660254037844386*Ghat[20]*dv10l; 
  outl[28] += -0.5*Ghat[28]*dv10l; 
  outl[29] += -0.8660254037844386*Ghat[23]*dv10l; 
  outl[30] += -0.8660254037844386*Ghat[24]*dv10l; 
  outl[31] += -0.8660254037844386*Ghat[28]*dv10l; 

  return std::abs(amid); 
} 

