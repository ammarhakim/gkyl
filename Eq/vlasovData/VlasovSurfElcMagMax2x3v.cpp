#include <VlasovModDecl.h> 
double VlasovSurfElcMag2x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[6]; 

  for(unsigned int i=0; i<6; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[6]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  double fjump[6]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  double alpha[6]; 

  alpha[0] = (-2.828427124746191*B1[0]*wv3)+2.828427124746191*B2[0]*wv2+2.828427124746191*E0[0]; 
  alpha[1] = (-2.828427124746191*B1[1]*wv3)+2.828427124746191*B2[1]*wv2+2.828427124746191*E0[1]; 
  alpha[2] = (-2.828427124746191*B1[2]*wv3)+2.828427124746191*B2[2]*wv2+2.828427124746191*E0[2]; 
  alpha[4] = 0.8164965809277261*B2[0]*dv2; 
  alpha[5] = -0.8164965809277261*B1[0]*dv3; 
  const double amid = 0.1767766952966368*alpha[0]; 
  Ghat[0] += 0.0883883476483184*(alpha[5]*favg[5]+alpha[4]*favg[4])-0.8660254037844386*fjump[3]+0.1530931089239486*alpha[0]*favg[3]+0.0883883476483184*(alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.0883883476483184*alpha[0]*favg[0]; 
  Ghat[1] += 0.1530931089239486*alpha[1]*favg[3]-0.5*fjump[1]+0.0883883476483184*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += 0.1530931089239486*alpha[2]*favg[3]-0.5*fjump[2]+0.0883883476483184*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[4] += (-0.5*fjump[4])+0.0883883476483184*alpha[0]*favg[4]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[4]; 
  Ghat[5] += (-0.5*fjump[5])+0.0883883476483184*alpha[0]*favg[5]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[5]; 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += 0.5*Ghat[2]*dv10r; 
  outr[3] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[4] += 0.5*Ghat[4]*dv10r; 
  outr[5] += 0.5*Ghat[5]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.5*Ghat[2]*dv10l; 
  outl[3] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[4] += -0.5*Ghat[4]*dv10l; 
  outl[5] += -0.5*Ghat[5]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double Ghat[21]; 

  for(unsigned int i=0; i<21; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[21]; 

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
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  double fjump[21]; 

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
  fjump[16] = amax*(1*fr[16]-fl[16]); 
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(1*fr[18]-fl[18]); 
  fjump[19] = amax*(1*fr[19]-fl[19]); 
  fjump[20] = amax*(1*fr[20]-fl[20]); 
  double alpha[21]; 

  alpha[0] = (-2.828427124746191*B1[0]*wv3)+2.828427124746191*B2[0]*wv2+2.828427124746191*E0[0]; 
  alpha[1] = (-2.828427124746191*B1[1]*wv3)+2.828427124746191*B2[1]*wv2+2.828427124746191*E0[1]; 
  alpha[2] = (-2.828427124746191*B1[2]*wv3)+2.828427124746191*B2[2]*wv2+2.828427124746191*E0[2]; 
  alpha[4] = 0.8164965809277261*B2[0]*dv2; 
  alpha[5] = -0.8164965809277261*B1[0]*dv3; 
  alpha[6] = (-2.828427124746191*B1[3]*wv3)+2.828427124746191*B2[3]*wv2+2.828427124746191*E0[3]; 
  alpha[9] = 0.8164965809277261*B2[1]*dv2; 
  alpha[10] = 0.8164965809277261*B2[2]*dv2; 
  alpha[12] = -0.8164965809277261*B1[1]*dv3; 
  alpha[13] = -0.8164965809277261*B1[2]*dv3; 
  alpha[16] = (-2.828427124746191*B1[4]*wv3)+2.828427124746191*B2[4]*wv2+2.828427124746191*E0[4]; 
  alpha[17] = (-2.828427124746191*B1[5]*wv3)+2.828427124746191*B2[5]*wv2+2.828427124746191*E0[5]; 
  const double amid = (-0.1976423537605236*alpha[17])-0.1976423537605236*alpha[16]+0.1767766952966368*alpha[0]; 
  Ghat[0] += (-1.118033988749895*fjump[18])+0.1976423537605236*alpha[0]*favg[18]+0.0883883476483184*(alpha[17]*favg[17]+alpha[16]*favg[16])+0.1530931089239486*alpha[5]*favg[14]+0.0883883476483184*(alpha[13]*favg[13]+alpha[12]*favg[12])+0.1530931089239486*alpha[4]*favg[11]+0.0883883476483184*(alpha[10]*favg[10]+alpha[9]*favg[9])+0.1530931089239486*(alpha[2]*favg[8]+alpha[1]*favg[7])+0.0883883476483184*(alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4])-0.8660254037844386*fjump[3]+0.1530931089239486*alpha[0]*favg[3]+0.0883883476483184*(alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.0883883476483184*alpha[0]*favg[0]; 
  Ghat[1] += alpha[1]*(0.1976423537605236*favg[18]+0.07905694150420944*favg[16])+(0.1369306393762915*favg[7]+0.07905694150420944*favg[1])*alpha[16]+0.1530931089239486*alpha[12]*favg[14]+0.0883883476483184*(alpha[5]*favg[12]+favg[5]*alpha[12])+0.1530931089239486*alpha[9]*favg[11]+0.0883883476483184*(alpha[4]*favg[9]+favg[4]*alpha[9])+0.1530931089239486*alpha[6]*favg[8]-0.8660254037844386*fjump[7]+0.1530931089239486*alpha[0]*favg[7]+0.0883883476483184*(alpha[2]*favg[6]+favg[2]*alpha[6])+0.1530931089239486*alpha[1]*favg[3]-0.5*fjump[1]+0.0883883476483184*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += alpha[2]*(0.1976423537605236*favg[18]+0.07905694150420944*favg[17])+(0.1369306393762915*favg[8]+0.07905694150420944*favg[2])*alpha[17]+0.1530931089239486*alpha[13]*favg[14]+0.0883883476483184*(alpha[5]*favg[13]+favg[5]*alpha[13])+0.1530931089239486*alpha[10]*favg[11]+0.0883883476483184*(alpha[4]*favg[10]+favg[4]*alpha[10])-0.8660254037844386*fjump[8]+0.1530931089239486*(alpha[0]*favg[8]+alpha[6]*favg[7])+0.0883883476483184*(alpha[1]*favg[6]+favg[1]*alpha[6])+0.1530931089239486*alpha[2]*favg[3]-0.5*fjump[2]+0.0883883476483184*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[4] += alpha[4]*(0.07905694150420944*favg[19]+0.1976423537605236*favg[18])+0.0883883476483184*alpha[5]*favg[15]-0.8660254037844386*fjump[11]+0.1530931089239486*alpha[0]*favg[11]+0.0883883476483184*alpha[2]*favg[10]+0.1530931089239486*favg[8]*alpha[10]+0.0883883476483184*(favg[2]*alpha[10]+alpha[1]*favg[9])+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[9]-0.5*fjump[4]+0.0883883476483184*alpha[0]*favg[4]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[4]; 
  Ghat[5] += alpha[5]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[18])+0.0883883476483184*alpha[4]*favg[15]-0.8660254037844386*fjump[14]+0.1530931089239486*alpha[0]*favg[14]+0.0883883476483184*alpha[2]*favg[13]+0.1530931089239486*favg[8]*alpha[13]+0.0883883476483184*(favg[2]*alpha[13]+alpha[1]*favg[12])+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[12]-0.5*fjump[5]+0.0883883476483184*alpha[0]*favg[5]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[5]; 
  Ghat[6] += 0.1976423537605236*alpha[6]*favg[18]+0.07905694150420944*(alpha[6]*favg[17]+favg[6]*alpha[17]+alpha[6]*favg[16]+favg[6]*alpha[16])+0.0883883476483184*(alpha[12]*favg[13]+favg[12]*alpha[13]+alpha[9]*favg[10]+favg[9]*alpha[10])+0.1530931089239486*(alpha[1]*favg[8]+alpha[2]*favg[7])-0.5*fjump[6]+0.0883883476483184*alpha[0]*favg[6]+0.1530931089239486*favg[3]*alpha[6]+0.0883883476483184*(favg[0]*alpha[6]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[9] += alpha[9]*(0.07905694150420944*favg[19]+0.1976423537605236*favg[18])+0.07905694150420944*(alpha[9]*favg[16]+favg[9]*alpha[16])+0.0883883476483184*alpha[12]*favg[15]+0.1530931089239486*alpha[1]*favg[11]+0.0883883476483184*(alpha[6]*favg[10]+favg[6]*alpha[10])-0.5*fjump[9]+0.0883883476483184*alpha[0]*favg[9]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[9]+0.1530931089239486*alpha[4]*favg[7]+0.0883883476483184*(alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[10] += alpha[10]*(0.07905694150420944*favg[19]+0.1976423537605236*favg[18])+0.07905694150420944*(alpha[10]*favg[17]+favg[10]*alpha[17])+0.0883883476483184*alpha[13]*favg[15]+0.1530931089239486*alpha[2]*favg[11]-0.5*fjump[10]+0.0883883476483184*alpha[0]*favg[10]+0.1530931089239486*favg[3]*alpha[10]+0.0883883476483184*(favg[0]*alpha[10]+alpha[6]*favg[9]+favg[6]*alpha[9])+0.1530931089239486*alpha[4]*favg[8]+0.0883883476483184*(alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[12] += alpha[12]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[18])+0.07905694150420944*(alpha[12]*favg[16]+favg[12]*alpha[16])+0.0883883476483184*alpha[9]*favg[15]+0.1530931089239486*alpha[1]*favg[14]+0.0883883476483184*(alpha[6]*favg[13]+favg[6]*alpha[13])-0.5*fjump[12]+0.0883883476483184*alpha[0]*favg[12]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[12]+0.1530931089239486*alpha[5]*favg[7]+0.0883883476483184*(alpha[1]*favg[5]+favg[1]*alpha[5]); 
  Ghat[13] += alpha[13]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[18])+0.07905694150420944*(alpha[13]*favg[17]+favg[13]*alpha[17])+0.0883883476483184*alpha[10]*favg[15]+0.1530931089239486*alpha[2]*favg[14]-0.5*fjump[13]+0.0883883476483184*alpha[0]*favg[13]+0.1530931089239486*favg[3]*alpha[13]+0.0883883476483184*(favg[0]*alpha[13]+alpha[6]*favg[12]+favg[6]*alpha[12])+0.1530931089239486*alpha[5]*favg[8]+0.0883883476483184*(alpha[2]*favg[5]+favg[2]*alpha[5]); 
  Ghat[15] += (-0.5*fjump[15])+0.0883883476483184*alpha[0]*favg[15]+0.1530931089239486*alpha[4]*favg[14]+0.0883883476483184*(alpha[10]*favg[13]+favg[10]*alpha[13]+alpha[9]*favg[12]+favg[9]*alpha[12])+0.1530931089239486*alpha[5]*favg[11]+0.0883883476483184*(alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[16] += 0.1976423537605236*alpha[16]*favg[18]-0.5*fjump[16]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[16]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[16]+0.07905694150420944*(alpha[12]*favg[12]+alpha[9]*favg[9])+0.1369306393762915*alpha[1]*favg[7]+0.07905694150420944*(alpha[6]*favg[6]+alpha[1]*favg[1]); 
  Ghat[17] += 0.1976423537605236*alpha[17]*favg[18]-0.5*fjump[17]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[17]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[17]+0.07905694150420944*(alpha[13]*favg[13]+alpha[10]*favg[10])+0.1369306393762915*alpha[2]*favg[8]+0.07905694150420944*(alpha[6]*favg[6]+alpha[2]*favg[2]); 
  Ghat[19] += (-0.5*fjump[19])+0.0883883476483184*alpha[0]*favg[19]+0.1369306393762915*alpha[4]*favg[11]+0.07905694150420944*(alpha[10]*favg[10]+alpha[9]*favg[9]+alpha[4]*favg[4]); 
  Ghat[20] += (-0.5*fjump[20])+0.0883883476483184*alpha[0]*favg[20]+0.1369306393762915*alpha[5]*favg[14]+0.07905694150420944*(alpha[13]*favg[13]+alpha[12]*favg[12]+alpha[5]*favg[5]); 

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
  outr[16] += 0.5*Ghat[16]*dv10r; 
  outr[17] += 0.5*Ghat[17]*dv10r; 
  outr[18] += 1.118033988749895*Ghat[0]*dv10r; 
  outr[19] += 0.5*Ghat[19]*dv10r; 
  outr[20] += 0.5*Ghat[20]*dv10r; 

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
  outl[16] += -0.5*Ghat[16]*dv10l; 
  outl[17] += -0.5*Ghat[17]*dv10l; 
  outl[18] += -1.118033988749895*Ghat[0]*dv10l; 
  outl[19] += -0.5*Ghat[19]*dv10l; 
  outl[20] += -0.5*Ghat[20]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  const double *B0 = &EM[30]; 
  const double *B1 = &EM[40]; 
  const double *B2 = &EM[50]; 

  double Ghat[56]; 

  for(unsigned int i=0; i<56; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[56]; 

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
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  favg[21] = -1*fr[21]+fl[21]; 
  favg[22] = 1*fr[22]+fl[22]; 
  favg[23] = -1*fr[23]+fl[23]; 
  favg[24] = -1*fr[24]+fl[24]; 
  favg[25] = 1*fr[25]+fl[25]; 
  favg[26] = -1*fr[26]+fl[26]; 
  favg[27] = -1*fr[27]+fl[27]; 
  favg[28] = 1*fr[28]+fl[28]; 
  favg[29] = 1*fr[29]+fl[29]; 
  favg[30] = -1*fr[30]+fl[30]; 
  favg[31] = 1*fr[31]+fl[31]; 
  favg[32] = 1*fr[32]+fl[32]; 
  favg[33] = -1*fr[33]+fl[33]; 
  favg[34] = -1*fr[34]+fl[34]; 
  favg[35] = 1*fr[35]+fl[35]; 
  favg[36] = 1*fr[36]+fl[36]; 
  favg[37] = 1*fr[37]+fl[37]; 
  favg[38] = 1*fr[38]+fl[38]; 
  favg[39] = 1*fr[39]+fl[39]; 
  favg[40] = 1*fr[40]+fl[40]; 
  favg[41] = 1*fr[41]+fl[41]; 
  favg[42] = -1*fr[42]+fl[42]; 
  favg[43] = 1*fr[43]+fl[43]; 
  favg[44] = 1*fr[44]+fl[44]; 
  favg[45] = 1*fr[45]+fl[45]; 
  favg[46] = 1*fr[46]+fl[46]; 
  favg[47] = 1*fr[47]+fl[47]; 
  favg[48] = 1*fr[48]+fl[48]; 
  favg[49] = -1*fr[49]+fl[49]; 
  favg[50] = 1*fr[50]+fl[50]; 
  favg[51] = 1*fr[51]+fl[51]; 
  favg[52] = 1*fr[52]+fl[52]; 
  favg[53] = -1*fr[53]+fl[53]; 
  favg[54] = 1*fr[54]+fl[54]; 
  favg[55] = 1*fr[55]+fl[55]; 
  double fjump[56]; 

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
  fjump[16] = amax*(1*fr[16]-fl[16]); 
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(1*fr[18]-fl[18]); 
  fjump[19] = amax*(1*fr[19]-fl[19]); 
  fjump[20] = amax*(1*fr[20]-fl[20]); 
  fjump[21] = amax*(-1*fr[21]-fl[21]); 
  fjump[22] = amax*(1*fr[22]-fl[22]); 
  fjump[23] = amax*(-1*fr[23]-fl[23]); 
  fjump[24] = amax*(-1*fr[24]-fl[24]); 
  fjump[25] = amax*(1*fr[25]-fl[25]); 
  fjump[26] = amax*(-1*fr[26]-fl[26]); 
  fjump[27] = amax*(-1*fr[27]-fl[27]); 
  fjump[28] = amax*(1*fr[28]-fl[28]); 
  fjump[29] = amax*(1*fr[29]-fl[29]); 
  fjump[30] = amax*(-1*fr[30]-fl[30]); 
  fjump[31] = amax*(1*fr[31]-fl[31]); 
  fjump[32] = amax*(1*fr[32]-fl[32]); 
  fjump[33] = amax*(-1*fr[33]-fl[33]); 
  fjump[34] = amax*(-1*fr[34]-fl[34]); 
  fjump[35] = amax*(1*fr[35]-fl[35]); 
  fjump[36] = amax*(1*fr[36]-fl[36]); 
  fjump[37] = amax*(1*fr[37]-fl[37]); 
  fjump[38] = amax*(1*fr[38]-fl[38]); 
  fjump[39] = amax*(1*fr[39]-fl[39]); 
  fjump[40] = amax*(1*fr[40]-fl[40]); 
  fjump[41] = amax*(1*fr[41]-fl[41]); 
  fjump[42] = amax*(-1*fr[42]-fl[42]); 
  fjump[43] = amax*(1*fr[43]-fl[43]); 
  fjump[44] = amax*(1*fr[44]-fl[44]); 
  fjump[45] = amax*(1*fr[45]-fl[45]); 
  fjump[46] = amax*(1*fr[46]-fl[46]); 
  fjump[47] = amax*(1*fr[47]-fl[47]); 
  fjump[48] = amax*(1*fr[48]-fl[48]); 
  fjump[49] = amax*(-1*fr[49]-fl[49]); 
  fjump[50] = amax*(1*fr[50]-fl[50]); 
  fjump[51] = amax*(1*fr[51]-fl[51]); 
  fjump[52] = amax*(1*fr[52]-fl[52]); 
  fjump[53] = amax*(-1*fr[53]-fl[53]); 
  fjump[54] = amax*(1*fr[54]-fl[54]); 
  fjump[55] = amax*(1*fr[55]-fl[55]); 
  double alpha[56]; 

  alpha[0] = (-2.828427124746191*B1[0]*wv3)+2.828427124746191*B2[0]*wv2+2.828427124746191*E0[0]; 
  alpha[1] = (-2.828427124746191*B1[1]*wv3)+2.828427124746191*B2[1]*wv2+2.828427124746191*E0[1]; 
  alpha[2] = (-2.828427124746191*B1[2]*wv3)+2.828427124746191*B2[2]*wv2+2.828427124746191*E0[2]; 
  alpha[4] = 0.8164965809277261*B2[0]*dv2; 
  alpha[5] = -0.8164965809277261*B1[0]*dv3; 
  alpha[6] = (-2.828427124746191*B1[3]*wv3)+2.828427124746191*B2[3]*wv2+2.828427124746191*E0[3]; 
  alpha[9] = 0.8164965809277261*B2[1]*dv2; 
  alpha[10] = 0.8164965809277261*B2[2]*dv2; 
  alpha[12] = -0.8164965809277261*B1[1]*dv3; 
  alpha[13] = -0.8164965809277261*B1[2]*dv3; 
  alpha[16] = (-2.828427124746191*B1[4]*wv3)+2.828427124746191*B2[4]*wv2+2.828427124746191*E0[4]; 
  alpha[17] = (-2.828427124746191*B1[5]*wv3)+2.828427124746191*B2[5]*wv2+2.828427124746191*E0[5]; 
  alpha[22] = 0.8164965809277261*B2[3]*dv2; 
  alpha[25] = -0.8164965809277261*B1[3]*dv3; 
  alpha[31] = (-2.828427124746191*B1[6]*wv3)+2.828427124746191*B2[6]*wv2+2.828427124746191*E0[6]; 
  alpha[32] = (-2.828427124746191*B1[7]*wv3)+2.828427124746191*B2[7]*wv2+2.828427124746191*E0[7]; 
  alpha[37] = 0.816496580927726*B2[4]*dv2; 
  alpha[38] = 0.816496580927726*B2[5]*dv2; 
  alpha[43] = -0.816496580927726*B1[4]*dv3; 
  alpha[44] = -0.816496580927726*B1[5]*dv3; 
  alpha[51] = (-2.828427124746191*B1[8]*wv3)+2.828427124746191*B2[8]*wv2+2.828427124746191*E0[8]; 
  alpha[52] = (-2.828427124746191*B1[9]*wv3)+2.828427124746191*B2[9]*wv2+2.828427124746191*E0[9]; 
  const double amid = (-0.1976423537605236*alpha[17])-0.1976423537605236*alpha[16]+0.1767766952966368*alpha[0]; 
  Ghat[0] += (-1.322875655532295*fjump[53])+0.2338535866733712*alpha[0]*favg[53]+0.0883883476483184*(alpha[52]*favg[52]+alpha[51]*favg[51])+0.1976423537605236*alpha[5]*favg[45]+0.0883883476483184*(alpha[44]*favg[44]+alpha[43]*favg[43])+0.1976423537605236*alpha[4]*favg[39]+0.0883883476483184*(alpha[38]*favg[38]+alpha[37]*favg[37])+0.1976423537605236*(alpha[2]*favg[36]+alpha[1]*favg[35])+0.1530931089239486*(alpha[17]*favg[34]+alpha[16]*favg[33])+0.0883883476483184*(alpha[32]*favg[32]+alpha[31]*favg[31])+0.1530931089239486*(alpha[13]*favg[27]+alpha[12]*favg[26])+0.0883883476483184*alpha[25]*favg[25]+0.1530931089239486*(alpha[10]*favg[24]+alpha[9]*favg[23])+0.0883883476483184*alpha[22]*favg[22]+0.1530931089239486*alpha[6]*favg[21]-1.118033988749895*fjump[18]+0.1976423537605236*alpha[0]*favg[18]+0.0883883476483184*(alpha[17]*favg[17]+alpha[16]*favg[16])+0.1530931089239486*alpha[5]*favg[14]+0.0883883476483184*(alpha[13]*favg[13]+alpha[12]*favg[12])+0.1530931089239486*alpha[4]*favg[11]+0.0883883476483184*(alpha[10]*favg[10]+alpha[9]*favg[9])+0.1530931089239486*(alpha[2]*favg[8]+alpha[1]*favg[7])+0.0883883476483184*(alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4])-0.8660254037844386*fjump[3]+0.1530931089239486*alpha[0]*favg[3]+0.0883883476483184*(alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.0883883476483184*alpha[0]*favg[0]; 
  Ghat[1] += 0.2338535866733712*alpha[1]*favg[53]+0.0776323754260148*alpha[16]*favg[51]+(0.1344632185501192*favg[33]+0.0776323754260148*favg[16])*alpha[51]+alpha[12]*(0.1976423537605236*favg[45]+0.07905694150420944*favg[43])+(0.1369306393762915*favg[26]+0.07905694150420944*favg[12])*alpha[43]+alpha[9]*(0.1976423537605236*favg[39]+0.07905694150420944*favg[37])+(0.1369306393762915*favg[23]+0.07905694150420944*favg[9])*alpha[37]+0.1976423537605236*alpha[6]*favg[36]-1.118033988749895*fjump[35]+(0.1767766952966368*alpha[16]+0.1976423537605236*alpha[0])*favg[35]+0.1530931089239486*alpha[32]*favg[34]+0.1369306393762915*alpha[1]*favg[33]+0.0883883476483184*(alpha[17]*favg[32]+favg[17]*alpha[32])+0.07905694150420944*alpha[6]*favg[31]+(0.1369306393762915*favg[21]+0.07905694150420944*favg[6])*alpha[31]+0.1530931089239486*(alpha[25]*favg[27]+alpha[5]*favg[26])+0.0883883476483184*(alpha[13]*favg[25]+favg[13]*alpha[25])+0.1530931089239486*(alpha[22]*favg[24]+alpha[4]*favg[23])+0.0883883476483184*(alpha[10]*favg[22]+favg[10]*alpha[22])+0.1530931089239486*alpha[2]*favg[21]+alpha[1]*(0.1976423537605236*favg[18]+0.07905694150420944*favg[16])+(0.1369306393762915*favg[7]+0.07905694150420944*favg[1])*alpha[16]+0.1530931089239486*alpha[12]*favg[14]+0.0883883476483184*(alpha[5]*favg[12]+favg[5]*alpha[12])+0.1530931089239486*alpha[9]*favg[11]+0.0883883476483184*(alpha[4]*favg[9]+favg[4]*alpha[9])+0.1530931089239486*alpha[6]*favg[8]-0.8660254037844386*fjump[7]+0.1530931089239486*alpha[0]*favg[7]+0.0883883476483184*(alpha[2]*favg[6]+favg[2]*alpha[6])+0.1530931089239486*alpha[1]*favg[3]-0.5*fjump[1]+0.0883883476483184*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += 0.2338535866733712*alpha[2]*favg[53]+0.0776323754260148*alpha[17]*favg[52]+(0.1344632185501192*favg[34]+0.0776323754260148*favg[17])*alpha[52]+alpha[13]*(0.1976423537605236*favg[45]+0.07905694150420944*favg[44])+(0.1369306393762915*favg[27]+0.07905694150420944*favg[13])*alpha[44]+alpha[10]*(0.1976423537605236*favg[39]+0.07905694150420944*favg[38])+(0.1369306393762915*favg[24]+0.07905694150420944*favg[10])*alpha[38]-1.118033988749895*fjump[36]+0.1767766952966368*alpha[17]*favg[36]+0.1976423537605236*(alpha[0]*favg[36]+alpha[6]*favg[35])+0.1369306393762915*alpha[2]*favg[34]+0.1530931089239486*alpha[31]*favg[33]+0.07905694150420944*alpha[6]*favg[32]+(0.1369306393762915*favg[21]+0.07905694150420944*favg[6])*alpha[32]+0.0883883476483184*(alpha[16]*favg[31]+favg[16]*alpha[31])+0.1530931089239486*(alpha[5]*favg[27]+alpha[25]*favg[26])+0.0883883476483184*(alpha[12]*favg[25]+favg[12]*alpha[25])+0.1530931089239486*(alpha[4]*favg[24]+alpha[22]*favg[23])+0.0883883476483184*(alpha[9]*favg[22]+favg[9]*alpha[22])+0.1530931089239486*alpha[1]*favg[21]+alpha[2]*(0.1976423537605236*favg[18]+0.07905694150420944*favg[17])+(0.1369306393762915*favg[8]+0.07905694150420944*favg[2])*alpha[17]+0.1530931089239486*alpha[13]*favg[14]+0.0883883476483184*(alpha[5]*favg[13]+favg[5]*alpha[13])+0.1530931089239486*alpha[10]*favg[11]+0.0883883476483184*(alpha[4]*favg[10]+favg[4]*alpha[10])-0.8660254037844386*fjump[8]+0.1530931089239486*(alpha[0]*favg[8]+alpha[6]*favg[7])+0.0883883476483184*(alpha[1]*favg[6]+favg[1]*alpha[6])+0.1530931089239486*alpha[2]*favg[3]-0.5*fjump[2]+0.0883883476483184*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[4] += alpha[4]*(0.2338535866733712*favg[53]+0.1369306393762915*favg[42])+0.07905694150420944*(alpha[10]*favg[41]+alpha[9]*favg[40])-1.118033988749895*fjump[39]+0.1976423537605236*alpha[0]*favg[39]+0.0883883476483184*alpha[17]*favg[38]+0.1530931089239486*favg[34]*alpha[38]+0.0883883476483184*(favg[17]*alpha[38]+alpha[16]*favg[37])+(0.1530931089239486*favg[33]+0.0883883476483184*favg[16])*alpha[37]+0.1976423537605236*(alpha[10]*favg[36]+alpha[9]*favg[35])+0.1530931089239486*alpha[5]*favg[30]+0.0883883476483184*(alpha[13]*favg[29]+alpha[12]*favg[28])+0.1530931089239486*(alpha[2]*favg[24]+alpha[1]*favg[23])+0.0883883476483184*alpha[6]*favg[22]+(0.1530931089239486*favg[21]+0.0883883476483184*favg[6])*alpha[22]+alpha[4]*(0.07905694150420944*favg[19]+0.1976423537605236*favg[18])+0.0883883476483184*alpha[5]*favg[15]-0.8660254037844386*fjump[11]+0.1530931089239486*alpha[0]*favg[11]+0.0883883476483184*alpha[2]*favg[10]+0.1530931089239486*favg[8]*alpha[10]+0.0883883476483184*(favg[2]*alpha[10]+alpha[1]*favg[9])+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[9]-0.5*fjump[4]+0.0883883476483184*alpha[0]*favg[4]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[4]; 
  Ghat[5] += alpha[5]*(0.2338535866733712*favg[53]+0.1369306393762915*favg[49])+0.07905694150420944*(alpha[13]*favg[48]+alpha[12]*favg[47])-1.118033988749895*fjump[45]+0.1976423537605236*alpha[0]*favg[45]+0.0883883476483184*alpha[17]*favg[44]+0.1530931089239486*favg[34]*alpha[44]+0.0883883476483184*(favg[17]*alpha[44]+alpha[16]*favg[43])+(0.1530931089239486*favg[33]+0.0883883476483184*favg[16])*alpha[43]+0.1976423537605236*(alpha[13]*favg[36]+alpha[12]*favg[35])+0.1530931089239486*alpha[4]*favg[30]+0.0883883476483184*(alpha[10]*favg[29]+alpha[9]*favg[28])+0.1530931089239486*(alpha[2]*favg[27]+alpha[1]*favg[26])+0.0883883476483184*alpha[6]*favg[25]+(0.1530931089239486*favg[21]+0.0883883476483184*favg[6])*alpha[25]+alpha[5]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[18])+0.0883883476483184*alpha[4]*favg[15]-0.8660254037844386*fjump[14]+0.1530931089239486*alpha[0]*favg[14]+0.0883883476483184*alpha[2]*favg[13]+0.1530931089239486*favg[8]*alpha[13]+0.0883883476483184*(favg[2]*alpha[13]+alpha[1]*favg[12])+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[12]-0.5*fjump[5]+0.0883883476483184*alpha[0]*favg[5]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[5]; 
  Ghat[6] += 0.2338535866733712*alpha[6]*favg[53]+0.0776323754260148*(alpha[32]*favg[52]+favg[32]*alpha[52]+alpha[31]*favg[51]+favg[31]*alpha[51])+0.1976423537605236*alpha[25]*favg[45]+0.07905694150420944*(alpha[25]*favg[44]+favg[25]*alpha[44]+alpha[25]*favg[43]+favg[25]*alpha[43])+0.1976423537605236*alpha[22]*favg[39]+0.07905694150420944*(alpha[22]*favg[38]+favg[22]*alpha[38]+alpha[22]*favg[37]+favg[22]*alpha[37])+(0.1767766952966368*alpha[32]+0.1976423537605236*alpha[1])*favg[36]+(0.1767766952966368*alpha[31]+0.1976423537605236*alpha[2])*favg[35]+0.1369306393762915*alpha[6]*(favg[34]+favg[33])+(0.07071067811865474*alpha[31]+0.07905694150420944*alpha[2])*favg[32]+(0.07071067811865474*favg[31]+0.1369306393762915*favg[8])*alpha[32]+0.07905694150420944*(favg[2]*alpha[32]+alpha[1]*favg[31])+(0.1369306393762915*favg[7]+0.07905694150420944*favg[1])*alpha[31]+0.1530931089239486*(alpha[12]*favg[27]+alpha[13]*favg[26])+0.0883883476483184*alpha[5]*favg[25]+(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])*alpha[25]+0.1530931089239486*(alpha[9]*favg[24]+alpha[10]*favg[23])+0.0883883476483184*alpha[4]*favg[22]+(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])*alpha[22]-0.8660254037844386*fjump[21]+(0.1369306393762915*(alpha[17]+alpha[16])+0.1530931089239486*alpha[0])*favg[21]+0.1976423537605236*alpha[6]*favg[18]+0.07905694150420944*(alpha[6]*favg[17]+favg[6]*alpha[17]+alpha[6]*favg[16]+favg[6]*alpha[16])+0.0883883476483184*(alpha[12]*favg[13]+favg[12]*alpha[13]+alpha[9]*favg[10]+favg[9]*alpha[10])+0.1530931089239486*(alpha[1]*favg[8]+alpha[2]*favg[7])-0.5*fjump[6]+0.0883883476483184*alpha[0]*favg[6]+0.1530931089239486*favg[3]*alpha[6]+0.0883883476483184*(favg[0]*alpha[6]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[9] += 0.2338535866733712*alpha[9]*favg[53]+0.0776323754260148*(alpha[37]*favg[51]+favg[37]*alpha[51])+0.07905694150420944*favg[28]*alpha[43]+0.1369306393762915*alpha[9]*favg[42]+0.07905694150420944*alpha[22]*favg[41]+(0.07071067811865474*alpha[37]+0.07905694150420944*alpha[4])*favg[40]+0.1976423537605236*alpha[1]*favg[39]+0.0883883476483184*(alpha[32]*favg[38]+favg[32]*alpha[38])+0.07905694150420944*alpha[1]*favg[37]+(0.1767766952966368*favg[35]+0.1369306393762915*favg[7]+0.07905694150420944*favg[1])*alpha[37]+0.1976423537605236*(alpha[22]*favg[36]+alpha[4]*favg[35])+0.1369306393762915*alpha[9]*favg[33]+0.07905694150420944*(alpha[22]*favg[31]+favg[22]*alpha[31])+0.1530931089239486*alpha[12]*favg[30]+0.0883883476483184*(alpha[25]*favg[29]+alpha[5]*favg[28])+0.1530931089239486*alpha[6]*favg[24]-0.8660254037844386*fjump[23]+(0.1369306393762915*alpha[16]+0.1530931089239486*alpha[0])*favg[23]+0.0883883476483184*alpha[2]*favg[22]+(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])*alpha[22]+0.1530931089239486*alpha[10]*favg[21]+alpha[9]*(0.07905694150420944*favg[19]+0.1976423537605236*favg[18])+0.07905694150420944*(alpha[9]*favg[16]+favg[9]*alpha[16])+0.0883883476483184*alpha[12]*favg[15]+0.1530931089239486*alpha[1]*favg[11]+0.0883883476483184*(alpha[6]*favg[10]+favg[6]*alpha[10])-0.5*fjump[9]+0.0883883476483184*alpha[0]*favg[9]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[9]+0.1530931089239486*alpha[4]*favg[7]+0.0883883476483184*(alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[10] += 0.2338535866733712*alpha[10]*favg[53]+0.0776323754260148*(alpha[38]*favg[52]+favg[38]*alpha[52])+0.07905694150420944*favg[29]*alpha[44]+0.1369306393762915*alpha[10]*favg[42]+0.07071067811865474*alpha[38]*favg[41]+0.07905694150420944*(alpha[4]*favg[41]+alpha[22]*favg[40])+alpha[2]*(0.1976423537605236*favg[39]+0.07905694150420944*favg[38])+(0.1767766952966368*favg[36]+0.1369306393762915*favg[8]+0.07905694150420944*favg[2])*alpha[38]+0.0883883476483184*(alpha[31]*favg[37]+favg[31]*alpha[37])+0.1976423537605236*(alpha[4]*favg[36]+alpha[22]*favg[35])+0.1369306393762915*alpha[10]*favg[34]+0.07905694150420944*(alpha[22]*favg[32]+favg[22]*alpha[32])+0.1530931089239486*alpha[13]*favg[30]+0.0883883476483184*(alpha[5]*favg[29]+alpha[25]*favg[28])-0.8660254037844386*fjump[24]+0.1369306393762915*alpha[17]*favg[24]+0.1530931089239486*(alpha[0]*favg[24]+alpha[6]*favg[23])+0.0883883476483184*alpha[1]*favg[22]+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[22]+0.1530931089239486*alpha[9]*favg[21]+alpha[10]*(0.07905694150420944*favg[19]+0.1976423537605236*favg[18])+0.07905694150420944*(alpha[10]*favg[17]+favg[10]*alpha[17])+0.0883883476483184*alpha[13]*favg[15]+0.1530931089239486*alpha[2]*favg[11]-0.5*fjump[10]+0.0883883476483184*alpha[0]*favg[10]+0.1530931089239486*favg[3]*alpha[10]+0.0883883476483184*(favg[0]*alpha[10]+alpha[6]*favg[9]+favg[6]*alpha[9])+0.1530931089239486*alpha[4]*favg[8]+0.0883883476483184*(alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[12] += 0.2338535866733712*alpha[12]*favg[53]+0.0776323754260148*(alpha[43]*favg[51]+favg[43]*alpha[51])+0.1369306393762915*alpha[12]*favg[49]+0.07905694150420944*alpha[25]*favg[48]+(0.07071067811865474*alpha[43]+0.07905694150420944*alpha[5])*favg[47]+0.1976423537605236*alpha[1]*favg[45]+0.0883883476483184*(alpha[32]*favg[44]+favg[32]*alpha[44])+0.07905694150420944*alpha[1]*favg[43]+(0.1767766952966368*favg[35]+0.1369306393762915*favg[7])*alpha[43]+0.07905694150420944*(favg[1]*alpha[43]+favg[28]*alpha[37])+0.1976423537605236*(alpha[25]*favg[36]+alpha[5]*favg[35])+0.1369306393762915*alpha[12]*favg[33]+0.07905694150420944*(alpha[25]*favg[31]+favg[25]*alpha[31])+0.1530931089239486*alpha[9]*favg[30]+0.0883883476483184*(alpha[22]*favg[29]+alpha[4]*favg[28])+0.1530931089239486*alpha[6]*favg[27]-0.8660254037844386*fjump[26]+(0.1369306393762915*alpha[16]+0.1530931089239486*alpha[0])*favg[26]+0.0883883476483184*alpha[2]*favg[25]+(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])*alpha[25]+0.1530931089239486*alpha[13]*favg[21]+alpha[12]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[18])+0.07905694150420944*(alpha[12]*favg[16]+favg[12]*alpha[16])+0.0883883476483184*alpha[9]*favg[15]+0.1530931089239486*alpha[1]*favg[14]+0.0883883476483184*(alpha[6]*favg[13]+favg[6]*alpha[13])-0.5*fjump[12]+0.0883883476483184*alpha[0]*favg[12]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[12]+0.1530931089239486*alpha[5]*favg[7]+0.0883883476483184*(alpha[1]*favg[5]+favg[1]*alpha[5]); 
  Ghat[13] += 0.2338535866733712*alpha[13]*favg[53]+0.0776323754260148*(alpha[44]*favg[52]+favg[44]*alpha[52])+0.1369306393762915*alpha[13]*favg[49]+0.07071067811865474*alpha[44]*favg[48]+0.07905694150420944*(alpha[5]*favg[48]+alpha[25]*favg[47])+alpha[2]*(0.1976423537605236*favg[45]+0.07905694150420944*favg[44])+(0.1767766952966368*favg[36]+0.1369306393762915*favg[8]+0.07905694150420944*favg[2])*alpha[44]+0.0883883476483184*(alpha[31]*favg[43]+favg[31]*alpha[43])+0.07905694150420944*favg[29]*alpha[38]+0.1976423537605236*(alpha[5]*favg[36]+alpha[25]*favg[35])+0.1369306393762915*alpha[13]*favg[34]+0.07905694150420944*(alpha[25]*favg[32]+favg[25]*alpha[32])+0.1530931089239486*alpha[10]*favg[30]+0.0883883476483184*(alpha[4]*favg[29]+alpha[22]*favg[28])-0.8660254037844386*fjump[27]+0.1369306393762915*alpha[17]*favg[27]+0.1530931089239486*(alpha[0]*favg[27]+alpha[6]*favg[26])+0.0883883476483184*alpha[1]*favg[25]+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[25]+0.1530931089239486*alpha[12]*favg[21]+alpha[13]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[18])+0.07905694150420944*(alpha[13]*favg[17]+favg[13]*alpha[17])+0.0883883476483184*alpha[10]*favg[15]+0.1530931089239486*alpha[2]*favg[14]-0.5*fjump[13]+0.0883883476483184*alpha[0]*favg[13]+0.1530931089239486*favg[3]*alpha[13]+0.0883883476483184*(favg[0]*alpha[13]+alpha[6]*favg[12]+favg[6]*alpha[12])+0.1530931089239486*alpha[5]*favg[8]+0.0883883476483184*(alpha[2]*favg[5]+favg[2]*alpha[5]); 
  Ghat[15] += 0.07905694150420944*alpha[5]*favg[50]+alpha[4]*(0.07905694150420944*favg[46]+0.1976423537605236*favg[45])+0.0883883476483184*(alpha[38]*favg[44]+favg[38]*alpha[44]+alpha[37]*favg[43]+favg[37]*alpha[43])+0.1976423537605236*alpha[5]*favg[39]-0.8660254037844386*fjump[30]+0.1530931089239486*alpha[0]*favg[30]+0.0883883476483184*(alpha[2]*favg[29]+alpha[1]*favg[28])+0.1530931089239486*(alpha[10]*favg[27]+alpha[9]*favg[26])+0.0883883476483184*(alpha[22]*favg[25]+favg[22]*alpha[25])+0.1530931089239486*(alpha[13]*favg[24]+alpha[12]*favg[23])-0.5*fjump[15]+0.0883883476483184*alpha[0]*favg[15]+0.1530931089239486*alpha[4]*favg[14]+0.0883883476483184*(alpha[10]*favg[13]+favg[10]*alpha[13]+alpha[9]*favg[12]+favg[9]*alpha[12])+0.1530931089239486*alpha[5]*favg[11]+0.0883883476483184*(alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[16] += 0.2338535866733712*alpha[16]*favg[53]+(0.05270462766947297*alpha[51]+0.0776323754260148*alpha[1])*favg[51]+(0.1735912687073533*favg[35]+0.1344632185501192*favg[7]+0.0776323754260148*favg[1])*alpha[51]+0.1976423537605236*alpha[43]*favg[45]+(0.05646924393157818*alpha[43]+0.0883883476483184*alpha[5])*favg[43]+(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])*alpha[43]+0.1976423537605236*alpha[37]*favg[39]+(0.05646924393157818*alpha[37]+0.0883883476483184*alpha[4])*favg[37]+(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])*alpha[37]+0.1976423537605236*alpha[31]*favg[36]+0.1767766952966368*alpha[1]*favg[35]-0.8660254037844386*fjump[33]+(0.09780759955449389*alpha[16]+0.1530931089239486*alpha[0])*favg[33]+0.07905694150420944*alpha[32]*favg[32]+(0.05646924393157818*alpha[31]+0.0883883476483184*alpha[2])*favg[31]+(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])*alpha[31]+0.1369306393762915*alpha[12]*favg[26]+0.07905694150420944*alpha[25]*favg[25]+0.1369306393762915*alpha[9]*favg[23]+0.07905694150420944*alpha[22]*favg[22]+0.1369306393762915*alpha[6]*favg[21]+0.1976423537605236*alpha[16]*favg[18]-0.5*fjump[16]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[16]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[16]+0.07905694150420944*(alpha[12]*favg[12]+alpha[9]*favg[9])+0.1369306393762915*alpha[1]*favg[7]+0.07905694150420944*(alpha[6]*favg[6]+alpha[1]*favg[1]); 
  Ghat[17] += 0.2338535866733712*alpha[17]*favg[53]+(0.05270462766947297*alpha[52]+0.0776323754260148*alpha[2])*favg[52]+(0.1735912687073533*favg[36]+0.1344632185501192*favg[8]+0.0776323754260148*favg[2])*alpha[52]+0.1976423537605236*alpha[44]*favg[45]+(0.05646924393157818*alpha[44]+0.0883883476483184*alpha[5])*favg[44]+(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])*alpha[44]+0.1976423537605236*alpha[38]*favg[39]+(0.05646924393157818*alpha[38]+0.0883883476483184*alpha[4])*favg[38]+(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])*alpha[38]+0.1767766952966368*alpha[2]*favg[36]+0.1976423537605236*alpha[32]*favg[35]-0.8660254037844386*fjump[34]+(0.09780759955449389*alpha[17]+0.1530931089239486*alpha[0])*favg[34]+(0.05646924393157818*alpha[32]+0.0883883476483184*alpha[1])*favg[32]+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[32]+0.07905694150420944*alpha[31]*favg[31]+0.1369306393762915*alpha[13]*favg[27]+0.07905694150420944*alpha[25]*favg[25]+0.1369306393762915*alpha[10]*favg[24]+0.07905694150420944*alpha[22]*favg[22]+0.1369306393762915*alpha[6]*favg[21]+0.1976423537605236*alpha[17]*favg[18]-0.5*fjump[17]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[17]+(0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[17]+0.07905694150420944*(alpha[13]*favg[13]+alpha[10]*favg[10])+0.1369306393762915*alpha[2]*favg[8]+0.07905694150420944*(alpha[6]*favg[6]+alpha[2]*favg[2]); 
  Ghat[19] += 0.0776323754260148*alpha[4]*favg[54]+0.0883883476483184*alpha[5]*favg[46]-0.8660254037844386*fjump[42]+0.1530931089239486*alpha[0]*favg[42]+0.0883883476483184*(alpha[2]*favg[41]+alpha[1]*favg[40])+0.1767766952966368*alpha[4]*favg[39]+0.07905694150420944*(alpha[38]*favg[38]+alpha[37]*favg[37])+0.1369306393762915*(alpha[10]*favg[24]+alpha[9]*favg[23])+0.07905694150420944*alpha[22]*favg[22]-0.5*fjump[19]+0.0883883476483184*alpha[0]*favg[19]+0.1369306393762915*alpha[4]*favg[11]+0.07905694150420944*(alpha[10]*favg[10]+alpha[9]*favg[9]+alpha[4]*favg[4]); 
  Ghat[20] += 0.0776323754260148*alpha[5]*favg[55]+0.0883883476483184*alpha[4]*favg[50]-0.8660254037844386*fjump[49]+0.1530931089239486*alpha[0]*favg[49]+0.0883883476483184*(alpha[2]*favg[48]+alpha[1]*favg[47])+0.1767766952966368*alpha[5]*favg[45]+0.07905694150420944*(alpha[44]*favg[44]+alpha[43]*favg[43])+0.1369306393762915*(alpha[13]*favg[27]+alpha[12]*favg[26])+0.07905694150420944*alpha[25]*favg[25]-0.5*fjump[20]+0.0883883476483184*alpha[0]*favg[20]+0.1369306393762915*alpha[5]*favg[14]+0.07905694150420944*(alpha[13]*favg[13]+alpha[12]*favg[12]+alpha[5]*favg[5]); 
  Ghat[22] += alpha[22]*(0.2338535866733712*favg[53]+0.1369306393762915*favg[42])+0.07905694150420944*(alpha[9]*favg[41]+alpha[10]*favg[40])+alpha[6]*(0.1976423537605236*favg[39]+0.07905694150420944*favg[38])+0.1369306393762915*favg[21]*alpha[38]+0.07905694150420944*(favg[6]*alpha[38]+alpha[6]*favg[37])+(0.1369306393762915*favg[21]+0.07905694150420944*favg[6])*alpha[37]+0.1976423537605236*(alpha[9]*favg[36]+alpha[10]*favg[35])+0.1369306393762915*alpha[22]*(favg[34]+favg[33])+0.07905694150420944*alpha[10]*favg[32]+0.1369306393762915*favg[24]*alpha[32]+0.07905694150420944*(favg[10]*alpha[32]+alpha[9]*favg[31])+(0.1369306393762915*favg[23]+0.07905694150420944*favg[9])*alpha[31]+0.1530931089239486*alpha[25]*favg[30]+0.0883883476483184*(alpha[12]*favg[29]+alpha[13]*favg[28]+favg[15]*alpha[25])+0.1530931089239486*(alpha[1]*favg[24]+alpha[2]*favg[23])-0.5*fjump[22]+(0.07905694150420944*(alpha[17]+alpha[16])+0.0883883476483184*alpha[0])*favg[22]+(0.07905694150420944*favg[19]+0.1976423537605236*favg[18]+0.07905694150420944*(favg[17]+favg[16])+0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[22]+0.1530931089239486*(alpha[4]*favg[21]+alpha[6]*favg[11])+0.0883883476483184*alpha[1]*favg[10]+0.1530931089239486*favg[7]*alpha[10]+0.0883883476483184*(favg[1]*alpha[10]+alpha[2]*favg[9])+0.1530931089239486*favg[8]*alpha[9]+0.0883883476483184*(favg[2]*alpha[9]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[25] += alpha[25]*(0.2338535866733712*favg[53]+0.1369306393762915*favg[49])+0.07905694150420944*(alpha[12]*favg[48]+alpha[13]*favg[47])+alpha[6]*(0.1976423537605236*favg[45]+0.07905694150420944*favg[44])+0.1369306393762915*favg[21]*alpha[44]+0.07905694150420944*(favg[6]*alpha[44]+alpha[6]*favg[43])+(0.1369306393762915*favg[21]+0.07905694150420944*favg[6])*alpha[43]+0.1976423537605236*(alpha[12]*favg[36]+alpha[13]*favg[35])+0.1369306393762915*alpha[25]*(favg[34]+favg[33])+0.07905694150420944*alpha[13]*favg[32]+0.1369306393762915*favg[27]*alpha[32]+0.07905694150420944*(favg[13]*alpha[32]+alpha[12]*favg[31])+(0.1369306393762915*favg[26]+0.07905694150420944*favg[12])*alpha[31]+0.1530931089239486*alpha[22]*favg[30]+0.0883883476483184*(alpha[9]*favg[29]+alpha[10]*favg[28])+0.1530931089239486*(alpha[1]*favg[27]+alpha[2]*favg[26])-0.5*fjump[25]+(0.07905694150420944*(alpha[17]+alpha[16])+0.0883883476483184*alpha[0])*favg[25]+(0.07905694150420944*favg[20]+0.1976423537605236*favg[18]+0.07905694150420944*(favg[17]+favg[16])+0.1530931089239486*favg[3])*alpha[25]+0.0883883476483184*(favg[0]*alpha[25]+favg[15]*alpha[22])+0.1530931089239486*(alpha[5]*favg[21]+alpha[6]*favg[14])+0.0883883476483184*alpha[1]*favg[13]+0.1530931089239486*favg[7]*alpha[13]+0.0883883476483184*(favg[1]*alpha[13]+alpha[2]*favg[12])+0.1530931089239486*favg[8]*alpha[12]+0.0883883476483184*(favg[2]*alpha[12]+alpha[5]*favg[6]+favg[5]*alpha[6]); 
  Ghat[28] += 0.07905694150420944*alpha[12]*favg[50]+alpha[9]*(0.07905694150420944*favg[46]+0.1976423537605236*favg[45]+0.07905694150420944*favg[43])+(0.1369306393762915*favg[23]+0.07905694150420944*favg[9])*alpha[43]+alpha[12]*(0.1976423537605236*favg[39]+0.07905694150420944*favg[37])+(0.1369306393762915*favg[26]+0.07905694150420944*favg[12])*alpha[37]+0.1530931089239486*alpha[1]*favg[30]+0.0883883476483184*alpha[6]*favg[29]-0.5*fjump[28]+(0.07905694150420944*alpha[16]+0.0883883476483184*alpha[0])*favg[28]+0.1530931089239486*(alpha[22]*favg[27]+alpha[4]*favg[26])+0.0883883476483184*alpha[10]*favg[25]+(0.1530931089239486*favg[24]+0.0883883476483184*favg[10])*alpha[25]+0.1530931089239486*alpha[5]*favg[23]+0.0883883476483184*(alpha[13]*favg[22]+favg[13]*alpha[22]+alpha[1]*favg[15])+0.1530931089239486*alpha[9]*favg[14]+0.0883883476483184*alpha[4]*favg[12]+0.1530931089239486*favg[11]*alpha[12]+0.0883883476483184*(favg[4]*alpha[12]+alpha[5]*favg[9]+favg[5]*alpha[9]); 
  Ghat[29] += 0.07905694150420944*alpha[13]*favg[50]+alpha[10]*(0.07905694150420944*favg[46]+0.1976423537605236*favg[45]+0.07905694150420944*favg[44])+(0.1369306393762915*favg[24]+0.07905694150420944*favg[10])*alpha[44]+alpha[13]*(0.1976423537605236*favg[39]+0.07905694150420944*favg[38])+(0.1369306393762915*favg[27]+0.07905694150420944*favg[13])*alpha[38]+0.1530931089239486*alpha[2]*favg[30]-0.5*fjump[29]+0.07905694150420944*alpha[17]*favg[29]+0.0883883476483184*(alpha[0]*favg[29]+alpha[6]*favg[28])+0.1530931089239486*(alpha[4]*favg[27]+alpha[22]*favg[26])+0.0883883476483184*alpha[9]*favg[25]+(0.1530931089239486*favg[23]+0.0883883476483184*favg[9])*alpha[25]+0.1530931089239486*alpha[5]*favg[24]+0.0883883476483184*(alpha[12]*favg[22]+favg[12]*alpha[22]+alpha[2]*favg[15])+0.1530931089239486*alpha[10]*favg[14]+0.0883883476483184*alpha[4]*favg[13]+0.1530931089239486*favg[11]*alpha[13]+0.0883883476483184*(favg[4]*alpha[13]+alpha[5]*favg[10]+favg[5]*alpha[10]); 
  Ghat[31] += 0.2338535866733712*alpha[31]*favg[53]+0.0776323754260148*alpha[6]*favg[51]+(0.1344632185501192*favg[21]+0.0776323754260148*favg[6])*alpha[51]+0.0883883476483184*alpha[13]*favg[43]+0.1530931089239486*favg[27]*alpha[43]+0.0883883476483184*(favg[13]*alpha[43]+alpha[10]*favg[37])+(0.1530931089239486*favg[24]+0.0883883476483184*favg[10])*alpha[37]+0.1976423537605236*alpha[16]*favg[36]+0.1767766952966368*alpha[6]*favg[35]+0.1369306393762915*alpha[31]*favg[34]+(0.09780759955449389*alpha[31]+0.1530931089239486*alpha[2])*favg[33]+0.07071067811865474*alpha[6]*favg[32]+(0.1224744871391589*favg[21]+0.07071067811865474*favg[6])*alpha[32]-0.5*fjump[31]+(0.07905694150420944*alpha[17]+0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[31]+(0.1976423537605236*favg[18]+0.07905694150420944*favg[17]+0.05646924393157818*favg[16]+0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[31]+0.1369306393762915*alpha[25]*favg[26]+0.07905694150420944*(alpha[12]*favg[25]+favg[12]*alpha[25])+0.1369306393762915*alpha[22]*favg[23]+0.07905694150420944*(alpha[9]*favg[22]+favg[9]*alpha[22])+0.1369306393762915*alpha[1]*favg[21]+0.0883883476483184*alpha[2]*favg[16]+(0.1530931089239486*favg[8]+0.0883883476483184*favg[2])*alpha[16]+0.1369306393762915*alpha[6]*favg[7]+0.07905694150420944*(alpha[1]*favg[6]+favg[1]*alpha[6]); 
  Ghat[32] += 0.2338535866733712*alpha[32]*favg[53]+0.0776323754260148*alpha[6]*favg[52]+(0.1344632185501192*favg[21]+0.0776323754260148*favg[6])*alpha[52]+0.0883883476483184*alpha[12]*favg[44]+0.1530931089239486*favg[26]*alpha[44]+0.0883883476483184*(favg[12]*alpha[44]+alpha[9]*favg[38])+(0.1530931089239486*favg[23]+0.0883883476483184*favg[9])*alpha[38]+0.1767766952966368*alpha[6]*favg[36]+0.1976423537605236*alpha[17]*favg[35]+(0.09780759955449389*alpha[32]+0.1530931089239486*alpha[1])*favg[34]+0.1369306393762915*alpha[32]*favg[33]-0.5*fjump[32]+(0.05646924393157818*alpha[17]+0.07905694150420944*alpha[16]+0.0883883476483184*alpha[0])*favg[32]+(0.1976423537605236*favg[18]+0.05646924393157818*favg[17]+0.07905694150420944*favg[16]+0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[32]+0.07071067811865474*alpha[6]*favg[31]+(0.1224744871391589*favg[21]+0.07071067811865474*favg[6])*alpha[31]+0.1369306393762915*alpha[25]*favg[27]+0.07905694150420944*(alpha[13]*favg[25]+favg[13]*alpha[25])+0.1369306393762915*alpha[22]*favg[24]+0.07905694150420944*(alpha[10]*favg[22]+favg[10]*alpha[22])+0.1369306393762915*alpha[2]*favg[21]+0.0883883476483184*alpha[1]*favg[17]+(0.1530931089239486*favg[7]+0.0883883476483184*favg[1])*alpha[17]+0.1369306393762915*alpha[6]*favg[8]+0.07905694150420944*(alpha[2]*favg[6]+favg[2]*alpha[6]); 
  Ghat[37] += 0.2338535866733712*alpha[37]*favg[53]+0.0776323754260148*alpha[9]*favg[51]+(0.1344632185501192*favg[23]+0.0776323754260148*favg[9])*alpha[51]+(0.1530931089239486*favg[30]+0.0883883476483184*favg[15])*alpha[43]+0.1369306393762915*alpha[37]*favg[42]+0.07071067811865474*alpha[9]*favg[40]+0.1976423537605236*alpha[16]*favg[39]-0.5*fjump[37]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[37]+(0.09780759955449389*favg[33]+0.07905694150420944*favg[19]+0.1976423537605236*favg[18]+0.05646924393157818*favg[16]+0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[37]+0.1767766952966368*alpha[9]*favg[35]+0.1530931089239486*alpha[4]*favg[33]+0.0883883476483184*alpha[10]*favg[31]+(0.1530931089239486*favg[24]+0.0883883476483184*favg[10])*alpha[31]+0.07905694150420944*alpha[12]*favg[28]+0.1369306393762915*alpha[1]*favg[23]+0.07905694150420944*alpha[6]*favg[22]+(0.1369306393762915*favg[21]+0.07905694150420944*favg[6])*alpha[22]+0.0883883476483184*alpha[4]*favg[16]+(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])*alpha[16]+0.07905694150420944*alpha[1]*favg[9]+(0.1369306393762915*favg[7]+0.07905694150420944*favg[1])*alpha[9]; 
  Ghat[38] += 0.2338535866733712*alpha[38]*favg[53]+0.0776323754260148*alpha[10]*favg[52]+(0.1344632185501192*favg[24]+0.0776323754260148*favg[10])*alpha[52]+(0.1530931089239486*favg[30]+0.0883883476483184*favg[15])*alpha[44]+0.1369306393762915*alpha[38]*favg[42]+0.07071067811865474*alpha[10]*favg[41]+0.1976423537605236*alpha[17]*favg[39]-0.5*fjump[38]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[38]+(0.09780759955449389*favg[34]+0.07905694150420944*favg[19]+0.1976423537605236*favg[18]+0.05646924393157818*favg[17]+0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[38]+0.1767766952966368*alpha[10]*favg[36]+0.1530931089239486*alpha[4]*favg[34]+0.0883883476483184*alpha[9]*favg[32]+(0.1530931089239486*favg[23]+0.0883883476483184*favg[9])*alpha[32]+0.07905694150420944*alpha[13]*favg[29]+0.1369306393762915*alpha[2]*favg[24]+0.07905694150420944*alpha[6]*favg[22]+(0.1369306393762915*favg[21]+0.07905694150420944*favg[6])*alpha[22]+0.0883883476483184*alpha[4]*favg[17]+(0.1530931089239486*favg[11]+0.0883883476483184*favg[4])*alpha[17]+0.07905694150420944*alpha[2]*favg[10]+(0.1369306393762915*favg[8]+0.07905694150420944*favg[2])*alpha[10]; 
  Ghat[40] += 0.0776323754260148*alpha[9]*favg[54]+0.0883883476483184*alpha[12]*favg[46]+0.1530931089239486*alpha[1]*favg[42]+0.0883883476483184*alpha[6]*favg[41]-0.5*fjump[40]+(0.07905694150420944*alpha[16]+0.0883883476483184*alpha[0])*favg[40]+alpha[9]*(0.1767766952966368*favg[39]+0.07071067811865474*favg[37])+(0.1224744871391589*favg[23]+0.07071067811865474*favg[9])*alpha[37]+0.1369306393762915*(alpha[22]*favg[24]+alpha[4]*favg[23])+0.07905694150420944*(alpha[10]*favg[22]+favg[10]*alpha[22])+0.0883883476483184*alpha[1]*favg[19]+0.1369306393762915*alpha[9]*favg[11]+0.07905694150420944*(alpha[4]*favg[9]+favg[4]*alpha[9]); 
  Ghat[41] += 0.0776323754260148*alpha[10]*favg[54]+0.0883883476483184*alpha[13]*favg[46]+0.1530931089239486*alpha[2]*favg[42]-0.5*fjump[41]+0.07905694150420944*alpha[17]*favg[41]+0.0883883476483184*(alpha[0]*favg[41]+alpha[6]*favg[40])+alpha[10]*(0.1767766952966368*favg[39]+0.07071067811865474*favg[38])+(0.1224744871391589*favg[24]+0.07071067811865474*favg[10])*alpha[38]+0.1369306393762915*(alpha[4]*favg[24]+alpha[22]*favg[23])+0.07905694150420944*(alpha[9]*favg[22]+favg[9]*alpha[22])+0.0883883476483184*alpha[2]*favg[19]+0.1369306393762915*alpha[10]*favg[11]+0.07905694150420944*(alpha[4]*favg[10]+favg[4]*alpha[10]); 
  Ghat[43] += 0.2338535866733712*alpha[43]*favg[53]+0.0776323754260148*alpha[12]*favg[51]+(0.1344632185501192*favg[26]+0.0776323754260148*favg[12])*alpha[51]+0.1369306393762915*alpha[43]*favg[49]+0.07071067811865474*alpha[12]*favg[47]+0.1976423537605236*alpha[16]*favg[45]-0.5*fjump[43]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[43]+(0.09780759955449389*favg[33]+0.07905694150420944*favg[20]+0.1976423537605236*favg[18]+0.05646924393157818*favg[16]+0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[43]+(0.1530931089239486*favg[30]+0.0883883476483184*favg[15])*alpha[37]+0.1767766952966368*alpha[12]*favg[35]+0.1530931089239486*alpha[5]*favg[33]+0.0883883476483184*alpha[13]*favg[31]+(0.1530931089239486*favg[27]+0.0883883476483184*favg[13])*alpha[31]+0.07905694150420944*alpha[9]*favg[28]+0.1369306393762915*alpha[1]*favg[26]+0.07905694150420944*alpha[6]*favg[25]+(0.1369306393762915*favg[21]+0.07905694150420944*favg[6])*alpha[25]+0.0883883476483184*alpha[5]*favg[16]+(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])*alpha[16]+0.07905694150420944*alpha[1]*favg[12]+(0.1369306393762915*favg[7]+0.07905694150420944*favg[1])*alpha[12]; 
  Ghat[44] += 0.2338535866733712*alpha[44]*favg[53]+0.0776323754260148*alpha[13]*favg[52]+(0.1344632185501192*favg[27]+0.0776323754260148*favg[13])*alpha[52]+0.1369306393762915*alpha[44]*favg[49]+0.07071067811865474*alpha[13]*favg[48]+0.1976423537605236*alpha[17]*favg[45]-0.5*fjump[44]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[44]+(0.09780759955449389*favg[34]+0.07905694150420944*favg[20]+0.1976423537605236*favg[18]+0.05646924393157818*favg[17]+0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[44]+(0.1530931089239486*favg[30]+0.0883883476483184*favg[15])*alpha[38]+0.1767766952966368*alpha[13]*favg[36]+0.1530931089239486*alpha[5]*favg[34]+0.0883883476483184*alpha[12]*favg[32]+(0.1530931089239486*favg[26]+0.0883883476483184*favg[12])*alpha[32]+0.07905694150420944*alpha[10]*favg[29]+0.1369306393762915*alpha[2]*favg[27]+0.07905694150420944*alpha[6]*favg[25]+(0.1369306393762915*favg[21]+0.07905694150420944*favg[6])*alpha[25]+0.0883883476483184*alpha[5]*favg[17]+(0.1530931089239486*favg[14]+0.0883883476483184*favg[5])*alpha[17]+0.07905694150420944*alpha[2]*favg[13]+(0.1369306393762915*favg[8]+0.07905694150420944*favg[2])*alpha[13]; 
  Ghat[46] += (-0.5*fjump[46])+0.0883883476483184*alpha[0]*favg[46]+0.1530931089239486*alpha[5]*favg[42]+0.0883883476483184*(alpha[13]*favg[41]+alpha[12]*favg[40])+0.1369306393762915*alpha[4]*favg[30]+0.07905694150420944*(alpha[10]*favg[29]+alpha[9]*favg[28])+0.0883883476483184*alpha[5]*favg[19]+0.07905694150420944*alpha[4]*favg[15]; 
  Ghat[47] += 0.0776323754260148*alpha[12]*favg[55]+0.0883883476483184*alpha[9]*favg[50]+0.1530931089239486*alpha[1]*favg[49]+0.0883883476483184*alpha[6]*favg[48]-0.5*fjump[47]+(0.07905694150420944*alpha[16]+0.0883883476483184*alpha[0])*favg[47]+alpha[12]*(0.1767766952966368*favg[45]+0.07071067811865474*favg[43])+(0.1224744871391589*favg[26]+0.07071067811865474*favg[12])*alpha[43]+0.1369306393762915*(alpha[25]*favg[27]+alpha[5]*favg[26])+0.07905694150420944*(alpha[13]*favg[25]+favg[13]*alpha[25])+0.0883883476483184*alpha[1]*favg[20]+0.1369306393762915*alpha[12]*favg[14]+0.07905694150420944*(alpha[5]*favg[12]+favg[5]*alpha[12]); 
  Ghat[48] += 0.0776323754260148*alpha[13]*favg[55]+0.0883883476483184*alpha[10]*favg[50]+0.1530931089239486*alpha[2]*favg[49]-0.5*fjump[48]+0.07905694150420944*alpha[17]*favg[48]+0.0883883476483184*(alpha[0]*favg[48]+alpha[6]*favg[47])+alpha[13]*(0.1767766952966368*favg[45]+0.07071067811865474*favg[44])+(0.1224744871391589*favg[27]+0.07071067811865474*favg[13])*alpha[44]+0.1369306393762915*(alpha[5]*favg[27]+alpha[25]*favg[26])+0.07905694150420944*(alpha[12]*favg[25]+favg[12]*alpha[25])+0.0883883476483184*alpha[2]*favg[20]+0.1369306393762915*alpha[13]*favg[14]+0.07905694150420944*(alpha[5]*favg[13]+favg[5]*alpha[13]); 
  Ghat[50] += (-0.5*fjump[50])+0.0883883476483184*alpha[0]*favg[50]+0.1530931089239486*alpha[4]*favg[49]+0.0883883476483184*(alpha[10]*favg[48]+alpha[9]*favg[47])+0.1369306393762915*alpha[5]*favg[30]+0.07905694150420944*(alpha[13]*favg[29]+alpha[12]*favg[28])+0.0883883476483184*alpha[4]*favg[20]+0.07905694150420944*alpha[5]*favg[15]; 
  Ghat[51] += 0.2338535866733712*alpha[51]*favg[53]-0.5*fjump[51]+(0.05270462766947297*alpha[16]+0.0883883476483184*alpha[0])*favg[51]+(0.09128709291752767*favg[33]+0.1976423537605236*favg[18]+0.05270462766947297*favg[16]+0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[51]+0.0776323754260148*alpha[12]*favg[43]+0.1344632185501192*favg[26]*alpha[43]+0.0776323754260148*(favg[12]*alpha[43]+alpha[9]*favg[37])+(0.1344632185501192*favg[23]+0.0776323754260148*favg[9])*alpha[37]+0.1735912687073533*alpha[16]*favg[35]+0.1344632185501192*alpha[1]*favg[33]+0.0776323754260148*alpha[6]*favg[31]+0.1344632185501192*favg[21]*alpha[31]+0.0776323754260148*(favg[6]*alpha[31]+alpha[1]*favg[16])+(0.1344632185501192*favg[7]+0.0776323754260148*favg[1])*alpha[16]; 
  Ghat[52] += 0.2338535866733712*alpha[52]*favg[53]-0.5*fjump[52]+(0.05270462766947297*alpha[17]+0.0883883476483184*alpha[0])*favg[52]+(0.09128709291752767*favg[34]+0.1976423537605236*favg[18]+0.05270462766947297*favg[17]+0.1530931089239486*favg[3]+0.0883883476483184*favg[0])*alpha[52]+0.0776323754260148*alpha[13]*favg[44]+0.1344632185501192*favg[27]*alpha[44]+0.0776323754260148*(favg[13]*alpha[44]+alpha[10]*favg[38])+(0.1344632185501192*favg[24]+0.0776323754260148*favg[10])*alpha[38]+0.1735912687073533*alpha[17]*favg[36]+0.1344632185501192*alpha[2]*favg[34]+0.0776323754260148*alpha[6]*favg[32]+0.1344632185501192*favg[21]*alpha[32]+0.0776323754260148*(favg[6]*alpha[32]+alpha[2]*favg[17])+(0.1344632185501192*favg[8]+0.0776323754260148*favg[2])*alpha[17]; 
  Ghat[54] += (-0.5*fjump[54])+0.0883883476483184*alpha[0]*favg[54]+0.1344632185501192*alpha[4]*favg[42]+0.0776323754260148*(alpha[10]*favg[41]+alpha[9]*favg[40]+alpha[4]*favg[19]); 
  Ghat[55] += (-0.5*fjump[55])+0.0883883476483184*alpha[0]*favg[55]+0.1344632185501192*alpha[5]*favg[49]+0.0776323754260148*(alpha[13]*favg[48]+alpha[12]*favg[47]+alpha[5]*favg[20]); 

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
  outr[16] += 0.5*Ghat[16]*dv10r; 
  outr[17] += 0.5*Ghat[17]*dv10r; 
  outr[18] += 1.118033988749895*Ghat[0]*dv10r; 
  outr[19] += 0.5*Ghat[19]*dv10r; 
  outr[20] += 0.5*Ghat[20]*dv10r; 
  outr[21] += -0.8660254037844386*Ghat[6]*dv10r; 
  outr[22] += 0.5*Ghat[22]*dv10r; 
  outr[23] += -0.8660254037844386*Ghat[9]*dv10r; 
  outr[24] += -0.8660254037844386*Ghat[10]*dv10r; 
  outr[25] += 0.5*Ghat[25]*dv10r; 
  outr[26] += -0.8660254037844386*Ghat[12]*dv10r; 
  outr[27] += -0.8660254037844386*Ghat[13]*dv10r; 
  outr[28] += 0.5*Ghat[28]*dv10r; 
  outr[29] += 0.5*Ghat[29]*dv10r; 
  outr[30] += -0.8660254037844386*Ghat[15]*dv10r; 
  outr[31] += 0.5*Ghat[31]*dv10r; 
  outr[32] += 0.5*Ghat[32]*dv10r; 
  outr[33] += -0.8660254037844387*Ghat[16]*dv10r; 
  outr[34] += -0.8660254037844387*Ghat[17]*dv10r; 
  outr[35] += 1.118033988749895*Ghat[1]*dv10r; 
  outr[36] += 1.118033988749895*Ghat[2]*dv10r; 
  outr[37] += 0.5*Ghat[37]*dv10r; 
  outr[38] += 0.5*Ghat[38]*dv10r; 
  outr[39] += 1.118033988749895*Ghat[4]*dv10r; 
  outr[40] += 0.5*Ghat[40]*dv10r; 
  outr[41] += 0.5*Ghat[41]*dv10r; 
  outr[42] += -0.8660254037844387*Ghat[19]*dv10r; 
  outr[43] += 0.5*Ghat[43]*dv10r; 
  outr[44] += 0.5*Ghat[44]*dv10r; 
  outr[45] += 1.118033988749895*Ghat[5]*dv10r; 
  outr[46] += 0.5*Ghat[46]*dv10r; 
  outr[47] += 0.5*Ghat[47]*dv10r; 
  outr[48] += 0.5*Ghat[48]*dv10r; 
  outr[49] += -0.8660254037844387*Ghat[20]*dv10r; 
  outr[50] += 0.5*Ghat[50]*dv10r; 
  outr[51] += 0.5*Ghat[51]*dv10r; 
  outr[52] += 0.5*Ghat[52]*dv10r; 
  outr[53] += -1.322875655532295*Ghat[0]*dv10r; 
  outr[54] += 0.5*Ghat[54]*dv10r; 
  outr[55] += 0.5*Ghat[55]*dv10r; 

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
  outl[16] += -0.5*Ghat[16]*dv10l; 
  outl[17] += -0.5*Ghat[17]*dv10l; 
  outl[18] += -1.118033988749895*Ghat[0]*dv10l; 
  outl[19] += -0.5*Ghat[19]*dv10l; 
  outl[20] += -0.5*Ghat[20]*dv10l; 
  outl[21] += -0.8660254037844386*Ghat[6]*dv10l; 
  outl[22] += -0.5*Ghat[22]*dv10l; 
  outl[23] += -0.8660254037844386*Ghat[9]*dv10l; 
  outl[24] += -0.8660254037844386*Ghat[10]*dv10l; 
  outl[25] += -0.5*Ghat[25]*dv10l; 
  outl[26] += -0.8660254037844386*Ghat[12]*dv10l; 
  outl[27] += -0.8660254037844386*Ghat[13]*dv10l; 
  outl[28] += -0.5*Ghat[28]*dv10l; 
  outl[29] += -0.5*Ghat[29]*dv10l; 
  outl[30] += -0.8660254037844386*Ghat[15]*dv10l; 
  outl[31] += -0.5*Ghat[31]*dv10l; 
  outl[32] += -0.5*Ghat[32]*dv10l; 
  outl[33] += -0.8660254037844387*Ghat[16]*dv10l; 
  outl[34] += -0.8660254037844387*Ghat[17]*dv10l; 
  outl[35] += -1.118033988749895*Ghat[1]*dv10l; 
  outl[36] += -1.118033988749895*Ghat[2]*dv10l; 
  outl[37] += -0.5*Ghat[37]*dv10l; 
  outl[38] += -0.5*Ghat[38]*dv10l; 
  outl[39] += -1.118033988749895*Ghat[4]*dv10l; 
  outl[40] += -0.5*Ghat[40]*dv10l; 
  outl[41] += -0.5*Ghat[41]*dv10l; 
  outl[42] += -0.8660254037844387*Ghat[19]*dv10l; 
  outl[43] += -0.5*Ghat[43]*dv10l; 
  outl[44] += -0.5*Ghat[44]*dv10l; 
  outl[45] += -1.118033988749895*Ghat[5]*dv10l; 
  outl[46] += -0.5*Ghat[46]*dv10l; 
  outl[47] += -0.5*Ghat[47]*dv10l; 
  outl[48] += -0.5*Ghat[48]*dv10l; 
  outl[49] += -0.8660254037844387*Ghat[20]*dv10l; 
  outl[50] += -0.5*Ghat[50]*dv10l; 
  outl[51] += -0.5*Ghat[51]*dv10l; 
  outl[52] += -0.5*Ghat[52]*dv10l; 
  outl[53] += -1.322875655532295*Ghat[0]*dv10l; 
  outl[54] += -0.5*Ghat[54]*dv10l; 
  outl[55] += -0.5*Ghat[55]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[3]; 
  double dv11r = 2/dxvr[3]; 
  const double *E1 = &EM[3]; 

  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[6]; 

  for(unsigned int i=0; i<6; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[6]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  double fjump[6]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(-1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  double alpha[6]; 

  alpha[0] = 2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]; 
  alpha[1] = 2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]; 
  alpha[2] = 2.828427124746191*B0[2]*wv3-2.828427124746191*B2[2]*wv1+2.828427124746191*E1[2]; 
  alpha[3] = -0.8164965809277261*B2[0]*dv1; 
  alpha[5] = 0.8164965809277261*B0[0]*dv3; 
  const double amid = 0.1767766952966368*alpha[0]; 
  Ghat[0] += 0.0883883476483184*alpha[5]*favg[5]-0.8660254037844386*fjump[4]+0.1530931089239486*alpha[0]*favg[4]+0.0883883476483184*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.0883883476483184*alpha[0]*favg[0]; 
  Ghat[1] += 0.1530931089239486*alpha[1]*favg[4]-0.5*fjump[1]+0.0883883476483184*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += 0.1530931089239486*alpha[2]*favg[4]-0.5*fjump[2]+0.0883883476483184*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] += 0.1530931089239486*alpha[3]*favg[4]-0.5*fjump[3]+0.0883883476483184*(alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[5] += (-0.5*fjump[5])+0.0883883476483184*alpha[0]*favg[5]+(0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[5]; 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += 0.5*Ghat[3]*dv11r; 
  outr[4] += -0.8660254037844386*Ghat[0]*dv11r; 
  outr[5] += 0.5*Ghat[5]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.5*Ghat[3]*dv11l; 
  outl[4] += -0.8660254037844386*Ghat[0]*dv11l; 
  outl[5] += -0.5*Ghat[5]*dv11l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[3]; 
  double dv11r = 2/dxvr[3]; 
  const double *E1 = &EM[6]; 

  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double Ghat[21]; 

  for(unsigned int i=0; i<21; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[21]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  double fjump[21]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(-1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(-1*fr[9]-fl[9]); 
  fjump[10] = amax*(-1*fr[10]-fl[10]); 
  fjump[11] = amax*(-1*fr[11]-fl[11]); 
  fjump[12] = amax*(1*fr[12]-fl[12]); 
  fjump[13] = amax*(1*fr[13]-fl[13]); 
  fjump[14] = amax*(1*fr[14]-fl[14]); 
  fjump[15] = amax*(-1*fr[15]-fl[15]); 
  fjump[16] = amax*(1*fr[16]-fl[16]); 
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(1*fr[18]-fl[18]); 
  fjump[19] = amax*(1*fr[19]-fl[19]); 
  fjump[20] = amax*(1*fr[20]-fl[20]); 
  double alpha[21]; 

  alpha[0] = 2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]; 
  alpha[1] = 2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]; 
  alpha[2] = 2.828427124746191*B0[2]*wv3-2.828427124746191*B2[2]*wv1+2.828427124746191*E1[2]; 
  alpha[3] = -0.8164965809277261*B2[0]*dv1; 
  alpha[5] = 0.8164965809277261*B0[0]*dv3; 
  alpha[6] = 2.828427124746191*B0[3]*wv3-2.828427124746191*B2[3]*wv1+2.828427124746191*E1[3]; 
  alpha[7] = -0.8164965809277261*B2[1]*dv1; 
  alpha[8] = -0.8164965809277261*B2[2]*dv1; 
  alpha[12] = 0.8164965809277261*B0[1]*dv3; 
  alpha[13] = 0.8164965809277261*B0[2]*dv3; 
  alpha[16] = 2.828427124746191*B0[4]*wv3-2.828427124746191*B2[4]*wv1+2.828427124746191*E1[4]; 
  alpha[17] = 2.828427124746191*B0[5]*wv3-2.828427124746191*B2[5]*wv1+2.828427124746191*E1[5]; 
  const double amid = (-0.1976423537605236*alpha[17])-0.1976423537605236*alpha[16]+0.1767766952966368*alpha[0]; 
  Ghat[0] += (-1.118033988749895*fjump[19])+0.1976423537605236*alpha[0]*favg[19]+0.0883883476483184*(alpha[17]*favg[17]+alpha[16]*favg[16])+0.1530931089239486*alpha[5]*favg[15]+0.0883883476483184*(alpha[13]*favg[13]+alpha[12]*favg[12])+0.1530931089239486*(alpha[3]*favg[11]+alpha[2]*favg[10]+alpha[1]*favg[9])+0.0883883476483184*(alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5])-0.8660254037844386*fjump[4]+0.1530931089239486*alpha[0]*favg[4]+0.0883883476483184*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.0883883476483184*alpha[0]*favg[0]; 
  Ghat[1] += alpha[1]*(0.1976423537605236*favg[19]+0.07905694150420944*favg[16])+(0.1369306393762915*favg[9]+0.07905694150420944*favg[1])*alpha[16]+0.1530931089239486*alpha[12]*favg[15]+0.0883883476483184*(alpha[5]*favg[12]+favg[5]*alpha[12])+0.1530931089239486*(alpha[7]*favg[11]+alpha[6]*favg[10])-0.8660254037844386*fjump[9]+0.1530931089239486*alpha[0]*favg[9]+0.0883883476483184*(alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6])+0.1530931089239486*alpha[1]*favg[4]-0.5*fjump[1]+0.0883883476483184*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += alpha[2]*(0.1976423537605236*favg[19]+0.07905694150420944*favg[17])+(0.1369306393762915*favg[10]+0.07905694150420944*favg[2])*alpha[17]+0.1530931089239486*alpha[13]*favg[15]+0.0883883476483184*(alpha[5]*favg[13]+favg[5]*alpha[13])+0.1530931089239486*alpha[8]*favg[11]-0.8660254037844386*fjump[10]+0.1530931089239486*(alpha[0]*favg[10]+alpha[6]*favg[9])+0.0883883476483184*(alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[1]*favg[6]+favg[1]*alpha[6])+0.1530931089239486*alpha[2]*favg[4]-0.5*fjump[2]+0.0883883476483184*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] += alpha[3]*(0.1976423537605236*favg[19]+0.07905694150420944*favg[18])+0.0883883476483184*alpha[5]*favg[14]-0.8660254037844386*fjump[11]+0.1530931089239486*(alpha[0]*favg[11]+alpha[8]*favg[10]+alpha[7]*favg[9])+0.0883883476483184*(alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[1]*favg[7]+favg[1]*alpha[7])+0.1530931089239486*alpha[3]*favg[4]-0.5*fjump[3]+0.0883883476483184*(alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[5] += alpha[5]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[19])-0.8660254037844386*fjump[15]+0.1530931089239486*alpha[0]*favg[15]+0.0883883476483184*(alpha[3]*favg[14]+alpha[2]*favg[13])+0.1530931089239486*favg[10]*alpha[13]+0.0883883476483184*(favg[2]*alpha[13]+alpha[1]*favg[12])+(0.1530931089239486*favg[9]+0.0883883476483184*favg[1])*alpha[12]-0.5*fjump[5]+0.0883883476483184*alpha[0]*favg[5]+(0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[5]; 
  Ghat[6] += 0.1976423537605236*alpha[6]*favg[19]+0.07905694150420944*(alpha[6]*favg[17]+favg[6]*alpha[17]+alpha[6]*favg[16]+favg[6]*alpha[16])+0.0883883476483184*(alpha[12]*favg[13]+favg[12]*alpha[13])+0.1530931089239486*(alpha[1]*favg[10]+alpha[2]*favg[9])+0.0883883476483184*(alpha[7]*favg[8]+favg[7]*alpha[8])-0.5*fjump[6]+0.0883883476483184*alpha[0]*favg[6]+0.1530931089239486*favg[4]*alpha[6]+0.0883883476483184*(favg[0]*alpha[6]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[7] += 0.1976423537605236*alpha[7]*favg[19]+0.07905694150420944*(alpha[7]*(favg[18]+favg[16])+favg[7]*alpha[16])+0.0883883476483184*alpha[12]*favg[14]+0.1530931089239486*(alpha[1]*favg[11]+alpha[3]*favg[9])+0.0883883476483184*(alpha[6]*favg[8]+favg[6]*alpha[8])-0.5*fjump[7]+0.0883883476483184*alpha[0]*favg[7]+0.1530931089239486*favg[4]*alpha[7]+0.0883883476483184*(favg[0]*alpha[7]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[8] += 0.1976423537605236*alpha[8]*favg[19]+0.07905694150420944*(alpha[8]*(favg[18]+favg[17])+favg[8]*alpha[17])+0.0883883476483184*alpha[13]*favg[14]+0.1530931089239486*(alpha[2]*favg[11]+alpha[3]*favg[10])-0.5*fjump[8]+0.0883883476483184*alpha[0]*favg[8]+0.1530931089239486*favg[4]*alpha[8]+0.0883883476483184*(favg[0]*alpha[8]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[12] += alpha[12]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[19])+0.07905694150420944*(alpha[12]*favg[16]+favg[12]*alpha[16])+0.1530931089239486*alpha[1]*favg[15]+0.0883883476483184*(alpha[7]*favg[14]+alpha[6]*favg[13]+favg[6]*alpha[13])-0.5*fjump[12]+0.0883883476483184*alpha[0]*favg[12]+(0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[12]+0.1530931089239486*alpha[5]*favg[9]+0.0883883476483184*(alpha[1]*favg[5]+favg[1]*alpha[5]); 
  Ghat[13] += alpha[13]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[19])+0.07905694150420944*(alpha[13]*favg[17]+favg[13]*alpha[17])+0.1530931089239486*alpha[2]*favg[15]+0.0883883476483184*alpha[8]*favg[14]-0.5*fjump[13]+0.0883883476483184*alpha[0]*favg[13]+0.1530931089239486*favg[4]*alpha[13]+0.0883883476483184*(favg[0]*alpha[13]+alpha[6]*favg[12]+favg[6]*alpha[12])+0.1530931089239486*alpha[5]*favg[10]+0.0883883476483184*(alpha[2]*favg[5]+favg[2]*alpha[5]); 
  Ghat[14] += 0.1530931089239486*alpha[3]*favg[15]-0.5*fjump[14]+0.0883883476483184*(alpha[0]*favg[14]+alpha[8]*favg[13]+favg[8]*alpha[13]+alpha[7]*favg[12]+favg[7]*alpha[12])+0.1530931089239486*alpha[5]*favg[11]+0.0883883476483184*(alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[16] += 0.1976423537605236*alpha[16]*favg[19]-0.5*fjump[16]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[16]+(0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[16]+0.07905694150420944*alpha[12]*favg[12]+0.1369306393762915*alpha[1]*favg[9]+0.07905694150420944*(alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[1]*favg[1]); 
  Ghat[17] += 0.1976423537605236*alpha[17]*favg[19]-0.5*fjump[17]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[17]+(0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[17]+0.07905694150420944*alpha[13]*favg[13]+0.1369306393762915*alpha[2]*favg[10]+0.07905694150420944*(alpha[8]*favg[8]+alpha[6]*favg[6]+alpha[2]*favg[2]); 
  Ghat[18] += (-0.5*fjump[18])+0.0883883476483184*alpha[0]*favg[18]+0.1369306393762915*alpha[3]*favg[11]+0.07905694150420944*(alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[3]*favg[3]); 
  Ghat[20] += (-0.5*fjump[20])+0.0883883476483184*alpha[0]*favg[20]+0.1369306393762915*alpha[5]*favg[15]+0.07905694150420944*(alpha[13]*favg[13]+alpha[12]*favg[12]+alpha[5]*favg[5]); 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += 0.5*Ghat[3]*dv11r; 
  outr[4] += -0.8660254037844386*Ghat[0]*dv11r; 
  outr[5] += 0.5*Ghat[5]*dv11r; 
  outr[6] += 0.5*Ghat[6]*dv11r; 
  outr[7] += 0.5*Ghat[7]*dv11r; 
  outr[8] += 0.5*Ghat[8]*dv11r; 
  outr[9] += -0.8660254037844386*Ghat[1]*dv11r; 
  outr[10] += -0.8660254037844386*Ghat[2]*dv11r; 
  outr[11] += -0.8660254037844386*Ghat[3]*dv11r; 
  outr[12] += 0.5*Ghat[12]*dv11r; 
  outr[13] += 0.5*Ghat[13]*dv11r; 
  outr[14] += 0.5*Ghat[14]*dv11r; 
  outr[15] += -0.8660254037844386*Ghat[5]*dv11r; 
  outr[16] += 0.5*Ghat[16]*dv11r; 
  outr[17] += 0.5*Ghat[17]*dv11r; 
  outr[18] += 0.5*Ghat[18]*dv11r; 
  outr[19] += 1.118033988749895*Ghat[0]*dv11r; 
  outr[20] += 0.5*Ghat[20]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.5*Ghat[3]*dv11l; 
  outl[4] += -0.8660254037844386*Ghat[0]*dv11l; 
  outl[5] += -0.5*Ghat[5]*dv11l; 
  outl[6] += -0.5*Ghat[6]*dv11l; 
  outl[7] += -0.5*Ghat[7]*dv11l; 
  outl[8] += -0.5*Ghat[8]*dv11l; 
  outl[9] += -0.8660254037844386*Ghat[1]*dv11l; 
  outl[10] += -0.8660254037844386*Ghat[2]*dv11l; 
  outl[11] += -0.8660254037844386*Ghat[3]*dv11l; 
  outl[12] += -0.5*Ghat[12]*dv11l; 
  outl[13] += -0.5*Ghat[13]*dv11l; 
  outl[14] += -0.5*Ghat[14]*dv11l; 
  outl[15] += -0.8660254037844386*Ghat[5]*dv11l; 
  outl[16] += -0.5*Ghat[16]*dv11l; 
  outl[17] += -0.5*Ghat[17]*dv11l; 
  outl[18] += -0.5*Ghat[18]*dv11l; 
  outl[19] += -1.118033988749895*Ghat[0]*dv11l; 
  outl[20] += -0.5*Ghat[20]*dv11l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[3]; 
  double dv11r = 2/dxvr[3]; 
  const double *E1 = &EM[10]; 

  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *B0 = &EM[30]; 
  const double *B1 = &EM[40]; 
  const double *B2 = &EM[50]; 

  double Ghat[56]; 

  for(unsigned int i=0; i<56; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[56]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  favg[21] = 1*fr[21]+fl[21]; 
  favg[22] = -1*fr[22]+fl[22]; 
  favg[23] = -1*fr[23]+fl[23]; 
  favg[24] = -1*fr[24]+fl[24]; 
  favg[25] = 1*fr[25]+fl[25]; 
  favg[26] = 1*fr[26]+fl[26]; 
  favg[27] = 1*fr[27]+fl[27]; 
  favg[28] = -1*fr[28]+fl[28]; 
  favg[29] = -1*fr[29]+fl[29]; 
  favg[30] = -1*fr[30]+fl[30]; 
  favg[31] = 1*fr[31]+fl[31]; 
  favg[32] = 1*fr[32]+fl[32]; 
  favg[33] = 1*fr[33]+fl[33]; 
  favg[34] = 1*fr[34]+fl[34]; 
  favg[35] = 1*fr[35]+fl[35]; 
  favg[36] = 1*fr[36]+fl[36]; 
  favg[37] = -1*fr[37]+fl[37]; 
  favg[38] = -1*fr[38]+fl[38]; 
  favg[39] = -1*fr[39]+fl[39]; 
  favg[40] = 1*fr[40]+fl[40]; 
  favg[41] = 1*fr[41]+fl[41]; 
  favg[42] = 1*fr[42]+fl[42]; 
  favg[43] = 1*fr[43]+fl[43]; 
  favg[44] = 1*fr[44]+fl[44]; 
  favg[45] = 1*fr[45]+fl[45]; 
  favg[46] = 1*fr[46]+fl[46]; 
  favg[47] = 1*fr[47]+fl[47]; 
  favg[48] = 1*fr[48]+fl[48]; 
  favg[49] = 1*fr[49]+fl[49]; 
  favg[50] = -1*fr[50]+fl[50]; 
  favg[51] = 1*fr[51]+fl[51]; 
  favg[52] = 1*fr[52]+fl[52]; 
  favg[53] = 1*fr[53]+fl[53]; 
  favg[54] = -1*fr[54]+fl[54]; 
  favg[55] = 1*fr[55]+fl[55]; 
  double fjump[56]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(-1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(-1*fr[9]-fl[9]); 
  fjump[10] = amax*(-1*fr[10]-fl[10]); 
  fjump[11] = amax*(-1*fr[11]-fl[11]); 
  fjump[12] = amax*(1*fr[12]-fl[12]); 
  fjump[13] = amax*(1*fr[13]-fl[13]); 
  fjump[14] = amax*(1*fr[14]-fl[14]); 
  fjump[15] = amax*(-1*fr[15]-fl[15]); 
  fjump[16] = amax*(1*fr[16]-fl[16]); 
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(1*fr[18]-fl[18]); 
  fjump[19] = amax*(1*fr[19]-fl[19]); 
  fjump[20] = amax*(1*fr[20]-fl[20]); 
  fjump[21] = amax*(1*fr[21]-fl[21]); 
  fjump[22] = amax*(-1*fr[22]-fl[22]); 
  fjump[23] = amax*(-1*fr[23]-fl[23]); 
  fjump[24] = amax*(-1*fr[24]-fl[24]); 
  fjump[25] = amax*(1*fr[25]-fl[25]); 
  fjump[26] = amax*(1*fr[26]-fl[26]); 
  fjump[27] = amax*(1*fr[27]-fl[27]); 
  fjump[28] = amax*(-1*fr[28]-fl[28]); 
  fjump[29] = amax*(-1*fr[29]-fl[29]); 
  fjump[30] = amax*(-1*fr[30]-fl[30]); 
  fjump[31] = amax*(1*fr[31]-fl[31]); 
  fjump[32] = amax*(1*fr[32]-fl[32]); 
  fjump[33] = amax*(1*fr[33]-fl[33]); 
  fjump[34] = amax*(1*fr[34]-fl[34]); 
  fjump[35] = amax*(1*fr[35]-fl[35]); 
  fjump[36] = amax*(1*fr[36]-fl[36]); 
  fjump[37] = amax*(-1*fr[37]-fl[37]); 
  fjump[38] = amax*(-1*fr[38]-fl[38]); 
  fjump[39] = amax*(-1*fr[39]-fl[39]); 
  fjump[40] = amax*(1*fr[40]-fl[40]); 
  fjump[41] = amax*(1*fr[41]-fl[41]); 
  fjump[42] = amax*(1*fr[42]-fl[42]); 
  fjump[43] = amax*(1*fr[43]-fl[43]); 
  fjump[44] = amax*(1*fr[44]-fl[44]); 
  fjump[45] = amax*(1*fr[45]-fl[45]); 
  fjump[46] = amax*(1*fr[46]-fl[46]); 
  fjump[47] = amax*(1*fr[47]-fl[47]); 
  fjump[48] = amax*(1*fr[48]-fl[48]); 
  fjump[49] = amax*(1*fr[49]-fl[49]); 
  fjump[50] = amax*(-1*fr[50]-fl[50]); 
  fjump[51] = amax*(1*fr[51]-fl[51]); 
  fjump[52] = amax*(1*fr[52]-fl[52]); 
  fjump[53] = amax*(1*fr[53]-fl[53]); 
  fjump[54] = amax*(-1*fr[54]-fl[54]); 
  fjump[55] = amax*(1*fr[55]-fl[55]); 
  double alpha[56]; 

  alpha[0] = 2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]; 
  alpha[1] = 2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]; 
  alpha[2] = 2.828427124746191*B0[2]*wv3-2.828427124746191*B2[2]*wv1+2.828427124746191*E1[2]; 
  alpha[3] = -0.8164965809277261*B2[0]*dv1; 
  alpha[5] = 0.8164965809277261*B0[0]*dv3; 
  alpha[6] = 2.828427124746191*B0[3]*wv3-2.828427124746191*B2[3]*wv1+2.828427124746191*E1[3]; 
  alpha[7] = -0.8164965809277261*B2[1]*dv1; 
  alpha[8] = -0.8164965809277261*B2[2]*dv1; 
  alpha[12] = 0.8164965809277261*B0[1]*dv3; 
  alpha[13] = 0.8164965809277261*B0[2]*dv3; 
  alpha[16] = 2.828427124746191*B0[4]*wv3-2.828427124746191*B2[4]*wv1+2.828427124746191*E1[4]; 
  alpha[17] = 2.828427124746191*B0[5]*wv3-2.828427124746191*B2[5]*wv1+2.828427124746191*E1[5]; 
  alpha[21] = -0.8164965809277261*B2[3]*dv1; 
  alpha[25] = 0.8164965809277261*B0[3]*dv3; 
  alpha[31] = 2.828427124746191*B0[6]*wv3-2.828427124746191*B2[6]*wv1+2.828427124746191*E1[6]; 
  alpha[32] = 2.828427124746191*B0[7]*wv3-2.828427124746191*B2[7]*wv1+2.828427124746191*E1[7]; 
  alpha[33] = -0.816496580927726*B2[4]*dv1; 
  alpha[34] = -0.816496580927726*B2[5]*dv1; 
  alpha[43] = 0.816496580927726*B0[4]*dv3; 
  alpha[44] = 0.816496580927726*B0[5]*dv3; 
  alpha[51] = 2.828427124746191*B0[8]*wv3-2.828427124746191*B2[8]*wv1+2.828427124746191*E1[8]; 
  alpha[52] = 2.828427124746191*B0[9]*wv3-2.828427124746191*B2[9]*wv1+2.828427124746191*E1[9]; 
  const double amid = (-0.1976423537605236*alpha[17])-0.1976423537605236*alpha[16]+0.1767766952966368*alpha[0]; 
  Ghat[0] += (-1.322875655532295*fjump[54])+0.2338535866733712*alpha[0]*favg[54]+0.0883883476483184*(alpha[52]*favg[52]+alpha[51]*favg[51])+0.1976423537605236*alpha[5]*favg[46]+0.0883883476483184*(alpha[44]*favg[44]+alpha[43]*favg[43])+0.1976423537605236*(alpha[3]*favg[42]+alpha[2]*favg[41]+alpha[1]*favg[40])+0.1530931089239486*(alpha[17]*favg[38]+alpha[16]*favg[37])+0.0883883476483184*(alpha[34]*favg[34]+alpha[33]*favg[33]+alpha[32]*favg[32]+alpha[31]*favg[31])+0.1530931089239486*(alpha[13]*favg[29]+alpha[12]*favg[28])+0.0883883476483184*alpha[25]*favg[25]+0.1530931089239486*(alpha[8]*favg[24]+alpha[7]*favg[23]+alpha[6]*favg[22])+0.0883883476483184*alpha[21]*favg[21]-1.118033988749895*fjump[19]+0.1976423537605236*alpha[0]*favg[19]+0.0883883476483184*(alpha[17]*favg[17]+alpha[16]*favg[16])+0.1530931089239486*alpha[5]*favg[15]+0.0883883476483184*(alpha[13]*favg[13]+alpha[12]*favg[12])+0.1530931089239486*(alpha[3]*favg[11]+alpha[2]*favg[10]+alpha[1]*favg[9])+0.0883883476483184*(alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5])-0.8660254037844386*fjump[4]+0.1530931089239486*alpha[0]*favg[4]+0.0883883476483184*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.0883883476483184*alpha[0]*favg[0]; 
  Ghat[1] += 0.2338535866733712*alpha[1]*favg[54]+0.0776323754260148*alpha[16]*favg[51]+(0.1344632185501192*favg[37]+0.0776323754260148*favg[16])*alpha[51]+alpha[12]*(0.1976423537605236*favg[46]+0.07905694150420944*favg[43])+(0.1369306393762915*favg[28]+0.07905694150420944*favg[12])*alpha[43]+0.1976423537605236*(alpha[7]*favg[42]+alpha[6]*favg[41])-1.118033988749895*fjump[40]+(0.1767766952966368*alpha[16]+0.1976423537605236*alpha[0])*favg[40]+0.1530931089239486*alpha[32]*favg[38]+0.1369306393762915*alpha[1]*favg[37]+0.07905694150420944*alpha[7]*favg[33]+(0.1369306393762915*favg[23]+0.07905694150420944*favg[7])*alpha[33]+0.0883883476483184*(alpha[17]*favg[32]+favg[17]*alpha[32])+0.07905694150420944*alpha[6]*favg[31]+(0.1369306393762915*favg[22]+0.07905694150420944*favg[6])*alpha[31]+0.1530931089239486*(alpha[25]*favg[29]+alpha[5]*favg[28])+0.0883883476483184*(alpha[13]*favg[25]+favg[13]*alpha[25])+0.1530931089239486*(alpha[21]*favg[24]+alpha[3]*favg[23]+alpha[2]*favg[22])+0.0883883476483184*(alpha[8]*favg[21]+favg[8]*alpha[21])+alpha[1]*(0.1976423537605236*favg[19]+0.07905694150420944*favg[16])+(0.1369306393762915*favg[9]+0.07905694150420944*favg[1])*alpha[16]+0.1530931089239486*alpha[12]*favg[15]+0.0883883476483184*(alpha[5]*favg[12]+favg[5]*alpha[12])+0.1530931089239486*(alpha[7]*favg[11]+alpha[6]*favg[10])-0.8660254037844386*fjump[9]+0.1530931089239486*alpha[0]*favg[9]+0.0883883476483184*(alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6])+0.1530931089239486*alpha[1]*favg[4]-0.5*fjump[1]+0.0883883476483184*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += 0.2338535866733712*alpha[2]*favg[54]+0.0776323754260148*alpha[17]*favg[52]+(0.1344632185501192*favg[38]+0.0776323754260148*favg[17])*alpha[52]+alpha[13]*(0.1976423537605236*favg[46]+0.07905694150420944*favg[44])+(0.1369306393762915*favg[29]+0.07905694150420944*favg[13])*alpha[44]+0.1976423537605236*alpha[8]*favg[42]-1.118033988749895*fjump[41]+0.1767766952966368*alpha[17]*favg[41]+0.1976423537605236*(alpha[0]*favg[41]+alpha[6]*favg[40])+0.1369306393762915*alpha[2]*favg[38]+0.1530931089239486*alpha[31]*favg[37]+0.07905694150420944*alpha[8]*favg[34]+0.1369306393762915*favg[24]*alpha[34]+0.07905694150420944*(favg[8]*alpha[34]+alpha[6]*favg[32])+(0.1369306393762915*favg[22]+0.07905694150420944*favg[6])*alpha[32]+0.0883883476483184*(alpha[16]*favg[31]+favg[16]*alpha[31])+0.1530931089239486*(alpha[5]*favg[29]+alpha[25]*favg[28])+0.0883883476483184*(alpha[12]*favg[25]+favg[12]*alpha[25])+0.1530931089239486*(alpha[3]*favg[24]+alpha[21]*favg[23]+alpha[1]*favg[22])+0.0883883476483184*(alpha[7]*favg[21]+favg[7]*alpha[21])+alpha[2]*(0.1976423537605236*favg[19]+0.07905694150420944*favg[17])+(0.1369306393762915*favg[10]+0.07905694150420944*favg[2])*alpha[17]+0.1530931089239486*alpha[13]*favg[15]+0.0883883476483184*(alpha[5]*favg[13]+favg[5]*alpha[13])+0.1530931089239486*alpha[8]*favg[11]-0.8660254037844386*fjump[10]+0.1530931089239486*(alpha[0]*favg[10]+alpha[6]*favg[9])+0.0883883476483184*(alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[1]*favg[6]+favg[1]*alpha[6])+0.1530931089239486*alpha[2]*favg[4]-0.5*fjump[2]+0.0883883476483184*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] += 0.2338535866733712*alpha[3]*favg[54]-1.118033988749895*fjump[42]+0.1976423537605236*(alpha[0]*favg[42]+alpha[8]*favg[41]+alpha[7]*favg[40])+0.1369306393762915*alpha[3]*favg[39]+0.1530931089239486*(alpha[34]*favg[38]+alpha[33]*favg[37])+0.07905694150420944*(alpha[8]*favg[36]+alpha[7]*favg[35])+0.0883883476483184*(alpha[17]*favg[34]+favg[17]*alpha[34]+alpha[16]*favg[33]+favg[16]*alpha[33])+0.1530931089239486*alpha[5]*favg[30]+0.0883883476483184*(alpha[13]*favg[27]+alpha[12]*favg[26])+0.1530931089239486*(alpha[2]*favg[24]+alpha[1]*favg[23]+alpha[21]*favg[22])+0.0883883476483184*(alpha[6]*favg[21]+favg[6]*alpha[21])+alpha[3]*(0.1976423537605236*favg[19]+0.07905694150420944*favg[18])+0.0883883476483184*alpha[5]*favg[14]-0.8660254037844386*fjump[11]+0.1530931089239486*(alpha[0]*favg[11]+alpha[8]*favg[10]+alpha[7]*favg[9])+0.0883883476483184*(alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[1]*favg[7]+favg[1]*alpha[7])+0.1530931089239486*alpha[3]*favg[4]-0.5*fjump[3]+0.0883883476483184*(alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[5] += alpha[5]*(0.2338535866733712*favg[54]+0.1369306393762915*favg[50])+0.07905694150420944*(alpha[13]*favg[48]+alpha[12]*favg[47])-1.118033988749895*fjump[46]+0.1976423537605236*alpha[0]*favg[46]+0.0883883476483184*alpha[17]*favg[44]+0.1530931089239486*favg[38]*alpha[44]+0.0883883476483184*(favg[17]*alpha[44]+alpha[16]*favg[43])+(0.1530931089239486*favg[37]+0.0883883476483184*favg[16])*alpha[43]+0.1976423537605236*(alpha[13]*favg[41]+alpha[12]*favg[40])+0.1530931089239486*(alpha[3]*favg[30]+alpha[2]*favg[29]+alpha[1]*favg[28])+0.0883883476483184*(alpha[8]*favg[27]+alpha[7]*favg[26]+alpha[6]*favg[25])+(0.1530931089239486*favg[22]+0.0883883476483184*favg[6])*alpha[25]+alpha[5]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[19])-0.8660254037844386*fjump[15]+0.1530931089239486*alpha[0]*favg[15]+0.0883883476483184*(alpha[3]*favg[14]+alpha[2]*favg[13])+0.1530931089239486*favg[10]*alpha[13]+0.0883883476483184*(favg[2]*alpha[13]+alpha[1]*favg[12])+(0.1530931089239486*favg[9]+0.0883883476483184*favg[1])*alpha[12]-0.5*fjump[5]+0.0883883476483184*alpha[0]*favg[5]+(0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[5]; 
  Ghat[6] += 0.2338535866733712*alpha[6]*favg[54]+0.0776323754260148*(alpha[32]*favg[52]+favg[32]*alpha[52]+alpha[31]*favg[51]+favg[31]*alpha[51])+0.1976423537605236*alpha[25]*favg[46]+0.07905694150420944*(alpha[25]*favg[44]+favg[25]*alpha[44]+alpha[25]*favg[43]+favg[25]*alpha[43])+0.1976423537605236*alpha[21]*favg[42]+(0.1767766952966368*alpha[32]+0.1976423537605236*alpha[1])*favg[41]+(0.1767766952966368*alpha[31]+0.1976423537605236*alpha[2])*favg[40]+0.1369306393762915*alpha[6]*(favg[38]+favg[37])+0.07905694150420944*(alpha[21]*favg[34]+favg[21]*alpha[34]+alpha[21]*favg[33]+favg[21]*alpha[33])+(0.07071067811865474*alpha[31]+0.07905694150420944*alpha[2])*favg[32]+(0.07071067811865474*favg[31]+0.1369306393762915*favg[10])*alpha[32]+0.07905694150420944*(favg[2]*alpha[32]+alpha[1]*favg[31])+(0.1369306393762915*favg[9]+0.07905694150420944*favg[1])*alpha[31]+0.1530931089239486*(alpha[12]*favg[29]+alpha[13]*favg[28])+0.0883883476483184*alpha[5]*favg[25]+(0.1530931089239486*favg[15]+0.0883883476483184*favg[5])*alpha[25]+0.1530931089239486*(alpha[7]*favg[24]+alpha[8]*favg[23])-0.8660254037844386*fjump[22]+(0.1369306393762915*(alpha[17]+alpha[16])+0.1530931089239486*alpha[0])*favg[22]+0.0883883476483184*alpha[3]*favg[21]+(0.1530931089239486*favg[11]+0.0883883476483184*favg[3])*alpha[21]+0.1976423537605236*alpha[6]*favg[19]+0.07905694150420944*(alpha[6]*favg[17]+favg[6]*alpha[17]+alpha[6]*favg[16]+favg[6]*alpha[16])+0.0883883476483184*(alpha[12]*favg[13]+favg[12]*alpha[13])+0.1530931089239486*(alpha[1]*favg[10]+alpha[2]*favg[9])+0.0883883476483184*(alpha[7]*favg[8]+favg[7]*alpha[8])-0.5*fjump[6]+0.0883883476483184*alpha[0]*favg[6]+0.1530931089239486*favg[4]*alpha[6]+0.0883883476483184*(favg[0]*alpha[6]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[7] += 0.2338535866733712*alpha[7]*favg[54]+0.0776323754260148*(alpha[33]*favg[51]+favg[33]*alpha[51])+0.07905694150420944*favg[26]*alpha[43]+0.1976423537605236*(alpha[1]*favg[42]+alpha[21]*favg[41])+(0.1767766952966368*alpha[33]+0.1976423537605236*alpha[3])*favg[40]+0.1369306393762915*alpha[7]*(favg[39]+favg[37])+0.07905694150420944*alpha[21]*favg[36]+(0.07071067811865474*alpha[33]+0.07905694150420944*alpha[3])*favg[35]+0.0883883476483184*(alpha[32]*favg[34]+favg[32]*alpha[34])+0.07905694150420944*alpha[1]*favg[33]+0.1369306393762915*favg[9]*alpha[33]+0.07905694150420944*(favg[1]*alpha[33]+alpha[21]*favg[31]+favg[21]*alpha[31])+0.1530931089239486*alpha[12]*favg[30]+0.0883883476483184*(alpha[25]*favg[27]+alpha[5]*favg[26])+0.1530931089239486*alpha[6]*favg[24]-0.8660254037844386*fjump[23]+0.1369306393762915*alpha[16]*favg[23]+0.1530931089239486*(alpha[0]*favg[23]+alpha[8]*favg[22])+0.0883883476483184*alpha[2]*favg[21]+(0.1530931089239486*favg[10]+0.0883883476483184*favg[2])*alpha[21]+0.1976423537605236*alpha[7]*favg[19]+0.07905694150420944*(alpha[7]*(favg[18]+favg[16])+favg[7]*alpha[16])+0.0883883476483184*alpha[12]*favg[14]+0.1530931089239486*(alpha[1]*favg[11]+alpha[3]*favg[9])+0.0883883476483184*(alpha[6]*favg[8]+favg[6]*alpha[8])-0.5*fjump[7]+0.0883883476483184*alpha[0]*favg[7]+0.1530931089239486*favg[4]*alpha[7]+0.0883883476483184*(favg[0]*alpha[7]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[8] += 0.2338535866733712*alpha[8]*favg[54]+0.0776323754260148*(alpha[34]*favg[52]+favg[34]*alpha[52])+0.07905694150420944*favg[27]*alpha[44]+0.1976423537605236*alpha[2]*favg[42]+0.1767766952966368*alpha[34]*favg[41]+0.1976423537605236*(alpha[3]*favg[41]+alpha[21]*favg[40])+0.1369306393762915*alpha[8]*(favg[39]+favg[38])+0.07071067811865474*alpha[34]*favg[36]+0.07905694150420944*(alpha[3]*favg[36]+alpha[21]*favg[35]+alpha[2]*favg[34])+(0.1369306393762915*favg[10]+0.07905694150420944*favg[2])*alpha[34]+0.0883883476483184*(alpha[31]*favg[33]+favg[31]*alpha[33])+0.07905694150420944*(alpha[21]*favg[32]+favg[21]*alpha[32])+0.1530931089239486*alpha[13]*favg[30]+0.0883883476483184*(alpha[5]*favg[27]+alpha[25]*favg[26])-0.8660254037844386*fjump[24]+0.1369306393762915*alpha[17]*favg[24]+0.1530931089239486*(alpha[0]*favg[24]+alpha[6]*favg[23]+alpha[7]*favg[22])+0.0883883476483184*alpha[1]*favg[21]+(0.1530931089239486*favg[9]+0.0883883476483184*favg[1])*alpha[21]+0.1976423537605236*alpha[8]*favg[19]+0.07905694150420944*(alpha[8]*(favg[18]+favg[17])+favg[8]*alpha[17])+0.0883883476483184*alpha[13]*favg[14]+0.1530931089239486*(alpha[2]*favg[11]+alpha[3]*favg[10])-0.5*fjump[8]+0.0883883476483184*alpha[0]*favg[8]+0.1530931089239486*favg[4]*alpha[8]+0.0883883476483184*(favg[0]*alpha[8]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[12] += 0.2338535866733712*alpha[12]*favg[54]+0.0776323754260148*(alpha[43]*favg[51]+favg[43]*alpha[51])+0.1369306393762915*alpha[12]*favg[50]+0.07905694150420944*alpha[25]*favg[48]+(0.07071067811865474*alpha[43]+0.07905694150420944*alpha[5])*favg[47]+0.1976423537605236*alpha[1]*favg[46]+0.0883883476483184*(alpha[32]*favg[44]+favg[32]*alpha[44])+0.07905694150420944*alpha[1]*favg[43]+(0.1767766952966368*favg[40]+0.1369306393762915*favg[9]+0.07905694150420944*favg[1])*alpha[43]+0.1976423537605236*(alpha[25]*favg[41]+alpha[5]*favg[40])+0.1369306393762915*alpha[12]*favg[37]+0.07905694150420944*(favg[26]*alpha[33]+alpha[25]*favg[31]+favg[25]*alpha[31])+0.1530931089239486*(alpha[7]*favg[30]+alpha[6]*favg[29])-0.8660254037844386*fjump[28]+(0.1369306393762915*alpha[16]+0.1530931089239486*alpha[0])*favg[28]+0.0883883476483184*(alpha[21]*favg[27]+alpha[3]*favg[26]+alpha[2]*favg[25])+(0.1530931089239486*favg[10]+0.0883883476483184*favg[2])*alpha[25]+0.1530931089239486*alpha[13]*favg[22]+alpha[12]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[19])+0.07905694150420944*(alpha[12]*favg[16]+favg[12]*alpha[16])+0.1530931089239486*alpha[1]*favg[15]+0.0883883476483184*(alpha[7]*favg[14]+alpha[6]*favg[13]+favg[6]*alpha[13])-0.5*fjump[12]+0.0883883476483184*alpha[0]*favg[12]+(0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[12]+0.1530931089239486*alpha[5]*favg[9]+0.0883883476483184*(alpha[1]*favg[5]+favg[1]*alpha[5]); 
  Ghat[13] += 0.2338535866733712*alpha[13]*favg[54]+0.0776323754260148*(alpha[44]*favg[52]+favg[44]*alpha[52])+0.1369306393762915*alpha[13]*favg[50]+0.07071067811865474*alpha[44]*favg[48]+0.07905694150420944*(alpha[5]*favg[48]+alpha[25]*favg[47])+alpha[2]*(0.1976423537605236*favg[46]+0.07905694150420944*favg[44])+(0.1767766952966368*favg[41]+0.1369306393762915*favg[10]+0.07905694150420944*favg[2])*alpha[44]+0.0883883476483184*(alpha[31]*favg[43]+favg[31]*alpha[43])+0.1976423537605236*(alpha[5]*favg[41]+alpha[25]*favg[40])+0.1369306393762915*alpha[13]*favg[38]+0.07905694150420944*(favg[27]*alpha[34]+alpha[25]*favg[32]+favg[25]*alpha[32])+0.1530931089239486*alpha[8]*favg[30]-0.8660254037844386*fjump[29]+0.1369306393762915*alpha[17]*favg[29]+0.1530931089239486*(alpha[0]*favg[29]+alpha[6]*favg[28])+0.0883883476483184*(alpha[3]*favg[27]+alpha[21]*favg[26]+alpha[1]*favg[25])+(0.1530931089239486*favg[9]+0.0883883476483184*favg[1])*alpha[25]+0.1530931089239486*alpha[12]*favg[22]+alpha[13]*(0.07905694150420944*favg[20]+0.1976423537605236*favg[19])+0.07905694150420944*(alpha[13]*favg[17]+favg[13]*alpha[17])+0.1530931089239486*alpha[2]*favg[15]+0.0883883476483184*alpha[8]*favg[14]-0.5*fjump[13]+0.0883883476483184*alpha[0]*favg[13]+0.1530931089239486*favg[4]*alpha[13]+0.0883883476483184*(favg[0]*alpha[13]+alpha[6]*favg[12]+favg[6]*alpha[12])+0.1530931089239486*alpha[5]*favg[10]+0.0883883476483184*(alpha[2]*favg[5]+favg[2]*alpha[5]); 
  Ghat[14] += 0.07905694150420944*alpha[5]*favg[49]+alpha[3]*(0.1976423537605236*favg[46]+0.07905694150420944*favg[45])+0.0883883476483184*(alpha[34]*favg[44]+favg[34]*alpha[44]+alpha[33]*favg[43]+favg[33]*alpha[43])+0.1976423537605236*alpha[5]*favg[42]-0.8660254037844386*fjump[30]+0.1530931089239486*(alpha[0]*favg[30]+alpha[8]*favg[29]+alpha[7]*favg[28])+0.0883883476483184*(alpha[2]*favg[27]+alpha[1]*favg[26]+alpha[21]*favg[25]+favg[21]*alpha[25])+0.1530931089239486*(alpha[13]*favg[24]+alpha[12]*favg[23]+alpha[3]*favg[15])-0.5*fjump[14]+0.0883883476483184*(alpha[0]*favg[14]+alpha[8]*favg[13]+favg[8]*alpha[13]+alpha[7]*favg[12]+favg[7]*alpha[12])+0.1530931089239486*alpha[5]*favg[11]+0.0883883476483184*(alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[16] += 0.2338535866733712*alpha[16]*favg[54]+(0.05270462766947297*alpha[51]+0.0776323754260148*alpha[1])*favg[51]+(0.1735912687073533*favg[40]+0.1344632185501192*favg[9]+0.0776323754260148*favg[1])*alpha[51]+0.1976423537605236*alpha[43]*favg[46]+(0.05646924393157818*alpha[43]+0.0883883476483184*alpha[5])*favg[43]+(0.1530931089239486*favg[15]+0.0883883476483184*favg[5])*alpha[43]+0.1976423537605236*(alpha[33]*favg[42]+alpha[31]*favg[41])+0.1767766952966368*alpha[1]*favg[40]-0.8660254037844386*fjump[37]+(0.09780759955449389*alpha[16]+0.1530931089239486*alpha[0])*favg[37]+(0.05646924393157818*alpha[33]+0.0883883476483184*alpha[3])*favg[33]+(0.1530931089239486*favg[11]+0.0883883476483184*favg[3])*alpha[33]+0.07905694150420944*alpha[32]*favg[32]+(0.05646924393157818*alpha[31]+0.0883883476483184*alpha[2])*favg[31]+(0.1530931089239486*favg[10]+0.0883883476483184*favg[2])*alpha[31]+0.1369306393762915*alpha[12]*favg[28]+0.07905694150420944*alpha[25]*favg[25]+0.1369306393762915*(alpha[7]*favg[23]+alpha[6]*favg[22])+0.07905694150420944*alpha[21]*favg[21]+0.1976423537605236*alpha[16]*favg[19]-0.5*fjump[16]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[16]+(0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[16]+0.07905694150420944*alpha[12]*favg[12]+0.1369306393762915*alpha[1]*favg[9]+0.07905694150420944*(alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[1]*favg[1]); 
  Ghat[17] += 0.2338535866733712*alpha[17]*favg[54]+(0.05270462766947297*alpha[52]+0.0776323754260148*alpha[2])*favg[52]+(0.1735912687073533*favg[41]+0.1344632185501192*favg[10]+0.0776323754260148*favg[2])*alpha[52]+0.1976423537605236*alpha[44]*favg[46]+(0.05646924393157818*alpha[44]+0.0883883476483184*alpha[5])*favg[44]+(0.1530931089239486*favg[15]+0.0883883476483184*favg[5])*alpha[44]+0.1976423537605236*alpha[34]*favg[42]+0.1767766952966368*alpha[2]*favg[41]+0.1976423537605236*alpha[32]*favg[40]-0.8660254037844386*fjump[38]+(0.09780759955449389*alpha[17]+0.1530931089239486*alpha[0])*favg[38]+(0.05646924393157818*alpha[34]+0.0883883476483184*alpha[3])*favg[34]+(0.1530931089239486*favg[11]+0.0883883476483184*favg[3])*alpha[34]+(0.05646924393157818*alpha[32]+0.0883883476483184*alpha[1])*favg[32]+(0.1530931089239486*favg[9]+0.0883883476483184*favg[1])*alpha[32]+0.07905694150420944*alpha[31]*favg[31]+0.1369306393762915*alpha[13]*favg[29]+0.07905694150420944*alpha[25]*favg[25]+0.1369306393762915*(alpha[8]*favg[24]+alpha[6]*favg[22])+0.07905694150420944*alpha[21]*favg[21]+0.1976423537605236*alpha[17]*favg[19]-0.5*fjump[17]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[17]+(0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[17]+0.07905694150420944*alpha[13]*favg[13]+0.1369306393762915*alpha[2]*favg[10]+0.07905694150420944*(alpha[8]*favg[8]+alpha[6]*favg[6]+alpha[2]*favg[2]); 
  Ghat[18] += 0.0776323754260148*alpha[3]*favg[53]+0.0883883476483184*alpha[5]*favg[45]+0.1767766952966368*alpha[3]*favg[42]-0.8660254037844386*fjump[39]+0.1530931089239486*alpha[0]*favg[39]+0.0883883476483184*(alpha[2]*favg[36]+alpha[1]*favg[35])+0.07905694150420944*(alpha[34]*favg[34]+alpha[33]*favg[33])+0.1369306393762915*(alpha[8]*favg[24]+alpha[7]*favg[23])+0.07905694150420944*alpha[21]*favg[21]-0.5*fjump[18]+0.0883883476483184*alpha[0]*favg[18]+0.1369306393762915*alpha[3]*favg[11]+0.07905694150420944*(alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[3]*favg[3]); 
  Ghat[20] += 0.0776323754260148*alpha[5]*favg[55]-0.8660254037844386*fjump[50]+0.1530931089239486*alpha[0]*favg[50]+0.0883883476483184*(alpha[3]*favg[49]+alpha[2]*favg[48]+alpha[1]*favg[47])+0.1767766952966368*alpha[5]*favg[46]+0.07905694150420944*(alpha[44]*favg[44]+alpha[43]*favg[43])+0.1369306393762915*(alpha[13]*favg[29]+alpha[12]*favg[28])+0.07905694150420944*alpha[25]*favg[25]-0.5*fjump[20]+0.0883883476483184*alpha[0]*favg[20]+0.1369306393762915*alpha[5]*favg[15]+0.07905694150420944*(alpha[13]*favg[13]+alpha[12]*favg[12]+alpha[5]*favg[5]); 
  Ghat[21] += 0.2338535866733712*alpha[21]*favg[54]+0.1976423537605236*(alpha[6]*favg[42]+alpha[7]*favg[41]+alpha[8]*favg[40])+0.1369306393762915*alpha[21]*(favg[39]+favg[38]+favg[37])+0.07905694150420944*(alpha[7]*favg[36]+alpha[8]*favg[35]+alpha[6]*favg[34])+0.1369306393762915*favg[22]*alpha[34]+0.07905694150420944*(favg[6]*alpha[34]+alpha[6]*favg[33])+0.1369306393762915*favg[22]*alpha[33]+0.07905694150420944*(favg[6]*alpha[33]+alpha[8]*favg[32])+0.1369306393762915*favg[24]*alpha[32]+0.07905694150420944*(favg[8]*alpha[32]+alpha[7]*favg[31])+(0.1369306393762915*favg[23]+0.07905694150420944*favg[7])*alpha[31]+0.1530931089239486*alpha[25]*favg[30]+0.0883883476483184*(alpha[12]*favg[27]+alpha[13]*favg[26]+favg[14]*alpha[25])+0.1530931089239486*(alpha[1]*favg[24]+alpha[2]*favg[23]+alpha[3]*favg[22])-0.5*fjump[21]+(0.07905694150420944*(alpha[17]+alpha[16])+0.0883883476483184*alpha[0])*favg[21]+(0.1976423537605236*favg[19]+0.07905694150420944*(favg[18]+favg[17]+favg[16])+0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[21]+0.1530931089239486*(alpha[6]*favg[11]+alpha[7]*favg[10]+alpha[8]*favg[9])+0.0883883476483184*(alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[3]*favg[6]+favg[3]*alpha[6]); 
  Ghat[25] += alpha[25]*(0.2338535866733712*favg[54]+0.1369306393762915*favg[50])+0.07905694150420944*(alpha[12]*favg[48]+alpha[13]*favg[47])+alpha[6]*(0.1976423537605236*favg[46]+0.07905694150420944*favg[44])+0.1369306393762915*favg[22]*alpha[44]+0.07905694150420944*(favg[6]*alpha[44]+alpha[6]*favg[43])+(0.1369306393762915*favg[22]+0.07905694150420944*favg[6])*alpha[43]+0.1976423537605236*(alpha[12]*favg[41]+alpha[13]*favg[40])+0.1369306393762915*alpha[25]*(favg[38]+favg[37])+0.07905694150420944*alpha[13]*favg[32]+0.1369306393762915*favg[29]*alpha[32]+0.07905694150420944*(favg[13]*alpha[32]+alpha[12]*favg[31])+(0.1369306393762915*favg[28]+0.07905694150420944*favg[12])*alpha[31]+0.1530931089239486*(alpha[21]*favg[30]+alpha[1]*favg[29]+alpha[2]*favg[28])+0.0883883476483184*(alpha[7]*favg[27]+alpha[8]*favg[26])-0.5*fjump[25]+(0.07905694150420944*(alpha[17]+alpha[16])+0.0883883476483184*alpha[0])*favg[25]+(0.07905694150420944*favg[20]+0.1976423537605236*favg[19]+0.07905694150420944*(favg[17]+favg[16])+0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[25]+0.1530931089239486*alpha[5]*favg[22]+0.0883883476483184*favg[14]*alpha[21]+0.1530931089239486*alpha[6]*favg[15]+0.0883883476483184*alpha[1]*favg[13]+0.1530931089239486*favg[9]*alpha[13]+0.0883883476483184*(favg[1]*alpha[13]+alpha[2]*favg[12])+0.1530931089239486*favg[10]*alpha[12]+0.0883883476483184*(favg[2]*alpha[12]+alpha[5]*favg[6]+favg[5]*alpha[6]); 
  Ghat[26] += 0.07905694150420944*alpha[12]*favg[49]+alpha[7]*(0.1976423537605236*favg[46]+0.07905694150420944*(favg[45]+favg[43]))+(0.1369306393762915*favg[23]+0.07905694150420944*favg[7])*alpha[43]+alpha[12]*(0.1976423537605236*favg[42]+0.07905694150420944*favg[33])+(0.1369306393762915*favg[28]+0.07905694150420944*favg[12])*alpha[33]+0.1530931089239486*(alpha[1]*favg[30]+alpha[21]*favg[29]+alpha[3]*favg[28])+0.0883883476483184*alpha[6]*favg[27]-0.5*fjump[26]+0.07905694150420944*alpha[16]*favg[26]+0.0883883476483184*(alpha[0]*favg[26]+alpha[8]*favg[25])+(0.1530931089239486*favg[24]+0.0883883476483184*favg[8])*alpha[25]+0.1530931089239486*alpha[5]*favg[23]+0.0883883476483184*(alpha[13]*favg[21]+favg[13]*alpha[21])+0.1530931089239486*alpha[7]*favg[15]+0.0883883476483184*(alpha[1]*favg[14]+alpha[3]*favg[12])+0.1530931089239486*favg[11]*alpha[12]+0.0883883476483184*(favg[3]*alpha[12]+alpha[5]*favg[7]+favg[5]*alpha[7]); 
  Ghat[27] += 0.07905694150420944*alpha[13]*favg[49]+alpha[8]*(0.1976423537605236*favg[46]+0.07905694150420944*(favg[45]+favg[44]))+(0.1369306393762915*favg[24]+0.07905694150420944*favg[8])*alpha[44]+alpha[13]*(0.1976423537605236*favg[42]+0.07905694150420944*favg[34])+(0.1369306393762915*favg[29]+0.07905694150420944*favg[13])*alpha[34]+0.1530931089239486*(alpha[2]*favg[30]+alpha[3]*favg[29]+alpha[21]*favg[28])-0.5*fjump[27]+0.07905694150420944*alpha[17]*favg[27]+0.0883883476483184*(alpha[0]*favg[27]+alpha[6]*favg[26]+alpha[7]*favg[25])+(0.1530931089239486*favg[23]+0.0883883476483184*favg[7])*alpha[25]+0.1530931089239486*alpha[5]*favg[24]+0.0883883476483184*(alpha[12]*favg[21]+favg[12]*alpha[21])+0.1530931089239486*alpha[8]*favg[15]+0.0883883476483184*(alpha[2]*favg[14]+alpha[3]*favg[13])+0.1530931089239486*favg[11]*alpha[13]+0.0883883476483184*(favg[3]*alpha[13]+alpha[5]*favg[8]+favg[5]*alpha[8]); 
  Ghat[31] += 0.2338535866733712*alpha[31]*favg[54]+0.0776323754260148*alpha[6]*favg[51]+(0.1344632185501192*favg[22]+0.0776323754260148*favg[6])*alpha[51]+0.0883883476483184*alpha[13]*favg[43]+(0.1530931089239486*favg[29]+0.0883883476483184*favg[13])*alpha[43]+0.1976423537605236*alpha[16]*favg[41]+0.1767766952966368*alpha[6]*favg[40]+0.1369306393762915*alpha[31]*favg[38]+(0.09780759955449389*alpha[31]+0.1530931089239486*alpha[2])*favg[37]+0.0883883476483184*alpha[8]*favg[33]+(0.1530931089239486*favg[24]+0.0883883476483184*favg[8])*alpha[33]+0.07071067811865474*alpha[6]*favg[32]+(0.1224744871391589*favg[22]+0.07071067811865474*favg[6])*alpha[32]-0.5*fjump[31]+(0.07905694150420944*alpha[17]+0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[31]+(0.1976423537605236*favg[19]+0.07905694150420944*favg[17]+0.05646924393157818*favg[16]+0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[31]+0.1369306393762915*alpha[25]*favg[28]+0.07905694150420944*(alpha[12]*favg[25]+favg[12]*alpha[25])+0.1369306393762915*(alpha[21]*favg[23]+alpha[1]*favg[22])+0.07905694150420944*(alpha[7]*favg[21]+favg[7]*alpha[21])+0.0883883476483184*alpha[2]*favg[16]+(0.1530931089239486*favg[10]+0.0883883476483184*favg[2])*alpha[16]+0.1369306393762915*alpha[6]*favg[9]+0.07905694150420944*(alpha[1]*favg[6]+favg[1]*alpha[6]); 
  Ghat[32] += 0.2338535866733712*alpha[32]*favg[54]+0.0776323754260148*alpha[6]*favg[52]+(0.1344632185501192*favg[22]+0.0776323754260148*favg[6])*alpha[52]+0.0883883476483184*alpha[12]*favg[44]+(0.1530931089239486*favg[28]+0.0883883476483184*favg[12])*alpha[44]+0.1767766952966368*alpha[6]*favg[41]+0.1976423537605236*alpha[17]*favg[40]+(0.09780759955449389*alpha[32]+0.1530931089239486*alpha[1])*favg[38]+0.1369306393762915*alpha[32]*favg[37]+0.0883883476483184*alpha[7]*favg[34]+(0.1530931089239486*favg[23]+0.0883883476483184*favg[7])*alpha[34]-0.5*fjump[32]+(0.05646924393157818*alpha[17]+0.07905694150420944*alpha[16]+0.0883883476483184*alpha[0])*favg[32]+(0.1976423537605236*favg[19]+0.05646924393157818*favg[17]+0.07905694150420944*favg[16]+0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[32]+0.07071067811865474*alpha[6]*favg[31]+(0.1224744871391589*favg[22]+0.07071067811865474*favg[6])*alpha[31]+0.1369306393762915*alpha[25]*favg[29]+0.07905694150420944*(alpha[13]*favg[25]+favg[13]*alpha[25])+0.1369306393762915*(alpha[21]*favg[24]+alpha[2]*favg[22])+0.07905694150420944*(alpha[8]*favg[21]+favg[8]*alpha[21])+0.0883883476483184*alpha[1]*favg[17]+(0.1530931089239486*favg[9]+0.0883883476483184*favg[1])*alpha[17]+0.1369306393762915*alpha[6]*favg[10]+0.07905694150420944*(alpha[2]*favg[6]+favg[2]*alpha[6]); 
  Ghat[33] += 0.2338535866733712*alpha[33]*favg[54]+0.0776323754260148*alpha[7]*favg[51]+(0.1344632185501192*favg[23]+0.0776323754260148*favg[7])*alpha[51]+(0.1530931089239486*favg[30]+0.0883883476483184*favg[14])*alpha[43]+0.1976423537605236*alpha[16]*favg[42]+0.1767766952966368*alpha[7]*favg[40]+0.1369306393762915*alpha[33]*favg[39]+(0.09780759955449389*alpha[33]+0.1530931089239486*alpha[3])*favg[37]+0.07071067811865474*alpha[7]*favg[35]-0.5*fjump[33]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[33]+(0.1976423537605236*favg[19]+0.07905694150420944*favg[18]+0.05646924393157818*favg[16]+0.1530931089239486*favg[4])*alpha[33]+0.0883883476483184*(favg[0]*alpha[33]+alpha[8]*favg[31])+(0.1530931089239486*favg[24]+0.0883883476483184*favg[8])*alpha[31]+0.07905694150420944*alpha[12]*favg[26]+0.1369306393762915*(alpha[1]*favg[23]+alpha[21]*favg[22])+0.07905694150420944*(alpha[6]*favg[21]+favg[6]*alpha[21])+0.0883883476483184*alpha[3]*favg[16]+(0.1530931089239486*favg[11]+0.0883883476483184*favg[3])*alpha[16]+0.1369306393762915*alpha[7]*favg[9]+0.07905694150420944*(alpha[1]*favg[7]+favg[1]*alpha[7]); 
  Ghat[34] += 0.2338535866733712*alpha[34]*favg[54]+0.0776323754260148*alpha[8]*favg[52]+(0.1344632185501192*favg[24]+0.0776323754260148*favg[8])*alpha[52]+(0.1530931089239486*favg[30]+0.0883883476483184*favg[14])*alpha[44]+0.1976423537605236*alpha[17]*favg[42]+0.1767766952966368*alpha[8]*favg[41]+0.1369306393762915*alpha[34]*favg[39]+(0.09780759955449389*alpha[34]+0.1530931089239486*alpha[3])*favg[38]+0.07071067811865474*alpha[8]*favg[36]-0.5*fjump[34]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[34]+(0.1976423537605236*favg[19]+0.07905694150420944*favg[18]+0.05646924393157818*favg[17]+0.1530931089239486*favg[4])*alpha[34]+0.0883883476483184*(favg[0]*alpha[34]+alpha[7]*favg[32])+(0.1530931089239486*favg[23]+0.0883883476483184*favg[7])*alpha[32]+0.07905694150420944*alpha[13]*favg[27]+0.1369306393762915*(alpha[2]*favg[24]+alpha[21]*favg[22])+0.07905694150420944*(alpha[6]*favg[21]+favg[6]*alpha[21])+0.0883883476483184*alpha[3]*favg[17]+(0.1530931089239486*favg[11]+0.0883883476483184*favg[3])*alpha[17]+0.1369306393762915*alpha[8]*favg[10]+0.07905694150420944*(alpha[2]*favg[8]+favg[2]*alpha[8]); 
  Ghat[35] += 0.0776323754260148*alpha[7]*favg[53]+0.0883883476483184*alpha[12]*favg[45]+0.1767766952966368*alpha[7]*favg[42]+0.1530931089239486*alpha[1]*favg[39]+0.0883883476483184*alpha[6]*favg[36]-0.5*fjump[35]+(0.07905694150420944*alpha[16]+0.0883883476483184*alpha[0])*favg[35]+0.07071067811865474*alpha[7]*favg[33]+(0.1224744871391589*favg[23]+0.07071067811865474*favg[7])*alpha[33]+0.1369306393762915*(alpha[21]*favg[24]+alpha[3]*favg[23])+0.07905694150420944*(alpha[8]*favg[21]+favg[8]*alpha[21])+0.0883883476483184*alpha[1]*favg[18]+0.1369306393762915*alpha[7]*favg[11]+0.07905694150420944*(alpha[3]*favg[7]+favg[3]*alpha[7]); 
  Ghat[36] += 0.0776323754260148*alpha[8]*favg[53]+0.0883883476483184*alpha[13]*favg[45]+0.1767766952966368*alpha[8]*favg[42]+0.1530931089239486*alpha[2]*favg[39]-0.5*fjump[36]+0.07905694150420944*alpha[17]*favg[36]+0.0883883476483184*(alpha[0]*favg[36]+alpha[6]*favg[35])+0.07071067811865474*alpha[8]*favg[34]+(0.1224744871391589*favg[24]+0.07071067811865474*favg[8])*alpha[34]+0.1369306393762915*(alpha[3]*favg[24]+alpha[21]*favg[23])+0.07905694150420944*(alpha[7]*favg[21]+favg[7]*alpha[21])+0.0883883476483184*alpha[2]*favg[18]+0.1369306393762915*alpha[8]*favg[11]+0.07905694150420944*(alpha[3]*favg[8]+favg[3]*alpha[8]); 
  Ghat[43] += 0.2338535866733712*alpha[43]*favg[54]+0.0776323754260148*alpha[12]*favg[51]+(0.1344632185501192*favg[28]+0.0776323754260148*favg[12])*alpha[51]+0.1369306393762915*alpha[43]*favg[50]+0.07071067811865474*alpha[12]*favg[47]+0.1976423537605236*alpha[16]*favg[46]-0.5*fjump[43]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[43]+(0.09780759955449389*favg[37]+0.07905694150420944*favg[20]+0.1976423537605236*favg[19]+0.05646924393157818*favg[16]+0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[43]+0.1767766952966368*alpha[12]*favg[40]+0.1530931089239486*(alpha[5]*favg[37]+favg[30]*alpha[33])+0.0883883476483184*(favg[14]*alpha[33]+alpha[13]*favg[31])+(0.1530931089239486*favg[29]+0.0883883476483184*favg[13])*alpha[31]+0.1369306393762915*alpha[1]*favg[28]+0.07905694150420944*(alpha[7]*favg[26]+alpha[6]*favg[25])+(0.1369306393762915*favg[22]+0.07905694150420944*favg[6])*alpha[25]+0.0883883476483184*alpha[5]*favg[16]+(0.1530931089239486*favg[15]+0.0883883476483184*favg[5])*alpha[16]+0.07905694150420944*alpha[1]*favg[12]+(0.1369306393762915*favg[9]+0.07905694150420944*favg[1])*alpha[12]; 
  Ghat[44] += 0.2338535866733712*alpha[44]*favg[54]+0.0776323754260148*alpha[13]*favg[52]+(0.1344632185501192*favg[29]+0.0776323754260148*favg[13])*alpha[52]+0.1369306393762915*alpha[44]*favg[50]+0.07071067811865474*alpha[13]*favg[48]+0.1976423537605236*alpha[17]*favg[46]-0.5*fjump[44]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[44]+(0.09780759955449389*favg[38]+0.07905694150420944*favg[20]+0.1976423537605236*favg[19]+0.05646924393157818*favg[17]+0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[44]+0.1767766952966368*alpha[13]*favg[41]+0.1530931089239486*(alpha[5]*favg[38]+favg[30]*alpha[34])+0.0883883476483184*(favg[14]*alpha[34]+alpha[12]*favg[32])+(0.1530931089239486*favg[28]+0.0883883476483184*favg[12])*alpha[32]+0.1369306393762915*alpha[2]*favg[29]+0.07905694150420944*(alpha[8]*favg[27]+alpha[6]*favg[25])+(0.1369306393762915*favg[22]+0.07905694150420944*favg[6])*alpha[25]+0.0883883476483184*alpha[5]*favg[17]+(0.1530931089239486*favg[15]+0.0883883476483184*favg[5])*alpha[17]+0.07905694150420944*alpha[2]*favg[13]+(0.1369306393762915*favg[10]+0.07905694150420944*favg[2])*alpha[13]; 
  Ghat[45] += (-0.5*fjump[45])+0.0883883476483184*alpha[0]*favg[45]+0.1530931089239486*alpha[5]*favg[39]+0.0883883476483184*(alpha[13]*favg[36]+alpha[12]*favg[35])+0.1369306393762915*alpha[3]*favg[30]+0.07905694150420944*(alpha[8]*favg[27]+alpha[7]*favg[26])+0.0883883476483184*alpha[5]*favg[18]+0.07905694150420944*alpha[3]*favg[14]; 
  Ghat[47] += 0.0776323754260148*alpha[12]*favg[55]+0.1530931089239486*alpha[1]*favg[50]+0.0883883476483184*(alpha[7]*favg[49]+alpha[6]*favg[48])-0.5*fjump[47]+(0.07905694150420944*alpha[16]+0.0883883476483184*alpha[0])*favg[47]+alpha[12]*(0.1767766952966368*favg[46]+0.07071067811865474*favg[43])+(0.1224744871391589*favg[28]+0.07071067811865474*favg[12])*alpha[43]+0.1369306393762915*(alpha[25]*favg[29]+alpha[5]*favg[28])+0.07905694150420944*(alpha[13]*favg[25]+favg[13]*alpha[25])+0.0883883476483184*alpha[1]*favg[20]+0.1369306393762915*alpha[12]*favg[15]+0.07905694150420944*(alpha[5]*favg[12]+favg[5]*alpha[12]); 
  Ghat[48] += 0.0776323754260148*alpha[13]*favg[55]+0.1530931089239486*alpha[2]*favg[50]+0.0883883476483184*alpha[8]*favg[49]-0.5*fjump[48]+0.07905694150420944*alpha[17]*favg[48]+0.0883883476483184*(alpha[0]*favg[48]+alpha[6]*favg[47])+alpha[13]*(0.1767766952966368*favg[46]+0.07071067811865474*favg[44])+(0.1224744871391589*favg[29]+0.07071067811865474*favg[13])*alpha[44]+0.1369306393762915*(alpha[5]*favg[29]+alpha[25]*favg[28])+0.07905694150420944*(alpha[12]*favg[25]+favg[12]*alpha[25])+0.0883883476483184*alpha[2]*favg[20]+0.1369306393762915*alpha[13]*favg[15]+0.07905694150420944*(alpha[5]*favg[13]+favg[5]*alpha[13]); 
  Ghat[49] += 0.1530931089239486*alpha[3]*favg[50]-0.5*fjump[49]+0.0883883476483184*(alpha[0]*favg[49]+alpha[8]*favg[48]+alpha[7]*favg[47])+0.1369306393762915*alpha[5]*favg[30]+0.07905694150420944*(alpha[13]*favg[27]+alpha[12]*favg[26])+0.0883883476483184*alpha[3]*favg[20]+0.07905694150420944*alpha[5]*favg[14]; 
  Ghat[51] += 0.2338535866733712*alpha[51]*favg[54]-0.5*fjump[51]+(0.05270462766947297*alpha[16]+0.0883883476483184*alpha[0])*favg[51]+(0.09128709291752767*favg[37]+0.1976423537605236*favg[19]+0.05270462766947297*favg[16]+0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[51]+0.0776323754260148*alpha[12]*favg[43]+(0.1344632185501192*favg[28]+0.0776323754260148*favg[12])*alpha[43]+0.1735912687073533*alpha[16]*favg[40]+0.1344632185501192*alpha[1]*favg[37]+0.0776323754260148*alpha[7]*favg[33]+0.1344632185501192*favg[23]*alpha[33]+0.0776323754260148*(favg[7]*alpha[33]+alpha[6]*favg[31])+0.1344632185501192*favg[22]*alpha[31]+0.0776323754260148*(favg[6]*alpha[31]+alpha[1]*favg[16])+(0.1344632185501192*favg[9]+0.0776323754260148*favg[1])*alpha[16]; 
  Ghat[52] += 0.2338535866733712*alpha[52]*favg[54]-0.5*fjump[52]+(0.05270462766947297*alpha[17]+0.0883883476483184*alpha[0])*favg[52]+(0.09128709291752767*favg[38]+0.1976423537605236*favg[19]+0.05270462766947297*favg[17]+0.1530931089239486*favg[4]+0.0883883476483184*favg[0])*alpha[52]+0.0776323754260148*alpha[13]*favg[44]+(0.1344632185501192*favg[29]+0.0776323754260148*favg[13])*alpha[44]+0.1735912687073533*alpha[17]*favg[41]+0.1344632185501192*alpha[2]*favg[38]+0.0776323754260148*alpha[8]*favg[34]+0.1344632185501192*favg[24]*alpha[34]+0.0776323754260148*(favg[8]*alpha[34]+alpha[6]*favg[32])+0.1344632185501192*favg[22]*alpha[32]+0.0776323754260148*(favg[6]*alpha[32]+alpha[2]*favg[17])+(0.1344632185501192*favg[10]+0.0776323754260148*favg[2])*alpha[17]; 
  Ghat[53] += (-0.5*fjump[53])+0.0883883476483184*alpha[0]*favg[53]+0.1344632185501192*alpha[3]*favg[39]+0.0776323754260148*(alpha[8]*favg[36]+alpha[7]*favg[35]+alpha[3]*favg[18]); 
  Ghat[55] += (-0.5*fjump[55])+0.0883883476483184*alpha[0]*favg[55]+0.1344632185501192*alpha[5]*favg[50]+0.0776323754260148*(alpha[13]*favg[48]+alpha[12]*favg[47]+alpha[5]*favg[20]); 

  outr[0] += 0.5*Ghat[0]*dv11r; 
  outr[1] += 0.5*Ghat[1]*dv11r; 
  outr[2] += 0.5*Ghat[2]*dv11r; 
  outr[3] += 0.5*Ghat[3]*dv11r; 
  outr[4] += -0.8660254037844386*Ghat[0]*dv11r; 
  outr[5] += 0.5*Ghat[5]*dv11r; 
  outr[6] += 0.5*Ghat[6]*dv11r; 
  outr[7] += 0.5*Ghat[7]*dv11r; 
  outr[8] += 0.5*Ghat[8]*dv11r; 
  outr[9] += -0.8660254037844386*Ghat[1]*dv11r; 
  outr[10] += -0.8660254037844386*Ghat[2]*dv11r; 
  outr[11] += -0.8660254037844386*Ghat[3]*dv11r; 
  outr[12] += 0.5*Ghat[12]*dv11r; 
  outr[13] += 0.5*Ghat[13]*dv11r; 
  outr[14] += 0.5*Ghat[14]*dv11r; 
  outr[15] += -0.8660254037844386*Ghat[5]*dv11r; 
  outr[16] += 0.5*Ghat[16]*dv11r; 
  outr[17] += 0.5*Ghat[17]*dv11r; 
  outr[18] += 0.5*Ghat[18]*dv11r; 
  outr[19] += 1.118033988749895*Ghat[0]*dv11r; 
  outr[20] += 0.5*Ghat[20]*dv11r; 
  outr[21] += 0.5*Ghat[21]*dv11r; 
  outr[22] += -0.8660254037844386*Ghat[6]*dv11r; 
  outr[23] += -0.8660254037844386*Ghat[7]*dv11r; 
  outr[24] += -0.8660254037844386*Ghat[8]*dv11r; 
  outr[25] += 0.5*Ghat[25]*dv11r; 
  outr[26] += 0.5*Ghat[26]*dv11r; 
  outr[27] += 0.5*Ghat[27]*dv11r; 
  outr[28] += -0.8660254037844386*Ghat[12]*dv11r; 
  outr[29] += -0.8660254037844386*Ghat[13]*dv11r; 
  outr[30] += -0.8660254037844386*Ghat[14]*dv11r; 
  outr[31] += 0.5*Ghat[31]*dv11r; 
  outr[32] += 0.5*Ghat[32]*dv11r; 
  outr[33] += 0.5*Ghat[33]*dv11r; 
  outr[34] += 0.5*Ghat[34]*dv11r; 
  outr[35] += 0.5*Ghat[35]*dv11r; 
  outr[36] += 0.5*Ghat[36]*dv11r; 
  outr[37] += -0.8660254037844387*Ghat[16]*dv11r; 
  outr[38] += -0.8660254037844387*Ghat[17]*dv11r; 
  outr[39] += -0.8660254037844387*Ghat[18]*dv11r; 
  outr[40] += 1.118033988749895*Ghat[1]*dv11r; 
  outr[41] += 1.118033988749895*Ghat[2]*dv11r; 
  outr[42] += 1.118033988749895*Ghat[3]*dv11r; 
  outr[43] += 0.5*Ghat[43]*dv11r; 
  outr[44] += 0.5*Ghat[44]*dv11r; 
  outr[45] += 0.5*Ghat[45]*dv11r; 
  outr[46] += 1.118033988749895*Ghat[5]*dv11r; 
  outr[47] += 0.5*Ghat[47]*dv11r; 
  outr[48] += 0.5*Ghat[48]*dv11r; 
  outr[49] += 0.5*Ghat[49]*dv11r; 
  outr[50] += -0.8660254037844387*Ghat[20]*dv11r; 
  outr[51] += 0.5*Ghat[51]*dv11r; 
  outr[52] += 0.5*Ghat[52]*dv11r; 
  outr[53] += 0.5*Ghat[53]*dv11r; 
  outr[54] += -1.322875655532295*Ghat[0]*dv11r; 
  outr[55] += 0.5*Ghat[55]*dv11r; 

  outl[0] += -0.5*Ghat[0]*dv11l; 
  outl[1] += -0.5*Ghat[1]*dv11l; 
  outl[2] += -0.5*Ghat[2]*dv11l; 
  outl[3] += -0.5*Ghat[3]*dv11l; 
  outl[4] += -0.8660254037844386*Ghat[0]*dv11l; 
  outl[5] += -0.5*Ghat[5]*dv11l; 
  outl[6] += -0.5*Ghat[6]*dv11l; 
  outl[7] += -0.5*Ghat[7]*dv11l; 
  outl[8] += -0.5*Ghat[8]*dv11l; 
  outl[9] += -0.8660254037844386*Ghat[1]*dv11l; 
  outl[10] += -0.8660254037844386*Ghat[2]*dv11l; 
  outl[11] += -0.8660254037844386*Ghat[3]*dv11l; 
  outl[12] += -0.5*Ghat[12]*dv11l; 
  outl[13] += -0.5*Ghat[13]*dv11l; 
  outl[14] += -0.5*Ghat[14]*dv11l; 
  outl[15] += -0.8660254037844386*Ghat[5]*dv11l; 
  outl[16] += -0.5*Ghat[16]*dv11l; 
  outl[17] += -0.5*Ghat[17]*dv11l; 
  outl[18] += -0.5*Ghat[18]*dv11l; 
  outl[19] += -1.118033988749895*Ghat[0]*dv11l; 
  outl[20] += -0.5*Ghat[20]*dv11l; 
  outl[21] += -0.5*Ghat[21]*dv11l; 
  outl[22] += -0.8660254037844386*Ghat[6]*dv11l; 
  outl[23] += -0.8660254037844386*Ghat[7]*dv11l; 
  outl[24] += -0.8660254037844386*Ghat[8]*dv11l; 
  outl[25] += -0.5*Ghat[25]*dv11l; 
  outl[26] += -0.5*Ghat[26]*dv11l; 
  outl[27] += -0.5*Ghat[27]*dv11l; 
  outl[28] += -0.8660254037844386*Ghat[12]*dv11l; 
  outl[29] += -0.8660254037844386*Ghat[13]*dv11l; 
  outl[30] += -0.8660254037844386*Ghat[14]*dv11l; 
  outl[31] += -0.5*Ghat[31]*dv11l; 
  outl[32] += -0.5*Ghat[32]*dv11l; 
  outl[33] += -0.5*Ghat[33]*dv11l; 
  outl[34] += -0.5*Ghat[34]*dv11l; 
  outl[35] += -0.5*Ghat[35]*dv11l; 
  outl[36] += -0.5*Ghat[36]*dv11l; 
  outl[37] += -0.8660254037844387*Ghat[16]*dv11l; 
  outl[38] += -0.8660254037844387*Ghat[17]*dv11l; 
  outl[39] += -0.8660254037844387*Ghat[18]*dv11l; 
  outl[40] += -1.118033988749895*Ghat[1]*dv11l; 
  outl[41] += -1.118033988749895*Ghat[2]*dv11l; 
  outl[42] += -1.118033988749895*Ghat[3]*dv11l; 
  outl[43] += -0.5*Ghat[43]*dv11l; 
  outl[44] += -0.5*Ghat[44]*dv11l; 
  outl[45] += -0.5*Ghat[45]*dv11l; 
  outl[46] += -1.118033988749895*Ghat[5]*dv11l; 
  outl[47] += -0.5*Ghat[47]*dv11l; 
  outl[48] += -0.5*Ghat[48]*dv11l; 
  outl[49] += -0.5*Ghat[49]*dv11l; 
  outl[50] += -0.8660254037844387*Ghat[20]*dv11l; 
  outl[51] += -0.5*Ghat[51]*dv11l; 
  outl[52] += -0.5*Ghat[52]*dv11l; 
  outl[53] += -0.5*Ghat[53]*dv11l; 
  outl[54] += -1.322875655532295*Ghat[0]*dv11l; 
  outl[55] += -0.5*Ghat[55]*dv11l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12l = 2/dxvl[4]; 
  double dv12r = 2/dxvr[4]; 
  const double *E2 = &EM[6]; 

  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double Ghat[6]; 

  for(unsigned int i=0; i<6; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[6]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  double fjump[6]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(-1*fr[5]-fl[5]); 
  double alpha[6]; 

  alpha[0] = (-2.828427124746191*B0[0]*wv2)+2.828427124746191*B1[0]*wv1+2.828427124746191*E2[0]; 
  alpha[1] = (-2.828427124746191*B0[1]*wv2)+2.828427124746191*B1[1]*wv1+2.828427124746191*E2[1]; 
  alpha[2] = (-2.828427124746191*B0[2]*wv2)+2.828427124746191*B1[2]*wv1+2.828427124746191*E2[2]; 
  alpha[3] = 0.8164965809277261*B1[0]*dv1; 
  alpha[4] = -0.8164965809277261*B0[0]*dv2; 
  const double amid = 0.1767766952966368*alpha[0]; 
  Ghat[0] += (-0.8660254037844386*fjump[5])+0.1530931089239486*alpha[0]*favg[5]+0.0883883476483184*(alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.0883883476483184*alpha[0]*favg[0]; 
  Ghat[1] += 0.1530931089239486*alpha[1]*favg[5]-0.5*fjump[1]+0.0883883476483184*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += 0.1530931089239486*alpha[2]*favg[5]-0.5*fjump[2]+0.0883883476483184*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] += 0.1530931089239486*alpha[3]*favg[5]-0.5*fjump[3]+0.0883883476483184*(alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] += 0.1530931089239486*alpha[4]*favg[5]-0.5*fjump[4]+0.0883883476483184*(alpha[0]*favg[4]+favg[0]*alpha[4]); 

  outr[0] += 0.5*Ghat[0]*dv12r; 
  outr[1] += 0.5*Ghat[1]*dv12r; 
  outr[2] += 0.5*Ghat[2]*dv12r; 
  outr[3] += 0.5*Ghat[3]*dv12r; 
  outr[4] += 0.5*Ghat[4]*dv12r; 
  outr[5] += -0.8660254037844386*Ghat[0]*dv12r; 

  outl[0] += -0.5*Ghat[0]*dv12l; 
  outl[1] += -0.5*Ghat[1]*dv12l; 
  outl[2] += -0.5*Ghat[2]*dv12l; 
  outl[3] += -0.5*Ghat[3]*dv12l; 
  outl[4] += -0.5*Ghat[4]*dv12l; 
  outl[5] += -0.8660254037844386*Ghat[0]*dv12l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12l = 2/dxvl[4]; 
  double dv12r = 2/dxvr[4]; 
  const double *E2 = &EM[12]; 

  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double Ghat[21]; 

  for(unsigned int i=0; i<21; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[21]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = -1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  double fjump[21]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(-1*fr[5]-fl[5]); 
  fjump[6] = amax*(1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(1*fr[9]-fl[9]); 
  fjump[10] = amax*(1*fr[10]-fl[10]); 
  fjump[11] = amax*(1*fr[11]-fl[11]); 
  fjump[12] = amax*(-1*fr[12]-fl[12]); 
  fjump[13] = amax*(-1*fr[13]-fl[13]); 
  fjump[14] = amax*(-1*fr[14]-fl[14]); 
  fjump[15] = amax*(-1*fr[15]-fl[15]); 
  fjump[16] = amax*(1*fr[16]-fl[16]); 
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(1*fr[18]-fl[18]); 
  fjump[19] = amax*(1*fr[19]-fl[19]); 
  fjump[20] = amax*(1*fr[20]-fl[20]); 
  double alpha[21]; 

  alpha[0] = (-2.828427124746191*B0[0]*wv2)+2.828427124746191*B1[0]*wv1+2.828427124746191*E2[0]; 
  alpha[1] = (-2.828427124746191*B0[1]*wv2)+2.828427124746191*B1[1]*wv1+2.828427124746191*E2[1]; 
  alpha[2] = (-2.828427124746191*B0[2]*wv2)+2.828427124746191*B1[2]*wv1+2.828427124746191*E2[2]; 
  alpha[3] = 0.8164965809277261*B1[0]*dv1; 
  alpha[4] = -0.8164965809277261*B0[0]*dv2; 
  alpha[6] = (-2.828427124746191*B0[3]*wv2)+2.828427124746191*B1[3]*wv1+2.828427124746191*E2[3]; 
  alpha[7] = 0.8164965809277261*B1[1]*dv1; 
  alpha[8] = 0.8164965809277261*B1[2]*dv1; 
  alpha[9] = -0.8164965809277261*B0[1]*dv2; 
  alpha[10] = -0.8164965809277261*B0[2]*dv2; 
  alpha[16] = (-2.828427124746191*B0[4]*wv2)+2.828427124746191*B1[4]*wv1+2.828427124746191*E2[4]; 
  alpha[17] = (-2.828427124746191*B0[5]*wv2)+2.828427124746191*B1[5]*wv1+2.828427124746191*E2[5]; 
  const double amid = (-0.1976423537605236*alpha[17])-0.1976423537605236*alpha[16]+0.1767766952966368*alpha[0]; 
  Ghat[0] += (-1.118033988749895*fjump[20])+0.1976423537605236*alpha[0]*favg[20]+0.0883883476483184*(alpha[17]*favg[17]+alpha[16]*favg[16])+0.1530931089239486*(alpha[4]*favg[15]+alpha[3]*favg[14]+alpha[2]*favg[13]+alpha[1]*favg[12])+0.0883883476483184*(alpha[10]*favg[10]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6])-0.8660254037844386*fjump[5]+0.1530931089239486*alpha[0]*favg[5]+0.0883883476483184*(alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.0883883476483184*alpha[0]*favg[0]; 
  Ghat[1] += alpha[1]*(0.1976423537605236*favg[20]+0.07905694150420944*favg[16])+(0.1369306393762915*favg[12]+0.07905694150420944*favg[1])*alpha[16]+0.1530931089239486*(alpha[9]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13])-0.8660254037844386*fjump[12]+0.1530931089239486*alpha[0]*favg[12]+0.0883883476483184*(alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6])+0.1530931089239486*alpha[1]*favg[5]-0.5*fjump[1]+0.0883883476483184*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += alpha[2]*(0.1976423537605236*favg[20]+0.07905694150420944*favg[17])+(0.1369306393762915*favg[13]+0.07905694150420944*favg[2])*alpha[17]+0.1530931089239486*(alpha[10]*favg[15]+alpha[8]*favg[14])-0.8660254037844386*fjump[13]+0.1530931089239486*(alpha[0]*favg[13]+alpha[6]*favg[12])+0.0883883476483184*(alpha[4]*favg[10]+favg[4]*alpha[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[1]*favg[6]+favg[1]*alpha[6])+0.1530931089239486*alpha[2]*favg[5]-0.5*fjump[2]+0.0883883476483184*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] += alpha[3]*(0.1976423537605236*favg[20]+0.07905694150420944*favg[18])-0.8660254037844386*fjump[14]+0.1530931089239486*(alpha[0]*favg[14]+alpha[8]*favg[13]+alpha[7]*favg[12])+0.0883883476483184*(alpha[4]*favg[11]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[1]*favg[7]+favg[1]*alpha[7])+0.1530931089239486*alpha[3]*favg[5]-0.5*fjump[3]+0.0883883476483184*(alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] += alpha[4]*(0.1976423537605236*favg[20]+0.07905694150420944*favg[19])-0.8660254037844386*fjump[15]+0.1530931089239486*(alpha[0]*favg[15]+alpha[10]*favg[13]+alpha[9]*favg[12])+0.0883883476483184*(alpha[3]*favg[11]+alpha[2]*favg[10]+favg[2]*alpha[10]+alpha[1]*favg[9]+favg[1]*alpha[9])+0.1530931089239486*alpha[4]*favg[5]-0.5*fjump[4]+0.0883883476483184*(alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[6] += 0.1976423537605236*alpha[6]*favg[20]+0.07905694150420944*(alpha[6]*favg[17]+favg[6]*alpha[17]+alpha[6]*favg[16]+favg[6]*alpha[16])+0.1530931089239486*(alpha[1]*favg[13]+alpha[2]*favg[12])+0.0883883476483184*(alpha[9]*favg[10]+favg[9]*alpha[10]+alpha[7]*favg[8]+favg[7]*alpha[8])-0.5*fjump[6]+0.0883883476483184*alpha[0]*favg[6]+0.1530931089239486*favg[5]*alpha[6]+0.0883883476483184*(favg[0]*alpha[6]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[7] += 0.1976423537605236*alpha[7]*favg[20]+0.07905694150420944*(alpha[7]*(favg[18]+favg[16])+favg[7]*alpha[16])+0.1530931089239486*(alpha[1]*favg[14]+alpha[3]*favg[12])+0.0883883476483184*(alpha[9]*favg[11]+alpha[6]*favg[8]+favg[6]*alpha[8])-0.5*fjump[7]+0.0883883476483184*alpha[0]*favg[7]+0.1530931089239486*favg[5]*alpha[7]+0.0883883476483184*(favg[0]*alpha[7]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[8] += 0.1976423537605236*alpha[8]*favg[20]+0.07905694150420944*(alpha[8]*(favg[18]+favg[17])+favg[8]*alpha[17])+0.1530931089239486*(alpha[2]*favg[14]+alpha[3]*favg[13])+0.0883883476483184*alpha[10]*favg[11]-0.5*fjump[8]+0.0883883476483184*alpha[0]*favg[8]+0.1530931089239486*favg[5]*alpha[8]+0.0883883476483184*(favg[0]*alpha[8]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[9] += 0.1976423537605236*alpha[9]*favg[20]+0.07905694150420944*(alpha[9]*(favg[19]+favg[16])+favg[9]*alpha[16])+0.1530931089239486*(alpha[1]*favg[15]+alpha[4]*favg[12])+0.0883883476483184*(alpha[7]*favg[11]+alpha[6]*favg[10]+favg[6]*alpha[10])-0.5*fjump[9]+0.0883883476483184*alpha[0]*favg[9]+0.1530931089239486*favg[5]*alpha[9]+0.0883883476483184*(favg[0]*alpha[9]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[10] += 0.1976423537605236*alpha[10]*favg[20]+0.07905694150420944*(alpha[10]*(favg[19]+favg[17])+favg[10]*alpha[17])+0.1530931089239486*(alpha[2]*favg[15]+alpha[4]*favg[13])+0.0883883476483184*alpha[8]*favg[11]-0.5*fjump[10]+0.0883883476483184*alpha[0]*favg[10]+0.1530931089239486*favg[5]*alpha[10]+0.0883883476483184*(favg[0]*alpha[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[11] += 0.1530931089239486*(alpha[3]*favg[15]+alpha[4]*favg[14])-0.5*fjump[11]+0.0883883476483184*(alpha[0]*favg[11]+alpha[8]*favg[10]+favg[8]*alpha[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[16] += 0.1976423537605236*alpha[16]*favg[20]-0.5*fjump[16]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[16]+(0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[16]+0.1369306393762915*alpha[1]*favg[12]+0.07905694150420944*(alpha[9]*favg[9]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[1]*favg[1]); 
  Ghat[17] += 0.1976423537605236*alpha[17]*favg[20]-0.5*fjump[17]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[17]+(0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[17]+0.1369306393762915*alpha[2]*favg[13]+0.07905694150420944*(alpha[10]*favg[10]+alpha[8]*favg[8]+alpha[6]*favg[6]+alpha[2]*favg[2]); 
  Ghat[18] += (-0.5*fjump[18])+0.0883883476483184*alpha[0]*favg[18]+0.1369306393762915*alpha[3]*favg[14]+0.07905694150420944*(alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[3]*favg[3]); 
  Ghat[19] += (-0.5*fjump[19])+0.0883883476483184*alpha[0]*favg[19]+0.1369306393762915*alpha[4]*favg[15]+0.07905694150420944*(alpha[10]*favg[10]+alpha[9]*favg[9]+alpha[4]*favg[4]); 

  outr[0] += 0.5*Ghat[0]*dv12r; 
  outr[1] += 0.5*Ghat[1]*dv12r; 
  outr[2] += 0.5*Ghat[2]*dv12r; 
  outr[3] += 0.5*Ghat[3]*dv12r; 
  outr[4] += 0.5*Ghat[4]*dv12r; 
  outr[5] += -0.8660254037844386*Ghat[0]*dv12r; 
  outr[6] += 0.5*Ghat[6]*dv12r; 
  outr[7] += 0.5*Ghat[7]*dv12r; 
  outr[8] += 0.5*Ghat[8]*dv12r; 
  outr[9] += 0.5*Ghat[9]*dv12r; 
  outr[10] += 0.5*Ghat[10]*dv12r; 
  outr[11] += 0.5*Ghat[11]*dv12r; 
  outr[12] += -0.8660254037844386*Ghat[1]*dv12r; 
  outr[13] += -0.8660254037844386*Ghat[2]*dv12r; 
  outr[14] += -0.8660254037844386*Ghat[3]*dv12r; 
  outr[15] += -0.8660254037844386*Ghat[4]*dv12r; 
  outr[16] += 0.5*Ghat[16]*dv12r; 
  outr[17] += 0.5*Ghat[17]*dv12r; 
  outr[18] += 0.5*Ghat[18]*dv12r; 
  outr[19] += 0.5*Ghat[19]*dv12r; 
  outr[20] += 1.118033988749895*Ghat[0]*dv12r; 

  outl[0] += -0.5*Ghat[0]*dv12l; 
  outl[1] += -0.5*Ghat[1]*dv12l; 
  outl[2] += -0.5*Ghat[2]*dv12l; 
  outl[3] += -0.5*Ghat[3]*dv12l; 
  outl[4] += -0.5*Ghat[4]*dv12l; 
  outl[5] += -0.8660254037844386*Ghat[0]*dv12l; 
  outl[6] += -0.5*Ghat[6]*dv12l; 
  outl[7] += -0.5*Ghat[7]*dv12l; 
  outl[8] += -0.5*Ghat[8]*dv12l; 
  outl[9] += -0.5*Ghat[9]*dv12l; 
  outl[10] += -0.5*Ghat[10]*dv12l; 
  outl[11] += -0.5*Ghat[11]*dv12l; 
  outl[12] += -0.8660254037844386*Ghat[1]*dv12l; 
  outl[13] += -0.8660254037844386*Ghat[2]*dv12l; 
  outl[14] += -0.8660254037844386*Ghat[3]*dv12l; 
  outl[15] += -0.8660254037844386*Ghat[4]*dv12l; 
  outl[16] += -0.5*Ghat[16]*dv12l; 
  outl[17] += -0.5*Ghat[17]*dv12l; 
  outl[18] += -0.5*Ghat[18]*dv12l; 
  outl[19] += -0.5*Ghat[19]*dv12l; 
  outl[20] += -1.118033988749895*Ghat[0]*dv12l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12l = 2/dxvl[4]; 
  double dv12r = 2/dxvr[4]; 
  const double *E2 = &EM[20]; 

  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *B0 = &EM[30]; 
  const double *B1 = &EM[40]; 
  const double *B2 = &EM[50]; 

  double Ghat[56]; 

  for(unsigned int i=0; i<56; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[56]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = -1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  favg[21] = 1*fr[21]+fl[21]; 
  favg[22] = 1*fr[22]+fl[22]; 
  favg[23] = 1*fr[23]+fl[23]; 
  favg[24] = 1*fr[24]+fl[24]; 
  favg[25] = -1*fr[25]+fl[25]; 
  favg[26] = -1*fr[26]+fl[26]; 
  favg[27] = -1*fr[27]+fl[27]; 
  favg[28] = -1*fr[28]+fl[28]; 
  favg[29] = -1*fr[29]+fl[29]; 
  favg[30] = -1*fr[30]+fl[30]; 
  favg[31] = 1*fr[31]+fl[31]; 
  favg[32] = 1*fr[32]+fl[32]; 
  favg[33] = 1*fr[33]+fl[33]; 
  favg[34] = 1*fr[34]+fl[34]; 
  favg[35] = 1*fr[35]+fl[35]; 
  favg[36] = 1*fr[36]+fl[36]; 
  favg[37] = 1*fr[37]+fl[37]; 
  favg[38] = 1*fr[38]+fl[38]; 
  favg[39] = 1*fr[39]+fl[39]; 
  favg[40] = 1*fr[40]+fl[40]; 
  favg[41] = 1*fr[41]+fl[41]; 
  favg[42] = 1*fr[42]+fl[42]; 
  favg[43] = -1*fr[43]+fl[43]; 
  favg[44] = -1*fr[44]+fl[44]; 
  favg[45] = -1*fr[45]+fl[45]; 
  favg[46] = -1*fr[46]+fl[46]; 
  favg[47] = 1*fr[47]+fl[47]; 
  favg[48] = 1*fr[48]+fl[48]; 
  favg[49] = 1*fr[49]+fl[49]; 
  favg[50] = 1*fr[50]+fl[50]; 
  favg[51] = 1*fr[51]+fl[51]; 
  favg[52] = 1*fr[52]+fl[52]; 
  favg[53] = 1*fr[53]+fl[53]; 
  favg[54] = 1*fr[54]+fl[54]; 
  favg[55] = -1*fr[55]+fl[55]; 
  double fjump[56]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(1*fr[2]-fl[2]); 
  fjump[3] = amax*(1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(-1*fr[5]-fl[5]); 
  fjump[6] = amax*(1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(1*fr[9]-fl[9]); 
  fjump[10] = amax*(1*fr[10]-fl[10]); 
  fjump[11] = amax*(1*fr[11]-fl[11]); 
  fjump[12] = amax*(-1*fr[12]-fl[12]); 
  fjump[13] = amax*(-1*fr[13]-fl[13]); 
  fjump[14] = amax*(-1*fr[14]-fl[14]); 
  fjump[15] = amax*(-1*fr[15]-fl[15]); 
  fjump[16] = amax*(1*fr[16]-fl[16]); 
  fjump[17] = amax*(1*fr[17]-fl[17]); 
  fjump[18] = amax*(1*fr[18]-fl[18]); 
  fjump[19] = amax*(1*fr[19]-fl[19]); 
  fjump[20] = amax*(1*fr[20]-fl[20]); 
  fjump[21] = amax*(1*fr[21]-fl[21]); 
  fjump[22] = amax*(1*fr[22]-fl[22]); 
  fjump[23] = amax*(1*fr[23]-fl[23]); 
  fjump[24] = amax*(1*fr[24]-fl[24]); 
  fjump[25] = amax*(-1*fr[25]-fl[25]); 
  fjump[26] = amax*(-1*fr[26]-fl[26]); 
  fjump[27] = amax*(-1*fr[27]-fl[27]); 
  fjump[28] = amax*(-1*fr[28]-fl[28]); 
  fjump[29] = amax*(-1*fr[29]-fl[29]); 
  fjump[30] = amax*(-1*fr[30]-fl[30]); 
  fjump[31] = amax*(1*fr[31]-fl[31]); 
  fjump[32] = amax*(1*fr[32]-fl[32]); 
  fjump[33] = amax*(1*fr[33]-fl[33]); 
  fjump[34] = amax*(1*fr[34]-fl[34]); 
  fjump[35] = amax*(1*fr[35]-fl[35]); 
  fjump[36] = amax*(1*fr[36]-fl[36]); 
  fjump[37] = amax*(1*fr[37]-fl[37]); 
  fjump[38] = amax*(1*fr[38]-fl[38]); 
  fjump[39] = amax*(1*fr[39]-fl[39]); 
  fjump[40] = amax*(1*fr[40]-fl[40]); 
  fjump[41] = amax*(1*fr[41]-fl[41]); 
  fjump[42] = amax*(1*fr[42]-fl[42]); 
  fjump[43] = amax*(-1*fr[43]-fl[43]); 
  fjump[44] = amax*(-1*fr[44]-fl[44]); 
  fjump[45] = amax*(-1*fr[45]-fl[45]); 
  fjump[46] = amax*(-1*fr[46]-fl[46]); 
  fjump[47] = amax*(1*fr[47]-fl[47]); 
  fjump[48] = amax*(1*fr[48]-fl[48]); 
  fjump[49] = amax*(1*fr[49]-fl[49]); 
  fjump[50] = amax*(1*fr[50]-fl[50]); 
  fjump[51] = amax*(1*fr[51]-fl[51]); 
  fjump[52] = amax*(1*fr[52]-fl[52]); 
  fjump[53] = amax*(1*fr[53]-fl[53]); 
  fjump[54] = amax*(1*fr[54]-fl[54]); 
  fjump[55] = amax*(-1*fr[55]-fl[55]); 
  double alpha[56]; 

  alpha[0] = (-2.828427124746191*B0[0]*wv2)+2.828427124746191*B1[0]*wv1+2.828427124746191*E2[0]; 
  alpha[1] = (-2.828427124746191*B0[1]*wv2)+2.828427124746191*B1[1]*wv1+2.828427124746191*E2[1]; 
  alpha[2] = (-2.828427124746191*B0[2]*wv2)+2.828427124746191*B1[2]*wv1+2.828427124746191*E2[2]; 
  alpha[3] = 0.8164965809277261*B1[0]*dv1; 
  alpha[4] = -0.8164965809277261*B0[0]*dv2; 
  alpha[6] = (-2.828427124746191*B0[3]*wv2)+2.828427124746191*B1[3]*wv1+2.828427124746191*E2[3]; 
  alpha[7] = 0.8164965809277261*B1[1]*dv1; 
  alpha[8] = 0.8164965809277261*B1[2]*dv1; 
  alpha[9] = -0.8164965809277261*B0[1]*dv2; 
  alpha[10] = -0.8164965809277261*B0[2]*dv2; 
  alpha[16] = (-2.828427124746191*B0[4]*wv2)+2.828427124746191*B1[4]*wv1+2.828427124746191*E2[4]; 
  alpha[17] = (-2.828427124746191*B0[5]*wv2)+2.828427124746191*B1[5]*wv1+2.828427124746191*E2[5]; 
  alpha[21] = 0.8164965809277261*B1[3]*dv1; 
  alpha[22] = -0.8164965809277261*B0[3]*dv2; 
  alpha[31] = (-2.828427124746191*B0[6]*wv2)+2.828427124746191*B1[6]*wv1+2.828427124746191*E2[6]; 
  alpha[32] = (-2.828427124746191*B0[7]*wv2)+2.828427124746191*B1[7]*wv1+2.828427124746191*E2[7]; 
  alpha[33] = 0.816496580927726*B1[4]*dv1; 
  alpha[34] = 0.816496580927726*B1[5]*dv1; 
  alpha[37] = -0.816496580927726*B0[4]*dv2; 
  alpha[38] = -0.816496580927726*B0[5]*dv2; 
  alpha[51] = (-2.828427124746191*B0[8]*wv2)+2.828427124746191*B1[8]*wv1+2.828427124746191*E2[8]; 
  alpha[52] = (-2.828427124746191*B0[9]*wv2)+2.828427124746191*B1[9]*wv1+2.828427124746191*E2[9]; 
  const double amid = (-0.1976423537605236*alpha[17])-0.1976423537605236*alpha[16]+0.1767766952966368*alpha[0]; 
  Ghat[0] += (-1.322875655532295*fjump[55])+0.2338535866733712*alpha[0]*favg[55]+0.0883883476483184*(alpha[52]*favg[52]+alpha[51]*favg[51])+0.1976423537605236*(alpha[4]*favg[50]+alpha[3]*favg[49]+alpha[2]*favg[48]+alpha[1]*favg[47])+0.1530931089239486*(alpha[17]*favg[44]+alpha[16]*favg[43])+0.0883883476483184*(alpha[38]*favg[38]+alpha[37]*favg[37]+alpha[34]*favg[34]+alpha[33]*favg[33]+alpha[32]*favg[32]+alpha[31]*favg[31])+0.1530931089239486*(alpha[10]*favg[29]+alpha[9]*favg[28]+alpha[8]*favg[27]+alpha[7]*favg[26]+alpha[6]*favg[25])+0.0883883476483184*(alpha[22]*favg[22]+alpha[21]*favg[21])-1.118033988749895*fjump[20]+0.1976423537605236*alpha[0]*favg[20]+0.0883883476483184*(alpha[17]*favg[17]+alpha[16]*favg[16])+0.1530931089239486*(alpha[4]*favg[15]+alpha[3]*favg[14]+alpha[2]*favg[13]+alpha[1]*favg[12])+0.0883883476483184*(alpha[10]*favg[10]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6])-0.8660254037844386*fjump[5]+0.1530931089239486*alpha[0]*favg[5]+0.0883883476483184*(alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1])-0.5*fjump[0]+0.0883883476483184*alpha[0]*favg[0]; 
  Ghat[1] += 0.2338535866733712*alpha[1]*favg[55]+0.0776323754260148*alpha[16]*favg[51]+(0.1344632185501192*favg[43]+0.0776323754260148*favg[16])*alpha[51]+0.1976423537605236*(alpha[9]*favg[50]+alpha[7]*favg[49]+alpha[6]*favg[48])-1.118033988749895*fjump[47]+(0.1767766952966368*alpha[16]+0.1976423537605236*alpha[0])*favg[47]+0.1530931089239486*alpha[32]*favg[44]+0.1369306393762915*alpha[1]*favg[43]+0.07905694150420944*alpha[9]*favg[37]+0.1369306393762915*favg[28]*alpha[37]+0.07905694150420944*(favg[9]*alpha[37]+alpha[7]*favg[33])+(0.1369306393762915*favg[26]+0.07905694150420944*favg[7])*alpha[33]+0.0883883476483184*(alpha[17]*favg[32]+favg[17]*alpha[32])+0.07905694150420944*alpha[6]*favg[31]+(0.1369306393762915*favg[25]+0.07905694150420944*favg[6])*alpha[31]+0.1530931089239486*(alpha[22]*favg[29]+alpha[4]*favg[28]+alpha[21]*favg[27]+alpha[3]*favg[26]+alpha[2]*favg[25])+0.0883883476483184*(alpha[10]*favg[22]+favg[10]*alpha[22]+alpha[8]*favg[21]+favg[8]*alpha[21])+alpha[1]*(0.1976423537605236*favg[20]+0.07905694150420944*favg[16])+(0.1369306393762915*favg[12]+0.07905694150420944*favg[1])*alpha[16]+0.1530931089239486*(alpha[9]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13])-0.8660254037844386*fjump[12]+0.1530931089239486*alpha[0]*favg[12]+0.0883883476483184*(alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6])+0.1530931089239486*alpha[1]*favg[5]-0.5*fjump[1]+0.0883883476483184*(alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] += 0.2338535866733712*alpha[2]*favg[55]+0.0776323754260148*alpha[17]*favg[52]+(0.1344632185501192*favg[44]+0.0776323754260148*favg[17])*alpha[52]+0.1976423537605236*(alpha[10]*favg[50]+alpha[8]*favg[49])-1.118033988749895*fjump[48]+0.1767766952966368*alpha[17]*favg[48]+0.1976423537605236*(alpha[0]*favg[48]+alpha[6]*favg[47])+0.1369306393762915*alpha[2]*favg[44]+0.1530931089239486*alpha[31]*favg[43]+0.07905694150420944*alpha[10]*favg[38]+0.1369306393762915*favg[29]*alpha[38]+0.07905694150420944*(favg[10]*alpha[38]+alpha[8]*favg[34])+0.1369306393762915*favg[27]*alpha[34]+0.07905694150420944*(favg[8]*alpha[34]+alpha[6]*favg[32])+(0.1369306393762915*favg[25]+0.07905694150420944*favg[6])*alpha[32]+0.0883883476483184*(alpha[16]*favg[31]+favg[16]*alpha[31])+0.1530931089239486*(alpha[4]*favg[29]+alpha[22]*favg[28]+alpha[3]*favg[27]+alpha[21]*favg[26]+alpha[1]*favg[25])+0.0883883476483184*(alpha[9]*favg[22]+favg[9]*alpha[22]+alpha[7]*favg[21]+favg[7]*alpha[21])+alpha[2]*(0.1976423537605236*favg[20]+0.07905694150420944*favg[17])+(0.1369306393762915*favg[13]+0.07905694150420944*favg[2])*alpha[17]+0.1530931089239486*(alpha[10]*favg[15]+alpha[8]*favg[14])-0.8660254037844386*fjump[13]+0.1530931089239486*(alpha[0]*favg[13]+alpha[6]*favg[12])+0.0883883476483184*(alpha[4]*favg[10]+favg[4]*alpha[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[1]*favg[6]+favg[1]*alpha[6])+0.1530931089239486*alpha[2]*favg[5]-0.5*fjump[2]+0.0883883476483184*(alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] += 0.2338535866733712*alpha[3]*favg[55]-1.118033988749895*fjump[49]+0.1976423537605236*(alpha[0]*favg[49]+alpha[8]*favg[48]+alpha[7]*favg[47])+0.1369306393762915*alpha[3]*favg[45]+0.1530931089239486*(alpha[34]*favg[44]+alpha[33]*favg[43])+0.07905694150420944*(alpha[8]*favg[36]+alpha[7]*favg[35])+0.0883883476483184*(alpha[17]*favg[34]+favg[17]*alpha[34]+alpha[16]*favg[33]+favg[16]*alpha[33])+0.1530931089239486*(alpha[4]*favg[30]+alpha[2]*favg[27]+alpha[1]*favg[26]+alpha[21]*favg[25])+0.0883883476483184*(alpha[10]*favg[24]+alpha[9]*favg[23]+alpha[6]*favg[21]+favg[6]*alpha[21])+alpha[3]*(0.1976423537605236*favg[20]+0.07905694150420944*favg[18])-0.8660254037844386*fjump[14]+0.1530931089239486*(alpha[0]*favg[14]+alpha[8]*favg[13]+alpha[7]*favg[12])+0.0883883476483184*(alpha[4]*favg[11]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[1]*favg[7]+favg[1]*alpha[7])+0.1530931089239486*alpha[3]*favg[5]-0.5*fjump[3]+0.0883883476483184*(alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] += 0.2338535866733712*alpha[4]*favg[55]-1.118033988749895*fjump[50]+0.1976423537605236*(alpha[0]*favg[50]+alpha[10]*favg[48]+alpha[9]*favg[47])+0.1369306393762915*alpha[4]*favg[46]+0.1530931089239486*(alpha[38]*favg[44]+alpha[37]*favg[43])+0.07905694150420944*(alpha[10]*favg[41]+alpha[9]*favg[40])+0.0883883476483184*(alpha[17]*favg[38]+favg[17]*alpha[38]+alpha[16]*favg[37]+favg[16]*alpha[37])+0.1530931089239486*(alpha[3]*favg[30]+alpha[2]*favg[29]+alpha[1]*favg[28]+alpha[22]*favg[25])+0.0883883476483184*(alpha[8]*favg[24]+alpha[7]*favg[23]+alpha[6]*favg[22]+favg[6]*alpha[22])+alpha[4]*(0.1976423537605236*favg[20]+0.07905694150420944*favg[19])-0.8660254037844386*fjump[15]+0.1530931089239486*(alpha[0]*favg[15]+alpha[10]*favg[13]+alpha[9]*favg[12])+0.0883883476483184*(alpha[3]*favg[11]+alpha[2]*favg[10]+favg[2]*alpha[10]+alpha[1]*favg[9]+favg[1]*alpha[9])+0.1530931089239486*alpha[4]*favg[5]-0.5*fjump[4]+0.0883883476483184*(alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[6] += 0.2338535866733712*alpha[6]*favg[55]+0.0776323754260148*(alpha[32]*favg[52]+favg[32]*alpha[52]+alpha[31]*favg[51]+favg[31]*alpha[51])+0.1976423537605236*(alpha[22]*favg[50]+alpha[21]*favg[49])+(0.1767766952966368*alpha[32]+0.1976423537605236*alpha[1])*favg[48]+(0.1767766952966368*alpha[31]+0.1976423537605236*alpha[2])*favg[47]+0.1369306393762915*alpha[6]*(favg[44]+favg[43])+0.07905694150420944*(alpha[22]*favg[38]+favg[22]*alpha[38]+alpha[22]*favg[37]+favg[22]*alpha[37]+alpha[21]*favg[34]+favg[21]*alpha[34]+alpha[21]*favg[33]+favg[21]*alpha[33])+(0.07071067811865474*alpha[31]+0.07905694150420944*alpha[2])*favg[32]+(0.07071067811865474*favg[31]+0.1369306393762915*favg[13])*alpha[32]+0.07905694150420944*(favg[2]*alpha[32]+alpha[1]*favg[31])+(0.1369306393762915*favg[12]+0.07905694150420944*favg[1])*alpha[31]+0.1530931089239486*(alpha[9]*favg[29]+alpha[10]*favg[28]+alpha[7]*favg[27]+alpha[8]*favg[26])-0.8660254037844386*fjump[25]+(0.1369306393762915*(alpha[17]+alpha[16])+0.1530931089239486*alpha[0])*favg[25]+0.0883883476483184*alpha[4]*favg[22]+0.1530931089239486*favg[15]*alpha[22]+0.0883883476483184*(favg[4]*alpha[22]+alpha[3]*favg[21])+(0.1530931089239486*favg[14]+0.0883883476483184*favg[3])*alpha[21]+0.1976423537605236*alpha[6]*favg[20]+0.07905694150420944*(alpha[6]*favg[17]+favg[6]*alpha[17]+alpha[6]*favg[16]+favg[6]*alpha[16])+0.1530931089239486*(alpha[1]*favg[13]+alpha[2]*favg[12])+0.0883883476483184*(alpha[9]*favg[10]+favg[9]*alpha[10]+alpha[7]*favg[8]+favg[7]*alpha[8])-0.5*fjump[6]+0.0883883476483184*alpha[0]*favg[6]+0.1530931089239486*favg[5]*alpha[6]+0.0883883476483184*(favg[0]*alpha[6]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[7] += 0.2338535866733712*alpha[7]*favg[55]+0.0776323754260148*(alpha[33]*favg[51]+favg[33]*alpha[51])+0.1976423537605236*(alpha[1]*favg[49]+alpha[21]*favg[48])+(0.1767766952966368*alpha[33]+0.1976423537605236*alpha[3])*favg[47]+0.1369306393762915*alpha[7]*(favg[45]+favg[43])+0.07905694150420944*(favg[23]*alpha[37]+alpha[21]*favg[36])+(0.07071067811865474*alpha[33]+0.07905694150420944*alpha[3])*favg[35]+0.0883883476483184*(alpha[32]*favg[34]+favg[32]*alpha[34])+0.07905694150420944*alpha[1]*favg[33]+0.1369306393762915*favg[12]*alpha[33]+0.07905694150420944*(favg[1]*alpha[33]+alpha[21]*favg[31]+favg[21]*alpha[31])+0.1530931089239486*(alpha[9]*favg[30]+alpha[6]*favg[27])-0.8660254037844386*fjump[26]+0.1369306393762915*alpha[16]*favg[26]+0.1530931089239486*(alpha[0]*favg[26]+alpha[8]*favg[25])+0.0883883476483184*(alpha[22]*favg[24]+alpha[4]*favg[23]+alpha[2]*favg[21])+(0.1530931089239486*favg[13]+0.0883883476483184*favg[2])*alpha[21]+0.1976423537605236*alpha[7]*favg[20]+0.07905694150420944*(alpha[7]*(favg[18]+favg[16])+favg[7]*alpha[16])+0.1530931089239486*(alpha[1]*favg[14]+alpha[3]*favg[12])+0.0883883476483184*(alpha[9]*favg[11]+alpha[6]*favg[8]+favg[6]*alpha[8])-0.5*fjump[7]+0.0883883476483184*alpha[0]*favg[7]+0.1530931089239486*favg[5]*alpha[7]+0.0883883476483184*(favg[0]*alpha[7]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[8] += 0.2338535866733712*alpha[8]*favg[55]+0.0776323754260148*(alpha[34]*favg[52]+favg[34]*alpha[52])+0.1976423537605236*alpha[2]*favg[49]+0.1767766952966368*alpha[34]*favg[48]+0.1976423537605236*(alpha[3]*favg[48]+alpha[21]*favg[47])+0.1369306393762915*alpha[8]*(favg[45]+favg[44])+0.07905694150420944*favg[24]*alpha[38]+0.07071067811865474*alpha[34]*favg[36]+0.07905694150420944*(alpha[3]*favg[36]+alpha[21]*favg[35]+alpha[2]*favg[34])+(0.1369306393762915*favg[13]+0.07905694150420944*favg[2])*alpha[34]+0.0883883476483184*(alpha[31]*favg[33]+favg[31]*alpha[33])+0.07905694150420944*(alpha[21]*favg[32]+favg[21]*alpha[32])+0.1530931089239486*alpha[10]*favg[30]-0.8660254037844386*fjump[27]+0.1369306393762915*alpha[17]*favg[27]+0.1530931089239486*(alpha[0]*favg[27]+alpha[6]*favg[26]+alpha[7]*favg[25])+0.0883883476483184*(alpha[4]*favg[24]+alpha[22]*favg[23]+alpha[1]*favg[21])+(0.1530931089239486*favg[12]+0.0883883476483184*favg[1])*alpha[21]+0.1976423537605236*alpha[8]*favg[20]+0.07905694150420944*(alpha[8]*(favg[18]+favg[17])+favg[8]*alpha[17])+0.1530931089239486*(alpha[2]*favg[14]+alpha[3]*favg[13])+0.0883883476483184*alpha[10]*favg[11]-0.5*fjump[8]+0.0883883476483184*alpha[0]*favg[8]+0.1530931089239486*favg[5]*alpha[8]+0.0883883476483184*(favg[0]*alpha[8]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[9] += 0.2338535866733712*alpha[9]*favg[55]+0.0776323754260148*(alpha[37]*favg[51]+favg[37]*alpha[51])+0.1976423537605236*(alpha[1]*favg[50]+alpha[22]*favg[48])+(0.1767766952966368*alpha[37]+0.1976423537605236*alpha[4])*favg[47]+0.1369306393762915*alpha[9]*(favg[46]+favg[43])+0.07905694150420944*alpha[22]*favg[41]+(0.07071067811865474*alpha[37]+0.07905694150420944*alpha[4])*favg[40]+0.0883883476483184*(alpha[32]*favg[38]+favg[32]*alpha[38])+0.07905694150420944*alpha[1]*favg[37]+0.1369306393762915*favg[12]*alpha[37]+0.07905694150420944*(favg[1]*alpha[37]+favg[23]*alpha[33]+alpha[22]*favg[31]+favg[22]*alpha[31])+0.1530931089239486*(alpha[7]*favg[30]+alpha[6]*favg[29])-0.8660254037844386*fjump[28]+0.1369306393762915*alpha[16]*favg[28]+0.1530931089239486*(alpha[0]*favg[28]+alpha[10]*favg[25])+0.0883883476483184*(alpha[21]*favg[24]+alpha[3]*favg[23]+alpha[2]*favg[22])+(0.1530931089239486*favg[13]+0.0883883476483184*favg[2])*alpha[22]+0.1976423537605236*alpha[9]*favg[20]+0.07905694150420944*(alpha[9]*(favg[19]+favg[16])+favg[9]*alpha[16])+0.1530931089239486*(alpha[1]*favg[15]+alpha[4]*favg[12])+0.0883883476483184*(alpha[7]*favg[11]+alpha[6]*favg[10]+favg[6]*alpha[10])-0.5*fjump[9]+0.0883883476483184*alpha[0]*favg[9]+0.1530931089239486*favg[5]*alpha[9]+0.0883883476483184*(favg[0]*alpha[9]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[10] += 0.2338535866733712*alpha[10]*favg[55]+0.0776323754260148*(alpha[38]*favg[52]+favg[38]*alpha[52])+0.1976423537605236*alpha[2]*favg[50]+0.1767766952966368*alpha[38]*favg[48]+0.1976423537605236*(alpha[4]*favg[48]+alpha[22]*favg[47])+0.1369306393762915*alpha[10]*(favg[46]+favg[44])+0.07071067811865474*alpha[38]*favg[41]+0.07905694150420944*(alpha[4]*favg[41]+alpha[22]*favg[40]+alpha[2]*favg[38])+(0.1369306393762915*favg[13]+0.07905694150420944*favg[2])*alpha[38]+0.0883883476483184*(alpha[31]*favg[37]+favg[31]*alpha[37])+0.07905694150420944*(favg[24]*alpha[34]+alpha[22]*favg[32]+favg[22]*alpha[32])+0.1530931089239486*alpha[8]*favg[30]-0.8660254037844386*fjump[29]+0.1369306393762915*alpha[17]*favg[29]+0.1530931089239486*(alpha[0]*favg[29]+alpha[6]*favg[28]+alpha[9]*favg[25])+0.0883883476483184*(alpha[3]*favg[24]+alpha[21]*favg[23]+alpha[1]*favg[22])+(0.1530931089239486*favg[12]+0.0883883476483184*favg[1])*alpha[22]+0.1976423537605236*alpha[10]*favg[20]+0.07905694150420944*(alpha[10]*(favg[19]+favg[17])+favg[10]*alpha[17])+0.1530931089239486*(alpha[2]*favg[15]+alpha[4]*favg[13])+0.0883883476483184*alpha[8]*favg[11]-0.5*fjump[10]+0.0883883476483184*alpha[0]*favg[10]+0.1530931089239486*favg[5]*alpha[10]+0.0883883476483184*(favg[0]*alpha[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[11] += 0.1976423537605236*(alpha[3]*favg[50]+alpha[4]*favg[49])+0.07905694150420944*(alpha[4]*favg[42]+alpha[3]*favg[39])+0.0883883476483184*(alpha[34]*favg[38]+favg[34]*alpha[38]+alpha[33]*favg[37]+favg[33]*alpha[37])-0.8660254037844386*fjump[30]+0.1530931089239486*(alpha[0]*favg[30]+alpha[8]*favg[29]+alpha[7]*favg[28]+alpha[10]*favg[27]+alpha[9]*favg[26])+0.0883883476483184*(alpha[2]*favg[24]+alpha[1]*favg[23]+alpha[21]*favg[22]+favg[21]*alpha[22])+0.1530931089239486*(alpha[3]*favg[15]+alpha[4]*favg[14])-0.5*fjump[11]+0.0883883476483184*(alpha[0]*favg[11]+alpha[8]*favg[10]+favg[8]*alpha[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[16] += 0.2338535866733712*alpha[16]*favg[55]+(0.05270462766947297*alpha[51]+0.0776323754260148*alpha[1])*favg[51]+(0.1735912687073533*favg[47]+0.1344632185501192*favg[12]+0.0776323754260148*favg[1])*alpha[51]+0.1976423537605236*(alpha[37]*favg[50]+alpha[33]*favg[49]+alpha[31]*favg[48])+0.1767766952966368*alpha[1]*favg[47]-0.8660254037844386*fjump[43]+(0.09780759955449389*alpha[16]+0.1530931089239486*alpha[0])*favg[43]+(0.05646924393157818*alpha[37]+0.0883883476483184*alpha[4])*favg[37]+(0.1530931089239486*favg[15]+0.0883883476483184*favg[4])*alpha[37]+(0.05646924393157818*alpha[33]+0.0883883476483184*alpha[3])*favg[33]+(0.1530931089239486*favg[14]+0.0883883476483184*favg[3])*alpha[33]+0.07905694150420944*alpha[32]*favg[32]+(0.05646924393157818*alpha[31]+0.0883883476483184*alpha[2])*favg[31]+(0.1530931089239486*favg[13]+0.0883883476483184*favg[2])*alpha[31]+0.1369306393762915*(alpha[9]*favg[28]+alpha[7]*favg[26]+alpha[6]*favg[25])+0.07905694150420944*(alpha[22]*favg[22]+alpha[21]*favg[21])+0.1976423537605236*alpha[16]*favg[20]-0.5*fjump[16]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[16]+(0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[16]+0.1369306393762915*alpha[1]*favg[12]+0.07905694150420944*(alpha[9]*favg[9]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[1]*favg[1]); 
  Ghat[17] += 0.2338535866733712*alpha[17]*favg[55]+(0.05270462766947297*alpha[52]+0.0776323754260148*alpha[2])*favg[52]+(0.1735912687073533*favg[48]+0.1344632185501192*favg[13]+0.0776323754260148*favg[2])*alpha[52]+0.1976423537605236*(alpha[38]*favg[50]+alpha[34]*favg[49])+0.1767766952966368*alpha[2]*favg[48]+0.1976423537605236*alpha[32]*favg[47]-0.8660254037844386*fjump[44]+(0.09780759955449389*alpha[17]+0.1530931089239486*alpha[0])*favg[44]+(0.05646924393157818*alpha[38]+0.0883883476483184*alpha[4])*favg[38]+(0.1530931089239486*favg[15]+0.0883883476483184*favg[4])*alpha[38]+(0.05646924393157818*alpha[34]+0.0883883476483184*alpha[3])*favg[34]+(0.1530931089239486*favg[14]+0.0883883476483184*favg[3])*alpha[34]+(0.05646924393157818*alpha[32]+0.0883883476483184*alpha[1])*favg[32]+(0.1530931089239486*favg[12]+0.0883883476483184*favg[1])*alpha[32]+0.07905694150420944*alpha[31]*favg[31]+0.1369306393762915*(alpha[10]*favg[29]+alpha[8]*favg[27]+alpha[6]*favg[25])+0.07905694150420944*(alpha[22]*favg[22]+alpha[21]*favg[21])+0.1976423537605236*alpha[17]*favg[20]-0.5*fjump[17]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[17]+(0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[17]+0.1369306393762915*alpha[2]*favg[13]+0.07905694150420944*(alpha[10]*favg[10]+alpha[8]*favg[8]+alpha[6]*favg[6]+alpha[2]*favg[2]); 
  Ghat[18] += alpha[3]*(0.0776323754260148*favg[53]+0.1767766952966368*favg[49])-0.8660254037844386*fjump[45]+0.1530931089239486*alpha[0]*favg[45]+0.0883883476483184*(alpha[4]*favg[39]+alpha[2]*favg[36]+alpha[1]*favg[35])+0.07905694150420944*(alpha[34]*favg[34]+alpha[33]*favg[33])+0.1369306393762915*(alpha[8]*favg[27]+alpha[7]*favg[26])+0.07905694150420944*alpha[21]*favg[21]-0.5*fjump[18]+0.0883883476483184*alpha[0]*favg[18]+0.1369306393762915*alpha[3]*favg[14]+0.07905694150420944*(alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[3]*favg[3]); 
  Ghat[19] += alpha[4]*(0.0776323754260148*favg[54]+0.1767766952966368*favg[50])-0.8660254037844386*fjump[46]+0.1530931089239486*alpha[0]*favg[46]+0.0883883476483184*(alpha[3]*favg[42]+alpha[2]*favg[41]+alpha[1]*favg[40])+0.07905694150420944*(alpha[38]*favg[38]+alpha[37]*favg[37])+0.1369306393762915*(alpha[10]*favg[29]+alpha[9]*favg[28])+0.07905694150420944*alpha[22]*favg[22]-0.5*fjump[19]+0.0883883476483184*alpha[0]*favg[19]+0.1369306393762915*alpha[4]*favg[15]+0.07905694150420944*(alpha[10]*favg[10]+alpha[9]*favg[9]+alpha[4]*favg[4]); 
  Ghat[21] += 0.2338535866733712*alpha[21]*favg[55]+0.1976423537605236*(alpha[6]*favg[49]+alpha[7]*favg[48]+alpha[8]*favg[47])+0.1369306393762915*alpha[21]*(favg[45]+favg[44]+favg[43])+0.07905694150420944*(alpha[7]*favg[36]+alpha[8]*favg[35]+alpha[6]*favg[34])+0.1369306393762915*favg[25]*alpha[34]+0.07905694150420944*(favg[6]*alpha[34]+alpha[6]*favg[33])+0.1369306393762915*favg[25]*alpha[33]+0.07905694150420944*(favg[6]*alpha[33]+alpha[8]*favg[32])+0.1369306393762915*favg[27]*alpha[32]+0.07905694150420944*(favg[8]*alpha[32]+alpha[7]*favg[31])+(0.1369306393762915*favg[26]+0.07905694150420944*favg[7])*alpha[31]+0.1530931089239486*(alpha[22]*favg[30]+alpha[1]*favg[27]+alpha[2]*favg[26]+alpha[3]*favg[25])+0.0883883476483184*(alpha[9]*favg[24]+alpha[10]*favg[23]+favg[11]*alpha[22])-0.5*fjump[21]+(0.07905694150420944*(alpha[17]+alpha[16])+0.0883883476483184*alpha[0])*favg[21]+(0.1976423537605236*favg[20]+0.07905694150420944*(favg[18]+favg[17]+favg[16])+0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[21]+0.1530931089239486*(alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[8]*favg[12])+0.0883883476483184*(alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[3]*favg[6]+favg[3]*alpha[6]); 
  Ghat[22] += 0.2338535866733712*alpha[22]*favg[55]+0.1976423537605236*(alpha[6]*favg[50]+alpha[9]*favg[48]+alpha[10]*favg[47])+0.1369306393762915*alpha[22]*(favg[46]+favg[44]+favg[43])+0.07905694150420944*(alpha[9]*favg[41]+alpha[10]*favg[40]+alpha[6]*favg[38])+0.1369306393762915*favg[25]*alpha[38]+0.07905694150420944*(favg[6]*alpha[38]+alpha[6]*favg[37])+0.1369306393762915*favg[25]*alpha[37]+0.07905694150420944*(favg[6]*alpha[37]+alpha[10]*favg[32])+0.1369306393762915*favg[29]*alpha[32]+0.07905694150420944*(favg[10]*alpha[32]+alpha[9]*favg[31])+(0.1369306393762915*favg[28]+0.07905694150420944*favg[9])*alpha[31]+0.1530931089239486*(alpha[21]*favg[30]+alpha[1]*favg[29]+alpha[2]*favg[28]+alpha[4]*favg[25])+0.0883883476483184*(alpha[7]*favg[24]+alpha[8]*favg[23])-0.5*fjump[22]+(0.07905694150420944*(alpha[17]+alpha[16])+0.0883883476483184*alpha[0])*favg[22]+(0.1976423537605236*favg[20]+0.07905694150420944*(favg[19]+favg[17]+favg[16])+0.1530931089239486*favg[5])*alpha[22]+0.0883883476483184*(favg[0]*alpha[22]+favg[11]*alpha[21])+0.1530931089239486*(alpha[6]*favg[15]+alpha[9]*favg[13]+alpha[10]*favg[12])+0.0883883476483184*(alpha[1]*favg[10]+favg[1]*alpha[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[23] += 0.1976423537605236*(alpha[7]*favg[50]+alpha[9]*favg[49])+0.07905694150420944*(alpha[9]*favg[42]+alpha[7]*(favg[39]+favg[37]))+0.1369306393762915*favg[26]*alpha[37]+0.07905694150420944*(favg[7]*alpha[37]+alpha[9]*favg[33])+(0.1369306393762915*favg[28]+0.07905694150420944*favg[9])*alpha[33]+0.1530931089239486*(alpha[1]*favg[30]+alpha[21]*favg[29]+alpha[3]*favg[28]+alpha[22]*favg[27]+alpha[4]*favg[26])+0.0883883476483184*alpha[6]*favg[24]-0.5*fjump[23]+0.07905694150420944*alpha[16]*favg[23]+0.0883883476483184*(alpha[0]*favg[23]+alpha[8]*favg[22]+favg[8]*alpha[22]+alpha[10]*favg[21]+favg[10]*alpha[21])+0.1530931089239486*(alpha[7]*favg[15]+alpha[9]*favg[14])+0.0883883476483184*(alpha[1]*favg[11]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[24] += 0.1976423537605236*(alpha[8]*favg[50]+alpha[10]*favg[49])+0.07905694150420944*(alpha[10]*favg[42]+alpha[8]*(favg[39]+favg[38]))+0.1369306393762915*favg[27]*alpha[38]+0.07905694150420944*(favg[8]*alpha[38]+alpha[10]*favg[34])+(0.1369306393762915*favg[29]+0.07905694150420944*favg[10])*alpha[34]+0.1530931089239486*(alpha[2]*favg[30]+alpha[3]*favg[29]+alpha[21]*favg[28]+alpha[4]*favg[27]+alpha[22]*favg[26])-0.5*fjump[24]+0.07905694150420944*alpha[17]*favg[24]+0.0883883476483184*(alpha[0]*favg[24]+alpha[6]*favg[23]+alpha[7]*favg[22]+favg[7]*alpha[22]+alpha[9]*favg[21]+favg[9]*alpha[21])+0.1530931089239486*(alpha[8]*favg[15]+alpha[10]*favg[14])+0.0883883476483184*(alpha[2]*favg[11]+alpha[3]*favg[10]+favg[3]*alpha[10]+alpha[4]*favg[8]+favg[4]*alpha[8]); 
  Ghat[31] += 0.2338535866733712*alpha[31]*favg[55]+0.0776323754260148*alpha[6]*favg[51]+(0.1344632185501192*favg[25]+0.0776323754260148*favg[6])*alpha[51]+0.1976423537605236*alpha[16]*favg[48]+0.1767766952966368*alpha[6]*favg[47]+0.1369306393762915*alpha[31]*favg[44]+(0.09780759955449389*alpha[31]+0.1530931089239486*alpha[2])*favg[43]+0.0883883476483184*alpha[10]*favg[37]+0.1530931089239486*favg[29]*alpha[37]+0.0883883476483184*(favg[10]*alpha[37]+alpha[8]*favg[33])+(0.1530931089239486*favg[27]+0.0883883476483184*favg[8])*alpha[33]+0.07071067811865474*alpha[6]*favg[32]+(0.1224744871391589*favg[25]+0.07071067811865474*favg[6])*alpha[32]-0.5*fjump[31]+(0.07905694150420944*alpha[17]+0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[31]+(0.1976423537605236*favg[20]+0.07905694150420944*favg[17]+0.05646924393157818*favg[16]+0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[31]+0.1369306393762915*(alpha[22]*favg[28]+alpha[21]*favg[26]+alpha[1]*favg[25])+0.07905694150420944*(alpha[9]*favg[22]+favg[9]*alpha[22]+alpha[7]*favg[21]+favg[7]*alpha[21])+0.0883883476483184*alpha[2]*favg[16]+(0.1530931089239486*favg[13]+0.0883883476483184*favg[2])*alpha[16]+0.1369306393762915*alpha[6]*favg[12]+0.07905694150420944*(alpha[1]*favg[6]+favg[1]*alpha[6]); 
  Ghat[32] += 0.2338535866733712*alpha[32]*favg[55]+0.0776323754260148*alpha[6]*favg[52]+(0.1344632185501192*favg[25]+0.0776323754260148*favg[6])*alpha[52]+0.1767766952966368*alpha[6]*favg[48]+0.1976423537605236*alpha[17]*favg[47]+(0.09780759955449389*alpha[32]+0.1530931089239486*alpha[1])*favg[44]+0.1369306393762915*alpha[32]*favg[43]+0.0883883476483184*alpha[9]*favg[38]+0.1530931089239486*favg[28]*alpha[38]+0.0883883476483184*(favg[9]*alpha[38]+alpha[7]*favg[34])+(0.1530931089239486*favg[26]+0.0883883476483184*favg[7])*alpha[34]-0.5*fjump[32]+(0.05646924393157818*alpha[17]+0.07905694150420944*alpha[16]+0.0883883476483184*alpha[0])*favg[32]+(0.1976423537605236*favg[20]+0.05646924393157818*favg[17]+0.07905694150420944*favg[16]+0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[32]+0.07071067811865474*alpha[6]*favg[31]+(0.1224744871391589*favg[25]+0.07071067811865474*favg[6])*alpha[31]+0.1369306393762915*(alpha[22]*favg[29]+alpha[21]*favg[27]+alpha[2]*favg[25])+0.07905694150420944*(alpha[10]*favg[22]+favg[10]*alpha[22]+alpha[8]*favg[21]+favg[8]*alpha[21])+0.0883883476483184*alpha[1]*favg[17]+(0.1530931089239486*favg[12]+0.0883883476483184*favg[1])*alpha[17]+0.1369306393762915*alpha[6]*favg[13]+0.07905694150420944*(alpha[2]*favg[6]+favg[2]*alpha[6]); 
  Ghat[33] += 0.2338535866733712*alpha[33]*favg[55]+0.0776323754260148*alpha[7]*favg[51]+(0.1344632185501192*favg[26]+0.0776323754260148*favg[7])*alpha[51]+0.1976423537605236*alpha[16]*favg[49]+0.1767766952966368*alpha[7]*favg[47]+0.1369306393762915*alpha[33]*favg[45]+(0.09780759955449389*alpha[33]+0.1530931089239486*alpha[3])*favg[43]+(0.1530931089239486*favg[30]+0.0883883476483184*favg[11])*alpha[37]+0.07071067811865474*alpha[7]*favg[35]-0.5*fjump[33]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[33]+(0.1976423537605236*favg[20]+0.07905694150420944*favg[18]+0.05646924393157818*favg[16]+0.1530931089239486*favg[5])*alpha[33]+0.0883883476483184*(favg[0]*alpha[33]+alpha[8]*favg[31])+(0.1530931089239486*favg[27]+0.0883883476483184*favg[8])*alpha[31]+0.1369306393762915*(alpha[1]*favg[26]+alpha[21]*favg[25])+0.07905694150420944*(alpha[9]*favg[23]+alpha[6]*favg[21]+favg[6]*alpha[21])+0.0883883476483184*alpha[3]*favg[16]+(0.1530931089239486*favg[14]+0.0883883476483184*favg[3])*alpha[16]+0.1369306393762915*alpha[7]*favg[12]+0.07905694150420944*(alpha[1]*favg[7]+favg[1]*alpha[7]); 
  Ghat[34] += 0.2338535866733712*alpha[34]*favg[55]+0.0776323754260148*alpha[8]*favg[52]+(0.1344632185501192*favg[27]+0.0776323754260148*favg[8])*alpha[52]+0.1976423537605236*alpha[17]*favg[49]+0.1767766952966368*alpha[8]*favg[48]+0.1369306393762915*alpha[34]*favg[45]+(0.09780759955449389*alpha[34]+0.1530931089239486*alpha[3])*favg[44]+(0.1530931089239486*favg[30]+0.0883883476483184*favg[11])*alpha[38]+0.07071067811865474*alpha[8]*favg[36]-0.5*fjump[34]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[34]+(0.1976423537605236*favg[20]+0.07905694150420944*favg[18]+0.05646924393157818*favg[17]+0.1530931089239486*favg[5])*alpha[34]+0.0883883476483184*(favg[0]*alpha[34]+alpha[7]*favg[32])+(0.1530931089239486*favg[26]+0.0883883476483184*favg[7])*alpha[32]+0.1369306393762915*(alpha[2]*favg[27]+alpha[21]*favg[25])+0.07905694150420944*(alpha[10]*favg[24]+alpha[6]*favg[21]+favg[6]*alpha[21])+0.0883883476483184*alpha[3]*favg[17]+(0.1530931089239486*favg[14]+0.0883883476483184*favg[3])*alpha[17]+0.1369306393762915*alpha[8]*favg[13]+0.07905694150420944*(alpha[2]*favg[8]+favg[2]*alpha[8]); 
  Ghat[35] += alpha[7]*(0.0776323754260148*favg[53]+0.1767766952966368*favg[49])+0.1530931089239486*alpha[1]*favg[45]+0.0883883476483184*(alpha[9]*favg[39]+alpha[6]*favg[36])-0.5*fjump[35]+(0.07905694150420944*alpha[16]+0.0883883476483184*alpha[0])*favg[35]+0.07071067811865474*alpha[7]*favg[33]+(0.1224744871391589*favg[26]+0.07071067811865474*favg[7])*alpha[33]+0.1369306393762915*(alpha[21]*favg[27]+alpha[3]*favg[26])+0.07905694150420944*(alpha[8]*favg[21]+favg[8]*alpha[21])+0.0883883476483184*alpha[1]*favg[18]+0.1369306393762915*alpha[7]*favg[14]+0.07905694150420944*(alpha[3]*favg[7]+favg[3]*alpha[7]); 
  Ghat[36] += alpha[8]*(0.0776323754260148*favg[53]+0.1767766952966368*favg[49])+0.1530931089239486*alpha[2]*favg[45]+0.0883883476483184*alpha[10]*favg[39]-0.5*fjump[36]+0.07905694150420944*alpha[17]*favg[36]+0.0883883476483184*(alpha[0]*favg[36]+alpha[6]*favg[35])+0.07071067811865474*alpha[8]*favg[34]+(0.1224744871391589*favg[27]+0.07071067811865474*favg[8])*alpha[34]+0.1369306393762915*(alpha[3]*favg[27]+alpha[21]*favg[26])+0.07905694150420944*(alpha[7]*favg[21]+favg[7]*alpha[21])+0.0883883476483184*alpha[2]*favg[18]+0.1369306393762915*alpha[8]*favg[14]+0.07905694150420944*(alpha[3]*favg[8]+favg[3]*alpha[8]); 
  Ghat[37] += 0.2338535866733712*alpha[37]*favg[55]+0.0776323754260148*alpha[9]*favg[51]+(0.1344632185501192*favg[28]+0.0776323754260148*favg[9])*alpha[51]+0.1976423537605236*alpha[16]*favg[50]+0.1767766952966368*alpha[9]*favg[47]+0.1369306393762915*alpha[37]*favg[46]+(0.09780759955449389*alpha[37]+0.1530931089239486*alpha[4])*favg[43]+0.07071067811865474*alpha[9]*favg[40]-0.5*fjump[37]+(0.05646924393157818*alpha[16]+0.0883883476483184*alpha[0])*favg[37]+(0.1976423537605236*favg[20]+0.07905694150420944*favg[19]+0.05646924393157818*favg[16]+0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[37]+0.1530931089239486*favg[30]*alpha[33]+0.0883883476483184*(favg[11]*alpha[33]+alpha[10]*favg[31])+(0.1530931089239486*favg[29]+0.0883883476483184*favg[10])*alpha[31]+0.1369306393762915*(alpha[1]*favg[28]+alpha[22]*favg[25])+0.07905694150420944*(alpha[7]*favg[23]+alpha[6]*favg[22]+favg[6]*alpha[22])+0.0883883476483184*alpha[4]*favg[16]+(0.1530931089239486*favg[15]+0.0883883476483184*favg[4])*alpha[16]+0.1369306393762915*alpha[9]*favg[12]+0.07905694150420944*(alpha[1]*favg[9]+favg[1]*alpha[9]); 
  Ghat[38] += 0.2338535866733712*alpha[38]*favg[55]+0.0776323754260148*alpha[10]*favg[52]+(0.1344632185501192*favg[29]+0.0776323754260148*favg[10])*alpha[52]+0.1976423537605236*alpha[17]*favg[50]+0.1767766952966368*alpha[10]*favg[48]+0.1369306393762915*alpha[38]*favg[46]+(0.09780759955449389*alpha[38]+0.1530931089239486*alpha[4])*favg[44]+0.07071067811865474*alpha[10]*favg[41]-0.5*fjump[38]+(0.05646924393157818*alpha[17]+0.0883883476483184*alpha[0])*favg[38]+(0.1976423537605236*favg[20]+0.07905694150420944*favg[19]+0.05646924393157818*favg[17]+0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[38]+0.1530931089239486*favg[30]*alpha[34]+0.0883883476483184*(favg[11]*alpha[34]+alpha[9]*favg[32])+(0.1530931089239486*favg[28]+0.0883883476483184*favg[9])*alpha[32]+0.1369306393762915*(alpha[2]*favg[29]+alpha[22]*favg[25])+0.07905694150420944*(alpha[8]*favg[24]+alpha[6]*favg[22]+favg[6]*alpha[22])+0.0883883476483184*alpha[4]*favg[17]+(0.1530931089239486*favg[15]+0.0883883476483184*favg[4])*alpha[17]+0.1369306393762915*alpha[10]*favg[13]+0.07905694150420944*(alpha[2]*favg[10]+favg[2]*alpha[10]); 
  Ghat[39] += 0.1530931089239486*alpha[4]*favg[45]-0.5*fjump[39]+0.0883883476483184*(alpha[0]*favg[39]+alpha[10]*favg[36]+alpha[9]*favg[35])+0.1369306393762915*alpha[3]*favg[30]+0.07905694150420944*(alpha[8]*favg[24]+alpha[7]*favg[23])+0.0883883476483184*alpha[4]*favg[18]+0.07905694150420944*alpha[3]*favg[11]; 
  Ghat[40] += alpha[9]*(0.0776323754260148*favg[54]+0.1767766952966368*favg[50])+0.1530931089239486*alpha[1]*favg[46]+0.0883883476483184*(alpha[7]*favg[42]+alpha[6]*favg[41])-0.5*fjump[40]+(0.07905694150420944*alpha[16]+0.0883883476483184*alpha[0])*favg[40]+0.07071067811865474*alpha[9]*favg[37]+(0.1224744871391589*favg[28]+0.07071067811865474*favg[9])*alpha[37]+0.1369306393762915*(alpha[22]*favg[29]+alpha[4]*favg[28])+0.07905694150420944*(alpha[10]*favg[22]+favg[10]*alpha[22])+0.0883883476483184*alpha[1]*favg[19]+0.1369306393762915*alpha[9]*favg[15]+0.07905694150420944*(alpha[4]*favg[9]+favg[4]*alpha[9]); 
  Ghat[41] += alpha[10]*(0.0776323754260148*favg[54]+0.1767766952966368*favg[50])+0.1530931089239486*alpha[2]*favg[46]+0.0883883476483184*alpha[8]*favg[42]-0.5*fjump[41]+0.07905694150420944*alpha[17]*favg[41]+0.0883883476483184*(alpha[0]*favg[41]+alpha[6]*favg[40])+0.07071067811865474*alpha[10]*favg[38]+(0.1224744871391589*favg[29]+0.07071067811865474*favg[10])*alpha[38]+0.1369306393762915*(alpha[4]*favg[29]+alpha[22]*favg[28])+0.07905694150420944*(alpha[9]*favg[22]+favg[9]*alpha[22])+0.0883883476483184*alpha[2]*favg[19]+0.1369306393762915*alpha[10]*favg[15]+0.07905694150420944*(alpha[4]*favg[10]+favg[4]*alpha[10]); 
  Ghat[42] += 0.1530931089239486*alpha[3]*favg[46]-0.5*fjump[42]+0.0883883476483184*(alpha[0]*favg[42]+alpha[8]*favg[41]+alpha[7]*favg[40])+0.1369306393762915*alpha[4]*favg[30]+0.07905694150420944*(alpha[10]*favg[24]+alpha[9]*favg[23])+0.0883883476483184*alpha[3]*favg[19]+0.07905694150420944*alpha[4]*favg[11]; 
  Ghat[51] += 0.2338535866733712*alpha[51]*favg[55]-0.5*fjump[51]+(0.05270462766947297*alpha[16]+0.0883883476483184*alpha[0])*favg[51]+(0.09128709291752767*favg[43]+0.1976423537605236*favg[20]+0.05270462766947297*favg[16]+0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[51]+0.1735912687073533*alpha[16]*favg[47]+0.1344632185501192*alpha[1]*favg[43]+0.0776323754260148*alpha[9]*favg[37]+0.1344632185501192*favg[28]*alpha[37]+0.0776323754260148*(favg[9]*alpha[37]+alpha[7]*favg[33])+0.1344632185501192*favg[26]*alpha[33]+0.0776323754260148*(favg[7]*alpha[33]+alpha[6]*favg[31])+0.1344632185501192*favg[25]*alpha[31]+0.0776323754260148*(favg[6]*alpha[31]+alpha[1]*favg[16])+(0.1344632185501192*favg[12]+0.0776323754260148*favg[1])*alpha[16]; 
  Ghat[52] += 0.2338535866733712*alpha[52]*favg[55]-0.5*fjump[52]+(0.05270462766947297*alpha[17]+0.0883883476483184*alpha[0])*favg[52]+(0.09128709291752767*favg[44]+0.1976423537605236*favg[20]+0.05270462766947297*favg[17]+0.1530931089239486*favg[5]+0.0883883476483184*favg[0])*alpha[52]+0.1735912687073533*alpha[17]*favg[48]+0.1344632185501192*alpha[2]*favg[44]+0.0776323754260148*alpha[10]*favg[38]+0.1344632185501192*favg[29]*alpha[38]+0.0776323754260148*(favg[10]*alpha[38]+alpha[8]*favg[34])+0.1344632185501192*favg[27]*alpha[34]+0.0776323754260148*(favg[8]*alpha[34]+alpha[6]*favg[32])+0.1344632185501192*favg[25]*alpha[32]+0.0776323754260148*(favg[6]*alpha[32]+alpha[2]*favg[17])+(0.1344632185501192*favg[13]+0.0776323754260148*favg[2])*alpha[17]; 
  Ghat[53] += (-0.5*fjump[53])+0.0883883476483184*alpha[0]*favg[53]+0.1344632185501192*alpha[3]*favg[45]+0.0776323754260148*(alpha[8]*favg[36]+alpha[7]*favg[35]+alpha[3]*favg[18]); 
  Ghat[54] += (-0.5*fjump[54])+0.0883883476483184*alpha[0]*favg[54]+0.1344632185501192*alpha[4]*favg[46]+0.0776323754260148*(alpha[10]*favg[41]+alpha[9]*favg[40]+alpha[4]*favg[19]); 

  outr[0] += 0.5*Ghat[0]*dv12r; 
  outr[1] += 0.5*Ghat[1]*dv12r; 
  outr[2] += 0.5*Ghat[2]*dv12r; 
  outr[3] += 0.5*Ghat[3]*dv12r; 
  outr[4] += 0.5*Ghat[4]*dv12r; 
  outr[5] += -0.8660254037844386*Ghat[0]*dv12r; 
  outr[6] += 0.5*Ghat[6]*dv12r; 
  outr[7] += 0.5*Ghat[7]*dv12r; 
  outr[8] += 0.5*Ghat[8]*dv12r; 
  outr[9] += 0.5*Ghat[9]*dv12r; 
  outr[10] += 0.5*Ghat[10]*dv12r; 
  outr[11] += 0.5*Ghat[11]*dv12r; 
  outr[12] += -0.8660254037844386*Ghat[1]*dv12r; 
  outr[13] += -0.8660254037844386*Ghat[2]*dv12r; 
  outr[14] += -0.8660254037844386*Ghat[3]*dv12r; 
  outr[15] += -0.8660254037844386*Ghat[4]*dv12r; 
  outr[16] += 0.5*Ghat[16]*dv12r; 
  outr[17] += 0.5*Ghat[17]*dv12r; 
  outr[18] += 0.5*Ghat[18]*dv12r; 
  outr[19] += 0.5*Ghat[19]*dv12r; 
  outr[20] += 1.118033988749895*Ghat[0]*dv12r; 
  outr[21] += 0.5*Ghat[21]*dv12r; 
  outr[22] += 0.5*Ghat[22]*dv12r; 
  outr[23] += 0.5*Ghat[23]*dv12r; 
  outr[24] += 0.5*Ghat[24]*dv12r; 
  outr[25] += -0.8660254037844386*Ghat[6]*dv12r; 
  outr[26] += -0.8660254037844386*Ghat[7]*dv12r; 
  outr[27] += -0.8660254037844386*Ghat[8]*dv12r; 
  outr[28] += -0.8660254037844386*Ghat[9]*dv12r; 
  outr[29] += -0.8660254037844386*Ghat[10]*dv12r; 
  outr[30] += -0.8660254037844386*Ghat[11]*dv12r; 
  outr[31] += 0.5*Ghat[31]*dv12r; 
  outr[32] += 0.5*Ghat[32]*dv12r; 
  outr[33] += 0.5*Ghat[33]*dv12r; 
  outr[34] += 0.5*Ghat[34]*dv12r; 
  outr[35] += 0.5*Ghat[35]*dv12r; 
  outr[36] += 0.5*Ghat[36]*dv12r; 
  outr[37] += 0.5*Ghat[37]*dv12r; 
  outr[38] += 0.5*Ghat[38]*dv12r; 
  outr[39] += 0.5*Ghat[39]*dv12r; 
  outr[40] += 0.5*Ghat[40]*dv12r; 
  outr[41] += 0.5*Ghat[41]*dv12r; 
  outr[42] += 0.5*Ghat[42]*dv12r; 
  outr[43] += -0.8660254037844387*Ghat[16]*dv12r; 
  outr[44] += -0.8660254037844387*Ghat[17]*dv12r; 
  outr[45] += -0.8660254037844387*Ghat[18]*dv12r; 
  outr[46] += -0.8660254037844387*Ghat[19]*dv12r; 
  outr[47] += 1.118033988749895*Ghat[1]*dv12r; 
  outr[48] += 1.118033988749895*Ghat[2]*dv12r; 
  outr[49] += 1.118033988749895*Ghat[3]*dv12r; 
  outr[50] += 1.118033988749895*Ghat[4]*dv12r; 
  outr[51] += 0.5*Ghat[51]*dv12r; 
  outr[52] += 0.5*Ghat[52]*dv12r; 
  outr[53] += 0.5*Ghat[53]*dv12r; 
  outr[54] += 0.5*Ghat[54]*dv12r; 
  outr[55] += -1.322875655532295*Ghat[0]*dv12r; 

  outl[0] += -0.5*Ghat[0]*dv12l; 
  outl[1] += -0.5*Ghat[1]*dv12l; 
  outl[2] += -0.5*Ghat[2]*dv12l; 
  outl[3] += -0.5*Ghat[3]*dv12l; 
  outl[4] += -0.5*Ghat[4]*dv12l; 
  outl[5] += -0.8660254037844386*Ghat[0]*dv12l; 
  outl[6] += -0.5*Ghat[6]*dv12l; 
  outl[7] += -0.5*Ghat[7]*dv12l; 
  outl[8] += -0.5*Ghat[8]*dv12l; 
  outl[9] += -0.5*Ghat[9]*dv12l; 
  outl[10] += -0.5*Ghat[10]*dv12l; 
  outl[11] += -0.5*Ghat[11]*dv12l; 
  outl[12] += -0.8660254037844386*Ghat[1]*dv12l; 
  outl[13] += -0.8660254037844386*Ghat[2]*dv12l; 
  outl[14] += -0.8660254037844386*Ghat[3]*dv12l; 
  outl[15] += -0.8660254037844386*Ghat[4]*dv12l; 
  outl[16] += -0.5*Ghat[16]*dv12l; 
  outl[17] += -0.5*Ghat[17]*dv12l; 
  outl[18] += -0.5*Ghat[18]*dv12l; 
  outl[19] += -0.5*Ghat[19]*dv12l; 
  outl[20] += -1.118033988749895*Ghat[0]*dv12l; 
  outl[21] += -0.5*Ghat[21]*dv12l; 
  outl[22] += -0.5*Ghat[22]*dv12l; 
  outl[23] += -0.5*Ghat[23]*dv12l; 
  outl[24] += -0.5*Ghat[24]*dv12l; 
  outl[25] += -0.8660254037844386*Ghat[6]*dv12l; 
  outl[26] += -0.8660254037844386*Ghat[7]*dv12l; 
  outl[27] += -0.8660254037844386*Ghat[8]*dv12l; 
  outl[28] += -0.8660254037844386*Ghat[9]*dv12l; 
  outl[29] += -0.8660254037844386*Ghat[10]*dv12l; 
  outl[30] += -0.8660254037844386*Ghat[11]*dv12l; 
  outl[31] += -0.5*Ghat[31]*dv12l; 
  outl[32] += -0.5*Ghat[32]*dv12l; 
  outl[33] += -0.5*Ghat[33]*dv12l; 
  outl[34] += -0.5*Ghat[34]*dv12l; 
  outl[35] += -0.5*Ghat[35]*dv12l; 
  outl[36] += -0.5*Ghat[36]*dv12l; 
  outl[37] += -0.5*Ghat[37]*dv12l; 
  outl[38] += -0.5*Ghat[38]*dv12l; 
  outl[39] += -0.5*Ghat[39]*dv12l; 
  outl[40] += -0.5*Ghat[40]*dv12l; 
  outl[41] += -0.5*Ghat[41]*dv12l; 
  outl[42] += -0.5*Ghat[42]*dv12l; 
  outl[43] += -0.8660254037844387*Ghat[16]*dv12l; 
  outl[44] += -0.8660254037844387*Ghat[17]*dv12l; 
  outl[45] += -0.8660254037844387*Ghat[18]*dv12l; 
  outl[46] += -0.8660254037844387*Ghat[19]*dv12l; 
  outl[47] += -1.118033988749895*Ghat[1]*dv12l; 
  outl[48] += -1.118033988749895*Ghat[2]*dv12l; 
  outl[49] += -1.118033988749895*Ghat[3]*dv12l; 
  outl[50] += -1.118033988749895*Ghat[4]*dv12l; 
  outl[51] += -0.5*Ghat[51]*dv12l; 
  outl[52] += -0.5*Ghat[52]*dv12l; 
  outl[53] += -0.5*Ghat[53]*dv12l; 
  outl[54] += -0.5*Ghat[54]*dv12l; 
  outl[55] += -1.322875655532295*Ghat[0]*dv12l; 
return std::abs(amid); 
} 
