#include <VlasovModDecl.h> 
double VlasovSurfElcMag1x1vSer_VX_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double abar0[2]; 


  abar0[0] = E0[0]; 
  abar0[1] = E0[1]; 

  double Ghat[4]; 

  for(unsigned int i=0; i<4; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[4]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  double fjump[4]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  const double amid = 0.7071067811865475*abar0[0]; 
  Ghat[0] = abar0[1]*(0.6123724356957946*favg[3]+0.3535533905932738*favg[1])-0.8660254037844386*fjump[2]+abar0[0]*(0.6123724356957946*favg[2]+0.3535533905932738*favg[0])-0.5*fjump[0]; 
  Ghat[1] = (-0.8660254037844386*fjump[3])+abar0[0]*(0.6123724356957946*favg[3]+0.3535533905932738*favg[1])+abar0[1]*(0.6123724356957946*favg[2]+0.3535533905932738*favg[0])-0.5*fjump[1]; 

  outr[0] += 0.5*Ghat[0]*dv10; 
  outr[1] += 0.5*Ghat[1]*dv10; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10; 
  outr[3] += -0.8660254037844386*Ghat[1]*dv10; 

  outl[0] += -0.5*Ghat[0]*dv10; 
  outl[1] += -0.5*Ghat[1]*dv10; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10; 
  outl[3] += -0.8660254037844386*Ghat[1]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x1vSer_VX_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 


  abar0[0] = E0[0]; 
  abar0[1] = E0[1]; 
  abar0[2] = E0[2]; 

  double Ghat[8]; 

  for(unsigned int i=0; i<8; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[8]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  double fjump[8]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  const double amid = 0.7071067811865475*abar0[0]-0.7905694150420947*abar0[2]; 
  Ghat[0] = abar0[1]*(0.7905694150420948*favg[7]+0.6123724356957945*favg[3]+0.3535533905932738*favg[1])+abar0[2]*(0.6123724356957947*favg[6]+0.3535533905932738*favg[4])-1.118033988749895*fjump[5]+abar0[0]*(0.7905694150420951*favg[5]+0.6123724356957945*favg[2]+0.3535533905932738*favg[0])-0.8660254037844385*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (-1.118033988749895*fjump[7])+abar0[0]*(0.7905694150420948*favg[7]+0.6123724356957945*favg[3]+0.3535533905932738*favg[1])+abar0[2]*(0.7071067811865477*favg[7]+0.5477225575051662*favg[3]+0.316227766016838*favg[1])+abar0[1]*(0.5477225575051662*favg[6]+0.7905694150420951*favg[5]+0.316227766016838*favg[4]+0.6123724356957945*favg[2]+0.3535533905932738*favg[0])-0.8660254037844385*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = abar0[1]*(0.7071067811865478*favg[7]+0.5477225575051663*favg[3]+0.316227766016838*favg[1])-0.8660254037844387*fjump[6]+abar0[0]*(0.6123724356957947*favg[6]+0.3535533905932738*favg[4])+abar0[2]*(0.3912303982179759*favg[6]+0.7905694150420951*favg[5]+0.2258769757263129*favg[4]+0.6123724356957946*favg[2]+0.3535533905932738*favg[0])-0.5*fjump[4]; 

  outr[0] += 0.5*Ghat[0]*dv10; 
  outr[1] += 0.5*Ghat[1]*dv10; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10; 
  outr[3] += -0.8660254037844386*Ghat[1]*dv10; 
  outr[4] += 0.5*Ghat[4]*dv10; 
  outr[5] += 1.118033988749895*Ghat[0]*dv10; 
  outr[6] += -0.8660254037844387*Ghat[4]*dv10; 
  outr[7] += 1.118033988749895*Ghat[1]*dv10; 

  outl[0] += -0.5*Ghat[0]*dv10; 
  outl[1] += -0.5*Ghat[1]*dv10; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10; 
  outl[3] += -0.8660254037844386*Ghat[1]*dv10; 
  outl[4] += -0.5*Ghat[4]*dv10; 
  outl[5] += -1.118033988749895*Ghat[0]*dv10; 
  outl[6] += -0.8660254037844387*Ghat[4]*dv10; 
  outl[7] += -1.118033988749895*Ghat[1]*dv10; 
return std::abs(amid); 
} 
