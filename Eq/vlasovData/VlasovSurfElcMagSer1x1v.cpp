#include <VlasovModDecl.h> 
double VlasovSurfElcMag1x1vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 

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

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[3] += -0.8660254037844386*Ghat[1]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[3] += -0.8660254037844386*Ghat[1]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x1vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 

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

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[3] += -0.8660254037844386*Ghat[1]*dv10r; 
  outr[4] += 0.5*Ghat[4]*dv10r; 
  outr[5] += 1.118033988749895*Ghat[0]*dv10r; 
  outr[6] += -0.8660254037844387*Ghat[4]*dv10r; 
  outr[7] += 1.118033988749895*Ghat[1]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[3] += -0.8660254037844386*Ghat[1]*dv10l; 
  outl[4] += -0.5*Ghat[4]*dv10l; 
  outl[5] += -1.118033988749895*Ghat[0]*dv10l; 
  outl[6] += -0.8660254037844387*Ghat[4]*dv10l; 
  outl[7] += -1.118033988749895*Ghat[1]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x1vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double abar0[4]; 


  abar0[0] = E0[0]; 
  abar0[1] = E0[1]; 
  abar0[2] = E0[2]; 
  abar0[3] = E0[3]; 

  double Ghat[12]; 

  for(unsigned int i=0; i<12; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[12]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  double fjump[12]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(-1*fr[9]-fl[9]); 
  fjump[10] = amax*(-1*fr[10]-fl[10]); 
  fjump[11] = amax*(-1*fr[11]-fl[11]); 
  const double amid = 0.7071067811865475*abar0[0]-0.7905694150420947*abar0[2]; 
  Ghat[0] = abar0[1]*(0.9354143466934854*favg[11]+0.790569415042095*favg[7]+0.6123724356957946*favg[3]+0.3535533905932738*favg[1])+abar0[3]*(0.6123724356957947*favg[10]+0.3535533905932738*favg[8])-1.322875655532296*fjump[9]+abar0[0]*(0.9354143466934857*favg[9]+0.7905694150420951*favg[5]+0.6123724356957946*favg[2]+0.3535533905932738*favg[0])+abar0[2]*(0.6123724356957947*favg[6]+0.3535533905932738*favg[4])-1.118033988749895*fjump[5]-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (-1.322875655532295*fjump[11])+abar0[0]*(0.9354143466934854*favg[11]+0.790569415042095*favg[7]+0.6123724356957946*favg[3]+0.3535533905932738*favg[1])+abar0[2]*(0.8366600265340758*favg[11]+0.5378528742004771*favg[10]+0.3105295017040594*favg[8]+0.7071067811865478*favg[7]+0.5477225575051663*favg[3]+0.316227766016838*favg[1])+abar0[1]*(0.9354143466934857*favg[9]+0.5477225575051662*favg[6]+0.7905694150420951*favg[5]+0.316227766016838*favg[4]+0.6123724356957946*favg[2]+0.3535533905932738*favg[0])-1.118033988749895*fjump[7]+abar0[3]*(0.5378528742004771*favg[6]+0.3105295017040594*favg[4])-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = abar0[1]*(0.8366600265340758*favg[11]+0.5378528742004771*favg[10]+0.3105295017040594*favg[8]+0.7071067811865478*favg[7]+0.5477225575051663*favg[3]+0.316227766016838*favg[1])+abar0[3]*(0.8215838362577493*favg[11]+0.365148371670111*favg[10]+0.2108185106778921*favg[8]+0.6943650748294138*favg[7]+0.5378528742004773*favg[3]+0.3105295017040594*favg[1])+abar0[2]*(0.9354143466934857*favg[9]+0.3912303982179759*favg[6]+0.7905694150420951*favg[5]+0.2258769757263129*favg[4]+0.6123724356957946*favg[2]+0.3535533905932738*favg[0])-0.8660254037844387*fjump[6]+abar0[0]*(0.6123724356957947*favg[6]+0.3535533905932738*favg[4])-0.5*fjump[4]; 
  Ghat[8] = abar0[2]*(0.8215838362577493*favg[11]+0.365148371670111*favg[10]+0.2108185106778921*favg[8]+0.6943650748294138*favg[7]+0.5378528742004773*favg[3]+0.3105295017040594*favg[1])-0.8660254037844388*fjump[10]+abar0[0]*(0.6123724356957947*favg[10]+0.3535533905932738*favg[8])+abar0[3]*(0.9354143466934857*favg[9]+0.365148371670111*favg[6]+0.7905694150420951*favg[5]+0.2108185106778921*favg[4]+0.6123724356957946*favg[2]+0.3535533905932738*favg[0])-0.5*fjump[8]+abar0[1]*(0.5378528742004771*favg[6]+0.3105295017040594*favg[4]); 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[3] += -0.8660254037844386*Ghat[1]*dv10r; 
  outr[4] += 0.5*Ghat[4]*dv10r; 
  outr[5] += 1.118033988749895*Ghat[0]*dv10r; 
  outr[6] += -0.8660254037844387*Ghat[4]*dv10r; 
  outr[7] += 1.118033988749895*Ghat[1]*dv10r; 
  outr[8] += 0.5*Ghat[8]*dv10r; 
  outr[9] += -1.322875655532295*Ghat[0]*dv10r; 
  outr[10] += -0.8660254037844386*Ghat[8]*dv10r; 
  outr[11] += -1.322875655532295*Ghat[1]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[3] += -0.8660254037844386*Ghat[1]*dv10l; 
  outl[4] += -0.5*Ghat[4]*dv10l; 
  outl[5] += -1.118033988749895*Ghat[0]*dv10l; 
  outl[6] += -0.8660254037844387*Ghat[4]*dv10l; 
  outl[7] += -1.118033988749895*Ghat[1]*dv10l; 
  outl[8] += -0.5*Ghat[8]*dv10l; 
  outl[9] += -1.322875655532295*Ghat[0]*dv10l; 
  outl[10] += -0.8660254037844386*Ghat[8]*dv10l; 
  outl[11] += -1.322875655532295*Ghat[1]*dv10l; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x1vSer_VX_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 

  const double *B0 = &EM[15]; 
  const double *B1 = &EM[20]; 
  const double *B2 = &EM[25]; 

  double abar0[5]; 


  abar0[0] = E0[0]; 
  abar0[1] = E0[1]; 
  abar0[2] = E0[2]; 
  abar0[3] = E0[3]; 
  abar0[4] = E0[4]; 

  double Ghat[17]; 

  for(unsigned int i=0; i<17; ++i){ 

    Ghat[i]=0.0; 

  }; 

  double favg[17]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = -1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  double fjump[17]; 

  fjump[0] = amax*(1*fr[0]-fl[0]); 
  fjump[1] = amax*(1*fr[1]-fl[1]); 
  fjump[2] = amax*(-1*fr[2]-fl[2]); 
  fjump[3] = amax*(-1*fr[3]-fl[3]); 
  fjump[4] = amax*(1*fr[4]-fl[4]); 
  fjump[5] = amax*(1*fr[5]-fl[5]); 
  fjump[6] = amax*(-1*fr[6]-fl[6]); 
  fjump[7] = amax*(1*fr[7]-fl[7]); 
  fjump[8] = amax*(1*fr[8]-fl[8]); 
  fjump[9] = amax*(-1*fr[9]-fl[9]); 
  fjump[10] = amax*(1*fr[10]-fl[10]); 
  fjump[11] = amax*(-1*fr[11]-fl[11]); 
  fjump[12] = amax*(-1*fr[12]-fl[12]); 
  fjump[13] = amax*(1*fr[13]-fl[13]); 
  fjump[14] = amax*(1*fr[14]-fl[14]); 
  fjump[15] = amax*(-1*fr[15]-fl[15]); 
  fjump[16] = amax*(1*fr[16]-fl[16]); 
  const double amid = 0.7954951288348656*abar0[4]-0.7905694150420947*abar0[2]+0.7071067811865475*abar0[0]; 
  Ghat[0] = abar0[1]*(1.060660171779821*favg[16]+0.9354143466934854*favg[12]+0.790569415042095*favg[7]+0.6123724356957946*favg[3]+0.3535533905932738*favg[1])+abar0[4]*(0.6123724356957946*favg[15]+0.3535533905932738*favg[13])-1.5*fjump[14]+abar0[0]*(1.060660171779821*favg[14]+0.9354143466934857*favg[9]+0.7905694150420951*favg[5]+0.6123724356957946*favg[2]+0.3535533905932738*favg[0])+abar0[3]*(0.6123724356957947*favg[11]+0.3535533905932738*favg[8])+abar0[2]*(0.7905694150420951*favg[10]+0.6123724356957947*favg[6]+0.3535533905932738*favg[4])-1.322875655532296*fjump[9]-1.118033988749895*fjump[5]-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (-1.5*fjump[16])+abar0[0]*(1.060660171779821*favg[16]+0.9354143466934854*favg[12]+0.790569415042095*favg[7]+0.6123724356957946*favg[3]+0.3535533905932738*favg[1])+abar0[2]*(0.9486832980505142*favg[16]+0.8366600265340758*favg[12]+0.5378528742004771*favg[11]+0.3105295017040594*favg[8]+0.7071067811865478*favg[7]+0.5477225575051663*favg[3]+0.316227766016838*favg[1])+abar0[3]*(0.534522483824849*favg[15]+0.308606699924184*favg[13]+0.6943650748294135*favg[10]+0.5378528742004771*favg[6]+0.3105295017040594*favg[4])+abar0[1]*(1.060660171779821*favg[14]+0.7071067811865477*favg[10]+0.9354143466934857*favg[9]+0.5477225575051662*favg[6]+0.7905694150420951*favg[5]+0.316227766016838*favg[4]+0.6123724356957946*favg[2]+0.3535533905932738*favg[0])-1.322875655532295*fjump[12]+abar0[4]*(0.534522483824849*favg[11]+0.308606699924184*favg[8])-1.118033988749895*fjump[7]-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = abar0[1]*(0.9486832980505141*favg[16]+0.8366600265340757*favg[12]+0.5378528742004771*favg[11]+0.3105295017040594*favg[8]+0.7071067811865477*favg[7]+0.5477225575051662*favg[3]+0.316227766016838*favg[1])+abar0[3]*(0.9315885051121781*favg[16]+0.8215838362577493*favg[12]+0.365148371670111*favg[11]+0.2108185106778921*favg[8]+0.6943650748294138*favg[7]+0.5378528742004771*favg[3]+0.3105295017040594*favg[1])+abar0[2]*(0.5248906591678237*favg[15]+1.060660171779821*favg[14]+0.3030457633656633*favg[13]+0.5050762722761055*favg[10]+0.9354143466934856*favg[9]+0.3912303982179759*favg[6]+0.7905694150420951*favg[5]+0.2258769757263129*favg[4]+0.6123724356957945*favg[2]+0.3535533905932738*favg[0])+abar0[4]*(0.3556639983799781*favg[15]+0.2053427052057391*favg[13]+0.6776309271789387*favg[10]+0.5248906591678241*favg[6]+0.3030457633656633*favg[4])-1.118033988749895*fjump[10]+abar0[0]*(0.7905694150420951*favg[10]+0.6123724356957947*favg[6]+0.3535533905932738*favg[4])-0.8660254037844387*fjump[6]-0.5*fjump[4]; 
  Ghat[8] = abar0[2]*(0.9315885051121781*favg[16]+0.8215838362577493*favg[12]+0.365148371670111*favg[11]+0.2108185106778921*favg[8]+0.6943650748294138*favg[7]+0.5378528742004771*favg[3]+0.3105295017040594*favg[1])+abar0[4]*(0.9258200997725516*favg[16]+0.8164965809277264*favg[12]+0.3340213285613425*favg[11]+0.1928473039599675*favg[8]+0.6900655593423546*favg[7]+0.534522483824849*favg[3]+0.3086066999241839*favg[1])+abar0[1]*(0.534522483824849*favg[15]+0.3086066999241839*favg[13]+0.6943650748294133*favg[10]+0.5378528742004771*favg[6]+0.3105295017040594*favg[4])+abar0[3]*(0.3340213285613425*favg[15]+1.060660171779821*favg[14]+0.1928473039599675*favg[13]+0.4714045207910318*favg[10]+0.9354143466934856*favg[9]+0.3651483716701109*favg[6]+0.7905694150420951*favg[5]+0.2108185106778921*favg[4]+0.6123724356957945*favg[2]+0.3535533905932738*favg[0])-0.8660254037844386*fjump[11]+abar0[0]*(0.6123724356957946*favg[11]+0.3535533905932738*favg[8])-0.5*fjump[8]; 
  Ghat[13] = abar0[3]*(0.9258200997725516*favg[16]+0.8164965809277264*favg[12]+0.3340213285613425*favg[11]+0.1928473039599675*favg[8]+0.6900655593423546*favg[7]+0.534522483824849*favg[3]+0.3086066999241839*favg[1])-0.8660254037844385*fjump[15]+abar0[0]*(0.6123724356957945*favg[15]+0.3535533905932738*favg[13])+abar0[2]*(0.3556639983799781*favg[15]+0.2053427052057391*favg[13]+0.6776309271789387*favg[10]+0.5248906591678241*favg[6]+0.3030457633656633*favg[4])+abar0[4]*(0.2973156880600959*favg[15]+1.060660171779821*favg[14]+0.1716552925357953*favg[13]+0.4591602475237324*favg[10]+0.9354143466934856*favg[9]+0.3556639983799782*favg[6]+0.7905694150420951*favg[5]+0.2053427052057391*favg[4]+0.6123724356957945*favg[2]+0.3535533905932738*favg[0])-0.5*fjump[13]+abar0[1]*(0.5345224838248489*favg[11]+0.3086066999241839*favg[8]); 

  outr[0] += 0.5*Ghat[0]*dv10r; 
  outr[1] += 0.5*Ghat[1]*dv10r; 
  outr[2] += -0.8660254037844386*Ghat[0]*dv10r; 
  outr[3] += -0.8660254037844386*Ghat[1]*dv10r; 
  outr[4] += 0.5*Ghat[4]*dv10r; 
  outr[5] += 1.118033988749895*Ghat[0]*dv10r; 
  outr[6] += -0.8660254037844387*Ghat[4]*dv10r; 
  outr[7] += 1.118033988749895*Ghat[1]*dv10r; 
  outr[8] += 0.5*Ghat[8]*dv10r; 
  outr[9] += -1.322875655532295*Ghat[0]*dv10r; 
  outr[10] += 1.118033988749895*Ghat[4]*dv10r; 
  outr[11] += -0.8660254037844386*Ghat[8]*dv10r; 
  outr[12] += -1.322875655532295*Ghat[1]*dv10r; 
  outr[13] += 0.5*Ghat[13]*dv10r; 
  outr[14] += 1.5*Ghat[0]*dv10r; 
  outr[15] += -0.8660254037844386*Ghat[13]*dv10r; 
  outr[16] += 1.5*Ghat[1]*dv10r; 

  outl[0] += -0.5*Ghat[0]*dv10l; 
  outl[1] += -0.5*Ghat[1]*dv10l; 
  outl[2] += -0.8660254037844386*Ghat[0]*dv10l; 
  outl[3] += -0.8660254037844386*Ghat[1]*dv10l; 
  outl[4] += -0.5*Ghat[4]*dv10l; 
  outl[5] += -1.118033988749895*Ghat[0]*dv10l; 
  outl[6] += -0.8660254037844387*Ghat[4]*dv10l; 
  outl[7] += -1.118033988749895*Ghat[1]*dv10l; 
  outl[8] += -0.5*Ghat[8]*dv10l; 
  outl[9] += -1.322875655532295*Ghat[0]*dv10l; 
  outl[10] += -1.118033988749895*Ghat[4]*dv10l; 
  outl[11] += -0.8660254037844386*Ghat[8]*dv10l; 
  outl[12] += -1.322875655532295*Ghat[1]*dv10l; 
  outl[13] += -0.5*Ghat[13]*dv10l; 
  outl[14] += -1.5*Ghat[0]*dv10l; 
  outl[15] += -0.8660254037844386*Ghat[13]*dv10l; 
  outl[16] += -1.5*Ghat[1]*dv10l; 
return std::abs(amid); 
} 
