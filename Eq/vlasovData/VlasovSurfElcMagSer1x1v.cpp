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

  double Ghat[4]; 

  double alpha[4]; 

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
  alpha[0] = 1.414213562373095*E0[0]; 
  alpha[1] = 1.414213562373095*E0[1]; 
  const double amid = 0.5*alpha[0]; 
  Ghat[0] = alpha[1]*(0.4330127018922193*favg[3]+0.25*favg[1])-0.8660254037844386*fjump[2]+alpha[0]*(0.4330127018922193*favg[2]+0.25*favg[0])-0.5*fjump[0]; 
  Ghat[1] = (-0.8660254037844386*fjump[3])+alpha[0]*(0.4330127018922193*favg[3]+0.25*favg[1])+alpha[1]*(0.4330127018922193*favg[2]+0.25*favg[0])-0.5*fjump[1]; 

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

  double Ghat[8]; 

  double alpha[8]; 

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
  alpha[0] = 1.414213562373095*E0[0]; 
  alpha[1] = 1.414213562373095*E0[1]; 
  alpha[4] = 1.414213562373095*E0[2]; 
  const double amid = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 
  Ghat[0] = alpha[1]*(0.5590169943749475*favg[7]+0.4330127018922193*favg[3]+0.25*favg[1])+alpha[4]*(0.4330127018922193*favg[6]+0.25*favg[4])-1.118033988749895*fjump[5]+alpha[0]*(0.5590169943749475*favg[5]+0.4330127018922193*favg[2]+0.25*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (-1.118033988749895*fjump[7])+alpha[0]*(0.5590169943749475*favg[7]+0.4330127018922193*favg[3]+0.25*favg[1])+alpha[4]*(0.5*favg[7]+0.3872983346207416*favg[3]+0.223606797749979*favg[1])+alpha[1]*(0.3872983346207416*favg[6]+0.5590169943749475*favg[5]+0.223606797749979*favg[4]+0.4330127018922193*favg[2]+0.25*favg[0])-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = alpha[1]*(0.5*favg[7]+0.3872983346207416*favg[3]+0.223606797749979*favg[1])-0.8660254037844386*fjump[6]+alpha[0]*(0.4330127018922193*favg[6]+0.25*favg[4])+alpha[4]*(0.276641667586244*favg[6]+0.5590169943749475*favg[5]+0.159719141249985*favg[4]+0.4330127018922193*favg[2]+0.25*favg[0])-0.5*fjump[4]; 

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

  double Ghat[12]; 

  double alpha[12]; 

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
  alpha[0] = 1.414213562373095*E0[0]; 
  alpha[1] = 1.414213562373095*E0[1]; 
  alpha[4] = 1.414213562373095*E0[2]; 
  alpha[8] = 1.414213562373095*E0[3]; 
  const double amid = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 
  Ghat[0] = alpha[1]*(0.6614378277661477*favg[11]+0.5590169943749475*favg[7]+0.4330127018922193*favg[3]+0.25*favg[1])+alpha[8]*(0.4330127018922193*favg[10]+0.25*favg[8])-1.322875655532295*fjump[9]+alpha[0]*(0.6614378277661477*favg[9]+0.5590169943749475*favg[5]+0.4330127018922193*favg[2]+0.25*favg[0])+alpha[4]*(0.4330127018922193*favg[6]+0.25*favg[4])-1.118033988749895*fjump[5]-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (-1.322875655532295*fjump[11])+alpha[0]*(0.6614378277661477*favg[11]+0.5590169943749475*favg[7]+0.4330127018922193*favg[3]+0.25*favg[1])+alpha[4]*(0.5916079783099616*favg[11]+0.3803194146278324*favg[10]+0.2195775164134199*favg[8]+0.5*favg[7]+0.3872983346207416*favg[3]+0.223606797749979*favg[1])+alpha[1]*(0.6614378277661477*favg[9]+0.3872983346207416*favg[6]+0.5590169943749475*favg[5]+0.223606797749979*favg[4]+0.4330127018922193*favg[2]+0.25*favg[0])+(0.3803194146278324*favg[6]+0.2195775164134199*favg[4])*alpha[8]-1.118033988749895*fjump[7]-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = alpha[1]*(0.5916079783099616*favg[11]+0.3803194146278324*favg[10]+0.2195775164134199*favg[8]+0.5*favg[7]+0.3872983346207416*favg[3]+0.223606797749979*favg[1])+alpha[8]*(0.5809475019311124*favg[11]+0.2581988897471612*favg[10]+0.149071198499986*favg[8]+0.4909902530309828*favg[7]+0.3803194146278324*favg[3]+0.2195775164134199*favg[1])+alpha[4]*(0.6614378277661477*favg[9]+0.276641667586244*favg[6]+0.5590169943749475*favg[5]+0.159719141249985*favg[4]+0.4330127018922193*favg[2]+0.25*favg[0])-0.8660254037844386*fjump[6]+alpha[0]*(0.4330127018922193*favg[6]+0.25*favg[4])-0.5*fjump[4]; 
  Ghat[8] = alpha[4]*(0.5809475019311124*favg[11]+0.2581988897471612*favg[10]+0.149071198499986*favg[8]+0.4909902530309828*favg[7]+0.3803194146278324*favg[3]+0.2195775164134199*favg[1])-0.8660254037844386*fjump[10]+alpha[0]*(0.4330127018922193*favg[10]+0.25*favg[8])+alpha[8]*(0.6614378277661477*favg[9]+0.2581988897471612*favg[6]+0.5590169943749475*favg[5]+0.149071198499986*favg[4]+0.4330127018922193*favg[2]+0.25*favg[0])-0.5*fjump[8]+alpha[1]*(0.3803194146278324*favg[6]+0.2195775164134199*favg[4]); 

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
