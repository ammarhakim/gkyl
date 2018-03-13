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

  double incr[4]; 

  double favg[4]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  double fjump[4]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = -1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  const double amid = 0.7071067811865475*abar0[0]; 
  incr[0] = (-0.4330127018922193*fjump[2]*amax)-0.25*fjump[0]*amax+0.3061862178478973*abar0[1]*favg[3]+0.3061862178478973*abar0[0]*favg[2]+0.1767766952966369*abar0[1]*favg[1]+0.1767766952966369*abar0[0]*favg[0]; 
  incr[1] = (-0.4330127018922193*fjump[3]*amax)-0.25*fjump[1]*amax+0.3061862178478973*abar0[0]*favg[3]+0.3061862178478973*abar0[1]*favg[2]+0.1767766952966369*abar0[0]*favg[1]+0.1767766952966369*favg[0]*abar0[1]; 
  incr[2] = 0.75*fjump[2]*amax+0.4330127018922193*fjump[0]*amax-0.5303300858899107*abar0[1]*favg[3]-0.5303300858899107*abar0[0]*favg[2]-0.3061862178478973*abar0[1]*favg[1]-0.3061862178478973*abar0[0]*favg[0]; 
  incr[3] = 0.75*fjump[3]*amax+0.4330127018922193*fjump[1]*amax-0.5303300858899107*abar0[0]*favg[3]-0.5303300858899107*abar0[1]*favg[2]-0.3061862178478973*abar0[0]*favg[1]-0.3061862178478973*favg[0]*abar0[1]; 

  outr[0] += incr[0]*dv10; 
  outr[1] += incr[1]*dv10; 
  outr[2] += incr[2]*dv10; 
  outr[3] += incr[3]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += incr[2]*dv10; 
  outl[3] += incr[3]*dv10; 
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

  double incr[8]; 

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

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = -1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  fjump[6] = -1*fr[6]-fl[6]; 
  fjump[7] = 1*fr[7]-fl[7]; 
  const double amid = 0.7071067811865475*abar0[0]-0.7905694150420947*abar0[2]; 
  incr[0] = (-0.5590169943749475*fjump[5]*amax)-0.4330127018922193*fjump[2]*amax-0.25*fjump[0]*amax+0.3952847075210474*abar0[1]*favg[7]+0.3061862178478973*abar0[2]*favg[6]+0.3952847075210476*abar0[0]*favg[5]+0.1767766952966369*abar0[2]*favg[4]+0.3061862178478972*abar0[1]*favg[3]+0.3061862178478972*abar0[0]*favg[2]+0.1767766952966369*abar0[1]*favg[1]+0.1767766952966369*abar0[0]*favg[0]; 
  incr[1] = (-0.5590169943749475*fjump[7]*amax)-0.4330127018922192*fjump[3]*amax-0.25*fjump[1]*amax+0.3535533905932738*abar0[2]*favg[7]+0.3952847075210474*abar0[0]*favg[7]+0.2738612787525831*abar0[1]*favg[6]+0.3952847075210476*abar0[1]*favg[5]+0.158113883008419*abar0[1]*favg[4]+0.2738612787525831*abar0[2]*favg[3]+0.3061862178478972*abar0[0]*favg[3]+0.3061862178478972*abar0[1]*favg[2]+0.158113883008419*favg[1]*abar0[2]+0.1767766952966369*abar0[0]*favg[1]+0.1767766952966369*favg[0]*abar0[1]; 
  incr[2] = 0.9682458365518543*fjump[5]*amax+0.75*fjump[2]*amax+0.4330127018922193*fjump[0]*amax-0.6846531968814578*abar0[1]*favg[7]-0.5303300858899107*abar0[2]*favg[6]-0.6846531968814579*abar0[0]*favg[5]-0.3061862178478973*abar0[2]*favg[4]-0.5303300858899107*abar0[1]*favg[3]-0.5303300858899107*abar0[0]*favg[2]-0.3061862178478973*abar0[1]*favg[1]-0.3061862178478973*abar0[0]*favg[0]; 
  incr[3] = 0.9682458365518543*fjump[7]*amax+0.75*fjump[3]*amax+0.4330127018922193*fjump[1]*amax-0.6123724356957947*abar0[2]*favg[7]-0.6846531968814578*abar0[0]*favg[7]-0.474341649025257*abar0[1]*favg[6]-0.6846531968814579*abar0[1]*favg[5]-0.2738612787525831*abar0[1]*favg[4]-0.4743416490252571*abar0[2]*favg[3]-0.5303300858899107*abar0[0]*favg[3]-0.5303300858899107*abar0[1]*favg[2]-0.2738612787525831*favg[1]*abar0[2]-0.3061862178478973*abar0[0]*favg[1]-0.3061862178478973*favg[0]*abar0[1]; 
  incr[4] = (-0.4330127018922194*fjump[6]*amax)-0.25*fjump[4]*amax+0.3535533905932739*abar0[1]*favg[7]+0.195615199108988*abar0[2]*favg[6]+0.3061862178478973*abar0[0]*favg[6]+0.3952847075210476*abar0[2]*favg[5]+0.1129384878631565*abar0[2]*favg[4]+0.1767766952966369*abar0[0]*favg[4]+0.2738612787525831*abar0[1]*favg[3]+0.3061862178478973*abar0[2]*favg[2]+0.1767766952966369*favg[0]*abar0[2]+0.158113883008419*abar0[1]*favg[1]; 
  incr[5] = (-1.25*fjump[5]*amax)-0.9682458365518543*fjump[2]*amax-0.5590169943749475*fjump[0]*amax+0.8838834764831843*abar0[1]*favg[7]+0.6846531968814578*abar0[2]*favg[6]+0.8838834764831844*abar0[0]*favg[5]+0.3952847075210474*abar0[2]*favg[4]+0.6846531968814576*abar0[1]*favg[3]+0.6846531968814576*abar0[0]*favg[2]+0.3952847075210474*abar0[1]*favg[1]+0.3952847075210474*abar0[0]*favg[0]; 
  incr[6] = 0.75*fjump[6]*amax+0.4330127018922194*fjump[4]*amax-0.6123724356957946*abar0[1]*favg[7]-0.3388154635894693*abar0[2]*favg[6]-0.5303300858899107*abar0[0]*favg[6]-0.6846531968814578*abar0[2]*favg[5]-0.1956151991089879*abar0[2]*favg[4]-0.3061862178478973*abar0[0]*favg[4]-0.474341649025257*abar0[1]*favg[3]-0.5303300858899106*abar0[2]*favg[2]-0.3061862178478973*favg[0]*abar0[2]-0.2738612787525831*abar0[1]*favg[1]; 
  incr[7] = (-1.25*fjump[7]*amax)-0.9682458365518543*fjump[3]*amax-0.5590169943749475*fjump[1]*amax+0.7905694150420951*abar0[2]*favg[7]+0.8838834764831844*abar0[0]*favg[7]+0.6123724356957946*abar0[1]*favg[6]+0.8838834764831844*abar0[1]*favg[5]+0.3535533905932738*abar0[1]*favg[4]+0.6123724356957947*abar0[2]*favg[3]+0.6846531968814578*abar0[0]*favg[3]+0.6846531968814578*abar0[1]*favg[2]+0.3535533905932738*favg[1]*abar0[2]+0.3952847075210474*abar0[0]*favg[1]+0.3952847075210474*favg[0]*abar0[1]; 

  outr[0] += incr[0]*dv10; 
  outr[1] += incr[1]*dv10; 
  outr[2] += incr[2]*dv10; 
  outr[3] += incr[3]*dv10; 
  outr[4] += incr[4]*dv10; 
  outr[5] += incr[5]*dv10; 
  outr[6] += incr[6]*dv10; 
  outr[7] += incr[7]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += incr[2]*dv10; 
  outl[3] += incr[3]*dv10; 
  outl[4] += -1.0*incr[4]*dv10; 
  outl[5] += -1.0*incr[5]*dv10; 
  outl[6] += incr[6]*dv10; 
  outl[7] += -1.0*incr[7]*dv10; 
return std::abs(amid); 
} 
