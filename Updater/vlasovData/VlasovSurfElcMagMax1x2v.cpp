#include <VlasovModDecl.h> 
double VlasovSurfElcMag1x2vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double abar0[2]; 

  double abar1[2]; 


  abar0[0] = E0[0]+wv2*B2[0]; 
  abar0[1] = E0[1]+wv2*B2[1]; 

  double incr[4]; 

  double favg[4]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  double fjump[4]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = -1*fr[2]-fl[2]; 
  fjump[3] = 1*fr[3]-fl[3]; 
  const double amid = 0.7071067811865475*abar0[0]; 
  incr[0] = 0.05103103630798285*B2[0]*favg[3]*dv2-0.4330127018922193*fjump[2]*amax-0.25*fjump[0]*amax+0.3061862178478971*abar0[0]*favg[2]+0.1767766952966368*abar0[1]*favg[1]+0.1767766952966368*abar0[0]*favg[0]; 
  incr[1] = 0.05103103630798285*B2[1]*favg[3]*dv2-0.25*fjump[1]*amax+0.3061862178478971*abar0[1]*favg[2]+0.1767766952966368*abar0[0]*favg[1]+0.1767766952966368*favg[0]*abar0[1]; 
  incr[2] = (-0.0883883476483184*B2[0]*favg[3]*dv2)+0.75*fjump[2]*amax+0.4330127018922193*fjump[0]*amax-0.5303300858899105*abar0[0]*favg[2]-0.3061862178478971*abar0[1]*favg[1]-0.3061862178478971*abar0[0]*favg[0]; 
  incr[3] = 0.0883883476483184*B2[0]*favg[2]*dv2+0.05103103630798285*favg[1]*B2[1]*dv2+0.05103103630798285*favg[0]*B2[0]*dv2-0.25*fjump[3]*amax+0.1767766952966368*abar0[0]*favg[3]; 

  outr[0] += incr[0]*dv10; 
  outr[1] += incr[1]*dv10; 
  outr[2] += incr[2]*dv10; 
  outr[3] += incr[3]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += incr[2]*dv10; 
  outl[3] += -1.0*incr[3]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 


  abar0[0] = E0[0]+wv2*B2[0]; 
  abar0[1] = E0[1]+wv2*B2[1]; 
  abar0[2] = E0[2]+wv2*B2[2]; 

  double incr[10]; 

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

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = -1*fr[2]-fl[2]; 
  fjump[3] = 1*fr[3]-fl[3]; 
  fjump[4] = -1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  fjump[6] = -1*fr[6]-fl[6]; 
  fjump[7] = 1*fr[7]-fl[7]; 
  fjump[8] = 1*fr[8]-fl[8]; 
  fjump[9] = 1*fr[9]-fl[9]; 
  const double amid = 0.7071067811865475*abar0[0]-0.7905694150420947*abar0[2]; 
  incr[0] = 0.0883883476483184*B2[0]*favg[6]*dv2+0.05103103630798285*B2[1]*favg[5]*dv2+0.05103103630798285*B2[0]*favg[3]*dv2-0.5590169943749475*fjump[8]*amax-0.4330127018922193*fjump[2]*amax-0.25*fjump[0]*amax+0.3952847075210473*abar0[0]*favg[8]+0.1767766952966368*abar0[2]*favg[7]+0.3061862178478971*abar0[1]*favg[4]+0.3061862178478971*abar0[0]*favg[2]+0.1767766952966368*abar0[1]*favg[1]+0.1767766952966368*abar0[0]*favg[0]; 
  incr[1] = 0.0883883476483184*B2[1]*favg[6]*dv2+0.04564354645876383*B2[2]*favg[5]*dv2+0.05103103630798285*B2[0]*favg[5]*dv2+0.05103103630798285*B2[1]*favg[3]*dv2-0.4330127018922193*fjump[4]*amax-0.25*fjump[1]*amax+0.3952847075210473*abar0[1]*favg[8]+0.1581138830084189*abar0[1]*favg[7]+0.2738612787525829*abar0[2]*favg[4]+0.3061862178478971*abar0[0]*favg[4]+0.3061862178478971*abar0[1]*favg[2]+0.1581138830084189*favg[1]*abar0[2]+0.1767766952966368*abar0[0]*favg[1]+0.1767766952966368*favg[0]*abar0[1]; 
  incr[2] = (-0.1530931089239486*B2[0]*favg[6]*dv2)-0.0883883476483184*B2[1]*favg[5]*dv2-0.0883883476483184*B2[0]*favg[3]*dv2+0.9682458365518543*fjump[8]*amax+0.75*fjump[2]*amax+0.4330127018922193*fjump[0]*amax-0.6846531968814573*abar0[0]*favg[8]-0.3061862178478971*abar0[2]*favg[7]-0.5303300858899105*abar0[1]*favg[4]-0.5303300858899105*abar0[0]*favg[2]-0.3061862178478971*abar0[1]*favg[1]-0.3061862178478971*abar0[0]*favg[0]; 
  incr[3] = 0.04564354645876383*B2[0]*favg[9]*dv2+0.1141088661469096*B2[0]*favg[8]*dv2+0.05103103630798285*B2[2]*favg[7]*dv2+0.0883883476483184*B2[1]*favg[4]*dv2+0.0883883476483184*B2[0]*favg[2]*dv2+0.05103103630798285*favg[1]*B2[1]*dv2+0.05103103630798285*favg[0]*B2[0]*dv2-0.4330127018922193*fjump[6]*amax-0.25*fjump[3]*amax+0.3061862178478971*abar0[0]*favg[6]+0.1767766952966368*abar0[1]*favg[5]+0.1767766952966368*abar0[0]*favg[3]; 
  incr[4] = (-0.1530931089239486*B2[1]*favg[6]*dv2)-0.07905694150420946*B2[2]*favg[5]*dv2-0.0883883476483184*B2[0]*favg[5]*dv2-0.0883883476483184*B2[1]*favg[3]*dv2+0.75*fjump[4]*amax+0.4330127018922193*fjump[1]*amax-0.6846531968814574*abar0[1]*favg[8]-0.2738612787525829*abar0[1]*favg[7]-0.4743416490252568*abar0[2]*favg[4]-0.5303300858899105*abar0[0]*favg[4]-0.5303300858899105*abar0[1]*favg[2]-0.2738612787525829*favg[1]*abar0[2]-0.3061862178478971*abar0[0]*favg[1]-0.3061862178478971*favg[0]*abar0[1]; 
  incr[5] = 0.04564354645876383*B2[1]*favg[9]*dv2+0.1141088661469096*B2[1]*favg[8]*dv2+0.04564354645876383*B2[1]*favg[7]*dv2+0.07905694150420946*B2[2]*favg[4]*dv2+0.0883883476483184*B2[0]*favg[4]*dv2+0.04564354645876383*favg[1]*B2[2]*dv2+0.0883883476483184*B2[1]*favg[2]*dv2+0.05103103630798285*favg[0]*B2[1]*dv2+0.05103103630798285*B2[0]*favg[1]*dv2-0.25*fjump[5]*amax+0.3061862178478971*abar0[1]*favg[6]+0.1581138830084189*abar0[2]*favg[5]+0.1767766952966368*abar0[0]*favg[5]+0.1767766952966368*abar0[1]*favg[3]; 
  incr[6] = (-0.07905694150420946*B2[0]*favg[9]*dv2)-0.1976423537605237*B2[0]*favg[8]*dv2-0.0883883476483184*B2[2]*favg[7]*dv2-0.1530931089239486*B2[1]*favg[4]*dv2-0.1530931089239486*B2[0]*favg[2]*dv2-0.0883883476483184*favg[1]*B2[1]*dv2-0.0883883476483184*favg[0]*B2[0]*dv2+0.75*fjump[6]*amax+0.4330127018922193*fjump[3]*amax-0.5303300858899105*abar0[0]*favg[6]-0.3061862178478971*abar0[1]*favg[5]-0.3061862178478971*abar0[0]*favg[3]; 
  incr[7] = 0.0883883476483184*B2[2]*favg[6]*dv2+0.04564354645876383*B2[1]*favg[5]*dv2+0.05103103630798286*B2[2]*favg[3]*dv2-0.25*fjump[7]*amax+0.3952847075210473*abar0[2]*favg[8]+0.1129384878631564*abar0[2]*favg[7]+0.1767766952966368*abar0[0]*favg[7]+0.273861278752583*abar0[1]*favg[4]+0.3061862178478971*abar0[2]*favg[2]+0.1767766952966368*favg[0]*abar0[2]+0.1581138830084189*abar0[1]*favg[1]; 
  incr[8] = 0.1976423537605236*B2[0]*favg[6]*dv2+0.1141088661469096*B2[1]*favg[5]*dv2+0.1141088661469096*B2[0]*favg[3]*dv2-1.25*fjump[8]*amax-0.9682458365518543*fjump[2]*amax-0.5590169943749475*fjump[0]*amax+0.883883476483184*abar0[0]*favg[8]+0.3952847075210473*abar0[2]*favg[7]+0.6846531968814572*abar0[1]*favg[4]+0.6846531968814572*abar0[0]*favg[2]+0.3952847075210473*abar0[1]*favg[1]+0.3952847075210473*abar0[0]*favg[0]; 
  incr[9] = 0.07905694150420947*B2[0]*favg[6]*dv2+0.04564354645876383*B2[1]*favg[5]*dv2+0.04564354645876383*B2[0]*favg[3]*dv2-0.25*fjump[9]*amax+0.1767766952966368*abar0[0]*favg[9]; 

  outr[0] += incr[0]*dv10; 
  outr[1] += incr[1]*dv10; 
  outr[2] += incr[2]*dv10; 
  outr[3] += incr[3]*dv10; 
  outr[4] += incr[4]*dv10; 
  outr[5] += incr[5]*dv10; 
  outr[6] += incr[6]*dv10; 
  outr[7] += incr[7]*dv10; 
  outr[8] += incr[8]*dv10; 
  outr[9] += incr[9]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += incr[2]*dv10; 
  outl[3] += -1.0*incr[3]*dv10; 
  outl[4] += incr[4]*dv10; 
  outl[5] += -1.0*incr[5]*dv10; 
  outl[6] += incr[6]*dv10; 
  outl[7] += -1.0*incr[7]*dv10; 
  outl[8] += -1.0*incr[8]*dv10; 
  outl[9] += -1.0*incr[9]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vMax_VY_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[2]; 

  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double abar0[2]; 

  double abar1[2]; 


  abar1[0] = E1[0]-wv1*B2[0]; 
  abar1[1] = E1[1]-wv1*B2[1]; 

  double incr[4]; 

  double favg[4]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  double fjump[4]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  const double amid = 0.7071067811865475*abar1[0]; 
  incr[0] = (-0.05103103630798285*B2[0]*favg[2]*dv1)-0.4330127018922193*fjump[3]*amax-0.25*fjump[0]*amax+0.3061862178478971*abar1[0]*favg[3]+0.1767766952966368*abar1[1]*favg[1]+0.1767766952966368*abar1[0]*favg[0]; 
  incr[1] = (-0.05103103630798285*B2[1]*favg[2]*dv1)-0.25*fjump[1]*amax+0.3061862178478971*abar1[1]*favg[3]+0.1767766952966368*abar1[0]*favg[1]+0.1767766952966368*favg[0]*abar1[1]; 
  incr[2] = (-0.0883883476483184*B2[0]*favg[3]*dv1)-0.05103103630798285*favg[1]*B2[1]*dv1-0.05103103630798285*favg[0]*B2[0]*dv1-0.25*fjump[2]*amax+0.1767766952966368*abar1[0]*favg[2]; 
  incr[3] = 0.0883883476483184*B2[0]*favg[2]*dv1+0.75*fjump[3]*amax+0.4330127018922193*fjump[0]*amax-0.5303300858899105*abar1[0]*favg[3]-0.3061862178478971*abar1[1]*favg[1]-0.3061862178478971*abar1[0]*favg[0]; 

  outr[0] += incr[0]*dv11; 
  outr[1] += incr[1]*dv11; 
  outr[2] += incr[2]*dv11; 
  outr[3] += incr[3]*dv11; 

  outl[0] += -1.0*incr[0]*dv11; 
  outl[1] += -1.0*incr[1]*dv11; 
  outl[2] += -1.0*incr[2]*dv11; 
  outl[3] += incr[3]*dv11; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x2vMax_VY_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[3]; 

  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 


  abar1[0] = E1[0]-wv1*B2[0]; 
  abar1[1] = E1[1]-wv1*B2[1]; 
  abar1[2] = E1[2]-wv1*B2[2]; 

  double incr[10]; 

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

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = -1*fr[5]-fl[5]; 
  fjump[6] = -1*fr[6]-fl[6]; 
  fjump[7] = 1*fr[7]-fl[7]; 
  fjump[8] = 1*fr[8]-fl[8]; 
  fjump[9] = 1*fr[9]-fl[9]; 
  const double amid = 0.7071067811865475*abar1[0]-0.7905694150420947*abar1[2]; 
  incr[0] = (-0.0883883476483184*B2[0]*favg[6]*dv1)-0.05103103630798285*B2[1]*favg[4]*dv1-0.05103103630798285*B2[0]*favg[2]*dv1-0.5590169943749475*fjump[9]*amax-0.4330127018922193*fjump[3]*amax-0.25*fjump[0]*amax+0.3952847075210473*abar1[0]*favg[9]+0.1767766952966368*abar1[2]*favg[7]+0.3061862178478971*abar1[1]*favg[5]+0.3061862178478971*abar1[0]*favg[3]+0.1767766952966368*abar1[1]*favg[1]+0.1767766952966368*abar1[0]*favg[0]; 
  incr[1] = (-0.0883883476483184*B2[1]*favg[6]*dv1)-0.04564354645876383*B2[2]*favg[4]*dv1-0.05103103630798285*B2[0]*favg[4]*dv1-0.05103103630798285*B2[1]*favg[2]*dv1-0.4330127018922193*fjump[5]*amax-0.25*fjump[1]*amax+0.3952847075210473*abar1[1]*favg[9]+0.1581138830084189*abar1[1]*favg[7]+0.2738612787525829*abar1[2]*favg[5]+0.3061862178478971*abar1[0]*favg[5]+0.3061862178478971*abar1[1]*favg[3]+0.1581138830084189*favg[1]*abar1[2]+0.1767766952966368*abar1[0]*favg[1]+0.1767766952966368*favg[0]*abar1[1]; 
  incr[2] = (-0.1141088661469096*B2[0]*favg[9]*dv1)-0.04564354645876383*B2[0]*favg[8]*dv1-0.05103103630798285*B2[2]*favg[7]*dv1-0.0883883476483184*B2[1]*favg[5]*dv1-0.0883883476483184*B2[0]*favg[3]*dv1-0.05103103630798285*favg[1]*B2[1]*dv1-0.05103103630798285*favg[0]*B2[0]*dv1-0.4330127018922193*fjump[6]*amax-0.25*fjump[2]*amax+0.3061862178478971*abar1[0]*favg[6]+0.1767766952966368*abar1[1]*favg[4]+0.1767766952966368*abar1[0]*favg[2]; 
  incr[3] = 0.1530931089239486*B2[0]*favg[6]*dv1+0.0883883476483184*B2[1]*favg[4]*dv1+0.0883883476483184*B2[0]*favg[2]*dv1+0.9682458365518543*fjump[9]*amax+0.75*fjump[3]*amax+0.4330127018922193*fjump[0]*amax-0.6846531968814573*abar1[0]*favg[9]-0.3061862178478971*abar1[2]*favg[7]-0.5303300858899105*abar1[1]*favg[5]-0.5303300858899105*abar1[0]*favg[3]-0.3061862178478971*abar1[1]*favg[1]-0.3061862178478971*abar1[0]*favg[0]; 
  incr[4] = (-0.1141088661469096*B2[1]*favg[9]*dv1)-0.04564354645876383*B2[1]*favg[8]*dv1-0.04564354645876383*B2[1]*favg[7]*dv1-0.07905694150420946*B2[2]*favg[5]*dv1-0.0883883476483184*B2[0]*favg[5]*dv1-0.0883883476483184*B2[1]*favg[3]*dv1-0.04564354645876383*favg[1]*B2[2]*dv1-0.05103103630798285*favg[0]*B2[1]*dv1-0.05103103630798285*B2[0]*favg[1]*dv1-0.25*fjump[4]*amax+0.3061862178478971*abar1[1]*favg[6]+0.1581138830084189*abar1[2]*favg[4]+0.1767766952966368*abar1[0]*favg[4]+0.1767766952966368*abar1[1]*favg[2]; 
  incr[5] = 0.1530931089239486*B2[1]*favg[6]*dv1+0.07905694150420946*B2[2]*favg[4]*dv1+0.0883883476483184*B2[0]*favg[4]*dv1+0.0883883476483184*B2[1]*favg[2]*dv1+0.75*fjump[5]*amax+0.4330127018922193*fjump[1]*amax-0.6846531968814574*abar1[1]*favg[9]-0.2738612787525829*abar1[1]*favg[7]-0.4743416490252568*abar1[2]*favg[5]-0.5303300858899105*abar1[0]*favg[5]-0.5303300858899105*abar1[1]*favg[3]-0.2738612787525829*favg[1]*abar1[2]-0.3061862178478971*abar1[0]*favg[1]-0.3061862178478971*favg[0]*abar1[1]; 
  incr[6] = 0.1976423537605237*B2[0]*favg[9]*dv1+0.07905694150420946*B2[0]*favg[8]*dv1+0.0883883476483184*B2[2]*favg[7]*dv1+0.1530931089239486*B2[1]*favg[5]*dv1+0.1530931089239486*B2[0]*favg[3]*dv1+0.0883883476483184*favg[1]*B2[1]*dv1+0.0883883476483184*favg[0]*B2[0]*dv1+0.75*fjump[6]*amax+0.4330127018922193*fjump[2]*amax-0.5303300858899105*abar1[0]*favg[6]-0.3061862178478971*abar1[1]*favg[4]-0.3061862178478971*abar1[0]*favg[2]; 
  incr[7] = (-0.0883883476483184*B2[2]*favg[6]*dv1)-0.04564354645876383*B2[1]*favg[4]*dv1-0.05103103630798286*favg[2]*B2[2]*dv1-0.25*fjump[7]*amax+0.3952847075210473*abar1[2]*favg[9]+0.1129384878631564*abar1[2]*favg[7]+0.1767766952966368*abar1[0]*favg[7]+0.273861278752583*abar1[1]*favg[5]+0.3061862178478971*abar1[2]*favg[3]+0.1767766952966368*favg[0]*abar1[2]+0.1581138830084189*abar1[1]*favg[1]; 
  incr[8] = (-0.07905694150420947*B2[0]*favg[6]*dv1)-0.04564354645876383*B2[1]*favg[4]*dv1-0.04564354645876383*B2[0]*favg[2]*dv1-0.25*fjump[8]*amax+0.1767766952966368*abar1[0]*favg[8]; 
  incr[9] = (-0.1976423537605236*B2[0]*favg[6]*dv1)-0.1141088661469096*B2[1]*favg[4]*dv1-0.1141088661469096*B2[0]*favg[2]*dv1-1.25*fjump[9]*amax-0.9682458365518543*fjump[3]*amax-0.5590169943749475*fjump[0]*amax+0.883883476483184*abar1[0]*favg[9]+0.3952847075210473*abar1[2]*favg[7]+0.6846531968814572*abar1[1]*favg[5]+0.6846531968814572*abar1[0]*favg[3]+0.3952847075210473*abar1[1]*favg[1]+0.3952847075210473*abar1[0]*favg[0]; 

  outr[0] += incr[0]*dv11; 
  outr[1] += incr[1]*dv11; 
  outr[2] += incr[2]*dv11; 
  outr[3] += incr[3]*dv11; 
  outr[4] += incr[4]*dv11; 
  outr[5] += incr[5]*dv11; 
  outr[6] += incr[6]*dv11; 
  outr[7] += incr[7]*dv11; 
  outr[8] += incr[8]*dv11; 
  outr[9] += incr[9]*dv11; 

  outl[0] += -1.0*incr[0]*dv11; 
  outl[1] += -1.0*incr[1]*dv11; 
  outl[2] += -1.0*incr[2]*dv11; 
  outl[3] += incr[3]*dv11; 
  outl[4] += -1.0*incr[4]*dv11; 
  outl[5] += incr[5]*dv11; 
  outl[6] += incr[6]*dv11; 
  outl[7] += -1.0*incr[7]*dv11; 
  outl[8] += -1.0*incr[8]*dv11; 
  outl[9] += -1.0*incr[9]*dv11; 
return std::abs(amid); 
} 
