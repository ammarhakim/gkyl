#include <VlasovModDecl.h> 
double VlasovSurfElcMag1x1vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr[3]; 

  double favg[3]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  double fjump[3]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = -1*fr[2]-fl[2]; 
  const double amid = 0.7071067811865475*abar0[0]; 
  incr[0] = (-0.4330127018922193*fjump[2]*amax)-0.25*fjump[0]*amax+0.3061862178478973*abar0[0]*favg[2]+0.1767766952966369*abar0[1]*favg[1]+0.1767766952966369*abar0[0]*favg[0]; 
  incr[1] = (-0.25*fjump[1]*amax)+0.3061862178478973*abar0[1]*favg[2]+0.1767766952966369*abar0[0]*favg[1]+0.1767766952966369*favg[0]*abar0[1]; 
  incr[2] = 0.75*fjump[2]*amax+0.4330127018922193*fjump[0]*amax-0.5303300858899107*abar0[0]*favg[2]-0.3061862178478973*abar0[1]*favg[1]-0.3061862178478973*abar0[0]*favg[0]; 

  outr[0] += incr[0]*dv10; 
  outr[1] += incr[1]*dv10; 
  outr[2] += incr[2]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += incr[2]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x1vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr[6]; 

  double favg[6]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  double fjump[6]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = -1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  const double amid = 0.7071067811865475*abar0[0]-0.7905694150420947*abar0[2]; 
  incr[0] = (-0.5590169943749475*fjump[5]*amax)-0.4330127018922193*fjump[2]*amax-0.25*fjump[0]*amax+0.3952847075210474*abar0[0]*favg[5]+0.1767766952966369*abar0[2]*favg[4]+0.3061862178478973*abar0[1]*favg[3]+0.3061862178478973*abar0[0]*favg[2]+0.1767766952966369*abar0[1]*favg[1]+0.1767766952966369*abar0[0]*favg[0]; 
  incr[1] = (-0.4330127018922193*fjump[3]*amax)-0.25*fjump[1]*amax+0.3952847075210476*abar0[1]*favg[5]+0.158113883008419*abar0[1]*favg[4]+0.2738612787525831*abar0[2]*favg[3]+0.3061862178478973*abar0[0]*favg[3]+0.3061862178478973*abar0[1]*favg[2]+0.158113883008419*favg[1]*abar0[2]+0.1767766952966369*abar0[0]*favg[1]+0.1767766952966369*favg[0]*abar0[1]; 
  incr[2] = 0.9682458365518543*fjump[5]*amax+0.75*fjump[2]*amax+0.4330127018922193*fjump[0]*amax-0.6846531968814578*abar0[0]*favg[5]-0.3061862178478973*abar0[2]*favg[4]-0.5303300858899107*abar0[1]*favg[3]-0.5303300858899107*abar0[0]*favg[2]-0.3061862178478973*abar0[1]*favg[1]-0.3061862178478973*abar0[0]*favg[0]; 
  incr[3] = 0.75*fjump[3]*amax+0.4330127018922193*fjump[1]*amax-0.6846531968814579*abar0[1]*favg[5]-0.2738612787525831*abar0[1]*favg[4]-0.4743416490252571*abar0[2]*favg[3]-0.5303300858899107*abar0[0]*favg[3]-0.5303300858899107*abar0[1]*favg[2]-0.2738612787525831*favg[1]*abar0[2]-0.3061862178478973*abar0[0]*favg[1]-0.3061862178478973*favg[0]*abar0[1]; 
  incr[4] = (-0.25*fjump[4]*amax)+0.3952847075210475*abar0[2]*favg[5]+0.1129384878631564*abar0[2]*favg[4]+0.1767766952966369*abar0[0]*favg[4]+0.2738612787525831*abar0[1]*favg[3]+0.3061862178478973*abar0[2]*favg[2]+0.1767766952966369*favg[0]*abar0[2]+0.158113883008419*abar0[1]*favg[1]; 
  incr[5] = (-1.25*fjump[5]*amax)-0.9682458365518543*fjump[2]*amax-0.5590169943749475*fjump[0]*amax+0.8838834764831844*abar0[0]*favg[5]+0.3952847075210474*abar0[2]*favg[4]+0.6846531968814578*abar0[1]*favg[3]+0.6846531968814578*abar0[0]*favg[2]+0.3952847075210474*abar0[1]*favg[1]+0.3952847075210474*abar0[0]*favg[0]; 

  outr[0] += incr[0]*dv10; 
  outr[1] += incr[1]*dv10; 
  outr[2] += incr[2]*dv10; 
  outr[3] += incr[3]*dv10; 
  outr[4] += incr[4]*dv10; 
  outr[5] += incr[5]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += incr[2]*dv10; 
  outl[3] += incr[3]*dv10; 
  outl[4] += -1.0*incr[4]*dv10; 
  outl[5] += -1.0*incr[5]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x1vMax_VX_P3(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double abar0[4]; 


  abar0[0] = E0[0]; 
  abar0[1] = E0[1]; 
  abar0[2] = E0[2]; 
  abar0[3] = E0[3]; 

  double incr[10]; 

  double favg[10]; 

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
  double fjump[10]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = -1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  fjump[6] = -1*fr[6]-fl[6]; 
  fjump[7] = 1*fr[7]-fl[7]; 
  fjump[8] = 1*fr[8]-fl[8]; 
  fjump[9] = -1*fr[9]-fl[9]; 
  const double amid = 0.7071067811865475*abar0[0]-0.7905694150420947*abar0[2]; 
  incr[0] = (-0.6614378277661477*fjump[9]*amax)-0.5590169943749475*fjump[5]*amax-0.4330127018922193*fjump[2]*amax-0.25*fjump[0]*amax+0.4677071733467427*abar0[0]*favg[9]+0.1767766952966369*abar0[3]*favg[8]+0.3952847075210474*abar0[1]*favg[7]+0.3061862178478973*abar0[2]*favg[6]+0.3952847075210476*abar0[0]*favg[5]+0.1767766952966369*abar0[2]*favg[4]+0.3061862178478972*abar0[1]*favg[3]+0.3061862178478972*abar0[0]*favg[2]+0.1767766952966369*abar0[1]*favg[1]+0.1767766952966369*abar0[0]*favg[0]; 
  incr[1] = (-0.5590169943749475*fjump[7]*amax)-0.4330127018922192*fjump[3]*amax-0.25*fjump[1]*amax+0.4677071733467428*abar0[1]*favg[9]+0.1552647508520297*abar0[2]*favg[8]+0.3535533905932739*abar0[2]*favg[7]+0.3952847075210475*abar0[0]*favg[7]+0.2689264371002386*abar0[3]*favg[6]+0.2738612787525831*abar0[1]*favg[6]+0.3952847075210476*abar0[1]*favg[5]+0.1552647508520297*abar0[3]*favg[4]+0.158113883008419*abar0[1]*favg[4]+0.2738612787525831*abar0[2]*favg[3]+0.3061862178478973*abar0[0]*favg[3]+0.3061862178478973*abar0[1]*favg[2]+0.158113883008419*favg[1]*abar0[2]+0.1767766952966369*abar0[0]*favg[1]+0.1767766952966369*favg[0]*abar0[1]; 
  incr[2] = 1.14564392373896*fjump[9]*amax+0.9682458365518543*fjump[5]*amax+0.75*fjump[2]*amax+0.4330127018922193*fjump[0]*amax-0.8100925873009827*abar0[0]*favg[9]-0.3061862178478973*abar0[3]*favg[8]-0.6846531968814578*abar0[1]*favg[7]-0.5303300858899107*abar0[2]*favg[6]-0.6846531968814579*abar0[0]*favg[5]-0.3061862178478973*abar0[2]*favg[4]-0.5303300858899107*abar0[1]*favg[3]-0.5303300858899107*abar0[0]*favg[2]-0.3061862178478973*abar0[1]*favg[1]-0.3061862178478973*abar0[0]*favg[0]; 
  incr[3] = 0.9682458365518543*fjump[7]*amax+0.75*fjump[3]*amax+0.4330127018922193*fjump[1]*amax-0.8100925873009827*abar0[1]*favg[9]-0.2689264371002386*abar0[2]*favg[8]-0.6123724356957947*abar0[2]*favg[7]-0.6846531968814578*abar0[0]*favg[7]-0.4657942525560891*abar0[3]*favg[6]-0.474341649025257*abar0[1]*favg[6]-0.6846531968814578*abar0[1]*favg[5]-0.2689264371002386*abar0[3]*favg[4]-0.2738612787525831*abar0[1]*favg[4]-0.4743416490252571*abar0[2]*favg[3]-0.5303300858899107*abar0[0]*favg[3]-0.5303300858899107*abar0[1]*favg[2]-0.2738612787525831*favg[1]*abar0[2]-0.3061862178478973*abar0[0]*favg[1]-0.3061862178478973*favg[0]*abar0[1]; 
  incr[4] = (-0.4330127018922194*fjump[6]*amax)-0.25*fjump[4]*amax+0.4677071733467428*abar0[2]*favg[9]+0.105409255338946*abar0[3]*favg[8]+0.1552647508520297*abar0[1]*favg[8]+0.3471825374147069*abar0[3]*favg[7]+0.3535533905932739*abar0[1]*favg[7]+0.195615199108988*abar0[2]*favg[6]+0.3061862178478973*abar0[0]*favg[6]+0.3952847075210476*abar0[2]*favg[5]+0.1129384878631565*abar0[2]*favg[4]+0.1767766952966369*abar0[0]*favg[4]+0.2689264371002387*abar0[3]*favg[3]+0.2738612787525831*abar0[1]*favg[3]+0.1552647508520297*favg[1]*abar0[3]+0.3061862178478973*abar0[2]*favg[2]+0.1767766952966369*favg[0]*abar0[2]+0.158113883008419*abar0[1]*favg[1]; 
  incr[5] = (-1.479019945774904*fjump[9]*amax)-1.25*fjump[5]*amax-0.9682458365518543*fjump[2]*amax-0.5590169943749475*fjump[0]*amax+1.045825033167594*abar0[0]*favg[9]+0.3952847075210474*abar0[3]*favg[8]+0.8838834764831843*abar0[1]*favg[7]+0.6846531968814578*abar0[2]*favg[6]+0.8838834764831844*abar0[0]*favg[5]+0.3952847075210474*abar0[2]*favg[4]+0.6846531968814576*abar0[1]*favg[3]+0.6846531968814576*abar0[0]*favg[2]+0.3952847075210474*abar0[1]*favg[1]+0.3952847075210474*abar0[0]*favg[0]; 
  incr[6] = 0.75*fjump[6]*amax+0.4330127018922194*fjump[4]*amax-0.8100925873009828*abar0[2]*favg[9]-0.1825741858350555*abar0[3]*favg[8]-0.2689264371002386*abar0[1]*favg[8]-0.601337794302955*abar0[3]*favg[7]-0.6123724356957947*abar0[1]*favg[7]-0.3388154635894693*abar0[2]*favg[6]-0.5303300858899107*abar0[0]*favg[6]-0.6846531968814579*abar0[2]*favg[5]-0.195615199108988*abar0[2]*favg[4]-0.3061862178478973*abar0[0]*favg[4]-0.4657942525560892*abar0[3]*favg[3]-0.4743416490252571*abar0[1]*favg[3]-0.2689264371002386*favg[1]*abar0[3]-0.5303300858899107*abar0[2]*favg[2]-0.3061862178478973*favg[0]*abar0[2]-0.2738612787525831*abar0[1]*favg[1]; 
  incr[7] = (-1.25*fjump[7]*amax)-0.9682458365518543*fjump[3]*amax-0.5590169943749475*fjump[1]*amax+1.045825033167595*abar0[1]*favg[9]+0.3471825374147069*abar0[2]*favg[8]+0.7905694150420951*abar0[2]*favg[7]+0.8838834764831847*abar0[0]*favg[7]+0.601337794302955*abar0[3]*favg[6]+0.6123724356957947*abar0[1]*favg[6]+0.8838834764831848*abar0[1]*favg[5]+0.3471825374147069*abar0[3]*favg[4]+0.3535533905932739*abar0[1]*favg[4]+0.6123724356957947*abar0[2]*favg[3]+0.6846531968814579*abar0[0]*favg[3]+0.6846531968814579*abar0[1]*favg[2]+0.3535533905932739*favg[1]*abar0[2]+0.3952847075210475*abar0[0]*favg[1]+0.3952847075210475*favg[0]*abar0[1]; 
  incr[8] = (-0.25*fjump[8]*amax)+0.4677071733467428*abar0[3]*favg[9]+0.105409255338946*abar0[2]*favg[8]+0.1767766952966369*abar0[0]*favg[8]+0.3471825374147069*abar0[2]*favg[7]+0.1825741858350555*abar0[3]*favg[6]+0.2689264371002386*abar0[1]*favg[6]+0.3952847075210476*abar0[3]*favg[5]+0.105409255338946*abar0[3]*favg[4]+0.1552647508520297*abar0[1]*favg[4]+0.2689264371002387*abar0[2]*favg[3]+0.3061862178478973*favg[2]*abar0[3]+0.1767766952966369*favg[0]*abar0[3]+0.1552647508520297*favg[1]*abar0[2]; 
  incr[9] = 1.75*fjump[9]*amax+1.479019945774904*fjump[5]*amax+1.14564392373896*fjump[2]*amax+0.6614378277661477*fjump[0]*amax-1.237436867076458*abar0[0]*favg[9]-0.4677071733467427*abar0[3]*favg[8]-1.045825033167594*abar0[1]*favg[7]-0.8100925873009828*abar0[2]*favg[6]-1.045825033167595*abar0[0]*favg[5]-0.4677071733467427*abar0[2]*favg[4]-0.8100925873009825*abar0[1]*favg[3]-0.8100925873009825*abar0[0]*favg[2]-0.4677071733467427*abar0[1]*favg[1]-0.4677071733467427*abar0[0]*favg[0]; 

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
  outl[3] += incr[3]*dv10; 
  outl[4] += -1.0*incr[4]*dv10; 
  outl[5] += -1.0*incr[5]*dv10; 
  outl[6] += incr[6]*dv10; 
  outl[7] += -1.0*incr[7]*dv10; 
  outl[8] += -1.0*incr[8]*dv10; 
  outl[9] += incr[9]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x1vMax_VX_P4(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[15]; 
  const double *B1 = &EM[20]; 
  const double *B2 = &EM[25]; 

  double abar0[5]; 


  abar0[0] = E0[0]; 
  abar0[1] = E0[1]; 
  abar0[2] = E0[2]; 
  abar0[3] = E0[3]; 
  abar0[4] = E0[4]; 

  double incr[15]; 

  double favg[15]; 

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
  double fjump[15]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = -1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  fjump[6] = -1*fr[6]-fl[6]; 
  fjump[7] = 1*fr[7]-fl[7]; 
  fjump[8] = 1*fr[8]-fl[8]; 
  fjump[9] = -1*fr[9]-fl[9]; 
  fjump[10] = 1*fr[10]-fl[10]; 
  fjump[11] = -1*fr[11]-fl[11]; 
  fjump[12] = -1*fr[12]-fl[12]; 
  fjump[13] = 1*fr[13]-fl[13]; 
  fjump[14] = 1*fr[14]-fl[14]; 
  const double amid = 0.7954951288348656*abar0[4]-0.7905694150420947*abar0[2]+0.7071067811865475*abar0[0]; 
  incr[0] = (-0.75*fjump[14]*amax)-0.6614378277661477*fjump[9]*amax-0.5590169943749475*fjump[5]*amax-0.4330127018922193*fjump[2]*amax-0.25*fjump[0]*amax+0.5303300858899107*abar0[0]*favg[14]+0.1767766952966369*abar0[4]*favg[13]+0.4677071733467427*abar0[1]*favg[12]+0.3061862178478973*abar0[3]*favg[11]+0.3952847075210476*abar0[2]*favg[10]+0.4677071733467428*abar0[0]*favg[9]+0.1767766952966369*abar0[3]*favg[8]+0.3952847075210475*abar0[1]*favg[7]+0.3061862178478973*abar0[2]*favg[6]+0.3952847075210476*abar0[0]*favg[5]+0.1767766952966369*abar0[2]*favg[4]+0.3061862178478973*abar0[1]*favg[3]+0.3061862178478973*abar0[0]*favg[2]+0.1767766952966369*abar0[1]*favg[1]+0.1767766952966369*abar0[0]*favg[0]; 
  incr[1] = (-0.6614378277661476*fjump[12]*amax)-0.5590169943749475*fjump[7]*amax-0.4330127018922192*fjump[3]*amax-0.25*fjump[1]*amax+0.5303300858899107*abar0[1]*favg[14]+0.154303349962092*abar0[3]*favg[13]+0.4183300132670379*abar0[2]*favg[12]+0.4677071733467427*abar0[0]*favg[12]+0.2672612419124245*abar0[4]*favg[11]+0.2689264371002386*abar0[2]*favg[11]+0.3471825374147067*abar0[3]*favg[10]+0.3535533905932738*abar0[1]*favg[10]+0.4677071733467428*abar0[1]*favg[9]+0.154303349962092*abar0[4]*favg[8]+0.1552647508520297*abar0[2]*favg[8]+0.3535533905932739*abar0[2]*favg[7]+0.3952847075210475*abar0[0]*favg[7]+0.2689264371002386*abar0[3]*favg[6]+0.2738612787525831*abar0[1]*favg[6]+0.3952847075210476*abar0[1]*favg[5]+0.1552647508520297*abar0[3]*favg[4]+0.158113883008419*abar0[1]*favg[4]+0.2738612787525831*abar0[2]*favg[3]+0.3061862178478973*abar0[0]*favg[3]+0.3061862178478973*abar0[1]*favg[2]+0.158113883008419*favg[1]*abar0[2]+0.1767766952966369*abar0[0]*favg[1]+0.1767766952966369*favg[0]*abar0[1]; 
  incr[2] = 1.299038105676658*fjump[14]*amax+1.14564392373896*fjump[9]*amax+0.9682458365518543*fjump[5]*amax+0.75*fjump[2]*amax+0.4330127018922193*fjump[0]*amax-0.9185586535436917*abar0[0]*favg[14]-0.3061862178478973*abar0[4]*favg[13]-0.8100925873009825*abar0[1]*favg[12]-0.5303300858899107*abar0[3]*favg[11]-0.6846531968814578*abar0[2]*favg[10]-0.8100925873009827*abar0[0]*favg[9]-0.3061862178478973*abar0[3]*favg[8]-0.6846531968814578*abar0[1]*favg[7]-0.5303300858899106*abar0[2]*favg[6]-0.6846531968814578*abar0[0]*favg[5]-0.3061862178478973*abar0[2]*favg[4]-0.5303300858899107*abar0[1]*favg[3]-0.5303300858899107*abar0[0]*favg[2]-0.3061862178478973*abar0[1]*favg[1]-0.3061862178478973*abar0[0]*favg[0]; 
  incr[3] = 1.14564392373896*fjump[12]*amax+0.9682458365518543*fjump[7]*amax+0.75*fjump[3]*amax+0.4330127018922193*fjump[1]*amax-0.9185586535436917*abar0[1]*favg[14]-0.2672612419124245*abar0[3]*favg[13]-0.7245688373094721*abar0[2]*favg[12]-0.8100925873009825*abar0[0]*favg[12]-0.4629100498862759*abar0[4]*favg[11]-0.4657942525560891*abar0[2]*favg[11]-0.6013377943029549*abar0[3]*favg[10]-0.6123724356957946*abar0[1]*favg[10]-0.8100925873009827*abar0[1]*favg[9]-0.2672612419124245*abar0[4]*favg[8]-0.2689264371002386*abar0[2]*favg[8]-0.6123724356957947*abar0[2]*favg[7]-0.6846531968814578*abar0[0]*favg[7]-0.4657942525560891*abar0[3]*favg[6]-0.474341649025257*abar0[1]*favg[6]-0.6846531968814578*abar0[1]*favg[5]-0.2689264371002386*abar0[3]*favg[4]-0.2738612787525831*abar0[1]*favg[4]-0.4743416490252571*abar0[2]*favg[3]-0.5303300858899107*abar0[0]*favg[3]-0.5303300858899107*abar0[1]*favg[2]-0.2738612787525831*favg[1]*abar0[2]-0.3061862178478973*abar0[0]*favg[1]-0.3061862178478973*favg[0]*abar0[1]; 
  incr[4] = (-0.5590169943749476*fjump[10]*amax)-0.4330127018922194*fjump[6]*amax-0.25*fjump[4]*amax+0.5303300858899107*abar0[2]*favg[14]+0.1026713526028695*abar0[4]*favg[13]+0.1515228816828317*abar0[2]*favg[13]+0.4107919181288747*abar0[3]*favg[12]+0.4183300132670378*abar0[1]*favg[12]+0.1825741858350555*abar0[3]*favg[11]+0.2689264371002386*abar0[1]*favg[11]+0.3388154635894693*abar0[4]*favg[10]+0.2525381361380528*abar0[2]*favg[10]+0.3952847075210476*abar0[0]*favg[10]+0.4677071733467428*abar0[2]*favg[9]+0.105409255338946*abar0[3]*favg[8]+0.1552647508520297*abar0[1]*favg[8]+0.3471825374147069*abar0[3]*favg[7]+0.3535533905932738*abar0[1]*favg[7]+0.262445329583912*abar0[4]*favg[6]+0.195615199108988*abar0[2]*favg[6]+0.3061862178478973*abar0[0]*favg[6]+0.3952847075210476*abar0[2]*favg[5]+0.1515228816828317*abar0[4]*favg[4]+0.1129384878631564*abar0[2]*favg[4]+0.1767766952966369*abar0[0]*favg[4]+0.2689264371002386*abar0[3]*favg[3]+0.2738612787525831*abar0[1]*favg[3]+0.1552647508520297*favg[1]*abar0[3]+0.3061862178478972*abar0[2]*favg[2]+0.1767766952966369*favg[0]*abar0[2]+0.158113883008419*abar0[1]*favg[1]; 
  incr[5] = (-1.677050983124842*fjump[14]*amax)-1.479019945774904*fjump[9]*amax-1.25*fjump[5]*amax-0.9682458365518543*fjump[2]*amax-0.5590169943749475*fjump[0]*amax+1.185854122563142*abar0[0]*favg[14]+0.3952847075210474*abar0[4]*favg[13]+1.045825033167594*abar0[1]*favg[12]+0.6846531968814576*abar0[3]*favg[11]+0.8838834764831844*abar0[2]*favg[10]+1.045825033167595*abar0[0]*favg[9]+0.3952847075210474*abar0[3]*favg[8]+0.8838834764831843*abar0[1]*favg[7]+0.6846531968814578*abar0[2]*favg[6]+0.8838834764831844*abar0[0]*favg[5]+0.3952847075210474*abar0[2]*favg[4]+0.6846531968814576*abar0[1]*favg[3]+0.6846531968814576*abar0[0]*favg[2]+0.3952847075210474*abar0[1]*favg[1]+0.3952847075210474*abar0[0]*favg[0]; 
  incr[6] = 0.9682458365518543*fjump[10]*amax+0.75*fjump[6]*amax+0.4330127018922194*fjump[4]*amax-0.9185586535436919*abar0[2]*favg[14]-0.1778319991899891*abar0[4]*favg[13]-0.262445329583912*abar0[2]*favg[13]-0.7115124735378855*abar0[3]*favg[12]-0.7245688373094721*abar0[1]*favg[12]-0.3162277660168381*abar0[3]*favg[11]-0.465794252556089*abar0[1]*favg[11]-0.5868455973269638*abar0[4]*favg[10]-0.4374088826398534*abar0[2]*favg[10]-0.6846531968814578*abar0[0]*favg[10]-0.8100925873009825*abar0[2]*favg[9]-0.1825741858350555*abar0[3]*favg[8]-0.2689264371002386*abar0[1]*favg[8]-0.601337794302955*abar0[3]*favg[7]-0.6123724356957946*abar0[1]*favg[7]-0.454568645048495*abar0[4]*favg[6]-0.3388154635894693*abar0[2]*favg[6]-0.5303300858899107*abar0[0]*favg[6]-0.6846531968814578*abar0[2]*favg[5]-0.262445329583912*abar0[4]*favg[4]-0.195615199108988*abar0[2]*favg[4]-0.3061862178478973*abar0[0]*favg[4]-0.4657942525560891*abar0[3]*favg[3]-0.474341649025257*abar0[1]*favg[3]-0.2689264371002386*favg[1]*abar0[3]-0.5303300858899106*abar0[2]*favg[2]-0.3061862178478973*favg[0]*abar0[2]-0.2738612787525831*abar0[1]*favg[1]; 
  incr[7] = (-1.479019945774904*fjump[12]*amax)-1.25*fjump[7]*amax-0.9682458365518543*fjump[3]*amax-0.5590169943749475*fjump[1]*amax+1.185854122563142*abar0[1]*favg[14]+0.3450327796711773*abar0[3]*favg[13]+0.9354143466934858*abar0[2]*favg[12]+1.045825033167595*abar0[0]*favg[12]+0.5976143046671971*abar0[4]*favg[11]+0.601337794302955*abar0[2]*favg[11]+0.7763237542601487*abar0[3]*favg[10]+0.7905694150420951*abar0[1]*favg[10]+1.045825033167595*abar0[1]*favg[9]+0.3450327796711773*abar0[4]*favg[8]+0.3471825374147069*abar0[2]*favg[8]+0.7905694150420951*abar0[2]*favg[7]+0.8838834764831847*abar0[0]*favg[7]+0.601337794302955*abar0[3]*favg[6]+0.6123724356957947*abar0[1]*favg[6]+0.8838834764831848*abar0[1]*favg[5]+0.3471825374147069*abar0[3]*favg[4]+0.3535533905932739*abar0[1]*favg[4]+0.6123724356957947*abar0[2]*favg[3]+0.6846531968814579*abar0[0]*favg[3]+0.6846531968814579*abar0[1]*favg[2]+0.3535533905932739*favg[1]*abar0[2]+0.3952847075210475*abar0[0]*favg[1]+0.3952847075210475*favg[0]*abar0[1]; 
  incr[8] = (-0.4330127018922193*fjump[11]*amax)-0.25*fjump[8]*amax+0.5303300858899107*abar0[3]*favg[14]+0.09642365197998377*abar0[3]*favg[13]+0.154303349962092*abar0[1]*favg[13]+0.4082482904638632*abar0[4]*favg[12]+0.4107919181288747*abar0[2]*favg[12]+0.1670106642806713*abar0[4]*favg[11]+0.1825741858350555*abar0[2]*favg[11]+0.3061862178478973*abar0[0]*favg[11]+0.2357022603955159*abar0[3]*favg[10]+0.3471825374147067*abar0[1]*favg[10]+0.4677071733467428*abar0[3]*favg[9]+0.09642365197998377*abar0[4]*favg[8]+0.105409255338946*abar0[2]*favg[8]+0.1767766952966369*abar0[0]*favg[8]+0.3450327796711773*abar0[4]*favg[7]+0.3471825374147069*abar0[2]*favg[7]+0.1825741858350555*abar0[3]*favg[6]+0.2689264371002386*abar0[1]*favg[6]+0.3952847075210476*abar0[3]*favg[5]+0.105409255338946*abar0[3]*favg[4]+0.1552647508520297*abar0[1]*favg[4]+0.2672612419124245*favg[3]*abar0[4]+0.154303349962092*favg[1]*abar0[4]+0.2689264371002386*abar0[2]*favg[3]+0.3061862178478972*favg[2]*abar0[3]+0.1767766952966369*favg[0]*abar0[3]+0.1552647508520297*favg[1]*abar0[2]; 
  incr[9] = 1.984313483298443*fjump[14]*amax+1.75*fjump[9]*amax+1.479019945774904*fjump[5]*amax+1.14564392373896*fjump[2]*amax+0.6614378277661477*fjump[0]*amax-1.403121520040228*abar0[0]*favg[14]-0.4677071733467427*abar0[4]*favg[13]-1.237436867076458*abar0[1]*favg[12]-0.8100925873009825*abar0[3]*favg[11]-1.045825033167595*abar0[2]*favg[10]-1.237436867076458*abar0[0]*favg[9]-0.4677071733467427*abar0[3]*favg[8]-1.045825033167594*abar0[1]*favg[7]-0.8100925873009828*abar0[2]*favg[6]-1.045825033167595*abar0[0]*favg[5]-0.4677071733467427*abar0[2]*favg[4]-0.8100925873009825*abar0[1]*favg[3]-0.8100925873009825*abar0[0]*favg[2]-0.4677071733467427*abar0[1]*favg[1]-0.4677071733467427*abar0[0]*favg[0]; 
  incr[10] = (-1.25*fjump[10]*amax)-0.9682458365518543*fjump[6]*amax-0.5590169943749475*fjump[4]*amax+1.185854122563142*abar0[2]*favg[14]+0.2295801237618662*abar0[4]*favg[13]+0.3388154635894693*abar0[2]*favg[13]+0.9185586535436919*abar0[3]*favg[12]+0.9354143466934854*abar0[1]*favg[12]+0.4082482904638632*abar0[3]*favg[11]+0.6013377943029545*abar0[1]*favg[11]+0.7576144084141582*abar0[4]*favg[10]+0.5646924393157823*abar0[2]*favg[10]+0.8838834764831844*abar0[0]*favg[10]+1.045825033167595*abar0[2]*favg[9]+0.2357022603955159*abar0[3]*favg[8]+0.3471825374147067*abar0[1]*favg[8]+0.7763237542601487*abar0[3]*favg[7]+0.790569415042095*abar0[1]*favg[7]+0.5868455973269638*abar0[4]*favg[6]+0.4374088826398534*abar0[2]*favg[6]+0.6846531968814578*abar0[0]*favg[6]+0.8838834764831844*abar0[2]*favg[5]+0.3388154635894693*abar0[4]*favg[4]+0.2525381361380528*abar0[2]*favg[4]+0.3952847075210474*abar0[0]*favg[4]+0.601337794302955*abar0[3]*favg[3]+0.6123724356957946*abar0[1]*favg[3]+0.3471825374147067*favg[1]*abar0[3]+0.6846531968814576*abar0[2]*favg[2]+0.3952847075210474*favg[0]*abar0[2]+0.3535533905932738*abar0[1]*favg[1]; 
  incr[11] = 0.75*fjump[11]*amax+0.4330127018922193*fjump[8]*amax-0.9185586535436919*abar0[3]*favg[14]-0.1670106642806713*abar0[3]*favg[13]-0.2672612419124245*abar0[1]*favg[13]-0.7071067811865478*abar0[4]*favg[12]-0.7115124735378855*abar0[2]*favg[12]-0.2892709559399513*abar0[4]*favg[11]-0.3162277660168381*abar0[2]*favg[11]-0.5303300858899107*abar0[0]*favg[11]-0.4082482904638632*abar0[3]*favg[10]-0.6013377943029545*abar0[1]*favg[10]-0.8100925873009823*abar0[3]*favg[9]-0.1670106642806713*abar0[4]*favg[8]-0.1825741858350555*abar0[2]*favg[8]-0.3061862178478973*abar0[0]*favg[8]-0.597614304667197*abar0[4]*favg[7]-0.601337794302955*abar0[2]*favg[7]-0.3162277660168381*abar0[3]*favg[6]-0.465794252556089*abar0[1]*favg[6]-0.6846531968814579*abar0[3]*favg[5]-0.1825741858350555*abar0[3]*favg[4]-0.2689264371002386*abar0[1]*favg[4]-0.4629100498862759*favg[3]*abar0[4]-0.2672612419124245*favg[1]*abar0[4]-0.4657942525560891*abar0[2]*favg[3]-0.5303300858899105*favg[2]*abar0[3]-0.3061862178478973*favg[0]*abar0[3]-0.2689264371002386*favg[1]*abar0[2]; 
  incr[12] = 1.75*fjump[12]*amax+1.479019945774904*fjump[7]*amax+1.14564392373896*fjump[3]*amax+0.6614378277661476*fjump[1]*amax-1.403121520040228*abar0[1]*favg[14]-0.4082482904638632*abar0[3]*favg[13]-1.106797181058933*abar0[2]*favg[12]-1.237436867076458*abar0[0]*favg[12]-0.7071067811865478*abar0[4]*favg[11]-0.7115124735378856*abar0[2]*favg[11]-0.9185586535436919*abar0[3]*favg[10]-0.9354143466934856*abar0[1]*favg[10]-1.237436867076458*abar0[1]*favg[9]-0.4082482904638632*abar0[4]*favg[8]-0.4107919181288747*abar0[2]*favg[8]-0.9354143466934858*abar0[2]*favg[7]-1.045825033167595*abar0[0]*favg[7]-0.7115124735378855*abar0[3]*favg[6]-0.7245688373094722*abar0[1]*favg[6]-1.045825033167595*abar0[1]*favg[5]-0.4107919181288747*abar0[3]*favg[4]-0.4183300132670379*abar0[1]*favg[4]-0.7245688373094721*abar0[2]*favg[3]-0.8100925873009825*abar0[0]*favg[3]-0.8100925873009825*abar0[1]*favg[2]-0.4183300132670379*favg[1]*abar0[2]-0.4677071733467427*abar0[0]*favg[1]-0.4677071733467427*favg[0]*abar0[1]; 
  incr[13] = (-0.25*fjump[13]*amax)+0.5303300858899107*abar0[4]*favg[14]+0.08582764626789764*abar0[4]*favg[13]+0.1026713526028695*abar0[2]*favg[13]+0.1767766952966369*abar0[0]*favg[13]+0.4082482904638632*abar0[3]*favg[12]+0.1670106642806713*abar0[3]*favg[11]+0.2672612419124245*abar0[1]*favg[11]+0.2295801237618662*abar0[4]*favg[10]+0.3388154635894693*abar0[2]*favg[10]+0.4677071733467428*abar0[4]*favg[9]+0.09642365197998377*abar0[3]*favg[8]+0.154303349962092*abar0[1]*favg[8]+0.3450327796711773*abar0[3]*favg[7]+0.1778319991899891*abar0[4]*favg[6]+0.262445329583912*abar0[2]*favg[6]+0.3952847075210476*abar0[4]*favg[5]+0.1026713526028695*abar0[4]*favg[4]+0.1515228816828317*abar0[2]*favg[4]+0.3061862178478972*favg[2]*abar0[4]+0.1767766952966369*favg[0]*abar0[4]+0.2672612419124245*abar0[3]*favg[3]+0.154303349962092*favg[1]*abar0[3]; 
  incr[14] = (-2.25*fjump[14]*amax)-1.984313483298443*fjump[9]*amax-1.677050983124842*fjump[5]*amax-1.299038105676658*fjump[2]*amax-0.75*fjump[0]*amax+1.590990257669732*abar0[0]*favg[14]+0.5303300858899107*abar0[4]*favg[13]+1.403121520040228*abar0[1]*favg[12]+0.9185586535436919*abar0[3]*favg[11]+1.185854122563142*abar0[2]*favg[10]+1.403121520040228*abar0[0]*favg[9]+0.5303300858899107*abar0[3]*favg[8]+1.185854122563142*abar0[1]*favg[7]+0.9185586535436919*abar0[2]*favg[6]+1.185854122563142*abar0[0]*favg[5]+0.5303300858899107*abar0[2]*favg[4]+0.9185586535436917*abar0[1]*favg[3]+0.9185586535436917*abar0[0]*favg[2]+0.5303300858899107*abar0[1]*favg[1]+0.5303300858899107*abar0[0]*favg[0]; 

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
  outr[10] += incr[10]*dv10; 
  outr[11] += incr[11]*dv10; 
  outr[12] += incr[12]*dv10; 
  outr[13] += incr[13]*dv10; 
  outr[14] += incr[14]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += incr[2]*dv10; 
  outl[3] += incr[3]*dv10; 
  outl[4] += -1.0*incr[4]*dv10; 
  outl[5] += -1.0*incr[5]*dv10; 
  outl[6] += incr[6]*dv10; 
  outl[7] += -1.0*incr[7]*dv10; 
  outl[8] += -1.0*incr[8]*dv10; 
  outl[9] += incr[9]*dv10; 
  outl[10] += -1.0*incr[10]*dv10; 
  outl[11] += incr[11]*dv10; 
  outl[12] += incr[12]*dv10; 
  outl[13] += -1.0*incr[13]*dv10; 
  outl[14] += -1.0*incr[14]*dv10; 
return std::abs(amid); 
} 
