#include <VlasovModDecl.h> 
void VlasovSurfElc1x2vMax_VX_P1(const double *w, const double *dxv, const double *Ein, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. E: electric field, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dv1 = 2/dxv[1]; 
  const double *E = &Ein[0]; 
  double incr[4]; 

  if (E[0]>0) { 
  incr[0] = 0.6123724356957944*E[0]*fl[2]*dv1+0.3535533905932737*E[1]*fl[1]*dv1+0.3535533905932737*E[0]*fl[0]*dv1; 
  incr[1] = 0.6123724356957944*E[1]*fl[2]*dv1+0.3535533905932737*E[0]*fl[1]*dv1+0.3535533905932737*fl[0]*E[1]*dv1; 
  incr[2] = (-1.060660171779821*E[0]*fl[2]*dv1)-0.6123724356957944*E[1]*fl[1]*dv1-0.6123724356957944*E[0]*fl[0]*dv1; 
  incr[3] = 0.3535533905932737*E[0]*fl[3]*dv1; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  } else { 
  incr[0] = (-0.6123724356957944*E[0]*fr[2]*dv1)+0.3535533905932737*E[1]*fr[1]*dv1+0.3535533905932737*E[0]*fr[0]*dv1; 
  incr[1] = (-0.6123724356957944*E[1]*fr[2]*dv1)+0.3535533905932737*E[0]*fr[1]*dv1+0.3535533905932737*fr[0]*E[1]*dv1; 
  incr[2] = 1.060660171779821*E[0]*fr[2]*dv1-0.6123724356957944*E[1]*fr[1]*dv1-0.6123724356957944*E[0]*fr[0]*dv1; 
  incr[3] = 0.3535533905932737*E[0]*fr[3]*dv1; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  } 
} 
void VlasovSurfElc1x2vMax_VX_P2(const double *w, const double *dxv, const double *Ein, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. E: electric field, fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
  double dv1 = 2/dxv[1]; 
  const double *E = &Ein[0]; 
  double incr[10]; 

  if (E[0]>0) { 
  incr[0] = 0.7905694150420947*E[0]*fl[8]*dv1+0.3535533905932737*E[2]*fl[7]*dv1+0.6123724356957944*E[1]*fl[4]*dv1+0.6123724356957944*E[0]*fl[2]*dv1+0.3535533905932737*E[1]*fl[1]*dv1+0.3535533905932737*E[0]*fl[0]*dv1; 
  incr[1] = 0.7905694150420947*E[1]*fl[8]*dv1+0.3162277660168379*E[1]*fl[7]*dv1+0.5477225575051661*E[2]*fl[4]*dv1+0.6123724356957944*E[0]*fl[4]*dv1+0.6123724356957944*E[1]*fl[2]*dv1+0.3162277660168379*fl[1]*E[2]*dv1+0.3535533905932737*E[0]*fl[1]*dv1+0.3535533905932737*fl[0]*E[1]*dv1; 
  incr[2] = (-1.369306393762915*E[0]*fl[8]*dv1)-0.6123724356957944*E[2]*fl[7]*dv1-1.060660171779821*E[1]*fl[4]*dv1-1.060660171779821*E[0]*fl[2]*dv1-0.6123724356957944*E[1]*fl[1]*dv1-0.6123724356957944*E[0]*fl[0]*dv1; 
  incr[3] = 0.6123724356957944*E[0]*fl[6]*dv1+0.3535533905932737*E[1]*fl[5]*dv1+0.3535533905932737*E[0]*fl[3]*dv1; 
  incr[4] = (-1.369306393762915*E[1]*fl[8]*dv1)-0.5477225575051661*E[1]*fl[7]*dv1-0.9486832980505137*E[2]*fl[4]*dv1-1.060660171779821*E[0]*fl[4]*dv1-1.060660171779821*E[1]*fl[2]*dv1-0.5477225575051661*fl[1]*E[2]*dv1-0.6123724356957944*E[0]*fl[1]*dv1-0.6123724356957944*fl[0]*E[1]*dv1; 
  incr[5] = 0.6123724356957944*E[1]*fl[6]*dv1+0.3162277660168379*E[2]*fl[5]*dv1+0.3535533905932737*E[0]*fl[5]*dv1+0.3535533905932737*E[1]*fl[3]*dv1; 
  incr[6] = (-1.060660171779821*E[0]*fl[6]*dv1)-0.6123724356957944*E[1]*fl[5]*dv1-0.6123724356957944*E[0]*fl[3]*dv1; 
  incr[7] = 0.7905694150420947*E[2]*fl[8]*dv1+0.2258769757263128*E[2]*fl[7]*dv1+0.3535533905932737*E[0]*fl[7]*dv1+0.5477225575051661*E[1]*fl[4]*dv1+0.6123724356957944*E[2]*fl[2]*dv1+0.3535533905932737*fl[0]*E[2]*dv1+0.3162277660168379*E[1]*fl[1]*dv1; 
  incr[8] = 1.767766952966368*E[0]*fl[8]*dv1+0.7905694150420947*E[2]*fl[7]*dv1+1.369306393762915*E[1]*fl[4]*dv1+1.369306393762915*E[0]*fl[2]*dv1+0.7905694150420947*E[1]*fl[1]*dv1+0.7905694150420947*E[0]*fl[0]*dv1; 
  incr[9] = 0.3535533905932737*E[0]*fl[9]*dv1; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  } else { 
  incr[0] = 0.7905694150420947*E[0]*fr[8]*dv1+0.3535533905932737*E[2]*fr[7]*dv1-0.6123724356957944*E[1]*fr[4]*dv1-0.6123724356957944*E[0]*fr[2]*dv1+0.3535533905932737*E[1]*fr[1]*dv1+0.3535533905932737*E[0]*fr[0]*dv1; 
  incr[1] = 0.7905694150420947*E[1]*fr[8]*dv1+0.3162277660168379*E[1]*fr[7]*dv1-0.5477225575051661*E[2]*fr[4]*dv1-0.6123724356957944*E[0]*fr[4]*dv1-0.6123724356957944*E[1]*fr[2]*dv1+0.3162277660168379*fr[1]*E[2]*dv1+0.3535533905932737*E[0]*fr[1]*dv1+0.3535533905932737*fr[0]*E[1]*dv1; 
  incr[2] = (-1.369306393762915*E[0]*fr[8]*dv1)-0.6123724356957944*E[2]*fr[7]*dv1+1.060660171779821*E[1]*fr[4]*dv1+1.060660171779821*E[0]*fr[2]*dv1-0.6123724356957944*E[1]*fr[1]*dv1-0.6123724356957944*E[0]*fr[0]*dv1; 
  incr[3] = (-0.6123724356957944*E[0]*fr[6]*dv1)+0.3535533905932737*E[1]*fr[5]*dv1+0.3535533905932737*E[0]*fr[3]*dv1; 
  incr[4] = (-1.369306393762915*E[1]*fr[8]*dv1)-0.5477225575051661*E[1]*fr[7]*dv1+0.9486832980505137*E[2]*fr[4]*dv1+1.060660171779821*E[0]*fr[4]*dv1+1.060660171779821*E[1]*fr[2]*dv1-0.5477225575051661*fr[1]*E[2]*dv1-0.6123724356957944*E[0]*fr[1]*dv1-0.6123724356957944*fr[0]*E[1]*dv1; 
  incr[5] = (-0.6123724356957944*E[1]*fr[6]*dv1)+0.3162277660168379*E[2]*fr[5]*dv1+0.3535533905932737*E[0]*fr[5]*dv1+0.3535533905932737*E[1]*fr[3]*dv1; 
  incr[6] = 1.060660171779821*E[0]*fr[6]*dv1-0.6123724356957944*E[1]*fr[5]*dv1-0.6123724356957944*E[0]*fr[3]*dv1; 
  incr[7] = 0.7905694150420947*E[2]*fr[8]*dv1+0.2258769757263128*E[2]*fr[7]*dv1+0.3535533905932737*E[0]*fr[7]*dv1-0.5477225575051661*E[1]*fr[4]*dv1-0.6123724356957944*E[2]*fr[2]*dv1+0.3535533905932737*fr[0]*E[2]*dv1+0.3162277660168379*E[1]*fr[1]*dv1; 
  incr[8] = 1.767766952966368*E[0]*fr[8]*dv1+0.7905694150420947*E[2]*fr[7]*dv1-1.369306393762915*E[1]*fr[4]*dv1-1.369306393762915*E[0]*fr[2]*dv1+0.7905694150420947*E[1]*fr[1]*dv1+0.7905694150420947*E[0]*fr[0]*dv1; 
  incr[9] = 0.3535533905932737*E[0]*fr[9]*dv1; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  } 
} 
