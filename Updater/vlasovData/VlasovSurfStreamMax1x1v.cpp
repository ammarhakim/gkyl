#include <VlasovModDecl.h> 
void VlasovSurfStream1x1vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[1]/dxv[0]; 
  double wx = w[1]*2/dxv[0]; 

  if (wx>0) { 
  outr[0] += 0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.2886751345948129*fl[2]*dvx; 
  outr[1] += (-1.5*fl[1]*wx)-0.8660254037844386*fl[0]*wx-0.5*fl[2]*dvx; 
  outr[2] += 0.5*fl[2]*wx+0.5*fl[1]*dvx+0.2886751345948129*fl[0]*dvx; 

  outl[0] += (-0.8660254037844386*fl[1]*wx)-0.5*fl[0]*wx-0.2886751345948129*fl[2]*dvx; 
  outl[1] += (-1.5*fl[1]*wx)-0.8660254037844386*fl[0]*wx-0.5*fl[2]*dvx; 
  outl[2] += (-0.5*fl[2]*wx)-0.5*fl[1]*dvx-0.2886751345948129*fl[0]*dvx; 
  } else { 
  outr[0] += (-0.8660254037844386*fr[1]*wx)+0.5*fr[0]*wx+0.2886751345948129*fr[2]*dvx; 
  outr[1] += 1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx-0.5*fr[2]*dvx; 
  outr[2] += 0.5*fr[2]*wx-0.5*fr[1]*dvx+0.2886751345948129*fr[0]*dvx; 

  outl[0] += 0.8660254037844386*fr[1]*wx-0.5*fr[0]*wx-0.2886751345948129*fr[2]*dvx; 
  outl[1] += 1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx-0.5*fr[2]*dvx; 
  outl[2] += (-0.5*fr[2]*wx)+0.5*fr[1]*dvx-0.2886751345948129*fr[0]*dvx; 
  } 
} 
void VlasovSurfStream1x1vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[1]/dxv[0]; 
  double wx = w[1]*2/dxv[0]; 

  if (wx>0) { 
  outr[0] += 1.118033988749895*fl[4]*wx+0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.5*fl[3]*dvx+0.2886751345948129*fl[2]*dvx; 
  outr[1] += (-1.936491673103709*fl[4]*wx)-1.5*fl[1]*wx-0.8660254037844386*fl[0]*wx-0.8660254037844386*fl[3]*dvx-0.5*fl[2]*dvx; 
  outr[2] += 0.8660254037844386*fl[3]*wx+0.5*fl[2]*wx+0.2581988897471612*fl[5]*dvx+0.6454972243679029*fl[4]*dvx+0.5*fl[1]*dvx+0.2886751345948129*fl[0]*dvx; 
  outr[3] += (-1.5*fl[3]*wx)-0.8660254037844386*fl[2]*wx-0.4472135954999579*fl[5]*dvx-1.118033988749895*fl[4]*dvx-0.8660254037844386*fl[1]*dvx-0.5*fl[0]*dvx; 
  outr[4] += 2.5*fl[4]*wx+1.936491673103709*fl[1]*wx+1.118033988749895*fl[0]*wx+1.118033988749895*fl[3]*dvx+0.6454972243679029*fl[2]*dvx; 
  outr[5] += 0.5*fl[5]*wx+0.4472135954999579*fl[3]*dvx+0.2581988897471612*fl[2]*dvx; 

  outl[0] += (-1.118033988749895*fl[4]*wx)-0.8660254037844386*fl[1]*wx-0.5*fl[0]*wx-0.5*fl[3]*dvx-0.2886751345948129*fl[2]*dvx; 
  outl[1] += (-1.936491673103709*fl[4]*wx)-1.5*fl[1]*wx-0.8660254037844386*fl[0]*wx-0.8660254037844386*fl[3]*dvx-0.5*fl[2]*dvx; 
  outl[2] += (-0.8660254037844386*fl[3]*wx)-0.5*fl[2]*wx-0.2581988897471612*fl[5]*dvx-0.6454972243679029*fl[4]*dvx-0.5*fl[1]*dvx-0.2886751345948129*fl[0]*dvx; 
  outl[3] += (-1.5*fl[3]*wx)-0.8660254037844386*fl[2]*wx-0.4472135954999579*fl[5]*dvx-1.118033988749895*fl[4]*dvx-0.8660254037844386*fl[1]*dvx-0.5*fl[0]*dvx; 
  outl[4] += (-2.5*fl[4]*wx)-1.936491673103709*fl[1]*wx-1.118033988749895*fl[0]*wx-1.118033988749895*fl[3]*dvx-0.6454972243679029*fl[2]*dvx; 
  outl[5] += (-0.5*fl[5]*wx)-0.4472135954999579*fl[3]*dvx-0.2581988897471612*fl[2]*dvx; 
  } else { 
  outr[0] += 1.118033988749895*fr[4]*wx-0.8660254037844386*fr[1]*wx+0.5*fr[0]*wx-0.5*fr[3]*dvx+0.2886751345948129*fr[2]*dvx; 
  outr[1] += (-1.936491673103709*fr[4]*wx)+1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx+0.8660254037844386*fr[3]*dvx-0.5*fr[2]*dvx; 
  outr[2] += (-0.8660254037844386*fr[3]*wx)+0.5*fr[2]*wx+0.2581988897471612*fr[5]*dvx+0.6454972243679029*fr[4]*dvx-0.5*fr[1]*dvx+0.2886751345948129*fr[0]*dvx; 
  outr[3] += 1.5*fr[3]*wx-0.8660254037844386*fr[2]*wx-0.4472135954999579*fr[5]*dvx-1.118033988749895*fr[4]*dvx+0.8660254037844386*fr[1]*dvx-0.5*fr[0]*dvx; 
  outr[4] += 2.5*fr[4]*wx-1.936491673103709*fr[1]*wx+1.118033988749895*fr[0]*wx-1.118033988749895*fr[3]*dvx+0.6454972243679029*fr[2]*dvx; 
  outr[5] += 0.5*fr[5]*wx-0.4472135954999579*fr[3]*dvx+0.2581988897471612*fr[2]*dvx; 

  outl[0] += (-1.118033988749895*fr[4]*wx)+0.8660254037844386*fr[1]*wx-0.5*fr[0]*wx+0.5*fr[3]*dvx-0.2886751345948129*fr[2]*dvx; 
  outl[1] += (-1.936491673103709*fr[4]*wx)+1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx+0.8660254037844386*fr[3]*dvx-0.5*fr[2]*dvx; 
  outl[2] += 0.8660254037844386*fr[3]*wx-0.5*fr[2]*wx-0.2581988897471612*fr[5]*dvx-0.6454972243679029*fr[4]*dvx+0.5*fr[1]*dvx-0.2886751345948129*fr[0]*dvx; 
  outl[3] += 1.5*fr[3]*wx-0.8660254037844386*fr[2]*wx-0.4472135954999579*fr[5]*dvx-1.118033988749895*fr[4]*dvx+0.8660254037844386*fr[1]*dvx-0.5*fr[0]*dvx; 
  outl[4] += (-2.5*fr[4]*wx)+1.936491673103709*fr[1]*wx-1.118033988749895*fr[0]*wx+1.118033988749895*fr[3]*dvx-0.6454972243679029*fr[2]*dvx; 
  outl[5] += (-0.5*fr[5]*wx)+0.4472135954999579*fr[3]*dvx-0.2581988897471612*fr[2]*dvx; 
  } 
} 
