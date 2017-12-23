#include <VlasovModDecl.h> 
void VlasovSurfStream1x1vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[1]/dxv[0]; 
  double wx = w[1]*2/dxv[0]; 

  double incr[4]; 

  if (wx>0) { 
  incr[0] = 0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.5*fl[3]*dvx+0.2886751345948129*fl[2]*dvx; 
  incr[1] = (-1.5*fl[1]*wx)-0.8660254037844386*fl[0]*wx-0.8660254037844386*fl[3]*dvx-0.5*fl[2]*dvx; 
  incr[2] = 0.8660254037844386*fl[3]*wx+0.5*fl[2]*wx+0.5*fl[1]*dvx+0.2886751345948129*fl[0]*dvx; 
  incr[3] = (-1.5*fl[3]*wx)-0.8660254037844386*fl[2]*wx-0.8660254037844386*fl[1]*dvx-0.5*fl[0]*dvx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = (-0.8660254037844386*fr[1]*wx)+0.5*fr[0]*wx-0.5*fr[3]*dvx+0.2886751345948129*fr[2]*dvx; 
  incr[1] = 1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx+0.8660254037844386*fr[3]*dvx-0.5*fr[2]*dvx; 
  incr[2] = (-0.8660254037844386*fr[3]*wx)+0.5*fr[2]*wx-0.5*fr[1]*dvx+0.2886751345948129*fr[0]*dvx; 
  incr[3] = 1.5*fr[3]*wx-0.8660254037844386*fr[2]*wx+0.8660254037844386*fr[1]*dvx-0.5*fr[0]*dvx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
} 
void VlasovSurfStream1x1vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[1]/dxv[0]; 
  double wx = w[1]*2/dxv[0]; 

  double incr[8]; 

  if (wx>0) { 
  incr[0] = 1.118033988749895*fl[4]*wx+0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.6454972243679029*fl[6]*dvx+0.5*fl[3]*dvx+0.2886751345948129*fl[2]*dvx; 
  incr[1] = (-1.936491673103709*fl[4]*wx)-1.5*fl[1]*wx-0.8660254037844386*fl[0]*wx-1.118033988749895*fl[6]*dvx-0.8660254037844386*fl[3]*dvx-0.5*fl[2]*dvx; 
  incr[2] = 1.118033988749895*fl[6]*wx+0.8660254037844386*fl[3]*wx+0.5*fl[2]*wx+0.4472135954999579*fl[7]*dvx+0.2581988897471612*fl[5]*dvx+0.6454972243679029*fl[4]*dvx+0.5*fl[1]*dvx+0.2886751345948129*fl[0]*dvx; 
  incr[3] = (-1.936491673103709*fl[6]*wx)-1.5*fl[3]*wx-0.8660254037844386*fl[2]*wx-0.7745966692414833*fl[7]*dvx-0.4472135954999579*fl[5]*dvx-1.118033988749895*fl[4]*dvx-0.8660254037844386*fl[1]*dvx-0.5*fl[0]*dvx; 
  incr[4] = 2.5*fl[4]*wx+1.936491673103709*fl[1]*wx+1.118033988749895*fl[0]*wx+1.443375672974065*fl[6]*dvx+1.118033988749895*fl[3]*dvx+0.6454972243679029*fl[2]*dvx; 
  incr[5] = 0.8660254037844386*fl[7]*wx+0.5*fl[5]*wx+0.5773502691896258*fl[6]*dvx+0.4472135954999579*fl[3]*dvx+0.2581988897471612*fl[2]*dvx; 
  incr[6] = 2.5*fl[6]*wx+1.936491673103709*fl[3]*wx+1.118033988749895*fl[2]*wx+fl[7]*dvx+0.5773502691896258*fl[5]*dvx+1.443375672974065*fl[4]*dvx+1.118033988749895*fl[1]*dvx+0.6454972243679029*fl[0]*dvx; 
  incr[7] = (-1.5*fl[7]*wx)-0.8660254037844386*fl[5]*wx-1.0*fl[6]*dvx-0.7745966692414833*fl[3]*dvx-0.4472135954999579*fl[2]*dvx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = 1.118033988749895*fr[4]*wx-0.8660254037844386*fr[1]*wx+0.5*fr[0]*wx+0.6454972243679029*fr[6]*dvx-0.5*fr[3]*dvx+0.2886751345948129*fr[2]*dvx; 
  incr[1] = (-1.936491673103709*fr[4]*wx)+1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx-1.118033988749895*fr[6]*dvx+0.8660254037844386*fr[3]*dvx-0.5*fr[2]*dvx; 
  incr[2] = 1.118033988749895*fr[6]*wx-0.8660254037844386*fr[3]*wx+0.5*fr[2]*wx-0.4472135954999579*fr[7]*dvx+0.2581988897471612*fr[5]*dvx+0.6454972243679029*fr[4]*dvx-0.5*fr[1]*dvx+0.2886751345948129*fr[0]*dvx; 
  incr[3] = (-1.936491673103709*fr[6]*wx)+1.5*fr[3]*wx-0.8660254037844386*fr[2]*wx+0.7745966692414833*fr[7]*dvx-0.4472135954999579*fr[5]*dvx-1.118033988749895*fr[4]*dvx+0.8660254037844386*fr[1]*dvx-0.5*fr[0]*dvx; 
  incr[4] = 2.5*fr[4]*wx-1.936491673103709*fr[1]*wx+1.118033988749895*fr[0]*wx+1.443375672974065*fr[6]*dvx-1.118033988749895*fr[3]*dvx+0.6454972243679029*fr[2]*dvx; 
  incr[5] = (-0.8660254037844386*fr[7]*wx)+0.5*fr[5]*wx+0.5773502691896258*fr[6]*dvx-0.4472135954999579*fr[3]*dvx+0.2581988897471612*fr[2]*dvx; 
  incr[6] = 2.5*fr[6]*wx-1.936491673103709*fr[3]*wx+1.118033988749895*fr[2]*wx-1.0*fr[7]*dvx+0.5773502691896258*fr[5]*dvx+1.443375672974065*fr[4]*dvx-1.118033988749895*fr[1]*dvx+0.6454972243679029*fr[0]*dvx; 
  incr[7] = 1.5*fr[7]*wx-0.8660254037844386*fr[5]*wx-1.0*fr[6]*dvx+0.7745966692414833*fr[3]*dvx-0.4472135954999579*fr[2]*dvx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } 
} 
