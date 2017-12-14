void VlasovSurfStream1x2vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  const unsigned int X=0, VX=1; 
  double dvx = dxv[VX]; 
  double wx = w[VX]; 

  if (wx>0) { 
  outr[0] += 0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.1443375672974065*fl[2]*dvx; 
  outr[1] += (-1.5*fl[1]*wx)-0.8660254037844386*fl[0]*wx-0.25*fl[2]*dvx; 
  outr[2] += 0.5*fl[2]*wx+0.25*fl[1]*dvx+0.1443375672974065*fl[0]*dvx; 
  outr[3] += 0.5*fl[3]*wx; 

  outl[0] += 0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.1443375672974065*fl[2]*dvx; 
  outl[1] += 1.5*fl[1]*wx+0.8660254037844386*fl[0]*wx+0.25*fl[2]*dvx; 
  outl[2] += 0.5*fl[2]*wx+0.25*fl[1]*dvx+0.1443375672974065*fl[0]*dvx; 
  outl[3] += 0.5*fl[3]*wx; 
  } else { 
  outr[0] += (-0.8660254037844386*fr[1]*wx)+0.5*fr[0]*wx+0.1443375672974065*fr[2]*dvx; 
  outr[1] += 1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx-0.25*fr[2]*dvx; 
  outr[2] += 0.5*fr[2]*wx-0.25*fr[1]*dvx+0.1443375672974065*fr[0]*dvx; 
  outr[3] += 0.5*fr[3]*wx; 

  outl[0] += (-0.8660254037844386*fr[1]*wx)+0.5*fr[0]*wx+0.1443375672974065*fr[2]*dvx; 
  outl[1] += (-1.5*fr[1]*wx)+0.8660254037844386*fr[0]*wx+0.25*fr[2]*dvx; 
  outl[2] += 0.5*fr[2]*wx-0.25*fr[1]*dvx+0.1443375672974065*fr[0]*dvx; 
  outl[3] += 0.5*fr[3]*wx; 
  } 
} 
void VlasovSurfStream1x2vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  const unsigned int X=0, VX=1; 
  double dvx = dxv[VX]; 
  double wx = w[VX]; 

  if (wx>0) { 
  outr[0] += 1.118033988749895*fl[7]*wx+0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.25*fl[4]*dvx+0.1443375672974065*fl[2]*dvx; 
  outr[1] += (-1.936491673103709*fl[7]*wx)-1.5*fl[1]*wx-0.8660254037844386*fl[0]*wx-0.4330127018922193*fl[4]*dvx-0.25*fl[2]*dvx; 
  outr[2] += 0.8660254037844386*fl[4]*wx+0.5*fl[2]*wx+0.1290994448735806*fl[8]*dvx+0.3227486121839515*fl[7]*dvx+0.25*fl[1]*dvx+0.1443375672974065*fl[0]*dvx; 
  outr[3] += 0.8660254037844386*fl[5]*wx+0.5*fl[3]*wx+0.1443375672974065*fl[6]*dvx; 
  outr[4] += (-1.5*fl[4]*wx)-0.8660254037844386*fl[2]*wx-0.223606797749979*fl[8]*dvx-0.5590169943749475*fl[7]*dvx-0.4330127018922193*fl[1]*dvx-0.25*fl[0]*dvx; 
  outr[5] += (-1.5*fl[5]*wx)-0.8660254037844386*fl[3]*wx-0.25*fl[6]*dvx; 
  outr[6] += 0.5*fl[6]*wx+0.25*fl[5]*dvx+0.1443375672974065*fl[3]*dvx; 
  outr[7] += 2.5*fl[7]*wx+1.936491673103709*fl[1]*wx+1.118033988749895*fl[0]*wx+0.5590169943749475*fl[4]*dvx+0.3227486121839515*fl[2]*dvx; 
  outr[8] += 0.5*fl[8]*wx+0.223606797749979*fl[4]*dvx+0.1290994448735806*fl[2]*dvx; 
  outr[9] += 0.5*fl[9]*wx; 

  outl[0] += 1.118033988749895*fl[7]*wx+0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.25*fl[4]*dvx+0.1443375672974065*fl[2]*dvx; 
  outl[1] += 1.936491673103709*fl[7]*wx+1.5*fl[1]*wx+0.8660254037844386*fl[0]*wx+0.4330127018922193*fl[4]*dvx+0.25*fl[2]*dvx; 
  outl[2] += 0.8660254037844386*fl[4]*wx+0.5*fl[2]*wx+0.1290994448735806*fl[8]*dvx+0.3227486121839515*fl[7]*dvx+0.25*fl[1]*dvx+0.1443375672974065*fl[0]*dvx; 
  outl[3] += 0.8660254037844386*fl[5]*wx+0.5*fl[3]*wx+0.1443375672974065*fl[6]*dvx; 
  outl[4] += 1.5*fl[4]*wx+0.8660254037844386*fl[2]*wx+0.223606797749979*fl[8]*dvx+0.5590169943749475*fl[7]*dvx+0.4330127018922193*fl[1]*dvx+0.25*fl[0]*dvx; 
  outl[5] += 1.5*fl[5]*wx+0.8660254037844386*fl[3]*wx+0.25*fl[6]*dvx; 
  outl[6] += 0.5*fl[6]*wx+0.25*fl[5]*dvx+0.1443375672974065*fl[3]*dvx; 
  outl[7] += 2.5*fl[7]*wx+1.936491673103709*fl[1]*wx+1.118033988749895*fl[0]*wx+0.5590169943749475*fl[4]*dvx+0.3227486121839515*fl[2]*dvx; 
  outl[8] += 0.5*fl[8]*wx+0.223606797749979*fl[4]*dvx+0.1290994448735806*fl[2]*dvx; 
  outl[9] += 0.5*fl[9]*wx; 
  } else { 
  outr[0] += 1.118033988749895*fr[7]*wx-0.8660254037844386*fr[1]*wx+0.5*fr[0]*wx-0.25*fr[4]*dvx+0.1443375672974065*fr[2]*dvx; 
  outr[1] += (-1.936491673103709*fr[7]*wx)+1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx+0.4330127018922193*fr[4]*dvx-0.25*fr[2]*dvx; 
  outr[2] += (-0.8660254037844386*fr[4]*wx)+0.5*fr[2]*wx+0.1290994448735806*fr[8]*dvx+0.3227486121839515*fr[7]*dvx-0.25*fr[1]*dvx+0.1443375672974065*fr[0]*dvx; 
  outr[3] += (-0.8660254037844386*fr[5]*wx)+0.5*fr[3]*wx+0.1443375672974065*fr[6]*dvx; 
  outr[4] += 1.5*fr[4]*wx-0.8660254037844386*fr[2]*wx-0.223606797749979*fr[8]*dvx-0.5590169943749475*fr[7]*dvx+0.4330127018922193*fr[1]*dvx-0.25*fr[0]*dvx; 
  outr[5] += 1.5*fr[5]*wx-0.8660254037844386*fr[3]*wx-0.25*fr[6]*dvx; 
  outr[6] += 0.5*fr[6]*wx-0.25*fr[5]*dvx+0.1443375672974065*fr[3]*dvx; 
  outr[7] += 2.5*fr[7]*wx-1.936491673103709*fr[1]*wx+1.118033988749895*fr[0]*wx-0.5590169943749475*fr[4]*dvx+0.3227486121839515*fr[2]*dvx; 
  outr[8] += 0.5*fr[8]*wx-0.223606797749979*fr[4]*dvx+0.1290994448735806*fr[2]*dvx; 
  outr[9] += 0.5*fr[9]*wx; 

  outl[0] += 1.118033988749895*fr[7]*wx-0.8660254037844386*fr[1]*wx+0.5*fr[0]*wx-0.25*fr[4]*dvx+0.1443375672974065*fr[2]*dvx; 
  outl[1] += 1.936491673103709*fr[7]*wx-1.5*fr[1]*wx+0.8660254037844386*fr[0]*wx-0.4330127018922193*fr[4]*dvx+0.25*fr[2]*dvx; 
  outl[2] += (-0.8660254037844386*fr[4]*wx)+0.5*fr[2]*wx+0.1290994448735806*fr[8]*dvx+0.3227486121839515*fr[7]*dvx-0.25*fr[1]*dvx+0.1443375672974065*fr[0]*dvx; 
  outl[3] += (-0.8660254037844386*fr[5]*wx)+0.5*fr[3]*wx+0.1443375672974065*fr[6]*dvx; 
  outl[4] += (-1.5*fr[4]*wx)+0.8660254037844386*fr[2]*wx+0.223606797749979*fr[8]*dvx+0.5590169943749475*fr[7]*dvx-0.4330127018922193*fr[1]*dvx+0.25*fr[0]*dvx; 
  outl[5] += (-1.5*fr[5]*wx)+0.8660254037844386*fr[3]*wx+0.25*fr[6]*dvx; 
  outl[6] += 0.5*fr[6]*wx-0.25*fr[5]*dvx+0.1443375672974065*fr[3]*dvx; 
  outl[7] += 2.5*fr[7]*wx-1.936491673103709*fr[1]*wx+1.118033988749895*fr[0]*wx-0.5590169943749475*fr[4]*dvx+0.3227486121839515*fr[2]*dvx; 
  outl[8] += 0.5*fr[8]*wx-0.223606797749979*fr[4]*dvx+0.1290994448735806*fr[2]*dvx; 
  outl[9] += 0.5*fr[9]*wx; 
  } 
} 
