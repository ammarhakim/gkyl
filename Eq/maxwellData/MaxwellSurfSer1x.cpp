#include <MaxwellModDecl.h> 
double MaxwellSurf1xSer_X_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[2]; 
  const double *ezl = &ql[4]; 
  const double *bxl = &ql[6]; 
  const double *byl = &ql[8]; 
  const double *bzl = &ql[10]; 
  const double *phl = &ql[12]; 
  const double *psl = &ql[14]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[2]; 
  double *outEzl = &outl[4]; 
  double *outBxl = &outl[6]; 
  double *outByl = &outl[8]; 
  double *outBzl = &outl[10]; 
  double *outPhl = &outl[12]; 
  double *outPsl = &outl[14]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[2]; 
  const double *ezr = &qr[4]; 
  const double *bxr = &qr[6]; 
  const double *byr = &qr[8]; 
  const double *bzr = &qr[10]; 
  const double *phr = &qr[12]; 
  const double *psr = &qr[14]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[2]; 
  double *outEzr = &outr[4]; 
  double *outBxr = &outr[6]; 
  double *outByr = &outr[8]; 
  double *outBzr = &outr[10]; 
  double *outPhr = &outr[12]; 
  double *outPsr = &outr[14]; 
 
  double incr[2]; 
 
  incr[0] = (-0.4330127018922193*phr[1]*c2*chi)+0.4330127018922193*phl[1]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi+0.4330127018922193*exr[1]*c*chi+0.4330127018922193*exl[1]*c*chi-0.25*exr[0]*c*chi+0.25*exl[0]*c*chi; 
  incr[1] = 0.75*phr[1]*c2*chi-0.75*phl[1]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi-0.75*exr[1]*c*chi-0.75*exl[1]*c*chi+0.4330127018922193*exr[0]*c*chi-0.4330127018922193*exl[0]*c*chi; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += incr[1]*dx1; 

 
  incr[0] = (-0.4330127018922193*bzr[1]*c2)+0.4330127018922193*bzl[1]*c2+0.25*bzr[0]*c2+0.25*bzl[0]*c2+0.4330127018922193*eyr[1]*c+0.4330127018922193*eyl[1]*c-0.25*eyr[0]*c+0.25*eyl[0]*c; 
  incr[1] = 0.75*bzr[1]*c2-0.75*bzl[1]*c2-0.4330127018922193*bzr[0]*c2-0.4330127018922193*bzl[0]*c2-0.75*eyr[1]*c-0.75*eyl[1]*c+0.4330127018922193*eyr[0]*c-0.4330127018922193*eyl[0]*c; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += incr[1]*dx1; 

 
  incr[0] = 0.4330127018922193*byr[1]*c2-0.4330127018922193*byl[1]*c2-0.25*byr[0]*c2-0.25*byl[0]*c2+0.4330127018922193*ezr[1]*c+0.4330127018922193*ezl[1]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = (-0.75*byr[1]*c2)+0.75*byl[1]*c2+0.4330127018922193*byr[0]*c2+0.4330127018922193*byl[0]*c2-0.75*ezr[1]*c-0.75*ezl[1]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += incr[1]*dx1; 

 
  incr[0] = 0.4330127018922193*bxr[1]*c*gamma+0.4330127018922193*bxl[1]*c*gamma-0.25*bxr[0]*c*gamma+0.25*bxl[0]*c*gamma-0.4330127018922193*psr[1]*gamma+0.4330127018922193*psl[1]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = (-0.75*bxr[1]*c*gamma)-0.75*bxl[1]*c*gamma+0.4330127018922193*bxr[0]*c*gamma-0.4330127018922193*bxl[0]*c*gamma+0.75*psr[1]*gamma-0.75*psl[1]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += incr[1]*dx1; 

 
  incr[0] = 0.4330127018922193*byr[1]*c+0.4330127018922193*byl[1]*c-0.25*byr[0]*c+0.25*byl[0]*c+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*ezr[0]-0.25*ezl[0]; 
  incr[1] = (-0.75*byr[1]*c)-0.75*byl[1]*c+0.4330127018922193*byr[0]*c-0.4330127018922193*byl[0]*c-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*ezr[0]+0.4330127018922193*ezl[0]; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += incr[1]*dx1; 

 
  incr[0] = 0.4330127018922193*bzr[1]*c+0.4330127018922193*bzl[1]*c-0.25*bzr[0]*c+0.25*bzl[0]*c-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*eyr[0]+0.25*eyl[0]; 
  incr[1] = (-0.75*bzr[1]*c)-0.75*bzl[1]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += incr[1]*dx1; 

 
  incr[0] = 0.4330127018922193*phr[1]*c*chi+0.4330127018922193*phl[1]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi-0.4330127018922193*exr[1]*chi+0.4330127018922193*exl[1]*chi+0.25*exr[0]*chi+0.25*exl[0]*chi; 
  incr[1] = (-0.75*phr[1]*c*chi)-0.75*phl[1]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi+0.75*exr[1]*chi-0.75*exl[1]*chi-0.4330127018922193*exr[0]*chi-0.4330127018922193*exl[0]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += incr[1]*dx1; 

 
  incr[0] = (-0.4330127018922193*bxr[1]*c2*gamma)+0.4330127018922193*bxl[1]*c2*gamma+0.25*bxr[0]*c2*gamma+0.25*bxl[0]*c2*gamma+0.4330127018922193*psr[1]*c*gamma+0.4330127018922193*psl[1]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = 0.75*bxr[1]*c2*gamma-0.75*bxl[1]*c2*gamma-0.4330127018922193*bxr[0]*c2*gamma-0.4330127018922193*bxl[0]*c2*gamma-0.75*psr[1]*c*gamma-0.75*psl[1]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += incr[1]*dx1; 

 
  return c; 
} 
double MaxwellSurf1xSer_X_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[3]; 
  const double *ezl = &ql[6]; 
  const double *bxl = &ql[9]; 
  const double *byl = &ql[12]; 
  const double *bzl = &ql[15]; 
  const double *phl = &ql[18]; 
  const double *psl = &ql[21]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[3]; 
  double *outEzl = &outl[6]; 
  double *outBxl = &outl[9]; 
  double *outByl = &outl[12]; 
  double *outBzl = &outl[15]; 
  double *outPhl = &outl[18]; 
  double *outPsl = &outl[21]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[3]; 
  const double *ezr = &qr[6]; 
  const double *bxr = &qr[9]; 
  const double *byr = &qr[12]; 
  const double *bzr = &qr[15]; 
  const double *phr = &qr[18]; 
  const double *psr = &qr[21]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[3]; 
  double *outEzr = &outr[6]; 
  double *outBxr = &outr[9]; 
  double *outByr = &outr[12]; 
  double *outBzr = &outr[15]; 
  double *outPhr = &outr[18]; 
  double *outPsr = &outr[21]; 
 
  double incr[3]; 
 
  incr[0] = 0.5590169943749475*phr[2]*c2*chi+0.5590169943749475*phl[2]*c2*chi-0.4330127018922193*phr[1]*c2*chi+0.4330127018922193*phl[1]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi-0.5590169943749475*exr[2]*c*chi+0.5590169943749475*exl[2]*c*chi+0.4330127018922193*exr[1]*c*chi+0.4330127018922193*exl[1]*c*chi-0.25*exr[0]*c*chi+0.25*exl[0]*c*chi; 
  incr[1] = (-0.9682458365518543*phr[2]*c2*chi)-0.9682458365518543*phl[2]*c2*chi+0.75*phr[1]*c2*chi-0.75*phl[1]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi+0.9682458365518543*exr[2]*c*chi-0.9682458365518543*exl[2]*c*chi-0.75*exr[1]*c*chi-0.75*exl[1]*c*chi+0.4330127018922193*exr[0]*c*chi-0.4330127018922193*exl[0]*c*chi; 
  incr[2] = 1.25*phr[2]*c2*chi+1.25*phl[2]*c2*chi-0.9682458365518543*phr[1]*c2*chi+0.9682458365518543*phl[1]*c2*chi+0.5590169943749475*phr[0]*c2*chi+0.5590169943749475*phl[0]*c2*chi-1.25*exr[2]*c*chi+1.25*exl[2]*c*chi+0.9682458365518543*exr[1]*c*chi+0.9682458365518543*exl[1]*c*chi-0.5590169943749475*exr[0]*c*chi+0.5590169943749475*exl[0]*c*chi; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += incr[1]*dx1; 
  outExl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = 0.5590169943749475*bzr[2]*c2+0.5590169943749475*bzl[2]*c2-0.4330127018922193*bzr[1]*c2+0.4330127018922193*bzl[1]*c2+0.25*bzr[0]*c2+0.25*bzl[0]*c2-0.5590169943749475*eyr[2]*c+0.5590169943749475*eyl[2]*c+0.4330127018922193*eyr[1]*c+0.4330127018922193*eyl[1]*c-0.25*eyr[0]*c+0.25*eyl[0]*c; 
  incr[1] = (-0.9682458365518543*bzr[2]*c2)-0.9682458365518543*bzl[2]*c2+0.75*bzr[1]*c2-0.75*bzl[1]*c2-0.4330127018922193*bzr[0]*c2-0.4330127018922193*bzl[0]*c2+0.9682458365518543*eyr[2]*c-0.9682458365518543*eyl[2]*c-0.75*eyr[1]*c-0.75*eyl[1]*c+0.4330127018922193*eyr[0]*c-0.4330127018922193*eyl[0]*c; 
  incr[2] = 1.25*bzr[2]*c2+1.25*bzl[2]*c2-0.9682458365518543*bzr[1]*c2+0.9682458365518543*bzl[1]*c2+0.5590169943749475*bzr[0]*c2+0.5590169943749475*bzl[0]*c2-1.25*eyr[2]*c+1.25*eyl[2]*c+0.9682458365518543*eyr[1]*c+0.9682458365518543*eyl[1]*c-0.5590169943749475*eyr[0]*c+0.5590169943749475*eyl[0]*c; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += incr[1]*dx1; 
  outEyl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = (-0.5590169943749475*byr[2]*c2)-0.5590169943749475*byl[2]*c2+0.4330127018922193*byr[1]*c2-0.4330127018922193*byl[1]*c2-0.25*byr[0]*c2-0.25*byl[0]*c2-0.5590169943749475*ezr[2]*c+0.5590169943749475*ezl[2]*c+0.4330127018922193*ezr[1]*c+0.4330127018922193*ezl[1]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = 0.9682458365518543*byr[2]*c2+0.9682458365518543*byl[2]*c2-0.75*byr[1]*c2+0.75*byl[1]*c2+0.4330127018922193*byr[0]*c2+0.4330127018922193*byl[0]*c2+0.9682458365518543*ezr[2]*c-0.9682458365518543*ezl[2]*c-0.75*ezr[1]*c-0.75*ezl[1]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 
  incr[2] = (-1.25*byr[2]*c2)-1.25*byl[2]*c2+0.9682458365518543*byr[1]*c2-0.9682458365518543*byl[1]*c2-0.5590169943749475*byr[0]*c2-0.5590169943749475*byl[0]*c2-1.25*ezr[2]*c+1.25*ezl[2]*c+0.9682458365518543*ezr[1]*c+0.9682458365518543*ezl[1]*c-0.5590169943749475*ezr[0]*c+0.5590169943749475*ezl[0]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += incr[1]*dx1; 
  outEzl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = (-0.5590169943749475*bxr[2]*c*gamma)+0.5590169943749475*bxl[2]*c*gamma+0.4330127018922193*bxr[1]*c*gamma+0.4330127018922193*bxl[1]*c*gamma-0.25*bxr[0]*c*gamma+0.25*bxl[0]*c*gamma+0.5590169943749475*psr[2]*gamma+0.5590169943749475*psl[2]*gamma-0.4330127018922193*psr[1]*gamma+0.4330127018922193*psl[1]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = 0.9682458365518543*bxr[2]*c*gamma-0.9682458365518543*bxl[2]*c*gamma-0.75*bxr[1]*c*gamma-0.75*bxl[1]*c*gamma+0.4330127018922193*bxr[0]*c*gamma-0.4330127018922193*bxl[0]*c*gamma-0.9682458365518543*psr[2]*gamma-0.9682458365518543*psl[2]*gamma+0.75*psr[1]*gamma-0.75*psl[1]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 
  incr[2] = (-1.25*bxr[2]*c*gamma)+1.25*bxl[2]*c*gamma+0.9682458365518543*bxr[1]*c*gamma+0.9682458365518543*bxl[1]*c*gamma-0.5590169943749475*bxr[0]*c*gamma+0.5590169943749475*bxl[0]*c*gamma+1.25*psr[2]*gamma+1.25*psl[2]*gamma-0.9682458365518543*psr[1]*gamma+0.9682458365518543*psl[1]*gamma+0.5590169943749475*psr[0]*gamma+0.5590169943749475*psl[0]*gamma; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += incr[1]*dx1; 
  outBxl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = (-0.5590169943749475*byr[2]*c)+0.5590169943749475*byl[2]*c+0.4330127018922193*byr[1]*c+0.4330127018922193*byl[1]*c-0.25*byr[0]*c+0.25*byl[0]*c-0.5590169943749475*ezr[2]-0.5590169943749475*ezl[2]+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*ezr[0]-0.25*ezl[0]; 
  incr[1] = 0.9682458365518543*byr[2]*c-0.9682458365518543*byl[2]*c-0.75*byr[1]*c-0.75*byl[1]*c+0.4330127018922193*byr[0]*c-0.4330127018922193*byl[0]*c+0.9682458365518543*ezr[2]+0.9682458365518543*ezl[2]-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*ezr[0]+0.4330127018922193*ezl[0]; 
  incr[2] = (-1.25*byr[2]*c)+1.25*byl[2]*c+0.9682458365518543*byr[1]*c+0.9682458365518543*byl[1]*c-0.5590169943749475*byr[0]*c+0.5590169943749475*byl[0]*c-1.25*ezr[2]-1.25*ezl[2]+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*ezr[0]-0.5590169943749475*ezl[0]; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += incr[1]*dx1; 
  outByl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = (-0.5590169943749475*bzr[2]*c)+0.5590169943749475*bzl[2]*c+0.4330127018922193*bzr[1]*c+0.4330127018922193*bzl[1]*c-0.25*bzr[0]*c+0.25*bzl[0]*c+0.5590169943749475*eyr[2]+0.5590169943749475*eyl[2]-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*eyr[0]+0.25*eyl[0]; 
  incr[1] = 0.9682458365518543*bzr[2]*c-0.9682458365518543*bzl[2]*c-0.75*bzr[1]*c-0.75*bzl[1]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c-0.9682458365518543*eyr[2]-0.9682458365518543*eyl[2]+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0]; 
  incr[2] = (-1.25*bzr[2]*c)+1.25*bzl[2]*c+0.9682458365518543*bzr[1]*c+0.9682458365518543*bzl[1]*c-0.5590169943749475*bzr[0]*c+0.5590169943749475*bzl[0]*c+1.25*eyr[2]+1.25*eyl[2]-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += incr[1]*dx1; 
  outBzl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = (-0.5590169943749475*phr[2]*c*chi)+0.5590169943749475*phl[2]*c*chi+0.4330127018922193*phr[1]*c*chi+0.4330127018922193*phl[1]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi+0.5590169943749475*exr[2]*chi+0.5590169943749475*exl[2]*chi-0.4330127018922193*exr[1]*chi+0.4330127018922193*exl[1]*chi+0.25*exr[0]*chi+0.25*exl[0]*chi; 
  incr[1] = 0.9682458365518543*phr[2]*c*chi-0.9682458365518543*phl[2]*c*chi-0.75*phr[1]*c*chi-0.75*phl[1]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi-0.9682458365518543*exr[2]*chi-0.9682458365518543*exl[2]*chi+0.75*exr[1]*chi-0.75*exl[1]*chi-0.4330127018922193*exr[0]*chi-0.4330127018922193*exl[0]*chi; 
  incr[2] = (-1.25*phr[2]*c*chi)+1.25*phl[2]*c*chi+0.9682458365518543*phr[1]*c*chi+0.9682458365518543*phl[1]*c*chi-0.5590169943749475*phr[0]*c*chi+0.5590169943749475*phl[0]*c*chi+1.25*exr[2]*chi+1.25*exl[2]*chi-0.9682458365518543*exr[1]*chi+0.9682458365518543*exl[1]*chi+0.5590169943749475*exr[0]*chi+0.5590169943749475*exl[0]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += incr[1]*dx1; 
  outPhl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = 0.5590169943749475*bxr[2]*c2*gamma+0.5590169943749475*bxl[2]*c2*gamma-0.4330127018922193*bxr[1]*c2*gamma+0.4330127018922193*bxl[1]*c2*gamma+0.25*bxr[0]*c2*gamma+0.25*bxl[0]*c2*gamma-0.5590169943749475*psr[2]*c*gamma+0.5590169943749475*psl[2]*c*gamma+0.4330127018922193*psr[1]*c*gamma+0.4330127018922193*psl[1]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = (-0.9682458365518543*bxr[2]*c2*gamma)-0.9682458365518543*bxl[2]*c2*gamma+0.75*bxr[1]*c2*gamma-0.75*bxl[1]*c2*gamma-0.4330127018922193*bxr[0]*c2*gamma-0.4330127018922193*bxl[0]*c2*gamma+0.9682458365518543*psr[2]*c*gamma-0.9682458365518543*psl[2]*c*gamma-0.75*psr[1]*c*gamma-0.75*psl[1]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 
  incr[2] = 1.25*bxr[2]*c2*gamma+1.25*bxl[2]*c2*gamma-0.9682458365518543*bxr[1]*c2*gamma+0.9682458365518543*bxl[1]*c2*gamma+0.5590169943749475*bxr[0]*c2*gamma+0.5590169943749475*bxl[0]*c2*gamma-1.25*psr[2]*c*gamma+1.25*psl[2]*c*gamma+0.9682458365518543*psr[1]*c*gamma+0.9682458365518543*psl[1]*c*gamma-0.5590169943749475*psr[0]*c*gamma+0.5590169943749475*psl[0]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += incr[1]*dx1; 
  outPsl[2] += -1.0*incr[2]*dx1; 

 
  return c; 
} 
double MaxwellSurf1xSer_X_P3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[4]; 
  const double *ezl = &ql[8]; 
  const double *bxl = &ql[12]; 
  const double *byl = &ql[16]; 
  const double *bzl = &ql[20]; 
  const double *phl = &ql[24]; 
  const double *psl = &ql[28]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[4]; 
  double *outEzl = &outl[8]; 
  double *outBxl = &outl[12]; 
  double *outByl = &outl[16]; 
  double *outBzl = &outl[20]; 
  double *outPhl = &outl[24]; 
  double *outPsl = &outl[28]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[4]; 
  const double *ezr = &qr[8]; 
  const double *bxr = &qr[12]; 
  const double *byr = &qr[16]; 
  const double *bzr = &qr[20]; 
  const double *phr = &qr[24]; 
  const double *psr = &qr[28]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[4]; 
  double *outEzr = &outr[8]; 
  double *outBxr = &outr[12]; 
  double *outByr = &outr[16]; 
  double *outBzr = &outr[20]; 
  double *outPhr = &outr[24]; 
  double *outPsr = &outr[28]; 
 
  double incr[4]; 
 
  incr[0] = (-0.6614378277661477*phr[3]*c2*chi)+0.6614378277661477*phl[3]*c2*chi+0.5590169943749475*phr[2]*c2*chi+0.5590169943749475*phl[2]*c2*chi-0.4330127018922193*phr[1]*c2*chi+0.4330127018922193*phl[1]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi+0.6614378277661477*exr[3]*c*chi+0.6614378277661477*exl[3]*c*chi-0.5590169943749475*exr[2]*c*chi+0.5590169943749475*exl[2]*c*chi+0.4330127018922193*exr[1]*c*chi+0.4330127018922193*exl[1]*c*chi-0.25*exr[0]*c*chi+0.25*exl[0]*c*chi; 
  incr[1] = 1.14564392373896*phr[3]*c2*chi-1.14564392373896*phl[3]*c2*chi-0.9682458365518543*phr[2]*c2*chi-0.9682458365518543*phl[2]*c2*chi+0.75*phr[1]*c2*chi-0.75*phl[1]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi-1.14564392373896*exr[3]*c*chi-1.14564392373896*exl[3]*c*chi+0.9682458365518543*exr[2]*c*chi-0.9682458365518543*exl[2]*c*chi-0.75*exr[1]*c*chi-0.75*exl[1]*c*chi+0.4330127018922193*exr[0]*c*chi-0.4330127018922193*exl[0]*c*chi; 
  incr[2] = (-1.479019945774904*phr[3]*c2*chi)+1.479019945774904*phl[3]*c2*chi+1.25*phr[2]*c2*chi+1.25*phl[2]*c2*chi-0.9682458365518543*phr[1]*c2*chi+0.9682458365518543*phl[1]*c2*chi+0.5590169943749475*phr[0]*c2*chi+0.5590169943749475*phl[0]*c2*chi+1.479019945774904*exr[3]*c*chi+1.479019945774904*exl[3]*c*chi-1.25*exr[2]*c*chi+1.25*exl[2]*c*chi+0.9682458365518543*exr[1]*c*chi+0.9682458365518543*exl[1]*c*chi-0.5590169943749475*exr[0]*c*chi+0.5590169943749475*exl[0]*c*chi; 
  incr[3] = 1.75*phr[3]*c2*chi-1.75*phl[3]*c2*chi-1.479019945774904*phr[2]*c2*chi-1.479019945774904*phl[2]*c2*chi+1.14564392373896*phr[1]*c2*chi-1.14564392373896*phl[1]*c2*chi-0.6614378277661477*phr[0]*c2*chi-0.6614378277661477*phl[0]*c2*chi-1.75*exr[3]*c*chi-1.75*exl[3]*c*chi+1.479019945774904*exr[2]*c*chi-1.479019945774904*exl[2]*c*chi-1.14564392373896*exr[1]*c*chi-1.14564392373896*exl[1]*c*chi+0.6614378277661477*exr[0]*c*chi-0.6614378277661477*exl[0]*c*chi; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 
  outExr[3] += incr[3]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += incr[1]*dx1; 
  outExl[2] += -1.0*incr[2]*dx1; 
  outExl[3] += incr[3]*dx1; 

 
  incr[0] = (-0.6614378277661477*bzr[3]*c2)+0.6614378277661477*bzl[3]*c2+0.5590169943749475*bzr[2]*c2+0.5590169943749475*bzl[2]*c2-0.4330127018922193*bzr[1]*c2+0.4330127018922193*bzl[1]*c2+0.25*bzr[0]*c2+0.25*bzl[0]*c2+0.6614378277661477*eyr[3]*c+0.6614378277661477*eyl[3]*c-0.5590169943749475*eyr[2]*c+0.5590169943749475*eyl[2]*c+0.4330127018922193*eyr[1]*c+0.4330127018922193*eyl[1]*c-0.25*eyr[0]*c+0.25*eyl[0]*c; 
  incr[1] = 1.14564392373896*bzr[3]*c2-1.14564392373896*bzl[3]*c2-0.9682458365518543*bzr[2]*c2-0.9682458365518543*bzl[2]*c2+0.75*bzr[1]*c2-0.75*bzl[1]*c2-0.4330127018922193*bzr[0]*c2-0.4330127018922193*bzl[0]*c2-1.14564392373896*eyr[3]*c-1.14564392373896*eyl[3]*c+0.9682458365518543*eyr[2]*c-0.9682458365518543*eyl[2]*c-0.75*eyr[1]*c-0.75*eyl[1]*c+0.4330127018922193*eyr[0]*c-0.4330127018922193*eyl[0]*c; 
  incr[2] = (-1.479019945774904*bzr[3]*c2)+1.479019945774904*bzl[3]*c2+1.25*bzr[2]*c2+1.25*bzl[2]*c2-0.9682458365518543*bzr[1]*c2+0.9682458365518543*bzl[1]*c2+0.5590169943749475*bzr[0]*c2+0.5590169943749475*bzl[0]*c2+1.479019945774904*eyr[3]*c+1.479019945774904*eyl[3]*c-1.25*eyr[2]*c+1.25*eyl[2]*c+0.9682458365518543*eyr[1]*c+0.9682458365518543*eyl[1]*c-0.5590169943749475*eyr[0]*c+0.5590169943749475*eyl[0]*c; 
  incr[3] = 1.75*bzr[3]*c2-1.75*bzl[3]*c2-1.479019945774904*bzr[2]*c2-1.479019945774904*bzl[2]*c2+1.14564392373896*bzr[1]*c2-1.14564392373896*bzl[1]*c2-0.6614378277661477*bzr[0]*c2-0.6614378277661477*bzl[0]*c2-1.75*eyr[3]*c-1.75*eyl[3]*c+1.479019945774904*eyr[2]*c-1.479019945774904*eyl[2]*c-1.14564392373896*eyr[1]*c-1.14564392373896*eyl[1]*c+0.6614378277661477*eyr[0]*c-0.6614378277661477*eyl[0]*c; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 
  outEyr[3] += incr[3]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += incr[1]*dx1; 
  outEyl[2] += -1.0*incr[2]*dx1; 
  outEyl[3] += incr[3]*dx1; 

 
  incr[0] = 0.6614378277661477*byr[3]*c2-0.6614378277661477*byl[3]*c2-0.5590169943749475*byr[2]*c2-0.5590169943749475*byl[2]*c2+0.4330127018922193*byr[1]*c2-0.4330127018922193*byl[1]*c2-0.25*byr[0]*c2-0.25*byl[0]*c2+0.6614378277661477*ezr[3]*c+0.6614378277661477*ezl[3]*c-0.5590169943749475*ezr[2]*c+0.5590169943749475*ezl[2]*c+0.4330127018922193*ezr[1]*c+0.4330127018922193*ezl[1]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = (-1.14564392373896*byr[3]*c2)+1.14564392373896*byl[3]*c2+0.9682458365518543*byr[2]*c2+0.9682458365518543*byl[2]*c2-0.75*byr[1]*c2+0.75*byl[1]*c2+0.4330127018922193*byr[0]*c2+0.4330127018922193*byl[0]*c2-1.14564392373896*ezr[3]*c-1.14564392373896*ezl[3]*c+0.9682458365518543*ezr[2]*c-0.9682458365518543*ezl[2]*c-0.75*ezr[1]*c-0.75*ezl[1]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 
  incr[2] = 1.479019945774904*byr[3]*c2-1.479019945774904*byl[3]*c2-1.25*byr[2]*c2-1.25*byl[2]*c2+0.9682458365518543*byr[1]*c2-0.9682458365518543*byl[1]*c2-0.5590169943749475*byr[0]*c2-0.5590169943749475*byl[0]*c2+1.479019945774904*ezr[3]*c+1.479019945774904*ezl[3]*c-1.25*ezr[2]*c+1.25*ezl[2]*c+0.9682458365518543*ezr[1]*c+0.9682458365518543*ezl[1]*c-0.5590169943749475*ezr[0]*c+0.5590169943749475*ezl[0]*c; 
  incr[3] = (-1.75*byr[3]*c2)+1.75*byl[3]*c2+1.479019945774904*byr[2]*c2+1.479019945774904*byl[2]*c2-1.14564392373896*byr[1]*c2+1.14564392373896*byl[1]*c2+0.6614378277661477*byr[0]*c2+0.6614378277661477*byl[0]*c2-1.75*ezr[3]*c-1.75*ezl[3]*c+1.479019945774904*ezr[2]*c-1.479019945774904*ezl[2]*c-1.14564392373896*ezr[1]*c-1.14564392373896*ezl[1]*c+0.6614378277661477*ezr[0]*c-0.6614378277661477*ezl[0]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 
  outEzr[3] += incr[3]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += incr[1]*dx1; 
  outEzl[2] += -1.0*incr[2]*dx1; 
  outEzl[3] += incr[3]*dx1; 

 
  incr[0] = 0.6614378277661477*bxr[3]*c*gamma+0.6614378277661477*bxl[3]*c*gamma-0.5590169943749475*bxr[2]*c*gamma+0.5590169943749475*bxl[2]*c*gamma+0.4330127018922193*bxr[1]*c*gamma+0.4330127018922193*bxl[1]*c*gamma-0.25*bxr[0]*c*gamma+0.25*bxl[0]*c*gamma-0.6614378277661477*psr[3]*gamma+0.6614378277661477*psl[3]*gamma+0.5590169943749475*psr[2]*gamma+0.5590169943749475*psl[2]*gamma-0.4330127018922193*psr[1]*gamma+0.4330127018922193*psl[1]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = (-1.14564392373896*bxr[3]*c*gamma)-1.14564392373896*bxl[3]*c*gamma+0.9682458365518543*bxr[2]*c*gamma-0.9682458365518543*bxl[2]*c*gamma-0.75*bxr[1]*c*gamma-0.75*bxl[1]*c*gamma+0.4330127018922193*bxr[0]*c*gamma-0.4330127018922193*bxl[0]*c*gamma+1.14564392373896*psr[3]*gamma-1.14564392373896*psl[3]*gamma-0.9682458365518543*psr[2]*gamma-0.9682458365518543*psl[2]*gamma+0.75*psr[1]*gamma-0.75*psl[1]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 
  incr[2] = 1.479019945774904*bxr[3]*c*gamma+1.479019945774904*bxl[3]*c*gamma-1.25*bxr[2]*c*gamma+1.25*bxl[2]*c*gamma+0.9682458365518543*bxr[1]*c*gamma+0.9682458365518543*bxl[1]*c*gamma-0.5590169943749475*bxr[0]*c*gamma+0.5590169943749475*bxl[0]*c*gamma-1.479019945774904*psr[3]*gamma+1.479019945774904*psl[3]*gamma+1.25*psr[2]*gamma+1.25*psl[2]*gamma-0.9682458365518543*psr[1]*gamma+0.9682458365518543*psl[1]*gamma+0.5590169943749475*psr[0]*gamma+0.5590169943749475*psl[0]*gamma; 
  incr[3] = (-1.75*bxr[3]*c*gamma)-1.75*bxl[3]*c*gamma+1.479019945774904*bxr[2]*c*gamma-1.479019945774904*bxl[2]*c*gamma-1.14564392373896*bxr[1]*c*gamma-1.14564392373896*bxl[1]*c*gamma+0.6614378277661477*bxr[0]*c*gamma-0.6614378277661477*bxl[0]*c*gamma+1.75*psr[3]*gamma-1.75*psl[3]*gamma-1.479019945774904*psr[2]*gamma-1.479019945774904*psl[2]*gamma+1.14564392373896*psr[1]*gamma-1.14564392373896*psl[1]*gamma-0.6614378277661477*psr[0]*gamma-0.6614378277661477*psl[0]*gamma; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 
  outBxr[3] += incr[3]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += incr[1]*dx1; 
  outBxl[2] += -1.0*incr[2]*dx1; 
  outBxl[3] += incr[3]*dx1; 

 
  incr[0] = 0.6614378277661477*byr[3]*c+0.6614378277661477*byl[3]*c-0.5590169943749475*byr[2]*c+0.5590169943749475*byl[2]*c+0.4330127018922193*byr[1]*c+0.4330127018922193*byl[1]*c-0.25*byr[0]*c+0.25*byl[0]*c+0.6614378277661477*ezr[3]-0.6614378277661477*ezl[3]-0.5590169943749475*ezr[2]-0.5590169943749475*ezl[2]+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*ezr[0]-0.25*ezl[0]; 
  incr[1] = (-1.14564392373896*byr[3]*c)-1.14564392373896*byl[3]*c+0.9682458365518543*byr[2]*c-0.9682458365518543*byl[2]*c-0.75*byr[1]*c-0.75*byl[1]*c+0.4330127018922193*byr[0]*c-0.4330127018922193*byl[0]*c-1.14564392373896*ezr[3]+1.14564392373896*ezl[3]+0.9682458365518543*ezr[2]+0.9682458365518543*ezl[2]-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*ezr[0]+0.4330127018922193*ezl[0]; 
  incr[2] = 1.479019945774904*byr[3]*c+1.479019945774904*byl[3]*c-1.25*byr[2]*c+1.25*byl[2]*c+0.9682458365518543*byr[1]*c+0.9682458365518543*byl[1]*c-0.5590169943749475*byr[0]*c+0.5590169943749475*byl[0]*c+1.479019945774904*ezr[3]-1.479019945774904*ezl[3]-1.25*ezr[2]-1.25*ezl[2]+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*ezr[0]-0.5590169943749475*ezl[0]; 
  incr[3] = (-1.75*byr[3]*c)-1.75*byl[3]*c+1.479019945774904*byr[2]*c-1.479019945774904*byl[2]*c-1.14564392373896*byr[1]*c-1.14564392373896*byl[1]*c+0.6614378277661477*byr[0]*c-0.6614378277661477*byl[0]*c-1.75*ezr[3]+1.75*ezl[3]+1.479019945774904*ezr[2]+1.479019945774904*ezl[2]-1.14564392373896*ezr[1]+1.14564392373896*ezl[1]+0.6614378277661477*ezr[0]+0.6614378277661477*ezl[0]; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 
  outByr[3] += incr[3]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += incr[1]*dx1; 
  outByl[2] += -1.0*incr[2]*dx1; 
  outByl[3] += incr[3]*dx1; 

 
  incr[0] = 0.6614378277661477*bzr[3]*c+0.6614378277661477*bzl[3]*c-0.5590169943749475*bzr[2]*c+0.5590169943749475*bzl[2]*c+0.4330127018922193*bzr[1]*c+0.4330127018922193*bzl[1]*c-0.25*bzr[0]*c+0.25*bzl[0]*c-0.6614378277661477*eyr[3]+0.6614378277661477*eyl[3]+0.5590169943749475*eyr[2]+0.5590169943749475*eyl[2]-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*eyr[0]+0.25*eyl[0]; 
  incr[1] = (-1.14564392373896*bzr[3]*c)-1.14564392373896*bzl[3]*c+0.9682458365518543*bzr[2]*c-0.9682458365518543*bzl[2]*c-0.75*bzr[1]*c-0.75*bzl[1]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c+1.14564392373896*eyr[3]-1.14564392373896*eyl[3]-0.9682458365518543*eyr[2]-0.9682458365518543*eyl[2]+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0]; 
  incr[2] = 1.479019945774904*bzr[3]*c+1.479019945774904*bzl[3]*c-1.25*bzr[2]*c+1.25*bzl[2]*c+0.9682458365518543*bzr[1]*c+0.9682458365518543*bzl[1]*c-0.5590169943749475*bzr[0]*c+0.5590169943749475*bzl[0]*c-1.479019945774904*eyr[3]+1.479019945774904*eyl[3]+1.25*eyr[2]+1.25*eyl[2]-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0]; 
  incr[3] = (-1.75*bzr[3]*c)-1.75*bzl[3]*c+1.479019945774904*bzr[2]*c-1.479019945774904*bzl[2]*c-1.14564392373896*bzr[1]*c-1.14564392373896*bzl[1]*c+0.6614378277661477*bzr[0]*c-0.6614378277661477*bzl[0]*c+1.75*eyr[3]-1.75*eyl[3]-1.479019945774904*eyr[2]-1.479019945774904*eyl[2]+1.14564392373896*eyr[1]-1.14564392373896*eyl[1]-0.6614378277661477*eyr[0]-0.6614378277661477*eyl[0]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 
  outBzr[3] += incr[3]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += incr[1]*dx1; 
  outBzl[2] += -1.0*incr[2]*dx1; 
  outBzl[3] += incr[3]*dx1; 

 
  incr[0] = 0.6614378277661477*phr[3]*c*chi+0.6614378277661477*phl[3]*c*chi-0.5590169943749475*phr[2]*c*chi+0.5590169943749475*phl[2]*c*chi+0.4330127018922193*phr[1]*c*chi+0.4330127018922193*phl[1]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi-0.6614378277661477*exr[3]*chi+0.6614378277661477*exl[3]*chi+0.5590169943749475*exr[2]*chi+0.5590169943749475*exl[2]*chi-0.4330127018922193*exr[1]*chi+0.4330127018922193*exl[1]*chi+0.25*exr[0]*chi+0.25*exl[0]*chi; 
  incr[1] = (-1.14564392373896*phr[3]*c*chi)-1.14564392373896*phl[3]*c*chi+0.9682458365518543*phr[2]*c*chi-0.9682458365518543*phl[2]*c*chi-0.75*phr[1]*c*chi-0.75*phl[1]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi+1.14564392373896*exr[3]*chi-1.14564392373896*exl[3]*chi-0.9682458365518543*exr[2]*chi-0.9682458365518543*exl[2]*chi+0.75*exr[1]*chi-0.75*exl[1]*chi-0.4330127018922193*exr[0]*chi-0.4330127018922193*exl[0]*chi; 
  incr[2] = 1.479019945774904*phr[3]*c*chi+1.479019945774904*phl[3]*c*chi-1.25*phr[2]*c*chi+1.25*phl[2]*c*chi+0.9682458365518543*phr[1]*c*chi+0.9682458365518543*phl[1]*c*chi-0.5590169943749475*phr[0]*c*chi+0.5590169943749475*phl[0]*c*chi-1.479019945774904*exr[3]*chi+1.479019945774904*exl[3]*chi+1.25*exr[2]*chi+1.25*exl[2]*chi-0.9682458365518543*exr[1]*chi+0.9682458365518543*exl[1]*chi+0.5590169943749475*exr[0]*chi+0.5590169943749475*exl[0]*chi; 
  incr[3] = (-1.75*phr[3]*c*chi)-1.75*phl[3]*c*chi+1.479019945774904*phr[2]*c*chi-1.479019945774904*phl[2]*c*chi-1.14564392373896*phr[1]*c*chi-1.14564392373896*phl[1]*c*chi+0.6614378277661477*phr[0]*c*chi-0.6614378277661477*phl[0]*c*chi+1.75*exr[3]*chi-1.75*exl[3]*chi-1.479019945774904*exr[2]*chi-1.479019945774904*exl[2]*chi+1.14564392373896*exr[1]*chi-1.14564392373896*exl[1]*chi-0.6614378277661477*exr[0]*chi-0.6614378277661477*exl[0]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 
  outPhr[3] += incr[3]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += incr[1]*dx1; 
  outPhl[2] += -1.0*incr[2]*dx1; 
  outPhl[3] += incr[3]*dx1; 

 
  incr[0] = (-0.6614378277661477*bxr[3]*c2*gamma)+0.6614378277661477*bxl[3]*c2*gamma+0.5590169943749475*bxr[2]*c2*gamma+0.5590169943749475*bxl[2]*c2*gamma-0.4330127018922193*bxr[1]*c2*gamma+0.4330127018922193*bxl[1]*c2*gamma+0.25*bxr[0]*c2*gamma+0.25*bxl[0]*c2*gamma+0.6614378277661477*psr[3]*c*gamma+0.6614378277661477*psl[3]*c*gamma-0.5590169943749475*psr[2]*c*gamma+0.5590169943749475*psl[2]*c*gamma+0.4330127018922193*psr[1]*c*gamma+0.4330127018922193*psl[1]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = 1.14564392373896*bxr[3]*c2*gamma-1.14564392373896*bxl[3]*c2*gamma-0.9682458365518543*bxr[2]*c2*gamma-0.9682458365518543*bxl[2]*c2*gamma+0.75*bxr[1]*c2*gamma-0.75*bxl[1]*c2*gamma-0.4330127018922193*bxr[0]*c2*gamma-0.4330127018922193*bxl[0]*c2*gamma-1.14564392373896*psr[3]*c*gamma-1.14564392373896*psl[3]*c*gamma+0.9682458365518543*psr[2]*c*gamma-0.9682458365518543*psl[2]*c*gamma-0.75*psr[1]*c*gamma-0.75*psl[1]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 
  incr[2] = (-1.479019945774904*bxr[3]*c2*gamma)+1.479019945774904*bxl[3]*c2*gamma+1.25*bxr[2]*c2*gamma+1.25*bxl[2]*c2*gamma-0.9682458365518543*bxr[1]*c2*gamma+0.9682458365518543*bxl[1]*c2*gamma+0.5590169943749475*bxr[0]*c2*gamma+0.5590169943749475*bxl[0]*c2*gamma+1.479019945774904*psr[3]*c*gamma+1.479019945774904*psl[3]*c*gamma-1.25*psr[2]*c*gamma+1.25*psl[2]*c*gamma+0.9682458365518543*psr[1]*c*gamma+0.9682458365518543*psl[1]*c*gamma-0.5590169943749475*psr[0]*c*gamma+0.5590169943749475*psl[0]*c*gamma; 
  incr[3] = 1.75*bxr[3]*c2*gamma-1.75*bxl[3]*c2*gamma-1.479019945774904*bxr[2]*c2*gamma-1.479019945774904*bxl[2]*c2*gamma+1.14564392373896*bxr[1]*c2*gamma-1.14564392373896*bxl[1]*c2*gamma-0.6614378277661477*bxr[0]*c2*gamma-0.6614378277661477*bxl[0]*c2*gamma-1.75*psr[3]*c*gamma-1.75*psl[3]*c*gamma+1.479019945774904*psr[2]*c*gamma-1.479019945774904*psl[2]*c*gamma-1.14564392373896*psr[1]*c*gamma-1.14564392373896*psl[1]*c*gamma+0.6614378277661477*psr[0]*c*gamma-0.6614378277661477*psl[0]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 
  outPsr[3] += incr[3]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += incr[1]*dx1; 
  outPsl[2] += -1.0*incr[2]*dx1; 
  outPsl[3] += incr[3]*dx1; 

 
  return c; 
} 
double MaxwellSurf1xSer_X_P4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[5]; 
  const double *ezl = &ql[10]; 
  const double *bxl = &ql[15]; 
  const double *byl = &ql[20]; 
  const double *bzl = &ql[25]; 
  const double *phl = &ql[30]; 
  const double *psl = &ql[35]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[5]; 
  double *outEzl = &outl[10]; 
  double *outBxl = &outl[15]; 
  double *outByl = &outl[20]; 
  double *outBzl = &outl[25]; 
  double *outPhl = &outl[30]; 
  double *outPsl = &outl[35]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[5]; 
  const double *ezr = &qr[10]; 
  const double *bxr = &qr[15]; 
  const double *byr = &qr[20]; 
  const double *bzr = &qr[25]; 
  const double *phr = &qr[30]; 
  const double *psr = &qr[35]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[5]; 
  double *outEzr = &outr[10]; 
  double *outBxr = &outr[15]; 
  double *outByr = &outr[20]; 
  double *outBzr = &outr[25]; 
  double *outPhr = &outr[30]; 
  double *outPsr = &outr[35]; 
 
  double incr[5]; 
 
  incr[0] = 0.75*phr[4]*c2*chi+0.75*phl[4]*c2*chi-0.6614378277661477*phr[3]*c2*chi+0.6614378277661477*phl[3]*c2*chi+0.5590169943749475*phr[2]*c2*chi+0.5590169943749475*phl[2]*c2*chi-0.4330127018922193*phr[1]*c2*chi+0.4330127018922193*phl[1]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi-0.75*exr[4]*c*chi+0.75*exl[4]*c*chi+0.6614378277661477*exr[3]*c*chi+0.6614378277661477*exl[3]*c*chi-0.5590169943749475*exr[2]*c*chi+0.5590169943749475*exl[2]*c*chi+0.4330127018922193*exr[1]*c*chi+0.4330127018922193*exl[1]*c*chi-0.25*exr[0]*c*chi+0.25*exl[0]*c*chi; 
  incr[1] = (-1.299038105676658*phr[4]*c2*chi)-1.299038105676658*phl[4]*c2*chi+1.14564392373896*phr[3]*c2*chi-1.14564392373896*phl[3]*c2*chi-0.9682458365518543*phr[2]*c2*chi-0.9682458365518543*phl[2]*c2*chi+0.75*phr[1]*c2*chi-0.75*phl[1]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi+1.299038105676658*exr[4]*c*chi-1.299038105676658*exl[4]*c*chi-1.14564392373896*exr[3]*c*chi-1.14564392373896*exl[3]*c*chi+0.9682458365518543*exr[2]*c*chi-0.9682458365518543*exl[2]*c*chi-0.75*exr[1]*c*chi-0.75*exl[1]*c*chi+0.4330127018922193*exr[0]*c*chi-0.4330127018922193*exl[0]*c*chi; 
  incr[2] = 1.677050983124842*phr[4]*c2*chi+1.677050983124842*phl[4]*c2*chi-1.479019945774904*phr[3]*c2*chi+1.479019945774904*phl[3]*c2*chi+1.25*phr[2]*c2*chi+1.25*phl[2]*c2*chi-0.9682458365518543*phr[1]*c2*chi+0.9682458365518543*phl[1]*c2*chi+0.5590169943749475*phr[0]*c2*chi+0.5590169943749475*phl[0]*c2*chi-1.677050983124842*exr[4]*c*chi+1.677050983124842*exl[4]*c*chi+1.479019945774904*exr[3]*c*chi+1.479019945774904*exl[3]*c*chi-1.25*exr[2]*c*chi+1.25*exl[2]*c*chi+0.9682458365518543*exr[1]*c*chi+0.9682458365518543*exl[1]*c*chi-0.5590169943749475*exr[0]*c*chi+0.5590169943749475*exl[0]*c*chi; 
  incr[3] = (-1.984313483298443*phr[4]*c2*chi)-1.984313483298443*phl[4]*c2*chi+1.75*phr[3]*c2*chi-1.75*phl[3]*c2*chi-1.479019945774904*phr[2]*c2*chi-1.479019945774904*phl[2]*c2*chi+1.14564392373896*phr[1]*c2*chi-1.14564392373896*phl[1]*c2*chi-0.6614378277661477*phr[0]*c2*chi-0.6614378277661477*phl[0]*c2*chi+1.984313483298443*exr[4]*c*chi-1.984313483298443*exl[4]*c*chi-1.75*exr[3]*c*chi-1.75*exl[3]*c*chi+1.479019945774904*exr[2]*c*chi-1.479019945774904*exl[2]*c*chi-1.14564392373896*exr[1]*c*chi-1.14564392373896*exl[1]*c*chi+0.6614378277661477*exr[0]*c*chi-0.6614378277661477*exl[0]*c*chi; 
  incr[4] = 2.25*phr[4]*c2*chi+2.25*phl[4]*c2*chi-1.984313483298443*phr[3]*c2*chi+1.984313483298443*phl[3]*c2*chi+1.677050983124842*phr[2]*c2*chi+1.677050983124842*phl[2]*c2*chi-1.299038105676658*phr[1]*c2*chi+1.299038105676658*phl[1]*c2*chi+0.75*phr[0]*c2*chi+0.75*phl[0]*c2*chi-2.25*exr[4]*c*chi+2.25*exl[4]*c*chi+1.984313483298443*exr[3]*c*chi+1.984313483298443*exl[3]*c*chi-1.677050983124842*exr[2]*c*chi+1.677050983124842*exl[2]*c*chi+1.299038105676658*exr[1]*c*chi+1.299038105676658*exl[1]*c*chi-0.75*exr[0]*c*chi+0.75*exl[0]*c*chi; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 
  outExr[3] += incr[3]*dx1; 
  outExr[4] += incr[4]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += incr[1]*dx1; 
  outExl[2] += -1.0*incr[2]*dx1; 
  outExl[3] += incr[3]*dx1; 
  outExl[4] += -1.0*incr[4]*dx1; 

 
  incr[0] = 0.75*bzr[4]*c2+0.75*bzl[4]*c2-0.6614378277661477*bzr[3]*c2+0.6614378277661477*bzl[3]*c2+0.5590169943749475*bzr[2]*c2+0.5590169943749475*bzl[2]*c2-0.4330127018922193*bzr[1]*c2+0.4330127018922193*bzl[1]*c2+0.25*bzr[0]*c2+0.25*bzl[0]*c2-0.75*eyr[4]*c+0.75*eyl[4]*c+0.6614378277661477*eyr[3]*c+0.6614378277661477*eyl[3]*c-0.5590169943749475*eyr[2]*c+0.5590169943749475*eyl[2]*c+0.4330127018922193*eyr[1]*c+0.4330127018922193*eyl[1]*c-0.25*eyr[0]*c+0.25*eyl[0]*c; 
  incr[1] = (-1.299038105676658*bzr[4]*c2)-1.299038105676658*bzl[4]*c2+1.14564392373896*bzr[3]*c2-1.14564392373896*bzl[3]*c2-0.9682458365518543*bzr[2]*c2-0.9682458365518543*bzl[2]*c2+0.75*bzr[1]*c2-0.75*bzl[1]*c2-0.4330127018922193*bzr[0]*c2-0.4330127018922193*bzl[0]*c2+1.299038105676658*eyr[4]*c-1.299038105676658*eyl[4]*c-1.14564392373896*eyr[3]*c-1.14564392373896*eyl[3]*c+0.9682458365518543*eyr[2]*c-0.9682458365518543*eyl[2]*c-0.75*eyr[1]*c-0.75*eyl[1]*c+0.4330127018922193*eyr[0]*c-0.4330127018922193*eyl[0]*c; 
  incr[2] = 1.677050983124842*bzr[4]*c2+1.677050983124842*bzl[4]*c2-1.479019945774904*bzr[3]*c2+1.479019945774904*bzl[3]*c2+1.25*bzr[2]*c2+1.25*bzl[2]*c2-0.9682458365518543*bzr[1]*c2+0.9682458365518543*bzl[1]*c2+0.5590169943749475*bzr[0]*c2+0.5590169943749475*bzl[0]*c2-1.677050983124842*eyr[4]*c+1.677050983124842*eyl[4]*c+1.479019945774904*eyr[3]*c+1.479019945774904*eyl[3]*c-1.25*eyr[2]*c+1.25*eyl[2]*c+0.9682458365518543*eyr[1]*c+0.9682458365518543*eyl[1]*c-0.5590169943749475*eyr[0]*c+0.5590169943749475*eyl[0]*c; 
  incr[3] = (-1.984313483298443*bzr[4]*c2)-1.984313483298443*bzl[4]*c2+1.75*bzr[3]*c2-1.75*bzl[3]*c2-1.479019945774904*bzr[2]*c2-1.479019945774904*bzl[2]*c2+1.14564392373896*bzr[1]*c2-1.14564392373896*bzl[1]*c2-0.6614378277661477*bzr[0]*c2-0.6614378277661477*bzl[0]*c2+1.984313483298443*eyr[4]*c-1.984313483298443*eyl[4]*c-1.75*eyr[3]*c-1.75*eyl[3]*c+1.479019945774904*eyr[2]*c-1.479019945774904*eyl[2]*c-1.14564392373896*eyr[1]*c-1.14564392373896*eyl[1]*c+0.6614378277661477*eyr[0]*c-0.6614378277661477*eyl[0]*c; 
  incr[4] = 2.25*bzr[4]*c2+2.25*bzl[4]*c2-1.984313483298443*bzr[3]*c2+1.984313483298443*bzl[3]*c2+1.677050983124842*bzr[2]*c2+1.677050983124842*bzl[2]*c2-1.299038105676658*bzr[1]*c2+1.299038105676658*bzl[1]*c2+0.75*bzr[0]*c2+0.75*bzl[0]*c2-2.25*eyr[4]*c+2.25*eyl[4]*c+1.984313483298443*eyr[3]*c+1.984313483298443*eyl[3]*c-1.677050983124842*eyr[2]*c+1.677050983124842*eyl[2]*c+1.299038105676658*eyr[1]*c+1.299038105676658*eyl[1]*c-0.75*eyr[0]*c+0.75*eyl[0]*c; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 
  outEyr[3] += incr[3]*dx1; 
  outEyr[4] += incr[4]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += incr[1]*dx1; 
  outEyl[2] += -1.0*incr[2]*dx1; 
  outEyl[3] += incr[3]*dx1; 
  outEyl[4] += -1.0*incr[4]*dx1; 

 
  incr[0] = (-0.75*byr[4]*c2)-0.75*byl[4]*c2+0.6614378277661477*byr[3]*c2-0.6614378277661477*byl[3]*c2-0.5590169943749475*byr[2]*c2-0.5590169943749475*byl[2]*c2+0.4330127018922193*byr[1]*c2-0.4330127018922193*byl[1]*c2-0.25*byr[0]*c2-0.25*byl[0]*c2-0.75*ezr[4]*c+0.75*ezl[4]*c+0.6614378277661477*ezr[3]*c+0.6614378277661477*ezl[3]*c-0.5590169943749475*ezr[2]*c+0.5590169943749475*ezl[2]*c+0.4330127018922193*ezr[1]*c+0.4330127018922193*ezl[1]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = 1.299038105676658*byr[4]*c2+1.299038105676658*byl[4]*c2-1.14564392373896*byr[3]*c2+1.14564392373896*byl[3]*c2+0.9682458365518543*byr[2]*c2+0.9682458365518543*byl[2]*c2-0.75*byr[1]*c2+0.75*byl[1]*c2+0.4330127018922193*byr[0]*c2+0.4330127018922193*byl[0]*c2+1.299038105676658*ezr[4]*c-1.299038105676658*ezl[4]*c-1.14564392373896*ezr[3]*c-1.14564392373896*ezl[3]*c+0.9682458365518543*ezr[2]*c-0.9682458365518543*ezl[2]*c-0.75*ezr[1]*c-0.75*ezl[1]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 
  incr[2] = (-1.677050983124842*byr[4]*c2)-1.677050983124842*byl[4]*c2+1.479019945774904*byr[3]*c2-1.479019945774904*byl[3]*c2-1.25*byr[2]*c2-1.25*byl[2]*c2+0.9682458365518543*byr[1]*c2-0.9682458365518543*byl[1]*c2-0.5590169943749475*byr[0]*c2-0.5590169943749475*byl[0]*c2-1.677050983124842*ezr[4]*c+1.677050983124842*ezl[4]*c+1.479019945774904*ezr[3]*c+1.479019945774904*ezl[3]*c-1.25*ezr[2]*c+1.25*ezl[2]*c+0.9682458365518543*ezr[1]*c+0.9682458365518543*ezl[1]*c-0.5590169943749475*ezr[0]*c+0.5590169943749475*ezl[0]*c; 
  incr[3] = 1.984313483298443*byr[4]*c2+1.984313483298443*byl[4]*c2-1.75*byr[3]*c2+1.75*byl[3]*c2+1.479019945774904*byr[2]*c2+1.479019945774904*byl[2]*c2-1.14564392373896*byr[1]*c2+1.14564392373896*byl[1]*c2+0.6614378277661477*byr[0]*c2+0.6614378277661477*byl[0]*c2+1.984313483298443*ezr[4]*c-1.984313483298443*ezl[4]*c-1.75*ezr[3]*c-1.75*ezl[3]*c+1.479019945774904*ezr[2]*c-1.479019945774904*ezl[2]*c-1.14564392373896*ezr[1]*c-1.14564392373896*ezl[1]*c+0.6614378277661477*ezr[0]*c-0.6614378277661477*ezl[0]*c; 
  incr[4] = (-2.25*byr[4]*c2)-2.25*byl[4]*c2+1.984313483298443*byr[3]*c2-1.984313483298443*byl[3]*c2-1.677050983124842*byr[2]*c2-1.677050983124842*byl[2]*c2+1.299038105676658*byr[1]*c2-1.299038105676658*byl[1]*c2-0.75*byr[0]*c2-0.75*byl[0]*c2-2.25*ezr[4]*c+2.25*ezl[4]*c+1.984313483298443*ezr[3]*c+1.984313483298443*ezl[3]*c-1.677050983124842*ezr[2]*c+1.677050983124842*ezl[2]*c+1.299038105676658*ezr[1]*c+1.299038105676658*ezl[1]*c-0.75*ezr[0]*c+0.75*ezl[0]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 
  outEzr[3] += incr[3]*dx1; 
  outEzr[4] += incr[4]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += incr[1]*dx1; 
  outEzl[2] += -1.0*incr[2]*dx1; 
  outEzl[3] += incr[3]*dx1; 
  outEzl[4] += -1.0*incr[4]*dx1; 

 
  incr[0] = (-0.75*bxr[4]*c*gamma)+0.75*bxl[4]*c*gamma+0.6614378277661477*bxr[3]*c*gamma+0.6614378277661477*bxl[3]*c*gamma-0.5590169943749475*bxr[2]*c*gamma+0.5590169943749475*bxl[2]*c*gamma+0.4330127018922193*bxr[1]*c*gamma+0.4330127018922193*bxl[1]*c*gamma-0.25*bxr[0]*c*gamma+0.25*bxl[0]*c*gamma+0.75*psr[4]*gamma+0.75*psl[4]*gamma-0.6614378277661477*psr[3]*gamma+0.6614378277661477*psl[3]*gamma+0.5590169943749475*psr[2]*gamma+0.5590169943749475*psl[2]*gamma-0.4330127018922193*psr[1]*gamma+0.4330127018922193*psl[1]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = 1.299038105676658*bxr[4]*c*gamma-1.299038105676658*bxl[4]*c*gamma-1.14564392373896*bxr[3]*c*gamma-1.14564392373896*bxl[3]*c*gamma+0.9682458365518543*bxr[2]*c*gamma-0.9682458365518543*bxl[2]*c*gamma-0.75*bxr[1]*c*gamma-0.75*bxl[1]*c*gamma+0.4330127018922193*bxr[0]*c*gamma-0.4330127018922193*bxl[0]*c*gamma-1.299038105676658*psr[4]*gamma-1.299038105676658*psl[4]*gamma+1.14564392373896*psr[3]*gamma-1.14564392373896*psl[3]*gamma-0.9682458365518543*psr[2]*gamma-0.9682458365518543*psl[2]*gamma+0.75*psr[1]*gamma-0.75*psl[1]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 
  incr[2] = (-1.677050983124842*bxr[4]*c*gamma)+1.677050983124842*bxl[4]*c*gamma+1.479019945774904*bxr[3]*c*gamma+1.479019945774904*bxl[3]*c*gamma-1.25*bxr[2]*c*gamma+1.25*bxl[2]*c*gamma+0.9682458365518543*bxr[1]*c*gamma+0.9682458365518543*bxl[1]*c*gamma-0.5590169943749475*bxr[0]*c*gamma+0.5590169943749475*bxl[0]*c*gamma+1.677050983124842*psr[4]*gamma+1.677050983124842*psl[4]*gamma-1.479019945774904*psr[3]*gamma+1.479019945774904*psl[3]*gamma+1.25*psr[2]*gamma+1.25*psl[2]*gamma-0.9682458365518543*psr[1]*gamma+0.9682458365518543*psl[1]*gamma+0.5590169943749475*psr[0]*gamma+0.5590169943749475*psl[0]*gamma; 
  incr[3] = 1.984313483298443*bxr[4]*c*gamma-1.984313483298443*bxl[4]*c*gamma-1.75*bxr[3]*c*gamma-1.75*bxl[3]*c*gamma+1.479019945774904*bxr[2]*c*gamma-1.479019945774904*bxl[2]*c*gamma-1.14564392373896*bxr[1]*c*gamma-1.14564392373896*bxl[1]*c*gamma+0.6614378277661477*bxr[0]*c*gamma-0.6614378277661477*bxl[0]*c*gamma-1.984313483298443*psr[4]*gamma-1.984313483298443*psl[4]*gamma+1.75*psr[3]*gamma-1.75*psl[3]*gamma-1.479019945774904*psr[2]*gamma-1.479019945774904*psl[2]*gamma+1.14564392373896*psr[1]*gamma-1.14564392373896*psl[1]*gamma-0.6614378277661477*psr[0]*gamma-0.6614378277661477*psl[0]*gamma; 
  incr[4] = (-2.25*bxr[4]*c*gamma)+2.25*bxl[4]*c*gamma+1.984313483298443*bxr[3]*c*gamma+1.984313483298443*bxl[3]*c*gamma-1.677050983124842*bxr[2]*c*gamma+1.677050983124842*bxl[2]*c*gamma+1.299038105676658*bxr[1]*c*gamma+1.299038105676658*bxl[1]*c*gamma-0.75*bxr[0]*c*gamma+0.75*bxl[0]*c*gamma+2.25*psr[4]*gamma+2.25*psl[4]*gamma-1.984313483298443*psr[3]*gamma+1.984313483298443*psl[3]*gamma+1.677050983124842*psr[2]*gamma+1.677050983124842*psl[2]*gamma-1.299038105676658*psr[1]*gamma+1.299038105676658*psl[1]*gamma+0.75*psr[0]*gamma+0.75*psl[0]*gamma; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 
  outBxr[3] += incr[3]*dx1; 
  outBxr[4] += incr[4]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += incr[1]*dx1; 
  outBxl[2] += -1.0*incr[2]*dx1; 
  outBxl[3] += incr[3]*dx1; 
  outBxl[4] += -1.0*incr[4]*dx1; 

 
  incr[0] = (-0.75*byr[4]*c)+0.75*byl[4]*c+0.6614378277661477*byr[3]*c+0.6614378277661477*byl[3]*c-0.5590169943749475*byr[2]*c+0.5590169943749475*byl[2]*c+0.4330127018922193*byr[1]*c+0.4330127018922193*byl[1]*c-0.25*byr[0]*c+0.25*byl[0]*c-0.75*ezr[4]-0.75*ezl[4]+0.6614378277661477*ezr[3]-0.6614378277661477*ezl[3]-0.5590169943749475*ezr[2]-0.5590169943749475*ezl[2]+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*ezr[0]-0.25*ezl[0]; 
  incr[1] = 1.299038105676658*byr[4]*c-1.299038105676658*byl[4]*c-1.14564392373896*byr[3]*c-1.14564392373896*byl[3]*c+0.9682458365518543*byr[2]*c-0.9682458365518543*byl[2]*c-0.75*byr[1]*c-0.75*byl[1]*c+0.4330127018922193*byr[0]*c-0.4330127018922193*byl[0]*c+1.299038105676658*ezr[4]+1.299038105676658*ezl[4]-1.14564392373896*ezr[3]+1.14564392373896*ezl[3]+0.9682458365518543*ezr[2]+0.9682458365518543*ezl[2]-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*ezr[0]+0.4330127018922193*ezl[0]; 
  incr[2] = (-1.677050983124842*byr[4]*c)+1.677050983124842*byl[4]*c+1.479019945774904*byr[3]*c+1.479019945774904*byl[3]*c-1.25*byr[2]*c+1.25*byl[2]*c+0.9682458365518543*byr[1]*c+0.9682458365518543*byl[1]*c-0.5590169943749475*byr[0]*c+0.5590169943749475*byl[0]*c-1.677050983124842*ezr[4]-1.677050983124842*ezl[4]+1.479019945774904*ezr[3]-1.479019945774904*ezl[3]-1.25*ezr[2]-1.25*ezl[2]+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*ezr[0]-0.5590169943749475*ezl[0]; 
  incr[3] = 1.984313483298443*byr[4]*c-1.984313483298443*byl[4]*c-1.75*byr[3]*c-1.75*byl[3]*c+1.479019945774904*byr[2]*c-1.479019945774904*byl[2]*c-1.14564392373896*byr[1]*c-1.14564392373896*byl[1]*c+0.6614378277661477*byr[0]*c-0.6614378277661477*byl[0]*c+1.984313483298443*ezr[4]+1.984313483298443*ezl[4]-1.75*ezr[3]+1.75*ezl[3]+1.479019945774904*ezr[2]+1.479019945774904*ezl[2]-1.14564392373896*ezr[1]+1.14564392373896*ezl[1]+0.6614378277661477*ezr[0]+0.6614378277661477*ezl[0]; 
  incr[4] = (-2.25*byr[4]*c)+2.25*byl[4]*c+1.984313483298443*byr[3]*c+1.984313483298443*byl[3]*c-1.677050983124842*byr[2]*c+1.677050983124842*byl[2]*c+1.299038105676658*byr[1]*c+1.299038105676658*byl[1]*c-0.75*byr[0]*c+0.75*byl[0]*c-2.25*ezr[4]-2.25*ezl[4]+1.984313483298443*ezr[3]-1.984313483298443*ezl[3]-1.677050983124842*ezr[2]-1.677050983124842*ezl[2]+1.299038105676658*ezr[1]-1.299038105676658*ezl[1]-0.75*ezr[0]-0.75*ezl[0]; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 
  outByr[3] += incr[3]*dx1; 
  outByr[4] += incr[4]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += incr[1]*dx1; 
  outByl[2] += -1.0*incr[2]*dx1; 
  outByl[3] += incr[3]*dx1; 
  outByl[4] += -1.0*incr[4]*dx1; 

 
  incr[0] = (-0.75*bzr[4]*c)+0.75*bzl[4]*c+0.6614378277661477*bzr[3]*c+0.6614378277661477*bzl[3]*c-0.5590169943749475*bzr[2]*c+0.5590169943749475*bzl[2]*c+0.4330127018922193*bzr[1]*c+0.4330127018922193*bzl[1]*c-0.25*bzr[0]*c+0.25*bzl[0]*c+0.75*eyr[4]+0.75*eyl[4]-0.6614378277661477*eyr[3]+0.6614378277661477*eyl[3]+0.5590169943749475*eyr[2]+0.5590169943749475*eyl[2]-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*eyr[0]+0.25*eyl[0]; 
  incr[1] = 1.299038105676658*bzr[4]*c-1.299038105676658*bzl[4]*c-1.14564392373896*bzr[3]*c-1.14564392373896*bzl[3]*c+0.9682458365518543*bzr[2]*c-0.9682458365518543*bzl[2]*c-0.75*bzr[1]*c-0.75*bzl[1]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c-1.299038105676658*eyr[4]-1.299038105676658*eyl[4]+1.14564392373896*eyr[3]-1.14564392373896*eyl[3]-0.9682458365518543*eyr[2]-0.9682458365518543*eyl[2]+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0]; 
  incr[2] = (-1.677050983124842*bzr[4]*c)+1.677050983124842*bzl[4]*c+1.479019945774904*bzr[3]*c+1.479019945774904*bzl[3]*c-1.25*bzr[2]*c+1.25*bzl[2]*c+0.9682458365518543*bzr[1]*c+0.9682458365518543*bzl[1]*c-0.5590169943749475*bzr[0]*c+0.5590169943749475*bzl[0]*c+1.677050983124842*eyr[4]+1.677050983124842*eyl[4]-1.479019945774904*eyr[3]+1.479019945774904*eyl[3]+1.25*eyr[2]+1.25*eyl[2]-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0]; 
  incr[3] = 1.984313483298443*bzr[4]*c-1.984313483298443*bzl[4]*c-1.75*bzr[3]*c-1.75*bzl[3]*c+1.479019945774904*bzr[2]*c-1.479019945774904*bzl[2]*c-1.14564392373896*bzr[1]*c-1.14564392373896*bzl[1]*c+0.6614378277661477*bzr[0]*c-0.6614378277661477*bzl[0]*c-1.984313483298443*eyr[4]-1.984313483298443*eyl[4]+1.75*eyr[3]-1.75*eyl[3]-1.479019945774904*eyr[2]-1.479019945774904*eyl[2]+1.14564392373896*eyr[1]-1.14564392373896*eyl[1]-0.6614378277661477*eyr[0]-0.6614378277661477*eyl[0]; 
  incr[4] = (-2.25*bzr[4]*c)+2.25*bzl[4]*c+1.984313483298443*bzr[3]*c+1.984313483298443*bzl[3]*c-1.677050983124842*bzr[2]*c+1.677050983124842*bzl[2]*c+1.299038105676658*bzr[1]*c+1.299038105676658*bzl[1]*c-0.75*bzr[0]*c+0.75*bzl[0]*c+2.25*eyr[4]+2.25*eyl[4]-1.984313483298443*eyr[3]+1.984313483298443*eyl[3]+1.677050983124842*eyr[2]+1.677050983124842*eyl[2]-1.299038105676658*eyr[1]+1.299038105676658*eyl[1]+0.75*eyr[0]+0.75*eyl[0]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 
  outBzr[3] += incr[3]*dx1; 
  outBzr[4] += incr[4]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += incr[1]*dx1; 
  outBzl[2] += -1.0*incr[2]*dx1; 
  outBzl[3] += incr[3]*dx1; 
  outBzl[4] += -1.0*incr[4]*dx1; 

 
  incr[0] = (-0.75*phr[4]*c*chi)+0.75*phl[4]*c*chi+0.6614378277661477*phr[3]*c*chi+0.6614378277661477*phl[3]*c*chi-0.5590169943749475*phr[2]*c*chi+0.5590169943749475*phl[2]*c*chi+0.4330127018922193*phr[1]*c*chi+0.4330127018922193*phl[1]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi+0.75*exr[4]*chi+0.75*exl[4]*chi-0.6614378277661477*exr[3]*chi+0.6614378277661477*exl[3]*chi+0.5590169943749475*exr[2]*chi+0.5590169943749475*exl[2]*chi-0.4330127018922193*exr[1]*chi+0.4330127018922193*exl[1]*chi+0.25*exr[0]*chi+0.25*exl[0]*chi; 
  incr[1] = 1.299038105676658*phr[4]*c*chi-1.299038105676658*phl[4]*c*chi-1.14564392373896*phr[3]*c*chi-1.14564392373896*phl[3]*c*chi+0.9682458365518543*phr[2]*c*chi-0.9682458365518543*phl[2]*c*chi-0.75*phr[1]*c*chi-0.75*phl[1]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi-1.299038105676658*exr[4]*chi-1.299038105676658*exl[4]*chi+1.14564392373896*exr[3]*chi-1.14564392373896*exl[3]*chi-0.9682458365518543*exr[2]*chi-0.9682458365518543*exl[2]*chi+0.75*exr[1]*chi-0.75*exl[1]*chi-0.4330127018922193*exr[0]*chi-0.4330127018922193*exl[0]*chi; 
  incr[2] = (-1.677050983124842*phr[4]*c*chi)+1.677050983124842*phl[4]*c*chi+1.479019945774904*phr[3]*c*chi+1.479019945774904*phl[3]*c*chi-1.25*phr[2]*c*chi+1.25*phl[2]*c*chi+0.9682458365518543*phr[1]*c*chi+0.9682458365518543*phl[1]*c*chi-0.5590169943749475*phr[0]*c*chi+0.5590169943749475*phl[0]*c*chi+1.677050983124842*exr[4]*chi+1.677050983124842*exl[4]*chi-1.479019945774904*exr[3]*chi+1.479019945774904*exl[3]*chi+1.25*exr[2]*chi+1.25*exl[2]*chi-0.9682458365518543*exr[1]*chi+0.9682458365518543*exl[1]*chi+0.5590169943749475*exr[0]*chi+0.5590169943749475*exl[0]*chi; 
  incr[3] = 1.984313483298443*phr[4]*c*chi-1.984313483298443*phl[4]*c*chi-1.75*phr[3]*c*chi-1.75*phl[3]*c*chi+1.479019945774904*phr[2]*c*chi-1.479019945774904*phl[2]*c*chi-1.14564392373896*phr[1]*c*chi-1.14564392373896*phl[1]*c*chi+0.6614378277661477*phr[0]*c*chi-0.6614378277661477*phl[0]*c*chi-1.984313483298443*exr[4]*chi-1.984313483298443*exl[4]*chi+1.75*exr[3]*chi-1.75*exl[3]*chi-1.479019945774904*exr[2]*chi-1.479019945774904*exl[2]*chi+1.14564392373896*exr[1]*chi-1.14564392373896*exl[1]*chi-0.6614378277661477*exr[0]*chi-0.6614378277661477*exl[0]*chi; 
  incr[4] = (-2.25*phr[4]*c*chi)+2.25*phl[4]*c*chi+1.984313483298443*phr[3]*c*chi+1.984313483298443*phl[3]*c*chi-1.677050983124842*phr[2]*c*chi+1.677050983124842*phl[2]*c*chi+1.299038105676658*phr[1]*c*chi+1.299038105676658*phl[1]*c*chi-0.75*phr[0]*c*chi+0.75*phl[0]*c*chi+2.25*exr[4]*chi+2.25*exl[4]*chi-1.984313483298443*exr[3]*chi+1.984313483298443*exl[3]*chi+1.677050983124842*exr[2]*chi+1.677050983124842*exl[2]*chi-1.299038105676658*exr[1]*chi+1.299038105676658*exl[1]*chi+0.75*exr[0]*chi+0.75*exl[0]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 
  outPhr[3] += incr[3]*dx1; 
  outPhr[4] += incr[4]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += incr[1]*dx1; 
  outPhl[2] += -1.0*incr[2]*dx1; 
  outPhl[3] += incr[3]*dx1; 
  outPhl[4] += -1.0*incr[4]*dx1; 

 
  incr[0] = 0.75*bxr[4]*c2*gamma+0.75*bxl[4]*c2*gamma-0.6614378277661477*bxr[3]*c2*gamma+0.6614378277661477*bxl[3]*c2*gamma+0.5590169943749475*bxr[2]*c2*gamma+0.5590169943749475*bxl[2]*c2*gamma-0.4330127018922193*bxr[1]*c2*gamma+0.4330127018922193*bxl[1]*c2*gamma+0.25*bxr[0]*c2*gamma+0.25*bxl[0]*c2*gamma-0.75*psr[4]*c*gamma+0.75*psl[4]*c*gamma+0.6614378277661477*psr[3]*c*gamma+0.6614378277661477*psl[3]*c*gamma-0.5590169943749475*psr[2]*c*gamma+0.5590169943749475*psl[2]*c*gamma+0.4330127018922193*psr[1]*c*gamma+0.4330127018922193*psl[1]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = (-1.299038105676658*bxr[4]*c2*gamma)-1.299038105676658*bxl[4]*c2*gamma+1.14564392373896*bxr[3]*c2*gamma-1.14564392373896*bxl[3]*c2*gamma-0.9682458365518543*bxr[2]*c2*gamma-0.9682458365518543*bxl[2]*c2*gamma+0.75*bxr[1]*c2*gamma-0.75*bxl[1]*c2*gamma-0.4330127018922193*bxr[0]*c2*gamma-0.4330127018922193*bxl[0]*c2*gamma+1.299038105676658*psr[4]*c*gamma-1.299038105676658*psl[4]*c*gamma-1.14564392373896*psr[3]*c*gamma-1.14564392373896*psl[3]*c*gamma+0.9682458365518543*psr[2]*c*gamma-0.9682458365518543*psl[2]*c*gamma-0.75*psr[1]*c*gamma-0.75*psl[1]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 
  incr[2] = 1.677050983124842*bxr[4]*c2*gamma+1.677050983124842*bxl[4]*c2*gamma-1.479019945774904*bxr[3]*c2*gamma+1.479019945774904*bxl[3]*c2*gamma+1.25*bxr[2]*c2*gamma+1.25*bxl[2]*c2*gamma-0.9682458365518543*bxr[1]*c2*gamma+0.9682458365518543*bxl[1]*c2*gamma+0.5590169943749475*bxr[0]*c2*gamma+0.5590169943749475*bxl[0]*c2*gamma-1.677050983124842*psr[4]*c*gamma+1.677050983124842*psl[4]*c*gamma+1.479019945774904*psr[3]*c*gamma+1.479019945774904*psl[3]*c*gamma-1.25*psr[2]*c*gamma+1.25*psl[2]*c*gamma+0.9682458365518543*psr[1]*c*gamma+0.9682458365518543*psl[1]*c*gamma-0.5590169943749475*psr[0]*c*gamma+0.5590169943749475*psl[0]*c*gamma; 
  incr[3] = (-1.984313483298443*bxr[4]*c2*gamma)-1.984313483298443*bxl[4]*c2*gamma+1.75*bxr[3]*c2*gamma-1.75*bxl[3]*c2*gamma-1.479019945774904*bxr[2]*c2*gamma-1.479019945774904*bxl[2]*c2*gamma+1.14564392373896*bxr[1]*c2*gamma-1.14564392373896*bxl[1]*c2*gamma-0.6614378277661477*bxr[0]*c2*gamma-0.6614378277661477*bxl[0]*c2*gamma+1.984313483298443*psr[4]*c*gamma-1.984313483298443*psl[4]*c*gamma-1.75*psr[3]*c*gamma-1.75*psl[3]*c*gamma+1.479019945774904*psr[2]*c*gamma-1.479019945774904*psl[2]*c*gamma-1.14564392373896*psr[1]*c*gamma-1.14564392373896*psl[1]*c*gamma+0.6614378277661477*psr[0]*c*gamma-0.6614378277661477*psl[0]*c*gamma; 
  incr[4] = 2.25*bxr[4]*c2*gamma+2.25*bxl[4]*c2*gamma-1.984313483298443*bxr[3]*c2*gamma+1.984313483298443*bxl[3]*c2*gamma+1.677050983124842*bxr[2]*c2*gamma+1.677050983124842*bxl[2]*c2*gamma-1.299038105676658*bxr[1]*c2*gamma+1.299038105676658*bxl[1]*c2*gamma+0.75*bxr[0]*c2*gamma+0.75*bxl[0]*c2*gamma-2.25*psr[4]*c*gamma+2.25*psl[4]*c*gamma+1.984313483298443*psr[3]*c*gamma+1.984313483298443*psl[3]*c*gamma-1.677050983124842*psr[2]*c*gamma+1.677050983124842*psl[2]*c*gamma+1.299038105676658*psr[1]*c*gamma+1.299038105676658*psl[1]*c*gamma-0.75*psr[0]*c*gamma+0.75*psl[0]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 
  outPsr[3] += incr[3]*dx1; 
  outPsr[4] += incr[4]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += incr[1]*dx1; 
  outPsl[2] += -1.0*incr[2]*dx1; 
  outPsl[3] += incr[3]*dx1; 
  outPsl[4] += -1.0*incr[4]*dx1; 

 
  return c; 
} 
