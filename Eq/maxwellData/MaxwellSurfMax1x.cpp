#include <MaxwellModDecl.h> 
void MaxwellSurf1xMax_X_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
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
 
  const double *exr = &ql[0]; 
  const double *eyr = &ql[2]; 
  const double *ezr = &ql[4]; 
  const double *bxr = &ql[6]; 
  const double *byr = &ql[8]; 
  const double *bzr = &ql[10]; 
  const double *phr = &ql[12]; 
  const double *psr = &ql[14]; 
 
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

 
} 
void MaxwellSurf1xMax_X_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
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
 
  const double *exr = &ql[0]; 
  const double *eyr = &ql[3]; 
  const double *ezr = &ql[6]; 
  const double *bxr = &ql[9]; 
  const double *byr = &ql[12]; 
  const double *bzr = &ql[15]; 
  const double *phr = &ql[18]; 
  const double *psr = &ql[21]; 
 
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

 
} 
