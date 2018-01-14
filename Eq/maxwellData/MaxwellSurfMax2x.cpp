#include <MaxwellModDecl.h> 
void MaxwellSurf2xMax_X_P0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[1]; 
  const double *ezl = &ql[2]; 
  const double *bxl = &ql[3]; 
  const double *byl = &ql[4]; 
  const double *bzl = &ql[5]; 
  const double *phl = &ql[6]; 
  const double *psl = &ql[7]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[1]; 
  double *outEzl = &outl[2]; 
  double *outBxl = &outl[3]; 
  double *outByl = &outl[4]; 
  double *outBzl = &outl[5]; 
  double *outPhl = &outl[6]; 
  double *outPsl = &outl[7]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[1]; 
  const double *ezr = &qr[2]; 
  const double *bxr = &qr[3]; 
  const double *byr = &qr[4]; 
  const double *bzr = &qr[5]; 
  const double *phr = &qr[6]; 
  const double *psr = &qr[7]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[1]; 
  double *outEzr = &outr[2]; 
  double *outBxr = &outr[3]; 
  double *outByr = &outr[4]; 
  double *outBzr = &outr[5]; 
  double *outPhr = &outr[6]; 
  double *outPsr = &outr[7]; 
 
  double incr[1]; 
 
  incr[0] = 0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi-0.25*exr[0]*c*chi+0.25*exl[0]*c*chi; 

  outExr[0] += incr[0]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = 0.25*bzr[0]*c2+0.25*bzl[0]*c2-0.25*eyr[0]*c+0.25*eyl[0]*c; 

  outEyr[0] += incr[0]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = (-0.25*byr[0]*c2)-0.25*byl[0]*c2-0.25*ezr[0]*c+0.25*ezl[0]*c; 

  outEzr[0] += incr[0]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = (-0.25*bxr[0]*c*gamma)+0.25*bxl[0]*c*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 

  outBxr[0] += incr[0]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = (-0.25*byr[0]*c)+0.25*byl[0]*c-0.25*ezr[0]-0.25*ezl[0]; 

  outByr[0] += incr[0]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = (-0.25*bzr[0]*c)+0.25*bzl[0]*c+0.25*eyr[0]+0.25*eyl[0]; 

  outBzr[0] += incr[0]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = (-0.25*phr[0]*c*chi)+0.25*phl[0]*c*chi+0.25*exr[0]*chi+0.25*exl[0]*chi; 

  outPhr[0] += incr[0]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = 0.25*bxr[0]*c2*gamma+0.25*bxl[0]*c2*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 

 
} 
void MaxwellSurf2xMax_X_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
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
 
  incr[0] = (-0.4330127018922193*phr[1]*c2*chi)+0.4330127018922193*phl[1]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi+0.4330127018922193*exr[1]*c*chi+0.4330127018922193*exl[1]*c*chi-0.25*exr[0]*c*chi+0.25*exl[0]*c*chi; 
  incr[1] = 0.75*phr[1]*c2*chi-0.75*phl[1]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi-0.75*exr[1]*c*chi-0.75*exl[1]*c*chi+0.4330127018922193*exr[0]*c*chi-0.4330127018922193*exl[0]*c*chi; 
  incr[2] = 0.25*phr[2]*c2*chi+0.25*phl[2]*c2*chi-0.25*exr[2]*c*chi+0.25*exl[2]*c*chi; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += incr[1]*dx1; 
  outExl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = (-0.4330127018922193*bzr[1]*c2)+0.4330127018922193*bzl[1]*c2+0.25*bzr[0]*c2+0.25*bzl[0]*c2+0.4330127018922193*eyr[1]*c+0.4330127018922193*eyl[1]*c-0.25*eyr[0]*c+0.25*eyl[0]*c; 
  incr[1] = 0.75*bzr[1]*c2-0.75*bzl[1]*c2-0.4330127018922193*bzr[0]*c2-0.4330127018922193*bzl[0]*c2-0.75*eyr[1]*c-0.75*eyl[1]*c+0.4330127018922193*eyr[0]*c-0.4330127018922193*eyl[0]*c; 
  incr[2] = 0.25*bzr[2]*c2+0.25*bzl[2]*c2-0.25*eyr[2]*c+0.25*eyl[2]*c; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += incr[1]*dx1; 
  outEyl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = 0.4330127018922193*byr[1]*c2-0.4330127018922193*byl[1]*c2-0.25*byr[0]*c2-0.25*byl[0]*c2+0.4330127018922193*ezr[1]*c+0.4330127018922193*ezl[1]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = (-0.75*byr[1]*c2)+0.75*byl[1]*c2+0.4330127018922193*byr[0]*c2+0.4330127018922193*byl[0]*c2-0.75*ezr[1]*c-0.75*ezl[1]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 
  incr[2] = (-0.25*byr[2]*c2)-0.25*byl[2]*c2-0.25*ezr[2]*c+0.25*ezl[2]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += incr[1]*dx1; 
  outEzl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = 0.4330127018922193*bxr[1]*c*gamma+0.4330127018922193*bxl[1]*c*gamma-0.25*bxr[0]*c*gamma+0.25*bxl[0]*c*gamma-0.4330127018922193*psr[1]*gamma+0.4330127018922193*psl[1]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = (-0.75*bxr[1]*c*gamma)-0.75*bxl[1]*c*gamma+0.4330127018922193*bxr[0]*c*gamma-0.4330127018922193*bxl[0]*c*gamma+0.75*psr[1]*gamma-0.75*psl[1]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 
  incr[2] = (-0.25*bxr[2]*c*gamma)+0.25*bxl[2]*c*gamma+0.25*psr[2]*gamma+0.25*psl[2]*gamma; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += incr[1]*dx1; 
  outBxl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = 0.4330127018922193*byr[1]*c+0.4330127018922193*byl[1]*c-0.25*byr[0]*c+0.25*byl[0]*c+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*ezr[0]-0.25*ezl[0]; 
  incr[1] = (-0.75*byr[1]*c)-0.75*byl[1]*c+0.4330127018922193*byr[0]*c-0.4330127018922193*byl[0]*c-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*ezr[0]+0.4330127018922193*ezl[0]; 
  incr[2] = (-0.25*byr[2]*c)+0.25*byl[2]*c-0.25*ezr[2]-0.25*ezl[2]; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += incr[1]*dx1; 
  outByl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = 0.4330127018922193*bzr[1]*c+0.4330127018922193*bzl[1]*c-0.25*bzr[0]*c+0.25*bzl[0]*c-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*eyr[0]+0.25*eyl[0]; 
  incr[1] = (-0.75*bzr[1]*c)-0.75*bzl[1]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0]; 
  incr[2] = (-0.25*bzr[2]*c)+0.25*bzl[2]*c+0.25*eyr[2]+0.25*eyl[2]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += incr[1]*dx1; 
  outBzl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = 0.4330127018922193*phr[1]*c*chi+0.4330127018922193*phl[1]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi-0.4330127018922193*exr[1]*chi+0.4330127018922193*exl[1]*chi+0.25*exr[0]*chi+0.25*exl[0]*chi; 
  incr[1] = (-0.75*phr[1]*c*chi)-0.75*phl[1]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi+0.75*exr[1]*chi-0.75*exl[1]*chi-0.4330127018922193*exr[0]*chi-0.4330127018922193*exl[0]*chi; 
  incr[2] = (-0.25*phr[2]*c*chi)+0.25*phl[2]*c*chi+0.25*exr[2]*chi+0.25*exl[2]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += incr[1]*dx1; 
  outPhl[2] += -1.0*incr[2]*dx1; 

 
  incr[0] = (-0.4330127018922193*bxr[1]*c2*gamma)+0.4330127018922193*bxl[1]*c2*gamma+0.25*bxr[0]*c2*gamma+0.25*bxl[0]*c2*gamma+0.4330127018922193*psr[1]*c*gamma+0.4330127018922193*psl[1]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = 0.75*bxr[1]*c2*gamma-0.75*bxl[1]*c2*gamma-0.4330127018922193*bxr[0]*c2*gamma-0.4330127018922193*bxl[0]*c2*gamma-0.75*psr[1]*c*gamma-0.75*psl[1]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 
  incr[2] = 0.25*bxr[2]*c2*gamma+0.25*bxl[2]*c2*gamma-0.25*psr[2]*c*gamma+0.25*psl[2]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += incr[1]*dx1; 
  outPsl[2] += -1.0*incr[2]*dx1; 

 
} 
void MaxwellSurf2xMax_X_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[6]; 
  const double *ezl = &ql[12]; 
  const double *bxl = &ql[18]; 
  const double *byl = &ql[24]; 
  const double *bzl = &ql[30]; 
  const double *phl = &ql[36]; 
  const double *psl = &ql[42]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[6]; 
  double *outEzl = &outl[12]; 
  double *outBxl = &outl[18]; 
  double *outByl = &outl[24]; 
  double *outBzl = &outl[30]; 
  double *outPhl = &outl[36]; 
  double *outPsl = &outl[42]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[6]; 
  const double *ezr = &qr[12]; 
  const double *bxr = &qr[18]; 
  const double *byr = &qr[24]; 
  const double *bzr = &qr[30]; 
  const double *phr = &qr[36]; 
  const double *psr = &qr[42]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[6]; 
  double *outEzr = &outr[12]; 
  double *outBxr = &outr[18]; 
  double *outByr = &outr[24]; 
  double *outBzr = &outr[30]; 
  double *outPhr = &outr[36]; 
  double *outPsr = &outr[42]; 
 
  double incr[6]; 
 
  incr[0] = 0.5590169943749475*phr[4]*c2*chi+0.5590169943749475*phl[4]*c2*chi-0.4330127018922193*phr[1]*c2*chi+0.4330127018922193*phl[1]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi-0.5590169943749475*exr[4]*c*chi+0.5590169943749475*exl[4]*c*chi+0.4330127018922193*exr[1]*c*chi+0.4330127018922193*exl[1]*c*chi-0.25*exr[0]*c*chi+0.25*exl[0]*c*chi; 
  incr[1] = (-0.9682458365518543*phr[4]*c2*chi)-0.9682458365518543*phl[4]*c2*chi+0.75*phr[1]*c2*chi-0.75*phl[1]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi+0.9682458365518543*exr[4]*c*chi-0.9682458365518543*exl[4]*c*chi-0.75*exr[1]*c*chi-0.75*exl[1]*c*chi+0.4330127018922193*exr[0]*c*chi-0.4330127018922193*exl[0]*c*chi; 
  incr[2] = (-0.4330127018922193*phr[3]*c2*chi)+0.4330127018922193*phl[3]*c2*chi+0.25*phr[2]*c2*chi+0.25*phl[2]*c2*chi+0.4330127018922193*exr[3]*c*chi+0.4330127018922193*exl[3]*c*chi-0.25*exr[2]*c*chi+0.25*exl[2]*c*chi; 
  incr[3] = 0.75*phr[3]*c2*chi-0.75*phl[3]*c2*chi-0.4330127018922193*phr[2]*c2*chi-0.4330127018922193*phl[2]*c2*chi-0.75*exr[3]*c*chi-0.75*exl[3]*c*chi+0.4330127018922193*exr[2]*c*chi-0.4330127018922193*exl[2]*c*chi; 
  incr[4] = 1.25*phr[4]*c2*chi+1.25*phl[4]*c2*chi-0.9682458365518543*phr[1]*c2*chi+0.9682458365518543*phl[1]*c2*chi+0.5590169943749475*phr[0]*c2*chi+0.5590169943749475*phl[0]*c2*chi-1.25*exr[4]*c*chi+1.25*exl[4]*c*chi+0.9682458365518543*exr[1]*c*chi+0.9682458365518543*exl[1]*c*chi-0.5590169943749475*exr[0]*c*chi+0.5590169943749475*exl[0]*c*chi; 
  incr[5] = 0.25*phr[5]*c2*chi+0.25*phl[5]*c2*chi-0.25*exr[5]*c*chi+0.25*exl[5]*c*chi; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 
  outExr[3] += incr[3]*dx1; 
  outExr[4] += incr[4]*dx1; 
  outExr[5] += incr[5]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += incr[1]*dx1; 
  outExl[2] += -1.0*incr[2]*dx1; 
  outExl[3] += incr[3]*dx1; 
  outExl[4] += -1.0*incr[4]*dx1; 
  outExl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = 0.5590169943749475*bzr[4]*c2+0.5590169943749475*bzl[4]*c2-0.4330127018922193*bzr[1]*c2+0.4330127018922193*bzl[1]*c2+0.25*bzr[0]*c2+0.25*bzl[0]*c2-0.5590169943749475*eyr[4]*c+0.5590169943749475*eyl[4]*c+0.4330127018922193*eyr[1]*c+0.4330127018922193*eyl[1]*c-0.25*eyr[0]*c+0.25*eyl[0]*c; 
  incr[1] = (-0.9682458365518543*bzr[4]*c2)-0.9682458365518543*bzl[4]*c2+0.75*bzr[1]*c2-0.75*bzl[1]*c2-0.4330127018922193*bzr[0]*c2-0.4330127018922193*bzl[0]*c2+0.9682458365518543*eyr[4]*c-0.9682458365518543*eyl[4]*c-0.75*eyr[1]*c-0.75*eyl[1]*c+0.4330127018922193*eyr[0]*c-0.4330127018922193*eyl[0]*c; 
  incr[2] = (-0.4330127018922193*bzr[3]*c2)+0.4330127018922193*bzl[3]*c2+0.25*bzr[2]*c2+0.25*bzl[2]*c2+0.4330127018922193*eyr[3]*c+0.4330127018922193*eyl[3]*c-0.25*eyr[2]*c+0.25*eyl[2]*c; 
  incr[3] = 0.75*bzr[3]*c2-0.75*bzl[3]*c2-0.4330127018922193*bzr[2]*c2-0.4330127018922193*bzl[2]*c2-0.75*eyr[3]*c-0.75*eyl[3]*c+0.4330127018922193*eyr[2]*c-0.4330127018922193*eyl[2]*c; 
  incr[4] = 1.25*bzr[4]*c2+1.25*bzl[4]*c2-0.9682458365518543*bzr[1]*c2+0.9682458365518543*bzl[1]*c2+0.5590169943749475*bzr[0]*c2+0.5590169943749475*bzl[0]*c2-1.25*eyr[4]*c+1.25*eyl[4]*c+0.9682458365518543*eyr[1]*c+0.9682458365518543*eyl[1]*c-0.5590169943749475*eyr[0]*c+0.5590169943749475*eyl[0]*c; 
  incr[5] = 0.25*bzr[5]*c2+0.25*bzl[5]*c2-0.25*eyr[5]*c+0.25*eyl[5]*c; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 
  outEyr[3] += incr[3]*dx1; 
  outEyr[4] += incr[4]*dx1; 
  outEyr[5] += incr[5]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += incr[1]*dx1; 
  outEyl[2] += -1.0*incr[2]*dx1; 
  outEyl[3] += incr[3]*dx1; 
  outEyl[4] += -1.0*incr[4]*dx1; 
  outEyl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = (-0.5590169943749475*byr[4]*c2)-0.5590169943749475*byl[4]*c2+0.4330127018922193*byr[1]*c2-0.4330127018922193*byl[1]*c2-0.25*byr[0]*c2-0.25*byl[0]*c2-0.5590169943749475*ezr[4]*c+0.5590169943749475*ezl[4]*c+0.4330127018922193*ezr[1]*c+0.4330127018922193*ezl[1]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = 0.9682458365518543*byr[4]*c2+0.9682458365518543*byl[4]*c2-0.75*byr[1]*c2+0.75*byl[1]*c2+0.4330127018922193*byr[0]*c2+0.4330127018922193*byl[0]*c2+0.9682458365518543*ezr[4]*c-0.9682458365518543*ezl[4]*c-0.75*ezr[1]*c-0.75*ezl[1]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 
  incr[2] = 0.4330127018922193*byr[3]*c2-0.4330127018922193*byl[3]*c2-0.25*byr[2]*c2-0.25*byl[2]*c2+0.4330127018922193*ezr[3]*c+0.4330127018922193*ezl[3]*c-0.25*ezr[2]*c+0.25*ezl[2]*c; 
  incr[3] = (-0.75*byr[3]*c2)+0.75*byl[3]*c2+0.4330127018922193*byr[2]*c2+0.4330127018922193*byl[2]*c2-0.75*ezr[3]*c-0.75*ezl[3]*c+0.4330127018922193*ezr[2]*c-0.4330127018922193*ezl[2]*c; 
  incr[4] = (-1.25*byr[4]*c2)-1.25*byl[4]*c2+0.9682458365518543*byr[1]*c2-0.9682458365518543*byl[1]*c2-0.5590169943749475*byr[0]*c2-0.5590169943749475*byl[0]*c2-1.25*ezr[4]*c+1.25*ezl[4]*c+0.9682458365518543*ezr[1]*c+0.9682458365518543*ezl[1]*c-0.5590169943749475*ezr[0]*c+0.5590169943749475*ezl[0]*c; 
  incr[5] = (-0.25*byr[5]*c2)-0.25*byl[5]*c2-0.25*ezr[5]*c+0.25*ezl[5]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 
  outEzr[3] += incr[3]*dx1; 
  outEzr[4] += incr[4]*dx1; 
  outEzr[5] += incr[5]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += incr[1]*dx1; 
  outEzl[2] += -1.0*incr[2]*dx1; 
  outEzl[3] += incr[3]*dx1; 
  outEzl[4] += -1.0*incr[4]*dx1; 
  outEzl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = (-0.5590169943749475*bxr[4]*c*gamma)+0.5590169943749475*bxl[4]*c*gamma+0.4330127018922193*bxr[1]*c*gamma+0.4330127018922193*bxl[1]*c*gamma-0.25*bxr[0]*c*gamma+0.25*bxl[0]*c*gamma+0.5590169943749475*psr[4]*gamma+0.5590169943749475*psl[4]*gamma-0.4330127018922193*psr[1]*gamma+0.4330127018922193*psl[1]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = 0.9682458365518543*bxr[4]*c*gamma-0.9682458365518543*bxl[4]*c*gamma-0.75*bxr[1]*c*gamma-0.75*bxl[1]*c*gamma+0.4330127018922193*bxr[0]*c*gamma-0.4330127018922193*bxl[0]*c*gamma-0.9682458365518543*psr[4]*gamma-0.9682458365518543*psl[4]*gamma+0.75*psr[1]*gamma-0.75*psl[1]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 
  incr[2] = 0.4330127018922193*bxr[3]*c*gamma+0.4330127018922193*bxl[3]*c*gamma-0.25*bxr[2]*c*gamma+0.25*bxl[2]*c*gamma-0.4330127018922193*psr[3]*gamma+0.4330127018922193*psl[3]*gamma+0.25*psr[2]*gamma+0.25*psl[2]*gamma; 
  incr[3] = (-0.75*bxr[3]*c*gamma)-0.75*bxl[3]*c*gamma+0.4330127018922193*bxr[2]*c*gamma-0.4330127018922193*bxl[2]*c*gamma+0.75*psr[3]*gamma-0.75*psl[3]*gamma-0.4330127018922193*psr[2]*gamma-0.4330127018922193*psl[2]*gamma; 
  incr[4] = (-1.25*bxr[4]*c*gamma)+1.25*bxl[4]*c*gamma+0.9682458365518543*bxr[1]*c*gamma+0.9682458365518543*bxl[1]*c*gamma-0.5590169943749475*bxr[0]*c*gamma+0.5590169943749475*bxl[0]*c*gamma+1.25*psr[4]*gamma+1.25*psl[4]*gamma-0.9682458365518543*psr[1]*gamma+0.9682458365518543*psl[1]*gamma+0.5590169943749475*psr[0]*gamma+0.5590169943749475*psl[0]*gamma; 
  incr[5] = (-0.25*bxr[5]*c*gamma)+0.25*bxl[5]*c*gamma+0.25*psr[5]*gamma+0.25*psl[5]*gamma; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 
  outBxr[3] += incr[3]*dx1; 
  outBxr[4] += incr[4]*dx1; 
  outBxr[5] += incr[5]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += incr[1]*dx1; 
  outBxl[2] += -1.0*incr[2]*dx1; 
  outBxl[3] += incr[3]*dx1; 
  outBxl[4] += -1.0*incr[4]*dx1; 
  outBxl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = (-0.5590169943749475*byr[4]*c)+0.5590169943749475*byl[4]*c+0.4330127018922193*byr[1]*c+0.4330127018922193*byl[1]*c-0.25*byr[0]*c+0.25*byl[0]*c-0.5590169943749475*ezr[4]-0.5590169943749475*ezl[4]+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*ezr[0]-0.25*ezl[0]; 
  incr[1] = 0.9682458365518543*byr[4]*c-0.9682458365518543*byl[4]*c-0.75*byr[1]*c-0.75*byl[1]*c+0.4330127018922193*byr[0]*c-0.4330127018922193*byl[0]*c+0.9682458365518543*ezr[4]+0.9682458365518543*ezl[4]-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*ezr[0]+0.4330127018922193*ezl[0]; 
  incr[2] = 0.4330127018922193*byr[3]*c+0.4330127018922193*byl[3]*c-0.25*byr[2]*c+0.25*byl[2]*c+0.4330127018922193*ezr[3]-0.4330127018922193*ezl[3]-0.25*ezr[2]-0.25*ezl[2]; 
  incr[3] = (-0.75*byr[3]*c)-0.75*byl[3]*c+0.4330127018922193*byr[2]*c-0.4330127018922193*byl[2]*c-0.75*ezr[3]+0.75*ezl[3]+0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]; 
  incr[4] = (-1.25*byr[4]*c)+1.25*byl[4]*c+0.9682458365518543*byr[1]*c+0.9682458365518543*byl[1]*c-0.5590169943749475*byr[0]*c+0.5590169943749475*byl[0]*c-1.25*ezr[4]-1.25*ezl[4]+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*ezr[0]-0.5590169943749475*ezl[0]; 
  incr[5] = (-0.25*byr[5]*c)+0.25*byl[5]*c-0.25*ezr[5]-0.25*ezl[5]; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 
  outByr[3] += incr[3]*dx1; 
  outByr[4] += incr[4]*dx1; 
  outByr[5] += incr[5]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += incr[1]*dx1; 
  outByl[2] += -1.0*incr[2]*dx1; 
  outByl[3] += incr[3]*dx1; 
  outByl[4] += -1.0*incr[4]*dx1; 
  outByl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = (-0.5590169943749475*bzr[4]*c)+0.5590169943749475*bzl[4]*c+0.4330127018922193*bzr[1]*c+0.4330127018922193*bzl[1]*c-0.25*bzr[0]*c+0.25*bzl[0]*c+0.5590169943749475*eyr[4]+0.5590169943749475*eyl[4]-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*eyr[0]+0.25*eyl[0]; 
  incr[1] = 0.9682458365518543*bzr[4]*c-0.9682458365518543*bzl[4]*c-0.75*bzr[1]*c-0.75*bzl[1]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c-0.9682458365518543*eyr[4]-0.9682458365518543*eyl[4]+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0]; 
  incr[2] = 0.4330127018922193*bzr[3]*c+0.4330127018922193*bzl[3]*c-0.25*bzr[2]*c+0.25*bzl[2]*c-0.4330127018922193*eyr[3]+0.4330127018922193*eyl[3]+0.25*eyr[2]+0.25*eyl[2]; 
  incr[3] = (-0.75*bzr[3]*c)-0.75*bzl[3]*c+0.4330127018922193*bzr[2]*c-0.4330127018922193*bzl[2]*c+0.75*eyr[3]-0.75*eyl[3]-0.4330127018922193*eyr[2]-0.4330127018922193*eyl[2]; 
  incr[4] = (-1.25*bzr[4]*c)+1.25*bzl[4]*c+0.9682458365518543*bzr[1]*c+0.9682458365518543*bzl[1]*c-0.5590169943749475*bzr[0]*c+0.5590169943749475*bzl[0]*c+1.25*eyr[4]+1.25*eyl[4]-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0]; 
  incr[5] = (-0.25*bzr[5]*c)+0.25*bzl[5]*c+0.25*eyr[5]+0.25*eyl[5]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 
  outBzr[3] += incr[3]*dx1; 
  outBzr[4] += incr[4]*dx1; 
  outBzr[5] += incr[5]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += incr[1]*dx1; 
  outBzl[2] += -1.0*incr[2]*dx1; 
  outBzl[3] += incr[3]*dx1; 
  outBzl[4] += -1.0*incr[4]*dx1; 
  outBzl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = (-0.5590169943749475*phr[4]*c*chi)+0.5590169943749475*phl[4]*c*chi+0.4330127018922193*phr[1]*c*chi+0.4330127018922193*phl[1]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi+0.5590169943749475*exr[4]*chi+0.5590169943749475*exl[4]*chi-0.4330127018922193*exr[1]*chi+0.4330127018922193*exl[1]*chi+0.25*exr[0]*chi+0.25*exl[0]*chi; 
  incr[1] = 0.9682458365518543*phr[4]*c*chi-0.9682458365518543*phl[4]*c*chi-0.75*phr[1]*c*chi-0.75*phl[1]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi-0.9682458365518543*exr[4]*chi-0.9682458365518543*exl[4]*chi+0.75*exr[1]*chi-0.75*exl[1]*chi-0.4330127018922193*exr[0]*chi-0.4330127018922193*exl[0]*chi; 
  incr[2] = 0.4330127018922193*phr[3]*c*chi+0.4330127018922193*phl[3]*c*chi-0.25*phr[2]*c*chi+0.25*phl[2]*c*chi-0.4330127018922193*exr[3]*chi+0.4330127018922193*exl[3]*chi+0.25*exr[2]*chi+0.25*exl[2]*chi; 
  incr[3] = (-0.75*phr[3]*c*chi)-0.75*phl[3]*c*chi+0.4330127018922193*phr[2]*c*chi-0.4330127018922193*phl[2]*c*chi+0.75*exr[3]*chi-0.75*exl[3]*chi-0.4330127018922193*exr[2]*chi-0.4330127018922193*exl[2]*chi; 
  incr[4] = (-1.25*phr[4]*c*chi)+1.25*phl[4]*c*chi+0.9682458365518543*phr[1]*c*chi+0.9682458365518543*phl[1]*c*chi-0.5590169943749475*phr[0]*c*chi+0.5590169943749475*phl[0]*c*chi+1.25*exr[4]*chi+1.25*exl[4]*chi-0.9682458365518543*exr[1]*chi+0.9682458365518543*exl[1]*chi+0.5590169943749475*exr[0]*chi+0.5590169943749475*exl[0]*chi; 
  incr[5] = (-0.25*phr[5]*c*chi)+0.25*phl[5]*c*chi+0.25*exr[5]*chi+0.25*exl[5]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 
  outPhr[3] += incr[3]*dx1; 
  outPhr[4] += incr[4]*dx1; 
  outPhr[5] += incr[5]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += incr[1]*dx1; 
  outPhl[2] += -1.0*incr[2]*dx1; 
  outPhl[3] += incr[3]*dx1; 
  outPhl[4] += -1.0*incr[4]*dx1; 
  outPhl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = 0.5590169943749475*bxr[4]*c2*gamma+0.5590169943749475*bxl[4]*c2*gamma-0.4330127018922193*bxr[1]*c2*gamma+0.4330127018922193*bxl[1]*c2*gamma+0.25*bxr[0]*c2*gamma+0.25*bxl[0]*c2*gamma-0.5590169943749475*psr[4]*c*gamma+0.5590169943749475*psl[4]*c*gamma+0.4330127018922193*psr[1]*c*gamma+0.4330127018922193*psl[1]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = (-0.9682458365518543*bxr[4]*c2*gamma)-0.9682458365518543*bxl[4]*c2*gamma+0.75*bxr[1]*c2*gamma-0.75*bxl[1]*c2*gamma-0.4330127018922193*bxr[0]*c2*gamma-0.4330127018922193*bxl[0]*c2*gamma+0.9682458365518543*psr[4]*c*gamma-0.9682458365518543*psl[4]*c*gamma-0.75*psr[1]*c*gamma-0.75*psl[1]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 
  incr[2] = (-0.4330127018922193*bxr[3]*c2*gamma)+0.4330127018922193*bxl[3]*c2*gamma+0.25*bxr[2]*c2*gamma+0.25*bxl[2]*c2*gamma+0.4330127018922193*psr[3]*c*gamma+0.4330127018922193*psl[3]*c*gamma-0.25*psr[2]*c*gamma+0.25*psl[2]*c*gamma; 
  incr[3] = 0.75*bxr[3]*c2*gamma-0.75*bxl[3]*c2*gamma-0.4330127018922193*bxr[2]*c2*gamma-0.4330127018922193*bxl[2]*c2*gamma-0.75*psr[3]*c*gamma-0.75*psl[3]*c*gamma+0.4330127018922193*psr[2]*c*gamma-0.4330127018922193*psl[2]*c*gamma; 
  incr[4] = 1.25*bxr[4]*c2*gamma+1.25*bxl[4]*c2*gamma-0.9682458365518543*bxr[1]*c2*gamma+0.9682458365518543*bxl[1]*c2*gamma+0.5590169943749475*bxr[0]*c2*gamma+0.5590169943749475*bxl[0]*c2*gamma-1.25*psr[4]*c*gamma+1.25*psl[4]*c*gamma+0.9682458365518543*psr[1]*c*gamma+0.9682458365518543*psl[1]*c*gamma-0.5590169943749475*psr[0]*c*gamma+0.5590169943749475*psl[0]*c*gamma; 
  incr[5] = 0.25*bxr[5]*c2*gamma+0.25*bxl[5]*c2*gamma-0.25*psr[5]*c*gamma+0.25*psl[5]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 
  outPsr[3] += incr[3]*dx1; 
  outPsr[4] += incr[4]*dx1; 
  outPsr[5] += incr[5]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += incr[1]*dx1; 
  outPsl[2] += -1.0*incr[2]*dx1; 
  outPsl[3] += incr[3]*dx1; 
  outPsl[4] += -1.0*incr[4]*dx1; 
  outPsl[5] += -1.0*incr[5]*dx1; 

 
} 
void MaxwellSurf2xMax_Y_P0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[1]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[1]; 
  const double *ezl = &ql[2]; 
  const double *bxl = &ql[3]; 
  const double *byl = &ql[4]; 
  const double *bzl = &ql[5]; 
  const double *phl = &ql[6]; 
  const double *psl = &ql[7]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[1]; 
  double *outEzl = &outl[2]; 
  double *outBxl = &outl[3]; 
  double *outByl = &outl[4]; 
  double *outBzl = &outl[5]; 
  double *outPhl = &outl[6]; 
  double *outPsl = &outl[7]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[1]; 
  const double *ezr = &qr[2]; 
  const double *bxr = &qr[3]; 
  const double *byr = &qr[4]; 
  const double *bzr = &qr[5]; 
  const double *phr = &qr[6]; 
  const double *psr = &qr[7]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[1]; 
  double *outEzr = &outr[2]; 
  double *outBxr = &outr[3]; 
  double *outByr = &outr[4]; 
  double *outBzr = &outr[5]; 
  double *outPhr = &outr[6]; 
  double *outPsr = &outr[7]; 
 
  double incr[1]; 
 
  incr[0] = (-0.25*bzr[0]*c2)-0.25*bzl[0]*c2-0.25*exr[0]*c+0.25*exl[0]*c; 

  outExr[0] += incr[0]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = 0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi-0.25*eyr[0]*c*chi+0.25*eyl[0]*c*chi; 

  outEyr[0] += incr[0]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = 0.25*bxr[0]*c2+0.25*bxl[0]*c2-0.25*ezr[0]*c+0.25*ezl[0]*c; 

  outEzr[0] += incr[0]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = (-0.25*bxr[0]*c)+0.25*bxl[0]*c+0.25*ezr[0]+0.25*ezl[0]; 

  outBxr[0] += incr[0]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = (-0.25*byr[0]*c*gamma)+0.25*byl[0]*c*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 

  outByr[0] += incr[0]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = (-0.25*bzr[0]*c)+0.25*bzl[0]*c-0.25*exr[0]-0.25*exl[0]; 

  outBzr[0] += incr[0]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = (-0.25*phr[0]*c*chi)+0.25*phl[0]*c*chi+0.25*eyr[0]*chi+0.25*eyl[0]*chi; 

  outPhr[0] += incr[0]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 

 
  incr[0] = 0.25*byr[0]*c2*gamma+0.25*byl[0]*c2*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 

 
} 
void MaxwellSurf2xMax_Y_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[1]; 
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
 
  incr[0] = 0.4330127018922193*bzr[2]*c2-0.4330127018922193*bzl[2]*c2-0.25*bzr[0]*c2-0.25*bzl[0]*c2+0.4330127018922193*exr[2]*c+0.4330127018922193*exl[2]*c-0.25*exr[0]*c+0.25*exl[0]*c; 
  incr[1] = (-0.25*bzr[1]*c2)-0.25*bzl[1]*c2-0.25*exr[1]*c+0.25*exl[1]*c; 
  incr[2] = (-0.75*bzr[2]*c2)+0.75*bzl[2]*c2+0.4330127018922193*bzr[0]*c2+0.4330127018922193*bzl[0]*c2-0.75*exr[2]*c-0.75*exl[2]*c+0.4330127018922193*exr[0]*c-0.4330127018922193*exl[0]*c; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += -1.0*incr[1]*dx1; 
  outExl[2] += incr[2]*dx1; 

 
  incr[0] = (-0.4330127018922193*phr[2]*c2*chi)+0.4330127018922193*phl[2]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi+0.4330127018922193*eyr[2]*c*chi+0.4330127018922193*eyl[2]*c*chi-0.25*eyr[0]*c*chi+0.25*eyl[0]*c*chi; 
  incr[1] = 0.25*phr[1]*c2*chi+0.25*phl[1]*c2*chi-0.25*eyr[1]*c*chi+0.25*eyl[1]*c*chi; 
  incr[2] = 0.75*phr[2]*c2*chi-0.75*phl[2]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi-0.75*eyr[2]*c*chi-0.75*eyl[2]*c*chi+0.4330127018922193*eyr[0]*c*chi-0.4330127018922193*eyl[0]*c*chi; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += -1.0*incr[1]*dx1; 
  outEyl[2] += incr[2]*dx1; 

 
  incr[0] = (-0.4330127018922193*bxr[2]*c2)+0.4330127018922193*bxl[2]*c2+0.25*bxr[0]*c2+0.25*bxl[0]*c2+0.4330127018922193*ezr[2]*c+0.4330127018922193*ezl[2]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = 0.25*bxr[1]*c2+0.25*bxl[1]*c2-0.25*ezr[1]*c+0.25*ezl[1]*c; 
  incr[2] = 0.75*bxr[2]*c2-0.75*bxl[2]*c2-0.4330127018922193*bxr[0]*c2-0.4330127018922193*bxl[0]*c2-0.75*ezr[2]*c-0.75*ezl[2]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += -1.0*incr[1]*dx1; 
  outEzl[2] += incr[2]*dx1; 

 
  incr[0] = 0.4330127018922193*bxr[2]*c+0.4330127018922193*bxl[2]*c-0.25*bxr[0]*c+0.25*bxl[0]*c-0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]+0.25*ezr[0]+0.25*ezl[0]; 
  incr[1] = (-0.25*bxr[1]*c)+0.25*bxl[1]*c+0.25*ezr[1]+0.25*ezl[1]; 
  incr[2] = (-0.75*bxr[2]*c)-0.75*bxl[2]*c+0.4330127018922193*bxr[0]*c-0.4330127018922193*bxl[0]*c+0.75*ezr[2]-0.75*ezl[2]-0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0]; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += -1.0*incr[1]*dx1; 
  outBxl[2] += incr[2]*dx1; 

 
  incr[0] = 0.4330127018922193*byr[2]*c*gamma+0.4330127018922193*byl[2]*c*gamma-0.25*byr[0]*c*gamma+0.25*byl[0]*c*gamma-0.4330127018922193*psr[2]*gamma+0.4330127018922193*psl[2]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = (-0.25*byr[1]*c*gamma)+0.25*byl[1]*c*gamma+0.25*psr[1]*gamma+0.25*psl[1]*gamma; 
  incr[2] = (-0.75*byr[2]*c*gamma)-0.75*byl[2]*c*gamma+0.4330127018922193*byr[0]*c*gamma-0.4330127018922193*byl[0]*c*gamma+0.75*psr[2]*gamma-0.75*psl[2]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += -1.0*incr[1]*dx1; 
  outByl[2] += incr[2]*dx1; 

 
  incr[0] = 0.4330127018922193*bzr[2]*c+0.4330127018922193*bzl[2]*c-0.25*bzr[0]*c+0.25*bzl[0]*c+0.4330127018922193*exr[2]-0.4330127018922193*exl[2]-0.25*exr[0]-0.25*exl[0]; 
  incr[1] = (-0.25*bzr[1]*c)+0.25*bzl[1]*c-0.25*exr[1]-0.25*exl[1]; 
  incr[2] = (-0.75*bzr[2]*c)-0.75*bzl[2]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c-0.75*exr[2]+0.75*exl[2]+0.4330127018922193*exr[0]+0.4330127018922193*exl[0]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += -1.0*incr[1]*dx1; 
  outBzl[2] += incr[2]*dx1; 

 
  incr[0] = 0.4330127018922193*phr[2]*c*chi+0.4330127018922193*phl[2]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi-0.4330127018922193*eyr[2]*chi+0.4330127018922193*eyl[2]*chi+0.25*eyr[0]*chi+0.25*eyl[0]*chi; 
  incr[1] = (-0.25*phr[1]*c*chi)+0.25*phl[1]*c*chi+0.25*eyr[1]*chi+0.25*eyl[1]*chi; 
  incr[2] = (-0.75*phr[2]*c*chi)-0.75*phl[2]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi+0.75*eyr[2]*chi-0.75*eyl[2]*chi-0.4330127018922193*eyr[0]*chi-0.4330127018922193*eyl[0]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += -1.0*incr[1]*dx1; 
  outPhl[2] += incr[2]*dx1; 

 
  incr[0] = (-0.4330127018922193*byr[2]*c2*gamma)+0.4330127018922193*byl[2]*c2*gamma+0.25*byr[0]*c2*gamma+0.25*byl[0]*c2*gamma+0.4330127018922193*psr[2]*c*gamma+0.4330127018922193*psl[2]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = 0.25*byr[1]*c2*gamma+0.25*byl[1]*c2*gamma-0.25*psr[1]*c*gamma+0.25*psl[1]*c*gamma; 
  incr[2] = 0.75*byr[2]*c2*gamma-0.75*byl[2]*c2*gamma-0.4330127018922193*byr[0]*c2*gamma-0.4330127018922193*byl[0]*c2*gamma-0.75*psr[2]*c*gamma-0.75*psl[2]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += -1.0*incr[1]*dx1; 
  outPsl[2] += incr[2]*dx1; 

 
} 
void MaxwellSurf2xMax_Y_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[1]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[6]; 
  const double *ezl = &ql[12]; 
  const double *bxl = &ql[18]; 
  const double *byl = &ql[24]; 
  const double *bzl = &ql[30]; 
  const double *phl = &ql[36]; 
  const double *psl = &ql[42]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[6]; 
  double *outEzl = &outl[12]; 
  double *outBxl = &outl[18]; 
  double *outByl = &outl[24]; 
  double *outBzl = &outl[30]; 
  double *outPhl = &outl[36]; 
  double *outPsl = &outl[42]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[6]; 
  const double *ezr = &qr[12]; 
  const double *bxr = &qr[18]; 
  const double *byr = &qr[24]; 
  const double *bzr = &qr[30]; 
  const double *phr = &qr[36]; 
  const double *psr = &qr[42]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[6]; 
  double *outEzr = &outr[12]; 
  double *outBxr = &outr[18]; 
  double *outByr = &outr[24]; 
  double *outBzr = &outr[30]; 
  double *outPhr = &outr[36]; 
  double *outPsr = &outr[42]; 
 
  double incr[6]; 
 
  incr[0] = (-0.5590169943749475*bzr[5]*c2)-0.5590169943749475*bzl[5]*c2+0.4330127018922193*bzr[2]*c2-0.4330127018922193*bzl[2]*c2-0.25*bzr[0]*c2-0.25*bzl[0]*c2-0.5590169943749475*exr[5]*c+0.5590169943749475*exl[5]*c+0.4330127018922193*exr[2]*c+0.4330127018922193*exl[2]*c-0.25*exr[0]*c+0.25*exl[0]*c; 
  incr[1] = 0.4330127018922193*bzr[3]*c2-0.4330127018922193*bzl[3]*c2-0.25*bzr[1]*c2-0.25*bzl[1]*c2+0.4330127018922193*exr[3]*c+0.4330127018922193*exl[3]*c-0.25*exr[1]*c+0.25*exl[1]*c; 
  incr[2] = 0.9682458365518543*bzr[5]*c2+0.9682458365518543*bzl[5]*c2-0.75*bzr[2]*c2+0.75*bzl[2]*c2+0.4330127018922193*bzr[0]*c2+0.4330127018922193*bzl[0]*c2+0.9682458365518543*exr[5]*c-0.9682458365518543*exl[5]*c-0.75*exr[2]*c-0.75*exl[2]*c+0.4330127018922193*exr[0]*c-0.4330127018922193*exl[0]*c; 
  incr[3] = (-0.75*bzr[3]*c2)+0.75*bzl[3]*c2+0.4330127018922193*bzr[1]*c2+0.4330127018922193*bzl[1]*c2-0.75*exr[3]*c-0.75*exl[3]*c+0.4330127018922193*exr[1]*c-0.4330127018922193*exl[1]*c; 
  incr[4] = (-0.25*bzr[4]*c2)-0.25*bzl[4]*c2-0.25*exr[4]*c+0.25*exl[4]*c; 
  incr[5] = (-1.25*bzr[5]*c2)-1.25*bzl[5]*c2+0.9682458365518543*bzr[2]*c2-0.9682458365518543*bzl[2]*c2-0.5590169943749475*bzr[0]*c2-0.5590169943749475*bzl[0]*c2-1.25*exr[5]*c+1.25*exl[5]*c+0.9682458365518543*exr[2]*c+0.9682458365518543*exl[2]*c-0.5590169943749475*exr[0]*c+0.5590169943749475*exl[0]*c; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 
  outExr[3] += incr[3]*dx1; 
  outExr[4] += incr[4]*dx1; 
  outExr[5] += incr[5]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += -1.0*incr[1]*dx1; 
  outExl[2] += incr[2]*dx1; 
  outExl[3] += incr[3]*dx1; 
  outExl[4] += -1.0*incr[4]*dx1; 
  outExl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = 0.5590169943749475*phr[5]*c2*chi+0.5590169943749475*phl[5]*c2*chi-0.4330127018922193*phr[2]*c2*chi+0.4330127018922193*phl[2]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi-0.5590169943749475*eyr[5]*c*chi+0.5590169943749475*eyl[5]*c*chi+0.4330127018922193*eyr[2]*c*chi+0.4330127018922193*eyl[2]*c*chi-0.25*eyr[0]*c*chi+0.25*eyl[0]*c*chi; 
  incr[1] = (-0.4330127018922193*phr[3]*c2*chi)+0.4330127018922193*phl[3]*c2*chi+0.25*phr[1]*c2*chi+0.25*phl[1]*c2*chi+0.4330127018922193*eyr[3]*c*chi+0.4330127018922193*eyl[3]*c*chi-0.25*eyr[1]*c*chi+0.25*eyl[1]*c*chi; 
  incr[2] = (-0.9682458365518543*phr[5]*c2*chi)-0.9682458365518543*phl[5]*c2*chi+0.75*phr[2]*c2*chi-0.75*phl[2]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi+0.9682458365518543*eyr[5]*c*chi-0.9682458365518543*eyl[5]*c*chi-0.75*eyr[2]*c*chi-0.75*eyl[2]*c*chi+0.4330127018922193*eyr[0]*c*chi-0.4330127018922193*eyl[0]*c*chi; 
  incr[3] = 0.75*phr[3]*c2*chi-0.75*phl[3]*c2*chi-0.4330127018922193*phr[1]*c2*chi-0.4330127018922193*phl[1]*c2*chi-0.75*eyr[3]*c*chi-0.75*eyl[3]*c*chi+0.4330127018922193*eyr[1]*c*chi-0.4330127018922193*eyl[1]*c*chi; 
  incr[4] = 0.25*phr[4]*c2*chi+0.25*phl[4]*c2*chi-0.25*eyr[4]*c*chi+0.25*eyl[4]*c*chi; 
  incr[5] = 1.25*phr[5]*c2*chi+1.25*phl[5]*c2*chi-0.9682458365518543*phr[2]*c2*chi+0.9682458365518543*phl[2]*c2*chi+0.5590169943749475*phr[0]*c2*chi+0.5590169943749475*phl[0]*c2*chi-1.25*eyr[5]*c*chi+1.25*eyl[5]*c*chi+0.9682458365518543*eyr[2]*c*chi+0.9682458365518543*eyl[2]*c*chi-0.5590169943749475*eyr[0]*c*chi+0.5590169943749475*eyl[0]*c*chi; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 
  outEyr[3] += incr[3]*dx1; 
  outEyr[4] += incr[4]*dx1; 
  outEyr[5] += incr[5]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += -1.0*incr[1]*dx1; 
  outEyl[2] += incr[2]*dx1; 
  outEyl[3] += incr[3]*dx1; 
  outEyl[4] += -1.0*incr[4]*dx1; 
  outEyl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = 0.5590169943749475*bxr[5]*c2+0.5590169943749475*bxl[5]*c2-0.4330127018922193*bxr[2]*c2+0.4330127018922193*bxl[2]*c2+0.25*bxr[0]*c2+0.25*bxl[0]*c2-0.5590169943749475*ezr[5]*c+0.5590169943749475*ezl[5]*c+0.4330127018922193*ezr[2]*c+0.4330127018922193*ezl[2]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = (-0.4330127018922193*bxr[3]*c2)+0.4330127018922193*bxl[3]*c2+0.25*bxr[1]*c2+0.25*bxl[1]*c2+0.4330127018922193*ezr[3]*c+0.4330127018922193*ezl[3]*c-0.25*ezr[1]*c+0.25*ezl[1]*c; 
  incr[2] = (-0.9682458365518543*bxr[5]*c2)-0.9682458365518543*bxl[5]*c2+0.75*bxr[2]*c2-0.75*bxl[2]*c2-0.4330127018922193*bxr[0]*c2-0.4330127018922193*bxl[0]*c2+0.9682458365518543*ezr[5]*c-0.9682458365518543*ezl[5]*c-0.75*ezr[2]*c-0.75*ezl[2]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 
  incr[3] = 0.75*bxr[3]*c2-0.75*bxl[3]*c2-0.4330127018922193*bxr[1]*c2-0.4330127018922193*bxl[1]*c2-0.75*ezr[3]*c-0.75*ezl[3]*c+0.4330127018922193*ezr[1]*c-0.4330127018922193*ezl[1]*c; 
  incr[4] = 0.25*bxr[4]*c2+0.25*bxl[4]*c2-0.25*ezr[4]*c+0.25*ezl[4]*c; 
  incr[5] = 1.25*bxr[5]*c2+1.25*bxl[5]*c2-0.9682458365518543*bxr[2]*c2+0.9682458365518543*bxl[2]*c2+0.5590169943749475*bxr[0]*c2+0.5590169943749475*bxl[0]*c2-1.25*ezr[5]*c+1.25*ezl[5]*c+0.9682458365518543*ezr[2]*c+0.9682458365518543*ezl[2]*c-0.5590169943749475*ezr[0]*c+0.5590169943749475*ezl[0]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 
  outEzr[3] += incr[3]*dx1; 
  outEzr[4] += incr[4]*dx1; 
  outEzr[5] += incr[5]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += -1.0*incr[1]*dx1; 
  outEzl[2] += incr[2]*dx1; 
  outEzl[3] += incr[3]*dx1; 
  outEzl[4] += -1.0*incr[4]*dx1; 
  outEzl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = (-0.5590169943749475*bxr[5]*c)+0.5590169943749475*bxl[5]*c+0.4330127018922193*bxr[2]*c+0.4330127018922193*bxl[2]*c-0.25*bxr[0]*c+0.25*bxl[0]*c+0.5590169943749475*ezr[5]+0.5590169943749475*ezl[5]-0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]+0.25*ezr[0]+0.25*ezl[0]; 
  incr[1] = 0.4330127018922193*bxr[3]*c+0.4330127018922193*bxl[3]*c-0.25*bxr[1]*c+0.25*bxl[1]*c-0.4330127018922193*ezr[3]+0.4330127018922193*ezl[3]+0.25*ezr[1]+0.25*ezl[1]; 
  incr[2] = 0.9682458365518543*bxr[5]*c-0.9682458365518543*bxl[5]*c-0.75*bxr[2]*c-0.75*bxl[2]*c+0.4330127018922193*bxr[0]*c-0.4330127018922193*bxl[0]*c-0.9682458365518543*ezr[5]-0.9682458365518543*ezl[5]+0.75*ezr[2]-0.75*ezl[2]-0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0]; 
  incr[3] = (-0.75*bxr[3]*c)-0.75*bxl[3]*c+0.4330127018922193*bxr[1]*c-0.4330127018922193*bxl[1]*c+0.75*ezr[3]-0.75*ezl[3]-0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]; 
  incr[4] = (-0.25*bxr[4]*c)+0.25*bxl[4]*c+0.25*ezr[4]+0.25*ezl[4]; 
  incr[5] = (-1.25*bxr[5]*c)+1.25*bxl[5]*c+0.9682458365518543*bxr[2]*c+0.9682458365518543*bxl[2]*c-0.5590169943749475*bxr[0]*c+0.5590169943749475*bxl[0]*c+1.25*ezr[5]+1.25*ezl[5]-0.9682458365518543*ezr[2]+0.9682458365518543*ezl[2]+0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0]; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 
  outBxr[3] += incr[3]*dx1; 
  outBxr[4] += incr[4]*dx1; 
  outBxr[5] += incr[5]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += -1.0*incr[1]*dx1; 
  outBxl[2] += incr[2]*dx1; 
  outBxl[3] += incr[3]*dx1; 
  outBxl[4] += -1.0*incr[4]*dx1; 
  outBxl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = (-0.5590169943749475*byr[5]*c*gamma)+0.5590169943749475*byl[5]*c*gamma+0.4330127018922193*byr[2]*c*gamma+0.4330127018922193*byl[2]*c*gamma-0.25*byr[0]*c*gamma+0.25*byl[0]*c*gamma+0.5590169943749475*psr[5]*gamma+0.5590169943749475*psl[5]*gamma-0.4330127018922193*psr[2]*gamma+0.4330127018922193*psl[2]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = 0.4330127018922193*byr[3]*c*gamma+0.4330127018922193*byl[3]*c*gamma-0.25*byr[1]*c*gamma+0.25*byl[1]*c*gamma-0.4330127018922193*psr[3]*gamma+0.4330127018922193*psl[3]*gamma+0.25*psr[1]*gamma+0.25*psl[1]*gamma; 
  incr[2] = 0.9682458365518543*byr[5]*c*gamma-0.9682458365518543*byl[5]*c*gamma-0.75*byr[2]*c*gamma-0.75*byl[2]*c*gamma+0.4330127018922193*byr[0]*c*gamma-0.4330127018922193*byl[0]*c*gamma-0.9682458365518543*psr[5]*gamma-0.9682458365518543*psl[5]*gamma+0.75*psr[2]*gamma-0.75*psl[2]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 
  incr[3] = (-0.75*byr[3]*c*gamma)-0.75*byl[3]*c*gamma+0.4330127018922193*byr[1]*c*gamma-0.4330127018922193*byl[1]*c*gamma+0.75*psr[3]*gamma-0.75*psl[3]*gamma-0.4330127018922193*psr[1]*gamma-0.4330127018922193*psl[1]*gamma; 
  incr[4] = (-0.25*byr[4]*c*gamma)+0.25*byl[4]*c*gamma+0.25*psr[4]*gamma+0.25*psl[4]*gamma; 
  incr[5] = (-1.25*byr[5]*c*gamma)+1.25*byl[5]*c*gamma+0.9682458365518543*byr[2]*c*gamma+0.9682458365518543*byl[2]*c*gamma-0.5590169943749475*byr[0]*c*gamma+0.5590169943749475*byl[0]*c*gamma+1.25*psr[5]*gamma+1.25*psl[5]*gamma-0.9682458365518543*psr[2]*gamma+0.9682458365518543*psl[2]*gamma+0.5590169943749475*psr[0]*gamma+0.5590169943749475*psl[0]*gamma; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 
  outByr[3] += incr[3]*dx1; 
  outByr[4] += incr[4]*dx1; 
  outByr[5] += incr[5]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += -1.0*incr[1]*dx1; 
  outByl[2] += incr[2]*dx1; 
  outByl[3] += incr[3]*dx1; 
  outByl[4] += -1.0*incr[4]*dx1; 
  outByl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = (-0.5590169943749475*bzr[5]*c)+0.5590169943749475*bzl[5]*c+0.4330127018922193*bzr[2]*c+0.4330127018922193*bzl[2]*c-0.25*bzr[0]*c+0.25*bzl[0]*c-0.5590169943749475*exr[5]-0.5590169943749475*exl[5]+0.4330127018922193*exr[2]-0.4330127018922193*exl[2]-0.25*exr[0]-0.25*exl[0]; 
  incr[1] = 0.4330127018922193*bzr[3]*c+0.4330127018922193*bzl[3]*c-0.25*bzr[1]*c+0.25*bzl[1]*c+0.4330127018922193*exr[3]-0.4330127018922193*exl[3]-0.25*exr[1]-0.25*exl[1]; 
  incr[2] = 0.9682458365518543*bzr[5]*c-0.9682458365518543*bzl[5]*c-0.75*bzr[2]*c-0.75*bzl[2]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c+0.9682458365518543*exr[5]+0.9682458365518543*exl[5]-0.75*exr[2]+0.75*exl[2]+0.4330127018922193*exr[0]+0.4330127018922193*exl[0]; 
  incr[3] = (-0.75*bzr[3]*c)-0.75*bzl[3]*c+0.4330127018922193*bzr[1]*c-0.4330127018922193*bzl[1]*c-0.75*exr[3]+0.75*exl[3]+0.4330127018922193*exr[1]+0.4330127018922193*exl[1]; 
  incr[4] = (-0.25*bzr[4]*c)+0.25*bzl[4]*c-0.25*exr[4]-0.25*exl[4]; 
  incr[5] = (-1.25*bzr[5]*c)+1.25*bzl[5]*c+0.9682458365518543*bzr[2]*c+0.9682458365518543*bzl[2]*c-0.5590169943749475*bzr[0]*c+0.5590169943749475*bzl[0]*c-1.25*exr[5]-1.25*exl[5]+0.9682458365518543*exr[2]-0.9682458365518543*exl[2]-0.5590169943749475*exr[0]-0.5590169943749475*exl[0]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 
  outBzr[3] += incr[3]*dx1; 
  outBzr[4] += incr[4]*dx1; 
  outBzr[5] += incr[5]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += -1.0*incr[1]*dx1; 
  outBzl[2] += incr[2]*dx1; 
  outBzl[3] += incr[3]*dx1; 
  outBzl[4] += -1.0*incr[4]*dx1; 
  outBzl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = (-0.5590169943749475*phr[5]*c*chi)+0.5590169943749475*phl[5]*c*chi+0.4330127018922193*phr[2]*c*chi+0.4330127018922193*phl[2]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi+0.5590169943749475*eyr[5]*chi+0.5590169943749475*eyl[5]*chi-0.4330127018922193*eyr[2]*chi+0.4330127018922193*eyl[2]*chi+0.25*eyr[0]*chi+0.25*eyl[0]*chi; 
  incr[1] = 0.4330127018922193*phr[3]*c*chi+0.4330127018922193*phl[3]*c*chi-0.25*phr[1]*c*chi+0.25*phl[1]*c*chi-0.4330127018922193*eyr[3]*chi+0.4330127018922193*eyl[3]*chi+0.25*eyr[1]*chi+0.25*eyl[1]*chi; 
  incr[2] = 0.9682458365518543*phr[5]*c*chi-0.9682458365518543*phl[5]*c*chi-0.75*phr[2]*c*chi-0.75*phl[2]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi-0.9682458365518543*eyr[5]*chi-0.9682458365518543*eyl[5]*chi+0.75*eyr[2]*chi-0.75*eyl[2]*chi-0.4330127018922193*eyr[0]*chi-0.4330127018922193*eyl[0]*chi; 
  incr[3] = (-0.75*phr[3]*c*chi)-0.75*phl[3]*c*chi+0.4330127018922193*phr[1]*c*chi-0.4330127018922193*phl[1]*c*chi+0.75*eyr[3]*chi-0.75*eyl[3]*chi-0.4330127018922193*eyr[1]*chi-0.4330127018922193*eyl[1]*chi; 
  incr[4] = (-0.25*phr[4]*c*chi)+0.25*phl[4]*c*chi+0.25*eyr[4]*chi+0.25*eyl[4]*chi; 
  incr[5] = (-1.25*phr[5]*c*chi)+1.25*phl[5]*c*chi+0.9682458365518543*phr[2]*c*chi+0.9682458365518543*phl[2]*c*chi-0.5590169943749475*phr[0]*c*chi+0.5590169943749475*phl[0]*c*chi+1.25*eyr[5]*chi+1.25*eyl[5]*chi-0.9682458365518543*eyr[2]*chi+0.9682458365518543*eyl[2]*chi+0.5590169943749475*eyr[0]*chi+0.5590169943749475*eyl[0]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 
  outPhr[3] += incr[3]*dx1; 
  outPhr[4] += incr[4]*dx1; 
  outPhr[5] += incr[5]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += -1.0*incr[1]*dx1; 
  outPhl[2] += incr[2]*dx1; 
  outPhl[3] += incr[3]*dx1; 
  outPhl[4] += -1.0*incr[4]*dx1; 
  outPhl[5] += -1.0*incr[5]*dx1; 

 
  incr[0] = 0.5590169943749475*byr[5]*c2*gamma+0.5590169943749475*byl[5]*c2*gamma-0.4330127018922193*byr[2]*c2*gamma+0.4330127018922193*byl[2]*c2*gamma+0.25*byr[0]*c2*gamma+0.25*byl[0]*c2*gamma-0.5590169943749475*psr[5]*c*gamma+0.5590169943749475*psl[5]*c*gamma+0.4330127018922193*psr[2]*c*gamma+0.4330127018922193*psl[2]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = (-0.4330127018922193*byr[3]*c2*gamma)+0.4330127018922193*byl[3]*c2*gamma+0.25*byr[1]*c2*gamma+0.25*byl[1]*c2*gamma+0.4330127018922193*psr[3]*c*gamma+0.4330127018922193*psl[3]*c*gamma-0.25*psr[1]*c*gamma+0.25*psl[1]*c*gamma; 
  incr[2] = (-0.9682458365518543*byr[5]*c2*gamma)-0.9682458365518543*byl[5]*c2*gamma+0.75*byr[2]*c2*gamma-0.75*byl[2]*c2*gamma-0.4330127018922193*byr[0]*c2*gamma-0.4330127018922193*byl[0]*c2*gamma+0.9682458365518543*psr[5]*c*gamma-0.9682458365518543*psl[5]*c*gamma-0.75*psr[2]*c*gamma-0.75*psl[2]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 
  incr[3] = 0.75*byr[3]*c2*gamma-0.75*byl[3]*c2*gamma-0.4330127018922193*byr[1]*c2*gamma-0.4330127018922193*byl[1]*c2*gamma-0.75*psr[3]*c*gamma-0.75*psl[3]*c*gamma+0.4330127018922193*psr[1]*c*gamma-0.4330127018922193*psl[1]*c*gamma; 
  incr[4] = 0.25*byr[4]*c2*gamma+0.25*byl[4]*c2*gamma-0.25*psr[4]*c*gamma+0.25*psl[4]*c*gamma; 
  incr[5] = 1.25*byr[5]*c2*gamma+1.25*byl[5]*c2*gamma-0.9682458365518543*byr[2]*c2*gamma+0.9682458365518543*byl[2]*c2*gamma+0.5590169943749475*byr[0]*c2*gamma+0.5590169943749475*byl[0]*c2*gamma-1.25*psr[5]*c*gamma+1.25*psl[5]*c*gamma+0.9682458365518543*psr[2]*c*gamma+0.9682458365518543*psl[2]*c*gamma-0.5590169943749475*psr[0]*c*gamma+0.5590169943749475*psl[0]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 
  outPsr[3] += incr[3]*dx1; 
  outPsr[4] += incr[4]*dx1; 
  outPsr[5] += incr[5]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += -1.0*incr[1]*dx1; 
  outPsl[2] += incr[2]*dx1; 
  outPsl[3] += incr[3]*dx1; 
  outPsl[4] += -1.0*incr[4]*dx1; 
  outPsl[5] += -1.0*incr[5]*dx1; 

 
} 
