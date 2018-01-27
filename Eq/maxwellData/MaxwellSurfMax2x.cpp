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
void MaxwellSurf2xMax_X_P3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[10]; 
  const double *ezl = &ql[20]; 
  const double *bxl = &ql[30]; 
  const double *byl = &ql[40]; 
  const double *bzl = &ql[50]; 
  const double *phl = &ql[60]; 
  const double *psl = &ql[70]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[10]; 
  double *outEzl = &outl[20]; 
  double *outBxl = &outl[30]; 
  double *outByl = &outl[40]; 
  double *outBzl = &outl[50]; 
  double *outPhl = &outl[60]; 
  double *outPsl = &outl[70]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[10]; 
  const double *ezr = &qr[20]; 
  const double *bxr = &qr[30]; 
  const double *byr = &qr[40]; 
  const double *bzr = &qr[50]; 
  const double *phr = &qr[60]; 
  const double *psr = &qr[70]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[10]; 
  double *outEzr = &outr[20]; 
  double *outBxr = &outr[30]; 
  double *outByr = &outr[40]; 
  double *outBzr = &outr[50]; 
  double *outPhr = &outr[60]; 
  double *outPsr = &outr[70]; 
 
  double incr[10]; 
 
  incr[0] = (-0.6614378277661477*phr[8]*c2*chi)+0.6614378277661477*phl[8]*c2*chi+0.5590169943749475*phr[4]*c2*chi+0.5590169943749475*phl[4]*c2*chi-0.4330127018922193*phr[1]*c2*chi+0.4330127018922193*phl[1]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi+0.6614378277661477*exr[8]*c*chi+0.6614378277661477*exl[8]*c*chi-0.5590169943749475*exr[4]*c*chi+0.5590169943749475*exl[4]*c*chi+0.4330127018922193*exr[1]*c*chi+0.4330127018922193*exl[1]*c*chi-0.25*exr[0]*c*chi+0.25*exl[0]*c*chi; 
  incr[1] = 1.14564392373896*phr[8]*c2*chi-1.14564392373896*phl[8]*c2*chi-0.9682458365518543*phr[4]*c2*chi-0.9682458365518543*phl[4]*c2*chi+0.75*phr[1]*c2*chi-0.75*phl[1]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi-1.14564392373896*exr[8]*c*chi-1.14564392373896*exl[8]*c*chi+0.9682458365518543*exr[4]*c*chi-0.9682458365518543*exl[4]*c*chi-0.75*exr[1]*c*chi-0.75*exl[1]*c*chi+0.4330127018922193*exr[0]*c*chi-0.4330127018922193*exl[0]*c*chi; 
  incr[2] = 0.5590169943749476*phr[6]*c2*chi+0.5590169943749476*phl[6]*c2*chi-0.4330127018922193*phr[3]*c2*chi+0.4330127018922193*phl[3]*c2*chi+0.25*phr[2]*c2*chi+0.25*phl[2]*c2*chi-0.5590169943749476*exr[6]*c*chi+0.5590169943749476*exl[6]*c*chi+0.4330127018922193*exr[3]*c*chi+0.4330127018922193*exl[3]*c*chi-0.25*exr[2]*c*chi+0.25*exl[2]*c*chi; 
  incr[3] = (-0.9682458365518543*phr[6]*c2*chi)-0.9682458365518543*phl[6]*c2*chi+0.75*phr[3]*c2*chi-0.75*phl[3]*c2*chi-0.4330127018922193*phr[2]*c2*chi-0.4330127018922193*phl[2]*c2*chi+0.9682458365518543*exr[6]*c*chi-0.9682458365518543*exl[6]*c*chi-0.75*exr[3]*c*chi-0.75*exl[3]*c*chi+0.4330127018922193*exr[2]*c*chi-0.4330127018922193*exl[2]*c*chi; 
  incr[4] = (-1.479019945774904*phr[8]*c2*chi)+1.479019945774904*phl[8]*c2*chi+1.25*phr[4]*c2*chi+1.25*phl[4]*c2*chi-0.9682458365518543*phr[1]*c2*chi+0.9682458365518543*phl[1]*c2*chi+0.5590169943749475*phr[0]*c2*chi+0.5590169943749475*phl[0]*c2*chi+1.479019945774904*exr[8]*c*chi+1.479019945774904*exl[8]*c*chi-1.25*exr[4]*c*chi+1.25*exl[4]*c*chi+0.9682458365518543*exr[1]*c*chi+0.9682458365518543*exl[1]*c*chi-0.5590169943749475*exr[0]*c*chi+0.5590169943749475*exl[0]*c*chi; 
  incr[5] = (-0.4330127018922194*phr[7]*c2*chi)+0.4330127018922194*phl[7]*c2*chi+0.25*phr[5]*c2*chi+0.25*phl[5]*c2*chi+0.4330127018922194*exr[7]*c*chi+0.4330127018922194*exl[7]*c*chi-0.25*exr[5]*c*chi+0.25*exl[5]*c*chi; 
  incr[6] = 1.25*phr[6]*c2*chi+1.25*phl[6]*c2*chi-0.9682458365518543*phr[3]*c2*chi+0.9682458365518543*phl[3]*c2*chi+0.5590169943749476*phr[2]*c2*chi+0.5590169943749476*phl[2]*c2*chi-1.25*exr[6]*c*chi+1.25*exl[6]*c*chi+0.9682458365518543*exr[3]*c*chi+0.9682458365518543*exl[3]*c*chi-0.5590169943749476*exr[2]*c*chi+0.5590169943749476*exl[2]*c*chi; 
  incr[7] = 0.75*phr[7]*c2*chi-0.75*phl[7]*c2*chi-0.4330127018922194*phr[5]*c2*chi-0.4330127018922194*phl[5]*c2*chi-0.75*exr[7]*c*chi-0.75*exl[7]*c*chi+0.4330127018922194*exr[5]*c*chi-0.4330127018922194*exl[5]*c*chi; 
  incr[8] = 1.75*phr[8]*c2*chi-1.75*phl[8]*c2*chi-1.479019945774904*phr[4]*c2*chi-1.479019945774904*phl[4]*c2*chi+1.14564392373896*phr[1]*c2*chi-1.14564392373896*phl[1]*c2*chi-0.6614378277661477*phr[0]*c2*chi-0.6614378277661477*phl[0]*c2*chi-1.75*exr[8]*c*chi-1.75*exl[8]*c*chi+1.479019945774904*exr[4]*c*chi-1.479019945774904*exl[4]*c*chi-1.14564392373896*exr[1]*c*chi-1.14564392373896*exl[1]*c*chi+0.6614378277661477*exr[0]*c*chi-0.6614378277661477*exl[0]*c*chi; 
  incr[9] = 0.25*phr[9]*c2*chi+0.25*phl[9]*c2*chi-0.25*exr[9]*c*chi+0.25*exl[9]*c*chi; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 
  outExr[3] += incr[3]*dx1; 
  outExr[4] += incr[4]*dx1; 
  outExr[5] += incr[5]*dx1; 
  outExr[6] += incr[6]*dx1; 
  outExr[7] += incr[7]*dx1; 
  outExr[8] += incr[8]*dx1; 
  outExr[9] += incr[9]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += incr[1]*dx1; 
  outExl[2] += -1.0*incr[2]*dx1; 
  outExl[3] += incr[3]*dx1; 
  outExl[4] += -1.0*incr[4]*dx1; 
  outExl[5] += -1.0*incr[5]*dx1; 
  outExl[6] += -1.0*incr[6]*dx1; 
  outExl[7] += incr[7]*dx1; 
  outExl[8] += incr[8]*dx1; 
  outExl[9] += -1.0*incr[9]*dx1; 

 
  incr[0] = (-0.6614378277661477*bzr[8]*c2)+0.6614378277661477*bzl[8]*c2+0.5590169943749475*bzr[4]*c2+0.5590169943749475*bzl[4]*c2-0.4330127018922193*bzr[1]*c2+0.4330127018922193*bzl[1]*c2+0.25*bzr[0]*c2+0.25*bzl[0]*c2+0.6614378277661477*eyr[8]*c+0.6614378277661477*eyl[8]*c-0.5590169943749475*eyr[4]*c+0.5590169943749475*eyl[4]*c+0.4330127018922193*eyr[1]*c+0.4330127018922193*eyl[1]*c-0.25*eyr[0]*c+0.25*eyl[0]*c; 
  incr[1] = 1.14564392373896*bzr[8]*c2-1.14564392373896*bzl[8]*c2-0.9682458365518543*bzr[4]*c2-0.9682458365518543*bzl[4]*c2+0.75*bzr[1]*c2-0.75*bzl[1]*c2-0.4330127018922193*bzr[0]*c2-0.4330127018922193*bzl[0]*c2-1.14564392373896*eyr[8]*c-1.14564392373896*eyl[8]*c+0.9682458365518543*eyr[4]*c-0.9682458365518543*eyl[4]*c-0.75*eyr[1]*c-0.75*eyl[1]*c+0.4330127018922193*eyr[0]*c-0.4330127018922193*eyl[0]*c; 
  incr[2] = 0.5590169943749476*bzr[6]*c2+0.5590169943749476*bzl[6]*c2-0.4330127018922193*bzr[3]*c2+0.4330127018922193*bzl[3]*c2+0.25*bzr[2]*c2+0.25*bzl[2]*c2-0.5590169943749476*eyr[6]*c+0.5590169943749476*eyl[6]*c+0.4330127018922193*eyr[3]*c+0.4330127018922193*eyl[3]*c-0.25*eyr[2]*c+0.25*eyl[2]*c; 
  incr[3] = (-0.9682458365518543*bzr[6]*c2)-0.9682458365518543*bzl[6]*c2+0.75*bzr[3]*c2-0.75*bzl[3]*c2-0.4330127018922193*bzr[2]*c2-0.4330127018922193*bzl[2]*c2+0.9682458365518543*eyr[6]*c-0.9682458365518543*eyl[6]*c-0.75*eyr[3]*c-0.75*eyl[3]*c+0.4330127018922193*eyr[2]*c-0.4330127018922193*eyl[2]*c; 
  incr[4] = (-1.479019945774904*bzr[8]*c2)+1.479019945774904*bzl[8]*c2+1.25*bzr[4]*c2+1.25*bzl[4]*c2-0.9682458365518543*bzr[1]*c2+0.9682458365518543*bzl[1]*c2+0.5590169943749475*bzr[0]*c2+0.5590169943749475*bzl[0]*c2+1.479019945774904*eyr[8]*c+1.479019945774904*eyl[8]*c-1.25*eyr[4]*c+1.25*eyl[4]*c+0.9682458365518543*eyr[1]*c+0.9682458365518543*eyl[1]*c-0.5590169943749475*eyr[0]*c+0.5590169943749475*eyl[0]*c; 
  incr[5] = (-0.4330127018922194*bzr[7]*c2)+0.4330127018922194*bzl[7]*c2+0.25*bzr[5]*c2+0.25*bzl[5]*c2+0.4330127018922194*eyr[7]*c+0.4330127018922194*eyl[7]*c-0.25*eyr[5]*c+0.25*eyl[5]*c; 
  incr[6] = 1.25*bzr[6]*c2+1.25*bzl[6]*c2-0.9682458365518543*bzr[3]*c2+0.9682458365518543*bzl[3]*c2+0.5590169943749476*bzr[2]*c2+0.5590169943749476*bzl[2]*c2-1.25*eyr[6]*c+1.25*eyl[6]*c+0.9682458365518543*eyr[3]*c+0.9682458365518543*eyl[3]*c-0.5590169943749476*eyr[2]*c+0.5590169943749476*eyl[2]*c; 
  incr[7] = 0.75*bzr[7]*c2-0.75*bzl[7]*c2-0.4330127018922194*bzr[5]*c2-0.4330127018922194*bzl[5]*c2-0.75*eyr[7]*c-0.75*eyl[7]*c+0.4330127018922194*eyr[5]*c-0.4330127018922194*eyl[5]*c; 
  incr[8] = 1.75*bzr[8]*c2-1.75*bzl[8]*c2-1.479019945774904*bzr[4]*c2-1.479019945774904*bzl[4]*c2+1.14564392373896*bzr[1]*c2-1.14564392373896*bzl[1]*c2-0.6614378277661477*bzr[0]*c2-0.6614378277661477*bzl[0]*c2-1.75*eyr[8]*c-1.75*eyl[8]*c+1.479019945774904*eyr[4]*c-1.479019945774904*eyl[4]*c-1.14564392373896*eyr[1]*c-1.14564392373896*eyl[1]*c+0.6614378277661477*eyr[0]*c-0.6614378277661477*eyl[0]*c; 
  incr[9] = 0.25*bzr[9]*c2+0.25*bzl[9]*c2-0.25*eyr[9]*c+0.25*eyl[9]*c; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 
  outEyr[3] += incr[3]*dx1; 
  outEyr[4] += incr[4]*dx1; 
  outEyr[5] += incr[5]*dx1; 
  outEyr[6] += incr[6]*dx1; 
  outEyr[7] += incr[7]*dx1; 
  outEyr[8] += incr[8]*dx1; 
  outEyr[9] += incr[9]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += incr[1]*dx1; 
  outEyl[2] += -1.0*incr[2]*dx1; 
  outEyl[3] += incr[3]*dx1; 
  outEyl[4] += -1.0*incr[4]*dx1; 
  outEyl[5] += -1.0*incr[5]*dx1; 
  outEyl[6] += -1.0*incr[6]*dx1; 
  outEyl[7] += incr[7]*dx1; 
  outEyl[8] += incr[8]*dx1; 
  outEyl[9] += -1.0*incr[9]*dx1; 

 
  incr[0] = 0.6614378277661477*byr[8]*c2-0.6614378277661477*byl[8]*c2-0.5590169943749475*byr[4]*c2-0.5590169943749475*byl[4]*c2+0.4330127018922193*byr[1]*c2-0.4330127018922193*byl[1]*c2-0.25*byr[0]*c2-0.25*byl[0]*c2+0.6614378277661477*ezr[8]*c+0.6614378277661477*ezl[8]*c-0.5590169943749475*ezr[4]*c+0.5590169943749475*ezl[4]*c+0.4330127018922193*ezr[1]*c+0.4330127018922193*ezl[1]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = (-1.14564392373896*byr[8]*c2)+1.14564392373896*byl[8]*c2+0.9682458365518543*byr[4]*c2+0.9682458365518543*byl[4]*c2-0.75*byr[1]*c2+0.75*byl[1]*c2+0.4330127018922193*byr[0]*c2+0.4330127018922193*byl[0]*c2-1.14564392373896*ezr[8]*c-1.14564392373896*ezl[8]*c+0.9682458365518543*ezr[4]*c-0.9682458365518543*ezl[4]*c-0.75*ezr[1]*c-0.75*ezl[1]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 
  incr[2] = (-0.5590169943749476*byr[6]*c2)-0.5590169943749476*byl[6]*c2+0.4330127018922193*byr[3]*c2-0.4330127018922193*byl[3]*c2-0.25*byr[2]*c2-0.25*byl[2]*c2-0.5590169943749476*ezr[6]*c+0.5590169943749476*ezl[6]*c+0.4330127018922193*ezr[3]*c+0.4330127018922193*ezl[3]*c-0.25*ezr[2]*c+0.25*ezl[2]*c; 
  incr[3] = 0.9682458365518543*byr[6]*c2+0.9682458365518543*byl[6]*c2-0.75*byr[3]*c2+0.75*byl[3]*c2+0.4330127018922193*byr[2]*c2+0.4330127018922193*byl[2]*c2+0.9682458365518543*ezr[6]*c-0.9682458365518543*ezl[6]*c-0.75*ezr[3]*c-0.75*ezl[3]*c+0.4330127018922193*ezr[2]*c-0.4330127018922193*ezl[2]*c; 
  incr[4] = 1.479019945774904*byr[8]*c2-1.479019945774904*byl[8]*c2-1.25*byr[4]*c2-1.25*byl[4]*c2+0.9682458365518543*byr[1]*c2-0.9682458365518543*byl[1]*c2-0.5590169943749475*byr[0]*c2-0.5590169943749475*byl[0]*c2+1.479019945774904*ezr[8]*c+1.479019945774904*ezl[8]*c-1.25*ezr[4]*c+1.25*ezl[4]*c+0.9682458365518543*ezr[1]*c+0.9682458365518543*ezl[1]*c-0.5590169943749475*ezr[0]*c+0.5590169943749475*ezl[0]*c; 
  incr[5] = 0.4330127018922194*byr[7]*c2-0.4330127018922194*byl[7]*c2-0.25*byr[5]*c2-0.25*byl[5]*c2+0.4330127018922194*ezr[7]*c+0.4330127018922194*ezl[7]*c-0.25*ezr[5]*c+0.25*ezl[5]*c; 
  incr[6] = (-1.25*byr[6]*c2)-1.25*byl[6]*c2+0.9682458365518543*byr[3]*c2-0.9682458365518543*byl[3]*c2-0.5590169943749476*byr[2]*c2-0.5590169943749476*byl[2]*c2-1.25*ezr[6]*c+1.25*ezl[6]*c+0.9682458365518543*ezr[3]*c+0.9682458365518543*ezl[3]*c-0.5590169943749476*ezr[2]*c+0.5590169943749476*ezl[2]*c; 
  incr[7] = (-0.75*byr[7]*c2)+0.75*byl[7]*c2+0.4330127018922194*byr[5]*c2+0.4330127018922194*byl[5]*c2-0.75*ezr[7]*c-0.75*ezl[7]*c+0.4330127018922194*ezr[5]*c-0.4330127018922194*ezl[5]*c; 
  incr[8] = (-1.75*byr[8]*c2)+1.75*byl[8]*c2+1.479019945774904*byr[4]*c2+1.479019945774904*byl[4]*c2-1.14564392373896*byr[1]*c2+1.14564392373896*byl[1]*c2+0.6614378277661477*byr[0]*c2+0.6614378277661477*byl[0]*c2-1.75*ezr[8]*c-1.75*ezl[8]*c+1.479019945774904*ezr[4]*c-1.479019945774904*ezl[4]*c-1.14564392373896*ezr[1]*c-1.14564392373896*ezl[1]*c+0.6614378277661477*ezr[0]*c-0.6614378277661477*ezl[0]*c; 
  incr[9] = (-0.25*byr[9]*c2)-0.25*byl[9]*c2-0.25*ezr[9]*c+0.25*ezl[9]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 
  outEzr[3] += incr[3]*dx1; 
  outEzr[4] += incr[4]*dx1; 
  outEzr[5] += incr[5]*dx1; 
  outEzr[6] += incr[6]*dx1; 
  outEzr[7] += incr[7]*dx1; 
  outEzr[8] += incr[8]*dx1; 
  outEzr[9] += incr[9]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += incr[1]*dx1; 
  outEzl[2] += -1.0*incr[2]*dx1; 
  outEzl[3] += incr[3]*dx1; 
  outEzl[4] += -1.0*incr[4]*dx1; 
  outEzl[5] += -1.0*incr[5]*dx1; 
  outEzl[6] += -1.0*incr[6]*dx1; 
  outEzl[7] += incr[7]*dx1; 
  outEzl[8] += incr[8]*dx1; 
  outEzl[9] += -1.0*incr[9]*dx1; 

 
  incr[0] = 0.6614378277661477*bxr[8]*c*gamma+0.6614378277661477*bxl[8]*c*gamma-0.5590169943749475*bxr[4]*c*gamma+0.5590169943749475*bxl[4]*c*gamma+0.4330127018922193*bxr[1]*c*gamma+0.4330127018922193*bxl[1]*c*gamma-0.25*bxr[0]*c*gamma+0.25*bxl[0]*c*gamma-0.6614378277661477*psr[8]*gamma+0.6614378277661477*psl[8]*gamma+0.5590169943749475*psr[4]*gamma+0.5590169943749475*psl[4]*gamma-0.4330127018922193*psr[1]*gamma+0.4330127018922193*psl[1]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = (-1.14564392373896*bxr[8]*c*gamma)-1.14564392373896*bxl[8]*c*gamma+0.9682458365518543*bxr[4]*c*gamma-0.9682458365518543*bxl[4]*c*gamma-0.75*bxr[1]*c*gamma-0.75*bxl[1]*c*gamma+0.4330127018922193*bxr[0]*c*gamma-0.4330127018922193*bxl[0]*c*gamma+1.14564392373896*psr[8]*gamma-1.14564392373896*psl[8]*gamma-0.9682458365518543*psr[4]*gamma-0.9682458365518543*psl[4]*gamma+0.75*psr[1]*gamma-0.75*psl[1]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 
  incr[2] = (-0.5590169943749476*bxr[6]*c*gamma)+0.5590169943749476*bxl[6]*c*gamma+0.4330127018922193*bxr[3]*c*gamma+0.4330127018922193*bxl[3]*c*gamma-0.25*bxr[2]*c*gamma+0.25*bxl[2]*c*gamma+0.5590169943749476*psr[6]*gamma+0.5590169943749476*psl[6]*gamma-0.4330127018922193*psr[3]*gamma+0.4330127018922193*psl[3]*gamma+0.25*psr[2]*gamma+0.25*psl[2]*gamma; 
  incr[3] = 0.9682458365518543*bxr[6]*c*gamma-0.9682458365518543*bxl[6]*c*gamma-0.75*bxr[3]*c*gamma-0.75*bxl[3]*c*gamma+0.4330127018922193*bxr[2]*c*gamma-0.4330127018922193*bxl[2]*c*gamma-0.9682458365518543*psr[6]*gamma-0.9682458365518543*psl[6]*gamma+0.75*psr[3]*gamma-0.75*psl[3]*gamma-0.4330127018922193*psr[2]*gamma-0.4330127018922193*psl[2]*gamma; 
  incr[4] = 1.479019945774904*bxr[8]*c*gamma+1.479019945774904*bxl[8]*c*gamma-1.25*bxr[4]*c*gamma+1.25*bxl[4]*c*gamma+0.9682458365518543*bxr[1]*c*gamma+0.9682458365518543*bxl[1]*c*gamma-0.5590169943749475*bxr[0]*c*gamma+0.5590169943749475*bxl[0]*c*gamma-1.479019945774904*psr[8]*gamma+1.479019945774904*psl[8]*gamma+1.25*psr[4]*gamma+1.25*psl[4]*gamma-0.9682458365518543*psr[1]*gamma+0.9682458365518543*psl[1]*gamma+0.5590169943749475*psr[0]*gamma+0.5590169943749475*psl[0]*gamma; 
  incr[5] = 0.4330127018922194*bxr[7]*c*gamma+0.4330127018922194*bxl[7]*c*gamma-0.25*bxr[5]*c*gamma+0.25*bxl[5]*c*gamma-0.4330127018922194*psr[7]*gamma+0.4330127018922194*psl[7]*gamma+0.25*psr[5]*gamma+0.25*psl[5]*gamma; 
  incr[6] = (-1.25*bxr[6]*c*gamma)+1.25*bxl[6]*c*gamma+0.9682458365518543*bxr[3]*c*gamma+0.9682458365518543*bxl[3]*c*gamma-0.5590169943749476*bxr[2]*c*gamma+0.5590169943749476*bxl[2]*c*gamma+1.25*psr[6]*gamma+1.25*psl[6]*gamma-0.9682458365518543*psr[3]*gamma+0.9682458365518543*psl[3]*gamma+0.5590169943749476*psr[2]*gamma+0.5590169943749476*psl[2]*gamma; 
  incr[7] = (-0.75*bxr[7]*c*gamma)-0.75*bxl[7]*c*gamma+0.4330127018922194*bxr[5]*c*gamma-0.4330127018922194*bxl[5]*c*gamma+0.75*psr[7]*gamma-0.75*psl[7]*gamma-0.4330127018922194*psr[5]*gamma-0.4330127018922194*psl[5]*gamma; 
  incr[8] = (-1.75*bxr[8]*c*gamma)-1.75*bxl[8]*c*gamma+1.479019945774904*bxr[4]*c*gamma-1.479019945774904*bxl[4]*c*gamma-1.14564392373896*bxr[1]*c*gamma-1.14564392373896*bxl[1]*c*gamma+0.6614378277661477*bxr[0]*c*gamma-0.6614378277661477*bxl[0]*c*gamma+1.75*psr[8]*gamma-1.75*psl[8]*gamma-1.479019945774904*psr[4]*gamma-1.479019945774904*psl[4]*gamma+1.14564392373896*psr[1]*gamma-1.14564392373896*psl[1]*gamma-0.6614378277661477*psr[0]*gamma-0.6614378277661477*psl[0]*gamma; 
  incr[9] = (-0.25*bxr[9]*c*gamma)+0.25*bxl[9]*c*gamma+0.25*psr[9]*gamma+0.25*psl[9]*gamma; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 
  outBxr[3] += incr[3]*dx1; 
  outBxr[4] += incr[4]*dx1; 
  outBxr[5] += incr[5]*dx1; 
  outBxr[6] += incr[6]*dx1; 
  outBxr[7] += incr[7]*dx1; 
  outBxr[8] += incr[8]*dx1; 
  outBxr[9] += incr[9]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += incr[1]*dx1; 
  outBxl[2] += -1.0*incr[2]*dx1; 
  outBxl[3] += incr[3]*dx1; 
  outBxl[4] += -1.0*incr[4]*dx1; 
  outBxl[5] += -1.0*incr[5]*dx1; 
  outBxl[6] += -1.0*incr[6]*dx1; 
  outBxl[7] += incr[7]*dx1; 
  outBxl[8] += incr[8]*dx1; 
  outBxl[9] += -1.0*incr[9]*dx1; 

 
  incr[0] = 0.6614378277661477*byr[8]*c+0.6614378277661477*byl[8]*c-0.5590169943749475*byr[4]*c+0.5590169943749475*byl[4]*c+0.4330127018922193*byr[1]*c+0.4330127018922193*byl[1]*c-0.25*byr[0]*c+0.25*byl[0]*c+0.6614378277661477*ezr[8]-0.6614378277661477*ezl[8]-0.5590169943749475*ezr[4]-0.5590169943749475*ezl[4]+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*ezr[0]-0.25*ezl[0]; 
  incr[1] = (-1.14564392373896*byr[8]*c)-1.14564392373896*byl[8]*c+0.9682458365518543*byr[4]*c-0.9682458365518543*byl[4]*c-0.75*byr[1]*c-0.75*byl[1]*c+0.4330127018922193*byr[0]*c-0.4330127018922193*byl[0]*c-1.14564392373896*ezr[8]+1.14564392373896*ezl[8]+0.9682458365518543*ezr[4]+0.9682458365518543*ezl[4]-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*ezr[0]+0.4330127018922193*ezl[0]; 
  incr[2] = (-0.5590169943749476*byr[6]*c)+0.5590169943749476*byl[6]*c+0.4330127018922193*byr[3]*c+0.4330127018922193*byl[3]*c-0.25*byr[2]*c+0.25*byl[2]*c-0.5590169943749476*ezr[6]-0.5590169943749476*ezl[6]+0.4330127018922193*ezr[3]-0.4330127018922193*ezl[3]-0.25*ezr[2]-0.25*ezl[2]; 
  incr[3] = 0.9682458365518543*byr[6]*c-0.9682458365518543*byl[6]*c-0.75*byr[3]*c-0.75*byl[3]*c+0.4330127018922193*byr[2]*c-0.4330127018922193*byl[2]*c+0.9682458365518543*ezr[6]+0.9682458365518543*ezl[6]-0.75*ezr[3]+0.75*ezl[3]+0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]; 
  incr[4] = 1.479019945774904*byr[8]*c+1.479019945774904*byl[8]*c-1.25*byr[4]*c+1.25*byl[4]*c+0.9682458365518543*byr[1]*c+0.9682458365518543*byl[1]*c-0.5590169943749475*byr[0]*c+0.5590169943749475*byl[0]*c+1.479019945774904*ezr[8]-1.479019945774904*ezl[8]-1.25*ezr[4]-1.25*ezl[4]+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*ezr[0]-0.5590169943749475*ezl[0]; 
  incr[5] = 0.4330127018922194*byr[7]*c+0.4330127018922194*byl[7]*c-0.25*byr[5]*c+0.25*byl[5]*c+0.4330127018922194*ezr[7]-0.4330127018922194*ezl[7]-0.25*ezr[5]-0.25*ezl[5]; 
  incr[6] = (-1.25*byr[6]*c)+1.25*byl[6]*c+0.9682458365518543*byr[3]*c+0.9682458365518543*byl[3]*c-0.5590169943749476*byr[2]*c+0.5590169943749476*byl[2]*c-1.25*ezr[6]-1.25*ezl[6]+0.9682458365518543*ezr[3]-0.9682458365518543*ezl[3]-0.5590169943749476*ezr[2]-0.5590169943749476*ezl[2]; 
  incr[7] = (-0.75*byr[7]*c)-0.75*byl[7]*c+0.4330127018922194*byr[5]*c-0.4330127018922194*byl[5]*c-0.75*ezr[7]+0.75*ezl[7]+0.4330127018922194*ezr[5]+0.4330127018922194*ezl[5]; 
  incr[8] = (-1.75*byr[8]*c)-1.75*byl[8]*c+1.479019945774904*byr[4]*c-1.479019945774904*byl[4]*c-1.14564392373896*byr[1]*c-1.14564392373896*byl[1]*c+0.6614378277661477*byr[0]*c-0.6614378277661477*byl[0]*c-1.75*ezr[8]+1.75*ezl[8]+1.479019945774904*ezr[4]+1.479019945774904*ezl[4]-1.14564392373896*ezr[1]+1.14564392373896*ezl[1]+0.6614378277661477*ezr[0]+0.6614378277661477*ezl[0]; 
  incr[9] = (-0.25*byr[9]*c)+0.25*byl[9]*c-0.25*ezr[9]-0.25*ezl[9]; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 
  outByr[3] += incr[3]*dx1; 
  outByr[4] += incr[4]*dx1; 
  outByr[5] += incr[5]*dx1; 
  outByr[6] += incr[6]*dx1; 
  outByr[7] += incr[7]*dx1; 
  outByr[8] += incr[8]*dx1; 
  outByr[9] += incr[9]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += incr[1]*dx1; 
  outByl[2] += -1.0*incr[2]*dx1; 
  outByl[3] += incr[3]*dx1; 
  outByl[4] += -1.0*incr[4]*dx1; 
  outByl[5] += -1.0*incr[5]*dx1; 
  outByl[6] += -1.0*incr[6]*dx1; 
  outByl[7] += incr[7]*dx1; 
  outByl[8] += incr[8]*dx1; 
  outByl[9] += -1.0*incr[9]*dx1; 

 
  incr[0] = 0.6614378277661477*bzr[8]*c+0.6614378277661477*bzl[8]*c-0.5590169943749475*bzr[4]*c+0.5590169943749475*bzl[4]*c+0.4330127018922193*bzr[1]*c+0.4330127018922193*bzl[1]*c-0.25*bzr[0]*c+0.25*bzl[0]*c-0.6614378277661477*eyr[8]+0.6614378277661477*eyl[8]+0.5590169943749475*eyr[4]+0.5590169943749475*eyl[4]-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*eyr[0]+0.25*eyl[0]; 
  incr[1] = (-1.14564392373896*bzr[8]*c)-1.14564392373896*bzl[8]*c+0.9682458365518543*bzr[4]*c-0.9682458365518543*bzl[4]*c-0.75*bzr[1]*c-0.75*bzl[1]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c+1.14564392373896*eyr[8]-1.14564392373896*eyl[8]-0.9682458365518543*eyr[4]-0.9682458365518543*eyl[4]+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0]; 
  incr[2] = (-0.5590169943749476*bzr[6]*c)+0.5590169943749476*bzl[6]*c+0.4330127018922193*bzr[3]*c+0.4330127018922193*bzl[3]*c-0.25*bzr[2]*c+0.25*bzl[2]*c+0.5590169943749476*eyr[6]+0.5590169943749476*eyl[6]-0.4330127018922193*eyr[3]+0.4330127018922193*eyl[3]+0.25*eyr[2]+0.25*eyl[2]; 
  incr[3] = 0.9682458365518543*bzr[6]*c-0.9682458365518543*bzl[6]*c-0.75*bzr[3]*c-0.75*bzl[3]*c+0.4330127018922193*bzr[2]*c-0.4330127018922193*bzl[2]*c-0.9682458365518543*eyr[6]-0.9682458365518543*eyl[6]+0.75*eyr[3]-0.75*eyl[3]-0.4330127018922193*eyr[2]-0.4330127018922193*eyl[2]; 
  incr[4] = 1.479019945774904*bzr[8]*c+1.479019945774904*bzl[8]*c-1.25*bzr[4]*c+1.25*bzl[4]*c+0.9682458365518543*bzr[1]*c+0.9682458365518543*bzl[1]*c-0.5590169943749475*bzr[0]*c+0.5590169943749475*bzl[0]*c-1.479019945774904*eyr[8]+1.479019945774904*eyl[8]+1.25*eyr[4]+1.25*eyl[4]-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0]; 
  incr[5] = 0.4330127018922194*bzr[7]*c+0.4330127018922194*bzl[7]*c-0.25*bzr[5]*c+0.25*bzl[5]*c-0.4330127018922194*eyr[7]+0.4330127018922194*eyl[7]+0.25*eyr[5]+0.25*eyl[5]; 
  incr[6] = (-1.25*bzr[6]*c)+1.25*bzl[6]*c+0.9682458365518543*bzr[3]*c+0.9682458365518543*bzl[3]*c-0.5590169943749476*bzr[2]*c+0.5590169943749476*bzl[2]*c+1.25*eyr[6]+1.25*eyl[6]-0.9682458365518543*eyr[3]+0.9682458365518543*eyl[3]+0.5590169943749476*eyr[2]+0.5590169943749476*eyl[2]; 
  incr[7] = (-0.75*bzr[7]*c)-0.75*bzl[7]*c+0.4330127018922194*bzr[5]*c-0.4330127018922194*bzl[5]*c+0.75*eyr[7]-0.75*eyl[7]-0.4330127018922194*eyr[5]-0.4330127018922194*eyl[5]; 
  incr[8] = (-1.75*bzr[8]*c)-1.75*bzl[8]*c+1.479019945774904*bzr[4]*c-1.479019945774904*bzl[4]*c-1.14564392373896*bzr[1]*c-1.14564392373896*bzl[1]*c+0.6614378277661477*bzr[0]*c-0.6614378277661477*bzl[0]*c+1.75*eyr[8]-1.75*eyl[8]-1.479019945774904*eyr[4]-1.479019945774904*eyl[4]+1.14564392373896*eyr[1]-1.14564392373896*eyl[1]-0.6614378277661477*eyr[0]-0.6614378277661477*eyl[0]; 
  incr[9] = (-0.25*bzr[9]*c)+0.25*bzl[9]*c+0.25*eyr[9]+0.25*eyl[9]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 
  outBzr[3] += incr[3]*dx1; 
  outBzr[4] += incr[4]*dx1; 
  outBzr[5] += incr[5]*dx1; 
  outBzr[6] += incr[6]*dx1; 
  outBzr[7] += incr[7]*dx1; 
  outBzr[8] += incr[8]*dx1; 
  outBzr[9] += incr[9]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += incr[1]*dx1; 
  outBzl[2] += -1.0*incr[2]*dx1; 
  outBzl[3] += incr[3]*dx1; 
  outBzl[4] += -1.0*incr[4]*dx1; 
  outBzl[5] += -1.0*incr[5]*dx1; 
  outBzl[6] += -1.0*incr[6]*dx1; 
  outBzl[7] += incr[7]*dx1; 
  outBzl[8] += incr[8]*dx1; 
  outBzl[9] += -1.0*incr[9]*dx1; 

 
  incr[0] = 0.6614378277661477*phr[8]*c*chi+0.6614378277661477*phl[8]*c*chi-0.5590169943749475*phr[4]*c*chi+0.5590169943749475*phl[4]*c*chi+0.4330127018922193*phr[1]*c*chi+0.4330127018922193*phl[1]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi-0.6614378277661477*exr[8]*chi+0.6614378277661477*exl[8]*chi+0.5590169943749475*exr[4]*chi+0.5590169943749475*exl[4]*chi-0.4330127018922193*exr[1]*chi+0.4330127018922193*exl[1]*chi+0.25*exr[0]*chi+0.25*exl[0]*chi; 
  incr[1] = (-1.14564392373896*phr[8]*c*chi)-1.14564392373896*phl[8]*c*chi+0.9682458365518543*phr[4]*c*chi-0.9682458365518543*phl[4]*c*chi-0.75*phr[1]*c*chi-0.75*phl[1]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi+1.14564392373896*exr[8]*chi-1.14564392373896*exl[8]*chi-0.9682458365518543*exr[4]*chi-0.9682458365518543*exl[4]*chi+0.75*exr[1]*chi-0.75*exl[1]*chi-0.4330127018922193*exr[0]*chi-0.4330127018922193*exl[0]*chi; 
  incr[2] = (-0.5590169943749476*phr[6]*c*chi)+0.5590169943749476*phl[6]*c*chi+0.4330127018922193*phr[3]*c*chi+0.4330127018922193*phl[3]*c*chi-0.25*phr[2]*c*chi+0.25*phl[2]*c*chi+0.5590169943749476*exr[6]*chi+0.5590169943749476*exl[6]*chi-0.4330127018922193*exr[3]*chi+0.4330127018922193*exl[3]*chi+0.25*exr[2]*chi+0.25*exl[2]*chi; 
  incr[3] = 0.9682458365518543*phr[6]*c*chi-0.9682458365518543*phl[6]*c*chi-0.75*phr[3]*c*chi-0.75*phl[3]*c*chi+0.4330127018922193*phr[2]*c*chi-0.4330127018922193*phl[2]*c*chi-0.9682458365518543*exr[6]*chi-0.9682458365518543*exl[6]*chi+0.75*exr[3]*chi-0.75*exl[3]*chi-0.4330127018922193*exr[2]*chi-0.4330127018922193*exl[2]*chi; 
  incr[4] = 1.479019945774904*phr[8]*c*chi+1.479019945774904*phl[8]*c*chi-1.25*phr[4]*c*chi+1.25*phl[4]*c*chi+0.9682458365518543*phr[1]*c*chi+0.9682458365518543*phl[1]*c*chi-0.5590169943749475*phr[0]*c*chi+0.5590169943749475*phl[0]*c*chi-1.479019945774904*exr[8]*chi+1.479019945774904*exl[8]*chi+1.25*exr[4]*chi+1.25*exl[4]*chi-0.9682458365518543*exr[1]*chi+0.9682458365518543*exl[1]*chi+0.5590169943749475*exr[0]*chi+0.5590169943749475*exl[0]*chi; 
  incr[5] = 0.4330127018922194*phr[7]*c*chi+0.4330127018922194*phl[7]*c*chi-0.25*phr[5]*c*chi+0.25*phl[5]*c*chi-0.4330127018922194*exr[7]*chi+0.4330127018922194*exl[7]*chi+0.25*exr[5]*chi+0.25*exl[5]*chi; 
  incr[6] = (-1.25*phr[6]*c*chi)+1.25*phl[6]*c*chi+0.9682458365518543*phr[3]*c*chi+0.9682458365518543*phl[3]*c*chi-0.5590169943749476*phr[2]*c*chi+0.5590169943749476*phl[2]*c*chi+1.25*exr[6]*chi+1.25*exl[6]*chi-0.9682458365518543*exr[3]*chi+0.9682458365518543*exl[3]*chi+0.5590169943749476*exr[2]*chi+0.5590169943749476*exl[2]*chi; 
  incr[7] = (-0.75*phr[7]*c*chi)-0.75*phl[7]*c*chi+0.4330127018922194*phr[5]*c*chi-0.4330127018922194*phl[5]*c*chi+0.75*exr[7]*chi-0.75*exl[7]*chi-0.4330127018922194*exr[5]*chi-0.4330127018922194*exl[5]*chi; 
  incr[8] = (-1.75*phr[8]*c*chi)-1.75*phl[8]*c*chi+1.479019945774904*phr[4]*c*chi-1.479019945774904*phl[4]*c*chi-1.14564392373896*phr[1]*c*chi-1.14564392373896*phl[1]*c*chi+0.6614378277661477*phr[0]*c*chi-0.6614378277661477*phl[0]*c*chi+1.75*exr[8]*chi-1.75*exl[8]*chi-1.479019945774904*exr[4]*chi-1.479019945774904*exl[4]*chi+1.14564392373896*exr[1]*chi-1.14564392373896*exl[1]*chi-0.6614378277661477*exr[0]*chi-0.6614378277661477*exl[0]*chi; 
  incr[9] = (-0.25*phr[9]*c*chi)+0.25*phl[9]*c*chi+0.25*exr[9]*chi+0.25*exl[9]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 
  outPhr[3] += incr[3]*dx1; 
  outPhr[4] += incr[4]*dx1; 
  outPhr[5] += incr[5]*dx1; 
  outPhr[6] += incr[6]*dx1; 
  outPhr[7] += incr[7]*dx1; 
  outPhr[8] += incr[8]*dx1; 
  outPhr[9] += incr[9]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += incr[1]*dx1; 
  outPhl[2] += -1.0*incr[2]*dx1; 
  outPhl[3] += incr[3]*dx1; 
  outPhl[4] += -1.0*incr[4]*dx1; 
  outPhl[5] += -1.0*incr[5]*dx1; 
  outPhl[6] += -1.0*incr[6]*dx1; 
  outPhl[7] += incr[7]*dx1; 
  outPhl[8] += incr[8]*dx1; 
  outPhl[9] += -1.0*incr[9]*dx1; 

 
  incr[0] = (-0.6614378277661477*bxr[8]*c2*gamma)+0.6614378277661477*bxl[8]*c2*gamma+0.5590169943749475*bxr[4]*c2*gamma+0.5590169943749475*bxl[4]*c2*gamma-0.4330127018922193*bxr[1]*c2*gamma+0.4330127018922193*bxl[1]*c2*gamma+0.25*bxr[0]*c2*gamma+0.25*bxl[0]*c2*gamma+0.6614378277661477*psr[8]*c*gamma+0.6614378277661477*psl[8]*c*gamma-0.5590169943749475*psr[4]*c*gamma+0.5590169943749475*psl[4]*c*gamma+0.4330127018922193*psr[1]*c*gamma+0.4330127018922193*psl[1]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = 1.14564392373896*bxr[8]*c2*gamma-1.14564392373896*bxl[8]*c2*gamma-0.9682458365518543*bxr[4]*c2*gamma-0.9682458365518543*bxl[4]*c2*gamma+0.75*bxr[1]*c2*gamma-0.75*bxl[1]*c2*gamma-0.4330127018922193*bxr[0]*c2*gamma-0.4330127018922193*bxl[0]*c2*gamma-1.14564392373896*psr[8]*c*gamma-1.14564392373896*psl[8]*c*gamma+0.9682458365518543*psr[4]*c*gamma-0.9682458365518543*psl[4]*c*gamma-0.75*psr[1]*c*gamma-0.75*psl[1]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 
  incr[2] = 0.5590169943749476*bxr[6]*c2*gamma+0.5590169943749476*bxl[6]*c2*gamma-0.4330127018922193*bxr[3]*c2*gamma+0.4330127018922193*bxl[3]*c2*gamma+0.25*bxr[2]*c2*gamma+0.25*bxl[2]*c2*gamma-0.5590169943749476*psr[6]*c*gamma+0.5590169943749476*psl[6]*c*gamma+0.4330127018922193*psr[3]*c*gamma+0.4330127018922193*psl[3]*c*gamma-0.25*psr[2]*c*gamma+0.25*psl[2]*c*gamma; 
  incr[3] = (-0.9682458365518543*bxr[6]*c2*gamma)-0.9682458365518543*bxl[6]*c2*gamma+0.75*bxr[3]*c2*gamma-0.75*bxl[3]*c2*gamma-0.4330127018922193*bxr[2]*c2*gamma-0.4330127018922193*bxl[2]*c2*gamma+0.9682458365518543*psr[6]*c*gamma-0.9682458365518543*psl[6]*c*gamma-0.75*psr[3]*c*gamma-0.75*psl[3]*c*gamma+0.4330127018922193*psr[2]*c*gamma-0.4330127018922193*psl[2]*c*gamma; 
  incr[4] = (-1.479019945774904*bxr[8]*c2*gamma)+1.479019945774904*bxl[8]*c2*gamma+1.25*bxr[4]*c2*gamma+1.25*bxl[4]*c2*gamma-0.9682458365518543*bxr[1]*c2*gamma+0.9682458365518543*bxl[1]*c2*gamma+0.5590169943749475*bxr[0]*c2*gamma+0.5590169943749475*bxl[0]*c2*gamma+1.479019945774904*psr[8]*c*gamma+1.479019945774904*psl[8]*c*gamma-1.25*psr[4]*c*gamma+1.25*psl[4]*c*gamma+0.9682458365518543*psr[1]*c*gamma+0.9682458365518543*psl[1]*c*gamma-0.5590169943749475*psr[0]*c*gamma+0.5590169943749475*psl[0]*c*gamma; 
  incr[5] = (-0.4330127018922194*bxr[7]*c2*gamma)+0.4330127018922194*bxl[7]*c2*gamma+0.25*bxr[5]*c2*gamma+0.25*bxl[5]*c2*gamma+0.4330127018922194*psr[7]*c*gamma+0.4330127018922194*psl[7]*c*gamma-0.25*psr[5]*c*gamma+0.25*psl[5]*c*gamma; 
  incr[6] = 1.25*bxr[6]*c2*gamma+1.25*bxl[6]*c2*gamma-0.9682458365518543*bxr[3]*c2*gamma+0.9682458365518543*bxl[3]*c2*gamma+0.5590169943749476*bxr[2]*c2*gamma+0.5590169943749476*bxl[2]*c2*gamma-1.25*psr[6]*c*gamma+1.25*psl[6]*c*gamma+0.9682458365518543*psr[3]*c*gamma+0.9682458365518543*psl[3]*c*gamma-0.5590169943749476*psr[2]*c*gamma+0.5590169943749476*psl[2]*c*gamma; 
  incr[7] = 0.75*bxr[7]*c2*gamma-0.75*bxl[7]*c2*gamma-0.4330127018922194*bxr[5]*c2*gamma-0.4330127018922194*bxl[5]*c2*gamma-0.75*psr[7]*c*gamma-0.75*psl[7]*c*gamma+0.4330127018922194*psr[5]*c*gamma-0.4330127018922194*psl[5]*c*gamma; 
  incr[8] = 1.75*bxr[8]*c2*gamma-1.75*bxl[8]*c2*gamma-1.479019945774904*bxr[4]*c2*gamma-1.479019945774904*bxl[4]*c2*gamma+1.14564392373896*bxr[1]*c2*gamma-1.14564392373896*bxl[1]*c2*gamma-0.6614378277661477*bxr[0]*c2*gamma-0.6614378277661477*bxl[0]*c2*gamma-1.75*psr[8]*c*gamma-1.75*psl[8]*c*gamma+1.479019945774904*psr[4]*c*gamma-1.479019945774904*psl[4]*c*gamma-1.14564392373896*psr[1]*c*gamma-1.14564392373896*psl[1]*c*gamma+0.6614378277661477*psr[0]*c*gamma-0.6614378277661477*psl[0]*c*gamma; 
  incr[9] = 0.25*bxr[9]*c2*gamma+0.25*bxl[9]*c2*gamma-0.25*psr[9]*c*gamma+0.25*psl[9]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 
  outPsr[3] += incr[3]*dx1; 
  outPsr[4] += incr[4]*dx1; 
  outPsr[5] += incr[5]*dx1; 
  outPsr[6] += incr[6]*dx1; 
  outPsr[7] += incr[7]*dx1; 
  outPsr[8] += incr[8]*dx1; 
  outPsr[9] += incr[9]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += incr[1]*dx1; 
  outPsl[2] += -1.0*incr[2]*dx1; 
  outPsl[3] += incr[3]*dx1; 
  outPsl[4] += -1.0*incr[4]*dx1; 
  outPsl[5] += -1.0*incr[5]*dx1; 
  outPsl[6] += -1.0*incr[6]*dx1; 
  outPsl[7] += incr[7]*dx1; 
  outPsl[8] += incr[8]*dx1; 
  outPsl[9] += -1.0*incr[9]*dx1; 

 
} 
void MaxwellSurf2xMax_X_P4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[15]; 
  const double *ezl = &ql[30]; 
  const double *bxl = &ql[45]; 
  const double *byl = &ql[60]; 
  const double *bzl = &ql[75]; 
  const double *phl = &ql[90]; 
  const double *psl = &ql[105]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[15]; 
  double *outEzl = &outl[30]; 
  double *outBxl = &outl[45]; 
  double *outByl = &outl[60]; 
  double *outBzl = &outl[75]; 
  double *outPhl = &outl[90]; 
  double *outPsl = &outl[105]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[15]; 
  const double *ezr = &qr[30]; 
  const double *bxr = &qr[45]; 
  const double *byr = &qr[60]; 
  const double *bzr = &qr[75]; 
  const double *phr = &qr[90]; 
  const double *psr = &qr[105]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[15]; 
  double *outEzr = &outr[30]; 
  double *outBxr = &outr[45]; 
  double *outByr = &outr[60]; 
  double *outBzr = &outr[75]; 
  double *outPhr = &outr[90]; 
  double *outPsr = &outr[105]; 
 
  double incr[15]; 
 
  incr[0] = 0.75*phr[13]*c2*chi+0.75*phl[13]*c2*chi-0.6614378277661477*phr[8]*c2*chi+0.6614378277661477*phl[8]*c2*chi+0.5590169943749475*phr[4]*c2*chi+0.5590169943749475*phl[4]*c2*chi-0.4330127018922193*phr[1]*c2*chi+0.4330127018922193*phl[1]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi-0.75*exr[13]*c*chi+0.75*exl[13]*c*chi+0.6614378277661477*exr[8]*c*chi+0.6614378277661477*exl[8]*c*chi-0.5590169943749475*exr[4]*c*chi+0.5590169943749475*exl[4]*c*chi+0.4330127018922193*exr[1]*c*chi+0.4330127018922193*exl[1]*c*chi-0.25*exr[0]*c*chi+0.25*exl[0]*c*chi; 
  incr[1] = (-1.299038105676658*phr[13]*c2*chi)-1.299038105676658*phl[13]*c2*chi+1.14564392373896*phr[8]*c2*chi-1.14564392373896*phl[8]*c2*chi-0.9682458365518543*phr[4]*c2*chi-0.9682458365518543*phl[4]*c2*chi+0.75*phr[1]*c2*chi-0.75*phl[1]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi+1.299038105676658*exr[13]*c*chi-1.299038105676658*exl[13]*c*chi-1.14564392373896*exr[8]*c*chi-1.14564392373896*exl[8]*c*chi+0.9682458365518543*exr[4]*c*chi-0.9682458365518543*exl[4]*c*chi-0.75*exr[1]*c*chi-0.75*exl[1]*c*chi+0.4330127018922193*exr[0]*c*chi-0.4330127018922193*exl[0]*c*chi; 
  incr[2] = (-0.6614378277661477*phr[11]*c2*chi)+0.6614378277661477*phl[11]*c2*chi+0.5590169943749476*phr[6]*c2*chi+0.5590169943749476*phl[6]*c2*chi-0.4330127018922193*phr[3]*c2*chi+0.4330127018922193*phl[3]*c2*chi+0.25*phr[2]*c2*chi+0.25*phl[2]*c2*chi+0.6614378277661477*exr[11]*c*chi+0.6614378277661477*exl[11]*c*chi-0.5590169943749476*exr[6]*c*chi+0.5590169943749476*exl[6]*c*chi+0.4330127018922193*exr[3]*c*chi+0.4330127018922193*exl[3]*c*chi-0.25*exr[2]*c*chi+0.25*exl[2]*c*chi; 
  incr[3] = 1.14564392373896*phr[11]*c2*chi-1.14564392373896*phl[11]*c2*chi-0.9682458365518543*phr[6]*c2*chi-0.9682458365518543*phl[6]*c2*chi+0.75*phr[3]*c2*chi-0.75*phl[3]*c2*chi-0.4330127018922193*phr[2]*c2*chi-0.4330127018922193*phl[2]*c2*chi-1.14564392373896*exr[11]*c*chi-1.14564392373896*exl[11]*c*chi+0.9682458365518543*exr[6]*c*chi-0.9682458365518543*exl[6]*c*chi-0.75*exr[3]*c*chi-0.75*exl[3]*c*chi+0.4330127018922193*exr[2]*c*chi-0.4330127018922193*exl[2]*c*chi; 
  incr[4] = 1.677050983124842*phr[13]*c2*chi+1.677050983124842*phl[13]*c2*chi-1.479019945774904*phr[8]*c2*chi+1.479019945774904*phl[8]*c2*chi+1.25*phr[4]*c2*chi+1.25*phl[4]*c2*chi-0.9682458365518543*phr[1]*c2*chi+0.9682458365518543*phl[1]*c2*chi+0.5590169943749475*phr[0]*c2*chi+0.5590169943749475*phl[0]*c2*chi-1.677050983124842*exr[13]*c*chi+1.677050983124842*exl[13]*c*chi+1.479019945774904*exr[8]*c*chi+1.479019945774904*exl[8]*c*chi-1.25*exr[4]*c*chi+1.25*exl[4]*c*chi+0.9682458365518543*exr[1]*c*chi+0.9682458365518543*exl[1]*c*chi-0.5590169943749475*exr[0]*c*chi+0.5590169943749475*exl[0]*c*chi; 
  incr[5] = 0.5590169943749475*phr[10]*c2*chi+0.5590169943749475*phl[10]*c2*chi-0.4330127018922194*phr[7]*c2*chi+0.4330127018922194*phl[7]*c2*chi+0.25*phr[5]*c2*chi+0.25*phl[5]*c2*chi-0.5590169943749475*exr[10]*c*chi+0.5590169943749475*exl[10]*c*chi+0.4330127018922194*exr[7]*c*chi+0.4330127018922194*exl[7]*c*chi-0.25*exr[5]*c*chi+0.25*exl[5]*c*chi; 
  incr[6] = (-1.479019945774904*phr[11]*c2*chi)+1.479019945774904*phl[11]*c2*chi+1.25*phr[6]*c2*chi+1.25*phl[6]*c2*chi-0.9682458365518543*phr[3]*c2*chi+0.9682458365518543*phl[3]*c2*chi+0.5590169943749476*phr[2]*c2*chi+0.5590169943749476*phl[2]*c2*chi+1.479019945774904*exr[11]*c*chi+1.479019945774904*exl[11]*c*chi-1.25*exr[6]*c*chi+1.25*exl[6]*c*chi+0.9682458365518543*exr[3]*c*chi+0.9682458365518543*exl[3]*c*chi-0.5590169943749476*exr[2]*c*chi+0.5590169943749476*exl[2]*c*chi; 
  incr[7] = (-0.9682458365518543*phr[10]*c2*chi)-0.9682458365518543*phl[10]*c2*chi+0.75*phr[7]*c2*chi-0.75*phl[7]*c2*chi-0.4330127018922194*phr[5]*c2*chi-0.4330127018922194*phl[5]*c2*chi+0.9682458365518543*exr[10]*c*chi-0.9682458365518543*exl[10]*c*chi-0.75*exr[7]*c*chi-0.75*exl[7]*c*chi+0.4330127018922194*exr[5]*c*chi-0.4330127018922194*exl[5]*c*chi; 
  incr[8] = (-1.984313483298443*phr[13]*c2*chi)-1.984313483298443*phl[13]*c2*chi+1.75*phr[8]*c2*chi-1.75*phl[8]*c2*chi-1.479019945774904*phr[4]*c2*chi-1.479019945774904*phl[4]*c2*chi+1.14564392373896*phr[1]*c2*chi-1.14564392373896*phl[1]*c2*chi-0.6614378277661477*phr[0]*c2*chi-0.6614378277661477*phl[0]*c2*chi+1.984313483298443*exr[13]*c*chi-1.984313483298443*exl[13]*c*chi-1.75*exr[8]*c*chi-1.75*exl[8]*c*chi+1.479019945774904*exr[4]*c*chi-1.479019945774904*exl[4]*c*chi-1.14564392373896*exr[1]*c*chi-1.14564392373896*exl[1]*c*chi+0.6614378277661477*exr[0]*c*chi-0.6614378277661477*exl[0]*c*chi; 
  incr[9] = (-0.4330127018922193*phr[12]*c2*chi)+0.4330127018922193*phl[12]*c2*chi+0.25*phr[9]*c2*chi+0.25*phl[9]*c2*chi+0.4330127018922193*exr[12]*c*chi+0.4330127018922193*exl[12]*c*chi-0.25*exr[9]*c*chi+0.25*exl[9]*c*chi; 
  incr[10] = 1.25*phr[10]*c2*chi+1.25*phl[10]*c2*chi-0.9682458365518543*phr[7]*c2*chi+0.9682458365518543*phl[7]*c2*chi+0.5590169943749475*phr[5]*c2*chi+0.5590169943749475*phl[5]*c2*chi-1.25*exr[10]*c*chi+1.25*exl[10]*c*chi+0.9682458365518543*exr[7]*c*chi+0.9682458365518543*exl[7]*c*chi-0.5590169943749475*exr[5]*c*chi+0.5590169943749475*exl[5]*c*chi; 
  incr[11] = 1.75*phr[11]*c2*chi-1.75*phl[11]*c2*chi-1.479019945774904*phr[6]*c2*chi-1.479019945774904*phl[6]*c2*chi+1.14564392373896*phr[3]*c2*chi-1.14564392373896*phl[3]*c2*chi-0.6614378277661477*phr[2]*c2*chi-0.6614378277661477*phl[2]*c2*chi-1.75*exr[11]*c*chi-1.75*exl[11]*c*chi+1.479019945774904*exr[6]*c*chi-1.479019945774904*exl[6]*c*chi-1.14564392373896*exr[3]*c*chi-1.14564392373896*exl[3]*c*chi+0.6614378277661477*exr[2]*c*chi-0.6614378277661477*exl[2]*c*chi; 
  incr[12] = 0.75*phr[12]*c2*chi-0.75*phl[12]*c2*chi-0.4330127018922193*phr[9]*c2*chi-0.4330127018922193*phl[9]*c2*chi-0.75*exr[12]*c*chi-0.75*exl[12]*c*chi+0.4330127018922193*exr[9]*c*chi-0.4330127018922193*exl[9]*c*chi; 
  incr[13] = 2.25*phr[13]*c2*chi+2.25*phl[13]*c2*chi-1.984313483298443*phr[8]*c2*chi+1.984313483298443*phl[8]*c2*chi+1.677050983124842*phr[4]*c2*chi+1.677050983124842*phl[4]*c2*chi-1.299038105676658*phr[1]*c2*chi+1.299038105676658*phl[1]*c2*chi+0.75*phr[0]*c2*chi+0.75*phl[0]*c2*chi-2.25*exr[13]*c*chi+2.25*exl[13]*c*chi+1.984313483298443*exr[8]*c*chi+1.984313483298443*exl[8]*c*chi-1.677050983124842*exr[4]*c*chi+1.677050983124842*exl[4]*c*chi+1.299038105676658*exr[1]*c*chi+1.299038105676658*exl[1]*c*chi-0.75*exr[0]*c*chi+0.75*exl[0]*c*chi; 
  incr[14] = 0.25*phr[14]*c2*chi+0.25*phl[14]*c2*chi-0.25*exr[14]*c*chi+0.25*exl[14]*c*chi; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 
  outExr[3] += incr[3]*dx1; 
  outExr[4] += incr[4]*dx1; 
  outExr[5] += incr[5]*dx1; 
  outExr[6] += incr[6]*dx1; 
  outExr[7] += incr[7]*dx1; 
  outExr[8] += incr[8]*dx1; 
  outExr[9] += incr[9]*dx1; 
  outExr[10] += incr[10]*dx1; 
  outExr[11] += incr[11]*dx1; 
  outExr[12] += incr[12]*dx1; 
  outExr[13] += incr[13]*dx1; 
  outExr[14] += incr[14]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += incr[1]*dx1; 
  outExl[2] += -1.0*incr[2]*dx1; 
  outExl[3] += incr[3]*dx1; 
  outExl[4] += -1.0*incr[4]*dx1; 
  outExl[5] += -1.0*incr[5]*dx1; 
  outExl[6] += -1.0*incr[6]*dx1; 
  outExl[7] += incr[7]*dx1; 
  outExl[8] += incr[8]*dx1; 
  outExl[9] += -1.0*incr[9]*dx1; 
  outExl[10] += -1.0*incr[10]*dx1; 
  outExl[11] += incr[11]*dx1; 
  outExl[12] += incr[12]*dx1; 
  outExl[13] += -1.0*incr[13]*dx1; 
  outExl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = 0.75*bzr[13]*c2+0.75*bzl[13]*c2-0.6614378277661477*bzr[8]*c2+0.6614378277661477*bzl[8]*c2+0.5590169943749475*bzr[4]*c2+0.5590169943749475*bzl[4]*c2-0.4330127018922193*bzr[1]*c2+0.4330127018922193*bzl[1]*c2+0.25*bzr[0]*c2+0.25*bzl[0]*c2-0.75*eyr[13]*c+0.75*eyl[13]*c+0.6614378277661477*eyr[8]*c+0.6614378277661477*eyl[8]*c-0.5590169943749475*eyr[4]*c+0.5590169943749475*eyl[4]*c+0.4330127018922193*eyr[1]*c+0.4330127018922193*eyl[1]*c-0.25*eyr[0]*c+0.25*eyl[0]*c; 
  incr[1] = (-1.299038105676658*bzr[13]*c2)-1.299038105676658*bzl[13]*c2+1.14564392373896*bzr[8]*c2-1.14564392373896*bzl[8]*c2-0.9682458365518543*bzr[4]*c2-0.9682458365518543*bzl[4]*c2+0.75*bzr[1]*c2-0.75*bzl[1]*c2-0.4330127018922193*bzr[0]*c2-0.4330127018922193*bzl[0]*c2+1.299038105676658*eyr[13]*c-1.299038105676658*eyl[13]*c-1.14564392373896*eyr[8]*c-1.14564392373896*eyl[8]*c+0.9682458365518543*eyr[4]*c-0.9682458365518543*eyl[4]*c-0.75*eyr[1]*c-0.75*eyl[1]*c+0.4330127018922193*eyr[0]*c-0.4330127018922193*eyl[0]*c; 
  incr[2] = (-0.6614378277661477*bzr[11]*c2)+0.6614378277661477*bzl[11]*c2+0.5590169943749476*bzr[6]*c2+0.5590169943749476*bzl[6]*c2-0.4330127018922193*bzr[3]*c2+0.4330127018922193*bzl[3]*c2+0.25*bzr[2]*c2+0.25*bzl[2]*c2+0.6614378277661477*eyr[11]*c+0.6614378277661477*eyl[11]*c-0.5590169943749476*eyr[6]*c+0.5590169943749476*eyl[6]*c+0.4330127018922193*eyr[3]*c+0.4330127018922193*eyl[3]*c-0.25*eyr[2]*c+0.25*eyl[2]*c; 
  incr[3] = 1.14564392373896*bzr[11]*c2-1.14564392373896*bzl[11]*c2-0.9682458365518543*bzr[6]*c2-0.9682458365518543*bzl[6]*c2+0.75*bzr[3]*c2-0.75*bzl[3]*c2-0.4330127018922193*bzr[2]*c2-0.4330127018922193*bzl[2]*c2-1.14564392373896*eyr[11]*c-1.14564392373896*eyl[11]*c+0.9682458365518543*eyr[6]*c-0.9682458365518543*eyl[6]*c-0.75*eyr[3]*c-0.75*eyl[3]*c+0.4330127018922193*eyr[2]*c-0.4330127018922193*eyl[2]*c; 
  incr[4] = 1.677050983124842*bzr[13]*c2+1.677050983124842*bzl[13]*c2-1.479019945774904*bzr[8]*c2+1.479019945774904*bzl[8]*c2+1.25*bzr[4]*c2+1.25*bzl[4]*c2-0.9682458365518543*bzr[1]*c2+0.9682458365518543*bzl[1]*c2+0.5590169943749475*bzr[0]*c2+0.5590169943749475*bzl[0]*c2-1.677050983124842*eyr[13]*c+1.677050983124842*eyl[13]*c+1.479019945774904*eyr[8]*c+1.479019945774904*eyl[8]*c-1.25*eyr[4]*c+1.25*eyl[4]*c+0.9682458365518543*eyr[1]*c+0.9682458365518543*eyl[1]*c-0.5590169943749475*eyr[0]*c+0.5590169943749475*eyl[0]*c; 
  incr[5] = 0.5590169943749475*bzr[10]*c2+0.5590169943749475*bzl[10]*c2-0.4330127018922194*bzr[7]*c2+0.4330127018922194*bzl[7]*c2+0.25*bzr[5]*c2+0.25*bzl[5]*c2-0.5590169943749475*eyr[10]*c+0.5590169943749475*eyl[10]*c+0.4330127018922194*eyr[7]*c+0.4330127018922194*eyl[7]*c-0.25*eyr[5]*c+0.25*eyl[5]*c; 
  incr[6] = (-1.479019945774904*bzr[11]*c2)+1.479019945774904*bzl[11]*c2+1.25*bzr[6]*c2+1.25*bzl[6]*c2-0.9682458365518543*bzr[3]*c2+0.9682458365518543*bzl[3]*c2+0.5590169943749476*bzr[2]*c2+0.5590169943749476*bzl[2]*c2+1.479019945774904*eyr[11]*c+1.479019945774904*eyl[11]*c-1.25*eyr[6]*c+1.25*eyl[6]*c+0.9682458365518543*eyr[3]*c+0.9682458365518543*eyl[3]*c-0.5590169943749476*eyr[2]*c+0.5590169943749476*eyl[2]*c; 
  incr[7] = (-0.9682458365518543*bzr[10]*c2)-0.9682458365518543*bzl[10]*c2+0.75*bzr[7]*c2-0.75*bzl[7]*c2-0.4330127018922194*bzr[5]*c2-0.4330127018922194*bzl[5]*c2+0.9682458365518543*eyr[10]*c-0.9682458365518543*eyl[10]*c-0.75*eyr[7]*c-0.75*eyl[7]*c+0.4330127018922194*eyr[5]*c-0.4330127018922194*eyl[5]*c; 
  incr[8] = (-1.984313483298443*bzr[13]*c2)-1.984313483298443*bzl[13]*c2+1.75*bzr[8]*c2-1.75*bzl[8]*c2-1.479019945774904*bzr[4]*c2-1.479019945774904*bzl[4]*c2+1.14564392373896*bzr[1]*c2-1.14564392373896*bzl[1]*c2-0.6614378277661477*bzr[0]*c2-0.6614378277661477*bzl[0]*c2+1.984313483298443*eyr[13]*c-1.984313483298443*eyl[13]*c-1.75*eyr[8]*c-1.75*eyl[8]*c+1.479019945774904*eyr[4]*c-1.479019945774904*eyl[4]*c-1.14564392373896*eyr[1]*c-1.14564392373896*eyl[1]*c+0.6614378277661477*eyr[0]*c-0.6614378277661477*eyl[0]*c; 
  incr[9] = (-0.4330127018922193*bzr[12]*c2)+0.4330127018922193*bzl[12]*c2+0.25*bzr[9]*c2+0.25*bzl[9]*c2+0.4330127018922193*eyr[12]*c+0.4330127018922193*eyl[12]*c-0.25*eyr[9]*c+0.25*eyl[9]*c; 
  incr[10] = 1.25*bzr[10]*c2+1.25*bzl[10]*c2-0.9682458365518543*bzr[7]*c2+0.9682458365518543*bzl[7]*c2+0.5590169943749475*bzr[5]*c2+0.5590169943749475*bzl[5]*c2-1.25*eyr[10]*c+1.25*eyl[10]*c+0.9682458365518543*eyr[7]*c+0.9682458365518543*eyl[7]*c-0.5590169943749475*eyr[5]*c+0.5590169943749475*eyl[5]*c; 
  incr[11] = 1.75*bzr[11]*c2-1.75*bzl[11]*c2-1.479019945774904*bzr[6]*c2-1.479019945774904*bzl[6]*c2+1.14564392373896*bzr[3]*c2-1.14564392373896*bzl[3]*c2-0.6614378277661477*bzr[2]*c2-0.6614378277661477*bzl[2]*c2-1.75*eyr[11]*c-1.75*eyl[11]*c+1.479019945774904*eyr[6]*c-1.479019945774904*eyl[6]*c-1.14564392373896*eyr[3]*c-1.14564392373896*eyl[3]*c+0.6614378277661477*eyr[2]*c-0.6614378277661477*eyl[2]*c; 
  incr[12] = 0.75*bzr[12]*c2-0.75*bzl[12]*c2-0.4330127018922193*bzr[9]*c2-0.4330127018922193*bzl[9]*c2-0.75*eyr[12]*c-0.75*eyl[12]*c+0.4330127018922193*eyr[9]*c-0.4330127018922193*eyl[9]*c; 
  incr[13] = 2.25*bzr[13]*c2+2.25*bzl[13]*c2-1.984313483298443*bzr[8]*c2+1.984313483298443*bzl[8]*c2+1.677050983124842*bzr[4]*c2+1.677050983124842*bzl[4]*c2-1.299038105676658*bzr[1]*c2+1.299038105676658*bzl[1]*c2+0.75*bzr[0]*c2+0.75*bzl[0]*c2-2.25*eyr[13]*c+2.25*eyl[13]*c+1.984313483298443*eyr[8]*c+1.984313483298443*eyl[8]*c-1.677050983124842*eyr[4]*c+1.677050983124842*eyl[4]*c+1.299038105676658*eyr[1]*c+1.299038105676658*eyl[1]*c-0.75*eyr[0]*c+0.75*eyl[0]*c; 
  incr[14] = 0.25*bzr[14]*c2+0.25*bzl[14]*c2-0.25*eyr[14]*c+0.25*eyl[14]*c; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 
  outEyr[3] += incr[3]*dx1; 
  outEyr[4] += incr[4]*dx1; 
  outEyr[5] += incr[5]*dx1; 
  outEyr[6] += incr[6]*dx1; 
  outEyr[7] += incr[7]*dx1; 
  outEyr[8] += incr[8]*dx1; 
  outEyr[9] += incr[9]*dx1; 
  outEyr[10] += incr[10]*dx1; 
  outEyr[11] += incr[11]*dx1; 
  outEyr[12] += incr[12]*dx1; 
  outEyr[13] += incr[13]*dx1; 
  outEyr[14] += incr[14]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += incr[1]*dx1; 
  outEyl[2] += -1.0*incr[2]*dx1; 
  outEyl[3] += incr[3]*dx1; 
  outEyl[4] += -1.0*incr[4]*dx1; 
  outEyl[5] += -1.0*incr[5]*dx1; 
  outEyl[6] += -1.0*incr[6]*dx1; 
  outEyl[7] += incr[7]*dx1; 
  outEyl[8] += incr[8]*dx1; 
  outEyl[9] += -1.0*incr[9]*dx1; 
  outEyl[10] += -1.0*incr[10]*dx1; 
  outEyl[11] += incr[11]*dx1; 
  outEyl[12] += incr[12]*dx1; 
  outEyl[13] += -1.0*incr[13]*dx1; 
  outEyl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = (-0.75*byr[13]*c2)-0.75*byl[13]*c2+0.6614378277661477*byr[8]*c2-0.6614378277661477*byl[8]*c2-0.5590169943749475*byr[4]*c2-0.5590169943749475*byl[4]*c2+0.4330127018922193*byr[1]*c2-0.4330127018922193*byl[1]*c2-0.25*byr[0]*c2-0.25*byl[0]*c2-0.75*ezr[13]*c+0.75*ezl[13]*c+0.6614378277661477*ezr[8]*c+0.6614378277661477*ezl[8]*c-0.5590169943749475*ezr[4]*c+0.5590169943749475*ezl[4]*c+0.4330127018922193*ezr[1]*c+0.4330127018922193*ezl[1]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = 1.299038105676658*byr[13]*c2+1.299038105676658*byl[13]*c2-1.14564392373896*byr[8]*c2+1.14564392373896*byl[8]*c2+0.9682458365518543*byr[4]*c2+0.9682458365518543*byl[4]*c2-0.75*byr[1]*c2+0.75*byl[1]*c2+0.4330127018922193*byr[0]*c2+0.4330127018922193*byl[0]*c2+1.299038105676658*ezr[13]*c-1.299038105676658*ezl[13]*c-1.14564392373896*ezr[8]*c-1.14564392373896*ezl[8]*c+0.9682458365518543*ezr[4]*c-0.9682458365518543*ezl[4]*c-0.75*ezr[1]*c-0.75*ezl[1]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 
  incr[2] = 0.6614378277661477*byr[11]*c2-0.6614378277661477*byl[11]*c2-0.5590169943749476*byr[6]*c2-0.5590169943749476*byl[6]*c2+0.4330127018922193*byr[3]*c2-0.4330127018922193*byl[3]*c2-0.25*byr[2]*c2-0.25*byl[2]*c2+0.6614378277661477*ezr[11]*c+0.6614378277661477*ezl[11]*c-0.5590169943749476*ezr[6]*c+0.5590169943749476*ezl[6]*c+0.4330127018922193*ezr[3]*c+0.4330127018922193*ezl[3]*c-0.25*ezr[2]*c+0.25*ezl[2]*c; 
  incr[3] = (-1.14564392373896*byr[11]*c2)+1.14564392373896*byl[11]*c2+0.9682458365518543*byr[6]*c2+0.9682458365518543*byl[6]*c2-0.75*byr[3]*c2+0.75*byl[3]*c2+0.4330127018922193*byr[2]*c2+0.4330127018922193*byl[2]*c2-1.14564392373896*ezr[11]*c-1.14564392373896*ezl[11]*c+0.9682458365518543*ezr[6]*c-0.9682458365518543*ezl[6]*c-0.75*ezr[3]*c-0.75*ezl[3]*c+0.4330127018922193*ezr[2]*c-0.4330127018922193*ezl[2]*c; 
  incr[4] = (-1.677050983124842*byr[13]*c2)-1.677050983124842*byl[13]*c2+1.479019945774904*byr[8]*c2-1.479019945774904*byl[8]*c2-1.25*byr[4]*c2-1.25*byl[4]*c2+0.9682458365518543*byr[1]*c2-0.9682458365518543*byl[1]*c2-0.5590169943749475*byr[0]*c2-0.5590169943749475*byl[0]*c2-1.677050983124842*ezr[13]*c+1.677050983124842*ezl[13]*c+1.479019945774904*ezr[8]*c+1.479019945774904*ezl[8]*c-1.25*ezr[4]*c+1.25*ezl[4]*c+0.9682458365518543*ezr[1]*c+0.9682458365518543*ezl[1]*c-0.5590169943749475*ezr[0]*c+0.5590169943749475*ezl[0]*c; 
  incr[5] = (-0.5590169943749475*byr[10]*c2)-0.5590169943749475*byl[10]*c2+0.4330127018922194*byr[7]*c2-0.4330127018922194*byl[7]*c2-0.25*byr[5]*c2-0.25*byl[5]*c2-0.5590169943749475*ezr[10]*c+0.5590169943749475*ezl[10]*c+0.4330127018922194*ezr[7]*c+0.4330127018922194*ezl[7]*c-0.25*ezr[5]*c+0.25*ezl[5]*c; 
  incr[6] = 1.479019945774904*byr[11]*c2-1.479019945774904*byl[11]*c2-1.25*byr[6]*c2-1.25*byl[6]*c2+0.9682458365518543*byr[3]*c2-0.9682458365518543*byl[3]*c2-0.5590169943749476*byr[2]*c2-0.5590169943749476*byl[2]*c2+1.479019945774904*ezr[11]*c+1.479019945774904*ezl[11]*c-1.25*ezr[6]*c+1.25*ezl[6]*c+0.9682458365518543*ezr[3]*c+0.9682458365518543*ezl[3]*c-0.5590169943749476*ezr[2]*c+0.5590169943749476*ezl[2]*c; 
  incr[7] = 0.9682458365518543*byr[10]*c2+0.9682458365518543*byl[10]*c2-0.75*byr[7]*c2+0.75*byl[7]*c2+0.4330127018922194*byr[5]*c2+0.4330127018922194*byl[5]*c2+0.9682458365518543*ezr[10]*c-0.9682458365518543*ezl[10]*c-0.75*ezr[7]*c-0.75*ezl[7]*c+0.4330127018922194*ezr[5]*c-0.4330127018922194*ezl[5]*c; 
  incr[8] = 1.984313483298443*byr[13]*c2+1.984313483298443*byl[13]*c2-1.75*byr[8]*c2+1.75*byl[8]*c2+1.479019945774904*byr[4]*c2+1.479019945774904*byl[4]*c2-1.14564392373896*byr[1]*c2+1.14564392373896*byl[1]*c2+0.6614378277661477*byr[0]*c2+0.6614378277661477*byl[0]*c2+1.984313483298443*ezr[13]*c-1.984313483298443*ezl[13]*c-1.75*ezr[8]*c-1.75*ezl[8]*c+1.479019945774904*ezr[4]*c-1.479019945774904*ezl[4]*c-1.14564392373896*ezr[1]*c-1.14564392373896*ezl[1]*c+0.6614378277661477*ezr[0]*c-0.6614378277661477*ezl[0]*c; 
  incr[9] = 0.4330127018922193*byr[12]*c2-0.4330127018922193*byl[12]*c2-0.25*byr[9]*c2-0.25*byl[9]*c2+0.4330127018922193*ezr[12]*c+0.4330127018922193*ezl[12]*c-0.25*ezr[9]*c+0.25*ezl[9]*c; 
  incr[10] = (-1.25*byr[10]*c2)-1.25*byl[10]*c2+0.9682458365518543*byr[7]*c2-0.9682458365518543*byl[7]*c2-0.5590169943749475*byr[5]*c2-0.5590169943749475*byl[5]*c2-1.25*ezr[10]*c+1.25*ezl[10]*c+0.9682458365518543*ezr[7]*c+0.9682458365518543*ezl[7]*c-0.5590169943749475*ezr[5]*c+0.5590169943749475*ezl[5]*c; 
  incr[11] = (-1.75*byr[11]*c2)+1.75*byl[11]*c2+1.479019945774904*byr[6]*c2+1.479019945774904*byl[6]*c2-1.14564392373896*byr[3]*c2+1.14564392373896*byl[3]*c2+0.6614378277661477*byr[2]*c2+0.6614378277661477*byl[2]*c2-1.75*ezr[11]*c-1.75*ezl[11]*c+1.479019945774904*ezr[6]*c-1.479019945774904*ezl[6]*c-1.14564392373896*ezr[3]*c-1.14564392373896*ezl[3]*c+0.6614378277661477*ezr[2]*c-0.6614378277661477*ezl[2]*c; 
  incr[12] = (-0.75*byr[12]*c2)+0.75*byl[12]*c2+0.4330127018922193*byr[9]*c2+0.4330127018922193*byl[9]*c2-0.75*ezr[12]*c-0.75*ezl[12]*c+0.4330127018922193*ezr[9]*c-0.4330127018922193*ezl[9]*c; 
  incr[13] = (-2.25*byr[13]*c2)-2.25*byl[13]*c2+1.984313483298443*byr[8]*c2-1.984313483298443*byl[8]*c2-1.677050983124842*byr[4]*c2-1.677050983124842*byl[4]*c2+1.299038105676658*byr[1]*c2-1.299038105676658*byl[1]*c2-0.75*byr[0]*c2-0.75*byl[0]*c2-2.25*ezr[13]*c+2.25*ezl[13]*c+1.984313483298443*ezr[8]*c+1.984313483298443*ezl[8]*c-1.677050983124842*ezr[4]*c+1.677050983124842*ezl[4]*c+1.299038105676658*ezr[1]*c+1.299038105676658*ezl[1]*c-0.75*ezr[0]*c+0.75*ezl[0]*c; 
  incr[14] = (-0.25*byr[14]*c2)-0.25*byl[14]*c2-0.25*ezr[14]*c+0.25*ezl[14]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 
  outEzr[3] += incr[3]*dx1; 
  outEzr[4] += incr[4]*dx1; 
  outEzr[5] += incr[5]*dx1; 
  outEzr[6] += incr[6]*dx1; 
  outEzr[7] += incr[7]*dx1; 
  outEzr[8] += incr[8]*dx1; 
  outEzr[9] += incr[9]*dx1; 
  outEzr[10] += incr[10]*dx1; 
  outEzr[11] += incr[11]*dx1; 
  outEzr[12] += incr[12]*dx1; 
  outEzr[13] += incr[13]*dx1; 
  outEzr[14] += incr[14]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += incr[1]*dx1; 
  outEzl[2] += -1.0*incr[2]*dx1; 
  outEzl[3] += incr[3]*dx1; 
  outEzl[4] += -1.0*incr[4]*dx1; 
  outEzl[5] += -1.0*incr[5]*dx1; 
  outEzl[6] += -1.0*incr[6]*dx1; 
  outEzl[7] += incr[7]*dx1; 
  outEzl[8] += incr[8]*dx1; 
  outEzl[9] += -1.0*incr[9]*dx1; 
  outEzl[10] += -1.0*incr[10]*dx1; 
  outEzl[11] += incr[11]*dx1; 
  outEzl[12] += incr[12]*dx1; 
  outEzl[13] += -1.0*incr[13]*dx1; 
  outEzl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = (-0.75*bxr[13]*c*gamma)+0.75*bxl[13]*c*gamma+0.6614378277661477*bxr[8]*c*gamma+0.6614378277661477*bxl[8]*c*gamma-0.5590169943749475*bxr[4]*c*gamma+0.5590169943749475*bxl[4]*c*gamma+0.4330127018922193*bxr[1]*c*gamma+0.4330127018922193*bxl[1]*c*gamma-0.25*bxr[0]*c*gamma+0.25*bxl[0]*c*gamma+0.75*psr[13]*gamma+0.75*psl[13]*gamma-0.6614378277661477*psr[8]*gamma+0.6614378277661477*psl[8]*gamma+0.5590169943749475*psr[4]*gamma+0.5590169943749475*psl[4]*gamma-0.4330127018922193*psr[1]*gamma+0.4330127018922193*psl[1]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = 1.299038105676658*bxr[13]*c*gamma-1.299038105676658*bxl[13]*c*gamma-1.14564392373896*bxr[8]*c*gamma-1.14564392373896*bxl[8]*c*gamma+0.9682458365518543*bxr[4]*c*gamma-0.9682458365518543*bxl[4]*c*gamma-0.75*bxr[1]*c*gamma-0.75*bxl[1]*c*gamma+0.4330127018922193*bxr[0]*c*gamma-0.4330127018922193*bxl[0]*c*gamma-1.299038105676658*psr[13]*gamma-1.299038105676658*psl[13]*gamma+1.14564392373896*psr[8]*gamma-1.14564392373896*psl[8]*gamma-0.9682458365518543*psr[4]*gamma-0.9682458365518543*psl[4]*gamma+0.75*psr[1]*gamma-0.75*psl[1]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 
  incr[2] = 0.6614378277661477*bxr[11]*c*gamma+0.6614378277661477*bxl[11]*c*gamma-0.5590169943749476*bxr[6]*c*gamma+0.5590169943749476*bxl[6]*c*gamma+0.4330127018922193*bxr[3]*c*gamma+0.4330127018922193*bxl[3]*c*gamma-0.25*bxr[2]*c*gamma+0.25*bxl[2]*c*gamma-0.6614378277661477*psr[11]*gamma+0.6614378277661477*psl[11]*gamma+0.5590169943749476*psr[6]*gamma+0.5590169943749476*psl[6]*gamma-0.4330127018922193*psr[3]*gamma+0.4330127018922193*psl[3]*gamma+0.25*psr[2]*gamma+0.25*psl[2]*gamma; 
  incr[3] = (-1.14564392373896*bxr[11]*c*gamma)-1.14564392373896*bxl[11]*c*gamma+0.9682458365518543*bxr[6]*c*gamma-0.9682458365518543*bxl[6]*c*gamma-0.75*bxr[3]*c*gamma-0.75*bxl[3]*c*gamma+0.4330127018922193*bxr[2]*c*gamma-0.4330127018922193*bxl[2]*c*gamma+1.14564392373896*psr[11]*gamma-1.14564392373896*psl[11]*gamma-0.9682458365518543*psr[6]*gamma-0.9682458365518543*psl[6]*gamma+0.75*psr[3]*gamma-0.75*psl[3]*gamma-0.4330127018922193*psr[2]*gamma-0.4330127018922193*psl[2]*gamma; 
  incr[4] = (-1.677050983124842*bxr[13]*c*gamma)+1.677050983124842*bxl[13]*c*gamma+1.479019945774904*bxr[8]*c*gamma+1.479019945774904*bxl[8]*c*gamma-1.25*bxr[4]*c*gamma+1.25*bxl[4]*c*gamma+0.9682458365518543*bxr[1]*c*gamma+0.9682458365518543*bxl[1]*c*gamma-0.5590169943749475*bxr[0]*c*gamma+0.5590169943749475*bxl[0]*c*gamma+1.677050983124842*psr[13]*gamma+1.677050983124842*psl[13]*gamma-1.479019945774904*psr[8]*gamma+1.479019945774904*psl[8]*gamma+1.25*psr[4]*gamma+1.25*psl[4]*gamma-0.9682458365518543*psr[1]*gamma+0.9682458365518543*psl[1]*gamma+0.5590169943749475*psr[0]*gamma+0.5590169943749475*psl[0]*gamma; 
  incr[5] = (-0.5590169943749475*bxr[10]*c*gamma)+0.5590169943749475*bxl[10]*c*gamma+0.4330127018922194*bxr[7]*c*gamma+0.4330127018922194*bxl[7]*c*gamma-0.25*bxr[5]*c*gamma+0.25*bxl[5]*c*gamma+0.5590169943749475*psr[10]*gamma+0.5590169943749475*psl[10]*gamma-0.4330127018922194*psr[7]*gamma+0.4330127018922194*psl[7]*gamma+0.25*psr[5]*gamma+0.25*psl[5]*gamma; 
  incr[6] = 1.479019945774904*bxr[11]*c*gamma+1.479019945774904*bxl[11]*c*gamma-1.25*bxr[6]*c*gamma+1.25*bxl[6]*c*gamma+0.9682458365518543*bxr[3]*c*gamma+0.9682458365518543*bxl[3]*c*gamma-0.5590169943749476*bxr[2]*c*gamma+0.5590169943749476*bxl[2]*c*gamma-1.479019945774904*psr[11]*gamma+1.479019945774904*psl[11]*gamma+1.25*psr[6]*gamma+1.25*psl[6]*gamma-0.9682458365518543*psr[3]*gamma+0.9682458365518543*psl[3]*gamma+0.5590169943749476*psr[2]*gamma+0.5590169943749476*psl[2]*gamma; 
  incr[7] = 0.9682458365518543*bxr[10]*c*gamma-0.9682458365518543*bxl[10]*c*gamma-0.75*bxr[7]*c*gamma-0.75*bxl[7]*c*gamma+0.4330127018922194*bxr[5]*c*gamma-0.4330127018922194*bxl[5]*c*gamma-0.9682458365518543*psr[10]*gamma-0.9682458365518543*psl[10]*gamma+0.75*psr[7]*gamma-0.75*psl[7]*gamma-0.4330127018922194*psr[5]*gamma-0.4330127018922194*psl[5]*gamma; 
  incr[8] = 1.984313483298443*bxr[13]*c*gamma-1.984313483298443*bxl[13]*c*gamma-1.75*bxr[8]*c*gamma-1.75*bxl[8]*c*gamma+1.479019945774904*bxr[4]*c*gamma-1.479019945774904*bxl[4]*c*gamma-1.14564392373896*bxr[1]*c*gamma-1.14564392373896*bxl[1]*c*gamma+0.6614378277661477*bxr[0]*c*gamma-0.6614378277661477*bxl[0]*c*gamma-1.984313483298443*psr[13]*gamma-1.984313483298443*psl[13]*gamma+1.75*psr[8]*gamma-1.75*psl[8]*gamma-1.479019945774904*psr[4]*gamma-1.479019945774904*psl[4]*gamma+1.14564392373896*psr[1]*gamma-1.14564392373896*psl[1]*gamma-0.6614378277661477*psr[0]*gamma-0.6614378277661477*psl[0]*gamma; 
  incr[9] = 0.4330127018922193*bxr[12]*c*gamma+0.4330127018922193*bxl[12]*c*gamma-0.25*bxr[9]*c*gamma+0.25*bxl[9]*c*gamma-0.4330127018922193*psr[12]*gamma+0.4330127018922193*psl[12]*gamma+0.25*psr[9]*gamma+0.25*psl[9]*gamma; 
  incr[10] = (-1.25*bxr[10]*c*gamma)+1.25*bxl[10]*c*gamma+0.9682458365518543*bxr[7]*c*gamma+0.9682458365518543*bxl[7]*c*gamma-0.5590169943749475*bxr[5]*c*gamma+0.5590169943749475*bxl[5]*c*gamma+1.25*psr[10]*gamma+1.25*psl[10]*gamma-0.9682458365518543*psr[7]*gamma+0.9682458365518543*psl[7]*gamma+0.5590169943749475*psr[5]*gamma+0.5590169943749475*psl[5]*gamma; 
  incr[11] = (-1.75*bxr[11]*c*gamma)-1.75*bxl[11]*c*gamma+1.479019945774904*bxr[6]*c*gamma-1.479019945774904*bxl[6]*c*gamma-1.14564392373896*bxr[3]*c*gamma-1.14564392373896*bxl[3]*c*gamma+0.6614378277661477*bxr[2]*c*gamma-0.6614378277661477*bxl[2]*c*gamma+1.75*psr[11]*gamma-1.75*psl[11]*gamma-1.479019945774904*psr[6]*gamma-1.479019945774904*psl[6]*gamma+1.14564392373896*psr[3]*gamma-1.14564392373896*psl[3]*gamma-0.6614378277661477*psr[2]*gamma-0.6614378277661477*psl[2]*gamma; 
  incr[12] = (-0.75*bxr[12]*c*gamma)-0.75*bxl[12]*c*gamma+0.4330127018922193*bxr[9]*c*gamma-0.4330127018922193*bxl[9]*c*gamma+0.75*psr[12]*gamma-0.75*psl[12]*gamma-0.4330127018922193*psr[9]*gamma-0.4330127018922193*psl[9]*gamma; 
  incr[13] = (-2.25*bxr[13]*c*gamma)+2.25*bxl[13]*c*gamma+1.984313483298443*bxr[8]*c*gamma+1.984313483298443*bxl[8]*c*gamma-1.677050983124842*bxr[4]*c*gamma+1.677050983124842*bxl[4]*c*gamma+1.299038105676658*bxr[1]*c*gamma+1.299038105676658*bxl[1]*c*gamma-0.75*bxr[0]*c*gamma+0.75*bxl[0]*c*gamma+2.25*psr[13]*gamma+2.25*psl[13]*gamma-1.984313483298443*psr[8]*gamma+1.984313483298443*psl[8]*gamma+1.677050983124842*psr[4]*gamma+1.677050983124842*psl[4]*gamma-1.299038105676658*psr[1]*gamma+1.299038105676658*psl[1]*gamma+0.75*psr[0]*gamma+0.75*psl[0]*gamma; 
  incr[14] = (-0.25*bxr[14]*c*gamma)+0.25*bxl[14]*c*gamma+0.25*psr[14]*gamma+0.25*psl[14]*gamma; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 
  outBxr[3] += incr[3]*dx1; 
  outBxr[4] += incr[4]*dx1; 
  outBxr[5] += incr[5]*dx1; 
  outBxr[6] += incr[6]*dx1; 
  outBxr[7] += incr[7]*dx1; 
  outBxr[8] += incr[8]*dx1; 
  outBxr[9] += incr[9]*dx1; 
  outBxr[10] += incr[10]*dx1; 
  outBxr[11] += incr[11]*dx1; 
  outBxr[12] += incr[12]*dx1; 
  outBxr[13] += incr[13]*dx1; 
  outBxr[14] += incr[14]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += incr[1]*dx1; 
  outBxl[2] += -1.0*incr[2]*dx1; 
  outBxl[3] += incr[3]*dx1; 
  outBxl[4] += -1.0*incr[4]*dx1; 
  outBxl[5] += -1.0*incr[5]*dx1; 
  outBxl[6] += -1.0*incr[6]*dx1; 
  outBxl[7] += incr[7]*dx1; 
  outBxl[8] += incr[8]*dx1; 
  outBxl[9] += -1.0*incr[9]*dx1; 
  outBxl[10] += -1.0*incr[10]*dx1; 
  outBxl[11] += incr[11]*dx1; 
  outBxl[12] += incr[12]*dx1; 
  outBxl[13] += -1.0*incr[13]*dx1; 
  outBxl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = (-0.75*byr[13]*c)+0.75*byl[13]*c+0.6614378277661477*byr[8]*c+0.6614378277661477*byl[8]*c-0.5590169943749475*byr[4]*c+0.5590169943749475*byl[4]*c+0.4330127018922193*byr[1]*c+0.4330127018922193*byl[1]*c-0.25*byr[0]*c+0.25*byl[0]*c-0.75*ezr[13]-0.75*ezl[13]+0.6614378277661477*ezr[8]-0.6614378277661477*ezl[8]-0.5590169943749475*ezr[4]-0.5590169943749475*ezl[4]+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*ezr[0]-0.25*ezl[0]; 
  incr[1] = 1.299038105676658*byr[13]*c-1.299038105676658*byl[13]*c-1.14564392373896*byr[8]*c-1.14564392373896*byl[8]*c+0.9682458365518543*byr[4]*c-0.9682458365518543*byl[4]*c-0.75*byr[1]*c-0.75*byl[1]*c+0.4330127018922193*byr[0]*c-0.4330127018922193*byl[0]*c+1.299038105676658*ezr[13]+1.299038105676658*ezl[13]-1.14564392373896*ezr[8]+1.14564392373896*ezl[8]+0.9682458365518543*ezr[4]+0.9682458365518543*ezl[4]-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*ezr[0]+0.4330127018922193*ezl[0]; 
  incr[2] = 0.6614378277661477*byr[11]*c+0.6614378277661477*byl[11]*c-0.5590169943749476*byr[6]*c+0.5590169943749476*byl[6]*c+0.4330127018922193*byr[3]*c+0.4330127018922193*byl[3]*c-0.25*byr[2]*c+0.25*byl[2]*c+0.6614378277661477*ezr[11]-0.6614378277661477*ezl[11]-0.5590169943749476*ezr[6]-0.5590169943749476*ezl[6]+0.4330127018922193*ezr[3]-0.4330127018922193*ezl[3]-0.25*ezr[2]-0.25*ezl[2]; 
  incr[3] = (-1.14564392373896*byr[11]*c)-1.14564392373896*byl[11]*c+0.9682458365518543*byr[6]*c-0.9682458365518543*byl[6]*c-0.75*byr[3]*c-0.75*byl[3]*c+0.4330127018922193*byr[2]*c-0.4330127018922193*byl[2]*c-1.14564392373896*ezr[11]+1.14564392373896*ezl[11]+0.9682458365518543*ezr[6]+0.9682458365518543*ezl[6]-0.75*ezr[3]+0.75*ezl[3]+0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]; 
  incr[4] = (-1.677050983124842*byr[13]*c)+1.677050983124842*byl[13]*c+1.479019945774904*byr[8]*c+1.479019945774904*byl[8]*c-1.25*byr[4]*c+1.25*byl[4]*c+0.9682458365518543*byr[1]*c+0.9682458365518543*byl[1]*c-0.5590169943749475*byr[0]*c+0.5590169943749475*byl[0]*c-1.677050983124842*ezr[13]-1.677050983124842*ezl[13]+1.479019945774904*ezr[8]-1.479019945774904*ezl[8]-1.25*ezr[4]-1.25*ezl[4]+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*ezr[0]-0.5590169943749475*ezl[0]; 
  incr[5] = (-0.5590169943749475*byr[10]*c)+0.5590169943749475*byl[10]*c+0.4330127018922194*byr[7]*c+0.4330127018922194*byl[7]*c-0.25*byr[5]*c+0.25*byl[5]*c-0.5590169943749475*ezr[10]-0.5590169943749475*ezl[10]+0.4330127018922194*ezr[7]-0.4330127018922194*ezl[7]-0.25*ezr[5]-0.25*ezl[5]; 
  incr[6] = 1.479019945774904*byr[11]*c+1.479019945774904*byl[11]*c-1.25*byr[6]*c+1.25*byl[6]*c+0.9682458365518543*byr[3]*c+0.9682458365518543*byl[3]*c-0.5590169943749476*byr[2]*c+0.5590169943749476*byl[2]*c+1.479019945774904*ezr[11]-1.479019945774904*ezl[11]-1.25*ezr[6]-1.25*ezl[6]+0.9682458365518543*ezr[3]-0.9682458365518543*ezl[3]-0.5590169943749476*ezr[2]-0.5590169943749476*ezl[2]; 
  incr[7] = 0.9682458365518543*byr[10]*c-0.9682458365518543*byl[10]*c-0.75*byr[7]*c-0.75*byl[7]*c+0.4330127018922194*byr[5]*c-0.4330127018922194*byl[5]*c+0.9682458365518543*ezr[10]+0.9682458365518543*ezl[10]-0.75*ezr[7]+0.75*ezl[7]+0.4330127018922194*ezr[5]+0.4330127018922194*ezl[5]; 
  incr[8] = 1.984313483298443*byr[13]*c-1.984313483298443*byl[13]*c-1.75*byr[8]*c-1.75*byl[8]*c+1.479019945774904*byr[4]*c-1.479019945774904*byl[4]*c-1.14564392373896*byr[1]*c-1.14564392373896*byl[1]*c+0.6614378277661477*byr[0]*c-0.6614378277661477*byl[0]*c+1.984313483298443*ezr[13]+1.984313483298443*ezl[13]-1.75*ezr[8]+1.75*ezl[8]+1.479019945774904*ezr[4]+1.479019945774904*ezl[4]-1.14564392373896*ezr[1]+1.14564392373896*ezl[1]+0.6614378277661477*ezr[0]+0.6614378277661477*ezl[0]; 
  incr[9] = 0.4330127018922193*byr[12]*c+0.4330127018922193*byl[12]*c-0.25*byr[9]*c+0.25*byl[9]*c+0.4330127018922193*ezr[12]-0.4330127018922193*ezl[12]-0.25*ezr[9]-0.25*ezl[9]; 
  incr[10] = (-1.25*byr[10]*c)+1.25*byl[10]*c+0.9682458365518543*byr[7]*c+0.9682458365518543*byl[7]*c-0.5590169943749475*byr[5]*c+0.5590169943749475*byl[5]*c-1.25*ezr[10]-1.25*ezl[10]+0.9682458365518543*ezr[7]-0.9682458365518543*ezl[7]-0.5590169943749475*ezr[5]-0.5590169943749475*ezl[5]; 
  incr[11] = (-1.75*byr[11]*c)-1.75*byl[11]*c+1.479019945774904*byr[6]*c-1.479019945774904*byl[6]*c-1.14564392373896*byr[3]*c-1.14564392373896*byl[3]*c+0.6614378277661477*byr[2]*c-0.6614378277661477*byl[2]*c-1.75*ezr[11]+1.75*ezl[11]+1.479019945774904*ezr[6]+1.479019945774904*ezl[6]-1.14564392373896*ezr[3]+1.14564392373896*ezl[3]+0.6614378277661477*ezr[2]+0.6614378277661477*ezl[2]; 
  incr[12] = (-0.75*byr[12]*c)-0.75*byl[12]*c+0.4330127018922193*byr[9]*c-0.4330127018922193*byl[9]*c-0.75*ezr[12]+0.75*ezl[12]+0.4330127018922193*ezr[9]+0.4330127018922193*ezl[9]; 
  incr[13] = (-2.25*byr[13]*c)+2.25*byl[13]*c+1.984313483298443*byr[8]*c+1.984313483298443*byl[8]*c-1.677050983124842*byr[4]*c+1.677050983124842*byl[4]*c+1.299038105676658*byr[1]*c+1.299038105676658*byl[1]*c-0.75*byr[0]*c+0.75*byl[0]*c-2.25*ezr[13]-2.25*ezl[13]+1.984313483298443*ezr[8]-1.984313483298443*ezl[8]-1.677050983124842*ezr[4]-1.677050983124842*ezl[4]+1.299038105676658*ezr[1]-1.299038105676658*ezl[1]-0.75*ezr[0]-0.75*ezl[0]; 
  incr[14] = (-0.25*byr[14]*c)+0.25*byl[14]*c-0.25*ezr[14]-0.25*ezl[14]; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 
  outByr[3] += incr[3]*dx1; 
  outByr[4] += incr[4]*dx1; 
  outByr[5] += incr[5]*dx1; 
  outByr[6] += incr[6]*dx1; 
  outByr[7] += incr[7]*dx1; 
  outByr[8] += incr[8]*dx1; 
  outByr[9] += incr[9]*dx1; 
  outByr[10] += incr[10]*dx1; 
  outByr[11] += incr[11]*dx1; 
  outByr[12] += incr[12]*dx1; 
  outByr[13] += incr[13]*dx1; 
  outByr[14] += incr[14]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += incr[1]*dx1; 
  outByl[2] += -1.0*incr[2]*dx1; 
  outByl[3] += incr[3]*dx1; 
  outByl[4] += -1.0*incr[4]*dx1; 
  outByl[5] += -1.0*incr[5]*dx1; 
  outByl[6] += -1.0*incr[6]*dx1; 
  outByl[7] += incr[7]*dx1; 
  outByl[8] += incr[8]*dx1; 
  outByl[9] += -1.0*incr[9]*dx1; 
  outByl[10] += -1.0*incr[10]*dx1; 
  outByl[11] += incr[11]*dx1; 
  outByl[12] += incr[12]*dx1; 
  outByl[13] += -1.0*incr[13]*dx1; 
  outByl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = (-0.75*bzr[13]*c)+0.75*bzl[13]*c+0.6614378277661477*bzr[8]*c+0.6614378277661477*bzl[8]*c-0.5590169943749475*bzr[4]*c+0.5590169943749475*bzl[4]*c+0.4330127018922193*bzr[1]*c+0.4330127018922193*bzl[1]*c-0.25*bzr[0]*c+0.25*bzl[0]*c+0.75*eyr[13]+0.75*eyl[13]-0.6614378277661477*eyr[8]+0.6614378277661477*eyl[8]+0.5590169943749475*eyr[4]+0.5590169943749475*eyl[4]-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*eyr[0]+0.25*eyl[0]; 
  incr[1] = 1.299038105676658*bzr[13]*c-1.299038105676658*bzl[13]*c-1.14564392373896*bzr[8]*c-1.14564392373896*bzl[8]*c+0.9682458365518543*bzr[4]*c-0.9682458365518543*bzl[4]*c-0.75*bzr[1]*c-0.75*bzl[1]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c-1.299038105676658*eyr[13]-1.299038105676658*eyl[13]+1.14564392373896*eyr[8]-1.14564392373896*eyl[8]-0.9682458365518543*eyr[4]-0.9682458365518543*eyl[4]+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0]; 
  incr[2] = 0.6614378277661477*bzr[11]*c+0.6614378277661477*bzl[11]*c-0.5590169943749476*bzr[6]*c+0.5590169943749476*bzl[6]*c+0.4330127018922193*bzr[3]*c+0.4330127018922193*bzl[3]*c-0.25*bzr[2]*c+0.25*bzl[2]*c-0.6614378277661477*eyr[11]+0.6614378277661477*eyl[11]+0.5590169943749476*eyr[6]+0.5590169943749476*eyl[6]-0.4330127018922193*eyr[3]+0.4330127018922193*eyl[3]+0.25*eyr[2]+0.25*eyl[2]; 
  incr[3] = (-1.14564392373896*bzr[11]*c)-1.14564392373896*bzl[11]*c+0.9682458365518543*bzr[6]*c-0.9682458365518543*bzl[6]*c-0.75*bzr[3]*c-0.75*bzl[3]*c+0.4330127018922193*bzr[2]*c-0.4330127018922193*bzl[2]*c+1.14564392373896*eyr[11]-1.14564392373896*eyl[11]-0.9682458365518543*eyr[6]-0.9682458365518543*eyl[6]+0.75*eyr[3]-0.75*eyl[3]-0.4330127018922193*eyr[2]-0.4330127018922193*eyl[2]; 
  incr[4] = (-1.677050983124842*bzr[13]*c)+1.677050983124842*bzl[13]*c+1.479019945774904*bzr[8]*c+1.479019945774904*bzl[8]*c-1.25*bzr[4]*c+1.25*bzl[4]*c+0.9682458365518543*bzr[1]*c+0.9682458365518543*bzl[1]*c-0.5590169943749475*bzr[0]*c+0.5590169943749475*bzl[0]*c+1.677050983124842*eyr[13]+1.677050983124842*eyl[13]-1.479019945774904*eyr[8]+1.479019945774904*eyl[8]+1.25*eyr[4]+1.25*eyl[4]-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0]; 
  incr[5] = (-0.5590169943749475*bzr[10]*c)+0.5590169943749475*bzl[10]*c+0.4330127018922194*bzr[7]*c+0.4330127018922194*bzl[7]*c-0.25*bzr[5]*c+0.25*bzl[5]*c+0.5590169943749475*eyr[10]+0.5590169943749475*eyl[10]-0.4330127018922194*eyr[7]+0.4330127018922194*eyl[7]+0.25*eyr[5]+0.25*eyl[5]; 
  incr[6] = 1.479019945774904*bzr[11]*c+1.479019945774904*bzl[11]*c-1.25*bzr[6]*c+1.25*bzl[6]*c+0.9682458365518543*bzr[3]*c+0.9682458365518543*bzl[3]*c-0.5590169943749476*bzr[2]*c+0.5590169943749476*bzl[2]*c-1.479019945774904*eyr[11]+1.479019945774904*eyl[11]+1.25*eyr[6]+1.25*eyl[6]-0.9682458365518543*eyr[3]+0.9682458365518543*eyl[3]+0.5590169943749476*eyr[2]+0.5590169943749476*eyl[2]; 
  incr[7] = 0.9682458365518543*bzr[10]*c-0.9682458365518543*bzl[10]*c-0.75*bzr[7]*c-0.75*bzl[7]*c+0.4330127018922194*bzr[5]*c-0.4330127018922194*bzl[5]*c-0.9682458365518543*eyr[10]-0.9682458365518543*eyl[10]+0.75*eyr[7]-0.75*eyl[7]-0.4330127018922194*eyr[5]-0.4330127018922194*eyl[5]; 
  incr[8] = 1.984313483298443*bzr[13]*c-1.984313483298443*bzl[13]*c-1.75*bzr[8]*c-1.75*bzl[8]*c+1.479019945774904*bzr[4]*c-1.479019945774904*bzl[4]*c-1.14564392373896*bzr[1]*c-1.14564392373896*bzl[1]*c+0.6614378277661477*bzr[0]*c-0.6614378277661477*bzl[0]*c-1.984313483298443*eyr[13]-1.984313483298443*eyl[13]+1.75*eyr[8]-1.75*eyl[8]-1.479019945774904*eyr[4]-1.479019945774904*eyl[4]+1.14564392373896*eyr[1]-1.14564392373896*eyl[1]-0.6614378277661477*eyr[0]-0.6614378277661477*eyl[0]; 
  incr[9] = 0.4330127018922193*bzr[12]*c+0.4330127018922193*bzl[12]*c-0.25*bzr[9]*c+0.25*bzl[9]*c-0.4330127018922193*eyr[12]+0.4330127018922193*eyl[12]+0.25*eyr[9]+0.25*eyl[9]; 
  incr[10] = (-1.25*bzr[10]*c)+1.25*bzl[10]*c+0.9682458365518543*bzr[7]*c+0.9682458365518543*bzl[7]*c-0.5590169943749475*bzr[5]*c+0.5590169943749475*bzl[5]*c+1.25*eyr[10]+1.25*eyl[10]-0.9682458365518543*eyr[7]+0.9682458365518543*eyl[7]+0.5590169943749475*eyr[5]+0.5590169943749475*eyl[5]; 
  incr[11] = (-1.75*bzr[11]*c)-1.75*bzl[11]*c+1.479019945774904*bzr[6]*c-1.479019945774904*bzl[6]*c-1.14564392373896*bzr[3]*c-1.14564392373896*bzl[3]*c+0.6614378277661477*bzr[2]*c-0.6614378277661477*bzl[2]*c+1.75*eyr[11]-1.75*eyl[11]-1.479019945774904*eyr[6]-1.479019945774904*eyl[6]+1.14564392373896*eyr[3]-1.14564392373896*eyl[3]-0.6614378277661477*eyr[2]-0.6614378277661477*eyl[2]; 
  incr[12] = (-0.75*bzr[12]*c)-0.75*bzl[12]*c+0.4330127018922193*bzr[9]*c-0.4330127018922193*bzl[9]*c+0.75*eyr[12]-0.75*eyl[12]-0.4330127018922193*eyr[9]-0.4330127018922193*eyl[9]; 
  incr[13] = (-2.25*bzr[13]*c)+2.25*bzl[13]*c+1.984313483298443*bzr[8]*c+1.984313483298443*bzl[8]*c-1.677050983124842*bzr[4]*c+1.677050983124842*bzl[4]*c+1.299038105676658*bzr[1]*c+1.299038105676658*bzl[1]*c-0.75*bzr[0]*c+0.75*bzl[0]*c+2.25*eyr[13]+2.25*eyl[13]-1.984313483298443*eyr[8]+1.984313483298443*eyl[8]+1.677050983124842*eyr[4]+1.677050983124842*eyl[4]-1.299038105676658*eyr[1]+1.299038105676658*eyl[1]+0.75*eyr[0]+0.75*eyl[0]; 
  incr[14] = (-0.25*bzr[14]*c)+0.25*bzl[14]*c+0.25*eyr[14]+0.25*eyl[14]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 
  outBzr[3] += incr[3]*dx1; 
  outBzr[4] += incr[4]*dx1; 
  outBzr[5] += incr[5]*dx1; 
  outBzr[6] += incr[6]*dx1; 
  outBzr[7] += incr[7]*dx1; 
  outBzr[8] += incr[8]*dx1; 
  outBzr[9] += incr[9]*dx1; 
  outBzr[10] += incr[10]*dx1; 
  outBzr[11] += incr[11]*dx1; 
  outBzr[12] += incr[12]*dx1; 
  outBzr[13] += incr[13]*dx1; 
  outBzr[14] += incr[14]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += incr[1]*dx1; 
  outBzl[2] += -1.0*incr[2]*dx1; 
  outBzl[3] += incr[3]*dx1; 
  outBzl[4] += -1.0*incr[4]*dx1; 
  outBzl[5] += -1.0*incr[5]*dx1; 
  outBzl[6] += -1.0*incr[6]*dx1; 
  outBzl[7] += incr[7]*dx1; 
  outBzl[8] += incr[8]*dx1; 
  outBzl[9] += -1.0*incr[9]*dx1; 
  outBzl[10] += -1.0*incr[10]*dx1; 
  outBzl[11] += incr[11]*dx1; 
  outBzl[12] += incr[12]*dx1; 
  outBzl[13] += -1.0*incr[13]*dx1; 
  outBzl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = (-0.75*phr[13]*c*chi)+0.75*phl[13]*c*chi+0.6614378277661477*phr[8]*c*chi+0.6614378277661477*phl[8]*c*chi-0.5590169943749475*phr[4]*c*chi+0.5590169943749475*phl[4]*c*chi+0.4330127018922193*phr[1]*c*chi+0.4330127018922193*phl[1]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi+0.75*exr[13]*chi+0.75*exl[13]*chi-0.6614378277661477*exr[8]*chi+0.6614378277661477*exl[8]*chi+0.5590169943749475*exr[4]*chi+0.5590169943749475*exl[4]*chi-0.4330127018922193*exr[1]*chi+0.4330127018922193*exl[1]*chi+0.25*exr[0]*chi+0.25*exl[0]*chi; 
  incr[1] = 1.299038105676658*phr[13]*c*chi-1.299038105676658*phl[13]*c*chi-1.14564392373896*phr[8]*c*chi-1.14564392373896*phl[8]*c*chi+0.9682458365518543*phr[4]*c*chi-0.9682458365518543*phl[4]*c*chi-0.75*phr[1]*c*chi-0.75*phl[1]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi-1.299038105676658*exr[13]*chi-1.299038105676658*exl[13]*chi+1.14564392373896*exr[8]*chi-1.14564392373896*exl[8]*chi-0.9682458365518543*exr[4]*chi-0.9682458365518543*exl[4]*chi+0.75*exr[1]*chi-0.75*exl[1]*chi-0.4330127018922193*exr[0]*chi-0.4330127018922193*exl[0]*chi; 
  incr[2] = 0.6614378277661477*phr[11]*c*chi+0.6614378277661477*phl[11]*c*chi-0.5590169943749476*phr[6]*c*chi+0.5590169943749476*phl[6]*c*chi+0.4330127018922193*phr[3]*c*chi+0.4330127018922193*phl[3]*c*chi-0.25*phr[2]*c*chi+0.25*phl[2]*c*chi-0.6614378277661477*exr[11]*chi+0.6614378277661477*exl[11]*chi+0.5590169943749476*exr[6]*chi+0.5590169943749476*exl[6]*chi-0.4330127018922193*exr[3]*chi+0.4330127018922193*exl[3]*chi+0.25*exr[2]*chi+0.25*exl[2]*chi; 
  incr[3] = (-1.14564392373896*phr[11]*c*chi)-1.14564392373896*phl[11]*c*chi+0.9682458365518543*phr[6]*c*chi-0.9682458365518543*phl[6]*c*chi-0.75*phr[3]*c*chi-0.75*phl[3]*c*chi+0.4330127018922193*phr[2]*c*chi-0.4330127018922193*phl[2]*c*chi+1.14564392373896*exr[11]*chi-1.14564392373896*exl[11]*chi-0.9682458365518543*exr[6]*chi-0.9682458365518543*exl[6]*chi+0.75*exr[3]*chi-0.75*exl[3]*chi-0.4330127018922193*exr[2]*chi-0.4330127018922193*exl[2]*chi; 
  incr[4] = (-1.677050983124842*phr[13]*c*chi)+1.677050983124842*phl[13]*c*chi+1.479019945774904*phr[8]*c*chi+1.479019945774904*phl[8]*c*chi-1.25*phr[4]*c*chi+1.25*phl[4]*c*chi+0.9682458365518543*phr[1]*c*chi+0.9682458365518543*phl[1]*c*chi-0.5590169943749475*phr[0]*c*chi+0.5590169943749475*phl[0]*c*chi+1.677050983124842*exr[13]*chi+1.677050983124842*exl[13]*chi-1.479019945774904*exr[8]*chi+1.479019945774904*exl[8]*chi+1.25*exr[4]*chi+1.25*exl[4]*chi-0.9682458365518543*exr[1]*chi+0.9682458365518543*exl[1]*chi+0.5590169943749475*exr[0]*chi+0.5590169943749475*exl[0]*chi; 
  incr[5] = (-0.5590169943749475*phr[10]*c*chi)+0.5590169943749475*phl[10]*c*chi+0.4330127018922194*phr[7]*c*chi+0.4330127018922194*phl[7]*c*chi-0.25*phr[5]*c*chi+0.25*phl[5]*c*chi+0.5590169943749475*exr[10]*chi+0.5590169943749475*exl[10]*chi-0.4330127018922194*exr[7]*chi+0.4330127018922194*exl[7]*chi+0.25*exr[5]*chi+0.25*exl[5]*chi; 
  incr[6] = 1.479019945774904*phr[11]*c*chi+1.479019945774904*phl[11]*c*chi-1.25*phr[6]*c*chi+1.25*phl[6]*c*chi+0.9682458365518543*phr[3]*c*chi+0.9682458365518543*phl[3]*c*chi-0.5590169943749476*phr[2]*c*chi+0.5590169943749476*phl[2]*c*chi-1.479019945774904*exr[11]*chi+1.479019945774904*exl[11]*chi+1.25*exr[6]*chi+1.25*exl[6]*chi-0.9682458365518543*exr[3]*chi+0.9682458365518543*exl[3]*chi+0.5590169943749476*exr[2]*chi+0.5590169943749476*exl[2]*chi; 
  incr[7] = 0.9682458365518543*phr[10]*c*chi-0.9682458365518543*phl[10]*c*chi-0.75*phr[7]*c*chi-0.75*phl[7]*c*chi+0.4330127018922194*phr[5]*c*chi-0.4330127018922194*phl[5]*c*chi-0.9682458365518543*exr[10]*chi-0.9682458365518543*exl[10]*chi+0.75*exr[7]*chi-0.75*exl[7]*chi-0.4330127018922194*exr[5]*chi-0.4330127018922194*exl[5]*chi; 
  incr[8] = 1.984313483298443*phr[13]*c*chi-1.984313483298443*phl[13]*c*chi-1.75*phr[8]*c*chi-1.75*phl[8]*c*chi+1.479019945774904*phr[4]*c*chi-1.479019945774904*phl[4]*c*chi-1.14564392373896*phr[1]*c*chi-1.14564392373896*phl[1]*c*chi+0.6614378277661477*phr[0]*c*chi-0.6614378277661477*phl[0]*c*chi-1.984313483298443*exr[13]*chi-1.984313483298443*exl[13]*chi+1.75*exr[8]*chi-1.75*exl[8]*chi-1.479019945774904*exr[4]*chi-1.479019945774904*exl[4]*chi+1.14564392373896*exr[1]*chi-1.14564392373896*exl[1]*chi-0.6614378277661477*exr[0]*chi-0.6614378277661477*exl[0]*chi; 
  incr[9] = 0.4330127018922193*phr[12]*c*chi+0.4330127018922193*phl[12]*c*chi-0.25*phr[9]*c*chi+0.25*phl[9]*c*chi-0.4330127018922193*exr[12]*chi+0.4330127018922193*exl[12]*chi+0.25*exr[9]*chi+0.25*exl[9]*chi; 
  incr[10] = (-1.25*phr[10]*c*chi)+1.25*phl[10]*c*chi+0.9682458365518543*phr[7]*c*chi+0.9682458365518543*phl[7]*c*chi-0.5590169943749475*phr[5]*c*chi+0.5590169943749475*phl[5]*c*chi+1.25*exr[10]*chi+1.25*exl[10]*chi-0.9682458365518543*exr[7]*chi+0.9682458365518543*exl[7]*chi+0.5590169943749475*exr[5]*chi+0.5590169943749475*exl[5]*chi; 
  incr[11] = (-1.75*phr[11]*c*chi)-1.75*phl[11]*c*chi+1.479019945774904*phr[6]*c*chi-1.479019945774904*phl[6]*c*chi-1.14564392373896*phr[3]*c*chi-1.14564392373896*phl[3]*c*chi+0.6614378277661477*phr[2]*c*chi-0.6614378277661477*phl[2]*c*chi+1.75*exr[11]*chi-1.75*exl[11]*chi-1.479019945774904*exr[6]*chi-1.479019945774904*exl[6]*chi+1.14564392373896*exr[3]*chi-1.14564392373896*exl[3]*chi-0.6614378277661477*exr[2]*chi-0.6614378277661477*exl[2]*chi; 
  incr[12] = (-0.75*phr[12]*c*chi)-0.75*phl[12]*c*chi+0.4330127018922193*phr[9]*c*chi-0.4330127018922193*phl[9]*c*chi+0.75*exr[12]*chi-0.75*exl[12]*chi-0.4330127018922193*exr[9]*chi-0.4330127018922193*exl[9]*chi; 
  incr[13] = (-2.25*phr[13]*c*chi)+2.25*phl[13]*c*chi+1.984313483298443*phr[8]*c*chi+1.984313483298443*phl[8]*c*chi-1.677050983124842*phr[4]*c*chi+1.677050983124842*phl[4]*c*chi+1.299038105676658*phr[1]*c*chi+1.299038105676658*phl[1]*c*chi-0.75*phr[0]*c*chi+0.75*phl[0]*c*chi+2.25*exr[13]*chi+2.25*exl[13]*chi-1.984313483298443*exr[8]*chi+1.984313483298443*exl[8]*chi+1.677050983124842*exr[4]*chi+1.677050983124842*exl[4]*chi-1.299038105676658*exr[1]*chi+1.299038105676658*exl[1]*chi+0.75*exr[0]*chi+0.75*exl[0]*chi; 
  incr[14] = (-0.25*phr[14]*c*chi)+0.25*phl[14]*c*chi+0.25*exr[14]*chi+0.25*exl[14]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 
  outPhr[3] += incr[3]*dx1; 
  outPhr[4] += incr[4]*dx1; 
  outPhr[5] += incr[5]*dx1; 
  outPhr[6] += incr[6]*dx1; 
  outPhr[7] += incr[7]*dx1; 
  outPhr[8] += incr[8]*dx1; 
  outPhr[9] += incr[9]*dx1; 
  outPhr[10] += incr[10]*dx1; 
  outPhr[11] += incr[11]*dx1; 
  outPhr[12] += incr[12]*dx1; 
  outPhr[13] += incr[13]*dx1; 
  outPhr[14] += incr[14]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += incr[1]*dx1; 
  outPhl[2] += -1.0*incr[2]*dx1; 
  outPhl[3] += incr[3]*dx1; 
  outPhl[4] += -1.0*incr[4]*dx1; 
  outPhl[5] += -1.0*incr[5]*dx1; 
  outPhl[6] += -1.0*incr[6]*dx1; 
  outPhl[7] += incr[7]*dx1; 
  outPhl[8] += incr[8]*dx1; 
  outPhl[9] += -1.0*incr[9]*dx1; 
  outPhl[10] += -1.0*incr[10]*dx1; 
  outPhl[11] += incr[11]*dx1; 
  outPhl[12] += incr[12]*dx1; 
  outPhl[13] += -1.0*incr[13]*dx1; 
  outPhl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = 0.75*bxr[13]*c2*gamma+0.75*bxl[13]*c2*gamma-0.6614378277661477*bxr[8]*c2*gamma+0.6614378277661477*bxl[8]*c2*gamma+0.5590169943749475*bxr[4]*c2*gamma+0.5590169943749475*bxl[4]*c2*gamma-0.4330127018922193*bxr[1]*c2*gamma+0.4330127018922193*bxl[1]*c2*gamma+0.25*bxr[0]*c2*gamma+0.25*bxl[0]*c2*gamma-0.75*psr[13]*c*gamma+0.75*psl[13]*c*gamma+0.6614378277661477*psr[8]*c*gamma+0.6614378277661477*psl[8]*c*gamma-0.5590169943749475*psr[4]*c*gamma+0.5590169943749475*psl[4]*c*gamma+0.4330127018922193*psr[1]*c*gamma+0.4330127018922193*psl[1]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = (-1.299038105676658*bxr[13]*c2*gamma)-1.299038105676658*bxl[13]*c2*gamma+1.14564392373896*bxr[8]*c2*gamma-1.14564392373896*bxl[8]*c2*gamma-0.9682458365518543*bxr[4]*c2*gamma-0.9682458365518543*bxl[4]*c2*gamma+0.75*bxr[1]*c2*gamma-0.75*bxl[1]*c2*gamma-0.4330127018922193*bxr[0]*c2*gamma-0.4330127018922193*bxl[0]*c2*gamma+1.299038105676658*psr[13]*c*gamma-1.299038105676658*psl[13]*c*gamma-1.14564392373896*psr[8]*c*gamma-1.14564392373896*psl[8]*c*gamma+0.9682458365518543*psr[4]*c*gamma-0.9682458365518543*psl[4]*c*gamma-0.75*psr[1]*c*gamma-0.75*psl[1]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 
  incr[2] = (-0.6614378277661477*bxr[11]*c2*gamma)+0.6614378277661477*bxl[11]*c2*gamma+0.5590169943749476*bxr[6]*c2*gamma+0.5590169943749476*bxl[6]*c2*gamma-0.4330127018922193*bxr[3]*c2*gamma+0.4330127018922193*bxl[3]*c2*gamma+0.25*bxr[2]*c2*gamma+0.25*bxl[2]*c2*gamma+0.6614378277661477*psr[11]*c*gamma+0.6614378277661477*psl[11]*c*gamma-0.5590169943749476*psr[6]*c*gamma+0.5590169943749476*psl[6]*c*gamma+0.4330127018922193*psr[3]*c*gamma+0.4330127018922193*psl[3]*c*gamma-0.25*psr[2]*c*gamma+0.25*psl[2]*c*gamma; 
  incr[3] = 1.14564392373896*bxr[11]*c2*gamma-1.14564392373896*bxl[11]*c2*gamma-0.9682458365518543*bxr[6]*c2*gamma-0.9682458365518543*bxl[6]*c2*gamma+0.75*bxr[3]*c2*gamma-0.75*bxl[3]*c2*gamma-0.4330127018922193*bxr[2]*c2*gamma-0.4330127018922193*bxl[2]*c2*gamma-1.14564392373896*psr[11]*c*gamma-1.14564392373896*psl[11]*c*gamma+0.9682458365518543*psr[6]*c*gamma-0.9682458365518543*psl[6]*c*gamma-0.75*psr[3]*c*gamma-0.75*psl[3]*c*gamma+0.4330127018922193*psr[2]*c*gamma-0.4330127018922193*psl[2]*c*gamma; 
  incr[4] = 1.677050983124842*bxr[13]*c2*gamma+1.677050983124842*bxl[13]*c2*gamma-1.479019945774904*bxr[8]*c2*gamma+1.479019945774904*bxl[8]*c2*gamma+1.25*bxr[4]*c2*gamma+1.25*bxl[4]*c2*gamma-0.9682458365518543*bxr[1]*c2*gamma+0.9682458365518543*bxl[1]*c2*gamma+0.5590169943749475*bxr[0]*c2*gamma+0.5590169943749475*bxl[0]*c2*gamma-1.677050983124842*psr[13]*c*gamma+1.677050983124842*psl[13]*c*gamma+1.479019945774904*psr[8]*c*gamma+1.479019945774904*psl[8]*c*gamma-1.25*psr[4]*c*gamma+1.25*psl[4]*c*gamma+0.9682458365518543*psr[1]*c*gamma+0.9682458365518543*psl[1]*c*gamma-0.5590169943749475*psr[0]*c*gamma+0.5590169943749475*psl[0]*c*gamma; 
  incr[5] = 0.5590169943749475*bxr[10]*c2*gamma+0.5590169943749475*bxl[10]*c2*gamma-0.4330127018922194*bxr[7]*c2*gamma+0.4330127018922194*bxl[7]*c2*gamma+0.25*bxr[5]*c2*gamma+0.25*bxl[5]*c2*gamma-0.5590169943749475*psr[10]*c*gamma+0.5590169943749475*psl[10]*c*gamma+0.4330127018922194*psr[7]*c*gamma+0.4330127018922194*psl[7]*c*gamma-0.25*psr[5]*c*gamma+0.25*psl[5]*c*gamma; 
  incr[6] = (-1.479019945774904*bxr[11]*c2*gamma)+1.479019945774904*bxl[11]*c2*gamma+1.25*bxr[6]*c2*gamma+1.25*bxl[6]*c2*gamma-0.9682458365518543*bxr[3]*c2*gamma+0.9682458365518543*bxl[3]*c2*gamma+0.5590169943749476*bxr[2]*c2*gamma+0.5590169943749476*bxl[2]*c2*gamma+1.479019945774904*psr[11]*c*gamma+1.479019945774904*psl[11]*c*gamma-1.25*psr[6]*c*gamma+1.25*psl[6]*c*gamma+0.9682458365518543*psr[3]*c*gamma+0.9682458365518543*psl[3]*c*gamma-0.5590169943749476*psr[2]*c*gamma+0.5590169943749476*psl[2]*c*gamma; 
  incr[7] = (-0.9682458365518543*bxr[10]*c2*gamma)-0.9682458365518543*bxl[10]*c2*gamma+0.75*bxr[7]*c2*gamma-0.75*bxl[7]*c2*gamma-0.4330127018922194*bxr[5]*c2*gamma-0.4330127018922194*bxl[5]*c2*gamma+0.9682458365518543*psr[10]*c*gamma-0.9682458365518543*psl[10]*c*gamma-0.75*psr[7]*c*gamma-0.75*psl[7]*c*gamma+0.4330127018922194*psr[5]*c*gamma-0.4330127018922194*psl[5]*c*gamma; 
  incr[8] = (-1.984313483298443*bxr[13]*c2*gamma)-1.984313483298443*bxl[13]*c2*gamma+1.75*bxr[8]*c2*gamma-1.75*bxl[8]*c2*gamma-1.479019945774904*bxr[4]*c2*gamma-1.479019945774904*bxl[4]*c2*gamma+1.14564392373896*bxr[1]*c2*gamma-1.14564392373896*bxl[1]*c2*gamma-0.6614378277661477*bxr[0]*c2*gamma-0.6614378277661477*bxl[0]*c2*gamma+1.984313483298443*psr[13]*c*gamma-1.984313483298443*psl[13]*c*gamma-1.75*psr[8]*c*gamma-1.75*psl[8]*c*gamma+1.479019945774904*psr[4]*c*gamma-1.479019945774904*psl[4]*c*gamma-1.14564392373896*psr[1]*c*gamma-1.14564392373896*psl[1]*c*gamma+0.6614378277661477*psr[0]*c*gamma-0.6614378277661477*psl[0]*c*gamma; 
  incr[9] = (-0.4330127018922193*bxr[12]*c2*gamma)+0.4330127018922193*bxl[12]*c2*gamma+0.25*bxr[9]*c2*gamma+0.25*bxl[9]*c2*gamma+0.4330127018922193*psr[12]*c*gamma+0.4330127018922193*psl[12]*c*gamma-0.25*psr[9]*c*gamma+0.25*psl[9]*c*gamma; 
  incr[10] = 1.25*bxr[10]*c2*gamma+1.25*bxl[10]*c2*gamma-0.9682458365518543*bxr[7]*c2*gamma+0.9682458365518543*bxl[7]*c2*gamma+0.5590169943749475*bxr[5]*c2*gamma+0.5590169943749475*bxl[5]*c2*gamma-1.25*psr[10]*c*gamma+1.25*psl[10]*c*gamma+0.9682458365518543*psr[7]*c*gamma+0.9682458365518543*psl[7]*c*gamma-0.5590169943749475*psr[5]*c*gamma+0.5590169943749475*psl[5]*c*gamma; 
  incr[11] = 1.75*bxr[11]*c2*gamma-1.75*bxl[11]*c2*gamma-1.479019945774904*bxr[6]*c2*gamma-1.479019945774904*bxl[6]*c2*gamma+1.14564392373896*bxr[3]*c2*gamma-1.14564392373896*bxl[3]*c2*gamma-0.6614378277661477*bxr[2]*c2*gamma-0.6614378277661477*bxl[2]*c2*gamma-1.75*psr[11]*c*gamma-1.75*psl[11]*c*gamma+1.479019945774904*psr[6]*c*gamma-1.479019945774904*psl[6]*c*gamma-1.14564392373896*psr[3]*c*gamma-1.14564392373896*psl[3]*c*gamma+0.6614378277661477*psr[2]*c*gamma-0.6614378277661477*psl[2]*c*gamma; 
  incr[12] = 0.75*bxr[12]*c2*gamma-0.75*bxl[12]*c2*gamma-0.4330127018922193*bxr[9]*c2*gamma-0.4330127018922193*bxl[9]*c2*gamma-0.75*psr[12]*c*gamma-0.75*psl[12]*c*gamma+0.4330127018922193*psr[9]*c*gamma-0.4330127018922193*psl[9]*c*gamma; 
  incr[13] = 2.25*bxr[13]*c2*gamma+2.25*bxl[13]*c2*gamma-1.984313483298443*bxr[8]*c2*gamma+1.984313483298443*bxl[8]*c2*gamma+1.677050983124842*bxr[4]*c2*gamma+1.677050983124842*bxl[4]*c2*gamma-1.299038105676658*bxr[1]*c2*gamma+1.299038105676658*bxl[1]*c2*gamma+0.75*bxr[0]*c2*gamma+0.75*bxl[0]*c2*gamma-2.25*psr[13]*c*gamma+2.25*psl[13]*c*gamma+1.984313483298443*psr[8]*c*gamma+1.984313483298443*psl[8]*c*gamma-1.677050983124842*psr[4]*c*gamma+1.677050983124842*psl[4]*c*gamma+1.299038105676658*psr[1]*c*gamma+1.299038105676658*psl[1]*c*gamma-0.75*psr[0]*c*gamma+0.75*psl[0]*c*gamma; 
  incr[14] = 0.25*bxr[14]*c2*gamma+0.25*bxl[14]*c2*gamma-0.25*psr[14]*c*gamma+0.25*psl[14]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 
  outPsr[3] += incr[3]*dx1; 
  outPsr[4] += incr[4]*dx1; 
  outPsr[5] += incr[5]*dx1; 
  outPsr[6] += incr[6]*dx1; 
  outPsr[7] += incr[7]*dx1; 
  outPsr[8] += incr[8]*dx1; 
  outPsr[9] += incr[9]*dx1; 
  outPsr[10] += incr[10]*dx1; 
  outPsr[11] += incr[11]*dx1; 
  outPsr[12] += incr[12]*dx1; 
  outPsr[13] += incr[13]*dx1; 
  outPsr[14] += incr[14]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += incr[1]*dx1; 
  outPsl[2] += -1.0*incr[2]*dx1; 
  outPsl[3] += incr[3]*dx1; 
  outPsl[4] += -1.0*incr[4]*dx1; 
  outPsl[5] += -1.0*incr[5]*dx1; 
  outPsl[6] += -1.0*incr[6]*dx1; 
  outPsl[7] += incr[7]*dx1; 
  outPsl[8] += incr[8]*dx1; 
  outPsl[9] += -1.0*incr[9]*dx1; 
  outPsl[10] += -1.0*incr[10]*dx1; 
  outPsl[11] += incr[11]*dx1; 
  outPsl[12] += incr[12]*dx1; 
  outPsl[13] += -1.0*incr[13]*dx1; 
  outPsl[14] += -1.0*incr[14]*dx1; 

 
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
void MaxwellSurf2xMax_Y_P3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[1]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[10]; 
  const double *ezl = &ql[20]; 
  const double *bxl = &ql[30]; 
  const double *byl = &ql[40]; 
  const double *bzl = &ql[50]; 
  const double *phl = &ql[60]; 
  const double *psl = &ql[70]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[10]; 
  double *outEzl = &outl[20]; 
  double *outBxl = &outl[30]; 
  double *outByl = &outl[40]; 
  double *outBzl = &outl[50]; 
  double *outPhl = &outl[60]; 
  double *outPsl = &outl[70]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[10]; 
  const double *ezr = &qr[20]; 
  const double *bxr = &qr[30]; 
  const double *byr = &qr[40]; 
  const double *bzr = &qr[50]; 
  const double *phr = &qr[60]; 
  const double *psr = &qr[70]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[10]; 
  double *outEzr = &outr[20]; 
  double *outBxr = &outr[30]; 
  double *outByr = &outr[40]; 
  double *outBzr = &outr[50]; 
  double *outPhr = &outr[60]; 
  double *outPsr = &outr[70]; 
 
  double incr[10]; 
 
  incr[0] = 0.6614378277661477*bzr[9]*c2-0.6614378277661477*bzl[9]*c2-0.5590169943749475*bzr[5]*c2-0.5590169943749475*bzl[5]*c2+0.4330127018922193*bzr[2]*c2-0.4330127018922193*bzl[2]*c2-0.25*bzr[0]*c2-0.25*bzl[0]*c2+0.6614378277661477*exr[9]*c+0.6614378277661477*exl[9]*c-0.5590169943749475*exr[5]*c+0.5590169943749475*exl[5]*c+0.4330127018922193*exr[2]*c+0.4330127018922193*exl[2]*c-0.25*exr[0]*c+0.25*exl[0]*c; 
  incr[1] = (-0.5590169943749476*bzr[7]*c2)-0.5590169943749476*bzl[7]*c2+0.4330127018922193*bzr[3]*c2-0.4330127018922193*bzl[3]*c2-0.25*bzr[1]*c2-0.25*bzl[1]*c2-0.5590169943749476*exr[7]*c+0.5590169943749476*exl[7]*c+0.4330127018922193*exr[3]*c+0.4330127018922193*exl[3]*c-0.25*exr[1]*c+0.25*exl[1]*c; 
  incr[2] = (-1.14564392373896*bzr[9]*c2)+1.14564392373896*bzl[9]*c2+0.9682458365518543*bzr[5]*c2+0.9682458365518543*bzl[5]*c2-0.75*bzr[2]*c2+0.75*bzl[2]*c2+0.4330127018922193*bzr[0]*c2+0.4330127018922193*bzl[0]*c2-1.14564392373896*exr[9]*c-1.14564392373896*exl[9]*c+0.9682458365518543*exr[5]*c-0.9682458365518543*exl[5]*c-0.75*exr[2]*c-0.75*exl[2]*c+0.4330127018922193*exr[0]*c-0.4330127018922193*exl[0]*c; 
  incr[3] = 0.9682458365518543*bzr[7]*c2+0.9682458365518543*bzl[7]*c2-0.75*bzr[3]*c2+0.75*bzl[3]*c2+0.4330127018922193*bzr[1]*c2+0.4330127018922193*bzl[1]*c2+0.9682458365518543*exr[7]*c-0.9682458365518543*exl[7]*c-0.75*exr[3]*c-0.75*exl[3]*c+0.4330127018922193*exr[1]*c-0.4330127018922193*exl[1]*c; 
  incr[4] = 0.4330127018922194*bzr[6]*c2-0.4330127018922194*bzl[6]*c2-0.25*bzr[4]*c2-0.25*bzl[4]*c2+0.4330127018922194*exr[6]*c+0.4330127018922194*exl[6]*c-0.25*exr[4]*c+0.25*exl[4]*c; 
  incr[5] = 1.479019945774904*bzr[9]*c2-1.479019945774904*bzl[9]*c2-1.25*bzr[5]*c2-1.25*bzl[5]*c2+0.9682458365518543*bzr[2]*c2-0.9682458365518543*bzl[2]*c2-0.5590169943749475*bzr[0]*c2-0.5590169943749475*bzl[0]*c2+1.479019945774904*exr[9]*c+1.479019945774904*exl[9]*c-1.25*exr[5]*c+1.25*exl[5]*c+0.9682458365518543*exr[2]*c+0.9682458365518543*exl[2]*c-0.5590169943749475*exr[0]*c+0.5590169943749475*exl[0]*c; 
  incr[6] = (-0.75*bzr[6]*c2)+0.75*bzl[6]*c2+0.4330127018922194*bzr[4]*c2+0.4330127018922194*bzl[4]*c2-0.75*exr[6]*c-0.75*exl[6]*c+0.4330127018922194*exr[4]*c-0.4330127018922194*exl[4]*c; 
  incr[7] = (-1.25*bzr[7]*c2)-1.25*bzl[7]*c2+0.9682458365518543*bzr[3]*c2-0.9682458365518543*bzl[3]*c2-0.5590169943749476*bzr[1]*c2-0.5590169943749476*bzl[1]*c2-1.25*exr[7]*c+1.25*exl[7]*c+0.9682458365518543*exr[3]*c+0.9682458365518543*exl[3]*c-0.5590169943749476*exr[1]*c+0.5590169943749476*exl[1]*c; 
  incr[8] = (-0.25*bzr[8]*c2)-0.25*bzl[8]*c2-0.25*exr[8]*c+0.25*exl[8]*c; 
  incr[9] = (-1.75*bzr[9]*c2)+1.75*bzl[9]*c2+1.479019945774904*bzr[5]*c2+1.479019945774904*bzl[5]*c2-1.14564392373896*bzr[2]*c2+1.14564392373896*bzl[2]*c2+0.6614378277661477*bzr[0]*c2+0.6614378277661477*bzl[0]*c2-1.75*exr[9]*c-1.75*exl[9]*c+1.479019945774904*exr[5]*c-1.479019945774904*exl[5]*c-1.14564392373896*exr[2]*c-1.14564392373896*exl[2]*c+0.6614378277661477*exr[0]*c-0.6614378277661477*exl[0]*c; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 
  outExr[3] += incr[3]*dx1; 
  outExr[4] += incr[4]*dx1; 
  outExr[5] += incr[5]*dx1; 
  outExr[6] += incr[6]*dx1; 
  outExr[7] += incr[7]*dx1; 
  outExr[8] += incr[8]*dx1; 
  outExr[9] += incr[9]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += -1.0*incr[1]*dx1; 
  outExl[2] += incr[2]*dx1; 
  outExl[3] += incr[3]*dx1; 
  outExl[4] += -1.0*incr[4]*dx1; 
  outExl[5] += -1.0*incr[5]*dx1; 
  outExl[6] += incr[6]*dx1; 
  outExl[7] += -1.0*incr[7]*dx1; 
  outExl[8] += -1.0*incr[8]*dx1; 
  outExl[9] += incr[9]*dx1; 

 
  incr[0] = (-0.6614378277661477*phr[9]*c2*chi)+0.6614378277661477*phl[9]*c2*chi+0.5590169943749475*phr[5]*c2*chi+0.5590169943749475*phl[5]*c2*chi-0.4330127018922193*phr[2]*c2*chi+0.4330127018922193*phl[2]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi+0.6614378277661477*eyr[9]*c*chi+0.6614378277661477*eyl[9]*c*chi-0.5590169943749475*eyr[5]*c*chi+0.5590169943749475*eyl[5]*c*chi+0.4330127018922193*eyr[2]*c*chi+0.4330127018922193*eyl[2]*c*chi-0.25*eyr[0]*c*chi+0.25*eyl[0]*c*chi; 
  incr[1] = 0.5590169943749476*phr[7]*c2*chi+0.5590169943749476*phl[7]*c2*chi-0.4330127018922193*phr[3]*c2*chi+0.4330127018922193*phl[3]*c2*chi+0.25*phr[1]*c2*chi+0.25*phl[1]*c2*chi-0.5590169943749476*eyr[7]*c*chi+0.5590169943749476*eyl[7]*c*chi+0.4330127018922193*eyr[3]*c*chi+0.4330127018922193*eyl[3]*c*chi-0.25*eyr[1]*c*chi+0.25*eyl[1]*c*chi; 
  incr[2] = 1.14564392373896*phr[9]*c2*chi-1.14564392373896*phl[9]*c2*chi-0.9682458365518543*phr[5]*c2*chi-0.9682458365518543*phl[5]*c2*chi+0.75*phr[2]*c2*chi-0.75*phl[2]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi-1.14564392373896*eyr[9]*c*chi-1.14564392373896*eyl[9]*c*chi+0.9682458365518543*eyr[5]*c*chi-0.9682458365518543*eyl[5]*c*chi-0.75*eyr[2]*c*chi-0.75*eyl[2]*c*chi+0.4330127018922193*eyr[0]*c*chi-0.4330127018922193*eyl[0]*c*chi; 
  incr[3] = (-0.9682458365518543*phr[7]*c2*chi)-0.9682458365518543*phl[7]*c2*chi+0.75*phr[3]*c2*chi-0.75*phl[3]*c2*chi-0.4330127018922193*phr[1]*c2*chi-0.4330127018922193*phl[1]*c2*chi+0.9682458365518543*eyr[7]*c*chi-0.9682458365518543*eyl[7]*c*chi-0.75*eyr[3]*c*chi-0.75*eyl[3]*c*chi+0.4330127018922193*eyr[1]*c*chi-0.4330127018922193*eyl[1]*c*chi; 
  incr[4] = (-0.4330127018922194*phr[6]*c2*chi)+0.4330127018922194*phl[6]*c2*chi+0.25*phr[4]*c2*chi+0.25*phl[4]*c2*chi+0.4330127018922194*eyr[6]*c*chi+0.4330127018922194*eyl[6]*c*chi-0.25*eyr[4]*c*chi+0.25*eyl[4]*c*chi; 
  incr[5] = (-1.479019945774904*phr[9]*c2*chi)+1.479019945774904*phl[9]*c2*chi+1.25*phr[5]*c2*chi+1.25*phl[5]*c2*chi-0.9682458365518543*phr[2]*c2*chi+0.9682458365518543*phl[2]*c2*chi+0.5590169943749475*phr[0]*c2*chi+0.5590169943749475*phl[0]*c2*chi+1.479019945774904*eyr[9]*c*chi+1.479019945774904*eyl[9]*c*chi-1.25*eyr[5]*c*chi+1.25*eyl[5]*c*chi+0.9682458365518543*eyr[2]*c*chi+0.9682458365518543*eyl[2]*c*chi-0.5590169943749475*eyr[0]*c*chi+0.5590169943749475*eyl[0]*c*chi; 
  incr[6] = 0.75*phr[6]*c2*chi-0.75*phl[6]*c2*chi-0.4330127018922194*phr[4]*c2*chi-0.4330127018922194*phl[4]*c2*chi-0.75*eyr[6]*c*chi-0.75*eyl[6]*c*chi+0.4330127018922194*eyr[4]*c*chi-0.4330127018922194*eyl[4]*c*chi; 
  incr[7] = 1.25*phr[7]*c2*chi+1.25*phl[7]*c2*chi-0.9682458365518543*phr[3]*c2*chi+0.9682458365518543*phl[3]*c2*chi+0.5590169943749476*phr[1]*c2*chi+0.5590169943749476*phl[1]*c2*chi-1.25*eyr[7]*c*chi+1.25*eyl[7]*c*chi+0.9682458365518543*eyr[3]*c*chi+0.9682458365518543*eyl[3]*c*chi-0.5590169943749476*eyr[1]*c*chi+0.5590169943749476*eyl[1]*c*chi; 
  incr[8] = 0.25*phr[8]*c2*chi+0.25*phl[8]*c2*chi-0.25*eyr[8]*c*chi+0.25*eyl[8]*c*chi; 
  incr[9] = 1.75*phr[9]*c2*chi-1.75*phl[9]*c2*chi-1.479019945774904*phr[5]*c2*chi-1.479019945774904*phl[5]*c2*chi+1.14564392373896*phr[2]*c2*chi-1.14564392373896*phl[2]*c2*chi-0.6614378277661477*phr[0]*c2*chi-0.6614378277661477*phl[0]*c2*chi-1.75*eyr[9]*c*chi-1.75*eyl[9]*c*chi+1.479019945774904*eyr[5]*c*chi-1.479019945774904*eyl[5]*c*chi-1.14564392373896*eyr[2]*c*chi-1.14564392373896*eyl[2]*c*chi+0.6614378277661477*eyr[0]*c*chi-0.6614378277661477*eyl[0]*c*chi; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 
  outEyr[3] += incr[3]*dx1; 
  outEyr[4] += incr[4]*dx1; 
  outEyr[5] += incr[5]*dx1; 
  outEyr[6] += incr[6]*dx1; 
  outEyr[7] += incr[7]*dx1; 
  outEyr[8] += incr[8]*dx1; 
  outEyr[9] += incr[9]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += -1.0*incr[1]*dx1; 
  outEyl[2] += incr[2]*dx1; 
  outEyl[3] += incr[3]*dx1; 
  outEyl[4] += -1.0*incr[4]*dx1; 
  outEyl[5] += -1.0*incr[5]*dx1; 
  outEyl[6] += incr[6]*dx1; 
  outEyl[7] += -1.0*incr[7]*dx1; 
  outEyl[8] += -1.0*incr[8]*dx1; 
  outEyl[9] += incr[9]*dx1; 

 
  incr[0] = (-0.6614378277661477*bxr[9]*c2)+0.6614378277661477*bxl[9]*c2+0.5590169943749475*bxr[5]*c2+0.5590169943749475*bxl[5]*c2-0.4330127018922193*bxr[2]*c2+0.4330127018922193*bxl[2]*c2+0.25*bxr[0]*c2+0.25*bxl[0]*c2+0.6614378277661477*ezr[9]*c+0.6614378277661477*ezl[9]*c-0.5590169943749475*ezr[5]*c+0.5590169943749475*ezl[5]*c+0.4330127018922193*ezr[2]*c+0.4330127018922193*ezl[2]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = 0.5590169943749476*bxr[7]*c2+0.5590169943749476*bxl[7]*c2-0.4330127018922193*bxr[3]*c2+0.4330127018922193*bxl[3]*c2+0.25*bxr[1]*c2+0.25*bxl[1]*c2-0.5590169943749476*ezr[7]*c+0.5590169943749476*ezl[7]*c+0.4330127018922193*ezr[3]*c+0.4330127018922193*ezl[3]*c-0.25*ezr[1]*c+0.25*ezl[1]*c; 
  incr[2] = 1.14564392373896*bxr[9]*c2-1.14564392373896*bxl[9]*c2-0.9682458365518543*bxr[5]*c2-0.9682458365518543*bxl[5]*c2+0.75*bxr[2]*c2-0.75*bxl[2]*c2-0.4330127018922193*bxr[0]*c2-0.4330127018922193*bxl[0]*c2-1.14564392373896*ezr[9]*c-1.14564392373896*ezl[9]*c+0.9682458365518543*ezr[5]*c-0.9682458365518543*ezl[5]*c-0.75*ezr[2]*c-0.75*ezl[2]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 
  incr[3] = (-0.9682458365518543*bxr[7]*c2)-0.9682458365518543*bxl[7]*c2+0.75*bxr[3]*c2-0.75*bxl[3]*c2-0.4330127018922193*bxr[1]*c2-0.4330127018922193*bxl[1]*c2+0.9682458365518543*ezr[7]*c-0.9682458365518543*ezl[7]*c-0.75*ezr[3]*c-0.75*ezl[3]*c+0.4330127018922193*ezr[1]*c-0.4330127018922193*ezl[1]*c; 
  incr[4] = (-0.4330127018922194*bxr[6]*c2)+0.4330127018922194*bxl[6]*c2+0.25*bxr[4]*c2+0.25*bxl[4]*c2+0.4330127018922194*ezr[6]*c+0.4330127018922194*ezl[6]*c-0.25*ezr[4]*c+0.25*ezl[4]*c; 
  incr[5] = (-1.479019945774904*bxr[9]*c2)+1.479019945774904*bxl[9]*c2+1.25*bxr[5]*c2+1.25*bxl[5]*c2-0.9682458365518543*bxr[2]*c2+0.9682458365518543*bxl[2]*c2+0.5590169943749475*bxr[0]*c2+0.5590169943749475*bxl[0]*c2+1.479019945774904*ezr[9]*c+1.479019945774904*ezl[9]*c-1.25*ezr[5]*c+1.25*ezl[5]*c+0.9682458365518543*ezr[2]*c+0.9682458365518543*ezl[2]*c-0.5590169943749475*ezr[0]*c+0.5590169943749475*ezl[0]*c; 
  incr[6] = 0.75*bxr[6]*c2-0.75*bxl[6]*c2-0.4330127018922194*bxr[4]*c2-0.4330127018922194*bxl[4]*c2-0.75*ezr[6]*c-0.75*ezl[6]*c+0.4330127018922194*ezr[4]*c-0.4330127018922194*ezl[4]*c; 
  incr[7] = 1.25*bxr[7]*c2+1.25*bxl[7]*c2-0.9682458365518543*bxr[3]*c2+0.9682458365518543*bxl[3]*c2+0.5590169943749476*bxr[1]*c2+0.5590169943749476*bxl[1]*c2-1.25*ezr[7]*c+1.25*ezl[7]*c+0.9682458365518543*ezr[3]*c+0.9682458365518543*ezl[3]*c-0.5590169943749476*ezr[1]*c+0.5590169943749476*ezl[1]*c; 
  incr[8] = 0.25*bxr[8]*c2+0.25*bxl[8]*c2-0.25*ezr[8]*c+0.25*ezl[8]*c; 
  incr[9] = 1.75*bxr[9]*c2-1.75*bxl[9]*c2-1.479019945774904*bxr[5]*c2-1.479019945774904*bxl[5]*c2+1.14564392373896*bxr[2]*c2-1.14564392373896*bxl[2]*c2-0.6614378277661477*bxr[0]*c2-0.6614378277661477*bxl[0]*c2-1.75*ezr[9]*c-1.75*ezl[9]*c+1.479019945774904*ezr[5]*c-1.479019945774904*ezl[5]*c-1.14564392373896*ezr[2]*c-1.14564392373896*ezl[2]*c+0.6614378277661477*ezr[0]*c-0.6614378277661477*ezl[0]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 
  outEzr[3] += incr[3]*dx1; 
  outEzr[4] += incr[4]*dx1; 
  outEzr[5] += incr[5]*dx1; 
  outEzr[6] += incr[6]*dx1; 
  outEzr[7] += incr[7]*dx1; 
  outEzr[8] += incr[8]*dx1; 
  outEzr[9] += incr[9]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += -1.0*incr[1]*dx1; 
  outEzl[2] += incr[2]*dx1; 
  outEzl[3] += incr[3]*dx1; 
  outEzl[4] += -1.0*incr[4]*dx1; 
  outEzl[5] += -1.0*incr[5]*dx1; 
  outEzl[6] += incr[6]*dx1; 
  outEzl[7] += -1.0*incr[7]*dx1; 
  outEzl[8] += -1.0*incr[8]*dx1; 
  outEzl[9] += incr[9]*dx1; 

 
  incr[0] = 0.6614378277661477*bxr[9]*c+0.6614378277661477*bxl[9]*c-0.5590169943749475*bxr[5]*c+0.5590169943749475*bxl[5]*c+0.4330127018922193*bxr[2]*c+0.4330127018922193*bxl[2]*c-0.25*bxr[0]*c+0.25*bxl[0]*c-0.6614378277661477*ezr[9]+0.6614378277661477*ezl[9]+0.5590169943749475*ezr[5]+0.5590169943749475*ezl[5]-0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]+0.25*ezr[0]+0.25*ezl[0]; 
  incr[1] = (-0.5590169943749476*bxr[7]*c)+0.5590169943749476*bxl[7]*c+0.4330127018922193*bxr[3]*c+0.4330127018922193*bxl[3]*c-0.25*bxr[1]*c+0.25*bxl[1]*c+0.5590169943749476*ezr[7]+0.5590169943749476*ezl[7]-0.4330127018922193*ezr[3]+0.4330127018922193*ezl[3]+0.25*ezr[1]+0.25*ezl[1]; 
  incr[2] = (-1.14564392373896*bxr[9]*c)-1.14564392373896*bxl[9]*c+0.9682458365518543*bxr[5]*c-0.9682458365518543*bxl[5]*c-0.75*bxr[2]*c-0.75*bxl[2]*c+0.4330127018922193*bxr[0]*c-0.4330127018922193*bxl[0]*c+1.14564392373896*ezr[9]-1.14564392373896*ezl[9]-0.9682458365518543*ezr[5]-0.9682458365518543*ezl[5]+0.75*ezr[2]-0.75*ezl[2]-0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0]; 
  incr[3] = 0.9682458365518543*bxr[7]*c-0.9682458365518543*bxl[7]*c-0.75*bxr[3]*c-0.75*bxl[3]*c+0.4330127018922193*bxr[1]*c-0.4330127018922193*bxl[1]*c-0.9682458365518543*ezr[7]-0.9682458365518543*ezl[7]+0.75*ezr[3]-0.75*ezl[3]-0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]; 
  incr[4] = 0.4330127018922194*bxr[6]*c+0.4330127018922194*bxl[6]*c-0.25*bxr[4]*c+0.25*bxl[4]*c-0.4330127018922194*ezr[6]+0.4330127018922194*ezl[6]+0.25*ezr[4]+0.25*ezl[4]; 
  incr[5] = 1.479019945774904*bxr[9]*c+1.479019945774904*bxl[9]*c-1.25*bxr[5]*c+1.25*bxl[5]*c+0.9682458365518543*bxr[2]*c+0.9682458365518543*bxl[2]*c-0.5590169943749475*bxr[0]*c+0.5590169943749475*bxl[0]*c-1.479019945774904*ezr[9]+1.479019945774904*ezl[9]+1.25*ezr[5]+1.25*ezl[5]-0.9682458365518543*ezr[2]+0.9682458365518543*ezl[2]+0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0]; 
  incr[6] = (-0.75*bxr[6]*c)-0.75*bxl[6]*c+0.4330127018922194*bxr[4]*c-0.4330127018922194*bxl[4]*c+0.75*ezr[6]-0.75*ezl[6]-0.4330127018922194*ezr[4]-0.4330127018922194*ezl[4]; 
  incr[7] = (-1.25*bxr[7]*c)+1.25*bxl[7]*c+0.9682458365518543*bxr[3]*c+0.9682458365518543*bxl[3]*c-0.5590169943749476*bxr[1]*c+0.5590169943749476*bxl[1]*c+1.25*ezr[7]+1.25*ezl[7]-0.9682458365518543*ezr[3]+0.9682458365518543*ezl[3]+0.5590169943749476*ezr[1]+0.5590169943749476*ezl[1]; 
  incr[8] = (-0.25*bxr[8]*c)+0.25*bxl[8]*c+0.25*ezr[8]+0.25*ezl[8]; 
  incr[9] = (-1.75*bxr[9]*c)-1.75*bxl[9]*c+1.479019945774904*bxr[5]*c-1.479019945774904*bxl[5]*c-1.14564392373896*bxr[2]*c-1.14564392373896*bxl[2]*c+0.6614378277661477*bxr[0]*c-0.6614378277661477*bxl[0]*c+1.75*ezr[9]-1.75*ezl[9]-1.479019945774904*ezr[5]-1.479019945774904*ezl[5]+1.14564392373896*ezr[2]-1.14564392373896*ezl[2]-0.6614378277661477*ezr[0]-0.6614378277661477*ezl[0]; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 
  outBxr[3] += incr[3]*dx1; 
  outBxr[4] += incr[4]*dx1; 
  outBxr[5] += incr[5]*dx1; 
  outBxr[6] += incr[6]*dx1; 
  outBxr[7] += incr[7]*dx1; 
  outBxr[8] += incr[8]*dx1; 
  outBxr[9] += incr[9]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += -1.0*incr[1]*dx1; 
  outBxl[2] += incr[2]*dx1; 
  outBxl[3] += incr[3]*dx1; 
  outBxl[4] += -1.0*incr[4]*dx1; 
  outBxl[5] += -1.0*incr[5]*dx1; 
  outBxl[6] += incr[6]*dx1; 
  outBxl[7] += -1.0*incr[7]*dx1; 
  outBxl[8] += -1.0*incr[8]*dx1; 
  outBxl[9] += incr[9]*dx1; 

 
  incr[0] = 0.6614378277661477*byr[9]*c*gamma+0.6614378277661477*byl[9]*c*gamma-0.5590169943749475*byr[5]*c*gamma+0.5590169943749475*byl[5]*c*gamma+0.4330127018922193*byr[2]*c*gamma+0.4330127018922193*byl[2]*c*gamma-0.25*byr[0]*c*gamma+0.25*byl[0]*c*gamma-0.6614378277661477*psr[9]*gamma+0.6614378277661477*psl[9]*gamma+0.5590169943749475*psr[5]*gamma+0.5590169943749475*psl[5]*gamma-0.4330127018922193*psr[2]*gamma+0.4330127018922193*psl[2]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = (-0.5590169943749476*byr[7]*c*gamma)+0.5590169943749476*byl[7]*c*gamma+0.4330127018922193*byr[3]*c*gamma+0.4330127018922193*byl[3]*c*gamma-0.25*byr[1]*c*gamma+0.25*byl[1]*c*gamma+0.5590169943749476*psr[7]*gamma+0.5590169943749476*psl[7]*gamma-0.4330127018922193*psr[3]*gamma+0.4330127018922193*psl[3]*gamma+0.25*psr[1]*gamma+0.25*psl[1]*gamma; 
  incr[2] = (-1.14564392373896*byr[9]*c*gamma)-1.14564392373896*byl[9]*c*gamma+0.9682458365518543*byr[5]*c*gamma-0.9682458365518543*byl[5]*c*gamma-0.75*byr[2]*c*gamma-0.75*byl[2]*c*gamma+0.4330127018922193*byr[0]*c*gamma-0.4330127018922193*byl[0]*c*gamma+1.14564392373896*psr[9]*gamma-1.14564392373896*psl[9]*gamma-0.9682458365518543*psr[5]*gamma-0.9682458365518543*psl[5]*gamma+0.75*psr[2]*gamma-0.75*psl[2]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 
  incr[3] = 0.9682458365518543*byr[7]*c*gamma-0.9682458365518543*byl[7]*c*gamma-0.75*byr[3]*c*gamma-0.75*byl[3]*c*gamma+0.4330127018922193*byr[1]*c*gamma-0.4330127018922193*byl[1]*c*gamma-0.9682458365518543*psr[7]*gamma-0.9682458365518543*psl[7]*gamma+0.75*psr[3]*gamma-0.75*psl[3]*gamma-0.4330127018922193*psr[1]*gamma-0.4330127018922193*psl[1]*gamma; 
  incr[4] = 0.4330127018922194*byr[6]*c*gamma+0.4330127018922194*byl[6]*c*gamma-0.25*byr[4]*c*gamma+0.25*byl[4]*c*gamma-0.4330127018922194*psr[6]*gamma+0.4330127018922194*psl[6]*gamma+0.25*psr[4]*gamma+0.25*psl[4]*gamma; 
  incr[5] = 1.479019945774904*byr[9]*c*gamma+1.479019945774904*byl[9]*c*gamma-1.25*byr[5]*c*gamma+1.25*byl[5]*c*gamma+0.9682458365518543*byr[2]*c*gamma+0.9682458365518543*byl[2]*c*gamma-0.5590169943749475*byr[0]*c*gamma+0.5590169943749475*byl[0]*c*gamma-1.479019945774904*psr[9]*gamma+1.479019945774904*psl[9]*gamma+1.25*psr[5]*gamma+1.25*psl[5]*gamma-0.9682458365518543*psr[2]*gamma+0.9682458365518543*psl[2]*gamma+0.5590169943749475*psr[0]*gamma+0.5590169943749475*psl[0]*gamma; 
  incr[6] = (-0.75*byr[6]*c*gamma)-0.75*byl[6]*c*gamma+0.4330127018922194*byr[4]*c*gamma-0.4330127018922194*byl[4]*c*gamma+0.75*psr[6]*gamma-0.75*psl[6]*gamma-0.4330127018922194*psr[4]*gamma-0.4330127018922194*psl[4]*gamma; 
  incr[7] = (-1.25*byr[7]*c*gamma)+1.25*byl[7]*c*gamma+0.9682458365518543*byr[3]*c*gamma+0.9682458365518543*byl[3]*c*gamma-0.5590169943749476*byr[1]*c*gamma+0.5590169943749476*byl[1]*c*gamma+1.25*psr[7]*gamma+1.25*psl[7]*gamma-0.9682458365518543*psr[3]*gamma+0.9682458365518543*psl[3]*gamma+0.5590169943749476*psr[1]*gamma+0.5590169943749476*psl[1]*gamma; 
  incr[8] = (-0.25*byr[8]*c*gamma)+0.25*byl[8]*c*gamma+0.25*psr[8]*gamma+0.25*psl[8]*gamma; 
  incr[9] = (-1.75*byr[9]*c*gamma)-1.75*byl[9]*c*gamma+1.479019945774904*byr[5]*c*gamma-1.479019945774904*byl[5]*c*gamma-1.14564392373896*byr[2]*c*gamma-1.14564392373896*byl[2]*c*gamma+0.6614378277661477*byr[0]*c*gamma-0.6614378277661477*byl[0]*c*gamma+1.75*psr[9]*gamma-1.75*psl[9]*gamma-1.479019945774904*psr[5]*gamma-1.479019945774904*psl[5]*gamma+1.14564392373896*psr[2]*gamma-1.14564392373896*psl[2]*gamma-0.6614378277661477*psr[0]*gamma-0.6614378277661477*psl[0]*gamma; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 
  outByr[3] += incr[3]*dx1; 
  outByr[4] += incr[4]*dx1; 
  outByr[5] += incr[5]*dx1; 
  outByr[6] += incr[6]*dx1; 
  outByr[7] += incr[7]*dx1; 
  outByr[8] += incr[8]*dx1; 
  outByr[9] += incr[9]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += -1.0*incr[1]*dx1; 
  outByl[2] += incr[2]*dx1; 
  outByl[3] += incr[3]*dx1; 
  outByl[4] += -1.0*incr[4]*dx1; 
  outByl[5] += -1.0*incr[5]*dx1; 
  outByl[6] += incr[6]*dx1; 
  outByl[7] += -1.0*incr[7]*dx1; 
  outByl[8] += -1.0*incr[8]*dx1; 
  outByl[9] += incr[9]*dx1; 

 
  incr[0] = 0.6614378277661477*bzr[9]*c+0.6614378277661477*bzl[9]*c-0.5590169943749475*bzr[5]*c+0.5590169943749475*bzl[5]*c+0.4330127018922193*bzr[2]*c+0.4330127018922193*bzl[2]*c-0.25*bzr[0]*c+0.25*bzl[0]*c+0.6614378277661477*exr[9]-0.6614378277661477*exl[9]-0.5590169943749475*exr[5]-0.5590169943749475*exl[5]+0.4330127018922193*exr[2]-0.4330127018922193*exl[2]-0.25*exr[0]-0.25*exl[0]; 
  incr[1] = (-0.5590169943749476*bzr[7]*c)+0.5590169943749476*bzl[7]*c+0.4330127018922193*bzr[3]*c+0.4330127018922193*bzl[3]*c-0.25*bzr[1]*c+0.25*bzl[1]*c-0.5590169943749476*exr[7]-0.5590169943749476*exl[7]+0.4330127018922193*exr[3]-0.4330127018922193*exl[3]-0.25*exr[1]-0.25*exl[1]; 
  incr[2] = (-1.14564392373896*bzr[9]*c)-1.14564392373896*bzl[9]*c+0.9682458365518543*bzr[5]*c-0.9682458365518543*bzl[5]*c-0.75*bzr[2]*c-0.75*bzl[2]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c-1.14564392373896*exr[9]+1.14564392373896*exl[9]+0.9682458365518543*exr[5]+0.9682458365518543*exl[5]-0.75*exr[2]+0.75*exl[2]+0.4330127018922193*exr[0]+0.4330127018922193*exl[0]; 
  incr[3] = 0.9682458365518543*bzr[7]*c-0.9682458365518543*bzl[7]*c-0.75*bzr[3]*c-0.75*bzl[3]*c+0.4330127018922193*bzr[1]*c-0.4330127018922193*bzl[1]*c+0.9682458365518543*exr[7]+0.9682458365518543*exl[7]-0.75*exr[3]+0.75*exl[3]+0.4330127018922193*exr[1]+0.4330127018922193*exl[1]; 
  incr[4] = 0.4330127018922194*bzr[6]*c+0.4330127018922194*bzl[6]*c-0.25*bzr[4]*c+0.25*bzl[4]*c+0.4330127018922194*exr[6]-0.4330127018922194*exl[6]-0.25*exr[4]-0.25*exl[4]; 
  incr[5] = 1.479019945774904*bzr[9]*c+1.479019945774904*bzl[9]*c-1.25*bzr[5]*c+1.25*bzl[5]*c+0.9682458365518543*bzr[2]*c+0.9682458365518543*bzl[2]*c-0.5590169943749475*bzr[0]*c+0.5590169943749475*bzl[0]*c+1.479019945774904*exr[9]-1.479019945774904*exl[9]-1.25*exr[5]-1.25*exl[5]+0.9682458365518543*exr[2]-0.9682458365518543*exl[2]-0.5590169943749475*exr[0]-0.5590169943749475*exl[0]; 
  incr[6] = (-0.75*bzr[6]*c)-0.75*bzl[6]*c+0.4330127018922194*bzr[4]*c-0.4330127018922194*bzl[4]*c-0.75*exr[6]+0.75*exl[6]+0.4330127018922194*exr[4]+0.4330127018922194*exl[4]; 
  incr[7] = (-1.25*bzr[7]*c)+1.25*bzl[7]*c+0.9682458365518543*bzr[3]*c+0.9682458365518543*bzl[3]*c-0.5590169943749476*bzr[1]*c+0.5590169943749476*bzl[1]*c-1.25*exr[7]-1.25*exl[7]+0.9682458365518543*exr[3]-0.9682458365518543*exl[3]-0.5590169943749476*exr[1]-0.5590169943749476*exl[1]; 
  incr[8] = (-0.25*bzr[8]*c)+0.25*bzl[8]*c-0.25*exr[8]-0.25*exl[8]; 
  incr[9] = (-1.75*bzr[9]*c)-1.75*bzl[9]*c+1.479019945774904*bzr[5]*c-1.479019945774904*bzl[5]*c-1.14564392373896*bzr[2]*c-1.14564392373896*bzl[2]*c+0.6614378277661477*bzr[0]*c-0.6614378277661477*bzl[0]*c-1.75*exr[9]+1.75*exl[9]+1.479019945774904*exr[5]+1.479019945774904*exl[5]-1.14564392373896*exr[2]+1.14564392373896*exl[2]+0.6614378277661477*exr[0]+0.6614378277661477*exl[0]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 
  outBzr[3] += incr[3]*dx1; 
  outBzr[4] += incr[4]*dx1; 
  outBzr[5] += incr[5]*dx1; 
  outBzr[6] += incr[6]*dx1; 
  outBzr[7] += incr[7]*dx1; 
  outBzr[8] += incr[8]*dx1; 
  outBzr[9] += incr[9]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += -1.0*incr[1]*dx1; 
  outBzl[2] += incr[2]*dx1; 
  outBzl[3] += incr[3]*dx1; 
  outBzl[4] += -1.0*incr[4]*dx1; 
  outBzl[5] += -1.0*incr[5]*dx1; 
  outBzl[6] += incr[6]*dx1; 
  outBzl[7] += -1.0*incr[7]*dx1; 
  outBzl[8] += -1.0*incr[8]*dx1; 
  outBzl[9] += incr[9]*dx1; 

 
  incr[0] = 0.6614378277661477*phr[9]*c*chi+0.6614378277661477*phl[9]*c*chi-0.5590169943749475*phr[5]*c*chi+0.5590169943749475*phl[5]*c*chi+0.4330127018922193*phr[2]*c*chi+0.4330127018922193*phl[2]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi-0.6614378277661477*eyr[9]*chi+0.6614378277661477*eyl[9]*chi+0.5590169943749475*eyr[5]*chi+0.5590169943749475*eyl[5]*chi-0.4330127018922193*eyr[2]*chi+0.4330127018922193*eyl[2]*chi+0.25*eyr[0]*chi+0.25*eyl[0]*chi; 
  incr[1] = (-0.5590169943749476*phr[7]*c*chi)+0.5590169943749476*phl[7]*c*chi+0.4330127018922193*phr[3]*c*chi+0.4330127018922193*phl[3]*c*chi-0.25*phr[1]*c*chi+0.25*phl[1]*c*chi+0.5590169943749476*eyr[7]*chi+0.5590169943749476*eyl[7]*chi-0.4330127018922193*eyr[3]*chi+0.4330127018922193*eyl[3]*chi+0.25*eyr[1]*chi+0.25*eyl[1]*chi; 
  incr[2] = (-1.14564392373896*phr[9]*c*chi)-1.14564392373896*phl[9]*c*chi+0.9682458365518543*phr[5]*c*chi-0.9682458365518543*phl[5]*c*chi-0.75*phr[2]*c*chi-0.75*phl[2]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi+1.14564392373896*eyr[9]*chi-1.14564392373896*eyl[9]*chi-0.9682458365518543*eyr[5]*chi-0.9682458365518543*eyl[5]*chi+0.75*eyr[2]*chi-0.75*eyl[2]*chi-0.4330127018922193*eyr[0]*chi-0.4330127018922193*eyl[0]*chi; 
  incr[3] = 0.9682458365518543*phr[7]*c*chi-0.9682458365518543*phl[7]*c*chi-0.75*phr[3]*c*chi-0.75*phl[3]*c*chi+0.4330127018922193*phr[1]*c*chi-0.4330127018922193*phl[1]*c*chi-0.9682458365518543*eyr[7]*chi-0.9682458365518543*eyl[7]*chi+0.75*eyr[3]*chi-0.75*eyl[3]*chi-0.4330127018922193*eyr[1]*chi-0.4330127018922193*eyl[1]*chi; 
  incr[4] = 0.4330127018922194*phr[6]*c*chi+0.4330127018922194*phl[6]*c*chi-0.25*phr[4]*c*chi+0.25*phl[4]*c*chi-0.4330127018922194*eyr[6]*chi+0.4330127018922194*eyl[6]*chi+0.25*eyr[4]*chi+0.25*eyl[4]*chi; 
  incr[5] = 1.479019945774904*phr[9]*c*chi+1.479019945774904*phl[9]*c*chi-1.25*phr[5]*c*chi+1.25*phl[5]*c*chi+0.9682458365518543*phr[2]*c*chi+0.9682458365518543*phl[2]*c*chi-0.5590169943749475*phr[0]*c*chi+0.5590169943749475*phl[0]*c*chi-1.479019945774904*eyr[9]*chi+1.479019945774904*eyl[9]*chi+1.25*eyr[5]*chi+1.25*eyl[5]*chi-0.9682458365518543*eyr[2]*chi+0.9682458365518543*eyl[2]*chi+0.5590169943749475*eyr[0]*chi+0.5590169943749475*eyl[0]*chi; 
  incr[6] = (-0.75*phr[6]*c*chi)-0.75*phl[6]*c*chi+0.4330127018922194*phr[4]*c*chi-0.4330127018922194*phl[4]*c*chi+0.75*eyr[6]*chi-0.75*eyl[6]*chi-0.4330127018922194*eyr[4]*chi-0.4330127018922194*eyl[4]*chi; 
  incr[7] = (-1.25*phr[7]*c*chi)+1.25*phl[7]*c*chi+0.9682458365518543*phr[3]*c*chi+0.9682458365518543*phl[3]*c*chi-0.5590169943749476*phr[1]*c*chi+0.5590169943749476*phl[1]*c*chi+1.25*eyr[7]*chi+1.25*eyl[7]*chi-0.9682458365518543*eyr[3]*chi+0.9682458365518543*eyl[3]*chi+0.5590169943749476*eyr[1]*chi+0.5590169943749476*eyl[1]*chi; 
  incr[8] = (-0.25*phr[8]*c*chi)+0.25*phl[8]*c*chi+0.25*eyr[8]*chi+0.25*eyl[8]*chi; 
  incr[9] = (-1.75*phr[9]*c*chi)-1.75*phl[9]*c*chi+1.479019945774904*phr[5]*c*chi-1.479019945774904*phl[5]*c*chi-1.14564392373896*phr[2]*c*chi-1.14564392373896*phl[2]*c*chi+0.6614378277661477*phr[0]*c*chi-0.6614378277661477*phl[0]*c*chi+1.75*eyr[9]*chi-1.75*eyl[9]*chi-1.479019945774904*eyr[5]*chi-1.479019945774904*eyl[5]*chi+1.14564392373896*eyr[2]*chi-1.14564392373896*eyl[2]*chi-0.6614378277661477*eyr[0]*chi-0.6614378277661477*eyl[0]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 
  outPhr[3] += incr[3]*dx1; 
  outPhr[4] += incr[4]*dx1; 
  outPhr[5] += incr[5]*dx1; 
  outPhr[6] += incr[6]*dx1; 
  outPhr[7] += incr[7]*dx1; 
  outPhr[8] += incr[8]*dx1; 
  outPhr[9] += incr[9]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += -1.0*incr[1]*dx1; 
  outPhl[2] += incr[2]*dx1; 
  outPhl[3] += incr[3]*dx1; 
  outPhl[4] += -1.0*incr[4]*dx1; 
  outPhl[5] += -1.0*incr[5]*dx1; 
  outPhl[6] += incr[6]*dx1; 
  outPhl[7] += -1.0*incr[7]*dx1; 
  outPhl[8] += -1.0*incr[8]*dx1; 
  outPhl[9] += incr[9]*dx1; 

 
  incr[0] = (-0.6614378277661477*byr[9]*c2*gamma)+0.6614378277661477*byl[9]*c2*gamma+0.5590169943749475*byr[5]*c2*gamma+0.5590169943749475*byl[5]*c2*gamma-0.4330127018922193*byr[2]*c2*gamma+0.4330127018922193*byl[2]*c2*gamma+0.25*byr[0]*c2*gamma+0.25*byl[0]*c2*gamma+0.6614378277661477*psr[9]*c*gamma+0.6614378277661477*psl[9]*c*gamma-0.5590169943749475*psr[5]*c*gamma+0.5590169943749475*psl[5]*c*gamma+0.4330127018922193*psr[2]*c*gamma+0.4330127018922193*psl[2]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = 0.5590169943749476*byr[7]*c2*gamma+0.5590169943749476*byl[7]*c2*gamma-0.4330127018922193*byr[3]*c2*gamma+0.4330127018922193*byl[3]*c2*gamma+0.25*byr[1]*c2*gamma+0.25*byl[1]*c2*gamma-0.5590169943749476*psr[7]*c*gamma+0.5590169943749476*psl[7]*c*gamma+0.4330127018922193*psr[3]*c*gamma+0.4330127018922193*psl[3]*c*gamma-0.25*psr[1]*c*gamma+0.25*psl[1]*c*gamma; 
  incr[2] = 1.14564392373896*byr[9]*c2*gamma-1.14564392373896*byl[9]*c2*gamma-0.9682458365518543*byr[5]*c2*gamma-0.9682458365518543*byl[5]*c2*gamma+0.75*byr[2]*c2*gamma-0.75*byl[2]*c2*gamma-0.4330127018922193*byr[0]*c2*gamma-0.4330127018922193*byl[0]*c2*gamma-1.14564392373896*psr[9]*c*gamma-1.14564392373896*psl[9]*c*gamma+0.9682458365518543*psr[5]*c*gamma-0.9682458365518543*psl[5]*c*gamma-0.75*psr[2]*c*gamma-0.75*psl[2]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 
  incr[3] = (-0.9682458365518543*byr[7]*c2*gamma)-0.9682458365518543*byl[7]*c2*gamma+0.75*byr[3]*c2*gamma-0.75*byl[3]*c2*gamma-0.4330127018922193*byr[1]*c2*gamma-0.4330127018922193*byl[1]*c2*gamma+0.9682458365518543*psr[7]*c*gamma-0.9682458365518543*psl[7]*c*gamma-0.75*psr[3]*c*gamma-0.75*psl[3]*c*gamma+0.4330127018922193*psr[1]*c*gamma-0.4330127018922193*psl[1]*c*gamma; 
  incr[4] = (-0.4330127018922194*byr[6]*c2*gamma)+0.4330127018922194*byl[6]*c2*gamma+0.25*byr[4]*c2*gamma+0.25*byl[4]*c2*gamma+0.4330127018922194*psr[6]*c*gamma+0.4330127018922194*psl[6]*c*gamma-0.25*psr[4]*c*gamma+0.25*psl[4]*c*gamma; 
  incr[5] = (-1.479019945774904*byr[9]*c2*gamma)+1.479019945774904*byl[9]*c2*gamma+1.25*byr[5]*c2*gamma+1.25*byl[5]*c2*gamma-0.9682458365518543*byr[2]*c2*gamma+0.9682458365518543*byl[2]*c2*gamma+0.5590169943749475*byr[0]*c2*gamma+0.5590169943749475*byl[0]*c2*gamma+1.479019945774904*psr[9]*c*gamma+1.479019945774904*psl[9]*c*gamma-1.25*psr[5]*c*gamma+1.25*psl[5]*c*gamma+0.9682458365518543*psr[2]*c*gamma+0.9682458365518543*psl[2]*c*gamma-0.5590169943749475*psr[0]*c*gamma+0.5590169943749475*psl[0]*c*gamma; 
  incr[6] = 0.75*byr[6]*c2*gamma-0.75*byl[6]*c2*gamma-0.4330127018922194*byr[4]*c2*gamma-0.4330127018922194*byl[4]*c2*gamma-0.75*psr[6]*c*gamma-0.75*psl[6]*c*gamma+0.4330127018922194*psr[4]*c*gamma-0.4330127018922194*psl[4]*c*gamma; 
  incr[7] = 1.25*byr[7]*c2*gamma+1.25*byl[7]*c2*gamma-0.9682458365518543*byr[3]*c2*gamma+0.9682458365518543*byl[3]*c2*gamma+0.5590169943749476*byr[1]*c2*gamma+0.5590169943749476*byl[1]*c2*gamma-1.25*psr[7]*c*gamma+1.25*psl[7]*c*gamma+0.9682458365518543*psr[3]*c*gamma+0.9682458365518543*psl[3]*c*gamma-0.5590169943749476*psr[1]*c*gamma+0.5590169943749476*psl[1]*c*gamma; 
  incr[8] = 0.25*byr[8]*c2*gamma+0.25*byl[8]*c2*gamma-0.25*psr[8]*c*gamma+0.25*psl[8]*c*gamma; 
  incr[9] = 1.75*byr[9]*c2*gamma-1.75*byl[9]*c2*gamma-1.479019945774904*byr[5]*c2*gamma-1.479019945774904*byl[5]*c2*gamma+1.14564392373896*byr[2]*c2*gamma-1.14564392373896*byl[2]*c2*gamma-0.6614378277661477*byr[0]*c2*gamma-0.6614378277661477*byl[0]*c2*gamma-1.75*psr[9]*c*gamma-1.75*psl[9]*c*gamma+1.479019945774904*psr[5]*c*gamma-1.479019945774904*psl[5]*c*gamma-1.14564392373896*psr[2]*c*gamma-1.14564392373896*psl[2]*c*gamma+0.6614378277661477*psr[0]*c*gamma-0.6614378277661477*psl[0]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 
  outPsr[3] += incr[3]*dx1; 
  outPsr[4] += incr[4]*dx1; 
  outPsr[5] += incr[5]*dx1; 
  outPsr[6] += incr[6]*dx1; 
  outPsr[7] += incr[7]*dx1; 
  outPsr[8] += incr[8]*dx1; 
  outPsr[9] += incr[9]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += -1.0*incr[1]*dx1; 
  outPsl[2] += incr[2]*dx1; 
  outPsl[3] += incr[3]*dx1; 
  outPsl[4] += -1.0*incr[4]*dx1; 
  outPsl[5] += -1.0*incr[5]*dx1; 
  outPsl[6] += incr[6]*dx1; 
  outPsl[7] += -1.0*incr[7]*dx1; 
  outPsl[8] += -1.0*incr[8]*dx1; 
  outPsl[9] += incr[9]*dx1; 

 
} 
void MaxwellSurf2xMax_Y_P4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double dx1 = 2.0/dx[1]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[15]; 
  const double *ezl = &ql[30]; 
  const double *bxl = &ql[45]; 
  const double *byl = &ql[60]; 
  const double *bzl = &ql[75]; 
  const double *phl = &ql[90]; 
  const double *psl = &ql[105]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[15]; 
  double *outEzl = &outl[30]; 
  double *outBxl = &outl[45]; 
  double *outByl = &outl[60]; 
  double *outBzl = &outl[75]; 
  double *outPhl = &outl[90]; 
  double *outPsl = &outl[105]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[15]; 
  const double *ezr = &qr[30]; 
  const double *bxr = &qr[45]; 
  const double *byr = &qr[60]; 
  const double *bzr = &qr[75]; 
  const double *phr = &qr[90]; 
  const double *psr = &qr[105]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[15]; 
  double *outEzr = &outr[30]; 
  double *outBxr = &outr[45]; 
  double *outByr = &outr[60]; 
  double *outBzr = &outr[75]; 
  double *outPhr = &outr[90]; 
  double *outPsr = &outr[105]; 
 
  double incr[15]; 
 
  incr[0] = (-0.75*bzr[14]*c2)-0.75*bzl[14]*c2+0.6614378277661477*bzr[9]*c2-0.6614378277661477*bzl[9]*c2-0.5590169943749475*bzr[5]*c2-0.5590169943749475*bzl[5]*c2+0.4330127018922193*bzr[2]*c2-0.4330127018922193*bzl[2]*c2-0.25*bzr[0]*c2-0.25*bzl[0]*c2-0.75*exr[14]*c+0.75*exl[14]*c+0.6614378277661477*exr[9]*c+0.6614378277661477*exl[9]*c-0.5590169943749475*exr[5]*c+0.5590169943749475*exl[5]*c+0.4330127018922193*exr[2]*c+0.4330127018922193*exl[2]*c-0.25*exr[0]*c+0.25*exl[0]*c; 
  incr[1] = 0.6614378277661477*bzr[12]*c2-0.6614378277661477*bzl[12]*c2-0.5590169943749476*bzr[7]*c2-0.5590169943749476*bzl[7]*c2+0.4330127018922193*bzr[3]*c2-0.4330127018922193*bzl[3]*c2-0.25*bzr[1]*c2-0.25*bzl[1]*c2+0.6614378277661477*exr[12]*c+0.6614378277661477*exl[12]*c-0.5590169943749476*exr[7]*c+0.5590169943749476*exl[7]*c+0.4330127018922193*exr[3]*c+0.4330127018922193*exl[3]*c-0.25*exr[1]*c+0.25*exl[1]*c; 
  incr[2] = 1.299038105676658*bzr[14]*c2+1.299038105676658*bzl[14]*c2-1.14564392373896*bzr[9]*c2+1.14564392373896*bzl[9]*c2+0.9682458365518543*bzr[5]*c2+0.9682458365518543*bzl[5]*c2-0.75*bzr[2]*c2+0.75*bzl[2]*c2+0.4330127018922193*bzr[0]*c2+0.4330127018922193*bzl[0]*c2+1.299038105676658*exr[14]*c-1.299038105676658*exl[14]*c-1.14564392373896*exr[9]*c-1.14564392373896*exl[9]*c+0.9682458365518543*exr[5]*c-0.9682458365518543*exl[5]*c-0.75*exr[2]*c-0.75*exl[2]*c+0.4330127018922193*exr[0]*c-0.4330127018922193*exl[0]*c; 
  incr[3] = (-1.14564392373896*bzr[12]*c2)+1.14564392373896*bzl[12]*c2+0.9682458365518543*bzr[7]*c2+0.9682458365518543*bzl[7]*c2-0.75*bzr[3]*c2+0.75*bzl[3]*c2+0.4330127018922193*bzr[1]*c2+0.4330127018922193*bzl[1]*c2-1.14564392373896*exr[12]*c-1.14564392373896*exl[12]*c+0.9682458365518543*exr[7]*c-0.9682458365518543*exl[7]*c-0.75*exr[3]*c-0.75*exl[3]*c+0.4330127018922193*exr[1]*c-0.4330127018922193*exl[1]*c; 
  incr[4] = (-0.5590169943749475*bzr[10]*c2)-0.5590169943749475*bzl[10]*c2+0.4330127018922194*bzr[6]*c2-0.4330127018922194*bzl[6]*c2-0.25*bzr[4]*c2-0.25*bzl[4]*c2-0.5590169943749475*exr[10]*c+0.5590169943749475*exl[10]*c+0.4330127018922194*exr[6]*c+0.4330127018922194*exl[6]*c-0.25*exr[4]*c+0.25*exl[4]*c; 
  incr[5] = (-1.677050983124842*bzr[14]*c2)-1.677050983124842*bzl[14]*c2+1.479019945774904*bzr[9]*c2-1.479019945774904*bzl[9]*c2-1.25*bzr[5]*c2-1.25*bzl[5]*c2+0.9682458365518543*bzr[2]*c2-0.9682458365518543*bzl[2]*c2-0.5590169943749475*bzr[0]*c2-0.5590169943749475*bzl[0]*c2-1.677050983124842*exr[14]*c+1.677050983124842*exl[14]*c+1.479019945774904*exr[9]*c+1.479019945774904*exl[9]*c-1.25*exr[5]*c+1.25*exl[5]*c+0.9682458365518543*exr[2]*c+0.9682458365518543*exl[2]*c-0.5590169943749475*exr[0]*c+0.5590169943749475*exl[0]*c; 
  incr[6] = 0.9682458365518543*bzr[10]*c2+0.9682458365518543*bzl[10]*c2-0.75*bzr[6]*c2+0.75*bzl[6]*c2+0.4330127018922194*bzr[4]*c2+0.4330127018922194*bzl[4]*c2+0.9682458365518543*exr[10]*c-0.9682458365518543*exl[10]*c-0.75*exr[6]*c-0.75*exl[6]*c+0.4330127018922194*exr[4]*c-0.4330127018922194*exl[4]*c; 
  incr[7] = 1.479019945774904*bzr[12]*c2-1.479019945774904*bzl[12]*c2-1.25*bzr[7]*c2-1.25*bzl[7]*c2+0.9682458365518543*bzr[3]*c2-0.9682458365518543*bzl[3]*c2-0.5590169943749476*bzr[1]*c2-0.5590169943749476*bzl[1]*c2+1.479019945774904*exr[12]*c+1.479019945774904*exl[12]*c-1.25*exr[7]*c+1.25*exl[7]*c+0.9682458365518543*exr[3]*c+0.9682458365518543*exl[3]*c-0.5590169943749476*exr[1]*c+0.5590169943749476*exl[1]*c; 
  incr[8] = 0.4330127018922193*bzr[11]*c2-0.4330127018922193*bzl[11]*c2-0.25*bzr[8]*c2-0.25*bzl[8]*c2+0.4330127018922193*exr[11]*c+0.4330127018922193*exl[11]*c-0.25*exr[8]*c+0.25*exl[8]*c; 
  incr[9] = 1.984313483298443*bzr[14]*c2+1.984313483298443*bzl[14]*c2-1.75*bzr[9]*c2+1.75*bzl[9]*c2+1.479019945774904*bzr[5]*c2+1.479019945774904*bzl[5]*c2-1.14564392373896*bzr[2]*c2+1.14564392373896*bzl[2]*c2+0.6614378277661477*bzr[0]*c2+0.6614378277661477*bzl[0]*c2+1.984313483298443*exr[14]*c-1.984313483298443*exl[14]*c-1.75*exr[9]*c-1.75*exl[9]*c+1.479019945774904*exr[5]*c-1.479019945774904*exl[5]*c-1.14564392373896*exr[2]*c-1.14564392373896*exl[2]*c+0.6614378277661477*exr[0]*c-0.6614378277661477*exl[0]*c; 
  incr[10] = (-1.25*bzr[10]*c2)-1.25*bzl[10]*c2+0.9682458365518543*bzr[6]*c2-0.9682458365518543*bzl[6]*c2-0.5590169943749475*bzr[4]*c2-0.5590169943749475*bzl[4]*c2-1.25*exr[10]*c+1.25*exl[10]*c+0.9682458365518543*exr[6]*c+0.9682458365518543*exl[6]*c-0.5590169943749475*exr[4]*c+0.5590169943749475*exl[4]*c; 
  incr[11] = (-0.75*bzr[11]*c2)+0.75*bzl[11]*c2+0.4330127018922193*bzr[8]*c2+0.4330127018922193*bzl[8]*c2-0.75*exr[11]*c-0.75*exl[11]*c+0.4330127018922193*exr[8]*c-0.4330127018922193*exl[8]*c; 
  incr[12] = (-1.75*bzr[12]*c2)+1.75*bzl[12]*c2+1.479019945774904*bzr[7]*c2+1.479019945774904*bzl[7]*c2-1.14564392373896*bzr[3]*c2+1.14564392373896*bzl[3]*c2+0.6614378277661477*bzr[1]*c2+0.6614378277661477*bzl[1]*c2-1.75*exr[12]*c-1.75*exl[12]*c+1.479019945774904*exr[7]*c-1.479019945774904*exl[7]*c-1.14564392373896*exr[3]*c-1.14564392373896*exl[3]*c+0.6614378277661477*exr[1]*c-0.6614378277661477*exl[1]*c; 
  incr[13] = (-0.25*bzr[13]*c2)-0.25*bzl[13]*c2-0.25*exr[13]*c+0.25*exl[13]*c; 
  incr[14] = (-2.25*bzr[14]*c2)-2.25*bzl[14]*c2+1.984313483298443*bzr[9]*c2-1.984313483298443*bzl[9]*c2-1.677050983124842*bzr[5]*c2-1.677050983124842*bzl[5]*c2+1.299038105676658*bzr[2]*c2-1.299038105676658*bzl[2]*c2-0.75*bzr[0]*c2-0.75*bzl[0]*c2-2.25*exr[14]*c+2.25*exl[14]*c+1.984313483298443*exr[9]*c+1.984313483298443*exl[9]*c-1.677050983124842*exr[5]*c+1.677050983124842*exl[5]*c+1.299038105676658*exr[2]*c+1.299038105676658*exl[2]*c-0.75*exr[0]*c+0.75*exl[0]*c; 

  outExr[0] += incr[0]*dx1; 
  outExr[1] += incr[1]*dx1; 
  outExr[2] += incr[2]*dx1; 
  outExr[3] += incr[3]*dx1; 
  outExr[4] += incr[4]*dx1; 
  outExr[5] += incr[5]*dx1; 
  outExr[6] += incr[6]*dx1; 
  outExr[7] += incr[7]*dx1; 
  outExr[8] += incr[8]*dx1; 
  outExr[9] += incr[9]*dx1; 
  outExr[10] += incr[10]*dx1; 
  outExr[11] += incr[11]*dx1; 
  outExr[12] += incr[12]*dx1; 
  outExr[13] += incr[13]*dx1; 
  outExr[14] += incr[14]*dx1; 

  outExl[0] += -1.0*incr[0]*dx1; 
  outExl[1] += -1.0*incr[1]*dx1; 
  outExl[2] += incr[2]*dx1; 
  outExl[3] += incr[3]*dx1; 
  outExl[4] += -1.0*incr[4]*dx1; 
  outExl[5] += -1.0*incr[5]*dx1; 
  outExl[6] += incr[6]*dx1; 
  outExl[7] += -1.0*incr[7]*dx1; 
  outExl[8] += -1.0*incr[8]*dx1; 
  outExl[9] += incr[9]*dx1; 
  outExl[10] += -1.0*incr[10]*dx1; 
  outExl[11] += incr[11]*dx1; 
  outExl[12] += incr[12]*dx1; 
  outExl[13] += -1.0*incr[13]*dx1; 
  outExl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = 0.75*phr[14]*c2*chi+0.75*phl[14]*c2*chi-0.6614378277661477*phr[9]*c2*chi+0.6614378277661477*phl[9]*c2*chi+0.5590169943749475*phr[5]*c2*chi+0.5590169943749475*phl[5]*c2*chi-0.4330127018922193*phr[2]*c2*chi+0.4330127018922193*phl[2]*c2*chi+0.25*phr[0]*c2*chi+0.25*phl[0]*c2*chi-0.75*eyr[14]*c*chi+0.75*eyl[14]*c*chi+0.6614378277661477*eyr[9]*c*chi+0.6614378277661477*eyl[9]*c*chi-0.5590169943749475*eyr[5]*c*chi+0.5590169943749475*eyl[5]*c*chi+0.4330127018922193*eyr[2]*c*chi+0.4330127018922193*eyl[2]*c*chi-0.25*eyr[0]*c*chi+0.25*eyl[0]*c*chi; 
  incr[1] = (-0.6614378277661477*phr[12]*c2*chi)+0.6614378277661477*phl[12]*c2*chi+0.5590169943749476*phr[7]*c2*chi+0.5590169943749476*phl[7]*c2*chi-0.4330127018922193*phr[3]*c2*chi+0.4330127018922193*phl[3]*c2*chi+0.25*phr[1]*c2*chi+0.25*phl[1]*c2*chi+0.6614378277661477*eyr[12]*c*chi+0.6614378277661477*eyl[12]*c*chi-0.5590169943749476*eyr[7]*c*chi+0.5590169943749476*eyl[7]*c*chi+0.4330127018922193*eyr[3]*c*chi+0.4330127018922193*eyl[3]*c*chi-0.25*eyr[1]*c*chi+0.25*eyl[1]*c*chi; 
  incr[2] = (-1.299038105676658*phr[14]*c2*chi)-1.299038105676658*phl[14]*c2*chi+1.14564392373896*phr[9]*c2*chi-1.14564392373896*phl[9]*c2*chi-0.9682458365518543*phr[5]*c2*chi-0.9682458365518543*phl[5]*c2*chi+0.75*phr[2]*c2*chi-0.75*phl[2]*c2*chi-0.4330127018922193*phr[0]*c2*chi-0.4330127018922193*phl[0]*c2*chi+1.299038105676658*eyr[14]*c*chi-1.299038105676658*eyl[14]*c*chi-1.14564392373896*eyr[9]*c*chi-1.14564392373896*eyl[9]*c*chi+0.9682458365518543*eyr[5]*c*chi-0.9682458365518543*eyl[5]*c*chi-0.75*eyr[2]*c*chi-0.75*eyl[2]*c*chi+0.4330127018922193*eyr[0]*c*chi-0.4330127018922193*eyl[0]*c*chi; 
  incr[3] = 1.14564392373896*phr[12]*c2*chi-1.14564392373896*phl[12]*c2*chi-0.9682458365518543*phr[7]*c2*chi-0.9682458365518543*phl[7]*c2*chi+0.75*phr[3]*c2*chi-0.75*phl[3]*c2*chi-0.4330127018922193*phr[1]*c2*chi-0.4330127018922193*phl[1]*c2*chi-1.14564392373896*eyr[12]*c*chi-1.14564392373896*eyl[12]*c*chi+0.9682458365518543*eyr[7]*c*chi-0.9682458365518543*eyl[7]*c*chi-0.75*eyr[3]*c*chi-0.75*eyl[3]*c*chi+0.4330127018922193*eyr[1]*c*chi-0.4330127018922193*eyl[1]*c*chi; 
  incr[4] = 0.5590169943749475*phr[10]*c2*chi+0.5590169943749475*phl[10]*c2*chi-0.4330127018922194*phr[6]*c2*chi+0.4330127018922194*phl[6]*c2*chi+0.25*phr[4]*c2*chi+0.25*phl[4]*c2*chi-0.5590169943749475*eyr[10]*c*chi+0.5590169943749475*eyl[10]*c*chi+0.4330127018922194*eyr[6]*c*chi+0.4330127018922194*eyl[6]*c*chi-0.25*eyr[4]*c*chi+0.25*eyl[4]*c*chi; 
  incr[5] = 1.677050983124842*phr[14]*c2*chi+1.677050983124842*phl[14]*c2*chi-1.479019945774904*phr[9]*c2*chi+1.479019945774904*phl[9]*c2*chi+1.25*phr[5]*c2*chi+1.25*phl[5]*c2*chi-0.9682458365518543*phr[2]*c2*chi+0.9682458365518543*phl[2]*c2*chi+0.5590169943749475*phr[0]*c2*chi+0.5590169943749475*phl[0]*c2*chi-1.677050983124842*eyr[14]*c*chi+1.677050983124842*eyl[14]*c*chi+1.479019945774904*eyr[9]*c*chi+1.479019945774904*eyl[9]*c*chi-1.25*eyr[5]*c*chi+1.25*eyl[5]*c*chi+0.9682458365518543*eyr[2]*c*chi+0.9682458365518543*eyl[2]*c*chi-0.5590169943749475*eyr[0]*c*chi+0.5590169943749475*eyl[0]*c*chi; 
  incr[6] = (-0.9682458365518543*phr[10]*c2*chi)-0.9682458365518543*phl[10]*c2*chi+0.75*phr[6]*c2*chi-0.75*phl[6]*c2*chi-0.4330127018922194*phr[4]*c2*chi-0.4330127018922194*phl[4]*c2*chi+0.9682458365518543*eyr[10]*c*chi-0.9682458365518543*eyl[10]*c*chi-0.75*eyr[6]*c*chi-0.75*eyl[6]*c*chi+0.4330127018922194*eyr[4]*c*chi-0.4330127018922194*eyl[4]*c*chi; 
  incr[7] = (-1.479019945774904*phr[12]*c2*chi)+1.479019945774904*phl[12]*c2*chi+1.25*phr[7]*c2*chi+1.25*phl[7]*c2*chi-0.9682458365518543*phr[3]*c2*chi+0.9682458365518543*phl[3]*c2*chi+0.5590169943749476*phr[1]*c2*chi+0.5590169943749476*phl[1]*c2*chi+1.479019945774904*eyr[12]*c*chi+1.479019945774904*eyl[12]*c*chi-1.25*eyr[7]*c*chi+1.25*eyl[7]*c*chi+0.9682458365518543*eyr[3]*c*chi+0.9682458365518543*eyl[3]*c*chi-0.5590169943749476*eyr[1]*c*chi+0.5590169943749476*eyl[1]*c*chi; 
  incr[8] = (-0.4330127018922193*phr[11]*c2*chi)+0.4330127018922193*phl[11]*c2*chi+0.25*phr[8]*c2*chi+0.25*phl[8]*c2*chi+0.4330127018922193*eyr[11]*c*chi+0.4330127018922193*eyl[11]*c*chi-0.25*eyr[8]*c*chi+0.25*eyl[8]*c*chi; 
  incr[9] = (-1.984313483298443*phr[14]*c2*chi)-1.984313483298443*phl[14]*c2*chi+1.75*phr[9]*c2*chi-1.75*phl[9]*c2*chi-1.479019945774904*phr[5]*c2*chi-1.479019945774904*phl[5]*c2*chi+1.14564392373896*phr[2]*c2*chi-1.14564392373896*phl[2]*c2*chi-0.6614378277661477*phr[0]*c2*chi-0.6614378277661477*phl[0]*c2*chi+1.984313483298443*eyr[14]*c*chi-1.984313483298443*eyl[14]*c*chi-1.75*eyr[9]*c*chi-1.75*eyl[9]*c*chi+1.479019945774904*eyr[5]*c*chi-1.479019945774904*eyl[5]*c*chi-1.14564392373896*eyr[2]*c*chi-1.14564392373896*eyl[2]*c*chi+0.6614378277661477*eyr[0]*c*chi-0.6614378277661477*eyl[0]*c*chi; 
  incr[10] = 1.25*phr[10]*c2*chi+1.25*phl[10]*c2*chi-0.9682458365518543*phr[6]*c2*chi+0.9682458365518543*phl[6]*c2*chi+0.5590169943749475*phr[4]*c2*chi+0.5590169943749475*phl[4]*c2*chi-1.25*eyr[10]*c*chi+1.25*eyl[10]*c*chi+0.9682458365518543*eyr[6]*c*chi+0.9682458365518543*eyl[6]*c*chi-0.5590169943749475*eyr[4]*c*chi+0.5590169943749475*eyl[4]*c*chi; 
  incr[11] = 0.75*phr[11]*c2*chi-0.75*phl[11]*c2*chi-0.4330127018922193*phr[8]*c2*chi-0.4330127018922193*phl[8]*c2*chi-0.75*eyr[11]*c*chi-0.75*eyl[11]*c*chi+0.4330127018922193*eyr[8]*c*chi-0.4330127018922193*eyl[8]*c*chi; 
  incr[12] = 1.75*phr[12]*c2*chi-1.75*phl[12]*c2*chi-1.479019945774904*phr[7]*c2*chi-1.479019945774904*phl[7]*c2*chi+1.14564392373896*phr[3]*c2*chi-1.14564392373896*phl[3]*c2*chi-0.6614378277661477*phr[1]*c2*chi-0.6614378277661477*phl[1]*c2*chi-1.75*eyr[12]*c*chi-1.75*eyl[12]*c*chi+1.479019945774904*eyr[7]*c*chi-1.479019945774904*eyl[7]*c*chi-1.14564392373896*eyr[3]*c*chi-1.14564392373896*eyl[3]*c*chi+0.6614378277661477*eyr[1]*c*chi-0.6614378277661477*eyl[1]*c*chi; 
  incr[13] = 0.25*phr[13]*c2*chi+0.25*phl[13]*c2*chi-0.25*eyr[13]*c*chi+0.25*eyl[13]*c*chi; 
  incr[14] = 2.25*phr[14]*c2*chi+2.25*phl[14]*c2*chi-1.984313483298443*phr[9]*c2*chi+1.984313483298443*phl[9]*c2*chi+1.677050983124842*phr[5]*c2*chi+1.677050983124842*phl[5]*c2*chi-1.299038105676658*phr[2]*c2*chi+1.299038105676658*phl[2]*c2*chi+0.75*phr[0]*c2*chi+0.75*phl[0]*c2*chi-2.25*eyr[14]*c*chi+2.25*eyl[14]*c*chi+1.984313483298443*eyr[9]*c*chi+1.984313483298443*eyl[9]*c*chi-1.677050983124842*eyr[5]*c*chi+1.677050983124842*eyl[5]*c*chi+1.299038105676658*eyr[2]*c*chi+1.299038105676658*eyl[2]*c*chi-0.75*eyr[0]*c*chi+0.75*eyl[0]*c*chi; 

  outEyr[0] += incr[0]*dx1; 
  outEyr[1] += incr[1]*dx1; 
  outEyr[2] += incr[2]*dx1; 
  outEyr[3] += incr[3]*dx1; 
  outEyr[4] += incr[4]*dx1; 
  outEyr[5] += incr[5]*dx1; 
  outEyr[6] += incr[6]*dx1; 
  outEyr[7] += incr[7]*dx1; 
  outEyr[8] += incr[8]*dx1; 
  outEyr[9] += incr[9]*dx1; 
  outEyr[10] += incr[10]*dx1; 
  outEyr[11] += incr[11]*dx1; 
  outEyr[12] += incr[12]*dx1; 
  outEyr[13] += incr[13]*dx1; 
  outEyr[14] += incr[14]*dx1; 

  outEyl[0] += -1.0*incr[0]*dx1; 
  outEyl[1] += -1.0*incr[1]*dx1; 
  outEyl[2] += incr[2]*dx1; 
  outEyl[3] += incr[3]*dx1; 
  outEyl[4] += -1.0*incr[4]*dx1; 
  outEyl[5] += -1.0*incr[5]*dx1; 
  outEyl[6] += incr[6]*dx1; 
  outEyl[7] += -1.0*incr[7]*dx1; 
  outEyl[8] += -1.0*incr[8]*dx1; 
  outEyl[9] += incr[9]*dx1; 
  outEyl[10] += -1.0*incr[10]*dx1; 
  outEyl[11] += incr[11]*dx1; 
  outEyl[12] += incr[12]*dx1; 
  outEyl[13] += -1.0*incr[13]*dx1; 
  outEyl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = 0.75*bxr[14]*c2+0.75*bxl[14]*c2-0.6614378277661477*bxr[9]*c2+0.6614378277661477*bxl[9]*c2+0.5590169943749475*bxr[5]*c2+0.5590169943749475*bxl[5]*c2-0.4330127018922193*bxr[2]*c2+0.4330127018922193*bxl[2]*c2+0.25*bxr[0]*c2+0.25*bxl[0]*c2-0.75*ezr[14]*c+0.75*ezl[14]*c+0.6614378277661477*ezr[9]*c+0.6614378277661477*ezl[9]*c-0.5590169943749475*ezr[5]*c+0.5590169943749475*ezl[5]*c+0.4330127018922193*ezr[2]*c+0.4330127018922193*ezl[2]*c-0.25*ezr[0]*c+0.25*ezl[0]*c; 
  incr[1] = (-0.6614378277661477*bxr[12]*c2)+0.6614378277661477*bxl[12]*c2+0.5590169943749476*bxr[7]*c2+0.5590169943749476*bxl[7]*c2-0.4330127018922193*bxr[3]*c2+0.4330127018922193*bxl[3]*c2+0.25*bxr[1]*c2+0.25*bxl[1]*c2+0.6614378277661477*ezr[12]*c+0.6614378277661477*ezl[12]*c-0.5590169943749476*ezr[7]*c+0.5590169943749476*ezl[7]*c+0.4330127018922193*ezr[3]*c+0.4330127018922193*ezl[3]*c-0.25*ezr[1]*c+0.25*ezl[1]*c; 
  incr[2] = (-1.299038105676658*bxr[14]*c2)-1.299038105676658*bxl[14]*c2+1.14564392373896*bxr[9]*c2-1.14564392373896*bxl[9]*c2-0.9682458365518543*bxr[5]*c2-0.9682458365518543*bxl[5]*c2+0.75*bxr[2]*c2-0.75*bxl[2]*c2-0.4330127018922193*bxr[0]*c2-0.4330127018922193*bxl[0]*c2+1.299038105676658*ezr[14]*c-1.299038105676658*ezl[14]*c-1.14564392373896*ezr[9]*c-1.14564392373896*ezl[9]*c+0.9682458365518543*ezr[5]*c-0.9682458365518543*ezl[5]*c-0.75*ezr[2]*c-0.75*ezl[2]*c+0.4330127018922193*ezr[0]*c-0.4330127018922193*ezl[0]*c; 
  incr[3] = 1.14564392373896*bxr[12]*c2-1.14564392373896*bxl[12]*c2-0.9682458365518543*bxr[7]*c2-0.9682458365518543*bxl[7]*c2+0.75*bxr[3]*c2-0.75*bxl[3]*c2-0.4330127018922193*bxr[1]*c2-0.4330127018922193*bxl[1]*c2-1.14564392373896*ezr[12]*c-1.14564392373896*ezl[12]*c+0.9682458365518543*ezr[7]*c-0.9682458365518543*ezl[7]*c-0.75*ezr[3]*c-0.75*ezl[3]*c+0.4330127018922193*ezr[1]*c-0.4330127018922193*ezl[1]*c; 
  incr[4] = 0.5590169943749475*bxr[10]*c2+0.5590169943749475*bxl[10]*c2-0.4330127018922194*bxr[6]*c2+0.4330127018922194*bxl[6]*c2+0.25*bxr[4]*c2+0.25*bxl[4]*c2-0.5590169943749475*ezr[10]*c+0.5590169943749475*ezl[10]*c+0.4330127018922194*ezr[6]*c+0.4330127018922194*ezl[6]*c-0.25*ezr[4]*c+0.25*ezl[4]*c; 
  incr[5] = 1.677050983124842*bxr[14]*c2+1.677050983124842*bxl[14]*c2-1.479019945774904*bxr[9]*c2+1.479019945774904*bxl[9]*c2+1.25*bxr[5]*c2+1.25*bxl[5]*c2-0.9682458365518543*bxr[2]*c2+0.9682458365518543*bxl[2]*c2+0.5590169943749475*bxr[0]*c2+0.5590169943749475*bxl[0]*c2-1.677050983124842*ezr[14]*c+1.677050983124842*ezl[14]*c+1.479019945774904*ezr[9]*c+1.479019945774904*ezl[9]*c-1.25*ezr[5]*c+1.25*ezl[5]*c+0.9682458365518543*ezr[2]*c+0.9682458365518543*ezl[2]*c-0.5590169943749475*ezr[0]*c+0.5590169943749475*ezl[0]*c; 
  incr[6] = (-0.9682458365518543*bxr[10]*c2)-0.9682458365518543*bxl[10]*c2+0.75*bxr[6]*c2-0.75*bxl[6]*c2-0.4330127018922194*bxr[4]*c2-0.4330127018922194*bxl[4]*c2+0.9682458365518543*ezr[10]*c-0.9682458365518543*ezl[10]*c-0.75*ezr[6]*c-0.75*ezl[6]*c+0.4330127018922194*ezr[4]*c-0.4330127018922194*ezl[4]*c; 
  incr[7] = (-1.479019945774904*bxr[12]*c2)+1.479019945774904*bxl[12]*c2+1.25*bxr[7]*c2+1.25*bxl[7]*c2-0.9682458365518543*bxr[3]*c2+0.9682458365518543*bxl[3]*c2+0.5590169943749476*bxr[1]*c2+0.5590169943749476*bxl[1]*c2+1.479019945774904*ezr[12]*c+1.479019945774904*ezl[12]*c-1.25*ezr[7]*c+1.25*ezl[7]*c+0.9682458365518543*ezr[3]*c+0.9682458365518543*ezl[3]*c-0.5590169943749476*ezr[1]*c+0.5590169943749476*ezl[1]*c; 
  incr[8] = (-0.4330127018922193*bxr[11]*c2)+0.4330127018922193*bxl[11]*c2+0.25*bxr[8]*c2+0.25*bxl[8]*c2+0.4330127018922193*ezr[11]*c+0.4330127018922193*ezl[11]*c-0.25*ezr[8]*c+0.25*ezl[8]*c; 
  incr[9] = (-1.984313483298443*bxr[14]*c2)-1.984313483298443*bxl[14]*c2+1.75*bxr[9]*c2-1.75*bxl[9]*c2-1.479019945774904*bxr[5]*c2-1.479019945774904*bxl[5]*c2+1.14564392373896*bxr[2]*c2-1.14564392373896*bxl[2]*c2-0.6614378277661477*bxr[0]*c2-0.6614378277661477*bxl[0]*c2+1.984313483298443*ezr[14]*c-1.984313483298443*ezl[14]*c-1.75*ezr[9]*c-1.75*ezl[9]*c+1.479019945774904*ezr[5]*c-1.479019945774904*ezl[5]*c-1.14564392373896*ezr[2]*c-1.14564392373896*ezl[2]*c+0.6614378277661477*ezr[0]*c-0.6614378277661477*ezl[0]*c; 
  incr[10] = 1.25*bxr[10]*c2+1.25*bxl[10]*c2-0.9682458365518543*bxr[6]*c2+0.9682458365518543*bxl[6]*c2+0.5590169943749475*bxr[4]*c2+0.5590169943749475*bxl[4]*c2-1.25*ezr[10]*c+1.25*ezl[10]*c+0.9682458365518543*ezr[6]*c+0.9682458365518543*ezl[6]*c-0.5590169943749475*ezr[4]*c+0.5590169943749475*ezl[4]*c; 
  incr[11] = 0.75*bxr[11]*c2-0.75*bxl[11]*c2-0.4330127018922193*bxr[8]*c2-0.4330127018922193*bxl[8]*c2-0.75*ezr[11]*c-0.75*ezl[11]*c+0.4330127018922193*ezr[8]*c-0.4330127018922193*ezl[8]*c; 
  incr[12] = 1.75*bxr[12]*c2-1.75*bxl[12]*c2-1.479019945774904*bxr[7]*c2-1.479019945774904*bxl[7]*c2+1.14564392373896*bxr[3]*c2-1.14564392373896*bxl[3]*c2-0.6614378277661477*bxr[1]*c2-0.6614378277661477*bxl[1]*c2-1.75*ezr[12]*c-1.75*ezl[12]*c+1.479019945774904*ezr[7]*c-1.479019945774904*ezl[7]*c-1.14564392373896*ezr[3]*c-1.14564392373896*ezl[3]*c+0.6614378277661477*ezr[1]*c-0.6614378277661477*ezl[1]*c; 
  incr[13] = 0.25*bxr[13]*c2+0.25*bxl[13]*c2-0.25*ezr[13]*c+0.25*ezl[13]*c; 
  incr[14] = 2.25*bxr[14]*c2+2.25*bxl[14]*c2-1.984313483298443*bxr[9]*c2+1.984313483298443*bxl[9]*c2+1.677050983124842*bxr[5]*c2+1.677050983124842*bxl[5]*c2-1.299038105676658*bxr[2]*c2+1.299038105676658*bxl[2]*c2+0.75*bxr[0]*c2+0.75*bxl[0]*c2-2.25*ezr[14]*c+2.25*ezl[14]*c+1.984313483298443*ezr[9]*c+1.984313483298443*ezl[9]*c-1.677050983124842*ezr[5]*c+1.677050983124842*ezl[5]*c+1.299038105676658*ezr[2]*c+1.299038105676658*ezl[2]*c-0.75*ezr[0]*c+0.75*ezl[0]*c; 

  outEzr[0] += incr[0]*dx1; 
  outEzr[1] += incr[1]*dx1; 
  outEzr[2] += incr[2]*dx1; 
  outEzr[3] += incr[3]*dx1; 
  outEzr[4] += incr[4]*dx1; 
  outEzr[5] += incr[5]*dx1; 
  outEzr[6] += incr[6]*dx1; 
  outEzr[7] += incr[7]*dx1; 
  outEzr[8] += incr[8]*dx1; 
  outEzr[9] += incr[9]*dx1; 
  outEzr[10] += incr[10]*dx1; 
  outEzr[11] += incr[11]*dx1; 
  outEzr[12] += incr[12]*dx1; 
  outEzr[13] += incr[13]*dx1; 
  outEzr[14] += incr[14]*dx1; 

  outEzl[0] += -1.0*incr[0]*dx1; 
  outEzl[1] += -1.0*incr[1]*dx1; 
  outEzl[2] += incr[2]*dx1; 
  outEzl[3] += incr[3]*dx1; 
  outEzl[4] += -1.0*incr[4]*dx1; 
  outEzl[5] += -1.0*incr[5]*dx1; 
  outEzl[6] += incr[6]*dx1; 
  outEzl[7] += -1.0*incr[7]*dx1; 
  outEzl[8] += -1.0*incr[8]*dx1; 
  outEzl[9] += incr[9]*dx1; 
  outEzl[10] += -1.0*incr[10]*dx1; 
  outEzl[11] += incr[11]*dx1; 
  outEzl[12] += incr[12]*dx1; 
  outEzl[13] += -1.0*incr[13]*dx1; 
  outEzl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = (-0.75*bxr[14]*c)+0.75*bxl[14]*c+0.6614378277661477*bxr[9]*c+0.6614378277661477*bxl[9]*c-0.5590169943749475*bxr[5]*c+0.5590169943749475*bxl[5]*c+0.4330127018922193*bxr[2]*c+0.4330127018922193*bxl[2]*c-0.25*bxr[0]*c+0.25*bxl[0]*c+0.75*ezr[14]+0.75*ezl[14]-0.6614378277661477*ezr[9]+0.6614378277661477*ezl[9]+0.5590169943749475*ezr[5]+0.5590169943749475*ezl[5]-0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]+0.25*ezr[0]+0.25*ezl[0]; 
  incr[1] = 0.6614378277661477*bxr[12]*c+0.6614378277661477*bxl[12]*c-0.5590169943749476*bxr[7]*c+0.5590169943749476*bxl[7]*c+0.4330127018922193*bxr[3]*c+0.4330127018922193*bxl[3]*c-0.25*bxr[1]*c+0.25*bxl[1]*c-0.6614378277661477*ezr[12]+0.6614378277661477*ezl[12]+0.5590169943749476*ezr[7]+0.5590169943749476*ezl[7]-0.4330127018922193*ezr[3]+0.4330127018922193*ezl[3]+0.25*ezr[1]+0.25*ezl[1]; 
  incr[2] = 1.299038105676658*bxr[14]*c-1.299038105676658*bxl[14]*c-1.14564392373896*bxr[9]*c-1.14564392373896*bxl[9]*c+0.9682458365518543*bxr[5]*c-0.9682458365518543*bxl[5]*c-0.75*bxr[2]*c-0.75*bxl[2]*c+0.4330127018922193*bxr[0]*c-0.4330127018922193*bxl[0]*c-1.299038105676658*ezr[14]-1.299038105676658*ezl[14]+1.14564392373896*ezr[9]-1.14564392373896*ezl[9]-0.9682458365518543*ezr[5]-0.9682458365518543*ezl[5]+0.75*ezr[2]-0.75*ezl[2]-0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0]; 
  incr[3] = (-1.14564392373896*bxr[12]*c)-1.14564392373896*bxl[12]*c+0.9682458365518543*bxr[7]*c-0.9682458365518543*bxl[7]*c-0.75*bxr[3]*c-0.75*bxl[3]*c+0.4330127018922193*bxr[1]*c-0.4330127018922193*bxl[1]*c+1.14564392373896*ezr[12]-1.14564392373896*ezl[12]-0.9682458365518543*ezr[7]-0.9682458365518543*ezl[7]+0.75*ezr[3]-0.75*ezl[3]-0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]; 
  incr[4] = (-0.5590169943749475*bxr[10]*c)+0.5590169943749475*bxl[10]*c+0.4330127018922194*bxr[6]*c+0.4330127018922194*bxl[6]*c-0.25*bxr[4]*c+0.25*bxl[4]*c+0.5590169943749475*ezr[10]+0.5590169943749475*ezl[10]-0.4330127018922194*ezr[6]+0.4330127018922194*ezl[6]+0.25*ezr[4]+0.25*ezl[4]; 
  incr[5] = (-1.677050983124842*bxr[14]*c)+1.677050983124842*bxl[14]*c+1.479019945774904*bxr[9]*c+1.479019945774904*bxl[9]*c-1.25*bxr[5]*c+1.25*bxl[5]*c+0.9682458365518543*bxr[2]*c+0.9682458365518543*bxl[2]*c-0.5590169943749475*bxr[0]*c+0.5590169943749475*bxl[0]*c+1.677050983124842*ezr[14]+1.677050983124842*ezl[14]-1.479019945774904*ezr[9]+1.479019945774904*ezl[9]+1.25*ezr[5]+1.25*ezl[5]-0.9682458365518543*ezr[2]+0.9682458365518543*ezl[2]+0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0]; 
  incr[6] = 0.9682458365518543*bxr[10]*c-0.9682458365518543*bxl[10]*c-0.75*bxr[6]*c-0.75*bxl[6]*c+0.4330127018922194*bxr[4]*c-0.4330127018922194*bxl[4]*c-0.9682458365518543*ezr[10]-0.9682458365518543*ezl[10]+0.75*ezr[6]-0.75*ezl[6]-0.4330127018922194*ezr[4]-0.4330127018922194*ezl[4]; 
  incr[7] = 1.479019945774904*bxr[12]*c+1.479019945774904*bxl[12]*c-1.25*bxr[7]*c+1.25*bxl[7]*c+0.9682458365518543*bxr[3]*c+0.9682458365518543*bxl[3]*c-0.5590169943749476*bxr[1]*c+0.5590169943749476*bxl[1]*c-1.479019945774904*ezr[12]+1.479019945774904*ezl[12]+1.25*ezr[7]+1.25*ezl[7]-0.9682458365518543*ezr[3]+0.9682458365518543*ezl[3]+0.5590169943749476*ezr[1]+0.5590169943749476*ezl[1]; 
  incr[8] = 0.4330127018922193*bxr[11]*c+0.4330127018922193*bxl[11]*c-0.25*bxr[8]*c+0.25*bxl[8]*c-0.4330127018922193*ezr[11]+0.4330127018922193*ezl[11]+0.25*ezr[8]+0.25*ezl[8]; 
  incr[9] = 1.984313483298443*bxr[14]*c-1.984313483298443*bxl[14]*c-1.75*bxr[9]*c-1.75*bxl[9]*c+1.479019945774904*bxr[5]*c-1.479019945774904*bxl[5]*c-1.14564392373896*bxr[2]*c-1.14564392373896*bxl[2]*c+0.6614378277661477*bxr[0]*c-0.6614378277661477*bxl[0]*c-1.984313483298443*ezr[14]-1.984313483298443*ezl[14]+1.75*ezr[9]-1.75*ezl[9]-1.479019945774904*ezr[5]-1.479019945774904*ezl[5]+1.14564392373896*ezr[2]-1.14564392373896*ezl[2]-0.6614378277661477*ezr[0]-0.6614378277661477*ezl[0]; 
  incr[10] = (-1.25*bxr[10]*c)+1.25*bxl[10]*c+0.9682458365518543*bxr[6]*c+0.9682458365518543*bxl[6]*c-0.5590169943749475*bxr[4]*c+0.5590169943749475*bxl[4]*c+1.25*ezr[10]+1.25*ezl[10]-0.9682458365518543*ezr[6]+0.9682458365518543*ezl[6]+0.5590169943749475*ezr[4]+0.5590169943749475*ezl[4]; 
  incr[11] = (-0.75*bxr[11]*c)-0.75*bxl[11]*c+0.4330127018922193*bxr[8]*c-0.4330127018922193*bxl[8]*c+0.75*ezr[11]-0.75*ezl[11]-0.4330127018922193*ezr[8]-0.4330127018922193*ezl[8]; 
  incr[12] = (-1.75*bxr[12]*c)-1.75*bxl[12]*c+1.479019945774904*bxr[7]*c-1.479019945774904*bxl[7]*c-1.14564392373896*bxr[3]*c-1.14564392373896*bxl[3]*c+0.6614378277661477*bxr[1]*c-0.6614378277661477*bxl[1]*c+1.75*ezr[12]-1.75*ezl[12]-1.479019945774904*ezr[7]-1.479019945774904*ezl[7]+1.14564392373896*ezr[3]-1.14564392373896*ezl[3]-0.6614378277661477*ezr[1]-0.6614378277661477*ezl[1]; 
  incr[13] = (-0.25*bxr[13]*c)+0.25*bxl[13]*c+0.25*ezr[13]+0.25*ezl[13]; 
  incr[14] = (-2.25*bxr[14]*c)+2.25*bxl[14]*c+1.984313483298443*bxr[9]*c+1.984313483298443*bxl[9]*c-1.677050983124842*bxr[5]*c+1.677050983124842*bxl[5]*c+1.299038105676658*bxr[2]*c+1.299038105676658*bxl[2]*c-0.75*bxr[0]*c+0.75*bxl[0]*c+2.25*ezr[14]+2.25*ezl[14]-1.984313483298443*ezr[9]+1.984313483298443*ezl[9]+1.677050983124842*ezr[5]+1.677050983124842*ezl[5]-1.299038105676658*ezr[2]+1.299038105676658*ezl[2]+0.75*ezr[0]+0.75*ezl[0]; 

  outBxr[0] += incr[0]*dx1; 
  outBxr[1] += incr[1]*dx1; 
  outBxr[2] += incr[2]*dx1; 
  outBxr[3] += incr[3]*dx1; 
  outBxr[4] += incr[4]*dx1; 
  outBxr[5] += incr[5]*dx1; 
  outBxr[6] += incr[6]*dx1; 
  outBxr[7] += incr[7]*dx1; 
  outBxr[8] += incr[8]*dx1; 
  outBxr[9] += incr[9]*dx1; 
  outBxr[10] += incr[10]*dx1; 
  outBxr[11] += incr[11]*dx1; 
  outBxr[12] += incr[12]*dx1; 
  outBxr[13] += incr[13]*dx1; 
  outBxr[14] += incr[14]*dx1; 

  outBxl[0] += -1.0*incr[0]*dx1; 
  outBxl[1] += -1.0*incr[1]*dx1; 
  outBxl[2] += incr[2]*dx1; 
  outBxl[3] += incr[3]*dx1; 
  outBxl[4] += -1.0*incr[4]*dx1; 
  outBxl[5] += -1.0*incr[5]*dx1; 
  outBxl[6] += incr[6]*dx1; 
  outBxl[7] += -1.0*incr[7]*dx1; 
  outBxl[8] += -1.0*incr[8]*dx1; 
  outBxl[9] += incr[9]*dx1; 
  outBxl[10] += -1.0*incr[10]*dx1; 
  outBxl[11] += incr[11]*dx1; 
  outBxl[12] += incr[12]*dx1; 
  outBxl[13] += -1.0*incr[13]*dx1; 
  outBxl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = (-0.75*byr[14]*c*gamma)+0.75*byl[14]*c*gamma+0.6614378277661477*byr[9]*c*gamma+0.6614378277661477*byl[9]*c*gamma-0.5590169943749475*byr[5]*c*gamma+0.5590169943749475*byl[5]*c*gamma+0.4330127018922193*byr[2]*c*gamma+0.4330127018922193*byl[2]*c*gamma-0.25*byr[0]*c*gamma+0.25*byl[0]*c*gamma+0.75*psr[14]*gamma+0.75*psl[14]*gamma-0.6614378277661477*psr[9]*gamma+0.6614378277661477*psl[9]*gamma+0.5590169943749475*psr[5]*gamma+0.5590169943749475*psl[5]*gamma-0.4330127018922193*psr[2]*gamma+0.4330127018922193*psl[2]*gamma+0.25*psr[0]*gamma+0.25*psl[0]*gamma; 
  incr[1] = 0.6614378277661477*byr[12]*c*gamma+0.6614378277661477*byl[12]*c*gamma-0.5590169943749476*byr[7]*c*gamma+0.5590169943749476*byl[7]*c*gamma+0.4330127018922193*byr[3]*c*gamma+0.4330127018922193*byl[3]*c*gamma-0.25*byr[1]*c*gamma+0.25*byl[1]*c*gamma-0.6614378277661477*psr[12]*gamma+0.6614378277661477*psl[12]*gamma+0.5590169943749476*psr[7]*gamma+0.5590169943749476*psl[7]*gamma-0.4330127018922193*psr[3]*gamma+0.4330127018922193*psl[3]*gamma+0.25*psr[1]*gamma+0.25*psl[1]*gamma; 
  incr[2] = 1.299038105676658*byr[14]*c*gamma-1.299038105676658*byl[14]*c*gamma-1.14564392373896*byr[9]*c*gamma-1.14564392373896*byl[9]*c*gamma+0.9682458365518543*byr[5]*c*gamma-0.9682458365518543*byl[5]*c*gamma-0.75*byr[2]*c*gamma-0.75*byl[2]*c*gamma+0.4330127018922193*byr[0]*c*gamma-0.4330127018922193*byl[0]*c*gamma-1.299038105676658*psr[14]*gamma-1.299038105676658*psl[14]*gamma+1.14564392373896*psr[9]*gamma-1.14564392373896*psl[9]*gamma-0.9682458365518543*psr[5]*gamma-0.9682458365518543*psl[5]*gamma+0.75*psr[2]*gamma-0.75*psl[2]*gamma-0.4330127018922193*psr[0]*gamma-0.4330127018922193*psl[0]*gamma; 
  incr[3] = (-1.14564392373896*byr[12]*c*gamma)-1.14564392373896*byl[12]*c*gamma+0.9682458365518543*byr[7]*c*gamma-0.9682458365518543*byl[7]*c*gamma-0.75*byr[3]*c*gamma-0.75*byl[3]*c*gamma+0.4330127018922193*byr[1]*c*gamma-0.4330127018922193*byl[1]*c*gamma+1.14564392373896*psr[12]*gamma-1.14564392373896*psl[12]*gamma-0.9682458365518543*psr[7]*gamma-0.9682458365518543*psl[7]*gamma+0.75*psr[3]*gamma-0.75*psl[3]*gamma-0.4330127018922193*psr[1]*gamma-0.4330127018922193*psl[1]*gamma; 
  incr[4] = (-0.5590169943749475*byr[10]*c*gamma)+0.5590169943749475*byl[10]*c*gamma+0.4330127018922194*byr[6]*c*gamma+0.4330127018922194*byl[6]*c*gamma-0.25*byr[4]*c*gamma+0.25*byl[4]*c*gamma+0.5590169943749475*psr[10]*gamma+0.5590169943749475*psl[10]*gamma-0.4330127018922194*psr[6]*gamma+0.4330127018922194*psl[6]*gamma+0.25*psr[4]*gamma+0.25*psl[4]*gamma; 
  incr[5] = (-1.677050983124842*byr[14]*c*gamma)+1.677050983124842*byl[14]*c*gamma+1.479019945774904*byr[9]*c*gamma+1.479019945774904*byl[9]*c*gamma-1.25*byr[5]*c*gamma+1.25*byl[5]*c*gamma+0.9682458365518543*byr[2]*c*gamma+0.9682458365518543*byl[2]*c*gamma-0.5590169943749475*byr[0]*c*gamma+0.5590169943749475*byl[0]*c*gamma+1.677050983124842*psr[14]*gamma+1.677050983124842*psl[14]*gamma-1.479019945774904*psr[9]*gamma+1.479019945774904*psl[9]*gamma+1.25*psr[5]*gamma+1.25*psl[5]*gamma-0.9682458365518543*psr[2]*gamma+0.9682458365518543*psl[2]*gamma+0.5590169943749475*psr[0]*gamma+0.5590169943749475*psl[0]*gamma; 
  incr[6] = 0.9682458365518543*byr[10]*c*gamma-0.9682458365518543*byl[10]*c*gamma-0.75*byr[6]*c*gamma-0.75*byl[6]*c*gamma+0.4330127018922194*byr[4]*c*gamma-0.4330127018922194*byl[4]*c*gamma-0.9682458365518543*psr[10]*gamma-0.9682458365518543*psl[10]*gamma+0.75*psr[6]*gamma-0.75*psl[6]*gamma-0.4330127018922194*psr[4]*gamma-0.4330127018922194*psl[4]*gamma; 
  incr[7] = 1.479019945774904*byr[12]*c*gamma+1.479019945774904*byl[12]*c*gamma-1.25*byr[7]*c*gamma+1.25*byl[7]*c*gamma+0.9682458365518543*byr[3]*c*gamma+0.9682458365518543*byl[3]*c*gamma-0.5590169943749476*byr[1]*c*gamma+0.5590169943749476*byl[1]*c*gamma-1.479019945774904*psr[12]*gamma+1.479019945774904*psl[12]*gamma+1.25*psr[7]*gamma+1.25*psl[7]*gamma-0.9682458365518543*psr[3]*gamma+0.9682458365518543*psl[3]*gamma+0.5590169943749476*psr[1]*gamma+0.5590169943749476*psl[1]*gamma; 
  incr[8] = 0.4330127018922193*byr[11]*c*gamma+0.4330127018922193*byl[11]*c*gamma-0.25*byr[8]*c*gamma+0.25*byl[8]*c*gamma-0.4330127018922193*psr[11]*gamma+0.4330127018922193*psl[11]*gamma+0.25*psr[8]*gamma+0.25*psl[8]*gamma; 
  incr[9] = 1.984313483298443*byr[14]*c*gamma-1.984313483298443*byl[14]*c*gamma-1.75*byr[9]*c*gamma-1.75*byl[9]*c*gamma+1.479019945774904*byr[5]*c*gamma-1.479019945774904*byl[5]*c*gamma-1.14564392373896*byr[2]*c*gamma-1.14564392373896*byl[2]*c*gamma+0.6614378277661477*byr[0]*c*gamma-0.6614378277661477*byl[0]*c*gamma-1.984313483298443*psr[14]*gamma-1.984313483298443*psl[14]*gamma+1.75*psr[9]*gamma-1.75*psl[9]*gamma-1.479019945774904*psr[5]*gamma-1.479019945774904*psl[5]*gamma+1.14564392373896*psr[2]*gamma-1.14564392373896*psl[2]*gamma-0.6614378277661477*psr[0]*gamma-0.6614378277661477*psl[0]*gamma; 
  incr[10] = (-1.25*byr[10]*c*gamma)+1.25*byl[10]*c*gamma+0.9682458365518543*byr[6]*c*gamma+0.9682458365518543*byl[6]*c*gamma-0.5590169943749475*byr[4]*c*gamma+0.5590169943749475*byl[4]*c*gamma+1.25*psr[10]*gamma+1.25*psl[10]*gamma-0.9682458365518543*psr[6]*gamma+0.9682458365518543*psl[6]*gamma+0.5590169943749475*psr[4]*gamma+0.5590169943749475*psl[4]*gamma; 
  incr[11] = (-0.75*byr[11]*c*gamma)-0.75*byl[11]*c*gamma+0.4330127018922193*byr[8]*c*gamma-0.4330127018922193*byl[8]*c*gamma+0.75*psr[11]*gamma-0.75*psl[11]*gamma-0.4330127018922193*psr[8]*gamma-0.4330127018922193*psl[8]*gamma; 
  incr[12] = (-1.75*byr[12]*c*gamma)-1.75*byl[12]*c*gamma+1.479019945774904*byr[7]*c*gamma-1.479019945774904*byl[7]*c*gamma-1.14564392373896*byr[3]*c*gamma-1.14564392373896*byl[3]*c*gamma+0.6614378277661477*byr[1]*c*gamma-0.6614378277661477*byl[1]*c*gamma+1.75*psr[12]*gamma-1.75*psl[12]*gamma-1.479019945774904*psr[7]*gamma-1.479019945774904*psl[7]*gamma+1.14564392373896*psr[3]*gamma-1.14564392373896*psl[3]*gamma-0.6614378277661477*psr[1]*gamma-0.6614378277661477*psl[1]*gamma; 
  incr[13] = (-0.25*byr[13]*c*gamma)+0.25*byl[13]*c*gamma+0.25*psr[13]*gamma+0.25*psl[13]*gamma; 
  incr[14] = (-2.25*byr[14]*c*gamma)+2.25*byl[14]*c*gamma+1.984313483298443*byr[9]*c*gamma+1.984313483298443*byl[9]*c*gamma-1.677050983124842*byr[5]*c*gamma+1.677050983124842*byl[5]*c*gamma+1.299038105676658*byr[2]*c*gamma+1.299038105676658*byl[2]*c*gamma-0.75*byr[0]*c*gamma+0.75*byl[0]*c*gamma+2.25*psr[14]*gamma+2.25*psl[14]*gamma-1.984313483298443*psr[9]*gamma+1.984313483298443*psl[9]*gamma+1.677050983124842*psr[5]*gamma+1.677050983124842*psl[5]*gamma-1.299038105676658*psr[2]*gamma+1.299038105676658*psl[2]*gamma+0.75*psr[0]*gamma+0.75*psl[0]*gamma; 

  outByr[0] += incr[0]*dx1; 
  outByr[1] += incr[1]*dx1; 
  outByr[2] += incr[2]*dx1; 
  outByr[3] += incr[3]*dx1; 
  outByr[4] += incr[4]*dx1; 
  outByr[5] += incr[5]*dx1; 
  outByr[6] += incr[6]*dx1; 
  outByr[7] += incr[7]*dx1; 
  outByr[8] += incr[8]*dx1; 
  outByr[9] += incr[9]*dx1; 
  outByr[10] += incr[10]*dx1; 
  outByr[11] += incr[11]*dx1; 
  outByr[12] += incr[12]*dx1; 
  outByr[13] += incr[13]*dx1; 
  outByr[14] += incr[14]*dx1; 

  outByl[0] += -1.0*incr[0]*dx1; 
  outByl[1] += -1.0*incr[1]*dx1; 
  outByl[2] += incr[2]*dx1; 
  outByl[3] += incr[3]*dx1; 
  outByl[4] += -1.0*incr[4]*dx1; 
  outByl[5] += -1.0*incr[5]*dx1; 
  outByl[6] += incr[6]*dx1; 
  outByl[7] += -1.0*incr[7]*dx1; 
  outByl[8] += -1.0*incr[8]*dx1; 
  outByl[9] += incr[9]*dx1; 
  outByl[10] += -1.0*incr[10]*dx1; 
  outByl[11] += incr[11]*dx1; 
  outByl[12] += incr[12]*dx1; 
  outByl[13] += -1.0*incr[13]*dx1; 
  outByl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = (-0.75*bzr[14]*c)+0.75*bzl[14]*c+0.6614378277661477*bzr[9]*c+0.6614378277661477*bzl[9]*c-0.5590169943749475*bzr[5]*c+0.5590169943749475*bzl[5]*c+0.4330127018922193*bzr[2]*c+0.4330127018922193*bzl[2]*c-0.25*bzr[0]*c+0.25*bzl[0]*c-0.75*exr[14]-0.75*exl[14]+0.6614378277661477*exr[9]-0.6614378277661477*exl[9]-0.5590169943749475*exr[5]-0.5590169943749475*exl[5]+0.4330127018922193*exr[2]-0.4330127018922193*exl[2]-0.25*exr[0]-0.25*exl[0]; 
  incr[1] = 0.6614378277661477*bzr[12]*c+0.6614378277661477*bzl[12]*c-0.5590169943749476*bzr[7]*c+0.5590169943749476*bzl[7]*c+0.4330127018922193*bzr[3]*c+0.4330127018922193*bzl[3]*c-0.25*bzr[1]*c+0.25*bzl[1]*c+0.6614378277661477*exr[12]-0.6614378277661477*exl[12]-0.5590169943749476*exr[7]-0.5590169943749476*exl[7]+0.4330127018922193*exr[3]-0.4330127018922193*exl[3]-0.25*exr[1]-0.25*exl[1]; 
  incr[2] = 1.299038105676658*bzr[14]*c-1.299038105676658*bzl[14]*c-1.14564392373896*bzr[9]*c-1.14564392373896*bzl[9]*c+0.9682458365518543*bzr[5]*c-0.9682458365518543*bzl[5]*c-0.75*bzr[2]*c-0.75*bzl[2]*c+0.4330127018922193*bzr[0]*c-0.4330127018922193*bzl[0]*c+1.299038105676658*exr[14]+1.299038105676658*exl[14]-1.14564392373896*exr[9]+1.14564392373896*exl[9]+0.9682458365518543*exr[5]+0.9682458365518543*exl[5]-0.75*exr[2]+0.75*exl[2]+0.4330127018922193*exr[0]+0.4330127018922193*exl[0]; 
  incr[3] = (-1.14564392373896*bzr[12]*c)-1.14564392373896*bzl[12]*c+0.9682458365518543*bzr[7]*c-0.9682458365518543*bzl[7]*c-0.75*bzr[3]*c-0.75*bzl[3]*c+0.4330127018922193*bzr[1]*c-0.4330127018922193*bzl[1]*c-1.14564392373896*exr[12]+1.14564392373896*exl[12]+0.9682458365518543*exr[7]+0.9682458365518543*exl[7]-0.75*exr[3]+0.75*exl[3]+0.4330127018922193*exr[1]+0.4330127018922193*exl[1]; 
  incr[4] = (-0.5590169943749475*bzr[10]*c)+0.5590169943749475*bzl[10]*c+0.4330127018922194*bzr[6]*c+0.4330127018922194*bzl[6]*c-0.25*bzr[4]*c+0.25*bzl[4]*c-0.5590169943749475*exr[10]-0.5590169943749475*exl[10]+0.4330127018922194*exr[6]-0.4330127018922194*exl[6]-0.25*exr[4]-0.25*exl[4]; 
  incr[5] = (-1.677050983124842*bzr[14]*c)+1.677050983124842*bzl[14]*c+1.479019945774904*bzr[9]*c+1.479019945774904*bzl[9]*c-1.25*bzr[5]*c+1.25*bzl[5]*c+0.9682458365518543*bzr[2]*c+0.9682458365518543*bzl[2]*c-0.5590169943749475*bzr[0]*c+0.5590169943749475*bzl[0]*c-1.677050983124842*exr[14]-1.677050983124842*exl[14]+1.479019945774904*exr[9]-1.479019945774904*exl[9]-1.25*exr[5]-1.25*exl[5]+0.9682458365518543*exr[2]-0.9682458365518543*exl[2]-0.5590169943749475*exr[0]-0.5590169943749475*exl[0]; 
  incr[6] = 0.9682458365518543*bzr[10]*c-0.9682458365518543*bzl[10]*c-0.75*bzr[6]*c-0.75*bzl[6]*c+0.4330127018922194*bzr[4]*c-0.4330127018922194*bzl[4]*c+0.9682458365518543*exr[10]+0.9682458365518543*exl[10]-0.75*exr[6]+0.75*exl[6]+0.4330127018922194*exr[4]+0.4330127018922194*exl[4]; 
  incr[7] = 1.479019945774904*bzr[12]*c+1.479019945774904*bzl[12]*c-1.25*bzr[7]*c+1.25*bzl[7]*c+0.9682458365518543*bzr[3]*c+0.9682458365518543*bzl[3]*c-0.5590169943749476*bzr[1]*c+0.5590169943749476*bzl[1]*c+1.479019945774904*exr[12]-1.479019945774904*exl[12]-1.25*exr[7]-1.25*exl[7]+0.9682458365518543*exr[3]-0.9682458365518543*exl[3]-0.5590169943749476*exr[1]-0.5590169943749476*exl[1]; 
  incr[8] = 0.4330127018922193*bzr[11]*c+0.4330127018922193*bzl[11]*c-0.25*bzr[8]*c+0.25*bzl[8]*c+0.4330127018922193*exr[11]-0.4330127018922193*exl[11]-0.25*exr[8]-0.25*exl[8]; 
  incr[9] = 1.984313483298443*bzr[14]*c-1.984313483298443*bzl[14]*c-1.75*bzr[9]*c-1.75*bzl[9]*c+1.479019945774904*bzr[5]*c-1.479019945774904*bzl[5]*c-1.14564392373896*bzr[2]*c-1.14564392373896*bzl[2]*c+0.6614378277661477*bzr[0]*c-0.6614378277661477*bzl[0]*c+1.984313483298443*exr[14]+1.984313483298443*exl[14]-1.75*exr[9]+1.75*exl[9]+1.479019945774904*exr[5]+1.479019945774904*exl[5]-1.14564392373896*exr[2]+1.14564392373896*exl[2]+0.6614378277661477*exr[0]+0.6614378277661477*exl[0]; 
  incr[10] = (-1.25*bzr[10]*c)+1.25*bzl[10]*c+0.9682458365518543*bzr[6]*c+0.9682458365518543*bzl[6]*c-0.5590169943749475*bzr[4]*c+0.5590169943749475*bzl[4]*c-1.25*exr[10]-1.25*exl[10]+0.9682458365518543*exr[6]-0.9682458365518543*exl[6]-0.5590169943749475*exr[4]-0.5590169943749475*exl[4]; 
  incr[11] = (-0.75*bzr[11]*c)-0.75*bzl[11]*c+0.4330127018922193*bzr[8]*c-0.4330127018922193*bzl[8]*c-0.75*exr[11]+0.75*exl[11]+0.4330127018922193*exr[8]+0.4330127018922193*exl[8]; 
  incr[12] = (-1.75*bzr[12]*c)-1.75*bzl[12]*c+1.479019945774904*bzr[7]*c-1.479019945774904*bzl[7]*c-1.14564392373896*bzr[3]*c-1.14564392373896*bzl[3]*c+0.6614378277661477*bzr[1]*c-0.6614378277661477*bzl[1]*c-1.75*exr[12]+1.75*exl[12]+1.479019945774904*exr[7]+1.479019945774904*exl[7]-1.14564392373896*exr[3]+1.14564392373896*exl[3]+0.6614378277661477*exr[1]+0.6614378277661477*exl[1]; 
  incr[13] = (-0.25*bzr[13]*c)+0.25*bzl[13]*c-0.25*exr[13]-0.25*exl[13]; 
  incr[14] = (-2.25*bzr[14]*c)+2.25*bzl[14]*c+1.984313483298443*bzr[9]*c+1.984313483298443*bzl[9]*c-1.677050983124842*bzr[5]*c+1.677050983124842*bzl[5]*c+1.299038105676658*bzr[2]*c+1.299038105676658*bzl[2]*c-0.75*bzr[0]*c+0.75*bzl[0]*c-2.25*exr[14]-2.25*exl[14]+1.984313483298443*exr[9]-1.984313483298443*exl[9]-1.677050983124842*exr[5]-1.677050983124842*exl[5]+1.299038105676658*exr[2]-1.299038105676658*exl[2]-0.75*exr[0]-0.75*exl[0]; 

  outBzr[0] += incr[0]*dx1; 
  outBzr[1] += incr[1]*dx1; 
  outBzr[2] += incr[2]*dx1; 
  outBzr[3] += incr[3]*dx1; 
  outBzr[4] += incr[4]*dx1; 
  outBzr[5] += incr[5]*dx1; 
  outBzr[6] += incr[6]*dx1; 
  outBzr[7] += incr[7]*dx1; 
  outBzr[8] += incr[8]*dx1; 
  outBzr[9] += incr[9]*dx1; 
  outBzr[10] += incr[10]*dx1; 
  outBzr[11] += incr[11]*dx1; 
  outBzr[12] += incr[12]*dx1; 
  outBzr[13] += incr[13]*dx1; 
  outBzr[14] += incr[14]*dx1; 

  outBzl[0] += -1.0*incr[0]*dx1; 
  outBzl[1] += -1.0*incr[1]*dx1; 
  outBzl[2] += incr[2]*dx1; 
  outBzl[3] += incr[3]*dx1; 
  outBzl[4] += -1.0*incr[4]*dx1; 
  outBzl[5] += -1.0*incr[5]*dx1; 
  outBzl[6] += incr[6]*dx1; 
  outBzl[7] += -1.0*incr[7]*dx1; 
  outBzl[8] += -1.0*incr[8]*dx1; 
  outBzl[9] += incr[9]*dx1; 
  outBzl[10] += -1.0*incr[10]*dx1; 
  outBzl[11] += incr[11]*dx1; 
  outBzl[12] += incr[12]*dx1; 
  outBzl[13] += -1.0*incr[13]*dx1; 
  outBzl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = (-0.75*phr[14]*c*chi)+0.75*phl[14]*c*chi+0.6614378277661477*phr[9]*c*chi+0.6614378277661477*phl[9]*c*chi-0.5590169943749475*phr[5]*c*chi+0.5590169943749475*phl[5]*c*chi+0.4330127018922193*phr[2]*c*chi+0.4330127018922193*phl[2]*c*chi-0.25*phr[0]*c*chi+0.25*phl[0]*c*chi+0.75*eyr[14]*chi+0.75*eyl[14]*chi-0.6614378277661477*eyr[9]*chi+0.6614378277661477*eyl[9]*chi+0.5590169943749475*eyr[5]*chi+0.5590169943749475*eyl[5]*chi-0.4330127018922193*eyr[2]*chi+0.4330127018922193*eyl[2]*chi+0.25*eyr[0]*chi+0.25*eyl[0]*chi; 
  incr[1] = 0.6614378277661477*phr[12]*c*chi+0.6614378277661477*phl[12]*c*chi-0.5590169943749476*phr[7]*c*chi+0.5590169943749476*phl[7]*c*chi+0.4330127018922193*phr[3]*c*chi+0.4330127018922193*phl[3]*c*chi-0.25*phr[1]*c*chi+0.25*phl[1]*c*chi-0.6614378277661477*eyr[12]*chi+0.6614378277661477*eyl[12]*chi+0.5590169943749476*eyr[7]*chi+0.5590169943749476*eyl[7]*chi-0.4330127018922193*eyr[3]*chi+0.4330127018922193*eyl[3]*chi+0.25*eyr[1]*chi+0.25*eyl[1]*chi; 
  incr[2] = 1.299038105676658*phr[14]*c*chi-1.299038105676658*phl[14]*c*chi-1.14564392373896*phr[9]*c*chi-1.14564392373896*phl[9]*c*chi+0.9682458365518543*phr[5]*c*chi-0.9682458365518543*phl[5]*c*chi-0.75*phr[2]*c*chi-0.75*phl[2]*c*chi+0.4330127018922193*phr[0]*c*chi-0.4330127018922193*phl[0]*c*chi-1.299038105676658*eyr[14]*chi-1.299038105676658*eyl[14]*chi+1.14564392373896*eyr[9]*chi-1.14564392373896*eyl[9]*chi-0.9682458365518543*eyr[5]*chi-0.9682458365518543*eyl[5]*chi+0.75*eyr[2]*chi-0.75*eyl[2]*chi-0.4330127018922193*eyr[0]*chi-0.4330127018922193*eyl[0]*chi; 
  incr[3] = (-1.14564392373896*phr[12]*c*chi)-1.14564392373896*phl[12]*c*chi+0.9682458365518543*phr[7]*c*chi-0.9682458365518543*phl[7]*c*chi-0.75*phr[3]*c*chi-0.75*phl[3]*c*chi+0.4330127018922193*phr[1]*c*chi-0.4330127018922193*phl[1]*c*chi+1.14564392373896*eyr[12]*chi-1.14564392373896*eyl[12]*chi-0.9682458365518543*eyr[7]*chi-0.9682458365518543*eyl[7]*chi+0.75*eyr[3]*chi-0.75*eyl[3]*chi-0.4330127018922193*eyr[1]*chi-0.4330127018922193*eyl[1]*chi; 
  incr[4] = (-0.5590169943749475*phr[10]*c*chi)+0.5590169943749475*phl[10]*c*chi+0.4330127018922194*phr[6]*c*chi+0.4330127018922194*phl[6]*c*chi-0.25*phr[4]*c*chi+0.25*phl[4]*c*chi+0.5590169943749475*eyr[10]*chi+0.5590169943749475*eyl[10]*chi-0.4330127018922194*eyr[6]*chi+0.4330127018922194*eyl[6]*chi+0.25*eyr[4]*chi+0.25*eyl[4]*chi; 
  incr[5] = (-1.677050983124842*phr[14]*c*chi)+1.677050983124842*phl[14]*c*chi+1.479019945774904*phr[9]*c*chi+1.479019945774904*phl[9]*c*chi-1.25*phr[5]*c*chi+1.25*phl[5]*c*chi+0.9682458365518543*phr[2]*c*chi+0.9682458365518543*phl[2]*c*chi-0.5590169943749475*phr[0]*c*chi+0.5590169943749475*phl[0]*c*chi+1.677050983124842*eyr[14]*chi+1.677050983124842*eyl[14]*chi-1.479019945774904*eyr[9]*chi+1.479019945774904*eyl[9]*chi+1.25*eyr[5]*chi+1.25*eyl[5]*chi-0.9682458365518543*eyr[2]*chi+0.9682458365518543*eyl[2]*chi+0.5590169943749475*eyr[0]*chi+0.5590169943749475*eyl[0]*chi; 
  incr[6] = 0.9682458365518543*phr[10]*c*chi-0.9682458365518543*phl[10]*c*chi-0.75*phr[6]*c*chi-0.75*phl[6]*c*chi+0.4330127018922194*phr[4]*c*chi-0.4330127018922194*phl[4]*c*chi-0.9682458365518543*eyr[10]*chi-0.9682458365518543*eyl[10]*chi+0.75*eyr[6]*chi-0.75*eyl[6]*chi-0.4330127018922194*eyr[4]*chi-0.4330127018922194*eyl[4]*chi; 
  incr[7] = 1.479019945774904*phr[12]*c*chi+1.479019945774904*phl[12]*c*chi-1.25*phr[7]*c*chi+1.25*phl[7]*c*chi+0.9682458365518543*phr[3]*c*chi+0.9682458365518543*phl[3]*c*chi-0.5590169943749476*phr[1]*c*chi+0.5590169943749476*phl[1]*c*chi-1.479019945774904*eyr[12]*chi+1.479019945774904*eyl[12]*chi+1.25*eyr[7]*chi+1.25*eyl[7]*chi-0.9682458365518543*eyr[3]*chi+0.9682458365518543*eyl[3]*chi+0.5590169943749476*eyr[1]*chi+0.5590169943749476*eyl[1]*chi; 
  incr[8] = 0.4330127018922193*phr[11]*c*chi+0.4330127018922193*phl[11]*c*chi-0.25*phr[8]*c*chi+0.25*phl[8]*c*chi-0.4330127018922193*eyr[11]*chi+0.4330127018922193*eyl[11]*chi+0.25*eyr[8]*chi+0.25*eyl[8]*chi; 
  incr[9] = 1.984313483298443*phr[14]*c*chi-1.984313483298443*phl[14]*c*chi-1.75*phr[9]*c*chi-1.75*phl[9]*c*chi+1.479019945774904*phr[5]*c*chi-1.479019945774904*phl[5]*c*chi-1.14564392373896*phr[2]*c*chi-1.14564392373896*phl[2]*c*chi+0.6614378277661477*phr[0]*c*chi-0.6614378277661477*phl[0]*c*chi-1.984313483298443*eyr[14]*chi-1.984313483298443*eyl[14]*chi+1.75*eyr[9]*chi-1.75*eyl[9]*chi-1.479019945774904*eyr[5]*chi-1.479019945774904*eyl[5]*chi+1.14564392373896*eyr[2]*chi-1.14564392373896*eyl[2]*chi-0.6614378277661477*eyr[0]*chi-0.6614378277661477*eyl[0]*chi; 
  incr[10] = (-1.25*phr[10]*c*chi)+1.25*phl[10]*c*chi+0.9682458365518543*phr[6]*c*chi+0.9682458365518543*phl[6]*c*chi-0.5590169943749475*phr[4]*c*chi+0.5590169943749475*phl[4]*c*chi+1.25*eyr[10]*chi+1.25*eyl[10]*chi-0.9682458365518543*eyr[6]*chi+0.9682458365518543*eyl[6]*chi+0.5590169943749475*eyr[4]*chi+0.5590169943749475*eyl[4]*chi; 
  incr[11] = (-0.75*phr[11]*c*chi)-0.75*phl[11]*c*chi+0.4330127018922193*phr[8]*c*chi-0.4330127018922193*phl[8]*c*chi+0.75*eyr[11]*chi-0.75*eyl[11]*chi-0.4330127018922193*eyr[8]*chi-0.4330127018922193*eyl[8]*chi; 
  incr[12] = (-1.75*phr[12]*c*chi)-1.75*phl[12]*c*chi+1.479019945774904*phr[7]*c*chi-1.479019945774904*phl[7]*c*chi-1.14564392373896*phr[3]*c*chi-1.14564392373896*phl[3]*c*chi+0.6614378277661477*phr[1]*c*chi-0.6614378277661477*phl[1]*c*chi+1.75*eyr[12]*chi-1.75*eyl[12]*chi-1.479019945774904*eyr[7]*chi-1.479019945774904*eyl[7]*chi+1.14564392373896*eyr[3]*chi-1.14564392373896*eyl[3]*chi-0.6614378277661477*eyr[1]*chi-0.6614378277661477*eyl[1]*chi; 
  incr[13] = (-0.25*phr[13]*c*chi)+0.25*phl[13]*c*chi+0.25*eyr[13]*chi+0.25*eyl[13]*chi; 
  incr[14] = (-2.25*phr[14]*c*chi)+2.25*phl[14]*c*chi+1.984313483298443*phr[9]*c*chi+1.984313483298443*phl[9]*c*chi-1.677050983124842*phr[5]*c*chi+1.677050983124842*phl[5]*c*chi+1.299038105676658*phr[2]*c*chi+1.299038105676658*phl[2]*c*chi-0.75*phr[0]*c*chi+0.75*phl[0]*c*chi+2.25*eyr[14]*chi+2.25*eyl[14]*chi-1.984313483298443*eyr[9]*chi+1.984313483298443*eyl[9]*chi+1.677050983124842*eyr[5]*chi+1.677050983124842*eyl[5]*chi-1.299038105676658*eyr[2]*chi+1.299038105676658*eyl[2]*chi+0.75*eyr[0]*chi+0.75*eyl[0]*chi; 

  outPhr[0] += incr[0]*dx1; 
  outPhr[1] += incr[1]*dx1; 
  outPhr[2] += incr[2]*dx1; 
  outPhr[3] += incr[3]*dx1; 
  outPhr[4] += incr[4]*dx1; 
  outPhr[5] += incr[5]*dx1; 
  outPhr[6] += incr[6]*dx1; 
  outPhr[7] += incr[7]*dx1; 
  outPhr[8] += incr[8]*dx1; 
  outPhr[9] += incr[9]*dx1; 
  outPhr[10] += incr[10]*dx1; 
  outPhr[11] += incr[11]*dx1; 
  outPhr[12] += incr[12]*dx1; 
  outPhr[13] += incr[13]*dx1; 
  outPhr[14] += incr[14]*dx1; 

  outPhl[0] += -1.0*incr[0]*dx1; 
  outPhl[1] += -1.0*incr[1]*dx1; 
  outPhl[2] += incr[2]*dx1; 
  outPhl[3] += incr[3]*dx1; 
  outPhl[4] += -1.0*incr[4]*dx1; 
  outPhl[5] += -1.0*incr[5]*dx1; 
  outPhl[6] += incr[6]*dx1; 
  outPhl[7] += -1.0*incr[7]*dx1; 
  outPhl[8] += -1.0*incr[8]*dx1; 
  outPhl[9] += incr[9]*dx1; 
  outPhl[10] += -1.0*incr[10]*dx1; 
  outPhl[11] += incr[11]*dx1; 
  outPhl[12] += incr[12]*dx1; 
  outPhl[13] += -1.0*incr[13]*dx1; 
  outPhl[14] += -1.0*incr[14]*dx1; 

 
  incr[0] = 0.75*byr[14]*c2*gamma+0.75*byl[14]*c2*gamma-0.6614378277661477*byr[9]*c2*gamma+0.6614378277661477*byl[9]*c2*gamma+0.5590169943749475*byr[5]*c2*gamma+0.5590169943749475*byl[5]*c2*gamma-0.4330127018922193*byr[2]*c2*gamma+0.4330127018922193*byl[2]*c2*gamma+0.25*byr[0]*c2*gamma+0.25*byl[0]*c2*gamma-0.75*psr[14]*c*gamma+0.75*psl[14]*c*gamma+0.6614378277661477*psr[9]*c*gamma+0.6614378277661477*psl[9]*c*gamma-0.5590169943749475*psr[5]*c*gamma+0.5590169943749475*psl[5]*c*gamma+0.4330127018922193*psr[2]*c*gamma+0.4330127018922193*psl[2]*c*gamma-0.25*psr[0]*c*gamma+0.25*psl[0]*c*gamma; 
  incr[1] = (-0.6614378277661477*byr[12]*c2*gamma)+0.6614378277661477*byl[12]*c2*gamma+0.5590169943749476*byr[7]*c2*gamma+0.5590169943749476*byl[7]*c2*gamma-0.4330127018922193*byr[3]*c2*gamma+0.4330127018922193*byl[3]*c2*gamma+0.25*byr[1]*c2*gamma+0.25*byl[1]*c2*gamma+0.6614378277661477*psr[12]*c*gamma+0.6614378277661477*psl[12]*c*gamma-0.5590169943749476*psr[7]*c*gamma+0.5590169943749476*psl[7]*c*gamma+0.4330127018922193*psr[3]*c*gamma+0.4330127018922193*psl[3]*c*gamma-0.25*psr[1]*c*gamma+0.25*psl[1]*c*gamma; 
  incr[2] = (-1.299038105676658*byr[14]*c2*gamma)-1.299038105676658*byl[14]*c2*gamma+1.14564392373896*byr[9]*c2*gamma-1.14564392373896*byl[9]*c2*gamma-0.9682458365518543*byr[5]*c2*gamma-0.9682458365518543*byl[5]*c2*gamma+0.75*byr[2]*c2*gamma-0.75*byl[2]*c2*gamma-0.4330127018922193*byr[0]*c2*gamma-0.4330127018922193*byl[0]*c2*gamma+1.299038105676658*psr[14]*c*gamma-1.299038105676658*psl[14]*c*gamma-1.14564392373896*psr[9]*c*gamma-1.14564392373896*psl[9]*c*gamma+0.9682458365518543*psr[5]*c*gamma-0.9682458365518543*psl[5]*c*gamma-0.75*psr[2]*c*gamma-0.75*psl[2]*c*gamma+0.4330127018922193*psr[0]*c*gamma-0.4330127018922193*psl[0]*c*gamma; 
  incr[3] = 1.14564392373896*byr[12]*c2*gamma-1.14564392373896*byl[12]*c2*gamma-0.9682458365518543*byr[7]*c2*gamma-0.9682458365518543*byl[7]*c2*gamma+0.75*byr[3]*c2*gamma-0.75*byl[3]*c2*gamma-0.4330127018922193*byr[1]*c2*gamma-0.4330127018922193*byl[1]*c2*gamma-1.14564392373896*psr[12]*c*gamma-1.14564392373896*psl[12]*c*gamma+0.9682458365518543*psr[7]*c*gamma-0.9682458365518543*psl[7]*c*gamma-0.75*psr[3]*c*gamma-0.75*psl[3]*c*gamma+0.4330127018922193*psr[1]*c*gamma-0.4330127018922193*psl[1]*c*gamma; 
  incr[4] = 0.5590169943749475*byr[10]*c2*gamma+0.5590169943749475*byl[10]*c2*gamma-0.4330127018922194*byr[6]*c2*gamma+0.4330127018922194*byl[6]*c2*gamma+0.25*byr[4]*c2*gamma+0.25*byl[4]*c2*gamma-0.5590169943749475*psr[10]*c*gamma+0.5590169943749475*psl[10]*c*gamma+0.4330127018922194*psr[6]*c*gamma+0.4330127018922194*psl[6]*c*gamma-0.25*psr[4]*c*gamma+0.25*psl[4]*c*gamma; 
  incr[5] = 1.677050983124842*byr[14]*c2*gamma+1.677050983124842*byl[14]*c2*gamma-1.479019945774904*byr[9]*c2*gamma+1.479019945774904*byl[9]*c2*gamma+1.25*byr[5]*c2*gamma+1.25*byl[5]*c2*gamma-0.9682458365518543*byr[2]*c2*gamma+0.9682458365518543*byl[2]*c2*gamma+0.5590169943749475*byr[0]*c2*gamma+0.5590169943749475*byl[0]*c2*gamma-1.677050983124842*psr[14]*c*gamma+1.677050983124842*psl[14]*c*gamma+1.479019945774904*psr[9]*c*gamma+1.479019945774904*psl[9]*c*gamma-1.25*psr[5]*c*gamma+1.25*psl[5]*c*gamma+0.9682458365518543*psr[2]*c*gamma+0.9682458365518543*psl[2]*c*gamma-0.5590169943749475*psr[0]*c*gamma+0.5590169943749475*psl[0]*c*gamma; 
  incr[6] = (-0.9682458365518543*byr[10]*c2*gamma)-0.9682458365518543*byl[10]*c2*gamma+0.75*byr[6]*c2*gamma-0.75*byl[6]*c2*gamma-0.4330127018922194*byr[4]*c2*gamma-0.4330127018922194*byl[4]*c2*gamma+0.9682458365518543*psr[10]*c*gamma-0.9682458365518543*psl[10]*c*gamma-0.75*psr[6]*c*gamma-0.75*psl[6]*c*gamma+0.4330127018922194*psr[4]*c*gamma-0.4330127018922194*psl[4]*c*gamma; 
  incr[7] = (-1.479019945774904*byr[12]*c2*gamma)+1.479019945774904*byl[12]*c2*gamma+1.25*byr[7]*c2*gamma+1.25*byl[7]*c2*gamma-0.9682458365518543*byr[3]*c2*gamma+0.9682458365518543*byl[3]*c2*gamma+0.5590169943749476*byr[1]*c2*gamma+0.5590169943749476*byl[1]*c2*gamma+1.479019945774904*psr[12]*c*gamma+1.479019945774904*psl[12]*c*gamma-1.25*psr[7]*c*gamma+1.25*psl[7]*c*gamma+0.9682458365518543*psr[3]*c*gamma+0.9682458365518543*psl[3]*c*gamma-0.5590169943749476*psr[1]*c*gamma+0.5590169943749476*psl[1]*c*gamma; 
  incr[8] = (-0.4330127018922193*byr[11]*c2*gamma)+0.4330127018922193*byl[11]*c2*gamma+0.25*byr[8]*c2*gamma+0.25*byl[8]*c2*gamma+0.4330127018922193*psr[11]*c*gamma+0.4330127018922193*psl[11]*c*gamma-0.25*psr[8]*c*gamma+0.25*psl[8]*c*gamma; 
  incr[9] = (-1.984313483298443*byr[14]*c2*gamma)-1.984313483298443*byl[14]*c2*gamma+1.75*byr[9]*c2*gamma-1.75*byl[9]*c2*gamma-1.479019945774904*byr[5]*c2*gamma-1.479019945774904*byl[5]*c2*gamma+1.14564392373896*byr[2]*c2*gamma-1.14564392373896*byl[2]*c2*gamma-0.6614378277661477*byr[0]*c2*gamma-0.6614378277661477*byl[0]*c2*gamma+1.984313483298443*psr[14]*c*gamma-1.984313483298443*psl[14]*c*gamma-1.75*psr[9]*c*gamma-1.75*psl[9]*c*gamma+1.479019945774904*psr[5]*c*gamma-1.479019945774904*psl[5]*c*gamma-1.14564392373896*psr[2]*c*gamma-1.14564392373896*psl[2]*c*gamma+0.6614378277661477*psr[0]*c*gamma-0.6614378277661477*psl[0]*c*gamma; 
  incr[10] = 1.25*byr[10]*c2*gamma+1.25*byl[10]*c2*gamma-0.9682458365518543*byr[6]*c2*gamma+0.9682458365518543*byl[6]*c2*gamma+0.5590169943749475*byr[4]*c2*gamma+0.5590169943749475*byl[4]*c2*gamma-1.25*psr[10]*c*gamma+1.25*psl[10]*c*gamma+0.9682458365518543*psr[6]*c*gamma+0.9682458365518543*psl[6]*c*gamma-0.5590169943749475*psr[4]*c*gamma+0.5590169943749475*psl[4]*c*gamma; 
  incr[11] = 0.75*byr[11]*c2*gamma-0.75*byl[11]*c2*gamma-0.4330127018922193*byr[8]*c2*gamma-0.4330127018922193*byl[8]*c2*gamma-0.75*psr[11]*c*gamma-0.75*psl[11]*c*gamma+0.4330127018922193*psr[8]*c*gamma-0.4330127018922193*psl[8]*c*gamma; 
  incr[12] = 1.75*byr[12]*c2*gamma-1.75*byl[12]*c2*gamma-1.479019945774904*byr[7]*c2*gamma-1.479019945774904*byl[7]*c2*gamma+1.14564392373896*byr[3]*c2*gamma-1.14564392373896*byl[3]*c2*gamma-0.6614378277661477*byr[1]*c2*gamma-0.6614378277661477*byl[1]*c2*gamma-1.75*psr[12]*c*gamma-1.75*psl[12]*c*gamma+1.479019945774904*psr[7]*c*gamma-1.479019945774904*psl[7]*c*gamma-1.14564392373896*psr[3]*c*gamma-1.14564392373896*psl[3]*c*gamma+0.6614378277661477*psr[1]*c*gamma-0.6614378277661477*psl[1]*c*gamma; 
  incr[13] = 0.25*byr[13]*c2*gamma+0.25*byl[13]*c2*gamma-0.25*psr[13]*c*gamma+0.25*psl[13]*c*gamma; 
  incr[14] = 2.25*byr[14]*c2*gamma+2.25*byl[14]*c2*gamma-1.984313483298443*byr[9]*c2*gamma+1.984313483298443*byl[9]*c2*gamma+1.677050983124842*byr[5]*c2*gamma+1.677050983124842*byl[5]*c2*gamma-1.299038105676658*byr[2]*c2*gamma+1.299038105676658*byl[2]*c2*gamma+0.75*byr[0]*c2*gamma+0.75*byl[0]*c2*gamma-2.25*psr[14]*c*gamma+2.25*psl[14]*c*gamma+1.984313483298443*psr[9]*c*gamma+1.984313483298443*psl[9]*c*gamma-1.677050983124842*psr[5]*c*gamma+1.677050983124842*psl[5]*c*gamma+1.299038105676658*psr[2]*c*gamma+1.299038105676658*psl[2]*c*gamma-0.75*psr[0]*c*gamma+0.75*psl[0]*c*gamma; 

  outPsr[0] += incr[0]*dx1; 
  outPsr[1] += incr[1]*dx1; 
  outPsr[2] += incr[2]*dx1; 
  outPsr[3] += incr[3]*dx1; 
  outPsr[4] += incr[4]*dx1; 
  outPsr[5] += incr[5]*dx1; 
  outPsr[6] += incr[6]*dx1; 
  outPsr[7] += incr[7]*dx1; 
  outPsr[8] += incr[8]*dx1; 
  outPsr[9] += incr[9]*dx1; 
  outPsr[10] += incr[10]*dx1; 
  outPsr[11] += incr[11]*dx1; 
  outPsr[12] += incr[12]*dx1; 
  outPsr[13] += incr[13]*dx1; 
  outPsr[14] += incr[14]*dx1; 

  outPsl[0] += -1.0*incr[0]*dx1; 
  outPsl[1] += -1.0*incr[1]*dx1; 
  outPsl[2] += incr[2]*dx1; 
  outPsl[3] += incr[3]*dx1; 
  outPsl[4] += -1.0*incr[4]*dx1; 
  outPsl[5] += -1.0*incr[5]*dx1; 
  outPsl[6] += incr[6]*dx1; 
  outPsl[7] += -1.0*incr[7]*dx1; 
  outPsl[8] += -1.0*incr[8]*dx1; 
  outPsl[9] += incr[9]*dx1; 
  outPsl[10] += -1.0*incr[10]*dx1; 
  outPsl[11] += incr[11]*dx1; 
  outPsl[12] += incr[12]*dx1; 
  outPsl[13] += -1.0*incr[13]*dx1; 
  outPsl[14] += -1.0*incr[14]*dx1; 

 
} 
