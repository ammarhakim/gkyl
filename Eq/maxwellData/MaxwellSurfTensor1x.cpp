#include <MaxwellModDecl.h> 
double MaxwellSurf1xTensor_X_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[0]; 
  const double dxr1 = 2.0/dxr[0]; 
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
 
  incr[0] = (0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+((-0.4330127018922193*phr[1])+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = ((-0.75*(exr[1]+exl[1]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+(0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 

 
  incr[0] = (0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+((-0.4330127018922193*bzr[1])+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = ((-0.75*(eyr[1]+eyl[1]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+(0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 

 
  incr[0] = (0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+(0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = ((-0.75*(ezr[1]+ezl[1]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+((-0.75*byr[1])+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (((-0.75*(bxr[1]+bxl[1]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = (((-0.75*(byr[1]+byl[1]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = (((-0.75*(bzr[1]+bzl[1]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = (((-0.75*(phr[1]+phl[1]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 

 
  incr[0] = (0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*bxr[1])+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = ((-0.75*(psr[1]+psl[1]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 

 
  return std::fmax(c, tau); 
} 
double MaxwellSurf1xTensor_X_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[0]; 
  const double dxr1 = 2.0/dxr[0]; 
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
 
  incr[0] = ((-0.5590169943749475*exr[2])+0.5590169943749475*exl[2]+0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+(0.5590169943749475*(phr[2]+phl[2])-0.4330127018922193*phr[1]+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.9682458365518543*exr[2]-0.9682458365518543*exl[2]-0.75*(exr[1]+exl[1])+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+((-0.9682458365518543*(phr[2]+phl[2]))+0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[2] = ((-1.25*exr[2])+1.25*exl[2]+0.9682458365518543*(exr[1]+exl[1])-0.5590169943749475*exr[0]+0.5590169943749475*exl[0])*c*chi+(1.25*(phr[2]+phl[2])-0.9682458365518543*phr[1]+0.9682458365518543*phl[1]+0.5590169943749475*(phr[0]+phl[0]))*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*eyr[2])+0.5590169943749475*eyl[2]+0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+(0.5590169943749475*(bzr[2]+bzl[2])-0.4330127018922193*bzr[1]+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = (0.9682458365518543*eyr[2]-0.9682458365518543*eyl[2]-0.75*(eyr[1]+eyl[1])+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+((-0.9682458365518543*(bzr[2]+bzl[2]))+0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[2] = ((-1.25*eyr[2])+1.25*eyl[2]+0.9682458365518543*(eyr[1]+eyl[1])-0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0])*tau+(1.25*(bzr[2]+bzl[2])-0.9682458365518543*bzr[1]+0.9682458365518543*bzl[1]+0.5590169943749475*(bzr[0]+bzl[0]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*ezr[2])+0.5590169943749475*ezl[2]+0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+((-0.5590169943749475*(byr[2]+byl[2]))+0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = (0.9682458365518543*ezr[2]-0.9682458365518543*ezl[2]-0.75*(ezr[1]+ezl[1])+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+(0.9682458365518543*(byr[2]+byl[2])-0.75*byr[1]+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[2] = ((-1.25*ezr[2])+1.25*ezl[2]+0.9682458365518543*(ezr[1]+ezl[1])-0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0])*tau+((-1.25*(byr[2]+byl[2]))+0.9682458365518543*byr[1]-0.9682458365518543*byl[1]-0.5590169943749475*(byr[0]+byl[0]))*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*bxr[2])+0.5590169943749475*bxl[2]+0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c+0.5590169943749475*(psr[2]+psl[2])-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = ((0.9682458365518543*bxr[2]-0.9682458365518543*bxl[2]-0.75*(bxr[1]+bxl[1])+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c-0.9682458365518543*(psr[2]+psl[2])+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[2] = (((-1.25*bxr[2])+1.25*bxl[2]+0.9682458365518543*(bxr[1]+bxl[1])-0.5590169943749475*bxr[0]+0.5590169943749475*bxl[0])*c+1.25*(psr[2]+psl[2])-0.9682458365518543*psr[1]+0.9682458365518543*psl[1]+0.5590169943749475*(psr[0]+psl[0]))*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*byr[2])+0.5590169943749475*byl[2]+0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau-0.5590169943749475*(ezr[2]+ezl[2])+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = ((0.9682458365518543*byr[2]-0.9682458365518543*byl[2]-0.75*(byr[1]+byl[1])+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau+0.9682458365518543*(ezr[2]+ezl[2])-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[2] = (((-1.25*byr[2])+1.25*byl[2]+0.9682458365518543*(byr[1]+byl[1])-0.5590169943749475*byr[0]+0.5590169943749475*byl[0])*c2)/tau-1.25*(ezr[2]+ezl[2])+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*(ezr[0]+ezl[0]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*bzr[2])+0.5590169943749475*bzl[2]+0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau+0.5590169943749475*(eyr[2]+eyl[2])-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = ((0.9682458365518543*bzr[2]-0.9682458365518543*bzl[2]-0.75*(bzr[1]+bzl[1])+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau-0.9682458365518543*(eyr[2]+eyl[2])+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[2] = (((-1.25*bzr[2])+1.25*bzl[2]+0.9682458365518543*(bzr[1]+bzl[1])-0.5590169943749475*bzr[0]+0.5590169943749475*bzl[0])*c2)/tau+1.25*(eyr[2]+eyl[2])-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*(eyr[0]+eyl[0]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*phr[2])+0.5590169943749475*phl[2]+0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c+0.5590169943749475*(exr[2]+exl[2])-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = ((0.9682458365518543*phr[2]-0.9682458365518543*phl[2]-0.75*(phr[1]+phl[1])+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c-0.9682458365518543*(exr[2]+exl[2])+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 
  incr[2] = (((-1.25*phr[2])+1.25*phl[2]+0.9682458365518543*(phr[1]+phl[1])-0.5590169943749475*phr[0]+0.5590169943749475*phl[0])*c+1.25*(exr[2]+exl[2])-0.9682458365518543*exr[1]+0.9682458365518543*exl[1]+0.5590169943749475*(exr[0]+exl[0]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*psr[2])+0.5590169943749475*psl[2]+0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+(0.5590169943749475*(bxr[2]+bxl[2])-0.4330127018922193*bxr[1]+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = (0.9682458365518543*psr[2]-0.9682458365518543*psl[2]-0.75*(psr[1]+psl[1])+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+((-0.9682458365518543*(bxr[2]+bxl[2]))+0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 
  incr[2] = ((-1.25*psr[2])+1.25*psl[2]+0.9682458365518543*(psr[1]+psl[1])-0.5590169943749475*psr[0]+0.5590169943749475*psl[0])*c*gamma+(1.25*(bxr[2]+bxl[2])-0.9682458365518543*bxr[1]+0.9682458365518543*bxl[1]+0.5590169943749475*(bxr[0]+bxl[0]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 

 
  return std::fmax(c, tau); 
} 
double MaxwellSurf1xTensor_X_P3(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[0]; 
  const double dxr1 = 2.0/dxr[0]; 
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
 
  incr[0] = (0.6614378277661477*(exr[3]+exl[3])-0.5590169943749475*exr[2]+0.5590169943749475*exl[2]+0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+((-0.6614378277661477*phr[3])+0.6614378277661477*phl[3]+0.5590169943749475*(phr[2]+phl[2])-0.4330127018922193*phr[1]+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = ((-1.14564392373896*(exr[3]+exl[3]))+0.9682458365518543*exr[2]-0.9682458365518543*exl[2]-0.75*(exr[1]+exl[1])+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+(1.14564392373896*phr[3]-1.14564392373896*phl[3]-0.9682458365518543*(phr[2]+phl[2])+0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[2] = (1.479019945774904*(exr[3]+exl[3])-1.25*exr[2]+1.25*exl[2]+0.9682458365518543*(exr[1]+exl[1])-0.5590169943749475*exr[0]+0.5590169943749475*exl[0])*c*chi+((-1.479019945774904*phr[3])+1.479019945774904*phl[3]+1.25*(phr[2]+phl[2])-0.9682458365518543*phr[1]+0.9682458365518543*phl[1]+0.5590169943749475*(phr[0]+phl[0]))*c2chi; 
  incr[3] = ((-1.75*(exr[3]+exl[3]))+1.479019945774904*exr[2]-1.479019945774904*exl[2]-1.14564392373896*(exr[1]+exl[1])+0.6614378277661477*exr[0]-0.6614378277661477*exl[0])*c*chi+(1.75*phr[3]-1.75*phl[3]-1.479019945774904*(phr[2]+phl[2])+1.14564392373896*phr[1]-1.14564392373896*phl[1]-0.6614378277661477*(phr[0]+phl[0]))*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 
  outExl[3] += incr[3]*dxl1; 

 
  incr[0] = (0.6614378277661477*(eyr[3]+eyl[3])-0.5590169943749475*eyr[2]+0.5590169943749475*eyl[2]+0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+((-0.6614378277661477*bzr[3])+0.6614378277661477*bzl[3]+0.5590169943749475*(bzr[2]+bzl[2])-0.4330127018922193*bzr[1]+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = ((-1.14564392373896*(eyr[3]+eyl[3]))+0.9682458365518543*eyr[2]-0.9682458365518543*eyl[2]-0.75*(eyr[1]+eyl[1])+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+(1.14564392373896*bzr[3]-1.14564392373896*bzl[3]-0.9682458365518543*(bzr[2]+bzl[2])+0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[2] = (1.479019945774904*(eyr[3]+eyl[3])-1.25*eyr[2]+1.25*eyl[2]+0.9682458365518543*(eyr[1]+eyl[1])-0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0])*tau+((-1.479019945774904*bzr[3])+1.479019945774904*bzl[3]+1.25*(bzr[2]+bzl[2])-0.9682458365518543*bzr[1]+0.9682458365518543*bzl[1]+0.5590169943749475*(bzr[0]+bzl[0]))*c2; 
  incr[3] = ((-1.75*(eyr[3]+eyl[3]))+1.479019945774904*eyr[2]-1.479019945774904*eyl[2]-1.14564392373896*(eyr[1]+eyl[1])+0.6614378277661477*eyr[0]-0.6614378277661477*eyl[0])*tau+(1.75*bzr[3]-1.75*bzl[3]-1.479019945774904*(bzr[2]+bzl[2])+1.14564392373896*bzr[1]-1.14564392373896*bzl[1]-0.6614378277661477*(bzr[0]+bzl[0]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 
  outEyl[3] += incr[3]*dxl1; 

 
  incr[0] = (0.6614378277661477*(ezr[3]+ezl[3])-0.5590169943749475*ezr[2]+0.5590169943749475*ezl[2]+0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+(0.6614378277661477*byr[3]-0.6614378277661477*byl[3]-0.5590169943749475*(byr[2]+byl[2])+0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = ((-1.14564392373896*(ezr[3]+ezl[3]))+0.9682458365518543*ezr[2]-0.9682458365518543*ezl[2]-0.75*(ezr[1]+ezl[1])+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+((-1.14564392373896*byr[3])+1.14564392373896*byl[3]+0.9682458365518543*(byr[2]+byl[2])-0.75*byr[1]+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[2] = (1.479019945774904*(ezr[3]+ezl[3])-1.25*ezr[2]+1.25*ezl[2]+0.9682458365518543*(ezr[1]+ezl[1])-0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0])*tau+(1.479019945774904*byr[3]-1.479019945774904*byl[3]-1.25*(byr[2]+byl[2])+0.9682458365518543*byr[1]-0.9682458365518543*byl[1]-0.5590169943749475*(byr[0]+byl[0]))*c2; 
  incr[3] = ((-1.75*(ezr[3]+ezl[3]))+1.479019945774904*ezr[2]-1.479019945774904*ezl[2]-1.14564392373896*(ezr[1]+ezl[1])+0.6614378277661477*ezr[0]-0.6614378277661477*ezl[0])*tau+((-1.75*byr[3])+1.75*byl[3]+1.479019945774904*(byr[2]+byl[2])-1.14564392373896*byr[1]+1.14564392373896*byl[1]+0.6614378277661477*(byr[0]+byl[0]))*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 
  outEzl[3] += incr[3]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(bxr[3]+bxl[3])-0.5590169943749475*bxr[2]+0.5590169943749475*bxl[2]+0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c-0.6614378277661477*psr[3]+0.6614378277661477*psl[3]+0.5590169943749475*(psr[2]+psl[2])-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (((-1.14564392373896*(bxr[3]+bxl[3]))+0.9682458365518543*bxr[2]-0.9682458365518543*bxl[2]-0.75*(bxr[1]+bxl[1])+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c+1.14564392373896*psr[3]-1.14564392373896*psl[3]-0.9682458365518543*(psr[2]+psl[2])+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[2] = ((1.479019945774904*(bxr[3]+bxl[3])-1.25*bxr[2]+1.25*bxl[2]+0.9682458365518543*(bxr[1]+bxl[1])-0.5590169943749475*bxr[0]+0.5590169943749475*bxl[0])*c-1.479019945774904*psr[3]+1.479019945774904*psl[3]+1.25*(psr[2]+psl[2])-0.9682458365518543*psr[1]+0.9682458365518543*psl[1]+0.5590169943749475*(psr[0]+psl[0]))*gamma; 
  incr[3] = (((-1.75*(bxr[3]+bxl[3]))+1.479019945774904*bxr[2]-1.479019945774904*bxl[2]-1.14564392373896*(bxr[1]+bxl[1])+0.6614378277661477*bxr[0]-0.6614378277661477*bxl[0])*c+1.75*psr[3]-1.75*psl[3]-1.479019945774904*(psr[2]+psl[2])+1.14564392373896*psr[1]-1.14564392373896*psl[1]-0.6614378277661477*(psr[0]+psl[0]))*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 
  outBxl[3] += incr[3]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(byr[3]+byl[3])-0.5590169943749475*byr[2]+0.5590169943749475*byl[2]+0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau+0.6614378277661477*ezr[3]-0.6614378277661477*ezl[3]-0.5590169943749475*(ezr[2]+ezl[2])+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = (((-1.14564392373896*(byr[3]+byl[3]))+0.9682458365518543*byr[2]-0.9682458365518543*byl[2]-0.75*(byr[1]+byl[1])+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau-1.14564392373896*ezr[3]+1.14564392373896*ezl[3]+0.9682458365518543*(ezr[2]+ezl[2])-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[2] = ((1.479019945774904*(byr[3]+byl[3])-1.25*byr[2]+1.25*byl[2]+0.9682458365518543*(byr[1]+byl[1])-0.5590169943749475*byr[0]+0.5590169943749475*byl[0])*c2)/tau+1.479019945774904*ezr[3]-1.479019945774904*ezl[3]-1.25*(ezr[2]+ezl[2])+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*(ezr[0]+ezl[0]); 
  incr[3] = (((-1.75*(byr[3]+byl[3]))+1.479019945774904*byr[2]-1.479019945774904*byl[2]-1.14564392373896*(byr[1]+byl[1])+0.6614378277661477*byr[0]-0.6614378277661477*byl[0])*c2)/tau-1.75*ezr[3]+1.75*ezl[3]+1.479019945774904*(ezr[2]+ezl[2])-1.14564392373896*ezr[1]+1.14564392373896*ezl[1]+0.6614378277661477*(ezr[0]+ezl[0]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 
  outByl[3] += incr[3]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(bzr[3]+bzl[3])-0.5590169943749475*bzr[2]+0.5590169943749475*bzl[2]+0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau-0.6614378277661477*eyr[3]+0.6614378277661477*eyl[3]+0.5590169943749475*(eyr[2]+eyl[2])-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = (((-1.14564392373896*(bzr[3]+bzl[3]))+0.9682458365518543*bzr[2]-0.9682458365518543*bzl[2]-0.75*(bzr[1]+bzl[1])+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau+1.14564392373896*eyr[3]-1.14564392373896*eyl[3]-0.9682458365518543*(eyr[2]+eyl[2])+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[2] = ((1.479019945774904*(bzr[3]+bzl[3])-1.25*bzr[2]+1.25*bzl[2]+0.9682458365518543*(bzr[1]+bzl[1])-0.5590169943749475*bzr[0]+0.5590169943749475*bzl[0])*c2)/tau-1.479019945774904*eyr[3]+1.479019945774904*eyl[3]+1.25*(eyr[2]+eyl[2])-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*(eyr[0]+eyl[0]); 
  incr[3] = (((-1.75*(bzr[3]+bzl[3]))+1.479019945774904*bzr[2]-1.479019945774904*bzl[2]-1.14564392373896*(bzr[1]+bzl[1])+0.6614378277661477*bzr[0]-0.6614378277661477*bzl[0])*c2)/tau+1.75*eyr[3]-1.75*eyl[3]-1.479019945774904*(eyr[2]+eyl[2])+1.14564392373896*eyr[1]-1.14564392373896*eyl[1]-0.6614378277661477*(eyr[0]+eyl[0]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 
  outBzl[3] += incr[3]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(phr[3]+phl[3])-0.5590169943749475*phr[2]+0.5590169943749475*phl[2]+0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c-0.6614378277661477*exr[3]+0.6614378277661477*exl[3]+0.5590169943749475*(exr[2]+exl[2])-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = (((-1.14564392373896*(phr[3]+phl[3]))+0.9682458365518543*phr[2]-0.9682458365518543*phl[2]-0.75*(phr[1]+phl[1])+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+1.14564392373896*exr[3]-1.14564392373896*exl[3]-0.9682458365518543*(exr[2]+exl[2])+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 
  incr[2] = ((1.479019945774904*(phr[3]+phl[3])-1.25*phr[2]+1.25*phl[2]+0.9682458365518543*(phr[1]+phl[1])-0.5590169943749475*phr[0]+0.5590169943749475*phl[0])*c-1.479019945774904*exr[3]+1.479019945774904*exl[3]+1.25*(exr[2]+exl[2])-0.9682458365518543*exr[1]+0.9682458365518543*exl[1]+0.5590169943749475*(exr[0]+exl[0]))*chi; 
  incr[3] = (((-1.75*(phr[3]+phl[3]))+1.479019945774904*phr[2]-1.479019945774904*phl[2]-1.14564392373896*(phr[1]+phl[1])+0.6614378277661477*phr[0]-0.6614378277661477*phl[0])*c+1.75*exr[3]-1.75*exl[3]-1.479019945774904*(exr[2]+exl[2])+1.14564392373896*exr[1]-1.14564392373896*exl[1]-0.6614378277661477*(exr[0]+exl[0]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 
  outPhl[3] += incr[3]*dxl1; 

 
  incr[0] = (0.6614378277661477*(psr[3]+psl[3])-0.5590169943749475*psr[2]+0.5590169943749475*psl[2]+0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.6614378277661477*bxr[3])+0.6614378277661477*bxl[3]+0.5590169943749475*(bxr[2]+bxl[2])-0.4330127018922193*bxr[1]+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = ((-1.14564392373896*(psr[3]+psl[3]))+0.9682458365518543*psr[2]-0.9682458365518543*psl[2]-0.75*(psr[1]+psl[1])+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(1.14564392373896*bxr[3]-1.14564392373896*bxl[3]-0.9682458365518543*(bxr[2]+bxl[2])+0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 
  incr[2] = (1.479019945774904*(psr[3]+psl[3])-1.25*psr[2]+1.25*psl[2]+0.9682458365518543*(psr[1]+psl[1])-0.5590169943749475*psr[0]+0.5590169943749475*psl[0])*c*gamma+((-1.479019945774904*bxr[3])+1.479019945774904*bxl[3]+1.25*(bxr[2]+bxl[2])-0.9682458365518543*bxr[1]+0.9682458365518543*bxl[1]+0.5590169943749475*(bxr[0]+bxl[0]))*c2gamma; 
  incr[3] = ((-1.75*(psr[3]+psl[3]))+1.479019945774904*psr[2]-1.479019945774904*psl[2]-1.14564392373896*(psr[1]+psl[1])+0.6614378277661477*psr[0]-0.6614378277661477*psl[0])*c*gamma+(1.75*bxr[3]-1.75*bxl[3]-1.479019945774904*(bxr[2]+bxl[2])+1.14564392373896*bxr[1]-1.14564392373896*bxl[1]-0.6614378277661477*(bxr[0]+bxl[0]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 
  outPsl[3] += incr[3]*dxl1; 

 
  return std::fmax(c, tau); 
} 
