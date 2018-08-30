#include <MaxwellModDecl.h> 
double MaxwellSurf3xMax_X_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
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
 
  incr[0] = (0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+((-0.4330127018922193*phr[1])+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = ((-0.75*(exr[1]+exl[1]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+(0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[2] = (0.25*exl[2]-0.25*exr[2])*c*chi+0.25*(phr[2]+phl[2])*c2chi; 
  incr[3] = (0.25*exl[3]-0.25*exr[3])*c*chi+0.25*(phr[3]+phl[3])*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 
  outExl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = (0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+((-0.4330127018922193*bzr[1])+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = ((-0.75*(eyr[1]+eyl[1]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+(0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[2] = 0.25*(eyl[2]*tau+(bzr[2]+bzl[2])*c2)-0.25*eyr[2]*tau; 
  incr[3] = 0.25*(eyl[3]*tau+(bzr[3]+bzl[3])*c2)-0.25*eyr[3]*tau; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 
  outEyl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = (0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+(0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = ((-0.75*(ezr[1]+ezl[1]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+((-0.75*byr[1])+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[2] = (0.25*ezl[2]-0.25*ezr[2])*tau-0.25*(byr[2]+byl[2])*c2; 
  incr[3] = (0.25*ezl[3]-0.25*ezr[3])*tau-0.25*(byr[3]+byl[3])*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 
  outEzl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (((-0.75*(bxr[1]+bxl[1]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[2] = (0.25*(bxl[2]*c+psr[2]+psl[2])-0.25*bxr[2]*c)*gamma; 
  incr[3] = (0.25*(bxl[3]*c+psr[3]+psl[3])-0.25*bxr[3]*c)*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 
  outBxl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = (((-0.75*(byr[1]+byl[1]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[2] = ((0.25*byl[2]-0.25*byr[2])*c2)/tau-0.25*(ezr[2]+ezl[2]); 
  incr[3] = ((0.25*byl[3]-0.25*byr[3])*c2)/tau-0.25*(ezr[3]+ezl[3]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 
  outByl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = (((-0.75*(bzr[1]+bzl[1]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[2] = ((0.25*bzl[2]-0.25*bzr[2])*c2)/tau+0.25*(eyr[2]+eyl[2]); 
  incr[3] = ((0.25*bzl[3]-0.25*bzr[3])*c2)/tau+0.25*(eyr[3]+eyl[3]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 
  outBzl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = (((-0.75*(phr[1]+phl[1]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 
  incr[2] = (0.25*(phl[2]*c+exr[2]+exl[2])-0.25*phr[2]*c)*chi; 
  incr[3] = (0.25*(phl[3]*c+exr[3]+exl[3])-0.25*phr[3]*c)*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 
  outPhl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = (0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*bxr[1])+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = ((-0.75*(psr[1]+psl[1]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 
  incr[2] = (0.25*psl[2]-0.25*psr[2])*c*gamma+0.25*(bxr[2]+bxl[2])*c2gamma; 
  incr[3] = (0.25*psl[3]-0.25*psr[3])*c*gamma+0.25*(bxr[3]+bxl[3])*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 
  outPsl[3] += -1.0*incr[3]*dxl1; 

 
  return c; 
} 
double MaxwellSurf3xMax_X_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[0]; 
  const double dxr1 = 2.0/dxr[0]; 
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
 
  incr[0] = ((-0.5590169943749475*exr[7])+0.5590169943749475*exl[7]+0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+(0.5590169943749475*(phr[7]+phl[7])-0.4330127018922193*phr[1]+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.9682458365518543*exr[7]-0.9682458365518543*exl[7]-0.75*(exr[1]+exl[1])+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+((-0.9682458365518543*(phr[7]+phl[7]))+0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[2] = (0.4330127018922193*(exr[4]+exl[4])-0.25*exr[2]+0.25*exl[2])*c*chi+((-0.4330127018922193*phr[4])+0.4330127018922193*phl[4]+0.25*(phr[2]+phl[2]))*c2chi; 
  incr[3] = (0.4330127018922193*(exr[5]+exl[5])-0.25*exr[3]+0.25*exl[3])*c*chi+((-0.4330127018922193*phr[5])+0.4330127018922193*phl[5]+0.25*(phr[3]+phl[3]))*c2chi; 
  incr[4] = ((-0.75*(exr[4]+exl[4]))+0.4330127018922193*exr[2]-0.4330127018922193*exl[2])*c*chi+(0.75*phr[4]-0.75*phl[4]-0.4330127018922193*(phr[2]+phl[2]))*c2chi; 
  incr[5] = ((-0.75*(exr[5]+exl[5]))+0.4330127018922193*exr[3]-0.4330127018922193*exl[3])*c*chi+(0.75*phr[5]-0.75*phl[5]-0.4330127018922193*(phr[3]+phl[3]))*c2chi; 
  incr[6] = (0.25*exl[6]-0.25*exr[6])*c*chi+0.25*(phr[6]+phl[6])*c2chi; 
  incr[7] = ((-1.25*exr[7])+1.25*exl[7]+0.9682458365518543*(exr[1]+exl[1])-0.5590169943749475*exr[0]+0.5590169943749475*exl[0])*c*chi+(1.25*(phr[7]+phl[7])-0.9682458365518543*phr[1]+0.9682458365518543*phl[1]+0.5590169943749475*(phr[0]+phl[0]))*c2chi; 
  incr[8] = (0.25*exl[8]-0.25*exr[8])*c*chi+0.25*(phr[8]+phl[8])*c2chi; 
  incr[9] = (0.25*exl[9]-0.25*exr[9])*c*chi+0.25*(phr[9]+phl[9])*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 
  outExr[4] += incr[4]*dxr1; 
  outExr[5] += incr[5]*dxr1; 
  outExr[6] += incr[6]*dxr1; 
  outExr[7] += incr[7]*dxr1; 
  outExr[8] += incr[8]*dxr1; 
  outExr[9] += incr[9]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 
  outExl[3] += -1.0*incr[3]*dxl1; 
  outExl[4] += incr[4]*dxl1; 
  outExl[5] += incr[5]*dxl1; 
  outExl[6] += -1.0*incr[6]*dxl1; 
  outExl[7] += -1.0*incr[7]*dxl1; 
  outExl[8] += -1.0*incr[8]*dxl1; 
  outExl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*eyr[7])+0.5590169943749475*eyl[7]+0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+(0.5590169943749475*(bzr[7]+bzl[7])-0.4330127018922193*bzr[1]+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = (0.9682458365518543*eyr[7]-0.9682458365518543*eyl[7]-0.75*(eyr[1]+eyl[1])+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+((-0.9682458365518543*(bzr[7]+bzl[7]))+0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[2] = (0.4330127018922193*(eyr[4]+eyl[4])-0.25*eyr[2]+0.25*eyl[2])*tau+((-0.4330127018922193*bzr[4])+0.4330127018922193*bzl[4]+0.25*(bzr[2]+bzl[2]))*c2; 
  incr[3] = (0.4330127018922193*(eyr[5]+eyl[5])-0.25*eyr[3]+0.25*eyl[3])*tau+((-0.4330127018922193*bzr[5])+0.4330127018922193*bzl[5]+0.25*(bzr[3]+bzl[3]))*c2; 
  incr[4] = ((-0.75*(eyr[4]+eyl[4]))+0.4330127018922193*eyr[2]-0.4330127018922193*eyl[2])*tau+(0.75*bzr[4]-0.75*bzl[4]-0.4330127018922193*(bzr[2]+bzl[2]))*c2; 
  incr[5] = ((-0.75*(eyr[5]+eyl[5]))+0.4330127018922193*eyr[3]-0.4330127018922193*eyl[3])*tau+(0.75*bzr[5]-0.75*bzl[5]-0.4330127018922193*(bzr[3]+bzl[3]))*c2; 
  incr[6] = 0.25*(eyl[6]*tau+(bzr[6]+bzl[6])*c2)-0.25*eyr[6]*tau; 
  incr[7] = ((-1.25*eyr[7])+1.25*eyl[7]+0.9682458365518543*(eyr[1]+eyl[1])-0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0])*tau+(1.25*(bzr[7]+bzl[7])-0.9682458365518543*bzr[1]+0.9682458365518543*bzl[1]+0.5590169943749475*(bzr[0]+bzl[0]))*c2; 
  incr[8] = 0.25*(eyl[8]*tau+(bzr[8]+bzl[8])*c2)-0.25*eyr[8]*tau; 
  incr[9] = 0.25*(eyl[9]*tau+(bzr[9]+bzl[9])*c2)-0.25*eyr[9]*tau; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 
  outEyr[4] += incr[4]*dxr1; 
  outEyr[5] += incr[5]*dxr1; 
  outEyr[6] += incr[6]*dxr1; 
  outEyr[7] += incr[7]*dxr1; 
  outEyr[8] += incr[8]*dxr1; 
  outEyr[9] += incr[9]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 
  outEyl[3] += -1.0*incr[3]*dxl1; 
  outEyl[4] += incr[4]*dxl1; 
  outEyl[5] += incr[5]*dxl1; 
  outEyl[6] += -1.0*incr[6]*dxl1; 
  outEyl[7] += -1.0*incr[7]*dxl1; 
  outEyl[8] += -1.0*incr[8]*dxl1; 
  outEyl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*ezr[7])+0.5590169943749475*ezl[7]+0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+((-0.5590169943749475*(byr[7]+byl[7]))+0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = (0.9682458365518543*ezr[7]-0.9682458365518543*ezl[7]-0.75*(ezr[1]+ezl[1])+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+(0.9682458365518543*(byr[7]+byl[7])-0.75*byr[1]+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[2] = (0.4330127018922193*(ezr[4]+ezl[4])-0.25*ezr[2]+0.25*ezl[2])*tau+(0.4330127018922193*byr[4]-0.4330127018922193*byl[4]-0.25*(byr[2]+byl[2]))*c2; 
  incr[3] = (0.4330127018922193*(ezr[5]+ezl[5])-0.25*ezr[3]+0.25*ezl[3])*tau+(0.4330127018922193*byr[5]-0.4330127018922193*byl[5]-0.25*(byr[3]+byl[3]))*c2; 
  incr[4] = ((-0.75*(ezr[4]+ezl[4]))+0.4330127018922193*ezr[2]-0.4330127018922193*ezl[2])*tau+((-0.75*byr[4])+0.75*byl[4]+0.4330127018922193*(byr[2]+byl[2]))*c2; 
  incr[5] = ((-0.75*(ezr[5]+ezl[5]))+0.4330127018922193*ezr[3]-0.4330127018922193*ezl[3])*tau+((-0.75*byr[5])+0.75*byl[5]+0.4330127018922193*(byr[3]+byl[3]))*c2; 
  incr[6] = (0.25*ezl[6]-0.25*ezr[6])*tau-0.25*(byr[6]+byl[6])*c2; 
  incr[7] = ((-1.25*ezr[7])+1.25*ezl[7]+0.9682458365518543*(ezr[1]+ezl[1])-0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0])*tau+((-1.25*(byr[7]+byl[7]))+0.9682458365518543*byr[1]-0.9682458365518543*byl[1]-0.5590169943749475*(byr[0]+byl[0]))*c2; 
  incr[8] = (0.25*ezl[8]-0.25*ezr[8])*tau-0.25*(byr[8]+byl[8])*c2; 
  incr[9] = (0.25*ezl[9]-0.25*ezr[9])*tau-0.25*(byr[9]+byl[9])*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 
  outEzr[4] += incr[4]*dxr1; 
  outEzr[5] += incr[5]*dxr1; 
  outEzr[6] += incr[6]*dxr1; 
  outEzr[7] += incr[7]*dxr1; 
  outEzr[8] += incr[8]*dxr1; 
  outEzr[9] += incr[9]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 
  outEzl[3] += -1.0*incr[3]*dxl1; 
  outEzl[4] += incr[4]*dxl1; 
  outEzl[5] += incr[5]*dxl1; 
  outEzl[6] += -1.0*incr[6]*dxl1; 
  outEzl[7] += -1.0*incr[7]*dxl1; 
  outEzl[8] += -1.0*incr[8]*dxl1; 
  outEzl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*bxr[7])+0.5590169943749475*bxl[7]+0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c+0.5590169943749475*(psr[7]+psl[7])-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = ((0.9682458365518543*bxr[7]-0.9682458365518543*bxl[7]-0.75*(bxr[1]+bxl[1])+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c-0.9682458365518543*(psr[7]+psl[7])+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[2] = ((0.4330127018922193*(bxr[4]+bxl[4])-0.25*bxr[2]+0.25*bxl[2])*c-0.4330127018922193*psr[4]+0.4330127018922193*psl[4]+0.25*(psr[2]+psl[2]))*gamma; 
  incr[3] = ((0.4330127018922193*(bxr[5]+bxl[5])-0.25*bxr[3]+0.25*bxl[3])*c-0.4330127018922193*psr[5]+0.4330127018922193*psl[5]+0.25*(psr[3]+psl[3]))*gamma; 
  incr[4] = (((-0.75*(bxr[4]+bxl[4]))+0.4330127018922193*bxr[2]-0.4330127018922193*bxl[2])*c+0.75*psr[4]-0.75*psl[4]-0.4330127018922193*(psr[2]+psl[2]))*gamma; 
  incr[5] = (((-0.75*(bxr[5]+bxl[5]))+0.4330127018922193*bxr[3]-0.4330127018922193*bxl[3])*c+0.75*psr[5]-0.75*psl[5]-0.4330127018922193*(psr[3]+psl[3]))*gamma; 
  incr[6] = (0.25*(bxl[6]*c+psr[6]+psl[6])-0.25*bxr[6]*c)*gamma; 
  incr[7] = (((-1.25*bxr[7])+1.25*bxl[7]+0.9682458365518543*(bxr[1]+bxl[1])-0.5590169943749475*bxr[0]+0.5590169943749475*bxl[0])*c+1.25*(psr[7]+psl[7])-0.9682458365518543*psr[1]+0.9682458365518543*psl[1]+0.5590169943749475*(psr[0]+psl[0]))*gamma; 
  incr[8] = (0.25*(bxl[8]*c+psr[8]+psl[8])-0.25*bxr[8]*c)*gamma; 
  incr[9] = (0.25*(bxl[9]*c+psr[9]+psl[9])-0.25*bxr[9]*c)*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 
  outBxr[4] += incr[4]*dxr1; 
  outBxr[5] += incr[5]*dxr1; 
  outBxr[6] += incr[6]*dxr1; 
  outBxr[7] += incr[7]*dxr1; 
  outBxr[8] += incr[8]*dxr1; 
  outBxr[9] += incr[9]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 
  outBxl[3] += -1.0*incr[3]*dxl1; 
  outBxl[4] += incr[4]*dxl1; 
  outBxl[5] += incr[5]*dxl1; 
  outBxl[6] += -1.0*incr[6]*dxl1; 
  outBxl[7] += -1.0*incr[7]*dxl1; 
  outBxl[8] += -1.0*incr[8]*dxl1; 
  outBxl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*byr[7])+0.5590169943749475*byl[7]+0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau-0.5590169943749475*(ezr[7]+ezl[7])+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = ((0.9682458365518543*byr[7]-0.9682458365518543*byl[7]-0.75*(byr[1]+byl[1])+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau+0.9682458365518543*(ezr[7]+ezl[7])-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[2] = ((0.4330127018922193*(byr[4]+byl[4])-0.25*byr[2]+0.25*byl[2])*c2)/tau+0.4330127018922193*ezr[4]-0.4330127018922193*ezl[4]-0.25*(ezr[2]+ezl[2]); 
  incr[3] = ((0.4330127018922193*(byr[5]+byl[5])-0.25*byr[3]+0.25*byl[3])*c2)/tau+0.4330127018922193*ezr[5]-0.4330127018922193*ezl[5]-0.25*(ezr[3]+ezl[3]); 
  incr[4] = (((-0.75*(byr[4]+byl[4]))+0.4330127018922193*byr[2]-0.4330127018922193*byl[2])*c2)/tau-0.75*ezr[4]+0.75*ezl[4]+0.4330127018922193*(ezr[2]+ezl[2]); 
  incr[5] = (((-0.75*(byr[5]+byl[5]))+0.4330127018922193*byr[3]-0.4330127018922193*byl[3])*c2)/tau-0.75*ezr[5]+0.75*ezl[5]+0.4330127018922193*(ezr[3]+ezl[3]); 
  incr[6] = ((0.25*byl[6]-0.25*byr[6])*c2)/tau-0.25*(ezr[6]+ezl[6]); 
  incr[7] = (((-1.25*byr[7])+1.25*byl[7]+0.9682458365518543*(byr[1]+byl[1])-0.5590169943749475*byr[0]+0.5590169943749475*byl[0])*c2)/tau-1.25*(ezr[7]+ezl[7])+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*(ezr[0]+ezl[0]); 
  incr[8] = ((0.25*byl[8]-0.25*byr[8])*c2)/tau-0.25*(ezr[8]+ezl[8]); 
  incr[9] = ((0.25*byl[9]-0.25*byr[9])*c2)/tau-0.25*(ezr[9]+ezl[9]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 
  outByr[4] += incr[4]*dxr1; 
  outByr[5] += incr[5]*dxr1; 
  outByr[6] += incr[6]*dxr1; 
  outByr[7] += incr[7]*dxr1; 
  outByr[8] += incr[8]*dxr1; 
  outByr[9] += incr[9]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 
  outByl[3] += -1.0*incr[3]*dxl1; 
  outByl[4] += incr[4]*dxl1; 
  outByl[5] += incr[5]*dxl1; 
  outByl[6] += -1.0*incr[6]*dxl1; 
  outByl[7] += -1.0*incr[7]*dxl1; 
  outByl[8] += -1.0*incr[8]*dxl1; 
  outByl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*bzr[7])+0.5590169943749475*bzl[7]+0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau+0.5590169943749475*(eyr[7]+eyl[7])-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = ((0.9682458365518543*bzr[7]-0.9682458365518543*bzl[7]-0.75*(bzr[1]+bzl[1])+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau-0.9682458365518543*(eyr[7]+eyl[7])+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[2] = ((0.4330127018922193*(bzr[4]+bzl[4])-0.25*bzr[2]+0.25*bzl[2])*c2)/tau-0.4330127018922193*eyr[4]+0.4330127018922193*eyl[4]+0.25*(eyr[2]+eyl[2]); 
  incr[3] = ((0.4330127018922193*(bzr[5]+bzl[5])-0.25*bzr[3]+0.25*bzl[3])*c2)/tau-0.4330127018922193*eyr[5]+0.4330127018922193*eyl[5]+0.25*(eyr[3]+eyl[3]); 
  incr[4] = (((-0.75*(bzr[4]+bzl[4]))+0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2])*c2)/tau+0.75*eyr[4]-0.75*eyl[4]-0.4330127018922193*(eyr[2]+eyl[2]); 
  incr[5] = (((-0.75*(bzr[5]+bzl[5]))+0.4330127018922193*bzr[3]-0.4330127018922193*bzl[3])*c2)/tau+0.75*eyr[5]-0.75*eyl[5]-0.4330127018922193*(eyr[3]+eyl[3]); 
  incr[6] = ((0.25*bzl[6]-0.25*bzr[6])*c2)/tau+0.25*(eyr[6]+eyl[6]); 
  incr[7] = (((-1.25*bzr[7])+1.25*bzl[7]+0.9682458365518543*(bzr[1]+bzl[1])-0.5590169943749475*bzr[0]+0.5590169943749475*bzl[0])*c2)/tau+1.25*(eyr[7]+eyl[7])-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*(eyr[0]+eyl[0]); 
  incr[8] = ((0.25*bzl[8]-0.25*bzr[8])*c2)/tau+0.25*(eyr[8]+eyl[8]); 
  incr[9] = ((0.25*bzl[9]-0.25*bzr[9])*c2)/tau+0.25*(eyr[9]+eyl[9]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 
  outBzr[4] += incr[4]*dxr1; 
  outBzr[5] += incr[5]*dxr1; 
  outBzr[6] += incr[6]*dxr1; 
  outBzr[7] += incr[7]*dxr1; 
  outBzr[8] += incr[8]*dxr1; 
  outBzr[9] += incr[9]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 
  outBzl[3] += -1.0*incr[3]*dxl1; 
  outBzl[4] += incr[4]*dxl1; 
  outBzl[5] += incr[5]*dxl1; 
  outBzl[6] += -1.0*incr[6]*dxl1; 
  outBzl[7] += -1.0*incr[7]*dxl1; 
  outBzl[8] += -1.0*incr[8]*dxl1; 
  outBzl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*phr[7])+0.5590169943749475*phl[7]+0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c+0.5590169943749475*(exr[7]+exl[7])-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = ((0.9682458365518543*phr[7]-0.9682458365518543*phl[7]-0.75*(phr[1]+phl[1])+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c-0.9682458365518543*(exr[7]+exl[7])+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 
  incr[2] = ((0.4330127018922193*(phr[4]+phl[4])-0.25*phr[2]+0.25*phl[2])*c-0.4330127018922193*exr[4]+0.4330127018922193*exl[4]+0.25*(exr[2]+exl[2]))*chi; 
  incr[3] = ((0.4330127018922193*(phr[5]+phl[5])-0.25*phr[3]+0.25*phl[3])*c-0.4330127018922193*exr[5]+0.4330127018922193*exl[5]+0.25*(exr[3]+exl[3]))*chi; 
  incr[4] = (((-0.75*(phr[4]+phl[4]))+0.4330127018922193*phr[2]-0.4330127018922193*phl[2])*c+0.75*exr[4]-0.75*exl[4]-0.4330127018922193*(exr[2]+exl[2]))*chi; 
  incr[5] = (((-0.75*(phr[5]+phl[5]))+0.4330127018922193*phr[3]-0.4330127018922193*phl[3])*c+0.75*exr[5]-0.75*exl[5]-0.4330127018922193*(exr[3]+exl[3]))*chi; 
  incr[6] = (0.25*(phl[6]*c+exr[6]+exl[6])-0.25*phr[6]*c)*chi; 
  incr[7] = (((-1.25*phr[7])+1.25*phl[7]+0.9682458365518543*(phr[1]+phl[1])-0.5590169943749475*phr[0]+0.5590169943749475*phl[0])*c+1.25*(exr[7]+exl[7])-0.9682458365518543*exr[1]+0.9682458365518543*exl[1]+0.5590169943749475*(exr[0]+exl[0]))*chi; 
  incr[8] = (0.25*(phl[8]*c+exr[8]+exl[8])-0.25*phr[8]*c)*chi; 
  incr[9] = (0.25*(phl[9]*c+exr[9]+exl[9])-0.25*phr[9]*c)*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 
  outPhr[4] += incr[4]*dxr1; 
  outPhr[5] += incr[5]*dxr1; 
  outPhr[6] += incr[6]*dxr1; 
  outPhr[7] += incr[7]*dxr1; 
  outPhr[8] += incr[8]*dxr1; 
  outPhr[9] += incr[9]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 
  outPhl[3] += -1.0*incr[3]*dxl1; 
  outPhl[4] += incr[4]*dxl1; 
  outPhl[5] += incr[5]*dxl1; 
  outPhl[6] += -1.0*incr[6]*dxl1; 
  outPhl[7] += -1.0*incr[7]*dxl1; 
  outPhl[8] += -1.0*incr[8]*dxl1; 
  outPhl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*psr[7])+0.5590169943749475*psl[7]+0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+(0.5590169943749475*(bxr[7]+bxl[7])-0.4330127018922193*bxr[1]+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = (0.9682458365518543*psr[7]-0.9682458365518543*psl[7]-0.75*(psr[1]+psl[1])+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+((-0.9682458365518543*(bxr[7]+bxl[7]))+0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 
  incr[2] = (0.4330127018922193*(psr[4]+psl[4])-0.25*psr[2]+0.25*psl[2])*c*gamma+((-0.4330127018922193*bxr[4])+0.4330127018922193*bxl[4]+0.25*(bxr[2]+bxl[2]))*c2gamma; 
  incr[3] = (0.4330127018922193*(psr[5]+psl[5])-0.25*psr[3]+0.25*psl[3])*c*gamma+((-0.4330127018922193*bxr[5])+0.4330127018922193*bxl[5]+0.25*(bxr[3]+bxl[3]))*c2gamma; 
  incr[4] = ((-0.75*(psr[4]+psl[4]))+0.4330127018922193*psr[2]-0.4330127018922193*psl[2])*c*gamma+(0.75*bxr[4]-0.75*bxl[4]-0.4330127018922193*(bxr[2]+bxl[2]))*c2gamma; 
  incr[5] = ((-0.75*(psr[5]+psl[5]))+0.4330127018922193*psr[3]-0.4330127018922193*psl[3])*c*gamma+(0.75*bxr[5]-0.75*bxl[5]-0.4330127018922193*(bxr[3]+bxl[3]))*c2gamma; 
  incr[6] = (0.25*psl[6]-0.25*psr[6])*c*gamma+0.25*(bxr[6]+bxl[6])*c2gamma; 
  incr[7] = ((-1.25*psr[7])+1.25*psl[7]+0.9682458365518543*(psr[1]+psl[1])-0.5590169943749475*psr[0]+0.5590169943749475*psl[0])*c*gamma+(1.25*(bxr[7]+bxl[7])-0.9682458365518543*bxr[1]+0.9682458365518543*bxl[1]+0.5590169943749475*(bxr[0]+bxl[0]))*c2gamma; 
  incr[8] = (0.25*psl[8]-0.25*psr[8])*c*gamma+0.25*(bxr[8]+bxl[8])*c2gamma; 
  incr[9] = (0.25*psl[9]-0.25*psr[9])*c*gamma+0.25*(bxr[9]+bxl[9])*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 
  outPsr[4] += incr[4]*dxr1; 
  outPsr[5] += incr[5]*dxr1; 
  outPsr[6] += incr[6]*dxr1; 
  outPsr[7] += incr[7]*dxr1; 
  outPsr[8] += incr[8]*dxr1; 
  outPsr[9] += incr[9]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 
  outPsl[3] += -1.0*incr[3]*dxl1; 
  outPsl[4] += incr[4]*dxl1; 
  outPsl[5] += incr[5]*dxl1; 
  outPsl[6] += -1.0*incr[6]*dxl1; 
  outPsl[7] += -1.0*incr[7]*dxl1; 
  outPsl[8] += -1.0*incr[8]*dxl1; 
  outPsl[9] += -1.0*incr[9]*dxl1; 

 
  return c; 
} 
double MaxwellSurf3xMax_X_P3(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[0]; 
  const double dxr1 = 2.0/dxr[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[20]; 
  const double *ezl = &ql[40]; 
  const double *bxl = &ql[60]; 
  const double *byl = &ql[80]; 
  const double *bzl = &ql[100]; 
  const double *phl = &ql[120]; 
  const double *psl = &ql[140]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[20]; 
  double *outEzl = &outl[40]; 
  double *outBxl = &outl[60]; 
  double *outByl = &outl[80]; 
  double *outBzl = &outl[100]; 
  double *outPhl = &outl[120]; 
  double *outPsl = &outl[140]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[20]; 
  const double *ezr = &qr[40]; 
  const double *bxr = &qr[60]; 
  const double *byr = &qr[80]; 
  const double *bzr = &qr[100]; 
  const double *phr = &qr[120]; 
  const double *psr = &qr[140]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[20]; 
  double *outEzr = &outr[40]; 
  double *outBxr = &outr[60]; 
  double *outByr = &outr[80]; 
  double *outBzr = &outr[100]; 
  double *outPhr = &outr[120]; 
  double *outPsr = &outr[140]; 
 
  double incr[20]; 
 
  incr[0] = (0.6614378277661477*(exr[17]+exl[17])-0.5590169943749475*exr[7]+0.5590169943749475*exl[7]+0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+((-0.6614378277661477*phr[17])+0.6614378277661477*phl[17]+0.5590169943749475*(phr[7]+phl[7])-0.4330127018922193*phr[1]+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = ((-1.14564392373896*(exr[17]+exl[17]))+0.9682458365518543*exr[7]-0.9682458365518543*exl[7]-0.75*(exr[1]+exl[1])+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+(1.14564392373896*phr[17]-1.14564392373896*phl[17]-0.9682458365518543*(phr[7]+phl[7])+0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[2] = ((-0.5590169943749476*exr[11])+0.5590169943749476*exl[11]+0.4330127018922193*(exr[4]+exl[4])-0.25*exr[2]+0.25*exl[2])*c*chi+(0.5590169943749476*(phr[11]+phl[11])-0.4330127018922193*phr[4]+0.4330127018922193*phl[4]+0.25*(phr[2]+phl[2]))*c2chi; 
  incr[3] = ((-0.5590169943749476*exr[13])+0.5590169943749476*exl[13]+0.4330127018922193*(exr[5]+exl[5])-0.25*exr[3]+0.25*exl[3])*c*chi+(0.5590169943749476*(phr[13]+phl[13])-0.4330127018922193*phr[5]+0.4330127018922193*phl[5]+0.25*(phr[3]+phl[3]))*c2chi; 
  incr[4] = (0.9682458365518543*exr[11]-0.9682458365518543*exl[11]-0.75*(exr[4]+exl[4])+0.4330127018922193*exr[2]-0.4330127018922193*exl[2])*c*chi+((-0.9682458365518543*(phr[11]+phl[11]))+0.75*phr[4]-0.75*phl[4]-0.4330127018922193*(phr[2]+phl[2]))*c2chi; 
  incr[5] = (0.9682458365518543*exr[13]-0.9682458365518543*exl[13]-0.75*(exr[5]+exl[5])+0.4330127018922193*exr[3]-0.4330127018922193*exl[3])*c*chi+((-0.9682458365518543*(phr[13]+phl[13]))+0.75*phr[5]-0.75*phl[5]-0.4330127018922193*(phr[3]+phl[3]))*c2chi; 
  incr[6] = (0.4330127018922193*(exr[10]+exl[10])-0.25*exr[6]+0.25*exl[6])*c*chi+((-0.4330127018922193*phr[10])+0.4330127018922193*phl[10]+0.25*(phr[6]+phl[6]))*c2chi; 
  incr[7] = (1.479019945774904*(exr[17]+exl[17])-1.25*exr[7]+1.25*exl[7]+0.9682458365518543*(exr[1]+exl[1])-0.5590169943749475*exr[0]+0.5590169943749475*exl[0])*c*chi+((-1.479019945774904*phr[17])+1.479019945774904*phl[17]+1.25*(phr[7]+phl[7])-0.9682458365518543*phr[1]+0.9682458365518543*phl[1]+0.5590169943749475*(phr[0]+phl[0]))*c2chi; 
  incr[8] = (0.4330127018922194*(exr[12]+exl[12])-0.25*exr[8]+0.25*exl[8])*c*chi+((-0.4330127018922194*phr[12])+0.4330127018922194*phl[12]+0.25*(phr[8]+phl[8]))*c2chi; 
  incr[9] = (0.4330127018922194*(exr[15]+exl[15])-0.25*exr[9]+0.25*exl[9])*c*chi+((-0.4330127018922194*phr[15])+0.4330127018922194*phl[15]+0.25*(phr[9]+phl[9]))*c2chi; 
  incr[10] = ((-0.75*(exr[10]+exl[10]))+0.4330127018922193*exr[6]-0.4330127018922193*exl[6])*c*chi+(0.75*phr[10]-0.75*phl[10]-0.4330127018922193*(phr[6]+phl[6]))*c2chi; 
  incr[11] = ((-1.25*exr[11])+1.25*exl[11]+0.9682458365518543*(exr[4]+exl[4])-0.5590169943749476*exr[2]+0.5590169943749476*exl[2])*c*chi+(1.25*(phr[11]+phl[11])-0.9682458365518543*phr[4]+0.9682458365518543*phl[4]+0.5590169943749476*(phr[2]+phl[2]))*c2chi; 
  incr[12] = ((-0.75*(exr[12]+exl[12]))+0.4330127018922194*exr[8]-0.4330127018922194*exl[8])*c*chi+(0.75*phr[12]-0.75*phl[12]-0.4330127018922194*(phr[8]+phl[8]))*c2chi; 
  incr[13] = ((-1.25*exr[13])+1.25*exl[13]+0.9682458365518543*(exr[5]+exl[5])-0.5590169943749476*exr[3]+0.5590169943749476*exl[3])*c*chi+(1.25*(phr[13]+phl[13])-0.9682458365518543*phr[5]+0.9682458365518543*phl[5]+0.5590169943749476*(phr[3]+phl[3]))*c2chi; 
  incr[14] = (0.25*exl[14]-0.25*exr[14])*c*chi+0.25*(phr[14]+phl[14])*c2chi; 
  incr[15] = ((-0.75*(exr[15]+exl[15]))+0.4330127018922194*exr[9]-0.4330127018922194*exl[9])*c*chi+(0.75*phr[15]-0.75*phl[15]-0.4330127018922194*(phr[9]+phl[9]))*c2chi; 
  incr[16] = (0.25*exl[16]-0.25*exr[16])*c*chi+0.25*(phr[16]+phl[16])*c2chi; 
  incr[17] = ((-1.75*(exr[17]+exl[17]))+1.479019945774904*exr[7]-1.479019945774904*exl[7]-1.14564392373896*(exr[1]+exl[1])+0.6614378277661477*exr[0]-0.6614378277661477*exl[0])*c*chi+(1.75*phr[17]-1.75*phl[17]-1.479019945774904*(phr[7]+phl[7])+1.14564392373896*phr[1]-1.14564392373896*phl[1]-0.6614378277661477*(phr[0]+phl[0]))*c2chi; 
  incr[18] = (0.25*exl[18]-0.25*exr[18])*c*chi+0.25*(phr[18]+phl[18])*c2chi; 
  incr[19] = (0.25*exl[19]-0.25*exr[19])*c*chi+0.25*(phr[19]+phl[19])*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 
  outExr[4] += incr[4]*dxr1; 
  outExr[5] += incr[5]*dxr1; 
  outExr[6] += incr[6]*dxr1; 
  outExr[7] += incr[7]*dxr1; 
  outExr[8] += incr[8]*dxr1; 
  outExr[9] += incr[9]*dxr1; 
  outExr[10] += incr[10]*dxr1; 
  outExr[11] += incr[11]*dxr1; 
  outExr[12] += incr[12]*dxr1; 
  outExr[13] += incr[13]*dxr1; 
  outExr[14] += incr[14]*dxr1; 
  outExr[15] += incr[15]*dxr1; 
  outExr[16] += incr[16]*dxr1; 
  outExr[17] += incr[17]*dxr1; 
  outExr[18] += incr[18]*dxr1; 
  outExr[19] += incr[19]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 
  outExl[3] += -1.0*incr[3]*dxl1; 
  outExl[4] += incr[4]*dxl1; 
  outExl[5] += incr[5]*dxl1; 
  outExl[6] += -1.0*incr[6]*dxl1; 
  outExl[7] += -1.0*incr[7]*dxl1; 
  outExl[8] += -1.0*incr[8]*dxl1; 
  outExl[9] += -1.0*incr[9]*dxl1; 
  outExl[10] += incr[10]*dxl1; 
  outExl[11] += -1.0*incr[11]*dxl1; 
  outExl[12] += incr[12]*dxl1; 
  outExl[13] += -1.0*incr[13]*dxl1; 
  outExl[14] += -1.0*incr[14]*dxl1; 
  outExl[15] += incr[15]*dxl1; 
  outExl[16] += -1.0*incr[16]*dxl1; 
  outExl[17] += incr[17]*dxl1; 
  outExl[18] += -1.0*incr[18]*dxl1; 
  outExl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = (0.6614378277661477*(eyr[17]+eyl[17])-0.5590169943749475*eyr[7]+0.5590169943749475*eyl[7]+0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+((-0.6614378277661477*bzr[17])+0.6614378277661477*bzl[17]+0.5590169943749475*(bzr[7]+bzl[7])-0.4330127018922193*bzr[1]+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = ((-1.14564392373896*(eyr[17]+eyl[17]))+0.9682458365518543*eyr[7]-0.9682458365518543*eyl[7]-0.75*(eyr[1]+eyl[1])+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+(1.14564392373896*bzr[17]-1.14564392373896*bzl[17]-0.9682458365518543*(bzr[7]+bzl[7])+0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[2] = ((-0.5590169943749476*eyr[11])+0.5590169943749476*eyl[11]+0.4330127018922193*(eyr[4]+eyl[4])-0.25*eyr[2]+0.25*eyl[2])*tau+(0.5590169943749476*(bzr[11]+bzl[11])-0.4330127018922193*bzr[4]+0.4330127018922193*bzl[4]+0.25*(bzr[2]+bzl[2]))*c2; 
  incr[3] = ((-0.5590169943749476*eyr[13])+0.5590169943749476*eyl[13]+0.4330127018922193*(eyr[5]+eyl[5])-0.25*eyr[3]+0.25*eyl[3])*tau+(0.5590169943749476*(bzr[13]+bzl[13])-0.4330127018922193*bzr[5]+0.4330127018922193*bzl[5]+0.25*(bzr[3]+bzl[3]))*c2; 
  incr[4] = (0.9682458365518543*eyr[11]-0.9682458365518543*eyl[11]-0.75*(eyr[4]+eyl[4])+0.4330127018922193*eyr[2]-0.4330127018922193*eyl[2])*tau+((-0.9682458365518543*(bzr[11]+bzl[11]))+0.75*bzr[4]-0.75*bzl[4]-0.4330127018922193*(bzr[2]+bzl[2]))*c2; 
  incr[5] = (0.9682458365518543*eyr[13]-0.9682458365518543*eyl[13]-0.75*(eyr[5]+eyl[5])+0.4330127018922193*eyr[3]-0.4330127018922193*eyl[3])*tau+((-0.9682458365518543*(bzr[13]+bzl[13]))+0.75*bzr[5]-0.75*bzl[5]-0.4330127018922193*(bzr[3]+bzl[3]))*c2; 
  incr[6] = (0.4330127018922193*(eyr[10]+eyl[10])-0.25*eyr[6]+0.25*eyl[6])*tau+((-0.4330127018922193*bzr[10])+0.4330127018922193*bzl[10]+0.25*(bzr[6]+bzl[6]))*c2; 
  incr[7] = (1.479019945774904*(eyr[17]+eyl[17])-1.25*eyr[7]+1.25*eyl[7]+0.9682458365518543*(eyr[1]+eyl[1])-0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0])*tau+((-1.479019945774904*bzr[17])+1.479019945774904*bzl[17]+1.25*(bzr[7]+bzl[7])-0.9682458365518543*bzr[1]+0.9682458365518543*bzl[1]+0.5590169943749475*(bzr[0]+bzl[0]))*c2; 
  incr[8] = (0.4330127018922194*(eyr[12]+eyl[12])-0.25*eyr[8]+0.25*eyl[8])*tau+((-0.4330127018922194*bzr[12])+0.4330127018922194*bzl[12]+0.25*(bzr[8]+bzl[8]))*c2; 
  incr[9] = (0.4330127018922194*(eyr[15]+eyl[15])-0.25*eyr[9]+0.25*eyl[9])*tau+((-0.4330127018922194*bzr[15])+0.4330127018922194*bzl[15]+0.25*(bzr[9]+bzl[9]))*c2; 
  incr[10] = ((-0.75*(eyr[10]+eyl[10]))+0.4330127018922193*eyr[6]-0.4330127018922193*eyl[6])*tau+(0.75*bzr[10]-0.75*bzl[10]-0.4330127018922193*(bzr[6]+bzl[6]))*c2; 
  incr[11] = ((-1.25*eyr[11])+1.25*eyl[11]+0.9682458365518543*(eyr[4]+eyl[4])-0.5590169943749476*eyr[2]+0.5590169943749476*eyl[2])*tau+(1.25*(bzr[11]+bzl[11])-0.9682458365518543*bzr[4]+0.9682458365518543*bzl[4]+0.5590169943749476*(bzr[2]+bzl[2]))*c2; 
  incr[12] = ((-0.75*(eyr[12]+eyl[12]))+0.4330127018922194*eyr[8]-0.4330127018922194*eyl[8])*tau+(0.75*bzr[12]-0.75*bzl[12]-0.4330127018922194*(bzr[8]+bzl[8]))*c2; 
  incr[13] = ((-1.25*eyr[13])+1.25*eyl[13]+0.9682458365518543*(eyr[5]+eyl[5])-0.5590169943749476*eyr[3]+0.5590169943749476*eyl[3])*tau+(1.25*(bzr[13]+bzl[13])-0.9682458365518543*bzr[5]+0.9682458365518543*bzl[5]+0.5590169943749476*(bzr[3]+bzl[3]))*c2; 
  incr[14] = 0.25*(eyl[14]*tau+(bzr[14]+bzl[14])*c2)-0.25*eyr[14]*tau; 
  incr[15] = ((-0.75*(eyr[15]+eyl[15]))+0.4330127018922194*eyr[9]-0.4330127018922194*eyl[9])*tau+(0.75*bzr[15]-0.75*bzl[15]-0.4330127018922194*(bzr[9]+bzl[9]))*c2; 
  incr[16] = 0.25*(eyl[16]*tau+(bzr[16]+bzl[16])*c2)-0.25*eyr[16]*tau; 
  incr[17] = ((-1.75*(eyr[17]+eyl[17]))+1.479019945774904*eyr[7]-1.479019945774904*eyl[7]-1.14564392373896*(eyr[1]+eyl[1])+0.6614378277661477*eyr[0]-0.6614378277661477*eyl[0])*tau+(1.75*bzr[17]-1.75*bzl[17]-1.479019945774904*(bzr[7]+bzl[7])+1.14564392373896*bzr[1]-1.14564392373896*bzl[1]-0.6614378277661477*(bzr[0]+bzl[0]))*c2; 
  incr[18] = 0.25*(eyl[18]*tau+(bzr[18]+bzl[18])*c2)-0.25*eyr[18]*tau; 
  incr[19] = 0.25*(eyl[19]*tau+(bzr[19]+bzl[19])*c2)-0.25*eyr[19]*tau; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 
  outEyr[4] += incr[4]*dxr1; 
  outEyr[5] += incr[5]*dxr1; 
  outEyr[6] += incr[6]*dxr1; 
  outEyr[7] += incr[7]*dxr1; 
  outEyr[8] += incr[8]*dxr1; 
  outEyr[9] += incr[9]*dxr1; 
  outEyr[10] += incr[10]*dxr1; 
  outEyr[11] += incr[11]*dxr1; 
  outEyr[12] += incr[12]*dxr1; 
  outEyr[13] += incr[13]*dxr1; 
  outEyr[14] += incr[14]*dxr1; 
  outEyr[15] += incr[15]*dxr1; 
  outEyr[16] += incr[16]*dxr1; 
  outEyr[17] += incr[17]*dxr1; 
  outEyr[18] += incr[18]*dxr1; 
  outEyr[19] += incr[19]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 
  outEyl[3] += -1.0*incr[3]*dxl1; 
  outEyl[4] += incr[4]*dxl1; 
  outEyl[5] += incr[5]*dxl1; 
  outEyl[6] += -1.0*incr[6]*dxl1; 
  outEyl[7] += -1.0*incr[7]*dxl1; 
  outEyl[8] += -1.0*incr[8]*dxl1; 
  outEyl[9] += -1.0*incr[9]*dxl1; 
  outEyl[10] += incr[10]*dxl1; 
  outEyl[11] += -1.0*incr[11]*dxl1; 
  outEyl[12] += incr[12]*dxl1; 
  outEyl[13] += -1.0*incr[13]*dxl1; 
  outEyl[14] += -1.0*incr[14]*dxl1; 
  outEyl[15] += incr[15]*dxl1; 
  outEyl[16] += -1.0*incr[16]*dxl1; 
  outEyl[17] += incr[17]*dxl1; 
  outEyl[18] += -1.0*incr[18]*dxl1; 
  outEyl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = (0.6614378277661477*(ezr[17]+ezl[17])-0.5590169943749475*ezr[7]+0.5590169943749475*ezl[7]+0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+(0.6614378277661477*byr[17]-0.6614378277661477*byl[17]-0.5590169943749475*(byr[7]+byl[7])+0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = ((-1.14564392373896*(ezr[17]+ezl[17]))+0.9682458365518543*ezr[7]-0.9682458365518543*ezl[7]-0.75*(ezr[1]+ezl[1])+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+((-1.14564392373896*byr[17])+1.14564392373896*byl[17]+0.9682458365518543*(byr[7]+byl[7])-0.75*byr[1]+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[2] = ((-0.5590169943749476*ezr[11])+0.5590169943749476*ezl[11]+0.4330127018922193*(ezr[4]+ezl[4])-0.25*ezr[2]+0.25*ezl[2])*tau+((-0.5590169943749476*(byr[11]+byl[11]))+0.4330127018922193*byr[4]-0.4330127018922193*byl[4]-0.25*(byr[2]+byl[2]))*c2; 
  incr[3] = ((-0.5590169943749476*ezr[13])+0.5590169943749476*ezl[13]+0.4330127018922193*(ezr[5]+ezl[5])-0.25*ezr[3]+0.25*ezl[3])*tau+((-0.5590169943749476*(byr[13]+byl[13]))+0.4330127018922193*byr[5]-0.4330127018922193*byl[5]-0.25*(byr[3]+byl[3]))*c2; 
  incr[4] = (0.9682458365518543*ezr[11]-0.9682458365518543*ezl[11]-0.75*(ezr[4]+ezl[4])+0.4330127018922193*ezr[2]-0.4330127018922193*ezl[2])*tau+(0.9682458365518543*(byr[11]+byl[11])-0.75*byr[4]+0.75*byl[4]+0.4330127018922193*(byr[2]+byl[2]))*c2; 
  incr[5] = (0.9682458365518543*ezr[13]-0.9682458365518543*ezl[13]-0.75*(ezr[5]+ezl[5])+0.4330127018922193*ezr[3]-0.4330127018922193*ezl[3])*tau+(0.9682458365518543*(byr[13]+byl[13])-0.75*byr[5]+0.75*byl[5]+0.4330127018922193*(byr[3]+byl[3]))*c2; 
  incr[6] = (0.4330127018922193*(ezr[10]+ezl[10])-0.25*ezr[6]+0.25*ezl[6])*tau+(0.4330127018922193*byr[10]-0.4330127018922193*byl[10]-0.25*(byr[6]+byl[6]))*c2; 
  incr[7] = (1.479019945774904*(ezr[17]+ezl[17])-1.25*ezr[7]+1.25*ezl[7]+0.9682458365518543*(ezr[1]+ezl[1])-0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0])*tau+(1.479019945774904*byr[17]-1.479019945774904*byl[17]-1.25*(byr[7]+byl[7])+0.9682458365518543*byr[1]-0.9682458365518543*byl[1]-0.5590169943749475*(byr[0]+byl[0]))*c2; 
  incr[8] = (0.4330127018922194*(ezr[12]+ezl[12])-0.25*ezr[8]+0.25*ezl[8])*tau+(0.4330127018922194*byr[12]-0.4330127018922194*byl[12]-0.25*(byr[8]+byl[8]))*c2; 
  incr[9] = (0.4330127018922194*(ezr[15]+ezl[15])-0.25*ezr[9]+0.25*ezl[9])*tau+(0.4330127018922194*byr[15]-0.4330127018922194*byl[15]-0.25*(byr[9]+byl[9]))*c2; 
  incr[10] = ((-0.75*(ezr[10]+ezl[10]))+0.4330127018922193*ezr[6]-0.4330127018922193*ezl[6])*tau+((-0.75*byr[10])+0.75*byl[10]+0.4330127018922193*(byr[6]+byl[6]))*c2; 
  incr[11] = ((-1.25*ezr[11])+1.25*ezl[11]+0.9682458365518543*(ezr[4]+ezl[4])-0.5590169943749476*ezr[2]+0.5590169943749476*ezl[2])*tau+((-1.25*(byr[11]+byl[11]))+0.9682458365518543*byr[4]-0.9682458365518543*byl[4]-0.5590169943749476*(byr[2]+byl[2]))*c2; 
  incr[12] = ((-0.75*(ezr[12]+ezl[12]))+0.4330127018922194*ezr[8]-0.4330127018922194*ezl[8])*tau+((-0.75*byr[12])+0.75*byl[12]+0.4330127018922194*(byr[8]+byl[8]))*c2; 
  incr[13] = ((-1.25*ezr[13])+1.25*ezl[13]+0.9682458365518543*(ezr[5]+ezl[5])-0.5590169943749476*ezr[3]+0.5590169943749476*ezl[3])*tau+((-1.25*(byr[13]+byl[13]))+0.9682458365518543*byr[5]-0.9682458365518543*byl[5]-0.5590169943749476*(byr[3]+byl[3]))*c2; 
  incr[14] = (0.25*ezl[14]-0.25*ezr[14])*tau-0.25*(byr[14]+byl[14])*c2; 
  incr[15] = ((-0.75*(ezr[15]+ezl[15]))+0.4330127018922194*ezr[9]-0.4330127018922194*ezl[9])*tau+((-0.75*byr[15])+0.75*byl[15]+0.4330127018922194*(byr[9]+byl[9]))*c2; 
  incr[16] = (0.25*ezl[16]-0.25*ezr[16])*tau-0.25*(byr[16]+byl[16])*c2; 
  incr[17] = ((-1.75*(ezr[17]+ezl[17]))+1.479019945774904*ezr[7]-1.479019945774904*ezl[7]-1.14564392373896*(ezr[1]+ezl[1])+0.6614378277661477*ezr[0]-0.6614378277661477*ezl[0])*tau+((-1.75*byr[17])+1.75*byl[17]+1.479019945774904*(byr[7]+byl[7])-1.14564392373896*byr[1]+1.14564392373896*byl[1]+0.6614378277661477*(byr[0]+byl[0]))*c2; 
  incr[18] = (0.25*ezl[18]-0.25*ezr[18])*tau-0.25*(byr[18]+byl[18])*c2; 
  incr[19] = (0.25*ezl[19]-0.25*ezr[19])*tau-0.25*(byr[19]+byl[19])*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 
  outEzr[4] += incr[4]*dxr1; 
  outEzr[5] += incr[5]*dxr1; 
  outEzr[6] += incr[6]*dxr1; 
  outEzr[7] += incr[7]*dxr1; 
  outEzr[8] += incr[8]*dxr1; 
  outEzr[9] += incr[9]*dxr1; 
  outEzr[10] += incr[10]*dxr1; 
  outEzr[11] += incr[11]*dxr1; 
  outEzr[12] += incr[12]*dxr1; 
  outEzr[13] += incr[13]*dxr1; 
  outEzr[14] += incr[14]*dxr1; 
  outEzr[15] += incr[15]*dxr1; 
  outEzr[16] += incr[16]*dxr1; 
  outEzr[17] += incr[17]*dxr1; 
  outEzr[18] += incr[18]*dxr1; 
  outEzr[19] += incr[19]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 
  outEzl[3] += -1.0*incr[3]*dxl1; 
  outEzl[4] += incr[4]*dxl1; 
  outEzl[5] += incr[5]*dxl1; 
  outEzl[6] += -1.0*incr[6]*dxl1; 
  outEzl[7] += -1.0*incr[7]*dxl1; 
  outEzl[8] += -1.0*incr[8]*dxl1; 
  outEzl[9] += -1.0*incr[9]*dxl1; 
  outEzl[10] += incr[10]*dxl1; 
  outEzl[11] += -1.0*incr[11]*dxl1; 
  outEzl[12] += incr[12]*dxl1; 
  outEzl[13] += -1.0*incr[13]*dxl1; 
  outEzl[14] += -1.0*incr[14]*dxl1; 
  outEzl[15] += incr[15]*dxl1; 
  outEzl[16] += -1.0*incr[16]*dxl1; 
  outEzl[17] += incr[17]*dxl1; 
  outEzl[18] += -1.0*incr[18]*dxl1; 
  outEzl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(bxr[17]+bxl[17])-0.5590169943749475*bxr[7]+0.5590169943749475*bxl[7]+0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c-0.6614378277661477*psr[17]+0.6614378277661477*psl[17]+0.5590169943749475*(psr[7]+psl[7])-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (((-1.14564392373896*(bxr[17]+bxl[17]))+0.9682458365518543*bxr[7]-0.9682458365518543*bxl[7]-0.75*(bxr[1]+bxl[1])+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c+1.14564392373896*psr[17]-1.14564392373896*psl[17]-0.9682458365518543*(psr[7]+psl[7])+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[2] = (((-0.5590169943749476*bxr[11])+0.5590169943749476*bxl[11]+0.4330127018922193*(bxr[4]+bxl[4])-0.25*bxr[2]+0.25*bxl[2])*c+0.5590169943749476*(psr[11]+psl[11])-0.4330127018922193*psr[4]+0.4330127018922193*psl[4]+0.25*(psr[2]+psl[2]))*gamma; 
  incr[3] = (((-0.5590169943749476*bxr[13])+0.5590169943749476*bxl[13]+0.4330127018922193*(bxr[5]+bxl[5])-0.25*bxr[3]+0.25*bxl[3])*c+0.5590169943749476*(psr[13]+psl[13])-0.4330127018922193*psr[5]+0.4330127018922193*psl[5]+0.25*(psr[3]+psl[3]))*gamma; 
  incr[4] = ((0.9682458365518543*bxr[11]-0.9682458365518543*bxl[11]-0.75*(bxr[4]+bxl[4])+0.4330127018922193*bxr[2]-0.4330127018922193*bxl[2])*c-0.9682458365518543*(psr[11]+psl[11])+0.75*psr[4]-0.75*psl[4]-0.4330127018922193*(psr[2]+psl[2]))*gamma; 
  incr[5] = ((0.9682458365518543*bxr[13]-0.9682458365518543*bxl[13]-0.75*(bxr[5]+bxl[5])+0.4330127018922193*bxr[3]-0.4330127018922193*bxl[3])*c-0.9682458365518543*(psr[13]+psl[13])+0.75*psr[5]-0.75*psl[5]-0.4330127018922193*(psr[3]+psl[3]))*gamma; 
  incr[6] = ((0.4330127018922193*(bxr[10]+bxl[10])-0.25*bxr[6]+0.25*bxl[6])*c-0.4330127018922193*psr[10]+0.4330127018922193*psl[10]+0.25*(psr[6]+psl[6]))*gamma; 
  incr[7] = ((1.479019945774904*(bxr[17]+bxl[17])-1.25*bxr[7]+1.25*bxl[7]+0.9682458365518543*(bxr[1]+bxl[1])-0.5590169943749475*bxr[0]+0.5590169943749475*bxl[0])*c-1.479019945774904*psr[17]+1.479019945774904*psl[17]+1.25*(psr[7]+psl[7])-0.9682458365518543*psr[1]+0.9682458365518543*psl[1]+0.5590169943749475*(psr[0]+psl[0]))*gamma; 
  incr[8] = ((0.4330127018922194*(bxr[12]+bxl[12])-0.25*bxr[8]+0.25*bxl[8])*c-0.4330127018922194*psr[12]+0.4330127018922194*psl[12]+0.25*(psr[8]+psl[8]))*gamma; 
  incr[9] = ((0.4330127018922194*(bxr[15]+bxl[15])-0.25*bxr[9]+0.25*bxl[9])*c-0.4330127018922194*psr[15]+0.4330127018922194*psl[15]+0.25*(psr[9]+psl[9]))*gamma; 
  incr[10] = (((-0.75*(bxr[10]+bxl[10]))+0.4330127018922193*bxr[6]-0.4330127018922193*bxl[6])*c+0.75*psr[10]-0.75*psl[10]-0.4330127018922193*(psr[6]+psl[6]))*gamma; 
  incr[11] = (((-1.25*bxr[11])+1.25*bxl[11]+0.9682458365518543*(bxr[4]+bxl[4])-0.5590169943749476*bxr[2]+0.5590169943749476*bxl[2])*c+1.25*(psr[11]+psl[11])-0.9682458365518543*psr[4]+0.9682458365518543*psl[4]+0.5590169943749476*(psr[2]+psl[2]))*gamma; 
  incr[12] = (((-0.75*(bxr[12]+bxl[12]))+0.4330127018922194*bxr[8]-0.4330127018922194*bxl[8])*c+0.75*psr[12]-0.75*psl[12]-0.4330127018922194*(psr[8]+psl[8]))*gamma; 
  incr[13] = (((-1.25*bxr[13])+1.25*bxl[13]+0.9682458365518543*(bxr[5]+bxl[5])-0.5590169943749476*bxr[3]+0.5590169943749476*bxl[3])*c+1.25*(psr[13]+psl[13])-0.9682458365518543*psr[5]+0.9682458365518543*psl[5]+0.5590169943749476*(psr[3]+psl[3]))*gamma; 
  incr[14] = (0.25*(bxl[14]*c+psr[14]+psl[14])-0.25*bxr[14]*c)*gamma; 
  incr[15] = (((-0.75*(bxr[15]+bxl[15]))+0.4330127018922194*bxr[9]-0.4330127018922194*bxl[9])*c+0.75*psr[15]-0.75*psl[15]-0.4330127018922194*(psr[9]+psl[9]))*gamma; 
  incr[16] = (0.25*(bxl[16]*c+psr[16]+psl[16])-0.25*bxr[16]*c)*gamma; 
  incr[17] = (((-1.75*(bxr[17]+bxl[17]))+1.479019945774904*bxr[7]-1.479019945774904*bxl[7]-1.14564392373896*(bxr[1]+bxl[1])+0.6614378277661477*bxr[0]-0.6614378277661477*bxl[0])*c+1.75*psr[17]-1.75*psl[17]-1.479019945774904*(psr[7]+psl[7])+1.14564392373896*psr[1]-1.14564392373896*psl[1]-0.6614378277661477*(psr[0]+psl[0]))*gamma; 
  incr[18] = (0.25*(bxl[18]*c+psr[18]+psl[18])-0.25*bxr[18]*c)*gamma; 
  incr[19] = (0.25*(bxl[19]*c+psr[19]+psl[19])-0.25*bxr[19]*c)*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 
  outBxr[4] += incr[4]*dxr1; 
  outBxr[5] += incr[5]*dxr1; 
  outBxr[6] += incr[6]*dxr1; 
  outBxr[7] += incr[7]*dxr1; 
  outBxr[8] += incr[8]*dxr1; 
  outBxr[9] += incr[9]*dxr1; 
  outBxr[10] += incr[10]*dxr1; 
  outBxr[11] += incr[11]*dxr1; 
  outBxr[12] += incr[12]*dxr1; 
  outBxr[13] += incr[13]*dxr1; 
  outBxr[14] += incr[14]*dxr1; 
  outBxr[15] += incr[15]*dxr1; 
  outBxr[16] += incr[16]*dxr1; 
  outBxr[17] += incr[17]*dxr1; 
  outBxr[18] += incr[18]*dxr1; 
  outBxr[19] += incr[19]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 
  outBxl[3] += -1.0*incr[3]*dxl1; 
  outBxl[4] += incr[4]*dxl1; 
  outBxl[5] += incr[5]*dxl1; 
  outBxl[6] += -1.0*incr[6]*dxl1; 
  outBxl[7] += -1.0*incr[7]*dxl1; 
  outBxl[8] += -1.0*incr[8]*dxl1; 
  outBxl[9] += -1.0*incr[9]*dxl1; 
  outBxl[10] += incr[10]*dxl1; 
  outBxl[11] += -1.0*incr[11]*dxl1; 
  outBxl[12] += incr[12]*dxl1; 
  outBxl[13] += -1.0*incr[13]*dxl1; 
  outBxl[14] += -1.0*incr[14]*dxl1; 
  outBxl[15] += incr[15]*dxl1; 
  outBxl[16] += -1.0*incr[16]*dxl1; 
  outBxl[17] += incr[17]*dxl1; 
  outBxl[18] += -1.0*incr[18]*dxl1; 
  outBxl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(byr[17]+byl[17])-0.5590169943749475*byr[7]+0.5590169943749475*byl[7]+0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau+0.6614378277661477*ezr[17]-0.6614378277661477*ezl[17]-0.5590169943749475*(ezr[7]+ezl[7])+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = (((-1.14564392373896*(byr[17]+byl[17]))+0.9682458365518543*byr[7]-0.9682458365518543*byl[7]-0.75*(byr[1]+byl[1])+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau-1.14564392373896*ezr[17]+1.14564392373896*ezl[17]+0.9682458365518543*(ezr[7]+ezl[7])-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[2] = (((-0.5590169943749476*byr[11])+0.5590169943749476*byl[11]+0.4330127018922193*(byr[4]+byl[4])-0.25*byr[2]+0.25*byl[2])*c2)/tau-0.5590169943749476*(ezr[11]+ezl[11])+0.4330127018922193*ezr[4]-0.4330127018922193*ezl[4]-0.25*(ezr[2]+ezl[2]); 
  incr[3] = (((-0.5590169943749476*byr[13])+0.5590169943749476*byl[13]+0.4330127018922193*(byr[5]+byl[5])-0.25*byr[3]+0.25*byl[3])*c2)/tau-0.5590169943749476*(ezr[13]+ezl[13])+0.4330127018922193*ezr[5]-0.4330127018922193*ezl[5]-0.25*(ezr[3]+ezl[3]); 
  incr[4] = ((0.9682458365518543*byr[11]-0.9682458365518543*byl[11]-0.75*(byr[4]+byl[4])+0.4330127018922193*byr[2]-0.4330127018922193*byl[2])*c2)/tau+0.9682458365518543*(ezr[11]+ezl[11])-0.75*ezr[4]+0.75*ezl[4]+0.4330127018922193*(ezr[2]+ezl[2]); 
  incr[5] = ((0.9682458365518543*byr[13]-0.9682458365518543*byl[13]-0.75*(byr[5]+byl[5])+0.4330127018922193*byr[3]-0.4330127018922193*byl[3])*c2)/tau+0.9682458365518543*(ezr[13]+ezl[13])-0.75*ezr[5]+0.75*ezl[5]+0.4330127018922193*(ezr[3]+ezl[3]); 
  incr[6] = ((0.4330127018922193*(byr[10]+byl[10])-0.25*byr[6]+0.25*byl[6])*c2)/tau+0.4330127018922193*ezr[10]-0.4330127018922193*ezl[10]-0.25*(ezr[6]+ezl[6]); 
  incr[7] = ((1.479019945774904*(byr[17]+byl[17])-1.25*byr[7]+1.25*byl[7]+0.9682458365518543*(byr[1]+byl[1])-0.5590169943749475*byr[0]+0.5590169943749475*byl[0])*c2)/tau+1.479019945774904*ezr[17]-1.479019945774904*ezl[17]-1.25*(ezr[7]+ezl[7])+0.9682458365518543*ezr[1]-0.9682458365518543*ezl[1]-0.5590169943749475*(ezr[0]+ezl[0]); 
  incr[8] = ((0.4330127018922194*(byr[12]+byl[12])-0.25*byr[8]+0.25*byl[8])*c2)/tau+0.4330127018922194*ezr[12]-0.4330127018922194*ezl[12]-0.25*(ezr[8]+ezl[8]); 
  incr[9] = ((0.4330127018922194*(byr[15]+byl[15])-0.25*byr[9]+0.25*byl[9])*c2)/tau+0.4330127018922194*ezr[15]-0.4330127018922194*ezl[15]-0.25*(ezr[9]+ezl[9]); 
  incr[10] = (((-0.75*(byr[10]+byl[10]))+0.4330127018922193*byr[6]-0.4330127018922193*byl[6])*c2)/tau-0.75*ezr[10]+0.75*ezl[10]+0.4330127018922193*(ezr[6]+ezl[6]); 
  incr[11] = (((-1.25*byr[11])+1.25*byl[11]+0.9682458365518543*(byr[4]+byl[4])-0.5590169943749476*byr[2]+0.5590169943749476*byl[2])*c2)/tau-1.25*(ezr[11]+ezl[11])+0.9682458365518543*ezr[4]-0.9682458365518543*ezl[4]-0.5590169943749476*(ezr[2]+ezl[2]); 
  incr[12] = (((-0.75*(byr[12]+byl[12]))+0.4330127018922194*byr[8]-0.4330127018922194*byl[8])*c2)/tau-0.75*ezr[12]+0.75*ezl[12]+0.4330127018922194*(ezr[8]+ezl[8]); 
  incr[13] = (((-1.25*byr[13])+1.25*byl[13]+0.9682458365518543*(byr[5]+byl[5])-0.5590169943749476*byr[3]+0.5590169943749476*byl[3])*c2)/tau-1.25*(ezr[13]+ezl[13])+0.9682458365518543*ezr[5]-0.9682458365518543*ezl[5]-0.5590169943749476*(ezr[3]+ezl[3]); 
  incr[14] = ((0.25*byl[14]-0.25*byr[14])*c2)/tau-0.25*(ezr[14]+ezl[14]); 
  incr[15] = (((-0.75*(byr[15]+byl[15]))+0.4330127018922194*byr[9]-0.4330127018922194*byl[9])*c2)/tau-0.75*ezr[15]+0.75*ezl[15]+0.4330127018922194*(ezr[9]+ezl[9]); 
  incr[16] = ((0.25*byl[16]-0.25*byr[16])*c2)/tau-0.25*(ezr[16]+ezl[16]); 
  incr[17] = (((-1.75*(byr[17]+byl[17]))+1.479019945774904*byr[7]-1.479019945774904*byl[7]-1.14564392373896*(byr[1]+byl[1])+0.6614378277661477*byr[0]-0.6614378277661477*byl[0])*c2)/tau-1.75*ezr[17]+1.75*ezl[17]+1.479019945774904*(ezr[7]+ezl[7])-1.14564392373896*ezr[1]+1.14564392373896*ezl[1]+0.6614378277661477*(ezr[0]+ezl[0]); 
  incr[18] = ((0.25*byl[18]-0.25*byr[18])*c2)/tau-0.25*(ezr[18]+ezl[18]); 
  incr[19] = ((0.25*byl[19]-0.25*byr[19])*c2)/tau-0.25*(ezr[19]+ezl[19]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 
  outByr[4] += incr[4]*dxr1; 
  outByr[5] += incr[5]*dxr1; 
  outByr[6] += incr[6]*dxr1; 
  outByr[7] += incr[7]*dxr1; 
  outByr[8] += incr[8]*dxr1; 
  outByr[9] += incr[9]*dxr1; 
  outByr[10] += incr[10]*dxr1; 
  outByr[11] += incr[11]*dxr1; 
  outByr[12] += incr[12]*dxr1; 
  outByr[13] += incr[13]*dxr1; 
  outByr[14] += incr[14]*dxr1; 
  outByr[15] += incr[15]*dxr1; 
  outByr[16] += incr[16]*dxr1; 
  outByr[17] += incr[17]*dxr1; 
  outByr[18] += incr[18]*dxr1; 
  outByr[19] += incr[19]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 
  outByl[3] += -1.0*incr[3]*dxl1; 
  outByl[4] += incr[4]*dxl1; 
  outByl[5] += incr[5]*dxl1; 
  outByl[6] += -1.0*incr[6]*dxl1; 
  outByl[7] += -1.0*incr[7]*dxl1; 
  outByl[8] += -1.0*incr[8]*dxl1; 
  outByl[9] += -1.0*incr[9]*dxl1; 
  outByl[10] += incr[10]*dxl1; 
  outByl[11] += -1.0*incr[11]*dxl1; 
  outByl[12] += incr[12]*dxl1; 
  outByl[13] += -1.0*incr[13]*dxl1; 
  outByl[14] += -1.0*incr[14]*dxl1; 
  outByl[15] += incr[15]*dxl1; 
  outByl[16] += -1.0*incr[16]*dxl1; 
  outByl[17] += incr[17]*dxl1; 
  outByl[18] += -1.0*incr[18]*dxl1; 
  outByl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(bzr[17]+bzl[17])-0.5590169943749475*bzr[7]+0.5590169943749475*bzl[7]+0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau-0.6614378277661477*eyr[17]+0.6614378277661477*eyl[17]+0.5590169943749475*(eyr[7]+eyl[7])-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = (((-1.14564392373896*(bzr[17]+bzl[17]))+0.9682458365518543*bzr[7]-0.9682458365518543*bzl[7]-0.75*(bzr[1]+bzl[1])+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau+1.14564392373896*eyr[17]-1.14564392373896*eyl[17]-0.9682458365518543*(eyr[7]+eyl[7])+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[2] = (((-0.5590169943749476*bzr[11])+0.5590169943749476*bzl[11]+0.4330127018922193*(bzr[4]+bzl[4])-0.25*bzr[2]+0.25*bzl[2])*c2)/tau+0.5590169943749476*(eyr[11]+eyl[11])-0.4330127018922193*eyr[4]+0.4330127018922193*eyl[4]+0.25*(eyr[2]+eyl[2]); 
  incr[3] = (((-0.5590169943749476*bzr[13])+0.5590169943749476*bzl[13]+0.4330127018922193*(bzr[5]+bzl[5])-0.25*bzr[3]+0.25*bzl[3])*c2)/tau+0.5590169943749476*(eyr[13]+eyl[13])-0.4330127018922193*eyr[5]+0.4330127018922193*eyl[5]+0.25*(eyr[3]+eyl[3]); 
  incr[4] = ((0.9682458365518543*bzr[11]-0.9682458365518543*bzl[11]-0.75*(bzr[4]+bzl[4])+0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2])*c2)/tau-0.9682458365518543*(eyr[11]+eyl[11])+0.75*eyr[4]-0.75*eyl[4]-0.4330127018922193*(eyr[2]+eyl[2]); 
  incr[5] = ((0.9682458365518543*bzr[13]-0.9682458365518543*bzl[13]-0.75*(bzr[5]+bzl[5])+0.4330127018922193*bzr[3]-0.4330127018922193*bzl[3])*c2)/tau-0.9682458365518543*(eyr[13]+eyl[13])+0.75*eyr[5]-0.75*eyl[5]-0.4330127018922193*(eyr[3]+eyl[3]); 
  incr[6] = ((0.4330127018922193*(bzr[10]+bzl[10])-0.25*bzr[6]+0.25*bzl[6])*c2)/tau-0.4330127018922193*eyr[10]+0.4330127018922193*eyl[10]+0.25*(eyr[6]+eyl[6]); 
  incr[7] = ((1.479019945774904*(bzr[17]+bzl[17])-1.25*bzr[7]+1.25*bzl[7]+0.9682458365518543*(bzr[1]+bzl[1])-0.5590169943749475*bzr[0]+0.5590169943749475*bzl[0])*c2)/tau-1.479019945774904*eyr[17]+1.479019945774904*eyl[17]+1.25*(eyr[7]+eyl[7])-0.9682458365518543*eyr[1]+0.9682458365518543*eyl[1]+0.5590169943749475*(eyr[0]+eyl[0]); 
  incr[8] = ((0.4330127018922194*(bzr[12]+bzl[12])-0.25*bzr[8]+0.25*bzl[8])*c2)/tau-0.4330127018922194*eyr[12]+0.4330127018922194*eyl[12]+0.25*(eyr[8]+eyl[8]); 
  incr[9] = ((0.4330127018922194*(bzr[15]+bzl[15])-0.25*bzr[9]+0.25*bzl[9])*c2)/tau-0.4330127018922194*eyr[15]+0.4330127018922194*eyl[15]+0.25*(eyr[9]+eyl[9]); 
  incr[10] = (((-0.75*(bzr[10]+bzl[10]))+0.4330127018922193*bzr[6]-0.4330127018922193*bzl[6])*c2)/tau+0.75*eyr[10]-0.75*eyl[10]-0.4330127018922193*(eyr[6]+eyl[6]); 
  incr[11] = (((-1.25*bzr[11])+1.25*bzl[11]+0.9682458365518543*(bzr[4]+bzl[4])-0.5590169943749476*bzr[2]+0.5590169943749476*bzl[2])*c2)/tau+1.25*(eyr[11]+eyl[11])-0.9682458365518543*eyr[4]+0.9682458365518543*eyl[4]+0.5590169943749476*(eyr[2]+eyl[2]); 
  incr[12] = (((-0.75*(bzr[12]+bzl[12]))+0.4330127018922194*bzr[8]-0.4330127018922194*bzl[8])*c2)/tau+0.75*eyr[12]-0.75*eyl[12]-0.4330127018922194*(eyr[8]+eyl[8]); 
  incr[13] = (((-1.25*bzr[13])+1.25*bzl[13]+0.9682458365518543*(bzr[5]+bzl[5])-0.5590169943749476*bzr[3]+0.5590169943749476*bzl[3])*c2)/tau+1.25*(eyr[13]+eyl[13])-0.9682458365518543*eyr[5]+0.9682458365518543*eyl[5]+0.5590169943749476*(eyr[3]+eyl[3]); 
  incr[14] = ((0.25*bzl[14]-0.25*bzr[14])*c2)/tau+0.25*(eyr[14]+eyl[14]); 
  incr[15] = (((-0.75*(bzr[15]+bzl[15]))+0.4330127018922194*bzr[9]-0.4330127018922194*bzl[9])*c2)/tau+0.75*eyr[15]-0.75*eyl[15]-0.4330127018922194*(eyr[9]+eyl[9]); 
  incr[16] = ((0.25*bzl[16]-0.25*bzr[16])*c2)/tau+0.25*(eyr[16]+eyl[16]); 
  incr[17] = (((-1.75*(bzr[17]+bzl[17]))+1.479019945774904*bzr[7]-1.479019945774904*bzl[7]-1.14564392373896*(bzr[1]+bzl[1])+0.6614378277661477*bzr[0]-0.6614378277661477*bzl[0])*c2)/tau+1.75*eyr[17]-1.75*eyl[17]-1.479019945774904*(eyr[7]+eyl[7])+1.14564392373896*eyr[1]-1.14564392373896*eyl[1]-0.6614378277661477*(eyr[0]+eyl[0]); 
  incr[18] = ((0.25*bzl[18]-0.25*bzr[18])*c2)/tau+0.25*(eyr[18]+eyl[18]); 
  incr[19] = ((0.25*bzl[19]-0.25*bzr[19])*c2)/tau+0.25*(eyr[19]+eyl[19]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 
  outBzr[4] += incr[4]*dxr1; 
  outBzr[5] += incr[5]*dxr1; 
  outBzr[6] += incr[6]*dxr1; 
  outBzr[7] += incr[7]*dxr1; 
  outBzr[8] += incr[8]*dxr1; 
  outBzr[9] += incr[9]*dxr1; 
  outBzr[10] += incr[10]*dxr1; 
  outBzr[11] += incr[11]*dxr1; 
  outBzr[12] += incr[12]*dxr1; 
  outBzr[13] += incr[13]*dxr1; 
  outBzr[14] += incr[14]*dxr1; 
  outBzr[15] += incr[15]*dxr1; 
  outBzr[16] += incr[16]*dxr1; 
  outBzr[17] += incr[17]*dxr1; 
  outBzr[18] += incr[18]*dxr1; 
  outBzr[19] += incr[19]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 
  outBzl[3] += -1.0*incr[3]*dxl1; 
  outBzl[4] += incr[4]*dxl1; 
  outBzl[5] += incr[5]*dxl1; 
  outBzl[6] += -1.0*incr[6]*dxl1; 
  outBzl[7] += -1.0*incr[7]*dxl1; 
  outBzl[8] += -1.0*incr[8]*dxl1; 
  outBzl[9] += -1.0*incr[9]*dxl1; 
  outBzl[10] += incr[10]*dxl1; 
  outBzl[11] += -1.0*incr[11]*dxl1; 
  outBzl[12] += incr[12]*dxl1; 
  outBzl[13] += -1.0*incr[13]*dxl1; 
  outBzl[14] += -1.0*incr[14]*dxl1; 
  outBzl[15] += incr[15]*dxl1; 
  outBzl[16] += -1.0*incr[16]*dxl1; 
  outBzl[17] += incr[17]*dxl1; 
  outBzl[18] += -1.0*incr[18]*dxl1; 
  outBzl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(phr[17]+phl[17])-0.5590169943749475*phr[7]+0.5590169943749475*phl[7]+0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c-0.6614378277661477*exr[17]+0.6614378277661477*exl[17]+0.5590169943749475*(exr[7]+exl[7])-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = (((-1.14564392373896*(phr[17]+phl[17]))+0.9682458365518543*phr[7]-0.9682458365518543*phl[7]-0.75*(phr[1]+phl[1])+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+1.14564392373896*exr[17]-1.14564392373896*exl[17]-0.9682458365518543*(exr[7]+exl[7])+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 
  incr[2] = (((-0.5590169943749476*phr[11])+0.5590169943749476*phl[11]+0.4330127018922193*(phr[4]+phl[4])-0.25*phr[2]+0.25*phl[2])*c+0.5590169943749476*(exr[11]+exl[11])-0.4330127018922193*exr[4]+0.4330127018922193*exl[4]+0.25*(exr[2]+exl[2]))*chi; 
  incr[3] = (((-0.5590169943749476*phr[13])+0.5590169943749476*phl[13]+0.4330127018922193*(phr[5]+phl[5])-0.25*phr[3]+0.25*phl[3])*c+0.5590169943749476*(exr[13]+exl[13])-0.4330127018922193*exr[5]+0.4330127018922193*exl[5]+0.25*(exr[3]+exl[3]))*chi; 
  incr[4] = ((0.9682458365518543*phr[11]-0.9682458365518543*phl[11]-0.75*(phr[4]+phl[4])+0.4330127018922193*phr[2]-0.4330127018922193*phl[2])*c-0.9682458365518543*(exr[11]+exl[11])+0.75*exr[4]-0.75*exl[4]-0.4330127018922193*(exr[2]+exl[2]))*chi; 
  incr[5] = ((0.9682458365518543*phr[13]-0.9682458365518543*phl[13]-0.75*(phr[5]+phl[5])+0.4330127018922193*phr[3]-0.4330127018922193*phl[3])*c-0.9682458365518543*(exr[13]+exl[13])+0.75*exr[5]-0.75*exl[5]-0.4330127018922193*(exr[3]+exl[3]))*chi; 
  incr[6] = ((0.4330127018922193*(phr[10]+phl[10])-0.25*phr[6]+0.25*phl[6])*c-0.4330127018922193*exr[10]+0.4330127018922193*exl[10]+0.25*(exr[6]+exl[6]))*chi; 
  incr[7] = ((1.479019945774904*(phr[17]+phl[17])-1.25*phr[7]+1.25*phl[7]+0.9682458365518543*(phr[1]+phl[1])-0.5590169943749475*phr[0]+0.5590169943749475*phl[0])*c-1.479019945774904*exr[17]+1.479019945774904*exl[17]+1.25*(exr[7]+exl[7])-0.9682458365518543*exr[1]+0.9682458365518543*exl[1]+0.5590169943749475*(exr[0]+exl[0]))*chi; 
  incr[8] = ((0.4330127018922194*(phr[12]+phl[12])-0.25*phr[8]+0.25*phl[8])*c-0.4330127018922194*exr[12]+0.4330127018922194*exl[12]+0.25*(exr[8]+exl[8]))*chi; 
  incr[9] = ((0.4330127018922194*(phr[15]+phl[15])-0.25*phr[9]+0.25*phl[9])*c-0.4330127018922194*exr[15]+0.4330127018922194*exl[15]+0.25*(exr[9]+exl[9]))*chi; 
  incr[10] = (((-0.75*(phr[10]+phl[10]))+0.4330127018922193*phr[6]-0.4330127018922193*phl[6])*c+0.75*exr[10]-0.75*exl[10]-0.4330127018922193*(exr[6]+exl[6]))*chi; 
  incr[11] = (((-1.25*phr[11])+1.25*phl[11]+0.9682458365518543*(phr[4]+phl[4])-0.5590169943749476*phr[2]+0.5590169943749476*phl[2])*c+1.25*(exr[11]+exl[11])-0.9682458365518543*exr[4]+0.9682458365518543*exl[4]+0.5590169943749476*(exr[2]+exl[2]))*chi; 
  incr[12] = (((-0.75*(phr[12]+phl[12]))+0.4330127018922194*phr[8]-0.4330127018922194*phl[8])*c+0.75*exr[12]-0.75*exl[12]-0.4330127018922194*(exr[8]+exl[8]))*chi; 
  incr[13] = (((-1.25*phr[13])+1.25*phl[13]+0.9682458365518543*(phr[5]+phl[5])-0.5590169943749476*phr[3]+0.5590169943749476*phl[3])*c+1.25*(exr[13]+exl[13])-0.9682458365518543*exr[5]+0.9682458365518543*exl[5]+0.5590169943749476*(exr[3]+exl[3]))*chi; 
  incr[14] = (0.25*(phl[14]*c+exr[14]+exl[14])-0.25*phr[14]*c)*chi; 
  incr[15] = (((-0.75*(phr[15]+phl[15]))+0.4330127018922194*phr[9]-0.4330127018922194*phl[9])*c+0.75*exr[15]-0.75*exl[15]-0.4330127018922194*(exr[9]+exl[9]))*chi; 
  incr[16] = (0.25*(phl[16]*c+exr[16]+exl[16])-0.25*phr[16]*c)*chi; 
  incr[17] = (((-1.75*(phr[17]+phl[17]))+1.479019945774904*phr[7]-1.479019945774904*phl[7]-1.14564392373896*(phr[1]+phl[1])+0.6614378277661477*phr[0]-0.6614378277661477*phl[0])*c+1.75*exr[17]-1.75*exl[17]-1.479019945774904*(exr[7]+exl[7])+1.14564392373896*exr[1]-1.14564392373896*exl[1]-0.6614378277661477*(exr[0]+exl[0]))*chi; 
  incr[18] = (0.25*(phl[18]*c+exr[18]+exl[18])-0.25*phr[18]*c)*chi; 
  incr[19] = (0.25*(phl[19]*c+exr[19]+exl[19])-0.25*phr[19]*c)*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 
  outPhr[4] += incr[4]*dxr1; 
  outPhr[5] += incr[5]*dxr1; 
  outPhr[6] += incr[6]*dxr1; 
  outPhr[7] += incr[7]*dxr1; 
  outPhr[8] += incr[8]*dxr1; 
  outPhr[9] += incr[9]*dxr1; 
  outPhr[10] += incr[10]*dxr1; 
  outPhr[11] += incr[11]*dxr1; 
  outPhr[12] += incr[12]*dxr1; 
  outPhr[13] += incr[13]*dxr1; 
  outPhr[14] += incr[14]*dxr1; 
  outPhr[15] += incr[15]*dxr1; 
  outPhr[16] += incr[16]*dxr1; 
  outPhr[17] += incr[17]*dxr1; 
  outPhr[18] += incr[18]*dxr1; 
  outPhr[19] += incr[19]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 
  outPhl[3] += -1.0*incr[3]*dxl1; 
  outPhl[4] += incr[4]*dxl1; 
  outPhl[5] += incr[5]*dxl1; 
  outPhl[6] += -1.0*incr[6]*dxl1; 
  outPhl[7] += -1.0*incr[7]*dxl1; 
  outPhl[8] += -1.0*incr[8]*dxl1; 
  outPhl[9] += -1.0*incr[9]*dxl1; 
  outPhl[10] += incr[10]*dxl1; 
  outPhl[11] += -1.0*incr[11]*dxl1; 
  outPhl[12] += incr[12]*dxl1; 
  outPhl[13] += -1.0*incr[13]*dxl1; 
  outPhl[14] += -1.0*incr[14]*dxl1; 
  outPhl[15] += incr[15]*dxl1; 
  outPhl[16] += -1.0*incr[16]*dxl1; 
  outPhl[17] += incr[17]*dxl1; 
  outPhl[18] += -1.0*incr[18]*dxl1; 
  outPhl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = (0.6614378277661477*(psr[17]+psl[17])-0.5590169943749475*psr[7]+0.5590169943749475*psl[7]+0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.6614378277661477*bxr[17])+0.6614378277661477*bxl[17]+0.5590169943749475*(bxr[7]+bxl[7])-0.4330127018922193*bxr[1]+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = ((-1.14564392373896*(psr[17]+psl[17]))+0.9682458365518543*psr[7]-0.9682458365518543*psl[7]-0.75*(psr[1]+psl[1])+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(1.14564392373896*bxr[17]-1.14564392373896*bxl[17]-0.9682458365518543*(bxr[7]+bxl[7])+0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 
  incr[2] = ((-0.5590169943749476*psr[11])+0.5590169943749476*psl[11]+0.4330127018922193*(psr[4]+psl[4])-0.25*psr[2]+0.25*psl[2])*c*gamma+(0.5590169943749476*(bxr[11]+bxl[11])-0.4330127018922193*bxr[4]+0.4330127018922193*bxl[4]+0.25*(bxr[2]+bxl[2]))*c2gamma; 
  incr[3] = ((-0.5590169943749476*psr[13])+0.5590169943749476*psl[13]+0.4330127018922193*(psr[5]+psl[5])-0.25*psr[3]+0.25*psl[3])*c*gamma+(0.5590169943749476*(bxr[13]+bxl[13])-0.4330127018922193*bxr[5]+0.4330127018922193*bxl[5]+0.25*(bxr[3]+bxl[3]))*c2gamma; 
  incr[4] = (0.9682458365518543*psr[11]-0.9682458365518543*psl[11]-0.75*(psr[4]+psl[4])+0.4330127018922193*psr[2]-0.4330127018922193*psl[2])*c*gamma+((-0.9682458365518543*(bxr[11]+bxl[11]))+0.75*bxr[4]-0.75*bxl[4]-0.4330127018922193*(bxr[2]+bxl[2]))*c2gamma; 
  incr[5] = (0.9682458365518543*psr[13]-0.9682458365518543*psl[13]-0.75*(psr[5]+psl[5])+0.4330127018922193*psr[3]-0.4330127018922193*psl[3])*c*gamma+((-0.9682458365518543*(bxr[13]+bxl[13]))+0.75*bxr[5]-0.75*bxl[5]-0.4330127018922193*(bxr[3]+bxl[3]))*c2gamma; 
  incr[6] = (0.4330127018922193*(psr[10]+psl[10])-0.25*psr[6]+0.25*psl[6])*c*gamma+((-0.4330127018922193*bxr[10])+0.4330127018922193*bxl[10]+0.25*(bxr[6]+bxl[6]))*c2gamma; 
  incr[7] = (1.479019945774904*(psr[17]+psl[17])-1.25*psr[7]+1.25*psl[7]+0.9682458365518543*(psr[1]+psl[1])-0.5590169943749475*psr[0]+0.5590169943749475*psl[0])*c*gamma+((-1.479019945774904*bxr[17])+1.479019945774904*bxl[17]+1.25*(bxr[7]+bxl[7])-0.9682458365518543*bxr[1]+0.9682458365518543*bxl[1]+0.5590169943749475*(bxr[0]+bxl[0]))*c2gamma; 
  incr[8] = (0.4330127018922194*(psr[12]+psl[12])-0.25*psr[8]+0.25*psl[8])*c*gamma+((-0.4330127018922194*bxr[12])+0.4330127018922194*bxl[12]+0.25*(bxr[8]+bxl[8]))*c2gamma; 
  incr[9] = (0.4330127018922194*(psr[15]+psl[15])-0.25*psr[9]+0.25*psl[9])*c*gamma+((-0.4330127018922194*bxr[15])+0.4330127018922194*bxl[15]+0.25*(bxr[9]+bxl[9]))*c2gamma; 
  incr[10] = ((-0.75*(psr[10]+psl[10]))+0.4330127018922193*psr[6]-0.4330127018922193*psl[6])*c*gamma+(0.75*bxr[10]-0.75*bxl[10]-0.4330127018922193*(bxr[6]+bxl[6]))*c2gamma; 
  incr[11] = ((-1.25*psr[11])+1.25*psl[11]+0.9682458365518543*(psr[4]+psl[4])-0.5590169943749476*psr[2]+0.5590169943749476*psl[2])*c*gamma+(1.25*(bxr[11]+bxl[11])-0.9682458365518543*bxr[4]+0.9682458365518543*bxl[4]+0.5590169943749476*(bxr[2]+bxl[2]))*c2gamma; 
  incr[12] = ((-0.75*(psr[12]+psl[12]))+0.4330127018922194*psr[8]-0.4330127018922194*psl[8])*c*gamma+(0.75*bxr[12]-0.75*bxl[12]-0.4330127018922194*(bxr[8]+bxl[8]))*c2gamma; 
  incr[13] = ((-1.25*psr[13])+1.25*psl[13]+0.9682458365518543*(psr[5]+psl[5])-0.5590169943749476*psr[3]+0.5590169943749476*psl[3])*c*gamma+(1.25*(bxr[13]+bxl[13])-0.9682458365518543*bxr[5]+0.9682458365518543*bxl[5]+0.5590169943749476*(bxr[3]+bxl[3]))*c2gamma; 
  incr[14] = (0.25*psl[14]-0.25*psr[14])*c*gamma+0.25*(bxr[14]+bxl[14])*c2gamma; 
  incr[15] = ((-0.75*(psr[15]+psl[15]))+0.4330127018922194*psr[9]-0.4330127018922194*psl[9])*c*gamma+(0.75*bxr[15]-0.75*bxl[15]-0.4330127018922194*(bxr[9]+bxl[9]))*c2gamma; 
  incr[16] = (0.25*psl[16]-0.25*psr[16])*c*gamma+0.25*(bxr[16]+bxl[16])*c2gamma; 
  incr[17] = ((-1.75*(psr[17]+psl[17]))+1.479019945774904*psr[7]-1.479019945774904*psl[7]-1.14564392373896*(psr[1]+psl[1])+0.6614378277661477*psr[0]-0.6614378277661477*psl[0])*c*gamma+(1.75*bxr[17]-1.75*bxl[17]-1.479019945774904*(bxr[7]+bxl[7])+1.14564392373896*bxr[1]-1.14564392373896*bxl[1]-0.6614378277661477*(bxr[0]+bxl[0]))*c2gamma; 
  incr[18] = (0.25*psl[18]-0.25*psr[18])*c*gamma+0.25*(bxr[18]+bxl[18])*c2gamma; 
  incr[19] = (0.25*psl[19]-0.25*psr[19])*c*gamma+0.25*(bxr[19]+bxl[19])*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 
  outPsr[4] += incr[4]*dxr1; 
  outPsr[5] += incr[5]*dxr1; 
  outPsr[6] += incr[6]*dxr1; 
  outPsr[7] += incr[7]*dxr1; 
  outPsr[8] += incr[8]*dxr1; 
  outPsr[9] += incr[9]*dxr1; 
  outPsr[10] += incr[10]*dxr1; 
  outPsr[11] += incr[11]*dxr1; 
  outPsr[12] += incr[12]*dxr1; 
  outPsr[13] += incr[13]*dxr1; 
  outPsr[14] += incr[14]*dxr1; 
  outPsr[15] += incr[15]*dxr1; 
  outPsr[16] += incr[16]*dxr1; 
  outPsr[17] += incr[17]*dxr1; 
  outPsr[18] += incr[18]*dxr1; 
  outPsr[19] += incr[19]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 
  outPsl[3] += -1.0*incr[3]*dxl1; 
  outPsl[4] += incr[4]*dxl1; 
  outPsl[5] += incr[5]*dxl1; 
  outPsl[6] += -1.0*incr[6]*dxl1; 
  outPsl[7] += -1.0*incr[7]*dxl1; 
  outPsl[8] += -1.0*incr[8]*dxl1; 
  outPsl[9] += -1.0*incr[9]*dxl1; 
  outPsl[10] += incr[10]*dxl1; 
  outPsl[11] += -1.0*incr[11]*dxl1; 
  outPsl[12] += incr[12]*dxl1; 
  outPsl[13] += -1.0*incr[13]*dxl1; 
  outPsl[14] += -1.0*incr[14]*dxl1; 
  outPsl[15] += incr[15]*dxl1; 
  outPsl[16] += -1.0*incr[16]*dxl1; 
  outPsl[17] += incr[17]*dxl1; 
  outPsl[18] += -1.0*incr[18]*dxl1; 
  outPsl[19] += -1.0*incr[19]*dxl1; 

 
  return c; 
} 
double MaxwellSurf3xMax_Y_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[1]; 
  const double dxr1 = 2.0/dxr[1]; 
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
 
  incr[0] = (0.4330127018922193*(exr[2]+exl[2])-0.25*exr[0]+0.25*exl[0])*tau+(0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2]-0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = (0.25*exl[1]-0.25*exr[1])*tau-0.25*(bzr[1]+bzl[1])*c2; 
  incr[2] = ((-0.75*(exr[2]+exl[2]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*tau+((-0.75*bzr[2])+0.75*bzl[2]+0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[3] = (0.25*exl[3]-0.25*exr[3])*tau-0.25*(bzr[3]+bzl[3])*c2; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += -1.0*incr[1]*dxl1; 
  outExl[2] += incr[2]*dxl1; 
  outExl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = (0.4330127018922193*(eyr[2]+eyl[2])-0.25*eyr[0]+0.25*eyl[0])*c*chi+((-0.4330127018922193*phr[2])+0.4330127018922193*phl[2]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.25*eyl[1]-0.25*eyr[1])*c*chi+0.25*(phr[1]+phl[1])*c2chi; 
  incr[2] = ((-0.75*(eyr[2]+eyl[2]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*c*chi+(0.75*phr[2]-0.75*phl[2]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[3] = (0.25*eyl[3]-0.25*eyr[3])*c*chi+0.25*(phr[3]+phl[3])*c2chi; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += -1.0*incr[1]*dxl1; 
  outEyl[2] += incr[2]*dxl1; 
  outEyl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = (0.4330127018922193*(ezr[2]+ezl[2])-0.25*ezr[0]+0.25*ezl[0])*tau+((-0.4330127018922193*bxr[2])+0.4330127018922193*bxl[2]+0.25*(bxr[0]+bxl[0]))*c2; 
  incr[1] = 0.25*(ezl[1]*tau+(bxr[1]+bxl[1])*c2)-0.25*ezr[1]*tau; 
  incr[2] = ((-0.75*(ezr[2]+ezl[2]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+(0.75*bxr[2]-0.75*bxl[2]-0.4330127018922193*(bxr[0]+bxl[0]))*c2; 
  incr[3] = 0.25*(ezl[3]*tau+(bxr[3]+bxl[3])*c2)-0.25*ezr[3]*tau; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += -1.0*incr[1]*dxl1; 
  outEzl[2] += incr[2]*dxl1; 
  outEzl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bxr[2]+bxl[2])-0.25*bxr[0]+0.25*bxl[0])*c2)/tau-0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]+0.25*(ezr[0]+ezl[0]); 
  incr[1] = ((0.25*bxl[1]-0.25*bxr[1])*c2)/tau+0.25*(ezr[1]+ezl[1]); 
  incr[2] = (((-0.75*(bxr[2]+bxl[2]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c2)/tau+0.75*ezr[2]-0.75*ezl[2]-0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[3] = ((0.25*bxl[3]-0.25*bxr[3])*c2)/tau+0.25*(ezr[3]+ezl[3]); 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += -1.0*incr[1]*dxl1; 
  outBxl[2] += incr[2]*dxl1; 
  outBxl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(byr[2]+byl[2])-0.25*byr[0]+0.25*byl[0])*c-0.4330127018922193*psr[2]+0.4330127018922193*psl[2]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (0.25*(byl[1]*c+psr[1]+psl[1])-0.25*byr[1]*c)*gamma; 
  incr[2] = (((-0.75*(byr[2]+byl[2]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c+0.75*psr[2]-0.75*psl[2]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[3] = (0.25*(byl[3]*c+psr[3]+psl[3])-0.25*byr[3]*c)*gamma; 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += -1.0*incr[1]*dxl1; 
  outByl[2] += incr[2]*dxl1; 
  outByl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bzr[2]+bzl[2])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau+0.4330127018922193*exr[2]-0.4330127018922193*exl[2]-0.25*(exr[0]+exl[0]); 
  incr[1] = ((0.25*bzl[1]-0.25*bzr[1])*c2)/tau-0.25*(exr[1]+exl[1]); 
  incr[2] = (((-0.75*(bzr[2]+bzl[2]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau-0.75*exr[2]+0.75*exl[2]+0.4330127018922193*(exr[0]+exl[0]); 
  incr[3] = ((0.25*bzl[3]-0.25*bzr[3])*c2)/tau-0.25*(exr[3]+exl[3]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += -1.0*incr[1]*dxl1; 
  outBzl[2] += incr[2]*dxl1; 
  outBzl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(phr[2]+phl[2])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*eyr[2]+0.4330127018922193*eyl[2]+0.25*(eyr[0]+eyl[0]))*chi; 
  incr[1] = (0.25*(phl[1]*c+eyr[1]+eyl[1])-0.25*phr[1]*c)*chi; 
  incr[2] = (((-0.75*(phr[2]+phl[2]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*eyr[2]-0.75*eyl[2]-0.4330127018922193*(eyr[0]+eyl[0]))*chi; 
  incr[3] = (0.25*(phl[3]*c+eyr[3]+eyl[3])-0.25*phr[3]*c)*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += -1.0*incr[1]*dxl1; 
  outPhl[2] += incr[2]*dxl1; 
  outPhl[3] += -1.0*incr[3]*dxl1; 

 
  incr[0] = (0.4330127018922193*(psr[2]+psl[2])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*byr[2])+0.4330127018922193*byl[2]+0.25*(byr[0]+byl[0]))*c2gamma; 
  incr[1] = (0.25*psl[1]-0.25*psr[1])*c*gamma+0.25*(byr[1]+byl[1])*c2gamma; 
  incr[2] = ((-0.75*(psr[2]+psl[2]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*byr[2]-0.75*byl[2]-0.4330127018922193*(byr[0]+byl[0]))*c2gamma; 
  incr[3] = (0.25*psl[3]-0.25*psr[3])*c*gamma+0.25*(byr[3]+byl[3])*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += -1.0*incr[1]*dxl1; 
  outPsl[2] += incr[2]*dxl1; 
  outPsl[3] += -1.0*incr[3]*dxl1; 

 
  return c; 
} 
double MaxwellSurf3xMax_Y_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[1]; 
  const double dxr1 = 2.0/dxr[1]; 
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
 
  incr[0] = ((-0.5590169943749475*exr[8])+0.5590169943749475*exl[8]+0.4330127018922193*(exr[2]+exl[2])-0.25*exr[0]+0.25*exl[0])*tau+((-0.5590169943749475*(bzr[8]+bzl[8]))+0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2]-0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = (0.4330127018922193*(exr[4]+exl[4])-0.25*exr[1]+0.25*exl[1])*tau+(0.4330127018922193*bzr[4]-0.4330127018922193*bzl[4]-0.25*(bzr[1]+bzl[1]))*c2; 
  incr[2] = (0.9682458365518543*exr[8]-0.9682458365518543*exl[8]-0.75*(exr[2]+exl[2])+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*tau+(0.9682458365518543*(bzr[8]+bzl[8])-0.75*bzr[2]+0.75*bzl[2]+0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[3] = (0.4330127018922193*(exr[6]+exl[6])-0.25*exr[3]+0.25*exl[3])*tau+(0.4330127018922193*bzr[6]-0.4330127018922193*bzl[6]-0.25*(bzr[3]+bzl[3]))*c2; 
  incr[4] = ((-0.75*(exr[4]+exl[4]))+0.4330127018922193*exr[1]-0.4330127018922193*exl[1])*tau+((-0.75*bzr[4])+0.75*bzl[4]+0.4330127018922193*(bzr[1]+bzl[1]))*c2; 
  incr[5] = (0.25*exl[5]-0.25*exr[5])*tau-0.25*(bzr[5]+bzl[5])*c2; 
  incr[6] = ((-0.75*(exr[6]+exl[6]))+0.4330127018922193*exr[3]-0.4330127018922193*exl[3])*tau+((-0.75*bzr[6])+0.75*bzl[6]+0.4330127018922193*(bzr[3]+bzl[3]))*c2; 
  incr[7] = (0.25*exl[7]-0.25*exr[7])*tau-0.25*(bzr[7]+bzl[7])*c2; 
  incr[8] = ((-1.25*exr[8])+1.25*exl[8]+0.9682458365518543*(exr[2]+exl[2])-0.5590169943749475*exr[0]+0.5590169943749475*exl[0])*tau+((-1.25*(bzr[8]+bzl[8]))+0.9682458365518543*bzr[2]-0.9682458365518543*bzl[2]-0.5590169943749475*(bzr[0]+bzl[0]))*c2; 
  incr[9] = (0.25*exl[9]-0.25*exr[9])*tau-0.25*(bzr[9]+bzl[9])*c2; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 
  outExr[4] += incr[4]*dxr1; 
  outExr[5] += incr[5]*dxr1; 
  outExr[6] += incr[6]*dxr1; 
  outExr[7] += incr[7]*dxr1; 
  outExr[8] += incr[8]*dxr1; 
  outExr[9] += incr[9]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += -1.0*incr[1]*dxl1; 
  outExl[2] += incr[2]*dxl1; 
  outExl[3] += -1.0*incr[3]*dxl1; 
  outExl[4] += incr[4]*dxl1; 
  outExl[5] += -1.0*incr[5]*dxl1; 
  outExl[6] += incr[6]*dxl1; 
  outExl[7] += -1.0*incr[7]*dxl1; 
  outExl[8] += -1.0*incr[8]*dxl1; 
  outExl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*eyr[8])+0.5590169943749475*eyl[8]+0.4330127018922193*(eyr[2]+eyl[2])-0.25*eyr[0]+0.25*eyl[0])*c*chi+(0.5590169943749475*(phr[8]+phl[8])-0.4330127018922193*phr[2]+0.4330127018922193*phl[2]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.4330127018922193*(eyr[4]+eyl[4])-0.25*eyr[1]+0.25*eyl[1])*c*chi+((-0.4330127018922193*phr[4])+0.4330127018922193*phl[4]+0.25*(phr[1]+phl[1]))*c2chi; 
  incr[2] = (0.9682458365518543*eyr[8]-0.9682458365518543*eyl[8]-0.75*(eyr[2]+eyl[2])+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*c*chi+((-0.9682458365518543*(phr[8]+phl[8]))+0.75*phr[2]-0.75*phl[2]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[3] = (0.4330127018922193*(eyr[6]+eyl[6])-0.25*eyr[3]+0.25*eyl[3])*c*chi+((-0.4330127018922193*phr[6])+0.4330127018922193*phl[6]+0.25*(phr[3]+phl[3]))*c2chi; 
  incr[4] = ((-0.75*(eyr[4]+eyl[4]))+0.4330127018922193*eyr[1]-0.4330127018922193*eyl[1])*c*chi+(0.75*phr[4]-0.75*phl[4]-0.4330127018922193*(phr[1]+phl[1]))*c2chi; 
  incr[5] = (0.25*eyl[5]-0.25*eyr[5])*c*chi+0.25*(phr[5]+phl[5])*c2chi; 
  incr[6] = ((-0.75*(eyr[6]+eyl[6]))+0.4330127018922193*eyr[3]-0.4330127018922193*eyl[3])*c*chi+(0.75*phr[6]-0.75*phl[6]-0.4330127018922193*(phr[3]+phl[3]))*c2chi; 
  incr[7] = (0.25*eyl[7]-0.25*eyr[7])*c*chi+0.25*(phr[7]+phl[7])*c2chi; 
  incr[8] = ((-1.25*eyr[8])+1.25*eyl[8]+0.9682458365518543*(eyr[2]+eyl[2])-0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0])*c*chi+(1.25*(phr[8]+phl[8])-0.9682458365518543*phr[2]+0.9682458365518543*phl[2]+0.5590169943749475*(phr[0]+phl[0]))*c2chi; 
  incr[9] = (0.25*eyl[9]-0.25*eyr[9])*c*chi+0.25*(phr[9]+phl[9])*c2chi; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 
  outEyr[4] += incr[4]*dxr1; 
  outEyr[5] += incr[5]*dxr1; 
  outEyr[6] += incr[6]*dxr1; 
  outEyr[7] += incr[7]*dxr1; 
  outEyr[8] += incr[8]*dxr1; 
  outEyr[9] += incr[9]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += -1.0*incr[1]*dxl1; 
  outEyl[2] += incr[2]*dxl1; 
  outEyl[3] += -1.0*incr[3]*dxl1; 
  outEyl[4] += incr[4]*dxl1; 
  outEyl[5] += -1.0*incr[5]*dxl1; 
  outEyl[6] += incr[6]*dxl1; 
  outEyl[7] += -1.0*incr[7]*dxl1; 
  outEyl[8] += -1.0*incr[8]*dxl1; 
  outEyl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*ezr[8])+0.5590169943749475*ezl[8]+0.4330127018922193*(ezr[2]+ezl[2])-0.25*ezr[0]+0.25*ezl[0])*tau+(0.5590169943749475*(bxr[8]+bxl[8])-0.4330127018922193*bxr[2]+0.4330127018922193*bxl[2]+0.25*(bxr[0]+bxl[0]))*c2; 
  incr[1] = (0.4330127018922193*(ezr[4]+ezl[4])-0.25*ezr[1]+0.25*ezl[1])*tau+((-0.4330127018922193*bxr[4])+0.4330127018922193*bxl[4]+0.25*(bxr[1]+bxl[1]))*c2; 
  incr[2] = (0.9682458365518543*ezr[8]-0.9682458365518543*ezl[8]-0.75*(ezr[2]+ezl[2])+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+((-0.9682458365518543*(bxr[8]+bxl[8]))+0.75*bxr[2]-0.75*bxl[2]-0.4330127018922193*(bxr[0]+bxl[0]))*c2; 
  incr[3] = (0.4330127018922193*(ezr[6]+ezl[6])-0.25*ezr[3]+0.25*ezl[3])*tau+((-0.4330127018922193*bxr[6])+0.4330127018922193*bxl[6]+0.25*(bxr[3]+bxl[3]))*c2; 
  incr[4] = ((-0.75*(ezr[4]+ezl[4]))+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1])*tau+(0.75*bxr[4]-0.75*bxl[4]-0.4330127018922193*(bxr[1]+bxl[1]))*c2; 
  incr[5] = 0.25*(ezl[5]*tau+(bxr[5]+bxl[5])*c2)-0.25*ezr[5]*tau; 
  incr[6] = ((-0.75*(ezr[6]+ezl[6]))+0.4330127018922193*ezr[3]-0.4330127018922193*ezl[3])*tau+(0.75*bxr[6]-0.75*bxl[6]-0.4330127018922193*(bxr[3]+bxl[3]))*c2; 
  incr[7] = 0.25*(ezl[7]*tau+(bxr[7]+bxl[7])*c2)-0.25*ezr[7]*tau; 
  incr[8] = ((-1.25*ezr[8])+1.25*ezl[8]+0.9682458365518543*(ezr[2]+ezl[2])-0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0])*tau+(1.25*(bxr[8]+bxl[8])-0.9682458365518543*bxr[2]+0.9682458365518543*bxl[2]+0.5590169943749475*(bxr[0]+bxl[0]))*c2; 
  incr[9] = 0.25*(ezl[9]*tau+(bxr[9]+bxl[9])*c2)-0.25*ezr[9]*tau; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 
  outEzr[4] += incr[4]*dxr1; 
  outEzr[5] += incr[5]*dxr1; 
  outEzr[6] += incr[6]*dxr1; 
  outEzr[7] += incr[7]*dxr1; 
  outEzr[8] += incr[8]*dxr1; 
  outEzr[9] += incr[9]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += -1.0*incr[1]*dxl1; 
  outEzl[2] += incr[2]*dxl1; 
  outEzl[3] += -1.0*incr[3]*dxl1; 
  outEzl[4] += incr[4]*dxl1; 
  outEzl[5] += -1.0*incr[5]*dxl1; 
  outEzl[6] += incr[6]*dxl1; 
  outEzl[7] += -1.0*incr[7]*dxl1; 
  outEzl[8] += -1.0*incr[8]*dxl1; 
  outEzl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*bxr[8])+0.5590169943749475*bxl[8]+0.4330127018922193*(bxr[2]+bxl[2])-0.25*bxr[0]+0.25*bxl[0])*c2)/tau+0.5590169943749475*(ezr[8]+ezl[8])-0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]+0.25*(ezr[0]+ezl[0]); 
  incr[1] = ((0.4330127018922193*(bxr[4]+bxl[4])-0.25*bxr[1]+0.25*bxl[1])*c2)/tau-0.4330127018922193*ezr[4]+0.4330127018922193*ezl[4]+0.25*(ezr[1]+ezl[1]); 
  incr[2] = ((0.9682458365518543*bxr[8]-0.9682458365518543*bxl[8]-0.75*(bxr[2]+bxl[2])+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c2)/tau-0.9682458365518543*(ezr[8]+ezl[8])+0.75*ezr[2]-0.75*ezl[2]-0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[3] = ((0.4330127018922193*(bxr[6]+bxl[6])-0.25*bxr[3]+0.25*bxl[3])*c2)/tau-0.4330127018922193*ezr[6]+0.4330127018922193*ezl[6]+0.25*(ezr[3]+ezl[3]); 
  incr[4] = (((-0.75*(bxr[4]+bxl[4]))+0.4330127018922193*bxr[1]-0.4330127018922193*bxl[1])*c2)/tau+0.75*ezr[4]-0.75*ezl[4]-0.4330127018922193*(ezr[1]+ezl[1]); 
  incr[5] = ((0.25*bxl[5]-0.25*bxr[5])*c2)/tau+0.25*(ezr[5]+ezl[5]); 
  incr[6] = (((-0.75*(bxr[6]+bxl[6]))+0.4330127018922193*bxr[3]-0.4330127018922193*bxl[3])*c2)/tau+0.75*ezr[6]-0.75*ezl[6]-0.4330127018922193*(ezr[3]+ezl[3]); 
  incr[7] = ((0.25*bxl[7]-0.25*bxr[7])*c2)/tau+0.25*(ezr[7]+ezl[7]); 
  incr[8] = (((-1.25*bxr[8])+1.25*bxl[8]+0.9682458365518543*(bxr[2]+bxl[2])-0.5590169943749475*bxr[0]+0.5590169943749475*bxl[0])*c2)/tau+1.25*(ezr[8]+ezl[8])-0.9682458365518543*ezr[2]+0.9682458365518543*ezl[2]+0.5590169943749475*(ezr[0]+ezl[0]); 
  incr[9] = ((0.25*bxl[9]-0.25*bxr[9])*c2)/tau+0.25*(ezr[9]+ezl[9]); 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 
  outBxr[4] += incr[4]*dxr1; 
  outBxr[5] += incr[5]*dxr1; 
  outBxr[6] += incr[6]*dxr1; 
  outBxr[7] += incr[7]*dxr1; 
  outBxr[8] += incr[8]*dxr1; 
  outBxr[9] += incr[9]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += -1.0*incr[1]*dxl1; 
  outBxl[2] += incr[2]*dxl1; 
  outBxl[3] += -1.0*incr[3]*dxl1; 
  outBxl[4] += incr[4]*dxl1; 
  outBxl[5] += -1.0*incr[5]*dxl1; 
  outBxl[6] += incr[6]*dxl1; 
  outBxl[7] += -1.0*incr[7]*dxl1; 
  outBxl[8] += -1.0*incr[8]*dxl1; 
  outBxl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*byr[8])+0.5590169943749475*byl[8]+0.4330127018922193*(byr[2]+byl[2])-0.25*byr[0]+0.25*byl[0])*c+0.5590169943749475*(psr[8]+psl[8])-0.4330127018922193*psr[2]+0.4330127018922193*psl[2]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = ((0.4330127018922193*(byr[4]+byl[4])-0.25*byr[1]+0.25*byl[1])*c-0.4330127018922193*psr[4]+0.4330127018922193*psl[4]+0.25*(psr[1]+psl[1]))*gamma; 
  incr[2] = ((0.9682458365518543*byr[8]-0.9682458365518543*byl[8]-0.75*(byr[2]+byl[2])+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c-0.9682458365518543*(psr[8]+psl[8])+0.75*psr[2]-0.75*psl[2]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[3] = ((0.4330127018922193*(byr[6]+byl[6])-0.25*byr[3]+0.25*byl[3])*c-0.4330127018922193*psr[6]+0.4330127018922193*psl[6]+0.25*(psr[3]+psl[3]))*gamma; 
  incr[4] = (((-0.75*(byr[4]+byl[4]))+0.4330127018922193*byr[1]-0.4330127018922193*byl[1])*c+0.75*psr[4]-0.75*psl[4]-0.4330127018922193*(psr[1]+psl[1]))*gamma; 
  incr[5] = (0.25*(byl[5]*c+psr[5]+psl[5])-0.25*byr[5]*c)*gamma; 
  incr[6] = (((-0.75*(byr[6]+byl[6]))+0.4330127018922193*byr[3]-0.4330127018922193*byl[3])*c+0.75*psr[6]-0.75*psl[6]-0.4330127018922193*(psr[3]+psl[3]))*gamma; 
  incr[7] = (0.25*(byl[7]*c+psr[7]+psl[7])-0.25*byr[7]*c)*gamma; 
  incr[8] = (((-1.25*byr[8])+1.25*byl[8]+0.9682458365518543*(byr[2]+byl[2])-0.5590169943749475*byr[0]+0.5590169943749475*byl[0])*c+1.25*(psr[8]+psl[8])-0.9682458365518543*psr[2]+0.9682458365518543*psl[2]+0.5590169943749475*(psr[0]+psl[0]))*gamma; 
  incr[9] = (0.25*(byl[9]*c+psr[9]+psl[9])-0.25*byr[9]*c)*gamma; 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 
  outByr[4] += incr[4]*dxr1; 
  outByr[5] += incr[5]*dxr1; 
  outByr[6] += incr[6]*dxr1; 
  outByr[7] += incr[7]*dxr1; 
  outByr[8] += incr[8]*dxr1; 
  outByr[9] += incr[9]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += -1.0*incr[1]*dxl1; 
  outByl[2] += incr[2]*dxl1; 
  outByl[3] += -1.0*incr[3]*dxl1; 
  outByl[4] += incr[4]*dxl1; 
  outByl[5] += -1.0*incr[5]*dxl1; 
  outByl[6] += incr[6]*dxl1; 
  outByl[7] += -1.0*incr[7]*dxl1; 
  outByl[8] += -1.0*incr[8]*dxl1; 
  outByl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*bzr[8])+0.5590169943749475*bzl[8]+0.4330127018922193*(bzr[2]+bzl[2])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau-0.5590169943749475*(exr[8]+exl[8])+0.4330127018922193*exr[2]-0.4330127018922193*exl[2]-0.25*(exr[0]+exl[0]); 
  incr[1] = ((0.4330127018922193*(bzr[4]+bzl[4])-0.25*bzr[1]+0.25*bzl[1])*c2)/tau+0.4330127018922193*exr[4]-0.4330127018922193*exl[4]-0.25*(exr[1]+exl[1]); 
  incr[2] = ((0.9682458365518543*bzr[8]-0.9682458365518543*bzl[8]-0.75*(bzr[2]+bzl[2])+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau+0.9682458365518543*(exr[8]+exl[8])-0.75*exr[2]+0.75*exl[2]+0.4330127018922193*(exr[0]+exl[0]); 
  incr[3] = ((0.4330127018922193*(bzr[6]+bzl[6])-0.25*bzr[3]+0.25*bzl[3])*c2)/tau+0.4330127018922193*exr[6]-0.4330127018922193*exl[6]-0.25*(exr[3]+exl[3]); 
  incr[4] = (((-0.75*(bzr[4]+bzl[4]))+0.4330127018922193*bzr[1]-0.4330127018922193*bzl[1])*c2)/tau-0.75*exr[4]+0.75*exl[4]+0.4330127018922193*(exr[1]+exl[1]); 
  incr[5] = ((0.25*bzl[5]-0.25*bzr[5])*c2)/tau-0.25*(exr[5]+exl[5]); 
  incr[6] = (((-0.75*(bzr[6]+bzl[6]))+0.4330127018922193*bzr[3]-0.4330127018922193*bzl[3])*c2)/tau-0.75*exr[6]+0.75*exl[6]+0.4330127018922193*(exr[3]+exl[3]); 
  incr[7] = ((0.25*bzl[7]-0.25*bzr[7])*c2)/tau-0.25*(exr[7]+exl[7]); 
  incr[8] = (((-1.25*bzr[8])+1.25*bzl[8]+0.9682458365518543*(bzr[2]+bzl[2])-0.5590169943749475*bzr[0]+0.5590169943749475*bzl[0])*c2)/tau-1.25*(exr[8]+exl[8])+0.9682458365518543*exr[2]-0.9682458365518543*exl[2]-0.5590169943749475*(exr[0]+exl[0]); 
  incr[9] = ((0.25*bzl[9]-0.25*bzr[9])*c2)/tau-0.25*(exr[9]+exl[9]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 
  outBzr[4] += incr[4]*dxr1; 
  outBzr[5] += incr[5]*dxr1; 
  outBzr[6] += incr[6]*dxr1; 
  outBzr[7] += incr[7]*dxr1; 
  outBzr[8] += incr[8]*dxr1; 
  outBzr[9] += incr[9]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += -1.0*incr[1]*dxl1; 
  outBzl[2] += incr[2]*dxl1; 
  outBzl[3] += -1.0*incr[3]*dxl1; 
  outBzl[4] += incr[4]*dxl1; 
  outBzl[5] += -1.0*incr[5]*dxl1; 
  outBzl[6] += incr[6]*dxl1; 
  outBzl[7] += -1.0*incr[7]*dxl1; 
  outBzl[8] += -1.0*incr[8]*dxl1; 
  outBzl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*phr[8])+0.5590169943749475*phl[8]+0.4330127018922193*(phr[2]+phl[2])-0.25*phr[0]+0.25*phl[0])*c+0.5590169943749475*(eyr[8]+eyl[8])-0.4330127018922193*eyr[2]+0.4330127018922193*eyl[2]+0.25*(eyr[0]+eyl[0]))*chi; 
  incr[1] = ((0.4330127018922193*(phr[4]+phl[4])-0.25*phr[1]+0.25*phl[1])*c-0.4330127018922193*eyr[4]+0.4330127018922193*eyl[4]+0.25*(eyr[1]+eyl[1]))*chi; 
  incr[2] = ((0.9682458365518543*phr[8]-0.9682458365518543*phl[8]-0.75*(phr[2]+phl[2])+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c-0.9682458365518543*(eyr[8]+eyl[8])+0.75*eyr[2]-0.75*eyl[2]-0.4330127018922193*(eyr[0]+eyl[0]))*chi; 
  incr[3] = ((0.4330127018922193*(phr[6]+phl[6])-0.25*phr[3]+0.25*phl[3])*c-0.4330127018922193*eyr[6]+0.4330127018922193*eyl[6]+0.25*(eyr[3]+eyl[3]))*chi; 
  incr[4] = (((-0.75*(phr[4]+phl[4]))+0.4330127018922193*phr[1]-0.4330127018922193*phl[1])*c+0.75*eyr[4]-0.75*eyl[4]-0.4330127018922193*(eyr[1]+eyl[1]))*chi; 
  incr[5] = (0.25*(phl[5]*c+eyr[5]+eyl[5])-0.25*phr[5]*c)*chi; 
  incr[6] = (((-0.75*(phr[6]+phl[6]))+0.4330127018922193*phr[3]-0.4330127018922193*phl[3])*c+0.75*eyr[6]-0.75*eyl[6]-0.4330127018922193*(eyr[3]+eyl[3]))*chi; 
  incr[7] = (0.25*(phl[7]*c+eyr[7]+eyl[7])-0.25*phr[7]*c)*chi; 
  incr[8] = (((-1.25*phr[8])+1.25*phl[8]+0.9682458365518543*(phr[2]+phl[2])-0.5590169943749475*phr[0]+0.5590169943749475*phl[0])*c+1.25*(eyr[8]+eyl[8])-0.9682458365518543*eyr[2]+0.9682458365518543*eyl[2]+0.5590169943749475*(eyr[0]+eyl[0]))*chi; 
  incr[9] = (0.25*(phl[9]*c+eyr[9]+eyl[9])-0.25*phr[9]*c)*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 
  outPhr[4] += incr[4]*dxr1; 
  outPhr[5] += incr[5]*dxr1; 
  outPhr[6] += incr[6]*dxr1; 
  outPhr[7] += incr[7]*dxr1; 
  outPhr[8] += incr[8]*dxr1; 
  outPhr[9] += incr[9]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += -1.0*incr[1]*dxl1; 
  outPhl[2] += incr[2]*dxl1; 
  outPhl[3] += -1.0*incr[3]*dxl1; 
  outPhl[4] += incr[4]*dxl1; 
  outPhl[5] += -1.0*incr[5]*dxl1; 
  outPhl[6] += incr[6]*dxl1; 
  outPhl[7] += -1.0*incr[7]*dxl1; 
  outPhl[8] += -1.0*incr[8]*dxl1; 
  outPhl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*psr[8])+0.5590169943749475*psl[8]+0.4330127018922193*(psr[2]+psl[2])-0.25*psr[0]+0.25*psl[0])*c*gamma+(0.5590169943749475*(byr[8]+byl[8])-0.4330127018922193*byr[2]+0.4330127018922193*byl[2]+0.25*(byr[0]+byl[0]))*c2gamma; 
  incr[1] = (0.4330127018922193*(psr[4]+psl[4])-0.25*psr[1]+0.25*psl[1])*c*gamma+((-0.4330127018922193*byr[4])+0.4330127018922193*byl[4]+0.25*(byr[1]+byl[1]))*c2gamma; 
  incr[2] = (0.9682458365518543*psr[8]-0.9682458365518543*psl[8]-0.75*(psr[2]+psl[2])+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+((-0.9682458365518543*(byr[8]+byl[8]))+0.75*byr[2]-0.75*byl[2]-0.4330127018922193*(byr[0]+byl[0]))*c2gamma; 
  incr[3] = (0.4330127018922193*(psr[6]+psl[6])-0.25*psr[3]+0.25*psl[3])*c*gamma+((-0.4330127018922193*byr[6])+0.4330127018922193*byl[6]+0.25*(byr[3]+byl[3]))*c2gamma; 
  incr[4] = ((-0.75*(psr[4]+psl[4]))+0.4330127018922193*psr[1]-0.4330127018922193*psl[1])*c*gamma+(0.75*byr[4]-0.75*byl[4]-0.4330127018922193*(byr[1]+byl[1]))*c2gamma; 
  incr[5] = (0.25*psl[5]-0.25*psr[5])*c*gamma+0.25*(byr[5]+byl[5])*c2gamma; 
  incr[6] = ((-0.75*(psr[6]+psl[6]))+0.4330127018922193*psr[3]-0.4330127018922193*psl[3])*c*gamma+(0.75*byr[6]-0.75*byl[6]-0.4330127018922193*(byr[3]+byl[3]))*c2gamma; 
  incr[7] = (0.25*psl[7]-0.25*psr[7])*c*gamma+0.25*(byr[7]+byl[7])*c2gamma; 
  incr[8] = ((-1.25*psr[8])+1.25*psl[8]+0.9682458365518543*(psr[2]+psl[2])-0.5590169943749475*psr[0]+0.5590169943749475*psl[0])*c*gamma+(1.25*(byr[8]+byl[8])-0.9682458365518543*byr[2]+0.9682458365518543*byl[2]+0.5590169943749475*(byr[0]+byl[0]))*c2gamma; 
  incr[9] = (0.25*psl[9]-0.25*psr[9])*c*gamma+0.25*(byr[9]+byl[9])*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 
  outPsr[4] += incr[4]*dxr1; 
  outPsr[5] += incr[5]*dxr1; 
  outPsr[6] += incr[6]*dxr1; 
  outPsr[7] += incr[7]*dxr1; 
  outPsr[8] += incr[8]*dxr1; 
  outPsr[9] += incr[9]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += -1.0*incr[1]*dxl1; 
  outPsl[2] += incr[2]*dxl1; 
  outPsl[3] += -1.0*incr[3]*dxl1; 
  outPsl[4] += incr[4]*dxl1; 
  outPsl[5] += -1.0*incr[5]*dxl1; 
  outPsl[6] += incr[6]*dxl1; 
  outPsl[7] += -1.0*incr[7]*dxl1; 
  outPsl[8] += -1.0*incr[8]*dxl1; 
  outPsl[9] += -1.0*incr[9]*dxl1; 

 
  return c; 
} 
double MaxwellSurf3xMax_Y_P3(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[1]; 
  const double dxr1 = 2.0/dxr[1]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[20]; 
  const double *ezl = &ql[40]; 
  const double *bxl = &ql[60]; 
  const double *byl = &ql[80]; 
  const double *bzl = &ql[100]; 
  const double *phl = &ql[120]; 
  const double *psl = &ql[140]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[20]; 
  double *outEzl = &outl[40]; 
  double *outBxl = &outl[60]; 
  double *outByl = &outl[80]; 
  double *outBzl = &outl[100]; 
  double *outPhl = &outl[120]; 
  double *outPsl = &outl[140]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[20]; 
  const double *ezr = &qr[40]; 
  const double *bxr = &qr[60]; 
  const double *byr = &qr[80]; 
  const double *bzr = &qr[100]; 
  const double *phr = &qr[120]; 
  const double *psr = &qr[140]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[20]; 
  double *outEzr = &outr[40]; 
  double *outBxr = &outr[60]; 
  double *outByr = &outr[80]; 
  double *outBzr = &outr[100]; 
  double *outPhr = &outr[120]; 
  double *outPsr = &outr[140]; 
 
  double incr[20]; 
 
  incr[0] = (0.6614378277661477*(exr[18]+exl[18])-0.5590169943749475*exr[8]+0.5590169943749475*exl[8]+0.4330127018922193*(exr[2]+exl[2])-0.25*exr[0]+0.25*exl[0])*tau+(0.6614378277661477*bzr[18]-0.6614378277661477*bzl[18]-0.5590169943749475*(bzr[8]+bzl[8])+0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2]-0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = ((-0.5590169943749476*exr[12])+0.5590169943749476*exl[12]+0.4330127018922193*(exr[4]+exl[4])-0.25*exr[1]+0.25*exl[1])*tau+((-0.5590169943749476*(bzr[12]+bzl[12]))+0.4330127018922193*bzr[4]-0.4330127018922193*bzl[4]-0.25*(bzr[1]+bzl[1]))*c2; 
  incr[2] = ((-1.14564392373896*(exr[18]+exl[18]))+0.9682458365518543*exr[8]-0.9682458365518543*exl[8]-0.75*(exr[2]+exl[2])+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*tau+((-1.14564392373896*bzr[18])+1.14564392373896*bzl[18]+0.9682458365518543*(bzr[8]+bzl[8])-0.75*bzr[2]+0.75*bzl[2]+0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[3] = ((-0.5590169943749476*exr[14])+0.5590169943749476*exl[14]+0.4330127018922193*(exr[6]+exl[6])-0.25*exr[3]+0.25*exl[3])*tau+((-0.5590169943749476*(bzr[14]+bzl[14]))+0.4330127018922193*bzr[6]-0.4330127018922193*bzl[6]-0.25*(bzr[3]+bzl[3]))*c2; 
  incr[4] = (0.9682458365518543*exr[12]-0.9682458365518543*exl[12]-0.75*(exr[4]+exl[4])+0.4330127018922193*exr[1]-0.4330127018922193*exl[1])*tau+(0.9682458365518543*(bzr[12]+bzl[12])-0.75*bzr[4]+0.75*bzl[4]+0.4330127018922193*(bzr[1]+bzl[1]))*c2; 
  incr[5] = (0.4330127018922193*(exr[10]+exl[10])-0.25*exr[5]+0.25*exl[5])*tau+(0.4330127018922193*bzr[10]-0.4330127018922193*bzl[10]-0.25*(bzr[5]+bzl[5]))*c2; 
  incr[6] = (0.9682458365518543*exr[14]-0.9682458365518543*exl[14]-0.75*(exr[6]+exl[6])+0.4330127018922193*exr[3]-0.4330127018922193*exl[3])*tau+(0.9682458365518543*(bzr[14]+bzl[14])-0.75*bzr[6]+0.75*bzl[6]+0.4330127018922193*(bzr[3]+bzl[3]))*c2; 
  incr[7] = (0.4330127018922194*(exr[11]+exl[11])-0.25*exr[7]+0.25*exl[7])*tau+(0.4330127018922194*bzr[11]-0.4330127018922194*bzl[11]-0.25*(bzr[7]+bzl[7]))*c2; 
  incr[8] = (1.479019945774904*(exr[18]+exl[18])-1.25*exr[8]+1.25*exl[8]+0.9682458365518543*(exr[2]+exl[2])-0.5590169943749475*exr[0]+0.5590169943749475*exl[0])*tau+(1.479019945774904*bzr[18]-1.479019945774904*bzl[18]-1.25*(bzr[8]+bzl[8])+0.9682458365518543*bzr[2]-0.9682458365518543*bzl[2]-0.5590169943749475*(bzr[0]+bzl[0]))*c2; 
  incr[9] = (0.4330127018922194*(exr[16]+exl[16])-0.25*exr[9]+0.25*exl[9])*tau+(0.4330127018922194*bzr[16]-0.4330127018922194*bzl[16]-0.25*(bzr[9]+bzl[9]))*c2; 
  incr[10] = ((-0.75*(exr[10]+exl[10]))+0.4330127018922193*exr[5]-0.4330127018922193*exl[5])*tau+((-0.75*bzr[10])+0.75*bzl[10]+0.4330127018922193*(bzr[5]+bzl[5]))*c2; 
  incr[11] = ((-0.75*(exr[11]+exl[11]))+0.4330127018922194*exr[7]-0.4330127018922194*exl[7])*tau+((-0.75*bzr[11])+0.75*bzl[11]+0.4330127018922194*(bzr[7]+bzl[7]))*c2; 
  incr[12] = ((-1.25*exr[12])+1.25*exl[12]+0.9682458365518543*(exr[4]+exl[4])-0.5590169943749476*exr[1]+0.5590169943749476*exl[1])*tau+((-1.25*(bzr[12]+bzl[12]))+0.9682458365518543*bzr[4]-0.9682458365518543*bzl[4]-0.5590169943749476*(bzr[1]+bzl[1]))*c2; 
  incr[13] = (0.25*exl[13]-0.25*exr[13])*tau-0.25*(bzr[13]+bzl[13])*c2; 
  incr[14] = ((-1.25*exr[14])+1.25*exl[14]+0.9682458365518543*(exr[6]+exl[6])-0.5590169943749476*exr[3]+0.5590169943749476*exl[3])*tau+((-1.25*(bzr[14]+bzl[14]))+0.9682458365518543*bzr[6]-0.9682458365518543*bzl[6]-0.5590169943749476*(bzr[3]+bzl[3]))*c2; 
  incr[15] = (0.25*exl[15]-0.25*exr[15])*tau-0.25*(bzr[15]+bzl[15])*c2; 
  incr[16] = ((-0.75*(exr[16]+exl[16]))+0.4330127018922194*exr[9]-0.4330127018922194*exl[9])*tau+((-0.75*bzr[16])+0.75*bzl[16]+0.4330127018922194*(bzr[9]+bzl[9]))*c2; 
  incr[17] = (0.25*exl[17]-0.25*exr[17])*tau-0.25*(bzr[17]+bzl[17])*c2; 
  incr[18] = ((-1.75*(exr[18]+exl[18]))+1.479019945774904*exr[8]-1.479019945774904*exl[8]-1.14564392373896*(exr[2]+exl[2])+0.6614378277661477*exr[0]-0.6614378277661477*exl[0])*tau+((-1.75*bzr[18])+1.75*bzl[18]+1.479019945774904*(bzr[8]+bzl[8])-1.14564392373896*bzr[2]+1.14564392373896*bzl[2]+0.6614378277661477*(bzr[0]+bzl[0]))*c2; 
  incr[19] = (0.25*exl[19]-0.25*exr[19])*tau-0.25*(bzr[19]+bzl[19])*c2; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 
  outExr[4] += incr[4]*dxr1; 
  outExr[5] += incr[5]*dxr1; 
  outExr[6] += incr[6]*dxr1; 
  outExr[7] += incr[7]*dxr1; 
  outExr[8] += incr[8]*dxr1; 
  outExr[9] += incr[9]*dxr1; 
  outExr[10] += incr[10]*dxr1; 
  outExr[11] += incr[11]*dxr1; 
  outExr[12] += incr[12]*dxr1; 
  outExr[13] += incr[13]*dxr1; 
  outExr[14] += incr[14]*dxr1; 
  outExr[15] += incr[15]*dxr1; 
  outExr[16] += incr[16]*dxr1; 
  outExr[17] += incr[17]*dxr1; 
  outExr[18] += incr[18]*dxr1; 
  outExr[19] += incr[19]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += -1.0*incr[1]*dxl1; 
  outExl[2] += incr[2]*dxl1; 
  outExl[3] += -1.0*incr[3]*dxl1; 
  outExl[4] += incr[4]*dxl1; 
  outExl[5] += -1.0*incr[5]*dxl1; 
  outExl[6] += incr[6]*dxl1; 
  outExl[7] += -1.0*incr[7]*dxl1; 
  outExl[8] += -1.0*incr[8]*dxl1; 
  outExl[9] += -1.0*incr[9]*dxl1; 
  outExl[10] += incr[10]*dxl1; 
  outExl[11] += incr[11]*dxl1; 
  outExl[12] += -1.0*incr[12]*dxl1; 
  outExl[13] += -1.0*incr[13]*dxl1; 
  outExl[14] += -1.0*incr[14]*dxl1; 
  outExl[15] += -1.0*incr[15]*dxl1; 
  outExl[16] += incr[16]*dxl1; 
  outExl[17] += -1.0*incr[17]*dxl1; 
  outExl[18] += incr[18]*dxl1; 
  outExl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = (0.6614378277661477*(eyr[18]+eyl[18])-0.5590169943749475*eyr[8]+0.5590169943749475*eyl[8]+0.4330127018922193*(eyr[2]+eyl[2])-0.25*eyr[0]+0.25*eyl[0])*c*chi+((-0.6614378277661477*phr[18])+0.6614378277661477*phl[18]+0.5590169943749475*(phr[8]+phl[8])-0.4330127018922193*phr[2]+0.4330127018922193*phl[2]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = ((-0.5590169943749476*eyr[12])+0.5590169943749476*eyl[12]+0.4330127018922193*(eyr[4]+eyl[4])-0.25*eyr[1]+0.25*eyl[1])*c*chi+(0.5590169943749476*(phr[12]+phl[12])-0.4330127018922193*phr[4]+0.4330127018922193*phl[4]+0.25*(phr[1]+phl[1]))*c2chi; 
  incr[2] = ((-1.14564392373896*(eyr[18]+eyl[18]))+0.9682458365518543*eyr[8]-0.9682458365518543*eyl[8]-0.75*(eyr[2]+eyl[2])+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*c*chi+(1.14564392373896*phr[18]-1.14564392373896*phl[18]-0.9682458365518543*(phr[8]+phl[8])+0.75*phr[2]-0.75*phl[2]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[3] = ((-0.5590169943749476*eyr[14])+0.5590169943749476*eyl[14]+0.4330127018922193*(eyr[6]+eyl[6])-0.25*eyr[3]+0.25*eyl[3])*c*chi+(0.5590169943749476*(phr[14]+phl[14])-0.4330127018922193*phr[6]+0.4330127018922193*phl[6]+0.25*(phr[3]+phl[3]))*c2chi; 
  incr[4] = (0.9682458365518543*eyr[12]-0.9682458365518543*eyl[12]-0.75*(eyr[4]+eyl[4])+0.4330127018922193*eyr[1]-0.4330127018922193*eyl[1])*c*chi+((-0.9682458365518543*(phr[12]+phl[12]))+0.75*phr[4]-0.75*phl[4]-0.4330127018922193*(phr[1]+phl[1]))*c2chi; 
  incr[5] = (0.4330127018922193*(eyr[10]+eyl[10])-0.25*eyr[5]+0.25*eyl[5])*c*chi+((-0.4330127018922193*phr[10])+0.4330127018922193*phl[10]+0.25*(phr[5]+phl[5]))*c2chi; 
  incr[6] = (0.9682458365518543*eyr[14]-0.9682458365518543*eyl[14]-0.75*(eyr[6]+eyl[6])+0.4330127018922193*eyr[3]-0.4330127018922193*eyl[3])*c*chi+((-0.9682458365518543*(phr[14]+phl[14]))+0.75*phr[6]-0.75*phl[6]-0.4330127018922193*(phr[3]+phl[3]))*c2chi; 
  incr[7] = (0.4330127018922194*(eyr[11]+eyl[11])-0.25*eyr[7]+0.25*eyl[7])*c*chi+((-0.4330127018922194*phr[11])+0.4330127018922194*phl[11]+0.25*(phr[7]+phl[7]))*c2chi; 
  incr[8] = (1.479019945774904*(eyr[18]+eyl[18])-1.25*eyr[8]+1.25*eyl[8]+0.9682458365518543*(eyr[2]+eyl[2])-0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0])*c*chi+((-1.479019945774904*phr[18])+1.479019945774904*phl[18]+1.25*(phr[8]+phl[8])-0.9682458365518543*phr[2]+0.9682458365518543*phl[2]+0.5590169943749475*(phr[0]+phl[0]))*c2chi; 
  incr[9] = (0.4330127018922194*(eyr[16]+eyl[16])-0.25*eyr[9]+0.25*eyl[9])*c*chi+((-0.4330127018922194*phr[16])+0.4330127018922194*phl[16]+0.25*(phr[9]+phl[9]))*c2chi; 
  incr[10] = ((-0.75*(eyr[10]+eyl[10]))+0.4330127018922193*eyr[5]-0.4330127018922193*eyl[5])*c*chi+(0.75*phr[10]-0.75*phl[10]-0.4330127018922193*(phr[5]+phl[5]))*c2chi; 
  incr[11] = ((-0.75*(eyr[11]+eyl[11]))+0.4330127018922194*eyr[7]-0.4330127018922194*eyl[7])*c*chi+(0.75*phr[11]-0.75*phl[11]-0.4330127018922194*(phr[7]+phl[7]))*c2chi; 
  incr[12] = ((-1.25*eyr[12])+1.25*eyl[12]+0.9682458365518543*(eyr[4]+eyl[4])-0.5590169943749476*eyr[1]+0.5590169943749476*eyl[1])*c*chi+(1.25*(phr[12]+phl[12])-0.9682458365518543*phr[4]+0.9682458365518543*phl[4]+0.5590169943749476*(phr[1]+phl[1]))*c2chi; 
  incr[13] = (0.25*eyl[13]-0.25*eyr[13])*c*chi+0.25*(phr[13]+phl[13])*c2chi; 
  incr[14] = ((-1.25*eyr[14])+1.25*eyl[14]+0.9682458365518543*(eyr[6]+eyl[6])-0.5590169943749476*eyr[3]+0.5590169943749476*eyl[3])*c*chi+(1.25*(phr[14]+phl[14])-0.9682458365518543*phr[6]+0.9682458365518543*phl[6]+0.5590169943749476*(phr[3]+phl[3]))*c2chi; 
  incr[15] = (0.25*eyl[15]-0.25*eyr[15])*c*chi+0.25*(phr[15]+phl[15])*c2chi; 
  incr[16] = ((-0.75*(eyr[16]+eyl[16]))+0.4330127018922194*eyr[9]-0.4330127018922194*eyl[9])*c*chi+(0.75*phr[16]-0.75*phl[16]-0.4330127018922194*(phr[9]+phl[9]))*c2chi; 
  incr[17] = (0.25*eyl[17]-0.25*eyr[17])*c*chi+0.25*(phr[17]+phl[17])*c2chi; 
  incr[18] = ((-1.75*(eyr[18]+eyl[18]))+1.479019945774904*eyr[8]-1.479019945774904*eyl[8]-1.14564392373896*(eyr[2]+eyl[2])+0.6614378277661477*eyr[0]-0.6614378277661477*eyl[0])*c*chi+(1.75*phr[18]-1.75*phl[18]-1.479019945774904*(phr[8]+phl[8])+1.14564392373896*phr[2]-1.14564392373896*phl[2]-0.6614378277661477*(phr[0]+phl[0]))*c2chi; 
  incr[19] = (0.25*eyl[19]-0.25*eyr[19])*c*chi+0.25*(phr[19]+phl[19])*c2chi; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 
  outEyr[4] += incr[4]*dxr1; 
  outEyr[5] += incr[5]*dxr1; 
  outEyr[6] += incr[6]*dxr1; 
  outEyr[7] += incr[7]*dxr1; 
  outEyr[8] += incr[8]*dxr1; 
  outEyr[9] += incr[9]*dxr1; 
  outEyr[10] += incr[10]*dxr1; 
  outEyr[11] += incr[11]*dxr1; 
  outEyr[12] += incr[12]*dxr1; 
  outEyr[13] += incr[13]*dxr1; 
  outEyr[14] += incr[14]*dxr1; 
  outEyr[15] += incr[15]*dxr1; 
  outEyr[16] += incr[16]*dxr1; 
  outEyr[17] += incr[17]*dxr1; 
  outEyr[18] += incr[18]*dxr1; 
  outEyr[19] += incr[19]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += -1.0*incr[1]*dxl1; 
  outEyl[2] += incr[2]*dxl1; 
  outEyl[3] += -1.0*incr[3]*dxl1; 
  outEyl[4] += incr[4]*dxl1; 
  outEyl[5] += -1.0*incr[5]*dxl1; 
  outEyl[6] += incr[6]*dxl1; 
  outEyl[7] += -1.0*incr[7]*dxl1; 
  outEyl[8] += -1.0*incr[8]*dxl1; 
  outEyl[9] += -1.0*incr[9]*dxl1; 
  outEyl[10] += incr[10]*dxl1; 
  outEyl[11] += incr[11]*dxl1; 
  outEyl[12] += -1.0*incr[12]*dxl1; 
  outEyl[13] += -1.0*incr[13]*dxl1; 
  outEyl[14] += -1.0*incr[14]*dxl1; 
  outEyl[15] += -1.0*incr[15]*dxl1; 
  outEyl[16] += incr[16]*dxl1; 
  outEyl[17] += -1.0*incr[17]*dxl1; 
  outEyl[18] += incr[18]*dxl1; 
  outEyl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = (0.6614378277661477*(ezr[18]+ezl[18])-0.5590169943749475*ezr[8]+0.5590169943749475*ezl[8]+0.4330127018922193*(ezr[2]+ezl[2])-0.25*ezr[0]+0.25*ezl[0])*tau+((-0.6614378277661477*bxr[18])+0.6614378277661477*bxl[18]+0.5590169943749475*(bxr[8]+bxl[8])-0.4330127018922193*bxr[2]+0.4330127018922193*bxl[2]+0.25*(bxr[0]+bxl[0]))*c2; 
  incr[1] = ((-0.5590169943749476*ezr[12])+0.5590169943749476*ezl[12]+0.4330127018922193*(ezr[4]+ezl[4])-0.25*ezr[1]+0.25*ezl[1])*tau+(0.5590169943749476*(bxr[12]+bxl[12])-0.4330127018922193*bxr[4]+0.4330127018922193*bxl[4]+0.25*(bxr[1]+bxl[1]))*c2; 
  incr[2] = ((-1.14564392373896*(ezr[18]+ezl[18]))+0.9682458365518543*ezr[8]-0.9682458365518543*ezl[8]-0.75*(ezr[2]+ezl[2])+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+(1.14564392373896*bxr[18]-1.14564392373896*bxl[18]-0.9682458365518543*(bxr[8]+bxl[8])+0.75*bxr[2]-0.75*bxl[2]-0.4330127018922193*(bxr[0]+bxl[0]))*c2; 
  incr[3] = ((-0.5590169943749476*ezr[14])+0.5590169943749476*ezl[14]+0.4330127018922193*(ezr[6]+ezl[6])-0.25*ezr[3]+0.25*ezl[3])*tau+(0.5590169943749476*(bxr[14]+bxl[14])-0.4330127018922193*bxr[6]+0.4330127018922193*bxl[6]+0.25*(bxr[3]+bxl[3]))*c2; 
  incr[4] = (0.9682458365518543*ezr[12]-0.9682458365518543*ezl[12]-0.75*(ezr[4]+ezl[4])+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1])*tau+((-0.9682458365518543*(bxr[12]+bxl[12]))+0.75*bxr[4]-0.75*bxl[4]-0.4330127018922193*(bxr[1]+bxl[1]))*c2; 
  incr[5] = (0.4330127018922193*(ezr[10]+ezl[10])-0.25*ezr[5]+0.25*ezl[5])*tau+((-0.4330127018922193*bxr[10])+0.4330127018922193*bxl[10]+0.25*(bxr[5]+bxl[5]))*c2; 
  incr[6] = (0.9682458365518543*ezr[14]-0.9682458365518543*ezl[14]-0.75*(ezr[6]+ezl[6])+0.4330127018922193*ezr[3]-0.4330127018922193*ezl[3])*tau+((-0.9682458365518543*(bxr[14]+bxl[14]))+0.75*bxr[6]-0.75*bxl[6]-0.4330127018922193*(bxr[3]+bxl[3]))*c2; 
  incr[7] = (0.4330127018922194*(ezr[11]+ezl[11])-0.25*ezr[7]+0.25*ezl[7])*tau+((-0.4330127018922194*bxr[11])+0.4330127018922194*bxl[11]+0.25*(bxr[7]+bxl[7]))*c2; 
  incr[8] = (1.479019945774904*(ezr[18]+ezl[18])-1.25*ezr[8]+1.25*ezl[8]+0.9682458365518543*(ezr[2]+ezl[2])-0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0])*tau+((-1.479019945774904*bxr[18])+1.479019945774904*bxl[18]+1.25*(bxr[8]+bxl[8])-0.9682458365518543*bxr[2]+0.9682458365518543*bxl[2]+0.5590169943749475*(bxr[0]+bxl[0]))*c2; 
  incr[9] = (0.4330127018922194*(ezr[16]+ezl[16])-0.25*ezr[9]+0.25*ezl[9])*tau+((-0.4330127018922194*bxr[16])+0.4330127018922194*bxl[16]+0.25*(bxr[9]+bxl[9]))*c2; 
  incr[10] = ((-0.75*(ezr[10]+ezl[10]))+0.4330127018922193*ezr[5]-0.4330127018922193*ezl[5])*tau+(0.75*bxr[10]-0.75*bxl[10]-0.4330127018922193*(bxr[5]+bxl[5]))*c2; 
  incr[11] = ((-0.75*(ezr[11]+ezl[11]))+0.4330127018922194*ezr[7]-0.4330127018922194*ezl[7])*tau+(0.75*bxr[11]-0.75*bxl[11]-0.4330127018922194*(bxr[7]+bxl[7]))*c2; 
  incr[12] = ((-1.25*ezr[12])+1.25*ezl[12]+0.9682458365518543*(ezr[4]+ezl[4])-0.5590169943749476*ezr[1]+0.5590169943749476*ezl[1])*tau+(1.25*(bxr[12]+bxl[12])-0.9682458365518543*bxr[4]+0.9682458365518543*bxl[4]+0.5590169943749476*(bxr[1]+bxl[1]))*c2; 
  incr[13] = 0.25*(ezl[13]*tau+(bxr[13]+bxl[13])*c2)-0.25*ezr[13]*tau; 
  incr[14] = ((-1.25*ezr[14])+1.25*ezl[14]+0.9682458365518543*(ezr[6]+ezl[6])-0.5590169943749476*ezr[3]+0.5590169943749476*ezl[3])*tau+(1.25*(bxr[14]+bxl[14])-0.9682458365518543*bxr[6]+0.9682458365518543*bxl[6]+0.5590169943749476*(bxr[3]+bxl[3]))*c2; 
  incr[15] = 0.25*(ezl[15]*tau+(bxr[15]+bxl[15])*c2)-0.25*ezr[15]*tau; 
  incr[16] = ((-0.75*(ezr[16]+ezl[16]))+0.4330127018922194*ezr[9]-0.4330127018922194*ezl[9])*tau+(0.75*bxr[16]-0.75*bxl[16]-0.4330127018922194*(bxr[9]+bxl[9]))*c2; 
  incr[17] = 0.25*(ezl[17]*tau+(bxr[17]+bxl[17])*c2)-0.25*ezr[17]*tau; 
  incr[18] = ((-1.75*(ezr[18]+ezl[18]))+1.479019945774904*ezr[8]-1.479019945774904*ezl[8]-1.14564392373896*(ezr[2]+ezl[2])+0.6614378277661477*ezr[0]-0.6614378277661477*ezl[0])*tau+(1.75*bxr[18]-1.75*bxl[18]-1.479019945774904*(bxr[8]+bxl[8])+1.14564392373896*bxr[2]-1.14564392373896*bxl[2]-0.6614378277661477*(bxr[0]+bxl[0]))*c2; 
  incr[19] = 0.25*(ezl[19]*tau+(bxr[19]+bxl[19])*c2)-0.25*ezr[19]*tau; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 
  outEzr[4] += incr[4]*dxr1; 
  outEzr[5] += incr[5]*dxr1; 
  outEzr[6] += incr[6]*dxr1; 
  outEzr[7] += incr[7]*dxr1; 
  outEzr[8] += incr[8]*dxr1; 
  outEzr[9] += incr[9]*dxr1; 
  outEzr[10] += incr[10]*dxr1; 
  outEzr[11] += incr[11]*dxr1; 
  outEzr[12] += incr[12]*dxr1; 
  outEzr[13] += incr[13]*dxr1; 
  outEzr[14] += incr[14]*dxr1; 
  outEzr[15] += incr[15]*dxr1; 
  outEzr[16] += incr[16]*dxr1; 
  outEzr[17] += incr[17]*dxr1; 
  outEzr[18] += incr[18]*dxr1; 
  outEzr[19] += incr[19]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += -1.0*incr[1]*dxl1; 
  outEzl[2] += incr[2]*dxl1; 
  outEzl[3] += -1.0*incr[3]*dxl1; 
  outEzl[4] += incr[4]*dxl1; 
  outEzl[5] += -1.0*incr[5]*dxl1; 
  outEzl[6] += incr[6]*dxl1; 
  outEzl[7] += -1.0*incr[7]*dxl1; 
  outEzl[8] += -1.0*incr[8]*dxl1; 
  outEzl[9] += -1.0*incr[9]*dxl1; 
  outEzl[10] += incr[10]*dxl1; 
  outEzl[11] += incr[11]*dxl1; 
  outEzl[12] += -1.0*incr[12]*dxl1; 
  outEzl[13] += -1.0*incr[13]*dxl1; 
  outEzl[14] += -1.0*incr[14]*dxl1; 
  outEzl[15] += -1.0*incr[15]*dxl1; 
  outEzl[16] += incr[16]*dxl1; 
  outEzl[17] += -1.0*incr[17]*dxl1; 
  outEzl[18] += incr[18]*dxl1; 
  outEzl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(bxr[18]+bxl[18])-0.5590169943749475*bxr[8]+0.5590169943749475*bxl[8]+0.4330127018922193*(bxr[2]+bxl[2])-0.25*bxr[0]+0.25*bxl[0])*c2)/tau-0.6614378277661477*ezr[18]+0.6614378277661477*ezl[18]+0.5590169943749475*(ezr[8]+ezl[8])-0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]+0.25*(ezr[0]+ezl[0]); 
  incr[1] = (((-0.5590169943749476*bxr[12])+0.5590169943749476*bxl[12]+0.4330127018922193*(bxr[4]+bxl[4])-0.25*bxr[1]+0.25*bxl[1])*c2)/tau+0.5590169943749476*(ezr[12]+ezl[12])-0.4330127018922193*ezr[4]+0.4330127018922193*ezl[4]+0.25*(ezr[1]+ezl[1]); 
  incr[2] = (((-1.14564392373896*(bxr[18]+bxl[18]))+0.9682458365518543*bxr[8]-0.9682458365518543*bxl[8]-0.75*(bxr[2]+bxl[2])+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c2)/tau+1.14564392373896*ezr[18]-1.14564392373896*ezl[18]-0.9682458365518543*(ezr[8]+ezl[8])+0.75*ezr[2]-0.75*ezl[2]-0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[3] = (((-0.5590169943749476*bxr[14])+0.5590169943749476*bxl[14]+0.4330127018922193*(bxr[6]+bxl[6])-0.25*bxr[3]+0.25*bxl[3])*c2)/tau+0.5590169943749476*(ezr[14]+ezl[14])-0.4330127018922193*ezr[6]+0.4330127018922193*ezl[6]+0.25*(ezr[3]+ezl[3]); 
  incr[4] = ((0.9682458365518543*bxr[12]-0.9682458365518543*bxl[12]-0.75*(bxr[4]+bxl[4])+0.4330127018922193*bxr[1]-0.4330127018922193*bxl[1])*c2)/tau-0.9682458365518543*(ezr[12]+ezl[12])+0.75*ezr[4]-0.75*ezl[4]-0.4330127018922193*(ezr[1]+ezl[1]); 
  incr[5] = ((0.4330127018922193*(bxr[10]+bxl[10])-0.25*bxr[5]+0.25*bxl[5])*c2)/tau-0.4330127018922193*ezr[10]+0.4330127018922193*ezl[10]+0.25*(ezr[5]+ezl[5]); 
  incr[6] = ((0.9682458365518543*bxr[14]-0.9682458365518543*bxl[14]-0.75*(bxr[6]+bxl[6])+0.4330127018922193*bxr[3]-0.4330127018922193*bxl[3])*c2)/tau-0.9682458365518543*(ezr[14]+ezl[14])+0.75*ezr[6]-0.75*ezl[6]-0.4330127018922193*(ezr[3]+ezl[3]); 
  incr[7] = ((0.4330127018922194*(bxr[11]+bxl[11])-0.25*bxr[7]+0.25*bxl[7])*c2)/tau-0.4330127018922194*ezr[11]+0.4330127018922194*ezl[11]+0.25*(ezr[7]+ezl[7]); 
  incr[8] = ((1.479019945774904*(bxr[18]+bxl[18])-1.25*bxr[8]+1.25*bxl[8]+0.9682458365518543*(bxr[2]+bxl[2])-0.5590169943749475*bxr[0]+0.5590169943749475*bxl[0])*c2)/tau-1.479019945774904*ezr[18]+1.479019945774904*ezl[18]+1.25*(ezr[8]+ezl[8])-0.9682458365518543*ezr[2]+0.9682458365518543*ezl[2]+0.5590169943749475*(ezr[0]+ezl[0]); 
  incr[9] = ((0.4330127018922194*(bxr[16]+bxl[16])-0.25*bxr[9]+0.25*bxl[9])*c2)/tau-0.4330127018922194*ezr[16]+0.4330127018922194*ezl[16]+0.25*(ezr[9]+ezl[9]); 
  incr[10] = (((-0.75*(bxr[10]+bxl[10]))+0.4330127018922193*bxr[5]-0.4330127018922193*bxl[5])*c2)/tau+0.75*ezr[10]-0.75*ezl[10]-0.4330127018922193*(ezr[5]+ezl[5]); 
  incr[11] = (((-0.75*(bxr[11]+bxl[11]))+0.4330127018922194*bxr[7]-0.4330127018922194*bxl[7])*c2)/tau+0.75*ezr[11]-0.75*ezl[11]-0.4330127018922194*(ezr[7]+ezl[7]); 
  incr[12] = (((-1.25*bxr[12])+1.25*bxl[12]+0.9682458365518543*(bxr[4]+bxl[4])-0.5590169943749476*bxr[1]+0.5590169943749476*bxl[1])*c2)/tau+1.25*(ezr[12]+ezl[12])-0.9682458365518543*ezr[4]+0.9682458365518543*ezl[4]+0.5590169943749476*(ezr[1]+ezl[1]); 
  incr[13] = ((0.25*bxl[13]-0.25*bxr[13])*c2)/tau+0.25*(ezr[13]+ezl[13]); 
  incr[14] = (((-1.25*bxr[14])+1.25*bxl[14]+0.9682458365518543*(bxr[6]+bxl[6])-0.5590169943749476*bxr[3]+0.5590169943749476*bxl[3])*c2)/tau+1.25*(ezr[14]+ezl[14])-0.9682458365518543*ezr[6]+0.9682458365518543*ezl[6]+0.5590169943749476*(ezr[3]+ezl[3]); 
  incr[15] = ((0.25*bxl[15]-0.25*bxr[15])*c2)/tau+0.25*(ezr[15]+ezl[15]); 
  incr[16] = (((-0.75*(bxr[16]+bxl[16]))+0.4330127018922194*bxr[9]-0.4330127018922194*bxl[9])*c2)/tau+0.75*ezr[16]-0.75*ezl[16]-0.4330127018922194*(ezr[9]+ezl[9]); 
  incr[17] = ((0.25*bxl[17]-0.25*bxr[17])*c2)/tau+0.25*(ezr[17]+ezl[17]); 
  incr[18] = (((-1.75*(bxr[18]+bxl[18]))+1.479019945774904*bxr[8]-1.479019945774904*bxl[8]-1.14564392373896*(bxr[2]+bxl[2])+0.6614378277661477*bxr[0]-0.6614378277661477*bxl[0])*c2)/tau+1.75*ezr[18]-1.75*ezl[18]-1.479019945774904*(ezr[8]+ezl[8])+1.14564392373896*ezr[2]-1.14564392373896*ezl[2]-0.6614378277661477*(ezr[0]+ezl[0]); 
  incr[19] = ((0.25*bxl[19]-0.25*bxr[19])*c2)/tau+0.25*(ezr[19]+ezl[19]); 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 
  outBxr[4] += incr[4]*dxr1; 
  outBxr[5] += incr[5]*dxr1; 
  outBxr[6] += incr[6]*dxr1; 
  outBxr[7] += incr[7]*dxr1; 
  outBxr[8] += incr[8]*dxr1; 
  outBxr[9] += incr[9]*dxr1; 
  outBxr[10] += incr[10]*dxr1; 
  outBxr[11] += incr[11]*dxr1; 
  outBxr[12] += incr[12]*dxr1; 
  outBxr[13] += incr[13]*dxr1; 
  outBxr[14] += incr[14]*dxr1; 
  outBxr[15] += incr[15]*dxr1; 
  outBxr[16] += incr[16]*dxr1; 
  outBxr[17] += incr[17]*dxr1; 
  outBxr[18] += incr[18]*dxr1; 
  outBxr[19] += incr[19]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += -1.0*incr[1]*dxl1; 
  outBxl[2] += incr[2]*dxl1; 
  outBxl[3] += -1.0*incr[3]*dxl1; 
  outBxl[4] += incr[4]*dxl1; 
  outBxl[5] += -1.0*incr[5]*dxl1; 
  outBxl[6] += incr[6]*dxl1; 
  outBxl[7] += -1.0*incr[7]*dxl1; 
  outBxl[8] += -1.0*incr[8]*dxl1; 
  outBxl[9] += -1.0*incr[9]*dxl1; 
  outBxl[10] += incr[10]*dxl1; 
  outBxl[11] += incr[11]*dxl1; 
  outBxl[12] += -1.0*incr[12]*dxl1; 
  outBxl[13] += -1.0*incr[13]*dxl1; 
  outBxl[14] += -1.0*incr[14]*dxl1; 
  outBxl[15] += -1.0*incr[15]*dxl1; 
  outBxl[16] += incr[16]*dxl1; 
  outBxl[17] += -1.0*incr[17]*dxl1; 
  outBxl[18] += incr[18]*dxl1; 
  outBxl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(byr[18]+byl[18])-0.5590169943749475*byr[8]+0.5590169943749475*byl[8]+0.4330127018922193*(byr[2]+byl[2])-0.25*byr[0]+0.25*byl[0])*c-0.6614378277661477*psr[18]+0.6614378277661477*psl[18]+0.5590169943749475*(psr[8]+psl[8])-0.4330127018922193*psr[2]+0.4330127018922193*psl[2]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (((-0.5590169943749476*byr[12])+0.5590169943749476*byl[12]+0.4330127018922193*(byr[4]+byl[4])-0.25*byr[1]+0.25*byl[1])*c+0.5590169943749476*(psr[12]+psl[12])-0.4330127018922193*psr[4]+0.4330127018922193*psl[4]+0.25*(psr[1]+psl[1]))*gamma; 
  incr[2] = (((-1.14564392373896*(byr[18]+byl[18]))+0.9682458365518543*byr[8]-0.9682458365518543*byl[8]-0.75*(byr[2]+byl[2])+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c+1.14564392373896*psr[18]-1.14564392373896*psl[18]-0.9682458365518543*(psr[8]+psl[8])+0.75*psr[2]-0.75*psl[2]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[3] = (((-0.5590169943749476*byr[14])+0.5590169943749476*byl[14]+0.4330127018922193*(byr[6]+byl[6])-0.25*byr[3]+0.25*byl[3])*c+0.5590169943749476*(psr[14]+psl[14])-0.4330127018922193*psr[6]+0.4330127018922193*psl[6]+0.25*(psr[3]+psl[3]))*gamma; 
  incr[4] = ((0.9682458365518543*byr[12]-0.9682458365518543*byl[12]-0.75*(byr[4]+byl[4])+0.4330127018922193*byr[1]-0.4330127018922193*byl[1])*c-0.9682458365518543*(psr[12]+psl[12])+0.75*psr[4]-0.75*psl[4]-0.4330127018922193*(psr[1]+psl[1]))*gamma; 
  incr[5] = ((0.4330127018922193*(byr[10]+byl[10])-0.25*byr[5]+0.25*byl[5])*c-0.4330127018922193*psr[10]+0.4330127018922193*psl[10]+0.25*(psr[5]+psl[5]))*gamma; 
  incr[6] = ((0.9682458365518543*byr[14]-0.9682458365518543*byl[14]-0.75*(byr[6]+byl[6])+0.4330127018922193*byr[3]-0.4330127018922193*byl[3])*c-0.9682458365518543*(psr[14]+psl[14])+0.75*psr[6]-0.75*psl[6]-0.4330127018922193*(psr[3]+psl[3]))*gamma; 
  incr[7] = ((0.4330127018922194*(byr[11]+byl[11])-0.25*byr[7]+0.25*byl[7])*c-0.4330127018922194*psr[11]+0.4330127018922194*psl[11]+0.25*(psr[7]+psl[7]))*gamma; 
  incr[8] = ((1.479019945774904*(byr[18]+byl[18])-1.25*byr[8]+1.25*byl[8]+0.9682458365518543*(byr[2]+byl[2])-0.5590169943749475*byr[0]+0.5590169943749475*byl[0])*c-1.479019945774904*psr[18]+1.479019945774904*psl[18]+1.25*(psr[8]+psl[8])-0.9682458365518543*psr[2]+0.9682458365518543*psl[2]+0.5590169943749475*(psr[0]+psl[0]))*gamma; 
  incr[9] = ((0.4330127018922194*(byr[16]+byl[16])-0.25*byr[9]+0.25*byl[9])*c-0.4330127018922194*psr[16]+0.4330127018922194*psl[16]+0.25*(psr[9]+psl[9]))*gamma; 
  incr[10] = (((-0.75*(byr[10]+byl[10]))+0.4330127018922193*byr[5]-0.4330127018922193*byl[5])*c+0.75*psr[10]-0.75*psl[10]-0.4330127018922193*(psr[5]+psl[5]))*gamma; 
  incr[11] = (((-0.75*(byr[11]+byl[11]))+0.4330127018922194*byr[7]-0.4330127018922194*byl[7])*c+0.75*psr[11]-0.75*psl[11]-0.4330127018922194*(psr[7]+psl[7]))*gamma; 
  incr[12] = (((-1.25*byr[12])+1.25*byl[12]+0.9682458365518543*(byr[4]+byl[4])-0.5590169943749476*byr[1]+0.5590169943749476*byl[1])*c+1.25*(psr[12]+psl[12])-0.9682458365518543*psr[4]+0.9682458365518543*psl[4]+0.5590169943749476*(psr[1]+psl[1]))*gamma; 
  incr[13] = (0.25*(byl[13]*c+psr[13]+psl[13])-0.25*byr[13]*c)*gamma; 
  incr[14] = (((-1.25*byr[14])+1.25*byl[14]+0.9682458365518543*(byr[6]+byl[6])-0.5590169943749476*byr[3]+0.5590169943749476*byl[3])*c+1.25*(psr[14]+psl[14])-0.9682458365518543*psr[6]+0.9682458365518543*psl[6]+0.5590169943749476*(psr[3]+psl[3]))*gamma; 
  incr[15] = (0.25*(byl[15]*c+psr[15]+psl[15])-0.25*byr[15]*c)*gamma; 
  incr[16] = (((-0.75*(byr[16]+byl[16]))+0.4330127018922194*byr[9]-0.4330127018922194*byl[9])*c+0.75*psr[16]-0.75*psl[16]-0.4330127018922194*(psr[9]+psl[9]))*gamma; 
  incr[17] = (0.25*(byl[17]*c+psr[17]+psl[17])-0.25*byr[17]*c)*gamma; 
  incr[18] = (((-1.75*(byr[18]+byl[18]))+1.479019945774904*byr[8]-1.479019945774904*byl[8]-1.14564392373896*(byr[2]+byl[2])+0.6614378277661477*byr[0]-0.6614378277661477*byl[0])*c+1.75*psr[18]-1.75*psl[18]-1.479019945774904*(psr[8]+psl[8])+1.14564392373896*psr[2]-1.14564392373896*psl[2]-0.6614378277661477*(psr[0]+psl[0]))*gamma; 
  incr[19] = (0.25*(byl[19]*c+psr[19]+psl[19])-0.25*byr[19]*c)*gamma; 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 
  outByr[4] += incr[4]*dxr1; 
  outByr[5] += incr[5]*dxr1; 
  outByr[6] += incr[6]*dxr1; 
  outByr[7] += incr[7]*dxr1; 
  outByr[8] += incr[8]*dxr1; 
  outByr[9] += incr[9]*dxr1; 
  outByr[10] += incr[10]*dxr1; 
  outByr[11] += incr[11]*dxr1; 
  outByr[12] += incr[12]*dxr1; 
  outByr[13] += incr[13]*dxr1; 
  outByr[14] += incr[14]*dxr1; 
  outByr[15] += incr[15]*dxr1; 
  outByr[16] += incr[16]*dxr1; 
  outByr[17] += incr[17]*dxr1; 
  outByr[18] += incr[18]*dxr1; 
  outByr[19] += incr[19]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += -1.0*incr[1]*dxl1; 
  outByl[2] += incr[2]*dxl1; 
  outByl[3] += -1.0*incr[3]*dxl1; 
  outByl[4] += incr[4]*dxl1; 
  outByl[5] += -1.0*incr[5]*dxl1; 
  outByl[6] += incr[6]*dxl1; 
  outByl[7] += -1.0*incr[7]*dxl1; 
  outByl[8] += -1.0*incr[8]*dxl1; 
  outByl[9] += -1.0*incr[9]*dxl1; 
  outByl[10] += incr[10]*dxl1; 
  outByl[11] += incr[11]*dxl1; 
  outByl[12] += -1.0*incr[12]*dxl1; 
  outByl[13] += -1.0*incr[13]*dxl1; 
  outByl[14] += -1.0*incr[14]*dxl1; 
  outByl[15] += -1.0*incr[15]*dxl1; 
  outByl[16] += incr[16]*dxl1; 
  outByl[17] += -1.0*incr[17]*dxl1; 
  outByl[18] += incr[18]*dxl1; 
  outByl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(bzr[18]+bzl[18])-0.5590169943749475*bzr[8]+0.5590169943749475*bzl[8]+0.4330127018922193*(bzr[2]+bzl[2])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau+0.6614378277661477*exr[18]-0.6614378277661477*exl[18]-0.5590169943749475*(exr[8]+exl[8])+0.4330127018922193*exr[2]-0.4330127018922193*exl[2]-0.25*(exr[0]+exl[0]); 
  incr[1] = (((-0.5590169943749476*bzr[12])+0.5590169943749476*bzl[12]+0.4330127018922193*(bzr[4]+bzl[4])-0.25*bzr[1]+0.25*bzl[1])*c2)/tau-0.5590169943749476*(exr[12]+exl[12])+0.4330127018922193*exr[4]-0.4330127018922193*exl[4]-0.25*(exr[1]+exl[1]); 
  incr[2] = (((-1.14564392373896*(bzr[18]+bzl[18]))+0.9682458365518543*bzr[8]-0.9682458365518543*bzl[8]-0.75*(bzr[2]+bzl[2])+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau-1.14564392373896*exr[18]+1.14564392373896*exl[18]+0.9682458365518543*(exr[8]+exl[8])-0.75*exr[2]+0.75*exl[2]+0.4330127018922193*(exr[0]+exl[0]); 
  incr[3] = (((-0.5590169943749476*bzr[14])+0.5590169943749476*bzl[14]+0.4330127018922193*(bzr[6]+bzl[6])-0.25*bzr[3]+0.25*bzl[3])*c2)/tau-0.5590169943749476*(exr[14]+exl[14])+0.4330127018922193*exr[6]-0.4330127018922193*exl[6]-0.25*(exr[3]+exl[3]); 
  incr[4] = ((0.9682458365518543*bzr[12]-0.9682458365518543*bzl[12]-0.75*(bzr[4]+bzl[4])+0.4330127018922193*bzr[1]-0.4330127018922193*bzl[1])*c2)/tau+0.9682458365518543*(exr[12]+exl[12])-0.75*exr[4]+0.75*exl[4]+0.4330127018922193*(exr[1]+exl[1]); 
  incr[5] = ((0.4330127018922193*(bzr[10]+bzl[10])-0.25*bzr[5]+0.25*bzl[5])*c2)/tau+0.4330127018922193*exr[10]-0.4330127018922193*exl[10]-0.25*(exr[5]+exl[5]); 
  incr[6] = ((0.9682458365518543*bzr[14]-0.9682458365518543*bzl[14]-0.75*(bzr[6]+bzl[6])+0.4330127018922193*bzr[3]-0.4330127018922193*bzl[3])*c2)/tau+0.9682458365518543*(exr[14]+exl[14])-0.75*exr[6]+0.75*exl[6]+0.4330127018922193*(exr[3]+exl[3]); 
  incr[7] = ((0.4330127018922194*(bzr[11]+bzl[11])-0.25*bzr[7]+0.25*bzl[7])*c2)/tau+0.4330127018922194*exr[11]-0.4330127018922194*exl[11]-0.25*(exr[7]+exl[7]); 
  incr[8] = ((1.479019945774904*(bzr[18]+bzl[18])-1.25*bzr[8]+1.25*bzl[8]+0.9682458365518543*(bzr[2]+bzl[2])-0.5590169943749475*bzr[0]+0.5590169943749475*bzl[0])*c2)/tau+1.479019945774904*exr[18]-1.479019945774904*exl[18]-1.25*(exr[8]+exl[8])+0.9682458365518543*exr[2]-0.9682458365518543*exl[2]-0.5590169943749475*(exr[0]+exl[0]); 
  incr[9] = ((0.4330127018922194*(bzr[16]+bzl[16])-0.25*bzr[9]+0.25*bzl[9])*c2)/tau+0.4330127018922194*exr[16]-0.4330127018922194*exl[16]-0.25*(exr[9]+exl[9]); 
  incr[10] = (((-0.75*(bzr[10]+bzl[10]))+0.4330127018922193*bzr[5]-0.4330127018922193*bzl[5])*c2)/tau-0.75*exr[10]+0.75*exl[10]+0.4330127018922193*(exr[5]+exl[5]); 
  incr[11] = (((-0.75*(bzr[11]+bzl[11]))+0.4330127018922194*bzr[7]-0.4330127018922194*bzl[7])*c2)/tau-0.75*exr[11]+0.75*exl[11]+0.4330127018922194*(exr[7]+exl[7]); 
  incr[12] = (((-1.25*bzr[12])+1.25*bzl[12]+0.9682458365518543*(bzr[4]+bzl[4])-0.5590169943749476*bzr[1]+0.5590169943749476*bzl[1])*c2)/tau-1.25*(exr[12]+exl[12])+0.9682458365518543*exr[4]-0.9682458365518543*exl[4]-0.5590169943749476*(exr[1]+exl[1]); 
  incr[13] = ((0.25*bzl[13]-0.25*bzr[13])*c2)/tau-0.25*(exr[13]+exl[13]); 
  incr[14] = (((-1.25*bzr[14])+1.25*bzl[14]+0.9682458365518543*(bzr[6]+bzl[6])-0.5590169943749476*bzr[3]+0.5590169943749476*bzl[3])*c2)/tau-1.25*(exr[14]+exl[14])+0.9682458365518543*exr[6]-0.9682458365518543*exl[6]-0.5590169943749476*(exr[3]+exl[3]); 
  incr[15] = ((0.25*bzl[15]-0.25*bzr[15])*c2)/tau-0.25*(exr[15]+exl[15]); 
  incr[16] = (((-0.75*(bzr[16]+bzl[16]))+0.4330127018922194*bzr[9]-0.4330127018922194*bzl[9])*c2)/tau-0.75*exr[16]+0.75*exl[16]+0.4330127018922194*(exr[9]+exl[9]); 
  incr[17] = ((0.25*bzl[17]-0.25*bzr[17])*c2)/tau-0.25*(exr[17]+exl[17]); 
  incr[18] = (((-1.75*(bzr[18]+bzl[18]))+1.479019945774904*bzr[8]-1.479019945774904*bzl[8]-1.14564392373896*(bzr[2]+bzl[2])+0.6614378277661477*bzr[0]-0.6614378277661477*bzl[0])*c2)/tau-1.75*exr[18]+1.75*exl[18]+1.479019945774904*(exr[8]+exl[8])-1.14564392373896*exr[2]+1.14564392373896*exl[2]+0.6614378277661477*(exr[0]+exl[0]); 
  incr[19] = ((0.25*bzl[19]-0.25*bzr[19])*c2)/tau-0.25*(exr[19]+exl[19]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 
  outBzr[4] += incr[4]*dxr1; 
  outBzr[5] += incr[5]*dxr1; 
  outBzr[6] += incr[6]*dxr1; 
  outBzr[7] += incr[7]*dxr1; 
  outBzr[8] += incr[8]*dxr1; 
  outBzr[9] += incr[9]*dxr1; 
  outBzr[10] += incr[10]*dxr1; 
  outBzr[11] += incr[11]*dxr1; 
  outBzr[12] += incr[12]*dxr1; 
  outBzr[13] += incr[13]*dxr1; 
  outBzr[14] += incr[14]*dxr1; 
  outBzr[15] += incr[15]*dxr1; 
  outBzr[16] += incr[16]*dxr1; 
  outBzr[17] += incr[17]*dxr1; 
  outBzr[18] += incr[18]*dxr1; 
  outBzr[19] += incr[19]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += -1.0*incr[1]*dxl1; 
  outBzl[2] += incr[2]*dxl1; 
  outBzl[3] += -1.0*incr[3]*dxl1; 
  outBzl[4] += incr[4]*dxl1; 
  outBzl[5] += -1.0*incr[5]*dxl1; 
  outBzl[6] += incr[6]*dxl1; 
  outBzl[7] += -1.0*incr[7]*dxl1; 
  outBzl[8] += -1.0*incr[8]*dxl1; 
  outBzl[9] += -1.0*incr[9]*dxl1; 
  outBzl[10] += incr[10]*dxl1; 
  outBzl[11] += incr[11]*dxl1; 
  outBzl[12] += -1.0*incr[12]*dxl1; 
  outBzl[13] += -1.0*incr[13]*dxl1; 
  outBzl[14] += -1.0*incr[14]*dxl1; 
  outBzl[15] += -1.0*incr[15]*dxl1; 
  outBzl[16] += incr[16]*dxl1; 
  outBzl[17] += -1.0*incr[17]*dxl1; 
  outBzl[18] += incr[18]*dxl1; 
  outBzl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(phr[18]+phl[18])-0.5590169943749475*phr[8]+0.5590169943749475*phl[8]+0.4330127018922193*(phr[2]+phl[2])-0.25*phr[0]+0.25*phl[0])*c-0.6614378277661477*eyr[18]+0.6614378277661477*eyl[18]+0.5590169943749475*(eyr[8]+eyl[8])-0.4330127018922193*eyr[2]+0.4330127018922193*eyl[2]+0.25*(eyr[0]+eyl[0]))*chi; 
  incr[1] = (((-0.5590169943749476*phr[12])+0.5590169943749476*phl[12]+0.4330127018922193*(phr[4]+phl[4])-0.25*phr[1]+0.25*phl[1])*c+0.5590169943749476*(eyr[12]+eyl[12])-0.4330127018922193*eyr[4]+0.4330127018922193*eyl[4]+0.25*(eyr[1]+eyl[1]))*chi; 
  incr[2] = (((-1.14564392373896*(phr[18]+phl[18]))+0.9682458365518543*phr[8]-0.9682458365518543*phl[8]-0.75*(phr[2]+phl[2])+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+1.14564392373896*eyr[18]-1.14564392373896*eyl[18]-0.9682458365518543*(eyr[8]+eyl[8])+0.75*eyr[2]-0.75*eyl[2]-0.4330127018922193*(eyr[0]+eyl[0]))*chi; 
  incr[3] = (((-0.5590169943749476*phr[14])+0.5590169943749476*phl[14]+0.4330127018922193*(phr[6]+phl[6])-0.25*phr[3]+0.25*phl[3])*c+0.5590169943749476*(eyr[14]+eyl[14])-0.4330127018922193*eyr[6]+0.4330127018922193*eyl[6]+0.25*(eyr[3]+eyl[3]))*chi; 
  incr[4] = ((0.9682458365518543*phr[12]-0.9682458365518543*phl[12]-0.75*(phr[4]+phl[4])+0.4330127018922193*phr[1]-0.4330127018922193*phl[1])*c-0.9682458365518543*(eyr[12]+eyl[12])+0.75*eyr[4]-0.75*eyl[4]-0.4330127018922193*(eyr[1]+eyl[1]))*chi; 
  incr[5] = ((0.4330127018922193*(phr[10]+phl[10])-0.25*phr[5]+0.25*phl[5])*c-0.4330127018922193*eyr[10]+0.4330127018922193*eyl[10]+0.25*(eyr[5]+eyl[5]))*chi; 
  incr[6] = ((0.9682458365518543*phr[14]-0.9682458365518543*phl[14]-0.75*(phr[6]+phl[6])+0.4330127018922193*phr[3]-0.4330127018922193*phl[3])*c-0.9682458365518543*(eyr[14]+eyl[14])+0.75*eyr[6]-0.75*eyl[6]-0.4330127018922193*(eyr[3]+eyl[3]))*chi; 
  incr[7] = ((0.4330127018922194*(phr[11]+phl[11])-0.25*phr[7]+0.25*phl[7])*c-0.4330127018922194*eyr[11]+0.4330127018922194*eyl[11]+0.25*(eyr[7]+eyl[7]))*chi; 
  incr[8] = ((1.479019945774904*(phr[18]+phl[18])-1.25*phr[8]+1.25*phl[8]+0.9682458365518543*(phr[2]+phl[2])-0.5590169943749475*phr[0]+0.5590169943749475*phl[0])*c-1.479019945774904*eyr[18]+1.479019945774904*eyl[18]+1.25*(eyr[8]+eyl[8])-0.9682458365518543*eyr[2]+0.9682458365518543*eyl[2]+0.5590169943749475*(eyr[0]+eyl[0]))*chi; 
  incr[9] = ((0.4330127018922194*(phr[16]+phl[16])-0.25*phr[9]+0.25*phl[9])*c-0.4330127018922194*eyr[16]+0.4330127018922194*eyl[16]+0.25*(eyr[9]+eyl[9]))*chi; 
  incr[10] = (((-0.75*(phr[10]+phl[10]))+0.4330127018922193*phr[5]-0.4330127018922193*phl[5])*c+0.75*eyr[10]-0.75*eyl[10]-0.4330127018922193*(eyr[5]+eyl[5]))*chi; 
  incr[11] = (((-0.75*(phr[11]+phl[11]))+0.4330127018922194*phr[7]-0.4330127018922194*phl[7])*c+0.75*eyr[11]-0.75*eyl[11]-0.4330127018922194*(eyr[7]+eyl[7]))*chi; 
  incr[12] = (((-1.25*phr[12])+1.25*phl[12]+0.9682458365518543*(phr[4]+phl[4])-0.5590169943749476*phr[1]+0.5590169943749476*phl[1])*c+1.25*(eyr[12]+eyl[12])-0.9682458365518543*eyr[4]+0.9682458365518543*eyl[4]+0.5590169943749476*(eyr[1]+eyl[1]))*chi; 
  incr[13] = (0.25*(phl[13]*c+eyr[13]+eyl[13])-0.25*phr[13]*c)*chi; 
  incr[14] = (((-1.25*phr[14])+1.25*phl[14]+0.9682458365518543*(phr[6]+phl[6])-0.5590169943749476*phr[3]+0.5590169943749476*phl[3])*c+1.25*(eyr[14]+eyl[14])-0.9682458365518543*eyr[6]+0.9682458365518543*eyl[6]+0.5590169943749476*(eyr[3]+eyl[3]))*chi; 
  incr[15] = (0.25*(phl[15]*c+eyr[15]+eyl[15])-0.25*phr[15]*c)*chi; 
  incr[16] = (((-0.75*(phr[16]+phl[16]))+0.4330127018922194*phr[9]-0.4330127018922194*phl[9])*c+0.75*eyr[16]-0.75*eyl[16]-0.4330127018922194*(eyr[9]+eyl[9]))*chi; 
  incr[17] = (0.25*(phl[17]*c+eyr[17]+eyl[17])-0.25*phr[17]*c)*chi; 
  incr[18] = (((-1.75*(phr[18]+phl[18]))+1.479019945774904*phr[8]-1.479019945774904*phl[8]-1.14564392373896*(phr[2]+phl[2])+0.6614378277661477*phr[0]-0.6614378277661477*phl[0])*c+1.75*eyr[18]-1.75*eyl[18]-1.479019945774904*(eyr[8]+eyl[8])+1.14564392373896*eyr[2]-1.14564392373896*eyl[2]-0.6614378277661477*(eyr[0]+eyl[0]))*chi; 
  incr[19] = (0.25*(phl[19]*c+eyr[19]+eyl[19])-0.25*phr[19]*c)*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 
  outPhr[4] += incr[4]*dxr1; 
  outPhr[5] += incr[5]*dxr1; 
  outPhr[6] += incr[6]*dxr1; 
  outPhr[7] += incr[7]*dxr1; 
  outPhr[8] += incr[8]*dxr1; 
  outPhr[9] += incr[9]*dxr1; 
  outPhr[10] += incr[10]*dxr1; 
  outPhr[11] += incr[11]*dxr1; 
  outPhr[12] += incr[12]*dxr1; 
  outPhr[13] += incr[13]*dxr1; 
  outPhr[14] += incr[14]*dxr1; 
  outPhr[15] += incr[15]*dxr1; 
  outPhr[16] += incr[16]*dxr1; 
  outPhr[17] += incr[17]*dxr1; 
  outPhr[18] += incr[18]*dxr1; 
  outPhr[19] += incr[19]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += -1.0*incr[1]*dxl1; 
  outPhl[2] += incr[2]*dxl1; 
  outPhl[3] += -1.0*incr[3]*dxl1; 
  outPhl[4] += incr[4]*dxl1; 
  outPhl[5] += -1.0*incr[5]*dxl1; 
  outPhl[6] += incr[6]*dxl1; 
  outPhl[7] += -1.0*incr[7]*dxl1; 
  outPhl[8] += -1.0*incr[8]*dxl1; 
  outPhl[9] += -1.0*incr[9]*dxl1; 
  outPhl[10] += incr[10]*dxl1; 
  outPhl[11] += incr[11]*dxl1; 
  outPhl[12] += -1.0*incr[12]*dxl1; 
  outPhl[13] += -1.0*incr[13]*dxl1; 
  outPhl[14] += -1.0*incr[14]*dxl1; 
  outPhl[15] += -1.0*incr[15]*dxl1; 
  outPhl[16] += incr[16]*dxl1; 
  outPhl[17] += -1.0*incr[17]*dxl1; 
  outPhl[18] += incr[18]*dxl1; 
  outPhl[19] += -1.0*incr[19]*dxl1; 

 
  incr[0] = (0.6614378277661477*(psr[18]+psl[18])-0.5590169943749475*psr[8]+0.5590169943749475*psl[8]+0.4330127018922193*(psr[2]+psl[2])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.6614378277661477*byr[18])+0.6614378277661477*byl[18]+0.5590169943749475*(byr[8]+byl[8])-0.4330127018922193*byr[2]+0.4330127018922193*byl[2]+0.25*(byr[0]+byl[0]))*c2gamma; 
  incr[1] = ((-0.5590169943749476*psr[12])+0.5590169943749476*psl[12]+0.4330127018922193*(psr[4]+psl[4])-0.25*psr[1]+0.25*psl[1])*c*gamma+(0.5590169943749476*(byr[12]+byl[12])-0.4330127018922193*byr[4]+0.4330127018922193*byl[4]+0.25*(byr[1]+byl[1]))*c2gamma; 
  incr[2] = ((-1.14564392373896*(psr[18]+psl[18]))+0.9682458365518543*psr[8]-0.9682458365518543*psl[8]-0.75*(psr[2]+psl[2])+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(1.14564392373896*byr[18]-1.14564392373896*byl[18]-0.9682458365518543*(byr[8]+byl[8])+0.75*byr[2]-0.75*byl[2]-0.4330127018922193*(byr[0]+byl[0]))*c2gamma; 
  incr[3] = ((-0.5590169943749476*psr[14])+0.5590169943749476*psl[14]+0.4330127018922193*(psr[6]+psl[6])-0.25*psr[3]+0.25*psl[3])*c*gamma+(0.5590169943749476*(byr[14]+byl[14])-0.4330127018922193*byr[6]+0.4330127018922193*byl[6]+0.25*(byr[3]+byl[3]))*c2gamma; 
  incr[4] = (0.9682458365518543*psr[12]-0.9682458365518543*psl[12]-0.75*(psr[4]+psl[4])+0.4330127018922193*psr[1]-0.4330127018922193*psl[1])*c*gamma+((-0.9682458365518543*(byr[12]+byl[12]))+0.75*byr[4]-0.75*byl[4]-0.4330127018922193*(byr[1]+byl[1]))*c2gamma; 
  incr[5] = (0.4330127018922193*(psr[10]+psl[10])-0.25*psr[5]+0.25*psl[5])*c*gamma+((-0.4330127018922193*byr[10])+0.4330127018922193*byl[10]+0.25*(byr[5]+byl[5]))*c2gamma; 
  incr[6] = (0.9682458365518543*psr[14]-0.9682458365518543*psl[14]-0.75*(psr[6]+psl[6])+0.4330127018922193*psr[3]-0.4330127018922193*psl[3])*c*gamma+((-0.9682458365518543*(byr[14]+byl[14]))+0.75*byr[6]-0.75*byl[6]-0.4330127018922193*(byr[3]+byl[3]))*c2gamma; 
  incr[7] = (0.4330127018922194*(psr[11]+psl[11])-0.25*psr[7]+0.25*psl[7])*c*gamma+((-0.4330127018922194*byr[11])+0.4330127018922194*byl[11]+0.25*(byr[7]+byl[7]))*c2gamma; 
  incr[8] = (1.479019945774904*(psr[18]+psl[18])-1.25*psr[8]+1.25*psl[8]+0.9682458365518543*(psr[2]+psl[2])-0.5590169943749475*psr[0]+0.5590169943749475*psl[0])*c*gamma+((-1.479019945774904*byr[18])+1.479019945774904*byl[18]+1.25*(byr[8]+byl[8])-0.9682458365518543*byr[2]+0.9682458365518543*byl[2]+0.5590169943749475*(byr[0]+byl[0]))*c2gamma; 
  incr[9] = (0.4330127018922194*(psr[16]+psl[16])-0.25*psr[9]+0.25*psl[9])*c*gamma+((-0.4330127018922194*byr[16])+0.4330127018922194*byl[16]+0.25*(byr[9]+byl[9]))*c2gamma; 
  incr[10] = ((-0.75*(psr[10]+psl[10]))+0.4330127018922193*psr[5]-0.4330127018922193*psl[5])*c*gamma+(0.75*byr[10]-0.75*byl[10]-0.4330127018922193*(byr[5]+byl[5]))*c2gamma; 
  incr[11] = ((-0.75*(psr[11]+psl[11]))+0.4330127018922194*psr[7]-0.4330127018922194*psl[7])*c*gamma+(0.75*byr[11]-0.75*byl[11]-0.4330127018922194*(byr[7]+byl[7]))*c2gamma; 
  incr[12] = ((-1.25*psr[12])+1.25*psl[12]+0.9682458365518543*(psr[4]+psl[4])-0.5590169943749476*psr[1]+0.5590169943749476*psl[1])*c*gamma+(1.25*(byr[12]+byl[12])-0.9682458365518543*byr[4]+0.9682458365518543*byl[4]+0.5590169943749476*(byr[1]+byl[1]))*c2gamma; 
  incr[13] = (0.25*psl[13]-0.25*psr[13])*c*gamma+0.25*(byr[13]+byl[13])*c2gamma; 
  incr[14] = ((-1.25*psr[14])+1.25*psl[14]+0.9682458365518543*(psr[6]+psl[6])-0.5590169943749476*psr[3]+0.5590169943749476*psl[3])*c*gamma+(1.25*(byr[14]+byl[14])-0.9682458365518543*byr[6]+0.9682458365518543*byl[6]+0.5590169943749476*(byr[3]+byl[3]))*c2gamma; 
  incr[15] = (0.25*psl[15]-0.25*psr[15])*c*gamma+0.25*(byr[15]+byl[15])*c2gamma; 
  incr[16] = ((-0.75*(psr[16]+psl[16]))+0.4330127018922194*psr[9]-0.4330127018922194*psl[9])*c*gamma+(0.75*byr[16]-0.75*byl[16]-0.4330127018922194*(byr[9]+byl[9]))*c2gamma; 
  incr[17] = (0.25*psl[17]-0.25*psr[17])*c*gamma+0.25*(byr[17]+byl[17])*c2gamma; 
  incr[18] = ((-1.75*(psr[18]+psl[18]))+1.479019945774904*psr[8]-1.479019945774904*psl[8]-1.14564392373896*(psr[2]+psl[2])+0.6614378277661477*psr[0]-0.6614378277661477*psl[0])*c*gamma+(1.75*byr[18]-1.75*byl[18]-1.479019945774904*(byr[8]+byl[8])+1.14564392373896*byr[2]-1.14564392373896*byl[2]-0.6614378277661477*(byr[0]+byl[0]))*c2gamma; 
  incr[19] = (0.25*psl[19]-0.25*psr[19])*c*gamma+0.25*(byr[19]+byl[19])*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 
  outPsr[4] += incr[4]*dxr1; 
  outPsr[5] += incr[5]*dxr1; 
  outPsr[6] += incr[6]*dxr1; 
  outPsr[7] += incr[7]*dxr1; 
  outPsr[8] += incr[8]*dxr1; 
  outPsr[9] += incr[9]*dxr1; 
  outPsr[10] += incr[10]*dxr1; 
  outPsr[11] += incr[11]*dxr1; 
  outPsr[12] += incr[12]*dxr1; 
  outPsr[13] += incr[13]*dxr1; 
  outPsr[14] += incr[14]*dxr1; 
  outPsr[15] += incr[15]*dxr1; 
  outPsr[16] += incr[16]*dxr1; 
  outPsr[17] += incr[17]*dxr1; 
  outPsr[18] += incr[18]*dxr1; 
  outPsr[19] += incr[19]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += -1.0*incr[1]*dxl1; 
  outPsl[2] += incr[2]*dxl1; 
  outPsl[3] += -1.0*incr[3]*dxl1; 
  outPsl[4] += incr[4]*dxl1; 
  outPsl[5] += -1.0*incr[5]*dxl1; 
  outPsl[6] += incr[6]*dxl1; 
  outPsl[7] += -1.0*incr[7]*dxl1; 
  outPsl[8] += -1.0*incr[8]*dxl1; 
  outPsl[9] += -1.0*incr[9]*dxl1; 
  outPsl[10] += incr[10]*dxl1; 
  outPsl[11] += incr[11]*dxl1; 
  outPsl[12] += -1.0*incr[12]*dxl1; 
  outPsl[13] += -1.0*incr[13]*dxl1; 
  outPsl[14] += -1.0*incr[14]*dxl1; 
  outPsl[15] += -1.0*incr[15]*dxl1; 
  outPsl[16] += incr[16]*dxl1; 
  outPsl[17] += -1.0*incr[17]*dxl1; 
  outPsl[18] += incr[18]*dxl1; 
  outPsl[19] += -1.0*incr[19]*dxl1; 

 
  return c; 
} 
double MaxwellSurf3xMax_Z_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[2]; 
  const double dxr1 = 2.0/dxr[2]; 
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
 
  incr[0] = (0.4330127018922193*(exr[3]+exl[3])-0.25*exr[0]+0.25*exl[0])*tau+((-0.4330127018922193*byr[3])+0.4330127018922193*byl[3]+0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = 0.25*(exl[1]*tau+(byr[1]+byl[1])*c2)-0.25*exr[1]*tau; 
  incr[2] = 0.25*(exl[2]*tau+(byr[2]+byl[2])*c2)-0.25*exr[2]*tau; 
  incr[3] = ((-0.75*(exr[3]+exl[3]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*tau+(0.75*byr[3]-0.75*byl[3]-0.4330127018922193*(byr[0]+byl[0]))*c2; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += -1.0*incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 
  outExl[3] += incr[3]*dxl1; 

 
  incr[0] = (0.4330127018922193*(eyr[3]+eyl[3])-0.25*eyr[0]+0.25*eyl[0])*tau+(0.4330127018922193*bxr[3]-0.4330127018922193*bxl[3]-0.25*(bxr[0]+bxl[0]))*c2; 
  incr[1] = (0.25*eyl[1]-0.25*eyr[1])*tau-0.25*(bxr[1]+bxl[1])*c2; 
  incr[2] = (0.25*eyl[2]-0.25*eyr[2])*tau-0.25*(bxr[2]+bxl[2])*c2; 
  incr[3] = ((-0.75*(eyr[3]+eyl[3]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+((-0.75*bxr[3])+0.75*bxl[3]+0.4330127018922193*(bxr[0]+bxl[0]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += -1.0*incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 
  outEyl[3] += incr[3]*dxl1; 

 
  incr[0] = (0.4330127018922193*(ezr[3]+ezl[3])-0.25*ezr[0]+0.25*ezl[0])*c*chi+((-0.4330127018922193*phr[3])+0.4330127018922193*phl[3]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.25*ezl[1]-0.25*ezr[1])*c*chi+0.25*(phr[1]+phl[1])*c2chi; 
  incr[2] = (0.25*ezl[2]-0.25*ezr[2])*c*chi+0.25*(phr[2]+phl[2])*c2chi; 
  incr[3] = ((-0.75*(ezr[3]+ezl[3]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*c*chi+(0.75*phr[3]-0.75*phl[3]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += -1.0*incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 
  outEzl[3] += incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bxr[3]+bxl[3])-0.25*bxr[0]+0.25*bxl[0])*c2)/tau+0.4330127018922193*eyr[3]-0.4330127018922193*eyl[3]-0.25*(eyr[0]+eyl[0]); 
  incr[1] = ((0.25*bxl[1]-0.25*bxr[1])*c2)/tau-0.25*(eyr[1]+eyl[1]); 
  incr[2] = ((0.25*bxl[2]-0.25*bxr[2])*c2)/tau-0.25*(eyr[2]+eyl[2]); 
  incr[3] = (((-0.75*(bxr[3]+bxl[3]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c2)/tau-0.75*eyr[3]+0.75*eyl[3]+0.4330127018922193*(eyr[0]+eyl[0]); 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += -1.0*incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 
  outBxl[3] += incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(byr[3]+byl[3])-0.25*byr[0]+0.25*byl[0])*c2)/tau-0.4330127018922193*exr[3]+0.4330127018922193*exl[3]+0.25*(exr[0]+exl[0]); 
  incr[1] = ((0.25*byl[1]-0.25*byr[1])*c2)/tau+0.25*(exr[1]+exl[1]); 
  incr[2] = ((0.25*byl[2]-0.25*byr[2])*c2)/tau+0.25*(exr[2]+exl[2]); 
  incr[3] = (((-0.75*(byr[3]+byl[3]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau+0.75*exr[3]-0.75*exl[3]-0.4330127018922193*(exr[0]+exl[0]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += -1.0*incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 
  outByl[3] += incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bzr[3]+bzl[3])-0.25*bzr[0]+0.25*bzl[0])*c-0.4330127018922193*psr[3]+0.4330127018922193*psl[3]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (0.25*(bzl[1]*c+psr[1]+psl[1])-0.25*bzr[1]*c)*gamma; 
  incr[2] = (0.25*(bzl[2]*c+psr[2]+psl[2])-0.25*bzr[2]*c)*gamma; 
  incr[3] = (((-0.75*(bzr[3]+bzl[3]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c+0.75*psr[3]-0.75*psl[3]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += -1.0*incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 
  outBzl[3] += incr[3]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(phr[3]+phl[3])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*ezr[3]+0.4330127018922193*ezl[3]+0.25*(ezr[0]+ezl[0]))*chi; 
  incr[1] = (0.25*(phl[1]*c+ezr[1]+ezl[1])-0.25*phr[1]*c)*chi; 
  incr[2] = (0.25*(phl[2]*c+ezr[2]+ezl[2])-0.25*phr[2]*c)*chi; 
  incr[3] = (((-0.75*(phr[3]+phl[3]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*ezr[3]-0.75*ezl[3]-0.4330127018922193*(ezr[0]+ezl[0]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += -1.0*incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 
  outPhl[3] += incr[3]*dxl1; 

 
  incr[0] = (0.4330127018922193*(psr[3]+psl[3])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*bzr[3])+0.4330127018922193*bzl[3]+0.25*(bzr[0]+bzl[0]))*c2gamma; 
  incr[1] = (0.25*psl[1]-0.25*psr[1])*c*gamma+0.25*(bzr[1]+bzl[1])*c2gamma; 
  incr[2] = (0.25*psl[2]-0.25*psr[2])*c*gamma+0.25*(bzr[2]+bzl[2])*c2gamma; 
  incr[3] = ((-0.75*(psr[3]+psl[3]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*bzr[3]-0.75*bzl[3]-0.4330127018922193*(bzr[0]+bzl[0]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += -1.0*incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 
  outPsl[3] += incr[3]*dxl1; 

 
  return c; 
} 
double MaxwellSurf3xMax_Z_P2(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[2]; 
  const double dxr1 = 2.0/dxr[2]; 
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
 
  incr[0] = ((-0.5590169943749475*exr[9])+0.5590169943749475*exl[9]+0.4330127018922193*(exr[3]+exl[3])-0.25*exr[0]+0.25*exl[0])*tau+(0.5590169943749475*(byr[9]+byl[9])-0.4330127018922193*byr[3]+0.4330127018922193*byl[3]+0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = (0.4330127018922193*(exr[5]+exl[5])-0.25*exr[1]+0.25*exl[1])*tau+((-0.4330127018922193*byr[5])+0.4330127018922193*byl[5]+0.25*(byr[1]+byl[1]))*c2; 
  incr[2] = (0.4330127018922193*(exr[6]+exl[6])-0.25*exr[2]+0.25*exl[2])*tau+((-0.4330127018922193*byr[6])+0.4330127018922193*byl[6]+0.25*(byr[2]+byl[2]))*c2; 
  incr[3] = (0.9682458365518543*exr[9]-0.9682458365518543*exl[9]-0.75*(exr[3]+exl[3])+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*tau+((-0.9682458365518543*(byr[9]+byl[9]))+0.75*byr[3]-0.75*byl[3]-0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[4] = 0.25*(exl[4]*tau+(byr[4]+byl[4])*c2)-0.25*exr[4]*tau; 
  incr[5] = ((-0.75*(exr[5]+exl[5]))+0.4330127018922193*exr[1]-0.4330127018922193*exl[1])*tau+(0.75*byr[5]-0.75*byl[5]-0.4330127018922193*(byr[1]+byl[1]))*c2; 
  incr[6] = ((-0.75*(exr[6]+exl[6]))+0.4330127018922193*exr[2]-0.4330127018922193*exl[2])*tau+(0.75*byr[6]-0.75*byl[6]-0.4330127018922193*(byr[2]+byl[2]))*c2; 
  incr[7] = 0.25*(exl[7]*tau+(byr[7]+byl[7])*c2)-0.25*exr[7]*tau; 
  incr[8] = 0.25*(exl[8]*tau+(byr[8]+byl[8])*c2)-0.25*exr[8]*tau; 
  incr[9] = ((-1.25*exr[9])+1.25*exl[9]+0.9682458365518543*(exr[3]+exl[3])-0.5590169943749475*exr[0]+0.5590169943749475*exl[0])*tau+(1.25*(byr[9]+byl[9])-0.9682458365518543*byr[3]+0.9682458365518543*byl[3]+0.5590169943749475*(byr[0]+byl[0]))*c2; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 
  outExr[4] += incr[4]*dxr1; 
  outExr[5] += incr[5]*dxr1; 
  outExr[6] += incr[6]*dxr1; 
  outExr[7] += incr[7]*dxr1; 
  outExr[8] += incr[8]*dxr1; 
  outExr[9] += incr[9]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += -1.0*incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 
  outExl[3] += incr[3]*dxl1; 
  outExl[4] += -1.0*incr[4]*dxl1; 
  outExl[5] += incr[5]*dxl1; 
  outExl[6] += incr[6]*dxl1; 
  outExl[7] += -1.0*incr[7]*dxl1; 
  outExl[8] += -1.0*incr[8]*dxl1; 
  outExl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*eyr[9])+0.5590169943749475*eyl[9]+0.4330127018922193*(eyr[3]+eyl[3])-0.25*eyr[0]+0.25*eyl[0])*tau+((-0.5590169943749475*(bxr[9]+bxl[9]))+0.4330127018922193*bxr[3]-0.4330127018922193*bxl[3]-0.25*(bxr[0]+bxl[0]))*c2; 
  incr[1] = (0.4330127018922193*(eyr[5]+eyl[5])-0.25*eyr[1]+0.25*eyl[1])*tau+(0.4330127018922193*bxr[5]-0.4330127018922193*bxl[5]-0.25*(bxr[1]+bxl[1]))*c2; 
  incr[2] = (0.4330127018922193*(eyr[6]+eyl[6])-0.25*eyr[2]+0.25*eyl[2])*tau+(0.4330127018922193*bxr[6]-0.4330127018922193*bxl[6]-0.25*(bxr[2]+bxl[2]))*c2; 
  incr[3] = (0.9682458365518543*eyr[9]-0.9682458365518543*eyl[9]-0.75*(eyr[3]+eyl[3])+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+(0.9682458365518543*(bxr[9]+bxl[9])-0.75*bxr[3]+0.75*bxl[3]+0.4330127018922193*(bxr[0]+bxl[0]))*c2; 
  incr[4] = (0.25*eyl[4]-0.25*eyr[4])*tau-0.25*(bxr[4]+bxl[4])*c2; 
  incr[5] = ((-0.75*(eyr[5]+eyl[5]))+0.4330127018922193*eyr[1]-0.4330127018922193*eyl[1])*tau+((-0.75*bxr[5])+0.75*bxl[5]+0.4330127018922193*(bxr[1]+bxl[1]))*c2; 
  incr[6] = ((-0.75*(eyr[6]+eyl[6]))+0.4330127018922193*eyr[2]-0.4330127018922193*eyl[2])*tau+((-0.75*bxr[6])+0.75*bxl[6]+0.4330127018922193*(bxr[2]+bxl[2]))*c2; 
  incr[7] = (0.25*eyl[7]-0.25*eyr[7])*tau-0.25*(bxr[7]+bxl[7])*c2; 
  incr[8] = (0.25*eyl[8]-0.25*eyr[8])*tau-0.25*(bxr[8]+bxl[8])*c2; 
  incr[9] = ((-1.25*eyr[9])+1.25*eyl[9]+0.9682458365518543*(eyr[3]+eyl[3])-0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0])*tau+((-1.25*(bxr[9]+bxl[9]))+0.9682458365518543*bxr[3]-0.9682458365518543*bxl[3]-0.5590169943749475*(bxr[0]+bxl[0]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 
  outEyr[4] += incr[4]*dxr1; 
  outEyr[5] += incr[5]*dxr1; 
  outEyr[6] += incr[6]*dxr1; 
  outEyr[7] += incr[7]*dxr1; 
  outEyr[8] += incr[8]*dxr1; 
  outEyr[9] += incr[9]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += -1.0*incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 
  outEyl[3] += incr[3]*dxl1; 
  outEyl[4] += -1.0*incr[4]*dxl1; 
  outEyl[5] += incr[5]*dxl1; 
  outEyl[6] += incr[6]*dxl1; 
  outEyl[7] += -1.0*incr[7]*dxl1; 
  outEyl[8] += -1.0*incr[8]*dxl1; 
  outEyl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*ezr[9])+0.5590169943749475*ezl[9]+0.4330127018922193*(ezr[3]+ezl[3])-0.25*ezr[0]+0.25*ezl[0])*c*chi+(0.5590169943749475*(phr[9]+phl[9])-0.4330127018922193*phr[3]+0.4330127018922193*phl[3]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.4330127018922193*(ezr[5]+ezl[5])-0.25*ezr[1]+0.25*ezl[1])*c*chi+((-0.4330127018922193*phr[5])+0.4330127018922193*phl[5]+0.25*(phr[1]+phl[1]))*c2chi; 
  incr[2] = (0.4330127018922193*(ezr[6]+ezl[6])-0.25*ezr[2]+0.25*ezl[2])*c*chi+((-0.4330127018922193*phr[6])+0.4330127018922193*phl[6]+0.25*(phr[2]+phl[2]))*c2chi; 
  incr[3] = (0.9682458365518543*ezr[9]-0.9682458365518543*ezl[9]-0.75*(ezr[3]+ezl[3])+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*c*chi+((-0.9682458365518543*(phr[9]+phl[9]))+0.75*phr[3]-0.75*phl[3]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[4] = (0.25*ezl[4]-0.25*ezr[4])*c*chi+0.25*(phr[4]+phl[4])*c2chi; 
  incr[5] = ((-0.75*(ezr[5]+ezl[5]))+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1])*c*chi+(0.75*phr[5]-0.75*phl[5]-0.4330127018922193*(phr[1]+phl[1]))*c2chi; 
  incr[6] = ((-0.75*(ezr[6]+ezl[6]))+0.4330127018922193*ezr[2]-0.4330127018922193*ezl[2])*c*chi+(0.75*phr[6]-0.75*phl[6]-0.4330127018922193*(phr[2]+phl[2]))*c2chi; 
  incr[7] = (0.25*ezl[7]-0.25*ezr[7])*c*chi+0.25*(phr[7]+phl[7])*c2chi; 
  incr[8] = (0.25*ezl[8]-0.25*ezr[8])*c*chi+0.25*(phr[8]+phl[8])*c2chi; 
  incr[9] = ((-1.25*ezr[9])+1.25*ezl[9]+0.9682458365518543*(ezr[3]+ezl[3])-0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0])*c*chi+(1.25*(phr[9]+phl[9])-0.9682458365518543*phr[3]+0.9682458365518543*phl[3]+0.5590169943749475*(phr[0]+phl[0]))*c2chi; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 
  outEzr[4] += incr[4]*dxr1; 
  outEzr[5] += incr[5]*dxr1; 
  outEzr[6] += incr[6]*dxr1; 
  outEzr[7] += incr[7]*dxr1; 
  outEzr[8] += incr[8]*dxr1; 
  outEzr[9] += incr[9]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += -1.0*incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 
  outEzl[3] += incr[3]*dxl1; 
  outEzl[4] += -1.0*incr[4]*dxl1; 
  outEzl[5] += incr[5]*dxl1; 
  outEzl[6] += incr[6]*dxl1; 
  outEzl[7] += -1.0*incr[7]*dxl1; 
  outEzl[8] += -1.0*incr[8]*dxl1; 
  outEzl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*bxr[9])+0.5590169943749475*bxl[9]+0.4330127018922193*(bxr[3]+bxl[3])-0.25*bxr[0]+0.25*bxl[0])*c2)/tau-0.5590169943749475*(eyr[9]+eyl[9])+0.4330127018922193*eyr[3]-0.4330127018922193*eyl[3]-0.25*(eyr[0]+eyl[0]); 
  incr[1] = ((0.4330127018922193*(bxr[5]+bxl[5])-0.25*bxr[1]+0.25*bxl[1])*c2)/tau+0.4330127018922193*eyr[5]-0.4330127018922193*eyl[5]-0.25*(eyr[1]+eyl[1]); 
  incr[2] = ((0.4330127018922193*(bxr[6]+bxl[6])-0.25*bxr[2]+0.25*bxl[2])*c2)/tau+0.4330127018922193*eyr[6]-0.4330127018922193*eyl[6]-0.25*(eyr[2]+eyl[2]); 
  incr[3] = ((0.9682458365518543*bxr[9]-0.9682458365518543*bxl[9]-0.75*(bxr[3]+bxl[3])+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c2)/tau+0.9682458365518543*(eyr[9]+eyl[9])-0.75*eyr[3]+0.75*eyl[3]+0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[4] = ((0.25*bxl[4]-0.25*bxr[4])*c2)/tau-0.25*(eyr[4]+eyl[4]); 
  incr[5] = (((-0.75*(bxr[5]+bxl[5]))+0.4330127018922193*bxr[1]-0.4330127018922193*bxl[1])*c2)/tau-0.75*eyr[5]+0.75*eyl[5]+0.4330127018922193*(eyr[1]+eyl[1]); 
  incr[6] = (((-0.75*(bxr[6]+bxl[6]))+0.4330127018922193*bxr[2]-0.4330127018922193*bxl[2])*c2)/tau-0.75*eyr[6]+0.75*eyl[6]+0.4330127018922193*(eyr[2]+eyl[2]); 
  incr[7] = ((0.25*bxl[7]-0.25*bxr[7])*c2)/tau-0.25*(eyr[7]+eyl[7]); 
  incr[8] = ((0.25*bxl[8]-0.25*bxr[8])*c2)/tau-0.25*(eyr[8]+eyl[8]); 
  incr[9] = (((-1.25*bxr[9])+1.25*bxl[9]+0.9682458365518543*(bxr[3]+bxl[3])-0.5590169943749475*bxr[0]+0.5590169943749475*bxl[0])*c2)/tau-1.25*(eyr[9]+eyl[9])+0.9682458365518543*eyr[3]-0.9682458365518543*eyl[3]-0.5590169943749475*(eyr[0]+eyl[0]); 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 
  outBxr[4] += incr[4]*dxr1; 
  outBxr[5] += incr[5]*dxr1; 
  outBxr[6] += incr[6]*dxr1; 
  outBxr[7] += incr[7]*dxr1; 
  outBxr[8] += incr[8]*dxr1; 
  outBxr[9] += incr[9]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += -1.0*incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 
  outBxl[3] += incr[3]*dxl1; 
  outBxl[4] += -1.0*incr[4]*dxl1; 
  outBxl[5] += incr[5]*dxl1; 
  outBxl[6] += incr[6]*dxl1; 
  outBxl[7] += -1.0*incr[7]*dxl1; 
  outBxl[8] += -1.0*incr[8]*dxl1; 
  outBxl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*byr[9])+0.5590169943749475*byl[9]+0.4330127018922193*(byr[3]+byl[3])-0.25*byr[0]+0.25*byl[0])*c2)/tau+0.5590169943749475*(exr[9]+exl[9])-0.4330127018922193*exr[3]+0.4330127018922193*exl[3]+0.25*(exr[0]+exl[0]); 
  incr[1] = ((0.4330127018922193*(byr[5]+byl[5])-0.25*byr[1]+0.25*byl[1])*c2)/tau-0.4330127018922193*exr[5]+0.4330127018922193*exl[5]+0.25*(exr[1]+exl[1]); 
  incr[2] = ((0.4330127018922193*(byr[6]+byl[6])-0.25*byr[2]+0.25*byl[2])*c2)/tau-0.4330127018922193*exr[6]+0.4330127018922193*exl[6]+0.25*(exr[2]+exl[2]); 
  incr[3] = ((0.9682458365518543*byr[9]-0.9682458365518543*byl[9]-0.75*(byr[3]+byl[3])+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau-0.9682458365518543*(exr[9]+exl[9])+0.75*exr[3]-0.75*exl[3]-0.4330127018922193*(exr[0]+exl[0]); 
  incr[4] = ((0.25*byl[4]-0.25*byr[4])*c2)/tau+0.25*(exr[4]+exl[4]); 
  incr[5] = (((-0.75*(byr[5]+byl[5]))+0.4330127018922193*byr[1]-0.4330127018922193*byl[1])*c2)/tau+0.75*exr[5]-0.75*exl[5]-0.4330127018922193*(exr[1]+exl[1]); 
  incr[6] = (((-0.75*(byr[6]+byl[6]))+0.4330127018922193*byr[2]-0.4330127018922193*byl[2])*c2)/tau+0.75*exr[6]-0.75*exl[6]-0.4330127018922193*(exr[2]+exl[2]); 
  incr[7] = ((0.25*byl[7]-0.25*byr[7])*c2)/tau+0.25*(exr[7]+exl[7]); 
  incr[8] = ((0.25*byl[8]-0.25*byr[8])*c2)/tau+0.25*(exr[8]+exl[8]); 
  incr[9] = (((-1.25*byr[9])+1.25*byl[9]+0.9682458365518543*(byr[3]+byl[3])-0.5590169943749475*byr[0]+0.5590169943749475*byl[0])*c2)/tau+1.25*(exr[9]+exl[9])-0.9682458365518543*exr[3]+0.9682458365518543*exl[3]+0.5590169943749475*(exr[0]+exl[0]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 
  outByr[4] += incr[4]*dxr1; 
  outByr[5] += incr[5]*dxr1; 
  outByr[6] += incr[6]*dxr1; 
  outByr[7] += incr[7]*dxr1; 
  outByr[8] += incr[8]*dxr1; 
  outByr[9] += incr[9]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += -1.0*incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 
  outByl[3] += incr[3]*dxl1; 
  outByl[4] += -1.0*incr[4]*dxl1; 
  outByl[5] += incr[5]*dxl1; 
  outByl[6] += incr[6]*dxl1; 
  outByl[7] += -1.0*incr[7]*dxl1; 
  outByl[8] += -1.0*incr[8]*dxl1; 
  outByl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*bzr[9])+0.5590169943749475*bzl[9]+0.4330127018922193*(bzr[3]+bzl[3])-0.25*bzr[0]+0.25*bzl[0])*c+0.5590169943749475*(psr[9]+psl[9])-0.4330127018922193*psr[3]+0.4330127018922193*psl[3]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = ((0.4330127018922193*(bzr[5]+bzl[5])-0.25*bzr[1]+0.25*bzl[1])*c-0.4330127018922193*psr[5]+0.4330127018922193*psl[5]+0.25*(psr[1]+psl[1]))*gamma; 
  incr[2] = ((0.4330127018922193*(bzr[6]+bzl[6])-0.25*bzr[2]+0.25*bzl[2])*c-0.4330127018922193*psr[6]+0.4330127018922193*psl[6]+0.25*(psr[2]+psl[2]))*gamma; 
  incr[3] = ((0.9682458365518543*bzr[9]-0.9682458365518543*bzl[9]-0.75*(bzr[3]+bzl[3])+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c-0.9682458365518543*(psr[9]+psl[9])+0.75*psr[3]-0.75*psl[3]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[4] = (0.25*(bzl[4]*c+psr[4]+psl[4])-0.25*bzr[4]*c)*gamma; 
  incr[5] = (((-0.75*(bzr[5]+bzl[5]))+0.4330127018922193*bzr[1]-0.4330127018922193*bzl[1])*c+0.75*psr[5]-0.75*psl[5]-0.4330127018922193*(psr[1]+psl[1]))*gamma; 
  incr[6] = (((-0.75*(bzr[6]+bzl[6]))+0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2])*c+0.75*psr[6]-0.75*psl[6]-0.4330127018922193*(psr[2]+psl[2]))*gamma; 
  incr[7] = (0.25*(bzl[7]*c+psr[7]+psl[7])-0.25*bzr[7]*c)*gamma; 
  incr[8] = (0.25*(bzl[8]*c+psr[8]+psl[8])-0.25*bzr[8]*c)*gamma; 
  incr[9] = (((-1.25*bzr[9])+1.25*bzl[9]+0.9682458365518543*(bzr[3]+bzl[3])-0.5590169943749475*bzr[0]+0.5590169943749475*bzl[0])*c+1.25*(psr[9]+psl[9])-0.9682458365518543*psr[3]+0.9682458365518543*psl[3]+0.5590169943749475*(psr[0]+psl[0]))*gamma; 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 
  outBzr[4] += incr[4]*dxr1; 
  outBzr[5] += incr[5]*dxr1; 
  outBzr[6] += incr[6]*dxr1; 
  outBzr[7] += incr[7]*dxr1; 
  outBzr[8] += incr[8]*dxr1; 
  outBzr[9] += incr[9]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += -1.0*incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 
  outBzl[3] += incr[3]*dxl1; 
  outBzl[4] += -1.0*incr[4]*dxl1; 
  outBzl[5] += incr[5]*dxl1; 
  outBzl[6] += incr[6]*dxl1; 
  outBzl[7] += -1.0*incr[7]*dxl1; 
  outBzl[8] += -1.0*incr[8]*dxl1; 
  outBzl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = (((-0.5590169943749475*phr[9])+0.5590169943749475*phl[9]+0.4330127018922193*(phr[3]+phl[3])-0.25*phr[0]+0.25*phl[0])*c+0.5590169943749475*(ezr[9]+ezl[9])-0.4330127018922193*ezr[3]+0.4330127018922193*ezl[3]+0.25*(ezr[0]+ezl[0]))*chi; 
  incr[1] = ((0.4330127018922193*(phr[5]+phl[5])-0.25*phr[1]+0.25*phl[1])*c-0.4330127018922193*ezr[5]+0.4330127018922193*ezl[5]+0.25*(ezr[1]+ezl[1]))*chi; 
  incr[2] = ((0.4330127018922193*(phr[6]+phl[6])-0.25*phr[2]+0.25*phl[2])*c-0.4330127018922193*ezr[6]+0.4330127018922193*ezl[6]+0.25*(ezr[2]+ezl[2]))*chi; 
  incr[3] = ((0.9682458365518543*phr[9]-0.9682458365518543*phl[9]-0.75*(phr[3]+phl[3])+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c-0.9682458365518543*(ezr[9]+ezl[9])+0.75*ezr[3]-0.75*ezl[3]-0.4330127018922193*(ezr[0]+ezl[0]))*chi; 
  incr[4] = (0.25*(phl[4]*c+ezr[4]+ezl[4])-0.25*phr[4]*c)*chi; 
  incr[5] = (((-0.75*(phr[5]+phl[5]))+0.4330127018922193*phr[1]-0.4330127018922193*phl[1])*c+0.75*ezr[5]-0.75*ezl[5]-0.4330127018922193*(ezr[1]+ezl[1]))*chi; 
  incr[6] = (((-0.75*(phr[6]+phl[6]))+0.4330127018922193*phr[2]-0.4330127018922193*phl[2])*c+0.75*ezr[6]-0.75*ezl[6]-0.4330127018922193*(ezr[2]+ezl[2]))*chi; 
  incr[7] = (0.25*(phl[7]*c+ezr[7]+ezl[7])-0.25*phr[7]*c)*chi; 
  incr[8] = (0.25*(phl[8]*c+ezr[8]+ezl[8])-0.25*phr[8]*c)*chi; 
  incr[9] = (((-1.25*phr[9])+1.25*phl[9]+0.9682458365518543*(phr[3]+phl[3])-0.5590169943749475*phr[0]+0.5590169943749475*phl[0])*c+1.25*(ezr[9]+ezl[9])-0.9682458365518543*ezr[3]+0.9682458365518543*ezl[3]+0.5590169943749475*(ezr[0]+ezl[0]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 
  outPhr[4] += incr[4]*dxr1; 
  outPhr[5] += incr[5]*dxr1; 
  outPhr[6] += incr[6]*dxr1; 
  outPhr[7] += incr[7]*dxr1; 
  outPhr[8] += incr[8]*dxr1; 
  outPhr[9] += incr[9]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += -1.0*incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 
  outPhl[3] += incr[3]*dxl1; 
  outPhl[4] += -1.0*incr[4]*dxl1; 
  outPhl[5] += incr[5]*dxl1; 
  outPhl[6] += incr[6]*dxl1; 
  outPhl[7] += -1.0*incr[7]*dxl1; 
  outPhl[8] += -1.0*incr[8]*dxl1; 
  outPhl[9] += -1.0*incr[9]*dxl1; 

 
  incr[0] = ((-0.5590169943749475*psr[9])+0.5590169943749475*psl[9]+0.4330127018922193*(psr[3]+psl[3])-0.25*psr[0]+0.25*psl[0])*c*gamma+(0.5590169943749475*(bzr[9]+bzl[9])-0.4330127018922193*bzr[3]+0.4330127018922193*bzl[3]+0.25*(bzr[0]+bzl[0]))*c2gamma; 
  incr[1] = (0.4330127018922193*(psr[5]+psl[5])-0.25*psr[1]+0.25*psl[1])*c*gamma+((-0.4330127018922193*bzr[5])+0.4330127018922193*bzl[5]+0.25*(bzr[1]+bzl[1]))*c2gamma; 
  incr[2] = (0.4330127018922193*(psr[6]+psl[6])-0.25*psr[2]+0.25*psl[2])*c*gamma+((-0.4330127018922193*bzr[6])+0.4330127018922193*bzl[6]+0.25*(bzr[2]+bzl[2]))*c2gamma; 
  incr[3] = (0.9682458365518543*psr[9]-0.9682458365518543*psl[9]-0.75*(psr[3]+psl[3])+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+((-0.9682458365518543*(bzr[9]+bzl[9]))+0.75*bzr[3]-0.75*bzl[3]-0.4330127018922193*(bzr[0]+bzl[0]))*c2gamma; 
  incr[4] = (0.25*psl[4]-0.25*psr[4])*c*gamma+0.25*(bzr[4]+bzl[4])*c2gamma; 
  incr[5] = ((-0.75*(psr[5]+psl[5]))+0.4330127018922193*psr[1]-0.4330127018922193*psl[1])*c*gamma+(0.75*bzr[5]-0.75*bzl[5]-0.4330127018922193*(bzr[1]+bzl[1]))*c2gamma; 
  incr[6] = ((-0.75*(psr[6]+psl[6]))+0.4330127018922193*psr[2]-0.4330127018922193*psl[2])*c*gamma+(0.75*bzr[6]-0.75*bzl[6]-0.4330127018922193*(bzr[2]+bzl[2]))*c2gamma; 
  incr[7] = (0.25*psl[7]-0.25*psr[7])*c*gamma+0.25*(bzr[7]+bzl[7])*c2gamma; 
  incr[8] = (0.25*psl[8]-0.25*psr[8])*c*gamma+0.25*(bzr[8]+bzl[8])*c2gamma; 
  incr[9] = ((-1.25*psr[9])+1.25*psl[9]+0.9682458365518543*(psr[3]+psl[3])-0.5590169943749475*psr[0]+0.5590169943749475*psl[0])*c*gamma+(1.25*(bzr[9]+bzl[9])-0.9682458365518543*bzr[3]+0.9682458365518543*bzl[3]+0.5590169943749475*(bzr[0]+bzl[0]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 
  outPsr[4] += incr[4]*dxr1; 
  outPsr[5] += incr[5]*dxr1; 
  outPsr[6] += incr[6]*dxr1; 
  outPsr[7] += incr[7]*dxr1; 
  outPsr[8] += incr[8]*dxr1; 
  outPsr[9] += incr[9]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += -1.0*incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 
  outPsl[3] += incr[3]*dxl1; 
  outPsl[4] += -1.0*incr[4]*dxl1; 
  outPsl[5] += incr[5]*dxl1; 
  outPsl[6] += incr[6]*dxl1; 
  outPsl[7] += -1.0*incr[7]*dxl1; 
  outPsl[8] += -1.0*incr[8]*dxl1; 
  outPsl[9] += -1.0*incr[9]*dxl1; 

 
  return c; 
} 
double MaxwellSurf3xMax_Z_P3(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[2]; 
  const double dxr1 = 2.0/dxr[2]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[20]; 
  const double *ezl = &ql[40]; 
  const double *bxl = &ql[60]; 
  const double *byl = &ql[80]; 
  const double *bzl = &ql[100]; 
  const double *phl = &ql[120]; 
  const double *psl = &ql[140]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[20]; 
  double *outEzl = &outl[40]; 
  double *outBxl = &outl[60]; 
  double *outByl = &outl[80]; 
  double *outBzl = &outl[100]; 
  double *outPhl = &outl[120]; 
  double *outPsl = &outl[140]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[20]; 
  const double *ezr = &qr[40]; 
  const double *bxr = &qr[60]; 
  const double *byr = &qr[80]; 
  const double *bzr = &qr[100]; 
  const double *phr = &qr[120]; 
  const double *psr = &qr[140]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[20]; 
  double *outEzr = &outr[40]; 
  double *outBxr = &outr[60]; 
  double *outByr = &outr[80]; 
  double *outBzr = &outr[100]; 
  double *outPhr = &outr[120]; 
  double *outPsr = &outr[140]; 
 
  double incr[20]; 
 
  incr[0] = (0.6614378277661477*(exr[19]+exl[19])-0.5590169943749475*exr[9]+0.5590169943749475*exl[9]+0.4330127018922193*(exr[3]+exl[3])-0.25*exr[0]+0.25*exl[0])*tau+((-0.6614378277661477*byr[19])+0.6614378277661477*byl[19]+0.5590169943749475*(byr[9]+byl[9])-0.4330127018922193*byr[3]+0.4330127018922193*byl[3]+0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = ((-0.5590169943749476*exr[15])+0.5590169943749476*exl[15]+0.4330127018922193*(exr[5]+exl[5])-0.25*exr[1]+0.25*exl[1])*tau+(0.5590169943749476*(byr[15]+byl[15])-0.4330127018922193*byr[5]+0.4330127018922193*byl[5]+0.25*(byr[1]+byl[1]))*c2; 
  incr[2] = ((-0.5590169943749476*exr[16])+0.5590169943749476*exl[16]+0.4330127018922193*(exr[6]+exl[6])-0.25*exr[2]+0.25*exl[2])*tau+(0.5590169943749476*(byr[16]+byl[16])-0.4330127018922193*byr[6]+0.4330127018922193*byl[6]+0.25*(byr[2]+byl[2]))*c2; 
  incr[3] = ((-1.14564392373896*(exr[19]+exl[19]))+0.9682458365518543*exr[9]-0.9682458365518543*exl[9]-0.75*(exr[3]+exl[3])+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*tau+(1.14564392373896*byr[19]-1.14564392373896*byl[19]-0.9682458365518543*(byr[9]+byl[9])+0.75*byr[3]-0.75*byl[3]-0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[4] = (0.4330127018922193*(exr[10]+exl[10])-0.25*exr[4]+0.25*exl[4])*tau+((-0.4330127018922193*byr[10])+0.4330127018922193*byl[10]+0.25*(byr[4]+byl[4]))*c2; 
  incr[5] = (0.9682458365518543*exr[15]-0.9682458365518543*exl[15]-0.75*(exr[5]+exl[5])+0.4330127018922193*exr[1]-0.4330127018922193*exl[1])*tau+((-0.9682458365518543*(byr[15]+byl[15]))+0.75*byr[5]-0.75*byl[5]-0.4330127018922193*(byr[1]+byl[1]))*c2; 
  incr[6] = (0.9682458365518543*exr[16]-0.9682458365518543*exl[16]-0.75*(exr[6]+exl[6])+0.4330127018922193*exr[2]-0.4330127018922193*exl[2])*tau+((-0.9682458365518543*(byr[16]+byl[16]))+0.75*byr[6]-0.75*byl[6]-0.4330127018922193*(byr[2]+byl[2]))*c2; 
  incr[7] = (0.4330127018922194*(exr[13]+exl[13])-0.25*exr[7]+0.25*exl[7])*tau+((-0.4330127018922194*byr[13])+0.4330127018922194*byl[13]+0.25*(byr[7]+byl[7]))*c2; 
  incr[8] = (0.4330127018922194*(exr[14]+exl[14])-0.25*exr[8]+0.25*exl[8])*tau+((-0.4330127018922194*byr[14])+0.4330127018922194*byl[14]+0.25*(byr[8]+byl[8]))*c2; 
  incr[9] = (1.479019945774904*(exr[19]+exl[19])-1.25*exr[9]+1.25*exl[9]+0.9682458365518543*(exr[3]+exl[3])-0.5590169943749475*exr[0]+0.5590169943749475*exl[0])*tau+((-1.479019945774904*byr[19])+1.479019945774904*byl[19]+1.25*(byr[9]+byl[9])-0.9682458365518543*byr[3]+0.9682458365518543*byl[3]+0.5590169943749475*(byr[0]+byl[0]))*c2; 
  incr[10] = ((-0.75*(exr[10]+exl[10]))+0.4330127018922193*exr[4]-0.4330127018922193*exl[4])*tau+(0.75*byr[10]-0.75*byl[10]-0.4330127018922193*(byr[4]+byl[4]))*c2; 
  incr[11] = 0.25*(exl[11]*tau+(byr[11]+byl[11])*c2)-0.25*exr[11]*tau; 
  incr[12] = 0.25*(exl[12]*tau+(byr[12]+byl[12])*c2)-0.25*exr[12]*tau; 
  incr[13] = ((-0.75*(exr[13]+exl[13]))+0.4330127018922194*exr[7]-0.4330127018922194*exl[7])*tau+(0.75*byr[13]-0.75*byl[13]-0.4330127018922194*(byr[7]+byl[7]))*c2; 
  incr[14] = ((-0.75*(exr[14]+exl[14]))+0.4330127018922194*exr[8]-0.4330127018922194*exl[8])*tau+(0.75*byr[14]-0.75*byl[14]-0.4330127018922194*(byr[8]+byl[8]))*c2; 
  incr[15] = ((-1.25*exr[15])+1.25*exl[15]+0.9682458365518543*(exr[5]+exl[5])-0.5590169943749476*exr[1]+0.5590169943749476*exl[1])*tau+(1.25*(byr[15]+byl[15])-0.9682458365518543*byr[5]+0.9682458365518543*byl[5]+0.5590169943749476*(byr[1]+byl[1]))*c2; 
  incr[16] = ((-1.25*exr[16])+1.25*exl[16]+0.9682458365518543*(exr[6]+exl[6])-0.5590169943749476*exr[2]+0.5590169943749476*exl[2])*tau+(1.25*(byr[16]+byl[16])-0.9682458365518543*byr[6]+0.9682458365518543*byl[6]+0.5590169943749476*(byr[2]+byl[2]))*c2; 
  incr[17] = 0.25*(exl[17]*tau+(byr[17]+byl[17])*c2)-0.25*exr[17]*tau; 
  incr[18] = 0.25*(exl[18]*tau+(byr[18]+byl[18])*c2)-0.25*exr[18]*tau; 
  incr[19] = ((-1.75*(exr[19]+exl[19]))+1.479019945774904*exr[9]-1.479019945774904*exl[9]-1.14564392373896*(exr[3]+exl[3])+0.6614378277661477*exr[0]-0.6614378277661477*exl[0])*tau+(1.75*byr[19]-1.75*byl[19]-1.479019945774904*(byr[9]+byl[9])+1.14564392373896*byr[3]-1.14564392373896*byl[3]-0.6614378277661477*(byr[0]+byl[0]))*c2; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 
  outExr[4] += incr[4]*dxr1; 
  outExr[5] += incr[5]*dxr1; 
  outExr[6] += incr[6]*dxr1; 
  outExr[7] += incr[7]*dxr1; 
  outExr[8] += incr[8]*dxr1; 
  outExr[9] += incr[9]*dxr1; 
  outExr[10] += incr[10]*dxr1; 
  outExr[11] += incr[11]*dxr1; 
  outExr[12] += incr[12]*dxr1; 
  outExr[13] += incr[13]*dxr1; 
  outExr[14] += incr[14]*dxr1; 
  outExr[15] += incr[15]*dxr1; 
  outExr[16] += incr[16]*dxr1; 
  outExr[17] += incr[17]*dxr1; 
  outExr[18] += incr[18]*dxr1; 
  outExr[19] += incr[19]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += -1.0*incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 
  outExl[3] += incr[3]*dxl1; 
  outExl[4] += -1.0*incr[4]*dxl1; 
  outExl[5] += incr[5]*dxl1; 
  outExl[6] += incr[6]*dxl1; 
  outExl[7] += -1.0*incr[7]*dxl1; 
  outExl[8] += -1.0*incr[8]*dxl1; 
  outExl[9] += -1.0*incr[9]*dxl1; 
  outExl[10] += incr[10]*dxl1; 
  outExl[11] += -1.0*incr[11]*dxl1; 
  outExl[12] += -1.0*incr[12]*dxl1; 
  outExl[13] += incr[13]*dxl1; 
  outExl[14] += incr[14]*dxl1; 
  outExl[15] += -1.0*incr[15]*dxl1; 
  outExl[16] += -1.0*incr[16]*dxl1; 
  outExl[17] += -1.0*incr[17]*dxl1; 
  outExl[18] += -1.0*incr[18]*dxl1; 
  outExl[19] += incr[19]*dxl1; 

 
  incr[0] = (0.6614378277661477*(eyr[19]+eyl[19])-0.5590169943749475*eyr[9]+0.5590169943749475*eyl[9]+0.4330127018922193*(eyr[3]+eyl[3])-0.25*eyr[0]+0.25*eyl[0])*tau+(0.6614378277661477*bxr[19]-0.6614378277661477*bxl[19]-0.5590169943749475*(bxr[9]+bxl[9])+0.4330127018922193*bxr[3]-0.4330127018922193*bxl[3]-0.25*(bxr[0]+bxl[0]))*c2; 
  incr[1] = ((-0.5590169943749476*eyr[15])+0.5590169943749476*eyl[15]+0.4330127018922193*(eyr[5]+eyl[5])-0.25*eyr[1]+0.25*eyl[1])*tau+((-0.5590169943749476*(bxr[15]+bxl[15]))+0.4330127018922193*bxr[5]-0.4330127018922193*bxl[5]-0.25*(bxr[1]+bxl[1]))*c2; 
  incr[2] = ((-0.5590169943749476*eyr[16])+0.5590169943749476*eyl[16]+0.4330127018922193*(eyr[6]+eyl[6])-0.25*eyr[2]+0.25*eyl[2])*tau+((-0.5590169943749476*(bxr[16]+bxl[16]))+0.4330127018922193*bxr[6]-0.4330127018922193*bxl[6]-0.25*(bxr[2]+bxl[2]))*c2; 
  incr[3] = ((-1.14564392373896*(eyr[19]+eyl[19]))+0.9682458365518543*eyr[9]-0.9682458365518543*eyl[9]-0.75*(eyr[3]+eyl[3])+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+((-1.14564392373896*bxr[19])+1.14564392373896*bxl[19]+0.9682458365518543*(bxr[9]+bxl[9])-0.75*bxr[3]+0.75*bxl[3]+0.4330127018922193*(bxr[0]+bxl[0]))*c2; 
  incr[4] = (0.4330127018922193*(eyr[10]+eyl[10])-0.25*eyr[4]+0.25*eyl[4])*tau+(0.4330127018922193*bxr[10]-0.4330127018922193*bxl[10]-0.25*(bxr[4]+bxl[4]))*c2; 
  incr[5] = (0.9682458365518543*eyr[15]-0.9682458365518543*eyl[15]-0.75*(eyr[5]+eyl[5])+0.4330127018922193*eyr[1]-0.4330127018922193*eyl[1])*tau+(0.9682458365518543*(bxr[15]+bxl[15])-0.75*bxr[5]+0.75*bxl[5]+0.4330127018922193*(bxr[1]+bxl[1]))*c2; 
  incr[6] = (0.9682458365518543*eyr[16]-0.9682458365518543*eyl[16]-0.75*(eyr[6]+eyl[6])+0.4330127018922193*eyr[2]-0.4330127018922193*eyl[2])*tau+(0.9682458365518543*(bxr[16]+bxl[16])-0.75*bxr[6]+0.75*bxl[6]+0.4330127018922193*(bxr[2]+bxl[2]))*c2; 
  incr[7] = (0.4330127018922194*(eyr[13]+eyl[13])-0.25*eyr[7]+0.25*eyl[7])*tau+(0.4330127018922194*bxr[13]-0.4330127018922194*bxl[13]-0.25*(bxr[7]+bxl[7]))*c2; 
  incr[8] = (0.4330127018922194*(eyr[14]+eyl[14])-0.25*eyr[8]+0.25*eyl[8])*tau+(0.4330127018922194*bxr[14]-0.4330127018922194*bxl[14]-0.25*(bxr[8]+bxl[8]))*c2; 
  incr[9] = (1.479019945774904*(eyr[19]+eyl[19])-1.25*eyr[9]+1.25*eyl[9]+0.9682458365518543*(eyr[3]+eyl[3])-0.5590169943749475*eyr[0]+0.5590169943749475*eyl[0])*tau+(1.479019945774904*bxr[19]-1.479019945774904*bxl[19]-1.25*(bxr[9]+bxl[9])+0.9682458365518543*bxr[3]-0.9682458365518543*bxl[3]-0.5590169943749475*(bxr[0]+bxl[0]))*c2; 
  incr[10] = ((-0.75*(eyr[10]+eyl[10]))+0.4330127018922193*eyr[4]-0.4330127018922193*eyl[4])*tau+((-0.75*bxr[10])+0.75*bxl[10]+0.4330127018922193*(bxr[4]+bxl[4]))*c2; 
  incr[11] = (0.25*eyl[11]-0.25*eyr[11])*tau-0.25*(bxr[11]+bxl[11])*c2; 
  incr[12] = (0.25*eyl[12]-0.25*eyr[12])*tau-0.25*(bxr[12]+bxl[12])*c2; 
  incr[13] = ((-0.75*(eyr[13]+eyl[13]))+0.4330127018922194*eyr[7]-0.4330127018922194*eyl[7])*tau+((-0.75*bxr[13])+0.75*bxl[13]+0.4330127018922194*(bxr[7]+bxl[7]))*c2; 
  incr[14] = ((-0.75*(eyr[14]+eyl[14]))+0.4330127018922194*eyr[8]-0.4330127018922194*eyl[8])*tau+((-0.75*bxr[14])+0.75*bxl[14]+0.4330127018922194*(bxr[8]+bxl[8]))*c2; 
  incr[15] = ((-1.25*eyr[15])+1.25*eyl[15]+0.9682458365518543*(eyr[5]+eyl[5])-0.5590169943749476*eyr[1]+0.5590169943749476*eyl[1])*tau+((-1.25*(bxr[15]+bxl[15]))+0.9682458365518543*bxr[5]-0.9682458365518543*bxl[5]-0.5590169943749476*(bxr[1]+bxl[1]))*c2; 
  incr[16] = ((-1.25*eyr[16])+1.25*eyl[16]+0.9682458365518543*(eyr[6]+eyl[6])-0.5590169943749476*eyr[2]+0.5590169943749476*eyl[2])*tau+((-1.25*(bxr[16]+bxl[16]))+0.9682458365518543*bxr[6]-0.9682458365518543*bxl[6]-0.5590169943749476*(bxr[2]+bxl[2]))*c2; 
  incr[17] = (0.25*eyl[17]-0.25*eyr[17])*tau-0.25*(bxr[17]+bxl[17])*c2; 
  incr[18] = (0.25*eyl[18]-0.25*eyr[18])*tau-0.25*(bxr[18]+bxl[18])*c2; 
  incr[19] = ((-1.75*(eyr[19]+eyl[19]))+1.479019945774904*eyr[9]-1.479019945774904*eyl[9]-1.14564392373896*(eyr[3]+eyl[3])+0.6614378277661477*eyr[0]-0.6614378277661477*eyl[0])*tau+((-1.75*bxr[19])+1.75*bxl[19]+1.479019945774904*(bxr[9]+bxl[9])-1.14564392373896*bxr[3]+1.14564392373896*bxl[3]+0.6614378277661477*(bxr[0]+bxl[0]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 
  outEyr[4] += incr[4]*dxr1; 
  outEyr[5] += incr[5]*dxr1; 
  outEyr[6] += incr[6]*dxr1; 
  outEyr[7] += incr[7]*dxr1; 
  outEyr[8] += incr[8]*dxr1; 
  outEyr[9] += incr[9]*dxr1; 
  outEyr[10] += incr[10]*dxr1; 
  outEyr[11] += incr[11]*dxr1; 
  outEyr[12] += incr[12]*dxr1; 
  outEyr[13] += incr[13]*dxr1; 
  outEyr[14] += incr[14]*dxr1; 
  outEyr[15] += incr[15]*dxr1; 
  outEyr[16] += incr[16]*dxr1; 
  outEyr[17] += incr[17]*dxr1; 
  outEyr[18] += incr[18]*dxr1; 
  outEyr[19] += incr[19]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += -1.0*incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 
  outEyl[3] += incr[3]*dxl1; 
  outEyl[4] += -1.0*incr[4]*dxl1; 
  outEyl[5] += incr[5]*dxl1; 
  outEyl[6] += incr[6]*dxl1; 
  outEyl[7] += -1.0*incr[7]*dxl1; 
  outEyl[8] += -1.0*incr[8]*dxl1; 
  outEyl[9] += -1.0*incr[9]*dxl1; 
  outEyl[10] += incr[10]*dxl1; 
  outEyl[11] += -1.0*incr[11]*dxl1; 
  outEyl[12] += -1.0*incr[12]*dxl1; 
  outEyl[13] += incr[13]*dxl1; 
  outEyl[14] += incr[14]*dxl1; 
  outEyl[15] += -1.0*incr[15]*dxl1; 
  outEyl[16] += -1.0*incr[16]*dxl1; 
  outEyl[17] += -1.0*incr[17]*dxl1; 
  outEyl[18] += -1.0*incr[18]*dxl1; 
  outEyl[19] += incr[19]*dxl1; 

 
  incr[0] = (0.6614378277661477*(ezr[19]+ezl[19])-0.5590169943749475*ezr[9]+0.5590169943749475*ezl[9]+0.4330127018922193*(ezr[3]+ezl[3])-0.25*ezr[0]+0.25*ezl[0])*c*chi+((-0.6614378277661477*phr[19])+0.6614378277661477*phl[19]+0.5590169943749475*(phr[9]+phl[9])-0.4330127018922193*phr[3]+0.4330127018922193*phl[3]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = ((-0.5590169943749476*ezr[15])+0.5590169943749476*ezl[15]+0.4330127018922193*(ezr[5]+ezl[5])-0.25*ezr[1]+0.25*ezl[1])*c*chi+(0.5590169943749476*(phr[15]+phl[15])-0.4330127018922193*phr[5]+0.4330127018922193*phl[5]+0.25*(phr[1]+phl[1]))*c2chi; 
  incr[2] = ((-0.5590169943749476*ezr[16])+0.5590169943749476*ezl[16]+0.4330127018922193*(ezr[6]+ezl[6])-0.25*ezr[2]+0.25*ezl[2])*c*chi+(0.5590169943749476*(phr[16]+phl[16])-0.4330127018922193*phr[6]+0.4330127018922193*phl[6]+0.25*(phr[2]+phl[2]))*c2chi; 
  incr[3] = ((-1.14564392373896*(ezr[19]+ezl[19]))+0.9682458365518543*ezr[9]-0.9682458365518543*ezl[9]-0.75*(ezr[3]+ezl[3])+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*c*chi+(1.14564392373896*phr[19]-1.14564392373896*phl[19]-0.9682458365518543*(phr[9]+phl[9])+0.75*phr[3]-0.75*phl[3]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[4] = (0.4330127018922193*(ezr[10]+ezl[10])-0.25*ezr[4]+0.25*ezl[4])*c*chi+((-0.4330127018922193*phr[10])+0.4330127018922193*phl[10]+0.25*(phr[4]+phl[4]))*c2chi; 
  incr[5] = (0.9682458365518543*ezr[15]-0.9682458365518543*ezl[15]-0.75*(ezr[5]+ezl[5])+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1])*c*chi+((-0.9682458365518543*(phr[15]+phl[15]))+0.75*phr[5]-0.75*phl[5]-0.4330127018922193*(phr[1]+phl[1]))*c2chi; 
  incr[6] = (0.9682458365518543*ezr[16]-0.9682458365518543*ezl[16]-0.75*(ezr[6]+ezl[6])+0.4330127018922193*ezr[2]-0.4330127018922193*ezl[2])*c*chi+((-0.9682458365518543*(phr[16]+phl[16]))+0.75*phr[6]-0.75*phl[6]-0.4330127018922193*(phr[2]+phl[2]))*c2chi; 
  incr[7] = (0.4330127018922194*(ezr[13]+ezl[13])-0.25*ezr[7]+0.25*ezl[7])*c*chi+((-0.4330127018922194*phr[13])+0.4330127018922194*phl[13]+0.25*(phr[7]+phl[7]))*c2chi; 
  incr[8] = (0.4330127018922194*(ezr[14]+ezl[14])-0.25*ezr[8]+0.25*ezl[8])*c*chi+((-0.4330127018922194*phr[14])+0.4330127018922194*phl[14]+0.25*(phr[8]+phl[8]))*c2chi; 
  incr[9] = (1.479019945774904*(ezr[19]+ezl[19])-1.25*ezr[9]+1.25*ezl[9]+0.9682458365518543*(ezr[3]+ezl[3])-0.5590169943749475*ezr[0]+0.5590169943749475*ezl[0])*c*chi+((-1.479019945774904*phr[19])+1.479019945774904*phl[19]+1.25*(phr[9]+phl[9])-0.9682458365518543*phr[3]+0.9682458365518543*phl[3]+0.5590169943749475*(phr[0]+phl[0]))*c2chi; 
  incr[10] = ((-0.75*(ezr[10]+ezl[10]))+0.4330127018922193*ezr[4]-0.4330127018922193*ezl[4])*c*chi+(0.75*phr[10]-0.75*phl[10]-0.4330127018922193*(phr[4]+phl[4]))*c2chi; 
  incr[11] = (0.25*ezl[11]-0.25*ezr[11])*c*chi+0.25*(phr[11]+phl[11])*c2chi; 
  incr[12] = (0.25*ezl[12]-0.25*ezr[12])*c*chi+0.25*(phr[12]+phl[12])*c2chi; 
  incr[13] = ((-0.75*(ezr[13]+ezl[13]))+0.4330127018922194*ezr[7]-0.4330127018922194*ezl[7])*c*chi+(0.75*phr[13]-0.75*phl[13]-0.4330127018922194*(phr[7]+phl[7]))*c2chi; 
  incr[14] = ((-0.75*(ezr[14]+ezl[14]))+0.4330127018922194*ezr[8]-0.4330127018922194*ezl[8])*c*chi+(0.75*phr[14]-0.75*phl[14]-0.4330127018922194*(phr[8]+phl[8]))*c2chi; 
  incr[15] = ((-1.25*ezr[15])+1.25*ezl[15]+0.9682458365518543*(ezr[5]+ezl[5])-0.5590169943749476*ezr[1]+0.5590169943749476*ezl[1])*c*chi+(1.25*(phr[15]+phl[15])-0.9682458365518543*phr[5]+0.9682458365518543*phl[5]+0.5590169943749476*(phr[1]+phl[1]))*c2chi; 
  incr[16] = ((-1.25*ezr[16])+1.25*ezl[16]+0.9682458365518543*(ezr[6]+ezl[6])-0.5590169943749476*ezr[2]+0.5590169943749476*ezl[2])*c*chi+(1.25*(phr[16]+phl[16])-0.9682458365518543*phr[6]+0.9682458365518543*phl[6]+0.5590169943749476*(phr[2]+phl[2]))*c2chi; 
  incr[17] = (0.25*ezl[17]-0.25*ezr[17])*c*chi+0.25*(phr[17]+phl[17])*c2chi; 
  incr[18] = (0.25*ezl[18]-0.25*ezr[18])*c*chi+0.25*(phr[18]+phl[18])*c2chi; 
  incr[19] = ((-1.75*(ezr[19]+ezl[19]))+1.479019945774904*ezr[9]-1.479019945774904*ezl[9]-1.14564392373896*(ezr[3]+ezl[3])+0.6614378277661477*ezr[0]-0.6614378277661477*ezl[0])*c*chi+(1.75*phr[19]-1.75*phl[19]-1.479019945774904*(phr[9]+phl[9])+1.14564392373896*phr[3]-1.14564392373896*phl[3]-0.6614378277661477*(phr[0]+phl[0]))*c2chi; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 
  outEzr[4] += incr[4]*dxr1; 
  outEzr[5] += incr[5]*dxr1; 
  outEzr[6] += incr[6]*dxr1; 
  outEzr[7] += incr[7]*dxr1; 
  outEzr[8] += incr[8]*dxr1; 
  outEzr[9] += incr[9]*dxr1; 
  outEzr[10] += incr[10]*dxr1; 
  outEzr[11] += incr[11]*dxr1; 
  outEzr[12] += incr[12]*dxr1; 
  outEzr[13] += incr[13]*dxr1; 
  outEzr[14] += incr[14]*dxr1; 
  outEzr[15] += incr[15]*dxr1; 
  outEzr[16] += incr[16]*dxr1; 
  outEzr[17] += incr[17]*dxr1; 
  outEzr[18] += incr[18]*dxr1; 
  outEzr[19] += incr[19]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += -1.0*incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 
  outEzl[3] += incr[3]*dxl1; 
  outEzl[4] += -1.0*incr[4]*dxl1; 
  outEzl[5] += incr[5]*dxl1; 
  outEzl[6] += incr[6]*dxl1; 
  outEzl[7] += -1.0*incr[7]*dxl1; 
  outEzl[8] += -1.0*incr[8]*dxl1; 
  outEzl[9] += -1.0*incr[9]*dxl1; 
  outEzl[10] += incr[10]*dxl1; 
  outEzl[11] += -1.0*incr[11]*dxl1; 
  outEzl[12] += -1.0*incr[12]*dxl1; 
  outEzl[13] += incr[13]*dxl1; 
  outEzl[14] += incr[14]*dxl1; 
  outEzl[15] += -1.0*incr[15]*dxl1; 
  outEzl[16] += -1.0*incr[16]*dxl1; 
  outEzl[17] += -1.0*incr[17]*dxl1; 
  outEzl[18] += -1.0*incr[18]*dxl1; 
  outEzl[19] += incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(bxr[19]+bxl[19])-0.5590169943749475*bxr[9]+0.5590169943749475*bxl[9]+0.4330127018922193*(bxr[3]+bxl[3])-0.25*bxr[0]+0.25*bxl[0])*c2)/tau+0.6614378277661477*eyr[19]-0.6614378277661477*eyl[19]-0.5590169943749475*(eyr[9]+eyl[9])+0.4330127018922193*eyr[3]-0.4330127018922193*eyl[3]-0.25*(eyr[0]+eyl[0]); 
  incr[1] = (((-0.5590169943749476*bxr[15])+0.5590169943749476*bxl[15]+0.4330127018922193*(bxr[5]+bxl[5])-0.25*bxr[1]+0.25*bxl[1])*c2)/tau-0.5590169943749476*(eyr[15]+eyl[15])+0.4330127018922193*eyr[5]-0.4330127018922193*eyl[5]-0.25*(eyr[1]+eyl[1]); 
  incr[2] = (((-0.5590169943749476*bxr[16])+0.5590169943749476*bxl[16]+0.4330127018922193*(bxr[6]+bxl[6])-0.25*bxr[2]+0.25*bxl[2])*c2)/tau-0.5590169943749476*(eyr[16]+eyl[16])+0.4330127018922193*eyr[6]-0.4330127018922193*eyl[6]-0.25*(eyr[2]+eyl[2]); 
  incr[3] = (((-1.14564392373896*(bxr[19]+bxl[19]))+0.9682458365518543*bxr[9]-0.9682458365518543*bxl[9]-0.75*(bxr[3]+bxl[3])+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c2)/tau-1.14564392373896*eyr[19]+1.14564392373896*eyl[19]+0.9682458365518543*(eyr[9]+eyl[9])-0.75*eyr[3]+0.75*eyl[3]+0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[4] = ((0.4330127018922193*(bxr[10]+bxl[10])-0.25*bxr[4]+0.25*bxl[4])*c2)/tau+0.4330127018922193*eyr[10]-0.4330127018922193*eyl[10]-0.25*(eyr[4]+eyl[4]); 
  incr[5] = ((0.9682458365518543*bxr[15]-0.9682458365518543*bxl[15]-0.75*(bxr[5]+bxl[5])+0.4330127018922193*bxr[1]-0.4330127018922193*bxl[1])*c2)/tau+0.9682458365518543*(eyr[15]+eyl[15])-0.75*eyr[5]+0.75*eyl[5]+0.4330127018922193*(eyr[1]+eyl[1]); 
  incr[6] = ((0.9682458365518543*bxr[16]-0.9682458365518543*bxl[16]-0.75*(bxr[6]+bxl[6])+0.4330127018922193*bxr[2]-0.4330127018922193*bxl[2])*c2)/tau+0.9682458365518543*(eyr[16]+eyl[16])-0.75*eyr[6]+0.75*eyl[6]+0.4330127018922193*(eyr[2]+eyl[2]); 
  incr[7] = ((0.4330127018922194*(bxr[13]+bxl[13])-0.25*bxr[7]+0.25*bxl[7])*c2)/tau+0.4330127018922194*eyr[13]-0.4330127018922194*eyl[13]-0.25*(eyr[7]+eyl[7]); 
  incr[8] = ((0.4330127018922194*(bxr[14]+bxl[14])-0.25*bxr[8]+0.25*bxl[8])*c2)/tau+0.4330127018922194*eyr[14]-0.4330127018922194*eyl[14]-0.25*(eyr[8]+eyl[8]); 
  incr[9] = ((1.479019945774904*(bxr[19]+bxl[19])-1.25*bxr[9]+1.25*bxl[9]+0.9682458365518543*(bxr[3]+bxl[3])-0.5590169943749475*bxr[0]+0.5590169943749475*bxl[0])*c2)/tau+1.479019945774904*eyr[19]-1.479019945774904*eyl[19]-1.25*(eyr[9]+eyl[9])+0.9682458365518543*eyr[3]-0.9682458365518543*eyl[3]-0.5590169943749475*(eyr[0]+eyl[0]); 
  incr[10] = (((-0.75*(bxr[10]+bxl[10]))+0.4330127018922193*bxr[4]-0.4330127018922193*bxl[4])*c2)/tau-0.75*eyr[10]+0.75*eyl[10]+0.4330127018922193*(eyr[4]+eyl[4]); 
  incr[11] = ((0.25*bxl[11]-0.25*bxr[11])*c2)/tau-0.25*(eyr[11]+eyl[11]); 
  incr[12] = ((0.25*bxl[12]-0.25*bxr[12])*c2)/tau-0.25*(eyr[12]+eyl[12]); 
  incr[13] = (((-0.75*(bxr[13]+bxl[13]))+0.4330127018922194*bxr[7]-0.4330127018922194*bxl[7])*c2)/tau-0.75*eyr[13]+0.75*eyl[13]+0.4330127018922194*(eyr[7]+eyl[7]); 
  incr[14] = (((-0.75*(bxr[14]+bxl[14]))+0.4330127018922194*bxr[8]-0.4330127018922194*bxl[8])*c2)/tau-0.75*eyr[14]+0.75*eyl[14]+0.4330127018922194*(eyr[8]+eyl[8]); 
  incr[15] = (((-1.25*bxr[15])+1.25*bxl[15]+0.9682458365518543*(bxr[5]+bxl[5])-0.5590169943749476*bxr[1]+0.5590169943749476*bxl[1])*c2)/tau-1.25*(eyr[15]+eyl[15])+0.9682458365518543*eyr[5]-0.9682458365518543*eyl[5]-0.5590169943749476*(eyr[1]+eyl[1]); 
  incr[16] = (((-1.25*bxr[16])+1.25*bxl[16]+0.9682458365518543*(bxr[6]+bxl[6])-0.5590169943749476*bxr[2]+0.5590169943749476*bxl[2])*c2)/tau-1.25*(eyr[16]+eyl[16])+0.9682458365518543*eyr[6]-0.9682458365518543*eyl[6]-0.5590169943749476*(eyr[2]+eyl[2]); 
  incr[17] = ((0.25*bxl[17]-0.25*bxr[17])*c2)/tau-0.25*(eyr[17]+eyl[17]); 
  incr[18] = ((0.25*bxl[18]-0.25*bxr[18])*c2)/tau-0.25*(eyr[18]+eyl[18]); 
  incr[19] = (((-1.75*(bxr[19]+bxl[19]))+1.479019945774904*bxr[9]-1.479019945774904*bxl[9]-1.14564392373896*(bxr[3]+bxl[3])+0.6614378277661477*bxr[0]-0.6614378277661477*bxl[0])*c2)/tau-1.75*eyr[19]+1.75*eyl[19]+1.479019945774904*(eyr[9]+eyl[9])-1.14564392373896*eyr[3]+1.14564392373896*eyl[3]+0.6614378277661477*(eyr[0]+eyl[0]); 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 
  outBxr[4] += incr[4]*dxr1; 
  outBxr[5] += incr[5]*dxr1; 
  outBxr[6] += incr[6]*dxr1; 
  outBxr[7] += incr[7]*dxr1; 
  outBxr[8] += incr[8]*dxr1; 
  outBxr[9] += incr[9]*dxr1; 
  outBxr[10] += incr[10]*dxr1; 
  outBxr[11] += incr[11]*dxr1; 
  outBxr[12] += incr[12]*dxr1; 
  outBxr[13] += incr[13]*dxr1; 
  outBxr[14] += incr[14]*dxr1; 
  outBxr[15] += incr[15]*dxr1; 
  outBxr[16] += incr[16]*dxr1; 
  outBxr[17] += incr[17]*dxr1; 
  outBxr[18] += incr[18]*dxr1; 
  outBxr[19] += incr[19]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += -1.0*incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 
  outBxl[3] += incr[3]*dxl1; 
  outBxl[4] += -1.0*incr[4]*dxl1; 
  outBxl[5] += incr[5]*dxl1; 
  outBxl[6] += incr[6]*dxl1; 
  outBxl[7] += -1.0*incr[7]*dxl1; 
  outBxl[8] += -1.0*incr[8]*dxl1; 
  outBxl[9] += -1.0*incr[9]*dxl1; 
  outBxl[10] += incr[10]*dxl1; 
  outBxl[11] += -1.0*incr[11]*dxl1; 
  outBxl[12] += -1.0*incr[12]*dxl1; 
  outBxl[13] += incr[13]*dxl1; 
  outBxl[14] += incr[14]*dxl1; 
  outBxl[15] += -1.0*incr[15]*dxl1; 
  outBxl[16] += -1.0*incr[16]*dxl1; 
  outBxl[17] += -1.0*incr[17]*dxl1; 
  outBxl[18] += -1.0*incr[18]*dxl1; 
  outBxl[19] += incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(byr[19]+byl[19])-0.5590169943749475*byr[9]+0.5590169943749475*byl[9]+0.4330127018922193*(byr[3]+byl[3])-0.25*byr[0]+0.25*byl[0])*c2)/tau-0.6614378277661477*exr[19]+0.6614378277661477*exl[19]+0.5590169943749475*(exr[9]+exl[9])-0.4330127018922193*exr[3]+0.4330127018922193*exl[3]+0.25*(exr[0]+exl[0]); 
  incr[1] = (((-0.5590169943749476*byr[15])+0.5590169943749476*byl[15]+0.4330127018922193*(byr[5]+byl[5])-0.25*byr[1]+0.25*byl[1])*c2)/tau+0.5590169943749476*(exr[15]+exl[15])-0.4330127018922193*exr[5]+0.4330127018922193*exl[5]+0.25*(exr[1]+exl[1]); 
  incr[2] = (((-0.5590169943749476*byr[16])+0.5590169943749476*byl[16]+0.4330127018922193*(byr[6]+byl[6])-0.25*byr[2]+0.25*byl[2])*c2)/tau+0.5590169943749476*(exr[16]+exl[16])-0.4330127018922193*exr[6]+0.4330127018922193*exl[6]+0.25*(exr[2]+exl[2]); 
  incr[3] = (((-1.14564392373896*(byr[19]+byl[19]))+0.9682458365518543*byr[9]-0.9682458365518543*byl[9]-0.75*(byr[3]+byl[3])+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau+1.14564392373896*exr[19]-1.14564392373896*exl[19]-0.9682458365518543*(exr[9]+exl[9])+0.75*exr[3]-0.75*exl[3]-0.4330127018922193*(exr[0]+exl[0]); 
  incr[4] = ((0.4330127018922193*(byr[10]+byl[10])-0.25*byr[4]+0.25*byl[4])*c2)/tau-0.4330127018922193*exr[10]+0.4330127018922193*exl[10]+0.25*(exr[4]+exl[4]); 
  incr[5] = ((0.9682458365518543*byr[15]-0.9682458365518543*byl[15]-0.75*(byr[5]+byl[5])+0.4330127018922193*byr[1]-0.4330127018922193*byl[1])*c2)/tau-0.9682458365518543*(exr[15]+exl[15])+0.75*exr[5]-0.75*exl[5]-0.4330127018922193*(exr[1]+exl[1]); 
  incr[6] = ((0.9682458365518543*byr[16]-0.9682458365518543*byl[16]-0.75*(byr[6]+byl[6])+0.4330127018922193*byr[2]-0.4330127018922193*byl[2])*c2)/tau-0.9682458365518543*(exr[16]+exl[16])+0.75*exr[6]-0.75*exl[6]-0.4330127018922193*(exr[2]+exl[2]); 
  incr[7] = ((0.4330127018922194*(byr[13]+byl[13])-0.25*byr[7]+0.25*byl[7])*c2)/tau-0.4330127018922194*exr[13]+0.4330127018922194*exl[13]+0.25*(exr[7]+exl[7]); 
  incr[8] = ((0.4330127018922194*(byr[14]+byl[14])-0.25*byr[8]+0.25*byl[8])*c2)/tau-0.4330127018922194*exr[14]+0.4330127018922194*exl[14]+0.25*(exr[8]+exl[8]); 
  incr[9] = ((1.479019945774904*(byr[19]+byl[19])-1.25*byr[9]+1.25*byl[9]+0.9682458365518543*(byr[3]+byl[3])-0.5590169943749475*byr[0]+0.5590169943749475*byl[0])*c2)/tau-1.479019945774904*exr[19]+1.479019945774904*exl[19]+1.25*(exr[9]+exl[9])-0.9682458365518543*exr[3]+0.9682458365518543*exl[3]+0.5590169943749475*(exr[0]+exl[0]); 
  incr[10] = (((-0.75*(byr[10]+byl[10]))+0.4330127018922193*byr[4]-0.4330127018922193*byl[4])*c2)/tau+0.75*exr[10]-0.75*exl[10]-0.4330127018922193*(exr[4]+exl[4]); 
  incr[11] = ((0.25*byl[11]-0.25*byr[11])*c2)/tau+0.25*(exr[11]+exl[11]); 
  incr[12] = ((0.25*byl[12]-0.25*byr[12])*c2)/tau+0.25*(exr[12]+exl[12]); 
  incr[13] = (((-0.75*(byr[13]+byl[13]))+0.4330127018922194*byr[7]-0.4330127018922194*byl[7])*c2)/tau+0.75*exr[13]-0.75*exl[13]-0.4330127018922194*(exr[7]+exl[7]); 
  incr[14] = (((-0.75*(byr[14]+byl[14]))+0.4330127018922194*byr[8]-0.4330127018922194*byl[8])*c2)/tau+0.75*exr[14]-0.75*exl[14]-0.4330127018922194*(exr[8]+exl[8]); 
  incr[15] = (((-1.25*byr[15])+1.25*byl[15]+0.9682458365518543*(byr[5]+byl[5])-0.5590169943749476*byr[1]+0.5590169943749476*byl[1])*c2)/tau+1.25*(exr[15]+exl[15])-0.9682458365518543*exr[5]+0.9682458365518543*exl[5]+0.5590169943749476*(exr[1]+exl[1]); 
  incr[16] = (((-1.25*byr[16])+1.25*byl[16]+0.9682458365518543*(byr[6]+byl[6])-0.5590169943749476*byr[2]+0.5590169943749476*byl[2])*c2)/tau+1.25*(exr[16]+exl[16])-0.9682458365518543*exr[6]+0.9682458365518543*exl[6]+0.5590169943749476*(exr[2]+exl[2]); 
  incr[17] = ((0.25*byl[17]-0.25*byr[17])*c2)/tau+0.25*(exr[17]+exl[17]); 
  incr[18] = ((0.25*byl[18]-0.25*byr[18])*c2)/tau+0.25*(exr[18]+exl[18]); 
  incr[19] = (((-1.75*(byr[19]+byl[19]))+1.479019945774904*byr[9]-1.479019945774904*byl[9]-1.14564392373896*(byr[3]+byl[3])+0.6614378277661477*byr[0]-0.6614378277661477*byl[0])*c2)/tau+1.75*exr[19]-1.75*exl[19]-1.479019945774904*(exr[9]+exl[9])+1.14564392373896*exr[3]-1.14564392373896*exl[3]-0.6614378277661477*(exr[0]+exl[0]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 
  outByr[4] += incr[4]*dxr1; 
  outByr[5] += incr[5]*dxr1; 
  outByr[6] += incr[6]*dxr1; 
  outByr[7] += incr[7]*dxr1; 
  outByr[8] += incr[8]*dxr1; 
  outByr[9] += incr[9]*dxr1; 
  outByr[10] += incr[10]*dxr1; 
  outByr[11] += incr[11]*dxr1; 
  outByr[12] += incr[12]*dxr1; 
  outByr[13] += incr[13]*dxr1; 
  outByr[14] += incr[14]*dxr1; 
  outByr[15] += incr[15]*dxr1; 
  outByr[16] += incr[16]*dxr1; 
  outByr[17] += incr[17]*dxr1; 
  outByr[18] += incr[18]*dxr1; 
  outByr[19] += incr[19]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += -1.0*incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 
  outByl[3] += incr[3]*dxl1; 
  outByl[4] += -1.0*incr[4]*dxl1; 
  outByl[5] += incr[5]*dxl1; 
  outByl[6] += incr[6]*dxl1; 
  outByl[7] += -1.0*incr[7]*dxl1; 
  outByl[8] += -1.0*incr[8]*dxl1; 
  outByl[9] += -1.0*incr[9]*dxl1; 
  outByl[10] += incr[10]*dxl1; 
  outByl[11] += -1.0*incr[11]*dxl1; 
  outByl[12] += -1.0*incr[12]*dxl1; 
  outByl[13] += incr[13]*dxl1; 
  outByl[14] += incr[14]*dxl1; 
  outByl[15] += -1.0*incr[15]*dxl1; 
  outByl[16] += -1.0*incr[16]*dxl1; 
  outByl[17] += -1.0*incr[17]*dxl1; 
  outByl[18] += -1.0*incr[18]*dxl1; 
  outByl[19] += incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(bzr[19]+bzl[19])-0.5590169943749475*bzr[9]+0.5590169943749475*bzl[9]+0.4330127018922193*(bzr[3]+bzl[3])-0.25*bzr[0]+0.25*bzl[0])*c-0.6614378277661477*psr[19]+0.6614378277661477*psl[19]+0.5590169943749475*(psr[9]+psl[9])-0.4330127018922193*psr[3]+0.4330127018922193*psl[3]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (((-0.5590169943749476*bzr[15])+0.5590169943749476*bzl[15]+0.4330127018922193*(bzr[5]+bzl[5])-0.25*bzr[1]+0.25*bzl[1])*c+0.5590169943749476*(psr[15]+psl[15])-0.4330127018922193*psr[5]+0.4330127018922193*psl[5]+0.25*(psr[1]+psl[1]))*gamma; 
  incr[2] = (((-0.5590169943749476*bzr[16])+0.5590169943749476*bzl[16]+0.4330127018922193*(bzr[6]+bzl[6])-0.25*bzr[2]+0.25*bzl[2])*c+0.5590169943749476*(psr[16]+psl[16])-0.4330127018922193*psr[6]+0.4330127018922193*psl[6]+0.25*(psr[2]+psl[2]))*gamma; 
  incr[3] = (((-1.14564392373896*(bzr[19]+bzl[19]))+0.9682458365518543*bzr[9]-0.9682458365518543*bzl[9]-0.75*(bzr[3]+bzl[3])+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c+1.14564392373896*psr[19]-1.14564392373896*psl[19]-0.9682458365518543*(psr[9]+psl[9])+0.75*psr[3]-0.75*psl[3]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[4] = ((0.4330127018922193*(bzr[10]+bzl[10])-0.25*bzr[4]+0.25*bzl[4])*c-0.4330127018922193*psr[10]+0.4330127018922193*psl[10]+0.25*(psr[4]+psl[4]))*gamma; 
  incr[5] = ((0.9682458365518543*bzr[15]-0.9682458365518543*bzl[15]-0.75*(bzr[5]+bzl[5])+0.4330127018922193*bzr[1]-0.4330127018922193*bzl[1])*c-0.9682458365518543*(psr[15]+psl[15])+0.75*psr[5]-0.75*psl[5]-0.4330127018922193*(psr[1]+psl[1]))*gamma; 
  incr[6] = ((0.9682458365518543*bzr[16]-0.9682458365518543*bzl[16]-0.75*(bzr[6]+bzl[6])+0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2])*c-0.9682458365518543*(psr[16]+psl[16])+0.75*psr[6]-0.75*psl[6]-0.4330127018922193*(psr[2]+psl[2]))*gamma; 
  incr[7] = ((0.4330127018922194*(bzr[13]+bzl[13])-0.25*bzr[7]+0.25*bzl[7])*c-0.4330127018922194*psr[13]+0.4330127018922194*psl[13]+0.25*(psr[7]+psl[7]))*gamma; 
  incr[8] = ((0.4330127018922194*(bzr[14]+bzl[14])-0.25*bzr[8]+0.25*bzl[8])*c-0.4330127018922194*psr[14]+0.4330127018922194*psl[14]+0.25*(psr[8]+psl[8]))*gamma; 
  incr[9] = ((1.479019945774904*(bzr[19]+bzl[19])-1.25*bzr[9]+1.25*bzl[9]+0.9682458365518543*(bzr[3]+bzl[3])-0.5590169943749475*bzr[0]+0.5590169943749475*bzl[0])*c-1.479019945774904*psr[19]+1.479019945774904*psl[19]+1.25*(psr[9]+psl[9])-0.9682458365518543*psr[3]+0.9682458365518543*psl[3]+0.5590169943749475*(psr[0]+psl[0]))*gamma; 
  incr[10] = (((-0.75*(bzr[10]+bzl[10]))+0.4330127018922193*bzr[4]-0.4330127018922193*bzl[4])*c+0.75*psr[10]-0.75*psl[10]-0.4330127018922193*(psr[4]+psl[4]))*gamma; 
  incr[11] = (0.25*(bzl[11]*c+psr[11]+psl[11])-0.25*bzr[11]*c)*gamma; 
  incr[12] = (0.25*(bzl[12]*c+psr[12]+psl[12])-0.25*bzr[12]*c)*gamma; 
  incr[13] = (((-0.75*(bzr[13]+bzl[13]))+0.4330127018922194*bzr[7]-0.4330127018922194*bzl[7])*c+0.75*psr[13]-0.75*psl[13]-0.4330127018922194*(psr[7]+psl[7]))*gamma; 
  incr[14] = (((-0.75*(bzr[14]+bzl[14]))+0.4330127018922194*bzr[8]-0.4330127018922194*bzl[8])*c+0.75*psr[14]-0.75*psl[14]-0.4330127018922194*(psr[8]+psl[8]))*gamma; 
  incr[15] = (((-1.25*bzr[15])+1.25*bzl[15]+0.9682458365518543*(bzr[5]+bzl[5])-0.5590169943749476*bzr[1]+0.5590169943749476*bzl[1])*c+1.25*(psr[15]+psl[15])-0.9682458365518543*psr[5]+0.9682458365518543*psl[5]+0.5590169943749476*(psr[1]+psl[1]))*gamma; 
  incr[16] = (((-1.25*bzr[16])+1.25*bzl[16]+0.9682458365518543*(bzr[6]+bzl[6])-0.5590169943749476*bzr[2]+0.5590169943749476*bzl[2])*c+1.25*(psr[16]+psl[16])-0.9682458365518543*psr[6]+0.9682458365518543*psl[6]+0.5590169943749476*(psr[2]+psl[2]))*gamma; 
  incr[17] = (0.25*(bzl[17]*c+psr[17]+psl[17])-0.25*bzr[17]*c)*gamma; 
  incr[18] = (0.25*(bzl[18]*c+psr[18]+psl[18])-0.25*bzr[18]*c)*gamma; 
  incr[19] = (((-1.75*(bzr[19]+bzl[19]))+1.479019945774904*bzr[9]-1.479019945774904*bzl[9]-1.14564392373896*(bzr[3]+bzl[3])+0.6614378277661477*bzr[0]-0.6614378277661477*bzl[0])*c+1.75*psr[19]-1.75*psl[19]-1.479019945774904*(psr[9]+psl[9])+1.14564392373896*psr[3]-1.14564392373896*psl[3]-0.6614378277661477*(psr[0]+psl[0]))*gamma; 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 
  outBzr[4] += incr[4]*dxr1; 
  outBzr[5] += incr[5]*dxr1; 
  outBzr[6] += incr[6]*dxr1; 
  outBzr[7] += incr[7]*dxr1; 
  outBzr[8] += incr[8]*dxr1; 
  outBzr[9] += incr[9]*dxr1; 
  outBzr[10] += incr[10]*dxr1; 
  outBzr[11] += incr[11]*dxr1; 
  outBzr[12] += incr[12]*dxr1; 
  outBzr[13] += incr[13]*dxr1; 
  outBzr[14] += incr[14]*dxr1; 
  outBzr[15] += incr[15]*dxr1; 
  outBzr[16] += incr[16]*dxr1; 
  outBzr[17] += incr[17]*dxr1; 
  outBzr[18] += incr[18]*dxr1; 
  outBzr[19] += incr[19]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += -1.0*incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 
  outBzl[3] += incr[3]*dxl1; 
  outBzl[4] += -1.0*incr[4]*dxl1; 
  outBzl[5] += incr[5]*dxl1; 
  outBzl[6] += incr[6]*dxl1; 
  outBzl[7] += -1.0*incr[7]*dxl1; 
  outBzl[8] += -1.0*incr[8]*dxl1; 
  outBzl[9] += -1.0*incr[9]*dxl1; 
  outBzl[10] += incr[10]*dxl1; 
  outBzl[11] += -1.0*incr[11]*dxl1; 
  outBzl[12] += -1.0*incr[12]*dxl1; 
  outBzl[13] += incr[13]*dxl1; 
  outBzl[14] += incr[14]*dxl1; 
  outBzl[15] += -1.0*incr[15]*dxl1; 
  outBzl[16] += -1.0*incr[16]*dxl1; 
  outBzl[17] += -1.0*incr[17]*dxl1; 
  outBzl[18] += -1.0*incr[18]*dxl1; 
  outBzl[19] += incr[19]*dxl1; 

 
  incr[0] = ((0.6614378277661477*(phr[19]+phl[19])-0.5590169943749475*phr[9]+0.5590169943749475*phl[9]+0.4330127018922193*(phr[3]+phl[3])-0.25*phr[0]+0.25*phl[0])*c-0.6614378277661477*ezr[19]+0.6614378277661477*ezl[19]+0.5590169943749475*(ezr[9]+ezl[9])-0.4330127018922193*ezr[3]+0.4330127018922193*ezl[3]+0.25*(ezr[0]+ezl[0]))*chi; 
  incr[1] = (((-0.5590169943749476*phr[15])+0.5590169943749476*phl[15]+0.4330127018922193*(phr[5]+phl[5])-0.25*phr[1]+0.25*phl[1])*c+0.5590169943749476*(ezr[15]+ezl[15])-0.4330127018922193*ezr[5]+0.4330127018922193*ezl[5]+0.25*(ezr[1]+ezl[1]))*chi; 
  incr[2] = (((-0.5590169943749476*phr[16])+0.5590169943749476*phl[16]+0.4330127018922193*(phr[6]+phl[6])-0.25*phr[2]+0.25*phl[2])*c+0.5590169943749476*(ezr[16]+ezl[16])-0.4330127018922193*ezr[6]+0.4330127018922193*ezl[6]+0.25*(ezr[2]+ezl[2]))*chi; 
  incr[3] = (((-1.14564392373896*(phr[19]+phl[19]))+0.9682458365518543*phr[9]-0.9682458365518543*phl[9]-0.75*(phr[3]+phl[3])+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+1.14564392373896*ezr[19]-1.14564392373896*ezl[19]-0.9682458365518543*(ezr[9]+ezl[9])+0.75*ezr[3]-0.75*ezl[3]-0.4330127018922193*(ezr[0]+ezl[0]))*chi; 
  incr[4] = ((0.4330127018922193*(phr[10]+phl[10])-0.25*phr[4]+0.25*phl[4])*c-0.4330127018922193*ezr[10]+0.4330127018922193*ezl[10]+0.25*(ezr[4]+ezl[4]))*chi; 
  incr[5] = ((0.9682458365518543*phr[15]-0.9682458365518543*phl[15]-0.75*(phr[5]+phl[5])+0.4330127018922193*phr[1]-0.4330127018922193*phl[1])*c-0.9682458365518543*(ezr[15]+ezl[15])+0.75*ezr[5]-0.75*ezl[5]-0.4330127018922193*(ezr[1]+ezl[1]))*chi; 
  incr[6] = ((0.9682458365518543*phr[16]-0.9682458365518543*phl[16]-0.75*(phr[6]+phl[6])+0.4330127018922193*phr[2]-0.4330127018922193*phl[2])*c-0.9682458365518543*(ezr[16]+ezl[16])+0.75*ezr[6]-0.75*ezl[6]-0.4330127018922193*(ezr[2]+ezl[2]))*chi; 
  incr[7] = ((0.4330127018922194*(phr[13]+phl[13])-0.25*phr[7]+0.25*phl[7])*c-0.4330127018922194*ezr[13]+0.4330127018922194*ezl[13]+0.25*(ezr[7]+ezl[7]))*chi; 
  incr[8] = ((0.4330127018922194*(phr[14]+phl[14])-0.25*phr[8]+0.25*phl[8])*c-0.4330127018922194*ezr[14]+0.4330127018922194*ezl[14]+0.25*(ezr[8]+ezl[8]))*chi; 
  incr[9] = ((1.479019945774904*(phr[19]+phl[19])-1.25*phr[9]+1.25*phl[9]+0.9682458365518543*(phr[3]+phl[3])-0.5590169943749475*phr[0]+0.5590169943749475*phl[0])*c-1.479019945774904*ezr[19]+1.479019945774904*ezl[19]+1.25*(ezr[9]+ezl[9])-0.9682458365518543*ezr[3]+0.9682458365518543*ezl[3]+0.5590169943749475*(ezr[0]+ezl[0]))*chi; 
  incr[10] = (((-0.75*(phr[10]+phl[10]))+0.4330127018922193*phr[4]-0.4330127018922193*phl[4])*c+0.75*ezr[10]-0.75*ezl[10]-0.4330127018922193*(ezr[4]+ezl[4]))*chi; 
  incr[11] = (0.25*(phl[11]*c+ezr[11]+ezl[11])-0.25*phr[11]*c)*chi; 
  incr[12] = (0.25*(phl[12]*c+ezr[12]+ezl[12])-0.25*phr[12]*c)*chi; 
  incr[13] = (((-0.75*(phr[13]+phl[13]))+0.4330127018922194*phr[7]-0.4330127018922194*phl[7])*c+0.75*ezr[13]-0.75*ezl[13]-0.4330127018922194*(ezr[7]+ezl[7]))*chi; 
  incr[14] = (((-0.75*(phr[14]+phl[14]))+0.4330127018922194*phr[8]-0.4330127018922194*phl[8])*c+0.75*ezr[14]-0.75*ezl[14]-0.4330127018922194*(ezr[8]+ezl[8]))*chi; 
  incr[15] = (((-1.25*phr[15])+1.25*phl[15]+0.9682458365518543*(phr[5]+phl[5])-0.5590169943749476*phr[1]+0.5590169943749476*phl[1])*c+1.25*(ezr[15]+ezl[15])-0.9682458365518543*ezr[5]+0.9682458365518543*ezl[5]+0.5590169943749476*(ezr[1]+ezl[1]))*chi; 
  incr[16] = (((-1.25*phr[16])+1.25*phl[16]+0.9682458365518543*(phr[6]+phl[6])-0.5590169943749476*phr[2]+0.5590169943749476*phl[2])*c+1.25*(ezr[16]+ezl[16])-0.9682458365518543*ezr[6]+0.9682458365518543*ezl[6]+0.5590169943749476*(ezr[2]+ezl[2]))*chi; 
  incr[17] = (0.25*(phl[17]*c+ezr[17]+ezl[17])-0.25*phr[17]*c)*chi; 
  incr[18] = (0.25*(phl[18]*c+ezr[18]+ezl[18])-0.25*phr[18]*c)*chi; 
  incr[19] = (((-1.75*(phr[19]+phl[19]))+1.479019945774904*phr[9]-1.479019945774904*phl[9]-1.14564392373896*(phr[3]+phl[3])+0.6614378277661477*phr[0]-0.6614378277661477*phl[0])*c+1.75*ezr[19]-1.75*ezl[19]-1.479019945774904*(ezr[9]+ezl[9])+1.14564392373896*ezr[3]-1.14564392373896*ezl[3]-0.6614378277661477*(ezr[0]+ezl[0]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 
  outPhr[4] += incr[4]*dxr1; 
  outPhr[5] += incr[5]*dxr1; 
  outPhr[6] += incr[6]*dxr1; 
  outPhr[7] += incr[7]*dxr1; 
  outPhr[8] += incr[8]*dxr1; 
  outPhr[9] += incr[9]*dxr1; 
  outPhr[10] += incr[10]*dxr1; 
  outPhr[11] += incr[11]*dxr1; 
  outPhr[12] += incr[12]*dxr1; 
  outPhr[13] += incr[13]*dxr1; 
  outPhr[14] += incr[14]*dxr1; 
  outPhr[15] += incr[15]*dxr1; 
  outPhr[16] += incr[16]*dxr1; 
  outPhr[17] += incr[17]*dxr1; 
  outPhr[18] += incr[18]*dxr1; 
  outPhr[19] += incr[19]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += -1.0*incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 
  outPhl[3] += incr[3]*dxl1; 
  outPhl[4] += -1.0*incr[4]*dxl1; 
  outPhl[5] += incr[5]*dxl1; 
  outPhl[6] += incr[6]*dxl1; 
  outPhl[7] += -1.0*incr[7]*dxl1; 
  outPhl[8] += -1.0*incr[8]*dxl1; 
  outPhl[9] += -1.0*incr[9]*dxl1; 
  outPhl[10] += incr[10]*dxl1; 
  outPhl[11] += -1.0*incr[11]*dxl1; 
  outPhl[12] += -1.0*incr[12]*dxl1; 
  outPhl[13] += incr[13]*dxl1; 
  outPhl[14] += incr[14]*dxl1; 
  outPhl[15] += -1.0*incr[15]*dxl1; 
  outPhl[16] += -1.0*incr[16]*dxl1; 
  outPhl[17] += -1.0*incr[17]*dxl1; 
  outPhl[18] += -1.0*incr[18]*dxl1; 
  outPhl[19] += incr[19]*dxl1; 

 
  incr[0] = (0.6614378277661477*(psr[19]+psl[19])-0.5590169943749475*psr[9]+0.5590169943749475*psl[9]+0.4330127018922193*(psr[3]+psl[3])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.6614378277661477*bzr[19])+0.6614378277661477*bzl[19]+0.5590169943749475*(bzr[9]+bzl[9])-0.4330127018922193*bzr[3]+0.4330127018922193*bzl[3]+0.25*(bzr[0]+bzl[0]))*c2gamma; 
  incr[1] = ((-0.5590169943749476*psr[15])+0.5590169943749476*psl[15]+0.4330127018922193*(psr[5]+psl[5])-0.25*psr[1]+0.25*psl[1])*c*gamma+(0.5590169943749476*(bzr[15]+bzl[15])-0.4330127018922193*bzr[5]+0.4330127018922193*bzl[5]+0.25*(bzr[1]+bzl[1]))*c2gamma; 
  incr[2] = ((-0.5590169943749476*psr[16])+0.5590169943749476*psl[16]+0.4330127018922193*(psr[6]+psl[6])-0.25*psr[2]+0.25*psl[2])*c*gamma+(0.5590169943749476*(bzr[16]+bzl[16])-0.4330127018922193*bzr[6]+0.4330127018922193*bzl[6]+0.25*(bzr[2]+bzl[2]))*c2gamma; 
  incr[3] = ((-1.14564392373896*(psr[19]+psl[19]))+0.9682458365518543*psr[9]-0.9682458365518543*psl[9]-0.75*(psr[3]+psl[3])+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(1.14564392373896*bzr[19]-1.14564392373896*bzl[19]-0.9682458365518543*(bzr[9]+bzl[9])+0.75*bzr[3]-0.75*bzl[3]-0.4330127018922193*(bzr[0]+bzl[0]))*c2gamma; 
  incr[4] = (0.4330127018922193*(psr[10]+psl[10])-0.25*psr[4]+0.25*psl[4])*c*gamma+((-0.4330127018922193*bzr[10])+0.4330127018922193*bzl[10]+0.25*(bzr[4]+bzl[4]))*c2gamma; 
  incr[5] = (0.9682458365518543*psr[15]-0.9682458365518543*psl[15]-0.75*(psr[5]+psl[5])+0.4330127018922193*psr[1]-0.4330127018922193*psl[1])*c*gamma+((-0.9682458365518543*(bzr[15]+bzl[15]))+0.75*bzr[5]-0.75*bzl[5]-0.4330127018922193*(bzr[1]+bzl[1]))*c2gamma; 
  incr[6] = (0.9682458365518543*psr[16]-0.9682458365518543*psl[16]-0.75*(psr[6]+psl[6])+0.4330127018922193*psr[2]-0.4330127018922193*psl[2])*c*gamma+((-0.9682458365518543*(bzr[16]+bzl[16]))+0.75*bzr[6]-0.75*bzl[6]-0.4330127018922193*(bzr[2]+bzl[2]))*c2gamma; 
  incr[7] = (0.4330127018922194*(psr[13]+psl[13])-0.25*psr[7]+0.25*psl[7])*c*gamma+((-0.4330127018922194*bzr[13])+0.4330127018922194*bzl[13]+0.25*(bzr[7]+bzl[7]))*c2gamma; 
  incr[8] = (0.4330127018922194*(psr[14]+psl[14])-0.25*psr[8]+0.25*psl[8])*c*gamma+((-0.4330127018922194*bzr[14])+0.4330127018922194*bzl[14]+0.25*(bzr[8]+bzl[8]))*c2gamma; 
  incr[9] = (1.479019945774904*(psr[19]+psl[19])-1.25*psr[9]+1.25*psl[9]+0.9682458365518543*(psr[3]+psl[3])-0.5590169943749475*psr[0]+0.5590169943749475*psl[0])*c*gamma+((-1.479019945774904*bzr[19])+1.479019945774904*bzl[19]+1.25*(bzr[9]+bzl[9])-0.9682458365518543*bzr[3]+0.9682458365518543*bzl[3]+0.5590169943749475*(bzr[0]+bzl[0]))*c2gamma; 
  incr[10] = ((-0.75*(psr[10]+psl[10]))+0.4330127018922193*psr[4]-0.4330127018922193*psl[4])*c*gamma+(0.75*bzr[10]-0.75*bzl[10]-0.4330127018922193*(bzr[4]+bzl[4]))*c2gamma; 
  incr[11] = (0.25*psl[11]-0.25*psr[11])*c*gamma+0.25*(bzr[11]+bzl[11])*c2gamma; 
  incr[12] = (0.25*psl[12]-0.25*psr[12])*c*gamma+0.25*(bzr[12]+bzl[12])*c2gamma; 
  incr[13] = ((-0.75*(psr[13]+psl[13]))+0.4330127018922194*psr[7]-0.4330127018922194*psl[7])*c*gamma+(0.75*bzr[13]-0.75*bzl[13]-0.4330127018922194*(bzr[7]+bzl[7]))*c2gamma; 
  incr[14] = ((-0.75*(psr[14]+psl[14]))+0.4330127018922194*psr[8]-0.4330127018922194*psl[8])*c*gamma+(0.75*bzr[14]-0.75*bzl[14]-0.4330127018922194*(bzr[8]+bzl[8]))*c2gamma; 
  incr[15] = ((-1.25*psr[15])+1.25*psl[15]+0.9682458365518543*(psr[5]+psl[5])-0.5590169943749476*psr[1]+0.5590169943749476*psl[1])*c*gamma+(1.25*(bzr[15]+bzl[15])-0.9682458365518543*bzr[5]+0.9682458365518543*bzl[5]+0.5590169943749476*(bzr[1]+bzl[1]))*c2gamma; 
  incr[16] = ((-1.25*psr[16])+1.25*psl[16]+0.9682458365518543*(psr[6]+psl[6])-0.5590169943749476*psr[2]+0.5590169943749476*psl[2])*c*gamma+(1.25*(bzr[16]+bzl[16])-0.9682458365518543*bzr[6]+0.9682458365518543*bzl[6]+0.5590169943749476*(bzr[2]+bzl[2]))*c2gamma; 
  incr[17] = (0.25*psl[17]-0.25*psr[17])*c*gamma+0.25*(bzr[17]+bzl[17])*c2gamma; 
  incr[18] = (0.25*psl[18]-0.25*psr[18])*c*gamma+0.25*(bzr[18]+bzl[18])*c2gamma; 
  incr[19] = ((-1.75*(psr[19]+psl[19]))+1.479019945774904*psr[9]-1.479019945774904*psl[9]-1.14564392373896*(psr[3]+psl[3])+0.6614378277661477*psr[0]-0.6614378277661477*psl[0])*c*gamma+(1.75*bzr[19]-1.75*bzl[19]-1.479019945774904*(bzr[9]+bzl[9])+1.14564392373896*bzr[3]-1.14564392373896*bzl[3]-0.6614378277661477*(bzr[0]+bzl[0]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 
  outPsr[4] += incr[4]*dxr1; 
  outPsr[5] += incr[5]*dxr1; 
  outPsr[6] += incr[6]*dxr1; 
  outPsr[7] += incr[7]*dxr1; 
  outPsr[8] += incr[8]*dxr1; 
  outPsr[9] += incr[9]*dxr1; 
  outPsr[10] += incr[10]*dxr1; 
  outPsr[11] += incr[11]*dxr1; 
  outPsr[12] += incr[12]*dxr1; 
  outPsr[13] += incr[13]*dxr1; 
  outPsr[14] += incr[14]*dxr1; 
  outPsr[15] += incr[15]*dxr1; 
  outPsr[16] += incr[16]*dxr1; 
  outPsr[17] += incr[17]*dxr1; 
  outPsr[18] += incr[18]*dxr1; 
  outPsr[19] += incr[19]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += -1.0*incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 
  outPsl[3] += incr[3]*dxl1; 
  outPsl[4] += -1.0*incr[4]*dxl1; 
  outPsl[5] += incr[5]*dxl1; 
  outPsl[6] += incr[6]*dxl1; 
  outPsl[7] += -1.0*incr[7]*dxl1; 
  outPsl[8] += -1.0*incr[8]*dxl1; 
  outPsl[9] += -1.0*incr[9]*dxl1; 
  outPsl[10] += incr[10]*dxl1; 
  outPsl[11] += -1.0*incr[11]*dxl1; 
  outPsl[12] += -1.0*incr[12]*dxl1; 
  outPsl[13] += incr[13]*dxl1; 
  outPsl[14] += incr[14]*dxl1; 
  outPsl[15] += -1.0*incr[15]*dxl1; 
  outPsl[16] += -1.0*incr[16]*dxl1; 
  outPsl[17] += -1.0*incr[17]*dxl1; 
  outPsl[18] += -1.0*incr[18]*dxl1; 
  outPsl[19] += incr[19]*dxl1; 

 
  return c; 
} 
