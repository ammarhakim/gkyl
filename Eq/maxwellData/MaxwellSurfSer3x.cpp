#include <MaxwellModDecl.h> 
double MaxwellSurf3xSer_X_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[0]; 
  const double dxr1 = 2.0/dxr[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[8]; 
  const double *ezl = &ql[16]; 
  const double *bxl = &ql[24]; 
  const double *byl = &ql[32]; 
  const double *bzl = &ql[40]; 
  const double *phl = &ql[48]; 
  const double *psl = &ql[56]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[8]; 
  double *outEzl = &outl[16]; 
  double *outBxl = &outl[24]; 
  double *outByl = &outl[32]; 
  double *outBzl = &outl[40]; 
  double *outPhl = &outl[48]; 
  double *outPsl = &outl[56]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[8]; 
  const double *ezr = &qr[16]; 
  const double *bxr = &qr[24]; 
  const double *byr = &qr[32]; 
  const double *bzr = &qr[40]; 
  const double *phr = &qr[48]; 
  const double *psr = &qr[56]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[8]; 
  double *outEzr = &outr[16]; 
  double *outBxr = &outr[24]; 
  double *outByr = &outr[32]; 
  double *outBzr = &outr[40]; 
  double *outPhr = &outr[48]; 
  double *outPsr = &outr[56]; 
 
  double incr[8]; 
 
  incr[0] = (0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+((-0.4330127018922193*phr[1])+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = ((-0.75*(exr[1]+exl[1]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+(0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[2] = (0.4330127018922193*(exr[4]+exl[4])-0.25*exr[2]+0.25*exl[2])*c*chi+((-0.4330127018922193*phr[4])+0.4330127018922193*phl[4]+0.25*(phr[2]+phl[2]))*c2chi; 
  incr[3] = (0.4330127018922193*(exr[5]+exl[5])-0.25*exr[3]+0.25*exl[3])*c*chi+((-0.4330127018922193*phr[5])+0.4330127018922193*phl[5]+0.25*(phr[3]+phl[3]))*c2chi; 
  incr[4] = ((-0.75*(exr[4]+exl[4]))+0.4330127018922193*exr[2]-0.4330127018922193*exl[2])*c*chi+(0.75*phr[4]-0.75*phl[4]-0.4330127018922193*(phr[2]+phl[2]))*c2chi; 
  incr[5] = ((-0.75*(exr[5]+exl[5]))+0.4330127018922193*exr[3]-0.4330127018922193*exl[3])*c*chi+(0.75*phr[5]-0.75*phl[5]-0.4330127018922193*(phr[3]+phl[3]))*c2chi; 
  incr[6] = (0.4330127018922193*(exr[7]+exl[7])-0.25*exr[6]+0.25*exl[6])*c*chi+((-0.4330127018922193*phr[7])+0.4330127018922193*phl[7]+0.25*(phr[6]+phl[6]))*c2chi; 
  incr[7] = ((-0.75*(exr[7]+exl[7]))+0.4330127018922193*exr[6]-0.4330127018922193*exl[6])*c*chi+(0.75*phr[7]-0.75*phl[7]-0.4330127018922193*(phr[6]+phl[6]))*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 
  outExr[4] += incr[4]*dxr1; 
  outExr[5] += incr[5]*dxr1; 
  outExr[6] += incr[6]*dxr1; 
  outExr[7] += incr[7]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 
  outExl[3] += -1.0*incr[3]*dxl1; 
  outExl[4] += incr[4]*dxl1; 
  outExl[5] += incr[5]*dxl1; 
  outExl[6] += -1.0*incr[6]*dxl1; 
  outExl[7] += incr[7]*dxl1; 

 
  incr[0] = (0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+((-0.4330127018922193*bzr[1])+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = ((-0.75*(eyr[1]+eyl[1]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+(0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[2] = (0.4330127018922193*(eyr[4]+eyl[4])-0.25*eyr[2]+0.25*eyl[2])*tau+((-0.4330127018922193*bzr[4])+0.4330127018922193*bzl[4]+0.25*(bzr[2]+bzl[2]))*c2; 
  incr[3] = (0.4330127018922193*(eyr[5]+eyl[5])-0.25*eyr[3]+0.25*eyl[3])*tau+((-0.4330127018922193*bzr[5])+0.4330127018922193*bzl[5]+0.25*(bzr[3]+bzl[3]))*c2; 
  incr[4] = ((-0.75*(eyr[4]+eyl[4]))+0.4330127018922193*eyr[2]-0.4330127018922193*eyl[2])*tau+(0.75*bzr[4]-0.75*bzl[4]-0.4330127018922193*(bzr[2]+bzl[2]))*c2; 
  incr[5] = ((-0.75*(eyr[5]+eyl[5]))+0.4330127018922193*eyr[3]-0.4330127018922193*eyl[3])*tau+(0.75*bzr[5]-0.75*bzl[5]-0.4330127018922193*(bzr[3]+bzl[3]))*c2; 
  incr[6] = (0.4330127018922193*(eyr[7]+eyl[7])-0.25*eyr[6]+0.25*eyl[6])*tau+((-0.4330127018922193*bzr[7])+0.4330127018922193*bzl[7]+0.25*(bzr[6]+bzl[6]))*c2; 
  incr[7] = ((-0.75*(eyr[7]+eyl[7]))+0.4330127018922193*eyr[6]-0.4330127018922193*eyl[6])*tau+(0.75*bzr[7]-0.75*bzl[7]-0.4330127018922193*(bzr[6]+bzl[6]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 
  outEyr[4] += incr[4]*dxr1; 
  outEyr[5] += incr[5]*dxr1; 
  outEyr[6] += incr[6]*dxr1; 
  outEyr[7] += incr[7]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 
  outEyl[3] += -1.0*incr[3]*dxl1; 
  outEyl[4] += incr[4]*dxl1; 
  outEyl[5] += incr[5]*dxl1; 
  outEyl[6] += -1.0*incr[6]*dxl1; 
  outEyl[7] += incr[7]*dxl1; 

 
  incr[0] = (0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+(0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = ((-0.75*(ezr[1]+ezl[1]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+((-0.75*byr[1])+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[2] = (0.4330127018922193*(ezr[4]+ezl[4])-0.25*ezr[2]+0.25*ezl[2])*tau+(0.4330127018922193*byr[4]-0.4330127018922193*byl[4]-0.25*(byr[2]+byl[2]))*c2; 
  incr[3] = (0.4330127018922193*(ezr[5]+ezl[5])-0.25*ezr[3]+0.25*ezl[3])*tau+(0.4330127018922193*byr[5]-0.4330127018922193*byl[5]-0.25*(byr[3]+byl[3]))*c2; 
  incr[4] = ((-0.75*(ezr[4]+ezl[4]))+0.4330127018922193*ezr[2]-0.4330127018922193*ezl[2])*tau+((-0.75*byr[4])+0.75*byl[4]+0.4330127018922193*(byr[2]+byl[2]))*c2; 
  incr[5] = ((-0.75*(ezr[5]+ezl[5]))+0.4330127018922193*ezr[3]-0.4330127018922193*ezl[3])*tau+((-0.75*byr[5])+0.75*byl[5]+0.4330127018922193*(byr[3]+byl[3]))*c2; 
  incr[6] = (0.4330127018922193*(ezr[7]+ezl[7])-0.25*ezr[6]+0.25*ezl[6])*tau+(0.4330127018922193*byr[7]-0.4330127018922193*byl[7]-0.25*(byr[6]+byl[6]))*c2; 
  incr[7] = ((-0.75*(ezr[7]+ezl[7]))+0.4330127018922193*ezr[6]-0.4330127018922193*ezl[6])*tau+((-0.75*byr[7])+0.75*byl[7]+0.4330127018922193*(byr[6]+byl[6]))*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 
  outEzr[4] += incr[4]*dxr1; 
  outEzr[5] += incr[5]*dxr1; 
  outEzr[6] += incr[6]*dxr1; 
  outEzr[7] += incr[7]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 
  outEzl[3] += -1.0*incr[3]*dxl1; 
  outEzl[4] += incr[4]*dxl1; 
  outEzl[5] += incr[5]*dxl1; 
  outEzl[6] += -1.0*incr[6]*dxl1; 
  outEzl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (((-0.75*(bxr[1]+bxl[1]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[2] = ((0.4330127018922193*(bxr[4]+bxl[4])-0.25*bxr[2]+0.25*bxl[2])*c-0.4330127018922193*psr[4]+0.4330127018922193*psl[4]+0.25*(psr[2]+psl[2]))*gamma; 
  incr[3] = ((0.4330127018922193*(bxr[5]+bxl[5])-0.25*bxr[3]+0.25*bxl[3])*c-0.4330127018922193*psr[5]+0.4330127018922193*psl[5]+0.25*(psr[3]+psl[3]))*gamma; 
  incr[4] = (((-0.75*(bxr[4]+bxl[4]))+0.4330127018922193*bxr[2]-0.4330127018922193*bxl[2])*c+0.75*psr[4]-0.75*psl[4]-0.4330127018922193*(psr[2]+psl[2]))*gamma; 
  incr[5] = (((-0.75*(bxr[5]+bxl[5]))+0.4330127018922193*bxr[3]-0.4330127018922193*bxl[3])*c+0.75*psr[5]-0.75*psl[5]-0.4330127018922193*(psr[3]+psl[3]))*gamma; 
  incr[6] = ((0.4330127018922193*(bxr[7]+bxl[7])-0.25*bxr[6]+0.25*bxl[6])*c-0.4330127018922193*psr[7]+0.4330127018922193*psl[7]+0.25*(psr[6]+psl[6]))*gamma; 
  incr[7] = (((-0.75*(bxr[7]+bxl[7]))+0.4330127018922193*bxr[6]-0.4330127018922193*bxl[6])*c+0.75*psr[7]-0.75*psl[7]-0.4330127018922193*(psr[6]+psl[6]))*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 
  outBxr[4] += incr[4]*dxr1; 
  outBxr[5] += incr[5]*dxr1; 
  outBxr[6] += incr[6]*dxr1; 
  outBxr[7] += incr[7]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 
  outBxl[3] += -1.0*incr[3]*dxl1; 
  outBxl[4] += incr[4]*dxl1; 
  outBxl[5] += incr[5]*dxl1; 
  outBxl[6] += -1.0*incr[6]*dxl1; 
  outBxl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = (((-0.75*(byr[1]+byl[1]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[2] = ((0.4330127018922193*(byr[4]+byl[4])-0.25*byr[2]+0.25*byl[2])*c2)/tau+0.4330127018922193*ezr[4]-0.4330127018922193*ezl[4]-0.25*(ezr[2]+ezl[2]); 
  incr[3] = ((0.4330127018922193*(byr[5]+byl[5])-0.25*byr[3]+0.25*byl[3])*c2)/tau+0.4330127018922193*ezr[5]-0.4330127018922193*ezl[5]-0.25*(ezr[3]+ezl[3]); 
  incr[4] = (((-0.75*(byr[4]+byl[4]))+0.4330127018922193*byr[2]-0.4330127018922193*byl[2])*c2)/tau-0.75*ezr[4]+0.75*ezl[4]+0.4330127018922193*(ezr[2]+ezl[2]); 
  incr[5] = (((-0.75*(byr[5]+byl[5]))+0.4330127018922193*byr[3]-0.4330127018922193*byl[3])*c2)/tau-0.75*ezr[5]+0.75*ezl[5]+0.4330127018922193*(ezr[3]+ezl[3]); 
  incr[6] = ((0.4330127018922193*(byr[7]+byl[7])-0.25*byr[6]+0.25*byl[6])*c2)/tau+0.4330127018922193*ezr[7]-0.4330127018922193*ezl[7]-0.25*(ezr[6]+ezl[6]); 
  incr[7] = (((-0.75*(byr[7]+byl[7]))+0.4330127018922193*byr[6]-0.4330127018922193*byl[6])*c2)/tau-0.75*ezr[7]+0.75*ezl[7]+0.4330127018922193*(ezr[6]+ezl[6]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 
  outByr[4] += incr[4]*dxr1; 
  outByr[5] += incr[5]*dxr1; 
  outByr[6] += incr[6]*dxr1; 
  outByr[7] += incr[7]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 
  outByl[3] += -1.0*incr[3]*dxl1; 
  outByl[4] += incr[4]*dxl1; 
  outByl[5] += incr[5]*dxl1; 
  outByl[6] += -1.0*incr[6]*dxl1; 
  outByl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = (((-0.75*(bzr[1]+bzl[1]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[2] = ((0.4330127018922193*(bzr[4]+bzl[4])-0.25*bzr[2]+0.25*bzl[2])*c2)/tau-0.4330127018922193*eyr[4]+0.4330127018922193*eyl[4]+0.25*(eyr[2]+eyl[2]); 
  incr[3] = ((0.4330127018922193*(bzr[5]+bzl[5])-0.25*bzr[3]+0.25*bzl[3])*c2)/tau-0.4330127018922193*eyr[5]+0.4330127018922193*eyl[5]+0.25*(eyr[3]+eyl[3]); 
  incr[4] = (((-0.75*(bzr[4]+bzl[4]))+0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2])*c2)/tau+0.75*eyr[4]-0.75*eyl[4]-0.4330127018922193*(eyr[2]+eyl[2]); 
  incr[5] = (((-0.75*(bzr[5]+bzl[5]))+0.4330127018922193*bzr[3]-0.4330127018922193*bzl[3])*c2)/tau+0.75*eyr[5]-0.75*eyl[5]-0.4330127018922193*(eyr[3]+eyl[3]); 
  incr[6] = ((0.4330127018922193*(bzr[7]+bzl[7])-0.25*bzr[6]+0.25*bzl[6])*c2)/tau-0.4330127018922193*eyr[7]+0.4330127018922193*eyl[7]+0.25*(eyr[6]+eyl[6]); 
  incr[7] = (((-0.75*(bzr[7]+bzl[7]))+0.4330127018922193*bzr[6]-0.4330127018922193*bzl[6])*c2)/tau+0.75*eyr[7]-0.75*eyl[7]-0.4330127018922193*(eyr[6]+eyl[6]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 
  outBzr[4] += incr[4]*dxr1; 
  outBzr[5] += incr[5]*dxr1; 
  outBzr[6] += incr[6]*dxr1; 
  outBzr[7] += incr[7]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 
  outBzl[3] += -1.0*incr[3]*dxl1; 
  outBzl[4] += incr[4]*dxl1; 
  outBzl[5] += incr[5]*dxl1; 
  outBzl[6] += -1.0*incr[6]*dxl1; 
  outBzl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = (((-0.75*(phr[1]+phl[1]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 
  incr[2] = ((0.4330127018922193*(phr[4]+phl[4])-0.25*phr[2]+0.25*phl[2])*c-0.4330127018922193*exr[4]+0.4330127018922193*exl[4]+0.25*(exr[2]+exl[2]))*chi; 
  incr[3] = ((0.4330127018922193*(phr[5]+phl[5])-0.25*phr[3]+0.25*phl[3])*c-0.4330127018922193*exr[5]+0.4330127018922193*exl[5]+0.25*(exr[3]+exl[3]))*chi; 
  incr[4] = (((-0.75*(phr[4]+phl[4]))+0.4330127018922193*phr[2]-0.4330127018922193*phl[2])*c+0.75*exr[4]-0.75*exl[4]-0.4330127018922193*(exr[2]+exl[2]))*chi; 
  incr[5] = (((-0.75*(phr[5]+phl[5]))+0.4330127018922193*phr[3]-0.4330127018922193*phl[3])*c+0.75*exr[5]-0.75*exl[5]-0.4330127018922193*(exr[3]+exl[3]))*chi; 
  incr[6] = ((0.4330127018922193*(phr[7]+phl[7])-0.25*phr[6]+0.25*phl[6])*c-0.4330127018922193*exr[7]+0.4330127018922193*exl[7]+0.25*(exr[6]+exl[6]))*chi; 
  incr[7] = (((-0.75*(phr[7]+phl[7]))+0.4330127018922193*phr[6]-0.4330127018922193*phl[6])*c+0.75*exr[7]-0.75*exl[7]-0.4330127018922193*(exr[6]+exl[6]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 
  outPhr[4] += incr[4]*dxr1; 
  outPhr[5] += incr[5]*dxr1; 
  outPhr[6] += incr[6]*dxr1; 
  outPhr[7] += incr[7]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 
  outPhl[3] += -1.0*incr[3]*dxl1; 
  outPhl[4] += incr[4]*dxl1; 
  outPhl[5] += incr[5]*dxl1; 
  outPhl[6] += -1.0*incr[6]*dxl1; 
  outPhl[7] += incr[7]*dxl1; 

 
  incr[0] = (0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*bxr[1])+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = ((-0.75*(psr[1]+psl[1]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 
  incr[2] = (0.4330127018922193*(psr[4]+psl[4])-0.25*psr[2]+0.25*psl[2])*c*gamma+((-0.4330127018922193*bxr[4])+0.4330127018922193*bxl[4]+0.25*(bxr[2]+bxl[2]))*c2gamma; 
  incr[3] = (0.4330127018922193*(psr[5]+psl[5])-0.25*psr[3]+0.25*psl[3])*c*gamma+((-0.4330127018922193*bxr[5])+0.4330127018922193*bxl[5]+0.25*(bxr[3]+bxl[3]))*c2gamma; 
  incr[4] = ((-0.75*(psr[4]+psl[4]))+0.4330127018922193*psr[2]-0.4330127018922193*psl[2])*c*gamma+(0.75*bxr[4]-0.75*bxl[4]-0.4330127018922193*(bxr[2]+bxl[2]))*c2gamma; 
  incr[5] = ((-0.75*(psr[5]+psl[5]))+0.4330127018922193*psr[3]-0.4330127018922193*psl[3])*c*gamma+(0.75*bxr[5]-0.75*bxl[5]-0.4330127018922193*(bxr[3]+bxl[3]))*c2gamma; 
  incr[6] = (0.4330127018922193*(psr[7]+psl[7])-0.25*psr[6]+0.25*psl[6])*c*gamma+((-0.4330127018922193*bxr[7])+0.4330127018922193*bxl[7]+0.25*(bxr[6]+bxl[6]))*c2gamma; 
  incr[7] = ((-0.75*(psr[7]+psl[7]))+0.4330127018922193*psr[6]-0.4330127018922193*psl[6])*c*gamma+(0.75*bxr[7]-0.75*bxl[7]-0.4330127018922193*(bxr[6]+bxl[6]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 
  outPsr[4] += incr[4]*dxr1; 
  outPsr[5] += incr[5]*dxr1; 
  outPsr[6] += incr[6]*dxr1; 
  outPsr[7] += incr[7]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 
  outPsl[3] += -1.0*incr[3]*dxl1; 
  outPsl[4] += incr[4]*dxl1; 
  outPsl[5] += incr[5]*dxl1; 
  outPsl[6] += -1.0*incr[6]*dxl1; 
  outPsl[7] += incr[7]*dxl1; 

 
  return std::fmax(c, tau); 
} 
double MaxwellSurf3xSer_Y_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[1]; 
  const double dxr1 = 2.0/dxr[1]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[8]; 
  const double *ezl = &ql[16]; 
  const double *bxl = &ql[24]; 
  const double *byl = &ql[32]; 
  const double *bzl = &ql[40]; 
  const double *phl = &ql[48]; 
  const double *psl = &ql[56]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[8]; 
  double *outEzl = &outl[16]; 
  double *outBxl = &outl[24]; 
  double *outByl = &outl[32]; 
  double *outBzl = &outl[40]; 
  double *outPhl = &outl[48]; 
  double *outPsl = &outl[56]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[8]; 
  const double *ezr = &qr[16]; 
  const double *bxr = &qr[24]; 
  const double *byr = &qr[32]; 
  const double *bzr = &qr[40]; 
  const double *phr = &qr[48]; 
  const double *psr = &qr[56]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[8]; 
  double *outEzr = &outr[16]; 
  double *outBxr = &outr[24]; 
  double *outByr = &outr[32]; 
  double *outBzr = &outr[40]; 
  double *outPhr = &outr[48]; 
  double *outPsr = &outr[56]; 
 
  double incr[8]; 
 
  incr[0] = (0.4330127018922193*(exr[2]+exl[2])-0.25*exr[0]+0.25*exl[0])*tau+(0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2]-0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = (0.4330127018922193*(exr[4]+exl[4])-0.25*exr[1]+0.25*exl[1])*tau+(0.4330127018922193*bzr[4]-0.4330127018922193*bzl[4]-0.25*(bzr[1]+bzl[1]))*c2; 
  incr[2] = ((-0.75*(exr[2]+exl[2]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*tau+((-0.75*bzr[2])+0.75*bzl[2]+0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[3] = (0.4330127018922193*(exr[6]+exl[6])-0.25*exr[3]+0.25*exl[3])*tau+(0.4330127018922193*bzr[6]-0.4330127018922193*bzl[6]-0.25*(bzr[3]+bzl[3]))*c2; 
  incr[4] = ((-0.75*(exr[4]+exl[4]))+0.4330127018922193*exr[1]-0.4330127018922193*exl[1])*tau+((-0.75*bzr[4])+0.75*bzl[4]+0.4330127018922193*(bzr[1]+bzl[1]))*c2; 
  incr[5] = (0.4330127018922193*(exr[7]+exl[7])-0.25*exr[5]+0.25*exl[5])*tau+(0.4330127018922193*bzr[7]-0.4330127018922193*bzl[7]-0.25*(bzr[5]+bzl[5]))*c2; 
  incr[6] = ((-0.75*(exr[6]+exl[6]))+0.4330127018922193*exr[3]-0.4330127018922193*exl[3])*tau+((-0.75*bzr[6])+0.75*bzl[6]+0.4330127018922193*(bzr[3]+bzl[3]))*c2; 
  incr[7] = ((-0.75*(exr[7]+exl[7]))+0.4330127018922193*exr[5]-0.4330127018922193*exl[5])*tau+((-0.75*bzr[7])+0.75*bzl[7]+0.4330127018922193*(bzr[5]+bzl[5]))*c2; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 
  outExr[4] += incr[4]*dxr1; 
  outExr[5] += incr[5]*dxr1; 
  outExr[6] += incr[6]*dxr1; 
  outExr[7] += incr[7]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += -1.0*incr[1]*dxl1; 
  outExl[2] += incr[2]*dxl1; 
  outExl[3] += -1.0*incr[3]*dxl1; 
  outExl[4] += incr[4]*dxl1; 
  outExl[5] += -1.0*incr[5]*dxl1; 
  outExl[6] += incr[6]*dxl1; 
  outExl[7] += incr[7]*dxl1; 

 
  incr[0] = (0.4330127018922193*(eyr[2]+eyl[2])-0.25*eyr[0]+0.25*eyl[0])*c*chi+((-0.4330127018922193*phr[2])+0.4330127018922193*phl[2]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.4330127018922193*(eyr[4]+eyl[4])-0.25*eyr[1]+0.25*eyl[1])*c*chi+((-0.4330127018922193*phr[4])+0.4330127018922193*phl[4]+0.25*(phr[1]+phl[1]))*c2chi; 
  incr[2] = ((-0.75*(eyr[2]+eyl[2]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*c*chi+(0.75*phr[2]-0.75*phl[2]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[3] = (0.4330127018922193*(eyr[6]+eyl[6])-0.25*eyr[3]+0.25*eyl[3])*c*chi+((-0.4330127018922193*phr[6])+0.4330127018922193*phl[6]+0.25*(phr[3]+phl[3]))*c2chi; 
  incr[4] = ((-0.75*(eyr[4]+eyl[4]))+0.4330127018922193*eyr[1]-0.4330127018922193*eyl[1])*c*chi+(0.75*phr[4]-0.75*phl[4]-0.4330127018922193*(phr[1]+phl[1]))*c2chi; 
  incr[5] = (0.4330127018922193*(eyr[7]+eyl[7])-0.25*eyr[5]+0.25*eyl[5])*c*chi+((-0.4330127018922193*phr[7])+0.4330127018922193*phl[7]+0.25*(phr[5]+phl[5]))*c2chi; 
  incr[6] = ((-0.75*(eyr[6]+eyl[6]))+0.4330127018922193*eyr[3]-0.4330127018922193*eyl[3])*c*chi+(0.75*phr[6]-0.75*phl[6]-0.4330127018922193*(phr[3]+phl[3]))*c2chi; 
  incr[7] = ((-0.75*(eyr[7]+eyl[7]))+0.4330127018922193*eyr[5]-0.4330127018922193*eyl[5])*c*chi+(0.75*phr[7]-0.75*phl[7]-0.4330127018922193*(phr[5]+phl[5]))*c2chi; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 
  outEyr[4] += incr[4]*dxr1; 
  outEyr[5] += incr[5]*dxr1; 
  outEyr[6] += incr[6]*dxr1; 
  outEyr[7] += incr[7]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += -1.0*incr[1]*dxl1; 
  outEyl[2] += incr[2]*dxl1; 
  outEyl[3] += -1.0*incr[3]*dxl1; 
  outEyl[4] += incr[4]*dxl1; 
  outEyl[5] += -1.0*incr[5]*dxl1; 
  outEyl[6] += incr[6]*dxl1; 
  outEyl[7] += incr[7]*dxl1; 

 
  incr[0] = (0.4330127018922193*(ezr[2]+ezl[2])-0.25*ezr[0]+0.25*ezl[0])*tau+((-0.4330127018922193*bxr[2])+0.4330127018922193*bxl[2]+0.25*(bxr[0]+bxl[0]))*c2; 
  incr[1] = (0.4330127018922193*(ezr[4]+ezl[4])-0.25*ezr[1]+0.25*ezl[1])*tau+((-0.4330127018922193*bxr[4])+0.4330127018922193*bxl[4]+0.25*(bxr[1]+bxl[1]))*c2; 
  incr[2] = ((-0.75*(ezr[2]+ezl[2]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+(0.75*bxr[2]-0.75*bxl[2]-0.4330127018922193*(bxr[0]+bxl[0]))*c2; 
  incr[3] = (0.4330127018922193*(ezr[6]+ezl[6])-0.25*ezr[3]+0.25*ezl[3])*tau+((-0.4330127018922193*bxr[6])+0.4330127018922193*bxl[6]+0.25*(bxr[3]+bxl[3]))*c2; 
  incr[4] = ((-0.75*(ezr[4]+ezl[4]))+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1])*tau+(0.75*bxr[4]-0.75*bxl[4]-0.4330127018922193*(bxr[1]+bxl[1]))*c2; 
  incr[5] = (0.4330127018922193*(ezr[7]+ezl[7])-0.25*ezr[5]+0.25*ezl[5])*tau+((-0.4330127018922193*bxr[7])+0.4330127018922193*bxl[7]+0.25*(bxr[5]+bxl[5]))*c2; 
  incr[6] = ((-0.75*(ezr[6]+ezl[6]))+0.4330127018922193*ezr[3]-0.4330127018922193*ezl[3])*tau+(0.75*bxr[6]-0.75*bxl[6]-0.4330127018922193*(bxr[3]+bxl[3]))*c2; 
  incr[7] = ((-0.75*(ezr[7]+ezl[7]))+0.4330127018922193*ezr[5]-0.4330127018922193*ezl[5])*tau+(0.75*bxr[7]-0.75*bxl[7]-0.4330127018922193*(bxr[5]+bxl[5]))*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 
  outEzr[4] += incr[4]*dxr1; 
  outEzr[5] += incr[5]*dxr1; 
  outEzr[6] += incr[6]*dxr1; 
  outEzr[7] += incr[7]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += -1.0*incr[1]*dxl1; 
  outEzl[2] += incr[2]*dxl1; 
  outEzl[3] += -1.0*incr[3]*dxl1; 
  outEzl[4] += incr[4]*dxl1; 
  outEzl[5] += -1.0*incr[5]*dxl1; 
  outEzl[6] += incr[6]*dxl1; 
  outEzl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bxr[2]+bxl[2])-0.25*bxr[0]+0.25*bxl[0])*c2)/tau-0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]+0.25*(ezr[0]+ezl[0]); 
  incr[1] = ((0.4330127018922193*(bxr[4]+bxl[4])-0.25*bxr[1]+0.25*bxl[1])*c2)/tau-0.4330127018922193*ezr[4]+0.4330127018922193*ezl[4]+0.25*(ezr[1]+ezl[1]); 
  incr[2] = (((-0.75*(bxr[2]+bxl[2]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c2)/tau+0.75*ezr[2]-0.75*ezl[2]-0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[3] = ((0.4330127018922193*(bxr[6]+bxl[6])-0.25*bxr[3]+0.25*bxl[3])*c2)/tau-0.4330127018922193*ezr[6]+0.4330127018922193*ezl[6]+0.25*(ezr[3]+ezl[3]); 
  incr[4] = (((-0.75*(bxr[4]+bxl[4]))+0.4330127018922193*bxr[1]-0.4330127018922193*bxl[1])*c2)/tau+0.75*ezr[4]-0.75*ezl[4]-0.4330127018922193*(ezr[1]+ezl[1]); 
  incr[5] = ((0.4330127018922193*(bxr[7]+bxl[7])-0.25*bxr[5]+0.25*bxl[5])*c2)/tau-0.4330127018922193*ezr[7]+0.4330127018922193*ezl[7]+0.25*(ezr[5]+ezl[5]); 
  incr[6] = (((-0.75*(bxr[6]+bxl[6]))+0.4330127018922193*bxr[3]-0.4330127018922193*bxl[3])*c2)/tau+0.75*ezr[6]-0.75*ezl[6]-0.4330127018922193*(ezr[3]+ezl[3]); 
  incr[7] = (((-0.75*(bxr[7]+bxl[7]))+0.4330127018922193*bxr[5]-0.4330127018922193*bxl[5])*c2)/tau+0.75*ezr[7]-0.75*ezl[7]-0.4330127018922193*(ezr[5]+ezl[5]); 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 
  outBxr[4] += incr[4]*dxr1; 
  outBxr[5] += incr[5]*dxr1; 
  outBxr[6] += incr[6]*dxr1; 
  outBxr[7] += incr[7]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += -1.0*incr[1]*dxl1; 
  outBxl[2] += incr[2]*dxl1; 
  outBxl[3] += -1.0*incr[3]*dxl1; 
  outBxl[4] += incr[4]*dxl1; 
  outBxl[5] += -1.0*incr[5]*dxl1; 
  outBxl[6] += incr[6]*dxl1; 
  outBxl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(byr[2]+byl[2])-0.25*byr[0]+0.25*byl[0])*c-0.4330127018922193*psr[2]+0.4330127018922193*psl[2]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = ((0.4330127018922193*(byr[4]+byl[4])-0.25*byr[1]+0.25*byl[1])*c-0.4330127018922193*psr[4]+0.4330127018922193*psl[4]+0.25*(psr[1]+psl[1]))*gamma; 
  incr[2] = (((-0.75*(byr[2]+byl[2]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c+0.75*psr[2]-0.75*psl[2]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[3] = ((0.4330127018922193*(byr[6]+byl[6])-0.25*byr[3]+0.25*byl[3])*c-0.4330127018922193*psr[6]+0.4330127018922193*psl[6]+0.25*(psr[3]+psl[3]))*gamma; 
  incr[4] = (((-0.75*(byr[4]+byl[4]))+0.4330127018922193*byr[1]-0.4330127018922193*byl[1])*c+0.75*psr[4]-0.75*psl[4]-0.4330127018922193*(psr[1]+psl[1]))*gamma; 
  incr[5] = ((0.4330127018922193*(byr[7]+byl[7])-0.25*byr[5]+0.25*byl[5])*c-0.4330127018922193*psr[7]+0.4330127018922193*psl[7]+0.25*(psr[5]+psl[5]))*gamma; 
  incr[6] = (((-0.75*(byr[6]+byl[6]))+0.4330127018922193*byr[3]-0.4330127018922193*byl[3])*c+0.75*psr[6]-0.75*psl[6]-0.4330127018922193*(psr[3]+psl[3]))*gamma; 
  incr[7] = (((-0.75*(byr[7]+byl[7]))+0.4330127018922193*byr[5]-0.4330127018922193*byl[5])*c+0.75*psr[7]-0.75*psl[7]-0.4330127018922193*(psr[5]+psl[5]))*gamma; 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 
  outByr[4] += incr[4]*dxr1; 
  outByr[5] += incr[5]*dxr1; 
  outByr[6] += incr[6]*dxr1; 
  outByr[7] += incr[7]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += -1.0*incr[1]*dxl1; 
  outByl[2] += incr[2]*dxl1; 
  outByl[3] += -1.0*incr[3]*dxl1; 
  outByl[4] += incr[4]*dxl1; 
  outByl[5] += -1.0*incr[5]*dxl1; 
  outByl[6] += incr[6]*dxl1; 
  outByl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bzr[2]+bzl[2])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau+0.4330127018922193*exr[2]-0.4330127018922193*exl[2]-0.25*(exr[0]+exl[0]); 
  incr[1] = ((0.4330127018922193*(bzr[4]+bzl[4])-0.25*bzr[1]+0.25*bzl[1])*c2)/tau+0.4330127018922193*exr[4]-0.4330127018922193*exl[4]-0.25*(exr[1]+exl[1]); 
  incr[2] = (((-0.75*(bzr[2]+bzl[2]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau-0.75*exr[2]+0.75*exl[2]+0.4330127018922193*(exr[0]+exl[0]); 
  incr[3] = ((0.4330127018922193*(bzr[6]+bzl[6])-0.25*bzr[3]+0.25*bzl[3])*c2)/tau+0.4330127018922193*exr[6]-0.4330127018922193*exl[6]-0.25*(exr[3]+exl[3]); 
  incr[4] = (((-0.75*(bzr[4]+bzl[4]))+0.4330127018922193*bzr[1]-0.4330127018922193*bzl[1])*c2)/tau-0.75*exr[4]+0.75*exl[4]+0.4330127018922193*(exr[1]+exl[1]); 
  incr[5] = ((0.4330127018922193*(bzr[7]+bzl[7])-0.25*bzr[5]+0.25*bzl[5])*c2)/tau+0.4330127018922193*exr[7]-0.4330127018922193*exl[7]-0.25*(exr[5]+exl[5]); 
  incr[6] = (((-0.75*(bzr[6]+bzl[6]))+0.4330127018922193*bzr[3]-0.4330127018922193*bzl[3])*c2)/tau-0.75*exr[6]+0.75*exl[6]+0.4330127018922193*(exr[3]+exl[3]); 
  incr[7] = (((-0.75*(bzr[7]+bzl[7]))+0.4330127018922193*bzr[5]-0.4330127018922193*bzl[5])*c2)/tau-0.75*exr[7]+0.75*exl[7]+0.4330127018922193*(exr[5]+exl[5]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 
  outBzr[4] += incr[4]*dxr1; 
  outBzr[5] += incr[5]*dxr1; 
  outBzr[6] += incr[6]*dxr1; 
  outBzr[7] += incr[7]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += -1.0*incr[1]*dxl1; 
  outBzl[2] += incr[2]*dxl1; 
  outBzl[3] += -1.0*incr[3]*dxl1; 
  outBzl[4] += incr[4]*dxl1; 
  outBzl[5] += -1.0*incr[5]*dxl1; 
  outBzl[6] += incr[6]*dxl1; 
  outBzl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(phr[2]+phl[2])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*eyr[2]+0.4330127018922193*eyl[2]+0.25*(eyr[0]+eyl[0]))*chi; 
  incr[1] = ((0.4330127018922193*(phr[4]+phl[4])-0.25*phr[1]+0.25*phl[1])*c-0.4330127018922193*eyr[4]+0.4330127018922193*eyl[4]+0.25*(eyr[1]+eyl[1]))*chi; 
  incr[2] = (((-0.75*(phr[2]+phl[2]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*eyr[2]-0.75*eyl[2]-0.4330127018922193*(eyr[0]+eyl[0]))*chi; 
  incr[3] = ((0.4330127018922193*(phr[6]+phl[6])-0.25*phr[3]+0.25*phl[3])*c-0.4330127018922193*eyr[6]+0.4330127018922193*eyl[6]+0.25*(eyr[3]+eyl[3]))*chi; 
  incr[4] = (((-0.75*(phr[4]+phl[4]))+0.4330127018922193*phr[1]-0.4330127018922193*phl[1])*c+0.75*eyr[4]-0.75*eyl[4]-0.4330127018922193*(eyr[1]+eyl[1]))*chi; 
  incr[5] = ((0.4330127018922193*(phr[7]+phl[7])-0.25*phr[5]+0.25*phl[5])*c-0.4330127018922193*eyr[7]+0.4330127018922193*eyl[7]+0.25*(eyr[5]+eyl[5]))*chi; 
  incr[6] = (((-0.75*(phr[6]+phl[6]))+0.4330127018922193*phr[3]-0.4330127018922193*phl[3])*c+0.75*eyr[6]-0.75*eyl[6]-0.4330127018922193*(eyr[3]+eyl[3]))*chi; 
  incr[7] = (((-0.75*(phr[7]+phl[7]))+0.4330127018922193*phr[5]-0.4330127018922193*phl[5])*c+0.75*eyr[7]-0.75*eyl[7]-0.4330127018922193*(eyr[5]+eyl[5]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 
  outPhr[4] += incr[4]*dxr1; 
  outPhr[5] += incr[5]*dxr1; 
  outPhr[6] += incr[6]*dxr1; 
  outPhr[7] += incr[7]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += -1.0*incr[1]*dxl1; 
  outPhl[2] += incr[2]*dxl1; 
  outPhl[3] += -1.0*incr[3]*dxl1; 
  outPhl[4] += incr[4]*dxl1; 
  outPhl[5] += -1.0*incr[5]*dxl1; 
  outPhl[6] += incr[6]*dxl1; 
  outPhl[7] += incr[7]*dxl1; 

 
  incr[0] = (0.4330127018922193*(psr[2]+psl[2])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*byr[2])+0.4330127018922193*byl[2]+0.25*(byr[0]+byl[0]))*c2gamma; 
  incr[1] = (0.4330127018922193*(psr[4]+psl[4])-0.25*psr[1]+0.25*psl[1])*c*gamma+((-0.4330127018922193*byr[4])+0.4330127018922193*byl[4]+0.25*(byr[1]+byl[1]))*c2gamma; 
  incr[2] = ((-0.75*(psr[2]+psl[2]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*byr[2]-0.75*byl[2]-0.4330127018922193*(byr[0]+byl[0]))*c2gamma; 
  incr[3] = (0.4330127018922193*(psr[6]+psl[6])-0.25*psr[3]+0.25*psl[3])*c*gamma+((-0.4330127018922193*byr[6])+0.4330127018922193*byl[6]+0.25*(byr[3]+byl[3]))*c2gamma; 
  incr[4] = ((-0.75*(psr[4]+psl[4]))+0.4330127018922193*psr[1]-0.4330127018922193*psl[1])*c*gamma+(0.75*byr[4]-0.75*byl[4]-0.4330127018922193*(byr[1]+byl[1]))*c2gamma; 
  incr[5] = (0.4330127018922193*(psr[7]+psl[7])-0.25*psr[5]+0.25*psl[5])*c*gamma+((-0.4330127018922193*byr[7])+0.4330127018922193*byl[7]+0.25*(byr[5]+byl[5]))*c2gamma; 
  incr[6] = ((-0.75*(psr[6]+psl[6]))+0.4330127018922193*psr[3]-0.4330127018922193*psl[3])*c*gamma+(0.75*byr[6]-0.75*byl[6]-0.4330127018922193*(byr[3]+byl[3]))*c2gamma; 
  incr[7] = ((-0.75*(psr[7]+psl[7]))+0.4330127018922193*psr[5]-0.4330127018922193*psl[5])*c*gamma+(0.75*byr[7]-0.75*byl[7]-0.4330127018922193*(byr[5]+byl[5]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 
  outPsr[4] += incr[4]*dxr1; 
  outPsr[5] += incr[5]*dxr1; 
  outPsr[6] += incr[6]*dxr1; 
  outPsr[7] += incr[7]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += -1.0*incr[1]*dxl1; 
  outPsl[2] += incr[2]*dxl1; 
  outPsl[3] += -1.0*incr[3]*dxl1; 
  outPsl[4] += incr[4]*dxl1; 
  outPsl[5] += -1.0*incr[5]*dxl1; 
  outPsl[6] += incr[6]*dxl1; 
  outPsl[7] += incr[7]*dxl1; 

 
  return std::fmax(c, tau); 
} 
double MaxwellSurf3xSer_Z_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[2]; 
  const double dxr1 = 2.0/dxr[2]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[8]; 
  const double *ezl = &ql[16]; 
  const double *bxl = &ql[24]; 
  const double *byl = &ql[32]; 
  const double *bzl = &ql[40]; 
  const double *phl = &ql[48]; 
  const double *psl = &ql[56]; 
 
  double *outExl = &outl[0]; 
  double *outEyl = &outl[8]; 
  double *outEzl = &outl[16]; 
  double *outBxl = &outl[24]; 
  double *outByl = &outl[32]; 
  double *outBzl = &outl[40]; 
  double *outPhl = &outl[48]; 
  double *outPsl = &outl[56]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[8]; 
  const double *ezr = &qr[16]; 
  const double *bxr = &qr[24]; 
  const double *byr = &qr[32]; 
  const double *bzr = &qr[40]; 
  const double *phr = &qr[48]; 
  const double *psr = &qr[56]; 
 
  double *outExr = &outr[0]; 
  double *outEyr = &outr[8]; 
  double *outEzr = &outr[16]; 
  double *outBxr = &outr[24]; 
  double *outByr = &outr[32]; 
  double *outBzr = &outr[40]; 
  double *outPhr = &outr[48]; 
  double *outPsr = &outr[56]; 
 
  double incr[8]; 
 
  incr[0] = (0.4330127018922193*(exr[3]+exl[3])-0.25*exr[0]+0.25*exl[0])*tau+((-0.4330127018922193*byr[3])+0.4330127018922193*byl[3]+0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = (0.4330127018922193*(exr[5]+exl[5])-0.25*exr[1]+0.25*exl[1])*tau+((-0.4330127018922193*byr[5])+0.4330127018922193*byl[5]+0.25*(byr[1]+byl[1]))*c2; 
  incr[2] = (0.4330127018922193*(exr[6]+exl[6])-0.25*exr[2]+0.25*exl[2])*tau+((-0.4330127018922193*byr[6])+0.4330127018922193*byl[6]+0.25*(byr[2]+byl[2]))*c2; 
  incr[3] = ((-0.75*(exr[3]+exl[3]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*tau+(0.75*byr[3]-0.75*byl[3]-0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[4] = (0.4330127018922193*(exr[7]+exl[7])-0.25*exr[4]+0.25*exl[4])*tau+((-0.4330127018922193*byr[7])+0.4330127018922193*byl[7]+0.25*(byr[4]+byl[4]))*c2; 
  incr[5] = ((-0.75*(exr[5]+exl[5]))+0.4330127018922193*exr[1]-0.4330127018922193*exl[1])*tau+(0.75*byr[5]-0.75*byl[5]-0.4330127018922193*(byr[1]+byl[1]))*c2; 
  incr[6] = ((-0.75*(exr[6]+exl[6]))+0.4330127018922193*exr[2]-0.4330127018922193*exl[2])*tau+(0.75*byr[6]-0.75*byl[6]-0.4330127018922193*(byr[2]+byl[2]))*c2; 
  incr[7] = ((-0.75*(exr[7]+exl[7]))+0.4330127018922193*exr[4]-0.4330127018922193*exl[4])*tau+(0.75*byr[7]-0.75*byl[7]-0.4330127018922193*(byr[4]+byl[4]))*c2; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 
  outExr[3] += incr[3]*dxr1; 
  outExr[4] += incr[4]*dxr1; 
  outExr[5] += incr[5]*dxr1; 
  outExr[6] += incr[6]*dxr1; 
  outExr[7] += incr[7]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += -1.0*incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 
  outExl[3] += incr[3]*dxl1; 
  outExl[4] += -1.0*incr[4]*dxl1; 
  outExl[5] += incr[5]*dxl1; 
  outExl[6] += incr[6]*dxl1; 
  outExl[7] += incr[7]*dxl1; 

 
  incr[0] = (0.4330127018922193*(eyr[3]+eyl[3])-0.25*eyr[0]+0.25*eyl[0])*tau+(0.4330127018922193*bxr[3]-0.4330127018922193*bxl[3]-0.25*(bxr[0]+bxl[0]))*c2; 
  incr[1] = (0.4330127018922193*(eyr[5]+eyl[5])-0.25*eyr[1]+0.25*eyl[1])*tau+(0.4330127018922193*bxr[5]-0.4330127018922193*bxl[5]-0.25*(bxr[1]+bxl[1]))*c2; 
  incr[2] = (0.4330127018922193*(eyr[6]+eyl[6])-0.25*eyr[2]+0.25*eyl[2])*tau+(0.4330127018922193*bxr[6]-0.4330127018922193*bxl[6]-0.25*(bxr[2]+bxl[2]))*c2; 
  incr[3] = ((-0.75*(eyr[3]+eyl[3]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+((-0.75*bxr[3])+0.75*bxl[3]+0.4330127018922193*(bxr[0]+bxl[0]))*c2; 
  incr[4] = (0.4330127018922193*(eyr[7]+eyl[7])-0.25*eyr[4]+0.25*eyl[4])*tau+(0.4330127018922193*bxr[7]-0.4330127018922193*bxl[7]-0.25*(bxr[4]+bxl[4]))*c2; 
  incr[5] = ((-0.75*(eyr[5]+eyl[5]))+0.4330127018922193*eyr[1]-0.4330127018922193*eyl[1])*tau+((-0.75*bxr[5])+0.75*bxl[5]+0.4330127018922193*(bxr[1]+bxl[1]))*c2; 
  incr[6] = ((-0.75*(eyr[6]+eyl[6]))+0.4330127018922193*eyr[2]-0.4330127018922193*eyl[2])*tau+((-0.75*bxr[6])+0.75*bxl[6]+0.4330127018922193*(bxr[2]+bxl[2]))*c2; 
  incr[7] = ((-0.75*(eyr[7]+eyl[7]))+0.4330127018922193*eyr[4]-0.4330127018922193*eyl[4])*tau+((-0.75*bxr[7])+0.75*bxl[7]+0.4330127018922193*(bxr[4]+bxl[4]))*c2; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 
  outEyr[3] += incr[3]*dxr1; 
  outEyr[4] += incr[4]*dxr1; 
  outEyr[5] += incr[5]*dxr1; 
  outEyr[6] += incr[6]*dxr1; 
  outEyr[7] += incr[7]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += -1.0*incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 
  outEyl[3] += incr[3]*dxl1; 
  outEyl[4] += -1.0*incr[4]*dxl1; 
  outEyl[5] += incr[5]*dxl1; 
  outEyl[6] += incr[6]*dxl1; 
  outEyl[7] += incr[7]*dxl1; 

 
  incr[0] = (0.4330127018922193*(ezr[3]+ezl[3])-0.25*ezr[0]+0.25*ezl[0])*c*chi+((-0.4330127018922193*phr[3])+0.4330127018922193*phl[3]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.4330127018922193*(ezr[5]+ezl[5])-0.25*ezr[1]+0.25*ezl[1])*c*chi+((-0.4330127018922193*phr[5])+0.4330127018922193*phl[5]+0.25*(phr[1]+phl[1]))*c2chi; 
  incr[2] = (0.4330127018922193*(ezr[6]+ezl[6])-0.25*ezr[2]+0.25*ezl[2])*c*chi+((-0.4330127018922193*phr[6])+0.4330127018922193*phl[6]+0.25*(phr[2]+phl[2]))*c2chi; 
  incr[3] = ((-0.75*(ezr[3]+ezl[3]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*c*chi+(0.75*phr[3]-0.75*phl[3]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[4] = (0.4330127018922193*(ezr[7]+ezl[7])-0.25*ezr[4]+0.25*ezl[4])*c*chi+((-0.4330127018922193*phr[7])+0.4330127018922193*phl[7]+0.25*(phr[4]+phl[4]))*c2chi; 
  incr[5] = ((-0.75*(ezr[5]+ezl[5]))+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1])*c*chi+(0.75*phr[5]-0.75*phl[5]-0.4330127018922193*(phr[1]+phl[1]))*c2chi; 
  incr[6] = ((-0.75*(ezr[6]+ezl[6]))+0.4330127018922193*ezr[2]-0.4330127018922193*ezl[2])*c*chi+(0.75*phr[6]-0.75*phl[6]-0.4330127018922193*(phr[2]+phl[2]))*c2chi; 
  incr[7] = ((-0.75*(ezr[7]+ezl[7]))+0.4330127018922193*ezr[4]-0.4330127018922193*ezl[4])*c*chi+(0.75*phr[7]-0.75*phl[7]-0.4330127018922193*(phr[4]+phl[4]))*c2chi; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 
  outEzr[3] += incr[3]*dxr1; 
  outEzr[4] += incr[4]*dxr1; 
  outEzr[5] += incr[5]*dxr1; 
  outEzr[6] += incr[6]*dxr1; 
  outEzr[7] += incr[7]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += -1.0*incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 
  outEzl[3] += incr[3]*dxl1; 
  outEzl[4] += -1.0*incr[4]*dxl1; 
  outEzl[5] += incr[5]*dxl1; 
  outEzl[6] += incr[6]*dxl1; 
  outEzl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bxr[3]+bxl[3])-0.25*bxr[0]+0.25*bxl[0])*c2)/tau+0.4330127018922193*eyr[3]-0.4330127018922193*eyl[3]-0.25*(eyr[0]+eyl[0]); 
  incr[1] = ((0.4330127018922193*(bxr[5]+bxl[5])-0.25*bxr[1]+0.25*bxl[1])*c2)/tau+0.4330127018922193*eyr[5]-0.4330127018922193*eyl[5]-0.25*(eyr[1]+eyl[1]); 
  incr[2] = ((0.4330127018922193*(bxr[6]+bxl[6])-0.25*bxr[2]+0.25*bxl[2])*c2)/tau+0.4330127018922193*eyr[6]-0.4330127018922193*eyl[6]-0.25*(eyr[2]+eyl[2]); 
  incr[3] = (((-0.75*(bxr[3]+bxl[3]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c2)/tau-0.75*eyr[3]+0.75*eyl[3]+0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[4] = ((0.4330127018922193*(bxr[7]+bxl[7])-0.25*bxr[4]+0.25*bxl[4])*c2)/tau+0.4330127018922193*eyr[7]-0.4330127018922193*eyl[7]-0.25*(eyr[4]+eyl[4]); 
  incr[5] = (((-0.75*(bxr[5]+bxl[5]))+0.4330127018922193*bxr[1]-0.4330127018922193*bxl[1])*c2)/tau-0.75*eyr[5]+0.75*eyl[5]+0.4330127018922193*(eyr[1]+eyl[1]); 
  incr[6] = (((-0.75*(bxr[6]+bxl[6]))+0.4330127018922193*bxr[2]-0.4330127018922193*bxl[2])*c2)/tau-0.75*eyr[6]+0.75*eyl[6]+0.4330127018922193*(eyr[2]+eyl[2]); 
  incr[7] = (((-0.75*(bxr[7]+bxl[7]))+0.4330127018922193*bxr[4]-0.4330127018922193*bxl[4])*c2)/tau-0.75*eyr[7]+0.75*eyl[7]+0.4330127018922193*(eyr[4]+eyl[4]); 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 
  outBxr[3] += incr[3]*dxr1; 
  outBxr[4] += incr[4]*dxr1; 
  outBxr[5] += incr[5]*dxr1; 
  outBxr[6] += incr[6]*dxr1; 
  outBxr[7] += incr[7]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += -1.0*incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 
  outBxl[3] += incr[3]*dxl1; 
  outBxl[4] += -1.0*incr[4]*dxl1; 
  outBxl[5] += incr[5]*dxl1; 
  outBxl[6] += incr[6]*dxl1; 
  outBxl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(byr[3]+byl[3])-0.25*byr[0]+0.25*byl[0])*c2)/tau-0.4330127018922193*exr[3]+0.4330127018922193*exl[3]+0.25*(exr[0]+exl[0]); 
  incr[1] = ((0.4330127018922193*(byr[5]+byl[5])-0.25*byr[1]+0.25*byl[1])*c2)/tau-0.4330127018922193*exr[5]+0.4330127018922193*exl[5]+0.25*(exr[1]+exl[1]); 
  incr[2] = ((0.4330127018922193*(byr[6]+byl[6])-0.25*byr[2]+0.25*byl[2])*c2)/tau-0.4330127018922193*exr[6]+0.4330127018922193*exl[6]+0.25*(exr[2]+exl[2]); 
  incr[3] = (((-0.75*(byr[3]+byl[3]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau+0.75*exr[3]-0.75*exl[3]-0.4330127018922193*(exr[0]+exl[0]); 
  incr[4] = ((0.4330127018922193*(byr[7]+byl[7])-0.25*byr[4]+0.25*byl[4])*c2)/tau-0.4330127018922193*exr[7]+0.4330127018922193*exl[7]+0.25*(exr[4]+exl[4]); 
  incr[5] = (((-0.75*(byr[5]+byl[5]))+0.4330127018922193*byr[1]-0.4330127018922193*byl[1])*c2)/tau+0.75*exr[5]-0.75*exl[5]-0.4330127018922193*(exr[1]+exl[1]); 
  incr[6] = (((-0.75*(byr[6]+byl[6]))+0.4330127018922193*byr[2]-0.4330127018922193*byl[2])*c2)/tau+0.75*exr[6]-0.75*exl[6]-0.4330127018922193*(exr[2]+exl[2]); 
  incr[7] = (((-0.75*(byr[7]+byl[7]))+0.4330127018922193*byr[4]-0.4330127018922193*byl[4])*c2)/tau+0.75*exr[7]-0.75*exl[7]-0.4330127018922193*(exr[4]+exl[4]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 
  outByr[3] += incr[3]*dxr1; 
  outByr[4] += incr[4]*dxr1; 
  outByr[5] += incr[5]*dxr1; 
  outByr[6] += incr[6]*dxr1; 
  outByr[7] += incr[7]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += -1.0*incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 
  outByl[3] += incr[3]*dxl1; 
  outByl[4] += -1.0*incr[4]*dxl1; 
  outByl[5] += incr[5]*dxl1; 
  outByl[6] += incr[6]*dxl1; 
  outByl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(bzr[3]+bzl[3])-0.25*bzr[0]+0.25*bzl[0])*c-0.4330127018922193*psr[3]+0.4330127018922193*psl[3]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = ((0.4330127018922193*(bzr[5]+bzl[5])-0.25*bzr[1]+0.25*bzl[1])*c-0.4330127018922193*psr[5]+0.4330127018922193*psl[5]+0.25*(psr[1]+psl[1]))*gamma; 
  incr[2] = ((0.4330127018922193*(bzr[6]+bzl[6])-0.25*bzr[2]+0.25*bzl[2])*c-0.4330127018922193*psr[6]+0.4330127018922193*psl[6]+0.25*(psr[2]+psl[2]))*gamma; 
  incr[3] = (((-0.75*(bzr[3]+bzl[3]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c+0.75*psr[3]-0.75*psl[3]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[4] = ((0.4330127018922193*(bzr[7]+bzl[7])-0.25*bzr[4]+0.25*bzl[4])*c-0.4330127018922193*psr[7]+0.4330127018922193*psl[7]+0.25*(psr[4]+psl[4]))*gamma; 
  incr[5] = (((-0.75*(bzr[5]+bzl[5]))+0.4330127018922193*bzr[1]-0.4330127018922193*bzl[1])*c+0.75*psr[5]-0.75*psl[5]-0.4330127018922193*(psr[1]+psl[1]))*gamma; 
  incr[6] = (((-0.75*(bzr[6]+bzl[6]))+0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2])*c+0.75*psr[6]-0.75*psl[6]-0.4330127018922193*(psr[2]+psl[2]))*gamma; 
  incr[7] = (((-0.75*(bzr[7]+bzl[7]))+0.4330127018922193*bzr[4]-0.4330127018922193*bzl[4])*c+0.75*psr[7]-0.75*psl[7]-0.4330127018922193*(psr[4]+psl[4]))*gamma; 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 
  outBzr[3] += incr[3]*dxr1; 
  outBzr[4] += incr[4]*dxr1; 
  outBzr[5] += incr[5]*dxr1; 
  outBzr[6] += incr[6]*dxr1; 
  outBzr[7] += incr[7]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += -1.0*incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 
  outBzl[3] += incr[3]*dxl1; 
  outBzl[4] += -1.0*incr[4]*dxl1; 
  outBzl[5] += incr[5]*dxl1; 
  outBzl[6] += incr[6]*dxl1; 
  outBzl[7] += incr[7]*dxl1; 

 
  incr[0] = ((0.4330127018922193*(phr[3]+phl[3])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*ezr[3]+0.4330127018922193*ezl[3]+0.25*(ezr[0]+ezl[0]))*chi; 
  incr[1] = ((0.4330127018922193*(phr[5]+phl[5])-0.25*phr[1]+0.25*phl[1])*c-0.4330127018922193*ezr[5]+0.4330127018922193*ezl[5]+0.25*(ezr[1]+ezl[1]))*chi; 
  incr[2] = ((0.4330127018922193*(phr[6]+phl[6])-0.25*phr[2]+0.25*phl[2])*c-0.4330127018922193*ezr[6]+0.4330127018922193*ezl[6]+0.25*(ezr[2]+ezl[2]))*chi; 
  incr[3] = (((-0.75*(phr[3]+phl[3]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*ezr[3]-0.75*ezl[3]-0.4330127018922193*(ezr[0]+ezl[0]))*chi; 
  incr[4] = ((0.4330127018922193*(phr[7]+phl[7])-0.25*phr[4]+0.25*phl[4])*c-0.4330127018922193*ezr[7]+0.4330127018922193*ezl[7]+0.25*(ezr[4]+ezl[4]))*chi; 
  incr[5] = (((-0.75*(phr[5]+phl[5]))+0.4330127018922193*phr[1]-0.4330127018922193*phl[1])*c+0.75*ezr[5]-0.75*ezl[5]-0.4330127018922193*(ezr[1]+ezl[1]))*chi; 
  incr[6] = (((-0.75*(phr[6]+phl[6]))+0.4330127018922193*phr[2]-0.4330127018922193*phl[2])*c+0.75*ezr[6]-0.75*ezl[6]-0.4330127018922193*(ezr[2]+ezl[2]))*chi; 
  incr[7] = (((-0.75*(phr[7]+phl[7]))+0.4330127018922193*phr[4]-0.4330127018922193*phl[4])*c+0.75*ezr[7]-0.75*ezl[7]-0.4330127018922193*(ezr[4]+ezl[4]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 
  outPhr[3] += incr[3]*dxr1; 
  outPhr[4] += incr[4]*dxr1; 
  outPhr[5] += incr[5]*dxr1; 
  outPhr[6] += incr[6]*dxr1; 
  outPhr[7] += incr[7]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += -1.0*incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 
  outPhl[3] += incr[3]*dxl1; 
  outPhl[4] += -1.0*incr[4]*dxl1; 
  outPhl[5] += incr[5]*dxl1; 
  outPhl[6] += incr[6]*dxl1; 
  outPhl[7] += incr[7]*dxl1; 

 
  incr[0] = (0.4330127018922193*(psr[3]+psl[3])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*bzr[3])+0.4330127018922193*bzl[3]+0.25*(bzr[0]+bzl[0]))*c2gamma; 
  incr[1] = (0.4330127018922193*(psr[5]+psl[5])-0.25*psr[1]+0.25*psl[1])*c*gamma+((-0.4330127018922193*bzr[5])+0.4330127018922193*bzl[5]+0.25*(bzr[1]+bzl[1]))*c2gamma; 
  incr[2] = (0.4330127018922193*(psr[6]+psl[6])-0.25*psr[2]+0.25*psl[2])*c*gamma+((-0.4330127018922193*bzr[6])+0.4330127018922193*bzl[6]+0.25*(bzr[2]+bzl[2]))*c2gamma; 
  incr[3] = ((-0.75*(psr[3]+psl[3]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*bzr[3]-0.75*bzl[3]-0.4330127018922193*(bzr[0]+bzl[0]))*c2gamma; 
  incr[4] = (0.4330127018922193*(psr[7]+psl[7])-0.25*psr[4]+0.25*psl[4])*c*gamma+((-0.4330127018922193*bzr[7])+0.4330127018922193*bzl[7]+0.25*(bzr[4]+bzl[4]))*c2gamma; 
  incr[5] = ((-0.75*(psr[5]+psl[5]))+0.4330127018922193*psr[1]-0.4330127018922193*psl[1])*c*gamma+(0.75*bzr[5]-0.75*bzl[5]-0.4330127018922193*(bzr[1]+bzl[1]))*c2gamma; 
  incr[6] = ((-0.75*(psr[6]+psl[6]))+0.4330127018922193*psr[2]-0.4330127018922193*psl[2])*c*gamma+(0.75*bzr[6]-0.75*bzl[6]-0.4330127018922193*(bzr[2]+bzl[2]))*c2gamma; 
  incr[7] = ((-0.75*(psr[7]+psl[7]))+0.4330127018922193*psr[4]-0.4330127018922193*psl[4])*c*gamma+(0.75*bzr[7]-0.75*bzl[7]-0.4330127018922193*(bzr[4]+bzl[4]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 
  outPsr[3] += incr[3]*dxr1; 
  outPsr[4] += incr[4]*dxr1; 
  outPsr[5] += incr[5]*dxr1; 
  outPsr[6] += incr[6]*dxr1; 
  outPsr[7] += incr[7]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += -1.0*incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 
  outPsl[3] += incr[3]*dxl1; 
  outPsl[4] += -1.0*incr[4]*dxl1; 
  outPsl[5] += incr[5]*dxl1; 
  outPsl[6] += incr[6]*dxl1; 
  outPsl[7] += incr[7]*dxl1; 

 
  return std::fmax(c, tau); 
} 
