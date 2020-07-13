#include <MaxwellModDecl.h> 
__host__ __device__ double MaxwellSurf2xMax_X_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
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
 
  incr[0] = (0.4330127018922193*(exr[1]+exl[1])-0.25*exr[0]+0.25*exl[0])*c*chi+((-0.4330127018922193*phr[1])+0.4330127018922193*phl[1]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = ((-0.75*(exr[1]+exl[1]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*c*chi+(0.75*phr[1]-0.75*phl[1]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 
  incr[2] = (0.25*exl[2]-0.25*exr[2])*c*chi+0.25*(phr[2]+phl[2])*c2chi; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += incr[1]*dxl1; 
  outExl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = (0.4330127018922193*(eyr[1]+eyl[1])-0.25*eyr[0]+0.25*eyl[0])*tau+((-0.4330127018922193*bzr[1])+0.4330127018922193*bzl[1]+0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = ((-0.75*(eyr[1]+eyl[1]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*tau+(0.75*bzr[1]-0.75*bzl[1]-0.4330127018922193*(bzr[0]+bzl[0]))*c2; 
  incr[2] = 0.25*(eyl[2]*tau+(bzr[2]+bzl[2])*c2)-0.25*eyr[2]*tau; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += incr[1]*dxl1; 
  outEyl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = (0.4330127018922193*(ezr[1]+ezl[1])-0.25*ezr[0]+0.25*ezl[0])*tau+(0.4330127018922193*byr[1]-0.4330127018922193*byl[1]-0.25*(byr[0]+byl[0]))*c2; 
  incr[1] = ((-0.75*(ezr[1]+ezl[1]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+((-0.75*byr[1])+0.75*byl[1]+0.4330127018922193*(byr[0]+byl[0]))*c2; 
  incr[2] = (0.25*ezl[2]-0.25*ezr[2])*tau-0.25*(byr[2]+byl[2])*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += incr[1]*dxl1; 
  outEzl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = ((0.4330127018922193*(bxr[1]+bxl[1])-0.25*bxr[0]+0.25*bxl[0])*c-0.4330127018922193*psr[1]+0.4330127018922193*psl[1]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (((-0.75*(bxr[1]+bxl[1]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c+0.75*psr[1]-0.75*psl[1]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 
  incr[2] = (0.25*(bxl[2]*c+psr[2]+psl[2])-0.25*bxr[2]*c)*gamma; 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += incr[1]*dxl1; 
  outBxl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = ((0.4330127018922193*(byr[1]+byl[1])-0.25*byr[0]+0.25*byl[0])*c2)/tau+0.4330127018922193*ezr[1]-0.4330127018922193*ezl[1]-0.25*(ezr[0]+ezl[0]); 
  incr[1] = (((-0.75*(byr[1]+byl[1]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c2)/tau-0.75*ezr[1]+0.75*ezl[1]+0.4330127018922193*(ezr[0]+ezl[0]); 
  incr[2] = ((0.25*byl[2]-0.25*byr[2])*c2)/tau-0.25*(ezr[2]+ezl[2]); 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += incr[1]*dxl1; 
  outByl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = ((0.4330127018922193*(bzr[1]+bzl[1])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau-0.4330127018922193*eyr[1]+0.4330127018922193*eyl[1]+0.25*(eyr[0]+eyl[0]); 
  incr[1] = (((-0.75*(bzr[1]+bzl[1]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau+0.75*eyr[1]-0.75*eyl[1]-0.4330127018922193*(eyr[0]+eyl[0]); 
  incr[2] = ((0.25*bzl[2]-0.25*bzr[2])*c2)/tau+0.25*(eyr[2]+eyl[2]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += incr[1]*dxl1; 
  outBzl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = ((0.4330127018922193*(phr[1]+phl[1])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*exr[1]+0.4330127018922193*exl[1]+0.25*(exr[0]+exl[0]))*chi; 
  incr[1] = (((-0.75*(phr[1]+phl[1]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*exr[1]-0.75*exl[1]-0.4330127018922193*(exr[0]+exl[0]))*chi; 
  incr[2] = (0.25*(phl[2]*c+exr[2]+exl[2])-0.25*phr[2]*c)*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += incr[1]*dxl1; 
  outPhl[2] += -1.0*incr[2]*dxl1; 

  incr[0] = (0.4330127018922193*(psr[1]+psl[1])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*bxr[1])+0.4330127018922193*bxl[1]+0.25*(bxr[0]+bxl[0]))*c2gamma; 
  incr[1] = ((-0.75*(psr[1]+psl[1]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*bxr[1]-0.75*bxl[1]-0.4330127018922193*(bxr[0]+bxl[0]))*c2gamma; 
  incr[2] = (0.25*psl[2]-0.25*psr[2])*c*gamma+0.25*(bxr[2]+bxl[2])*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += incr[1]*dxl1; 
  outPsl[2] += -1.0*incr[2]*dxl1; 

  return std::fmax(c, tau); 
} 
__host__ __device__ double MaxwellSurf2xMax_Y_P1(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dxl1 = 2.0/dxl[1]; 
  const double dxr1 = 2.0/dxr[1]; 
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
 
  incr[0] = (0.4330127018922193*(exr[2]+exl[2])-0.25*exr[0]+0.25*exl[0])*tau+(0.4330127018922193*bzr[2]-0.4330127018922193*bzl[2]-0.25*(bzr[0]+bzl[0]))*c2; 
  incr[1] = (0.25*exl[1]-0.25*exr[1])*tau-0.25*(bzr[1]+bzl[1])*c2; 
  incr[2] = ((-0.75*(exr[2]+exl[2]))+0.4330127018922193*exr[0]-0.4330127018922193*exl[0])*tau+((-0.75*bzr[2])+0.75*bzl[2]+0.4330127018922193*(bzr[0]+bzl[0]))*c2; 

  outExr[0] += incr[0]*dxr1; 
  outExr[1] += incr[1]*dxr1; 
  outExr[2] += incr[2]*dxr1; 

  outExl[0] += -1.0*incr[0]*dxl1; 
  outExl[1] += -1.0*incr[1]*dxl1; 
  outExl[2] += incr[2]*dxl1; 

  incr[0] = (0.4330127018922193*(eyr[2]+eyl[2])-0.25*eyr[0]+0.25*eyl[0])*c*chi+((-0.4330127018922193*phr[2])+0.4330127018922193*phl[2]+0.25*(phr[0]+phl[0]))*c2chi; 
  incr[1] = (0.25*eyl[1]-0.25*eyr[1])*c*chi+0.25*(phr[1]+phl[1])*c2chi; 
  incr[2] = ((-0.75*(eyr[2]+eyl[2]))+0.4330127018922193*eyr[0]-0.4330127018922193*eyl[0])*c*chi+(0.75*phr[2]-0.75*phl[2]-0.4330127018922193*(phr[0]+phl[0]))*c2chi; 

  outEyr[0] += incr[0]*dxr1; 
  outEyr[1] += incr[1]*dxr1; 
  outEyr[2] += incr[2]*dxr1; 

  outEyl[0] += -1.0*incr[0]*dxl1; 
  outEyl[1] += -1.0*incr[1]*dxl1; 
  outEyl[2] += incr[2]*dxl1; 

  incr[0] = (0.4330127018922193*(ezr[2]+ezl[2])-0.25*ezr[0]+0.25*ezl[0])*tau+((-0.4330127018922193*bxr[2])+0.4330127018922193*bxl[2]+0.25*(bxr[0]+bxl[0]))*c2; 
  incr[1] = 0.25*(ezl[1]*tau+(bxr[1]+bxl[1])*c2)-0.25*ezr[1]*tau; 
  incr[2] = ((-0.75*(ezr[2]+ezl[2]))+0.4330127018922193*ezr[0]-0.4330127018922193*ezl[0])*tau+(0.75*bxr[2]-0.75*bxl[2]-0.4330127018922193*(bxr[0]+bxl[0]))*c2; 

  outEzr[0] += incr[0]*dxr1; 
  outEzr[1] += incr[1]*dxr1; 
  outEzr[2] += incr[2]*dxr1; 

  outEzl[0] += -1.0*incr[0]*dxl1; 
  outEzl[1] += -1.0*incr[1]*dxl1; 
  outEzl[2] += incr[2]*dxl1; 

  incr[0] = ((0.4330127018922193*(bxr[2]+bxl[2])-0.25*bxr[0]+0.25*bxl[0])*c2)/tau-0.4330127018922193*ezr[2]+0.4330127018922193*ezl[2]+0.25*(ezr[0]+ezl[0]); 
  incr[1] = ((0.25*bxl[1]-0.25*bxr[1])*c2)/tau+0.25*(ezr[1]+ezl[1]); 
  incr[2] = (((-0.75*(bxr[2]+bxl[2]))+0.4330127018922193*bxr[0]-0.4330127018922193*bxl[0])*c2)/tau+0.75*ezr[2]-0.75*ezl[2]-0.4330127018922193*(ezr[0]+ezl[0]); 

  outBxr[0] += incr[0]*dxr1; 
  outBxr[1] += incr[1]*dxr1; 
  outBxr[2] += incr[2]*dxr1; 

  outBxl[0] += -1.0*incr[0]*dxl1; 
  outBxl[1] += -1.0*incr[1]*dxl1; 
  outBxl[2] += incr[2]*dxl1; 

  incr[0] = ((0.4330127018922193*(byr[2]+byl[2])-0.25*byr[0]+0.25*byl[0])*c-0.4330127018922193*psr[2]+0.4330127018922193*psl[2]+0.25*(psr[0]+psl[0]))*gamma; 
  incr[1] = (0.25*(byl[1]*c+psr[1]+psl[1])-0.25*byr[1]*c)*gamma; 
  incr[2] = (((-0.75*(byr[2]+byl[2]))+0.4330127018922193*byr[0]-0.4330127018922193*byl[0])*c+0.75*psr[2]-0.75*psl[2]-0.4330127018922193*(psr[0]+psl[0]))*gamma; 

  outByr[0] += incr[0]*dxr1; 
  outByr[1] += incr[1]*dxr1; 
  outByr[2] += incr[2]*dxr1; 

  outByl[0] += -1.0*incr[0]*dxl1; 
  outByl[1] += -1.0*incr[1]*dxl1; 
  outByl[2] += incr[2]*dxl1; 

  incr[0] = ((0.4330127018922193*(bzr[2]+bzl[2])-0.25*bzr[0]+0.25*bzl[0])*c2)/tau+0.4330127018922193*exr[2]-0.4330127018922193*exl[2]-0.25*(exr[0]+exl[0]); 
  incr[1] = ((0.25*bzl[1]-0.25*bzr[1])*c2)/tau-0.25*(exr[1]+exl[1]); 
  incr[2] = (((-0.75*(bzr[2]+bzl[2]))+0.4330127018922193*bzr[0]-0.4330127018922193*bzl[0])*c2)/tau-0.75*exr[2]+0.75*exl[2]+0.4330127018922193*(exr[0]+exl[0]); 

  outBzr[0] += incr[0]*dxr1; 
  outBzr[1] += incr[1]*dxr1; 
  outBzr[2] += incr[2]*dxr1; 

  outBzl[0] += -1.0*incr[0]*dxl1; 
  outBzl[1] += -1.0*incr[1]*dxl1; 
  outBzl[2] += incr[2]*dxl1; 

  incr[0] = ((0.4330127018922193*(phr[2]+phl[2])-0.25*phr[0]+0.25*phl[0])*c-0.4330127018922193*eyr[2]+0.4330127018922193*eyl[2]+0.25*(eyr[0]+eyl[0]))*chi; 
  incr[1] = (0.25*(phl[1]*c+eyr[1]+eyl[1])-0.25*phr[1]*c)*chi; 
  incr[2] = (((-0.75*(phr[2]+phl[2]))+0.4330127018922193*phr[0]-0.4330127018922193*phl[0])*c+0.75*eyr[2]-0.75*eyl[2]-0.4330127018922193*(eyr[0]+eyl[0]))*chi; 

  outPhr[0] += incr[0]*dxr1; 
  outPhr[1] += incr[1]*dxr1; 
  outPhr[2] += incr[2]*dxr1; 

  outPhl[0] += -1.0*incr[0]*dxl1; 
  outPhl[1] += -1.0*incr[1]*dxl1; 
  outPhl[2] += incr[2]*dxl1; 

  incr[0] = (0.4330127018922193*(psr[2]+psl[2])-0.25*psr[0]+0.25*psl[0])*c*gamma+((-0.4330127018922193*byr[2])+0.4330127018922193*byl[2]+0.25*(byr[0]+byl[0]))*c2gamma; 
  incr[1] = (0.25*psl[1]-0.25*psr[1])*c*gamma+0.25*(byr[1]+byl[1])*c2gamma; 
  incr[2] = ((-0.75*(psr[2]+psl[2]))+0.4330127018922193*psr[0]-0.4330127018922193*psl[0])*c*gamma+(0.75*byr[2]-0.75*byl[2]-0.4330127018922193*(byr[0]+byl[0]))*c2gamma; 

  outPsr[0] += incr[0]*dxr1; 
  outPsr[1] += incr[1]*dxr1; 
  outPsr[2] += incr[2]*dxr1; 

  outPsl[0] += -1.0*incr[0]*dxl1; 
  outPsl[1] += -1.0*incr[1]*dxl1; 
  outPsl[2] += incr[2]*dxl1; 

  return std::fmax(c, tau); 
} 
