#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream1x2vSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double Ghat[4]; 
  if (wr[1]>0) { 
  Ghat[0] = dxvl[1]*(0.3535533905932737*fl[4]+0.2041241452319315*fl[2])+(1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wl[1]; 
  Ghat[1] = wl[1]*(1.224744871391589*fl[4]+0.7071067811865475*fl[2])+dxvl[1]*(0.3535533905932737*fl[1]+0.2041241452319315*fl[0]); 
  Ghat[2] = dxvl[1]*(0.3535533905932737*fl[7]+0.2041241452319315*fl[6])+wl[1]*(1.224744871391589*fl[5]+0.7071067811865475*fl[3]); 
  Ghat[3] = wl[1]*(1.224744871391589*fl[7]+0.7071067811865475*fl[6])+dxvl[1]*(0.3535533905932737*fl[5]+0.2041241452319315*fl[3]); 

  outr[0] += 0.7071067811865475*Ghat[0]*rdxr2; 
  outr[1] += -1.224744871391589*Ghat[0]*rdxr2; 
  outr[2] += 0.7071067811865475*Ghat[1]*rdxr2; 
  outr[3] += 0.7071067811865475*Ghat[2]*rdxr2; 
  outr[4] += -1.224744871391589*Ghat[1]*rdxr2; 
  outr[5] += -1.224744871391589*Ghat[2]*rdxr2; 
  outr[6] += 0.7071067811865475*Ghat[3]*rdxr2; 
  outr[7] += -1.224744871391589*Ghat[3]*rdxr2; 

  outl[0] += -0.7071067811865475*Ghat[0]*rdxl2; 
  outl[1] += -1.224744871391589*Ghat[0]*rdxl2; 
  outl[2] += -0.7071067811865475*Ghat[1]*rdxl2; 
  outl[3] += -0.7071067811865475*Ghat[2]*rdxl2; 
  outl[4] += -1.224744871391589*Ghat[1]*rdxl2; 
  outl[5] += -1.224744871391589*Ghat[2]*rdxl2; 
  outl[6] += -0.7071067811865475*Ghat[3]*rdxl2; 
  outl[7] += -1.224744871391589*Ghat[3]*rdxl2; 
  } else { 
  Ghat[0] = dxvr[1]*(0.2041241452319315*fr[2]-0.3535533905932737*fr[4])+(0.7071067811865475*fr[0]-1.224744871391589*fr[1])*wr[1]; 
  Ghat[1] = wr[1]*(0.7071067811865475*fr[2]-1.224744871391589*fr[4])+dxvr[1]*(0.2041241452319315*fr[0]-0.3535533905932737*fr[1]); 
  Ghat[2] = dxvr[1]*(0.2041241452319315*fr[6]-0.3535533905932737*fr[7])+wr[1]*(0.7071067811865475*fr[3]-1.224744871391589*fr[5]); 
  Ghat[3] = wr[1]*(0.7071067811865475*fr[6]-1.224744871391589*fr[7])+dxvr[1]*(0.2041241452319315*fr[3]-0.3535533905932737*fr[5]); 

  outr[0] += 0.7071067811865475*Ghat[0]*rdxr2; 
  outr[1] += -1.224744871391589*Ghat[0]*rdxr2; 
  outr[2] += 0.7071067811865475*Ghat[1]*rdxr2; 
  outr[3] += 0.7071067811865475*Ghat[2]*rdxr2; 
  outr[4] += -1.224744871391589*Ghat[1]*rdxr2; 
  outr[5] += -1.224744871391589*Ghat[2]*rdxr2; 
  outr[6] += 0.7071067811865475*Ghat[3]*rdxr2; 
  outr[7] += -1.224744871391589*Ghat[3]*rdxr2; 

  outl[0] += -0.7071067811865475*Ghat[0]*rdxl2; 
  outl[1] += -1.224744871391589*Ghat[0]*rdxl2; 
  outl[2] += -0.7071067811865475*Ghat[1]*rdxl2; 
  outl[3] += -0.7071067811865475*Ghat[2]*rdxl2; 
  outl[4] += -1.224744871391589*Ghat[1]*rdxl2; 
  outl[5] += -1.224744871391589*Ghat[2]*rdxl2; 
  outl[6] += -0.7071067811865475*Ghat[3]*rdxl2; 
  outl[7] += -1.224744871391589*Ghat[3]*rdxl2; 
  } 
} 
