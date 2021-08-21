#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream1x1vSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double Ghat[3]; 
  if (wr[1]>0) { 
  Ghat[0] = dxvl[1]*(0.4564354645876383*fl[6]+0.3535533905932737*fl[3]+0.2041241452319315*fl[2])+wl[1]*(1.58113883008419*fl[4]+1.224744871391589*fl[1]+0.7071067811865475*fl[0]); 
  Ghat[1] = dxvl[1]*(0.3162277660168379*fl[7]+0.1825741858350554*fl[5]+0.4564354645876384*fl[4]+0.3535533905932737*fl[1]+0.2041241452319315*fl[0])+wl[1]*(1.58113883008419*fl[6]+1.224744871391589*fl[3]+0.7071067811865475*fl[2]); 
  Ghat[2] = wl[1]*(1.224744871391589*fl[7]+0.7071067811865475*fl[5])+dxvl[1]*(0.408248290463863*fl[6]+0.3162277660168379*fl[3]+0.1825741858350554*fl[2]); 

  outr[0] += 0.7071067811865475*Ghat[0]*rdxr2; 
  outr[1] += -1.224744871391589*Ghat[0]*rdxr2; 
  outr[2] += 0.7071067811865475*Ghat[1]*rdxr2; 
  outr[3] += -1.224744871391589*Ghat[1]*rdxr2; 
  outr[4] += 1.58113883008419*Ghat[0]*rdxr2; 
  outr[5] += 0.7071067811865475*Ghat[2]*rdxr2; 
  outr[6] += 1.58113883008419*Ghat[1]*rdxr2; 
  outr[7] += -1.224744871391589*Ghat[2]*rdxr2; 

  outl[0] += -0.7071067811865475*Ghat[0]*rdxl2; 
  outl[1] += -1.224744871391589*Ghat[0]*rdxl2; 
  outl[2] += -0.7071067811865475*Ghat[1]*rdxl2; 
  outl[3] += -1.224744871391589*Ghat[1]*rdxl2; 
  outl[4] += -1.58113883008419*Ghat[0]*rdxl2; 
  outl[5] += -0.7071067811865475*Ghat[2]*rdxl2; 
  outl[6] += -1.58113883008419*Ghat[1]*rdxl2; 
  outl[7] += -1.224744871391589*Ghat[2]*rdxl2; 
  } else { 
  Ghat[0] = dxvr[1]*(0.4564354645876383*fr[6]-0.3535533905932737*fr[3]+0.2041241452319315*fr[2])+wr[1]*(1.58113883008419*fr[4]-1.224744871391589*fr[1]+0.7071067811865475*fr[0]); 
  Ghat[1] = dxvr[1]*((-0.3162277660168379*fr[7])+0.1825741858350554*fr[5]+0.4564354645876384*fr[4]-0.3535533905932737*fr[1]+0.2041241452319315*fr[0])+wr[1]*(1.58113883008419*fr[6]-1.224744871391589*fr[3]+0.7071067811865475*fr[2]); 
  Ghat[2] = wr[1]*(0.7071067811865475*fr[5]-1.224744871391589*fr[7])+dxvr[1]*(0.408248290463863*fr[6]-0.3162277660168379*fr[3]+0.1825741858350554*fr[2]); 

  outr[0] += 0.7071067811865475*Ghat[0]*rdxr2; 
  outr[1] += -1.224744871391589*Ghat[0]*rdxr2; 
  outr[2] += 0.7071067811865475*Ghat[1]*rdxr2; 
  outr[3] += -1.224744871391589*Ghat[1]*rdxr2; 
  outr[4] += 1.58113883008419*Ghat[0]*rdxr2; 
  outr[5] += 0.7071067811865475*Ghat[2]*rdxr2; 
  outr[6] += 1.58113883008419*Ghat[1]*rdxr2; 
  outr[7] += -1.224744871391589*Ghat[2]*rdxr2; 

  outl[0] += -0.7071067811865475*Ghat[0]*rdxl2; 
  outl[1] += -1.224744871391589*Ghat[0]*rdxl2; 
  outl[2] += -0.7071067811865475*Ghat[1]*rdxl2; 
  outl[3] += -1.224744871391589*Ghat[1]*rdxl2; 
  outl[4] += -1.58113883008419*Ghat[0]*rdxl2; 
  outl[5] += -0.7071067811865475*Ghat[2]*rdxl2; 
  outl[6] += -1.58113883008419*Ghat[1]*rdxl2; 
  outl[7] += -1.224744871391589*Ghat[2]*rdxl2; 
  } 
} 
