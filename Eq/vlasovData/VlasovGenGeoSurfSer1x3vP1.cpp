#include <VlasovModDecl.h>
__host__ __device__ void VlasovGenGeoSurf1x3vSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *alphaGeo, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // alphaGeo:  General geometry alph.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.

  const double *ag0 = &alphaGeo[0]; 

  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double alpha[8]; 
  double Ghat[8]; 

  alpha[0] = 1.224744871391589*ag0[1]+0.7071067811865475*ag0[0]; 
  alpha[1] = 1.224744871391589*ag0[5]+0.7071067811865475*ag0[2]; 
  alpha[2] = 1.224744871391589*ag0[6]+0.7071067811865475*ag0[3]; 
  alpha[3] = 1.224744871391589*ag0[8]+0.7071067811865475*ag0[4]; 
  alpha[4] = 1.224744871391589*ag0[11]+0.7071067811865475*ag0[7]; 
  alpha[5] = 1.224744871391589*ag0[12]+0.7071067811865475*ag0[9]; 
  alpha[6] = 1.224744871391589*ag0[13]+0.7071067811865475*ag0[10]; 
  alpha[7] = 1.224744871391589*ag0[15]+0.7071067811865475*ag0[14]; 

  if (wr[1]>0) { 
  Ghat[0] = 0.4330127018922193*alpha[7]*fl[15]+0.25*alpha[7]*fl[14]+0.4330127018922193*alpha[6]*fl[13]+0.4330127018922193*alpha[5]*fl[12]+0.4330127018922193*alpha[4]*fl[11]+0.25*alpha[6]*fl[10]+0.25*alpha[5]*fl[9]+0.4330127018922193*alpha[3]*fl[8]+0.25*alpha[4]*fl[7]+0.4330127018922193*alpha[2]*fl[6]+0.4330127018922193*alpha[1]*fl[5]+0.25*alpha[3]*fl[4]+0.25*alpha[2]*fl[3]+0.25*alpha[1]*fl[2]+0.4330127018922193*alpha[0]*fl[1]+0.25*alpha[0]*fl[0]; 
  Ghat[1] = 0.4330127018922193*alpha[6]*fl[15]+0.25*alpha[6]*fl[14]+0.4330127018922193*alpha[7]*fl[13]+0.4330127018922193*alpha[3]*fl[12]+0.4330127018922193*alpha[2]*fl[11]+0.25*alpha[7]*fl[10]+0.25*alpha[3]*fl[9]+0.4330127018922193*alpha[5]*fl[8]+0.25*alpha[2]*fl[7]+0.4330127018922193*alpha[4]*fl[6]+0.4330127018922193*alpha[0]*fl[5]+0.25*fl[4]*alpha[5]+0.25*fl[3]*alpha[4]+0.25*alpha[0]*fl[2]+0.4330127018922193*alpha[1]*fl[1]+0.25*fl[0]*alpha[1]; 
  Ghat[2] = 0.4330127018922193*alpha[5]*fl[15]+0.25*alpha[5]*fl[14]+0.4330127018922193*alpha[3]*fl[13]+0.4330127018922193*alpha[7]*fl[12]+0.4330127018922193*alpha[1]*fl[11]+0.25*alpha[3]*fl[10]+0.25*alpha[7]*fl[9]+0.4330127018922193*alpha[6]*fl[8]+0.25*alpha[1]*fl[7]+0.4330127018922193*alpha[0]*fl[6]+0.25*fl[4]*alpha[6]+0.4330127018922193*alpha[4]*fl[5]+0.25*fl[2]*alpha[4]+0.25*alpha[0]*fl[3]+0.4330127018922193*fl[1]*alpha[2]+0.25*fl[0]*alpha[2]; 
  Ghat[3] = 0.4330127018922193*alpha[4]*fl[15]+0.25*alpha[4]*fl[14]+0.4330127018922193*alpha[2]*fl[13]+0.4330127018922193*alpha[1]*fl[12]+0.4330127018922193*alpha[7]*fl[11]+0.25*alpha[2]*fl[10]+0.25*alpha[1]*fl[9]+0.4330127018922193*alpha[0]*fl[8]+0.25*alpha[7]*fl[7]+0.4330127018922193*alpha[6]*fl[6]+0.25*fl[3]*alpha[6]+0.4330127018922193*alpha[5]*fl[5]+0.25*fl[2]*alpha[5]+0.25*alpha[0]*fl[4]+0.4330127018922193*fl[1]*alpha[3]+0.25*fl[0]*alpha[3]; 
  Ghat[4] = 0.4330127018922193*alpha[3]*fl[15]+0.25*alpha[3]*fl[14]+0.4330127018922193*alpha[5]*fl[13]+0.4330127018922193*alpha[6]*fl[12]+0.4330127018922193*alpha[0]*fl[11]+0.25*alpha[5]*fl[10]+0.25*alpha[6]*fl[9]+0.4330127018922193*alpha[7]*fl[8]+0.25*alpha[0]*fl[7]+0.25*fl[4]*alpha[7]+0.4330127018922193*alpha[1]*fl[6]+0.4330127018922193*alpha[2]*fl[5]+0.4330127018922193*fl[1]*alpha[4]+0.25*fl[0]*alpha[4]+0.25*alpha[1]*fl[3]+0.25*alpha[2]*fl[2]; 
  Ghat[5] = 0.4330127018922193*alpha[2]*fl[15]+0.25*alpha[2]*fl[14]+0.4330127018922193*alpha[4]*fl[13]+0.4330127018922193*alpha[0]*fl[12]+0.4330127018922193*alpha[6]*fl[11]+0.25*alpha[4]*fl[10]+0.25*alpha[0]*fl[9]+0.4330127018922193*alpha[1]*fl[8]+0.25*alpha[6]*fl[7]+0.4330127018922193*fl[6]*alpha[7]+0.25*fl[3]*alpha[7]+0.4330127018922193*alpha[3]*fl[5]+0.4330127018922193*fl[1]*alpha[5]+0.25*fl[0]*alpha[5]+0.25*alpha[1]*fl[4]+0.25*fl[2]*alpha[3]; 
  Ghat[6] = 0.4330127018922193*alpha[1]*fl[15]+0.25*alpha[1]*fl[14]+0.4330127018922193*alpha[0]*fl[13]+0.4330127018922193*alpha[4]*fl[12]+0.4330127018922193*alpha[5]*fl[11]+0.25*alpha[0]*fl[10]+0.25*alpha[4]*fl[9]+0.4330127018922193*alpha[2]*fl[8]+0.25*alpha[5]*fl[7]+0.4330127018922193*fl[5]*alpha[7]+0.25*fl[2]*alpha[7]+0.4330127018922193*alpha[3]*fl[6]+0.4330127018922193*fl[1]*alpha[6]+0.25*fl[0]*alpha[6]+0.25*alpha[2]*fl[4]+0.25*alpha[3]*fl[3]; 
  Ghat[7] = 0.4330127018922193*alpha[0]*fl[15]+0.25*alpha[0]*fl[14]+0.4330127018922193*alpha[1]*fl[13]+0.4330127018922193*alpha[2]*fl[12]+0.4330127018922193*alpha[3]*fl[11]+0.25*alpha[1]*fl[10]+0.25*alpha[2]*fl[9]+0.4330127018922193*alpha[4]*fl[8]+0.25*alpha[3]*fl[7]+0.4330127018922193*fl[1]*alpha[7]+0.25*fl[0]*alpha[7]+0.4330127018922193*alpha[5]*fl[6]+0.4330127018922193*fl[5]*alpha[6]+0.25*fl[2]*alpha[6]+0.25*fl[3]*alpha[5]+0.25*alpha[4]*fl[4]; 

  outr[0] += 0.7071067811865475*Ghat[0]*rdxr2; 
  outr[1] += -1.224744871391589*Ghat[0]*rdxr2; 
  outr[2] += 0.7071067811865475*Ghat[1]*rdxr2; 
  outr[3] += 0.7071067811865475*Ghat[2]*rdxr2; 
  outr[4] += 0.7071067811865475*Ghat[3]*rdxr2; 
  outr[5] += -1.224744871391589*Ghat[1]*rdxr2; 
  outr[6] += -1.224744871391589*Ghat[2]*rdxr2; 
  outr[7] += 0.7071067811865475*Ghat[4]*rdxr2; 
  outr[8] += -1.224744871391589*Ghat[3]*rdxr2; 
  outr[9] += 0.7071067811865475*Ghat[5]*rdxr2; 
  outr[10] += 0.7071067811865475*Ghat[6]*rdxr2; 
  outr[11] += -1.224744871391589*Ghat[4]*rdxr2; 
  outr[12] += -1.224744871391589*Ghat[5]*rdxr2; 
  outr[13] += -1.224744871391589*Ghat[6]*rdxr2; 
  outr[14] += 0.7071067811865475*Ghat[7]*rdxr2; 
  outr[15] += -1.224744871391589*Ghat[7]*rdxr2; 

  outl[0] += -0.7071067811865475*Ghat[0]*rdxl2; 
  outl[1] += -1.224744871391589*Ghat[0]*rdxl2; 
  outl[2] += -0.7071067811865475*Ghat[1]*rdxl2; 
  outl[3] += -0.7071067811865475*Ghat[2]*rdxl2; 
  outl[4] += -0.7071067811865475*Ghat[3]*rdxl2; 
  outl[5] += -1.224744871391589*Ghat[1]*rdxl2; 
  outl[6] += -1.224744871391589*Ghat[2]*rdxl2; 
  outl[7] += -0.7071067811865475*Ghat[4]*rdxl2; 
  outl[8] += -1.224744871391589*Ghat[3]*rdxl2; 
  outl[9] += -0.7071067811865475*Ghat[5]*rdxl2; 
  outl[10] += -0.7071067811865475*Ghat[6]*rdxl2; 
  outl[11] += -1.224744871391589*Ghat[4]*rdxl2; 
  outl[12] += -1.224744871391589*Ghat[5]*rdxl2; 
  outl[13] += -1.224744871391589*Ghat[6]*rdxl2; 
  outl[14] += -0.7071067811865475*Ghat[7]*rdxl2; 
  outl[15] += -1.224744871391589*Ghat[7]*rdxl2; 
  } else { 
  Ghat[0] = (-0.4330127018922193*alpha[7]*fr[15])+0.25*alpha[7]*fr[14]-0.4330127018922193*alpha[6]*fr[13]-0.4330127018922193*alpha[5]*fr[12]-0.4330127018922193*alpha[4]*fr[11]+0.25*alpha[6]*fr[10]+0.25*alpha[5]*fr[9]-0.4330127018922193*alpha[3]*fr[8]+0.25*alpha[4]*fr[7]-0.4330127018922193*alpha[2]*fr[6]-0.4330127018922193*alpha[1]*fr[5]+0.25*alpha[3]*fr[4]+0.25*alpha[2]*fr[3]+0.25*alpha[1]*fr[2]-0.4330127018922193*alpha[0]*fr[1]+0.25*alpha[0]*fr[0]; 
  Ghat[1] = (-0.4330127018922193*alpha[6]*fr[15])+0.25*alpha[6]*fr[14]-0.4330127018922193*alpha[7]*fr[13]-0.4330127018922193*alpha[3]*fr[12]-0.4330127018922193*alpha[2]*fr[11]+0.25*alpha[7]*fr[10]+0.25*alpha[3]*fr[9]-0.4330127018922193*alpha[5]*fr[8]+0.25*alpha[2]*fr[7]-0.4330127018922193*alpha[4]*fr[6]-0.4330127018922193*alpha[0]*fr[5]+0.25*fr[4]*alpha[5]+0.25*fr[3]*alpha[4]+0.25*alpha[0]*fr[2]-0.4330127018922193*alpha[1]*fr[1]+0.25*fr[0]*alpha[1]; 
  Ghat[2] = (-0.4330127018922193*alpha[5]*fr[15])+0.25*alpha[5]*fr[14]-0.4330127018922193*alpha[3]*fr[13]-0.4330127018922193*alpha[7]*fr[12]-0.4330127018922193*alpha[1]*fr[11]+0.25*alpha[3]*fr[10]+0.25*alpha[7]*fr[9]-0.4330127018922193*alpha[6]*fr[8]+0.25*alpha[1]*fr[7]-0.4330127018922193*alpha[0]*fr[6]+0.25*fr[4]*alpha[6]-0.4330127018922193*alpha[4]*fr[5]+0.25*fr[2]*alpha[4]+0.25*alpha[0]*fr[3]-0.4330127018922193*fr[1]*alpha[2]+0.25*fr[0]*alpha[2]; 
  Ghat[3] = (-0.4330127018922193*alpha[4]*fr[15])+0.25*alpha[4]*fr[14]-0.4330127018922193*alpha[2]*fr[13]-0.4330127018922193*alpha[1]*fr[12]-0.4330127018922193*alpha[7]*fr[11]+0.25*alpha[2]*fr[10]+0.25*alpha[1]*fr[9]-0.4330127018922193*alpha[0]*fr[8]+0.25*alpha[7]*fr[7]-0.4330127018922193*alpha[6]*fr[6]+0.25*fr[3]*alpha[6]-0.4330127018922193*alpha[5]*fr[5]+0.25*fr[2]*alpha[5]+0.25*alpha[0]*fr[4]-0.4330127018922193*fr[1]*alpha[3]+0.25*fr[0]*alpha[3]; 
  Ghat[4] = (-0.4330127018922193*alpha[3]*fr[15])+0.25*alpha[3]*fr[14]-0.4330127018922193*alpha[5]*fr[13]-0.4330127018922193*alpha[6]*fr[12]-0.4330127018922193*alpha[0]*fr[11]+0.25*alpha[5]*fr[10]+0.25*alpha[6]*fr[9]-0.4330127018922193*alpha[7]*fr[8]+0.25*alpha[0]*fr[7]+0.25*fr[4]*alpha[7]-0.4330127018922193*alpha[1]*fr[6]-0.4330127018922193*alpha[2]*fr[5]-0.4330127018922193*fr[1]*alpha[4]+0.25*fr[0]*alpha[4]+0.25*alpha[1]*fr[3]+0.25*alpha[2]*fr[2]; 
  Ghat[5] = (-0.4330127018922193*alpha[2]*fr[15])+0.25*alpha[2]*fr[14]-0.4330127018922193*alpha[4]*fr[13]-0.4330127018922193*alpha[0]*fr[12]-0.4330127018922193*alpha[6]*fr[11]+0.25*alpha[4]*fr[10]+0.25*alpha[0]*fr[9]-0.4330127018922193*alpha[1]*fr[8]+0.25*alpha[6]*fr[7]-0.4330127018922193*fr[6]*alpha[7]+0.25*fr[3]*alpha[7]-0.4330127018922193*alpha[3]*fr[5]-0.4330127018922193*fr[1]*alpha[5]+0.25*fr[0]*alpha[5]+0.25*alpha[1]*fr[4]+0.25*fr[2]*alpha[3]; 
  Ghat[6] = (-0.4330127018922193*alpha[1]*fr[15])+0.25*alpha[1]*fr[14]-0.4330127018922193*alpha[0]*fr[13]-0.4330127018922193*alpha[4]*fr[12]-0.4330127018922193*alpha[5]*fr[11]+0.25*alpha[0]*fr[10]+0.25*alpha[4]*fr[9]-0.4330127018922193*alpha[2]*fr[8]+0.25*alpha[5]*fr[7]-0.4330127018922193*fr[5]*alpha[7]+0.25*fr[2]*alpha[7]-0.4330127018922193*alpha[3]*fr[6]-0.4330127018922193*fr[1]*alpha[6]+0.25*fr[0]*alpha[6]+0.25*alpha[2]*fr[4]+0.25*alpha[3]*fr[3]; 
  Ghat[7] = (-0.4330127018922193*alpha[0]*fr[15])+0.25*alpha[0]*fr[14]-0.4330127018922193*alpha[1]*fr[13]-0.4330127018922193*alpha[2]*fr[12]-0.4330127018922193*alpha[3]*fr[11]+0.25*alpha[1]*fr[10]+0.25*alpha[2]*fr[9]-0.4330127018922193*alpha[4]*fr[8]+0.25*alpha[3]*fr[7]-0.4330127018922193*fr[1]*alpha[7]+0.25*fr[0]*alpha[7]-0.4330127018922193*alpha[5]*fr[6]-0.4330127018922193*fr[5]*alpha[6]+0.25*fr[2]*alpha[6]+0.25*fr[3]*alpha[5]+0.25*alpha[4]*fr[4]; 

  outr[0] += 0.7071067811865475*Ghat[0]*rdxr2; 
  outr[1] += -1.224744871391589*Ghat[0]*rdxr2; 
  outr[2] += 0.7071067811865475*Ghat[1]*rdxr2; 
  outr[3] += 0.7071067811865475*Ghat[2]*rdxr2; 
  outr[4] += 0.7071067811865475*Ghat[3]*rdxr2; 
  outr[5] += -1.224744871391589*Ghat[1]*rdxr2; 
  outr[6] += -1.224744871391589*Ghat[2]*rdxr2; 
  outr[7] += 0.7071067811865475*Ghat[4]*rdxr2; 
  outr[8] += -1.224744871391589*Ghat[3]*rdxr2; 
  outr[9] += 0.7071067811865475*Ghat[5]*rdxr2; 
  outr[10] += 0.7071067811865475*Ghat[6]*rdxr2; 
  outr[11] += -1.224744871391589*Ghat[4]*rdxr2; 
  outr[12] += -1.224744871391589*Ghat[5]*rdxr2; 
  outr[13] += -1.224744871391589*Ghat[6]*rdxr2; 
  outr[14] += 0.7071067811865475*Ghat[7]*rdxr2; 
  outr[15] += -1.224744871391589*Ghat[7]*rdxr2; 

  outl[0] += -0.7071067811865475*Ghat[0]*rdxl2; 
  outl[1] += -1.224744871391589*Ghat[0]*rdxl2; 
  outl[2] += -0.7071067811865475*Ghat[1]*rdxl2; 
  outl[3] += -0.7071067811865475*Ghat[2]*rdxl2; 
  outl[4] += -0.7071067811865475*Ghat[3]*rdxl2; 
  outl[5] += -1.224744871391589*Ghat[1]*rdxl2; 
  outl[6] += -1.224744871391589*Ghat[2]*rdxl2; 
  outl[7] += -0.7071067811865475*Ghat[4]*rdxl2; 
  outl[8] += -1.224744871391589*Ghat[3]*rdxl2; 
  outl[9] += -0.7071067811865475*Ghat[5]*rdxl2; 
  outl[10] += -0.7071067811865475*Ghat[6]*rdxl2; 
  outl[11] += -1.224744871391589*Ghat[4]*rdxl2; 
  outl[12] += -1.224744871391589*Ghat[5]*rdxl2; 
  outl[13] += -1.224744871391589*Ghat[6]*rdxl2; 
  outl[14] += -0.7071067811865475*Ghat[7]*rdxl2; 
  outl[15] += -1.224744871391589*Ghat[7]*rdxl2; 
  } 
} 
