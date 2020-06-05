#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag2x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[2]; 
  double dv10r = 2/dxvr[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 
  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double Ghat[15]; 
  double favg[15]; 
  double alpha[15]; 

  favg[0] = 1.58113883008419*fr[18]+1.58113883008419*fl[18]-1.224744871391589*fr[3]+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[7])+1.224744871391589*fl[7]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[8])+1.224744871391589*fl[8]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = 0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = 0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[7] = 0.7071067811865475*fr[10]+0.7071067811865475*fl[10]; 
  favg[8] = 0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = 0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 
  favg[10] = 0.7071067811865475*fr[15]+0.7071067811865475*fl[15]; 
  favg[11] = 0.7071067811865475*fr[16]+0.7071067811865475*fl[16]; 
  favg[12] = 0.7071067811865475*fr[17]+0.7071067811865475*fl[17]; 
  favg[13] = 0.7071067811865475*fr[19]+0.7071067811865475*fl[19]; 
  favg[14] = 0.7071067811865475*fr[20]+0.7071067811865475*fl[20]; 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 2.0*(B2[2]*wv2+E0[2])-2.0*B1[2]*wv3; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[4] = -0.5773502691896258*B1[0]*dv3; 
  alpha[5] = 2.0*(B2[3]*wv2+E0[3])-2.0*B1[3]*wv3; 
  alpha[6] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 0.5773502691896258*B2[2]*dv2; 
  alpha[8] = -0.5773502691896258*B1[1]*dv3; 
  alpha[9] = -0.5773502691896258*B1[2]*dv3; 
  alpha[11] = 2.0*(B2[4]*wv2+E0[4])-2.0*B1[4]*wv3; 
  alpha[12] = 2.0*(B2[5]*wv2+E0[5])-2.0*B1[5]*wv3; 

  const double amid = (-0.2795084971874737*alpha[12])-0.2795084971874737*alpha[11]+0.25*alpha[0]; 

  Ghat[0] = 0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[18]-1.0*(2.23606797749979*fl[18]+1.732050807568877*(fr[3]+fl[3]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[1]+fl[1])*amax+0.025*(4.47213595499958*(alpha[1]*favg[11]+favg[1]*alpha[11])+5.0*(alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[8]+fl[8])-1.0*fr[2]+fl[2])*amax+0.025*(4.47213595499958*(alpha[2]*favg[12]+favg[2]*alpha[12])+5.0*(alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[4]+fl[4])*amax+0.025*(4.47213595499958*alpha[3]*favg[13]+5.0*(alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3])); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[5]+fl[5])*amax+0.025*(4.47213595499958*alpha[4]*favg[14]+5.0*(alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4])); 
  Ghat[5] = 0.025*(4.47213595499958*(alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[5]*favg[11]+favg[5]*alpha[11])+5.0*(alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]))-0.3535533905932737*(fr[6]-1.0*fl[6])*amax; 
  Ghat[6] = 0.025*(4.47213595499958*(alpha[6]*(favg[13]+favg[11])+favg[6]*alpha[11])+5.0*(alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]))-0.3535533905932737*(fr[9]-1.0*fl[9])*amax; 
  Ghat[7] = 0.025*(4.47213595499958*(alpha[7]*(favg[13]+favg[12])+favg[7]*alpha[12])+5.0*(alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]))-0.3535533905932737*(fr[10]-1.0*fl[10])*amax; 
  Ghat[8] = 0.025*(4.47213595499958*(alpha[8]*(favg[14]+favg[11])+favg[8]*alpha[11])+5.0*(alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]))-0.3535533905932737*(fr[12]-1.0*fl[12])*amax; 
  Ghat[9] = 0.025*(4.47213595499958*(alpha[9]*(favg[14]+favg[12])+favg[9]*alpha[12])+5.0*(alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]))-0.3535533905932737*(fr[13]-1.0*fl[13])*amax; 
  Ghat[10] = 0.125*(alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4])-0.3535533905932737*(fr[15]-1.0*fl[15])*amax; 
  Ghat[11] = 0.003571428571428571*((22.3606797749979*alpha[11]+35.0*alpha[0])*favg[11]+35.0*favg[0]*alpha[11]+31.30495168499706*(alpha[8]*favg[8]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[1]*favg[1]))-0.3535533905932737*(fr[16]-1.0*fl[16])*amax; 
  Ghat[12] = 0.003571428571428571*((22.3606797749979*alpha[12]+35.0*alpha[0])*favg[12]+35.0*favg[0]*alpha[12]+31.30495168499706*(alpha[9]*favg[9]+alpha[7]*favg[7]+alpha[5]*favg[5]+alpha[2]*favg[2]))-0.3535533905932737*(fr[17]-1.0*fl[17])*amax; 
  Ghat[13] = 0.025*(5.0*alpha[0]*favg[13]+4.47213595499958*(alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[3]*favg[3]))-0.3535533905932737*(fr[19]-1.0*fl[19])*amax; 
  Ghat[14] = 0.025*(5.0*alpha[0]*favg[14]+4.47213595499958*(alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[4]*favg[4]))-0.3535533905932737*(fr[20]-1.0*fl[20])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[3] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv10r; 
  outr[6] += 0.7071067811865475*Ghat[5]*dv10r; 
  outr[7] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[8] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[9] += 0.7071067811865475*Ghat[6]*dv10r; 
  outr[10] += 0.7071067811865475*Ghat[7]*dv10r; 
  outr[11] += -1.224744871391589*Ghat[3]*dv10r; 
  outr[12] += 0.7071067811865475*Ghat[8]*dv10r; 
  outr[13] += 0.7071067811865475*Ghat[9]*dv10r; 
  outr[14] += -1.224744871391589*Ghat[4]*dv10r; 
  outr[15] += 0.7071067811865475*Ghat[10]*dv10r; 
  outr[16] += 0.7071067811865475*Ghat[11]*dv10r; 
  outr[17] += 0.7071067811865475*Ghat[12]*dv10r; 
  outr[18] += 1.58113883008419*Ghat[0]*dv10r; 
  outr[19] += 0.7071067811865475*Ghat[13]*dv10r; 
  outr[20] += 0.7071067811865475*Ghat[14]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[3] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv10l; 
  outl[6] += -0.7071067811865475*Ghat[5]*dv10l; 
  outl[7] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[8] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[9] += -0.7071067811865475*Ghat[6]*dv10l; 
  outl[10] += -0.7071067811865475*Ghat[7]*dv10l; 
  outl[11] += -1.224744871391589*Ghat[3]*dv10l; 
  outl[12] += -0.7071067811865475*Ghat[8]*dv10l; 
  outl[13] += -0.7071067811865475*Ghat[9]*dv10l; 
  outl[14] += -1.224744871391589*Ghat[4]*dv10l; 
  outl[15] += -0.7071067811865475*Ghat[10]*dv10l; 
  outl[16] += -0.7071067811865475*Ghat[11]*dv10l; 
  outl[17] += -0.7071067811865475*Ghat[12]*dv10l; 
  outl[18] += -1.58113883008419*Ghat[0]*dv10l; 
  outl[19] += -0.7071067811865475*Ghat[13]*dv10l; 
  outl[20] += -0.7071067811865475*Ghat[14]*dv10l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag2x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[3]; 
  double dv11r = 2/dxvr[3]; 
  const double *E1 = &EM[6]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 
  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double Ghat[15]; 
  double favg[15]; 
  double alpha[15]; 

  favg[0] = 1.58113883008419*fr[19]+1.58113883008419*fl[19]-1.224744871391589*fr[4]+1.224744871391589*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[9])+1.224744871391589*fl[9]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = 0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = 0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = 0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[8] = 0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = 0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 
  favg[10] = 0.7071067811865475*fr[14]+0.7071067811865475*fl[14]; 
  favg[11] = 0.7071067811865475*fr[16]+0.7071067811865475*fl[16]; 
  favg[12] = 0.7071067811865475*fr[17]+0.7071067811865475*fl[17]; 
  favg[13] = 0.7071067811865475*fr[18]+0.7071067811865475*fl[18]; 
  favg[14] = 0.7071067811865475*fr[20]+0.7071067811865475*fl[20]; 

  alpha[0] = 2.0*B0[0]*wv3-2.0*B2[0]*wv1+2.0*E1[0]; 
  alpha[1] = 2.0*B0[1]*wv3-2.0*B2[1]*wv1+2.0*E1[1]; 
  alpha[2] = 2.0*B0[2]*wv3-2.0*B2[2]*wv1+2.0*E1[2]; 
  alpha[3] = -0.5773502691896258*B2[0]*dv1; 
  alpha[4] = 0.5773502691896258*B0[0]*dv3; 
  alpha[5] = 2.0*B0[3]*wv3-2.0*B2[3]*wv1+2.0*E1[3]; 
  alpha[6] = -0.5773502691896258*B2[1]*dv1; 
  alpha[7] = -0.5773502691896258*B2[2]*dv1; 
  alpha[8] = 0.5773502691896258*B0[1]*dv3; 
  alpha[9] = 0.5773502691896258*B0[2]*dv3; 
  alpha[11] = 2.0*B0[4]*wv3-2.0*B2[4]*wv1+2.0*E1[4]; 
  alpha[12] = 2.0*B0[5]*wv3-2.0*B2[5]*wv1+2.0*E1[5]; 

  const double amid = (-0.2795084971874737*alpha[12])-0.2795084971874737*alpha[11]+0.25*alpha[0]; 

  Ghat[0] = 0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[19]-1.0*(2.23606797749979*fl[19]+1.732050807568877*(fr[4]+fl[4]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[9]+fl[9])-1.0*fr[1]+fl[1])*amax+0.025*(4.47213595499958*(alpha[1]*favg[11]+favg[1]*alpha[11])+5.0*(alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[10]+fl[10])-1.0*fr[2]+fl[2])*amax+0.025*(4.47213595499958*(alpha[2]*favg[12]+favg[2]*alpha[12])+5.0*(alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[3]+fl[3])*amax+0.025*(4.47213595499958*alpha[3]*favg[13]+5.0*(alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3])); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[5]+fl[5])*amax+0.025*(4.47213595499958*alpha[4]*favg[14]+5.0*(alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4])); 
  Ghat[5] = 0.025*(4.47213595499958*(alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[5]*favg[11]+favg[5]*alpha[11])+5.0*(alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]))-0.3535533905932737*(fr[6]-1.0*fl[6])*amax; 
  Ghat[6] = 0.025*(4.47213595499958*(alpha[6]*(favg[13]+favg[11])+favg[6]*alpha[11])+5.0*(alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]))-0.3535533905932737*(fr[7]-1.0*fl[7])*amax; 
  Ghat[7] = 0.025*(4.47213595499958*(alpha[7]*(favg[13]+favg[12])+favg[7]*alpha[12])+5.0*(alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]))-0.3535533905932737*(fr[8]-1.0*fl[8])*amax; 
  Ghat[8] = 0.025*(4.47213595499958*(alpha[8]*(favg[14]+favg[11])+favg[8]*alpha[11])+5.0*(alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]))-0.3535533905932737*(fr[12]-1.0*fl[12])*amax; 
  Ghat[9] = 0.025*(4.47213595499958*(alpha[9]*(favg[14]+favg[12])+favg[9]*alpha[12])+5.0*(alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]))-0.3535533905932737*(fr[13]-1.0*fl[13])*amax; 
  Ghat[10] = 0.125*(alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4])-0.3535533905932737*(fr[14]-1.0*fl[14])*amax; 
  Ghat[11] = 0.003571428571428571*((22.3606797749979*alpha[11]+35.0*alpha[0])*favg[11]+35.0*favg[0]*alpha[11]+31.30495168499706*(alpha[8]*favg[8]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[1]*favg[1]))-0.3535533905932737*(fr[16]-1.0*fl[16])*amax; 
  Ghat[12] = 0.003571428571428571*((22.3606797749979*alpha[12]+35.0*alpha[0])*favg[12]+35.0*favg[0]*alpha[12]+31.30495168499706*(alpha[9]*favg[9]+alpha[7]*favg[7]+alpha[5]*favg[5]+alpha[2]*favg[2]))-0.3535533905932737*(fr[17]-1.0*fl[17])*amax; 
  Ghat[13] = 0.025*(5.0*alpha[0]*favg[13]+4.47213595499958*(alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[3]*favg[3]))-0.3535533905932737*(fr[18]-1.0*fl[18])*amax; 
  Ghat[14] = 0.025*(5.0*alpha[0]*favg[14]+4.47213595499958*(alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[4]*favg[4]))-0.3535533905932737*(fr[20]-1.0*fl[20])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv11r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv11r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv11r; 
  outr[3] += 0.7071067811865475*Ghat[3]*dv11r; 
  outr[4] += -1.224744871391589*Ghat[0]*dv11r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv11r; 
  outr[6] += 0.7071067811865475*Ghat[5]*dv11r; 
  outr[7] += 0.7071067811865475*Ghat[6]*dv11r; 
  outr[8] += 0.7071067811865475*Ghat[7]*dv11r; 
  outr[9] += -1.224744871391589*Ghat[1]*dv11r; 
  outr[10] += -1.224744871391589*Ghat[2]*dv11r; 
  outr[11] += -1.224744871391589*Ghat[3]*dv11r; 
  outr[12] += 0.7071067811865475*Ghat[8]*dv11r; 
  outr[13] += 0.7071067811865475*Ghat[9]*dv11r; 
  outr[14] += 0.7071067811865475*Ghat[10]*dv11r; 
  outr[15] += -1.224744871391589*Ghat[4]*dv11r; 
  outr[16] += 0.7071067811865475*Ghat[11]*dv11r; 
  outr[17] += 0.7071067811865475*Ghat[12]*dv11r; 
  outr[18] += 0.7071067811865475*Ghat[13]*dv11r; 
  outr[19] += 1.58113883008419*Ghat[0]*dv11r; 
  outr[20] += 0.7071067811865475*Ghat[14]*dv11r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv11l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv11l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv11l; 
  outl[3] += -0.7071067811865475*Ghat[3]*dv11l; 
  outl[4] += -1.224744871391589*Ghat[0]*dv11l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv11l; 
  outl[6] += -0.7071067811865475*Ghat[5]*dv11l; 
  outl[7] += -0.7071067811865475*Ghat[6]*dv11l; 
  outl[8] += -0.7071067811865475*Ghat[7]*dv11l; 
  outl[9] += -1.224744871391589*Ghat[1]*dv11l; 
  outl[10] += -1.224744871391589*Ghat[2]*dv11l; 
  outl[11] += -1.224744871391589*Ghat[3]*dv11l; 
  outl[12] += -0.7071067811865475*Ghat[8]*dv11l; 
  outl[13] += -0.7071067811865475*Ghat[9]*dv11l; 
  outl[14] += -0.7071067811865475*Ghat[10]*dv11l; 
  outl[15] += -1.224744871391589*Ghat[4]*dv11l; 
  outl[16] += -0.7071067811865475*Ghat[11]*dv11l; 
  outl[17] += -0.7071067811865475*Ghat[12]*dv11l; 
  outl[18] += -0.7071067811865475*Ghat[13]*dv11l; 
  outl[19] += -1.58113883008419*Ghat[0]*dv11l; 
  outl[20] += -0.7071067811865475*Ghat[14]*dv11l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag2x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12l = 2/dxvl[4]; 
  double dv12r = 2/dxvr[4]; 
  const double *E2 = &EM[12]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 
  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double Ghat[15]; 
  double favg[15]; 
  double alpha[15]; 

  favg[0] = 1.58113883008419*fr[20]+1.58113883008419*fl[20]-1.224744871391589*fr[5]+1.224744871391589*fl[5]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[12])+1.224744871391589*fl[12]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[13])+1.224744871391589*fl[13]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[5] = 0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = 0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = 0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[8] = 0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[9] = 0.7071067811865475*fr[10]+0.7071067811865475*fl[10]; 
  favg[10] = 0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 
  favg[11] = 0.7071067811865475*fr[16]+0.7071067811865475*fl[16]; 
  favg[12] = 0.7071067811865475*fr[17]+0.7071067811865475*fl[17]; 
  favg[13] = 0.7071067811865475*fr[18]+0.7071067811865475*fl[18]; 
  favg[14] = 0.7071067811865475*fr[19]+0.7071067811865475*fl[19]; 

  alpha[0] = 2.0*(B1[0]*wv1+E2[0])-2.0*B0[0]*wv2; 
  alpha[1] = 2.0*(B1[1]*wv1+E2[1])-2.0*B0[1]*wv2; 
  alpha[2] = 2.0*(B1[2]*wv1+E2[2])-2.0*B0[2]*wv2; 
  alpha[3] = 0.5773502691896258*B1[0]*dv1; 
  alpha[4] = -0.5773502691896258*B0[0]*dv2; 
  alpha[5] = 2.0*(B1[3]*wv1+E2[3])-2.0*B0[3]*wv2; 
  alpha[6] = 0.5773502691896258*B1[1]*dv1; 
  alpha[7] = 0.5773502691896258*B1[2]*dv1; 
  alpha[8] = -0.5773502691896258*B0[1]*dv2; 
  alpha[9] = -0.5773502691896258*B0[2]*dv2; 
  alpha[11] = 2.0*(B1[4]*wv1+E2[4])-2.0*B0[4]*wv2; 
  alpha[12] = 2.0*(B1[5]*wv1+E2[5])-2.0*B0[5]*wv2; 

  const double amid = (-0.2795084971874737*alpha[12])-0.2795084971874737*alpha[11]+0.25*alpha[0]; 

  Ghat[0] = 0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[20]-1.0*(2.23606797749979*fl[20]+1.732050807568877*(fr[5]+fl[5]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[12]+fl[12])-1.0*fr[1]+fl[1])*amax+0.025*(4.47213595499958*(alpha[1]*favg[11]+favg[1]*alpha[11])+5.0*(alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[13]+fl[13])-1.0*fr[2]+fl[2])*amax+0.025*(4.47213595499958*(alpha[2]*favg[12]+favg[2]*alpha[12])+5.0*(alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[3]+fl[3])*amax+0.025*(4.47213595499958*alpha[3]*favg[13]+5.0*(alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3])); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[4]+fl[4])*amax+0.025*(4.47213595499958*alpha[4]*favg[14]+5.0*(alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4])); 
  Ghat[5] = 0.025*(4.47213595499958*(alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[5]*favg[11]+favg[5]*alpha[11])+5.0*(alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]))-0.3535533905932737*(fr[6]-1.0*fl[6])*amax; 
  Ghat[6] = 0.025*(4.47213595499958*(alpha[6]*(favg[13]+favg[11])+favg[6]*alpha[11])+5.0*(alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]))-0.3535533905932737*(fr[7]-1.0*fl[7])*amax; 
  Ghat[7] = 0.025*(4.47213595499958*(alpha[7]*(favg[13]+favg[12])+favg[7]*alpha[12])+5.0*(alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]))-0.3535533905932737*(fr[8]-1.0*fl[8])*amax; 
  Ghat[8] = 0.025*(4.47213595499958*(alpha[8]*(favg[14]+favg[11])+favg[8]*alpha[11])+5.0*(alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]))-0.3535533905932737*(fr[9]-1.0*fl[9])*amax; 
  Ghat[9] = 0.025*(4.47213595499958*(alpha[9]*(favg[14]+favg[12])+favg[9]*alpha[12])+5.0*(alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]))-0.3535533905932737*(fr[10]-1.0*fl[10])*amax; 
  Ghat[10] = 0.125*(alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4])-0.3535533905932737*(fr[11]-1.0*fl[11])*amax; 
  Ghat[11] = 0.003571428571428571*((22.3606797749979*alpha[11]+35.0*alpha[0])*favg[11]+35.0*favg[0]*alpha[11]+31.30495168499706*(alpha[8]*favg[8]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[1]*favg[1]))-0.3535533905932737*(fr[16]-1.0*fl[16])*amax; 
  Ghat[12] = 0.003571428571428571*((22.3606797749979*alpha[12]+35.0*alpha[0])*favg[12]+35.0*favg[0]*alpha[12]+31.30495168499706*(alpha[9]*favg[9]+alpha[7]*favg[7]+alpha[5]*favg[5]+alpha[2]*favg[2]))-0.3535533905932737*(fr[17]-1.0*fl[17])*amax; 
  Ghat[13] = 0.025*(5.0*alpha[0]*favg[13]+4.47213595499958*(alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[3]*favg[3]))-0.3535533905932737*(fr[18]-1.0*fl[18])*amax; 
  Ghat[14] = 0.025*(5.0*alpha[0]*favg[14]+4.47213595499958*(alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[4]*favg[4]))-0.3535533905932737*(fr[19]-1.0*fl[19])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv12r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv12r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv12r; 
  outr[3] += 0.7071067811865475*Ghat[3]*dv12r; 
  outr[4] += 0.7071067811865475*Ghat[4]*dv12r; 
  outr[5] += -1.224744871391589*Ghat[0]*dv12r; 
  outr[6] += 0.7071067811865475*Ghat[5]*dv12r; 
  outr[7] += 0.7071067811865475*Ghat[6]*dv12r; 
  outr[8] += 0.7071067811865475*Ghat[7]*dv12r; 
  outr[9] += 0.7071067811865475*Ghat[8]*dv12r; 
  outr[10] += 0.7071067811865475*Ghat[9]*dv12r; 
  outr[11] += 0.7071067811865475*Ghat[10]*dv12r; 
  outr[12] += -1.224744871391589*Ghat[1]*dv12r; 
  outr[13] += -1.224744871391589*Ghat[2]*dv12r; 
  outr[14] += -1.224744871391589*Ghat[3]*dv12r; 
  outr[15] += -1.224744871391589*Ghat[4]*dv12r; 
  outr[16] += 0.7071067811865475*Ghat[11]*dv12r; 
  outr[17] += 0.7071067811865475*Ghat[12]*dv12r; 
  outr[18] += 0.7071067811865475*Ghat[13]*dv12r; 
  outr[19] += 0.7071067811865475*Ghat[14]*dv12r; 
  outr[20] += 1.58113883008419*Ghat[0]*dv12r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv12l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv12l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv12l; 
  outl[3] += -0.7071067811865475*Ghat[3]*dv12l; 
  outl[4] += -0.7071067811865475*Ghat[4]*dv12l; 
  outl[5] += -1.224744871391589*Ghat[0]*dv12l; 
  outl[6] += -0.7071067811865475*Ghat[5]*dv12l; 
  outl[7] += -0.7071067811865475*Ghat[6]*dv12l; 
  outl[8] += -0.7071067811865475*Ghat[7]*dv12l; 
  outl[9] += -0.7071067811865475*Ghat[8]*dv12l; 
  outl[10] += -0.7071067811865475*Ghat[9]*dv12l; 
  outl[11] += -0.7071067811865475*Ghat[10]*dv12l; 
  outl[12] += -1.224744871391589*Ghat[1]*dv12l; 
  outl[13] += -1.224744871391589*Ghat[2]*dv12l; 
  outl[14] += -1.224744871391589*Ghat[3]*dv12l; 
  outl[15] += -1.224744871391589*Ghat[4]*dv12l; 
  outl[16] += -0.7071067811865475*Ghat[11]*dv12l; 
  outl[17] += -0.7071067811865475*Ghat[12]*dv12l; 
  outl[18] += -0.7071067811865475*Ghat[13]*dv12l; 
  outl[19] += -0.7071067811865475*Ghat[14]*dv12l; 
  outl[20] += -1.58113883008419*Ghat[0]*dv12l; 

  return std::abs(amid); 
} 
