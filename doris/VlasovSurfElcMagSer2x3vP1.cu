//#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag2x3vSer_VX_P1(const double* __restrict__ wl, const double* __restrict__ wr, const double* __restrict__ dxvl, const double* __restrict__ dxvr, const double amax, const double* __restrict__ EM, const double* __restrict__ fl, const double* __restrict__ fr, double *outl, double *outr) 
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
  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double Ghat[16]; 
  double favg[16]; 
  double alpha[16]; 

  favg[0] = (-1.224744871391589*fr[3])+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[7])+1.224744871391589*fl[7]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[8])+1.224744871391589*fl[8]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[16])+1.224744871391589*fl[16]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[18])+1.224744871391589*fl[18]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[7] = (-1.224744871391589*fr[19])+1.224744871391589*fl[19]+0.7071067811865475*fr[10]+0.7071067811865475*fl[10]; 
  favg[8] = (-1.224744871391589*fr[21])+1.224744871391589*fl[21]+0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = (-1.224744871391589*fr[22])+1.224744871391589*fl[22]+0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 
  favg[10] = (-1.224744871391589*fr[25])+1.224744871391589*fl[25]+0.7071067811865475*fr[15]+0.7071067811865475*fl[15]; 
  favg[11] = (-1.224744871391589*fr[26])+1.224744871391589*fl[26]+0.7071067811865475*fr[17]+0.7071067811865475*fl[17]; 
  favg[12] = (-1.224744871391589*fr[27])+1.224744871391589*fl[27]+0.7071067811865475*fr[20]+0.7071067811865475*fl[20]; 
  favg[13] = (-1.224744871391589*fr[29])+1.224744871391589*fl[29]+0.7071067811865475*fr[23]+0.7071067811865475*fl[23]; 
  favg[14] = (-1.224744871391589*fr[30])+1.224744871391589*fl[30]+0.7071067811865475*fr[24]+0.7071067811865475*fl[24]; 
  favg[15] = (-1.224744871391589*fr[31])+1.224744871391589*fl[31]+0.7071067811865475*fr[28]+0.7071067811865475*fl[28]; 

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
  alpha[11] = 0.5773502691896258*B2[3]*dv2; 
  alpha[12] = -0.5773502691896258*B1[3]*dv3; 

  const double amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[3]+fl[3])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[1]+fl[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[8]+fl[8])-1.0*fr[2]+fl[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[4]+fl[4])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[5]+fl[5])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[16]+fl[16])-1.0*fr[6]+fl[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[18]+fl[18])-1.0*fr[9]+fl[9])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[19]+fl[19])-1.0*fr[10]+fl[10])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = 0.3535533905932737*(1.732050807568877*(fr[21]+fl[21])-1.0*fr[12]+fl[12])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = 0.3535533905932737*(1.732050807568877*(fr[22]+fl[22])-1.0*fr[13]+fl[13])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = 0.3535533905932737*(1.732050807568877*(fr[25]+fl[25])-1.0*fr[15]+fl[15])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = 0.3535533905932737*(1.732050807568877*(fr[26]+fl[26])-1.0*fr[17]+fl[17])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = 0.3535533905932737*(1.732050807568877*(fr[27]+fl[27])-1.0*fr[20]+fl[20])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = 0.3535533905932737*(1.732050807568877*(fr[29]+fl[29])-1.0*fr[23]+fl[23])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = 0.3535533905932737*(1.732050807568877*(fr[30]+fl[30])-1.0*fr[24]+fl[24])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = 0.3535533905932737*(1.732050807568877*(fr[31]+fl[31])-1.0*fr[28]+fl[28])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

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
  outr[16] += -1.224744871391589*Ghat[5]*dv10r; 
  outr[17] += 0.7071067811865475*Ghat[11]*dv10r; 
  outr[18] += -1.224744871391589*Ghat[6]*dv10r; 
  outr[19] += -1.224744871391589*Ghat[7]*dv10r; 
  outr[20] += 0.7071067811865475*Ghat[12]*dv10r; 
  outr[21] += -1.224744871391589*Ghat[8]*dv10r; 
  outr[22] += -1.224744871391589*Ghat[9]*dv10r; 
  outr[23] += 0.7071067811865475*Ghat[13]*dv10r; 
  outr[24] += 0.7071067811865475*Ghat[14]*dv10r; 
  outr[25] += -1.224744871391589*Ghat[10]*dv10r; 
  outr[26] += -1.224744871391589*Ghat[11]*dv10r; 
  outr[27] += -1.224744871391589*Ghat[12]*dv10r; 
  outr[28] += 0.7071067811865475*Ghat[15]*dv10r; 
  outr[29] += -1.224744871391589*Ghat[13]*dv10r; 
  outr[30] += -1.224744871391589*Ghat[14]*dv10r; 
  outr[31] += -1.224744871391589*Ghat[15]*dv10r; 

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
  outl[16] += -1.224744871391589*Ghat[5]*dv10l; 
  outl[17] += -0.7071067811865475*Ghat[11]*dv10l; 
  outl[18] += -1.224744871391589*Ghat[6]*dv10l; 
  outl[19] += -1.224744871391589*Ghat[7]*dv10l; 
  outl[20] += -0.7071067811865475*Ghat[12]*dv10l; 
  outl[21] += -1.224744871391589*Ghat[8]*dv10l; 
  outl[22] += -1.224744871391589*Ghat[9]*dv10l; 
  outl[23] += -0.7071067811865475*Ghat[13]*dv10l; 
  outl[24] += -0.7071067811865475*Ghat[14]*dv10l; 
  outl[25] += -1.224744871391589*Ghat[10]*dv10l; 
  outl[26] += -1.224744871391589*Ghat[11]*dv10l; 
  outl[27] += -1.224744871391589*Ghat[12]*dv10l; 
  outl[28] += -0.7071067811865475*Ghat[15]*dv10l; 
  outl[29] += -1.224744871391589*Ghat[13]*dv10l; 
  outl[30] += -1.224744871391589*Ghat[14]*dv10l; 
  outl[31] += -1.224744871391589*Ghat[15]*dv10l; 

  return std::abs(amid); 
} 
