#include <VlasovModDecl.h> 

__host__ __device__ double VlasovPhiSurf2x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double qDm, const double *phi, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qDm:       Species charge (q) divided by its mass (m).
  // phi:       electrostatic potential.
  // EM:        external EM field vectors.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double rdx2qDm = 2.*qDm/dxvl[0]; 
  const double rdy2qDm = 2.*qDm/dxvl[1]; 
  double dv10l = 2./dxvl[2]; 
  double dv10r = 2./dxvr[2]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *E0 = &EM[0]; 

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

  alpha[0] = 2.0*E0[0]-3.464101615137754*phi[1]*rdx2qDm; 
  alpha[1] = 2.0*E0[1]; 
  alpha[2] = 2.0*E0[2]-3.464101615137754*phi[3]*rdx2qDm; 
  alpha[5] = 2.0*E0[3]; 

  const double amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[3]+fl[3])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[5]*favg[5]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[7]+fl[7])-1.0*fr[1]+fl[1])*amax+0.125*(alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[8]+fl[8])-1.0*fr[2]+fl[2])*amax+0.125*(alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[4]+fl[4])*amax+0.125*(alpha[5]*favg[11]+alpha[2]*favg[7]+alpha[1]*favg[6]+alpha[0]*favg[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[5]+fl[5])*amax+0.125*(alpha[5]*favg[12]+alpha[2]*favg[9]+alpha[1]*favg[8]+alpha[0]*favg[4]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[16]+fl[16])-1.0*fr[6]+fl[6])*amax+0.125*(alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[18]+fl[18])-1.0*fr[9]+fl[9])*amax+0.125*(alpha[2]*favg[11]+alpha[5]*favg[7]+alpha[0]*favg[6]+alpha[1]*favg[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[19]+fl[19])-1.0*fr[10]+fl[10])*amax+0.125*(alpha[1]*favg[11]+alpha[0]*favg[7]+alpha[5]*favg[6]+alpha[2]*favg[3]); 
  Ghat[8] = 0.3535533905932737*(1.732050807568877*(fr[21]+fl[21])-1.0*fr[12]+fl[12])*amax+0.125*(alpha[2]*favg[12]+alpha[5]*favg[9]+alpha[0]*favg[8]+alpha[1]*favg[4]); 
  Ghat[9] = 0.3535533905932737*(1.732050807568877*(fr[22]+fl[22])-1.0*fr[13]+fl[13])*amax+0.125*(alpha[1]*favg[12]+alpha[0]*favg[9]+alpha[5]*favg[8]+alpha[2]*favg[4]); 
  Ghat[10] = 0.3535533905932737*(1.732050807568877*(fr[25]+fl[25])-1.0*fr[15]+fl[15])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[0]*favg[10]); 
  Ghat[11] = 0.3535533905932737*(1.732050807568877*(fr[26]+fl[26])-1.0*fr[17]+fl[17])*amax+0.125*(alpha[0]*favg[11]+alpha[1]*favg[7]+alpha[2]*favg[6]+favg[3]*alpha[5]); 
  Ghat[12] = 0.3535533905932737*(1.732050807568877*(fr[27]+fl[27])-1.0*fr[20]+fl[20])*amax+0.125*(alpha[0]*favg[12]+alpha[1]*favg[9]+alpha[2]*favg[8]+favg[4]*alpha[5]); 
  Ghat[13] = 0.3535533905932737*(1.732050807568877*(fr[29]+fl[29])-1.0*fr[23]+fl[23])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[1]*favg[10]); 
  Ghat[14] = 0.3535533905932737*(1.732050807568877*(fr[30]+fl[30])-1.0*fr[24]+fl[24])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[2]*favg[10]); 
  Ghat[15] = 0.3535533905932737*(1.732050807568877*(fr[31]+fl[31])-1.0*fr[28]+fl[28])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[5]*favg[10]); 

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

__host__ __device__ double VlasovPhiSurf2x3vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double qDm, const double *phi, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qDm:       Species charge (q) divided by its mass (m).
  // phi:       electrostatic potential.
  // EM:        external EM field vectors.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double rdx2qDm = 2.*qDm/dxvl[0]; 
  const double rdy2qDm = 2.*qDm/dxvl[1]; 
  double dv11l = 2./dxvl[3]; 
  double dv11r = 2./dxvr[3]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *E1 = &EM[4]; 

  double Ghat[16]; 
  double favg[16]; 
  double alpha[16]; 

  favg[0] = (-1.224744871391589*fr[4])+1.224744871391589*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[9])+1.224744871391589*fl[9]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[17])+1.224744871391589*fl[17]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[18])+1.224744871391589*fl[18]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = (-1.224744871391589*fr[19])+1.224744871391589*fl[19]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[8] = (-1.224744871391589*fr[23])+1.224744871391589*fl[23]+0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = (-1.224744871391589*fr[24])+1.224744871391589*fl[24]+0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 
  favg[10] = (-1.224744871391589*fr[25])+1.224744871391589*fl[25]+0.7071067811865475*fr[14]+0.7071067811865475*fl[14]; 
  favg[11] = (-1.224744871391589*fr[26])+1.224744871391589*fl[26]+0.7071067811865475*fr[16]+0.7071067811865475*fl[16]; 
  favg[12] = (-1.224744871391589*fr[28])+1.224744871391589*fl[28]+0.7071067811865475*fr[20]+0.7071067811865475*fl[20]; 
  favg[13] = (-1.224744871391589*fr[29])+1.224744871391589*fl[29]+0.7071067811865475*fr[21]+0.7071067811865475*fl[21]; 
  favg[14] = (-1.224744871391589*fr[30])+1.224744871391589*fl[30]+0.7071067811865475*fr[22]+0.7071067811865475*fl[22]; 
  favg[15] = (-1.224744871391589*fr[31])+1.224744871391589*fl[31]+0.7071067811865475*fr[27]+0.7071067811865475*fl[27]; 

  alpha[0] = 2.0*E1[0]-3.464101615137754*phi[2]*rdy2qDm; 
  alpha[1] = 2.0*E1[1]-3.464101615137754*phi[3]*rdy2qDm; 
  alpha[2] = 2.0*E1[2]; 
  alpha[5] = 2.0*E1[3]; 

  const double amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[4]+fl[4])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[5]*favg[5]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[9]+fl[9])-1.0*fr[1]+fl[1])*amax+0.125*(alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[10]+fl[10])-1.0*fr[2]+fl[2])*amax+0.125*(alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[3]+fl[3])*amax+0.125*(alpha[5]*favg[11]+alpha[2]*favg[7]+alpha[1]*favg[6]+alpha[0]*favg[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[5]+fl[5])*amax+0.125*(alpha[5]*favg[12]+alpha[2]*favg[9]+alpha[1]*favg[8]+alpha[0]*favg[4]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[17]+fl[17])-1.0*fr[6]+fl[6])*amax+0.125*(alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[18]+fl[18])-1.0*fr[7]+fl[7])*amax+0.125*(alpha[2]*favg[11]+alpha[5]*favg[7]+alpha[0]*favg[6]+alpha[1]*favg[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[19]+fl[19])-1.0*fr[8]+fl[8])*amax+0.125*(alpha[1]*favg[11]+alpha[0]*favg[7]+alpha[5]*favg[6]+alpha[2]*favg[3]); 
  Ghat[8] = 0.3535533905932737*(1.732050807568877*(fr[23]+fl[23])-1.0*fr[12]+fl[12])*amax+0.125*(alpha[2]*favg[12]+alpha[5]*favg[9]+alpha[0]*favg[8]+alpha[1]*favg[4]); 
  Ghat[9] = 0.3535533905932737*(1.732050807568877*(fr[24]+fl[24])-1.0*fr[13]+fl[13])*amax+0.125*(alpha[1]*favg[12]+alpha[0]*favg[9]+alpha[5]*favg[8]+alpha[2]*favg[4]); 
  Ghat[10] = 0.3535533905932737*(1.732050807568877*(fr[25]+fl[25])-1.0*fr[14]+fl[14])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[0]*favg[10]); 
  Ghat[11] = 0.3535533905932737*(1.732050807568877*(fr[26]+fl[26])-1.0*fr[16]+fl[16])*amax+0.125*(alpha[0]*favg[11]+alpha[1]*favg[7]+alpha[2]*favg[6]+favg[3]*alpha[5]); 
  Ghat[12] = 0.3535533905932737*(1.732050807568877*(fr[28]+fl[28])-1.0*fr[20]+fl[20])*amax+0.125*(alpha[0]*favg[12]+alpha[1]*favg[9]+alpha[2]*favg[8]+favg[4]*alpha[5]); 
  Ghat[13] = 0.3535533905932737*(1.732050807568877*(fr[29]+fl[29])-1.0*fr[21]+fl[21])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[1]*favg[10]); 
  Ghat[14] = 0.3535533905932737*(1.732050807568877*(fr[30]+fl[30])-1.0*fr[22]+fl[22])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[2]*favg[10]); 
  Ghat[15] = 0.3535533905932737*(1.732050807568877*(fr[31]+fl[31])-1.0*fr[27]+fl[27])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[5]*favg[10]); 

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
  outr[17] += -1.224744871391589*Ghat[5]*dv11r; 
  outr[18] += -1.224744871391589*Ghat[6]*dv11r; 
  outr[19] += -1.224744871391589*Ghat[7]*dv11r; 
  outr[20] += 0.7071067811865475*Ghat[12]*dv11r; 
  outr[21] += 0.7071067811865475*Ghat[13]*dv11r; 
  outr[22] += 0.7071067811865475*Ghat[14]*dv11r; 
  outr[23] += -1.224744871391589*Ghat[8]*dv11r; 
  outr[24] += -1.224744871391589*Ghat[9]*dv11r; 
  outr[25] += -1.224744871391589*Ghat[10]*dv11r; 
  outr[26] += -1.224744871391589*Ghat[11]*dv11r; 
  outr[27] += 0.7071067811865475*Ghat[15]*dv11r; 
  outr[28] += -1.224744871391589*Ghat[12]*dv11r; 
  outr[29] += -1.224744871391589*Ghat[13]*dv11r; 
  outr[30] += -1.224744871391589*Ghat[14]*dv11r; 
  outr[31] += -1.224744871391589*Ghat[15]*dv11r; 

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
  outl[17] += -1.224744871391589*Ghat[5]*dv11l; 
  outl[18] += -1.224744871391589*Ghat[6]*dv11l; 
  outl[19] += -1.224744871391589*Ghat[7]*dv11l; 
  outl[20] += -0.7071067811865475*Ghat[12]*dv11l; 
  outl[21] += -0.7071067811865475*Ghat[13]*dv11l; 
  outl[22] += -0.7071067811865475*Ghat[14]*dv11l; 
  outl[23] += -1.224744871391589*Ghat[8]*dv11l; 
  outl[24] += -1.224744871391589*Ghat[9]*dv11l; 
  outl[25] += -1.224744871391589*Ghat[10]*dv11l; 
  outl[26] += -1.224744871391589*Ghat[11]*dv11l; 
  outl[27] += -0.7071067811865475*Ghat[15]*dv11l; 
  outl[28] += -1.224744871391589*Ghat[12]*dv11l; 
  outl[29] += -1.224744871391589*Ghat[13]*dv11l; 
  outl[30] += -1.224744871391589*Ghat[14]*dv11l; 
  outl[31] += -1.224744871391589*Ghat[15]*dv11l; 

  return std::abs(amid); 
} 

__host__ __device__ double VlasovPhiSurf2x3vSer_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double qDm, const double *phi, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qDm:       Species charge (q) divided by its mass (m).
  // phi:       electrostatic potential.
  // EM:        external EM field vectors.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double rdx2qDm = 2.*qDm/dxvl[0]; 
  const double rdy2qDm = 2.*qDm/dxvl[1]; 
  double dv12l = 2./dxvl[4]; 
  double dv12r = 2./dxvr[4]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *E2 = &EM[8]; 

  double Ghat[16]; 
  double favg[16]; 
  double alpha[16]; 

  favg[0] = (-1.224744871391589*fr[5])+1.224744871391589*fl[5]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[12])+1.224744871391589*fl[12]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[13])+1.224744871391589*fl[13]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[5] = (-1.224744871391589*fr[20])+1.224744871391589*fl[20]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[21])+1.224744871391589*fl[21]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = (-1.224744871391589*fr[22])+1.224744871391589*fl[22]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[8] = (-1.224744871391589*fr[23])+1.224744871391589*fl[23]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[9] = (-1.224744871391589*fr[24])+1.224744871391589*fl[24]+0.7071067811865475*fr[10]+0.7071067811865475*fl[10]; 
  favg[10] = (-1.224744871391589*fr[25])+1.224744871391589*fl[25]+0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 
  favg[11] = (-1.224744871391589*fr[27])+1.224744871391589*fl[27]+0.7071067811865475*fr[16]+0.7071067811865475*fl[16]; 
  favg[12] = (-1.224744871391589*fr[28])+1.224744871391589*fl[28]+0.7071067811865475*fr[17]+0.7071067811865475*fl[17]; 
  favg[13] = (-1.224744871391589*fr[29])+1.224744871391589*fl[29]+0.7071067811865475*fr[18]+0.7071067811865475*fl[18]; 
  favg[14] = (-1.224744871391589*fr[30])+1.224744871391589*fl[30]+0.7071067811865475*fr[19]+0.7071067811865475*fl[19]; 
  favg[15] = (-1.224744871391589*fr[31])+1.224744871391589*fl[31]+0.7071067811865475*fr[26]+0.7071067811865475*fl[26]; 

  alpha[0] = 2.0*E2[0]; 
  alpha[1] = 2.0*E2[1]; 
  alpha[2] = 2.0*E2[2]; 
  alpha[5] = 2.0*E2[3]; 

  const double amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[5]+fl[5])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[5]*favg[5]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[12]+fl[12])-1.0*fr[1]+fl[1])*amax+0.125*(alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[13]+fl[13])-1.0*fr[2]+fl[2])*amax+0.125*(alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[3]+fl[3])*amax+0.125*(alpha[5]*favg[11]+alpha[2]*favg[7]+alpha[1]*favg[6]+alpha[0]*favg[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[4]+fl[4])*amax+0.125*(alpha[5]*favg[12]+alpha[2]*favg[9]+alpha[1]*favg[8]+alpha[0]*favg[4]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[20]+fl[20])-1.0*fr[6]+fl[6])*amax+0.125*(alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[21]+fl[21])-1.0*fr[7]+fl[7])*amax+0.125*(alpha[2]*favg[11]+alpha[5]*favg[7]+alpha[0]*favg[6]+alpha[1]*favg[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[22]+fl[22])-1.0*fr[8]+fl[8])*amax+0.125*(alpha[1]*favg[11]+alpha[0]*favg[7]+alpha[5]*favg[6]+alpha[2]*favg[3]); 
  Ghat[8] = 0.3535533905932737*(1.732050807568877*(fr[23]+fl[23])-1.0*fr[9]+fl[9])*amax+0.125*(alpha[2]*favg[12]+alpha[5]*favg[9]+alpha[0]*favg[8]+alpha[1]*favg[4]); 
  Ghat[9] = 0.3535533905932737*(1.732050807568877*(fr[24]+fl[24])-1.0*fr[10]+fl[10])*amax+0.125*(alpha[1]*favg[12]+alpha[0]*favg[9]+alpha[5]*favg[8]+alpha[2]*favg[4]); 
  Ghat[10] = 0.3535533905932737*(1.732050807568877*(fr[25]+fl[25])-1.0*fr[11]+fl[11])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[0]*favg[10]); 
  Ghat[11] = 0.3535533905932737*(1.732050807568877*(fr[27]+fl[27])-1.0*fr[16]+fl[16])*amax+0.125*(alpha[0]*favg[11]+alpha[1]*favg[7]+alpha[2]*favg[6]+favg[3]*alpha[5]); 
  Ghat[12] = 0.3535533905932737*(1.732050807568877*(fr[28]+fl[28])-1.0*fr[17]+fl[17])*amax+0.125*(alpha[0]*favg[12]+alpha[1]*favg[9]+alpha[2]*favg[8]+favg[4]*alpha[5]); 
  Ghat[13] = 0.3535533905932737*(1.732050807568877*(fr[29]+fl[29])-1.0*fr[18]+fl[18])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[1]*favg[10]); 
  Ghat[14] = 0.3535533905932737*(1.732050807568877*(fr[30]+fl[30])-1.0*fr[19]+fl[19])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[2]*favg[10]); 
  Ghat[15] = 0.3535533905932737*(1.732050807568877*(fr[31]+fl[31])-1.0*fr[26]+fl[26])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[5]*favg[10]); 

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
  outr[20] += -1.224744871391589*Ghat[5]*dv12r; 
  outr[21] += -1.224744871391589*Ghat[6]*dv12r; 
  outr[22] += -1.224744871391589*Ghat[7]*dv12r; 
  outr[23] += -1.224744871391589*Ghat[8]*dv12r; 
  outr[24] += -1.224744871391589*Ghat[9]*dv12r; 
  outr[25] += -1.224744871391589*Ghat[10]*dv12r; 
  outr[26] += 0.7071067811865475*Ghat[15]*dv12r; 
  outr[27] += -1.224744871391589*Ghat[11]*dv12r; 
  outr[28] += -1.224744871391589*Ghat[12]*dv12r; 
  outr[29] += -1.224744871391589*Ghat[13]*dv12r; 
  outr[30] += -1.224744871391589*Ghat[14]*dv12r; 
  outr[31] += -1.224744871391589*Ghat[15]*dv12r; 

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
  outl[20] += -1.224744871391589*Ghat[5]*dv12l; 
  outl[21] += -1.224744871391589*Ghat[6]*dv12l; 
  outl[22] += -1.224744871391589*Ghat[7]*dv12l; 
  outl[23] += -1.224744871391589*Ghat[8]*dv12l; 
  outl[24] += -1.224744871391589*Ghat[9]*dv12l; 
  outl[25] += -1.224744871391589*Ghat[10]*dv12l; 
  outl[26] += -0.7071067811865475*Ghat[15]*dv12l; 
  outl[27] += -1.224744871391589*Ghat[11]*dv12l; 
  outl[28] += -1.224744871391589*Ghat[12]*dv12l; 
  outl[29] += -1.224744871391589*Ghat[13]*dv12l; 
  outl[30] += -1.224744871391589*Ghat[14]*dv12l; 
  outl[31] += -1.224744871391589*Ghat[15]*dv12l; 

  return std::abs(amid); 
} 

__host__ __device__ double VlasovPhiBextSurf2x3vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double qDm, const double *phi, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qDm:       Species charge (q) divided by its mass (m).
  // phi:       electrostatic potential.
  // EM:        external EM field vectors.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double rdx2qDm = 2.*qDm/dxvl[0]; 
  const double rdy2qDm = 2.*qDm/dxvl[1]; 
  double dv10l = 2./dxvl[2]; 
  double dv10r = 2./dxvr[2]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *E0 = &EM[0]; 
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

  alpha[0] = (-2.0*B1[0]*wv3)+2.0*B2[0]*wv2-3.464101615137754*phi[1]*rdx2qDm+2.0*E0[0]; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = (-2.0*B1[2]*wv3)+2.0*B2[2]*wv2-3.464101615137754*phi[3]*rdx2qDm+2.0*E0[2]; 
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

__host__ __device__ double VlasovPhiBextSurf2x3vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double qDm, const double *phi, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qDm:       Species charge (q) divided by its mass (m).
  // phi:       electrostatic potential.
  // EM:        external EM field vectors.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double rdx2qDm = 2.*qDm/dxvl[0]; 
  const double rdy2qDm = 2.*qDm/dxvl[1]; 
  double dv11l = 2./dxvl[3]; 
  double dv11r = 2./dxvr[3]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *E1 = &EM[4]; 
  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double Ghat[16]; 
  double favg[16]; 
  double alpha[16]; 

  favg[0] = (-1.224744871391589*fr[4])+1.224744871391589*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[9])+1.224744871391589*fl[9]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[10])+1.224744871391589*fl[10]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[11])+1.224744871391589*fl[11]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[17])+1.224744871391589*fl[17]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[18])+1.224744871391589*fl[18]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = (-1.224744871391589*fr[19])+1.224744871391589*fl[19]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[8] = (-1.224744871391589*fr[23])+1.224744871391589*fl[23]+0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = (-1.224744871391589*fr[24])+1.224744871391589*fl[24]+0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 
  favg[10] = (-1.224744871391589*fr[25])+1.224744871391589*fl[25]+0.7071067811865475*fr[14]+0.7071067811865475*fl[14]; 
  favg[11] = (-1.224744871391589*fr[26])+1.224744871391589*fl[26]+0.7071067811865475*fr[16]+0.7071067811865475*fl[16]; 
  favg[12] = (-1.224744871391589*fr[28])+1.224744871391589*fl[28]+0.7071067811865475*fr[20]+0.7071067811865475*fl[20]; 
  favg[13] = (-1.224744871391589*fr[29])+1.224744871391589*fl[29]+0.7071067811865475*fr[21]+0.7071067811865475*fl[21]; 
  favg[14] = (-1.224744871391589*fr[30])+1.224744871391589*fl[30]+0.7071067811865475*fr[22]+0.7071067811865475*fl[22]; 
  favg[15] = (-1.224744871391589*fr[31])+1.224744871391589*fl[31]+0.7071067811865475*fr[27]+0.7071067811865475*fl[27]; 

  alpha[0] = 2.0*B0[0]*wv3-2.0*B2[0]*wv1-3.464101615137754*phi[2]*rdy2qDm+2.0*E1[0]; 
  alpha[1] = 2.0*B0[1]*wv3-2.0*B2[1]*wv1-3.464101615137754*phi[3]*rdy2qDm+2.0*E1[1]; 
  alpha[2] = 2.0*B0[2]*wv3-2.0*B2[2]*wv1+2.0*E1[2]; 
  alpha[3] = -0.5773502691896258*B2[0]*dv1; 
  alpha[4] = 0.5773502691896258*B0[0]*dv3; 
  alpha[5] = 2.0*B0[3]*wv3-2.0*B2[3]*wv1+2.0*E1[3]; 
  alpha[6] = -0.5773502691896258*B2[1]*dv1; 
  alpha[7] = -0.5773502691896258*B2[2]*dv1; 
  alpha[8] = 0.5773502691896258*B0[1]*dv3; 
  alpha[9] = 0.5773502691896258*B0[2]*dv3; 
  alpha[11] = -0.5773502691896258*B2[3]*dv1; 
  alpha[12] = 0.5773502691896258*B0[3]*dv3; 

  const double amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[4]+fl[4])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[9]+fl[9])-1.0*fr[1]+fl[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[10]+fl[10])-1.0*fr[2]+fl[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[11]+fl[11])-1.0*fr[3]+fl[3])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[5]+fl[5])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[17]+fl[17])-1.0*fr[6]+fl[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[18]+fl[18])-1.0*fr[7]+fl[7])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[19]+fl[19])-1.0*fr[8]+fl[8])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = 0.3535533905932737*(1.732050807568877*(fr[23]+fl[23])-1.0*fr[12]+fl[12])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = 0.3535533905932737*(1.732050807568877*(fr[24]+fl[24])-1.0*fr[13]+fl[13])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = 0.3535533905932737*(1.732050807568877*(fr[25]+fl[25])-1.0*fr[14]+fl[14])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = 0.3535533905932737*(1.732050807568877*(fr[26]+fl[26])-1.0*fr[16]+fl[16])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = 0.3535533905932737*(1.732050807568877*(fr[28]+fl[28])-1.0*fr[20]+fl[20])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = 0.3535533905932737*(1.732050807568877*(fr[29]+fl[29])-1.0*fr[21]+fl[21])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = 0.3535533905932737*(1.732050807568877*(fr[30]+fl[30])-1.0*fr[22]+fl[22])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = 0.3535533905932737*(1.732050807568877*(fr[31]+fl[31])-1.0*fr[27]+fl[27])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

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
  outr[17] += -1.224744871391589*Ghat[5]*dv11r; 
  outr[18] += -1.224744871391589*Ghat[6]*dv11r; 
  outr[19] += -1.224744871391589*Ghat[7]*dv11r; 
  outr[20] += 0.7071067811865475*Ghat[12]*dv11r; 
  outr[21] += 0.7071067811865475*Ghat[13]*dv11r; 
  outr[22] += 0.7071067811865475*Ghat[14]*dv11r; 
  outr[23] += -1.224744871391589*Ghat[8]*dv11r; 
  outr[24] += -1.224744871391589*Ghat[9]*dv11r; 
  outr[25] += -1.224744871391589*Ghat[10]*dv11r; 
  outr[26] += -1.224744871391589*Ghat[11]*dv11r; 
  outr[27] += 0.7071067811865475*Ghat[15]*dv11r; 
  outr[28] += -1.224744871391589*Ghat[12]*dv11r; 
  outr[29] += -1.224744871391589*Ghat[13]*dv11r; 
  outr[30] += -1.224744871391589*Ghat[14]*dv11r; 
  outr[31] += -1.224744871391589*Ghat[15]*dv11r; 

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
  outl[17] += -1.224744871391589*Ghat[5]*dv11l; 
  outl[18] += -1.224744871391589*Ghat[6]*dv11l; 
  outl[19] += -1.224744871391589*Ghat[7]*dv11l; 
  outl[20] += -0.7071067811865475*Ghat[12]*dv11l; 
  outl[21] += -0.7071067811865475*Ghat[13]*dv11l; 
  outl[22] += -0.7071067811865475*Ghat[14]*dv11l; 
  outl[23] += -1.224744871391589*Ghat[8]*dv11l; 
  outl[24] += -1.224744871391589*Ghat[9]*dv11l; 
  outl[25] += -1.224744871391589*Ghat[10]*dv11l; 
  outl[26] += -1.224744871391589*Ghat[11]*dv11l; 
  outl[27] += -0.7071067811865475*Ghat[15]*dv11l; 
  outl[28] += -1.224744871391589*Ghat[12]*dv11l; 
  outl[29] += -1.224744871391589*Ghat[13]*dv11l; 
  outl[30] += -1.224744871391589*Ghat[14]*dv11l; 
  outl[31] += -1.224744871391589*Ghat[15]*dv11l; 

  return std::abs(amid); 
} 

__host__ __device__ double VlasovPhiBextSurf2x3vSer_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double qDm, const double *phi, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qDm:       Species charge (q) divided by its mass (m).
  // phi:       electrostatic potential.
  // EM:        external EM field vectors.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double rdx2qDm = 2.*qDm/dxvl[0]; 
  const double rdy2qDm = 2.*qDm/dxvl[1]; 
  double dv12l = 2./dxvl[4]; 
  double dv12r = 2./dxvr[4]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double dv3 = dxvr[4], wv3 = wr[4]; 

  const double *E2 = &EM[8]; 
  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double Ghat[16]; 
  double favg[16]; 
  double alpha[16]; 

  favg[0] = (-1.224744871391589*fr[5])+1.224744871391589*fl[5]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = (-1.224744871391589*fr[12])+1.224744871391589*fl[12]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = (-1.224744871391589*fr[13])+1.224744871391589*fl[13]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = (-1.224744871391589*fr[14])+1.224744871391589*fl[14]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[5] = (-1.224744871391589*fr[20])+1.224744871391589*fl[20]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[21])+1.224744871391589*fl[21]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = (-1.224744871391589*fr[22])+1.224744871391589*fl[22]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[8] = (-1.224744871391589*fr[23])+1.224744871391589*fl[23]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[9] = (-1.224744871391589*fr[24])+1.224744871391589*fl[24]+0.7071067811865475*fr[10]+0.7071067811865475*fl[10]; 
  favg[10] = (-1.224744871391589*fr[25])+1.224744871391589*fl[25]+0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 
  favg[11] = (-1.224744871391589*fr[27])+1.224744871391589*fl[27]+0.7071067811865475*fr[16]+0.7071067811865475*fl[16]; 
  favg[12] = (-1.224744871391589*fr[28])+1.224744871391589*fl[28]+0.7071067811865475*fr[17]+0.7071067811865475*fl[17]; 
  favg[13] = (-1.224744871391589*fr[29])+1.224744871391589*fl[29]+0.7071067811865475*fr[18]+0.7071067811865475*fl[18]; 
  favg[14] = (-1.224744871391589*fr[30])+1.224744871391589*fl[30]+0.7071067811865475*fr[19]+0.7071067811865475*fl[19]; 
  favg[15] = (-1.224744871391589*fr[31])+1.224744871391589*fl[31]+0.7071067811865475*fr[26]+0.7071067811865475*fl[26]; 

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
  alpha[11] = 0.5773502691896258*B1[3]*dv1; 
  alpha[12] = -0.5773502691896258*B0[3]*dv2; 

  const double amid = 0.25*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(1.732050807568877*(fr[5]+fl[5])-1.0*fr[0]+fl[0])*amax+0.125*(alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[9]*favg[9]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.3535533905932737*(1.732050807568877*(fr[12]+fl[12])-1.0*fr[1]+fl[1])*amax+0.125*(alpha[9]*favg[12]+favg[9]*alpha[12]+alpha[7]*favg[11]+favg[7]*alpha[11]+alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[0]*favg[1]+favg[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(1.732050807568877*(fr[13]+fl[13])-1.0*fr[2]+fl[2])*amax+0.125*(alpha[8]*favg[12]+favg[8]*alpha[12]+alpha[6]*favg[11]+favg[6]*alpha[11]+alpha[4]*favg[9]+favg[4]*alpha[9]+alpha[3]*favg[7]+favg[3]*alpha[7]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[2]+favg[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(1.732050807568877*(fr[14]+fl[14])-1.0*fr[3]+fl[3])*amax+0.125*(alpha[12]*favg[15]+alpha[9]*favg[14]+alpha[8]*favg[13]+alpha[5]*favg[11]+favg[5]*alpha[11]+alpha[4]*favg[10]+alpha[2]*favg[7]+favg[2]*alpha[7]+alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[0]*favg[3]+favg[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[4]+fl[4])*amax+0.125*(alpha[11]*favg[15]+alpha[7]*favg[14]+alpha[6]*favg[13]+alpha[5]*favg[12]+favg[5]*alpha[12]+alpha[3]*favg[10]+alpha[2]*favg[9]+favg[2]*alpha[9]+alpha[1]*favg[8]+favg[1]*alpha[8]+alpha[0]*favg[4]+favg[0]*alpha[4]); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[20]+fl[20])-1.0*fr[6]+fl[6])*amax+0.125*(alpha[4]*favg[12]+favg[4]*alpha[12]+alpha[3]*favg[11]+favg[3]*alpha[11]+alpha[8]*favg[9]+favg[8]*alpha[9]+alpha[6]*favg[7]+favg[6]*alpha[7]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[2]+favg[1]*alpha[2]); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[21]+fl[21])-1.0*fr[7]+fl[7])*amax+0.125*(alpha[9]*favg[15]+alpha[12]*favg[14]+alpha[4]*favg[13]+alpha[2]*favg[11]+favg[2]*alpha[11]+alpha[8]*favg[10]+alpha[5]*favg[7]+favg[5]*alpha[7]+alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[1]*favg[3]+favg[1]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(1.732050807568877*(fr[22]+fl[22])-1.0*fr[8]+fl[8])*amax+0.125*(alpha[8]*favg[15]+alpha[4]*favg[14]+alpha[12]*favg[13]+alpha[1]*favg[11]+favg[1]*alpha[11]+alpha[9]*favg[10]+alpha[0]*favg[7]+favg[0]*alpha[7]+alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[2]*favg[3]+favg[2]*alpha[3]); 
  Ghat[8] = 0.3535533905932737*(1.732050807568877*(fr[23]+fl[23])-1.0*fr[9]+fl[9])*amax+0.125*(alpha[7]*favg[15]+alpha[11]*favg[14]+alpha[3]*favg[13]+alpha[2]*favg[12]+favg[2]*alpha[12]+alpha[6]*favg[10]+alpha[5]*favg[9]+favg[5]*alpha[9]+alpha[0]*favg[8]+favg[0]*alpha[8]+alpha[1]*favg[4]+favg[1]*alpha[4]); 
  Ghat[9] = 0.3535533905932737*(1.732050807568877*(fr[24]+fl[24])-1.0*fr[10]+fl[10])*amax+0.125*(alpha[6]*favg[15]+alpha[3]*favg[14]+alpha[11]*favg[13]+alpha[1]*favg[12]+favg[1]*alpha[12]+alpha[7]*favg[10]+alpha[0]*favg[9]+favg[0]*alpha[9]+alpha[5]*favg[8]+favg[5]*alpha[8]+alpha[2]*favg[4]+favg[2]*alpha[4]); 
  Ghat[10] = 0.3535533905932737*(1.732050807568877*(fr[25]+fl[25])-1.0*fr[11]+fl[11])*amax+0.125*(alpha[5]*favg[15]+alpha[2]*favg[14]+alpha[1]*favg[13]+alpha[11]*favg[12]+favg[11]*alpha[12]+alpha[0]*favg[10]+alpha[7]*favg[9]+favg[7]*alpha[9]+alpha[6]*favg[8]+favg[6]*alpha[8]+alpha[3]*favg[4]+favg[3]*alpha[4]); 
  Ghat[11] = 0.3535533905932737*(1.732050807568877*(fr[27]+fl[27])-1.0*fr[16]+fl[16])*amax+0.125*(alpha[4]*favg[15]+alpha[8]*favg[14]+alpha[9]*favg[13]+favg[10]*alpha[12]+alpha[0]*favg[11]+favg[0]*alpha[11]+alpha[1]*favg[7]+favg[1]*alpha[7]+alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[3]*favg[5]+favg[3]*alpha[5]); 
  Ghat[12] = 0.3535533905932737*(1.732050807568877*(fr[28]+fl[28])-1.0*fr[17]+fl[17])*amax+0.125*(alpha[3]*favg[15]+alpha[6]*favg[14]+alpha[7]*favg[13]+alpha[0]*favg[12]+favg[0]*alpha[12]+favg[10]*alpha[11]+alpha[1]*favg[9]+favg[1]*alpha[9]+alpha[2]*favg[8]+favg[2]*alpha[8]+alpha[4]*favg[5]+favg[4]*alpha[5]); 
  Ghat[13] = 0.3535533905932737*(1.732050807568877*(fr[29]+fl[29])-1.0*fr[18]+fl[18])*amax+0.125*(alpha[2]*favg[15]+alpha[5]*favg[14]+alpha[0]*favg[13]+alpha[7]*favg[12]+favg[7]*alpha[12]+alpha[9]*favg[11]+favg[9]*alpha[11]+alpha[1]*favg[10]+alpha[3]*favg[8]+favg[3]*alpha[8]+alpha[4]*favg[6]+favg[4]*alpha[6]); 
  Ghat[14] = 0.3535533905932737*(1.732050807568877*(fr[30]+fl[30])-1.0*fr[19]+fl[19])*amax+0.125*(alpha[1]*favg[15]+alpha[0]*favg[14]+alpha[5]*favg[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[8]*favg[11]+favg[8]*alpha[11]+alpha[2]*favg[10]+alpha[3]*favg[9]+favg[3]*alpha[9]+alpha[4]*favg[7]+favg[4]*alpha[7]); 
  Ghat[15] = 0.3535533905932737*(1.732050807568877*(fr[31]+fl[31])-1.0*fr[26]+fl[26])*amax+0.125*(alpha[0]*favg[15]+alpha[1]*favg[14]+alpha[2]*favg[13]+alpha[3]*favg[12]+favg[3]*alpha[12]+alpha[4]*favg[11]+favg[4]*alpha[11]+alpha[5]*favg[10]+alpha[6]*favg[9]+favg[6]*alpha[9]+alpha[7]*favg[8]+favg[7]*alpha[8]); 

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
  outr[20] += -1.224744871391589*Ghat[5]*dv12r; 
  outr[21] += -1.224744871391589*Ghat[6]*dv12r; 
  outr[22] += -1.224744871391589*Ghat[7]*dv12r; 
  outr[23] += -1.224744871391589*Ghat[8]*dv12r; 
  outr[24] += -1.224744871391589*Ghat[9]*dv12r; 
  outr[25] += -1.224744871391589*Ghat[10]*dv12r; 
  outr[26] += 0.7071067811865475*Ghat[15]*dv12r; 
  outr[27] += -1.224744871391589*Ghat[11]*dv12r; 
  outr[28] += -1.224744871391589*Ghat[12]*dv12r; 
  outr[29] += -1.224744871391589*Ghat[13]*dv12r; 
  outr[30] += -1.224744871391589*Ghat[14]*dv12r; 
  outr[31] += -1.224744871391589*Ghat[15]*dv12r; 

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
  outl[20] += -1.224744871391589*Ghat[5]*dv12l; 
  outl[21] += -1.224744871391589*Ghat[6]*dv12l; 
  outl[22] += -1.224744871391589*Ghat[7]*dv12l; 
  outl[23] += -1.224744871391589*Ghat[8]*dv12l; 
  outl[24] += -1.224744871391589*Ghat[9]*dv12l; 
  outl[25] += -1.224744871391589*Ghat[10]*dv12l; 
  outl[26] += -0.7071067811865475*Ghat[15]*dv12l; 
  outl[27] += -1.224744871391589*Ghat[11]*dv12l; 
  outl[28] += -1.224744871391589*Ghat[12]*dv12l; 
  outl[29] += -1.224744871391589*Ghat[13]*dv12l; 
  outl[30] += -1.224744871391589*Ghat[14]*dv12l; 
  outl[31] += -1.224744871391589*Ghat[15]*dv12l; 

  return std::abs(amid); 
} 

