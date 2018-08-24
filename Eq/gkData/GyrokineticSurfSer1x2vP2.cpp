#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x2vSer_X_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.5*wv; 

  double alpha[20]; 
  alpha[0] = 2.828427124746191*wv; 
  alpha[2] = 1.632993161855453/dfac_v; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  double favg[20]; 

  favg[0] = 1.0*fr[0]+fl[0]; 
  favg[1] = -1.0*fr[1]+fl[1]; 
  favg[2] = 1.0*fr[2]+fl[2]; 
  favg[3] = 1.0*fr[3]+fl[3]; 
  favg[4] = -1.0*fr[4]+fl[4]; 
  favg[5] = -1.0*fr[5]+fl[5]; 
  favg[6] = 1.0*fr[6]+fl[6]; 
  favg[7] = 1.0*fr[7]+fl[7]; 
  favg[8] = 1.0*fr[8]+fl[8]; 
  favg[9] = 1.0*fr[9]+fl[9]; 
  favg[10] = -1.0*fr[10]+fl[10]; 
  favg[11] = 1.0*fr[11]+fl[11]; 
  favg[12] = -1.0*fr[12]+fl[12]; 
  favg[13] = 1.0*fr[13]+fl[13]; 
  favg[14] = 1.0*fr[14]+fl[14]; 
  favg[15] = -1.0*fr[15]+fl[15]; 
  favg[16] = 1.0*fr[16]+fl[16]; 
  favg[17] = 1.0*fr[17]+fl[17]; 
  favg[18] = -1.0*fr[18]+fl[18]; 
  favg[19] = -1.0*fr[19]+fl[19]; 
  double fjump[20]; 

  fjump[0] = amax*(1.0*fr[0]-fl[0]); 
  fjump[1] = amax*(-1.0*fr[1]-fl[1]); 
  fjump[2] = amax*(1.0*fr[2]-fl[2]); 
  fjump[3] = amax*(1.0*fr[3]-fl[3]); 
  fjump[4] = amax*(-1.0*fr[4]-fl[4]); 
  fjump[5] = amax*(-1.0*fr[5]-fl[5]); 
  fjump[6] = amax*(1.0*fr[6]-fl[6]); 
  fjump[7] = amax*(1.0*fr[7]-fl[7]); 
  fjump[8] = amax*(1.0*fr[8]-fl[8]); 
  fjump[9] = amax*(1.0*fr[9]-fl[9]); 
  fjump[10] = amax*(-1.0*fr[10]-fl[10]); 
  fjump[11] = amax*(1.0*fr[11]-fl[11]); 
  fjump[12] = amax*(-1.0*fr[12]-fl[12]); 
  fjump[13] = amax*(1.0*fr[13]-fl[13]); 
  fjump[14] = amax*(1.0*fr[14]-fl[14]); 
  fjump[15] = amax*(-1.0*fr[15]-fl[15]); 
  fjump[16] = amax*(1.0*fr[16]-fl[16]); 
  fjump[17] = amax*(1.0*fr[17]-fl[17]); 
  fjump[18] = amax*(-1.0*fr[18]-fl[18]); 
  fjump[19] = amax*(-1.0*fr[19]-fl[19]); 
  Ghat[0] = alpha[2]*(0.3952847075210473*favg[11]+0.3061862178478971*favg[4]+0.1767766952966368*favg[2])-1.118033988749895*fjump[7]+alpha[0]*(0.3952847075210473*favg[7]+0.3061862178478971*favg[1]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[1]-0.5*fjump[0]; 
  Ghat[2] = alpha[2]*(0.273861278752583*favg[12]+0.1581138830084189*favg[8]+0.3952847075210473*favg[7]+0.3061862178478971*favg[1]+0.1767766952966368*favg[0])-1.118033988749895*fjump[11]+alpha[0]*(0.3952847075210473*favg[11]+0.3061862178478971*favg[4]+0.1767766952966368*favg[2])-0.8660254037844386*fjump[4]-0.5*fjump[2]; 
  Ghat[3] = alpha[2]*(0.3952847075210473*favg[17]+0.3061862178478971*favg[10]+0.1767766952966368*favg[6])-1.118033988749895*fjump[13]+alpha[0]*(0.3952847075210473*favg[13]+0.3061862178478971*favg[5]+0.1767766952966368*favg[3])-0.8660254037844386*fjump[5]-0.5*fjump[3]; 
  Ghat[6] = alpha[2]*(0.273861278752583*favg[18]+0.1581138830084189*favg[14]+0.3952847075210473*favg[13]+0.3061862178478971*favg[5]+0.1767766952966368*favg[3])-1.118033988749895*fjump[17]+alpha[0]*(0.3952847075210473*favg[17]+0.3061862178478971*favg[10]+0.1767766952966368*favg[6])-0.8660254037844386*fjump[10]-0.5*fjump[6]; 
  Ghat[8] = (-0.8660254037844387*fjump[12])+alpha[0]*(0.3061862178478971*favg[12]+0.1767766952966368*favg[8])+alpha[2]*(0.3535533905932737*favg[11]+0.273861278752583*favg[4]+0.1581138830084189*favg[2])-0.5*fjump[8]; 
  Ghat[9] = alpha[2]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[16])-0.8660254037844387*fjump[15]+alpha[0]*(0.3061862178478971*favg[15]+0.1767766952966368*favg[9])-0.5*fjump[9]; 
  Ghat[14] = (-0.8660254037844387*fjump[18])+alpha[0]*(0.3061862178478971*favg[18]+0.1767766952966368*favg[14])+alpha[2]*(0.3535533905932737*favg[17]+0.273861278752583*favg[10]+0.1581138830084189*favg[6])-0.5*fjump[14]; 
  Ghat[16] = (-0.8660254037844387*fjump[19])+alpha[0]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[16])-0.5*fjump[16]+alpha[2]*(0.3061862178478971*favg[15]+0.1767766952966368*favg[9]); 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = 0.5*Ghat[3]*dfac_x; 
  incr[4] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[5] = -0.8660254037844386*Ghat[3]*dfac_x; 
  incr[6] = 0.5*Ghat[6]*dfac_x; 
  incr[7] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[8] = 0.5*Ghat[8]*dfac_x; 
  incr[9] = 0.5*Ghat[9]*dfac_x; 
  incr[10] = -0.8660254037844386*Ghat[6]*dfac_x; 
  incr[11] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[12] = -0.8660254037844387*Ghat[8]*dfac_x; 
  incr[13] = 1.118033988749895*Ghat[3]*dfac_x; 
  incr[14] = 0.5*Ghat[14]*dfac_x; 
  incr[15] = -0.8660254037844387*Ghat[9]*dfac_x; 
  incr[16] = 0.5*Ghat[16]*dfac_x; 
  incr[17] = 1.118033988749895*Ghat[6]*dfac_x; 
  incr[18] = -0.8660254037844387*Ghat[14]*dfac_x; 
  incr[19] = -0.8660254037844387*Ghat[16]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.6123724356957944*Phi[1]*dfac_x*q_)/m_; 

  double alpha[20]; 
  alpha[0] = -(3.464101615137754*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(7.745966692414834*Phi[2]*dfac_x*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  double favg[20]; 

  favg[0] = 1.0*fr[0]+fl[0]; 
  favg[1] = 1.0*fr[1]+fl[1]; 
  favg[2] = -1.0*fr[2]+fl[2]; 
  favg[3] = 1.0*fr[3]+fl[3]; 
  favg[4] = -1.0*fr[4]+fl[4]; 
  favg[5] = 1.0*fr[5]+fl[5]; 
  favg[6] = -1.0*fr[6]+fl[6]; 
  favg[7] = 1.0*fr[7]+fl[7]; 
  favg[8] = 1.0*fr[8]+fl[8]; 
  favg[9] = 1.0*fr[9]+fl[9]; 
  favg[10] = -1.0*fr[10]+fl[10]; 
  favg[11] = -1.0*fr[11]+fl[11]; 
  favg[12] = 1.0*fr[12]+fl[12]; 
  favg[13] = 1.0*fr[13]+fl[13]; 
  favg[14] = 1.0*fr[14]+fl[14]; 
  favg[15] = 1.0*fr[15]+fl[15]; 
  favg[16] = -1.0*fr[16]+fl[16]; 
  favg[17] = -1.0*fr[17]+fl[17]; 
  favg[18] = 1.0*fr[18]+fl[18]; 
  favg[19] = -1.0*fr[19]+fl[19]; 
  double fjump[20]; 

  fjump[0] = amax*(1.0*fr[0]-fl[0]); 
  fjump[1] = amax*(1.0*fr[1]-fl[1]); 
  fjump[2] = amax*(-1.0*fr[2]-fl[2]); 
  fjump[3] = amax*(1.0*fr[3]-fl[3]); 
  fjump[4] = amax*(-1.0*fr[4]-fl[4]); 
  fjump[5] = amax*(1.0*fr[5]-fl[5]); 
  fjump[6] = amax*(-1.0*fr[6]-fl[6]); 
  fjump[7] = amax*(1.0*fr[7]-fl[7]); 
  fjump[8] = amax*(1.0*fr[8]-fl[8]); 
  fjump[9] = amax*(1.0*fr[9]-fl[9]); 
  fjump[10] = amax*(-1.0*fr[10]-fl[10]); 
  fjump[11] = amax*(-1.0*fr[11]-fl[11]); 
  fjump[12] = amax*(1.0*fr[12]-fl[12]); 
  fjump[13] = amax*(1.0*fr[13]-fl[13]); 
  fjump[14] = amax*(1.0*fr[14]-fl[14]); 
  fjump[15] = amax*(1.0*fr[15]-fl[15]); 
  fjump[16] = amax*(-1.0*fr[16]-fl[16]); 
  fjump[17] = amax*(-1.0*fr[17]-fl[17]); 
  fjump[18] = amax*(1.0*fr[18]-fl[18]); 
  fjump[19] = amax*(-1.0*fr[19]-fl[19]); 
  Ghat[0] = alpha[1]*(0.3952847075210473*favg[12]+0.3061862178478971*favg[4]+0.1767766952966368*favg[1])-1.118033988749895*fjump[8]+alpha[0]*(0.3952847075210473*favg[8]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (-1.118033988749895*fjump[12])+alpha[0]*(0.3952847075210473*favg[12]+0.3061862178478971*favg[4]+0.1767766952966368*favg[1])+alpha[1]*(0.273861278752583*favg[11]+0.3952847075210473*favg[8]+0.1581138830084189*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[4]-0.5*fjump[1]; 
  Ghat[3] = alpha[1]*(0.3952847075210473*favg[18]+0.3061862178478971*favg[10]+0.1767766952966368*favg[5])-1.118033988749895*fjump[14]+alpha[0]*(0.3952847075210473*favg[14]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])-0.8660254037844386*fjump[6]-0.5*fjump[3]; 
  Ghat[5] = (-1.118033988749895*fjump[18])+alpha[0]*(0.3952847075210473*favg[18]+0.3061862178478971*favg[10]+0.1767766952966368*favg[5])+alpha[1]*(0.273861278752583*favg[17]+0.3952847075210473*favg[14]+0.1581138830084189*favg[13]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])-0.8660254037844386*fjump[10]-0.5*fjump[5]; 
  Ghat[7] = alpha[1]*(0.3535533905932737*favg[12]+0.273861278752583*favg[4]+0.1581138830084189*favg[1])-0.8660254037844387*fjump[11]+alpha[0]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[7])-0.5*fjump[7]; 
  Ghat[9] = alpha[1]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[15])-0.8660254037844387*fjump[16]+alpha[0]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[9])-0.5*fjump[9]; 
  Ghat[13] = alpha[1]*(0.3535533905932737*favg[18]+0.273861278752583*favg[10]+0.1581138830084189*favg[5])-0.8660254037844387*fjump[17]+alpha[0]*(0.3061862178478971*favg[17]+0.1767766952966368*favg[13])-0.5*fjump[13]; 
  Ghat[15] = (-0.8660254037844387*fjump[19])+alpha[0]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[15])+alpha[1]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[9])-0.5*fjump[15]; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = 0.5*Ghat[3]*dfac_v; 
  incr[4] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[5] = 0.5*Ghat[5]*dfac_v; 
  incr[6] = -0.8660254037844386*Ghat[3]*dfac_v; 
  incr[7] = 0.5*Ghat[7]*dfac_v; 
  incr[8] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[9] = 0.5*Ghat[9]*dfac_v; 
  incr[10] = -0.8660254037844386*Ghat[5]*dfac_v; 
  incr[11] = -0.8660254037844387*Ghat[7]*dfac_v; 
  incr[12] = 1.118033988749895*Ghat[1]*dfac_v; 
  incr[13] = 0.5*Ghat[13]*dfac_v; 
  incr[14] = 1.118033988749895*Ghat[3]*dfac_v; 
  incr[15] = 0.5*Ghat[15]*dfac_v; 
  incr[16] = -0.8660254037844387*Ghat[9]*dfac_v; 
  incr[17] = -0.8660254037844387*Ghat[13]*dfac_v; 
  incr[18] = 1.118033988749895*Ghat[5]*dfac_v; 
  incr[19] = -0.8660254037844387*Ghat[15]*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_X_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.5*wv; 

  double alpha[20]; 
  alpha[0] = 2.828427124746191*wv; 
  alpha[2] = 1.632993161855453/dfac_v; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  double favg[20]; 

  favg[0] = 1.0*fr[0]+fl[0]; 
  favg[1] = -1.0*fr[1]+fl[1]; 
  favg[2] = 1.0*fr[2]+fl[2]; 
  favg[3] = 1.0*fr[3]+fl[3]; 
  favg[4] = -1.0*fr[4]+fl[4]; 
  favg[5] = -1.0*fr[5]+fl[5]; 
  favg[6] = 1.0*fr[6]+fl[6]; 
  favg[7] = 1.0*fr[7]+fl[7]; 
  favg[8] = 1.0*fr[8]+fl[8]; 
  favg[9] = 1.0*fr[9]+fl[9]; 
  favg[10] = -1.0*fr[10]+fl[10]; 
  favg[11] = 1.0*fr[11]+fl[11]; 
  favg[12] = -1.0*fr[12]+fl[12]; 
  favg[13] = 1.0*fr[13]+fl[13]; 
  favg[14] = 1.0*fr[14]+fl[14]; 
  favg[15] = -1.0*fr[15]+fl[15]; 
  favg[16] = 1.0*fr[16]+fl[16]; 
  favg[17] = 1.0*fr[17]+fl[17]; 
  favg[18] = -1.0*fr[18]+fl[18]; 
  favg[19] = -1.0*fr[19]+fl[19]; 
  double fjump[20]; 

  fjump[0] = amax*(1.0*fr[0]-fl[0]); 
  fjump[1] = amax*(-1.0*fr[1]-fl[1]); 
  fjump[2] = amax*(1.0*fr[2]-fl[2]); 
  fjump[3] = amax*(1.0*fr[3]-fl[3]); 
  fjump[4] = amax*(-1.0*fr[4]-fl[4]); 
  fjump[5] = amax*(-1.0*fr[5]-fl[5]); 
  fjump[6] = amax*(1.0*fr[6]-fl[6]); 
  fjump[7] = amax*(1.0*fr[7]-fl[7]); 
  fjump[8] = amax*(1.0*fr[8]-fl[8]); 
  fjump[9] = amax*(1.0*fr[9]-fl[9]); 
  fjump[10] = amax*(-1.0*fr[10]-fl[10]); 
  fjump[11] = amax*(1.0*fr[11]-fl[11]); 
  fjump[12] = amax*(-1.0*fr[12]-fl[12]); 
  fjump[13] = amax*(1.0*fr[13]-fl[13]); 
  fjump[14] = amax*(1.0*fr[14]-fl[14]); 
  fjump[15] = amax*(-1.0*fr[15]-fl[15]); 
  fjump[16] = amax*(1.0*fr[16]-fl[16]); 
  fjump[17] = amax*(1.0*fr[17]-fl[17]); 
  fjump[18] = amax*(-1.0*fr[18]-fl[18]); 
  fjump[19] = amax*(-1.0*fr[19]-fl[19]); 
  Ghat[0] = alpha[2]*(0.3952847075210473*favg[11]+0.3061862178478971*favg[4]+0.1767766952966368*favg[2])-1.118033988749895*fjump[7]+alpha[0]*(0.3952847075210473*favg[7]+0.3061862178478971*favg[1]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[1]-0.5*fjump[0]; 
  Ghat[2] = alpha[2]*(0.273861278752583*favg[12]+0.1581138830084189*favg[8]+0.3952847075210473*favg[7]+0.3061862178478971*favg[1]+0.1767766952966368*favg[0])-1.118033988749895*fjump[11]+alpha[0]*(0.3952847075210473*favg[11]+0.3061862178478971*favg[4]+0.1767766952966368*favg[2])-0.8660254037844386*fjump[4]-0.5*fjump[2]; 
  Ghat[3] = alpha[2]*(0.3952847075210473*favg[17]+0.3061862178478971*favg[10]+0.1767766952966368*favg[6])-1.118033988749895*fjump[13]+alpha[0]*(0.3952847075210473*favg[13]+0.3061862178478971*favg[5]+0.1767766952966368*favg[3])-0.8660254037844386*fjump[5]-0.5*fjump[3]; 
  Ghat[6] = alpha[2]*(0.273861278752583*favg[18]+0.1581138830084189*favg[14]+0.3952847075210473*favg[13]+0.3061862178478971*favg[5]+0.1767766952966368*favg[3])-1.118033988749895*fjump[17]+alpha[0]*(0.3952847075210473*favg[17]+0.3061862178478971*favg[10]+0.1767766952966368*favg[6])-0.8660254037844386*fjump[10]-0.5*fjump[6]; 
  Ghat[8] = (-0.8660254037844387*fjump[12])+alpha[0]*(0.3061862178478971*favg[12]+0.1767766952966368*favg[8])+alpha[2]*(0.3535533905932737*favg[11]+0.273861278752583*favg[4]+0.1581138830084189*favg[2])-0.5*fjump[8]; 
  Ghat[9] = alpha[2]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[16])-0.8660254037844387*fjump[15]+alpha[0]*(0.3061862178478971*favg[15]+0.1767766952966368*favg[9])-0.5*fjump[9]; 
  Ghat[14] = (-0.8660254037844387*fjump[18])+alpha[0]*(0.3061862178478971*favg[18]+0.1767766952966368*favg[14])+alpha[2]*(0.3535533905932737*favg[17]+0.273861278752583*favg[10]+0.1581138830084189*favg[6])-0.5*fjump[14]; 
  Ghat[16] = (-0.8660254037844387*fjump[19])+alpha[0]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[16])-0.5*fjump[16]+alpha[2]*(0.3061862178478971*favg[15]+0.1767766952966368*favg[9]); 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = 0.5*Ghat[3]*dfac_x; 
  incr[4] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[5] = -0.8660254037844386*Ghat[3]*dfac_x; 
  incr[6] = 0.5*Ghat[6]*dfac_x; 
  incr[7] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[8] = 0.5*Ghat[8]*dfac_x; 
  incr[9] = 0.5*Ghat[9]*dfac_x; 
  incr[10] = -0.8660254037844386*Ghat[6]*dfac_x; 
  incr[11] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[12] = -0.8660254037844387*Ghat[8]*dfac_x; 
  incr[13] = 1.118033988749895*Ghat[3]*dfac_x; 
  incr[14] = 0.5*Ghat[14]*dfac_x; 
  incr[15] = -0.8660254037844387*Ghat[9]*dfac_x; 
  incr[16] = 0.5*Ghat[16]*dfac_x; 
  incr[17] = 1.118033988749895*Ghat[6]*dfac_x; 
  incr[18] = -0.8660254037844387*Ghat[14]*dfac_x; 
  incr[19] = -0.8660254037844387*Ghat[16]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.6123724356957944*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 

  double alpha[20]; 
  alpha[0] = (-(3.464101615137754*Bmag[1]*dfac_x*wm)/m_)-(3.464101615137754*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = (-(7.745966692414834*Bmag[2]*dfac_x*wm)/m_)-(7.745966692414834*Phi[2]*dfac_x*q_)/m_; 
  alpha[3] = -(2.0*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[5] = -(4.47213595499958*Bmag[2]*dfac_x)/(dfac_m*m_); 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[20]; 
  double favg[20]; 

  favg[0] = 1.0*fr[0]+fl[0]; 
  favg[1] = 1.0*fr[1]+fl[1]; 
  favg[2] = -1.0*fr[2]+fl[2]; 
  favg[3] = 1.0*fr[3]+fl[3]; 
  favg[4] = -1.0*fr[4]+fl[4]; 
  favg[5] = 1.0*fr[5]+fl[5]; 
  favg[6] = -1.0*fr[6]+fl[6]; 
  favg[7] = 1.0*fr[7]+fl[7]; 
  favg[8] = 1.0*fr[8]+fl[8]; 
  favg[9] = 1.0*fr[9]+fl[9]; 
  favg[10] = -1.0*fr[10]+fl[10]; 
  favg[11] = -1.0*fr[11]+fl[11]; 
  favg[12] = 1.0*fr[12]+fl[12]; 
  favg[13] = 1.0*fr[13]+fl[13]; 
  favg[14] = 1.0*fr[14]+fl[14]; 
  favg[15] = 1.0*fr[15]+fl[15]; 
  favg[16] = -1.0*fr[16]+fl[16]; 
  favg[17] = -1.0*fr[17]+fl[17]; 
  favg[18] = 1.0*fr[18]+fl[18]; 
  favg[19] = -1.0*fr[19]+fl[19]; 
  double fjump[20]; 

  fjump[0] = amax*(1.0*fr[0]-fl[0]); 
  fjump[1] = amax*(1.0*fr[1]-fl[1]); 
  fjump[2] = amax*(-1.0*fr[2]-fl[2]); 
  fjump[3] = amax*(1.0*fr[3]-fl[3]); 
  fjump[4] = amax*(-1.0*fr[4]-fl[4]); 
  fjump[5] = amax*(1.0*fr[5]-fl[5]); 
  fjump[6] = amax*(-1.0*fr[6]-fl[6]); 
  fjump[7] = amax*(1.0*fr[7]-fl[7]); 
  fjump[8] = amax*(1.0*fr[8]-fl[8]); 
  fjump[9] = amax*(1.0*fr[9]-fl[9]); 
  fjump[10] = amax*(-1.0*fr[10]-fl[10]); 
  fjump[11] = amax*(-1.0*fr[11]-fl[11]); 
  fjump[12] = amax*(1.0*fr[12]-fl[12]); 
  fjump[13] = amax*(1.0*fr[13]-fl[13]); 
  fjump[14] = amax*(1.0*fr[14]-fl[14]); 
  fjump[15] = amax*(1.0*fr[15]-fl[15]); 
  fjump[16] = amax*(-1.0*fr[16]-fl[16]); 
  fjump[17] = amax*(-1.0*fr[17]-fl[17]); 
  fjump[18] = amax*(1.0*fr[18]-fl[18]); 
  fjump[19] = amax*(-1.0*fr[19]-fl[19]); 
  Ghat[0] = alpha[5]*(0.3952847075210473*favg[18]+0.3061862178478971*favg[10]+0.1767766952966368*favg[5])+alpha[3]*(0.3952847075210473*favg[14]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+alpha[1]*(0.3952847075210473*favg[12]+0.3061862178478971*favg[4]+0.1767766952966368*favg[1])-1.118033988749895*fjump[8]+alpha[0]*(0.3952847075210473*favg[8]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = alpha[3]*(0.3952847075210473*favg[18]+0.3061862178478971*favg[10]+0.1767766952966368*favg[5])+alpha[5]*(0.273861278752583*favg[17]+0.3952847075210473*favg[14]+0.1581138830084189*favg[13]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])-1.118033988749895*fjump[12]+alpha[0]*(0.3952847075210473*favg[12]+0.3061862178478971*favg[4]+0.1767766952966368*favg[1])+alpha[1]*(0.273861278752583*favg[11]+0.3952847075210473*favg[8]+0.1581138830084189*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[4]-0.5*fjump[1]; 
  Ghat[3] = alpha[5]*(0.273861278752583*favg[19]+0.1581138830084189*favg[15]+0.3952847075210473*favg[12]+0.3061862178478971*favg[4]+0.1767766952966368*favg[1])+alpha[1]*(0.3952847075210473*favg[18]+0.3061862178478971*favg[10]+0.1767766952966368*favg[5])+alpha[3]*(0.273861278752583*favg[16]+0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-1.118033988749895*fjump[14]+alpha[0]*(0.3952847075210473*favg[14]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])-0.8660254037844386*fjump[6]-0.5*fjump[3]; 
  Ghat[5] = alpha[3]*(0.273861278752583*favg[19]+0.1581138830084189*favg[15]+0.3952847075210473*favg[12]+0.3061862178478971*favg[4]+0.1767766952966368*favg[1])-1.118033988749895*fjump[18]+alpha[0]*(0.3952847075210473*favg[18]+0.3061862178478971*favg[10]+0.1767766952966368*favg[5])+alpha[1]*(0.273861278752583*favg[17]+0.3952847075210473*favg[14]+0.1581138830084189*favg[13]+0.3061862178478971*favg[6]+0.1767766952966368*favg[3])+alpha[5]*(0.273861278752583*favg[16]+0.273861278752583*favg[11]+0.1581138830084189*favg[9]+0.3952847075210473*favg[8]+0.1581138830084189*favg[7]+0.3061862178478971*favg[2]+0.1767766952966368*favg[0])-0.8660254037844386*fjump[10]-0.5*fjump[5]; 
  Ghat[7] = alpha[5]*(0.3535533905932737*favg[18]+0.273861278752583*favg[10]+0.1581138830084189*favg[5])+alpha[3]*(0.3061862178478971*favg[17]+0.1767766952966368*favg[13])+alpha[1]*(0.3535533905932737*favg[12]+0.273861278752583*favg[4]+0.1581138830084189*favg[1])-0.8660254037844387*fjump[11]+alpha[0]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[7])-0.5*fjump[7]; 
  Ghat[9] = alpha[1]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[15])+alpha[5]*(0.3535533905932737*favg[18]+0.273861278752583*favg[10]+0.1581138830084189*favg[5])-0.8660254037844387*fjump[16]+alpha[0]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[9])+alpha[3]*(0.3535533905932737*favg[14]+0.273861278752583*favg[6]+0.1581138830084189*favg[3])-0.5*fjump[9]; 
  Ghat[13] = alpha[5]*(0.2449489742783177*favg[19]+0.1414213562373095*favg[15]+0.3535533905932737*favg[12]+0.273861278752583*favg[4]+0.1581138830084189*favg[1])+alpha[1]*(0.3535533905932737*favg[18]+0.273861278752583*favg[10]+0.1581138830084189*favg[5])-0.8660254037844387*fjump[17]+alpha[0]*(0.3061862178478971*favg[17]+0.1767766952966368*favg[13])-0.5*fjump[13]+alpha[3]*(0.3061862178478971*favg[11]+0.1767766952966368*favg[7]); 
  Ghat[15] = (-0.8660254037844387*fjump[19])+alpha[0]*(0.3061862178478971*favg[19]+0.1767766952966368*favg[15])+alpha[3]*(0.3535533905932737*favg[18]+0.273861278752583*favg[10]+0.1581138830084189*favg[5])+alpha[5]*(0.2449489742783177*favg[17]+0.3535533905932737*favg[14]+0.1414213562373095*favg[13]+0.273861278752583*favg[6]+0.1581138830084189*favg[3])+alpha[1]*(0.3061862178478971*favg[16]+0.1767766952966368*favg[9])-0.5*fjump[15]; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = 0.5*Ghat[3]*dfac_v; 
  incr[4] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[5] = 0.5*Ghat[5]*dfac_v; 
  incr[6] = -0.8660254037844386*Ghat[3]*dfac_v; 
  incr[7] = 0.5*Ghat[7]*dfac_v; 
  incr[8] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[9] = 0.5*Ghat[9]*dfac_v; 
  incr[10] = -0.8660254037844386*Ghat[5]*dfac_v; 
  incr[11] = -0.8660254037844387*Ghat[7]*dfac_v; 
  incr[12] = 1.118033988749895*Ghat[1]*dfac_v; 
  incr[13] = 0.5*Ghat[13]*dfac_v; 
  incr[14] = 1.118033988749895*Ghat[3]*dfac_v; 
  incr[15] = 0.5*Ghat[15]*dfac_v; 
  incr[16] = -0.8660254037844387*Ghat[9]*dfac_v; 
  incr[17] = -0.8660254037844387*Ghat[13]*dfac_v; 
  incr[18] = 1.118033988749895*Ghat[5]*dfac_v; 
  incr[19] = -0.8660254037844387*Ghat[15]*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += incr[19]; 
return std::abs(alpha0); 
} 
