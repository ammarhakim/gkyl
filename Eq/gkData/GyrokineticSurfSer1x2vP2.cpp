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

  double alpha[8]; 
  alpha[0] = 2.0*wv; 
  alpha[1] = 1.154700538379252/dfac_v; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  double favg[8]; 
  favg[0] = 0.7071067811865475*(2.23606797749979*(fr[7]+fl[7])+1.732050807568877*(fl[1]-1.0*fr[1])+fr[0]+fl[0]); 
  favg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[11]+fl[11])+3.0*(fl[4]-1.0*fr[4]))+3.0*(fr[2]+fl[2])); 
  favg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[13]+fl[13])+3.0*(fl[5]-1.0*fr[5]))+3.0*(fr[3]+fl[3])); 
  favg[3] = 0.7071067811865475*(2.23606797749979*(fr[17]+fl[17])+1.732050807568877*(fl[10]-1.0*fr[10])+fr[6]+fl[6]); 
  favg[4] = -0.1414213562373095*(8.660254037844387*fr[12]-1.0*(8.660254037844387*fl[12]+5.0*(fr[8]+fl[8]))); 
  favg[5] = -0.1414213562373095*(8.660254037844387*fr[15]-1.0*(8.660254037844387*fl[15]+5.0*(fr[9]+fl[9]))); 
  favg[6] = -0.1414213562373095*(8.660254037844387*fr[18]-1.0*(8.660254037844387*fl[18]+5.0*(fr[14]+fl[14]))); 
  favg[7] = -0.1414213562373095*(8.660254037844387*fr[19]-1.0*(8.660254037844387*fl[19]+5.0*(fr[16]+fl[16]))); 
  Ghat[0] = 0.25*(alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[7]-1.0*(2.23606797749979*fl[7]+1.732050807568877*(fr[1]+fl[1]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.05*(4.47213595499958*alpha[1]*favg[4]+5.0*(alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[11]-1.0*(3.872983346207417*fl[11]+3.0*(fr[4]+fl[4])))+3.0*(fr[2]-1.0*fl[2]))*amax; 
  Ghat[2] = 0.25*(alpha[1]*favg[3]+alpha[0]*favg[2])-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[13]-1.0*(3.872983346207417*fl[13]+3.0*(fr[5]+fl[5])))+3.0*(fr[3]-1.0*fl[3]))*amax; 
  Ghat[3] = 0.01666666666666667*(13.41640786499874*alpha[1]*favg[6]+15.0*(alpha[0]*favg[3]+alpha[1]*favg[2]))-0.3535533905932737*(2.23606797749979*fr[17]-1.0*(2.23606797749979*fl[17]+1.732050807568877*(fr[10]+fl[10]))+fr[6]-1.0*fl[6])*amax; 
  Ghat[4] = 0.07071067811865474*(8.660254037844387*(fr[12]+fl[12])+5.0*(fl[8]-1.0*fr[8]))*amax+0.05*(5.0*alpha[0]*favg[4]+4.47213595499958*alpha[1]*favg[1]); 
  Ghat[5] = 0.07071067811865474*(8.660254037844387*(fr[15]+fl[15])+5.0*(fl[9]-1.0*fr[9]))*amax+0.01666666666666667*(15.0*alpha[1]*favg[7]+15.0*alpha[0]*favg[5]); 
  Ghat[6] = 0.07071067811865474*(8.660254037844387*(fr[18]+fl[18])+5.0*(fl[14]-1.0*fr[14]))*amax+0.01666666666666667*(15.0*alpha[0]*favg[6]+13.41640786499874*alpha[1]*favg[3]); 
  Ghat[7] = 0.07071067811865474*(8.660254037844387*(fr[19]+fl[19])+5.0*(fl[16]-1.0*fr[16]))*amax+0.01666666666666667*(15.0*alpha[0]*favg[7]+15.0*alpha[1]*favg[5]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_x; 
  incr[1] = -1.224744871391589*Ghat[0]*dfac_x; 
  incr[2] = 0.7071067811865475*Ghat[1]*dfac_x; 
  incr[3] = 0.7071067811865475*Ghat[2]*dfac_x; 
  incr[4] = -1.224744871391589*Ghat[1]*dfac_x; 
  incr[5] = -1.224744871391589*Ghat[2]*dfac_x; 
  incr[6] = 0.7071067811865475*Ghat[3]*dfac_x; 
  incr[7] = 1.58113883008419*Ghat[0]*dfac_x; 
  incr[8] = 0.7071067811865475*Ghat[4]*dfac_x; 
  incr[9] = 0.7071067811865475*Ghat[5]*dfac_x; 
  incr[10] = -1.224744871391589*Ghat[3]*dfac_x; 
  incr[11] = 1.58113883008419*Ghat[1]*dfac_x; 
  incr[12] = -1.224744871391589*Ghat[4]*dfac_x; 
  incr[13] = 1.58113883008419*Ghat[2]*dfac_x; 
  incr[14] = 0.7071067811865475*Ghat[6]*dfac_x; 
  incr[15] = -1.224744871391589*Ghat[5]*dfac_x; 
  incr[16] = 0.7071067811865475*Ghat[7]*dfac_x; 
  incr[17] = 1.58113883008419*Ghat[3]*dfac_x; 
  incr[18] = -1.224744871391589*Ghat[6]*dfac_x; 
  incr[19] = -1.224744871391589*Ghat[7]*dfac_x; 

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

  double alpha[8]; 
  alpha[0] = -(2.449489742783178*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(5.477225575051662*Phi[2]*dfac_x*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  double favg[8]; 
  favg[0] = 0.7071067811865475*(2.23606797749979*(fr[8]+fl[8])+1.732050807568877*(fl[2]-1.0*fr[2])+fr[0]+fl[0]); 
  favg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[12]+fl[12])+3.0*(fl[4]-1.0*fr[4]))+3.0*(fr[1]+fl[1])); 
  favg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[14]+fl[14])+3.0*(fl[6]-1.0*fr[6]))+3.0*(fr[3]+fl[3])); 
  favg[3] = 0.7071067811865475*(2.23606797749979*(fr[18]+fl[18])+1.732050807568877*(fl[10]-1.0*fr[10])+fr[5]+fl[5]); 
  favg[4] = -0.1414213562373095*(8.660254037844387*fr[11]-1.0*(8.660254037844387*fl[11]+5.0*(fr[7]+fl[7]))); 
  favg[5] = -0.1414213562373095*(8.660254037844387*fr[16]-1.0*(8.660254037844387*fl[16]+5.0*(fr[9]+fl[9]))); 
  favg[6] = -0.1414213562373095*(8.660254037844387*fr[17]-1.0*(8.660254037844387*fl[17]+5.0*(fr[13]+fl[13]))); 
  favg[7] = -0.1414213562373095*(8.660254037844387*fr[19]-1.0*(8.660254037844387*fl[19]+5.0*(fr[15]+fl[15]))); 
  Ghat[0] = 0.25*(alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[8]-1.0*(2.23606797749979*fl[8]+1.732050807568877*(fr[2]+fl[2]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.05*(4.47213595499958*alpha[1]*favg[4]+5.0*(alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[12]-1.0*(3.872983346207417*fl[12]+3.0*(fr[4]+fl[4])))+3.0*(fr[1]-1.0*fl[1]))*amax; 
  Ghat[2] = 0.25*(alpha[1]*favg[3]+alpha[0]*favg[2])-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[14]-1.0*(3.872983346207417*fl[14]+3.0*(fr[6]+fl[6])))+3.0*(fr[3]-1.0*fl[3]))*amax; 
  Ghat[3] = 0.01666666666666667*(13.41640786499874*alpha[1]*favg[6]+15.0*(alpha[0]*favg[3]+alpha[1]*favg[2]))-0.3535533905932737*(2.23606797749979*fr[18]-1.0*(2.23606797749979*fl[18]+1.732050807568877*(fr[10]+fl[10]))+fr[5]-1.0*fl[5])*amax; 
  Ghat[4] = 0.07071067811865474*(8.660254037844387*(fr[11]+fl[11])+5.0*(fl[7]-1.0*fr[7]))*amax+0.05*(5.0*alpha[0]*favg[4]+4.47213595499958*alpha[1]*favg[1]); 
  Ghat[5] = 0.07071067811865474*(8.660254037844387*(fr[16]+fl[16])+5.0*(fl[9]-1.0*fr[9]))*amax+0.01666666666666667*(15.0*alpha[1]*favg[7]+15.0*alpha[0]*favg[5]); 
  Ghat[6] = 0.07071067811865474*(8.660254037844387*(fr[17]+fl[17])+5.0*(fl[13]-1.0*fr[13]))*amax+0.01666666666666667*(15.0*alpha[0]*favg[6]+13.41640786499874*alpha[1]*favg[3]); 
  Ghat[7] = 0.07071067811865474*(8.660254037844387*(fr[19]+fl[19])+5.0*(fl[15]-1.0*fr[15]))*amax+0.01666666666666667*(15.0*alpha[0]*favg[7]+15.0*alpha[1]*favg[5]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_v; 
  incr[1] = 0.7071067811865475*Ghat[1]*dfac_v; 
  incr[2] = -1.224744871391589*Ghat[0]*dfac_v; 
  incr[3] = 0.7071067811865475*Ghat[2]*dfac_v; 
  incr[4] = -1.224744871391589*Ghat[1]*dfac_v; 
  incr[5] = 0.7071067811865475*Ghat[3]*dfac_v; 
  incr[6] = -1.224744871391589*Ghat[2]*dfac_v; 
  incr[7] = 0.7071067811865475*Ghat[4]*dfac_v; 
  incr[8] = 1.58113883008419*Ghat[0]*dfac_v; 
  incr[9] = 0.7071067811865475*Ghat[5]*dfac_v; 
  incr[10] = -1.224744871391589*Ghat[3]*dfac_v; 
  incr[11] = -1.224744871391589*Ghat[4]*dfac_v; 
  incr[12] = 1.58113883008419*Ghat[1]*dfac_v; 
  incr[13] = 0.7071067811865475*Ghat[6]*dfac_v; 
  incr[14] = 1.58113883008419*Ghat[2]*dfac_v; 
  incr[15] = 0.7071067811865475*Ghat[7]*dfac_v; 
  incr[16] = -1.224744871391589*Ghat[5]*dfac_v; 
  incr[17] = -1.224744871391589*Ghat[6]*dfac_v; 
  incr[18] = 1.58113883008419*Ghat[3]*dfac_v; 
  incr[19] = -1.224744871391589*Ghat[7]*dfac_v; 

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

  double alpha[8]; 
  alpha[0] = 2.0*wv; 
  alpha[1] = 1.154700538379252/dfac_v; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  double favg[8]; 
  favg[0] = 0.7071067811865475*(2.23606797749979*(fr[7]+fl[7])+1.732050807568877*(fl[1]-1.0*fr[1])+fr[0]+fl[0]); 
  favg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[11]+fl[11])+3.0*(fl[4]-1.0*fr[4]))+3.0*(fr[2]+fl[2])); 
  favg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[13]+fl[13])+3.0*(fl[5]-1.0*fr[5]))+3.0*(fr[3]+fl[3])); 
  favg[3] = 0.7071067811865475*(2.23606797749979*(fr[17]+fl[17])+1.732050807568877*(fl[10]-1.0*fr[10])+fr[6]+fl[6]); 
  favg[4] = -0.1414213562373095*(8.660254037844387*fr[12]-1.0*(8.660254037844387*fl[12]+5.0*(fr[8]+fl[8]))); 
  favg[5] = -0.1414213562373095*(8.660254037844387*fr[15]-1.0*(8.660254037844387*fl[15]+5.0*(fr[9]+fl[9]))); 
  favg[6] = -0.1414213562373095*(8.660254037844387*fr[18]-1.0*(8.660254037844387*fl[18]+5.0*(fr[14]+fl[14]))); 
  favg[7] = -0.1414213562373095*(8.660254037844387*fr[19]-1.0*(8.660254037844387*fl[19]+5.0*(fr[16]+fl[16]))); 
  Ghat[0] = 0.25*(alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[7]-1.0*(2.23606797749979*fl[7]+1.732050807568877*(fr[1]+fl[1]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.05*(4.47213595499958*alpha[1]*favg[4]+5.0*(alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[11]-1.0*(3.872983346207417*fl[11]+3.0*(fr[4]+fl[4])))+3.0*(fr[2]-1.0*fl[2]))*amax; 
  Ghat[2] = 0.25*(alpha[1]*favg[3]+alpha[0]*favg[2])-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[13]-1.0*(3.872983346207417*fl[13]+3.0*(fr[5]+fl[5])))+3.0*(fr[3]-1.0*fl[3]))*amax; 
  Ghat[3] = 0.01666666666666667*(13.41640786499874*alpha[1]*favg[6]+15.0*(alpha[0]*favg[3]+alpha[1]*favg[2]))-0.3535533905932737*(2.23606797749979*fr[17]-1.0*(2.23606797749979*fl[17]+1.732050807568877*(fr[10]+fl[10]))+fr[6]-1.0*fl[6])*amax; 
  Ghat[4] = 0.07071067811865474*(8.660254037844387*(fr[12]+fl[12])+5.0*(fl[8]-1.0*fr[8]))*amax+0.05*(5.0*alpha[0]*favg[4]+4.47213595499958*alpha[1]*favg[1]); 
  Ghat[5] = 0.07071067811865474*(8.660254037844387*(fr[15]+fl[15])+5.0*(fl[9]-1.0*fr[9]))*amax+0.01666666666666667*(15.0*alpha[1]*favg[7]+15.0*alpha[0]*favg[5]); 
  Ghat[6] = 0.07071067811865474*(8.660254037844387*(fr[18]+fl[18])+5.0*(fl[14]-1.0*fr[14]))*amax+0.01666666666666667*(15.0*alpha[0]*favg[6]+13.41640786499874*alpha[1]*favg[3]); 
  Ghat[7] = 0.07071067811865474*(8.660254037844387*(fr[19]+fl[19])+5.0*(fl[16]-1.0*fr[16]))*amax+0.01666666666666667*(15.0*alpha[0]*favg[7]+15.0*alpha[1]*favg[5]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_x; 
  incr[1] = -1.224744871391589*Ghat[0]*dfac_x; 
  incr[2] = 0.7071067811865475*Ghat[1]*dfac_x; 
  incr[3] = 0.7071067811865475*Ghat[2]*dfac_x; 
  incr[4] = -1.224744871391589*Ghat[1]*dfac_x; 
  incr[5] = -1.224744871391589*Ghat[2]*dfac_x; 
  incr[6] = 0.7071067811865475*Ghat[3]*dfac_x; 
  incr[7] = 1.58113883008419*Ghat[0]*dfac_x; 
  incr[8] = 0.7071067811865475*Ghat[4]*dfac_x; 
  incr[9] = 0.7071067811865475*Ghat[5]*dfac_x; 
  incr[10] = -1.224744871391589*Ghat[3]*dfac_x; 
  incr[11] = 1.58113883008419*Ghat[1]*dfac_x; 
  incr[12] = -1.224744871391589*Ghat[4]*dfac_x; 
  incr[13] = 1.58113883008419*Ghat[2]*dfac_x; 
  incr[14] = 0.7071067811865475*Ghat[6]*dfac_x; 
  incr[15] = -1.224744871391589*Ghat[5]*dfac_x; 
  incr[16] = 0.7071067811865475*Ghat[7]*dfac_x; 
  incr[17] = 1.58113883008419*Ghat[3]*dfac_x; 
  incr[18] = -1.224744871391589*Ghat[6]*dfac_x; 
  incr[19] = -1.224744871391589*Ghat[7]*dfac_x; 

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

  double alpha[8]; 
  alpha[0] = (-(2.449489742783178*Bmag[1]*dfac_x*wm)/m_)-(2.449489742783178*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = (-(5.477225575051662*Bmag[2]*dfac_x*wm)/m_)-(5.477225575051662*Phi[2]*dfac_x*q_)/m_; 
  alpha[2] = -(1.414213562373095*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[3] = -(3.16227766016838*Bmag[2]*dfac_x)/(dfac_m*m_); 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  double favg[8]; 
  favg[0] = 0.7071067811865475*(2.23606797749979*(fr[8]+fl[8])+1.732050807568877*(fl[2]-1.0*fr[2])+fr[0]+fl[0]); 
  favg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[12]+fl[12])+3.0*(fl[4]-1.0*fr[4]))+3.0*(fr[1]+fl[1])); 
  favg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[14]+fl[14])+3.0*(fl[6]-1.0*fr[6]))+3.0*(fr[3]+fl[3])); 
  favg[3] = 0.7071067811865475*(2.23606797749979*(fr[18]+fl[18])+1.732050807568877*(fl[10]-1.0*fr[10])+fr[5]+fl[5]); 
  favg[4] = -0.1414213562373095*(8.660254037844387*fr[11]-1.0*(8.660254037844387*fl[11]+5.0*(fr[7]+fl[7]))); 
  favg[5] = -0.1414213562373095*(8.660254037844387*fr[16]-1.0*(8.660254037844387*fl[16]+5.0*(fr[9]+fl[9]))); 
  favg[6] = -0.1414213562373095*(8.660254037844387*fr[17]-1.0*(8.660254037844387*fl[17]+5.0*(fr[13]+fl[13]))); 
  favg[7] = -0.1414213562373095*(8.660254037844387*fr[19]-1.0*(8.660254037844387*fl[19]+5.0*(fr[15]+fl[15]))); 
  Ghat[0] = 0.25*(alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])-0.3535533905932737*(2.23606797749979*fr[8]-1.0*(2.23606797749979*fl[8]+1.732050807568877*(fr[2]+fl[2]))+fr[0]-1.0*fl[0])*amax; 
  Ghat[1] = 0.01666666666666667*(13.41640786499874*alpha[3]*favg[6]+3.0*(4.47213595499958*alpha[1]*favg[4]+5.0*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1])))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[12]-1.0*(3.872983346207417*fl[12]+3.0*(fr[4]+fl[4])))+3.0*(fr[1]-1.0*fl[1]))*amax; 
  Ghat[2] = 0.01666666666666667*(13.41640786499874*alpha[3]*favg[7]+3.0*(4.47213595499958*alpha[2]*favg[5]+5.0*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2])))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[14]-1.0*(3.872983346207417*fl[14]+3.0*(fr[6]+fl[6])))+3.0*(fr[3]-1.0*fl[3]))*amax; 
  Ghat[3] = 0.01666666666666667*(13.41640786499874*(alpha[2]*favg[7]+alpha[1]*favg[6])+3.0*(4.47213595499958*alpha[3]*(favg[5]+favg[4])+5.0*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2])))-0.3535533905932737*(2.23606797749979*fr[18]-1.0*(2.23606797749979*fl[18]+1.732050807568877*(fr[10]+fl[10]))+fr[5]-1.0*fl[5])*amax; 
  Ghat[4] = 0.07071067811865474*(8.660254037844387*(fr[11]+fl[11])+5.0*(fl[7]-1.0*fr[7]))*amax+0.01666666666666667*(15.0*alpha[2]*favg[6]+3.0*(5.0*alpha[0]*favg[4]+4.47213595499958*(alpha[3]*favg[3]+alpha[1]*favg[1]))); 
  Ghat[5] = 0.07071067811865474*(8.660254037844387*(fr[16]+fl[16])+5.0*(fl[9]-1.0*fr[9]))*amax+0.01666666666666667*(15.0*alpha[1]*favg[7]+3.0*(5.0*alpha[0]*favg[5]+4.47213595499958*(alpha[3]*favg[3]+alpha[2]*favg[2]))); 
  Ghat[6] = 0.07071067811865474*(8.660254037844387*(fr[17]+fl[17])+5.0*(fl[13]-1.0*fr[13]))*amax+0.01666666666666667*(3.0*(4.0*alpha[3]*favg[7]+5.0*alpha[0]*favg[6])+6.708203932499369*(2.23606797749979*alpha[2]*favg[4]+2.0*(alpha[1]*favg[3]+favg[1]*alpha[3]))); 
  Ghat[7] = 0.07071067811865474*(8.660254037844387*(fr[19]+fl[19])+5.0*(fl[15]-1.0*fr[15]))*amax+0.01666666666666667*(3.0*(5.0*alpha[0]*favg[7]+4.0*alpha[3]*favg[6])+6.708203932499369*(2.23606797749979*alpha[1]*favg[5]+2.0*(alpha[2]*favg[3]+favg[2]*alpha[3]))); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_v; 
  incr[1] = 0.7071067811865475*Ghat[1]*dfac_v; 
  incr[2] = -1.224744871391589*Ghat[0]*dfac_v; 
  incr[3] = 0.7071067811865475*Ghat[2]*dfac_v; 
  incr[4] = -1.224744871391589*Ghat[1]*dfac_v; 
  incr[5] = 0.7071067811865475*Ghat[3]*dfac_v; 
  incr[6] = -1.224744871391589*Ghat[2]*dfac_v; 
  incr[7] = 0.7071067811865475*Ghat[4]*dfac_v; 
  incr[8] = 1.58113883008419*Ghat[0]*dfac_v; 
  incr[9] = 0.7071067811865475*Ghat[5]*dfac_v; 
  incr[10] = -1.224744871391589*Ghat[3]*dfac_v; 
  incr[11] = -1.224744871391589*Ghat[4]*dfac_v; 
  incr[12] = 1.58113883008419*Ghat[1]*dfac_v; 
  incr[13] = 0.7071067811865475*Ghat[6]*dfac_v; 
  incr[14] = 1.58113883008419*Ghat[2]*dfac_v; 
  incr[15] = 0.7071067811865475*Ghat[7]*dfac_v; 
  incr[16] = -1.224744871391589*Ghat[5]*dfac_v; 
  incr[17] = -1.224744871391589*Ghat[6]*dfac_v; 
  incr[18] = 1.58113883008419*Ghat[3]*dfac_v; 
  incr[19] = -1.224744871391589*Ghat[7]*dfac_v; 

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
