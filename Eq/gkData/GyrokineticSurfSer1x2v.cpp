#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double alpha[8]; 
  alpha[0] = 2.828427124746191*wv; 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_X_P2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[20]; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double alpha[20]; 
  alpha[0] = 2.828427124746191*wv; 
  alpha[2] = 1.632993161855453/dfac_v; 
  if (alpha0>0) { 
  incr[0] = 0.05892556509887893*(6.708203932499369*alpha[2]*fl[11]+6.708203932499369*alpha[0]*fl[7]+5.196152422706631*alpha[2]*fl[4]+3.0*alpha[2]*fl[2]+5.196152422706631*alpha[0]*fl[1]+3.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*(3.872983346207417*alpha[2]*fl[11]+3.872983346207417*alpha[0]*fl[7]+3.0*alpha[2]*fl[4]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[12]+33.54101966249684*alpha[0]*fl[11]+13.41640786499874*alpha[2]*fl[8]+33.54101966249685*alpha[2]*fl[7]+25.98076211353316*alpha[0]*fl[4]+15.0*alpha[0]*fl[2]+25.98076211353316*fl[1]*alpha[2]+15.0*fl[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05892556509887893*(6.708203932499369*alpha[2]*fl[17]+6.708203932499369*alpha[0]*fl[13]+5.196152422706631*alpha[2]*fl[10]+3.0*alpha[2]*fl[6]+5.196152422706631*alpha[0]*fl[5]+3.0*alpha[0]*fl[3])*dfac_x; 
  incr[4] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[12]+19.36491673103708*alpha[0]*fl[11]+7.745966692414834*alpha[2]*fl[8]+19.36491673103709*alpha[2]*fl[7]+15.0*alpha[0]*fl[4]+8.660254037844386*alpha[0]*fl[2]+15.0*fl[1]*alpha[2]+8.660254037844386*fl[0]*alpha[2])*dfac_x; 
  incr[5] = -0.1767766952966368*(3.872983346207417*alpha[2]*fl[17]+3.872983346207417*alpha[0]*fl[13]+3.0*alpha[2]*fl[10]+1.732050807568877*alpha[2]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[3])*dfac_x; 
  incr[6] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[18]+33.54101966249685*alpha[0]*fl[17]+13.41640786499874*alpha[2]*fl[14]+33.54101966249684*alpha[2]*fl[13]+25.98076211353316*alpha[0]*fl[10]+15.0*alpha[0]*fl[6]+25.98076211353316*alpha[2]*fl[5]+15.0*alpha[2]*fl[3])*dfac_x; 
  incr[7] = 0.05892556509887893*(15.0*alpha[2]*fl[11]+15.0*alpha[0]*fl[7]+11.61895003862225*alpha[2]*fl[4]+6.708203932499369*alpha[2]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[8] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[12]+30.0*alpha[2]*fl[11]+15.0*alpha[0]*fl[8]+23.2379000772445*alpha[2]*fl[4]+13.41640786499874*alpha[2]*fl[2])*dfac_x; 
  incr[9] = 0.01178511301977579*(25.98076211353316*alpha[2]*fl[19]+15.0*alpha[2]*fl[16]+25.98076211353316*alpha[0]*fl[15]+15.0*alpha[0]*fl[9])*dfac_x; 
  incr[10] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[18]+19.36491673103709*alpha[0]*fl[17]+7.745966692414834*alpha[2]*fl[14]+19.36491673103708*alpha[2]*fl[13]+15.0*alpha[0]*fl[10]+8.660254037844386*alpha[0]*fl[6]+15.0*alpha[2]*fl[5]+8.660254037844386*alpha[2]*fl[3])*dfac_x; 
  incr[11] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[12]+75.0*alpha[0]*fl[11]+30.0*alpha[2]*fl[8]+75.00000000000001*alpha[2]*fl[7]+58.09475019311126*alpha[0]*fl[4]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*fl[1]*alpha[2]+33.54101966249684*fl[0]*alpha[2])*dfac_x; 
  incr[12] = -0.03535533905932736*(15.0*alpha[0]*fl[12]+17.32050807568877*alpha[2]*fl[11]+8.660254037844387*alpha[0]*fl[8]+13.41640786499874*alpha[2]*fl[4]+7.745966692414834*alpha[2]*fl[2])*dfac_x; 
  incr[13] = 0.05892556509887893*(15.0*alpha[2]*fl[17]+15.0*alpha[0]*fl[13]+11.61895003862225*alpha[2]*fl[10]+6.708203932499369*alpha[2]*fl[6]+11.61895003862225*alpha[0]*fl[5]+6.708203932499369*alpha[0]*fl[3])*dfac_x; 
  incr[14] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[18]+30.0*alpha[2]*fl[17]+15.0*alpha[0]*fl[14]+23.2379000772445*alpha[2]*fl[10]+13.41640786499874*alpha[2]*fl[6])*dfac_x; 
  incr[15] = -0.03535533905932736*(15.0*alpha[2]*fl[19]+8.660254037844386*alpha[2]*fl[16]+15.0*alpha[0]*fl[15]+8.660254037844387*alpha[0]*fl[9])*dfac_x; 
  incr[16] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[19]+15.0*alpha[0]*fl[16]+25.98076211353316*alpha[2]*fl[15]+15.0*alpha[2]*fl[9])*dfac_x; 
  incr[17] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[18]+75.0*alpha[0]*fl[17]+30.0*alpha[2]*fl[14]+75.00000000000001*alpha[2]*fl[13]+58.09475019311126*alpha[0]*fl[10]+33.54101966249685*alpha[0]*fl[6]+58.09475019311126*alpha[2]*fl[5]+33.54101966249685*alpha[2]*fl[3])*dfac_x; 
  incr[18] = -0.03535533905932736*(15.0*alpha[0]*fl[18]+17.32050807568877*alpha[2]*fl[17]+8.660254037844387*alpha[0]*fl[14]+13.41640786499874*alpha[2]*fl[10]+7.745966692414834*alpha[2]*fl[6])*dfac_x; 
  incr[19] = -0.03535533905932736*(15.0*alpha[0]*fl[19]+8.660254037844387*alpha[0]*fl[16]+15.0*alpha[2]*fl[15]+8.660254037844386*alpha[2]*fl[9])*dfac_x; 

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
  } else { 
  incr[0] = 0.05892556509887893*(6.708203932499369*alpha[2]*fr[11]+6.708203932499369*alpha[0]*fr[7]-5.196152422706631*alpha[2]*fr[4]+3.0*alpha[2]*fr[2]-5.196152422706631*alpha[0]*fr[1]+3.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = -0.1767766952966368*(3.872983346207417*alpha[2]*fr[11]+3.872983346207417*alpha[0]*fr[7]-3.0*alpha[2]*fr[4]+1.732050807568877*alpha[2]*fr[2]-3.0*alpha[0]*fr[1]+1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[12]-33.54101966249684*alpha[0]*fr[11]-13.41640786499874*alpha[2]*fr[8]-33.54101966249685*alpha[2]*fr[7]+25.98076211353316*alpha[0]*fr[4]-15.0*alpha[0]*fr[2]+25.98076211353316*fr[1]*alpha[2]-15.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05892556509887893*(6.708203932499369*alpha[2]*fr[17]+6.708203932499369*alpha[0]*fr[13]-5.196152422706631*alpha[2]*fr[10]+3.0*alpha[2]*fr[6]-5.196152422706631*alpha[0]*fr[5]+3.0*alpha[0]*fr[3])*dfac_x; 
  incr[4] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[12]-19.36491673103708*alpha[0]*fr[11]-7.745966692414834*alpha[2]*fr[8]-19.36491673103709*alpha[2]*fr[7]+15.0*alpha[0]*fr[4]-8.660254037844386*alpha[0]*fr[2]+15.0*fr[1]*alpha[2]-8.660254037844386*fr[0]*alpha[2])*dfac_x; 
  incr[5] = -0.1767766952966368*(3.872983346207417*alpha[2]*fr[17]+3.872983346207417*alpha[0]*fr[13]-3.0*alpha[2]*fr[10]+1.732050807568877*alpha[2]*fr[6]-3.0*alpha[0]*fr[5]+1.732050807568877*alpha[0]*fr[3])*dfac_x; 
  incr[6] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[18]-33.54101966249685*alpha[0]*fr[17]-13.41640786499874*alpha[2]*fr[14]-33.54101966249684*alpha[2]*fr[13]+25.98076211353316*alpha[0]*fr[10]-15.0*alpha[0]*fr[6]+25.98076211353316*alpha[2]*fr[5]-15.0*alpha[2]*fr[3])*dfac_x; 
  incr[7] = 0.05892556509887893*(15.0*alpha[2]*fr[11]+15.0*alpha[0]*fr[7]-11.61895003862225*alpha[2]*fr[4]+6.708203932499369*alpha[2]*fr[2]-11.61895003862225*alpha[0]*fr[1]+6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[8] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[12]-30.0*alpha[2]*fr[11]-15.0*alpha[0]*fr[8]+23.2379000772445*alpha[2]*fr[4]-13.41640786499874*alpha[2]*fr[2])*dfac_x; 
  incr[9] = -0.01178511301977579*(25.98076211353316*alpha[2]*fr[19]-15.0*alpha[2]*fr[16]+25.98076211353316*alpha[0]*fr[15]-15.0*alpha[0]*fr[9])*dfac_x; 
  incr[10] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[18]-19.36491673103709*alpha[0]*fr[17]-7.745966692414834*alpha[2]*fr[14]-19.36491673103708*alpha[2]*fr[13]+15.0*alpha[0]*fr[10]-8.660254037844386*alpha[0]*fr[6]+15.0*alpha[2]*fr[5]-8.660254037844386*alpha[2]*fr[3])*dfac_x; 
  incr[11] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[12]-75.0*alpha[0]*fr[11]-30.0*alpha[2]*fr[8]-75.00000000000001*alpha[2]*fr[7]+58.09475019311126*alpha[0]*fr[4]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*fr[1]*alpha[2]-33.54101966249684*fr[0]*alpha[2])*dfac_x; 
  incr[12] = 0.03535533905932736*(15.0*alpha[0]*fr[12]-17.32050807568877*alpha[2]*fr[11]-8.660254037844387*alpha[0]*fr[8]+13.41640786499874*alpha[2]*fr[4]-7.745966692414834*alpha[2]*fr[2])*dfac_x; 
  incr[13] = 0.05892556509887893*(15.0*alpha[2]*fr[17]+15.0*alpha[0]*fr[13]-11.61895003862225*alpha[2]*fr[10]+6.708203932499369*alpha[2]*fr[6]-11.61895003862225*alpha[0]*fr[5]+6.708203932499369*alpha[0]*fr[3])*dfac_x; 
  incr[14] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[18]-30.0*alpha[2]*fr[17]-15.0*alpha[0]*fr[14]+23.2379000772445*alpha[2]*fr[10]-13.41640786499874*alpha[2]*fr[6])*dfac_x; 
  incr[15] = 0.03535533905932736*(15.0*alpha[2]*fr[19]-8.660254037844386*alpha[2]*fr[16]+15.0*alpha[0]*fr[15]-8.660254037844387*alpha[0]*fr[9])*dfac_x; 
  incr[16] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[19]-15.0*alpha[0]*fr[16]+25.98076211353316*alpha[2]*fr[15]-15.0*alpha[2]*fr[9])*dfac_x; 
  incr[17] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[18]-75.0*alpha[0]*fr[17]-30.0*alpha[2]*fr[14]-75.00000000000001*alpha[2]*fr[13]+58.09475019311126*alpha[0]*fr[10]-33.54101966249685*alpha[0]*fr[6]+58.09475019311126*alpha[2]*fr[5]-33.54101966249685*alpha[2]*fr[3])*dfac_x; 
  incr[18] = 0.03535533905932736*(15.0*alpha[0]*fr[18]-17.32050807568877*alpha[2]*fr[17]-8.660254037844387*alpha[0]*fr[14]+13.41640786499874*alpha[2]*fr[10]-7.745966692414834*alpha[2]*fr[6])*dfac_x; 
  incr[19] = 0.03535533905932736*(15.0*alpha[0]*fr[19]-8.660254037844387*alpha[0]*fr[16]+15.0*alpha[2]*fr[15]-8.660254037844386*alpha[2]*fr[9])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.732050807568877*Phi[1]*dfac_x*q_)/m_; 

  double alpha[8]; 
  alpha[0] = -(3.464101615137754*Phi[1]*dfac_x*q_)/m_; 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[4]+fl[1])*dfac_v; 
  incr[2] = -0.1767766952966368*alpha[0]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[6]+fl[3])*dfac_v; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[1])*dfac_v; 
  incr[5] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[7]+fl[5])*dfac_v; 
  incr[6] = -0.1767766952966368*alpha[0]*(3.0*fl[6]+1.732050807568877*fl[3])*dfac_v; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[5])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  incr[1] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[1])*dfac_v; 
  incr[2] = 0.1767766952966368*alpha[0]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[6]-1.0*fr[3])*dfac_v; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[1])*dfac_v; 
  incr[5] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[5])*dfac_v; 
  incr[6] = 0.1767766952966368*alpha[0]*(3.0*fr[6]-1.732050807568877*fr[3])*dfac_v; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[5])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[20]; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.732050807568877*Phi[1]*dfac_x*q_)/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-(3.061862178478972*Phi[2]*fr[12]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[8]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[8]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-1.118033988749895*fr[8]*amax+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Phi[1]*fr[12]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[12]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[11]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[11]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[8]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[8]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[7]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[4]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-1.118033988749895*fr[12]*amax+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[3] = (-(3.061862178478972*Phi[2]*fr[18]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[14]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[14]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[10]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[3]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[3]*dfac_x*q_)/m_-1.118033988749895*fr[14]*amax+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[5] = (-(1.369306393762915*Phi[1]*fr[18]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[18]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[17]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[17]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[14]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[14]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[13]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[13]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[10]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[10]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[3]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[3]*dfac_x*q_)/m_-1.118033988749895*fr[18]*amax+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 
  Ghat[7] = (-(2.738612787525831*Phi[2]*fr[12]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[12]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[11]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[11]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[7]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[4]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_+0.8660254037844387*fr[11]*amax+0.8660254037844387*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax; 
  Ghat[9] = (2.371708245126284*Phi[2]*fr[19]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[19]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[16]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[16]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[9]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[9]*dfac_x*q_)/m_+0.8660254037844387*fr[16]*amax+0.8660254037844387*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[13] = (-(2.738612787525831*Phi[2]*fr[18]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[18]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[17]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[17]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[13]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[13]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[10]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[5]*dfac_x*q_)/m_+0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax; 
  Ghat[15] = (1.060660171779821*Phi[1]*fr[19]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[19]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[16]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[16]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[9]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[9]*dfac_x*q_)/m_+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax; 

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
double EmGyrokineticSurf1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double alpha[8]; 
  alpha[0] = 2.828427124746191*wv; 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_X_P2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[20]; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double alpha[20]; 
  alpha[0] = 2.828427124746191*wv; 
  alpha[2] = 1.632993161855453/dfac_v; 
  if (alpha0>0) { 
  incr[0] = 0.05892556509887893*(6.708203932499369*alpha[2]*fl[11]+6.708203932499369*alpha[0]*fl[7]+5.196152422706631*alpha[2]*fl[4]+3.0*alpha[2]*fl[2]+5.196152422706631*alpha[0]*fl[1]+3.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*(3.872983346207417*alpha[2]*fl[11]+3.872983346207417*alpha[0]*fl[7]+3.0*alpha[2]*fl[4]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[12]+33.54101966249684*alpha[0]*fl[11]+13.41640786499874*alpha[2]*fl[8]+33.54101966249685*alpha[2]*fl[7]+25.98076211353316*alpha[0]*fl[4]+15.0*alpha[0]*fl[2]+25.98076211353316*fl[1]*alpha[2]+15.0*fl[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05892556509887893*(6.708203932499369*alpha[2]*fl[17]+6.708203932499369*alpha[0]*fl[13]+5.196152422706631*alpha[2]*fl[10]+3.0*alpha[2]*fl[6]+5.196152422706631*alpha[0]*fl[5]+3.0*alpha[0]*fl[3])*dfac_x; 
  incr[4] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[12]+19.36491673103708*alpha[0]*fl[11]+7.745966692414834*alpha[2]*fl[8]+19.36491673103709*alpha[2]*fl[7]+15.0*alpha[0]*fl[4]+8.660254037844386*alpha[0]*fl[2]+15.0*fl[1]*alpha[2]+8.660254037844386*fl[0]*alpha[2])*dfac_x; 
  incr[5] = -0.1767766952966368*(3.872983346207417*alpha[2]*fl[17]+3.872983346207417*alpha[0]*fl[13]+3.0*alpha[2]*fl[10]+1.732050807568877*alpha[2]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[3])*dfac_x; 
  incr[6] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[18]+33.54101966249685*alpha[0]*fl[17]+13.41640786499874*alpha[2]*fl[14]+33.54101966249684*alpha[2]*fl[13]+25.98076211353316*alpha[0]*fl[10]+15.0*alpha[0]*fl[6]+25.98076211353316*alpha[2]*fl[5]+15.0*alpha[2]*fl[3])*dfac_x; 
  incr[7] = 0.05892556509887893*(15.0*alpha[2]*fl[11]+15.0*alpha[0]*fl[7]+11.61895003862225*alpha[2]*fl[4]+6.708203932499369*alpha[2]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[8] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[12]+30.0*alpha[2]*fl[11]+15.0*alpha[0]*fl[8]+23.2379000772445*alpha[2]*fl[4]+13.41640786499874*alpha[2]*fl[2])*dfac_x; 
  incr[9] = 0.01178511301977579*(25.98076211353316*alpha[2]*fl[19]+15.0*alpha[2]*fl[16]+25.98076211353316*alpha[0]*fl[15]+15.0*alpha[0]*fl[9])*dfac_x; 
  incr[10] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[18]+19.36491673103709*alpha[0]*fl[17]+7.745966692414834*alpha[2]*fl[14]+19.36491673103708*alpha[2]*fl[13]+15.0*alpha[0]*fl[10]+8.660254037844386*alpha[0]*fl[6]+15.0*alpha[2]*fl[5]+8.660254037844386*alpha[2]*fl[3])*dfac_x; 
  incr[11] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[12]+75.0*alpha[0]*fl[11]+30.0*alpha[2]*fl[8]+75.00000000000001*alpha[2]*fl[7]+58.09475019311126*alpha[0]*fl[4]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*fl[1]*alpha[2]+33.54101966249684*fl[0]*alpha[2])*dfac_x; 
  incr[12] = -0.03535533905932736*(15.0*alpha[0]*fl[12]+17.32050807568877*alpha[2]*fl[11]+8.660254037844387*alpha[0]*fl[8]+13.41640786499874*alpha[2]*fl[4]+7.745966692414834*alpha[2]*fl[2])*dfac_x; 
  incr[13] = 0.05892556509887893*(15.0*alpha[2]*fl[17]+15.0*alpha[0]*fl[13]+11.61895003862225*alpha[2]*fl[10]+6.708203932499369*alpha[2]*fl[6]+11.61895003862225*alpha[0]*fl[5]+6.708203932499369*alpha[0]*fl[3])*dfac_x; 
  incr[14] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[18]+30.0*alpha[2]*fl[17]+15.0*alpha[0]*fl[14]+23.2379000772445*alpha[2]*fl[10]+13.41640786499874*alpha[2]*fl[6])*dfac_x; 
  incr[15] = -0.03535533905932736*(15.0*alpha[2]*fl[19]+8.660254037844386*alpha[2]*fl[16]+15.0*alpha[0]*fl[15]+8.660254037844387*alpha[0]*fl[9])*dfac_x; 
  incr[16] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[19]+15.0*alpha[0]*fl[16]+25.98076211353316*alpha[2]*fl[15]+15.0*alpha[2]*fl[9])*dfac_x; 
  incr[17] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[18]+75.0*alpha[0]*fl[17]+30.0*alpha[2]*fl[14]+75.00000000000001*alpha[2]*fl[13]+58.09475019311126*alpha[0]*fl[10]+33.54101966249685*alpha[0]*fl[6]+58.09475019311126*alpha[2]*fl[5]+33.54101966249685*alpha[2]*fl[3])*dfac_x; 
  incr[18] = -0.03535533905932736*(15.0*alpha[0]*fl[18]+17.32050807568877*alpha[2]*fl[17]+8.660254037844387*alpha[0]*fl[14]+13.41640786499874*alpha[2]*fl[10]+7.745966692414834*alpha[2]*fl[6])*dfac_x; 
  incr[19] = -0.03535533905932736*(15.0*alpha[0]*fl[19]+8.660254037844387*alpha[0]*fl[16]+15.0*alpha[2]*fl[15]+8.660254037844386*alpha[2]*fl[9])*dfac_x; 

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
  } else { 
  incr[0] = 0.05892556509887893*(6.708203932499369*alpha[2]*fr[11]+6.708203932499369*alpha[0]*fr[7]-5.196152422706631*alpha[2]*fr[4]+3.0*alpha[2]*fr[2]-5.196152422706631*alpha[0]*fr[1]+3.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = -0.1767766952966368*(3.872983346207417*alpha[2]*fr[11]+3.872983346207417*alpha[0]*fr[7]-3.0*alpha[2]*fr[4]+1.732050807568877*alpha[2]*fr[2]-3.0*alpha[0]*fr[1]+1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[12]-33.54101966249684*alpha[0]*fr[11]-13.41640786499874*alpha[2]*fr[8]-33.54101966249685*alpha[2]*fr[7]+25.98076211353316*alpha[0]*fr[4]-15.0*alpha[0]*fr[2]+25.98076211353316*fr[1]*alpha[2]-15.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05892556509887893*(6.708203932499369*alpha[2]*fr[17]+6.708203932499369*alpha[0]*fr[13]-5.196152422706631*alpha[2]*fr[10]+3.0*alpha[2]*fr[6]-5.196152422706631*alpha[0]*fr[5]+3.0*alpha[0]*fr[3])*dfac_x; 
  incr[4] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[12]-19.36491673103708*alpha[0]*fr[11]-7.745966692414834*alpha[2]*fr[8]-19.36491673103709*alpha[2]*fr[7]+15.0*alpha[0]*fr[4]-8.660254037844386*alpha[0]*fr[2]+15.0*fr[1]*alpha[2]-8.660254037844386*fr[0]*alpha[2])*dfac_x; 
  incr[5] = -0.1767766952966368*(3.872983346207417*alpha[2]*fr[17]+3.872983346207417*alpha[0]*fr[13]-3.0*alpha[2]*fr[10]+1.732050807568877*alpha[2]*fr[6]-3.0*alpha[0]*fr[5]+1.732050807568877*alpha[0]*fr[3])*dfac_x; 
  incr[6] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[18]-33.54101966249685*alpha[0]*fr[17]-13.41640786499874*alpha[2]*fr[14]-33.54101966249684*alpha[2]*fr[13]+25.98076211353316*alpha[0]*fr[10]-15.0*alpha[0]*fr[6]+25.98076211353316*alpha[2]*fr[5]-15.0*alpha[2]*fr[3])*dfac_x; 
  incr[7] = 0.05892556509887893*(15.0*alpha[2]*fr[11]+15.0*alpha[0]*fr[7]-11.61895003862225*alpha[2]*fr[4]+6.708203932499369*alpha[2]*fr[2]-11.61895003862225*alpha[0]*fr[1]+6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[8] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[12]-30.0*alpha[2]*fr[11]-15.0*alpha[0]*fr[8]+23.2379000772445*alpha[2]*fr[4]-13.41640786499874*alpha[2]*fr[2])*dfac_x; 
  incr[9] = -0.01178511301977579*(25.98076211353316*alpha[2]*fr[19]-15.0*alpha[2]*fr[16]+25.98076211353316*alpha[0]*fr[15]-15.0*alpha[0]*fr[9])*dfac_x; 
  incr[10] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[18]-19.36491673103709*alpha[0]*fr[17]-7.745966692414834*alpha[2]*fr[14]-19.36491673103708*alpha[2]*fr[13]+15.0*alpha[0]*fr[10]-8.660254037844386*alpha[0]*fr[6]+15.0*alpha[2]*fr[5]-8.660254037844386*alpha[2]*fr[3])*dfac_x; 
  incr[11] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[12]-75.0*alpha[0]*fr[11]-30.0*alpha[2]*fr[8]-75.00000000000001*alpha[2]*fr[7]+58.09475019311126*alpha[0]*fr[4]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*fr[1]*alpha[2]-33.54101966249684*fr[0]*alpha[2])*dfac_x; 
  incr[12] = 0.03535533905932736*(15.0*alpha[0]*fr[12]-17.32050807568877*alpha[2]*fr[11]-8.660254037844387*alpha[0]*fr[8]+13.41640786499874*alpha[2]*fr[4]-7.745966692414834*alpha[2]*fr[2])*dfac_x; 
  incr[13] = 0.05892556509887893*(15.0*alpha[2]*fr[17]+15.0*alpha[0]*fr[13]-11.61895003862225*alpha[2]*fr[10]+6.708203932499369*alpha[2]*fr[6]-11.61895003862225*alpha[0]*fr[5]+6.708203932499369*alpha[0]*fr[3])*dfac_x; 
  incr[14] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[18]-30.0*alpha[2]*fr[17]-15.0*alpha[0]*fr[14]+23.2379000772445*alpha[2]*fr[10]-13.41640786499874*alpha[2]*fr[6])*dfac_x; 
  incr[15] = 0.03535533905932736*(15.0*alpha[2]*fr[19]-8.660254037844386*alpha[2]*fr[16]+15.0*alpha[0]*fr[15]-8.660254037844387*alpha[0]*fr[9])*dfac_x; 
  incr[16] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[19]-15.0*alpha[0]*fr[16]+25.98076211353316*alpha[2]*fr[15]-15.0*alpha[2]*fr[9])*dfac_x; 
  incr[17] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[18]-75.0*alpha[0]*fr[17]-30.0*alpha[2]*fr[14]-75.00000000000001*alpha[2]*fr[13]+58.09475019311126*alpha[0]*fr[10]-33.54101966249685*alpha[0]*fr[6]+58.09475019311126*alpha[2]*fr[5]-33.54101966249685*alpha[2]*fr[3])*dfac_x; 
  incr[18] = 0.03535533905932736*(15.0*alpha[0]*fr[18]-17.32050807568877*alpha[2]*fr[17]-8.660254037844387*alpha[0]*fr[14]+13.41640786499874*alpha[2]*fr[10]-7.745966692414834*alpha[2]*fr[6])*dfac_x; 
  incr[19] = 0.03535533905932736*(15.0*alpha[0]*fr[19]-8.660254037844387*alpha[0]*fr[16]+15.0*alpha[2]*fr[15]-8.660254037844386*alpha[2]*fr[9])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(1.732050807568877*Phi[1]*dfac_x+dApardt[0])*q_)/m_; 

  double alpha[8]; 
  alpha[0] = (-(3.464101615137754*Phi[1]*dfac_x*q_)/m_)-(2.0*dApardt[0]*q_)/m_; 
  alpha[1] = -(2.0*dApardt[1]*q_)/m_; 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[1]*fl[4]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[4]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.1767766952966368*(3.0*alpha[1]*fl[4]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[1]*fl[7]+1.732050807568877*alpha[0]*fl[6]+alpha[1]*fl[5]+alpha[0]*fl[3])*dfac_v; 
  incr[4] = -0.1767766952966368*(3.0*alpha[0]*fl[4]+3.0*alpha[1]*fl[2]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[6]+alpha[0]*fl[5]+alpha[1]*fl[3])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*alpha[1]*fl[7]+3.0*alpha[0]*fl[6]+1.732050807568877*alpha[1]*fl[5]+1.732050807568877*alpha[0]*fl[3])*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fl[7]+3.0*alpha[1]*fl[6]+1.732050807568877*alpha[0]*fl[5]+1.732050807568877*alpha[1]*fl[3])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[1]*fr[4]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[4]+1.732050807568877*alpha[1]*fr[2]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(3.0*alpha[1]*fr[4]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[1]*fr[7]+1.732050807568877*alpha[0]*fr[6]-1.0*alpha[1]*fr[5]-1.0*alpha[0]*fr[3])*dfac_v; 
  incr[4] = 0.1767766952966368*(3.0*alpha[0]*fr[4]+3.0*alpha[1]*fr[2]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[7]+1.732050807568877*alpha[1]*fr[6]-1.0*alpha[0]*fr[5]-1.0*alpha[1]*fr[3])*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*alpha[1]*fr[7]+3.0*alpha[0]*fr[6]-1.732050807568877*alpha[1]*fr[5]-1.732050807568877*alpha[0]*fr[3])*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fr[7]+3.0*alpha[1]*fr[6]-1.732050807568877*alpha[0]*fr[5]-1.732050807568877*alpha[1]*fr[3])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[20]; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(1.732050807568877*Phi[1]*dfac_x+dApardt[0])*q_)/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-(3.061862178478972*Phi[2]*fr[12]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[8]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[8]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-(0.7905694150420948*dApardt[1]*fr[12]*q_)/m_-(0.7905694150420948*dApardt[1]*fl[12]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[11]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[0]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[0]*fl[8]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[7]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[7]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[4]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[2]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[0]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[0]*q_)/m_-1.118033988749895*fr[8]*amax+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Phi[1]*fr[12]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[12]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[11]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[11]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[8]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[8]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[7]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[4]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[2]*fr[12]*q_)/m_-(0.7905694150420948*dApardt[0]*fr[12]*q_)/m_-(0.7071067811865475*dApardt[2]*fl[12]*q_)/m_-(0.7905694150420948*dApardt[0]*fl[12]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[11]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[1]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[1]*fl[8]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[7]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[7]*q_)/m_+(0.5477225575051661*dApardt[2]*fr[4]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[4]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[4]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[2]*q_)/m_-(0.3162277660168379*fr[1]*dApardt[2]*q_)/m_-(0.3162277660168379*fl[1]*dApardt[2]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[1]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[1]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[1]*q_)/m_-1.118033988749895*fr[12]*amax+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[3] = (-(3.061862178478972*Phi[2]*fr[18]*dfac_x*q_)/m_)-(3.061862178478972*Phi[2]*fl[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[14]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[14]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[10]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[3]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[3]*dfac_x*q_)/m_-(0.7905694150420947*dApardt[1]*fr[18]*q_)/m_-(0.7905694150420947*dApardt[1]*fl[18]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[17]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[17]*q_)/m_-(0.7905694150420948*dApardt[0]*fr[14]*q_)/m_-(0.7905694150420948*dApardt[0]*fl[14]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[13]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[13]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[10]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[6]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[5]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[3]*q_)/m_-1.118033988749895*fr[14]*amax+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[5] = (-(1.369306393762915*Phi[1]*fr[18]*dfac_x*q_)/m_)-(1.369306393762915*Phi[1]*fl[18]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[17]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[17]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[14]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[14]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[13]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[13]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[10]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[10]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[3]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[3]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[2]*fr[18]*q_)/m_-(0.7905694150420947*dApardt[0]*fr[18]*q_)/m_-(0.7071067811865475*dApardt[2]*fl[18]*q_)/m_-(0.7905694150420947*dApardt[0]*fl[18]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[17]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[17]*q_)/m_-(0.7905694150420948*dApardt[1]*fr[14]*q_)/m_-(0.7905694150420948*dApardt[1]*fl[14]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[13]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[13]*q_)/m_+(0.5477225575051661*dApardt[2]*fr[10]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[10]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[10]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[6]*q_)/m_-(0.3162277660168379*dApardt[2]*fr[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[5]*q_)/m_-(0.3162277660168379*dApardt[2]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[3]*q_)/m_-1.118033988749895*fr[18]*amax+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 
  Ghat[7] = (-(2.738612787525831*Phi[2]*fr[12]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[12]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[11]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[11]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[7]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[4]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[1]*fr[12]*q_)/m_-(0.7071067811865475*dApardt[1]*fl[12]*q_)/m_+(0.3912303982179757*dApardt[2]*fr[11]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[11]*q_)/m_-(0.3912303982179757*dApardt[2]*fl[11]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[2]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[2]*fl[8]*q_)/m_-(0.2258769757263128*dApardt[2]*fr[7]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[7]*q_)/m_-(0.2258769757263128*dApardt[2]*fl[7]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[7]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[4]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[2]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[2]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[2]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[1]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[1]*q_)/m_+0.8660254037844387*fr[11]*amax+0.8660254037844387*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax; 
  Ghat[9] = (2.371708245126284*Phi[2]*fr[19]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[19]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[16]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[16]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[9]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[9]*dfac_x*q_)/m_+(0.6123724356957944*dApardt[1]*fr[19]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[19]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[16]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[16]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[15]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[9]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[9]*q_)/m_+0.8660254037844387*fr[16]*amax+0.8660254037844387*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[13] = (-(2.738612787525831*Phi[2]*fr[18]*dfac_x*q_)/m_)-(2.738612787525831*Phi[2]*fl[18]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[17]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[17]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[13]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[13]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[10]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[1]*fr[18]*q_)/m_-(0.7071067811865475*dApardt[1]*fl[18]*q_)/m_+(0.3912303982179757*dApardt[2]*fr[17]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[17]*q_)/m_-(0.3912303982179757*dApardt[2]*fl[17]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[17]*q_)/m_-(0.7905694150420947*dApardt[2]*fr[14]*q_)/m_-(0.7905694150420947*dApardt[2]*fl[14]*q_)/m_-(0.2258769757263128*dApardt[2]*fr[13]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[13]*q_)/m_-(0.2258769757263128*dApardt[2]*fl[13]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[13]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[10]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[6]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[5]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[3]*q_)/m_+0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax; 
  Ghat[15] = (1.060660171779821*Phi[1]*fr[19]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[19]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[16]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[16]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[9]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[9]*dfac_x*q_)/m_+(0.5477225575051661*dApardt[2]*fr[19]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[19]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[19]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[19]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[16]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[16]*q_)/m_-(0.3162277660168379*dApardt[2]*fr[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[15]*q_)/m_-(0.3162277660168379*dApardt[2]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[9]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[9]*q_)/m_+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax; 

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
double GyrokineticSurf1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double alpha[8]; 
  alpha[0] = 2.828427124746191*wv; 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_X_P2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[20]; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double alpha[20]; 
  alpha[0] = 2.828427124746191*wv; 
  alpha[2] = 1.632993161855453/dfac_v; 
  if (alpha0>0) { 
  incr[0] = 0.05892556509887893*(6.708203932499369*alpha[2]*fl[11]+6.708203932499369*alpha[0]*fl[7]+5.196152422706631*alpha[2]*fl[4]+3.0*alpha[2]*fl[2]+5.196152422706631*alpha[0]*fl[1]+3.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*(3.872983346207417*alpha[2]*fl[11]+3.872983346207417*alpha[0]*fl[7]+3.0*alpha[2]*fl[4]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[12]+33.54101966249684*alpha[0]*fl[11]+13.41640786499874*alpha[2]*fl[8]+33.54101966249685*alpha[2]*fl[7]+25.98076211353316*alpha[0]*fl[4]+15.0*alpha[0]*fl[2]+25.98076211353316*fl[1]*alpha[2]+15.0*fl[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05892556509887893*(6.708203932499369*alpha[2]*fl[17]+6.708203932499369*alpha[0]*fl[13]+5.196152422706631*alpha[2]*fl[10]+3.0*alpha[2]*fl[6]+5.196152422706631*alpha[0]*fl[5]+3.0*alpha[0]*fl[3])*dfac_x; 
  incr[4] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[12]+19.36491673103708*alpha[0]*fl[11]+7.745966692414834*alpha[2]*fl[8]+19.36491673103709*alpha[2]*fl[7]+15.0*alpha[0]*fl[4]+8.660254037844386*alpha[0]*fl[2]+15.0*fl[1]*alpha[2]+8.660254037844386*fl[0]*alpha[2])*dfac_x; 
  incr[5] = -0.1767766952966368*(3.872983346207417*alpha[2]*fl[17]+3.872983346207417*alpha[0]*fl[13]+3.0*alpha[2]*fl[10]+1.732050807568877*alpha[2]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[3])*dfac_x; 
  incr[6] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[18]+33.54101966249685*alpha[0]*fl[17]+13.41640786499874*alpha[2]*fl[14]+33.54101966249684*alpha[2]*fl[13]+25.98076211353316*alpha[0]*fl[10]+15.0*alpha[0]*fl[6]+25.98076211353316*alpha[2]*fl[5]+15.0*alpha[2]*fl[3])*dfac_x; 
  incr[7] = 0.05892556509887893*(15.0*alpha[2]*fl[11]+15.0*alpha[0]*fl[7]+11.61895003862225*alpha[2]*fl[4]+6.708203932499369*alpha[2]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[8] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[12]+30.0*alpha[2]*fl[11]+15.0*alpha[0]*fl[8]+23.2379000772445*alpha[2]*fl[4]+13.41640786499874*alpha[2]*fl[2])*dfac_x; 
  incr[9] = 0.01178511301977579*(25.98076211353316*alpha[2]*fl[19]+15.0*alpha[2]*fl[16]+25.98076211353316*alpha[0]*fl[15]+15.0*alpha[0]*fl[9])*dfac_x; 
  incr[10] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[18]+19.36491673103709*alpha[0]*fl[17]+7.745966692414834*alpha[2]*fl[14]+19.36491673103708*alpha[2]*fl[13]+15.0*alpha[0]*fl[10]+8.660254037844386*alpha[0]*fl[6]+15.0*alpha[2]*fl[5]+8.660254037844386*alpha[2]*fl[3])*dfac_x; 
  incr[11] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[12]+75.0*alpha[0]*fl[11]+30.0*alpha[2]*fl[8]+75.00000000000001*alpha[2]*fl[7]+58.09475019311126*alpha[0]*fl[4]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*fl[1]*alpha[2]+33.54101966249684*fl[0]*alpha[2])*dfac_x; 
  incr[12] = -0.03535533905932736*(15.0*alpha[0]*fl[12]+17.32050807568877*alpha[2]*fl[11]+8.660254037844387*alpha[0]*fl[8]+13.41640786499874*alpha[2]*fl[4]+7.745966692414834*alpha[2]*fl[2])*dfac_x; 
  incr[13] = 0.05892556509887893*(15.0*alpha[2]*fl[17]+15.0*alpha[0]*fl[13]+11.61895003862225*alpha[2]*fl[10]+6.708203932499369*alpha[2]*fl[6]+11.61895003862225*alpha[0]*fl[5]+6.708203932499369*alpha[0]*fl[3])*dfac_x; 
  incr[14] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[18]+30.0*alpha[2]*fl[17]+15.0*alpha[0]*fl[14]+23.2379000772445*alpha[2]*fl[10]+13.41640786499874*alpha[2]*fl[6])*dfac_x; 
  incr[15] = -0.03535533905932736*(15.0*alpha[2]*fl[19]+8.660254037844386*alpha[2]*fl[16]+15.0*alpha[0]*fl[15]+8.660254037844387*alpha[0]*fl[9])*dfac_x; 
  incr[16] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[19]+15.0*alpha[0]*fl[16]+25.98076211353316*alpha[2]*fl[15]+15.0*alpha[2]*fl[9])*dfac_x; 
  incr[17] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[18]+75.0*alpha[0]*fl[17]+30.0*alpha[2]*fl[14]+75.00000000000001*alpha[2]*fl[13]+58.09475019311126*alpha[0]*fl[10]+33.54101966249685*alpha[0]*fl[6]+58.09475019311126*alpha[2]*fl[5]+33.54101966249685*alpha[2]*fl[3])*dfac_x; 
  incr[18] = -0.03535533905932736*(15.0*alpha[0]*fl[18]+17.32050807568877*alpha[2]*fl[17]+8.660254037844387*alpha[0]*fl[14]+13.41640786499874*alpha[2]*fl[10]+7.745966692414834*alpha[2]*fl[6])*dfac_x; 
  incr[19] = -0.03535533905932736*(15.0*alpha[0]*fl[19]+8.660254037844387*alpha[0]*fl[16]+15.0*alpha[2]*fl[15]+8.660254037844386*alpha[2]*fl[9])*dfac_x; 

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
  } else { 
  incr[0] = 0.05892556509887893*(6.708203932499369*alpha[2]*fr[11]+6.708203932499369*alpha[0]*fr[7]-5.196152422706631*alpha[2]*fr[4]+3.0*alpha[2]*fr[2]-5.196152422706631*alpha[0]*fr[1]+3.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = -0.1767766952966368*(3.872983346207417*alpha[2]*fr[11]+3.872983346207417*alpha[0]*fr[7]-3.0*alpha[2]*fr[4]+1.732050807568877*alpha[2]*fr[2]-3.0*alpha[0]*fr[1]+1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[12]-33.54101966249684*alpha[0]*fr[11]-13.41640786499874*alpha[2]*fr[8]-33.54101966249685*alpha[2]*fr[7]+25.98076211353316*alpha[0]*fr[4]-15.0*alpha[0]*fr[2]+25.98076211353316*fr[1]*alpha[2]-15.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05892556509887893*(6.708203932499369*alpha[2]*fr[17]+6.708203932499369*alpha[0]*fr[13]-5.196152422706631*alpha[2]*fr[10]+3.0*alpha[2]*fr[6]-5.196152422706631*alpha[0]*fr[5]+3.0*alpha[0]*fr[3])*dfac_x; 
  incr[4] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[12]-19.36491673103708*alpha[0]*fr[11]-7.745966692414834*alpha[2]*fr[8]-19.36491673103709*alpha[2]*fr[7]+15.0*alpha[0]*fr[4]-8.660254037844386*alpha[0]*fr[2]+15.0*fr[1]*alpha[2]-8.660254037844386*fr[0]*alpha[2])*dfac_x; 
  incr[5] = -0.1767766952966368*(3.872983346207417*alpha[2]*fr[17]+3.872983346207417*alpha[0]*fr[13]-3.0*alpha[2]*fr[10]+1.732050807568877*alpha[2]*fr[6]-3.0*alpha[0]*fr[5]+1.732050807568877*alpha[0]*fr[3])*dfac_x; 
  incr[6] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[18]-33.54101966249685*alpha[0]*fr[17]-13.41640786499874*alpha[2]*fr[14]-33.54101966249684*alpha[2]*fr[13]+25.98076211353316*alpha[0]*fr[10]-15.0*alpha[0]*fr[6]+25.98076211353316*alpha[2]*fr[5]-15.0*alpha[2]*fr[3])*dfac_x; 
  incr[7] = 0.05892556509887893*(15.0*alpha[2]*fr[11]+15.0*alpha[0]*fr[7]-11.61895003862225*alpha[2]*fr[4]+6.708203932499369*alpha[2]*fr[2]-11.61895003862225*alpha[0]*fr[1]+6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[8] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[12]-30.0*alpha[2]*fr[11]-15.0*alpha[0]*fr[8]+23.2379000772445*alpha[2]*fr[4]-13.41640786499874*alpha[2]*fr[2])*dfac_x; 
  incr[9] = -0.01178511301977579*(25.98076211353316*alpha[2]*fr[19]-15.0*alpha[2]*fr[16]+25.98076211353316*alpha[0]*fr[15]-15.0*alpha[0]*fr[9])*dfac_x; 
  incr[10] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[18]-19.36491673103709*alpha[0]*fr[17]-7.745966692414834*alpha[2]*fr[14]-19.36491673103708*alpha[2]*fr[13]+15.0*alpha[0]*fr[10]-8.660254037844386*alpha[0]*fr[6]+15.0*alpha[2]*fr[5]-8.660254037844386*alpha[2]*fr[3])*dfac_x; 
  incr[11] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[12]-75.0*alpha[0]*fr[11]-30.0*alpha[2]*fr[8]-75.00000000000001*alpha[2]*fr[7]+58.09475019311126*alpha[0]*fr[4]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*fr[1]*alpha[2]-33.54101966249684*fr[0]*alpha[2])*dfac_x; 
  incr[12] = 0.03535533905932736*(15.0*alpha[0]*fr[12]-17.32050807568877*alpha[2]*fr[11]-8.660254037844387*alpha[0]*fr[8]+13.41640786499874*alpha[2]*fr[4]-7.745966692414834*alpha[2]*fr[2])*dfac_x; 
  incr[13] = 0.05892556509887893*(15.0*alpha[2]*fr[17]+15.0*alpha[0]*fr[13]-11.61895003862225*alpha[2]*fr[10]+6.708203932499369*alpha[2]*fr[6]-11.61895003862225*alpha[0]*fr[5]+6.708203932499369*alpha[0]*fr[3])*dfac_x; 
  incr[14] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[18]-30.0*alpha[2]*fr[17]-15.0*alpha[0]*fr[14]+23.2379000772445*alpha[2]*fr[10]-13.41640786499874*alpha[2]*fr[6])*dfac_x; 
  incr[15] = 0.03535533905932736*(15.0*alpha[2]*fr[19]-8.660254037844386*alpha[2]*fr[16]+15.0*alpha[0]*fr[15]-8.660254037844387*alpha[0]*fr[9])*dfac_x; 
  incr[16] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[19]-15.0*alpha[0]*fr[16]+25.98076211353316*alpha[2]*fr[15]-15.0*alpha[2]*fr[9])*dfac_x; 
  incr[17] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[18]-75.0*alpha[0]*fr[17]-30.0*alpha[2]*fr[14]-75.00000000000001*alpha[2]*fr[13]+58.09475019311126*alpha[0]*fr[10]-33.54101966249685*alpha[0]*fr[6]+58.09475019311126*alpha[2]*fr[5]-33.54101966249685*alpha[2]*fr[3])*dfac_x; 
  incr[18] = 0.03535533905932736*(15.0*alpha[0]*fr[18]-17.32050807568877*alpha[2]*fr[17]-8.660254037844387*alpha[0]*fr[14]+13.41640786499874*alpha[2]*fr[10]-7.745966692414834*alpha[2]*fr[6])*dfac_x; 
  incr[19] = 0.03535533905932736*(15.0*alpha[0]*fr[19]-8.660254037844387*alpha[0]*fr[16]+15.0*alpha[2]*fr[15]-8.660254037844386*alpha[2]*fr[9])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(1.732050807568877*Bmag[1]*dfac_x*wm+1.732050807568877*Phi[1]*dfac_x*q_))/m_; 

  double alpha[8]; 
  alpha[0] = (-(3.464101615137754*Bmag[1]*dfac_x*wm)/m_)-(3.464101615137754*Phi[1]*dfac_x*q_)/m_; 
  alpha[3] = -(2.0*Bmag[1]*dfac_x)/(dfac_m*m_); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[3]*fl[6]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*alpha[3]*fl[7]+alpha[3]*fl[5]+1.732050807568877*alpha[0]*fl[4]+alpha[0]*fl[1])*dfac_v; 
  incr[2] = -0.1767766952966368*(3.0*alpha[3]*fl[6]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[6]+alpha[0]*fl[3]+1.732050807568877*fl[2]*alpha[3]+fl[0]*alpha[3])*dfac_v; 
  incr[4] = -0.1767766952966368*(3.0*alpha[3]*fl[7]+1.732050807568877*alpha[3]*fl[5]+3.0*alpha[0]*fl[4]+1.732050807568877*alpha[0]*fl[1])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[7]+alpha[0]*fl[5]+1.732050807568877*alpha[3]*fl[4]+fl[1]*alpha[3])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*alpha[0]*fl[6]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fl[7]+1.732050807568877*alpha[0]*fl[5]+3.0*alpha[3]*fl[4]+1.732050807568877*fl[1]*alpha[3])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[3]*fr[6]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*alpha[3]*fr[7]-1.0*alpha[3]*fr[5]+1.732050807568877*alpha[0]*fr[4]-1.0*alpha[0]*fr[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(3.0*alpha[3]*fr[6]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[6]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_v; 
  incr[4] = 0.1767766952966368*(3.0*alpha[3]*fr[7]-1.732050807568877*alpha[3]*fr[5]+3.0*alpha[0]*fr[4]-1.732050807568877*alpha[0]*fr[1])*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[7]-1.0*alpha[0]*fr[5]+1.732050807568877*alpha[3]*fr[4]-1.0*fr[1]*alpha[3])*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*alpha[0]*fr[6]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fr[7]-1.732050807568877*alpha[0]*fr[5]+3.0*alpha[3]*fr[4]-1.732050807568877*fr[1]*alpha[3])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[20]; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(1.732050807568877*Bmag[1]*dfac_x*wm+1.732050807568877*Phi[1]*dfac_x*q_))/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-(3.061862178478972*Bmag[2]*fr[12]*dfac_x*wm)/m_)-(3.061862178478972*Bmag[2]*fl[12]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fr[8]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fl[8]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[4]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[4]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[2]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[2]*dfac_x*wm)/m_-(1.369306393762915*fr[1]*Bmag[2]*dfac_x*wm)/m_-(1.369306393762915*fl[1]*Bmag[2]*dfac_x*wm)/m_-(0.6123724356957944*fr[0]*Bmag[1]*dfac_x*wm)/m_-(0.6123724356957944*fl[0]*Bmag[1]*dfac_x*wm)/m_-(3.061862178478972*Phi[2]*fr[12]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[8]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[8]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-(1.767766952966368*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fr[14]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fl[14]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[6]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[6]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[3]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[3]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[8]*amax+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Bmag[1]*fr[12]*dfac_x*wm)/m_)-(1.369306393762915*Bmag[1]*fl[12]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[11]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[11]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fr[8]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fl[8]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[7]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[7]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[4]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[4]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[2]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[2]*dfac_x*wm)/m_-(1.369306393762915*fr[0]*Bmag[2]*dfac_x*wm)/m_-(1.369306393762915*fl[0]*Bmag[2]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[1]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[1]*dfac_x*wm)/m_-(1.369306393762915*Phi[1]*fr[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[12]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[11]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[11]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[8]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[8]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[7]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[4]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-(0.7905694150420947*Bmag[1]*fr[18]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fl[18]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[17]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[17]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fr[14]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fl[14]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[13]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[13]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[10]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[10]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[6]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[6]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[5]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[5]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fr[3]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fl[3]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[12]*amax+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[3] = (-(3.061862178478972*Bmag[2]*fr[18]*dfac_x*wm)/m_)-(3.061862178478972*Bmag[2]*fl[18]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fr[14]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fl[14]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[10]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[10]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[6]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[6]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[5]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[5]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[3]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[3]*dfac_x*wm)/m_-(3.061862178478972*Phi[2]*fr[18]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[14]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[14]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[10]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[3]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[3]*dfac_x*q_)/m_+(1.224744871391589*Bmag[2]*fr[19]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[19]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[16]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[16]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[15]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[15]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fr[12]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fl[12]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[9]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[9]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fr[8]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fl[8]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[4]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[4]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[2]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fr[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fl[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*fr[0]*Bmag[1]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*fl[0]*Bmag[1]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[14]*amax+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[5] = (-(1.369306393762915*Bmag[1]*fr[18]*dfac_x*wm)/m_)-(1.369306393762915*Bmag[1]*fl[18]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[17]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[17]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fr[14]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fl[14]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[13]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[13]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[10]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[10]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[6]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[6]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[5]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[5]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[3]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[3]*dfac_x*wm)/m_-(1.369306393762915*Phi[1]*fr[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[18]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[17]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[17]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[14]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[14]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[13]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[13]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[10]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[10]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[3]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[3]*dfac_x*q_)/m_+(0.5477225575051661*Bmag[1]*fr[19]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[19]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[16]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[16]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[15]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[15]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fr[12]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fl[12]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[11]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[11]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[9]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[9]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fr[8]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fl[8]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[7]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[7]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[4]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[4]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[2]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fr[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fl[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[1]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[1]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[18]*amax+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 
  Ghat[7] = (-(2.738612787525831*Bmag[2]*fr[12]*dfac_x*wm)/m_)-(2.738612787525831*Bmag[2]*fl[12]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[11]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[11]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[7]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[7]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[4]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[4]*dfac_x*wm)/m_-(1.224744871391589*fr[1]*Bmag[2]*dfac_x*wm)/m_-(1.224744871391589*fl[1]*Bmag[2]*dfac_x*wm)/m_-(2.738612787525831*Phi[2]*fr[12]*dfac_x*q_)/m_-(2.738612787525831*Phi[2]*fl[12]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[11]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[11]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[7]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[4]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_-(1.58113883008419*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[17]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[17]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[13]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[13]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[11]*amax+0.8660254037844387*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax; 
  Ghat[9] = (2.371708245126284*Bmag[2]*fr[19]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[19]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[16]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[16]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[15]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[15]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[9]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[9]*dfac_x*wm)/m_+(2.371708245126284*Phi[2]*fr[19]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[19]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[16]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[16]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[9]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[9]*dfac_x*q_)/m_-(1.58113883008419*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fr[14]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fl[14]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[6]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[6]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[3]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[3]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[16]*amax+0.8660254037844387*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[13] = (-(2.738612787525831*Bmag[2]*fr[18]*dfac_x*wm)/m_)-(2.738612787525831*Bmag[2]*fl[18]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[17]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[17]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[13]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[13]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[10]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[10]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[5]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[5]*dfac_x*wm)/m_-(2.738612787525831*Phi[2]*fr[18]*dfac_x*q_)/m_-(2.738612787525831*Phi[2]*fl[18]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[17]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[17]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[13]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[13]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[10]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[5]*dfac_x*q_)/m_+(1.095445115010332*Bmag[2]*fr[19]*dfac_x)/(dfac_m*m_)-(1.095445115010332*Bmag[2]*fl[19]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fr[15]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fl[15]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fr[12]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[12]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[11]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[11]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[7]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[7]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[4]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[4]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*fr[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*fl[1]*Bmag[2]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax; 
  Ghat[15] = (1.060660171779821*Bmag[1]*fr[19]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[19]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[16]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[16]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[15]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[15]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[9]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[9]*dfac_x*wm)/m_+(1.060660171779821*Phi[1]*fr[19]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[19]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[16]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[16]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[9]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[9]*dfac_x*q_)/m_-(0.7071067811865475*Bmag[1]*fr[18]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fl[18]*dfac_x)/(dfac_m*m_)+(1.095445115010332*Bmag[2]*fr[17]*dfac_x)/(dfac_m*m_)-(1.095445115010332*Bmag[2]*fl[17]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fr[14]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[14]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fr[13]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fl[13]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[10]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[10]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[6]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[6]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[5]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[3]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[3]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax; 

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
double EmGyrokineticSurf1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double alpha[8]; 
  alpha[0] = 2.828427124746191*wv; 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_X_P2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[20]; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 1.414213562373095*wv; 

  double alpha[20]; 
  alpha[0] = 2.828427124746191*wv; 
  alpha[2] = 1.632993161855453/dfac_v; 
  if (alpha0>0) { 
  incr[0] = 0.05892556509887893*(6.708203932499369*alpha[2]*fl[11]+6.708203932499369*alpha[0]*fl[7]+5.196152422706631*alpha[2]*fl[4]+3.0*alpha[2]*fl[2]+5.196152422706631*alpha[0]*fl[1]+3.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*(3.872983346207417*alpha[2]*fl[11]+3.872983346207417*alpha[0]*fl[7]+3.0*alpha[2]*fl[4]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[12]+33.54101966249684*alpha[0]*fl[11]+13.41640786499874*alpha[2]*fl[8]+33.54101966249685*alpha[2]*fl[7]+25.98076211353316*alpha[0]*fl[4]+15.0*alpha[0]*fl[2]+25.98076211353316*fl[1]*alpha[2]+15.0*fl[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05892556509887893*(6.708203932499369*alpha[2]*fl[17]+6.708203932499369*alpha[0]*fl[13]+5.196152422706631*alpha[2]*fl[10]+3.0*alpha[2]*fl[6]+5.196152422706631*alpha[0]*fl[5]+3.0*alpha[0]*fl[3])*dfac_x; 
  incr[4] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[12]+19.36491673103708*alpha[0]*fl[11]+7.745966692414834*alpha[2]*fl[8]+19.36491673103709*alpha[2]*fl[7]+15.0*alpha[0]*fl[4]+8.660254037844386*alpha[0]*fl[2]+15.0*fl[1]*alpha[2]+8.660254037844386*fl[0]*alpha[2])*dfac_x; 
  incr[5] = -0.1767766952966368*(3.872983346207417*alpha[2]*fl[17]+3.872983346207417*alpha[0]*fl[13]+3.0*alpha[2]*fl[10]+1.732050807568877*alpha[2]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[3])*dfac_x; 
  incr[6] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[18]+33.54101966249685*alpha[0]*fl[17]+13.41640786499874*alpha[2]*fl[14]+33.54101966249684*alpha[2]*fl[13]+25.98076211353316*alpha[0]*fl[10]+15.0*alpha[0]*fl[6]+25.98076211353316*alpha[2]*fl[5]+15.0*alpha[2]*fl[3])*dfac_x; 
  incr[7] = 0.05892556509887893*(15.0*alpha[2]*fl[11]+15.0*alpha[0]*fl[7]+11.61895003862225*alpha[2]*fl[4]+6.708203932499369*alpha[2]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[8] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[12]+30.0*alpha[2]*fl[11]+15.0*alpha[0]*fl[8]+23.2379000772445*alpha[2]*fl[4]+13.41640786499874*alpha[2]*fl[2])*dfac_x; 
  incr[9] = 0.01178511301977579*(25.98076211353316*alpha[2]*fl[19]+15.0*alpha[2]*fl[16]+25.98076211353316*alpha[0]*fl[15]+15.0*alpha[0]*fl[9])*dfac_x; 
  incr[10] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[18]+19.36491673103709*alpha[0]*fl[17]+7.745966692414834*alpha[2]*fl[14]+19.36491673103708*alpha[2]*fl[13]+15.0*alpha[0]*fl[10]+8.660254037844386*alpha[0]*fl[6]+15.0*alpha[2]*fl[5]+8.660254037844386*alpha[2]*fl[3])*dfac_x; 
  incr[11] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[12]+75.0*alpha[0]*fl[11]+30.0*alpha[2]*fl[8]+75.00000000000001*alpha[2]*fl[7]+58.09475019311126*alpha[0]*fl[4]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*fl[1]*alpha[2]+33.54101966249684*fl[0]*alpha[2])*dfac_x; 
  incr[12] = -0.03535533905932736*(15.0*alpha[0]*fl[12]+17.32050807568877*alpha[2]*fl[11]+8.660254037844387*alpha[0]*fl[8]+13.41640786499874*alpha[2]*fl[4]+7.745966692414834*alpha[2]*fl[2])*dfac_x; 
  incr[13] = 0.05892556509887893*(15.0*alpha[2]*fl[17]+15.0*alpha[0]*fl[13]+11.61895003862225*alpha[2]*fl[10]+6.708203932499369*alpha[2]*fl[6]+11.61895003862225*alpha[0]*fl[5]+6.708203932499369*alpha[0]*fl[3])*dfac_x; 
  incr[14] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[18]+30.0*alpha[2]*fl[17]+15.0*alpha[0]*fl[14]+23.2379000772445*alpha[2]*fl[10]+13.41640786499874*alpha[2]*fl[6])*dfac_x; 
  incr[15] = -0.03535533905932736*(15.0*alpha[2]*fl[19]+8.660254037844386*alpha[2]*fl[16]+15.0*alpha[0]*fl[15]+8.660254037844387*alpha[0]*fl[9])*dfac_x; 
  incr[16] = 0.01178511301977579*(25.98076211353316*alpha[0]*fl[19]+15.0*alpha[0]*fl[16]+25.98076211353316*alpha[2]*fl[15]+15.0*alpha[2]*fl[9])*dfac_x; 
  incr[17] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[18]+75.0*alpha[0]*fl[17]+30.0*alpha[2]*fl[14]+75.00000000000001*alpha[2]*fl[13]+58.09475019311126*alpha[0]*fl[10]+33.54101966249685*alpha[0]*fl[6]+58.09475019311126*alpha[2]*fl[5]+33.54101966249685*alpha[2]*fl[3])*dfac_x; 
  incr[18] = -0.03535533905932736*(15.0*alpha[0]*fl[18]+17.32050807568877*alpha[2]*fl[17]+8.660254037844387*alpha[0]*fl[14]+13.41640786499874*alpha[2]*fl[10]+7.745966692414834*alpha[2]*fl[6])*dfac_x; 
  incr[19] = -0.03535533905932736*(15.0*alpha[0]*fl[19]+8.660254037844387*alpha[0]*fl[16]+15.0*alpha[2]*fl[15]+8.660254037844386*alpha[2]*fl[9])*dfac_x; 

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
  } else { 
  incr[0] = 0.05892556509887893*(6.708203932499369*alpha[2]*fr[11]+6.708203932499369*alpha[0]*fr[7]-5.196152422706631*alpha[2]*fr[4]+3.0*alpha[2]*fr[2]-5.196152422706631*alpha[0]*fr[1]+3.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = -0.1767766952966368*(3.872983346207417*alpha[2]*fr[11]+3.872983346207417*alpha[0]*fr[7]-3.0*alpha[2]*fr[4]+1.732050807568877*alpha[2]*fr[2]-3.0*alpha[0]*fr[1]+1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[12]-33.54101966249684*alpha[0]*fr[11]-13.41640786499874*alpha[2]*fr[8]-33.54101966249685*alpha[2]*fr[7]+25.98076211353316*alpha[0]*fr[4]-15.0*alpha[0]*fr[2]+25.98076211353316*fr[1]*alpha[2]-15.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05892556509887893*(6.708203932499369*alpha[2]*fr[17]+6.708203932499369*alpha[0]*fr[13]-5.196152422706631*alpha[2]*fr[10]+3.0*alpha[2]*fr[6]-5.196152422706631*alpha[0]*fr[5]+3.0*alpha[0]*fr[3])*dfac_x; 
  incr[4] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[12]-19.36491673103708*alpha[0]*fr[11]-7.745966692414834*alpha[2]*fr[8]-19.36491673103709*alpha[2]*fr[7]+15.0*alpha[0]*fr[4]-8.660254037844386*alpha[0]*fr[2]+15.0*fr[1]*alpha[2]-8.660254037844386*fr[0]*alpha[2])*dfac_x; 
  incr[5] = -0.1767766952966368*(3.872983346207417*alpha[2]*fr[17]+3.872983346207417*alpha[0]*fr[13]-3.0*alpha[2]*fr[10]+1.732050807568877*alpha[2]*fr[6]-3.0*alpha[0]*fr[5]+1.732050807568877*alpha[0]*fr[3])*dfac_x; 
  incr[6] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[18]-33.54101966249685*alpha[0]*fr[17]-13.41640786499874*alpha[2]*fr[14]-33.54101966249684*alpha[2]*fr[13]+25.98076211353316*alpha[0]*fr[10]-15.0*alpha[0]*fr[6]+25.98076211353316*alpha[2]*fr[5]-15.0*alpha[2]*fr[3])*dfac_x; 
  incr[7] = 0.05892556509887893*(15.0*alpha[2]*fr[11]+15.0*alpha[0]*fr[7]-11.61895003862225*alpha[2]*fr[4]+6.708203932499369*alpha[2]*fr[2]-11.61895003862225*alpha[0]*fr[1]+6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[8] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[12]-30.0*alpha[2]*fr[11]-15.0*alpha[0]*fr[8]+23.2379000772445*alpha[2]*fr[4]-13.41640786499874*alpha[2]*fr[2])*dfac_x; 
  incr[9] = -0.01178511301977579*(25.98076211353316*alpha[2]*fr[19]-15.0*alpha[2]*fr[16]+25.98076211353316*alpha[0]*fr[15]-15.0*alpha[0]*fr[9])*dfac_x; 
  incr[10] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[18]-19.36491673103709*alpha[0]*fr[17]-7.745966692414834*alpha[2]*fr[14]-19.36491673103708*alpha[2]*fr[13]+15.0*alpha[0]*fr[10]-8.660254037844386*alpha[0]*fr[6]+15.0*alpha[2]*fr[5]-8.660254037844386*alpha[2]*fr[3])*dfac_x; 
  incr[11] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[12]-75.0*alpha[0]*fr[11]-30.0*alpha[2]*fr[8]-75.00000000000001*alpha[2]*fr[7]+58.09475019311126*alpha[0]*fr[4]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*fr[1]*alpha[2]-33.54101966249684*fr[0]*alpha[2])*dfac_x; 
  incr[12] = 0.03535533905932736*(15.0*alpha[0]*fr[12]-17.32050807568877*alpha[2]*fr[11]-8.660254037844387*alpha[0]*fr[8]+13.41640786499874*alpha[2]*fr[4]-7.745966692414834*alpha[2]*fr[2])*dfac_x; 
  incr[13] = 0.05892556509887893*(15.0*alpha[2]*fr[17]+15.0*alpha[0]*fr[13]-11.61895003862225*alpha[2]*fr[10]+6.708203932499369*alpha[2]*fr[6]-11.61895003862225*alpha[0]*fr[5]+6.708203932499369*alpha[0]*fr[3])*dfac_x; 
  incr[14] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[18]-30.0*alpha[2]*fr[17]-15.0*alpha[0]*fr[14]+23.2379000772445*alpha[2]*fr[10]-13.41640786499874*alpha[2]*fr[6])*dfac_x; 
  incr[15] = 0.03535533905932736*(15.0*alpha[2]*fr[19]-8.660254037844386*alpha[2]*fr[16]+15.0*alpha[0]*fr[15]-8.660254037844387*alpha[0]*fr[9])*dfac_x; 
  incr[16] = -0.01178511301977579*(25.98076211353316*alpha[0]*fr[19]-15.0*alpha[0]*fr[16]+25.98076211353316*alpha[2]*fr[15]-15.0*alpha[2]*fr[9])*dfac_x; 
  incr[17] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[18]-75.0*alpha[0]*fr[17]-30.0*alpha[2]*fr[14]-75.00000000000001*alpha[2]*fr[13]+58.09475019311126*alpha[0]*fr[10]-33.54101966249685*alpha[0]*fr[6]+58.09475019311126*alpha[2]*fr[5]-33.54101966249685*alpha[2]*fr[3])*dfac_x; 
  incr[18] = 0.03535533905932736*(15.0*alpha[0]*fr[18]-17.32050807568877*alpha[2]*fr[17]-8.660254037844387*alpha[0]*fr[14]+13.41640786499874*alpha[2]*fr[10]-7.745966692414834*alpha[2]*fr[6])*dfac_x; 
  incr[19] = 0.03535533905932736*(15.0*alpha[0]*fr[19]-8.660254037844387*alpha[0]*fr[16]+15.0*alpha[2]*fr[15]-8.660254037844386*alpha[2]*fr[9])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(1.732050807568877*Bmag[1]*dfac_x*wm+(1.732050807568877*Phi[1]*dfac_x+dApardt[0])*q_))/m_; 

  double alpha[8]; 
  alpha[0] = (-(3.464101615137754*Bmag[1]*dfac_x*wm)/m_)-(3.464101615137754*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[0]*q_)/m_; 
  alpha[1] = -(2.0*dApardt[1]*q_)/m_; 
  alpha[3] = -(2.0*Bmag[1]*dfac_x)/(dfac_m*m_); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[3]*fl[6]+1.732050807568877*alpha[1]*fl[4]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*alpha[3]*fl[7]+alpha[3]*fl[5]+1.732050807568877*alpha[0]*fl[4]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.1767766952966368*(3.0*alpha[3]*fl[6]+3.0*alpha[1]*fl[4]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[1]*fl[7]+1.732050807568877*alpha[0]*fl[6]+alpha[1]*fl[5]+alpha[0]*fl[3]+1.732050807568877*fl[2]*alpha[3]+fl[0]*alpha[3])*dfac_v; 
  incr[4] = -0.1767766952966368*(3.0*alpha[3]*fl[7]+1.732050807568877*alpha[3]*fl[5]+3.0*alpha[0]*fl[4]+3.0*alpha[1]*fl[2]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[6]+alpha[0]*fl[5]+1.732050807568877*alpha[3]*fl[4]+alpha[1]*fl[3]+fl[1]*alpha[3])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*alpha[1]*fl[7]+3.0*alpha[0]*fl[6]+1.732050807568877*alpha[1]*fl[5]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fl[7]+3.0*alpha[1]*fl[6]+1.732050807568877*alpha[0]*fl[5]+3.0*alpha[3]*fl[4]+1.732050807568877*alpha[1]*fl[3]+1.732050807568877*fl[1]*alpha[3])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[3]*fr[6]+1.732050807568877*alpha[1]*fr[4]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*alpha[3]*fr[7]-1.0*alpha[3]*fr[5]+1.732050807568877*alpha[0]*fr[4]+1.732050807568877*alpha[1]*fr[2]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(3.0*alpha[3]*fr[6]+3.0*alpha[1]*fr[4]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[1]*fr[7]+1.732050807568877*alpha[0]*fr[6]-1.0*alpha[1]*fr[5]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_v; 
  incr[4] = 0.1767766952966368*(3.0*alpha[3]*fr[7]-1.732050807568877*alpha[3]*fr[5]+3.0*alpha[0]*fr[4]+3.0*alpha[1]*fr[2]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[7]+1.732050807568877*alpha[1]*fr[6]-1.0*alpha[0]*fr[5]+1.732050807568877*alpha[3]*fr[4]-1.0*alpha[1]*fr[3]-1.0*fr[1]*alpha[3])*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*alpha[1]*fr[7]+3.0*alpha[0]*fr[6]-1.732050807568877*alpha[1]*fr[5]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fr[7]+3.0*alpha[1]*fr[6]-1.732050807568877*alpha[0]*fr[5]+3.0*alpha[3]*fr[4]-1.732050807568877*alpha[1]*fr[3]-1.732050807568877*fr[1]*alpha[3])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x2vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[20]; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(1.0*(1.732050807568877*Bmag[1]*dfac_x*wm+(1.732050807568877*Phi[1]*dfac_x+dApardt[0])*q_))/m_; 

  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-(3.061862178478972*Bmag[2]*fr[12]*dfac_x*wm)/m_)-(3.061862178478972*Bmag[2]*fl[12]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fr[8]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fl[8]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[4]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[4]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[2]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[2]*dfac_x*wm)/m_-(1.369306393762915*fr[1]*Bmag[2]*dfac_x*wm)/m_-(1.369306393762915*fl[1]*Bmag[2]*dfac_x*wm)/m_-(0.6123724356957944*fr[0]*Bmag[1]*dfac_x*wm)/m_-(0.6123724356957944*fl[0]*Bmag[1]*dfac_x*wm)/m_-(3.061862178478972*Phi[2]*fr[12]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[8]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[8]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[4]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[2]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*fr[0]*Phi[1]*dfac_x*q_)/m_-(0.6123724356957944*fl[0]*Phi[1]*dfac_x*q_)/m_-(0.7905694150420948*dApardt[1]*fr[12]*q_)/m_-(0.7905694150420948*dApardt[1]*fl[12]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[11]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[0]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[0]*fl[8]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[7]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[7]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[4]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[2]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[0]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[0]*q_)/m_-(1.767766952966368*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fr[14]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fl[14]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[6]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[6]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[3]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[3]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[8]*amax+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax; 
  Ghat[1] = (-(1.369306393762915*Bmag[1]*fr[12]*dfac_x*wm)/m_)-(1.369306393762915*Bmag[1]*fl[12]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[11]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[11]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fr[8]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fl[8]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[7]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[7]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[4]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[4]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[2]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[2]*dfac_x*wm)/m_-(1.369306393762915*fr[0]*Bmag[2]*dfac_x*wm)/m_-(1.369306393762915*fl[0]*Bmag[2]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[1]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[1]*dfac_x*wm)/m_-(1.369306393762915*Phi[1]*fr[12]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[12]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[11]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[11]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[8]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[8]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[7]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[7]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[4]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[4]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[2]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[2]*dfac_x*q_)/m_-(1.369306393762915*fr[0]*Phi[2]*dfac_x*q_)/m_-(1.369306393762915*fl[0]*Phi[2]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[1]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[1]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[2]*fr[12]*q_)/m_-(0.7905694150420948*dApardt[0]*fr[12]*q_)/m_-(0.7071067811865475*dApardt[2]*fl[12]*q_)/m_-(0.7905694150420948*dApardt[0]*fl[12]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[11]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[1]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[1]*fl[8]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[7]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[7]*q_)/m_+(0.5477225575051661*dApardt[2]*fr[4]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[4]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[4]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[2]*q_)/m_-(0.3162277660168379*fr[1]*dApardt[2]*q_)/m_-(0.3162277660168379*fl[1]*dApardt[2]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[1]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[1]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[1]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[1]*q_)/m_-(0.7905694150420947*Bmag[1]*fr[18]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fl[18]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[17]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[17]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fr[14]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fl[14]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[13]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[13]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[10]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[10]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[6]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[6]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[5]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[5]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fr[3]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[2]*fl[3]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[12]*amax+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax; 
  Ghat[3] = (-(3.061862178478972*Bmag[2]*fr[18]*dfac_x*wm)/m_)-(3.061862178478972*Bmag[2]*fl[18]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fr[14]*dfac_x*wm)/m_-(1.369306393762915*Bmag[1]*fl[14]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[10]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[10]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[6]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[6]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[5]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[5]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[3]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[3]*dfac_x*wm)/m_-(3.061862178478972*Phi[2]*fr[18]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fr[14]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[14]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[10]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[6]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[6]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[3]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[3]*dfac_x*q_)/m_-(0.7905694150420947*dApardt[1]*fr[18]*q_)/m_-(0.7905694150420947*dApardt[1]*fl[18]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[17]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[17]*q_)/m_-(0.7905694150420948*dApardt[0]*fr[14]*q_)/m_-(0.7905694150420948*dApardt[0]*fl[14]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[13]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[13]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[10]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[6]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[5]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[3]*q_)/m_+(1.224744871391589*Bmag[2]*fr[19]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[19]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[16]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[16]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[15]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[15]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fr[12]*dfac_x)/(dfac_m*m_)-(1.767766952966369*Bmag[2]*fl[12]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[9]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[9]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fr[8]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*Bmag[1]*fl[8]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[4]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[4]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[2]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fr[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fl[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*fr[0]*Bmag[1]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*fl[0]*Bmag[1]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[14]*amax+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax; 
  Ghat[5] = (-(1.369306393762915*Bmag[1]*fr[18]*dfac_x*wm)/m_)-(1.369306393762915*Bmag[1]*fl[18]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[17]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[17]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fr[14]*dfac_x*wm)/m_-(3.061862178478972*Bmag[2]*fl[14]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[13]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[13]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[10]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[10]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[6]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[6]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[5]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[5]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[3]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[3]*dfac_x*wm)/m_-(1.369306393762915*Phi[1]*fr[18]*dfac_x*q_)/m_-(1.369306393762915*Phi[1]*fl[18]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[17]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[17]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fr[14]*dfac_x*q_)/m_-(3.061862178478972*Phi[2]*fl[14]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[13]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[13]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[10]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[10]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[6]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[6]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[5]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[5]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[3]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[3]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[2]*fr[18]*q_)/m_-(0.7905694150420947*dApardt[0]*fr[18]*q_)/m_-(0.7071067811865475*dApardt[2]*fl[18]*q_)/m_-(0.7905694150420947*dApardt[0]*fl[18]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[17]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[17]*q_)/m_-(0.7905694150420948*dApardt[1]*fr[14]*q_)/m_-(0.7905694150420948*dApardt[1]*fl[14]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[13]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[13]*q_)/m_+(0.5477225575051661*dApardt[2]*fr[10]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[10]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[10]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[6]*q_)/m_-(0.3162277660168379*dApardt[2]*fr[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[5]*q_)/m_-(0.3162277660168379*dApardt[2]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[3]*q_)/m_+(0.5477225575051661*Bmag[1]*fr[19]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[19]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[16]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[16]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[15]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[15]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fr[12]*dfac_x)/(dfac_m*m_)-(0.7905694150420948*Bmag[1]*fl[12]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[11]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[11]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[9]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[9]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fr[8]*dfac_x)/(dfac_m*m_)-(1.767766952966368*Bmag[2]*fl[8]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[7]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[7]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[4]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[4]*dfac_x)/(dfac_m*m_)+(1.369306393762915*Bmag[2]*fr[2]*dfac_x)/(dfac_m*m_)-(1.369306393762915*Bmag[2]*fl[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fr[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7905694150420947*fl[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[1]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[1]*dfac_x)/(dfac_m*m_)-1.118033988749895*fr[18]*amax+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax; 
  Ghat[7] = (-(2.738612787525831*Bmag[2]*fr[12]*dfac_x*wm)/m_)-(2.738612787525831*Bmag[2]*fl[12]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[11]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[11]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[7]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[7]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[4]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[4]*dfac_x*wm)/m_-(1.224744871391589*fr[1]*Bmag[2]*dfac_x*wm)/m_-(1.224744871391589*fl[1]*Bmag[2]*dfac_x*wm)/m_-(2.738612787525831*Phi[2]*fr[12]*dfac_x*q_)/m_-(2.738612787525831*Phi[2]*fl[12]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[11]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[11]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[7]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[7]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[4]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[4]*dfac_x*q_)/m_-(1.224744871391589*fr[1]*Phi[2]*dfac_x*q_)/m_-(1.224744871391589*fl[1]*Phi[2]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[1]*fr[12]*q_)/m_-(0.7071067811865475*dApardt[1]*fl[12]*q_)/m_+(0.3912303982179757*dApardt[2]*fr[11]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[11]*q_)/m_-(0.3912303982179757*dApardt[2]*fl[11]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[11]*q_)/m_-(0.7905694150420947*dApardt[2]*fr[8]*q_)/m_-(0.7905694150420947*dApardt[2]*fl[8]*q_)/m_-(0.2258769757263128*dApardt[2]*fr[7]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[7]*q_)/m_-(0.2258769757263128*dApardt[2]*fl[7]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[7]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[4]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[4]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[2]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[2]*q_)/m_-(0.3535533905932737*fr[0]*dApardt[2]*q_)/m_-(0.3535533905932737*fl[0]*dApardt[2]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[1]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[1]*q_)/m_-(1.58113883008419*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[17]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[17]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[13]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[13]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[11]*amax+0.8660254037844387*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax; 
  Ghat[9] = (2.371708245126284*Bmag[2]*fr[19]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[19]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[16]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[16]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[15]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[15]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[9]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[9]*dfac_x*wm)/m_+(2.371708245126284*Phi[2]*fr[19]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[19]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[16]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[16]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[9]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[9]*dfac_x*q_)/m_+(0.6123724356957944*dApardt[1]*fr[19]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[19]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[16]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[16]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[15]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[9]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[9]*q_)/m_-(1.58113883008419*Bmag[2]*fr[18]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[18]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fr[14]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fl[14]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[10]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[10]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[6]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[6]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[5]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[3]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[3]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[16]*amax+0.8660254037844387*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax; 
  Ghat[13] = (-(2.738612787525831*Bmag[2]*fr[18]*dfac_x*wm)/m_)-(2.738612787525831*Bmag[2]*fl[18]*dfac_x*wm)/m_+(1.060660171779821*Bmag[1]*fr[17]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[17]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[13]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[13]*dfac_x*wm)/m_+(2.121320343559642*Bmag[2]*fr[10]*dfac_x*wm)/m_-(2.121320343559642*Bmag[2]*fl[10]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fr[5]*dfac_x*wm)/m_-(1.224744871391589*Bmag[2]*fl[5]*dfac_x*wm)/m_-(2.738612787525831*Phi[2]*fr[18]*dfac_x*q_)/m_-(2.738612787525831*Phi[2]*fl[18]*dfac_x*q_)/m_+(1.060660171779821*Phi[1]*fr[17]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[17]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[13]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[13]*dfac_x*q_)/m_+(2.121320343559642*Phi[2]*fr[10]*dfac_x*q_)/m_-(2.121320343559642*Phi[2]*fl[10]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fr[5]*dfac_x*q_)/m_-(1.224744871391589*Phi[2]*fl[5]*dfac_x*q_)/m_-(0.7071067811865475*dApardt[1]*fr[18]*q_)/m_-(0.7071067811865475*dApardt[1]*fl[18]*q_)/m_+(0.3912303982179757*dApardt[2]*fr[17]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[17]*q_)/m_-(0.3912303982179757*dApardt[2]*fl[17]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[17]*q_)/m_-(0.7905694150420947*dApardt[2]*fr[14]*q_)/m_-(0.7905694150420947*dApardt[2]*fl[14]*q_)/m_-(0.2258769757263128*dApardt[2]*fr[13]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[13]*q_)/m_-(0.2258769757263128*dApardt[2]*fl[13]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[13]*q_)/m_+(0.5477225575051661*dApardt[1]*fr[10]*q_)/m_-(0.5477225575051661*dApardt[1]*fl[10]*q_)/m_+(0.6123724356957944*dApardt[2]*fr[6]*q_)/m_-(0.6123724356957944*dApardt[2]*fl[6]*q_)/m_-(0.3162277660168379*dApardt[1]*fr[5]*q_)/m_-(0.3162277660168379*dApardt[1]*fl[5]*q_)/m_-(0.3535533905932737*dApardt[2]*fr[3]*q_)/m_-(0.3535533905932737*dApardt[2]*fl[3]*q_)/m_+(1.095445115010332*Bmag[2]*fr[19]*dfac_x)/(dfac_m*m_)-(1.095445115010332*Bmag[2]*fl[19]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fr[15]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fl[15]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fr[12]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[12]*dfac_x)/(dfac_m*m_)+(0.6123724356957944*Bmag[1]*fr[11]*dfac_x)/(dfac_m*m_)-(0.6123724356957944*Bmag[1]*fl[11]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fr[7]*dfac_x)/(dfac_m*m_)-(0.3535533905932737*Bmag[1]*fl[7]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[4]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[4]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*fr[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*fl[1]*Bmag[2]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax; 
  Ghat[15] = (1.060660171779821*Bmag[1]*fr[19]*dfac_x*wm)/m_-(1.060660171779821*Bmag[1]*fl[19]*dfac_x*wm)/m_+(2.371708245126284*Bmag[2]*fr[16]*dfac_x*wm)/m_-(2.371708245126284*Bmag[2]*fl[16]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fr[15]*dfac_x*wm)/m_-(0.6123724356957944*Bmag[1]*fl[15]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fr[9]*dfac_x*wm)/m_-(1.369306393762915*Bmag[2]*fl[9]*dfac_x*wm)/m_+(1.060660171779821*Phi[1]*fr[19]*dfac_x*q_)/m_-(1.060660171779821*Phi[1]*fl[19]*dfac_x*q_)/m_+(2.371708245126284*Phi[2]*fr[16]*dfac_x*q_)/m_-(2.371708245126284*Phi[2]*fl[16]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fr[15]*dfac_x*q_)/m_-(0.6123724356957944*Phi[1]*fl[15]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fr[9]*dfac_x*q_)/m_-(1.369306393762915*Phi[2]*fl[9]*dfac_x*q_)/m_+(0.5477225575051661*dApardt[2]*fr[19]*q_)/m_+(0.6123724356957944*dApardt[0]*fr[19]*q_)/m_-(0.5477225575051661*dApardt[2]*fl[19]*q_)/m_-(0.6123724356957944*dApardt[0]*fl[19]*q_)/m_+(0.6123724356957944*dApardt[1]*fr[16]*q_)/m_-(0.6123724356957944*dApardt[1]*fl[16]*q_)/m_-(0.3162277660168379*dApardt[2]*fr[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fr[15]*q_)/m_-(0.3162277660168379*dApardt[2]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[0]*fl[15]*q_)/m_-(0.3535533905932737*dApardt[1]*fr[9]*q_)/m_-(0.3535533905932737*dApardt[1]*fl[9]*q_)/m_-(0.7071067811865475*Bmag[1]*fr[18]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[1]*fl[18]*dfac_x)/(dfac_m*m_)+(1.095445115010332*Bmag[2]*fr[17]*dfac_x)/(dfac_m*m_)-(1.095445115010332*Bmag[2]*fl[17]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fr[14]*dfac_x)/(dfac_m*m_)-(1.58113883008419*Bmag[2]*fl[14]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fr[13]*dfac_x)/(dfac_m*m_)-(0.6324555320336759*Bmag[2]*fl[13]*dfac_x)/(dfac_m*m_)+(0.5477225575051661*Bmag[1]*fr[10]*dfac_x)/(dfac_m*m_)-(0.5477225575051661*Bmag[1]*fl[10]*dfac_x)/(dfac_m*m_)+(1.224744871391589*Bmag[2]*fr[6]*dfac_x)/(dfac_m*m_)-(1.224744871391589*Bmag[2]*fl[6]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fr[5]*dfac_x)/(dfac_m*m_)-(0.3162277660168379*Bmag[1]*fl[5]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fr[3]*dfac_x)/(dfac_m*m_)-(0.7071067811865475*Bmag[2]*fl[3]*dfac_x)/(dfac_m*m_)+0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax; 

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
