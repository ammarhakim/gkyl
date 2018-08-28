#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = Gradpar[0]*wv; 

  double alpha[8]; 
  alpha[0] = 2.0*Gradpar[0]*wv; 
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
double GyrokineticSurf1x2vSer_X_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = Gradpar[0]*wv; 

  double alpha[20]; 
  alpha[0] = 2.0*Gradpar[0]*wv; 
  alpha[2] = (1.154700538379252*Gradpar[0])/dfac_v; 
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
double GyrokineticSurf1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.224744871391589*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[8]; 
  alpha[0] = -(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
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
double GyrokineticSurf1x2vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.224744871391589*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[20]; 
  alpha[0] = -(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(5.477225575051662*Gradpar[0]*Phi[2]*dfac_x*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-1.118033988749895*fr[8]*amax)+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax+0.3952847075210473*alpha[1]*fr[12]+0.3952847075210473*alpha[1]*fl[12]+0.3952847075210473*alpha[0]*fr[8]+0.3952847075210473*alpha[0]*fl[8]-0.3061862178478971*alpha[1]*fr[4]+0.3061862178478971*alpha[1]*fl[4]-0.3061862178478971*alpha[0]*fr[2]+0.3061862178478971*alpha[0]*fl[2]+0.1767766952966368*alpha[1]*fr[1]+0.1767766952966368*alpha[1]*fl[1]+0.1767766952966368*alpha[0]*fr[0]+0.1767766952966368*alpha[0]*fl[0]; 
  Ghat[1] = (-1.118033988749895*fr[12]*amax)+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax+0.3952847075210473*alpha[0]*fr[12]+0.3952847075210473*alpha[0]*fl[12]-0.273861278752583*alpha[1]*fr[11]+0.273861278752583*alpha[1]*fl[11]+0.3952847075210473*alpha[1]*fr[8]+0.3952847075210473*alpha[1]*fl[8]+0.1581138830084189*alpha[1]*fr[7]+0.1581138830084189*alpha[1]*fl[7]-0.3061862178478971*alpha[0]*fr[4]+0.3061862178478971*alpha[0]*fl[4]-0.3061862178478971*alpha[1]*fr[2]+0.3061862178478971*alpha[1]*fl[2]+0.1767766952966368*alpha[0]*fr[1]+0.1767766952966368*alpha[0]*fl[1]+0.1767766952966368*fr[0]*alpha[1]+0.1767766952966368*fl[0]*alpha[1]; 
  Ghat[3] = (-1.118033988749895*fr[14]*amax)+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax+0.3952847075210473*alpha[1]*fr[18]+0.3952847075210473*alpha[1]*fl[18]+0.3952847075210473*alpha[0]*fr[14]+0.3952847075210473*alpha[0]*fl[14]-0.3061862178478971*alpha[1]*fr[10]+0.3061862178478971*alpha[1]*fl[10]-0.3061862178478971*alpha[0]*fr[6]+0.3061862178478971*alpha[0]*fl[6]+0.1767766952966368*alpha[1]*fr[5]+0.1767766952966368*alpha[1]*fl[5]+0.1767766952966368*alpha[0]*fr[3]+0.1767766952966368*alpha[0]*fl[3]; 
  Ghat[5] = (-1.118033988749895*fr[18]*amax)+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax+0.3952847075210473*alpha[0]*fr[18]+0.3952847075210473*alpha[0]*fl[18]-0.273861278752583*alpha[1]*fr[17]+0.273861278752583*alpha[1]*fl[17]+0.3952847075210473*alpha[1]*fr[14]+0.3952847075210473*alpha[1]*fl[14]+0.1581138830084189*alpha[1]*fr[13]+0.1581138830084189*alpha[1]*fl[13]-0.3061862178478971*alpha[0]*fr[10]+0.3061862178478971*alpha[0]*fl[10]-0.3061862178478971*alpha[1]*fr[6]+0.3061862178478971*alpha[1]*fl[6]+0.1767766952966368*alpha[0]*fr[5]+0.1767766952966368*alpha[0]*fl[5]+0.1767766952966368*alpha[1]*fr[3]+0.1767766952966368*alpha[1]*fl[3]; 
  Ghat[7] = 0.8660254037844386*fr[11]*amax+0.8660254037844386*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax+0.3535533905932737*alpha[1]*fr[12]+0.3535533905932737*alpha[1]*fl[12]-0.3061862178478972*alpha[0]*fr[11]+0.3061862178478972*alpha[0]*fl[11]+0.1767766952966368*alpha[0]*fr[7]+0.1767766952966368*alpha[0]*fl[7]-0.273861278752583*alpha[1]*fr[4]+0.273861278752583*alpha[1]*fl[4]+0.1581138830084189*alpha[1]*fr[1]+0.1581138830084189*alpha[1]*fl[1]; 
  Ghat[9] = 0.8660254037844386*fr[16]*amax+0.8660254037844386*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax-0.3061862178478971*alpha[1]*fr[19]+0.3061862178478971*alpha[1]*fl[19]-0.3061862178478972*alpha[0]*fr[16]+0.3061862178478972*alpha[0]*fl[16]+0.1767766952966368*alpha[1]*fr[15]+0.1767766952966368*alpha[1]*fl[15]+0.1767766952966368*alpha[0]*fr[9]+0.1767766952966368*alpha[0]*fl[9]; 
  Ghat[13] = 0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax+0.3535533905932737*alpha[1]*fr[18]+0.3535533905932737*alpha[1]*fl[18]-0.3061862178478971*alpha[0]*fr[17]+0.3061862178478971*alpha[0]*fl[17]+0.1767766952966368*alpha[0]*fr[13]+0.1767766952966368*alpha[0]*fl[13]-0.2738612787525829*alpha[1]*fr[10]+0.2738612787525829*alpha[1]*fl[10]+0.1581138830084189*alpha[1]*fr[5]+0.1581138830084189*alpha[1]*fl[5]; 
  Ghat[15] = 0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax-0.3061862178478971*alpha[0]*fr[19]+0.3061862178478971*alpha[0]*fl[19]-0.3061862178478971*alpha[1]*fr[16]+0.3061862178478971*alpha[1]*fl[16]+0.1767766952966368*alpha[0]*fr[15]+0.1767766952966368*alpha[0]*fl[15]+0.1767766952966368*alpha[1]*fr[9]+0.1767766952966368*alpha[1]*fl[9]; 

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
double EmGyrokineticSurf1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = 0.7071067811865475*(1.732050807568877*geoY[0]*Apar[1]*dfac_x+1.414213562373095*Gradpar[0])*wv; 

  double alpha[8]; 
  alpha[0] = 2.449489742783178*geoY[0]*Apar[1]*dfac_x*wv+2.0*Gradpar[0]*wv; 
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
double EmGyrokineticSurf1x2vSer_X_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -0.5*((9.48683298050514*geoY[0]*Apar[2]-2.449489742783178*geoY[0]*Apar[1])*dfac_x-2.0*Gradpar[0])*wv; 

  double alpha[20]; 
  alpha[0] = (-9.486832980505138*geoY[0]*Apar[2]*dfac_x*wv)+2.449489742783178*geoY[0]*Apar[1]*dfac_x*wv+2.0*Gradpar[0]*wv; 
  alpha[2] = (-(5.477225575051661*geoY[0]*Apar[2]*dfac_x)/dfac_v)+(1.414213562373095*geoY[0]*Apar[1]*dfac_x)/dfac_v+(1.154700538379251*Gradpar[0])/dfac_v; 
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
double EmGyrokineticSurf1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.3535533905932737*(4.242640687119286*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x+3.464101615137754*Gradpar[0]*Phi[1]*dfac_x+2.828427124746191*dApardt[0])*q_)/m_; 

  double alpha[8]; 
  alpha[0] = (-(3.0*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x*q_)/m_)-(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[0]*q_)/m_; 
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
double EmGyrokineticSurf1x2vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.3535533905932737*((21.21320343559643*geoY[0]*Apar[2]*Phi[2]+4.242640687119286*geoY[0]*Apar[1]*Phi[1])*dfac_x*dfac_x+3.464101615137754*Gradpar[0]*Phi[1]*dfac_x+2.828427124746191*dApardt[0])*q_)/m_; 

  double alpha[20]; 
  alpha[0] = (-(15.0*geoY[0]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_)-(3.0*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[0]*q_)/m_; 
  alpha[1] = (-(6.708203932499369*geoY[0]*Apar[1]*Phi[2]*dfac_x*dfac_x*q_)/m_)-(6.708203932499369*geoY[0]*Phi[1]*Apar[2]*dfac_x*dfac_x*q_)/m_-(5.477225575051662*Gradpar[0]*Phi[2]*dfac_x*q_)/m_-(2.0*dApardt[1]*q_)/m_; 
  alpha[7] = (-(13.41640786499874*geoY[0]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_)-(2.0*dApardt[2]*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-1.118033988749895*fr[8]*amax)+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax+0.3952847075210473*alpha[1]*fr[12]+0.3952847075210473*alpha[1]*fl[12]-0.3061862178478971*alpha[7]*fr[11]+0.3061862178478971*alpha[7]*fl[11]+0.3952847075210473*alpha[0]*fr[8]+0.3952847075210473*alpha[0]*fl[8]+0.1767766952966368*alpha[7]*fr[7]+0.1767766952966368*alpha[7]*fl[7]-0.3061862178478971*alpha[1]*fr[4]+0.3061862178478971*alpha[1]*fl[4]-0.3061862178478971*alpha[0]*fr[2]+0.3061862178478971*alpha[0]*fl[2]+0.1767766952966368*alpha[1]*fr[1]+0.1767766952966368*alpha[1]*fl[1]+0.1767766952966368*alpha[0]*fr[0]+0.1767766952966368*alpha[0]*fl[0]; 
  Ghat[1] = (-1.118033988749895*fr[12]*amax)+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax+0.3535533905932737*alpha[7]*fr[12]+0.3952847075210473*alpha[0]*fr[12]+0.3535533905932737*alpha[7]*fl[12]+0.3952847075210473*alpha[0]*fl[12]-0.273861278752583*alpha[1]*fr[11]+0.273861278752583*alpha[1]*fl[11]+0.3952847075210473*alpha[1]*fr[8]+0.3952847075210473*alpha[1]*fl[8]+0.1581138830084189*alpha[1]*fr[7]+0.1581138830084189*alpha[1]*fl[7]-0.273861278752583*fr[4]*alpha[7]+0.273861278752583*fl[4]*alpha[7]+0.1581138830084189*fr[1]*alpha[7]+0.1581138830084189*fl[1]*alpha[7]-0.3061862178478971*alpha[0]*fr[4]+0.3061862178478971*alpha[0]*fl[4]-0.3061862178478971*alpha[1]*fr[2]+0.3061862178478971*alpha[1]*fl[2]+0.1767766952966368*alpha[0]*fr[1]+0.1767766952966368*alpha[0]*fl[1]+0.1767766952966368*fr[0]*alpha[1]+0.1767766952966368*fl[0]*alpha[1]; 
  Ghat[3] = (-1.118033988749895*fr[14]*amax)+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax+0.3952847075210473*alpha[1]*fr[18]+0.3952847075210473*alpha[1]*fl[18]-0.3061862178478971*alpha[7]*fr[17]+0.3061862178478971*alpha[7]*fl[17]+0.3952847075210473*alpha[0]*fr[14]+0.3952847075210473*alpha[0]*fl[14]+0.1767766952966368*alpha[7]*fr[13]+0.1767766952966368*alpha[7]*fl[13]-0.3061862178478971*alpha[1]*fr[10]+0.3061862178478971*alpha[1]*fl[10]-0.3061862178478971*alpha[0]*fr[6]+0.3061862178478971*alpha[0]*fl[6]+0.1767766952966368*alpha[1]*fr[5]+0.1767766952966368*alpha[1]*fl[5]+0.1767766952966368*alpha[0]*fr[3]+0.1767766952966368*alpha[0]*fl[3]; 
  Ghat[5] = (-1.118033988749895*fr[18]*amax)+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax+0.3535533905932737*alpha[7]*fr[18]+0.3952847075210473*alpha[0]*fr[18]+0.3535533905932737*alpha[7]*fl[18]+0.3952847075210473*alpha[0]*fl[18]-0.273861278752583*alpha[1]*fr[17]+0.273861278752583*alpha[1]*fl[17]+0.3952847075210473*alpha[1]*fr[14]+0.3952847075210473*alpha[1]*fl[14]+0.1581138830084189*alpha[1]*fr[13]+0.1581138830084189*alpha[1]*fl[13]-0.273861278752583*alpha[7]*fr[10]-0.3061862178478971*alpha[0]*fr[10]+0.273861278752583*alpha[7]*fl[10]+0.3061862178478971*alpha[0]*fl[10]+0.1581138830084189*fr[5]*alpha[7]+0.1581138830084189*fl[5]*alpha[7]-0.3061862178478971*alpha[1]*fr[6]+0.3061862178478971*alpha[1]*fl[6]+0.1767766952966368*alpha[0]*fr[5]+0.1767766952966368*alpha[0]*fl[5]+0.1767766952966368*alpha[1]*fr[3]+0.1767766952966368*alpha[1]*fl[3]; 
  Ghat[7] = 0.8660254037844386*fr[11]*amax+0.8660254037844386*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax+0.3535533905932737*alpha[1]*fr[12]+0.3535533905932737*alpha[1]*fl[12]-0.1956151991089878*alpha[7]*fr[11]-0.3061862178478972*alpha[0]*fr[11]+0.1956151991089878*alpha[7]*fl[11]+0.3061862178478972*alpha[0]*fl[11]+0.3952847075210473*alpha[7]*fr[8]+0.3952847075210473*alpha[7]*fl[8]+0.1129384878631564*alpha[7]*fr[7]+0.1767766952966368*alpha[0]*fr[7]+0.1129384878631564*alpha[7]*fl[7]+0.1767766952966368*alpha[0]*fl[7]-0.3061862178478971*fr[2]*alpha[7]+0.3061862178478971*fl[2]*alpha[7]+0.1767766952966368*fr[0]*alpha[7]+0.1767766952966368*fl[0]*alpha[7]-0.273861278752583*alpha[1]*fr[4]+0.273861278752583*alpha[1]*fl[4]+0.1581138830084189*alpha[1]*fr[1]+0.1581138830084189*alpha[1]*fl[1]; 
  Ghat[9] = 0.8660254037844386*fr[16]*amax+0.8660254037844386*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax-0.3061862178478971*alpha[1]*fr[19]+0.3061862178478971*alpha[1]*fl[19]-0.3061862178478972*alpha[0]*fr[16]+0.3061862178478972*alpha[0]*fl[16]+0.1767766952966368*alpha[1]*fr[15]+0.1767766952966368*alpha[1]*fl[15]+0.1767766952966368*alpha[0]*fr[9]+0.1767766952966368*alpha[0]*fl[9]; 
  Ghat[13] = 0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax+0.3535533905932737*alpha[1]*fr[18]+0.3535533905932737*alpha[1]*fl[18]-0.1956151991089878*alpha[7]*fr[17]-0.3061862178478971*alpha[0]*fr[17]+0.1956151991089878*alpha[7]*fl[17]+0.3061862178478971*alpha[0]*fl[17]+0.3952847075210473*alpha[7]*fr[14]+0.3952847075210473*alpha[7]*fl[14]+0.1129384878631564*alpha[7]*fr[13]+0.1767766952966368*alpha[0]*fr[13]+0.1129384878631564*alpha[7]*fl[13]+0.1767766952966368*alpha[0]*fl[13]-0.2738612787525829*alpha[1]*fr[10]+0.2738612787525829*alpha[1]*fl[10]-0.3061862178478971*fr[6]*alpha[7]+0.3061862178478971*fl[6]*alpha[7]+0.1767766952966368*fr[3]*alpha[7]+0.1767766952966368*fl[3]*alpha[7]+0.1581138830084189*alpha[1]*fr[5]+0.1581138830084189*alpha[1]*fl[5]; 
  Ghat[15] = 0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax-0.2738612787525829*alpha[7]*fr[19]-0.3061862178478971*alpha[0]*fr[19]+0.2738612787525829*alpha[7]*fl[19]+0.3061862178478971*alpha[0]*fl[19]-0.3061862178478971*alpha[1]*fr[16]+0.3061862178478971*alpha[1]*fl[16]+0.1581138830084189*alpha[7]*fr[15]+0.1767766952966368*alpha[0]*fr[15]+0.1581138830084189*alpha[7]*fl[15]+0.1767766952966368*alpha[0]*fl[15]+0.1767766952966368*alpha[1]*fr[9]+0.1767766952966368*alpha[1]*fl[9]; 

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
double GyrokineticSurf1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.5*(((5.196152422706631*Bmag[1]*BmagInv[1]-3.0*BmagInv[0]*Bmag[1])*geoY[1]-3.0*geoY[0]*Bmag[1]*BmagInv[1]+1.732050807568877*BmagInv[0]*geoY[0]*Bmag[1])*dfac_x*m_*wv2+(3.464101615137754*Gradpar[1]-2.0*Gradpar[0])*q_*wv))/q_; 

  double alpha[8]; 
  alpha[0] = (-(5.196152422706631*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv2)/q_)+(3.0*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv2)/q_+(3.0*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv2)/q_-(1.732050807568877*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv2)/q_-3.464101615137754*Gradpar[1]*wv+2.0*Gradpar[0]*wv; 
  alpha[2] = (-(3.0*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_))+(1.732050807568877*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(1.732050807568877*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv)/(dfac_v*q_)-(1.0*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[2]*fl[4]+alpha[2]*fl[2]+1.732050807568877*alpha[0]*fl[1]+alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*(3.0*alpha[2]*fl[4]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[4]+alpha[0]*fl[2]+1.732050807568877*fl[1]*alpha[2]+fl[0]*alpha[2])*dfac_x; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[2]*fl[7]+alpha[2]*fl[6]+1.732050807568877*alpha[0]*fl[5]+alpha[0]*fl[3])*dfac_x; 
  incr[4] = -0.1767766952966368*(3.0*alpha[0]*fl[4]+1.732050807568877*alpha[0]*fl[2]+3.0*fl[1]*alpha[2]+1.732050807568877*fl[0]*alpha[2])*dfac_x; 
  incr[5] = -0.1767766952966368*(3.0*alpha[2]*fl[7]+1.732050807568877*alpha[2]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[3])*dfac_x; 
  incr[6] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[7]+alpha[0]*fl[6]+1.732050807568877*alpha[2]*fl[5]+alpha[2]*fl[3])*dfac_x; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fl[7]+1.732050807568877*alpha[0]*fl[6]+3.0*alpha[2]*fl[5]+1.732050807568877*alpha[2]*fl[3])*dfac_x; 

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
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[2]*fr[4]-1.0*alpha[2]*fr[2]+1.732050807568877*alpha[0]*fr[1]-1.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.1767766952966368*(3.0*alpha[2]*fr[4]-1.732050807568877*alpha[2]*fr[2]+3.0*alpha[0]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[4]-1.0*alpha[0]*fr[2]+1.732050807568877*fr[1]*alpha[2]-1.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[2]*fr[7]-1.0*alpha[2]*fr[6]+1.732050807568877*alpha[0]*fr[5]-1.0*alpha[0]*fr[3])*dfac_x; 
  incr[4] = 0.1767766952966368*(3.0*alpha[0]*fr[4]-1.732050807568877*alpha[0]*fr[2]+3.0*fr[1]*alpha[2]-1.732050807568877*fr[0]*alpha[2])*dfac_x; 
  incr[5] = 0.1767766952966368*(3.0*alpha[2]*fr[7]-1.732050807568877*alpha[2]*fr[6]+3.0*alpha[0]*fr[5]-1.732050807568877*alpha[0]*fr[3])*dfac_x; 
  incr[6] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[7]-1.0*alpha[0]*fr[6]+1.732050807568877*alpha[2]*fr[5]-1.0*alpha[2]*fr[3])*dfac_x; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fr[7]-1.732050807568877*alpha[0]*fr[6]+3.0*alpha[2]*fr[5]-1.732050807568877*alpha[2]*fr[3])*dfac_x; 

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
double GyrokineticSurf1x2vSer_X_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.1666666666666667*((((100.6230589874906*Bmag[2]-25.98076211353316*Bmag[1])*BmagInv[2]+(45.0*BmagInv[0]-77.94228634059945*BmagInv[1])*Bmag[2]+20.12461179749811*Bmag[1]*BmagInv[1]-11.61895003862225*BmagInv[0]*Bmag[1])*geoY[2]+((45.0*geoY[0]-77.94228634059945*geoY[1])*Bmag[2]+20.12461179749811*Bmag[1]*geoY[1]-11.61895003862225*geoY[0]*Bmag[1])*BmagInv[2]+((60.37383539249433*BmagInv[1]-34.85685011586674*BmagInv[0])*geoY[1]-34.85685011586674*geoY[0]*BmagInv[1]+20.12461179749811*BmagInv[0]*geoY[0])*Bmag[2]+(9.0*BmagInv[0]*Bmag[1]-15.58845726811989*Bmag[1]*BmagInv[1])*geoY[1]+9.0*geoY[0]*Bmag[1]*BmagInv[1]-5.196152422706631*BmagInv[0]*geoY[0]*Bmag[1])*dfac_v2*dfac_x*m_*wv2+(13.41640786499874*Gradpar[2]-10.39230484541326*Gradpar[1]+6.0*Gradpar[0])*dfac_v2*q_*wv+(((33.54101966249685*Bmag[2]-8.660254037844386*Bmag[1])*BmagInv[2]+(15.0*BmagInv[0]-25.98076211353316*BmagInv[1])*Bmag[2]+6.708203932499369*Bmag[1]*BmagInv[1]-3.872983346207417*BmagInv[0]*Bmag[1])*geoY[2]+((15.0*geoY[0]-25.98076211353316*geoY[1])*Bmag[2]+6.708203932499369*Bmag[1]*geoY[1]-3.872983346207417*geoY[0]*Bmag[1])*BmagInv[2]+((20.12461179749811*BmagInv[1]-11.61895003862225*BmagInv[0])*geoY[1]-11.61895003862225*geoY[0]*BmagInv[1]+6.708203932499369*BmagInv[0]*geoY[0])*Bmag[2]+(3.0*BmagInv[0]*Bmag[1]-5.196152422706631*Bmag[1]*BmagInv[1])*geoY[1]+3.0*geoY[0]*Bmag[1]*BmagInv[1]-1.732050807568877*BmagInv[0]*geoY[0]*Bmag[1])*dfac_x*m_))/(dfac_v2*q_); 

  double alpha[20]; 
  alpha[0] = (33.54101966249685*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_*wv2)/q_-(8.660254037844386*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_*wv2)/q_-(25.98076211353315*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_*wv2)/q_+(15.0*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_*wv2)/q_+(6.708203932499369*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_*wv2)/q_-(3.872983346207417*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_*wv2)/q_-(25.98076211353315*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv2)/q_+(15.0*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv2)/q_+(6.708203932499369*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_*wv2)/q_-(3.872983346207417*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_*wv2)/q_+(20.12461179749811*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_*wv2)/q_-(11.61895003862225*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_*wv2)/q_-(11.61895003862225*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_*wv2)/q_+(6.708203932499369*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_*wv2)/q_-(5.19615242270663*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv2)/q_+(3.0*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv2)/q_+(3.0*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv2)/q_-(1.732050807568877*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv2)/q_+4.47213595499958*Gradpar[2]*wv-3.464101615137754*Gradpar[1]*wv+2.0*Gradpar[0]*wv+(11.18033988749895*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(2.886751345948128*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(8.660254037844386*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(5.0*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(2.23606797749979*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(1.290994448735806*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(8.660254037844386*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(5.0*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(2.23606797749979*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)-(1.290994448735806*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(6.708203932499369*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(3.872983346207417*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(3.872983346207417*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)+(2.23606797749979*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(1.732050807568877*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(1.0*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(1.0*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_)/(dfac_v2*q_)-(0.5773502691896257*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_)/(dfac_v2*q_); 
  alpha[2] = (38.72983346207418*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(10.0*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(30.0*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)+(17.32050807568877*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)+(7.745966692414834*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(4.47213595499958*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(30.0*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(17.32050807568877*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(7.745966692414834*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)-(4.47213595499958*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(23.2379000772445*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(13.41640786499874*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(13.41640786499874*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)+(7.745966692414834*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(6.0*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(3.464101615137754*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(3.464101615137754*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv)/(dfac_v*q_)-(2.0*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv)/(dfac_v*q_)+(2.581988897471611*Gradpar[2])/dfac_v-(2.0*Gradpar[1])/dfac_v+(1.154700538379251*Gradpar[0])/dfac_v; 
  alpha[8] = (10.0*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(2.581988897471612*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(7.745966692414835*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(4.47213595499958*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(2.0*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(1.154700538379251*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(7.745966692414835*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(4.47213595499958*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(2.0*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)-(1.154700538379251*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(6.0*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(3.464101615137754*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(3.464101615137754*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)+(2.0*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(1.549193338482967*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.8944271909999159*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.8944271909999159*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_)/(dfac_v2*q_)-(0.5163977794943223*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_)/(dfac_v2*q_); 
  if (alpha0>0) { 
  incr[0] = 0.01178511301977579*(25.98076211353316*alpha[8]*fl[12]+33.54101966249684*alpha[2]*fl[11]+15.0*alpha[8]*fl[8]+33.54101966249685*alpha[0]*fl[7]+25.98076211353316*alpha[2]*fl[4]+15.0*alpha[2]*fl[2]+25.98076211353316*alpha[0]*fl[1]+15.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.03535533905932736*(15.0*alpha[8]*fl[12]+19.36491673103708*alpha[2]*fl[11]+8.660254037844386*alpha[8]*fl[8]+19.36491673103709*alpha[0]*fl[7]+15.0*alpha[2]*fl[4]+8.660254037844386*alpha[2]*fl[2]+15.0*alpha[0]*fl[1]+8.660254037844386*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[12]+30.0*alpha[8]*fl[11]+33.54101966249684*alpha[0]*fl[11]+13.41640786499874*alpha[2]*fl[8]+23.2379000772445*fl[4]*alpha[8]+13.41640786499874*fl[2]*alpha[8]+33.54101966249685*alpha[2]*fl[7]+25.98076211353316*alpha[0]*fl[4]+15.0*alpha[0]*fl[2]+25.98076211353316*fl[1]*alpha[2]+15.0*fl[0]*alpha[2])*dfac_x; 
  incr[3] = 0.01178511301977579*(25.98076211353316*alpha[8]*fl[18]+33.54101966249685*alpha[2]*fl[17]+15.0*alpha[8]*fl[14]+33.54101966249684*alpha[0]*fl[13]+25.98076211353316*alpha[2]*fl[10]+15.0*alpha[2]*fl[6]+25.98076211353316*alpha[0]*fl[5]+15.0*alpha[0]*fl[3])*dfac_x; 
  incr[4] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[12]+17.32050807568877*alpha[8]*fl[11]+19.36491673103708*alpha[0]*fl[11]+7.745966692414834*alpha[2]*fl[8]+13.41640786499874*fl[4]*alpha[8]+7.745966692414834*fl[2]*alpha[8]+19.36491673103709*alpha[2]*fl[7]+15.0*alpha[0]*fl[4]+8.660254037844386*alpha[0]*fl[2]+15.0*fl[1]*alpha[2]+8.660254037844386*fl[0]*alpha[2])*dfac_x; 
  incr[5] = -0.03535533905932736*(15.0*alpha[8]*fl[18]+19.36491673103709*alpha[2]*fl[17]+8.660254037844387*alpha[8]*fl[14]+19.36491673103708*alpha[0]*fl[13]+15.0*alpha[2]*fl[10]+8.660254037844386*alpha[2]*fl[6]+15.0*alpha[0]*fl[5]+8.660254037844386*alpha[0]*fl[3])*dfac_x; 
  incr[6] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[18]+30.0*alpha[8]*fl[17]+33.54101966249685*alpha[0]*fl[17]+13.41640786499874*alpha[2]*fl[14]+33.54101966249684*alpha[2]*fl[13]+23.2379000772445*alpha[8]*fl[10]+25.98076211353316*alpha[0]*fl[10]+13.41640786499874*fl[6]*alpha[8]+15.0*alpha[0]*fl[6]+25.98076211353316*alpha[2]*fl[5]+15.0*alpha[2]*fl[3])*dfac_x; 
  incr[7] = 0.05892556509887893*(11.61895003862225*alpha[8]*fl[12]+15.0*alpha[2]*fl[11]+6.708203932499369*alpha[8]*fl[8]+15.0*alpha[0]*fl[7]+11.61895003862225*alpha[2]*fl[4]+6.708203932499369*alpha[2]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[8] = 0.001683587574253684*(116.1895003862225*alpha[8]*fl[12]+181.8653347947321*alpha[0]*fl[12]+210.0*alpha[2]*fl[11]+67.0820393249937*alpha[8]*fl[8]+105.0*alpha[0]*fl[8]+234.787137637478*fl[7]*alpha[8]+181.8653347947321*fl[1]*alpha[8]+105.0*fl[0]*alpha[8]+162.6653005407115*alpha[2]*fl[4]+93.91485505499116*alpha[2]*fl[2])*dfac_x; 
  incr[9] = 0.01178511301977579*(25.98076211353316*alpha[2]*fl[19]+15.0*alpha[2]*fl[16]+25.98076211353316*alpha[0]*fl[15]+15.0*alpha[0]*fl[9])*dfac_x; 
  incr[10] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[18]+17.32050807568877*alpha[8]*fl[17]+19.36491673103709*alpha[0]*fl[17]+7.745966692414834*alpha[2]*fl[14]+19.36491673103708*alpha[2]*fl[13]+13.41640786499874*alpha[8]*fl[10]+15.0*alpha[0]*fl[10]+7.745966692414834*fl[6]*alpha[8]+8.660254037844386*alpha[0]*fl[6]+15.0*alpha[2]*fl[5]+8.660254037844386*alpha[2]*fl[3])*dfac_x; 
  incr[11] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[12]+67.0820393249937*alpha[8]*fl[11]+75.0*alpha[0]*fl[11]+30.0*alpha[2]*fl[8]+51.96152422706632*fl[4]*alpha[8]+30.0*fl[2]*alpha[8]+75.00000000000001*alpha[2]*fl[7]+58.09475019311126*alpha[0]*fl[4]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*fl[1]*alpha[2]+33.54101966249684*fl[0]*alpha[2])*dfac_x; 
  incr[12] = -0.005050762722761052*(67.0820393249937*alpha[8]*fl[12]+105.0*alpha[0]*fl[12]+121.2435565298214*alpha[2]*fl[11]+38.72983346207417*alpha[8]*fl[8]+60.62177826491071*alpha[0]*fl[8]+135.5544171172596*fl[7]*alpha[8]+105.0*fl[1]*alpha[8]+60.62177826491071*fl[0]*alpha[8]+93.91485505499116*alpha[2]*fl[4]+54.22176684690384*alpha[2]*fl[2])*dfac_x; 
  incr[13] = 0.05892556509887893*(11.61895003862225*alpha[8]*fl[18]+15.0*alpha[2]*fl[17]+6.708203932499369*alpha[8]*fl[14]+15.0*alpha[0]*fl[13]+11.61895003862225*alpha[2]*fl[10]+6.708203932499369*alpha[2]*fl[6]+11.61895003862225*alpha[0]*fl[5]+6.708203932499369*alpha[0]*fl[3])*dfac_x; 
  incr[14] = 0.001683587574253684*(116.1895003862225*alpha[8]*fl[18]+181.8653347947321*alpha[0]*fl[18]+210.0*alpha[2]*fl[17]+67.0820393249937*alpha[8]*fl[14]+105.0*alpha[0]*fl[14]+234.787137637478*alpha[8]*fl[13]+162.6653005407115*alpha[2]*fl[10]+181.8653347947321*fl[5]*alpha[8]+105.0*fl[3]*alpha[8]+93.91485505499116*alpha[2]*fl[6])*dfac_x; 
  incr[15] = -0.03535533905932736*(15.0*alpha[2]*fl[19]+8.660254037844386*alpha[2]*fl[16]+15.0*alpha[0]*fl[15]+8.660254037844387*alpha[0]*fl[9])*dfac_x; 
  incr[16] = 0.01178511301977579*(23.2379000772445*alpha[8]*fl[19]+25.98076211353316*alpha[0]*fl[19]+13.41640786499874*alpha[8]*fl[16]+15.0*alpha[0]*fl[16]+25.98076211353316*alpha[2]*fl[15]+15.0*alpha[2]*fl[9])*dfac_x; 
  incr[17] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[18]+67.0820393249937*alpha[8]*fl[17]+75.0*alpha[0]*fl[17]+30.0*alpha[2]*fl[14]+75.00000000000001*alpha[2]*fl[13]+51.96152422706631*alpha[8]*fl[10]+58.09475019311126*alpha[0]*fl[10]+30.0*fl[6]*alpha[8]+33.54101966249685*alpha[0]*fl[6]+58.09475019311126*alpha[2]*fl[5]+33.54101966249685*alpha[2]*fl[3])*dfac_x; 
  incr[18] = -0.005050762722761052*(67.0820393249937*alpha[8]*fl[18]+105.0*alpha[0]*fl[18]+121.2435565298214*alpha[2]*fl[17]+38.72983346207417*alpha[8]*fl[14]+60.62177826491071*alpha[0]*fl[14]+135.5544171172596*alpha[8]*fl[13]+93.91485505499116*alpha[2]*fl[10]+105.0*fl[5]*alpha[8]+60.6217782649107*fl[3]*alpha[8]+54.22176684690384*alpha[2]*fl[6])*dfac_x; 
  incr[19] = -0.03535533905932736*(13.41640786499874*alpha[8]*fl[19]+15.0*alpha[0]*fl[19]+7.745966692414834*alpha[8]*fl[16]+8.660254037844387*alpha[0]*fl[16]+15.0*alpha[2]*fl[15]+8.660254037844386*alpha[2]*fl[9])*dfac_x; 

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
  incr[0] = -0.01178511301977579*(25.98076211353316*alpha[8]*fr[12]-33.54101966249684*alpha[2]*fr[11]-15.0*alpha[8]*fr[8]-33.54101966249685*alpha[0]*fr[7]+25.98076211353316*alpha[2]*fr[4]-15.0*alpha[2]*fr[2]+25.98076211353316*alpha[0]*fr[1]-15.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.03535533905932736*(15.0*alpha[8]*fr[12]-19.36491673103708*alpha[2]*fr[11]-8.660254037844386*alpha[8]*fr[8]-19.36491673103709*alpha[0]*fr[7]+15.0*alpha[2]*fr[4]-8.660254037844386*alpha[2]*fr[2]+15.0*alpha[0]*fr[1]-8.660254037844386*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[12]-30.0*alpha[8]*fr[11]-33.54101966249684*alpha[0]*fr[11]-13.41640786499874*alpha[2]*fr[8]+23.2379000772445*fr[4]*alpha[8]-13.41640786499874*fr[2]*alpha[8]-33.54101966249685*alpha[2]*fr[7]+25.98076211353316*alpha[0]*fr[4]-15.0*alpha[0]*fr[2]+25.98076211353316*fr[1]*alpha[2]-15.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = -0.01178511301977579*(25.98076211353316*alpha[8]*fr[18]-33.54101966249685*alpha[2]*fr[17]-15.0*alpha[8]*fr[14]-33.54101966249684*alpha[0]*fr[13]+25.98076211353316*alpha[2]*fr[10]-15.0*alpha[2]*fr[6]+25.98076211353316*alpha[0]*fr[5]-15.0*alpha[0]*fr[3])*dfac_x; 
  incr[4] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[12]-17.32050807568877*alpha[8]*fr[11]-19.36491673103708*alpha[0]*fr[11]-7.745966692414834*alpha[2]*fr[8]+13.41640786499874*fr[4]*alpha[8]-7.745966692414834*fr[2]*alpha[8]-19.36491673103709*alpha[2]*fr[7]+15.0*alpha[0]*fr[4]-8.660254037844386*alpha[0]*fr[2]+15.0*fr[1]*alpha[2]-8.660254037844386*fr[0]*alpha[2])*dfac_x; 
  incr[5] = 0.03535533905932736*(15.0*alpha[8]*fr[18]-19.36491673103709*alpha[2]*fr[17]-8.660254037844387*alpha[8]*fr[14]-19.36491673103708*alpha[0]*fr[13]+15.0*alpha[2]*fr[10]-8.660254037844386*alpha[2]*fr[6]+15.0*alpha[0]*fr[5]-8.660254037844386*alpha[0]*fr[3])*dfac_x; 
  incr[6] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[18]-30.0*alpha[8]*fr[17]-33.54101966249685*alpha[0]*fr[17]-13.41640786499874*alpha[2]*fr[14]-33.54101966249684*alpha[2]*fr[13]+23.2379000772445*alpha[8]*fr[10]+25.98076211353316*alpha[0]*fr[10]-13.41640786499874*fr[6]*alpha[8]-15.0*alpha[0]*fr[6]+25.98076211353316*alpha[2]*fr[5]-15.0*alpha[2]*fr[3])*dfac_x; 
  incr[7] = -0.05892556509887893*(11.61895003862225*alpha[8]*fr[12]-15.0*alpha[2]*fr[11]-6.708203932499369*alpha[8]*fr[8]-15.0*alpha[0]*fr[7]+11.61895003862225*alpha[2]*fr[4]-6.708203932499369*alpha[2]*fr[2]+11.61895003862225*alpha[0]*fr[1]-6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[8] = -0.001683587574253684*(116.1895003862225*alpha[8]*fr[12]+181.8653347947321*alpha[0]*fr[12]-210.0*alpha[2]*fr[11]-67.0820393249937*alpha[8]*fr[8]-105.0*alpha[0]*fr[8]-234.787137637478*fr[7]*alpha[8]+181.8653347947321*fr[1]*alpha[8]-105.0*fr[0]*alpha[8]+162.6653005407115*alpha[2]*fr[4]-93.91485505499116*alpha[2]*fr[2])*dfac_x; 
  incr[9] = -0.01178511301977579*(25.98076211353316*alpha[2]*fr[19]-15.0*alpha[2]*fr[16]+25.98076211353316*alpha[0]*fr[15]-15.0*alpha[0]*fr[9])*dfac_x; 
  incr[10] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[18]-17.32050807568877*alpha[8]*fr[17]-19.36491673103709*alpha[0]*fr[17]-7.745966692414834*alpha[2]*fr[14]-19.36491673103708*alpha[2]*fr[13]+13.41640786499874*alpha[8]*fr[10]+15.0*alpha[0]*fr[10]-7.745966692414834*fr[6]*alpha[8]-8.660254037844386*alpha[0]*fr[6]+15.0*alpha[2]*fr[5]-8.660254037844386*alpha[2]*fr[3])*dfac_x; 
  incr[11] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[12]-67.0820393249937*alpha[8]*fr[11]-75.0*alpha[0]*fr[11]-30.0*alpha[2]*fr[8]+51.96152422706632*fr[4]*alpha[8]-30.0*fr[2]*alpha[8]-75.00000000000001*alpha[2]*fr[7]+58.09475019311126*alpha[0]*fr[4]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*fr[1]*alpha[2]-33.54101966249684*fr[0]*alpha[2])*dfac_x; 
  incr[12] = 0.005050762722761052*(67.0820393249937*alpha[8]*fr[12]+105.0*alpha[0]*fr[12]-121.2435565298214*alpha[2]*fr[11]-38.72983346207417*alpha[8]*fr[8]-60.62177826491071*alpha[0]*fr[8]-135.5544171172596*fr[7]*alpha[8]+105.0*fr[1]*alpha[8]-60.62177826491071*fr[0]*alpha[8]+93.91485505499116*alpha[2]*fr[4]-54.22176684690384*alpha[2]*fr[2])*dfac_x; 
  incr[13] = -0.05892556509887893*(11.61895003862225*alpha[8]*fr[18]-15.0*alpha[2]*fr[17]-6.708203932499369*alpha[8]*fr[14]-15.0*alpha[0]*fr[13]+11.61895003862225*alpha[2]*fr[10]-6.708203932499369*alpha[2]*fr[6]+11.61895003862225*alpha[0]*fr[5]-6.708203932499369*alpha[0]*fr[3])*dfac_x; 
  incr[14] = -0.001683587574253684*(116.1895003862225*alpha[8]*fr[18]+181.8653347947321*alpha[0]*fr[18]-210.0*alpha[2]*fr[17]-67.0820393249937*alpha[8]*fr[14]-105.0*alpha[0]*fr[14]-234.787137637478*alpha[8]*fr[13]+162.6653005407115*alpha[2]*fr[10]+181.8653347947321*fr[5]*alpha[8]-105.0*fr[3]*alpha[8]-93.91485505499116*alpha[2]*fr[6])*dfac_x; 
  incr[15] = 0.03535533905932736*(15.0*alpha[2]*fr[19]-8.660254037844386*alpha[2]*fr[16]+15.0*alpha[0]*fr[15]-8.660254037844387*alpha[0]*fr[9])*dfac_x; 
  incr[16] = -0.01178511301977579*(23.2379000772445*alpha[8]*fr[19]+25.98076211353316*alpha[0]*fr[19]-13.41640786499874*alpha[8]*fr[16]-15.0*alpha[0]*fr[16]+25.98076211353316*alpha[2]*fr[15]-15.0*alpha[2]*fr[9])*dfac_x; 
  incr[17] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[18]-67.0820393249937*alpha[8]*fr[17]-75.0*alpha[0]*fr[17]-30.0*alpha[2]*fr[14]-75.00000000000001*alpha[2]*fr[13]+51.96152422706631*alpha[8]*fr[10]+58.09475019311126*alpha[0]*fr[10]-30.0*fr[6]*alpha[8]-33.54101966249685*alpha[0]*fr[6]+58.09475019311126*alpha[2]*fr[5]-33.54101966249685*alpha[2]*fr[3])*dfac_x; 
  incr[18] = 0.005050762722761052*(67.0820393249937*alpha[8]*fr[18]+105.0*alpha[0]*fr[18]-121.2435565298214*alpha[2]*fr[17]-38.72983346207417*alpha[8]*fr[14]-60.62177826491071*alpha[0]*fr[14]-135.5544171172596*alpha[8]*fr[13]+93.91485505499116*alpha[2]*fr[10]+105.0*fr[5]*alpha[8]-60.6217782649107*fr[3]*alpha[8]-54.22176684690384*alpha[2]*fr[6])*dfac_x; 
  incr[19] = 0.03535533905932736*(13.41640786499874*alpha[8]*fr[19]+15.0*alpha[0]*fr[19]-7.745966692414834*alpha[8]*fr[16]-8.660254037844387*alpha[0]*fr[16]+15.0*alpha[2]*fr[15]-8.660254037844386*alpha[2]*fr[9])*dfac_x; 

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
double GyrokineticSurf1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.25*(((4.242640687119286*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]+4.242640687119286*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1])*dfac_v*dfac_x*dfac_x*m_*wm+(4.242640687119286*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]+4.242640687119286*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_v*dfac_x*dfac_x*m_*q_)*wv+(((-4.242640687119286*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1])-4.242640687119286*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1])*dfac_x*dfac_x*m_-4.898979485566357*Gradpar[0]*Bmag[1]*dfac_v*dfac_x*q_)*wm+((-4.242640687119286*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1])-4.242640687119286*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_x*dfac_x*m_*q_-4.898979485566357*Gradpar[0]*Phi[1]*dfac_v*dfac_x*q2))/(dfac_v*m_*q_); 

  double alpha[8]; 
  alpha[0] = (2.121320343559643*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559643*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm*wv)/q_+2.121320343559643*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+2.121320343559643*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*wv-(2.121320343559643*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559643*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.449489742783178*Gradpar[0]*Bmag[1]*dfac_x*wm)/m_-(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(2.121320343559643*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(2.121320343559643*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[1] = (2.121320343559643*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559643*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm*wv)/q_+2.121320343559643*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+2.121320343559643*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*wv-(2.121320343559643*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559643*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.449489742783178*Bmag[1]*Gradpar[1]*dfac_x*wm)/m_-(2.449489742783178*Gradpar[1]*Phi[1]*dfac_x*q_)/m_-(2.121320343559643*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(2.121320343559643*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[3] = (1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)-(1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.414213562373095*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[5] = (1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)-(1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.414213562373095*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[5]*fl[7]+1.732050807568877*alpha[3]*fl[6]+alpha[5]*fl[5]+1.732050807568877*alpha[1]*fl[4]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*alpha[3]*fl[7]+1.732050807568877*alpha[5]*fl[6]+alpha[3]*fl[5]+fl[3]*alpha[5]+1.732050807568877*alpha[0]*fl[4]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.1767766952966368*(3.0*alpha[5]*fl[7]+3.0*alpha[3]*fl[6]+1.732050807568877*alpha[5]*fl[5]+3.0*alpha[1]*fl[4]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[1]*fl[7]+1.732050807568877*alpha[0]*fl[6]+alpha[1]*fl[5]+1.732050807568877*fl[4]*alpha[5]+fl[1]*alpha[5]+alpha[0]*fl[3]+1.732050807568877*fl[2]*alpha[3]+fl[0]*alpha[3])*dfac_v; 
  incr[4] = -0.1767766952966368*(3.0*alpha[3]*fl[7]+3.0*alpha[5]*fl[6]+1.732050807568877*alpha[3]*fl[5]+1.732050807568877*fl[3]*alpha[5]+3.0*alpha[0]*fl[4]+3.0*alpha[1]*fl[2]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[6]+alpha[0]*fl[5]+1.732050807568877*fl[2]*alpha[5]+fl[0]*alpha[5]+1.732050807568877*alpha[3]*fl[4]+alpha[1]*fl[3]+fl[1]*alpha[3])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*alpha[1]*fl[7]+3.0*alpha[0]*fl[6]+1.732050807568877*alpha[1]*fl[5]+3.0*fl[4]*alpha[5]+1.732050807568877*fl[1]*alpha[5]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fl[7]+3.0*alpha[1]*fl[6]+1.732050807568877*alpha[0]*fl[5]+3.0*fl[2]*alpha[5]+1.732050807568877*fl[0]*alpha[5]+3.0*alpha[3]*fl[4]+1.732050807568877*alpha[1]*fl[3]+1.732050807568877*fl[1]*alpha[3])*dfac_v; 

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
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[5]*fr[7]+1.732050807568877*alpha[3]*fr[6]-1.0*alpha[5]*fr[5]+1.732050807568877*alpha[1]*fr[4]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*alpha[3]*fr[7]+1.732050807568877*alpha[5]*fr[6]-1.0*alpha[3]*fr[5]-1.0*fr[3]*alpha[5]+1.732050807568877*alpha[0]*fr[4]+1.732050807568877*alpha[1]*fr[2]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(3.0*alpha[5]*fr[7]+3.0*alpha[3]*fr[6]-1.732050807568877*alpha[5]*fr[5]+3.0*alpha[1]*fr[4]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[1]*fr[7]+1.732050807568877*alpha[0]*fr[6]-1.0*alpha[1]*fr[5]+1.732050807568877*fr[4]*alpha[5]-1.0*fr[1]*alpha[5]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_v; 
  incr[4] = 0.1767766952966368*(3.0*alpha[3]*fr[7]+3.0*alpha[5]*fr[6]-1.732050807568877*alpha[3]*fr[5]-1.732050807568877*fr[3]*alpha[5]+3.0*alpha[0]*fr[4]+3.0*alpha[1]*fr[2]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[7]+1.732050807568877*alpha[1]*fr[6]-1.0*alpha[0]*fr[5]+1.732050807568877*fr[2]*alpha[5]-1.0*fr[0]*alpha[5]+1.732050807568877*alpha[3]*fr[4]-1.0*alpha[1]*fr[3]-1.0*fr[1]*alpha[3])*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*alpha[1]*fr[7]+3.0*alpha[0]*fr[6]-1.732050807568877*alpha[1]*fr[5]+3.0*fr[4]*alpha[5]-1.732050807568877*fr[1]*alpha[5]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fr[7]+3.0*alpha[1]*fr[6]-1.732050807568877*alpha[0]*fr[5]+3.0*fr[2]*alpha[5]-1.732050807568877*fr[0]*alpha[5]+3.0*alpha[3]*fr[4]-1.732050807568877*alpha[1]*fr[3]-1.732050807568877*fr[1]*alpha[3])*dfac_v; 

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
double GyrokineticSurf1x2vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.05050762722761052*(((((165.0*Bmag[2]*Bmag[2]+21.0*Bmag[1]*Bmag[1])*BmagInv[2]+93.91485505499116*BmagInv[0]*Bmag[2]*Bmag[2]+84.0*Bmag[1]*BmagInv[1]*Bmag[2])*geoY[2]+(93.91485505499116*geoY[0]*Bmag[2]*Bmag[2]+84.0*Bmag[1]*geoY[1]*Bmag[2])*BmagInv[2]+(189.0*BmagInv[1]*geoY[1]+105.0*BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+(93.91485505499116*BmagInv[0]*Bmag[1]*geoY[1]+93.91485505499116*geoY[0]*Bmag[1]*BmagInv[1])*Bmag[2]+21.0*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]+21.0*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1])*dfac_v*dfac_x*dfac_x*m_*wm+(((165.0*Bmag[2]*BmagInv[2]+93.91485505499116*BmagInv[0]*Bmag[2]+42.0*Bmag[1]*BmagInv[1])*Phi[2]+21.0*Bmag[1]*Phi[1]*BmagInv[2]+42.0*BmagInv[1]*Phi[1]*Bmag[2])*geoY[2]+((93.91485505499116*geoY[0]*Bmag[2]+42.0*Bmag[1]*geoY[1])*BmagInv[2]+(189.0*BmagInv[1]*geoY[1]+105.0*BmagInv[0]*geoY[0])*Bmag[2]+46.95742752749558*BmagInv[0]*Bmag[1]*geoY[1]+46.95742752749558*geoY[0]*Bmag[1]*BmagInv[1])*Phi[2]+42.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]+(46.95742752749558*BmagInv[0]*Phi[1]*geoY[1]+46.95742752749558*geoY[0]*BmagInv[1]*Phi[1])*Bmag[2]+21.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]+21.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_v*dfac_x*dfac_x*m_*q_)*wv+(((-54.22176684690384*Gradpar[1]*Bmag[2])-24.24871130596428*Gradpar[0]*Bmag[1])*dfac_v*dfac_x*q_+((((-165.0*Bmag[2]*Bmag[2])-21.0*Bmag[1]*Bmag[1])*BmagInv[2]-93.91485505499116*BmagInv[0]*Bmag[2]*Bmag[2]-84.0*Bmag[1]*BmagInv[1]*Bmag[2])*geoY[2]+((-93.91485505499116*geoY[0]*Bmag[2]*Bmag[2])-84.0*Bmag[1]*geoY[1]*Bmag[2])*BmagInv[2]+((-189.0*BmagInv[1]*geoY[1])-105.0*BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+((-93.91485505499116*BmagInv[0]*Bmag[1]*geoY[1])-93.91485505499116*geoY[0]*Bmag[1]*BmagInv[1])*Bmag[2]-21.0*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]-21.0*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1])*dfac_x*dfac_x*m_)*wm+((((-165.0*Bmag[2]*BmagInv[2])-93.91485505499116*BmagInv[0]*Bmag[2]-42.0*Bmag[1]*BmagInv[1])*Phi[2]-21.0*Bmag[1]*Phi[1]*BmagInv[2]-42.0*BmagInv[1]*Phi[1]*Bmag[2])*geoY[2]+(((-93.91485505499116*geoY[0]*Bmag[2])-42.0*Bmag[1]*geoY[1])*BmagInv[2]+((-189.0*BmagInv[1]*geoY[1])-105.0*BmagInv[0]*geoY[0])*Bmag[2]-46.95742752749558*BmagInv[0]*Bmag[1]*geoY[1]-46.95742752749558*geoY[0]*Bmag[1]*BmagInv[1])*Phi[2]-42.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]+((-46.95742752749558*BmagInv[0]*Phi[1]*geoY[1])-46.95742752749558*geoY[0]*BmagInv[1]*Phi[1])*Bmag[2]-21.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]-21.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_x*dfac_x*m_*q_+((-54.22176684690384*Gradpar[1]*Phi[2])-24.24871130596428*Gradpar[0]*Phi[1])*dfac_v*dfac_x*q2))/(dfac_v*m_*q_); 

  double alpha[20]; 
  alpha[0] = (16.66751698511147*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(9.486832980505136*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(8.485281374238568*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(9.486832980505136*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(8.485281374238568*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(19.09188309203678*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(10.60660171779821*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(9.486832980505136*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(9.486832980505136*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm*wv)/q_+16.66751698511147*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+9.486832980505136*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+4.242640687119284*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+2.121320343559642*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+4.242640687119284*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+9.486832980505136*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119284*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+19.09188309203678*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+10.60660171779821*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+4.743416490252568*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+4.743416490252568*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119284*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+4.743416490252568*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+4.743416490252568*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+2.121320343559642*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*wv-(16.66751698511147*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505136*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.485281374238568*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505136*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.485281374238568*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(19.09188309203678*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(10.60660171779821*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505136*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505136*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(5.47722557505166*Gradpar[1]*Bmag[2]*dfac_x*wm)/m_-(2.449489742783177*Gradpar[0]*Bmag[1]*dfac_x*wm)/m_-(5.47722557505166*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(2.449489742783177*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(16.66751698511147*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(9.486832980505136*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(9.486832980505136*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(19.09188309203678*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(10.60660171779821*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.743416490252568*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.743416490252568*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(4.743416490252568*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(4.743416490252568*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[1] = (14.90788039793665*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(24.39471337844179*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(8.48528137423857*BmagInv[0]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(1.897366596101028*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(24.39471337844179*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(8.48528137423857*geoY[0]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(1.897366596101028*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(19.09188309203678*BmagInv[0]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(19.09188309203678*geoY[0]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(17.07629936490925*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(9.486832980505142*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm*wv)/q_+7.453940198968325*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+24.39471337844179*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+4.242640687119285*BmagInv[0]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+7.453940198968325*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+4.242640687119285*BmagInv[0]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+1.897366596101028*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*wv+24.39471337844179*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119285*geoY[0]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+19.09188309203678*BmagInv[0]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+19.09188309203678*geoY[0]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+8.538149682454625*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+4.743416490252571*BmagInv[0]*geoY[0]*Bmag[1]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119285*geoY[0]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+1.897366596101028*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wv+8.538149682454625*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+4.743416490252571*BmagInv[0]*geoY[0]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+2.121320343559642*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+2.121320343559642*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*wv-(14.90788039793665*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(24.39471337844179*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.48528137423857*BmagInv[0]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(1.897366596101028*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(24.39471337844179*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.48528137423857*geoY[0]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(1.897366596101028*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(19.09188309203678*BmagInv[0]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(19.09188309203678*geoY[0]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(17.07629936490925*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505142*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(4.898979485566356*Bmag[2]*Gradpar[2]*dfac_x*wm)/m_-(5.477225575051661*Gradpar[0]*Bmag[2]*dfac_x*wm)/m_-(2.449489742783178*Bmag[1]*Gradpar[1]*dfac_x*wm)/m_-(4.898979485566356*Gradpar[2]*Phi[2]*dfac_x*q_)/m_-(5.477225575051661*Gradpar[0]*Phi[2]*dfac_x*q_)/m_-(2.449489742783178*Gradpar[1]*Phi[1]*dfac_x*q_)/m_-(7.453940198968325*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(24.39471337844179*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119285*BmagInv[0]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(7.453940198968325*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119285*BmagInv[0]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(1.897366596101028*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x)/dfac_v-(24.39471337844179*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119285*geoY[0]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(19.09188309203678*BmagInv[0]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(19.09188309203678*geoY[0]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(8.538149682454625*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.743416490252571*BmagInv[0]*geoY[0]*Bmag[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119285*geoY[0]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(1.897366596101028*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(8.538149682454625*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(4.743416490252571*BmagInv[0]*geoY[0]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[3] = (9.62299541807677*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051659*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566355*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051659*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566355*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(11.0227038425243*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(6.123724356957943*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051659*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051659*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)-(9.62299541807677*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.477225575051659*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566355*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.477225575051659*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566355*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(11.0227038425243*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(6.123724356957943*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.477225575051659*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.477225575051659*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(3.162277660168379*Gradpar[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(1.414213562373095*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[5] = (8.607068760795467*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(14.08429433584713*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566356*BmagInv[0]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.095445115010332*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(14.08429433584713*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566356*geoY[0]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.095445115010332*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(11.0227038425243*BmagInv[0]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(11.0227038425243*geoY[0]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(9.859006035092989*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051662*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)-(8.607068760795467*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(14.08429433584713*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566356*BmagInv[0]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.095445115010332*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(14.08429433584713*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566356*geoY[0]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.095445115010332*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(11.0227038425243*BmagInv[0]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(11.0227038425243*geoY[0]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(9.859006035092989*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.477225575051662*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(2.82842712474619*Bmag[2]*Gradpar[2]*dfac_x)/(dfac_m*m_)-(3.16227766016838*Gradpar[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(1.414213562373095*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  alpha[7] = (27.10523708715754*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(1.355261854357877*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(16.66751698511148*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(14.90788039793665*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(16.66751698511148*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(14.90788039793665*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(24.39471337844179*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(9.486832980505142*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(8.48528137423857*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(8.48528137423857*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(1.897366596101028*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm*wv)/q_+27.10523708715754*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+16.66751698511148*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+7.453940198968325*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+1.355261854357877*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+7.453940198968325*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+2.121320343559642*BmagInv[0]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*wv+16.66751698511148*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+7.453940198968325*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+24.39471337844179*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+9.486832980505142*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119285*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119285*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*wv+7.453940198968325*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+2.121320343559642*geoY[0]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x*wv+4.242640687119285*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+4.242640687119285*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+1.897366596101028*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv-(27.10523708715754*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(1.355261854357877*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(16.66751698511148*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(14.90788039793665*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(16.66751698511148*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(14.90788039793665*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(24.39471337844179*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505142*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.48528137423857*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.48528137423857*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(1.897366596101028*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.449489742783178*Bmag[1]*Gradpar[2]*dfac_x*wm)/m_-(4.898979485566356*Gradpar[1]*Bmag[2]*dfac_x*wm)/m_-(4.898979485566356*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(2.449489742783178*Phi[1]*Gradpar[2]*dfac_x*q_)/m_-(27.10523708715754*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(16.66751698511148*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(7.453940198968325*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(1.355261854357877*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(7.453940198968325*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*BmagInv[0]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x)/dfac_v-(16.66751698511148*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(7.453940198968325*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(24.39471337844179*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(9.486832980505142*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119285*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119285*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(7.453940198968325*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*geoY[0]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119285*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119285*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(1.897366596101028*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[13] = (15.64921592871903*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(0.7824607964359517*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(9.622995418076774*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(8.60706876079547*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(9.622995418076774*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(8.60706876079547*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(14.08429433584713*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051661*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566358*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566358*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.095445115010332*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)-(15.64921592871903*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(0.7824607964359517*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(9.622995418076774*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(8.60706876079547*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(9.622995418076774*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(8.60706876079547*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(14.08429433584713*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.477225575051661*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566358*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566358*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.095445115010332*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.414213562373095*Bmag[1]*Gradpar[2]*dfac_x)/(dfac_m*m_)-(2.82842712474619*Gradpar[1]*Bmag[2]*dfac_x)/(dfac_m*m_); 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-1.118033988749895*fr[8]*amax)+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax+0.3952847075210473*alpha[5]*fr[18]+0.3952847075210473*alpha[5]*fl[18]-0.3061862178478971*alpha[13]*fr[17]+0.3061862178478971*alpha[13]*fl[17]+0.3952847075210473*alpha[3]*fr[14]+0.3952847075210473*alpha[3]*fl[14]+0.1767766952966368*alpha[13]*fr[13]+0.1767766952966368*alpha[13]*fl[13]+0.3952847075210473*alpha[1]*fr[12]+0.3952847075210473*alpha[1]*fl[12]-0.3061862178478971*alpha[7]*fr[11]+0.3061862178478971*alpha[7]*fl[11]-0.3061862178478971*alpha[5]*fr[10]+0.3061862178478971*alpha[5]*fl[10]+0.3952847075210473*alpha[0]*fr[8]+0.3952847075210473*alpha[0]*fl[8]+0.1767766952966368*alpha[7]*fr[7]+0.1767766952966368*alpha[7]*fl[7]-0.3061862178478971*alpha[3]*fr[6]+0.3061862178478971*alpha[3]*fl[6]+0.1767766952966368*alpha[5]*fr[5]+0.1767766952966368*alpha[5]*fl[5]-0.3061862178478971*alpha[1]*fr[4]+0.3061862178478971*alpha[1]*fl[4]+0.1767766952966368*alpha[3]*fr[3]+0.1767766952966368*alpha[3]*fl[3]-0.3061862178478971*alpha[0]*fr[2]+0.3061862178478971*alpha[0]*fl[2]+0.1767766952966368*alpha[1]*fr[1]+0.1767766952966368*alpha[1]*fl[1]+0.1767766952966368*alpha[0]*fr[0]+0.1767766952966368*alpha[0]*fl[0]; 
  Ghat[1] = (-1.118033988749895*fr[12]*amax)+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax+0.3535533905932737*alpha[13]*fr[18]+0.3952847075210473*alpha[3]*fr[18]+0.3535533905932737*alpha[13]*fl[18]+0.3952847075210473*alpha[3]*fl[18]-0.273861278752583*alpha[5]*fr[17]+0.273861278752583*alpha[5]*fl[17]+0.3952847075210473*alpha[5]*fr[14]+0.3952847075210473*alpha[5]*fl[14]+0.1581138830084189*alpha[5]*fr[13]+0.1581138830084189*alpha[5]*fl[13]-0.273861278752583*fr[10]*alpha[13]+0.273861278752583*fl[10]*alpha[13]+0.1581138830084189*fr[5]*alpha[13]+0.1581138830084189*fl[5]*alpha[13]+0.3535533905932737*alpha[7]*fr[12]+0.3952847075210473*alpha[0]*fr[12]+0.3535533905932737*alpha[7]*fl[12]+0.3952847075210473*alpha[0]*fl[12]-0.273861278752583*alpha[1]*fr[11]+0.273861278752583*alpha[1]*fl[11]-0.3061862178478971*alpha[3]*fr[10]+0.3061862178478971*alpha[3]*fl[10]+0.3952847075210473*alpha[1]*fr[8]+0.3952847075210473*alpha[1]*fl[8]+0.1581138830084189*alpha[1]*fr[7]+0.1581138830084189*alpha[1]*fl[7]-0.273861278752583*fr[4]*alpha[7]+0.273861278752583*fl[4]*alpha[7]+0.1581138830084189*fr[1]*alpha[7]+0.1581138830084189*fl[1]*alpha[7]-0.3061862178478971*alpha[5]*fr[6]+0.3061862178478971*alpha[5]*fl[6]+0.1767766952966368*alpha[3]*fr[5]+0.1767766952966368*alpha[3]*fl[5]+0.1767766952966368*fr[3]*alpha[5]+0.1767766952966368*fl[3]*alpha[5]-0.3061862178478971*alpha[0]*fr[4]+0.3061862178478971*alpha[0]*fl[4]-0.3061862178478971*alpha[1]*fr[2]+0.3061862178478971*alpha[1]*fl[2]+0.1767766952966368*alpha[0]*fr[1]+0.1767766952966368*alpha[0]*fl[1]+0.1767766952966368*fr[0]*alpha[1]+0.1767766952966368*fl[0]*alpha[1]; 
  Ghat[3] = (-1.118033988749895*fr[14]*amax)+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax-0.273861278752583*alpha[5]*fr[19]+0.273861278752583*alpha[5]*fl[19]+0.3952847075210473*alpha[1]*fr[18]+0.3952847075210473*alpha[1]*fl[18]-0.3061862178478971*alpha[7]*fr[17]+0.3061862178478971*alpha[7]*fl[17]-0.273861278752583*alpha[3]*fr[16]+0.273861278752583*alpha[3]*fl[16]+0.1581138830084189*alpha[5]*fr[15]+0.1581138830084189*alpha[5]*fl[15]+0.3952847075210473*alpha[0]*fr[14]+0.3952847075210473*alpha[0]*fl[14]+0.1767766952966368*alpha[7]*fr[13]+0.1767766952966368*alpha[7]*fl[13]-0.3061862178478971*fr[11]*alpha[13]+0.3061862178478971*fl[11]*alpha[13]+0.1767766952966368*fr[7]*alpha[13]+0.1767766952966368*fl[7]*alpha[13]+0.3952847075210473*alpha[5]*fr[12]+0.3952847075210473*alpha[5]*fl[12]-0.3061862178478971*alpha[1]*fr[10]+0.3061862178478971*alpha[1]*fl[10]+0.1581138830084189*alpha[3]*fr[9]+0.1581138830084189*alpha[3]*fl[9]+0.3952847075210473*alpha[3]*fr[8]+0.3952847075210473*alpha[3]*fl[8]-0.3061862178478971*alpha[0]*fr[6]+0.3061862178478971*alpha[0]*fl[6]+0.1767766952966368*alpha[1]*fr[5]+0.1767766952966368*alpha[1]*fl[5]-0.3061862178478971*fr[4]*alpha[5]+0.3061862178478971*fl[4]*alpha[5]+0.1767766952966368*fr[1]*alpha[5]+0.1767766952966368*fl[1]*alpha[5]+0.1767766952966368*alpha[0]*fr[3]+0.1767766952966368*alpha[0]*fl[3]-0.3061862178478971*fr[2]*alpha[3]+0.3061862178478971*fl[2]*alpha[3]+0.1767766952966368*fr[0]*alpha[3]+0.1767766952966368*fl[0]*alpha[3]; 
  Ghat[5] = (-1.118033988749895*fr[18]*amax)+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax-0.2449489742783177*alpha[13]*fr[19]-0.273861278752583*alpha[3]*fr[19]+0.2449489742783177*alpha[13]*fl[19]+0.273861278752583*alpha[3]*fl[19]+0.3535533905932737*alpha[7]*fr[18]+0.3952847075210473*alpha[0]*fr[18]+0.3535533905932737*alpha[7]*fl[18]+0.3952847075210473*alpha[0]*fl[18]-0.273861278752583*alpha[1]*fr[17]+0.273861278752583*alpha[1]*fl[17]-0.2738612787525829*alpha[5]*fr[16]+0.2738612787525829*alpha[5]*fl[16]+0.1414213562373095*alpha[13]*fr[15]+0.1581138830084189*alpha[3]*fr[15]+0.1414213562373095*alpha[13]*fl[15]+0.1581138830084189*alpha[3]*fl[15]+0.3952847075210473*alpha[1]*fr[14]+0.3952847075210473*alpha[1]*fl[14]+0.1581138830084189*alpha[1]*fr[13]+0.1581138830084189*alpha[1]*fl[13]+0.3535533905932737*fr[12]*alpha[13]+0.3535533905932737*fl[12]*alpha[13]-0.2738612787525829*fr[4]*alpha[13]+0.2738612787525829*fl[4]*alpha[13]+0.1581138830084189*fr[1]*alpha[13]+0.1581138830084189*fl[1]*alpha[13]+0.3952847075210473*alpha[3]*fr[12]+0.3952847075210473*alpha[3]*fl[12]-0.2738612787525829*alpha[5]*fr[11]+0.2738612787525829*alpha[5]*fl[11]-0.273861278752583*alpha[7]*fr[10]-0.3061862178478971*alpha[0]*fr[10]+0.273861278752583*alpha[7]*fl[10]+0.3061862178478971*alpha[0]*fl[10]+0.1581138830084189*alpha[5]*fr[9]+0.1581138830084189*alpha[5]*fl[9]+0.3952847075210473*alpha[5]*fr[8]+0.3952847075210473*alpha[5]*fl[8]+0.1581138830084189*alpha[5]*fr[7]+0.1581138830084189*alpha[5]*fl[7]+0.1581138830084189*fr[5]*alpha[7]+0.1581138830084189*fl[5]*alpha[7]-0.3061862178478971*alpha[1]*fr[6]+0.3061862178478971*alpha[1]*fl[6]+0.1767766952966368*alpha[0]*fr[5]+0.1767766952966368*alpha[0]*fl[5]-0.3061862178478971*fr[2]*alpha[5]+0.3061862178478971*fl[2]*alpha[5]+0.1767766952966368*fr[0]*alpha[5]+0.1767766952966368*fl[0]*alpha[5]-0.3061862178478971*alpha[3]*fr[4]+0.3061862178478971*alpha[3]*fl[4]+0.1767766952966368*alpha[1]*fr[3]+0.1767766952966368*alpha[1]*fl[3]+0.1767766952966368*fr[1]*alpha[3]+0.1767766952966368*fl[1]*alpha[3]; 
  Ghat[7] = 0.8660254037844386*fr[11]*amax+0.8660254037844386*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax+0.3535533905932737*alpha[5]*fr[18]+0.3535533905932737*alpha[5]*fl[18]-0.1956151991089878*alpha[13]*fr[17]-0.3061862178478971*alpha[3]*fr[17]+0.1956151991089878*alpha[13]*fl[17]+0.3061862178478971*alpha[3]*fl[17]+0.3952847075210473*alpha[13]*fr[14]+0.3952847075210473*alpha[13]*fl[14]+0.1129384878631564*alpha[13]*fr[13]+0.1767766952966368*alpha[3]*fr[13]+0.1129384878631564*alpha[13]*fl[13]+0.1767766952966368*alpha[3]*fl[13]-0.3061862178478972*fr[6]*alpha[13]+0.3061862178478972*fl[6]*alpha[13]+0.1767766952966368*fr[3]*alpha[13]+0.1767766952966368*fl[3]*alpha[13]+0.3535533905932737*alpha[1]*fr[12]+0.3535533905932737*alpha[1]*fl[12]-0.1956151991089878*alpha[7]*fr[11]-0.3061862178478972*alpha[0]*fr[11]+0.1956151991089878*alpha[7]*fl[11]+0.3061862178478972*alpha[0]*fl[11]-0.273861278752583*alpha[5]*fr[10]+0.273861278752583*alpha[5]*fl[10]+0.3952847075210473*alpha[7]*fr[8]+0.3952847075210473*alpha[7]*fl[8]+0.1129384878631564*alpha[7]*fr[7]+0.1767766952966368*alpha[0]*fr[7]+0.1129384878631564*alpha[7]*fl[7]+0.1767766952966368*alpha[0]*fl[7]-0.3061862178478971*fr[2]*alpha[7]+0.3061862178478971*fl[2]*alpha[7]+0.1767766952966368*fr[0]*alpha[7]+0.1767766952966368*fl[0]*alpha[7]+0.1581138830084189*alpha[5]*fr[5]+0.1581138830084189*alpha[5]*fl[5]-0.273861278752583*alpha[1]*fr[4]+0.273861278752583*alpha[1]*fl[4]+0.1581138830084189*alpha[1]*fr[1]+0.1581138830084189*alpha[1]*fl[1]; 
  Ghat[9] = 0.8660254037844386*fr[16]*amax+0.8660254037844386*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax-0.3061862178478971*alpha[1]*fr[19]+0.3061862178478971*alpha[1]*fl[19]+0.3535533905932737*alpha[5]*fr[18]+0.3535533905932737*alpha[5]*fl[18]-0.273861278752583*alpha[13]*fr[17]+0.273861278752583*alpha[13]*fl[17]-0.3061862178478972*alpha[0]*fr[16]+0.3061862178478972*alpha[0]*fl[16]+0.1767766952966368*alpha[1]*fr[15]+0.1767766952966368*alpha[1]*fl[15]+0.3535533905932737*alpha[3]*fr[14]+0.3535533905932737*alpha[3]*fl[14]+0.1581138830084189*alpha[13]*fr[13]+0.1581138830084189*alpha[13]*fl[13]-0.273861278752583*alpha[5]*fr[10]+0.273861278752583*alpha[5]*fl[10]+0.1767766952966368*alpha[0]*fr[9]+0.1767766952966368*alpha[0]*fl[9]-0.273861278752583*alpha[3]*fr[6]+0.273861278752583*alpha[3]*fl[6]+0.1581138830084189*alpha[5]*fr[5]+0.1581138830084189*alpha[5]*fl[5]+0.1581138830084189*alpha[3]*fr[3]+0.1581138830084189*alpha[3]*fl[3]; 
  Ghat[13] = 0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax-0.2449489742783177*alpha[5]*fr[19]+0.2449489742783177*alpha[5]*fl[19]+0.3535533905932737*alpha[1]*fr[18]+0.3535533905932737*alpha[1]*fl[18]-0.1956151991089878*alpha[7]*fr[17]-0.3061862178478971*alpha[0]*fr[17]+0.1956151991089878*alpha[7]*fl[17]+0.3061862178478971*alpha[0]*fl[17]-0.2738612787525829*alpha[13]*fr[16]+0.2738612787525829*alpha[13]*fl[16]+0.1414213562373095*alpha[5]*fr[15]+0.1414213562373095*alpha[5]*fl[15]+0.3952847075210473*alpha[7]*fr[14]+0.3952847075210473*alpha[7]*fl[14]+0.1129384878631564*alpha[7]*fr[13]+0.1767766952966368*alpha[0]*fr[13]+0.1129384878631564*alpha[7]*fl[13]+0.1767766952966368*alpha[0]*fl[13]-0.1956151991089878*fr[11]*alpha[13]+0.1956151991089878*fl[11]*alpha[13]+0.1581138830084189*fr[9]*alpha[13]+0.1581138830084189*fl[9]*alpha[13]+0.3952847075210473*fr[8]*alpha[13]+0.3952847075210473*fl[8]*alpha[13]+0.1129384878631564*fr[7]*alpha[13]+0.1129384878631564*fl[7]*alpha[13]-0.3061862178478971*fr[2]*alpha[13]+0.3061862178478971*fl[2]*alpha[13]+0.1767766952966368*fr[0]*alpha[13]+0.1767766952966368*fl[0]*alpha[13]+0.3535533905932737*alpha[5]*fr[12]+0.3535533905932737*alpha[5]*fl[12]-0.3061862178478971*alpha[3]*fr[11]+0.3061862178478971*alpha[3]*fl[11]-0.2738612787525829*alpha[1]*fr[10]+0.2738612787525829*alpha[1]*fl[10]+0.1767766952966368*alpha[3]*fr[7]+0.1767766952966368*alpha[3]*fl[7]-0.3061862178478971*fr[6]*alpha[7]+0.3061862178478971*fl[6]*alpha[7]+0.1767766952966368*fr[3]*alpha[7]+0.1767766952966368*fl[3]*alpha[7]+0.1581138830084189*alpha[1]*fr[5]+0.1581138830084189*alpha[1]*fl[5]-0.2738612787525829*fr[4]*alpha[5]+0.2738612787525829*fl[4]*alpha[5]+0.1581138830084189*fr[1]*alpha[5]+0.1581138830084189*fl[1]*alpha[5]; 
  Ghat[15] = 0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax-0.2738612787525829*alpha[7]*fr[19]-0.3061862178478971*alpha[0]*fr[19]+0.2738612787525829*alpha[7]*fl[19]+0.3061862178478971*alpha[0]*fl[19]+0.3162277660168379*alpha[13]*fr[18]+0.3535533905932737*alpha[3]*fr[18]+0.3162277660168379*alpha[13]*fl[18]+0.3535533905932737*alpha[3]*fl[18]-0.2449489742783177*alpha[5]*fr[17]+0.2449489742783177*alpha[5]*fl[17]-0.3061862178478971*alpha[1]*fr[16]+0.3061862178478971*alpha[1]*fl[16]+0.1581138830084189*alpha[7]*fr[15]+0.1767766952966368*alpha[0]*fr[15]+0.1581138830084189*alpha[7]*fl[15]+0.1767766952966368*alpha[0]*fl[15]+0.3535533905932737*alpha[5]*fr[14]+0.3535533905932737*alpha[5]*fl[14]+0.1414213562373095*alpha[5]*fr[13]+0.1414213562373095*alpha[5]*fl[13]-0.2449489742783178*fr[10]*alpha[13]+0.2449489742783178*fl[10]*alpha[13]+0.1414213562373095*fr[5]*alpha[13]+0.1414213562373095*fl[5]*alpha[13]-0.2738612787525829*alpha[3]*fr[10]+0.2738612787525829*alpha[3]*fl[10]+0.1767766952966368*alpha[1]*fr[9]+0.1767766952966368*alpha[1]*fl[9]-0.2738612787525829*alpha[5]*fr[6]+0.2738612787525829*alpha[5]*fl[6]+0.1581138830084189*alpha[3]*fr[5]+0.1581138830084189*alpha[3]*fl[5]+0.1581138830084189*fr[3]*alpha[5]+0.1581138830084189*fl[3]*alpha[5]; 

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
double EmGyrokineticSurf1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.3535533905932737*(((7.348469228349534*Bmag[1]*BmagInv[1]-4.242640687119286*BmagInv[0]*Bmag[1])*geoY[1]-4.242640687119286*geoY[0]*Bmag[1]*BmagInv[1]+2.449489742783178*BmagInv[0]*geoY[0]*Bmag[1])*dfac_x*m_*wv2+((((5.196152422706631*Apar[0]-9.0*Apar[1])*Bmag[1]*BmagInv[1]+(5.196152422706631*BmagInv[0]*Apar[1]-3.0*Apar[0]*BmagInv[0])*Bmag[1]+6.0*Apar[1])*geoY[1]+(5.196152422706631*geoY[0]*Apar[1]-3.0*Apar[0]*geoY[0])*Bmag[1]*BmagInv[1]+(1.732050807568877*Apar[0]*BmagInv[0]*geoY[0]-3.0*BmagInv[0]*geoY[0]*Apar[1])*Bmag[1]-3.464101615137754*geoY[0]*Apar[1])*dfac_x+4.898979485566357*Gradpar[1]-2.828427124746191*Gradpar[0])*q_*wv))/q_; 

  double alpha[8]; 
  alpha[0] = (-(5.196152422706631*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv2)/q_)+(3.0*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv2)/q_+(3.0*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv2)/q_-(1.732050807568877*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv2)/q_+6.363961030678928*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*wv-3.674234614174767*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*wv-3.674234614174767*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*dfac_x*wv+2.121320343559643*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*wv-4.242640687119286*Apar[1]*geoY[1]*dfac_x*wv-3.674234614174767*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*dfac_x*wv+2.121320343559643*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*wv+2.121320343559643*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*dfac_x*wv-1.224744871391589*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*wv+2.449489742783178*geoY[0]*Apar[1]*dfac_x*wv-3.464101615137754*Gradpar[1]*wv+2.0*Gradpar[0]*wv; 
  alpha[2] = (-(3.0*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_))+(1.732050807568877*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(1.732050807568877*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv)/(dfac_v*q_)-(1.0*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[2]*fl[4]+alpha[2]*fl[2]+1.732050807568877*alpha[0]*fl[1]+alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*(3.0*alpha[2]*fl[4]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[4]+alpha[0]*fl[2]+1.732050807568877*fl[1]*alpha[2]+fl[0]*alpha[2])*dfac_x; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[2]*fl[7]+alpha[2]*fl[6]+1.732050807568877*alpha[0]*fl[5]+alpha[0]*fl[3])*dfac_x; 
  incr[4] = -0.1767766952966368*(3.0*alpha[0]*fl[4]+1.732050807568877*alpha[0]*fl[2]+3.0*fl[1]*alpha[2]+1.732050807568877*fl[0]*alpha[2])*dfac_x; 
  incr[5] = -0.1767766952966368*(3.0*alpha[2]*fl[7]+1.732050807568877*alpha[2]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[3])*dfac_x; 
  incr[6] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[7]+alpha[0]*fl[6]+1.732050807568877*alpha[2]*fl[5]+alpha[2]*fl[3])*dfac_x; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fl[7]+1.732050807568877*alpha[0]*fl[6]+3.0*alpha[2]*fl[5]+1.732050807568877*alpha[2]*fl[3])*dfac_x; 

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
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[2]*fr[4]-1.0*alpha[2]*fr[2]+1.732050807568877*alpha[0]*fr[1]-1.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.1767766952966368*(3.0*alpha[2]*fr[4]-1.732050807568877*alpha[2]*fr[2]+3.0*alpha[0]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[4]-1.0*alpha[0]*fr[2]+1.732050807568877*fr[1]*alpha[2]-1.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[2]*fr[7]-1.0*alpha[2]*fr[6]+1.732050807568877*alpha[0]*fr[5]-1.0*alpha[0]*fr[3])*dfac_x; 
  incr[4] = 0.1767766952966368*(3.0*alpha[0]*fr[4]-1.732050807568877*alpha[0]*fr[2]+3.0*fr[1]*alpha[2]-1.732050807568877*fr[0]*alpha[2])*dfac_x; 
  incr[5] = 0.1767766952966368*(3.0*alpha[2]*fr[7]-1.732050807568877*alpha[2]*fr[6]+3.0*alpha[0]*fr[5]-1.732050807568877*alpha[0]*fr[3])*dfac_x; 
  incr[6] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[7]-1.0*alpha[0]*fr[6]+1.732050807568877*alpha[2]*fr[5]-1.0*alpha[2]*fr[3])*dfac_x; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fr[7]-1.732050807568877*alpha[0]*fr[6]+3.0*alpha[2]*fr[5]-1.732050807568877*alpha[2]*fr[3])*dfac_x; 

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
double EmGyrokineticSurf1x2vSer_X_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.1178511301977579*((((142.3024947075771*Bmag[2]-36.74234614174767*Bmag[1])*BmagInv[2]+(63.63961030678928*BmagInv[0]-110.227038425243*BmagInv[1])*Bmag[2]+28.46049894151542*Bmag[1]*BmagInv[1]-16.43167672515498*BmagInv[0]*Bmag[1])*geoY[2]+((63.63961030678928*geoY[0]-110.227038425243*geoY[1])*Bmag[2]+28.46049894151542*Bmag[1]*geoY[1]-16.43167672515498*geoY[0]*Bmag[1])*BmagInv[2]+((85.38149682454625*BmagInv[1]-49.29503017546494*BmagInv[0])*geoY[1]-49.29503017546494*geoY[0]*BmagInv[1]+28.46049894151542*BmagInv[0]*geoY[0])*Bmag[2]+(12.72792206135786*BmagInv[0]*Bmag[1]-22.0454076850486*Bmag[1]*BmagInv[1])*geoY[1]+12.72792206135786*geoY[0]*Bmag[1]*BmagInv[1]-7.348469228349534*BmagInv[0]*geoY[0]*Bmag[1])*dfac_v2*dfac_x*m_*wv2+(((((225.0*Apar[2]-174.2842505793337*Apar[1]+100.6230589874906*Apar[0])*Bmag[2]-58.09475019311126*Bmag[1]*Apar[2]+(45.0*Apar[1]-25.98076211353316*Apar[0])*Bmag[1])*BmagInv[2]+((100.6230589874906*BmagInv[0]-174.2842505793337*BmagInv[1])*Apar[2]+(135.0*Apar[1]-77.94228634059945*Apar[0])*BmagInv[1]-77.94228634059945*BmagInv[0]*Apar[1]+45.0*Apar[0]*BmagInv[0])*Bmag[2]+(45.0*Bmag[1]*BmagInv[1]-25.98076211353316*BmagInv[0]*Bmag[1]-90.0)*Apar[2]+(20.12461179749811*Apar[0]-34.85685011586674*Apar[1])*Bmag[1]*BmagInv[1]+(20.12461179749811*BmagInv[0]*Apar[1]-11.61895003862225*Apar[0]*BmagInv[0])*Bmag[1]+23.2379000772445*Apar[1])*geoY[2]+(((100.6230589874906*geoY[0]-174.2842505793337*geoY[1])*Apar[2]+(135.0*Apar[1]-77.94228634059945*Apar[0])*geoY[1]-77.94228634059945*geoY[0]*Apar[1]+45.0*Apar[0]*geoY[0])*Bmag[2]+(45.0*Bmag[1]*geoY[1]-25.98076211353316*geoY[0]*Bmag[1])*Apar[2]+(20.12461179749811*Apar[0]-34.85685011586674*Apar[1])*Bmag[1]*geoY[1]+(20.12461179749811*geoY[0]*Apar[1]-11.61895003862225*Apar[0]*geoY[0])*Bmag[1])*BmagInv[2]+(((135.0*BmagInv[1]-77.94228634059945*BmagInv[0])*geoY[1]-77.94228634059945*geoY[0]*BmagInv[1]+45.0*BmagInv[0]*geoY[0])*Apar[2]+((60.37383539249433*Apar[0]-104.5705503476002*Apar[1])*BmagInv[1]+60.37383539249433*BmagInv[0]*Apar[1]-34.85685011586674*Apar[0]*BmagInv[0])*geoY[1]+(60.37383539249433*geoY[0]*Apar[1]-34.85685011586674*Apar[0]*geoY[0])*BmagInv[1]-34.85685011586674*BmagInv[0]*geoY[0]*Apar[1]+20.12461179749811*Apar[0]*BmagInv[0]*geoY[0])*Bmag[2]+(((-34.85685011586674*Bmag[1]*BmagInv[1])+20.12461179749811*BmagInv[0]*Bmag[1]+69.71370023173348)*geoY[1]+20.12461179749811*geoY[0]*Bmag[1]*BmagInv[1]-11.61895003862225*BmagInv[0]*geoY[0]*Bmag[1]-40.24922359499622*geoY[0])*Apar[2]+((27.0*Apar[1]-15.58845726811989*Apar[0])*Bmag[1]*BmagInv[1]+(9.0*Apar[0]*BmagInv[0]-15.58845726811989*BmagInv[0]*Apar[1])*Bmag[1]-18.0*Apar[1])*geoY[1]+(9.0*Apar[0]*geoY[0]-15.58845726811989*geoY[0]*Apar[1])*Bmag[1]*BmagInv[1]+(9.0*BmagInv[0]*geoY[0]*Apar[1]-5.196152422706631*Apar[0]*BmagInv[0]*geoY[0])*Bmag[1]+10.39230484541326*geoY[0]*Apar[1])*dfac_v2*dfac_x+(18.97366596101028*Gradpar[2]-14.69693845669907*Gradpar[1]+8.485281374238571*Gradpar[0])*dfac_v2)*q_*wv+(((47.43416490252571*Bmag[2]-12.24744871391589*Bmag[1])*BmagInv[2]+(21.21320343559643*BmagInv[0]-36.74234614174767*BmagInv[1])*Bmag[2]+9.48683298050514*Bmag[1]*BmagInv[1]-5.477225575051662*BmagInv[0]*Bmag[1])*geoY[2]+((21.21320343559643*geoY[0]-36.74234614174767*geoY[1])*Bmag[2]+9.48683298050514*Bmag[1]*geoY[1]-5.477225575051662*geoY[0]*Bmag[1])*BmagInv[2]+((28.46049894151542*BmagInv[1]-16.43167672515498*BmagInv[0])*geoY[1]-16.43167672515498*geoY[0]*BmagInv[1]+9.48683298050514*BmagInv[0]*geoY[0])*Bmag[2]+(4.242640687119286*BmagInv[0]*Bmag[1]-7.348469228349534*Bmag[1]*BmagInv[1])*geoY[1]+4.242640687119286*geoY[0]*Bmag[1]*BmagInv[1]-2.449489742783178*BmagInv[0]*geoY[0]*Bmag[1])*dfac_x*m_))/(dfac_v2*q_); 

  double alpha[20]; 
  alpha[0] = (33.54101966249685*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_*wv2)/q_-(8.660254037844386*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_*wv2)/q_-(25.98076211353315*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_*wv2)/q_+(15.0*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_*wv2)/q_+(6.708203932499369*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_*wv2)/q_-(3.872983346207417*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_*wv2)/q_-(25.98076211353315*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv2)/q_+(15.0*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv2)/q_+(6.708203932499369*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_*wv2)/q_-(3.872983346207417*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_*wv2)/q_+(20.12461179749811*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_*wv2)/q_-(11.61895003862225*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_*wv2)/q_-(11.61895003862225*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_*wv2)/q_+(6.708203932499369*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_*wv2)/q_-(5.19615242270663*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv2)/q_+(3.0*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv2)/q_+(3.0*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv2)/q_-(1.732050807568877*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv2)/q_+53.03300858899107*Apar[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*wv-41.07919181288746*Apar[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*wv+23.71708245126285*Apar[0]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*wv-13.69306393762916*Bmag[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x*wv+10.60660171779821*Apar[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*wv-6.123724356957945*Apar[0]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*wv-41.07919181288746*BmagInv[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*wv+23.71708245126285*BmagInv[0]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*wv+31.81980515339464*Apar[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*wv-18.37117307087383*Apar[0]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*wv-18.37117307087383*BmagInv[0]*Apar[1]*Bmag[2]*geoY[2]*dfac_x*wv+10.60660171779821*Apar[0]*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*wv+10.60660171779821*Bmag[1]*BmagInv[1]*Apar[2]*geoY[2]*dfac_x*wv-6.123724356957945*BmagInv[0]*Bmag[1]*Apar[2]*geoY[2]*dfac_x*wv-21.21320343559643*Apar[2]*geoY[2]*dfac_x*wv-8.21583836257749*Apar[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*wv+4.743416490252569*Apar[0]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*wv+4.743416490252569*BmagInv[0]*Apar[1]*Bmag[1]*geoY[2]*dfac_x*wv-2.738612787525831*Apar[0]*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*wv+5.477225575051662*Apar[1]*geoY[2]*dfac_x*wv-41.07919181288746*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*wv+23.71708245126285*geoY[0]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*wv+31.81980515339464*Apar[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*wv-18.37117307087383*Apar[0]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*wv-18.37117307087383*geoY[0]*Apar[1]*Bmag[2]*BmagInv[2]*dfac_x*wv+10.60660171779821*Apar[0]*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*wv+10.60660171779821*Bmag[1]*geoY[1]*Apar[2]*BmagInv[2]*dfac_x*wv-6.123724356957945*geoY[0]*Bmag[1]*Apar[2]*BmagInv[2]*dfac_x*wv-8.21583836257749*Apar[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*wv+4.743416490252569*Apar[0]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*wv+4.743416490252569*geoY[0]*Apar[1]*Bmag[1]*BmagInv[2]*dfac_x*wv-2.738612787525831*Apar[0]*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*wv+31.81980515339464*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*wv-18.37117307087383*BmagInv[0]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*wv-18.37117307087383*geoY[0]*BmagInv[1]*Apar[2]*Bmag[2]*dfac_x*wv+10.60660171779821*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*dfac_x*wv-24.64751508773247*Apar[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*wv+14.23024947075771*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*wv+14.23024947075771*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*dfac_x*wv-8.21583836257749*Apar[0]*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*wv+14.23024947075771*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*dfac_x*wv-8.21583836257749*Apar[0]*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*wv-8.21583836257749*BmagInv[0]*geoY[0]*Apar[1]*Bmag[2]*dfac_x*wv+4.743416490252569*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*wv-8.21583836257749*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*dfac_x*wv+4.743416490252569*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*dfac_x*wv+16.43167672515499*geoY[1]*Apar[2]*dfac_x*wv+4.743416490252569*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*dfac_x*wv-2.738612787525831*BmagInv[0]*geoY[0]*Bmag[1]*Apar[2]*dfac_x*wv-9.486832980505138*geoY[0]*Apar[2]*dfac_x*wv+6.363961030678928*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*wv-3.674234614174767*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*wv-3.674234614174767*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*dfac_x*wv+2.121320343559643*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*wv-4.242640687119286*Apar[1]*geoY[1]*dfac_x*wv-3.674234614174767*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*dfac_x*wv+2.121320343559643*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*wv+2.121320343559643*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*dfac_x*wv-1.224744871391589*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*wv+2.449489742783178*geoY[0]*Apar[1]*dfac_x*wv+4.47213595499958*Gradpar[2]*wv-3.464101615137754*Gradpar[1]*wv+2.0*Gradpar[0]*wv+(11.18033988749895*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(2.886751345948128*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(8.660254037844386*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(5.0*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(2.23606797749979*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(1.290994448735806*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(8.660254037844386*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(5.0*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(2.23606797749979*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)-(1.290994448735806*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(6.708203932499369*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(3.872983346207417*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(3.872983346207417*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)+(2.23606797749979*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(1.732050807568877*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(1.0*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(1.0*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_)/(dfac_v2*q_)-(0.5773502691896257*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_)/(dfac_v2*q_); 
  alpha[2] = (38.72983346207418*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(10.0*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(30.0*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)+(17.32050807568877*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)+(7.745966692414834*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(4.47213595499958*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(30.0*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(17.32050807568877*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(7.745966692414834*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)-(4.47213595499958*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(23.2379000772445*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(13.41640786499874*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(13.41640786499874*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)+(7.745966692414834*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(6.0*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(3.464101615137754*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(3.464101615137754*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv)/(dfac_v*q_)-(2.0*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv)/(dfac_v*q_)+(30.61862178478972*Apar[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v-(23.71708245126285*Apar[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v+(13.69306393762916*Apar[0]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v-(7.905694150420951*Bmag[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v+(6.123724356957945*Apar[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v-(3.535533905932738*Apar[0]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v-(23.71708245126285*BmagInv[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x)/dfac_v+(13.69306393762916*BmagInv[0]*Apar[2]*Bmag[2]*geoY[2]*dfac_x)/dfac_v+(18.37117307087383*Apar[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x)/dfac_v-(10.60660171779821*Apar[0]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x)/dfac_v-(10.60660171779821*BmagInv[0]*Apar[1]*Bmag[2]*geoY[2]*dfac_x)/dfac_v+(6.123724356957945*Apar[0]*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x)/dfac_v+(6.123724356957945*Bmag[1]*BmagInv[1]*Apar[2]*geoY[2]*dfac_x)/dfac_v-(3.535533905932738*BmagInv[0]*Bmag[1]*Apar[2]*geoY[2]*dfac_x)/dfac_v-(12.24744871391589*Apar[2]*geoY[2]*dfac_x)/dfac_v-(4.743416490252569*Apar[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x)/dfac_v+(2.738612787525831*Apar[0]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x)/dfac_v+(2.738612787525831*BmagInv[0]*Apar[1]*Bmag[1]*geoY[2]*dfac_x)/dfac_v-(1.58113883008419*Apar[0]*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x)/dfac_v+(3.16227766016838*Apar[1]*geoY[2]*dfac_x)/dfac_v-(23.71708245126285*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v+(13.69306393762916*geoY[0]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v+(18.37117307087383*Apar[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v-(10.60660171779821*Apar[0]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v-(10.60660171779821*geoY[0]*Apar[1]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v+(6.123724356957945*Apar[0]*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v+(6.123724356957945*Bmag[1]*geoY[1]*Apar[2]*BmagInv[2]*dfac_x)/dfac_v-(3.535533905932738*geoY[0]*Bmag[1]*Apar[2]*BmagInv[2]*dfac_x)/dfac_v-(4.743416490252569*Apar[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x)/dfac_v+(2.738612787525831*Apar[0]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x)/dfac_v+(2.738612787525831*geoY[0]*Apar[1]*Bmag[1]*BmagInv[2]*dfac_x)/dfac_v-(1.58113883008419*Apar[0]*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x)/dfac_v+(18.37117307087383*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x)/dfac_v-(10.60660171779821*BmagInv[0]*geoY[1]*Apar[2]*Bmag[2]*dfac_x)/dfac_v-(10.60660171779821*geoY[0]*BmagInv[1]*Apar[2]*Bmag[2]*dfac_x)/dfac_v+(6.123724356957945*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*dfac_x)/dfac_v-(14.23024947075771*Apar[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x)/dfac_v+(8.21583836257749*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x)/dfac_v+(8.21583836257749*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*dfac_x)/dfac_v-(4.743416490252569*Apar[0]*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x)/dfac_v+(8.21583836257749*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*dfac_x)/dfac_v-(4.743416490252569*Apar[0]*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x)/dfac_v-(4.743416490252569*BmagInv[0]*geoY[0]*Apar[1]*Bmag[2]*dfac_x)/dfac_v+(2.738612787525831*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x)/dfac_v-(4.743416490252569*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*dfac_x)/dfac_v+(2.738612787525831*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*dfac_x)/dfac_v+(9.486832980505138*geoY[1]*Apar[2]*dfac_x)/dfac_v+(2.738612787525831*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*dfac_x)/dfac_v-(1.58113883008419*BmagInv[0]*geoY[0]*Bmag[1]*Apar[2]*dfac_x)/dfac_v-(5.477225575051662*geoY[0]*Apar[2]*dfac_x)/dfac_v+(3.674234614174767*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x)/dfac_v-(2.121320343559643*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x)/dfac_v-(2.121320343559643*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*dfac_x)/dfac_v+(1.224744871391589*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x)/dfac_v-(2.449489742783178*Apar[1]*geoY[1]*dfac_x)/dfac_v-(2.121320343559643*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*dfac_x)/dfac_v+(1.224744871391589*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x)/dfac_v+(1.224744871391589*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*dfac_x)/dfac_v-(0.7071067811865476*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x)/dfac_v+(1.414213562373095*geoY[0]*Apar[1]*dfac_x)/dfac_v+(2.581988897471611*Gradpar[2])/dfac_v-(2.0*Gradpar[1])/dfac_v+(1.154700538379251*Gradpar[0])/dfac_v; 
  alpha[8] = (10.0*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(2.581988897471612*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(7.745966692414835*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(4.47213595499958*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(2.0*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(1.154700538379251*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(7.745966692414835*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(4.47213595499958*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(2.0*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)-(1.154700538379251*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(6.0*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(3.464101615137754*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(3.464101615137754*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)+(2.0*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(1.549193338482967*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.8944271909999159*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.8944271909999159*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_)/(dfac_v2*q_)-(0.5163977794943223*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_)/(dfac_v2*q_); 
  if (alpha0>0) { 
  incr[0] = 0.01178511301977579*(25.98076211353316*alpha[8]*fl[12]+33.54101966249684*alpha[2]*fl[11]+15.0*alpha[8]*fl[8]+33.54101966249685*alpha[0]*fl[7]+25.98076211353316*alpha[2]*fl[4]+15.0*alpha[2]*fl[2]+25.98076211353316*alpha[0]*fl[1]+15.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.03535533905932736*(15.0*alpha[8]*fl[12]+19.36491673103708*alpha[2]*fl[11]+8.660254037844386*alpha[8]*fl[8]+19.36491673103709*alpha[0]*fl[7]+15.0*alpha[2]*fl[4]+8.660254037844386*alpha[2]*fl[2]+15.0*alpha[0]*fl[1]+8.660254037844386*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[12]+30.0*alpha[8]*fl[11]+33.54101966249684*alpha[0]*fl[11]+13.41640786499874*alpha[2]*fl[8]+23.2379000772445*fl[4]*alpha[8]+13.41640786499874*fl[2]*alpha[8]+33.54101966249685*alpha[2]*fl[7]+25.98076211353316*alpha[0]*fl[4]+15.0*alpha[0]*fl[2]+25.98076211353316*fl[1]*alpha[2]+15.0*fl[0]*alpha[2])*dfac_x; 
  incr[3] = 0.01178511301977579*(25.98076211353316*alpha[8]*fl[18]+33.54101966249685*alpha[2]*fl[17]+15.0*alpha[8]*fl[14]+33.54101966249684*alpha[0]*fl[13]+25.98076211353316*alpha[2]*fl[10]+15.0*alpha[2]*fl[6]+25.98076211353316*alpha[0]*fl[5]+15.0*alpha[0]*fl[3])*dfac_x; 
  incr[4] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[12]+17.32050807568877*alpha[8]*fl[11]+19.36491673103708*alpha[0]*fl[11]+7.745966692414834*alpha[2]*fl[8]+13.41640786499874*fl[4]*alpha[8]+7.745966692414834*fl[2]*alpha[8]+19.36491673103709*alpha[2]*fl[7]+15.0*alpha[0]*fl[4]+8.660254037844386*alpha[0]*fl[2]+15.0*fl[1]*alpha[2]+8.660254037844386*fl[0]*alpha[2])*dfac_x; 
  incr[5] = -0.03535533905932736*(15.0*alpha[8]*fl[18]+19.36491673103709*alpha[2]*fl[17]+8.660254037844387*alpha[8]*fl[14]+19.36491673103708*alpha[0]*fl[13]+15.0*alpha[2]*fl[10]+8.660254037844386*alpha[2]*fl[6]+15.0*alpha[0]*fl[5]+8.660254037844386*alpha[0]*fl[3])*dfac_x; 
  incr[6] = 0.01178511301977579*(23.2379000772445*alpha[2]*fl[18]+30.0*alpha[8]*fl[17]+33.54101966249685*alpha[0]*fl[17]+13.41640786499874*alpha[2]*fl[14]+33.54101966249684*alpha[2]*fl[13]+23.2379000772445*alpha[8]*fl[10]+25.98076211353316*alpha[0]*fl[10]+13.41640786499874*fl[6]*alpha[8]+15.0*alpha[0]*fl[6]+25.98076211353316*alpha[2]*fl[5]+15.0*alpha[2]*fl[3])*dfac_x; 
  incr[7] = 0.05892556509887893*(11.61895003862225*alpha[8]*fl[12]+15.0*alpha[2]*fl[11]+6.708203932499369*alpha[8]*fl[8]+15.0*alpha[0]*fl[7]+11.61895003862225*alpha[2]*fl[4]+6.708203932499369*alpha[2]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[8] = 0.001683587574253684*(116.1895003862225*alpha[8]*fl[12]+181.8653347947321*alpha[0]*fl[12]+210.0*alpha[2]*fl[11]+67.0820393249937*alpha[8]*fl[8]+105.0*alpha[0]*fl[8]+234.787137637478*fl[7]*alpha[8]+181.8653347947321*fl[1]*alpha[8]+105.0*fl[0]*alpha[8]+162.6653005407115*alpha[2]*fl[4]+93.91485505499116*alpha[2]*fl[2])*dfac_x; 
  incr[9] = 0.01178511301977579*(25.98076211353316*alpha[2]*fl[19]+15.0*alpha[2]*fl[16]+25.98076211353316*alpha[0]*fl[15]+15.0*alpha[0]*fl[9])*dfac_x; 
  incr[10] = -0.03535533905932736*(13.41640786499874*alpha[2]*fl[18]+17.32050807568877*alpha[8]*fl[17]+19.36491673103709*alpha[0]*fl[17]+7.745966692414834*alpha[2]*fl[14]+19.36491673103708*alpha[2]*fl[13]+13.41640786499874*alpha[8]*fl[10]+15.0*alpha[0]*fl[10]+7.745966692414834*fl[6]*alpha[8]+8.660254037844386*alpha[0]*fl[6]+15.0*alpha[2]*fl[5]+8.660254037844386*alpha[2]*fl[3])*dfac_x; 
  incr[11] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[12]+67.0820393249937*alpha[8]*fl[11]+75.0*alpha[0]*fl[11]+30.0*alpha[2]*fl[8]+51.96152422706632*fl[4]*alpha[8]+30.0*fl[2]*alpha[8]+75.00000000000001*alpha[2]*fl[7]+58.09475019311126*alpha[0]*fl[4]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*fl[1]*alpha[2]+33.54101966249684*fl[0]*alpha[2])*dfac_x; 
  incr[12] = -0.005050762722761052*(67.0820393249937*alpha[8]*fl[12]+105.0*alpha[0]*fl[12]+121.2435565298214*alpha[2]*fl[11]+38.72983346207417*alpha[8]*fl[8]+60.62177826491071*alpha[0]*fl[8]+135.5544171172596*fl[7]*alpha[8]+105.0*fl[1]*alpha[8]+60.62177826491071*fl[0]*alpha[8]+93.91485505499116*alpha[2]*fl[4]+54.22176684690384*alpha[2]*fl[2])*dfac_x; 
  incr[13] = 0.05892556509887893*(11.61895003862225*alpha[8]*fl[18]+15.0*alpha[2]*fl[17]+6.708203932499369*alpha[8]*fl[14]+15.0*alpha[0]*fl[13]+11.61895003862225*alpha[2]*fl[10]+6.708203932499369*alpha[2]*fl[6]+11.61895003862225*alpha[0]*fl[5]+6.708203932499369*alpha[0]*fl[3])*dfac_x; 
  incr[14] = 0.001683587574253684*(116.1895003862225*alpha[8]*fl[18]+181.8653347947321*alpha[0]*fl[18]+210.0*alpha[2]*fl[17]+67.0820393249937*alpha[8]*fl[14]+105.0*alpha[0]*fl[14]+234.787137637478*alpha[8]*fl[13]+162.6653005407115*alpha[2]*fl[10]+181.8653347947321*fl[5]*alpha[8]+105.0*fl[3]*alpha[8]+93.91485505499116*alpha[2]*fl[6])*dfac_x; 
  incr[15] = -0.03535533905932736*(15.0*alpha[2]*fl[19]+8.660254037844386*alpha[2]*fl[16]+15.0*alpha[0]*fl[15]+8.660254037844387*alpha[0]*fl[9])*dfac_x; 
  incr[16] = 0.01178511301977579*(23.2379000772445*alpha[8]*fl[19]+25.98076211353316*alpha[0]*fl[19]+13.41640786499874*alpha[8]*fl[16]+15.0*alpha[0]*fl[16]+25.98076211353316*alpha[2]*fl[15]+15.0*alpha[2]*fl[9])*dfac_x; 
  incr[17] = 0.01178511301977579*(51.96152422706631*alpha[2]*fl[18]+67.0820393249937*alpha[8]*fl[17]+75.0*alpha[0]*fl[17]+30.0*alpha[2]*fl[14]+75.00000000000001*alpha[2]*fl[13]+51.96152422706631*alpha[8]*fl[10]+58.09475019311126*alpha[0]*fl[10]+30.0*fl[6]*alpha[8]+33.54101966249685*alpha[0]*fl[6]+58.09475019311126*alpha[2]*fl[5]+33.54101966249685*alpha[2]*fl[3])*dfac_x; 
  incr[18] = -0.005050762722761052*(67.0820393249937*alpha[8]*fl[18]+105.0*alpha[0]*fl[18]+121.2435565298214*alpha[2]*fl[17]+38.72983346207417*alpha[8]*fl[14]+60.62177826491071*alpha[0]*fl[14]+135.5544171172596*alpha[8]*fl[13]+93.91485505499116*alpha[2]*fl[10]+105.0*fl[5]*alpha[8]+60.6217782649107*fl[3]*alpha[8]+54.22176684690384*alpha[2]*fl[6])*dfac_x; 
  incr[19] = -0.03535533905932736*(13.41640786499874*alpha[8]*fl[19]+15.0*alpha[0]*fl[19]+7.745966692414834*alpha[8]*fl[16]+8.660254037844387*alpha[0]*fl[16]+15.0*alpha[2]*fl[15]+8.660254037844386*alpha[2]*fl[9])*dfac_x; 

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
  incr[0] = -0.01178511301977579*(25.98076211353316*alpha[8]*fr[12]-33.54101966249684*alpha[2]*fr[11]-15.0*alpha[8]*fr[8]-33.54101966249685*alpha[0]*fr[7]+25.98076211353316*alpha[2]*fr[4]-15.0*alpha[2]*fr[2]+25.98076211353316*alpha[0]*fr[1]-15.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.03535533905932736*(15.0*alpha[8]*fr[12]-19.36491673103708*alpha[2]*fr[11]-8.660254037844386*alpha[8]*fr[8]-19.36491673103709*alpha[0]*fr[7]+15.0*alpha[2]*fr[4]-8.660254037844386*alpha[2]*fr[2]+15.0*alpha[0]*fr[1]-8.660254037844386*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[12]-30.0*alpha[8]*fr[11]-33.54101966249684*alpha[0]*fr[11]-13.41640786499874*alpha[2]*fr[8]+23.2379000772445*fr[4]*alpha[8]-13.41640786499874*fr[2]*alpha[8]-33.54101966249685*alpha[2]*fr[7]+25.98076211353316*alpha[0]*fr[4]-15.0*alpha[0]*fr[2]+25.98076211353316*fr[1]*alpha[2]-15.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = -0.01178511301977579*(25.98076211353316*alpha[8]*fr[18]-33.54101966249685*alpha[2]*fr[17]-15.0*alpha[8]*fr[14]-33.54101966249684*alpha[0]*fr[13]+25.98076211353316*alpha[2]*fr[10]-15.0*alpha[2]*fr[6]+25.98076211353316*alpha[0]*fr[5]-15.0*alpha[0]*fr[3])*dfac_x; 
  incr[4] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[12]-17.32050807568877*alpha[8]*fr[11]-19.36491673103708*alpha[0]*fr[11]-7.745966692414834*alpha[2]*fr[8]+13.41640786499874*fr[4]*alpha[8]-7.745966692414834*fr[2]*alpha[8]-19.36491673103709*alpha[2]*fr[7]+15.0*alpha[0]*fr[4]-8.660254037844386*alpha[0]*fr[2]+15.0*fr[1]*alpha[2]-8.660254037844386*fr[0]*alpha[2])*dfac_x; 
  incr[5] = 0.03535533905932736*(15.0*alpha[8]*fr[18]-19.36491673103709*alpha[2]*fr[17]-8.660254037844387*alpha[8]*fr[14]-19.36491673103708*alpha[0]*fr[13]+15.0*alpha[2]*fr[10]-8.660254037844386*alpha[2]*fr[6]+15.0*alpha[0]*fr[5]-8.660254037844386*alpha[0]*fr[3])*dfac_x; 
  incr[6] = -0.01178511301977579*(23.2379000772445*alpha[2]*fr[18]-30.0*alpha[8]*fr[17]-33.54101966249685*alpha[0]*fr[17]-13.41640786499874*alpha[2]*fr[14]-33.54101966249684*alpha[2]*fr[13]+23.2379000772445*alpha[8]*fr[10]+25.98076211353316*alpha[0]*fr[10]-13.41640786499874*fr[6]*alpha[8]-15.0*alpha[0]*fr[6]+25.98076211353316*alpha[2]*fr[5]-15.0*alpha[2]*fr[3])*dfac_x; 
  incr[7] = -0.05892556509887893*(11.61895003862225*alpha[8]*fr[12]-15.0*alpha[2]*fr[11]-6.708203932499369*alpha[8]*fr[8]-15.0*alpha[0]*fr[7]+11.61895003862225*alpha[2]*fr[4]-6.708203932499369*alpha[2]*fr[2]+11.61895003862225*alpha[0]*fr[1]-6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[8] = -0.001683587574253684*(116.1895003862225*alpha[8]*fr[12]+181.8653347947321*alpha[0]*fr[12]-210.0*alpha[2]*fr[11]-67.0820393249937*alpha[8]*fr[8]-105.0*alpha[0]*fr[8]-234.787137637478*fr[7]*alpha[8]+181.8653347947321*fr[1]*alpha[8]-105.0*fr[0]*alpha[8]+162.6653005407115*alpha[2]*fr[4]-93.91485505499116*alpha[2]*fr[2])*dfac_x; 
  incr[9] = -0.01178511301977579*(25.98076211353316*alpha[2]*fr[19]-15.0*alpha[2]*fr[16]+25.98076211353316*alpha[0]*fr[15]-15.0*alpha[0]*fr[9])*dfac_x; 
  incr[10] = 0.03535533905932736*(13.41640786499874*alpha[2]*fr[18]-17.32050807568877*alpha[8]*fr[17]-19.36491673103709*alpha[0]*fr[17]-7.745966692414834*alpha[2]*fr[14]-19.36491673103708*alpha[2]*fr[13]+13.41640786499874*alpha[8]*fr[10]+15.0*alpha[0]*fr[10]-7.745966692414834*fr[6]*alpha[8]-8.660254037844386*alpha[0]*fr[6]+15.0*alpha[2]*fr[5]-8.660254037844386*alpha[2]*fr[3])*dfac_x; 
  incr[11] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[12]-67.0820393249937*alpha[8]*fr[11]-75.0*alpha[0]*fr[11]-30.0*alpha[2]*fr[8]+51.96152422706632*fr[4]*alpha[8]-30.0*fr[2]*alpha[8]-75.00000000000001*alpha[2]*fr[7]+58.09475019311126*alpha[0]*fr[4]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*fr[1]*alpha[2]-33.54101966249684*fr[0]*alpha[2])*dfac_x; 
  incr[12] = 0.005050762722761052*(67.0820393249937*alpha[8]*fr[12]+105.0*alpha[0]*fr[12]-121.2435565298214*alpha[2]*fr[11]-38.72983346207417*alpha[8]*fr[8]-60.62177826491071*alpha[0]*fr[8]-135.5544171172596*fr[7]*alpha[8]+105.0*fr[1]*alpha[8]-60.62177826491071*fr[0]*alpha[8]+93.91485505499116*alpha[2]*fr[4]-54.22176684690384*alpha[2]*fr[2])*dfac_x; 
  incr[13] = -0.05892556509887893*(11.61895003862225*alpha[8]*fr[18]-15.0*alpha[2]*fr[17]-6.708203932499369*alpha[8]*fr[14]-15.0*alpha[0]*fr[13]+11.61895003862225*alpha[2]*fr[10]-6.708203932499369*alpha[2]*fr[6]+11.61895003862225*alpha[0]*fr[5]-6.708203932499369*alpha[0]*fr[3])*dfac_x; 
  incr[14] = -0.001683587574253684*(116.1895003862225*alpha[8]*fr[18]+181.8653347947321*alpha[0]*fr[18]-210.0*alpha[2]*fr[17]-67.0820393249937*alpha[8]*fr[14]-105.0*alpha[0]*fr[14]-234.787137637478*alpha[8]*fr[13]+162.6653005407115*alpha[2]*fr[10]+181.8653347947321*fr[5]*alpha[8]-105.0*fr[3]*alpha[8]-93.91485505499116*alpha[2]*fr[6])*dfac_x; 
  incr[15] = 0.03535533905932736*(15.0*alpha[2]*fr[19]-8.660254037844386*alpha[2]*fr[16]+15.0*alpha[0]*fr[15]-8.660254037844387*alpha[0]*fr[9])*dfac_x; 
  incr[16] = -0.01178511301977579*(23.2379000772445*alpha[8]*fr[19]+25.98076211353316*alpha[0]*fr[19]-13.41640786499874*alpha[8]*fr[16]-15.0*alpha[0]*fr[16]+25.98076211353316*alpha[2]*fr[15]-15.0*alpha[2]*fr[9])*dfac_x; 
  incr[17] = -0.01178511301977579*(51.96152422706631*alpha[2]*fr[18]-67.0820393249937*alpha[8]*fr[17]-75.0*alpha[0]*fr[17]-30.0*alpha[2]*fr[14]-75.00000000000001*alpha[2]*fr[13]+51.96152422706631*alpha[8]*fr[10]+58.09475019311126*alpha[0]*fr[10]-30.0*fr[6]*alpha[8]-33.54101966249685*alpha[0]*fr[6]+58.09475019311126*alpha[2]*fr[5]-33.54101966249685*alpha[2]*fr[3])*dfac_x; 
  incr[18] = 0.005050762722761052*(67.0820393249937*alpha[8]*fr[18]+105.0*alpha[0]*fr[18]-121.2435565298214*alpha[2]*fr[17]-38.72983346207417*alpha[8]*fr[14]-60.62177826491071*alpha[0]*fr[14]-135.5544171172596*alpha[8]*fr[13]+93.91485505499116*alpha[2]*fr[10]+105.0*fr[5]*alpha[8]-60.6217782649107*fr[3]*alpha[8]-54.22176684690384*alpha[2]*fr[6])*dfac_x; 
  incr[19] = 0.03535533905932736*(13.41640786499874*alpha[8]*fr[19]+15.0*alpha[0]*fr[19]-7.745966692414834*alpha[8]*fr[16]-8.660254037844387*alpha[0]*fr[16]+15.0*alpha[2]*fr[15]-8.660254037844386*alpha[2]*fr[9])*dfac_x; 

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
double EmGyrokineticSurf1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.1767766952966368*(((6.0*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]+6.0*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1])*dfac_v*dfac_x*dfac_x*m_*wm+(6.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]+6.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_v*dfac_x*dfac_x*m_*q_)*wv+((((4.242640687119286*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[1]+4.242640687119286*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[1])*geoY[1]+4.242640687119286*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]+4.242640687119286*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]-8.485281374238571*geoY[0]*Apar[1]*Bmag[1])*dfac_v*dfac_x*dfac_x-6.928203230275509*Gradpar[0]*Bmag[1]*dfac_v*dfac_x)*q_+((-6.0*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1])-6.0*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1])*dfac_x*dfac_x*m_)*wm+((-6.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1])-6.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_x*dfac_x*m_*q_+(((4.242640687119286*Apar[0]*Bmag[1]*BmagInv[1]+4.242640687119286*BmagInv[0]*Apar[1]*Bmag[1])*Phi[1]*geoY[1]+(4.242640687119286*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]+4.242640687119286*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]-8.485281374238571*geoY[0]*Apar[1])*Phi[1])*dfac_v*dfac_x*dfac_x-6.928203230275509*Gradpar[0]*Phi[1]*dfac_v*dfac_x-5.656854249492382*dApardt[0]*dfac_v)*q2))/(dfac_v*m_*q_); 

  double alpha[8]; 
  alpha[0] = (2.121320343559642*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm*wv)/q_+2.121320343559642*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*wv-(2.121320343559642*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)+(1.5*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/m_+(1.5*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/m_+(1.5*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm)/m_+(1.5*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm)/m_-(3.0*geoY[0]*Apar[1]*Bmag[1]*dfac_x*dfac_x*wm)/m_-(2.449489742783177*Gradpar[0]*Bmag[1]*dfac_x*wm)/m_+(1.5*Apar[0]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.5*BmagInv[0]*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.5*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_+(1.5*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(3.0*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(2.449489742783177*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[0]*q_)/m_-(2.121320343559642*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[1] = (2.121320343559642*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm*wv)/q_+2.121320343559642*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+2.121320343559642*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*wv-(2.121320343559642*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)+(2.7*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/m_+(1.5*Apar[0]*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/m_-(3.0*Apar[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/m_+(1.5*Apar[0]*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm)/m_+(1.5*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm)/m_-(2.449489742783177*Bmag[1]*Gradpar[1]*dfac_x*wm)/m_+(2.7*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.5*Apar[0]*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_-(3.0*Apar[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.5*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_+(1.5*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(2.449489742783177*Gradpar[1]*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[1]*q_)/m_-(2.121320343559642*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[3] = (1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)-(1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)+(0.8660254037844385*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x)/(dfac_m*m_)-(1.732050807568877*geoY[0]*Apar[1]*Bmag[1]*dfac_x*dfac_x)/(dfac_m*m_)-(1.414213562373095*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[5] = (1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)-(1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)+(1.55884572681199*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*Apar[0]*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)-(1.732050807568877*Apar[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*Apar[0]*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*dfac_x*dfac_x)/(dfac_m*m_)-(1.414213562373095*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[5]*fl[7]+1.732050807568877*alpha[3]*fl[6]+alpha[5]*fl[5]+1.732050807568877*alpha[1]*fl[4]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*alpha[3]*fl[7]+1.732050807568877*alpha[5]*fl[6]+alpha[3]*fl[5]+fl[3]*alpha[5]+1.732050807568877*alpha[0]*fl[4]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.1767766952966368*(3.0*alpha[5]*fl[7]+3.0*alpha[3]*fl[6]+1.732050807568877*alpha[5]*fl[5]+3.0*alpha[1]*fl[4]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[1]*fl[7]+1.732050807568877*alpha[0]*fl[6]+alpha[1]*fl[5]+1.732050807568877*fl[4]*alpha[5]+fl[1]*alpha[5]+alpha[0]*fl[3]+1.732050807568877*fl[2]*alpha[3]+fl[0]*alpha[3])*dfac_v; 
  incr[4] = -0.1767766952966368*(3.0*alpha[3]*fl[7]+3.0*alpha[5]*fl[6]+1.732050807568877*alpha[3]*fl[5]+1.732050807568877*fl[3]*alpha[5]+3.0*alpha[0]*fl[4]+3.0*alpha[1]*fl[2]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[6]+alpha[0]*fl[5]+1.732050807568877*fl[2]*alpha[5]+fl[0]*alpha[5]+1.732050807568877*alpha[3]*fl[4]+alpha[1]*fl[3]+fl[1]*alpha[3])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*alpha[1]*fl[7]+3.0*alpha[0]*fl[6]+1.732050807568877*alpha[1]*fl[5]+3.0*fl[4]*alpha[5]+1.732050807568877*fl[1]*alpha[5]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fl[7]+3.0*alpha[1]*fl[6]+1.732050807568877*alpha[0]*fl[5]+3.0*fl[2]*alpha[5]+1.732050807568877*fl[0]*alpha[5]+3.0*alpha[3]*fl[4]+1.732050807568877*alpha[1]*fl[3]+1.732050807568877*fl[1]*alpha[3])*dfac_v; 

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
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[5]*fr[7]+1.732050807568877*alpha[3]*fr[6]-1.0*alpha[5]*fr[5]+1.732050807568877*alpha[1]*fr[4]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*alpha[3]*fr[7]+1.732050807568877*alpha[5]*fr[6]-1.0*alpha[3]*fr[5]-1.0*fr[3]*alpha[5]+1.732050807568877*alpha[0]*fr[4]+1.732050807568877*alpha[1]*fr[2]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(3.0*alpha[5]*fr[7]+3.0*alpha[3]*fr[6]-1.732050807568877*alpha[5]*fr[5]+3.0*alpha[1]*fr[4]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[1]*fr[7]+1.732050807568877*alpha[0]*fr[6]-1.0*alpha[1]*fr[5]+1.732050807568877*fr[4]*alpha[5]-1.0*fr[1]*alpha[5]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_v; 
  incr[4] = 0.1767766952966368*(3.0*alpha[3]*fr[7]+3.0*alpha[5]*fr[6]-1.732050807568877*alpha[3]*fr[5]-1.732050807568877*fr[3]*alpha[5]+3.0*alpha[0]*fr[4]+3.0*alpha[1]*fr[2]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[7]+1.732050807568877*alpha[1]*fr[6]-1.0*alpha[0]*fr[5]+1.732050807568877*fr[2]*alpha[5]-1.0*fr[0]*alpha[5]+1.732050807568877*alpha[3]*fr[4]-1.0*alpha[1]*fr[3]-1.0*fr[1]*alpha[3])*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*alpha[1]*fr[7]+3.0*alpha[0]*fr[6]-1.732050807568877*alpha[1]*fr[5]+3.0*fr[4]*alpha[5]-1.732050807568877*fr[1]*alpha[5]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fr[7]+3.0*alpha[1]*fr[6]-1.732050807568877*alpha[0]*fr[5]+3.0*fr[2]*alpha[5]-1.732050807568877*fr[0]*alpha[5]+3.0*alpha[3]*fr[4]-1.732050807568877*alpha[1]*fr[3]-1.732050807568877*fr[1]*alpha[3])*dfac_v; 

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
double EmGyrokineticSurf1x2vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.005050762722761052*(((((1650.0*Bmag[2]*Bmag[2]+210.0*Bmag[1]*Bmag[1])*BmagInv[2]+939.1485505499119*BmagInv[0]*Bmag[2]*Bmag[2]+840.0*Bmag[1]*BmagInv[1]*Bmag[2])*geoY[2]+(939.1485505499119*geoY[0]*Bmag[2]*Bmag[2]+840.0*Bmag[1]*geoY[1]*Bmag[2])*BmagInv[2]+(1890.0*BmagInv[1]*geoY[1]+1050.0*BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+(939.1485505499119*BmagInv[0]*Bmag[1]*geoY[1]+939.1485505499119*geoY[0]*Bmag[1]*BmagInv[1])*Bmag[2]+210.0*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]+210.0*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1])*dfac_v*dfac_x*dfac_x*m_*wm+(((1650.0*Bmag[2]*BmagInv[2]+939.1485505499119*BmagInv[0]*Bmag[2]+420.0*Bmag[1]*BmagInv[1])*Phi[2]+210.0*Bmag[1]*Phi[1]*BmagInv[2]+420.0*BmagInv[1]*Phi[1]*Bmag[2])*geoY[2]+((939.1485505499119*geoY[0]*Bmag[2]+420.0*Bmag[1]*geoY[1])*BmagInv[2]+(1890.0*BmagInv[1]*geoY[1]+1050.0*BmagInv[0]*geoY[0])*Bmag[2]+469.5742752749559*BmagInv[0]*Bmag[1]*geoY[1]+469.5742752749559*geoY[0]*Bmag[1]*BmagInv[1])*Phi[2]+420.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]+(469.5742752749559*BmagInv[0]*Phi[1]*geoY[1]+469.5742752749559*geoY[0]*BmagInv[1]*Phi[1])*Bmag[2]+210.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]+210.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_v*dfac_x*dfac_x*m_*q_)*wv+((((((1897.366596101029*Apar[2]+1166.726188957804*Apar[0])*Bmag[2]*Bmag[2]+1043.551627855566*Apar[1]*Bmag[1]*Bmag[2]+94.86832980505142*Bmag[1]*Bmag[1]*Apar[2]+148.492424049175*Apar[0]*Bmag[1]*Bmag[1])*BmagInv[2]+(1166.726188957804*BmagInv[0]*Apar[2]+1707.629936490926*Apar[1]*BmagInv[1]+664.0783086353599*Apar[0]*BmagInv[0])*Bmag[2]*Bmag[2]+((1043.551627855566*Bmag[1]*BmagInv[1]-1328.15661727072)*Apar[2]+593.9696961967002*Apar[0]*Bmag[1]*BmagInv[1]+593.9696961967002*BmagInv[0]*Apar[1]*Bmag[1])*Bmag[2]+148.492424049175*BmagInv[0]*Bmag[1]*Bmag[1]*Apar[2]+132.815661727072*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1])*geoY[2]+((1166.726188957804*geoY[0]*Apar[2]+1707.629936490926*Apar[1]*geoY[1]+664.0783086353599*Apar[0]*geoY[0])*Bmag[2]*Bmag[2]+(1043.551627855566*Bmag[1]*geoY[1]*Apar[2]+593.9696961967002*Apar[0]*Bmag[1]*geoY[1]+593.9696961967002*geoY[0]*Apar[1]*Bmag[1])*Bmag[2]+148.492424049175*geoY[0]*Bmag[1]*Bmag[1]*Apar[2]+132.815661727072*Apar[1]*Bmag[1]*Bmag[1]*geoY[1])*BmagInv[2]+((1707.629936490926*BmagInv[1]*geoY[1]+664.0783086353599*BmagInv[0]*geoY[0])*Apar[2]+(1336.431816442575*Apar[0]*BmagInv[1]+1336.431816442575*BmagInv[0]*Apar[1])*geoY[1]+1336.431816442575*geoY[0]*Apar[1]*BmagInv[1]+742.462120245875*Apar[0]*BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+((593.9696961967002*BmagInv[0]*Bmag[1]*geoY[1]+593.9696961967002*geoY[0]*Bmag[1]*BmagInv[1]-1484.92424049175*geoY[0])*Apar[2]+(1195.340955543648*Apar[1]*Bmag[1]*BmagInv[1]+664.0783086353599*Apar[0]*BmagInv[0]*Bmag[1]-664.0783086353599*Apar[1])*geoY[1]+664.0783086353599*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]+664.0783086353599*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1])*Bmag[2]+(132.815661727072*Bmag[1]*Bmag[1]*BmagInv[1]-664.0783086353599*Bmag[1])*geoY[1]*Apar[2]+(148.492424049175*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[1]+148.492424049175*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[1])*geoY[1]+148.492424049175*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]+148.492424049175*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]-296.98484809835*geoY[0]*Apar[1]*Bmag[1])*dfac_v*dfac_x*dfac_x+((-542.2176684690385*Gradpar[1]*Bmag[2])-242.4871130596428*Gradpar[0]*Bmag[1])*dfac_v*dfac_x)*q_+((((-1650.0*Bmag[2]*Bmag[2])-210.0*Bmag[1]*Bmag[1])*BmagInv[2]-939.1485505499119*BmagInv[0]*Bmag[2]*Bmag[2]-840.0*Bmag[1]*BmagInv[1]*Bmag[2])*geoY[2]+((-939.1485505499119*geoY[0]*Bmag[2]*Bmag[2])-840.0*Bmag[1]*geoY[1]*Bmag[2])*BmagInv[2]+((-1890.0*BmagInv[1]*geoY[1])-1050.0*BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+((-939.1485505499119*BmagInv[0]*Bmag[1]*geoY[1])-939.1485505499119*geoY[0]*Bmag[1]*BmagInv[1])*Bmag[2]-210.0*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]-210.0*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1])*dfac_x*dfac_x*m_)*wm+((((-1650.0*Bmag[2]*BmagInv[2])-939.1485505499119*BmagInv[0]*Bmag[2]-420.0*Bmag[1]*BmagInv[1])*Phi[2]-210.0*Bmag[1]*Phi[1]*BmagInv[2]-420.0*BmagInv[1]*Phi[1]*Bmag[2])*geoY[2]+(((-939.1485505499119*geoY[0]*Bmag[2])-420.0*Bmag[1]*geoY[1])*BmagInv[2]+((-1890.0*BmagInv[1]*geoY[1])-1050.0*BmagInv[0]*geoY[0])*Bmag[2]-469.5742752749559*BmagInv[0]*Bmag[1]*geoY[1]-469.5742752749559*geoY[0]*Bmag[1]*BmagInv[1])*Phi[2]-420.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]+((-469.5742752749559*BmagInv[0]*Phi[1]*geoY[1])-469.5742752749559*geoY[0]*BmagInv[1]*Phi[1])*Bmag[2]-210.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]-210.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_x*dfac_x*m_*q_+((((((1897.366596101029*Apar[2]+1166.726188957804*Apar[0])*Bmag[2]+521.7758139277827*Apar[1]*Bmag[1])*BmagInv[2]+(1166.726188957804*BmagInv[0]*Apar[2]+1707.629936490926*Apar[1]*BmagInv[1]+664.0783086353599*Apar[0]*BmagInv[0])*Bmag[2]+(521.7758139277827*Bmag[1]*BmagInv[1]-1328.15661727072)*Apar[2]+296.98484809835*Apar[0]*Bmag[1]*BmagInv[1]+296.98484809835*BmagInv[0]*Apar[1]*Bmag[1])*Phi[2]+(521.7758139277827*Apar[1]*Phi[1]*Bmag[2]+94.86832980505142*Bmag[1]*Phi[1]*Apar[2]+148.492424049175*Apar[0]*Bmag[1]*Phi[1])*BmagInv[2]+(521.7758139277827*BmagInv[1]*Phi[1]*Apar[2]+(296.98484809835*Apar[0]*BmagInv[1]+296.98484809835*BmagInv[0]*Apar[1])*Phi[1])*Bmag[2]+148.492424049175*BmagInv[0]*Bmag[1]*Phi[1]*Apar[2]+132.815661727072*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1])*geoY[2]+(((1166.726188957804*geoY[0]*Apar[2]+1707.629936490926*Apar[1]*geoY[1]+664.0783086353599*Apar[0]*geoY[0])*Bmag[2]+521.7758139277827*Bmag[1]*geoY[1]*Apar[2]+296.98484809835*Apar[0]*Bmag[1]*geoY[1]+296.98484809835*geoY[0]*Apar[1]*Bmag[1])*BmagInv[2]+((1707.629936490926*BmagInv[1]*geoY[1]+664.0783086353599*BmagInv[0]*geoY[0])*Apar[2]+(1336.431816442575*Apar[0]*BmagInv[1]+1336.431816442575*BmagInv[0]*Apar[1])*geoY[1]+1336.431816442575*geoY[0]*Apar[1]*BmagInv[1]+742.462120245875*Apar[0]*BmagInv[0]*geoY[0])*Bmag[2]+(296.98484809835*BmagInv[0]*Bmag[1]*geoY[1]+296.98484809835*geoY[0]*Bmag[1]*BmagInv[1]-1484.92424049175*geoY[0])*Apar[2]+(597.6704777718237*Apar[1]*Bmag[1]*BmagInv[1]+332.0391543176799*Apar[0]*BmagInv[0]*Bmag[1]-664.0783086353599*Apar[1])*geoY[1]+332.0391543176799*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]+332.0391543176799*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1])*Phi[2]+((521.7758139277827*Phi[1]*geoY[1]*Apar[2]+296.98484809835*Apar[0]*Phi[1]*geoY[1]+296.98484809835*geoY[0]*Apar[1]*Phi[1])*Bmag[2]+148.492424049175*geoY[0]*Bmag[1]*Phi[1]*Apar[2]+132.815661727072*Apar[1]*Bmag[1]*Phi[1]*geoY[1])*BmagInv[2]+((296.98484809835*BmagInv[0]*Phi[1]*geoY[1]+296.98484809835*geoY[0]*BmagInv[1]*Phi[1])*Apar[2]+(597.6704777718237*Apar[1]*BmagInv[1]+332.0391543176799*Apar[0]*BmagInv[0])*Phi[1]*geoY[1]+(332.0391543176799*Apar[0]*geoY[0]*BmagInv[1]+332.0391543176799*BmagInv[0]*geoY[0]*Apar[1])*Phi[1])*Bmag[2]+(132.815661727072*Bmag[1]*BmagInv[1]-664.0783086353599)*Phi[1]*geoY[1]*Apar[2]+(148.492424049175*Apar[0]*Bmag[1]*BmagInv[1]+148.492424049175*BmagInv[0]*Apar[1]*Bmag[1])*Phi[1]*geoY[1]+(148.492424049175*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]+148.492424049175*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]-296.98484809835*geoY[0]*Apar[1])*Phi[1])*dfac_v*dfac_x*dfac_x+((-542.2176684690385*Gradpar[1]*Phi[2])-242.4871130596428*Gradpar[0]*Phi[1])*dfac_v*dfac_x-197.9898987322334*dApardt[0]*dfac_v)*q2))/(dfac_v*m_*q_); 

  double alpha[20]; 
  alpha[0] = (16.66751698511147*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(9.48683298050514*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(8.485281374238568*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(9.48683298050514*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(8.485281374238568*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(19.09188309203678*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(10.60660171779821*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(9.48683298050514*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(9.48683298050514*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm*wv)/q_+16.66751698511147*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+9.48683298050514*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+4.242640687119284*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+2.121320343559642*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+4.242640687119284*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+9.48683298050514*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119284*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+19.09188309203678*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+10.60660171779821*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+4.74341649025257*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+4.74341649025257*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119284*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+4.74341649025257*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+4.74341649025257*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+2.121320343559642*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*wv-(16.66751698511147*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505136*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.485281374238568*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505136*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.485281374238568*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(19.09188309203678*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(10.60660171779821*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505136*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505136*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)+(19.16629694999821*Apar[2]*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(11.78571428571428*Apar[0]*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*Apar[1]*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(0.95831484749991*Bmag[1]*Bmag[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(1.5*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(11.78571428571428*BmagInv[0]*Apar[2]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499839*Apar[1]*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(6.70820393249937*Apar[0]*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*Bmag[1]*BmagInv[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_-(13.41640786499874*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*Apar[0]*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(1.5*BmagInv[0]*Bmag[1]*Bmag[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x*wm)/m_+(11.78571428571428*geoY[0]*Apar[2]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499839*Apar[1]*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(6.70820393249937*Apar[0]*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*Bmag[1]*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*Apar[0]*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*geoY[0]*Apar[1]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(1.5*geoY[0]*Bmag[1]*Bmag[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*Apar[1]*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499839*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.70820393249937*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(13.5*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(13.5*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(13.5*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(7.499999999999999*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_-(15.0*geoY[0]*Apar[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(12.07476707849886*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.70820393249937*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_-(6.70820393249937*Apar[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.70820393249937*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.70820393249937*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*wm)/m_-(6.70820393249937*Bmag[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*wm)/m_+(1.5*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/m_+(1.5*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/m_+(1.5*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm)/m_+(1.5*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm)/m_-(3.0*geoY[0]*Apar[1]*Bmag[1]*dfac_x*dfac_x*wm)/m_-(5.47722557505166*Gradpar[1]*Bmag[2]*dfac_x*wm)/m_-(2.449489742783177*Gradpar[0]*Bmag[1]*dfac_x*wm)/m_+(19.16629694999821*Apar[2]*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(11.78571428571428*Apar[0]*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Apar[1]*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(11.78571428571428*BmagInv[0]*Apar[2]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499839*Apar[1]*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(6.70820393249937*Apar[0]*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Bmag[1]*BmagInv[1]*Apar[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_-(13.41640786499874*Apar[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.0*BmagInv[0]*Apar[1]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Apar[1]*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(0.95831484749991*Bmag[1]*Phi[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.5*Apar[0]*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*BmagInv[1]*Phi[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.0*BmagInv[0]*Apar[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.5*BmagInv[0]*Bmag[1]*Phi[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_+(11.78571428571428*geoY[0]*Apar[2]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499839*Apar[1]*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(6.70820393249937*Apar[0]*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Bmag[1]*geoY[1]*Apar[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.0*geoY[0]*Apar[1]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499839*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(6.70820393249937*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(13.5*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(13.5*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(13.5*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(7.499999999999999*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.0*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.0*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_-(15.0*geoY[0]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(6.037383539249431*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.354101966249685*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_-(6.70820393249937*Apar[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.354101966249685*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.354101966249685*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Phi[1]*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(3.0*geoY[0]*Apar[1]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(1.5*geoY[0]*Bmag[1]*Phi[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(3.0*BmagInv[0]*Phi[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(3.0*geoY[0]*BmagInv[1]*Phi[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(6.037383539249431*Apar[1]*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(3.354101966249685*Apar[0]*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(3.354101966249685*Apar[0]*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(3.354101966249685*BmagInv[0]*geoY[0]*Apar[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*q_)/m_-(6.70820393249937*Phi[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*q_)/m_+(1.5*Apar[0]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.5*BmagInv[0]*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.5*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_+(1.5*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(3.0*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(5.47722557505166*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(2.449489742783177*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[0]*q_)/m_-(16.66751698511147*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(9.486832980505136*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(9.486832980505136*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(19.09188309203678*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(10.60660171779821*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.743416490252568*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.743416490252568*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(4.743416490252568*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(4.743416490252568*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[1] = (14.90788039793665*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(24.39471337844178*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(8.485281374238568*BmagInv[0]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(1.897366596101027*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(24.39471337844178*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(8.485281374238568*geoY[0]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(1.897366596101027*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(19.09188309203678*BmagInv[0]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(19.09188309203678*geoY[0]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(17.07629936490925*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(9.48683298050514*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm*wv)/q_+7.453940198968323*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+24.39471337844178*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+4.242640687119284*BmagInv[0]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+7.453940198968323*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+4.242640687119284*BmagInv[0]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+1.897366596101027*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*wv+24.39471337844178*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119284*geoY[0]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+19.09188309203678*BmagInv[0]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+19.09188309203678*geoY[0]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+8.538149682454623*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+4.74341649025257*BmagInv[0]*geoY[0]*Bmag[1]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119284*geoY[0]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+1.897366596101027*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wv+8.538149682454623*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+4.74341649025257*BmagInv[0]*geoY[0]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+2.121320343559642*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+2.121320343559642*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*wv-(14.90788039793665*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(24.39471337844178*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.485281374238568*BmagInv[0]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(1.897366596101027*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(24.39471337844178*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.485281374238568*geoY[0]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(1.897366596101027*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(19.09188309203678*BmagInv[0]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(19.09188309203678*geoY[0]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(17.07629936490925*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.48683298050514*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)+(28.92857142857142*Apar[1]*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(17.14285714285715*Bmag[1]*Apar[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*Apar[0]*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(2.357142857142857*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(28.92857142857142*BmagInv[1]*Apar[2]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499839*Apar[0]*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499839*BmagInv[0]*Apar[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*BmagInv[0]*Bmag[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(15.42857142857143*Apar[1]*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*Apar[0]*BmagInv[0]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_-(6.000000000000001*Apar[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(2.357142857142857*Bmag[1]*Bmag[1]*BmagInv[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*wm)/m_-(6.000000000000001*Bmag[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x*wm)/m_+(28.92857142857142*geoY[1]*Apar[2]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499839*Apar[0]*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499839*geoY[0]*Apar[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*geoY[0]*Bmag[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(15.42857142857143*Apar[1]*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*Apar[0]*geoY[0]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(2.357142857142857*Bmag[1]*Bmag[1]*geoY[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*Apar[0]*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499839*BmagInv[0]*geoY[1]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499839*geoY[0]*BmagInv[1]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(28.92857142857142*Apar[1]*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(13.5*Apar[0]*BmagInv[0]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(13.5*Apar[0]*geoY[0]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(13.5*BmagInv[0]*geoY[0]*Apar[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(15.42857142857143*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_-(27.0*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*BmagInv[0]*geoY[0]*Bmag[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(12.07476707849886*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(12.07476707849886*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(12.07476707849886*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.70820393249937*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_-(6.70820393249937*geoY[0]*Apar[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*Apar[2]*dfac_x*dfac_x*wm)/m_-(6.70820393249937*geoY[0]*Bmag[1]*Apar[2]*dfac_x*dfac_x*wm)/m_+(2.699999999999999*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/m_+(1.5*Apar[0]*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/m_-(3.0*Apar[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/m_+(1.5*Apar[0]*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm)/m_+(1.5*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wm)/m_-(4.898979485566355*Bmag[2]*Gradpar[2]*dfac_x*wm)/m_-(5.47722557505166*Gradpar[0]*Bmag[2]*dfac_x*wm)/m_-(2.449489742783177*Bmag[1]*Gradpar[1]*dfac_x*wm)/m_+(28.92857142857142*Apar[1]*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(8.571428571428573*Bmag[1]*Apar[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Apar[0]*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(28.92857142857142*BmagInv[1]*Apar[2]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499839*Apar[0]*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499839*BmagInv[0]*Apar[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*BmagInv[0]*Bmag[1]*Apar[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(7.714285714285715*Apar[1]*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*BmagInv[0]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_-(6.000000000000001*Apar[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(8.571428571428573*Phi[1]*Apar[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Apar[0]*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(2.357142857142857*Apar[1]*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*BmagInv[0]*Phi[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(7.714285714285715*Apar[1]*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*BmagInv[0]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(2.357142857142857*Bmag[1]*BmagInv[1]*Phi[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*q_)/m_-(6.000000000000001*Phi[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*Apar[0]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*BmagInv[0]*Apar[1]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_+(28.92857142857142*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499839*Apar[0]*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499839*geoY[0]*Apar[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*geoY[0]*Bmag[1]*Apar[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(7.714285714285715*Apar[1]*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*geoY[0]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499839*BmagInv[0]*geoY[1]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499839*geoY[0]*BmagInv[1]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(28.92857142857142*Apar[1]*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(13.5*Apar[0]*BmagInv[0]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(13.5*Apar[0]*geoY[0]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(13.5*BmagInv[0]*geoY[0]*Apar[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(7.714285714285715*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_-(27.0*geoY[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.0*BmagInv[0]*geoY[0]*Bmag[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(6.037383539249431*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(6.037383539249431*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(6.037383539249431*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.354101966249685*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Phi[2]*dfac_x*dfac_x*q_)/m_-(6.70820393249937*geoY[0]*Apar[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*geoY[0]*Phi[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(7.714285714285715*Apar[1]*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*geoY[0]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(2.357142857142857*Bmag[1]*Phi[1]*geoY[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*Apar[0]*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*geoY[0]*Apar[1]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(7.714285714285715*BmagInv[1]*Phi[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(3.0*BmagInv[0]*geoY[0]*Phi[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(6.037383539249431*Apar[0]*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(6.037383539249431*BmagInv[0]*Apar[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(6.037383539249431*geoY[0]*Apar[1]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(3.354101966249685*Apar[0]*BmagInv[0]*geoY[0]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*Apar[2]*dfac_x*dfac_x*q_)/m_-(6.70820393249937*geoY[0]*Phi[1]*Apar[2]*dfac_x*dfac_x*q_)/m_+(2.699999999999999*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.5*Apar[0]*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_-(3.0*Apar[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.5*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_+(1.5*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(4.898979485566355*Gradpar[2]*Phi[2]*dfac_x*q_)/m_-(5.47722557505166*Gradpar[0]*Phi[2]*dfac_x*q_)/m_-(2.449489742783177*Gradpar[1]*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[1]*q_)/m_-(7.453940198968323*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(24.39471337844178*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*BmagInv[0]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(7.453940198968323*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*BmagInv[0]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(1.897366596101027*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x)/dfac_v-(24.39471337844178*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*geoY[0]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(19.09188309203678*BmagInv[0]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(19.09188309203678*geoY[0]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(8.538149682454623*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.74341649025257*BmagInv[0]*geoY[0]*Bmag[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*geoY[0]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(1.897366596101027*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(8.538149682454623*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(4.74341649025257*BmagInv[0]*geoY[0]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[3] = (9.62299541807677*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051661*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566355*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051661*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566355*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(11.0227038425243*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(6.123724356957943*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051661*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051661*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)-(9.622995418076767*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.47722557505166*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566355*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.47722557505166*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566355*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(11.0227038425243*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(6.123724356957943*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.47722557505166*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.47722557505166*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)+(11.06566670344977*Apar[2]*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.80448531544916*Apar[0]*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*Apar[1]*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.5532833351724883*Bmag[1]*Bmag[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.80448531544916*BmagInv[0]*Apar[2]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.95910003310479*Apar[1]*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.872983346207417*Apar[0]*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*Bmag[1]*BmagInv[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)-(7.745966692414837*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137754*Apar[0]*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137754*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*BmagInv[0]*Bmag[1]*Bmag[1]*Apar[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.80448531544916*geoY[0]*Apar[2]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.95910003310479*Apar[1]*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.872983346207417*Apar[0]*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*Bmag[1]*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137754*Apar[0]*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137754*geoY[0]*Apar[1]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*geoY[0]*Bmag[1]*Bmag[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*Apar[1]*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.95910003310479*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.872983346207417*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(7.794228634059946*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(7.794228634059946*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(7.794228634059946*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(4.330127018922192*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137754*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137754*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)-(8.660254037844386*geoY[0]*Apar[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.971370023173351*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.872983346207417*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)-(3.872983346207417*Apar[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.872983346207417*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.872983346207417*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*dfac_x*dfac_x)/(dfac_m*m_)-(3.872983346207417*Bmag[1]*geoY[1]*Apar[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*dfac_x*dfac_x)/(dfac_m*m_)-(1.732050807568877*geoY[0]*Apar[1]*Bmag[1]*dfac_x*dfac_x)/(dfac_m*m_)-(3.162277660168379*Gradpar[1]*Bmag[2]*dfac_x)/(dfac_m*m_)-(1.414213562373095*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[5] = (8.607068760795467*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(14.08429433584713*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566355*BmagInv[0]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.095445115010332*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(14.08429433584713*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566355*geoY[0]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.095445115010332*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(11.0227038425243*BmagInv[0]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(11.0227038425243*geoY[0]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(9.859006035092987*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.477225575051661*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)-(8.607068760795467*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(14.08429433584713*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566355*BmagInv[0]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.095445115010332*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(14.08429433584713*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566355*geoY[0]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.095445115010332*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(11.0227038425243*BmagInv[0]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(11.0227038425243*geoY[0]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(9.859006035092987*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.477225575051661*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)+(16.70191850155703*Apar[1]*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.897433186107872*Bmag[1]*Apar[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*Apar[0]*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(1.360897063089832*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(16.70191850155703*BmagInv[1]*Apar[2]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.95910003310479*Apar[0]*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.95910003310479*BmagInv[0]*Apar[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*BmagInv[0]*Bmag[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(8.907689867497085*Apar[1]*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137754*Apar[0]*BmagInv[0]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)-(3.464101615137754*Apar[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(1.360897063089832*Bmag[1]*Bmag[1]*BmagInv[1]*Apar[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)-(3.464101615137754*Bmag[1]*Apar[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(16.70191850155703*geoY[1]*Apar[2]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.95910003310479*Apar[0]*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.95910003310479*geoY[0]*Apar[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*geoY[0]*Bmag[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(8.907689867497085*Apar[1]*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137754*Apar[0]*geoY[0]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(1.360897063089832*Bmag[1]*Bmag[1]*geoY[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*Apar[0]*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.95910003310479*BmagInv[0]*geoY[1]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.95910003310479*geoY[0]*BmagInv[1]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(16.70191850155703*Apar[1]*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(7.794228634059946*Apar[0]*BmagInv[0]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(7.794228634059946*Apar[0]*geoY[0]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(7.794228634059946*BmagInv[0]*geoY[0]*Apar[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(8.907689867497085*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)-(15.58845726811989*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137754*BmagInv[0]*geoY[0]*Bmag[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.971370023173351*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.971370023173351*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.971370023173351*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.872983346207417*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)-(3.872983346207417*geoY[0]*Apar[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*Apar[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*Apar[2]*dfac_x*dfac_x)/(dfac_m*m_)-(3.872983346207417*geoY[0]*Bmag[1]*Apar[2]*dfac_x*dfac_x)/(dfac_m*m_)+(1.558845726811989*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*Apar[0]*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)-(1.732050807568877*Apar[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*Apar[0]*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844385*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*dfac_x*dfac_x)/(dfac_m*m_)-(2.828427124746189*Bmag[2]*Gradpar[2]*dfac_x)/(dfac_m*m_)-(3.162277660168379*Gradpar[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(1.414213562373095*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  alpha[7] = (27.10523708715754*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(1.355261854357877*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(16.66751698511147*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(14.90788039793664*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x*wm*wv)/q_+(16.66751698511147*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(14.90788039793664*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(2.121320343559642*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x*wm*wv)/q_+(24.39471337844178*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(9.486832980505138*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(8.485281374238568*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(8.485281374238568*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm*wv)/q_+(1.897366596101027*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm*wv)/q_+27.10523708715754*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+16.66751698511147*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+7.453940198968322*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+1.355261854357877*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+7.453940198968322*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+2.121320343559642*BmagInv[0]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*wv+16.66751698511147*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+7.453940198968322*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+24.39471337844178*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+9.486832980505138*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119284*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+4.242640687119284*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*wv+7.453940198968322*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+2.121320343559642*geoY[0]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x*wv+4.242640687119284*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+4.242640687119284*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+1.897366596101027*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv-(27.10523708715754*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(1.355261854357877*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(16.66751698511147*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(14.90788039793664*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(16.66751698511147*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(14.90788039793664*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(2.121320343559642*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(24.39471337844178*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(9.486832980505138*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.485281374238568*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(8.485281374238568*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm)/(dfac_v*q_)-(1.897366596101027*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/(dfac_v*q_)+(35.55194805194805*Apar[2]*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(19.16629694999821*Apar[0]*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(17.14285714285715*Apar[1]*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(3.214285714285714*Bmag[1]*Bmag[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(0.95831484749991*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(19.16629694999821*BmagInv[0]*Apar[2]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(28.92857142857143*Apar[1]*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(11.78571428571428*Apar[0]*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(17.14285714285715*Bmag[1]*BmagInv[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_-(23.57142857142857*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*Apar[0]*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(0.95831484749991*BmagInv[0]*Bmag[1]*Bmag[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*wm)/m_+(2.357142857142857*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x*wm)/m_+(1.5*Apar[0]*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x*wm)/m_-(3.0*Apar[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x*wm)/m_+(19.16629694999821*geoY[0]*Apar[2]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(28.92857142857143*Apar[1]*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(11.78571428571428*Apar[0]*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(17.14285714285715*Bmag[1]*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*Apar[0]*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*geoY[0]*Apar[1]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(0.95831484749991*geoY[0]*Bmag[1]*Bmag[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(2.357142857142857*Apar[1]*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(1.5*Apar[0]*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x*wm)/m_+(28.92857142857143*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(11.78571428571428*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499838*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499838*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(17.24966725499838*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.70820393249937*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(10.54146332249901*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_-(13.41640786499874*geoY[0]*Apar[2]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(15.42857142857143*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_-(6.000000000000001*Apar[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(6.000000000000001*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Bmag[2]*dfac_x*dfac_x*wm)/m_+(2.357142857142857*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*wm)/m_-(6.000000000000001*Bmag[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*wm)/m_+(1.5*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*Apar[2]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x*wm)/m_+(1.341640786499874*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x*wm)/m_-(2.449489742783177*Bmag[1]*Gradpar[2]*dfac_x*wm)/m_-(4.898979485566355*Gradpar[1]*Bmag[2]*dfac_x*wm)/m_+(35.55194805194805*Apar[2]*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(19.16629694999821*Apar[0]*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(8.571428571428573*Apar[1]*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(19.16629694999821*BmagInv[0]*Apar[2]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(28.92857142857143*Apar[1]*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(11.78571428571428*Apar[0]*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(8.571428571428573*Bmag[1]*BmagInv[1]*Apar[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_-(23.57142857142857*Apar[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Apar[0]*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*BmagInv[0]*Apar[1]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(8.571428571428573*Apar[1]*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.214285714285714*Bmag[1]*Phi[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(0.95831484749991*Apar[0]*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(8.571428571428573*BmagInv[1]*Phi[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Apar[0]*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*BmagInv[0]*Apar[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(0.95831484749991*BmagInv[0]*Bmag[1]*Phi[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(2.357142857142857*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.5*Apar[0]*BmagInv[0]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_-(3.0*Apar[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_+(19.16629694999821*geoY[0]*Apar[2]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(28.92857142857143*Apar[1]*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(11.78571428571428*Apar[0]*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(8.571428571428573*Bmag[1]*geoY[1]*Apar[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Apar[0]*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*geoY[0]*Apar[1]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(28.92857142857143*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(11.78571428571428*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499838*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499838*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(17.24966725499838*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(6.70820393249937*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_-(13.41640786499874*geoY[0]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(7.714285714285715*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_-(6.000000000000001*Apar[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.0*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(8.571428571428573*Phi[1]*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*Apar[0]*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*geoY[0]*Apar[1]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(0.95831484749991*geoY[0]*Bmag[1]*Phi[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(2.357142857142857*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(1.5*Apar[0]*geoY[0]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*BmagInv[0]*Phi[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(5.270731661249505*geoY[0]*BmagInv[1]*Phi[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(7.714285714285715*Apar[1]*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(3.0*Apar[0]*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(3.0*BmagInv[0]*geoY[0]*Apar[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(2.357142857142857*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*q_)/m_-(6.000000000000001*Phi[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*q_)/m_+(1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*Apar[2]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*Apar[0]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*BmagInv[0]*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.341640786499874*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(4.898979485566355*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(2.449489742783177*Phi[1]*Gradpar[2]*dfac_x*q_)/m_-(2.0*dApardt[2]*q_)/m_-(27.10523708715754*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(16.66751698511147*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(7.453940198968322*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(1.355261854357877*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(7.453940198968322*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*BmagInv[0]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x)/dfac_v-(16.66751698511147*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(7.453940198968322*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(24.39471337844178*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(9.486832980505138*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(7.453940198968322*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(2.121320343559642*geoY[0]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(4.242640687119284*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(1.897366596101027*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[13] = (15.64921592871903*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(0.7824607964359515*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(9.622995418076771*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(8.607068760795466*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(9.622995418076771*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(8.607068760795466*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(14.08429433584713*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(5.47722557505166*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566355*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(4.898979485566355*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x*wv)/(dfac_m*q_)+(1.095445115010332*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x*wv)/(dfac_m*q_)-(15.64921592871903*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(0.7824607964359515*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(9.622995418076771*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(8.607068760795466*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(9.622995418076771*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(8.607068760795466*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.224744871391589*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(14.08429433584713*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(5.477225575051659*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566355*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(4.898979485566355*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)-(1.095445115010332*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*dfac_v*q_)+(20.5259267780078*Apar[2]*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(11.06566670344976*Apar[0]*Bmag[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.897433186107877*Apar[1]*Bmag[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(1.855768722395226*Bmag[1]*Bmag[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.5532833351724882*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(11.06566670344976*BmagInv[0]*Apar[2]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(16.70191850155703*Apar[1]*BmagInv[1]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.804485315449162*Apar[0]*BmagInv[0]*Bmag[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.897433186107877*Bmag[1]*BmagInv[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)-(13.60897063089833*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*Apar[0]*Bmag[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.5532833351724882*BmagInv[0]*Bmag[1]*Bmag[1]*Apar[2]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(1.360897063089832*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844386*Apar[0]*BmagInv[0]*Bmag[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)-(1.732050807568878*Apar[1]*Bmag[1]*geoY[2]*dfac_x*dfac_x)/(dfac_m*m_)+(11.06566670344976*geoY[0]*Apar[2]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(16.70191850155703*Apar[1]*geoY[1]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.804485315449162*Apar[0]*geoY[0]*Bmag[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.897433186107877*Bmag[1]*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*Apar[0]*Bmag[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*geoY[0]*Apar[1]*Bmag[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.5532833351724882*geoY[0]*Bmag[1]*Bmag[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(1.360897063089832*Apar[1]*Bmag[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844386*Apar[0]*geoY[0]*Bmag[1]*Bmag[1]*BmagInv[2]*dfac_x*dfac_x)/(dfac_m*m_)+(16.70191850155703*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.804485315449162*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.959100033104788*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.959100033104788*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(9.959100033104788*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.872983346207417*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(6.08611668689737*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)-(7.745966692414835*geoY[0]*Apar[2]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(8.907689867497087*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137756*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)-(3.464101615137756*Apar[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137756*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(3.464101615137756*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Bmag[2]*dfac_x*dfac_x)/(dfac_m*m_)+(1.360897063089832*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*dfac_x*dfac_x)/(dfac_m*m_)-(3.464101615137756*Bmag[1]*geoY[1]*Apar[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.8660254037844386*BmagInv[0]*geoY[0]*Bmag[1]*Bmag[1]*Apar[2]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*Apar[0]*Bmag[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*BmagInv[0]*Apar[1]*Bmag[1]*Bmag[1]*geoY[1]*dfac_x*dfac_x)/(dfac_m*m_)+(0.7745966692414834*geoY[0]*Apar[1]*Bmag[1]*Bmag[1]*BmagInv[1]*dfac_x*dfac_x)/(dfac_m*m_)-(1.414213562373095*Bmag[1]*Gradpar[2]*dfac_x)/(dfac_m*m_)-(2.828427124746189*Gradpar[1]*Bmag[2]*dfac_x)/(dfac_m*m_); 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-1.118033988749895*fr[8]*amax)+1.118033988749895*fl[8]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax+0.3952847075210473*alpha[5]*fr[18]+0.3952847075210473*alpha[5]*fl[18]-0.3061862178478971*alpha[13]*fr[17]+0.3061862178478971*alpha[13]*fl[17]+0.3952847075210473*alpha[3]*fr[14]+0.3952847075210473*alpha[3]*fl[14]+0.1767766952966368*alpha[13]*fr[13]+0.1767766952966368*alpha[13]*fl[13]+0.3952847075210473*alpha[1]*fr[12]+0.3952847075210473*alpha[1]*fl[12]-0.3061862178478971*alpha[7]*fr[11]+0.3061862178478971*alpha[7]*fl[11]-0.3061862178478971*alpha[5]*fr[10]+0.3061862178478971*alpha[5]*fl[10]+0.3952847075210473*alpha[0]*fr[8]+0.3952847075210473*alpha[0]*fl[8]+0.1767766952966368*alpha[7]*fr[7]+0.1767766952966368*alpha[7]*fl[7]-0.3061862178478971*alpha[3]*fr[6]+0.3061862178478971*alpha[3]*fl[6]+0.1767766952966368*alpha[5]*fr[5]+0.1767766952966368*alpha[5]*fl[5]-0.3061862178478971*alpha[1]*fr[4]+0.3061862178478971*alpha[1]*fl[4]+0.1767766952966368*alpha[3]*fr[3]+0.1767766952966368*alpha[3]*fl[3]-0.3061862178478971*alpha[0]*fr[2]+0.3061862178478971*alpha[0]*fl[2]+0.1767766952966368*alpha[1]*fr[1]+0.1767766952966368*alpha[1]*fl[1]+0.1767766952966368*alpha[0]*fr[0]+0.1767766952966368*alpha[0]*fl[0]; 
  Ghat[1] = (-1.118033988749895*fr[12]*amax)+1.118033988749895*fl[12]*amax+0.8660254037844386*fr[4]*amax+0.8660254037844386*fl[4]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax+0.3535533905932737*alpha[13]*fr[18]+0.3952847075210473*alpha[3]*fr[18]+0.3535533905932737*alpha[13]*fl[18]+0.3952847075210473*alpha[3]*fl[18]-0.273861278752583*alpha[5]*fr[17]+0.273861278752583*alpha[5]*fl[17]+0.3952847075210473*alpha[5]*fr[14]+0.3952847075210473*alpha[5]*fl[14]+0.1581138830084189*alpha[5]*fr[13]+0.1581138830084189*alpha[5]*fl[13]-0.273861278752583*fr[10]*alpha[13]+0.273861278752583*fl[10]*alpha[13]+0.1581138830084189*fr[5]*alpha[13]+0.1581138830084189*fl[5]*alpha[13]+0.3535533905932737*alpha[7]*fr[12]+0.3952847075210473*alpha[0]*fr[12]+0.3535533905932737*alpha[7]*fl[12]+0.3952847075210473*alpha[0]*fl[12]-0.273861278752583*alpha[1]*fr[11]+0.273861278752583*alpha[1]*fl[11]-0.3061862178478971*alpha[3]*fr[10]+0.3061862178478971*alpha[3]*fl[10]+0.3952847075210473*alpha[1]*fr[8]+0.3952847075210473*alpha[1]*fl[8]+0.1581138830084189*alpha[1]*fr[7]+0.1581138830084189*alpha[1]*fl[7]-0.273861278752583*fr[4]*alpha[7]+0.273861278752583*fl[4]*alpha[7]+0.1581138830084189*fr[1]*alpha[7]+0.1581138830084189*fl[1]*alpha[7]-0.3061862178478971*alpha[5]*fr[6]+0.3061862178478971*alpha[5]*fl[6]+0.1767766952966368*alpha[3]*fr[5]+0.1767766952966368*alpha[3]*fl[5]+0.1767766952966368*fr[3]*alpha[5]+0.1767766952966368*fl[3]*alpha[5]-0.3061862178478971*alpha[0]*fr[4]+0.3061862178478971*alpha[0]*fl[4]-0.3061862178478971*alpha[1]*fr[2]+0.3061862178478971*alpha[1]*fl[2]+0.1767766952966368*alpha[0]*fr[1]+0.1767766952966368*alpha[0]*fl[1]+0.1767766952966368*fr[0]*alpha[1]+0.1767766952966368*fl[0]*alpha[1]; 
  Ghat[3] = (-1.118033988749895*fr[14]*amax)+1.118033988749895*fl[14]*amax+0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[3]*amax+0.5*fl[3]*amax-0.273861278752583*alpha[5]*fr[19]+0.273861278752583*alpha[5]*fl[19]+0.3952847075210473*alpha[1]*fr[18]+0.3952847075210473*alpha[1]*fl[18]-0.3061862178478971*alpha[7]*fr[17]+0.3061862178478971*alpha[7]*fl[17]-0.273861278752583*alpha[3]*fr[16]+0.273861278752583*alpha[3]*fl[16]+0.1581138830084189*alpha[5]*fr[15]+0.1581138830084189*alpha[5]*fl[15]+0.3952847075210473*alpha[0]*fr[14]+0.3952847075210473*alpha[0]*fl[14]+0.1767766952966368*alpha[7]*fr[13]+0.1767766952966368*alpha[7]*fl[13]-0.3061862178478971*fr[11]*alpha[13]+0.3061862178478971*fl[11]*alpha[13]+0.1767766952966368*fr[7]*alpha[13]+0.1767766952966368*fl[7]*alpha[13]+0.3952847075210473*alpha[5]*fr[12]+0.3952847075210473*alpha[5]*fl[12]-0.3061862178478971*alpha[1]*fr[10]+0.3061862178478971*alpha[1]*fl[10]+0.1581138830084189*alpha[3]*fr[9]+0.1581138830084189*alpha[3]*fl[9]+0.3952847075210473*alpha[3]*fr[8]+0.3952847075210473*alpha[3]*fl[8]-0.3061862178478971*alpha[0]*fr[6]+0.3061862178478971*alpha[0]*fl[6]+0.1767766952966368*alpha[1]*fr[5]+0.1767766952966368*alpha[1]*fl[5]-0.3061862178478971*fr[4]*alpha[5]+0.3061862178478971*fl[4]*alpha[5]+0.1767766952966368*fr[1]*alpha[5]+0.1767766952966368*fl[1]*alpha[5]+0.1767766952966368*alpha[0]*fr[3]+0.1767766952966368*alpha[0]*fl[3]-0.3061862178478971*fr[2]*alpha[3]+0.3061862178478971*fl[2]*alpha[3]+0.1767766952966368*fr[0]*alpha[3]+0.1767766952966368*fl[0]*alpha[3]; 
  Ghat[5] = (-1.118033988749895*fr[18]*amax)+1.118033988749895*fl[18]*amax+0.8660254037844386*fr[10]*amax+0.8660254037844386*fl[10]*amax-0.5*fr[5]*amax+0.5*fl[5]*amax-0.2449489742783177*alpha[13]*fr[19]-0.273861278752583*alpha[3]*fr[19]+0.2449489742783177*alpha[13]*fl[19]+0.273861278752583*alpha[3]*fl[19]+0.3535533905932737*alpha[7]*fr[18]+0.3952847075210473*alpha[0]*fr[18]+0.3535533905932737*alpha[7]*fl[18]+0.3952847075210473*alpha[0]*fl[18]-0.273861278752583*alpha[1]*fr[17]+0.273861278752583*alpha[1]*fl[17]-0.2738612787525829*alpha[5]*fr[16]+0.2738612787525829*alpha[5]*fl[16]+0.1414213562373095*alpha[13]*fr[15]+0.1581138830084189*alpha[3]*fr[15]+0.1414213562373095*alpha[13]*fl[15]+0.1581138830084189*alpha[3]*fl[15]+0.3952847075210473*alpha[1]*fr[14]+0.3952847075210473*alpha[1]*fl[14]+0.1581138830084189*alpha[1]*fr[13]+0.1581138830084189*alpha[1]*fl[13]+0.3535533905932737*fr[12]*alpha[13]+0.3535533905932737*fl[12]*alpha[13]-0.2738612787525829*fr[4]*alpha[13]+0.2738612787525829*fl[4]*alpha[13]+0.1581138830084189*fr[1]*alpha[13]+0.1581138830084189*fl[1]*alpha[13]+0.3952847075210473*alpha[3]*fr[12]+0.3952847075210473*alpha[3]*fl[12]-0.2738612787525829*alpha[5]*fr[11]+0.2738612787525829*alpha[5]*fl[11]-0.273861278752583*alpha[7]*fr[10]-0.3061862178478971*alpha[0]*fr[10]+0.273861278752583*alpha[7]*fl[10]+0.3061862178478971*alpha[0]*fl[10]+0.1581138830084189*alpha[5]*fr[9]+0.1581138830084189*alpha[5]*fl[9]+0.3952847075210473*alpha[5]*fr[8]+0.3952847075210473*alpha[5]*fl[8]+0.1581138830084189*alpha[5]*fr[7]+0.1581138830084189*alpha[5]*fl[7]+0.1581138830084189*fr[5]*alpha[7]+0.1581138830084189*fl[5]*alpha[7]-0.3061862178478971*alpha[1]*fr[6]+0.3061862178478971*alpha[1]*fl[6]+0.1767766952966368*alpha[0]*fr[5]+0.1767766952966368*alpha[0]*fl[5]-0.3061862178478971*fr[2]*alpha[5]+0.3061862178478971*fl[2]*alpha[5]+0.1767766952966368*fr[0]*alpha[5]+0.1767766952966368*fl[0]*alpha[5]-0.3061862178478971*alpha[3]*fr[4]+0.3061862178478971*alpha[3]*fl[4]+0.1767766952966368*alpha[1]*fr[3]+0.1767766952966368*alpha[1]*fl[3]+0.1767766952966368*fr[1]*alpha[3]+0.1767766952966368*fl[1]*alpha[3]; 
  Ghat[7] = 0.8660254037844386*fr[11]*amax+0.8660254037844386*fl[11]*amax-0.5*fr[7]*amax+0.5*fl[7]*amax+0.3535533905932737*alpha[5]*fr[18]+0.3535533905932737*alpha[5]*fl[18]-0.1956151991089878*alpha[13]*fr[17]-0.3061862178478971*alpha[3]*fr[17]+0.1956151991089878*alpha[13]*fl[17]+0.3061862178478971*alpha[3]*fl[17]+0.3952847075210473*alpha[13]*fr[14]+0.3952847075210473*alpha[13]*fl[14]+0.1129384878631564*alpha[13]*fr[13]+0.1767766952966368*alpha[3]*fr[13]+0.1129384878631564*alpha[13]*fl[13]+0.1767766952966368*alpha[3]*fl[13]-0.3061862178478972*fr[6]*alpha[13]+0.3061862178478972*fl[6]*alpha[13]+0.1767766952966368*fr[3]*alpha[13]+0.1767766952966368*fl[3]*alpha[13]+0.3535533905932737*alpha[1]*fr[12]+0.3535533905932737*alpha[1]*fl[12]-0.1956151991089878*alpha[7]*fr[11]-0.3061862178478972*alpha[0]*fr[11]+0.1956151991089878*alpha[7]*fl[11]+0.3061862178478972*alpha[0]*fl[11]-0.273861278752583*alpha[5]*fr[10]+0.273861278752583*alpha[5]*fl[10]+0.3952847075210473*alpha[7]*fr[8]+0.3952847075210473*alpha[7]*fl[8]+0.1129384878631564*alpha[7]*fr[7]+0.1767766952966368*alpha[0]*fr[7]+0.1129384878631564*alpha[7]*fl[7]+0.1767766952966368*alpha[0]*fl[7]-0.3061862178478971*fr[2]*alpha[7]+0.3061862178478971*fl[2]*alpha[7]+0.1767766952966368*fr[0]*alpha[7]+0.1767766952966368*fl[0]*alpha[7]+0.1581138830084189*alpha[5]*fr[5]+0.1581138830084189*alpha[5]*fl[5]-0.273861278752583*alpha[1]*fr[4]+0.273861278752583*alpha[1]*fl[4]+0.1581138830084189*alpha[1]*fr[1]+0.1581138830084189*alpha[1]*fl[1]; 
  Ghat[9] = 0.8660254037844386*fr[16]*amax+0.8660254037844386*fl[16]*amax-0.5*fr[9]*amax+0.5*fl[9]*amax-0.3061862178478971*alpha[1]*fr[19]+0.3061862178478971*alpha[1]*fl[19]+0.3535533905932737*alpha[5]*fr[18]+0.3535533905932737*alpha[5]*fl[18]-0.273861278752583*alpha[13]*fr[17]+0.273861278752583*alpha[13]*fl[17]-0.3061862178478972*alpha[0]*fr[16]+0.3061862178478972*alpha[0]*fl[16]+0.1767766952966368*alpha[1]*fr[15]+0.1767766952966368*alpha[1]*fl[15]+0.3535533905932737*alpha[3]*fr[14]+0.3535533905932737*alpha[3]*fl[14]+0.1581138830084189*alpha[13]*fr[13]+0.1581138830084189*alpha[13]*fl[13]-0.273861278752583*alpha[5]*fr[10]+0.273861278752583*alpha[5]*fl[10]+0.1767766952966368*alpha[0]*fr[9]+0.1767766952966368*alpha[0]*fl[9]-0.273861278752583*alpha[3]*fr[6]+0.273861278752583*alpha[3]*fl[6]+0.1581138830084189*alpha[5]*fr[5]+0.1581138830084189*alpha[5]*fl[5]+0.1581138830084189*alpha[3]*fr[3]+0.1581138830084189*alpha[3]*fl[3]; 
  Ghat[13] = 0.8660254037844387*fr[17]*amax+0.8660254037844387*fl[17]*amax-0.5*fr[13]*amax+0.5*fl[13]*amax-0.2449489742783177*alpha[5]*fr[19]+0.2449489742783177*alpha[5]*fl[19]+0.3535533905932737*alpha[1]*fr[18]+0.3535533905932737*alpha[1]*fl[18]-0.1956151991089878*alpha[7]*fr[17]-0.3061862178478971*alpha[0]*fr[17]+0.1956151991089878*alpha[7]*fl[17]+0.3061862178478971*alpha[0]*fl[17]-0.2738612787525829*alpha[13]*fr[16]+0.2738612787525829*alpha[13]*fl[16]+0.1414213562373095*alpha[5]*fr[15]+0.1414213562373095*alpha[5]*fl[15]+0.3952847075210473*alpha[7]*fr[14]+0.3952847075210473*alpha[7]*fl[14]+0.1129384878631564*alpha[7]*fr[13]+0.1767766952966368*alpha[0]*fr[13]+0.1129384878631564*alpha[7]*fl[13]+0.1767766952966368*alpha[0]*fl[13]-0.1956151991089878*fr[11]*alpha[13]+0.1956151991089878*fl[11]*alpha[13]+0.1581138830084189*fr[9]*alpha[13]+0.1581138830084189*fl[9]*alpha[13]+0.3952847075210473*fr[8]*alpha[13]+0.3952847075210473*fl[8]*alpha[13]+0.1129384878631564*fr[7]*alpha[13]+0.1129384878631564*fl[7]*alpha[13]-0.3061862178478971*fr[2]*alpha[13]+0.3061862178478971*fl[2]*alpha[13]+0.1767766952966368*fr[0]*alpha[13]+0.1767766952966368*fl[0]*alpha[13]+0.3535533905932737*alpha[5]*fr[12]+0.3535533905932737*alpha[5]*fl[12]-0.3061862178478971*alpha[3]*fr[11]+0.3061862178478971*alpha[3]*fl[11]-0.2738612787525829*alpha[1]*fr[10]+0.2738612787525829*alpha[1]*fl[10]+0.1767766952966368*alpha[3]*fr[7]+0.1767766952966368*alpha[3]*fl[7]-0.3061862178478971*fr[6]*alpha[7]+0.3061862178478971*fl[6]*alpha[7]+0.1767766952966368*fr[3]*alpha[7]+0.1767766952966368*fl[3]*alpha[7]+0.1581138830084189*alpha[1]*fr[5]+0.1581138830084189*alpha[1]*fl[5]-0.2738612787525829*fr[4]*alpha[5]+0.2738612787525829*fl[4]*alpha[5]+0.1581138830084189*fr[1]*alpha[5]+0.1581138830084189*fl[1]*alpha[5]; 
  Ghat[15] = 0.8660254037844387*fr[19]*amax+0.8660254037844387*fl[19]*amax-0.5*fr[15]*amax+0.5*fl[15]*amax-0.2738612787525829*alpha[7]*fr[19]-0.3061862178478971*alpha[0]*fr[19]+0.2738612787525829*alpha[7]*fl[19]+0.3061862178478971*alpha[0]*fl[19]+0.3162277660168379*alpha[13]*fr[18]+0.3535533905932737*alpha[3]*fr[18]+0.3162277660168379*alpha[13]*fl[18]+0.3535533905932737*alpha[3]*fl[18]-0.2449489742783177*alpha[5]*fr[17]+0.2449489742783177*alpha[5]*fl[17]-0.3061862178478971*alpha[1]*fr[16]+0.3061862178478971*alpha[1]*fl[16]+0.1581138830084189*alpha[7]*fr[15]+0.1767766952966368*alpha[0]*fr[15]+0.1581138830084189*alpha[7]*fl[15]+0.1767766952966368*alpha[0]*fl[15]+0.3535533905932737*alpha[5]*fr[14]+0.3535533905932737*alpha[5]*fl[14]+0.1414213562373095*alpha[5]*fr[13]+0.1414213562373095*alpha[5]*fl[13]-0.2449489742783178*fr[10]*alpha[13]+0.2449489742783178*fl[10]*alpha[13]+0.1414213562373095*fr[5]*alpha[13]+0.1414213562373095*fl[5]*alpha[13]-0.2738612787525829*alpha[3]*fr[10]+0.2738612787525829*alpha[3]*fl[10]+0.1767766952966368*alpha[1]*fr[9]+0.1767766952966368*alpha[1]*fl[9]-0.2738612787525829*alpha[5]*fr[6]+0.2738612787525829*alpha[5]*fl[6]+0.1581138830084189*alpha[3]*fr[5]+0.1581138830084189*alpha[3]*fl[5]+0.1581138830084189*fr[3]*alpha[5]+0.1581138830084189*fl[3]*alpha[5]; 

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
