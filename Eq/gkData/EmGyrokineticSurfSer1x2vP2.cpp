#include <GyrokineticModDecl.h> 
double EmGyrokineticSurf1x2vSer_X_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  if (alpha0>0) { 
  incr[0] = 0.08333333333333333*(6.708203932499369*alpha[1]*fl[11]+6.708203932499369*alpha[0]*fl[7]+5.196152422706631*alpha[1]*fl[4]+3.0*alpha[1]*fl[2]+5.196152422706631*alpha[0]*fl[1]+3.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.25*(3.872983346207417*alpha[1]*fl[11]+3.872983346207417*alpha[0]*fl[7]+3.0*alpha[1]*fl[4]+1.732050807568877*alpha[1]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01666666666666667*(23.2379000772445*alpha[1]*fl[12]+33.54101966249684*alpha[0]*fl[11]+13.41640786499874*alpha[1]*fl[8]+33.54101966249685*alpha[1]*fl[7]+25.98076211353316*alpha[0]*fl[4]+15.0*alpha[0]*fl[2]+25.98076211353316*alpha[1]*fl[1]+15.0*fl[0]*alpha[1])*dfac_x; 
  incr[3] = 0.08333333333333333*(6.708203932499369*alpha[1]*fl[17]+6.708203932499369*alpha[0]*fl[13]+5.196152422706631*alpha[1]*fl[10]+3.0*alpha[1]*fl[6]+5.196152422706631*alpha[0]*fl[5]+3.0*alpha[0]*fl[3])*dfac_x; 
  incr[4] = -0.05*(13.41640786499874*alpha[1]*fl[12]+19.36491673103708*alpha[0]*fl[11]+7.745966692414834*alpha[1]*fl[8]+19.36491673103709*alpha[1]*fl[7]+15.0*alpha[0]*fl[4]+8.660254037844386*alpha[0]*fl[2]+15.0*alpha[1]*fl[1]+8.660254037844386*fl[0]*alpha[1])*dfac_x; 
  incr[5] = -0.25*(3.872983346207417*alpha[1]*fl[17]+3.872983346207417*alpha[0]*fl[13]+3.0*alpha[1]*fl[10]+1.732050807568877*alpha[1]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[3])*dfac_x; 
  incr[6] = 0.01666666666666667*(23.2379000772445*alpha[1]*fl[18]+33.54101966249685*alpha[0]*fl[17]+13.41640786499874*alpha[1]*fl[14]+33.54101966249684*alpha[1]*fl[13]+25.98076211353316*alpha[0]*fl[10]+15.0*alpha[0]*fl[6]+25.98076211353316*alpha[1]*fl[5]+15.0*alpha[1]*fl[3])*dfac_x; 
  incr[7] = 0.08333333333333333*(15.0*alpha[1]*fl[11]+15.0*alpha[0]*fl[7]+11.61895003862225*alpha[1]*fl[4]+6.708203932499369*alpha[1]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[8] = 0.01666666666666667*(25.98076211353316*alpha[0]*fl[12]+30.0*alpha[1]*fl[11]+15.0*alpha[0]*fl[8]+23.2379000772445*alpha[1]*fl[4]+13.41640786499874*alpha[1]*fl[2])*dfac_x; 
  incr[9] = 0.01666666666666667*(25.98076211353316*alpha[1]*fl[19]+15.0*alpha[1]*fl[16]+25.98076211353316*alpha[0]*fl[15]+15.0*alpha[0]*fl[9])*dfac_x; 
  incr[10] = -0.05*(13.41640786499874*alpha[1]*fl[18]+19.36491673103709*alpha[0]*fl[17]+7.745966692414834*alpha[1]*fl[14]+19.36491673103708*alpha[1]*fl[13]+15.0*alpha[0]*fl[10]+8.660254037844386*alpha[0]*fl[6]+15.0*alpha[1]*fl[5]+8.660254037844386*alpha[1]*fl[3])*dfac_x; 
  incr[11] = 0.01666666666666667*(51.96152422706631*alpha[1]*fl[12]+75.0*alpha[0]*fl[11]+30.0*alpha[1]*fl[8]+75.00000000000001*alpha[1]*fl[7]+58.09475019311126*alpha[0]*fl[4]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*alpha[1]*fl[1]+33.54101966249684*fl[0]*alpha[1])*dfac_x; 
  incr[12] = -0.05*(15.0*alpha[0]*fl[12]+17.32050807568877*alpha[1]*fl[11]+8.660254037844387*alpha[0]*fl[8]+13.41640786499874*alpha[1]*fl[4]+7.745966692414834*alpha[1]*fl[2])*dfac_x; 
  incr[13] = 0.08333333333333333*(15.0*alpha[1]*fl[17]+15.0*alpha[0]*fl[13]+11.61895003862225*alpha[1]*fl[10]+6.708203932499369*alpha[1]*fl[6]+11.61895003862225*alpha[0]*fl[5]+6.708203932499369*alpha[0]*fl[3])*dfac_x; 
  incr[14] = 0.01666666666666667*(25.98076211353316*alpha[0]*fl[18]+30.0*alpha[1]*fl[17]+15.0*alpha[0]*fl[14]+23.2379000772445*alpha[1]*fl[10]+13.41640786499874*alpha[1]*fl[6])*dfac_x; 
  incr[15] = -0.05*(15.0*alpha[1]*fl[19]+8.660254037844386*alpha[1]*fl[16]+15.0*alpha[0]*fl[15]+8.660254037844387*alpha[0]*fl[9])*dfac_x; 
  incr[16] = 0.01666666666666667*(25.98076211353316*alpha[0]*fl[19]+15.0*alpha[0]*fl[16]+25.98076211353316*alpha[1]*fl[15]+15.0*alpha[1]*fl[9])*dfac_x; 
  incr[17] = 0.01666666666666667*(51.96152422706631*alpha[1]*fl[18]+75.0*alpha[0]*fl[17]+30.0*alpha[1]*fl[14]+75.00000000000001*alpha[1]*fl[13]+58.09475019311126*alpha[0]*fl[10]+33.54101966249685*alpha[0]*fl[6]+58.09475019311126*alpha[1]*fl[5]+33.54101966249685*alpha[1]*fl[3])*dfac_x; 
  incr[18] = -0.05*(15.0*alpha[0]*fl[18]+17.32050807568877*alpha[1]*fl[17]+8.660254037844387*alpha[0]*fl[14]+13.41640786499874*alpha[1]*fl[10]+7.745966692414834*alpha[1]*fl[6])*dfac_x; 
  incr[19] = -0.05*(15.0*alpha[0]*fl[19]+8.660254037844387*alpha[0]*fl[16]+15.0*alpha[1]*fl[15]+8.660254037844386*alpha[1]*fl[9])*dfac_x; 

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
  incr[0] = 0.08333333333333333*(6.708203932499369*alpha[1]*fr[11]+6.708203932499369*alpha[0]*fr[7]-5.196152422706631*alpha[1]*fr[4]+3.0*alpha[1]*fr[2]-5.196152422706631*alpha[0]*fr[1]+3.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = -0.25*(3.872983346207417*alpha[1]*fr[11]+3.872983346207417*alpha[0]*fr[7]-3.0*alpha[1]*fr[4]+1.732050807568877*alpha[1]*fr[2]-3.0*alpha[0]*fr[1]+1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01666666666666667*(23.2379000772445*alpha[1]*fr[12]-33.54101966249684*alpha[0]*fr[11]-13.41640786499874*alpha[1]*fr[8]-33.54101966249685*alpha[1]*fr[7]+25.98076211353316*alpha[0]*fr[4]-15.0*alpha[0]*fr[2]+25.98076211353316*alpha[1]*fr[1]-15.0*fr[0]*alpha[1])*dfac_x; 
  incr[3] = 0.08333333333333333*(6.708203932499369*alpha[1]*fr[17]+6.708203932499369*alpha[0]*fr[13]-5.196152422706631*alpha[1]*fr[10]+3.0*alpha[1]*fr[6]-5.196152422706631*alpha[0]*fr[5]+3.0*alpha[0]*fr[3])*dfac_x; 
  incr[4] = 0.05*(13.41640786499874*alpha[1]*fr[12]-19.36491673103708*alpha[0]*fr[11]-7.745966692414834*alpha[1]*fr[8]-19.36491673103709*alpha[1]*fr[7]+15.0*alpha[0]*fr[4]-8.660254037844386*alpha[0]*fr[2]+15.0*alpha[1]*fr[1]-8.660254037844386*fr[0]*alpha[1])*dfac_x; 
  incr[5] = -0.25*(3.872983346207417*alpha[1]*fr[17]+3.872983346207417*alpha[0]*fr[13]-3.0*alpha[1]*fr[10]+1.732050807568877*alpha[1]*fr[6]-3.0*alpha[0]*fr[5]+1.732050807568877*alpha[0]*fr[3])*dfac_x; 
  incr[6] = -0.01666666666666667*(23.2379000772445*alpha[1]*fr[18]-33.54101966249685*alpha[0]*fr[17]-13.41640786499874*alpha[1]*fr[14]-33.54101966249684*alpha[1]*fr[13]+25.98076211353316*alpha[0]*fr[10]-15.0*alpha[0]*fr[6]+25.98076211353316*alpha[1]*fr[5]-15.0*alpha[1]*fr[3])*dfac_x; 
  incr[7] = 0.08333333333333333*(15.0*alpha[1]*fr[11]+15.0*alpha[0]*fr[7]-11.61895003862225*alpha[1]*fr[4]+6.708203932499369*alpha[1]*fr[2]-11.61895003862225*alpha[0]*fr[1]+6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[8] = -0.01666666666666667*(25.98076211353316*alpha[0]*fr[12]-30.0*alpha[1]*fr[11]-15.0*alpha[0]*fr[8]+23.2379000772445*alpha[1]*fr[4]-13.41640786499874*alpha[1]*fr[2])*dfac_x; 
  incr[9] = -0.01666666666666667*(25.98076211353316*alpha[1]*fr[19]-15.0*alpha[1]*fr[16]+25.98076211353316*alpha[0]*fr[15]-15.0*alpha[0]*fr[9])*dfac_x; 
  incr[10] = 0.05*(13.41640786499874*alpha[1]*fr[18]-19.36491673103709*alpha[0]*fr[17]-7.745966692414834*alpha[1]*fr[14]-19.36491673103708*alpha[1]*fr[13]+15.0*alpha[0]*fr[10]-8.660254037844386*alpha[0]*fr[6]+15.0*alpha[1]*fr[5]-8.660254037844386*alpha[1]*fr[3])*dfac_x; 
  incr[11] = -0.01666666666666667*(51.96152422706631*alpha[1]*fr[12]-75.0*alpha[0]*fr[11]-30.0*alpha[1]*fr[8]-75.00000000000001*alpha[1]*fr[7]+58.09475019311126*alpha[0]*fr[4]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*alpha[1]*fr[1]-33.54101966249684*fr[0]*alpha[1])*dfac_x; 
  incr[12] = 0.05*(15.0*alpha[0]*fr[12]-17.32050807568877*alpha[1]*fr[11]-8.660254037844387*alpha[0]*fr[8]+13.41640786499874*alpha[1]*fr[4]-7.745966692414834*alpha[1]*fr[2])*dfac_x; 
  incr[13] = 0.08333333333333333*(15.0*alpha[1]*fr[17]+15.0*alpha[0]*fr[13]-11.61895003862225*alpha[1]*fr[10]+6.708203932499369*alpha[1]*fr[6]-11.61895003862225*alpha[0]*fr[5]+6.708203932499369*alpha[0]*fr[3])*dfac_x; 
  incr[14] = -0.01666666666666667*(25.98076211353316*alpha[0]*fr[18]-30.0*alpha[1]*fr[17]-15.0*alpha[0]*fr[14]+23.2379000772445*alpha[1]*fr[10]-13.41640786499874*alpha[1]*fr[6])*dfac_x; 
  incr[15] = 0.05*(15.0*alpha[1]*fr[19]-8.660254037844386*alpha[1]*fr[16]+15.0*alpha[0]*fr[15]-8.660254037844387*alpha[0]*fr[9])*dfac_x; 
  incr[16] = -0.01666666666666667*(25.98076211353316*alpha[0]*fr[19]-15.0*alpha[0]*fr[16]+25.98076211353316*alpha[1]*fr[15]-15.0*alpha[1]*fr[9])*dfac_x; 
  incr[17] = -0.01666666666666667*(51.96152422706631*alpha[1]*fr[18]-75.0*alpha[0]*fr[17]-30.0*alpha[1]*fr[14]-75.00000000000001*alpha[1]*fr[13]+58.09475019311126*alpha[0]*fr[10]-33.54101966249685*alpha[0]*fr[6]+58.09475019311126*alpha[1]*fr[5]-33.54101966249685*alpha[1]*fr[3])*dfac_x; 
  incr[18] = 0.05*(15.0*alpha[0]*fr[18]-17.32050807568877*alpha[1]*fr[17]-8.660254037844387*alpha[0]*fr[14]+13.41640786499874*alpha[1]*fr[10]-7.745966692414834*alpha[1]*fr[6])*dfac_x; 
  incr[19] = 0.05*(15.0*alpha[0]*fr[19]-8.660254037844387*alpha[0]*fr[16]+15.0*alpha[1]*fr[15]-8.660254037844386*alpha[1]*fr[9])*dfac_x; 

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
