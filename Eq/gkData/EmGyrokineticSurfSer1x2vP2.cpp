#include <GyrokineticModDecl.h> 
double EmGyrokineticSurf1x2vSer_X_P2_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  incr[0] = 0.08333333333333333*(6.708203932499369*alpha[1]*fl[11]+6.708203932499369*alpha[0]*fl[7]+alpha[1]*(5.196152422706631*fl[4]+3.0*fl[2])+alpha[0]*(5.196152422706631*fl[1]+3.0*fl[0]))*dfac_x; 
  incr[1] = -0.25*(3.872983346207417*(alpha[1]*fl[11]+alpha[0]*fl[7])+alpha[1]*(3.0*fl[4]+1.732050807568877*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.01666666666666667*(23.2379000772445*alpha[1]*fl[12]+33.54101966249684*alpha[0]*fl[11]+alpha[1]*(13.41640786499874*fl[8]+33.54101966249685*fl[7])+alpha[0]*(25.98076211353316*fl[4]+15.0*fl[2])+alpha[1]*(25.98076211353316*fl[1]+15.0*fl[0]))*dfac_x; 
  incr[3] = 0.08333333333333333*(6.708203932499369*alpha[1]*fl[17]+6.708203932499369*alpha[0]*fl[13]+alpha[1]*(5.196152422706631*fl[10]+3.0*fl[6])+alpha[0]*(5.196152422706631*fl[5]+3.0*fl[3]))*dfac_x; 
  incr[4] = -0.05*(13.41640786499874*alpha[1]*fl[12]+19.36491673103708*alpha[0]*fl[11]+alpha[1]*(7.745966692414834*fl[8]+19.36491673103709*fl[7])+alpha[0]*(15.0*fl[4]+8.660254037844386*fl[2])+alpha[1]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[5] = -0.25*(3.872983346207417*(alpha[1]*fl[17]+alpha[0]*fl[13])+alpha[1]*(3.0*fl[10]+1.732050807568877*fl[6])+alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3]))*dfac_x; 
  incr[6] = 0.01666666666666667*(23.2379000772445*alpha[1]*fl[18]+33.54101966249685*alpha[0]*fl[17]+alpha[1]*(13.41640786499874*fl[14]+33.54101966249684*fl[13])+alpha[0]*(25.98076211353316*fl[10]+15.0*fl[6])+alpha[1]*(25.98076211353316*fl[5]+15.0*fl[3]))*dfac_x; 
  incr[7] = 0.08333333333333333*(15.0*alpha[1]*fl[11]+15.0*alpha[0]*fl[7]+alpha[1]*(11.61895003862225*fl[4]+6.708203932499369*fl[2])+alpha[0]*(11.61895003862225*fl[1]+6.708203932499369*fl[0]))*dfac_x; 
  incr[8] = 0.01666666666666667*(25.98076211353316*alpha[0]*fl[12]+30.0*alpha[1]*fl[11]+15.0*alpha[0]*fl[8]+alpha[1]*(23.2379000772445*fl[4]+13.41640786499874*fl[2]))*dfac_x; 
  incr[9] = 0.01666666666666667*(alpha[1]*(25.98076211353316*fl[19]+15.0*fl[16])+alpha[0]*(25.98076211353316*fl[15]+15.0*fl[9]))*dfac_x; 
  incr[10] = -0.05*(13.41640786499874*alpha[1]*fl[18]+19.36491673103709*alpha[0]*fl[17]+alpha[1]*(7.745966692414834*fl[14]+19.36491673103708*fl[13])+alpha[0]*(15.0*fl[10]+8.660254037844386*fl[6])+alpha[1]*(15.0*fl[5]+8.660254037844386*fl[3]))*dfac_x; 
  incr[11] = 0.01666666666666667*(51.96152422706631*alpha[1]*fl[12]+75.0*alpha[0]*fl[11]+alpha[1]*(30.0*fl[8]+75.00000000000001*fl[7])+alpha[0]*(58.09475019311126*fl[4]+33.54101966249684*fl[2])+alpha[1]*(58.09475019311126*fl[1]+33.54101966249684*fl[0]))*dfac_x; 
  incr[12] = -0.05*(15.0*alpha[0]*fl[12]+17.32050807568877*alpha[1]*fl[11]+8.660254037844387*alpha[0]*fl[8]+alpha[1]*(13.41640786499874*fl[4]+7.745966692414834*fl[2]))*dfac_x; 
  incr[13] = 0.08333333333333333*(15.0*alpha[1]*fl[17]+15.0*alpha[0]*fl[13]+alpha[1]*(11.61895003862225*fl[10]+6.708203932499369*fl[6])+alpha[0]*(11.61895003862225*fl[5]+6.708203932499369*fl[3]))*dfac_x; 
  incr[14] = 0.01666666666666667*(25.98076211353316*alpha[0]*fl[18]+30.0*alpha[1]*fl[17]+15.0*alpha[0]*fl[14]+alpha[1]*(23.2379000772445*fl[10]+13.41640786499874*fl[6]))*dfac_x; 
  incr[15] = -0.05*(alpha[1]*(15.0*fl[19]+8.660254037844386*fl[16])+alpha[0]*(15.0*fl[15]+8.660254037844387*fl[9]))*dfac_x; 
  incr[16] = 0.01666666666666667*(alpha[0]*(25.98076211353316*fl[19]+15.0*fl[16])+alpha[1]*(25.98076211353316*fl[15]+15.0*fl[9]))*dfac_x; 
  incr[17] = 0.01666666666666667*(51.96152422706631*alpha[1]*fl[18]+75.0*alpha[0]*fl[17]+alpha[1]*(30.0*fl[14]+75.00000000000001*fl[13])+alpha[0]*(58.09475019311126*fl[10]+33.54101966249685*fl[6])+alpha[1]*(58.09475019311126*fl[5]+33.54101966249685*fl[3]))*dfac_x; 
  incr[18] = -0.05*(15.0*alpha[0]*fl[18]+17.32050807568877*alpha[1]*fl[17]+8.660254037844387*alpha[0]*fl[14]+alpha[1]*(13.41640786499874*fl[10]+7.745966692414834*fl[6]))*dfac_x; 
  incr[19] = -0.05*(alpha[0]*(15.0*fl[19]+8.660254037844387*fl[16])+alpha[1]*(15.0*fl[15]+8.660254037844386*fl[9]))*dfac_x; 

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
  incr[0] = 0.08333333333333333*(6.708203932499369*alpha[1]*fr[11]+6.708203932499369*alpha[0]*fr[7]+alpha[1]*(3.0*fr[2]-5.196152422706631*fr[4])+alpha[0]*(3.0*fr[0]-5.196152422706631*fr[1]))*dfac_x; 
  incr[1] = -0.25*(3.872983346207417*(alpha[1]*fr[11]+alpha[0]*fr[7])+alpha[1]*(1.732050807568877*fr[2]-3.0*fr[4])+alpha[0]*(1.732050807568877*fr[0]-3.0*fr[1]))*dfac_x; 
  incr[2] = -0.01666666666666667*(23.2379000772445*alpha[1]*fr[12]-33.54101966249684*alpha[0]*fr[11]+alpha[1]*((-13.41640786499874*fr[8])-33.54101966249685*fr[7])+alpha[0]*(25.98076211353316*fr[4]-15.0*fr[2])+alpha[1]*(25.98076211353316*fr[1]-15.0*fr[0]))*dfac_x; 
  incr[3] = 0.08333333333333333*(6.708203932499369*alpha[1]*fr[17]+6.708203932499369*alpha[0]*fr[13]+alpha[1]*(3.0*fr[6]-5.196152422706631*fr[10])+alpha[0]*(3.0*fr[3]-5.196152422706631*fr[5]))*dfac_x; 
  incr[4] = 0.05*(13.41640786499874*alpha[1]*fr[12]-19.36491673103708*alpha[0]*fr[11]+alpha[1]*((-7.745966692414834*fr[8])-19.36491673103709*fr[7])+alpha[0]*(15.0*fr[4]-8.660254037844386*fr[2])+alpha[1]*(15.0*fr[1]-8.660254037844386*fr[0]))*dfac_x; 
  incr[5] = -0.25*(3.872983346207417*(alpha[1]*fr[17]+alpha[0]*fr[13])+alpha[1]*(1.732050807568877*fr[6]-3.0*fr[10])+alpha[0]*(1.732050807568877*fr[3]-3.0*fr[5]))*dfac_x; 
  incr[6] = -0.01666666666666667*(23.2379000772445*alpha[1]*fr[18]-33.54101966249685*alpha[0]*fr[17]+alpha[1]*((-13.41640786499874*fr[14])-33.54101966249684*fr[13])+alpha[0]*(25.98076211353316*fr[10]-15.0*fr[6])+alpha[1]*(25.98076211353316*fr[5]-15.0*fr[3]))*dfac_x; 
  incr[7] = 0.08333333333333333*(15.0*alpha[1]*fr[11]+15.0*alpha[0]*fr[7]+alpha[1]*(6.708203932499369*fr[2]-11.61895003862225*fr[4])+alpha[0]*(6.708203932499369*fr[0]-11.61895003862225*fr[1]))*dfac_x; 
  incr[8] = -0.01666666666666667*(25.98076211353316*alpha[0]*fr[12]-30.0*alpha[1]*fr[11]-15.0*alpha[0]*fr[8]+alpha[1]*(23.2379000772445*fr[4]-13.41640786499874*fr[2]))*dfac_x; 
  incr[9] = -0.01666666666666667*(alpha[1]*(25.98076211353316*fr[19]-15.0*fr[16])+alpha[0]*(25.98076211353316*fr[15]-15.0*fr[9]))*dfac_x; 
  incr[10] = 0.05*(13.41640786499874*alpha[1]*fr[18]-19.36491673103709*alpha[0]*fr[17]+alpha[1]*((-7.745966692414834*fr[14])-19.36491673103708*fr[13])+alpha[0]*(15.0*fr[10]-8.660254037844386*fr[6])+alpha[1]*(15.0*fr[5]-8.660254037844386*fr[3]))*dfac_x; 
  incr[11] = -0.01666666666666667*(51.96152422706631*alpha[1]*fr[12]-75.0*alpha[0]*fr[11]+alpha[1]*((-30.0*fr[8])-75.00000000000001*fr[7])+alpha[0]*(58.09475019311126*fr[4]-33.54101966249684*fr[2])+alpha[1]*(58.09475019311126*fr[1]-33.54101966249684*fr[0]))*dfac_x; 
  incr[12] = 0.05*(15.0*alpha[0]*fr[12]-17.32050807568877*alpha[1]*fr[11]-8.660254037844387*alpha[0]*fr[8]+alpha[1]*(13.41640786499874*fr[4]-7.745966692414834*fr[2]))*dfac_x; 
  incr[13] = 0.08333333333333333*(15.0*alpha[1]*fr[17]+15.0*alpha[0]*fr[13]+alpha[1]*(6.708203932499369*fr[6]-11.61895003862225*fr[10])+alpha[0]*(6.708203932499369*fr[3]-11.61895003862225*fr[5]))*dfac_x; 
  incr[14] = -0.01666666666666667*(25.98076211353316*alpha[0]*fr[18]-30.0*alpha[1]*fr[17]-15.0*alpha[0]*fr[14]+alpha[1]*(23.2379000772445*fr[10]-13.41640786499874*fr[6]))*dfac_x; 
  incr[15] = 0.05*(alpha[1]*(15.0*fr[19]-8.660254037844386*fr[16])+alpha[0]*(15.0*fr[15]-8.660254037844387*fr[9]))*dfac_x; 
  incr[16] = -0.01666666666666667*(alpha[0]*(25.98076211353316*fr[19]-15.0*fr[16])+alpha[1]*(25.98076211353316*fr[15]-15.0*fr[9]))*dfac_x; 
  incr[17] = -0.01666666666666667*(51.96152422706631*alpha[1]*fr[18]-75.0*alpha[0]*fr[17]+alpha[1]*((-30.0*fr[14])-75.00000000000001*fr[13])+alpha[0]*(58.09475019311126*fr[10]-33.54101966249685*fr[6])+alpha[1]*(58.09475019311126*fr[5]-33.54101966249685*fr[3]))*dfac_x; 
  incr[18] = 0.05*(15.0*alpha[0]*fr[18]-17.32050807568877*alpha[1]*fr[17]-8.660254037844387*alpha[0]*fr[14]+alpha[1]*(13.41640786499874*fr[10]-7.745966692414834*fr[6]))*dfac_x; 
  incr[19] = 0.05*(alpha[0]*(15.0*fr[19]-8.660254037844387*fr[16])+alpha[1]*(15.0*fr[15]-8.660254037844386*fr[9]))*dfac_x; 

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
double EmGyrokineticSurf1x2vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.3535533905932737*(1.732050807568877*Phi[1]*dfac_x+dApardt[0])*q_)/m_; 

  double alpha[8]; 
  alpha[0] = -(1.0*(2.449489742783178*Phi[1]*dfac_x+1.414213562373095*dApardt[0])*q_)/m_; 
  alpha[1] = -(1.0*(5.477225575051662*Phi[2]*dfac_x+1.414213562373095*dApardt[1])*q_)/m_; 
  alpha[4] = -(1.414213562373095*dApardt[2]*q_)/m_; 
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

//
// //
  Ghat[0] = -0.1767766952966368*((4.47213595499958*fr[8]-4.47213595499958*fl[8]-3.464101615137754*(fr[2]+fl[2])+2.0*fr[0]-2.0*fl[0])*amax-1.414213562373095*(alpha[4]*favg[4]+alpha[1]*favg[1]+alpha[0]*favg[0])); 
  Ghat[1] = -0.01178511301977579*((67.08203932499369*fr[12]-67.08203932499369*fl[12]-51.96152422706631*(fr[4]+fl[4])+30.0*fr[1]-30.0*fl[1])*amax-18.97366596101028*(alpha[1]*favg[4]+favg[1]*alpha[4])-21.21320343559643*(alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = -0.01178511301977579*((67.08203932499369*fr[14]-67.08203932499369*fl[14]-51.96152422706631*(fr[6]+fl[6])+30.0*fr[3]-30.0*fl[3])*amax-21.21320343559643*alpha[4]*favg[6]-21.21320343559643*(alpha[1]*favg[3]+alpha[0]*favg[2])); 
  Ghat[3] = -0.01178511301977579*((67.0820393249937*fr[18]-67.0820393249937*fl[18]-51.96152422706631*(fr[10]+fl[10])+30.0*fr[5]-30.0*fl[5])*amax-18.97366596101028*(alpha[1]*favg[6]+favg[3]*alpha[4])-21.21320343559643*(alpha[0]*favg[3]+alpha[1]*favg[2])); 
  Ghat[4] = 0.005050762722761052*((121.2435565298214*(fr[11]+fl[11])-70.0*fr[7]+70.0*fl[7])*amax+(31.62277660168381*alpha[4]+49.49747468305833*alpha[0])*favg[4]+49.49747468305833*favg[0]*alpha[4]+44.27188724235732*alpha[1]*favg[1]); 
  Ghat[5] = 0.02635231383473648*((23.2379000772445*(fr[16]+fl[16])-13.41640786499874*fr[9]+13.41640786499874*fl[9])*amax+9.48683298050514*(alpha[1]*favg[7]+alpha[0]*favg[5])); 
  Ghat[6] = 0.001683587574253684*((363.7306695894642*(fr[17]+fl[17])-210.0*fr[13]+210.0*fl[13])*amax+(94.86832980505142*alpha[4]+148.492424049175*alpha[0])*favg[6]+148.492424049175*favg[2]*alpha[4]+132.815661727072*alpha[1]*favg[3]); 
  Ghat[7] = 0.009128709291752763*((67.0820393249937*(fr[19]+fl[19])-38.72983346207417*fr[15]+38.72983346207417*fl[15])*amax+(24.49489742783179*alpha[4]+27.38612787525831*alpha[0])*favg[7]+27.38612787525831*alpha[1]*favg[5]); 
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
double EmGyrokineticSurf1x2vSer_X_P2_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  incr[0] = 0.08333333333333333*(6.708203932499369*alpha[1]*fl[11]+6.708203932499369*alpha[0]*fl[7]+alpha[1]*(5.196152422706631*fl[4]+3.0*fl[2])+alpha[0]*(5.196152422706631*fl[1]+3.0*fl[0]))*dfac_x; 
  incr[1] = -0.25*(3.872983346207417*(alpha[1]*fl[11]+alpha[0]*fl[7])+alpha[1]*(3.0*fl[4]+1.732050807568877*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.01666666666666667*(23.2379000772445*alpha[1]*fl[12]+33.54101966249684*alpha[0]*fl[11]+alpha[1]*(13.41640786499874*fl[8]+33.54101966249685*fl[7])+alpha[0]*(25.98076211353316*fl[4]+15.0*fl[2])+alpha[1]*(25.98076211353316*fl[1]+15.0*fl[0]))*dfac_x; 
  incr[3] = 0.08333333333333333*(6.708203932499369*alpha[1]*fl[17]+6.708203932499369*alpha[0]*fl[13]+alpha[1]*(5.196152422706631*fl[10]+3.0*fl[6])+alpha[0]*(5.196152422706631*fl[5]+3.0*fl[3]))*dfac_x; 
  incr[4] = -0.05*(13.41640786499874*alpha[1]*fl[12]+19.36491673103708*alpha[0]*fl[11]+alpha[1]*(7.745966692414834*fl[8]+19.36491673103709*fl[7])+alpha[0]*(15.0*fl[4]+8.660254037844386*fl[2])+alpha[1]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[5] = -0.25*(3.872983346207417*(alpha[1]*fl[17]+alpha[0]*fl[13])+alpha[1]*(3.0*fl[10]+1.732050807568877*fl[6])+alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3]))*dfac_x; 
  incr[6] = 0.01666666666666667*(23.2379000772445*alpha[1]*fl[18]+33.54101966249685*alpha[0]*fl[17]+alpha[1]*(13.41640786499874*fl[14]+33.54101966249684*fl[13])+alpha[0]*(25.98076211353316*fl[10]+15.0*fl[6])+alpha[1]*(25.98076211353316*fl[5]+15.0*fl[3]))*dfac_x; 
  incr[7] = 0.08333333333333333*(15.0*alpha[1]*fl[11]+15.0*alpha[0]*fl[7]+alpha[1]*(11.61895003862225*fl[4]+6.708203932499369*fl[2])+alpha[0]*(11.61895003862225*fl[1]+6.708203932499369*fl[0]))*dfac_x; 
  incr[8] = 0.01666666666666667*(25.98076211353316*alpha[0]*fl[12]+30.0*alpha[1]*fl[11]+15.0*alpha[0]*fl[8]+alpha[1]*(23.2379000772445*fl[4]+13.41640786499874*fl[2]))*dfac_x; 
  incr[9] = 0.01666666666666667*(alpha[1]*(25.98076211353316*fl[19]+15.0*fl[16])+alpha[0]*(25.98076211353316*fl[15]+15.0*fl[9]))*dfac_x; 
  incr[10] = -0.05*(13.41640786499874*alpha[1]*fl[18]+19.36491673103709*alpha[0]*fl[17]+alpha[1]*(7.745966692414834*fl[14]+19.36491673103708*fl[13])+alpha[0]*(15.0*fl[10]+8.660254037844386*fl[6])+alpha[1]*(15.0*fl[5]+8.660254037844386*fl[3]))*dfac_x; 
  incr[11] = 0.01666666666666667*(51.96152422706631*alpha[1]*fl[12]+75.0*alpha[0]*fl[11]+alpha[1]*(30.0*fl[8]+75.00000000000001*fl[7])+alpha[0]*(58.09475019311126*fl[4]+33.54101966249684*fl[2])+alpha[1]*(58.09475019311126*fl[1]+33.54101966249684*fl[0]))*dfac_x; 
  incr[12] = -0.05*(15.0*alpha[0]*fl[12]+17.32050807568877*alpha[1]*fl[11]+8.660254037844387*alpha[0]*fl[8]+alpha[1]*(13.41640786499874*fl[4]+7.745966692414834*fl[2]))*dfac_x; 
  incr[13] = 0.08333333333333333*(15.0*alpha[1]*fl[17]+15.0*alpha[0]*fl[13]+alpha[1]*(11.61895003862225*fl[10]+6.708203932499369*fl[6])+alpha[0]*(11.61895003862225*fl[5]+6.708203932499369*fl[3]))*dfac_x; 
  incr[14] = 0.01666666666666667*(25.98076211353316*alpha[0]*fl[18]+30.0*alpha[1]*fl[17]+15.0*alpha[0]*fl[14]+alpha[1]*(23.2379000772445*fl[10]+13.41640786499874*fl[6]))*dfac_x; 
  incr[15] = -0.05*(alpha[1]*(15.0*fl[19]+8.660254037844386*fl[16])+alpha[0]*(15.0*fl[15]+8.660254037844387*fl[9]))*dfac_x; 
  incr[16] = 0.01666666666666667*(alpha[0]*(25.98076211353316*fl[19]+15.0*fl[16])+alpha[1]*(25.98076211353316*fl[15]+15.0*fl[9]))*dfac_x; 
  incr[17] = 0.01666666666666667*(51.96152422706631*alpha[1]*fl[18]+75.0*alpha[0]*fl[17]+alpha[1]*(30.0*fl[14]+75.00000000000001*fl[13])+alpha[0]*(58.09475019311126*fl[10]+33.54101966249685*fl[6])+alpha[1]*(58.09475019311126*fl[5]+33.54101966249685*fl[3]))*dfac_x; 
  incr[18] = -0.05*(15.0*alpha[0]*fl[18]+17.32050807568877*alpha[1]*fl[17]+8.660254037844387*alpha[0]*fl[14]+alpha[1]*(13.41640786499874*fl[10]+7.745966692414834*fl[6]))*dfac_x; 
  incr[19] = -0.05*(alpha[0]*(15.0*fl[19]+8.660254037844387*fl[16])+alpha[1]*(15.0*fl[15]+8.660254037844386*fl[9]))*dfac_x; 

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
  incr[0] = 0.08333333333333333*(6.708203932499369*alpha[1]*fr[11]+6.708203932499369*alpha[0]*fr[7]+alpha[1]*(3.0*fr[2]-5.196152422706631*fr[4])+alpha[0]*(3.0*fr[0]-5.196152422706631*fr[1]))*dfac_x; 
  incr[1] = -0.25*(3.872983346207417*(alpha[1]*fr[11]+alpha[0]*fr[7])+alpha[1]*(1.732050807568877*fr[2]-3.0*fr[4])+alpha[0]*(1.732050807568877*fr[0]-3.0*fr[1]))*dfac_x; 
  incr[2] = -0.01666666666666667*(23.2379000772445*alpha[1]*fr[12]-33.54101966249684*alpha[0]*fr[11]+alpha[1]*((-13.41640786499874*fr[8])-33.54101966249685*fr[7])+alpha[0]*(25.98076211353316*fr[4]-15.0*fr[2])+alpha[1]*(25.98076211353316*fr[1]-15.0*fr[0]))*dfac_x; 
  incr[3] = 0.08333333333333333*(6.708203932499369*alpha[1]*fr[17]+6.708203932499369*alpha[0]*fr[13]+alpha[1]*(3.0*fr[6]-5.196152422706631*fr[10])+alpha[0]*(3.0*fr[3]-5.196152422706631*fr[5]))*dfac_x; 
  incr[4] = 0.05*(13.41640786499874*alpha[1]*fr[12]-19.36491673103708*alpha[0]*fr[11]+alpha[1]*((-7.745966692414834*fr[8])-19.36491673103709*fr[7])+alpha[0]*(15.0*fr[4]-8.660254037844386*fr[2])+alpha[1]*(15.0*fr[1]-8.660254037844386*fr[0]))*dfac_x; 
  incr[5] = -0.25*(3.872983346207417*(alpha[1]*fr[17]+alpha[0]*fr[13])+alpha[1]*(1.732050807568877*fr[6]-3.0*fr[10])+alpha[0]*(1.732050807568877*fr[3]-3.0*fr[5]))*dfac_x; 
  incr[6] = -0.01666666666666667*(23.2379000772445*alpha[1]*fr[18]-33.54101966249685*alpha[0]*fr[17]+alpha[1]*((-13.41640786499874*fr[14])-33.54101966249684*fr[13])+alpha[0]*(25.98076211353316*fr[10]-15.0*fr[6])+alpha[1]*(25.98076211353316*fr[5]-15.0*fr[3]))*dfac_x; 
  incr[7] = 0.08333333333333333*(15.0*alpha[1]*fr[11]+15.0*alpha[0]*fr[7]+alpha[1]*(6.708203932499369*fr[2]-11.61895003862225*fr[4])+alpha[0]*(6.708203932499369*fr[0]-11.61895003862225*fr[1]))*dfac_x; 
  incr[8] = -0.01666666666666667*(25.98076211353316*alpha[0]*fr[12]-30.0*alpha[1]*fr[11]-15.0*alpha[0]*fr[8]+alpha[1]*(23.2379000772445*fr[4]-13.41640786499874*fr[2]))*dfac_x; 
  incr[9] = -0.01666666666666667*(alpha[1]*(25.98076211353316*fr[19]-15.0*fr[16])+alpha[0]*(25.98076211353316*fr[15]-15.0*fr[9]))*dfac_x; 
  incr[10] = 0.05*(13.41640786499874*alpha[1]*fr[18]-19.36491673103709*alpha[0]*fr[17]+alpha[1]*((-7.745966692414834*fr[14])-19.36491673103708*fr[13])+alpha[0]*(15.0*fr[10]-8.660254037844386*fr[6])+alpha[1]*(15.0*fr[5]-8.660254037844386*fr[3]))*dfac_x; 
  incr[11] = -0.01666666666666667*(51.96152422706631*alpha[1]*fr[12]-75.0*alpha[0]*fr[11]+alpha[1]*((-30.0*fr[8])-75.00000000000001*fr[7])+alpha[0]*(58.09475019311126*fr[4]-33.54101966249684*fr[2])+alpha[1]*(58.09475019311126*fr[1]-33.54101966249684*fr[0]))*dfac_x; 
  incr[12] = 0.05*(15.0*alpha[0]*fr[12]-17.32050807568877*alpha[1]*fr[11]-8.660254037844387*alpha[0]*fr[8]+alpha[1]*(13.41640786499874*fr[4]-7.745966692414834*fr[2]))*dfac_x; 
  incr[13] = 0.08333333333333333*(15.0*alpha[1]*fr[17]+15.0*alpha[0]*fr[13]+alpha[1]*(6.708203932499369*fr[6]-11.61895003862225*fr[10])+alpha[0]*(6.708203932499369*fr[3]-11.61895003862225*fr[5]))*dfac_x; 
  incr[14] = -0.01666666666666667*(25.98076211353316*alpha[0]*fr[18]-30.0*alpha[1]*fr[17]-15.0*alpha[0]*fr[14]+alpha[1]*(23.2379000772445*fr[10]-13.41640786499874*fr[6]))*dfac_x; 
  incr[15] = 0.05*(alpha[1]*(15.0*fr[19]-8.660254037844386*fr[16])+alpha[0]*(15.0*fr[15]-8.660254037844387*fr[9]))*dfac_x; 
  incr[16] = -0.01666666666666667*(alpha[0]*(25.98076211353316*fr[19]-15.0*fr[16])+alpha[1]*(25.98076211353316*fr[15]-15.0*fr[9]))*dfac_x; 
  incr[17] = -0.01666666666666667*(51.96152422706631*alpha[1]*fr[18]-75.0*alpha[0]*fr[17]+alpha[1]*((-30.0*fr[14])-75.00000000000001*fr[13])+alpha[0]*(58.09475019311126*fr[10]-33.54101966249685*fr[6])+alpha[1]*(58.09475019311126*fr[5]-33.54101966249685*fr[3]))*dfac_x; 
  incr[18] = 0.05*(15.0*alpha[0]*fr[18]-17.32050807568877*alpha[1]*fr[17]-8.660254037844387*alpha[0]*fr[14]+alpha[1]*(13.41640786499874*fr[10]-7.745966692414834*fr[6]))*dfac_x; 
  incr[19] = 0.05*(alpha[0]*(15.0*fr[19]-8.660254037844387*fr[16])+alpha[1]*(15.0*fr[15]-8.660254037844386*fr[9]))*dfac_x; 

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
double EmGyrokineticSurf1x2vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.3535533905932737*(1.732050807568877*Bmag[1]*dfac_x*wm+(1.732050807568877*Phi[1]*dfac_x+dApardt[0])*q_))/m_; 

  double alpha[8]; 
  alpha[0] = -(1.0*(2.449489742783178*Bmag[1]*dfac_x*wm+(2.449489742783178*Phi[1]*dfac_x+1.414213562373095*dApardt[0])*q_))/m_; 
  alpha[1] = -(1.0*(5.477225575051662*Bmag[2]*dfac_x*wm+(5.477225575051662*Phi[2]*dfac_x+1.414213562373095*dApardt[1])*q_))/m_; 
  alpha[2] = -(1.414213562373095*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[3] = -(3.16227766016838*Bmag[2]*dfac_x)/(dfac_m*m_); 
  alpha[4] = -(1.414213562373095*dApardt[2]*q_)/m_; 
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

//
// //
  Ghat[0] = -0.1767766952966368*((4.47213595499958*fr[8]-4.47213595499958*fl[8]-3.464101615137754*(fr[2]+fl[2])+2.0*fr[0]-2.0*fl[0])*amax-1.414213562373095*(alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])); 
  Ghat[1] = -0.01178511301977579*((67.08203932499369*fr[12]-67.08203932499369*fl[12]-51.96152422706631*(fr[4]+fl[4])+30.0*fr[1]-30.0*fl[1])*amax-18.97366596101028*(alpha[3]*favg[6]+alpha[1]*favg[4]+favg[1]*alpha[4])-21.21320343559643*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = -0.02041241452319315*((38.72983346207417*fr[14]-38.72983346207417*fl[14]-30.0*(fr[6]+fl[6])+17.32050807568877*fr[3]-17.32050807568877*fl[3])*amax-10.95445115010333*alpha[3]*favg[7]-12.24744871391589*alpha[4]*favg[6]-10.95445115010332*alpha[2]*favg[5]-12.24744871391589*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = -0.01178511301977579*((67.0820393249937*fr[18]-67.0820393249937*fl[18]-51.96152422706631*(fr[10]+fl[10])+30.0*fr[5]-30.0*fl[5])*amax-18.97366596101028*(alpha[2]*favg[7]+alpha[1]*favg[6]+alpha[3]*(favg[5]+favg[4])+favg[3]*alpha[4])-21.21320343559643*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2])); 
  Ghat[4] = 0.001683587574253684*((363.7306695894642*(fr[11]+fl[11])-210.0*fr[7]+210.0*fl[7])*amax+148.492424049175*alpha[2]*favg[6]+(94.86832980505142*alpha[4]+148.492424049175*alpha[0])*favg[4]+148.492424049175*favg[0]*alpha[4]+132.815661727072*(alpha[3]*favg[3]+alpha[1]*favg[1])); 
  Ghat[5] = 0.01178511301977579*((51.96152422706632*(fr[16]+fl[16])-30.0*fr[9]+30.0*fl[9])*amax+21.21320343559643*alpha[1]*favg[7]+21.21320343559643*alpha[0]*favg[5]+18.97366596101028*(alpha[3]*favg[3]+alpha[2]*favg[2])); 
  Ghat[6] = 0.001304101327393252*((469.5742752749559*(fr[17]+fl[17])-271.1088342345192*fr[13]+271.1088342345192*fl[13])*amax+153.3623161014466*alpha[3]*favg[7]+(122.474487139159*alpha[4]+191.7028951268082*alpha[0])*favg[6]+191.7028951268082*(alpha[2]*favg[4]+favg[2]*alpha[4])+171.4642819948225*(alpha[1]*favg[3]+favg[1]*alpha[3])); 
  Ghat[7] = 0.009128709291752763*((67.0820393249937*(fr[19]+fl[19])-38.72983346207417*fr[15]+38.72983346207417*fl[15])*amax+(24.49489742783179*alpha[4]+27.38612787525831*alpha[0])*favg[7]+21.90890230020666*alpha[3]*favg[6]+27.38612787525831*alpha[1]*favg[5]+24.49489742783179*(alpha[2]*favg[3]+favg[2]*alpha[3])); 
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
