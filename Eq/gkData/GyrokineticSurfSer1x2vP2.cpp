#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x2vSer_X_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[8]; 
  alpha[0] = 1.414213562373095*Gradpar[0]*wv; 
  alpha[1] = (0.8164965809277261*Gradpar[0])/dfac_v; 
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
double GyrokineticSurf1x2vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[8]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(3.872983346207417*Gradpar[0]*Phi[2]*dfac_x*q_)/m_; 
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
  Ghat[0] = -0.1767766952966368*((4.47213595499958*fr[8]-4.47213595499958*fl[8]-3.464101615137754*(fr[2]+fl[2])+2.0*fr[0]-2.0*fl[0])*amax-1.414213562373095*(alpha[1]*favg[1]+alpha[0]*favg[0])); 
  Ghat[1] = -0.01178511301977579*((67.08203932499369*fr[12]-67.08203932499369*fl[12]-51.96152422706631*(fr[4]+fl[4])+30.0*fr[1]-30.0*fl[1])*amax-18.97366596101028*alpha[1]*favg[4]-21.21320343559643*(alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = -0.05892556509887893*((13.41640786499874*fr[14]-13.41640786499874*fl[14]-10.39230484541326*(fr[6]+fl[6])+6.0*fr[3]-6.0*fl[3])*amax-4.242640687119286*(alpha[1]*favg[3]+alpha[0]*favg[2])); 
  Ghat[3] = -0.01178511301977579*((67.0820393249937*fr[18]-67.0820393249937*fl[18]-51.96152422706631*(fr[10]+fl[10])+30.0*fr[5]-30.0*fl[5])*amax-18.97366596101028*alpha[1]*favg[6]-21.21320343559643*(alpha[0]*favg[3]+alpha[1]*favg[2])); 
  Ghat[4] = 0.03535533905932736*((17.32050807568877*(fr[11]+fl[11])-10.0*fr[7]+10.0*fl[7])*amax+7.071067811865476*alpha[0]*favg[4]+6.324555320336761*alpha[1]*favg[1]); 
  Ghat[5] = 0.02635231383473648*((23.2379000772445*(fr[16]+fl[16])-13.41640786499874*fr[9]+13.41640786499874*fl[9])*amax+9.48683298050514*(alpha[1]*favg[7]+alpha[0]*favg[5])); 
  Ghat[6] = 0.01178511301977579*((51.96152422706632*(fr[17]+fl[17])-30.0*fr[13]+30.0*fl[13])*amax+21.21320343559643*alpha[0]*favg[6]+18.97366596101028*alpha[1]*favg[3]); 
  Ghat[7] = 0.04564354645876382*((13.41640786499874*(fr[19]+fl[19])-7.745966692414834*fr[15]+7.745966692414834*fl[15])*amax+5.477225575051662*(alpha[0]*favg[7]+alpha[1]*favg[5])); 
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
double GyrokineticSurf1x2vSer_X_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.05892556509887893*(dfac_v2*wv*((((100.6230589874906*Bmag[2]-25.98076211353316*Bmag[1])*BmagInv[2]+5.0*(9.0*BmagInv[0]-15.58845726811989*BmagInv[1])*Bmag[2]+2.23606797749979*Bmag[1]*(9.0*BmagInv[1]-5.196152422706631*BmagInv[0]))*geoY[2]+(5.0*(9.0*geoY[0]-15.58845726811989*geoY[1])*Bmag[2]+2.23606797749979*Bmag[1]*(9.0*geoY[1]-5.196152422706631*geoY[0]))*BmagInv[2]+2.23606797749979*((27.0*BmagInv[1]-15.58845726811989*BmagInv[0])*geoY[1]+geoY[0]*(9.0*BmagInv[0]-15.58845726811989*BmagInv[1]))*Bmag[2]+Bmag[1]*((9.0*BmagInv[0]-15.58845726811989*BmagInv[1])*geoY[1]+geoY[0]*(9.0*BmagInv[1]-5.196152422706631*BmagInv[0])))*dfac_x*m_*wv+2.0*(6.708203932499369*Gradpar[2]-5.196152422706631*Gradpar[1]+3.0*Gradpar[0])*q_)+(((33.54101966249685*Bmag[2]-8.660254037844386*Bmag[1])*BmagInv[2]+5.0*(3.0*BmagInv[0]-5.196152422706631*BmagInv[1])*Bmag[2]+2.23606797749979*Bmag[1]*(3.0*BmagInv[1]-1.732050807568877*BmagInv[0]))*geoY[2]+(5.0*(3.0*geoY[0]-5.196152422706631*geoY[1])*Bmag[2]+2.23606797749979*Bmag[1]*(3.0*geoY[1]-1.732050807568877*geoY[0]))*BmagInv[2]+2.23606797749979*((9.0*BmagInv[1]-5.196152422706631*BmagInv[0])*geoY[1]+geoY[0]*(3.0*BmagInv[0]-5.196152422706631*BmagInv[1]))*Bmag[2]+Bmag[1]*((3.0*BmagInv[0]-5.196152422706631*BmagInv[1])*geoY[1]+geoY[0]*(3.0*BmagInv[1]-1.732050807568877*BmagInv[0])))*dfac_x*m_))/(dfac_v2*q_); 

  double alpha[8]; 
  alpha[0] = (0.1666666666666667*(dfac_v2*((((142.3024947075771*Bmag[2]-36.74234614174767*Bmag[1])*BmagInv[2]+(63.63961030678928*BmagInv[0]-110.227038425243*BmagInv[1])*Bmag[2]+Bmag[1]*(28.46049894151542*BmagInv[1]-16.43167672515498*BmagInv[0]))*geoY[2]+((63.63961030678928*geoY[0]-110.227038425243*geoY[1])*Bmag[2]+Bmag[1]*(28.46049894151542*geoY[1]-16.43167672515498*geoY[0]))*BmagInv[2]+((85.38149682454625*BmagInv[1]-49.29503017546494*BmagInv[0])*geoY[1]+geoY[0]*(28.46049894151542*BmagInv[0]-49.29503017546494*BmagInv[1]))*Bmag[2]+Bmag[1]*((12.72792206135786*BmagInv[0]-22.0454076850486*BmagInv[1])*geoY[1]+geoY[0]*(12.72792206135786*BmagInv[1]-7.348469228349534*BmagInv[0])))*dfac_x*m_*wv2+(18.97366596101028*Gradpar[2]-14.69693845669907*Gradpar[1]+8.485281374238571*Gradpar[0])*q_*wv)+(((47.43416490252571*Bmag[2]-12.24744871391589*Bmag[1])*BmagInv[2]+(21.21320343559643*BmagInv[0]-36.74234614174767*BmagInv[1])*Bmag[2]+Bmag[1]*(9.48683298050514*BmagInv[1]-5.477225575051662*BmagInv[0]))*geoY[2]+((21.21320343559643*geoY[0]-36.74234614174767*geoY[1])*Bmag[2]+Bmag[1]*(9.48683298050514*geoY[1]-5.477225575051662*geoY[0]))*BmagInv[2]+((28.46049894151542*BmagInv[1]-16.43167672515498*BmagInv[0])*geoY[1]+geoY[0]*(9.48683298050514*BmagInv[0]-16.43167672515498*BmagInv[1]))*Bmag[2]+Bmag[1]*((4.242640687119286*BmagInv[0]-7.348469228349534*BmagInv[1])*geoY[1]+geoY[0]*(4.242640687119286*BmagInv[1]-2.449489742783178*BmagInv[0])))*dfac_x*m_))/(dfac_v2*q_); 
  alpha[1] = (0.3333333333333333*((((82.15838362577493*Bmag[2]-21.21320343559643*Bmag[1])*BmagInv[2]+(36.74234614174767*BmagInv[0]-63.63961030678928*BmagInv[1])*Bmag[2]+Bmag[1]*(16.43167672515498*BmagInv[1]-9.48683298050514*BmagInv[0]))*geoY[2]+((36.74234614174767*geoY[0]-63.63961030678928*geoY[1])*Bmag[2]+Bmag[1]*(16.43167672515498*geoY[1]-9.48683298050514*geoY[0]))*BmagInv[2]+((49.29503017546494*BmagInv[1]-28.46049894151542*BmagInv[0])*geoY[1]+geoY[0]*(16.43167672515498*BmagInv[0]-28.46049894151542*BmagInv[1]))*Bmag[2]+Bmag[1]*((7.348469228349534*BmagInv[0]-12.72792206135786*BmagInv[1])*geoY[1]+geoY[0]*(7.348469228349534*BmagInv[1]-4.242640687119286*BmagInv[0])))*dfac_x*m_*wv+(5.477225575051662*Gradpar[2]-4.242640687119286*Gradpar[1]+2.449489742783178*Gradpar[0])*q_))/(dfac_v*q_); 
  alpha[4] = (0.06666666666666667*(((106.0660171779821*Bmag[2]-27.38612787525831*Bmag[1])*BmagInv[2]+(47.43416490252571*BmagInv[0]-82.15838362577493*BmagInv[1])*Bmag[2]+Bmag[1]*(21.21320343559643*BmagInv[1]-12.24744871391589*BmagInv[0]))*geoY[2]+((47.43416490252571*geoY[0]-82.15838362577493*geoY[1])*Bmag[2]+Bmag[1]*(21.21320343559643*geoY[1]-12.24744871391589*geoY[0]))*BmagInv[2]+((63.63961030678928*BmagInv[1]-36.74234614174767*BmagInv[0])*geoY[1]+geoY[0]*(21.21320343559643*BmagInv[0]-36.74234614174767*BmagInv[1]))*Bmag[2]+Bmag[1]*((9.48683298050514*BmagInv[0]-16.43167672515498*BmagInv[1])*geoY[1]+geoY[0]*(9.48683298050514*BmagInv[1]-5.477225575051662*BmagInv[0])))*dfac_x*m_)/(dfac_v2*q_); 
  if (alpha0>0) { 
  incr[0] = 0.01666666666666667*(25.98076211353316*alpha[4]*fl[12]+33.54101966249684*alpha[1]*fl[11]+15.0*alpha[4]*fl[8]+33.54101966249685*alpha[0]*fl[7]+alpha[1]*(25.98076211353316*fl[4]+15.0*fl[2])+alpha[0]*(25.98076211353316*fl[1]+15.0*fl[0]))*dfac_x; 
  incr[1] = -0.05*(15.0*alpha[4]*fl[12]+19.36491673103708*alpha[1]*fl[11]+8.660254037844386*alpha[4]*fl[8]+19.36491673103709*alpha[0]*fl[7]+alpha[1]*(15.0*fl[4]+8.660254037844386*fl[2])+alpha[0]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[2] = 0.01666666666666667*(23.2379000772445*alpha[1]*fl[12]+(30.0*alpha[4]+33.54101966249684*alpha[0])*fl[11]+alpha[1]*(13.41640786499874*fl[8]+33.54101966249685*fl[7])+(23.2379000772445*alpha[4]+25.98076211353316*alpha[0])*fl[4]+fl[2]*(13.41640786499874*alpha[4]+15.0*alpha[0])+alpha[1]*(25.98076211353316*fl[1]+15.0*fl[0]))*dfac_x; 
  incr[3] = 0.01666666666666667*(25.98076211353316*alpha[4]*fl[18]+33.54101966249685*alpha[1]*fl[17]+15.0*alpha[4]*fl[14]+33.54101966249684*alpha[0]*fl[13]+alpha[1]*(25.98076211353316*fl[10]+15.0*fl[6])+alpha[0]*(25.98076211353316*fl[5]+15.0*fl[3]))*dfac_x; 
  incr[4] = -0.05*(13.41640786499874*alpha[1]*fl[12]+(17.32050807568877*alpha[4]+19.36491673103708*alpha[0])*fl[11]+alpha[1]*(7.745966692414834*fl[8]+19.36491673103709*fl[7])+(13.41640786499874*alpha[4]+15.0*alpha[0])*fl[4]+fl[2]*(7.745966692414834*alpha[4]+8.660254037844386*alpha[0])+alpha[1]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[5] = -0.05*(15.0*alpha[4]*fl[18]+19.36491673103709*alpha[1]*fl[17]+8.660254037844387*alpha[4]*fl[14]+19.36491673103708*alpha[0]*fl[13]+alpha[1]*(15.0*fl[10]+8.660254037844386*fl[6])+alpha[0]*(15.0*fl[5]+8.660254037844386*fl[3]))*dfac_x; 
  incr[6] = 0.01666666666666667*(23.2379000772445*alpha[1]*fl[18]+(30.0*alpha[4]+33.54101966249685*alpha[0])*fl[17]+alpha[1]*(13.41640786499874*fl[14]+33.54101966249684*fl[13])+(23.2379000772445*alpha[4]+25.98076211353316*alpha[0])*fl[10]+(13.41640786499874*alpha[4]+15.0*alpha[0])*fl[6]+alpha[1]*(25.98076211353316*fl[5]+15.0*fl[3]))*dfac_x; 
  incr[7] = 0.08333333333333333*(11.61895003862225*alpha[4]*fl[12]+15.0*alpha[1]*fl[11]+6.708203932499369*alpha[4]*fl[8]+15.0*alpha[0]*fl[7]+alpha[1]*(11.61895003862225*fl[4]+6.708203932499369*fl[2])+alpha[0]*(11.61895003862225*fl[1]+6.708203932499369*fl[0]))*dfac_x; 
  incr[8] = 0.002380952380952381*((116.1895003862225*alpha[4]+181.8653347947321*alpha[0])*fl[12]+210.0*alpha[1]*fl[11]+(67.0820393249937*alpha[4]+105.0*alpha[0])*fl[8]+234.787137637478*alpha[4]*fl[7]+162.6653005407115*alpha[1]*fl[4]+(181.8653347947321*fl[1]+105.0*fl[0])*alpha[4]+93.91485505499116*alpha[1]*fl[2])*dfac_x; 
  incr[9] = 0.01666666666666667*(alpha[1]*(25.98076211353316*fl[19]+15.0*fl[16])+alpha[0]*(25.98076211353316*fl[15]+15.0*fl[9]))*dfac_x; 
  incr[10] = -0.05*(13.41640786499874*alpha[1]*fl[18]+(17.32050807568877*alpha[4]+19.36491673103709*alpha[0])*fl[17]+alpha[1]*(7.745966692414834*fl[14]+19.36491673103708*fl[13])+(13.41640786499874*alpha[4]+15.0*alpha[0])*fl[10]+(7.745966692414834*alpha[4]+8.660254037844386*alpha[0])*fl[6]+alpha[1]*(15.0*fl[5]+8.660254037844386*fl[3]))*dfac_x; 
  incr[11] = 0.01666666666666667*(51.96152422706631*alpha[1]*fl[12]+(67.0820393249937*alpha[4]+75.0*alpha[0])*fl[11]+alpha[1]*(30.0*fl[8]+75.00000000000001*fl[7])+(51.96152422706632*alpha[4]+58.09475019311126*alpha[0])*fl[4]+fl[2]*(30.0*alpha[4]+33.54101966249684*alpha[0])+alpha[1]*(58.09475019311126*fl[1]+33.54101966249684*fl[0]))*dfac_x; 
  incr[12] = -0.007142857142857143*((67.0820393249937*alpha[4]+105.0*alpha[0])*fl[12]+121.2435565298214*alpha[1]*fl[11]+(38.72983346207417*alpha[4]+60.62177826491071*alpha[0])*fl[8]+135.5544171172596*alpha[4]*fl[7]+93.91485505499116*alpha[1]*fl[4]+(105.0*fl[1]+60.62177826491071*fl[0])*alpha[4]+54.22176684690384*alpha[1]*fl[2])*dfac_x; 
  incr[13] = 0.08333333333333333*(11.61895003862225*alpha[4]*fl[18]+15.0*alpha[1]*fl[17]+6.708203932499369*alpha[4]*fl[14]+15.0*alpha[0]*fl[13]+alpha[1]*(11.61895003862225*fl[10]+6.708203932499369*fl[6])+alpha[0]*(11.61895003862225*fl[5]+6.708203932499369*fl[3]))*dfac_x; 
  incr[14] = 0.002380952380952381*((116.1895003862225*alpha[4]+181.8653347947321*alpha[0])*fl[18]+210.0*alpha[1]*fl[17]+(67.0820393249937*alpha[4]+105.0*alpha[0])*fl[14]+234.787137637478*alpha[4]*fl[13]+alpha[1]*(162.6653005407115*fl[10]+93.91485505499116*fl[6])+alpha[4]*(181.8653347947321*fl[5]+105.0*fl[3]))*dfac_x; 
  incr[15] = -0.05*(alpha[1]*(15.0*fl[19]+8.660254037844386*fl[16])+alpha[0]*(15.0*fl[15]+8.660254037844387*fl[9]))*dfac_x; 
  incr[16] = 0.01666666666666667*((23.2379000772445*alpha[4]+25.98076211353316*alpha[0])*fl[19]+(13.41640786499874*alpha[4]+15.0*alpha[0])*fl[16]+alpha[1]*(25.98076211353316*fl[15]+15.0*fl[9]))*dfac_x; 
  incr[17] = 0.01666666666666667*(51.96152422706631*alpha[1]*fl[18]+(67.0820393249937*alpha[4]+75.0*alpha[0])*fl[17]+alpha[1]*(30.0*fl[14]+75.00000000000001*fl[13])+(51.96152422706631*alpha[4]+58.09475019311126*alpha[0])*fl[10]+(30.0*alpha[4]+33.54101966249685*alpha[0])*fl[6]+alpha[1]*(58.09475019311126*fl[5]+33.54101966249685*fl[3]))*dfac_x; 
  incr[18] = -0.007142857142857143*((67.0820393249937*alpha[4]+105.0*alpha[0])*fl[18]+121.2435565298214*alpha[1]*fl[17]+(38.72983346207417*alpha[4]+60.62177826491071*alpha[0])*fl[14]+135.5544171172596*alpha[4]*fl[13]+alpha[1]*(93.91485505499116*fl[10]+54.22176684690384*fl[6])+alpha[4]*(105.0*fl[5]+60.6217782649107*fl[3]))*dfac_x; 
  incr[19] = -0.05*((13.41640786499874*alpha[4]+15.0*alpha[0])*fl[19]+(7.745966692414834*alpha[4]+8.660254037844387*alpha[0])*fl[16]+alpha[1]*(15.0*fl[15]+8.660254037844386*fl[9]))*dfac_x; 

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
  incr[0] = -0.01666666666666667*(25.98076211353316*alpha[4]*fr[12]-33.54101966249684*alpha[1]*fr[11]-15.0*alpha[4]*fr[8]-33.54101966249685*alpha[0]*fr[7]+alpha[1]*(25.98076211353316*fr[4]-15.0*fr[2])+alpha[0]*(25.98076211353316*fr[1]-15.0*fr[0]))*dfac_x; 
  incr[1] = 0.05*(15.0*alpha[4]*fr[12]-19.36491673103708*alpha[1]*fr[11]-8.660254037844386*alpha[4]*fr[8]-19.36491673103709*alpha[0]*fr[7]+alpha[1]*(15.0*fr[4]-8.660254037844386*fr[2])+alpha[0]*(15.0*fr[1]-8.660254037844386*fr[0]))*dfac_x; 
  incr[2] = -0.01666666666666667*(23.2379000772445*alpha[1]*fr[12]+((-30.0*alpha[4])-33.54101966249684*alpha[0])*fr[11]+alpha[1]*((-13.41640786499874*fr[8])-33.54101966249685*fr[7])+(23.2379000772445*alpha[4]+25.98076211353316*alpha[0])*fr[4]+fr[2]*((-13.41640786499874*alpha[4])-15.0*alpha[0])+alpha[1]*(25.98076211353316*fr[1]-15.0*fr[0]))*dfac_x; 
  incr[3] = -0.01666666666666667*(25.98076211353316*alpha[4]*fr[18]-33.54101966249685*alpha[1]*fr[17]-15.0*alpha[4]*fr[14]-33.54101966249684*alpha[0]*fr[13]+alpha[1]*(25.98076211353316*fr[10]-15.0*fr[6])+alpha[0]*(25.98076211353316*fr[5]-15.0*fr[3]))*dfac_x; 
  incr[4] = 0.05*(13.41640786499874*alpha[1]*fr[12]+((-17.32050807568877*alpha[4])-19.36491673103708*alpha[0])*fr[11]+alpha[1]*((-7.745966692414834*fr[8])-19.36491673103709*fr[7])+(13.41640786499874*alpha[4]+15.0*alpha[0])*fr[4]+fr[2]*((-7.745966692414834*alpha[4])-8.660254037844386*alpha[0])+alpha[1]*(15.0*fr[1]-8.660254037844386*fr[0]))*dfac_x; 
  incr[5] = 0.05*(15.0*alpha[4]*fr[18]-19.36491673103709*alpha[1]*fr[17]-8.660254037844387*alpha[4]*fr[14]-19.36491673103708*alpha[0]*fr[13]+alpha[1]*(15.0*fr[10]-8.660254037844386*fr[6])+alpha[0]*(15.0*fr[5]-8.660254037844386*fr[3]))*dfac_x; 
  incr[6] = -0.01666666666666667*(23.2379000772445*alpha[1]*fr[18]+((-30.0*alpha[4])-33.54101966249685*alpha[0])*fr[17]+alpha[1]*((-13.41640786499874*fr[14])-33.54101966249684*fr[13])+(23.2379000772445*alpha[4]+25.98076211353316*alpha[0])*fr[10]+((-13.41640786499874*alpha[4])-15.0*alpha[0])*fr[6]+alpha[1]*(25.98076211353316*fr[5]-15.0*fr[3]))*dfac_x; 
  incr[7] = -0.08333333333333333*(11.61895003862225*alpha[4]*fr[12]-15.0*alpha[1]*fr[11]-6.708203932499369*alpha[4]*fr[8]-15.0*alpha[0]*fr[7]+alpha[1]*(11.61895003862225*fr[4]-6.708203932499369*fr[2])+alpha[0]*(11.61895003862225*fr[1]-6.708203932499369*fr[0]))*dfac_x; 
  incr[8] = -0.002380952380952381*((116.1895003862225*alpha[4]+181.8653347947321*alpha[0])*fr[12]-210.0*alpha[1]*fr[11]+((-67.0820393249937*alpha[4])-105.0*alpha[0])*fr[8]-234.787137637478*alpha[4]*fr[7]+162.6653005407115*alpha[1]*fr[4]+(181.8653347947321*fr[1]-105.0*fr[0])*alpha[4]-93.91485505499116*alpha[1]*fr[2])*dfac_x; 
  incr[9] = -0.01666666666666667*(alpha[1]*(25.98076211353316*fr[19]-15.0*fr[16])+alpha[0]*(25.98076211353316*fr[15]-15.0*fr[9]))*dfac_x; 
  incr[10] = 0.05*(13.41640786499874*alpha[1]*fr[18]+((-17.32050807568877*alpha[4])-19.36491673103709*alpha[0])*fr[17]+alpha[1]*((-7.745966692414834*fr[14])-19.36491673103708*fr[13])+(13.41640786499874*alpha[4]+15.0*alpha[0])*fr[10]+((-7.745966692414834*alpha[4])-8.660254037844386*alpha[0])*fr[6]+alpha[1]*(15.0*fr[5]-8.660254037844386*fr[3]))*dfac_x; 
  incr[11] = -0.01666666666666667*(51.96152422706631*alpha[1]*fr[12]+((-67.0820393249937*alpha[4])-75.0*alpha[0])*fr[11]+alpha[1]*((-30.0*fr[8])-75.00000000000001*fr[7])+(51.96152422706632*alpha[4]+58.09475019311126*alpha[0])*fr[4]+fr[2]*((-30.0*alpha[4])-33.54101966249684*alpha[0])+alpha[1]*(58.09475019311126*fr[1]-33.54101966249684*fr[0]))*dfac_x; 
  incr[12] = 0.007142857142857143*((67.0820393249937*alpha[4]+105.0*alpha[0])*fr[12]-121.2435565298214*alpha[1]*fr[11]+((-38.72983346207417*alpha[4])-60.62177826491071*alpha[0])*fr[8]-135.5544171172596*alpha[4]*fr[7]+93.91485505499116*alpha[1]*fr[4]+(105.0*fr[1]-60.62177826491071*fr[0])*alpha[4]-54.22176684690384*alpha[1]*fr[2])*dfac_x; 
  incr[13] = -0.08333333333333333*(11.61895003862225*alpha[4]*fr[18]-15.0*alpha[1]*fr[17]-6.708203932499369*alpha[4]*fr[14]-15.0*alpha[0]*fr[13]+alpha[1]*(11.61895003862225*fr[10]-6.708203932499369*fr[6])+alpha[0]*(11.61895003862225*fr[5]-6.708203932499369*fr[3]))*dfac_x; 
  incr[14] = -0.002380952380952381*((116.1895003862225*alpha[4]+181.8653347947321*alpha[0])*fr[18]-210.0*alpha[1]*fr[17]+((-67.0820393249937*alpha[4])-105.0*alpha[0])*fr[14]-234.787137637478*alpha[4]*fr[13]+alpha[1]*(162.6653005407115*fr[10]-93.91485505499116*fr[6])+alpha[4]*(181.8653347947321*fr[5]-105.0*fr[3]))*dfac_x; 
  incr[15] = 0.05*(alpha[1]*(15.0*fr[19]-8.660254037844386*fr[16])+alpha[0]*(15.0*fr[15]-8.660254037844387*fr[9]))*dfac_x; 
  incr[16] = -0.01666666666666667*((23.2379000772445*alpha[4]+25.98076211353316*alpha[0])*fr[19]+((-13.41640786499874*alpha[4])-15.0*alpha[0])*fr[16]+alpha[1]*(25.98076211353316*fr[15]-15.0*fr[9]))*dfac_x; 
  incr[17] = -0.01666666666666667*(51.96152422706631*alpha[1]*fr[18]+((-67.0820393249937*alpha[4])-75.0*alpha[0])*fr[17]+alpha[1]*((-30.0*fr[14])-75.00000000000001*fr[13])+(51.96152422706631*alpha[4]+58.09475019311126*alpha[0])*fr[10]+((-30.0*alpha[4])-33.54101966249685*alpha[0])*fr[6]+alpha[1]*(58.09475019311126*fr[5]-33.54101966249685*fr[3]))*dfac_x; 
  incr[18] = 0.007142857142857143*((67.0820393249937*alpha[4]+105.0*alpha[0])*fr[18]-121.2435565298214*alpha[1]*fr[17]+((-38.72983346207417*alpha[4])-60.62177826491071*alpha[0])*fr[14]-135.5544171172596*alpha[4]*fr[13]+alpha[1]*(93.91485505499116*fr[10]-54.22176684690384*fr[6])+alpha[4]*(105.0*fr[5]-60.6217782649107*fr[3]))*dfac_x; 
  incr[19] = 0.05*((13.41640786499874*alpha[4]+15.0*alpha[0])*fr[19]+((-7.745966692414834*alpha[4])-8.660254037844387*alpha[0])*fr[16]+alpha[1]*(15.0*fr[15]-8.660254037844386*fr[9]))*dfac_x; 

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
double GyrokineticSurf1x2vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[20]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.01785714285714286*dfac_x*(3.0*dfac_v*dfac_x*m_*((((55.0*Bmag[2]*Bmag[2]+7.0*Bmag[1]*Bmag[1])*BmagInv[2]+14.0*Bmag[2]*(2.23606797749979*BmagInv[0]*Bmag[2]+2.0*Bmag[1]*BmagInv[1]))*geoY[2]+7.0*(Bmag[2]*(2.0*(2.23606797749979*geoY[0]*Bmag[2]+2.0*Bmag[1]*geoY[1])*BmagInv[2]+(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2])+Bmag[1]*(4.47213595499958*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))*wm+(((55.0*Bmag[2]*BmagInv[2]+14.0*(2.23606797749979*BmagInv[0]*Bmag[2]+Bmag[1]*BmagInv[1]))*Phi[2]+7.0*Phi[1]*(Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2]))*geoY[2]+7.0*((2.0*(2.23606797749979*geoY[0]*Bmag[2]+Bmag[1]*geoY[1])*BmagInv[2]+(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+2.23606797749979*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*(2.0*geoY[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))*q_)*wv-1.0*((24.24871130596428*(2.23606797749979*Gradpar[1]*Bmag[2]+Gradpar[0]*Bmag[1])*dfac_v*q_+3.0*(((55.0*Bmag[2]*Bmag[2]+7.0*Bmag[1]*Bmag[1])*BmagInv[2]+14.0*Bmag[2]*(2.23606797749979*BmagInv[0]*Bmag[2]+2.0*Bmag[1]*BmagInv[1]))*geoY[2]+7.0*(Bmag[2]*(2.0*(2.23606797749979*geoY[0]*Bmag[2]+2.0*Bmag[1]*geoY[1])*BmagInv[2]+(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2])+Bmag[1]*(4.47213595499958*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))*dfac_x*m_)*wm+q_*(24.24871130596428*(2.23606797749979*Gradpar[1]*Phi[2]+Gradpar[0]*Phi[1])*dfac_v*q_+3.0*(((55.0*Bmag[2]*BmagInv[2]+14.0*(2.23606797749979*BmagInv[0]*Bmag[2]+Bmag[1]*BmagInv[1]))*Phi[2]+7.0*Phi[1]*(Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2]))*geoY[2]+7.0*((2.0*(2.23606797749979*geoY[0]*Bmag[2]+Bmag[1]*geoY[1])*BmagInv[2]+(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+2.23606797749979*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*(2.0*geoY[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))*dfac_x*m_))))/(dfac_v*m_*q_); 

  double alpha[8]; 
  alpha[0] = (0.07142857142857142*(dfac_v*dfac_x2*m_*((((165.0*Bmag[2]*Bmag[2]+21.0*Bmag[1]*Bmag[1])*BmagInv[2]+Bmag[2]*(93.91485505499116*BmagInv[0]*Bmag[2]+84.0*Bmag[1]*BmagInv[1]))*geoY[2]+Bmag[2]*((93.91485505499116*geoY[0]*Bmag[2]+84.0*Bmag[1]*geoY[1])*BmagInv[2]+(189.0*BmagInv[1]*geoY[1]+105.0*BmagInv[0]*geoY[0])*Bmag[2])+Bmag[1]*(93.91485505499116*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+21.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*wm+(((Bmag[2]*(165.0*BmagInv[2]+93.91485505499116*BmagInv[0])+42.0*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*(21.0*Bmag[1]*BmagInv[2]+42.0*BmagInv[1]*Bmag[2]))*geoY[2]+((93.91485505499116*geoY[0]*Bmag[2]+42.0*Bmag[1]*geoY[1])*BmagInv[2]+(189.0*BmagInv[1]*geoY[1]+105.0*BmagInv[0]*geoY[0])*Bmag[2]+46.95742752749558*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*(42.0*geoY[1]*BmagInv[2]+46.95742752749558*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+21.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*q_)*wv+(((-54.22176684690384*Gradpar[1]*Bmag[2])-24.24871130596428*Gradpar[0]*Bmag[1])*dfac_v*dfac_x*q_+((((-165.0*Bmag[2]*Bmag[2])-21.0*Bmag[1]*Bmag[1])*BmagInv[2]+Bmag[2]*((-93.91485505499116*BmagInv[0]*Bmag[2])-84.0*Bmag[1]*BmagInv[1]))*geoY[2]+Bmag[2]*(((-93.91485505499116*geoY[0]*Bmag[2])-84.0*Bmag[1]*geoY[1])*BmagInv[2]+((-189.0*BmagInv[1]*geoY[1])-105.0*BmagInv[0]*geoY[0])*Bmag[2])+Bmag[1]*((-93.91485505499116*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2])-21.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*dfac_x2*m_)*wm+(((Bmag[2]*((-165.0*BmagInv[2])-93.91485505499116*BmagInv[0])-42.0*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*((-21.0*Bmag[1]*BmagInv[2])-42.0*BmagInv[1]*Bmag[2]))*geoY[2]+(((-93.91485505499116*geoY[0]*Bmag[2])-42.0*Bmag[1]*geoY[1])*BmagInv[2]+((-189.0*BmagInv[1]*geoY[1])-105.0*BmagInv[0]*geoY[0])*Bmag[2]-46.95742752749558*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*((-42.0*geoY[1]*BmagInv[2])-46.95742752749558*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))-21.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*dfac_x2*m_*q_+((-54.22176684690384*Gradpar[1]*Phi[2])-24.24871130596428*Gradpar[0]*Phi[1])*dfac_v*dfac_x*q2))/(dfac_v*m_*q_); 
  alpha[1] = (0.01428571428571429*(dfac_v*dfac_x2*m_*(((Bmag[2]*(737.9024325749308*Bmag[1]*BmagInv[2]+1207.476707849887*BmagInv[1]*Bmag[2])+Bmag[1]*(420.0*BmagInv[0]*Bmag[2]+93.91485505499116*Bmag[1]*BmagInv[1]))*geoY[2]+(1207.476707849887*geoY[1]*Bmag[2]*Bmag[2]+Bmag[1]*(420.0*geoY[0]*Bmag[2]+93.91485505499116*Bmag[1]*geoY[1]))*BmagInv[2]+945.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]*Bmag[2]+Bmag[1]*((845.2336954949205*BmagInv[1]*geoY[1]+469.5742752749559*BmagInv[0]*geoY[0])*Bmag[2]+105.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*wm+(((368.9512162874654*Bmag[1]*BmagInv[2]+1207.476707849887*BmagInv[1]*Bmag[2]+210.0*BmagInv[0]*Bmag[1])*Phi[2]+Phi[1]*(Bmag[2]*(368.9512162874654*BmagInv[2]+210.0*BmagInv[0])+93.91485505499116*Bmag[1]*BmagInv[1]))*geoY[2]+((1207.476707849887*geoY[1]*Bmag[2]+210.0*geoY[0]*Bmag[1])*BmagInv[2]+945.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(422.6168477474603*BmagInv[1]*geoY[1]+234.787137637478*BmagInv[0]*geoY[0]))*Phi[2]+Phi[1]*((210.0*geoY[0]*Bmag[2]+93.91485505499116*Bmag[1]*geoY[1])*BmagInv[2]+(422.6168477474603*BmagInv[1]*geoY[1]+234.787137637478*BmagInv[0]*geoY[0])*Bmag[2]+105.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*q_)*wv+((Bmag[2]*((-242.4871130596428*Gradpar[2])-271.1088342345192*Gradpar[0])-121.2435565298214*Bmag[1]*Gradpar[1])*dfac_v*dfac_x*q_+((Bmag[2]*((-737.9024325749308*Bmag[1]*BmagInv[2])-1207.476707849887*BmagInv[1]*Bmag[2])+Bmag[1]*((-420.0*BmagInv[0]*Bmag[2])-93.91485505499116*Bmag[1]*BmagInv[1]))*geoY[2]+(Bmag[1]*((-420.0*geoY[0]*Bmag[2])-93.91485505499116*Bmag[1]*geoY[1])-1207.476707849887*geoY[1]*Bmag[2]*Bmag[2])*BmagInv[2]-945.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]*Bmag[2]+Bmag[1]*(((-845.2336954949205*BmagInv[1]*geoY[1])-469.5742752749559*BmagInv[0]*geoY[0])*Bmag[2]-105.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_x2*m_)*wm+((((-368.9512162874654*Bmag[1]*BmagInv[2])-1207.476707849887*BmagInv[1]*Bmag[2]-210.0*BmagInv[0]*Bmag[1])*Phi[2]+Phi[1]*(Bmag[2]*((-368.9512162874654*BmagInv[2])-210.0*BmagInv[0])-93.91485505499116*Bmag[1]*BmagInv[1]))*geoY[2]+(((-1207.476707849887*geoY[1]*Bmag[2])-210.0*geoY[0]*Bmag[1])*BmagInv[2]-945.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*((-422.6168477474603*BmagInv[1]*geoY[1])-234.787137637478*BmagInv[0]*geoY[0]))*Phi[2]+Phi[1]*(((-210.0*geoY[0]*Bmag[2])-93.91485505499116*Bmag[1]*geoY[1])*BmagInv[2]+((-422.6168477474603*BmagInv[1]*geoY[1])-234.787137637478*BmagInv[0]*geoY[0])*Bmag[2]-105.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_x2*m_*q_+(((-242.4871130596428*Gradpar[2])-271.1088342345192*Gradpar[0])*Phi[2]-121.2435565298214*Gradpar[1]*Phi[1])*dfac_v*dfac_x*q2))/(dfac_v*m_*q_); 
  alpha[2] = (0.07142857142857142*(dfac_v*((((95.26279441628824*Bmag[2]*Bmag[2]+12.12435565298214*Bmag[1]*Bmag[1])*BmagInv[2]+Bmag[2]*(54.22176684690384*BmagInv[0]*Bmag[2]+48.49742261192856*Bmag[1]*BmagInv[1]))*geoY[2]+Bmag[2]*((54.22176684690384*geoY[0]*Bmag[2]+48.49742261192856*Bmag[1]*geoY[1])*BmagInv[2]+(109.1192008768392*BmagInv[1]*geoY[1]+60.6217782649107*BmagInv[0]*geoY[0])*Bmag[2])+Bmag[1]*(54.22176684690384*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+12.12435565298214*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*dfac_x2*m_*wv+((-31.30495168499706*Gradpar[1]*Bmag[2])-14.0*Gradpar[0]*Bmag[1])*dfac_x*q_)+((((-95.26279441628824*Bmag[2]*Bmag[2])-12.12435565298214*Bmag[1]*Bmag[1])*BmagInv[2]+Bmag[2]*((-54.22176684690384*BmagInv[0]*Bmag[2])-48.49742261192856*Bmag[1]*BmagInv[1]))*geoY[2]+Bmag[2]*(((-54.22176684690384*geoY[0]*Bmag[2])-48.49742261192856*Bmag[1]*geoY[1])*BmagInv[2]+((-109.1192008768392*BmagInv[1]*geoY[1])-60.6217782649107*BmagInv[0]*geoY[0])*Bmag[2])+Bmag[1]*((-54.22176684690384*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2])-12.12435565298214*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*dfac_x2*m_))/(dfac_m*dfac_v*m_*q_); 
  alpha[3] = (0.01428571428571429*(dfac_v*(((Bmag[2]*(426.0281680828159*Bmag[1]*BmagInv[2]+697.1370023173349*BmagInv[1]*Bmag[2])+Bmag[1]*(242.4871130596428*BmagInv[0]*Bmag[2]+54.22176684690384*Bmag[1]*BmagInv[1]))*geoY[2]+(697.1370023173349*geoY[1]*Bmag[2]*Bmag[2]+Bmag[1]*(242.4871130596428*geoY[0]*Bmag[2]+54.22176684690384*Bmag[1]*geoY[1]))*BmagInv[2]+545.5960043841961*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]*Bmag[2]+Bmag[1]*((487.9959016221344*BmagInv[1]*geoY[1]+271.1088342345192*BmagInv[0]*geoY[0])*Bmag[2]+60.6217782649107*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_x2*m_*wv+(Bmag[2]*((-140.0*Gradpar[2])-156.5247584249853*Gradpar[0])-70.0*Bmag[1]*Gradpar[1])*dfac_x*q_)+((Bmag[2]*((-426.0281680828159*Bmag[1]*BmagInv[2])-697.1370023173349*BmagInv[1]*Bmag[2])+Bmag[1]*((-242.4871130596428*BmagInv[0]*Bmag[2])-54.22176684690384*Bmag[1]*BmagInv[1]))*geoY[2]+(Bmag[1]*((-242.4871130596428*geoY[0]*Bmag[2])-54.22176684690384*Bmag[1]*geoY[1])-697.1370023173349*geoY[1]*Bmag[2]*Bmag[2])*BmagInv[2]-545.5960043841961*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]*Bmag[2]+Bmag[1]*(((-487.9959016221344*BmagInv[1]*geoY[1])-271.1088342345192*BmagInv[0]*geoY[0])*Bmag[2]-60.6217782649107*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_x2*m_))/(dfac_m*dfac_v*m_*q_); 
  alpha[4] = (0.01428571428571429*(dfac_v*dfac_x2*m_*((((1341.640786499874*Bmag[2]*Bmag[2]+67.0820393249937*Bmag[1]*Bmag[1])*BmagInv[2]+825.0*BmagInv[0]*Bmag[2]*Bmag[2]+Bmag[1]*(737.9024325749308*BmagInv[1]*Bmag[2]+105.0*BmagInv[0]*Bmag[1]))*geoY[2]+(825.0*geoY[0]*Bmag[2]*Bmag[2]+Bmag[1]*(737.9024325749308*geoY[1]*Bmag[2]+105.0*geoY[0]*Bmag[1]))*BmagInv[2]+(1207.476707849887*BmagInv[1]*geoY[1]+469.5742752749559*BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+Bmag[1]*(420.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+93.91485505499116*Bmag[1]*BmagInv[1]*geoY[1]))*wm+(((Bmag[2]*(1341.640786499874*BmagInv[2]+825.0*BmagInv[0])+368.9512162874654*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*(67.0820393249937*Bmag[1]*BmagInv[2]+368.9512162874654*BmagInv[1]*Bmag[2]+105.0*BmagInv[0]*Bmag[1]))*geoY[2]+((825.0*geoY[0]*Bmag[2]+368.9512162874654*Bmag[1]*geoY[1])*BmagInv[2]+(1207.476707849887*BmagInv[1]*geoY[1]+469.5742752749559*BmagInv[0]*geoY[0])*Bmag[2]+210.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*((368.9512162874654*geoY[1]*Bmag[2]+105.0*geoY[0]*Bmag[1])*BmagInv[2]+210.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+93.91485505499116*Bmag[1]*BmagInv[1]*geoY[1]))*q_)*wv+(((-121.2435565298214*Bmag[1]*Gradpar[2])-242.4871130596428*Gradpar[1]*Bmag[2])*dfac_v*dfac_x*q_+((((-1341.640786499874*Bmag[2]*Bmag[2])-67.0820393249937*Bmag[1]*Bmag[1])*BmagInv[2]-825.0*BmagInv[0]*Bmag[2]*Bmag[2]+Bmag[1]*((-737.9024325749308*BmagInv[1]*Bmag[2])-105.0*BmagInv[0]*Bmag[1]))*geoY[2]+(Bmag[1]*((-737.9024325749308*geoY[1]*Bmag[2])-105.0*geoY[0]*Bmag[1])-825.0*geoY[0]*Bmag[2]*Bmag[2])*BmagInv[2]+((-1207.476707849887*BmagInv[1]*geoY[1])-469.5742752749559*BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+Bmag[1]*((-420.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2])-93.91485505499116*Bmag[1]*BmagInv[1]*geoY[1]))*dfac_x2*m_)*wm+(((Bmag[2]*((-1341.640786499874*BmagInv[2])-825.0*BmagInv[0])-368.9512162874654*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*((-67.0820393249937*Bmag[1]*BmagInv[2])-368.9512162874654*BmagInv[1]*Bmag[2]-105.0*BmagInv[0]*Bmag[1]))*geoY[2]+(((-825.0*geoY[0]*Bmag[2])-368.9512162874654*Bmag[1]*geoY[1])*BmagInv[2]+((-1207.476707849887*BmagInv[1]*geoY[1])-469.5742752749559*BmagInv[0]*geoY[0])*Bmag[2]-210.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(((-368.9512162874654*geoY[1]*Bmag[2])-105.0*geoY[0]*Bmag[1])*BmagInv[2]-210.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]-93.91485505499116*Bmag[1]*BmagInv[1]*geoY[1]))*dfac_x2*m_*q_+((-242.4871130596428*Gradpar[1]*Phi[2])-121.2435565298214*Phi[1]*Gradpar[2])*dfac_v*dfac_x*q2))/(dfac_v*m_*q_); 
  alpha[6] = (0.004761904761904762*(dfac_v*((((2323.79000772445*Bmag[2]*Bmag[2]+116.1895003862225*Bmag[1]*Bmag[1])*BmagInv[2]+1428.941916244324*BmagInv[0]*Bmag[2]*Bmag[2]+Bmag[1]*(1278.084504248448*BmagInv[1]*Bmag[2]+181.8653347947321*BmagInv[0]*Bmag[1]))*geoY[2]+(1428.941916244324*geoY[0]*Bmag[2]*Bmag[2]+Bmag[1]*(1278.084504248448*geoY[1]*Bmag[2]+181.8653347947321*geoY[0]*Bmag[1]))*BmagInv[2]+(2091.411006952005*BmagInv[1]*geoY[1]+813.3265027035576*BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+Bmag[1]*(727.4613391789285*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+162.6653005407115*Bmag[1]*BmagInv[1]*geoY[1]))*dfac_x2*m_*wv+((-210.0*Bmag[1]*Gradpar[2])-420.0000000000001*Gradpar[1]*Bmag[2])*dfac_x*q_)+((((-2323.79000772445*Bmag[2]*Bmag[2])-116.1895003862225*Bmag[1]*Bmag[1])*BmagInv[2]-1428.941916244324*BmagInv[0]*Bmag[2]*Bmag[2]+Bmag[1]*((-1278.084504248448*BmagInv[1]*Bmag[2])-181.8653347947321*BmagInv[0]*Bmag[1]))*geoY[2]+(Bmag[1]*((-1278.084504248448*geoY[1]*Bmag[2])-181.8653347947321*geoY[0]*Bmag[1])-1428.941916244324*geoY[0]*Bmag[2]*Bmag[2])*BmagInv[2]+((-2091.411006952005*BmagInv[1]*geoY[1])-813.3265027035576*BmagInv[0]*geoY[0])*Bmag[2]*Bmag[2]+Bmag[1]*((-727.4613391789285*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2])-162.6653005407115*Bmag[1]*BmagInv[1]*geoY[1]))*dfac_x2*m_))/(dfac_m*dfac_v*m_*q_); 
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
  Ghat[0] = -0.1767766952966368*((4.47213595499958*fr[8]-4.47213595499958*fl[8]-3.464101615137754*(fr[2]+fl[2])+2.0*fr[0]-2.0*fl[0])*amax-1.414213562373095*(alpha[6]*favg[6]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])); 
  Ghat[1] = -0.01178511301977579*((67.08203932499369*fr[12]-67.08203932499369*fl[12]-51.96152422706631*(fr[4]+fl[4])+30.0*fr[1]-30.0*fl[1])*amax-18.97366596101028*(alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[1]*favg[4]+favg[1]*alpha[4])-21.21320343559643*(alpha[2]*favg[3]+favg[2]*alpha[3]+alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = -0.02041241452319315*((38.72983346207417*fr[14]-38.72983346207417*fl[14]-30.0*(fr[6]+fl[6])+17.32050807568877*fr[3]-17.32050807568877*fl[3])*amax-10.95445115010333*alpha[3]*favg[7]-12.24744871391589*(alpha[4]*favg[6]+favg[4]*alpha[6])-10.95445115010332*alpha[2]*favg[5]-12.24744871391589*(alpha[1]*favg[3]+favg[1]*alpha[3]+alpha[0]*favg[2]+favg[0]*alpha[2])); 
  Ghat[3] = -7.365695637359863e-4*((1073.312629199899*fr[18]-1073.312629199899*fl[18]-831.384387633061*(fr[10]+fl[10])+480.0*fr[5]-480.0*fl[5])*amax+((-271.5290039756345*alpha[6])-303.5786553761646*alpha[2]-82.15838362577493*alpha[0])*favg[7]-303.5786553761646*(alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[3]*(favg[5]+favg[4])+favg[3]*alpha[4])-339.411254969543*(alpha[0]*favg[3]+favg[0]*alpha[3]+alpha[1]*favg[2]+favg[1]*alpha[2])); 
  Ghat[4] = 0.001683587574253684*((363.7306695894642*(fr[11]+fl[11])-210.0*fr[7]+210.0*fl[7])*amax+(94.86832980505142*alpha[6]+148.492424049175*alpha[2])*favg[6]+148.492424049175*favg[2]*alpha[6]+(94.86832980505142*alpha[4]+148.492424049175*alpha[0])*favg[4]+148.492424049175*favg[0]*alpha[4]+132.815661727072*(alpha[3]*favg[3]+alpha[1]*favg[1])); 
  Ghat[5] = 0.01178511301977579*((51.96152422706632*(fr[16]+fl[16])-30.0*fr[9]+30.0*fl[9])*amax+21.21320343559643*alpha[1]*favg[7]+18.97366596101028*alpha[6]*favg[6]+21.21320343559643*alpha[0]*favg[5]+18.97366596101028*(alpha[3]*favg[3]+alpha[2]*favg[2])); 
  Ghat[6] = 0.001304101327393252*((469.5742752749559*(fr[17]+fl[17])-271.1088342345192*fr[13]+271.1088342345192*fl[13])*amax+153.3623161014466*alpha[3]*favg[7]+(122.474487139159*alpha[4]+191.7028951268082*alpha[0])*favg[6]+(171.4642819948225*favg[5]+122.474487139159*favg[4]+191.7028951268082*favg[0])*alpha[6]+191.7028951268082*(alpha[2]*favg[4]+favg[2]*alpha[4])+171.4642819948225*(alpha[1]*favg[3]+favg[1]*alpha[3])); 
  Ghat[7] = 0.01178511301977579*((51.96152422706632*(fr[19]+fl[19])-30.0*fr[15]+30.0*fl[15])*amax+(18.97366596101028*alpha[4]+21.21320343559643*alpha[0])*favg[7]+16.97056274847715*(alpha[3]*favg[6]+favg[3]*alpha[6])+21.21320343559643*alpha[1]*favg[5]+18.97366596101028*(alpha[2]*favg[3]+favg[2]*alpha[3])); 
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
