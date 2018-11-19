#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x1vSer_X_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[3]; 
  alpha[0] = Gradpar[0]*wv; 
  alpha[1] = (0.5773502691896258*Gradpar[0])/dfac_v; 
  if (alpha0>0) { 
  incr[0] = 0.1178511301977579*(6.708203932499369*alpha[1]*fl[6]+6.708203932499369*alpha[0]*fl[4]+alpha[1]*(5.196152422706631*fl[3]+3.0*fl[2])+alpha[0]*(5.196152422706631*fl[1]+3.0*fl[0]))*dfac_x; 
  incr[1] = -0.3535533905932737*(3.872983346207417*(alpha[1]*fl[6]+alpha[0]*fl[4])+alpha[1]*(3.0*fl[3]+1.732050807568877*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.02357022603955158*(23.2379000772445*alpha[1]*fl[7]+33.54101966249684*alpha[0]*fl[6]+alpha[1]*(13.41640786499874*fl[5]+33.54101966249685*fl[4])+alpha[0]*(25.98076211353316*fl[3]+15.0*fl[2])+alpha[1]*(25.98076211353316*fl[1]+15.0*fl[0]))*dfac_x; 
  incr[3] = -0.07071067811865474*(13.41640786499874*alpha[1]*fl[7]+19.36491673103708*alpha[0]*fl[6]+alpha[1]*(7.745966692414834*fl[5]+19.36491673103709*fl[4])+alpha[0]*(15.0*fl[3]+8.660254037844386*fl[2])+alpha[1]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[4] = 0.1178511301977579*(15.0*alpha[1]*fl[6]+15.0*alpha[0]*fl[4]+alpha[1]*(11.61895003862225*fl[3]+6.708203932499369*fl[2])+alpha[0]*(11.61895003862225*fl[1]+6.708203932499369*fl[0]))*dfac_x; 
  incr[5] = 0.02357022603955158*(25.98076211353316*alpha[0]*fl[7]+30.0*alpha[1]*fl[6]+15.0*alpha[0]*fl[5]+alpha[1]*(23.2379000772445*fl[3]+13.41640786499874*fl[2]))*dfac_x; 
  incr[6] = 0.02357022603955158*(51.96152422706631*alpha[1]*fl[7]+75.0*alpha[0]*fl[6]+alpha[1]*(30.0*fl[5]+75.00000000000001*fl[4])+alpha[0]*(58.09475019311126*fl[3]+33.54101966249684*fl[2])+alpha[1]*(58.09475019311126*fl[1]+33.54101966249684*fl[0]))*dfac_x; 
  incr[7] = -0.07071067811865474*(15.0*alpha[0]*fl[7]+17.32050807568877*alpha[1]*fl[6]+8.660254037844387*alpha[0]*fl[5]+alpha[1]*(13.41640786499874*fl[3]+7.745966692414834*fl[2]))*dfac_x; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = 0.1178511301977579*(6.708203932499369*alpha[1]*fr[6]+6.708203932499369*alpha[0]*fr[4]+alpha[1]*(3.0*fr[2]-5.196152422706631*fr[3])+alpha[0]*(3.0*fr[0]-5.196152422706631*fr[1]))*dfac_x; 
  incr[1] = -0.3535533905932737*(3.872983346207417*(alpha[1]*fr[6]+alpha[0]*fr[4])+alpha[1]*(1.732050807568877*fr[2]-3.0*fr[3])+alpha[0]*(1.732050807568877*fr[0]-3.0*fr[1]))*dfac_x; 
  incr[2] = -0.02357022603955158*(23.2379000772445*alpha[1]*fr[7]-33.54101966249684*alpha[0]*fr[6]+alpha[1]*((-13.41640786499874*fr[5])-33.54101966249685*fr[4])+alpha[0]*(25.98076211353316*fr[3]-15.0*fr[2])+alpha[1]*(25.98076211353316*fr[1]-15.0*fr[0]))*dfac_x; 
  incr[3] = 0.07071067811865474*(13.41640786499874*alpha[1]*fr[7]-19.36491673103708*alpha[0]*fr[6]+alpha[1]*((-7.745966692414834*fr[5])-19.36491673103709*fr[4])+alpha[0]*(15.0*fr[3]-8.660254037844386*fr[2])+alpha[1]*(15.0*fr[1]-8.660254037844386*fr[0]))*dfac_x; 
  incr[4] = 0.1178511301977579*(15.0*alpha[1]*fr[6]+15.0*alpha[0]*fr[4]+alpha[1]*(6.708203932499369*fr[2]-11.61895003862225*fr[3])+alpha[0]*(6.708203932499369*fr[0]-11.61895003862225*fr[1]))*dfac_x; 
  incr[5] = -0.02357022603955158*(25.98076211353316*alpha[0]*fr[7]-30.0*alpha[1]*fr[6]-15.0*alpha[0]*fr[5]+alpha[1]*(23.2379000772445*fr[3]-13.41640786499874*fr[2]))*dfac_x; 
  incr[6] = -0.02357022603955158*(51.96152422706631*alpha[1]*fr[7]-75.0*alpha[0]*fr[6]+alpha[1]*((-30.0*fr[5])-75.00000000000001*fr[4])+alpha[0]*(58.09475019311126*fr[3]-33.54101966249684*fr[2])+alpha[1]*(58.09475019311126*fr[1]-33.54101966249684*fr[0]))*dfac_x; 
  incr[7] = 0.07071067811865474*(15.0*alpha[0]*fr[7]-17.32050807568877*alpha[1]*fr[6]-8.660254037844387*alpha[0]*fr[5]+alpha[1]*(13.41640786499874*fr[3]-7.745966692414834*fr[2]))*dfac_x; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[3]; 
  alpha[0] = -(1.224744871391589*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(2.738612787525831*Gradpar[0]*Phi[2]*dfac_x*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[3]; 
  double favg[3]; 
  favg[0] = 0.7071067811865475*(2.23606797749979*(fr[5]+fl[5])+1.732050807568877*(fl[2]-1.0*fr[2])+fr[0]+fl[0]); 
  favg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[7]+fl[7])+3.0*(fl[3]-1.0*fr[3]))+3.0*(fr[1]+fl[1])); 
  favg[2] = -0.1414213562373095*(8.660254037844387*fr[6]-1.0*(8.660254037844387*fl[6]+5.0*(fr[4]+fl[4]))); 

//
// //
  Ghat[0] = -0.25*((3.16227766016838*fr[5]-3.16227766016838*fl[5]-2.449489742783178*(fr[2]+fl[2])+1.414213562373095*fr[0]-1.414213562373095*fl[0])*amax-1.414213562373095*(alpha[1]*favg[1]+alpha[0]*favg[0])); 
  Ghat[1] = -0.01666666666666667*((47.43416490252569*fr[7]-47.43416490252569*fl[7]-36.74234614174767*(fr[3]+fl[3])+21.21320343559643*fr[1]-21.21320343559643*fl[1])*amax-18.97366596101028*alpha[1]*favg[2]-21.21320343559643*(alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.05*((12.24744871391589*(fr[6]+fl[6])-7.071067811865476*fr[4]+7.071067811865476*fl[4])*amax+7.071067811865476*alpha[0]*favg[2]+6.324555320336761*alpha[1]*favg[1]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_v; 
  incr[1] = 0.7071067811865475*Ghat[1]*dfac_v; 
  incr[2] = -1.224744871391589*Ghat[0]*dfac_v; 
  incr[3] = -1.224744871391589*Ghat[1]*dfac_v; 
  incr[4] = 0.7071067811865475*Ghat[2]*dfac_v; 
  incr[5] = 1.58113883008419*Ghat[0]*dfac_v; 
  incr[6] = -1.224744871391589*Ghat[2]*dfac_v; 
  incr[7] = 1.58113883008419*Ghat[1]*dfac_v; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_X_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.05892556509887893*(dfac_v2*wv*((((100.6230589874906*Bmag[2]-25.98076211353316*Bmag[1])*BmagInv[2]+5.0*(9.0*BmagInv[0]-15.58845726811989*BmagInv[1])*Bmag[2]+2.23606797749979*Bmag[1]*(9.0*BmagInv[1]-5.196152422706631*BmagInv[0]))*geoY[2]+(5.0*(9.0*geoY[0]-15.58845726811989*geoY[1])*Bmag[2]+2.23606797749979*Bmag[1]*(9.0*geoY[1]-5.196152422706631*geoY[0]))*BmagInv[2]+2.23606797749979*((27.0*BmagInv[1]-15.58845726811989*BmagInv[0])*geoY[1]+geoY[0]*(9.0*BmagInv[0]-15.58845726811989*BmagInv[1]))*Bmag[2]+Bmag[1]*((9.0*BmagInv[0]-15.58845726811989*BmagInv[1])*geoY[1]+geoY[0]*(9.0*BmagInv[1]-5.196152422706631*BmagInv[0])))*dfac_x*m_*wv+2.0*(6.708203932499369*Gradpar[2]-5.196152422706631*Gradpar[1]+3.0*Gradpar[0])*q_)+(((33.54101966249685*Bmag[2]-8.660254037844386*Bmag[1])*BmagInv[2]+5.0*(3.0*BmagInv[0]-5.196152422706631*BmagInv[1])*Bmag[2]+2.23606797749979*Bmag[1]*(3.0*BmagInv[1]-1.732050807568877*BmagInv[0]))*geoY[2]+(5.0*(3.0*geoY[0]-5.196152422706631*geoY[1])*Bmag[2]+2.23606797749979*Bmag[1]*(3.0*geoY[1]-1.732050807568877*geoY[0]))*BmagInv[2]+2.23606797749979*((9.0*BmagInv[1]-5.196152422706631*BmagInv[0])*geoY[1]+geoY[0]*(3.0*BmagInv[0]-5.196152422706631*BmagInv[1]))*Bmag[2]+Bmag[1]*((3.0*BmagInv[0]-5.196152422706631*BmagInv[1])*geoY[1]+geoY[0]*(3.0*BmagInv[1]-1.732050807568877*BmagInv[0])))*dfac_x*m_))/(dfac_v2*q_); 

  double alpha[3]; 
  alpha[0] = (0.1666666666666667*(dfac_v2*((((100.6230589874906*Bmag[2]-25.98076211353316*Bmag[1])*BmagInv[2]+(45.0*BmagInv[0]-77.94228634059945*BmagInv[1])*Bmag[2]+Bmag[1]*(20.12461179749811*BmagInv[1]-11.61895003862225*BmagInv[0]))*geoY[2]+((45.0*geoY[0]-77.94228634059945*geoY[1])*Bmag[2]+Bmag[1]*(20.12461179749811*geoY[1]-11.61895003862225*geoY[0]))*BmagInv[2]+((60.37383539249433*BmagInv[1]-34.85685011586674*BmagInv[0])*geoY[1]+geoY[0]*(20.12461179749811*BmagInv[0]-34.85685011586674*BmagInv[1]))*Bmag[2]+Bmag[1]*((9.0*BmagInv[0]-15.58845726811989*BmagInv[1])*geoY[1]+geoY[0]*(9.0*BmagInv[1]-5.196152422706631*BmagInv[0])))*dfac_x*m_*wv2+(13.41640786499874*Gradpar[2]-10.39230484541326*Gradpar[1]+6.0*Gradpar[0])*q_*wv)+(((33.54101966249685*Bmag[2]-8.660254037844386*Bmag[1])*BmagInv[2]+(15.0*BmagInv[0]-25.98076211353316*BmagInv[1])*Bmag[2]+Bmag[1]*(6.708203932499369*BmagInv[1]-3.872983346207417*BmagInv[0]))*geoY[2]+((15.0*geoY[0]-25.98076211353316*geoY[1])*Bmag[2]+Bmag[1]*(6.708203932499369*geoY[1]-3.872983346207417*geoY[0]))*BmagInv[2]+((20.12461179749811*BmagInv[1]-11.61895003862225*BmagInv[0])*geoY[1]+geoY[0]*(6.708203932499369*BmagInv[0]-11.61895003862225*BmagInv[1]))*Bmag[2]+Bmag[1]*((3.0*BmagInv[0]-5.196152422706631*BmagInv[1])*geoY[1]+geoY[0]*(3.0*BmagInv[1]-1.732050807568877*BmagInv[0])))*dfac_x*m_))/(dfac_v2*q_); 
  alpha[1] = (0.3333333333333333*((((58.09475019311126*Bmag[2]-15.0*Bmag[1])*BmagInv[2]+(25.98076211353316*BmagInv[0]-45.0*BmagInv[1])*Bmag[2]+Bmag[1]*(11.61895003862225*BmagInv[1]-6.708203932499369*BmagInv[0]))*geoY[2]+((25.98076211353316*geoY[0]-45.0*geoY[1])*Bmag[2]+Bmag[1]*(11.61895003862225*geoY[1]-6.708203932499369*geoY[0]))*BmagInv[2]+((34.85685011586674*BmagInv[1]-20.12461179749811*BmagInv[0])*geoY[1]+geoY[0]*(11.61895003862225*BmagInv[0]-20.12461179749811*BmagInv[1]))*Bmag[2]+Bmag[1]*((5.196152422706631*BmagInv[0]-9.0*BmagInv[1])*geoY[1]+geoY[0]*(5.196152422706631*BmagInv[1]-3.0*BmagInv[0])))*dfac_x*m_*wv+(3.872983346207417*Gradpar[2]-3.0*Gradpar[1]+1.732050807568877*Gradpar[0])*q_))/(dfac_v*q_); 
  alpha[2] = (0.06666666666666667*(((75.0*Bmag[2]-19.36491673103709*Bmag[1])*BmagInv[2]+(33.54101966249685*BmagInv[0]-58.09475019311126*BmagInv[1])*Bmag[2]+Bmag[1]*(15.0*BmagInv[1]-8.660254037844386*BmagInv[0]))*geoY[2]+((33.54101966249685*geoY[0]-58.09475019311126*geoY[1])*Bmag[2]+Bmag[1]*(15.0*geoY[1]-8.660254037844386*geoY[0]))*BmagInv[2]+((45.0*BmagInv[1]-25.98076211353316*BmagInv[0])*geoY[1]+geoY[0]*(15.0*BmagInv[0]-25.98076211353316*BmagInv[1]))*Bmag[2]+Bmag[1]*((6.708203932499369*BmagInv[0]-11.61895003862225*BmagInv[1])*geoY[1]+geoY[0]*(6.708203932499369*BmagInv[1]-3.872983346207417*BmagInv[0])))*dfac_x*m_)/(dfac_v2*q_); 
  if (alpha0>0) { 
  incr[0] = 0.02357022603955158*(25.98076211353316*alpha[2]*fl[7]+33.54101966249684*alpha[1]*fl[6]+15.0*alpha[2]*fl[5]+33.54101966249685*alpha[0]*fl[4]+alpha[1]*(25.98076211353316*fl[3]+15.0*fl[2])+alpha[0]*(25.98076211353316*fl[1]+15.0*fl[0]))*dfac_x; 
  incr[1] = -0.07071067811865474*(15.0*alpha[2]*fl[7]+19.36491673103708*alpha[1]*fl[6]+8.660254037844386*alpha[2]*fl[5]+19.36491673103709*alpha[0]*fl[4]+alpha[1]*(15.0*fl[3]+8.660254037844386*fl[2])+alpha[0]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[2] = 0.02357022603955158*(23.2379000772445*alpha[1]*fl[7]+(30.0*alpha[2]+33.54101966249684*alpha[0])*fl[6]+alpha[1]*(13.41640786499874*fl[5]+33.54101966249685*fl[4])+(23.2379000772445*alpha[2]+25.98076211353316*alpha[0])*fl[3]+(13.41640786499874*alpha[2]+15.0*alpha[0])*fl[2]+alpha[1]*(25.98076211353316*fl[1]+15.0*fl[0]))*dfac_x; 
  incr[3] = -0.07071067811865474*(13.41640786499874*alpha[1]*fl[7]+(17.32050807568877*alpha[2]+19.36491673103708*alpha[0])*fl[6]+alpha[1]*(7.745966692414834*fl[5]+19.36491673103709*fl[4])+(13.41640786499874*alpha[2]+15.0*alpha[0])*fl[3]+(7.745966692414834*alpha[2]+8.660254037844386*alpha[0])*fl[2]+alpha[1]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[4] = 0.1178511301977579*(11.61895003862225*alpha[2]*fl[7]+15.0*alpha[1]*fl[6]+6.708203932499369*alpha[2]*fl[5]+15.0*alpha[0]*fl[4]+alpha[1]*(11.61895003862225*fl[3]+6.708203932499369*fl[2])+alpha[0]*(11.61895003862225*fl[1]+6.708203932499369*fl[0]))*dfac_x; 
  incr[5] = 0.003367175148507369*((116.1895003862225*alpha[2]+181.8653347947321*alpha[0])*fl[7]+210.0*alpha[1]*fl[6]+(67.0820393249937*alpha[2]+105.0*alpha[0])*fl[5]+234.787137637478*alpha[2]*fl[4]+alpha[1]*(162.6653005407115*fl[3]+93.91485505499116*fl[2])+(181.8653347947321*fl[1]+105.0*fl[0])*alpha[2])*dfac_x; 
  incr[6] = 0.02357022603955158*(51.96152422706631*alpha[1]*fl[7]+(67.0820393249937*alpha[2]+75.0*alpha[0])*fl[6]+alpha[1]*(30.0*fl[5]+75.00000000000001*fl[4])+(51.96152422706632*alpha[2]+58.09475019311126*alpha[0])*fl[3]+(30.0*alpha[2]+33.54101966249684*alpha[0])*fl[2]+alpha[1]*(58.09475019311126*fl[1]+33.54101966249684*fl[0]))*dfac_x; 
  incr[7] = -0.01010152544552211*((67.0820393249937*alpha[2]+105.0*alpha[0])*fl[7]+121.2435565298214*alpha[1]*fl[6]+(38.72983346207417*alpha[2]+60.62177826491071*alpha[0])*fl[5]+135.5544171172596*alpha[2]*fl[4]+alpha[1]*(93.91485505499116*fl[3]+54.22176684690384*fl[2])+(105.0*fl[1]+60.62177826491071*fl[0])*alpha[2])*dfac_x; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } else { 
  incr[0] = -0.02357022603955158*(25.98076211353316*alpha[2]*fr[7]-33.54101966249684*alpha[1]*fr[6]-15.0*alpha[2]*fr[5]-33.54101966249685*alpha[0]*fr[4]+alpha[1]*(25.98076211353316*fr[3]-15.0*fr[2])+alpha[0]*(25.98076211353316*fr[1]-15.0*fr[0]))*dfac_x; 
  incr[1] = 0.07071067811865474*(15.0*alpha[2]*fr[7]-19.36491673103708*alpha[1]*fr[6]-8.660254037844386*alpha[2]*fr[5]-19.36491673103709*alpha[0]*fr[4]+alpha[1]*(15.0*fr[3]-8.660254037844386*fr[2])+alpha[0]*(15.0*fr[1]-8.660254037844386*fr[0]))*dfac_x; 
  incr[2] = -0.02357022603955158*(23.2379000772445*alpha[1]*fr[7]+((-30.0*alpha[2])-33.54101966249684*alpha[0])*fr[6]+alpha[1]*((-13.41640786499874*fr[5])-33.54101966249685*fr[4])+(23.2379000772445*alpha[2]+25.98076211353316*alpha[0])*fr[3]+((-13.41640786499874*alpha[2])-15.0*alpha[0])*fr[2]+alpha[1]*(25.98076211353316*fr[1]-15.0*fr[0]))*dfac_x; 
  incr[3] = 0.07071067811865474*(13.41640786499874*alpha[1]*fr[7]+((-17.32050807568877*alpha[2])-19.36491673103708*alpha[0])*fr[6]+alpha[1]*((-7.745966692414834*fr[5])-19.36491673103709*fr[4])+(13.41640786499874*alpha[2]+15.0*alpha[0])*fr[3]+((-7.745966692414834*alpha[2])-8.660254037844386*alpha[0])*fr[2]+alpha[1]*(15.0*fr[1]-8.660254037844386*fr[0]))*dfac_x; 
  incr[4] = -0.1178511301977579*(11.61895003862225*alpha[2]*fr[7]-15.0*alpha[1]*fr[6]-6.708203932499369*alpha[2]*fr[5]-15.0*alpha[0]*fr[4]+alpha[1]*(11.61895003862225*fr[3]-6.708203932499369*fr[2])+alpha[0]*(11.61895003862225*fr[1]-6.708203932499369*fr[0]))*dfac_x; 
  incr[5] = -0.003367175148507369*((116.1895003862225*alpha[2]+181.8653347947321*alpha[0])*fr[7]-210.0*alpha[1]*fr[6]+((-67.0820393249937*alpha[2])-105.0*alpha[0])*fr[5]-234.787137637478*alpha[2]*fr[4]+alpha[1]*(162.6653005407115*fr[3]-93.91485505499116*fr[2])+(181.8653347947321*fr[1]-105.0*fr[0])*alpha[2])*dfac_x; 
  incr[6] = -0.02357022603955158*(51.96152422706631*alpha[1]*fr[7]+((-67.0820393249937*alpha[2])-75.0*alpha[0])*fr[6]+alpha[1]*((-30.0*fr[5])-75.00000000000001*fr[4])+(51.96152422706632*alpha[2]+58.09475019311126*alpha[0])*fr[3]+((-30.0*alpha[2])-33.54101966249684*alpha[0])*fr[2]+alpha[1]*(58.09475019311126*fr[1]-33.54101966249684*fr[0]))*dfac_x; 
  incr[7] = 0.01010152544552211*((67.0820393249937*alpha[2]+105.0*alpha[0])*fr[7]-121.2435565298214*alpha[1]*fr[6]+((-38.72983346207417*alpha[2])-60.62177826491071*alpha[0])*fr[5]-135.5544171172596*alpha[2]*fr[4]+alpha[1]*(93.91485505499116*fr[3]-54.22176684690384*fr[2])+(105.0*fr[1]-60.62177826491071*fr[0])*alpha[2])*dfac_x; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.01785714285714286*dfac_x*(3.0*(((55.0*Bmag[2]*BmagInv[2]+14.0*(2.23606797749979*BmagInv[0]*Bmag[2]+Bmag[1]*BmagInv[1]))*Phi[2]+7.0*Phi[1]*(Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2]))*geoY[2]+7.0*((2.0*(2.23606797749979*geoY[0]*Bmag[2]+Bmag[1]*geoY[1])*BmagInv[2]+(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+2.23606797749979*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*(2.0*geoY[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))*dfac_v*dfac_x*m_*wv-1.0*(24.24871130596428*(2.23606797749979*Gradpar[1]*Phi[2]+Gradpar[0]*Phi[1])*dfac_v*q_+3.0*(((55.0*Bmag[2]*BmagInv[2]+14.0*(2.23606797749979*BmagInv[0]*Bmag[2]+Bmag[1]*BmagInv[1]))*Phi[2]+7.0*Phi[1]*(Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2]))*geoY[2]+7.0*((2.0*(2.23606797749979*geoY[0]*Bmag[2]+Bmag[1]*geoY[1])*BmagInv[2]+(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+2.23606797749979*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*(2.0*geoY[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))*dfac_x*m_)))/(dfac_v*m_); 

  double alpha[3]; 
  alpha[0] = (0.05050762722761052*(dfac_v*((((Bmag[2]*(165.0*BmagInv[2]+93.91485505499116*BmagInv[0])+42.0*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*(21.0*Bmag[1]*BmagInv[2]+42.0*BmagInv[1]*Bmag[2]))*geoY[2]+((93.91485505499116*geoY[0]*Bmag[2]+42.0*Bmag[1]*geoY[1])*BmagInv[2]+(189.0*BmagInv[1]*geoY[1]+105.0*BmagInv[0]*geoY[0])*Bmag[2]+46.95742752749558*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*(42.0*geoY[1]*BmagInv[2]+46.95742752749558*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+21.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*dfac_x2*m_*wv+((-54.22176684690384*Gradpar[1]*Phi[2])-24.24871130596428*Gradpar[0]*Phi[1])*dfac_x*q_)+(((Bmag[2]*((-165.0*BmagInv[2])-93.91485505499116*BmagInv[0])-42.0*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*((-21.0*Bmag[1]*BmagInv[2])-42.0*BmagInv[1]*Bmag[2]))*geoY[2]+(((-93.91485505499116*geoY[0]*Bmag[2])-42.0*Bmag[1]*geoY[1])*BmagInv[2]+((-189.0*BmagInv[1]*geoY[1])-105.0*BmagInv[0]*geoY[0])*Bmag[2]-46.95742752749558*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*((-42.0*geoY[1]*BmagInv[2])-46.95742752749558*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))-21.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*dfac_x2*m_))/(dfac_v*m_); 
  alpha[1] = (0.01010152544552211*(dfac_v*((((368.9512162874654*Bmag[1]*BmagInv[2]+1207.476707849887*BmagInv[1]*Bmag[2]+210.0*BmagInv[0]*Bmag[1])*Phi[2]+Phi[1]*(Bmag[2]*(368.9512162874654*BmagInv[2]+210.0*BmagInv[0])+93.91485505499116*Bmag[1]*BmagInv[1]))*geoY[2]+((1207.476707849887*geoY[1]*Bmag[2]+210.0*geoY[0]*Bmag[1])*BmagInv[2]+945.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(422.6168477474603*BmagInv[1]*geoY[1]+234.787137637478*BmagInv[0]*geoY[0]))*Phi[2]+Phi[1]*((210.0*geoY[0]*Bmag[2]+93.91485505499116*Bmag[1]*geoY[1])*BmagInv[2]+(422.6168477474603*BmagInv[1]*geoY[1]+234.787137637478*BmagInv[0]*geoY[0])*Bmag[2]+105.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_x2*m_*wv+(((-242.4871130596428*Gradpar[2])-271.1088342345192*Gradpar[0])*Phi[2]-121.2435565298214*Gradpar[1]*Phi[1])*dfac_x*q_)+((((-368.9512162874654*Bmag[1]*BmagInv[2])-1207.476707849887*BmagInv[1]*Bmag[2]-210.0*BmagInv[0]*Bmag[1])*Phi[2]+Phi[1]*(Bmag[2]*((-368.9512162874654*BmagInv[2])-210.0*BmagInv[0])-93.91485505499116*Bmag[1]*BmagInv[1]))*geoY[2]+(((-1207.476707849887*geoY[1]*Bmag[2])-210.0*geoY[0]*Bmag[1])*BmagInv[2]-945.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*((-422.6168477474603*BmagInv[1]*geoY[1])-234.787137637478*BmagInv[0]*geoY[0]))*Phi[2]+Phi[1]*(((-210.0*geoY[0]*Bmag[2])-93.91485505499116*Bmag[1]*geoY[1])*BmagInv[2]+((-422.6168477474603*BmagInv[1]*geoY[1])-234.787137637478*BmagInv[0]*geoY[0])*Bmag[2]-105.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_x2*m_))/(dfac_v*m_); 
  alpha[2] = (0.01010152544552211*(dfac_v*((((Bmag[2]*(1341.640786499874*BmagInv[2]+825.0*BmagInv[0])+368.9512162874654*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*(67.0820393249937*Bmag[1]*BmagInv[2]+368.9512162874654*BmagInv[1]*Bmag[2]+105.0*BmagInv[0]*Bmag[1]))*geoY[2]+((825.0*geoY[0]*Bmag[2]+368.9512162874654*Bmag[1]*geoY[1])*BmagInv[2]+(1207.476707849887*BmagInv[1]*geoY[1]+469.5742752749559*BmagInv[0]*geoY[0])*Bmag[2]+210.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*((368.9512162874654*geoY[1]*Bmag[2]+105.0*geoY[0]*Bmag[1])*BmagInv[2]+210.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+93.91485505499116*Bmag[1]*BmagInv[1]*geoY[1]))*dfac_x2*m_*wv+((-242.4871130596428*Gradpar[1]*Phi[2])-121.2435565298214*Phi[1]*Gradpar[2])*dfac_x*q_)+(((Bmag[2]*((-1341.640786499874*BmagInv[2])-825.0*BmagInv[0])-368.9512162874654*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*((-67.0820393249937*Bmag[1]*BmagInv[2])-368.9512162874654*BmagInv[1]*Bmag[2]-105.0*BmagInv[0]*Bmag[1]))*geoY[2]+(((-825.0*geoY[0]*Bmag[2])-368.9512162874654*Bmag[1]*geoY[1])*BmagInv[2]+((-1207.476707849887*BmagInv[1]*geoY[1])-469.5742752749559*BmagInv[0]*geoY[0])*Bmag[2]-210.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(((-368.9512162874654*geoY[1]*Bmag[2])-105.0*geoY[0]*Bmag[1])*BmagInv[2]-210.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]-93.91485505499116*Bmag[1]*BmagInv[1]*geoY[1]))*dfac_x2*m_))/(dfac_v*m_); 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[3]; 
  double favg[3]; 
  favg[0] = 0.7071067811865475*(2.23606797749979*(fr[5]+fl[5])+1.732050807568877*(fl[2]-1.0*fr[2])+fr[0]+fl[0]); 
  favg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fr[7]+fl[7])+3.0*(fl[3]-1.0*fr[3]))+3.0*(fr[1]+fl[1])); 
  favg[2] = -0.1414213562373095*(8.660254037844387*fr[6]-1.0*(8.660254037844387*fl[6]+5.0*(fr[4]+fl[4]))); 

//
// //
  Ghat[0] = -0.25*((3.16227766016838*fr[5]-3.16227766016838*fl[5]-2.449489742783178*(fr[2]+fl[2])+1.414213562373095*fr[0]-1.414213562373095*fl[0])*amax-1.414213562373095*(alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0])); 
  Ghat[1] = -0.01666666666666667*((47.43416490252569*fr[7]-47.43416490252569*fl[7]-36.74234614174767*(fr[3]+fl[3])+21.21320343559643*fr[1]-21.21320343559643*fl[1])*amax-18.97366596101028*(alpha[1]*favg[2]+favg[1]*alpha[2])-21.21320343559643*(alpha[0]*favg[1]+favg[0]*alpha[1])); 
  Ghat[2] = 0.0159719141249985*((38.34057902536163*(fr[6]+fl[6])-22.13594362117866*fr[4]+22.13594362117866*fl[4])*amax+(14.14213562373095*alpha[2]+22.13594362117866*alpha[0])*favg[2]+22.13594362117866*favg[0]*alpha[2]+19.79898987322333*alpha[1]*favg[1]); 
  incr[0] = 0.7071067811865475*Ghat[0]*dfac_v; 
  incr[1] = 0.7071067811865475*Ghat[1]*dfac_v; 
  incr[2] = -1.224744871391589*Ghat[0]*dfac_v; 
  incr[3] = -1.224744871391589*Ghat[1]*dfac_v; 
  incr[4] = 0.7071067811865475*Ghat[2]*dfac_v; 
  incr[5] = 1.58113883008419*Ghat[0]*dfac_v; 
  incr[6] = -1.224744871391589*Ghat[2]*dfac_v; 
  incr[7] = 1.58113883008419*Ghat[1]*dfac_v; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
return std::abs(alpha0); 
} 
