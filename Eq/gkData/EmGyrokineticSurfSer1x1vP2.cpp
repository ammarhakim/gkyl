#include <GyrokineticModDecl.h> 
double EmGyrokineticSurf1x1vSer_X_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -0.1767766952966368*(BmagInv[0]*geoY[0]*(6.708203932499369*Apar[2]-1.732050807568877*Apar[1])*dfac_x-2.0*Gradpar[0])*wv; 

  double alpha[3]; 
  alpha[0] = -0.5*(BmagInv[0]*geoY[0]*(6.708203932499369*Apar[2]-1.732050807568877*Apar[1])*dfac_x-2.0*Gradpar[0])*wv; 
  alpha[1] = -(0.1666666666666667*(BmagInv[0]*geoY[0]*(11.61895003862225*Apar[2]-3.0*Apar[1])*dfac_x-3.464101615137754*Gradpar[0]))/dfac_v; 
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
double EmGyrokineticSurf1x1vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.125*(dfac_x*(3.0*BmagInv[0]*geoY[0]*(5.0*Apar[2]*Phi[2]+Apar[1]*Phi[1])*dfac_x+3.464101615137754*Gradpar[0]*Phi[1])+2.828427124746191*dApardt[0])*q_)/m_; 

  double alpha[3]; 
  alpha[0] = -(0.3535533905932737*(BmagInv[0]*geoY[0]*(15.0*Apar[2]*Phi[2]+3.0*Apar[1]*Phi[1])*dfac_x2+3.464101615137754*Gradpar[0]*Phi[1]*dfac_x+2.828427124746191*dApardt[0])*q_)/m_; 
  alpha[1] = -(0.3535533905932737*(6.708203932499369*BmagInv[0]*geoY[0]*(Apar[1]*Phi[2]+Phi[1]*Apar[2])*dfac_x2+7.745966692414834*Gradpar[0]*Phi[2]*dfac_x+2.828427124746191*dApardt[1])*q_)/m_; 
  alpha[2] = -(0.7071067811865475*(6.708203932499369*BmagInv[0]*geoY[0]*Apar[2]*Phi[2]*dfac_x2+1.414213562373095*dApardt[2])*q_)/m_; 
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
double EmGyrokineticSurf1x1vSer_X_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.04166666666666666*(dfac_v2*wv*(1.414213562373095*(((100.6230589874906*Bmag[2]-25.98076211353316*Bmag[1])*BmagInv[2]+5.0*(9.0*BmagInv[0]-15.58845726811989*BmagInv[1])*Bmag[2]+2.23606797749979*Bmag[1]*(9.0*BmagInv[1]-5.196152422706631*BmagInv[0]))*geoY[2]+(5.0*(9.0*geoY[0]-15.58845726811989*geoY[1])*Bmag[2]+2.23606797749979*Bmag[1]*(9.0*geoY[1]-5.196152422706631*geoY[0]))*BmagInv[2]+2.23606797749979*((27.0*BmagInv[1]-15.58845726811989*BmagInv[0])*geoY[1]+geoY[0]*(9.0*BmagInv[0]-15.58845726811989*BmagInv[1]))*Bmag[2]+Bmag[1]*((9.0*BmagInv[0]-15.58845726811989*BmagInv[1])*geoY[1]+geoY[0]*(9.0*BmagInv[1]-5.196152422706631*BmagInv[0])))*dfac_x*m_*wv+(((((225.0*Apar[2]+11.18033988749895*(9.0*Apar[0]-15.58845726811989*Apar[1]))*Bmag[2]-11.18033988749895*(5.196152422706631*Bmag[1]+12.72792206135786)*Apar[2]+5.0*((9.0*Apar[1]-5.196152422706631*Apar[0])*Bmag[1]+7.348469228349534*Apar[1]))*BmagInv[2]+(11.18033988749895*(9.0*BmagInv[0]-15.58845726811989*BmagInv[1])*Apar[2]+5.0*((27.0*Apar[1]-15.58845726811989*Apar[0])*BmagInv[1]+BmagInv[0]*(9.0*Apar[0]-15.58845726811989*Apar[1])))*Bmag[2]+5.0*((9.0*Bmag[1]+22.0454076850486)*BmagInv[1]-1.0*BmagInv[0]*(5.196152422706631*Bmag[1]+12.72792206135786))*Apar[2]+2.23606797749979*(((9.0*Apar[0]-15.58845726811989*Apar[1])*Bmag[1]-12.72792206135786*Apar[1])*BmagInv[1]+BmagInv[0]*((9.0*Apar[1]-5.196152422706631*Apar[0])*Bmag[1]+7.348469228349534*Apar[1])))*geoY[2]+((11.18033988749895*(9.0*geoY[0]-15.58845726811989*geoY[1])*Apar[2]+5.0*((27.0*Apar[1]-15.58845726811989*Apar[0])*geoY[1]+geoY[0]*(9.0*Apar[0]-15.58845726811989*Apar[1])))*Bmag[2]+5.0*((9.0*Bmag[1]+22.0454076850486)*geoY[1]-1.0*geoY[0]*(5.196152422706631*Bmag[1]+12.72792206135786))*Apar[2]+2.23606797749979*(((9.0*Apar[0]-15.58845726811989*Apar[1])*Bmag[1]-12.72792206135786*Apar[1])*geoY[1]+geoY[0]*((9.0*Apar[1]-5.196152422706631*Apar[0])*Bmag[1]+7.348469228349534*Apar[1])))*BmagInv[2]+(5.0*((27.0*BmagInv[1]-15.58845726811989*BmagInv[0])*geoY[1]+geoY[0]*(9.0*BmagInv[0]-15.58845726811989*BmagInv[1]))*Apar[2]+2.23606797749979*(((27.0*Apar[0]-46.76537180435967*Apar[1])*BmagInv[1]+BmagInv[0]*(27.0*Apar[1]-15.58845726811989*Apar[0]))*geoY[1]+geoY[0]*((27.0*Apar[1]-15.58845726811989*Apar[0])*BmagInv[1]+BmagInv[0]*(9.0*Apar[0]-15.58845726811989*Apar[1]))))*Bmag[2]-2.23606797749979*(((15.58845726811989*Bmag[1]+38.18376618407357)*BmagInv[1]-1.0*BmagInv[0]*(9.0*Bmag[1]+22.0454076850486))*geoY[1]+geoY[0]*(BmagInv[0]*(5.196152422706631*Bmag[1]+12.72792206135786)-1.0*(9.0*Bmag[1]+22.0454076850486)*BmagInv[1]))*Apar[2]+(((27.0*Apar[1]-15.58845726811989*Apar[0])*Bmag[1]+22.0454076850486*Apar[1])*BmagInv[1]+BmagInv[0]*((9.0*Apar[0]-15.58845726811989*Apar[1])*Bmag[1]-12.72792206135786*Apar[1]))*geoY[1]+geoY[0]*(((9.0*Apar[0]-15.58845726811989*Apar[1])*Bmag[1]-12.72792206135786*Apar[1])*BmagInv[1]+BmagInv[0]*((9.0*Apar[1]-5.196152422706631*Apar[0])*Bmag[1]+7.348469228349534*Apar[1])))*dfac_x+2.828427124746191*(6.708203932499369*Gradpar[2]-5.196152422706631*Gradpar[1]+3.0*Gradpar[0]))*q_)+1.414213562373095*(((33.54101966249685*Bmag[2]-8.660254037844386*Bmag[1])*BmagInv[2]+5.0*(3.0*BmagInv[0]-5.196152422706631*BmagInv[1])*Bmag[2]+2.23606797749979*Bmag[1]*(3.0*BmagInv[1]-1.732050807568877*BmagInv[0]))*geoY[2]+(5.0*(3.0*geoY[0]-5.196152422706631*geoY[1])*Bmag[2]+2.23606797749979*Bmag[1]*(3.0*geoY[1]-1.732050807568877*geoY[0]))*BmagInv[2]+2.23606797749979*((9.0*BmagInv[1]-5.196152422706631*BmagInv[0])*geoY[1]+geoY[0]*(3.0*BmagInv[0]-5.196152422706631*BmagInv[1]))*Bmag[2]+Bmag[1]*((3.0*BmagInv[0]-5.196152422706631*BmagInv[1])*geoY[1]+geoY[0]*(3.0*BmagInv[1]-1.732050807568877*BmagInv[0])))*dfac_x*m_))/(dfac_v2*q_); 

  double alpha[3]; 
  alpha[0] = (0.1178511301977579*(dfac_v2*((((142.3024947075771*Bmag[2]-36.74234614174767*Bmag[1])*BmagInv[2]+(63.63961030678928*BmagInv[0]-110.227038425243*BmagInv[1])*Bmag[2]+Bmag[1]*(28.46049894151542*BmagInv[1]-16.43167672515498*BmagInv[0]))*geoY[2]+((63.63961030678928*geoY[0]-110.227038425243*geoY[1])*Bmag[2]+Bmag[1]*(28.46049894151542*geoY[1]-16.43167672515498*geoY[0]))*BmagInv[2]+((85.38149682454625*BmagInv[1]-49.29503017546494*BmagInv[0])*geoY[1]+geoY[0]*(28.46049894151542*BmagInv[0]-49.29503017546494*BmagInv[1]))*Bmag[2]+Bmag[1]*((12.72792206135786*BmagInv[0]-22.0454076850486*BmagInv[1])*geoY[1]+geoY[0]*(12.72792206135786*BmagInv[1]-7.348469228349534*BmagInv[0])))*dfac_x*m_*wv2+(((((225.0*Apar[2]-174.2842505793337*Apar[1]+100.6230589874906*Apar[0])*Bmag[2]+((-58.09475019311126*Bmag[1])-142.3024947075771)*Apar[2]+(45.0*Apar[1]-25.98076211353316*Apar[0])*Bmag[1]+36.74234614174767*Apar[1])*BmagInv[2]+((100.6230589874906*BmagInv[0]-174.2842505793337*BmagInv[1])*Apar[2]+(135.0*Apar[1]-77.94228634059945*Apar[0])*BmagInv[1]+BmagInv[0]*(45.0*Apar[0]-77.94228634059945*Apar[1]))*Bmag[2]+((45.0*Bmag[1]+110.227038425243)*BmagInv[1]+BmagInv[0]*((-25.98076211353316*Bmag[1])-63.63961030678928))*Apar[2]+((20.12461179749811*Apar[0]-34.85685011586674*Apar[1])*Bmag[1]-28.46049894151542*Apar[1])*BmagInv[1]+BmagInv[0]*((20.12461179749811*Apar[1]-11.61895003862225*Apar[0])*Bmag[1]+16.43167672515498*Apar[1]))*geoY[2]+(((100.6230589874906*geoY[0]-174.2842505793337*geoY[1])*Apar[2]+(135.0*Apar[1]-77.94228634059945*Apar[0])*geoY[1]+geoY[0]*(45.0*Apar[0]-77.94228634059945*Apar[1]))*Bmag[2]+((45.0*Bmag[1]+110.227038425243)*geoY[1]+geoY[0]*((-25.98076211353316*Bmag[1])-63.63961030678928))*Apar[2]+((20.12461179749811*Apar[0]-34.85685011586674*Apar[1])*Bmag[1]-28.46049894151542*Apar[1])*geoY[1]+geoY[0]*((20.12461179749811*Apar[1]-11.61895003862225*Apar[0])*Bmag[1]+16.43167672515498*Apar[1]))*BmagInv[2]+(((135.0*BmagInv[1]-77.94228634059945*BmagInv[0])*geoY[1]+geoY[0]*(45.0*BmagInv[0]-77.94228634059945*BmagInv[1]))*Apar[2]+((60.37383539249433*Apar[0]-104.5705503476002*Apar[1])*BmagInv[1]+BmagInv[0]*(60.37383539249433*Apar[1]-34.85685011586674*Apar[0]))*geoY[1]+geoY[0]*((60.37383539249433*Apar[1]-34.85685011586674*Apar[0])*BmagInv[1]+BmagInv[0]*(20.12461179749811*Apar[0]-34.85685011586674*Apar[1])))*Bmag[2]+((((-34.85685011586674*Bmag[1])-85.38149682454625)*BmagInv[1]+BmagInv[0]*(20.12461179749811*Bmag[1]+49.29503017546494))*geoY[1]+geoY[0]*((20.12461179749811*Bmag[1]+49.29503017546494)*BmagInv[1]+BmagInv[0]*((-11.61895003862225*Bmag[1])-28.46049894151542)))*Apar[2]+(((27.0*Apar[1]-15.58845726811989*Apar[0])*Bmag[1]+22.0454076850486*Apar[1])*BmagInv[1]+BmagInv[0]*((9.0*Apar[0]-15.58845726811989*Apar[1])*Bmag[1]-12.72792206135786*Apar[1]))*geoY[1]+geoY[0]*(((9.0*Apar[0]-15.58845726811989*Apar[1])*Bmag[1]-12.72792206135786*Apar[1])*BmagInv[1]+BmagInv[0]*((9.0*Apar[1]-5.196152422706631*Apar[0])*Bmag[1]+7.348469228349534*Apar[1])))*dfac_x+18.97366596101028*Gradpar[2]-14.69693845669907*Gradpar[1]+8.485281374238571*Gradpar[0])*q_*wv)+(((47.43416490252571*Bmag[2]-12.24744871391589*Bmag[1])*BmagInv[2]+(21.21320343559643*BmagInv[0]-36.74234614174767*BmagInv[1])*Bmag[2]+Bmag[1]*(9.48683298050514*BmagInv[1]-5.477225575051662*BmagInv[0]))*geoY[2]+((21.21320343559643*geoY[0]-36.74234614174767*geoY[1])*Bmag[2]+Bmag[1]*(9.48683298050514*geoY[1]-5.477225575051662*geoY[0]))*BmagInv[2]+((28.46049894151542*BmagInv[1]-16.43167672515498*BmagInv[0])*geoY[1]+geoY[0]*(9.48683298050514*BmagInv[0]-16.43167672515498*BmagInv[1]))*Bmag[2]+Bmag[1]*((4.242640687119286*BmagInv[0]-7.348469228349534*BmagInv[1])*geoY[1]+geoY[0]*(4.242640687119286*BmagInv[1]-2.449489742783178*BmagInv[0])))*dfac_x*m_))/(dfac_v2*q_); 
  alpha[1] = (0.1178511301977579*((((164.3167672515499*Bmag[2]-42.42640687119286*Bmag[1])*BmagInv[2]+(73.48469228349535*BmagInv[0]-127.2792206135786*BmagInv[1])*Bmag[2]+Bmag[1]*(32.86335345030997*BmagInv[1]-18.97366596101028*BmagInv[0]))*geoY[2]+((73.48469228349535*geoY[0]-127.2792206135786*geoY[1])*Bmag[2]+Bmag[1]*(32.86335345030997*geoY[1]-18.97366596101028*geoY[0]))*BmagInv[2]+((98.5900603509299*BmagInv[1]-56.92099788303084*BmagInv[0])*geoY[1]+geoY[0]*(32.86335345030997*BmagInv[0]-56.92099788303084*BmagInv[1]))*Bmag[2]+Bmag[1]*((14.69693845669907*BmagInv[0]-25.45584412271572*BmagInv[1])*geoY[1]+geoY[0]*(14.69693845669907*BmagInv[1]-8.485281374238571*BmagInv[0])))*dfac_x*m_*wv+(((((129.9038105676658*Apar[2]-100.6230589874906*Apar[1]+58.09475019311126*Apar[0])*Bmag[2]+((-33.54101966249685*Bmag[1])-82.15838362577493)*Apar[2]+(25.98076211353316*Apar[1]-15.0*Apar[0])*Bmag[1]+21.21320343559643*Apar[1])*BmagInv[2]+((58.09475019311126*BmagInv[0]-100.6230589874906*BmagInv[1])*Apar[2]+(77.94228634059945*Apar[1]-45.0*Apar[0])*BmagInv[1]+BmagInv[0]*(25.98076211353316*Apar[0]-45.0*Apar[1]))*Bmag[2]+((25.98076211353316*Bmag[1]+63.63961030678928)*BmagInv[1]+BmagInv[0]*((-15.0*Bmag[1])-36.74234614174767))*Apar[2]+((11.61895003862225*Apar[0]-20.12461179749811*Apar[1])*Bmag[1]-16.43167672515498*Apar[1])*BmagInv[1]+BmagInv[0]*((11.61895003862225*Apar[1]-6.708203932499369*Apar[0])*Bmag[1]+9.48683298050514*Apar[1]))*geoY[2]+(((58.09475019311126*geoY[0]-100.6230589874906*geoY[1])*Apar[2]+(77.94228634059945*Apar[1]-45.0*Apar[0])*geoY[1]+geoY[0]*(25.98076211353316*Apar[0]-45.0*Apar[1]))*Bmag[2]+((25.98076211353316*Bmag[1]+63.63961030678928)*geoY[1]+geoY[0]*((-15.0*Bmag[1])-36.74234614174767))*Apar[2]+((11.61895003862225*Apar[0]-20.12461179749811*Apar[1])*Bmag[1]-16.43167672515498*Apar[1])*geoY[1]+geoY[0]*((11.61895003862225*Apar[1]-6.708203932499369*Apar[0])*Bmag[1]+9.48683298050514*Apar[1]))*BmagInv[2]+(((77.94228634059945*BmagInv[1]-45.0*BmagInv[0])*geoY[1]+geoY[0]*(25.98076211353316*BmagInv[0]-45.0*BmagInv[1]))*Apar[2]+((34.85685011586674*Apar[0]-60.37383539249433*Apar[1])*BmagInv[1]+BmagInv[0]*(34.85685011586674*Apar[1]-20.12461179749811*Apar[0]))*geoY[1]+geoY[0]*((34.85685011586674*Apar[1]-20.12461179749811*Apar[0])*BmagInv[1]+BmagInv[0]*(11.61895003862225*Apar[0]-20.12461179749811*Apar[1])))*Bmag[2]+((((-20.12461179749811*Bmag[1])-49.29503017546494)*BmagInv[1]+BmagInv[0]*(11.61895003862225*Bmag[1]+28.46049894151542))*geoY[1]+geoY[0]*((11.61895003862225*Bmag[1]+28.46049894151542)*BmagInv[1]+BmagInv[0]*((-6.708203932499369*Bmag[1])-16.43167672515498)))*Apar[2]+(((15.58845726811989*Apar[1]-9.0*Apar[0])*Bmag[1]+12.72792206135786*Apar[1])*BmagInv[1]+BmagInv[0]*((5.196152422706631*Apar[0]-9.0*Apar[1])*Bmag[1]-7.348469228349534*Apar[1]))*geoY[1]+geoY[0]*(((5.196152422706631*Apar[0]-9.0*Apar[1])*Bmag[1]-7.348469228349534*Apar[1])*BmagInv[1]+BmagInv[0]*((5.196152422706631*Apar[1]-3.0*Apar[0])*Bmag[1]+4.242640687119286*Apar[1])))*dfac_x+10.95445115010332*Gradpar[2]-8.485281374238571*Gradpar[1]+4.898979485566357*Gradpar[0])*q_))/(dfac_v*q_); 
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
double EmGyrokineticSurf1x1vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.001785714285714286*(dfac_v*(6.0*(((275.0*Bmag[2]*BmagInv[2]+14.0*(11.18033988749895*BmagInv[0]*Bmag[2]+5.0*Bmag[1]*BmagInv[1]))*Phi[2]+35.0*Phi[1]*(Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2]))*geoY[2]+7.0*((2.0*(11.18033988749895*geoY[0]*Bmag[2]+5.0*Bmag[1]*geoY[1])*BmagInv[2]+5.0*(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+11.18033988749895*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*(10.0*geoY[1]*BmagInv[2]+11.18033988749895*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+5.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))*dfac_x2*m_*wv+(3.0*((((1.414213562373095*(447.213595499958*Apar[2]+275.0*Apar[0])*Bmag[2]+11.0*(15.8113883008419*Apar[1]*Bmag[1]-50.0*Apar[2]))*BmagInv[2]+1.414213562373095*(275.0*BmagInv[0]*Apar[2]+11.18033988749895*(36.0*Apar[1]*BmagInv[1]+14.0*Apar[0]*BmagInv[0]))*Bmag[2]+11.18033988749895*(15.55634918610405*Bmag[1]*BmagInv[1]-28.0*BmagInv[0])*Apar[2]+35.0*((2.828427124746191*Apar[0]*Bmag[1]-4.0*Apar[1])*BmagInv[1]+2.828427124746191*BmagInv[0]*Apar[1]*Bmag[1]))*Phi[2]+Phi[1]*((15.8113883008419*(11.0*Apar[1]*Bmag[2]+2.0*Bmag[1]*Apar[2])+35.0*(1.414213562373095*Apar[0]*Bmag[1]-2.0*Apar[1]))*BmagInv[2]+1.414213562373095*(122.9837387624885*BmagInv[1]*Apar[2]+70.0*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Bmag[2]+7.0*(5.0*(1.414213562373095*BmagInv[0]*Bmag[1]-4.0*BmagInv[1])*Apar[2]+6.324555320336761*Apar[1]*Bmag[1]*BmagInv[1])))*geoY[2]+((1.414213562373095*(275.0*geoY[0]*Apar[2]+11.18033988749895*(36.0*Apar[1]*geoY[1]+14.0*Apar[0]*geoY[0]))*Bmag[2]+11.18033988749895*(15.55634918610405*Bmag[1]*geoY[1]-28.0*geoY[0])*Apar[2]+35.0*((2.828427124746191*Apar[0]*Bmag[1]-4.0*Apar[1])*geoY[1]+2.828427124746191*geoY[0]*Apar[1]*Bmag[1]))*BmagInv[2]+1.414213562373095*(22.3606797749979*(18.0*BmagInv[1]*geoY[1]+7.0*BmagInv[0]*geoY[0])*Apar[2]+35.0*(9.0*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1])*geoY[1]+geoY[0]*(9.0*Apar[1]*BmagInv[1]+5.0*Apar[0]*BmagInv[0])))*Bmag[2]+7.0*(5.0*((2.828427124746191*BmagInv[0]*Bmag[1]-18.0*BmagInv[1])*geoY[1]+geoY[0]*(2.828427124746191*Bmag[1]*BmagInv[1]-10.0*BmagInv[0]))*Apar[2]+2.23606797749979*((12.72792206135786*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(7.071067811865476*Apar[0]*Bmag[1]-10.0*Apar[1]))*geoY[1]+5.0*geoY[0]*((1.414213562373095*Apar[0]*Bmag[1]-2.0*Apar[1])*BmagInv[1]+1.414213562373095*BmagInv[0]*Apar[1]*Bmag[1]))))*Phi[2]+Phi[1]*((1.414213562373095*(122.9837387624885*geoY[1]*Apar[2]+70.0*(Apar[0]*geoY[1]+geoY[0]*Apar[1]))*Bmag[2]+7.0*(5.0*(1.414213562373095*geoY[0]*Bmag[1]-4.0*geoY[1])*Apar[2]+6.324555320336761*Apar[1]*Bmag[1]*geoY[1]))*BmagInv[2]+7.0*(1.414213562373095*(10.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Apar[2]+2.23606797749979*((9.0*Apar[1]*BmagInv[1]+5.0*Apar[0]*BmagInv[0])*geoY[1]+5.0*geoY[0]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1])))*Bmag[2]+2.23606797749979*((2.828427124746191*Bmag[1]*BmagInv[1]-10.0*BmagInv[0])*geoY[1]-10.0*geoY[0]*BmagInv[1])*Apar[2]+5.0*(((1.414213562373095*Apar[0]*Bmag[1]-2.0*Apar[1])*BmagInv[1]+1.414213562373095*BmagInv[0]*Apar[1]*Bmag[1])*geoY[1]+geoY[0]*(1.414213562373095*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(1.414213562373095*Apar[0]*Bmag[1]-2.0*Apar[1]))))))*dfac_x2-7.0*(6.928203230275509*(11.18033988749895*Gradpar[1]*Phi[2]+5.0*Gradpar[0]*Phi[1])*dfac_x+28.28427124746191*dApardt[0]))*q_)-6.0*(((275.0*Bmag[2]*BmagInv[2]+14.0*(11.18033988749895*BmagInv[0]*Bmag[2]+5.0*Bmag[1]*BmagInv[1]))*Phi[2]+35.0*Phi[1]*(Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2]))*geoY[2]+7.0*((2.0*(11.18033988749895*geoY[0]*Bmag[2]+5.0*Bmag[1]*geoY[1])*BmagInv[2]+5.0*(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+11.18033988749895*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*(10.0*geoY[1]*BmagInv[2]+11.18033988749895*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+5.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))))*dfac_x2*m_))/(dfac_v*m_); 

  double alpha[3]; 
  alpha[0] = (0.005050762722761052*(dfac_v*((((Bmag[2]*(1650.0*BmagInv[2]+939.1485505499119*BmagInv[0])+420.0*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*(210.0*Bmag[1]*BmagInv[2]+420.0*BmagInv[1]*Bmag[2]))*geoY[2]+((939.1485505499119*geoY[0]*Bmag[2]+420.0*Bmag[1]*geoY[1])*BmagInv[2]+(1890.0*BmagInv[1]*geoY[1]+1050.0*BmagInv[0]*geoY[0])*Bmag[2]+469.5742752749559*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*(420.0*geoY[1]*BmagInv[2]+469.5742752749559*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+210.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*dfac_x2*m_*wv+((((((1897.366596101029*Apar[2]+1166.726188957804*Apar[0])*Bmag[2]-1650.0*Apar[2]+521.7758139277827*Apar[1]*Bmag[1])*BmagInv[2]+(1166.726188957804*BmagInv[0]*Apar[2]+1707.629936490926*Apar[1]*BmagInv[1]+664.0783086353599*Apar[0]*BmagInv[0])*Bmag[2]+(521.7758139277827*Bmag[1]*BmagInv[1]-939.1485505499119*BmagInv[0])*Apar[2]+(296.98484809835*Apar[0]*Bmag[1]-420.0*Apar[1])*BmagInv[1]+296.98484809835*BmagInv[0]*Apar[1]*Bmag[1])*Phi[2]+Phi[1]*((521.7758139277827*Apar[1]*Bmag[2]+Bmag[1]*(94.86832980505142*Apar[2]+148.492424049175*Apar[0])-210.0*Apar[1])*BmagInv[2]+(521.7758139277827*BmagInv[1]*Apar[2]+296.98484809835*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Bmag[2]+(148.492424049175*BmagInv[0]*Bmag[1]-420.0*BmagInv[1])*Apar[2]+132.815661727072*Apar[1]*Bmag[1]*BmagInv[1]))*geoY[2]+(((1166.726188957804*geoY[0]*Apar[2]+1707.629936490926*Apar[1]*geoY[1]+664.0783086353599*Apar[0]*geoY[0])*Bmag[2]+(521.7758139277827*Bmag[1]*geoY[1]-939.1485505499119*geoY[0])*Apar[2]+(296.98484809835*Apar[0]*Bmag[1]-420.0*Apar[1])*geoY[1]+296.98484809835*geoY[0]*Apar[1]*Bmag[1])*BmagInv[2]+((1707.629936490926*BmagInv[1]*geoY[1]+664.0783086353599*BmagInv[0]*geoY[0])*Apar[2]+1336.431816442575*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1])*geoY[1]+geoY[0]*(1336.431816442575*Apar[1]*BmagInv[1]+742.462120245875*Apar[0]*BmagInv[0]))*Bmag[2]+((296.98484809835*BmagInv[0]*Bmag[1]-1890.0*BmagInv[1])*geoY[1]+geoY[0]*(296.98484809835*Bmag[1]*BmagInv[1]-1050.0*BmagInv[0]))*Apar[2]+(597.6704777718237*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(332.0391543176799*Apar[0]*Bmag[1]-469.5742752749559*Apar[1]))*geoY[1]+geoY[0]*((332.0391543176799*Apar[0]*Bmag[1]-469.5742752749559*Apar[1])*BmagInv[1]+332.0391543176799*BmagInv[0]*Apar[1]*Bmag[1]))*Phi[2]+Phi[1]*(((521.7758139277827*geoY[1]*Apar[2]+296.98484809835*(Apar[0]*geoY[1]+geoY[0]*Apar[1]))*Bmag[2]+(148.492424049175*geoY[0]*Bmag[1]-420.0*geoY[1])*Apar[2]+132.815661727072*Apar[1]*Bmag[1]*geoY[1])*BmagInv[2]+(296.98484809835*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Apar[2]+(597.6704777718237*Apar[1]*BmagInv[1]+332.0391543176799*Apar[0]*BmagInv[0])*geoY[1]+332.0391543176799*geoY[0]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Bmag[2]+((132.815661727072*Bmag[1]*BmagInv[1]-469.5742752749559*BmagInv[0])*geoY[1]-469.5742752749559*geoY[0]*BmagInv[1])*Apar[2]+((148.492424049175*Apar[0]*Bmag[1]-210.0*Apar[1])*BmagInv[1]+148.492424049175*BmagInv[0]*Apar[1]*Bmag[1])*geoY[1]+geoY[0]*(148.492424049175*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(148.492424049175*Apar[0]*Bmag[1]-210.0*Apar[1]))))*dfac_x2+((-542.2176684690385*Gradpar[1]*Phi[2])-242.4871130596428*Gradpar[0]*Phi[1])*dfac_x-197.9898987322334*dApardt[0])*q_)+(((Bmag[2]*((-1650.0*BmagInv[2])-939.1485505499119*BmagInv[0])-420.0*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*((-210.0*Bmag[1]*BmagInv[2])-420.0*BmagInv[1]*Bmag[2]))*geoY[2]+(((-939.1485505499119*geoY[0]*Bmag[2])-420.0*Bmag[1]*geoY[1])*BmagInv[2]+((-1890.0*BmagInv[1]*geoY[1])-1050.0*BmagInv[0]*geoY[0])*Bmag[2]-469.5742752749559*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(Bmag[2]*((-420.0*geoY[1]*BmagInv[2])-469.5742752749559*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))-210.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])))*dfac_x2*m_))/(dfac_v*m_); 
  alpha[1] = (0.005050762722761052*(dfac_v*((((737.9024325749308*Bmag[1]*BmagInv[2]+2414.953415699773*BmagInv[1]*Bmag[2]+420.0*BmagInv[0]*Bmag[1])*Phi[2]+Phi[1]*(Bmag[2]*(737.9024325749308*BmagInv[2]+420.0*BmagInv[0])+187.8297101099823*Bmag[1]*BmagInv[1]))*geoY[2]+((2414.953415699773*geoY[1]*Bmag[2]+420.0*geoY[0]*Bmag[1])*BmagInv[2]+1890.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*(845.2336954949205*BmagInv[1]*geoY[1]+469.5742752749559*BmagInv[0]*geoY[0]))*Phi[2]+Phi[1]*((420.0*geoY[0]*Bmag[2]+187.8297101099823*Bmag[1]*geoY[1])*BmagInv[2]+(845.2336954949205*BmagInv[1]*geoY[1]+469.5742752749559*BmagInv[0]*geoY[0])*Bmag[2]+210.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_x2*m_*wv+(((((2863.782463805518*Apar[1]*Bmag[2]+Bmag[1]*(848.5281374238575*Apar[2]+521.7758139277827*Apar[0])-737.9024325749308*Apar[1])*BmagInv[2]+(2863.782463805518*BmagInv[1]*Apar[2]+1707.629936490926*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Bmag[2]+(521.7758139277827*BmagInv[0]*Bmag[1]-2414.953415699773*BmagInv[1])*Apar[2]+763.6753236814716*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(296.98484809835*Apar[0]*Bmag[1]-420.0*Apar[1]))*Phi[2]+Phi[1]*(((848.5281374238575*Apar[2]+521.7758139277827*Apar[0])*Bmag[2]-737.9024325749308*Apar[2]+233.3452377915607*Apar[1]*Bmag[1])*BmagInv[2]+(521.7758139277827*BmagInv[0]*Apar[2]+763.6753236814716*Apar[1]*BmagInv[1]+296.98484809835*Apar[0]*BmagInv[0])*Bmag[2]+(233.3452377915607*Bmag[1]*BmagInv[1]-420.0*BmagInv[0])*Apar[2]+(132.815661727072*Apar[0]*Bmag[1]-187.8297101099823*Apar[1])*BmagInv[1]+132.815661727072*BmagInv[0]*Apar[1]*Bmag[1]))*geoY[2]+(((2863.782463805518*geoY[1]*Apar[2]+1707.629936490926*(Apar[0]*geoY[1]+geoY[0]*Apar[1]))*Bmag[2]+(521.7758139277827*geoY[0]*Bmag[1]-2414.953415699773*geoY[1])*Apar[2]+763.6753236814716*Apar[1]*Bmag[1]*geoY[1]+geoY[0]*(296.98484809835*Apar[0]*Bmag[1]-420.0*Apar[1]))*BmagInv[2]+(1707.629936490926*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Apar[2]+(2863.782463805518*Apar[1]*BmagInv[1]+1336.431816442575*Apar[0]*BmagInv[0])*geoY[1]+1336.431816442575*geoY[0]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Bmag[2]+((763.6753236814716*Bmag[1]*BmagInv[1]-1890.0*BmagInv[0])*geoY[1]+geoY[0]*(296.98484809835*BmagInv[0]*Bmag[1]-1890.0*BmagInv[1]))*Apar[2]+((597.6704777718237*Apar[0]*Bmag[1]-845.2336954949205*Apar[1])*BmagInv[1]+597.6704777718237*BmagInv[0]*Apar[1]*Bmag[1])*geoY[1]+geoY[0]*(597.6704777718237*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(332.0391543176799*Apar[0]*Bmag[1]-469.5742752749559*Apar[1])))*Phi[2]+Phi[1]*(((521.7758139277827*geoY[0]*Apar[2]+763.6753236814716*Apar[1]*geoY[1]+296.98484809835*Apar[0]*geoY[0])*Bmag[2]+(233.3452377915607*Bmag[1]*geoY[1]-420.0*geoY[0])*Apar[2]+(132.815661727072*Apar[0]*Bmag[1]-187.8297101099823*Apar[1])*geoY[1]+132.815661727072*geoY[0]*Apar[1]*Bmag[1])*BmagInv[2]+((763.6753236814716*BmagInv[1]*geoY[1]+296.98484809835*BmagInv[0]*geoY[0])*Apar[2]+597.6704777718237*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1])*geoY[1]+geoY[0]*(597.6704777718237*Apar[1]*BmagInv[1]+332.0391543176799*Apar[0]*BmagInv[0]))*Bmag[2]+((132.815661727072*BmagInv[0]*Bmag[1]-845.2336954949205*BmagInv[1])*geoY[1]+geoY[0]*(132.815661727072*Bmag[1]*BmagInv[1]-469.5742752749559*BmagInv[0]))*Apar[2]+(267.286363288515*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(148.492424049175*Apar[0]*Bmag[1]-210.0*Apar[1]))*geoY[1]+geoY[0]*((148.492424049175*Apar[0]*Bmag[1]-210.0*Apar[1])*BmagInv[1]+148.492424049175*BmagInv[0]*Apar[1]*Bmag[1])))*dfac_x2+(((-484.9742261192856*Gradpar[2])-542.2176684690385*Gradpar[0])*Phi[2]-242.4871130596428*Gradpar[1]*Phi[1])*dfac_x-197.9898987322334*dApardt[1])*q_)+((((-737.9024325749308*Bmag[1]*BmagInv[2])-2414.953415699773*BmagInv[1]*Bmag[2]-420.0*BmagInv[0]*Bmag[1])*Phi[2]+Phi[1]*(Bmag[2]*((-737.9024325749308*BmagInv[2])-420.0*BmagInv[0])-187.8297101099823*Bmag[1]*BmagInv[1]))*geoY[2]+(((-2414.953415699773*geoY[1]*Bmag[2])-420.0*geoY[0]*Bmag[1])*BmagInv[2]-1890.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*((-845.2336954949205*BmagInv[1]*geoY[1])-469.5742752749559*BmagInv[0]*geoY[0]))*Phi[2]+Phi[1]*(((-420.0*geoY[0]*Bmag[2])-187.8297101099823*Bmag[1]*geoY[1])*BmagInv[2]+((-845.2336954949205*BmagInv[1]*geoY[1])-469.5742752749559*BmagInv[0]*geoY[0])*Bmag[2]-210.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_x2*m_))/(dfac_v*m_); 
  alpha[2] = (4.59160247523732e-4*(dfac_v*((((Bmag[2]*(29516.09730299723*BmagInv[2]+18150.0*BmagInv[0])+8116.926758324238*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*(1475.804865149862*Bmag[1]*BmagInv[2]+8116.926758324238*BmagInv[1]*Bmag[2]+2310.0*BmagInv[0]*Bmag[1]))*geoY[2]+((18150.0*geoY[0]*Bmag[2]+8116.926758324238*Bmag[1]*geoY[1])*BmagInv[2]+(26564.48757269751*BmagInv[1]*geoY[1]+10330.63405604903*BmagInv[0]*geoY[0])*Bmag[2]+4620.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*((8116.926758324238*geoY[1]*Bmag[2]+2310.0*geoY[0]*Bmag[1])*BmagInv[2]+4620.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+2066.126811209806*Bmag[1]*BmagInv[1]*geoY[1]))*dfac_x2*m_*wv+((((((38714.09626996348*Apar[2]+20871.03255711132*Apar[0])*Bmag[2]-29516.09730299723*Apar[2]+9333.809511662432*Apar[1]*Bmag[1])*BmagInv[2]+(20871.03255711132*BmagInv[0]*Apar[2]+31501.60710186069*Apar[1]*BmagInv[1]+12833.98807853584*Apar[0]*BmagInv[0])*Bmag[2]+(9333.809511662432*Bmag[1]*BmagInv[1]-18150.0*BmagInv[0])*Apar[2]+(5739.53395320561*Apar[0]*Bmag[1]-8116.926758324238*Apar[1])*BmagInv[1]+5739.53395320561*BmagInv[0]*Apar[1]*Bmag[1])*Phi[2]+Phi[1]*((9333.809511662432*Apar[1]*Bmag[2]+Bmag[1]*(3500.178566873411*Apar[2]+1043.551627855566*Apar[0])-1475.804865149862*Apar[1])*BmagInv[2]+(9333.809511662432*BmagInv[1]*Apar[2]+5739.53395320561*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Bmag[2]+(1043.551627855566*BmagInv[0]*Bmag[1]-8116.926758324238*BmagInv[1])*Apar[2]+2566.797615707168*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(1633.416664540925*Apar[0]*Bmag[1]-2310.0*Apar[1])))*geoY[2]+(((20871.03255711132*geoY[0]*Apar[2]+31501.60710186069*Apar[1]*geoY[1]+12833.98807853584*Apar[0]*geoY[0])*Bmag[2]+(9333.809511662432*Bmag[1]*geoY[1]-18150.0*geoY[0])*Apar[2]+(5739.53395320561*Apar[0]*Bmag[1]-8116.926758324238*Apar[1])*geoY[1]+5739.53395320561*geoY[0]*Apar[1]*Bmag[1])*BmagInv[2]+((31501.60710186069*BmagInv[1]*geoY[1]+12833.98807853584*BmagInv[0]*geoY[0])*Apar[2]+18783.92930140019*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1])*geoY[1]+geoY[0]*(18783.92930140019*Apar[1]*BmagInv[1]+7304.86139498896*Apar[0]*BmagInv[0]))*Bmag[2]+((5739.53395320561*BmagInv[0]*Bmag[1]-26564.48757269751*BmagInv[1])*geoY[1]+geoY[0]*(5739.53395320561*Bmag[1]*BmagInv[1]-10330.63405604903*BmagInv[0]))*Apar[2]+(8400.428560496188*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(3266.83332908185*Apar[0]*Bmag[1]-4620.0*Apar[1]))*geoY[1]+geoY[0]*((3266.83332908185*Apar[0]*Bmag[1]-4620.0*Apar[1])*BmagInv[1]+3266.83332908185*BmagInv[0]*Apar[1]*Bmag[1]))*Phi[2]+Phi[1]*(((9333.809511662432*geoY[1]*Apar[2]+5739.53395320561*(Apar[0]*geoY[1]+geoY[0]*Apar[1]))*Bmag[2]+(1043.551627855566*geoY[0]*Bmag[1]-8116.926758324238*geoY[1])*Apar[2]+2566.797615707168*Apar[1]*Bmag[1]*geoY[1]+geoY[0]*(1633.416664540925*Apar[0]*Bmag[1]-2310.0*Apar[1]))*BmagInv[2]+(5739.53395320561*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Apar[2]+(8400.428560496188*Apar[1]*BmagInv[1]+3266.83332908185*Apar[0]*BmagInv[0])*geoY[1]+3266.83332908185*geoY[0]*(Apar[0]*BmagInv[1]+BmagInv[0]*Apar[1]))*Bmag[2]+((2566.797615707168*Bmag[1]*BmagInv[1]-4620.0*BmagInv[0])*geoY[1]+geoY[0]*(1633.416664540925*BmagInv[0]*Bmag[1]-4620.0*BmagInv[1]))*Apar[2]+((1460.972278997792*Apar[0]*Bmag[1]-2066.126811209806*Apar[1])*BmagInv[1]+1460.972278997792*BmagInv[0]*Apar[1]*Bmag[1])*geoY[1]+1460.972278997792*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]))*dfac_x2+((-5334.716487312142*Gradpar[1]*Phi[2])-2667.358243656071*Phi[1]*Gradpar[2])*dfac_x-2177.888886054567*dApardt[2])*q_)+(((Bmag[2]*((-29516.09730299723*BmagInv[2])-18150.0*BmagInv[0])-8116.926758324238*Bmag[1]*BmagInv[1])*Phi[2]+Phi[1]*((-1475.804865149862*Bmag[1]*BmagInv[2])-8116.926758324238*BmagInv[1]*Bmag[2]-2310.0*BmagInv[0]*Bmag[1]))*geoY[2]+(((-18150.0*geoY[0]*Bmag[2])-8116.926758324238*Bmag[1]*geoY[1])*BmagInv[2]+((-26564.48757269751*BmagInv[1]*geoY[1])-10330.63405604903*BmagInv[0]*geoY[0])*Bmag[2]-4620.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*Phi[2]+Phi[1]*(((-8116.926758324238*geoY[1]*Bmag[2])-2310.0*geoY[0]*Bmag[1])*BmagInv[2]-4620.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]-2066.126811209806*Bmag[1]*BmagInv[1]*geoY[1]))*dfac_x2*m_))/(dfac_v*m_); 
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
