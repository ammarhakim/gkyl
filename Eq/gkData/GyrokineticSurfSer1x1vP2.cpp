#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x1vSer_X_P2_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[3]; 
  alpha[0] = Gradpar[0]*wv; 
  alpha[1] = (0.5773502691896258*Gradpar[0])/dfac_v; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1178511301977579*(6.708203932499369*alpha[1]*fl[6]+6.708203932499369*alpha[0]*fl[4]+alpha[1]*(5.196152422706631*fl[3]+3.0*fl[2])+alpha[0]*(5.196152422706631*fl[1]+3.0*fl[0]))*dfac_x; 
  incr[1] = -0.3535533905932737*(3.872983346207417*(alpha[1]*fl[6]+alpha[0]*fl[4])+alpha[1]*(3.0*fl[3]+1.732050807568877*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.02357022603955158*(23.2379000772445*alpha[1]*fl[7]+33.54101966249684*alpha[0]*fl[6]+alpha[1]*(13.41640786499874*fl[5]+33.54101966249685*fl[4])+alpha[0]*(25.98076211353316*fl[3]+15.0*fl[2])+alpha[1]*(25.98076211353316*fl[1]+15.0*fl[0]))*dfac_x; 
  incr[3] = -0.07071067811865474*(13.41640786499874*alpha[1]*fl[7]+19.36491673103708*alpha[0]*fl[6]+alpha[1]*(7.745966692414834*fl[5]+19.36491673103709*fl[4])+alpha[0]*(15.0*fl[3]+8.660254037844386*fl[2])+alpha[1]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[4] = 0.1178511301977579*(15.0*alpha[1]*fl[6]+15.0*alpha[0]*fl[4]+alpha[1]*(11.61895003862225*fl[3]+6.708203932499369*fl[2])+alpha[0]*(11.61895003862225*fl[1]+6.708203932499369*fl[0]))*dfac_x; 
  incr[5] = 0.02357022603955158*(25.98076211353316*alpha[0]*fl[7]+30.0*alpha[1]*fl[6]+15.0*alpha[0]*fl[5]+alpha[1]*(23.2379000772445*fl[3]+13.41640786499874*fl[2]))*dfac_x; 
  incr[6] = 0.02357022603955158*(51.96152422706631*alpha[1]*fl[7]+75.0*alpha[0]*fl[6]+alpha[1]*(30.0*fl[5]+75.00000000000001*fl[4])+alpha[0]*(58.09475019311126*fl[3]+33.54101966249684*fl[2])+alpha[1]*(58.09475019311126*fl[1]+33.54101966249684*fl[0]))*dfac_x; 
  incr[7] = -0.07071067811865474*(15.0*alpha[0]*fl[7]+17.32050807568877*alpha[1]*fl[6]+8.660254037844387*alpha[0]*fl[5]+alpha[1]*(13.41640786499874*fl[3]+7.745966692414834*fl[2]))*dfac_x; 
  } else { 
  incr[0] = 0.1178511301977579*(6.708203932499369*alpha[1]*fr[6]+6.708203932499369*alpha[0]*fr[4]+alpha[1]*(3.0*fr[2]-5.196152422706631*fr[3])+alpha[0]*(3.0*fr[0]-5.196152422706631*fr[1]))*dfac_x; 
  incr[1] = -0.3535533905932737*(3.872983346207417*(alpha[1]*fr[6]+alpha[0]*fr[4])+alpha[1]*(1.732050807568877*fr[2]-3.0*fr[3])+alpha[0]*(1.732050807568877*fr[0]-3.0*fr[1]))*dfac_x; 
  incr[2] = -0.02357022603955158*(23.2379000772445*alpha[1]*fr[7]-33.54101966249684*alpha[0]*fr[6]+alpha[1]*((-13.41640786499874*fr[5])-33.54101966249685*fr[4])+alpha[0]*(25.98076211353316*fr[3]-15.0*fr[2])+alpha[1]*(25.98076211353316*fr[1]-15.0*fr[0]))*dfac_x; 
  incr[3] = 0.07071067811865474*(13.41640786499874*alpha[1]*fr[7]-19.36491673103708*alpha[0]*fr[6]+alpha[1]*((-7.745966692414834*fr[5])-19.36491673103709*fr[4])+alpha[0]*(15.0*fr[3]-8.660254037844386*fr[2])+alpha[1]*(15.0*fr[1]-8.660254037844386*fr[0]))*dfac_x; 
  incr[4] = 0.1178511301977579*(15.0*alpha[1]*fr[6]+15.0*alpha[0]*fr[4]+alpha[1]*(6.708203932499369*fr[2]-11.61895003862225*fr[3])+alpha[0]*(6.708203932499369*fr[0]-11.61895003862225*fr[1]))*dfac_x; 
  incr[5] = -0.02357022603955158*(25.98076211353316*alpha[0]*fr[7]-30.0*alpha[1]*fr[6]-15.0*alpha[0]*fr[5]+alpha[1]*(23.2379000772445*fr[3]-13.41640786499874*fr[2]))*dfac_x; 
  incr[6] = -0.02357022603955158*(51.96152422706631*alpha[1]*fr[7]-75.0*alpha[0]*fr[6]+alpha[1]*((-30.0*fr[5])-75.00000000000001*fr[4])+alpha[0]*(58.09475019311126*fr[3]-33.54101966249684*fr[2])+alpha[1]*(58.09475019311126*fr[1]-33.54101966249684*fr[0]))*dfac_x; 
  incr[7] = 0.07071067811865474*(15.0*alpha[0]*fr[7]-17.32050807568877*alpha[1]*fr[6]-8.660254037844387*alpha[0]*fr[5]+alpha[1]*(13.41640786499874*fr[3]-7.745966692414834*fr[2]))*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[3];
double fupwindQuad[3];
double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]-0.9486832980505137*alpha[1]; 
  fupwindQuad[0] = 0.5*((0.7745966692414833*(fr[7]+fl[7])+1.5*fr[6]-1.5*fl[6]-0.4472135954999579*fr[5]+0.4472135954999579*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]-1.161895003862225*(fr[3]+fl[3])+0.6708203932499369*fr[2]-0.6708203932499369*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.7745966692414833*fr[7]+0.7745966692414833*fl[7]-1.5*(fr[6]+fl[6])+0.4472135954999579*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])+1.161895003862225*fr[3]-1.161895003862225*fl[3]-0.6708203932499369*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-0.9682458365518543*(fr[7]+fl[7]))+0.5590169943749475*fr[5]-0.5590169943749475*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)+0.9682458365518543*fr[7]-0.9682458365518543*fl[7]-0.5590169943749475*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0]; 
  fupwindQuad[2] = 0.5*((0.7745966692414833*(fr[7]+fl[7])-1.5*fr[6]+1.5*fl[6]-0.4472135954999579*fr[5]+0.4472135954999579*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]+1.161895003862225*(fr[3]+fl[3])-0.6708203932499369*fr[2]+0.6708203932499369*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.7745966692414833*fr[7]+0.7745966692414833*fl[7]+1.5*(fr[6]+fl[6])+0.4472135954999579*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])-1.161895003862225*fr[3]+1.161895003862225*fl[3]+0.6708203932499369*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  fupwind[0] = 0.3928371006591931*fupwindQuad[2]+0.6285393610547092*fupwindQuad[1]+0.3928371006591931*fupwindQuad[0]; 
  fupwind[1] = 0.52704627669473*fupwindQuad[2]-0.52704627669473*fupwindQuad[0]; 
  fupwind[2] = 0.3513641844631533*fupwindQuad[2]-0.7027283689263066*fupwindQuad[1]+0.3513641844631533*fupwindQuad[0]; 
  incr[0] = 0.5*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_x; 
  incr[1] = -0.8660254037844386*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_x; 
  incr[2] = 0.1*(4.47213595499958*alpha[1]*fupwind[2]+5.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_x; 
  incr[3] = -0.1732050807568877*(4.47213595499958*alpha[1]*fupwind[2]+5.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_x; 
  incr[4] = 1.118033988749895*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_x; 
  incr[5] = 0.1*(5.0*alpha[0]*fupwind[2]+4.47213595499958*alpha[1]*fupwind[1])*dfac_x; 
  incr[6] = 0.223606797749979*(4.47213595499958*alpha[1]*fupwind[2]+5.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_x; 
  incr[7] = -0.3872983346207417*(2.23606797749979*alpha[0]*fupwind[2]+2.0*alpha[1]*fupwind[1])*dfac_x; 

#endif 
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
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
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
double GyrokineticSurf1x1vSer_X_P2_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*(2.23606797749979*Gradpar[2]-1.732050807568877*Gradpar[1]+Gradpar[0])*wv; 

  double alpha[3]; 
  alpha[0] = (2.23606797749979*Gradpar[2]-1.732050807568877*Gradpar[1]+Gradpar[0])*wv; 
  alpha[1] = (0.3333333333333333*(3.872983346207417*Gradpar[2]-3.0*Gradpar[1]+1.732050807568877*Gradpar[0]))/dfac_v; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  incr[0] = 0.1178511301977579*(6.708203932499369*alpha[1]*fl[6]+6.708203932499369*alpha[0]*fl[4]+alpha[1]*(5.196152422706631*fl[3]+3.0*fl[2])+alpha[0]*(5.196152422706631*fl[1]+3.0*fl[0]))*dfac_x; 
  incr[1] = -0.3535533905932737*(3.872983346207417*(alpha[1]*fl[6]+alpha[0]*fl[4])+alpha[1]*(3.0*fl[3]+1.732050807568877*fl[2])+alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0]))*dfac_x; 
  incr[2] = 0.02357022603955158*(23.2379000772445*alpha[1]*fl[7]+33.54101966249684*alpha[0]*fl[6]+alpha[1]*(13.41640786499874*fl[5]+33.54101966249685*fl[4])+alpha[0]*(25.98076211353316*fl[3]+15.0*fl[2])+alpha[1]*(25.98076211353316*fl[1]+15.0*fl[0]))*dfac_x; 
  incr[3] = -0.07071067811865474*(13.41640786499874*alpha[1]*fl[7]+19.36491673103708*alpha[0]*fl[6]+alpha[1]*(7.745966692414834*fl[5]+19.36491673103709*fl[4])+alpha[0]*(15.0*fl[3]+8.660254037844386*fl[2])+alpha[1]*(15.0*fl[1]+8.660254037844386*fl[0]))*dfac_x; 
  incr[4] = 0.1178511301977579*(15.0*alpha[1]*fl[6]+15.0*alpha[0]*fl[4]+alpha[1]*(11.61895003862225*fl[3]+6.708203932499369*fl[2])+alpha[0]*(11.61895003862225*fl[1]+6.708203932499369*fl[0]))*dfac_x; 
  incr[5] = 0.02357022603955158*(25.98076211353316*alpha[0]*fl[7]+30.0*alpha[1]*fl[6]+15.0*alpha[0]*fl[5]+alpha[1]*(23.2379000772445*fl[3]+13.41640786499874*fl[2]))*dfac_x; 
  incr[6] = 0.02357022603955158*(51.96152422706631*alpha[1]*fl[7]+75.0*alpha[0]*fl[6]+alpha[1]*(30.0*fl[5]+75.00000000000001*fl[4])+alpha[0]*(58.09475019311126*fl[3]+33.54101966249684*fl[2])+alpha[1]*(58.09475019311126*fl[1]+33.54101966249684*fl[0]))*dfac_x; 
  incr[7] = -0.07071067811865474*(15.0*alpha[0]*fl[7]+17.32050807568877*alpha[1]*fl[6]+8.660254037844387*alpha[0]*fl[5]+alpha[1]*(13.41640786499874*fl[3]+7.745966692414834*fl[2]))*dfac_x; 
  } else { 
  incr[0] = 0.1178511301977579*(6.708203932499369*alpha[1]*fr[6]+6.708203932499369*alpha[0]*fr[4]+alpha[1]*(3.0*fr[2]-5.196152422706631*fr[3])+alpha[0]*(3.0*fr[0]-5.196152422706631*fr[1]))*dfac_x; 
  incr[1] = -0.3535533905932737*(3.872983346207417*(alpha[1]*fr[6]+alpha[0]*fr[4])+alpha[1]*(1.732050807568877*fr[2]-3.0*fr[3])+alpha[0]*(1.732050807568877*fr[0]-3.0*fr[1]))*dfac_x; 
  incr[2] = -0.02357022603955158*(23.2379000772445*alpha[1]*fr[7]-33.54101966249684*alpha[0]*fr[6]+alpha[1]*((-13.41640786499874*fr[5])-33.54101966249685*fr[4])+alpha[0]*(25.98076211353316*fr[3]-15.0*fr[2])+alpha[1]*(25.98076211353316*fr[1]-15.0*fr[0]))*dfac_x; 
  incr[3] = 0.07071067811865474*(13.41640786499874*alpha[1]*fr[7]-19.36491673103708*alpha[0]*fr[6]+alpha[1]*((-7.745966692414834*fr[5])-19.36491673103709*fr[4])+alpha[0]*(15.0*fr[3]-8.660254037844386*fr[2])+alpha[1]*(15.0*fr[1]-8.660254037844386*fr[0]))*dfac_x; 
  incr[4] = 0.1178511301977579*(15.0*alpha[1]*fr[6]+15.0*alpha[0]*fr[4]+alpha[1]*(6.708203932499369*fr[2]-11.61895003862225*fr[3])+alpha[0]*(6.708203932499369*fr[0]-11.61895003862225*fr[1]))*dfac_x; 
  incr[5] = -0.02357022603955158*(25.98076211353316*alpha[0]*fr[7]-30.0*alpha[1]*fr[6]-15.0*alpha[0]*fr[5]+alpha[1]*(23.2379000772445*fr[3]-13.41640786499874*fr[2]))*dfac_x; 
  incr[6] = -0.02357022603955158*(51.96152422706631*alpha[1]*fr[7]-75.0*alpha[0]*fr[6]+alpha[1]*((-30.0*fr[5])-75.00000000000001*fr[4])+alpha[0]*(58.09475019311126*fr[3]-33.54101966249684*fr[2])+alpha[1]*(58.09475019311126*fr[1]-33.54101966249684*fr[0]))*dfac_x; 
  incr[7] = 0.07071067811865474*(15.0*alpha[0]*fr[7]-17.32050807568877*alpha[1]*fr[6]-8.660254037844387*alpha[0]*fr[5]+alpha[1]*(13.41640786499874*fr[3]-7.745966692414834*fr[2]))*dfac_x; 
  }
#elif upwindType == QUAD 
double fupwind[3];
double fupwindQuad[3];
double alphaQuad;
  alphaQuad = 0.7071067811865475*alpha[0]-0.9486832980505137*alpha[1]; 
  fupwindQuad[0] = 0.5*((0.7745966692414833*(fr[7]+fl[7])+1.5*fr[6]-1.5*fl[6]-0.4472135954999579*fr[5]+0.4472135954999579*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]-1.161895003862225*(fr[3]+fl[3])+0.6708203932499369*fr[2]-0.6708203932499369*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.7745966692414833*fr[7]+0.7745966692414833*fl[7]-1.5*(fr[6]+fl[6])+0.4472135954999579*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])+1.161895003862225*fr[3]-1.161895003862225*fl[3]-0.6708203932499369*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  fupwindQuad[1] = 0.5*(((-0.9682458365518543*(fr[7]+fl[7]))+0.5590169943749475*fr[5]-0.5590169943749475*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)+0.9682458365518543*fr[7]-0.9682458365518543*fl[7]-0.5590169943749475*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  alphaQuad = 0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0]; 
  fupwindQuad[2] = 0.5*((0.7745966692414833*(fr[7]+fl[7])-1.5*fr[6]+1.5*fl[6]-0.4472135954999579*fr[5]+0.4472135954999579*fl[5]-1.118033988749895*fr[4]+1.118033988749895*fl[4]+1.161895003862225*(fr[3]+fl[3])-0.6708203932499369*fr[2]+0.6708203932499369*fl[2]+0.8660254037844386*(fr[1]+fl[1])-0.5*fr[0]+0.5*fl[0])*sgn(alphaQuad)-0.7745966692414833*fr[7]+0.7745966692414833*fl[7]+1.5*(fr[6]+fl[6])+0.4472135954999579*(fr[5]+fl[5])+1.118033988749895*(fr[4]+fl[4])-1.161895003862225*fr[3]+1.161895003862225*fl[3]+0.6708203932499369*(fr[2]+fl[2])-0.8660254037844386*fr[1]+0.8660254037844386*fl[1]+0.5*(fr[0]+fl[0])); 
  fupwind[0] = 0.3928371006591931*fupwindQuad[2]+0.6285393610547092*fupwindQuad[1]+0.3928371006591931*fupwindQuad[0]; 
  fupwind[1] = 0.52704627669473*fupwindQuad[2]-0.52704627669473*fupwindQuad[0]; 
  fupwind[2] = 0.3513641844631533*fupwindQuad[2]-0.7027283689263066*fupwindQuad[1]+0.3513641844631533*fupwindQuad[0]; 
  incr[0] = 0.5*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_x; 
  incr[1] = -0.8660254037844386*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_x; 
  incr[2] = 0.1*(4.47213595499958*alpha[1]*fupwind[2]+5.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_x; 
  incr[3] = -0.1732050807568877*(4.47213595499958*alpha[1]*fupwind[2]+5.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_x; 
  incr[4] = 1.118033988749895*(alpha[1]*fupwind[1]+alpha[0]*fupwind[0])*dfac_x; 
  incr[5] = 0.1*(5.0*alpha[0]*fupwind[2]+4.47213595499958*alpha[1]*fupwind[1])*dfac_x; 
  incr[6] = 0.223606797749979*(4.47213595499958*alpha[1]*fupwind[2]+5.0*(alpha[0]*fupwind[1]+fupwind[0]*alpha[1]))*dfac_x; 
  incr[7] = -0.3872983346207417*(2.23606797749979*alpha[0]*fupwind[2]+2.0*alpha[1]*fupwind[1])*dfac_x; 

#endif 
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
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*(2.23606797749979*Gradpar[1]*Phi[2]+Gradpar[0]*Phi[1])*dfac_x*q_)/m_; 

  double alpha[3]; 
  alpha[0] = -(0.7071067811865475*(3.872983346207417*Gradpar[1]*Phi[2]+1.732050807568877*Gradpar[0]*Phi[1])*dfac_x*q_)/m_; 
  alpha[1] = -(0.7071067811865475*((3.464101615137754*Gradpar[2]+3.872983346207417*Gradpar[0])*Phi[2]+1.732050807568877*Gradpar[1]*Phi[1])*dfac_x*q_)/m_; 
  alpha[2] = -(0.7071067811865475*(3.464101615137754*Gradpar[1]*Phi[2]+1.732050807568877*Phi[1]*Gradpar[2])*dfac_x*q_)/m_; 
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
