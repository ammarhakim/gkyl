#include <GyrokineticModDecl.h>
double EmGyrokineticSimpleHelicalSurf1x2vSer_x_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[20]; 
  hamilR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(bmag[0]*wmuR+phi[0]*q_))+2.0*m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamilR[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[7] = 2.0*phi[2]*q_; 
  hamilR[8] = (0.421637021355784*m_)/rdvpar2SqR; 

  double BstarZdBmagR[20]; 
  BstarZdBmagR[0] = 2.0*cmag[0]; 

  double alphaR[8]; 
  alphaR[0] = (0.4330127018922193*BstarZdBmagR[0]*hamilR[2]*rdvpar2R)/m_; 
  alphaR[1] = (0.9682458365518543*BstarZdBmagR[0]*hamilR[8]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.1082531754730548*BstarZdBmagR[0]*hamilR[2]*rdvpar2R)/m_; 

  double incr[20]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.08333333333333333*(6.708203932499369*alphaR[1]*fL[11]+6.708203932499369*alphaR[0]*fL[7]+alphaR[1]*(5.196152422706631*fL[4]+3.0*fL[2])+alphaR[0]*(5.196152422706631*fL[1]+3.0*fL[0])); 
  incr[1] = -0.25*(3.872983346207417*(alphaR[1]*fL[11]+alphaR[0]*fL[7])+alphaR[1]*(3.0*fL[4]+1.732050807568877*fL[2])+alphaR[0]*(3.0*fL[1]+1.732050807568877*fL[0])); 
  incr[2] = 0.01666666666666667*(23.2379000772445*alphaR[1]*fL[12]+33.54101966249684*alphaR[0]*fL[11]+alphaR[1]*(13.41640786499874*fL[8]+33.54101966249685*fL[7])+alphaR[0]*(25.98076211353316*fL[4]+15.0*fL[2])+alphaR[1]*(25.98076211353316*fL[1]+15.0*fL[0])); 
  incr[3] = 0.08333333333333333*(6.708203932499369*alphaR[1]*fL[17]+6.708203932499369*alphaR[0]*fL[13]+alphaR[1]*(5.196152422706631*fL[10]+3.0*fL[6])+alphaR[0]*(5.196152422706631*fL[5]+3.0*fL[3])); 
  incr[4] = -0.05*(13.41640786499874*alphaR[1]*fL[12]+19.36491673103708*alphaR[0]*fL[11]+alphaR[1]*(7.745966692414834*fL[8]+19.36491673103709*fL[7])+alphaR[0]*(15.0*fL[4]+8.660254037844386*fL[2])+alphaR[1]*(15.0*fL[1]+8.660254037844386*fL[0])); 
  incr[5] = -0.25*(3.872983346207417*(alphaR[1]*fL[17]+alphaR[0]*fL[13])+alphaR[1]*(3.0*fL[10]+1.732050807568877*fL[6])+alphaR[0]*(3.0*fL[5]+1.732050807568877*fL[3])); 
  incr[6] = 0.01666666666666667*(23.2379000772445*alphaR[1]*fL[18]+33.54101966249685*alphaR[0]*fL[17]+alphaR[1]*(13.41640786499874*fL[14]+33.54101966249684*fL[13])+alphaR[0]*(25.98076211353316*fL[10]+15.0*fL[6])+alphaR[1]*(25.98076211353316*fL[5]+15.0*fL[3])); 
  incr[7] = 0.08333333333333333*(15.0*alphaR[1]*fL[11]+15.0*alphaR[0]*fL[7]+alphaR[1]*(11.61895003862225*fL[4]+6.708203932499369*fL[2])+alphaR[0]*(11.61895003862225*fL[1]+6.708203932499369*fL[0])); 
  incr[8] = 0.01666666666666667*(25.98076211353316*alphaR[0]*fL[12]+30.0*alphaR[1]*fL[11]+15.0*alphaR[0]*fL[8]+alphaR[1]*(23.2379000772445*fL[4]+13.41640786499874*fL[2])); 
  incr[9] = 0.01666666666666667*(alphaR[1]*(25.98076211353316*fL[19]+15.0*fL[16])+alphaR[0]*(25.98076211353316*fL[15]+15.0*fL[9])); 
  incr[10] = -0.05*(13.41640786499874*alphaR[1]*fL[18]+19.36491673103709*alphaR[0]*fL[17]+alphaR[1]*(7.745966692414834*fL[14]+19.36491673103708*fL[13])+alphaR[0]*(15.0*fL[10]+8.660254037844386*fL[6])+alphaR[1]*(15.0*fL[5]+8.660254037844386*fL[3])); 
  incr[11] = 0.01666666666666667*(51.96152422706631*alphaR[1]*fL[12]+75.0*alphaR[0]*fL[11]+alphaR[1]*(30.0*fL[8]+75.00000000000001*fL[7])+alphaR[0]*(58.09475019311126*fL[4]+33.54101966249684*fL[2])+alphaR[1]*(58.09475019311126*fL[1]+33.54101966249684*fL[0])); 
  incr[12] = -0.05*(15.0*alphaR[0]*fL[12]+17.32050807568877*alphaR[1]*fL[11]+8.660254037844387*alphaR[0]*fL[8]+alphaR[1]*(13.41640786499874*fL[4]+7.745966692414834*fL[2])); 
  incr[13] = 0.08333333333333333*(15.0*alphaR[1]*fL[17]+15.0*alphaR[0]*fL[13]+alphaR[1]*(11.61895003862225*fL[10]+6.708203932499369*fL[6])+alphaR[0]*(11.61895003862225*fL[5]+6.708203932499369*fL[3])); 
  incr[14] = 0.01666666666666667*(25.98076211353316*alphaR[0]*fL[18]+30.0*alphaR[1]*fL[17]+15.0*alphaR[0]*fL[14]+alphaR[1]*(23.2379000772445*fL[10]+13.41640786499874*fL[6])); 
  incr[15] = -0.05*(alphaR[1]*(15.0*fL[19]+8.660254037844386*fL[16])+alphaR[0]*(15.0*fL[15]+8.660254037844387*fL[9])); 
  incr[16] = 0.01666666666666667*(alphaR[0]*(25.98076211353316*fL[19]+15.0*fL[16])+alphaR[1]*(25.98076211353316*fL[15]+15.0*fL[9])); 
  incr[17] = 0.01666666666666667*(51.96152422706631*alphaR[1]*fL[18]+75.0*alphaR[0]*fL[17]+alphaR[1]*(30.0*fL[14]+75.00000000000001*fL[13])+alphaR[0]*(58.09475019311126*fL[10]+33.54101966249685*fL[6])+alphaR[1]*(58.09475019311126*fL[5]+33.54101966249685*fL[3])); 
  incr[18] = -0.05*(15.0*alphaR[0]*fL[18]+17.32050807568877*alphaR[1]*fL[17]+8.660254037844387*alphaR[0]*fL[14]+alphaR[1]*(13.41640786499874*fL[10]+7.745966692414834*fL[6])); 
  incr[19] = -0.05*(alphaR[0]*(15.0*fL[19]+8.660254037844387*fL[16])+alphaR[1]*(15.0*fL[15]+8.660254037844386*fL[9])); 
  } else { 
  incr[0] = 0.08333333333333333*(6.708203932499369*alphaR[1]*fR[11]+6.708203932499369*alphaR[0]*fR[7]+alphaR[1]*(3.0*fR[2]-5.196152422706631*fR[4])+alphaR[0]*(3.0*fR[0]-5.196152422706631*fR[1])); 
  incr[1] = -0.25*(3.872983346207417*(alphaR[1]*fR[11]+alphaR[0]*fR[7])+alphaR[1]*(1.732050807568877*fR[2]-3.0*fR[4])+alphaR[0]*(1.732050807568877*fR[0]-3.0*fR[1])); 
  incr[2] = -0.01666666666666667*(23.2379000772445*alphaR[1]*fR[12]-33.54101966249684*alphaR[0]*fR[11]+alphaR[1]*((-13.41640786499874*fR[8])-33.54101966249685*fR[7])+alphaR[0]*(25.98076211353316*fR[4]-15.0*fR[2])+alphaR[1]*(25.98076211353316*fR[1]-15.0*fR[0])); 
  incr[3] = 0.08333333333333333*(6.708203932499369*alphaR[1]*fR[17]+6.708203932499369*alphaR[0]*fR[13]+alphaR[1]*(3.0*fR[6]-5.196152422706631*fR[10])+alphaR[0]*(3.0*fR[3]-5.196152422706631*fR[5])); 
  incr[4] = 0.05*(13.41640786499874*alphaR[1]*fR[12]-19.36491673103708*alphaR[0]*fR[11]+alphaR[1]*((-7.745966692414834*fR[8])-19.36491673103709*fR[7])+alphaR[0]*(15.0*fR[4]-8.660254037844386*fR[2])+alphaR[1]*(15.0*fR[1]-8.660254037844386*fR[0])); 
  incr[5] = -0.25*(3.872983346207417*(alphaR[1]*fR[17]+alphaR[0]*fR[13])+alphaR[1]*(1.732050807568877*fR[6]-3.0*fR[10])+alphaR[0]*(1.732050807568877*fR[3]-3.0*fR[5])); 
  incr[6] = -0.01666666666666667*(23.2379000772445*alphaR[1]*fR[18]-33.54101966249685*alphaR[0]*fR[17]+alphaR[1]*((-13.41640786499874*fR[14])-33.54101966249684*fR[13])+alphaR[0]*(25.98076211353316*fR[10]-15.0*fR[6])+alphaR[1]*(25.98076211353316*fR[5]-15.0*fR[3])); 
  incr[7] = 0.08333333333333333*(15.0*alphaR[1]*fR[11]+15.0*alphaR[0]*fR[7]+alphaR[1]*(6.708203932499369*fR[2]-11.61895003862225*fR[4])+alphaR[0]*(6.708203932499369*fR[0]-11.61895003862225*fR[1])); 
  incr[8] = -0.01666666666666667*(25.98076211353316*alphaR[0]*fR[12]-30.0*alphaR[1]*fR[11]-15.0*alphaR[0]*fR[8]+alphaR[1]*(23.2379000772445*fR[4]-13.41640786499874*fR[2])); 
  incr[9] = -0.01666666666666667*(alphaR[1]*(25.98076211353316*fR[19]-15.0*fR[16])+alphaR[0]*(25.98076211353316*fR[15]-15.0*fR[9])); 
  incr[10] = 0.05*(13.41640786499874*alphaR[1]*fR[18]-19.36491673103709*alphaR[0]*fR[17]+alphaR[1]*((-7.745966692414834*fR[14])-19.36491673103708*fR[13])+alphaR[0]*(15.0*fR[10]-8.660254037844386*fR[6])+alphaR[1]*(15.0*fR[5]-8.660254037844386*fR[3])); 
  incr[11] = -0.01666666666666667*(51.96152422706631*alphaR[1]*fR[12]-75.0*alphaR[0]*fR[11]+alphaR[1]*((-30.0*fR[8])-75.00000000000001*fR[7])+alphaR[0]*(58.09475019311126*fR[4]-33.54101966249684*fR[2])+alphaR[1]*(58.09475019311126*fR[1]-33.54101966249684*fR[0])); 
  incr[12] = 0.05*(15.0*alphaR[0]*fR[12]-17.32050807568877*alphaR[1]*fR[11]-8.660254037844387*alphaR[0]*fR[8]+alphaR[1]*(13.41640786499874*fR[4]-7.745966692414834*fR[2])); 
  incr[13] = 0.08333333333333333*(15.0*alphaR[1]*fR[17]+15.0*alphaR[0]*fR[13]+alphaR[1]*(6.708203932499369*fR[6]-11.61895003862225*fR[10])+alphaR[0]*(6.708203932499369*fR[3]-11.61895003862225*fR[5])); 
  incr[14] = -0.01666666666666667*(25.98076211353316*alphaR[0]*fR[18]-30.0*alphaR[1]*fR[17]-15.0*alphaR[0]*fR[14]+alphaR[1]*(23.2379000772445*fR[10]-13.41640786499874*fR[6])); 
  incr[15] = 0.05*(alphaR[1]*(15.0*fR[19]-8.660254037844386*fR[16])+alphaR[0]*(15.0*fR[15]-8.660254037844387*fR[9])); 
  incr[16] = -0.01666666666666667*(alphaR[0]*(25.98076211353316*fR[19]-15.0*fR[16])+alphaR[1]*(25.98076211353316*fR[15]-15.0*fR[9])); 
  incr[17] = -0.01666666666666667*(51.96152422706631*alphaR[1]*fR[18]-75.0*alphaR[0]*fR[17]+alphaR[1]*((-30.0*fR[14])-75.00000000000001*fR[13])+alphaR[0]*(58.09475019311126*fR[10]-33.54101966249685*fR[6])+alphaR[1]*(58.09475019311126*fR[5]-33.54101966249685*fR[3])); 
  incr[18] = 0.05*(15.0*alphaR[0]*fR[18]-17.32050807568877*alphaR[1]*fR[17]-8.660254037844387*alphaR[0]*fR[14]+alphaR[1]*(13.41640786499874*fR[10]-7.745966692414834*fR[6])); 
  incr[19] = 0.05*(alphaR[0]*(15.0*fR[19]-8.660254037844387*fR[16])+alphaR[1]*(15.0*fR[15]-8.660254037844386*fR[9])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.5*alphaR[0]-0.6708203932499369*alphaR[1]; 
  fUpOrd[0] = 0.5*(((-0.7348469228349525*(fR[19]+fL[19]))-0.7348469228349533*(fR[18]+fL[18])-1.42302494707577*fR[17]+1.42302494707577*fL[17]+0.4242640687119281*fR[16]-0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])+0.4242640687119285*fR[14]-0.4242640687119285*fL[14]+1.060660171779821*fR[13]-1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])+1.060660171779821*fR[11]-1.060660171779821*fL[11]+1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6363961030678926*fR[6]+0.6363961030678926*fL[6]-0.8215838362577489*(fR[5]+fL[5]+fR[4]+fL[4])+0.4743416490252568*(fR[3]+fR[2])-0.4743416490252568*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.7348469228349525*fR[19]-0.7348469228349525*fL[19]+0.7348469228349533*fR[18]-0.7348469228349533*fL[18]+1.42302494707577*(fR[17]+fL[17])-0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]-0.4242640687119285*(fR[14]+fL[14])-1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]-1.060660171779821*(fR[11]+fL[11])-1.10227038425243*fR[10]+1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])+0.6363961030678926*(fR[6]+fL[6])+0.8215838362577489*(fR[5]+fR[4])-0.8215838362577489*(fL[5]+fL[4])-0.4743416490252568*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[1] = 0.5*((0.6846531968814573*(fR[18]+fL[18])-0.3952847075210473*fR[14]+0.3952847075210473*fL[14]+0.7905694150420948*fR[13]-0.7905694150420948*fL[13]-0.6846531968814574*(fR[12]+fL[12])+0.3952847075210473*fR[8]-0.3952847075210473*fL[8]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6123724356957944*(fR[5]+fL[5])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6846531968814573*fR[18]+0.6846531968814573*fL[18]+0.3952847075210473*(fR[14]+fL[14])-0.7905694150420948*(fR[13]+fL[13])+0.6846531968814574*fR[12]-0.6846531968814574*fL[12]-0.3952847075210473*(fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])+0.6123724356957944*fR[5]-0.6123724356957944*fL[5]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[2] = 0.5*((0.7348469228349525*(fR[19]+fL[19])-0.7348469228349533*(fR[18]+fL[18])+1.42302494707577*fR[17]-1.42302494707577*fL[17]-0.4242640687119281*fR[16]+0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])+0.4242640687119285*fR[14]-0.4242640687119285*fL[14]+1.060660171779821*fR[13]-1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])-1.060660171779821*fR[11]+1.060660171779821*fL[11]-1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6363961030678926*fR[6]-0.6363961030678926*fL[6]-0.8215838362577489*(fR[5]+fL[5])+0.8215838362577489*(fR[4]+fL[4])+0.4743416490252568*fR[3]-0.4743416490252568*(fL[3]+fR[2])+0.4743416490252568*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.7348469228349525*fR[19]+0.7348469228349525*fL[19]+0.7348469228349533*fR[18]-0.7348469228349533*fL[18]-1.42302494707577*(fR[17]+fL[17])+0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]-0.4242640687119285*(fR[14]+fL[14])-1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]+1.060660171779821*(fR[11]+fL[11])+1.10227038425243*fR[10]-1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6363961030678926*(fR[6]+fL[6])+0.8215838362577489*fR[5]-0.8215838362577489*(fL[5]+fR[4])+0.8215838362577489*fL[4]-0.4743416490252568*(fR[3]+fL[3])+0.4743416490252568*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.5*alphaR[1]; 
  fUpOrd[3] = 0.5*((0.6846531968814573*(fR[19]+fL[19])-0.3952847075210473*fR[16]+0.3952847075210473*fL[16]-0.6846531968814574*(fR[15]+fL[15])+0.7905694150420948*fR[11]-0.7905694150420948*fL[11]+0.3952847075210473*fR[9]-0.3952847075210473*fL[9]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[2]-0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6846531968814573*fR[19]+0.6846531968814573*fL[19]+0.3952847075210473*(fR[16]+fL[16])+0.6846531968814574*fR[15]-0.6846531968814574*fL[15]-0.7905694150420948*(fR[11]+fL[11])-0.3952847075210473*(fR[9]+fL[9])+0.7905694150420947*(fR[7]+fL[7])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]-0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[1]+alphaR[0]); 
  fUpOrd[4] = 0.5*(((-0.6846531968814573*(fR[19]+fL[19]))+0.3952847075210473*fR[16]-0.3952847075210473*fL[16]-0.6846531968814574*(fR[15]+fL[15])-0.7905694150420948*fR[11]+0.7905694150420948*fL[11]+0.3952847075210473*fR[9]-0.3952847075210473*fL[9]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[2]+0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6846531968814573*fR[19]-0.6846531968814573*fL[19]-0.3952847075210473*(fR[16]+fL[16])+0.6846531968814574*fR[15]-0.6846531968814574*fL[15]+0.7905694150420948*(fR[11]+fL[11])-0.3952847075210473*(fR[9]+fL[9])+0.7905694150420947*(fR[7]+fL[7])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]+0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.6708203932499369*alphaR[1]; 
  fUpOrd[5] = 0.5*(((-0.7348469228349525*(fR[19]+fL[19]))+0.7348469228349533*(fR[18]+fL[18])+1.42302494707577*fR[17]-1.42302494707577*fL[17]+0.4242640687119281*fR[16]-0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])-0.4242640687119285*fR[14]+0.4242640687119285*fL[14]-1.060660171779821*fR[13]+1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])+1.060660171779821*fR[11]-1.060660171779821*fL[11]-1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6363961030678926*fR[6]-0.6363961030678926*fL[6]+0.8215838362577489*(fR[5]+fL[5])-0.8215838362577489*(fR[4]+fL[4])-0.4743416490252568*fR[3]+0.4743416490252568*(fL[3]+fR[2])-0.4743416490252568*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.7348469228349525*fR[19]-0.7348469228349525*fL[19]-0.7348469228349533*fR[18]+0.7348469228349533*fL[18]-1.42302494707577*(fR[17]+fL[17])-0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]+0.4242640687119285*(fR[14]+fL[14])+1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]-1.060660171779821*(fR[11]+fL[11])+1.10227038425243*fR[10]-1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6363961030678926*(fR[6]+fL[6])-0.8215838362577489*fR[5]+0.8215838362577489*(fL[5]+fR[4])-0.8215838362577489*fL[4]+0.4743416490252568*(fR[3]+fL[3])-0.4743416490252568*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[6] = 0.5*(((-0.6846531968814573*(fR[18]+fL[18]))+0.3952847075210473*fR[14]-0.3952847075210473*fL[14]-0.7905694150420948*fR[13]+0.7905694150420948*fL[13]-0.6846531968814574*(fR[12]+fL[12])+0.3952847075210473*fR[8]-0.3952847075210473*fL[8]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6123724356957944*(fR[5]+fL[5])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6846531968814573*fR[18]-0.6846531968814573*fL[18]-0.3952847075210473*(fR[14]+fL[14])+0.7905694150420948*(fR[13]+fL[13])+0.6846531968814574*fR[12]-0.6846531968814574*fL[12]-0.3952847075210473*(fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6123724356957944*fR[5]+0.6123724356957944*fL[5]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[7] = 0.5*((0.7348469228349525*(fR[19]+fL[19])+0.7348469228349533*(fR[18]+fL[18])-1.42302494707577*fR[17]+1.42302494707577*fL[17]-0.4242640687119281*fR[16]+0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])-0.4242640687119285*fR[14]+0.4242640687119285*fL[14]-1.060660171779821*fR[13]+1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])-1.060660171779821*fR[11]+1.060660171779821*fL[11]+1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6363961030678926*fR[6]+0.6363961030678926*fL[6]+0.8215838362577489*(fR[5]+fL[5]+fR[4]+fL[4])-0.4743416490252568*(fR[3]+fR[2])+0.4743416490252568*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.7348469228349525*fR[19]+0.7348469228349525*fL[19]-0.7348469228349533*fR[18]+0.7348469228349533*fL[18]+1.42302494707577*(fR[17]+fL[17])+0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]+0.4242640687119285*(fR[14]+fL[14])+1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]+1.060660171779821*(fR[11]+fL[11])-1.10227038425243*fR[10]+1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])+0.6363961030678926*(fR[6]+fL[6])-0.8215838362577489*(fR[5]+fR[4])+0.8215838362577489*(fL[5]+fL[4])+0.4743416490252568*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.03846153846153846*(5.0*fUpOrd[7]+8.0*fUpOrd[6]+5.0*fUpOrd[5]+8.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]+8.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[1] = 0.00828173324999922*(25.0*(fUpOrd[7]-1.0*fUpOrd[5])+53.66563145999496*(fUpOrd[4]-1.0*fUpOrd[3])+25.0*(fUpOrd[2]-1.0*fUpOrd[0])); 
  fUp[2] = 0.01851851851851852*(11.18033988749895*fUpOrd[7]+24.0*fUpOrd[6]+11.18033988749895*fUpOrd[5]-1.0*(11.18033988749895*fUpOrd[2]+24.0*fUpOrd[1]+11.18033988749895*fUpOrd[0])); 
  fUp[3] = 0.2777777777777778*(fUpOrd[7]-1.0*(fUpOrd[5]+fUpOrd[2])+fUpOrd[0]); 
  fUp[4] = 0.03440104580768907*(5.0*fUpOrd[7]-18.0*fUpOrd[6]+5.0*fUpOrd[5]+8.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]-18.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[5] = 0.03440104580768907*(5.0*fUpOrd[7]+8.0*fUpOrd[6]+5.0*fUpOrd[5]-18.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]+8.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[6] = 0.08281733249999224*(2.23606797749979*fUpOrd[7]-6.0*fUpOrd[6]+2.23606797749979*(fUpOrd[5]-1.0*fUpOrd[2])+6.0*fUpOrd[1]-2.23606797749979*fUpOrd[0]); 
  fUp[7] = 0.01656346649999844*(11.18033988749895*(fUpOrd[7]-1.0*fUpOrd[5])+10.0*(3.0*fUpOrd[3]-3.0*fUpOrd[4])+11.18033988749895*(fUpOrd[2]-1.0*fUpOrd[0])); 

  incr[0] = 0.3535533905932737*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = -0.6123724356957944*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[2] = 0.07071067811865474*(4.47213595499958*alphaR[1]*fUp[4]+5.0*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[3] = 0.3535533905932737*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2]); 
  incr[4] = -0.07071067811865474*(7.745966692414834*alphaR[1]*fUp[4]+8.660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[5] = -0.6123724356957944*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2]); 
  incr[6] = 0.02357022603955158*(13.41640786499874*alphaR[1]*fUp[6]+15.0*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[7] = 0.7905694150420947*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[8] = 0.07071067811865474*(5.0*alphaR[0]*fUp[4]+4.47213595499958*alphaR[1]*fUp[1]); 
  incr[9] = 0.02357022603955158*(15.0*alphaR[1]*fUp[7]+15.0*alphaR[0]*fUp[5]); 
  incr[10] = -0.07071067811865474*(7.745966692414834*alphaR[1]*fUp[6]+8.660254037844386*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[11] = 0.02357022603955158*(30.0*alphaR[1]*fUp[4]+33.54101966249684*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[12] = -0.07071067811865474*(8.660254037844387*alphaR[0]*fUp[4]+7.745966692414834*alphaR[1]*fUp[1]); 
  incr[13] = 0.7905694150420945*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2]); 
  incr[14] = 0.02357022603955158*(15.0*alphaR[0]*fUp[6]+13.41640786499874*alphaR[1]*fUp[3]); 
  incr[15] = -0.07071067811865474*(8.660254037844386*alphaR[1]*fUp[7]+8.660254037844387*alphaR[0]*fUp[5]); 
  incr[16] = 0.02357022603955158*(15.0*alphaR[0]*fUp[7]+15.0*alphaR[1]*fUp[5]); 
  incr[17] = 0.02357022603955158*(30.0*alphaR[1]*fUp[6]+33.54101966249685*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[18] = -0.07071067811865474*(8.660254037844387*alphaR[0]*fUp[6]+7.745966692414834*alphaR[1]*fUp[3]); 
  incr[19] = -0.07071067811865474*(8.660254037844387*alphaR[0]*fUp[7]+8.660254037844386*alphaR[1]*fUp[5]); 

#endif 
  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 
  outR[8] += incr[8]*rdx2R; 
  outR[9] += incr[9]*rdx2R; 
  outR[10] += incr[10]*rdx2R; 
  outR[11] += incr[11]*rdx2R; 
  outR[12] += incr[12]*rdx2R; 
  outR[13] += incr[13]*rdx2R; 
  outR[14] += incr[14]*rdx2R; 
  outR[15] += incr[15]*rdx2R; 
  outR[16] += incr[16]*rdx2R; 
  outR[17] += incr[17]*rdx2R; 
  outR[18] += incr[18]*rdx2R; 
  outR[19] += incr[19]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += -1.0*incr[7]*rdx2L; 
  outL[8] += -1.0*incr[8]*rdx2L; 
  outL[9] += -1.0*incr[9]*rdx2L; 
  outL[10] += incr[10]*rdx2L; 
  outL[11] += -1.0*incr[11]*rdx2L; 
  outL[12] += incr[12]*rdx2L; 
  outL[13] += -1.0*incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 
  outL[16] += -1.0*incr[16]*rdx2L; 
  outL[17] += -1.0*incr[17]*rdx2L; 
  outL[18] += incr[18]*rdx2L; 
  outL[19] += incr[19]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmGyrokineticSimpleHelicalSurf1x2vSer_vpar_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[20]; 
  hamilR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(bmag[0]*wmuR+phi[0]*q_))+2.0*m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamilR[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[7] = 2.0*phi[2]*q_; 
  hamilR[8] = (0.421637021355784*m_)/rdvpar2SqR; 

  double BstarZdBmagR[20]; 
  BstarZdBmagR[0] = 2.0*cmag[0]; 

  double alphaR[8]; 
  alphaR[0] = -(0.25*(1.732050807568877*BstarZdBmagR[0]*hamilR[1]*rdx2R+5.656854249492382*dApardt[0]*q_))/m_; 
  alphaR[1] = -(0.25*(3.872983346207417*BstarZdBmagR[0]*hamilR[7]*rdx2R+5.656854249492382*dApardt[1]*q_))/m_; 
  alphaR[4] = -(1.414213562373095*dApardt[2]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(1.732050807568877*BstarZdBmagR[0]*hamilR[1]*rdx2R+5.656854249492382*dApardt[0]*q_))/m_; 

  double incr[20]; 
  double amax = amax_in; 

  double fAvg[8]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[8]+fL[8])+1.732050807568877*(fL[2]-1.0*fR[2])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[12]+fL[12])+3.0*(fL[4]-1.0*fR[4]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[14]+fL[14])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.7071067811865475*(2.23606797749979*(fR[18]+fL[18])+1.732050807568877*(fL[10]-1.0*fR[10])+fR[5]+fL[5]); 
  fAvg[4] = -0.1414213562373095*(8.660254037844387*fR[11]-1.0*(8.660254037844387*fL[11]+5.0*(fR[7]+fL[7]))); 
  fAvg[5] = -0.1414213562373095*(8.660254037844387*fR[16]-1.0*(8.660254037844387*fL[16]+5.0*(fR[9]+fL[9]))); 
  fAvg[6] = -0.1414213562373095*(8.660254037844387*fR[17]-1.0*(8.660254037844387*fL[17]+5.0*(fR[13]+fL[13]))); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[19]-1.0*(8.660254037844387*fL[19]+5.0*(fR[15]+fL[15]))); 

  double Ghat[8]; 
  Ghat[0] = -0.1767766952966368*((4.47213595499958*fR[8]-4.47213595499958*fL[8]-3.464101615137754*(fR[2]+fL[2])+2.0*fR[0]-2.0*fL[0])*amax-1.414213562373095*(alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.01178511301977579*((67.08203932499369*fR[12]-67.08203932499369*fL[12]-51.96152422706631*(fR[4]+fL[4])+30.0*fR[1]-30.0*fL[1])*amax-18.97366596101028*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])-21.21320343559643*(alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.01178511301977579*((67.08203932499369*fR[14]-67.08203932499369*fL[14]-51.96152422706631*(fR[6]+fL[6])+30.0*fR[3]-30.0*fL[3])*amax-21.21320343559643*alphaR[4]*fAvg[6]-21.21320343559643*(alphaR[1]*fAvg[3]+alphaR[0]*fAvg[2])); 
  Ghat[3] = -0.01178511301977579*((67.0820393249937*fR[18]-67.0820393249937*fL[18]-51.96152422706631*(fR[10]+fL[10])+30.0*fR[5]-30.0*fL[5])*amax-18.97366596101028*(alphaR[1]*fAvg[6]+fAvg[3]*alphaR[4])-21.21320343559643*(alphaR[0]*fAvg[3]+alphaR[1]*fAvg[2])); 
  Ghat[4] = 0.005050762722761052*((121.2435565298214*(fR[11]+fL[11])-70.0*fR[7]+70.0*fL[7])*amax+(31.62277660168381*alphaR[4]+49.49747468305833*alphaR[0])*fAvg[4]+49.49747468305833*fAvg[0]*alphaR[4]+44.27188724235732*alphaR[1]*fAvg[1]); 
  Ghat[5] = 0.01178511301977579*((51.96152422706632*(fR[16]+fL[16])-30.0*fR[9]+30.0*fL[9])*amax+21.21320343559643*alphaR[1]*fAvg[7]+21.21320343559643*alphaR[0]*fAvg[5]); 
  Ghat[6] = 0.001683587574253684*((363.7306695894642*(fR[17]+fL[17])-210.0*fR[13]+210.0*fL[13])*amax+(94.86832980505142*alphaR[4]+148.492424049175*alphaR[0])*fAvg[6]+148.492424049175*fAvg[2]*alphaR[4]+132.815661727072*alphaR[1]*fAvg[3]); 
  Ghat[7] = 0.01178511301977579*((51.96152422706632*(fR[19]+fL[19])-30.0*fR[15]+30.0*fL[15])*amax+(18.97366596101028*alphaR[4]+21.21320343559643*alphaR[0])*fAvg[7]+21.21320343559643*alphaR[1]*fAvg[5]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = -1.224744871391589*Ghat[0]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = -1.224744871391589*Ghat[1]; 
  incr[5] = 0.7071067811865475*Ghat[3]; 
  incr[6] = -1.224744871391589*Ghat[2]; 
  incr[7] = 0.7071067811865475*Ghat[4]; 
  incr[8] = 1.58113883008419*Ghat[0]; 
  incr[9] = 0.7071067811865475*Ghat[5]; 
  incr[10] = -1.224744871391589*Ghat[3]; 
  incr[11] = -1.224744871391589*Ghat[4]; 
  incr[12] = 1.58113883008419*Ghat[1]; 
  incr[13] = 0.7071067811865475*Ghat[6]; 
  incr[14] = 1.58113883008419*Ghat[2]; 
  incr[15] = 0.7071067811865475*Ghat[7]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = -1.224744871391589*Ghat[6]; 
  incr[18] = 1.58113883008419*Ghat[3]; 
  incr[19] = -1.224744871391589*Ghat[7]; 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 
  outR[8] += incr[8]*rdvpar2R; 
  outR[9] += incr[9]*rdvpar2R; 
  outR[10] += incr[10]*rdvpar2R; 
  outR[11] += incr[11]*rdvpar2R; 
  outR[12] += incr[12]*rdvpar2R; 
  outR[13] += incr[13]*rdvpar2R; 
  outR[14] += incr[14]*rdvpar2R; 
  outR[15] += incr[15]*rdvpar2R; 
  outR[16] += incr[16]*rdvpar2R; 
  outR[17] += incr[17]*rdvpar2R; 
  outR[18] += incr[18]*rdvpar2R; 
  outR[19] += incr[19]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += -1.0*incr[7]*rdvpar2L; 
  outL[8] += -1.0*incr[8]*rdvpar2L; 
  outL[9] += -1.0*incr[9]*rdvpar2L; 
  outL[10] += incr[10]*rdvpar2L; 
  outL[11] += incr[11]*rdvpar2L; 
  outL[12] += -1.0*incr[12]*rdvpar2L; 
  outL[13] += -1.0*incr[13]*rdvpar2L; 
  outL[14] += -1.0*incr[14]*rdvpar2L; 
  outL[15] += -1.0*incr[15]*rdvpar2L; 
  outL[16] += incr[16]*rdvpar2L; 
  outL[17] += incr[17]*rdvpar2L; 
  outL[18] += -1.0*incr[18]*rdvpar2L; 
  outL[19] += incr[19]*rdvpar2L; 
return std::abs(alphaSurfAvgR); 
} 
double EmGyrokineticSimpleHelicalSurf1x2vSer_x_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[20]; 
  hamilR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(bmag[0]*wmuR+phi[0]*q_))+2.0*m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamilR[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[7] = 2.0*(bmag[2]*wmuR+phi[2]*q_); 
  hamilR[8] = (0.421637021355784*m_)/rdvpar2SqR; 
  hamilR[13] = (1.154700538379251*bmag[2])/rdmu2R; 

  double BstarZdBmagR[20]; 
  BstarZdBmagR[0] = 2.0*cmag[0]; 
  BstarZdBmagR[1] = 2.0*cmag[1]; 
  BstarZdBmagR[7] = 2.0*cmag[2]; 

  double alphaR[8]; 
  alphaR[0] = (0.25*hamilR[2]*(3.872983346207417*BstarZdBmagR[7]-3.0*BstarZdBmagR[1]+1.732050807568877*BstarZdBmagR[0])*rdvpar2R)/m_; 
  alphaR[1] = (0.25*(8.660254037844386*BstarZdBmagR[7]-6.708203932499369*BstarZdBmagR[1]+3.872983346207417*BstarZdBmagR[0])*hamilR[8]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(3.872983346207417*hamilR[2]*BstarZdBmagR[7]+(1.732050807568877*BstarZdBmagR[0]-3.0*BstarZdBmagR[1])*hamilR[2])*rdvpar2R)/m_; 

  double incr[20]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.08333333333333333*(6.708203932499369*alphaR[1]*fL[11]+6.708203932499369*alphaR[0]*fL[7]+alphaR[1]*(5.196152422706631*fL[4]+3.0*fL[2])+alphaR[0]*(5.196152422706631*fL[1]+3.0*fL[0])); 
  incr[1] = -0.25*(3.872983346207417*(alphaR[1]*fL[11]+alphaR[0]*fL[7])+alphaR[1]*(3.0*fL[4]+1.732050807568877*fL[2])+alphaR[0]*(3.0*fL[1]+1.732050807568877*fL[0])); 
  incr[2] = 0.01666666666666667*(23.2379000772445*alphaR[1]*fL[12]+33.54101966249684*alphaR[0]*fL[11]+alphaR[1]*(13.41640786499874*fL[8]+33.54101966249685*fL[7])+alphaR[0]*(25.98076211353316*fL[4]+15.0*fL[2])+alphaR[1]*(25.98076211353316*fL[1]+15.0*fL[0])); 
  incr[3] = 0.08333333333333333*(6.708203932499369*alphaR[1]*fL[17]+6.708203932499369*alphaR[0]*fL[13]+alphaR[1]*(5.196152422706631*fL[10]+3.0*fL[6])+alphaR[0]*(5.196152422706631*fL[5]+3.0*fL[3])); 
  incr[4] = -0.05*(13.41640786499874*alphaR[1]*fL[12]+19.36491673103708*alphaR[0]*fL[11]+alphaR[1]*(7.745966692414834*fL[8]+19.36491673103709*fL[7])+alphaR[0]*(15.0*fL[4]+8.660254037844386*fL[2])+alphaR[1]*(15.0*fL[1]+8.660254037844386*fL[0])); 
  incr[5] = -0.25*(3.872983346207417*(alphaR[1]*fL[17]+alphaR[0]*fL[13])+alphaR[1]*(3.0*fL[10]+1.732050807568877*fL[6])+alphaR[0]*(3.0*fL[5]+1.732050807568877*fL[3])); 
  incr[6] = 0.01666666666666667*(23.2379000772445*alphaR[1]*fL[18]+33.54101966249685*alphaR[0]*fL[17]+alphaR[1]*(13.41640786499874*fL[14]+33.54101966249684*fL[13])+alphaR[0]*(25.98076211353316*fL[10]+15.0*fL[6])+alphaR[1]*(25.98076211353316*fL[5]+15.0*fL[3])); 
  incr[7] = 0.08333333333333333*(15.0*alphaR[1]*fL[11]+15.0*alphaR[0]*fL[7]+alphaR[1]*(11.61895003862225*fL[4]+6.708203932499369*fL[2])+alphaR[0]*(11.61895003862225*fL[1]+6.708203932499369*fL[0])); 
  incr[8] = 0.01666666666666667*(25.98076211353316*alphaR[0]*fL[12]+30.0*alphaR[1]*fL[11]+15.0*alphaR[0]*fL[8]+alphaR[1]*(23.2379000772445*fL[4]+13.41640786499874*fL[2])); 
  incr[9] = 0.01666666666666667*(alphaR[1]*(25.98076211353316*fL[19]+15.0*fL[16])+alphaR[0]*(25.98076211353316*fL[15]+15.0*fL[9])); 
  incr[10] = -0.05*(13.41640786499874*alphaR[1]*fL[18]+19.36491673103709*alphaR[0]*fL[17]+alphaR[1]*(7.745966692414834*fL[14]+19.36491673103708*fL[13])+alphaR[0]*(15.0*fL[10]+8.660254037844386*fL[6])+alphaR[1]*(15.0*fL[5]+8.660254037844386*fL[3])); 
  incr[11] = 0.01666666666666667*(51.96152422706631*alphaR[1]*fL[12]+75.0*alphaR[0]*fL[11]+alphaR[1]*(30.0*fL[8]+75.00000000000001*fL[7])+alphaR[0]*(58.09475019311126*fL[4]+33.54101966249684*fL[2])+alphaR[1]*(58.09475019311126*fL[1]+33.54101966249684*fL[0])); 
  incr[12] = -0.05*(15.0*alphaR[0]*fL[12]+17.32050807568877*alphaR[1]*fL[11]+8.660254037844387*alphaR[0]*fL[8]+alphaR[1]*(13.41640786499874*fL[4]+7.745966692414834*fL[2])); 
  incr[13] = 0.08333333333333333*(15.0*alphaR[1]*fL[17]+15.0*alphaR[0]*fL[13]+alphaR[1]*(11.61895003862225*fL[10]+6.708203932499369*fL[6])+alphaR[0]*(11.61895003862225*fL[5]+6.708203932499369*fL[3])); 
  incr[14] = 0.01666666666666667*(25.98076211353316*alphaR[0]*fL[18]+30.0*alphaR[1]*fL[17]+15.0*alphaR[0]*fL[14]+alphaR[1]*(23.2379000772445*fL[10]+13.41640786499874*fL[6])); 
  incr[15] = -0.05*(alphaR[1]*(15.0*fL[19]+8.660254037844386*fL[16])+alphaR[0]*(15.0*fL[15]+8.660254037844387*fL[9])); 
  incr[16] = 0.01666666666666667*(alphaR[0]*(25.98076211353316*fL[19]+15.0*fL[16])+alphaR[1]*(25.98076211353316*fL[15]+15.0*fL[9])); 
  incr[17] = 0.01666666666666667*(51.96152422706631*alphaR[1]*fL[18]+75.0*alphaR[0]*fL[17]+alphaR[1]*(30.0*fL[14]+75.00000000000001*fL[13])+alphaR[0]*(58.09475019311126*fL[10]+33.54101966249685*fL[6])+alphaR[1]*(58.09475019311126*fL[5]+33.54101966249685*fL[3])); 
  incr[18] = -0.05*(15.0*alphaR[0]*fL[18]+17.32050807568877*alphaR[1]*fL[17]+8.660254037844387*alphaR[0]*fL[14]+alphaR[1]*(13.41640786499874*fL[10]+7.745966692414834*fL[6])); 
  incr[19] = -0.05*(alphaR[0]*(15.0*fL[19]+8.660254037844387*fL[16])+alphaR[1]*(15.0*fL[15]+8.660254037844386*fL[9])); 
  } else { 
  incr[0] = 0.08333333333333333*(6.708203932499369*alphaR[1]*fR[11]+6.708203932499369*alphaR[0]*fR[7]+alphaR[1]*(3.0*fR[2]-5.196152422706631*fR[4])+alphaR[0]*(3.0*fR[0]-5.196152422706631*fR[1])); 
  incr[1] = -0.25*(3.872983346207417*(alphaR[1]*fR[11]+alphaR[0]*fR[7])+alphaR[1]*(1.732050807568877*fR[2]-3.0*fR[4])+alphaR[0]*(1.732050807568877*fR[0]-3.0*fR[1])); 
  incr[2] = -0.01666666666666667*(23.2379000772445*alphaR[1]*fR[12]-33.54101966249684*alphaR[0]*fR[11]+alphaR[1]*((-13.41640786499874*fR[8])-33.54101966249685*fR[7])+alphaR[0]*(25.98076211353316*fR[4]-15.0*fR[2])+alphaR[1]*(25.98076211353316*fR[1]-15.0*fR[0])); 
  incr[3] = 0.08333333333333333*(6.708203932499369*alphaR[1]*fR[17]+6.708203932499369*alphaR[0]*fR[13]+alphaR[1]*(3.0*fR[6]-5.196152422706631*fR[10])+alphaR[0]*(3.0*fR[3]-5.196152422706631*fR[5])); 
  incr[4] = 0.05*(13.41640786499874*alphaR[1]*fR[12]-19.36491673103708*alphaR[0]*fR[11]+alphaR[1]*((-7.745966692414834*fR[8])-19.36491673103709*fR[7])+alphaR[0]*(15.0*fR[4]-8.660254037844386*fR[2])+alphaR[1]*(15.0*fR[1]-8.660254037844386*fR[0])); 
  incr[5] = -0.25*(3.872983346207417*(alphaR[1]*fR[17]+alphaR[0]*fR[13])+alphaR[1]*(1.732050807568877*fR[6]-3.0*fR[10])+alphaR[0]*(1.732050807568877*fR[3]-3.0*fR[5])); 
  incr[6] = -0.01666666666666667*(23.2379000772445*alphaR[1]*fR[18]-33.54101966249685*alphaR[0]*fR[17]+alphaR[1]*((-13.41640786499874*fR[14])-33.54101966249684*fR[13])+alphaR[0]*(25.98076211353316*fR[10]-15.0*fR[6])+alphaR[1]*(25.98076211353316*fR[5]-15.0*fR[3])); 
  incr[7] = 0.08333333333333333*(15.0*alphaR[1]*fR[11]+15.0*alphaR[0]*fR[7]+alphaR[1]*(6.708203932499369*fR[2]-11.61895003862225*fR[4])+alphaR[0]*(6.708203932499369*fR[0]-11.61895003862225*fR[1])); 
  incr[8] = -0.01666666666666667*(25.98076211353316*alphaR[0]*fR[12]-30.0*alphaR[1]*fR[11]-15.0*alphaR[0]*fR[8]+alphaR[1]*(23.2379000772445*fR[4]-13.41640786499874*fR[2])); 
  incr[9] = -0.01666666666666667*(alphaR[1]*(25.98076211353316*fR[19]-15.0*fR[16])+alphaR[0]*(25.98076211353316*fR[15]-15.0*fR[9])); 
  incr[10] = 0.05*(13.41640786499874*alphaR[1]*fR[18]-19.36491673103709*alphaR[0]*fR[17]+alphaR[1]*((-7.745966692414834*fR[14])-19.36491673103708*fR[13])+alphaR[0]*(15.0*fR[10]-8.660254037844386*fR[6])+alphaR[1]*(15.0*fR[5]-8.660254037844386*fR[3])); 
  incr[11] = -0.01666666666666667*(51.96152422706631*alphaR[1]*fR[12]-75.0*alphaR[0]*fR[11]+alphaR[1]*((-30.0*fR[8])-75.00000000000001*fR[7])+alphaR[0]*(58.09475019311126*fR[4]-33.54101966249684*fR[2])+alphaR[1]*(58.09475019311126*fR[1]-33.54101966249684*fR[0])); 
  incr[12] = 0.05*(15.0*alphaR[0]*fR[12]-17.32050807568877*alphaR[1]*fR[11]-8.660254037844387*alphaR[0]*fR[8]+alphaR[1]*(13.41640786499874*fR[4]-7.745966692414834*fR[2])); 
  incr[13] = 0.08333333333333333*(15.0*alphaR[1]*fR[17]+15.0*alphaR[0]*fR[13]+alphaR[1]*(6.708203932499369*fR[6]-11.61895003862225*fR[10])+alphaR[0]*(6.708203932499369*fR[3]-11.61895003862225*fR[5])); 
  incr[14] = -0.01666666666666667*(25.98076211353316*alphaR[0]*fR[18]-30.0*alphaR[1]*fR[17]-15.0*alphaR[0]*fR[14]+alphaR[1]*(23.2379000772445*fR[10]-13.41640786499874*fR[6])); 
  incr[15] = 0.05*(alphaR[1]*(15.0*fR[19]-8.660254037844386*fR[16])+alphaR[0]*(15.0*fR[15]-8.660254037844387*fR[9])); 
  incr[16] = -0.01666666666666667*(alphaR[0]*(25.98076211353316*fR[19]-15.0*fR[16])+alphaR[1]*(25.98076211353316*fR[15]-15.0*fR[9])); 
  incr[17] = -0.01666666666666667*(51.96152422706631*alphaR[1]*fR[18]-75.0*alphaR[0]*fR[17]+alphaR[1]*((-30.0*fR[14])-75.00000000000001*fR[13])+alphaR[0]*(58.09475019311126*fR[10]-33.54101966249685*fR[6])+alphaR[1]*(58.09475019311126*fR[5]-33.54101966249685*fR[3])); 
  incr[18] = 0.05*(15.0*alphaR[0]*fR[18]-17.32050807568877*alphaR[1]*fR[17]-8.660254037844387*alphaR[0]*fR[14]+alphaR[1]*(13.41640786499874*fR[10]-7.745966692414834*fR[6])); 
  incr[19] = 0.05*(alphaR[0]*(15.0*fR[19]-8.660254037844387*fR[16])+alphaR[1]*(15.0*fR[15]-8.660254037844386*fR[9])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.5*alphaR[0]-0.6708203932499369*alphaR[1]; 
  fUpOrd[0] = 0.5*(((-0.7348469228349525*(fR[19]+fL[19]))-0.7348469228349533*(fR[18]+fL[18])-1.42302494707577*fR[17]+1.42302494707577*fL[17]+0.4242640687119281*fR[16]-0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])+0.4242640687119285*fR[14]-0.4242640687119285*fL[14]+1.060660171779821*fR[13]-1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])+1.060660171779821*fR[11]-1.060660171779821*fL[11]+1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6363961030678926*fR[6]+0.6363961030678926*fL[6]-0.8215838362577489*(fR[5]+fL[5]+fR[4]+fL[4])+0.4743416490252568*(fR[3]+fR[2])-0.4743416490252568*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.7348469228349525*fR[19]-0.7348469228349525*fL[19]+0.7348469228349533*fR[18]-0.7348469228349533*fL[18]+1.42302494707577*(fR[17]+fL[17])-0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]-0.4242640687119285*(fR[14]+fL[14])-1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]-1.060660171779821*(fR[11]+fL[11])-1.10227038425243*fR[10]+1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])+0.6363961030678926*(fR[6]+fL[6])+0.8215838362577489*(fR[5]+fR[4])-0.8215838362577489*(fL[5]+fL[4])-0.4743416490252568*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[1] = 0.5*((0.6846531968814573*(fR[18]+fL[18])-0.3952847075210473*fR[14]+0.3952847075210473*fL[14]+0.7905694150420948*fR[13]-0.7905694150420948*fL[13]-0.6846531968814574*(fR[12]+fL[12])+0.3952847075210473*fR[8]-0.3952847075210473*fL[8]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6123724356957944*(fR[5]+fL[5])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6846531968814573*fR[18]+0.6846531968814573*fL[18]+0.3952847075210473*(fR[14]+fL[14])-0.7905694150420948*(fR[13]+fL[13])+0.6846531968814574*fR[12]-0.6846531968814574*fL[12]-0.3952847075210473*(fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])+0.6123724356957944*fR[5]-0.6123724356957944*fL[5]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[2] = 0.5*((0.7348469228349525*(fR[19]+fL[19])-0.7348469228349533*(fR[18]+fL[18])+1.42302494707577*fR[17]-1.42302494707577*fL[17]-0.4242640687119281*fR[16]+0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])+0.4242640687119285*fR[14]-0.4242640687119285*fL[14]+1.060660171779821*fR[13]-1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])-1.060660171779821*fR[11]+1.060660171779821*fL[11]-1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6363961030678926*fR[6]-0.6363961030678926*fL[6]-0.8215838362577489*(fR[5]+fL[5])+0.8215838362577489*(fR[4]+fL[4])+0.4743416490252568*fR[3]-0.4743416490252568*(fL[3]+fR[2])+0.4743416490252568*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.7348469228349525*fR[19]+0.7348469228349525*fL[19]+0.7348469228349533*fR[18]-0.7348469228349533*fL[18]-1.42302494707577*(fR[17]+fL[17])+0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]-0.4242640687119285*(fR[14]+fL[14])-1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]+1.060660171779821*(fR[11]+fL[11])+1.10227038425243*fR[10]-1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6363961030678926*(fR[6]+fL[6])+0.8215838362577489*fR[5]-0.8215838362577489*(fL[5]+fR[4])+0.8215838362577489*fL[4]-0.4743416490252568*(fR[3]+fL[3])+0.4743416490252568*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.5*alphaR[1]; 
  fUpOrd[3] = 0.5*((0.6846531968814573*(fR[19]+fL[19])-0.3952847075210473*fR[16]+0.3952847075210473*fL[16]-0.6846531968814574*(fR[15]+fL[15])+0.7905694150420948*fR[11]-0.7905694150420948*fL[11]+0.3952847075210473*fR[9]-0.3952847075210473*fL[9]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[2]-0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6846531968814573*fR[19]+0.6846531968814573*fL[19]+0.3952847075210473*(fR[16]+fL[16])+0.6846531968814574*fR[15]-0.6846531968814574*fL[15]-0.7905694150420948*(fR[11]+fL[11])-0.3952847075210473*(fR[9]+fL[9])+0.7905694150420947*(fR[7]+fL[7])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]-0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[1]+alphaR[0]); 
  fUpOrd[4] = 0.5*(((-0.6846531968814573*(fR[19]+fL[19]))+0.3952847075210473*fR[16]-0.3952847075210473*fL[16]-0.6846531968814574*(fR[15]+fL[15])-0.7905694150420948*fR[11]+0.7905694150420948*fL[11]+0.3952847075210473*fR[9]-0.3952847075210473*fL[9]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[2]+0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6846531968814573*fR[19]-0.6846531968814573*fL[19]-0.3952847075210473*(fR[16]+fL[16])+0.6846531968814574*fR[15]-0.6846531968814574*fL[15]+0.7905694150420948*(fR[11]+fL[11])-0.3952847075210473*(fR[9]+fL[9])+0.7905694150420947*(fR[7]+fL[7])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]+0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.6708203932499369*alphaR[1]; 
  fUpOrd[5] = 0.5*(((-0.7348469228349525*(fR[19]+fL[19]))+0.7348469228349533*(fR[18]+fL[18])+1.42302494707577*fR[17]-1.42302494707577*fL[17]+0.4242640687119281*fR[16]-0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])-0.4242640687119285*fR[14]+0.4242640687119285*fL[14]-1.060660171779821*fR[13]+1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])+1.060660171779821*fR[11]-1.060660171779821*fL[11]-1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6363961030678926*fR[6]-0.6363961030678926*fL[6]+0.8215838362577489*(fR[5]+fL[5])-0.8215838362577489*(fR[4]+fL[4])-0.4743416490252568*fR[3]+0.4743416490252568*(fL[3]+fR[2])-0.4743416490252568*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.7348469228349525*fR[19]-0.7348469228349525*fL[19]-0.7348469228349533*fR[18]+0.7348469228349533*fL[18]-1.42302494707577*(fR[17]+fL[17])-0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]+0.4242640687119285*(fR[14]+fL[14])+1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]-1.060660171779821*(fR[11]+fL[11])+1.10227038425243*fR[10]-1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6363961030678926*(fR[6]+fL[6])-0.8215838362577489*fR[5]+0.8215838362577489*(fL[5]+fR[4])-0.8215838362577489*fL[4]+0.4743416490252568*(fR[3]+fL[3])-0.4743416490252568*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[6] = 0.5*(((-0.6846531968814573*(fR[18]+fL[18]))+0.3952847075210473*fR[14]-0.3952847075210473*fL[14]-0.7905694150420948*fR[13]+0.7905694150420948*fL[13]-0.6846531968814574*(fR[12]+fL[12])+0.3952847075210473*fR[8]-0.3952847075210473*fL[8]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6123724356957944*(fR[5]+fL[5])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6846531968814573*fR[18]-0.6846531968814573*fL[18]-0.3952847075210473*(fR[14]+fL[14])+0.7905694150420948*(fR[13]+fL[13])+0.6846531968814574*fR[12]-0.6846531968814574*fL[12]-0.3952847075210473*(fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6123724356957944*fR[5]+0.6123724356957944*fL[5]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[7] = 0.5*((0.7348469228349525*(fR[19]+fL[19])+0.7348469228349533*(fR[18]+fL[18])-1.42302494707577*fR[17]+1.42302494707577*fL[17]-0.4242640687119281*fR[16]+0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])-0.4242640687119285*fR[14]+0.4242640687119285*fL[14]-1.060660171779821*fR[13]+1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])-1.060660171779821*fR[11]+1.060660171779821*fL[11]+1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6363961030678926*fR[6]+0.6363961030678926*fL[6]+0.8215838362577489*(fR[5]+fL[5]+fR[4]+fL[4])-0.4743416490252568*(fR[3]+fR[2])+0.4743416490252568*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.7348469228349525*fR[19]+0.7348469228349525*fL[19]-0.7348469228349533*fR[18]+0.7348469228349533*fL[18]+1.42302494707577*(fR[17]+fL[17])+0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]+0.4242640687119285*(fR[14]+fL[14])+1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]+1.060660171779821*(fR[11]+fL[11])-1.10227038425243*fR[10]+1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])+0.6363961030678926*(fR[6]+fL[6])-0.8215838362577489*(fR[5]+fR[4])+0.8215838362577489*(fL[5]+fL[4])+0.4743416490252568*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.03846153846153846*(5.0*fUpOrd[7]+8.0*fUpOrd[6]+5.0*fUpOrd[5]+8.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]+8.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[1] = 0.00828173324999922*(25.0*(fUpOrd[7]-1.0*fUpOrd[5])+53.66563145999496*(fUpOrd[4]-1.0*fUpOrd[3])+25.0*(fUpOrd[2]-1.0*fUpOrd[0])); 
  fUp[2] = 0.01851851851851852*(11.18033988749895*fUpOrd[7]+24.0*fUpOrd[6]+11.18033988749895*fUpOrd[5]-1.0*(11.18033988749895*fUpOrd[2]+24.0*fUpOrd[1]+11.18033988749895*fUpOrd[0])); 
  fUp[3] = 0.2777777777777778*(fUpOrd[7]-1.0*(fUpOrd[5]+fUpOrd[2])+fUpOrd[0]); 
  fUp[4] = 0.03440104580768907*(5.0*fUpOrd[7]-18.0*fUpOrd[6]+5.0*fUpOrd[5]+8.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]-18.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[5] = 0.03440104580768907*(5.0*fUpOrd[7]+8.0*fUpOrd[6]+5.0*fUpOrd[5]-18.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]+8.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[6] = 0.08281733249999224*(2.23606797749979*fUpOrd[7]-6.0*fUpOrd[6]+2.23606797749979*(fUpOrd[5]-1.0*fUpOrd[2])+6.0*fUpOrd[1]-2.23606797749979*fUpOrd[0]); 
  fUp[7] = 0.01656346649999844*(11.18033988749895*(fUpOrd[7]-1.0*fUpOrd[5])+10.0*(3.0*fUpOrd[3]-3.0*fUpOrd[4])+11.18033988749895*(fUpOrd[2]-1.0*fUpOrd[0])); 

  incr[0] = 0.3535533905932737*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = -0.6123724356957944*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[2] = 0.07071067811865474*(4.47213595499958*alphaR[1]*fUp[4]+5.0*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[3] = 0.3535533905932737*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2]); 
  incr[4] = -0.07071067811865474*(7.745966692414834*alphaR[1]*fUp[4]+8.660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[5] = -0.6123724356957944*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2]); 
  incr[6] = 0.02357022603955158*(13.41640786499874*alphaR[1]*fUp[6]+15.0*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[7] = 0.7905694150420947*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[8] = 0.07071067811865474*(5.0*alphaR[0]*fUp[4]+4.47213595499958*alphaR[1]*fUp[1]); 
  incr[9] = 0.02357022603955158*(15.0*alphaR[1]*fUp[7]+15.0*alphaR[0]*fUp[5]); 
  incr[10] = -0.07071067811865474*(7.745966692414834*alphaR[1]*fUp[6]+8.660254037844386*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[11] = 0.02357022603955158*(30.0*alphaR[1]*fUp[4]+33.54101966249684*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[12] = -0.07071067811865474*(8.660254037844387*alphaR[0]*fUp[4]+7.745966692414834*alphaR[1]*fUp[1]); 
  incr[13] = 0.7905694150420945*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2]); 
  incr[14] = 0.02357022603955158*(15.0*alphaR[0]*fUp[6]+13.41640786499874*alphaR[1]*fUp[3]); 
  incr[15] = -0.07071067811865474*(8.660254037844386*alphaR[1]*fUp[7]+8.660254037844387*alphaR[0]*fUp[5]); 
  incr[16] = 0.02357022603955158*(15.0*alphaR[0]*fUp[7]+15.0*alphaR[1]*fUp[5]); 
  incr[17] = 0.02357022603955158*(30.0*alphaR[1]*fUp[6]+33.54101966249685*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[18] = -0.07071067811865474*(8.660254037844387*alphaR[0]*fUp[6]+7.745966692414834*alphaR[1]*fUp[3]); 
  incr[19] = -0.07071067811865474*(8.660254037844387*alphaR[0]*fUp[7]+8.660254037844386*alphaR[1]*fUp[5]); 

#endif 
  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 
  outR[8] += incr[8]*rdx2R; 
  outR[9] += incr[9]*rdx2R; 
  outR[10] += incr[10]*rdx2R; 
  outR[11] += incr[11]*rdx2R; 
  outR[12] += incr[12]*rdx2R; 
  outR[13] += incr[13]*rdx2R; 
  outR[14] += incr[14]*rdx2R; 
  outR[15] += incr[15]*rdx2R; 
  outR[16] += incr[16]*rdx2R; 
  outR[17] += incr[17]*rdx2R; 
  outR[18] += incr[18]*rdx2R; 
  outR[19] += incr[19]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += -1.0*incr[7]*rdx2L; 
  outL[8] += -1.0*incr[8]*rdx2L; 
  outL[9] += -1.0*incr[9]*rdx2L; 
  outL[10] += incr[10]*rdx2L; 
  outL[11] += -1.0*incr[11]*rdx2L; 
  outL[12] += incr[12]*rdx2L; 
  outL[13] += -1.0*incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 
  outL[16] += -1.0*incr[16]*rdx2L; 
  outL[17] += -1.0*incr[17]*rdx2L; 
  outL[18] += incr[18]*rdx2L; 
  outL[19] += incr[19]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmGyrokineticSimpleHelicalSurf1x2vSer_vpar_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[20]; 
  hamilR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(bmag[0]*wmuR+phi[0]*q_))+2.0*m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamilR[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[7] = 2.0*(bmag[2]*wmuR+phi[2]*q_); 
  hamilR[8] = (0.421637021355784*m_)/rdvpar2SqR; 
  hamilR[13] = (1.154700538379251*bmag[2])/rdmu2R; 

  double BstarZdBmagR[20]; 
  BstarZdBmagR[0] = 2.0*cmag[0]; 
  BstarZdBmagR[1] = 2.0*cmag[1]; 
  BstarZdBmagR[7] = 2.0*cmag[2]; 

  double alphaR[8]; 
  alphaR[0] = -(0.25*((3.872983346207417*BstarZdBmagR[1]*hamilR[7]+1.732050807568877*BstarZdBmagR[0]*hamilR[1])*rdx2R+5.656854249492382*dApardt[0]*q_))/m_; 
  alphaR[1] = -(0.25*(((3.464101615137754*BstarZdBmagR[7]+3.872983346207417*BstarZdBmagR[0])*hamilR[7]+1.732050807568877*BstarZdBmagR[1]*hamilR[1])*rdx2R+5.656854249492382*dApardt[1]*q_))/m_; 
  alphaR[2] = -(0.25*(3.872983346207417*BstarZdBmagR[1]*hamilR[13]+1.732050807568877*BstarZdBmagR[0]*hamilR[5])*rdx2R)/m_; 
  alphaR[3] = -(0.05*((17.32050807568877*BstarZdBmagR[7]+19.36491673103708*BstarZdBmagR[0])*hamilR[13]+8.660254037844386*BstarZdBmagR[1]*hamilR[5])*rdx2R)/m_; 
  alphaR[4] = -(0.25*((3.464101615137754*BstarZdBmagR[1]*hamilR[7]+1.732050807568877*hamilR[1]*BstarZdBmagR[7])*rdx2R+5.656854249492382*dApardt[2]*q_))/m_; 
  alphaR[6] = -(0.05*(17.32050807568877*BstarZdBmagR[1]*hamilR[13]+8.660254037844387*hamilR[5]*BstarZdBmagR[7])*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*((3.872983346207417*BstarZdBmagR[1]*hamilR[7]+1.732050807568877*BstarZdBmagR[0]*hamilR[1])*rdx2R+5.656854249492382*dApardt[0]*q_))/m_; 

  double incr[20]; 
  double amax = amax_in; 

  double fAvg[8]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[8]+fL[8])+1.732050807568877*(fL[2]-1.0*fR[2])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[12]+fL[12])+3.0*(fL[4]-1.0*fR[4]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[14]+fL[14])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.7071067811865475*(2.23606797749979*(fR[18]+fL[18])+1.732050807568877*(fL[10]-1.0*fR[10])+fR[5]+fL[5]); 
  fAvg[4] = -0.1414213562373095*(8.660254037844387*fR[11]-1.0*(8.660254037844387*fL[11]+5.0*(fR[7]+fL[7]))); 
  fAvg[5] = -0.1414213562373095*(8.660254037844387*fR[16]-1.0*(8.660254037844387*fL[16]+5.0*(fR[9]+fL[9]))); 
  fAvg[6] = -0.1414213562373095*(8.660254037844387*fR[17]-1.0*(8.660254037844387*fL[17]+5.0*(fR[13]+fL[13]))); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[19]-1.0*(8.660254037844387*fL[19]+5.0*(fR[15]+fL[15]))); 

  double Ghat[8]; 
  Ghat[0] = -0.1767766952966368*((4.47213595499958*fR[8]-4.47213595499958*fL[8]-3.464101615137754*(fR[2]+fL[2])+2.0*fR[0]-2.0*fL[0])*amax-1.414213562373095*(alphaR[6]*fAvg[6]+alphaR[4]*fAvg[4]+alphaR[3]*fAvg[3]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.01178511301977579*((67.08203932499369*fR[12]-67.08203932499369*fL[12]-51.96152422706631*(fR[4]+fL[4])+30.0*fR[1]-30.0*fL[1])*amax-18.97366596101028*(alphaR[3]*fAvg[6]+fAvg[3]*alphaR[6]+alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])-21.21320343559643*(alphaR[2]*fAvg[3]+fAvg[2]*alphaR[3]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.01178511301977579*((67.08203932499369*fR[14]-67.08203932499369*fL[14]-51.96152422706631*(fR[6]+fL[6])+30.0*fR[3]-30.0*fL[3])*amax-18.97366596101028*alphaR[3]*fAvg[7]-21.21320343559643*(alphaR[4]*fAvg[6]+fAvg[4]*alphaR[6])-18.97366596101028*alphaR[2]*fAvg[5]-21.21320343559643*(alphaR[1]*fAvg[3]+fAvg[1]*alphaR[3]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.01178511301977579*((67.0820393249937*fR[18]-67.0820393249937*fL[18]-51.96152422706631*(fR[10]+fL[10])+30.0*fR[5]-30.0*fL[5])*amax+((-16.97056274847715*alphaR[6])-18.97366596101028*alphaR[2])*fAvg[7]-18.97366596101028*(alphaR[1]*fAvg[6]+fAvg[1]*alphaR[6]+alphaR[3]*(fAvg[5]+fAvg[4])+fAvg[3]*alphaR[4])-21.21320343559643*(alphaR[0]*fAvg[3]+fAvg[0]*alphaR[3]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[4] = 0.001683587574253684*((363.7306695894642*(fR[11]+fL[11])-210.0*fR[7]+210.0*fL[7])*amax+(94.86832980505142*alphaR[6]+148.492424049175*alphaR[2])*fAvg[6]+148.492424049175*fAvg[2]*alphaR[6]+(94.86832980505142*alphaR[4]+148.492424049175*alphaR[0])*fAvg[4]+148.492424049175*fAvg[0]*alphaR[4]+132.815661727072*(alphaR[3]*fAvg[3]+alphaR[1]*fAvg[1])); 
  Ghat[5] = 0.01178511301977579*((51.96152422706632*(fR[16]+fL[16])-30.0*fR[9]+30.0*fL[9])*amax+21.21320343559643*alphaR[1]*fAvg[7]+18.97366596101028*alphaR[6]*fAvg[6]+21.21320343559643*alphaR[0]*fAvg[5]+18.97366596101028*(alphaR[3]*fAvg[3]+alphaR[2]*fAvg[2])); 
  Ghat[6] = 0.001683587574253684*((363.7306695894642*(fR[17]+fL[17])-210.0*fR[13]+210.0*fL[13])*amax+118.79393923934*alphaR[3]*fAvg[7]+(94.86832980505142*alphaR[4]+148.492424049175*alphaR[0])*fAvg[6]+(132.815661727072*fAvg[5]+94.86832980505142*fAvg[4]+148.492424049175*fAvg[0])*alphaR[6]+148.492424049175*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])+132.815661727072*(alphaR[1]*fAvg[3]+fAvg[1]*alphaR[3])); 
  Ghat[7] = 0.01178511301977579*((51.96152422706632*(fR[19]+fL[19])-30.0*fR[15]+30.0*fL[15])*amax+(18.97366596101028*alphaR[4]+21.21320343559643*alphaR[0])*fAvg[7]+16.97056274847715*(alphaR[3]*fAvg[6]+fAvg[3]*alphaR[6])+21.21320343559643*alphaR[1]*fAvg[5]+18.97366596101028*(alphaR[2]*fAvg[3]+fAvg[2]*alphaR[3])); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = -1.224744871391589*Ghat[0]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = -1.224744871391589*Ghat[1]; 
  incr[5] = 0.7071067811865475*Ghat[3]; 
  incr[6] = -1.224744871391589*Ghat[2]; 
  incr[7] = 0.7071067811865475*Ghat[4]; 
  incr[8] = 1.58113883008419*Ghat[0]; 
  incr[9] = 0.7071067811865475*Ghat[5]; 
  incr[10] = -1.224744871391589*Ghat[3]; 
  incr[11] = -1.224744871391589*Ghat[4]; 
  incr[12] = 1.58113883008419*Ghat[1]; 
  incr[13] = 0.7071067811865475*Ghat[6]; 
  incr[14] = 1.58113883008419*Ghat[2]; 
  incr[15] = 0.7071067811865475*Ghat[7]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = -1.224744871391589*Ghat[6]; 
  incr[18] = 1.58113883008419*Ghat[3]; 
  incr[19] = -1.224744871391589*Ghat[7]; 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 
  outR[8] += incr[8]*rdvpar2R; 
  outR[9] += incr[9]*rdvpar2R; 
  outR[10] += incr[10]*rdvpar2R; 
  outR[11] += incr[11]*rdvpar2R; 
  outR[12] += incr[12]*rdvpar2R; 
  outR[13] += incr[13]*rdvpar2R; 
  outR[14] += incr[14]*rdvpar2R; 
  outR[15] += incr[15]*rdvpar2R; 
  outR[16] += incr[16]*rdvpar2R; 
  outR[17] += incr[17]*rdvpar2R; 
  outR[18] += incr[18]*rdvpar2R; 
  outR[19] += incr[19]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += -1.0*incr[7]*rdvpar2L; 
  outL[8] += -1.0*incr[8]*rdvpar2L; 
  outL[9] += -1.0*incr[9]*rdvpar2L; 
  outL[10] += incr[10]*rdvpar2L; 
  outL[11] += incr[11]*rdvpar2L; 
  outL[12] += -1.0*incr[12]*rdvpar2L; 
  outL[13] += -1.0*incr[13]*rdvpar2L; 
  outL[14] += -1.0*incr[14]*rdvpar2L; 
  outL[15] += -1.0*incr[15]*rdvpar2L; 
  outL[16] += incr[16]*rdvpar2L; 
  outL[17] += incr[17]*rdvpar2L; 
  outL[18] += -1.0*incr[18]*rdvpar2L; 
  outL[19] += incr[19]*rdvpar2L; 
return std::abs(alphaSurfAvgR); 
} 
