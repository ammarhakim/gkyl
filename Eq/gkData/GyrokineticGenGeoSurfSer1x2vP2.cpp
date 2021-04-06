#include <GyrokineticModDecl.h>
double GyrokineticGenGeoSurf1x2vSer_x_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

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
double GyrokineticGenGeoSurf1x2vSer_vpar_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double alphaR[8]; 
  alphaR[0] = -(0.4330127018922193*BstarZdBmagR[0]*hamilR[1]*rdx2R)/m_; 
  alphaR[1] = -(0.9682458365518543*BstarZdBmagR[0]*hamilR[7]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.1082531754730548*BstarZdBmagR[0]*hamilR[1]*rdx2R)/m_; 

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
  Ghat[0] = -0.1767766952966368*((4.47213595499958*fR[8]-4.47213595499958*fL[8]-3.464101615137754*(fR[2]+fL[2])+2.0*fR[0]-2.0*fL[0])*amax-1.414213562373095*(alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.01178511301977579*((67.08203932499369*fR[12]-67.08203932499369*fL[12]-51.96152422706631*(fR[4]+fL[4])+30.0*fR[1]-30.0*fL[1])*amax-18.97366596101028*alphaR[1]*fAvg[4]-21.21320343559643*(alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.05892556509887893*((13.41640786499874*fR[14]-13.41640786499874*fL[14]-10.39230484541326*(fR[6]+fL[6])+6.0*fR[3]-6.0*fL[3])*amax-4.242640687119286*(alphaR[1]*fAvg[3]+alphaR[0]*fAvg[2])); 
  Ghat[3] = -0.01178511301977579*((67.0820393249937*fR[18]-67.0820393249937*fL[18]-51.96152422706631*(fR[10]+fL[10])+30.0*fR[5]-30.0*fL[5])*amax-18.97366596101028*alphaR[1]*fAvg[6]-21.21320343559643*(alphaR[0]*fAvg[3]+alphaR[1]*fAvg[2])); 
  Ghat[4] = 0.03535533905932736*((17.32050807568877*(fR[11]+fL[11])-10.0*fR[7]+10.0*fL[7])*amax+7.071067811865476*alphaR[0]*fAvg[4]+6.324555320336761*alphaR[1]*fAvg[1]); 
  Ghat[5] = 0.01178511301977579*((51.96152422706632*(fR[16]+fL[16])-30.0*fR[9]+30.0*fL[9])*amax+21.21320343559643*alphaR[1]*fAvg[7]+21.21320343559643*alphaR[0]*fAvg[5]); 
  Ghat[6] = 0.01178511301977579*((51.96152422706632*(fR[17]+fL[17])-30.0*fR[13]+30.0*fL[13])*amax+21.21320343559643*alphaR[0]*fAvg[6]+18.97366596101028*alphaR[1]*fAvg[3]); 
  Ghat[7] = 0.01178511301977579*((51.96152422706632*(fR[19]+fL[19])-30.0*fR[15]+30.0*fL[15])*amax+21.21320343559643*alphaR[0]*fAvg[7]+21.21320343559643*alphaR[1]*fAvg[5]); 

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
double GyrokineticGenGeoSurf1x2vSer_x_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R*wvparR+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (0.2828427124746191*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2R*wvparR+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*(b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[7] = (0.04040610178208843*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R*wvparR+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmagR[11] = (1.414213562373095*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = (0.25*(hamilR[8]*(8.660254037844387*BstarZdBmagR[11]-6.708203932499369*BstarZdBmagR[4]+3.872983346207417*BstarZdBmagR[2])+hamilR[2]*(3.872983346207417*BstarZdBmagR[7]-3.0*BstarZdBmagR[1]+1.732050807568877*BstarZdBmagR[0]))*rdvpar2R)/m_; 
  alphaR[1] = (0.25*(3.872983346207417*hamilR[2]*BstarZdBmagR[11]+(8.660254037844386*BstarZdBmagR[7]-6.708203932499369*BstarZdBmagR[1]+3.872983346207417*BstarZdBmagR[0])*hamilR[8]+hamilR[2]*(1.732050807568877*BstarZdBmagR[2]-3.0*BstarZdBmagR[4]))*rdvpar2R)/m_; 
  alphaR[4] = (0.5*hamilR[8]*(3.872983346207417*BstarZdBmagR[11]-3.0*BstarZdBmagR[4]+1.732050807568877*BstarZdBmagR[2])*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(8.660254037844387*hamilR[8]*BstarZdBmagR[11]+(3.872983346207417*BstarZdBmagR[2]-6.708203932499369*BstarZdBmagR[4])*hamilR[8]+3.872983346207417*hamilR[2]*BstarZdBmagR[7]+(1.732050807568877*BstarZdBmagR[0]-3.0*BstarZdBmagR[1])*hamilR[2])*rdvpar2R)/m_; 

  double incr[20]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.01666666666666667*(25.98076211353316*alphaR[4]*fL[12]+33.54101966249684*alphaR[1]*fL[11]+15.0*alphaR[4]*fL[8]+33.54101966249685*alphaR[0]*fL[7]+alphaR[1]*(25.98076211353316*fL[4]+15.0*fL[2])+alphaR[0]*(25.98076211353316*fL[1]+15.0*fL[0])); 
  incr[1] = -0.05*(15.0*alphaR[4]*fL[12]+19.36491673103708*alphaR[1]*fL[11]+8.660254037844386*alphaR[4]*fL[8]+19.36491673103709*alphaR[0]*fL[7]+alphaR[1]*(15.0*fL[4]+8.660254037844386*fL[2])+alphaR[0]*(15.0*fL[1]+8.660254037844386*fL[0])); 
  incr[2] = 0.01666666666666667*(23.2379000772445*alphaR[1]*fL[12]+(30.0*alphaR[4]+33.54101966249684*alphaR[0])*fL[11]+alphaR[1]*(13.41640786499874*fL[8]+33.54101966249685*fL[7])+(23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fL[4]+fL[2]*(13.41640786499874*alphaR[4]+15.0*alphaR[0])+alphaR[1]*(25.98076211353316*fL[1]+15.0*fL[0])); 
  incr[3] = 0.01666666666666667*(25.98076211353316*alphaR[4]*fL[18]+33.54101966249685*alphaR[1]*fL[17]+15.0*alphaR[4]*fL[14]+33.54101966249684*alphaR[0]*fL[13]+alphaR[1]*(25.98076211353316*fL[10]+15.0*fL[6])+alphaR[0]*(25.98076211353316*fL[5]+15.0*fL[3])); 
  incr[4] = -0.05*(13.41640786499874*alphaR[1]*fL[12]+(17.32050807568877*alphaR[4]+19.36491673103708*alphaR[0])*fL[11]+alphaR[1]*(7.745966692414834*fL[8]+19.36491673103709*fL[7])+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fL[4]+fL[2]*(7.745966692414834*alphaR[4]+8.660254037844386*alphaR[0])+alphaR[1]*(15.0*fL[1]+8.660254037844386*fL[0])); 
  incr[5] = -0.05*(15.0*alphaR[4]*fL[18]+19.36491673103709*alphaR[1]*fL[17]+8.660254037844387*alphaR[4]*fL[14]+19.36491673103708*alphaR[0]*fL[13]+alphaR[1]*(15.0*fL[10]+8.660254037844386*fL[6])+alphaR[0]*(15.0*fL[5]+8.660254037844386*fL[3])); 
  incr[6] = 0.01666666666666667*(23.2379000772445*alphaR[1]*fL[18]+(30.0*alphaR[4]+33.54101966249685*alphaR[0])*fL[17]+alphaR[1]*(13.41640786499874*fL[14]+33.54101966249684*fL[13])+(23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fL[10]+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fL[6]+alphaR[1]*(25.98076211353316*fL[5]+15.0*fL[3])); 
  incr[7] = 0.08333333333333333*(11.61895003862225*alphaR[4]*fL[12]+15.0*alphaR[1]*fL[11]+6.708203932499369*alphaR[4]*fL[8]+15.0*alphaR[0]*fL[7]+alphaR[1]*(11.61895003862225*fL[4]+6.708203932499369*fL[2])+alphaR[0]*(11.61895003862225*fL[1]+6.708203932499369*fL[0])); 
  incr[8] = 0.002380952380952381*((116.1895003862225*alphaR[4]+181.8653347947321*alphaR[0])*fL[12]+210.0*alphaR[1]*fL[11]+(67.0820393249937*alphaR[4]+105.0*alphaR[0])*fL[8]+234.787137637478*alphaR[4]*fL[7]+162.6653005407115*alphaR[1]*fL[4]+(181.8653347947321*fL[1]+105.0*fL[0])*alphaR[4]+93.91485505499116*alphaR[1]*fL[2]); 
  incr[9] = 0.01666666666666667*(alphaR[1]*(25.98076211353316*fL[19]+15.0*fL[16])+alphaR[0]*(25.98076211353316*fL[15]+15.0*fL[9])); 
  incr[10] = -0.05*(13.41640786499874*alphaR[1]*fL[18]+(17.32050807568877*alphaR[4]+19.36491673103709*alphaR[0])*fL[17]+alphaR[1]*(7.745966692414834*fL[14]+19.36491673103708*fL[13])+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fL[10]+(7.745966692414834*alphaR[4]+8.660254037844386*alphaR[0])*fL[6]+alphaR[1]*(15.0*fL[5]+8.660254037844386*fL[3])); 
  incr[11] = 0.01666666666666667*(51.96152422706631*alphaR[1]*fL[12]+(67.0820393249937*alphaR[4]+75.0*alphaR[0])*fL[11]+alphaR[1]*(30.0*fL[8]+75.00000000000001*fL[7])+(51.96152422706632*alphaR[4]+58.09475019311126*alphaR[0])*fL[4]+fL[2]*(30.0*alphaR[4]+33.54101966249684*alphaR[0])+alphaR[1]*(58.09475019311126*fL[1]+33.54101966249684*fL[0])); 
  incr[12] = -0.007142857142857143*((67.0820393249937*alphaR[4]+105.0*alphaR[0])*fL[12]+121.2435565298214*alphaR[1]*fL[11]+(38.72983346207417*alphaR[4]+60.62177826491071*alphaR[0])*fL[8]+135.5544171172596*alphaR[4]*fL[7]+93.91485505499116*alphaR[1]*fL[4]+(105.0*fL[1]+60.62177826491071*fL[0])*alphaR[4]+54.22176684690384*alphaR[1]*fL[2]); 
  incr[13] = 0.08333333333333333*(11.61895003862225*alphaR[4]*fL[18]+15.0*alphaR[1]*fL[17]+6.708203932499369*alphaR[4]*fL[14]+15.0*alphaR[0]*fL[13]+alphaR[1]*(11.61895003862225*fL[10]+6.708203932499369*fL[6])+alphaR[0]*(11.61895003862225*fL[5]+6.708203932499369*fL[3])); 
  incr[14] = 0.002380952380952381*((116.1895003862225*alphaR[4]+181.8653347947321*alphaR[0])*fL[18]+210.0*alphaR[1]*fL[17]+(67.0820393249937*alphaR[4]+105.0*alphaR[0])*fL[14]+234.787137637478*alphaR[4]*fL[13]+alphaR[1]*(162.6653005407115*fL[10]+93.91485505499116*fL[6])+alphaR[4]*(181.8653347947321*fL[5]+105.0*fL[3])); 
  incr[15] = -0.05*(alphaR[1]*(15.0*fL[19]+8.660254037844386*fL[16])+alphaR[0]*(15.0*fL[15]+8.660254037844387*fL[9])); 
  incr[16] = 0.01666666666666667*((23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fL[19]+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fL[16]+alphaR[1]*(25.98076211353316*fL[15]+15.0*fL[9])); 
  incr[17] = 0.01666666666666667*(51.96152422706631*alphaR[1]*fL[18]+(67.0820393249937*alphaR[4]+75.0*alphaR[0])*fL[17]+alphaR[1]*(30.0*fL[14]+75.00000000000001*fL[13])+(51.96152422706631*alphaR[4]+58.09475019311126*alphaR[0])*fL[10]+(30.0*alphaR[4]+33.54101966249685*alphaR[0])*fL[6]+alphaR[1]*(58.09475019311126*fL[5]+33.54101966249685*fL[3])); 
  incr[18] = -0.007142857142857143*((67.0820393249937*alphaR[4]+105.0*alphaR[0])*fL[18]+121.2435565298214*alphaR[1]*fL[17]+(38.72983346207417*alphaR[4]+60.62177826491071*alphaR[0])*fL[14]+135.5544171172596*alphaR[4]*fL[13]+alphaR[1]*(93.91485505499116*fL[10]+54.22176684690384*fL[6])+alphaR[4]*(105.0*fL[5]+60.6217782649107*fL[3])); 
  incr[19] = -0.05*((13.41640786499874*alphaR[4]+15.0*alphaR[0])*fL[19]+(7.745966692414834*alphaR[4]+8.660254037844387*alphaR[0])*fL[16]+alphaR[1]*(15.0*fL[15]+8.660254037844386*fL[9])); 
  } else { 
  incr[0] = -0.01666666666666667*(25.98076211353316*alphaR[4]*fR[12]-33.54101966249684*alphaR[1]*fR[11]-15.0*alphaR[4]*fR[8]-33.54101966249685*alphaR[0]*fR[7]+alphaR[1]*(25.98076211353316*fR[4]-15.0*fR[2])+alphaR[0]*(25.98076211353316*fR[1]-15.0*fR[0])); 
  incr[1] = 0.05*(15.0*alphaR[4]*fR[12]-19.36491673103708*alphaR[1]*fR[11]-8.660254037844386*alphaR[4]*fR[8]-19.36491673103709*alphaR[0]*fR[7]+alphaR[1]*(15.0*fR[4]-8.660254037844386*fR[2])+alphaR[0]*(15.0*fR[1]-8.660254037844386*fR[0])); 
  incr[2] = -0.01666666666666667*(23.2379000772445*alphaR[1]*fR[12]+((-30.0*alphaR[4])-33.54101966249684*alphaR[0])*fR[11]+alphaR[1]*((-13.41640786499874*fR[8])-33.54101966249685*fR[7])+(23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fR[4]+fR[2]*((-13.41640786499874*alphaR[4])-15.0*alphaR[0])+alphaR[1]*(25.98076211353316*fR[1]-15.0*fR[0])); 
  incr[3] = -0.01666666666666667*(25.98076211353316*alphaR[4]*fR[18]-33.54101966249685*alphaR[1]*fR[17]-15.0*alphaR[4]*fR[14]-33.54101966249684*alphaR[0]*fR[13]+alphaR[1]*(25.98076211353316*fR[10]-15.0*fR[6])+alphaR[0]*(25.98076211353316*fR[5]-15.0*fR[3])); 
  incr[4] = 0.05*(13.41640786499874*alphaR[1]*fR[12]+((-17.32050807568877*alphaR[4])-19.36491673103708*alphaR[0])*fR[11]+alphaR[1]*((-7.745966692414834*fR[8])-19.36491673103709*fR[7])+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fR[4]+fR[2]*((-7.745966692414834*alphaR[4])-8.660254037844386*alphaR[0])+alphaR[1]*(15.0*fR[1]-8.660254037844386*fR[0])); 
  incr[5] = 0.05*(15.0*alphaR[4]*fR[18]-19.36491673103709*alphaR[1]*fR[17]-8.660254037844387*alphaR[4]*fR[14]-19.36491673103708*alphaR[0]*fR[13]+alphaR[1]*(15.0*fR[10]-8.660254037844386*fR[6])+alphaR[0]*(15.0*fR[5]-8.660254037844386*fR[3])); 
  incr[6] = -0.01666666666666667*(23.2379000772445*alphaR[1]*fR[18]+((-30.0*alphaR[4])-33.54101966249685*alphaR[0])*fR[17]+alphaR[1]*((-13.41640786499874*fR[14])-33.54101966249684*fR[13])+(23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fR[10]+((-13.41640786499874*alphaR[4])-15.0*alphaR[0])*fR[6]+alphaR[1]*(25.98076211353316*fR[5]-15.0*fR[3])); 
  incr[7] = -0.08333333333333333*(11.61895003862225*alphaR[4]*fR[12]-15.0*alphaR[1]*fR[11]-6.708203932499369*alphaR[4]*fR[8]-15.0*alphaR[0]*fR[7]+alphaR[1]*(11.61895003862225*fR[4]-6.708203932499369*fR[2])+alphaR[0]*(11.61895003862225*fR[1]-6.708203932499369*fR[0])); 
  incr[8] = -0.002380952380952381*((116.1895003862225*alphaR[4]+181.8653347947321*alphaR[0])*fR[12]-210.0*alphaR[1]*fR[11]+((-67.0820393249937*alphaR[4])-105.0*alphaR[0])*fR[8]-234.787137637478*alphaR[4]*fR[7]+162.6653005407115*alphaR[1]*fR[4]+(181.8653347947321*fR[1]-105.0*fR[0])*alphaR[4]-93.91485505499116*alphaR[1]*fR[2]); 
  incr[9] = -0.01666666666666667*(alphaR[1]*(25.98076211353316*fR[19]-15.0*fR[16])+alphaR[0]*(25.98076211353316*fR[15]-15.0*fR[9])); 
  incr[10] = 0.05*(13.41640786499874*alphaR[1]*fR[18]+((-17.32050807568877*alphaR[4])-19.36491673103709*alphaR[0])*fR[17]+alphaR[1]*((-7.745966692414834*fR[14])-19.36491673103708*fR[13])+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fR[10]+((-7.745966692414834*alphaR[4])-8.660254037844386*alphaR[0])*fR[6]+alphaR[1]*(15.0*fR[5]-8.660254037844386*fR[3])); 
  incr[11] = -0.01666666666666667*(51.96152422706631*alphaR[1]*fR[12]+((-67.0820393249937*alphaR[4])-75.0*alphaR[0])*fR[11]+alphaR[1]*((-30.0*fR[8])-75.00000000000001*fR[7])+(51.96152422706632*alphaR[4]+58.09475019311126*alphaR[0])*fR[4]+fR[2]*((-30.0*alphaR[4])-33.54101966249684*alphaR[0])+alphaR[1]*(58.09475019311126*fR[1]-33.54101966249684*fR[0])); 
  incr[12] = 0.007142857142857143*((67.0820393249937*alphaR[4]+105.0*alphaR[0])*fR[12]-121.2435565298214*alphaR[1]*fR[11]+((-38.72983346207417*alphaR[4])-60.62177826491071*alphaR[0])*fR[8]-135.5544171172596*alphaR[4]*fR[7]+93.91485505499116*alphaR[1]*fR[4]+(105.0*fR[1]-60.62177826491071*fR[0])*alphaR[4]-54.22176684690384*alphaR[1]*fR[2]); 
  incr[13] = -0.08333333333333333*(11.61895003862225*alphaR[4]*fR[18]-15.0*alphaR[1]*fR[17]-6.708203932499369*alphaR[4]*fR[14]-15.0*alphaR[0]*fR[13]+alphaR[1]*(11.61895003862225*fR[10]-6.708203932499369*fR[6])+alphaR[0]*(11.61895003862225*fR[5]-6.708203932499369*fR[3])); 
  incr[14] = -0.002380952380952381*((116.1895003862225*alphaR[4]+181.8653347947321*alphaR[0])*fR[18]-210.0*alphaR[1]*fR[17]+((-67.0820393249937*alphaR[4])-105.0*alphaR[0])*fR[14]-234.787137637478*alphaR[4]*fR[13]+alphaR[1]*(162.6653005407115*fR[10]-93.91485505499116*fR[6])+alphaR[4]*(181.8653347947321*fR[5]-105.0*fR[3])); 
  incr[15] = 0.05*(alphaR[1]*(15.0*fR[19]-8.660254037844386*fR[16])+alphaR[0]*(15.0*fR[15]-8.660254037844387*fR[9])); 
  incr[16] = -0.01666666666666667*((23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fR[19]+((-13.41640786499874*alphaR[4])-15.0*alphaR[0])*fR[16]+alphaR[1]*(25.98076211353316*fR[15]-15.0*fR[9])); 
  incr[17] = -0.01666666666666667*(51.96152422706631*alphaR[1]*fR[18]+((-67.0820393249937*alphaR[4])-75.0*alphaR[0])*fR[17]+alphaR[1]*((-30.0*fR[14])-75.00000000000001*fR[13])+(51.96152422706631*alphaR[4]+58.09475019311126*alphaR[0])*fR[10]+((-30.0*alphaR[4])-33.54101966249685*alphaR[0])*fR[6]+alphaR[1]*(58.09475019311126*fR[5]-33.54101966249685*fR[3])); 
  incr[18] = 0.007142857142857143*((67.0820393249937*alphaR[4]+105.0*alphaR[0])*fR[18]-121.2435565298214*alphaR[1]*fR[17]+((-38.72983346207417*alphaR[4])-60.62177826491071*alphaR[0])*fR[14]-135.5544171172596*alphaR[4]*fR[13]+alphaR[1]*(93.91485505499116*fR[10]-54.22176684690384*fR[6])+alphaR[4]*(105.0*fR[5]-60.6217782649107*fR[3])); 
  incr[19] = 0.05*((13.41640786499874*alphaR[4]+15.0*alphaR[0])*fR[19]+((-7.745966692414834*alphaR[4])-8.660254037844387*alphaR[0])*fR[16]+alphaR[1]*(15.0*fR[15]-8.660254037844386*fR[9])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.4472135954999579*alphaR[4]-0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[0] = 0.5*(((-0.7348469228349525*(fR[19]+fL[19]))-0.7348469228349533*(fR[18]+fL[18])-1.42302494707577*fR[17]+1.42302494707577*fL[17]+0.4242640687119281*fR[16]-0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])+0.4242640687119285*fR[14]-0.4242640687119285*fL[14]+1.060660171779821*fR[13]-1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])+1.060660171779821*fR[11]-1.060660171779821*fL[11]+1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6363961030678926*fR[6]+0.6363961030678926*fL[6]-0.8215838362577489*(fR[5]+fL[5]+fR[4]+fL[4])+0.4743416490252568*(fR[3]+fR[2])-0.4743416490252568*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.7348469228349525*fR[19]-0.7348469228349525*fL[19]+0.7348469228349533*fR[18]-0.7348469228349533*fL[18]+1.42302494707577*(fR[17]+fL[17])-0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]-0.4242640687119285*(fR[14]+fL[14])-1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]-1.060660171779821*(fR[11]+fL[11])-1.10227038425243*fR[10]+1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])+0.6363961030678926*(fR[6]+fL[6])+0.8215838362577489*(fR[5]+fR[4])-0.8215838362577489*(fL[5]+fL[4])-0.4743416490252568*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.5590169943749475*alphaR[4]; 
  fUpOrd[1] = 0.5*((0.6846531968814573*(fR[18]+fL[18])-0.3952847075210473*fR[14]+0.3952847075210473*fL[14]+0.7905694150420948*fR[13]-0.7905694150420948*fL[13]-0.6846531968814574*(fR[12]+fL[12])+0.3952847075210473*fR[8]-0.3952847075210473*fL[8]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6123724356957944*(fR[5]+fL[5])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6846531968814573*fR[18]+0.6846531968814573*fL[18]+0.3952847075210473*(fR[14]+fL[14])-0.7905694150420948*(fR[13]+fL[13])+0.6846531968814574*fR[12]-0.6846531968814574*fL[12]-0.3952847075210473*(fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])+0.6123724356957944*fR[5]-0.6123724356957944*fL[5]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.4472135954999579*alphaR[4]+0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[2] = 0.5*((0.7348469228349525*(fR[19]+fL[19])-0.7348469228349533*(fR[18]+fL[18])+1.42302494707577*fR[17]-1.42302494707577*fL[17]-0.4242640687119281*fR[16]+0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])+0.4242640687119285*fR[14]-0.4242640687119285*fL[14]+1.060660171779821*fR[13]-1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])-1.060660171779821*fR[11]+1.060660171779821*fL[11]-1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6363961030678926*fR[6]-0.6363961030678926*fL[6]-0.8215838362577489*(fR[5]+fL[5])+0.8215838362577489*(fR[4]+fL[4])+0.4743416490252568*fR[3]-0.4743416490252568*(fL[3]+fR[2])+0.4743416490252568*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.7348469228349525*fR[19]+0.7348469228349525*fL[19]+0.7348469228349533*fR[18]-0.7348469228349533*fL[18]-1.42302494707577*(fR[17]+fL[17])+0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]-0.4242640687119285*(fR[14]+fL[14])-1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]+1.060660171779821*(fR[11]+fL[11])+1.10227038425243*fR[10]-1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6363961030678926*(fR[6]+fL[6])+0.8215838362577489*fR[5]-0.8215838362577489*(fL[5]+fR[4])+0.8215838362577489*fL[4]-0.4743416490252568*(fR[3]+fL[3])+0.4743416490252568*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.5*alphaR[1]; 
  fUpOrd[3] = 0.5*((0.6846531968814573*(fR[19]+fL[19])-0.3952847075210473*fR[16]+0.3952847075210473*fL[16]-0.6846531968814574*(fR[15]+fL[15])+0.7905694150420948*fR[11]-0.7905694150420948*fL[11]+0.3952847075210473*fR[9]-0.3952847075210473*fL[9]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[2]-0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6846531968814573*fR[19]+0.6846531968814573*fL[19]+0.3952847075210473*(fR[16]+fL[16])+0.6846531968814574*fR[15]-0.6846531968814574*fL[15]-0.7905694150420948*(fR[11]+fL[11])-0.3952847075210473*(fR[9]+fL[9])+0.7905694150420947*(fR[7]+fL[7])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]-0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[1]+alphaR[0]); 
  fUpOrd[4] = 0.5*(((-0.6846531968814573*(fR[19]+fL[19]))+0.3952847075210473*fR[16]-0.3952847075210473*fL[16]-0.6846531968814574*(fR[15]+fL[15])-0.7905694150420948*fR[11]+0.7905694150420948*fL[11]+0.3952847075210473*fR[9]-0.3952847075210473*fL[9]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[2]+0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6846531968814573*fR[19]-0.6846531968814573*fL[19]-0.3952847075210473*(fR[16]+fL[16])+0.6846531968814574*fR[15]-0.6846531968814574*fL[15]+0.7905694150420948*(fR[11]+fL[11])-0.3952847075210473*(fR[9]+fL[9])+0.7905694150420947*(fR[7]+fL[7])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]+0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.4472135954999579*alphaR[4]-0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[5] = 0.5*(((-0.7348469228349525*(fR[19]+fL[19]))+0.7348469228349533*(fR[18]+fL[18])+1.42302494707577*fR[17]-1.42302494707577*fL[17]+0.4242640687119281*fR[16]-0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])-0.4242640687119285*fR[14]+0.4242640687119285*fL[14]-1.060660171779821*fR[13]+1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])+1.060660171779821*fR[11]-1.060660171779821*fL[11]-1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6363961030678926*fR[6]-0.6363961030678926*fL[6]+0.8215838362577489*(fR[5]+fL[5])-0.8215838362577489*(fR[4]+fL[4])-0.4743416490252568*fR[3]+0.4743416490252568*(fL[3]+fR[2])-0.4743416490252568*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.7348469228349525*fR[19]-0.7348469228349525*fL[19]-0.7348469228349533*fR[18]+0.7348469228349533*fL[18]-1.42302494707577*(fR[17]+fL[17])-0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]+0.4242640687119285*(fR[14]+fL[14])+1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]-1.060660171779821*(fR[11]+fL[11])+1.10227038425243*fR[10]-1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6363961030678926*(fR[6]+fL[6])-0.8215838362577489*fR[5]+0.8215838362577489*(fL[5]+fR[4])-0.8215838362577489*fL[4]+0.4743416490252568*(fR[3]+fL[3])-0.4743416490252568*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.5590169943749475*alphaR[4]; 
  fUpOrd[6] = 0.5*(((-0.6846531968814573*(fR[18]+fL[18]))+0.3952847075210473*fR[14]-0.3952847075210473*fL[14]-0.7905694150420948*fR[13]+0.7905694150420948*fL[13]-0.6846531968814574*(fR[12]+fL[12])+0.3952847075210473*fR[8]-0.3952847075210473*fL[8]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6123724356957944*(fR[5]+fL[5])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6846531968814573*fR[18]-0.6846531968814573*fL[18]-0.3952847075210473*(fR[14]+fL[14])+0.7905694150420948*(fR[13]+fL[13])+0.6846531968814574*fR[12]-0.6846531968814574*fL[12]-0.3952847075210473*(fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6123724356957944*fR[5]+0.6123724356957944*fL[5]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.4472135954999579*alphaR[4]+0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
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

  incr[0] = 0.3535533905932737*(alphaR[4]*fUp[4]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = -0.6123724356957944*(alphaR[4]*fUp[4]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[2] = 0.07071067811865474*(4.47213595499958*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4])+5.0*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[3] = 0.02357022603955158*(15.0*alphaR[4]*fUp[6]+15.0*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2])); 
  incr[4] = -0.07071067811865474*(7.745966692414834*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4])+8.660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[5] = -0.07071067811865474*(8.660254037844387*alphaR[4]*fUp[6]+8.660254037844386*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2])); 
  incr[6] = 0.02357022603955158*(13.41640786499874*alphaR[1]*fUp[6]+13.41640786499874*fUp[3]*alphaR[4]+15.0*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[7] = 0.7905694150420947*(alphaR[4]*fUp[4]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[8] = 0.01010152544552211*((22.3606797749979*alphaR[4]+35.0*alphaR[0])*fUp[4]+35.0*fUp[0]*alphaR[4]+31.30495168499706*alphaR[1]*fUp[1]); 
  incr[9] = 0.02357022603955158*(15.0*alphaR[1]*fUp[7]+15.0*alphaR[0]*fUp[5]); 
  incr[10] = -0.07071067811865474*(7.745966692414834*(alphaR[1]*fUp[6]+fUp[3]*alphaR[4])+8.660254037844386*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[11] = 0.02357022603955158*(30.0*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4])+33.54101966249684*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[12] = -0.01010152544552211*((38.72983346207417*alphaR[4]+60.62177826491071*alphaR[0])*fUp[4]+60.62177826491071*fUp[0]*alphaR[4]+54.22176684690384*alphaR[1]*fUp[1]); 
  incr[13] = 0.1178511301977579*(6.708203932499369*alphaR[4]*fUp[6]+6.708203932499369*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2])); 
  incr[14] = 0.003367175148507369*((67.0820393249937*alphaR[4]+105.0*alphaR[0])*fUp[6]+105.0*fUp[2]*alphaR[4]+93.91485505499116*alphaR[1]*fUp[3]); 
  incr[15] = -0.07071067811865474*(8.660254037844386*alphaR[1]*fUp[7]+8.660254037844387*alphaR[0]*fUp[5]); 
  incr[16] = 0.02357022603955158*((13.41640786499874*alphaR[4]+15.0*alphaR[0])*fUp[7]+15.0*alphaR[1]*fUp[5]); 
  incr[17] = 0.02357022603955158*(30.0*alphaR[1]*fUp[6]+30.0*fUp[3]*alphaR[4]+33.54101966249685*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[18] = -0.01010152544552211*((38.72983346207417*alphaR[4]+60.62177826491071*alphaR[0])*fUp[6]+60.6217782649107*fUp[2]*alphaR[4]+54.22176684690384*alphaR[1]*fUp[3]); 
  incr[19] = -0.07071067811865474*((7.745966692414834*alphaR[4]+8.660254037844387*alphaR[0])*fUp[7]+8.660254037844386*alphaR[1]*fUp[5]); 

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
double GyrokineticGenGeoSurf1x2vSer_vpar_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R*wvparR+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (0.2828427124746191*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2R*wvparR+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*(b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[7] = (0.04040610178208843*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R*wvparR+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmagR[11] = (1.414213562373095*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = (0.25*((6.708203932499369*BstarZdBmagR[4]-3.872983346207417*BstarZdBmagR[1])*hamilR[7]+hamilR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]))*rdx2R)/m_; 
  alphaR[1] = (0.05*(hamilR[7]*(30.0*BstarZdBmagR[11]-17.32050807568877*BstarZdBmagR[7]+33.54101966249685*BstarZdBmagR[2]-19.36491673103709*BstarZdBmagR[0])+hamilR[1]*(15.0*BstarZdBmagR[4]-8.660254037844386*BstarZdBmagR[1]))*rdx2R)/m_; 
  alphaR[2] = (0.25*((6.708203932499369*BstarZdBmagR[4]-3.872983346207417*BstarZdBmagR[1])*hamilR[13]+(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*hamilR[5])*rdx2R)/m_; 
  alphaR[3] = (0.05*((30.0*BstarZdBmagR[11]-17.32050807568877*BstarZdBmagR[7]+33.54101966249684*BstarZdBmagR[2]-19.36491673103708*BstarZdBmagR[0])*hamilR[13]+(15.0*BstarZdBmagR[4]-8.660254037844386*BstarZdBmagR[1])*hamilR[5])*rdx2R)/m_; 
  alphaR[4] = (0.05*(15.0*hamilR[1]*BstarZdBmagR[11]+(30.0*BstarZdBmagR[4]-17.32050807568877*BstarZdBmagR[1])*hamilR[7]-8.660254037844386*hamilR[1]*BstarZdBmagR[7])*rdx2R)/m_; 
  alphaR[6] = (0.05*((30.0*BstarZdBmagR[4]-17.32050807568877*BstarZdBmagR[1])*hamilR[13]+hamilR[5]*(15.0*BstarZdBmagR[11]-8.660254037844387*BstarZdBmagR[7]))*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*((6.708203932499369*BstarZdBmagR[4]-3.872983346207417*BstarZdBmagR[1])*hamilR[7]+3.0*hamilR[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamilR[1])*rdx2R)/m_; 

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
double GyrokineticGenGeoSurf1x2vSer_x_P2_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

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
double GyrokineticGenGeoSurf1x2vSer_vpar_P2_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double alphaR[8]; 
  alphaR[0] = -(0.4330127018922193*BstarZdBmagR[0]*hamilR[1]*rdx2R)/m_; 
  alphaR[1] = -(0.9682458365518543*BstarZdBmagR[0]*hamilR[7]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.1082531754730548*BstarZdBmagR[0]*hamilR[1]*rdx2R)/m_; 

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
  Ghat[0] = -0.1767766952966368*((4.47213595499958*fR[8]-4.47213595499958*fL[8]-3.464101615137754*(fR[2]+fL[2])+2.0*fR[0]-2.0*fL[0])*amax-1.414213562373095*(alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.01178511301977579*((67.08203932499369*fR[12]-67.08203932499369*fL[12]-51.96152422706631*(fR[4]+fL[4])+30.0*fR[1]-30.0*fL[1])*amax-18.97366596101028*alphaR[1]*fAvg[4]-21.21320343559643*(alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.05892556509887893*((13.41640786499874*fR[14]-13.41640786499874*fL[14]-10.39230484541326*(fR[6]+fL[6])+6.0*fR[3]-6.0*fL[3])*amax-4.242640687119286*(alphaR[1]*fAvg[3]+alphaR[0]*fAvg[2])); 
  Ghat[3] = -0.01178511301977579*((67.0820393249937*fR[18]-67.0820393249937*fL[18]-51.96152422706631*(fR[10]+fL[10])+30.0*fR[5]-30.0*fL[5])*amax-18.97366596101028*alphaR[1]*fAvg[6]-21.21320343559643*(alphaR[0]*fAvg[3]+alphaR[1]*fAvg[2])); 
  Ghat[4] = 0.03535533905932736*((17.32050807568877*(fR[11]+fL[11])-10.0*fR[7]+10.0*fL[7])*amax+7.071067811865476*alphaR[0]*fAvg[4]+6.324555320336761*alphaR[1]*fAvg[1]); 
  Ghat[5] = 0.01178511301977579*((51.96152422706632*(fR[16]+fL[16])-30.0*fR[9]+30.0*fL[9])*amax+21.21320343559643*alphaR[1]*fAvg[7]+21.21320343559643*alphaR[0]*fAvg[5]); 
  Ghat[6] = 0.01178511301977579*((51.96152422706632*(fR[17]+fL[17])-30.0*fR[13]+30.0*fL[13])*amax+21.21320343559643*alphaR[0]*fAvg[6]+18.97366596101028*alphaR[1]*fAvg[3]); 
  Ghat[7] = 0.01178511301977579*((51.96152422706632*(fR[19]+fL[19])-30.0*fR[15]+30.0*fL[15])*amax+21.21320343559643*alphaR[0]*fAvg[7]+21.21320343559643*alphaR[1]*fAvg[5]); 

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
double GyrokineticGenGeoSurf1x2vSer_x_P2_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R*wvparR+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (0.2828427124746191*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2R*wvparR+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*(b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[7] = (0.04040610178208843*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R*wvparR+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmagR[11] = (1.414213562373095*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = (0.25*(hamilR[8]*(8.660254037844387*BstarZdBmagR[11]-6.708203932499369*BstarZdBmagR[4]+3.872983346207417*BstarZdBmagR[2])+hamilR[2]*(3.872983346207417*BstarZdBmagR[7]-3.0*BstarZdBmagR[1]+1.732050807568877*BstarZdBmagR[0]))*rdvpar2R)/m_; 
  alphaR[1] = (0.25*(3.872983346207417*hamilR[2]*BstarZdBmagR[11]+(8.660254037844386*BstarZdBmagR[7]-6.708203932499369*BstarZdBmagR[1]+3.872983346207417*BstarZdBmagR[0])*hamilR[8]+hamilR[2]*(1.732050807568877*BstarZdBmagR[2]-3.0*BstarZdBmagR[4]))*rdvpar2R)/m_; 
  alphaR[4] = (0.5*hamilR[8]*(3.872983346207417*BstarZdBmagR[11]-3.0*BstarZdBmagR[4]+1.732050807568877*BstarZdBmagR[2])*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(8.660254037844387*hamilR[8]*BstarZdBmagR[11]+(3.872983346207417*BstarZdBmagR[2]-6.708203932499369*BstarZdBmagR[4])*hamilR[8]+3.872983346207417*hamilR[2]*BstarZdBmagR[7]+(1.732050807568877*BstarZdBmagR[0]-3.0*BstarZdBmagR[1])*hamilR[2])*rdvpar2R)/m_; 

  double incr[20]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.01666666666666667*(25.98076211353316*alphaR[4]*fL[12]+33.54101966249684*alphaR[1]*fL[11]+15.0*alphaR[4]*fL[8]+33.54101966249685*alphaR[0]*fL[7]+alphaR[1]*(25.98076211353316*fL[4]+15.0*fL[2])+alphaR[0]*(25.98076211353316*fL[1]+15.0*fL[0])); 
  incr[1] = -0.05*(15.0*alphaR[4]*fL[12]+19.36491673103708*alphaR[1]*fL[11]+8.660254037844386*alphaR[4]*fL[8]+19.36491673103709*alphaR[0]*fL[7]+alphaR[1]*(15.0*fL[4]+8.660254037844386*fL[2])+alphaR[0]*(15.0*fL[1]+8.660254037844386*fL[0])); 
  incr[2] = 0.01666666666666667*(23.2379000772445*alphaR[1]*fL[12]+(30.0*alphaR[4]+33.54101966249684*alphaR[0])*fL[11]+alphaR[1]*(13.41640786499874*fL[8]+33.54101966249685*fL[7])+(23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fL[4]+fL[2]*(13.41640786499874*alphaR[4]+15.0*alphaR[0])+alphaR[1]*(25.98076211353316*fL[1]+15.0*fL[0])); 
  incr[3] = 0.01666666666666667*(25.98076211353316*alphaR[4]*fL[18]+33.54101966249685*alphaR[1]*fL[17]+15.0*alphaR[4]*fL[14]+33.54101966249684*alphaR[0]*fL[13]+alphaR[1]*(25.98076211353316*fL[10]+15.0*fL[6])+alphaR[0]*(25.98076211353316*fL[5]+15.0*fL[3])); 
  incr[4] = -0.05*(13.41640786499874*alphaR[1]*fL[12]+(17.32050807568877*alphaR[4]+19.36491673103708*alphaR[0])*fL[11]+alphaR[1]*(7.745966692414834*fL[8]+19.36491673103709*fL[7])+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fL[4]+fL[2]*(7.745966692414834*alphaR[4]+8.660254037844386*alphaR[0])+alphaR[1]*(15.0*fL[1]+8.660254037844386*fL[0])); 
  incr[5] = -0.05*(15.0*alphaR[4]*fL[18]+19.36491673103709*alphaR[1]*fL[17]+8.660254037844387*alphaR[4]*fL[14]+19.36491673103708*alphaR[0]*fL[13]+alphaR[1]*(15.0*fL[10]+8.660254037844386*fL[6])+alphaR[0]*(15.0*fL[5]+8.660254037844386*fL[3])); 
  incr[6] = 0.01666666666666667*(23.2379000772445*alphaR[1]*fL[18]+(30.0*alphaR[4]+33.54101966249685*alphaR[0])*fL[17]+alphaR[1]*(13.41640786499874*fL[14]+33.54101966249684*fL[13])+(23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fL[10]+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fL[6]+alphaR[1]*(25.98076211353316*fL[5]+15.0*fL[3])); 
  incr[7] = 0.08333333333333333*(11.61895003862225*alphaR[4]*fL[12]+15.0*alphaR[1]*fL[11]+6.708203932499369*alphaR[4]*fL[8]+15.0*alphaR[0]*fL[7]+alphaR[1]*(11.61895003862225*fL[4]+6.708203932499369*fL[2])+alphaR[0]*(11.61895003862225*fL[1]+6.708203932499369*fL[0])); 
  incr[8] = 0.002380952380952381*((116.1895003862225*alphaR[4]+181.8653347947321*alphaR[0])*fL[12]+210.0*alphaR[1]*fL[11]+(67.0820393249937*alphaR[4]+105.0*alphaR[0])*fL[8]+234.787137637478*alphaR[4]*fL[7]+162.6653005407115*alphaR[1]*fL[4]+(181.8653347947321*fL[1]+105.0*fL[0])*alphaR[4]+93.91485505499116*alphaR[1]*fL[2]); 
  incr[9] = 0.01666666666666667*(alphaR[1]*(25.98076211353316*fL[19]+15.0*fL[16])+alphaR[0]*(25.98076211353316*fL[15]+15.0*fL[9])); 
  incr[10] = -0.05*(13.41640786499874*alphaR[1]*fL[18]+(17.32050807568877*alphaR[4]+19.36491673103709*alphaR[0])*fL[17]+alphaR[1]*(7.745966692414834*fL[14]+19.36491673103708*fL[13])+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fL[10]+(7.745966692414834*alphaR[4]+8.660254037844386*alphaR[0])*fL[6]+alphaR[1]*(15.0*fL[5]+8.660254037844386*fL[3])); 
  incr[11] = 0.01666666666666667*(51.96152422706631*alphaR[1]*fL[12]+(67.0820393249937*alphaR[4]+75.0*alphaR[0])*fL[11]+alphaR[1]*(30.0*fL[8]+75.00000000000001*fL[7])+(51.96152422706632*alphaR[4]+58.09475019311126*alphaR[0])*fL[4]+fL[2]*(30.0*alphaR[4]+33.54101966249684*alphaR[0])+alphaR[1]*(58.09475019311126*fL[1]+33.54101966249684*fL[0])); 
  incr[12] = -0.007142857142857143*((67.0820393249937*alphaR[4]+105.0*alphaR[0])*fL[12]+121.2435565298214*alphaR[1]*fL[11]+(38.72983346207417*alphaR[4]+60.62177826491071*alphaR[0])*fL[8]+135.5544171172596*alphaR[4]*fL[7]+93.91485505499116*alphaR[1]*fL[4]+(105.0*fL[1]+60.62177826491071*fL[0])*alphaR[4]+54.22176684690384*alphaR[1]*fL[2]); 
  incr[13] = 0.08333333333333333*(11.61895003862225*alphaR[4]*fL[18]+15.0*alphaR[1]*fL[17]+6.708203932499369*alphaR[4]*fL[14]+15.0*alphaR[0]*fL[13]+alphaR[1]*(11.61895003862225*fL[10]+6.708203932499369*fL[6])+alphaR[0]*(11.61895003862225*fL[5]+6.708203932499369*fL[3])); 
  incr[14] = 0.002380952380952381*((116.1895003862225*alphaR[4]+181.8653347947321*alphaR[0])*fL[18]+210.0*alphaR[1]*fL[17]+(67.0820393249937*alphaR[4]+105.0*alphaR[0])*fL[14]+234.787137637478*alphaR[4]*fL[13]+alphaR[1]*(162.6653005407115*fL[10]+93.91485505499116*fL[6])+alphaR[4]*(181.8653347947321*fL[5]+105.0*fL[3])); 
  incr[15] = -0.05*(alphaR[1]*(15.0*fL[19]+8.660254037844386*fL[16])+alphaR[0]*(15.0*fL[15]+8.660254037844387*fL[9])); 
  incr[16] = 0.01666666666666667*((23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fL[19]+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fL[16]+alphaR[1]*(25.98076211353316*fL[15]+15.0*fL[9])); 
  incr[17] = 0.01666666666666667*(51.96152422706631*alphaR[1]*fL[18]+(67.0820393249937*alphaR[4]+75.0*alphaR[0])*fL[17]+alphaR[1]*(30.0*fL[14]+75.00000000000001*fL[13])+(51.96152422706631*alphaR[4]+58.09475019311126*alphaR[0])*fL[10]+(30.0*alphaR[4]+33.54101966249685*alphaR[0])*fL[6]+alphaR[1]*(58.09475019311126*fL[5]+33.54101966249685*fL[3])); 
  incr[18] = -0.007142857142857143*((67.0820393249937*alphaR[4]+105.0*alphaR[0])*fL[18]+121.2435565298214*alphaR[1]*fL[17]+(38.72983346207417*alphaR[4]+60.62177826491071*alphaR[0])*fL[14]+135.5544171172596*alphaR[4]*fL[13]+alphaR[1]*(93.91485505499116*fL[10]+54.22176684690384*fL[6])+alphaR[4]*(105.0*fL[5]+60.6217782649107*fL[3])); 
  incr[19] = -0.05*((13.41640786499874*alphaR[4]+15.0*alphaR[0])*fL[19]+(7.745966692414834*alphaR[4]+8.660254037844387*alphaR[0])*fL[16]+alphaR[1]*(15.0*fL[15]+8.660254037844386*fL[9])); 
  } else { 
  incr[0] = -0.01666666666666667*(25.98076211353316*alphaR[4]*fR[12]-33.54101966249684*alphaR[1]*fR[11]-15.0*alphaR[4]*fR[8]-33.54101966249685*alphaR[0]*fR[7]+alphaR[1]*(25.98076211353316*fR[4]-15.0*fR[2])+alphaR[0]*(25.98076211353316*fR[1]-15.0*fR[0])); 
  incr[1] = 0.05*(15.0*alphaR[4]*fR[12]-19.36491673103708*alphaR[1]*fR[11]-8.660254037844386*alphaR[4]*fR[8]-19.36491673103709*alphaR[0]*fR[7]+alphaR[1]*(15.0*fR[4]-8.660254037844386*fR[2])+alphaR[0]*(15.0*fR[1]-8.660254037844386*fR[0])); 
  incr[2] = -0.01666666666666667*(23.2379000772445*alphaR[1]*fR[12]+((-30.0*alphaR[4])-33.54101966249684*alphaR[0])*fR[11]+alphaR[1]*((-13.41640786499874*fR[8])-33.54101966249685*fR[7])+(23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fR[4]+fR[2]*((-13.41640786499874*alphaR[4])-15.0*alphaR[0])+alphaR[1]*(25.98076211353316*fR[1]-15.0*fR[0])); 
  incr[3] = -0.01666666666666667*(25.98076211353316*alphaR[4]*fR[18]-33.54101966249685*alphaR[1]*fR[17]-15.0*alphaR[4]*fR[14]-33.54101966249684*alphaR[0]*fR[13]+alphaR[1]*(25.98076211353316*fR[10]-15.0*fR[6])+alphaR[0]*(25.98076211353316*fR[5]-15.0*fR[3])); 
  incr[4] = 0.05*(13.41640786499874*alphaR[1]*fR[12]+((-17.32050807568877*alphaR[4])-19.36491673103708*alphaR[0])*fR[11]+alphaR[1]*((-7.745966692414834*fR[8])-19.36491673103709*fR[7])+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fR[4]+fR[2]*((-7.745966692414834*alphaR[4])-8.660254037844386*alphaR[0])+alphaR[1]*(15.0*fR[1]-8.660254037844386*fR[0])); 
  incr[5] = 0.05*(15.0*alphaR[4]*fR[18]-19.36491673103709*alphaR[1]*fR[17]-8.660254037844387*alphaR[4]*fR[14]-19.36491673103708*alphaR[0]*fR[13]+alphaR[1]*(15.0*fR[10]-8.660254037844386*fR[6])+alphaR[0]*(15.0*fR[5]-8.660254037844386*fR[3])); 
  incr[6] = -0.01666666666666667*(23.2379000772445*alphaR[1]*fR[18]+((-30.0*alphaR[4])-33.54101966249685*alphaR[0])*fR[17]+alphaR[1]*((-13.41640786499874*fR[14])-33.54101966249684*fR[13])+(23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fR[10]+((-13.41640786499874*alphaR[4])-15.0*alphaR[0])*fR[6]+alphaR[1]*(25.98076211353316*fR[5]-15.0*fR[3])); 
  incr[7] = -0.08333333333333333*(11.61895003862225*alphaR[4]*fR[12]-15.0*alphaR[1]*fR[11]-6.708203932499369*alphaR[4]*fR[8]-15.0*alphaR[0]*fR[7]+alphaR[1]*(11.61895003862225*fR[4]-6.708203932499369*fR[2])+alphaR[0]*(11.61895003862225*fR[1]-6.708203932499369*fR[0])); 
  incr[8] = -0.002380952380952381*((116.1895003862225*alphaR[4]+181.8653347947321*alphaR[0])*fR[12]-210.0*alphaR[1]*fR[11]+((-67.0820393249937*alphaR[4])-105.0*alphaR[0])*fR[8]-234.787137637478*alphaR[4]*fR[7]+162.6653005407115*alphaR[1]*fR[4]+(181.8653347947321*fR[1]-105.0*fR[0])*alphaR[4]-93.91485505499116*alphaR[1]*fR[2]); 
  incr[9] = -0.01666666666666667*(alphaR[1]*(25.98076211353316*fR[19]-15.0*fR[16])+alphaR[0]*(25.98076211353316*fR[15]-15.0*fR[9])); 
  incr[10] = 0.05*(13.41640786499874*alphaR[1]*fR[18]+((-17.32050807568877*alphaR[4])-19.36491673103709*alphaR[0])*fR[17]+alphaR[1]*((-7.745966692414834*fR[14])-19.36491673103708*fR[13])+(13.41640786499874*alphaR[4]+15.0*alphaR[0])*fR[10]+((-7.745966692414834*alphaR[4])-8.660254037844386*alphaR[0])*fR[6]+alphaR[1]*(15.0*fR[5]-8.660254037844386*fR[3])); 
  incr[11] = -0.01666666666666667*(51.96152422706631*alphaR[1]*fR[12]+((-67.0820393249937*alphaR[4])-75.0*alphaR[0])*fR[11]+alphaR[1]*((-30.0*fR[8])-75.00000000000001*fR[7])+(51.96152422706632*alphaR[4]+58.09475019311126*alphaR[0])*fR[4]+fR[2]*((-30.0*alphaR[4])-33.54101966249684*alphaR[0])+alphaR[1]*(58.09475019311126*fR[1]-33.54101966249684*fR[0])); 
  incr[12] = 0.007142857142857143*((67.0820393249937*alphaR[4]+105.0*alphaR[0])*fR[12]-121.2435565298214*alphaR[1]*fR[11]+((-38.72983346207417*alphaR[4])-60.62177826491071*alphaR[0])*fR[8]-135.5544171172596*alphaR[4]*fR[7]+93.91485505499116*alphaR[1]*fR[4]+(105.0*fR[1]-60.62177826491071*fR[0])*alphaR[4]-54.22176684690384*alphaR[1]*fR[2]); 
  incr[13] = -0.08333333333333333*(11.61895003862225*alphaR[4]*fR[18]-15.0*alphaR[1]*fR[17]-6.708203932499369*alphaR[4]*fR[14]-15.0*alphaR[0]*fR[13]+alphaR[1]*(11.61895003862225*fR[10]-6.708203932499369*fR[6])+alphaR[0]*(11.61895003862225*fR[5]-6.708203932499369*fR[3])); 
  incr[14] = -0.002380952380952381*((116.1895003862225*alphaR[4]+181.8653347947321*alphaR[0])*fR[18]-210.0*alphaR[1]*fR[17]+((-67.0820393249937*alphaR[4])-105.0*alphaR[0])*fR[14]-234.787137637478*alphaR[4]*fR[13]+alphaR[1]*(162.6653005407115*fR[10]-93.91485505499116*fR[6])+alphaR[4]*(181.8653347947321*fR[5]-105.0*fR[3])); 
  incr[15] = 0.05*(alphaR[1]*(15.0*fR[19]-8.660254037844386*fR[16])+alphaR[0]*(15.0*fR[15]-8.660254037844387*fR[9])); 
  incr[16] = -0.01666666666666667*((23.2379000772445*alphaR[4]+25.98076211353316*alphaR[0])*fR[19]+((-13.41640786499874*alphaR[4])-15.0*alphaR[0])*fR[16]+alphaR[1]*(25.98076211353316*fR[15]-15.0*fR[9])); 
  incr[17] = -0.01666666666666667*(51.96152422706631*alphaR[1]*fR[18]+((-67.0820393249937*alphaR[4])-75.0*alphaR[0])*fR[17]+alphaR[1]*((-30.0*fR[14])-75.00000000000001*fR[13])+(51.96152422706631*alphaR[4]+58.09475019311126*alphaR[0])*fR[10]+((-30.0*alphaR[4])-33.54101966249685*alphaR[0])*fR[6]+alphaR[1]*(58.09475019311126*fR[5]-33.54101966249685*fR[3])); 
  incr[18] = 0.007142857142857143*((67.0820393249937*alphaR[4]+105.0*alphaR[0])*fR[18]-121.2435565298214*alphaR[1]*fR[17]+((-38.72983346207417*alphaR[4])-60.62177826491071*alphaR[0])*fR[14]-135.5544171172596*alphaR[4]*fR[13]+alphaR[1]*(93.91485505499116*fR[10]-54.22176684690384*fR[6])+alphaR[4]*(105.0*fR[5]-60.6217782649107*fR[3])); 
  incr[19] = 0.05*((13.41640786499874*alphaR[4]+15.0*alphaR[0])*fR[19]+((-7.745966692414834*alphaR[4])-8.660254037844387*alphaR[0])*fR[16]+alphaR[1]*(15.0*fR[15]-8.660254037844386*fR[9])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.4472135954999579*alphaR[4]-0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[0] = 0.5*(((-0.7348469228349525*(fR[19]+fL[19]))-0.7348469228349533*(fR[18]+fL[18])-1.42302494707577*fR[17]+1.42302494707577*fL[17]+0.4242640687119281*fR[16]-0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])+0.4242640687119285*fR[14]-0.4242640687119285*fL[14]+1.060660171779821*fR[13]-1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])+1.060660171779821*fR[11]-1.060660171779821*fL[11]+1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6363961030678926*fR[6]+0.6363961030678926*fL[6]-0.8215838362577489*(fR[5]+fL[5]+fR[4]+fL[4])+0.4743416490252568*(fR[3]+fR[2])-0.4743416490252568*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.7348469228349525*fR[19]-0.7348469228349525*fL[19]+0.7348469228349533*fR[18]-0.7348469228349533*fL[18]+1.42302494707577*(fR[17]+fL[17])-0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]-0.4242640687119285*(fR[14]+fL[14])-1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]-1.060660171779821*(fR[11]+fL[11])-1.10227038425243*fR[10]+1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])+0.6363961030678926*(fR[6]+fL[6])+0.8215838362577489*(fR[5]+fR[4])-0.8215838362577489*(fL[5]+fL[4])-0.4743416490252568*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.5590169943749475*alphaR[4]; 
  fUpOrd[1] = 0.5*((0.6846531968814573*(fR[18]+fL[18])-0.3952847075210473*fR[14]+0.3952847075210473*fL[14]+0.7905694150420948*fR[13]-0.7905694150420948*fL[13]-0.6846531968814574*(fR[12]+fL[12])+0.3952847075210473*fR[8]-0.3952847075210473*fL[8]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6123724356957944*(fR[5]+fL[5])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6846531968814573*fR[18]+0.6846531968814573*fL[18]+0.3952847075210473*(fR[14]+fL[14])-0.7905694150420948*(fR[13]+fL[13])+0.6846531968814574*fR[12]-0.6846531968814574*fL[12]-0.3952847075210473*(fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])+0.6123724356957944*fR[5]-0.6123724356957944*fL[5]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.4472135954999579*alphaR[4]+0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[2] = 0.5*((0.7348469228349525*(fR[19]+fL[19])-0.7348469228349533*(fR[18]+fL[18])+1.42302494707577*fR[17]-1.42302494707577*fL[17]-0.4242640687119281*fR[16]+0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])+0.4242640687119285*fR[14]-0.4242640687119285*fL[14]+1.060660171779821*fR[13]-1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])-1.060660171779821*fR[11]+1.060660171779821*fL[11]-1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6363961030678926*fR[6]-0.6363961030678926*fL[6]-0.8215838362577489*(fR[5]+fL[5])+0.8215838362577489*(fR[4]+fL[4])+0.4743416490252568*fR[3]-0.4743416490252568*(fL[3]+fR[2])+0.4743416490252568*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.7348469228349525*fR[19]+0.7348469228349525*fL[19]+0.7348469228349533*fR[18]-0.7348469228349533*fL[18]-1.42302494707577*(fR[17]+fL[17])+0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]-0.4242640687119285*(fR[14]+fL[14])-1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]+1.060660171779821*(fR[11]+fL[11])+1.10227038425243*fR[10]-1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6363961030678926*(fR[6]+fL[6])+0.8215838362577489*fR[5]-0.8215838362577489*(fL[5]+fR[4])+0.8215838362577489*fL[4]-0.4743416490252568*(fR[3]+fL[3])+0.4743416490252568*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.5*alphaR[1]; 
  fUpOrd[3] = 0.5*((0.6846531968814573*(fR[19]+fL[19])-0.3952847075210473*fR[16]+0.3952847075210473*fL[16]-0.6846531968814574*(fR[15]+fL[15])+0.7905694150420948*fR[11]-0.7905694150420948*fL[11]+0.3952847075210473*fR[9]-0.3952847075210473*fL[9]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]-0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[2]-0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6846531968814573*fR[19]+0.6846531968814573*fL[19]+0.3952847075210473*(fR[16]+fL[16])+0.6846531968814574*fR[15]-0.6846531968814574*fL[15]-0.7905694150420948*(fR[11]+fL[11])-0.3952847075210473*(fR[9]+fL[9])+0.7905694150420947*(fR[7]+fL[7])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]-0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[1]+alphaR[0]); 
  fUpOrd[4] = 0.5*(((-0.6846531968814573*(fR[19]+fL[19]))+0.3952847075210473*fR[16]-0.3952847075210473*fL[16]-0.6846531968814574*(fR[15]+fL[15])-0.7905694150420948*fR[11]+0.7905694150420948*fL[11]+0.3952847075210473*fR[9]-0.3952847075210473*fL[9]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[2]+0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6846531968814573*fR[19]-0.6846531968814573*fL[19]-0.3952847075210473*(fR[16]+fL[16])+0.6846531968814574*fR[15]-0.6846531968814574*fL[15]+0.7905694150420948*(fR[11]+fL[11])-0.3952847075210473*(fR[9]+fL[9])+0.7905694150420947*(fR[7]+fL[7])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]+0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.4472135954999579*alphaR[4]-0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[5] = 0.5*(((-0.7348469228349525*(fR[19]+fL[19]))+0.7348469228349533*(fR[18]+fL[18])+1.42302494707577*fR[17]-1.42302494707577*fL[17]+0.4242640687119281*fR[16]-0.4242640687119281*fL[16]+0.5477225575051661*(fR[15]+fL[15])-0.4242640687119285*fR[14]+0.4242640687119285*fL[14]-1.060660171779821*fR[13]+1.060660171779821*fL[13]+0.5477225575051661*(fR[12]+fL[12])+1.060660171779821*fR[11]-1.060660171779821*fL[11]-1.10227038425243*(fR[10]+fL[10])-0.3162277660168379*(fR[9]+fR[8])+0.3162277660168379*(fL[9]+fL[8])-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6363961030678926*fR[6]-0.6363961030678926*fL[6]+0.8215838362577489*(fR[5]+fL[5])-0.8215838362577489*(fR[4]+fL[4])-0.4743416490252568*fR[3]+0.4743416490252568*(fL[3]+fR[2])-0.4743416490252568*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.7348469228349525*fR[19]-0.7348469228349525*fL[19]-0.7348469228349533*fR[18]+0.7348469228349533*fL[18]-1.42302494707577*(fR[17]+fL[17])-0.4242640687119281*(fR[16]+fL[16])-0.5477225575051661*fR[15]+0.5477225575051661*fL[15]+0.4242640687119285*(fR[14]+fL[14])+1.060660171779821*(fR[13]+fL[13])-0.5477225575051661*fR[12]+0.5477225575051661*fL[12]-1.060660171779821*(fR[11]+fL[11])+1.10227038425243*fR[10]-1.10227038425243*fL[10]+0.3162277660168379*(fR[9]+fL[9]+fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6363961030678926*(fR[6]+fL[6])-0.8215838362577489*fR[5]+0.8215838362577489*(fL[5]+fR[4])-0.8215838362577489*fL[4]+0.4743416490252568*(fR[3]+fL[3])-0.4743416490252568*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.5590169943749475*alphaR[4]; 
  fUpOrd[6] = 0.5*(((-0.6846531968814573*(fR[18]+fL[18]))+0.3952847075210473*fR[14]-0.3952847075210473*fL[14]-0.7905694150420948*fR[13]+0.7905694150420948*fL[13]-0.6846531968814574*(fR[12]+fL[12])+0.3952847075210473*fR[8]-0.3952847075210473*fL[8]-0.7905694150420947*fR[7]+0.7905694150420947*fL[7]+0.6123724356957944*(fR[5]+fL[5])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6846531968814573*fR[18]-0.6846531968814573*fL[18]-0.3952847075210473*(fR[14]+fL[14])+0.7905694150420948*(fR[13]+fL[13])+0.6846531968814574*fR[12]-0.6846531968814574*fL[12]-0.3952847075210473*(fR[8]+fL[8])+0.7905694150420947*(fR[7]+fL[7])-0.6123724356957944*fR[5]+0.6123724356957944*fL[5]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.4472135954999579*alphaR[4]+0.6708203932499369*alphaR[1]+0.5*alphaR[0]; 
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

  incr[0] = 0.3535533905932737*(alphaR[4]*fUp[4]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = -0.6123724356957944*(alphaR[4]*fUp[4]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[2] = 0.07071067811865474*(4.47213595499958*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4])+5.0*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[3] = 0.02357022603955158*(15.0*alphaR[4]*fUp[6]+15.0*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2])); 
  incr[4] = -0.07071067811865474*(7.745966692414834*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4])+8.660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[5] = -0.07071067811865474*(8.660254037844387*alphaR[4]*fUp[6]+8.660254037844386*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2])); 
  incr[6] = 0.02357022603955158*(13.41640786499874*alphaR[1]*fUp[6]+13.41640786499874*fUp[3]*alphaR[4]+15.0*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[7] = 0.7905694150420947*(alphaR[4]*fUp[4]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[8] = 0.01010152544552211*((22.3606797749979*alphaR[4]+35.0*alphaR[0])*fUp[4]+35.0*fUp[0]*alphaR[4]+31.30495168499706*alphaR[1]*fUp[1]); 
  incr[9] = 0.02357022603955158*(15.0*alphaR[1]*fUp[7]+15.0*alphaR[0]*fUp[5]); 
  incr[10] = -0.07071067811865474*(7.745966692414834*(alphaR[1]*fUp[6]+fUp[3]*alphaR[4])+8.660254037844386*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[11] = 0.02357022603955158*(30.0*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4])+33.54101966249684*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1])); 
  incr[12] = -0.01010152544552211*((38.72983346207417*alphaR[4]+60.62177826491071*alphaR[0])*fUp[4]+60.62177826491071*fUp[0]*alphaR[4]+54.22176684690384*alphaR[1]*fUp[1]); 
  incr[13] = 0.1178511301977579*(6.708203932499369*alphaR[4]*fUp[6]+6.708203932499369*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2])); 
  incr[14] = 0.003367175148507369*((67.0820393249937*alphaR[4]+105.0*alphaR[0])*fUp[6]+105.0*fUp[2]*alphaR[4]+93.91485505499116*alphaR[1]*fUp[3]); 
  incr[15] = -0.07071067811865474*(8.660254037844386*alphaR[1]*fUp[7]+8.660254037844387*alphaR[0]*fUp[5]); 
  incr[16] = 0.02357022603955158*((13.41640786499874*alphaR[4]+15.0*alphaR[0])*fUp[7]+15.0*alphaR[1]*fUp[5]); 
  incr[17] = 0.02357022603955158*(30.0*alphaR[1]*fUp[6]+30.0*fUp[3]*alphaR[4]+33.54101966249685*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2])); 
  incr[18] = -0.01010152544552211*((38.72983346207417*alphaR[4]+60.62177826491071*alphaR[0])*fUp[6]+60.6217782649107*fUp[2]*alphaR[4]+54.22176684690384*alphaR[1]*fUp[3]); 
  incr[19] = -0.07071067811865474*((7.745966692414834*alphaR[4]+8.660254037844387*alphaR[0])*fUp[7]+8.660254037844386*alphaR[1]*fUp[5]); 

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
double GyrokineticGenGeoSurf1x2vSer_vpar_P2_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
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
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R*wvparR+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (0.2828427124746191*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2R*wvparR+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*(b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[7] = (0.04040610178208843*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R*wvparR+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmagR[11] = (1.414213562373095*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = (0.25*((6.708203932499369*BstarZdBmagR[4]-3.872983346207417*BstarZdBmagR[1])*hamilR[7]+hamilR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]))*rdx2R)/m_; 
  alphaR[1] = (0.05*(hamilR[7]*(30.0*BstarZdBmagR[11]-17.32050807568877*BstarZdBmagR[7]+33.54101966249685*BstarZdBmagR[2]-19.36491673103709*BstarZdBmagR[0])+hamilR[1]*(15.0*BstarZdBmagR[4]-8.660254037844386*BstarZdBmagR[1]))*rdx2R)/m_; 
  alphaR[2] = (0.25*((6.708203932499369*BstarZdBmagR[4]-3.872983346207417*BstarZdBmagR[1])*hamilR[13]+(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*hamilR[5])*rdx2R)/m_; 
  alphaR[3] = (0.05*((30.0*BstarZdBmagR[11]-17.32050807568877*BstarZdBmagR[7]+33.54101966249684*BstarZdBmagR[2]-19.36491673103708*BstarZdBmagR[0])*hamilR[13]+(15.0*BstarZdBmagR[4]-8.660254037844386*BstarZdBmagR[1])*hamilR[5])*rdx2R)/m_; 
  alphaR[4] = (0.05*(15.0*hamilR[1]*BstarZdBmagR[11]+(30.0*BstarZdBmagR[4]-17.32050807568877*BstarZdBmagR[1])*hamilR[7]-8.660254037844386*hamilR[1]*BstarZdBmagR[7])*rdx2R)/m_; 
  alphaR[6] = (0.05*((30.0*BstarZdBmagR[4]-17.32050807568877*BstarZdBmagR[1])*hamilR[13]+hamilR[5]*(15.0*BstarZdBmagR[11]-8.660254037844387*BstarZdBmagR[7]))*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*((6.708203932499369*BstarZdBmagR[4]-3.872983346207417*BstarZdBmagR[1])*hamilR[7]+3.0*hamilR[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamilR[1])*rdx2R)/m_; 

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
