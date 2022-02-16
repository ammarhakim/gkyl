#include <VlasovModDecl.h> 
__host__ __device__ double VlasovUpwindSurfNeutral1x2vSer_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double *Fo0 = &EM[0]; 

  double alpha[8]; 
  double incr[20]; 

  alpha[0] = 1.414213562373095*Fo0[0]; 
  alpha[1] = 1.414213562373095*Fo0[1]; 
  alpha[4] = 1.414213562373095*Fo0[2]; 

  const double amid = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 

  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0]; 
  fUpOrd[0] = 0.5*(((-0.7348469228349525*(fr[19]+fl[19]))-1.42302494707577*fr[18]+1.42302494707577*fl[18]-0.7348469228349533*(fr[17]+fl[17])+0.5477225575051661*(fr[16]+fl[16])+0.4242640687119281*fr[15]-0.4242640687119281*fl[15]+1.060660171779821*fr[14]-1.060660171779821*fl[14]+0.4242640687119285*fr[13]-0.4242640687119285*fl[13]+1.060660171779821*fr[12]-1.060660171779821*fl[12]+0.5477225575051661*(fr[11]+fl[11])+1.10227038425243*(fr[10]+fl[10])-0.3162277660168379*fr[9]+0.3162277660168379*fl[9]-0.7905694150420947*fr[8]+0.7905694150420947*fl[8]-0.3162277660168379*fr[7]+0.3162277660168379*fl[7]-0.8215838362577489*(fr[6]+fl[6])-0.6363961030678926*fr[5]+0.6363961030678926*fl[5]-0.8215838362577489*(fr[4]+fl[4])+0.4743416490252568*fr[3]-0.4743416490252568*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.4743416490252568*fr[1]-0.4743416490252568*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaOrdR)+0.7348469228349525*fr[19]-0.7348469228349525*fl[19]+1.42302494707577*(fr[18]+fl[18])+0.7348469228349533*fr[17]-0.7348469228349533*fl[17]-0.5477225575051661*fr[16]+0.5477225575051661*fl[16]-0.4242640687119281*(fr[15]+fl[15])-1.060660171779821*(fr[14]+fl[14])-0.4242640687119285*(fr[13]+fl[13])-1.060660171779821*(fr[12]+fl[12])-0.5477225575051661*fr[11]+0.5477225575051661*fl[11]-1.10227038425243*fr[10]+1.10227038425243*fl[10]+0.3162277660168379*(fr[9]+fl[9])+0.7905694150420947*(fr[8]+fl[8])+0.3162277660168379*(fr[7]+fl[7])+0.8215838362577489*fr[6]-0.8215838362577489*fl[6]+0.6363961030678926*(fr[5]+fl[5])+0.8215838362577489*fr[4]-0.8215838362577489*fl[4]-0.4743416490252568*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.4743416490252568*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 
  fUpOrd[1] = 0.5*((0.6846531968814573*(fr[17]+fl[17])+0.7905694150420948*fr[14]-0.7905694150420948*fl[14]-0.3952847075210473*fr[13]+0.3952847075210473*fl[13]-0.6846531968814574*(fr[11]+fl[11])-0.7905694150420947*fr[8]+0.7905694150420947*fl[8]+0.3952847075210473*fr[7]-0.3952847075210473*fl[7]-0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaOrdR)-0.6846531968814573*fr[17]+0.6846531968814573*fl[17]-0.7905694150420948*(fr[14]+fl[14])+0.3952847075210473*(fr[13]+fl[13])+0.6846531968814574*fr[11]-0.6846531968814574*fl[11]+0.7905694150420947*(fr[8]+fl[8])-0.3952847075210473*(fr[7]+fl[7])+0.6123724356957944*fr[6]-0.6123724356957944*fl[6]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0]; 
  fUpOrd[2] = 0.5*((0.7348469228349525*(fr[19]+fl[19])+1.42302494707577*fr[18]-1.42302494707577*fl[18]-0.7348469228349533*(fr[17]+fl[17])+0.5477225575051661*(fr[16]+fl[16])-0.4242640687119281*fr[15]+0.4242640687119281*fl[15]+1.060660171779821*fr[14]-1.060660171779821*fl[14]+0.4242640687119285*fr[13]-0.4242640687119285*fl[13]-1.060660171779821*fr[12]+1.060660171779821*fl[12]+0.5477225575051661*(fr[11]+fl[11])-1.10227038425243*(fr[10]+fl[10])-0.3162277660168379*fr[9]+0.3162277660168379*fl[9]-0.7905694150420947*fr[8]+0.7905694150420947*fl[8]-0.3162277660168379*fr[7]+0.3162277660168379*fl[7]-0.8215838362577489*(fr[6]+fl[6])+0.6363961030678926*fr[5]-0.6363961030678926*fl[5]+0.8215838362577489*(fr[4]+fl[4])+0.4743416490252568*fr[3]-0.4743416490252568*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.4743416490252568*fr[1]+0.4743416490252568*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaOrdR)-0.7348469228349525*fr[19]+0.7348469228349525*fl[19]-1.42302494707577*(fr[18]+fl[18])+0.7348469228349533*fr[17]-0.7348469228349533*fl[17]-0.5477225575051661*fr[16]+0.5477225575051661*fl[16]+0.4242640687119281*(fr[15]+fl[15])-1.060660171779821*(fr[14]+fl[14])-0.4242640687119285*(fr[13]+fl[13])+1.060660171779821*(fr[12]+fl[12])-0.5477225575051661*fr[11]+0.5477225575051661*fl[11]+1.10227038425243*fr[10]-1.10227038425243*fl[10]+0.3162277660168379*(fr[9]+fl[9])+0.7905694150420947*(fr[8]+fl[8])+0.3162277660168379*(fr[7]+fl[7])+0.8215838362577489*fr[6]-0.8215838362577489*fl[6]-0.6363961030678926*(fr[5]+fl[5])-0.8215838362577489*fr[4]+0.8215838362577489*fl[4]-0.4743416490252568*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.4743416490252568*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*alpha[0]-0.5*alpha[1]; 
  fUpOrd[3] = 0.5*((0.6846531968814573*(fr[19]+fl[19])-0.6846531968814574*(fr[16]+fl[16])-0.3952847075210473*fr[15]+0.3952847075210473*fl[15]+0.7905694150420948*fr[12]-0.7905694150420948*fl[12]+0.3952847075210473*fr[9]-0.3952847075210473*fl[9]-0.7905694150420947*fr[8]+0.7905694150420947*fl[8]-0.6123724356957944*(fr[4]+fl[4])+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaOrdR)-0.6846531968814573*fr[19]+0.6846531968814573*fl[19]+0.6846531968814574*fr[16]-0.6846531968814574*fl[16]+0.3952847075210473*(fr[15]+fl[15])-0.7905694150420948*(fr[12]+fl[12])-0.3952847075210473*(fr[9]+fl[9])+0.7905694150420947*(fr[8]+fl[8])+0.6123724356957944*fr[4]-0.6123724356957944*(fl[4]+fr[2])+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*(alpha[1]+alpha[0]); 
  fUpOrd[4] = 0.5*(((-0.6846531968814573*(fr[19]+fl[19]))-0.6846531968814574*(fr[16]+fl[16])+0.3952847075210473*fr[15]-0.3952847075210473*fl[15]-0.7905694150420948*fr[12]+0.7905694150420948*fl[12]+0.3952847075210473*fr[9]-0.3952847075210473*fl[9]-0.7905694150420947*fr[8]+0.7905694150420947*fl[8]+0.6123724356957944*(fr[4]+fl[4]+fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaOrdR)+0.6846531968814573*fr[19]-0.6846531968814573*fl[19]+0.6846531968814574*fr[16]-0.6846531968814574*fl[16]-0.3952847075210473*(fr[15]+fl[15])+0.7905694150420948*(fr[12]+fl[12])-0.3952847075210473*(fr[9]+fl[9])+0.7905694150420947*(fr[8]+fl[8])-0.6123724356957944*(fr[4]+fr[2])+0.6123724356957944*(fl[4]+fl[2])+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaOrdR = 0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0]; 
  fUpOrd[5] = 0.5*(((-0.7348469228349525*(fr[19]+fl[19]))+1.42302494707577*fr[18]-1.42302494707577*fl[18]+0.7348469228349533*(fr[17]+fl[17])+0.5477225575051661*(fr[16]+fl[16])+0.4242640687119281*fr[15]-0.4242640687119281*fl[15]-1.060660171779821*fr[14]+1.060660171779821*fl[14]-0.4242640687119285*fr[13]+0.4242640687119285*fl[13]+1.060660171779821*fr[12]-1.060660171779821*fl[12]+0.5477225575051661*(fr[11]+fl[11])-1.10227038425243*(fr[10]+fl[10])-0.3162277660168379*fr[9]+0.3162277660168379*fl[9]-0.7905694150420947*fr[8]+0.7905694150420947*fl[8]-0.3162277660168379*fr[7]+0.3162277660168379*fl[7]+0.8215838362577489*(fr[6]+fl[6])+0.6363961030678926*fr[5]-0.6363961030678926*fl[5]-0.8215838362577489*(fr[4]+fl[4])-0.4743416490252568*fr[3]+0.4743416490252568*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.4743416490252568*fr[1]-0.4743416490252568*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaOrdR)+0.7348469228349525*fr[19]-0.7348469228349525*fl[19]-1.42302494707577*(fr[18]+fl[18])-0.7348469228349533*fr[17]+0.7348469228349533*fl[17]-0.5477225575051661*fr[16]+0.5477225575051661*fl[16]-0.4242640687119281*(fr[15]+fl[15])+1.060660171779821*(fr[14]+fl[14])+0.4242640687119285*(fr[13]+fl[13])-1.060660171779821*(fr[12]+fl[12])-0.5477225575051661*fr[11]+0.5477225575051661*fl[11]+1.10227038425243*fr[10]-1.10227038425243*fl[10]+0.3162277660168379*(fr[9]+fl[9])+0.7905694150420947*(fr[8]+fl[8])+0.3162277660168379*(fr[7]+fl[7])-0.8215838362577489*fr[6]+0.8215838362577489*fl[6]-0.6363961030678926*(fr[5]+fl[5])+0.8215838362577489*fr[4]-0.8215838362577489*fl[4]+0.4743416490252568*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.4743416490252568*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 
  fUpOrd[6] = 0.5*(((-0.6846531968814573*(fr[17]+fl[17]))-0.7905694150420948*fr[14]+0.7905694150420948*fl[14]+0.3952847075210473*fr[13]-0.3952847075210473*fl[13]-0.6846531968814574*(fr[11]+fl[11])-0.7905694150420947*fr[8]+0.7905694150420947*fl[8]+0.3952847075210473*fr[7]-0.3952847075210473*fl[7]+0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaOrdR)+0.6846531968814573*fr[17]-0.6846531968814573*fl[17]+0.7905694150420948*(fr[14]+fl[14])-0.3952847075210473*(fr[13]+fl[13])+0.6846531968814574*fr[11]-0.6846531968814574*fl[11]+0.7905694150420947*(fr[8]+fl[8])-0.3952847075210473*(fr[7]+fl[7])-0.6123724356957944*fr[6]+0.6123724356957944*fl[6]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0]; 
  fUpOrd[7] = 0.5*((0.7348469228349525*(fr[19]+fl[19])-1.42302494707577*fr[18]+1.42302494707577*fl[18]+0.7348469228349533*(fr[17]+fl[17])+0.5477225575051661*(fr[16]+fl[16])-0.4242640687119281*fr[15]+0.4242640687119281*fl[15]-1.060660171779821*fr[14]+1.060660171779821*fl[14]-0.4242640687119285*fr[13]+0.4242640687119285*fl[13]-1.060660171779821*fr[12]+1.060660171779821*fl[12]+0.5477225575051661*(fr[11]+fl[11])+1.10227038425243*(fr[10]+fl[10])-0.3162277660168379*fr[9]+0.3162277660168379*fl[9]-0.7905694150420947*fr[8]+0.7905694150420947*fl[8]-0.3162277660168379*fr[7]+0.3162277660168379*fl[7]+0.8215838362577489*(fr[6]+fl[6])-0.6363961030678926*fr[5]+0.6363961030678926*fl[5]+0.8215838362577489*(fr[4]+fl[4])-0.4743416490252568*fr[3]+0.4743416490252568*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.4743416490252568*fr[1]+0.4743416490252568*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaOrdR)-0.7348469228349525*fr[19]+0.7348469228349525*fl[19]+1.42302494707577*(fr[18]+fl[18])-0.7348469228349533*fr[17]+0.7348469228349533*fl[17]-0.5477225575051661*fr[16]+0.5477225575051661*fl[16]+0.4242640687119281*(fr[15]+fl[15])+1.060660171779821*(fr[14]+fl[14])+0.4242640687119285*(fr[13]+fl[13])+1.060660171779821*(fr[12]+fl[12])-0.5477225575051661*fr[11]+0.5477225575051661*fl[11]-1.10227038425243*fr[10]+1.10227038425243*fl[10]+0.3162277660168379*(fr[9]+fl[9])+0.7905694150420947*(fr[8]+fl[8])+0.3162277660168379*(fr[7]+fl[7])-0.8215838362577489*fr[6]+0.8215838362577489*fl[6]+0.6363961030678926*(fr[5]+fl[5])-0.8215838362577489*fr[4]+0.8215838362577489*fl[4]+0.4743416490252568*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.4743416490252568*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 

  double fUp[8];
  fUp[0] = 0.03846153846153846*(5.0*fUpOrd[7]+8.0*fUpOrd[6]+5.0*fUpOrd[5]+8.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]+8.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[1] = 0.00828173324999922*(25.0*(fUpOrd[7]-1.0*fUpOrd[5])+53.66563145999496*(fUpOrd[4]-1.0*fUpOrd[3])+25.0*(fUpOrd[2]-1.0*fUpOrd[0])); 
  fUp[2] = 0.01851851851851852*(11.18033988749895*fUpOrd[7]+24.0*fUpOrd[6]+11.18033988749895*fUpOrd[5]-1.0*(11.18033988749895*fUpOrd[2]+24.0*fUpOrd[1]+11.18033988749895*fUpOrd[0])); 
  fUp[3] = 0.2777777777777778*(fUpOrd[7]-1.0*(fUpOrd[5]+fUpOrd[2])+fUpOrd[0]); 
  fUp[4] = 0.03440104580768907*(5.0*fUpOrd[7]-18.0*fUpOrd[6]+5.0*fUpOrd[5]+8.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]-18.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[5] = 0.03440104580768907*(5.0*fUpOrd[7]+8.0*fUpOrd[6]+5.0*fUpOrd[5]-18.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]+8.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[6] = 0.08281733249999224*(2.23606797749979*fUpOrd[7]-6.0*fUpOrd[6]+2.23606797749979*(fUpOrd[5]-1.0*fUpOrd[2])+6.0*fUpOrd[1]-2.23606797749979*fUpOrd[0]); 
  fUp[7] = 0.01656346649999844*(11.18033988749895*(fUpOrd[7]-1.0*fUpOrd[5])+10.0*(3.0*fUpOrd[3]-3.0*fUpOrd[4])+11.18033988749895*(fUpOrd[2]-1.0*fUpOrd[0])); 

  incr[0] = 0.3535533905932737*(alpha[4]*fUp[4]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[1] = 0.07071067811865474*(4.47213595499958*(alpha[1]*fUp[4]+fUp[1]*alpha[4])+5.0*(alpha[0]*fUp[1]+fUp[0]*alpha[1])); 
  incr[2] = -0.6123724356957944*(alpha[4]*fUp[4]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[3] = 0.02357022603955158*(15.0*alpha[4]*fUp[6]+15.0*(alpha[1]*fUp[3]+alpha[0]*fUp[2])); 
  incr[4] = -0.07071067811865474*(7.745966692414834*(alpha[1]*fUp[4]+fUp[1]*alpha[4])+8.660254037844386*(alpha[0]*fUp[1]+fUp[0]*alpha[1])); 
  incr[5] = 0.02357022603955158*(13.41640786499874*alpha[1]*fUp[6]+13.41640786499874*fUp[3]*alpha[4]+15.0*(alpha[0]*fUp[3]+alpha[1]*fUp[2])); 
  incr[6] = -0.07071067811865474*(8.660254037844387*alpha[4]*fUp[6]+8.660254037844386*(alpha[1]*fUp[3]+alpha[0]*fUp[2])); 
  incr[7] = 0.01010152544552211*((22.3606797749979*alpha[4]+35.0*alpha[0])*fUp[4]+35.0*fUp[0]*alpha[4]+31.30495168499706*alpha[1]*fUp[1]); 
  incr[8] = 0.7905694150420947*(alpha[4]*fUp[4]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[9] = 0.02357022603955158*(15.0*alpha[1]*fUp[7]+15.0*alpha[0]*fUp[5]); 
  incr[10] = -0.07071067811865474*(7.745966692414834*(alpha[1]*fUp[6]+fUp[3]*alpha[4])+8.660254037844386*(alpha[0]*fUp[3]+alpha[1]*fUp[2])); 
  incr[11] = -0.01010152544552211*((38.72983346207417*alpha[4]+60.62177826491071*alpha[0])*fUp[4]+60.62177826491071*fUp[0]*alpha[4]+54.22176684690384*alpha[1]*fUp[1]); 
  incr[12] = 0.02357022603955158*(30.0*(alpha[1]*fUp[4]+fUp[1]*alpha[4])+33.54101966249684*(alpha[0]*fUp[1]+fUp[0]*alpha[1])); 
  incr[13] = 0.003367175148507369*((67.0820393249937*alpha[4]+105.0*alpha[0])*fUp[6]+105.0*fUp[2]*alpha[4]+93.91485505499116*alpha[1]*fUp[3]); 
  incr[14] = 0.1178511301977579*(6.708203932499369*alpha[4]*fUp[6]+6.708203932499369*(alpha[1]*fUp[3]+alpha[0]*fUp[2])); 
  incr[15] = 0.02357022603955158*((13.41640786499874*alpha[4]+15.0*alpha[0])*fUp[7]+15.0*alpha[1]*fUp[5]); 
  incr[16] = -0.07071067811865474*(8.660254037844386*alpha[1]*fUp[7]+8.660254037844387*alpha[0]*fUp[5]); 
  incr[17] = -0.01010152544552211*((38.72983346207417*alpha[4]+60.62177826491071*alpha[0])*fUp[6]+60.6217782649107*fUp[2]*alpha[4]+54.22176684690384*alpha[1]*fUp[3]); 
  incr[18] = 0.02357022603955158*(30.0*alpha[1]*fUp[6]+30.0*fUp[3]*alpha[4]+33.54101966249685*(alpha[0]*fUp[3]+alpha[1]*fUp[2])); 
  incr[19] = -0.07071067811865474*((7.745966692414834*alpha[4]+8.660254037844387*alpha[0])*fUp[7]+8.660254037844386*alpha[1]*fUp[5]); 

  outr[0] += incr[0]*dv10r; 
  outr[1] += incr[1]*dv10r; 
  outr[2] += incr[2]*dv10r; 
  outr[3] += incr[3]*dv10r; 
  outr[4] += incr[4]*dv10r; 
  outr[5] += incr[5]*dv10r; 
  outr[6] += incr[6]*dv10r; 
  outr[7] += incr[7]*dv10r; 
  outr[8] += incr[8]*dv10r; 
  outr[9] += incr[9]*dv10r; 
  outr[10] += incr[10]*dv10r; 
  outr[11] += incr[11]*dv10r; 
  outr[12] += incr[12]*dv10r; 
  outr[13] += incr[13]*dv10r; 
  outr[14] += incr[14]*dv10r; 
  outr[15] += incr[15]*dv10r; 
  outr[16] += incr[16]*dv10r; 
  outr[17] += incr[17]*dv10r; 
  outr[18] += incr[18]*dv10r; 
  outr[19] += incr[19]*dv10r; 

  outl[0] += -1.0*incr[0]*dv10l; 
  outl[1] += -1.0*incr[1]*dv10l; 
  outl[2] += incr[2]*dv10l; 
  outl[3] += -1.0*incr[3]*dv10l; 
  outl[4] += incr[4]*dv10l; 
  outl[5] += -1.0*incr[5]*dv10l; 
  outl[6] += incr[6]*dv10l; 
  outl[7] += -1.0*incr[7]*dv10l; 
  outl[8] += -1.0*incr[8]*dv10l; 
  outl[9] += -1.0*incr[9]*dv10l; 
  outl[10] += incr[10]*dv10l; 
  outl[11] += incr[11]*dv10l; 
  outl[12] += -1.0*incr[12]*dv10l; 
  outl[13] += -1.0*incr[13]*dv10l; 
  outl[14] += -1.0*incr[14]*dv10l; 
  outl[15] += -1.0*incr[15]*dv10l; 
  outl[16] += incr[16]*dv10l; 
  outl[17] += incr[17]*dv10l; 
  outl[18] += -1.0*incr[18]*dv10l; 
  outl[19] += incr[19]*dv10l; 


  return std::abs(amid); 
} 
__host__ __device__ double VlasovUpwindSurfNeutral1x2vSer_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[2]; 
  double dv11r = 2/dxvr[2]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double *Fo1 = &EM[3]; 

  double alpha[8]; 
  double incr[20]; 

  alpha[0] = 1.414213562373095*Fo1[0]; 
  alpha[1] = 1.414213562373095*Fo1[1]; 
  alpha[4] = 1.414213562373095*Fo1[2]; 

  const double amid = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 

  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0]; 
  fUpOrd[0] = 0.5*(((-1.42302494707577*fr[19])+1.42302494707577*fl[19]-0.7348469228349525*(fr[18]+fl[18])-0.7348469228349533*(fr[17]+fl[17])+1.060660171779821*(fr[16]+fr[15])-1.060660171779821*(fl[16]+fl[15])+0.5477225575051661*(fr[14]+fl[14]+fr[13]+fl[13])+0.4242640687119281*fr[12]-0.4242640687119281*fl[12]+0.4242640687119285*fr[11]-0.4242640687119285*fl[11]+1.10227038425243*(fr[10]+fl[10])-0.7905694150420947*fr[9]+0.7905694150420947*fl[9]-0.3162277660168379*(fr[8]+fr[7])+0.3162277660168379*(fl[8]+fl[7])-0.8215838362577489*(fr[6]+fl[6]+fr[5]+fl[5])-0.6363961030678926*fr[4]+0.6363961030678926*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.4743416490252568*(fr[2]+fr[1])-0.4743416490252568*(fl[2]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaOrdR)+1.42302494707577*(fr[19]+fl[19])+0.7348469228349525*fr[18]-0.7348469228349525*fl[18]+0.7348469228349533*fr[17]-0.7348469228349533*fl[17]-1.060660171779821*(fr[16]+fl[16]+fr[15]+fl[15])-0.5477225575051661*(fr[14]+fr[13])+0.5477225575051661*(fl[14]+fl[13])-0.4242640687119281*(fr[12]+fl[12])-0.4242640687119285*(fr[11]+fl[11])-1.10227038425243*fr[10]+1.10227038425243*fl[10]+0.7905694150420947*(fr[9]+fl[9])+0.3162277660168379*(fr[8]+fl[8]+fr[7]+fl[7])+0.8215838362577489*(fr[6]+fr[5])-0.8215838362577489*(fl[6]+fl[5])+0.6363961030678926*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.4743416490252568*(fr[2]+fl[2]+fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 
  fUpOrd[1] = 0.5*((0.6846531968814573*(fr[17]+fl[17])+0.7905694150420948*fr[16]-0.7905694150420948*fl[16]-0.6846531968814574*(fr[13]+fl[13])-0.3952847075210473*fr[11]+0.3952847075210473*fl[11]-0.7905694150420947*fr[9]+0.7905694150420947*fl[9]+0.3952847075210473*fr[7]-0.3952847075210473*fl[7]-0.6123724356957944*(fr[6]+fl[6])+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*fr[2]-0.3535533905932737*(fl[2]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaOrdR)-0.6846531968814573*fr[17]+0.6846531968814573*fl[17]-0.7905694150420948*(fr[16]+fl[16])+0.6846531968814574*fr[13]-0.6846531968814574*fl[13]+0.3952847075210473*(fr[11]+fl[11])+0.7905694150420947*(fr[9]+fl[9])-0.3952847075210473*(fr[7]+fl[7])+0.6123724356957944*fr[6]-0.6123724356957944*(fl[6]+fr[3])+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0]; 
  fUpOrd[2] = 0.5*((1.42302494707577*fr[19]-1.42302494707577*fl[19]+0.7348469228349525*(fr[18]+fl[18])-0.7348469228349533*(fr[17]+fl[17])+1.060660171779821*fr[16]-1.060660171779821*(fl[16]+fr[15])+1.060660171779821*fl[15]+0.5477225575051661*(fr[14]+fl[14]+fr[13]+fl[13])-0.4242640687119281*fr[12]+0.4242640687119281*fl[12]+0.4242640687119285*fr[11]-0.4242640687119285*fl[11]-1.10227038425243*(fr[10]+fl[10])-0.7905694150420947*fr[9]+0.7905694150420947*fl[9]-0.3162277660168379*(fr[8]+fr[7])+0.3162277660168379*(fl[8]+fl[7])-0.8215838362577489*(fr[6]+fl[6])+0.8215838362577489*(fr[5]+fl[5])+0.6363961030678926*fr[4]-0.6363961030678926*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.4743416490252568*fr[2]-0.4743416490252568*(fl[2]+fr[1])+0.4743416490252568*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaOrdR)-1.42302494707577*(fr[19]+fl[19])-0.7348469228349525*fr[18]+0.7348469228349525*fl[18]+0.7348469228349533*fr[17]-0.7348469228349533*fl[17]-1.060660171779821*(fr[16]+fl[16])+1.060660171779821*(fr[15]+fl[15])-0.5477225575051661*(fr[14]+fr[13])+0.5477225575051661*(fl[14]+fl[13])+0.4242640687119281*(fr[12]+fl[12])-0.4242640687119285*(fr[11]+fl[11])+1.10227038425243*fr[10]-1.10227038425243*fl[10]+0.7905694150420947*(fr[9]+fl[9])+0.3162277660168379*(fr[8]+fl[8]+fr[7]+fl[7])+0.8215838362577489*fr[6]-0.8215838362577489*(fl[6]+fr[5])+0.8215838362577489*fl[5]-0.6363961030678926*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.4743416490252568*(fr[2]+fl[2])+0.4743416490252568*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*alpha[0]-0.5*alpha[1]; 
  fUpOrd[3] = 0.5*((0.6846531968814573*(fr[18]+fl[18])+0.7905694150420948*fr[15]-0.7905694150420948*fl[15]-0.6846531968814574*(fr[14]+fl[14])-0.3952847075210473*fr[12]+0.3952847075210473*fl[12]-0.7905694150420947*fr[9]+0.7905694150420947*fl[9]+0.3952847075210473*fr[8]-0.3952847075210473*fl[8]-0.6123724356957944*(fr[5]+fl[5])+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaOrdR)-0.6846531968814573*fr[18]+0.6846531968814573*fl[18]-0.7905694150420948*(fr[15]+fl[15])+0.6846531968814574*fr[14]-0.6846531968814574*fl[14]+0.3952847075210473*(fr[12]+fl[12])+0.7905694150420947*(fr[9]+fl[9])-0.3952847075210473*(fr[8]+fl[8])+0.6123724356957944*fr[5]-0.6123724356957944*(fl[5]+fr[3])+0.6123724356957944*fl[3]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*(alpha[1]+alpha[0]); 
  fUpOrd[4] = 0.5*(((-0.6846531968814573*(fr[18]+fl[18]))-0.7905694150420948*fr[15]+0.7905694150420948*fl[15]-0.6846531968814574*(fr[14]+fl[14])+0.3952847075210473*fr[12]-0.3952847075210473*fl[12]-0.7905694150420947*fr[9]+0.7905694150420947*fl[9]+0.3952847075210473*fr[8]-0.3952847075210473*fl[8]+0.6123724356957944*(fr[5]+fl[5]+fr[3]+fl[3])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaOrdR)+0.6846531968814573*fr[18]-0.6846531968814573*fl[18]+0.7905694150420948*(fr[15]+fl[15])+0.6846531968814574*fr[14]-0.6846531968814574*fl[14]-0.3952847075210473*(fr[12]+fl[12])+0.7905694150420947*(fr[9]+fl[9])-0.3952847075210473*(fr[8]+fl[8])-0.6123724356957944*(fr[5]+fr[3])+0.6123724356957944*(fl[5]+fl[3])+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaOrdR = 0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0]; 
  fUpOrd[5] = 0.5*((1.42302494707577*fr[19]-1.42302494707577*fl[19]-0.7348469228349525*(fr[18]+fl[18])+0.7348469228349533*(fr[17]+fl[17])-1.060660171779821*fr[16]+1.060660171779821*(fl[16]+fr[15])-1.060660171779821*fl[15]+0.5477225575051661*(fr[14]+fl[14]+fr[13]+fl[13])+0.4242640687119281*fr[12]-0.4242640687119281*fl[12]-0.4242640687119285*fr[11]+0.4242640687119285*fl[11]-1.10227038425243*(fr[10]+fl[10])-0.7905694150420947*fr[9]+0.7905694150420947*fl[9]-0.3162277660168379*(fr[8]+fr[7])+0.3162277660168379*(fl[8]+fl[7])+0.8215838362577489*(fr[6]+fl[6])-0.8215838362577489*(fr[5]+fl[5])+0.6363961030678926*fr[4]-0.6363961030678926*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.4743416490252568*fr[2]+0.4743416490252568*(fl[2]+fr[1])-0.4743416490252568*fl[1]-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaOrdR)-1.42302494707577*(fr[19]+fl[19])+0.7348469228349525*fr[18]-0.7348469228349525*fl[18]-0.7348469228349533*fr[17]+0.7348469228349533*fl[17]+1.060660171779821*(fr[16]+fl[16])-1.060660171779821*(fr[15]+fl[15])-0.5477225575051661*(fr[14]+fr[13])+0.5477225575051661*(fl[14]+fl[13])-0.4242640687119281*(fr[12]+fl[12])+0.4242640687119285*(fr[11]+fl[11])+1.10227038425243*fr[10]-1.10227038425243*fl[10]+0.7905694150420947*(fr[9]+fl[9])+0.3162277660168379*(fr[8]+fl[8]+fr[7]+fl[7])-0.8215838362577489*fr[6]+0.8215838362577489*(fl[6]+fr[5])-0.8215838362577489*fl[5]-0.6363961030678926*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.4743416490252568*(fr[2]+fl[2])-0.4743416490252568*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*alpha[0]-0.5590169943749475*alpha[4]; 
  fUpOrd[6] = 0.5*(((-0.6846531968814573*(fr[17]+fl[17]))-0.7905694150420948*fr[16]+0.7905694150420948*fl[16]-0.6846531968814574*(fr[13]+fl[13])+0.3952847075210473*fr[11]-0.3952847075210473*fl[11]-0.7905694150420947*fr[9]+0.7905694150420947*fl[9]+0.3952847075210473*fr[7]-0.3952847075210473*fl[7]+0.6123724356957944*(fr[6]+fl[6]+fr[3]+fl[3])-0.3535533905932737*(fr[2]+fr[0])+0.3535533905932737*(fl[2]+fl[0]))*sgn(alphaOrdR)+0.6846531968814573*fr[17]-0.6846531968814573*fl[17]+0.7905694150420948*(fr[16]+fl[16])+0.6846531968814574*fr[13]-0.6846531968814574*fl[13]-0.3952847075210473*(fr[11]+fl[11])+0.7905694150420947*(fr[9]+fl[9])-0.3952847075210473*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fr[3])+0.6123724356957944*(fl[6]+fl[3])+0.3535533905932737*(fr[2]+fl[2]+fr[0]+fl[0])); 
  alphaOrdR = 0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0]; 
  fUpOrd[7] = 0.5*(((-1.42302494707577*fr[19])+1.42302494707577*fl[19]+0.7348469228349525*(fr[18]+fl[18])+0.7348469228349533*(fr[17]+fl[17])-1.060660171779821*(fr[16]+fr[15])+1.060660171779821*(fl[16]+fl[15])+0.5477225575051661*(fr[14]+fl[14]+fr[13]+fl[13])-0.4242640687119281*fr[12]+0.4242640687119281*fl[12]-0.4242640687119285*fr[11]+0.4242640687119285*fl[11]+1.10227038425243*(fr[10]+fl[10])-0.7905694150420947*fr[9]+0.7905694150420947*fl[9]-0.3162277660168379*(fr[8]+fr[7])+0.3162277660168379*(fl[8]+fl[7])+0.8215838362577489*(fr[6]+fl[6]+fr[5]+fl[5])-0.6363961030678926*fr[4]+0.6363961030678926*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.4743416490252568*(fr[2]+fr[1])+0.4743416490252568*(fl[2]+fl[1])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaOrdR)+1.42302494707577*(fr[19]+fl[19])-0.7348469228349525*fr[18]+0.7348469228349525*fl[18]-0.7348469228349533*fr[17]+0.7348469228349533*fl[17]+1.060660171779821*(fr[16]+fl[16]+fr[15]+fl[15])-0.5477225575051661*(fr[14]+fr[13])+0.5477225575051661*(fl[14]+fl[13])+0.4242640687119281*(fr[12]+fl[12])+0.4242640687119285*(fr[11]+fl[11])-1.10227038425243*fr[10]+1.10227038425243*fl[10]+0.7905694150420947*(fr[9]+fl[9])+0.3162277660168379*(fr[8]+fl[8]+fr[7]+fl[7])-0.8215838362577489*(fr[6]+fr[5])+0.8215838362577489*(fl[6]+fl[5])+0.6363961030678926*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.4743416490252568*(fr[2]+fl[2]+fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 

  double fUp[8];
  fUp[0] = 0.03846153846153846*(5.0*fUpOrd[7]+8.0*fUpOrd[6]+5.0*fUpOrd[5]+8.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]+8.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[1] = 0.00828173324999922*(25.0*(fUpOrd[7]-1.0*fUpOrd[5])+53.66563145999496*(fUpOrd[4]-1.0*fUpOrd[3])+25.0*(fUpOrd[2]-1.0*fUpOrd[0])); 
  fUp[2] = 0.01851851851851852*(11.18033988749895*fUpOrd[7]+24.0*fUpOrd[6]+11.18033988749895*fUpOrd[5]-1.0*(11.18033988749895*fUpOrd[2]+24.0*fUpOrd[1]+11.18033988749895*fUpOrd[0])); 
  fUp[3] = 0.2777777777777778*(fUpOrd[7]-1.0*(fUpOrd[5]+fUpOrd[2])+fUpOrd[0]); 
  fUp[4] = 0.03440104580768907*(5.0*fUpOrd[7]-18.0*fUpOrd[6]+5.0*fUpOrd[5]+8.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]-18.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[5] = 0.03440104580768907*(5.0*fUpOrd[7]+8.0*fUpOrd[6]+5.0*fUpOrd[5]-18.0*(fUpOrd[4]+fUpOrd[3])+5.0*fUpOrd[2]+8.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[6] = 0.08281733249999224*(2.23606797749979*fUpOrd[7]-6.0*fUpOrd[6]+2.23606797749979*(fUpOrd[5]-1.0*fUpOrd[2])+6.0*fUpOrd[1]-2.23606797749979*fUpOrd[0]); 
  fUp[7] = 0.01656346649999844*(11.18033988749895*(fUpOrd[7]-1.0*fUpOrd[5])+10.0*(3.0*fUpOrd[3]-3.0*fUpOrd[4])+11.18033988749895*(fUpOrd[2]-1.0*fUpOrd[0])); 

  incr[0] = 0.3535533905932737*(alpha[4]*fUp[4]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[1] = 0.07071067811865474*(4.47213595499958*(alpha[1]*fUp[4]+fUp[1]*alpha[4])+5.0*(alpha[0]*fUp[1]+fUp[0]*alpha[1])); 
  incr[2] = 0.02357022603955158*(15.0*alpha[4]*fUp[6]+15.0*(alpha[1]*fUp[3]+alpha[0]*fUp[2])); 
  incr[3] = -0.6123724356957944*(alpha[4]*fUp[4]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[4] = 0.02357022603955158*(13.41640786499874*alpha[1]*fUp[6]+13.41640786499874*fUp[3]*alpha[4]+15.0*(alpha[0]*fUp[3]+alpha[1]*fUp[2])); 
  incr[5] = -0.07071067811865474*(7.745966692414834*(alpha[1]*fUp[4]+fUp[1]*alpha[4])+8.660254037844386*(alpha[0]*fUp[1]+fUp[0]*alpha[1])); 
  incr[6] = -0.07071067811865474*(8.660254037844387*alpha[4]*fUp[6]+8.660254037844386*(alpha[1]*fUp[3]+alpha[0]*fUp[2])); 
  incr[7] = 0.01010152544552211*((22.3606797749979*alpha[4]+35.0*alpha[0])*fUp[4]+35.0*fUp[0]*alpha[4]+31.30495168499706*alpha[1]*fUp[1]); 
  incr[8] = 0.02357022603955158*(15.0*alpha[1]*fUp[7]+15.0*alpha[0]*fUp[5]); 
  incr[9] = 0.7905694150420947*(alpha[4]*fUp[4]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[10] = -0.07071067811865474*(7.745966692414834*(alpha[1]*fUp[6]+fUp[3]*alpha[4])+8.660254037844386*(alpha[0]*fUp[3]+alpha[1]*fUp[2])); 
  incr[11] = 0.003367175148507369*((67.0820393249937*alpha[4]+105.0*alpha[0])*fUp[6]+105.0*fUp[2]*alpha[4]+93.91485505499116*alpha[1]*fUp[3]); 
  incr[12] = 0.02357022603955158*((13.41640786499874*alpha[4]+15.0*alpha[0])*fUp[7]+15.0*alpha[1]*fUp[5]); 
  incr[13] = -0.01010152544552211*((38.72983346207417*alpha[4]+60.62177826491071*alpha[0])*fUp[4]+60.62177826491071*fUp[0]*alpha[4]+54.22176684690384*alpha[1]*fUp[1]); 
  incr[14] = -0.07071067811865474*(8.660254037844386*alpha[1]*fUp[7]+8.660254037844387*alpha[0]*fUp[5]); 
  incr[15] = 0.02357022603955158*(30.0*(alpha[1]*fUp[4]+fUp[1]*alpha[4])+33.54101966249684*(alpha[0]*fUp[1]+fUp[0]*alpha[1])); 
  incr[16] = 0.1178511301977579*(6.708203932499369*alpha[4]*fUp[6]+6.708203932499369*(alpha[1]*fUp[3]+alpha[0]*fUp[2])); 
  incr[17] = -0.01010152544552211*((38.72983346207417*alpha[4]+60.62177826491071*alpha[0])*fUp[6]+60.6217782649107*fUp[2]*alpha[4]+54.22176684690384*alpha[1]*fUp[3]); 
  incr[18] = -0.07071067811865474*((7.745966692414834*alpha[4]+8.660254037844387*alpha[0])*fUp[7]+8.660254037844386*alpha[1]*fUp[5]); 
  incr[19] = 0.02357022603955158*(30.0*alpha[1]*fUp[6]+30.0*fUp[3]*alpha[4]+33.54101966249685*(alpha[0]*fUp[3]+alpha[1]*fUp[2])); 

  outr[0] += incr[0]*dv11r; 
  outr[1] += incr[1]*dv11r; 
  outr[2] += incr[2]*dv11r; 
  outr[3] += incr[3]*dv11r; 
  outr[4] += incr[4]*dv11r; 
  outr[5] += incr[5]*dv11r; 
  outr[6] += incr[6]*dv11r; 
  outr[7] += incr[7]*dv11r; 
  outr[8] += incr[8]*dv11r; 
  outr[9] += incr[9]*dv11r; 
  outr[10] += incr[10]*dv11r; 
  outr[11] += incr[11]*dv11r; 
  outr[12] += incr[12]*dv11r; 
  outr[13] += incr[13]*dv11r; 
  outr[14] += incr[14]*dv11r; 
  outr[15] += incr[15]*dv11r; 
  outr[16] += incr[16]*dv11r; 
  outr[17] += incr[17]*dv11r; 
  outr[18] += incr[18]*dv11r; 
  outr[19] += incr[19]*dv11r; 

  outl[0] += -1.0*incr[0]*dv11l; 
  outl[1] += -1.0*incr[1]*dv11l; 
  outl[2] += -1.0*incr[2]*dv11l; 
  outl[3] += incr[3]*dv11l; 
  outl[4] += -1.0*incr[4]*dv11l; 
  outl[5] += incr[5]*dv11l; 
  outl[6] += incr[6]*dv11l; 
  outl[7] += -1.0*incr[7]*dv11l; 
  outl[8] += -1.0*incr[8]*dv11l; 
  outl[9] += -1.0*incr[9]*dv11l; 
  outl[10] += incr[10]*dv11l; 
  outl[11] += -1.0*incr[11]*dv11l; 
  outl[12] += -1.0*incr[12]*dv11l; 
  outl[13] += incr[13]*dv11l; 
  outl[14] += incr[14]*dv11l; 
  outl[15] += -1.0*incr[15]*dv11l; 
  outl[16] += -1.0*incr[16]*dv11l; 
  outl[17] += incr[17]*dv11l; 
  outl[18] += incr[18]*dv11l; 
  outl[19] += -1.0*incr[19]*dv11l; 


  return std::abs(amid); 
} 
