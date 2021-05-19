#include <VmLBOModDecl.h> 
double VmLBOconstNuSurfNonUniform1x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[10]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 

  double fjump[10]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 

  double Ghat[10]; 
  for(unsigned short int i=0; i<10; ++i){ 
    Ghat[i]=0.0; 
  }; 

  const double dxvl1R2 = pow(dxvl[1],2);
  const double dxvl1R3 = pow(dxvl[1],3);
  const double dxvl1R4 = pow(dxvl[1],4);
  const double dxvl1R5 = pow(dxvl[1],5);
  const double dxvl1R6 = pow(dxvl[1],6);
  const double dxvr1R2 = pow(dxvr[1],2);
  const double dxvr1R3 = pow(dxvr[1],3);
  const double dxvr1R4 = pow(dxvr[1],4);
  const double dxvr1R5 = pow(dxvr[1],5);
  const double dxvr1R6 = pow(dxvr[1],6);

  Ghat[0] = ((375.6594202199647*nuVtSqSum[0]*dxvl1R3*dxvr1R3+214.6625258399798*nuVtSqSum[0]*dxvl1R4*dxvr1R2-93.91485505499116*nuVtSqSum[0]*dxvl1R5*dxvr[1]-67.0820393249937*nuVtSqSum[0]*dxvl1R6)*fr[8]+(67.0820393249937*nuVtSqSum[0]*dxvr1R6+93.91485505499116*nuVtSqSum[0]*dxvl[1]*dxvr1R5-214.6625258399798*nuVtSqSum[0]*dxvl1R2*dxvr1R4-375.6594202199647*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[8]+600.0*dxvl1R3*dxvr1R3*nuVtSqSum[2]*fr[7]-600.0*dxvl1R3*dxvr1R3*nuVtSqSum[2]*fl[7]+((-727.4613391789284*dxvl1R3*dxvr1R3)-138.5640646055102*dxvl1R4*dxvr1R2+86.60254037844386*dxvl1R5*dxvr[1]+17.32050807568877*dxvl1R6)*nuVtSqSum[1]*fr[4]+(17.32050807568877*dxvr1R6+86.60254037844386*dxvl[1]*dxvr1R5-138.5640646055102*dxvl1R2*dxvr1R4-727.4613391789284*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[4]+((-727.4613391789284*nuVtSqSum[0]*dxvl1R3*dxvr1R3)-138.5640646055102*nuVtSqSum[0]*dxvl1R4*dxvr1R2+86.60254037844386*nuVtSqSum[0]*dxvl1R5*dxvr[1]+17.32050807568877*nuVtSqSum[0]*dxvl1R6)*fr[2]+(17.32050807568877*nuVtSqSum[0]*dxvr1R6+86.60254037844386*nuVtSqSum[0]*dxvl[1]*dxvr1R5-138.5640646055102*nuVtSqSum[0]*dxvl1R2*dxvr1R4-727.4613391789284*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[2]+(600.0*dxvl1R3*dxvr1R3*fr[1]-600.0*dxvl1R3*dxvr1R3*fl[1])*nuVtSqSum[1]+(600.0*fr[0]-600.0*fl[0])*nuVtSqSum[0]*dxvl1R3*dxvr1R3)/(7.071067811865476*dxvl[1]*dxvr1R6+35.35533905932738*dxvl1R2*dxvr1R5+70.71067811865477*dxvl1R3*dxvr1R4+70.71067811865477*dxvl1R4*dxvr1R3+35.35533905932738*dxvl1R5*dxvr1R2+7.071067811865476*dxvl1R6*dxvr[1])-0.5*(2.23606797749979*fjump[8]+1.732050807568877*fjump[2]+fjump[0])+0.3535533905932737*(2.23606797749979*alphaDrag[0]*favg[8]+alphaDrag[2]*favg[7]+1.732050807568877*alphaDrag[1]*favg[4]+1.732050807568877*alphaDrag[0]*favg[2]+alphaDrag[1]*favg[1]+alphaDrag[0]*favg[0]); 
  Ghat[1] = ((375.6594202199647*dxvl1R3*dxvr1R3+214.6625258399798*dxvl1R4*dxvr1R2-93.91485505499116*dxvl1R5*dxvr[1]-67.0820393249937*dxvl1R6)*nuVtSqSum[1]*fr[8]+(67.0820393249937*dxvr1R6+93.91485505499116*dxvl[1]*dxvr1R5-214.6625258399798*dxvl1R2*dxvr1R4-375.6594202199647*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[8]+536.6563145999496*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fr[7]-536.6563145999496*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fl[7]+(((-650.6612021628459*dxvl1R3*dxvr1R3)-123.9354670786373*dxvl1R4*dxvr1R2+77.45966692414835*dxvl1R5*dxvr[1]+15.49193338482967*dxvl1R6)*nuVtSqSum[2]-727.4613391789284*nuVtSqSum[0]*dxvl1R3*dxvr1R3-138.5640646055102*nuVtSqSum[0]*dxvl1R4*dxvr1R2+86.60254037844386*nuVtSqSum[0]*dxvl1R5*dxvr[1]+17.32050807568877*nuVtSqSum[0]*dxvl1R6)*fr[4]+((15.49193338482967*dxvr1R6+77.45966692414835*dxvl[1]*dxvr1R5-123.9354670786373*dxvl1R2*dxvr1R4-650.6612021628459*dxvl1R3*dxvr1R3)*nuVtSqSum[2]+17.32050807568877*nuVtSqSum[0]*dxvr1R6+86.60254037844386*nuVtSqSum[0]*dxvl[1]*dxvr1R5-138.5640646055102*nuVtSqSum[0]*dxvl1R2*dxvr1R4-727.4613391789284*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[4]+(536.6563145999496*dxvl1R3*dxvr1R3*fr[1]-536.6563145999496*dxvl1R3*dxvr1R3*fl[1])*nuVtSqSum[2]+((-727.4613391789284*dxvl1R3*dxvr1R3)-138.5640646055102*dxvl1R4*dxvr1R2+86.60254037844386*dxvl1R5*dxvr[1]+17.32050807568877*dxvl1R6)*nuVtSqSum[1]*fr[2]+(17.32050807568877*dxvr1R6+86.60254037844386*dxvl[1]*dxvr1R5-138.5640646055102*dxvl1R2*dxvr1R4-727.4613391789284*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[2]+(600.0*fr[0]-600.0*fl[0])*dxvl1R3*dxvr1R3*nuVtSqSum[1]+600.0*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fr[1]-600.0*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fl[1])/(7.071067811865476*dxvl[1]*dxvr1R6+35.35533905932738*dxvl1R2*dxvr1R5+70.71067811865477*dxvl1R3*dxvr1R4+70.71067811865477*dxvl1R4*dxvr1R3+35.35533905932738*dxvl1R5*dxvr1R2+7.071067811865476*dxvl1R6*dxvr[1])+0.07071067811865474*(11.18033988749895*alphaDrag[1]*favg[8]+4.47213595499958*alphaDrag[1]*favg[7]+(7.745966692414834*alphaDrag[2]+8.660254037844386*alphaDrag[0])*favg[4]+8.660254037844386*alphaDrag[1]*favg[2]+4.47213595499958*favg[1]*alphaDrag[2]+5.0*alphaDrag[0]*favg[1]+5.0*favg[0]*alphaDrag[1])-0.5*(1.732050807568877*fjump[4]+fjump[1]); 
  Ghat[3] = (-(1.0*((145.4922678357857*nuVtSqSum[0]*dxvl1R3*dxvr1R3+27.71281292110204*nuVtSqSum[0]*dxvl1R4*dxvr1R2-17.32050807568877*nuVtSqSum[0]*dxvl1R5*dxvr[1]-3.464101615137754*nuVtSqSum[0]*dxvl1R6)*fr[6]+((-3.464101615137754*nuVtSqSum[0]*dxvr1R6)-17.32050807568877*nuVtSqSum[0]*dxvl[1]*dxvr1R5+27.71281292110204*nuVtSqSum[0]*dxvl1R2*dxvr1R4+145.4922678357857*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[6]-120.0*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fr[5]+120.0*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fl[5]-120.0*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fr[3]+120.0*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fl[3]))/(1.414213562373095*dxvl[1]*dxvr1R6+7.071067811865476*dxvl1R2*dxvr1R5+14.14213562373095*dxvl1R3*dxvr1R4+14.14213562373095*dxvl1R4*dxvr1R3+7.071067811865476*dxvl1R5*dxvr1R2+1.414213562373095*dxvl1R6*dxvr[1]))-0.5*(1.732050807568877*fjump[6]+fjump[3])+0.3535533905932737*(1.732050807568877*alphaDrag[0]*favg[6]+alphaDrag[1]*favg[5]+alphaDrag[0]*favg[3]); 
  Ghat[5] = (-(1.0*((145.4922678357857*dxvl1R3*dxvr1R3+27.71281292110204*dxvl1R4*dxvr1R2-17.32050807568877*dxvl1R5*dxvr[1]-3.464101615137754*dxvl1R6)*nuVtSqSum[1]*fr[6]+((-3.464101615137754*dxvr1R6)-17.32050807568877*dxvl[1]*dxvr1R5+27.71281292110204*dxvl1R2*dxvr1R4+145.4922678357857*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[6]+((-107.3312629199899*dxvl1R3*dxvr1R3*nuVtSqSum[2])-120.0*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fr[5]+(107.3312629199899*dxvl1R3*dxvr1R3*nuVtSqSum[2]+120.0*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[5]-120.0*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fr[3]+120.0*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fl[3]))/(1.414213562373095*dxvl[1]*dxvr1R6+7.071067811865476*dxvl1R2*dxvr1R5+14.14213562373095*dxvl1R3*dxvr1R4+14.14213562373095*dxvl1R4*dxvr1R3+7.071067811865476*dxvl1R5*dxvr1R2+1.414213562373095*dxvl1R6*dxvr[1]))+0.07071067811865474*(8.660254037844386*alphaDrag[1]*favg[6]+(4.47213595499958*alphaDrag[2]+5.0*alphaDrag[0])*favg[5]+5.0*alphaDrag[1]*favg[3])-0.5*fjump[5]; 
  Ghat[7] = ((2629.615941539753*dxvl1R3*dxvr1R3+1502.637680879859*dxvl1R4*dxvr1R2-657.4039853849382*dxvl1R5*dxvr[1]-469.5742752749559*dxvl1R6)*nuVtSqSum[2]*fr[8]+(469.5742752749559*dxvr1R6+657.4039853849382*dxvl[1]*dxvr1R5-1502.637680879859*dxvl1R2*dxvr1R4-2629.615941539753*dxvl1R3*dxvr1R3)*nuVtSqSum[2]*fl[8]+(2683.281572999748*dxvl1R3*dxvr1R3*nuVtSqSum[2]+4200.0*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fr[7]+((-2683.281572999748*dxvl1R3*dxvr1R3*nuVtSqSum[2])-4200.0*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[7]+((-4554.628415139921*dxvl1R3*dxvr1R3)-867.5482695504614*dxvl1R4*dxvr1R2+542.2176684690385*dxvl1R5*dxvr[1]+108.4435336938077*dxvl1R6)*nuVtSqSum[1]*fr[4]+(108.4435336938077*dxvr1R6+542.2176684690385*dxvl[1]*dxvr1R5-867.5482695504614*dxvl1R2*dxvr1R4-4554.628415139921*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[4]+(((-5092.229374252499*dxvl1R3*dxvr1R3)-969.9484522385712*dxvl1R4*dxvr1R2+606.217782649107*dxvl1R5*dxvr[1]+121.2435565298214*dxvl1R6)*fr[2]+(121.2435565298214*dxvr1R6+606.217782649107*dxvl[1]*dxvr1R5-969.9484522385712*dxvl1R2*dxvr1R4-5092.229374252499*dxvl1R3*dxvr1R3)*fl[2]+(4200.0*fr[0]-4200.0*fl[0])*dxvl1R3*dxvr1R3)*nuVtSqSum[2]+(3756.594202199647*dxvl1R3*dxvr1R3*fr[1]-3756.594202199647*dxvl1R3*dxvr1R3*fl[1])*nuVtSqSum[1])/(49.49747468305833*dxvl[1]*dxvr1R6+247.4873734152916*dxvl1R2*dxvr1R5+494.9747468305834*dxvl1R3*dxvr1R4+494.9747468305834*dxvl1R4*dxvr1R3+247.4873734152916*dxvl1R5*dxvr1R2+49.49747468305833*dxvl1R6*dxvr[1])+0.01010152544552211*(78.26237921249266*alphaDrag[2]*favg[8]+(22.3606797749979*alphaDrag[2]+35.0*alphaDrag[0])*favg[7]+54.22176684690384*alphaDrag[1]*favg[4]+60.6217782649107*alphaDrag[2]*favg[2]+35.0*favg[0]*alphaDrag[2]+31.30495168499706*alphaDrag[1]*favg[1])-0.5*fjump[7]; 
  Ghat[9] = (120.0*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[9]-120.0*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[9])/(1.414213562373095*dxvr1R5+7.071067811865476*dxvl[1]*dxvr1R4+14.14213562373095*dxvl1R2*dxvr1R3+14.14213562373095*dxvl1R3*dxvr1R2+7.071067811865476*dxvl1R4*dxvr[1]+1.414213562373095*dxvl1R5)-0.5*fjump[9]+0.3535533905932737*alphaDrag[0]*favg[9]; 

  double incr1[10]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[1]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[3]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = -1.118033988749895*Ghat[0]; 
  incr1[9] = -0.5*Ghat[9]; 

  double incr2[10]; 

  incr2[2] = ((54.22176684690384*nuVtSqSum[0]*dxvl1R3*dxvr1R2+61.96773353931867*nuVtSqSum[0]*dxvl1R4*dxvr[1]+19.36491673103709*nuVtSqSum[0]*dxvl1R5)*fr[8]+(19.36491673103709*nuVtSqSum[0]*dxvr1R5+61.96773353931867*nuVtSqSum[0]*dxvl[1]*dxvr1R4+54.22176684690384*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[8]+(86.60254037844386*dxvl1R3*dxvr1R2+43.30127018922193*dxvl1R4*dxvr[1]+8.660254037844386*dxvl1R5)*nuVtSqSum[2]*fr[7]+(8.660254037844386*dxvr1R5+43.30127018922193*dxvl[1]*dxvr1R4+86.60254037844386*dxvl1R2*dxvr1R3)*nuVtSqSum[2]*fl[7]+((-105.0*dxvl1R3*dxvr1R2)-75.0*dxvl1R4*dxvr[1]-15.0*dxvl1R5)*nuVtSqSum[1]*fr[4]+(15.0*dxvr1R5+75.0*dxvl[1]*dxvr1R4+105.0*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[4]+((-105.0*nuVtSqSum[0]*dxvl1R3*dxvr1R2)-75.0*nuVtSqSum[0]*dxvl1R4*dxvr[1]-15.0*nuVtSqSum[0]*dxvl1R5)*fr[2]+(15.0*nuVtSqSum[0]*dxvr1R5+75.0*nuVtSqSum[0]*dxvl[1]*dxvr1R4+105.0*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[2]+((86.60254037844386*dxvl1R3*dxvr1R2+43.30127018922193*dxvl1R4*dxvr[1]+8.660254037844386*dxvl1R5)*fr[1]+(8.660254037844386*dxvr1R5+43.30127018922193*dxvl[1]*dxvr1R4+86.60254037844386*dxvl1R2*dxvr1R3)*fl[1])*nuVtSqSum[1]+8.660254037844386*fl[0]*nuVtSqSum[0]*dxvr1R5+43.30127018922193*fl[0]*nuVtSqSum[0]*dxvl[1]*dxvr1R4+86.60254037844386*fl[0]*nuVtSqSum[0]*dxvl1R2*dxvr1R3+86.60254037844386*fr[0]*nuVtSqSum[0]*dxvl1R3*dxvr1R2+43.30127018922193*fr[0]*nuVtSqSum[0]*dxvl1R4*dxvr[1]+8.660254037844386*fr[0]*nuVtSqSum[0]*dxvl1R5)/(14.14213562373095*dxvr1R5+70.71067811865477*dxvl[1]*dxvr1R4+141.4213562373096*dxvl1R2*dxvr1R3+141.4213562373096*dxvl1R3*dxvr1R2+70.71067811865477*dxvl1R4*dxvr[1]+14.14213562373095*dxvl1R5); 
  incr2[4] = ((54.22176684690384*dxvl1R3*dxvr1R2+61.96773353931867*dxvl1R4*dxvr[1]+19.36491673103709*dxvl1R5)*nuVtSqSum[1]*fr[8]+(19.36491673103709*dxvr1R5+61.96773353931867*dxvl[1]*dxvr1R4+54.22176684690384*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[8]+(77.45966692414835*dxvl1R3*dxvr1R2+38.72983346207418*dxvl1R4*dxvr[1]+7.745966692414834*dxvl1R5)*nuVtSqSum[1]*fr[7]+(7.745966692414834*dxvr1R5+38.72983346207418*dxvl[1]*dxvr1R4+77.45966692414835*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[7]+(((-93.91485505499116*dxvl1R3*dxvr1R2)-67.0820393249937*dxvl1R4*dxvr[1]-13.41640786499874*dxvl1R5)*nuVtSqSum[2]-105.0*nuVtSqSum[0]*dxvl1R3*dxvr1R2-75.0*nuVtSqSum[0]*dxvl1R4*dxvr[1]-15.0*nuVtSqSum[0]*dxvl1R5)*fr[4]+((13.41640786499874*dxvr1R5+67.0820393249937*dxvl[1]*dxvr1R4+93.91485505499116*dxvl1R2*dxvr1R3)*nuVtSqSum[2]+15.0*nuVtSqSum[0]*dxvr1R5+75.0*nuVtSqSum[0]*dxvl[1]*dxvr1R4+105.0*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[4]+((77.45966692414835*dxvl1R3*dxvr1R2+38.72983346207418*dxvl1R4*dxvr[1]+7.745966692414834*dxvl1R5)*fr[1]+(7.745966692414834*dxvr1R5+38.72983346207418*dxvl[1]*dxvr1R4+77.45966692414835*dxvl1R2*dxvr1R3)*fl[1])*nuVtSqSum[2]+((-105.0*dxvl1R3*dxvr1R2)-75.0*dxvl1R4*dxvr[1]-15.0*dxvl1R5)*nuVtSqSum[1]*fr[2]+(15.0*dxvr1R5+75.0*dxvl[1]*dxvr1R4+105.0*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[2]+(8.660254037844386*fl[0]*dxvr1R5+43.30127018922193*fl[0]*dxvl[1]*dxvr1R4+86.60254037844386*fl[0]*dxvl1R2*dxvr1R3+86.60254037844386*fr[0]*dxvl1R3*dxvr1R2+43.30127018922193*fr[0]*dxvl1R4*dxvr[1]+8.660254037844386*fr[0]*dxvl1R5)*nuVtSqSum[1]+(86.60254037844386*nuVtSqSum[0]*dxvl1R3*dxvr1R2+43.30127018922193*nuVtSqSum[0]*dxvl1R4*dxvr[1]+8.660254037844386*nuVtSqSum[0]*dxvl1R5)*fr[1]+(8.660254037844386*nuVtSqSum[0]*dxvr1R5+43.30127018922193*nuVtSqSum[0]*dxvl[1]*dxvr1R4+86.60254037844386*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[1])/(14.14213562373095*dxvr1R5+70.71067811865477*dxvl[1]*dxvr1R4+141.4213562373096*dxvl1R2*dxvr1R3+141.4213562373096*dxvl1R3*dxvr1R2+70.71067811865477*dxvl1R4*dxvr[1]+14.14213562373095*dxvl1R5); 
  incr2[6] = -(1.0*((21.0*nuVtSqSum[0]*dxvl1R3*dxvr1R2+15.0*nuVtSqSum[0]*dxvl1R4*dxvr[1]+3.0*nuVtSqSum[0]*dxvl1R5)*fr[6]+((-3.0*nuVtSqSum[0]*dxvr1R5)-15.0*nuVtSqSum[0]*dxvl[1]*dxvr1R4-21.0*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[6]+((-17.32050807568877*dxvl1R3*dxvr1R2)-8.660254037844386*dxvl1R4*dxvr[1]-1.732050807568877*dxvl1R5)*nuVtSqSum[1]*fr[5]+((-1.732050807568877*dxvr1R5)-8.660254037844386*dxvl[1]*dxvr1R4-17.32050807568877*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[5]+((-17.32050807568877*nuVtSqSum[0]*dxvl1R3*dxvr1R2)-8.660254037844386*nuVtSqSum[0]*dxvl1R4*dxvr[1]-1.732050807568877*nuVtSqSum[0]*dxvl1R5)*fr[3]+((-1.732050807568877*nuVtSqSum[0]*dxvr1R5)-8.660254037844386*nuVtSqSum[0]*dxvl[1]*dxvr1R4-17.32050807568877*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[3]))/(2.828427124746191*dxvr1R5+14.14213562373095*dxvl[1]*dxvr1R4+28.28427124746191*dxvl1R2*dxvr1R3+28.28427124746191*dxvl1R3*dxvr1R2+14.14213562373095*dxvl1R4*dxvr[1]+2.828427124746191*dxvl1R5); 
  incr2[8] = -(1.0*((42.0*nuVtSqSum[0]*dxvl1R3*dxvr1R2+48.0*nuVtSqSum[0]*dxvl1R4*dxvr[1]+15.0*nuVtSqSum[0]*dxvl1R5)*fr[8]+(15.0*nuVtSqSum[0]*dxvr1R5+48.0*nuVtSqSum[0]*dxvl[1]*dxvr1R4+42.0*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[8]+(67.0820393249937*dxvl1R3*dxvr1R2+33.54101966249685*dxvl1R4*dxvr[1]+6.708203932499369*dxvl1R5)*nuVtSqSum[2]*fr[7]+(6.708203932499369*dxvr1R5+33.54101966249685*dxvl[1]*dxvr1R4+67.0820393249937*dxvl1R2*dxvr1R3)*nuVtSqSum[2]*fl[7]+((-81.33265027035574*dxvl1R3*dxvr1R2)-58.09475019311126*dxvl1R4*dxvr[1]-11.61895003862225*dxvl1R5)*nuVtSqSum[1]*fr[4]+(11.61895003862225*dxvr1R5+58.09475019311126*dxvl[1]*dxvr1R4+81.33265027035574*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[4]+((-81.33265027035574*nuVtSqSum[0]*dxvl1R3*dxvr1R2)-58.09475019311126*nuVtSqSum[0]*dxvl1R4*dxvr[1]-11.61895003862225*nuVtSqSum[0]*dxvl1R5)*fr[2]+(11.61895003862225*nuVtSqSum[0]*dxvr1R5+58.09475019311126*nuVtSqSum[0]*dxvl[1]*dxvr1R4+81.33265027035574*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[2]+((67.0820393249937*dxvl1R3*dxvr1R2+33.54101966249685*dxvl1R4*dxvr[1]+6.708203932499369*dxvl1R5)*fr[1]+(6.708203932499369*dxvr1R5+33.54101966249685*dxvl[1]*dxvr1R4+67.0820393249937*dxvl1R2*dxvr1R3)*fl[1])*nuVtSqSum[1]+6.708203932499369*fl[0]*nuVtSqSum[0]*dxvr1R5+33.54101966249685*fl[0]*nuVtSqSum[0]*dxvl[1]*dxvr1R4+67.0820393249937*fl[0]*nuVtSqSum[0]*dxvl1R2*dxvr1R3+67.0820393249937*fr[0]*nuVtSqSum[0]*dxvl1R3*dxvr1R2+33.54101966249685*fr[0]*nuVtSqSum[0]*dxvl1R4*dxvr[1]+6.708203932499369*fr[0]*nuVtSqSum[0]*dxvl1R5))/(2.828427124746191*dxvr1R5+14.14213562373095*dxvl[1]*dxvr1R4+28.28427124746191*dxvl1R2*dxvr1R3+28.28427124746191*dxvl1R3*dxvr1R2+14.14213562373095*dxvl1R4*dxvr[1]+2.828427124746191*dxvl1R5); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr2[8]*rdvSq4R+incr1[8]*rdv2R; 
  outr[9] += incr1[9]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += incr2[8]*rdvSq4L-1.0*incr1[8]*rdv2L; 
  outl[9] += -1.0*incr1[9]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUx[2])/nuSum-(0.7071067811865475*sumNuUx[0])/nuSum+wl[1]); 
} 
double VmLBOconstNuSurfNonUniform1x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUy = &nuUSum[3]; 

  double favg[10]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 

  double fjump[10]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[2]*nuSum+0.7071067811865475*dxvl[2]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 

  double Ghat[10]; 
  for(unsigned short int i=0; i<10; ++i){ 
    Ghat[i]=0.0; 
  }; 

  const double dxvl2R2 = pow(dxvl[2],2);
  const double dxvl2R3 = pow(dxvl[2],3);
  const double dxvl2R4 = pow(dxvl[2],4);
  const double dxvl2R5 = pow(dxvl[2],5);
  const double dxvl2R6 = pow(dxvl[2],6);
  const double dxvr2R2 = pow(dxvr[2],2);
  const double dxvr2R3 = pow(dxvr[2],3);
  const double dxvr2R4 = pow(dxvr[2],4);
  const double dxvr2R5 = pow(dxvr[2],5);
  const double dxvr2R6 = pow(dxvr[2],6);

  Ghat[0] = ((375.6594202199647*nuVtSqSum[0]*dxvl2R3*dxvr2R3+214.6625258399798*nuVtSqSum[0]*dxvl2R4*dxvr2R2-93.91485505499116*nuVtSqSum[0]*dxvl2R5*dxvr[2]-67.0820393249937*nuVtSqSum[0]*dxvl2R6)*fr[9]+(67.0820393249937*nuVtSqSum[0]*dxvr2R6+93.91485505499116*nuVtSqSum[0]*dxvl[2]*dxvr2R5-214.6625258399798*nuVtSqSum[0]*dxvl2R2*dxvr2R4-375.6594202199647*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[9]+600.0*dxvl2R3*dxvr2R3*nuVtSqSum[2]*fr[7]-600.0*dxvl2R3*dxvr2R3*nuVtSqSum[2]*fl[7]+((-727.4613391789284*nuVtSqSum[1]*dxvl2R3*dxvr2R3)-138.5640646055102*nuVtSqSum[1]*dxvl2R4*dxvr2R2+86.60254037844386*nuVtSqSum[1]*dxvl2R5*dxvr[2]+17.32050807568877*nuVtSqSum[1]*dxvl2R6)*fr[5]+(17.32050807568877*nuVtSqSum[1]*dxvr2R6+86.60254037844386*nuVtSqSum[1]*dxvl[2]*dxvr2R5-138.5640646055102*nuVtSqSum[1]*dxvl2R2*dxvr2R4-727.4613391789284*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[5]+((-727.4613391789284*nuVtSqSum[0]*dxvl2R3*dxvr2R3)-138.5640646055102*nuVtSqSum[0]*dxvl2R4*dxvr2R2+86.60254037844386*nuVtSqSum[0]*dxvl2R5*dxvr[2]+17.32050807568877*nuVtSqSum[0]*dxvl2R6)*fr[3]+(17.32050807568877*nuVtSqSum[0]*dxvr2R6+86.60254037844386*nuVtSqSum[0]*dxvl[2]*dxvr2R5-138.5640646055102*nuVtSqSum[0]*dxvl2R2*dxvr2R4-727.4613391789284*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[3]+((600.0*fr[1]-600.0*fl[1])*nuVtSqSum[1]+(600.0*fr[0]-600.0*fl[0])*nuVtSqSum[0])*dxvl2R3*dxvr2R3)/(7.071067811865476*dxvl[2]*dxvr2R6+35.35533905932738*dxvl2R2*dxvr2R5+70.71067811865477*dxvl2R3*dxvr2R4+70.71067811865477*dxvl2R4*dxvr2R3+35.35533905932738*dxvl2R5*dxvr2R2+7.071067811865476*dxvl2R6*dxvr[2])-0.5*(2.23606797749979*fjump[9]+1.732050807568877*fjump[3]+fjump[0])+0.3535533905932737*(2.23606797749979*alphaDrag[0]*favg[9]+alphaDrag[2]*favg[7]+1.732050807568877*alphaDrag[1]*favg[5]+1.732050807568877*alphaDrag[0]*favg[3]+alphaDrag[1]*favg[1]+alphaDrag[0]*favg[0]); 
  Ghat[1] = ((375.6594202199647*nuVtSqSum[1]*dxvl2R3*dxvr2R3+214.6625258399798*nuVtSqSum[1]*dxvl2R4*dxvr2R2-93.91485505499116*nuVtSqSum[1]*dxvl2R5*dxvr[2]-67.0820393249937*nuVtSqSum[1]*dxvl2R6)*fr[9]+(67.0820393249937*nuVtSqSum[1]*dxvr2R6+93.91485505499116*nuVtSqSum[1]*dxvl[2]*dxvr2R5-214.6625258399798*nuVtSqSum[1]*dxvl2R2*dxvr2R4-375.6594202199647*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[9]+536.6563145999496*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fr[7]-536.6563145999496*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fl[7]+(((-650.6612021628459*dxvl2R3*dxvr2R3)-123.9354670786373*dxvl2R4*dxvr2R2+77.45966692414835*dxvl2R5*dxvr[2]+15.49193338482967*dxvl2R6)*nuVtSqSum[2]-727.4613391789284*nuVtSqSum[0]*dxvl2R3*dxvr2R3-138.5640646055102*nuVtSqSum[0]*dxvl2R4*dxvr2R2+86.60254037844386*nuVtSqSum[0]*dxvl2R5*dxvr[2]+17.32050807568877*nuVtSqSum[0]*dxvl2R6)*fr[5]+((15.49193338482967*dxvr2R6+77.45966692414835*dxvl[2]*dxvr2R5-123.9354670786373*dxvl2R2*dxvr2R4-650.6612021628459*dxvl2R3*dxvr2R3)*nuVtSqSum[2]+17.32050807568877*nuVtSqSum[0]*dxvr2R6+86.60254037844386*nuVtSqSum[0]*dxvl[2]*dxvr2R5-138.5640646055102*nuVtSqSum[0]*dxvl2R2*dxvr2R4-727.4613391789284*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[5]+((-727.4613391789284*nuVtSqSum[1]*dxvl2R3*dxvr2R3)-138.5640646055102*nuVtSqSum[1]*dxvl2R4*dxvr2R2+86.60254037844386*nuVtSqSum[1]*dxvl2R5*dxvr[2]+17.32050807568877*nuVtSqSum[1]*dxvl2R6)*fr[3]+(17.32050807568877*nuVtSqSum[1]*dxvr2R6+86.60254037844386*nuVtSqSum[1]*dxvl[2]*dxvr2R5-138.5640646055102*nuVtSqSum[1]*dxvl2R2*dxvr2R4-727.4613391789284*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[3]+(536.6563145999496*fr[1]-536.6563145999496*fl[1])*dxvl2R3*dxvr2R3*nuVtSqSum[2]+((600.0*fr[0]-600.0*fl[0])*nuVtSqSum[1]+600.0*nuVtSqSum[0]*fr[1]-600.0*nuVtSqSum[0]*fl[1])*dxvl2R3*dxvr2R3)/(7.071067811865476*dxvl[2]*dxvr2R6+35.35533905932738*dxvl2R2*dxvr2R5+70.71067811865477*dxvl2R3*dxvr2R4+70.71067811865477*dxvl2R4*dxvr2R3+35.35533905932738*dxvl2R5*dxvr2R2+7.071067811865476*dxvl2R6*dxvr[2])+0.07071067811865474*(11.18033988749895*alphaDrag[1]*favg[9]+4.47213595499958*alphaDrag[1]*favg[7]+(7.745966692414834*alphaDrag[2]+8.660254037844386*alphaDrag[0])*favg[5]+8.660254037844386*alphaDrag[1]*favg[3]+4.47213595499958*favg[1]*alphaDrag[2]+5.0*alphaDrag[0]*favg[1]+5.0*favg[0]*alphaDrag[1])-0.5*(1.732050807568877*fjump[5]+fjump[1]); 
  Ghat[2] = (-(1.0*((145.4922678357857*nuVtSqSum[0]*dxvl2R3*dxvr2R3+27.71281292110204*nuVtSqSum[0]*dxvl2R4*dxvr2R2-17.32050807568877*nuVtSqSum[0]*dxvl2R5*dxvr[2]-3.464101615137754*nuVtSqSum[0]*dxvl2R6)*fr[6]+((-3.464101615137754*nuVtSqSum[0]*dxvr2R6)-17.32050807568877*nuVtSqSum[0]*dxvl[2]*dxvr2R5+27.71281292110204*nuVtSqSum[0]*dxvl2R2*dxvr2R4+145.4922678357857*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[6]-120.0*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fr[4]+120.0*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fl[4]-120.0*nuVtSqSum[0]*dxvl2R3*dxvr2R3*fr[2]+120.0*nuVtSqSum[0]*dxvl2R3*dxvr2R3*fl[2]))/(1.414213562373095*dxvl[2]*dxvr2R6+7.071067811865476*dxvl2R2*dxvr2R5+14.14213562373095*dxvl2R3*dxvr2R4+14.14213562373095*dxvl2R4*dxvr2R3+7.071067811865476*dxvl2R5*dxvr2R2+1.414213562373095*dxvl2R6*dxvr[2]))-0.5*(1.732050807568877*fjump[6]+fjump[2])+0.3535533905932737*(1.732050807568877*alphaDrag[0]*favg[6]+alphaDrag[1]*favg[4]+alphaDrag[0]*favg[2]); 
  Ghat[4] = (-(1.0*((145.4922678357857*nuVtSqSum[1]*dxvl2R3*dxvr2R3+27.71281292110204*nuVtSqSum[1]*dxvl2R4*dxvr2R2-17.32050807568877*nuVtSqSum[1]*dxvl2R5*dxvr[2]-3.464101615137754*nuVtSqSum[1]*dxvl2R6)*fr[6]+((-3.464101615137754*nuVtSqSum[1]*dxvr2R6)-17.32050807568877*nuVtSqSum[1]*dxvl[2]*dxvr2R5+27.71281292110204*nuVtSqSum[1]*dxvl2R2*dxvr2R4+145.4922678357857*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[6]+((-107.3312629199899*dxvl2R3*dxvr2R3*nuVtSqSum[2])-120.0*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fr[4]+(107.3312629199899*dxvl2R3*dxvr2R3*nuVtSqSum[2]+120.0*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[4]-120.0*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fr[2]+120.0*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fl[2]))/(1.414213562373095*dxvl[2]*dxvr2R6+7.071067811865476*dxvl2R2*dxvr2R5+14.14213562373095*dxvl2R3*dxvr2R4+14.14213562373095*dxvl2R4*dxvr2R3+7.071067811865476*dxvl2R5*dxvr2R2+1.414213562373095*dxvl2R6*dxvr[2]))+0.07071067811865474*(8.660254037844386*alphaDrag[1]*favg[6]+(4.47213595499958*alphaDrag[2]+5.0*alphaDrag[0])*favg[4]+5.0*alphaDrag[1]*favg[2])-0.5*fjump[4]; 
  Ghat[7] = ((2629.615941539753*dxvl2R3*dxvr2R3+1502.637680879859*dxvl2R4*dxvr2R2-657.4039853849382*dxvl2R5*dxvr[2]-469.5742752749559*dxvl2R6)*nuVtSqSum[2]*fr[9]+(469.5742752749559*dxvr2R6+657.4039853849382*dxvl[2]*dxvr2R5-1502.637680879859*dxvl2R2*dxvr2R4-2629.615941539753*dxvl2R3*dxvr2R3)*nuVtSqSum[2]*fl[9]+(2683.281572999748*dxvl2R3*dxvr2R3*nuVtSqSum[2]+4200.0*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fr[7]+((-2683.281572999748*dxvl2R3*dxvr2R3*nuVtSqSum[2])-4200.0*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[7]+((-4554.628415139921*nuVtSqSum[1]*dxvl2R3*dxvr2R3)-867.5482695504614*nuVtSqSum[1]*dxvl2R4*dxvr2R2+542.2176684690385*nuVtSqSum[1]*dxvl2R5*dxvr[2]+108.4435336938077*nuVtSqSum[1]*dxvl2R6)*fr[5]+(108.4435336938077*nuVtSqSum[1]*dxvr2R6+542.2176684690385*nuVtSqSum[1]*dxvl[2]*dxvr2R5-867.5482695504614*nuVtSqSum[1]*dxvl2R2*dxvr2R4-4554.628415139921*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[5]+((-5092.229374252499*dxvl2R3*dxvr2R3)-969.9484522385712*dxvl2R4*dxvr2R2+606.217782649107*dxvl2R5*dxvr[2]+121.2435565298214*dxvl2R6)*nuVtSqSum[2]*fr[3]+(121.2435565298214*dxvr2R6+606.217782649107*dxvl[2]*dxvr2R5-969.9484522385712*dxvl2R2*dxvr2R4-5092.229374252499*dxvl2R3*dxvr2R3)*nuVtSqSum[2]*fl[3]+(4200.0*fr[0]-4200.0*fl[0])*dxvl2R3*dxvr2R3*nuVtSqSum[2]+(3756.594202199647*fr[1]-3756.594202199647*fl[1])*nuVtSqSum[1]*dxvl2R3*dxvr2R3)/(49.49747468305833*dxvl[2]*dxvr2R6+247.4873734152916*dxvl2R2*dxvr2R5+494.9747468305834*dxvl2R3*dxvr2R4+494.9747468305834*dxvl2R4*dxvr2R3+247.4873734152916*dxvl2R5*dxvr2R2+49.49747468305833*dxvl2R6*dxvr[2])+0.01010152544552211*(78.26237921249266*alphaDrag[2]*favg[9]+(22.3606797749979*alphaDrag[2]+35.0*alphaDrag[0])*favg[7]+54.22176684690384*alphaDrag[1]*favg[5]+60.6217782649107*alphaDrag[2]*favg[3]+35.0*favg[0]*alphaDrag[2]+31.30495168499706*alphaDrag[1]*favg[1])-0.5*fjump[7]; 
  Ghat[8] = (120.0*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[8]-120.0*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[8])/(1.414213562373095*dxvr2R5+7.071067811865476*dxvl[2]*dxvr2R4+14.14213562373095*dxvl2R2*dxvr2R3+14.14213562373095*dxvl2R3*dxvr2R2+7.071067811865476*dxvl2R4*dxvr[2]+1.414213562373095*dxvl2R5)-0.5*fjump[8]+0.3535533905932737*alphaDrag[0]*favg[8]; 

  double incr1[10]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = 0.8660254037844386*Ghat[2]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = -1.118033988749895*Ghat[0]; 

  double incr2[10]; 

  incr2[3] = ((54.22176684690384*nuVtSqSum[0]*dxvl2R3*dxvr2R2+61.96773353931867*nuVtSqSum[0]*dxvl2R4*dxvr[2]+19.36491673103709*nuVtSqSum[0]*dxvl2R5)*fr[9]+(19.36491673103709*nuVtSqSum[0]*dxvr2R5+61.96773353931867*nuVtSqSum[0]*dxvl[2]*dxvr2R4+54.22176684690384*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[9]+(86.60254037844386*dxvl2R3*dxvr2R2+43.30127018922193*dxvl2R4*dxvr[2]+8.660254037844386*dxvl2R5)*nuVtSqSum[2]*fr[7]+(8.660254037844386*dxvr2R5+43.30127018922193*dxvl[2]*dxvr2R4+86.60254037844386*dxvl2R2*dxvr2R3)*nuVtSqSum[2]*fl[7]+((-105.0*nuVtSqSum[1]*dxvl2R3*dxvr2R2)-75.0*nuVtSqSum[1]*dxvl2R4*dxvr[2]-15.0*nuVtSqSum[1]*dxvl2R5)*fr[5]+(15.0*nuVtSqSum[1]*dxvr2R5+75.0*nuVtSqSum[1]*dxvl[2]*dxvr2R4+105.0*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[5]+((-105.0*nuVtSqSum[0]*dxvl2R3*dxvr2R2)-75.0*nuVtSqSum[0]*dxvl2R4*dxvr[2]-15.0*nuVtSqSum[0]*dxvl2R5)*fr[3]+(15.0*nuVtSqSum[0]*dxvr2R5+75.0*nuVtSqSum[0]*dxvl[2]*dxvr2R4+105.0*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[3]+(8.660254037844386*fl[1]*nuVtSqSum[1]+8.660254037844386*fl[0]*nuVtSqSum[0])*dxvr2R5+(43.30127018922193*fl[1]*nuVtSqSum[1]+43.30127018922193*fl[0]*nuVtSqSum[0])*dxvl[2]*dxvr2R4+(86.60254037844386*fl[1]*nuVtSqSum[1]+86.60254037844386*fl[0]*nuVtSqSum[0])*dxvl2R2*dxvr2R3+(86.60254037844386*fr[1]*nuVtSqSum[1]+86.60254037844386*fr[0]*nuVtSqSum[0])*dxvl2R3*dxvr2R2+(43.30127018922193*fr[1]*nuVtSqSum[1]+43.30127018922193*fr[0]*nuVtSqSum[0])*dxvl2R4*dxvr[2]+(8.660254037844386*fr[1]*nuVtSqSum[1]+8.660254037844386*fr[0]*nuVtSqSum[0])*dxvl2R5)/(14.14213562373095*dxvr2R5+70.71067811865477*dxvl[2]*dxvr2R4+141.4213562373096*dxvl2R2*dxvr2R3+141.4213562373096*dxvl2R3*dxvr2R2+70.71067811865477*dxvl2R4*dxvr[2]+14.14213562373095*dxvl2R5); 
  incr2[5] = ((54.22176684690384*nuVtSqSum[1]*dxvl2R3*dxvr2R2+61.96773353931867*nuVtSqSum[1]*dxvl2R4*dxvr[2]+19.36491673103709*nuVtSqSum[1]*dxvl2R5)*fr[9]+(19.36491673103709*nuVtSqSum[1]*dxvr2R5+61.96773353931867*nuVtSqSum[1]*dxvl[2]*dxvr2R4+54.22176684690384*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[9]+(77.45966692414835*nuVtSqSum[1]*dxvl2R3*dxvr2R2+38.72983346207418*nuVtSqSum[1]*dxvl2R4*dxvr[2]+7.745966692414834*nuVtSqSum[1]*dxvl2R5)*fr[7]+(7.745966692414834*nuVtSqSum[1]*dxvr2R5+38.72983346207418*nuVtSqSum[1]*dxvl[2]*dxvr2R4+77.45966692414835*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[7]+(((-93.91485505499116*dxvl2R3*dxvr2R2)-67.0820393249937*dxvl2R4*dxvr[2]-13.41640786499874*dxvl2R5)*nuVtSqSum[2]-105.0*nuVtSqSum[0]*dxvl2R3*dxvr2R2-75.0*nuVtSqSum[0]*dxvl2R4*dxvr[2]-15.0*nuVtSqSum[0]*dxvl2R5)*fr[5]+((13.41640786499874*dxvr2R5+67.0820393249937*dxvl[2]*dxvr2R4+93.91485505499116*dxvl2R2*dxvr2R3)*nuVtSqSum[2]+15.0*nuVtSqSum[0]*dxvr2R5+75.0*nuVtSqSum[0]*dxvl[2]*dxvr2R4+105.0*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[5]+((-105.0*nuVtSqSum[1]*dxvl2R3*dxvr2R2)-75.0*nuVtSqSum[1]*dxvl2R4*dxvr[2]-15.0*nuVtSqSum[1]*dxvl2R5)*fr[3]+(15.0*nuVtSqSum[1]*dxvr2R5+75.0*nuVtSqSum[1]*dxvl[2]*dxvr2R4+105.0*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[3]+(7.745966692414834*fl[1]*dxvr2R5+38.72983346207418*fl[1]*dxvl[2]*dxvr2R4+77.45966692414835*fl[1]*dxvl2R2*dxvr2R3+77.45966692414835*fr[1]*dxvl2R3*dxvr2R2+38.72983346207418*fr[1]*dxvl2R4*dxvr[2]+7.745966692414834*fr[1]*dxvl2R5)*nuVtSqSum[2]+(8.660254037844386*fl[0]*nuVtSqSum[1]+8.660254037844386*nuVtSqSum[0]*fl[1])*dxvr2R5+(43.30127018922193*fl[0]*nuVtSqSum[1]+43.30127018922193*nuVtSqSum[0]*fl[1])*dxvl[2]*dxvr2R4+(86.60254037844386*fl[0]*nuVtSqSum[1]+86.60254037844386*nuVtSqSum[0]*fl[1])*dxvl2R2*dxvr2R3+(86.60254037844386*fr[0]*nuVtSqSum[1]+86.60254037844386*nuVtSqSum[0]*fr[1])*dxvl2R3*dxvr2R2+(43.30127018922193*fr[0]*nuVtSqSum[1]+43.30127018922193*nuVtSqSum[0]*fr[1])*dxvl2R4*dxvr[2]+(8.660254037844386*fr[0]*nuVtSqSum[1]+8.660254037844386*nuVtSqSum[0]*fr[1])*dxvl2R5)/(14.14213562373095*dxvr2R5+70.71067811865477*dxvl[2]*dxvr2R4+141.4213562373096*dxvl2R2*dxvr2R3+141.4213562373096*dxvl2R3*dxvr2R2+70.71067811865477*dxvl2R4*dxvr[2]+14.14213562373095*dxvl2R5); 
  incr2[6] = -(1.0*((21.0*nuVtSqSum[0]*dxvl2R3*dxvr2R2+15.0*nuVtSqSum[0]*dxvl2R4*dxvr[2]+3.0*nuVtSqSum[0]*dxvl2R5)*fr[6]+((-3.0*nuVtSqSum[0]*dxvr2R5)-15.0*nuVtSqSum[0]*dxvl[2]*dxvr2R4-21.0*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[6]+((-17.32050807568877*nuVtSqSum[1]*dxvl2R3*dxvr2R2)-8.660254037844386*nuVtSqSum[1]*dxvl2R4*dxvr[2]-1.732050807568877*nuVtSqSum[1]*dxvl2R5)*fr[4]+((-1.732050807568877*nuVtSqSum[1]*dxvr2R5)-8.660254037844386*nuVtSqSum[1]*dxvl[2]*dxvr2R4-17.32050807568877*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[4]+((-17.32050807568877*nuVtSqSum[0]*dxvl2R3*dxvr2R2)-8.660254037844386*nuVtSqSum[0]*dxvl2R4*dxvr[2]-1.732050807568877*nuVtSqSum[0]*dxvl2R5)*fr[2]+((-1.732050807568877*nuVtSqSum[0]*dxvr2R5)-8.660254037844386*nuVtSqSum[0]*dxvl[2]*dxvr2R4-17.32050807568877*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[2]))/(2.828427124746191*dxvr2R5+14.14213562373095*dxvl[2]*dxvr2R4+28.28427124746191*dxvl2R2*dxvr2R3+28.28427124746191*dxvl2R3*dxvr2R2+14.14213562373095*dxvl2R4*dxvr[2]+2.828427124746191*dxvl2R5); 
  incr2[9] = -(1.0*((42.0*nuVtSqSum[0]*dxvl2R3*dxvr2R2+48.0*nuVtSqSum[0]*dxvl2R4*dxvr[2]+15.0*nuVtSqSum[0]*dxvl2R5)*fr[9]+(15.0*nuVtSqSum[0]*dxvr2R5+48.0*nuVtSqSum[0]*dxvl[2]*dxvr2R4+42.0*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[9]+(67.0820393249937*dxvl2R3*dxvr2R2+33.54101966249685*dxvl2R4*dxvr[2]+6.708203932499369*dxvl2R5)*nuVtSqSum[2]*fr[7]+(6.708203932499369*dxvr2R5+33.54101966249685*dxvl[2]*dxvr2R4+67.0820393249937*dxvl2R2*dxvr2R3)*nuVtSqSum[2]*fl[7]+((-81.33265027035574*nuVtSqSum[1]*dxvl2R3*dxvr2R2)-58.09475019311126*nuVtSqSum[1]*dxvl2R4*dxvr[2]-11.61895003862225*nuVtSqSum[1]*dxvl2R5)*fr[5]+(11.61895003862225*nuVtSqSum[1]*dxvr2R5+58.09475019311126*nuVtSqSum[1]*dxvl[2]*dxvr2R4+81.33265027035574*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[5]+((-81.33265027035574*nuVtSqSum[0]*dxvl2R3*dxvr2R2)-58.09475019311126*nuVtSqSum[0]*dxvl2R4*dxvr[2]-11.61895003862225*nuVtSqSum[0]*dxvl2R5)*fr[3]+(11.61895003862225*nuVtSqSum[0]*dxvr2R5+58.09475019311126*nuVtSqSum[0]*dxvl[2]*dxvr2R4+81.33265027035574*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[3]+(6.708203932499369*fl[1]*nuVtSqSum[1]+6.708203932499369*fl[0]*nuVtSqSum[0])*dxvr2R5+(33.54101966249685*fl[1]*nuVtSqSum[1]+33.54101966249685*fl[0]*nuVtSqSum[0])*dxvl[2]*dxvr2R4+(67.0820393249937*fl[1]*nuVtSqSum[1]+67.0820393249937*fl[0]*nuVtSqSum[0])*dxvl2R2*dxvr2R3+(67.0820393249937*fr[1]*nuVtSqSum[1]+67.0820393249937*fr[0]*nuVtSqSum[0])*dxvl2R3*dxvr2R2+(33.54101966249685*fr[1]*nuVtSqSum[1]+33.54101966249685*fr[0]*nuVtSqSum[0])*dxvl2R4*dxvr[2]+(6.708203932499369*fr[1]*nuVtSqSum[1]+6.708203932499369*fr[0]*nuVtSqSum[0])*dxvl2R5))/(2.828427124746191*dxvr2R5+14.14213562373095*dxvl[2]*dxvr2R4+28.28427124746191*dxvl2R2*dxvr2R3+28.28427124746191*dxvl2R3*dxvr2R2+14.14213562373095*dxvl2R4*dxvr[2]+2.828427124746191*dxvl2R5); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += incr2[9]*rdvSq4L-1.0*incr1[9]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUy[2])/nuSum-(0.7071067811865475*sumNuUy[0])/nuSum+wl[2]); 
} 
