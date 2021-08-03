#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x2vTensor_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
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

  double favg[27]; 
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
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = -1*fr[16]+fl[16]; 
  favg[17] = -1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = -1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  favg[21] = 1*fr[21]+fl[21]; 
  favg[22] = 1*fr[22]+fl[22]; 
  favg[23] = 1*fr[23]+fl[23]; 
  favg[24] = -1*fr[24]+fl[24]; 
  favg[25] = 1*fr[25]+fl[25]; 
  favg[26] = 1*fr[26]+fl[26]; 

  double fjump[27]; 
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
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(-1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(-1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(-1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(-1*fr[19])); 
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(1*fr[20])); 
  fjump[21] = nuSum*vMuMidMax*(fl[21]-(1*fr[21])); 
  fjump[22] = nuSum*vMuMidMax*(fl[22]-(1*fr[22])); 
  fjump[23] = nuSum*vMuMidMax*(fl[23]-(1*fr[23])); 
  fjump[24] = nuSum*vMuMidMax*(fl[24]-(-1*fr[24])); 
  fjump[25] = nuSum*vMuMidMax*(fl[25]-(1*fr[25])); 
  fjump[26] = nuSum*vMuMidMax*(fl[26]-(1*fr[26])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 

  double Gdiff[27]; 
  double Ghat[27]; 
  double incr2[27]; 


  incr2[2] = 0.0110485434560398*(27.11088342345192*nuVtSqSum[2]*fr[20]+27.11088342345192*nuVtSqSum[2]*fl[20]+27.11088342345192*nuVtSqSum[1]*fr[12]+27.11088342345192*nuVtSqSum[1]*fl[12]-39.0*nuVtSqSum[2]*fr[11]+39.0*nuVtSqSum[2]*fl[11]+27.11088342345192*nuVtSqSum[0]*fr[8]+27.11088342345192*nuVtSqSum[0]*fl[8]+27.71281292110204*nuVtSqSum[2]*fr[7]+27.71281292110204*nuVtSqSum[2]*fl[7]-39.0*nuVtSqSum[1]*fr[4]+39.0*nuVtSqSum[1]*fl[4]-39.0*nuVtSqSum[0]*fr[2]+39.0*nuVtSqSum[0]*fl[2]+(27.71281292110204*fr[1]+27.71281292110204*fl[1])*nuVtSqSum[1]+(27.71281292110204*fr[0]+27.71281292110204*fl[0])*nuVtSqSum[0]); 
  incr2[4] = 0.002209708691207959*(121.2435565298214*nuVtSqSum[1]*fr[20]+121.2435565298214*nuVtSqSum[1]*fl[20]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fr[12]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fl[12]-174.4133022449836*nuVtSqSum[1]*fr[11]+174.4133022449836*nuVtSqSum[1]*fl[11]+135.5544171172596*nuVtSqSum[1]*fr[8]+135.5544171172596*nuVtSqSum[1]*fl[8]+123.9354670786373*nuVtSqSum[1]*fr[7]+123.9354670786373*nuVtSqSum[1]*fl[7]+((-174.4133022449836*nuVtSqSum[2])-195.0*nuVtSqSum[0])*fr[4]+(174.4133022449836*nuVtSqSum[2]+195.0*nuVtSqSum[0])*fl[4]+(123.9354670786373*fr[1]+123.9354670786373*fl[1])*nuVtSqSum[2]-195.0*nuVtSqSum[1]*fr[2]+195.0*nuVtSqSum[1]*fl[2]+(138.5640646055102*fr[0]+138.5640646055102*fl[0])*nuVtSqSum[1]+138.5640646055102*nuVtSqSum[0]*fr[1]+138.5640646055102*nuVtSqSum[0]*fl[1]); 
  incr2[6] = 0.0110485434560398*(27.11088342345192*nuVtSqSum[2]*fr[23]+27.11088342345192*nuVtSqSum[2]*fl[23]+27.11088342345192*nuVtSqSum[1]*fr[18]+27.11088342345192*nuVtSqSum[1]*fl[18]-39.0*nuVtSqSum[2]*fr[17]+39.0*nuVtSqSum[2]*fl[17]+27.11088342345192*nuVtSqSum[0]*fr[14]+27.11088342345192*nuVtSqSum[0]*fl[14]+27.71281292110204*nuVtSqSum[2]*fr[13]+27.71281292110204*nuVtSqSum[2]*fl[13]-39.0*nuVtSqSum[1]*fr[10]+39.0*nuVtSqSum[1]*fl[10]-39.0*nuVtSqSum[0]*fr[6]+39.0*nuVtSqSum[0]*fl[6]+27.71281292110204*nuVtSqSum[1]*fr[5]+27.71281292110204*nuVtSqSum[1]*fl[5]+27.71281292110204*nuVtSqSum[0]*fr[3]+27.71281292110204*nuVtSqSum[0]*fl[3]); 
  incr2[8] = -0.0110485434560398*(105.0*nuVtSqSum[2]*fr[20]+105.0*nuVtSqSum[2]*fl[20]+105.0*nuVtSqSum[1]*fr[12]+105.0*nuVtSqSum[1]*fl[12]-151.0463505020892*nuVtSqSum[2]*fr[11]+151.0463505020892*nuVtSqSum[2]*fl[11]+105.0*nuVtSqSum[0]*fr[8]+105.0*nuVtSqSum[0]*fl[8]+107.3312629199899*nuVtSqSum[2]*fr[7]+107.3312629199899*nuVtSqSum[2]*fl[7]-151.0463505020892*nuVtSqSum[1]*fr[4]+151.0463505020892*nuVtSqSum[1]*fl[4]-151.0463505020892*nuVtSqSum[0]*fr[2]+151.0463505020892*nuVtSqSum[0]*fl[2]+(107.3312629199899*fr[1]+107.3312629199899*fl[1])*nuVtSqSum[1]+(107.3312629199899*fr[0]+107.3312629199899*fl[0])*nuVtSqSum[0]); 
  incr2[10] = 0.002209708691207959*(121.2435565298214*nuVtSqSum[1]*fr[23]+121.2435565298214*nuVtSqSum[1]*fl[23]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fr[18]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fl[18]-174.4133022449836*nuVtSqSum[1]*fr[17]+174.4133022449836*nuVtSqSum[1]*fl[17]+135.5544171172596*nuVtSqSum[1]*fr[14]+135.5544171172596*nuVtSqSum[1]*fl[14]+123.9354670786373*nuVtSqSum[1]*fr[13]+123.9354670786373*nuVtSqSum[1]*fl[13]+((-174.4133022449836*nuVtSqSum[2])-195.0*nuVtSqSum[0])*fr[10]+(174.4133022449836*nuVtSqSum[2]+195.0*nuVtSqSum[0])*fl[10]-195.0*nuVtSqSum[1]*fr[6]+195.0*nuVtSqSum[1]*fl[6]+(123.9354670786373*nuVtSqSum[2]+138.5640646055102*nuVtSqSum[0])*fr[5]+(123.9354670786373*nuVtSqSum[2]+138.5640646055102*nuVtSqSum[0])*fl[5]+138.5640646055102*nuVtSqSum[1]*fr[3]+138.5640646055102*nuVtSqSum[1]*fl[3]); 
  incr2[11] = 3.156726701725656e-4*((606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fr[20]+(606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fl[20]+848.7048957087499*nuVtSqSum[1]*fr[12]+848.7048957087499*nuVtSqSum[1]*fl[12]+((-872.0665112249181*nuVtSqSum[2])-1365.0*nuVtSqSum[0])*fr[11]+(872.0665112249181*nuVtSqSum[2]+1365.0*nuVtSqSum[0])*fl[11]+948.8809198208172*nuVtSqSum[2]*fr[8]+948.8809198208172*nuVtSqSum[2]*fl[8]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fr[7]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fl[7]-1220.893115714885*nuVtSqSum[1]*fr[4]+1220.893115714885*nuVtSqSum[1]*fl[4]+((-1365.0*fr[2])+1365.0*fl[2]+969.9484522385712*fr[0]+969.9484522385712*fl[0])*nuVtSqSum[2]+(867.5482695504614*fr[1]+867.5482695504614*fl[1])*nuVtSqSum[1]); 
  incr2[12] = -0.0110485434560398*(93.91485505499116*nuVtSqSum[1]*fr[20]+93.91485505499116*nuVtSqSum[1]*fl[20]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fr[12]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fl[12]-135.0999629903724*nuVtSqSum[1]*fr[11]+135.0999629903724*nuVtSqSum[1]*fl[11]+105.0*nuVtSqSum[1]*fr[8]+105.0*nuVtSqSum[1]*fl[8]+96.0*nuVtSqSum[1]*fr[7]+96.0*nuVtSqSum[1]*fl[7]+((-135.0999629903724*nuVtSqSum[2])-151.0463505020892*nuVtSqSum[0])*fr[4]+(135.0999629903724*nuVtSqSum[2]+151.0463505020892*nuVtSqSum[0])*fl[4]+(96.0*fr[1]+96.0*fl[1])*nuVtSqSum[2]-151.0463505020892*nuVtSqSum[1]*fr[2]+151.0463505020892*nuVtSqSum[1]*fl[2]+(107.3312629199899*fr[0]+107.3312629199899*fl[0])*nuVtSqSum[1]+107.3312629199899*nuVtSqSum[0]*fr[1]+107.3312629199899*nuVtSqSum[0]*fl[1]); 
  incr2[14] = -0.0110485434560398*(105.0*nuVtSqSum[2]*fr[23]+105.0*nuVtSqSum[2]*fl[23]+105.0*nuVtSqSum[1]*fr[18]+105.0*nuVtSqSum[1]*fl[18]-151.0463505020892*nuVtSqSum[2]*fr[17]+151.0463505020892*nuVtSqSum[2]*fl[17]+105.0*nuVtSqSum[0]*fr[14]+105.0*nuVtSqSum[0]*fl[14]+107.3312629199899*nuVtSqSum[2]*fr[13]+107.3312629199899*nuVtSqSum[2]*fl[13]-151.0463505020892*nuVtSqSum[1]*fr[10]+151.0463505020892*nuVtSqSum[1]*fl[10]-151.0463505020892*nuVtSqSum[0]*fr[6]+151.0463505020892*nuVtSqSum[0]*fl[6]+107.3312629199899*nuVtSqSum[1]*fr[5]+107.3312629199899*nuVtSqSum[1]*fl[5]+107.3312629199899*nuVtSqSum[0]*fr[3]+107.3312629199899*nuVtSqSum[0]*fl[3]); 
  incr2[16] = 0.0110485434560398*(27.11088342345192*nuVtSqSum[2]*fr[26]+27.11088342345192*nuVtSqSum[2]*fl[26]+27.11088342345192*nuVtSqSum[1]*fr[25]+27.11088342345192*nuVtSqSum[1]*fl[25]-39.0*nuVtSqSum[2]*fr[24]+39.0*nuVtSqSum[2]*fl[24]+27.11088342345192*nuVtSqSum[0]*fr[22]+27.11088342345192*nuVtSqSum[0]*fl[22]+27.71281292110204*nuVtSqSum[2]*fr[21]+27.71281292110204*nuVtSqSum[2]*fl[21]-39.0*nuVtSqSum[1]*fr[19]+39.0*nuVtSqSum[1]*fl[19]-39.0*nuVtSqSum[0]*fr[16]+39.0*nuVtSqSum[0]*fl[16]+27.71281292110204*nuVtSqSum[1]*fr[15]+27.71281292110204*nuVtSqSum[1]*fl[15]+27.71281292110204*nuVtSqSum[0]*fr[9]+27.71281292110204*nuVtSqSum[0]*fl[9]); 
  incr2[17] = 3.156726701725656e-4*((606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fr[23]+(606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fl[23]+848.7048957087499*nuVtSqSum[1]*fr[18]+848.7048957087499*nuVtSqSum[1]*fl[18]+((-872.0665112249181*nuVtSqSum[2])-1365.0*nuVtSqSum[0])*fr[17]+(872.0665112249181*nuVtSqSum[2]+1365.0*nuVtSqSum[0])*fl[17]+948.8809198208172*nuVtSqSum[2]*fr[14]+948.8809198208172*nuVtSqSum[2]*fl[14]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fr[13]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fl[13]-1220.893115714885*nuVtSqSum[1]*fr[10]+1220.893115714885*nuVtSqSum[1]*fl[10]-1365.0*nuVtSqSum[2]*fr[6]+1365.0*nuVtSqSum[2]*fl[6]+867.5482695504614*nuVtSqSum[1]*fr[5]+867.5482695504614*nuVtSqSum[1]*fl[5]+969.9484522385712*nuVtSqSum[2]*fr[3]+969.9484522385712*nuVtSqSum[2]*fl[3]); 
  incr2[18] = -0.0110485434560398*(93.91485505499116*nuVtSqSum[1]*fr[23]+93.91485505499116*nuVtSqSum[1]*fl[23]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fr[18]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fl[18]-135.0999629903724*nuVtSqSum[1]*fr[17]+135.0999629903724*nuVtSqSum[1]*fl[17]+105.0*nuVtSqSum[1]*fr[14]+105.0*nuVtSqSum[1]*fl[14]+96.0*nuVtSqSum[1]*fr[13]+96.0*nuVtSqSum[1]*fl[13]+((-135.0999629903724*nuVtSqSum[2])-151.0463505020892*nuVtSqSum[0])*fr[10]+(135.0999629903724*nuVtSqSum[2]+151.0463505020892*nuVtSqSum[0])*fl[10]-151.0463505020892*nuVtSqSum[1]*fr[6]+151.0463505020892*nuVtSqSum[1]*fl[6]+(96.0*nuVtSqSum[2]+107.3312629199899*nuVtSqSum[0])*fr[5]+(96.0*nuVtSqSum[2]+107.3312629199899*nuVtSqSum[0])*fl[5]+107.3312629199899*nuVtSqSum[1]*fr[3]+107.3312629199899*nuVtSqSum[1]*fl[3]); 
  incr2[19] = 0.002209708691207959*(121.2435565298214*nuVtSqSum[1]*fr[26]+121.2435565298214*nuVtSqSum[1]*fl[26]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fr[25]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fl[25]-174.4133022449836*nuVtSqSum[1]*fr[24]+174.4133022449836*nuVtSqSum[1]*fl[24]+135.5544171172596*nuVtSqSum[1]*fr[22]+135.5544171172596*nuVtSqSum[1]*fl[22]+123.9354670786373*nuVtSqSum[1]*fr[21]+123.9354670786373*nuVtSqSum[1]*fl[21]+((-174.4133022449836*nuVtSqSum[2])-195.0*nuVtSqSum[0])*fr[19]+(174.4133022449836*nuVtSqSum[2]+195.0*nuVtSqSum[0])*fl[19]-195.0*nuVtSqSum[1]*fr[16]+195.0*nuVtSqSum[1]*fl[16]+(123.9354670786373*nuVtSqSum[2]+138.5640646055102*nuVtSqSum[0])*fr[15]+(123.9354670786373*nuVtSqSum[2]+138.5640646055102*nuVtSqSum[0])*fl[15]+138.5640646055102*nuVtSqSum[1]*fr[9]+138.5640646055102*nuVtSqSum[1]*fl[9]); 
  incr2[20] = -0.001578363350862828*((469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fr[20]+(469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fl[20]+657.4039853849382*nuVtSqSum[1]*fr[12]+657.4039853849382*nuVtSqSum[1]*fl[12]+((-675.499814951862*nuVtSqSum[2])-1057.324453514625*nuVtSqSum[0])*fr[11]+(675.499814951862*nuVtSqSum[2]+1057.324453514625*nuVtSqSum[0])*fl[11]+735.0*nuVtSqSum[2]*fr[8]+735.0*nuVtSqSum[2]*fl[8]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fr[7]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fl[7]-945.6997409326069*nuVtSqSum[1]*fr[4]+945.6997409326069*nuVtSqSum[1]*fl[4]+((-1057.324453514625*fr[2])+1057.324453514625*fl[2]+751.3188404399293*fr[0]+751.3188404399293*fl[0])*nuVtSqSum[2]+(672.0*fr[1]+672.0*fl[1])*nuVtSqSum[1]); 
  incr2[22] = -0.0110485434560398*(105.0*nuVtSqSum[2]*fr[26]+105.0*nuVtSqSum[2]*fl[26]+105.0*nuVtSqSum[1]*fr[25]+105.0*nuVtSqSum[1]*fl[25]-151.0463505020892*nuVtSqSum[2]*fr[24]+151.0463505020892*nuVtSqSum[2]*fl[24]+105.0*nuVtSqSum[0]*fr[22]+105.0*nuVtSqSum[0]*fl[22]+107.3312629199899*nuVtSqSum[2]*fr[21]+107.3312629199899*nuVtSqSum[2]*fl[21]-151.0463505020892*nuVtSqSum[1]*fr[19]+151.0463505020892*nuVtSqSum[1]*fl[19]-151.0463505020892*nuVtSqSum[0]*fr[16]+151.0463505020892*nuVtSqSum[0]*fl[16]+107.3312629199899*nuVtSqSum[1]*fr[15]+107.3312629199899*nuVtSqSum[1]*fl[15]+107.3312629199899*nuVtSqSum[0]*fr[9]+107.3312629199899*nuVtSqSum[0]*fl[9]); 
  incr2[23] = -0.001578363350862828*((469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fr[23]+(469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fl[23]+657.4039853849382*nuVtSqSum[1]*fr[18]+657.4039853849382*nuVtSqSum[1]*fl[18]+((-675.499814951862*nuVtSqSum[2])-1057.324453514625*nuVtSqSum[0])*fr[17]+(675.499814951862*nuVtSqSum[2]+1057.324453514625*nuVtSqSum[0])*fl[17]+735.0*nuVtSqSum[2]*fr[14]+735.0*nuVtSqSum[2]*fl[14]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fr[13]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fl[13]-945.6997409326069*nuVtSqSum[1]*fr[10]+945.6997409326069*nuVtSqSum[1]*fl[10]-1057.324453514625*nuVtSqSum[2]*fr[6]+1057.324453514625*nuVtSqSum[2]*fl[6]+672.0*nuVtSqSum[1]*fr[5]+672.0*nuVtSqSum[1]*fl[5]+751.3188404399293*nuVtSqSum[2]*fr[3]+751.3188404399293*nuVtSqSum[2]*fl[3]); 
  incr2[24] = 3.156726701725656e-4*((606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fr[26]+(606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fl[26]+848.7048957087499*nuVtSqSum[1]*fr[25]+848.7048957087499*nuVtSqSum[1]*fl[25]+((-872.0665112249181*nuVtSqSum[2])-1365.0*nuVtSqSum[0])*fr[24]+(872.0665112249181*nuVtSqSum[2]+1365.0*nuVtSqSum[0])*fl[24]+948.8809198208172*nuVtSqSum[2]*fr[22]+948.8809198208172*nuVtSqSum[2]*fl[22]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fr[21]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fl[21]-1220.893115714885*nuVtSqSum[1]*fr[19]+1220.893115714885*nuVtSqSum[1]*fl[19]-1365.0*nuVtSqSum[2]*fr[16]+1365.0*nuVtSqSum[2]*fl[16]+867.5482695504614*nuVtSqSum[1]*fr[15]+867.5482695504614*nuVtSqSum[1]*fl[15]+969.9484522385712*nuVtSqSum[2]*fr[9]+969.9484522385712*nuVtSqSum[2]*fl[9]); 
  incr2[25] = -0.0110485434560398*(93.91485505499116*nuVtSqSum[1]*fr[26]+93.91485505499116*nuVtSqSum[1]*fl[26]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fr[25]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fl[25]-135.0999629903724*nuVtSqSum[1]*fr[24]+135.0999629903724*nuVtSqSum[1]*fl[24]+105.0*nuVtSqSum[1]*fr[22]+105.0*nuVtSqSum[1]*fl[22]+96.0*nuVtSqSum[1]*fr[21]+96.0*nuVtSqSum[1]*fl[21]+((-135.0999629903724*nuVtSqSum[2])-151.0463505020892*nuVtSqSum[0])*fr[19]+(135.0999629903724*nuVtSqSum[2]+151.0463505020892*nuVtSqSum[0])*fl[19]-151.0463505020892*nuVtSqSum[1]*fr[16]+151.0463505020892*nuVtSqSum[1]*fl[16]+(96.0*nuVtSqSum[2]+107.3312629199899*nuVtSqSum[0])*fr[15]+(96.0*nuVtSqSum[2]+107.3312629199899*nuVtSqSum[0])*fl[15]+107.3312629199899*nuVtSqSum[1]*fr[9]+107.3312629199899*nuVtSqSum[1]*fl[9]); 
  incr2[26] = -0.001578363350862828*((469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fr[26]+(469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fl[26]+657.4039853849382*nuVtSqSum[1]*fr[25]+657.4039853849382*nuVtSqSum[1]*fl[25]+((-675.499814951862*nuVtSqSum[2])-1057.324453514625*nuVtSqSum[0])*fr[24]+(675.499814951862*nuVtSqSum[2]+1057.324453514625*nuVtSqSum[0])*fl[24]+735.0*nuVtSqSum[2]*fr[22]+735.0*nuVtSqSum[2]*fl[22]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fr[21]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fl[21]-945.6997409326069*nuVtSqSum[1]*fr[19]+945.6997409326069*nuVtSqSum[1]*fl[19]-1057.324453514625*nuVtSqSum[2]*fr[16]+1057.324453514625*nuVtSqSum[2]*fl[16]+672.0*nuVtSqSum[1]*fr[15]+672.0*nuVtSqSum[1]*fl[15]+751.3188404399293*nuVtSqSum[2]*fr[9]+751.3188404399293*nuVtSqSum[2]*fl[9]); 


  Gdiff[0] = 0.01767766952966368*(53.66563145999496*nuVtSqSum[2]*fr[20]-53.66563145999496*nuVtSqSum[2]*fl[20]+53.66563145999495*nuVtSqSum[1]*fr[12]-53.66563145999495*nuVtSqSum[1]*fl[12]-95.26279441628826*nuVtSqSum[2]*fr[11]-95.26279441628826*nuVtSqSum[2]*fl[11]+53.66563145999496*nuVtSqSum[0]*fr[8]-53.66563145999496*nuVtSqSum[0]*fl[8]+75.0*nuVtSqSum[2]*fr[7]-75.0*nuVtSqSum[2]*fl[7]-95.26279441628824*nuVtSqSum[1]*fr[4]-95.26279441628824*nuVtSqSum[1]*fl[4]-95.26279441628824*nuVtSqSum[0]*fr[2]-95.26279441628824*nuVtSqSum[0]*fl[2]+(75.0*fr[1]-75.0*fl[1])*nuVtSqSum[1]+(75.0*fr[0]-75.0*fl[0])*nuVtSqSum[0]); 
  Gdiff[1] = 0.003535533905932736*(240.0*nuVtSqSum[1]*fr[20]-240.0*nuVtSqSum[1]*fl[20]+(240.0*nuVtSqSum[2]+268.3281572999747*nuVtSqSum[0])*fr[12]+((-240.0*nuVtSqSum[2])-268.3281572999747*nuVtSqSum[0])*fl[12]-426.0281680828159*nuVtSqSum[1]*fr[11]-426.0281680828159*nuVtSqSum[1]*fl[11]+268.3281572999748*nuVtSqSum[1]*fr[8]-268.3281572999748*nuVtSqSum[1]*fl[8]+335.4101966249685*nuVtSqSum[1]*fr[7]-335.4101966249685*nuVtSqSum[1]*fl[7]+((-426.0281680828159*nuVtSqSum[2])-476.3139720814412*nuVtSqSum[0])*fr[4]+((-426.0281680828159*nuVtSqSum[2])-476.3139720814412*nuVtSqSum[0])*fl[4]+(335.4101966249685*fr[1]-335.4101966249685*fl[1])*nuVtSqSum[2]-476.3139720814412*nuVtSqSum[1]*fr[2]-476.3139720814412*nuVtSqSum[1]*fl[2]+(375.0*fr[0]-375.0*fl[0])*nuVtSqSum[1]+375.0*nuVtSqSum[0]*fr[1]-375.0*nuVtSqSum[0]*fl[1]); 
  Gdiff[3] = 0.01767766952966368*(53.66563145999496*nuVtSqSum[2]*fr[23]-53.66563145999496*nuVtSqSum[2]*fl[23]+53.66563145999496*nuVtSqSum[1]*fr[18]-53.66563145999496*nuVtSqSum[1]*fl[18]-95.26279441628824*nuVtSqSum[2]*fr[17]-95.26279441628824*nuVtSqSum[2]*fl[17]+53.66563145999495*nuVtSqSum[0]*fr[14]-53.66563145999495*nuVtSqSum[0]*fl[14]+75.00000000000001*nuVtSqSum[2]*fr[13]-75.00000000000001*nuVtSqSum[2]*fl[13]-95.26279441628824*nuVtSqSum[1]*fr[10]-95.26279441628824*nuVtSqSum[1]*fl[10]-95.26279441628824*nuVtSqSum[0]*fr[6]-95.26279441628824*nuVtSqSum[0]*fl[6]+75.0*nuVtSqSum[1]*fr[5]-75.0*nuVtSqSum[1]*fl[5]+75.0*nuVtSqSum[0]*fr[3]-75.0*nuVtSqSum[0]*fl[3]); 
  Gdiff[5] = 0.01767766952966368*(48.0*nuVtSqSum[1]*fr[23]-48.0*nuVtSqSum[1]*fl[23]+(48.0*nuVtSqSum[2]+53.66563145999496*nuVtSqSum[0])*fr[18]+((-48.0*nuVtSqSum[2])-53.66563145999496*nuVtSqSum[0])*fl[18]-85.20563361656318*nuVtSqSum[1]*fr[17]-85.20563361656318*nuVtSqSum[1]*fl[17]+53.66563145999495*nuVtSqSum[1]*fr[14]-53.66563145999495*nuVtSqSum[1]*fl[14]+67.08203932499369*nuVtSqSum[1]*fr[13]-67.08203932499369*nuVtSqSum[1]*fl[13]+((-85.20563361656318*nuVtSqSum[2])-95.26279441628824*nuVtSqSum[0])*fr[10]+((-85.20563361656318*nuVtSqSum[2])-95.26279441628824*nuVtSqSum[0])*fl[10]-95.26279441628824*nuVtSqSum[1]*fr[6]-95.26279441628824*nuVtSqSum[1]*fl[6]+(67.0820393249937*nuVtSqSum[2]+75.0*nuVtSqSum[0])*fr[5]+((-67.0820393249937*nuVtSqSum[2])-75.0*nuVtSqSum[0])*fl[5]+75.0*nuVtSqSum[1]*fr[3]-75.0*nuVtSqSum[1]*fl[3]); 
  Gdiff[7] = 5.050762722761051e-4*((1200.0*nuVtSqSum[2]+1878.297101099824*nuVtSqSum[0])*fr[20]+((-1200.0*nuVtSqSum[2])-1878.297101099824*nuVtSqSum[0])*fl[20]+1680.0*nuVtSqSum[1]*fr[12]-1680.0*nuVtSqSum[1]*fl[12]+((-2130.140840414079*nuVtSqSum[2])-3334.19780457009*nuVtSqSum[0])*fr[11]+((-2130.140840414079*nuVtSqSum[2])-3334.19780457009*nuVtSqSum[0])*fl[11]+1878.297101099824*nuVtSqSum[2]*fr[8]-1878.297101099824*nuVtSqSum[2]*fl[8]+(1677.050983124843*nuVtSqSum[2]+2625.0*nuVtSqSum[0])*fr[7]+((-1677.050983124843*nuVtSqSum[2])-2625.0*nuVtSqSum[0])*fl[7]-2982.197176579711*nuVtSqSum[1]*fr[4]-2982.197176579711*nuVtSqSum[1]*fl[4]+((-3334.197804570088*fr[2])-3334.197804570088*fl[2]+2625.0*fr[0]-2625.0*fl[0])*nuVtSqSum[2]+(2347.87137637478*fr[1]-2347.87137637478*fl[1])*nuVtSqSum[1]); 
  Gdiff[9] = 0.01767766952966368*(53.66563145999496*nuVtSqSum[2]*fr[26]-53.66563145999496*nuVtSqSum[2]*fl[26]+53.66563145999496*nuVtSqSum[1]*fr[25]-53.66563145999496*nuVtSqSum[1]*fl[25]-95.26279441628824*nuVtSqSum[2]*fr[24]-95.26279441628824*nuVtSqSum[2]*fl[24]+53.66563145999496*nuVtSqSum[0]*fr[22]-53.66563145999496*nuVtSqSum[0]*fl[22]+75.0*nuVtSqSum[2]*fr[21]-75.0*nuVtSqSum[2]*fl[21]-95.26279441628824*nuVtSqSum[1]*fr[19]-95.26279441628824*nuVtSqSum[1]*fl[19]-95.26279441628826*nuVtSqSum[0]*fr[16]-95.26279441628826*nuVtSqSum[0]*fl[16]+75.00000000000001*nuVtSqSum[1]*fr[15]-75.00000000000001*nuVtSqSum[1]*fl[15]+75.0*nuVtSqSum[0]*fr[9]-75.0*nuVtSqSum[0]*fl[9]); 
  Gdiff[13] = 5.050762722761051e-4*((1200.0*nuVtSqSum[2]+1878.297101099823*nuVtSqSum[0])*fr[23]+((-1200.0*nuVtSqSum[2])-1878.297101099823*nuVtSqSum[0])*fl[23]+1680.0*nuVtSqSum[1]*fr[18]-1680.0*nuVtSqSum[1]*fl[18]+((-2130.140840414079*nuVtSqSum[2])-3334.19780457009*nuVtSqSum[0])*fr[17]+((-2130.140840414079*nuVtSqSum[2])-3334.19780457009*nuVtSqSum[0])*fl[17]+1878.297101099824*nuVtSqSum[2]*fr[14]-1878.297101099824*nuVtSqSum[2]*fl[14]+(1677.050983124843*nuVtSqSum[2]+2625.0*nuVtSqSum[0])*fr[13]+((-1677.050983124843*nuVtSqSum[2])-2625.0*nuVtSqSum[0])*fl[13]-2982.197176579711*nuVtSqSum[1]*fr[10]-2982.197176579711*nuVtSqSum[1]*fl[10]-3334.19780457009*nuVtSqSum[2]*fr[6]-3334.19780457009*nuVtSqSum[2]*fl[6]+2347.871376374779*nuVtSqSum[1]*fr[5]-2347.871376374779*nuVtSqSum[1]*fl[5]+2625.0*nuVtSqSum[2]*fr[3]-2625.0*nuVtSqSum[2]*fl[3]); 
  Gdiff[15] = 0.003535533905932736*(240.0*nuVtSqSum[1]*fr[26]-240.0*nuVtSqSum[1]*fl[26]+(240.0*nuVtSqSum[2]+268.3281572999747*nuVtSqSum[0])*fr[25]+((-240.0*nuVtSqSum[2])-268.3281572999747*nuVtSqSum[0])*fl[25]-426.0281680828159*nuVtSqSum[1]*fr[24]-426.0281680828159*nuVtSqSum[1]*fl[24]+268.3281572999747*nuVtSqSum[1]*fr[22]-268.3281572999747*nuVtSqSum[1]*fl[22]+335.4101966249685*nuVtSqSum[1]*fr[21]-335.4101966249685*nuVtSqSum[1]*fl[21]+((-426.0281680828159*nuVtSqSum[2])-476.3139720814414*nuVtSqSum[0])*fr[19]+((-426.0281680828159*nuVtSqSum[2])-476.3139720814414*nuVtSqSum[0])*fl[19]-476.3139720814412*nuVtSqSum[1]*fr[16]-476.3139720814412*nuVtSqSum[1]*fl[16]+(335.4101966249685*nuVtSqSum[2]+375.0*nuVtSqSum[0])*fr[15]+((-335.4101966249685*nuVtSqSum[2])-375.0*nuVtSqSum[0])*fl[15]+375.0000000000001*nuVtSqSum[1]*fr[9]-375.0000000000001*nuVtSqSum[1]*fl[9]); 
  Gdiff[21] = 0.002525381361380526*((240.0*nuVtSqSum[2]+375.6594202199647*nuVtSqSum[0])*fr[26]+((-240.0*nuVtSqSum[2])-375.6594202199647*nuVtSqSum[0])*fl[26]+336.0*nuVtSqSum[1]*fr[25]-336.0*nuVtSqSum[1]*fl[25]+((-426.0281680828159*nuVtSqSum[2])-666.8395609140177*nuVtSqSum[0])*fr[24]+((-426.0281680828159*nuVtSqSum[2])-666.8395609140177*nuVtSqSum[0])*fl[24]+375.6594202199647*nuVtSqSum[2]*fr[22]-375.6594202199647*nuVtSqSum[2]*fl[22]+(335.4101966249685*nuVtSqSum[2]+525.0*nuVtSqSum[0])*fr[21]+((-335.4101966249685*nuVtSqSum[2])-525.0*nuVtSqSum[0])*fl[21]-596.4394353159422*nuVtSqSum[1]*fr[19]-596.4394353159422*nuVtSqSum[1]*fl[19]-666.8395609140179*nuVtSqSum[2]*fr[16]-666.8395609140179*nuVtSqSum[2]*fl[16]+469.5742752749558*nuVtSqSum[1]*fr[15]-469.5742752749558*nuVtSqSum[1]*fl[15]+525.0*nuVtSqSum[2]*fr[9]-525.0*nuVtSqSum[2]*fl[9]); 

  Ghat[0] = Gdiff[0]*rdv2L+alphaDrag[2]*(0.7905694150420947*favg[20]+0.6123724356957944*favg[11]+0.3535533905932737*favg[7])+alphaDrag[1]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])-1.118033988749895*fjump[8]+alphaDrag[0]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[20]+0.5477225575051661*favg[11]+0.7905694150420947*favg[8]+0.3162277660168379*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-1.118033988749895*fjump[12]+alphaDrag[0]*(0.7905694150420948*favg[12]+0.6123724356957944*favg[4]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])-0.8660254037844386*fjump[4]-0.5*fjump[1]; 
  Ghat[3] = Gdiff[3]*rdv2L+alphaDrag[2]*(0.7905694150420947*favg[23]+0.6123724356957944*favg[17]+0.3535533905932737*favg[13])+alphaDrag[1]*(0.7905694150420947*favg[18]+0.6123724356957944*favg[10]+0.3535533905932737*favg[5])-1.118033988749895*fjump[14]+alphaDrag[0]*(0.7905694150420948*favg[14]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-0.8660254037844386*fjump[6]-0.5*fjump[3]; 
  Ghat[5] = Gdiff[5]*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[23]+0.5477225575051661*favg[17]+0.7905694150420948*favg[14]+0.3162277660168379*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])-1.118033988749895*fjump[18]+alphaDrag[0]*(0.7905694150420947*favg[18]+0.6123724356957944*favg[10]+0.3535533905932737*favg[5])+alphaDrag[2]*(0.7071067811865475*favg[18]+0.5477225575051661*favg[10]+0.3162277660168379*favg[5])-0.8660254037844386*fjump[10]-0.5*fjump[5]; 
  Ghat[7] = Gdiff[7]*rdv2L-1.118033988749895*fjump[20]+alphaDrag[0]*(0.7905694150420947*favg[20]+0.6123724356957944*favg[11]+0.3535533905932737*favg[7])+alphaDrag[2]*(0.5050762722761053*favg[20]+0.3912303982179757*favg[11]+0.7905694150420947*favg[8]+0.2258769757263128*favg[7]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+alphaDrag[1]*(0.7071067811865475*favg[12]+0.5477225575051661*favg[4]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[11]-0.5*fjump[7]; 
  Ghat[9] = Gdiff[9]*rdv2L+alphaDrag[2]*(0.7905694150420947*favg[26]+0.6123724356957944*favg[24]+0.3535533905932737*favg[21])+alphaDrag[1]*(0.7905694150420947*favg[25]+0.6123724356957944*favg[19]+0.3535533905932737*favg[15])-1.118033988749895*fjump[22]+alphaDrag[0]*(0.7905694150420947*favg[22]+0.6123724356957944*favg[16]+0.3535533905932737*favg[9])-0.8660254037844387*fjump[16]-0.5*fjump[9]; 
  Ghat[13] = Gdiff[13]*rdv2L-1.118033988749895*fjump[23]+alphaDrag[0]*(0.7905694150420948*favg[23]+0.6123724356957944*favg[17]+0.3535533905932737*favg[13])+alphaDrag[2]*(0.5050762722761054*favg[23]+0.3912303982179757*favg[17]+0.7905694150420947*favg[14]+0.2258769757263128*favg[13]+0.6123724356957944*favg[6]+0.3535533905932737*favg[3])+alphaDrag[1]*(0.7071067811865475*favg[18]+0.5477225575051661*favg[10]+0.3162277660168379*favg[5])-0.8660254037844387*fjump[17]-0.5*fjump[13]; 
  Ghat[15] = Gdiff[15]*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[26]+0.5477225575051661*favg[24]+0.7905694150420948*favg[22]+0.3162277660168379*favg[21]+0.6123724356957944*favg[16]+0.3535533905932737*favg[9])-1.118033988749895*fjump[25]+alphaDrag[0]*(0.7905694150420948*favg[25]+0.6123724356957944*favg[19]+0.3535533905932737*favg[15])+alphaDrag[2]*(0.7071067811865475*favg[25]+0.5477225575051661*favg[19]+0.3162277660168379*favg[15])-0.8660254037844387*fjump[19]-0.5*fjump[15]; 
  Ghat[21] = Gdiff[21]*rdv2L-1.118033988749895*fjump[26]+alphaDrag[0]*(0.7905694150420947*favg[26]+0.6123724356957944*favg[24]+0.3535533905932737*favg[21])+alphaDrag[2]*(0.5050762722761053*favg[26]+0.3912303982179757*favg[24]+0.7905694150420947*favg[22]+0.2258769757263128*favg[21]+0.6123724356957944*favg[16]+0.3535533905932737*favg[9])+alphaDrag[1]*(0.7071067811865475*favg[25]+0.5477225575051661*favg[19]+0.3162277660168379*favg[15])-0.8660254037844386*fjump[24]-0.5*fjump[21]; 

  double incr1[27]; 
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
  incr1[10] = 0.8660254037844386*Ghat[5]; 
  incr1[11] = 0.8660254037844387*Ghat[7]; 
  incr1[12] = -1.118033988749895*Ghat[1]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = -1.118033988749895*Ghat[3]; 
  incr1[15] = -0.5*Ghat[15]; 
  incr1[16] = 0.8660254037844387*Ghat[9]; 
  incr1[17] = 0.8660254037844387*Ghat[13]; 
  incr1[18] = -1.118033988749895*Ghat[5]; 
  incr1[19] = 0.8660254037844387*Ghat[15]; 
  incr1[20] = -1.118033988749895*Ghat[7]; 
  incr1[21] = -0.5*Ghat[21]; 
  incr1[22] = -1.118033988749895*Ghat[9]; 
  incr1[23] = -1.118033988749895*Ghat[13]; 
  incr1[24] = 0.8660254037844386*Ghat[21]; 
  incr1[25] = -1.118033988749895*Ghat[15]; 
  incr1[26] = -1.118033988749895*Ghat[21]; 

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
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr2[11]*rdvSq4R+incr1[11]*rdv2R; 
  outr[12] += incr2[12]*rdvSq4R+incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr1[15]*rdv2R; 
  outr[16] += incr2[16]*rdvSq4R+incr1[16]*rdv2R; 
  outr[17] += incr2[17]*rdvSq4R+incr1[17]*rdv2R; 
  outr[18] += incr2[18]*rdvSq4R+incr1[18]*rdv2R; 
  outr[19] += incr2[19]*rdvSq4R+incr1[19]*rdv2R; 
  outr[20] += incr2[20]*rdvSq4R+incr1[20]*rdv2R; 
  outr[21] += incr1[21]*rdv2R; 
  outr[22] += incr2[22]*rdvSq4R+incr1[22]*rdv2R; 
  outr[23] += incr2[23]*rdvSq4R+incr1[23]*rdv2R; 
  outr[24] += incr2[24]*rdvSq4R+incr1[24]*rdv2R; 
  outr[25] += incr2[25]*rdvSq4R+incr1[25]*rdv2R; 
  outr[26] += incr2[26]*rdvSq4R+incr1[26]*rdv2R; 

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
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += incr1[11]*rdv2L-1.0*incr2[11]*rdvSq4L; 
  outl[12] += incr2[12]*rdvSq4L-1.0*incr1[12]*rdv2L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += incr2[14]*rdvSq4L-1.0*incr1[14]*rdv2L; 
  outl[15] += -1.0*incr1[15]*rdv2L; 
  outl[16] += incr1[16]*rdv2L-1.0*incr2[16]*rdvSq4L; 
  outl[17] += incr1[17]*rdv2L-1.0*incr2[17]*rdvSq4L; 
  outl[18] += incr2[18]*rdvSq4L-1.0*incr1[18]*rdv2L; 
  outl[19] += incr1[19]*rdv2L-1.0*incr2[19]*rdvSq4L; 
  outl[20] += incr2[20]*rdvSq4L-1.0*incr1[20]*rdv2L; 
  outl[21] += -1.0*incr1[21]*rdv2L; 
  outl[22] += incr2[22]*rdvSq4L-1.0*incr1[22]*rdv2L; 
  outl[23] += incr2[23]*rdvSq4L-1.0*incr1[23]*rdv2L; 
  outl[24] += incr1[24]*rdv2L-1.0*incr2[24]*rdvSq4L; 
  outl[25] += incr2[25]*rdvSq4L-1.0*incr1[25]*rdv2L; 
  outl[26] += incr2[26]*rdvSq4L-1.0*incr1[26]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUx[2])/nuSum-(0.7071067811865475*sumNuUx[0])/nuSum+wl[1]); 
} 
double VmLBOconstNuSurf1x2vTensor_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
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

  double favg[27]; 
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
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = -1*fr[17]+fl[17]; 
  favg[18] = -1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  favg[21] = 1*fr[21]+fl[21]; 
  favg[22] = 1*fr[22]+fl[22]; 
  favg[23] = -1*fr[23]+fl[23]; 
  favg[24] = 1*fr[24]+fl[24]; 
  favg[25] = 1*fr[25]+fl[25]; 
  favg[26] = 1*fr[26]+fl[26]; 

  double fjump[27]; 
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
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(-1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(-1*fr[14])); 
  fjump[15] = nuSum*vMuMidMax*(fl[15]-(1*fr[15])); 
  fjump[16] = nuSum*vMuMidMax*(fl[16]-(1*fr[16])); 
  fjump[17] = nuSum*vMuMidMax*(fl[17]-(-1*fr[17])); 
  fjump[18] = nuSum*vMuMidMax*(fl[18]-(-1*fr[18])); 
  fjump[19] = nuSum*vMuMidMax*(fl[19]-(1*fr[19])); 
  fjump[20] = nuSum*vMuMidMax*(fl[20]-(1*fr[20])); 
  fjump[21] = nuSum*vMuMidMax*(fl[21]-(1*fr[21])); 
  fjump[22] = nuSum*vMuMidMax*(fl[22]-(1*fr[22])); 
  fjump[23] = nuSum*vMuMidMax*(fl[23]-(-1*fr[23])); 
  fjump[24] = nuSum*vMuMidMax*(fl[24]-(1*fr[24])); 
  fjump[25] = nuSum*vMuMidMax*(fl[25]-(1*fr[25])); 
  fjump[26] = nuSum*vMuMidMax*(fl[26]-(1*fr[26])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[2]*nuSum+0.7071067811865475*dxvl[2]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 

  double Gdiff[27]; 
  double Ghat[27]; 
  double incr2[27]; 


  incr2[3] = 0.0110485434560398*(27.11088342345192*nuVtSqSum[2]*fr[21]+27.11088342345192*nuVtSqSum[2]*fl[21]+27.11088342345192*nuVtSqSum[1]*fr[15]+27.11088342345192*nuVtSqSum[1]*fl[15]-39.0*nuVtSqSum[2]*fr[13]+39.0*nuVtSqSum[2]*fl[13]+27.11088342345192*nuVtSqSum[0]*fr[9]+27.11088342345192*nuVtSqSum[0]*fl[9]+27.71281292110204*nuVtSqSum[2]*fr[7]+27.71281292110204*nuVtSqSum[2]*fl[7]-39.0*nuVtSqSum[1]*fr[5]+39.0*nuVtSqSum[1]*fl[5]-39.0*nuVtSqSum[0]*fr[3]+39.0*nuVtSqSum[0]*fl[3]+(27.71281292110204*fr[1]+27.71281292110204*fl[1])*nuVtSqSum[1]+(27.71281292110204*fr[0]+27.71281292110204*fl[0])*nuVtSqSum[0]); 
  incr2[5] = 0.002209708691207959*(121.2435565298214*nuVtSqSum[1]*fr[21]+121.2435565298214*nuVtSqSum[1]*fl[21]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fr[15]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fl[15]-174.4133022449836*nuVtSqSum[1]*fr[13]+174.4133022449836*nuVtSqSum[1]*fl[13]+135.5544171172596*nuVtSqSum[1]*fr[9]+135.5544171172596*nuVtSqSum[1]*fl[9]+123.9354670786373*nuVtSqSum[1]*fr[7]+123.9354670786373*nuVtSqSum[1]*fl[7]+((-174.4133022449836*nuVtSqSum[2])-195.0*nuVtSqSum[0])*fr[5]+(174.4133022449836*nuVtSqSum[2]+195.0*nuVtSqSum[0])*fl[5]-195.0*nuVtSqSum[1]*fr[3]+195.0*nuVtSqSum[1]*fl[3]+(123.9354670786373*fr[1]+123.9354670786373*fl[1])*nuVtSqSum[2]+(138.5640646055102*fr[0]+138.5640646055102*fl[0])*nuVtSqSum[1]+138.5640646055102*nuVtSqSum[0]*fr[1]+138.5640646055102*nuVtSqSum[0]*fl[1]); 
  incr2[6] = 0.0110485434560398*(27.11088342345192*nuVtSqSum[2]*fr[24]+27.11088342345192*nuVtSqSum[2]*fl[24]+27.11088342345192*nuVtSqSum[1]*fr[19]+27.11088342345192*nuVtSqSum[1]*fl[19]-39.0*nuVtSqSum[2]*fr[17]+39.0*nuVtSqSum[2]*fl[17]+27.11088342345192*nuVtSqSum[0]*fr[16]+27.11088342345192*nuVtSqSum[0]*fl[16]+27.71281292110204*nuVtSqSum[2]*fr[11]+27.71281292110204*nuVtSqSum[2]*fl[11]-39.0*nuVtSqSum[1]*fr[10]+39.0*nuVtSqSum[1]*fl[10]-39.0*nuVtSqSum[0]*fr[6]+39.0*nuVtSqSum[0]*fl[6]+27.71281292110204*nuVtSqSum[1]*fr[4]+27.71281292110204*nuVtSqSum[1]*fl[4]+27.71281292110204*nuVtSqSum[0]*fr[2]+27.71281292110204*nuVtSqSum[0]*fl[2]); 
  incr2[9] = -0.0110485434560398*(105.0*nuVtSqSum[2]*fr[21]+105.0*nuVtSqSum[2]*fl[21]+105.0*nuVtSqSum[1]*fr[15]+105.0*nuVtSqSum[1]*fl[15]-151.0463505020892*nuVtSqSum[2]*fr[13]+151.0463505020892*nuVtSqSum[2]*fl[13]+105.0*nuVtSqSum[0]*fr[9]+105.0*nuVtSqSum[0]*fl[9]+107.3312629199899*nuVtSqSum[2]*fr[7]+107.3312629199899*nuVtSqSum[2]*fl[7]-151.0463505020892*nuVtSqSum[1]*fr[5]+151.0463505020892*nuVtSqSum[1]*fl[5]-151.0463505020892*nuVtSqSum[0]*fr[3]+151.0463505020892*nuVtSqSum[0]*fl[3]+(107.3312629199899*fr[1]+107.3312629199899*fl[1])*nuVtSqSum[1]+(107.3312629199899*fr[0]+107.3312629199899*fl[0])*nuVtSqSum[0]); 
  incr2[10] = 0.002209708691207959*(121.2435565298214*nuVtSqSum[1]*fr[24]+121.2435565298214*nuVtSqSum[1]*fl[24]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fr[19]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fl[19]-174.4133022449836*nuVtSqSum[1]*fr[17]+174.4133022449836*nuVtSqSum[1]*fl[17]+135.5544171172596*nuVtSqSum[1]*fr[16]+135.5544171172596*nuVtSqSum[1]*fl[16]+123.9354670786373*nuVtSqSum[1]*fr[11]+123.9354670786373*nuVtSqSum[1]*fl[11]+((-174.4133022449836*nuVtSqSum[2])-195.0*nuVtSqSum[0])*fr[10]+(174.4133022449836*nuVtSqSum[2]+195.0*nuVtSqSum[0])*fl[10]-195.0*nuVtSqSum[1]*fr[6]+195.0*nuVtSqSum[1]*fl[6]+(123.9354670786373*nuVtSqSum[2]+138.5640646055102*nuVtSqSum[0])*fr[4]+(123.9354670786373*nuVtSqSum[2]+138.5640646055102*nuVtSqSum[0])*fl[4]+138.5640646055102*nuVtSqSum[1]*fr[2]+138.5640646055102*nuVtSqSum[1]*fl[2]); 
  incr2[13] = 3.156726701725656e-4*((606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fr[21]+(606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fl[21]+848.7048957087499*nuVtSqSum[1]*fr[15]+848.7048957087499*nuVtSqSum[1]*fl[15]+((-872.0665112249181*nuVtSqSum[2])-1365.0*nuVtSqSum[0])*fr[13]+(872.0665112249181*nuVtSqSum[2]+1365.0*nuVtSqSum[0])*fl[13]+948.8809198208172*nuVtSqSum[2]*fr[9]+948.8809198208172*nuVtSqSum[2]*fl[9]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fr[7]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fl[7]-1220.893115714885*nuVtSqSum[1]*fr[5]+1220.893115714885*nuVtSqSum[1]*fl[5]-1365.0*nuVtSqSum[2]*fr[3]+1365.0*nuVtSqSum[2]*fl[3]+(969.9484522385712*fr[0]+969.9484522385712*fl[0])*nuVtSqSum[2]+(867.5482695504614*fr[1]+867.5482695504614*fl[1])*nuVtSqSum[1]); 
  incr2[14] = 0.0110485434560398*(27.11088342345192*nuVtSqSum[2]*fr[26]+27.11088342345192*nuVtSqSum[2]*fl[26]+27.11088342345192*nuVtSqSum[1]*fr[25]+27.11088342345192*nuVtSqSum[1]*fl[25]-39.0*nuVtSqSum[2]*fr[23]+39.0*nuVtSqSum[2]*fl[23]+27.11088342345192*nuVtSqSum[0]*fr[22]+27.11088342345192*nuVtSqSum[0]*fl[22]+27.71281292110204*nuVtSqSum[2]*fr[20]+27.71281292110204*nuVtSqSum[2]*fl[20]-39.0*nuVtSqSum[1]*fr[18]+39.0*nuVtSqSum[1]*fl[18]-39.0*nuVtSqSum[0]*fr[14]+39.0*nuVtSqSum[0]*fl[14]+27.71281292110204*nuVtSqSum[1]*fr[12]+27.71281292110204*nuVtSqSum[1]*fl[12]+27.71281292110204*nuVtSqSum[0]*fr[8]+27.71281292110204*nuVtSqSum[0]*fl[8]); 
  incr2[15] = -0.0110485434560398*(93.91485505499116*nuVtSqSum[1]*fr[21]+93.91485505499116*nuVtSqSum[1]*fl[21]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fr[15]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fl[15]-135.0999629903724*nuVtSqSum[1]*fr[13]+135.0999629903724*nuVtSqSum[1]*fl[13]+105.0*nuVtSqSum[1]*fr[9]+105.0*nuVtSqSum[1]*fl[9]+96.0*nuVtSqSum[1]*fr[7]+96.0*nuVtSqSum[1]*fl[7]+((-135.0999629903724*nuVtSqSum[2])-151.0463505020892*nuVtSqSum[0])*fr[5]+(135.0999629903724*nuVtSqSum[2]+151.0463505020892*nuVtSqSum[0])*fl[5]-151.0463505020892*nuVtSqSum[1]*fr[3]+151.0463505020892*nuVtSqSum[1]*fl[3]+(96.0*fr[1]+96.0*fl[1])*nuVtSqSum[2]+(107.3312629199899*fr[0]+107.3312629199899*fl[0])*nuVtSqSum[1]+107.3312629199899*nuVtSqSum[0]*fr[1]+107.3312629199899*nuVtSqSum[0]*fl[1]); 
  incr2[16] = -0.0110485434560398*(105.0*nuVtSqSum[2]*fr[24]+105.0*nuVtSqSum[2]*fl[24]+105.0*nuVtSqSum[1]*fr[19]+105.0*nuVtSqSum[1]*fl[19]-151.0463505020892*nuVtSqSum[2]*fr[17]+151.0463505020892*nuVtSqSum[2]*fl[17]+105.0*nuVtSqSum[0]*fr[16]+105.0*nuVtSqSum[0]*fl[16]+107.3312629199899*nuVtSqSum[2]*fr[11]+107.3312629199899*nuVtSqSum[2]*fl[11]-151.0463505020892*nuVtSqSum[1]*fr[10]+151.0463505020892*nuVtSqSum[1]*fl[10]-151.0463505020892*nuVtSqSum[0]*fr[6]+151.0463505020892*nuVtSqSum[0]*fl[6]+107.3312629199899*nuVtSqSum[1]*fr[4]+107.3312629199899*nuVtSqSum[1]*fl[4]+107.3312629199899*nuVtSqSum[0]*fr[2]+107.3312629199899*nuVtSqSum[0]*fl[2]); 
  incr2[17] = 3.156726701725656e-4*((606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fr[24]+(606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fl[24]+848.7048957087499*nuVtSqSum[1]*fr[19]+848.7048957087499*nuVtSqSum[1]*fl[19]+((-872.0665112249181*nuVtSqSum[2])-1365.0*nuVtSqSum[0])*fr[17]+(872.0665112249181*nuVtSqSum[2]+1365.0*nuVtSqSum[0])*fl[17]+948.8809198208172*nuVtSqSum[2]*fr[16]+948.8809198208172*nuVtSqSum[2]*fl[16]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fr[11]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fl[11]-1220.893115714885*nuVtSqSum[1]*fr[10]+1220.893115714885*nuVtSqSum[1]*fl[10]-1365.0*nuVtSqSum[2]*fr[6]+1365.0*nuVtSqSum[2]*fl[6]+867.5482695504614*nuVtSqSum[1]*fr[4]+867.5482695504614*nuVtSqSum[1]*fl[4]+(969.9484522385712*fr[2]+969.9484522385712*fl[2])*nuVtSqSum[2]); 
  incr2[18] = 0.002209708691207959*(121.2435565298214*nuVtSqSum[1]*fr[26]+121.2435565298214*nuVtSqSum[1]*fl[26]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fr[25]+(121.2435565298214*nuVtSqSum[2]+135.5544171172596*nuVtSqSum[0])*fl[25]-174.4133022449836*nuVtSqSum[1]*fr[23]+174.4133022449836*nuVtSqSum[1]*fl[23]+135.5544171172596*nuVtSqSum[1]*fr[22]+135.5544171172596*nuVtSqSum[1]*fl[22]+123.9354670786373*nuVtSqSum[1]*fr[20]+123.9354670786373*nuVtSqSum[1]*fl[20]+((-174.4133022449836*nuVtSqSum[2])-195.0*nuVtSqSum[0])*fr[18]+(174.4133022449836*nuVtSqSum[2]+195.0*nuVtSqSum[0])*fl[18]-195.0*nuVtSqSum[1]*fr[14]+195.0*nuVtSqSum[1]*fl[14]+(123.9354670786373*nuVtSqSum[2]+138.5640646055102*nuVtSqSum[0])*fr[12]+(123.9354670786373*nuVtSqSum[2]+138.5640646055102*nuVtSqSum[0])*fl[12]+138.5640646055102*nuVtSqSum[1]*fr[8]+138.5640646055102*nuVtSqSum[1]*fl[8]); 
  incr2[19] = -0.0110485434560398*(93.91485505499116*nuVtSqSum[1]*fr[24]+93.91485505499116*nuVtSqSum[1]*fl[24]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fr[19]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fl[19]-135.0999629903724*nuVtSqSum[1]*fr[17]+135.0999629903724*nuVtSqSum[1]*fl[17]+105.0*nuVtSqSum[1]*fr[16]+105.0*nuVtSqSum[1]*fl[16]+96.0*nuVtSqSum[1]*fr[11]+96.0*nuVtSqSum[1]*fl[11]+((-135.0999629903724*nuVtSqSum[2])-151.0463505020892*nuVtSqSum[0])*fr[10]+(135.0999629903724*nuVtSqSum[2]+151.0463505020892*nuVtSqSum[0])*fl[10]-151.0463505020892*nuVtSqSum[1]*fr[6]+151.0463505020892*nuVtSqSum[1]*fl[6]+(96.0*nuVtSqSum[2]+107.3312629199899*nuVtSqSum[0])*fr[4]+(96.0*nuVtSqSum[2]+107.3312629199899*nuVtSqSum[0])*fl[4]+107.3312629199899*nuVtSqSum[1]*fr[2]+107.3312629199899*nuVtSqSum[1]*fl[2]); 
  incr2[21] = -0.001578363350862828*((469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fr[21]+(469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fl[21]+657.4039853849382*nuVtSqSum[1]*fr[15]+657.4039853849382*nuVtSqSum[1]*fl[15]+((-675.499814951862*nuVtSqSum[2])-1057.324453514625*nuVtSqSum[0])*fr[13]+(675.499814951862*nuVtSqSum[2]+1057.324453514625*nuVtSqSum[0])*fl[13]+735.0*nuVtSqSum[2]*fr[9]+735.0*nuVtSqSum[2]*fl[9]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fr[7]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fl[7]-945.6997409326069*nuVtSqSum[1]*fr[5]+945.6997409326069*nuVtSqSum[1]*fl[5]-1057.324453514625*nuVtSqSum[2]*fr[3]+1057.324453514625*nuVtSqSum[2]*fl[3]+(751.3188404399293*fr[0]+751.3188404399293*fl[0])*nuVtSqSum[2]+(672.0*fr[1]+672.0*fl[1])*nuVtSqSum[1]); 
  incr2[22] = -0.0110485434560398*(105.0*nuVtSqSum[2]*fr[26]+105.0*nuVtSqSum[2]*fl[26]+105.0*nuVtSqSum[1]*fr[25]+105.0*nuVtSqSum[1]*fl[25]-151.0463505020892*nuVtSqSum[2]*fr[23]+151.0463505020892*nuVtSqSum[2]*fl[23]+105.0*nuVtSqSum[0]*fr[22]+105.0*nuVtSqSum[0]*fl[22]+107.3312629199899*nuVtSqSum[2]*fr[20]+107.3312629199899*nuVtSqSum[2]*fl[20]-151.0463505020892*nuVtSqSum[1]*fr[18]+151.0463505020892*nuVtSqSum[1]*fl[18]-151.0463505020892*nuVtSqSum[0]*fr[14]+151.0463505020892*nuVtSqSum[0]*fl[14]+107.3312629199899*nuVtSqSum[1]*fr[12]+107.3312629199899*nuVtSqSum[1]*fl[12]+107.3312629199899*nuVtSqSum[0]*fr[8]+107.3312629199899*nuVtSqSum[0]*fl[8]); 
  incr2[23] = 3.156726701725656e-4*((606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fr[26]+(606.217782649107*nuVtSqSum[2]+948.8809198208172*nuVtSqSum[0])*fl[26]+848.7048957087499*nuVtSqSum[1]*fr[25]+848.7048957087499*nuVtSqSum[1]*fl[25]+((-872.0665112249181*nuVtSqSum[2])-1365.0*nuVtSqSum[0])*fr[23]+(872.0665112249181*nuVtSqSum[2]+1365.0*nuVtSqSum[0])*fl[23]+948.8809198208172*nuVtSqSum[2]*fr[22]+948.8809198208172*nuVtSqSum[2]*fl[22]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fr[20]+(619.6773353931868*nuVtSqSum[2]+969.9484522385712*nuVtSqSum[0])*fl[20]-1220.893115714885*nuVtSqSum[1]*fr[18]+1220.893115714885*nuVtSqSum[1]*fl[18]-1365.0*nuVtSqSum[2]*fr[14]+1365.0*nuVtSqSum[2]*fl[14]+867.5482695504614*nuVtSqSum[1]*fr[12]+867.5482695504614*nuVtSqSum[1]*fl[12]+969.9484522385712*nuVtSqSum[2]*fr[8]+969.9484522385712*nuVtSqSum[2]*fl[8]); 
  incr2[24] = -0.001578363350862828*((469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fr[24]+(469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fl[24]+657.4039853849382*nuVtSqSum[1]*fr[19]+657.4039853849382*nuVtSqSum[1]*fl[19]+((-675.499814951862*nuVtSqSum[2])-1057.324453514625*nuVtSqSum[0])*fr[17]+(675.499814951862*nuVtSqSum[2]+1057.324453514625*nuVtSqSum[0])*fl[17]+735.0*nuVtSqSum[2]*fr[16]+735.0*nuVtSqSum[2]*fl[16]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fr[11]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fl[11]-945.6997409326069*nuVtSqSum[1]*fr[10]+945.6997409326069*nuVtSqSum[1]*fl[10]-1057.324453514625*nuVtSqSum[2]*fr[6]+1057.324453514625*nuVtSqSum[2]*fl[6]+672.0*nuVtSqSum[1]*fr[4]+672.0*nuVtSqSum[1]*fl[4]+(751.3188404399293*fr[2]+751.3188404399293*fl[2])*nuVtSqSum[2]); 
  incr2[25] = -0.0110485434560398*(93.91485505499116*nuVtSqSum[1]*fr[26]+93.91485505499116*nuVtSqSum[1]*fl[26]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fr[25]+(93.91485505499116*nuVtSqSum[2]+105.0*nuVtSqSum[0])*fl[25]-135.0999629903724*nuVtSqSum[1]*fr[23]+135.0999629903724*nuVtSqSum[1]*fl[23]+105.0*nuVtSqSum[1]*fr[22]+105.0*nuVtSqSum[1]*fl[22]+96.0*nuVtSqSum[1]*fr[20]+96.0*nuVtSqSum[1]*fl[20]+((-135.0999629903724*nuVtSqSum[2])-151.0463505020892*nuVtSqSum[0])*fr[18]+(135.0999629903724*nuVtSqSum[2]+151.0463505020892*nuVtSqSum[0])*fl[18]-151.0463505020892*nuVtSqSum[1]*fr[14]+151.0463505020892*nuVtSqSum[1]*fl[14]+(96.0*nuVtSqSum[2]+107.3312629199899*nuVtSqSum[0])*fr[12]+(96.0*nuVtSqSum[2]+107.3312629199899*nuVtSqSum[0])*fl[12]+107.3312629199899*nuVtSqSum[1]*fr[8]+107.3312629199899*nuVtSqSum[1]*fl[8]); 
  incr2[26] = -0.001578363350862828*((469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fr[26]+(469.5742752749559*nuVtSqSum[2]+735.0*nuVtSqSum[0])*fl[26]+657.4039853849382*nuVtSqSum[1]*fr[25]+657.4039853849382*nuVtSqSum[1]*fl[25]+((-675.499814951862*nuVtSqSum[2])-1057.324453514625*nuVtSqSum[0])*fr[23]+(675.499814951862*nuVtSqSum[2]+1057.324453514625*nuVtSqSum[0])*fl[23]+735.0*nuVtSqSum[2]*fr[22]+735.0*nuVtSqSum[2]*fl[22]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fr[20]+(480.0*nuVtSqSum[2]+751.3188404399293*nuVtSqSum[0])*fl[20]-945.6997409326069*nuVtSqSum[1]*fr[18]+945.6997409326069*nuVtSqSum[1]*fl[18]-1057.324453514625*nuVtSqSum[2]*fr[14]+1057.324453514625*nuVtSqSum[2]*fl[14]+672.0*nuVtSqSum[1]*fr[12]+672.0*nuVtSqSum[1]*fl[12]+751.3188404399293*nuVtSqSum[2]*fr[8]+751.3188404399293*nuVtSqSum[2]*fl[8]); 


  Gdiff[0] = 0.01767766952966368*(53.66563145999496*nuVtSqSum[2]*fr[21]-53.66563145999496*nuVtSqSum[2]*fl[21]+53.66563145999495*nuVtSqSum[1]*fr[15]-53.66563145999495*nuVtSqSum[1]*fl[15]-95.26279441628826*nuVtSqSum[2]*fr[13]-95.26279441628826*nuVtSqSum[2]*fl[13]+53.66563145999496*nuVtSqSum[0]*fr[9]-53.66563145999496*nuVtSqSum[0]*fl[9]+75.0*nuVtSqSum[2]*fr[7]-75.0*nuVtSqSum[2]*fl[7]-95.26279441628824*nuVtSqSum[1]*fr[5]-95.26279441628824*nuVtSqSum[1]*fl[5]-95.26279441628824*nuVtSqSum[0]*fr[3]-95.26279441628824*nuVtSqSum[0]*fl[3]+(75.0*fr[1]-75.0*fl[1])*nuVtSqSum[1]+(75.0*fr[0]-75.0*fl[0])*nuVtSqSum[0]); 
  Gdiff[1] = 0.003535533905932736*(240.0*nuVtSqSum[1]*fr[21]-240.0*nuVtSqSum[1]*fl[21]+(240.0*nuVtSqSum[2]+268.3281572999747*nuVtSqSum[0])*fr[15]+((-240.0*nuVtSqSum[2])-268.3281572999747*nuVtSqSum[0])*fl[15]-426.0281680828159*nuVtSqSum[1]*fr[13]-426.0281680828159*nuVtSqSum[1]*fl[13]+268.3281572999748*nuVtSqSum[1]*fr[9]-268.3281572999748*nuVtSqSum[1]*fl[9]+335.4101966249685*nuVtSqSum[1]*fr[7]-335.4101966249685*nuVtSqSum[1]*fl[7]+((-426.0281680828159*nuVtSqSum[2])-476.3139720814412*nuVtSqSum[0])*fr[5]+((-426.0281680828159*nuVtSqSum[2])-476.3139720814412*nuVtSqSum[0])*fl[5]-476.3139720814412*nuVtSqSum[1]*fr[3]-476.3139720814412*nuVtSqSum[1]*fl[3]+(335.4101966249685*fr[1]-335.4101966249685*fl[1])*nuVtSqSum[2]+(375.0*fr[0]-375.0*fl[0])*nuVtSqSum[1]+375.0*nuVtSqSum[0]*fr[1]-375.0*nuVtSqSum[0]*fl[1]); 
  Gdiff[2] = 0.01767766952966368*(53.66563145999496*nuVtSqSum[2]*fr[24]-53.66563145999496*nuVtSqSum[2]*fl[24]+53.66563145999496*nuVtSqSum[1]*fr[19]-53.66563145999496*nuVtSqSum[1]*fl[19]-95.26279441628824*nuVtSqSum[2]*fr[17]-95.26279441628824*nuVtSqSum[2]*fl[17]+53.66563145999495*nuVtSqSum[0]*fr[16]-53.66563145999495*nuVtSqSum[0]*fl[16]+75.00000000000001*nuVtSqSum[2]*fr[11]-75.00000000000001*nuVtSqSum[2]*fl[11]-95.26279441628824*nuVtSqSum[1]*fr[10]-95.26279441628824*nuVtSqSum[1]*fl[10]-95.26279441628824*nuVtSqSum[0]*fr[6]-95.26279441628824*nuVtSqSum[0]*fl[6]+75.0*nuVtSqSum[1]*fr[4]-75.0*nuVtSqSum[1]*fl[4]+75.0*nuVtSqSum[0]*fr[2]-75.0*nuVtSqSum[0]*fl[2]); 
  Gdiff[4] = 0.01767766952966368*(48.0*nuVtSqSum[1]*fr[24]-48.0*nuVtSqSum[1]*fl[24]+(48.0*nuVtSqSum[2]+53.66563145999496*nuVtSqSum[0])*fr[19]+((-48.0*nuVtSqSum[2])-53.66563145999496*nuVtSqSum[0])*fl[19]-85.20563361656318*nuVtSqSum[1]*fr[17]-85.20563361656318*nuVtSqSum[1]*fl[17]+53.66563145999495*nuVtSqSum[1]*fr[16]-53.66563145999495*nuVtSqSum[1]*fl[16]+67.08203932499369*nuVtSqSum[1]*fr[11]-67.08203932499369*nuVtSqSum[1]*fl[11]+((-85.20563361656318*nuVtSqSum[2])-95.26279441628824*nuVtSqSum[0])*fr[10]+((-85.20563361656318*nuVtSqSum[2])-95.26279441628824*nuVtSqSum[0])*fl[10]-95.26279441628824*nuVtSqSum[1]*fr[6]-95.26279441628824*nuVtSqSum[1]*fl[6]+(67.0820393249937*nuVtSqSum[2]+75.0*nuVtSqSum[0])*fr[4]+((-67.0820393249937*nuVtSqSum[2])-75.0*nuVtSqSum[0])*fl[4]+75.0*nuVtSqSum[1]*fr[2]-75.0*nuVtSqSum[1]*fl[2]); 
  Gdiff[7] = 5.050762722761051e-4*((1200.0*nuVtSqSum[2]+1878.297101099824*nuVtSqSum[0])*fr[21]+((-1200.0*nuVtSqSum[2])-1878.297101099824*nuVtSqSum[0])*fl[21]+1680.0*nuVtSqSum[1]*fr[15]-1680.0*nuVtSqSum[1]*fl[15]+((-2130.140840414079*nuVtSqSum[2])-3334.19780457009*nuVtSqSum[0])*fr[13]+((-2130.140840414079*nuVtSqSum[2])-3334.19780457009*nuVtSqSum[0])*fl[13]+1878.297101099824*nuVtSqSum[2]*fr[9]-1878.297101099824*nuVtSqSum[2]*fl[9]+(1677.050983124843*nuVtSqSum[2]+2625.0*nuVtSqSum[0])*fr[7]+((-1677.050983124843*nuVtSqSum[2])-2625.0*nuVtSqSum[0])*fl[7]-2982.197176579711*nuVtSqSum[1]*fr[5]-2982.197176579711*nuVtSqSum[1]*fl[5]-3334.197804570088*nuVtSqSum[2]*fr[3]-3334.197804570088*nuVtSqSum[2]*fl[3]+(2625.0*fr[0]-2625.0*fl[0])*nuVtSqSum[2]+(2347.87137637478*fr[1]-2347.87137637478*fl[1])*nuVtSqSum[1]); 
  Gdiff[8] = 0.01767766952966368*(53.66563145999496*nuVtSqSum[2]*fr[26]-53.66563145999496*nuVtSqSum[2]*fl[26]+53.66563145999496*nuVtSqSum[1]*fr[25]-53.66563145999496*nuVtSqSum[1]*fl[25]-95.26279441628824*nuVtSqSum[2]*fr[23]-95.26279441628824*nuVtSqSum[2]*fl[23]+53.66563145999496*nuVtSqSum[0]*fr[22]-53.66563145999496*nuVtSqSum[0]*fl[22]+75.0*nuVtSqSum[2]*fr[20]-75.0*nuVtSqSum[2]*fl[20]-95.26279441628824*nuVtSqSum[1]*fr[18]-95.26279441628824*nuVtSqSum[1]*fl[18]-95.26279441628826*nuVtSqSum[0]*fr[14]-95.26279441628826*nuVtSqSum[0]*fl[14]+75.00000000000001*nuVtSqSum[1]*fr[12]-75.00000000000001*nuVtSqSum[1]*fl[12]+75.0*nuVtSqSum[0]*fr[8]-75.0*nuVtSqSum[0]*fl[8]); 
  Gdiff[11] = 5.050762722761051e-4*((1200.0*nuVtSqSum[2]+1878.297101099823*nuVtSqSum[0])*fr[24]+((-1200.0*nuVtSqSum[2])-1878.297101099823*nuVtSqSum[0])*fl[24]+1680.0*nuVtSqSum[1]*fr[19]-1680.0*nuVtSqSum[1]*fl[19]+((-2130.140840414079*nuVtSqSum[2])-3334.19780457009*nuVtSqSum[0])*fr[17]+((-2130.140840414079*nuVtSqSum[2])-3334.19780457009*nuVtSqSum[0])*fl[17]+1878.297101099824*nuVtSqSum[2]*fr[16]-1878.297101099824*nuVtSqSum[2]*fl[16]+(1677.050983124843*nuVtSqSum[2]+2625.0*nuVtSqSum[0])*fr[11]+((-1677.050983124843*nuVtSqSum[2])-2625.0*nuVtSqSum[0])*fl[11]-2982.197176579711*nuVtSqSum[1]*fr[10]-2982.197176579711*nuVtSqSum[1]*fl[10]-3334.19780457009*nuVtSqSum[2]*fr[6]-3334.19780457009*nuVtSqSum[2]*fl[6]+2347.871376374779*nuVtSqSum[1]*fr[4]-2347.871376374779*nuVtSqSum[1]*fl[4]+(2625.0*fr[2]-2625.0*fl[2])*nuVtSqSum[2]); 
  Gdiff[12] = 0.003535533905932736*(240.0*nuVtSqSum[1]*fr[26]-240.0*nuVtSqSum[1]*fl[26]+(240.0*nuVtSqSum[2]+268.3281572999747*nuVtSqSum[0])*fr[25]+((-240.0*nuVtSqSum[2])-268.3281572999747*nuVtSqSum[0])*fl[25]-426.0281680828159*nuVtSqSum[1]*fr[23]-426.0281680828159*nuVtSqSum[1]*fl[23]+268.3281572999747*nuVtSqSum[1]*fr[22]-268.3281572999747*nuVtSqSum[1]*fl[22]+335.4101966249685*nuVtSqSum[1]*fr[20]-335.4101966249685*nuVtSqSum[1]*fl[20]+((-426.0281680828159*nuVtSqSum[2])-476.3139720814414*nuVtSqSum[0])*fr[18]+((-426.0281680828159*nuVtSqSum[2])-476.3139720814414*nuVtSqSum[0])*fl[18]-476.3139720814412*nuVtSqSum[1]*fr[14]-476.3139720814412*nuVtSqSum[1]*fl[14]+(335.4101966249685*nuVtSqSum[2]+375.0*nuVtSqSum[0])*fr[12]+((-335.4101966249685*nuVtSqSum[2])-375.0*nuVtSqSum[0])*fl[12]+375.0000000000001*nuVtSqSum[1]*fr[8]-375.0000000000001*nuVtSqSum[1]*fl[8]); 
  Gdiff[20] = 0.002525381361380526*((240.0*nuVtSqSum[2]+375.6594202199647*nuVtSqSum[0])*fr[26]+((-240.0*nuVtSqSum[2])-375.6594202199647*nuVtSqSum[0])*fl[26]+336.0*nuVtSqSum[1]*fr[25]-336.0*nuVtSqSum[1]*fl[25]+((-426.0281680828159*nuVtSqSum[2])-666.8395609140177*nuVtSqSum[0])*fr[23]+((-426.0281680828159*nuVtSqSum[2])-666.8395609140177*nuVtSqSum[0])*fl[23]+375.6594202199647*nuVtSqSum[2]*fr[22]-375.6594202199647*nuVtSqSum[2]*fl[22]+(335.4101966249685*nuVtSqSum[2]+525.0*nuVtSqSum[0])*fr[20]+((-335.4101966249685*nuVtSqSum[2])-525.0*nuVtSqSum[0])*fl[20]-596.4394353159422*nuVtSqSum[1]*fr[18]-596.4394353159422*nuVtSqSum[1]*fl[18]-666.8395609140179*nuVtSqSum[2]*fr[14]-666.8395609140179*nuVtSqSum[2]*fl[14]+469.5742752749558*nuVtSqSum[1]*fr[12]-469.5742752749558*nuVtSqSum[1]*fl[12]+525.0*nuVtSqSum[2]*fr[8]-525.0*nuVtSqSum[2]*fl[8]); 

  Ghat[0] = Gdiff[0]*rdv2L+alphaDrag[2]*(0.7905694150420947*favg[21]+0.6123724356957944*favg[13]+0.3535533905932737*favg[7])+alphaDrag[1]*(0.7905694150420948*favg[15]+0.6123724356957944*favg[5]+0.3535533905932737*favg[1])-1.118033988749895*fjump[9]+alphaDrag[0]*(0.7905694150420947*favg[9]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[3]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[21]+0.5477225575051661*favg[13]+0.7905694150420947*favg[9]+0.3162277660168379*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])-1.118033988749895*fjump[15]+alphaDrag[0]*(0.7905694150420948*favg[15]+0.6123724356957944*favg[5]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.7071067811865475*favg[15]+0.5477225575051661*favg[5]+0.3162277660168379*favg[1])-0.8660254037844386*fjump[5]-0.5*fjump[1]; 
  Ghat[2] = Gdiff[2]*rdv2L+alphaDrag[2]*(0.7905694150420947*favg[24]+0.6123724356957944*favg[17]+0.3535533905932737*favg[11])+alphaDrag[1]*(0.7905694150420947*favg[19]+0.6123724356957944*favg[10]+0.3535533905932737*favg[4])-1.118033988749895*fjump[16]+alphaDrag[0]*(0.7905694150420948*favg[16]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-0.8660254037844386*fjump[6]-0.5*fjump[2]; 
  Ghat[4] = Gdiff[4]*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[24]+0.5477225575051661*favg[17]+0.7905694150420948*favg[16]+0.3162277660168379*favg[11]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])-1.118033988749895*fjump[19]+alphaDrag[0]*(0.7905694150420947*favg[19]+0.6123724356957944*favg[10]+0.3535533905932737*favg[4])+alphaDrag[2]*(0.7071067811865475*favg[19]+0.5477225575051661*favg[10]+0.3162277660168379*favg[4])-0.8660254037844386*fjump[10]-0.5*fjump[4]; 
  Ghat[7] = Gdiff[7]*rdv2L-1.118033988749895*fjump[21]+alphaDrag[0]*(0.7905694150420947*favg[21]+0.6123724356957944*favg[13]+0.3535533905932737*favg[7])+alphaDrag[2]*(0.5050762722761053*favg[21]+0.3912303982179757*favg[13]+0.7905694150420947*favg[9]+0.2258769757263128*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[0])+alphaDrag[1]*(0.7071067811865475*favg[15]+0.5477225575051661*favg[5]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[13]-0.5*fjump[7]; 
  Ghat[8] = Gdiff[8]*rdv2L+alphaDrag[2]*(0.7905694150420947*favg[26]+0.6123724356957944*favg[23]+0.3535533905932737*favg[20])+alphaDrag[1]*(0.7905694150420947*favg[25]+0.6123724356957944*favg[18]+0.3535533905932737*favg[12])-1.118033988749895*fjump[22]+alphaDrag[0]*(0.7905694150420947*favg[22]+0.6123724356957944*favg[14]+0.3535533905932737*favg[8])-0.8660254037844387*fjump[14]-0.5*fjump[8]; 
  Ghat[11] = Gdiff[11]*rdv2L-1.118033988749895*fjump[24]+alphaDrag[0]*(0.7905694150420948*favg[24]+0.6123724356957944*favg[17]+0.3535533905932737*favg[11])+alphaDrag[2]*(0.5050762722761054*favg[24]+0.3912303982179757*favg[17]+0.7905694150420947*favg[16]+0.2258769757263128*favg[11]+0.6123724356957944*favg[6]+0.3535533905932737*favg[2])+alphaDrag[1]*(0.7071067811865475*favg[19]+0.5477225575051661*favg[10]+0.3162277660168379*favg[4])-0.8660254037844387*fjump[17]-0.5*fjump[11]; 
  Ghat[12] = Gdiff[12]*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[26]+0.5477225575051661*favg[23]+0.7905694150420948*favg[22]+0.3162277660168379*favg[20]+0.6123724356957944*favg[14]+0.3535533905932737*favg[8])-1.118033988749895*fjump[25]+alphaDrag[0]*(0.7905694150420948*favg[25]+0.6123724356957944*favg[18]+0.3535533905932737*favg[12])+alphaDrag[2]*(0.7071067811865475*favg[25]+0.5477225575051661*favg[18]+0.3162277660168379*favg[12])-0.8660254037844387*fjump[18]-0.5*fjump[12]; 
  Ghat[20] = Gdiff[20]*rdv2L-1.118033988749895*fjump[26]+alphaDrag[0]*(0.7905694150420947*favg[26]+0.6123724356957944*favg[23]+0.3535533905932737*favg[20])+alphaDrag[2]*(0.5050762722761053*favg[26]+0.3912303982179757*favg[23]+0.7905694150420947*favg[22]+0.2258769757263128*favg[20]+0.6123724356957944*favg[14]+0.3535533905932737*favg[8])+alphaDrag[1]*(0.7071067811865475*favg[25]+0.5477225575051661*favg[18]+0.3162277660168379*favg[12])-0.8660254037844386*fjump[23]-0.5*fjump[20]; 

  double incr1[27]; 
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
  incr1[10] = 0.8660254037844386*Ghat[4]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = 0.8660254037844387*Ghat[7]; 
  incr1[14] = 0.8660254037844387*Ghat[8]; 
  incr1[15] = -1.118033988749895*Ghat[1]; 
  incr1[16] = -1.118033988749895*Ghat[2]; 
  incr1[17] = 0.8660254037844387*Ghat[11]; 
  incr1[18] = 0.8660254037844387*Ghat[12]; 
  incr1[19] = -1.118033988749895*Ghat[4]; 
  incr1[20] = -0.5*Ghat[20]; 
  incr1[21] = -1.118033988749895*Ghat[7]; 
  incr1[22] = -1.118033988749895*Ghat[8]; 
  incr1[23] = 0.8660254037844386*Ghat[20]; 
  incr1[24] = -1.118033988749895*Ghat[11]; 
  incr1[25] = -1.118033988749895*Ghat[12]; 
  incr1[26] = -1.118033988749895*Ghat[20]; 

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
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 
  outr[16] += incr2[16]*rdvSq4R+incr1[16]*rdv2R; 
  outr[17] += incr2[17]*rdvSq4R+incr1[17]*rdv2R; 
  outr[18] += incr2[18]*rdvSq4R+incr1[18]*rdv2R; 
  outr[19] += incr2[19]*rdvSq4R+incr1[19]*rdv2R; 
  outr[20] += incr1[20]*rdv2R; 
  outr[21] += incr2[21]*rdvSq4R+incr1[21]*rdv2R; 
  outr[22] += incr2[22]*rdvSq4R+incr1[22]*rdv2R; 
  outr[23] += incr2[23]*rdvSq4R+incr1[23]*rdv2R; 
  outr[24] += incr2[24]*rdvSq4R+incr1[24]*rdv2R; 
  outr[25] += incr2[25]*rdvSq4R+incr1[25]*rdv2R; 
  outr[26] += incr2[26]*rdvSq4R+incr1[26]*rdv2R; 

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
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += incr1[13]*rdv2L-1.0*incr2[13]*rdvSq4L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr2[15]*rdvSq4L-1.0*incr1[15]*rdv2L; 
  outl[16] += incr2[16]*rdvSq4L-1.0*incr1[16]*rdv2L; 
  outl[17] += incr1[17]*rdv2L-1.0*incr2[17]*rdvSq4L; 
  outl[18] += incr1[18]*rdv2L-1.0*incr2[18]*rdvSq4L; 
  outl[19] += incr2[19]*rdvSq4L-1.0*incr1[19]*rdv2L; 
  outl[20] += -1.0*incr1[20]*rdv2L; 
  outl[21] += incr2[21]*rdvSq4L-1.0*incr1[21]*rdv2L; 
  outl[22] += incr2[22]*rdvSq4L-1.0*incr1[22]*rdv2L; 
  outl[23] += incr1[23]*rdv2L-1.0*incr2[23]*rdvSq4L; 
  outl[24] += incr2[24]*rdvSq4L-1.0*incr1[24]*rdv2L; 
  outl[25] += incr2[25]*rdvSq4L-1.0*incr1[25]*rdv2L; 
  outl[26] += incr2[26]*rdvSq4L-1.0*incr1[26]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUy[2])/nuSum-(0.7071067811865475*sumNuUy[0])/nuSum+wl[2]); 
} 
