#include <VmLBOModDecl.h> 
double VmLBOconstNuUpwindSurf1x1vTensor_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:          Cell-center coordinates. 
  // dxv[2]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[3]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf[3]; 
  alphaDrSurf[0] = -0.7071067811865475*(((-2.0*wl[1])-1.0*dxvl[1])*nuSum+1.414213562373095*sumNuUx[0]); 
  alphaDrSurf[1] = -1.0*sumNuUx[1]; 
  alphaDrSurf[2] = -1.0*sumNuUx[2]; 

  double fUpwindQuad[3];
  fUpwindQuad[0] = 0.5*((fr[8]-(fl[8]-(1.5*fl[7]-1.5*fr[7])+0.7745966692414833*(fr[6]+fl[6])-(1.118033988749895*fr[5]-1.118033988749895*fl[5])-0.4472135954999579*fr[4]+0.4472135954999579*fl[4]-(1.161895003862225*(fr[3]+fl[3])-0.8660254037844386*(fr[2]+fl[2]))+0.6708203932499369*fr[1]-0.6708203932499369*fl[1]-(0.5*fr[0]-0.5*fl[0])))*sgn(0.6324555320336759*alphaDrSurf[2]-0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0])+fr[8]+fl[8]-1.5*(fr[7]+fl[7])-0.7745966692414833*fr[6]+0.7745966692414833*fl[6]+1.118033988749895*(fr[5]+fl[5])+0.4472135954999579*(fr[4]+fl[4])+1.161895003862225*fr[3]-1.161895003862225*fl[3]-0.8660254037844386*fr[2]+0.8660254037844386*fl[2]-0.6708203932499369*(fr[1]+fl[1])+0.5*(fr[0]+fl[0])); 
  fUpwindQuad[1] = 0.5*((-1.25*(fr[8]+fl[8]))+0.9682458365518543*fr[6]-0.9682458365518543*fl[6]+1.118033988749895*(fr[5]+fl[5])-0.5590169943749475*(fr[4]+fl[4])-0.8660254037844386*fr[2]+0.8660254037844386*fl[2]+0.5*(fr[0]+fl[0]))-0.5*(1.25*fr[8]-1.25*fl[8]-0.9682458365518543*(fr[6]+fl[6])-1.118033988749895*fr[5]+1.118033988749895*fl[5]+0.5590169943749475*fr[4]-0.5590169943749475*fl[4]+0.8660254037844386*(fr[2]+fl[2])-0.5*fr[0]+0.5*fl[0])*sgn(0.7071067811865475*alphaDrSurf[0]-0.7905694150420947*alphaDrSurf[2]); 
  fUpwindQuad[2] = 0.5*((fr[8]-(fl[8]-(1.5*fr[7]-1.5*fl[7])+0.7745966692414833*(fr[6]+fl[6])-(1.118033988749895*fr[5]-1.118033988749895*fl[5])-0.4472135954999579*fr[4]+0.4472135954999579*fl[4]-((-1.161895003862225*(fr[3]+fl[3]))-0.8660254037844386*(fr[2]+fl[2]))-0.6708203932499369*fr[1]+0.6708203932499369*fl[1]-(0.5*fr[0]-0.5*fl[0])))*sgn(0.6324555320336759*alphaDrSurf[2]+0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0])+fr[8]+fl[8]+1.5*(fr[7]+fl[7])-0.7745966692414833*fr[6]+0.7745966692414833*fl[6]+1.118033988749895*(fr[5]+fl[5])+0.4472135954999579*(fr[4]+fl[4])-1.161895003862225*fr[3]+1.161895003862225*fl[3]-0.8660254037844386*fr[2]+0.8660254037844386*fl[2]+0.6708203932499369*(fr[1]+fl[1])+0.5*(fr[0]+fl[0])); 

  double fUpwind[3];
  fUpwind[0] = 0.07856742013183861*(5.0*fUpwindQuad[2]+8.0*fUpwindQuad[1]+5.0*fUpwindQuad[0]); 
  fUpwind[1] = 0.5270462766947298*(fUpwindQuad[2]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.3513641844631533*(fUpwindQuad[2]-2.0*fUpwindQuad[1]+fUpwindQuad[0]); 

  double Gdiff[9]; 
  double Ghat[9]; 
  double incr2[9]; 


  incr2[2] = 0.0078125*(38.34057902536163*nuVtSqSum[2]*fr[8]+38.34057902536163*nuVtSqSum[2]*fl[8]+38.34057902536163*nuVtSqSum[1]*fr[7]+38.34057902536163*nuVtSqSum[1]*fl[7]-55.15432893255071*nuVtSqSum[2]*fr[6]+55.15432893255071*nuVtSqSum[2]*fl[6]+38.34057902536163*nuVtSqSum[0]*fr[5]+38.34057902536163*nuVtSqSum[0]*fl[5]+39.19183588453087*nuVtSqSum[2]*fr[4]+39.19183588453087*nuVtSqSum[2]*fl[4]-55.15432893255071*nuVtSqSum[1]*fr[3]+55.15432893255071*nuVtSqSum[1]*fl[3]-55.15432893255071*nuVtSqSum[0]*fr[2]+55.15432893255071*nuVtSqSum[0]*fl[2]+(39.19183588453087*fr[1]+39.19183588453087*fl[1])*nuVtSqSum[1]+(39.19183588453087*fr[0]+39.19183588453087*fl[0])*nuVtSqSum[0]); 
  incr2[3] = 0.0015625*(171.4642819948225*nuVtSqSum[1]*fr[8]+171.4642819948225*nuVtSqSum[1]*fl[8]+(171.4642819948225*nuVtSqSum[2]+191.7028951268082*nuVtSqSum[0])*fr[7]+(171.4642819948225*nuVtSqSum[2]+191.7028951268082*nuVtSqSum[0])*fl[7]-246.6576574931337*nuVtSqSum[1]*fr[6]+246.6576574931337*nuVtSqSum[1]*fl[6]+191.7028951268082*nuVtSqSum[1]*fr[5]+191.7028951268082*nuVtSqSum[1]*fl[5]+175.2712184016533*nuVtSqSum[1]*fr[4]+175.2712184016533*nuVtSqSum[1]*fl[4]+((-246.6576574931337*nuVtSqSum[2])-275.7716446627535*nuVtSqSum[0])*fr[3]+(246.6576574931337*nuVtSqSum[2]+275.7716446627535*nuVtSqSum[0])*fl[3]+(175.2712184016533*fr[1]+175.2712184016533*fl[1])*nuVtSqSum[2]-275.7716446627535*nuVtSqSum[1]*fr[2]+275.7716446627535*nuVtSqSum[1]*fl[2]+(195.9591794226543*fr[0]+195.9591794226543*fl[0])*nuVtSqSum[1]+195.9591794226543*nuVtSqSum[0]*fr[1]+195.9591794226543*nuVtSqSum[0]*fl[1]); 
  incr2[5] = -0.0078125*(148.492424049175*nuVtSqSum[2]*fr[8]+148.492424049175*nuVtSqSum[2]*fl[8]+148.492424049175*nuVtSqSum[1]*fr[7]+148.492424049175*nuVtSqSum[1]*fl[7]-213.6117974270148*nuVtSqSum[2]*fr[6]+213.6117974270148*nuVtSqSum[2]*fl[6]+148.492424049175*nuVtSqSum[0]*fr[5]+148.492424049175*nuVtSqSum[0]*fl[5]+151.7893276880823*nuVtSqSum[2]*fr[4]+151.7893276880823*nuVtSqSum[2]*fl[4]-213.6117974270148*nuVtSqSum[1]*fr[3]+213.6117974270148*nuVtSqSum[1]*fl[3]-213.6117974270148*nuVtSqSum[0]*fr[2]+213.6117974270148*nuVtSqSum[0]*fl[2]+(151.7893276880823*fr[1]+151.7893276880823*fl[1])*nuVtSqSum[1]+(151.7893276880823*fr[0]+151.7893276880823*fl[0])*nuVtSqSum[0]); 
  incr2[6] = 2.232142857142857e-4*((857.3214099741125*nuVtSqSum[2]+1341.920265887657*nuVtSqSum[0])*fr[8]+(857.3214099741125*nuVtSqSum[2]+1341.920265887657*nuVtSqSum[0])*fl[8]+1200.249973963757*nuVtSqSum[1]*fr[7]+1200.249973963757*nuVtSqSum[1]*fl[7]+((-1233.288287465668*nuVtSqSum[2])-1930.401512639275*nuVtSqSum[0])*fr[6]+(1233.288287465668*nuVtSqSum[2]+1930.401512639275*nuVtSqSum[0])*fl[6]+1341.920265887657*nuVtSqSum[2]*fr[5]+1341.920265887657*nuVtSqSum[2]*fl[5]+(876.3560920082665*nuVtSqSum[2]+1371.71425595858*nuVtSqSum[0])*fr[4]+(876.3560920082665*nuVtSqSum[2]+1371.71425595858*nuVtSqSum[0])*fl[4]-1726.603602451936*nuVtSqSum[1]*fr[3]+1726.603602451936*nuVtSqSum[1]*fl[3]+((-1930.401512639275*fr[2])+1930.401512639275*fl[2]+1371.71425595858*fr[0]+1371.71425595858*fl[0])*nuVtSqSum[2]+(1226.898528811573*fr[1]+1226.898528811573*fl[1])*nuVtSqSum[1]); 
  incr2[7] = -0.0078125*(132.815661727072*nuVtSqSum[1]*fr[8]+132.815661727072*nuVtSqSum[1]*fl[8]+(132.815661727072*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*fr[7]+(132.815661727072*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*fl[7]-191.0601999370879*nuVtSqSum[1]*fr[6]+191.0601999370879*nuVtSqSum[1]*fl[6]+148.492424049175*nuVtSqSum[1]*fr[5]+148.492424049175*nuVtSqSum[1]*fl[5]+135.7645019878172*nuVtSqSum[1]*fr[4]+135.7645019878172*nuVtSqSum[1]*fl[4]+((-191.0601999370879*nuVtSqSum[2])-213.6117974270148*nuVtSqSum[0])*fr[3]+(191.0601999370879*nuVtSqSum[2]+213.6117974270148*nuVtSqSum[0])*fl[3]+(135.7645019878172*fr[1]+135.7645019878172*fl[1])*nuVtSqSum[2]-213.6117974270148*nuVtSqSum[1]*fr[2]+213.6117974270148*nuVtSqSum[1]*fl[2]+(151.7893276880823*fr[0]+151.7893276880823*fl[0])*nuVtSqSum[1]+151.7893276880823*nuVtSqSum[0]*fr[1]+151.7893276880823*nuVtSqSum[0]*fl[1]); 
  incr2[8] = -0.001116071428571429*((664.0783086353599*nuVtSqSum[2]+1039.446968344225*nuVtSqSum[0])*fr[8]+(664.0783086353599*nuVtSqSum[2]+1039.446968344225*nuVtSqSum[0])*fl[8]+929.7096320895039*nuVtSqSum[1]*fr[7]+929.7096320895039*nuVtSqSum[1]*fl[7]+((-955.3009996854395*nuVtSqSum[2])-1495.282581989103*nuVtSqSum[0])*fr[6]+(955.3009996854395*nuVtSqSum[2]+1495.282581989103*nuVtSqSum[0])*fl[6]+1039.446968344225*nuVtSqSum[2]*fr[5]+1039.446968344225*nuVtSqSum[2]*fl[5]+(678.8225099390861*nuVtSqSum[2]+1062.525293816576*nuVtSqSum[0])*fr[4]+(678.8225099390861*nuVtSqSum[2]+1062.525293816576*nuVtSqSum[0])*fl[4]-1337.421399559615*nuVtSqSum[1]*fr[3]+1337.421399559615*nuVtSqSum[1]*fl[3]+((-1495.282581989103*fr[2])+1495.282581989103*fl[2]+1062.525293816576*fr[0]+1062.525293816576*fl[0])*nuVtSqSum[2]+(950.3515139147205*fr[1]+950.3515139147205*fl[1])*nuVtSqSum[1]); 


  Gdiff[0] = -0.0125*(3.16227766016838*((-24.0*nuVtSqSum[2]*fr[8])+24.0*nuVtSqSum[2]*fl[8]-24.0*nuVtSqSum[0]*fr[5]+24.0*nuVtSqSum[0]*fl[5])+3.872983346207417*(2.449489742783178*(8.0*nuVtSqSum[1]*fl[7]-8.0*nuVtSqSum[1]*fr[7])+3.16227766016838*(11.0*nuVtSqSum[2]*fr[6]+11.0*nuVtSqSum[2]*fl[6]))+1.414213562373095*((-75.0*nuVtSqSum[2]*fr[4])+75.0*nuVtSqSum[2]*fl[4]+(75.0*fl[1]-75.0*fr[1])*nuVtSqSum[1]+(75.0*fl[0]-75.0*fr[0])*nuVtSqSum[0])+2.449489742783178*(55.0*nuVtSqSum[1]*fr[3]+55.0*nuVtSqSum[1]*fl[3]+55.0*nuVtSqSum[0]*fr[2]+55.0*nuVtSqSum[0]*fl[2])); 
  Gdiff[1] = 0.0025*(1.414213562373095*(240.0*nuVtSqSum[1]*fr[8]-240.0*nuVtSqSum[1]*fl[8]+(375.0*fr[0]-375.0*fl[0])*nuVtSqSum[1]+375.0*nuVtSqSum[0]*fr[1]-375.0*nuVtSqSum[0]*fl[1])+3.872983346207417*(5.477225575051662*(16.0*nuVtSqSum[2]*fr[7]-16.0*nuVtSqSum[2]*fl[7])+2.449489742783178*(40.0*nuVtSqSum[0]*fr[7]-40.0*nuVtSqSum[0]*fl[7])+1.414213562373095*((-110.0*nuVtSqSum[1]*fr[6])-110.0*nuVtSqSum[1]*fl[6]))+2.23606797749979*(1.414213562373095*(120.0*nuVtSqSum[1]*fr[5]-120.0*nuVtSqSum[1]*fl[5]+150.0*nuVtSqSum[1]*fr[4]-150.0*nuVtSqSum[1]*fl[4]+(150.0*fr[1]-150.0*fl[1])*nuVtSqSum[2])+2.449489742783178*((-110.0*nuVtSqSum[2]*fr[3])-110.0*nuVtSqSum[2]*fl[3]))+2.449489742783178*((-275.0*nuVtSqSum[0]*fr[3])-275.0*nuVtSqSum[0]*fl[3]-275.0*nuVtSqSum[1]*fr[2]-275.0*nuVtSqSum[1]*fl[2])); 
  Gdiff[4] = 3.571428571428572e-4*(2.23606797749979*(1.414213562373095*(840.0*nuVtSqSum[0]*fr[8]-840.0*nuVtSqSum[0]*fl[8]+840.0*nuVtSqSum[2]*fr[5]-840.0*nuVtSqSum[2]*fl[5]+750.0*nuVtSqSum[2]*fr[4]-750.0*nuVtSqSum[2]*fl[4]+(1050.0*fr[1]-1050.0*fl[1])*nuVtSqSum[1])+2.449489742783178*((-770.0*nuVtSqSum[1]*fr[3])-770.0*nuVtSqSum[1]*fl[3]))+1.414213562373095*(1200.0*nuVtSqSum[2]*fr[8]-1200.0*nuVtSqSum[2]*fl[8]+2625.0*nuVtSqSum[0]*fr[4]-2625.0*nuVtSqSum[0]*fl[4]+(2625.0*fr[0]-2625.0*fl[0])*nuVtSqSum[2])+3.872983346207417*(2.23606797749979*(2.449489742783178*(112.0*nuVtSqSum[1]*fr[7]-112.0*nuVtSqSum[1]*fl[7])+1.414213562373095*((-385.0*nuVtSqSum[0]*fr[6])-385.0*nuVtSqSum[0]*fl[6]))+1.414213562373095*((-550.0*nuVtSqSum[2]*fr[6])-550.0*nuVtSqSum[2]*fl[6]))+2.449489742783178*((-1925.0*fr[2])-1925.0*fl[2])*nuVtSqSum[2]); 

  Ghat[0] = Gdiff[0]*rdv2L+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+0.8944271909999159*alphaDrSurf[1]*fUpwind[2]+0.8944271909999159*fUpwind[1]*alphaDrSurf[2]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]; 
  Ghat[4] = Gdiff[4]*rdv2L+0.6388765649999399*alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]+0.8944271909999159*alphaDrSurf[1]*fUpwind[1]; 

  double incr1[9]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -1.118033988749895*Ghat[0]; 
  incr1[6] = 0.8660254037844387*Ghat[4]; 
  incr1[7] = -1.118033988749895*Ghat[1]; 
  incr1[8] = -1.118033988749895*Ghat[4]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr2[8]*rdvSq4R+incr1[8]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr2[5]*rdvSq4L-1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr2[7]*rdvSq4L-1.0*incr1[7]*rdv2L; 
  outl[8] += incr2[8]*rdvSq4L-1.0*incr1[8]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUx[2])/nuSum-(0.7071067811865475*sumNuUx[0])/nuSum+wl[1]); 
} 
