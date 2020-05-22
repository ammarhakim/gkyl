#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x1vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:          Cell-center coordinates. 
  // dxv[2]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[4]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[10]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 

  double fjump[10]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(-1*fr[9])); 

  double alphaDrag[4]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 
  alphaDrag[3] = -1.0*sumNuUx[3]; 

  double Ghat[10]; 
  for(unsigned short int i=0; i<10; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.929292090055312*nuVtSqSum[0]*fr[9])-1.929292090055312*nuVtSqSum[0]*fl[9]+3.866990209613929*nuVtSqSum[3]*fr[8]-3.866990209613929*nuVtSqSum[3]*fl[8]+4.39754237117165*nuVtSqSum[1]*fr[7]-4.39754237117165*nuVtSqSum[1]*fl[7]-5.473078644031159*nuVtSqSum[2]*fr[6]-5.473078644031159*nuVtSqSum[2]*fl[6]+4.397542371171649*nuVtSqSum[0]*fr[5]-4.397542371171649*nuVtSqSum[0]*fl[5]+3.866990209613929*nuVtSqSum[2]*fr[4]-3.866990209613929*nuVtSqSum[2]*fl[4]-5.473078644031159*nuVtSqSum[1]*fr[3]-5.473078644031159*nuVtSqSum[1]*fl[3]-5.473078644031159*nuVtSqSum[0]*fr[2]-5.473078644031159*nuVtSqSum[0]*fl[2]+3.866990209613929*fr[1]*nuVtSqSum[1]-3.866990209613929*fl[1]*nuVtSqSum[1]+3.866990209613929*fr[0]*nuVtSqSum[0]-3.866990209613929*fl[0]*nuVtSqSum[0])*rdv-1.322875655532295*fjump[9]+alphaDrag[0]*(0.9354143466934851*favg[9]+0.7905694150420947*favg[5]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+0.3535533905932737*alphaDrag[3]*favg[8]+alphaDrag[1]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])-1.118033988749895*fjump[5]-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = ((-1.929292090055312*nuVtSqSum[1]*fr[9])-1.929292090055312*nuVtSqSum[1]*fl[9]+3.396416424888148*nuVtSqSum[2]*fr[8]-3.396416424888148*nuVtSqSum[2]*fl[8]+3.933281470350169*nuVtSqSum[2]*fr[7]+4.39754237117165*nuVtSqSum[0]*fr[7]-3.933281470350169*nuVtSqSum[2]*fl[7]-4.39754237117165*nuVtSqSum[0]*fl[7]-4.807060063166761*nuVtSqSum[3]*fr[6]-4.89527035770242*nuVtSqSum[1]*fr[6]-4.807060063166761*nuVtSqSum[3]*fl[6]-4.89527035770242*nuVtSqSum[1]*fl[6]+4.397542371171649*nuVtSqSum[1]*fr[5]-4.397542371171649*nuVtSqSum[1]*fl[5]+3.396416424888148*nuVtSqSum[3]*fr[4]+3.458741190809164*nuVtSqSum[1]*fr[4]-3.396416424888148*nuVtSqSum[3]*fl[4]-3.458741190809164*nuVtSqSum[1]*fl[4]-4.89527035770242*nuVtSqSum[2]*fr[3]-5.473078644031159*nuVtSqSum[0]*fr[3]-4.89527035770242*nuVtSqSum[2]*fl[3]-5.473078644031159*nuVtSqSum[0]*fl[3]+3.458741190809164*fr[1]*nuVtSqSum[2]-3.458741190809164*fl[1]*nuVtSqSum[2]-5.473078644031159*nuVtSqSum[1]*fr[2]-5.473078644031159*nuVtSqSum[1]*fl[2]+3.866990209613929*fr[0]*nuVtSqSum[1]-3.866990209613929*fl[0]*nuVtSqSum[1]+3.866990209613929*nuVtSqSum[0]*fr[1]-3.866990209613929*nuVtSqSum[0]*fl[1])*rdv+alphaDrag[1]*(0.9354143466934851*favg[9]+0.5477225575051661*favg[6]+0.7905694150420947*favg[5]+0.3162277660168379*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+alphaDrag[2]*(0.3105295017040592*favg[8]+0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])-1.118033988749895*fjump[7]+alphaDrag[0]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+alphaDrag[3]*(0.5378528742004769*favg[6]+0.3105295017040592*favg[4])-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = ((-1.929292090055312*nuVtSqSum[2]*fr[9])-1.929292090055312*nuVtSqSum[2]*fl[9]+2.305827460539443*nuVtSqSum[3]*fr[8]+3.396416424888148*nuVtSqSum[1]*fr[8]-2.305827460539443*nuVtSqSum[3]*fl[8]-3.396416424888148*nuVtSqSum[1]*fl[8]+3.86240572873861*nuVtSqSum[3]*fr[7]+3.933281470350169*nuVtSqSum[1]*fr[7]-3.86240572873861*nuVtSqSum[3]*fl[7]-3.933281470350169*nuVtSqSum[1]*fl[7]-3.496621684073157*nuVtSqSum[2]*fr[6]-5.473078644031159*nuVtSqSum[0]*fr[6]-3.496621684073157*nuVtSqSum[2]*fl[6]-5.473078644031159*nuVtSqSum[0]*fl[6]+4.397542371171649*nuVtSqSum[2]*fr[5]-4.397542371171649*nuVtSqSum[2]*fl[5]+2.470529422006546*nuVtSqSum[2]*fr[4]+3.866990209613929*nuVtSqSum[0]*fr[4]-2.470529422006546*nuVtSqSum[2]*fl[4]-3.866990209613929*nuVtSqSum[0]*fl[4]-4.807060063166761*fr[3]*nuVtSqSum[3]-4.807060063166761*fl[3]*nuVtSqSum[3]+3.396416424888148*fr[1]*nuVtSqSum[3]-3.396416424888148*fl[1]*nuVtSqSum[3]-4.89527035770242*nuVtSqSum[1]*fr[3]-4.89527035770242*nuVtSqSum[1]*fl[3]-5.473078644031159*fr[2]*nuVtSqSum[2]-5.473078644031159*fl[2]*nuVtSqSum[2]+3.866990209613929*fr[0]*nuVtSqSum[2]-3.866990209613929*fl[0]*nuVtSqSum[2]+3.458741190809164*fr[1]*nuVtSqSum[1]-3.458741190809164*fl[1]*nuVtSqSum[1])*rdv+alphaDrag[2]*(0.9354143466934851*favg[9]+0.3912303982179757*favg[6]+0.7905694150420947*favg[5]+0.2258769757263128*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+alphaDrag[1]*(0.3105295017040592*favg[8]+0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])+alphaDrag[3]*(0.210818510677892*favg[8]+0.6943650748294133*favg[7]+0.537852874200477*favg[3]+0.3105295017040592*favg[1])-0.8660254037844387*fjump[6]+alphaDrag[0]*(0.6123724356957944*favg[6]+0.3535533905932737*favg[4])-0.5*fjump[4]; 
  Ghat[8] = ((-1.929292090055312*nuVtSqSum[3]*fr[9])-1.929292090055312*nuVtSqSum[3]*fl[9]+2.305827460539443*nuVtSqSum[2]*fr[8]+3.866990209613929*nuVtSqSum[0]*fr[8]-2.305827460539443*nuVtSqSum[2]*fl[8]-3.866990209613929*nuVtSqSum[0]*fl[8]+3.86240572873861*nuVtSqSum[2]*fr[7]-3.86240572873861*nuVtSqSum[2]*fl[7]-3.263513571801613*nuVtSqSum[3]*fr[6]-4.807060063166761*nuVtSqSum[1]*fr[6]-3.263513571801613*nuVtSqSum[3]*fl[6]-4.807060063166761*nuVtSqSum[1]*fl[6]+4.397542371171649*nuVtSqSum[3]*fr[5]-4.397542371171649*nuVtSqSum[3]*fl[5]+2.305827460539443*nuVtSqSum[3]*fr[4]+3.396416424888148*nuVtSqSum[1]*fr[4]-2.305827460539443*nuVtSqSum[3]*fl[4]-3.396416424888148*nuVtSqSum[1]*fl[4]-5.473078644031159*fr[2]*nuVtSqSum[3]-5.473078644031159*fl[2]*nuVtSqSum[3]+3.866990209613929*fr[0]*nuVtSqSum[3]-3.866990209613929*fl[0]*nuVtSqSum[3]-4.807060063166761*nuVtSqSum[2]*fr[3]-4.807060063166761*nuVtSqSum[2]*fl[3]+3.396416424888148*fr[1]*nuVtSqSum[2]-3.396416424888148*fl[1]*nuVtSqSum[2])*rdv+alphaDrag[3]*(0.9354143466934851*favg[9]+0.3651483716701107*favg[6]+0.7905694150420947*favg[5]+0.210818510677892*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.5*fjump[8]+0.3535533905932737*alphaDrag[0]*favg[8]+alphaDrag[2]*(0.210818510677892*favg[8]+0.6943650748294133*favg[7]+0.537852874200477*favg[3]+0.3105295017040592*favg[1])+alphaDrag[1]*(0.5378528742004769*favg[6]+0.3105295017040592*favg[4]); 

  double incr1[10]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -1.118033988749895*Ghat[0]; 
  incr1[6] = 0.8660254037844387*Ghat[4]; 
  incr1[7] = -1.118033988749895*Ghat[1]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = 1.322875655532295*Ghat[0]; 

  double incr2[10]; 
  incr2[2] = nuVtSqSum[0]*((-0.2314550249431378*fr[9])+0.2314550249431378*fl[9]+0.4279082480509108*(fr[5]+fl[5])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+0.3061862178478971*nuVtSqSum[3]*(fr[8]+fl[8])+nuVtSqSum[1]*(0.4279082480509107*(fr[7]+fl[7])-0.4640388251536715*fr[3]+0.4640388251536715*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.4640388251536715*fr[6])+0.4640388251536715*fl[6]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[3] = nuVtSqSum[1]*((-0.2314550249431378*fr[9])+0.2314550249431378*fl[9]-0.4150489428970995*fr[6]+0.4150489428970995*fl[6]+0.4279082480509108*(fr[5]+fl[5])+0.273861278752583*(fr[4]+fl[4])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[2]*(0.2689264371002384*(fr[8]+fl[8])+0.3827327723098714*(fr[7]+fl[7])-0.4150489428970995*fr[3]+0.4150489428970995*fl[3]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[0]*(0.4279082480509107*(fr[7]+fl[7])-0.4640388251536715*fr[3]+0.4640388251536715*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[3]*((-0.4075699709865777*fr[6])+0.4075699709865777*fl[6]+0.2689264371002384*(fr[4]+fl[4])); 
  incr2[5] = nuVtSqSum[0]*(0.896421457000795*fr[9]-0.896421457000795*fl[9]-1.657281518405969*(fr[5]+fl[5])+1.797214641813825*fr[2]-1.797214641813825*fl[2]-1.185854122563142*(fr[0]+fl[0]))-1.185854122563142*nuVtSqSum[3]*(fr[8]+fl[8])+nuVtSqSum[1]*((-1.65728151840597*(fr[7]+fl[7]))+1.797214641813825*fr[3]-1.797214641813825*fl[3]-1.185854122563142*(fr[1]+fl[1]))+nuVtSqSum[2]*(1.797214641813825*fr[6]-1.797214641813825*fl[6]-1.185854122563142*(fr[4]+fl[4])); 
  incr2[6] = nuVtSqSum[2]*((-0.2314550249431378*fr[9])+0.2314550249431378*fl[9]-0.2964635306407854*fr[6]+0.2964635306407854*fl[6]+0.4279082480509108*(fr[5]+fl[5])+0.1956151991089878*(fr[4]+fl[4])-0.4640388251536715*fr[2]+0.4640388251536715*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[1]*(0.2689264371002384*(fr[8]+fl[8])+0.3827327723098713*(fr[7]+fl[7])-0.4150489428970995*fr[3]+0.4150489428970995*fl[3]+0.273861278752583*(fr[1]+fl[1]))+nuVtSqSum[3]*(0.1825741858350553*(fr[8]+fl[8])+0.3758361214393465*(fr[7]+fl[7])-0.4075699709865777*fr[3]+0.4075699709865777*fl[3]+0.2689264371002384*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.4640388251536715*fr[6])+0.4640388251536715*fl[6]+0.3061862178478971*(fr[4]+fl[4])); 
  incr2[7] = nuVtSqSum[1]*(0.8964214570007949*fr[9]-0.8964214570007949*fl[9]+1.60747764370146*fr[6]-1.60747764370146*fl[6]-1.65728151840597*(fr[5]+fl[5])-1.060660171779821*(fr[4]+fl[4])+1.797214641813825*fr[2]-1.797214641813825*fl[2]-1.185854122563142*(fr[0]+fl[0]))+nuVtSqSum[2]*((-1.04154761224412*(fr[8]+fl[8]))-1.482317653203927*(fr[7]+fl[7])+1.60747764370146*fr[3]-1.60747764370146*fl[3]-1.060660171779821*(fr[1]+fl[1]))+nuVtSqSum[0]*((-1.657281518405969*(fr[7]+fl[7]))+1.797214641813825*fr[3]-1.797214641813825*fl[3]-1.185854122563142*(fr[1]+fl[1]))+nuVtSqSum[3]*(1.578511710045255*fr[6]-1.578511710045255*fl[6]-1.04154761224412*(fr[4]+fl[4])); 
  incr2[9] = nuVtSqSum[0]*((-2.121320343559642*fr[9])+2.121320343559642*fl[9]+3.921843874378478*(fr[5]+fl[5])-4.252986083330156*fr[2]+4.252986083330156*fl[2]+2.806243040080455*(fr[0]+fl[0]))+2.806243040080455*nuVtSqSum[3]*(fr[8]+fl[8])+nuVtSqSum[1]*(3.921843874378477*(fr[7]+fl[7])-4.252986083330156*fr[3]+4.252986083330156*fl[3]+2.806243040080455*(fr[1]+fl[1]))+nuVtSqSum[2]*((-4.252986083330157*fr[6])+4.252986083330157*fl[6]+2.806243040080455*(fr[4]+fl[4])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr2[5]*rdvSq4L-1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr2[7]*rdvSq4L-1.0*incr1[7]*rdv2L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += incr1[9]*rdv2L-1.0*incr2[9]*rdvSq4L; 

  return std::abs((0.7905694150420947*sumNuUx[2])/nuSum-(0.7071067811865475*sumNuUx[0])/nuSum+wl[1]); 
} 
