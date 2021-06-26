#include <VmLBOModDecl.h> 
double VmLBOconstNuSurfNonUniform1x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[9]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[15]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 

  double fjump[15]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(-1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(-1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(-1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 

  double Ghat[15]; 
  for(unsigned short int i=0; i<15; ++i){ 
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

  Ghat[0] = ((265.631323454144*nuVtSqSum[0]*dxvl1R3*dxvr1R3+151.7893276880823*nuVtSqSum[0]*dxvl1R4*dxvr1R2-66.40783086353598*nuVtSqSum[0]*dxvl1R5*dxvr[1]-47.43416490252571*nuVtSqSum[0]*dxvl1R6)*fr[12]+(47.43416490252571*nuVtSqSum[0]*dxvr1R6+66.40783086353598*nuVtSqSum[0]*dxvl[1]*dxvr1R5-151.7893276880823*nuVtSqSum[0]*dxvl1R2*dxvr1R4-265.631323454144*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[12]+424.2640687119287*dxvl1R3*dxvr1R3*nuVtSqSum[2]*fr[11]-424.2640687119287*dxvl1R3*dxvr1R3*nuVtSqSum[2]*fl[11]+((-514.3928459844674*dxvl1R3*dxvr1R3)-97.97958971132716*dxvl1R4*dxvr1R2+61.23724356957945*dxvl1R5*dxvr[1]+12.24744871391589*dxvl1R6)*nuVtSqSum[1]*fr[5]+(12.24744871391589*dxvr1R6+61.23724356957945*dxvl[1]*dxvr1R5-97.97958971132716*dxvl1R2*dxvr1R4-514.3928459844674*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[5]+((-514.3928459844674*nuVtSqSum[0]*dxvl1R3*dxvr1R3)-97.97958971132716*nuVtSqSum[0]*dxvl1R4*dxvr1R2+61.23724356957945*nuVtSqSum[0]*dxvl1R5*dxvr[1]+12.24744871391589*nuVtSqSum[0]*dxvl1R6)*fr[2]+(12.24744871391589*nuVtSqSum[0]*dxvr1R6+61.23724356957945*nuVtSqSum[0]*dxvl[1]*dxvr1R5-97.97958971132716*nuVtSqSum[0]*dxvl1R2*dxvr1R4-514.3928459844674*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[2]+(424.2640687119287*dxvl1R3*dxvr1R3*fr[1]-424.2640687119287*dxvl1R3*dxvr1R3*fl[1])*nuVtSqSum[1]+(424.2640687119287*fr[0]-424.2640687119287*fl[0])*nuVtSqSum[0]*dxvl1R3*dxvr1R3)/(5.0*dxvl[1]*dxvr1R6+25.0*dxvl1R2*dxvr1R5+50.0*dxvl1R3*dxvr1R4+50.0*dxvl1R4*dxvr1R3+25.0*dxvl1R5*dxvr1R2+5.0*dxvl1R6*dxvr[1])-0.5*(2.23606797749979*fjump[12]+1.732050807568877*fjump[2]+fjump[0])+0.25*(3.16227766016838*alphaDrag[0]*favg[12]+1.414213562373095*alphaDrag[2]*favg[11]+2.449489742783178*alphaDrag[1]*favg[5]+2.449489742783178*alphaDrag[0]*favg[2]+1.414213562373095*alphaDrag[1]*favg[1]+1.414213562373095*alphaDrag[0]*favg[0]); 
  Ghat[1] = ((265.631323454144*dxvl1R3*dxvr1R3+151.7893276880823*dxvl1R4*dxvr1R2-66.40783086353598*dxvl1R5*dxvr[1]-47.43416490252571*dxvl1R6)*nuVtSqSum[1]*fr[12]+(47.43416490252571*dxvr1R6+66.40783086353598*dxvl[1]*dxvr1R5-151.7893276880823*dxvl1R2*dxvr1R4-265.631323454144*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[12]+379.4733192202058*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fr[11]-379.4733192202058*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fl[11]+(((-460.0869483043397*dxvl1R3*dxvr1R3)-87.63560920082662*dxvl1R4*dxvr1R2+54.77225575051663*dxvl1R5*dxvr[1]+10.95445115010332*dxvl1R6)*nuVtSqSum[2]-514.3928459844674*nuVtSqSum[0]*dxvl1R3*dxvr1R3-97.97958971132716*nuVtSqSum[0]*dxvl1R4*dxvr1R2+61.23724356957945*nuVtSqSum[0]*dxvl1R5*dxvr[1]+12.24744871391589*nuVtSqSum[0]*dxvl1R6)*fr[5]+((10.95445115010332*dxvr1R6+54.77225575051663*dxvl[1]*dxvr1R5-87.63560920082662*dxvl1R2*dxvr1R4-460.0869483043397*dxvl1R3*dxvr1R3)*nuVtSqSum[2]+12.24744871391589*nuVtSqSum[0]*dxvr1R6+61.23724356957945*nuVtSqSum[0]*dxvl[1]*dxvr1R5-97.97958971132716*nuVtSqSum[0]*dxvl1R2*dxvr1R4-514.3928459844674*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[5]+(379.4733192202058*dxvl1R3*dxvr1R3*fr[1]-379.4733192202058*dxvl1R3*dxvr1R3*fl[1])*nuVtSqSum[2]+((-514.3928459844674*dxvl1R3*dxvr1R3)-97.97958971132716*dxvl1R4*dxvr1R2+61.23724356957945*dxvl1R5*dxvr[1]+12.24744871391589*dxvl1R6)*nuVtSqSum[1]*fr[2]+(12.24744871391589*dxvr1R6+61.23724356957945*dxvl[1]*dxvr1R5-97.97958971132716*dxvl1R2*dxvr1R4-514.3928459844674*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[2]+(424.2640687119287*fr[0]-424.2640687119287*fl[0])*dxvl1R3*dxvr1R3*nuVtSqSum[1]+424.2640687119287*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fr[1]-424.2640687119287*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fl[1])/(5.0*dxvl[1]*dxvr1R6+25.0*dxvl1R2*dxvr1R5+50.0*dxvl1R3*dxvr1R4+50.0*dxvl1R4*dxvr1R3+25.0*dxvl1R5*dxvr1R2+5.0*dxvl1R6*dxvr[1])+0.05*(15.8113883008419*alphaDrag[1]*favg[12]+6.324555320336761*alphaDrag[1]*favg[11]+(10.95445115010332*alphaDrag[2]+12.24744871391589*alphaDrag[0])*favg[5]+12.24744871391589*alphaDrag[1]*favg[2]+6.324555320336761*favg[1]*alphaDrag[2]+7.071067811865476*alphaDrag[0]*favg[1]+7.071067811865476*favg[0]*alphaDrag[1])-0.5*(1.732050807568877*fjump[5]+fjump[1]); 
  Ghat[3] = (-(1.0*((102.8785691968935*nuVtSqSum[0]*dxvl1R3*dxvr1R3+19.59591794226543*nuVtSqSum[0]*dxvl1R4*dxvr1R2-12.24744871391589*nuVtSqSum[0]*dxvl1R5*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R6)*fr[7]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R6)-12.24744871391589*nuVtSqSum[0]*dxvl[1]*dxvr1R5+19.59591794226543*nuVtSqSum[0]*dxvl1R2*dxvr1R4+102.8785691968935*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[7]-84.85281374238573*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fr[6]+84.85281374238573*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fl[6]-84.85281374238573*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fr[3]+84.85281374238573*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fl[3]))/(dxvl[1]*dxvr1R6+5.0*dxvl1R2*dxvr1R5+10.0*dxvl1R3*dxvr1R4+10.0*dxvl1R4*dxvr1R3+5.0*dxvl1R5*dxvr1R2+dxvl1R6*dxvr[1]))-0.5*(1.732050807568877*fjump[7]+fjump[3])+0.25*(2.449489742783178*alphaDrag[0]*favg[7]+1.414213562373095*alphaDrag[1]*favg[6]+1.414213562373095*alphaDrag[0]*favg[3]); 
  Ghat[4] = (-(1.0*((102.8785691968935*nuVtSqSum[0]*dxvl1R3*dxvr1R3+19.59591794226543*nuVtSqSum[0]*dxvl1R4*dxvr1R2-12.24744871391589*nuVtSqSum[0]*dxvl1R5*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R6)*fr[9]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R6)-12.24744871391589*nuVtSqSum[0]*dxvl[1]*dxvr1R5+19.59591794226543*nuVtSqSum[0]*dxvl1R2*dxvr1R4+102.8785691968935*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[9]-84.85281374238573*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fr[8]+84.85281374238573*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fl[8]-84.85281374238573*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fr[4]+84.85281374238573*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fl[4]))/(dxvl[1]*dxvr1R6+5.0*dxvl1R2*dxvr1R5+10.0*dxvl1R3*dxvr1R4+10.0*dxvl1R4*dxvr1R3+5.0*dxvl1R5*dxvr1R2+dxvl1R6*dxvr[1]))-0.5*(1.732050807568877*fjump[9]+fjump[4])+0.25*(2.449489742783178*alphaDrag[0]*favg[9]+1.414213562373095*alphaDrag[1]*favg[8]+1.414213562373095*alphaDrag[0]*favg[4]); 
  Ghat[6] = (-(1.0*((102.8785691968935*dxvl1R3*dxvr1R3+19.59591794226543*dxvl1R4*dxvr1R2-12.24744871391589*dxvl1R5*dxvr[1]-2.449489742783178*dxvl1R6)*nuVtSqSum[1]*fr[7]+((-2.449489742783178*dxvr1R6)-12.24744871391589*dxvl[1]*dxvr1R5+19.59591794226543*dxvl1R2*dxvr1R4+102.8785691968935*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[7]+((-75.89466384404115*dxvl1R3*dxvr1R3*nuVtSqSum[2])-84.85281374238573*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fr[6]+(75.89466384404115*dxvl1R3*dxvr1R3*nuVtSqSum[2]+84.85281374238573*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[6]-84.85281374238573*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fr[3]+84.85281374238573*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fl[3]))/(dxvl[1]*dxvr1R6+5.0*dxvl1R2*dxvr1R5+10.0*dxvl1R3*dxvr1R4+10.0*dxvl1R4*dxvr1R3+5.0*dxvl1R5*dxvr1R2+dxvl1R6*dxvr[1]))+0.05*(12.24744871391589*alphaDrag[1]*favg[7]+(6.324555320336761*alphaDrag[2]+7.071067811865476*alphaDrag[0])*favg[6]+7.071067811865476*alphaDrag[1]*favg[3])-0.5*fjump[6]; 
  Ghat[8] = (-(1.0*((102.8785691968935*dxvl1R3*dxvr1R3+19.59591794226543*dxvl1R4*dxvr1R2-12.24744871391589*dxvl1R5*dxvr[1]-2.449489742783178*dxvl1R6)*nuVtSqSum[1]*fr[9]+((-2.449489742783178*dxvr1R6)-12.24744871391589*dxvl[1]*dxvr1R5+19.59591794226543*dxvl1R2*dxvr1R4+102.8785691968935*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[9]+((-75.89466384404115*dxvl1R3*dxvr1R3*nuVtSqSum[2])-84.85281374238573*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fr[8]+(75.89466384404115*dxvl1R3*dxvr1R3*nuVtSqSum[2]+84.85281374238573*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[8]-84.85281374238573*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fr[4]+84.85281374238573*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fl[4]))/(dxvl[1]*dxvr1R6+5.0*dxvl1R2*dxvr1R5+10.0*dxvl1R3*dxvr1R4+10.0*dxvl1R4*dxvr1R3+5.0*dxvl1R5*dxvr1R2+dxvl1R6*dxvr[1]))+0.05*(12.24744871391589*alphaDrag[1]*favg[9]+(6.324555320336761*alphaDrag[2]+7.071067811865476*alphaDrag[0])*favg[8]+7.071067811865476*alphaDrag[1]*favg[4])-0.5*fjump[8]; 
  Ghat[10] = (84.85281374238573*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[10]-84.85281374238573*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[10])/(dxvr1R5+5.0*dxvl[1]*dxvr1R4+10.0*dxvl1R2*dxvr1R3+10.0*dxvl1R3*dxvr1R2+5.0*dxvl1R4*dxvr[1]+dxvl1R5)-0.5*fjump[10]+0.3535533905932737*alphaDrag[0]*favg[10]; 
  Ghat[11] = ((1859.419264179008*dxvl1R3*dxvr1R3+1062.525293816576*dxvl1R4*dxvr1R2-464.8548160447518*dxvl1R5*dxvr[1]-332.0391543176799*dxvl1R6)*nuVtSqSum[2]*fr[12]+(332.0391543176799*dxvr1R6+464.8548160447518*dxvl[1]*dxvr1R5-1062.525293816576*dxvl1R2*dxvr1R4-1859.419264179008*dxvl1R3*dxvr1R3)*nuVtSqSum[2]*fl[12]+(1897.366596101029*dxvl1R3*dxvr1R3*nuVtSqSum[2]+2969.848480983501*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fr[11]+((-1897.366596101029*dxvl1R3*dxvr1R3*nuVtSqSum[2])-2969.848480983501*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[11]+((-3220.608638130377*dxvl1R3*dxvr1R3)-613.4492644057864*dxvl1R4*dxvr1R2+383.4057902536164*dxvl1R5*dxvr[1]+76.68115805072327*dxvl1R6)*nuVtSqSum[1]*fr[5]+(76.68115805072327*dxvr1R6+383.4057902536164*dxvl[1]*dxvr1R5-613.4492644057864*dxvl1R2*dxvr1R4-3220.608638130377*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[5]+(((-3600.749921891272*dxvl1R3*dxvr1R3)-685.8571279792902*dxvl1R4*dxvr1R2+428.6607049870561*dxvl1R5*dxvr[1]+85.73214099741124*dxvl1R6)*fr[2]+(85.73214099741124*dxvr1R6+428.6607049870561*dxvl[1]*dxvr1R5-685.8571279792902*dxvl1R2*dxvr1R4-3600.749921891272*dxvl1R3*dxvr1R3)*fl[2]+(2969.848480983501*fr[0]-2969.848480983501*fl[0])*dxvl1R3*dxvr1R3)*nuVtSqSum[2]+(2656.313234541441*dxvl1R3*dxvr1R3*fr[1]-2656.313234541441*dxvl1R3*dxvr1R3*fl[1])*nuVtSqSum[1])/(35.0*dxvl[1]*dxvr1R6+175.0*dxvl1R2*dxvr1R5+350.0*dxvl1R3*dxvr1R4+350.0*dxvl1R4*dxvr1R3+175.0*dxvl1R5*dxvr1R2+35.0*dxvl1R6*dxvr[1])+0.007142857142857143*(110.6797181058933*alphaDrag[2]*favg[12]+(31.62277660168381*alphaDrag[2]+49.49747468305833*alphaDrag[0])*favg[11]+76.68115805072327*alphaDrag[1]*favg[5]+85.73214099741124*alphaDrag[2]*favg[2]+49.49747468305833*favg[0]*alphaDrag[2]+44.27188724235732*alphaDrag[1]*favg[1])-0.5*fjump[11]; 
  Ghat[13] = (84.85281374238573*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[13]-84.85281374238573*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[13])/(dxvr1R5+5.0*dxvl[1]*dxvr1R4+10.0*dxvl1R2*dxvr1R3+10.0*dxvl1R3*dxvr1R2+5.0*dxvl1R4*dxvr[1]+dxvl1R5)-0.5*fjump[13]+0.3535533905932737*alphaDrag[0]*favg[13]; 
  Ghat[14] = (84.85281374238573*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fr[14]-84.85281374238573*nuVtSqSum[0]*dxvl1R2*dxvr1R2*fl[14])/(dxvr1R5+5.0*dxvl[1]*dxvr1R4+10.0*dxvl1R2*dxvr1R3+10.0*dxvl1R3*dxvr1R2+5.0*dxvl1R4*dxvr[1]+dxvl1R5)-0.5*fjump[14]+0.3535533905932737*alphaDrag[0]*favg[14]; 

  double incr1[15]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = 0.8660254037844386*Ghat[3]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = 0.8660254037844386*Ghat[4]; 
  incr1[10] = -0.5*Ghat[10]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = -1.118033988749895*Ghat[0]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = -0.5*Ghat[14]; 

  double incr2[15]; 

  incr2[2] = ((76.68115805072327*nuVtSqSum[0]*dxvl1R3*dxvr1R2+87.63560920082662*nuVtSqSum[0]*dxvl1R4*dxvr[1]+27.38612787525831*nuVtSqSum[0]*dxvl1R5)*fr[12]+(27.38612787525831*nuVtSqSum[0]*dxvr1R5+87.63560920082662*nuVtSqSum[0]*dxvl[1]*dxvr1R4+76.68115805072327*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[12]+(122.4744871391589*dxvl1R3*dxvr1R2+61.23724356957945*dxvl1R4*dxvr[1]+12.24744871391589*dxvl1R5)*nuVtSqSum[2]*fr[11]+(12.24744871391589*dxvr1R5+61.23724356957945*dxvl[1]*dxvr1R4+122.4744871391589*dxvl1R2*dxvr1R3)*nuVtSqSum[2]*fl[11]+((-148.492424049175*dxvl1R3*dxvr1R2)-106.0660171779821*dxvl1R4*dxvr[1]-21.21320343559643*dxvl1R5)*nuVtSqSum[1]*fr[5]+(21.21320343559643*dxvr1R5+106.0660171779821*dxvl[1]*dxvr1R4+148.492424049175*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[5]+((-148.492424049175*nuVtSqSum[0]*dxvl1R3*dxvr1R2)-106.0660171779821*nuVtSqSum[0]*dxvl1R4*dxvr[1]-21.21320343559643*nuVtSqSum[0]*dxvl1R5)*fr[2]+(21.21320343559643*nuVtSqSum[0]*dxvr1R5+106.0660171779821*nuVtSqSum[0]*dxvl[1]*dxvr1R4+148.492424049175*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[2]+((122.4744871391589*dxvl1R3*dxvr1R2+61.23724356957945*dxvl1R4*dxvr[1]+12.24744871391589*dxvl1R5)*fr[1]+(12.24744871391589*dxvr1R5+61.23724356957945*dxvl[1]*dxvr1R4+122.4744871391589*dxvl1R2*dxvr1R3)*fl[1])*nuVtSqSum[1]+12.24744871391589*fl[0]*nuVtSqSum[0]*dxvr1R5+61.23724356957945*fl[0]*nuVtSqSum[0]*dxvl[1]*dxvr1R4+122.4744871391589*fl[0]*nuVtSqSum[0]*dxvl1R2*dxvr1R3+122.4744871391589*fr[0]*nuVtSqSum[0]*dxvl1R3*dxvr1R2+61.23724356957945*fr[0]*nuVtSqSum[0]*dxvl1R4*dxvr[1]+12.24744871391589*fr[0]*nuVtSqSum[0]*dxvl1R5)/(20.0*dxvr1R5+100.0*dxvl[1]*dxvr1R4+200.0*dxvl1R2*dxvr1R3+200.0*dxvl1R3*dxvr1R2+100.0*dxvl1R4*dxvr[1]+20.0*dxvl1R5); 
  incr2[5] = ((76.68115805072327*dxvl1R3*dxvr1R2+87.63560920082662*dxvl1R4*dxvr[1]+27.38612787525831*dxvl1R5)*nuVtSqSum[1]*fr[12]+(27.38612787525831*dxvr1R5+87.63560920082662*dxvl[1]*dxvr1R4+76.68115805072327*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[12]+(109.5445115010333*dxvl1R3*dxvr1R2+54.77225575051663*dxvl1R4*dxvr[1]+10.95445115010332*dxvl1R5)*nuVtSqSum[1]*fr[11]+(10.95445115010332*dxvr1R5+54.77225575051663*dxvl[1]*dxvr1R4+109.5445115010333*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[11]+(((-132.815661727072*dxvl1R3*dxvr1R2)-94.86832980505142*dxvl1R4*dxvr[1]-18.97366596101028*dxvl1R5)*nuVtSqSum[2]-148.492424049175*nuVtSqSum[0]*dxvl1R3*dxvr1R2-106.0660171779821*nuVtSqSum[0]*dxvl1R4*dxvr[1]-21.21320343559643*nuVtSqSum[0]*dxvl1R5)*fr[5]+((18.97366596101028*dxvr1R5+94.86832980505142*dxvl[1]*dxvr1R4+132.815661727072*dxvl1R2*dxvr1R3)*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0]*dxvr1R5+106.0660171779821*nuVtSqSum[0]*dxvl[1]*dxvr1R4+148.492424049175*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[5]+((109.5445115010333*dxvl1R3*dxvr1R2+54.77225575051663*dxvl1R4*dxvr[1]+10.95445115010332*dxvl1R5)*fr[1]+(10.95445115010332*dxvr1R5+54.77225575051663*dxvl[1]*dxvr1R4+109.5445115010333*dxvl1R2*dxvr1R3)*fl[1])*nuVtSqSum[2]+((-148.492424049175*dxvl1R3*dxvr1R2)-106.0660171779821*dxvl1R4*dxvr[1]-21.21320343559643*dxvl1R5)*nuVtSqSum[1]*fr[2]+(21.21320343559643*dxvr1R5+106.0660171779821*dxvl[1]*dxvr1R4+148.492424049175*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[2]+(12.24744871391589*fl[0]*dxvr1R5+61.23724356957945*fl[0]*dxvl[1]*dxvr1R4+122.4744871391589*fl[0]*dxvl1R2*dxvr1R3+122.4744871391589*fr[0]*dxvl1R3*dxvr1R2+61.23724356957945*fr[0]*dxvl1R4*dxvr[1]+12.24744871391589*fr[0]*dxvl1R5)*nuVtSqSum[1]+(122.4744871391589*nuVtSqSum[0]*dxvl1R3*dxvr1R2+61.23724356957945*nuVtSqSum[0]*dxvl1R4*dxvr[1]+12.24744871391589*nuVtSqSum[0]*dxvl1R5)*fr[1]+(12.24744871391589*nuVtSqSum[0]*dxvr1R5+61.23724356957945*nuVtSqSum[0]*dxvl[1]*dxvr1R4+122.4744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[1])/(20.0*dxvr1R5+100.0*dxvl[1]*dxvr1R4+200.0*dxvl1R2*dxvr1R3+200.0*dxvl1R3*dxvr1R2+100.0*dxvl1R4*dxvr[1]+20.0*dxvl1R5); 
  incr2[7] = -(1.0*((29.698484809835*nuVtSqSum[0]*dxvl1R3*dxvr1R2+21.21320343559643*nuVtSqSum[0]*dxvl1R4*dxvr[1]+4.242640687119286*nuVtSqSum[0]*dxvl1R5)*fr[7]+((-4.242640687119286*nuVtSqSum[0]*dxvr1R5)-21.21320343559643*nuVtSqSum[0]*dxvl[1]*dxvr1R4-29.698484809835*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[7]+((-24.49489742783179*dxvl1R3*dxvr1R2)-12.24744871391589*dxvl1R4*dxvr[1]-2.449489742783178*dxvl1R5)*nuVtSqSum[1]*fr[6]+((-2.449489742783178*dxvr1R5)-12.24744871391589*dxvl[1]*dxvr1R4-24.49489742783179*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[6]+((-24.49489742783179*nuVtSqSum[0]*dxvl1R3*dxvr1R2)-12.24744871391589*nuVtSqSum[0]*dxvl1R4*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R5)*fr[3]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R5)-12.24744871391589*nuVtSqSum[0]*dxvl[1]*dxvr1R4-24.49489742783179*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[3]))/(4.0*dxvr1R5+20.0*dxvl[1]*dxvr1R4+40.0*dxvl1R2*dxvr1R3+40.0*dxvl1R3*dxvr1R2+20.0*dxvl1R4*dxvr[1]+4.0*dxvl1R5); 
  incr2[9] = -(1.0*((29.698484809835*nuVtSqSum[0]*dxvl1R3*dxvr1R2+21.21320343559643*nuVtSqSum[0]*dxvl1R4*dxvr[1]+4.242640687119286*nuVtSqSum[0]*dxvl1R5)*fr[9]+((-4.242640687119286*nuVtSqSum[0]*dxvr1R5)-21.21320343559643*nuVtSqSum[0]*dxvl[1]*dxvr1R4-29.698484809835*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[9]+((-24.49489742783179*dxvl1R3*dxvr1R2)-12.24744871391589*dxvl1R4*dxvr[1]-2.449489742783178*dxvl1R5)*nuVtSqSum[1]*fr[8]+((-2.449489742783178*dxvr1R5)-12.24744871391589*dxvl[1]*dxvr1R4-24.49489742783179*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[8]+((-24.49489742783179*nuVtSqSum[0]*dxvl1R3*dxvr1R2)-12.24744871391589*nuVtSqSum[0]*dxvl1R4*dxvr[1]-2.449489742783178*nuVtSqSum[0]*dxvl1R5)*fr[4]+((-2.449489742783178*nuVtSqSum[0]*dxvr1R5)-12.24744871391589*nuVtSqSum[0]*dxvl[1]*dxvr1R4-24.49489742783179*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[4]))/(4.0*dxvr1R5+20.0*dxvl[1]*dxvr1R4+40.0*dxvl1R2*dxvr1R3+40.0*dxvl1R3*dxvr1R2+20.0*dxvl1R4*dxvr[1]+4.0*dxvl1R5); 
  incr2[12] = -(1.0*((59.39696961967*nuVtSqSum[0]*dxvl1R3*dxvr1R2+67.8822509939086*nuVtSqSum[0]*dxvl1R4*dxvr[1]+21.21320343559643*nuVtSqSum[0]*dxvl1R5)*fr[12]+(21.21320343559643*nuVtSqSum[0]*dxvr1R5+67.8822509939086*nuVtSqSum[0]*dxvl[1]*dxvr1R4+59.39696961967*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[12]+(94.86832980505142*dxvl1R3*dxvr1R2+47.43416490252571*dxvl1R4*dxvr[1]+9.48683298050514*dxvl1R5)*nuVtSqSum[2]*fr[11]+(9.48683298050514*dxvr1R5+47.43416490252571*dxvl[1]*dxvr1R4+94.86832980505142*dxvl1R2*dxvr1R3)*nuVtSqSum[2]*fl[11]+((-115.0217370760849*dxvl1R3*dxvr1R2)-82.15838362577493*dxvl1R4*dxvr[1]-16.43167672515498*dxvl1R5)*nuVtSqSum[1]*fr[5]+(16.43167672515498*dxvr1R5+82.15838362577493*dxvl[1]*dxvr1R4+115.0217370760849*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[5]+((-115.0217370760849*nuVtSqSum[0]*dxvl1R3*dxvr1R2)-82.15838362577493*nuVtSqSum[0]*dxvl1R4*dxvr[1]-16.43167672515498*nuVtSqSum[0]*dxvl1R5)*fr[2]+(16.43167672515498*nuVtSqSum[0]*dxvr1R5+82.15838362577493*nuVtSqSum[0]*dxvl[1]*dxvr1R4+115.0217370760849*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[2]+((94.86832980505142*dxvl1R3*dxvr1R2+47.43416490252571*dxvl1R4*dxvr[1]+9.48683298050514*dxvl1R5)*fr[1]+(9.48683298050514*dxvr1R5+47.43416490252571*dxvl[1]*dxvr1R4+94.86832980505142*dxvl1R2*dxvr1R3)*fl[1])*nuVtSqSum[1]+9.48683298050514*fl[0]*nuVtSqSum[0]*dxvr1R5+47.43416490252571*fl[0]*nuVtSqSum[0]*dxvl[1]*dxvr1R4+94.86832980505142*fl[0]*nuVtSqSum[0]*dxvl1R2*dxvr1R3+94.86832980505142*fr[0]*nuVtSqSum[0]*dxvl1R3*dxvr1R2+47.43416490252571*fr[0]*nuVtSqSum[0]*dxvl1R4*dxvr[1]+9.48683298050514*fr[0]*nuVtSqSum[0]*dxvl1R5))/(4.0*dxvr1R5+20.0*dxvl[1]*dxvr1R4+40.0*dxvl1R2*dxvr1R3+40.0*dxvl1R3*dxvr1R2+20.0*dxvl1R4*dxvr[1]+4.0*dxvl1R5); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 
  outr[10] += incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr2[12]*rdvSq4R+incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr1[14]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += incr1[9]*rdv2L-1.0*incr2[9]*rdvSq4L; 
  outl[10] += -1.0*incr1[10]*rdv2L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += incr2[12]*rdvSq4L-1.0*incr1[12]*rdv2L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += -1.0*incr1[14]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUx[2])/nuSum-(0.7071067811865475*sumNuUx[0])/nuSum+wl[1]); 
} 
double VmLBOconstNuSurfNonUniform1x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[9]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUy = &nuUSum[3]; 

  double favg[15]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 

  double fjump[15]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(-1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[2]*nuSum+0.7071067811865475*dxvl[2]*nuSum-1.0*sumNuUy[0]; 
  alphaDrag[1] = -1.0*sumNuUy[1]; 
  alphaDrag[2] = -1.0*sumNuUy[2]; 

  double Ghat[15]; 
  for(unsigned short int i=0; i<15; ++i){ 
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

  Ghat[0] = ((265.631323454144*nuVtSqSum[0]*dxvl2R3*dxvr2R3+151.7893276880823*nuVtSqSum[0]*dxvl2R4*dxvr2R2-66.40783086353598*nuVtSqSum[0]*dxvl2R5*dxvr[2]-47.43416490252571*nuVtSqSum[0]*dxvl2R6)*fr[13]+(47.43416490252571*nuVtSqSum[0]*dxvr2R6+66.40783086353598*nuVtSqSum[0]*dxvl[2]*dxvr2R5-151.7893276880823*nuVtSqSum[0]*dxvl2R2*dxvr2R4-265.631323454144*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[13]+424.2640687119287*dxvl2R3*dxvr2R3*nuVtSqSum[2]*fr[11]-424.2640687119287*dxvl2R3*dxvr2R3*nuVtSqSum[2]*fl[11]+((-514.3928459844674*nuVtSqSum[1]*dxvl2R3*dxvr2R3)-97.97958971132716*nuVtSqSum[1]*dxvl2R4*dxvr2R2+61.23724356957945*nuVtSqSum[1]*dxvl2R5*dxvr[2]+12.24744871391589*nuVtSqSum[1]*dxvl2R6)*fr[6]+(12.24744871391589*nuVtSqSum[1]*dxvr2R6+61.23724356957945*nuVtSqSum[1]*dxvl[2]*dxvr2R5-97.97958971132716*nuVtSqSum[1]*dxvl2R2*dxvr2R4-514.3928459844674*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[6]+((-514.3928459844674*nuVtSqSum[0]*dxvl2R3*dxvr2R3)-97.97958971132716*nuVtSqSum[0]*dxvl2R4*dxvr2R2+61.23724356957945*nuVtSqSum[0]*dxvl2R5*dxvr[2]+12.24744871391589*nuVtSqSum[0]*dxvl2R6)*fr[3]+(12.24744871391589*nuVtSqSum[0]*dxvr2R6+61.23724356957945*nuVtSqSum[0]*dxvl[2]*dxvr2R5-97.97958971132716*nuVtSqSum[0]*dxvl2R2*dxvr2R4-514.3928459844674*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[3]+((424.2640687119287*fr[1]-424.2640687119287*fl[1])*nuVtSqSum[1]+(424.2640687119287*fr[0]-424.2640687119287*fl[0])*nuVtSqSum[0])*dxvl2R3*dxvr2R3)/(5.0*dxvl[2]*dxvr2R6+25.0*dxvl2R2*dxvr2R5+50.0*dxvl2R3*dxvr2R4+50.0*dxvl2R4*dxvr2R3+25.0*dxvl2R5*dxvr2R2+5.0*dxvl2R6*dxvr[2])-0.5*(2.23606797749979*fjump[13]+1.732050807568877*fjump[3]+fjump[0])+0.25*(3.16227766016838*alphaDrag[0]*favg[13]+1.414213562373095*alphaDrag[2]*favg[11]+2.449489742783178*alphaDrag[1]*favg[6]+2.449489742783178*alphaDrag[0]*favg[3]+1.414213562373095*alphaDrag[1]*favg[1]+1.414213562373095*alphaDrag[0]*favg[0]); 
  Ghat[1] = ((265.631323454144*nuVtSqSum[1]*dxvl2R3*dxvr2R3+151.7893276880823*nuVtSqSum[1]*dxvl2R4*dxvr2R2-66.40783086353598*nuVtSqSum[1]*dxvl2R5*dxvr[2]-47.43416490252571*nuVtSqSum[1]*dxvl2R6)*fr[13]+(47.43416490252571*nuVtSqSum[1]*dxvr2R6+66.40783086353598*nuVtSqSum[1]*dxvl[2]*dxvr2R5-151.7893276880823*nuVtSqSum[1]*dxvl2R2*dxvr2R4-265.631323454144*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[13]+379.4733192202058*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fr[11]-379.4733192202058*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fl[11]+(((-460.0869483043397*dxvl2R3*dxvr2R3)-87.63560920082662*dxvl2R4*dxvr2R2+54.77225575051663*dxvl2R5*dxvr[2]+10.95445115010332*dxvl2R6)*nuVtSqSum[2]-514.3928459844674*nuVtSqSum[0]*dxvl2R3*dxvr2R3-97.97958971132716*nuVtSqSum[0]*dxvl2R4*dxvr2R2+61.23724356957945*nuVtSqSum[0]*dxvl2R5*dxvr[2]+12.24744871391589*nuVtSqSum[0]*dxvl2R6)*fr[6]+((10.95445115010332*dxvr2R6+54.77225575051663*dxvl[2]*dxvr2R5-87.63560920082662*dxvl2R2*dxvr2R4-460.0869483043397*dxvl2R3*dxvr2R3)*nuVtSqSum[2]+12.24744871391589*nuVtSqSum[0]*dxvr2R6+61.23724356957945*nuVtSqSum[0]*dxvl[2]*dxvr2R5-97.97958971132716*nuVtSqSum[0]*dxvl2R2*dxvr2R4-514.3928459844674*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[6]+((-514.3928459844674*nuVtSqSum[1]*dxvl2R3*dxvr2R3)-97.97958971132716*nuVtSqSum[1]*dxvl2R4*dxvr2R2+61.23724356957945*nuVtSqSum[1]*dxvl2R5*dxvr[2]+12.24744871391589*nuVtSqSum[1]*dxvl2R6)*fr[3]+(12.24744871391589*nuVtSqSum[1]*dxvr2R6+61.23724356957945*nuVtSqSum[1]*dxvl[2]*dxvr2R5-97.97958971132716*nuVtSqSum[1]*dxvl2R2*dxvr2R4-514.3928459844674*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[3]+(379.4733192202058*fr[1]-379.4733192202058*fl[1])*dxvl2R3*dxvr2R3*nuVtSqSum[2]+((424.2640687119287*fr[0]-424.2640687119287*fl[0])*nuVtSqSum[1]+424.2640687119287*nuVtSqSum[0]*fr[1]-424.2640687119287*nuVtSqSum[0]*fl[1])*dxvl2R3*dxvr2R3)/(5.0*dxvl[2]*dxvr2R6+25.0*dxvl2R2*dxvr2R5+50.0*dxvl2R3*dxvr2R4+50.0*dxvl2R4*dxvr2R3+25.0*dxvl2R5*dxvr2R2+5.0*dxvl2R6*dxvr[2])+0.05*(15.8113883008419*alphaDrag[1]*favg[13]+6.324555320336761*alphaDrag[1]*favg[11]+(10.95445115010332*alphaDrag[2]+12.24744871391589*alphaDrag[0])*favg[6]+12.24744871391589*alphaDrag[1]*favg[3]+6.324555320336761*favg[1]*alphaDrag[2]+7.071067811865476*alphaDrag[0]*favg[1]+7.071067811865476*favg[0]*alphaDrag[1])-0.5*(1.732050807568877*fjump[6]+fjump[1]); 
  Ghat[2] = (-(1.0*((102.8785691968935*nuVtSqSum[0]*dxvl2R3*dxvr2R3+19.59591794226543*nuVtSqSum[0]*dxvl2R4*dxvr2R2-12.24744871391589*nuVtSqSum[0]*dxvl2R5*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R6)*fr[7]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R6)-12.24744871391589*nuVtSqSum[0]*dxvl[2]*dxvr2R5+19.59591794226543*nuVtSqSum[0]*dxvl2R2*dxvr2R4+102.8785691968935*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[7]-84.85281374238573*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fr[5]+84.85281374238573*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fl[5]-84.85281374238573*nuVtSqSum[0]*dxvl2R3*dxvr2R3*fr[2]+84.85281374238573*nuVtSqSum[0]*dxvl2R3*dxvr2R3*fl[2]))/(dxvl[2]*dxvr2R6+5.0*dxvl2R2*dxvr2R5+10.0*dxvl2R3*dxvr2R4+10.0*dxvl2R4*dxvr2R3+5.0*dxvl2R5*dxvr2R2+dxvl2R6*dxvr[2]))-0.5*(1.732050807568877*fjump[7]+fjump[2])+0.25*(2.449489742783178*alphaDrag[0]*favg[7]+1.414213562373095*alphaDrag[1]*favg[5]+1.414213562373095*alphaDrag[0]*favg[2]); 
  Ghat[4] = (-(1.0*((102.8785691968935*nuVtSqSum[0]*dxvl2R3*dxvr2R3+19.59591794226543*nuVtSqSum[0]*dxvl2R4*dxvr2R2-12.24744871391589*nuVtSqSum[0]*dxvl2R5*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R6)*fr[10]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R6)-12.24744871391589*nuVtSqSum[0]*dxvl[2]*dxvr2R5+19.59591794226543*nuVtSqSum[0]*dxvl2R2*dxvr2R4+102.8785691968935*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[10]-84.85281374238573*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fr[8]+84.85281374238573*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fl[8]-84.85281374238573*nuVtSqSum[0]*dxvl2R3*dxvr2R3*fr[4]+84.85281374238573*nuVtSqSum[0]*dxvl2R3*dxvr2R3*fl[4]))/(dxvl[2]*dxvr2R6+5.0*dxvl2R2*dxvr2R5+10.0*dxvl2R3*dxvr2R4+10.0*dxvl2R4*dxvr2R3+5.0*dxvl2R5*dxvr2R2+dxvl2R6*dxvr[2]))-0.5*(1.732050807568877*fjump[10]+fjump[4])+0.25*(2.449489742783178*alphaDrag[0]*favg[10]+1.414213562373095*alphaDrag[1]*favg[8]+1.414213562373095*alphaDrag[0]*favg[4]); 
  Ghat[5] = (-(1.0*((102.8785691968935*nuVtSqSum[1]*dxvl2R3*dxvr2R3+19.59591794226543*nuVtSqSum[1]*dxvl2R4*dxvr2R2-12.24744871391589*nuVtSqSum[1]*dxvl2R5*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R6)*fr[7]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R6)-12.24744871391589*nuVtSqSum[1]*dxvl[2]*dxvr2R5+19.59591794226543*nuVtSqSum[1]*dxvl2R2*dxvr2R4+102.8785691968935*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[7]+((-75.89466384404115*dxvl2R3*dxvr2R3*nuVtSqSum[2])-84.85281374238573*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fr[5]+(75.89466384404115*dxvl2R3*dxvr2R3*nuVtSqSum[2]+84.85281374238573*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[5]-84.85281374238573*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fr[2]+84.85281374238573*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fl[2]))/(dxvl[2]*dxvr2R6+5.0*dxvl2R2*dxvr2R5+10.0*dxvl2R3*dxvr2R4+10.0*dxvl2R4*dxvr2R3+5.0*dxvl2R5*dxvr2R2+dxvl2R6*dxvr[2]))+0.05*(12.24744871391589*alphaDrag[1]*favg[7]+(6.324555320336761*alphaDrag[2]+7.071067811865476*alphaDrag[0])*favg[5]+7.071067811865476*alphaDrag[1]*favg[2])-0.5*fjump[5]; 
  Ghat[8] = (-(1.0*((102.8785691968935*nuVtSqSum[1]*dxvl2R3*dxvr2R3+19.59591794226543*nuVtSqSum[1]*dxvl2R4*dxvr2R2-12.24744871391589*nuVtSqSum[1]*dxvl2R5*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R6)*fr[10]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R6)-12.24744871391589*nuVtSqSum[1]*dxvl[2]*dxvr2R5+19.59591794226543*nuVtSqSum[1]*dxvl2R2*dxvr2R4+102.8785691968935*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[10]+((-75.89466384404115*dxvl2R3*dxvr2R3*nuVtSqSum[2])-84.85281374238573*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fr[8]+(75.89466384404115*dxvl2R3*dxvr2R3*nuVtSqSum[2]+84.85281374238573*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[8]-84.85281374238573*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fr[4]+84.85281374238573*nuVtSqSum[1]*dxvl2R3*dxvr2R3*fl[4]))/(dxvl[2]*dxvr2R6+5.0*dxvl2R2*dxvr2R5+10.0*dxvl2R3*dxvr2R4+10.0*dxvl2R4*dxvr2R3+5.0*dxvl2R5*dxvr2R2+dxvl2R6*dxvr[2]))+0.05*(12.24744871391589*alphaDrag[1]*favg[10]+(6.324555320336761*alphaDrag[2]+7.071067811865476*alphaDrag[0])*favg[8]+7.071067811865476*alphaDrag[1]*favg[4])-0.5*fjump[8]; 
  Ghat[9] = (84.85281374238573*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[9]-84.85281374238573*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[9])/(dxvr2R5+5.0*dxvl[2]*dxvr2R4+10.0*dxvl2R2*dxvr2R3+10.0*dxvl2R3*dxvr2R2+5.0*dxvl2R4*dxvr[2]+dxvl2R5)-0.5*fjump[9]+0.3535533905932737*alphaDrag[0]*favg[9]; 
  Ghat[11] = ((1859.419264179008*dxvl2R3*dxvr2R3+1062.525293816576*dxvl2R4*dxvr2R2-464.8548160447518*dxvl2R5*dxvr[2]-332.0391543176799*dxvl2R6)*nuVtSqSum[2]*fr[13]+(332.0391543176799*dxvr2R6+464.8548160447518*dxvl[2]*dxvr2R5-1062.525293816576*dxvl2R2*dxvr2R4-1859.419264179008*dxvl2R3*dxvr2R3)*nuVtSqSum[2]*fl[13]+(1897.366596101029*dxvl2R3*dxvr2R3*nuVtSqSum[2]+2969.848480983501*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fr[11]+((-1897.366596101029*dxvl2R3*dxvr2R3*nuVtSqSum[2])-2969.848480983501*nuVtSqSum[0]*dxvl2R3*dxvr2R3)*fl[11]+((-3220.608638130377*nuVtSqSum[1]*dxvl2R3*dxvr2R3)-613.4492644057864*nuVtSqSum[1]*dxvl2R4*dxvr2R2+383.4057902536164*nuVtSqSum[1]*dxvl2R5*dxvr[2]+76.68115805072327*nuVtSqSum[1]*dxvl2R6)*fr[6]+(76.68115805072327*nuVtSqSum[1]*dxvr2R6+383.4057902536164*nuVtSqSum[1]*dxvl[2]*dxvr2R5-613.4492644057864*nuVtSqSum[1]*dxvl2R2*dxvr2R4-3220.608638130377*nuVtSqSum[1]*dxvl2R3*dxvr2R3)*fl[6]+((-3600.749921891272*dxvl2R3*dxvr2R3)-685.8571279792902*dxvl2R4*dxvr2R2+428.6607049870561*dxvl2R5*dxvr[2]+85.73214099741124*dxvl2R6)*nuVtSqSum[2]*fr[3]+(85.73214099741124*dxvr2R6+428.6607049870561*dxvl[2]*dxvr2R5-685.8571279792902*dxvl2R2*dxvr2R4-3600.749921891272*dxvl2R3*dxvr2R3)*nuVtSqSum[2]*fl[3]+(2969.848480983501*fr[0]-2969.848480983501*fl[0])*dxvl2R3*dxvr2R3*nuVtSqSum[2]+(2656.313234541441*fr[1]-2656.313234541441*fl[1])*nuVtSqSum[1]*dxvl2R3*dxvr2R3)/(35.0*dxvl[2]*dxvr2R6+175.0*dxvl2R2*dxvr2R5+350.0*dxvl2R3*dxvr2R4+350.0*dxvl2R4*dxvr2R3+175.0*dxvl2R5*dxvr2R2+35.0*dxvl2R6*dxvr[2])+0.007142857142857143*(110.6797181058933*alphaDrag[2]*favg[13]+(31.62277660168381*alphaDrag[2]+49.49747468305833*alphaDrag[0])*favg[11]+76.68115805072327*alphaDrag[1]*favg[6]+85.73214099741124*alphaDrag[2]*favg[3]+49.49747468305833*favg[0]*alphaDrag[2]+44.27188724235732*alphaDrag[1]*favg[1])-0.5*fjump[11]; 
  Ghat[12] = (84.85281374238573*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[12]-84.85281374238573*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[12])/(dxvr2R5+5.0*dxvl[2]*dxvr2R4+10.0*dxvl2R2*dxvr2R3+10.0*dxvl2R3*dxvr2R2+5.0*dxvl2R4*dxvr[2]+dxvl2R5)-0.5*fjump[12]+0.3535533905932737*alphaDrag[0]*favg[12]; 
  Ghat[14] = (84.85281374238573*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fr[14]-84.85281374238573*nuVtSqSum[0]*dxvl2R2*dxvr2R2*fl[14])/(dxvr2R5+5.0*dxvl[2]*dxvr2R4+10.0*dxvl2R2*dxvr2R3+10.0*dxvl2R3*dxvr2R2+5.0*dxvl2R4*dxvr[2]+dxvl2R5)-0.5*fjump[14]+0.3535533905932737*alphaDrag[0]*favg[14]; 

  double incr1[15]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[1]; 
  incr1[7] = 0.8660254037844386*Ghat[2]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = -0.5*Ghat[9]; 
  incr1[10] = 0.8660254037844386*Ghat[4]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = -1.118033988749895*Ghat[0]; 
  incr1[14] = -0.5*Ghat[14]; 

  double incr2[15]; 

  incr2[3] = ((76.68115805072327*nuVtSqSum[0]*dxvl2R3*dxvr2R2+87.63560920082662*nuVtSqSum[0]*dxvl2R4*dxvr[2]+27.38612787525831*nuVtSqSum[0]*dxvl2R5)*fr[13]+(27.38612787525831*nuVtSqSum[0]*dxvr2R5+87.63560920082662*nuVtSqSum[0]*dxvl[2]*dxvr2R4+76.68115805072327*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[13]+(122.4744871391589*dxvl2R3*dxvr2R2+61.23724356957945*dxvl2R4*dxvr[2]+12.24744871391589*dxvl2R5)*nuVtSqSum[2]*fr[11]+(12.24744871391589*dxvr2R5+61.23724356957945*dxvl[2]*dxvr2R4+122.4744871391589*dxvl2R2*dxvr2R3)*nuVtSqSum[2]*fl[11]+((-148.492424049175*nuVtSqSum[1]*dxvl2R3*dxvr2R2)-106.0660171779821*nuVtSqSum[1]*dxvl2R4*dxvr[2]-21.21320343559643*nuVtSqSum[1]*dxvl2R5)*fr[6]+(21.21320343559643*nuVtSqSum[1]*dxvr2R5+106.0660171779821*nuVtSqSum[1]*dxvl[2]*dxvr2R4+148.492424049175*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[6]+((-148.492424049175*nuVtSqSum[0]*dxvl2R3*dxvr2R2)-106.0660171779821*nuVtSqSum[0]*dxvl2R4*dxvr[2]-21.21320343559643*nuVtSqSum[0]*dxvl2R5)*fr[3]+(21.21320343559643*nuVtSqSum[0]*dxvr2R5+106.0660171779821*nuVtSqSum[0]*dxvl[2]*dxvr2R4+148.492424049175*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[3]+(12.24744871391589*fl[1]*nuVtSqSum[1]+12.24744871391589*fl[0]*nuVtSqSum[0])*dxvr2R5+(61.23724356957945*fl[1]*nuVtSqSum[1]+61.23724356957945*fl[0]*nuVtSqSum[0])*dxvl[2]*dxvr2R4+(122.4744871391589*fl[1]*nuVtSqSum[1]+122.4744871391589*fl[0]*nuVtSqSum[0])*dxvl2R2*dxvr2R3+(122.4744871391589*fr[1]*nuVtSqSum[1]+122.4744871391589*fr[0]*nuVtSqSum[0])*dxvl2R3*dxvr2R2+(61.23724356957945*fr[1]*nuVtSqSum[1]+61.23724356957945*fr[0]*nuVtSqSum[0])*dxvl2R4*dxvr[2]+(12.24744871391589*fr[1]*nuVtSqSum[1]+12.24744871391589*fr[0]*nuVtSqSum[0])*dxvl2R5)/(20.0*dxvr2R5+100.0*dxvl[2]*dxvr2R4+200.0*dxvl2R2*dxvr2R3+200.0*dxvl2R3*dxvr2R2+100.0*dxvl2R4*dxvr[2]+20.0*dxvl2R5); 
  incr2[6] = ((76.68115805072327*nuVtSqSum[1]*dxvl2R3*dxvr2R2+87.63560920082662*nuVtSqSum[1]*dxvl2R4*dxvr[2]+27.38612787525831*nuVtSqSum[1]*dxvl2R5)*fr[13]+(27.38612787525831*nuVtSqSum[1]*dxvr2R5+87.63560920082662*nuVtSqSum[1]*dxvl[2]*dxvr2R4+76.68115805072327*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[13]+(109.5445115010333*nuVtSqSum[1]*dxvl2R3*dxvr2R2+54.77225575051663*nuVtSqSum[1]*dxvl2R4*dxvr[2]+10.95445115010332*nuVtSqSum[1]*dxvl2R5)*fr[11]+(10.95445115010332*nuVtSqSum[1]*dxvr2R5+54.77225575051663*nuVtSqSum[1]*dxvl[2]*dxvr2R4+109.5445115010333*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[11]+(((-132.815661727072*dxvl2R3*dxvr2R2)-94.86832980505142*dxvl2R4*dxvr[2]-18.97366596101028*dxvl2R5)*nuVtSqSum[2]-148.492424049175*nuVtSqSum[0]*dxvl2R3*dxvr2R2-106.0660171779821*nuVtSqSum[0]*dxvl2R4*dxvr[2]-21.21320343559643*nuVtSqSum[0]*dxvl2R5)*fr[6]+((18.97366596101028*dxvr2R5+94.86832980505142*dxvl[2]*dxvr2R4+132.815661727072*dxvl2R2*dxvr2R3)*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0]*dxvr2R5+106.0660171779821*nuVtSqSum[0]*dxvl[2]*dxvr2R4+148.492424049175*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[6]+((-148.492424049175*nuVtSqSum[1]*dxvl2R3*dxvr2R2)-106.0660171779821*nuVtSqSum[1]*dxvl2R4*dxvr[2]-21.21320343559643*nuVtSqSum[1]*dxvl2R5)*fr[3]+(21.21320343559643*nuVtSqSum[1]*dxvr2R5+106.0660171779821*nuVtSqSum[1]*dxvl[2]*dxvr2R4+148.492424049175*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[3]+(10.95445115010332*fl[1]*dxvr2R5+54.77225575051663*fl[1]*dxvl[2]*dxvr2R4+109.5445115010333*fl[1]*dxvl2R2*dxvr2R3+109.5445115010333*fr[1]*dxvl2R3*dxvr2R2+54.77225575051663*fr[1]*dxvl2R4*dxvr[2]+10.95445115010332*fr[1]*dxvl2R5)*nuVtSqSum[2]+(12.24744871391589*fl[0]*nuVtSqSum[1]+12.24744871391589*nuVtSqSum[0]*fl[1])*dxvr2R5+(61.23724356957945*fl[0]*nuVtSqSum[1]+61.23724356957945*nuVtSqSum[0]*fl[1])*dxvl[2]*dxvr2R4+(122.4744871391589*fl[0]*nuVtSqSum[1]+122.4744871391589*nuVtSqSum[0]*fl[1])*dxvl2R2*dxvr2R3+(122.4744871391589*fr[0]*nuVtSqSum[1]+122.4744871391589*nuVtSqSum[0]*fr[1])*dxvl2R3*dxvr2R2+(61.23724356957945*fr[0]*nuVtSqSum[1]+61.23724356957945*nuVtSqSum[0]*fr[1])*dxvl2R4*dxvr[2]+(12.24744871391589*fr[0]*nuVtSqSum[1]+12.24744871391589*nuVtSqSum[0]*fr[1])*dxvl2R5)/(20.0*dxvr2R5+100.0*dxvl[2]*dxvr2R4+200.0*dxvl2R2*dxvr2R3+200.0*dxvl2R3*dxvr2R2+100.0*dxvl2R4*dxvr[2]+20.0*dxvl2R5); 
  incr2[7] = -(1.0*((29.698484809835*nuVtSqSum[0]*dxvl2R3*dxvr2R2+21.21320343559643*nuVtSqSum[0]*dxvl2R4*dxvr[2]+4.242640687119286*nuVtSqSum[0]*dxvl2R5)*fr[7]+((-4.242640687119286*nuVtSqSum[0]*dxvr2R5)-21.21320343559643*nuVtSqSum[0]*dxvl[2]*dxvr2R4-29.698484809835*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[7]+((-24.49489742783179*nuVtSqSum[1]*dxvl2R3*dxvr2R2)-12.24744871391589*nuVtSqSum[1]*dxvl2R4*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R5)*fr[5]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R5)-12.24744871391589*nuVtSqSum[1]*dxvl[2]*dxvr2R4-24.49489742783179*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[5]+((-24.49489742783179*nuVtSqSum[0]*dxvl2R3*dxvr2R2)-12.24744871391589*nuVtSqSum[0]*dxvl2R4*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R5)*fr[2]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R5)-12.24744871391589*nuVtSqSum[0]*dxvl[2]*dxvr2R4-24.49489742783179*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[2]))/(4.0*dxvr2R5+20.0*dxvl[2]*dxvr2R4+40.0*dxvl2R2*dxvr2R3+40.0*dxvl2R3*dxvr2R2+20.0*dxvl2R4*dxvr[2]+4.0*dxvl2R5); 
  incr2[10] = -(1.0*((29.698484809835*nuVtSqSum[0]*dxvl2R3*dxvr2R2+21.21320343559643*nuVtSqSum[0]*dxvl2R4*dxvr[2]+4.242640687119286*nuVtSqSum[0]*dxvl2R5)*fr[10]+((-4.242640687119286*nuVtSqSum[0]*dxvr2R5)-21.21320343559643*nuVtSqSum[0]*dxvl[2]*dxvr2R4-29.698484809835*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[10]+((-24.49489742783179*nuVtSqSum[1]*dxvl2R3*dxvr2R2)-12.24744871391589*nuVtSqSum[1]*dxvl2R4*dxvr[2]-2.449489742783178*nuVtSqSum[1]*dxvl2R5)*fr[8]+((-2.449489742783178*nuVtSqSum[1]*dxvr2R5)-12.24744871391589*nuVtSqSum[1]*dxvl[2]*dxvr2R4-24.49489742783179*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[8]+((-24.49489742783179*nuVtSqSum[0]*dxvl2R3*dxvr2R2)-12.24744871391589*nuVtSqSum[0]*dxvl2R4*dxvr[2]-2.449489742783178*nuVtSqSum[0]*dxvl2R5)*fr[4]+((-2.449489742783178*nuVtSqSum[0]*dxvr2R5)-12.24744871391589*nuVtSqSum[0]*dxvl[2]*dxvr2R4-24.49489742783179*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[4]))/(4.0*dxvr2R5+20.0*dxvl[2]*dxvr2R4+40.0*dxvl2R2*dxvr2R3+40.0*dxvl2R3*dxvr2R2+20.0*dxvl2R4*dxvr[2]+4.0*dxvl2R5); 
  incr2[13] = -(1.0*((59.39696961967*nuVtSqSum[0]*dxvl2R3*dxvr2R2+67.8822509939086*nuVtSqSum[0]*dxvl2R4*dxvr[2]+21.21320343559643*nuVtSqSum[0]*dxvl2R5)*fr[13]+(21.21320343559643*nuVtSqSum[0]*dxvr2R5+67.8822509939086*nuVtSqSum[0]*dxvl[2]*dxvr2R4+59.39696961967*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[13]+(94.86832980505142*dxvl2R3*dxvr2R2+47.43416490252571*dxvl2R4*dxvr[2]+9.48683298050514*dxvl2R5)*nuVtSqSum[2]*fr[11]+(9.48683298050514*dxvr2R5+47.43416490252571*dxvl[2]*dxvr2R4+94.86832980505142*dxvl2R2*dxvr2R3)*nuVtSqSum[2]*fl[11]+((-115.0217370760849*nuVtSqSum[1]*dxvl2R3*dxvr2R2)-82.15838362577493*nuVtSqSum[1]*dxvl2R4*dxvr[2]-16.43167672515498*nuVtSqSum[1]*dxvl2R5)*fr[6]+(16.43167672515498*nuVtSqSum[1]*dxvr2R5+82.15838362577493*nuVtSqSum[1]*dxvl[2]*dxvr2R4+115.0217370760849*nuVtSqSum[1]*dxvl2R2*dxvr2R3)*fl[6]+((-115.0217370760849*nuVtSqSum[0]*dxvl2R3*dxvr2R2)-82.15838362577493*nuVtSqSum[0]*dxvl2R4*dxvr[2]-16.43167672515498*nuVtSqSum[0]*dxvl2R5)*fr[3]+(16.43167672515498*nuVtSqSum[0]*dxvr2R5+82.15838362577493*nuVtSqSum[0]*dxvl[2]*dxvr2R4+115.0217370760849*nuVtSqSum[0]*dxvl2R2*dxvr2R3)*fl[3]+(9.48683298050514*fl[1]*nuVtSqSum[1]+9.48683298050514*fl[0]*nuVtSqSum[0])*dxvr2R5+(47.43416490252571*fl[1]*nuVtSqSum[1]+47.43416490252571*fl[0]*nuVtSqSum[0])*dxvl[2]*dxvr2R4+(94.86832980505142*fl[1]*nuVtSqSum[1]+94.86832980505142*fl[0]*nuVtSqSum[0])*dxvl2R2*dxvr2R3+(94.86832980505142*fr[1]*nuVtSqSum[1]+94.86832980505142*fr[0]*nuVtSqSum[0])*dxvl2R3*dxvr2R2+(47.43416490252571*fr[1]*nuVtSqSum[1]+47.43416490252571*fr[0]*nuVtSqSum[0])*dxvl2R4*dxvr[2]+(9.48683298050514*fr[1]*nuVtSqSum[1]+9.48683298050514*fr[0]*nuVtSqSum[0])*dxvl2R5))/(4.0*dxvr2R5+20.0*dxvl[2]*dxvr2R4+40.0*dxvl2R2*dxvr2R3+40.0*dxvl2R3*dxvr2R2+20.0*dxvl2R4*dxvr[2]+4.0*dxvl2R5); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr1[9]*rdv2R; 
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr1[14]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += -1.0*incr1[9]*rdv2L; 
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += incr2[13]*rdvSq4L-1.0*incr1[13]*rdv2L; 
  outl[14] += -1.0*incr1[14]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUy[2])/nuSum-(0.7071067811865475*sumNuUy[0])/nuSum+wl[2]); 
} 
double VmLBOconstNuSurfNonUniform1x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[9]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  const double *sumNuUz = &nuUSum[6]; 

  double favg[15]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = -1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 

  double fjump[15]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(-1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(-1*fr[8])); 
  fjump[9] = nuSum*vMuMidMax*(fl[9]-(-1*fr[9])); 
  fjump[10] = nuSum*vMuMidMax*(fl[10]-(-1*fr[10])); 
  fjump[11] = nuSum*vMuMidMax*(fl[11]-(1*fr[11])); 
  fjump[12] = nuSum*vMuMidMax*(fl[12]-(1*fr[12])); 
  fjump[13] = nuSum*vMuMidMax*(fl[13]-(1*fr[13])); 
  fjump[14] = nuSum*vMuMidMax*(fl[14]-(1*fr[14])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[3]*nuSum+0.7071067811865475*dxvl[3]*nuSum-1.0*sumNuUz[0]; 
  alphaDrag[1] = -1.0*sumNuUz[1]; 
  alphaDrag[2] = -1.0*sumNuUz[2]; 

  double Ghat[15]; 
  for(unsigned short int i=0; i<15; ++i){ 
    Ghat[i]=0.0; 
  }; 

  const double dxvl3R2 = pow(dxvl[3],2);
  const double dxvl3R3 = pow(dxvl[3],3);
  const double dxvl3R4 = pow(dxvl[3],4);
  const double dxvl3R5 = pow(dxvl[3],5);
  const double dxvr3R2 = pow(dxvr[3],2);
  const double dxvr3R3 = pow(dxvr[3],3);
  const double dxvr3R4 = pow(dxvr[3],4);
  const double dxvr3R5 = pow(dxvr[3],5);

  Ghat[0] = ((739.9729724794009*nuVtSqSum[0]*dxvl3R2*dxvr3R4+151.7893276880823*nuVtSqSum[0]*dxvl3R3*dxvr3R3+((-66.40783086353598*nuVtSqSum[0]*dxvl3R4)-1897.366596101029*nuVtSqSum[0]*dxvl3R2)*dxvr3R2-47.43416490252571*nuVtSqSum[0]*dxvl3R5*dxvr[3])*fr[14]+(47.43416490252571*nuVtSqSum[0]*dxvl[3]*dxvr3R5+66.40783086353598*nuVtSqSum[0]*dxvl3R2*dxvr3R4-151.7893276880823*nuVtSqSum[0]*dxvl3R3*dxvr3R3+(1897.366596101029*nuVtSqSum[0]*dxvl3R2-739.9729724794009*nuVtSqSum[0]*dxvl3R4)*dxvr3R2)*fl[14]+1697.056274847715*nuVtSqSum[2]*dxvl3R2*dxvr3R2*fr[11]-1697.056274847715*nuVtSqSum[2]*dxvl3R2*dxvr3R2*fl[11]+((-1028.785691968935*nuVtSqSum[1]*dxvl3R2*dxvr3R3)-195.9591794226543*nuVtSqSum[1]*dxvl3R3*dxvr3R2+122.4744871391589*nuVtSqSum[1]*dxvl3R4*dxvr[3]+24.49489742783179*nuVtSqSum[1]*dxvl3R5)*fr[8]+(24.49489742783179*nuVtSqSum[1]*dxvr3R5+122.4744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R4-195.9591794226543*nuVtSqSum[1]*dxvl3R2*dxvr3R3-1028.785691968935*nuVtSqSum[1]*dxvl3R3*dxvr3R2)*fl[8]+((-1028.785691968935*nuVtSqSum[0]*dxvl3R2*dxvr3R3)-195.9591794226543*nuVtSqSum[0]*dxvl3R3*dxvr3R2+122.4744871391589*nuVtSqSum[0]*dxvl3R4*dxvr[3]+24.49489742783179*nuVtSqSum[0]*dxvl3R5)*fr[4]+(24.49489742783179*nuVtSqSum[0]*dxvr3R5+122.4744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R4-195.9591794226543*nuVtSqSum[0]*dxvl3R2*dxvr3R3-1028.785691968935*nuVtSqSum[0]*dxvl3R3*dxvr3R2)*fl[4]+((1697.056274847715*fr[1]-1697.056274847715*fl[1])*nuVtSqSum[1]+(1697.056274847715*fr[0]-1697.056274847715*fl[0])*nuVtSqSum[0])*dxvl3R2*dxvr3R2)/(20.0*dxvr3R5+100.0*dxvl[3]*dxvr3R4+200.0*dxvl3R2*dxvr3R3+200.0*dxvl3R3*dxvr3R2+100.0*dxvl3R4*dxvr[3]+20.0*dxvl3R5)-0.5*(2.23606797749979*fjump[14]+1.732050807568877*fjump[4]+fjump[0])+0.25*(3.16227766016838*alphaDrag[0]*favg[14]+1.414213562373095*alphaDrag[2]*favg[11]+2.449489742783178*alphaDrag[1]*favg[8]+2.449489742783178*alphaDrag[0]*favg[4]+1.414213562373095*alphaDrag[1]*favg[1]+1.414213562373095*alphaDrag[0]*favg[0]); 
  Ghat[1] = ((739.9729724794009*nuVtSqSum[1]*dxvl3R2*dxvr3R4+151.7893276880823*nuVtSqSum[1]*dxvl3R3*dxvr3R3+((-66.40783086353598*nuVtSqSum[1]*dxvl3R4)-1897.366596101029*nuVtSqSum[1]*dxvl3R2)*dxvr3R2-47.43416490252571*nuVtSqSum[1]*dxvl3R5*dxvr[3])*fr[14]+(47.43416490252571*nuVtSqSum[1]*dxvl[3]*dxvr3R5+66.40783086353598*nuVtSqSum[1]*dxvl3R2*dxvr3R4-151.7893276880823*nuVtSqSum[1]*dxvl3R3*dxvr3R3+(1897.366596101029*nuVtSqSum[1]*dxvl3R2-739.9729724794009*nuVtSqSum[1]*dxvl3R4)*dxvr3R2)*fl[14]+1517.893276880824*nuVtSqSum[1]*dxvl3R2*dxvr3R2*fr[11]-1517.893276880824*nuVtSqSum[1]*dxvl3R2*dxvr3R2*fl[11]+(((-920.1738966086795*nuVtSqSum[2])-1028.785691968935*nuVtSqSum[0])*dxvl3R2*dxvr3R3+((-175.2712184016533*nuVtSqSum[2])-195.9591794226543*nuVtSqSum[0])*dxvl3R3*dxvr3R2+(109.5445115010333*nuVtSqSum[2]+122.4744871391589*nuVtSqSum[0])*dxvl3R4*dxvr[3]+(21.90890230020665*nuVtSqSum[2]+24.49489742783179*nuVtSqSum[0])*dxvl3R5)*fr[8]+((21.90890230020665*nuVtSqSum[2]+24.49489742783179*nuVtSqSum[0])*dxvr3R5+(109.5445115010333*nuVtSqSum[2]+122.4744871391589*nuVtSqSum[0])*dxvl[3]*dxvr3R4+((-175.2712184016533*nuVtSqSum[2])-195.9591794226543*nuVtSqSum[0])*dxvl3R2*dxvr3R3+((-920.1738966086795*nuVtSqSum[2])-1028.785691968935*nuVtSqSum[0])*dxvl3R3*dxvr3R2)*fl[8]+((-1028.785691968935*nuVtSqSum[1]*dxvl3R2*dxvr3R3)-195.9591794226543*nuVtSqSum[1]*dxvl3R3*dxvr3R2+122.4744871391589*nuVtSqSum[1]*dxvl3R4*dxvr[3]+24.49489742783179*nuVtSqSum[1]*dxvl3R5)*fr[4]+(24.49489742783179*nuVtSqSum[1]*dxvr3R5+122.4744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R4-195.9591794226543*nuVtSqSum[1]*dxvl3R2*dxvr3R3-1028.785691968935*nuVtSqSum[1]*dxvl3R3*dxvr3R2)*fl[4]+((1517.893276880824*fr[1]-1517.893276880824*fl[1])*nuVtSqSum[2]+(1697.056274847715*fr[0]-1697.056274847715*fl[0])*nuVtSqSum[1]+1697.056274847715*nuVtSqSum[0]*fr[1]-1697.056274847715*nuVtSqSum[0]*fl[1])*dxvl3R2*dxvr3R2)/(20.0*dxvr3R5+100.0*dxvl[3]*dxvr3R4+200.0*dxvl3R2*dxvr3R3+200.0*dxvl3R3*dxvr3R2+100.0*dxvl3R4*dxvr[3]+20.0*dxvl3R5)+0.05*(15.8113883008419*alphaDrag[1]*favg[14]+6.324555320336761*alphaDrag[1]*favg[11]+(10.95445115010332*alphaDrag[2]+12.24744871391589*alphaDrag[0])*favg[8]+12.24744871391589*alphaDrag[1]*favg[4]+6.324555320336761*favg[1]*alphaDrag[2]+7.071067811865476*alphaDrag[0]*favg[1]+7.071067811865476*favg[0]*alphaDrag[1])-0.5*(1.732050807568877*fjump[8]+fjump[1]); 
  Ghat[2] = (-(1.0*((102.8785691968935*nuVtSqSum[0]*dxvl3R2*dxvr3R3+19.59591794226543*nuVtSqSum[0]*dxvl3R3*dxvr3R2-12.24744871391589*nuVtSqSum[0]*dxvl3R4*dxvr[3]-2.449489742783178*nuVtSqSum[0]*dxvl3R5)*fr[9]+((-2.449489742783178*nuVtSqSum[0]*dxvr3R5)-12.24744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R4+19.59591794226543*nuVtSqSum[0]*dxvl3R2*dxvr3R3+102.8785691968935*nuVtSqSum[0]*dxvl3R3*dxvr3R2)*fl[9]-169.7056274847715*nuVtSqSum[1]*dxvl3R2*dxvr3R2*fr[5]+169.7056274847715*nuVtSqSum[1]*dxvl3R2*dxvr3R2*fl[5]+(169.7056274847715*nuVtSqSum[0]*fl[2]-169.7056274847715*nuVtSqSum[0]*fr[2])*dxvl3R2*dxvr3R2))/(2.0*dxvr3R5+10.0*dxvl[3]*dxvr3R4+20.0*dxvl3R2*dxvr3R3+20.0*dxvl3R3*dxvr3R2+10.0*dxvl3R4*dxvr[3]+2.0*dxvl3R5))-0.5*(1.732050807568877*fjump[9]+fjump[2])+0.25*(2.449489742783178*alphaDrag[0]*favg[9]+1.414213562373095*alphaDrag[1]*favg[5]+1.414213562373095*alphaDrag[0]*favg[2]); 
  Ghat[3] = (-(1.0*((102.8785691968935*nuVtSqSum[0]*dxvl3R2*dxvr3R3+19.59591794226543*nuVtSqSum[0]*dxvl3R3*dxvr3R2-12.24744871391589*nuVtSqSum[0]*dxvl3R4*dxvr[3]-2.449489742783178*nuVtSqSum[0]*dxvl3R5)*fr[10]+((-2.449489742783178*nuVtSqSum[0]*dxvr3R5)-12.24744871391589*nuVtSqSum[0]*dxvl[3]*dxvr3R4+19.59591794226543*nuVtSqSum[0]*dxvl3R2*dxvr3R3+102.8785691968935*nuVtSqSum[0]*dxvl3R3*dxvr3R2)*fl[10]-169.7056274847715*nuVtSqSum[1]*dxvl3R2*dxvr3R2*fr[6]+169.7056274847715*nuVtSqSum[1]*dxvl3R2*dxvr3R2*fl[6]-169.7056274847715*nuVtSqSum[0]*dxvl3R2*dxvr3R2*fr[3]+169.7056274847715*nuVtSqSum[0]*dxvl3R2*dxvr3R2*fl[3]))/(2.0*dxvr3R5+10.0*dxvl[3]*dxvr3R4+20.0*dxvl3R2*dxvr3R3+20.0*dxvl3R3*dxvr3R2+10.0*dxvl3R4*dxvr[3]+2.0*dxvl3R5))-0.5*(1.732050807568877*fjump[10]+fjump[3])+0.25*(2.449489742783178*alphaDrag[0]*favg[10]+1.414213562373095*alphaDrag[1]*favg[6]+1.414213562373095*alphaDrag[0]*favg[3]); 
  Ghat[5] = (-(1.0*((102.8785691968935*nuVtSqSum[1]*dxvl3R2*dxvr3R3+19.59591794226543*nuVtSqSum[1]*dxvl3R3*dxvr3R2-12.24744871391589*nuVtSqSum[1]*dxvl3R4*dxvr[3]-2.449489742783178*nuVtSqSum[1]*dxvl3R5)*fr[9]+((-2.449489742783178*nuVtSqSum[1]*dxvr3R5)-12.24744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R4+19.59591794226543*nuVtSqSum[1]*dxvl3R2*dxvr3R3+102.8785691968935*nuVtSqSum[1]*dxvl3R3*dxvr3R2)*fl[9]+((-151.7893276880823*nuVtSqSum[2])-169.7056274847715*nuVtSqSum[0])*dxvl3R2*dxvr3R2*fr[5]+(151.7893276880823*nuVtSqSum[2]+169.7056274847715*nuVtSqSum[0])*dxvl3R2*dxvr3R2*fl[5]+(169.7056274847715*nuVtSqSum[1]*fl[2]-169.7056274847715*nuVtSqSum[1]*fr[2])*dxvl3R2*dxvr3R2))/(2.0*dxvr3R5+10.0*dxvl[3]*dxvr3R4+20.0*dxvl3R2*dxvr3R3+20.0*dxvl3R3*dxvr3R2+10.0*dxvl3R4*dxvr[3]+2.0*dxvl3R5))+0.05*(12.24744871391589*alphaDrag[1]*favg[9]+(6.324555320336761*alphaDrag[2]+7.071067811865476*alphaDrag[0])*favg[5]+7.071067811865476*alphaDrag[1]*favg[2])-0.5*fjump[5]; 
  Ghat[6] = (-(1.0*((102.8785691968935*nuVtSqSum[1]*dxvl3R2*dxvr3R3+19.59591794226543*nuVtSqSum[1]*dxvl3R3*dxvr3R2-12.24744871391589*nuVtSqSum[1]*dxvl3R4*dxvr[3]-2.449489742783178*nuVtSqSum[1]*dxvl3R5)*fr[10]+((-2.449489742783178*nuVtSqSum[1]*dxvr3R5)-12.24744871391589*nuVtSqSum[1]*dxvl[3]*dxvr3R4+19.59591794226543*nuVtSqSum[1]*dxvl3R2*dxvr3R3+102.8785691968935*nuVtSqSum[1]*dxvl3R3*dxvr3R2)*fl[10]+((-151.7893276880823*nuVtSqSum[2])-169.7056274847715*nuVtSqSum[0])*dxvl3R2*dxvr3R2*fr[6]+(151.7893276880823*nuVtSqSum[2]+169.7056274847715*nuVtSqSum[0])*dxvl3R2*dxvr3R2*fl[6]-169.7056274847715*nuVtSqSum[1]*dxvl3R2*dxvr3R2*fr[3]+169.7056274847715*nuVtSqSum[1]*dxvl3R2*dxvr3R2*fl[3]))/(2.0*dxvr3R5+10.0*dxvl[3]*dxvr3R4+20.0*dxvl3R2*dxvr3R3+20.0*dxvl3R3*dxvr3R2+10.0*dxvl3R4*dxvr[3]+2.0*dxvl3R5))+0.05*(12.24744871391589*alphaDrag[1]*favg[10]+(6.324555320336761*alphaDrag[2]+7.071067811865476*alphaDrag[0])*favg[6]+7.071067811865476*alphaDrag[1]*favg[3])-0.5*fjump[6]; 
  Ghat[7] = (84.85281374238573*nuVtSqSum[0]*dxvl3R2*dxvr3R2*fr[7]-84.85281374238573*nuVtSqSum[0]*dxvl3R2*dxvr3R2*fl[7])/(dxvr3R5+5.0*dxvl[3]*dxvr3R4+10.0*dxvl3R2*dxvr3R3+10.0*dxvl3R3*dxvr3R2+5.0*dxvl3R4*dxvr[3]+dxvl3R5)-0.5*fjump[7]+0.3535533905932737*alphaDrag[0]*favg[7]; 
  Ghat[11] = ((5179.810807355807*nuVtSqSum[2]*dxvl3R2*dxvr3R4+1062.525293816576*nuVtSqSum[2]*dxvl3R3*dxvr3R3+((-464.8548160447518*nuVtSqSum[2]*dxvl3R4)-13281.5661727072*nuVtSqSum[2]*dxvl3R2)*dxvr3R2-332.0391543176799*nuVtSqSum[2]*dxvl3R5*dxvr[3])*fr[14]+(332.0391543176799*nuVtSqSum[2]*dxvl[3]*dxvr3R5+464.8548160447518*nuVtSqSum[2]*dxvl3R2*dxvr3R4-1062.525293816576*nuVtSqSum[2]*dxvl3R3*dxvr3R3+(13281.5661727072*nuVtSqSum[2]*dxvl3R2-5179.810807355807*nuVtSqSum[2]*dxvl3R4)*dxvr3R2)*fl[14]+(7589.466384404118*nuVtSqSum[2]+11879.39392393401*nuVtSqSum[0])*dxvl3R2*dxvr3R2*fr[11]+((-7589.466384404118*nuVtSqSum[2])-11879.39392393401*nuVtSqSum[0])*dxvl3R2*dxvr3R2*fl[11]+((-6441.217276260756*nuVtSqSum[1]*dxvl3R2*dxvr3R3)-1226.898528811573*nuVtSqSum[1]*dxvl3R3*dxvr3R2+766.811580507233*nuVtSqSum[1]*dxvl3R4*dxvr[3]+153.3623161014466*nuVtSqSum[1]*dxvl3R5)*fr[8]+(153.3623161014466*nuVtSqSum[1]*dxvr3R5+766.811580507233*nuVtSqSum[1]*dxvl[3]*dxvr3R4-1226.898528811573*nuVtSqSum[1]*dxvl3R2*dxvr3R3-6441.217276260756*nuVtSqSum[1]*dxvl3R3*dxvr3R2)*fl[8]+((-7201.499843782545*nuVtSqSum[2]*dxvl3R2*dxvr3R3)-1371.71425595858*nuVtSqSum[2]*dxvl3R3*dxvr3R2+857.3214099741125*nuVtSqSum[2]*dxvl3R4*dxvr[3]+171.4642819948225*nuVtSqSum[2]*dxvl3R5)*fr[4]+(171.4642819948225*nuVtSqSum[2]*dxvr3R5+857.3214099741125*nuVtSqSum[2]*dxvl[3]*dxvr3R4-1371.71425595858*nuVtSqSum[2]*dxvl3R2*dxvr3R3-7201.499843782545*nuVtSqSum[2]*dxvl3R3*dxvr3R2)*fl[4]+((11879.39392393401*fr[0]-11879.39392393401*fl[0])*nuVtSqSum[2]+(10625.25293816576*fr[1]-10625.25293816576*fl[1])*nuVtSqSum[1])*dxvl3R2*dxvr3R2)/(140.0*dxvr3R5+700.0*dxvl[3]*dxvr3R4+1400.0*dxvl3R2*dxvr3R3+1400.0*dxvl3R3*dxvr3R2+700.0*dxvl3R4*dxvr[3]+140.0*dxvl3R5)+0.007142857142857143*(110.6797181058933*alphaDrag[2]*favg[14]+(31.62277660168381*alphaDrag[2]+49.49747468305833*alphaDrag[0])*favg[11]+76.68115805072327*alphaDrag[1]*favg[8]+85.73214099741124*alphaDrag[2]*favg[4]+49.49747468305833*favg[0]*alphaDrag[2]+44.27188724235732*alphaDrag[1]*favg[1])-0.5*fjump[11]; 
  Ghat[12] = (84.85281374238573*nuVtSqSum[0]*dxvl3R2*dxvr3R2*fr[12]-84.85281374238573*nuVtSqSum[0]*dxvl3R2*dxvr3R2*fl[12])/(dxvr3R5+5.0*dxvl[3]*dxvr3R4+10.0*dxvl3R2*dxvr3R3+10.0*dxvl3R3*dxvr3R2+5.0*dxvl3R4*dxvr[3]+dxvl3R5)-0.5*fjump[12]+0.3535533905932737*alphaDrag[0]*favg[12]; 
  Ghat[13] = (84.85281374238573*nuVtSqSum[0]*dxvl3R2*dxvr3R2*fr[13]-84.85281374238573*nuVtSqSum[0]*dxvl3R2*dxvr3R2*fl[13])/(dxvr3R5+5.0*dxvl[3]*dxvr3R4+10.0*dxvl3R2*dxvr3R3+10.0*dxvl3R3*dxvr3R2+5.0*dxvl3R4*dxvr[3]+dxvl3R5)-0.5*fjump[13]+0.3535533905932737*alphaDrag[0]*favg[13]; 

  double incr1[15]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[0]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = 0.8660254037844386*Ghat[1]; 
  incr1[9] = 0.8660254037844386*Ghat[2]; 
  incr1[10] = 0.8660254037844386*Ghat[3]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = -0.5*Ghat[13]; 
  incr1[14] = -1.118033988749895*Ghat[0]; 

  double incr2[15]; 

  incr2[4] = ((427.2235948540296*nuVtSqSum[0]*dxvl3R3*dxvr3R4+312.2018577779447*nuVtSqSum[0]*dxvl3R4*dxvr3R3+(82.15838362577493*nuVtSqSum[0]*dxvl3R5-1095.445115010333*nuVtSqSum[0]*dxvl3R3)*dxvr3R2-547.7225575051664*nuVtSqSum[0]*dxvl3R4*dxvr[3]-109.5445115010333*nuVtSqSum[0]*dxvl3R5)*fr[14]+((82.15838362577493*nuVtSqSum[0]*dxvl3R2-109.5445115010333*nuVtSqSum[0])*dxvr3R5+(312.2018577779447*nuVtSqSum[0]*dxvl3R3-547.7225575051664*nuVtSqSum[0]*dxvl[3])*dxvr3R4+(427.2235948540296*nuVtSqSum[0]*dxvl3R4-1095.445115010333*nuVtSqSum[0]*dxvl3R2)*dxvr3R3)*fl[14]+(979.7958971132716*nuVtSqSum[2]*dxvl3R3*dxvr3R2+489.8979485566358*nuVtSqSum[2]*dxvl3R4*dxvr[3]+97.97958971132716*nuVtSqSum[2]*dxvl3R5)*fr[11]+(97.97958971132716*nuVtSqSum[2]*dxvr3R5+489.8979485566358*nuVtSqSum[2]*dxvl[3]*dxvr3R4+979.7958971132716*nuVtSqSum[2]*dxvl3R2*dxvr3R3)*fl[11]+((-593.9696961967002*nuVtSqSum[1]*dxvl3R3*dxvr3R3)-424.2640687119287*nuVtSqSum[1]*dxvl3R4*dxvr3R2-84.85281374238573*nuVtSqSum[1]*dxvl3R5*dxvr[3])*fr[8]+(84.85281374238573*nuVtSqSum[1]*dxvl[3]*dxvr3R5+424.2640687119287*nuVtSqSum[1]*dxvl3R2*dxvr3R4+593.9696961967002*nuVtSqSum[1]*dxvl3R3*dxvr3R3)*fl[8]+((-593.9696961967002*nuVtSqSum[0]*dxvl3R3*dxvr3R3)-424.2640687119287*nuVtSqSum[0]*dxvl3R4*dxvr3R2-84.85281374238573*nuVtSqSum[0]*dxvl3R5*dxvr[3])*fr[4]+(84.85281374238573*nuVtSqSum[0]*dxvl[3]*dxvr3R5+424.2640687119287*nuVtSqSum[0]*dxvl3R2*dxvr3R4+593.9696961967002*nuVtSqSum[0]*dxvl3R3*dxvr3R3)*fl[4]+(97.97958971132716*fl[1]*nuVtSqSum[1]+97.97958971132716*fl[0]*nuVtSqSum[0])*dxvr3R5+(489.8979485566358*fl[1]*nuVtSqSum[1]+489.8979485566358*fl[0]*nuVtSqSum[0])*dxvl[3]*dxvr3R4+(979.7958971132716*fl[1]*nuVtSqSum[1]+979.7958971132716*fl[0]*nuVtSqSum[0])*dxvl3R2*dxvr3R3+(979.7958971132716*fr[1]*nuVtSqSum[1]+979.7958971132716*fr[0]*nuVtSqSum[0])*dxvl3R3*dxvr3R2+(489.8979485566358*fr[1]*nuVtSqSum[1]+489.8979485566358*fr[0]*nuVtSqSum[0])*dxvl3R4*dxvr[3]+(97.97958971132716*fr[1]*nuVtSqSum[1]+97.97958971132716*fr[0]*nuVtSqSum[0])*dxvl3R5)/(160.0*dxvr3R5+800.0*dxvl[3]*dxvr3R4+1600.0*dxvl3R2*dxvr3R3+1600.0*dxvl3R3*dxvr3R2+800.0*dxvl3R4*dxvr[3]+160.0*dxvl3R5); 
  incr2[8] = ((427.2235948540296*nuVtSqSum[1]*dxvl3R3*dxvr3R4+312.2018577779447*nuVtSqSum[1]*dxvl3R4*dxvr3R3+(82.15838362577493*nuVtSqSum[1]*dxvl3R5-1095.445115010333*nuVtSqSum[1]*dxvl3R3)*dxvr3R2-547.7225575051664*nuVtSqSum[1]*dxvl3R4*dxvr[3]-109.5445115010333*nuVtSqSum[1]*dxvl3R5)*fr[14]+((82.15838362577493*nuVtSqSum[1]*dxvl3R2-109.5445115010333*nuVtSqSum[1])*dxvr3R5+(312.2018577779447*nuVtSqSum[1]*dxvl3R3-547.7225575051664*nuVtSqSum[1]*dxvl[3])*dxvr3R4+(427.2235948540296*nuVtSqSum[1]*dxvl3R4-1095.445115010333*nuVtSqSum[1]*dxvl3R2)*dxvr3R3)*fl[14]+(876.3560920082665*nuVtSqSum[1]*dxvl3R3*dxvr3R2+438.1780460041332*nuVtSqSum[1]*dxvl3R4*dxvr[3]+87.63560920082662*nuVtSqSum[1]*dxvl3R5)*fr[11]+(87.63560920082662*nuVtSqSum[1]*dxvr3R5+438.1780460041332*nuVtSqSum[1]*dxvl[3]*dxvr3R4+876.3560920082665*nuVtSqSum[1]*dxvl3R2*dxvr3R3)*fl[11]+(((-531.262646908288*nuVtSqSum[2])-593.9696961967002*nuVtSqSum[0])*dxvl3R3*dxvr3R3+((-379.4733192202058*nuVtSqSum[2])-424.2640687119287*nuVtSqSum[0])*dxvl3R4*dxvr3R2+((-75.89466384404115*nuVtSqSum[2])-84.85281374238573*nuVtSqSum[0])*dxvl3R5*dxvr[3])*fr[8]+((75.89466384404115*nuVtSqSum[2]+84.85281374238573*nuVtSqSum[0])*dxvl[3]*dxvr3R5+(379.4733192202058*nuVtSqSum[2]+424.2640687119287*nuVtSqSum[0])*dxvl3R2*dxvr3R4+(531.262646908288*nuVtSqSum[2]+593.9696961967002*nuVtSqSum[0])*dxvl3R3*dxvr3R3)*fl[8]+((-593.9696961967002*nuVtSqSum[1]*dxvl3R3*dxvr3R3)-424.2640687119287*nuVtSqSum[1]*dxvl3R4*dxvr3R2-84.85281374238573*nuVtSqSum[1]*dxvl3R5*dxvr[3])*fr[4]+(84.85281374238573*nuVtSqSum[1]*dxvl[3]*dxvr3R5+424.2640687119287*nuVtSqSum[1]*dxvl3R2*dxvr3R4+593.9696961967002*nuVtSqSum[1]*dxvl3R3*dxvr3R3)*fl[4]+(87.63560920082662*fl[1]*nuVtSqSum[2]+97.97958971132716*fl[0]*nuVtSqSum[1]+97.97958971132716*nuVtSqSum[0]*fl[1])*dxvr3R5+(438.1780460041332*fl[1]*nuVtSqSum[2]+489.8979485566358*fl[0]*nuVtSqSum[1]+489.8979485566358*nuVtSqSum[0]*fl[1])*dxvl[3]*dxvr3R4+(876.3560920082665*fl[1]*nuVtSqSum[2]+979.7958971132716*fl[0]*nuVtSqSum[1]+979.7958971132716*nuVtSqSum[0]*fl[1])*dxvl3R2*dxvr3R3+(876.3560920082665*fr[1]*nuVtSqSum[2]+979.7958971132716*fr[0]*nuVtSqSum[1]+979.7958971132716*nuVtSqSum[0]*fr[1])*dxvl3R3*dxvr3R2+(438.1780460041332*fr[1]*nuVtSqSum[2]+489.8979485566358*fr[0]*nuVtSqSum[1]+489.8979485566358*nuVtSqSum[0]*fr[1])*dxvl3R4*dxvr[3]+(87.63560920082662*fr[1]*nuVtSqSum[2]+97.97958971132716*fr[0]*nuVtSqSum[1]+97.97958971132716*nuVtSqSum[0]*fr[1])*dxvl3R5)/(160.0*dxvr3R5+800.0*dxvl[3]*dxvr3R4+1600.0*dxvl3R2*dxvr3R3+1600.0*dxvl3R3*dxvr3R2+800.0*dxvl3R4*dxvr[3]+160.0*dxvl3R5); 
  incr2[9] = -(1.0*((29.698484809835*nuVtSqSum[0]*dxvl3R3*dxvr3R3+21.21320343559643*nuVtSqSum[0]*dxvl3R4*dxvr3R2+4.242640687119286*nuVtSqSum[0]*dxvl3R5*dxvr[3])*fr[9]+((-4.242640687119286*nuVtSqSum[0]*dxvl[3]*dxvr3R5)-21.21320343559643*nuVtSqSum[0]*dxvl3R2*dxvr3R4-29.698484809835*nuVtSqSum[0]*dxvl3R3*dxvr3R3)*fl[9]+((-48.98979485566358*nuVtSqSum[1]*dxvl3R3*dxvr3R2)-24.49489742783179*nuVtSqSum[1]*dxvl3R4*dxvr[3]-4.898979485566357*nuVtSqSum[1]*dxvl3R5)*fr[5]+((-4.898979485566357*nuVtSqSum[1]*dxvr3R5)-24.49489742783179*nuVtSqSum[1]*dxvl[3]*dxvr3R4-48.98979485566358*nuVtSqSum[1]*dxvl3R2*dxvr3R3)*fl[5]-4.898979485566357*nuVtSqSum[0]*fl[2]*dxvr3R5-24.49489742783179*nuVtSqSum[0]*fl[2]*dxvl[3]*dxvr3R4-48.98979485566358*nuVtSqSum[0]*fl[2]*dxvl3R2*dxvr3R3-48.98979485566358*nuVtSqSum[0]*fr[2]*dxvl3R3*dxvr3R2-24.49489742783179*nuVtSqSum[0]*fr[2]*dxvl3R4*dxvr[3]-4.898979485566357*nuVtSqSum[0]*fr[2]*dxvl3R5))/(8.0*dxvr3R5+40.0*dxvl[3]*dxvr3R4+80.0*dxvl3R2*dxvr3R3+80.0*dxvl3R3*dxvr3R2+40.0*dxvl3R4*dxvr[3]+8.0*dxvl3R5); 
  incr2[10] = -(1.0*((29.698484809835*nuVtSqSum[0]*dxvl3R3*dxvr3R3+21.21320343559643*nuVtSqSum[0]*dxvl3R4*dxvr3R2+4.242640687119286*nuVtSqSum[0]*dxvl3R5*dxvr[3])*fr[10]+((-4.242640687119286*nuVtSqSum[0]*dxvl[3]*dxvr3R5)-21.21320343559643*nuVtSqSum[0]*dxvl3R2*dxvr3R4-29.698484809835*nuVtSqSum[0]*dxvl3R3*dxvr3R3)*fl[10]+((-48.98979485566358*nuVtSqSum[1]*dxvl3R3*dxvr3R2)-24.49489742783179*nuVtSqSum[1]*dxvl3R4*dxvr[3]-4.898979485566357*nuVtSqSum[1]*dxvl3R5)*fr[6]+((-4.898979485566357*nuVtSqSum[1]*dxvr3R5)-24.49489742783179*nuVtSqSum[1]*dxvl[3]*dxvr3R4-48.98979485566358*nuVtSqSum[1]*dxvl3R2*dxvr3R3)*fl[6]+((-48.98979485566358*nuVtSqSum[0]*dxvl3R3*dxvr3R2)-24.49489742783179*nuVtSqSum[0]*dxvl3R4*dxvr[3]-4.898979485566357*nuVtSqSum[0]*dxvl3R5)*fr[3]+((-4.898979485566357*nuVtSqSum[0]*dxvr3R5)-24.49489742783179*nuVtSqSum[0]*dxvl[3]*dxvr3R4-48.98979485566358*nuVtSqSum[0]*dxvl3R2*dxvr3R3)*fl[3]))/(8.0*dxvr3R5+40.0*dxvl[3]*dxvr3R4+80.0*dxvl3R2*dxvr3R3+80.0*dxvl3R3*dxvr3R2+40.0*dxvl3R4*dxvr[3]+8.0*dxvl3R5); 
  incr2[14] = -(1.0*((330.9259735953043*nuVtSqSum[0]*dxvl3R3*dxvr3R4+241.8305191657993*nuVtSqSum[0]*dxvl3R4*dxvr3R3+(63.63961030678928*nuVtSqSum[0]*dxvl3R5-848.5281374238575*nuVtSqSum[0]*dxvl3R3)*dxvr3R2-424.2640687119287*nuVtSqSum[0]*dxvl3R4*dxvr[3]-84.85281374238573*nuVtSqSum[0]*dxvl3R5)*fr[14]+((63.63961030678928*nuVtSqSum[0]*dxvl3R2-84.85281374238573*nuVtSqSum[0])*dxvr3R5+(241.8305191657993*nuVtSqSum[0]*dxvl3R3-424.2640687119287*nuVtSqSum[0]*dxvl[3])*dxvr3R4+(330.9259735953043*nuVtSqSum[0]*dxvl3R4-848.5281374238575*nuVtSqSum[0]*dxvl3R2)*dxvr3R3)*fl[14]+(758.9466384404116*nuVtSqSum[2]*dxvl3R3*dxvr3R2+379.4733192202058*nuVtSqSum[2]*dxvl3R4*dxvr[3]+75.89466384404115*nuVtSqSum[2]*dxvl3R5)*fr[11]+(75.89466384404115*nuVtSqSum[2]*dxvr3R5+379.4733192202058*nuVtSqSum[2]*dxvl[3]*dxvr3R4+758.9466384404116*nuVtSqSum[2]*dxvl3R2*dxvr3R3)*fl[11]+((-460.0869483043397*nuVtSqSum[1]*dxvl3R3*dxvr3R3)-328.6335345030998*nuVtSqSum[1]*dxvl3R4*dxvr3R2-65.72670690061996*nuVtSqSum[1]*dxvl3R5*dxvr[3])*fr[8]+(65.72670690061996*nuVtSqSum[1]*dxvl[3]*dxvr3R5+328.6335345030998*nuVtSqSum[1]*dxvl3R2*dxvr3R4+460.0869483043397*nuVtSqSum[1]*dxvl3R3*dxvr3R3)*fl[8]+((-460.0869483043397*nuVtSqSum[0]*dxvl3R3*dxvr3R3)-328.6335345030998*nuVtSqSum[0]*dxvl3R4*dxvr3R2-65.72670690061996*nuVtSqSum[0]*dxvl3R5*dxvr[3])*fr[4]+(65.72670690061996*nuVtSqSum[0]*dxvl[3]*dxvr3R5+328.6335345030998*nuVtSqSum[0]*dxvl3R2*dxvr3R4+460.0869483043397*nuVtSqSum[0]*dxvl3R3*dxvr3R3)*fl[4]+(75.89466384404115*fl[1]*nuVtSqSum[1]+75.89466384404115*fl[0]*nuVtSqSum[0])*dxvr3R5+(379.4733192202058*fl[1]*nuVtSqSum[1]+379.4733192202058*fl[0]*nuVtSqSum[0])*dxvl[3]*dxvr3R4+(758.9466384404116*fl[1]*nuVtSqSum[1]+758.9466384404116*fl[0]*nuVtSqSum[0])*dxvl3R2*dxvr3R3+(758.9466384404116*fr[1]*nuVtSqSum[1]+758.9466384404116*fr[0]*nuVtSqSum[0])*dxvl3R3*dxvr3R2+(379.4733192202058*fr[1]*nuVtSqSum[1]+379.4733192202058*fr[0]*nuVtSqSum[0])*dxvl3R4*dxvr[3]+(75.89466384404115*fr[1]*nuVtSqSum[1]+75.89466384404115*fr[0]*nuVtSqSum[0])*dxvl3R5))/(32.0*dxvr3R5+160.0*dxvl[3]*dxvr3R4+320.0*dxvl3R2*dxvr3R3+320.0*dxvl3R3*dxvr3R2+160.0*dxvl3R4*dxvr[3]+32.0*dxvl3R5); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr2[8]*rdvSq4R+incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += incr1[8]*rdv2L-1.0*incr2[8]*rdvSq4L; 
  outl[9] += incr1[9]*rdv2L-1.0*incr2[9]*rdvSq4L; 
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += -1.0*incr1[13]*rdv2L; 
  outl[14] += incr2[14]*rdvSq4L-1.0*incr1[14]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUz[2])/nuSum-(0.7071067811865475*sumNuUz[0])/nuSum+wl[3]); 
} 
