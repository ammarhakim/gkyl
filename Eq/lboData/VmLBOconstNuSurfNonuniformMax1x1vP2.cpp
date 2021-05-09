#include <VmLBOModDecl.h> 
double VmLBOconstNuSurfNonUniform1x1vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
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

  double favg[6]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 

  double fjump[6]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 

  double Ghat[6]; 
  for(unsigned short int i=0; i<6; ++i){ 
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

  Ghat[0] = ((265.631323454144*nuVtSqSum[0]*dxvl1R3*dxvr1R3+151.7893276880823*nuVtSqSum[0]*dxvl1R4*dxvr1R2-66.40783086353598*nuVtSqSum[0]*dxvl1R5*dxvr[1]-47.43416490252571*nuVtSqSum[0]*dxvl1R6)*fr[5]+(47.43416490252571*nuVtSqSum[0]*dxvr1R6+66.40783086353598*nuVtSqSum[0]*dxvl[1]*dxvr1R5-151.7893276880823*nuVtSqSum[0]*dxvl1R2*dxvr1R4-265.631323454144*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[5]+424.2640687119287*dxvl1R3*dxvr1R3*nuVtSqSum[2]*fr[4]-424.2640687119287*dxvl1R3*dxvr1R3*nuVtSqSum[2]*fl[4]+((-514.3928459844674*dxvl1R3*dxvr1R3)-97.97958971132716*dxvl1R4*dxvr1R2+61.23724356957945*dxvl1R5*dxvr[1]+12.24744871391589*dxvl1R6)*nuVtSqSum[1]*fr[3]+(12.24744871391589*dxvr1R6+61.23724356957945*dxvl[1]*dxvr1R5-97.97958971132716*dxvl1R2*dxvr1R4-514.3928459844674*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[3]+((-514.3928459844674*nuVtSqSum[0]*dxvl1R3*dxvr1R3)-97.97958971132716*nuVtSqSum[0]*dxvl1R4*dxvr1R2+61.23724356957945*nuVtSqSum[0]*dxvl1R5*dxvr[1]+12.24744871391589*nuVtSqSum[0]*dxvl1R6)*fr[2]+(12.24744871391589*nuVtSqSum[0]*dxvr1R6+61.23724356957945*nuVtSqSum[0]*dxvl[1]*dxvr1R5-97.97958971132716*nuVtSqSum[0]*dxvl1R2*dxvr1R4-514.3928459844674*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[2]+(424.2640687119287*dxvl1R3*dxvr1R3*fr[1]-424.2640687119287*dxvl1R3*dxvr1R3*fl[1])*nuVtSqSum[1]+(424.2640687119287*fr[0]-424.2640687119287*fl[0])*nuVtSqSum[0]*dxvl1R3*dxvr1R3)/(5.0*dxvl[1]*dxvr1R6+25.0*dxvl1R2*dxvr1R5+50.0*dxvl1R3*dxvr1R4+50.0*dxvl1R4*dxvr1R3+25.0*dxvl1R5*dxvr1R2+5.0*dxvl1R6*dxvr[1])-0.5*(2.23606797749979*fjump[5]+1.732050807568877*fjump[2]+fjump[0])+0.25*(3.16227766016838*alphaDrag[0]*favg[5]+1.414213562373095*alphaDrag[2]*favg[4]+2.449489742783178*alphaDrag[1]*favg[3]+2.449489742783178*alphaDrag[0]*favg[2]+1.414213562373095*alphaDrag[1]*favg[1]+1.414213562373095*alphaDrag[0]*favg[0]); 
  Ghat[1] = ((265.631323454144*dxvl1R3*dxvr1R3+151.7893276880823*dxvl1R4*dxvr1R2-66.40783086353598*dxvl1R5*dxvr[1]-47.43416490252571*dxvl1R6)*nuVtSqSum[1]*fr[5]+(47.43416490252571*dxvr1R6+66.40783086353598*dxvl[1]*dxvr1R5-151.7893276880823*dxvl1R2*dxvr1R4-265.631323454144*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[5]+379.4733192202058*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fr[4]-379.4733192202058*dxvl1R3*dxvr1R3*nuVtSqSum[1]*fl[4]+(((-460.0869483043397*dxvl1R3*dxvr1R3)-87.63560920082662*dxvl1R4*dxvr1R2+54.77225575051663*dxvl1R5*dxvr[1]+10.95445115010332*dxvl1R6)*nuVtSqSum[2]-514.3928459844674*nuVtSqSum[0]*dxvl1R3*dxvr1R3-97.97958971132716*nuVtSqSum[0]*dxvl1R4*dxvr1R2+61.23724356957945*nuVtSqSum[0]*dxvl1R5*dxvr[1]+12.24744871391589*nuVtSqSum[0]*dxvl1R6)*fr[3]+((10.95445115010332*dxvr1R6+54.77225575051663*dxvl[1]*dxvr1R5-87.63560920082662*dxvl1R2*dxvr1R4-460.0869483043397*dxvl1R3*dxvr1R3)*nuVtSqSum[2]+12.24744871391589*nuVtSqSum[0]*dxvr1R6+61.23724356957945*nuVtSqSum[0]*dxvl[1]*dxvr1R5-97.97958971132716*nuVtSqSum[0]*dxvl1R2*dxvr1R4-514.3928459844674*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[3]+(379.4733192202058*dxvl1R3*dxvr1R3*fr[1]-379.4733192202058*dxvl1R3*dxvr1R3*fl[1])*nuVtSqSum[2]+((-514.3928459844674*dxvl1R3*dxvr1R3)-97.97958971132716*dxvl1R4*dxvr1R2+61.23724356957945*dxvl1R5*dxvr[1]+12.24744871391589*dxvl1R6)*nuVtSqSum[1]*fr[2]+(12.24744871391589*dxvr1R6+61.23724356957945*dxvl[1]*dxvr1R5-97.97958971132716*dxvl1R2*dxvr1R4-514.3928459844674*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[2]+(424.2640687119287*fr[0]-424.2640687119287*fl[0])*dxvl1R3*dxvr1R3*nuVtSqSum[1]+424.2640687119287*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fr[1]-424.2640687119287*nuVtSqSum[0]*dxvl1R3*dxvr1R3*fl[1])/(5.0*dxvl[1]*dxvr1R6+25.0*dxvl1R2*dxvr1R5+50.0*dxvl1R3*dxvr1R4+50.0*dxvl1R4*dxvr1R3+25.0*dxvl1R5*dxvr1R2+5.0*dxvl1R6*dxvr[1])+0.05*(15.8113883008419*alphaDrag[1]*favg[5]+6.324555320336761*alphaDrag[1]*favg[4]+(10.95445115010332*alphaDrag[2]+12.24744871391589*alphaDrag[0])*favg[3]+12.24744871391589*alphaDrag[1]*favg[2]+6.324555320336761*favg[1]*alphaDrag[2]+7.071067811865476*alphaDrag[0]*favg[1]+7.071067811865476*favg[0]*alphaDrag[1])-0.5*(1.732050807568877*fjump[3]+fjump[1]); 
  Ghat[4] = ((1859.419264179008*dxvl1R3*dxvr1R3+1062.525293816576*dxvl1R4*dxvr1R2-464.8548160447518*dxvl1R5*dxvr[1]-332.0391543176799*dxvl1R6)*nuVtSqSum[2]*fr[5]+(332.0391543176799*dxvr1R6+464.8548160447518*dxvl[1]*dxvr1R5-1062.525293816576*dxvl1R2*dxvr1R4-1859.419264179008*dxvl1R3*dxvr1R3)*nuVtSqSum[2]*fl[5]+(1897.366596101029*dxvl1R3*dxvr1R3*nuVtSqSum[2]+2969.848480983501*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fr[4]+((-1897.366596101029*dxvl1R3*dxvr1R3*nuVtSqSum[2])-2969.848480983501*nuVtSqSum[0]*dxvl1R3*dxvr1R3)*fl[4]+((-3220.608638130377*dxvl1R3*dxvr1R3)-613.4492644057864*dxvl1R4*dxvr1R2+383.4057902536164*dxvl1R5*dxvr[1]+76.68115805072327*dxvl1R6)*nuVtSqSum[1]*fr[3]+(76.68115805072327*dxvr1R6+383.4057902536164*dxvl[1]*dxvr1R5-613.4492644057864*dxvl1R2*dxvr1R4-3220.608638130377*dxvl1R3*dxvr1R3)*nuVtSqSum[1]*fl[3]+(((-3600.749921891272*dxvl1R3*dxvr1R3)-685.8571279792902*dxvl1R4*dxvr1R2+428.6607049870561*dxvl1R5*dxvr[1]+85.73214099741124*dxvl1R6)*fr[2]+(85.73214099741124*dxvr1R6+428.6607049870561*dxvl[1]*dxvr1R5-685.8571279792902*dxvl1R2*dxvr1R4-3600.749921891272*dxvl1R3*dxvr1R3)*fl[2]+(2969.848480983501*fr[0]-2969.848480983501*fl[0])*dxvl1R3*dxvr1R3)*nuVtSqSum[2]+(2656.313234541441*dxvl1R3*dxvr1R3*fr[1]-2656.313234541441*dxvl1R3*dxvr1R3*fl[1])*nuVtSqSum[1])/(35.0*dxvl[1]*dxvr1R6+175.0*dxvl1R2*dxvr1R5+350.0*dxvl1R3*dxvr1R4+350.0*dxvl1R4*dxvr1R3+175.0*dxvl1R5*dxvr1R2+35.0*dxvl1R6*dxvr[1])+0.007142857142857143*(110.6797181058933*alphaDrag[2]*favg[5]+(31.62277660168381*alphaDrag[2]+49.49747468305833*alphaDrag[0])*favg[4]+76.68115805072327*alphaDrag[1]*favg[3]+85.73214099741124*alphaDrag[2]*favg[2]+49.49747468305833*favg[0]*alphaDrag[2]+44.27188724235732*alphaDrag[1]*favg[1])-0.5*fjump[4]; 

  double incr1[6]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -1.118033988749895*Ghat[0]; 

  double incr2[6]; 

  incr2[2] = ((76.68115805072327*nuVtSqSum[0]*dxvl1R3*dxvr1R2+87.63560920082662*nuVtSqSum[0]*dxvl1R4*dxvr[1]+27.38612787525831*nuVtSqSum[0]*dxvl1R5)*fr[5]+(27.38612787525831*nuVtSqSum[0]*dxvr1R5+87.63560920082662*nuVtSqSum[0]*dxvl[1]*dxvr1R4+76.68115805072327*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[5]+(122.4744871391589*dxvl1R3*dxvr1R2+61.23724356957945*dxvl1R4*dxvr[1]+12.24744871391589*dxvl1R5)*nuVtSqSum[2]*fr[4]+(12.24744871391589*dxvr1R5+61.23724356957945*dxvl[1]*dxvr1R4+122.4744871391589*dxvl1R2*dxvr1R3)*nuVtSqSum[2]*fl[4]+((-148.492424049175*dxvl1R3*dxvr1R2)-106.0660171779821*dxvl1R4*dxvr[1]-21.21320343559643*dxvl1R5)*nuVtSqSum[1]*fr[3]+(21.21320343559643*dxvr1R5+106.0660171779821*dxvl[1]*dxvr1R4+148.492424049175*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[3]+((-148.492424049175*nuVtSqSum[0]*dxvl1R3*dxvr1R2)-106.0660171779821*nuVtSqSum[0]*dxvl1R4*dxvr[1]-21.21320343559643*nuVtSqSum[0]*dxvl1R5)*fr[2]+(21.21320343559643*nuVtSqSum[0]*dxvr1R5+106.0660171779821*nuVtSqSum[0]*dxvl[1]*dxvr1R4+148.492424049175*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[2]+((122.4744871391589*dxvl1R3*dxvr1R2+61.23724356957945*dxvl1R4*dxvr[1]+12.24744871391589*dxvl1R5)*fr[1]+(12.24744871391589*dxvr1R5+61.23724356957945*dxvl[1]*dxvr1R4+122.4744871391589*dxvl1R2*dxvr1R3)*fl[1])*nuVtSqSum[1]+12.24744871391589*fl[0]*nuVtSqSum[0]*dxvr1R5+61.23724356957945*fl[0]*nuVtSqSum[0]*dxvl[1]*dxvr1R4+122.4744871391589*fl[0]*nuVtSqSum[0]*dxvl1R2*dxvr1R3+122.4744871391589*fr[0]*nuVtSqSum[0]*dxvl1R3*dxvr1R2+61.23724356957945*fr[0]*nuVtSqSum[0]*dxvl1R4*dxvr[1]+12.24744871391589*fr[0]*nuVtSqSum[0]*dxvl1R5)/(20.0*dxvr1R5+100.0*dxvl[1]*dxvr1R4+200.0*dxvl1R2*dxvr1R3+200.0*dxvl1R3*dxvr1R2+100.0*dxvl1R4*dxvr[1]+20.0*dxvl1R5); 
  incr2[3] = ((76.68115805072327*dxvl1R3*dxvr1R2+87.63560920082662*dxvl1R4*dxvr[1]+27.38612787525831*dxvl1R5)*nuVtSqSum[1]*fr[5]+(27.38612787525831*dxvr1R5+87.63560920082662*dxvl[1]*dxvr1R4+76.68115805072327*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[5]+(109.5445115010333*dxvl1R3*dxvr1R2+54.77225575051663*dxvl1R4*dxvr[1]+10.95445115010332*dxvl1R5)*nuVtSqSum[1]*fr[4]+(10.95445115010332*dxvr1R5+54.77225575051663*dxvl[1]*dxvr1R4+109.5445115010333*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[4]+(((-132.815661727072*dxvl1R3*dxvr1R2)-94.86832980505142*dxvl1R4*dxvr[1]-18.97366596101028*dxvl1R5)*nuVtSqSum[2]-148.492424049175*nuVtSqSum[0]*dxvl1R3*dxvr1R2-106.0660171779821*nuVtSqSum[0]*dxvl1R4*dxvr[1]-21.21320343559643*nuVtSqSum[0]*dxvl1R5)*fr[3]+((18.97366596101028*dxvr1R5+94.86832980505142*dxvl[1]*dxvr1R4+132.815661727072*dxvl1R2*dxvr1R3)*nuVtSqSum[2]+21.21320343559643*nuVtSqSum[0]*dxvr1R5+106.0660171779821*nuVtSqSum[0]*dxvl[1]*dxvr1R4+148.492424049175*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[3]+((109.5445115010333*dxvl1R3*dxvr1R2+54.77225575051663*dxvl1R4*dxvr[1]+10.95445115010332*dxvl1R5)*fr[1]+(10.95445115010332*dxvr1R5+54.77225575051663*dxvl[1]*dxvr1R4+109.5445115010333*dxvl1R2*dxvr1R3)*fl[1])*nuVtSqSum[2]+((-148.492424049175*dxvl1R3*dxvr1R2)-106.0660171779821*dxvl1R4*dxvr[1]-21.21320343559643*dxvl1R5)*nuVtSqSum[1]*fr[2]+(21.21320343559643*dxvr1R5+106.0660171779821*dxvl[1]*dxvr1R4+148.492424049175*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[2]+(12.24744871391589*fl[0]*dxvr1R5+61.23724356957945*fl[0]*dxvl[1]*dxvr1R4+122.4744871391589*fl[0]*dxvl1R2*dxvr1R3+122.4744871391589*fr[0]*dxvl1R3*dxvr1R2+61.23724356957945*fr[0]*dxvl1R4*dxvr[1]+12.24744871391589*fr[0]*dxvl1R5)*nuVtSqSum[1]+(122.4744871391589*nuVtSqSum[0]*dxvl1R3*dxvr1R2+61.23724356957945*nuVtSqSum[0]*dxvl1R4*dxvr[1]+12.24744871391589*nuVtSqSum[0]*dxvl1R5)*fr[1]+(12.24744871391589*nuVtSqSum[0]*dxvr1R5+61.23724356957945*nuVtSqSum[0]*dxvl[1]*dxvr1R4+122.4744871391589*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[1])/(20.0*dxvr1R5+100.0*dxvl[1]*dxvr1R4+200.0*dxvl1R2*dxvr1R3+200.0*dxvl1R3*dxvr1R2+100.0*dxvl1R4*dxvr[1]+20.0*dxvl1R5); 
  incr2[5] = -(1.0*((59.39696961967*nuVtSqSum[0]*dxvl1R3*dxvr1R2+67.8822509939086*nuVtSqSum[0]*dxvl1R4*dxvr[1]+21.21320343559643*nuVtSqSum[0]*dxvl1R5)*fr[5]+(21.21320343559643*nuVtSqSum[0]*dxvr1R5+67.8822509939086*nuVtSqSum[0]*dxvl[1]*dxvr1R4+59.39696961967*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[5]+(94.86832980505142*dxvl1R3*dxvr1R2+47.43416490252571*dxvl1R4*dxvr[1]+9.48683298050514*dxvl1R5)*nuVtSqSum[2]*fr[4]+(9.48683298050514*dxvr1R5+47.43416490252571*dxvl[1]*dxvr1R4+94.86832980505142*dxvl1R2*dxvr1R3)*nuVtSqSum[2]*fl[4]+((-115.0217370760849*dxvl1R3*dxvr1R2)-82.15838362577493*dxvl1R4*dxvr[1]-16.43167672515498*dxvl1R5)*nuVtSqSum[1]*fr[3]+(16.43167672515498*dxvr1R5+82.15838362577493*dxvl[1]*dxvr1R4+115.0217370760849*dxvl1R2*dxvr1R3)*nuVtSqSum[1]*fl[3]+((-115.0217370760849*nuVtSqSum[0]*dxvl1R3*dxvr1R2)-82.15838362577493*nuVtSqSum[0]*dxvl1R4*dxvr[1]-16.43167672515498*nuVtSqSum[0]*dxvl1R5)*fr[2]+(16.43167672515498*nuVtSqSum[0]*dxvr1R5+82.15838362577493*nuVtSqSum[0]*dxvl[1]*dxvr1R4+115.0217370760849*nuVtSqSum[0]*dxvl1R2*dxvr1R3)*fl[2]+((94.86832980505142*dxvl1R3*dxvr1R2+47.43416490252571*dxvl1R4*dxvr[1]+9.48683298050514*dxvl1R5)*fr[1]+(9.48683298050514*dxvr1R5+47.43416490252571*dxvl[1]*dxvr1R4+94.86832980505142*dxvl1R2*dxvr1R3)*fl[1])*nuVtSqSum[1]+9.48683298050514*fl[0]*nuVtSqSum[0]*dxvr1R5+47.43416490252571*fl[0]*nuVtSqSum[0]*dxvl[1]*dxvr1R4+94.86832980505142*fl[0]*nuVtSqSum[0]*dxvl1R2*dxvr1R3+94.86832980505142*fr[0]*nuVtSqSum[0]*dxvl1R3*dxvr1R2+47.43416490252571*fr[0]*nuVtSqSum[0]*dxvl1R4*dxvr[1]+9.48683298050514*fr[0]*nuVtSqSum[0]*dxvl1R5))/(4.0*dxvr1R5+20.0*dxvl[1]*dxvr1R4+40.0*dxvl1R2*dxvr1R3+40.0*dxvl1R3*dxvr1R2+20.0*dxvl1R4*dxvr[1]+4.0*dxvl1R5); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr2[5]*rdvSq4L-1.0*incr1[5]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUx[2])/nuSum-(0.7071067811865475*sumNuUx[0])/nuSum+wl[1]); 
} 
