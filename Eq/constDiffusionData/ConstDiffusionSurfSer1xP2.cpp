#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf1xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = (-0.6708203932499369*fr[2])+0.6708203932499369*fl[2]+1.190784930203603*fr[1]+1.190784930203603*fl[1]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.161895003862225*fr[2]-1.161895003862225*fl[2]-2.0625*fr[1]-2.0625*fl[1]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[2] = (-1.5*fr[2])+1.5*fl[2]+2.662676050517599*fr[1]+2.662676050517599*fl[1]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 

  incr2[1] = 0.4236075534914363*fr[2]+0.4236075534914363*fl[2]-0.609375*fr[1]+0.609375*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[2] = (-1.640625*fr[2])-1.640625*fl[2]+2.360099226595144*fr[1]-2.360099226595144*fl[1]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul-1.0*incr1[2]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf1xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 
  double incr3[3]; 
  double incr4[3]; 

  incr1[0] = 6.708203932499369*fr[2]-6.708203932499369*fl[2]-8.11898816047911*fr[1]-8.11898816047911*fl[1]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = (-11.61895003862225*fr[2])+11.61895003862225*fl[2]+14.0625*fr[1]+14.0625*fl[1]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[2] = 15.0*fr[2]-15.0*fl[2]-18.15460943534727*fr[1]-18.15460943534727*fl[1]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 

  incr2[1] = (-2.541645320948617*fr[2])-2.541645320948617*fl[2]+1.40625*fr[1]-1.40625*fl[1]; 
  incr2[2] = 9.84375*fr[2]+9.84375*fl[2]-5.446382830604179*fr[1]+5.446382830604179*fl[1]; 

  incr3[2] = (-4.5*fr[2])+4.5*fl[2]+7.988028151552797*fr[1]+7.988028151552797*fl[1]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 


  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += (-1.0*incr2[1]*rdxFnur)-1.0*incr1[1]*rdxFnur; 
  outr[2] += (-1.0*incr3[2]*rdxFnur)-1.0*incr2[2]*rdxFnur-1.0*incr1[2]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr2[1]*rdxFnul-1.0*incr1[1]*rdxFnul; 
  outl[2] += incr3[2]*rdxFnul-1.0*incr2[2]*rdxFnul+incr1[2]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf1xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 
  double incr3[3]; 
  double incr4[3]; 
  double incr5[3]; 
  double incr6[3]; 

  incr1[0] = (-35.21807064562169*fr[2])+35.21807064562169*fl[2]+34.09975027401226*fr[1]+34.09975027401226*fl[1]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = 60.9994877027668*fr[2]-60.9994877027668*fl[2]-59.0625*fr[1]-59.0625*fl[1]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[2] = (-78.75*fr[2])+78.75*fl[2]+76.2493596284585*fr[1]+76.2493596284585*fl[1]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 

  incr2[1] = 9.531169953557313*fr[2]+9.531169953557313*fl[2]-2.4609375*fr[1]+2.4609375*fl[1]; 
  incr2[2] = (-36.9140625*fr[2])-36.9140625*fl[2]+9.531169953557313*fr[1]-9.531169953557313*fl[1]; 

  incr3[2] = 45.0*fr[2]-45.0*fl[2]-54.4638283060418*fr[1]-54.4638283060418*fl[1]+31.4447059335908*fr[0]-31.4447059335908*fl[0]; 




  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr3[2]*rdxFnur+incr2[2]*rdxFnur+incr1[2]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += (-1.0*incr3[2]*rdxFnul)+incr2[2]*rdxFnul-1.0*incr1[2]*rdxFnul; 

} 
void ConstDiffusionVarCoeffSurf1xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[3]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = (-1.060660171779821*fr[2]*nul[2])+1.060660171779821*fl[2]*nul[2]+1.882796291424007*fr[1]*nul[2]+1.882796291424007*fl[1]*nul[2]-1.482317653203927*fr[0]*nul[2]+1.482317653203927*fl[0]*nul[2]-0.8215838362577489*nul[1]*fr[2]-0.4743416490252568*nul[0]*fr[2]+0.8215838362577489*nul[1]*fl[2]+0.4743416490252568*nul[0]*fl[2]+1.458407736197253*fr[1]*nul[1]+1.458407736197253*fl[1]*nul[1]-1.148198316929614*fr[0]*nul[1]+1.148198316929614*fl[0]*nul[1]+0.8420120990817169*nul[0]*fr[1]+0.8420120990817169*nul[0]*fl[1]-0.6629126073623879*fr[0]*nul[0]+0.6629126073623879*fl[0]*nul[0]; 
  incr1[1] = 1.837117307087383*fr[2]*nul[2]-1.837117307087383*fl[2]*nul[2]-3.261098837048639*fr[1]*nul[2]-3.261098837048639*fl[1]*nul[2]+2.567449488305465*fr[0]*nul[2]-2.567449488305465*fl[0]*nul[2]+1.42302494707577*nul[1]*fr[2]+0.8215838362577489*nul[0]*fr[2]-1.42302494707577*nul[1]*fl[2]-0.8215838362577489*nul[0]*fl[2]-2.526036297245151*fr[1]*nul[1]-2.526036297245151*fl[1]*nul[1]+1.988737822087164*fr[0]*nul[1]-1.988737822087164*fl[0]*nul[1]-1.458407736197253*nul[0]*fr[1]-1.458407736197253*nul[0]*fl[1]+1.148198316929614*fr[0]*nul[0]-1.148198316929614*fl[0]*nul[0]; 
  incr1[2] = (-2.371708245126284*fr[2]*nul[2])+2.371708245126284*fl[2]*nul[2]+4.210060495408585*fr[1]*nul[2]+4.210060495408585*fl[1]*nul[2]-3.31456303681194*fr[0]*nul[2]+3.31456303681194*fl[0]*nul[2]-1.837117307087383*nul[1]*fr[2]-1.060660171779821*nul[0]*fr[2]+1.837117307087383*nul[1]*fl[2]+1.060660171779821*nul[0]*fl[2]+3.261098837048639*fr[1]*nul[1]+3.261098837048639*fl[1]*nul[1]-2.567449488305465*fr[0]*nul[1]+2.567449488305465*fl[0]*nul[1]+1.882796291424007*nul[0]*fr[1]+1.882796291424007*nul[0]*fl[1]-1.482317653203927*fr[0]*nul[0]+1.482317653203927*fl[0]*nul[0]; 

  incr2[1] = 0.6697823515422747*fr[2]*nul[2]+0.6697823515422747*fl[2]*nul[2]-0.9635064745825523*fr[1]*nul[2]+0.9635064745825523*fl[1]*nul[2]+0.6846531968814573*fr[0]*nul[2]+0.6846531968814573*fl[0]*nul[2]+0.5188111786213743*nul[1]*fr[2]+0.2995357736356374*nul[0]*fr[2]+0.5188111786213743*nul[1]*fl[2]+0.2995357736356374*nul[0]*fl[2]-0.7463289060042488*fr[1]*nul[1]+0.7463289060042488*fl[1]*nul[1]+0.5303300858899105*fr[0]*nul[1]+0.5303300858899105*fl[0]*nul[1]-0.430893194785552*nul[0]*fr[1]+0.430893194785552*nul[0]*fl[1]+0.3061862178478971*fr[0]*nul[0]+0.3061862178478971*fl[0]*nul[0]; 
  incr2[2] = (-2.594055893106872*fr[2]*nul[2])-2.594055893106872*fl[2]*nul[2]+3.731644530021244*fr[1]*nul[2]-3.731644530021244*fl[1]*nul[2]-2.651650429449552*fr[0]*nul[2]-2.651650429449552*fl[0]*nul[2]-2.009347054626824*nul[1]*fr[2]-1.160097062884178*nul[0]*fr[2]-2.009347054626824*nul[1]*fl[2]-1.160097062884178*nul[0]*fl[2]+2.890519423747657*fr[1]*nul[1]-2.890519423747657*fl[1]*nul[1]-2.053959590644372*fr[0]*nul[1]-2.053959590644372*fl[0]*nul[1]+1.668842167398551*nul[0]*fr[1]-1.668842167398551*nul[0]*fl[1]-1.185854122563142*fr[0]*nul[0]-1.185854122563142*fl[0]*nul[0]; 

  outr[0] += incr1[0]*rdxFr; 
  outr[1] += incr2[1]*rdxFr+incr1[1]*rdxFr; 
  outr[2] += incr2[2]*rdxFr+incr1[2]*rdxFr; 

  outl[0] += -1.0*incr1[0]*rdxFl; 
  outl[1] += incr1[1]*rdxFl-1.0*incr2[1]*rdxFl; 
  outl[2] += incr2[2]*rdxFl-1.0*incr1[2]*rdxFl; 

} 
