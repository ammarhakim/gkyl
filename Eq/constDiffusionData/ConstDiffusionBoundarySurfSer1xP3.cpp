#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf1xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 

  if (edge < 0) {

  incr2[1] = (-2.29128784747792*fr[3])+1.936491673103708*fr[2]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[2] = 8.874119674649426*fr[3]-7.5*fr[2]+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 
  incr2[3] = (-21.0*fr[3])+17.74823934929885*fr[2]-13.74772708486752*fr[1]+7.937253933193772*fr[0]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 

  } else {

  incr2[1] = 2.29128784747792*fl[3]+1.936491673103708*fl[2]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[2] = (-8.874119674649426*fl[3])-7.5*fl[2]-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 
  incr2[3] = 21.0*fl[3]+17.74823934929885*fl[2]+13.74772708486752*fl[1]+7.937253933193772*fl[0]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf1xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 
  double incr3[4]; 
  double incr4[4]; 

  if (edge < 0) {


  incr2[1] = 11.78376607274358*fr[3]-11.61895003862225*fr[2]+4.5*fr[1]; 
  incr2[2] = (-45.6383297553399*fr[3])+45.0*fr[2]-17.42842505793337*fr[1]; 
  incr2[3] = 108.0*fr[3]-106.4894360957931*fr[2]+41.24318125460255*fr[1]; 


  incr4[3] = (-16.2*fr[3])+28.39718295887816*fr[2]-30.24499958670854*fr[1]+19.84313483298443*fr[0]; 

  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[3] += (-1.0*incr4[3]*rdxFnur)-1.0*incr2[3]*rdxFnur; 

  } else {


  incr2[1] = 11.78376607274358*fl[3]+11.61895003862225*fl[2]+4.5*fl[1]; 
  incr2[2] = 45.6383297553399*fl[3]+45.0*fl[2]+17.42842505793337*fl[1]; 
  incr2[3] = 108.0*fl[3]+106.4894360957931*fl[2]+41.24318125460255*fl[1]; 


  incr4[3] = (-16.2*fl[3])-28.39718295887816*fl[2]-30.24499958670854*fl[1]-19.84313483298443*fl[0]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += (-1.0*incr4[3]*rdxFnul)-1.0*incr2[3]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf1xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 
  double incr3[4]; 
  double incr4[4]; 
  double incr5[4]; 
  double incr6[4]; 

  if (edge < 0) {


  incr2[1] = (-103.1079531365064*fr[3])+76.2493596284585*fr[2]-19.6875*fr[1]; 
  incr2[2] = 399.3353853592242*fr[3]-295.3125*fr[2]+76.2493596284585*fr[1]; 
  incr2[3] = (-945.0*fr[3])+698.8369243786424*fr[2]-180.4389179888861*fr[1]; 


  incr4[3] = 225.0*fr[3]-241.2651286545313*fr[2]+96.66370606547471*fr[1]; 



  outr[1] += incr2[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr4[3]*rdxFnur+incr2[3]*rdxFnur; 

  } else {


  incr2[1] = (-103.1079531365064*fl[3])-76.2493596284585*fl[2]-19.6875*fl[1]; 
  incr2[2] = (-399.3353853592242*fl[3])-295.3125*fl[2]-76.2493596284585*fl[1]; 
  incr2[3] = (-945.0*fl[3])-698.8369243786424*fl[2]-180.4389179888861*fl[1]; 


  incr4[3] = 225.0*fl[3]+241.2651286545313*fl[2]+96.66370606547471*fl[1]; 



  outl[1] += incr2[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul; 
  outl[3] += incr4[3]*rdxFnul+incr2[3]*rdxFnul; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf1xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 

  if (edge < 0) {

  incr2[1] = (-4.28660704987056*fr[3]*nul[3])+3.622844186547359*fr[2]*nul[3]-2.806243040080455*fr[1]*nul[3]+1.620185174601965*fr[0]*nul[3]-3.62284418654736*nul[2]*fr[3]-2.806243040080455*nul[1]*fr[3]-1.620185174601965*nul[0]*fr[3]+3.06186217847897*fr[2]*nul[2]-2.371708245126284*fr[1]*nul[2]+1.369306393762915*fr[0]*nul[2]+2.371708245126284*nul[1]*fr[2]+1.369306393762915*nul[0]*fr[2]-1.837117307087383*fr[1]*nul[1]+1.060660171779821*fr[0]*nul[1]-1.060660171779821*nul[0]*fr[1]+0.6123724356957944*fr[0]*nul[0]; 
  incr2[2] = 16.60195771588399*fr[3]*nul[3]-14.03121520040228*fr[2]*nul[3]+10.86853255964208*fr[1]*nul[3]-6.274950199005565*fr[0]*nul[3]+14.03121520040228*nul[2]*fr[3]+10.86853255964208*nul[1]*fr[3]+6.274950199005568*nul[0]*fr[3]-11.85854122563142*fr[2]*nul[2]+9.185586535436915*fr[1]*nul[2]-5.303300858899105*fr[0]*nul[2]-9.18558653543691*nul[1]*fr[2]-5.303300858899105*nul[0]*fr[2]+7.115124735378852*fr[1]*nul[1]-4.107919181288745*fr[0]*nul[1]+4.107919181288745*nul[0]*fr[1]-2.371708245126284*fr[0]*nul[0]; 
  incr2[3] = (-39.28740256112638*fr[3]*nul[3])+33.20391543176798*fr[2]*nul[3]-25.71964229922336*fr[1]*nul[3]+14.8492424049175*fr[0]*nul[3]-33.20391543176798*nul[2]*fr[3]-25.71964229922336*nul[1]*fr[3]-14.8492424049175*nul[0]*fr[3]+28.06243040080456*fr[2]*nul[2]-21.73706511928415*fr[1]*nul[2]+12.54990039801113*fr[0]*nul[2]+21.73706511928415*nul[1]*fr[2]+12.54990039801113*nul[0]*fr[2]-16.83745824048274*fr[1]*nul[1]+9.721111047611789*fr[0]*nul[1]-9.721111047611789*nul[0]*fr[1]+5.612486080160912*fr[0]*nul[0]; 

  outr[1] += incr2[1]*rdxFr; 
  outr[2] += incr2[2]*rdxFr; 
  outr[3] += incr2[3]*rdxFr; 

  } else {

  incr2[1] = 4.28660704987056*fl[3]*nul[3]+3.622844186547359*fl[2]*nul[3]+2.806243040080455*fl[1]*nul[3]+1.620185174601965*fl[0]*nul[3]+3.62284418654736*nul[2]*fl[3]+2.806243040080455*nul[1]*fl[3]+1.620185174601965*nul[0]*fl[3]+3.06186217847897*fl[2]*nul[2]+2.371708245126284*fl[1]*nul[2]+1.369306393762915*fl[0]*nul[2]+2.371708245126284*nul[1]*fl[2]+1.369306393762915*nul[0]*fl[2]+1.837117307087383*fl[1]*nul[1]+1.060660171779821*fl[0]*nul[1]+1.060660171779821*nul[0]*fl[1]+0.6123724356957944*fl[0]*nul[0]; 
  incr2[2] = (-16.60195771588399*fl[3]*nul[3])-14.03121520040228*fl[2]*nul[3]-10.86853255964208*fl[1]*nul[3]-6.274950199005565*fl[0]*nul[3]-14.03121520040228*nul[2]*fl[3]-10.86853255964208*nul[1]*fl[3]-6.274950199005568*nul[0]*fl[3]-11.85854122563142*fl[2]*nul[2]-9.185586535436915*fl[1]*nul[2]-5.303300858899105*fl[0]*nul[2]-9.18558653543691*nul[1]*fl[2]-5.303300858899105*nul[0]*fl[2]-7.115124735378852*fl[1]*nul[1]-4.107919181288745*fl[0]*nul[1]-4.107919181288745*nul[0]*fl[1]-2.371708245126284*fl[0]*nul[0]; 
  incr2[3] = 39.28740256112638*fl[3]*nul[3]+33.20391543176798*fl[2]*nul[3]+25.71964229922336*fl[1]*nul[3]+14.8492424049175*fl[0]*nul[3]+33.20391543176798*nul[2]*fl[3]+25.71964229922336*nul[1]*fl[3]+14.8492424049175*nul[0]*fl[3]+28.06243040080456*fl[2]*nul[2]+21.73706511928415*fl[1]*nul[2]+12.54990039801113*fl[0]*nul[2]+21.73706511928415*nul[1]*fl[2]+12.54990039801113*nul[0]*fl[2]+16.83745824048274*fl[1]*nul[1]+9.721111047611789*fl[0]*nul[1]+9.721111047611789*nul[0]*fl[1]+5.612486080160912*fl[0]*nul[0]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 
  outl[2] += incr2[2]*rdxFl; 
  outl[3] += -1.0*incr2[3]*rdxFl; 

  }

} 
