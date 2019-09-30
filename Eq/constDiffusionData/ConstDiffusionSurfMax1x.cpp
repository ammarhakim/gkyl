#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf1xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[2]; 
  incr1[0] = 1.082531754730548*fr[1]+1.082531754730548*fl[1]-1.125*fr[0]+1.125*fl[0]; 
  incr1[1] = (-1.875*fr[1])-1.875*fl[1]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 

  double incr2[2]; 
  incr2[1] = (-1.0*fr[1])+fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 

} 
void ConstDiffusionSurf1xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[3]; 
  incr1[0] = (-1.341640786499874*fr[2])+1.341640786499874*fl[2]+2.381569860407206*fr[1]+2.381569860407206*fl[1]-1.875*fr[0]+1.875*fl[0]; 
  incr1[1] = 2.32379000772445*fr[2]-2.32379000772445*fl[2]-4.125*fr[1]-4.125*fl[1]+3.247595264191645*fr[0]-3.247595264191645*fl[0]; 
  incr1[2] = (-3.0*fr[2])+3.0*fl[2]+5.325352101035199*fr[1]+5.325352101035199*fl[1]-4.192627457812106*fr[0]+4.192627457812106*fl[0]; 

  double incr2[3]; 
  incr2[1] = 0.8472151069828725*fr[2]+0.8472151069828725*fl[2]-1.21875*fr[1]+1.21875*fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
  incr2[2] = (-3.28125*fr[2])-3.28125*fl[2]+4.720198453190289*fr[1]-4.720198453190289*fl[1]-3.354101966249685*fr[0]-3.354101966249685*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 
  outr[2] += incr2[2]*rdxSq2nur+incr1[2]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 
  outl[2] += incr2[2]*rdxSq2nul-1.0*incr1[2]*rdxSq2nul; 

} 
void ConstDiffusionSurf1xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[4]; 
  incr1[0] = 1.36421551976768*fr[3]+1.36421551976768*fl[3]-3.109532031210645*fr[2]+3.109532031210645*fl[2]+3.87005102316171*fr[1]+3.87005102316171*fl[1]-2.734375*fr[0]+2.734375*fl[0]; 
  incr1[1] = (-2.362890592711605*fr[3])-2.362890592711605*fl[3]+5.385867465819689*fr[2]-5.385867465819689*fl[2]-6.703125*fr[1]-6.703125*fl[1]+4.736076426946148*fr[0]-4.736076426946148*fl[0]; 
  incr1[2] = 3.05047863816074*fr[3]+3.05047863816074*fl[3]-6.953125*fr[2]+6.953125*fl[2]+8.653697164182198*fr[1]+8.653697164182198*fl[1]-6.114248375975989*fr[0]+6.114248375975989*fl[0]; 
  incr1[3] = (-3.609375*fr[3])-3.609375*fl[3]+8.227048448372905*fr[2]-8.227048448372905*fl[2]-10.23919256841696*fr[1]-10.23919256841696*fl[1]+7.23447624119224*fr[0]-7.23447624119224*fl[0]; 

  double incr2[4]; 
  incr2[1] = (-0.6546536707079771*fr[3])+0.6546536707079771*fl[3]+1.210307295689818*fr[2]+1.210307295689818*fl[2]-1.3125*fr[1]+1.3125*fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
  incr2[2] = 2.53546276418555*fr[3]-2.53546276418555*fl[3]-4.6875*fr[2]-4.6875*fl[2]+5.083290641897234*fr[1]-5.083290641897234*fl[1]-3.354101966249685*fr[0]-3.354101966249685*fl[0]; 
  incr2[3] = (-6.0*fr[3])+6.0*fl[3]+11.09264959331178*fr[2]+11.09264959331178*fl[2]-12.02926119925908*fr[1]+12.02926119925908*fl[1]+7.937253933193772*fr[0]+7.937253933193772*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 
  outr[2] += incr2[2]*rdxSq2nur+incr1[2]*rdxSq2nur; 
  outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 
  outl[2] += incr2[2]*rdxSq2nul-1.0*incr1[2]*rdxSq2nul; 
  outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 

} 
void ConstDiffusionSurf1xMax_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[5]; 
  incr1[0] = (-1.25*fr[4])+1.25*fl[4]+3.379533901242661*fr[3]+3.379533901242661*fl[3]-5.162172557431155*fr[2]+5.162172557431155*fl[2]+5.527677772592862*fr[1]+5.527677772592862*fl[1]-3.69140625*fr[0]+3.69140625*fl[0]; 
  incr1[1] = 2.165063509461096*fr[4]-2.165063509461096*fl[4]-5.853524422853748*fr[3]-5.853524422853748*fl[3]+8.941145146908527*fr[2]-8.941145146908527*fl[2]-9.57421875*fr[1]-9.57421875*fl[1]+6.393703176377298*fr[0]-6.393703176377298*fl[0]; 
  incr1[2] = (-2.795084971874738*fr[4])+2.795084971874738*fl[4]+7.556867535443652*fr[3]+7.556867535443652*fl[3]-11.54296875*fr[2]+11.54296875*fl[2]+12.36026325723227*fr[1]+12.36026325723227*fl[1]-8.254235307567585*fr[0]+8.254235307567585*fl[0]; 
  incr1[3] = 3.307189138830738*fr[4]-3.307189138830738*fl[4]-8.94140625*fr[3]-8.94140625*fl[3]+13.65782481176513*fr[2]-13.65782481176513*fl[2]-14.62486071398016*fr[1]-14.62486071398016*fl[1]+9.766542925609524*fr[0]-9.766542925609524*fl[0]; 
  incr1[4] = (-3.75*fr[4])+3.75*fl[4]+10.13860170372798*fr[3]+10.13860170372798*fl[3]-15.48651767229347*fr[2]+15.48651767229347*fl[2]+16.58303331777859*fr[1]+16.58303331777859*fl[1]-11.07421875*fr[0]+11.07421875*fl[0]; 

  double incr2[5]; 
  incr2[1] = 0.4837563778952138*fr[4]+0.4837563778952138*fl[4]-1.07020531715347*fr[3]+1.07020531715347*fl[3]+1.406982231239413*fr[2]+1.406982231239413*fl[2]-1.36328125*fr[1]+1.36328125*fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
  incr2[2] = (-1.873580395209785*fr[4])-1.873580395209785*fl[4]+4.144887370358018*fr[3]-4.144887370358018*fl[3]-5.44921875*fr[2]-5.44921875*fl[2]+5.27996557744683*fr[1]-5.27996557744683*fl[1]-3.354101966249685*fr[0]-3.354101966249685*fl[0]; 
  incr2[3] = 4.433700439244959*fr[4]+4.433700439244959*fl[4]-9.80859375*fr[3]+9.80859375*fl[3]+12.89520515222495*fr[2]+12.89520515222495*fl[2]-12.49467904327803*fr[1]+12.49467904327803*fl[1]+7.937253933193772*fr[0]+7.937253933193772*fl[0]; 
  incr2[4] = (-8.37890625*fr[4])-8.37890625*fl[4]+18.53649983840175*fr[3]-18.53649983840175*fl[3]-24.36964709853287*fr[2]-24.36964709853287*fl[2]+23.61272390006008*fr[1]-23.61272390006008*fl[1]-15.0*fr[0]-15.0*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 
  outr[2] += incr2[2]*rdxSq2nur+incr1[2]*rdxSq2nur; 
  outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 
  outr[4] += incr2[4]*rdxSq2nur+incr1[4]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 
  outl[2] += incr2[2]*rdxSq2nul-1.0*incr1[2]*rdxSq2nul; 
  outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  outl[4] += incr2[4]*rdxSq2nul-1.0*incr1[4]*rdxSq2nul; 

} 
