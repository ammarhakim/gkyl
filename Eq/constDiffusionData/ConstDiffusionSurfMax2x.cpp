#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf2xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[3]; 
  incr1[0] = 1.082531754730548*fr[1]+1.082531754730548*fl[1]-1.125*fr[0]+1.125*fl[0]; 
  incr1[1] = (-1.875*fr[1])-1.875*fl[1]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 
  incr1[2] = 1.125*fl[2]-1.125*fr[2]; 

  double incr2[3]; 
  incr2[1] = (-1.0*fr[1])+fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 
  outr[2] += incr1[2]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 
  outl[2] += -1.0*incr1[2]*rdxSq2nul; 

} 
void ConstDiffusionSurf2xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[6]; 
  incr1[0] = (-1.341640786499874*fr[4])+1.341640786499874*fl[4]+2.381569860407206*fr[1]+2.381569860407206*fl[1]-1.875*fr[0]+1.875*fl[0]; 
  incr1[1] = 2.32379000772445*fr[4]-2.32379000772445*fl[4]-4.125*fr[1]-4.125*fl[1]+3.247595264191645*fr[0]-3.247595264191645*fl[0]; 
  incr1[2] = 2.381569860407206*fr[3]+2.381569860407206*fl[3]-1.875*fr[2]+1.875*fl[2]; 
  incr1[3] = (-4.125*fr[3])-4.125*fl[3]+3.247595264191645*fr[2]-3.247595264191645*fl[2]; 
  incr1[4] = (-3.0*fr[4])+3.0*fl[4]+5.325352101035199*fr[1]+5.325352101035199*fl[1]-4.192627457812106*fr[0]+4.192627457812106*fl[0]; 
  incr1[5] = 1.875*fl[5]-1.875*fr[5]; 

  double incr2[6]; 
  incr2[1] = 0.8472151069828725*fr[4]+0.8472151069828725*fl[4]-1.21875*fr[1]+1.21875*fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
  incr2[3] = (-1.21875*fr[3])+1.21875*fl[3]+0.8660254037844386*fr[2]+0.8660254037844386*fl[2]; 
  incr2[4] = (-3.28125*fr[4])-3.28125*fl[4]+4.720198453190289*fr[1]-4.720198453190289*fl[1]-3.354101966249685*fr[0]-3.354101966249685*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 
  outr[2] += incr1[2]*rdxSq2nur; 
  outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 
  outr[4] += incr2[4]*rdxSq2nur+incr1[4]*rdxSq2nur; 
  outr[5] += incr1[5]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 
  outl[2] += -1.0*incr1[2]*rdxSq2nul; 
  outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  outl[4] += incr2[4]*rdxSq2nul-1.0*incr1[4]*rdxSq2nul; 
  outl[5] += -1.0*incr1[5]*rdxSq2nul; 

} 
void ConstDiffusionSurf2xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[10]; 
  incr1[0] = 1.36421551976768*fr[8]+1.36421551976768*fl[8]-3.109532031210645*fr[4]+3.109532031210645*fl[4]+3.87005102316171*fr[1]+3.87005102316171*fl[1]-2.734375*fr[0]+2.734375*fl[0]; 
  incr1[1] = (-2.362890592711605*fr[8])-2.362890592711605*fl[8]+5.385867465819689*fr[4]-5.385867465819689*fl[4]-6.703125*fr[1]-6.703125*fl[1]+4.736076426946148*fr[0]-4.736076426946148*fl[0]; 
  incr1[2] = (-3.109532031210645*fr[6])+3.109532031210645*fl[6]+3.87005102316171*fr[3]+3.87005102316171*fl[3]-2.734375*fr[2]+2.734375*fl[2]; 
  incr1[3] = 5.38586746581969*fr[6]-5.38586746581969*fl[6]-6.703125*fr[3]-6.703125*fl[3]+4.736076426946148*fr[2]-4.736076426946148*fl[2]; 
  incr1[4] = 3.05047863816074*fr[8]+3.05047863816074*fl[8]-6.953125*fr[4]+6.953125*fl[4]+8.653697164182198*fr[1]+8.653697164182198*fl[1]-6.114248375975989*fr[0]+6.114248375975989*fl[0]; 
  incr1[5] = 3.87005102316171*fr[7]+3.87005102316171*fl[7]-2.734375*fr[5]+2.734375*fl[5]; 
  incr1[6] = (-6.953125*fr[6])+6.953125*fl[6]+8.653697164182198*fr[3]+8.653697164182198*fl[3]-6.114248375975989*fr[2]+6.114248375975989*fl[2]; 
  incr1[7] = (-6.703125*fr[7])-6.703125*fl[7]+4.73607642694615*fr[5]-4.73607642694615*fl[5]; 
  incr1[8] = (-3.609375*fr[8])-3.609375*fl[8]+8.227048448372905*fr[4]-8.227048448372905*fl[4]-10.23919256841696*fr[1]-10.23919256841696*fl[1]+7.23447624119224*fr[0]-7.23447624119224*fl[0]; 
  incr1[9] = 2.734375*fl[9]-2.734375*fr[9]; 

  double incr2[10]; 
  incr2[1] = (-0.6546536707079771*fr[8])+0.6546536707079771*fl[8]+1.210307295689818*fr[4]+1.210307295689818*fl[4]-1.3125*fr[1]+1.3125*fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
  incr2[3] = 1.210307295689818*fr[6]+1.210307295689818*fl[6]-1.3125*fr[3]+1.3125*fl[3]+0.8660254037844386*fr[2]+0.8660254037844386*fl[2]; 
  incr2[4] = 2.53546276418555*fr[8]-2.53546276418555*fl[8]-4.6875*fr[4]-4.6875*fl[4]+5.083290641897234*fr[1]-5.083290641897234*fl[1]-3.354101966249685*fr[0]-3.354101966249685*fl[0]; 
  incr2[6] = (-4.6875*fr[6])-4.6875*fl[6]+5.083290641897235*fr[3]-5.083290641897235*fl[3]-3.354101966249684*fr[2]-3.354101966249684*fl[2]; 
  incr2[7] = (-1.3125*fr[7])+1.3125*fl[7]+0.8660254037844387*fr[5]+0.8660254037844387*fl[5]; 
  incr2[8] = (-6.0*fr[8])+6.0*fl[8]+11.09264959331178*fr[4]+11.09264959331178*fl[4]-12.02926119925908*fr[1]+12.02926119925908*fl[1]+7.937253933193772*fr[0]+7.937253933193772*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 
  outr[2] += incr1[2]*rdxSq2nur; 
  outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 
  outr[4] += incr2[4]*rdxSq2nur+incr1[4]*rdxSq2nur; 
  outr[5] += incr1[5]*rdxSq2nur; 
  outr[6] += incr2[6]*rdxSq2nur+incr1[6]*rdxSq2nur; 
  outr[7] += incr2[7]*rdxSq2nur+incr1[7]*rdxSq2nur; 
  outr[8] += incr2[8]*rdxSq2nur+incr1[8]*rdxSq2nur; 
  outr[9] += incr1[9]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 
  outl[2] += -1.0*incr1[2]*rdxSq2nul; 
  outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  outl[4] += incr2[4]*rdxSq2nul-1.0*incr1[4]*rdxSq2nul; 
  outl[5] += -1.0*incr1[5]*rdxSq2nul; 
  outl[6] += incr2[6]*rdxSq2nul-1.0*incr1[6]*rdxSq2nul; 
  outl[7] += incr1[7]*rdxSq2nul-1.0*incr2[7]*rdxSq2nul; 
  outl[8] += incr1[8]*rdxSq2nul-1.0*incr2[8]*rdxSq2nul; 
  outl[9] += -1.0*incr1[9]*rdxSq2nul; 

} 
void ConstDiffusionSurf2xMax_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[15]; 
  incr1[0] = (-1.25*fr[13])+1.25*fl[13]+3.379533901242661*fr[8]+3.379533901242661*fl[8]-5.162172557431155*fr[4]+5.162172557431155*fl[4]+5.527677772592862*fr[1]+5.527677772592862*fl[1]-3.69140625*fr[0]+3.69140625*fl[0]; 
  incr1[1] = 2.165063509461096*fr[13]-2.165063509461096*fl[13]-5.853524422853748*fr[8]-5.853524422853748*fl[8]+8.941145146908527*fr[4]-8.941145146908527*fl[4]-9.57421875*fr[1]-9.57421875*fl[1]+6.393703176377298*fr[0]-6.393703176377298*fl[0]; 
  incr1[2] = 3.37953390124266*fr[11]+3.37953390124266*fl[11]-5.162172557431155*fr[6]+5.162172557431155*fl[6]+5.527677772592862*fr[3]+5.527677772592862*fl[3]-3.69140625*fr[2]+3.69140625*fl[2]; 
  incr1[3] = (-5.853524422853749*fr[11])-5.853524422853749*fl[11]+8.941145146908529*fr[6]-8.941145146908529*fl[6]-9.57421875*fr[3]-9.57421875*fl[3]+6.393703176377298*fr[2]-6.393703176377298*fl[2]; 
  incr1[4] = (-2.795084971874738*fr[13])+2.795084971874738*fl[13]+7.556867535443652*fr[8]+7.556867535443652*fl[8]-11.54296875*fr[4]+11.54296875*fl[4]+12.36026325723227*fr[1]+12.36026325723227*fl[1]-8.254235307567585*fr[0]+8.254235307567585*fl[0]; 
  incr1[5] = (-5.162172557431155*fr[10])+5.162172557431155*fl[10]+5.527677772592862*fr[7]+5.527677772592862*fl[7]-3.69140625*fr[5]+3.69140625*fl[5]; 
  incr1[6] = 7.556867535443651*fr[11]+7.556867535443651*fl[11]-11.54296875*fr[6]+11.54296875*fl[6]+12.36026325723226*fr[3]+12.36026325723226*fl[3]-8.254235307567583*fr[2]+8.254235307567583*fl[2]; 
  incr1[7] = 8.941145146908529*fr[10]-8.941145146908529*fl[10]-9.57421875*fr[7]-9.57421875*fl[7]+6.393703176377302*fr[5]-6.393703176377302*fl[5]; 
  incr1[8] = 3.307189138830738*fr[13]-3.307189138830738*fl[13]-8.94140625*fr[8]-8.94140625*fl[8]+13.65782481176513*fr[4]-13.65782481176513*fl[4]-14.62486071398016*fr[1]-14.62486071398016*fl[1]+9.766542925609524*fr[0]-9.766542925609524*fl[0]; 
  incr1[9] = 5.527677772592861*fr[12]+5.527677772592861*fl[12]-3.69140625*fr[9]+3.69140625*fl[9]; 
  incr1[10] = (-11.54296875*fr[10])+11.54296875*fl[10]+12.36026325723226*fr[7]+12.36026325723226*fl[7]-8.254235307567585*fr[5]+8.254235307567585*fl[5]; 
  incr1[11] = (-8.94140625*fr[11])-8.94140625*fl[11]+13.65782481176513*fr[6]-13.65782481176513*fl[6]-14.62486071398016*fr[3]-14.62486071398016*fl[3]+9.76654292560952*fr[2]-9.76654292560952*fl[2]; 
  incr1[12] = (-9.57421875*fr[12])-9.57421875*fl[12]+6.393703176377302*fr[9]-6.393703176377302*fl[9]; 
  incr1[13] = (-3.75*fr[13])+3.75*fl[13]+10.13860170372798*fr[8]+10.13860170372798*fl[8]-15.48651767229347*fr[4]+15.48651767229347*fl[4]+16.58303331777859*fr[1]+16.58303331777859*fl[1]-11.07421875*fr[0]+11.07421875*fl[0]; 
  incr1[14] = 3.69140625*fl[14]-3.69140625*fr[14]; 

  double incr2[15]; 
  incr2[1] = 0.4837563778952138*fr[13]+0.4837563778952138*fl[13]-1.07020531715347*fr[8]+1.07020531715347*fl[8]+1.406982231239413*fr[4]+1.406982231239413*fl[4]-1.36328125*fr[1]+1.36328125*fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
  incr2[3] = (-1.07020531715347*fr[11])+1.07020531715347*fl[11]+1.406982231239413*fr[6]+1.406982231239413*fl[6]-1.36328125*fr[3]+1.36328125*fl[3]+0.8660254037844386*fr[2]+0.8660254037844386*fl[2]; 
  incr2[4] = (-1.873580395209785*fr[13])-1.873580395209785*fl[13]+4.144887370358018*fr[8]-4.144887370358018*fl[8]-5.44921875*fr[4]-5.44921875*fl[4]+5.27996557744683*fr[1]-5.27996557744683*fl[1]-3.354101966249685*fr[0]-3.354101966249685*fl[0]; 
  incr2[6] = 4.144887370358018*fr[11]-4.144887370358018*fl[11]-5.44921875*fr[6]-5.44921875*fl[6]+5.279965577446831*fr[3]-5.279965577446831*fl[3]-3.354101966249684*fr[2]-3.354101966249684*fl[2]; 
  incr2[7] = 1.406982231239413*fr[10]+1.406982231239413*fl[10]-1.36328125*fr[7]+1.36328125*fl[7]+0.8660254037844387*fr[5]+0.8660254037844387*fl[5]; 
  incr2[8] = 4.433700439244959*fr[13]+4.433700439244959*fl[13]-9.80859375*fr[8]+9.80859375*fl[8]+12.89520515222495*fr[4]+12.89520515222495*fl[4]-12.49467904327803*fr[1]+12.49467904327803*fl[1]+7.937253933193772*fr[0]+7.937253933193772*fl[0]; 
  incr2[10] = (-5.44921875*fr[10])-5.44921875*fl[10]+5.279965577446831*fr[7]-5.279965577446831*fl[7]-3.354101966249685*fr[5]-3.354101966249685*fl[5]; 
  incr2[11] = (-9.80859375*fr[11])+9.80859375*fl[11]+12.89520515222494*fr[6]+12.89520515222494*fl[6]-12.49467904327803*fr[3]+12.49467904327803*fl[3]+7.93725393319377*fr[2]+7.93725393319377*fl[2]; 
  incr2[12] = (-1.36328125*fr[12])+1.36328125*fl[12]+0.8660254037844387*fr[9]+0.8660254037844387*fl[9]; 
  incr2[13] = (-8.37890625*fr[13])-8.37890625*fl[13]+18.53649983840175*fr[8]-18.53649983840175*fl[8]-24.36964709853287*fr[4]-24.36964709853287*fl[4]+23.61272390006008*fr[1]-23.61272390006008*fl[1]-15.0*fr[0]-15.0*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 
  outr[2] += incr1[2]*rdxSq2nur; 
  outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 
  outr[4] += incr2[4]*rdxSq2nur+incr1[4]*rdxSq2nur; 
  outr[5] += incr1[5]*rdxSq2nur; 
  outr[6] += incr2[6]*rdxSq2nur+incr1[6]*rdxSq2nur; 
  outr[7] += incr2[7]*rdxSq2nur+incr1[7]*rdxSq2nur; 
  outr[8] += incr2[8]*rdxSq2nur+incr1[8]*rdxSq2nur; 
  outr[9] += incr1[9]*rdxSq2nur; 
  outr[10] += incr2[10]*rdxSq2nur+incr1[10]*rdxSq2nur; 
  outr[11] += incr2[11]*rdxSq2nur+incr1[11]*rdxSq2nur; 
  outr[12] += incr2[12]*rdxSq2nur+incr1[12]*rdxSq2nur; 
  outr[13] += incr2[13]*rdxSq2nur+incr1[13]*rdxSq2nur; 
  outr[14] += incr1[14]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 
  outl[2] += -1.0*incr1[2]*rdxSq2nul; 
  outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  outl[4] += incr2[4]*rdxSq2nul-1.0*incr1[4]*rdxSq2nul; 
  outl[5] += -1.0*incr1[5]*rdxSq2nul; 
  outl[6] += incr2[6]*rdxSq2nul-1.0*incr1[6]*rdxSq2nul; 
  outl[7] += incr1[7]*rdxSq2nul-1.0*incr2[7]*rdxSq2nul; 
  outl[8] += incr1[8]*rdxSq2nul-1.0*incr2[8]*rdxSq2nul; 
  outl[9] += -1.0*incr1[9]*rdxSq2nul; 
  outl[10] += incr2[10]*rdxSq2nul-1.0*incr1[10]*rdxSq2nul; 
  outl[11] += incr1[11]*rdxSq2nul-1.0*incr2[11]*rdxSq2nul; 
  outl[12] += incr1[12]*rdxSq2nul-1.0*incr2[12]*rdxSq2nul; 
  outl[13] += incr2[13]*rdxSq2nul-1.0*incr1[13]*rdxSq2nul; 
  outl[14] += -1.0*incr1[14]*rdxSq2nul; 

} 
void ConstDiffusionSurf2xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 

  double incr1[3]; 
  incr1[0] = 1.082531754730548*fr[2]+1.082531754730548*fl[2]-1.125*fr[0]+1.125*fl[0]; 
  incr1[1] = 1.125*fl[1]-1.125*fr[1]; 
  incr1[2] = (-1.875*fr[2])-1.875*fl[2]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 

  double incr2[3]; 
  incr2[2] = (-1.0*fr[2])+fl[2]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr1[1]*rdxSq2nur; 
  outr[2] += incr2[2]*rdxSq2nur+incr1[2]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += -1.0*incr1[1]*rdxSq2nul; 
  outl[2] += incr1[2]*rdxSq2nul-1.0*incr2[2]*rdxSq2nul; 

} 
void ConstDiffusionSurf2xMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 

  double incr1[6]; 
  incr1[0] = (-1.341640786499874*fr[5])+1.341640786499874*fl[5]+2.381569860407206*fr[2]+2.381569860407206*fl[2]-1.875*fr[0]+1.875*fl[0]; 
  incr1[1] = 2.381569860407206*fr[3]+2.381569860407206*fl[3]-1.875*fr[1]+1.875*fl[1]; 
  incr1[2] = 2.32379000772445*fr[5]-2.32379000772445*fl[5]-4.125*fr[2]-4.125*fl[2]+3.247595264191645*fr[0]-3.247595264191645*fl[0]; 
  incr1[3] = (-4.125*fr[3])-4.125*fl[3]+3.247595264191645*fr[1]-3.247595264191645*fl[1]; 
  incr1[4] = 1.875*fl[4]-1.875*fr[4]; 
  incr1[5] = (-3.0*fr[5])+3.0*fl[5]+5.325352101035199*fr[2]+5.325352101035199*fl[2]-4.192627457812106*fr[0]+4.192627457812106*fl[0]; 

  double incr2[6]; 
  incr2[2] = 0.8472151069828725*fr[5]+0.8472151069828725*fl[5]-1.21875*fr[2]+1.21875*fl[2]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
  incr2[3] = (-1.21875*fr[3])+1.21875*fl[3]+0.8660254037844386*fr[1]+0.8660254037844386*fl[1]; 
  incr2[5] = (-3.28125*fr[5])-3.28125*fl[5]+4.720198453190289*fr[2]-4.720198453190289*fl[2]-3.354101966249685*fr[0]-3.354101966249685*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr1[1]*rdxSq2nur; 
  outr[2] += incr2[2]*rdxSq2nur+incr1[2]*rdxSq2nur; 
  outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 
  outr[4] += incr1[4]*rdxSq2nur; 
  outr[5] += incr2[5]*rdxSq2nur+incr1[5]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += -1.0*incr1[1]*rdxSq2nul; 
  outl[2] += incr1[2]*rdxSq2nul-1.0*incr2[2]*rdxSq2nul; 
  outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  outl[4] += -1.0*incr1[4]*rdxSq2nul; 
  outl[5] += incr2[5]*rdxSq2nul-1.0*incr1[5]*rdxSq2nul; 

} 
void ConstDiffusionSurf2xMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 

  double incr1[10]; 
  incr1[0] = 1.36421551976768*fr[9]+1.36421551976768*fl[9]-3.109532031210645*fr[5]+3.109532031210645*fl[5]+3.87005102316171*fr[2]+3.87005102316171*fl[2]-2.734375*fr[0]+2.734375*fl[0]; 
  incr1[1] = (-3.109532031210645*fr[7])+3.109532031210645*fl[7]+3.87005102316171*fr[3]+3.87005102316171*fl[3]-2.734375*fr[1]+2.734375*fl[1]; 
  incr1[2] = (-2.362890592711605*fr[9])-2.362890592711605*fl[9]+5.385867465819689*fr[5]-5.385867465819689*fl[5]-6.703125*fr[2]-6.703125*fl[2]+4.736076426946148*fr[0]-4.736076426946148*fl[0]; 
  incr1[3] = 5.38586746581969*fr[7]-5.38586746581969*fl[7]-6.703125*fr[3]-6.703125*fl[3]+4.736076426946148*fr[1]-4.736076426946148*fl[1]; 
  incr1[4] = 3.87005102316171*fr[6]+3.87005102316171*fl[6]-2.734375*fr[4]+2.734375*fl[4]; 
  incr1[5] = 3.05047863816074*fr[9]+3.05047863816074*fl[9]-6.953125*fr[5]+6.953125*fl[5]+8.653697164182198*fr[2]+8.653697164182198*fl[2]-6.114248375975989*fr[0]+6.114248375975989*fl[0]; 
  incr1[6] = (-6.703125*fr[6])-6.703125*fl[6]+4.73607642694615*fr[4]-4.73607642694615*fl[4]; 
  incr1[7] = (-6.953125*fr[7])+6.953125*fl[7]+8.653697164182198*fr[3]+8.653697164182198*fl[3]-6.114248375975989*fr[1]+6.114248375975989*fl[1]; 
  incr1[8] = 2.734375*fl[8]-2.734375*fr[8]; 
  incr1[9] = (-3.609375*fr[9])-3.609375*fl[9]+8.227048448372905*fr[5]-8.227048448372905*fl[5]-10.23919256841696*fr[2]-10.23919256841696*fl[2]+7.23447624119224*fr[0]-7.23447624119224*fl[0]; 

  double incr2[10]; 
  incr2[2] = (-0.6546536707079771*fr[9])+0.6546536707079771*fl[9]+1.210307295689818*fr[5]+1.210307295689818*fl[5]-1.3125*fr[2]+1.3125*fl[2]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
  incr2[3] = 1.210307295689818*fr[7]+1.210307295689818*fl[7]-1.3125*fr[3]+1.3125*fl[3]+0.8660254037844386*fr[1]+0.8660254037844386*fl[1]; 
  incr2[5] = 2.53546276418555*fr[9]-2.53546276418555*fl[9]-4.6875*fr[5]-4.6875*fl[5]+5.083290641897234*fr[2]-5.083290641897234*fl[2]-3.354101966249685*fr[0]-3.354101966249685*fl[0]; 
  incr2[6] = (-1.3125*fr[6])+1.3125*fl[6]+0.8660254037844387*fr[4]+0.8660254037844387*fl[4]; 
  incr2[7] = (-4.6875*fr[7])-4.6875*fl[7]+5.083290641897235*fr[3]-5.083290641897235*fl[3]-3.354101966249684*fr[1]-3.354101966249684*fl[1]; 
  incr2[9] = (-6.0*fr[9])+6.0*fl[9]+11.09264959331178*fr[5]+11.09264959331178*fl[5]-12.02926119925908*fr[2]+12.02926119925908*fl[2]+7.937253933193772*fr[0]+7.937253933193772*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr1[1]*rdxSq2nur; 
  outr[2] += incr2[2]*rdxSq2nur+incr1[2]*rdxSq2nur; 
  outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 
  outr[4] += incr1[4]*rdxSq2nur; 
  outr[5] += incr2[5]*rdxSq2nur+incr1[5]*rdxSq2nur; 
  outr[6] += incr2[6]*rdxSq2nur+incr1[6]*rdxSq2nur; 
  outr[7] += incr2[7]*rdxSq2nur+incr1[7]*rdxSq2nur; 
  outr[8] += incr1[8]*rdxSq2nur; 
  outr[9] += incr2[9]*rdxSq2nur+incr1[9]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += -1.0*incr1[1]*rdxSq2nul; 
  outl[2] += incr1[2]*rdxSq2nul-1.0*incr2[2]*rdxSq2nul; 
  outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  outl[4] += -1.0*incr1[4]*rdxSq2nul; 
  outl[5] += incr2[5]*rdxSq2nul-1.0*incr1[5]*rdxSq2nul; 
  outl[6] += incr1[6]*rdxSq2nul-1.0*incr2[6]*rdxSq2nul; 
  outl[7] += incr2[7]*rdxSq2nul-1.0*incr1[7]*rdxSq2nul; 
  outl[8] += -1.0*incr1[8]*rdxSq2nul; 
  outl[9] += incr1[9]*rdxSq2nul-1.0*incr2[9]*rdxSq2nul; 

} 
void ConstDiffusionSurf2xMax_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 

  double incr1[15]; 
  incr1[0] = (-1.25*fr[14])+1.25*fl[14]+3.379533901242661*fr[9]+3.379533901242661*fl[9]-5.162172557431155*fr[5]+5.162172557431155*fl[5]+5.527677772592862*fr[2]+5.527677772592862*fl[2]-3.69140625*fr[0]+3.69140625*fl[0]; 
  incr1[1] = 3.37953390124266*fr[12]+3.37953390124266*fl[12]-5.162172557431155*fr[7]+5.162172557431155*fl[7]+5.527677772592862*fr[3]+5.527677772592862*fl[3]-3.69140625*fr[1]+3.69140625*fl[1]; 
  incr1[2] = 2.165063509461096*fr[14]-2.165063509461096*fl[14]-5.853524422853748*fr[9]-5.853524422853748*fl[9]+8.941145146908527*fr[5]-8.941145146908527*fl[5]-9.57421875*fr[2]-9.57421875*fl[2]+6.393703176377298*fr[0]-6.393703176377298*fl[0]; 
  incr1[3] = (-5.853524422853749*fr[12])-5.853524422853749*fl[12]+8.941145146908529*fr[7]-8.941145146908529*fl[7]-9.57421875*fr[3]-9.57421875*fl[3]+6.393703176377298*fr[1]-6.393703176377298*fl[1]; 
  incr1[4] = (-5.162172557431155*fr[10])+5.162172557431155*fl[10]+5.527677772592862*fr[6]+5.527677772592862*fl[6]-3.69140625*fr[4]+3.69140625*fl[4]; 
  incr1[5] = (-2.795084971874738*fr[14])+2.795084971874738*fl[14]+7.556867535443652*fr[9]+7.556867535443652*fl[9]-11.54296875*fr[5]+11.54296875*fl[5]+12.36026325723227*fr[2]+12.36026325723227*fl[2]-8.254235307567585*fr[0]+8.254235307567585*fl[0]; 
  incr1[6] = 8.941145146908529*fr[10]-8.941145146908529*fl[10]-9.57421875*fr[6]-9.57421875*fl[6]+6.393703176377302*fr[4]-6.393703176377302*fl[4]; 
  incr1[7] = 7.556867535443651*fr[12]+7.556867535443651*fl[12]-11.54296875*fr[7]+11.54296875*fl[7]+12.36026325723226*fr[3]+12.36026325723226*fl[3]-8.254235307567583*fr[1]+8.254235307567583*fl[1]; 
  incr1[8] = 5.527677772592861*fr[11]+5.527677772592861*fl[11]-3.69140625*fr[8]+3.69140625*fl[8]; 
  incr1[9] = 3.307189138830738*fr[14]-3.307189138830738*fl[14]-8.94140625*fr[9]-8.94140625*fl[9]+13.65782481176513*fr[5]-13.65782481176513*fl[5]-14.62486071398016*fr[2]-14.62486071398016*fl[2]+9.766542925609524*fr[0]-9.766542925609524*fl[0]; 
  incr1[10] = (-11.54296875*fr[10])+11.54296875*fl[10]+12.36026325723226*fr[6]+12.36026325723226*fl[6]-8.254235307567585*fr[4]+8.254235307567585*fl[4]; 
  incr1[11] = (-9.57421875*fr[11])-9.57421875*fl[11]+6.393703176377302*fr[8]-6.393703176377302*fl[8]; 
  incr1[12] = (-8.94140625*fr[12])-8.94140625*fl[12]+13.65782481176513*fr[7]-13.65782481176513*fl[7]-14.62486071398016*fr[3]-14.62486071398016*fl[3]+9.76654292560952*fr[1]-9.76654292560952*fl[1]; 
  incr1[13] = 3.69140625*fl[13]-3.69140625*fr[13]; 
  incr1[14] = (-3.75*fr[14])+3.75*fl[14]+10.13860170372798*fr[9]+10.13860170372798*fl[9]-15.48651767229347*fr[5]+15.48651767229347*fl[5]+16.58303331777859*fr[2]+16.58303331777859*fl[2]-11.07421875*fr[0]+11.07421875*fl[0]; 

  double incr2[15]; 
  incr2[2] = 0.4837563778952138*fr[14]+0.4837563778952138*fl[14]-1.07020531715347*fr[9]+1.07020531715347*fl[9]+1.406982231239413*fr[5]+1.406982231239413*fl[5]-1.36328125*fr[2]+1.36328125*fl[2]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
  incr2[3] = (-1.07020531715347*fr[12])+1.07020531715347*fl[12]+1.406982231239413*fr[7]+1.406982231239413*fl[7]-1.36328125*fr[3]+1.36328125*fl[3]+0.8660254037844386*fr[1]+0.8660254037844386*fl[1]; 
  incr2[5] = (-1.873580395209785*fr[14])-1.873580395209785*fl[14]+4.144887370358018*fr[9]-4.144887370358018*fl[9]-5.44921875*fr[5]-5.44921875*fl[5]+5.27996557744683*fr[2]-5.27996557744683*fl[2]-3.354101966249685*fr[0]-3.354101966249685*fl[0]; 
  incr2[6] = 1.406982231239413*fr[10]+1.406982231239413*fl[10]-1.36328125*fr[6]+1.36328125*fl[6]+0.8660254037844387*fr[4]+0.8660254037844387*fl[4]; 
  incr2[7] = 4.144887370358018*fr[12]-4.144887370358018*fl[12]-5.44921875*fr[7]-5.44921875*fl[7]+5.279965577446831*fr[3]-5.279965577446831*fl[3]-3.354101966249684*fr[1]-3.354101966249684*fl[1]; 
  incr2[9] = 4.433700439244959*fr[14]+4.433700439244959*fl[14]-9.80859375*fr[9]+9.80859375*fl[9]+12.89520515222495*fr[5]+12.89520515222495*fl[5]-12.49467904327803*fr[2]+12.49467904327803*fl[2]+7.937253933193772*fr[0]+7.937253933193772*fl[0]; 
  incr2[10] = (-5.44921875*fr[10])-5.44921875*fl[10]+5.279965577446831*fr[6]-5.279965577446831*fl[6]-3.354101966249685*fr[4]-3.354101966249685*fl[4]; 
  incr2[11] = (-1.36328125*fr[11])+1.36328125*fl[11]+0.8660254037844387*fr[8]+0.8660254037844387*fl[8]; 
  incr2[12] = (-9.80859375*fr[12])+9.80859375*fl[12]+12.89520515222494*fr[7]+12.89520515222494*fl[7]-12.49467904327803*fr[3]+12.49467904327803*fl[3]+7.93725393319377*fr[1]+7.93725393319377*fl[1]; 
  incr2[14] = (-8.37890625*fr[14])-8.37890625*fl[14]+18.53649983840175*fr[9]-18.53649983840175*fl[9]-24.36964709853287*fr[5]-24.36964709853287*fl[5]+23.61272390006008*fr[2]-23.61272390006008*fl[2]-15.0*fr[0]-15.0*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr1[1]*rdxSq2nur; 
  outr[2] += incr2[2]*rdxSq2nur+incr1[2]*rdxSq2nur; 
  outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 
  outr[4] += incr1[4]*rdxSq2nur; 
  outr[5] += incr2[5]*rdxSq2nur+incr1[5]*rdxSq2nur; 
  outr[6] += incr2[6]*rdxSq2nur+incr1[6]*rdxSq2nur; 
  outr[7] += incr2[7]*rdxSq2nur+incr1[7]*rdxSq2nur; 
  outr[8] += incr1[8]*rdxSq2nur; 
  outr[9] += incr2[9]*rdxSq2nur+incr1[9]*rdxSq2nur; 
  outr[10] += incr2[10]*rdxSq2nur+incr1[10]*rdxSq2nur; 
  outr[11] += incr2[11]*rdxSq2nur+incr1[11]*rdxSq2nur; 
  outr[12] += incr2[12]*rdxSq2nur+incr1[12]*rdxSq2nur; 
  outr[13] += incr1[13]*rdxSq2nur; 
  outr[14] += incr2[14]*rdxSq2nur+incr1[14]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += -1.0*incr1[1]*rdxSq2nul; 
  outl[2] += incr1[2]*rdxSq2nul-1.0*incr2[2]*rdxSq2nul; 
  outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  outl[4] += -1.0*incr1[4]*rdxSq2nul; 
  outl[5] += incr2[5]*rdxSq2nul-1.0*incr1[5]*rdxSq2nul; 
  outl[6] += incr1[6]*rdxSq2nul-1.0*incr2[6]*rdxSq2nul; 
  outl[7] += incr2[7]*rdxSq2nul-1.0*incr1[7]*rdxSq2nul; 
  outl[8] += -1.0*incr1[8]*rdxSq2nul; 
  outl[9] += incr1[9]*rdxSq2nul-1.0*incr2[9]*rdxSq2nul; 
  outl[10] += incr2[10]*rdxSq2nul-1.0*incr1[10]*rdxSq2nul; 
  outl[11] += incr1[11]*rdxSq2nul-1.0*incr2[11]*rdxSq2nul; 
  outl[12] += incr1[12]*rdxSq2nul-1.0*incr2[12]*rdxSq2nul; 
  outl[13] += -1.0*incr1[13]*rdxSq2nul; 
  outl[14] += incr2[14]*rdxSq2nul-1.0*incr1[14]*rdxSq2nul; 

} 
