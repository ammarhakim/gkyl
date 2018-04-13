#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf1xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dxl = 4.0/(dxvl[0]*dxvl[0]); 
  double dxr = 4.0/(dxvr[0]*dxvr[0]); 
  double incr[2]; 

  incr[0] = 1.082531754730548*fr[1]+1.082531754730548*fl[1]-1.125*fr[0]+1.125*fl[0]; 
  incr[1] = (-1.875*fr[1])-1.875*fl[1]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 

  outr[0] += incr[0]*dxr*nu; 
  outr[1] += incr[1]*dxr*nu; 

  outl[0] += -1.0*incr[0]*dxl*nu; 
  outl[1] += incr[1]*dxl*nu; 

} 
void ConstDiffusionSurf1xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dxl = 4.0/(dxvl[0]*dxvl[0]); 
  double dxr = 4.0/(dxvr[0]*dxvr[0]); 
  double incr[3]; 

  incr[0] = (-1.341640786499874*fr[2])+1.341640786499874*fl[2]+2.381569860407206*fr[1]+2.381569860407206*fl[1]-1.875*fr[0]+1.875*fl[0]; 
  incr[1] = 2.32379000772445*fr[2]-2.32379000772445*fl[2]-4.125*fr[1]-4.125*fl[1]+3.247595264191645*fr[0]-3.247595264191645*fl[0]; 
  incr[2] = (-3.0*fr[2])+3.0*fl[2]+5.325352101035199*fr[1]+5.325352101035199*fl[1]-4.192627457812106*fr[0]+4.192627457812106*fl[0]; 

  outr[0] += incr[0]*dxr*nu; 
  outr[1] += incr[1]*dxr*nu; 
  outr[2] += incr[2]*dxr*nu; 

  outl[0] += -1.0*incr[0]*dxl*nu; 
  outl[1] += incr[1]*dxl*nu; 
  outl[2] += -1.0*incr[2]*dxl*nu; 

} 
void ConstDiffusionSurf1xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dxl = 4.0/(dxvl[0]*dxvl[0]); 
  double dxr = 4.0/(dxvr[0]*dxvr[0]); 
  double incr[4]; 

  incr[0] = 1.36421551976768*fr[3]+1.36421551976768*fl[3]-3.109532031210645*fr[2]+3.109532031210645*fl[2]+3.87005102316171*fr[1]+3.87005102316171*fl[1]-2.734375*fr[0]+2.734375*fl[0]; 
  incr[1] = (-2.362890592711605*fr[3])-2.362890592711605*fl[3]+5.385867465819689*fr[2]-5.385867465819689*fl[2]-6.703125*fr[1]-6.703125*fl[1]+4.736076426946148*fr[0]-4.736076426946148*fl[0]; 
  incr[2] = 3.05047863816074*fr[3]+3.05047863816074*fl[3]-6.953125*fr[2]+6.953125*fl[2]+8.653697164182198*fr[1]+8.653697164182198*fl[1]-6.114248375975989*fr[0]+6.114248375975989*fl[0]; 
  incr[3] = (-3.609375*fr[3])-3.609375*fl[3]+8.227048448372905*fr[2]-8.227048448372905*fl[2]-10.23919256841696*fr[1]-10.23919256841696*fl[1]+7.23447624119224*fr[0]-7.23447624119224*fl[0]; 

  outr[0] += incr[0]*dxr*nu; 
  outr[1] += incr[1]*dxr*nu; 
  outr[2] += incr[2]*dxr*nu; 
  outr[3] += incr[3]*dxr*nu; 

  outl[0] += -1.0*incr[0]*dxl*nu; 
  outl[1] += incr[1]*dxl*nu; 
  outl[2] += -1.0*incr[2]*dxl*nu; 
  outl[3] += incr[3]*dxl*nu; 

} 
void ConstDiffusionSurf1xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dxl = 4.0/(dxvl[0]*dxvl[0]); 
  double dxr = 4.0/(dxvr[0]*dxvr[0]); 
  double incr[5]; 

  incr[0] = (-1.25*fr[4])+1.25*fl[4]+3.379533901242661*fr[3]+3.379533901242661*fl[3]-5.162172557431155*fr[2]+5.162172557431155*fl[2]+5.527677772592862*fr[1]+5.527677772592862*fl[1]-3.69140625*fr[0]+3.69140625*fl[0]; 
  incr[1] = 2.165063509461096*fr[4]-2.165063509461096*fl[4]-5.853524422853748*fr[3]-5.853524422853748*fl[3]+8.941145146908527*fr[2]-8.941145146908527*fl[2]-9.57421875*fr[1]-9.57421875*fl[1]+6.393703176377298*fr[0]-6.393703176377298*fl[0]; 
  incr[2] = (-2.795084971874738*fr[4])+2.795084971874738*fl[4]+7.556867535443652*fr[3]+7.556867535443652*fl[3]-11.54296875*fr[2]+11.54296875*fl[2]+12.36026325723227*fr[1]+12.36026325723227*fl[1]-8.254235307567585*fr[0]+8.254235307567585*fl[0]; 
  incr[3] = 3.307189138830738*fr[4]-3.307189138830738*fl[4]-8.94140625*fr[3]-8.94140625*fl[3]+13.65782481176513*fr[2]-13.65782481176513*fl[2]-14.62486071398016*fr[1]-14.62486071398016*fl[1]+9.766542925609524*fr[0]-9.766542925609524*fl[0]; 
  incr[4] = (-3.75*fr[4])+3.75*fl[4]+10.13860170372798*fr[3]+10.13860170372798*fl[3]-15.48651767229347*fr[2]+15.48651767229347*fl[2]+16.58303331777859*fr[1]+16.58303331777859*fl[1]-11.07421875*fr[0]+11.07421875*fl[0]; 

  outr[0] += incr[0]*dxr*nu; 
  outr[1] += incr[1]*dxr*nu; 
  outr[2] += incr[2]*dxr*nu; 
  outr[3] += incr[3]*dxr*nu; 
  outr[4] += incr[4]*dxr*nu; 

  outl[0] += -1.0*incr[0]*dxl*nu; 
  outl[1] += incr[1]*dxl*nu; 
  outl[2] += -1.0*incr[2]*dxl*nu; 
  outl[3] += incr[3]*dxl*nu; 
  outl[4] += -1.0*incr[4]*dxl*nu; 

} 
