#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf3xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 
  double incr[4]; 

  incr[0] = 1.082531754730548*fr[1]+1.082531754730548*fl[1]-1.125*fr[0]+1.125*fl[0]; 
  incr[1] = (-1.875*fr[1])-1.875*fl[1]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 
  incr[2] = 1.125*fl[2]-1.125*fr[2]; 
  incr[3] = 1.125*fl[3]-1.125*fr[3]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += incr[1]*rdxSq2nul; 
  outl[2] += -1.0*incr[2]*rdxSq2nul; 
  outl[3] += -1.0*incr[3]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 
  double incr[10]; 

  incr[0] = (-1.341640786499874*fr[7])+1.341640786499874*fl[7]+2.381569860407206*fr[1]+2.381569860407206*fl[1]-1.875*fr[0]+1.875*fl[0]; 
  incr[1] = 2.32379000772445*fr[7]-2.32379000772445*fl[7]-4.125*fr[1]-4.125*fl[1]+3.247595264191645*fr[0]-3.247595264191645*fl[0]; 
  incr[2] = 2.381569860407206*fr[4]+2.381569860407206*fl[4]-1.875*fr[2]+1.875*fl[2]; 
  incr[3] = 2.381569860407206*fr[5]+2.381569860407206*fl[5]-1.875*fr[3]+1.875*fl[3]; 
  incr[4] = (-4.125*fr[4])-4.125*fl[4]+3.247595264191645*fr[2]-3.247595264191645*fl[2]; 
  incr[5] = (-4.125*fr[5])-4.125*fl[5]+3.247595264191645*fr[3]-3.247595264191645*fl[3]; 
  incr[6] = 1.875*fl[6]-1.875*fr[6]; 
  incr[7] = (-3.0*fr[7])+3.0*fl[7]+5.325352101035199*fr[1]+5.325352101035199*fl[1]-4.192627457812106*fr[0]+4.192627457812106*fl[0]; 
  incr[8] = 1.875*fl[8]-1.875*fr[8]; 
  incr[9] = 1.875*fl[9]-1.875*fr[9]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 
  outr[4] += incr[4]*rdxSq2nur; 
  outr[5] += incr[5]*rdxSq2nur; 
  outr[6] += incr[6]*rdxSq2nur; 
  outr[7] += incr[7]*rdxSq2nur; 
  outr[8] += incr[8]*rdxSq2nur; 
  outr[9] += incr[9]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += incr[1]*rdxSq2nul; 
  outl[2] += -1.0*incr[2]*rdxSq2nul; 
  outl[3] += -1.0*incr[3]*rdxSq2nul; 
  outl[4] += incr[4]*rdxSq2nul; 
  outl[5] += incr[5]*rdxSq2nul; 
  outl[6] += -1.0*incr[6]*rdxSq2nul; 
  outl[7] += -1.0*incr[7]*rdxSq2nul; 
  outl[8] += -1.0*incr[8]*rdxSq2nul; 
  outl[9] += -1.0*incr[9]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 
  double incr[20]; 

  incr[0] = 1.36421551976768*fr[17]+1.36421551976768*fl[17]-3.109532031210645*fr[7]+3.109532031210645*fl[7]+3.87005102316171*fr[1]+3.87005102316171*fl[1]-2.734375*fr[0]+2.734375*fl[0]; 
  incr[1] = (-2.362890592711605*fr[17])-2.362890592711605*fl[17]+5.385867465819689*fr[7]-5.385867465819689*fl[7]-6.703125*fr[1]-6.703125*fl[1]+4.736076426946148*fr[0]-4.736076426946148*fl[0]; 
  incr[2] = (-3.109532031210645*fr[11])+3.109532031210645*fl[11]+3.87005102316171*fr[4]+3.87005102316171*fl[4]-2.734375*fr[2]+2.734375*fl[2]; 
  incr[3] = (-3.109532031210645*fr[13])+3.109532031210645*fl[13]+3.87005102316171*fr[5]+3.87005102316171*fl[5]-2.734375*fr[3]+2.734375*fl[3]; 
  incr[4] = 5.38586746581969*fr[11]-5.38586746581969*fl[11]-6.703125*fr[4]-6.703125*fl[4]+4.736076426946148*fr[2]-4.736076426946148*fl[2]; 
  incr[5] = 5.38586746581969*fr[13]-5.38586746581969*fl[13]-6.703125*fr[5]-6.703125*fl[5]+4.736076426946148*fr[3]-4.736076426946148*fl[3]; 
  incr[6] = 3.87005102316171*fr[10]+3.87005102316171*fl[10]-2.734375*fr[6]+2.734375*fl[6]; 
  incr[7] = 3.05047863816074*fr[17]+3.05047863816074*fl[17]-6.953125*fr[7]+6.953125*fl[7]+8.653697164182198*fr[1]+8.653697164182198*fl[1]-6.114248375975989*fr[0]+6.114248375975989*fl[0]; 
  incr[8] = 3.87005102316171*fr[12]+3.87005102316171*fl[12]-2.734375*fr[8]+2.734375*fl[8]; 
  incr[9] = 3.87005102316171*fr[15]+3.87005102316171*fl[15]-2.734375*fr[9]+2.734375*fl[9]; 
  incr[10] = (-6.703125*fr[10])-6.703125*fl[10]+4.736076426946148*fr[6]-4.736076426946148*fl[6]; 
  incr[11] = (-6.953125*fr[11])+6.953125*fl[11]+8.653697164182198*fr[4]+8.653697164182198*fl[4]-6.114248375975989*fr[2]+6.114248375975989*fl[2]; 
  incr[12] = (-6.703125*fr[12])-6.703125*fl[12]+4.73607642694615*fr[8]-4.73607642694615*fl[8]; 
  incr[13] = (-6.953125*fr[13])+6.953125*fl[13]+8.653697164182198*fr[5]+8.653697164182198*fl[5]-6.114248375975989*fr[3]+6.114248375975989*fl[3]; 
  incr[14] = 2.734375*fl[14]-2.734375*fr[14]; 
  incr[15] = (-6.703125*fr[15])-6.703125*fl[15]+4.73607642694615*fr[9]-4.73607642694615*fl[9]; 
  incr[16] = 2.734375*fl[16]-2.734375*fr[16]; 
  incr[17] = (-3.609375*fr[17])-3.609375*fl[17]+8.227048448372905*fr[7]-8.227048448372905*fl[7]-10.23919256841696*fr[1]-10.23919256841696*fl[1]+7.23447624119224*fr[0]-7.23447624119224*fl[0]; 
  incr[18] = 2.734375*fl[18]-2.734375*fr[18]; 
  incr[19] = 2.734375*fl[19]-2.734375*fr[19]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 
  outr[4] += incr[4]*rdxSq2nur; 
  outr[5] += incr[5]*rdxSq2nur; 
  outr[6] += incr[6]*rdxSq2nur; 
  outr[7] += incr[7]*rdxSq2nur; 
  outr[8] += incr[8]*rdxSq2nur; 
  outr[9] += incr[9]*rdxSq2nur; 
  outr[10] += incr[10]*rdxSq2nur; 
  outr[11] += incr[11]*rdxSq2nur; 
  outr[12] += incr[12]*rdxSq2nur; 
  outr[13] += incr[13]*rdxSq2nur; 
  outr[14] += incr[14]*rdxSq2nur; 
  outr[15] += incr[15]*rdxSq2nur; 
  outr[16] += incr[16]*rdxSq2nur; 
  outr[17] += incr[17]*rdxSq2nur; 
  outr[18] += incr[18]*rdxSq2nur; 
  outr[19] += incr[19]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += incr[1]*rdxSq2nul; 
  outl[2] += -1.0*incr[2]*rdxSq2nul; 
  outl[3] += -1.0*incr[3]*rdxSq2nul; 
  outl[4] += incr[4]*rdxSq2nul; 
  outl[5] += incr[5]*rdxSq2nul; 
  outl[6] += -1.0*incr[6]*rdxSq2nul; 
  outl[7] += -1.0*incr[7]*rdxSq2nul; 
  outl[8] += -1.0*incr[8]*rdxSq2nul; 
  outl[9] += -1.0*incr[9]*rdxSq2nul; 
  outl[10] += incr[10]*rdxSq2nul; 
  outl[11] += -1.0*incr[11]*rdxSq2nul; 
  outl[12] += incr[12]*rdxSq2nul; 
  outl[13] += -1.0*incr[13]*rdxSq2nul; 
  outl[14] += -1.0*incr[14]*rdxSq2nul; 
  outl[15] += incr[15]*rdxSq2nul; 
  outl[16] += -1.0*incr[16]*rdxSq2nul; 
  outl[17] += incr[17]*rdxSq2nul; 
  outl[18] += -1.0*incr[18]*rdxSq2nul; 
  outl[19] += -1.0*incr[19]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 
  double incr[35]; 

  incr[0] = (-1.25*fr[32])+1.25*fl[32]+3.379533901242661*fr[17]+3.379533901242661*fl[17]-5.162172557431155*fr[7]+5.162172557431155*fl[7]+5.527677772592862*fr[1]+5.527677772592862*fl[1]-3.69140625*fr[0]+3.69140625*fl[0]; 
  incr[1] = 2.165063509461096*fr[32]-2.165063509461096*fl[32]-5.853524422853748*fr[17]-5.853524422853748*fl[17]+8.941145146908527*fr[7]-8.941145146908527*fl[7]-9.57421875*fr[1]-9.57421875*fl[1]+6.393703176377298*fr[0]-6.393703176377298*fl[0]; 
  incr[2] = 3.37953390124266*fr[26]+3.37953390124266*fl[26]-5.162172557431155*fr[11]+5.162172557431155*fl[11]+5.527677772592862*fr[4]+5.527677772592862*fl[4]-3.69140625*fr[2]+3.69140625*fl[2]; 
  incr[3] = 3.37953390124266*fr[28]+3.37953390124266*fl[28]-5.162172557431155*fr[13]+5.162172557431155*fl[13]+5.527677772592862*fr[5]+5.527677772592862*fl[5]-3.69140625*fr[3]+3.69140625*fl[3]; 
  incr[4] = (-5.853524422853749*fr[26])-5.853524422853749*fl[26]+8.941145146908529*fr[11]-8.941145146908529*fl[11]-9.57421875*fr[4]-9.57421875*fl[4]+6.393703176377298*fr[2]-6.393703176377298*fl[2]; 
  incr[5] = (-5.853524422853749*fr[28])-5.853524422853749*fl[28]+8.941145146908529*fr[13]-8.941145146908529*fl[13]-9.57421875*fr[5]-9.57421875*fl[5]+6.393703176377298*fr[3]-6.393703176377298*fl[3]; 
  incr[6] = (-5.162172557431155*fr[20])+5.162172557431155*fl[20]+5.527677772592862*fr[10]+5.527677772592862*fl[10]-3.69140625*fr[6]+3.69140625*fl[6]; 
  incr[7] = (-2.795084971874738*fr[32])+2.795084971874738*fl[32]+7.556867535443652*fr[17]+7.556867535443652*fl[17]-11.54296875*fr[7]+11.54296875*fl[7]+12.36026325723227*fr[1]+12.36026325723227*fl[1]-8.254235307567585*fr[0]+8.254235307567585*fl[0]; 
  incr[8] = (-5.162172557431155*fr[23])+5.162172557431155*fl[23]+5.527677772592862*fr[12]+5.527677772592862*fl[12]-3.69140625*fr[8]+3.69140625*fl[8]; 
  incr[9] = (-5.162172557431155*fr[24])+5.162172557431155*fl[24]+5.527677772592862*fr[15]+5.527677772592862*fl[15]-3.69140625*fr[9]+3.69140625*fl[9]; 
  incr[10] = 8.941145146908527*fr[20]-8.941145146908527*fl[20]-9.57421875*fr[10]-9.57421875*fl[10]+6.393703176377298*fr[6]-6.393703176377298*fl[6]; 
  incr[11] = 7.556867535443651*fr[26]+7.556867535443651*fl[26]-11.54296875*fr[11]+11.54296875*fl[11]+12.36026325723226*fr[4]+12.36026325723226*fl[4]-8.254235307567583*fr[2]+8.254235307567583*fl[2]; 
  incr[12] = 8.941145146908529*fr[23]-8.941145146908529*fl[23]-9.57421875*fr[12]-9.57421875*fl[12]+6.393703176377302*fr[8]-6.393703176377302*fl[8]; 
  incr[13] = 7.556867535443651*fr[28]+7.556867535443651*fl[28]-11.54296875*fr[13]+11.54296875*fl[13]+12.36026325723226*fr[5]+12.36026325723226*fl[5]-8.254235307567583*fr[3]+8.254235307567583*fl[3]; 
  incr[14] = 5.527677772592862*fr[21]+5.527677772592862*fl[21]-3.69140625*fr[14]+3.69140625*fl[14]; 
  incr[15] = 8.941145146908529*fr[24]-8.941145146908529*fl[24]-9.57421875*fr[15]-9.57421875*fl[15]+6.393703176377302*fr[9]-6.393703176377302*fl[9]; 
  incr[16] = 5.527677772592862*fr[22]+5.527677772592862*fl[22]-3.69140625*fr[16]+3.69140625*fl[16]; 
  incr[17] = 3.307189138830738*fr[32]-3.307189138830738*fl[32]-8.94140625*fr[17]-8.94140625*fl[17]+13.65782481176513*fr[7]-13.65782481176513*fl[7]-14.62486071398016*fr[1]-14.62486071398016*fl[1]+9.766542925609524*fr[0]-9.766542925609524*fl[0]; 
  incr[18] = 5.527677772592861*fr[27]+5.527677772592861*fl[27]-3.69140625*fr[18]+3.69140625*fl[18]; 
  incr[19] = 5.527677772592861*fr[30]+5.527677772592861*fl[30]-3.69140625*fr[19]+3.69140625*fl[19]; 
  incr[20] = (-11.54296875*fr[20])+11.54296875*fl[20]+12.36026325723227*fr[10]+12.36026325723227*fl[10]-8.254235307567585*fr[6]+8.254235307567585*fl[6]; 
  incr[21] = (-9.57421875*fr[21])-9.57421875*fl[21]+6.393703176377302*fr[14]-6.393703176377302*fl[14]; 
  incr[22] = (-9.57421875*fr[22])-9.57421875*fl[22]+6.393703176377302*fr[16]-6.393703176377302*fl[16]; 
  incr[23] = (-11.54296875*fr[23])+11.54296875*fl[23]+12.36026325723226*fr[12]+12.36026325723226*fl[12]-8.254235307567585*fr[8]+8.254235307567585*fl[8]; 
  incr[24] = (-11.54296875*fr[24])+11.54296875*fl[24]+12.36026325723226*fr[15]+12.36026325723226*fl[15]-8.254235307567585*fr[9]+8.254235307567585*fl[9]; 
  incr[25] = 3.69140625*fl[25]-3.69140625*fr[25]; 
  incr[26] = (-8.94140625*fr[26])-8.94140625*fl[26]+13.65782481176513*fr[11]-13.65782481176513*fl[11]-14.62486071398016*fr[4]-14.62486071398016*fl[4]+9.76654292560952*fr[2]-9.76654292560952*fl[2]; 
  incr[27] = (-9.57421875*fr[27])-9.57421875*fl[27]+6.393703176377302*fr[18]-6.393703176377302*fl[18]; 
  incr[28] = (-8.94140625*fr[28])-8.94140625*fl[28]+13.65782481176513*fr[13]-13.65782481176513*fl[13]-14.62486071398016*fr[5]-14.62486071398016*fl[5]+9.76654292560952*fr[3]-9.76654292560952*fl[3]; 
  incr[29] = 3.69140625*fl[29]-3.69140625*fr[29]; 
  incr[30] = (-9.57421875*fr[30])-9.57421875*fl[30]+6.393703176377302*fr[19]-6.393703176377302*fl[19]; 
  incr[31] = 3.69140625*fl[31]-3.69140625*fr[31]; 
  incr[32] = (-3.75*fr[32])+3.75*fl[32]+10.13860170372798*fr[17]+10.13860170372798*fl[17]-15.48651767229347*fr[7]+15.48651767229347*fl[7]+16.58303331777859*fr[1]+16.58303331777859*fl[1]-11.07421875*fr[0]+11.07421875*fl[0]; 
  incr[33] = 3.69140625*fl[33]-3.69140625*fr[33]; 
  incr[34] = 3.69140625*fl[34]-3.69140625*fr[34]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 
  outr[4] += incr[4]*rdxSq2nur; 
  outr[5] += incr[5]*rdxSq2nur; 
  outr[6] += incr[6]*rdxSq2nur; 
  outr[7] += incr[7]*rdxSq2nur; 
  outr[8] += incr[8]*rdxSq2nur; 
  outr[9] += incr[9]*rdxSq2nur; 
  outr[10] += incr[10]*rdxSq2nur; 
  outr[11] += incr[11]*rdxSq2nur; 
  outr[12] += incr[12]*rdxSq2nur; 
  outr[13] += incr[13]*rdxSq2nur; 
  outr[14] += incr[14]*rdxSq2nur; 
  outr[15] += incr[15]*rdxSq2nur; 
  outr[16] += incr[16]*rdxSq2nur; 
  outr[17] += incr[17]*rdxSq2nur; 
  outr[18] += incr[18]*rdxSq2nur; 
  outr[19] += incr[19]*rdxSq2nur; 
  outr[20] += incr[20]*rdxSq2nur; 
  outr[21] += incr[21]*rdxSq2nur; 
  outr[22] += incr[22]*rdxSq2nur; 
  outr[23] += incr[23]*rdxSq2nur; 
  outr[24] += incr[24]*rdxSq2nur; 
  outr[25] += incr[25]*rdxSq2nur; 
  outr[26] += incr[26]*rdxSq2nur; 
  outr[27] += incr[27]*rdxSq2nur; 
  outr[28] += incr[28]*rdxSq2nur; 
  outr[29] += incr[29]*rdxSq2nur; 
  outr[30] += incr[30]*rdxSq2nur; 
  outr[31] += incr[31]*rdxSq2nur; 
  outr[32] += incr[32]*rdxSq2nur; 
  outr[33] += incr[33]*rdxSq2nur; 
  outr[34] += incr[34]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += incr[1]*rdxSq2nul; 
  outl[2] += -1.0*incr[2]*rdxSq2nul; 
  outl[3] += -1.0*incr[3]*rdxSq2nul; 
  outl[4] += incr[4]*rdxSq2nul; 
  outl[5] += incr[5]*rdxSq2nul; 
  outl[6] += -1.0*incr[6]*rdxSq2nul; 
  outl[7] += -1.0*incr[7]*rdxSq2nul; 
  outl[8] += -1.0*incr[8]*rdxSq2nul; 
  outl[9] += -1.0*incr[9]*rdxSq2nul; 
  outl[10] += incr[10]*rdxSq2nul; 
  outl[11] += -1.0*incr[11]*rdxSq2nul; 
  outl[12] += incr[12]*rdxSq2nul; 
  outl[13] += -1.0*incr[13]*rdxSq2nul; 
  outl[14] += -1.0*incr[14]*rdxSq2nul; 
  outl[15] += incr[15]*rdxSq2nul; 
  outl[16] += -1.0*incr[16]*rdxSq2nul; 
  outl[17] += incr[17]*rdxSq2nul; 
  outl[18] += -1.0*incr[18]*rdxSq2nul; 
  outl[19] += -1.0*incr[19]*rdxSq2nul; 
  outl[20] += -1.0*incr[20]*rdxSq2nul; 
  outl[21] += incr[21]*rdxSq2nul; 
  outl[22] += incr[22]*rdxSq2nul; 
  outl[23] += -1.0*incr[23]*rdxSq2nul; 
  outl[24] += -1.0*incr[24]*rdxSq2nul; 
  outl[25] += -1.0*incr[25]*rdxSq2nul; 
  outl[26] += incr[26]*rdxSq2nul; 
  outl[27] += incr[27]*rdxSq2nul; 
  outl[28] += incr[28]*rdxSq2nul; 
  outl[29] += -1.0*incr[29]*rdxSq2nul; 
  outl[30] += incr[30]*rdxSq2nul; 
  outl[31] += -1.0*incr[31]*rdxSq2nul; 
  outl[32] += -1.0*incr[32]*rdxSq2nul; 
  outl[33] += -1.0*incr[33]*rdxSq2nul; 
  outl[34] += -1.0*incr[34]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 
  double incr[4]; 

  incr[0] = 1.082531754730548*fr[2]+1.082531754730548*fl[2]-1.125*fr[0]+1.125*fl[0]; 
  incr[1] = 1.125*fl[1]-1.125*fr[1]; 
  incr[2] = (-1.875*fr[2])-1.875*fl[2]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 
  incr[3] = 1.125*fl[3]-1.125*fr[3]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += -1.0*incr[1]*rdxSq2nul; 
  outl[2] += incr[2]*rdxSq2nul; 
  outl[3] += -1.0*incr[3]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 
  double incr[10]; 

  incr[0] = (-1.341640786499874*fr[8])+1.341640786499874*fl[8]+2.381569860407206*fr[2]+2.381569860407206*fl[2]-1.875*fr[0]+1.875*fl[0]; 
  incr[1] = 2.381569860407206*fr[4]+2.381569860407206*fl[4]-1.875*fr[1]+1.875*fl[1]; 
  incr[2] = 2.32379000772445*fr[8]-2.32379000772445*fl[8]-4.125*fr[2]-4.125*fl[2]+3.247595264191645*fr[0]-3.247595264191645*fl[0]; 
  incr[3] = 2.381569860407206*fr[6]+2.381569860407206*fl[6]-1.875*fr[3]+1.875*fl[3]; 
  incr[4] = (-4.125*fr[4])-4.125*fl[4]+3.247595264191645*fr[1]-3.247595264191645*fl[1]; 
  incr[5] = 1.875*fl[5]-1.875*fr[5]; 
  incr[6] = (-4.125*fr[6])-4.125*fl[6]+3.247595264191645*fr[3]-3.247595264191645*fl[3]; 
  incr[7] = 1.875*fl[7]-1.875*fr[7]; 
  incr[8] = (-3.0*fr[8])+3.0*fl[8]+5.325352101035199*fr[2]+5.325352101035199*fl[2]-4.192627457812106*fr[0]+4.192627457812106*fl[0]; 
  incr[9] = 1.875*fl[9]-1.875*fr[9]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 
  outr[4] += incr[4]*rdxSq2nur; 
  outr[5] += incr[5]*rdxSq2nur; 
  outr[6] += incr[6]*rdxSq2nur; 
  outr[7] += incr[7]*rdxSq2nur; 
  outr[8] += incr[8]*rdxSq2nur; 
  outr[9] += incr[9]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += -1.0*incr[1]*rdxSq2nul; 
  outl[2] += incr[2]*rdxSq2nul; 
  outl[3] += -1.0*incr[3]*rdxSq2nul; 
  outl[4] += incr[4]*rdxSq2nul; 
  outl[5] += -1.0*incr[5]*rdxSq2nul; 
  outl[6] += incr[6]*rdxSq2nul; 
  outl[7] += -1.0*incr[7]*rdxSq2nul; 
  outl[8] += -1.0*incr[8]*rdxSq2nul; 
  outl[9] += -1.0*incr[9]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 
  double incr[20]; 

  incr[0] = 1.36421551976768*fr[18]+1.36421551976768*fl[18]-3.109532031210645*fr[8]+3.109532031210645*fl[8]+3.87005102316171*fr[2]+3.87005102316171*fl[2]-2.734375*fr[0]+2.734375*fl[0]; 
  incr[1] = (-3.109532031210645*fr[12])+3.109532031210645*fl[12]+3.87005102316171*fr[4]+3.87005102316171*fl[4]-2.734375*fr[1]+2.734375*fl[1]; 
  incr[2] = (-2.362890592711605*fr[18])-2.362890592711605*fl[18]+5.385867465819689*fr[8]-5.385867465819689*fl[8]-6.703125*fr[2]-6.703125*fl[2]+4.736076426946148*fr[0]-4.736076426946148*fl[0]; 
  incr[3] = (-3.109532031210645*fr[14])+3.109532031210645*fl[14]+3.87005102316171*fr[6]+3.87005102316171*fl[6]-2.734375*fr[3]+2.734375*fl[3]; 
  incr[4] = 5.38586746581969*fr[12]-5.38586746581969*fl[12]-6.703125*fr[4]-6.703125*fl[4]+4.736076426946148*fr[1]-4.736076426946148*fl[1]; 
  incr[5] = 3.87005102316171*fr[10]+3.87005102316171*fl[10]-2.734375*fr[5]+2.734375*fl[5]; 
  incr[6] = 5.38586746581969*fr[14]-5.38586746581969*fl[14]-6.703125*fr[6]-6.703125*fl[6]+4.736076426946148*fr[3]-4.736076426946148*fl[3]; 
  incr[7] = 3.87005102316171*fr[11]+3.87005102316171*fl[11]-2.734375*fr[7]+2.734375*fl[7]; 
  incr[8] = 3.05047863816074*fr[18]+3.05047863816074*fl[18]-6.953125*fr[8]+6.953125*fl[8]+8.653697164182198*fr[2]+8.653697164182198*fl[2]-6.114248375975989*fr[0]+6.114248375975989*fl[0]; 
  incr[9] = 3.87005102316171*fr[16]+3.87005102316171*fl[16]-2.734375*fr[9]+2.734375*fl[9]; 
  incr[10] = (-6.703125*fr[10])-6.703125*fl[10]+4.736076426946148*fr[5]-4.736076426946148*fl[5]; 
  incr[11] = (-6.703125*fr[11])-6.703125*fl[11]+4.73607642694615*fr[7]-4.73607642694615*fl[7]; 
  incr[12] = (-6.953125*fr[12])+6.953125*fl[12]+8.653697164182198*fr[4]+8.653697164182198*fl[4]-6.114248375975989*fr[1]+6.114248375975989*fl[1]; 
  incr[13] = 2.734375*fl[13]-2.734375*fr[13]; 
  incr[14] = (-6.953125*fr[14])+6.953125*fl[14]+8.653697164182198*fr[6]+8.653697164182198*fl[6]-6.114248375975989*fr[3]+6.114248375975989*fl[3]; 
  incr[15] = 2.734375*fl[15]-2.734375*fr[15]; 
  incr[16] = (-6.703125*fr[16])-6.703125*fl[16]+4.73607642694615*fr[9]-4.73607642694615*fl[9]; 
  incr[17] = 2.734375*fl[17]-2.734375*fr[17]; 
  incr[18] = (-3.609375*fr[18])-3.609375*fl[18]+8.227048448372905*fr[8]-8.227048448372905*fl[8]-10.23919256841696*fr[2]-10.23919256841696*fl[2]+7.23447624119224*fr[0]-7.23447624119224*fl[0]; 
  incr[19] = 2.734375*fl[19]-2.734375*fr[19]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 
  outr[4] += incr[4]*rdxSq2nur; 
  outr[5] += incr[5]*rdxSq2nur; 
  outr[6] += incr[6]*rdxSq2nur; 
  outr[7] += incr[7]*rdxSq2nur; 
  outr[8] += incr[8]*rdxSq2nur; 
  outr[9] += incr[9]*rdxSq2nur; 
  outr[10] += incr[10]*rdxSq2nur; 
  outr[11] += incr[11]*rdxSq2nur; 
  outr[12] += incr[12]*rdxSq2nur; 
  outr[13] += incr[13]*rdxSq2nur; 
  outr[14] += incr[14]*rdxSq2nur; 
  outr[15] += incr[15]*rdxSq2nur; 
  outr[16] += incr[16]*rdxSq2nur; 
  outr[17] += incr[17]*rdxSq2nur; 
  outr[18] += incr[18]*rdxSq2nur; 
  outr[19] += incr[19]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += -1.0*incr[1]*rdxSq2nul; 
  outl[2] += incr[2]*rdxSq2nul; 
  outl[3] += -1.0*incr[3]*rdxSq2nul; 
  outl[4] += incr[4]*rdxSq2nul; 
  outl[5] += -1.0*incr[5]*rdxSq2nul; 
  outl[6] += incr[6]*rdxSq2nul; 
  outl[7] += -1.0*incr[7]*rdxSq2nul; 
  outl[8] += -1.0*incr[8]*rdxSq2nul; 
  outl[9] += -1.0*incr[9]*rdxSq2nul; 
  outl[10] += incr[10]*rdxSq2nul; 
  outl[11] += incr[11]*rdxSq2nul; 
  outl[12] += -1.0*incr[12]*rdxSq2nul; 
  outl[13] += -1.0*incr[13]*rdxSq2nul; 
  outl[14] += -1.0*incr[14]*rdxSq2nul; 
  outl[15] += -1.0*incr[15]*rdxSq2nul; 
  outl[16] += incr[16]*rdxSq2nul; 
  outl[17] += -1.0*incr[17]*rdxSq2nul; 
  outl[18] += incr[18]*rdxSq2nul; 
  outl[19] += -1.0*incr[19]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 
  double incr[35]; 

  incr[0] = (-1.25*fr[33])+1.25*fl[33]+3.379533901242661*fr[18]+3.379533901242661*fl[18]-5.162172557431155*fr[8]+5.162172557431155*fl[8]+5.527677772592862*fr[2]+5.527677772592862*fl[2]-3.69140625*fr[0]+3.69140625*fl[0]; 
  incr[1] = 3.37953390124266*fr[27]+3.37953390124266*fl[27]-5.162172557431155*fr[12]+5.162172557431155*fl[12]+5.527677772592862*fr[4]+5.527677772592862*fl[4]-3.69140625*fr[1]+3.69140625*fl[1]; 
  incr[2] = 2.165063509461096*fr[33]-2.165063509461096*fl[33]-5.853524422853748*fr[18]-5.853524422853748*fl[18]+8.941145146908527*fr[8]-8.941145146908527*fl[8]-9.57421875*fr[2]-9.57421875*fl[2]+6.393703176377298*fr[0]-6.393703176377298*fl[0]; 
  incr[3] = 3.37953390124266*fr[29]+3.37953390124266*fl[29]-5.162172557431155*fr[14]+5.162172557431155*fl[14]+5.527677772592862*fr[6]+5.527677772592862*fl[6]-3.69140625*fr[3]+3.69140625*fl[3]; 
  incr[4] = (-5.853524422853749*fr[27])-5.853524422853749*fl[27]+8.941145146908529*fr[12]-8.941145146908529*fl[12]-9.57421875*fr[4]-9.57421875*fl[4]+6.393703176377298*fr[1]-6.393703176377298*fl[1]; 
  incr[5] = (-5.162172557431155*fr[21])+5.162172557431155*fl[21]+5.527677772592862*fr[10]+5.527677772592862*fl[10]-3.69140625*fr[5]+3.69140625*fl[5]; 
  incr[6] = (-5.853524422853749*fr[29])-5.853524422853749*fl[29]+8.941145146908529*fr[14]-8.941145146908529*fl[14]-9.57421875*fr[6]-9.57421875*fl[6]+6.393703176377298*fr[3]-6.393703176377298*fl[3]; 
  incr[7] = (-5.162172557431155*fr[23])+5.162172557431155*fl[23]+5.527677772592862*fr[11]+5.527677772592862*fl[11]-3.69140625*fr[7]+3.69140625*fl[7]; 
  incr[8] = (-2.795084971874738*fr[33])+2.795084971874738*fl[33]+7.556867535443652*fr[18]+7.556867535443652*fl[18]-11.54296875*fr[8]+11.54296875*fl[8]+12.36026325723227*fr[2]+12.36026325723227*fl[2]-8.254235307567585*fr[0]+8.254235307567585*fl[0]; 
  incr[9] = (-5.162172557431155*fr[25])+5.162172557431155*fl[25]+5.527677772592862*fr[16]+5.527677772592862*fl[16]-3.69140625*fr[9]+3.69140625*fl[9]; 
  incr[10] = 8.941145146908527*fr[21]-8.941145146908527*fl[21]-9.57421875*fr[10]-9.57421875*fl[10]+6.393703176377298*fr[5]-6.393703176377298*fl[5]; 
  incr[11] = 8.941145146908529*fr[23]-8.941145146908529*fl[23]-9.57421875*fr[11]-9.57421875*fl[11]+6.393703176377302*fr[7]-6.393703176377302*fl[7]; 
  incr[12] = 7.556867535443651*fr[27]+7.556867535443651*fl[27]-11.54296875*fr[12]+11.54296875*fl[12]+12.36026325723226*fr[4]+12.36026325723226*fl[4]-8.254235307567583*fr[1]+8.254235307567583*fl[1]; 
  incr[13] = 5.527677772592862*fr[20]+5.527677772592862*fl[20]-3.69140625*fr[13]+3.69140625*fl[13]; 
  incr[14] = 7.556867535443651*fr[29]+7.556867535443651*fl[29]-11.54296875*fr[14]+11.54296875*fl[14]+12.36026325723226*fr[6]+12.36026325723226*fl[6]-8.254235307567583*fr[3]+8.254235307567583*fl[3]; 
  incr[15] = 5.527677772592862*fr[22]+5.527677772592862*fl[22]-3.69140625*fr[15]+3.69140625*fl[15]; 
  incr[16] = 8.941145146908529*fr[25]-8.941145146908529*fl[25]-9.57421875*fr[16]-9.57421875*fl[16]+6.393703176377302*fr[9]-6.393703176377302*fl[9]; 
  incr[17] = 5.527677772592861*fr[26]+5.527677772592861*fl[26]-3.69140625*fr[17]+3.69140625*fl[17]; 
  incr[18] = 3.307189138830738*fr[33]-3.307189138830738*fl[33]-8.94140625*fr[18]-8.94140625*fl[18]+13.65782481176513*fr[8]-13.65782481176513*fl[8]-14.62486071398016*fr[2]-14.62486071398016*fl[2]+9.766542925609524*fr[0]-9.766542925609524*fl[0]; 
  incr[19] = 5.527677772592861*fr[31]+5.527677772592861*fl[31]-3.69140625*fr[19]+3.69140625*fl[19]; 
  incr[20] = (-9.57421875*fr[20])-9.57421875*fl[20]+6.393703176377302*fr[13]-6.393703176377302*fl[13]; 
  incr[21] = (-11.54296875*fr[21])+11.54296875*fl[21]+12.36026325723227*fr[10]+12.36026325723227*fl[10]-8.254235307567585*fr[5]+8.254235307567585*fl[5]; 
  incr[22] = (-9.57421875*fr[22])-9.57421875*fl[22]+6.393703176377302*fr[15]-6.393703176377302*fl[15]; 
  incr[23] = (-11.54296875*fr[23])+11.54296875*fl[23]+12.36026325723226*fr[11]+12.36026325723226*fl[11]-8.254235307567585*fr[7]+8.254235307567585*fl[7]; 
  incr[24] = 3.69140625*fl[24]-3.69140625*fr[24]; 
  incr[25] = (-11.54296875*fr[25])+11.54296875*fl[25]+12.36026325723226*fr[16]+12.36026325723226*fl[16]-8.254235307567585*fr[9]+8.254235307567585*fl[9]; 
  incr[26] = (-9.57421875*fr[26])-9.57421875*fl[26]+6.393703176377302*fr[17]-6.393703176377302*fl[17]; 
  incr[27] = (-8.94140625*fr[27])-8.94140625*fl[27]+13.65782481176513*fr[12]-13.65782481176513*fl[12]-14.62486071398016*fr[4]-14.62486071398016*fl[4]+9.76654292560952*fr[1]-9.76654292560952*fl[1]; 
  incr[28] = 3.69140625*fl[28]-3.69140625*fr[28]; 
  incr[29] = (-8.94140625*fr[29])-8.94140625*fl[29]+13.65782481176513*fr[14]-13.65782481176513*fl[14]-14.62486071398016*fr[6]-14.62486071398016*fl[6]+9.76654292560952*fr[3]-9.76654292560952*fl[3]; 
  incr[30] = 3.69140625*fl[30]-3.69140625*fr[30]; 
  incr[31] = (-9.57421875*fr[31])-9.57421875*fl[31]+6.393703176377302*fr[19]-6.393703176377302*fl[19]; 
  incr[32] = 3.69140625*fl[32]-3.69140625*fr[32]; 
  incr[33] = (-3.75*fr[33])+3.75*fl[33]+10.13860170372798*fr[18]+10.13860170372798*fl[18]-15.48651767229347*fr[8]+15.48651767229347*fl[8]+16.58303331777859*fr[2]+16.58303331777859*fl[2]-11.07421875*fr[0]+11.07421875*fl[0]; 
  incr[34] = 3.69140625*fl[34]-3.69140625*fr[34]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 
  outr[4] += incr[4]*rdxSq2nur; 
  outr[5] += incr[5]*rdxSq2nur; 
  outr[6] += incr[6]*rdxSq2nur; 
  outr[7] += incr[7]*rdxSq2nur; 
  outr[8] += incr[8]*rdxSq2nur; 
  outr[9] += incr[9]*rdxSq2nur; 
  outr[10] += incr[10]*rdxSq2nur; 
  outr[11] += incr[11]*rdxSq2nur; 
  outr[12] += incr[12]*rdxSq2nur; 
  outr[13] += incr[13]*rdxSq2nur; 
  outr[14] += incr[14]*rdxSq2nur; 
  outr[15] += incr[15]*rdxSq2nur; 
  outr[16] += incr[16]*rdxSq2nur; 
  outr[17] += incr[17]*rdxSq2nur; 
  outr[18] += incr[18]*rdxSq2nur; 
  outr[19] += incr[19]*rdxSq2nur; 
  outr[20] += incr[20]*rdxSq2nur; 
  outr[21] += incr[21]*rdxSq2nur; 
  outr[22] += incr[22]*rdxSq2nur; 
  outr[23] += incr[23]*rdxSq2nur; 
  outr[24] += incr[24]*rdxSq2nur; 
  outr[25] += incr[25]*rdxSq2nur; 
  outr[26] += incr[26]*rdxSq2nur; 
  outr[27] += incr[27]*rdxSq2nur; 
  outr[28] += incr[28]*rdxSq2nur; 
  outr[29] += incr[29]*rdxSq2nur; 
  outr[30] += incr[30]*rdxSq2nur; 
  outr[31] += incr[31]*rdxSq2nur; 
  outr[32] += incr[32]*rdxSq2nur; 
  outr[33] += incr[33]*rdxSq2nur; 
  outr[34] += incr[34]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += -1.0*incr[1]*rdxSq2nul; 
  outl[2] += incr[2]*rdxSq2nul; 
  outl[3] += -1.0*incr[3]*rdxSq2nul; 
  outl[4] += incr[4]*rdxSq2nul; 
  outl[5] += -1.0*incr[5]*rdxSq2nul; 
  outl[6] += incr[6]*rdxSq2nul; 
  outl[7] += -1.0*incr[7]*rdxSq2nul; 
  outl[8] += -1.0*incr[8]*rdxSq2nul; 
  outl[9] += -1.0*incr[9]*rdxSq2nul; 
  outl[10] += incr[10]*rdxSq2nul; 
  outl[11] += incr[11]*rdxSq2nul; 
  outl[12] += -1.0*incr[12]*rdxSq2nul; 
  outl[13] += -1.0*incr[13]*rdxSq2nul; 
  outl[14] += -1.0*incr[14]*rdxSq2nul; 
  outl[15] += -1.0*incr[15]*rdxSq2nul; 
  outl[16] += incr[16]*rdxSq2nul; 
  outl[17] += -1.0*incr[17]*rdxSq2nul; 
  outl[18] += incr[18]*rdxSq2nul; 
  outl[19] += -1.0*incr[19]*rdxSq2nul; 
  outl[20] += incr[20]*rdxSq2nul; 
  outl[21] += -1.0*incr[21]*rdxSq2nul; 
  outl[22] += incr[22]*rdxSq2nul; 
  outl[23] += -1.0*incr[23]*rdxSq2nul; 
  outl[24] += -1.0*incr[24]*rdxSq2nul; 
  outl[25] += -1.0*incr[25]*rdxSq2nul; 
  outl[26] += incr[26]*rdxSq2nul; 
  outl[27] += incr[27]*rdxSq2nul; 
  outl[28] += -1.0*incr[28]*rdxSq2nul; 
  outl[29] += incr[29]*rdxSq2nul; 
  outl[30] += -1.0*incr[30]*rdxSq2nul; 
  outl[31] += incr[31]*rdxSq2nul; 
  outl[32] += -1.0*incr[32]*rdxSq2nul; 
  outl[33] += -1.0*incr[33]*rdxSq2nul; 
  outl[34] += -1.0*incr[34]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[2]/(dxvl[2]*dxvl[2]); 
  double rdxSq2nur = 2.0*nu[2]/(dxvr[2]*dxvr[2]); 
  double incr[4]; 

  incr[0] = 1.082531754730548*fr[3]+1.082531754730548*fl[3]-1.125*fr[0]+1.125*fl[0]; 
  incr[1] = 1.125*fl[1]-1.125*fr[1]; 
  incr[2] = 1.125*fl[2]-1.125*fr[2]; 
  incr[3] = (-1.875*fr[3])-1.875*fl[3]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += -1.0*incr[1]*rdxSq2nul; 
  outl[2] += -1.0*incr[2]*rdxSq2nul; 
  outl[3] += incr[3]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[2]/(dxvl[2]*dxvl[2]); 
  double rdxSq2nur = 2.0*nu[2]/(dxvr[2]*dxvr[2]); 
  double incr[10]; 

  incr[0] = (-1.341640786499874*fr[9])+1.341640786499874*fl[9]+2.381569860407206*fr[3]+2.381569860407206*fl[3]-1.875*fr[0]+1.875*fl[0]; 
  incr[1] = 2.381569860407206*fr[5]+2.381569860407206*fl[5]-1.875*fr[1]+1.875*fl[1]; 
  incr[2] = 2.381569860407206*fr[6]+2.381569860407206*fl[6]-1.875*fr[2]+1.875*fl[2]; 
  incr[3] = 2.32379000772445*fr[9]-2.32379000772445*fl[9]-4.125*fr[3]-4.125*fl[3]+3.247595264191645*fr[0]-3.247595264191645*fl[0]; 
  incr[4] = 1.875*fl[4]-1.875*fr[4]; 
  incr[5] = (-4.125*fr[5])-4.125*fl[5]+3.247595264191645*fr[1]-3.247595264191645*fl[1]; 
  incr[6] = (-4.125*fr[6])-4.125*fl[6]+3.247595264191645*fr[2]-3.247595264191645*fl[2]; 
  incr[7] = 1.875*fl[7]-1.875*fr[7]; 
  incr[8] = 1.875*fl[8]-1.875*fr[8]; 
  incr[9] = (-3.0*fr[9])+3.0*fl[9]+5.325352101035199*fr[3]+5.325352101035199*fl[3]-4.192627457812106*fr[0]+4.192627457812106*fl[0]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 
  outr[4] += incr[4]*rdxSq2nur; 
  outr[5] += incr[5]*rdxSq2nur; 
  outr[6] += incr[6]*rdxSq2nur; 
  outr[7] += incr[7]*rdxSq2nur; 
  outr[8] += incr[8]*rdxSq2nur; 
  outr[9] += incr[9]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += -1.0*incr[1]*rdxSq2nul; 
  outl[2] += -1.0*incr[2]*rdxSq2nul; 
  outl[3] += incr[3]*rdxSq2nul; 
  outl[4] += -1.0*incr[4]*rdxSq2nul; 
  outl[5] += incr[5]*rdxSq2nul; 
  outl[6] += incr[6]*rdxSq2nul; 
  outl[7] += -1.0*incr[7]*rdxSq2nul; 
  outl[8] += -1.0*incr[8]*rdxSq2nul; 
  outl[9] += -1.0*incr[9]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_Z_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[2]/(dxvl[2]*dxvl[2]); 
  double rdxSq2nur = 2.0*nu[2]/(dxvr[2]*dxvr[2]); 
  double incr[20]; 

  incr[0] = 1.36421551976768*fr[19]+1.36421551976768*fl[19]-3.109532031210645*fr[9]+3.109532031210645*fl[9]+3.87005102316171*fr[3]+3.87005102316171*fl[3]-2.734375*fr[0]+2.734375*fl[0]; 
  incr[1] = (-3.109532031210645*fr[15])+3.109532031210645*fl[15]+3.87005102316171*fr[5]+3.87005102316171*fl[5]-2.734375*fr[1]+2.734375*fl[1]; 
  incr[2] = (-3.109532031210645*fr[16])+3.109532031210645*fl[16]+3.87005102316171*fr[6]+3.87005102316171*fl[6]-2.734375*fr[2]+2.734375*fl[2]; 
  incr[3] = (-2.362890592711605*fr[19])-2.362890592711605*fl[19]+5.385867465819689*fr[9]-5.385867465819689*fl[9]-6.703125*fr[3]-6.703125*fl[3]+4.736076426946148*fr[0]-4.736076426946148*fl[0]; 
  incr[4] = 3.87005102316171*fr[10]+3.87005102316171*fl[10]-2.734375*fr[4]+2.734375*fl[4]; 
  incr[5] = 5.38586746581969*fr[15]-5.38586746581969*fl[15]-6.703125*fr[5]-6.703125*fl[5]+4.736076426946148*fr[1]-4.736076426946148*fl[1]; 
  incr[6] = 5.38586746581969*fr[16]-5.38586746581969*fl[16]-6.703125*fr[6]-6.703125*fl[6]+4.736076426946148*fr[2]-4.736076426946148*fl[2]; 
  incr[7] = 3.87005102316171*fr[13]+3.87005102316171*fl[13]-2.734375*fr[7]+2.734375*fl[7]; 
  incr[8] = 3.87005102316171*fr[14]+3.87005102316171*fl[14]-2.734375*fr[8]+2.734375*fl[8]; 
  incr[9] = 3.05047863816074*fr[19]+3.05047863816074*fl[19]-6.953125*fr[9]+6.953125*fl[9]+8.653697164182198*fr[3]+8.653697164182198*fl[3]-6.114248375975989*fr[0]+6.114248375975989*fl[0]; 
  incr[10] = (-6.703125*fr[10])-6.703125*fl[10]+4.736076426946148*fr[4]-4.736076426946148*fl[4]; 
  incr[11] = 2.734375*fl[11]-2.734375*fr[11]; 
  incr[12] = 2.734375*fl[12]-2.734375*fr[12]; 
  incr[13] = (-6.703125*fr[13])-6.703125*fl[13]+4.73607642694615*fr[7]-4.73607642694615*fl[7]; 
  incr[14] = (-6.703125*fr[14])-6.703125*fl[14]+4.73607642694615*fr[8]-4.73607642694615*fl[8]; 
  incr[15] = (-6.953125*fr[15])+6.953125*fl[15]+8.653697164182198*fr[5]+8.653697164182198*fl[5]-6.114248375975989*fr[1]+6.114248375975989*fl[1]; 
  incr[16] = (-6.953125*fr[16])+6.953125*fl[16]+8.653697164182198*fr[6]+8.653697164182198*fl[6]-6.114248375975989*fr[2]+6.114248375975989*fl[2]; 
  incr[17] = 2.734375*fl[17]-2.734375*fr[17]; 
  incr[18] = 2.734375*fl[18]-2.734375*fr[18]; 
  incr[19] = (-3.609375*fr[19])-3.609375*fl[19]+8.227048448372905*fr[9]-8.227048448372905*fl[9]-10.23919256841696*fr[3]-10.23919256841696*fl[3]+7.23447624119224*fr[0]-7.23447624119224*fl[0]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 
  outr[4] += incr[4]*rdxSq2nur; 
  outr[5] += incr[5]*rdxSq2nur; 
  outr[6] += incr[6]*rdxSq2nur; 
  outr[7] += incr[7]*rdxSq2nur; 
  outr[8] += incr[8]*rdxSq2nur; 
  outr[9] += incr[9]*rdxSq2nur; 
  outr[10] += incr[10]*rdxSq2nur; 
  outr[11] += incr[11]*rdxSq2nur; 
  outr[12] += incr[12]*rdxSq2nur; 
  outr[13] += incr[13]*rdxSq2nur; 
  outr[14] += incr[14]*rdxSq2nur; 
  outr[15] += incr[15]*rdxSq2nur; 
  outr[16] += incr[16]*rdxSq2nur; 
  outr[17] += incr[17]*rdxSq2nur; 
  outr[18] += incr[18]*rdxSq2nur; 
  outr[19] += incr[19]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += -1.0*incr[1]*rdxSq2nul; 
  outl[2] += -1.0*incr[2]*rdxSq2nul; 
  outl[3] += incr[3]*rdxSq2nul; 
  outl[4] += -1.0*incr[4]*rdxSq2nul; 
  outl[5] += incr[5]*rdxSq2nul; 
  outl[6] += incr[6]*rdxSq2nul; 
  outl[7] += -1.0*incr[7]*rdxSq2nul; 
  outl[8] += -1.0*incr[8]*rdxSq2nul; 
  outl[9] += -1.0*incr[9]*rdxSq2nul; 
  outl[10] += incr[10]*rdxSq2nul; 
  outl[11] += -1.0*incr[11]*rdxSq2nul; 
  outl[12] += -1.0*incr[12]*rdxSq2nul; 
  outl[13] += incr[13]*rdxSq2nul; 
  outl[14] += incr[14]*rdxSq2nul; 
  outl[15] += -1.0*incr[15]*rdxSq2nul; 
  outl[16] += -1.0*incr[16]*rdxSq2nul; 
  outl[17] += -1.0*incr[17]*rdxSq2nul; 
  outl[18] += -1.0*incr[18]*rdxSq2nul; 
  outl[19] += incr[19]*rdxSq2nul; 

} 
void ConstDiffusionSurf3xMax_Z_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[2]/(dxvl[2]*dxvl[2]); 
  double rdxSq2nur = 2.0*nu[2]/(dxvr[2]*dxvr[2]); 
  double incr[35]; 

  incr[0] = (-1.25*fr[34])+1.25*fl[34]+3.379533901242661*fr[19]+3.379533901242661*fl[19]-5.162172557431155*fr[9]+5.162172557431155*fl[9]+5.527677772592862*fr[3]+5.527677772592862*fl[3]-3.69140625*fr[0]+3.69140625*fl[0]; 
  incr[1] = 3.37953390124266*fr[30]+3.37953390124266*fl[30]-5.162172557431155*fr[15]+5.162172557431155*fl[15]+5.527677772592862*fr[5]+5.527677772592862*fl[5]-3.69140625*fr[1]+3.69140625*fl[1]; 
  incr[2] = 3.37953390124266*fr[31]+3.37953390124266*fl[31]-5.162172557431155*fr[16]+5.162172557431155*fl[16]+5.527677772592862*fr[6]+5.527677772592862*fl[6]-3.69140625*fr[2]+3.69140625*fl[2]; 
  incr[3] = 2.165063509461096*fr[34]-2.165063509461096*fl[34]-5.853524422853748*fr[19]-5.853524422853748*fl[19]+8.941145146908527*fr[9]-8.941145146908527*fl[9]-9.57421875*fr[3]-9.57421875*fl[3]+6.393703176377298*fr[0]-6.393703176377298*fl[0]; 
  incr[4] = (-5.162172557431155*fr[22])+5.162172557431155*fl[22]+5.527677772592862*fr[10]+5.527677772592862*fl[10]-3.69140625*fr[4]+3.69140625*fl[4]; 
  incr[5] = (-5.853524422853749*fr[30])-5.853524422853749*fl[30]+8.941145146908529*fr[15]-8.941145146908529*fl[15]-9.57421875*fr[5]-9.57421875*fl[5]+6.393703176377298*fr[1]-6.393703176377298*fl[1]; 
  incr[6] = (-5.853524422853749*fr[31])-5.853524422853749*fl[31]+8.941145146908529*fr[16]-8.941145146908529*fl[16]-9.57421875*fr[6]-9.57421875*fl[6]+6.393703176377298*fr[2]-6.393703176377298*fl[2]; 
  incr[7] = (-5.162172557431155*fr[24])+5.162172557431155*fl[24]+5.527677772592862*fr[13]+5.527677772592862*fl[13]-3.69140625*fr[7]+3.69140625*fl[7]; 
  incr[8] = (-5.162172557431155*fr[25])+5.162172557431155*fl[25]+5.527677772592862*fr[14]+5.527677772592862*fl[14]-3.69140625*fr[8]+3.69140625*fl[8]; 
  incr[9] = (-2.795084971874738*fr[34])+2.795084971874738*fl[34]+7.556867535443652*fr[19]+7.556867535443652*fl[19]-11.54296875*fr[9]+11.54296875*fl[9]+12.36026325723227*fr[3]+12.36026325723227*fl[3]-8.254235307567585*fr[0]+8.254235307567585*fl[0]; 
  incr[10] = 8.941145146908527*fr[22]-8.941145146908527*fl[22]-9.57421875*fr[10]-9.57421875*fl[10]+6.393703176377298*fr[4]-6.393703176377298*fl[4]; 
  incr[11] = 5.527677772592862*fr[20]+5.527677772592862*fl[20]-3.69140625*fr[11]+3.69140625*fl[11]; 
  incr[12] = 5.527677772592862*fr[21]+5.527677772592862*fl[21]-3.69140625*fr[12]+3.69140625*fl[12]; 
  incr[13] = 8.941145146908529*fr[24]-8.941145146908529*fl[24]-9.57421875*fr[13]-9.57421875*fl[13]+6.393703176377302*fr[7]-6.393703176377302*fl[7]; 
  incr[14] = 8.941145146908529*fr[25]-8.941145146908529*fl[25]-9.57421875*fr[14]-9.57421875*fl[14]+6.393703176377302*fr[8]-6.393703176377302*fl[8]; 
  incr[15] = 7.556867535443651*fr[30]+7.556867535443651*fl[30]-11.54296875*fr[15]+11.54296875*fl[15]+12.36026325723226*fr[5]+12.36026325723226*fl[5]-8.254235307567583*fr[1]+8.254235307567583*fl[1]; 
  incr[16] = 7.556867535443651*fr[31]+7.556867535443651*fl[31]-11.54296875*fr[16]+11.54296875*fl[16]+12.36026325723226*fr[6]+12.36026325723226*fl[6]-8.254235307567583*fr[2]+8.254235307567583*fl[2]; 
  incr[17] = 5.527677772592861*fr[28]+5.527677772592861*fl[28]-3.69140625*fr[17]+3.69140625*fl[17]; 
  incr[18] = 5.527677772592861*fr[29]+5.527677772592861*fl[29]-3.69140625*fr[18]+3.69140625*fl[18]; 
  incr[19] = 3.307189138830738*fr[34]-3.307189138830738*fl[34]-8.94140625*fr[19]-8.94140625*fl[19]+13.65782481176513*fr[9]-13.65782481176513*fl[9]-14.62486071398016*fr[3]-14.62486071398016*fl[3]+9.766542925609524*fr[0]-9.766542925609524*fl[0]; 
  incr[20] = (-9.57421875*fr[20])-9.57421875*fl[20]+6.393703176377302*fr[11]-6.393703176377302*fl[11]; 
  incr[21] = (-9.57421875*fr[21])-9.57421875*fl[21]+6.393703176377302*fr[12]-6.393703176377302*fl[12]; 
  incr[22] = (-11.54296875*fr[22])+11.54296875*fl[22]+12.36026325723227*fr[10]+12.36026325723227*fl[10]-8.254235307567585*fr[4]+8.254235307567585*fl[4]; 
  incr[23] = 3.69140625*fl[23]-3.69140625*fr[23]; 
  incr[24] = (-11.54296875*fr[24])+11.54296875*fl[24]+12.36026325723226*fr[13]+12.36026325723226*fl[13]-8.254235307567585*fr[7]+8.254235307567585*fl[7]; 
  incr[25] = (-11.54296875*fr[25])+11.54296875*fl[25]+12.36026325723226*fr[14]+12.36026325723226*fl[14]-8.254235307567585*fr[8]+8.254235307567585*fl[8]; 
  incr[26] = 3.69140625*fl[26]-3.69140625*fr[26]; 
  incr[27] = 3.69140625*fl[27]-3.69140625*fr[27]; 
  incr[28] = (-9.57421875*fr[28])-9.57421875*fl[28]+6.393703176377302*fr[17]-6.393703176377302*fl[17]; 
  incr[29] = (-9.57421875*fr[29])-9.57421875*fl[29]+6.393703176377302*fr[18]-6.393703176377302*fl[18]; 
  incr[30] = (-8.94140625*fr[30])-8.94140625*fl[30]+13.65782481176513*fr[15]-13.65782481176513*fl[15]-14.62486071398016*fr[5]-14.62486071398016*fl[5]+9.76654292560952*fr[1]-9.76654292560952*fl[1]; 
  incr[31] = (-8.94140625*fr[31])-8.94140625*fl[31]+13.65782481176513*fr[16]-13.65782481176513*fl[16]-14.62486071398016*fr[6]-14.62486071398016*fl[6]+9.76654292560952*fr[2]-9.76654292560952*fl[2]; 
  incr[32] = 3.69140625*fl[32]-3.69140625*fr[32]; 
  incr[33] = 3.69140625*fl[33]-3.69140625*fr[33]; 
  incr[34] = (-3.75*fr[34])+3.75*fl[34]+10.13860170372798*fr[19]+10.13860170372798*fl[19]-15.48651767229347*fr[9]+15.48651767229347*fl[9]+16.58303331777859*fr[3]+16.58303331777859*fl[3]-11.07421875*fr[0]+11.07421875*fl[0]; 

  outr[0] += incr[0]*rdxSq2nur; 
  outr[1] += incr[1]*rdxSq2nur; 
  outr[2] += incr[2]*rdxSq2nur; 
  outr[3] += incr[3]*rdxSq2nur; 
  outr[4] += incr[4]*rdxSq2nur; 
  outr[5] += incr[5]*rdxSq2nur; 
  outr[6] += incr[6]*rdxSq2nur; 
  outr[7] += incr[7]*rdxSq2nur; 
  outr[8] += incr[8]*rdxSq2nur; 
  outr[9] += incr[9]*rdxSq2nur; 
  outr[10] += incr[10]*rdxSq2nur; 
  outr[11] += incr[11]*rdxSq2nur; 
  outr[12] += incr[12]*rdxSq2nur; 
  outr[13] += incr[13]*rdxSq2nur; 
  outr[14] += incr[14]*rdxSq2nur; 
  outr[15] += incr[15]*rdxSq2nur; 
  outr[16] += incr[16]*rdxSq2nur; 
  outr[17] += incr[17]*rdxSq2nur; 
  outr[18] += incr[18]*rdxSq2nur; 
  outr[19] += incr[19]*rdxSq2nur; 
  outr[20] += incr[20]*rdxSq2nur; 
  outr[21] += incr[21]*rdxSq2nur; 
  outr[22] += incr[22]*rdxSq2nur; 
  outr[23] += incr[23]*rdxSq2nur; 
  outr[24] += incr[24]*rdxSq2nur; 
  outr[25] += incr[25]*rdxSq2nur; 
  outr[26] += incr[26]*rdxSq2nur; 
  outr[27] += incr[27]*rdxSq2nur; 
  outr[28] += incr[28]*rdxSq2nur; 
  outr[29] += incr[29]*rdxSq2nur; 
  outr[30] += incr[30]*rdxSq2nur; 
  outr[31] += incr[31]*rdxSq2nur; 
  outr[32] += incr[32]*rdxSq2nur; 
  outr[33] += incr[33]*rdxSq2nur; 
  outr[34] += incr[34]*rdxSq2nur; 

  outl[0] += -1.0*incr[0]*rdxSq2nul; 
  outl[1] += -1.0*incr[1]*rdxSq2nul; 
  outl[2] += -1.0*incr[2]*rdxSq2nul; 
  outl[3] += incr[3]*rdxSq2nul; 
  outl[4] += -1.0*incr[4]*rdxSq2nul; 
  outl[5] += incr[5]*rdxSq2nul; 
  outl[6] += incr[6]*rdxSq2nul; 
  outl[7] += -1.0*incr[7]*rdxSq2nul; 
  outl[8] += -1.0*incr[8]*rdxSq2nul; 
  outl[9] += -1.0*incr[9]*rdxSq2nul; 
  outl[10] += incr[10]*rdxSq2nul; 
  outl[11] += -1.0*incr[11]*rdxSq2nul; 
  outl[12] += -1.0*incr[12]*rdxSq2nul; 
  outl[13] += incr[13]*rdxSq2nul; 
  outl[14] += incr[14]*rdxSq2nul; 
  outl[15] += -1.0*incr[15]*rdxSq2nul; 
  outl[16] += -1.0*incr[16]*rdxSq2nul; 
  outl[17] += -1.0*incr[17]*rdxSq2nul; 
  outl[18] += -1.0*incr[18]*rdxSq2nul; 
  outl[19] += incr[19]*rdxSq2nul; 
  outl[20] += incr[20]*rdxSq2nul; 
  outl[21] += incr[21]*rdxSq2nul; 
  outl[22] += -1.0*incr[22]*rdxSq2nul; 
  outl[23] += -1.0*incr[23]*rdxSq2nul; 
  outl[24] += -1.0*incr[24]*rdxSq2nul; 
  outl[25] += -1.0*incr[25]*rdxSq2nul; 
  outl[26] += -1.0*incr[26]*rdxSq2nul; 
  outl[27] += -1.0*incr[27]*rdxSq2nul; 
  outl[28] += incr[28]*rdxSq2nul; 
  outl[29] += incr[29]*rdxSq2nul; 
  outl[30] += incr[30]*rdxSq2nul; 
  outl[31] += incr[31]*rdxSq2nul; 
  outl[32] += -1.0*incr[32]*rdxSq2nul; 
  outl[33] += -1.0*incr[33]*rdxSq2nul; 
  outl[34] += -1.0*incr[34]*rdxSq2nul; 

} 
