#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf1xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 

  incr1[0] = 0.6821077598838398*fr[3]+0.6821077598838398*fl[3]-1.554766015605323*fr[2]+1.554766015605323*fl[2]+1.935025511580855*fr[1]+1.935025511580855*fl[1]-1.3671875*fr[0]+1.3671875*fl[0]; 
  incr1[1] = (-1.181445296355802*fr[3])-1.181445296355802*fl[3]+2.692933732909844*fr[2]-2.692933732909844*fl[2]-3.3515625*fr[1]-3.3515625*fl[1]+2.368038213473074*fr[0]-2.368038213473074*fl[0]; 
  incr1[2] = 1.52523931908037*fr[3]+1.52523931908037*fl[3]-3.4765625*fr[2]+3.4765625*fl[2]+4.326848582091099*fr[1]+4.326848582091099*fl[1]-3.057124187987994*fr[0]+3.057124187987994*fl[0]; 
  incr1[3] = (-1.8046875*fr[3])-1.8046875*fl[3]+4.113524224186452*fr[2]-4.113524224186452*fl[2]-5.119596284208478*fr[1]-5.119596284208478*fl[1]+3.61723812059612*fr[0]-3.61723812059612*fl[0]; 

  incr2[1] = (-0.3273268353539885*fr[3])+0.3273268353539885*fl[3]+0.605153647844909*fr[2]+0.605153647844909*fl[2]-0.65625*fr[1]+0.65625*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[2] = 1.267731382092775*fr[3]-1.267731382092775*fl[3]-2.34375*fr[2]-2.34375*fl[2]+2.541645320948617*fr[1]-2.541645320948617*fl[1]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 
  incr2[3] = (-3.0*fr[3])+3.0*fl[3]+5.546324796655891*fr[2]+5.546324796655891*fl[2]-6.01463059962954*fr[1]+6.01463059962954*fl[1]+3.968626966596886*fr[0]+3.968626966596886*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul-1.0*incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf1xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 
  double incr3[4]; 
  double incr4[4]; 

  incr1[0] = (-15.34742459738639*fr[3])-15.34742459738639*fl[3]+26.59698043549555*fr[2]-26.59698043549555*fl[2]-24.86440124146728*fr[1]-24.86440124146728*fl[1]+14.35546875*fr[0]-14.35546875*fl[0]; 
  incr1[1] = 26.58251916800555*fr[3]+26.58251916800555*fl[3]-46.06732144219369*fr[2]+46.06732144219369*fl[2]+43.06640625*fr[1]+43.06640625*fl[1]-24.86440124146728*fr[0]+24.86440124146728*fl[0]; 
  incr1[2] = (-34.31788467930833*fr[3])-34.31788467930833*fl[3]+59.47265625*fr[2]-59.47265625*fl[2]-55.598491395751*fr[1]-55.598491395751*fl[1]+32.09980397387394*fr[0]-32.09980397387394*fl[0]; 
  incr1[3] = 40.60546875*fr[3]+40.60546875*fl[3]-70.36899585757163*fr[2]+70.36899585757163*fl[2]+65.78502218344809*fr[1]+65.78502218344809*fl[1]-37.98100026625927*fr[0]+37.98100026625927*fl[0]; 

  incr2[1] = 4.909902530309828*fr[3]-4.909902530309828*fl[3]-5.264836736250706*fr[2]-5.264836736250706*fl[2]+2.109375*fr[1]-2.109375*fl[1]; 
  incr2[2] = (-19.01597073139163*fr[3])+19.01597073139163*fl[3]+20.390625*fr[2]+20.390625*fl[2]-8.169574245906269*fr[1]+8.169574245906269*fl[1]; 
  incr2[3] = 45.0*fr[3]-45.0*fl[3]-48.25302573090624*fr[2]-48.25302573090624*fl[2]+19.33274121309494*fr[1]-19.33274121309494*fl[1]; 

  incr3[2] = 4.57571795724111*fr[3]+4.57571795724111*fl[3]-10.4296875*fr[2]+10.4296875*fl[2]+12.98054574627329*fr[1]+12.98054574627329*fl[1]-9.171372563963983*fr[0]+9.171372563963983*fl[0]; 
  incr3[3] = (-27.0703125*fr[3])-27.0703125*fl[3]+61.70286336279679*fr[2]-61.70286336279679*fl[2]-76.79394426312716*fr[1]-76.79394426312716*fl[1]+54.2585718089418*fr[0]-54.2585718089418*fl[0]; 

  incr4[3] = (-7.5*fr[3])+7.5*fl[3]+13.86581199163973*fr[2]+13.86581199163973*fl[2]-15.03657649907385*fr[1]+15.03657649907385*fl[1]+9.921567416492215*fr[0]+9.921567416492215*fl[0]; 

  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += (-1.0*incr2[1]*rdxFnur)-1.0*incr1[1]*rdxFnur; 
  outr[2] += (-1.0*incr3[2]*rdxFnur)-1.0*incr2[2]*rdxFnur-1.0*incr1[2]*rdxFnur; 
  outr[3] += (-1.0*incr4[3]*rdxFnur)-1.0*incr3[3]*rdxFnur-1.0*incr2[3]*rdxFnur-1.0*incr1[3]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr2[1]*rdxFnul-1.0*incr1[1]*rdxFnul; 
  outl[2] += incr3[2]*rdxFnul-1.0*incr2[2]*rdxFnul+incr1[2]*rdxFnul; 
  outl[3] += incr4[3]*rdxFnul-1.0*incr3[3]*rdxFnul+incr2[3]*rdxFnul-1.0*incr1[3]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf1xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
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

  incr1[0] = 214.8639443634096*fr[3]+214.8639443634096*fl[3]-313.6609416875682*fr[2]+313.6609416875682*fl[2]+268.5355334078465*fr[1]+268.5355334078465*fl[1]-155.0390625*fr[0]+155.0390625*fl[0]; 
  incr1[1] = (-372.1552683520777*fr[3])-372.1552683520777*fl[3]+543.2766873527669*fr[2]-543.2766873527669*fl[2]-465.1171875*fr[1]-465.1171875*fl[1]+268.5355334078465*fr[0]-268.5355334078465*fl[0]; 
  incr1[2] = 480.4503855103166*fr[3]+480.4503855103166*fl[3]-701.3671875*fr[2]+701.3671875*fl[2]+600.4637070741106*fr[1]+600.4637070741106*fl[1]-346.6778829178386*fr[0]+346.6778829178386*fl[0]; 
  incr1[3] = (-568.4765625*fr[3])-568.4765625*fl[3]+829.8688476996379*fr[2]-829.8688476996379*fl[2]-710.4782395812391*fr[1]-710.4782395812391*fl[1]+410.1948028756001*fr[0]-410.1948028756001*fl[0]; 

  incr2[1] = (-51.55397656825318*fr[3])+51.55397656825318*fl[3]+38.12467981422925*fr[2]+38.12467981422925*fl[2]-9.84375*fr[1]+9.84375*fl[1]; 
  incr2[2] = 199.6676926796121*fr[3]-199.6676926796121*fl[3]-147.65625*fr[2]-147.65625*fl[2]+38.12467981422925*fr[1]-38.12467981422925*fl[1]; 
  incr2[3] = (-472.5*fr[3])+472.5*fl[3]+349.4184621893212*fr[2]+349.4184621893212*fl[2]-90.21945899444307*fr[1]+90.21945899444307*fl[1]; 

  incr3[2] = (-102.953654037925*fr[3])-102.953654037925*fl[3]+178.41796875*fr[2]-178.41796875*fl[2]-166.795474187253*fr[1]-166.795474187253*fl[1]+96.29941192162183*fr[0]-96.29941192162183*fl[0]; 
  incr3[3] = 609.08203125*fr[3]+609.08203125*fl[3]-1055.534937863574*fr[2]+1055.534937863574*fl[2]+986.7753327517214*fr[1]+986.7753327517214*fl[1]-569.715003993889*fr[0]+569.715003993889*fl[0]; 

  incr4[3] = 112.5*fr[3]-112.5*fl[3]-120.6325643272656*fr[2]-120.6325643272656*fl[2]+48.33185303273736*fr[1]-48.33185303273736*fl[1]; 



  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr3[2]*rdxFnur+incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
  outr[3] += incr4[3]*rdxFnur+incr3[3]*rdxFnur+incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += (-1.0*incr3[2]*rdxFnul)+incr2[2]*rdxFnul-1.0*incr1[2]*rdxFnul; 
  outl[3] += (-1.0*incr4[3]*rdxFnul)+incr3[3]*rdxFnul-1.0*incr2[3]*rdxFnul+incr1[3]*rdxFnul; 

} 
