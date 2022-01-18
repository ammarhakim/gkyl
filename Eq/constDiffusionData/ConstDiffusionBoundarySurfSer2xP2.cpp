#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf2xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 

  if (edge < 0) {

  incr2[1] = 1.936491673103708*fr[4]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[3] = 1.936491673103709*fr[6]-1.5*fr[3]+0.8660254037844386*fr[2]; 
  incr2[4] = (-7.5*fr[4])+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 
  incr2[6] = (-7.5*fr[6])+5.809475019311126*fr[3]-3.354101966249684*fr[2]; 
  incr2[7] = 0.8660254037844387*fr[5]-1.5*fr[7]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 

  } else {

  incr2[1] = 1.936491673103708*fl[4]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[3] = 1.936491673103709*fl[6]+1.5*fl[3]+0.8660254037844386*fl[2]; 
  incr2[4] = (-7.5*fl[4])-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 
  incr2[6] = (-7.5*fl[6])-5.809475019311126*fl[3]-3.354101966249684*fl[2]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844387*fl[5]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf2xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 

  if (edge < 0) {

  incr2[2] = 1.936491673103708*fr[5]-1.5*fr[2]+0.8660254037844386*fr[0]; 
  incr2[3] = 1.936491673103709*fr[7]-1.5*fr[3]+0.8660254037844386*fr[1]; 
  incr2[5] = (-7.5*fr[5])+5.809475019311125*fr[2]-3.354101966249685*fr[0]; 
  incr2[6] = 0.8660254037844387*fr[4]-1.5*fr[6]; 
  incr2[7] = (-7.5*fr[7])+5.809475019311126*fr[3]-3.354101966249684*fr[1]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 

  } else {

  incr2[2] = 1.936491673103708*fl[5]+1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[3] = 1.936491673103709*fl[7]+1.5*fl[3]+0.8660254037844386*fl[1]; 
  incr2[5] = (-7.5*fl[5])-5.809475019311125*fl[2]-3.354101966249685*fl[0]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844387*fl[4]; 
  incr2[7] = (-7.5*fl[7])-5.809475019311126*fl[3]-3.354101966249684*fl[1]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 
  double incr3[8]; 
  double incr4[8]; 

  if (edge < 0) {


  incr2[1] = 2.8125*fr[1]-5.083290641897234*fr[4]; 
  incr2[3] = 2.8125*fr[3]-5.083290641897235*fr[6]; 
  incr2[4] = 19.6875*fr[4]-10.89276566120836*fr[1]; 
  incr2[6] = 19.6875*fr[6]-10.89276566120836*fr[3]; 
  incr2[7] = 2.8125*fr[7]; 



  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 

  } else {


  incr2[1] = 5.083290641897234*fl[4]+2.8125*fl[1]; 
  incr2[3] = 5.083290641897235*fl[6]+2.8125*fl[3]; 
  incr2[4] = 19.6875*fl[4]+10.89276566120836*fl[1]; 
  incr2[6] = 19.6875*fl[6]+10.89276566120836*fl[3]; 
  incr2[7] = 2.8125*fl[7]; 



  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 
  double incr3[8]; 
  double incr4[8]; 

  if (edge < 0) {


  incr2[2] = 2.8125*fr[2]-5.083290641897234*fr[5]; 
  incr2[3] = 2.8125*fr[3]-5.083290641897235*fr[7]; 
  incr2[5] = 19.6875*fr[5]-10.89276566120836*fr[2]; 
  incr2[6] = 2.8125*fr[6]; 
  incr2[7] = 19.6875*fr[7]-10.89276566120836*fr[3]; 



  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 

  } else {


  incr2[2] = 5.083290641897234*fl[5]+2.8125*fl[2]; 
  incr2[3] = 5.083290641897235*fl[7]+2.8125*fl[3]; 
  incr2[5] = 19.6875*fl[5]+10.89276566120836*fl[2]; 
  incr2[6] = 2.8125*fl[6]; 
  incr2[7] = 19.6875*fl[7]+10.89276566120836*fl[3]; 



  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf2xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 
  double incr3[8]; 
  double incr4[8]; 
  double incr5[8]; 
  double incr6[8]; 

  if (edge < 0) {


  incr2[1] = 19.06233990711463*fr[4]-4.921875*fr[1]; 
  incr2[3] = 19.06233990711463*fr[6]-4.921875*fr[3]; 
  incr2[4] = 19.06233990711463*fr[1]-73.828125*fr[4]; 
  incr2[6] = 19.06233990711463*fr[3]-73.828125*fr[6]; 
  incr2[7] = -4.921875*fr[7]; 





  outr[1] += incr2[1]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 

  } else {


  incr2[1] = (-19.06233990711463*fl[4])-4.921875*fl[1]; 
  incr2[3] = (-19.06233990711463*fl[6])-4.921875*fl[3]; 
  incr2[4] = (-73.828125*fl[4])-19.06233990711463*fl[1]; 
  incr2[6] = (-73.828125*fl[6])-19.06233990711463*fl[3]; 
  incr2[7] = -4.921875*fl[7]; 





  outl[1] += incr2[1]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf2xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 
  double incr3[8]; 
  double incr4[8]; 
  double incr5[8]; 
  double incr6[8]; 

  if (edge < 0) {


  incr2[2] = 19.06233990711463*fr[5]-4.921875*fr[2]; 
  incr2[3] = 19.06233990711463*fr[7]-4.921875*fr[3]; 
  incr2[5] = 19.06233990711463*fr[2]-73.828125*fr[5]; 
  incr2[6] = -4.921875*fr[6]; 
  incr2[7] = 19.06233990711463*fr[3]-73.828125*fr[7]; 





  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 

  } else {


  incr2[2] = (-19.06233990711463*fl[5])-4.921875*fl[2]; 
  incr2[3] = (-19.06233990711463*fl[7])-4.921875*fl[3]; 
  incr2[5] = (-73.828125*fl[5])-19.06233990711463*fl[2]; 
  incr2[6] = -4.921875*fl[6]; 
  incr2[7] = (-73.828125*fl[7])-19.06233990711463*fl[3]; 





  outl[2] += incr2[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf2xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[16]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 

  if (edge < 0) {

  incr2[1] = 2.165063509461096*fr[4]*nul[4]-1.677050983124842*fr[1]*nul[4]+0.9682458365518543*fr[0]*nul[4]+1.677050983124842*nul[1]*fr[4]+0.9682458365518541*nul[0]*fr[4]-1.299038105676658*fr[1]*nul[1]+0.75*fr[0]*nul[1]-0.75*nul[0]*fr[1]+0.4330127018922193*fr[0]*nul[0]; 
  incr2[3] = 2.165063509461097*nul[4]*fr[6]+1.677050983124842*nul[1]*fr[6]+0.9682458365518543*nul[0]*fr[6]-1.677050983124842*fr[3]*nul[4]+0.9682458365518543*fr[2]*nul[4]-1.299038105676658*nul[1]*fr[3]-0.75*nul[0]*fr[3]+0.75*nul[1]*fr[2]+0.4330127018922193*nul[0]*fr[2]; 
  incr2[4] = (-8.385254915624213*fr[4]*nul[4])+6.495190528383289*fr[1]*nul[4]-3.75*fr[0]*nul[4]-6.495190528383286*nul[1]*fr[4]-3.75*nul[0]*fr[4]+5.031152949374527*fr[1]*nul[1]-2.904737509655563*fr[0]*nul[1]+2.904737509655563*nul[0]*fr[1]-1.677050983124842*fr[0]*nul[0]; 
  incr2[6] = (-8.385254915624213*nul[4]*fr[6])-6.495190528383286*nul[1]*fr[6]-3.75*nul[0]*fr[6]+6.49519052838329*fr[3]*nul[4]-3.75*fr[2]*nul[4]+5.031152949374526*nul[1]*fr[3]+2.904737509655563*nul[0]*fr[3]-2.904737509655563*nul[1]*fr[2]-1.677050983124842*nul[0]*fr[2]; 
  incr2[7] = (-1.677050983124842*nul[4]*fr[7])-1.299038105676658*nul[1]*fr[7]-0.75*nul[0]*fr[7]+0.9682458365518543*nul[4]*fr[5]+0.75*nul[1]*fr[5]+0.4330127018922194*nul[0]*fr[5]; 

  outr[1] += incr2[1]*rdxFr; 
  outr[3] += incr2[3]*rdxFr; 
  outr[4] += incr2[4]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 

  } else {

  incr2[1] = 2.165063509461096*fl[4]*nul[4]+1.677050983124842*fl[1]*nul[4]+0.9682458365518543*fl[0]*nul[4]+1.677050983124842*nul[1]*fl[4]+0.9682458365518541*nul[0]*fl[4]+1.299038105676658*fl[1]*nul[1]+0.75*fl[0]*nul[1]+0.75*nul[0]*fl[1]+0.4330127018922193*fl[0]*nul[0]; 
  incr2[3] = 2.165063509461097*nul[4]*fl[6]+1.677050983124842*nul[1]*fl[6]+0.9682458365518543*nul[0]*fl[6]+1.677050983124842*fl[3]*nul[4]+0.9682458365518543*fl[2]*nul[4]+1.299038105676658*nul[1]*fl[3]+0.75*nul[0]*fl[3]+0.75*nul[1]*fl[2]+0.4330127018922193*nul[0]*fl[2]; 
  incr2[4] = (-8.385254915624213*fl[4]*nul[4])-6.495190528383289*fl[1]*nul[4]-3.75*fl[0]*nul[4]-6.495190528383286*nul[1]*fl[4]-3.75*nul[0]*fl[4]-5.031152949374527*fl[1]*nul[1]-2.904737509655563*fl[0]*nul[1]-2.904737509655563*nul[0]*fl[1]-1.677050983124842*fl[0]*nul[0]; 
  incr2[6] = (-8.385254915624213*nul[4]*fl[6])-6.495190528383286*nul[1]*fl[6]-3.75*nul[0]*fl[6]-6.49519052838329*fl[3]*nul[4]-3.75*fl[2]*nul[4]-5.031152949374526*nul[1]*fl[3]-2.904737509655563*nul[0]*fl[3]-2.904737509655563*nul[1]*fl[2]-1.677050983124842*nul[0]*fl[2]; 
  incr2[7] = 1.677050983124842*nul[4]*fl[7]+1.299038105676658*nul[1]*fl[7]+0.75*nul[0]*fl[7]+0.9682458365518543*nul[4]*fl[5]+0.75*nul[1]*fl[5]+0.4330127018922194*nul[0]*fl[5]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 
  outl[3] += -1.0*incr2[3]*rdxFl; 
  outl[4] += incr2[4]*rdxFl; 
  outl[6] += incr2[6]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf2xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[16]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[1]*dxl[1]); 
  double rdxFr = 4.0/(dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 

  if (edge < 0) {

  incr2[2] = 2.165063509461096*fr[5]*nul[13]-1.677050983124842*fr[2]*nul[13]+0.9682458365518543*fr[0]*nul[13]+1.677050983124842*fr[5]*nul[10]-1.299038105676658*fr[2]*nul[10]+0.75*fr[0]*nul[10]+0.9682458365518541*fr[5]*nul[8]-0.75*fr[2]*nul[8]+0.4330127018922193*fr[0]*nul[8]; 
  incr2[3] = 2.165063509461097*fr[7]*nul[13]-1.677050983124842*fr[3]*nul[13]+0.9682458365518543*fr[1]*nul[13]+1.677050983124842*fr[7]*nul[10]-1.299038105676658*fr[3]*nul[10]+0.75*fr[1]*nul[10]+0.9682458365518543*fr[7]*nul[8]-0.75*fr[3]*nul[8]+0.4330127018922193*fr[1]*nul[8]; 
  incr2[5] = (-8.385254915624213*fr[5]*nul[13])+6.495190528383289*fr[2]*nul[13]-3.75*fr[0]*nul[13]-6.495190528383286*fr[5]*nul[10]+5.031152949374527*fr[2]*nul[10]-2.904737509655563*fr[0]*nul[10]-3.75*fr[5]*nul[8]+2.904737509655563*fr[2]*nul[8]-1.677050983124842*fr[0]*nul[8]; 
  incr2[6] = (-1.677050983124842*fr[6]*nul[13])+0.9682458365518543*fr[4]*nul[13]-1.299038105676658*fr[6]*nul[10]+0.75*fr[4]*nul[10]-0.75*fr[6]*nul[8]+0.4330127018922194*fr[4]*nul[8]; 
  incr2[7] = (-8.385254915624213*fr[7]*nul[13])+6.49519052838329*fr[3]*nul[13]-3.75*fr[1]*nul[13]-6.495190528383286*fr[7]*nul[10]+5.031152949374526*fr[3]*nul[10]-2.904737509655563*fr[1]*nul[10]-3.75*fr[7]*nul[8]+2.904737509655563*fr[3]*nul[8]-1.677050983124842*fr[1]*nul[8]; 

  outr[2] += incr2[2]*rdxFr; 
  outr[3] += incr2[3]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 

  } else {

  incr2[2] = 2.165063509461096*fl[5]*nul[13]+1.677050983124842*fl[2]*nul[13]+0.9682458365518543*fl[0]*nul[13]+1.677050983124842*fl[5]*nul[10]+1.299038105676658*fl[2]*nul[10]+0.75*fl[0]*nul[10]+0.9682458365518541*fl[5]*nul[8]+0.75*fl[2]*nul[8]+0.4330127018922193*fl[0]*nul[8]; 
  incr2[3] = 2.165063509461097*fl[7]*nul[13]+1.677050983124842*fl[3]*nul[13]+0.9682458365518543*fl[1]*nul[13]+1.677050983124842*fl[7]*nul[10]+1.299038105676658*fl[3]*nul[10]+0.75*fl[1]*nul[10]+0.9682458365518543*fl[7]*nul[8]+0.75*fl[3]*nul[8]+0.4330127018922193*fl[1]*nul[8]; 
  incr2[5] = (-8.385254915624213*fl[5]*nul[13])-6.495190528383289*fl[2]*nul[13]-3.75*fl[0]*nul[13]-6.495190528383286*fl[5]*nul[10]-5.031152949374527*fl[2]*nul[10]-2.904737509655563*fl[0]*nul[10]-3.75*fl[5]*nul[8]-2.904737509655563*fl[2]*nul[8]-1.677050983124842*fl[0]*nul[8]; 
  incr2[6] = 1.677050983124842*fl[6]*nul[13]+0.9682458365518543*fl[4]*nul[13]+1.299038105676658*fl[6]*nul[10]+0.75*fl[4]*nul[10]+0.75*fl[6]*nul[8]+0.4330127018922194*fl[4]*nul[8]; 
  incr2[7] = (-8.385254915624213*fl[7]*nul[13])-6.49519052838329*fl[3]*nul[13]-3.75*fl[1]*nul[13]-6.495190528383286*fl[7]*nul[10]-5.031152949374526*fl[3]*nul[10]-2.904737509655563*fl[1]*nul[10]-3.75*fl[7]*nul[8]-2.904737509655563*fl[3]*nul[8]-1.677050983124842*fl[1]*nul[8]; 

  outl[2] += -1.0*incr2[2]*rdxFl; 
  outl[3] += -1.0*incr2[3]*rdxFl; 
  outl[5] += incr2[5]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[7] += incr2[7]*rdxFl; 

  }

} 
