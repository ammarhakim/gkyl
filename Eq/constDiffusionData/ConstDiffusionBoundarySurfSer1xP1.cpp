#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf1xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[1]:    current grid index.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[2]; 
  double incr2[2]; 

  if (idxr[0] == 1) {

  incr2[1] = 0.8660254037844386*fr[0]-1.5*fr[1]; 

  outr[1] += incr2[1]*rdxFnur; 

  } else {

  incr2[1] = 1.5*fl[1]+0.8660254037844386*fl[0]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf1xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[1]:    current grid index.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[2]; 
  double incr2[2]; 
  double incr3[2]; 
  double incr4[2]; 

  if (idxr[0] == 1) {


  incr2[1] = 1.5*fr[1]; 



  outr[1] += -1.0*incr2[1]*rdxFnur; 

  } else {


  incr2[1] = 1.5*fl[1]; 



  outl[1] += -1.0*incr2[1]*rdxFnul; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf1xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[1]:    current grid index.
  // nu[2]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[2]; 
  double incr2[2]; 

  if (idxr[0] == 1) {

  incr2[1] = (-1.837117307087383*fr[1]*nul[1])+1.060660171779821*fr[0]*nul[1]-1.060660171779821*nul[0]*fr[1]+0.6123724356957944*fr[0]*nul[0]; 

  outr[1] += incr2[1]*rdxFr; 

  } else {

  incr2[1] = 1.837117307087383*fl[1]*nul[1]+1.060660171779821*fl[0]*nul[1]+1.060660171779821*nul[0]*fl[1]+0.6123724356957944*fl[0]*nul[0]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 

  }

} 
