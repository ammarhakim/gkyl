#include <CartFieldDifferentialOpModDecl.h> 
void CartFieldDifferentialOpDxxRecoverySurf1xMax_X_P1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:       Cell-center coordinates. 
  // dx[1]:      Cell spacing. 
  // idx[1]:     current grid index.
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdx = 1.0/dxl[0]; 
  double rdx2L = 2.0/dxl[0]; 
  double rdx2R = 2.0/dxr[0]; 
  double rdxSq4L = 4.0/(dxl[0]*dxl[0]); 
  double rdxSq4R = 4.0/(dxr[0]*dxr[0]); 

  double Gdiff[2]; 
  Gdiff[0] = (-2.165063509461096*(fr[1]+fl[1]))+2.25*fr[0]-2.25*fl[0]; 

  double Ghat[2]; 
  for(unsigned short int i=0; i<2; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = Gdiff[0]*rdx; 

  double incr1[2]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = 0.8660254037844386*Ghat[0]; 

  double incr2[2]; 
  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*(fr[0]+fl[0]); 

  outr[0] += incr1[0]*rdx2R; 
  outr[1] += incr2[1]*rdxSq4R+incr1[1]*rdx2R; 

  outl[0] += -1.0*incr1[0]*rdx2L; 
  outl[1] += incr1[1]*rdx2L-1.0*incr2[1]*rdxSq4L; 

} 
