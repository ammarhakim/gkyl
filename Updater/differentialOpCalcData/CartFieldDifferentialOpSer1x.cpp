#include <CartFieldDifferentialOpModDecl.h> 
void CartFieldDifferentialOpDxxVol1xSer_P1(const double *w, const double *dx, const int *idx, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates. 
  // dx[1]:  Cell spacing. 
  // idx[1]: Current grid index.
  // f:      Input function. 
  // out:    Incremented output 
  double rdx2[1]; 
  double rdxSq4[1]; 
  rdx2[0] = 2.0/dx[0]; 
  rdxSq4[0] = 4.0/(dx[0]*dx[0]); 


} 
