#include <VmLBOModDecl.h> 
void VmLBOconstNuCFLfreqMin1x1vTensorP2(const double *w, const double *dxv, const double *Lv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, double *cflFreqDrag, double *cflFreqDiff) 
{ 
  // w[2]:        cell-center coordinates. 
  // dxv[2]:      cell spacing. 
  // nuSum:       collisionalities added (self and cross species collisionalities). 
  // nuUSum:      sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum:   sum of thermal speeds squared time their respective collisionalities. 
  // Lv[1]:       domain length in velocity space.
  // cflFreqDrag: min CFL frequency supported by drag terms.
  // cflFreqDiff: min CFL frequency supported by diffusion terms.

  cflFreqDrag[0] = std::fmax(cflFreqDrag[0], std::fabs(((-3.141592653589793*w[1]*nuSum)-2.483647066449025*nuUSum[2]+2.221441469079183*nuUSum[0])/Lv[0])); 
  cflFreqDiff[0] = std::fmax(cflFreqDiff[0], std::fabs((50.24782223739994*nuVtSqSum[0]-56.17877312207591*nuVtSqSum[2])/pow(Lv[0],2.0))); 

} 
