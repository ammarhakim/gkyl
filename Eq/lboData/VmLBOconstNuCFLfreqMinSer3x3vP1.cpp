#include <VmLBOModDecl.h> 
void VmLBOconstNuCFLfreqMin3x3vSerP1(const double *w, const double *dxv, const double *Lv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, double *cflFreqDrag, double *cflFreqDiff) 
{ 
  // w[6]:        cell-center coordinates. 
  // dxv[6]:      cell spacing. 
  // nuSum:       collisionalities added (self and cross species collisionalities). 
  // nuUSum:      sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum:   sum of thermal speeds squared time their respective collisionalities. 
  // Lv[3]:       domain length in velocity space.
  // cflFreqDrag: min CFL frequency supported by drag terms.
  // cflFreqDiff: min CFL frequency supported by diffusion terms.

  cflFreqDrag[0] = std::fmax(cflFreqDrag[0], std::fabs((1.110720734539592*nuUSum[0]-3.141592653589793*w[3]*nuSum)/Lv[0])+std::fabs((1.110720734539592*nuUSum[8]-3.141592653589793*w[4]*nuSum)/Lv[1])+std::fabs((1.110720734539592*nuUSum[16]-3.141592653589793*w[5]*nuSum)/Lv[2])); 
  cflFreqDiff[0] = std::fmax(cflFreqDiff[0], std::fabs((18.61030453237035*nuVtSqSum[0])/pow(Lv[0],2.0))+std::fabs((18.61030453237035*nuVtSqSum[0])/pow(Lv[1],2.0))+std::fabs((18.61030453237035*nuVtSqSum[0])/pow(Lv[2],2.0))); 

} 
