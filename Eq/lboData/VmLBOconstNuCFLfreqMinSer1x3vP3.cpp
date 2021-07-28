#include <VmLBOModDecl.h> 
void VmLBOconstNuCFLfreqMin1x3vSerP3(const double *w, const double *dxv, const double *Lv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, double *cflFreqDrag, double *cflFreqDiff) 
{ 
  // w[4]:        cell-center coordinates. 
  // dxv[4]:      cell spacing. 
  // nuSum:       collisionalities added (self and cross species collisionalities). 
  // nuUSum:      sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum:   sum of thermal speeds squared time their respective collisionalities. 
  // Lv[3]:       domain length in velocity space.
  // cflFreqDrag: min CFL frequency supported by drag terms.
  // cflFreqDiff: min CFL frequency supported by diffusion terms.

  cflFreqDrag[0] = std::fmax(cflFreqDrag[0], std::fabs(((-3.141592653589793*w[1]*nuSum)-2.483647066449025*nuUSum[2]+2.221441469079183*nuUSum[0])/Lv[0])+std::fabs(((-3.141592653589793*w[2]*nuSum)-2.483647066449025*nuUSum[6]+2.221441469079183*nuUSum[4])/Lv[1])+std::fabs(((-3.141592653589793*w[3]*nuSum)-2.483647066449025*nuUSum[10]+2.221441469079183*nuUSum[8])/Lv[2])); 
  cflFreqDiff[0] = std::fmax(cflFreqDiff[0], std::fabs((63.80675839669836*nuVtSqSum[0]-71.33812459946151*nuVtSqSum[2])/pow(Lv[0],2.0))+std::fabs((63.80675839669836*nuVtSqSum[0]-71.33812459946151*nuVtSqSum[2])/pow(Lv[1],2.0))+std::fabs((63.80675839669836*nuVtSqSum[0]-71.33812459946151*nuVtSqSum[2])/pow(Lv[2],2.0))); 

} 
