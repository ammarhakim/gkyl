#include <VmLBOModDecl.h> 
void VmLBOconstNuCFLfreqMin2x3vSerP1(const double *w, const double *dxv, const double *Lv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, double *cflFreqDrag, double *cflFreqDiff) 
{ 
  // w[5]:        cell-center coordinates. 
  // dxv[5]:      cell spacing. 
  // nuSum:       collisionalities added (self and cross species collisionalities). 
  // nuUSum:      sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum:   sum of thermal speeds squared time their respective collisionalities. 
  // Lv[3]:       domain length in velocity space.
  // cflFreqDrag: min CFL frequency supported by drag terms.
  // cflFreqDiff: min CFL frequency supported by diffusion terms.

  cflFreqDrag[0] = std::fmax(cflFreqDrag[0], std::fabs((1.570796326794897*nuUSum[0]-3.141592653589793*w[2]*nuSum)/Lv[0])+std::fabs((1.570796326794897*nuUSum[4]-3.141592653589793*w[3]*nuSum)/Lv[1])+std::fabs((1.570796326794897*nuUSum[8]-3.141592653589793*w[4]*nuSum)/Lv[2])); 
  cflFreqDiff[0] = std::fmax(cflFreqDiff[0], std::fabs((26.31894506957162*nuVtSqSum[0])/pow(Lv[0],2.0))+std::fabs((26.31894506957162*nuVtSqSum[0])/pow(Lv[1],2.0))+std::fabs((26.31894506957162*nuVtSqSum[0])/pow(Lv[2],2.0))); 

} 
