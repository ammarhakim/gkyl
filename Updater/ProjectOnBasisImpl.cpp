#include <ProjectOnBasisImpl.h>
void projectF(double* f, double* weights, double* basisAtOrdinates, double* fv, int numVal, int numBasis, int numOrd)
{
  int offset = 0;
  for(int n=0; n<numVal; n++) {
    for(int k=0; k<numBasis; k++) {
      f[offset+k] = 0.0;
    }

    for(int imu=0; imu<numOrd; imu++) {
      double tmp = weights[imu]*fv[n+numVal*imu];
      for(int k=0; k<numBasis; k++) {
        f[offset+k] += tmp*basisAtOrdinates[k+numBasis*imu];
      } 
    }
    offset += numBasis;
  }
}
