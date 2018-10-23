#include <math.h>
extern "C" {
  void MaxwellianInnerLoop(double * n, double * u, double * vth2,
			   double * fItr,
			   double * weights, double * dz, double * zc,
			   double * ordinates, double * basisAtOrdinates,
			   int numBasis, int numConfOrds, int numVelOrds,
			   int numConfDims, int numVelDims);
}
