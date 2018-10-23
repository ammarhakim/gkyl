#include <MaxwellianOnBasisImpl.h>

void MaxwellianInnerLoop(double * n, double * u, double * vth2,
			 double * fItr,
			 double * weights, double * dz, double * zc,
			 double * ordinates, double * basisAtOrdinates,
			 int numBasis, int numConfOrds, int numVelOrds,
			 int numConfDims, int numVelDims) {
  double denom, maxwellian, v, v2;
  int phaseIdx;
  int numPhaseDims = numConfDims + numVelDims;

  for (int k = 0; k < numBasis; ++k)
    fItr[k] = 0;

  // Loop over configuration space ordinates
  for (int confIdx = 0; confIdx < numConfOrds; ++confIdx) {
    // Prepare the Maxvellian normalization
    denom = 1.0;
    for (int d = 0; d < numVelDims; ++d)
      denom *= 2 * M_PI * vth2[confIdx];
    denom = 1 / sqrt(denom);

    // Loop over velocity space ordinates
    for (int velIdx = 0; velIdx < numVelOrds; ++velIdx) {
      maxwellian = n[confIdx] * denom;
      phaseIdx = confIdx*numVelOrds + velIdx;

      v2 = 0;
      for (int d = 0; d < numVelDims; ++d) {
	v = ordinates[phaseIdx*numPhaseDims + numConfDims + d];
	// convert logical to physical coordinates
	v = 0.5*dz[numConfDims + d]*v + zc[numConfDims + d];
	v2 += (v - u[confIdx*numVelDims + d]) * (v - u[confIdx*numVelDims + d]);
      }
      maxwellian *= exp(-0.5 * v2 / vth2[confIdx]);

      for (int k = 0; k < numBasis; ++k) 
	fItr[k] +=
	  weights[phaseIdx] * maxwellian * basisAtOrdinates[phaseIdx*numBasis + k];
    }
  }
}
