#include <MaxwellianOnBasisImpl.h>

class Indexer {
private:
  int nr, nc;

public:
  Indexer(int numRows, int numColumns) 
    : nr(numRows), nc(numColumns) {}

  int get(int row, int column) {
    return row*nc + column;
  }
  
};

void MaxwellianInnerLoop(/* Number density, drift speed, and thermal velocity squared
			    n[numConfOrds], u[numConfOrds, numVelDims], vth2[numConfOrds] */
			 double * n, double * u, double * vth2,
			 /* Pointer to output Maxwellian f[numPhaseBasis] */
			 double * fItr,
			 /* weights[numPhaseOrds], dz[numPhaseDims], zc[numPhaseDims] */
			 double * weights, double * dz, double * zc,
			 /* ordinates[numPhaseOrds, numPhaseDims] */
			 double * ordinates, 
			 /* basisAtOrdinates[numPhaseOrds, numPhaseBasis] */
			 double * basisAtOrdinates,
			 /* phaseToConfOrdMap[numPhaseOrds] */
			 double * phaseToConfOrdMap,
			 int numPhaseBasis,
			 int numConfOrds, int numPhaseOrds,
			 int numConfDims, int numPhaseDims) {

  int numVelDims = numPhaseDims - numConfDims;
  Indexer uIdx(numConfOrds, numVelDims);
  Indexer oIdx(numPhaseOrds, numPhaseDims);
  Indexer bIdx(numPhaseOrds, numPhaseBasis);

  for (int k = 0; k < numPhaseBasis; ++k)
    fItr[k] = 0;

  double maxwellNorm[numConfOrds];
  for (int confOrdIdx = 0; confOrdIdx < numConfOrds; ++confOrdIdx) {
    // Prepare the Maxvellian normalization
    double denom = 1.0;
    for (int d = 0; d < numVelDims; ++d)
      denom *= 2*M_PI*vth2[confOrdIdx];
    maxwellNorm[confOrdIdx] = n[confOrdIdx]/std::sqrt(denom);
  }
  
  for (int phaseOrdIdx = 0; phaseOrdIdx < numPhaseOrds; ++phaseOrdIdx) {
    int confOrdIdx = phaseToConfOrdMap[phaseOrdIdx]-1;
    double maxwellian = maxwellNorm[confOrdIdx];
    double v2 = 0;
    for (int d = 0; d < numVelDims; ++d) {
      double v = ordinates[oIdx.get(phaseOrdIdx, numConfDims+d)];
      // convert logical to physical coordinates
      v = 0.5*v*dz[numConfDims+d] + zc[numConfDims+d] - u[uIdx.get(confOrdIdx, d)];
      v2 += v*v;
    }
    maxwellian *= exp(-0.5*v2/vth2[confOrdIdx]);
    
    for (int k = 0; k < numPhaseBasis; ++k) 
      fItr[k] +=
	weights[phaseOrdIdx]*maxwellian*basisAtOrdinates[bIdx.get(phaseOrdIdx, k)];
  }
}
