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

void GkMaxwellianInnerLoop(/* Number density, drift speed, and thermal velocity squared
			    n[numConfOrds], uPar[numConfOrds], vtSq[numConfOrds],]*/
			   double * n, double * uPar, double * vtSq,
			   /* Magnetic field bmag[bumConfOrds and particle mass*/
			   double * bmag, double m_,
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
  Indexer oIdx(numPhaseOrds, numPhaseDims);
  Indexer bIdx(numPhaseOrds, numPhaseBasis);

  for (int k = 0; k < numPhaseBasis; ++k)
    fItr[k] = 0;

  double maxwellNorm[numConfOrds];
  for (int confOrdIdx = 0; confOrdIdx < numConfOrds; ++confOrdIdx) {
    // Prepare the Maxvellian normalization
    double denom = 1.0;
    if (numVelDims == 1)
      denom *= 2*M_PI*vtSq[confOrdIdx];
    else
      denom *= std::pow(2*M_PI*vtSq[confOrdIdx],3.0);
    maxwellNorm[confOrdIdx] = n[confOrdIdx]/std::sqrt(denom); 
  }
  
  for (int phaseOrdIdx = 0; phaseOrdIdx < numPhaseOrds; ++phaseOrdIdx) {
    int confOrdIdx = phaseToConfOrdMap[phaseOrdIdx]-1;
    double maxwellian = maxwellNorm[confOrdIdx];
    double v2 = 0;

    if (numVelDims == 1) {
      double v = ordinates[oIdx.get(phaseOrdIdx, numConfDims)];
      // convert logical to physical coordinates
      v = 0.5*v*dz[numConfDims] + zc[numConfDims] - uPar[confOrdIdx];
      v2 += v*v;

      if (vtSq[confOrdIdx] < 0)
      	//printf("GkMaxwellian: vtSq is less than 0! \n");
      	maxwellian *= 0; 
      else
	maxwellian *= exp(-0.5*v2/vtSq[confOrdIdx]);
    
      for (int k = 0; k < numPhaseBasis; ++k) 
	fItr[k] +=
	  weights[phaseOrdIdx]*(maxwellian+1e-40)*basisAtOrdinates[bIdx.get(phaseOrdIdx, k)];
    }

    else {
      double v = ordinates[oIdx.get(phaseOrdIdx, numConfDims)];
      double mu = ordinates[oIdx.get(phaseOrdIdx, numConfDims+1)];
      // convert logical to physical coordinates
      v = 0.5*v*dz[numConfDims] + zc[numConfDims] - uPar[confOrdIdx];
      v2 += v*v;
      mu = 0.5*mu*dz[numConfDims+1] + zc[numConfDims+1];
      // multiply by jacobian, the magnetic field
      maxwellian *= bmag[confOrdIdx]*exp((-0.5*v2-bmag[confOrdIdx]*mu/m_)/vtSq[confOrdIdx]);
    
      for (int k = 0; k < numPhaseBasis; ++k) 
	fItr[k] +=
	  weights[phaseOrdIdx]*(maxwellian+1e-40)*basisAtOrdinates[bIdx.get(phaseOrdIdx, k)];
    }

  }  
}
