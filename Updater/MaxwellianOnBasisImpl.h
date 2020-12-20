#include <cmath>
#include <stdio.h>
extern "C" {
  void MaxwellianInnerLoopOrg(double * n, double * u, double * vtSq,
			   double * fItr,
			   double * weights, double * dz, double * zc,
			   double * ordinates,
			   double * basisAtOrdinates,
			   double * phaseToConfOrdMap,
			   int numPhaseBasis,
			   int numConfOrds, int numPhaseOrds,
			   int numConfDims, int numPhaseDims);

  void GkMaxwellianInnerLoopOrg(double * n, double * u, double * vth2, double * bmag, double m_,
			     double * fItr,
			     double * weights, double * dz, double * zc,
			     double * ordinates,
			     double * basisAtOrdinates,
			     double * phaseToConfOrdMap,
			     int numPhaseBasis,
			     int numConfOrds, int numPhaseOrds,
			     int numConfDims, int numPhaseDims);
 
}
