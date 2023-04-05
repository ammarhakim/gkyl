#ifndef BIMAXWELLIANONBASIS_MOD_DECL_H 
#define BIMAXWELLIANONBASIS_MOD_DECL_H 
#include <cmath> 

extern "C" { 

void GkBiMaxwellianOnBasisGauss1x1vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtparSq, const double *vtperpSq, const double *bmag, double *flowUOrd, double *vtparSqOrd, double *vtperpSqOrd, double *fMFacOrd, double *bmagOrd);
void GkBiMaxwellianOnBasisGauss1x1vSer_P1_phaseQuad(const double *flowUOrd, const double *vtparSqOrd, const double *vtperpSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkBiMaxwellianOnBasisGauss1x1vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtparSq, const double *vtperpSq, const double *bmag, double *flowUOrd, double *vtparSqOrd, double *vtperpSqOrd, double *fMFacOrd, double *bmagOrd);
void GkBiMaxwellianOnBasisGauss1x1vSer_P2_phaseQuad(const double *flowUOrd, const double *vtparSqOrd, const double *vtperpSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkBiMaxwellianOnBasisGauss1x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtparSq, const double *vtperpSq, const double *bmag, double *flowUOrd, double *vtparSqOrd, double *vtperpSqOrd, double *fMFacOrd, double *bmagOrd);
void GkBiMaxwellianOnBasisGauss1x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtparSqOrd, const double *vtperpSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkBiMaxwellianOnBasisGauss1x2vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtparSq, const double *vtperpSq, const double *bmag, double *flowUOrd, double *vtparSqOrd, double *vtperpSqOrd, double *fMFacOrd, double *bmagOrd);
void GkBiMaxwellianOnBasisGauss1x2vSer_P2_phaseQuad(const double *flowUOrd, const double *vtparSqOrd, const double *vtperpSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkBiMaxwellianOnBasisGauss1x1vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtparSq, const double *vtperpSq, const double *bmag, double *flowUOrd, double *vtparSqOrd, double *vtperpSqOrd, double *fMFacOrd, double *bmagOrd);
void GkBiMaxwellianOnBasisGauss1x1vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtparSqOrd, const double *vtperpSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkBiMaxwellianOnBasisGauss1x1vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtparSq, const double *vtperpSq, const double *bmag, double *flowUOrd, double *vtparSqOrd, double *vtperpSqOrd, double *fMFacOrd, double *bmagOrd);
void GkBiMaxwellianOnBasisGauss1x1vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtparSqOrd, const double *vtperpSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkBiMaxwellianOnBasisGauss1x2vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtparSq, const double *vtperpSq, const double *bmag, double *flowUOrd, double *vtparSqOrd, double *vtperpSqOrd, double *fMFacOrd, double *bmagOrd);
void GkBiMaxwellianOnBasisGauss1x2vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtparSqOrd, const double *vtperpSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkBiMaxwellianOnBasisGauss1x2vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtparSq, const double *vtperpSq, const double *bmag, double *flowUOrd, double *vtparSqOrd, double *vtperpSqOrd, double *fMFacOrd, double *bmagOrd);
void GkBiMaxwellianOnBasisGauss1x2vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtparSqOrd, const double *vtperpSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);

void MaxwellianInnerLoop(double *n, double *u, double *vtparSq, double *vtperpSq, double *bmag, double m_, double *fItr, double *weights, double *dz, double *zc, double *ordinates, double *basisAtOrdinates, double *phaseToConfOrdMap, int numPhaseBasis, int numConfOrds, int numPhaseOrds, int numConfDims, int numPhaseDims);
void GkMaxwellianInnerLoop(double *n, double *u, double *vtparSq, double *vtperpSq, double *bmag, double m_, double *fItr, double *weights, double *dz, double *zc, double *ordinates, double *basisAtOrdinates, double *phaseToConfOrdMap, int numPhaseBasis, int numConfOrds, int numPhaseOrds, int numConfDims, int numPhaseDims);
} 

#endif 
