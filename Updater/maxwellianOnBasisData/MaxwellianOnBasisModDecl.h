#ifndef MAXWELLIANONBASIS_MOD_DECL_H 
#define MAXWELLIANONBASIS_MOD_DECL_H 
#include <cmath> 

extern "C" { 

void MaxwellianOnBasisGauss1x1vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x1vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x1vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x1vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x1vSer_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x1vSer_P3_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x2vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x2vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x2vSer_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x2vSer_P3_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x3vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x3vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x3vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x3vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x3vSer_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x3vSer_P3_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss2x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss2x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss2x2vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss2x2vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss2x2vSer_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss2x2vSer_P3_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss2x3vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss2x3vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss2x3vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss2x3vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss3x3vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss3x3vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x1vMax_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x1vMax_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x1vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x1vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x2vMax_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x2vMax_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss2x2vMax_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss2x2vMax_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss2x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss2x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x1vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x1vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x1vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x1vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x2vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x2vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss1x2vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss1x2vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss2x2vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss2x2vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);
void MaxwellianOnBasisGauss2x2vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd);
void MaxwellianOnBasisGauss2x2vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dvx, double *fMOut);

} 

#endif 
