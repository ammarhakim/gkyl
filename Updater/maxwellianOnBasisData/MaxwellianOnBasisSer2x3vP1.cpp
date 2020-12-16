#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss2x3vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[4];
  m0Ord[0] = 0.4999999999999999*den[3]-0.4999999999999999*(den[2]+den[1])+0.5*den[0]; 
  m0Ord[1] = (-0.4999999999999999*den[3])+0.4999999999999999*den[2]-0.4999999999999999*den[1]+0.5*den[0]; 
  m0Ord[2] = (-0.4999999999999999*(den[3]+den[2]))+0.4999999999999999*den[1]+0.5*den[0]; 
  m0Ord[3] = 0.4999999999999999*(den[3]+den[2]+den[1])+0.5*den[0]; 

  flowUOrd[0] = 0.4999999999999999*flowU[3]-0.4999999999999999*(flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[1] = (-0.4999999999999999*flowU[3])+0.4999999999999999*flowU[2]-0.4999999999999999*flowU[1]+0.5*flowU[0]; 
  flowUOrd[2] = (-0.4999999999999999*(flowU[3]+flowU[2]))+0.4999999999999999*flowU[1]+0.5*flowU[0]; 
  flowUOrd[3] = 0.4999999999999999*(flowU[3]+flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[4] = 0.4999999999999999*flowU[7]-0.4999999999999999*(flowU[6]+flowU[5])+0.5*flowU[4]; 
  flowUOrd[5] = (-0.4999999999999999*flowU[7])+0.4999999999999999*flowU[6]-0.4999999999999999*flowU[5]+0.5*flowU[4]; 
  flowUOrd[6] = (-0.4999999999999999*(flowU[7]+flowU[6]))+0.4999999999999999*flowU[5]+0.5*flowU[4]; 
  flowUOrd[7] = 0.4999999999999999*(flowU[7]+flowU[6]+flowU[5])+0.5*flowU[4]; 
  flowUOrd[8] = 0.4999999999999999*flowU[11]-0.4999999999999999*(flowU[10]+flowU[9])+0.5*flowU[8]; 
  flowUOrd[9] = (-0.4999999999999999*flowU[11])+0.4999999999999999*flowU[10]-0.4999999999999999*flowU[9]+0.5*flowU[8]; 
  flowUOrd[10] = (-0.4999999999999999*(flowU[11]+flowU[10]))+0.4999999999999999*flowU[9]+0.5*flowU[8]; 
  flowUOrd[11] = 0.4999999999999999*(flowU[11]+flowU[10]+flowU[9])+0.5*flowU[8]; 

  vtSqOrd[0] = 0.4999999999999999*vtSq[3]-0.4999999999999999*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 
  vtSqOrd[1] = (-0.4999999999999999*vtSq[3])+0.4999999999999999*vtSq[2]-0.4999999999999999*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[2] = (-0.4999999999999999*(vtSq[3]+vtSq[2]))+0.4999999999999999*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[3] = 0.4999999999999999*(vtSq[3]+vtSq[2]+vtSq[1])+0.5*vtSq[0]; 

  if ((vtSqOrd[0] <= 0.0) || (m0Ord[0] <= 0.0))
    fMFacOrd[0] = 0.;
  else
    fMFacOrd[0] = m0Ord[0]/std::pow(2.506628274631001*sqrt(vtSqOrd[0]),3.0); 
  if ((vtSqOrd[1] <= 0.0) || (m0Ord[1] <= 0.0))
    fMFacOrd[1] = 0.;
  else
    fMFacOrd[1] = m0Ord[1]/std::pow(2.506628274631001*sqrt(vtSqOrd[1]),3.0); 
  if ((vtSqOrd[2] <= 0.0) || (m0Ord[2] <= 0.0))
    fMFacOrd[2] = 0.;
  else
    fMFacOrd[2] = m0Ord[2]/std::pow(2.506628274631001*sqrt(vtSqOrd[2]),3.0); 
  if ((vtSqOrd[3] <= 0.0) || (m0Ord[3] <= 0.0))
    fMFacOrd[3] = 0.;
  else
    fMFacOrd[3] = m0Ord[3]/std::pow(2.506628274631001*sqrt(vtSqOrd[3]),3.0); 

}

void MaxwellianOnBasisGauss2x3vSerUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[4];
  m0Ord[0] = 0.4999999999999999*den[3]-0.4999999999999999*(den[2]+den[1])+0.5*den[0]; 
  m0Ord[1] = (-0.4999999999999999*den[3])+0.4999999999999999*den[2]-0.4999999999999999*den[1]+0.5*den[0]; 
  m0Ord[2] = (-0.4999999999999999*(den[3]+den[2]))+0.4999999999999999*den[1]+0.5*den[0]; 
  m0Ord[3] = 0.4999999999999999*(den[3]+den[2]+den[1])+0.5*den[0]; 

  flowUOrd[0] = 0.0; 
  flowUOrd[1] = 0.0; 
  flowUOrd[2] = 0.0; 
  flowUOrd[3] = 0.0; 
  flowUOrd[4] = 0.0; 
  flowUOrd[5] = 0.0; 
  flowUOrd[6] = 0.0; 
  flowUOrd[7] = 0.0; 
  flowUOrd[8] = 0.4999999999999999*flowU[3]-0.4999999999999999*(flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[9] = (-0.4999999999999999*flowU[3])+0.4999999999999999*flowU[2]-0.4999999999999999*flowU[1]+0.5*flowU[0]; 
  flowUOrd[10] = (-0.4999999999999999*(flowU[3]+flowU[2]))+0.4999999999999999*flowU[1]+0.5*flowU[0]; 
  flowUOrd[11] = 0.4999999999999999*(flowU[3]+flowU[2]+flowU[1])+0.5*flowU[0]; 

  vtSqOrd[0] = 0.4999999999999999*vtSq[3]-0.4999999999999999*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 
  vtSqOrd[1] = (-0.4999999999999999*vtSq[3])+0.4999999999999999*vtSq[2]-0.4999999999999999*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[2] = (-0.4999999999999999*(vtSq[3]+vtSq[2]))+0.4999999999999999*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[3] = 0.4999999999999999*(vtSq[3]+vtSq[2]+vtSq[1])+0.5*vtSq[0]; 

  if ((vtSqOrd[0] <= 0.0) || (m0Ord[0] <= 0.0))
    fMFacOrd[0] = 0.;
  else
    fMFacOrd[0] = m0Ord[0]/std::pow(2.506628274631001*sqrt(vtSqOrd[0]),3.0); 
  if ((vtSqOrd[1] <= 0.0) || (m0Ord[1] <= 0.0))
    fMFacOrd[1] = 0.;
  else
    fMFacOrd[1] = m0Ord[1]/std::pow(2.506628274631001*sqrt(vtSqOrd[1]),3.0); 
  if ((vtSqOrd[2] <= 0.0) || (m0Ord[2] <= 0.0))
    fMFacOrd[2] = 0.;
  else
    fMFacOrd[2] = m0Ord[2]/std::pow(2.506628274631001*sqrt(vtSqOrd[2]),3.0); 
  if ((vtSqOrd[3] <= 0.0) || (m0Ord[3] <= 0.0))
    fMFacOrd[3] = 0.;
  else
    fMFacOrd[3] = m0Ord[3]/std::pow(2.506628274631001*sqrt(vtSqOrd[3]),3.0); 

}

void MaxwellianOnBasisGauss2x3vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[32];
  if ((vtSqOrd[0] <= 0.0) || (fMFacOrd[0] <= 0.0)) {
    fMquad[0] = 0;
    fMquad[1] = 0;
    fMquad[2] = 0;
    fMquad[3] = 0;
    fMquad[4] = 0;
    fMquad[5] = 0;
    fMquad[6] = 0;
    fMquad[7] = 0;
  } else {
    fMquad[0] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[8])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[4])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[1] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[8])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[4])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[2] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[8])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[4])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[3] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[8])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[4])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[4] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[8])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[4])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[5] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[8])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[4])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[6] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[8])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[4])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[7] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[8])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[4])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  };
  if ((vtSqOrd[1] <= 0.0) || (fMFacOrd[1] <= 0.0)) {
    fMquad[8] = 0;
    fMquad[9] = 0;
    fMquad[10] = 0;
    fMquad[11] = 0;
    fMquad[12] = 0;
    fMquad[13] = 0;
    fMquad[14] = 0;
    fMquad[15] = 0;
  } else {
    fMquad[8] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[5])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[9] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[5])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[10] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[5])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[11] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[5])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[12] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[5])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[13] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[5])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[14] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[5])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[15] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[5])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  };
  if ((vtSqOrd[2] <= 0.0) || (fMFacOrd[2] <= 0.0)) {
    fMquad[16] = 0;
    fMquad[17] = 0;
    fMquad[18] = 0;
    fMquad[19] = 0;
    fMquad[20] = 0;
    fMquad[21] = 0;
    fMquad[22] = 0;
    fMquad[23] = 0;
  } else {
    fMquad[16] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[6])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[17] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[6])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[18] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[6])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[19] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[6])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[20] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[6])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[21] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[6])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[22] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[6])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[23] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[6])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
  };
  if ((vtSqOrd[3] <= 0.0) || (fMFacOrd[3] <= 0.0)) {
    fMquad[24] = 0;
    fMquad[25] = 0;
    fMquad[26] = 0;
    fMquad[27] = 0;
    fMquad[28] = 0;
    fMquad[29] = 0;
    fMquad[30] = 0;
    fMquad[31] = 0;
  } else {
    fMquad[24] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[7])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[25] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[7])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[26] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[7])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[27] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[7])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[28] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[7])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[29] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[7])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[30] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[4]-0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[7])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[31] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[4]+0.2886751345948129*dxv[4],2.0)+std::pow((-1.0*flowUOrd[7])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
  };

  fMOut[0] = 0.1767766952966368*(fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28]+fMquad[27]+fMquad[26]+fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[1] = 0.1767766952966368*(fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28]+fMquad[27]+fMquad[26]+fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16])-0.1767766952966368*(fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[2] = 0.1767766952966368*(fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28]+fMquad[27]+fMquad[26]+fMquad[25]+fMquad[24])-0.1767766952966368*(fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16])+0.1767766952966368*(fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8])-0.1767766952966368*(fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[3] = 0.1767766952966368*(fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28])-0.1767766952966368*(fMquad[27]+fMquad[26]+fMquad[25]+fMquad[24])+0.1767766952966368*(fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20])-0.1767766952966368*(fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16])+0.1767766952966368*(fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12])-0.1767766952966368*(fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8])+0.1767766952966368*(fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4])-0.1767766952966368*(fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[4] = 0.1767766952966368*(fMquad[31]+fMquad[30])-0.1767766952966368*(fMquad[29]+fMquad[28])+0.1767766952966368*(fMquad[27]+fMquad[26])-0.1767766952966368*(fMquad[25]+fMquad[24])+0.1767766952966368*(fMquad[23]+fMquad[22])-0.1767766952966368*(fMquad[21]+fMquad[20])+0.1767766952966368*(fMquad[19]+fMquad[18])-0.1767766952966368*(fMquad[17]+fMquad[16])+0.1767766952966368*(fMquad[15]+fMquad[14])-0.1767766952966368*(fMquad[13]+fMquad[12])+0.1767766952966368*(fMquad[11]+fMquad[10])-0.1767766952966368*(fMquad[9]+fMquad[8])+0.1767766952966368*(fMquad[7]+fMquad[6])-0.1767766952966368*(fMquad[5]+fMquad[4])+0.1767766952966368*(fMquad[3]+fMquad[2])-0.1767766952966368*(fMquad[1]+fMquad[0]); 
  fMOut[5] = 0.1767766952966368*fMquad[31]-0.1767766952966368*fMquad[30]+0.1767766952966368*fMquad[29]-0.1767766952966368*fMquad[28]+0.1767766952966368*fMquad[27]-0.1767766952966368*fMquad[26]+0.1767766952966368*fMquad[25]-0.1767766952966368*fMquad[24]+0.1767766952966368*fMquad[23]-0.1767766952966368*fMquad[22]+0.1767766952966368*fMquad[21]-0.1767766952966368*fMquad[20]+0.1767766952966368*fMquad[19]-0.1767766952966368*fMquad[18]+0.1767766952966368*fMquad[17]-0.1767766952966368*fMquad[16]+0.1767766952966368*fMquad[15]-0.1767766952966368*fMquad[14]+0.1767766952966368*fMquad[13]-0.1767766952966368*fMquad[12]+0.1767766952966368*fMquad[11]-0.1767766952966368*fMquad[10]+0.1767766952966368*fMquad[9]-0.1767766952966368*fMquad[8]+0.1767766952966368*fMquad[7]-0.1767766952966368*fMquad[6]+0.1767766952966368*fMquad[5]-0.1767766952966368*fMquad[4]+0.1767766952966368*fMquad[3]-0.1767766952966368*fMquad[2]+0.1767766952966368*fMquad[1]-0.1767766952966368*fMquad[0]; 
  fMOut[6] = 0.1767766952966368*(fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28]+fMquad[27]+fMquad[26]+fMquad[25]+fMquad[24])-0.1767766952966368*(fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8])+0.1767766952966368*(fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[7] = 0.1767766952966368*(fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28])-0.1767766952966368*(fMquad[27]+fMquad[26]+fMquad[25]+fMquad[24])+0.1767766952966368*(fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20])-0.1767766952966368*(fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12])+0.1767766952966368*(fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8])-0.1767766952966368*(fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4])+0.1767766952966368*(fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[8] = 0.1767766952966368*(fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28])-0.1767766952966368*(fMquad[27]+fMquad[26]+fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20])+0.1767766952966368*(fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12])-0.1767766952966368*(fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4])+0.1767766952966368*(fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[9] = 0.1767766952966368*(fMquad[31]+fMquad[30])-0.1767766952966368*(fMquad[29]+fMquad[28])+0.1767766952966368*(fMquad[27]+fMquad[26])-0.1767766952966368*(fMquad[25]+fMquad[24])+0.1767766952966368*(fMquad[23]+fMquad[22])-0.1767766952966368*(fMquad[21]+fMquad[20])+0.1767766952966368*(fMquad[19]+fMquad[18])-0.1767766952966368*(fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14])+0.1767766952966368*(fMquad[13]+fMquad[12])-0.1767766952966368*(fMquad[11]+fMquad[10])+0.1767766952966368*(fMquad[9]+fMquad[8])-0.1767766952966368*(fMquad[7]+fMquad[6])+0.1767766952966368*(fMquad[5]+fMquad[4])-0.1767766952966368*(fMquad[3]+fMquad[2])+0.1767766952966368*(fMquad[1]+fMquad[0]); 
  fMOut[10] = 0.1767766952966368*(fMquad[31]+fMquad[30])-0.1767766952966368*(fMquad[29]+fMquad[28])+0.1767766952966368*(fMquad[27]+fMquad[26])-0.1767766952966368*(fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22])+0.1767766952966368*(fMquad[21]+fMquad[20])-0.1767766952966368*(fMquad[19]+fMquad[18])+0.1767766952966368*(fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14])-0.1767766952966368*(fMquad[13]+fMquad[12])+0.1767766952966368*(fMquad[11]+fMquad[10])-0.1767766952966368*(fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6])+0.1767766952966368*(fMquad[5]+fMquad[4])-0.1767766952966368*(fMquad[3]+fMquad[2])+0.1767766952966368*(fMquad[1]+fMquad[0]); 
  fMOut[11] = 0.1767766952966368*(fMquad[31]+fMquad[30])-0.1767766952966368*(fMquad[29]+fMquad[28]+fMquad[27]+fMquad[26])+0.1767766952966368*(fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22])-0.1767766952966368*(fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18])+0.1767766952966368*(fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14])-0.1767766952966368*(fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10])+0.1767766952966368*(fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6])-0.1767766952966368*(fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2])+0.1767766952966368*(fMquad[1]+fMquad[0]); 
  fMOut[12] = 0.1767766952966368*fMquad[31]-0.1767766952966368*fMquad[30]+0.1767766952966368*fMquad[29]-0.1767766952966368*fMquad[28]+0.1767766952966368*fMquad[27]-0.1767766952966368*fMquad[26]+0.1767766952966368*fMquad[25]-0.1767766952966368*fMquad[24]+0.1767766952966368*fMquad[23]-0.1767766952966368*fMquad[22]+0.1767766952966368*fMquad[21]-0.1767766952966368*fMquad[20]+0.1767766952966368*fMquad[19]-0.1767766952966368*fMquad[18]+0.1767766952966368*fMquad[17]-0.1767766952966368*(fMquad[16]+fMquad[15])+0.1767766952966368*fMquad[14]-0.1767766952966368*fMquad[13]+0.1767766952966368*fMquad[12]-0.1767766952966368*fMquad[11]+0.1767766952966368*fMquad[10]-0.1767766952966368*fMquad[9]+0.1767766952966368*fMquad[8]-0.1767766952966368*fMquad[7]+0.1767766952966368*fMquad[6]-0.1767766952966368*fMquad[5]+0.1767766952966368*fMquad[4]-0.1767766952966368*fMquad[3]+0.1767766952966368*fMquad[2]-0.1767766952966368*fMquad[1]+0.1767766952966368*fMquad[0]; 
  fMOut[13] = 0.1767766952966368*fMquad[31]-0.1767766952966368*fMquad[30]+0.1767766952966368*fMquad[29]-0.1767766952966368*fMquad[28]+0.1767766952966368*fMquad[27]-0.1767766952966368*fMquad[26]+0.1767766952966368*fMquad[25]-0.1767766952966368*(fMquad[24]+fMquad[23])+0.1767766952966368*fMquad[22]-0.1767766952966368*fMquad[21]+0.1767766952966368*fMquad[20]-0.1767766952966368*fMquad[19]+0.1767766952966368*fMquad[18]-0.1767766952966368*fMquad[17]+0.1767766952966368*(fMquad[16]+fMquad[15])-0.1767766952966368*fMquad[14]+0.1767766952966368*fMquad[13]-0.1767766952966368*fMquad[12]+0.1767766952966368*fMquad[11]-0.1767766952966368*fMquad[10]+0.1767766952966368*fMquad[9]-0.1767766952966368*(fMquad[8]+fMquad[7])+0.1767766952966368*fMquad[6]-0.1767766952966368*fMquad[5]+0.1767766952966368*fMquad[4]-0.1767766952966368*fMquad[3]+0.1767766952966368*fMquad[2]-0.1767766952966368*fMquad[1]+0.1767766952966368*fMquad[0]; 
  fMOut[14] = 0.1767766952966368*fMquad[31]-0.1767766952966368*fMquad[30]+0.1767766952966368*fMquad[29]-0.1767766952966368*(fMquad[28]+fMquad[27])+0.1767766952966368*fMquad[26]-0.1767766952966368*fMquad[25]+0.1767766952966368*(fMquad[24]+fMquad[23])-0.1767766952966368*fMquad[22]+0.1767766952966368*fMquad[21]-0.1767766952966368*(fMquad[20]+fMquad[19])+0.1767766952966368*fMquad[18]-0.1767766952966368*fMquad[17]+0.1767766952966368*(fMquad[16]+fMquad[15])-0.1767766952966368*fMquad[14]+0.1767766952966368*fMquad[13]-0.1767766952966368*(fMquad[12]+fMquad[11])+0.1767766952966368*fMquad[10]-0.1767766952966368*fMquad[9]+0.1767766952966368*(fMquad[8]+fMquad[7])-0.1767766952966368*fMquad[6]+0.1767766952966368*fMquad[5]-0.1767766952966368*(fMquad[4]+fMquad[3])+0.1767766952966368*fMquad[2]-0.1767766952966368*fMquad[1]+0.1767766952966368*fMquad[0]; 
  fMOut[15] = 0.1767766952966368*fMquad[31]-0.1767766952966368*(fMquad[30]+fMquad[29])+0.1767766952966368*(fMquad[28]+fMquad[27])-0.1767766952966368*(fMquad[26]+fMquad[25])+0.1767766952966368*(fMquad[24]+fMquad[23])-0.1767766952966368*(fMquad[22]+fMquad[21])+0.1767766952966368*(fMquad[20]+fMquad[19])-0.1767766952966368*(fMquad[18]+fMquad[17])+0.1767766952966368*(fMquad[16]+fMquad[15])-0.1767766952966368*(fMquad[14]+fMquad[13])+0.1767766952966368*(fMquad[12]+fMquad[11])-0.1767766952966368*(fMquad[10]+fMquad[9])+0.1767766952966368*(fMquad[8]+fMquad[7])-0.1767766952966368*(fMquad[6]+fMquad[5])+0.1767766952966368*(fMquad[4]+fMquad[3])-0.1767766952966368*(fMquad[2]+fMquad[1])+0.1767766952966368*fMquad[0]; 
  fMOut[16] = 0.1767766952966368*(fMquad[31]+fMquad[30]+fMquad[29]+fMquad[28])-0.1767766952966368*(fMquad[27]+fMquad[26]+fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22]+fMquad[21]+fMquad[20])+0.1767766952966368*(fMquad[19]+fMquad[18]+fMquad[17]+fMquad[16])-0.1767766952966368*(fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12])+0.1767766952966368*(fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4])-0.1767766952966368*(fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[17] = 0.1767766952966368*(fMquad[31]+fMquad[30])-0.1767766952966368*(fMquad[29]+fMquad[28])+0.1767766952966368*(fMquad[27]+fMquad[26])-0.1767766952966368*(fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22])+0.1767766952966368*(fMquad[21]+fMquad[20])-0.1767766952966368*(fMquad[19]+fMquad[18])+0.1767766952966368*(fMquad[17]+fMquad[16])-0.1767766952966368*(fMquad[15]+fMquad[14])+0.1767766952966368*(fMquad[13]+fMquad[12])-0.1767766952966368*(fMquad[11]+fMquad[10])+0.1767766952966368*(fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6])-0.1767766952966368*(fMquad[5]+fMquad[4])+0.1767766952966368*(fMquad[3]+fMquad[2])-0.1767766952966368*(fMquad[1]+fMquad[0]); 
  fMOut[18] = 0.1767766952966368*(fMquad[31]+fMquad[30])-0.1767766952966368*(fMquad[29]+fMquad[28]+fMquad[27]+fMquad[26])+0.1767766952966368*(fMquad[25]+fMquad[24]+fMquad[23]+fMquad[22])-0.1767766952966368*(fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18])+0.1767766952966368*(fMquad[17]+fMquad[16])-0.1767766952966368*(fMquad[15]+fMquad[14])+0.1767766952966368*(fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10])-0.1767766952966368*(fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6])+0.1767766952966368*(fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2])-0.1767766952966368*(fMquad[1]+fMquad[0]); 
  fMOut[19] = 0.1767766952966368*(fMquad[31]+fMquad[30])-0.1767766952966368*(fMquad[29]+fMquad[28]+fMquad[27]+fMquad[26])+0.1767766952966368*(fMquad[25]+fMquad[24])-0.1767766952966368*(fMquad[23]+fMquad[22])+0.1767766952966368*(fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18])-0.1767766952966368*(fMquad[17]+fMquad[16])+0.1767766952966368*(fMquad[15]+fMquad[14])-0.1767766952966368*(fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10])+0.1767766952966368*(fMquad[9]+fMquad[8])-0.1767766952966368*(fMquad[7]+fMquad[6])+0.1767766952966368*(fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2])-0.1767766952966368*(fMquad[1]+fMquad[0]); 
  fMOut[20] = 0.1767766952966368*fMquad[31]-0.1767766952966368*fMquad[30]+0.1767766952966368*fMquad[29]-0.1767766952966368*fMquad[28]+0.1767766952966368*fMquad[27]-0.1767766952966368*fMquad[26]+0.1767766952966368*fMquad[25]-0.1767766952966368*(fMquad[24]+fMquad[23])+0.1767766952966368*fMquad[22]-0.1767766952966368*fMquad[21]+0.1767766952966368*fMquad[20]-0.1767766952966368*fMquad[19]+0.1767766952966368*fMquad[18]-0.1767766952966368*fMquad[17]+0.1767766952966368*fMquad[16]-0.1767766952966368*fMquad[15]+0.1767766952966368*fMquad[14]-0.1767766952966368*fMquad[13]+0.1767766952966368*fMquad[12]-0.1767766952966368*fMquad[11]+0.1767766952966368*fMquad[10]-0.1767766952966368*fMquad[9]+0.1767766952966368*(fMquad[8]+fMquad[7])-0.1767766952966368*fMquad[6]+0.1767766952966368*fMquad[5]-0.1767766952966368*fMquad[4]+0.1767766952966368*fMquad[3]-0.1767766952966368*fMquad[2]+0.1767766952966368*fMquad[1]-0.1767766952966368*fMquad[0]; 
  fMOut[21] = 0.1767766952966368*fMquad[31]-0.1767766952966368*fMquad[30]+0.1767766952966368*fMquad[29]-0.1767766952966368*(fMquad[28]+fMquad[27])+0.1767766952966368*fMquad[26]-0.1767766952966368*fMquad[25]+0.1767766952966368*(fMquad[24]+fMquad[23])-0.1767766952966368*fMquad[22]+0.1767766952966368*fMquad[21]-0.1767766952966368*(fMquad[20]+fMquad[19])+0.1767766952966368*fMquad[18]-0.1767766952966368*fMquad[17]+0.1767766952966368*fMquad[16]-0.1767766952966368*fMquad[15]+0.1767766952966368*fMquad[14]-0.1767766952966368*fMquad[13]+0.1767766952966368*(fMquad[12]+fMquad[11])-0.1767766952966368*fMquad[10]+0.1767766952966368*fMquad[9]-0.1767766952966368*(fMquad[8]+fMquad[7])+0.1767766952966368*fMquad[6]-0.1767766952966368*fMquad[5]+0.1767766952966368*(fMquad[4]+fMquad[3])-0.1767766952966368*fMquad[2]+0.1767766952966368*fMquad[1]-0.1767766952966368*fMquad[0]; 
  fMOut[22] = 0.1767766952966368*fMquad[31]-0.1767766952966368*fMquad[30]+0.1767766952966368*fMquad[29]-0.1767766952966368*(fMquad[28]+fMquad[27])+0.1767766952966368*fMquad[26]-0.1767766952966368*fMquad[25]+0.1767766952966368*fMquad[24]-0.1767766952966368*fMquad[23]+0.1767766952966368*fMquad[22]-0.1767766952966368*fMquad[21]+0.1767766952966368*(fMquad[20]+fMquad[19])-0.1767766952966368*fMquad[18]+0.1767766952966368*fMquad[17]-0.1767766952966368*fMquad[16]+0.1767766952966368*fMquad[15]-0.1767766952966368*fMquad[14]+0.1767766952966368*fMquad[13]-0.1767766952966368*(fMquad[12]+fMquad[11])+0.1767766952966368*fMquad[10]-0.1767766952966368*fMquad[9]+0.1767766952966368*fMquad[8]-0.1767766952966368*fMquad[7]+0.1767766952966368*fMquad[6]-0.1767766952966368*fMquad[5]+0.1767766952966368*(fMquad[4]+fMquad[3])-0.1767766952966368*fMquad[2]+0.1767766952966368*fMquad[1]-0.1767766952966368*fMquad[0]; 
  fMOut[23] = 0.1767766952966368*fMquad[31]-0.1767766952966368*(fMquad[30]+fMquad[29])+0.1767766952966368*(fMquad[28]+fMquad[27])-0.1767766952966368*(fMquad[26]+fMquad[25])+0.1767766952966368*(fMquad[24]+fMquad[23])-0.1767766952966368*(fMquad[22]+fMquad[21])+0.1767766952966368*(fMquad[20]+fMquad[19])-0.1767766952966368*(fMquad[18]+fMquad[17])+0.1767766952966368*fMquad[16]-0.1767766952966368*fMquad[15]+0.1767766952966368*(fMquad[14]+fMquad[13])-0.1767766952966368*(fMquad[12]+fMquad[11])+0.1767766952966368*(fMquad[10]+fMquad[9])-0.1767766952966368*(fMquad[8]+fMquad[7])+0.1767766952966368*(fMquad[6]+fMquad[5])-0.1767766952966368*(fMquad[4]+fMquad[3])+0.1767766952966368*(fMquad[2]+fMquad[1])-0.1767766952966368*fMquad[0]; 
  fMOut[24] = 0.1767766952966368*fMquad[31]-0.1767766952966368*(fMquad[30]+fMquad[29])+0.1767766952966368*(fMquad[28]+fMquad[27])-0.1767766952966368*(fMquad[26]+fMquad[25])+0.1767766952966368*fMquad[24]-0.1767766952966368*fMquad[23]+0.1767766952966368*(fMquad[22]+fMquad[21])-0.1767766952966368*(fMquad[20]+fMquad[19])+0.1767766952966368*(fMquad[18]+fMquad[17])-0.1767766952966368*fMquad[16]+0.1767766952966368*fMquad[15]-0.1767766952966368*(fMquad[14]+fMquad[13])+0.1767766952966368*(fMquad[12]+fMquad[11])-0.1767766952966368*(fMquad[10]+fMquad[9])+0.1767766952966368*fMquad[8]-0.1767766952966368*fMquad[7]+0.1767766952966368*(fMquad[6]+fMquad[5])-0.1767766952966368*(fMquad[4]+fMquad[3])+0.1767766952966368*(fMquad[2]+fMquad[1])-0.1767766952966368*fMquad[0]; 
  fMOut[25] = 0.1767766952966368*fMquad[31]-0.1767766952966368*(fMquad[30]+fMquad[29])+0.1767766952966368*fMquad[28]-0.1767766952966368*fMquad[27]+0.1767766952966368*(fMquad[26]+fMquad[25])-0.1767766952966368*fMquad[24]+0.1767766952966368*fMquad[23]-0.1767766952966368*(fMquad[22]+fMquad[21])+0.1767766952966368*fMquad[20]-0.1767766952966368*fMquad[19]+0.1767766952966368*(fMquad[18]+fMquad[17])-0.1767766952966368*fMquad[16]+0.1767766952966368*fMquad[15]-0.1767766952966368*(fMquad[14]+fMquad[13])+0.1767766952966368*fMquad[12]-0.1767766952966368*fMquad[11]+0.1767766952966368*(fMquad[10]+fMquad[9])-0.1767766952966368*fMquad[8]+0.1767766952966368*fMquad[7]-0.1767766952966368*(fMquad[6]+fMquad[5])+0.1767766952966368*fMquad[4]-0.1767766952966368*fMquad[3]+0.1767766952966368*(fMquad[2]+fMquad[1])-0.1767766952966368*fMquad[0]; 
  fMOut[26] = 0.1767766952966368*(fMquad[31]+fMquad[30])-0.1767766952966368*(fMquad[29]+fMquad[28]+fMquad[27]+fMquad[26])+0.1767766952966368*(fMquad[25]+fMquad[24])-0.1767766952966368*(fMquad[23]+fMquad[22])+0.1767766952966368*(fMquad[21]+fMquad[20]+fMquad[19]+fMquad[18])-0.1767766952966368*(fMquad[17]+fMquad[16]+fMquad[15]+fMquad[14])+0.1767766952966368*(fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10])-0.1767766952966368*(fMquad[9]+fMquad[8])+0.1767766952966368*(fMquad[7]+fMquad[6])-0.1767766952966368*(fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2])+0.1767766952966368*(fMquad[1]+fMquad[0]); 
  fMOut[27] = 0.1767766952966368*fMquad[31]-0.1767766952966368*fMquad[30]+0.1767766952966368*fMquad[29]-0.1767766952966368*(fMquad[28]+fMquad[27])+0.1767766952966368*fMquad[26]-0.1767766952966368*fMquad[25]+0.1767766952966368*fMquad[24]-0.1767766952966368*fMquad[23]+0.1767766952966368*fMquad[22]-0.1767766952966368*fMquad[21]+0.1767766952966368*(fMquad[20]+fMquad[19])-0.1767766952966368*fMquad[18]+0.1767766952966368*fMquad[17]-0.1767766952966368*(fMquad[16]+fMquad[15])+0.1767766952966368*fMquad[14]-0.1767766952966368*fMquad[13]+0.1767766952966368*(fMquad[12]+fMquad[11])-0.1767766952966368*fMquad[10]+0.1767766952966368*fMquad[9]-0.1767766952966368*fMquad[8]+0.1767766952966368*fMquad[7]-0.1767766952966368*fMquad[6]+0.1767766952966368*fMquad[5]-0.1767766952966368*(fMquad[4]+fMquad[3])+0.1767766952966368*fMquad[2]-0.1767766952966368*fMquad[1]+0.1767766952966368*fMquad[0]; 
  fMOut[28] = 0.1767766952966368*fMquad[31]-0.1767766952966368*(fMquad[30]+fMquad[29])+0.1767766952966368*(fMquad[28]+fMquad[27])-0.1767766952966368*(fMquad[26]+fMquad[25])+0.1767766952966368*fMquad[24]-0.1767766952966368*fMquad[23]+0.1767766952966368*(fMquad[22]+fMquad[21])-0.1767766952966368*(fMquad[20]+fMquad[19])+0.1767766952966368*(fMquad[18]+fMquad[17])-0.1767766952966368*(fMquad[16]+fMquad[15])+0.1767766952966368*(fMquad[14]+fMquad[13])-0.1767766952966368*(fMquad[12]+fMquad[11])+0.1767766952966368*(fMquad[10]+fMquad[9])-0.1767766952966368*fMquad[8]+0.1767766952966368*fMquad[7]-0.1767766952966368*(fMquad[6]+fMquad[5])+0.1767766952966368*(fMquad[4]+fMquad[3])-0.1767766952966368*(fMquad[2]+fMquad[1])+0.1767766952966368*fMquad[0]; 
  fMOut[29] = 0.1767766952966368*fMquad[31]-0.1767766952966368*(fMquad[30]+fMquad[29])+0.1767766952966368*fMquad[28]-0.1767766952966368*fMquad[27]+0.1767766952966368*(fMquad[26]+fMquad[25])-0.1767766952966368*fMquad[24]+0.1767766952966368*fMquad[23]-0.1767766952966368*(fMquad[22]+fMquad[21])+0.1767766952966368*fMquad[20]-0.1767766952966368*fMquad[19]+0.1767766952966368*(fMquad[18]+fMquad[17])-0.1767766952966368*(fMquad[16]+fMquad[15])+0.1767766952966368*(fMquad[14]+fMquad[13])-0.1767766952966368*fMquad[12]+0.1767766952966368*fMquad[11]-0.1767766952966368*(fMquad[10]+fMquad[9])+0.1767766952966368*fMquad[8]-0.1767766952966368*fMquad[7]+0.1767766952966368*(fMquad[6]+fMquad[5])-0.1767766952966368*fMquad[4]+0.1767766952966368*fMquad[3]-0.1767766952966368*(fMquad[2]+fMquad[1])+0.1767766952966368*fMquad[0]; 
  fMOut[30] = 0.1767766952966368*fMquad[31]-0.1767766952966368*(fMquad[30]+fMquad[29])+0.1767766952966368*fMquad[28]-0.1767766952966368*fMquad[27]+0.1767766952966368*(fMquad[26]+fMquad[25])-0.1767766952966368*(fMquad[24]+fMquad[23])+0.1767766952966368*(fMquad[22]+fMquad[21])-0.1767766952966368*fMquad[20]+0.1767766952966368*fMquad[19]-0.1767766952966368*(fMquad[18]+fMquad[17])+0.1767766952966368*(fMquad[16]+fMquad[15])-0.1767766952966368*(fMquad[14]+fMquad[13])+0.1767766952966368*fMquad[12]-0.1767766952966368*fMquad[11]+0.1767766952966368*(fMquad[10]+fMquad[9])-0.1767766952966368*(fMquad[8]+fMquad[7])+0.1767766952966368*(fMquad[6]+fMquad[5])-0.1767766952966368*fMquad[4]+0.1767766952966368*fMquad[3]-0.1767766952966368*(fMquad[2]+fMquad[1])+0.1767766952966368*fMquad[0]; 
  fMOut[31] = 0.1767766952966367*fMquad[31]-0.1767766952966367*(fMquad[30]+fMquad[29])+0.1767766952966367*fMquad[28]-0.1767766952966367*fMquad[27]+0.1767766952966367*(fMquad[26]+fMquad[25])-0.1767766952966367*(fMquad[24]+fMquad[23])+0.1767766952966367*(fMquad[22]+fMquad[21])-0.1767766952966367*fMquad[20]+0.1767766952966367*fMquad[19]-0.1767766952966367*(fMquad[18]+fMquad[17])+0.1767766952966367*fMquad[16]-0.1767766952966367*fMquad[15]+0.1767766952966367*(fMquad[14]+fMquad[13])-0.1767766952966367*fMquad[12]+0.1767766952966367*fMquad[11]-0.1767766952966367*(fMquad[10]+fMquad[9])+0.1767766952966367*(fMquad[8]+fMquad[7])-0.1767766952966367*(fMquad[6]+fMquad[5])+0.1767766952966367*fMquad[4]-0.1767766952966367*fMquad[3]+0.1767766952966367*(fMquad[2]+fMquad[1])-0.1767766952966367*fMquad[0]; 

}
