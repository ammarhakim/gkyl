#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss2x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[9];
  m0Ord[0] = 0.5*den[0]-0.5590169943749475*(den[5]+den[4]); 
  m0Ord[1] = 0.4472135954999581*den[5]-0.5590169943749475*den[4]-0.6708203932499369*den[2]+0.5*den[0]; 
  m0Ord[2] = 0.4472135954999581*den[5]-0.5590169943749475*den[4]+0.6708203932499369*den[2]+0.5*den[0]; 
  m0Ord[3] = (-0.5590169943749475*den[5])+0.4472135954999581*den[4]-0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[4] = 0.4472135954999581*(den[5]+den[4])+0.9*den[3]-0.6708203932499369*(den[2]+den[1])+0.5*den[0]; 
  m0Ord[5] = 0.4472135954999581*(den[5]+den[4])-0.9*den[3]+0.6708203932499369*den[2]-0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[6] = (-0.5590169943749475*den[5])+0.4472135954999581*den[4]+0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[7] = 0.4472135954999581*(den[5]+den[4])-0.9*den[3]-0.6708203932499369*den[2]+0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[8] = 0.4472135954999581*(den[5]+den[4])+0.9*den[3]+0.6708203932499369*(den[2]+den[1])+0.5*den[0]; 

  flowUOrd[0] = 0.5*flowU[0]-0.5590169943749475*(flowU[5]+flowU[4]); 
  flowUOrd[1] = 0.4472135954999581*flowU[5]-0.5590169943749475*flowU[4]-0.6708203932499369*flowU[2]+0.5*flowU[0]; 
  flowUOrd[2] = 0.4472135954999581*flowU[5]-0.5590169943749475*flowU[4]+0.6708203932499369*flowU[2]+0.5*flowU[0]; 
  flowUOrd[3] = (-0.5590169943749475*flowU[5])+0.4472135954999581*flowU[4]-0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[4] = 0.4472135954999581*(flowU[5]+flowU[4])+0.9*flowU[3]-0.6708203932499369*(flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[5] = 0.4472135954999581*(flowU[5]+flowU[4])-0.9*flowU[3]+0.6708203932499369*flowU[2]-0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[6] = (-0.5590169943749475*flowU[5])+0.4472135954999581*flowU[4]+0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[7] = 0.4472135954999581*(flowU[5]+flowU[4])-0.9*flowU[3]-0.6708203932499369*flowU[2]+0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[8] = 0.4472135954999581*(flowU[5]+flowU[4])+0.9*flowU[3]+0.6708203932499369*(flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[9] = 0.5*flowU[6]-0.5590169943749475*(flowU[11]+flowU[10]); 
  flowUOrd[10] = 0.4472135954999581*flowU[11]-0.5590169943749475*flowU[10]-0.6708203932499369*flowU[8]+0.5*flowU[6]; 
  flowUOrd[11] = 0.4472135954999581*flowU[11]-0.5590169943749475*flowU[10]+0.6708203932499369*flowU[8]+0.5*flowU[6]; 
  flowUOrd[12] = (-0.5590169943749475*flowU[11])+0.4472135954999581*flowU[10]-0.6708203932499369*flowU[7]+0.5*flowU[6]; 
  flowUOrd[13] = 0.4472135954999581*(flowU[11]+flowU[10])+0.9*flowU[9]-0.6708203932499369*(flowU[8]+flowU[7])+0.5*flowU[6]; 
  flowUOrd[14] = 0.4472135954999581*(flowU[11]+flowU[10])-0.9*flowU[9]+0.6708203932499369*flowU[8]-0.6708203932499369*flowU[7]+0.5*flowU[6]; 
  flowUOrd[15] = (-0.5590169943749475*flowU[11])+0.4472135954999581*flowU[10]+0.6708203932499369*flowU[7]+0.5*flowU[6]; 
  flowUOrd[16] = 0.4472135954999581*(flowU[11]+flowU[10])-0.9*flowU[9]-0.6708203932499369*flowU[8]+0.6708203932499369*flowU[7]+0.5*flowU[6]; 
  flowUOrd[17] = 0.4472135954999581*(flowU[11]+flowU[10])+0.9*flowU[9]+0.6708203932499369*(flowU[8]+flowU[7])+0.5*flowU[6]; 

  vtSqOrd[0] = 0.5*vtSq[0]-0.5590169943749475*(vtSq[5]+vtSq[4]); 
  vtSqOrd[1] = 0.4472135954999581*vtSq[5]-0.5590169943749475*vtSq[4]-0.6708203932499369*vtSq[2]+0.5*vtSq[0]; 
  vtSqOrd[2] = 0.4472135954999581*vtSq[5]-0.5590169943749475*vtSq[4]+0.6708203932499369*vtSq[2]+0.5*vtSq[0]; 
  vtSqOrd[3] = (-0.5590169943749475*vtSq[5])+0.4472135954999581*vtSq[4]-0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[4] = 0.4472135954999581*(vtSq[5]+vtSq[4])+0.9*vtSq[3]-0.6708203932499369*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 
  vtSqOrd[5] = 0.4472135954999581*(vtSq[5]+vtSq[4])-0.9*vtSq[3]+0.6708203932499369*vtSq[2]-0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[6] = (-0.5590169943749475*vtSq[5])+0.4472135954999581*vtSq[4]+0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[7] = 0.4472135954999581*(vtSq[5]+vtSq[4])-0.9*vtSq[3]-0.6708203932499369*vtSq[2]+0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[8] = 0.4472135954999581*(vtSq[5]+vtSq[4])+0.9*vtSq[3]+0.6708203932499369*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (0.1591549430918953*m0Ord[0])/vtSqOrd[0]; 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (0.1591549430918953*m0Ord[1])/vtSqOrd[1]; 
  else
    fMFacOrd[1] = 0.0;
  if ((vtSqOrd[2] > 0.) && (m0Ord[2] > 0.))
    fMFacOrd[2] = (0.1591549430918953*m0Ord[2])/vtSqOrd[2]; 
  else
    fMFacOrd[2] = 0.0;
  if ((vtSqOrd[3] > 0.) && (m0Ord[3] > 0.))
    fMFacOrd[3] = (0.1591549430918953*m0Ord[3])/vtSqOrd[3]; 
  else
    fMFacOrd[3] = 0.0;
  if ((vtSqOrd[4] > 0.) && (m0Ord[4] > 0.))
    fMFacOrd[4] = (0.1591549430918953*m0Ord[4])/vtSqOrd[4]; 
  else
    fMFacOrd[4] = 0.0;
  if ((vtSqOrd[5] > 0.) && (m0Ord[5] > 0.))
    fMFacOrd[5] = (0.1591549430918953*m0Ord[5])/vtSqOrd[5]; 
  else
    fMFacOrd[5] = 0.0;
  if ((vtSqOrd[6] > 0.) && (m0Ord[6] > 0.))
    fMFacOrd[6] = (0.1591549430918953*m0Ord[6])/vtSqOrd[6]; 
  else
    fMFacOrd[6] = 0.0;
  if ((vtSqOrd[7] > 0.) && (m0Ord[7] > 0.))
    fMFacOrd[7] = (0.1591549430918953*m0Ord[7])/vtSqOrd[7]; 
  else
    fMFacOrd[7] = 0.0;
  if ((vtSqOrd[8] > 0.) && (m0Ord[8] > 0.))
    fMFacOrd[8] = (0.1591549430918953*m0Ord[8])/vtSqOrd[8]; 
  else
    fMFacOrd[8] = 0.0;

}

void MaxwellianOnBasisGauss2x2vMaxUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[9];
  m0Ord[0] = 0.5*den[0]-0.5590169943749475*(den[5]+den[4]); 
  m0Ord[1] = 0.4472135954999581*den[5]-0.5590169943749475*den[4]-0.6708203932499369*den[2]+0.5*den[0]; 
  m0Ord[2] = 0.4472135954999581*den[5]-0.5590169943749475*den[4]+0.6708203932499369*den[2]+0.5*den[0]; 
  m0Ord[3] = (-0.5590169943749475*den[5])+0.4472135954999581*den[4]-0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[4] = 0.4472135954999581*(den[5]+den[4])+0.9*den[3]-0.6708203932499369*(den[2]+den[1])+0.5*den[0]; 
  m0Ord[5] = 0.4472135954999581*(den[5]+den[4])-0.9*den[3]+0.6708203932499369*den[2]-0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[6] = (-0.5590169943749475*den[5])+0.4472135954999581*den[4]+0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[7] = 0.4472135954999581*(den[5]+den[4])-0.9*den[3]-0.6708203932499369*den[2]+0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[8] = 0.4472135954999581*(den[5]+den[4])+0.9*den[3]+0.6708203932499369*(den[2]+den[1])+0.5*den[0]; 

  flowUOrd[0] = 0.0; 
  flowUOrd[1] = 0.0; 
  flowUOrd[2] = 0.0; 
  flowUOrd[3] = 0.0; 
  flowUOrd[4] = 0.0; 
  flowUOrd[5] = 0.0; 
  flowUOrd[6] = 0.0; 
  flowUOrd[7] = 0.0; 
  flowUOrd[8] = 0.0; 
  flowUOrd[9] = 0.5*flowU[0]-0.5590169943749475*(flowU[5]+flowU[4]); 
  flowUOrd[10] = 0.4472135954999581*flowU[5]-0.5590169943749475*flowU[4]-0.6708203932499369*flowU[2]+0.5*flowU[0]; 
  flowUOrd[11] = 0.4472135954999581*flowU[5]-0.5590169943749475*flowU[4]+0.6708203932499369*flowU[2]+0.5*flowU[0]; 
  flowUOrd[12] = (-0.5590169943749475*flowU[5])+0.4472135954999581*flowU[4]-0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[13] = 0.4472135954999581*(flowU[5]+flowU[4])+0.9*flowU[3]-0.6708203932499369*(flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[14] = 0.4472135954999581*(flowU[5]+flowU[4])-0.9*flowU[3]+0.6708203932499369*flowU[2]-0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[15] = (-0.5590169943749475*flowU[5])+0.4472135954999581*flowU[4]+0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[16] = 0.4472135954999581*(flowU[5]+flowU[4])-0.9*flowU[3]-0.6708203932499369*flowU[2]+0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[17] = 0.4472135954999581*(flowU[5]+flowU[4])+0.9*flowU[3]+0.6708203932499369*(flowU[2]+flowU[1])+0.5*flowU[0]; 

  vtSqOrd[0] = 0.5*vtSq[0]-0.5590169943749475*(vtSq[5]+vtSq[4]); 
  vtSqOrd[1] = 0.4472135954999581*vtSq[5]-0.5590169943749475*vtSq[4]-0.6708203932499369*vtSq[2]+0.5*vtSq[0]; 
  vtSqOrd[2] = 0.4472135954999581*vtSq[5]-0.5590169943749475*vtSq[4]+0.6708203932499369*vtSq[2]+0.5*vtSq[0]; 
  vtSqOrd[3] = (-0.5590169943749475*vtSq[5])+0.4472135954999581*vtSq[4]-0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[4] = 0.4472135954999581*(vtSq[5]+vtSq[4])+0.9*vtSq[3]-0.6708203932499369*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 
  vtSqOrd[5] = 0.4472135954999581*(vtSq[5]+vtSq[4])-0.9*vtSq[3]+0.6708203932499369*vtSq[2]-0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[6] = (-0.5590169943749475*vtSq[5])+0.4472135954999581*vtSq[4]+0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[7] = 0.4472135954999581*(vtSq[5]+vtSq[4])-0.9*vtSq[3]-0.6708203932499369*vtSq[2]+0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[8] = 0.4472135954999581*(vtSq[5]+vtSq[4])+0.9*vtSq[3]+0.6708203932499369*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (0.1591549430918953*m0Ord[0])/vtSqOrd[0]; 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (0.1591549430918953*m0Ord[1])/vtSqOrd[1]; 
  else
    fMFacOrd[1] = 0.0;
  if ((vtSqOrd[2] > 0.) && (m0Ord[2] > 0.))
    fMFacOrd[2] = (0.1591549430918953*m0Ord[2])/vtSqOrd[2]; 
  else
    fMFacOrd[2] = 0.0;
  if ((vtSqOrd[3] > 0.) && (m0Ord[3] > 0.))
    fMFacOrd[3] = (0.1591549430918953*m0Ord[3])/vtSqOrd[3]; 
  else
    fMFacOrd[3] = 0.0;
  if ((vtSqOrd[4] > 0.) && (m0Ord[4] > 0.))
    fMFacOrd[4] = (0.1591549430918953*m0Ord[4])/vtSqOrd[4]; 
  else
    fMFacOrd[4] = 0.0;
  if ((vtSqOrd[5] > 0.) && (m0Ord[5] > 0.))
    fMFacOrd[5] = (0.1591549430918953*m0Ord[5])/vtSqOrd[5]; 
  else
    fMFacOrd[5] = 0.0;
  if ((vtSqOrd[6] > 0.) && (m0Ord[6] > 0.))
    fMFacOrd[6] = (0.1591549430918953*m0Ord[6])/vtSqOrd[6]; 
  else
    fMFacOrd[6] = 0.0;
  if ((vtSqOrd[7] > 0.) && (m0Ord[7] > 0.))
    fMFacOrd[7] = (0.1591549430918953*m0Ord[7])/vtSqOrd[7]; 
  else
    fMFacOrd[7] = 0.0;
  if ((vtSqOrd[8] > 0.) && (m0Ord[8] > 0.))
    fMFacOrd[8] = (0.1591549430918953*m0Ord[8])/vtSqOrd[8]; 
  else
    fMFacOrd[8] = 0.0;

}

void MaxwellianOnBasisGauss2x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[81];
  if ((vtSqOrd[0] > 0.) && (fMFacOrd[0] > 0.)) {
    fMquad[0] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[9],2.0)+std::pow(wc[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[1] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[2] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[3] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[9],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[4] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[5] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[6] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[9],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[7] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[8] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  } else {
    fMquad[0] = 0.0;
    fMquad[1] = 0.0;
    fMquad[2] = 0.0;
    fMquad[3] = 0.0;
    fMquad[4] = 0.0;
    fMquad[5] = 0.0;
    fMquad[6] = 0.0;
    fMquad[7] = 0.0;
    fMquad[8] = 0.0;
  };
  if ((vtSqOrd[1] > 0.) && (fMFacOrd[1] > 0.)) {
    fMquad[9] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[10],2.0)+std::pow(wc[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[10] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[11] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[12] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[10],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[13] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[14] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[15] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[10],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[16] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[17] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  } else {
    fMquad[9] = 0.0;
    fMquad[10] = 0.0;
    fMquad[11] = 0.0;
    fMquad[12] = 0.0;
    fMquad[13] = 0.0;
    fMquad[14] = 0.0;
    fMquad[15] = 0.0;
    fMquad[16] = 0.0;
    fMquad[17] = 0.0;
  };
  if ((vtSqOrd[2] > 0.) && (fMFacOrd[2] > 0.)) {
    fMquad[18] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[11],2.0)+std::pow(wc[2]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
    fMquad[19] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
    fMquad[20] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
    fMquad[21] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[11],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[22] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[23] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[24] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[11],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[25] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
    fMquad[26] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  } else {
    fMquad[18] = 0.0;
    fMquad[19] = 0.0;
    fMquad[20] = 0.0;
    fMquad[21] = 0.0;
    fMquad[22] = 0.0;
    fMquad[23] = 0.0;
    fMquad[24] = 0.0;
    fMquad[25] = 0.0;
    fMquad[26] = 0.0;
  };
  if ((vtSqOrd[3] > 0.) && (fMFacOrd[3] > 0.)) {
    fMquad[27] = fMFacOrd[3]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[12],2.0)+std::pow(wc[2]-1.0*flowUOrd[3],2.0)))/vtSqOrd[3]); 
    fMquad[28] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[3],2.0)))/vtSqOrd[3]); 
    fMquad[29] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[3],2.0)))/vtSqOrd[3]); 
    fMquad[30] = fMFacOrd[3]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[12],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[31] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[32] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[33] = fMFacOrd[3]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[12],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[34] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
    fMquad[35] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  } else {
    fMquad[27] = 0.0;
    fMquad[28] = 0.0;
    fMquad[29] = 0.0;
    fMquad[30] = 0.0;
    fMquad[31] = 0.0;
    fMquad[32] = 0.0;
    fMquad[33] = 0.0;
    fMquad[34] = 0.0;
    fMquad[35] = 0.0;
  };
  if ((vtSqOrd[4] > 0.) && (fMFacOrd[4] > 0.)) {
    fMquad[36] = fMFacOrd[4]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[13],2.0)+std::pow(wc[2]-1.0*flowUOrd[4],2.0)))/vtSqOrd[4]); 
    fMquad[37] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[4],2.0)))/vtSqOrd[4]); 
    fMquad[38] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[4],2.0)))/vtSqOrd[4]); 
    fMquad[39] = fMFacOrd[4]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[13],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
    fMquad[40] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
    fMquad[41] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
    fMquad[42] = fMFacOrd[4]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[13],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
    fMquad[43] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
    fMquad[44] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  } else {
    fMquad[36] = 0.0;
    fMquad[37] = 0.0;
    fMquad[38] = 0.0;
    fMquad[39] = 0.0;
    fMquad[40] = 0.0;
    fMquad[41] = 0.0;
    fMquad[42] = 0.0;
    fMquad[43] = 0.0;
    fMquad[44] = 0.0;
  };
  if ((vtSqOrd[5] > 0.) && (fMFacOrd[5] > 0.)) {
    fMquad[45] = fMFacOrd[5]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[14],2.0)+std::pow(wc[2]-1.0*flowUOrd[5],2.0)))/vtSqOrd[5]); 
    fMquad[46] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[5],2.0)))/vtSqOrd[5]); 
    fMquad[47] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[5],2.0)))/vtSqOrd[5]); 
    fMquad[48] = fMFacOrd[5]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[14],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
    fMquad[49] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
    fMquad[50] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
    fMquad[51] = fMFacOrd[5]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[14],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
    fMquad[52] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
    fMquad[53] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  } else {
    fMquad[45] = 0.0;
    fMquad[46] = 0.0;
    fMquad[47] = 0.0;
    fMquad[48] = 0.0;
    fMquad[49] = 0.0;
    fMquad[50] = 0.0;
    fMquad[51] = 0.0;
    fMquad[52] = 0.0;
    fMquad[53] = 0.0;
  };
  if ((vtSqOrd[6] > 0.) && (fMFacOrd[6] > 0.)) {
    fMquad[54] = fMFacOrd[6]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[15],2.0)+std::pow(wc[2]-1.0*flowUOrd[6],2.0)))/vtSqOrd[6]); 
    fMquad[55] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[6],2.0)))/vtSqOrd[6]); 
    fMquad[56] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[6],2.0)))/vtSqOrd[6]); 
    fMquad[57] = fMFacOrd[6]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[15],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
    fMquad[58] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
    fMquad[59] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
    fMquad[60] = fMFacOrd[6]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[15],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
    fMquad[61] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
    fMquad[62] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  } else {
    fMquad[54] = 0.0;
    fMquad[55] = 0.0;
    fMquad[56] = 0.0;
    fMquad[57] = 0.0;
    fMquad[58] = 0.0;
    fMquad[59] = 0.0;
    fMquad[60] = 0.0;
    fMquad[61] = 0.0;
    fMquad[62] = 0.0;
  };
  if ((vtSqOrd[7] > 0.) && (fMFacOrd[7] > 0.)) {
    fMquad[63] = fMFacOrd[7]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[16],2.0)+std::pow(wc[2]-1.0*flowUOrd[7],2.0)))/vtSqOrd[7]); 
    fMquad[64] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[7],2.0)))/vtSqOrd[7]); 
    fMquad[65] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[7],2.0)))/vtSqOrd[7]); 
    fMquad[66] = fMFacOrd[7]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[16],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
    fMquad[67] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
    fMquad[68] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
    fMquad[69] = fMFacOrd[7]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[16],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
    fMquad[70] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
    fMquad[71] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  } else {
    fMquad[63] = 0.0;
    fMquad[64] = 0.0;
    fMquad[65] = 0.0;
    fMquad[66] = 0.0;
    fMquad[67] = 0.0;
    fMquad[68] = 0.0;
    fMquad[69] = 0.0;
    fMquad[70] = 0.0;
    fMquad[71] = 0.0;
  };
  if ((vtSqOrd[8] > 0.) && (fMFacOrd[8] > 0.)) {
    fMquad[72] = fMFacOrd[8]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[17],2.0)+std::pow(wc[2]-1.0*flowUOrd[8],2.0)))/vtSqOrd[8]); 
    fMquad[73] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[8],2.0)))/vtSqOrd[8]); 
    fMquad[74] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[8],2.0)))/vtSqOrd[8]); 
    fMquad[75] = fMFacOrd[8]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[17],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
    fMquad[76] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
    fMquad[77] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
    fMquad[78] = fMFacOrd[8]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[17],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
    fMquad[79] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
    fMquad[80] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  } else {
    fMquad[72] = 0.0;
    fMquad[73] = 0.0;
    fMquad[74] = 0.0;
    fMquad[75] = 0.0;
    fMquad[76] = 0.0;
    fMquad[77] = 0.0;
    fMquad[78] = 0.0;
    fMquad[79] = 0.0;
    fMquad[80] = 0.0;
  };

  fMOut[0] = 0.02381496723060508*(fMquad[80]+fMquad[79])+0.03810394756896813*fMquad[78]+0.02381496723060508*(fMquad[77]+fMquad[76])+0.03810394756896813*fMquad[75]+0.03810394756896812*(fMquad[74]+fMquad[73])+0.06096631611034901*fMquad[72]+0.02381496723060508*(fMquad[71]+fMquad[70])+0.03810394756896813*fMquad[69]+0.02381496723060508*(fMquad[68]+fMquad[67])+0.03810394756896813*fMquad[66]+0.03810394756896812*(fMquad[65]+fMquad[64])+0.06096631611034901*fMquad[63]+0.03810394756896812*(fMquad[62]+fMquad[61])+0.06096631611034901*fMquad[60]+0.03810394756896812*(fMquad[59]+fMquad[58])+0.06096631611034901*fMquad[57]+0.06096631611034901*(fMquad[56]+fMquad[55])+0.09754610577655842*fMquad[54]+0.02381496723060508*(fMquad[53]+fMquad[52])+0.03810394756896813*fMquad[51]+0.02381496723060508*(fMquad[50]+fMquad[49])+0.03810394756896813*fMquad[48]+0.03810394756896812*(fMquad[47]+fMquad[46])+0.06096631611034901*fMquad[45]+0.02381496723060508*(fMquad[44]+fMquad[43])+0.03810394756896813*fMquad[42]+0.02381496723060508*(fMquad[41]+fMquad[40])+0.03810394756896813*fMquad[39]+0.03810394756896812*(fMquad[38]+fMquad[37])+0.06096631611034901*fMquad[36]+0.03810394756896812*(fMquad[35]+fMquad[34])+0.06096631611034901*fMquad[33]+0.03810394756896812*(fMquad[32]+fMquad[31])+0.06096631611034901*fMquad[30]+0.06096631611034901*(fMquad[29]+fMquad[28])+0.09754610577655842*fMquad[27]+0.03810394756896812*(fMquad[26]+fMquad[25])+0.06096631611034901*fMquad[24]+0.03810394756896812*(fMquad[23]+fMquad[22])+0.06096631611034901*fMquad[21]+0.06096631611034901*(fMquad[20]+fMquad[19])+0.09754610577655842*fMquad[18]+0.03810394756896812*(fMquad[17]+fMquad[16])+0.06096631611034901*fMquad[15]+0.03810394756896812*(fMquad[14]+fMquad[13])+0.06096631611034901*fMquad[12]+0.06096631611034901*(fMquad[11]+fMquad[10])+0.09754610577655842*fMquad[9]+0.06096631611034901*(fMquad[8]+fMquad[7])+0.09754610577655842*fMquad[6]+0.06096631611034901*(fMquad[5]+fMquad[4])+0.09754610577655842*(fMquad[3]+fMquad[2]+fMquad[1])+0.1560737692424935*fMquad[0]; 
  fMOut[1] = 0.03195113136573771*(fMquad[80]+fMquad[79])+0.05112181018518035*fMquad[78]+0.03195113136573771*(fMquad[77]+fMquad[76])+0.05112181018518035*fMquad[75]+0.05112181018518035*(fMquad[74]+fMquad[73])+0.08179489629628857*fMquad[72]+0.03195113136573771*(fMquad[71]+fMquad[70])+0.05112181018518035*fMquad[69]+0.03195113136573771*(fMquad[68]+fMquad[67])+0.05112181018518035*fMquad[66]+0.05112181018518035*(fMquad[65]+fMquad[64])+0.08179489629628857*fMquad[63]+0.05112181018518035*(fMquad[62]+fMquad[61])+0.08179489629628857*fMquad[60]+0.05112181018518035*(fMquad[59]+fMquad[58])+0.08179489629628857*(fMquad[57]+fMquad[56]+fMquad[55])+0.1308718340740617*fMquad[54]-0.03195113136573771*(fMquad[53]+fMquad[52])-0.05112181018518035*fMquad[51]-0.03195113136573771*(fMquad[50]+fMquad[49])-0.05112181018518035*fMquad[48]-0.05112181018518035*(fMquad[47]+fMquad[46])-0.08179489629628857*fMquad[45]-0.03195113136573771*(fMquad[44]+fMquad[43])-0.05112181018518035*fMquad[42]-0.03195113136573771*(fMquad[41]+fMquad[40])-0.05112181018518035*fMquad[39]-0.05112181018518035*(fMquad[38]+fMquad[37])-0.08179489629628857*fMquad[36]-0.05112181018518035*(fMquad[35]+fMquad[34])-0.08179489629628857*fMquad[33]-0.05112181018518035*(fMquad[32]+fMquad[31])-0.08179489629628857*(fMquad[30]+fMquad[29]+fMquad[28])-0.1308718340740617*fMquad[27]; 
  fMOut[2] = 0.03195113136573771*(fMquad[80]+fMquad[79])+0.05112181018518035*fMquad[78]+0.03195113136573771*(fMquad[77]+fMquad[76])+0.05112181018518035*fMquad[75]+0.05112181018518035*(fMquad[74]+fMquad[73])+0.08179489629628857*fMquad[72]-0.03195113136573771*(fMquad[71]+fMquad[70])-0.05112181018518035*fMquad[69]-0.03195113136573771*(fMquad[68]+fMquad[67])-0.05112181018518035*fMquad[66]-0.05112181018518035*(fMquad[65]+fMquad[64])-0.08179489629628857*fMquad[63]+0.03195113136573771*(fMquad[53]+fMquad[52])+0.05112181018518035*fMquad[51]+0.03195113136573771*(fMquad[50]+fMquad[49])+0.05112181018518035*fMquad[48]+0.05112181018518035*(fMquad[47]+fMquad[46])+0.08179489629628857*fMquad[45]-0.03195113136573771*(fMquad[44]+fMquad[43])-0.05112181018518035*fMquad[42]-0.03195113136573771*(fMquad[41]+fMquad[40])-0.05112181018518035*fMquad[39]-0.05112181018518035*(fMquad[38]+fMquad[37])-0.08179489629628857*fMquad[36]+0.05112181018518035*(fMquad[26]+fMquad[25])+0.08179489629628857*fMquad[24]+0.05112181018518035*(fMquad[23]+fMquad[22])+0.08179489629628857*(fMquad[21]+fMquad[20]+fMquad[19])+0.1308718340740617*fMquad[18]-0.05112181018518035*(fMquad[17]+fMquad[16])-0.08179489629628857*fMquad[15]-0.05112181018518035*(fMquad[14]+fMquad[13])-0.08179489629628857*(fMquad[12]+fMquad[11]+fMquad[10])-0.1308718340740617*fMquad[9]; 
  fMOut[3] = 0.03195113136573771*(fMquad[80]+fMquad[79])+0.05112181018518035*fMquad[78]-0.03195113136573771*(fMquad[77]+fMquad[76])-0.05112181018518035*fMquad[75]+0.03195113136573771*(fMquad[71]+fMquad[70])+0.05112181018518035*fMquad[69]-0.03195113136573771*(fMquad[68]+fMquad[67])-0.05112181018518035*fMquad[66]+0.05112181018518035*(fMquad[62]+fMquad[61])+0.08179489629628857*fMquad[60]-0.05112181018518035*(fMquad[59]+fMquad[58])-0.08179489629628857*fMquad[57]+0.03195113136573771*(fMquad[53]+fMquad[52])+0.05112181018518035*fMquad[51]-0.03195113136573771*(fMquad[50]+fMquad[49])-0.05112181018518035*fMquad[48]+0.03195113136573771*(fMquad[44]+fMquad[43])+0.05112181018518035*fMquad[42]-0.03195113136573771*(fMquad[41]+fMquad[40])-0.05112181018518035*fMquad[39]+0.05112181018518035*(fMquad[35]+fMquad[34])+0.08179489629628857*fMquad[33]-0.05112181018518035*(fMquad[32]+fMquad[31])-0.08179489629628857*fMquad[30]+0.05112181018518035*(fMquad[26]+fMquad[25])+0.08179489629628857*fMquad[24]-0.05112181018518035*(fMquad[23]+fMquad[22])-0.08179489629628857*fMquad[21]+0.05112181018518035*(fMquad[17]+fMquad[16])+0.08179489629628857*fMquad[15]-0.05112181018518035*(fMquad[14]+fMquad[13])-0.08179489629628857*fMquad[12]+0.08179489629628857*(fMquad[8]+fMquad[7])+0.1308718340740617*fMquad[6]-0.08179489629628857*(fMquad[5]+fMquad[4])-0.1308718340740617*fMquad[3]; 
  fMOut[4] = 0.03195113136573771*fMquad[80]-0.03195113136573771*fMquad[79]+0.03195113136573771*fMquad[77]-0.03195113136573771*fMquad[76]+0.05112181018518035*fMquad[74]-0.05112181018518035*fMquad[73]+0.03195113136573771*fMquad[71]-0.03195113136573771*fMquad[70]+0.03195113136573771*fMquad[68]-0.03195113136573771*fMquad[67]+0.05112181018518035*fMquad[65]-0.05112181018518035*fMquad[64]+0.05112181018518035*fMquad[62]-0.05112181018518035*fMquad[61]+0.05112181018518035*fMquad[59]-0.05112181018518035*fMquad[58]+0.08179489629628857*fMquad[56]-0.08179489629628857*fMquad[55]+0.03195113136573771*fMquad[53]-0.03195113136573771*fMquad[52]+0.03195113136573771*fMquad[50]-0.03195113136573771*fMquad[49]+0.05112181018518035*fMquad[47]-0.05112181018518035*fMquad[46]+0.03195113136573771*fMquad[44]-0.03195113136573771*fMquad[43]+0.03195113136573771*fMquad[41]-0.03195113136573771*fMquad[40]+0.05112181018518035*fMquad[38]-0.05112181018518035*fMquad[37]+0.05112181018518035*fMquad[35]-0.05112181018518035*fMquad[34]+0.05112181018518035*fMquad[32]-0.05112181018518035*fMquad[31]+0.08179489629628857*fMquad[29]-0.08179489629628857*fMquad[28]+0.05112181018518035*fMquad[26]-0.05112181018518035*fMquad[25]+0.05112181018518035*fMquad[23]-0.05112181018518035*fMquad[22]+0.08179489629628857*fMquad[20]-0.08179489629628857*fMquad[19]+0.05112181018518035*fMquad[17]-0.05112181018518035*fMquad[16]+0.05112181018518035*fMquad[14]-0.05112181018518035*fMquad[13]+0.08179489629628857*fMquad[11]-0.08179489629628857*fMquad[10]+0.08179489629628857*fMquad[8]-0.08179489629628857*fMquad[7]+0.08179489629628857*fMquad[5]-0.08179489629628857*fMquad[4]+0.1308718340740617*fMquad[2]-0.1308718340740617*fMquad[1]; 
  fMOut[5] = 0.04286694101508914*(fMquad[80]+fMquad[79])+0.06858710562414264*fMquad[78]+0.04286694101508914*(fMquad[77]+fMquad[76])+0.06858710562414264*fMquad[75]+0.06858710562414262*(fMquad[74]+fMquad[73])+0.1097393689986282*fMquad[72]-0.04286694101508914*(fMquad[71]+fMquad[70])-0.06858710562414264*fMquad[69]-0.04286694101508914*(fMquad[68]+fMquad[67])-0.06858710562414264*fMquad[66]-0.06858710562414262*(fMquad[65]+fMquad[64])-0.1097393689986282*fMquad[63]-0.04286694101508914*(fMquad[53]+fMquad[52])-0.06858710562414264*fMquad[51]-0.04286694101508914*(fMquad[50]+fMquad[49])-0.06858710562414264*fMquad[48]-0.06858710562414262*(fMquad[47]+fMquad[46])-0.1097393689986282*fMquad[45]+0.04286694101508914*(fMquad[44]+fMquad[43])+0.06858710562414264*fMquad[42]+0.04286694101508914*(fMquad[41]+fMquad[40])+0.06858710562414264*fMquad[39]+0.06858710562414262*(fMquad[38]+fMquad[37])+0.1097393689986282*fMquad[36]; 
  fMOut[6] = 0.04286694101508914*(fMquad[80]+fMquad[79])+0.06858710562414264*fMquad[78]-0.04286694101508914*(fMquad[77]+fMquad[76])-0.06858710562414264*fMquad[75]+0.04286694101508914*(fMquad[71]+fMquad[70])+0.06858710562414264*fMquad[69]-0.04286694101508914*(fMquad[68]+fMquad[67])-0.06858710562414264*fMquad[66]+0.06858710562414262*(fMquad[62]+fMquad[61])+0.1097393689986282*fMquad[60]-0.06858710562414262*(fMquad[59]+fMquad[58])-0.1097393689986282*fMquad[57]-0.04286694101508914*(fMquad[53]+fMquad[52])-0.06858710562414264*fMquad[51]+0.04286694101508914*(fMquad[50]+fMquad[49])+0.06858710562414264*fMquad[48]-0.04286694101508914*(fMquad[44]+fMquad[43])-0.06858710562414264*fMquad[42]+0.04286694101508914*(fMquad[41]+fMquad[40])+0.06858710562414264*fMquad[39]-0.06858710562414262*(fMquad[35]+fMquad[34])-0.1097393689986282*fMquad[33]+0.06858710562414262*(fMquad[32]+fMquad[31])+0.1097393689986282*fMquad[30]; 
  fMOut[7] = 0.04286694101508914*(fMquad[80]+fMquad[79])+0.06858710562414264*fMquad[78]-0.04286694101508914*(fMquad[77]+fMquad[76])-0.06858710562414264*fMquad[75]-0.04286694101508914*(fMquad[71]+fMquad[70])-0.06858710562414264*fMquad[69]+0.04286694101508914*(fMquad[68]+fMquad[67])+0.06858710562414264*fMquad[66]+0.04286694101508914*(fMquad[53]+fMquad[52])+0.06858710562414264*fMquad[51]-0.04286694101508914*(fMquad[50]+fMquad[49])-0.06858710562414264*fMquad[48]-0.04286694101508914*(fMquad[44]+fMquad[43])-0.06858710562414264*fMquad[42]+0.04286694101508914*(fMquad[41]+fMquad[40])+0.06858710562414264*fMquad[39]+0.06858710562414262*(fMquad[26]+fMquad[25])+0.1097393689986282*fMquad[24]-0.06858710562414262*(fMquad[23]+fMquad[22])-0.1097393689986282*fMquad[21]-0.06858710562414262*(fMquad[17]+fMquad[16])-0.1097393689986282*fMquad[15]+0.06858710562414262*(fMquad[14]+fMquad[13])+0.1097393689986282*fMquad[12]; 
  fMOut[8] = 0.04286694101508914*fMquad[80]-0.04286694101508914*fMquad[79]+0.04286694101508914*fMquad[77]-0.04286694101508914*fMquad[76]+0.06858710562414262*fMquad[74]-0.06858710562414262*fMquad[73]+0.04286694101508914*fMquad[71]-0.04286694101508914*fMquad[70]+0.04286694101508914*fMquad[68]-0.04286694101508914*fMquad[67]+0.06858710562414262*fMquad[65]-0.06858710562414262*fMquad[64]+0.06858710562414262*fMquad[62]-0.06858710562414262*fMquad[61]+0.06858710562414262*fMquad[59]-0.06858710562414262*fMquad[58]+0.1097393689986282*fMquad[56]-0.1097393689986282*fMquad[55]-0.04286694101508914*fMquad[53]+0.04286694101508914*fMquad[52]-0.04286694101508914*fMquad[50]+0.04286694101508914*fMquad[49]-0.06858710562414262*fMquad[47]+0.06858710562414262*fMquad[46]-0.04286694101508914*fMquad[44]+0.04286694101508914*fMquad[43]-0.04286694101508914*fMquad[41]+0.04286694101508914*fMquad[40]-0.06858710562414262*fMquad[38]+0.06858710562414262*fMquad[37]-0.06858710562414262*fMquad[35]+0.06858710562414262*fMquad[34]-0.06858710562414262*fMquad[32]+0.06858710562414262*fMquad[31]-0.1097393689986282*fMquad[29]+0.1097393689986282*fMquad[28]; 
  fMOut[9] = 0.04286694101508914*fMquad[80]-0.04286694101508914*fMquad[79]+0.04286694101508914*fMquad[77]-0.04286694101508914*fMquad[76]+0.06858710562414262*fMquad[74]-0.06858710562414262*fMquad[73]-0.04286694101508914*fMquad[71]+0.04286694101508914*fMquad[70]-0.04286694101508914*fMquad[68]+0.04286694101508914*fMquad[67]-0.06858710562414262*fMquad[65]+0.06858710562414262*fMquad[64]+0.04286694101508914*fMquad[53]-0.04286694101508914*fMquad[52]+0.04286694101508914*fMquad[50]-0.04286694101508914*fMquad[49]+0.06858710562414262*fMquad[47]-0.06858710562414262*fMquad[46]-0.04286694101508914*fMquad[44]+0.04286694101508914*fMquad[43]-0.04286694101508914*fMquad[41]+0.04286694101508914*fMquad[40]-0.06858710562414262*fMquad[38]+0.06858710562414262*(fMquad[37]+fMquad[26])-0.06858710562414262*fMquad[25]+0.06858710562414262*fMquad[23]-0.06858710562414262*fMquad[22]+0.1097393689986282*fMquad[20]-0.1097393689986282*fMquad[19]-0.06858710562414262*fMquad[17]+0.06858710562414262*fMquad[16]-0.06858710562414262*fMquad[14]+0.06858710562414262*fMquad[13]-0.1097393689986282*fMquad[11]+0.1097393689986282*fMquad[10]; 
  fMOut[10] = 0.04286694101508914*fMquad[80]-0.04286694101508914*(fMquad[79]+fMquad[77])+0.04286694101508914*(fMquad[76]+fMquad[71])-0.04286694101508914*(fMquad[70]+fMquad[68])+0.04286694101508914*fMquad[67]+0.06858710562414262*fMquad[62]-0.06858710562414262*(fMquad[61]+fMquad[59])+0.06858710562414262*fMquad[58]+0.04286694101508914*fMquad[53]-0.04286694101508914*(fMquad[52]+fMquad[50])+0.04286694101508914*(fMquad[49]+fMquad[44])-0.04286694101508914*(fMquad[43]+fMquad[41])+0.04286694101508914*fMquad[40]+0.06858710562414262*fMquad[35]-0.06858710562414262*(fMquad[34]+fMquad[32])+0.06858710562414262*(fMquad[31]+fMquad[26])-0.06858710562414262*(fMquad[25]+fMquad[23])+0.06858710562414262*(fMquad[22]+fMquad[17])-0.06858710562414262*(fMquad[16]+fMquad[14])+0.06858710562414262*fMquad[13]+0.1097393689986282*fMquad[8]-0.1097393689986282*(fMquad[7]+fMquad[5])+0.1097393689986282*fMquad[4]; 
  fMOut[11] = 0.02130075424382515*(fMquad[80]+fMquad[79])+0.03408120679012025*fMquad[78]+0.02130075424382515*(fMquad[77]+fMquad[76])+0.03408120679012025*fMquad[75]+0.03408120679012024*(fMquad[74]+fMquad[73])+0.0545299308641924*fMquad[72]+0.02130075424382515*(fMquad[71]+fMquad[70])+0.03408120679012025*fMquad[69]+0.02130075424382515*(fMquad[68]+fMquad[67])+0.03408120679012025*fMquad[66]+0.03408120679012024*(fMquad[65]+fMquad[64])+0.0545299308641924*fMquad[63]+0.03408120679012024*(fMquad[62]+fMquad[61])+0.0545299308641924*fMquad[60]+0.03408120679012024*(fMquad[59]+fMquad[58])+0.0545299308641924*(fMquad[57]+fMquad[56]+fMquad[55])+0.08724788938270785*fMquad[54]+0.02130075424382515*(fMquad[53]+fMquad[52])+0.03408120679012025*fMquad[51]+0.02130075424382515*(fMquad[50]+fMquad[49])+0.03408120679012025*fMquad[48]+0.03408120679012024*(fMquad[47]+fMquad[46])+0.0545299308641924*fMquad[45]+0.02130075424382515*(fMquad[44]+fMquad[43])+0.03408120679012025*fMquad[42]+0.02130075424382515*(fMquad[41]+fMquad[40])+0.03408120679012025*fMquad[39]+0.03408120679012024*(fMquad[38]+fMquad[37])+0.0545299308641924*fMquad[36]+0.03408120679012024*(fMquad[35]+fMquad[34])+0.0545299308641924*fMquad[33]+0.03408120679012024*(fMquad[32]+fMquad[31])+0.0545299308641924*(fMquad[30]+fMquad[29]+fMquad[28])+0.08724788938270785*fMquad[27]-0.04260150848765029*(fMquad[26]+fMquad[25])-0.06816241358024047*fMquad[24]-0.04260150848765029*(fMquad[23]+fMquad[22])-0.06816241358024047*fMquad[21]-0.06816241358024049*(fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]-0.04260150848765029*(fMquad[17]+fMquad[16])-0.06816241358024047*fMquad[15]-0.04260150848765029*(fMquad[14]+fMquad[13])-0.06816241358024047*fMquad[12]-0.06816241358024049*(fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]-0.06816241358024049*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]-0.06816241358024049*(fMquad[5]+fMquad[4])-0.1090598617283848*(fMquad[3]+fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[12] = 0.02130075424382515*(fMquad[80]+fMquad[79])+0.03408120679012025*fMquad[78]+0.02130075424382515*(fMquad[77]+fMquad[76])+0.03408120679012025*fMquad[75]+0.03408120679012024*(fMquad[74]+fMquad[73])+0.0545299308641924*fMquad[72]+0.02130075424382515*(fMquad[71]+fMquad[70])+0.03408120679012025*fMquad[69]+0.02130075424382515*(fMquad[68]+fMquad[67])+0.03408120679012025*fMquad[66]+0.03408120679012024*(fMquad[65]+fMquad[64])+0.0545299308641924*fMquad[63]-0.04260150848765029*(fMquad[62]+fMquad[61])-0.06816241358024047*fMquad[60]-0.04260150848765029*(fMquad[59]+fMquad[58])-0.06816241358024047*fMquad[57]-0.06816241358024049*(fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382515*(fMquad[53]+fMquad[52])+0.03408120679012025*fMquad[51]+0.02130075424382515*(fMquad[50]+fMquad[49])+0.03408120679012025*fMquad[48]+0.03408120679012024*(fMquad[47]+fMquad[46])+0.0545299308641924*fMquad[45]+0.02130075424382515*(fMquad[44]+fMquad[43])+0.03408120679012025*fMquad[42]+0.02130075424382515*(fMquad[41]+fMquad[40])+0.03408120679012025*fMquad[39]+0.03408120679012024*(fMquad[38]+fMquad[37])+0.0545299308641924*fMquad[36]-0.04260150848765029*(fMquad[35]+fMquad[34])-0.06816241358024047*fMquad[33]-0.04260150848765029*(fMquad[32]+fMquad[31])-0.06816241358024047*fMquad[30]-0.06816241358024049*(fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012024*(fMquad[26]+fMquad[25])+0.0545299308641924*fMquad[24]+0.03408120679012024*(fMquad[23]+fMquad[22])+0.0545299308641924*(fMquad[21]+fMquad[20]+fMquad[19])+0.08724788938270785*fMquad[18]+0.03408120679012024*(fMquad[17]+fMquad[16])+0.0545299308641924*fMquad[15]+0.03408120679012024*(fMquad[14]+fMquad[13])+0.0545299308641924*(fMquad[12]+fMquad[11]+fMquad[10])+0.08724788938270785*fMquad[9]-0.06816241358024049*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]-0.06816241358024049*(fMquad[5]+fMquad[4])-0.1090598617283848*(fMquad[3]+fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[13] = 0.02130075424382515*(fMquad[80]+fMquad[79])+0.03408120679012025*fMquad[78]+0.02130075424382515*(fMquad[77]+fMquad[76])+0.03408120679012025*fMquad[75]-0.04260150848765029*(fMquad[74]+fMquad[73])-0.06816241358024047*fMquad[72]+0.02130075424382515*(fMquad[71]+fMquad[70])+0.03408120679012025*fMquad[69]+0.02130075424382515*(fMquad[68]+fMquad[67])+0.03408120679012025*fMquad[66]-0.04260150848765029*(fMquad[65]+fMquad[64])-0.06816241358024047*fMquad[63]+0.03408120679012024*(fMquad[62]+fMquad[61])+0.0545299308641924*fMquad[60]+0.03408120679012024*(fMquad[59]+fMquad[58])+0.0545299308641924*fMquad[57]-0.06816241358024049*(fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382515*(fMquad[53]+fMquad[52])+0.03408120679012025*fMquad[51]+0.02130075424382515*(fMquad[50]+fMquad[49])+0.03408120679012025*fMquad[48]-0.04260150848765029*(fMquad[47]+fMquad[46])-0.06816241358024047*fMquad[45]+0.02130075424382515*(fMquad[44]+fMquad[43])+0.03408120679012025*fMquad[42]+0.02130075424382515*(fMquad[41]+fMquad[40])+0.03408120679012025*fMquad[39]-0.04260150848765029*(fMquad[38]+fMquad[37])-0.06816241358024047*fMquad[36]+0.03408120679012024*(fMquad[35]+fMquad[34])+0.0545299308641924*fMquad[33]+0.03408120679012024*(fMquad[32]+fMquad[31])+0.0545299308641924*fMquad[30]-0.06816241358024049*(fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012024*(fMquad[26]+fMquad[25])+0.0545299308641924*fMquad[24]+0.03408120679012024*(fMquad[23]+fMquad[22])+0.0545299308641924*fMquad[21]-0.06816241358024049*(fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]+0.03408120679012024*(fMquad[17]+fMquad[16])+0.0545299308641924*fMquad[15]+0.03408120679012024*(fMquad[14]+fMquad[13])+0.0545299308641924*fMquad[12]-0.06816241358024049*(fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]+0.0545299308641924*(fMquad[8]+fMquad[7])+0.08724788938270785*fMquad[6]+0.0545299308641924*(fMquad[5]+fMquad[4])+0.08724788938270785*fMquad[3]-0.1090598617283848*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[14] = 0.02130075424382515*(fMquad[80]+fMquad[79])-0.0426015084876503*fMquad[78]+0.02130075424382515*(fMquad[77]+fMquad[76])-0.0426015084876503*fMquad[75]+0.03408120679012024*(fMquad[74]+fMquad[73])-0.06816241358024047*fMquad[72]+0.02130075424382515*(fMquad[71]+fMquad[70])-0.0426015084876503*fMquad[69]+0.02130075424382515*(fMquad[68]+fMquad[67])-0.0426015084876503*fMquad[66]+0.03408120679012024*(fMquad[65]+fMquad[64])-0.06816241358024047*fMquad[63]+0.03408120679012024*(fMquad[62]+fMquad[61])-0.06816241358024047*fMquad[60]+0.03408120679012024*(fMquad[59]+fMquad[58])-0.06816241358024047*fMquad[57]+0.0545299308641924*(fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382515*(fMquad[53]+fMquad[52])-0.0426015084876503*fMquad[51]+0.02130075424382515*(fMquad[50]+fMquad[49])-0.0426015084876503*fMquad[48]+0.03408120679012024*(fMquad[47]+fMquad[46])-0.06816241358024047*fMquad[45]+0.02130075424382515*(fMquad[44]+fMquad[43])-0.0426015084876503*fMquad[42]+0.02130075424382515*(fMquad[41]+fMquad[40])-0.0426015084876503*fMquad[39]+0.03408120679012024*(fMquad[38]+fMquad[37])-0.06816241358024047*fMquad[36]+0.03408120679012024*(fMquad[35]+fMquad[34])-0.06816241358024047*fMquad[33]+0.03408120679012024*(fMquad[32]+fMquad[31])-0.06816241358024047*fMquad[30]+0.0545299308641924*(fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012024*(fMquad[26]+fMquad[25])-0.06816241358024047*fMquad[24]+0.03408120679012024*(fMquad[23]+fMquad[22])-0.06816241358024047*fMquad[21]+0.0545299308641924*(fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]+0.03408120679012024*(fMquad[17]+fMquad[16])-0.06816241358024047*fMquad[15]+0.03408120679012024*(fMquad[14]+fMquad[13])-0.06816241358024047*fMquad[12]+0.0545299308641924*(fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]+0.0545299308641924*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]+0.0545299308641924*(fMquad[5]+fMquad[4])-0.1090598617283848*fMquad[3]+0.08724788938270785*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 

}
void GkMaxwellianOnBasisGauss2x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[9];
  m0Ord[0] = 0.5*den[0]-0.5590169943749475*(den[5]+den[4]); 
  m0Ord[1] = 0.4472135954999581*den[5]-0.5590169943749475*den[4]-0.6708203932499369*den[2]+0.5*den[0]; 
  m0Ord[2] = 0.4472135954999581*den[5]-0.5590169943749475*den[4]+0.6708203932499369*den[2]+0.5*den[0]; 
  m0Ord[3] = (-0.5590169943749475*den[5])+0.4472135954999581*den[4]-0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[4] = 0.4472135954999581*(den[5]+den[4])+0.9*den[3]-0.6708203932499369*(den[2]+den[1])+0.5*den[0]; 
  m0Ord[5] = 0.4472135954999581*(den[5]+den[4])-0.9*den[3]+0.6708203932499369*den[2]-0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[6] = (-0.5590169943749475*den[5])+0.4472135954999581*den[4]+0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[7] = 0.4472135954999581*(den[5]+den[4])-0.9*den[3]-0.6708203932499369*den[2]+0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[8] = 0.4472135954999581*(den[5]+den[4])+0.9*den[3]+0.6708203932499369*(den[2]+den[1])+0.5*den[0]; 

  flowUOrd[0] = 0.5*flowU[0]-0.5590169943749475*(flowU[5]+flowU[4]); 
  flowUOrd[1] = 0.4472135954999581*flowU[5]-0.5590169943749475*flowU[4]-0.6708203932499369*flowU[2]+0.5*flowU[0]; 
  flowUOrd[2] = 0.4472135954999581*flowU[5]-0.5590169943749475*flowU[4]+0.6708203932499369*flowU[2]+0.5*flowU[0]; 
  flowUOrd[3] = (-0.5590169943749475*flowU[5])+0.4472135954999581*flowU[4]-0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[4] = 0.4472135954999581*(flowU[5]+flowU[4])+0.9*flowU[3]-0.6708203932499369*(flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[5] = 0.4472135954999581*(flowU[5]+flowU[4])-0.9*flowU[3]+0.6708203932499369*flowU[2]-0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[6] = (-0.5590169943749475*flowU[5])+0.4472135954999581*flowU[4]+0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[7] = 0.4472135954999581*(flowU[5]+flowU[4])-0.9*flowU[3]-0.6708203932499369*flowU[2]+0.6708203932499369*flowU[1]+0.5*flowU[0]; 
  flowUOrd[8] = 0.4472135954999581*(flowU[5]+flowU[4])+0.9*flowU[3]+0.6708203932499369*(flowU[2]+flowU[1])+0.5*flowU[0]; 

  vtSqOrd[0] = 0.5*vtSq[0]-0.5590169943749475*(vtSq[5]+vtSq[4]); 
  vtSqOrd[1] = 0.4472135954999581*vtSq[5]-0.5590169943749475*vtSq[4]-0.6708203932499369*vtSq[2]+0.5*vtSq[0]; 
  vtSqOrd[2] = 0.4472135954999581*vtSq[5]-0.5590169943749475*vtSq[4]+0.6708203932499369*vtSq[2]+0.5*vtSq[0]; 
  vtSqOrd[3] = (-0.5590169943749475*vtSq[5])+0.4472135954999581*vtSq[4]-0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[4] = 0.4472135954999581*(vtSq[5]+vtSq[4])+0.9*vtSq[3]-0.6708203932499369*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 
  vtSqOrd[5] = 0.4472135954999581*(vtSq[5]+vtSq[4])-0.9*vtSq[3]+0.6708203932499369*vtSq[2]-0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[6] = (-0.5590169943749475*vtSq[5])+0.4472135954999581*vtSq[4]+0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[7] = 0.4472135954999581*(vtSq[5]+vtSq[4])-0.9*vtSq[3]-0.6708203932499369*vtSq[2]+0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[8] = 0.4472135954999581*(vtSq[5]+vtSq[4])+0.9*vtSq[3]+0.6708203932499369*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 

  bmagOrd[0] = 0.5*bmag[0]-0.5590169943749475*(bmag[5]+bmag[4]); 
  bmagOrd[1] = 0.4472135954999581*bmag[5]-0.5590169943749475*bmag[4]-0.6708203932499369*bmag[2]+0.5*bmag[0]; 
  bmagOrd[2] = 0.4472135954999581*bmag[5]-0.5590169943749475*bmag[4]+0.6708203932499369*bmag[2]+0.5*bmag[0]; 
  bmagOrd[3] = (-0.5590169943749475*bmag[5])+0.4472135954999581*bmag[4]-0.6708203932499369*bmag[1]+0.5*bmag[0]; 
  bmagOrd[4] = 0.4472135954999581*(bmag[5]+bmag[4])+0.9*bmag[3]-0.6708203932499369*(bmag[2]+bmag[1])+0.5*bmag[0]; 
  bmagOrd[5] = 0.4472135954999581*(bmag[5]+bmag[4])-0.9*bmag[3]+0.6708203932499369*bmag[2]-0.6708203932499369*bmag[1]+0.5*bmag[0]; 
  bmagOrd[6] = (-0.5590169943749475*bmag[5])+0.4472135954999581*bmag[4]+0.6708203932499369*bmag[1]+0.5*bmag[0]; 
  bmagOrd[7] = 0.4472135954999581*(bmag[5]+bmag[4])-0.9*bmag[3]-0.6708203932499369*bmag[2]+0.6708203932499369*bmag[1]+0.5*bmag[0]; 
  bmagOrd[8] = 0.4472135954999581*(bmag[5]+bmag[4])+0.9*bmag[3]+0.6708203932499369*(bmag[2]+bmag[1])+0.5*bmag[0]; 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (bmagOrd[0]*m0Ord[0])/std::pow(2.506628274631001*sqrt(vtSqOrd[0]),3.0); 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (bmagOrd[1]*m0Ord[1])/std::pow(2.506628274631001*sqrt(vtSqOrd[1]),3.0); 
  else
    fMFacOrd[1] = 0.0;
  if ((vtSqOrd[2] > 0.) && (m0Ord[2] > 0.))
    fMFacOrd[2] = (bmagOrd[2]*m0Ord[2])/std::pow(2.506628274631001*sqrt(vtSqOrd[2]),3.0); 
  else
    fMFacOrd[2] = 0.0;
  if ((vtSqOrd[3] > 0.) && (m0Ord[3] > 0.))
    fMFacOrd[3] = (bmagOrd[3]*m0Ord[3])/std::pow(2.506628274631001*sqrt(vtSqOrd[3]),3.0); 
  else
    fMFacOrd[3] = 0.0;
  if ((vtSqOrd[4] > 0.) && (m0Ord[4] > 0.))
    fMFacOrd[4] = (bmagOrd[4]*m0Ord[4])/std::pow(2.506628274631001*sqrt(vtSqOrd[4]),3.0); 
  else
    fMFacOrd[4] = 0.0;
  if ((vtSqOrd[5] > 0.) && (m0Ord[5] > 0.))
    fMFacOrd[5] = (bmagOrd[5]*m0Ord[5])/std::pow(2.506628274631001*sqrt(vtSqOrd[5]),3.0); 
  else
    fMFacOrd[5] = 0.0;
  if ((vtSqOrd[6] > 0.) && (m0Ord[6] > 0.))
    fMFacOrd[6] = (bmagOrd[6]*m0Ord[6])/std::pow(2.506628274631001*sqrt(vtSqOrd[6]),3.0); 
  else
    fMFacOrd[6] = 0.0;
  if ((vtSqOrd[7] > 0.) && (m0Ord[7] > 0.))
    fMFacOrd[7] = (bmagOrd[7]*m0Ord[7])/std::pow(2.506628274631001*sqrt(vtSqOrd[7]),3.0); 
  else
    fMFacOrd[7] = 0.0;
  if ((vtSqOrd[8] > 0.) && (m0Ord[8] > 0.))
    fMFacOrd[8] = (bmagOrd[8]*m0Ord[8])/std::pow(2.506628274631001*sqrt(vtSqOrd[8]),3.0); 
  else
    fMFacOrd[8] = 0.0;

}

void GkMaxwellianOnBasisGauss2x2vMaxUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[9];
  m0Ord[0] = 0.5*den[0]-0.5590169943749475*(den[5]+den[4]); 
  m0Ord[1] = 0.4472135954999581*den[5]-0.5590169943749475*den[4]-0.6708203932499369*den[2]+0.5*den[0]; 
  m0Ord[2] = 0.4472135954999581*den[5]-0.5590169943749475*den[4]+0.6708203932499369*den[2]+0.5*den[0]; 
  m0Ord[3] = (-0.5590169943749475*den[5])+0.4472135954999581*den[4]-0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[4] = 0.4472135954999581*(den[5]+den[4])+0.9*den[3]-0.6708203932499369*(den[2]+den[1])+0.5*den[0]; 
  m0Ord[5] = 0.4472135954999581*(den[5]+den[4])-0.9*den[3]+0.6708203932499369*den[2]-0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[6] = (-0.5590169943749475*den[5])+0.4472135954999581*den[4]+0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[7] = 0.4472135954999581*(den[5]+den[4])-0.9*den[3]-0.6708203932499369*den[2]+0.6708203932499369*den[1]+0.5*den[0]; 
  m0Ord[8] = 0.4472135954999581*(den[5]+den[4])+0.9*den[3]+0.6708203932499369*(den[2]+den[1])+0.5*den[0]; 

  flowUOrd[0] = 0.5*flowU[12]-0.5590169943749475*(flowU[17]+flowU[16]); 
  flowUOrd[1] = 0.4472135954999581*flowU[17]-0.5590169943749475*flowU[16]-0.6708203932499369*flowU[14]+0.5*flowU[12]; 
  flowUOrd[2] = 0.4472135954999581*flowU[17]-0.5590169943749475*flowU[16]+0.6708203932499369*flowU[14]+0.5*flowU[12]; 
  flowUOrd[3] = (-0.5590169943749475*flowU[17])+0.4472135954999581*flowU[16]-0.6708203932499369*flowU[13]+0.5*flowU[12]; 
  flowUOrd[4] = 0.4472135954999581*(flowU[17]+flowU[16])+0.9*flowU[15]-0.6708203932499369*(flowU[14]+flowU[13])+0.5*flowU[12]; 
  flowUOrd[5] = 0.4472135954999581*(flowU[17]+flowU[16])-0.9*flowU[15]+0.6708203932499369*flowU[14]-0.6708203932499369*flowU[13]+0.5*flowU[12]; 
  flowUOrd[6] = (-0.5590169943749475*flowU[17])+0.4472135954999581*flowU[16]+0.6708203932499369*flowU[13]+0.5*flowU[12]; 
  flowUOrd[7] = 0.4472135954999581*(flowU[17]+flowU[16])-0.9*flowU[15]-0.6708203932499369*flowU[14]+0.6708203932499369*flowU[13]+0.5*flowU[12]; 
  flowUOrd[8] = 0.4472135954999581*(flowU[17]+flowU[16])+0.9*flowU[15]+0.6708203932499369*(flowU[14]+flowU[13])+0.5*flowU[12]; 

  vtSqOrd[0] = 0.5*vtSq[0]-0.5590169943749475*(vtSq[5]+vtSq[4]); 
  vtSqOrd[1] = 0.4472135954999581*vtSq[5]-0.5590169943749475*vtSq[4]-0.6708203932499369*vtSq[2]+0.5*vtSq[0]; 
  vtSqOrd[2] = 0.4472135954999581*vtSq[5]-0.5590169943749475*vtSq[4]+0.6708203932499369*vtSq[2]+0.5*vtSq[0]; 
  vtSqOrd[3] = (-0.5590169943749475*vtSq[5])+0.4472135954999581*vtSq[4]-0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[4] = 0.4472135954999581*(vtSq[5]+vtSq[4])+0.9*vtSq[3]-0.6708203932499369*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 
  vtSqOrd[5] = 0.4472135954999581*(vtSq[5]+vtSq[4])-0.9*vtSq[3]+0.6708203932499369*vtSq[2]-0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[6] = (-0.5590169943749475*vtSq[5])+0.4472135954999581*vtSq[4]+0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[7] = 0.4472135954999581*(vtSq[5]+vtSq[4])-0.9*vtSq[3]-0.6708203932499369*vtSq[2]+0.6708203932499369*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[8] = 0.4472135954999581*(vtSq[5]+vtSq[4])+0.9*vtSq[3]+0.6708203932499369*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 

  bmagOrd[0] = 0.5*bmag[0]-0.5590169943749475*(bmag[5]+bmag[4]); 
  bmagOrd[1] = 0.4472135954999581*bmag[5]-0.5590169943749475*bmag[4]-0.6708203932499369*bmag[2]+0.5*bmag[0]; 
  bmagOrd[2] = 0.4472135954999581*bmag[5]-0.5590169943749475*bmag[4]+0.6708203932499369*bmag[2]+0.5*bmag[0]; 
  bmagOrd[3] = (-0.5590169943749475*bmag[5])+0.4472135954999581*bmag[4]-0.6708203932499369*bmag[1]+0.5*bmag[0]; 
  bmagOrd[4] = 0.4472135954999581*(bmag[5]+bmag[4])+0.9*bmag[3]-0.6708203932499369*(bmag[2]+bmag[1])+0.5*bmag[0]; 
  bmagOrd[5] = 0.4472135954999581*(bmag[5]+bmag[4])-0.9*bmag[3]+0.6708203932499369*bmag[2]-0.6708203932499369*bmag[1]+0.5*bmag[0]; 
  bmagOrd[6] = (-0.5590169943749475*bmag[5])+0.4472135954999581*bmag[4]+0.6708203932499369*bmag[1]+0.5*bmag[0]; 
  bmagOrd[7] = 0.4472135954999581*(bmag[5]+bmag[4])-0.9*bmag[3]-0.6708203932499369*bmag[2]+0.6708203932499369*bmag[1]+0.5*bmag[0]; 
  bmagOrd[8] = 0.4472135954999581*(bmag[5]+bmag[4])+0.9*bmag[3]+0.6708203932499369*(bmag[2]+bmag[1])+0.5*bmag[0]; 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (bmagOrd[0]*m0Ord[0])/std::pow(2.506628274631001*sqrt(vtSqOrd[0]),3.0); 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (bmagOrd[1]*m0Ord[1])/std::pow(2.506628274631001*sqrt(vtSqOrd[1]),3.0); 
  else
    fMFacOrd[1] = 0.0;
  if ((vtSqOrd[2] > 0.) && (m0Ord[2] > 0.))
    fMFacOrd[2] = (bmagOrd[2]*m0Ord[2])/std::pow(2.506628274631001*sqrt(vtSqOrd[2]),3.0); 
  else
    fMFacOrd[2] = 0.0;
  if ((vtSqOrd[3] > 0.) && (m0Ord[3] > 0.))
    fMFacOrd[3] = (bmagOrd[3]*m0Ord[3])/std::pow(2.506628274631001*sqrt(vtSqOrd[3]),3.0); 
  else
    fMFacOrd[3] = 0.0;
  if ((vtSqOrd[4] > 0.) && (m0Ord[4] > 0.))
    fMFacOrd[4] = (bmagOrd[4]*m0Ord[4])/std::pow(2.506628274631001*sqrt(vtSqOrd[4]),3.0); 
  else
    fMFacOrd[4] = 0.0;
  if ((vtSqOrd[5] > 0.) && (m0Ord[5] > 0.))
    fMFacOrd[5] = (bmagOrd[5]*m0Ord[5])/std::pow(2.506628274631001*sqrt(vtSqOrd[5]),3.0); 
  else
    fMFacOrd[5] = 0.0;
  if ((vtSqOrd[6] > 0.) && (m0Ord[6] > 0.))
    fMFacOrd[6] = (bmagOrd[6]*m0Ord[6])/std::pow(2.506628274631001*sqrt(vtSqOrd[6]),3.0); 
  else
    fMFacOrd[6] = 0.0;
  if ((vtSqOrd[7] > 0.) && (m0Ord[7] > 0.))
    fMFacOrd[7] = (bmagOrd[7]*m0Ord[7])/std::pow(2.506628274631001*sqrt(vtSqOrd[7]),3.0); 
  else
    fMFacOrd[7] = 0.0;
  if ((vtSqOrd[8] > 0.) && (m0Ord[8] > 0.))
    fMFacOrd[8] = (bmagOrd[8]*m0Ord[8])/std::pow(2.506628274631001*sqrt(vtSqOrd[8]),3.0); 
  else
    fMFacOrd[8] = 0.0;

}

void GkMaxwellianOnBasisGauss2x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[81];
  if ((vtSqOrd[0] > 0.) && (fMFacOrd[0] > 0.)) {
    fMquad[0] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[1] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[2] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[3] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[4] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[5] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[6] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[7] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[8] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
  } else {
    fMquad[0] = 9.999999999999999e-41;
    fMquad[1] = 9.999999999999999e-41;
    fMquad[2] = 9.999999999999999e-41;
    fMquad[3] = 9.999999999999999e-41;
    fMquad[4] = 9.999999999999999e-41;
    fMquad[5] = 9.999999999999999e-41;
    fMquad[6] = 9.999999999999999e-41;
    fMquad[7] = 9.999999999999999e-41;
    fMquad[8] = 9.999999999999999e-41;
  };
  if ((vtSqOrd[1] > 0.) && (fMFacOrd[1] > 0.)) {
    fMquad[9] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[10] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[11] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[12] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[13] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[14] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[15] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[16] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[17] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
  } else {
    fMquad[9] = 9.999999999999999e-41;
    fMquad[10] = 9.999999999999999e-41;
    fMquad[11] = 9.999999999999999e-41;
    fMquad[12] = 9.999999999999999e-41;
    fMquad[13] = 9.999999999999999e-41;
    fMquad[14] = 9.999999999999999e-41;
    fMquad[15] = 9.999999999999999e-41;
    fMquad[16] = 9.999999999999999e-41;
    fMquad[17] = 9.999999999999999e-41;
  };
  if ((vtSqOrd[2] > 0.) && (fMFacOrd[2] > 0.)) {
    fMquad[18] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[19] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[20] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[21] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[22] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[23] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[24] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[25] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[26] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
  } else {
    fMquad[18] = 9.999999999999999e-41;
    fMquad[19] = 9.999999999999999e-41;
    fMquad[20] = 9.999999999999999e-41;
    fMquad[21] = 9.999999999999999e-41;
    fMquad[22] = 9.999999999999999e-41;
    fMquad[23] = 9.999999999999999e-41;
    fMquad[24] = 9.999999999999999e-41;
    fMquad[25] = 9.999999999999999e-41;
    fMquad[26] = 9.999999999999999e-41;
  };
  if ((vtSqOrd[3] > 0.) && (fMFacOrd[3] > 0.)) {
    fMquad[27] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[3],2.0))/vtSqOrd[3])+9.999999999999999e-41; 
    fMquad[28] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[3],2.0))/vtSqOrd[3])+9.999999999999999e-41; 
    fMquad[29] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[3],2.0))/vtSqOrd[3])+9.999999999999999e-41; 
    fMquad[30] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[3])+9.999999999999999e-41; 
    fMquad[31] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[3])+9.999999999999999e-41; 
    fMquad[32] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[3])+9.999999999999999e-41; 
    fMquad[33] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[3])+9.999999999999999e-41; 
    fMquad[34] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[3])+9.999999999999999e-41; 
    fMquad[35] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[3])+9.999999999999999e-41; 
  } else {
    fMquad[27] = 9.999999999999999e-41;
    fMquad[28] = 9.999999999999999e-41;
    fMquad[29] = 9.999999999999999e-41;
    fMquad[30] = 9.999999999999999e-41;
    fMquad[31] = 9.999999999999999e-41;
    fMquad[32] = 9.999999999999999e-41;
    fMquad[33] = 9.999999999999999e-41;
    fMquad[34] = 9.999999999999999e-41;
    fMquad[35] = 9.999999999999999e-41;
  };
  if ((vtSqOrd[4] > 0.) && (fMFacOrd[4] > 0.)) {
    fMquad[36] = fMFacOrd[4]*exp(((-(1.0*bmagOrd[4]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[4],2.0))/vtSqOrd[4])+9.999999999999999e-41; 
    fMquad[37] = fMFacOrd[4]*exp(((-(1.0*bmagOrd[4]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[4],2.0))/vtSqOrd[4])+9.999999999999999e-41; 
    fMquad[38] = fMFacOrd[4]*exp(((-(1.0*bmagOrd[4]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[4],2.0))/vtSqOrd[4])+9.999999999999999e-41; 
    fMquad[39] = fMFacOrd[4]*exp(((-(1.0*bmagOrd[4]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[4])+9.999999999999999e-41; 
    fMquad[40] = fMFacOrd[4]*exp(((-(1.0*bmagOrd[4]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[4])+9.999999999999999e-41; 
    fMquad[41] = fMFacOrd[4]*exp(((-(1.0*bmagOrd[4]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[4])+9.999999999999999e-41; 
    fMquad[42] = fMFacOrd[4]*exp(((-(1.0*bmagOrd[4]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[4])+9.999999999999999e-41; 
    fMquad[43] = fMFacOrd[4]*exp(((-(1.0*bmagOrd[4]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[4])+9.999999999999999e-41; 
    fMquad[44] = fMFacOrd[4]*exp(((-(1.0*bmagOrd[4]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[4])+9.999999999999999e-41; 
  } else {
    fMquad[36] = 9.999999999999999e-41;
    fMquad[37] = 9.999999999999999e-41;
    fMquad[38] = 9.999999999999999e-41;
    fMquad[39] = 9.999999999999999e-41;
    fMquad[40] = 9.999999999999999e-41;
    fMquad[41] = 9.999999999999999e-41;
    fMquad[42] = 9.999999999999999e-41;
    fMquad[43] = 9.999999999999999e-41;
    fMquad[44] = 9.999999999999999e-41;
  };
  if ((vtSqOrd[5] > 0.) && (fMFacOrd[5] > 0.)) {
    fMquad[45] = fMFacOrd[5]*exp(((-(1.0*bmagOrd[5]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[5],2.0))/vtSqOrd[5])+9.999999999999999e-41; 
    fMquad[46] = fMFacOrd[5]*exp(((-(1.0*bmagOrd[5]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[5],2.0))/vtSqOrd[5])+9.999999999999999e-41; 
    fMquad[47] = fMFacOrd[5]*exp(((-(1.0*bmagOrd[5]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[5],2.0))/vtSqOrd[5])+9.999999999999999e-41; 
    fMquad[48] = fMFacOrd[5]*exp(((-(1.0*bmagOrd[5]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[5])+9.999999999999999e-41; 
    fMquad[49] = fMFacOrd[5]*exp(((-(1.0*bmagOrd[5]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[5])+9.999999999999999e-41; 
    fMquad[50] = fMFacOrd[5]*exp(((-(1.0*bmagOrd[5]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[5])+9.999999999999999e-41; 
    fMquad[51] = fMFacOrd[5]*exp(((-(1.0*bmagOrd[5]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[5])+9.999999999999999e-41; 
    fMquad[52] = fMFacOrd[5]*exp(((-(1.0*bmagOrd[5]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[5])+9.999999999999999e-41; 
    fMquad[53] = fMFacOrd[5]*exp(((-(1.0*bmagOrd[5]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[5])+9.999999999999999e-41; 
  } else {
    fMquad[45] = 9.999999999999999e-41;
    fMquad[46] = 9.999999999999999e-41;
    fMquad[47] = 9.999999999999999e-41;
    fMquad[48] = 9.999999999999999e-41;
    fMquad[49] = 9.999999999999999e-41;
    fMquad[50] = 9.999999999999999e-41;
    fMquad[51] = 9.999999999999999e-41;
    fMquad[52] = 9.999999999999999e-41;
    fMquad[53] = 9.999999999999999e-41;
  };
  if ((vtSqOrd[6] > 0.) && (fMFacOrd[6] > 0.)) {
    fMquad[54] = fMFacOrd[6]*exp(((-(1.0*bmagOrd[6]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[6],2.0))/vtSqOrd[6])+9.999999999999999e-41; 
    fMquad[55] = fMFacOrd[6]*exp(((-(1.0*bmagOrd[6]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[6],2.0))/vtSqOrd[6])+9.999999999999999e-41; 
    fMquad[56] = fMFacOrd[6]*exp(((-(1.0*bmagOrd[6]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[6],2.0))/vtSqOrd[6])+9.999999999999999e-41; 
    fMquad[57] = fMFacOrd[6]*exp(((-(1.0*bmagOrd[6]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[6])+9.999999999999999e-41; 
    fMquad[58] = fMFacOrd[6]*exp(((-(1.0*bmagOrd[6]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[6])+9.999999999999999e-41; 
    fMquad[59] = fMFacOrd[6]*exp(((-(1.0*bmagOrd[6]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[6])+9.999999999999999e-41; 
    fMquad[60] = fMFacOrd[6]*exp(((-(1.0*bmagOrd[6]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[6])+9.999999999999999e-41; 
    fMquad[61] = fMFacOrd[6]*exp(((-(1.0*bmagOrd[6]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[6])+9.999999999999999e-41; 
    fMquad[62] = fMFacOrd[6]*exp(((-(1.0*bmagOrd[6]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[6])+9.999999999999999e-41; 
  } else {
    fMquad[54] = 9.999999999999999e-41;
    fMquad[55] = 9.999999999999999e-41;
    fMquad[56] = 9.999999999999999e-41;
    fMquad[57] = 9.999999999999999e-41;
    fMquad[58] = 9.999999999999999e-41;
    fMquad[59] = 9.999999999999999e-41;
    fMquad[60] = 9.999999999999999e-41;
    fMquad[61] = 9.999999999999999e-41;
    fMquad[62] = 9.999999999999999e-41;
  };
  if ((vtSqOrd[7] > 0.) && (fMFacOrd[7] > 0.)) {
    fMquad[63] = fMFacOrd[7]*exp(((-(1.0*bmagOrd[7]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[7],2.0))/vtSqOrd[7])+9.999999999999999e-41; 
    fMquad[64] = fMFacOrd[7]*exp(((-(1.0*bmagOrd[7]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[7],2.0))/vtSqOrd[7])+9.999999999999999e-41; 
    fMquad[65] = fMFacOrd[7]*exp(((-(1.0*bmagOrd[7]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[7],2.0))/vtSqOrd[7])+9.999999999999999e-41; 
    fMquad[66] = fMFacOrd[7]*exp(((-(1.0*bmagOrd[7]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[7])+9.999999999999999e-41; 
    fMquad[67] = fMFacOrd[7]*exp(((-(1.0*bmagOrd[7]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[7])+9.999999999999999e-41; 
    fMquad[68] = fMFacOrd[7]*exp(((-(1.0*bmagOrd[7]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[7])+9.999999999999999e-41; 
    fMquad[69] = fMFacOrd[7]*exp(((-(1.0*bmagOrd[7]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[7])+9.999999999999999e-41; 
    fMquad[70] = fMFacOrd[7]*exp(((-(1.0*bmagOrd[7]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[7])+9.999999999999999e-41; 
    fMquad[71] = fMFacOrd[7]*exp(((-(1.0*bmagOrd[7]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[7])+9.999999999999999e-41; 
  } else {
    fMquad[63] = 9.999999999999999e-41;
    fMquad[64] = 9.999999999999999e-41;
    fMquad[65] = 9.999999999999999e-41;
    fMquad[66] = 9.999999999999999e-41;
    fMquad[67] = 9.999999999999999e-41;
    fMquad[68] = 9.999999999999999e-41;
    fMquad[69] = 9.999999999999999e-41;
    fMquad[70] = 9.999999999999999e-41;
    fMquad[71] = 9.999999999999999e-41;
  };
  if ((vtSqOrd[8] > 0.) && (fMFacOrd[8] > 0.)) {
    fMquad[72] = fMFacOrd[8]*exp(((-(1.0*bmagOrd[8]*std::abs(wc[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[8],2.0))/vtSqOrd[8])+9.999999999999999e-41; 
    fMquad[73] = fMFacOrd[8]*exp(((-(1.0*bmagOrd[8]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[8],2.0))/vtSqOrd[8])+9.999999999999999e-41; 
    fMquad[74] = fMFacOrd[8]*exp(((-(1.0*bmagOrd[8]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[8],2.0))/vtSqOrd[8])+9.999999999999999e-41; 
    fMquad[75] = fMFacOrd[8]*exp(((-(1.0*bmagOrd[8]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[8])+9.999999999999999e-41; 
    fMquad[76] = fMFacOrd[8]*exp(((-(1.0*bmagOrd[8]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[8])+9.999999999999999e-41; 
    fMquad[77] = fMFacOrd[8]*exp(((-(1.0*bmagOrd[8]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[8])+9.999999999999999e-41; 
    fMquad[78] = fMFacOrd[8]*exp(((-(1.0*bmagOrd[8]*std::abs(wc[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[8])+9.999999999999999e-41; 
    fMquad[79] = fMFacOrd[8]*exp(((-(1.0*bmagOrd[8]*std::abs(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[8])+9.999999999999999e-41; 
    fMquad[80] = fMFacOrd[8]*exp(((-(1.0*bmagOrd[8]*std::abs(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[8])+9.999999999999999e-41; 
  } else {
    fMquad[72] = 9.999999999999999e-41;
    fMquad[73] = 9.999999999999999e-41;
    fMquad[74] = 9.999999999999999e-41;
    fMquad[75] = 9.999999999999999e-41;
    fMquad[76] = 9.999999999999999e-41;
    fMquad[77] = 9.999999999999999e-41;
    fMquad[78] = 9.999999999999999e-41;
    fMquad[79] = 9.999999999999999e-41;
    fMquad[80] = 9.999999999999999e-41;
  };

  fMOut[0] = 0.02381496723060508*(fMquad[80]+fMquad[79])+0.03810394756896813*fMquad[78]+0.02381496723060508*(fMquad[77]+fMquad[76])+0.03810394756896813*fMquad[75]+0.03810394756896812*(fMquad[74]+fMquad[73])+0.06096631611034901*fMquad[72]+0.02381496723060508*(fMquad[71]+fMquad[70])+0.03810394756896813*fMquad[69]+0.02381496723060508*(fMquad[68]+fMquad[67])+0.03810394756896813*fMquad[66]+0.03810394756896812*(fMquad[65]+fMquad[64])+0.06096631611034901*fMquad[63]+0.03810394756896812*(fMquad[62]+fMquad[61])+0.06096631611034901*fMquad[60]+0.03810394756896812*(fMquad[59]+fMquad[58])+0.06096631611034901*fMquad[57]+0.06096631611034901*(fMquad[56]+fMquad[55])+0.09754610577655842*fMquad[54]+0.02381496723060508*(fMquad[53]+fMquad[52])+0.03810394756896813*fMquad[51]+0.02381496723060508*(fMquad[50]+fMquad[49])+0.03810394756896813*fMquad[48]+0.03810394756896812*(fMquad[47]+fMquad[46])+0.06096631611034901*fMquad[45]+0.02381496723060508*(fMquad[44]+fMquad[43])+0.03810394756896813*fMquad[42]+0.02381496723060508*(fMquad[41]+fMquad[40])+0.03810394756896813*fMquad[39]+0.03810394756896812*(fMquad[38]+fMquad[37])+0.06096631611034901*fMquad[36]+0.03810394756896812*(fMquad[35]+fMquad[34])+0.06096631611034901*fMquad[33]+0.03810394756896812*(fMquad[32]+fMquad[31])+0.06096631611034901*fMquad[30]+0.06096631611034901*(fMquad[29]+fMquad[28])+0.09754610577655842*fMquad[27]+0.03810394756896812*(fMquad[26]+fMquad[25])+0.06096631611034901*fMquad[24]+0.03810394756896812*(fMquad[23]+fMquad[22])+0.06096631611034901*fMquad[21]+0.06096631611034901*(fMquad[20]+fMquad[19])+0.09754610577655842*fMquad[18]+0.03810394756896812*(fMquad[17]+fMquad[16])+0.06096631611034901*fMquad[15]+0.03810394756896812*(fMquad[14]+fMquad[13])+0.06096631611034901*fMquad[12]+0.06096631611034901*(fMquad[11]+fMquad[10])+0.09754610577655842*fMquad[9]+0.06096631611034901*(fMquad[8]+fMquad[7])+0.09754610577655842*fMquad[6]+0.06096631611034901*(fMquad[5]+fMquad[4])+0.09754610577655842*(fMquad[3]+fMquad[2]+fMquad[1])+0.1560737692424935*fMquad[0]; 
  fMOut[1] = 0.03195113136573771*(fMquad[80]+fMquad[79])+0.05112181018518035*fMquad[78]+0.03195113136573771*(fMquad[77]+fMquad[76])+0.05112181018518035*fMquad[75]+0.05112181018518035*(fMquad[74]+fMquad[73])+0.08179489629628857*fMquad[72]+0.03195113136573771*(fMquad[71]+fMquad[70])+0.05112181018518035*fMquad[69]+0.03195113136573771*(fMquad[68]+fMquad[67])+0.05112181018518035*fMquad[66]+0.05112181018518035*(fMquad[65]+fMquad[64])+0.08179489629628857*fMquad[63]+0.05112181018518035*(fMquad[62]+fMquad[61])+0.08179489629628857*fMquad[60]+0.05112181018518035*(fMquad[59]+fMquad[58])+0.08179489629628857*(fMquad[57]+fMquad[56]+fMquad[55])+0.1308718340740617*fMquad[54]-0.03195113136573771*(fMquad[53]+fMquad[52])-0.05112181018518035*fMquad[51]-0.03195113136573771*(fMquad[50]+fMquad[49])-0.05112181018518035*fMquad[48]-0.05112181018518035*(fMquad[47]+fMquad[46])-0.08179489629628857*fMquad[45]-0.03195113136573771*(fMquad[44]+fMquad[43])-0.05112181018518035*fMquad[42]-0.03195113136573771*(fMquad[41]+fMquad[40])-0.05112181018518035*fMquad[39]-0.05112181018518035*(fMquad[38]+fMquad[37])-0.08179489629628857*fMquad[36]-0.05112181018518035*(fMquad[35]+fMquad[34])-0.08179489629628857*fMquad[33]-0.05112181018518035*(fMquad[32]+fMquad[31])-0.08179489629628857*(fMquad[30]+fMquad[29]+fMquad[28])-0.1308718340740617*fMquad[27]; 
  fMOut[2] = 0.03195113136573771*(fMquad[80]+fMquad[79])+0.05112181018518035*fMquad[78]+0.03195113136573771*(fMquad[77]+fMquad[76])+0.05112181018518035*fMquad[75]+0.05112181018518035*(fMquad[74]+fMquad[73])+0.08179489629628857*fMquad[72]-0.03195113136573771*(fMquad[71]+fMquad[70])-0.05112181018518035*fMquad[69]-0.03195113136573771*(fMquad[68]+fMquad[67])-0.05112181018518035*fMquad[66]-0.05112181018518035*(fMquad[65]+fMquad[64])-0.08179489629628857*fMquad[63]+0.03195113136573771*(fMquad[53]+fMquad[52])+0.05112181018518035*fMquad[51]+0.03195113136573771*(fMquad[50]+fMquad[49])+0.05112181018518035*fMquad[48]+0.05112181018518035*(fMquad[47]+fMquad[46])+0.08179489629628857*fMquad[45]-0.03195113136573771*(fMquad[44]+fMquad[43])-0.05112181018518035*fMquad[42]-0.03195113136573771*(fMquad[41]+fMquad[40])-0.05112181018518035*fMquad[39]-0.05112181018518035*(fMquad[38]+fMquad[37])-0.08179489629628857*fMquad[36]+0.05112181018518035*(fMquad[26]+fMquad[25])+0.08179489629628857*fMquad[24]+0.05112181018518035*(fMquad[23]+fMquad[22])+0.08179489629628857*(fMquad[21]+fMquad[20]+fMquad[19])+0.1308718340740617*fMquad[18]-0.05112181018518035*(fMquad[17]+fMquad[16])-0.08179489629628857*fMquad[15]-0.05112181018518035*(fMquad[14]+fMquad[13])-0.08179489629628857*(fMquad[12]+fMquad[11]+fMquad[10])-0.1308718340740617*fMquad[9]; 
  fMOut[3] = 0.03195113136573771*(fMquad[80]+fMquad[79])+0.05112181018518035*fMquad[78]-0.03195113136573771*(fMquad[77]+fMquad[76])-0.05112181018518035*fMquad[75]+0.03195113136573771*(fMquad[71]+fMquad[70])+0.05112181018518035*fMquad[69]-0.03195113136573771*(fMquad[68]+fMquad[67])-0.05112181018518035*fMquad[66]+0.05112181018518035*(fMquad[62]+fMquad[61])+0.08179489629628857*fMquad[60]-0.05112181018518035*(fMquad[59]+fMquad[58])-0.08179489629628857*fMquad[57]+0.03195113136573771*(fMquad[53]+fMquad[52])+0.05112181018518035*fMquad[51]-0.03195113136573771*(fMquad[50]+fMquad[49])-0.05112181018518035*fMquad[48]+0.03195113136573771*(fMquad[44]+fMquad[43])+0.05112181018518035*fMquad[42]-0.03195113136573771*(fMquad[41]+fMquad[40])-0.05112181018518035*fMquad[39]+0.05112181018518035*(fMquad[35]+fMquad[34])+0.08179489629628857*fMquad[33]-0.05112181018518035*(fMquad[32]+fMquad[31])-0.08179489629628857*fMquad[30]+0.05112181018518035*(fMquad[26]+fMquad[25])+0.08179489629628857*fMquad[24]-0.05112181018518035*(fMquad[23]+fMquad[22])-0.08179489629628857*fMquad[21]+0.05112181018518035*(fMquad[17]+fMquad[16])+0.08179489629628857*fMquad[15]-0.05112181018518035*(fMquad[14]+fMquad[13])-0.08179489629628857*fMquad[12]+0.08179489629628857*(fMquad[8]+fMquad[7])+0.1308718340740617*fMquad[6]-0.08179489629628857*(fMquad[5]+fMquad[4])-0.1308718340740617*fMquad[3]; 
  fMOut[4] = 0.03195113136573771*fMquad[80]-0.03195113136573771*fMquad[79]+0.03195113136573771*fMquad[77]-0.03195113136573771*fMquad[76]+0.05112181018518035*fMquad[74]-0.05112181018518035*fMquad[73]+0.03195113136573771*fMquad[71]-0.03195113136573771*fMquad[70]+0.03195113136573771*fMquad[68]-0.03195113136573771*fMquad[67]+0.05112181018518035*fMquad[65]-0.05112181018518035*fMquad[64]+0.05112181018518035*fMquad[62]-0.05112181018518035*fMquad[61]+0.05112181018518035*fMquad[59]-0.05112181018518035*fMquad[58]+0.08179489629628857*fMquad[56]-0.08179489629628857*fMquad[55]+0.03195113136573771*fMquad[53]-0.03195113136573771*fMquad[52]+0.03195113136573771*fMquad[50]-0.03195113136573771*fMquad[49]+0.05112181018518035*fMquad[47]-0.05112181018518035*fMquad[46]+0.03195113136573771*fMquad[44]-0.03195113136573771*fMquad[43]+0.03195113136573771*fMquad[41]-0.03195113136573771*fMquad[40]+0.05112181018518035*fMquad[38]-0.05112181018518035*fMquad[37]+0.05112181018518035*fMquad[35]-0.05112181018518035*fMquad[34]+0.05112181018518035*fMquad[32]-0.05112181018518035*fMquad[31]+0.08179489629628857*fMquad[29]-0.08179489629628857*fMquad[28]+0.05112181018518035*fMquad[26]-0.05112181018518035*fMquad[25]+0.05112181018518035*fMquad[23]-0.05112181018518035*fMquad[22]+0.08179489629628857*fMquad[20]-0.08179489629628857*fMquad[19]+0.05112181018518035*fMquad[17]-0.05112181018518035*fMquad[16]+0.05112181018518035*fMquad[14]-0.05112181018518035*fMquad[13]+0.08179489629628857*fMquad[11]-0.08179489629628857*fMquad[10]+0.08179489629628857*fMquad[8]-0.08179489629628857*fMquad[7]+0.08179489629628857*fMquad[5]-0.08179489629628857*fMquad[4]+0.1308718340740617*fMquad[2]-0.1308718340740617*fMquad[1]; 
  fMOut[5] = 0.04286694101508914*(fMquad[80]+fMquad[79])+0.06858710562414264*fMquad[78]+0.04286694101508914*(fMquad[77]+fMquad[76])+0.06858710562414264*fMquad[75]+0.06858710562414262*(fMquad[74]+fMquad[73])+0.1097393689986282*fMquad[72]-0.04286694101508914*(fMquad[71]+fMquad[70])-0.06858710562414264*fMquad[69]-0.04286694101508914*(fMquad[68]+fMquad[67])-0.06858710562414264*fMquad[66]-0.06858710562414262*(fMquad[65]+fMquad[64])-0.1097393689986282*fMquad[63]-0.04286694101508914*(fMquad[53]+fMquad[52])-0.06858710562414264*fMquad[51]-0.04286694101508914*(fMquad[50]+fMquad[49])-0.06858710562414264*fMquad[48]-0.06858710562414262*(fMquad[47]+fMquad[46])-0.1097393689986282*fMquad[45]+0.04286694101508914*(fMquad[44]+fMquad[43])+0.06858710562414264*fMquad[42]+0.04286694101508914*(fMquad[41]+fMquad[40])+0.06858710562414264*fMquad[39]+0.06858710562414262*(fMquad[38]+fMquad[37])+0.1097393689986282*fMquad[36]; 
  fMOut[6] = 0.04286694101508914*(fMquad[80]+fMquad[79])+0.06858710562414264*fMquad[78]-0.04286694101508914*(fMquad[77]+fMquad[76])-0.06858710562414264*fMquad[75]+0.04286694101508914*(fMquad[71]+fMquad[70])+0.06858710562414264*fMquad[69]-0.04286694101508914*(fMquad[68]+fMquad[67])-0.06858710562414264*fMquad[66]+0.06858710562414262*(fMquad[62]+fMquad[61])+0.1097393689986282*fMquad[60]-0.06858710562414262*(fMquad[59]+fMquad[58])-0.1097393689986282*fMquad[57]-0.04286694101508914*(fMquad[53]+fMquad[52])-0.06858710562414264*fMquad[51]+0.04286694101508914*(fMquad[50]+fMquad[49])+0.06858710562414264*fMquad[48]-0.04286694101508914*(fMquad[44]+fMquad[43])-0.06858710562414264*fMquad[42]+0.04286694101508914*(fMquad[41]+fMquad[40])+0.06858710562414264*fMquad[39]-0.06858710562414262*(fMquad[35]+fMquad[34])-0.1097393689986282*fMquad[33]+0.06858710562414262*(fMquad[32]+fMquad[31])+0.1097393689986282*fMquad[30]; 
  fMOut[7] = 0.04286694101508914*(fMquad[80]+fMquad[79])+0.06858710562414264*fMquad[78]-0.04286694101508914*(fMquad[77]+fMquad[76])-0.06858710562414264*fMquad[75]-0.04286694101508914*(fMquad[71]+fMquad[70])-0.06858710562414264*fMquad[69]+0.04286694101508914*(fMquad[68]+fMquad[67])+0.06858710562414264*fMquad[66]+0.04286694101508914*(fMquad[53]+fMquad[52])+0.06858710562414264*fMquad[51]-0.04286694101508914*(fMquad[50]+fMquad[49])-0.06858710562414264*fMquad[48]-0.04286694101508914*(fMquad[44]+fMquad[43])-0.06858710562414264*fMquad[42]+0.04286694101508914*(fMquad[41]+fMquad[40])+0.06858710562414264*fMquad[39]+0.06858710562414262*(fMquad[26]+fMquad[25])+0.1097393689986282*fMquad[24]-0.06858710562414262*(fMquad[23]+fMquad[22])-0.1097393689986282*fMquad[21]-0.06858710562414262*(fMquad[17]+fMquad[16])-0.1097393689986282*fMquad[15]+0.06858710562414262*(fMquad[14]+fMquad[13])+0.1097393689986282*fMquad[12]; 
  fMOut[8] = 0.04286694101508914*fMquad[80]-0.04286694101508914*fMquad[79]+0.04286694101508914*fMquad[77]-0.04286694101508914*fMquad[76]+0.06858710562414262*fMquad[74]-0.06858710562414262*fMquad[73]+0.04286694101508914*fMquad[71]-0.04286694101508914*fMquad[70]+0.04286694101508914*fMquad[68]-0.04286694101508914*fMquad[67]+0.06858710562414262*fMquad[65]-0.06858710562414262*fMquad[64]+0.06858710562414262*fMquad[62]-0.06858710562414262*fMquad[61]+0.06858710562414262*fMquad[59]-0.06858710562414262*fMquad[58]+0.1097393689986282*fMquad[56]-0.1097393689986282*fMquad[55]-0.04286694101508914*fMquad[53]+0.04286694101508914*fMquad[52]-0.04286694101508914*fMquad[50]+0.04286694101508914*fMquad[49]-0.06858710562414262*fMquad[47]+0.06858710562414262*fMquad[46]-0.04286694101508914*fMquad[44]+0.04286694101508914*fMquad[43]-0.04286694101508914*fMquad[41]+0.04286694101508914*fMquad[40]-0.06858710562414262*fMquad[38]+0.06858710562414262*fMquad[37]-0.06858710562414262*fMquad[35]+0.06858710562414262*fMquad[34]-0.06858710562414262*fMquad[32]+0.06858710562414262*fMquad[31]-0.1097393689986282*fMquad[29]+0.1097393689986282*fMquad[28]; 
  fMOut[9] = 0.04286694101508914*fMquad[80]-0.04286694101508914*fMquad[79]+0.04286694101508914*fMquad[77]-0.04286694101508914*fMquad[76]+0.06858710562414262*fMquad[74]-0.06858710562414262*fMquad[73]-0.04286694101508914*fMquad[71]+0.04286694101508914*fMquad[70]-0.04286694101508914*fMquad[68]+0.04286694101508914*fMquad[67]-0.06858710562414262*fMquad[65]+0.06858710562414262*fMquad[64]+0.04286694101508914*fMquad[53]-0.04286694101508914*fMquad[52]+0.04286694101508914*fMquad[50]-0.04286694101508914*fMquad[49]+0.06858710562414262*fMquad[47]-0.06858710562414262*fMquad[46]-0.04286694101508914*fMquad[44]+0.04286694101508914*fMquad[43]-0.04286694101508914*fMquad[41]+0.04286694101508914*fMquad[40]-0.06858710562414262*fMquad[38]+0.06858710562414262*(fMquad[37]+fMquad[26])-0.06858710562414262*fMquad[25]+0.06858710562414262*fMquad[23]-0.06858710562414262*fMquad[22]+0.1097393689986282*fMquad[20]-0.1097393689986282*fMquad[19]-0.06858710562414262*fMquad[17]+0.06858710562414262*fMquad[16]-0.06858710562414262*fMquad[14]+0.06858710562414262*fMquad[13]-0.1097393689986282*fMquad[11]+0.1097393689986282*fMquad[10]; 
  fMOut[10] = 0.04286694101508914*fMquad[80]-0.04286694101508914*(fMquad[79]+fMquad[77])+0.04286694101508914*(fMquad[76]+fMquad[71])-0.04286694101508914*(fMquad[70]+fMquad[68])+0.04286694101508914*fMquad[67]+0.06858710562414262*fMquad[62]-0.06858710562414262*(fMquad[61]+fMquad[59])+0.06858710562414262*fMquad[58]+0.04286694101508914*fMquad[53]-0.04286694101508914*(fMquad[52]+fMquad[50])+0.04286694101508914*(fMquad[49]+fMquad[44])-0.04286694101508914*(fMquad[43]+fMquad[41])+0.04286694101508914*fMquad[40]+0.06858710562414262*fMquad[35]-0.06858710562414262*(fMquad[34]+fMquad[32])+0.06858710562414262*(fMquad[31]+fMquad[26])-0.06858710562414262*(fMquad[25]+fMquad[23])+0.06858710562414262*(fMquad[22]+fMquad[17])-0.06858710562414262*(fMquad[16]+fMquad[14])+0.06858710562414262*fMquad[13]+0.1097393689986282*fMquad[8]-0.1097393689986282*(fMquad[7]+fMquad[5])+0.1097393689986282*fMquad[4]; 
  fMOut[11] = 0.02130075424382515*(fMquad[80]+fMquad[79])+0.03408120679012025*fMquad[78]+0.02130075424382515*(fMquad[77]+fMquad[76])+0.03408120679012025*fMquad[75]+0.03408120679012024*(fMquad[74]+fMquad[73])+0.0545299308641924*fMquad[72]+0.02130075424382515*(fMquad[71]+fMquad[70])+0.03408120679012025*fMquad[69]+0.02130075424382515*(fMquad[68]+fMquad[67])+0.03408120679012025*fMquad[66]+0.03408120679012024*(fMquad[65]+fMquad[64])+0.0545299308641924*fMquad[63]+0.03408120679012024*(fMquad[62]+fMquad[61])+0.0545299308641924*fMquad[60]+0.03408120679012024*(fMquad[59]+fMquad[58])+0.0545299308641924*(fMquad[57]+fMquad[56]+fMquad[55])+0.08724788938270785*fMquad[54]+0.02130075424382515*(fMquad[53]+fMquad[52])+0.03408120679012025*fMquad[51]+0.02130075424382515*(fMquad[50]+fMquad[49])+0.03408120679012025*fMquad[48]+0.03408120679012024*(fMquad[47]+fMquad[46])+0.0545299308641924*fMquad[45]+0.02130075424382515*(fMquad[44]+fMquad[43])+0.03408120679012025*fMquad[42]+0.02130075424382515*(fMquad[41]+fMquad[40])+0.03408120679012025*fMquad[39]+0.03408120679012024*(fMquad[38]+fMquad[37])+0.0545299308641924*fMquad[36]+0.03408120679012024*(fMquad[35]+fMquad[34])+0.0545299308641924*fMquad[33]+0.03408120679012024*(fMquad[32]+fMquad[31])+0.0545299308641924*(fMquad[30]+fMquad[29]+fMquad[28])+0.08724788938270785*fMquad[27]-0.04260150848765029*(fMquad[26]+fMquad[25])-0.06816241358024047*fMquad[24]-0.04260150848765029*(fMquad[23]+fMquad[22])-0.06816241358024047*fMquad[21]-0.06816241358024049*(fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]-0.04260150848765029*(fMquad[17]+fMquad[16])-0.06816241358024047*fMquad[15]-0.04260150848765029*(fMquad[14]+fMquad[13])-0.06816241358024047*fMquad[12]-0.06816241358024049*(fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]-0.06816241358024049*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]-0.06816241358024049*(fMquad[5]+fMquad[4])-0.1090598617283848*(fMquad[3]+fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[12] = 0.02130075424382515*(fMquad[80]+fMquad[79])+0.03408120679012025*fMquad[78]+0.02130075424382515*(fMquad[77]+fMquad[76])+0.03408120679012025*fMquad[75]+0.03408120679012024*(fMquad[74]+fMquad[73])+0.0545299308641924*fMquad[72]+0.02130075424382515*(fMquad[71]+fMquad[70])+0.03408120679012025*fMquad[69]+0.02130075424382515*(fMquad[68]+fMquad[67])+0.03408120679012025*fMquad[66]+0.03408120679012024*(fMquad[65]+fMquad[64])+0.0545299308641924*fMquad[63]-0.04260150848765029*(fMquad[62]+fMquad[61])-0.06816241358024047*fMquad[60]-0.04260150848765029*(fMquad[59]+fMquad[58])-0.06816241358024047*fMquad[57]-0.06816241358024049*(fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382515*(fMquad[53]+fMquad[52])+0.03408120679012025*fMquad[51]+0.02130075424382515*(fMquad[50]+fMquad[49])+0.03408120679012025*fMquad[48]+0.03408120679012024*(fMquad[47]+fMquad[46])+0.0545299308641924*fMquad[45]+0.02130075424382515*(fMquad[44]+fMquad[43])+0.03408120679012025*fMquad[42]+0.02130075424382515*(fMquad[41]+fMquad[40])+0.03408120679012025*fMquad[39]+0.03408120679012024*(fMquad[38]+fMquad[37])+0.0545299308641924*fMquad[36]-0.04260150848765029*(fMquad[35]+fMquad[34])-0.06816241358024047*fMquad[33]-0.04260150848765029*(fMquad[32]+fMquad[31])-0.06816241358024047*fMquad[30]-0.06816241358024049*(fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012024*(fMquad[26]+fMquad[25])+0.0545299308641924*fMquad[24]+0.03408120679012024*(fMquad[23]+fMquad[22])+0.0545299308641924*(fMquad[21]+fMquad[20]+fMquad[19])+0.08724788938270785*fMquad[18]+0.03408120679012024*(fMquad[17]+fMquad[16])+0.0545299308641924*fMquad[15]+0.03408120679012024*(fMquad[14]+fMquad[13])+0.0545299308641924*(fMquad[12]+fMquad[11]+fMquad[10])+0.08724788938270785*fMquad[9]-0.06816241358024049*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]-0.06816241358024049*(fMquad[5]+fMquad[4])-0.1090598617283848*(fMquad[3]+fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[13] = 0.02130075424382515*(fMquad[80]+fMquad[79])+0.03408120679012025*fMquad[78]+0.02130075424382515*(fMquad[77]+fMquad[76])+0.03408120679012025*fMquad[75]-0.04260150848765029*(fMquad[74]+fMquad[73])-0.06816241358024047*fMquad[72]+0.02130075424382515*(fMquad[71]+fMquad[70])+0.03408120679012025*fMquad[69]+0.02130075424382515*(fMquad[68]+fMquad[67])+0.03408120679012025*fMquad[66]-0.04260150848765029*(fMquad[65]+fMquad[64])-0.06816241358024047*fMquad[63]+0.03408120679012024*(fMquad[62]+fMquad[61])+0.0545299308641924*fMquad[60]+0.03408120679012024*(fMquad[59]+fMquad[58])+0.0545299308641924*fMquad[57]-0.06816241358024049*(fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382515*(fMquad[53]+fMquad[52])+0.03408120679012025*fMquad[51]+0.02130075424382515*(fMquad[50]+fMquad[49])+0.03408120679012025*fMquad[48]-0.04260150848765029*(fMquad[47]+fMquad[46])-0.06816241358024047*fMquad[45]+0.02130075424382515*(fMquad[44]+fMquad[43])+0.03408120679012025*fMquad[42]+0.02130075424382515*(fMquad[41]+fMquad[40])+0.03408120679012025*fMquad[39]-0.04260150848765029*(fMquad[38]+fMquad[37])-0.06816241358024047*fMquad[36]+0.03408120679012024*(fMquad[35]+fMquad[34])+0.0545299308641924*fMquad[33]+0.03408120679012024*(fMquad[32]+fMquad[31])+0.0545299308641924*fMquad[30]-0.06816241358024049*(fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012024*(fMquad[26]+fMquad[25])+0.0545299308641924*fMquad[24]+0.03408120679012024*(fMquad[23]+fMquad[22])+0.0545299308641924*fMquad[21]-0.06816241358024049*(fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]+0.03408120679012024*(fMquad[17]+fMquad[16])+0.0545299308641924*fMquad[15]+0.03408120679012024*(fMquad[14]+fMquad[13])+0.0545299308641924*fMquad[12]-0.06816241358024049*(fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]+0.0545299308641924*(fMquad[8]+fMquad[7])+0.08724788938270785*fMquad[6]+0.0545299308641924*(fMquad[5]+fMquad[4])+0.08724788938270785*fMquad[3]-0.1090598617283848*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[14] = 0.02130075424382515*(fMquad[80]+fMquad[79])-0.0426015084876503*fMquad[78]+0.02130075424382515*(fMquad[77]+fMquad[76])-0.0426015084876503*fMquad[75]+0.03408120679012024*(fMquad[74]+fMquad[73])-0.06816241358024047*fMquad[72]+0.02130075424382515*(fMquad[71]+fMquad[70])-0.0426015084876503*fMquad[69]+0.02130075424382515*(fMquad[68]+fMquad[67])-0.0426015084876503*fMquad[66]+0.03408120679012024*(fMquad[65]+fMquad[64])-0.06816241358024047*fMquad[63]+0.03408120679012024*(fMquad[62]+fMquad[61])-0.06816241358024047*fMquad[60]+0.03408120679012024*(fMquad[59]+fMquad[58])-0.06816241358024047*fMquad[57]+0.0545299308641924*(fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382515*(fMquad[53]+fMquad[52])-0.0426015084876503*fMquad[51]+0.02130075424382515*(fMquad[50]+fMquad[49])-0.0426015084876503*fMquad[48]+0.03408120679012024*(fMquad[47]+fMquad[46])-0.06816241358024047*fMquad[45]+0.02130075424382515*(fMquad[44]+fMquad[43])-0.0426015084876503*fMquad[42]+0.02130075424382515*(fMquad[41]+fMquad[40])-0.0426015084876503*fMquad[39]+0.03408120679012024*(fMquad[38]+fMquad[37])-0.06816241358024047*fMquad[36]+0.03408120679012024*(fMquad[35]+fMquad[34])-0.06816241358024047*fMquad[33]+0.03408120679012024*(fMquad[32]+fMquad[31])-0.06816241358024047*fMquad[30]+0.0545299308641924*(fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012024*(fMquad[26]+fMquad[25])-0.06816241358024047*fMquad[24]+0.03408120679012024*(fMquad[23]+fMquad[22])-0.06816241358024047*fMquad[21]+0.0545299308641924*(fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]+0.03408120679012024*(fMquad[17]+fMquad[16])-0.06816241358024047*fMquad[15]+0.03408120679012024*(fMquad[14]+fMquad[13])-0.06816241358024047*fMquad[12]+0.0545299308641924*(fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]+0.0545299308641924*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]+0.0545299308641924*(fMquad[5]+fMquad[4])-0.1090598617283848*fMquad[3]+0.08724788938270785*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 

}
