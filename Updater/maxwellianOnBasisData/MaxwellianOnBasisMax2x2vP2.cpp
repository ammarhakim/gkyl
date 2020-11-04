#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss2x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

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

  if (vtSqOrd[0] <= 0.0)
    fMFacOrd[0] = 0;
  else
    fMFacOrd[0] = (0.1591549430918953*(0.5*den[0]-0.5590169943749475*(den[5]+den[4])))/vtSqOrd[0]; 
  if (vtSqOrd[1] <= 0.0)
    fMFacOrd[1] = 0;
  else
    fMFacOrd[1] = (0.1591549430918953*(0.4472135954999581*den[5]-0.5590169943749475*den[4]-0.6708203932499369*den[2]+0.5*den[0]))/vtSqOrd[1]; 
  if (vtSqOrd[2] <= 0.0)
    fMFacOrd[2] = 0;
  else
    fMFacOrd[2] = (0.1591549430918953*(0.4472135954999581*den[5]-0.5590169943749475*den[4]+0.6708203932499369*den[2]+0.5*den[0]))/vtSqOrd[2]; 
  if (vtSqOrd[3] <= 0.0)
    fMFacOrd[3] = 0;
  else
    fMFacOrd[3] = (0.1591549430918953*((-0.5590169943749475*den[5])+0.4472135954999581*den[4]-0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[3]; 
  if (vtSqOrd[4] <= 0.0)
    fMFacOrd[4] = 0;
  else
    fMFacOrd[4] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])+0.9*den[3]-0.6708203932499369*(den[2]+den[1])+0.5*den[0]))/vtSqOrd[4]; 
  if (vtSqOrd[5] <= 0.0)
    fMFacOrd[5] = 0;
  else
    fMFacOrd[5] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])-0.9*den[3]+0.6708203932499369*den[2]-0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[5]; 
  if (vtSqOrd[6] <= 0.0)
    fMFacOrd[6] = 0;
  else
    fMFacOrd[6] = (0.1591549430918953*((-0.5590169943749475*den[5])+0.4472135954999581*den[4]+0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[6]; 
  if (vtSqOrd[7] <= 0.0)
    fMFacOrd[7] = 0;
  else
    fMFacOrd[7] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])-0.9*den[3]-0.6708203932499369*den[2]+0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[7]; 
  if (vtSqOrd[8] <= 0.0)
    fMFacOrd[8] = 0;
  else
    fMFacOrd[8] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])+0.9*den[3]+0.6708203932499369*(den[2]+den[1])+0.5*den[0]))/vtSqOrd[8]; 

}

void MaxwellianOnBasisGauss2x2vMaxUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

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

  if (vtSqOrd[0] <= 0.0)
    fMFacOrd[0] = 0;
  else
    fMFacOrd[0] = (0.1591549430918953*(0.5*den[0]-0.5590169943749475*(den[5]+den[4])))/vtSqOrd[0]; 
  if (vtSqOrd[1] <= 0.0)
    fMFacOrd[1] = 0;
  else
    fMFacOrd[1] = (0.1591549430918953*(0.4472135954999581*den[5]-0.5590169943749475*den[4]-0.6708203932499369*den[2]+0.5*den[0]))/vtSqOrd[1]; 
  if (vtSqOrd[2] <= 0.0)
    fMFacOrd[2] = 0;
  else
    fMFacOrd[2] = (0.1591549430918953*(0.4472135954999581*den[5]-0.5590169943749475*den[4]+0.6708203932499369*den[2]+0.5*den[0]))/vtSqOrd[2]; 
  if (vtSqOrd[3] <= 0.0)
    fMFacOrd[3] = 0;
  else
    fMFacOrd[3] = (0.1591549430918953*((-0.5590169943749475*den[5])+0.4472135954999581*den[4]-0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[3]; 
  if (vtSqOrd[4] <= 0.0)
    fMFacOrd[4] = 0;
  else
    fMFacOrd[4] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])+0.9*den[3]-0.6708203932499369*(den[2]+den[1])+0.5*den[0]))/vtSqOrd[4]; 
  if (vtSqOrd[5] <= 0.0)
    fMFacOrd[5] = 0;
  else
    fMFacOrd[5] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])-0.9*den[3]+0.6708203932499369*den[2]-0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[5]; 
  if (vtSqOrd[6] <= 0.0)
    fMFacOrd[6] = 0;
  else
    fMFacOrd[6] = (0.1591549430918953*((-0.5590169943749475*den[5])+0.4472135954999581*den[4]+0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[6]; 
  if (vtSqOrd[7] <= 0.0)
    fMFacOrd[7] = 0;
  else
    fMFacOrd[7] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])-0.9*den[3]-0.6708203932499369*den[2]+0.6708203932499369*den[1]+0.5*den[0]))/vtSqOrd[7]; 
  if (vtSqOrd[8] <= 0.0)
    fMFacOrd[8] = 0;
  else
    fMFacOrd[8] = (0.1591549430918953*(0.4472135954999581*(den[5]+den[4])+0.9*den[3]+0.6708203932499369*(den[2]+den[1])+0.5*den[0]))/vtSqOrd[8]; 

}

void MaxwellianOnBasisGauss2x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[81];
  fMquad[0] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[9],2.0)+std::pow(wc[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[1] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[2] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[3] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[9],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[4] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[5] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[6] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[9],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[7] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[8] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[9])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[9] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[10],2.0)+std::pow(wc[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[10] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[11] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[12] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[10],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[13] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[14] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[15] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[10],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[16] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[17] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[10])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[18] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[11],2.0)+std::pow(wc[2]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
  fMquad[19] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
  fMquad[20] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
  fMquad[21] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[11],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[22] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[23] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[24] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[11],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[25] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[26] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[11])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[27] = fMFacOrd[3]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[12],2.0)+std::pow(wc[2]-1.0*flowUOrd[3],2.0)))/vtSqOrd[3]); 
  fMquad[28] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[3],2.0)))/vtSqOrd[3]); 
  fMquad[29] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[3],2.0)))/vtSqOrd[3]); 
  fMquad[30] = fMFacOrd[3]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[12],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[31] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[32] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[33] = fMFacOrd[3]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[12],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[34] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[35] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[12])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[36] = fMFacOrd[4]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[13],2.0)+std::pow(wc[2]-1.0*flowUOrd[4],2.0)))/vtSqOrd[4]); 
  fMquad[37] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[4],2.0)))/vtSqOrd[4]); 
  fMquad[38] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[4],2.0)))/vtSqOrd[4]); 
  fMquad[39] = fMFacOrd[4]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[13],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[40] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[41] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[42] = fMFacOrd[4]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[13],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[43] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[44] = fMFacOrd[4]*exp(-(0.5*(std::pow((-1.0*flowUOrd[13])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[4]); 
  fMquad[45] = fMFacOrd[5]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[14],2.0)+std::pow(wc[2]-1.0*flowUOrd[5],2.0)))/vtSqOrd[5]); 
  fMquad[46] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[5],2.0)))/vtSqOrd[5]); 
  fMquad[47] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[5],2.0)))/vtSqOrd[5]); 
  fMquad[48] = fMFacOrd[5]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[14],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[49] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[50] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[51] = fMFacOrd[5]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[14],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[52] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[53] = fMFacOrd[5]*exp(-(0.5*(std::pow((-1.0*flowUOrd[14])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[5]); 
  fMquad[54] = fMFacOrd[6]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[15],2.0)+std::pow(wc[2]-1.0*flowUOrd[6],2.0)))/vtSqOrd[6]); 
  fMquad[55] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[6],2.0)))/vtSqOrd[6]); 
  fMquad[56] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[6],2.0)))/vtSqOrd[6]); 
  fMquad[57] = fMFacOrd[6]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[15],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[58] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[59] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[60] = fMFacOrd[6]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[15],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[61] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[62] = fMFacOrd[6]*exp(-(0.5*(std::pow((-1.0*flowUOrd[15])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[6]); 
  fMquad[63] = fMFacOrd[7]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[16],2.0)+std::pow(wc[2]-1.0*flowUOrd[7],2.0)))/vtSqOrd[7]); 
  fMquad[64] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[7],2.0)))/vtSqOrd[7]); 
  fMquad[65] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[7],2.0)))/vtSqOrd[7]); 
  fMquad[66] = fMFacOrd[7]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[16],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[67] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[68] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[69] = fMFacOrd[7]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[16],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[70] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[71] = fMFacOrd[7]*exp(-(0.5*(std::pow((-1.0*flowUOrd[16])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[7]); 
  fMquad[72] = fMFacOrd[8]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[17],2.0)+std::pow(wc[2]-1.0*flowUOrd[8],2.0)))/vtSqOrd[8]); 
  fMquad[73] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[8],2.0)))/vtSqOrd[8]); 
  fMquad[74] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[8],2.0)))/vtSqOrd[8]); 
  fMquad[75] = fMFacOrd[8]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[17],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  fMquad[76] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  fMquad[77] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  fMquad[78] = fMFacOrd[8]*exp(-(0.5*(std::pow(wc[3]-1.0*flowUOrd[17],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  fMquad[79] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]-0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 
  fMquad[80] = fMFacOrd[8]*exp(-(0.5*(std::pow((-1.0*flowUOrd[17])+wc[3]+0.3872983346207417*dxv[3],2.0)+std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0)))/vtSqOrd[8]); 

  fMOut[0] = 0.0238149672306051*(fMquad[80]+fMquad[79])+0.03810394756896815*fMquad[78]+0.0238149672306051*(fMquad[77]+fMquad[76])+0.03810394756896815*(fMquad[75]+fMquad[74]+fMquad[73])+0.06096631611034903*fMquad[72]+0.0238149672306051*(fMquad[71]+fMquad[70])+0.03810394756896815*fMquad[69]+0.0238149672306051*(fMquad[68]+fMquad[67])+0.03810394756896815*(fMquad[66]+fMquad[65]+fMquad[64])+0.06096631611034903*fMquad[63]+0.03810394756896815*(fMquad[62]+fMquad[61])+0.06096631611034903*fMquad[60]+0.03810394756896815*(fMquad[59]+fMquad[58])+0.06096631611034903*(fMquad[57]+fMquad[56]+fMquad[55])+0.09754610577655842*fMquad[54]+0.0238149672306051*(fMquad[53]+fMquad[52])+0.03810394756896815*fMquad[51]+0.0238149672306051*(fMquad[50]+fMquad[49])+0.03810394756896815*(fMquad[48]+fMquad[47]+fMquad[46])+0.06096631611034903*fMquad[45]+0.0238149672306051*(fMquad[44]+fMquad[43])+0.03810394756896815*fMquad[42]+0.0238149672306051*(fMquad[41]+fMquad[40])+0.03810394756896815*(fMquad[39]+fMquad[38]+fMquad[37])+0.06096631611034903*fMquad[36]+0.03810394756896815*(fMquad[35]+fMquad[34])+0.06096631611034903*fMquad[33]+0.03810394756896815*(fMquad[32]+fMquad[31])+0.06096631611034903*(fMquad[30]+fMquad[29]+fMquad[28])+0.09754610577655842*fMquad[27]+0.03810394756896815*(fMquad[26]+fMquad[25])+0.06096631611034903*fMquad[24]+0.03810394756896815*(fMquad[23]+fMquad[22])+0.06096631611034903*(fMquad[21]+fMquad[20]+fMquad[19])+0.09754610577655842*fMquad[18]+0.03810394756896815*(fMquad[17]+fMquad[16])+0.06096631611034903*fMquad[15]+0.03810394756896815*(fMquad[14]+fMquad[13])+0.06096631611034903*(fMquad[12]+fMquad[11]+fMquad[10])+0.09754610577655842*fMquad[9]+0.06096631611034904*(fMquad[8]+fMquad[7])+0.09754610577655844*fMquad[6]+0.06096631611034904*(fMquad[5]+fMquad[4])+0.09754610577655844*fMquad[3]+0.09754610577655845*(fMquad[2]+fMquad[1])+0.1560737692424935*fMquad[0]; 
  fMOut[1] = 0.03195113136573775*(fMquad[80]+fMquad[79])+0.05112181018518038*fMquad[78]+0.03195113136573775*(fMquad[77]+fMquad[76])+0.05112181018518038*(fMquad[75]+fMquad[74]+fMquad[73])+0.0817948962962886*fMquad[72]+0.03195113136573775*(fMquad[71]+fMquad[70])+0.05112181018518038*fMquad[69]+0.03195113136573775*(fMquad[68]+fMquad[67])+0.05112181018518038*(fMquad[66]+fMquad[65]+fMquad[64])+0.0817948962962886*fMquad[63]+0.05112181018518038*(fMquad[62]+fMquad[61])+0.0817948962962886*fMquad[60]+0.05112181018518038*(fMquad[59]+fMquad[58])+0.0817948962962886*(fMquad[57]+fMquad[56]+fMquad[55])+0.1308718340740617*fMquad[54]-0.03195113136573775*(fMquad[53]+fMquad[52])-0.05112181018518038*fMquad[51]-0.03195113136573775*(fMquad[50]+fMquad[49])-0.05112181018518038*(fMquad[48]+fMquad[47]+fMquad[46])-0.0817948962962886*fMquad[45]-0.03195113136573775*(fMquad[44]+fMquad[43])-0.05112181018518038*fMquad[42]-0.03195113136573775*(fMquad[41]+fMquad[40])-0.05112181018518038*(fMquad[39]+fMquad[38]+fMquad[37])-0.0817948962962886*fMquad[36]-0.05112181018518038*(fMquad[35]+fMquad[34])-0.0817948962962886*fMquad[33]-0.05112181018518038*(fMquad[32]+fMquad[31])-0.0817948962962886*(fMquad[30]+fMquad[29]+fMquad[28])-0.1308718340740617*fMquad[27]; 
  fMOut[2] = 0.03195113136573775*(fMquad[80]+fMquad[79])+0.05112181018518038*fMquad[78]+0.03195113136573775*(fMquad[77]+fMquad[76])+0.05112181018518038*(fMquad[75]+fMquad[74]+fMquad[73])+0.0817948962962886*fMquad[72]-0.03195113136573775*(fMquad[71]+fMquad[70])-0.05112181018518038*fMquad[69]-0.03195113136573775*(fMquad[68]+fMquad[67])-0.05112181018518038*(fMquad[66]+fMquad[65]+fMquad[64])-0.0817948962962886*fMquad[63]+0.03195113136573775*(fMquad[53]+fMquad[52])+0.05112181018518038*fMquad[51]+0.03195113136573775*(fMquad[50]+fMquad[49])+0.05112181018518038*(fMquad[48]+fMquad[47]+fMquad[46])+0.0817948962962886*fMquad[45]-0.03195113136573775*(fMquad[44]+fMquad[43])-0.05112181018518038*fMquad[42]-0.03195113136573775*(fMquad[41]+fMquad[40])-0.05112181018518038*(fMquad[39]+fMquad[38]+fMquad[37])-0.0817948962962886*fMquad[36]+0.05112181018518038*(fMquad[26]+fMquad[25])+0.0817948962962886*fMquad[24]+0.05112181018518038*(fMquad[23]+fMquad[22])+0.0817948962962886*(fMquad[21]+fMquad[20]+fMquad[19])+0.1308718340740617*fMquad[18]-0.05112181018518038*(fMquad[17]+fMquad[16])-0.0817948962962886*fMquad[15]-0.05112181018518038*(fMquad[14]+fMquad[13])-0.0817948962962886*(fMquad[12]+fMquad[11]+fMquad[10])-0.1308718340740617*fMquad[9]; 
  fMOut[3] = 0.03195113136573775*(fMquad[80]+fMquad[79])+0.05112181018518038*fMquad[78]-0.03195113136573775*(fMquad[77]+fMquad[76])-0.05112181018518038*fMquad[75]+0.03195113136573775*(fMquad[71]+fMquad[70])+0.05112181018518038*fMquad[69]-0.03195113136573775*(fMquad[68]+fMquad[67])-0.05112181018518038*fMquad[66]+0.05112181018518038*(fMquad[62]+fMquad[61])+0.0817948962962886*fMquad[60]-0.05112181018518038*(fMquad[59]+fMquad[58])-0.0817948962962886*fMquad[57]+0.03195113136573775*(fMquad[53]+fMquad[52])+0.05112181018518038*fMquad[51]-0.03195113136573775*(fMquad[50]+fMquad[49])-0.05112181018518038*fMquad[48]+0.03195113136573775*(fMquad[44]+fMquad[43])+0.05112181018518038*fMquad[42]-0.03195113136573775*(fMquad[41]+fMquad[40])-0.05112181018518038*fMquad[39]+0.05112181018518038*(fMquad[35]+fMquad[34])+0.0817948962962886*fMquad[33]-0.05112181018518038*(fMquad[32]+fMquad[31])-0.0817948962962886*fMquad[30]+0.05112181018518038*(fMquad[26]+fMquad[25])+0.0817948962962886*fMquad[24]-0.05112181018518038*(fMquad[23]+fMquad[22])-0.0817948962962886*fMquad[21]+0.05112181018518038*(fMquad[17]+fMquad[16])+0.0817948962962886*fMquad[15]-0.05112181018518038*(fMquad[14]+fMquad[13])-0.0817948962962886*fMquad[12]+0.0817948962962886*(fMquad[8]+fMquad[7])+0.1308718340740617*fMquad[6]-0.0817948962962886*(fMquad[5]+fMquad[4])-0.1308718340740617*fMquad[3]; 
  fMOut[4] = 0.03195113136573775*fMquad[80]-0.03195113136573775*fMquad[79]+0.03195113136573775*fMquad[77]-0.03195113136573775*fMquad[76]+0.05112181018518038*fMquad[74]-0.05112181018518038*fMquad[73]+0.03195113136573775*fMquad[71]-0.03195113136573775*fMquad[70]+0.03195113136573775*fMquad[68]-0.03195113136573775*fMquad[67]+0.05112181018518038*fMquad[65]-0.05112181018518038*fMquad[64]+0.05112181018518038*fMquad[62]-0.05112181018518038*fMquad[61]+0.05112181018518038*fMquad[59]-0.05112181018518038*fMquad[58]+0.0817948962962886*fMquad[56]-0.0817948962962886*fMquad[55]+0.03195113136573775*fMquad[53]-0.03195113136573775*fMquad[52]+0.03195113136573775*fMquad[50]-0.03195113136573775*fMquad[49]+0.05112181018518038*fMquad[47]-0.05112181018518038*fMquad[46]+0.03195113136573775*fMquad[44]-0.03195113136573775*fMquad[43]+0.03195113136573775*fMquad[41]-0.03195113136573775*fMquad[40]+0.05112181018518038*fMquad[38]-0.05112181018518038*fMquad[37]+0.05112181018518038*fMquad[35]-0.05112181018518038*fMquad[34]+0.05112181018518038*fMquad[32]-0.05112181018518038*fMquad[31]+0.0817948962962886*fMquad[29]-0.0817948962962886*fMquad[28]+0.05112181018518038*fMquad[26]-0.05112181018518038*fMquad[25]+0.05112181018518038*fMquad[23]-0.05112181018518038*fMquad[22]+0.0817948962962886*fMquad[20]-0.0817948962962886*fMquad[19]+0.05112181018518038*fMquad[17]-0.05112181018518038*fMquad[16]+0.05112181018518038*fMquad[14]-0.05112181018518038*fMquad[13]+0.0817948962962886*fMquad[11]-0.0817948962962886*fMquad[10]+0.0817948962962886*fMquad[8]-0.0817948962962886*fMquad[7]+0.0817948962962886*fMquad[5]-0.0817948962962886*fMquad[4]+0.1308718340740617*fMquad[2]-0.1308718340740617*fMquad[1]; 
  fMOut[5] = 0.04286694101508918*(fMquad[80]+fMquad[79])+0.06858710562414268*fMquad[78]+0.04286694101508918*(fMquad[77]+fMquad[76])+0.06858710562414268*(fMquad[75]+fMquad[74]+fMquad[73])+0.1097393689986283*fMquad[72]-0.04286694101508918*(fMquad[71]+fMquad[70])-0.06858710562414268*fMquad[69]-0.04286694101508918*(fMquad[68]+fMquad[67])-0.06858710562414268*(fMquad[66]+fMquad[65]+fMquad[64])-0.1097393689986283*fMquad[63]-0.04286694101508918*(fMquad[53]+fMquad[52])-0.06858710562414268*fMquad[51]-0.04286694101508918*(fMquad[50]+fMquad[49])-0.06858710562414268*(fMquad[48]+fMquad[47]+fMquad[46])-0.1097393689986283*fMquad[45]+0.04286694101508918*(fMquad[44]+fMquad[43])+0.06858710562414268*fMquad[42]+0.04286694101508918*(fMquad[41]+fMquad[40])+0.06858710562414268*(fMquad[39]+fMquad[38]+fMquad[37])+0.1097393689986283*fMquad[36]; 
  fMOut[6] = 0.04286694101508918*(fMquad[80]+fMquad[79])+0.06858710562414268*fMquad[78]-0.04286694101508918*(fMquad[77]+fMquad[76])-0.06858710562414268*fMquad[75]+0.04286694101508918*(fMquad[71]+fMquad[70])+0.06858710562414268*fMquad[69]-0.04286694101508918*(fMquad[68]+fMquad[67])-0.06858710562414268*fMquad[66]+0.06858710562414268*(fMquad[62]+fMquad[61])+0.1097393689986283*fMquad[60]-0.06858710562414268*(fMquad[59]+fMquad[58])-0.1097393689986283*fMquad[57]-0.04286694101508918*(fMquad[53]+fMquad[52])-0.06858710562414268*fMquad[51]+0.04286694101508918*(fMquad[50]+fMquad[49])+0.06858710562414268*fMquad[48]-0.04286694101508918*(fMquad[44]+fMquad[43])-0.06858710562414268*fMquad[42]+0.04286694101508918*(fMquad[41]+fMquad[40])+0.06858710562414268*fMquad[39]-0.06858710562414268*(fMquad[35]+fMquad[34])-0.1097393689986283*fMquad[33]+0.06858710562414268*(fMquad[32]+fMquad[31])+0.1097393689986283*fMquad[30]; 
  fMOut[7] = 0.04286694101508918*(fMquad[80]+fMquad[79])+0.06858710562414268*fMquad[78]-0.04286694101508918*(fMquad[77]+fMquad[76])-0.06858710562414268*fMquad[75]-0.04286694101508918*(fMquad[71]+fMquad[70])-0.06858710562414268*fMquad[69]+0.04286694101508918*(fMquad[68]+fMquad[67])+0.06858710562414268*fMquad[66]+0.04286694101508918*(fMquad[53]+fMquad[52])+0.06858710562414268*fMquad[51]-0.04286694101508918*(fMquad[50]+fMquad[49])-0.06858710562414268*fMquad[48]-0.04286694101508918*(fMquad[44]+fMquad[43])-0.06858710562414268*fMquad[42]+0.04286694101508918*(fMquad[41]+fMquad[40])+0.06858710562414268*(fMquad[39]+fMquad[26]+fMquad[25])+0.1097393689986283*fMquad[24]-0.06858710562414268*(fMquad[23]+fMquad[22])-0.1097393689986283*fMquad[21]-0.06858710562414268*(fMquad[17]+fMquad[16])-0.1097393689986283*fMquad[15]+0.06858710562414268*(fMquad[14]+fMquad[13])+0.1097393689986283*fMquad[12]; 
  fMOut[8] = 0.04286694101508918*fMquad[80]-0.04286694101508918*fMquad[79]+0.04286694101508918*fMquad[77]-0.04286694101508918*fMquad[76]+0.06858710562414268*fMquad[74]-0.06858710562414268*fMquad[73]+0.04286694101508918*fMquad[71]-0.04286694101508918*fMquad[70]+0.04286694101508918*fMquad[68]-0.04286694101508918*fMquad[67]+0.06858710562414268*fMquad[65]-0.06858710562414268*fMquad[64]+0.06858710562414268*fMquad[62]-0.06858710562414268*fMquad[61]+0.06858710562414268*fMquad[59]-0.06858710562414268*fMquad[58]+0.1097393689986283*fMquad[56]-0.1097393689986283*fMquad[55]-0.04286694101508918*fMquad[53]+0.04286694101508918*fMquad[52]-0.04286694101508918*fMquad[50]+0.04286694101508918*fMquad[49]-0.06858710562414268*fMquad[47]+0.06858710562414268*fMquad[46]-0.04286694101508918*fMquad[44]+0.04286694101508918*fMquad[43]-0.04286694101508918*fMquad[41]+0.04286694101508918*fMquad[40]-0.06858710562414268*fMquad[38]+0.06858710562414268*fMquad[37]-0.06858710562414268*fMquad[35]+0.06858710562414268*fMquad[34]-0.06858710562414268*fMquad[32]+0.06858710562414268*fMquad[31]-0.1097393689986283*fMquad[29]+0.1097393689986283*fMquad[28]; 
  fMOut[9] = 0.04286694101508918*fMquad[80]-0.04286694101508918*fMquad[79]+0.04286694101508918*fMquad[77]-0.04286694101508918*fMquad[76]+0.06858710562414268*fMquad[74]-0.06858710562414268*fMquad[73]-0.04286694101508918*fMquad[71]+0.04286694101508918*fMquad[70]-0.04286694101508918*fMquad[68]+0.04286694101508918*fMquad[67]-0.06858710562414268*fMquad[65]+0.06858710562414268*fMquad[64]+0.04286694101508918*fMquad[53]-0.04286694101508918*fMquad[52]+0.04286694101508918*fMquad[50]-0.04286694101508918*fMquad[49]+0.06858710562414268*fMquad[47]-0.06858710562414268*fMquad[46]-0.04286694101508918*fMquad[44]+0.04286694101508918*fMquad[43]-0.04286694101508918*fMquad[41]+0.04286694101508918*fMquad[40]-0.06858710562414268*fMquad[38]+0.06858710562414268*(fMquad[37]+fMquad[26])-0.06858710562414268*fMquad[25]+0.06858710562414268*fMquad[23]-0.06858710562414268*fMquad[22]+0.1097393689986283*fMquad[20]-0.1097393689986283*fMquad[19]-0.06858710562414268*fMquad[17]+0.06858710562414268*fMquad[16]-0.06858710562414268*fMquad[14]+0.06858710562414268*fMquad[13]-0.1097393689986283*fMquad[11]+0.1097393689986283*fMquad[10]; 
  fMOut[10] = 0.04286694101508918*fMquad[80]-0.04286694101508918*(fMquad[79]+fMquad[77])+0.04286694101508918*(fMquad[76]+fMquad[71])-0.04286694101508918*(fMquad[70]+fMquad[68])+0.04286694101508918*fMquad[67]+0.06858710562414268*fMquad[62]-0.06858710562414268*(fMquad[61]+fMquad[59])+0.06858710562414268*fMquad[58]+0.04286694101508918*fMquad[53]-0.04286694101508918*(fMquad[52]+fMquad[50])+0.04286694101508918*(fMquad[49]+fMquad[44])-0.04286694101508918*(fMquad[43]+fMquad[41])+0.04286694101508918*fMquad[40]+0.06858710562414268*fMquad[35]-0.06858710562414268*(fMquad[34]+fMquad[32])+0.06858710562414268*(fMquad[31]+fMquad[26])-0.06858710562414268*(fMquad[25]+fMquad[23])+0.06858710562414268*(fMquad[22]+fMquad[17])-0.06858710562414268*(fMquad[16]+fMquad[14])+0.06858710562414268*fMquad[13]+0.1097393689986283*fMquad[8]-0.1097393689986283*(fMquad[7]+fMquad[5])+0.1097393689986283*fMquad[4]; 
  fMOut[11] = 0.02130075424382517*(fMquad[80]+fMquad[79])+0.03408120679012027*fMquad[78]+0.02130075424382517*(fMquad[77]+fMquad[76])+0.03408120679012027*(fMquad[75]+fMquad[74]+fMquad[73])+0.05452993086419242*fMquad[72]+0.02130075424382517*(fMquad[71]+fMquad[70])+0.03408120679012027*fMquad[69]+0.02130075424382517*(fMquad[68]+fMquad[67])+0.03408120679012027*(fMquad[66]+fMquad[65]+fMquad[64])+0.05452993086419242*fMquad[63]+0.03408120679012027*(fMquad[62]+fMquad[61])+0.05452993086419242*fMquad[60]+0.03408120679012027*(fMquad[59]+fMquad[58])+0.05452993086419242*(fMquad[57]+fMquad[56]+fMquad[55])+0.08724788938270785*fMquad[54]+0.02130075424382517*(fMquad[53]+fMquad[52])+0.03408120679012027*fMquad[51]+0.02130075424382517*(fMquad[50]+fMquad[49])+0.03408120679012027*(fMquad[48]+fMquad[47]+fMquad[46])+0.05452993086419242*fMquad[45]+0.02130075424382517*(fMquad[44]+fMquad[43])+0.03408120679012027*fMquad[42]+0.02130075424382517*(fMquad[41]+fMquad[40])+0.03408120679012027*(fMquad[39]+fMquad[38]+fMquad[37])+0.05452993086419242*fMquad[36]+0.03408120679012027*(fMquad[35]+fMquad[34])+0.05452993086419242*fMquad[33]+0.03408120679012027*(fMquad[32]+fMquad[31])+0.05452993086419242*(fMquad[30]+fMquad[29]+fMquad[28])+0.08724788938270785*fMquad[27]-0.04260150848765032*(fMquad[26]+fMquad[25])-0.0681624135802405*fMquad[24]-0.04260150848765032*(fMquad[23]+fMquad[22])-0.0681624135802405*(fMquad[21]+fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]-0.04260150848765032*(fMquad[17]+fMquad[16])-0.0681624135802405*fMquad[15]-0.04260150848765032*(fMquad[14]+fMquad[13])-0.0681624135802405*(fMquad[12]+fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]-0.06816241358024051*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]-0.06816241358024051*(fMquad[5]+fMquad[4])-0.1090598617283848*fMquad[3]-0.1090598617283848*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[12] = 0.02130075424382517*(fMquad[80]+fMquad[79])+0.03408120679012027*fMquad[78]+0.02130075424382517*(fMquad[77]+fMquad[76])+0.03408120679012027*(fMquad[75]+fMquad[74]+fMquad[73])+0.05452993086419242*fMquad[72]+0.02130075424382517*(fMquad[71]+fMquad[70])+0.03408120679012027*fMquad[69]+0.02130075424382517*(fMquad[68]+fMquad[67])+0.03408120679012027*(fMquad[66]+fMquad[65]+fMquad[64])+0.05452993086419242*fMquad[63]-0.04260150848765032*(fMquad[62]+fMquad[61])-0.0681624135802405*fMquad[60]-0.04260150848765032*(fMquad[59]+fMquad[58])-0.0681624135802405*(fMquad[57]+fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382517*(fMquad[53]+fMquad[52])+0.03408120679012027*fMquad[51]+0.02130075424382517*(fMquad[50]+fMquad[49])+0.03408120679012027*(fMquad[48]+fMquad[47]+fMquad[46])+0.05452993086419242*fMquad[45]+0.02130075424382517*(fMquad[44]+fMquad[43])+0.03408120679012027*fMquad[42]+0.02130075424382517*(fMquad[41]+fMquad[40])+0.03408120679012027*(fMquad[39]+fMquad[38]+fMquad[37])+0.05452993086419242*fMquad[36]-0.04260150848765032*(fMquad[35]+fMquad[34])-0.0681624135802405*fMquad[33]-0.04260150848765032*(fMquad[32]+fMquad[31])-0.0681624135802405*(fMquad[30]+fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012027*(fMquad[26]+fMquad[25])+0.05452993086419242*fMquad[24]+0.03408120679012027*(fMquad[23]+fMquad[22])+0.05452993086419242*(fMquad[21]+fMquad[20]+fMquad[19])+0.08724788938270785*fMquad[18]+0.03408120679012027*(fMquad[17]+fMquad[16])+0.05452993086419242*fMquad[15]+0.03408120679012027*(fMquad[14]+fMquad[13])+0.05452993086419242*(fMquad[12]+fMquad[11]+fMquad[10])+0.08724788938270785*fMquad[9]-0.06816241358024051*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]-0.06816241358024051*(fMquad[5]+fMquad[4])-0.1090598617283848*fMquad[3]-0.1090598617283848*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[13] = 0.02130075424382517*(fMquad[80]+fMquad[79])+0.03408120679012027*fMquad[78]+0.02130075424382517*(fMquad[77]+fMquad[76])+0.03408120679012027*fMquad[75]-0.04260150848765032*(fMquad[74]+fMquad[73])-0.0681624135802405*fMquad[72]+0.02130075424382517*(fMquad[71]+fMquad[70])+0.03408120679012027*fMquad[69]+0.02130075424382517*(fMquad[68]+fMquad[67])+0.03408120679012027*fMquad[66]-0.04260150848765032*(fMquad[65]+fMquad[64])-0.0681624135802405*fMquad[63]+0.03408120679012027*(fMquad[62]+fMquad[61])+0.05452993086419242*fMquad[60]+0.03408120679012027*(fMquad[59]+fMquad[58])+0.05452993086419242*fMquad[57]-0.0681624135802405*(fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382517*(fMquad[53]+fMquad[52])+0.03408120679012027*fMquad[51]+0.02130075424382517*(fMquad[50]+fMquad[49])+0.03408120679012027*fMquad[48]-0.04260150848765032*(fMquad[47]+fMquad[46])-0.0681624135802405*fMquad[45]+0.02130075424382517*(fMquad[44]+fMquad[43])+0.03408120679012027*fMquad[42]+0.02130075424382517*(fMquad[41]+fMquad[40])+0.03408120679012027*fMquad[39]-0.04260150848765032*(fMquad[38]+fMquad[37])-0.0681624135802405*fMquad[36]+0.03408120679012027*(fMquad[35]+fMquad[34])+0.05452993086419242*fMquad[33]+0.03408120679012027*(fMquad[32]+fMquad[31])+0.05452993086419242*fMquad[30]-0.0681624135802405*(fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012027*(fMquad[26]+fMquad[25])+0.05452993086419242*fMquad[24]+0.03408120679012027*(fMquad[23]+fMquad[22])+0.05452993086419242*fMquad[21]-0.0681624135802405*(fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]+0.03408120679012027*(fMquad[17]+fMquad[16])+0.05452993086419242*fMquad[15]+0.03408120679012027*(fMquad[14]+fMquad[13])+0.05452993086419242*fMquad[12]-0.0681624135802405*(fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]+0.05452993086419243*(fMquad[8]+fMquad[7])+0.08724788938270786*fMquad[6]+0.05452993086419243*(fMquad[5]+fMquad[4])+0.08724788938270786*fMquad[3]-0.1090598617283848*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[14] = 0.02130075424382517*(fMquad[80]+fMquad[79])-0.04260150848765032*fMquad[78]+0.02130075424382517*(fMquad[77]+fMquad[76])-0.04260150848765032*fMquad[75]+0.03408120679012027*(fMquad[74]+fMquad[73])-0.0681624135802405*fMquad[72]+0.02130075424382517*(fMquad[71]+fMquad[70])-0.04260150848765032*fMquad[69]+0.02130075424382517*(fMquad[68]+fMquad[67])-0.04260150848765032*fMquad[66]+0.03408120679012027*(fMquad[65]+fMquad[64])-0.0681624135802405*fMquad[63]+0.03408120679012027*(fMquad[62]+fMquad[61])-0.0681624135802405*fMquad[60]+0.03408120679012027*(fMquad[59]+fMquad[58])-0.0681624135802405*fMquad[57]+0.05452993086419242*(fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382517*(fMquad[53]+fMquad[52])-0.04260150848765032*fMquad[51]+0.02130075424382517*(fMquad[50]+fMquad[49])-0.04260150848765032*fMquad[48]+0.03408120679012027*(fMquad[47]+fMquad[46])-0.0681624135802405*fMquad[45]+0.02130075424382517*(fMquad[44]+fMquad[43])-0.04260150848765032*fMquad[42]+0.02130075424382517*(fMquad[41]+fMquad[40])-0.04260150848765032*fMquad[39]+0.03408120679012027*(fMquad[38]+fMquad[37])-0.0681624135802405*fMquad[36]+0.03408120679012027*(fMquad[35]+fMquad[34])-0.0681624135802405*fMquad[33]+0.03408120679012027*(fMquad[32]+fMquad[31])-0.0681624135802405*fMquad[30]+0.05452993086419242*(fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012027*(fMquad[26]+fMquad[25])-0.0681624135802405*fMquad[24]+0.03408120679012027*(fMquad[23]+fMquad[22])-0.0681624135802405*fMquad[21]+0.05452993086419242*(fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]+0.03408120679012027*(fMquad[17]+fMquad[16])-0.0681624135802405*fMquad[15]+0.03408120679012027*(fMquad[14]+fMquad[13])-0.0681624135802405*fMquad[12]+0.05452993086419242*(fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]+0.05452993086419243*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]+0.05452993086419243*(fMquad[5]+fMquad[4])-0.1090598617283848*fMquad[3]+0.08724788938270787*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 

}
void GkMaxwellianOnBasisGauss2x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

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

  if (vtSqOrd[0] <= 0.0)
    fMFacOrd[0] = 9.999999999999999e-41;
  else
    fMFacOrd[0] = (bmagOrd[0]*(0.5*den[0]-0.5590169943749475*(den[5]+den[4])))/std::pow(2.506628274631001*sqrt(vtSqOrd[0]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[1] <= 0.0)
    fMFacOrd[1] = 9.999999999999999e-41;
  else
    fMFacOrd[1] = (bmagOrd[1]*(0.4472135954999581*den[5]-0.5590169943749475*den[4]-0.6708203932499369*den[2]+0.5*den[0]))/std::pow(2.506628274631001*sqrt(vtSqOrd[1]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[2] <= 0.0)
    fMFacOrd[2] = 9.999999999999999e-41;
  else
    fMFacOrd[2] = (bmagOrd[2]*(0.4472135954999581*den[5]-0.5590169943749475*den[4]+0.6708203932499369*den[2]+0.5*den[0]))/std::pow(2.506628274631001*sqrt(vtSqOrd[2]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[3] <= 0.0)
    fMFacOrd[3] = 9.999999999999999e-41;
  else
    fMFacOrd[3] = (bmagOrd[3]*((-0.5590169943749475*den[5])+0.4472135954999581*den[4]-0.6708203932499369*den[1]+0.5*den[0]))/std::pow(2.506628274631001*sqrt(vtSqOrd[3]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[4] <= 0.0)
    fMFacOrd[4] = 9.999999999999999e-41;
  else
    fMFacOrd[4] = (bmagOrd[4]*(0.4472135954999581*(den[5]+den[4])+0.9*den[3]-0.6708203932499369*(den[2]+den[1])+0.5*den[0]))/std::pow(2.506628274631001*sqrt(vtSqOrd[4]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[5] <= 0.0)
    fMFacOrd[5] = 9.999999999999999e-41;
  else
    fMFacOrd[5] = (bmagOrd[5]*(0.4472135954999581*(den[5]+den[4])-0.9*den[3]+0.6708203932499369*den[2]-0.6708203932499369*den[1]+0.5*den[0]))/std::pow(2.506628274631001*sqrt(vtSqOrd[5]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[6] <= 0.0)
    fMFacOrd[6] = 9.999999999999999e-41;
  else
    fMFacOrd[6] = (((-0.5590169943749475*den[5])+0.4472135954999581*den[4]+0.6708203932499369*den[1]+0.5*den[0])*bmagOrd[6])/std::pow(2.506628274631001*sqrt(vtSqOrd[6]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[7] <= 0.0)
    fMFacOrd[7] = 9.999999999999999e-41;
  else
    fMFacOrd[7] = ((0.4472135954999581*(den[5]+den[4])-0.9*den[3]-0.6708203932499369*den[2]+0.6708203932499369*den[1]+0.5*den[0])*bmagOrd[7])/std::pow(2.506628274631001*sqrt(vtSqOrd[7]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[8] <= 0.0)
    fMFacOrd[8] = 9.999999999999999e-41;
  else
    fMFacOrd[8] = ((0.4472135954999581*(den[5]+den[4])+0.9*den[3]+0.6708203932499369*(den[2]+den[1])+0.5*den[0])*bmagOrd[8])/std::pow(2.506628274631001*sqrt(vtSqOrd[8]),3.0)+9.999999999999999e-41; 

}

void GkMaxwellianOnBasisGauss2x2vMaxUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

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

  if (vtSqOrd[0] <= 0.0)
    fMFacOrd[0] = 9.999999999999999e-41;
  else
    fMFacOrd[0] = (bmagOrd[0]*(0.5*den[0]-0.5590169943749475*(den[5]+den[4])))/std::pow(2.506628274631001*sqrt(vtSqOrd[0]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[1] <= 0.0)
    fMFacOrd[1] = 9.999999999999999e-41;
  else
    fMFacOrd[1] = (bmagOrd[1]*(0.4472135954999581*den[5]-0.5590169943749475*den[4]-0.6708203932499369*den[2]+0.5*den[0]))/std::pow(2.506628274631001*sqrt(vtSqOrd[1]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[2] <= 0.0)
    fMFacOrd[2] = 9.999999999999999e-41;
  else
    fMFacOrd[2] = (bmagOrd[2]*(0.4472135954999581*den[5]-0.5590169943749475*den[4]+0.6708203932499369*den[2]+0.5*den[0]))/std::pow(2.506628274631001*sqrt(vtSqOrd[2]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[3] <= 0.0)
    fMFacOrd[3] = 9.999999999999999e-41;
  else
    fMFacOrd[3] = (bmagOrd[3]*((-0.5590169943749475*den[5])+0.4472135954999581*den[4]-0.6708203932499369*den[1]+0.5*den[0]))/std::pow(2.506628274631001*sqrt(vtSqOrd[3]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[4] <= 0.0)
    fMFacOrd[4] = 9.999999999999999e-41;
  else
    fMFacOrd[4] = (bmagOrd[4]*(0.4472135954999581*(den[5]+den[4])+0.9*den[3]-0.6708203932499369*(den[2]+den[1])+0.5*den[0]))/std::pow(2.506628274631001*sqrt(vtSqOrd[4]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[5] <= 0.0)
    fMFacOrd[5] = 9.999999999999999e-41;
  else
    fMFacOrd[5] = (bmagOrd[5]*(0.4472135954999581*(den[5]+den[4])-0.9*den[3]+0.6708203932499369*den[2]-0.6708203932499369*den[1]+0.5*den[0]))/std::pow(2.506628274631001*sqrt(vtSqOrd[5]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[6] <= 0.0)
    fMFacOrd[6] = 9.999999999999999e-41;
  else
    fMFacOrd[6] = (((-0.5590169943749475*den[5])+0.4472135954999581*den[4]+0.6708203932499369*den[1]+0.5*den[0])*bmagOrd[6])/std::pow(2.506628274631001*sqrt(vtSqOrd[6]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[7] <= 0.0)
    fMFacOrd[7] = 9.999999999999999e-41;
  else
    fMFacOrd[7] = ((0.4472135954999581*(den[5]+den[4])-0.9*den[3]-0.6708203932499369*den[2]+0.6708203932499369*den[1]+0.5*den[0])*bmagOrd[7])/std::pow(2.506628274631001*sqrt(vtSqOrd[7]),3.0)+9.999999999999999e-41; 
  if (vtSqOrd[8] <= 0.0)
    fMFacOrd[8] = 9.999999999999999e-41;
  else
    fMFacOrd[8] = ((0.4472135954999581*(den[5]+den[4])+0.9*den[3]+0.6708203932499369*(den[2]+den[1])+0.5*den[0])*bmagOrd[8])/std::pow(2.506628274631001*sqrt(vtSqOrd[8]),3.0)+9.999999999999999e-41; 

}

void GkMaxwellianOnBasisGauss2x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[81];
  fMquad[0] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*wc[3])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[1] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[2] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[3] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*wc[3])/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[4] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[5] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[6] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*wc[3])/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[7] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[8] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[9] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*wc[3])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1]); 
  fMquad[10] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1]); 
  fMquad[11] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1]); 
  fMquad[12] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*wc[3])/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1]); 
  fMquad[13] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1]); 
  fMquad[14] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1]); 
  fMquad[15] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*wc[3])/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1]); 
  fMquad[16] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1]); 
  fMquad[17] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]+0.3872983346207417*dxv[2]-1.0*flowUOrd[1],2.0))/vtSqOrd[1]); 
  fMquad[18] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*wc[3])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2],2.0))/vtSqOrd[2]); 
  fMquad[19] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2],2.0))/vtSqOrd[2]); 
  fMquad[20] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2],2.0))/vtSqOrd[2]); 
  fMquad[21] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*wc[3])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[2]); 
  fMquad[22] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[2]); 
  fMquad[23] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[2]); 
  fMquad[24] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*wc[3])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[2]); 
  fMquad[25] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[2]); 
  fMquad[26] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[2]); 
  fMquad[27] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*wc[3])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[3],2.0))/vtSqOrd[3]); 
  fMquad[28] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[3],2.0))/vtSqOrd[3]); 
  fMquad[29] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[3],2.0))/vtSqOrd[3]); 
  fMquad[30] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*wc[3])/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[3]); 
  fMquad[31] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[3]); 
  fMquad[32] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[3]); 
  fMquad[33] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*wc[3])/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[3]); 
  fMquad[34] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*(wc[3]-0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[3]); 
  fMquad[35] = fMFacOrd[3]*exp(((-(1.0*bmagOrd[3]*(wc[3]+0.3872983346207417*dxv[3]))/m_)-0.5*std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[3]); 
  fMquad[36] = fMFacOrd[4]*exp(((-(1.0*wc[3]*bmagOrd[4])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[4],2.0))/vtSqOrd[4]); 
  fMquad[37] = fMFacOrd[4]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[4])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[4],2.0))/vtSqOrd[4]); 
  fMquad[38] = fMFacOrd[4]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[4])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[4],2.0))/vtSqOrd[4]); 
  fMquad[39] = fMFacOrd[4]*exp(((-(1.0*wc[3]*bmagOrd[4])/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[4]); 
  fMquad[40] = fMFacOrd[4]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[4])/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[4]); 
  fMquad[41] = fMFacOrd[4]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[4])/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[4]); 
  fMquad[42] = fMFacOrd[4]*exp(((-(1.0*wc[3]*bmagOrd[4])/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[4]); 
  fMquad[43] = fMFacOrd[4]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[4])/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[4]); 
  fMquad[44] = fMFacOrd[4]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[4])/m_)-0.5*std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[4]); 
  fMquad[45] = fMFacOrd[5]*exp(((-(1.0*wc[3]*bmagOrd[5])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[5],2.0))/vtSqOrd[5]); 
  fMquad[46] = fMFacOrd[5]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[5])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[5],2.0))/vtSqOrd[5]); 
  fMquad[47] = fMFacOrd[5]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[5])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[5],2.0))/vtSqOrd[5]); 
  fMquad[48] = fMFacOrd[5]*exp(((-(1.0*wc[3]*bmagOrd[5])/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[5]); 
  fMquad[49] = fMFacOrd[5]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[5])/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[5]); 
  fMquad[50] = fMFacOrd[5]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[5])/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[5]); 
  fMquad[51] = fMFacOrd[5]*exp(((-(1.0*wc[3]*bmagOrd[5])/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[5]); 
  fMquad[52] = fMFacOrd[5]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[5])/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[5]); 
  fMquad[53] = fMFacOrd[5]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[5])/m_)-0.5*std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[5]); 
  fMquad[54] = fMFacOrd[6]*exp(((-(1.0*wc[3]*bmagOrd[6])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[6],2.0))/vtSqOrd[6]); 
  fMquad[55] = fMFacOrd[6]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[6])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[6],2.0))/vtSqOrd[6]); 
  fMquad[56] = fMFacOrd[6]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[6])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[6],2.0))/vtSqOrd[6]); 
  fMquad[57] = fMFacOrd[6]*exp(((-(1.0*wc[3]*bmagOrd[6])/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[6]); 
  fMquad[58] = fMFacOrd[6]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[6])/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[6]); 
  fMquad[59] = fMFacOrd[6]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[6])/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[6]); 
  fMquad[60] = fMFacOrd[6]*exp(((-(1.0*wc[3]*bmagOrd[6])/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[6]); 
  fMquad[61] = fMFacOrd[6]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[6])/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[6]); 
  fMquad[62] = fMFacOrd[6]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[6])/m_)-0.5*std::pow((-1.0*flowUOrd[6])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[6]); 
  fMquad[63] = fMFacOrd[7]*exp(((-(1.0*wc[3]*bmagOrd[7])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[7],2.0))/vtSqOrd[7]); 
  fMquad[64] = fMFacOrd[7]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[7])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[7],2.0))/vtSqOrd[7]); 
  fMquad[65] = fMFacOrd[7]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[7])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[7],2.0))/vtSqOrd[7]); 
  fMquad[66] = fMFacOrd[7]*exp(((-(1.0*wc[3]*bmagOrd[7])/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[7]); 
  fMquad[67] = fMFacOrd[7]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[7])/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[7]); 
  fMquad[68] = fMFacOrd[7]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[7])/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[7]); 
  fMquad[69] = fMFacOrd[7]*exp(((-(1.0*wc[3]*bmagOrd[7])/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[7]); 
  fMquad[70] = fMFacOrd[7]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[7])/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[7]); 
  fMquad[71] = fMFacOrd[7]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[7])/m_)-0.5*std::pow((-1.0*flowUOrd[7])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[7]); 
  fMquad[72] = fMFacOrd[8]*exp(((-(1.0*wc[3]*bmagOrd[8])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[8],2.0))/vtSqOrd[8]); 
  fMquad[73] = fMFacOrd[8]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[8])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[8],2.0))/vtSqOrd[8]); 
  fMquad[74] = fMFacOrd[8]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[8])/m_)-0.5*std::pow(wc[2]-1.0*flowUOrd[8],2.0))/vtSqOrd[8]); 
  fMquad[75] = fMFacOrd[8]*exp(((-(1.0*wc[3]*bmagOrd[8])/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[8]); 
  fMquad[76] = fMFacOrd[8]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[8])/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[8]); 
  fMquad[77] = fMFacOrd[8]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[8])/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]-0.3872983346207417*dxv[2],2.0))/vtSqOrd[8]); 
  fMquad[78] = fMFacOrd[8]*exp(((-(1.0*wc[3]*bmagOrd[8])/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[8]); 
  fMquad[79] = fMFacOrd[8]*exp(((-(1.0*(wc[3]-0.3872983346207417*dxv[3])*bmagOrd[8])/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[8]); 
  fMquad[80] = fMFacOrd[8]*exp(((-(1.0*(wc[3]+0.3872983346207417*dxv[3])*bmagOrd[8])/m_)-0.5*std::pow((-1.0*flowUOrd[8])+wc[2]+0.3872983346207417*dxv[2],2.0))/vtSqOrd[8]); 

  fMOut[0] = 0.0238149672306051*(fMquad[80]+fMquad[79])+0.03810394756896815*fMquad[78]+0.0238149672306051*(fMquad[77]+fMquad[76])+0.03810394756896815*(fMquad[75]+fMquad[74]+fMquad[73])+0.06096631611034903*fMquad[72]+0.0238149672306051*(fMquad[71]+fMquad[70])+0.03810394756896815*fMquad[69]+0.0238149672306051*(fMquad[68]+fMquad[67])+0.03810394756896815*(fMquad[66]+fMquad[65]+fMquad[64])+0.06096631611034903*fMquad[63]+0.03810394756896815*(fMquad[62]+fMquad[61])+0.06096631611034903*fMquad[60]+0.03810394756896815*(fMquad[59]+fMquad[58])+0.06096631611034903*(fMquad[57]+fMquad[56]+fMquad[55])+0.09754610577655842*fMquad[54]+0.0238149672306051*(fMquad[53]+fMquad[52])+0.03810394756896815*fMquad[51]+0.0238149672306051*(fMquad[50]+fMquad[49])+0.03810394756896815*(fMquad[48]+fMquad[47]+fMquad[46])+0.06096631611034903*fMquad[45]+0.0238149672306051*(fMquad[44]+fMquad[43])+0.03810394756896815*fMquad[42]+0.0238149672306051*(fMquad[41]+fMquad[40])+0.03810394756896815*(fMquad[39]+fMquad[38]+fMquad[37])+0.06096631611034903*fMquad[36]+0.03810394756896815*(fMquad[35]+fMquad[34])+0.06096631611034903*fMquad[33]+0.03810394756896815*(fMquad[32]+fMquad[31])+0.06096631611034903*(fMquad[30]+fMquad[29]+fMquad[28])+0.09754610577655842*fMquad[27]+0.03810394756896815*(fMquad[26]+fMquad[25])+0.06096631611034903*fMquad[24]+0.03810394756896815*(fMquad[23]+fMquad[22])+0.06096631611034903*(fMquad[21]+fMquad[20]+fMquad[19])+0.09754610577655842*fMquad[18]+0.03810394756896815*(fMquad[17]+fMquad[16])+0.06096631611034903*fMquad[15]+0.03810394756896815*(fMquad[14]+fMquad[13])+0.06096631611034903*(fMquad[12]+fMquad[11]+fMquad[10])+0.09754610577655842*fMquad[9]+0.06096631611034904*(fMquad[8]+fMquad[7])+0.09754610577655844*fMquad[6]+0.06096631611034904*(fMquad[5]+fMquad[4])+0.09754610577655844*fMquad[3]+0.09754610577655845*(fMquad[2]+fMquad[1])+0.1560737692424935*fMquad[0]; 
  fMOut[1] = 0.03195113136573775*(fMquad[80]+fMquad[79])+0.05112181018518038*fMquad[78]+0.03195113136573775*(fMquad[77]+fMquad[76])+0.05112181018518038*(fMquad[75]+fMquad[74]+fMquad[73])+0.0817948962962886*fMquad[72]+0.03195113136573775*(fMquad[71]+fMquad[70])+0.05112181018518038*fMquad[69]+0.03195113136573775*(fMquad[68]+fMquad[67])+0.05112181018518038*(fMquad[66]+fMquad[65]+fMquad[64])+0.0817948962962886*fMquad[63]+0.05112181018518038*(fMquad[62]+fMquad[61])+0.0817948962962886*fMquad[60]+0.05112181018518038*(fMquad[59]+fMquad[58])+0.0817948962962886*(fMquad[57]+fMquad[56]+fMquad[55])+0.1308718340740617*fMquad[54]-0.03195113136573775*(fMquad[53]+fMquad[52])-0.05112181018518038*fMquad[51]-0.03195113136573775*(fMquad[50]+fMquad[49])-0.05112181018518038*(fMquad[48]+fMquad[47]+fMquad[46])-0.0817948962962886*fMquad[45]-0.03195113136573775*(fMquad[44]+fMquad[43])-0.05112181018518038*fMquad[42]-0.03195113136573775*(fMquad[41]+fMquad[40])-0.05112181018518038*(fMquad[39]+fMquad[38]+fMquad[37])-0.0817948962962886*fMquad[36]-0.05112181018518038*(fMquad[35]+fMquad[34])-0.0817948962962886*fMquad[33]-0.05112181018518038*(fMquad[32]+fMquad[31])-0.0817948962962886*(fMquad[30]+fMquad[29]+fMquad[28])-0.1308718340740617*fMquad[27]; 
  fMOut[2] = 0.03195113136573775*(fMquad[80]+fMquad[79])+0.05112181018518038*fMquad[78]+0.03195113136573775*(fMquad[77]+fMquad[76])+0.05112181018518038*(fMquad[75]+fMquad[74]+fMquad[73])+0.0817948962962886*fMquad[72]-0.03195113136573775*(fMquad[71]+fMquad[70])-0.05112181018518038*fMquad[69]-0.03195113136573775*(fMquad[68]+fMquad[67])-0.05112181018518038*(fMquad[66]+fMquad[65]+fMquad[64])-0.0817948962962886*fMquad[63]+0.03195113136573775*(fMquad[53]+fMquad[52])+0.05112181018518038*fMquad[51]+0.03195113136573775*(fMquad[50]+fMquad[49])+0.05112181018518038*(fMquad[48]+fMquad[47]+fMquad[46])+0.0817948962962886*fMquad[45]-0.03195113136573775*(fMquad[44]+fMquad[43])-0.05112181018518038*fMquad[42]-0.03195113136573775*(fMquad[41]+fMquad[40])-0.05112181018518038*(fMquad[39]+fMquad[38]+fMquad[37])-0.0817948962962886*fMquad[36]+0.05112181018518038*(fMquad[26]+fMquad[25])+0.0817948962962886*fMquad[24]+0.05112181018518038*(fMquad[23]+fMquad[22])+0.0817948962962886*(fMquad[21]+fMquad[20]+fMquad[19])+0.1308718340740617*fMquad[18]-0.05112181018518038*(fMquad[17]+fMquad[16])-0.0817948962962886*fMquad[15]-0.05112181018518038*(fMquad[14]+fMquad[13])-0.0817948962962886*(fMquad[12]+fMquad[11]+fMquad[10])-0.1308718340740617*fMquad[9]; 
  fMOut[3] = 0.03195113136573775*(fMquad[80]+fMquad[79])+0.05112181018518038*fMquad[78]-0.03195113136573775*(fMquad[77]+fMquad[76])-0.05112181018518038*fMquad[75]+0.03195113136573775*(fMquad[71]+fMquad[70])+0.05112181018518038*fMquad[69]-0.03195113136573775*(fMquad[68]+fMquad[67])-0.05112181018518038*fMquad[66]+0.05112181018518038*(fMquad[62]+fMquad[61])+0.0817948962962886*fMquad[60]-0.05112181018518038*(fMquad[59]+fMquad[58])-0.0817948962962886*fMquad[57]+0.03195113136573775*(fMquad[53]+fMquad[52])+0.05112181018518038*fMquad[51]-0.03195113136573775*(fMquad[50]+fMquad[49])-0.05112181018518038*fMquad[48]+0.03195113136573775*(fMquad[44]+fMquad[43])+0.05112181018518038*fMquad[42]-0.03195113136573775*(fMquad[41]+fMquad[40])-0.05112181018518038*fMquad[39]+0.05112181018518038*(fMquad[35]+fMquad[34])+0.0817948962962886*fMquad[33]-0.05112181018518038*(fMquad[32]+fMquad[31])-0.0817948962962886*fMquad[30]+0.05112181018518038*(fMquad[26]+fMquad[25])+0.0817948962962886*fMquad[24]-0.05112181018518038*(fMquad[23]+fMquad[22])-0.0817948962962886*fMquad[21]+0.05112181018518038*(fMquad[17]+fMquad[16])+0.0817948962962886*fMquad[15]-0.05112181018518038*(fMquad[14]+fMquad[13])-0.0817948962962886*fMquad[12]+0.0817948962962886*(fMquad[8]+fMquad[7])+0.1308718340740617*fMquad[6]-0.0817948962962886*(fMquad[5]+fMquad[4])-0.1308718340740617*fMquad[3]; 
  fMOut[4] = 0.03195113136573775*fMquad[80]-0.03195113136573775*fMquad[79]+0.03195113136573775*fMquad[77]-0.03195113136573775*fMquad[76]+0.05112181018518038*fMquad[74]-0.05112181018518038*fMquad[73]+0.03195113136573775*fMquad[71]-0.03195113136573775*fMquad[70]+0.03195113136573775*fMquad[68]-0.03195113136573775*fMquad[67]+0.05112181018518038*fMquad[65]-0.05112181018518038*fMquad[64]+0.05112181018518038*fMquad[62]-0.05112181018518038*fMquad[61]+0.05112181018518038*fMquad[59]-0.05112181018518038*fMquad[58]+0.0817948962962886*fMquad[56]-0.0817948962962886*fMquad[55]+0.03195113136573775*fMquad[53]-0.03195113136573775*fMquad[52]+0.03195113136573775*fMquad[50]-0.03195113136573775*fMquad[49]+0.05112181018518038*fMquad[47]-0.05112181018518038*fMquad[46]+0.03195113136573775*fMquad[44]-0.03195113136573775*fMquad[43]+0.03195113136573775*fMquad[41]-0.03195113136573775*fMquad[40]+0.05112181018518038*fMquad[38]-0.05112181018518038*fMquad[37]+0.05112181018518038*fMquad[35]-0.05112181018518038*fMquad[34]+0.05112181018518038*fMquad[32]-0.05112181018518038*fMquad[31]+0.0817948962962886*fMquad[29]-0.0817948962962886*fMquad[28]+0.05112181018518038*fMquad[26]-0.05112181018518038*fMquad[25]+0.05112181018518038*fMquad[23]-0.05112181018518038*fMquad[22]+0.0817948962962886*fMquad[20]-0.0817948962962886*fMquad[19]+0.05112181018518038*fMquad[17]-0.05112181018518038*fMquad[16]+0.05112181018518038*fMquad[14]-0.05112181018518038*fMquad[13]+0.0817948962962886*fMquad[11]-0.0817948962962886*fMquad[10]+0.0817948962962886*fMquad[8]-0.0817948962962886*fMquad[7]+0.0817948962962886*fMquad[5]-0.0817948962962886*fMquad[4]+0.1308718340740617*fMquad[2]-0.1308718340740617*fMquad[1]; 
  fMOut[5] = 0.04286694101508918*(fMquad[80]+fMquad[79])+0.06858710562414268*fMquad[78]+0.04286694101508918*(fMquad[77]+fMquad[76])+0.06858710562414268*(fMquad[75]+fMquad[74]+fMquad[73])+0.1097393689986283*fMquad[72]-0.04286694101508918*(fMquad[71]+fMquad[70])-0.06858710562414268*fMquad[69]-0.04286694101508918*(fMquad[68]+fMquad[67])-0.06858710562414268*(fMquad[66]+fMquad[65]+fMquad[64])-0.1097393689986283*fMquad[63]-0.04286694101508918*(fMquad[53]+fMquad[52])-0.06858710562414268*fMquad[51]-0.04286694101508918*(fMquad[50]+fMquad[49])-0.06858710562414268*(fMquad[48]+fMquad[47]+fMquad[46])-0.1097393689986283*fMquad[45]+0.04286694101508918*(fMquad[44]+fMquad[43])+0.06858710562414268*fMquad[42]+0.04286694101508918*(fMquad[41]+fMquad[40])+0.06858710562414268*(fMquad[39]+fMquad[38]+fMquad[37])+0.1097393689986283*fMquad[36]; 
  fMOut[6] = 0.04286694101508918*(fMquad[80]+fMquad[79])+0.06858710562414268*fMquad[78]-0.04286694101508918*(fMquad[77]+fMquad[76])-0.06858710562414268*fMquad[75]+0.04286694101508918*(fMquad[71]+fMquad[70])+0.06858710562414268*fMquad[69]-0.04286694101508918*(fMquad[68]+fMquad[67])-0.06858710562414268*fMquad[66]+0.06858710562414268*(fMquad[62]+fMquad[61])+0.1097393689986283*fMquad[60]-0.06858710562414268*(fMquad[59]+fMquad[58])-0.1097393689986283*fMquad[57]-0.04286694101508918*(fMquad[53]+fMquad[52])-0.06858710562414268*fMquad[51]+0.04286694101508918*(fMquad[50]+fMquad[49])+0.06858710562414268*fMquad[48]-0.04286694101508918*(fMquad[44]+fMquad[43])-0.06858710562414268*fMquad[42]+0.04286694101508918*(fMquad[41]+fMquad[40])+0.06858710562414268*fMquad[39]-0.06858710562414268*(fMquad[35]+fMquad[34])-0.1097393689986283*fMquad[33]+0.06858710562414268*(fMquad[32]+fMquad[31])+0.1097393689986283*fMquad[30]; 
  fMOut[7] = 0.04286694101508918*(fMquad[80]+fMquad[79])+0.06858710562414268*fMquad[78]-0.04286694101508918*(fMquad[77]+fMquad[76])-0.06858710562414268*fMquad[75]-0.04286694101508918*(fMquad[71]+fMquad[70])-0.06858710562414268*fMquad[69]+0.04286694101508918*(fMquad[68]+fMquad[67])+0.06858710562414268*fMquad[66]+0.04286694101508918*(fMquad[53]+fMquad[52])+0.06858710562414268*fMquad[51]-0.04286694101508918*(fMquad[50]+fMquad[49])-0.06858710562414268*fMquad[48]-0.04286694101508918*(fMquad[44]+fMquad[43])-0.06858710562414268*fMquad[42]+0.04286694101508918*(fMquad[41]+fMquad[40])+0.06858710562414268*(fMquad[39]+fMquad[26]+fMquad[25])+0.1097393689986283*fMquad[24]-0.06858710562414268*(fMquad[23]+fMquad[22])-0.1097393689986283*fMquad[21]-0.06858710562414268*(fMquad[17]+fMquad[16])-0.1097393689986283*fMquad[15]+0.06858710562414268*(fMquad[14]+fMquad[13])+0.1097393689986283*fMquad[12]; 
  fMOut[8] = 0.04286694101508918*fMquad[80]-0.04286694101508918*fMquad[79]+0.04286694101508918*fMquad[77]-0.04286694101508918*fMquad[76]+0.06858710562414268*fMquad[74]-0.06858710562414268*fMquad[73]+0.04286694101508918*fMquad[71]-0.04286694101508918*fMquad[70]+0.04286694101508918*fMquad[68]-0.04286694101508918*fMquad[67]+0.06858710562414268*fMquad[65]-0.06858710562414268*fMquad[64]+0.06858710562414268*fMquad[62]-0.06858710562414268*fMquad[61]+0.06858710562414268*fMquad[59]-0.06858710562414268*fMquad[58]+0.1097393689986283*fMquad[56]-0.1097393689986283*fMquad[55]-0.04286694101508918*fMquad[53]+0.04286694101508918*fMquad[52]-0.04286694101508918*fMquad[50]+0.04286694101508918*fMquad[49]-0.06858710562414268*fMquad[47]+0.06858710562414268*fMquad[46]-0.04286694101508918*fMquad[44]+0.04286694101508918*fMquad[43]-0.04286694101508918*fMquad[41]+0.04286694101508918*fMquad[40]-0.06858710562414268*fMquad[38]+0.06858710562414268*fMquad[37]-0.06858710562414268*fMquad[35]+0.06858710562414268*fMquad[34]-0.06858710562414268*fMquad[32]+0.06858710562414268*fMquad[31]-0.1097393689986283*fMquad[29]+0.1097393689986283*fMquad[28]; 
  fMOut[9] = 0.04286694101508918*fMquad[80]-0.04286694101508918*fMquad[79]+0.04286694101508918*fMquad[77]-0.04286694101508918*fMquad[76]+0.06858710562414268*fMquad[74]-0.06858710562414268*fMquad[73]-0.04286694101508918*fMquad[71]+0.04286694101508918*fMquad[70]-0.04286694101508918*fMquad[68]+0.04286694101508918*fMquad[67]-0.06858710562414268*fMquad[65]+0.06858710562414268*fMquad[64]+0.04286694101508918*fMquad[53]-0.04286694101508918*fMquad[52]+0.04286694101508918*fMquad[50]-0.04286694101508918*fMquad[49]+0.06858710562414268*fMquad[47]-0.06858710562414268*fMquad[46]-0.04286694101508918*fMquad[44]+0.04286694101508918*fMquad[43]-0.04286694101508918*fMquad[41]+0.04286694101508918*fMquad[40]-0.06858710562414268*fMquad[38]+0.06858710562414268*(fMquad[37]+fMquad[26])-0.06858710562414268*fMquad[25]+0.06858710562414268*fMquad[23]-0.06858710562414268*fMquad[22]+0.1097393689986283*fMquad[20]-0.1097393689986283*fMquad[19]-0.06858710562414268*fMquad[17]+0.06858710562414268*fMquad[16]-0.06858710562414268*fMquad[14]+0.06858710562414268*fMquad[13]-0.1097393689986283*fMquad[11]+0.1097393689986283*fMquad[10]; 
  fMOut[10] = 0.04286694101508918*fMquad[80]-0.04286694101508918*(fMquad[79]+fMquad[77])+0.04286694101508918*(fMquad[76]+fMquad[71])-0.04286694101508918*(fMquad[70]+fMquad[68])+0.04286694101508918*fMquad[67]+0.06858710562414268*fMquad[62]-0.06858710562414268*(fMquad[61]+fMquad[59])+0.06858710562414268*fMquad[58]+0.04286694101508918*fMquad[53]-0.04286694101508918*(fMquad[52]+fMquad[50])+0.04286694101508918*(fMquad[49]+fMquad[44])-0.04286694101508918*(fMquad[43]+fMquad[41])+0.04286694101508918*fMquad[40]+0.06858710562414268*fMquad[35]-0.06858710562414268*(fMquad[34]+fMquad[32])+0.06858710562414268*(fMquad[31]+fMquad[26])-0.06858710562414268*(fMquad[25]+fMquad[23])+0.06858710562414268*(fMquad[22]+fMquad[17])-0.06858710562414268*(fMquad[16]+fMquad[14])+0.06858710562414268*fMquad[13]+0.1097393689986283*fMquad[8]-0.1097393689986283*(fMquad[7]+fMquad[5])+0.1097393689986283*fMquad[4]; 
  fMOut[11] = 0.02130075424382517*(fMquad[80]+fMquad[79])+0.03408120679012027*fMquad[78]+0.02130075424382517*(fMquad[77]+fMquad[76])+0.03408120679012027*(fMquad[75]+fMquad[74]+fMquad[73])+0.05452993086419242*fMquad[72]+0.02130075424382517*(fMquad[71]+fMquad[70])+0.03408120679012027*fMquad[69]+0.02130075424382517*(fMquad[68]+fMquad[67])+0.03408120679012027*(fMquad[66]+fMquad[65]+fMquad[64])+0.05452993086419242*fMquad[63]+0.03408120679012027*(fMquad[62]+fMquad[61])+0.05452993086419242*fMquad[60]+0.03408120679012027*(fMquad[59]+fMquad[58])+0.05452993086419242*(fMquad[57]+fMquad[56]+fMquad[55])+0.08724788938270785*fMquad[54]+0.02130075424382517*(fMquad[53]+fMquad[52])+0.03408120679012027*fMquad[51]+0.02130075424382517*(fMquad[50]+fMquad[49])+0.03408120679012027*(fMquad[48]+fMquad[47]+fMquad[46])+0.05452993086419242*fMquad[45]+0.02130075424382517*(fMquad[44]+fMquad[43])+0.03408120679012027*fMquad[42]+0.02130075424382517*(fMquad[41]+fMquad[40])+0.03408120679012027*(fMquad[39]+fMquad[38]+fMquad[37])+0.05452993086419242*fMquad[36]+0.03408120679012027*(fMquad[35]+fMquad[34])+0.05452993086419242*fMquad[33]+0.03408120679012027*(fMquad[32]+fMquad[31])+0.05452993086419242*(fMquad[30]+fMquad[29]+fMquad[28])+0.08724788938270785*fMquad[27]-0.04260150848765032*(fMquad[26]+fMquad[25])-0.0681624135802405*fMquad[24]-0.04260150848765032*(fMquad[23]+fMquad[22])-0.0681624135802405*(fMquad[21]+fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]-0.04260150848765032*(fMquad[17]+fMquad[16])-0.0681624135802405*fMquad[15]-0.04260150848765032*(fMquad[14]+fMquad[13])-0.0681624135802405*(fMquad[12]+fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]-0.06816241358024051*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]-0.06816241358024051*(fMquad[5]+fMquad[4])-0.1090598617283848*fMquad[3]-0.1090598617283848*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[12] = 0.02130075424382517*(fMquad[80]+fMquad[79])+0.03408120679012027*fMquad[78]+0.02130075424382517*(fMquad[77]+fMquad[76])+0.03408120679012027*(fMquad[75]+fMquad[74]+fMquad[73])+0.05452993086419242*fMquad[72]+0.02130075424382517*(fMquad[71]+fMquad[70])+0.03408120679012027*fMquad[69]+0.02130075424382517*(fMquad[68]+fMquad[67])+0.03408120679012027*(fMquad[66]+fMquad[65]+fMquad[64])+0.05452993086419242*fMquad[63]-0.04260150848765032*(fMquad[62]+fMquad[61])-0.0681624135802405*fMquad[60]-0.04260150848765032*(fMquad[59]+fMquad[58])-0.0681624135802405*(fMquad[57]+fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382517*(fMquad[53]+fMquad[52])+0.03408120679012027*fMquad[51]+0.02130075424382517*(fMquad[50]+fMquad[49])+0.03408120679012027*(fMquad[48]+fMquad[47]+fMquad[46])+0.05452993086419242*fMquad[45]+0.02130075424382517*(fMquad[44]+fMquad[43])+0.03408120679012027*fMquad[42]+0.02130075424382517*(fMquad[41]+fMquad[40])+0.03408120679012027*(fMquad[39]+fMquad[38]+fMquad[37])+0.05452993086419242*fMquad[36]-0.04260150848765032*(fMquad[35]+fMquad[34])-0.0681624135802405*fMquad[33]-0.04260150848765032*(fMquad[32]+fMquad[31])-0.0681624135802405*(fMquad[30]+fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012027*(fMquad[26]+fMquad[25])+0.05452993086419242*fMquad[24]+0.03408120679012027*(fMquad[23]+fMquad[22])+0.05452993086419242*(fMquad[21]+fMquad[20]+fMquad[19])+0.08724788938270785*fMquad[18]+0.03408120679012027*(fMquad[17]+fMquad[16])+0.05452993086419242*fMquad[15]+0.03408120679012027*(fMquad[14]+fMquad[13])+0.05452993086419242*(fMquad[12]+fMquad[11]+fMquad[10])+0.08724788938270785*fMquad[9]-0.06816241358024051*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]-0.06816241358024051*(fMquad[5]+fMquad[4])-0.1090598617283848*fMquad[3]-0.1090598617283848*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[13] = 0.02130075424382517*(fMquad[80]+fMquad[79])+0.03408120679012027*fMquad[78]+0.02130075424382517*(fMquad[77]+fMquad[76])+0.03408120679012027*fMquad[75]-0.04260150848765032*(fMquad[74]+fMquad[73])-0.0681624135802405*fMquad[72]+0.02130075424382517*(fMquad[71]+fMquad[70])+0.03408120679012027*fMquad[69]+0.02130075424382517*(fMquad[68]+fMquad[67])+0.03408120679012027*fMquad[66]-0.04260150848765032*(fMquad[65]+fMquad[64])-0.0681624135802405*fMquad[63]+0.03408120679012027*(fMquad[62]+fMquad[61])+0.05452993086419242*fMquad[60]+0.03408120679012027*(fMquad[59]+fMquad[58])+0.05452993086419242*fMquad[57]-0.0681624135802405*(fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382517*(fMquad[53]+fMquad[52])+0.03408120679012027*fMquad[51]+0.02130075424382517*(fMquad[50]+fMquad[49])+0.03408120679012027*fMquad[48]-0.04260150848765032*(fMquad[47]+fMquad[46])-0.0681624135802405*fMquad[45]+0.02130075424382517*(fMquad[44]+fMquad[43])+0.03408120679012027*fMquad[42]+0.02130075424382517*(fMquad[41]+fMquad[40])+0.03408120679012027*fMquad[39]-0.04260150848765032*(fMquad[38]+fMquad[37])-0.0681624135802405*fMquad[36]+0.03408120679012027*(fMquad[35]+fMquad[34])+0.05452993086419242*fMquad[33]+0.03408120679012027*(fMquad[32]+fMquad[31])+0.05452993086419242*fMquad[30]-0.0681624135802405*(fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012027*(fMquad[26]+fMquad[25])+0.05452993086419242*fMquad[24]+0.03408120679012027*(fMquad[23]+fMquad[22])+0.05452993086419242*fMquad[21]-0.0681624135802405*(fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]+0.03408120679012027*(fMquad[17]+fMquad[16])+0.05452993086419242*fMquad[15]+0.03408120679012027*(fMquad[14]+fMquad[13])+0.05452993086419242*fMquad[12]-0.0681624135802405*(fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]+0.05452993086419243*(fMquad[8]+fMquad[7])+0.08724788938270786*fMquad[6]+0.05452993086419243*(fMquad[5]+fMquad[4])+0.08724788938270786*fMquad[3]-0.1090598617283848*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 
  fMOut[14] = 0.02130075424382517*(fMquad[80]+fMquad[79])-0.04260150848765032*fMquad[78]+0.02130075424382517*(fMquad[77]+fMquad[76])-0.04260150848765032*fMquad[75]+0.03408120679012027*(fMquad[74]+fMquad[73])-0.0681624135802405*fMquad[72]+0.02130075424382517*(fMquad[71]+fMquad[70])-0.04260150848765032*fMquad[69]+0.02130075424382517*(fMquad[68]+fMquad[67])-0.04260150848765032*fMquad[66]+0.03408120679012027*(fMquad[65]+fMquad[64])-0.0681624135802405*fMquad[63]+0.03408120679012027*(fMquad[62]+fMquad[61])-0.0681624135802405*fMquad[60]+0.03408120679012027*(fMquad[59]+fMquad[58])-0.0681624135802405*fMquad[57]+0.05452993086419242*(fMquad[56]+fMquad[55])-0.1090598617283848*fMquad[54]+0.02130075424382517*(fMquad[53]+fMquad[52])-0.04260150848765032*fMquad[51]+0.02130075424382517*(fMquad[50]+fMquad[49])-0.04260150848765032*fMquad[48]+0.03408120679012027*(fMquad[47]+fMquad[46])-0.0681624135802405*fMquad[45]+0.02130075424382517*(fMquad[44]+fMquad[43])-0.04260150848765032*fMquad[42]+0.02130075424382517*(fMquad[41]+fMquad[40])-0.04260150848765032*fMquad[39]+0.03408120679012027*(fMquad[38]+fMquad[37])-0.0681624135802405*fMquad[36]+0.03408120679012027*(fMquad[35]+fMquad[34])-0.0681624135802405*fMquad[33]+0.03408120679012027*(fMquad[32]+fMquad[31])-0.0681624135802405*fMquad[30]+0.05452993086419242*(fMquad[29]+fMquad[28])-0.1090598617283848*fMquad[27]+0.03408120679012027*(fMquad[26]+fMquad[25])-0.0681624135802405*fMquad[24]+0.03408120679012027*(fMquad[23]+fMquad[22])-0.0681624135802405*fMquad[21]+0.05452993086419242*(fMquad[20]+fMquad[19])-0.1090598617283848*fMquad[18]+0.03408120679012027*(fMquad[17]+fMquad[16])-0.0681624135802405*fMquad[15]+0.03408120679012027*(fMquad[14]+fMquad[13])-0.0681624135802405*fMquad[12]+0.05452993086419242*(fMquad[11]+fMquad[10])-0.1090598617283848*fMquad[9]+0.05452993086419243*(fMquad[8]+fMquad[7])-0.1090598617283848*fMquad[6]+0.05452993086419243*(fMquad[5]+fMquad[4])-0.1090598617283848*fMquad[3]+0.08724788938270787*(fMquad[2]+fMquad[1])-0.1744957787654157*fMquad[0]; 

}
