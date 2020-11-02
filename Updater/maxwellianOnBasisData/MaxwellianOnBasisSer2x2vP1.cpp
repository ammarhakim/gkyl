#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss2x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd) {

  flowUOrd[0] = 0.4999999999999999*flowU[3]-0.4999999999999999*(flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[1] = (-0.4999999999999999*flowU[3])+0.4999999999999999*flowU[2]-0.4999999999999999*flowU[1]+0.5*flowU[0]; 
  flowUOrd[2] = (-0.4999999999999999*(flowU[3]+flowU[2]))+0.4999999999999999*flowU[1]+0.5*flowU[0]; 
  flowUOrd[3] = 0.4999999999999999*(flowU[3]+flowU[2]+flowU[1])+0.5*flowU[0]; 
  flowUOrd[4] = 0.4999999999999999*flowU[7]-0.4999999999999999*(flowU[6]+flowU[5])+0.5*flowU[4]; 
  flowUOrd[5] = (-0.4999999999999999*flowU[7])+0.4999999999999999*flowU[6]-0.4999999999999999*flowU[5]+0.5*flowU[4]; 
  flowUOrd[6] = (-0.4999999999999999*(flowU[7]+flowU[6]))+0.4999999999999999*flowU[5]+0.5*flowU[4]; 
  flowUOrd[7] = 0.4999999999999999*(flowU[7]+flowU[6]+flowU[5])+0.5*flowU[4]; 

  vtSqOrd[0] = 0.4999999999999999*vtSq[3]-0.4999999999999999*(vtSq[2]+vtSq[1])+0.5*vtSq[0]; 
  vtSqOrd[1] = (-0.4999999999999999*vtSq[3])+0.4999999999999999*vtSq[2]-0.4999999999999999*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[2] = (-0.4999999999999999*(vtSq[3]+vtSq[2]))+0.4999999999999999*vtSq[1]+0.5*vtSq[0]; 
  vtSqOrd[3] = 0.4999999999999999*(vtSq[3]+vtSq[2]+vtSq[1])+0.5*vtSq[0]; 

  fMFacOrd[0] = (0.1591549430918953*(0.4999999999999999*den[3]-0.4999999999999999*(den[2]+den[1])+0.5*den[0]))/(0.4999999999999999*vtSq[3]-0.4999999999999999*(vtSq[2]+vtSq[1])+0.5*vtSq[0]); 
  fMFacOrd[1] = (0.1591549430918953*((-0.4999999999999999*den[3])+0.4999999999999999*den[2]-0.4999999999999999*den[1]+0.5*den[0]))/((-0.4999999999999999*vtSq[3])+0.4999999999999999*vtSq[2]-0.4999999999999999*vtSq[1]+0.5*vtSq[0]); 
  fMFacOrd[2] = (0.1591549430918953*((-0.4999999999999999*(den[3]+den[2]))+0.4999999999999999*den[1]+0.5*den[0]))/((-0.4999999999999999*(vtSq[3]+vtSq[2]))+0.4999999999999999*vtSq[1]+0.5*vtSq[0]); 
  fMFacOrd[3] = (0.1591549430918953*(0.4999999999999999*(den[3]+den[2]+den[1])+0.5*den[0]))/(0.4999999999999999*(vtSq[3]+vtSq[2]+vtSq[1])+0.5*vtSq[0]); 

}
void MaxwellianOnBasisGauss2x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[16];
  fMquad[0] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[1] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[2] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[3] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[4] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[5] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[6] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[7] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]+0.2886751345948129*dxv[2]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
  fMquad[8] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[6])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[9] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[6])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[10] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[6])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[11] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[6])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[2]); 
  fMquad[12] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[7])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[13] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[7])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[14] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[7])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 
  fMquad[15] = fMFacOrd[3]*exp(-(0.5*(std::pow((-1.0*flowUOrd[7])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.2886751345948129*dxv[2],2.0)))/vtSqOrd[3]); 

  fMOut[0] = 0.25*(fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[1] = 0.25*(fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8])-0.25*(fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[2] = 0.25*(fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12])-0.25*(fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8])+0.25*(fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4])-0.25*(fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[3] = 0.25*(fMquad[15]+fMquad[14])-0.25*(fMquad[13]+fMquad[12])+0.25*(fMquad[11]+fMquad[10])-0.25*(fMquad[9]+fMquad[8])+0.25*(fMquad[7]+fMquad[6])-0.25*(fMquad[5]+fMquad[4])+0.25*(fMquad[3]+fMquad[2])-0.25*(fMquad[1]+fMquad[0]); 
  fMOut[4] = 0.25*fMquad[15]-0.25*fMquad[14]+0.25*fMquad[13]-0.25*fMquad[12]+0.25*fMquad[11]-0.25*fMquad[10]+0.25*fMquad[9]-0.25*fMquad[8]+0.25*fMquad[7]-0.25*fMquad[6]+0.25*fMquad[5]-0.25*fMquad[4]+0.25*fMquad[3]-0.25*fMquad[2]+0.25*fMquad[1]-0.25*fMquad[0]; 
  fMOut[5] = 0.25*(fMquad[15]+fMquad[14]+fMquad[13]+fMquad[12])-0.25*(fMquad[11]+fMquad[10]+fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4])+0.25*(fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[6] = 0.25*(fMquad[15]+fMquad[14])-0.25*(fMquad[13]+fMquad[12])+0.25*(fMquad[11]+fMquad[10])-0.25*(fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6])+0.25*(fMquad[5]+fMquad[4])-0.25*(fMquad[3]+fMquad[2])+0.25*(fMquad[1]+fMquad[0]); 
  fMOut[7] = 0.25*(fMquad[15]+fMquad[14])-0.25*(fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10])+0.25*(fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6])-0.25*(fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2])+0.25*(fMquad[1]+fMquad[0]); 
  fMOut[8] = 0.25*fMquad[15]-0.25*fMquad[14]+0.25*fMquad[13]-0.25*fMquad[12]+0.25*fMquad[11]-0.25*fMquad[10]+0.25*fMquad[9]-0.25*(fMquad[8]+fMquad[7])+0.25*fMquad[6]-0.25*fMquad[5]+0.25*fMquad[4]-0.25*fMquad[3]+0.25*fMquad[2]-0.25*fMquad[1]+0.25*fMquad[0]; 
  fMOut[9] = 0.25*fMquad[15]-0.25*fMquad[14]+0.25*fMquad[13]-0.25*(fMquad[12]+fMquad[11])+0.25*fMquad[10]-0.25*fMquad[9]+0.25*(fMquad[8]+fMquad[7])-0.25*fMquad[6]+0.25*fMquad[5]-0.25*(fMquad[4]+fMquad[3])+0.25*fMquad[2]-0.25*fMquad[1]+0.25*fMquad[0]; 
  fMOut[10] = 0.25*fMquad[15]-0.25*(fMquad[14]+fMquad[13])+0.25*(fMquad[12]+fMquad[11])-0.25*(fMquad[10]+fMquad[9])+0.25*(fMquad[8]+fMquad[7])-0.25*(fMquad[6]+fMquad[5])+0.25*(fMquad[4]+fMquad[3])-0.25*(fMquad[2]+fMquad[1])+0.25*fMquad[0]; 
  fMOut[11] = 0.2499999999999999*(fMquad[15]+fMquad[14])-0.2499999999999999*(fMquad[13]+fMquad[12]+fMquad[11]+fMquad[10])+0.2499999999999999*(fMquad[9]+fMquad[8])-0.2499999999999999*(fMquad[7]+fMquad[6])+0.2499999999999999*(fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2])-0.2499999999999999*(fMquad[1]+fMquad[0]); 
  fMOut[12] = 0.2499999999999999*fMquad[15]-0.2499999999999999*fMquad[14]+0.2499999999999999*fMquad[13]-0.2499999999999999*(fMquad[12]+fMquad[11])+0.2499999999999999*fMquad[10]-0.2499999999999999*fMquad[9]+0.2499999999999999*fMquad[8]-0.2499999999999999*fMquad[7]+0.2499999999999999*fMquad[6]-0.2499999999999999*fMquad[5]+0.2499999999999999*(fMquad[4]+fMquad[3])-0.2499999999999999*fMquad[2]+0.2499999999999999*fMquad[1]-0.2499999999999999*fMquad[0]; 
  fMOut[13] = 0.2499999999999999*fMquad[15]-0.2499999999999999*(fMquad[14]+fMquad[13])+0.2499999999999999*(fMquad[12]+fMquad[11])-0.2499999999999999*(fMquad[10]+fMquad[9])+0.2499999999999999*fMquad[8]-0.2499999999999999*fMquad[7]+0.2499999999999999*(fMquad[6]+fMquad[5])-0.2499999999999999*(fMquad[4]+fMquad[3])+0.2499999999999999*(fMquad[2]+fMquad[1])-0.2499999999999999*fMquad[0]; 
  fMOut[14] = 0.2499999999999999*fMquad[15]-0.2499999999999999*(fMquad[14]+fMquad[13])+0.2499999999999999*fMquad[12]-0.2499999999999999*fMquad[11]+0.2499999999999999*(fMquad[10]+fMquad[9])-0.2499999999999999*fMquad[8]+0.2499999999999999*fMquad[7]-0.2499999999999999*(fMquad[6]+fMquad[5])+0.2499999999999999*fMquad[4]-0.2499999999999999*fMquad[3]+0.2499999999999999*(fMquad[2]+fMquad[1])-0.2499999999999999*fMquad[0]; 
  fMOut[15] = 0.25*fMquad[15]-0.25*(fMquad[14]+fMquad[13])+0.25*fMquad[12]-0.25*fMquad[11]+0.25*(fMquad[10]+fMquad[9])-0.25*(fMquad[8]+fMquad[7])+0.25*(fMquad[6]+fMquad[5])-0.25*fMquad[4]+0.25*fMquad[3]-0.25*(fMquad[2]+fMquad[1])+0.25*fMquad[0]; 

}
