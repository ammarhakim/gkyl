#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss1x3vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd) {

  flowUOrd[0] = 0.7071067811865475*flowU[0]-0.7071067811865474*flowU[1]; 
  flowUOrd[1] = 0.7071067811865474*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[2] = 0.7071067811865475*flowU[2]-0.7071067811865474*flowU[3]; 
  flowUOrd[3] = 0.7071067811865474*flowU[3]+0.7071067811865475*flowU[2]; 
  flowUOrd[4] = 0.7071067811865475*flowU[4]-0.7071067811865474*flowU[5]; 
  flowUOrd[5] = 0.7071067811865474*flowU[5]+0.7071067811865475*flowU[4]; 

  vtSqOrd[0] = 0.7071067811865475*vtSq[0]-0.7071067811865474*vtSq[1]; 
  vtSqOrd[1] = 0.7071067811865474*vtSq[1]+0.7071067811865475*vtSq[0]; 

  fMFacOrd[0] = (0.7071067811865475*den[0]-0.7071067811865474*den[1])/std::pow(2.506628274631001*sqrt(0.7071067811865475*vtSq[0]-0.7071067811865474*vtSq[1]),3.0); 
  fMFacOrd[1] = (0.7071067811865474*den[1]+0.7071067811865475*den[0])/std::pow(2.506628274631001*sqrt(0.7071067811865474*vtSq[1]+0.7071067811865475*vtSq[0]),3.0); 

}
void MaxwellianOnBasisGauss1x3vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[16];
  fMquad[0] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[1] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[2] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[3] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[4] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]+0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[5] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]-0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]+0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[6] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]+0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[7] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow(wc[2]-1.0*flowUOrd[2]+0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]+0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
  fMquad[8] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]-0.2886751345948129*dxv[1],2.0)))/vtSqOrd[1]); 
  fMquad[9] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]-0.2886751345948129*dxv[1],2.0)))/vtSqOrd[1]); 
  fMquad[10] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]-0.2886751345948129*dxv[1],2.0)))/vtSqOrd[1]); 
  fMquad[11] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]-0.2886751345948129*dxv[1],2.0)))/vtSqOrd[1]); 
  fMquad[12] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]+0.2886751345948129*dxv[1],2.0)))/vtSqOrd[1]); 
  fMquad[13] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]-0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]+0.2886751345948129*dxv[1],2.0)))/vtSqOrd[1]); 
  fMquad[14] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]-0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]+0.2886751345948129*dxv[1],2.0)))/vtSqOrd[1]); 
  fMquad[15] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[3]+0.2886751345948129*dxv[3],2.0)+std::pow((-1.0*flowUOrd[3])+wc[2]+0.2886751345948129*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]+0.2886751345948129*dxv[1],2.0)))/vtSqOrd[1]); 

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
