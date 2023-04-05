#include <BiMaxwellianOnBasisModDecl.h>

void GkBiMaxwellianOnBasisGauss1x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtparSq, const double *vtperpSq, const double *bmag, double *flowUOrd, double *vtparSqOrd, double *vtperpSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*den[0]-0.7071067811865474*den[1]; 
  m0Ord[1] = 0.7071067811865474*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7071067811865475*flowU[0]-0.7071067811865474*flowU[1]; 
  flowUOrd[1] = 0.7071067811865474*flowU[1]+0.7071067811865475*flowU[0]; 

  vtparSqOrd[0] = 0.7071067811865475*vtparSq[0]-0.7071067811865474*vtparSq[1]; 
  vtparSqOrd[1] = 0.7071067811865474*vtparSq[1]+0.7071067811865475*vtparSq[0]; 

  vtperpSqOrd[0] = 0.7071067811865475*vtperpSq[0]-0.7071067811865474*vtperpSq[1]; 
  vtperpSqOrd[1] = 0.7071067811865474*vtperpSq[1]+0.7071067811865475*vtperpSq[0]; 

  bmagOrd[0] = 0.7071067811865475*bmag[0]-0.7071067811865474*bmag[1]; 
  bmagOrd[1] = 0.7071067811865474*bmag[1]+0.7071067811865475*bmag[0]; 

  if ((vtparSqOrd[0] > 0.) && (vtperpSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (0.06349363593424097*bmagOrd[0]*m0Ord[0])/(sqrt(vtparSqOrd[0])*vtperpSqOrd[0]); 
  else
    fMFacOrd[0] = 0.0;
  if ((vtparSqOrd[1] > 0.) && (vtperpSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (0.06349363593424097*bmagOrd[1]*m0Ord[1])/(sqrt(vtparSqOrd[1])*vtperpSqOrd[1]); 
  else
    fMFacOrd[1] = 0.0;

}

void GkBiMaxwellianOnBasisGauss1x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtparSqOrd, const double *vtperpSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[8];
  if ((vtparSqOrd[0] > 0.) && (vtperpSqOrd[0] > 0.) && (fMFacOrd[0] > 0.)) {
    fMquad[0] = fMFacOrd[0]*exp((-(1.0*bmagOrd[0]*std::abs(wc[2]-0.2886751345948129*dxv[2]))/(vtperpSqOrd[0]*m_))-(0.5*std::pow(wc[1]-0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0))/vtparSqOrd[0])+9.999999999999999e-41; 
    fMquad[1] = fMFacOrd[0]*exp((-(1.0*bmagOrd[0]*std::abs(wc[2]+0.2886751345948129*dxv[2]))/(vtperpSqOrd[0]*m_))-(0.5*std::pow(wc[1]-0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0))/vtparSqOrd[0])+9.999999999999999e-41; 
    fMquad[2] = fMFacOrd[0]*exp((-(1.0*bmagOrd[0]*std::abs(wc[2]-0.2886751345948129*dxv[2]))/(vtperpSqOrd[0]*m_))-(0.5*std::pow(wc[1]+0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0))/vtparSqOrd[0])+9.999999999999999e-41; 
    fMquad[3] = fMFacOrd[0]*exp((-(1.0*bmagOrd[0]*std::abs(wc[2]+0.2886751345948129*dxv[2]))/(vtperpSqOrd[0]*m_))-(0.5*std::pow(wc[1]+0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0))/vtparSqOrd[0])+9.999999999999999e-41; 
  } else {
    fMquad[0] = 9.999999999999999e-41;
    fMquad[1] = 9.999999999999999e-41;
    fMquad[2] = 9.999999999999999e-41;
    fMquad[3] = 9.999999999999999e-41;
  };
  if ((vtparSqOrd[1] > 0.) && (vtperpSqOrd[1] > 0.) && (fMFacOrd[1] > 0.)) {
    fMquad[4] = fMFacOrd[1]*exp((-(1.0*bmagOrd[1]*std::abs(wc[2]-0.2886751345948129*dxv[2]))/(vtperpSqOrd[1]*m_))-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.2886751345948129*dxv[1],2.0))/vtparSqOrd[1])+9.999999999999999e-41; 
    fMquad[5] = fMFacOrd[1]*exp((-(1.0*bmagOrd[1]*std::abs(wc[2]+0.2886751345948129*dxv[2]))/(vtperpSqOrd[1]*m_))-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.2886751345948129*dxv[1],2.0))/vtparSqOrd[1])+9.999999999999999e-41; 
    fMquad[6] = fMFacOrd[1]*exp((-(1.0*bmagOrd[1]*std::abs(wc[2]-0.2886751345948129*dxv[2]))/(vtperpSqOrd[1]*m_))-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.2886751345948129*dxv[1],2.0))/vtparSqOrd[1])+9.999999999999999e-41; 
    fMquad[7] = fMFacOrd[1]*exp((-(1.0*bmagOrd[1]*std::abs(wc[2]+0.2886751345948129*dxv[2]))/(vtperpSqOrd[1]*m_))-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.2886751345948129*dxv[1],2.0))/vtparSqOrd[1])+9.999999999999999e-41; 
  } else {
    fMquad[4] = 9.999999999999999e-41;
    fMquad[5] = 9.999999999999999e-41;
    fMquad[6] = 9.999999999999999e-41;
    fMquad[7] = 9.999999999999999e-41;
  };

  fMOut[0] = 0.3535533905932737*(fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[1] = 0.3535533905932736*(fMquad[7]+fMquad[6]+fMquad[5]+fMquad[4])-0.3535533905932736*(fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[2] = 0.3535533905932736*(fMquad[7]+fMquad[6])-0.3535533905932736*(fMquad[5]+fMquad[4])+0.3535533905932736*(fMquad[3]+fMquad[2])-0.3535533905932736*(fMquad[1]+fMquad[0]); 
  fMOut[3] = 0.3535533905932736*fMquad[7]-0.3535533905932736*fMquad[6]+0.3535533905932736*fMquad[5]-0.3535533905932736*fMquad[4]+0.3535533905932736*fMquad[3]-0.3535533905932736*fMquad[2]+0.3535533905932736*fMquad[1]-0.3535533905932736*fMquad[0]; 
  fMOut[4] = 0.3535533905932736*(fMquad[7]+fMquad[6])-0.3535533905932736*(fMquad[5]+fMquad[4]+fMquad[3]+fMquad[2])+0.3535533905932736*(fMquad[1]+fMquad[0]); 
  fMOut[5] = 0.3535533905932736*fMquad[7]-0.3535533905932736*fMquad[6]+0.3535533905932736*fMquad[5]-0.3535533905932736*(fMquad[4]+fMquad[3])+0.3535533905932736*fMquad[2]-0.3535533905932736*fMquad[1]+0.3535533905932736*fMquad[0]; 
  fMOut[6] = 0.3535533905932736*fMquad[7]-0.3535533905932736*(fMquad[6]+fMquad[5])+0.3535533905932736*(fMquad[4]+fMquad[3])-0.3535533905932736*(fMquad[2]+fMquad[1])+0.3535533905932736*fMquad[0]; 
  fMOut[7] = 0.3535533905932736*fMquad[7]-0.3535533905932736*(fMquad[6]+fMquad[5])+0.3535533905932736*fMquad[4]-0.3535533905932736*fMquad[3]+0.3535533905932736*(fMquad[2]+fMquad[1])-0.3535533905932736*fMquad[0]; 

}
