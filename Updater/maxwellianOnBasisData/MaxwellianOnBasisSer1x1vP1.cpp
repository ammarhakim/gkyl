#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss1x1vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*den[0]-0.7071067811865474*den[1]; 
  m0Ord[1] = 0.7071067811865474*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7071067811865475*flowU[0]-0.7071067811865474*flowU[1]; 
  flowUOrd[1] = 0.7071067811865474*flowU[1]+0.7071067811865475*flowU[0]; 

  vtSqOrd[0] = 0.7071067811865475*vtSq[0]-0.7071067811865474*vtSq[1]; 
  vtSqOrd[1] = 0.7071067811865474*vtSq[1]+0.7071067811865475*vtSq[0]; 

  if ((vtSqOrd[0] <= 0.0) || (m0Ord[0] <= 0.0))
    fMFacOrd[0] = 0.;
  else
    fMFacOrd[0] = (0.3989422804014326*m0Ord[0])/sqrt(vtSqOrd[0]); 
  if ((vtSqOrd[1] <= 0.0) || (m0Ord[1] <= 0.0))
    fMFacOrd[1] = 0.;
  else
    fMFacOrd[1] = (0.3989422804014326*m0Ord[1])/sqrt(vtSqOrd[1]); 

}

void MaxwellianOnBasisGauss1x1vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[4];
  if ((vtSqOrd[0] <= 0.0) || (fMFacOrd[0] <= 0.0)) {
    fMquad[0] = 0;
    fMquad[1] = 0;
  } else {
    fMquad[0] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
    fMquad[1] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]+0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  };
  if ((vtSqOrd[1] <= 0.0) || (fMFacOrd[1] <= 0.0)) {
    fMquad[2] = 0;
    fMquad[3] = 0;
  } else {
    fMquad[2] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.2886751345948129*dxv[1],2.0))/vtSqOrd[1]); 
    fMquad[3] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.2886751345948129*dxv[1],2.0))/vtSqOrd[1]); 
  };

  fMOut[0] = 0.5*(fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[1] = 0.4999999999999999*(fMquad[3]+fMquad[2])-0.4999999999999999*(fMquad[1]+fMquad[0]); 
  fMOut[2] = 0.4999999999999999*fMquad[3]-0.4999999999999999*fMquad[2]+0.4999999999999999*fMquad[1]-0.4999999999999999*fMquad[0]; 
  fMOut[3] = 0.4999999999999999*fMquad[3]-0.4999999999999999*(fMquad[2]+fMquad[1])+0.4999999999999999*fMquad[0]; 

}
void GkMaxwellianOnBasisGauss1x1vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*den[0]-0.7071067811865474*den[1]; 
  m0Ord[1] = 0.7071067811865474*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7071067811865475*flowU[0]-0.7071067811865474*flowU[1]; 
  flowUOrd[1] = 0.7071067811865474*flowU[1]+0.7071067811865475*flowU[0]; 

  vtSqOrd[0] = 0.7071067811865475*vtSq[0]-0.7071067811865474*vtSq[1]; 
  vtSqOrd[1] = 0.7071067811865474*vtSq[1]+0.7071067811865475*vtSq[0]; 

  if ((vtSqOrd[0] <= 0.0) || (m0Ord[0] <= 0.0))
    fMFacOrd[0] = 0.;
  else
    fMFacOrd[0] = (0.3989422804014326*m0Ord[0]*(0.7071067811865475*bmag[0]-0.7071067811865474*bmag[1]))/sqrt(vtSqOrd[0])+9.999999999999999e-41; 
  if ((vtSqOrd[1] <= 0.0) || (m0Ord[1] <= 0.0))
    fMFacOrd[1] = 0.;
  else
    fMFacOrd[1] = (0.3989422804014326*(0.7071067811865474*bmag[1]+0.7071067811865475*bmag[0])*m0Ord[1])/sqrt(vtSqOrd[1])+9.999999999999999e-41; 

}

void GkMaxwellianOnBasisGauss1x1vSerUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*den[0]-0.7071067811865474*den[1]; 
  m0Ord[1] = 0.7071067811865474*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7071067811865475*flowU[4]-0.7071067811865474*flowU[5]; 
  flowUOrd[1] = 0.7071067811865474*flowU[5]+0.7071067811865475*flowU[4]; 

  vtSqOrd[0] = 0.7071067811865475*vtSq[0]-0.7071067811865474*vtSq[1]; 
  vtSqOrd[1] = 0.7071067811865474*vtSq[1]+0.7071067811865475*vtSq[0]; 

  if ((vtSqOrd[0] <= 0.0) || (m0Ord[0] <= 0.0))
    fMFacOrd[0] = 0.;
  else
    fMFacOrd[0] = (0.3989422804014326*m0Ord[0]*(0.7071067811865475*bmag[0]-0.7071067811865474*bmag[1]))/sqrt(vtSqOrd[0])+9.999999999999999e-41; 
  if ((vtSqOrd[1] <= 0.0) || (m0Ord[1] <= 0.0))
    fMFacOrd[1] = 0.;
  else
    fMFacOrd[1] = (0.3989422804014326*(0.7071067811865474*bmag[1]+0.7071067811865475*bmag[0])*m0Ord[1])/sqrt(vtSqOrd[1])+9.999999999999999e-41; 

}

void GkMaxwellianOnBasisGauss1x1vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[4];
  if ((vtSqOrd[0] <= 0.0) || (fMFacOrd[0] <= 0.0)) {
    fMquad[0] = 9.999999999999999e-41;
    fMquad[1] = 9.999999999999999e-41;
  } else {
    fMquad[0] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
    fMquad[1] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]+0.2886751345948129*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  };
  if ((vtSqOrd[1] <= 0.0) || (fMFacOrd[1] <= 0.0)) {
    fMquad[2] = 9.999999999999999e-41;
    fMquad[3] = 9.999999999999999e-41;
  } else {
    fMquad[2] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.2886751345948129*dxv[1],2.0))/vtSqOrd[1]); 
    fMquad[3] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.2886751345948129*dxv[1],2.0))/vtSqOrd[1]); 
  };

  fMOut[0] = 0.5*(fMquad[3]+fMquad[2]+fMquad[1]+fMquad[0]); 
  fMOut[1] = 0.4999999999999999*(fMquad[3]+fMquad[2])-0.4999999999999999*(fMquad[1]+fMquad[0]); 
  fMOut[2] = 0.4999999999999999*fMquad[3]-0.4999999999999999*fMquad[2]+0.4999999999999999*fMquad[1]-0.4999999999999999*fMquad[0]; 
  fMOut[3] = 0.4999999999999999*fMquad[3]-0.4999999999999999*(fMquad[2]+fMquad[1])+0.4999999999999999*fMquad[0]; 

}
