#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss1x1vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*(den[0]-den[1]); 
  m0Ord[1] = 0.7071067811865475*(den[1]+den[0]); 

  flowUOrd[0] = 0.7071067811865475*(flowU[0]-flowU[1]); 
  flowUOrd[1] = 0.7071067811865475*(flowU[1]+flowU[0]); 

  vtSqOrd[0] = 0.7071067811865475*(vtSq[0]-vtSq[1]); 
  vtSqOrd[1] = 0.7071067811865475*(vtSq[1]+vtSq[0]); 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (0.3989422804014326*m0Ord[0])/sqrt(vtSqOrd[0]); 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (0.3989422804014326*m0Ord[1])/sqrt(vtSqOrd[1]); 
  else
    fMFacOrd[1] = 0.0;

}

void MaxwellianOnBasisGauss1x1vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[6];
  if ((vtSqOrd[0] > 0.) && (fMFacOrd[0] > 0.)) {
    fMquad[0] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-0.3872983346207416*dxv[1]-flowUOrd[0],2.0))/vtSqOrd[0]); 
    fMquad[1] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-flowUOrd[0],2.0))/vtSqOrd[0]); 
    fMquad[2] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]+0.3872983346207416*dxv[1]-flowUOrd[0],2.0))/vtSqOrd[0]); 
  } else {
    fMquad[0] = 0.0;
    fMquad[1] = 0.0;
    fMquad[2] = 0.0;
  };
  if ((vtSqOrd[1] > 0.) && (fMFacOrd[1] > 0.)) {
    fMquad[3] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-flowUOrd[1]-0.3872983346207416*dxv[1],2.0))/vtSqOrd[1]); 
    fMquad[4] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-flowUOrd[1],2.0))/vtSqOrd[1]); 
    fMquad[5] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-flowUOrd[1]+0.3872983346207416*dxv[1],2.0))/vtSqOrd[1]); 
  } else {
    fMquad[3] = 0.0;
    fMquad[4] = 0.0;
    fMquad[5] = 0.0;
  };

  fMOut[0] = 0.2777777777777778*fMquad[5]+0.4444444444444444*fMquad[4]+0.2777777777777778*(fMquad[3]+fMquad[2])+0.4444444444444444*fMquad[1]+0.2777777777777778*fMquad[0]; 
  fMOut[1] = 1.732050807568877*(0.1603750747748961*fMquad[5]+0.2566001196398337*fMquad[4]+0.1603750747748961*fMquad[3]-0.1603750747748961*fMquad[2]-0.2566001196398337*fMquad[1]-0.1603750747748961*fMquad[0]); 
  fMOut[2] = 1.732050807568877*(0.2151657414559676*fMquad[5]-0.2151657414559676*fMquad[3]+0.2151657414559676*fMquad[2]-0.2151657414559676*fMquad[0]); 
  fMOut[3] = 0.3726779962499649*fMquad[5]-0.3726779962499649*(fMquad[3]+fMquad[2])+0.3726779962499649*fMquad[0]; 
  fMOut[4] = 2.23606797749979*(0.1111111111111111*fMquad[5]-0.2222222222222222*fMquad[4]+0.1111111111111111*(fMquad[3]+fMquad[2])-0.2222222222222222*fMquad[1]+0.1111111111111111*fMquad[0]); 
  fMOut[5] = 3.872983346207417*(0.0641500299099584*fMquad[5]-0.1283000598199168*fMquad[4]+0.0641500299099584*fMquad[3]-0.0641500299099584*fMquad[2]+0.1283000598199168*fMquad[1]-0.0641500299099584*fMquad[0]); 

}
void GkMaxwellianOnBasisGauss1x1vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*(den[0]-den[1]); 
  m0Ord[1] = 0.7071067811865475*(den[1]+den[0]); 

  flowUOrd[0] = 0.7071067811865475*(flowU[0]-flowU[1]); 
  flowUOrd[1] = 0.7071067811865475*(flowU[1]+flowU[0]); 

  vtSqOrd[0] = 0.7071067811865475*(vtSq[0]-vtSq[1]); 
  vtSqOrd[1] = 0.7071067811865475*(vtSq[1]+vtSq[0]); 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (0.3989422804014326*m0Ord[0]*(0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1]))/sqrt(vtSqOrd[0]); 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (0.2820947917738781*(bmag[1]+bmag[0])*m0Ord[1])/sqrt(vtSqOrd[1]); 
  else
    fMFacOrd[1] = 0.0;

}

void GkMaxwellianOnBasisGauss1x1vSerUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*(den[0]-den[1]); 
  m0Ord[1] = 0.7071067811865475*(den[1]+den[0]); 

  flowUOrd[0] = 0.7071067811865475*(flowU[0]-flowU[1]); 
  flowUOrd[1] = 0.7071067811865475*(flowU[1]+flowU[0]); 

  vtSqOrd[0] = 0.7071067811865475*(vtSq[0]-vtSq[1]); 
  vtSqOrd[1] = 0.7071067811865475*(vtSq[1]+vtSq[0]); 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (0.3989422804014326*m0Ord[0]*(0.7071067811865475*bmag[0]-0.7071067811865475*bmag[1]))/sqrt(vtSqOrd[0]); 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (0.2820947917738781*(bmag[1]+bmag[0])*m0Ord[1])/sqrt(vtSqOrd[1]); 
  else
    fMFacOrd[1] = 0.0;

}

void GkMaxwellianOnBasisGauss1x1vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[6];
  if ((vtSqOrd[0] > 0.) && (fMFacOrd[0] > 0.)) {
    fMquad[0] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-0.3872983346207416*dxv[1]-flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[1] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[2] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]+0.3872983346207416*dxv[1]-flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
  } else {
    fMquad[0] = 9.999999999999999e-41;
    fMquad[1] = 9.999999999999999e-41;
    fMquad[2] = 9.999999999999999e-41;
  };
  if ((vtSqOrd[1] > 0.) && (fMFacOrd[1] > 0.)) {
    fMquad[3] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-flowUOrd[1]-0.3872983346207416*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[4] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[5] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-flowUOrd[1]+0.3872983346207416*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
  } else {
    fMquad[3] = 9.999999999999999e-41;
    fMquad[4] = 9.999999999999999e-41;
    fMquad[5] = 9.999999999999999e-41;
  };

  fMOut[0] = 0.2777777777777778*fMquad[5]+0.4444444444444444*fMquad[4]+0.2777777777777778*(fMquad[3]+fMquad[2])+0.4444444444444444*fMquad[1]+0.2777777777777778*fMquad[0]; 
  fMOut[1] = 1.732050807568877*(0.1603750747748961*fMquad[5]+0.2566001196398337*fMquad[4]+0.1603750747748961*fMquad[3]-0.1603750747748961*fMquad[2]-0.2566001196398337*fMquad[1]-0.1603750747748961*fMquad[0]); 
  fMOut[2] = 1.732050807568877*(0.2151657414559676*fMquad[5]-0.2151657414559676*fMquad[3]+0.2151657414559676*fMquad[2]-0.2151657414559676*fMquad[0]); 
  fMOut[3] = 0.3726779962499649*fMquad[5]-0.3726779962499649*(fMquad[3]+fMquad[2])+0.3726779962499649*fMquad[0]; 
  fMOut[4] = 2.23606797749979*(0.1111111111111111*fMquad[5]-0.2222222222222222*fMquad[4]+0.1111111111111111*(fMquad[3]+fMquad[2])-0.2222222222222222*fMquad[1]+0.1111111111111111*fMquad[0]); 
  fMOut[5] = 3.872983346207417*(0.0641500299099584*fMquad[5]-0.1283000598199168*fMquad[4]+0.0641500299099584*fMquad[3]-0.0641500299099584*fMquad[2]+0.1283000598199168*fMquad[1]-0.0641500299099584*fMquad[0]); 

}
