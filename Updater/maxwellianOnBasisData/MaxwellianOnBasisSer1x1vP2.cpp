#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss1x1vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, double *flowUOrd, double *vtSqOrd, double *fMFacOrd) {

  flowUOrd[0] = 0.7071067811865475*flowU[0]-0.7905694150420947*flowU[2]; 
  flowUOrd[1] = 0.6324555320336759*flowU[2]-0.9486832980505137*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[2] = 0.6324555320336759*flowU[2]+0.9486832980505137*flowU[1]+0.7071067811865475*flowU[0]; 

  vtSqOrd[0] = 0.7071067811865475*vtSq[0]-0.7905694150420947*vtSq[2]; 
  vtSqOrd[1] = 0.6324555320336759*vtSq[2]-0.9486832980505137*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[2] = 0.6324555320336759*vtSq[2]+0.9486832980505137*vtSq[1]+0.7071067811865475*vtSq[0]; 

  if (vtSqOrd[0] <= 0.0)
    fMFacOrd[0] = 0.;
  else
    fMFacOrd[0] = (0.3989422804014326*(0.7071067811865475*den[0]-0.7905694150420947*den[2]))/sqrt(vtSqOrd[0]); 
  if (vtSqOrd[1] <= 0.0)
    fMFacOrd[1] = 0.;
  else
    fMFacOrd[1] = (0.3989422804014326*(0.6324555320336759*den[2]-0.9486832980505137*den[1]+0.7071067811865475*den[0]))/sqrt(vtSqOrd[1]); 
  if (vtSqOrd[2] <= 0.0)
    fMFacOrd[2] = 0.;
  else
    fMFacOrd[2] = (0.3989422804014326*(0.6324555320336759*den[2]+0.9486832980505137*den[1]+0.7071067811865475*den[0]))/sqrt(vtSqOrd[2]); 

}

void MaxwellianOnBasisGauss1x1vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[9];
  fMquad[0] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[1] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[2] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]+0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  fMquad[3] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1],2.0))/vtSqOrd[1]); 
  fMquad[4] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.3872983346207417*dxv[1],2.0))/vtSqOrd[1]); 
  fMquad[5] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.3872983346207417*dxv[1],2.0))/vtSqOrd[1]); 
  fMquad[6] = fMFacOrd[2]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[2],2.0))/vtSqOrd[2]); 
  fMquad[7] = fMFacOrd[2]*exp(-(0.5*std::pow((-1.0*flowUOrd[2])+wc[1]-0.3872983346207417*dxv[1],2.0))/vtSqOrd[2]); 
  fMquad[8] = fMFacOrd[2]*exp(-(0.5*std::pow((-1.0*flowUOrd[2])+wc[1]+0.3872983346207417*dxv[1],2.0))/vtSqOrd[2]); 

  fMOut[0] = 0.154320987654321*(fMquad[8]+fMquad[7])+0.2469135802469136*fMquad[6]+0.154320987654321*(fMquad[5]+fMquad[4])+0.2469135802469136*(fMquad[3]+fMquad[2]+fMquad[1])+0.3950617283950617*fMquad[0]; 
  fMOut[1] = 0.2070433312499806*(fMquad[8]+fMquad[7])+0.3312693299999688*fMquad[6]-0.2070433312499806*(fMquad[5]+fMquad[4])-0.3312693299999688*fMquad[3]; 
  fMOut[2] = 0.2070433312499806*fMquad[8]-0.2070433312499806*fMquad[7]+0.2070433312499806*fMquad[5]-0.2070433312499806*fMquad[4]+0.3312693299999688*fMquad[2]-0.3312693299999688*fMquad[1]; 
  fMOut[3] = 0.2777777777777778*fMquad[8]-0.2777777777777778*(fMquad[7]+fMquad[5])+0.2777777777777778*fMquad[4]; 
  fMOut[4] = 0.1380288874999871*(fMquad[8]+fMquad[7])+0.2208462199999793*fMquad[6]+0.1380288874999871*(fMquad[5]+fMquad[4])+0.2208462199999793*fMquad[3]-0.276057774999974*(fMquad[2]+fMquad[1])-0.4416924399999584*fMquad[0]; 
  fMOut[5] = 0.1380288874999871*(fMquad[8]+fMquad[7])-0.276057774999974*fMquad[6]+0.1380288874999871*(fMquad[5]+fMquad[4])-0.276057774999974*fMquad[3]+0.2208462199999793*(fMquad[2]+fMquad[1])-0.4416924399999584*fMquad[0]; 
  fMOut[6] = 0.1851851851851853*fMquad[8]-0.1851851851851853*fMquad[7]+0.1851851851851853*fMquad[5]-0.1851851851851853*fMquad[4]-0.3703703703703703*fMquad[2]+0.3703703703703703*fMquad[1]; 
  fMOut[7] = 0.1851851851851853*(fMquad[8]+fMquad[7])-0.3703703703703703*fMquad[6]-0.1851851851851853*(fMquad[5]+fMquad[4])+0.3703703703703703*fMquad[3]; 

}
