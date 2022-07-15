#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss1x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*(den[0]-den[1]); 
  m0Ord[1] = 0.7071067811865475*(den[1]+den[0]); 

  flowUOrd[0] = 0.7071067811865475*(flowU[0]-flowU[1]); 
  flowUOrd[1] = 0.7071067811865475*(flowU[1]+flowU[0]); 
  flowUOrd[2] = 0.7071067811865475*(flowU[2]-flowU[3]); 
  flowUOrd[3] = 0.7071067811865475*(flowU[3]+flowU[2]); 

  vtSqOrd[0] = 0.7071067811865475*(vtSq[0]-vtSq[1]); 
  vtSqOrd[1] = 0.7071067811865475*(vtSq[1]+vtSq[0]); 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (0.1591549430918953*m0Ord[0])/vtSqOrd[0]; 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (0.1591549430918953*m0Ord[1])/vtSqOrd[1]; 
  else
    fMFacOrd[1] = 0.0;

}

void MaxwellianOnBasisGauss1x2vSerUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*(den[0]-den[1]); 
  m0Ord[1] = 0.7071067811865475*(den[1]+den[0]); 

  flowUOrd[0] = 0.7071067811865475*(flowU[0]-flowU[1]); 
  flowUOrd[1] = 0.7071067811865475*(flowU[1]+flowU[0]); 
  flowUOrd[2] = 0.0; 
  flowUOrd[3] = 0.0; 

  vtSqOrd[0] = 0.7071067811865475*(vtSq[0]-vtSq[1]); 
  vtSqOrd[1] = 0.7071067811865475*(vtSq[1]+vtSq[0]); 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (0.1591549430918953*m0Ord[0])/vtSqOrd[0]; 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (0.1591549430918953*m0Ord[1])/vtSqOrd[1]; 
  else
    fMFacOrd[1] = 0.0;

}

void MaxwellianOnBasisGauss1x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[18];
  if ((vtSqOrd[0] > 0.) && (fMFacOrd[0] > 0.)) {
    fMquad[0] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[2]-0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]-0.3872983346207416*dxv[1]-flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[1] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[2],2.0)+std::pow(wc[1]-0.3872983346207416*dxv[1]-flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[2] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[2]+0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]-0.3872983346207416*dxv[1]-flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[3] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[2]-0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]-flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[4] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[2],2.0)+std::pow(wc[1]-flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[5] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[2]+0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]-flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[6] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[2]-0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]+0.3872983346207416*dxv[1]-flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[7] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[2],2.0)+std::pow(wc[1]+0.3872983346207416*dxv[1]-flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[8] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[2]+0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]+0.3872983346207416*dxv[1]-flowUOrd[0],2.0)))/vtSqOrd[0]); 
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
    fMquad[9] = fMFacOrd[1]*exp(-(0.5*(std::pow((-flowUOrd[3])+wc[2]-0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]-flowUOrd[1]-0.3872983346207416*dxv[1],2.0)))/vtSqOrd[1]); 
    fMquad[10] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[3],2.0)+std::pow(wc[1]-flowUOrd[1]-0.3872983346207416*dxv[1],2.0)))/vtSqOrd[1]); 
    fMquad[11] = fMFacOrd[1]*exp(-(0.5*(std::pow((-flowUOrd[3])+wc[2]+0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]-flowUOrd[1]-0.3872983346207416*dxv[1],2.0)))/vtSqOrd[1]); 
    fMquad[12] = fMFacOrd[1]*exp(-(0.5*(std::pow((-flowUOrd[3])+wc[2]-0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]-flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[13] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[3],2.0)+std::pow(wc[1]-flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[14] = fMFacOrd[1]*exp(-(0.5*(std::pow((-flowUOrd[3])+wc[2]+0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]-flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[15] = fMFacOrd[1]*exp(-(0.5*(std::pow((-flowUOrd[3])+wc[2]-0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]-flowUOrd[1]+0.3872983346207416*dxv[1],2.0)))/vtSqOrd[1]); 
    fMquad[16] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[2]-flowUOrd[3],2.0)+std::pow(wc[1]-flowUOrd[1]+0.3872983346207416*dxv[1],2.0)))/vtSqOrd[1]); 
    fMquad[17] = fMFacOrd[1]*exp(-(0.5*(std::pow((-flowUOrd[3])+wc[2]+0.3872983346207416*dxv[2],2.0)+std::pow(wc[1]-flowUOrd[1]+0.3872983346207416*dxv[1],2.0)))/vtSqOrd[1]); 
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

  fMOut[0] = 0.3535533905932737*(0.308641975308642*fMquad[17]+0.4938271604938271*fMquad[16]+0.308641975308642*fMquad[15]+0.4938271604938271*fMquad[14]+0.7901234567901234*fMquad[13]+0.4938271604938271*fMquad[12]+0.308641975308642*fMquad[11]+0.4938271604938271*fMquad[10]+0.308641975308642*(fMquad[9]+fMquad[8])+0.4938271604938271*fMquad[7]+0.308641975308642*fMquad[6]+0.4938271604938271*fMquad[5]+0.7901234567901234*fMquad[4]+0.4938271604938271*fMquad[3]+0.308641975308642*fMquad[2]+0.4938271604938271*fMquad[1]+0.308641975308642*fMquad[0]); 
  fMOut[1] = 0.6123724356957944*(0.1781945275276623*fMquad[17]+0.2851112440442596*fMquad[16]+0.1781945275276623*fMquad[15]+0.2851112440442596*fMquad[14]+0.4561779904708155*fMquad[13]+0.2851112440442596*fMquad[12]+0.1781945275276623*fMquad[11]+0.2851112440442596*fMquad[10]+0.1781945275276623*fMquad[9]-0.1781945275276623*fMquad[8]-0.2851112440442596*fMquad[7]-0.1781945275276623*fMquad[6]-0.2851112440442596*fMquad[5]-0.4561779904708155*fMquad[4]-0.2851112440442596*fMquad[3]-0.1781945275276623*fMquad[2]-0.2851112440442596*fMquad[1]-0.1781945275276623*fMquad[0]); 
  fMOut[2] = 0.6123724356957944*(0.2390730460621862*fMquad[17]+0.3825168736994979*fMquad[16]+0.2390730460621862*fMquad[15]-0.2390730460621862*fMquad[11]-0.3825168736994979*fMquad[10]-0.2390730460621862*fMquad[9]+0.2390730460621862*fMquad[8]+0.3825168736994979*fMquad[7]+0.2390730460621862*fMquad[6]-0.2390730460621862*fMquad[2]-0.3825168736994979*fMquad[1]-0.2390730460621862*fMquad[0]); 
  fMOut[3] = 0.6123724356957944*(0.2390730460621862*fMquad[17]-0.2390730460621862*fMquad[15]+0.3825168736994979*fMquad[14]-0.3825168736994979*fMquad[12]+0.2390730460621862*fMquad[11]-0.2390730460621862*fMquad[9]+0.2390730460621862*fMquad[8]-0.2390730460621862*fMquad[6]+0.3825168736994979*fMquad[5]-0.3825168736994979*fMquad[3]+0.2390730460621862*fMquad[2]-0.2390730460621862*fMquad[0]); 
  fMOut[4] = 0.3535533905932737*(0.4140866624999611*fMquad[17]+0.6625386599999377*fMquad[16]+0.4140866624999611*fMquad[15]-0.4140866624999611*fMquad[11]-0.6625386599999377*fMquad[10]-0.4140866624999611*(fMquad[9]+fMquad[8])-0.6625386599999377*fMquad[7]-0.4140866624999611*fMquad[6]+0.4140866624999611*fMquad[2]+0.6625386599999377*fMquad[1]+0.4140866624999611*fMquad[0]); 
  fMOut[5] = 0.3535533905932737*(0.4140866624999611*fMquad[17]-0.4140866624999611*fMquad[15]+0.6625386599999377*fMquad[14]-0.6625386599999377*fMquad[12]+0.4140866624999611*fMquad[11]-0.4140866624999611*(fMquad[9]+fMquad[8])+0.4140866624999611*fMquad[6]-0.6625386599999377*fMquad[5]+0.6625386599999377*fMquad[3]-0.4140866624999611*fMquad[2]+0.4140866624999611*fMquad[0]); 
  fMOut[6] = 0.3535533905932737*(0.5555555555555555*fMquad[17]-0.5555555555555555*(fMquad[15]+fMquad[11])+0.5555555555555555*(fMquad[9]+fMquad[8])-0.5555555555555555*(fMquad[6]+fMquad[2])+0.5555555555555555*fMquad[0]); 
  fMOut[7] = 1.837117307087383*(0.1069167165165974*fMquad[17]-0.1069167165165974*(fMquad[15]+fMquad[11])+0.1069167165165974*fMquad[9]-0.1069167165165974*fMquad[8]+0.1069167165165974*(fMquad[6]+fMquad[2])-0.1069167165165974*fMquad[0]); 
  fMOut[8] = 0.3952847075210473*(0.2469135802469135*fMquad[17]+0.3950617283950615*fMquad[16]+0.2469135802469135*fMquad[15]-0.4938271604938271*fMquad[14]-0.7901234567901234*fMquad[13]-0.4938271604938271*fMquad[12]+0.2469135802469135*fMquad[11]+0.3950617283950615*fMquad[10]+0.2469135802469135*(fMquad[9]+fMquad[8])+0.3950617283950615*fMquad[7]+0.2469135802469135*fMquad[6]-0.4938271604938271*fMquad[5]-0.7901234567901234*fMquad[4]-0.4938271604938271*fMquad[3]+0.2469135802469135*fMquad[2]+0.3950617283950615*fMquad[1]+0.2469135802469135*fMquad[0]); 
  fMOut[9] = 0.6846531968814574*(0.1425556220221298*fMquad[17]+0.2280889952354076*fMquad[16]+0.1425556220221298*fMquad[15]-0.2851112440442596*fMquad[14]-0.4561779904708154*fMquad[13]-0.2851112440442596*fMquad[12]+0.1425556220221298*fMquad[11]+0.2280889952354076*fMquad[10]+0.1425556220221298*fMquad[9]-0.1425556220221298*fMquad[8]-0.2280889952354076*fMquad[7]-0.1425556220221298*fMquad[6]+0.2851112440442596*fMquad[5]+0.4561779904708154*fMquad[4]+0.2851112440442596*fMquad[3]-0.1425556220221298*fMquad[2]-0.2280889952354076*fMquad[1]-0.1425556220221298*fMquad[0]); 
  fMOut[10] = 0.6846531968814574*(0.1912584368497489*fMquad[17]-0.1912584368497489*fMquad[15]-0.3825168736994979*fMquad[14]+0.3825168736994979*fMquad[12]+0.1912584368497489*fMquad[11]-0.1912584368497489*fMquad[9]+0.1912584368497489*fMquad[8]-0.1912584368497489*fMquad[6]-0.3825168736994979*fMquad[5]+0.3825168736994979*fMquad[3]+0.1912584368497489*fMquad[2]-0.1912584368497489*fMquad[0]); 
  fMOut[11] = 0.3952847075210473*(0.3312693299999688*fMquad[17]-0.3312693299999688*fMquad[15]-0.6625386599999377*fMquad[14]+0.6625386599999377*fMquad[12]+0.3312693299999688*fMquad[11]-0.3312693299999688*(fMquad[9]+fMquad[8])+0.3312693299999688*fMquad[6]+0.6625386599999377*fMquad[5]-0.6625386599999377*fMquad[3]-0.3312693299999688*fMquad[2]+0.3312693299999688*fMquad[0]); 
  fMOut[12] = 0.3952847075210473*(0.2469135802469135*fMquad[17]-0.4938271604938271*fMquad[16]+0.2469135802469135*fMquad[15]+0.3950617283950615*fMquad[14]-0.7901234567901234*fMquad[13]+0.3950617283950615*fMquad[12]+0.2469135802469135*fMquad[11]-0.4938271604938271*fMquad[10]+0.2469135802469135*(fMquad[9]+fMquad[8])-0.4938271604938271*fMquad[7]+0.2469135802469135*fMquad[6]+0.3950617283950615*fMquad[5]-0.7901234567901234*fMquad[4]+0.3950617283950615*fMquad[3]+0.2469135802469135*fMquad[2]-0.4938271604938271*fMquad[1]+0.2469135802469135*fMquad[0]); 
  fMOut[13] = 0.6846531968814574*(0.1425556220221298*fMquad[17]-0.2851112440442596*fMquad[16]+0.1425556220221298*fMquad[15]+0.2280889952354076*fMquad[14]-0.4561779904708154*fMquad[13]+0.2280889952354076*fMquad[12]+0.1425556220221298*fMquad[11]-0.2851112440442596*fMquad[10]+0.1425556220221298*fMquad[9]-0.1425556220221298*fMquad[8]+0.2851112440442596*fMquad[7]-0.1425556220221298*fMquad[6]-0.2280889952354076*fMquad[5]+0.4561779904708154*fMquad[4]-0.2280889952354076*fMquad[3]-0.1425556220221298*fMquad[2]+0.2851112440442596*fMquad[1]-0.1425556220221298*fMquad[0]); 
  fMOut[14] = 0.6846531968814574*(0.1912584368497489*fMquad[17]-0.3825168736994979*fMquad[16]+0.1912584368497489*fMquad[15]-0.1912584368497489*fMquad[11]+0.3825168736994979*fMquad[10]-0.1912584368497489*fMquad[9]+0.1912584368497489*fMquad[8]-0.3825168736994979*fMquad[7]+0.1912584368497489*fMquad[6]-0.1912584368497489*fMquad[2]+0.3825168736994979*fMquad[1]-0.1912584368497489*fMquad[0]); 
  fMOut[15] = 0.3952847075210473*(0.3312693299999688*fMquad[17]-0.6625386599999376*fMquad[16]+0.3312693299999688*fMquad[15]-0.3312693299999688*fMquad[11]+0.6625386599999376*fMquad[10]-0.3312693299999688*(fMquad[9]+fMquad[8])+0.6625386599999376*fMquad[7]-0.3312693299999688*fMquad[6]+0.3312693299999688*fMquad[2]-0.6625386599999376*fMquad[1]+0.3312693299999688*fMquad[0]); 

}
void GkMaxwellianOnBasisGauss1x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*(den[0]-den[1]); 
  m0Ord[1] = 0.7071067811865475*(den[1]+den[0]); 

  flowUOrd[0] = 0.7071067811865475*(flowU[0]-flowU[1]); 
  flowUOrd[1] = 0.7071067811865475*(flowU[1]+flowU[0]); 

  vtSqOrd[0] = 0.7071067811865475*(vtSq[0]-vtSq[1]); 
  vtSqOrd[1] = 0.7071067811865475*(vtSq[1]+vtSq[0]); 

  bmagOrd[0] = 0.7071067811865475*(bmag[0]-bmag[1]); 
  bmagOrd[1] = 0.7071067811865475*(bmag[1]+bmag[0]); 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (bmagOrd[0]*m0Ord[0])/std::pow(2.506628274631001*sqrt(vtSqOrd[0]),3.0); 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (bmagOrd[1]*m0Ord[1])/std::pow(2.506628274631001*sqrt(vtSqOrd[1]),3.0); 
  else
    fMFacOrd[1] = 0.0;

}

void GkMaxwellianOnBasisGauss1x2vSerUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[2];
  m0Ord[0] = 0.7071067811865475*(den[0]-den[1]); 
  m0Ord[1] = 0.7071067811865475*(den[1]+den[0]); 

  flowUOrd[0] = 0.7071067811865475*(flowU[0]-flowU[1]); 
  flowUOrd[1] = 0.7071067811865475*(flowU[1]+flowU[0]); 

  vtSqOrd[0] = 0.7071067811865475*(vtSq[0]-vtSq[1]); 
  vtSqOrd[1] = 0.7071067811865475*(vtSq[1]+vtSq[0]); 

  bmagOrd[0] = 0.7071067811865475*(bmag[0]-bmag[1]); 
  bmagOrd[1] = 0.7071067811865475*(bmag[1]+bmag[0]); 

  if ((vtSqOrd[0] > 0.) && (m0Ord[0] > 0.))
    fMFacOrd[0] = (bmagOrd[0]*m0Ord[0])/std::pow(2.506628274631001*sqrt(vtSqOrd[0]),3.0); 
  else
    fMFacOrd[0] = 0.0;
  if ((vtSqOrd[1] > 0.) && (m0Ord[1] > 0.))
    fMFacOrd[1] = (bmagOrd[1]*m0Ord[1])/std::pow(2.506628274631001*sqrt(vtSqOrd[1]),3.0); 
  else
    fMFacOrd[1] = 0.0;

}

void GkMaxwellianOnBasisGauss1x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[12];
  if ((vtSqOrd[0] > 0.) && (fMFacOrd[0] > 0.)) {
    fMquad[0] = fMFacOrd[0]*exp(((-(bmagOrd[0]*std::abs(wc[2]-0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]-0.3872983346207416*dxv[1]-flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[1] = fMFacOrd[0]*exp(((-(bmagOrd[0]*std::abs(wc[2]-0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]-flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[2] = fMFacOrd[0]*exp(((-(bmagOrd[0]*std::abs(wc[2]-0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]+0.3872983346207416*dxv[1]-flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[3] = fMFacOrd[0]*exp(((-(bmagOrd[0]*std::abs(wc[2]+0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]-0.3872983346207416*dxv[1]-flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[4] = fMFacOrd[0]*exp(((-(bmagOrd[0]*std::abs(wc[2]+0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]-flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[5] = fMFacOrd[0]*exp(((-(bmagOrd[0]*std::abs(wc[2]+0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]+0.3872983346207416*dxv[1]-flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
  } else {
    fMquad[0] = 9.999999999999999e-41;
    fMquad[1] = 9.999999999999999e-41;
    fMquad[2] = 9.999999999999999e-41;
    fMquad[3] = 9.999999999999999e-41;
    fMquad[4] = 9.999999999999999e-41;
    fMquad[5] = 9.999999999999999e-41;
  };
  if ((vtSqOrd[1] > 0.) && (fMFacOrd[1] > 0.)) {
    fMquad[6] = fMFacOrd[1]*exp(((-(bmagOrd[1]*std::abs(wc[2]-0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]-flowUOrd[1]-0.3872983346207416*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[7] = fMFacOrd[1]*exp(((-(bmagOrd[1]*std::abs(wc[2]-0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]-flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[8] = fMFacOrd[1]*exp(((-(bmagOrd[1]*std::abs(wc[2]-0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]-flowUOrd[1]+0.3872983346207416*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[9] = fMFacOrd[1]*exp(((-(bmagOrd[1]*std::abs(wc[2]+0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]-flowUOrd[1]-0.3872983346207416*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[10] = fMFacOrd[1]*exp(((-(bmagOrd[1]*std::abs(wc[2]+0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]-flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[11] = fMFacOrd[1]*exp(((-(bmagOrd[1]*std::abs(wc[2]+0.2886751345948129*dxv[2]))/m_)-0.5*std::pow(wc[1]-flowUOrd[1]+0.3872983346207416*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
  } else {
    fMquad[6] = 9.999999999999999e-41;
    fMquad[7] = 9.999999999999999e-41;
    fMquad[8] = 9.999999999999999e-41;
    fMquad[9] = 9.999999999999999e-41;
    fMquad[10] = 9.999999999999999e-41;
    fMquad[11] = 9.999999999999999e-41;
  };

  fMOut[0] = 0.3535533905932737*(0.5555555555555556*fMquad[11]+0.8888888888888888*fMquad[10]+0.5555555555555556*(fMquad[9]+fMquad[8])+0.8888888888888888*fMquad[7]+0.5555555555555556*(fMquad[6]+fMquad[5])+0.8888888888888888*fMquad[4]+0.5555555555555556*(fMquad[3]+fMquad[2])+0.8888888888888888*fMquad[1]+0.5555555555555556*fMquad[0]); 
  fMOut[1] = 0.6123724356957944*(0.3207501495497921*fMquad[11]+0.5132002392796674*fMquad[10]+0.3207501495497921*(fMquad[9]+fMquad[8])+0.5132002392796674*fMquad[7]+0.3207501495497921*fMquad[6]-0.3207501495497921*fMquad[5]-0.5132002392796674*fMquad[4]-0.3207501495497921*(fMquad[3]+fMquad[2])-0.5132002392796674*fMquad[1]-0.3207501495497921*fMquad[0]); 
  fMOut[2] = 0.6123724356957944*(0.4303314829119352*fMquad[11]-0.4303314829119352*fMquad[9]+0.4303314829119352*fMquad[8]-0.4303314829119352*fMquad[6]+0.4303314829119352*fMquad[5]-0.4303314829119352*fMquad[3]+0.4303314829119352*fMquad[2]-0.4303314829119352*fMquad[0]); 
  fMOut[3] = 0.6123724356957944*(0.3207501495497921*fMquad[11]+0.5132002392796674*fMquad[10]+0.3207501495497921*fMquad[9]-0.3207501495497921*fMquad[8]-0.5132002392796674*fMquad[7]-0.3207501495497921*fMquad[6]+0.3207501495497921*fMquad[5]+0.5132002392796674*fMquad[4]+0.3207501495497921*fMquad[3]-0.3207501495497921*fMquad[2]-0.5132002392796674*fMquad[1]-0.3207501495497921*fMquad[0]); 
  fMOut[4] = 0.3535533905932737*(0.7453559924999299*fMquad[11]-0.7453559924999299*fMquad[9]+0.7453559924999299*fMquad[8]-0.7453559924999299*(fMquad[6]+fMquad[5])+0.7453559924999299*fMquad[3]-0.7453559924999299*fMquad[2]+0.7453559924999299*fMquad[0]); 
  fMOut[5] = 0.3535533905932737*(0.5555555555555557*fMquad[11]+0.8888888888888891*fMquad[10]+0.5555555555555557*fMquad[9]-0.5555555555555557*fMquad[8]-0.8888888888888891*fMquad[7]-0.5555555555555557*(fMquad[6]+fMquad[5])-0.8888888888888891*fMquad[4]-0.5555555555555557*fMquad[3]+0.5555555555555557*fMquad[2]+0.8888888888888891*fMquad[1]+0.5555555555555557*fMquad[0]); 
  fMOut[6] = 0.3535533905932737*(0.7453559924999299*fMquad[11]-0.7453559924999299*(fMquad[9]+fMquad[8])+0.7453559924999299*(fMquad[6]+fMquad[5])-0.7453559924999299*(fMquad[3]+fMquad[2])+0.7453559924999299*fMquad[0]); 
  fMOut[7] = 1.837117307087383*(0.1434438276373118*fMquad[11]-0.1434438276373118*(fMquad[9]+fMquad[8])+0.1434438276373118*fMquad[6]-0.1434438276373118*fMquad[5]+0.1434438276373118*(fMquad[3]+fMquad[2])-0.1434438276373118*fMquad[0]); 
  fMOut[8] = 0.3952847075210473*(0.4444444444444443*fMquad[11]-0.8888888888888888*fMquad[10]+0.4444444444444443*(fMquad[9]+fMquad[8])-0.8888888888888888*fMquad[7]+0.4444444444444443*(fMquad[6]+fMquad[5])-0.8888888888888888*fMquad[4]+0.4444444444444443*(fMquad[3]+fMquad[2])-0.8888888888888888*fMquad[1]+0.4444444444444443*fMquad[0]); 
  fMOut[9] = 0.6846531968814574*(0.2566001196398336*fMquad[11]-0.5132002392796673*fMquad[10]+0.2566001196398336*(fMquad[9]+fMquad[8])-0.5132002392796673*fMquad[7]+0.2566001196398336*fMquad[6]-0.2566001196398336*fMquad[5]+0.5132002392796673*fMquad[4]-0.2566001196398336*(fMquad[3]+fMquad[2])+0.5132002392796673*fMquad[1]-0.2566001196398336*fMquad[0]); 
  fMOut[10] = 0.6846531968814574*(0.2566001196398336*fMquad[11]-0.5132002392796674*fMquad[10]+0.2566001196398336*fMquad[9]-0.2566001196398336*fMquad[8]+0.5132002392796674*fMquad[7]-0.2566001196398336*fMquad[6]+0.2566001196398336*fMquad[5]-0.5132002392796674*fMquad[4]+0.2566001196398336*fMquad[3]-0.2566001196398336*fMquad[2]+0.5132002392796674*fMquad[1]-0.2566001196398336*fMquad[0]); 
  fMOut[11] = 0.3952847075210473*(0.4444444444444444*fMquad[11]-0.8888888888888891*fMquad[10]+0.4444444444444444*fMquad[9]-0.4444444444444444*fMquad[8]+0.8888888888888891*fMquad[7]-0.4444444444444444*(fMquad[6]+fMquad[5])+0.8888888888888891*fMquad[4]-0.4444444444444444*fMquad[3]+0.4444444444444444*fMquad[2]-0.8888888888888891*fMquad[1]+0.4444444444444444*fMquad[0]); 

}
