#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss1x2vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[3];
  m0Ord[0] = 0.7071067811865475*den[0]-0.7905694150420947*den[2]; 
  m0Ord[1] = 0.6324555320336759*den[2]-0.9486832980505137*den[1]+0.7071067811865475*den[0]; 
  m0Ord[2] = 0.6324555320336759*den[2]+0.9486832980505137*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7071067811865475*flowU[0]-0.7905694150420947*flowU[2]; 
  flowUOrd[1] = 0.6324555320336759*flowU[2]-0.9486832980505137*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[2] = 0.6324555320336759*flowU[2]+0.9486832980505137*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[3] = 0.7071067811865475*flowU[3]-0.7905694150420947*flowU[5]; 
  flowUOrd[4] = 0.6324555320336759*flowU[5]-0.9486832980505137*flowU[4]+0.7071067811865475*flowU[3]; 
  flowUOrd[5] = 0.6324555320336759*flowU[5]+0.9486832980505137*flowU[4]+0.7071067811865475*flowU[3]; 

  vtSqOrd[0] = 0.7071067811865475*vtSq[0]-0.7905694150420947*vtSq[2]; 
  vtSqOrd[1] = 0.6324555320336759*vtSq[2]-0.9486832980505137*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[2] = 0.6324555320336759*vtSq[2]+0.9486832980505137*vtSq[1]+0.7071067811865475*vtSq[0]; 

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

}

void MaxwellianOnBasisGauss1x2vSerUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[3];
  m0Ord[0] = 0.7071067811865475*den[0]-0.7905694150420947*den[2]; 
  m0Ord[1] = 0.6324555320336759*den[2]-0.9486832980505137*den[1]+0.7071067811865475*den[0]; 
  m0Ord[2] = 0.6324555320336759*den[2]+0.9486832980505137*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7071067811865475*flowU[0]-0.7905694150420947*flowU[2]; 
  flowUOrd[1] = 0.6324555320336759*flowU[2]-0.9486832980505137*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[2] = 0.6324555320336759*flowU[2]+0.9486832980505137*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[3] = 0.0; 
  flowUOrd[4] = 0.0; 
  flowUOrd[5] = 0.0; 

  vtSqOrd[0] = 0.7071067811865475*vtSq[0]-0.7905694150420947*vtSq[2]; 
  vtSqOrd[1] = 0.6324555320336759*vtSq[2]-0.9486832980505137*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[2] = 0.6324555320336759*vtSq[2]+0.9486832980505137*vtSq[1]+0.7071067811865475*vtSq[0]; 

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

}

void MaxwellianOnBasisGauss1x2vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[27];
  if ((vtSqOrd[0] > 0.) && (fMFacOrd[0] > 0.)) {
    fMquad[0] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-1.0*flowUOrd[3],2.0)+std::pow(wc[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[1] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[2] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[3] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-1.0*flowUOrd[3],2.0)+std::pow(wc[1]-0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[4] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[5] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[6] = fMFacOrd[0]*exp(-(0.5*(std::pow(wc[2]-1.0*flowUOrd[3],2.0)+std::pow(wc[1]+0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[7] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[3])+wc[2]-0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]+0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
    fMquad[8] = fMFacOrd[0]*exp(-(0.5*(std::pow((-1.0*flowUOrd[3])+wc[2]+0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]+0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0)))/vtSqOrd[0]); 
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
    fMquad[9] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[2]-1.0*flowUOrd[4],2.0)+std::pow(wc[1]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[10] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[11] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1],2.0)))/vtSqOrd[1]); 
    fMquad[12] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[2]-1.0*flowUOrd[4],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]-0.3872983346207417*dxv[1],2.0)))/vtSqOrd[1]); 
    fMquad[13] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]-0.3872983346207417*dxv[1],2.0)))/vtSqOrd[1]); 
    fMquad[14] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]-0.3872983346207417*dxv[1],2.0)))/vtSqOrd[1]); 
    fMquad[15] = fMFacOrd[1]*exp(-(0.5*(std::pow(wc[2]-1.0*flowUOrd[4],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]+0.3872983346207417*dxv[1],2.0)))/vtSqOrd[1]); 
    fMquad[16] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[2]-0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]+0.3872983346207417*dxv[1],2.0)))/vtSqOrd[1]); 
    fMquad[17] = fMFacOrd[1]*exp(-(0.5*(std::pow((-1.0*flowUOrd[4])+wc[2]+0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[1]+0.3872983346207417*dxv[1],2.0)))/vtSqOrd[1]); 
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
    fMquad[18] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[2]-1.0*flowUOrd[5],2.0)+std::pow(wc[1]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
    fMquad[19] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
    fMquad[20] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)+std::pow(wc[1]-1.0*flowUOrd[2],2.0)))/vtSqOrd[2]); 
    fMquad[21] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[2]-1.0*flowUOrd[5],2.0)+std::pow((-1.0*flowUOrd[2])+wc[1]-0.3872983346207417*dxv[1],2.0)))/vtSqOrd[2]); 
    fMquad[22] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)+std::pow((-1.0*flowUOrd[2])+wc[1]-0.3872983346207417*dxv[1],2.0)))/vtSqOrd[2]); 
    fMquad[23] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)+std::pow((-1.0*flowUOrd[2])+wc[1]-0.3872983346207417*dxv[1],2.0)))/vtSqOrd[2]); 
    fMquad[24] = fMFacOrd[2]*exp(-(0.5*(std::pow(wc[2]-1.0*flowUOrd[5],2.0)+std::pow((-1.0*flowUOrd[2])+wc[1]+0.3872983346207417*dxv[1],2.0)))/vtSqOrd[2]); 
    fMquad[25] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[2]-0.3872983346207417*dxv[2],2.0)+std::pow((-1.0*flowUOrd[2])+wc[1]+0.3872983346207417*dxv[1],2.0)))/vtSqOrd[2]); 
    fMquad[26] = fMFacOrd[2]*exp(-(0.5*(std::pow((-1.0*flowUOrd[5])+wc[2]+0.3872983346207417*dxv[2],2.0)+std::pow((-1.0*flowUOrd[2])+wc[1]+0.3872983346207417*dxv[1],2.0)))/vtSqOrd[2]); 
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

  fMOut[0] = 0.06062300936098654*(fMquad[26]+fMquad[25])+0.09699681497757848*fMquad[24]+0.06062300936098654*(fMquad[23]+fMquad[22])+0.09699681497757848*(fMquad[21]+fMquad[20]+fMquad[19])+0.1551949039641256*fMquad[18]+0.06062300936098654*(fMquad[17]+fMquad[16])+0.09699681497757848*fMquad[15]+0.06062300936098654*(fMquad[14]+fMquad[13])+0.09699681497757848*(fMquad[12]+fMquad[11]+fMquad[10])+0.1551949039641256*fMquad[9]+0.09699681497757848*(fMquad[8]+fMquad[7])+0.1551949039641256*fMquad[6]+0.09699681497757848*(fMquad[5]+fMquad[4])+0.1551949039641256*(fMquad[3]+fMquad[2]+fMquad[1])+0.248311846342601*fMquad[0]; 
  fMOut[1] = 0.0813343019590632*(fMquad[26]+fMquad[25])+0.1301348831345011*fMquad[24]+0.0813343019590632*(fMquad[23]+fMquad[22])+0.1301348831345011*(fMquad[21]+fMquad[20]+fMquad[19])+0.2082158130152018*fMquad[18]-0.0813343019590632*(fMquad[17]+fMquad[16])-0.1301348831345011*fMquad[15]-0.0813343019590632*(fMquad[14]+fMquad[13])-0.1301348831345011*(fMquad[12]+fMquad[11]+fMquad[10])-0.2082158130152018*fMquad[9]; 
  fMOut[2] = 0.0813343019590632*(fMquad[26]+fMquad[25])+0.1301348831345011*fMquad[24]-0.0813343019590632*(fMquad[23]+fMquad[22])-0.1301348831345011*fMquad[21]+0.0813343019590632*(fMquad[17]+fMquad[16])+0.1301348831345011*fMquad[15]-0.0813343019590632*(fMquad[14]+fMquad[13])-0.1301348831345011*fMquad[12]+0.1301348831345011*(fMquad[8]+fMquad[7])+0.2082158130152018*fMquad[6]-0.1301348831345011*(fMquad[5]+fMquad[4])-0.2082158130152018*fMquad[3]; 
  fMOut[3] = 0.0813343019590632*fMquad[26]-0.0813343019590632*fMquad[25]+0.0813343019590632*fMquad[23]-0.0813343019590632*fMquad[22]+0.1301348831345011*fMquad[20]-0.1301348831345011*fMquad[19]+0.0813343019590632*fMquad[17]-0.0813343019590632*fMquad[16]+0.0813343019590632*fMquad[14]-0.0813343019590632*fMquad[13]+0.1301348831345011*fMquad[11]-0.1301348831345011*fMquad[10]+0.1301348831345011*fMquad[8]-0.1301348831345011*fMquad[7]+0.1301348831345011*fMquad[5]-0.1301348831345011*fMquad[4]+0.2082158130152018*fMquad[2]-0.2082158130152018*fMquad[1]; 
  fMOut[4] = 0.1091214168497758*(fMquad[26]+fMquad[25])+0.1745942669596413*fMquad[24]-0.1091214168497758*(fMquad[23]+fMquad[22])-0.1745942669596413*fMquad[21]-0.1091214168497758*(fMquad[17]+fMquad[16])-0.1745942669596413*fMquad[15]+0.1091214168497758*(fMquad[14]+fMquad[13])+0.1745942669596413*fMquad[12]; 
  fMOut[5] = 0.1091214168497758*fMquad[26]-0.1091214168497758*fMquad[25]+0.1091214168497758*fMquad[23]-0.1091214168497758*fMquad[22]+0.1745942669596413*fMquad[20]-0.1745942669596413*fMquad[19]-0.1091214168497758*fMquad[17]+0.1091214168497758*fMquad[16]-0.1091214168497758*fMquad[14]+0.1091214168497758*fMquad[13]-0.1745942669596413*fMquad[11]+0.1745942669596413*fMquad[10]; 
  fMOut[6] = 0.1091214168497758*fMquad[26]-0.1091214168497758*(fMquad[25]+fMquad[23])+0.1091214168497758*(fMquad[22]+fMquad[17])-0.1091214168497758*(fMquad[16]+fMquad[14])+0.1091214168497758*fMquad[13]+0.1745942669596413*fMquad[8]-0.1745942669596413*(fMquad[7]+fMquad[5])+0.1745942669596413*fMquad[4]; 
  fMOut[7] = 0.0542228679727088*(fMquad[26]+fMquad[25])+0.0867565887563341*fMquad[24]+0.0542228679727088*(fMquad[23]+fMquad[22])+0.0867565887563341*(fMquad[21]+fMquad[20]+fMquad[19])+0.1388105420101346*fMquad[18]+0.0542228679727088*(fMquad[17]+fMquad[16])+0.0867565887563341*fMquad[15]+0.0542228679727088*(fMquad[14]+fMquad[13])+0.0867565887563341*(fMquad[12]+fMquad[11]+fMquad[10])+0.1388105420101346*fMquad[9]-0.1084457359454176*(fMquad[8]+fMquad[7])-0.1735131775126682*fMquad[6]-0.1084457359454176*(fMquad[5]+fMquad[4])-0.1735131775126682*(fMquad[3]+fMquad[2]+fMquad[1])-0.2776210840202691*fMquad[0]; 
  fMOut[8] = 0.0542228679727088*(fMquad[26]+fMquad[25])+0.0867565887563341*fMquad[24]+0.0542228679727088*(fMquad[23]+fMquad[22])+0.0867565887563341*fMquad[21]-0.1084457359454176*(fMquad[20]+fMquad[19])-0.1735131775126682*fMquad[18]+0.0542228679727088*(fMquad[17]+fMquad[16])+0.0867565887563341*fMquad[15]+0.0542228679727088*(fMquad[14]+fMquad[13])+0.0867565887563341*fMquad[12]-0.1084457359454176*(fMquad[11]+fMquad[10])-0.1735131775126682*fMquad[9]+0.0867565887563341*(fMquad[8]+fMquad[7])+0.1388105420101346*fMquad[6]+0.0867565887563341*(fMquad[5]+fMquad[4])+0.1388105420101346*fMquad[3]-0.1735131775126682*(fMquad[2]+fMquad[1])-0.2776210840202691*fMquad[0]; 
  fMOut[9] = 0.0542228679727088*(fMquad[26]+fMquad[25])-0.1084457359454176*fMquad[24]+0.0542228679727088*(fMquad[23]+fMquad[22])-0.1084457359454176*fMquad[21]+0.0867565887563341*(fMquad[20]+fMquad[19])-0.1735131775126682*fMquad[18]+0.0542228679727088*(fMquad[17]+fMquad[16])-0.1084457359454176*fMquad[15]+0.0542228679727088*(fMquad[14]+fMquad[13])-0.1084457359454176*fMquad[12]+0.0867565887563341*(fMquad[11]+fMquad[10])-0.1735131775126682*fMquad[9]+0.0867565887563341*(fMquad[8]+fMquad[7])-0.1735131775126682*fMquad[6]+0.0867565887563341*(fMquad[5]+fMquad[4])-0.1735131775126682*fMquad[3]+0.1388105420101346*(fMquad[2]+fMquad[1])-0.2776210840202691*fMquad[0]; 
  fMOut[10] = 0.1464017435263137*fMquad[26]-0.1464017435263137*(fMquad[25]+fMquad[23])+0.1464017435263137*fMquad[22]-0.1464017435263137*fMquad[17]+0.1464017435263137*(fMquad[16]+fMquad[14])-0.1464017435263137*fMquad[13]; 
  fMOut[11] = 0.07274761123318388*(fMquad[26]+fMquad[25])+0.1163961779730942*fMquad[24]-0.07274761123318388*(fMquad[23]+fMquad[22])-0.1163961779730942*fMquad[21]+0.07274761123318388*(fMquad[17]+fMquad[16])+0.1163961779730942*fMquad[15]-0.07274761123318388*(fMquad[14]+fMquad[13])-0.1163961779730942*fMquad[12]-0.1454952224663677*(fMquad[8]+fMquad[7])-0.2327923559461884*fMquad[6]+0.1454952224663677*(fMquad[5]+fMquad[4])+0.2327923559461884*fMquad[3]; 
  fMOut[12] = 0.07274761123318388*(fMquad[26]+fMquad[25])+0.1163961779730942*fMquad[24]+0.07274761123318388*(fMquad[23]+fMquad[22])+0.1163961779730942*fMquad[21]-0.1454952224663677*(fMquad[20]+fMquad[19])-0.2327923559461884*fMquad[18]-0.07274761123318388*(fMquad[17]+fMquad[16])-0.1163961779730942*fMquad[15]-0.07274761123318388*(fMquad[14]+fMquad[13])-0.1163961779730942*fMquad[12]+0.1454952224663677*(fMquad[11]+fMquad[10])+0.2327923559461884*fMquad[9]; 
  fMOut[13] = 0.07274761123318388*fMquad[26]-0.07274761123318388*fMquad[25]+0.07274761123318388*fMquad[23]-0.07274761123318388*fMquad[22]+0.1163961779730942*fMquad[20]-0.1163961779730942*fMquad[19]+0.07274761123318388*fMquad[17]-0.07274761123318388*fMquad[16]+0.07274761123318388*fMquad[14]-0.07274761123318388*fMquad[13]+0.1163961779730942*fMquad[11]-0.1163961779730942*fMquad[10]-0.1454952224663677*fMquad[8]+0.1454952224663677*fMquad[7]-0.1454952224663677*fMquad[5]+0.1454952224663677*fMquad[4]-0.2327923559461884*fMquad[2]+0.2327923559461884*fMquad[1]; 
  fMOut[14] = 0.07274761123318388*fMquad[26]-0.07274761123318388*fMquad[25]+0.07274761123318388*fMquad[23]-0.07274761123318388*fMquad[22]-0.1454952224663677*fMquad[20]+0.1454952224663677*fMquad[19]+0.07274761123318388*fMquad[17]-0.07274761123318388*fMquad[16]+0.07274761123318388*fMquad[14]-0.07274761123318388*fMquad[13]-0.1454952224663677*fMquad[11]+0.1454952224663677*fMquad[10]+0.1163961779730942*fMquad[8]-0.1163961779730942*fMquad[7]+0.1163961779730942*fMquad[5]-0.1163961779730942*fMquad[4]-0.2327923559461884*fMquad[2]+0.2327923559461884*fMquad[1]; 
  fMOut[15] = 0.07274761123318388*(fMquad[26]+fMquad[25])-0.1454952224663677*fMquad[24]+0.07274761123318388*(fMquad[23]+fMquad[22])-0.1454952224663677*fMquad[21]+0.1163961779730942*(fMquad[20]+fMquad[19])-0.2327923559461884*fMquad[18]-0.07274761123318388*(fMquad[17]+fMquad[16])+0.1454952224663677*fMquad[15]-0.07274761123318388*(fMquad[14]+fMquad[13])+0.1454952224663677*fMquad[12]-0.1163961779730942*(fMquad[11]+fMquad[10])+0.2327923559461884*fMquad[9]; 
  fMOut[16] = 0.07274761123318388*(fMquad[26]+fMquad[25])-0.1454952224663677*fMquad[24]-0.07274761123318388*(fMquad[23]+fMquad[22])+0.1454952224663677*fMquad[21]+0.07274761123318388*(fMquad[17]+fMquad[16])-0.1454952224663677*fMquad[15]-0.07274761123318388*(fMquad[14]+fMquad[13])+0.1454952224663677*fMquad[12]+0.1163961779730942*(fMquad[8]+fMquad[7])-0.2327923559461884*fMquad[6]-0.1163961779730942*(fMquad[5]+fMquad[4])+0.2327923559461884*fMquad[3]; 
  fMOut[17] = 0.09760116235087586*fMquad[26]-0.09760116235087586*(fMquad[25]+fMquad[23])+0.09760116235087586*(fMquad[22]+fMquad[17])-0.09760116235087586*(fMquad[16]+fMquad[14])+0.09760116235087586*fMquad[13]-0.1952023247017517*fMquad[8]+0.1952023247017517*(fMquad[7]+fMquad[5])-0.1952023247017517*fMquad[4]; 
  fMOut[18] = 0.09760116235087589*fMquad[26]-0.09760116235087589*fMquad[25]+0.09760116235087589*fMquad[23]-0.09760116235087589*fMquad[22]-0.1952023247017516*fMquad[20]+0.1952023247017516*fMquad[19]-0.09760116235087589*fMquad[17]+0.09760116235087589*fMquad[16]-0.09760116235087589*fMquad[14]+0.09760116235087589*fMquad[13]+0.1952023247017516*fMquad[11]-0.1952023247017516*fMquad[10]; 
  fMOut[19] = 0.09760116235087589*(fMquad[26]+fMquad[25])-0.1952023247017516*fMquad[24]-0.09760116235087589*(fMquad[23]+fMquad[22])+0.1952023247017516*fMquad[21]-0.09760116235087589*(fMquad[17]+fMquad[16])+0.1952023247017516*fMquad[15]+0.09760116235087589*(fMquad[14]+fMquad[13])-0.1952023247017516*fMquad[12]; 

}
void GkMaxwellianOnBasisGauss1x2vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[3];
  m0Ord[0] = 0.7071067811865475*den[0]-0.7905694150420947*den[2]; 
  m0Ord[1] = 0.6324555320336759*den[2]-0.9486832980505137*den[1]+0.7071067811865475*den[0]; 
  m0Ord[2] = 0.6324555320336759*den[2]+0.9486832980505137*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7071067811865475*flowU[0]-0.7905694150420947*flowU[2]; 
  flowUOrd[1] = 0.6324555320336759*flowU[2]-0.9486832980505137*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[2] = 0.6324555320336759*flowU[2]+0.9486832980505137*flowU[1]+0.7071067811865475*flowU[0]; 

  vtSqOrd[0] = 0.7071067811865475*vtSq[0]-0.7905694150420947*vtSq[2]; 
  vtSqOrd[1] = 0.6324555320336759*vtSq[2]-0.9486832980505137*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[2] = 0.6324555320336759*vtSq[2]+0.9486832980505137*vtSq[1]+0.7071067811865475*vtSq[0]; 

  bmagOrd[0] = 0.7071067811865475*bmag[0]-0.7905694150420947*bmag[2]; 
  bmagOrd[1] = 0.6324555320336759*bmag[2]-0.9486832980505137*bmag[1]+0.7071067811865475*bmag[0]; 
  bmagOrd[2] = 0.6324555320336759*bmag[2]+0.9486832980505137*bmag[1]+0.7071067811865475*bmag[0]; 

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

}

void GkMaxwellianOnBasisGauss1x2vSerUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[3];
  m0Ord[0] = 0.7071067811865475*den[0]-0.7905694150420947*den[2]; 
  m0Ord[1] = 0.6324555320336759*den[2]-0.9486832980505137*den[1]+0.7071067811865475*den[0]; 
  m0Ord[2] = 0.6324555320336759*den[2]+0.9486832980505137*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7071067811865475*flowU[0]-0.7905694150420947*flowU[2]; 
  flowUOrd[1] = 0.6324555320336759*flowU[2]-0.9486832980505137*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[2] = 0.6324555320336759*flowU[2]+0.9486832980505137*flowU[1]+0.7071067811865475*flowU[0]; 

  vtSqOrd[0] = 0.7071067811865475*vtSq[0]-0.7905694150420947*vtSq[2]; 
  vtSqOrd[1] = 0.6324555320336759*vtSq[2]-0.9486832980505137*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[2] = 0.6324555320336759*vtSq[2]+0.9486832980505137*vtSq[1]+0.7071067811865475*vtSq[0]; 

  bmagOrd[0] = 0.7071067811865475*bmag[0]-0.7905694150420947*bmag[2]; 
  bmagOrd[1] = 0.6324555320336759*bmag[2]-0.9486832980505137*bmag[1]+0.7071067811865475*bmag[0]; 
  bmagOrd[2] = 0.6324555320336759*bmag[2]+0.9486832980505137*bmag[1]+0.7071067811865475*bmag[0]; 

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

}

void GkMaxwellianOnBasisGauss1x2vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[27];
  if ((vtSqOrd[0] > 0.) && (fMFacOrd[0] > 0.)) {
    fMquad[0] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[1] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[2]-0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[2] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[2]+0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[3] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[2]))/m_)-0.5*std::pow(wc[1]-0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[4] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[2]-0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[5] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[2]+0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[6] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[2]))/m_)-0.5*std::pow(wc[1]+0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[7] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[2]-0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]+0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
    fMquad[8] = fMFacOrd[0]*exp(((-(1.0*bmagOrd[0]*std::abs(wc[2]+0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]+0.3872983346207417*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0])+9.999999999999999e-41; 
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
    fMquad[9] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[10] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[2]-0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[11] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[2]+0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[12] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.3872983346207417*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[13] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[2]-0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.3872983346207417*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[14] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[2]+0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.3872983346207417*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[15] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.3872983346207417*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[16] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[2]-0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.3872983346207417*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
    fMquad[17] = fMFacOrd[1]*exp(((-(1.0*bmagOrd[1]*std::abs(wc[2]+0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.3872983346207417*dxv[1],2.0))/vtSqOrd[1])+9.999999999999999e-41; 
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
    fMquad[18] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[19] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[2]-0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[20] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[2]+0.3872983346207417*dxv[2]))/m_)-0.5*std::pow(wc[1]-1.0*flowUOrd[2],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[21] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[2]))/m_)-0.5*std::pow((-1.0*flowUOrd[2])+wc[1]-0.3872983346207417*dxv[1],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[22] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[2]-0.3872983346207417*dxv[2]))/m_)-0.5*std::pow((-1.0*flowUOrd[2])+wc[1]-0.3872983346207417*dxv[1],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[23] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[2]+0.3872983346207417*dxv[2]))/m_)-0.5*std::pow((-1.0*flowUOrd[2])+wc[1]-0.3872983346207417*dxv[1],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[24] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[2]))/m_)-0.5*std::pow((-1.0*flowUOrd[2])+wc[1]+0.3872983346207417*dxv[1],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[25] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[2]-0.3872983346207417*dxv[2]))/m_)-0.5*std::pow((-1.0*flowUOrd[2])+wc[1]+0.3872983346207417*dxv[1],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
    fMquad[26] = fMFacOrd[2]*exp(((-(1.0*bmagOrd[2]*std::abs(wc[2]+0.3872983346207417*dxv[2]))/m_)-0.5*std::pow((-1.0*flowUOrd[2])+wc[1]+0.3872983346207417*dxv[1],2.0))/vtSqOrd[2])+9.999999999999999e-41; 
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

  fMOut[0] = 0.06062300936098654*(fMquad[26]+fMquad[25])+0.09699681497757848*fMquad[24]+0.06062300936098654*(fMquad[23]+fMquad[22])+0.09699681497757848*(fMquad[21]+fMquad[20]+fMquad[19])+0.1551949039641256*fMquad[18]+0.06062300936098654*(fMquad[17]+fMquad[16])+0.09699681497757848*fMquad[15]+0.06062300936098654*(fMquad[14]+fMquad[13])+0.09699681497757848*(fMquad[12]+fMquad[11]+fMquad[10])+0.1551949039641256*fMquad[9]+0.09699681497757848*(fMquad[8]+fMquad[7])+0.1551949039641256*fMquad[6]+0.09699681497757848*(fMquad[5]+fMquad[4])+0.1551949039641256*(fMquad[3]+fMquad[2]+fMquad[1])+0.248311846342601*fMquad[0]; 
  fMOut[1] = 0.0813343019590632*(fMquad[26]+fMquad[25])+0.1301348831345011*fMquad[24]+0.0813343019590632*(fMquad[23]+fMquad[22])+0.1301348831345011*(fMquad[21]+fMquad[20]+fMquad[19])+0.2082158130152018*fMquad[18]-0.0813343019590632*(fMquad[17]+fMquad[16])-0.1301348831345011*fMquad[15]-0.0813343019590632*(fMquad[14]+fMquad[13])-0.1301348831345011*(fMquad[12]+fMquad[11]+fMquad[10])-0.2082158130152018*fMquad[9]; 
  fMOut[2] = 0.0813343019590632*(fMquad[26]+fMquad[25])+0.1301348831345011*fMquad[24]-0.0813343019590632*(fMquad[23]+fMquad[22])-0.1301348831345011*fMquad[21]+0.0813343019590632*(fMquad[17]+fMquad[16])+0.1301348831345011*fMquad[15]-0.0813343019590632*(fMquad[14]+fMquad[13])-0.1301348831345011*fMquad[12]+0.1301348831345011*(fMquad[8]+fMquad[7])+0.2082158130152018*fMquad[6]-0.1301348831345011*(fMquad[5]+fMquad[4])-0.2082158130152018*fMquad[3]; 
  fMOut[3] = 0.0813343019590632*fMquad[26]-0.0813343019590632*fMquad[25]+0.0813343019590632*fMquad[23]-0.0813343019590632*fMquad[22]+0.1301348831345011*fMquad[20]-0.1301348831345011*fMquad[19]+0.0813343019590632*fMquad[17]-0.0813343019590632*fMquad[16]+0.0813343019590632*fMquad[14]-0.0813343019590632*fMquad[13]+0.1301348831345011*fMquad[11]-0.1301348831345011*fMquad[10]+0.1301348831345011*fMquad[8]-0.1301348831345011*fMquad[7]+0.1301348831345011*fMquad[5]-0.1301348831345011*fMquad[4]+0.2082158130152018*fMquad[2]-0.2082158130152018*fMquad[1]; 
  fMOut[4] = 0.1091214168497758*(fMquad[26]+fMquad[25])+0.1745942669596413*fMquad[24]-0.1091214168497758*(fMquad[23]+fMquad[22])-0.1745942669596413*fMquad[21]-0.1091214168497758*(fMquad[17]+fMquad[16])-0.1745942669596413*fMquad[15]+0.1091214168497758*(fMquad[14]+fMquad[13])+0.1745942669596413*fMquad[12]; 
  fMOut[5] = 0.1091214168497758*fMquad[26]-0.1091214168497758*fMquad[25]+0.1091214168497758*fMquad[23]-0.1091214168497758*fMquad[22]+0.1745942669596413*fMquad[20]-0.1745942669596413*fMquad[19]-0.1091214168497758*fMquad[17]+0.1091214168497758*fMquad[16]-0.1091214168497758*fMquad[14]+0.1091214168497758*fMquad[13]-0.1745942669596413*fMquad[11]+0.1745942669596413*fMquad[10]; 
  fMOut[6] = 0.1091214168497758*fMquad[26]-0.1091214168497758*(fMquad[25]+fMquad[23])+0.1091214168497758*(fMquad[22]+fMquad[17])-0.1091214168497758*(fMquad[16]+fMquad[14])+0.1091214168497758*fMquad[13]+0.1745942669596413*fMquad[8]-0.1745942669596413*(fMquad[7]+fMquad[5])+0.1745942669596413*fMquad[4]; 
  fMOut[7] = 0.0542228679727088*(fMquad[26]+fMquad[25])+0.0867565887563341*fMquad[24]+0.0542228679727088*(fMquad[23]+fMquad[22])+0.0867565887563341*(fMquad[21]+fMquad[20]+fMquad[19])+0.1388105420101346*fMquad[18]+0.0542228679727088*(fMquad[17]+fMquad[16])+0.0867565887563341*fMquad[15]+0.0542228679727088*(fMquad[14]+fMquad[13])+0.0867565887563341*(fMquad[12]+fMquad[11]+fMquad[10])+0.1388105420101346*fMquad[9]-0.1084457359454176*(fMquad[8]+fMquad[7])-0.1735131775126682*fMquad[6]-0.1084457359454176*(fMquad[5]+fMquad[4])-0.1735131775126682*(fMquad[3]+fMquad[2]+fMquad[1])-0.2776210840202691*fMquad[0]; 
  fMOut[8] = 0.0542228679727088*(fMquad[26]+fMquad[25])+0.0867565887563341*fMquad[24]+0.0542228679727088*(fMquad[23]+fMquad[22])+0.0867565887563341*fMquad[21]-0.1084457359454176*(fMquad[20]+fMquad[19])-0.1735131775126682*fMquad[18]+0.0542228679727088*(fMquad[17]+fMquad[16])+0.0867565887563341*fMquad[15]+0.0542228679727088*(fMquad[14]+fMquad[13])+0.0867565887563341*fMquad[12]-0.1084457359454176*(fMquad[11]+fMquad[10])-0.1735131775126682*fMquad[9]+0.0867565887563341*(fMquad[8]+fMquad[7])+0.1388105420101346*fMquad[6]+0.0867565887563341*(fMquad[5]+fMquad[4])+0.1388105420101346*fMquad[3]-0.1735131775126682*(fMquad[2]+fMquad[1])-0.2776210840202691*fMquad[0]; 
  fMOut[9] = 0.0542228679727088*(fMquad[26]+fMquad[25])-0.1084457359454176*fMquad[24]+0.0542228679727088*(fMquad[23]+fMquad[22])-0.1084457359454176*fMquad[21]+0.0867565887563341*(fMquad[20]+fMquad[19])-0.1735131775126682*fMquad[18]+0.0542228679727088*(fMquad[17]+fMquad[16])-0.1084457359454176*fMquad[15]+0.0542228679727088*(fMquad[14]+fMquad[13])-0.1084457359454176*fMquad[12]+0.0867565887563341*(fMquad[11]+fMquad[10])-0.1735131775126682*fMquad[9]+0.0867565887563341*(fMquad[8]+fMquad[7])-0.1735131775126682*fMquad[6]+0.0867565887563341*(fMquad[5]+fMquad[4])-0.1735131775126682*fMquad[3]+0.1388105420101346*(fMquad[2]+fMquad[1])-0.2776210840202691*fMquad[0]; 
  fMOut[10] = 0.1464017435263137*fMquad[26]-0.1464017435263137*(fMquad[25]+fMquad[23])+0.1464017435263137*fMquad[22]-0.1464017435263137*fMquad[17]+0.1464017435263137*(fMquad[16]+fMquad[14])-0.1464017435263137*fMquad[13]; 
  fMOut[11] = 0.07274761123318388*(fMquad[26]+fMquad[25])+0.1163961779730942*fMquad[24]-0.07274761123318388*(fMquad[23]+fMquad[22])-0.1163961779730942*fMquad[21]+0.07274761123318388*(fMquad[17]+fMquad[16])+0.1163961779730942*fMquad[15]-0.07274761123318388*(fMquad[14]+fMquad[13])-0.1163961779730942*fMquad[12]-0.1454952224663677*(fMquad[8]+fMquad[7])-0.2327923559461884*fMquad[6]+0.1454952224663677*(fMquad[5]+fMquad[4])+0.2327923559461884*fMquad[3]; 
  fMOut[12] = 0.07274761123318388*(fMquad[26]+fMquad[25])+0.1163961779730942*fMquad[24]+0.07274761123318388*(fMquad[23]+fMquad[22])+0.1163961779730942*fMquad[21]-0.1454952224663677*(fMquad[20]+fMquad[19])-0.2327923559461884*fMquad[18]-0.07274761123318388*(fMquad[17]+fMquad[16])-0.1163961779730942*fMquad[15]-0.07274761123318388*(fMquad[14]+fMquad[13])-0.1163961779730942*fMquad[12]+0.1454952224663677*(fMquad[11]+fMquad[10])+0.2327923559461884*fMquad[9]; 
  fMOut[13] = 0.07274761123318388*fMquad[26]-0.07274761123318388*fMquad[25]+0.07274761123318388*fMquad[23]-0.07274761123318388*fMquad[22]+0.1163961779730942*fMquad[20]-0.1163961779730942*fMquad[19]+0.07274761123318388*fMquad[17]-0.07274761123318388*fMquad[16]+0.07274761123318388*fMquad[14]-0.07274761123318388*fMquad[13]+0.1163961779730942*fMquad[11]-0.1163961779730942*fMquad[10]-0.1454952224663677*fMquad[8]+0.1454952224663677*fMquad[7]-0.1454952224663677*fMquad[5]+0.1454952224663677*fMquad[4]-0.2327923559461884*fMquad[2]+0.2327923559461884*fMquad[1]; 
  fMOut[14] = 0.07274761123318388*fMquad[26]-0.07274761123318388*fMquad[25]+0.07274761123318388*fMquad[23]-0.07274761123318388*fMquad[22]-0.1454952224663677*fMquad[20]+0.1454952224663677*fMquad[19]+0.07274761123318388*fMquad[17]-0.07274761123318388*fMquad[16]+0.07274761123318388*fMquad[14]-0.07274761123318388*fMquad[13]-0.1454952224663677*fMquad[11]+0.1454952224663677*fMquad[10]+0.1163961779730942*fMquad[8]-0.1163961779730942*fMquad[7]+0.1163961779730942*fMquad[5]-0.1163961779730942*fMquad[4]-0.2327923559461884*fMquad[2]+0.2327923559461884*fMquad[1]; 
  fMOut[15] = 0.07274761123318388*(fMquad[26]+fMquad[25])-0.1454952224663677*fMquad[24]+0.07274761123318388*(fMquad[23]+fMquad[22])-0.1454952224663677*fMquad[21]+0.1163961779730942*(fMquad[20]+fMquad[19])-0.2327923559461884*fMquad[18]-0.07274761123318388*(fMquad[17]+fMquad[16])+0.1454952224663677*fMquad[15]-0.07274761123318388*(fMquad[14]+fMquad[13])+0.1454952224663677*fMquad[12]-0.1163961779730942*(fMquad[11]+fMquad[10])+0.2327923559461884*fMquad[9]; 
  fMOut[16] = 0.07274761123318388*(fMquad[26]+fMquad[25])-0.1454952224663677*fMquad[24]-0.07274761123318388*(fMquad[23]+fMquad[22])+0.1454952224663677*fMquad[21]+0.07274761123318388*(fMquad[17]+fMquad[16])-0.1454952224663677*fMquad[15]-0.07274761123318388*(fMquad[14]+fMquad[13])+0.1454952224663677*fMquad[12]+0.1163961779730942*(fMquad[8]+fMquad[7])-0.2327923559461884*fMquad[6]-0.1163961779730942*(fMquad[5]+fMquad[4])+0.2327923559461884*fMquad[3]; 
  fMOut[17] = 0.09760116235087586*fMquad[26]-0.09760116235087586*(fMquad[25]+fMquad[23])+0.09760116235087586*(fMquad[22]+fMquad[17])-0.09760116235087586*(fMquad[16]+fMquad[14])+0.09760116235087586*fMquad[13]-0.1952023247017517*fMquad[8]+0.1952023247017517*(fMquad[7]+fMquad[5])-0.1952023247017517*fMquad[4]; 
  fMOut[18] = 0.09760116235087589*fMquad[26]-0.09760116235087589*fMquad[25]+0.09760116235087589*fMquad[23]-0.09760116235087589*fMquad[22]-0.1952023247017516*fMquad[20]+0.1952023247017516*fMquad[19]-0.09760116235087589*fMquad[17]+0.09760116235087589*fMquad[16]-0.09760116235087589*fMquad[14]+0.09760116235087589*fMquad[13]+0.1952023247017516*fMquad[11]-0.1952023247017516*fMquad[10]; 
  fMOut[19] = 0.09760116235087589*(fMquad[26]+fMquad[25])-0.1952023247017516*fMquad[24]-0.09760116235087589*(fMquad[23]+fMquad[22])+0.1952023247017516*fMquad[21]-0.09760116235087589*(fMquad[17]+fMquad[16])+0.1952023247017516*fMquad[15]+0.09760116235087589*(fMquad[14]+fMquad[13])-0.1952023247017516*fMquad[12]; 

}
