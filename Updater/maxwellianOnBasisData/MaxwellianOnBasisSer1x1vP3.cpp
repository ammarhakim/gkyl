#include <MaxwellianOnBasisModDecl.h>

void MaxwellianOnBasisGauss1x1vSer_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[4];
  m0Ord[0] = 0.7702725556588816*den[3]-0.5164305132317774*den[2]-0.416390039500913*den[1]+0.7071067811865475*den[0]; 
  m0Ord[1] = (-0.7702725556588816*den[3])-0.5164305132317774*den[2]+0.416390039500913*den[1]+0.7071067811865475*den[0]; 
  m0Ord[2] = 0.5701294036773671*den[3]+0.9681844646844028*den[2]+1.054672281193885*den[1]+0.7071067811865475*den[0]; 
  m0Ord[3] = (-0.5701294036773671*den[3])+0.9681844646844028*den[2]-1.054672281193885*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7702725556588816*flowU[3]-0.5164305132317774*flowU[2]-0.416390039500913*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[1] = (-0.7702725556588816*flowU[3])-0.5164305132317774*flowU[2]+0.416390039500913*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[2] = 0.5701294036773671*flowU[3]+0.9681844646844028*flowU[2]+1.054672281193885*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[3] = (-0.5701294036773671*flowU[3])+0.9681844646844028*flowU[2]-1.054672281193885*flowU[1]+0.7071067811865475*flowU[0]; 

  vtSqOrd[0] = 0.7702725556588816*vtSq[3]-0.5164305132317774*vtSq[2]-0.416390039500913*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[1] = (-0.7702725556588816*vtSq[3])-0.5164305132317774*vtSq[2]+0.416390039500913*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[2] = 0.5701294036773671*vtSq[3]+0.9681844646844028*vtSq[2]+1.054672281193885*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[3] = (-0.5701294036773671*vtSq[3])+0.9681844646844028*vtSq[2]-1.054672281193885*vtSq[1]+0.7071067811865475*vtSq[0]; 

  if ((vtSqOrd[0] <= 0.0) || (m0Ord[0] <= 0.0))
    fMFacOrd[0] = 0.;
  else
    fMFacOrd[0] = (0.3989422804014326*m0Ord[0])/sqrt(vtSqOrd[0]); 
  if ((vtSqOrd[1] <= 0.0) || (m0Ord[1] <= 0.0))
    fMFacOrd[1] = 0.;
  else
    fMFacOrd[1] = (0.3989422804014326*m0Ord[1])/sqrt(vtSqOrd[1]); 
  if ((vtSqOrd[2] <= 0.0) || (m0Ord[2] <= 0.0))
    fMFacOrd[2] = 0.;
  else
    fMFacOrd[2] = (0.3989422804014326*m0Ord[2])/sqrt(vtSqOrd[2]); 
  if ((vtSqOrd[3] <= 0.0) || (m0Ord[3] <= 0.0))
    fMFacOrd[3] = 0.;
  else
    fMFacOrd[3] = (0.3989422804014326*m0Ord[3])/sqrt(vtSqOrd[3]); 

}

void MaxwellianOnBasisGauss1x1vSer_P3_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[16];
  if ((vtSqOrd[0] <= 0.0) || (fMFacOrd[0] <= 0.0)) {
    fMquad[0] = 0;
    fMquad[1] = 0;
    fMquad[2] = 0;
    fMquad[3] = 0;
  } else {
    fMquad[0] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-0.1699905217924281*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
    fMquad[1] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]+0.1699905217924281*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
    fMquad[2] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]+0.4305681557970263*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
    fMquad[3] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-0.4305681557970263*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  };
  if ((vtSqOrd[1] <= 0.0) || (fMFacOrd[1] <= 0.0)) {
    fMquad[4] = 0;
    fMquad[5] = 0;
    fMquad[6] = 0;
    fMquad[7] = 0;
  } else {
    fMquad[4] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.1699905217924281*dxv[1],2.0))/vtSqOrd[1]); 
    fMquad[5] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.1699905217924281*dxv[1],2.0))/vtSqOrd[1]); 
    fMquad[6] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.4305681557970263*dxv[1],2.0))/vtSqOrd[1]); 
    fMquad[7] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.4305681557970263*dxv[1],2.0))/vtSqOrd[1]); 
  };
  if ((vtSqOrd[2] <= 0.0) || (fMFacOrd[2] <= 0.0)) {
    fMquad[8] = 0;
    fMquad[9] = 0;
    fMquad[10] = 0;
    fMquad[11] = 0;
  } else {
    fMquad[8] = fMFacOrd[2]*exp(-(0.5*std::pow((-1.0*flowUOrd[2])+wc[1]-0.1699905217924281*dxv[1],2.0))/vtSqOrd[2]); 
    fMquad[9] = fMFacOrd[2]*exp(-(0.5*std::pow((-1.0*flowUOrd[2])+wc[1]+0.1699905217924281*dxv[1],2.0))/vtSqOrd[2]); 
    fMquad[10] = fMFacOrd[2]*exp(-(0.5*std::pow((-1.0*flowUOrd[2])+wc[1]+0.4305681557970263*dxv[1],2.0))/vtSqOrd[2]); 
    fMquad[11] = fMFacOrd[2]*exp(-(0.5*std::pow((-1.0*flowUOrd[2])+wc[1]-0.4305681557970263*dxv[1],2.0))/vtSqOrd[2]); 
  };
  if ((vtSqOrd[3] <= 0.0) || (fMFacOrd[3] <= 0.0)) {
    fMquad[12] = 0;
    fMquad[13] = 0;
    fMquad[14] = 0;
    fMquad[15] = 0;
  } else {
    fMquad[12] = fMFacOrd[3]*exp(-(0.5*std::pow((-1.0*flowUOrd[3])+wc[1]-0.1699905217924281*dxv[1],2.0))/vtSqOrd[3]); 
    fMquad[13] = fMFacOrd[3]*exp(-(0.5*std::pow((-1.0*flowUOrd[3])+wc[1]+0.1699905217924281*dxv[1],2.0))/vtSqOrd[3]); 
    fMquad[14] = fMFacOrd[3]*exp(-(0.5*std::pow((-1.0*flowUOrd[3])+wc[1]+0.4305681557970263*dxv[1],2.0))/vtSqOrd[3]); 
    fMquad[15] = fMFacOrd[3]*exp(-(0.5*std::pow((-1.0*flowUOrd[3])+wc[1]-0.4305681557970263*dxv[1],2.0))/vtSqOrd[3]); 
  };

  fMOut[0] = 0.06050149664280098*(fMquad[15]+fMquad[14])+0.1134259259259259*(fMquad[13]+fMquad[12])+0.06050149664280098*(fMquad[11]+fMquad[10])+0.1134259259259259*(fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6])+0.2126466515053471*(fMquad[5]+fMquad[4])+0.1134259259259259*(fMquad[3]+fMquad[2])+0.2126466515053471*(fMquad[1]+fMquad[0]); 
  fMOut[1] = (-0.0902399088477601*(fMquad[15]+fMquad[14]))-0.169178380445011*(fMquad[13]+fMquad[12])+0.0902399088477601*(fMquad[11]+fMquad[10])+0.169178380445011*(fMquad[9]+fMquad[8])+0.06679249447653642*(fMquad[7]+fMquad[6])+0.1252200515903253*(fMquad[5]+fMquad[4])-0.06679249447653642*(fMquad[3]+fMquad[2])-0.1252200515903253*(fMquad[1]+fMquad[0]); 
  fMOut[2] = (-0.0902399088477601*fMquad[15])+0.0902399088477601*fMquad[14]+0.06679249447653642*fMquad[13]-0.06679249447653642*fMquad[12]-0.0902399088477601*fMquad[11]+0.0902399088477601*fMquad[10]+0.06679249447653642*fMquad[9]-0.06679249447653642*fMquad[8]-0.169178380445011*fMquad[7]+0.169178380445011*fMquad[6]+0.1252200515903253*fMquad[5]-0.1252200515903253*fMquad[4]-0.169178380445011*fMquad[3]+0.169178380445011*fMquad[2]+0.1252200515903253*fMquad[1]-0.1252200515903253*fMquad[0]; 
  fMOut[3] = 0.1345956976391758*fMquad[15]-0.1345956976391758*fMquad[14]-0.09962313244682941*fMquad[13]+0.09962313244682941*fMquad[12]-0.1345956976391758*fMquad[11]+0.1345956976391758*fMquad[10]+0.09962313244682941*fMquad[9]-0.09962313244682941*(fMquad[8]+fMquad[7])+0.09962313244682941*fMquad[6]+0.07373763569415742*fMquad[5]-0.07373763569415742*fMquad[4]+0.09962313244682941*fMquad[3]-0.09962313244682941*fMquad[2]-0.07373763569415742*fMquad[1]+0.07373763569415742*fMquad[0]; 
  fMOut[4] = 0.08283983508321342*(fMquad[15]+fMquad[14])+0.1553050010207066*(fMquad[13]+fMquad[12])+0.08283983508321342*(fMquad[11]+fMquad[10])+0.1553050010207066*(fMquad[9]+fMquad[8])-0.08283983508321344*(fMquad[7]+fMquad[6])-0.1553050010207066*(fMquad[5]+fMquad[4])-0.08283983508321344*(fMquad[3]+fMquad[2])-0.1553050010207066*(fMquad[1]+fMquad[0]); 
  fMOut[5] = 0.08283983508321342*(fMquad[15]+fMquad[14])-0.08283983508321344*(fMquad[13]+fMquad[12])+0.08283983508321342*(fMquad[11]+fMquad[10])-0.08283983508321344*(fMquad[9]+fMquad[8])+0.1553050010207066*(fMquad[7]+fMquad[6])-0.1553050010207066*(fMquad[5]+fMquad[4])+0.1553050010207066*(fMquad[3]+fMquad[2])-0.1553050010207066*(fMquad[1]+fMquad[0]); 
  fMOut[6] = (-0.1235582519719726*fMquad[15])+0.1235582519719726*fMquad[14]+0.0914535926259784*fMquad[13]-0.0914535926259784*fMquad[12]-0.1235582519719726*fMquad[11]+0.1235582519719726*fMquad[10]+0.0914535926259784*fMquad[9]-0.0914535926259784*fMquad[8]+0.1235582519719726*fMquad[7]-0.1235582519719726*fMquad[6]-0.09145359262597841*fMquad[5]+0.09145359262597841*fMquad[4]+0.1235582519719726*fMquad[3]-0.1235582519719726*fMquad[2]-0.09145359262597841*fMquad[1]+0.09145359262597841*fMquad[0]; 
  fMOut[7] = (-0.1235582519719726*(fMquad[15]+fMquad[14]))+0.1235582519719726*(fMquad[13]+fMquad[12])+0.1235582519719726*(fMquad[11]+fMquad[10])-0.1235582519719726*(fMquad[9]+fMquad[8])+0.09145359262597838*(fMquad[7]+fMquad[6])-0.09145359262597841*(fMquad[5]+fMquad[4])-0.09145359262597838*(fMquad[3]+fMquad[2])+0.09145359262597841*(fMquad[1]+fMquad[0]); 
  fMOut[8] = (-0.04878143318703135*(fMquad[15]+fMquad[14]))-0.0914535926259784*(fMquad[13]+fMquad[12])+0.04878143318703135*(fMquad[11]+fMquad[10])+0.0914535926259784*(fMquad[9]+fMquad[8])-0.1235582519719726*(fMquad[7]+fMquad[6])-0.2316423545429344*(fMquad[5]+fMquad[4])+0.1235582519719726*(fMquad[3]+fMquad[2])+0.2316423545429344*(fMquad[1]+fMquad[0]); 
  fMOut[9] = (-0.04878143318703135*fMquad[15])+0.04878143318703135*fMquad[14]-0.1235582519719726*fMquad[13]+0.1235582519719726*fMquad[12]-0.04878143318703135*fMquad[11]+0.04878143318703135*fMquad[10]-0.1235582519719726*fMquad[9]+0.1235582519719726*fMquad[8]-0.0914535926259784*fMquad[7]+0.0914535926259784*fMquad[6]-0.2316423545429344*fMquad[5]+0.2316423545429344*fMquad[4]-0.0914535926259784*fMquad[3]+0.0914535926259784*fMquad[2]-0.2316423545429344*fMquad[1]+0.2316423545429344*fMquad[0]; 
  fMOut[10] = 0.07275906099067718*fMquad[15]-0.07275906099067718*fMquad[14]-0.05385376870821617*fMquad[13]+0.05385376870821617*fMquad[12]-0.07275906099067718*fMquad[11]+0.07275906099067718*fMquad[10]+0.05385376870821617*fMquad[9]-0.05385376870821617*fMquad[8]+0.1842910673957039*fMquad[7]-0.1842910673957039*fMquad[6]-0.1364059456428416*fMquad[5]+0.1364059456428416*fMquad[4]-0.1842910673957039*fMquad[3]+0.1842910673957039*fMquad[2]+0.1364059456428416*fMquad[1]-0.1364059456428416*fMquad[0]; 
  fMOut[11] = 0.07275906099067721*fMquad[15]-0.07275906099067721*fMquad[14]+0.1842910673957038*fMquad[13]-0.1842910673957038*fMquad[12]-0.07275906099067721*fMquad[11]+0.07275906099067721*fMquad[10]-0.1842910673957038*fMquad[9]+0.1842910673957038*fMquad[8]-0.05385376870821612*fMquad[7]+0.05385376870821612*fMquad[6]-0.1364059456428416*fMquad[5]+0.1364059456428416*fMquad[4]+0.05385376870821612*fMquad[3]-0.05385376870821612*fMquad[2]+0.1364059456428416*fMquad[1]-0.1364059456428416*fMquad[0]; 

}
void GkMaxwellianOnBasisGauss1x1vSer_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[4];
  m0Ord[0] = 0.7702725556588816*den[3]-0.5164305132317774*den[2]-0.416390039500913*den[1]+0.7071067811865475*den[0]; 
  m0Ord[1] = (-0.7702725556588816*den[3])-0.5164305132317774*den[2]+0.416390039500913*den[1]+0.7071067811865475*den[0]; 
  m0Ord[2] = 0.5701294036773671*den[3]+0.9681844646844028*den[2]+1.054672281193885*den[1]+0.7071067811865475*den[0]; 
  m0Ord[3] = (-0.5701294036773671*den[3])+0.9681844646844028*den[2]-1.054672281193885*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7702725556588816*flowU[3]-0.5164305132317774*flowU[2]-0.416390039500913*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[1] = (-0.7702725556588816*flowU[3])-0.5164305132317774*flowU[2]+0.416390039500913*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[2] = 0.5701294036773671*flowU[3]+0.9681844646844028*flowU[2]+1.054672281193885*flowU[1]+0.7071067811865475*flowU[0]; 
  flowUOrd[3] = (-0.5701294036773671*flowU[3])+0.9681844646844028*flowU[2]-1.054672281193885*flowU[1]+0.7071067811865475*flowU[0]; 

  vtSqOrd[0] = 0.7702725556588816*vtSq[3]-0.5164305132317774*vtSq[2]-0.416390039500913*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[1] = (-0.7702725556588816*vtSq[3])-0.5164305132317774*vtSq[2]+0.416390039500913*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[2] = 0.5701294036773671*vtSq[3]+0.9681844646844028*vtSq[2]+1.054672281193885*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[3] = (-0.5701294036773671*vtSq[3])+0.9681844646844028*vtSq[2]-1.054672281193885*vtSq[1]+0.7071067811865475*vtSq[0]; 

  if ((vtSqOrd[0] <= 0.0) || (m0Ord[0] <= 0.0))
    fMFacOrd[0] = 0.;
  else
    fMFacOrd[0] = (0.3989422804014326*m0Ord[0]*(0.7702725556588816*bmag[3]-0.5164305132317774*bmag[2]-0.416390039500913*bmag[1]+0.7071067811865475*bmag[0]))/sqrt(vtSqOrd[0]); 
  if ((vtSqOrd[1] <= 0.0) || (m0Ord[1] <= 0.0))
    fMFacOrd[1] = 0.;
  else
    fMFacOrd[1] = (0.3989422804014326*m0Ord[1]*((-0.7702725556588816*bmag[3])-0.5164305132317774*bmag[2]+0.416390039500913*bmag[1]+0.7071067811865475*bmag[0]))/sqrt(vtSqOrd[1]); 
  if ((vtSqOrd[2] <= 0.0) || (m0Ord[2] <= 0.0))
    fMFacOrd[2] = 0.;
  else
    fMFacOrd[2] = (0.3989422804014326*m0Ord[2]*(0.5701294036773671*bmag[3]+0.9681844646844028*bmag[2]+1.054672281193885*bmag[1]+0.7071067811865475*bmag[0]))/sqrt(vtSqOrd[2]); 
  if ((vtSqOrd[3] <= 0.0) || (m0Ord[3] <= 0.0))
    fMFacOrd[3] = 0.;
  else
    fMFacOrd[3] = (0.3989422804014326*((-0.5701294036773671*bmag[3])+0.9681844646844028*bmag[2]-1.054672281193885*bmag[1]+0.7071067811865475*bmag[0])*m0Ord[3])/sqrt(vtSqOrd[3]); 

}

void GkMaxwellianOnBasisGauss1x1vSerUz_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd) {

  double m0Ord[4];
  m0Ord[0] = 0.7702725556588816*den[3]-0.5164305132317774*den[2]-0.416390039500913*den[1]+0.7071067811865475*den[0]; 
  m0Ord[1] = (-0.7702725556588816*den[3])-0.5164305132317774*den[2]+0.416390039500913*den[1]+0.7071067811865475*den[0]; 
  m0Ord[2] = 0.5701294036773671*den[3]+0.9681844646844028*den[2]+1.054672281193885*den[1]+0.7071067811865475*den[0]; 
  m0Ord[3] = (-0.5701294036773671*den[3])+0.9681844646844028*den[2]-1.054672281193885*den[1]+0.7071067811865475*den[0]; 

  flowUOrd[0] = 0.7702725556588816*flowU[11]-0.5164305132317774*flowU[10]-0.416390039500913*flowU[9]+0.7071067811865475*flowU[8]; 
  flowUOrd[1] = (-0.7702725556588816*flowU[11])-0.5164305132317774*flowU[10]+0.416390039500913*flowU[9]+0.7071067811865475*flowU[8]; 
  flowUOrd[2] = 0.5701294036773671*flowU[11]+0.9681844646844028*flowU[10]+1.054672281193885*flowU[9]+0.7071067811865475*flowU[8]; 
  flowUOrd[3] = (-0.5701294036773671*flowU[11])+0.9681844646844028*flowU[10]-1.054672281193885*flowU[9]+0.7071067811865475*flowU[8]; 

  vtSqOrd[0] = 0.7702725556588816*vtSq[3]-0.5164305132317774*vtSq[2]-0.416390039500913*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[1] = (-0.7702725556588816*vtSq[3])-0.5164305132317774*vtSq[2]+0.416390039500913*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[2] = 0.5701294036773671*vtSq[3]+0.9681844646844028*vtSq[2]+1.054672281193885*vtSq[1]+0.7071067811865475*vtSq[0]; 
  vtSqOrd[3] = (-0.5701294036773671*vtSq[3])+0.9681844646844028*vtSq[2]-1.054672281193885*vtSq[1]+0.7071067811865475*vtSq[0]; 

  if ((vtSqOrd[0] <= 0.0) || (m0Ord[0] <= 0.0))
    fMFacOrd[0] = 0.;
  else
    fMFacOrd[0] = (0.3989422804014326*m0Ord[0]*(0.7702725556588816*bmag[3]-0.5164305132317774*bmag[2]-0.416390039500913*bmag[1]+0.7071067811865475*bmag[0]))/sqrt(vtSqOrd[0]); 
  if ((vtSqOrd[1] <= 0.0) || (m0Ord[1] <= 0.0))
    fMFacOrd[1] = 0.;
  else
    fMFacOrd[1] = (0.3989422804014326*m0Ord[1]*((-0.7702725556588816*bmag[3])-0.5164305132317774*bmag[2]+0.416390039500913*bmag[1]+0.7071067811865475*bmag[0]))/sqrt(vtSqOrd[1]); 
  if ((vtSqOrd[2] <= 0.0) || (m0Ord[2] <= 0.0))
    fMFacOrd[2] = 0.;
  else
    fMFacOrd[2] = (0.3989422804014326*m0Ord[2]*(0.5701294036773671*bmag[3]+0.9681844646844028*bmag[2]+1.054672281193885*bmag[1]+0.7071067811865475*bmag[0]))/sqrt(vtSqOrd[2]); 
  if ((vtSqOrd[3] <= 0.0) || (m0Ord[3] <= 0.0))
    fMFacOrd[3] = 0.;
  else
    fMFacOrd[3] = (0.3989422804014326*((-0.5701294036773671*bmag[3])+0.9681844646844028*bmag[2]-1.054672281193885*bmag[1]+0.7071067811865475*bmag[0])*m0Ord[3])/sqrt(vtSqOrd[3]); 

}

void GkMaxwellianOnBasisGauss1x1vSer_P3_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut) {

  double fMquad[16];
  if ((vtSqOrd[0] <= 0.0) || (fMFacOrd[0] <= 0.0)) {
    fMquad[0] = 9.999999999999999e-41;
    fMquad[1] = 9.999999999999999e-41;
    fMquad[2] = 9.999999999999999e-41;
    fMquad[3] = 9.999999999999999e-41;
  } else {
    fMquad[0] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-0.1699905217924281*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
    fMquad[1] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]+0.1699905217924281*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
    fMquad[2] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]+0.4305681557970263*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
    fMquad[3] = fMFacOrd[0]*exp(-(0.5*std::pow(wc[1]-0.4305681557970263*dxv[1]-1.0*flowUOrd[0],2.0))/vtSqOrd[0]); 
  };
  if ((vtSqOrd[1] <= 0.0) || (fMFacOrd[1] <= 0.0)) {
    fMquad[4] = 9.999999999999999e-41;
    fMquad[5] = 9.999999999999999e-41;
    fMquad[6] = 9.999999999999999e-41;
    fMquad[7] = 9.999999999999999e-41;
  } else {
    fMquad[4] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.1699905217924281*dxv[1],2.0))/vtSqOrd[1]); 
    fMquad[5] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.1699905217924281*dxv[1],2.0))/vtSqOrd[1]); 
    fMquad[6] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]+0.4305681557970263*dxv[1],2.0))/vtSqOrd[1]); 
    fMquad[7] = fMFacOrd[1]*exp(-(0.5*std::pow(wc[1]-1.0*flowUOrd[1]-0.4305681557970263*dxv[1],2.0))/vtSqOrd[1]); 
  };
  if ((vtSqOrd[2] <= 0.0) || (fMFacOrd[2] <= 0.0)) {
    fMquad[8] = 9.999999999999999e-41;
    fMquad[9] = 9.999999999999999e-41;
    fMquad[10] = 9.999999999999999e-41;
    fMquad[11] = 9.999999999999999e-41;
  } else {
    fMquad[8] = fMFacOrd[2]*exp(-(0.5*std::pow((-1.0*flowUOrd[2])+wc[1]-0.1699905217924281*dxv[1],2.0))/vtSqOrd[2]); 
    fMquad[9] = fMFacOrd[2]*exp(-(0.5*std::pow((-1.0*flowUOrd[2])+wc[1]+0.1699905217924281*dxv[1],2.0))/vtSqOrd[2]); 
    fMquad[10] = fMFacOrd[2]*exp(-(0.5*std::pow((-1.0*flowUOrd[2])+wc[1]+0.4305681557970263*dxv[1],2.0))/vtSqOrd[2]); 
    fMquad[11] = fMFacOrd[2]*exp(-(0.5*std::pow((-1.0*flowUOrd[2])+wc[1]-0.4305681557970263*dxv[1],2.0))/vtSqOrd[2]); 
  };
  if ((vtSqOrd[3] <= 0.0) || (fMFacOrd[3] <= 0.0)) {
    fMquad[12] = 9.999999999999999e-41;
    fMquad[13] = 9.999999999999999e-41;
    fMquad[14] = 9.999999999999999e-41;
    fMquad[15] = 9.999999999999999e-41;
  } else {
    fMquad[12] = fMFacOrd[3]*exp(-(0.5*std::pow((-1.0*flowUOrd[3])+wc[1]-0.1699905217924281*dxv[1],2.0))/vtSqOrd[3]); 
    fMquad[13] = fMFacOrd[3]*exp(-(0.5*std::pow((-1.0*flowUOrd[3])+wc[1]+0.1699905217924281*dxv[1],2.0))/vtSqOrd[3]); 
    fMquad[14] = fMFacOrd[3]*exp(-(0.5*std::pow((-1.0*flowUOrd[3])+wc[1]+0.4305681557970263*dxv[1],2.0))/vtSqOrd[3]); 
    fMquad[15] = fMFacOrd[3]*exp(-(0.5*std::pow((-1.0*flowUOrd[3])+wc[1]-0.4305681557970263*dxv[1],2.0))/vtSqOrd[3]); 
  };

  fMOut[0] = 0.06050149664280098*(fMquad[15]+fMquad[14])+0.1134259259259259*(fMquad[13]+fMquad[12])+0.06050149664280098*(fMquad[11]+fMquad[10])+0.1134259259259259*(fMquad[9]+fMquad[8]+fMquad[7]+fMquad[6])+0.2126466515053471*(fMquad[5]+fMquad[4])+0.1134259259259259*(fMquad[3]+fMquad[2])+0.2126466515053471*(fMquad[1]+fMquad[0]); 
  fMOut[1] = (-0.0902399088477601*(fMquad[15]+fMquad[14]))-0.169178380445011*(fMquad[13]+fMquad[12])+0.0902399088477601*(fMquad[11]+fMquad[10])+0.169178380445011*(fMquad[9]+fMquad[8])+0.06679249447653642*(fMquad[7]+fMquad[6])+0.1252200515903253*(fMquad[5]+fMquad[4])-0.06679249447653642*(fMquad[3]+fMquad[2])-0.1252200515903253*(fMquad[1]+fMquad[0]); 
  fMOut[2] = (-0.0902399088477601*fMquad[15])+0.0902399088477601*fMquad[14]+0.06679249447653642*fMquad[13]-0.06679249447653642*fMquad[12]-0.0902399088477601*fMquad[11]+0.0902399088477601*fMquad[10]+0.06679249447653642*fMquad[9]-0.06679249447653642*fMquad[8]-0.169178380445011*fMquad[7]+0.169178380445011*fMquad[6]+0.1252200515903253*fMquad[5]-0.1252200515903253*fMquad[4]-0.169178380445011*fMquad[3]+0.169178380445011*fMquad[2]+0.1252200515903253*fMquad[1]-0.1252200515903253*fMquad[0]; 
  fMOut[3] = 0.1345956976391758*fMquad[15]-0.1345956976391758*fMquad[14]-0.09962313244682941*fMquad[13]+0.09962313244682941*fMquad[12]-0.1345956976391758*fMquad[11]+0.1345956976391758*fMquad[10]+0.09962313244682941*fMquad[9]-0.09962313244682941*(fMquad[8]+fMquad[7])+0.09962313244682941*fMquad[6]+0.07373763569415742*fMquad[5]-0.07373763569415742*fMquad[4]+0.09962313244682941*fMquad[3]-0.09962313244682941*fMquad[2]-0.07373763569415742*fMquad[1]+0.07373763569415742*fMquad[0]; 
  fMOut[4] = 0.08283983508321342*(fMquad[15]+fMquad[14])+0.1553050010207066*(fMquad[13]+fMquad[12])+0.08283983508321342*(fMquad[11]+fMquad[10])+0.1553050010207066*(fMquad[9]+fMquad[8])-0.08283983508321344*(fMquad[7]+fMquad[6])-0.1553050010207066*(fMquad[5]+fMquad[4])-0.08283983508321344*(fMquad[3]+fMquad[2])-0.1553050010207066*(fMquad[1]+fMquad[0]); 
  fMOut[5] = 0.08283983508321342*(fMquad[15]+fMquad[14])-0.08283983508321344*(fMquad[13]+fMquad[12])+0.08283983508321342*(fMquad[11]+fMquad[10])-0.08283983508321344*(fMquad[9]+fMquad[8])+0.1553050010207066*(fMquad[7]+fMquad[6])-0.1553050010207066*(fMquad[5]+fMquad[4])+0.1553050010207066*(fMquad[3]+fMquad[2])-0.1553050010207066*(fMquad[1]+fMquad[0]); 
  fMOut[6] = (-0.1235582519719726*fMquad[15])+0.1235582519719726*fMquad[14]+0.0914535926259784*fMquad[13]-0.0914535926259784*fMquad[12]-0.1235582519719726*fMquad[11]+0.1235582519719726*fMquad[10]+0.0914535926259784*fMquad[9]-0.0914535926259784*fMquad[8]+0.1235582519719726*fMquad[7]-0.1235582519719726*fMquad[6]-0.09145359262597841*fMquad[5]+0.09145359262597841*fMquad[4]+0.1235582519719726*fMquad[3]-0.1235582519719726*fMquad[2]-0.09145359262597841*fMquad[1]+0.09145359262597841*fMquad[0]; 
  fMOut[7] = (-0.1235582519719726*(fMquad[15]+fMquad[14]))+0.1235582519719726*(fMquad[13]+fMquad[12])+0.1235582519719726*(fMquad[11]+fMquad[10])-0.1235582519719726*(fMquad[9]+fMquad[8])+0.09145359262597838*(fMquad[7]+fMquad[6])-0.09145359262597841*(fMquad[5]+fMquad[4])-0.09145359262597838*(fMquad[3]+fMquad[2])+0.09145359262597841*(fMquad[1]+fMquad[0]); 
  fMOut[8] = (-0.04878143318703135*(fMquad[15]+fMquad[14]))-0.0914535926259784*(fMquad[13]+fMquad[12])+0.04878143318703135*(fMquad[11]+fMquad[10])+0.0914535926259784*(fMquad[9]+fMquad[8])-0.1235582519719726*(fMquad[7]+fMquad[6])-0.2316423545429344*(fMquad[5]+fMquad[4])+0.1235582519719726*(fMquad[3]+fMquad[2])+0.2316423545429344*(fMquad[1]+fMquad[0]); 
  fMOut[9] = (-0.04878143318703135*fMquad[15])+0.04878143318703135*fMquad[14]-0.1235582519719726*fMquad[13]+0.1235582519719726*fMquad[12]-0.04878143318703135*fMquad[11]+0.04878143318703135*fMquad[10]-0.1235582519719726*fMquad[9]+0.1235582519719726*fMquad[8]-0.0914535926259784*fMquad[7]+0.0914535926259784*fMquad[6]-0.2316423545429344*fMquad[5]+0.2316423545429344*fMquad[4]-0.0914535926259784*fMquad[3]+0.0914535926259784*fMquad[2]-0.2316423545429344*fMquad[1]+0.2316423545429344*fMquad[0]; 
  fMOut[10] = 0.07275906099067718*fMquad[15]-0.07275906099067718*fMquad[14]-0.05385376870821617*fMquad[13]+0.05385376870821617*fMquad[12]-0.07275906099067718*fMquad[11]+0.07275906099067718*fMquad[10]+0.05385376870821617*fMquad[9]-0.05385376870821617*fMquad[8]+0.1842910673957039*fMquad[7]-0.1842910673957039*fMquad[6]-0.1364059456428416*fMquad[5]+0.1364059456428416*fMquad[4]-0.1842910673957039*fMquad[3]+0.1842910673957039*fMquad[2]+0.1364059456428416*fMquad[1]-0.1364059456428416*fMquad[0]; 
  fMOut[11] = 0.07275906099067721*fMquad[15]-0.07275906099067721*fMquad[14]+0.1842910673957038*fMquad[13]-0.1842910673957038*fMquad[12]-0.07275906099067721*fMquad[11]+0.07275906099067721*fMquad[10]-0.1842910673957038*fMquad[9]+0.1842910673957038*fMquad[8]-0.05385376870821612*fMquad[7]+0.05385376870821612*fMquad[6]-0.1364059456428416*fMquad[5]+0.1364059456428416*fMquad[4]+0.05385376870821612*fMquad[3]-0.05385376870821612*fMquad[2]+0.1364059456428416*fMquad[1]-0.1364059456428416*fMquad[0]; 

}
