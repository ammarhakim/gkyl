#include <DeltaFGyrokineticModDecl.h>
double DeltaFGyrokineticGenGeoSurf1x1vSer_x_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 
  hamil0R[5] = (0.2981423969999719*m_)/rdvpar2SqR; 

  double hamil1R[8]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 
  hamil1R[4] = 1.414213562373095*phi[2]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = cmag[0]*jacobTotInv[0]; 

  double alpha0R[3]; 
  alpha0R[0] = (0.6123724356957944*BstarZdBmagR[0]*hamil0R[2]*rdvpar2R)/m_; 
  alpha0R[1] = (1.369306393762915*BstarZdBmagR[0]*hamil0R[5]*rdvpar2R)/m_; 

  double alpha1R[3]; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.2165063509461096*BstarZdBmagR[0]*hamil0R[2]*rdvpar2R)/m_; 

  double incr[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double fUpOrd[3];
  if (0.7071067811865475*alpha0R[0]-0.9486832980505137*alpha0R[1] > 0) { 
    fUpOrd[0] = 0.7745966692414833*f1L[7]-1.5*f1L[6]+0.4472135954999579*f1L[5]+1.118033988749895*f1L[4]-1.161895003862225*f1L[3]-0.6708203932499369*f1L[2]+0.8660254037844386*f1L[1]+0.5*f1L[0]; 
  } else { 
    fUpOrd[0] = (-0.7745966692414833*f1R[7])-1.5*f1R[6]+0.4472135954999579*f1R[5]+1.118033988749895*f1R[4]+1.161895003862225*f1R[3]-0.6708203932499369*f1R[2]-0.8660254037844386*f1R[1]+0.5*f1R[0]; 
  } 
  if (0.7071067811865475*alpha0R[0] > 0) { 
    fUpOrd[1] = (-0.9682458365518543*f1L[7])-0.5590169943749475*f1L[5]+1.118033988749895*f1L[4]+0.8660254037844386*f1L[1]+0.5*f1L[0]; 
  } else { 
    fUpOrd[1] = 0.9682458365518543*f1R[7]-0.5590169943749475*f1R[5]+1.118033988749895*f1R[4]-0.8660254037844386*f1R[1]+0.5*f1R[0]; 
  } 
  if (0.9486832980505137*alpha0R[1]+0.7071067811865475*alpha0R[0] > 0) { 
    fUpOrd[2] = 0.7745966692414833*f1L[7]+1.5*f1L[6]+0.4472135954999579*f1L[5]+1.118033988749895*f1L[4]+1.161895003862225*f1L[3]+0.6708203932499369*f1L[2]+0.8660254037844386*f1L[1]+0.5*f1L[0]; 
  } else { 
    fUpOrd[2] = (-0.7745966692414833*f1R[7])+1.5*f1R[6]+0.4472135954999579*f1R[5]+1.118033988749895*f1R[4]-1.161895003862225*f1R[3]+0.6708203932499369*f1R[2]-0.8660254037844386*f1R[1]+0.5*f1R[0]; 
  } 

  double fUp[3];
  fUp[0] = 0.3928371006591929*fUpOrd[2]+0.6285393610547088*fUpOrd[1]+0.3928371006591929*fUpOrd[0]; 
  fUp[1] = 0.5270462766947305*fUpOrd[2]-0.5270462766947305*fUpOrd[0]; 
  fUp[2] = 0.3513641844631533*fUpOrd[2]-0.7027283689263067*fUpOrd[1]+0.3513641844631533*fUpOrd[0]; 

  incr[0] = 0.5*(alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[1] = -0.8660254037844386*(alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[2] = 0.1*(4.47213595499958*alpha0R[1]*fUp[2]+5.0*(alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1])); 
  incr[3] = -0.1*(7.745966692414834*alpha0R[1]*fUp[2]+8.660254037844386*(alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1])); 
  incr[4] = 1.118033988749895*(alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[5] = 0.1*(5.0*alpha0R[0]*fUp[2]+4.47213595499958*alpha0R[1]*fUp[1]); 
  incr[6] = 0.03333333333333333*(30.0*alpha0R[1]*fUp[2]+33.54101966249684*(alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1])); 
  incr[7] = -0.1*(8.660254037844387*alpha0R[0]*fUp[2]+7.745966692414834*alpha0R[1]*fUp[1]); 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += incr[3]*rdx2L; 
  outL[4] += -1.0*incr[4]*rdx2L; 
  outL[5] += -1.0*incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += incr[7]*rdx2L; 

  // alpha1 == 0, so nothing to do for alpha1*f0 and alpha1*f1 terms.
  return std::abs(alphaSurfAvgR); 
} 
double DeltaFGyrokineticGenGeoSurf1x1vSer_vpar_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 
  hamil0R[5] = (0.2981423969999719*m_)/rdvpar2SqR; 

  double hamil1R[8]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 
  hamil1R[4] = 1.414213562373095*phi[2]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = cmag[0]*jacobTotInv[0]; 

  double alpha0R[3]; 

  double alpha1R[3]; 
  alpha1R[0] = -(0.6123724356957944*BstarZdBmagR[0]*hamil1R[1]*rdx2R)/m_; 
  alpha1R[1] = -(1.369306393762915*BstarZdBmagR[0]*hamil1R[4]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.2165063509461096*BstarZdBmagR[0]*hamil1R[1]*rdx2R)/m_; 

  double incr[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double fUpOrd[3];
  if (0.7071067811865475*alpha1R[0]-0.9486832980505137*alpha1R[1] > 0) { 
    fUpOrd[0] = (-1.5*f1L[7])+0.7745966692414833*f1L[6]+1.118033988749895*f1L[5]+0.4472135954999579*f1L[4]-1.161895003862225*f1L[3]+0.8660254037844386*f1L[2]-0.6708203932499369*f1L[1]+0.5*f1L[0]; 
  } else { 
    fUpOrd[0] = (-1.5*f1R[7])-0.7745966692414833*f1R[6]+1.118033988749895*f1R[5]+0.4472135954999579*f1R[4]+1.161895003862225*f1R[3]-0.8660254037844386*f1R[2]-0.6708203932499369*f1R[1]+0.5*f1R[0]; 
  } 
  if (0.7071067811865475*alpha1R[0] > 0) { 
    fUpOrd[1] = (-0.9682458365518543*f1L[6])+1.118033988749895*f1L[5]-0.5590169943749475*f1L[4]+0.8660254037844386*f1L[2]+0.5*f1L[0]; 
  } else { 
    fUpOrd[1] = 0.9682458365518543*f1R[6]+1.118033988749895*f1R[5]-0.5590169943749475*f1R[4]-0.8660254037844386*f1R[2]+0.5*f1R[0]; 
  } 
  if (0.9486832980505137*alpha1R[1]+0.7071067811865475*alpha1R[0] > 0) { 
    fUpOrd[2] = 1.5*f1L[7]+0.7745966692414833*f1L[6]+1.118033988749895*f1L[5]+0.4472135954999579*f1L[4]+1.161895003862225*f1L[3]+0.8660254037844386*f1L[2]+0.6708203932499369*f1L[1]+0.5*f1L[0]; 
  } else { 
    fUpOrd[2] = 1.5*f1R[7]-0.7745966692414833*f1R[6]+1.118033988749895*f1R[5]+0.4472135954999579*f1R[4]-1.161895003862225*f1R[3]-0.8660254037844386*f1R[2]+0.6708203932499369*f1R[1]+0.5*f1R[0]; 
  } 

  double fUp[3];
  fUp[0] = 0.3928371006591929*fUpOrd[2]+0.6285393610547088*fUpOrd[1]+0.3928371006591929*fUpOrd[0]; 
  fUp[1] = 0.5270462766947305*fUpOrd[2]-0.5270462766947305*fUpOrd[0]; 
  fUp[2] = 0.3513641844631533*fUpOrd[2]-0.7027283689263067*fUpOrd[1]+0.3513641844631533*fUpOrd[0]; 

  incr[0] = 0.5*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.1*(4.47213595499958*alpha1R[1]*fUp[2]+5.0*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 
  incr[2] = -0.8660254037844386*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = -0.1*(7.745966692414834*alpha1R[1]*fUp[2]+8.660254037844386*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 
  incr[4] = 0.1*(5.0*alpha1R[0]*fUp[2]+4.47213595499958*alpha1R[1]*fUp[1]); 
  incr[5] = 1.118033988749895*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[6] = -0.1*(8.660254037844387*alpha1R[0]*fUp[2]+7.745966692414834*alpha1R[1]*fUp[1]); 
  incr[7] = 0.03333333333333333*(30.0*alpha1R[1]*fUp[2]+33.54101966249684*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 
  outL[4] += -1.0*incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += -1.0*incr[7]*rdvpar2L; 

  // linear term alpha1*f0 
  if (0.7071067811865475*alpha1R[0]-0.9486832980505137*alpha1R[1] > 0) { 
    fUpOrd[0] = (-1.5*f0L[7])+0.7745966692414833*f0L[6]+1.118033988749895*f0L[5]+0.4472135954999579*f0L[4]-1.161895003862225*f0L[3]+0.8660254037844386*f0L[2]-0.6708203932499369*f0L[1]+0.5*f0L[0]; 
  } else { 
    fUpOrd[0] = (-1.5*f0R[7])-0.7745966692414833*f0R[6]+1.118033988749895*f0R[5]+0.4472135954999579*f0R[4]+1.161895003862225*f0R[3]-0.8660254037844386*f0R[2]-0.6708203932499369*f0R[1]+0.5*f0R[0]; 
  } 
  if (0.7071067811865475*alpha1R[0] > 0) { 
    fUpOrd[1] = (-0.9682458365518543*f0L[6])+1.118033988749895*f0L[5]-0.5590169943749475*f0L[4]+0.8660254037844386*f0L[2]+0.5*f0L[0]; 
  } else { 
    fUpOrd[1] = 0.9682458365518543*f0R[6]+1.118033988749895*f0R[5]-0.5590169943749475*f0R[4]-0.8660254037844386*f0R[2]+0.5*f0R[0]; 
  } 
  if (0.9486832980505137*alpha1R[1]+0.7071067811865475*alpha1R[0] > 0) { 
    fUpOrd[2] = 1.5*f0L[7]+0.7745966692414833*f0L[6]+1.118033988749895*f0L[5]+0.4472135954999579*f0L[4]+1.161895003862225*f0L[3]+0.8660254037844386*f0L[2]+0.6708203932499369*f0L[1]+0.5*f0L[0]; 
  } else { 
    fUpOrd[2] = 1.5*f0R[7]-0.7745966692414833*f0R[6]+1.118033988749895*f0R[5]+0.4472135954999579*f0R[4]-1.161895003862225*f0R[3]-0.8660254037844386*f0R[2]+0.6708203932499369*f0R[1]+0.5*f0R[0]; 
  } 

  fUp[0] = 0.3928371006591929*fUpOrd[2]+0.6285393610547088*fUpOrd[1]+0.3928371006591929*fUpOrd[0]; 
  fUp[1] = 0.5270462766947305*fUpOrd[2]-0.5270462766947305*fUpOrd[0]; 
  fUp[2] = 0.3513641844631533*fUpOrd[2]-0.7027283689263067*fUpOrd[1]+0.3513641844631533*fUpOrd[0]; 

  incr[0] = 0.5*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.1*(4.47213595499958*alpha1R[1]*fUp[2]+5.0*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 
  incr[2] = -0.8660254037844386*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = -0.1*(7.745966692414834*alpha1R[1]*fUp[2]+8.660254037844386*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 
  incr[4] = 0.1*(5.0*alpha1R[0]*fUp[2]+4.47213595499958*alpha1R[1]*fUp[1]); 
  incr[5] = 1.118033988749895*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[6] = -0.1*(8.660254037844387*alpha1R[0]*fUp[2]+7.745966692414834*alpha1R[1]*fUp[1]); 
  incr[7] = 0.03333333333333333*(30.0*alpha1R[1]*fUp[2]+33.54101966249684*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 
  outL[4] += -1.0*incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += -1.0*incr[7]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double DeltaFGyrokineticGenGeoSurf1x1vSer_x_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 
  hamil0R[5] = (0.2981423969999719*m_)/rdvpar2SqR; 

  double hamil1R[8]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 
  hamil1R[4] = 1.414213562373095*phi[2]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R*wvparR+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmagR[1] = (0.2*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2R*wvparR+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmagR[2] = ((2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = ((b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (0.02857142857142857*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R*wvparR+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmagR[6] = (1.0*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R)/(q_*rdvpar2R); 

  double alpha0R[3]; 
  alpha0R[0] = (0.3535533905932737*(hamil0R[5]*(8.660254037844387*BstarZdBmagR[6]-6.708203932499369*BstarZdBmagR[3]+3.872983346207417*BstarZdBmagR[2])+hamil0R[2]*(3.872983346207417*BstarZdBmagR[4]-3.0*BstarZdBmagR[1]+1.732050807568877*BstarZdBmagR[0]))*rdvpar2R)/m_; 
  alpha0R[1] = (0.3535533905932737*(3.872983346207417*hamil0R[2]*BstarZdBmagR[6]+(8.660254037844386*BstarZdBmagR[4]-6.708203932499369*BstarZdBmagR[1]+3.872983346207417*BstarZdBmagR[0])*hamil0R[5]+hamil0R[2]*(1.732050807568877*BstarZdBmagR[2]-3.0*BstarZdBmagR[3]))*rdvpar2R)/m_; 
  alpha0R[2] = (0.7071067811865475*hamil0R[5]*(3.872983346207417*BstarZdBmagR[6]-3.0*BstarZdBmagR[3]+1.732050807568877*BstarZdBmagR[2])*rdvpar2R)/m_; 

  double alpha1R[3]; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.125*(8.660254037844387*hamil0R[5]*BstarZdBmagR[6]+(3.872983346207417*BstarZdBmagR[2]-6.708203932499369*BstarZdBmagR[3])*hamil0R[5]+3.872983346207417*hamil0R[2]*BstarZdBmagR[4]+(1.732050807568877*BstarZdBmagR[0]-3.0*BstarZdBmagR[1])*hamil0R[2])*rdvpar2R)/m_; 

  double incr[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double fUpOrd[3];
  if (0.6324555320336759*alpha0R[2]-0.9486832980505137*alpha0R[1]+0.7071067811865475*alpha0R[0] > 0) { 
    fUpOrd[0] = 0.7745966692414833*f1L[7]-1.5*f1L[6]+0.4472135954999579*f1L[5]+1.118033988749895*f1L[4]-1.161895003862225*f1L[3]-0.6708203932499369*f1L[2]+0.8660254037844386*f1L[1]+0.5*f1L[0]; 
  } else { 
    fUpOrd[0] = (-0.7745966692414833*f1R[7])-1.5*f1R[6]+0.4472135954999579*f1R[5]+1.118033988749895*f1R[4]+1.161895003862225*f1R[3]-0.6708203932499369*f1R[2]-0.8660254037844386*f1R[1]+0.5*f1R[0]; 
  } 
  if (0.7071067811865475*alpha0R[0]-0.7905694150420947*alpha0R[2] > 0) { 
    fUpOrd[1] = (-0.9682458365518543*f1L[7])-0.5590169943749475*f1L[5]+1.118033988749895*f1L[4]+0.8660254037844386*f1L[1]+0.5*f1L[0]; 
  } else { 
    fUpOrd[1] = 0.9682458365518543*f1R[7]-0.5590169943749475*f1R[5]+1.118033988749895*f1R[4]-0.8660254037844386*f1R[1]+0.5*f1R[0]; 
  } 
  if (0.6324555320336759*alpha0R[2]+0.9486832980505137*alpha0R[1]+0.7071067811865475*alpha0R[0] > 0) { 
    fUpOrd[2] = 0.7745966692414833*f1L[7]+1.5*f1L[6]+0.4472135954999579*f1L[5]+1.118033988749895*f1L[4]+1.161895003862225*f1L[3]+0.6708203932499369*f1L[2]+0.8660254037844386*f1L[1]+0.5*f1L[0]; 
  } else { 
    fUpOrd[2] = (-0.7745966692414833*f1R[7])+1.5*f1R[6]+0.4472135954999579*f1R[5]+1.118033988749895*f1R[4]-1.161895003862225*f1R[3]+0.6708203932499369*f1R[2]-0.8660254037844386*f1R[1]+0.5*f1R[0]; 
  } 

  double fUp[3];
  fUp[0] = 0.3928371006591929*fUpOrd[2]+0.6285393610547088*fUpOrd[1]+0.3928371006591929*fUpOrd[0]; 
  fUp[1] = 0.5270462766947305*fUpOrd[2]-0.5270462766947305*fUpOrd[0]; 
  fUp[2] = 0.3513641844631533*fUpOrd[2]-0.7027283689263067*fUpOrd[1]+0.3513641844631533*fUpOrd[0]; 

  incr[0] = 0.5*(alpha0R[2]*fUp[2]+alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[1] = -0.8660254037844386*(alpha0R[2]*fUp[2]+alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[2] = 0.1*(4.47213595499958*(alpha0R[1]*fUp[2]+fUp[1]*alpha0R[2])+5.0*(alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1])); 
  incr[3] = -0.1*(7.745966692414834*(alpha0R[1]*fUp[2]+fUp[1]*alpha0R[2])+8.660254037844386*(alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1])); 
  incr[4] = 1.118033988749895*(alpha0R[2]*fUp[2]+alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[5] = 0.01428571428571429*((22.3606797749979*alpha0R[2]+35.0*alpha0R[0])*fUp[2]+35.0*fUp[0]*alpha0R[2]+31.30495168499706*alpha0R[1]*fUp[1]); 
  incr[6] = 0.03333333333333333*(30.0*(alpha0R[1]*fUp[2]+fUp[1]*alpha0R[2])+33.54101966249684*(alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1])); 
  incr[7] = -0.01428571428571429*((38.72983346207417*alpha0R[2]+60.62177826491071*alpha0R[0])*fUp[2]+60.62177826491071*fUp[0]*alpha0R[2]+54.22176684690384*alpha0R[1]*fUp[1]); 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += incr[3]*rdx2L; 
  outL[4] += -1.0*incr[4]*rdx2L; 
  outL[5] += -1.0*incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += incr[7]*rdx2L; 

  // alpha1 == 0, so nothing to do for alpha1*f0 and alpha1*f1 terms.
  return std::abs(alphaSurfAvgR); 
} 
double DeltaFGyrokineticGenGeoSurf1x1vSer_vpar_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 
  hamil0R[5] = (0.2981423969999719*m_)/rdvpar2SqR; 

  double hamil1R[8]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 
  hamil1R[4] = 1.414213562373095*phi[2]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R*wvparR+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmagR[1] = (0.2*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2R*wvparR+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmagR[2] = ((2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = ((b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (0.02857142857142857*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R*wvparR+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmagR[6] = (1.0*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2R)/(q_*rdvpar2R); 

  double alpha0R[3]; 

  double alpha1R[3]; 
  alpha1R[0] = (0.3535533905932737*((6.708203932499369*BstarZdBmagR[3]-3.872983346207417*BstarZdBmagR[1])*hamil1R[4]+hamil1R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]))*rdx2R)/m_; 
  alpha1R[1] = (0.07071067811865474*(hamil1R[4]*(30.0*BstarZdBmagR[6]-17.32050807568877*BstarZdBmagR[4]+33.54101966249685*BstarZdBmagR[2]-19.36491673103709*BstarZdBmagR[0])+hamil1R[1]*(15.0*BstarZdBmagR[3]-8.660254037844386*BstarZdBmagR[1]))*rdx2R)/m_; 
  alpha1R[2] = (0.07071067811865474*(15.0*hamil1R[1]*BstarZdBmagR[6]+(30.0*BstarZdBmagR[3]-17.32050807568877*BstarZdBmagR[1])*hamil1R[4]-8.660254037844386*hamil1R[1]*BstarZdBmagR[4])*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.125*((6.708203932499369*BstarZdBmagR[3]-3.872983346207417*BstarZdBmagR[1])*hamil1R[4]+3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamil1R[1])*rdx2R)/m_; 

  double incr[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double fUpOrd[3];
  if (0.6324555320336759*alpha1R[2]-0.9486832980505137*alpha1R[1]+0.7071067811865475*alpha1R[0] > 0) { 
    fUpOrd[0] = (-1.5*f1L[7])+0.7745966692414833*f1L[6]+1.118033988749895*f1L[5]+0.4472135954999579*f1L[4]-1.161895003862225*f1L[3]+0.8660254037844386*f1L[2]-0.6708203932499369*f1L[1]+0.5*f1L[0]; 
  } else { 
    fUpOrd[0] = (-1.5*f1R[7])-0.7745966692414833*f1R[6]+1.118033988749895*f1R[5]+0.4472135954999579*f1R[4]+1.161895003862225*f1R[3]-0.8660254037844386*f1R[2]-0.6708203932499369*f1R[1]+0.5*f1R[0]; 
  } 
  if (0.7071067811865475*alpha1R[0]-0.7905694150420947*alpha1R[2] > 0) { 
    fUpOrd[1] = (-0.9682458365518543*f1L[6])+1.118033988749895*f1L[5]-0.5590169943749475*f1L[4]+0.8660254037844386*f1L[2]+0.5*f1L[0]; 
  } else { 
    fUpOrd[1] = 0.9682458365518543*f1R[6]+1.118033988749895*f1R[5]-0.5590169943749475*f1R[4]-0.8660254037844386*f1R[2]+0.5*f1R[0]; 
  } 
  if (0.6324555320336759*alpha1R[2]+0.9486832980505137*alpha1R[1]+0.7071067811865475*alpha1R[0] > 0) { 
    fUpOrd[2] = 1.5*f1L[7]+0.7745966692414833*f1L[6]+1.118033988749895*f1L[5]+0.4472135954999579*f1L[4]+1.161895003862225*f1L[3]+0.8660254037844386*f1L[2]+0.6708203932499369*f1L[1]+0.5*f1L[0]; 
  } else { 
    fUpOrd[2] = 1.5*f1R[7]-0.7745966692414833*f1R[6]+1.118033988749895*f1R[5]+0.4472135954999579*f1R[4]-1.161895003862225*f1R[3]-0.8660254037844386*f1R[2]+0.6708203932499369*f1R[1]+0.5*f1R[0]; 
  } 

  double fUp[3];
  fUp[0] = 0.3928371006591929*fUpOrd[2]+0.6285393610547088*fUpOrd[1]+0.3928371006591929*fUpOrd[0]; 
  fUp[1] = 0.5270462766947305*fUpOrd[2]-0.5270462766947305*fUpOrd[0]; 
  fUp[2] = 0.3513641844631533*fUpOrd[2]-0.7027283689263067*fUpOrd[1]+0.3513641844631533*fUpOrd[0]; 

  incr[0] = 0.5*(alpha1R[2]*fUp[2]+alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.1*(4.47213595499958*(alpha1R[1]*fUp[2]+fUp[1]*alpha1R[2])+5.0*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 
  incr[2] = -0.8660254037844386*(alpha1R[2]*fUp[2]+alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = -0.1*(7.745966692414834*(alpha1R[1]*fUp[2]+fUp[1]*alpha1R[2])+8.660254037844386*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 
  incr[4] = 0.01428571428571429*((22.3606797749979*alpha1R[2]+35.0*alpha1R[0])*fUp[2]+35.0*fUp[0]*alpha1R[2]+31.30495168499706*alpha1R[1]*fUp[1]); 
  incr[5] = 1.118033988749895*(alpha1R[2]*fUp[2]+alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[6] = -0.01428571428571429*((38.72983346207417*alpha1R[2]+60.62177826491071*alpha1R[0])*fUp[2]+60.62177826491071*fUp[0]*alpha1R[2]+54.22176684690384*alpha1R[1]*fUp[1]); 
  incr[7] = 0.03333333333333333*(30.0*(alpha1R[1]*fUp[2]+fUp[1]*alpha1R[2])+33.54101966249684*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 
  outL[4] += -1.0*incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += -1.0*incr[7]*rdvpar2L; 

  // linear term alpha1*f0 
  if (0.6324555320336759*alpha1R[2]-0.9486832980505137*alpha1R[1]+0.7071067811865475*alpha1R[0] > 0) { 
    fUpOrd[0] = (-1.5*f0L[7])+0.7745966692414833*f0L[6]+1.118033988749895*f0L[5]+0.4472135954999579*f0L[4]-1.161895003862225*f0L[3]+0.8660254037844386*f0L[2]-0.6708203932499369*f0L[1]+0.5*f0L[0]; 
  } else { 
    fUpOrd[0] = (-1.5*f0R[7])-0.7745966692414833*f0R[6]+1.118033988749895*f0R[5]+0.4472135954999579*f0R[4]+1.161895003862225*f0R[3]-0.8660254037844386*f0R[2]-0.6708203932499369*f0R[1]+0.5*f0R[0]; 
  } 
  if (0.7071067811865475*alpha1R[0]-0.7905694150420947*alpha1R[2] > 0) { 
    fUpOrd[1] = (-0.9682458365518543*f0L[6])+1.118033988749895*f0L[5]-0.5590169943749475*f0L[4]+0.8660254037844386*f0L[2]+0.5*f0L[0]; 
  } else { 
    fUpOrd[1] = 0.9682458365518543*f0R[6]+1.118033988749895*f0R[5]-0.5590169943749475*f0R[4]-0.8660254037844386*f0R[2]+0.5*f0R[0]; 
  } 
  if (0.6324555320336759*alpha1R[2]+0.9486832980505137*alpha1R[1]+0.7071067811865475*alpha1R[0] > 0) { 
    fUpOrd[2] = 1.5*f0L[7]+0.7745966692414833*f0L[6]+1.118033988749895*f0L[5]+0.4472135954999579*f0L[4]+1.161895003862225*f0L[3]+0.8660254037844386*f0L[2]+0.6708203932499369*f0L[1]+0.5*f0L[0]; 
  } else { 
    fUpOrd[2] = 1.5*f0R[7]-0.7745966692414833*f0R[6]+1.118033988749895*f0R[5]+0.4472135954999579*f0R[4]-1.161895003862225*f0R[3]-0.8660254037844386*f0R[2]+0.6708203932499369*f0R[1]+0.5*f0R[0]; 
  } 

  fUp[0] = 0.3928371006591929*fUpOrd[2]+0.6285393610547088*fUpOrd[1]+0.3928371006591929*fUpOrd[0]; 
  fUp[1] = 0.5270462766947305*fUpOrd[2]-0.5270462766947305*fUpOrd[0]; 
  fUp[2] = 0.3513641844631533*fUpOrd[2]-0.7027283689263067*fUpOrd[1]+0.3513641844631533*fUpOrd[0]; 

  incr[0] = 0.5*(alpha1R[2]*fUp[2]+alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.1*(4.47213595499958*(alpha1R[1]*fUp[2]+fUp[1]*alpha1R[2])+5.0*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 
  incr[2] = -0.8660254037844386*(alpha1R[2]*fUp[2]+alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = -0.1*(7.745966692414834*(alpha1R[1]*fUp[2]+fUp[1]*alpha1R[2])+8.660254037844386*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 
  incr[4] = 0.01428571428571429*((22.3606797749979*alpha1R[2]+35.0*alpha1R[0])*fUp[2]+35.0*fUp[0]*alpha1R[2]+31.30495168499706*alpha1R[1]*fUp[1]); 
  incr[5] = 1.118033988749895*(alpha1R[2]*fUp[2]+alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[6] = -0.01428571428571429*((38.72983346207417*alpha1R[2]+60.62177826491071*alpha1R[0])*fUp[2]+60.62177826491071*fUp[0]*alpha1R[2]+54.22176684690384*alpha1R[1]*fUp[1]); 
  incr[7] = 0.03333333333333333*(30.0*(alpha1R[1]*fUp[2]+fUp[1]*alpha1R[2])+33.54101966249684*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1])); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 
  outL[4] += -1.0*incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += -1.0*incr[7]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
