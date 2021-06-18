#include <gyrofluid_mod_decl.h>

double gyrofluid_surf_1x_p2_ser_x(const double q_, const double m_, const double *wL1, const double *dxL1, const double *wR1, const double *dxR1, const double cMaxIn, const double *jacL, const double *rBmagL, const double *jacDbmagL, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, const double *csL1, const double *csR1, double *outL, double *outR) 
{ 
  // q_,m_:              species charge and mass.
  // wL,wR:              cell-center in left and right cells.
  // dxL,dxR:            cell length in left and right cells.
  // cMaxIn:             maximum sound speed (or some factor like it).
  // jac:                jacobian.
  // rBmag:              reciprocal of magnetic field magnitude (1/B).
  // jacDbmag:           jacobian divided by B (J/B).
  // sMomL,sMomR:        stepped moments (times Jacobian) in left and right cells.
  // phiL,phiR:          electrostatic potential in left and right cells.
  // primMomL,primMomR:  primitive moments (upar, Tpar, Tperp) in left and right cells.
  // csL,csR:            sound speed in left and right cells.
  // outL/outR:          output increment in left and right cells.

  double wxL = wL1[0];
  double rdx2L = 2.0/dxL1[0];
  double rdxSq4L = rdx2L*rdx2L;
  double wxR = wR1[0];
  double rdx2R = 2.0/dxR1[0];
  double rdxSq4R = rdx2R*rdx2R;

  double uparL[1]; 
  uparL[0] = 0.7071067811865475*(2.23606797749979*primMomL1[2]+1.732050807568877*primMomL1[1]+primMomL1[0]); 

  double uparR[1]; 
  uparR[0] = 0.7071067811865475*(2.23606797749979*primMomR1[2]-1.732050807568877*primMomR1[1]+primMomR1[0]); 

  double cMax = cMaxIn + fabs(fmax(1.58113883008419*primMomL1[2]+1.224744871391589*primMomL1[1]+0.7071067811865475*primMomL1[0],1.58113883008419*primMomR1[2]-1.224744871391589*primMomR1[1]+0.7071067811865475*primMomR1[0]));
  double sMom1Favg[1];
  sMom1Favg[0] = 0.25*(3.16227766016838*(sMomR1[5]+sMomL1[5])-2.449489742783178*sMomR1[4]+2.449489742783178*sMomL1[4]+1.414213562373095*(sMomR1[3]+sMomL1[3])); 

  double momHat1[1];
  momHat1[0] = -0.3535533905932737*((2.23606797749979*sMomR1[2]-2.23606797749979*sMomL1[2]-1.732050807568877*(sMomR1[1]+sMomL1[1])+sMomR1[0]-1.0*sMomL1[0])*cMax-2.828427124746191*sMom1Favg[0]); 

  double sMom2Favg[1];
  sMom2Favg[0] = (0.01428571428571429*((110.6797181058933*(sMomR1[8]+sMomL1[8])-85.73214099741124*sMomR1[7]+85.73214099741124*sMomL1[7]+49.49747468305833*(sMomR1[6]+sMomL1[6]))*m_+((-85.0*sMomR1[2])+54.22176684690384*sMomR1[1]-78.26237921249266*sMomR1[0])*primMomR1[8]+((-85.0*sMomL1[2])-54.22176684690384*sMomL1[1]-78.26237921249266*sMomL1[0])*primMomL1[8]+(54.22176684690384*sMomR1[2]-105.0*sMomR1[1]+60.6217782649107*sMomR1[0])*primMomR1[7]+((-54.22176684690384*sMomL1[2])-105.0*sMomL1[1]-60.6217782649107*sMomL1[0])*primMomL1[7]+((-78.26237921249266*sMomR1[2])+60.6217782649107*sMomR1[1]-35.0*sMomR1[0])*primMomR1[6]+((-78.26237921249266*sMomL1[2])-60.6217782649107*sMomL1[1]-35.0*sMomL1[0])*primMomL1[6]))/m_; 

  double momHat2[1];
  momHat2[0] = -0.3535533905932737*((2.23606797749979*sMomR1[5]-2.23606797749979*sMomL1[5]-1.732050807568877*(sMomR1[4]+sMomL1[4])+sMomR1[3]-1.0*sMomL1[3])*cMax-2.828427124746191*sMom2Favg[0]); 

  double sMom3Favg[1];
  sMom3Favg[0] = (5.102040816326531e-4*(((1190.0*primMomR1[2]-759.1047358566536*primMomR1[1]+1095.673308974897*primMomR1[0])*sMomR1[8]+(1190.0*primMomL1[2]+759.1047358566536*primMomL1[1]+1095.673308974897*primMomL1[0])*sMomL1[8]+((-759.1047358566536*primMomR1[2])+1470.0*primMomR1[1]-848.7048957087499*primMomR1[0])*sMomR1[7]+(759.1047358566536*primMomL1[2]+1470.0*primMomL1[1]+848.7048957087499*primMomL1[0])*sMomL1[7]+(1095.673308974897*primMomR1[2]-848.7048957087499*primMomR1[1]+490.0*primMomR1[0])*sMomR1[6]+(1095.673308974897*primMomL1[2]+848.7048957087499*primMomL1[1]+490.0*primMomL1[0])*sMomL1[6])*m_+((1312.345228969878*primMomR1[2]-943.0535509715237*primMomR1[1]+841.4570696119916*primMomR1[0])*sMomR1[2]+(841.4570696119916*sMomR1[0]-480.099989585503*sMomR1[1])*primMomR1[2]+(929.7096320895039*primMomR1[1]-536.768106355063*primMomR1[0])*sMomR1[1]+sMomR1[0]*(774.7580267412532*primMomR1[0]-536.768106355063*primMomR1[1]))*primMomR1[5]+((1312.345228969878*primMomL1[2]+943.0535509715237*primMomL1[1]+841.4570696119916*primMomL1[0])*sMomL1[2]+(480.099989585503*sMomL1[1]+841.4570696119916*sMomL1[0])*primMomL1[2]+(929.7096320895039*primMomL1[1]+536.768106355063*primMomL1[0])*sMomL1[1]+sMomL1[0]*(536.768106355063*primMomL1[1]+774.7580267412532*primMomL1[0]))*primMomL1[5]+(((-480.099989585503*primMomR1[2])+929.7096320895039*primMomR1[1]-536.768106355063*primMomR1[0])*sMomR1[2]+(1527.380109861327*sMomR1[1]-536.768106355063*sMomR1[0])*primMomR1[2]+(1039.446968344225*primMomR1[0]-1080.224976567381*primMomR1[1])*sMomR1[1]+sMomR1[0]*(1039.446968344225*primMomR1[1]-600.1249869818786*primMomR1[0]))*primMomR1[4]+((480.099989585503*primMomL1[2]+929.7096320895039*primMomL1[1]+536.768106355063*primMomL1[0])*sMomL1[2]+(1527.380109861327*sMomL1[1]+536.768106355063*sMomL1[0])*primMomL1[2]+(1080.224976567381*primMomL1[1]+1039.446968344225*primMomL1[0])*sMomL1[1]+sMomL1[0]*(1039.446968344225*primMomL1[1]+600.1249869818786*primMomL1[0]))*primMomL1[4]+((841.4570696119916*primMomR1[2]-536.768106355063*primMomR1[1]+774.7580267412532*primMomR1[0])*sMomR1[2]+(774.7580267412532*sMomR1[0]-536.768106355063*sMomR1[1])*primMomR1[2]+(1039.446968344225*primMomR1[1]-600.1249869818786*primMomR1[0])*sMomR1[1]+sMomR1[0]*(346.4823227814083*primMomR1[0]-600.1249869818786*primMomR1[1]))*primMomR1[3]+((841.4570696119916*primMomL1[2]+536.768106355063*primMomL1[1]+774.7580267412532*primMomL1[0])*sMomL1[2]+(536.768106355063*sMomL1[1]+774.7580267412532*sMomL1[0])*primMomL1[2]+(1039.446968344225*primMomL1[1]+600.1249869818786*primMomL1[0])*sMomL1[1]+sMomL1[0]*(600.1249869818786*primMomL1[1]+346.4823227814083*primMomL1[0]))*primMomL1[3]))/m_; 

  double momHat3[1];
  momHat3[0] = -0.3535533905932737*((2.23606797749979*sMomR1[8]-2.23606797749979*sMomL1[8]-1.732050807568877*(sMomR1[7]+sMomL1[7])+sMomR1[6]-1.0*sMomL1[6])*cMax-2.828427124746191*sMom3Favg[0]); 

  double sMom4Favg[1];
  sMom4Favg[0] = 0.007142857142857143*((85.0*primMomR1[2]-54.22176684690384*primMomR1[1]+78.26237921249266*primMomR1[0])*sMomR1[11]+(85.0*primMomL1[2]+54.22176684690384*primMomL1[1]+78.26237921249266*primMomL1[0])*sMomL1[11]+((-54.22176684690384*primMomR1[2])+105.0*primMomR1[1]-60.6217782649107*primMomR1[0])*sMomR1[10]+(54.22176684690384*primMomL1[2]+105.0*primMomL1[1]+60.6217782649107*primMomL1[0])*sMomL1[10]+(78.26237921249266*primMomR1[2]-60.6217782649107*primMomR1[1]+35.0*primMomR1[0])*sMomR1[9]+(78.26237921249266*primMomL1[2]+60.6217782649107*primMomL1[1]+35.0*primMomL1[0])*sMomL1[9]); 

  double momHat4[1];
  momHat4[0] = -0.3535533905932737*((2.23606797749979*sMomR1[11]-2.23606797749979*sMomL1[11]-1.732050807568877*(sMomR1[10]+sMomL1[10])+sMomR1[9]-1.0*sMomL1[9])*cMax-2.828427124746191*sMom4Favg[0]); 

  double incr1[3];
  incr1[0] = 0.7071067811865475*momHat1[0]; 
  incr1[1] = -1.224744871391589*momHat1[0]; 
  incr1[2] = 1.58113883008419*momHat1[0]; 

  double incr2[3];
  incr2[0] = 0.7071067811865475*momHat2[0]; 
  incr2[1] = -1.224744871391589*momHat2[0]; 
  incr2[2] = 1.58113883008419*momHat2[0]; 

  double incr3[3];
  incr3[0] = 0.7071067811865475*momHat3[0]; 
  incr3[1] = -1.224744871391589*momHat3[0]; 
  incr3[2] = 1.58113883008419*momHat3[0]; 

  double incr4[3];
  incr4[0] = 0.7071067811865475*momHat4[0]; 
  incr4[1] = -1.224744871391589*momHat4[0]; 
  incr4[2] = 1.58113883008419*momHat4[0]; 

  outR[0] += incr1[0]*rdx2R; 
  outR[1] += incr1[1]*rdx2R; 
  outR[2] += incr1[2]*rdx2R; 

  outL[0] += -1.0*incr1[0]*rdx2L; 
  outL[1] += incr1[1]*rdx2L; 
  outL[2] += -1.0*incr1[2]*rdx2L; 

  outR[3] += incr2[0]*rdx2R; 
  outR[4] += incr2[1]*rdx2R; 
  outR[5] += incr2[2]*rdx2R; 

  outL[3] += -1.0*incr2[0]*rdx2L; 
  outL[4] += incr2[1]*rdx2L; 
  outL[5] += -1.0*incr2[2]*rdx2L; 

  outR[6] += incr3[0]*rdx2R; 
  outR[7] += incr3[1]*rdx2R; 
  outR[8] += incr3[2]*rdx2R; 

  outL[6] += -1.0*incr3[0]*rdx2L; 
  outL[7] += incr3[1]*rdx2L; 
  outL[8] += -1.0*incr3[2]*rdx2L; 

  outR[9] += incr4[0]*rdx2R; 
  outR[10] += incr4[1]*rdx2R; 
  outR[11] += incr4[2]*rdx2R; 

  outL[9] += -1.0*incr4[0]*rdx2L; 
  outL[10] += incr4[1]*rdx2L; 
  outL[11] += -1.0*incr4[2]*rdx2L; 

  return 0.7071067811865475*csL1[0]-0.7905694150420947*csL1[2]; 
}
