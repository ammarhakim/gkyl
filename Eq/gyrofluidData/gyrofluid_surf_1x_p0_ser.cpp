#include <gyrofluid_mod_decl.h>

double gyrofluid_surf_1x_p0_ser_x(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *dxL1, const double *wR1, const double *dxR1, const double uMaxIn, const double *jacL, const double *rBmagL, const double *jacDbmagL, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR) 
{ 
  // q_,m_:              species charge and mass.
  // kappaPar,kappaPerp: heat conductivity coefficients.
  // kperpSq:            k_perp^2.
  // wL,wR:              cell-center in left and right cells.
  // dxL,dxR:            cell length in left and right cells.
  // uMaxIn:             maximum speed.
  // jac:                jacobian.
  // rBmag:              reciprocal of magnetic field magnitude (1/B).
  // jacDbmag:           jacobian divided by B (J/B).
  // sMomL,sMomR:        stepped moments (times Jacobian) in left and right cells.
  // phiL,phiR:          electrostatic potential in left and right cells.
  // primMomL,primMomR:  primitive moments (upar, Tpar, Tperp) in left and right cells.
  // outL/outR:          output increment in left and right cells.

  double wxL = wL1[0];
  double rdx2L = 2.0/dxL1[0];
  double rdxSq4L = rdx2L*rdx2L;
  double wxR = wR1[0];
  double rdx2R = 2.0/dxR1[0];
  double rdxSq4R = rdx2R*rdx2R;

  double uparL[1]; 
  uparL[0] = 0.7071067811865475*primMomL1[0]; 

  double uparR[1]; 
  uparR[0] = 0.7071067811865475*primMomR1[0]; 

  double sMom1Favg[1];
  sMom1Favg[0] = 0.3535533905932738*(sMomR1[1]+sMomL1[1]); 

  double momHat1[1];
  momHat1[0] = -0.3535533905932737*((sMomR1[0]-1.0*sMomL1[0])*uMaxIn-2.828427124746191*sMom1Favg[0]); 

  double sMom2Favg[1];
  sMom2Favg[0] = (0.5*(1.414213562373095*(sMomR1[2]+sMomL1[2])*m_-1.0*(sMomR1[0]*primMomR1[2]+sMomL1[0]*primMomL1[2])))/m_; 

  double momHat2[1];
  momHat2[0] = -0.3535533905932737*((sMomR1[1]-1.0*sMomL1[1])*uMaxIn-2.828427124746191*sMom2Favg[0]); 

  double sMom3Favg[1];
  sMom3Favg[0] = (0.125*(2.0*(primMomR1[0]*sMomR1[2]+primMomL1[0]*sMomL1[2])*m_+1.414213562373095*(primMomR1[0]*sMomR1[0]*primMomR1[1]+primMomL1[0]*sMomL1[0]*primMomL1[1])))/m_; 

  double momHat3[1];
  momHat3[0] = -0.3535533905932737*((sMomR1[2]-1.0*sMomL1[2])*uMaxIn-2.828427124746191*sMom3Favg[0]); 

  double sMom4Favg[1];
  sMom4Favg[0] = 0.25*(primMomR1[0]*sMomR1[3]+primMomL1[0]*sMomL1[3]); 

  double momHat4[1];
  momHat4[0] = -0.3535533905932737*((sMomR1[3]-1.0*sMomL1[3])*uMaxIn-2.828427124746191*sMom4Favg[0]); 

  double GheatF3[1];

  GheatF3[0] = -0.125*((1.414213562373095*jacL[0]*phiR1[0]-1.414213562373095*jacL[0]*phiL1[0])*kappaPerp*kperpSq*q_+(2.828427124746191*jacL[0]*primMomL1[2]-2.828427124746191*jacL[0]*primMomR1[2])*kappaPerp+(1.414213562373095*jacL[0]*primMomL1[1]-1.414213562373095*jacL[0]*primMomR1[1])*kappaPar); 

  double GheatF4[1];

  GheatF4[0] = -0.125*((1.414213562373095*jacDbmagL[0]*phiR1[0]-1.414213562373095*jacDbmagL[0]*phiL1[0])*kappaPerp*kperpSq*q_+(2.828427124746191*jacDbmagL[0]*primMomL1[2]-2.828427124746191*jacDbmagL[0]*primMomR1[2])*kappaPerp); 

  double incr1[1];
  incr1[0] = 0.7071067811865475*momHat1[0]; 

  double incrNonFlux1[1];

  double incr2[1];
  incr2[0] = 0.7071067811865475*momHat2[0]; 

  double incrNonFlux2[1];

  double incr3[1];
  incr3[0] = 0.7071067811865475*momHat3[0]-0.5*GheatF3[0]; 

  double incrNonFlux3[1];

  double incr4[1];
  incr4[0] = 0.7071067811865475*momHat4[0]-0.5*GheatF4[0]; 

  double incrNonFlux4[1];

  outR[0] += incr1[0]*rdx2R; 

  outL[0] += -1.0*incr1[0]*rdx2L; 

  outR[1] += incr2[0]*rdx2R; 

  outL[1] += -1.0*incr2[0]*rdx2L; 

  outR[2] += incr3[0]*rdx2R; 

  outL[2] += -1.0*incr3[0]*rdx2L; 

  outR[3] += incr4[0]*rdx2R; 

  outL[3] += -1.0*incr4[0]*rdx2L; 

  return fabs(0.4854917717073234*sqrt((2.0*primMomL1[2])/m_+primMomL1[1]/m_)+0.7071067811865475*primMomL1[0]); 
}
