#include <gyrofluid_mod_decl.h>

double gyrofluid_surf_1x_p1_ser_x(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *dxL1, const double *wR1, const double *dxR1, const double uMaxIn, const double *jacL, const double *rBmagL, const double *jacDbmagL, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR) 
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
  uparL[0] = 0.5*(2.449489742783178*primMomL1[1]+1.414213562373095*primMomL1[0]); 

  double uparR[1]; 
  uparR[0] = -0.5*(2.449489742783178*primMomR1[1]-1.414213562373095*primMomR1[0]); 

  double sMom1Favg[1];
  sMom1Favg[0] = -0.25*(2.449489742783178*sMomR1[3]-2.449489742783178*sMomL1[3]-1.414213562373095*(sMomR1[2]+sMomL1[2])); 

  double momHat1[1];
  momHat1[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[1]+sMomL1[1])-1.0*sMomR1[0]+sMomL1[0])*uMaxIn+2.828427124746191*sMom1Favg[0]); 

  double sMom2Favg[1];
  sMom2Favg[0] = -(0.5*((2.449489742783178*sMomR1[5]-2.449489742783178*sMomL1[5]-1.414213562373095*(sMomR1[4]+sMomL1[4]))*m_+(sMomR1[1]-1.732050807568877*sMomR1[0])*primMomR1[5]+(sMomL1[1]+1.732050807568877*sMomL1[0])*primMomL1[5]+(sMomR1[0]-1.732050807568877*sMomR1[1])*primMomR1[4]+(1.732050807568877*sMomL1[1]+sMomL1[0])*primMomL1[4]))/m_; 

  double momHat2[1];
  momHat2[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[3]+sMomL1[3])-1.0*sMomR1[2]+sMomL1[2])*uMaxIn+2.828427124746191*sMom2Favg[0]); 

  double sMom3Favg[1];
  sMom3Favg[0] = (0.125*(((2.0*primMomR1[1]-3.464101615137754*primMomR1[0])*sMomR1[5]+(2.0*primMomL1[1]+3.464101615137754*primMomL1[0])*sMomL1[5]+(2.0*primMomR1[0]-3.464101615137754*primMomR1[1])*sMomR1[4]+(3.464101615137754*primMomL1[1]+2.0*primMomL1[0])*sMomL1[4])*m_+((1.414213562373095*primMomR1[0]-2.449489742783178*primMomR1[1])*sMomR1[1]+sMomR1[0]*(1.414213562373095*primMomR1[1]-2.449489742783178*primMomR1[0]))*primMomR1[3]+((2.449489742783178*primMomL1[1]+1.414213562373095*primMomL1[0])*sMomL1[1]+sMomL1[0]*(1.414213562373095*primMomL1[1]+2.449489742783178*primMomL1[0]))*primMomL1[3]+((1.414213562373095*primMomR1[1]-2.449489742783178*primMomR1[0])*sMomR1[1]+sMomR1[0]*(1.414213562373095*primMomR1[0]-2.449489742783178*primMomR1[1]))*primMomR1[2]+((1.414213562373095*primMomL1[1]+2.449489742783178*primMomL1[0])*sMomL1[1]+sMomL1[0]*(2.449489742783178*primMomL1[1]+1.414213562373095*primMomL1[0]))*primMomL1[2]))/m_; 

  double momHat3[1];
  momHat3[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[5]+sMomL1[5])-1.0*sMomR1[4]+sMomL1[4])*uMaxIn+2.828427124746191*sMom3Favg[0]); 

  double sMom4Favg[1];
  sMom4Favg[0] = 0.25*((primMomR1[1]-1.732050807568877*primMomR1[0])*sMomR1[7]+(primMomL1[1]+1.732050807568877*primMomL1[0])*sMomL1[7]+(primMomR1[0]-1.732050807568877*primMomR1[1])*sMomR1[6]+(1.732050807568877*primMomL1[1]+primMomL1[0])*sMomL1[6]); 

  double momHat4[1];
  momHat4[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[7]+sMomL1[7])-1.0*sMomR1[6]+sMomL1[6])*uMaxIn+2.828427124746191*sMom4Favg[0]); 

  double GheatF3[2];

  GheatF3[0] = 0.03125*(((21.21320343559643*jacL[1]+12.24744871391589*jacL[0])*phiR1[1]+(21.21320343559643*jacL[1]+12.24744871391589*jacL[0])*phiL1[1]+(22.0454076850486*phiL1[0]-22.0454076850486*phiR1[0])*jacL[1]-12.72792206135786*jacL[0]*phiR1[0]+12.72792206135786*jacL[0]*phiL1[0])*kappaPerp*kperpSq*q_+(((-42.42640687119286*jacL[1])-24.49489742783179*jacL[0])*primMomR1[5]+((-42.42640687119286*jacL[1])-24.49489742783179*jacL[0])*primMomL1[5]+(44.0908153700972*jacL[1]+25.45584412271572*jacL[0])*primMomR1[4]+((-44.0908153700972*jacL[1])-25.45584412271572*jacL[0])*primMomL1[4])*kappaPerp+(((-21.21320343559643*jacL[1])-12.24744871391589*jacL[0])*primMomR1[3]+((-21.21320343559643*jacL[1])-12.24744871391589*jacL[0])*primMomL1[3]+(22.0454076850486*jacL[1]+12.72792206135786*jacL[0])*primMomR1[2]+((-22.0454076850486*jacL[1])-12.72792206135786*jacL[0])*primMomL1[2])*kappaPar); 

  double GheatF4[2];

  GheatF4[0] = 0.03125*(((21.21320343559643*jacDbmagL[1]+12.24744871391589*jacDbmagL[0])*phiR1[1]+(21.21320343559643*jacDbmagL[1]+12.24744871391589*jacDbmagL[0])*phiL1[1]+(22.0454076850486*phiL1[0]-22.0454076850486*phiR1[0])*jacDbmagL[1]-12.72792206135786*jacDbmagL[0]*phiR1[0]+12.72792206135786*jacDbmagL[0]*phiL1[0])*kappaPerp*kperpSq*q_+(((-42.42640687119286*jacDbmagL[1])-24.49489742783179*jacDbmagL[0])*primMomR1[5]+((-42.42640687119286*jacDbmagL[1])-24.49489742783179*jacDbmagL[0])*primMomL1[5]+(44.0908153700972*jacDbmagL[1]+25.45584412271572*jacDbmagL[0])*primMomR1[4]+((-44.0908153700972*jacDbmagL[1])-25.45584412271572*jacDbmagL[0])*primMomL1[4])*kappaPerp); 

  double incr1[2];
  incr1[0] = 0.7071067811865475*momHat1[0]; 
  incr1[1] = -1.224744871391589*momHat1[0]; 

  double incrNonFlux1[2];

  double incr2[2];
  incr2[0] = 0.7071067811865475*momHat2[0]; 
  incr2[1] = -1.224744871391589*momHat2[0]; 

  double incrNonFlux2[2];

  double incr3[2];
  incr3[0] = 0.7071067811865475*momHat3[0]-0.5*GheatF3[0]; 
  incr3[1] = 0.8660254037844386*GheatF3[0]-1.224744871391589*momHat3[0]; 

  double incrNonFlux3[2];
  incrNonFlux3[1] = 0.0625*(kappaPerp*(((4.898979485566357*jacL[1]+2.828427124746191*jacL[0])*phiR1[1]+((-4.898979485566357*jacL[1])-2.828427124746191*jacL[0])*phiL1[1]+(phiR1[0]+phiL1[0])*((-4.242640687119286*jacL[1])-2.449489742783178*jacL[0]))*kperpSq*q_+((-9.797958971132715*jacL[1])-5.656854249492382*jacL[0])*primMomR1[5]+(9.797958971132715*jacL[1]+5.656854249492382*jacL[0])*primMomL1[5]+(8.485281374238571*jacL[1]+4.898979485566357*jacL[0])*(primMomR1[4]+primMomL1[4]))+(((-4.898979485566357*jacL[1])-2.828427124746191*jacL[0])*primMomR1[3]+(4.898979485566357*jacL[1]+2.828427124746191*jacL[0])*primMomL1[3]+(4.242640687119286*jacL[1]+2.449489742783178*jacL[0])*(primMomR1[2]+primMomL1[2]))*kappaPar); 

  double incr4[2];
  incr4[0] = 0.7071067811865475*momHat4[0]-0.5*GheatF4[0]; 
  incr4[1] = 0.8660254037844386*GheatF4[0]-1.224744871391589*momHat4[0]; 

  double incrNonFlux4[2];
  incrNonFlux4[1] = 0.0625*kappaPerp*(((4.898979485566357*jacDbmagL[1]+2.828427124746191*jacDbmagL[0])*phiR1[1]+((-4.898979485566357*jacDbmagL[1])-2.828427124746191*jacDbmagL[0])*phiL1[1]+(phiR1[0]+phiL1[0])*((-4.242640687119286*jacDbmagL[1])-2.449489742783178*jacDbmagL[0]))*kperpSq*q_+((-9.797958971132715*jacDbmagL[1])-5.656854249492382*jacDbmagL[0])*primMomR1[5]+(9.797958971132715*jacDbmagL[1]+5.656854249492382*jacDbmagL[0])*primMomL1[5]+(8.485281374238571*jacDbmagL[1]+4.898979485566357*jacDbmagL[0])*(primMomR1[4]+primMomL1[4])); 

  outR[0] += incr1[0]*rdx2R; 
  outR[1] += incr1[1]*rdx2R; 

  outL[0] += -1.0*incr1[0]*rdx2L; 
  outL[1] += incr1[1]*rdx2L; 

  outR[2] += incr2[0]*rdx2R; 
  outR[3] += incr2[1]*rdx2R; 

  outL[2] += -1.0*incr2[0]*rdx2L; 
  outL[3] += incr2[1]*rdx2L; 

  outR[4] += incr3[0]*rdx2R; 
  outR[5] += incrNonFlux3[1]*rdxSq4R+incr3[1]*rdx2R; 

  outL[4] += -1.0*incr3[0]*rdx2L; 
  outL[5] += incr3[1]*rdx2L-1.0*incrNonFlux3[1]*rdxSq4L; 

  outR[6] += incr4[0]*rdx2R; 
  outR[7] += incrNonFlux4[1]*rdxSq4R+incr4[1]*rdx2R; 

  outL[6] += -1.0*incr4[0]*rdx2L; 
  outL[7] += incr4[1]*rdx2L-1.0*incrNonFlux4[1]*rdxSq4L; 

  return fabs(0.4854917717073234*sqrt((2.0*primMomL1[4])/m_+primMomL1[2]/m_)+0.7071067811865475*primMomL1[0]); 
}
