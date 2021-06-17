#include <gyrofluid_heatflux_mod_decl.h>

double gyrofluid_heatflux_surf_1x_p1_ser_x(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *dxL1, const double *wR1, const double *dxR1, const double cMaxIn, const double *jacL, const double *jacDbmagL, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR) 
{ 
  // q_,m_:              species charge and mass.
  // kappaPar,kappaPerp: heat conductivity coefficients.
  // kperpSq:            k_perp^2.
  // wL,wR:              cell-center in left and right cells.
  // dxL,dxR:            cell length in left and right cells.
  // cMaxIn:             maximum sound speed (or some factor like it).
  // jac:                jacobian.
  // jacDbmag:           jacobian divided by B (J/B).
  // sMomL,sMomR:        stepped moments (times Jacobian) in left and right cells.
  // phiL,phiR:          electrostatic potential in left and right cells.
  // primMomL,primMomR:  primitive moments (upar, Tpar, Tperp) in left and right cells.
  // outL/outR:          output increment in left and right cells.

  double rdx2L = 2.0/dxL1[0];
  double rdxSq4L = rdx2L*rdx2L;
  double rdx2R = 2.0/dxR1[0];
  double rdxSq4R = rdx2R*rdx2R;

  double GheatF3[2];

  GheatF3[0] = 0.03125*(((21.21320343559643*jacL[1]+12.24744871391589*jacL[0])*phiR1[1]+(21.21320343559643*jacL[1]+12.24744871391589*jacL[0])*phiL1[1]+(22.0454076850486*phiL1[0]-22.0454076850486*phiR1[0])*jacL[1]-12.72792206135786*jacL[0]*phiR1[0]+12.72792206135786*jacL[0]*phiL1[0])*kappaPerp*kperpSq*q_+(((-42.42640687119286*jacL[1])-24.49489742783179*jacL[0])*primMomR1[5]+((-42.42640687119286*jacL[1])-24.49489742783179*jacL[0])*primMomL1[5]+(44.0908153700972*jacL[1]+25.45584412271572*jacL[0])*primMomR1[4]+((-44.0908153700972*jacL[1])-25.45584412271572*jacL[0])*primMomL1[4])*kappaPerp+(((-21.21320343559643*jacL[1])-12.24744871391589*jacL[0])*primMomR1[3]+((-21.21320343559643*jacL[1])-12.24744871391589*jacL[0])*primMomL1[3]+(22.0454076850486*jacL[1]+12.72792206135786*jacL[0])*primMomR1[2]+((-22.0454076850486*jacL[1])-12.72792206135786*jacL[0])*primMomL1[2])*kappaPar)*rdx2L; 

  double GheatF4[2];

  GheatF4[0] = 0.03125*(((21.21320343559643*jacDbmagL[1]+12.24744871391589*jacDbmagL[0])*phiR1[1]+(21.21320343559643*jacDbmagL[1]+12.24744871391589*jacDbmagL[0])*phiL1[1]+(22.0454076850486*phiL1[0]-22.0454076850486*phiR1[0])*jacDbmagL[1]-12.72792206135786*jacDbmagL[0]*phiR1[0]+12.72792206135786*jacDbmagL[0]*phiL1[0])*kappaPerp*kperpSq*q_+(((-42.42640687119286*jacDbmagL[1])-24.49489742783179*jacDbmagL[0])*primMomR1[5]+((-42.42640687119286*jacDbmagL[1])-24.49489742783179*jacDbmagL[0])*primMomL1[5]+(44.0908153700972*jacDbmagL[1]+25.45584412271572*jacDbmagL[0])*primMomR1[4]+((-44.0908153700972*jacDbmagL[1])-25.45584412271572*jacDbmagL[0])*primMomL1[4])*kappaPerp)*rdx2L; 

  double incr3[2];
  incr3[0] = -0.5*GheatF3[0]; 
  incr3[1] = 0.8660254037844386*GheatF3[0]; 

  double incrNonFlux3[2];
  incrNonFlux3[1] = 0.0625*(kappaPerp*(((4.898979485566357*jacL[1]+2.828427124746191*jacL[0])*phiR1[1]+((-4.898979485566357*jacL[1])-2.828427124746191*jacL[0])*phiL1[1]+(phiR1[0]+phiL1[0])*((-4.242640687119286*jacL[1])-2.449489742783178*jacL[0]))*kperpSq*q_+((-9.797958971132715*jacL[1])-5.656854249492382*jacL[0])*primMomR1[5]+(9.797958971132715*jacL[1]+5.656854249492382*jacL[0])*primMomL1[5]+(8.485281374238571*jacL[1]+4.898979485566357*jacL[0])*(primMomR1[4]+primMomL1[4]))+(((-4.898979485566357*jacL[1])-2.828427124746191*jacL[0])*primMomR1[3]+(4.898979485566357*jacL[1]+2.828427124746191*jacL[0])*primMomL1[3]+(4.242640687119286*jacL[1]+2.449489742783178*jacL[0])*(primMomR1[2]+primMomL1[2]))*kappaPar); 

  double incr4[2];
  incr4[0] = -0.5*GheatF4[0]; 
  incr4[1] = 0.8660254037844386*GheatF4[0]; 

  double incrNonFlux4[2];
  incrNonFlux4[1] = 0.0625*kappaPerp*(((4.898979485566357*jacDbmagL[1]+2.828427124746191*jacDbmagL[0])*phiR1[1]+((-4.898979485566357*jacDbmagL[1])-2.828427124746191*jacDbmagL[0])*phiL1[1]+(phiR1[0]+phiL1[0])*((-4.242640687119286*jacDbmagL[1])-2.449489742783178*jacDbmagL[0]))*kperpSq*q_+((-9.797958971132715*jacDbmagL[1])-5.656854249492382*jacDbmagL[0])*primMomR1[5]+(9.797958971132715*jacDbmagL[1]+5.656854249492382*jacDbmagL[0])*primMomL1[5]+(8.485281374238571*jacDbmagL[1]+4.898979485566357*jacDbmagL[0])*(primMomR1[4]+primMomL1[4])); 

  outR[4] += incr3[0]*rdx2R; 
  outR[5] += incrNonFlux3[1]*rdxSq4R+incr3[1]*rdx2R; 

  outL[4] += -1.0*incr3[0]*rdx2L; 
  outL[5] += incr3[1]*rdx2L-1.0*incrNonFlux3[1]*rdxSq4L; 

  outR[6] += incr4[0]*rdx2R; 
  outR[7] += incrNonFlux4[1]*rdxSq4R+incr4[1]*rdx2R; 

  outL[6] += -1.0*incr4[0]*rdx2L; 
  outL[7] += incr4[1]*rdx2L-1.0*incrNonFlux4[1]*rdxSq4L; 

}
