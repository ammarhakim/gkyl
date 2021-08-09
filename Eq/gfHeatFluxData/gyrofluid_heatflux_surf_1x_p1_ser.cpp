#include <gyrofluid_heatflux_mod_decl.h>

void gyrofluid_heatflux_surf_1x_p1_ser_x(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *wR1, const double *dxL1, const double *dxR1, const double cMaxIn, const double *rBmagL1, const double *rBmagR1, const double *rBmagSqL1, const double *rBmagSqR1, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR) 
{ 
  // q_,m_:              species charge and mass.
  // kappaPar,kappaPerp: heat conductivity coefficients.
  // kperpSq:            k_perp^2.
  // wL,wR:              cell-center in left and right cells.
  // dxL,dxR:            cell length in left and right cells.
  // cMaxIn:             maximum sound speed (or some factor like it).
  // rBmag:              reciprocal of magnetic field magnitude (1/B).
  // rBmagSq:            rBmag^2.
  // sMomL,sMomR:        stepped moments (times Jacobian) in left and right cells.
  // phiL,phiR:          electrostatic potential in left and right cells.
  // primMomL,primMomR:  primitive moments (upar, Tpar, Tperp) in left and right cells.
  // outL/outR:          output increment in left and right cells.

  double rdx2L = 2.0/dxL1[0];
  double rdxSq4L = rdx2L*rdx2L;
  double rdx2R = 2.0/dxR1[0];
  double rdxSq4R = rdx2R*rdx2R;

  double GheatF3[2];

  GheatF3[0] = 0.03125*(((21.21320343559643*phiR1[1]+21.21320343559643*phiL1[1]-22.0454076850486*phiR1[0]+22.0454076850486*phiL1[0])*rBmagL1[1]+12.24744871391589*rBmagL1[0]*phiR1[1]+12.24744871391589*rBmagL1[0]*phiL1[1]+(12.72792206135786*phiL1[0]-12.72792206135786*phiR1[0])*rBmagL1[0])*kappaPerp*kperpSq*q_+(((-42.42640687119286*rBmagL1[1])-24.49489742783179*rBmagL1[0])*primMomR1[5]+((-42.42640687119286*rBmagL1[1])-24.49489742783179*rBmagL1[0])*primMomL1[5]+(44.0908153700972*rBmagL1[1]+25.45584412271572*rBmagL1[0])*primMomR1[4]+((-44.0908153700972*rBmagL1[1])-25.45584412271572*rBmagL1[0])*primMomL1[4])*kappaPerp+(((-21.21320343559643*rBmagL1[1])-12.24744871391589*rBmagL1[0])*primMomR1[3]+((-21.21320343559643*rBmagL1[1])-12.24744871391589*rBmagL1[0])*primMomL1[3]+(22.0454076850486*rBmagL1[1]+12.72792206135786*rBmagL1[0])*primMomR1[2]+((-22.0454076850486*rBmagL1[1])-12.72792206135786*rBmagL1[0])*primMomL1[2])*kappaPar)*rdx2L; 

  double GheatF4[2];

  GheatF4[0] = 0.03125*(((21.21320343559643*phiR1[1]+21.21320343559643*phiL1[1]-22.0454076850486*phiR1[0]+22.0454076850486*phiL1[0])*rBmagSqL1[1]+12.24744871391589*rBmagSqL1[0]*phiR1[1]+12.24744871391589*rBmagSqL1[0]*phiL1[1]+(12.72792206135786*phiL1[0]-12.72792206135786*phiR1[0])*rBmagSqL1[0])*kappaPerp*kperpSq*q_+(((-42.42640687119286*rBmagSqL1[1])-24.49489742783179*rBmagSqL1[0])*primMomR1[5]+((-42.42640687119286*rBmagSqL1[1])-24.49489742783179*rBmagSqL1[0])*primMomL1[5]+(44.0908153700972*rBmagSqL1[1]+25.45584412271572*rBmagSqL1[0])*primMomR1[4]+((-44.0908153700972*rBmagSqL1[1])-25.45584412271572*rBmagSqL1[0])*primMomL1[4])*kappaPerp)*rdx2L; 

  double incr3[2];
  incr3[0] = -0.5*GheatF3[0]; 
  incr3[1] = 0.8660254037844386*GheatF3[0]; 

  double incrNonFlux3[2];
  incrNonFlux3[1] = 0.0625*(kappaPerp*(((4.898979485566357*phiR1[1]-4.898979485566357*phiL1[1]-4.242640687119286*(phiR1[0]+phiL1[0]))*rBmagL1[1]+rBmagL1[0]*(2.828427124746191*phiR1[1]-2.828427124746191*phiL1[1]-2.449489742783178*(phiR1[0]+phiL1[0])))*kperpSq*q_+((-9.797958971132715*rBmagL1[1])-5.656854249492382*rBmagL1[0])*primMomR1[5]+(9.797958971132715*rBmagL1[1]+5.656854249492382*rBmagL1[0])*primMomL1[5]+(8.485281374238571*rBmagL1[1]+4.898979485566357*rBmagL1[0])*(primMomR1[4]+primMomL1[4]))+(((-4.898979485566357*rBmagL1[1])-2.828427124746191*rBmagL1[0])*primMomR1[3]+(4.898979485566357*rBmagL1[1]+2.828427124746191*rBmagL1[0])*primMomL1[3]+(4.242640687119286*rBmagL1[1]+2.449489742783178*rBmagL1[0])*(primMomR1[2]+primMomL1[2]))*kappaPar); 

  double incr4[2];
  incr4[0] = -0.5*GheatF4[0]; 
  incr4[1] = 0.8660254037844386*GheatF4[0]; 

  double incrNonFlux4[2];
  incrNonFlux4[1] = 0.0625*kappaPerp*(((4.898979485566357*phiR1[1]-4.898979485566357*phiL1[1]-4.242640687119286*(phiR1[0]+phiL1[0]))*rBmagSqL1[1]+rBmagSqL1[0]*(2.828427124746191*phiR1[1]-2.828427124746191*phiL1[1]-2.449489742783178*(phiR1[0]+phiL1[0])))*kperpSq*q_+((-9.797958971132715*rBmagSqL1[1])-5.656854249492382*rBmagSqL1[0])*primMomR1[5]+(9.797958971132715*rBmagSqL1[1]+5.656854249492382*rBmagSqL1[0])*primMomL1[5]+(8.485281374238571*rBmagSqL1[1]+4.898979485566357*rBmagSqL1[0])*(primMomR1[4]+primMomL1[4])); 

  outR[4] += incr3[0]*rdx2R; 
  outR[5] += incrNonFlux3[1]*rdxSq4R+incr3[1]*rdx2R; 

  outL[4] += -1.0*incr3[0]*rdx2L; 
  outL[5] += incr3[1]*rdx2L-1.0*incrNonFlux3[1]*rdxSq4L; 

  outR[6] += incr4[0]*rdx2R; 
  outR[7] += incrNonFlux4[1]*rdxSq4R+incr4[1]*rdx2R; 

  outL[6] += -1.0*incr4[0]*rdx2L; 
  outL[7] += incr4[1]*rdx2L-1.0*incrNonFlux4[1]*rdxSq4L; 

}
