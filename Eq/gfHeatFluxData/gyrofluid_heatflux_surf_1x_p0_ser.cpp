#include <gyrofluid_heatflux_mod_decl.h>

void gyrofluid_heatflux_surf_1x_p0_ser_x(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *wR1, const double *dxL1, const double *dxR1, const double cMaxIn, const double *rBmagL1, const double *rBmagR1, const double *rBmagSqL1, const double *rBmagSqR1, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR) 
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

  double GheatF3[1];

  GheatF3[0] = -0.125*((1.414213562373095*phiR1[0]-1.414213562373095*phiL1[0])*rBmagL1[0]*kappaPerp*kperpSq*q_+(2.828427124746191*rBmagL1[0]*primMomL1[2]-2.828427124746191*rBmagL1[0]*primMomR1[2])*kappaPerp+(1.414213562373095*rBmagL1[0]*primMomL1[1]-1.414213562373095*rBmagL1[0]*primMomR1[1])*kappaPar)*rdx2L; 

  double GheatF4[1];

  GheatF4[0] = -0.125*((1.414213562373095*phiR1[0]-1.414213562373095*phiL1[0])*rBmagSqL1[0]*kappaPerp*kperpSq*q_+(2.828427124746191*rBmagSqL1[0]*primMomL1[2]-2.828427124746191*rBmagSqL1[0]*primMomR1[2])*kappaPerp)*rdx2L; 

  double incr3[1];
  incr3[0] = -0.5*GheatF3[0]; 

  double incrNonFlux3[1];

  double incr4[1];
  incr4[0] = -0.5*GheatF4[0]; 

  double incrNonFlux4[1];

  outR[2] += incr3[0]*rdx2R; 

  outL[2] += -1.0*incr3[0]*rdx2L; 

  outR[3] += incr4[0]*rdx2R; 

  outL[3] += -1.0*incr4[0]*rdx2L; 

}
