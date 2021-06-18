#include <gyrofluid_mod_decl.h>

double gyrofluid_surf_1x_p1_ser_x(const double q_, const double m_, const double *wL1, const double *dxL1, const double *wR1, const double *dxR1, const double cMaxIn, const double *jacL, const double *rBmagL, const double *jacDbmagL, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, const double *csL1, const double *csR1, double *outL, double *outR) 
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
  uparL[0] = 0.5*(2.449489742783178*primMomL1[1]+1.414213562373095*primMomL1[0]); 

  double uparR[1]; 
  uparR[0] = -0.5*(2.449489742783178*primMomR1[1]-1.414213562373095*primMomR1[0]); 

  double cMax = cMaxIn + fabs(fmax(1.224744871391589*primMomL1[1]+0.7071067811865475*primMomL1[0],0.7071067811865475*primMomR1[0]-1.224744871391589*primMomR1[1]));
  double sMom1Favg[1];
  sMom1Favg[0] = -0.25*(2.449489742783178*sMomR1[3]-2.449489742783178*sMomL1[3]-1.414213562373095*(sMomR1[2]+sMomL1[2])); 

  double momHat1[1];
  momHat1[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[1]+sMomL1[1])-1.0*sMomR1[0]+sMomL1[0])*cMax+2.828427124746191*sMom1Favg[0]); 

  double sMom2Favg[1];
  sMom2Favg[0] = -(0.5*((2.449489742783178*sMomR1[5]-2.449489742783178*sMomL1[5]-1.414213562373095*(sMomR1[4]+sMomL1[4]))*m_+(sMomR1[1]-1.732050807568877*sMomR1[0])*primMomR1[5]+(sMomL1[1]+1.732050807568877*sMomL1[0])*primMomL1[5]+(sMomR1[0]-1.732050807568877*sMomR1[1])*primMomR1[4]+(1.732050807568877*sMomL1[1]+sMomL1[0])*primMomL1[4]))/m_; 

  double momHat2[1];
  momHat2[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[3]+sMomL1[3])-1.0*sMomR1[2]+sMomL1[2])*cMax+2.828427124746191*sMom2Favg[0]); 

  double sMom3Favg[1];
  sMom3Favg[0] = (0.125*(((2.0*primMomR1[1]-3.464101615137754*primMomR1[0])*sMomR1[5]+(2.0*primMomL1[1]+3.464101615137754*primMomL1[0])*sMomL1[5]+(2.0*primMomR1[0]-3.464101615137754*primMomR1[1])*sMomR1[4]+(3.464101615137754*primMomL1[1]+2.0*primMomL1[0])*sMomL1[4])*m_+((1.414213562373095*primMomR1[0]-2.449489742783178*primMomR1[1])*sMomR1[1]+sMomR1[0]*(1.414213562373095*primMomR1[1]-2.449489742783178*primMomR1[0]))*primMomR1[3]+((2.449489742783178*primMomL1[1]+1.414213562373095*primMomL1[0])*sMomL1[1]+sMomL1[0]*(1.414213562373095*primMomL1[1]+2.449489742783178*primMomL1[0]))*primMomL1[3]+((1.414213562373095*primMomR1[1]-2.449489742783178*primMomR1[0])*sMomR1[1]+sMomR1[0]*(1.414213562373095*primMomR1[0]-2.449489742783178*primMomR1[1]))*primMomR1[2]+((1.414213562373095*primMomL1[1]+2.449489742783178*primMomL1[0])*sMomL1[1]+sMomL1[0]*(2.449489742783178*primMomL1[1]+1.414213562373095*primMomL1[0]))*primMomL1[2]))/m_; 

  double momHat3[1];
  momHat3[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[5]+sMomL1[5])-1.0*sMomR1[4]+sMomL1[4])*cMax+2.828427124746191*sMom3Favg[0]); 

  double sMom4Favg[1];
  sMom4Favg[0] = 0.25*((primMomR1[1]-1.732050807568877*primMomR1[0])*sMomR1[7]+(primMomL1[1]+1.732050807568877*primMomL1[0])*sMomL1[7]+(primMomR1[0]-1.732050807568877*primMomR1[1])*sMomR1[6]+(1.732050807568877*primMomL1[1]+primMomL1[0])*sMomL1[6]); 

  double momHat4[1];
  momHat4[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[7]+sMomL1[7])-1.0*sMomR1[6]+sMomL1[6])*cMax+2.828427124746191*sMom4Favg[0]); 

  double incr1[2];
  incr1[0] = 0.7071067811865475*momHat1[0]; 
  incr1[1] = -1.224744871391589*momHat1[0]; 

  double incr2[2];
  incr2[0] = 0.7071067811865475*momHat2[0]; 
  incr2[1] = -1.224744871391589*momHat2[0]; 

  double incr3[2];
  incr3[0] = 0.7071067811865475*momHat3[0]; 
  incr3[1] = -1.224744871391589*momHat3[0]; 

  double incr4[2];
  incr4[0] = 0.7071067811865475*momHat4[0]; 
  incr4[1] = -1.224744871391589*momHat4[0]; 

  outR[0] += incr1[0]*rdx2R; 
  outR[1] += incr1[1]*rdx2R; 

  outL[0] += -1.0*incr1[0]*rdx2L; 
  outL[1] += incr1[1]*rdx2L; 

  outR[2] += incr2[0]*rdx2R; 
  outR[3] += incr2[1]*rdx2R; 

  outL[2] += -1.0*incr2[0]*rdx2L; 
  outL[3] += incr2[1]*rdx2L; 

  outR[4] += incr3[0]*rdx2R; 
  outR[5] += incr3[1]*rdx2R; 

  outL[4] += -1.0*incr3[0]*rdx2L; 
  outL[5] += incr3[1]*rdx2L; 

  outR[6] += incr4[0]*rdx2R; 
  outR[7] += incr4[1]*rdx2R; 

  outL[6] += -1.0*incr4[0]*rdx2L; 
  outL[7] += incr4[1]*rdx2L; 

  return 0.7071067811865475*csL1[0]; 
}
