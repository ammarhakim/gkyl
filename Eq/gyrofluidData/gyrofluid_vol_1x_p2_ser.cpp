#include <gyrofluid_mod_decl.h>

double gyrofluid_vol_1x_p2_ser(const double q_, const double m_, const double *w, const double *dx, const double *jac, const double *rBmag, const double *jacDbmag, const double *dBdz, const double *sMom, const double *phi, const double *primMom, const double *cs, double *out) 
{ 
  // q_,m_:   species charge and mass.
  // w:       cell-center.
  // dx:      cell length.
  // uMaxIn:  maximum speed.
  // jac:     jacobian.
  // rBmag:   reciprocal of magnetic field magnitude (1/B).
  // jacDbmag:jacobian divided by B (J/B).
  // sMom:    stepped moments (times Jacobian).
  // phi:     electrostatic potential.
  // primMom: primitive moments (upar, Tpar, Tperp).
  // cs:      sound speed.
  // out:     output increment.

  double wx = w[0];
  double rdx2 = 2.0/dx[0];
  double rdxSq4 = rdx2*rdx2;

  out[1] += 1.732050807568877*sMom[3]*rdx2; 
  out[2] += 3.872983346207417*sMom[4]*rdx2; 

  out[3] += (((-2.738612787525831*sMom[1]*phi[2])-1.224744871391589*sMom[0]*phi[1])*q_*rdx2)/m_-0.7071067811865475*(dBdz[2]*sMom[11]+dBdz[1]*sMom[10]+dBdz[0]*sMom[9]); 
  out[4] += (((phi[2]*((-2.449489742783178*sMom[2])-2.738612787525831*sMom[0])-1.224744871391589*phi[1]*sMom[1])*q_-2.449489742783178*(sMom[2]*primMom[8]+sMom[1]*primMom[7]+sMom[0]*primMom[6]))/m_+3.464101615137754*sMom[6])*rdx2-0.6324555320336759*(dBdz[1]*sMom[11]+dBdz[2]*sMom[10])-0.7071067811865475*(dBdz[0]*sMom[10]+dBdz[1]*sMom[9]); 
  out[5] += ((sMom[1]*((-2.449489742783178*phi[2]*q_)-4.898979485566357*primMom[8])-1.224744871391589*phi[1]*sMom[2]*q_-5.477225575051662*(sMom[0]*primMom[7]+sMom[1]*primMom[6])-4.898979485566357*sMom[2]*primMom[7])/m_+7.745966692414834*sMom[7])*rdx2+((-0.4517539514526256*dBdz[2])-0.7071067811865475*dBdz[0])*sMom[11]-0.6324555320336759*dBdz[1]*sMom[10]-0.7071067811865475*dBdz[2]*sMom[9]; 

  out[6] += (((-2.738612787525831*phi[2]*sMom[4])-1.224744871391589*phi[1]*sMom[3])*q_*rdx2)/m_; 
  out[7] += ((((-2.449489742783178*phi[2]*sMom[5])-1.224744871391589*phi[1]*sMom[4]-2.738612787525831*phi[2]*sMom[3])*q_+(0.8660254037844386*(primMom[0]*sMom[2]+sMom[0]*primMom[2])+0.5532833351724881*primMom[2]*sMom[2]+0.7745966692414833*primMom[1]*sMom[1])*primMom[5]+0.8660254037844386*((primMom[0]*sMom[1]+sMom[0]*primMom[1])*primMom[4]+(primMom[2]*sMom[2]+primMom[1]*sMom[1]+primMom[0]*sMom[0])*primMom[3])+0.7745966692414833*(primMom[1]*sMom[2]+sMom[1]*primMom[2])*primMom[4])/m_+1.224744871391589*(primMom[2]*sMom[8]+primMom[1]*sMom[7]+primMom[0]*sMom[6]))*rdx2; 
  out[8] += ((((-1.224744871391589*phi[1]*sMom[5])-2.449489742783178*phi[2]*sMom[4])*q_+(3.043058343448684*primMom[1]*sMom[2]+1.549193338482967*sMom[1]*primMom[2]+1.732050807568877*(primMom[0]*sMom[1]+sMom[0]*primMom[1]))*primMom[5]+(1.732050807568877*(primMom[0]*sMom[2]+sMom[0]*primMom[2])+1.549193338482967*primMom[2]*sMom[2]+3.485685011586674*primMom[1]*sMom[1]+1.936491673103709*primMom[0]*sMom[0])*primMom[4]+(1.732050807568877*(primMom[1]*sMom[2]+sMom[1]*primMom[2])+1.936491673103709*(primMom[0]*sMom[1]+sMom[0]*primMom[1]))*primMom[3])/m_+2.449489742783178*(primMom[1]*sMom[8]+primMom[2]*sMom[7])+2.738612787525831*(primMom[0]*sMom[7]+primMom[1]*sMom[6]))*rdx2; 

  out[10] += 1.224744871391589*(primMom[2]*sMom[11]+primMom[1]*sMom[10]+primMom[0]*sMom[9])*rdx2; 
  out[11] += (2.449489742783178*(primMom[1]*sMom[11]+primMom[2]*sMom[10])+2.738612787525831*(primMom[0]*sMom[10]+primMom[1]*sMom[9]))*rdx2; 

  return rdx2*(fabs(0.7071067811865475*primMom[0]-0.7905694150420947*primMom[2]) + 0.7071067811865475*cs[0]-0.7905694150420947*cs[2]); 
}
