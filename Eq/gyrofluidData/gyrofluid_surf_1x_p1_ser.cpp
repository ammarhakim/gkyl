#include <gyrofluid_mod_decl.h>

double gyrofluid_surf_1x_p1_ser_x(const double q_, const double m_, const double *wL1, const double *wR1, const double *dxL1, const double *dxR1, const double cMaxIn, const double *rJacL1, const double *rJacR1, const double *rBmagL1, const double *rBmagR1, const double *rBmagSqL1, const double *rBmagSqR1, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, const double *cRusL1, const double *cRusR1, double *outL, double *outR) 
{ 
  // q_,m_:              species charge and mass.
  // wL,wR:              cell-center in left and right cells.
  // dxL,dxR:            cell length in left and right cells.
  // cMaxIn:             maximum phase speed for numerical fluxes.
  // rJac:               reciprocal of jacobian (1/J).
  // rBmag:              reciprocal of magnetic field magnitude (1/B).
  // rBmagSq:            1/B^2.
  // sMomL,sMomR:        stepped moments (times Jacobian) in left and right cells.
  // phiL,phiR:          electrostatic potential in left and right cells.
  // primMomL,primMomR:  primitive moments (upar, Tpar, Tperp) in left and right cells.
  // cRus:               phase speed in Rusanov fluxes.
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

  double cMax = cMaxIn + fmax(fabs(1.224744871391589*primMomL1[1]+0.7071067811865475*primMomL1[0]),fabs(0.7071067811865475*primMomR1[0]-1.224744871391589*primMomR1[1]));
  double sMom1Favg[1];
  sMom1Favg[0] = -0.025*(((22.0454076850486*rBmagR1[1]-7.071067811865476*rBmagR1[0])*rJacR1[1]+rJacR1[0]*(12.24744871391589*rBmagR1[0]-7.071067811865476*rBmagR1[1]))*sMomR1[3]+(((-22.0454076850486*rBmagL1[1])-7.071067811865476*rBmagL1[0])*rJacL1[1]+rJacL1[0]*((-7.071067811865476*rBmagL1[1])-12.24744871391589*rBmagL1[0]))*sMomL1[3]+((12.24744871391589*rBmagR1[0]-7.071067811865476*rBmagR1[1])*rJacR1[1]+rJacR1[0]*(12.24744871391589*rBmagR1[1]-7.071067811865476*rBmagR1[0]))*sMomR1[2]+(((-7.071067811865476*rBmagL1[1])-12.24744871391589*rBmagL1[0])*rJacL1[1]+rJacL1[0]*((-12.24744871391589*rBmagL1[1])-7.071067811865476*rBmagL1[0]))*sMomL1[2]); 

  double momHat1[1];
  momHat1[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[1]+sMomL1[1])-1.0*sMomR1[0]+sMomL1[0])*cMax+2.828427124746191*sMom1Favg[0]); 

  double sMom2Favg[1];
  sMom2Favg[0] = -(0.05*((((22.0454076850486*rBmagR1[1]-7.071067811865476*rBmagR1[0])*rJacR1[1]+rJacR1[0]*(12.24744871391589*rBmagR1[0]-7.071067811865476*rBmagR1[1]))*sMomR1[5]+(((-22.0454076850486*rBmagL1[1])-7.071067811865476*rBmagL1[0])*rJacL1[1]+rJacL1[0]*((-7.071067811865476*rBmagL1[1])-12.24744871391589*rBmagL1[0]))*sMomL1[5]+((12.24744871391589*rBmagR1[0]-7.071067811865476*rBmagR1[1])*rJacR1[1]+rJacR1[0]*(12.24744871391589*rBmagR1[1]-7.071067811865476*rBmagR1[0]))*sMomR1[4]+(((-7.071067811865476*rBmagL1[1])-12.24744871391589*rBmagL1[0])*rJacL1[1]+rJacL1[0]*((-12.24744871391589*rBmagL1[1])-7.071067811865476*rBmagL1[0]))*sMomL1[4])*m_+(((5.0*rBmagR1[1]-8.660254037844386*rBmagR1[0])*rJacR1[1]+rJacR1[0]*(5.0*rBmagR1[0]-8.660254037844386*rBmagR1[1]))*sMomR1[1]+sMomR1[0]*((5.0*rBmagR1[0]-15.58845726811989*rBmagR1[1])*rJacR1[1]+rJacR1[0]*(5.0*rBmagR1[1]-8.660254037844386*rBmagR1[0])))*primMomR1[5]+(((5.0*rBmagL1[1]+8.660254037844386*rBmagL1[0])*rJacL1[1]+rJacL1[0]*(8.660254037844386*rBmagL1[1]+5.0*rBmagL1[0]))*sMomL1[1]+sMomL1[0]*((15.58845726811989*rBmagL1[1]+5.0*rBmagL1[0])*rJacL1[1]+rJacL1[0]*(5.0*rBmagL1[1]+8.660254037844386*rBmagL1[0])))*primMomL1[5]+(((5.0*rBmagR1[0]-15.58845726811989*rBmagR1[1])*rJacR1[1]+rJacR1[0]*(5.0*rBmagR1[1]-8.660254037844386*rBmagR1[0]))*sMomR1[1]+sMomR1[0]*((5.0*rBmagR1[1]-8.660254037844386*rBmagR1[0])*rJacR1[1]+rJacR1[0]*(5.0*rBmagR1[0]-8.660254037844386*rBmagR1[1])))*primMomR1[4]+(((15.58845726811989*rBmagL1[1]+5.0*rBmagL1[0])*rJacL1[1]+rJacL1[0]*(5.0*rBmagL1[1]+8.660254037844386*rBmagL1[0]))*sMomL1[1]+sMomL1[0]*((5.0*rBmagL1[1]+8.660254037844386*rBmagL1[0])*rJacL1[1]+rJacL1[0]*(8.660254037844386*rBmagL1[1]+5.0*rBmagL1[0])))*primMomL1[4]))/m_; 

  double momHat2[1];
  momHat2[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[3]+sMomL1[3])-1.0*sMomR1[2]+sMomL1[2])*cMax+2.828427124746191*sMom2Favg[0]); 

  double sMom3Favg[1];
  sMom3Favg[0] = (0.0125*(((((18.0*primMomR1[1]-31.17691453623978*primMomR1[0])*rBmagR1[1]+rBmagR1[0]*(10.0*primMomR1[0]-31.17691453623978*primMomR1[1]))*rJacR1[1]+rJacR1[0]*((10.0*primMomR1[0]-31.17691453623978*primMomR1[1])*rBmagR1[1]+rBmagR1[0]*(10.0*primMomR1[1]-17.32050807568877*primMomR1[0])))*sMomR1[5]+(((18.0*primMomL1[1]+31.17691453623978*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(31.17691453623978*primMomL1[1]+10.0*primMomL1[0]))*rJacL1[1]+rJacL1[0]*((31.17691453623978*primMomL1[1]+10.0*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(10.0*primMomL1[1]+17.32050807568877*primMomL1[0])))*sMomL1[5]+(((10.0*primMomR1[0]-31.17691453623978*primMomR1[1])*rBmagR1[1]+rBmagR1[0]*(10.0*primMomR1[1]-17.32050807568877*primMomR1[0]))*rJacR1[1]+rJacR1[0]*((10.0*primMomR1[1]-17.32050807568877*primMomR1[0])*rBmagR1[1]+rBmagR1[0]*(10.0*primMomR1[0]-17.32050807568877*primMomR1[1])))*sMomR1[4]+(((31.17691453623978*primMomL1[1]+10.0*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(10.0*primMomL1[1]+17.32050807568877*primMomL1[0]))*rJacL1[1]+rJacL1[0]*((10.0*primMomL1[1]+17.32050807568877*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(17.32050807568877*primMomL1[1]+10.0*primMomL1[0])))*sMomL1[4])*m_+((((7.071067811865476*primMomR1[0]-22.0454076850486*primMomR1[1])*rBmagR1[1]+rBmagR1[0]*(7.071067811865476*primMomR1[1]-12.24744871391589*primMomR1[0]))*rJacR1[1]+rJacR1[0]*((7.071067811865476*primMomR1[1]-12.24744871391589*primMomR1[0])*rBmagR1[1]+rBmagR1[0]*(7.071067811865476*primMomR1[0]-12.24744871391589*primMomR1[1])))*sMomR1[1]+sMomR1[0]*(((12.72792206135786*primMomR1[1]-22.0454076850486*primMomR1[0])*rBmagR1[1]+rBmagR1[0]*(7.071067811865476*primMomR1[0]-22.0454076850486*primMomR1[1]))*rJacR1[1]+rJacR1[0]*((7.071067811865476*primMomR1[0]-22.0454076850486*primMomR1[1])*rBmagR1[1]+rBmagR1[0]*(7.071067811865476*primMomR1[1]-12.24744871391589*primMomR1[0]))))*primMomR1[3]+((((22.0454076850486*primMomL1[1]+7.071067811865476*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(7.071067811865476*primMomL1[1]+12.24744871391589*primMomL1[0]))*rJacL1[1]+rJacL1[0]*((7.071067811865476*primMomL1[1]+12.24744871391589*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(12.24744871391589*primMomL1[1]+7.071067811865476*primMomL1[0])))*sMomL1[1]+sMomL1[0]*(((12.72792206135786*primMomL1[1]+22.0454076850486*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(22.0454076850486*primMomL1[1]+7.071067811865476*primMomL1[0]))*rJacL1[1]+rJacL1[0]*((22.0454076850486*primMomL1[1]+7.071067811865476*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(7.071067811865476*primMomL1[1]+12.24744871391589*primMomL1[0]))))*primMomL1[3]+((((12.72792206135786*primMomR1[1]-22.0454076850486*primMomR1[0])*rBmagR1[1]+rBmagR1[0]*(7.071067811865476*primMomR1[0]-22.0454076850486*primMomR1[1]))*rJacR1[1]+rJacR1[0]*((7.071067811865476*primMomR1[0]-22.0454076850486*primMomR1[1])*rBmagR1[1]+rBmagR1[0]*(7.071067811865476*primMomR1[1]-12.24744871391589*primMomR1[0])))*sMomR1[1]+sMomR1[0]*(((7.071067811865476*primMomR1[0]-22.0454076850486*primMomR1[1])*rBmagR1[1]+rBmagR1[0]*(7.071067811865476*primMomR1[1]-12.24744871391589*primMomR1[0]))*rJacR1[1]+rJacR1[0]*((7.071067811865476*primMomR1[1]-12.24744871391589*primMomR1[0])*rBmagR1[1]+rBmagR1[0]*(7.071067811865476*primMomR1[0]-12.24744871391589*primMomR1[1]))))*primMomR1[2]+((((12.72792206135786*primMomL1[1]+22.0454076850486*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(22.0454076850486*primMomL1[1]+7.071067811865476*primMomL1[0]))*rJacL1[1]+rJacL1[0]*((22.0454076850486*primMomL1[1]+7.071067811865476*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(7.071067811865476*primMomL1[1]+12.24744871391589*primMomL1[0])))*sMomL1[1]+sMomL1[0]*(((22.0454076850486*primMomL1[1]+7.071067811865476*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(7.071067811865476*primMomL1[1]+12.24744871391589*primMomL1[0]))*rJacL1[1]+rJacL1[0]*((7.071067811865476*primMomL1[1]+12.24744871391589*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(12.24744871391589*primMomL1[1]+7.071067811865476*primMomL1[0]))))*primMomL1[2]))/m_; 

  double momHat3[1];
  momHat3[0] = 0.3535533905932737*((1.732050807568877*(sMomR1[5]+sMomL1[5])-1.0*sMomR1[4]+sMomL1[4])*cMax+2.828427124746191*sMom3Favg[0]); 

  double sMom4Favg[1];
  sMom4Favg[0] = 0.025*((((9.0*primMomR1[1]-15.58845726811989*primMomR1[0])*rBmagR1[1]+rBmagR1[0]*(5.0*primMomR1[0]-15.58845726811989*primMomR1[1]))*rJacR1[1]+rJacR1[0]*((5.0*primMomR1[0]-15.58845726811989*primMomR1[1])*rBmagR1[1]+rBmagR1[0]*(5.0*primMomR1[1]-8.660254037844386*primMomR1[0])))*sMomR1[7]+(((9.0*primMomL1[1]+15.58845726811989*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(15.58845726811989*primMomL1[1]+5.0*primMomL1[0]))*rJacL1[1]+rJacL1[0]*((15.58845726811989*primMomL1[1]+5.0*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(5.0*primMomL1[1]+8.660254037844386*primMomL1[0])))*sMomL1[7]+(((5.0*primMomR1[0]-15.58845726811989*primMomR1[1])*rBmagR1[1]+rBmagR1[0]*(5.0*primMomR1[1]-8.660254037844386*primMomR1[0]))*rJacR1[1]+rJacR1[0]*((5.0*primMomR1[1]-8.660254037844386*primMomR1[0])*rBmagR1[1]+rBmagR1[0]*(5.0*primMomR1[0]-8.660254037844386*primMomR1[1])))*sMomR1[6]+(((15.58845726811989*primMomL1[1]+5.0*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(5.0*primMomL1[1]+8.660254037844386*primMomL1[0]))*rJacL1[1]+rJacL1[0]*((5.0*primMomL1[1]+8.660254037844386*primMomL1[0])*rBmagL1[1]+rBmagL1[0]*(8.660254037844386*primMomL1[1]+5.0*primMomL1[0])))*sMomL1[6]); 

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

  return 0.7071067811865475*cRusL1[0]; 
}
