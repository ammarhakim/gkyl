#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x1vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[1]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[1]*dxv[1]); 

  double alpha_mid = 0.0; 
  double alpha_drag[3]; 
  double alpha_diffusion[3]; 

  // Expand rdv2nu*(vx-ux) in phase basis.
  alpha_drag[0] = (1.414213562373095*u[0]-2.0*w[1])*rdvx2nu; 
  alpha_drag[1] = 1.414213562373095*u[1]*rdvx2nu; 
  alpha_drag[2] = -0.5773502691896258*dxv[1]*rdvx2nu; 

  // x-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.25*alpha_drag[0]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 1.414213562373095*vtSq[0]; 
  alpha_diffusion[1] = 1.414213562373095*vtSq[1]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.6666666666666666*alpha_diffusion[0]*rdvxSq4nu); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol1x1vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[1]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[1]*dxv[1]); 

  double alpha_mid = 0.0; 
  double alpha_drag[6]; 
  double alpha_diffusion[6]; 

  // Expand rdv2nu*(vx-ux) in phase basis.
  alpha_drag[0] = (1.414213562373095*u[0]-2.0*w[1])*rdvx2nu; 
  alpha_drag[1] = 1.414213562373095*u[1]*rdvx2nu; 
  alpha_drag[2] = -0.5773502691896258*dxv[1]*rdvx2nu; 
  alpha_drag[4] = 1.414213562373095*u[2]*rdvx2nu; 

  // x-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.25*alpha_drag[0]-0.2795084971874737*alpha_drag[4]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 1.414213562373095*vtSq[0]; 
  alpha_diffusion[1] = 1.414213562373095*vtSq[1]; 
  alpha_diffusion[4] = 1.414213562373095*vtSq[2]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((0.9*alpha_diffusion[0]-1.006230589874905*alpha_diffusion[4])*rdvxSq4nu); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alpha_drag[4]*f[4]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.7745966692414833*(alpha_drag[1]*f[4]+f[1]*alpha_drag[4])+0.8660254037844386*(alpha_drag[2]*f[3]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[5] += 3.354101966249685*(alpha_diffusion[4]*f[4]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4nu+1.732050807568877*alpha_drag[2]*f[5]+1.936491673103709*(alpha_drag[1]*f[3]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol1x1vMaxP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[1]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[1]*dxv[1]); 

  double alpha_mid = 0.0; 
  double alpha_drag[10]; 
  double alpha_diffusion[10]; 

  // Expand rdv2nu*(vx-ux) in phase basis.
  alpha_drag[0] = (1.414213562373095*u[0]-2.0*w[1])*rdvx2nu; 
  alpha_drag[1] = 1.414213562373095*u[1]*rdvx2nu; 
  alpha_drag[2] = -0.5773502691896258*dxv[1]*rdvx2nu; 
  alpha_drag[4] = 1.414213562373095*u[2]*rdvx2nu; 
  alpha_drag[8] = 1.414213562373095*u[3]*rdvx2nu; 

  // x-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.25*alpha_drag[0]-0.2795084971874737*alpha_drag[4]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 1.414213562373095*vtSq[0]; 
  alpha_diffusion[1] = 1.414213562373095*vtSq[1]; 
  alpha_diffusion[4] = 1.414213562373095*vtSq[2]; 
  alpha_diffusion[8] = 1.414213562373095*vtSq[3]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((1.142857142857143*alpha_diffusion[0]-1.27775312999988*alpha_diffusion[4])*rdvxSq4nu); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alpha_drag[8]*f[8]+alpha_drag[4]*f[4]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.7606388292556648*(alpha_drag[4]*f[8]+f[4]*alpha_drag[8])+0.7745966692414833*(alpha_drag[1]*f[4]+f[1]*alpha_drag[4])+0.8660254037844386*(alpha_drag[2]*f[3]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[5] += 3.354101966249685*(alpha_diffusion[8]*f[8]+alpha_diffusion[4]*f[4]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4nu+1.936491673103709*alpha_drag[4]*f[6]+1.732050807568877*alpha_drag[2]*f[5]+1.936491673103709*(alpha_drag[1]*f[3]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[6] += 0.5163977794943223*alpha_drag[8]*f[8]+0.7606388292556648*(alpha_drag[1]*f[8]+f[1]*alpha_drag[8])+0.8660254037844386*alpha_drag[2]*f[6]+0.5532833351724881*alpha_drag[4]*f[4]+0.8660254037844386*(alpha_drag[0]*f[4]+f[0]*alpha_drag[4])+0.7745966692414833*alpha_drag[1]*f[1]; 
  out[7] += (2.945941518185896*(alpha_diffusion[4]*f[8]+f[4]*alpha_diffusion[8])+3.0*(alpha_diffusion[1]*f[4]+f[1]*alpha_diffusion[4])+3.354101966249685*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvxSq4nu+1.700840128541522*f[6]*alpha_drag[8]+1.732050807568877*(alpha_drag[2]*f[7]+alpha_drag[1]*f[6]+f[3]*alpha_drag[4])+1.936491673103709*(alpha_drag[0]*f[3]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 
  out[9] += 11.4564392373896*(alpha_diffusion[4]*f[6]+alpha_diffusion[1]*f[3]+alpha_diffusion[0]*f[2])*rdvxSq4nu+2.598076211353316*alpha_drag[2]*f[9]+1.322875655532295*alpha_drag[8]*f[8]+2.958039891549809*(alpha_drag[1]*f[7]+alpha_drag[0]*f[5])+1.322875655532295*alpha_drag[4]*f[4]+3.968626966596886*alpha_drag[2]*f[2]+1.322875655532295*(alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 

  return alpha_mid; 

} 
