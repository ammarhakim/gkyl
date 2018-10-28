#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x1vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[2]:    Cell-center coordinates. 
  // dxv[2]:  Cell spacing. 
  // nu:      diffusion coefficient (collisionality). 
  // u[2]:    bulk velocity. 
  // vtSq[2]: thermal speed squared. 
  // f[4]:    Input distribution function. 
  // out[4]:  Incremented output 
  double rdv2nu[1]; 
  double rdvSq4nu[1]; 
  rdv2nu[0] = 2.0/dxv[1]; 
  rdvSq4nu[0] = nu*rdv2nu[0]*rdv2nu[0]; 
  rdv2nu[0] = nu*rdv2nu[0]; 

  double alpha0[4]; 
  alpha0[0] = rdv2nu[0]*(1.414213562373095*u[0]-2.0*w[1]); 
  alpha0[1] = 1.414213562373095*rdv2nu[0]*u[1]; 
  alpha0[2] = -1.154700538379252*nu; 

  double alpha1[4]; 
  alpha1[0] = 1.414213562373095*rdvSq4nu[0]*vtSq[0]; 
  alpha1[1] = 1.414213562373095*rdvSq4nu[0]*vtSq[1]; 

  out[2] += 0.8660254037844386*(alpha0[2]*f[2]+alpha0[1]*f[1]+alpha0[0]*f[0]); 
  out[3] += 0.8660254037844386*(alpha0[2]*f[3]+alpha0[0]*f[1]+f[0]*alpha0[1]); 

  const double alpha0Mid = 0.5*alpha0[0]; 
  const double alpha1Mid = 0.6666666666666666*alpha1[0]; 
  return std::abs(alpha0Mid) + alpha1Mid; 

} 
double GkLBOconstNuVol1x1vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[2]:    Cell-center coordinates. 
  // dxv[2]:  Cell spacing. 
  // nu:      diffusion coefficient (collisionality). 
  // u[3]:    bulk velocity. 
  // vtSq[3]: thermal speed squared. 
  // f[8]:    Input distribution function. 
  // out[8]:  Incremented output 
  double rdv2nu[1]; 
  double rdvSq4nu[1]; 
  rdv2nu[0] = 2.0/dxv[1]; 
  rdvSq4nu[0] = nu*rdv2nu[0]*rdv2nu[0]; 
  rdv2nu[0] = nu*rdv2nu[0]; 

  double alpha0[8]; 
  alpha0[0] = rdv2nu[0]*(1.414213562373095*u[0]-2.0*w[1]); 
  alpha0[1] = 1.414213562373095*rdv2nu[0]*u[1]; 
  alpha0[2] = -1.154700538379252*nu; 
  alpha0[4] = 1.414213562373095*rdv2nu[0]*u[2]; 

  double alpha1[8]; 
  alpha1[0] = 1.414213562373095*rdvSq4nu[0]*vtSq[0]; 
  alpha1[1] = 1.414213562373095*rdvSq4nu[0]*vtSq[1]; 
  alpha1[4] = 1.414213562373095*rdvSq4nu[0]*vtSq[2]; 

  out[2] += 0.8660254037844386*(alpha0[4]*f[4]+alpha0[2]*f[2]+alpha0[1]*f[1]+alpha0[0]*f[0]); 
  out[3] += 0.7745966692414833*(alpha0[1]*f[4]+f[1]*alpha0[4])+0.8660254037844386*(alpha0[2]*f[3]+alpha0[0]*f[1]+f[0]*alpha0[1]); 
  out[5] += 1.936491673103709*alpha0[4]*f[6]+1.732050807568877*alpha0[2]*f[5]+3.354101966249685*alpha1[4]*f[4]+1.936491673103709*(alpha0[1]*f[3]+alpha0[0]*f[2]+f[0]*alpha0[2])+3.354101966249685*(alpha1[1]*f[1]+alpha1[0]*f[0]); 
  out[6] += 0.8660254037844386*alpha0[2]*f[6]+0.5532833351724881*alpha0[4]*f[4]+0.8660254037844386*(alpha0[0]*f[4]+f[0]*alpha0[4])+0.7745966692414833*alpha0[1]*f[1]; 
  out[7] += 1.732050807568877*(alpha0[2]*f[7]+alpha0[1]*f[6])+3.0*(alpha1[1]*f[4]+f[1]*alpha1[4])+1.732050807568877*f[3]*alpha0[4]+1.936491673103709*(alpha0[0]*f[3]+alpha0[1]*f[2]+f[1]*alpha0[2])+3.354101966249685*(alpha1[0]*f[1]+f[0]*alpha1[1]); 

  const double alpha0Mid = 0.5*alpha0[0]-0.5590169943749475*alpha0[4]; 
  const double alpha1Mid = (9*(0.5*alpha1[0]-0.5590169943749475*alpha1[4]))/5; 
  return std::abs(alpha0Mid) + alpha1Mid; 

} 
double GkLBOconstNuVol1x1vSerP3(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[2]:    Cell-center coordinates. 
  // dxv[2]:  Cell spacing. 
  // nu:      diffusion coefficient (collisionality). 
  // u[4]:    bulk velocity. 
  // vtSq[4]: thermal speed squared. 
  // f[12]:    Input distribution function. 
  // out[12]:  Incremented output 
  double rdv2nu[1]; 
  double rdvSq4nu[1]; 
  rdv2nu[0] = 2.0/dxv[1]; 
  rdvSq4nu[0] = nu*rdv2nu[0]*rdv2nu[0]; 
  rdv2nu[0] = nu*rdv2nu[0]; 

  double alpha0[12]; 
  alpha0[0] = rdv2nu[0]*(1.414213562373095*u[0]-2.0*w[1]); 
  alpha0[1] = 1.414213562373095*rdv2nu[0]*u[1]; 
  alpha0[2] = -1.154700538379252*nu; 
  alpha0[4] = 1.414213562373095*rdv2nu[0]*u[2]; 
  alpha0[8] = 1.414213562373095*rdv2nu[0]*u[3]; 

  double alpha1[12]; 
  alpha1[0] = 1.414213562373095*rdvSq4nu[0]*vtSq[0]; 
  alpha1[1] = 1.414213562373095*rdvSq4nu[0]*vtSq[1]; 
  alpha1[4] = 1.414213562373095*rdvSq4nu[0]*vtSq[2]; 
  alpha1[8] = 1.414213562373095*rdvSq4nu[0]*vtSq[3]; 

  out[2] += 0.8660254037844386*(alpha0[8]*f[8]+alpha0[4]*f[4]+alpha0[2]*f[2]+alpha0[1]*f[1]+alpha0[0]*f[0]); 
  out[3] += 0.7606388292556648*(alpha0[4]*f[8]+f[4]*alpha0[8])+0.7745966692414833*(alpha0[1]*f[4]+f[1]*alpha0[4])+0.8660254037844386*(alpha0[2]*f[3]+alpha0[0]*f[1]+f[0]*alpha0[1]); 
  out[5] += 1.936491673103709*alpha0[8]*f[10]+3.354101966249685*alpha1[8]*f[8]+1.936491673103709*alpha0[4]*f[6]+1.732050807568877*alpha0[2]*f[5]+3.354101966249685*alpha1[4]*f[4]+1.936491673103709*(alpha0[1]*f[3]+alpha0[0]*f[2]+f[0]*alpha0[2])+3.354101966249685*(alpha1[1]*f[1]+alpha1[0]*f[0]); 
  out[6] += 0.5163977794943223*alpha0[8]*f[8]+0.7606388292556648*(alpha0[1]*f[8]+f[1]*alpha0[8])+0.8660254037844386*alpha0[2]*f[6]+0.5532833351724881*alpha0[4]*f[4]+0.8660254037844386*(alpha0[0]*f[4]+f[0]*alpha0[4])+0.7745966692414833*alpha0[1]*f[1]; 
  out[7] += 1.700840128541522*alpha0[4]*f[10]+2.945941518185896*(alpha1[4]*f[8]+f[4]*alpha1[8])+1.700840128541522*f[6]*alpha0[8]+1.732050807568877*(alpha0[2]*f[7]+alpha0[1]*f[6])+3.0*(alpha1[1]*f[4]+f[1]*alpha1[4])+1.732050807568877*f[3]*alpha0[4]+1.936491673103709*(alpha0[0]*f[3]+alpha0[1]*f[2]+f[1]*alpha0[2])+3.354101966249685*(alpha1[0]*f[1]+f[0]*alpha1[1]); 
  out[9] += 11.4564392373896*alpha1[8]*f[10]+2.598076211353316*alpha0[2]*f[9]+1.322875655532295*alpha0[8]*f[8]+2.958039891549809*alpha0[1]*f[7]+11.4564392373896*alpha1[4]*f[6]+2.958039891549809*alpha0[0]*f[5]+1.322875655532295*alpha0[4]*f[4]+11.4564392373896*alpha1[1]*f[3]+(3.968626966596886*alpha0[2]+11.4564392373896*alpha1[0])*f[2]+1.322875655532295*(alpha0[1]*f[1]+alpha0[0]*f[0]); 
  out[10] += 0.8660254037844386*alpha0[2]*f[10]+(0.5163977794943223*alpha0[4]+0.8660254037844386*alpha0[0])*f[8]+(0.5163977794943223*f[4]+0.8660254037844386*f[0])*alpha0[8]+0.7606388292556648*(alpha0[1]*f[4]+f[1]*alpha0[4]); 
  out[11] += 2.598076211353316*alpha0[2]*f[11]+10.06230589874905*alpha1[4]*f[10]+1.161895003862225*alpha0[4]*f[8]+10.06230589874905*f[6]*alpha1[8]+1.161895003862225*f[4]*alpha0[8]+(2.645751311064591*alpha0[4]+2.958039891549809*alpha0[0])*f[7]+10.2469507659596*alpha1[1]*f[6]+alpha0[1]*(2.958039891549809*f[5]+1.183215956619923*f[4])+10.2469507659596*f[3]*alpha1[4]+1.183215956619923*f[1]*alpha0[4]+3.968626966596886*alpha0[2]*f[3]+11.4564392373896*(alpha1[0]*f[3]+alpha1[1]*f[2])+1.322875655532295*(alpha0[0]*f[1]+f[0]*alpha0[1]); 

  const double alpha0Mid = 0.5*alpha0[0]-0.5590169943749475*alpha0[4]; 
  const double alpha1Mid = (16*(0.5*alpha1[0]-0.5590169943749475*alpha1[4]))/7; 
  return std::abs(alpha0Mid) + alpha1Mid; 

} 
