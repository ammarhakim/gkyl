#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x1vSerP1(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdvx2 = 2/dxv[1]; 
  const double rdvxSq4 = 4/(dxv[1]*dxv[1]); 

  double alpha_mid = 0.0; 
  double alpha_drag[4]; 
  double alpha_diffusion[4]; 

  // Expand rdv2*nu*(vx-ux) in phase basis.
  alpha_drag[0] = (1.414213562373095*nuU[0]-2.0*w[1]*nu)*rdvx2; 
  alpha_drag[1] = 1.414213562373095*nuU[1]*rdvx2; 
  alpha_drag[2] = -0.5773502691896258*dxv[1]*nu*rdvx2; 

  // cvard[1]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.25*alpha_drag[0]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 1.414213562373095*nuVtSq[0]; 
  alpha_diffusion[1] = 1.414213562373095*nuVtSq[1]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.6666666666666666*alpha_diffusion[0]*rdvxSq4); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.8660254037844386*(alpha_drag[2]*f[3]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol1x1vSerP2(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdvx2 = 2/dxv[1]; 
  const double rdvxSq4 = 4/(dxv[1]*dxv[1]); 

  double alpha_mid = 0.0; 
  double alpha_drag[8]; 
  double alpha_diffusion[8]; 

  // Expand rdv2*nu*(vx-ux) in phase basis.
  alpha_drag[0] = (1.414213562373095*nuU[0]-2.0*w[1]*nu)*rdvx2; 
  alpha_drag[1] = 1.414213562373095*nuU[1]*rdvx2; 
  alpha_drag[2] = -0.5773502691896258*dxv[1]*nu*rdvx2; 
  alpha_drag[4] = 1.414213562373095*nuU[2]*rdvx2; 

  // cvard[1]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.25*alpha_drag[0]-0.2795084971874737*alpha_drag[4]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 1.414213562373095*nuVtSq[0]; 
  alpha_diffusion[1] = 1.414213562373095*nuVtSq[1]; 
  alpha_diffusion[4] = 1.414213562373095*nuVtSq[2]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((0.9*alpha_diffusion[0]-1.006230589874905*alpha_diffusion[4])*rdvxSq4); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alpha_drag[4]*f[4]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.7745966692414833*(alpha_drag[1]*f[4]+f[1]*alpha_drag[4])+0.8660254037844386*(alpha_drag[2]*f[3]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[5] += 3.354101966249685*(alpha_diffusion[4]*f[4]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4+1.936491673103709*alpha_drag[4]*f[6]+1.732050807568877*alpha_drag[2]*f[5]+1.936491673103709*(alpha_drag[1]*f[3]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[6] += 0.8660254037844386*alpha_drag[2]*f[6]+0.5532833351724881*alpha_drag[4]*f[4]+0.8660254037844386*(alpha_drag[0]*f[4]+f[0]*alpha_drag[4])+0.7745966692414833*alpha_drag[1]*f[1]; 
  out[7] += (3.0*(alpha_diffusion[1]*f[4]+f[1]*alpha_diffusion[4])+3.354101966249685*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvxSq4+1.732050807568877*(alpha_drag[2]*f[7]+alpha_drag[1]*f[6]+f[3]*alpha_drag[4])+1.936491673103709*(alpha_drag[0]*f[3]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol1x1vSerP3(const double *w, const double *dxv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  const double rdvx2 = 2/dxv[1]; 
  const double rdvxSq4 = 4/(dxv[1]*dxv[1]); 

  double alpha_mid = 0.0; 
  double alpha_drag[12]; 
  double alpha_diffusion[12]; 

  // Expand rdv2*nu*(vx-ux) in phase basis.
  alpha_drag[0] = (1.414213562373095*nuU[0]-2.0*w[1]*nu)*rdvx2; 
  alpha_drag[1] = 1.414213562373095*nuU[1]*rdvx2; 
  alpha_drag[2] = -0.5773502691896258*dxv[1]*nu*rdvx2; 
  alpha_drag[4] = 1.414213562373095*nuU[2]*rdvx2; 
  alpha_drag[8] = 1.414213562373095*nuU[3]*rdvx2; 

  // cvard[1]-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.25*alpha_drag[0]-0.2795084971874737*alpha_drag[4]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 1.414213562373095*nuVtSq[0]; 
  alpha_diffusion[1] = 1.414213562373095*nuVtSq[1]; 
  alpha_diffusion[4] = 1.414213562373095*nuVtSq[2]; 
  alpha_diffusion[8] = 1.414213562373095*nuVtSq[3]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((1.142857142857143*alpha_diffusion[0]-1.27775312999988*alpha_diffusion[4])*rdvxSq4); 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alpha_drag[8]*f[8]+alpha_drag[4]*f[4]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[3] += 0.7606388292556648*(alpha_drag[4]*f[8]+f[4]*alpha_drag[8])+0.7745966692414833*(alpha_drag[1]*f[4]+f[1]*alpha_drag[4])+0.8660254037844386*(alpha_drag[2]*f[3]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[5] += 3.354101966249685*(alpha_diffusion[8]*f[8]+alpha_diffusion[4]*f[4]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4+1.936491673103709*(alpha_drag[8]*f[10]+alpha_drag[4]*f[6])+1.732050807568877*alpha_drag[2]*f[5]+1.936491673103709*(alpha_drag[1]*f[3]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[6] += 0.5163977794943223*alpha_drag[8]*f[8]+0.7606388292556648*(alpha_drag[1]*f[8]+f[1]*alpha_drag[8])+0.8660254037844386*alpha_drag[2]*f[6]+0.5532833351724881*alpha_drag[4]*f[4]+0.8660254037844386*(alpha_drag[0]*f[4]+f[0]*alpha_drag[4])+0.7745966692414833*alpha_drag[1]*f[1]; 
  out[7] += (2.945941518185896*(alpha_diffusion[4]*f[8]+f[4]*alpha_diffusion[8])+3.0*(alpha_diffusion[1]*f[4]+f[1]*alpha_diffusion[4])+3.354101966249685*(alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvxSq4+1.700840128541522*(alpha_drag[4]*f[10]+f[6]*alpha_drag[8])+1.732050807568877*(alpha_drag[2]*f[7]+alpha_drag[1]*f[6]+f[3]*alpha_drag[4])+1.936491673103709*(alpha_drag[0]*f[3]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 
  out[9] += 11.4564392373896*(alpha_diffusion[8]*f[10]+alpha_diffusion[4]*f[6]+alpha_diffusion[1]*f[3]+alpha_diffusion[0]*f[2])*rdvxSq4+2.598076211353316*alpha_drag[2]*f[9]+1.322875655532295*alpha_drag[8]*f[8]+2.958039891549809*(alpha_drag[1]*f[7]+alpha_drag[0]*f[5])+1.322875655532295*alpha_drag[4]*f[4]+3.968626966596886*alpha_drag[2]*f[2]+1.322875655532295*(alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[10] += 0.8660254037844386*alpha_drag[2]*f[10]+(0.5163977794943223*alpha_drag[4]+0.8660254037844386*alpha_drag[0])*f[8]+(0.5163977794943223*f[4]+0.8660254037844386*f[0])*alpha_drag[8]+0.7606388292556648*(alpha_drag[1]*f[4]+f[1]*alpha_drag[4]); 
  out[11] += (10.06230589874905*(alpha_diffusion[4]*f[10]+f[6]*alpha_diffusion[8])+10.2469507659596*(alpha_diffusion[1]*f[6]+f[3]*alpha_diffusion[4])+11.4564392373896*(alpha_diffusion[0]*f[3]+alpha_diffusion[1]*f[2]))*rdvxSq4+2.598076211353316*alpha_drag[2]*f[11]+1.161895003862225*(alpha_drag[4]*f[8]+f[4]*alpha_drag[8])+2.645751311064591*alpha_drag[4]*f[7]+2.958039891549809*(alpha_drag[0]*f[7]+alpha_drag[1]*f[5])+1.183215956619923*(alpha_drag[1]*f[4]+f[1]*alpha_drag[4])+3.968626966596886*alpha_drag[2]*f[3]+1.322875655532295*(alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 

  return alpha_mid; 

} 
