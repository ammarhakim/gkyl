#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x1vSerP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
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

  double alphaDrag[4]; 
  // Expand rdv2nu*(vx-ux) in phase basis.
  alphaDrag[0] = (1.414213562373095*u[0]-2.0*w[1])*rdvx2nu; 
  alphaDrag[1] = 1.414213562373095*u[1]*rdvx2nu; 
  alphaDrag[2] = -0.5773502691896258*dxv[1]*rdvx2nu; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaDrag[2]*f[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 

  return std::abs(0.25*alphaDrag[0])+std::abs(0.9428090415820636*vtSq[0]*rdvxSq4nu); 

} 
double VmLBOconstNuVol1x1vSerP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
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

  double alphaDrag[8]; 
  // Expand rdv2nu*(vx-ux) in phase basis.
  alphaDrag[0] = (1.414213562373095*u[0]-2.0*w[1])*rdvx2nu; 
  alphaDrag[1] = 1.414213562373095*u[1]*rdvx2nu; 
  alphaDrag[2] = -0.5773502691896258*dxv[1]*rdvx2nu; 
  alphaDrag[4] = 1.414213562373095*u[2]*rdvx2nu; 

  double facDiff[3]; 
  // Expand nu*vthSq in phase basis.
  facDiff[0] = vtSq[0]; 
  facDiff[1] = vtSq[1]; 
  facDiff[2] = vtSq[2]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alphaDrag[4]*f[4]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.7745966692414833*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+0.8660254037844386*(alphaDrag[2]*f[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 4.743416490252569*(facDiff[2]*f[4]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4nu+1.936491673103709*alphaDrag[4]*f[6]+1.732050807568877*alphaDrag[2]*f[5]+1.936491673103709*(alphaDrag[1]*f[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[6] += 0.8660254037844386*alphaDrag[2]*f[6]+0.5532833351724881*alphaDrag[4]*f[4]+0.8660254037844386*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+0.7745966692414833*alphaDrag[1]*f[1]; 
  out[7] += (4.242640687119286*(facDiff[1]*f[4]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4nu+1.732050807568877*(alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+f[3]*alphaDrag[4])+1.936491673103709*(alphaDrag[0]*f[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 

  return std::abs(0.25*alphaDrag[0]-0.2795084971874737*alphaDrag[4])+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvxSq4nu); 

} 
double VmLBOconstNuVol1x1vSerP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
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

  double alphaDrag[12]; 
  // Expand rdv2nu*(vx-ux) in phase basis.
  alphaDrag[0] = (1.414213562373095*u[0]-2.0*w[1])*rdvx2nu; 
  alphaDrag[1] = 1.414213562373095*u[1]*rdvx2nu; 
  alphaDrag[2] = -0.5773502691896258*dxv[1]*rdvx2nu; 
  alphaDrag[4] = 1.414213562373095*u[2]*rdvx2nu; 
  alphaDrag[8] = 1.414213562373095*u[3]*rdvx2nu; 

  double facDiff[4]; 
  // Expand nu*vthSq in phase basis.
  facDiff[0] = vtSq[0]; 
  facDiff[1] = vtSq[1]; 
  facDiff[2] = vtSq[2]; 
  facDiff[3] = vtSq[3]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alphaDrag[8]*f[8]+alphaDrag[4]*f[4]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.7606388292556648*(alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+0.7745966692414833*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+0.8660254037844386*(alphaDrag[2]*f[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 4.743416490252569*(facDiff[3]*f[8]+facDiff[2]*f[4]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4nu+1.936491673103709*(alphaDrag[8]*f[10]+alphaDrag[4]*f[6])+1.732050807568877*alphaDrag[2]*f[5]+1.936491673103709*(alphaDrag[1]*f[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[6] += 0.5163977794943223*alphaDrag[8]*f[8]+0.7606388292556648*(alphaDrag[1]*f[8]+f[1]*alphaDrag[8])+0.8660254037844386*alphaDrag[2]*f[6]+0.5532833351724881*alphaDrag[4]*f[4]+0.8660254037844386*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+0.7745966692414833*alphaDrag[1]*f[1]; 
  out[7] += (4.166190448976479*facDiff[2]*f[8]+4.242640687119286*(facDiff[1]*f[4]+f[1]*facDiff[2])+4.166190448976479*facDiff[3]*f[4]+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4nu+1.700840128541522*(alphaDrag[4]*f[10]+f[6]*alphaDrag[8])+1.732050807568877*(alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+f[3]*alphaDrag[4])+1.936491673103709*(alphaDrag[0]*f[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[9] += 16.20185174601965*(facDiff[3]*f[10]+facDiff[2]*f[6]+facDiff[1]*f[3]+facDiff[0]*f[2])*rdvxSq4nu+2.598076211353316*alphaDrag[2]*f[9]+1.322875655532295*alphaDrag[8]*f[8]+2.958039891549809*(alphaDrag[1]*f[7]+alphaDrag[0]*f[5])+1.322875655532295*alphaDrag[4]*f[4]+3.968626966596886*alphaDrag[2]*f[2]+1.322875655532295*(alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[10] += 0.8660254037844386*alphaDrag[2]*f[10]+(0.5163977794943223*alphaDrag[4]+0.8660254037844386*alphaDrag[0])*f[8]+(0.5163977794943223*f[4]+0.8660254037844386*f[0])*alphaDrag[8]+0.7606388292556648*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4]); 
  out[11] += (14.23024947075771*facDiff[2]*f[10]+(14.23024947075771*facDiff[3]+14.49137674618944*facDiff[1])*f[6]+16.20185174601965*(facDiff[0]*f[3]+facDiff[1]*f[2])+14.49137674618944*facDiff[2]*f[3])*rdvxSq4nu+2.598076211353316*alphaDrag[2]*f[11]+1.161895003862225*(alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+2.645751311064591*alphaDrag[4]*f[7]+2.958039891549809*(alphaDrag[0]*f[7]+alphaDrag[1]*f[5])+1.183215956619923*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+3.968626966596886*alphaDrag[2]*f[3]+1.322875655532295*(alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 

  return std::abs(0.25*alphaDrag[0]-0.2795084971874737*alphaDrag[4])+std::abs((1.616244071283538*facDiff[0]-1.807015805810503*facDiff[2])*rdvxSq4nu); 

} 
