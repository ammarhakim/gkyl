#include <VmLBOModDecl.h> 
double VmLBOconstNuVol2x3vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates. 
  // dxv[5]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[2]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[2]*dxv[2]); 
  const double rdvy2nu = 2.0*nu/dxv[3]; 
  const double rdvySq4nu = 4.0*nu/(dxv[3]*dxv[3]); 
  const double rdvz2nu = 2.0*nu/dxv[4]; 
  const double rdvzSq4nu = 4.0*nu/(dxv[4]*dxv[4]); 

  double alphaDrag[18]; 
  // Expand rdv2nu*(vx-ux) in phase basis.
  alphaDrag[0] = (2.828427124746191*u[0]-5.656854249492382*w[2])*rdvx2nu; 
  alphaDrag[1] = 2.828427124746191*u[1]*rdvx2nu; 
  alphaDrag[2] = 2.828427124746191*u[2]*rdvx2nu; 
  alphaDrag[3] = -1.632993161855453*dxv[2]*rdvx2nu; 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alphaDrag[6] = (2.828427124746191*u[3]-5.656854249492382*w[3])*rdvy2nu; 
  alphaDrag[7] = 2.828427124746191*u[4]*rdvy2nu; 
  alphaDrag[8] = 2.828427124746191*u[5]*rdvy2nu; 
  alphaDrag[10] = -1.632993161855453*dxv[3]*rdvy2nu; 

  // Expand rdv2nu*(vz-uz) in phase basis.
  alphaDrag[12] = (2.828427124746191*u[6]-5.656854249492382*w[4])*rdvz2nu; 
  alphaDrag[13] = 2.828427124746191*u[7]*rdvz2nu; 
  alphaDrag[14] = 2.828427124746191*u[8]*rdvz2nu; 
  alphaDrag[17] = -1.632993161855453*dxv[4]*rdvz2nu; 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.3061862178478971*(alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[4]*alphaDrag[10]+f[2]*alphaDrag[8]+f[1]*alphaDrag[7]+f[0]*alphaDrag[6]); 
  out[5] += 0.3061862178478971*(f[5]*alphaDrag[17]+f[2]*alphaDrag[14]+f[1]*alphaDrag[13]+f[0]*alphaDrag[12]); 

  return std::abs(0.0883883476483184*alphaDrag[0])+std::abs(0.0883883476483184*alphaDrag[6])+std::abs(0.0883883476483184*alphaDrag[12])+std::abs(0.6666666666666666*vtSq[0]*rdvxSq4nu)+std::abs(0.6666666666666666*vtSq[0]*rdvySq4nu)+std::abs(0.6666666666666666*vtSq[0]*rdvzSq4nu); 

} 
double VmLBOconstNuVol2x3vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[5]:   Cell-center coordinates. 
  // dxv[5]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[2]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[2]*dxv[2]); 
  const double rdvy2nu = 2.0*nu/dxv[3]; 
  const double rdvySq4nu = 4.0*nu/(dxv[3]*dxv[3]); 
  const double rdvz2nu = 2.0*nu/dxv[4]; 
  const double rdvzSq4nu = 4.0*nu/(dxv[4]*dxv[4]); 

  double alphaDrag[63]; 
  // Expand rdv2nu*(vx-ux) in phase basis.
  alphaDrag[0] = (2.828427124746191*u[0]-5.656854249492382*w[2])*rdvx2nu; 
  alphaDrag[1] = 2.828427124746191*u[1]*rdvx2nu; 
  alphaDrag[2] = 2.828427124746191*u[2]*rdvx2nu; 
  alphaDrag[3] = -1.632993161855453*dxv[2]*rdvx2nu; 
  alphaDrag[6] = 2.828427124746191*u[3]*rdvx2nu; 
  alphaDrag[16] = 2.828427124746191*u[4]*rdvx2nu; 
  alphaDrag[17] = 2.828427124746191*u[5]*rdvx2nu; 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alphaDrag[21] = (2.828427124746191*u[6]-5.656854249492382*w[3])*rdvy2nu; 
  alphaDrag[22] = 2.828427124746191*u[7]*rdvy2nu; 
  alphaDrag[23] = 2.828427124746191*u[8]*rdvy2nu; 
  alphaDrag[25] = -1.632993161855453*dxv[3]*rdvy2nu; 
  alphaDrag[27] = 2.828427124746191*u[9]*rdvy2nu; 
  alphaDrag[37] = 2.828427124746191*u[10]*rdvy2nu; 
  alphaDrag[38] = 2.828427124746191*u[11]*rdvy2nu; 

  // Expand rdv2nu*(vz-uz) in phase basis.
  alphaDrag[42] = (2.828427124746191*u[12]-5.656854249492382*w[4])*rdvz2nu; 
  alphaDrag[43] = 2.828427124746191*u[13]*rdvz2nu; 
  alphaDrag[44] = 2.828427124746191*u[14]*rdvz2nu; 
  alphaDrag[47] = -1.632993161855453*dxv[4]*rdvz2nu; 
  alphaDrag[48] = 2.828427124746191*u[15]*rdvz2nu; 
  alphaDrag[58] = 2.828427124746191*u[16]*rdvz2nu; 
  alphaDrag[59] = 2.828427124746191*u[17]*rdvz2nu; 

  double facDiff[6]; 
  // Expand nu*vthSq in phase basis.
  facDiff[0] = vtSq[0]; 
  facDiff[1] = vtSq[1]; 
  facDiff[2] = vtSq[2]; 
  facDiff[3] = vtSq[3]; 
  facDiff[4] = vtSq[4]; 
  facDiff[5] = vtSq[5]; 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.3061862178478971*(alphaDrag[17]*f[17]+alphaDrag[16]*f[16]+alphaDrag[6]*f[6]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[17]*alphaDrag[38]+f[16]*alphaDrag[37]+f[6]*alphaDrag[27]+f[4]*alphaDrag[25]+f[2]*alphaDrag[23]+f[1]*alphaDrag[22]+f[0]*alphaDrag[21]); 
  out[5] += 0.3061862178478971*(f[17]*alphaDrag[59]+f[16]*alphaDrag[58]+f[6]*alphaDrag[48]+f[5]*alphaDrag[47]+f[2]*alphaDrag[44]+f[1]*alphaDrag[43]+f[0]*alphaDrag[42]); 
  out[7] += 0.273861278752583*(alphaDrag[1]*f[16]+f[1]*alphaDrag[16])+0.3061862178478971*(alphaDrag[3]*f[7]+alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[8] += 0.273861278752583*(alphaDrag[2]*f[17]+f[2]*alphaDrag[17])+0.3061862178478971*(alphaDrag[3]*f[8]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 0.273861278752583*f[1]*alphaDrag[37]+0.3061862178478971*(f[2]*alphaDrag[27]+f[9]*alphaDrag[25]+f[6]*alphaDrag[23])+0.273861278752583*f[16]*alphaDrag[22]+0.3061862178478971*(f[0]*alphaDrag[22]+f[1]*alphaDrag[21]); 
  out[10] += 0.273861278752583*f[2]*alphaDrag[38]+0.3061862178478971*(f[1]*alphaDrag[27]+f[10]*alphaDrag[25])+0.273861278752583*f[17]*alphaDrag[23]+0.3061862178478971*(f[0]*alphaDrag[23]+f[6]*alphaDrag[22]+f[2]*alphaDrag[21]); 
  out[11] += 0.3061862178478971*(f[11]*alphaDrag[25]+f[8]*alphaDrag[23]+f[7]*alphaDrag[22]+f[3]*alphaDrag[21]+alphaDrag[3]*f[11]+alphaDrag[2]*f[10]+alphaDrag[1]*f[9]+alphaDrag[0]*f[4]); 
  out[12] += 0.273861278752583*f[1]*alphaDrag[58]+0.3061862178478971*(f[2]*alphaDrag[48]+f[12]*alphaDrag[47]+f[6]*alphaDrag[44])+0.273861278752583*f[16]*alphaDrag[43]+0.3061862178478971*(f[0]*alphaDrag[43]+f[1]*alphaDrag[42]); 
  out[13] += 0.273861278752583*f[2]*alphaDrag[59]+0.3061862178478971*(f[1]*alphaDrag[48]+f[13]*alphaDrag[47])+0.273861278752583*f[17]*alphaDrag[44]+0.3061862178478971*(f[0]*alphaDrag[44]+f[6]*alphaDrag[43]+f[2]*alphaDrag[42]); 
  out[14] += 0.3061862178478971*(f[14]*alphaDrag[47]+f[8]*alphaDrag[44]+f[7]*alphaDrag[43]+f[3]*alphaDrag[42]+alphaDrag[3]*f[14]+alphaDrag[2]*f[13]+alphaDrag[1]*f[12]+alphaDrag[0]*f[5]); 
  out[15] += 0.3061862178478971*(f[15]*alphaDrag[47]+f[10]*alphaDrag[44]+f[9]*alphaDrag[43]+f[4]*alphaDrag[42]+f[15]*alphaDrag[25]+f[13]*alphaDrag[23]+f[12]*alphaDrag[22]+f[5]*alphaDrag[21]); 
  out[18] += 3.354101966249685*(facDiff[5]*f[17]+facDiff[4]*f[16]+facDiff[3]*f[6]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4nu+0.6123724356957944*alphaDrag[3]*f[18]+0.6846531968814573*(alphaDrag[2]*f[8]+alphaDrag[1]*f[7]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[19] += 3.354101966249685*(facDiff[5]*f[17]+facDiff[4]*f[16]+facDiff[3]*f[6]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4nu+0.6123724356957944*f[19]*alphaDrag[25]+0.6846531968814573*(f[0]*alphaDrag[25]+f[10]*alphaDrag[23]+f[9]*alphaDrag[22]+f[4]*alphaDrag[21]); 
  out[20] += 3.354101966249685*(facDiff[5]*f[17]+facDiff[4]*f[16]+facDiff[3]*f[6]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvzSq4nu+0.6123724356957944*f[20]*alphaDrag[47]+0.6846531968814573*(f[0]*alphaDrag[47]+f[13]*alphaDrag[44]+f[12]*alphaDrag[43]+f[5]*alphaDrag[42]); 

  return std::abs(0.0883883476483184*alphaDrag[0]-0.09882117688026182*(alphaDrag[17]+alphaDrag[16]))+std::abs(0.0883883476483184*alphaDrag[21]-0.09882117688026182*(alphaDrag[38]+alphaDrag[37]))+std::abs(0.0883883476483184*alphaDrag[42]-0.09882117688026182*(alphaDrag[59]+alphaDrag[58]))+std::abs((0.9*facDiff[0]-1.006230589874905*(facDiff[5]+facDiff[4]))*rdvxSq4nu)+std::abs((0.9*facDiff[0]-1.006230589874905*(facDiff[5]+facDiff[4]))*rdvySq4nu)+std::abs((0.9*facDiff[0]-1.006230589874905*(facDiff[5]+facDiff[4]))*rdvzSq4nu); 

} 
