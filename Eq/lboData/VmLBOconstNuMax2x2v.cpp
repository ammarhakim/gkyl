#include <VmLBOModDecl.h> 
double VmLBOconstNuVol2x2vMaxP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[2]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[2]*dxv[2]); 
  const double rdvy2nu = 2.0*nu/dxv[3]; 
  const double rdvySq4nu = 4.0*nu/(dxv[3]*dxv[3]); 

  double alpha_mid = 0.0; 
  double alpha_drag[10]; 
  double alpha_diffusion[5]; 

  // Expand rdv2nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.0*u[0]-4.0*w[2])*rdvx2nu; 
  alpha_drag[1] = 2.0*u[1]*rdvx2nu; 
  alpha_drag[2] = 2.0*u[2]*rdvx2nu; 
  alpha_drag[3] = -1.154700538379252*dxv[2]*rdvx2nu; 

  // x-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[0]); 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alpha_drag[5] = (2.0*u[3]-4.0*w[3])*rdvy2nu; 
  alpha_drag[6] = 2.0*u[4]*rdvy2nu; 
  alpha_drag[7] = 2.0*u[5]*rdvy2nu; 
  alpha_drag[9] = -1.154700538379252*dxv[3]*rdvy2nu; 

  // y-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[5]); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.0*vtSq[0]; 
  alpha_diffusion[1] = 2.0*vtSq[1]; 
  alpha_diffusion[2] = 2.0*vtSq[2]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.3333333333333333*alpha_diffusion[0]*rdvxSq4nu); 
  alpha_mid += std::abs(0.3333333333333333*alpha_diffusion[0]*rdvySq4nu); 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.4330127018922193*(alpha_drag[3]*f[3]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[4]*alpha_drag[9]+f[2]*alpha_drag[7]+f[1]*alpha_drag[6]+f[0]*alpha_drag[5]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol2x2vMaxP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[2]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[2]*dxv[2]); 
  const double rdvy2nu = 2.0*nu/dxv[3]; 
  const double rdvySq4nu = 4.0*nu/(dxv[3]*dxv[3]); 

  double alpha_mid = 0.0; 
  double alpha_drag[30]; 
  double alpha_diffusion[15]; 

  // Expand rdv2nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.0*u[0]-4.0*w[2])*rdvx2nu; 
  alpha_drag[1] = 2.0*u[1]*rdvx2nu; 
  alpha_drag[2] = 2.0*u[2]*rdvx2nu; 
  alpha_drag[3] = -1.154700538379252*dxv[2]*rdvx2nu; 
  alpha_drag[5] = 2.0*u[3]*rdvx2nu; 
  alpha_drag[11] = 2.0*u[4]*rdvx2nu; 
  alpha_drag[12] = 2.0*u[5]*rdvx2nu; 

  // x-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[0]-0.1397542485937369*(alpha_drag[12]+alpha_drag[11])); 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alpha_drag[15] = (2.0*u[6]-4.0*w[3])*rdvy2nu; 
  alpha_drag[16] = 2.0*u[7]*rdvy2nu; 
  alpha_drag[17] = 2.0*u[8]*rdvy2nu; 
  alpha_drag[19] = -1.154700538379252*dxv[3]*rdvy2nu; 
  alpha_drag[20] = 2.0*u[9]*rdvy2nu; 
  alpha_drag[26] = 2.0*u[10]*rdvy2nu; 
  alpha_drag[27] = 2.0*u[11]*rdvy2nu; 

  // y-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[15]-0.1397542485937369*(alpha_drag[27]+alpha_drag[26])); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.0*vtSq[0]; 
  alpha_diffusion[1] = 2.0*vtSq[1]; 
  alpha_diffusion[2] = 2.0*vtSq[2]; 
  alpha_diffusion[5] = 2.0*vtSq[3]; 
  alpha_diffusion[11] = 2.0*vtSq[4]; 
  alpha_diffusion[12] = 2.0*vtSq[5]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((0.45*alpha_diffusion[0]-0.5031152949374527*(alpha_diffusion[12]+alpha_diffusion[11]))*rdvxSq4nu); 
  alpha_mid += std::abs((0.45*alpha_diffusion[0]-0.5031152949374527*(alpha_diffusion[12]+alpha_diffusion[11]))*rdvySq4nu); 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.4330127018922193*(alpha_drag[12]*f[12]+alpha_drag[11]*f[11]+alpha_drag[5]*f[5]+alpha_drag[3]*f[3]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[12]*alpha_drag[27]+f[11]*alpha_drag[26]+f[5]*alpha_drag[20]+f[4]*alpha_drag[19]+f[2]*alpha_drag[17]+f[1]*alpha_drag[16]+f[0]*alpha_drag[15]); 
  out[6] += 0.3872983346207416*(alpha_drag[1]*f[11]+f[1]*alpha_drag[11])+0.4330127018922193*(alpha_drag[3]*f[6]+alpha_drag[2]*f[5]+f[2]*alpha_drag[5]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[7] += 0.3872983346207416*(alpha_drag[2]*f[12]+f[2]*alpha_drag[12])+0.4330127018922193*(alpha_drag[3]*f[7]+alpha_drag[1]*f[5]+f[1]*alpha_drag[5]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[8] += 0.3872983346207416*f[1]*alpha_drag[26]+0.4330127018922193*(f[2]*alpha_drag[20]+f[8]*alpha_drag[19]+f[5]*alpha_drag[17])+0.3872983346207416*f[11]*alpha_drag[16]+0.4330127018922193*(f[0]*alpha_drag[16]+f[1]*alpha_drag[15]); 
  out[9] += 0.3872983346207416*f[2]*alpha_drag[27]+0.4330127018922193*(f[1]*alpha_drag[20]+f[9]*alpha_drag[19])+0.3872983346207416*f[12]*alpha_drag[17]+0.4330127018922193*(f[0]*alpha_drag[17]+f[5]*alpha_drag[16]+f[2]*alpha_drag[15]); 
  out[10] += 0.4330127018922193*(f[10]*alpha_drag[19]+f[7]*alpha_drag[17]+f[6]*alpha_drag[16]+f[3]*alpha_drag[15]+alpha_drag[3]*f[10]+alpha_drag[2]*f[9]+alpha_drag[1]*f[8]+alpha_drag[0]*f[4]); 
  out[13] += 1.677050983124842*(alpha_diffusion[12]*f[12]+alpha_diffusion[11]*f[11]+alpha_diffusion[5]*f[5]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4nu+0.8660254037844386*alpha_drag[3]*f[13]+0.9682458365518543*(alpha_drag[2]*f[7]+alpha_drag[1]*f[6]+alpha_drag[0]*f[3]+f[0]*alpha_drag[3]); 
  out[14] += 1.677050983124842*(alpha_diffusion[12]*f[12]+alpha_diffusion[11]*f[11]+alpha_diffusion[5]*f[5]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvySq4nu+0.8660254037844386*f[14]*alpha_drag[19]+0.9682458365518543*(f[0]*alpha_drag[19]+f[9]*alpha_drag[17]+f[8]*alpha_drag[16]+f[4]*alpha_drag[15]); 

  return alpha_mid; 

} 
double VmLBOconstNuVol2x2vMaxP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nu:     diffusion coefficient (collisionality). 
  // u:    bulk velocity. 
  // vtSq: thermal speed squared. 
  // f:      Input distribution function.
  // out:    Incremented output 
  const double rdvx2nu = 2.0*nu/dxv[2]; 
  const double rdvxSq4nu = 4.0*nu/(dxv[2]*dxv[2]); 
  const double rdvy2nu = 2.0*nu/dxv[3]; 
  const double rdvySq4nu = 4.0*nu/(dxv[3]*dxv[3]); 

  double alpha_mid = 0.0; 
  double alpha_drag[70]; 
  double alpha_diffusion[35]; 

  // Expand rdv2nu*(vx-ux) in phase basis.
  alpha_drag[0] = (2.0*u[0]-4.0*w[2])*rdvx2nu; 
  alpha_drag[1] = 2.0*u[1]*rdvx2nu; 
  alpha_drag[2] = 2.0*u[2]*rdvx2nu; 
  alpha_drag[3] = -1.154700538379252*dxv[2]*rdvx2nu; 
  alpha_drag[5] = 2.0*u[3]*rdvx2nu; 
  alpha_drag[11] = 2.0*u[4]*rdvx2nu; 
  alpha_drag[12] = 2.0*u[5]*rdvx2nu; 
  alpha_drag[19] = 2.0*u[6]*rdvx2nu; 
  alpha_drag[20] = 2.0*u[7]*rdvx2nu; 
  alpha_drag[31] = 2.0*u[8]*rdvx2nu; 
  alpha_drag[32] = 2.0*u[9]*rdvx2nu; 

  // x-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[0]-0.1397542485937369*(alpha_drag[12]+alpha_drag[11])); 

  // Expand rdv2nu*(vy-uy) in phase basis.
  alpha_drag[35] = (2.0*u[10]-4.0*w[3])*rdvy2nu; 
  alpha_drag[36] = 2.0*u[11]*rdvy2nu; 
  alpha_drag[37] = 2.0*u[12]*rdvy2nu; 
  alpha_drag[39] = -1.154700538379252*dxv[3]*rdvy2nu; 
  alpha_drag[40] = 2.0*u[13]*rdvy2nu; 
  alpha_drag[46] = 2.0*u[14]*rdvy2nu; 
  alpha_drag[47] = 2.0*u[15]*rdvy2nu; 
  alpha_drag[54] = 2.0*u[16]*rdvy2nu; 
  alpha_drag[55] = 2.0*u[17]*rdvy2nu; 
  alpha_drag[66] = 2.0*u[18]*rdvy2nu; 
  alpha_drag[67] = 2.0*u[19]*rdvy2nu; 

  // y-drag contribution to the midpoint value used for CFL.
  alpha_mid += std::abs(0.125*alpha_drag[35]-0.1397542485937369*(alpha_drag[47]+alpha_drag[46])); 

  // Expand nu*vthSq in phase basis.
  alpha_diffusion[0] = 2.0*vtSq[0]; 
  alpha_diffusion[1] = 2.0*vtSq[1]; 
  alpha_diffusion[2] = 2.0*vtSq[2]; 
  alpha_diffusion[5] = 2.0*vtSq[3]; 
  alpha_diffusion[11] = 2.0*vtSq[4]; 
  alpha_diffusion[12] = 2.0*vtSq[5]; 
  alpha_diffusion[19] = 2.0*vtSq[6]; 
  alpha_diffusion[20] = 2.0*vtSq[7]; 
  alpha_diffusion[31] = 2.0*vtSq[8]; 
  alpha_diffusion[32] = 2.0*vtSq[9]; 

  // Diffusion contribution to the midpoint value used for CFL.
  alpha_mid += std::abs((0.5714285714285714*alpha_diffusion[0]-0.6388765649999399*(alpha_diffusion[12]+alpha_diffusion[11]))*rdvxSq4nu); 
  alpha_mid += std::abs((0.5714285714285714*alpha_diffusion[0]-0.6388765649999399*(alpha_diffusion[12]+alpha_diffusion[11]))*rdvySq4nu); 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.4330127018922193*(alpha_drag[32]*f[32]+alpha_drag[31]*f[31]+alpha_drag[20]*f[20]+alpha_drag[19]*f[19]+alpha_drag[12]*f[12]+alpha_drag[11]*f[11]+alpha_drag[5]*f[5]+alpha_drag[3]*f[3]+alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[32]*alpha_drag[67]+f[31]*alpha_drag[66]+f[20]*alpha_drag[55]+f[19]*alpha_drag[54]+f[12]*alpha_drag[47]+f[11]*alpha_drag[46]+f[5]*alpha_drag[40]+f[4]*alpha_drag[39]+f[2]*alpha_drag[37]+f[1]*alpha_drag[36]+f[0]*alpha_drag[35]); 
  out[6] += 0.3803194146278324*(alpha_drag[11]*f[31]+f[11]*alpha_drag[31])+0.4330127018922193*(alpha_drag[12]*f[20]+f[12]*alpha_drag[20])+0.3872983346207416*(alpha_drag[5]*f[19]+f[5]*alpha_drag[19]+alpha_drag[1]*f[11]+f[1]*alpha_drag[11])+0.4330127018922193*(alpha_drag[3]*f[6]+alpha_drag[2]*f[5]+f[2]*alpha_drag[5]+alpha_drag[0]*f[1]+f[0]*alpha_drag[1]); 
  out[7] += 0.3803194146278324*(alpha_drag[12]*f[32]+f[12]*alpha_drag[32])+0.3872983346207416*(alpha_drag[5]*f[20]+f[5]*alpha_drag[20])+0.4330127018922193*(alpha_drag[11]*f[19]+f[11]*alpha_drag[19])+0.3872983346207416*(alpha_drag[2]*f[12]+f[2]*alpha_drag[12])+0.4330127018922193*(alpha_drag[3]*f[7]+alpha_drag[1]*f[5]+f[1]*alpha_drag[5]+alpha_drag[0]*f[2]+f[0]*alpha_drag[2]); 
  out[8] += 0.3803194146278324*f[11]*alpha_drag[66]+0.4330127018922193*f[12]*alpha_drag[55]+0.3872983346207416*f[5]*alpha_drag[54]+0.4330127018922193*f[20]*alpha_drag[47]+0.3803194146278324*f[31]*alpha_drag[46]+0.3872983346207416*(f[1]*alpha_drag[46]+f[19]*alpha_drag[40])+0.4330127018922193*(f[2]*alpha_drag[40]+f[8]*alpha_drag[39]+f[5]*alpha_drag[37])+0.3872983346207416*f[11]*alpha_drag[36]+0.4330127018922193*(f[0]*alpha_drag[36]+f[1]*alpha_drag[35]); 
  out[9] += 0.3803194146278324*f[12]*alpha_drag[67]+0.3872983346207416*f[5]*alpha_drag[55]+0.4330127018922193*f[11]*alpha_drag[54]+(0.3803194146278324*f[32]+0.3872983346207416*f[2])*alpha_drag[47]+0.4330127018922193*f[19]*alpha_drag[46]+0.3872983346207416*f[20]*alpha_drag[40]+0.4330127018922193*(f[1]*alpha_drag[40]+f[9]*alpha_drag[39])+0.3872983346207416*f[12]*alpha_drag[37]+0.4330127018922193*(f[0]*alpha_drag[37]+f[5]*alpha_drag[36]+f[2]*alpha_drag[35]); 
  out[10] += 0.4330127018922193*(f[22]*alpha_drag[47]+f[21]*alpha_drag[46]+f[15]*alpha_drag[40]+f[10]*alpha_drag[39]+f[7]*alpha_drag[37]+f[6]*alpha_drag[36]+f[3]*alpha_drag[35]+alpha_drag[12]*f[26]+alpha_drag[11]*f[25]+alpha_drag[5]*f[16]+alpha_drag[3]*f[10]+alpha_drag[2]*f[9]+alpha_drag[1]*f[8]+alpha_drag[0]*f[4]); 
  out[13] += 1.677050983124842*(alpha_diffusion[32]*f[32]+alpha_diffusion[31]*f[31]+alpha_diffusion[20]*f[20]+alpha_diffusion[19]*f[19]+alpha_diffusion[12]*f[12]+alpha_diffusion[11]*f[11]+alpha_diffusion[5]*f[5]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvxSq4nu+0.9682458365518543*(alpha_drag[12]*f[22]+alpha_drag[11]*f[21]+alpha_drag[5]*f[15])+0.8660254037844386*alpha_drag[3]*f[13]+0.9682458365518543*(alpha_drag[2]*f[7]+alpha_drag[1]*f[6]+alpha_drag[0]*f[3]+f[0]*alpha_drag[3]); 
  out[14] += 1.677050983124842*(alpha_diffusion[32]*f[32]+alpha_diffusion[31]*f[31]+alpha_diffusion[20]*f[20]+alpha_diffusion[19]*f[19]+alpha_diffusion[12]*f[12]+alpha_diffusion[11]*f[11]+alpha_diffusion[5]*f[5]+alpha_diffusion[2]*f[2]+alpha_diffusion[1]*f[1]+alpha_diffusion[0]*f[0])*rdvySq4nu+0.9682458365518543*(f[26]*alpha_drag[47]+f[25]*alpha_drag[46]+f[16]*alpha_drag[40])+0.8660254037844386*f[14]*alpha_drag[39]+0.9682458365518543*(f[0]*alpha_drag[39]+f[9]*alpha_drag[37]+f[8]*alpha_drag[36]+f[4]*alpha_drag[35]); 
  out[15] += 0.3803194146278324*(alpha_drag[20]*f[32]+f[20]*alpha_drag[32]+alpha_drag[19]*f[31]+f[19]*alpha_drag[31])+(0.3464101615137755*alpha_drag[19]+0.3872983346207416*alpha_drag[2])*f[20]+0.3464101615137755*f[19]*alpha_drag[20]+0.3872983346207416*(f[2]*alpha_drag[20]+alpha_drag[1]*f[19]+f[1]*alpha_drag[19])+0.4330127018922193*alpha_drag[3]*f[15]+0.3872983346207416*(alpha_drag[5]*f[12]+f[5]*alpha_drag[12]+alpha_drag[5]*f[11]+f[5]*alpha_drag[11])+0.4330127018922193*(alpha_drag[0]*f[5]+f[0]*alpha_drag[5]+alpha_drag[1]*f[2]+f[1]*alpha_drag[2]); 
  out[16] += 0.3803194146278324*(f[20]*alpha_drag[67]+f[19]*alpha_drag[66])+(0.3803194146278324*f[32]+0.3464101615137755*f[19]+0.3872983346207416*f[2])*alpha_drag[55]+(0.3803194146278324*f[31]+0.3464101615137755*f[20])*alpha_drag[54]+0.3872983346207416*(f[1]*alpha_drag[54]+f[5]*(alpha_drag[47]+alpha_drag[46])+(f[12]+f[11])*alpha_drag[40])+0.4330127018922193*(f[0]*alpha_drag[40]+f[16]*alpha_drag[39])+(0.3872983346207416*f[20]+0.4330127018922193*f[1])*alpha_drag[37]+0.3872983346207416*f[19]*alpha_drag[36]+0.4330127018922193*(f[2]*alpha_drag[36]+f[5]*alpha_drag[35]); 
  out[17] += 0.3803194146278324*f[21]*alpha_drag[66]+0.4330127018922193*f[22]*alpha_drag[55]+0.3872983346207416*(f[15]*alpha_drag[54]+f[6]*alpha_drag[46])+0.4330127018922193*(f[7]*alpha_drag[40]+f[17]*alpha_drag[39]+f[15]*alpha_drag[37])+0.3872983346207416*f[21]*alpha_drag[36]+0.4330127018922193*(f[3]*alpha_drag[36]+f[6]*alpha_drag[35])+0.3803194146278324*f[25]*alpha_drag[31]+0.4330127018922193*alpha_drag[20]*f[26]+0.3872983346207416*(alpha_drag[1]*f[25]+f[16]*alpha_drag[19])+0.4330127018922193*(alpha_drag[3]*f[17]+alpha_drag[2]*f[16])+0.3872983346207416*f[8]*alpha_drag[11]+0.4330127018922193*(alpha_drag[5]*f[9]+alpha_drag[0]*f[8]+alpha_drag[1]*f[4]); 
  out[18] += 0.3803194146278324*f[22]*alpha_drag[67]+0.3872983346207416*f[15]*alpha_drag[55]+0.4330127018922193*f[21]*alpha_drag[54]+0.3872983346207416*f[7]*alpha_drag[47]+0.4330127018922193*(f[6]*alpha_drag[40]+f[18]*alpha_drag[39])+0.3872983346207416*f[22]*alpha_drag[37]+0.4330127018922193*(f[3]*alpha_drag[37]+f[15]*alpha_drag[36]+f[7]*alpha_drag[35])+f[26]*(0.3803194146278324*alpha_drag[32]+0.3872983346207416*alpha_drag[2])+0.4330127018922193*alpha_drag[19]*f[25]+0.3872983346207416*f[16]*alpha_drag[20]+0.4330127018922193*(alpha_drag[3]*f[18]+alpha_drag[1]*f[16])+0.3872983346207416*f[9]*alpha_drag[12]+0.4330127018922193*(alpha_drag[0]*f[9]+alpha_drag[5]*f[8]+alpha_drag[2]*f[4]); 
  out[21] += 0.2581988897471612*alpha_drag[31]*f[31]+0.3803194146278324*(alpha_drag[1]*f[31]+f[1]*alpha_drag[31])+0.4330127018922193*alpha_drag[3]*f[21]+0.3872983346207416*alpha_drag[20]*f[20]+0.276641667586244*alpha_drag[19]*f[19]+0.4330127018922193*(alpha_drag[2]*f[19]+f[2]*alpha_drag[19])+0.276641667586244*alpha_drag[11]*f[11]+0.4330127018922193*(alpha_drag[0]*f[11]+f[0]*alpha_drag[11])+0.3872983346207416*(alpha_drag[5]*f[5]+alpha_drag[1]*f[1]); 
  out[22] += 0.2581988897471612*alpha_drag[32]*f[32]+0.3803194146278324*(alpha_drag[2]*f[32]+f[2]*alpha_drag[32])+0.4330127018922193*alpha_drag[3]*f[22]+0.276641667586244*alpha_drag[20]*f[20]+0.4330127018922193*(alpha_drag[1]*f[20]+f[1]*alpha_drag[20])+0.3872983346207416*alpha_drag[19]*f[19]+0.276641667586244*alpha_drag[12]*f[12]+0.4330127018922193*(alpha_drag[0]*f[12]+f[0]*alpha_drag[12])+0.3872983346207416*(alpha_drag[5]*f[5]+alpha_drag[2]*f[2]); 
  out[23] += (1.472970759092948*(alpha_diffusion[11]*f[31]+f[11]*alpha_diffusion[31])+1.677050983124842*(alpha_diffusion[12]*f[20]+f[12]*alpha_diffusion[20])+1.5*(alpha_diffusion[5]*f[19]+f[5]*alpha_diffusion[19]+alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[2]*f[5]+f[2]*alpha_diffusion[5]+alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvxSq4nu+0.8504200642707612*f[21]*alpha_drag[31]+0.8660254037844386*alpha_drag[3]*f[23]+0.9682458365518543*alpha_drag[20]*f[22]+0.8660254037844386*alpha_drag[1]*f[21]+f[15]*(0.8660254037844386*alpha_drag[19]+0.9682458365518543*alpha_drag[2])+0.8660254037844386*f[6]*alpha_drag[11]+0.9682458365518543*(alpha_drag[5]*f[7]+alpha_drag[0]*f[6]+alpha_drag[1]*f[3]+f[1]*alpha_drag[3]); 
  out[24] += (1.472970759092948*(alpha_diffusion[12]*f[32]+f[12]*alpha_diffusion[32])+1.5*(alpha_diffusion[5]*f[20]+f[5]*alpha_diffusion[20])+1.677050983124842*(alpha_diffusion[11]*f[19]+f[11]*alpha_diffusion[19])+1.5*(alpha_diffusion[2]*f[12]+f[2]*alpha_diffusion[12])+1.677050983124842*(alpha_diffusion[1]*f[5]+f[1]*alpha_diffusion[5]+alpha_diffusion[0]*f[2]+f[0]*alpha_diffusion[2]))*rdvxSq4nu+0.8504200642707612*f[22]*alpha_drag[32]+0.8660254037844386*(alpha_drag[3]*f[24]+alpha_drag[2]*f[22])+0.9682458365518543*alpha_drag[19]*f[21]+f[15]*(0.8660254037844386*alpha_drag[20]+0.9682458365518543*alpha_drag[1])+0.8660254037844386*f[7]*alpha_drag[12]+0.9682458365518543*(alpha_drag[0]*f[7]+alpha_drag[5]*f[6]+alpha_drag[2]*f[3]+f[2]*alpha_drag[3]); 
  out[25] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_drag[66]+0.3872983346207416*f[20]*alpha_drag[55]+(0.276641667586244*f[19]+0.4330127018922193*f[2])*alpha_drag[54]+(0.276641667586244*f[11]+0.4330127018922193*f[0])*alpha_drag[46]+0.3872983346207416*f[5]*alpha_drag[40]+0.4330127018922193*(f[25]*alpha_drag[39]+f[19]*alpha_drag[37])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_drag[36]+0.4330127018922193*f[11]*alpha_drag[35]; 
  out[26] += (0.2581988897471612*f[32]+0.3803194146278324*f[2])*alpha_drag[67]+(0.276641667586244*f[20]+0.4330127018922193*f[1])*alpha_drag[55]+0.3872983346207416*f[19]*alpha_drag[54]+(0.276641667586244*f[12]+0.4330127018922193*f[0])*alpha_drag[47]+0.3872983346207416*f[5]*alpha_drag[40]+0.4330127018922193*f[26]*alpha_drag[39]+(0.3803194146278324*f[32]+0.3872983346207416*f[2])*alpha_drag[37]+0.4330127018922193*(f[20]*alpha_drag[36]+f[12]*alpha_drag[35]); 
  out[27] += 1.677050983124842*(alpha_diffusion[12]*f[26]+alpha_diffusion[11]*f[25]+alpha_diffusion[5]*f[16]+alpha_diffusion[2]*f[9]+alpha_diffusion[1]*f[8]+alpha_diffusion[0]*f[4])*rdvxSq4nu+0.4330127018922193*(f[27]*alpha_drag[39]+f[24]*alpha_drag[37]+f[23]*alpha_drag[36]+f[13]*alpha_drag[35])+0.8660254037844386*alpha_drag[3]*f[27]+0.9682458365518543*(alpha_drag[2]*f[18]+alpha_drag[1]*f[17]+alpha_drag[0]*f[10]+alpha_drag[3]*f[4]); 
  out[28] += (1.472970759092948*(alpha_diffusion[11]*f[31]+f[11]*alpha_diffusion[31])+1.677050983124842*(alpha_diffusion[12]*f[20]+f[12]*alpha_diffusion[20])+1.5*(alpha_diffusion[5]*f[19]+f[5]*alpha_diffusion[19]+alpha_diffusion[1]*f[11]+f[1]*alpha_diffusion[11])+1.677050983124842*(alpha_diffusion[2]*f[5]+f[2]*alpha_diffusion[5]+alpha_diffusion[0]*f[1]+f[0]*alpha_diffusion[1]))*rdvySq4nu+0.8504200642707612*f[25]*alpha_drag[66]+0.9682458365518543*f[26]*alpha_drag[55]+0.8660254037844386*(f[16]*alpha_drag[54]+f[8]*alpha_drag[46])+0.9682458365518543*f[9]*alpha_drag[40]+0.8660254037844386*f[28]*alpha_drag[39]+0.9682458365518543*(f[1]*alpha_drag[39]+f[16]*alpha_drag[37])+0.8660254037844386*f[25]*alpha_drag[36]+0.9682458365518543*(f[4]*alpha_drag[36]+f[8]*alpha_drag[35]); 
  out[29] += (1.472970759092948*(alpha_diffusion[12]*f[32]+f[12]*alpha_diffusion[32])+1.5*(alpha_diffusion[5]*f[20]+f[5]*alpha_diffusion[20])+1.677050983124842*(alpha_diffusion[11]*f[19]+f[11]*alpha_diffusion[19])+1.5*(alpha_diffusion[2]*f[12]+f[2]*alpha_diffusion[12])+1.677050983124842*(alpha_diffusion[1]*f[5]+f[1]*alpha_diffusion[5]+alpha_diffusion[0]*f[2]+f[0]*alpha_diffusion[2]))*rdvySq4nu+0.8504200642707612*f[26]*alpha_drag[67]+0.8660254037844386*f[16]*alpha_drag[55]+0.9682458365518543*f[25]*alpha_drag[54]+0.8660254037844386*f[9]*alpha_drag[47]+0.9682458365518543*f[8]*alpha_drag[40]+(0.8660254037844386*f[29]+0.9682458365518543*f[2])*alpha_drag[39]+0.8660254037844386*f[26]*alpha_drag[37]+0.9682458365518543*(f[4]*alpha_drag[37]+f[16]*alpha_drag[36]+f[9]*alpha_drag[35]); 
  out[30] += 1.677050983124842*(alpha_diffusion[12]*f[22]+alpha_diffusion[11]*f[21]+alpha_diffusion[5]*f[15]+alpha_diffusion[2]*f[7]+alpha_diffusion[1]*f[6]+alpha_diffusion[0]*f[3])*rdvySq4nu+0.8660254037844386*f[30]*alpha_drag[39]+0.9682458365518543*(f[3]*alpha_drag[39]+f[18]*alpha_drag[37]+f[17]*alpha_drag[36]+f[10]*alpha_drag[35])+0.4330127018922193*(alpha_drag[3]*f[30]+alpha_drag[2]*f[29]+alpha_drag[1]*f[28]+alpha_drag[0]*f[14]); 
  out[33] += 5.7282196186948*(alpha_diffusion[12]*f[22]+alpha_diffusion[11]*f[21]+alpha_diffusion[5]*f[15]+alpha_diffusion[2]*f[7]+alpha_diffusion[1]*f[6]+alpha_diffusion[0]*f[3])*rdvxSq4nu+1.299038105676658*alpha_drag[3]*f[33]+0.6614378277661477*(alpha_drag[32]*f[32]+alpha_drag[31]*f[31])+1.479019945774904*(alpha_drag[2]*f[24]+alpha_drag[1]*f[23])+0.6614378277661477*(alpha_drag[20]*f[20]+alpha_drag[19]*f[19])+1.479019945774904*alpha_drag[0]*f[13]+0.6614378277661477*(alpha_drag[12]*f[12]+alpha_drag[11]*f[11]+alpha_drag[5]*f[5])+1.984313483298443*alpha_drag[3]*f[3]+0.6614378277661477*(alpha_drag[2]*f[2]+alpha_drag[1]*f[1]+alpha_drag[0]*f[0]); 
  out[34] += 5.7282196186948*(alpha_diffusion[12]*f[26]+alpha_diffusion[11]*f[25]+alpha_diffusion[5]*f[16]+alpha_diffusion[2]*f[9]+alpha_diffusion[1]*f[8]+alpha_diffusion[0]*f[4])*rdvySq4nu+0.6614378277661477*(f[32]*alpha_drag[67]+f[31]*alpha_drag[66]+f[20]*alpha_drag[55]+f[19]*alpha_drag[54]+f[12]*alpha_drag[47]+f[11]*alpha_drag[46]+f[5]*alpha_drag[40])+(1.299038105676658*f[34]+1.984313483298443*f[4])*alpha_drag[39]+(1.479019945774904*f[29]+0.6614378277661477*f[2])*alpha_drag[37]+(1.479019945774904*f[28]+0.6614378277661477*f[1])*alpha_drag[36]+(1.479019945774904*f[14]+0.6614378277661477*f[0])*alpha_drag[35]; 

  return alpha_mid; 

} 
