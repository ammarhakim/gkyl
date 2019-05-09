#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x2vSerP1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 

  double alphaDrag[16]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.828427124746191*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.0*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -0.8164965809277261*dxv[1]*nuSum*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[8] = (2.0*nuUSum[2]-2.828427124746191*w[2]*nuSum)*rdvy2; 
  alphaDrag[9] = 2.0*nuUSum[3]*rdvy2; 
  alphaDrag[11] = -0.8164965809277261*dxv[2]*nuSum*rdvy2; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[3]*alphaDrag[11]+f[1]*alphaDrag[9]+f[0]*alphaDrag[8]); 
  out[4] += 0.6123724356957944*(alphaDrag[2]*f[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.6123724356957944*(f[5]*alphaDrag[11]+f[0]*alphaDrag[9]+f[1]*alphaDrag[8]); 
  out[6] += 0.6123724356957944*(f[6]*alphaDrag[11]+f[4]*alphaDrag[9]+f[2]*alphaDrag[8]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[7] += 0.6123724356957944*(f[7]*alphaDrag[11]+f[2]*alphaDrag[9]+f[4]*alphaDrag[8]+alphaDrag[2]*f[7]+alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 

  return std::abs(0.1767766952966368*alphaDrag[0])+std::abs(0.1767766952966368*alphaDrag[8])+std::abs(0.9428090415820636*nuVtSqSum[0]*rdvxSq4)+std::abs(0.9428090415820636*nuVtSqSum[0]*rdvySq4); 

} 
double VmLBOconstNuVol1x2vSerP2(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 

  double alphaDrag[40]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.828427124746191*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.0*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -0.8164965809277261*dxv[1]*nuSum*rdvx2; 
  alphaDrag[7] = 2.0*nuUSum[2]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[20] = (2.0*nuUSum[3]-2.828427124746191*w[2]*nuSum)*rdvy2; 
  alphaDrag[21] = 2.0*nuUSum[4]*rdvy2; 
  alphaDrag[23] = -0.8164965809277261*dxv[2]*nuSum*rdvy2; 
  alphaDrag[27] = 2.0*nuUSum[5]*rdvy2; 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alphaDrag[7]*f[7]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*alphaDrag[27]+f[3]*alphaDrag[23]+f[1]*alphaDrag[21]+f[0]*alphaDrag[20]); 
  out[4] += 0.5477225575051661*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.6123724356957944*(alphaDrag[2]*f[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.5477225575051661*f[1]*alphaDrag[27]+0.6123724356957944*f[5]*alphaDrag[23]+0.5477225575051661*f[7]*alphaDrag[21]+0.6123724356957944*(f[0]*alphaDrag[21]+f[1]*alphaDrag[20]); 
  out[6] += 0.6123724356957944*(f[11]*alphaDrag[27]+f[6]*alphaDrag[23]+f[4]*alphaDrag[21]+f[2]*alphaDrag[20]+alphaDrag[7]*f[13]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[8] += 4.743416490252569*(facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+1.369306393762915*alphaDrag[7]*f[11]+1.224744871391589*alphaDrag[2]*f[8]+1.369306393762915*(alphaDrag[1]*f[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 4.743416490252569*(facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+1.369306393762915*f[13]*alphaDrag[27]+1.224744871391589*f[9]*alphaDrag[23]+1.369306393762915*(f[0]*alphaDrag[23]+f[5]*alphaDrag[21]+f[3]*alphaDrag[20]); 
  out[10] += 0.5477225575051661*f[4]*alphaDrag[27]+0.6123724356957944*f[10]*alphaDrag[23]+0.5477225575051661*f[11]*alphaDrag[21]+0.6123724356957944*(f[2]*alphaDrag[21]+f[4]*alphaDrag[20])+0.5477225575051661*alphaDrag[1]*f[13]+0.6123724356957944*alphaDrag[2]*f[10]+0.5477225575051661*f[5]*alphaDrag[7]+0.6123724356957944*(alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[11] += 0.6123724356957944*alphaDrag[2]*f[11]+0.3912303982179757*alphaDrag[7]*f[7]+0.6123724356957944*(alphaDrag[0]*f[7]+f[0]*alphaDrag[7])+0.5477225575051661*alphaDrag[1]*f[1]; 
  out[12] += (4.242640687119286*(facDiff[1]*f[7]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4+1.224744871391589*(alphaDrag[2]*f[12]+alphaDrag[1]*f[11]+f[4]*alphaDrag[7])+1.369306393762915*(alphaDrag[0]*f[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[13] += 0.3912303982179757*f[7]*alphaDrag[27]+0.6123724356957944*(f[0]*alphaDrag[27]+f[13]*alphaDrag[23])+0.5477225575051661*f[1]*alphaDrag[21]+0.6123724356957944*f[7]*alphaDrag[20]; 
  out[14] += 4.743416490252569*(facDiff[2]*f[13]+facDiff[1]*f[5]+facDiff[0]*f[3])*rdvxSq4+0.6123724356957944*(f[14]*alphaDrag[23]+f[12]*alphaDrag[21]+f[8]*alphaDrag[20])+1.369306393762915*alphaDrag[7]*f[17]+1.224744871391589*alphaDrag[2]*f[14]+1.369306393762915*(alphaDrag[1]*f[10]+alphaDrag[0]*f[6]+alphaDrag[2]*f[3]); 
  out[15] += (4.242640687119286*(facDiff[1]*f[7]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvySq4+1.224744871391589*f[5]*alphaDrag[27]+(1.224744871391589*f[15]+1.369306393762915*f[1])*alphaDrag[23]+1.224744871391589*f[13]*alphaDrag[21]+1.369306393762915*(f[3]*alphaDrag[21]+f[5]*alphaDrag[20]); 
  out[16] += 4.743416490252569*(facDiff[2]*f[11]+facDiff[1]*f[4]+facDiff[0]*f[2])*rdvySq4+1.369306393762915*f[17]*alphaDrag[27]+1.224744871391589*f[16]*alphaDrag[23]+1.369306393762915*(f[2]*alphaDrag[23]+f[10]*alphaDrag[21]+f[6]*alphaDrag[20])+0.6123724356957944*(alphaDrag[2]*f[16]+alphaDrag[1]*f[15]+alphaDrag[0]*f[9]); 
  out[17] += 0.3912303982179757*f[11]*alphaDrag[27]+0.6123724356957944*(f[2]*alphaDrag[27]+f[17]*alphaDrag[23])+0.5477225575051661*f[4]*alphaDrag[21]+0.6123724356957944*(f[11]*alphaDrag[20]+alphaDrag[2]*f[17])+0.3912303982179757*alphaDrag[7]*f[13]+0.6123724356957944*(alphaDrag[0]*f[13]+f[3]*alphaDrag[7])+0.5477225575051661*alphaDrag[1]*f[5]; 
  out[18] += (4.242640687119286*facDiff[1]*f[13]+4.743416490252569*(facDiff[0]*f[5]+facDiff[1]*f[3])+4.242640687119286*facDiff[2]*f[5])*rdvxSq4+0.5477225575051661*f[12]*alphaDrag[27]+0.6123724356957944*(f[18]*alphaDrag[23]+f[8]*alphaDrag[21]+f[12]*alphaDrag[20])+1.224744871391589*(alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[7]*f[10])+1.369306393762915*(alphaDrag[0]*f[10]+alphaDrag[1]*f[6]+alphaDrag[2]*f[5]); 
  out[19] += (4.242640687119286*facDiff[1]*f[11]+4.743416490252569*(facDiff[0]*f[4]+facDiff[1]*f[2])+4.242640687119286*facDiff[2]*f[4])*rdvySq4+1.224744871391589*f[10]*alphaDrag[27]+(1.224744871391589*f[19]+1.369306393762915*f[4])*alphaDrag[23]+1.224744871391589*f[17]*alphaDrag[21]+1.369306393762915*(f[6]*alphaDrag[21]+f[10]*alphaDrag[20])+0.6123724356957944*alphaDrag[2]*f[19]+0.5477225575051661*alphaDrag[7]*f[15]+0.6123724356957944*(alphaDrag[0]*f[15]+alphaDrag[1]*f[9]); 

  return std::abs(0.1767766952966368*alphaDrag[0]-0.1976423537605236*alphaDrag[7])+std::abs(0.1767766952966368*alphaDrag[20]-0.1976423537605236*alphaDrag[27])+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvxSq4)+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvySq4); 

} 
double VmLBOconstNuVol1x2vSerP3(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 

  double alphaDrag[64]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.828427124746191*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.0*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -0.8164965809277261*dxv[1]*nuSum*rdvx2; 
  alphaDrag[7] = 2.0*nuUSum[2]*rdvx2; 
  alphaDrag[17] = 2.0*nuUSum[3]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[32] = (2.0*nuUSum[4]-2.828427124746191*w[2]*nuSum)*rdvy2; 
  alphaDrag[33] = 2.0*nuUSum[5]*rdvy2; 
  alphaDrag[35] = -0.8164965809277261*dxv[2]*nuSum*rdvy2; 
  alphaDrag[39] = 2.0*nuUSum[6]*rdvy2; 
  alphaDrag[49] = 2.0*nuUSum[7]*rdvy2; 

  double facDiff[4]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 
  facDiff[3] = nuVtSqSum[3]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alphaDrag[17]*f[17]+alphaDrag[7]*f[7]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[17]*alphaDrag[49]+f[7]*alphaDrag[39]+f[3]*alphaDrag[35]+f[1]*alphaDrag[33]+f[0]*alphaDrag[32]); 
  out[4] += 0.537852874200477*(alphaDrag[7]*f[17]+f[7]*alphaDrag[17])+0.5477225575051661*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.6123724356957944*(alphaDrag[2]*f[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.537852874200477*f[7]*alphaDrag[49]+(0.537852874200477*f[17]+0.5477225575051661*f[1])*alphaDrag[39]+0.6123724356957944*f[5]*alphaDrag[35]+0.5477225575051661*f[7]*alphaDrag[33]+0.6123724356957944*(f[0]*alphaDrag[33]+f[1]*alphaDrag[32]); 
  out[6] += 0.6123724356957944*(f[23]*alphaDrag[49]+f[11]*alphaDrag[39]+f[6]*alphaDrag[35]+f[4]*alphaDrag[33]+f[2]*alphaDrag[32]+alphaDrag[17]*f[25]+alphaDrag[7]*f[13]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[8] += 4.743416490252569*(facDiff[3]*f[17]+facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+1.369306393762915*(alphaDrag[17]*f[23]+alphaDrag[7]*f[11])+1.224744871391589*alphaDrag[2]*f[8]+1.369306393762915*(alphaDrag[1]*f[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 4.743416490252569*(facDiff[3]*f[17]+facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+1.369306393762915*(f[25]*alphaDrag[49]+f[13]*alphaDrag[39])+1.224744871391589*f[9]*alphaDrag[35]+1.369306393762915*(f[0]*alphaDrag[35]+f[5]*alphaDrag[33]+f[3]*alphaDrag[32]); 
  out[10] += 0.537852874200477*f[11]*alphaDrag[49]+(0.537852874200477*f[23]+0.5477225575051661*f[4])*alphaDrag[39]+0.6123724356957944*f[10]*alphaDrag[35]+0.5477225575051661*f[11]*alphaDrag[33]+0.6123724356957944*(f[2]*alphaDrag[33]+f[4]*alphaDrag[32])+0.537852874200477*alphaDrag[7]*f[25]+f[13]*(0.537852874200477*alphaDrag[17]+0.5477225575051661*alphaDrag[1])+0.6123724356957944*alphaDrag[2]*f[10]+0.5477225575051661*f[5]*alphaDrag[7]+0.6123724356957944*(alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[11] += 0.3651483716701108*alphaDrag[17]*f[17]+0.537852874200477*(alphaDrag[1]*f[17]+f[1]*alphaDrag[17])+0.6123724356957944*alphaDrag[2]*f[11]+0.3912303982179757*alphaDrag[7]*f[7]+0.6123724356957944*(alphaDrag[0]*f[7]+f[0]*alphaDrag[7])+0.5477225575051661*alphaDrag[1]*f[1]; 
  out[12] += (4.166190448976479*facDiff[2]*f[17]+4.242640687119286*(facDiff[1]*f[7]+f[1]*facDiff[2])+4.166190448976479*facDiff[3]*f[7]+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4+1.202675588605909*(alphaDrag[7]*f[23]+f[11]*alphaDrag[17])+1.224744871391589*(alphaDrag[2]*f[12]+alphaDrag[1]*f[11]+f[4]*alphaDrag[7])+1.369306393762915*(alphaDrag[0]*f[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[13] += (0.3651483716701108*f[17]+0.537852874200477*f[1])*alphaDrag[49]+0.3912303982179757*f[7]*alphaDrag[39]+0.6123724356957944*(f[0]*alphaDrag[39]+f[13]*alphaDrag[35])+(0.537852874200477*f[17]+0.5477225575051661*f[1])*alphaDrag[33]+0.6123724356957944*f[7]*alphaDrag[32]; 
  out[14] += 4.743416490252569*(facDiff[3]*f[25]+facDiff[2]*f[13]+facDiff[1]*f[5]+facDiff[0]*f[3])*rdvxSq4+0.6123724356957944*(f[14]*alphaDrag[35]+f[12]*alphaDrag[33]+f[8]*alphaDrag[32])+1.369306393762915*(alphaDrag[17]*f[29]+alphaDrag[7]*f[20])+1.224744871391589*alphaDrag[2]*f[14]+1.369306393762915*(alphaDrag[1]*f[10]+alphaDrag[0]*f[6]+alphaDrag[2]*f[3]); 
  out[15] += (4.166190448976479*facDiff[2]*f[17]+4.242640687119286*(facDiff[1]*f[7]+f[1]*facDiff[2])+4.166190448976479*facDiff[3]*f[7]+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvySq4+1.202675588605909*f[13]*alphaDrag[49]+(1.202675588605909*f[25]+1.224744871391589*f[5])*alphaDrag[39]+(1.224744871391589*f[15]+1.369306393762915*f[1])*alphaDrag[35]+1.224744871391589*f[13]*alphaDrag[33]+1.369306393762915*(f[3]*alphaDrag[33]+f[5]*alphaDrag[32]); 
  out[16] += 4.743416490252569*(facDiff[3]*f[23]+facDiff[2]*f[11]+facDiff[1]*f[4]+facDiff[0]*f[2])*rdvySq4+1.369306393762915*(f[29]*alphaDrag[49]+f[20]*alphaDrag[39])+1.224744871391589*f[16]*alphaDrag[35]+1.369306393762915*(f[2]*alphaDrag[35]+f[10]*alphaDrag[33]+f[6]*alphaDrag[32])+0.6123724356957944*(alphaDrag[2]*f[16]+alphaDrag[1]*f[15]+alphaDrag[0]*f[9]); 
  out[18] += 16.20185174601965*(facDiff[3]*f[23]+facDiff[2]*f[11]+facDiff[1]*f[4]+facDiff[0]*f[2])*rdvxSq4+1.837117307087383*alphaDrag[2]*f[18]+0.9354143466934851*alphaDrag[17]*f[17]+2.091650066335188*(alphaDrag[1]*f[12]+alphaDrag[0]*f[8])+0.9354143466934851*alphaDrag[7]*f[7]+2.806243040080455*alphaDrag[2]*f[2]+0.9354143466934851*(alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[19] += 16.20185174601965*(facDiff[3]*f[25]+facDiff[2]*f[13]+facDiff[1]*f[5]+facDiff[0]*f[3])*rdvySq4+0.9354143466934851*(f[17]*alphaDrag[49]+f[7]*alphaDrag[39])+(1.837117307087383*f[19]+2.806243040080455*f[3])*alphaDrag[35]+(2.091650066335188*f[15]+0.9354143466934851*f[1])*alphaDrag[33]+(2.091650066335188*f[9]+0.9354143466934851*f[0])*alphaDrag[32]; 
  out[20] += (0.3651483716701108*f[23]+0.537852874200477*f[4])*alphaDrag[49]+0.3912303982179757*f[11]*alphaDrag[39]+0.6123724356957944*(f[2]*alphaDrag[39]+f[20]*alphaDrag[35])+(0.537852874200477*f[23]+0.5477225575051661*f[4])*alphaDrag[33]+0.6123724356957944*f[11]*alphaDrag[32]+(0.3651483716701108*alphaDrag[17]+0.537852874200477*alphaDrag[1])*f[25]+0.6123724356957944*alphaDrag[2]*f[20]+0.537852874200477*f[5]*alphaDrag[17]+0.3912303982179757*alphaDrag[7]*f[13]+0.6123724356957944*(alphaDrag[0]*f[13]+f[3]*alphaDrag[7])+0.5477225575051661*alphaDrag[1]*f[5]; 
  out[21] += (4.166190448976479*facDiff[2]*f[25]+(4.166190448976479*facDiff[3]+4.242640687119286*facDiff[1])*f[13]+4.743416490252569*(facDiff[0]*f[5]+facDiff[1]*f[3])+4.242640687119286*facDiff[2]*f[5])*rdvxSq4+0.5477225575051661*f[12]*alphaDrag[39]+0.6123724356957944*(f[21]*alphaDrag[35]+f[8]*alphaDrag[33]+f[12]*alphaDrag[32])+1.202675588605909*alphaDrag[7]*f[29]+1.224744871391589*alphaDrag[2]*f[21]+1.202675588605909*alphaDrag[17]*f[20]+1.224744871391589*(alphaDrag[1]*f[20]+alphaDrag[7]*f[10])+1.369306393762915*(alphaDrag[0]*f[10]+alphaDrag[1]*f[6]+alphaDrag[2]*f[5]); 
  out[22] += (4.166190448976479*facDiff[2]*f[23]+(4.166190448976479*facDiff[3]+4.242640687119286*facDiff[1])*f[11]+4.743416490252569*(facDiff[0]*f[4]+facDiff[1]*f[2])+4.242640687119286*facDiff[2]*f[4])*rdvySq4+1.202675588605909*f[20]*alphaDrag[49]+(1.202675588605909*f[29]+1.224744871391589*f[10])*alphaDrag[39]+(1.224744871391589*f[22]+1.369306393762915*f[4])*alphaDrag[35]+1.224744871391589*f[20]*alphaDrag[33]+1.369306393762915*(f[6]*alphaDrag[33]+f[10]*alphaDrag[32])+0.6123724356957944*alphaDrag[2]*f[22]+0.5477225575051661*alphaDrag[7]*f[15]+0.6123724356957944*(alphaDrag[0]*f[15]+alphaDrag[1]*f[9]); 
  out[23] += 0.6123724356957944*alphaDrag[2]*f[23]+(0.3651483716701108*alphaDrag[7]+0.6123724356957944*alphaDrag[0])*f[17]+(0.3651483716701108*f[7]+0.6123724356957944*f[0])*alphaDrag[17]+0.537852874200477*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7]); 
  out[24] += (14.23024947075771*facDiff[2]*f[23]+(14.23024947075771*facDiff[3]+14.49137674618944*facDiff[1])*f[11]+16.20185174601965*(facDiff[0]*f[4]+facDiff[1]*f[2])+14.49137674618944*facDiff[2]*f[4])*rdvxSq4+1.837117307087383*alphaDrag[2]*f[24]+0.8215838362577489*(alphaDrag[7]*f[17]+f[7]*alphaDrag[17])+1.870828693386971*alphaDrag[7]*f[12]+2.091650066335188*(alphaDrag[0]*f[12]+alphaDrag[1]*f[8])+0.8366600265340755*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+2.806243040080455*alphaDrag[2]*f[4]+0.9354143466934851*(alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[25] += (0.3651483716701108*f[7]+0.6123724356957944*f[0])*alphaDrag[49]+(0.3651483716701108*f[17]+0.537852874200477*f[1])*alphaDrag[39]+0.6123724356957944*f[25]*alphaDrag[35]+0.537852874200477*f[7]*alphaDrag[33]+0.6123724356957944*f[17]*alphaDrag[32]; 
  out[26] += 16.20185174601965*(facDiff[3]*f[29]+facDiff[2]*f[20]+facDiff[1]*f[10]+facDiff[0]*f[6])*rdvxSq4+0.6123724356957944*(f[26]*alphaDrag[35]+f[24]*alphaDrag[33]+f[18]*alphaDrag[32])+1.837117307087383*alphaDrag[2]*f[26]+0.9354143466934851*alphaDrag[17]*f[25]+2.091650066335188*(alphaDrag[1]*f[21]+alphaDrag[0]*f[14])+0.9354143466934851*alphaDrag[7]*f[13]+2.806243040080455*alphaDrag[2]*f[6]+0.9354143466934851*(alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[27] += (14.23024947075771*facDiff[2]*f[25]+(14.23024947075771*facDiff[3]+14.49137674618944*facDiff[1])*f[13]+16.20185174601965*(facDiff[0]*f[5]+facDiff[1]*f[3])+14.49137674618944*facDiff[2]*f[5])*rdvySq4+0.8215838362577489*f[7]*alphaDrag[49]+(0.8215838362577489*f[17]+1.870828693386971*f[15]+0.8366600265340755*f[1])*alphaDrag[39]+(1.837117307087383*f[27]+2.806243040080455*f[5])*alphaDrag[35]+(2.091650066335188*f[9]+0.8366600265340755*f[7]+0.9354143466934851*f[0])*alphaDrag[33]+(2.091650066335188*f[15]+0.9354143466934851*f[1])*alphaDrag[32]; 
  out[28] += 16.20185174601965*(facDiff[3]*f[29]+facDiff[2]*f[20]+facDiff[1]*f[10]+facDiff[0]*f[6])*rdvySq4+0.9354143466934851*(f[23]*alphaDrag[49]+f[11]*alphaDrag[39])+(1.837117307087383*f[28]+2.806243040080455*f[6])*alphaDrag[35]+(2.091650066335188*f[22]+0.9354143466934851*f[4])*alphaDrag[33]+(2.091650066335188*f[16]+0.9354143466934851*f[2])*alphaDrag[32]+0.6123724356957944*(alphaDrag[2]*f[28]+alphaDrag[1]*f[27]+alphaDrag[0]*f[19]); 
  out[29] += (0.3651483716701108*f[11]+0.6123724356957944*f[2])*alphaDrag[49]+(0.3651483716701108*f[23]+0.537852874200477*f[4])*alphaDrag[39]+0.6123724356957944*f[29]*alphaDrag[35]+0.537852874200477*f[11]*alphaDrag[33]+0.6123724356957944*(f[23]*alphaDrag[32]+alphaDrag[2]*f[29])+(0.3651483716701108*alphaDrag[7]+0.6123724356957944*alphaDrag[0])*f[25]+(0.3651483716701108*f[13]+0.6123724356957944*f[3])*alphaDrag[17]+0.537852874200477*(alphaDrag[1]*f[13]+f[5]*alphaDrag[7]); 
  out[30] += (14.23024947075771*facDiff[2]*f[29]+(14.23024947075771*facDiff[3]+14.49137674618944*facDiff[1])*f[20]+16.20185174601965*(facDiff[0]*f[10]+facDiff[1]*f[6])+14.49137674618944*facDiff[2]*f[10])*rdvxSq4+0.5477225575051661*f[24]*alphaDrag[39]+0.6123724356957944*(f[30]*alphaDrag[35]+f[18]*alphaDrag[33]+f[24]*alphaDrag[32])+1.837117307087383*alphaDrag[2]*f[30]+0.8215838362577489*alphaDrag[7]*f[25]+(1.870828693386971*alphaDrag[7]+2.091650066335188*alphaDrag[0])*f[21]+0.8215838362577489*f[13]*alphaDrag[17]+alphaDrag[1]*(2.091650066335188*f[14]+0.8366600265340755*f[13])+2.806243040080455*alphaDrag[2]*f[10]+0.8366600265340755*f[5]*alphaDrag[7]+0.9354143466934851*(alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[31] += (14.23024947075771*facDiff[2]*f[29]+(14.23024947075771*facDiff[3]+14.49137674618944*facDiff[1])*f[20]+16.20185174601965*(facDiff[0]*f[10]+facDiff[1]*f[6])+14.49137674618944*facDiff[2]*f[10])*rdvySq4+0.8215838362577489*f[11]*alphaDrag[49]+(0.8215838362577489*f[23]+1.870828693386971*f[22]+0.8366600265340755*f[4])*alphaDrag[39]+(1.837117307087383*f[31]+2.806243040080455*f[10])*alphaDrag[35]+(2.091650066335188*f[16]+0.8366600265340755*f[11]+0.9354143466934851*f[2])*alphaDrag[33]+(2.091650066335188*f[22]+0.9354143466934851*f[4])*alphaDrag[32]+0.6123724356957944*alphaDrag[2]*f[31]+0.5477225575051661*alphaDrag[7]*f[27]+0.6123724356957944*(alphaDrag[0]*f[27]+alphaDrag[1]*f[19]); 

  return std::abs(0.1767766952966368*alphaDrag[0]-0.1976423537605236*alphaDrag[7])+std::abs(0.1767766952966368*alphaDrag[32]-0.1976423537605236*alphaDrag[39])+std::abs((1.616244071283538*facDiff[0]-1.807015805810503*facDiff[2])*rdvxSq4)+std::abs((1.616244071283538*facDiff[0]-1.807015805810503*facDiff[2])*rdvySq4); 

} 
