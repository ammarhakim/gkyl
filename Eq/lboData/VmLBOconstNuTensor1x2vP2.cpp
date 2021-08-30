#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x2vTensorP2(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
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

  double alphaDrag[54]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.828427124746191*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.0*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -0.8164965809277261*dxv[1]*nuSum*rdvx2; 
  alphaDrag[7] = 2.0*nuUSum[2]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[27] = (2.0*nuUSum[3]-2.828427124746191*w[2]*nuSum)*rdvy2; 
  alphaDrag[28] = 2.0*nuUSum[4]*rdvy2; 
  alphaDrag[30] = -0.8164965809277261*dxv[2]*nuSum*rdvy2; 
  alphaDrag[34] = 2.0*nuUSum[5]*rdvy2; 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alphaDrag[7]*f[7]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*alphaDrag[34]+f[3]*alphaDrag[30]+f[1]*alphaDrag[28]+f[0]*alphaDrag[27]); 
  out[4] += 0.5477225575051661*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.6123724356957944*(alphaDrag[2]*f[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.5477225575051661*f[1]*alphaDrag[34]+0.6123724356957944*f[5]*alphaDrag[30]+0.5477225575051661*f[7]*alphaDrag[28]+0.6123724356957944*(f[0]*alphaDrag[28]+f[1]*alphaDrag[27]); 
  out[6] += 0.6123724356957944*(f[11]*alphaDrag[34]+f[6]*alphaDrag[30]+f[4]*alphaDrag[28]+f[2]*alphaDrag[27]+alphaDrag[7]*f[13]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[8] += 4.743416490252569*(facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+1.369306393762915*alphaDrag[7]*f[11]+1.224744871391589*alphaDrag[2]*f[8]+1.369306393762915*(alphaDrag[1]*f[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 4.743416490252569*(facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+1.369306393762915*f[13]*alphaDrag[34]+1.224744871391589*f[9]*alphaDrag[30]+1.369306393762915*(f[0]*alphaDrag[30]+f[5]*alphaDrag[28]+f[3]*alphaDrag[27]); 
  out[10] += 0.5477225575051661*f[4]*alphaDrag[34]+0.6123724356957944*f[10]*alphaDrag[30]+0.5477225575051661*f[11]*alphaDrag[28]+0.6123724356957944*(f[2]*alphaDrag[28]+f[4]*alphaDrag[27])+0.5477225575051661*alphaDrag[1]*f[13]+0.6123724356957944*alphaDrag[2]*f[10]+0.5477225575051661*f[5]*alphaDrag[7]+0.6123724356957944*(alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[11] += 0.6123724356957944*alphaDrag[2]*f[11]+0.3912303982179757*alphaDrag[7]*f[7]+0.6123724356957944*(alphaDrag[0]*f[7]+f[0]*alphaDrag[7])+0.5477225575051661*alphaDrag[1]*f[1]; 
  out[12] += (4.242640687119286*(facDiff[1]*f[7]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4+1.224744871391589*(alphaDrag[2]*f[12]+alphaDrag[1]*f[11]+f[4]*alphaDrag[7])+1.369306393762915*(alphaDrag[0]*f[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[13] += 0.3912303982179757*f[7]*alphaDrag[34]+0.6123724356957944*(f[0]*alphaDrag[34]+f[13]*alphaDrag[30])+0.5477225575051661*f[1]*alphaDrag[28]+0.6123724356957944*f[7]*alphaDrag[27]; 
  out[14] += 4.743416490252569*(facDiff[2]*f[13]+facDiff[1]*f[5]+facDiff[0]*f[3])*rdvxSq4+0.6123724356957944*(f[20]*alphaDrag[34]+f[14]*alphaDrag[30]+f[12]*alphaDrag[28]+f[8]*alphaDrag[27])+1.369306393762915*alphaDrag[7]*f[17]+1.224744871391589*alphaDrag[2]*f[14]+1.369306393762915*(alphaDrag[1]*f[10]+alphaDrag[0]*f[6]+alphaDrag[2]*f[3]); 
  out[15] += (4.242640687119286*(facDiff[1]*f[7]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvySq4+1.224744871391589*f[5]*alphaDrag[34]+(1.224744871391589*f[15]+1.369306393762915*f[1])*alphaDrag[30]+1.224744871391589*f[13]*alphaDrag[28]+1.369306393762915*(f[3]*alphaDrag[28]+f[5]*alphaDrag[27]); 
  out[16] += 4.743416490252569*(facDiff[2]*f[11]+facDiff[1]*f[4]+facDiff[0]*f[2])*rdvySq4+1.369306393762915*f[17]*alphaDrag[34]+1.224744871391589*f[16]*alphaDrag[30]+1.369306393762915*(f[2]*alphaDrag[30]+f[10]*alphaDrag[28]+f[6]*alphaDrag[27])+0.6123724356957944*(alphaDrag[7]*f[21]+alphaDrag[2]*f[16]+alphaDrag[1]*f[15]+alphaDrag[0]*f[9]); 
  out[17] += 0.3912303982179757*f[11]*alphaDrag[34]+0.6123724356957944*(f[2]*alphaDrag[34]+f[17]*alphaDrag[30])+0.5477225575051661*f[4]*alphaDrag[28]+0.6123724356957944*(f[11]*alphaDrag[27]+alphaDrag[2]*f[17])+0.3912303982179757*alphaDrag[7]*f[13]+0.6123724356957944*(alphaDrag[0]*f[13]+f[3]*alphaDrag[7])+0.5477225575051661*alphaDrag[1]*f[5]; 
  out[18] += (4.242640687119286*facDiff[1]*f[13]+4.743416490252569*(facDiff[0]*f[5]+facDiff[1]*f[3])+4.242640687119286*facDiff[2]*f[5])*rdvxSq4+0.5477225575051661*f[12]*alphaDrag[34]+0.6123724356957944*f[18]*alphaDrag[30]+0.5477225575051661*f[20]*alphaDrag[28]+0.6123724356957944*(f[8]*alphaDrag[28]+f[12]*alphaDrag[27])+1.224744871391589*(alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[7]*f[10])+1.369306393762915*(alphaDrag[0]*f[10]+alphaDrag[1]*f[6]+alphaDrag[2]*f[5]); 
  out[19] += (4.242640687119286*facDiff[1]*f[11]+4.743416490252569*(facDiff[0]*f[4]+facDiff[1]*f[2])+4.242640687119286*facDiff[2]*f[4])*rdvySq4+1.224744871391589*f[10]*alphaDrag[34]+(1.224744871391589*f[19]+1.369306393762915*f[4])*alphaDrag[30]+1.224744871391589*f[17]*alphaDrag[28]+1.369306393762915*(f[6]*alphaDrag[28]+f[10]*alphaDrag[27])+0.5477225575051661*alphaDrag[1]*f[21]+0.6123724356957944*alphaDrag[2]*f[19]+0.5477225575051661*alphaDrag[7]*f[15]+0.6123724356957944*(alphaDrag[0]*f[15]+alphaDrag[1]*f[9]); 
  out[20] += (4.743416490252569*(facDiff[0]*f[7]+f[0]*facDiff[2])+3.030457633656632*facDiff[2]*f[7]+4.242640687119286*f[1]*facDiff[1])*rdvxSq4+1.224744871391589*alphaDrag[2]*f[20]+0.8748177652797062*alphaDrag[7]*f[11]+1.369306393762915*(alphaDrag[0]*f[11]+alphaDrag[2]*f[7]+f[2]*alphaDrag[7])+1.224744871391589*alphaDrag[1]*f[4]; 
  out[21] += (4.743416490252569*(facDiff[0]*f[7]+f[0]*facDiff[2])+3.030457633656632*facDiff[2]*f[7]+4.242640687119286*f[1]*facDiff[1])*rdvySq4+(0.8748177652797062*f[13]+1.369306393762915*f[3])*alphaDrag[34]+(1.224744871391589*f[21]+1.369306393762915*f[7])*alphaDrag[30]+1.224744871391589*f[5]*alphaDrag[28]+1.369306393762915*f[13]*alphaDrag[27]; 
  out[22] += 4.743416490252569*((facDiff[2]*f[20]+facDiff[1]*f[12]+facDiff[0]*f[8])*rdvySq4+(facDiff[2]*f[21]+facDiff[1]*f[15]+facDiff[0]*f[9])*rdvxSq4)+1.369306393762915*f[23]*alphaDrag[34]+1.224744871391589*f[22]*alphaDrag[30]+1.369306393762915*(f[8]*alphaDrag[30]+f[18]*alphaDrag[28]+f[14]*alphaDrag[27]+alphaDrag[7]*f[24])+1.224744871391589*alphaDrag[2]*f[22]+1.369306393762915*(alphaDrag[1]*f[19]+alphaDrag[0]*f[16]+alphaDrag[2]*f[9]); 
  out[23] += ((3.030457633656632*facDiff[2]+4.743416490252569*facDiff[0])*f[13]+4.242640687119286*facDiff[1]*f[5]+4.743416490252569*facDiff[2]*f[3])*rdvxSq4+0.3912303982179757*f[20]*alphaDrag[34]+0.6123724356957944*(f[8]*alphaDrag[34]+f[23]*alphaDrag[30])+0.5477225575051661*f[12]*alphaDrag[28]+0.6123724356957944*f[20]*alphaDrag[27]+1.224744871391589*alphaDrag[2]*f[23]+0.8748177652797062*alphaDrag[7]*f[17]+1.369306393762915*(alphaDrag[0]*f[17]+alphaDrag[2]*f[13])+1.224744871391589*alphaDrag[1]*f[10]+1.369306393762915*f[6]*alphaDrag[7]; 
  out[24] += ((3.030457633656632*facDiff[2]+4.743416490252569*facDiff[0])*f[11]+4.242640687119286*facDiff[1]*f[4]+4.743416490252569*f[2]*facDiff[2])*rdvySq4+(0.8748177652797062*f[17]+1.369306393762915*f[6])*alphaDrag[34]+(1.224744871391589*f[24]+1.369306393762915*f[11])*alphaDrag[30]+1.224744871391589*f[10]*alphaDrag[28]+1.369306393762915*f[17]*alphaDrag[27]+0.6123724356957944*alphaDrag[2]*f[24]+(0.3912303982179757*alphaDrag[7]+0.6123724356957944*alphaDrag[0])*f[21]+0.5477225575051661*alphaDrag[1]*f[15]+0.6123724356957944*alphaDrag[7]*f[9]; 
  out[25] += (4.242640687119286*facDiff[1]*f[20]+4.743416490252569*(facDiff[0]*f[12]+facDiff[1]*f[8])+4.242640687119286*facDiff[2]*f[12])*rdvySq4+(4.242640687119286*facDiff[1]*f[21]+4.743416490252569*(facDiff[0]*f[15]+facDiff[1]*f[9])+4.242640687119286*facDiff[2]*f[15])*rdvxSq4+1.224744871391589*f[18]*alphaDrag[34]+(1.224744871391589*f[25]+1.369306393762915*f[12])*alphaDrag[30]+1.224744871391589*f[23]*alphaDrag[28]+1.369306393762915*(f[14]*alphaDrag[28]+f[18]*alphaDrag[27])+1.224744871391589*(alphaDrag[2]*f[25]+alphaDrag[1]*f[24]+alphaDrag[7]*f[19])+1.369306393762915*(alphaDrag[0]*f[19]+alphaDrag[1]*f[16]+alphaDrag[2]*f[15]); 
  out[26] += ((3.030457633656632*facDiff[2]+4.743416490252569*facDiff[0])*f[20]+4.242640687119286*facDiff[1]*f[12]+4.743416490252569*facDiff[2]*f[8])*rdvySq4+((3.030457633656632*facDiff[2]+4.743416490252569*facDiff[0])*f[21]+4.242640687119286*facDiff[1]*f[15]+4.743416490252569*facDiff[2]*f[9])*rdvxSq4+(0.8748177652797062*f[23]+1.369306393762915*f[14])*alphaDrag[34]+(1.224744871391589*f[26]+1.369306393762915*f[20])*alphaDrag[30]+1.224744871391589*f[18]*alphaDrag[28]+1.369306393762915*f[23]*alphaDrag[27]+1.224744871391589*alphaDrag[2]*f[26]+0.8748177652797062*alphaDrag[7]*f[24]+1.369306393762915*(alphaDrag[0]*f[24]+alphaDrag[2]*f[21])+1.224744871391589*alphaDrag[1]*f[19]+1.369306393762915*alphaDrag[7]*f[16]; 

  return std::abs(0.1767766952966368*alphaDrag[0]-0.1976423537605236*alphaDrag[7])+std::abs(0.1767766952966368*alphaDrag[27]-0.1976423537605236*alphaDrag[34])+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvxSq4)+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvySq4); 

} 
