#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionVol3xSerP1(const double *w, const double *dxv, double *positivityWeightByDir, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing.
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  double dfac3 = 2.0/dxv[2]; 
  double w3 = w[2]; 
  const double *v1 = &f[8]; 
  const double *v2 = &f[16]; 
  const double *v3 = &f[24]; 
  double cflRate = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  positivityWeightByDir[1] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(2.449489742783178*v1[1]-1.414213562373095*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(2.449489742783178*v1[1]+1.414213562373095*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*((-1.837117307087383*v1[7])+1.060660171779821*(v1[6]+v1[5]+v1[4])-0.6123724356957944*(v1[3]+v1[2]+v1[1])+0.3535533905932737*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.837117307087383*v1[7]-1.060660171779821*v1[6]+1.060660171779821*v1[5]-1.060660171779821*v1[4]-0.6123724356957944*v1[3]+0.6123724356957944*v1[2]-0.6123724356957944*v1[1]+0.3535533905932737*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.837117307087383*v1[7]-1.060660171779821*(v1[6]+v1[5])+1.060660171779821*v1[4]+0.6123724356957944*v1[3]-0.6123724356957944*(v1[2]+v1[1])+0.3535533905932737*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-1.837117307087383*v1[7])+1.060660171779821*v1[6]-1.060660171779821*(v1[5]+v1[4])+0.6123724356957944*(v1[3]+v1[2])-0.6123724356957944*v1[1]+0.3535533905932737*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(1.837117307087383*v1[7]+1.060660171779821*v1[6]-1.060660171779821*(v1[5]+v1[4])-0.6123724356957944*(v1[3]+v1[2])+0.6123724356957944*v1[1]+0.3535533905932737*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-1.837117307087383*v1[7])-1.060660171779821*(v1[6]+v1[5])+1.060660171779821*v1[4]-0.6123724356957944*v1[3]+0.6123724356957944*(v1[2]+v1[1])+0.3535533905932737*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-1.837117307087383*v1[7])-1.060660171779821*v1[6]+1.060660171779821*v1[5]-1.060660171779821*v1[4]+0.6123724356957944*v1[3]-0.6123724356957944*v1[2]+0.6123724356957944*v1[1]+0.3535533905932737*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(1.837117307087383*v1[7]+1.060660171779821*(v1[6]+v1[5]+v1[4])+0.6123724356957944*(v1[3]+v1[2]+v1[1])+0.3535533905932737*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  positivityWeightByDir[2] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(2.449489742783178*v2[2]-1.414213562373095*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(2.449489742783178*v2[2]+1.414213562373095*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*((-1.837117307087383*v2[7])+1.060660171779821*(v2[6]+v2[5]+v2[4])-0.6123724356957944*(v2[3]+v2[2]+v2[1])+0.3535533905932737*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.837117307087383*v2[7]+1.060660171779821*v2[6]-1.060660171779821*(v2[5]+v2[4])-0.6123724356957944*(v2[3]+v2[2])+0.6123724356957944*v2[1]+0.3535533905932737*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.837117307087383*v2[7]-1.060660171779821*(v2[6]+v2[5])+1.060660171779821*v2[4]+0.6123724356957944*v2[3]-0.6123724356957944*(v2[2]+v2[1])+0.3535533905932737*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-1.837117307087383*v2[7])-1.060660171779821*v2[6]+1.060660171779821*v2[5]-1.060660171779821*v2[4]+0.6123724356957944*v2[3]-0.6123724356957944*v2[2]+0.6123724356957944*v2[1]+0.3535533905932737*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(1.837117307087383*v2[7]-1.060660171779821*v2[6]+1.060660171779821*v2[5]-1.060660171779821*v2[4]-0.6123724356957944*v2[3]+0.6123724356957944*v2[2]-0.6123724356957944*v2[1]+0.3535533905932737*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-1.837117307087383*v2[7])-1.060660171779821*(v2[6]+v2[5])+1.060660171779821*v2[4]-0.6123724356957944*v2[3]+0.6123724356957944*(v2[2]+v2[1])+0.3535533905932737*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-1.837117307087383*v2[7])+1.060660171779821*v2[6]-1.060660171779821*(v2[5]+v2[4])+0.6123724356957944*(v2[3]+v2[2])-0.6123724356957944*v2[1]+0.3535533905932737*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(1.837117307087383*v2[7]+1.060660171779821*(v2[6]+v2[5]+v2[4])+0.6123724356957944*(v2[3]+v2[2]+v2[1])+0.3535533905932737*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  positivityWeightByDir[3] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(2.449489742783178*v3[3]-1.414213562373095*v3[0])*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[3] += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(2.449489742783178*v3[3]+1.414213562373095*v3[0])*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[3] += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*((-1.837117307087383*v3[7])+1.060660171779821*(v3[6]+v3[5]+v3[4])-0.6123724356957944*(v3[3]+v3[2]+v3[1])+0.3535533905932737*v3[0])*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[3] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.837117307087383*v3[7]+1.060660171779821*v3[6]-1.060660171779821*(v3[5]+v3[4])-0.6123724356957944*(v3[3]+v3[2])+0.6123724356957944*v3[1]+0.3535533905932737*v3[0])*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[3] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.837117307087383*v3[7]-1.060660171779821*v3[6]+1.060660171779821*v3[5]-1.060660171779821*v3[4]-0.6123724356957944*v3[3]+0.6123724356957944*v3[2]-0.6123724356957944*v3[1]+0.3535533905932737*v3[0])*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[3] += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-1.837117307087383*v3[7])-1.060660171779821*(v3[6]+v3[5])+1.060660171779821*v3[4]-0.6123724356957944*v3[3]+0.6123724356957944*(v3[2]+v3[1])+0.3535533905932737*v3[0])*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[3] += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(1.837117307087383*v3[7]-1.060660171779821*(v3[6]+v3[5])+1.060660171779821*v3[4]+0.6123724356957944*v3[3]-0.6123724356957944*(v3[2]+v3[1])+0.3535533905932737*v3[0])*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[3] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-1.837117307087383*v3[7])-1.060660171779821*v3[6]+1.060660171779821*v3[5]-1.060660171779821*v3[4]+0.6123724356957944*v3[3]-0.6123724356957944*v3[2]+0.6123724356957944*v3[1]+0.3535533905932737*v3[0])*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[3] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-1.837117307087383*v3[7])+1.060660171779821*v3[6]-1.060660171779821*(v3[5]+v3[4])+0.6123724356957944*(v3[3]+v3[2])-0.6123724356957944*v3[1]+0.3535533905932737*v3[0])*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[3] += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(1.837117307087383*v3[7]+1.060660171779821*(v3[6]+v3[5]+v3[4])+0.6123724356957944*(v3[3]+v3[2]+v3[1])+0.3535533905932737*v3[0])*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[3] += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*(f[7]*v1[7]+f[6]*v1[6]+f[5]*v1[5]+f[4]*v1[4]+f[3]*v1[3]+f[2]*v1[2]+f[1]*v1[1]+f[0]*v1[0])*dfac1; 
  out[2] += 0.6123724356957944*(f[7]*v2[7]+f[6]*v2[6]+f[5]*v2[5]+f[4]*v2[4]+f[3]*v2[3]+f[2]*v2[2]+f[1]*v2[1]+f[0]*v2[0])*dfac2; 
  out[3] += 0.6123724356957944*(f[7]*v3[7]+f[6]*v3[6]+f[5]*v3[5]+f[4]*v3[4]+f[3]*v3[3]+f[2]*v3[2]+f[1]*v3[1]+f[0]*v3[0])*dfac3; 
  out[4] += 0.6123724356957946*((f[6]*v2[7]+v2[6]*f[7]+f[3]*v2[5]+v2[3]*f[5]+f[2]*v2[4]+v2[2]*f[4]+f[0]*v2[1]+v2[0]*f[1])*dfac2+(f[5]*v1[7]+v1[5]*f[7]+f[3]*v1[6]+v1[3]*f[6]+f[1]*v1[4]+v1[1]*f[4]+f[0]*v1[2]+v1[0]*f[2])*dfac1); 
  out[5] += 0.6123724356957946*((f[6]*v3[7]+v3[6]*f[7]+f[3]*v3[5]+v3[3]*f[5]+f[2]*v3[4]+v3[2]*f[4]+f[0]*v3[1]+v3[0]*f[1])*dfac3+(f[4]*v1[7]+v1[4]*f[7]+f[2]*v1[6]+v1[2]*f[6]+f[1]*v1[5]+v1[1]*f[5]+f[0]*v1[3]+v1[0]*f[3])*dfac1); 
  out[6] += 0.6123724356957946*((f[5]*v3[7]+v3[5]*f[7]+f[3]*v3[6]+v3[3]*f[6]+f[1]*v3[4]+v3[1]*f[4]+f[0]*v3[2]+v3[0]*f[2])*dfac3+(f[4]*v2[7]+v2[4]*f[7]+f[2]*v2[6]+v2[2]*f[6]+f[1]*v2[5]+v2[1]*f[5]+f[0]*v2[3]+v2[0]*f[3])*dfac2); 
  out[7] += 0.6123724356957946*((f[3]*v3[7]+v3[3]*f[7]+f[5]*v3[6]+v3[5]*f[6]+f[0]*v3[4]+v3[0]*f[4]+f[1]*v3[2]+v3[1]*f[2])*dfac3+(f[2]*v2[7]+v2[2]*f[7]+f[4]*v2[6]+v2[4]*f[6]+f[0]*v2[5]+v2[0]*f[5]+f[1]*v2[3]+v2[1]*f[3])*dfac2+(f[1]*v1[7]+v1[1]*f[7]+f[0]*v1[6]+v1[0]*f[6]+f[4]*v1[5]+v1[4]*f[5]+f[2]*v1[3]+v1[2]*f[3])*dfac1); 
  return cflRate; 
} 
