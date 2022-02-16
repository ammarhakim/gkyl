#include <VlasovModDecl.h> 
__host__ __device__ double VlasovUpwindSurfNeutral1x2vSer_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double *Fo0 = &EM[0]; 

  double alpha[4]; 
  double incr[8]; 

  alpha[0] = 1.414213562373095*Fo0[0]; 
  alpha[1] = 1.414213562373095*Fo0[1]; 

  const double amid = 0.5*alpha[0]; 

  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alpha[0]-0.5*alpha[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaOrdR)-0.6123724356957944*fr[7]+0.6123724356957944*(fl[7]+fr[6])-0.6123724356957944*fl[6]+0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*(alpha[1]+alpha[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaOrdR)+0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6])-0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaOrdR = 0.5*alpha[0]-0.5*alpha[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaOrdR)+0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*fl[6]-0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*(alpha[1]+alpha[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaOrdR)-0.6123724356957944*(fr[7]+fr[6])+0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[1] = 0.3535533905932737*(alpha[0]*fUp[1]+fUp[0]*alpha[1]); 
  incr[2] = -0.6123724356957944*(alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[3] = 0.3535533905932737*(alpha[1]*fUp[3]+alpha[0]*fUp[2]); 
  incr[4] = -0.6123724356957944*(alpha[0]*fUp[1]+fUp[0]*alpha[1]); 
  incr[5] = 0.3535533905932737*(alpha[0]*fUp[3]+alpha[1]*fUp[2]); 
  incr[6] = -0.6123724356957944*(alpha[1]*fUp[3]+alpha[0]*fUp[2]); 
  incr[7] = -0.6123724356957944*(alpha[0]*fUp[3]+alpha[1]*fUp[2]); 

  outr[0] += incr[0]*dv10r; 
  outr[1] += incr[1]*dv10r; 
  outr[2] += incr[2]*dv10r; 
  outr[3] += incr[3]*dv10r; 
  outr[4] += incr[4]*dv10r; 
  outr[5] += incr[5]*dv10r; 
  outr[6] += incr[6]*dv10r; 
  outr[7] += incr[7]*dv10r; 

  outl[0] += -1.0*incr[0]*dv10l; 
  outl[1] += -1.0*incr[1]*dv10l; 
  outl[2] += incr[2]*dv10l; 
  outl[3] += -1.0*incr[3]*dv10l; 
  outl[4] += incr[4]*dv10l; 
  outl[5] += -1.0*incr[5]*dv10l; 
  outl[6] += incr[6]*dv10l; 
  outl[7] += incr[7]*dv10l; 


  return std::abs(amid); 
} 
__host__ __device__ double VlasovUpwindSurfNeutral1x2vSer_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // EM:        EM field.
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: output distribution function in left/right cells 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[2]; 
  double dv11r = 2/dxvr[2]; 

  const double dv1 = dxvr[1], wv1 = wr[1]; 
  const double dv2 = dxvr[2], wv2 = wr[2]; 
  const double *Fo1 = &EM[2]; 

  double alpha[4]; 
  double incr[8]; 

  alpha[0] = 1.414213562373095*Fo1[0]; 
  alpha[1] = 1.414213562373095*Fo1[1]; 

  const double amid = 0.5*alpha[0]; 

  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alpha[0]-0.5*alpha[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6]+fr[5]+fl[5])-0.3535533905932737*fr[4]+0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*(fr[2]+fr[1])-0.3535533905932737*(fl[2]+fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaOrdR)-0.6123724356957944*fr[7]+0.6123724356957944*(fl[7]+fr[6]+fr[5])-0.6123724356957944*(fl[6]+fl[5])+0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2]+fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*(alpha[1]+alpha[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.6123724356957944*(fr[5]+fl[5])+0.3535533905932737*fr[4]-0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*fr[2]-0.3535533905932737*(fl[2]+fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaOrdR)+0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6]+fr[5])+0.6123724356957944*fl[5]-0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2])+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0])); 
  alphaOrdR = 0.5*alpha[0]-0.5*alpha[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])-0.6123724356957944*(fr[5]+fl[5])+0.3535533905932737*fr[4]-0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*fr[2]+0.3535533905932737*(fl[2]+fr[1])-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaOrdR)+0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*(fl[6]+fr[5])-0.6123724356957944*fl[5]-0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0])); 
  alphaOrdR = 0.5*(alpha[1]+alpha[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.3535533905932737*fr[4]+0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*(fr[2]+fr[1]+fr[0])+0.3535533905932737*(fl[2]+fl[1]+fl[0]))*sgn(alphaOrdR)-0.6123724356957944*(fr[7]+fr[6]+fr[5])+0.6123724356957944*(fl[7]+fl[6]+fl[5])+0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[1] = 0.3535533905932737*(alpha[0]*fUp[1]+fUp[0]*alpha[1]); 
  incr[2] = 0.3535533905932737*(alpha[1]*fUp[3]+alpha[0]*fUp[2]); 
  incr[3] = -0.6123724356957944*(alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[4] = 0.3535533905932737*(alpha[0]*fUp[3]+alpha[1]*fUp[2]); 
  incr[5] = -0.6123724356957944*(alpha[0]*fUp[1]+fUp[0]*alpha[1]); 
  incr[6] = -0.6123724356957944*(alpha[1]*fUp[3]+alpha[0]*fUp[2]); 
  incr[7] = -0.6123724356957944*(alpha[0]*fUp[3]+alpha[1]*fUp[2]); 

  outr[0] += incr[0]*dv11r; 
  outr[1] += incr[1]*dv11r; 
  outr[2] += incr[2]*dv11r; 
  outr[3] += incr[3]*dv11r; 
  outr[4] += incr[4]*dv11r; 
  outr[5] += incr[5]*dv11r; 
  outr[6] += incr[6]*dv11r; 
  outr[7] += incr[7]*dv11r; 

  outl[0] += -1.0*incr[0]*dv11l; 
  outl[1] += -1.0*incr[1]*dv11l; 
  outl[2] += -1.0*incr[2]*dv11l; 
  outl[3] += incr[3]*dv11l; 
  outl[4] += -1.0*incr[4]*dv11l; 
  outl[5] += incr[5]*dv11l; 
  outl[6] += incr[6]*dv11l; 
  outl[7] += incr[7]*dv11l; 


  return std::abs(amid); 
} 
