#include <VlasovModDecl.h> 
__host__ __device__ double VlasovRecoverySurfElcMag1x1vTensor_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 

  double Gleft[4]; 
  double Gright[4]; 
  double incr[4]; 

  Gleft[0] = 0.7071067811865475*(E0[1]*fl[1]+E0[0]*fl[0]); 
  Gleft[1] = 0.7071067811865475*(E0[0]*fl[1]+fl[0]*E0[1]); 
  Gleft[2] = 0.7071067811865475*(E0[1]*fl[3]+E0[0]*fl[2]); 
  Gleft[3] = 0.7071067811865475*(E0[0]*fl[3]+E0[1]*fl[2]); 

  Gright[0] = 0.7071067811865475*(E0[1]*fr[1]+E0[0]*fr[0]); 
  Gright[1] = 0.7071067811865475*(E0[0]*fr[1]+fr[0]*E0[1]); 
  Gright[2] = 0.7071067811865475*(E0[1]*fr[3]+E0[0]*fr[2]); 
  Gright[3] = 0.7071067811865475*(E0[0]*fr[3]+E0[1]*fr[2]); 

  incr[0] = (-0.2886751345948129*Gright[2])+0.2886751345948129*Gleft[2]+0.25*(Gright[0]+Gleft[0]); 
  incr[1] = (-0.2886751345948129*Gright[3])+0.2886751345948129*Gleft[3]+0.25*(Gright[1]+Gleft[1]); 
  incr[2] = 0.5*Gright[2]-0.5*Gleft[2]-0.4330127018922193*(Gright[0]+Gleft[0]); 
  incr[3] = 0.5*Gright[3]-0.5*Gleft[3]-0.4330127018922193*(Gright[1]+Gleft[1]); 

  outr[0] += incr[0]*dv10r; 
  outr[1] += incr[1]*dv10r; 
  outr[2] += incr[2]*dv10r; 
  outr[3] += incr[3]*dv10r; 
  outl[0] += -1.0*incr[0]*dv10l; 
  outl[1] += -1.0*incr[1]*dv10l; 
  outl[2] += incr[2]*dv10l; 
  outl[3] += incr[3]*dv10l; 
  const double amid = 0.7071067811865475*E0[0]; 

  return std::abs(amid); 
} 
