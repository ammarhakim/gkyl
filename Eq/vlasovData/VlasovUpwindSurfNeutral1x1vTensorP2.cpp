#include <VlasovModDecl.h> 
__host__ __device__ double VlasovUpwindSurfNeutral1x1vTensor_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
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
  const double *Fo0 = &EM[0]; 

  double alpha[3]; 
  double incr[9]; 

  alpha[0] = Fo0[0]; 
  alpha[1] = Fo0[1]; 
  alpha[2] = Fo0[2]; 

  const double amid = 0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2]; 

  double alphaOrdR;
  double fUpOrd[3];
  alphaOrdR = 0.6324555320336759*alpha[2]-0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0]; 
  fUpOrd[0] = 0.5*(fr[8]+fl[8]-1.5*(fr[7]+fl[7])-0.7745966692414833*fr[6]+0.7745966692414833*fl[6]+1.118033988749895*(fr[5]+fl[5])+0.4472135954999579*(fr[4]+fl[4])+1.161895003862225*fr[3]-1.161895003862225*fl[3]-0.8660254037844386*fr[2]+0.8660254037844386*fl[2]-0.6708203932499369*(fr[1]+fl[1])+0.5*(fr[0]+fl[0]))-0.5*(fr[8]-(fl[8]-(1.5*fl[7]-1.5*fr[7])+0.7745966692414833*(fr[6]+fl[6])-(1.118033988749895*fr[5]-1.118033988749895*fl[5])-0.4472135954999579*fr[4]+0.4472135954999579*fl[4]-(1.161895003862225*(fr[3]+fl[3])-0.8660254037844386*(fr[2]+fl[2]))+0.6708203932499369*fr[1]-0.6708203932499369*fl[1]-(0.5*fr[0]-0.5*fl[0])))*sgn(alphaOrdR); 
  alphaOrdR = 0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2]; 
  fUpOrd[1] = 0.5*((1.25*fr[8]-1.25*fl[8]-0.9682458365518543*(fr[6]+fl[6])-1.118033988749895*fr[5]+1.118033988749895*fl[5]+0.5590169943749475*fr[4]-0.5590169943749475*fl[4]+0.8660254037844386*(fr[2]+fl[2])-0.5*fr[0]+0.5*fl[0])*sgn(alphaOrdR)-1.25*(fr[8]+fl[8])+0.9682458365518543*fr[6]-0.9682458365518543*fl[6]+1.118033988749895*(fr[5]+fl[5])-0.5590169943749475*(fr[4]+fl[4])-0.8660254037844386*fr[2]+0.8660254037844386*fl[2]+0.5*(fr[0]+fl[0])); 
  alphaOrdR = 0.6324555320336759*alpha[2]+0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0]; 
  fUpOrd[2] = 0.5*(fr[8]+fl[8]+1.5*(fr[7]+fl[7])-0.7745966692414833*fr[6]+0.7745966692414833*fl[6]+1.118033988749895*(fr[5]+fl[5])+0.4472135954999579*(fr[4]+fl[4])-1.161895003862225*fr[3]+1.161895003862225*fl[3]-0.8660254037844386*fr[2]+0.8660254037844386*fl[2]+0.6708203932499369*(fr[1]+fl[1])+0.5*(fr[0]+fl[0]))-0.5*(fr[8]-(fl[8]-(1.5*fr[7]-1.5*fl[7])+0.7745966692414833*(fr[6]+fl[6])-(1.118033988749895*fr[5]-1.118033988749895*fl[5])-0.4472135954999579*fr[4]+0.4472135954999579*fl[4]-((-1.161895003862225*(fr[3]+fl[3]))-0.8660254037844386*(fr[2]+fl[2]))-0.6708203932499369*fr[1]+0.6708203932499369*fl[1]-(0.5*fr[0]-0.5*fl[0])))*sgn(alphaOrdR); 

  double fUp[3];
  fUp[0] = 0.07856742013183861*(5.0*fUpOrd[2]+8.0*fUpOrd[1]+5.0*fUpOrd[0]); 
  fUp[1] = 0.5270462766947298*(fUpOrd[2]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3513641844631533*(fUpOrd[2]-2.0*fUpOrd[1]+fUpOrd[0]); 

  incr[0] = 0.5*(alpha[2]*fUp[2]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[1] = 0.1*(4.47213595499958*(alpha[1]*fUp[2]+fUp[1]*alpha[2])+5.0*(alpha[0]*fUp[1]+fUp[0]*alpha[1])); 
  incr[2] = -0.8660254037844386*(alpha[2]*fUp[2]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[3] = -0.1*(7.745966692414834*(alpha[1]*fUp[2]+fUp[1]*alpha[2])+8.660254037844386*(alpha[0]*fUp[1]+fUp[0]*alpha[1])); 
  incr[4] = 0.01428571428571429*(2.23606797749979*(10.0*alpha[2]*fUp[2]+14.0*alpha[1]*fUp[1])+35.0*(alpha[0]*fUp[2]+fUp[0]*alpha[2])); 
  incr[5] = 1.118033988749895*(alpha[2]*fUp[2]+alpha[1]*fUp[1]+alpha[0]*fUp[0]); 
  incr[6] = -0.05532833351724881*(15.65247584249853*(alpha[0]*fUp[2]+fUp[0]*alpha[2])+10.0*alpha[2]*fUp[2]+14.0*alpha[1]*fUp[1]); 
  incr[7] = 0.1290994448735805*(7.745966692414834*(alpha[1]*fUp[2]+fUp[1]*alpha[2])+8.660254037844386*(alpha[0]*fUp[1]+fUp[0]*alpha[1])); 
  incr[8] = 0.07142857142857142*(15.65247584249853*(alpha[0]*fUp[2]+fUp[0]*alpha[2])+10.0*alpha[2]*fUp[2]+14.0*alpha[1]*fUp[1]); 

  outr[0] += incr[0]*dv10r; 
  outr[1] += incr[1]*dv10r; 
  outr[2] += incr[2]*dv10r; 
  outr[3] += incr[3]*dv10r; 
  outr[4] += incr[4]*dv10r; 
  outr[5] += incr[5]*dv10r; 
  outr[6] += incr[6]*dv10r; 
  outr[7] += incr[7]*dv10r; 
  outr[8] += incr[8]*dv10r; 

  outl[0] += -1.0*incr[0]*dv10l; 
  outl[1] += -1.0*incr[1]*dv10l; 
  outl[2] += incr[2]*dv10l; 
  outl[3] += incr[3]*dv10l; 
  outl[4] += -1.0*incr[4]*dv10l; 
  outl[5] += -1.0*incr[5]*dv10l; 
  outl[6] += incr[6]*dv10l; 
  outl[7] += -1.0*incr[7]*dv10l; 
  outl[8] += -1.0*incr[8]*dv10l; 


  return std::abs(amid); 
} 
