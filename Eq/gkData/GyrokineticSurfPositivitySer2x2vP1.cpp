#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity2x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.125*geoZ[0]*(1.732050807568877*Phi[2]-3.0*Phi[3])*dfac_y; 

  double alpha[8]; 
  alpha[0] = -0.5*geoZ[0]*(2.449489742783178*Phi[2]-4.242640687119286*Phi[3])*dfac_y; 
  double f0Quad[8]; 
  double f1Quad[8]; 
  double limQuad[8]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[0] = 0.1767766952966368*((fr[14]-1.0*(fl[14]+fr[10])+fl[10]-1.0*fr[9]+fl[9]-1.0*fr[7]+fl[7]+fr[4]-1.0*fl[4]+fr[3]-1.0*fl[3]+fr[2]-1.0*(fl[2]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14])+fr[10]+fl[10]+fr[9]+fl[9]+fr[7]+fl[7]-1.0*(fr[4]+fl[4]+fr[3]+fl[3]+fr[2]+fl[2])+fr[0]+fl[0]); 
  f1Quad[0] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[13]+fl[13]+fr[12]+fl[12]+fr[11]+fl[11])+fr[8]+fl[8]+fr[6]+fl[6]+fr[5]+fl[5]-1.0*(fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[13]-1.0*fl[13]+fr[12]-1.0*fl[12]+fr[11]-1.0*(fl[11]+fr[8])+fl[8]-1.0*fr[6]+fl[6]-1.0*fr[5]+fl[5]+fr[1]-1.0*fl[1]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[1] = -0.1767766952966368*((fr[14]-1.0*fl[14]+fr[10]-1.0*(fl[10]+fr[9])+fl[9]-1.0*fr[7]+fl[7]-1.0*fr[4]+fl[4]-1.0*fr[3]+fl[3]+fr[2]-1.0*fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14]+fr[10]+fl[10])+fr[9]+fl[9]+fr[7]+fl[7]+fr[4]+fl[4]+fr[3]+fl[3]-1.0*(fr[2]+fl[2]+fr[0]+fl[0])); 
  f1Quad[1] = -0.1767766952966368*((fr[15]+fl[15]+fr[13]+fl[13]-1.0*(fr[12]+fl[12]+fr[11]+fl[11]+fr[8]+fl[8]+fr[6]+fl[6])+fr[5]+fl[5]+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[13]+fl[13]+fr[12]-1.0*fl[12]+fr[11]-1.0*fl[11]+fr[8]-1.0*fl[8]+fr[6]-1.0*(fl[6]+fr[5])+fl[5]-1.0*fr[1]+fl[1]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[2] = -0.1767766952966368*((fr[14]-1.0*(fl[14]+fr[10])+fl[10]+fr[9]-1.0*(fl[9]+fr[7])+fl[7]-1.0*fr[4]+fl[4]+fr[3]-1.0*(fl[3]+fr[2])+fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14])+fr[10]+fl[10]-1.0*(fr[9]+fl[9])+fr[7]+fl[7]+fr[4]+fl[4]-1.0*(fr[3]+fl[3])+fr[2]+fl[2]-1.0*(fr[0]+fl[0])); 
  f1Quad[2] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[13]+fl[13])+fr[12]+fl[12]-1.0*(fr[11]+fl[11]+fr[8]+fl[8])+fr[6]+fl[6]-1.0*(fr[5]+fl[5])+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[13]-1.0*(fl[13]+fr[12])+fl[12]+fr[11]-1.0*fl[11]+fr[8]-1.0*(fl[8]+fr[6])+fl[6]+fr[5]-1.0*(fl[5]+fr[1])+fl[1]); 
  limQuad[2] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[3] = 0.1767766952966368*((fr[14]-1.0*fl[14]+fr[10]-1.0*fl[10]+fr[9]-1.0*(fl[9]+fr[7])+fl[7]+fr[4]-1.0*(fl[4]+fr[3])+fl[3]-1.0*fr[2]+fl[2]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14]+fr[10]+fl[10]+fr[9]+fl[9])+fr[7]+fl[7]-1.0*(fr[4]+fl[4])+fr[3]+fl[3]+fr[2]+fl[2]+fr[0]+fl[0]); 
  f1Quad[3] = 0.1767766952966368*((fr[15]+fl[15]+fr[13]+fl[13]+fr[12]+fl[12]-1.0*(fr[11]+fl[11])+fr[8]+fl[8]-1.0*(fr[6]+fl[6]+fr[5]+fl[5]+fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[13]+fl[13]-1.0*fr[12]+fl[12]+fr[11]-1.0*(fl[11]+fr[8])+fl[8]+fr[6]-1.0*fl[6]+fr[5]-1.0*fl[5]+fr[1]-1.0*fl[1]); 
  limQuad[3] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[4] = -0.1767766952966368*((fr[14]-1.0*(fl[14]+fr[10])+fl[10]-1.0*fr[9]+fl[9]+fr[7]-1.0*fl[7]+fr[4]-1.0*(fl[4]+fr[3])+fl[3]-1.0*fr[2]+fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14])+fr[10]+fl[10]+fr[9]+fl[9]-1.0*(fr[7]+fl[7]+fr[4]+fl[4])+fr[3]+fl[3]+fr[2]+fl[2]-1.0*(fr[0]+fl[0])); 
  f1Quad[4] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[13]+fl[13]+fr[12]+fl[12])+fr[11]+fl[11]+fr[8]+fl[8]-1.0*(fr[6]+fl[6]+fr[5]+fl[5])+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[13]-1.0*fl[13]+fr[12]-1.0*(fl[12]+fr[11])+fl[11]-1.0*fr[8]+fl[8]+fr[6]-1.0*fl[6]+fr[5]-1.0*(fl[5]+fr[1])+fl[1]); 
  limQuad[4] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[5] = 0.1767766952966368*((fr[14]-1.0*fl[14]+fr[10]-1.0*(fl[10]+fr[9])+fl[9]+fr[7]-1.0*(fl[7]+fr[4])+fl[4]+fr[3]-1.0*(fl[3]+fr[2])+fl[2]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14]+fr[10]+fl[10])+fr[9]+fl[9]-1.0*(fr[7]+fl[7])+fr[4]+fl[4]-1.0*(fr[3]+fl[3])+fr[2]+fl[2]+fr[0]+fl[0]); 
  f1Quad[5] = 0.1767766952966368*((fr[15]+fl[15]+fr[13]+fl[13]-1.0*(fr[12]+fl[12])+fr[11]+fl[11]-1.0*(fr[8]+fl[8])+fr[6]+fl[6]-1.0*(fr[5]+fl[5]+fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[13]+fl[13]+fr[12]-1.0*(fl[12]+fr[11])+fl[11]+fr[8]-1.0*(fl[8]+fr[6])+fl[6]+fr[5]-1.0*fl[5]+fr[1]-1.0*fl[1]); 
  limQuad[5] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[6] = 0.1767766952966368*((fr[14]-1.0*(fl[14]+fr[10])+fl[10]+fr[9]-1.0*fl[9]+fr[7]-1.0*(fl[7]+fr[4])+fl[4]-1.0*fr[3]+fl[3]+fr[2]-1.0*(fl[2]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14])+fr[10]+fl[10]-1.0*(fr[9]+fl[9]+fr[7]+fl[7])+fr[4]+fl[4]+fr[3]+fl[3]-1.0*(fr[2]+fl[2])+fr[0]+fl[0]); 
  f1Quad[6] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[13]+fl[13])+fr[12]+fl[12]+fr[11]+fl[11]-1.0*(fr[8]+fl[8]+fr[6]+fl[6])+fr[5]+fl[5]-1.0*(fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[13]-1.0*(fl[13]+fr[12])+fl[12]-1.0*fr[11]+fl[11]+fr[8]-1.0*fl[8]+fr[6]-1.0*(fl[6]+fr[5])+fl[5]+fr[1]-1.0*fl[1]); 
  limQuad[6] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[7] = -0.1767766952966368*((fr[14]-1.0*fl[14]+fr[10]-1.0*fl[10]+fr[9]-1.0*fl[9]+fr[7]-1.0*fl[7]+fr[4]-1.0*fl[4]+fr[3]-1.0*fl[3]+fr[2]-1.0*fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14]+fr[10]+fl[10]+fr[9]+fl[9]+fr[7]+fl[7]+fr[4]+fl[4]+fr[3]+fl[3]+fr[2]+fl[2]+fr[0]+fl[0])); 
  f1Quad[7] = -0.1767766952966368*((fr[15]+fl[15]+fr[13]+fl[13]+fr[12]+fl[12]+fr[11]+fl[11]+fr[8]+fl[8]+fr[6]+fl[6]+fr[5]+fl[5]+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[13]+fl[13]-1.0*fr[12]+fl[12]-1.0*fr[11]+fl[11]-1.0*fr[8]+fl[8]-1.0*fr[6]+fl[6]-1.0*fr[5]+fl[5]-1.0*fr[1]+fl[1]); 
  limQuad[7] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  double fhat[16]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[2] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[3] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[6] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[7] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[8] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[9] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[10] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[11] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[12] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[13] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[14] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[15] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[8];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrl[0] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[13]+fhat[12]+fhat[11]+3.0*fhat[1])+9.0*(fhat[8]+fhat[6]+fhat[5])))/(36.0*EPSILON-1.732050807568877*(fhat[14]+3.0*(fhat[4]+fhat[3]+fhat[2]))+3.0*(fhat[10]+fhat[9]+fhat[7]+3.0*fhat[0])); 
  rCtrl[1] = (3.0*fhat[15]+5.196152422706631*(fhat[13]-1.0*(fhat[12]+fhat[11])+3.0*fhat[1])+9.0*(fhat[5]-1.0*(fhat[8]+fhat[6])))/(36.0*EPSILON+1.732050807568877*(fhat[14]+3.0*(fhat[2]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0]))); 
  rCtrl[2] = (3.0*fhat[15]-5.196152422706631*(fhat[13]-1.0*fhat[12]+fhat[11]-3.0*fhat[1])+9.0*((-1.0*fhat[8])+fhat[6]-1.0*fhat[5]))/(36.0*EPSILON+1.732050807568877*(fhat[14]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[2]))+3.0*((-1.0*fhat[10])+fhat[9]-1.0*fhat[7]+3.0*fhat[0])); 
  rCtrl[3] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[13]+fhat[12]-1.0*(fhat[11]+3.0*fhat[1]))+9.0*(fhat[8]-1.0*(fhat[6]+fhat[5]))))/(36.0*EPSILON-1.732050807568877*(fhat[14]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[2])))+3.0*((-1.0*(fhat[10]+fhat[9]))+fhat[7]+3.0*fhat[0])); 
  rCtrl[4] = (3.0*fhat[15]-5.196152422706631*(fhat[13]+fhat[12]-1.0*(fhat[11]+3.0*fhat[1]))+9.0*(fhat[8]-1.0*(fhat[6]+fhat[5])))/(36.0*EPSILON+1.732050807568877*(fhat[14]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[2])))+3.0*((-1.0*(fhat[10]+fhat[9]))+fhat[7]+3.0*fhat[0])); 
  rCtrl[5] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[13]-1.0*fhat[12]+fhat[11]-3.0*fhat[1])+9.0*((-1.0*fhat[8])+fhat[6]-1.0*fhat[5])))/(36.0*EPSILON-1.732050807568877*(fhat[14]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[2]))+3.0*((-1.0*fhat[10])+fhat[9]-1.0*fhat[7]+3.0*fhat[0])); 
  rCtrl[6] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[13]-1.0*(fhat[12]+fhat[11])+3.0*fhat[1])+9.0*(fhat[5]-1.0*(fhat[8]+fhat[6]))))/(36.0*EPSILON-1.732050807568877*(fhat[14]+3.0*(fhat[2]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0]))); 
  rCtrl[7] = (3.0*fhat[15]+5.196152422706631*(fhat[13]+fhat[12]+fhat[11]+3.0*fhat[1])+9.0*(fhat[8]+fhat[6]+fhat[5]))/(36.0*EPSILON+1.732050807568877*(fhat[14]+3.0*(fhat[4]+fhat[3]+fhat[2]))+3.0*(fhat[10]+fhat[9]+fhat[7]+3.0*fhat[0])); 
  double fhatCtrl[8];  // fhatCtrl = anti-limited fhat evaluated at each control node on x surface 
  fhatCtrl[0] = -0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[4]+fhat[3]+fhat[2]))-3.0*(fhat[10]+fhat[9]+fhat[7]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[2]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[9]+fhat[7])+3.0*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[2]))-3.0*(fhat[10]-1.0*fhat[9]+fhat[7]-3.0*fhat[0]))*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = -0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[2])))+3.0*(fhat[10]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0])))*limTheta(rCtrl[3],-1.0); 
  fhatCtrl[4] = 0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[2])))-3.0*(fhat[10]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0])))*limTheta(rCtrl[4],-1.0); 
  fhatCtrl[5] = -0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[2]))+3.0*(fhat[10]-1.0*fhat[9]+fhat[7]-3.0*fhat[0]))*limTheta(rCtrl[5],-1.0); 
  fhatCtrl[6] = -0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[2]-1.0*(fhat[4]+fhat[3])))-3.0*(fhat[10]-1.0*(fhat[9]+fhat[7])+3.0*fhat[0]))*limTheta(rCtrl[6],-1.0); 
  fhatCtrl[7] = 0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[4]+fhat[3]+fhat[2]))+3.0*(fhat[10]+fhat[9]+fhat[7]+3.0*fhat[0]))*limTheta(rCtrl[7],-1.0); 
  double fhatAL[8];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.3535533905932737*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.6123724356957944*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 1.060660171779821*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[5] = 1.060660171779821*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[6] = 1.060660171779821*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 1.837117307087383*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[8]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[6]+1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])-1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[7]+fhatAL[1])-0.5773502691896258*(fhatAL[3]+fhatAL[2]))+1.414213562373095*(fhatAL[6]-1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*fhatAL[2]-0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[1]))+1.414213562373095*((-1.0*fhatAL[6])-0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[5])+fhatAL[0])), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[1]-1.0*fhatAL[7])-0.5773502691896258*fhatAL[3]+0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])+0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[5])+fhatAL[0])), limQuad[3])); 
  fhatALQuad[4] = std::max(0.0, std::min(0.25*(2.449489742783178*((-0.5773502691896258*(fhatAL[1]-1.0*fhatAL[7]))+0.5773502691896258*fhatAL[3]-0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])-0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[4])+fhatAL[0])), limQuad[4])); 
  fhatALQuad[5] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[1])-0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])+0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[4])+fhatAL[0])), limQuad[5])); 
  fhatALQuad[6] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[3]+fhatAL[2])-0.5773502691896258*(fhatAL[7]+fhatAL[1]))+1.414213562373095*(fhatAL[6]-1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[6])); 
  fhatALQuad[7] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])+1.414213562373095*(fhatAL[6]+1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[7])); 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.25*alpha[0]*fhatAL[0]*dfac_x; 
  incr[1] = -0.4330127018922193*alpha[0]*fhatAL[0]*dfac_x; 
  incr[2] = 0.25*alpha[0]*fhatAL[1]*dfac_x; 
  incr[3] = 0.25*alpha[0]*fhatAL[2]*dfac_x; 
  incr[4] = 0.25*alpha[0]*fhatAL[3]*dfac_x; 
  incr[5] = -0.4330127018922193*alpha[0]*fhatAL[1]*dfac_x; 
  incr[6] = -0.4330127018922193*alpha[0]*fhatAL[2]*dfac_x; 
  incr[7] = 0.25*alpha[0]*fhatAL[4]*dfac_x; 
  incr[8] = -0.4330127018922193*alpha[0]*fhatAL[3]*dfac_x; 
  incr[9] = 0.25*alpha[0]*fhatAL[5]*dfac_x; 
  incr[10] = 0.25*alpha[0]*fhatAL[6]*dfac_x; 
  incr[11] = -0.4330127018922193*alpha[0]*fhatAL[4]*dfac_x; 
  incr[12] = -0.4330127018922193*alpha[0]*fhatAL[5]*dfac_x; 
  incr[13] = -0.4330127018922193*alpha[0]*fhatAL[6]*dfac_x; 
  incr[14] = 0.25*alpha[0]*fhatAL[7]*dfac_x; 
  incr[15] = -0.4330127018922193*alpha[0]*fhatAL[7]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Y_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.125*geoZ[0]*(1.732050807568877*Phi[1]-3.0*Phi[3])*dfac_x; 

  double alpha[8]; 
  alpha[0] = 0.5*geoZ[0]*(2.449489742783178*Phi[1]-4.242640687119286*Phi[3])*dfac_x; 
  double f0Quad[8]; 
  double f1Quad[8]; 
  double limQuad[8]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[0] = 0.1767766952966368*((fr[13]-1.0*(fl[13]+fr[10])+fl[10]-1.0*fr[8]+fl[8]-1.0*fr[6]+fl[6]+fr[4]-1.0*fl[4]+fr[3]-1.0*fl[3]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13])+fr[10]+fl[10]+fr[8]+fl[8]+fr[6]+fl[6]-1.0*(fr[4]+fl[4]+fr[3]+fl[3]+fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[0] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14]+fr[12]+fl[12]+fr[11]+fl[11])+fr[9]+fl[9]+fr[7]+fl[7]+fr[5]+fl[5]-1.0*(fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*fl[14]+fr[12]-1.0*fl[12]+fr[11]-1.0*(fl[11]+fr[9])+fl[9]-1.0*fr[7]+fl[7]-1.0*fr[5]+fl[5]+fr[2]-1.0*fl[2]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[1] = -0.1767766952966368*((fr[13]-1.0*fl[13]+fr[10]-1.0*(fl[10]+fr[8])+fl[8]-1.0*fr[6]+fl[6]-1.0*fr[4]+fl[4]-1.0*fr[3]+fl[3]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13]+fr[10]+fl[10])+fr[8]+fl[8]+fr[6]+fl[6]+fr[4]+fl[4]+fr[3]+fl[3]-1.0*(fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[1] = -0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]-1.0*(fr[12]+fl[12]+fr[11]+fl[11]+fr[9]+fl[9]+fr[7]+fl[7])+fr[5]+fl[5]+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]+fr[12]-1.0*fl[12]+fr[11]-1.0*fl[11]+fr[9]-1.0*fl[9]+fr[7]-1.0*(fl[7]+fr[5])+fl[5]-1.0*fr[2]+fl[2]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[2] = -0.1767766952966368*((fr[13]-1.0*(fl[13]+fr[10])+fl[10]+fr[8]-1.0*(fl[8]+fr[6])+fl[6]-1.0*fr[4]+fl[4]+fr[3]-1.0*(fl[3]+fr[1])+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13])+fr[10]+fl[10]-1.0*(fr[8]+fl[8])+fr[6]+fl[6]+fr[4]+fl[4]-1.0*(fr[3]+fl[3])+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[2] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14])+fr[12]+fl[12]-1.0*(fr[11]+fl[11]+fr[9]+fl[9])+fr[7]+fl[7]-1.0*(fr[5]+fl[5])+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*(fl[14]+fr[12])+fl[12]+fr[11]-1.0*fl[11]+fr[9]-1.0*(fl[9]+fr[7])+fl[7]+fr[5]-1.0*(fl[5]+fr[2])+fl[2]); 
  limQuad[2] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[3] = 0.1767766952966368*((fr[13]-1.0*fl[13]+fr[10]-1.0*fl[10]+fr[8]-1.0*(fl[8]+fr[6])+fl[6]+fr[4]-1.0*(fl[4]+fr[3])+fl[3]-1.0*fr[1]+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13]+fr[10]+fl[10]+fr[8]+fl[8])+fr[6]+fl[6]-1.0*(fr[4]+fl[4])+fr[3]+fl[3]+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[3] = 0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]+fr[12]+fl[12]-1.0*(fr[11]+fl[11])+fr[9]+fl[9]-1.0*(fr[7]+fl[7]+fr[5]+fl[5]+fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]-1.0*fr[12]+fl[12]+fr[11]-1.0*(fl[11]+fr[9])+fl[9]+fr[7]-1.0*fl[7]+fr[5]-1.0*fl[5]+fr[2]-1.0*fl[2]); 
  limQuad[3] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[4] = -0.1767766952966368*((fr[13]-1.0*(fl[13]+fr[10])+fl[10]-1.0*fr[8]+fl[8]+fr[6]-1.0*fl[6]+fr[4]-1.0*(fl[4]+fr[3])+fl[3]-1.0*fr[1]+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13])+fr[10]+fl[10]+fr[8]+fl[8]-1.0*(fr[6]+fl[6]+fr[4]+fl[4])+fr[3]+fl[3]+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[4] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14]+fr[12]+fl[12])+fr[11]+fl[11]+fr[9]+fl[9]-1.0*(fr[7]+fl[7]+fr[5]+fl[5])+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*fl[14]+fr[12]-1.0*(fl[12]+fr[11])+fl[11]-1.0*fr[9]+fl[9]+fr[7]-1.0*fl[7]+fr[5]-1.0*(fl[5]+fr[2])+fl[2]); 
  limQuad[4] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[5] = 0.1767766952966368*((fr[13]-1.0*fl[13]+fr[10]-1.0*(fl[10]+fr[8])+fl[8]+fr[6]-1.0*(fl[6]+fr[4])+fl[4]+fr[3]-1.0*(fl[3]+fr[1])+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13]+fr[10]+fl[10])+fr[8]+fl[8]-1.0*(fr[6]+fl[6])+fr[4]+fl[4]-1.0*(fr[3]+fl[3])+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[5] = 0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]-1.0*(fr[12]+fl[12])+fr[11]+fl[11]-1.0*(fr[9]+fl[9])+fr[7]+fl[7]-1.0*(fr[5]+fl[5]+fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]+fr[12]-1.0*(fl[12]+fr[11])+fl[11]+fr[9]-1.0*(fl[9]+fr[7])+fl[7]+fr[5]-1.0*fl[5]+fr[2]-1.0*fl[2]); 
  limQuad[5] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[6] = 0.1767766952966368*((fr[13]-1.0*(fl[13]+fr[10])+fl[10]+fr[8]-1.0*fl[8]+fr[6]-1.0*(fl[6]+fr[4])+fl[4]-1.0*fr[3]+fl[3]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13])+fr[10]+fl[10]-1.0*(fr[8]+fl[8]+fr[6]+fl[6])+fr[4]+fl[4]+fr[3]+fl[3]-1.0*(fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[6] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14])+fr[12]+fl[12]+fr[11]+fl[11]-1.0*(fr[9]+fl[9]+fr[7]+fl[7])+fr[5]+fl[5]-1.0*(fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*(fl[14]+fr[12])+fl[12]-1.0*fr[11]+fl[11]+fr[9]-1.0*fl[9]+fr[7]-1.0*(fl[7]+fr[5])+fl[5]+fr[2]-1.0*fl[2]); 
  limQuad[6] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[7] = -0.1767766952966368*((fr[13]-1.0*fl[13]+fr[10]-1.0*fl[10]+fr[8]-1.0*fl[8]+fr[6]-1.0*fl[6]+fr[4]-1.0*fl[4]+fr[3]-1.0*fl[3]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13]+fr[10]+fl[10]+fr[8]+fl[8]+fr[6]+fl[6]+fr[4]+fl[4]+fr[3]+fl[3]+fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[7] = -0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]+fr[12]+fl[12]+fr[11]+fl[11]+fr[9]+fl[9]+fr[7]+fl[7]+fr[5]+fl[5]+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]-1.0*fr[12]+fl[12]-1.0*fr[11]+fl[11]-1.0*fr[9]+fl[9]-1.0*fr[7]+fl[7]-1.0*fr[5]+fl[5]-1.0*fr[2]+fl[2]); 
  limQuad[7] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  double fhat[16]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[3] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[6] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[7] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[8] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[9] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[10] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[11] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[12] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[13] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[14] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[15] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[8];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrl[0] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[14]+fhat[12]+fhat[11]+3.0*fhat[2])+9.0*(fhat[9]+fhat[7]+fhat[5])))/(36.0*EPSILON-1.732050807568877*(fhat[13]+3.0*(fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[10]+fhat[8]+fhat[6]+3.0*fhat[0])); 
  rCtrl[1] = (3.0*fhat[15]+5.196152422706631*(fhat[14]-1.0*(fhat[12]+fhat[11])+3.0*fhat[2])+9.0*(fhat[5]-1.0*(fhat[9]+fhat[7])))/(36.0*EPSILON+1.732050807568877*(fhat[13]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[8]+fhat[6]-3.0*fhat[0]))); 
  rCtrl[2] = (3.0*fhat[15]-5.196152422706631*(fhat[14]-1.0*fhat[12]+fhat[11]-3.0*fhat[2])+9.0*((-1.0*fhat[9])+fhat[7]-1.0*fhat[5]))/(36.0*EPSILON+1.732050807568877*(fhat[13]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[1]))+3.0*((-1.0*fhat[10])+fhat[8]-1.0*fhat[6]+3.0*fhat[0])); 
  rCtrl[3] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[14]+fhat[12]-1.0*(fhat[11]+3.0*fhat[2]))+9.0*(fhat[9]-1.0*(fhat[7]+fhat[5]))))/(36.0*EPSILON-1.732050807568877*(fhat[13]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[1])))+3.0*((-1.0*(fhat[10]+fhat[8]))+fhat[6]+3.0*fhat[0])); 
  rCtrl[4] = (3.0*fhat[15]-5.196152422706631*(fhat[14]+fhat[12]-1.0*(fhat[11]+3.0*fhat[2]))+9.0*(fhat[9]-1.0*(fhat[7]+fhat[5])))/(36.0*EPSILON+1.732050807568877*(fhat[13]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[1])))+3.0*((-1.0*(fhat[10]+fhat[8]))+fhat[6]+3.0*fhat[0])); 
  rCtrl[5] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[14]-1.0*fhat[12]+fhat[11]-3.0*fhat[2])+9.0*((-1.0*fhat[9])+fhat[7]-1.0*fhat[5])))/(36.0*EPSILON-1.732050807568877*(fhat[13]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[1]))+3.0*((-1.0*fhat[10])+fhat[8]-1.0*fhat[6]+3.0*fhat[0])); 
  rCtrl[6] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[14]-1.0*(fhat[12]+fhat[11])+3.0*fhat[2])+9.0*(fhat[5]-1.0*(fhat[9]+fhat[7]))))/(36.0*EPSILON-1.732050807568877*(fhat[13]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[8]+fhat[6]-3.0*fhat[0]))); 
  rCtrl[7] = (3.0*fhat[15]+5.196152422706631*(fhat[14]+fhat[12]+fhat[11]+3.0*fhat[2])+9.0*(fhat[9]+fhat[7]+fhat[5]))/(36.0*EPSILON+1.732050807568877*(fhat[13]+3.0*(fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[10]+fhat[8]+fhat[6]+3.0*fhat[0])); 
  double fhatCtrl[8];  // fhatCtrl = anti-limited fhat evaluated at each control node on y surface 
  fhatCtrl[0] = -0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[4]+fhat[3]+fhat[1]))-3.0*(fhat[10]+fhat[8]+fhat[6]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[8]+fhat[6])+3.0*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[1]))-3.0*(fhat[10]-1.0*fhat[8]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = -0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[1])))+3.0*(fhat[10]+fhat[8]-1.0*(fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[3],-1.0); 
  fhatCtrl[4] = 0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[1])))-3.0*(fhat[10]+fhat[8]-1.0*(fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[4],-1.0); 
  fhatCtrl[5] = -0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[1]))+3.0*(fhat[10]-1.0*fhat[8]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[5],-1.0); 
  fhatCtrl[6] = -0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[3])))-3.0*(fhat[10]-1.0*(fhat[8]+fhat[6])+3.0*fhat[0]))*limTheta(rCtrl[6],-1.0); 
  fhatCtrl[7] = 0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[10]+fhat[8]+fhat[6]+3.0*fhat[0]))*limTheta(rCtrl[7],-1.0); 
  double fhatAL[8];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.3535533905932737*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.6123724356957944*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 1.060660171779821*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[5] = 1.060660171779821*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[6] = 1.060660171779821*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 1.837117307087383*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[8]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[6]+1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])-1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[7]+fhatAL[1])-0.5773502691896258*(fhatAL[3]+fhatAL[2]))+1.414213562373095*(fhatAL[6]-1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*fhatAL[2]-0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[1]))+1.414213562373095*((-1.0*fhatAL[6])-0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[5])+fhatAL[0])), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[1]-1.0*fhatAL[7])-0.5773502691896258*fhatAL[3]+0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])+0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[5])+fhatAL[0])), limQuad[3])); 
  fhatALQuad[4] = std::max(0.0, std::min(0.25*(2.449489742783178*((-0.5773502691896258*(fhatAL[1]-1.0*fhatAL[7]))+0.5773502691896258*fhatAL[3]-0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])-0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[4])+fhatAL[0])), limQuad[4])); 
  fhatALQuad[5] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[1])-0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])+0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[4])+fhatAL[0])), limQuad[5])); 
  fhatALQuad[6] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[3]+fhatAL[2])-0.5773502691896258*(fhatAL[7]+fhatAL[1]))+1.414213562373095*(fhatAL[6]-1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[6])); 
  fhatALQuad[7] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])+1.414213562373095*(fhatAL[6]+1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[7])); 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.25*alpha[0]*fhatAL[0]*dfac_y; 
  incr[1] = 0.25*alpha[0]*fhatAL[1]*dfac_y; 
  incr[2] = -0.4330127018922193*alpha[0]*fhatAL[0]*dfac_y; 
  incr[3] = 0.25*alpha[0]*fhatAL[2]*dfac_y; 
  incr[4] = 0.25*alpha[0]*fhatAL[3]*dfac_y; 
  incr[5] = -0.4330127018922193*alpha[0]*fhatAL[1]*dfac_y; 
  incr[6] = 0.25*alpha[0]*fhatAL[4]*dfac_y; 
  incr[7] = -0.4330127018922193*alpha[0]*fhatAL[2]*dfac_y; 
  incr[8] = 0.25*alpha[0]*fhatAL[5]*dfac_y; 
  incr[9] = -0.4330127018922193*alpha[0]*fhatAL[3]*dfac_y; 
  incr[10] = 0.25*alpha[0]*fhatAL[6]*dfac_y; 
  incr[11] = -0.4330127018922193*alpha[0]*fhatAL[4]*dfac_y; 
  incr[12] = -0.4330127018922193*alpha[0]*fhatAL[5]*dfac_y; 
  incr[13] = 0.25*alpha[0]*fhatAL[7]*dfac_y; 
  incr[14] = -0.4330127018922193*alpha[0]*fhatAL[6]*dfac_y; 
  incr[15] = -0.4330127018922193*alpha[0]*fhatAL[7]*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.0; 

  double alpha[8]; 
  // alpha == 0, so nothing to do 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.125*(1.732050807568877*(3.0*geoZ[1]*Phi[3]+geoZ[0]*Phi[2])-3.0*(geoZ[0]*Phi[3]+geoZ[1]*Phi[2]))*dfac_y; 

  double alpha[8]; 
  alpha[0] = -0.5*(2.449489742783178*(3.0*geoZ[1]*Phi[3]+geoZ[0]*Phi[2])-4.242640687119286*(geoZ[0]*Phi[3]+geoZ[1]*Phi[2]))*dfac_y; 
  double f0Quad[8]; 
  double f1Quad[8]; 
  double limQuad[8]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[0] = 0.1767766952966368*((fr[14]-1.0*(fl[14]+fr[10])+fl[10]-1.0*fr[9]+fl[9]-1.0*fr[7]+fl[7]+fr[4]-1.0*fl[4]+fr[3]-1.0*fl[3]+fr[2]-1.0*(fl[2]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14])+fr[10]+fl[10]+fr[9]+fl[9]+fr[7]+fl[7]-1.0*(fr[4]+fl[4]+fr[3]+fl[3]+fr[2]+fl[2])+fr[0]+fl[0]); 
  f1Quad[0] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[13]+fl[13]+fr[12]+fl[12]+fr[11]+fl[11])+fr[8]+fl[8]+fr[6]+fl[6]+fr[5]+fl[5]-1.0*(fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[13]-1.0*fl[13]+fr[12]-1.0*fl[12]+fr[11]-1.0*(fl[11]+fr[8])+fl[8]-1.0*fr[6]+fl[6]-1.0*fr[5]+fl[5]+fr[1]-1.0*fl[1]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[1] = -0.1767766952966368*((fr[14]-1.0*fl[14]+fr[10]-1.0*(fl[10]+fr[9])+fl[9]-1.0*fr[7]+fl[7]-1.0*fr[4]+fl[4]-1.0*fr[3]+fl[3]+fr[2]-1.0*fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14]+fr[10]+fl[10])+fr[9]+fl[9]+fr[7]+fl[7]+fr[4]+fl[4]+fr[3]+fl[3]-1.0*(fr[2]+fl[2]+fr[0]+fl[0])); 
  f1Quad[1] = -0.1767766952966368*((fr[15]+fl[15]+fr[13]+fl[13]-1.0*(fr[12]+fl[12]+fr[11]+fl[11]+fr[8]+fl[8]+fr[6]+fl[6])+fr[5]+fl[5]+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[13]+fl[13]+fr[12]-1.0*fl[12]+fr[11]-1.0*fl[11]+fr[8]-1.0*fl[8]+fr[6]-1.0*(fl[6]+fr[5])+fl[5]-1.0*fr[1]+fl[1]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[2] = -0.1767766952966368*((fr[14]-1.0*(fl[14]+fr[10])+fl[10]+fr[9]-1.0*(fl[9]+fr[7])+fl[7]-1.0*fr[4]+fl[4]+fr[3]-1.0*(fl[3]+fr[2])+fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14])+fr[10]+fl[10]-1.0*(fr[9]+fl[9])+fr[7]+fl[7]+fr[4]+fl[4]-1.0*(fr[3]+fl[3])+fr[2]+fl[2]-1.0*(fr[0]+fl[0])); 
  f1Quad[2] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[13]+fl[13])+fr[12]+fl[12]-1.0*(fr[11]+fl[11]+fr[8]+fl[8])+fr[6]+fl[6]-1.0*(fr[5]+fl[5])+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[13]-1.0*(fl[13]+fr[12])+fl[12]+fr[11]-1.0*fl[11]+fr[8]-1.0*(fl[8]+fr[6])+fl[6]+fr[5]-1.0*(fl[5]+fr[1])+fl[1]); 
  limQuad[2] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[3] = 0.1767766952966368*((fr[14]-1.0*fl[14]+fr[10]-1.0*fl[10]+fr[9]-1.0*(fl[9]+fr[7])+fl[7]+fr[4]-1.0*(fl[4]+fr[3])+fl[3]-1.0*fr[2]+fl[2]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14]+fr[10]+fl[10]+fr[9]+fl[9])+fr[7]+fl[7]-1.0*(fr[4]+fl[4])+fr[3]+fl[3]+fr[2]+fl[2]+fr[0]+fl[0]); 
  f1Quad[3] = 0.1767766952966368*((fr[15]+fl[15]+fr[13]+fl[13]+fr[12]+fl[12]-1.0*(fr[11]+fl[11])+fr[8]+fl[8]-1.0*(fr[6]+fl[6]+fr[5]+fl[5]+fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[13]+fl[13]-1.0*fr[12]+fl[12]+fr[11]-1.0*(fl[11]+fr[8])+fl[8]+fr[6]-1.0*fl[6]+fr[5]-1.0*fl[5]+fr[1]-1.0*fl[1]); 
  limQuad[3] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[4] = -0.1767766952966368*((fr[14]-1.0*(fl[14]+fr[10])+fl[10]-1.0*fr[9]+fl[9]+fr[7]-1.0*fl[7]+fr[4]-1.0*(fl[4]+fr[3])+fl[3]-1.0*fr[2]+fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14])+fr[10]+fl[10]+fr[9]+fl[9]-1.0*(fr[7]+fl[7]+fr[4]+fl[4])+fr[3]+fl[3]+fr[2]+fl[2]-1.0*(fr[0]+fl[0])); 
  f1Quad[4] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[13]+fl[13]+fr[12]+fl[12])+fr[11]+fl[11]+fr[8]+fl[8]-1.0*(fr[6]+fl[6]+fr[5]+fl[5])+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[13]-1.0*fl[13]+fr[12]-1.0*(fl[12]+fr[11])+fl[11]-1.0*fr[8]+fl[8]+fr[6]-1.0*fl[6]+fr[5]-1.0*(fl[5]+fr[1])+fl[1]); 
  limQuad[4] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[5] = 0.1767766952966368*((fr[14]-1.0*fl[14]+fr[10]-1.0*(fl[10]+fr[9])+fl[9]+fr[7]-1.0*(fl[7]+fr[4])+fl[4]+fr[3]-1.0*(fl[3]+fr[2])+fl[2]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14]+fr[10]+fl[10])+fr[9]+fl[9]-1.0*(fr[7]+fl[7])+fr[4]+fl[4]-1.0*(fr[3]+fl[3])+fr[2]+fl[2]+fr[0]+fl[0]); 
  f1Quad[5] = 0.1767766952966368*((fr[15]+fl[15]+fr[13]+fl[13]-1.0*(fr[12]+fl[12])+fr[11]+fl[11]-1.0*(fr[8]+fl[8])+fr[6]+fl[6]-1.0*(fr[5]+fl[5]+fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[13]+fl[13]+fr[12]-1.0*(fl[12]+fr[11])+fl[11]+fr[8]-1.0*(fl[8]+fr[6])+fl[6]+fr[5]-1.0*fl[5]+fr[1]-1.0*fl[1]); 
  limQuad[5] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[6] = 0.1767766952966368*((fr[14]-1.0*(fl[14]+fr[10])+fl[10]+fr[9]-1.0*fl[9]+fr[7]-1.0*(fl[7]+fr[4])+fl[4]-1.0*fr[3]+fl[3]+fr[2]-1.0*(fl[2]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14])+fr[10]+fl[10]-1.0*(fr[9]+fl[9]+fr[7]+fl[7])+fr[4]+fl[4]+fr[3]+fl[3]-1.0*(fr[2]+fl[2])+fr[0]+fl[0]); 
  f1Quad[6] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[13]+fl[13])+fr[12]+fl[12]+fr[11]+fl[11]-1.0*(fr[8]+fl[8]+fr[6]+fl[6])+fr[5]+fl[5]-1.0*(fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[13]-1.0*(fl[13]+fr[12])+fl[12]-1.0*fr[11]+fl[11]+fr[8]-1.0*fl[8]+fr[6]-1.0*(fl[6]+fr[5])+fl[5]+fr[1]-1.0*fl[1]); 
  limQuad[6] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]; 
  f0Quad[7] = -0.1767766952966368*((fr[14]-1.0*fl[14]+fr[10]-1.0*fl[10]+fr[9]-1.0*fl[9]+fr[7]-1.0*fl[7]+fr[4]-1.0*fl[4]+fr[3]-1.0*fl[3]+fr[2]-1.0*fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[14]+fl[14]+fr[10]+fl[10]+fr[9]+fl[9]+fr[7]+fl[7]+fr[4]+fl[4]+fr[3]+fl[3]+fr[2]+fl[2]+fr[0]+fl[0])); 
  f1Quad[7] = -0.1767766952966368*((fr[15]+fl[15]+fr[13]+fl[13]+fr[12]+fl[12]+fr[11]+fl[11]+fr[8]+fl[8]+fr[6]+fl[6]+fr[5]+fl[5]+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[13]+fl[13]-1.0*fr[12]+fl[12]-1.0*fr[11]+fl[11]-1.0*fr[8]+fl[8]-1.0*fr[6]+fl[6]-1.0*fr[5]+fl[5]-1.0*fr[1]+fl[1]); 
  limQuad[7] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  double fhat[16]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[2] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[3] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[6] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[7] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[8] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[9] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[10] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[11] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[12] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[13] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[14] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[15] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[8];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrl[0] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[13]+fhat[12]+fhat[11]+3.0*fhat[1])+9.0*(fhat[8]+fhat[6]+fhat[5])))/(36.0*EPSILON-1.732050807568877*(fhat[14]+3.0*(fhat[4]+fhat[3]+fhat[2]))+3.0*(fhat[10]+fhat[9]+fhat[7]+3.0*fhat[0])); 
  rCtrl[1] = (3.0*fhat[15]+5.196152422706631*(fhat[13]-1.0*(fhat[12]+fhat[11])+3.0*fhat[1])+9.0*(fhat[5]-1.0*(fhat[8]+fhat[6])))/(36.0*EPSILON+1.732050807568877*(fhat[14]+3.0*(fhat[2]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0]))); 
  rCtrl[2] = (3.0*fhat[15]-5.196152422706631*(fhat[13]-1.0*fhat[12]+fhat[11]-3.0*fhat[1])+9.0*((-1.0*fhat[8])+fhat[6]-1.0*fhat[5]))/(36.0*EPSILON+1.732050807568877*(fhat[14]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[2]))+3.0*((-1.0*fhat[10])+fhat[9]-1.0*fhat[7]+3.0*fhat[0])); 
  rCtrl[3] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[13]+fhat[12]-1.0*(fhat[11]+3.0*fhat[1]))+9.0*(fhat[8]-1.0*(fhat[6]+fhat[5]))))/(36.0*EPSILON-1.732050807568877*(fhat[14]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[2])))+3.0*((-1.0*(fhat[10]+fhat[9]))+fhat[7]+3.0*fhat[0])); 
  rCtrl[4] = (3.0*fhat[15]-5.196152422706631*(fhat[13]+fhat[12]-1.0*(fhat[11]+3.0*fhat[1]))+9.0*(fhat[8]-1.0*(fhat[6]+fhat[5])))/(36.0*EPSILON+1.732050807568877*(fhat[14]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[2])))+3.0*((-1.0*(fhat[10]+fhat[9]))+fhat[7]+3.0*fhat[0])); 
  rCtrl[5] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[13]-1.0*fhat[12]+fhat[11]-3.0*fhat[1])+9.0*((-1.0*fhat[8])+fhat[6]-1.0*fhat[5])))/(36.0*EPSILON-1.732050807568877*(fhat[14]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[2]))+3.0*((-1.0*fhat[10])+fhat[9]-1.0*fhat[7]+3.0*fhat[0])); 
  rCtrl[6] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[13]-1.0*(fhat[12]+fhat[11])+3.0*fhat[1])+9.0*(fhat[5]-1.0*(fhat[8]+fhat[6]))))/(36.0*EPSILON-1.732050807568877*(fhat[14]+3.0*(fhat[2]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[9]+fhat[7]-3.0*fhat[0]))); 
  rCtrl[7] = (3.0*fhat[15]+5.196152422706631*(fhat[13]+fhat[12]+fhat[11]+3.0*fhat[1])+9.0*(fhat[8]+fhat[6]+fhat[5]))/(36.0*EPSILON+1.732050807568877*(fhat[14]+3.0*(fhat[4]+fhat[3]+fhat[2]))+3.0*(fhat[10]+fhat[9]+fhat[7]+3.0*fhat[0])); 
  double fhatCtrl[8];  // fhatCtrl = anti-limited fhat evaluated at each control node on x surface 
  fhatCtrl[0] = -0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[4]+fhat[3]+fhat[2]))-3.0*(fhat[10]+fhat[9]+fhat[7]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[2]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[9]+fhat[7])+3.0*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[2]))-3.0*(fhat[10]-1.0*fhat[9]+fhat[7]-3.0*fhat[0]))*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = -0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[2])))+3.0*(fhat[10]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0])))*limTheta(rCtrl[3],-1.0); 
  fhatCtrl[4] = 0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[2])))-3.0*(fhat[10]+fhat[9]-1.0*(fhat[7]+3.0*fhat[0])))*limTheta(rCtrl[4],-1.0); 
  fhatCtrl[5] = -0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[2]))+3.0*(fhat[10]-1.0*fhat[9]+fhat[7]-3.0*fhat[0]))*limTheta(rCtrl[5],-1.0); 
  fhatCtrl[6] = -0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[2]-1.0*(fhat[4]+fhat[3])))-3.0*(fhat[10]-1.0*(fhat[9]+fhat[7])+3.0*fhat[0]))*limTheta(rCtrl[6],-1.0); 
  fhatCtrl[7] = 0.02777777777777778*(1.732050807568877*(fhat[14]+3.0*(fhat[4]+fhat[3]+fhat[2]))+3.0*(fhat[10]+fhat[9]+fhat[7]+3.0*fhat[0]))*limTheta(rCtrl[7],-1.0); 
  double fhatAL[8];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.3535533905932737*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.6123724356957944*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 1.060660171779821*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[5] = 1.060660171779821*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[6] = 1.060660171779821*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 1.837117307087383*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[8]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[6]+1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])-1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[7]+fhatAL[1])-0.5773502691896258*(fhatAL[3]+fhatAL[2]))+1.414213562373095*(fhatAL[6]-1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*fhatAL[2]-0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[1]))+1.414213562373095*((-1.0*fhatAL[6])-0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[5])+fhatAL[0])), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[1]-1.0*fhatAL[7])-0.5773502691896258*fhatAL[3]+0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])+0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[5])+fhatAL[0])), limQuad[3])); 
  fhatALQuad[4] = std::max(0.0, std::min(0.25*(2.449489742783178*((-0.5773502691896258*(fhatAL[1]-1.0*fhatAL[7]))+0.5773502691896258*fhatAL[3]-0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])-0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[4])+fhatAL[0])), limQuad[4])); 
  fhatALQuad[5] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[1])-0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])+0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[4])+fhatAL[0])), limQuad[5])); 
  fhatALQuad[6] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[3]+fhatAL[2])-0.5773502691896258*(fhatAL[7]+fhatAL[1]))+1.414213562373095*(fhatAL[6]-1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[6])); 
  fhatALQuad[7] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])+1.414213562373095*(fhatAL[6]+1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[7])); 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.25*alpha[0]*fhatAL[0]*dfac_x; 
  incr[1] = -0.4330127018922193*alpha[0]*fhatAL[0]*dfac_x; 
  incr[2] = 0.25*alpha[0]*fhatAL[1]*dfac_x; 
  incr[3] = 0.25*alpha[0]*fhatAL[2]*dfac_x; 
  incr[4] = 0.25*alpha[0]*fhatAL[3]*dfac_x; 
  incr[5] = -0.4330127018922193*alpha[0]*fhatAL[1]*dfac_x; 
  incr[6] = -0.4330127018922193*alpha[0]*fhatAL[2]*dfac_x; 
  incr[7] = 0.25*alpha[0]*fhatAL[4]*dfac_x; 
  incr[8] = -0.4330127018922193*alpha[0]*fhatAL[3]*dfac_x; 
  incr[9] = 0.25*alpha[0]*fhatAL[5]*dfac_x; 
  incr[10] = 0.25*alpha[0]*fhatAL[6]*dfac_x; 
  incr[11] = -0.4330127018922193*alpha[0]*fhatAL[4]*dfac_x; 
  incr[12] = -0.4330127018922193*alpha[0]*fhatAL[5]*dfac_x; 
  incr[13] = -0.4330127018922193*alpha[0]*fhatAL[6]*dfac_x; 
  incr[14] = 0.25*alpha[0]*fhatAL[7]*dfac_x; 
  incr[15] = -0.4330127018922193*alpha[0]*fhatAL[7]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Y_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.0625*dfac_x*(1.732050807568877*(Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*m_*wv2+2.0*geoZ[0]*(Bmag[1]*wm+Phi[1]*q_))-6.0*geoZ[0]*Phi[3]*q_))/q_; 

  double alpha[8]; 
  alpha[0] = (0.3535533905932737*dfac_x*(1.732050807568877*(Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*m_*wv2+2.0*geoZ[0]*(Bmag[1]*wm+Phi[1]*q_))-6.0*geoZ[0]*Phi[3]*q_))/q_; 
  alpha[1] = (0.3535533905932737*dfac_x*(1.732050807568877*(Bmag[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*m_*wv2+2.0*geoZ[1]*(Bmag[1]*wm+Phi[1]*q_))-6.0*geoZ[1]*Phi[3]*q_))/q_; 
  alpha[2] = (0.3535533905932737*Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*dfac_x*m_*wv)/(dfac_v*q_); 
  alpha[3] = (0.7071067811865475*geoZ[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alpha[4] = (0.3535533905932737*Bmag[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*dfac_x*m_*wv)/(dfac_v*q_); 
  alpha[5] = (0.7071067811865475*Bmag[1]*geoZ[1]*dfac_x)/(dfac_m*q_); 
  double f0Quad[8]; 
  double f1Quad[8]; 
  double limQuad[8]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.3535533905932737*(alpha[5]+alpha[4])-0.3535533905932737*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  f0Quad[0] = 0.1767766952966368*((fr[13]-1.0*(fl[13]+fr[10])+fl[10]-1.0*fr[8]+fl[8]-1.0*fr[6]+fl[6]+fr[4]-1.0*fl[4]+fr[3]-1.0*fl[3]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13])+fr[10]+fl[10]+fr[8]+fl[8]+fr[6]+fl[6]-1.0*(fr[4]+fl[4]+fr[3]+fl[3]+fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[0] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14]+fr[12]+fl[12]+fr[11]+fl[11])+fr[9]+fl[9]+fr[7]+fl[7]+fr[5]+fl[5]-1.0*(fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*fl[14]+fr[12]-1.0*fl[12]+fr[11]-1.0*(fl[11]+fr[9])+fl[9]-1.0*fr[7]+fl[7]-1.0*fr[5]+fl[5]+fr[2]-1.0*fl[2]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0])-0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2]); 
  f0Quad[1] = -0.1767766952966368*((fr[13]-1.0*fl[13]+fr[10]-1.0*(fl[10]+fr[8])+fl[8]-1.0*fr[6]+fl[6]-1.0*fr[4]+fl[4]-1.0*fr[3]+fl[3]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13]+fr[10]+fl[10])+fr[8]+fl[8]+fr[6]+fl[6]+fr[4]+fl[4]+fr[3]+fl[3]-1.0*(fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[1] = -0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]-1.0*(fr[12]+fl[12]+fr[11]+fl[11]+fr[9]+fl[9]+fr[7]+fl[7])+fr[5]+fl[5]+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]+fr[12]-1.0*fl[12]+fr[11]-1.0*fl[11]+fr[9]-1.0*fl[9]+fr[7]-1.0*(fl[7]+fr[5])+fl[5]-1.0*fr[2]+fl[2]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[5]-0.3535533905932737*(alpha[4]+alpha[3])+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  f0Quad[2] = -0.1767766952966368*((fr[13]-1.0*(fl[13]+fr[10])+fl[10]+fr[8]-1.0*(fl[8]+fr[6])+fl[6]-1.0*fr[4]+fl[4]+fr[3]-1.0*(fl[3]+fr[1])+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13])+fr[10]+fl[10]-1.0*(fr[8]+fl[8])+fr[6]+fl[6]+fr[4]+fl[4]-1.0*(fr[3]+fl[3])+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[2] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14])+fr[12]+fl[12]-1.0*(fr[11]+fl[11]+fr[9]+fl[9])+fr[7]+fl[7]-1.0*(fr[5]+fl[5])+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*(fl[14]+fr[12])+fl[12]+fr[11]-1.0*fl[11]+fr[9]-1.0*(fl[9]+fr[7])+fl[7]+fr[5]-1.0*(fl[5]+fr[2])+fl[2]); 
  limQuad[2] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = (-0.3535533905932737*alpha[5])+0.3535533905932737*alpha[4]-0.3535533905932737*alpha[3]+0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]); 
  f0Quad[3] = 0.1767766952966368*((fr[13]-1.0*fl[13]+fr[10]-1.0*fl[10]+fr[8]-1.0*(fl[8]+fr[6])+fl[6]+fr[4]-1.0*(fl[4]+fr[3])+fl[3]-1.0*fr[1]+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13]+fr[10]+fl[10]+fr[8]+fl[8])+fr[6]+fl[6]-1.0*(fr[4]+fl[4])+fr[3]+fl[3]+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[3] = 0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]+fr[12]+fl[12]-1.0*(fr[11]+fl[11])+fr[9]+fl[9]-1.0*(fr[7]+fl[7]+fr[5]+fl[5]+fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]-1.0*fr[12]+fl[12]+fr[11]-1.0*(fl[11]+fr[9])+fl[9]+fr[7]-1.0*fl[7]+fr[5]-1.0*fl[5]+fr[2]-1.0*fl[2]); 
  limQuad[3] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = (-0.3535533905932737*alpha[5])+0.3535533905932737*(alpha[4]+alpha[3])-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0]; 
  f0Quad[4] = -0.1767766952966368*((fr[13]-1.0*(fl[13]+fr[10])+fl[10]-1.0*fr[8]+fl[8]+fr[6]-1.0*fl[6]+fr[4]-1.0*(fl[4]+fr[3])+fl[3]-1.0*fr[1]+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13])+fr[10]+fl[10]+fr[8]+fl[8]-1.0*(fr[6]+fl[6]+fr[4]+fl[4])+fr[3]+fl[3]+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[4] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14]+fr[12]+fl[12])+fr[11]+fl[11]+fr[9]+fl[9]-1.0*(fr[7]+fl[7]+fr[5]+fl[5])+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*fl[14]+fr[12]-1.0*(fl[12]+fr[11])+fl[11]-1.0*fr[9]+fl[9]+fr[7]-1.0*fl[7]+fr[5]-1.0*(fl[5]+fr[2])+fl[2]); 
  limQuad[4] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[5]-0.3535533905932737*alpha[4]+0.3535533905932737*alpha[3]-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]); 
  f0Quad[5] = 0.1767766952966368*((fr[13]-1.0*fl[13]+fr[10]-1.0*(fl[10]+fr[8])+fl[8]+fr[6]-1.0*(fl[6]+fr[4])+fl[4]+fr[3]-1.0*(fl[3]+fr[1])+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13]+fr[10]+fl[10])+fr[8]+fl[8]-1.0*(fr[6]+fl[6])+fr[4]+fl[4]-1.0*(fr[3]+fl[3])+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[5] = 0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]-1.0*(fr[12]+fl[12])+fr[11]+fl[11]-1.0*(fr[9]+fl[9])+fr[7]+fl[7]-1.0*(fr[5]+fl[5]+fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]+fr[12]-1.0*(fl[12]+fr[11])+fl[11]+fr[9]-1.0*(fl[9]+fr[7])+fl[7]+fr[5]-1.0*fl[5]+fr[2]-1.0*fl[2]); 
  limQuad[5] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = (-0.3535533905932737*(alpha[5]+alpha[4]))+0.3535533905932737*(alpha[3]+alpha[2])-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0]; 
  f0Quad[6] = 0.1767766952966368*((fr[13]-1.0*(fl[13]+fr[10])+fl[10]+fr[8]-1.0*fl[8]+fr[6]-1.0*(fl[6]+fr[4])+fl[4]-1.0*fr[3]+fl[3]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13])+fr[10]+fl[10]-1.0*(fr[8]+fl[8]+fr[6]+fl[6])+fr[4]+fl[4]+fr[3]+fl[3]-1.0*(fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[6] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14])+fr[12]+fl[12]+fr[11]+fl[11]-1.0*(fr[9]+fl[9]+fr[7]+fl[7])+fr[5]+fl[5]-1.0*(fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*(fl[14]+fr[12])+fl[12]-1.0*fr[11]+fl[11]+fr[9]-1.0*fl[9]+fr[7]-1.0*(fl[7]+fr[5])+fl[5]+fr[2]-1.0*fl[2]); 
  limQuad[6] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*(alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  f0Quad[7] = -0.1767766952966368*((fr[13]-1.0*fl[13]+fr[10]-1.0*fl[10]+fr[8]-1.0*fl[8]+fr[6]-1.0*fl[6]+fr[4]-1.0*fl[4]+fr[3]-1.0*fl[3]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[13]+fl[13]+fr[10]+fl[10]+fr[8]+fl[8]+fr[6]+fl[6]+fr[4]+fl[4]+fr[3]+fl[3]+fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[7] = -0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]+fr[12]+fl[12]+fr[11]+fl[11]+fr[9]+fl[9]+fr[7]+fl[7]+fr[5]+fl[5]+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]-1.0*fr[12]+fl[12]-1.0*fr[11]+fl[11]-1.0*fr[9]+fl[9]-1.0*fr[7]+fl[7]-1.0*fr[5]+fl[5]-1.0*fr[2]+fl[2]); 
  limQuad[7] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  double fhat[16]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[3] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[6] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[7] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[8] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[9] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[10] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[11] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[12] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[13] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[14] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[15] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[8];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrl[0] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[14]+fhat[12]+fhat[11]+3.0*fhat[2])+9.0*(fhat[9]+fhat[7]+fhat[5])))/(36.0*EPSILON-1.732050807568877*(fhat[13]+3.0*(fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[10]+fhat[8]+fhat[6]+3.0*fhat[0])); 
  rCtrl[1] = (3.0*fhat[15]+5.196152422706631*(fhat[14]-1.0*(fhat[12]+fhat[11])+3.0*fhat[2])+9.0*(fhat[5]-1.0*(fhat[9]+fhat[7])))/(36.0*EPSILON+1.732050807568877*(fhat[13]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[8]+fhat[6]-3.0*fhat[0]))); 
  rCtrl[2] = (3.0*fhat[15]-5.196152422706631*(fhat[14]-1.0*fhat[12]+fhat[11]-3.0*fhat[2])+9.0*((-1.0*fhat[9])+fhat[7]-1.0*fhat[5]))/(36.0*EPSILON+1.732050807568877*(fhat[13]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[1]))+3.0*((-1.0*fhat[10])+fhat[8]-1.0*fhat[6]+3.0*fhat[0])); 
  rCtrl[3] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[14]+fhat[12]-1.0*(fhat[11]+3.0*fhat[2]))+9.0*(fhat[9]-1.0*(fhat[7]+fhat[5]))))/(36.0*EPSILON-1.732050807568877*(fhat[13]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[1])))+3.0*((-1.0*(fhat[10]+fhat[8]))+fhat[6]+3.0*fhat[0])); 
  rCtrl[4] = (3.0*fhat[15]-5.196152422706631*(fhat[14]+fhat[12]-1.0*(fhat[11]+3.0*fhat[2]))+9.0*(fhat[9]-1.0*(fhat[7]+fhat[5])))/(36.0*EPSILON+1.732050807568877*(fhat[13]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[1])))+3.0*((-1.0*(fhat[10]+fhat[8]))+fhat[6]+3.0*fhat[0])); 
  rCtrl[5] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[14]-1.0*fhat[12]+fhat[11]-3.0*fhat[2])+9.0*((-1.0*fhat[9])+fhat[7]-1.0*fhat[5])))/(36.0*EPSILON-1.732050807568877*(fhat[13]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[1]))+3.0*((-1.0*fhat[10])+fhat[8]-1.0*fhat[6]+3.0*fhat[0])); 
  rCtrl[6] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[14]-1.0*(fhat[12]+fhat[11])+3.0*fhat[2])+9.0*(fhat[5]-1.0*(fhat[9]+fhat[7]))))/(36.0*EPSILON-1.732050807568877*(fhat[13]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[8]+fhat[6]-3.0*fhat[0]))); 
  rCtrl[7] = (3.0*fhat[15]+5.196152422706631*(fhat[14]+fhat[12]+fhat[11]+3.0*fhat[2])+9.0*(fhat[9]+fhat[7]+fhat[5]))/(36.0*EPSILON+1.732050807568877*(fhat[13]+3.0*(fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[10]+fhat[8]+fhat[6]+3.0*fhat[0])); 
  double fhatCtrl[8];  // fhatCtrl = anti-limited fhat evaluated at each control node on y surface 
  fhatCtrl[0] = -0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[4]+fhat[3]+fhat[1]))-3.0*(fhat[10]+fhat[8]+fhat[6]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[3])))+3.0*(fhat[10]-1.0*(fhat[8]+fhat[6])+3.0*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[1]))-3.0*(fhat[10]-1.0*fhat[8]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = -0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[1])))+3.0*(fhat[10]+fhat[8]-1.0*(fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[3],-1.0); 
  fhatCtrl[4] = 0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[4]-1.0*(fhat[3]+fhat[1])))-3.0*(fhat[10]+fhat[8]-1.0*(fhat[6]+3.0*fhat[0])))*limTheta(rCtrl[4],-1.0); 
  fhatCtrl[5] = -0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*((-1.0*fhat[4])+fhat[3]-1.0*fhat[1]))+3.0*(fhat[10]-1.0*fhat[8]+fhat[6]-3.0*fhat[0]))*limTheta(rCtrl[5],-1.0); 
  fhatCtrl[6] = -0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[3])))-3.0*(fhat[10]-1.0*(fhat[8]+fhat[6])+3.0*fhat[0]))*limTheta(rCtrl[6],-1.0); 
  fhatCtrl[7] = 0.02777777777777778*(1.732050807568877*(fhat[13]+3.0*(fhat[4]+fhat[3]+fhat[1]))+3.0*(fhat[10]+fhat[8]+fhat[6]+3.0*fhat[0]))*limTheta(rCtrl[7],-1.0); 
  double fhatAL[8];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.3535533905932737*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.6123724356957944*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 1.060660171779821*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[5] = 1.060660171779821*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[6] = 1.060660171779821*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 1.837117307087383*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[8]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[6]+1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])-1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[7]+fhatAL[1])-0.5773502691896258*(fhatAL[3]+fhatAL[2]))+1.414213562373095*(fhatAL[6]-1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*fhatAL[2]-0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[1]))+1.414213562373095*((-1.0*fhatAL[6])-0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[5])+fhatAL[0])), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[1]-1.0*fhatAL[7])-0.5773502691896258*fhatAL[3]+0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])+0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[5])+fhatAL[0])), limQuad[3])); 
  fhatALQuad[4] = std::max(0.0, std::min(0.25*(2.449489742783178*((-0.5773502691896258*(fhatAL[1]-1.0*fhatAL[7]))+0.5773502691896258*fhatAL[3]-0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])-0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[4])+fhatAL[0])), limQuad[4])); 
  fhatALQuad[5] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[1])-0.5773502691896258*fhatAL[2])+1.414213562373095*((-1.0*fhatAL[6])+0.5773502691896258*(1.732050807568877*fhatAL[5]-1.732050807568877*fhatAL[4])+fhatAL[0])), limQuad[5])); 
  fhatALQuad[6] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[3]+fhatAL[2])-0.5773502691896258*(fhatAL[7]+fhatAL[1]))+1.414213562373095*(fhatAL[6]-1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[6])); 
  fhatALQuad[7] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])+1.414213562373095*(fhatAL[6]+1.0*(fhatAL[5]+fhatAL[4])+fhatAL[0])), limQuad[7])); 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.25*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_y; 
  incr[1] = 0.25*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_y; 
  incr[2] = -0.4330127018922193*(alpha[5]*fhatAL[5]+alpha[4]*fhatAL[4]+alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_y; 
  incr[3] = 0.25*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_y; 
  incr[4] = 0.25*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_y; 
  incr[5] = -0.4330127018922193*(alpha[3]*fhatAL[5]+fhatAL[3]*alpha[5]+alpha[2]*fhatAL[4]+fhatAL[2]*alpha[4]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_y; 
  incr[6] = 0.25*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_y; 
  incr[7] = -0.4330127018922193*(alpha[5]*fhatAL[7]+alpha[3]*fhatAL[6]+alpha[1]*fhatAL[4]+fhatAL[1]*alpha[4]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_y; 
  incr[8] = 0.25*(alpha[2]*fhatAL[7]+alpha[4]*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_y; 
  incr[9] = -0.4330127018922193*(alpha[4]*fhatAL[7]+alpha[2]*fhatAL[6]+alpha[1]*fhatAL[5]+fhatAL[1]*alpha[5]+alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3])*dfac_y; 
  incr[10] = 0.25*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_y; 
  incr[11] = -0.4330127018922193*(alpha[3]*fhatAL[7]+alpha[5]*fhatAL[6]+alpha[0]*fhatAL[4]+fhatAL[0]*alpha[4]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_y; 
  incr[12] = -0.4330127018922193*(alpha[2]*fhatAL[7]+alpha[4]*fhatAL[6]+alpha[0]*fhatAL[5]+fhatAL[0]*alpha[5]+alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3])*dfac_y; 
  incr[13] = 0.25*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])*dfac_y; 
  incr[14] = -0.4330127018922193*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6]+alpha[4]*fhatAL[5]+fhatAL[4]*alpha[5]+alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3])*dfac_y; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6]+alpha[2]*fhatAL[5]+fhatAL[2]*alpha[5]+alpha[3]*fhatAL[4]+fhatAL[3]*alpha[4])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.09375*Bmag[1]*((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Phi[2])*dfac_x*dfac_y*(dfac_v*wv-1.0))/dfac_v; 

  double alpha[8]; 
  alpha[0] = -(0.1767766952966368*Bmag[1]*((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Phi[2])*dfac_x*dfac_y*(3.0*dfac_v*wv-3.0))/dfac_v; 
  alpha[1] = -(0.03535533905932736*Bmag[1]*dfac_x*dfac_y*(((27.0*BmagInv[1]*geoZ[1]+15.0*BmagInv[0]*geoZ[0])*Phi[3]+15.0*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[2])*dfac_v*wv+((-27.0*BmagInv[1]*geoZ[1])-15.0*BmagInv[0]*geoZ[0])*Phi[3]-15.0*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[2]))/dfac_v; 
  double f0Quad[8]; 
  double f1Quad[8]; 
  double limQuad[8]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1]; 
  f0Quad[0] = 0.1767766952966368*((fr[12]-1.0*(fl[12]+fr[9])+fl[9]-1.0*fr[8]+fl[8]-1.0*fr[5]+fl[5]+fr[4]-1.0*fl[4]+fr[2]-1.0*fl[2]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12])+fr[9]+fl[9]+fr[8]+fl[8]+fr[5]+fl[5]-1.0*(fr[4]+fl[4]+fr[2]+fl[2]+fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[0] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14]+fr[13]+fl[13]+fr[11]+fl[11])+fr[10]+fl[10]+fr[7]+fl[7]+fr[6]+fl[6]-1.0*(fr[3]+fl[3]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*fl[14]+fr[13]-1.0*fl[13]+fr[11]-1.0*(fl[11]+fr[10])+fl[10]-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]+fr[3]-1.0*fl[3]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0]); 
  f0Quad[1] = -0.1767766952966368*((fr[12]-1.0*fl[12]+fr[9]-1.0*(fl[9]+fr[8])+fl[8]-1.0*fr[5]+fl[5]-1.0*fr[4]+fl[4]-1.0*fr[2]+fl[2]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12]+fr[9]+fl[9])+fr[8]+fl[8]+fr[5]+fl[5]+fr[4]+fl[4]+fr[2]+fl[2]-1.0*(fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[1] = -0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]-1.0*(fr[13]+fl[13]+fr[11]+fl[11]+fr[10]+fl[10]+fr[7]+fl[7])+fr[6]+fl[6]+fr[3]+fl[3])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]+fr[13]-1.0*fl[13]+fr[11]-1.0*fl[11]+fr[10]-1.0*fl[10]+fr[7]-1.0*(fl[7]+fr[6])+fl[6]-1.0*fr[3]+fl[3]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1]; 
  f0Quad[2] = -0.1767766952966368*((fr[12]-1.0*(fl[12]+fr[9])+fl[9]+fr[8]-1.0*(fl[8]+fr[5])+fl[5]-1.0*fr[4]+fl[4]+fr[2]-1.0*(fl[2]+fr[1])+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12])+fr[9]+fl[9]-1.0*(fr[8]+fl[8])+fr[5]+fl[5]+fr[4]+fl[4]-1.0*(fr[2]+fl[2])+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[2] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14])+fr[13]+fl[13]-1.0*(fr[11]+fl[11]+fr[10]+fl[10])+fr[7]+fl[7]-1.0*(fr[6]+fl[6])+fr[3]+fl[3])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*(fl[14]+fr[13])+fl[13]+fr[11]-1.0*fl[11]+fr[10]-1.0*(fl[10]+fr[7])+fl[7]+fr[6]-1.0*(fl[6]+fr[3])+fl[3]); 
  limQuad[2] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0]); 
  f0Quad[3] = 0.1767766952966368*((fr[12]-1.0*fl[12]+fr[9]-1.0*fl[9]+fr[8]-1.0*(fl[8]+fr[5])+fl[5]+fr[4]-1.0*(fl[4]+fr[2])+fl[2]-1.0*fr[1]+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12]+fr[9]+fl[9]+fr[8]+fl[8])+fr[5]+fl[5]-1.0*(fr[4]+fl[4])+fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[3] = 0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]-1.0*(fr[11]+fl[11])+fr[10]+fl[10]-1.0*(fr[7]+fl[7]+fr[6]+fl[6]+fr[3]+fl[3]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]-1.0*fr[13]+fl[13]+fr[11]-1.0*(fl[11]+fr[10])+fl[10]+fr[7]-1.0*fl[7]+fr[6]-1.0*fl[6]+fr[3]-1.0*fl[3]); 
  limQuad[3] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1]; 
  f0Quad[4] = -0.1767766952966368*((fr[12]-1.0*(fl[12]+fr[9])+fl[9]-1.0*fr[8]+fl[8]+fr[5]-1.0*fl[5]+fr[4]-1.0*(fl[4]+fr[2])+fl[2]-1.0*fr[1]+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12])+fr[9]+fl[9]+fr[8]+fl[8]-1.0*(fr[5]+fl[5]+fr[4]+fl[4])+fr[2]+fl[2]+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[4] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14]+fr[13]+fl[13])+fr[11]+fl[11]+fr[10]+fl[10]-1.0*(fr[7]+fl[7]+fr[6]+fl[6])+fr[3]+fl[3])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*fl[14]+fr[13]-1.0*(fl[13]+fr[11])+fl[11]-1.0*fr[10]+fl[10]+fr[7]-1.0*fl[7]+fr[6]-1.0*(fl[6]+fr[3])+fl[3]); 
  limQuad[4] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0]); 
  f0Quad[5] = 0.1767766952966368*((fr[12]-1.0*fl[12]+fr[9]-1.0*(fl[9]+fr[8])+fl[8]+fr[5]-1.0*(fl[5]+fr[4])+fl[4]+fr[2]-1.0*(fl[2]+fr[1])+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12]+fr[9]+fl[9])+fr[8]+fl[8]-1.0*(fr[5]+fl[5])+fr[4]+fl[4]-1.0*(fr[2]+fl[2])+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[5] = 0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]-1.0*(fr[13]+fl[13])+fr[11]+fl[11]-1.0*(fr[10]+fl[10])+fr[7]+fl[7]-1.0*(fr[6]+fl[6]+fr[3]+fl[3]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]+fr[13]-1.0*(fl[13]+fr[11])+fl[11]+fr[10]-1.0*(fl[10]+fr[7])+fl[7]+fr[6]-1.0*fl[6]+fr[3]-1.0*fl[3]); 
  limQuad[5] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*alpha[0]-0.3535533905932737*alpha[1]; 
  f0Quad[6] = 0.1767766952966368*((fr[12]-1.0*(fl[12]+fr[9])+fl[9]+fr[8]-1.0*fl[8]+fr[5]-1.0*(fl[5]+fr[4])+fl[4]-1.0*fr[2]+fl[2]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12])+fr[9]+fl[9]-1.0*(fr[8]+fl[8]+fr[5]+fl[5])+fr[4]+fl[4]+fr[2]+fl[2]-1.0*(fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[6] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14])+fr[13]+fl[13]+fr[11]+fl[11]-1.0*(fr[10]+fl[10]+fr[7]+fl[7])+fr[6]+fl[6]-1.0*(fr[3]+fl[3]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*(fl[14]+fr[13])+fl[13]-1.0*fr[11]+fl[11]+fr[10]-1.0*fl[10]+fr[7]-1.0*(fl[7]+fr[6])+fl[6]+fr[3]-1.0*fl[3]); 
  limQuad[6] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  alphaQuad = 0.3535533905932737*(alpha[1]+alpha[0]); 
  f0Quad[7] = -0.1767766952966368*((fr[12]-1.0*fl[12]+fr[9]-1.0*fl[9]+fr[8]-1.0*fl[8]+fr[5]-1.0*fl[5]+fr[4]-1.0*fl[4]+fr[2]-1.0*fl[2]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12]+fr[9]+fl[9]+fr[8]+fl[8]+fr[5]+fl[5]+fr[4]+fl[4]+fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[7] = -0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]+fr[11]+fl[11]+fr[10]+fl[10]+fr[7]+fl[7]+fr[6]+fl[6]+fr[3]+fl[3])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]-1.0*fr[13]+fl[13]-1.0*fr[11]+fl[11]-1.0*fr[10]+fl[10]-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]-1.0*fr[3]+fl[3]); 
  limQuad[7] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.25; 
  double fhat[16]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[3] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[4] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[6] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[7] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[8] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[9] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[10] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[11] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[12] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[13] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[14] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[15] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[8];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrl[0] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[14]+fhat[13]+fhat[11]+3.0*fhat[3])+9.0*(fhat[10]+fhat[7]+fhat[6])))/(36.0*EPSILON-1.732050807568877*(fhat[12]+3.0*(fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[9]+fhat[8]+fhat[5]+3.0*fhat[0])); 
  rCtrl[1] = (3.0*fhat[15]+5.196152422706631*(fhat[14]-1.0*(fhat[13]+fhat[11])+3.0*fhat[3])+9.0*(fhat[6]-1.0*(fhat[10]+fhat[7])))/(36.0*EPSILON+1.732050807568877*(fhat[12]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[2])))+3.0*(fhat[9]-1.0*(fhat[8]+fhat[5]-3.0*fhat[0]))); 
  rCtrl[2] = (3.0*fhat[15]-5.196152422706631*(fhat[14]-1.0*fhat[13]+fhat[11]-3.0*fhat[3])+9.0*((-1.0*fhat[10])+fhat[7]-1.0*fhat[6]))/(36.0*EPSILON+1.732050807568877*(fhat[12]+3.0*((-1.0*fhat[4])+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*fhat[9])+fhat[8]-1.0*fhat[5]+3.0*fhat[0])); 
  rCtrl[3] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[14]+fhat[13]-1.0*(fhat[11]+3.0*fhat[3]))+9.0*(fhat[10]-1.0*(fhat[7]+fhat[6]))))/(36.0*EPSILON-1.732050807568877*(fhat[12]+3.0*(fhat[4]-1.0*(fhat[2]+fhat[1])))+3.0*((-1.0*(fhat[9]+fhat[8]))+fhat[5]+3.0*fhat[0])); 
  rCtrl[4] = (3.0*fhat[15]-5.196152422706631*(fhat[14]+fhat[13]-1.0*(fhat[11]+3.0*fhat[3]))+9.0*(fhat[10]-1.0*(fhat[7]+fhat[6])))/(36.0*EPSILON+1.732050807568877*(fhat[12]+3.0*(fhat[4]-1.0*(fhat[2]+fhat[1])))+3.0*((-1.0*(fhat[9]+fhat[8]))+fhat[5]+3.0*fhat[0])); 
  rCtrl[5] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[14]-1.0*fhat[13]+fhat[11]-3.0*fhat[3])+9.0*((-1.0*fhat[10])+fhat[7]-1.0*fhat[6])))/(36.0*EPSILON-1.732050807568877*(fhat[12]+3.0*((-1.0*fhat[4])+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*fhat[9])+fhat[8]-1.0*fhat[5]+3.0*fhat[0])); 
  rCtrl[6] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[14]-1.0*(fhat[13]+fhat[11])+3.0*fhat[3])+9.0*(fhat[6]-1.0*(fhat[10]+fhat[7]))))/(36.0*EPSILON-1.732050807568877*(fhat[12]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[2])))+3.0*(fhat[9]-1.0*(fhat[8]+fhat[5]-3.0*fhat[0]))); 
  rCtrl[7] = (3.0*fhat[15]+5.196152422706631*(fhat[14]+fhat[13]+fhat[11]+3.0*fhat[3])+9.0*(fhat[10]+fhat[7]+fhat[6]))/(36.0*EPSILON+1.732050807568877*(fhat[12]+3.0*(fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[9]+fhat[8]+fhat[5]+3.0*fhat[0])); 
  double fhatCtrl[8];  // fhatCtrl = anti-limited fhat evaluated at each control node on vx surface 
  fhatCtrl[0] = -0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[4]+fhat[2]+fhat[1]))-3.0*(fhat[9]+fhat[8]+fhat[5]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[2])))+3.0*(fhat[9]-1.0*(fhat[8]+fhat[5])+3.0*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*((-1.0*fhat[4])+fhat[2]-1.0*fhat[1]))-3.0*(fhat[9]-1.0*fhat[8]+fhat[5]-3.0*fhat[0]))*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = -0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[4]-1.0*(fhat[2]+fhat[1])))+3.0*(fhat[9]+fhat[8]-1.0*(fhat[5]+3.0*fhat[0])))*limTheta(rCtrl[3],-1.0); 
  fhatCtrl[4] = 0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[4]-1.0*(fhat[2]+fhat[1])))-3.0*(fhat[9]+fhat[8]-1.0*(fhat[5]+3.0*fhat[0])))*limTheta(rCtrl[4],-1.0); 
  fhatCtrl[5] = -0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*((-1.0*fhat[4])+fhat[2]-1.0*fhat[1]))+3.0*(fhat[9]-1.0*fhat[8]+fhat[5]-3.0*fhat[0]))*limTheta(rCtrl[5],-1.0); 
  fhatCtrl[6] = -0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[2])))-3.0*(fhat[9]-1.0*(fhat[8]+fhat[5])+3.0*fhat[0]))*limTheta(rCtrl[6],-1.0); 
  fhatCtrl[7] = 0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[9]+fhat[8]+fhat[5]+3.0*fhat[0]))*limTheta(rCtrl[7],-1.0); 
  double fhatAL[8];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.3535533905932737*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.6123724356957944*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 1.060660171779821*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[5] = 1.060660171779821*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[6] = 1.060660171779821*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 1.837117307087383*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[8]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.25*(1.414213562373095*(1.0*(fhatAL[6]+fhatAL[4])+fhatAL[5]+fhatAL[0])-1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*fhatAL[1]-0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[2]))+1.414213562373095*((-0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[6]))-1.0*fhatAL[5]+fhatAL[0])), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[7]+fhatAL[2])-0.5773502691896258*(fhatAL[3]+fhatAL[1]))+1.414213562373095*((-1.0*(fhatAL[6]+fhatAL[4]))+fhatAL[5]+fhatAL[0])), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[2]-1.0*fhatAL[7])-0.5773502691896258*fhatAL[3]+0.5773502691896258*fhatAL[1])+1.414213562373095*(0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[6])-1.0*fhatAL[5]+fhatAL[0])), limQuad[3])); 
  fhatALQuad[4] = std::max(0.0, std::min(0.25*(2.449489742783178*((-0.5773502691896258*(fhatAL[2]-1.0*fhatAL[7]))+0.5773502691896258*fhatAL[3]-0.5773502691896258*fhatAL[1])+1.414213562373095*((-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[4]))-1.0*fhatAL[5]+fhatAL[0])), limQuad[4])); 
  fhatALQuad[5] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[3]+fhatAL[1])-0.5773502691896258*(fhatAL[7]+fhatAL[2]))+1.414213562373095*((-1.0*(fhatAL[6]+fhatAL[4]))+fhatAL[5]+fhatAL[0])), limQuad[5])); 
  fhatALQuad[6] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[2])-0.5773502691896258*fhatAL[1])+1.414213562373095*(0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[4])-1.0*fhatAL[5]+fhatAL[0])), limQuad[6])); 
  fhatALQuad[7] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])+1.414213562373095*(1.0*(fhatAL[6]+fhatAL[4])+fhatAL[5]+fhatAL[0])), limQuad[7])); 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.25*(alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[1] = 0.25*(alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(alpha[1]*fhatAL[4]+alpha[0]*fhatAL[2])*dfac_v; 
  incr[3] = -0.4330127018922193*(alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[4] = 0.25*(alpha[1]*fhatAL[5]+alpha[0]*fhatAL[3])*dfac_v; 
  incr[5] = 0.25*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2])*dfac_v; 
  incr[6] = -0.4330127018922193*(alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[7] = -0.4330127018922193*(alpha[1]*fhatAL[4]+alpha[0]*fhatAL[2])*dfac_v; 
  incr[8] = 0.25*(alpha[0]*fhatAL[5]+alpha[1]*fhatAL[3])*dfac_v; 
  incr[9] = 0.25*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6])*dfac_v; 
  incr[10] = -0.4330127018922193*(alpha[1]*fhatAL[5]+alpha[0]*fhatAL[3])*dfac_v; 
  incr[11] = -0.4330127018922193*(alpha[0]*fhatAL[4]+alpha[1]*fhatAL[2])*dfac_v; 
  incr[12] = 0.25*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6])*dfac_v; 
  incr[13] = -0.4330127018922193*(alpha[0]*fhatAL[5]+alpha[1]*fhatAL[3])*dfac_v; 
  incr[14] = -0.4330127018922193*(alpha[1]*fhatAL[7]+alpha[0]*fhatAL[6])*dfac_v; 
  incr[15] = -0.4330127018922193*(alpha[0]*fhatAL[7]+alpha[1]*fhatAL[6])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  return std::abs(alpha0); 
} 
