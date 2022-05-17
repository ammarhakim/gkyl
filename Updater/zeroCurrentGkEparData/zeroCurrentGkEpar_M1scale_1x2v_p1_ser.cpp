#include <zeroCurrentGkEpar_mod_decl.h>

void zeroCurrentGkEpar_M1scale_1x2v_p1_ser(const double *xvc, const double *dxv, const double vparMax, const double *alpha, double *distf) 
{ 
  // xvc: cell center.
  // dxv: cell length.
  // vparmax: maximum vpar of the grid.
  // alpha: vpar scaling factor.
  // distf: distribution function (electron).
  if ((fabs(0.7071067811865475*alpha[0]-1.224744871391589*alpha[1]) < 1./vparMax) && (fabs(0.7071067811865475*alpha[0]) < 1./vparMax) && (fabs(1.224744871391589*alpha[1]+0.7071067811865475*alpha[0]) < 1./vparMax)) {

  distf[0] = 0.2041241452319315*alpha[1]*dxv[1]*distf[4]+0.2041241452319315*alpha[0]*dxv[1]*distf[2]+0.7071067811865475*alpha[1]*distf[1]*xvc[1]+0.7071067811865475*alpha[0]*distf[0]*xvc[1]+distf[0]; 
  distf[1] = 0.2041241452319315*alpha[0]*dxv[1]*distf[4]+0.2041241452319315*alpha[1]*dxv[1]*distf[2]+0.7071067811865475*alpha[0]*distf[1]*xvc[1]+0.7071067811865475*distf[0]*alpha[1]*xvc[1]+distf[1]; 
  distf[2] = 0.7071067811865475*alpha[1]*xvc[1]*distf[4]+0.7071067811865475*alpha[0]*xvc[1]*distf[2]+distf[2]+0.2041241452319315*alpha[1]*distf[1]*dxv[1]+0.2041241452319315*alpha[0]*distf[0]*dxv[1]; 
  distf[3] = 0.2041241452319315*alpha[1]*dxv[1]*distf[7]+0.2041241452319315*alpha[0]*dxv[1]*distf[6]+0.7071067811865475*alpha[1]*xvc[1]*distf[5]+0.7071067811865475*alpha[0]*xvc[1]*distf[3]+distf[3]; 
  distf[4] = 0.7071067811865475*alpha[0]*xvc[1]*distf[4]+distf[4]+0.7071067811865475*alpha[1]*xvc[1]*distf[2]+0.2041241452319315*alpha[0]*distf[1]*dxv[1]+0.2041241452319315*distf[0]*alpha[1]*dxv[1]; 
  distf[5] = 0.2041241452319315*alpha[0]*dxv[1]*distf[7]+0.2041241452319315*alpha[1]*dxv[1]*distf[6]+0.7071067811865475*alpha[0]*xvc[1]*distf[5]+distf[5]+0.7071067811865475*alpha[1]*xvc[1]*distf[3]; 
  distf[6] = 0.7071067811865475*alpha[1]*xvc[1]*distf[7]+0.7071067811865475*alpha[0]*xvc[1]*distf[6]+distf[6]+0.2041241452319315*alpha[1]*dxv[1]*distf[5]+0.2041241452319315*alpha[0]*dxv[1]*distf[3]; 
  distf[7] = 0.7071067811865475*alpha[0]*xvc[1]*distf[7]+distf[7]+0.7071067811865475*alpha[1]*xvc[1]*distf[6]+0.2041241452319315*alpha[0]*dxv[1]*distf[5]+0.2041241452319315*alpha[1]*dxv[1]*distf[3]; 

  }

}

