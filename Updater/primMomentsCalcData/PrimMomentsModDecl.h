#ifndef PRIM_MOMENTS_MOD_DECL_H 
#define PRIM_MOMENTS_MOD_DECL_H 
 
// Eigen include statements. 
#include <Eigen/Dense> 
#include <CartFieldBinOpModDecl.h> 
 
extern "C" { 
void VmSelfPrimMoments1x1vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star1x1vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star1x1vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral1x1vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x1vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x1vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x1vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x2vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star1x2vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star1x2vSer_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star1x2vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral1x2vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vSer_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vSer_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x3vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star1x3vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star1x3vSer_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star1x3vSer_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star1x3vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral1x3vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x3vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x3vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x1vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral1x1vSer_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x1vSer_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x1vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x1vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x2vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral1x2vSer_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vSer_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vSer_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x2vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star2x2vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star2x2vSer_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star2x2vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral2x2vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vSer_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vSer_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x3vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star2x3vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star2x3vSer_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star2x3vSer_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star2x3vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral2x3vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x3vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x3vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments2x2vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral2x2vSer_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vSer_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vSer_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene2x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons2x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments3x3vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star3x3vSer_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star3x3vSer_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star3x3vSer_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star3x3vSer(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral3x3vSer_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene3x3vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons3x3vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments3x2vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral3x2vSer_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vSer_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vSer_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene3x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons3x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 


void VmSelfPrimMoments1x1vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x1vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x1vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x1vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x1vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x2vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x2vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vSer_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vSer_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x3vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x3vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x3vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x3vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x1vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral1x1vSer_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x1vSer_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x1vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x1vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x2vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral1x2vSer_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vSer_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vSer_vF_VY_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x2vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral2x2vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vSer_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vSer_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x3vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral2x3vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x3vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x3vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments2x2vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral2x2vSer_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vSer_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vSer_vF_VY_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene2x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons2x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments3x3vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral3x3vSer_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene3x3vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons3x3vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments3x2vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral3x2vSer_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vSer_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vSer_vF_VY_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene3x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons3x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 


void VmSelfPrimMoments1x1vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x1vSer_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x1vSer_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x1vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x1vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x2vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x2vSer_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vSer_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vSer_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vSer_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x2vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x2vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x3vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x3vSer_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_F_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vSer_vF_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x3vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x3vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x1vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral1x1vSer_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x1vSer_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x1vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x1vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x2vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral1x2vSer_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vSer_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vSer_vF_VY_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x2vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x2vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x2vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral2x2vSer_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vSer_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vSer_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vSer_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x2vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x2vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x3vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral2x3vSer_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_F_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vSer_vF_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x3vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x3vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments2x2vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral2x2vSer_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vSer_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vSer_vF_VY_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene2x2vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons2x2vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments3x3vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral3x3vSer_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_F_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vSer_vF_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene3x3vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons3x3vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments3x2vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral3x2vSer_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vSer_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vSer_vF_VY_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene3x2vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons3x2vSer_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 


void VmSelfPrimMoments1x1vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star1x1vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star1x1vMax(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral1x1vMax_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x1vMax_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x1vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x1vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x2vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star1x2vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star1x2vMax_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star1x2vMax(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral1x2vMax_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vMax_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vMax_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vMax_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x2vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x2vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x3vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star1x3vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star1x3vMax_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star1x3vMax_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star1x3vMax(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral1x3vMax_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x3vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x3vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x1vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void GkM0Star1x1vMax_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void GkM1iM2Star1x1vMax(const double *w, const double *dxv, const double intFac, const double m_, const double *Bmag, const double *f, double *outM1i, double *outM2); 

void GkBoundaryIntegral1x1vMax_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x1vMax_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x1vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x1vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x2vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void GkM0Star1x2vMax_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void GkM1iM2Star1x2vMax(const double *w, const double *dxv, const double intFac, const double m_, const double *Bmag, const double *f, double *outM1i, double *outM2); 

void GkBoundaryIntegral1x2vMax_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vMax_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vMax_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x2vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x2vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x2vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star2x2vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star2x2vMax_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star2x2vMax(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral2x2vMax_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vMax_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vMax_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vMax_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x2vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x2vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x3vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star2x3vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star2x3vMax_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star2x3vMax_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star2x3vMax(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral2x3vMax_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x3vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x3vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments2x2vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void GkM0Star2x2vMax_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void GkM1iM2Star2x2vMax(const double *w, const double *dxv, const double intFac, const double m_, const double *Bmag, const double *f, double *outM1i, double *outM2); 

void GkBoundaryIntegral2x2vMax_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vMax_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vMax_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene2x2vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons2x2vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments3x3vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void VmM0Star3x3vMax_VX(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star3x3vMax_VY(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM0Star3x3vMax_VZ(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void VmM1iM2Star3x3vMax(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2); 

void VmBoundaryIntegral3x3vMax_F_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_vF_VX_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_F_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_vF_VY_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_F_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_vF_VZ_P1(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene3x3vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons3x3vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments3x2vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq); 

void GkM0Star3x2vMax_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out); 

void GkM1iM2Star3x2vMax(const double *w, const double *dxv, const double intFac, const double m_, const double *Bmag, const double *f, double *outM1i, double *outM2); 

void GkBoundaryIntegral3x2vMax_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vMax_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vMax_vF_VY_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene3x2vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons3x2vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 


void VmSelfPrimMoments1x1vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x1vMax_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x1vMax_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x1vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x1vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x2vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x2vMax_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vMax_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vMax_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vMax_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x2vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x2vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x3vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x3vMax_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x3vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x3vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x1vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral1x1vMax_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x1vMax_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x1vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x1vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x2vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral1x2vMax_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vMax_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vMax_vF_VY_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x2vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x2vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x2vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral2x2vMax_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vMax_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vMax_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vMax_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x2vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x2vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x3vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral2x3vMax_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x3vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x3vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments2x2vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral2x2vMax_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vMax_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vMax_vF_VY_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene2x2vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons2x2vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments3x3vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral3x3vMax_F_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_vF_VX_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_F_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_vF_VY_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_F_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_vF_VZ_P2(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene3x3vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons3x3vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments3x2vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral3x2vMax_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vMax_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vMax_vF_VY_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene3x2vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons3x2vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 


void VmSelfPrimMoments1x1vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x1vMax_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x1vMax_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x1vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x1vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x2vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x2vMax_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vMax_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vMax_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x2vMax_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x2vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x2vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments1x3vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral1x3vMax_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_F_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral1x3vMax_vF_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene1x3vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons1x3vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x1vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral1x1vMax_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x1vMax_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x1vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x1vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments1x2vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral1x2vMax_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vMax_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral1x2vMax_vF_VY_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene1x2vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons1x2vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x2vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral2x2vMax_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vMax_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vMax_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x2vMax_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x2vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x2vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments2x3vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral2x3vMax_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_F_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral2x3vMax_vF_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene2x3vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons2x3vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments2x2vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral2x2vMax_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vMax_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral2x2vMax_vF_VY_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene2x2vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons2x2vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmSelfPrimMoments3x3vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void VmBoundaryIntegral3x3vMax_F_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_vF_VX_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_F_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_vF_VY_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_F_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmBoundaryIntegral3x3vMax_vF_VZ_P3(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void VmCrossPrimMomentsGreene3x3vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void VmCrossPrimMomentsHeavyIons3x3vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkSelfPrimMoments3x2vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq); 

void GkBoundaryIntegral3x2vMax_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vMax_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkBoundaryIntegral3x2vMax_vF_VY_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out); 

void GkCrossPrimMomentsGreene3x2vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 

void GkCrossPrimMomentsHeavyIons3x2vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross); 


} 
#endif 
