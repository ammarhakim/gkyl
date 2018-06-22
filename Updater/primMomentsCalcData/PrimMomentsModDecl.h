#ifndef PRIM_MOMENTS_MOD_DECL_H 
#define PRIM_MOMENTS_MOD_DECL_H 
 
// Eigen include statements. 
#include <Eigen/Dense> 
 
extern "C" { 
void SelfPrimMoments1x1vSer_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x1vSer_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x1vSer_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x1vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x1vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x1vSer_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x1vSer_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x1vSer_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x1vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x1vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x1vSer_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x1vSer_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x1vSer_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x1vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x1vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x2vSer_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x2vSer_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vSer_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vSer_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vSer_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x2vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x2vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x2vSer_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x2vSer_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vSer_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vSer_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vSer_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x2vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x2vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x2vSer_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x2vSer_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vSer_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vSer_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vSer_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x2vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x2vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x3vSer_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x3vSer_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_F_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_vF_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x3vSer_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x3vSer_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_F_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_vF_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x3vSer_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x3vSer_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_F_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vSer_vF_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 


void SelfPrimMoments2x2vSer_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x2vSer_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vSer_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vSer_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vSer_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x2vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x2vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x2vSer_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x2vSer_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vSer_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vSer_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vSer_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x2vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x2vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x2vSer_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x2vSer_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vSer_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vSer_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vSer_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x2vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x2vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x3vSer_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x3vSer_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_F_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_vF_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x3vSer_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x3vSer_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_F_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_vF_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x3vSer_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x3vSer_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_F_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vSer_vF_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 


void SelfPrimMoments3x3vSer_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral3x3vSer_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_F_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_vF_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_3x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_3x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments3x3vSer_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral3x3vSer_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_F_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_vF_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_3x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_3x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments3x3vSer_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral3x3vSer_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_F_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vSer_vF_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_3x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_3x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 


void SelfPrimMoments1x1vMax_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x1vMax_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x1vMax_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x1vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x1vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x1vMax_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x1vMax_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x1vMax_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x1vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x1vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x1vMax_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x1vMax_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x1vMax_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x1vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x1vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x2vMax_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x2vMax_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vMax_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vMax_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vMax_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x2vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x2vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x2vMax_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x2vMax_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vMax_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vMax_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vMax_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x2vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x2vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x2vMax_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x2vMax_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vMax_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vMax_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x2vMax_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x2vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x2vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x3vMax_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x3vMax_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_F_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_vF_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x3vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x3vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x3vMax_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x3vMax_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_F_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_vF_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x3vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x3vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x3vMax_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral1x3vMax_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_F_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral1x3vMax_vF_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_1x3vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x3vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 


void SelfPrimMoments2x2vMax_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x2vMax_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vMax_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vMax_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vMax_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x2vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x2vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x2vMax_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x2vMax_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vMax_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vMax_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vMax_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x2vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x2vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x2vMax_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x2vMax_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vMax_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vMax_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x2vMax_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x2vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x2vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x3vMax_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x3vMax_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_F_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_vF_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x3vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x3vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x3vMax_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x3vMax_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_F_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_vF_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x3vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x3vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x3vMax_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral2x3vMax_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_F_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral2x3vMax_vF_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_2x3vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x3vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 


void SelfPrimMoments3x3vMax_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral3x3vMax_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_F_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_vF_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_3x3vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_3x3vMax_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments3x3vMax_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral3x3vMax_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_F_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_vF_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_3x3vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_3x3vMax_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments3x3vMax_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0B, const double *m1B, double *u, double *vtSq); 

void BoundaryIntegral3x3vMax_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_F_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void BoundaryIntegral3x3vMax_vF_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out); 

void CrossPrimMoments_VmeiLBO_3x3vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_3x3vMax_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 


} 
#endif 
