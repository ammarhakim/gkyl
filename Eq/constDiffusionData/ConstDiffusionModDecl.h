#ifndef CONSTDIFFUSION_MOD_DELC_H 
#define CONSTDIFFUSION_MOD_DELC_H 
#include <cmath> 
extern "C" { 
double ConstDiffusionVol1xMaxP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity1xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC1xSerDirichlet_Xlower_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC1xSerDirichlet_Xupper_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC1xSerNeumann_Xlower_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC1xSerNeumann_Xupper_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);

double ConstDiffusionVol1xMaxP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC1xSerDirichlet_Xlower_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC1xSerDirichlet_Xupper_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC1xSerNeumann_Xlower_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC1xSerNeumann_Xupper_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);

double ConstDiffusionVol1xMaxP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC1xSerDirichlet_Xlower_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC1xSerDirichlet_Xupper_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC1xSerNeumann_Xlower_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC1xSerNeumann_Xupper_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);


double ConstDiffusionVol2xMaxP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity2xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC2xSerDirichlet_Xlower_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerDirichlet_Xupper_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Xlower_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Xupper_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionSurf2xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity2xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC2xSerDirichlet_Ylower_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerDirichlet_Yupper_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Ylower_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Yupper_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);

double ConstDiffusionVol2xMaxP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC2xSerDirichlet_Xlower_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerDirichlet_Xupper_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Xlower_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Xupper_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionSurf2xMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC2xSerDirichlet_Ylower_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerDirichlet_Yupper_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Ylower_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Yupper_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);

double ConstDiffusionVol2xMaxP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC2xSerDirichlet_Xlower_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerDirichlet_Xupper_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Xlower_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Xupper_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionSurf2xMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC2xSerDirichlet_Ylower_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerDirichlet_Yupper_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Ylower_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC2xSerNeumann_Yupper_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);


double ConstDiffusionVol3xMaxP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC3xSerDirichlet_Xlower_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerDirichlet_Xupper_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Xlower_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Xupper_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionSurf3xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC3xSerDirichlet_Ylower_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerDirichlet_Yupper_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Ylower_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Yupper_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionSurf3xMax_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xMax_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC3xSerDirichlet_Zlower_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerDirichlet_Zupper_P1(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Zlower_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Zupper_P1(const double dx, const double *fSkin, const double fpBC, double *fGhost);

double ConstDiffusionVol3xMaxP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC3xSerDirichlet_Xlower_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerDirichlet_Xupper_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Xlower_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Xupper_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionSurf3xMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC3xSerDirichlet_Ylower_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerDirichlet_Yupper_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Ylower_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Yupper_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionSurf3xMax_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC3xSerDirichlet_Zlower_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerDirichlet_Zupper_P2(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Zlower_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Zupper_P2(const double dx, const double *fSkin, const double fpBC, double *fGhost);

double ConstDiffusionVol3xMaxP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC3xSerDirichlet_Xlower_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerDirichlet_Xupper_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Xlower_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Xupper_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionSurf3xMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC3xSerDirichlet_Ylower_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerDirichlet_Yupper_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Ylower_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Yupper_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionSurf3xMax_Z_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xMax_Z_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBC3xSerDirichlet_Zlower_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerDirichlet_Zupper_P3(const double dx, const double *fSkin, const double fBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Zlower_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);
void ConstDiffusionBC3xSerNeumann_Zupper_P3(const double dx, const double *fSkin, const double fpBC, double *fGhost);



 
double ConstDiffusionVol1xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity1xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol1xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol1xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol1xSerP4(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf1xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf1xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 


double ConstDiffusionVol2xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity2xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity2xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol2xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol2xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xSer_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol2xSerP4(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf2xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf2xSer_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf2xSer_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 


double ConstDiffusionVol3xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurfPositivity3xSer_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol3xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol3xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Z_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Z_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 

double ConstDiffusionVol3xSerP4(const double *w, const double *dxv, const double *nu, const double *f, double *out); 
void ConstDiffusionSurf3xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_X_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Y_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionSurf3xSer_Z_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 
void ConstDiffusionBoundarySurf3xSer_Z_P4(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr); 



 
} 
#endif 
