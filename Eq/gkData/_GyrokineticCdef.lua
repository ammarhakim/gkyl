local ffi = require "ffi" 

ffi.cdef [[
double GyrokineticVol1x1vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *f, double *out); 
double EmGyrokineticVol1x1vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *f, double *out); 
double dAdtVol1x1vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *dApardt, const double *f, double *out); 
double GyrokineticSurf1x1vSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double EmGyrokineticSurf1x1vSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf1x1vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double EmGyrokineticSurf1x1vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr); 
double dAdtSurf1x1vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr); 
double calcSheathDeltaPhi1xSer_P1(const double *phi, const double *phiWall, const double zVal);
void calcSheathPartialReflection1x1vSer_P1(const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);


double GyrokineticVol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *f, double *out); 
double EmGyrokineticVol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *f, double *out); 
double dAdtVol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *dApardt, const double *f, double *out); 
double GyrokineticSurf2x2vSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double EmGyrokineticSurf2x2vSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf2x2vSer_Y_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double EmGyrokineticSurf2x2vSer_Y_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf2x2vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double EmGyrokineticSurf2x2vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr); 
double dAdtSurf2x2vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr); 
double calcSheathDeltaPhi2xSer_P1(const double *phi, const double *phiWall, const double zVal);
void calcSheathPartialReflection2x2vSer_P1(const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);


double GyrokineticVol3x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *f, double *out); 
double EmGyrokineticVol3x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *f, double *out); 
double dAdtVol3x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *dApardt, const double *f, double *out); 
double GyrokineticSurf3x2vSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double EmGyrokineticSurf3x2vSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf3x2vSer_Y_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double EmGyrokineticSurf3x2vSer_Y_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf3x2vSer_Z_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double EmGyrokineticSurf3x2vSer_Z_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr); 
double GyrokineticSurf3x2vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr); 
double EmGyrokineticSurf3x2vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr); 
double dAdtSurf3x2vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr); 
double calcSheathDeltaPhi3xSer_P1(const double *phi, const double *phiWall, const double zVal);
void calcSheathPartialReflection3x2vSer_P1(const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);


]]
