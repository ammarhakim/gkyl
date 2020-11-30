local ffi = require "ffi" 

ffi.cdef [[
double IncompEulerVol1xSerP1(const double q_, const double m_, const double *w, const double *dxv, double *positivityWeightByDir, const double *phi, const double *f, double *out); 
double IncompEulerSurf1xSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double *phi, const double *fl, const double *fr, double *outl, double *outr); 
double IncompEulerSurfPositivity1xSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *phi, const double *fl, const double *fr, double *outl, double *outr); 

double IncompEulerVol2xSerP1(const double q_, const double m_, const double *w, const double *dxv, double *positivityWeightByDir, const double *phi, const double *f, double *out); 
double IncompEulerSurf2xSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double *phi, const double *fl, const double *fr, double *outl, double *outr); 
double IncompEulerSurfPositivity2xSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *phi, const double *fl, const double *fr, double *outl, double *outr); 
double IncompEulerSurf2xSer_Y_P1(const double q_, const double m_, const double *w, const double *dxv, const double *phi, const double *fl, const double *fr, double *outl, double *outr); 
double IncompEulerSurfPositivity2xSer_Y_P1(const double q_, const double m_, const double *w, const double *dxv, const double dtApprox, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *phi, const double *fl, const double *fr, double *outl, double *outr); 

]]
