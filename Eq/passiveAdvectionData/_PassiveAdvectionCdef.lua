local ffi = require "ffi" 

ffi.cdef [[
double PassiveAdvectionVol1xSerP1(const double *w, const double *dxv, const double *f, double *out); 
double PassiveAdvectionSurf1xSer_X1_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol1xSerP2(const double *w, const double *dxv, const double *f, double *out); 
double PassiveAdvectionSurf1xSer_X1_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol2xSerP1(const double *w, const double *dxv, const double *f, double *out); 
double PassiveAdvectionSurf2xSer_X1_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf2xSer_X2_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol2xSerP2(const double *w, const double *dxv, const double *f, double *out); 
double PassiveAdvectionSurf2xSer_X1_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf2xSer_X2_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol3xSerP1(const double *w, const double *dxv, const double *f, double *out); 
double PassiveAdvectionSurf3xSer_X1_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf3xSer_X2_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf3xSer_X3_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol3xSerP2(const double *w, const double *dxv, const double *f, double *out); 
double PassiveAdvectionSurf3xSer_X1_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf3xSer_X2_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf3xSer_X3_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol4xSerP1(const double *w, const double *dxv, const double *f, double *out); 
double PassiveAdvectionSurf4xSer_X1_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf4xSer_X2_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf4xSer_X3_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf4xSer_X4_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol4xSerP2(const double *w, const double *dxv, const double *f, double *out); 
double PassiveAdvectionSurf4xSer_X1_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf4xSer_X2_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf4xSer_X3_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf4xSer_X4_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol5xSerP1(const double *w, const double *dxv, const double *f, double *out); 
double PassiveAdvectionSurf5xSer_X1_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf5xSer_X2_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf5xSer_X3_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf5xSer_X4_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf5xSer_X5_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

double PassiveAdvectionVol5xSerP2(const double *w, const double *dxv, const double *f, double *out); 
double PassiveAdvectionSurf5xSer_X1_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf5xSer_X2_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf5xSer_X3_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf5xSer_X4_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
double PassiveAdvectionSurf5xSer_X5_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

]]
