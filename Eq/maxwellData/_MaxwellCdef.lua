-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the Maxwell equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

typedef struct { double c, chi, gamma; } MaxwellEq_t; 
 
void MaxwellVol1xMaxP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xMax_X_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol1xMaxP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xMax_X_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 


void MaxwellVol2xMaxP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xMax_X_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xMax_Y_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol2xMaxP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xMax_X_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xMax_Y_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 


void MaxwellVol3xMaxP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xMax_X_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Y_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Z_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol3xMaxP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xMax_X_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Y_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Z_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 



 
void MaxwellVol1xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xSer_X_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol1xSerP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xSer_X_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 


void MaxwellVol2xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xSer_X_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xSer_Y_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol2xSerP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xSer_X_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xSer_Y_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 


void MaxwellVol3xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xSer_X_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Y_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Z_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol3xSerP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xSer_X_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Y_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Z_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr);

]]
