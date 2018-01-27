#ifndef MAXWELL_MOD_DELC_H 
#define MAXWELL_MOD_DELC_H 
extern "C" { 
typedef struct { double c, chi, gamma; } MaxwellEq_t; 
 
void MaxwellVol1xMaxP0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xMax_X_P0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol1xMaxP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xMax_X_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol1xMaxP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xMax_X_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol1xMaxP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xMax_X_P3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol1xMaxP4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xMax_X_P4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 


void MaxwellVol2xMaxP0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xMax_X_P0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xMax_Y_P0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol2xMaxP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xMax_X_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xMax_Y_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol2xMaxP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xMax_X_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xMax_Y_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol2xMaxP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xMax_X_P3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xMax_Y_P3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol2xMaxP4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xMax_X_P4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xMax_Y_P4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 


void MaxwellVol3xMaxP0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xMax_X_P0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Y_P0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Z_P0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol3xMaxP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xMax_X_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Y_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Z_P1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol3xMaxP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xMax_X_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Y_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Z_P2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol3xMaxP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xMax_X_P3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Y_P3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Z_P3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol3xMaxP4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xMax_X_P4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Y_P4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xMax_Z_P4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 



 
void MaxwellVol1xSerP0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xSer_X_P0(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol1xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xSer_X_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol1xSerP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xSer_X_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol1xSerP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xSer_X_P3(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol1xSerP4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf1xSer_X_P4(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 


void MaxwellVol2xSerP0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xSer_X_P0(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xSer_Y_P0(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol2xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xSer_X_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xSer_Y_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol2xSerP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xSer_X_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xSer_Y_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol2xSerP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xSer_X_P3(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xSer_Y_P3(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol2xSerP4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf2xSer_X_P4(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf2xSer_Y_P4(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 


void MaxwellVol3xSerP0(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xSer_X_P0(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Y_P0(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Z_P0(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol3xSerP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xSer_X_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Y_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Z_P1(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol3xSerP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xSer_X_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Y_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Z_P2(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol3xSerP3(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xSer_X_P3(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Y_P3(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Z_P3(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 

void MaxwellVol3xSerP4(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out); 
void MaxwellSurf3xSer_X_P4(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Y_P4(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 
void MaxwellSurf3xSer_Z_P4(const MaxwellEq_t * meq, const double *w, const double *dx, const double *ql, const double *qr, double *outl, double *outr); 



 
} 
#endif 
