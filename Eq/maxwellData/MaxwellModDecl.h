#ifndef MAXWELL_MOD_DELC_H 
#define MAXWELL_MOD_DELC_H 
extern "C" { 
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



 
} 
#endif 
