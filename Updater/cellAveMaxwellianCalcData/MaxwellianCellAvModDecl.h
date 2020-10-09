#ifndef MAXWELLIANCELLAV_MOD_DELC_H 
#define MAXWELLIANCELLAV_MOD_DELC_H 
#include <cmath> 

#include <algorithm> 

extern "C" { 
void MaxwellianCellAvSer1x1v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvSer1x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvSer1x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvSer1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 

void GkMaxwellianCellAvSer1x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvSer2x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvSer2x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvSer2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvSer3x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvSer3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvSer1x1v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvSer1x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvSer1x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvSer1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 

void GkMaxwellianCellAvSer1x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvSer2x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvSer2x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvSer2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvSer3x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvSer3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvSer1x1v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvSer1x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvSer1x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvSer1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 

void GkMaxwellianCellAvSer1x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvSer2x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvSer2x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvSer2x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvSer3x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvSer3x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvMax1x1v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvMax1x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvMax1x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvMax1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 

void GkMaxwellianCellAvMax1x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvMax2x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvMax2x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvMax2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvMax3x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvMax3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvMax1x1v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvMax1x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvMax1x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvMax1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 

void GkMaxwellianCellAvMax1x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvMax2x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvMax2x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvMax2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvMax3x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvMax3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvMax1x1v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvMax1x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvMax1x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvMax1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 

void GkMaxwellianCellAvMax1x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvMax2x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void MaxwellianCellAvMax2x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvMax2x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 


void MaxwellianCellAvMax3x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax); 

void GkMaxwellianCellAvMax3x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax); 



 
} 
#endif 
