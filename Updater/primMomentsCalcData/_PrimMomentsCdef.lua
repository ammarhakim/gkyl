-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in calculating primitive moments.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

void SelfPrimMoments1x1vSer_P1(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x1vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x1vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x1vSer_P2(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x1vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x1vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x1vSer_P3(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x1vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x1vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x1vSer_P4(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x1vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x1vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x2vSer_P1(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x2vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x2vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x2vSer_P2(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x2vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x2vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x2vSer_P3(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x2vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x2vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x2vSer_P4(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x2vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x2vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x3vSer_P1(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x3vSer_P2(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x3vSer_P3(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments1x3vSer_P4(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_1x3vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_1x3vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 


void SelfPrimMoments2x2vSer_P1(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_2x2vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x2vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x2vSer_P2(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_2x2vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x2vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x2vSer_P3(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_2x2vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x2vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x2vSer_P4(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_2x2vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x2vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x3vSer_P1(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_2x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x3vSer_P2(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_2x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x3vSer_P3(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_2x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments2x3vSer_P4(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_2x3vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_2x3vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 


void SelfPrimMoments3x3vSer_P1(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_3x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_3x3vSer_P1(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments3x3vSer_P2(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_3x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_3x3vSer_P2(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments3x3vSer_P3(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_3x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_3x3vSer_P3(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

void SelfPrimMoments3x3vSer_P4(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq); 

void CrossPrimMoments_VmeiLBO_3x3vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 
void CrossPrimMoments_VmieLBO_3x3vSer_P4(const double mRat, const double *ne, const double *ue, const double *vtSqe, const double *ni, const double *ui, const double *vtSqi, double *uie, double *vtSqie); 

]]

