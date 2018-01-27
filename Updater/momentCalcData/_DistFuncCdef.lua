-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the moment calc
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[
void MomentCalc1x1vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x1vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x1vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x1vSer_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vSer_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vSer_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x1vSer_P4(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x2vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x2vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x2vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x2vSer_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vSer_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vSer_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x2vSer_P4(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x3vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x3vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x3vSer_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vSer_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vSer_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x3vSer_P4(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x2vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x2vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x2vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x2vSer_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vSer_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vSer_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x2vSer_P4(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x3vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x3vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x3vSer_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vSer_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vSer_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x3vSer_P4(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc3x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc3x3vSer_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc3x3vSer_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc3x3vSer_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vSer_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vSer_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc3x3vSer_P4(const double *w, const double *dxv, const double *f, double *out); 



 
void MomentCalc1x1vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x1vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x1vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x1vMax_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x1vMax_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x1vMax_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x1vMax_P4(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x2vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x2vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x2vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x2vMax_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x2vMax_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x2vMax_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x2vMax_P4(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc1x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x3vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x3vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x3vMax_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc1x3vMax_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc1x3vMax_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc1x3vMax_P4(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x2vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x2vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x2vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x2vMax_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x2vMax_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x2vMax_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x2vMax_P4(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc2x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x3vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x3vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x3vMax_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc2x3vMax_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc2x3vMax_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc2x3vMax_P4(const double *w, const double *dxv, const double *f, double *out); 


void MomentCalc3x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc3x3vMax_P1(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc3x3vMax_P2(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vMax_M0_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M1i_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2ij_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2_P3(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M3i_P3(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc3x3vMax_P3(const double *w, const double *dxv, const double *f, double *out); 

void MomentCalc3x3vMax_M0_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M1i_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2ij_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M2_P4(const double *w, const double *dxv, const double *f, double *out); 
void MomentCalc3x3vMax_M3i_P4(const double *w, const double *dxv, const double *f, double *out); 
void IntMomentCalc3x3vMax_P4(const double *w, const double *dxv, const double *f, double *out);  
]]
