-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the moment calc
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[
void accumulateConfToPhase1x1vSer_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x1vSer_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x1vMax_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x1vMax_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x1vSer_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x1vSer_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x1vMax_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x1vMax_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x1vSer_P3(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x1vSer_P3(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x1vMax_P3(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x1vMax_P3(const double fact, const double *fconf, double *fphase); 


void accumulateConfToPhase1x2vSer_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x2vSer_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x2vMax_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x2vMax_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x2vSer_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x2vSer_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x2vMax_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x2vMax_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x2vSer_P3(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x2vSer_P3(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x2vMax_P3(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x2vMax_P3(const double fact, const double *fconf, double *fphase); 


void accumulateConfToPhase1x3vSer_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x3vSer_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x3vMax_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x3vMax_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x3vSer_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x3vSer_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x3vMax_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x3vMax_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x3vSer_P3(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x3vSer_P3(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase1x3vMax_P3(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase1x3vMax_P3(const double fact, const double *fconf, double *fphase); 


void accumulateConfToPhase2x2vSer_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x2vSer_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase2x2vMax_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x2vMax_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase2x2vSer_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x2vSer_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase2x2vMax_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x2vMax_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase2x2vSer_P3(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x2vSer_P3(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase2x2vMax_P3(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x2vMax_P3(const double fact, const double *fconf, double *fphase); 


void accumulateConfToPhase2x3vSer_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x3vSer_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase2x3vMax_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x3vMax_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase2x3vSer_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x3vSer_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase2x3vMax_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x3vMax_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase2x3vMax_P3(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase2x3vMax_P3(const double fact, const double *fconf, double *fphase); 


void accumulateConfToPhase3x2vSer_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase3x2vSer_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase3x2vMax_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase3x2vMax_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase3x2vSer_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase3x2vSer_P2(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase3x2vMax_P2(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase3x2vMax_P2(const double fact, const double *fconf, double *fphase);  

void accumulateConfToPhase3x2vMax_P3(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase3x2vMax_P3(const double fact, const double *fconf, double *fphase); 


void accumulateConfToPhase3x3vSer_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase3x3vSer_P1(const double fact, const double *fconf, double *fphase); 

void accumulateConfToPhase3x3vMax_P1(const double fact, const double *fconf, double *fphase); 
void assignConfToPhase3x3vMax_P1(const double fact, const double *fconf, double *fphase);
]]
