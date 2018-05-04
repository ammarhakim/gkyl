-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the binary operator calc.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[
void BinOpMultiply1xSer_P1(const double *A, const double *B, const short int vDim, double *out); 

void BinOpMultiply1xSer_P2(const double *A, const double *B, const short int vDim, double *out); 

void BinOpDivide1xSer_P1(const double *A, const double *B, const short int vDim, double *out); 

void BinOpDivide1xSer_P2(const double *A, const double *B, const short int vDim, double *out); 

]]
