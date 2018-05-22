-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the binary operator calc.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[
void CartFieldBinOpMultiply1xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply1xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpMultiply2xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply2xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpMultiply3xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply3xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpDivide1xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide1xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpDivide2xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide2xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpDivide3xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide3xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpDotProduct1xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out);

void CartFieldBinOpDotProduct1xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out);


void CartFieldBinOpDotProduct2xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out);

void CartFieldBinOpDotProduct2xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out);


void CartFieldBinOpDotProduct3xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out);

void CartFieldBinOpDotProduct3xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out);


void CartFieldBinOpMultiply1x1vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply1x1vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply1x2vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply1x2vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply1x3vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply1x3vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpMultiply2x2vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply2x2vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply2x3vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply2x3vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpMultiply3x3vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply3x3vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpDivide1x1vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide1x1vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide1x2vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide1x2vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide1x3vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide1x3vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpDivide2x2vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide2x2vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide2x3vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide2x3vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpDivide3x3vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide3x3vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 


void CartFieldBinOpMultiply3x2vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpMultiply3x2vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide3x2vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

void CartFieldBinOpDivide3x2vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out); 

]]
