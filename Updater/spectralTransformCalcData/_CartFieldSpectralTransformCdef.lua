-- Gkyl -------------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the CartFieldSpectralTransform updater.
--    _______     ___
-- + 6 @ |||| # P ||| +
---------------------------------------------------------------------------------------

local ffi = require "ffi" 

ffi.cdef [[

typedef struct spectralTransform spectralTransform;
spectralTransform* new_spectralTransform(const int nModes, const int nCells, const int pOrder, const int nSurfB);
void assignLHSMatrixSer(spectralTransform *sTransObj, const int pOrder, const int cellIdx, const int spectralIdx, const double *spectralBasisIn);
void getLHSMatrixInverse(spectralTransform *sTransObj);
void assignRHSMatrix1x1vSer_P1OpDir1(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
void assignRHSMatrix1x1vSer_P1OpDir2(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
void assignRHSMatrix1x1vSer_P2OpDir1(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
void assignRHSMatrix1x1vSer_P2OpDir2(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
void solveTransform(spectralTransform *sTransObj);
void getSolution1x1vSer_P1(spectralTransform *sTransObj, const int cellIdx, double *fSpectral);
void getSolution1x1vSer_P2(spectralTransform *sTransObj, const int cellIdx, double *fSpectral);

]]
