local ffi = require "ffi"

ffi.cdef [[

void calcMom1x1vSer_M0_P1(unsigned int numBlocks, unsigned int numThreads, unsigned int *nCells, double *w, double *dxv, double *fIn, double *out)
 
]]
