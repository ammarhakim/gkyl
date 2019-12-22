local ffi = require "ffi"

ffi.cdef [[

void calcMom1x1vSer_M0_P1(RectCart_t *grid, Range_t *range, GkDeviceProp *prop, int numBlocks, int numThreads, double *fIn, double *out)
 
]]
