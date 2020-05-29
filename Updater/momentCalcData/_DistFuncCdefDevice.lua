local ffi = require "ffi"

ffi.cdef [[

void cuda_MomentCalc1x1vSer_M0_P1(GkylRectCart_t *grid, GkylRange_t *pRange, GkylRange_t *cRange, GkDeviceProp *prop, int numBlocks, int numThreads, const double *fIn, double *out);
 
]]
