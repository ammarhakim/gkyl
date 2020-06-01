-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for device kernels used by the moment calculation updater.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi      = require "ffi"
local RectCart = require "Grid.RectCart"
local Range    = require "Lib.Range"

ffi.cdef [[

void cuda_MomentCalc1x1vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out);
void cuda_MomentCalc1x1vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out);
void cuda_MomentCalc1x1vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out);
void cuda_MomentCalc1x1vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out);
 
]]
