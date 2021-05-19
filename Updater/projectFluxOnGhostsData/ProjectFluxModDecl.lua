-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for projecting flux onto ghosts.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Kernel function to compute projection of flux onto ghost cells. 
function _M.projectFlux(basisNm, CDIM, VDIM, dir, polyOrder)
   local funcType = "double"
   local funcNm = string.format("ProjectFluxOnGhosts%dx%dvDir%d%s_P%d", CDIM, VDIM, dir, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double wv, const double dv, const double zVal, const double *fIn, double *fHat)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
