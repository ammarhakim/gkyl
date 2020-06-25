-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into CartFieldInterpolate C++ kernel functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select interpolation operator kernel.
function _M.selectInterpolation(basisNm, dim, polyOrder)
   local funcType = "void"
   local funcNm = string.format("CartFieldInterp%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select interpolation kernel for the cases that don't interpolate in one direction.
function _M.selectInterpolationExcept1dir(basisNm, dim, polyOrder, sameDir)
   local pDirs  = {"X","Y","Z","VX","VY","VZ"}
   local funcType = "void"
   local funcNm = string.format("CartFieldInterp%dx%s_%s_P%d", dim, basisNmMap[basisNm], pDirs[sameDir], polyOrder)
   local funcSign = "(const double *wC, const double *wF, const double *dxC, const double *dxF, const double *fldC, double *fldF)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end


return _M
