-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into CartFieldInterpolate C++ kernel functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Updater.interpolateCalcData._CartFieldInterpolateCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select interpolation operator kernel.
function _M.selectInterpolation(basisNm, dim, polyOrder)
   local funcNm = string.format("CartFieldInterp%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select interpolation kernel for the cases that don't interpolate in one direction.
function _M.selectInterpolationExcept1dir(basisNm, dim, polyOrder, sameDir)
   local pDirs  = {"X","Y","Z","VX","VY","VZ"}
   local funcNm = string.format("CartFieldInterp%dx%s_%s_P%d", dim, basisNmMap[basisNm], pDirs[sameDir], polyOrder)
   return ffi.C[funcNm]
end


return _M
