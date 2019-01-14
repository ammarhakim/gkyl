-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing the spitzer collisionality.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.spitzerNuCalcData._SpitzerNuCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Select kernel function to compute expansion coefficients of nu. 
function _M.selectSpitzerNuScale(basisNm, CDIM, polyOrder)
   local funcNm = string.format("SpitzerNuScale%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select kernel function to compute cell-wise constant nu. 
function _M.selectCellAvSpitzerNuScale(basisNm, CDIM, polyOrder)
   local funcNm = string.format("SpitzerNuCellAvScale%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select kernel function to build expansion coefficients of nu from Spitzer formula. 
function _M.selectSpitzerNuBuild(basisNm, CDIM, polyOrder)
   local funcNm = string.format("SpitzerNuBuild%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select kernel function to compute cell-wise constant nu from Spitzer formula. 
function _M.selectCellAvSpitzerNuBuild(basisNm, CDIM, polyOrder)
   local funcNm = string.format("SpitzerNuCellAvBuild%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
