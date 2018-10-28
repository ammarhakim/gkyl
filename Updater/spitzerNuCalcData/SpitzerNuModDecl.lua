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
function _M.selectSpitzerNuCalc(basisNm, CDIM, polyOrder)
   local funcNm = string.format("SpitzerNu%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select kernel function to compute cell-wise constant nu. 
function _M.selectCellAvSpitzerNuCalc(basisNm, CDIM, polyOrder)
   local funcNm = string.format("SpitzerNuCellAv%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
