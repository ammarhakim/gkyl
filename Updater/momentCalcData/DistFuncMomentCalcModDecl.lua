-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into moment calculation C++ kernel functions based on
-- CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.momentCalcData._DistFuncCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select function to compute specified moment
function _M.selectMomCalc(mom, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("MomentCalc%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], mom, polyOrder)
   return ffi.C[funcNm]
end

-- select function to compute integrated moments
function _M.selectIntMomCalc(basisNm, VDIM, polyOrder)

end

return _M
