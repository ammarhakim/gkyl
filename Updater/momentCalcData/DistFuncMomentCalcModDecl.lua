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

-- select function to compute specified gyrokinetic moment
function _M.selectGkMomCalc(mom, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("GkMomentCalc%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], mom, polyOrder)
   return ffi.C[funcNm]
end

-- select function to compute integrated moments
function _M.selectIntMomCalc(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("IntMomentCalc%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
