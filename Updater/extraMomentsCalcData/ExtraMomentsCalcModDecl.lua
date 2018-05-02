-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Extra moments calculation C++ kernel functions based on
-- CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.extraMomentsCalcData._ExtraMomentsCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select function to compute specified moment
function _M.selectMomCalc(mom, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("ExtraMomentsCalc%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], mom, polyOrder)
   return ffi.C[funcNm]
end

-- select function to compute specified gyrokinetic moment
function _M.selectGkMomCalc(mom, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("GkExtraMomentsCalc%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], string.sub(mom,3), polyOrder)
   return ffi.C[funcNm]
end

return _M
