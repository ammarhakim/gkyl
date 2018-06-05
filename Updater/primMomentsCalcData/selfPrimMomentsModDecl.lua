-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing primitive moments u and vtSq
-- for self collisions based on CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.primMomentsCalcData._VmLBOprimMomentsCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select kernel function to 
function _M.selectSelfPrimMomentsCalc(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("VmLBOconstNuPrimMoments%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
