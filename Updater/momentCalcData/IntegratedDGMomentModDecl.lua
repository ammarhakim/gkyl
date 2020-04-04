-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into integrated moment calculation C++ kernel functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Updater.momentCalcData._IntegratedDGMomentCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

function _M.selectIntMomCalc(moment, basisNm, DIM, pOrder)
   funcNm = string.format("IntDGMoment%dx%s_%s_P%d", DIM, basisNmMap[basisNm], moment, pOrder)
   return ffi.C[funcNm]
end

function _M.selectIntVspaceMomCalc(moment, basisNm, CDIM, VDIM, pOrder)
   funcNm = string.format("IntDGMoment%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], moment, pOrder)
   return ffi.C[funcNm]
end


return _M
