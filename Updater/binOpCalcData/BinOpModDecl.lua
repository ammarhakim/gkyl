-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into binary operator calculation C++ kernel functions based on
-- CDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.binOpCalcData._BinOpCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select function to compute specified moment
function _M.selectBinOpCalc(op, basisNm, CDIM, polyOrder)
   local funcNm = string.format("BinOp%s%dx%s_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
