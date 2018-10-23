-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into binary operator calculation C++ kernel functions based on
-- CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.binOpCalcData._CartFieldBinOpCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- select kernel function to compute specified operation
-- between inputs with same dimensionality.
function _M.selectBinOpCalcS(op, basisNm, CDIM, polyOrder)
   local funcNm = string.format("CartFieldBinOp%s%dx%s_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- select kernel function to compute specified operation
-- between inputs of different dimensionality.
function _M.selectBinOpCalcD(op, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("CartFieldBinOp%s%dx%dv%s_P%d", op, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
