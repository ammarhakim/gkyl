-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into binary operator calculation C++ kernel functions based on
-- CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Updater.binOpCalcData._CartFieldBinOpCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select kernel function to compute specified operation
-- between inputs with same dimensionality.
function _M.selectBinOpCalcS(op, basisNm, CDIM, VDIM, polyOrder, applyPos, inclVx2)
   local posString = ""
   if (applyPos and op=="Divide") then posString = "Positivity" end
   local funcNm
   if VDIM then 
      local inclVx2str = ""
      if inclVx2 and polyOrder==1 then inclVx2str = "InclVx2" end
      funcNm = string.format("CartFieldBinOp%s%s%dx%dv%s_P%d", op, posString, CDIM, VDIM, basisNmMap[basisNm], inclVx2str, polyOrder)
   else
      funcNm = string.format("CartFieldBinOp%s%s%dx%s_P%d", op, posString, CDIM, basisNmMap[basisNm], polyOrder)
   end
   return ffi.C[funcNm]
end

-- Select kernel function to compute specified operation
-- between inputs of different dimensionality.
function _M.selectBinOpCalcD(op, basisNm, CDIM, VDIM, polyOrder, inclVx2)
   local inclVx2str = ""
   if inclVx2 and polyOrder==1 then inclVx2str = "InclVx2" end
   local funcNm = string.format("CartFieldBinOpConfPhase%s%dx%dv%s_P%d", op, CDIM, VDIM, basisNmMap[basisNm], inclVx2str, polyOrder)
   return ffi.C[funcNm]
end

return _M
