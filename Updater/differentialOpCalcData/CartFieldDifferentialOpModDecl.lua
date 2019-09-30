-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into binary operator calculation C++ kernel functions based on
-- CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.differentialOpCalcData._CartFieldDifferentialOpCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select volume kernel function to compute specified operation.
function _M.selectOperatorVol(op, basisNm, CDIM, polyOrder)
   local funcNm = string.format("CartFieldDifferentialOp%sVol%dx%s_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select surface kernel functions to compute surface contributions (output is a table of functions).
function _M.selectOperatorSurf(op, basisNm, CDIM, polyOrder)
   if CDIM == 1 then
      local funcNmX = string.format("CartFieldDifferentialOp%sSurf%dx%s_X_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], nullFunc, nullFunc }
   elseif CDIM == 2 then
      local funcNmX = string.format("CartFieldDifferentialOp%sSurf%dx%s_X_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("CartFieldDifferentialOp%sSurf%dx%s_Y_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], nullFunc }
   elseif CDIM == 3 then
      local funcNmX = string.format("CartFieldDifferentialOp%sSurf%dx%s_X_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("CartFieldDifferentialOp%sSurf%dx%s_Y_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("CartFieldDifferentialOp%sSurf%dx%s_Z_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end

-- Select boundary suface function kernels to compute boundary surface terms (output is a table of functions).
function _M.selectOperatorBoundarySurf(op, basisNm, CDIM, polyOrder)
   if CDIM == 1 then
      local funcNmX = string.format("CartFieldDifferentialOp%sBoundarySurf%dx%s_X_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], nullFunc, nullFunc }
   elseif CDIM == 2 then
      local funcNmX = string.format("CartFieldDifferentialOp%sBoundarySurf%dx%s_X_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("CartFieldDifferentialOp%sBoundarySurf%dx%s_Y_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], nullFunc }
   elseif CDIM == 3 then
      local funcNmX = string.format("CartFieldDifferentialOp%sBoundarySurf%dx%s_X_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("CartFieldDifferentialOp%sBoundarySurf%dx%s_Y_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("CartFieldDifferentialOp%sBoundarySurf%dx%s_Z_P%d", op, CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end

return _M
