-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing primitive moments u and vtSq
-- for cross collisions based on CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.primMomentsCalcData._PrimMomentsCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select kernel function to compute the primitive moments for self-collision terms. 
function _M.selectSelfPrimMomentsCalc(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("SelfPrimMoments%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- select kernel function to compute the cross primitive moments. 
function _M.selectCrossPrimMomentsCalc(collOp, collide, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("CrossPrimMoments_%s%s_%dx%dv%s_P%d", collide, operator, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
