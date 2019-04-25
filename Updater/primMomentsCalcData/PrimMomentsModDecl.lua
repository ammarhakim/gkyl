-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing primitive moments u and vtSq
-- for cross collisions based on CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.primMomentsCalcData._PrimMomentsCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Select kernel function to compute the primitive moments for self-collision terms. 
function _M.selectSelfPrimMomentsCalc(kinSpecies, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("%sSelfPrimMoments%dx%dv%s_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select kernel function to compute the cross-collision primitive moments. 
function _M.selectCrossPrimMomentsCalc(kinSpecies, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("%sCrossPrimMoments%dx%dv%s_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
