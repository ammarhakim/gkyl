-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into CartFieldInterpolate C++ kernel functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Updater.interpolateCalcData._CartFieldInterpolateCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select restriction operator kernel.
function _M.selectRestriction(basisNm, dim, polyOrder)
   local funcNm = string.format("CartFieldInterpRestrict%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select prolongation operator kernel.
function _M.selectProlongation(basisNm, dim, polyOrder)
   local funcNm = string.format("CartFieldInterpProlong%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end


return _M
