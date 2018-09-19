-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing the Lagrange
-- multiplier fix for a distribution function
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.lagrangeFixData._VlasovLagrangeFixCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select kernel function
function _M.selectLagrangeFixFunction(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("VlasovLagrangeFix%s%dx%dv%dp",
				basisNmMap[basisNm], CDIM, VDIM, polyOrder)
   return ffi.C[funcNm]
end

return _M
