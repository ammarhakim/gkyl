-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing charge exchange quantities.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.chargeExchangeCalcData._RelativeVelocityCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Kernel function to compute relative velocity for Pauls CX model. 
function _M.GkVrelProdCX(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("GkProdCXcellAv%s%dx%dv_P%d", basisNmMap[basisNm], CDIM, VDIM, polyOrder)
   return ffi.C[funcNm]
end

-- Kernel function to compute relative velocity for Pauls CX model. 
function _M.VmVrelProdCX(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("VmProdCXcellAv%s%dx%dv_P%d", basisNmMap[basisNm], CDIM, VDIM, polyOrder)
   return ffi.C[funcNm]
end

return _M
