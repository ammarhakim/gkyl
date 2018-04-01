-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions based on
-- CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.confToPhaseData._ConfToPhaseCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select function to compute specified operation
function _M.selectConfToPhase(op, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("%sConfToPhase%dx%dv%s_P%d", op, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
