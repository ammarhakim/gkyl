-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing the twit-shift modification to the phi smoother BCs.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

--select Kernel to use for modifying the BC
function _M.selectBcModifier(edge, basisNm, CDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("skin_ghost_avg_%s_%dx_%s_p%d",edge, CDIM,basisNmMap[basisNm],polyOrder)
   local funcSign = "(double *skinField, double *ghostField)"
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
