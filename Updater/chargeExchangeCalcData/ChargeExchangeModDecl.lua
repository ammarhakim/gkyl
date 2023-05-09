-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing charge exchange quantities.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser" }

local _M = {}

-- Kernel function to compute charge exchange cross section based on fitting function.
function _M.SigmaCX(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = string.format("SigmaCXcellAv%s%dx%dv_P%d", basisNmMap[basisNm], CDIM, VDIM, polyOrder)
   local funcSign = "(const double a, const double b, const double *m0, const double *uIon, const double *uNeut, const double *vtSqIon, double vtSqIonMin, const double *vtSqNeut, double vtSqNeutMin, double *vSigmaCX)"
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
