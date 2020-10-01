-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing the Voronov reaction rate.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Kernel function to compute Voronov reaction rate. 
function _M.voronovCoef(basisNm, CDIM, polyOrder)
   local funcType = "double"
   local funcNm = string.format("VoronovReactRateCellAv%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double elemCharge, const double m_, const double *m0, const double *vtSqNeut, const double *vtSqElc, const double E, const double A, const double K, const double P, const double X, double *coefIz)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.ionizationTemp(basisNm, CDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("IonizationTemp%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
