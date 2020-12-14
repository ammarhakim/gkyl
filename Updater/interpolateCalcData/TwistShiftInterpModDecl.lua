-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into TwistShiftInterp C++ kernel functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select twist-shift interpolation operator kernel.
function _M.selectTwistShiftInterp(basisNm, dim, polyOrder)
   local funcType = "void"
   local funcNm = string.format("TwistShiftInterp%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *xLimLo, const double *xLimUp, const double *fldSrc, double *fldDest)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectTwistShiftInterpGen(basisNm, dim, polyOrder)
   local funcType = "void"
   local funcNm = string.format("TwistShiftInterp_limXvarYfixed%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double *fldSrc, double *fldDest)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectTwistShiftInterpGenSub(basisNm, dim, polyOrder)
   local funcType = "void"
   local funcNm   = string.format("TwistShiftInterp_xLimDG%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyLim, const double ycLim, const double *fldSrc, double *fldDest)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectTwistShiftInterpGenSubY(basisNm, dim, polyOrder)
   local funcType = "void"
   local funcNm   = string.format("TwistShiftInterp_yLimDG%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dxLim, const double xcLim, const double *fldSrc, double *fldDest)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectTwistShiftInterpGenSubMinusTwoCorners(basisNm, dim, polyOrder)
   local funcType = "void"
   local funcNm   = string.format("TwistShiftInterp_mTwoCorners%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *xLimLoL, const double *xLimUpL, const double yLimLoL, const double yLimUpL, const double dyLimL, const double ycLimL, const double *xLimLoR, const double *xLimUpR, const double yLimLoR, const double yLimUpR, const double dyLimR, const double ycLimR, const double *fldSrc, double *fldDest)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M