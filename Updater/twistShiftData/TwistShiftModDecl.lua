-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into TwistShiftInterp C++ kernel functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- C Struct needed by all kernels/C functions.
ffi.cdef[[
  typedef struct tsStruct tsStruct;
]]

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select function that allocates the vectors and (temporary) matrix used by the twist-shift operation.
function _M.selectTwistShiftAlloc()
   local funcType = "tsStruct*"
   local funcNm   = "twistShift_alloc" 
   local funcSign = "(const int numMats, const int matN)" 
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select function that allocates matrices that multiply donor field DG coefficients in twist-shift operation.
function _M.selectTwistShiftAllocCellMat()
   local funcType = "void" 
   local funcNm   = "twistShift_allocCellMat" 
   local funcSign = "(tsStruct *tsMatVecs, const int cellIdx, const int numMats)" 
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select twist-shift interpolation operator kernel.
function _M.selectTwistShiftIntSubX(dim, basisNm, polyOrder, yShPolyOrder)
   local funcType = "void"
   local funcNm   = string.format("twistShift_xLimDG%dx%sP%d_yShP%d", dim, basisNmMap[basisNm], polyOrder, yShPolyOrder)
   local funcSign = "(const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const int matIdx)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectTwistShiftIntSubY(dim, basisNm, polyOrder, yShPolyOrder)
   local funcType = "void"
   local funcNm   = string.format("twistShift_yLimDG%dx%sP%d_yShP%d", dim, basisNmMap[basisNm], polyOrder, yShPolyOrder)
   local funcSign = "(const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const int matIdx)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectTwistShiftIntSubM2Corners(dim, basisNm, polyOrder, yShPolyOrder)
   local funcType = "void"
   local funcNm   = string.format("twistShift_M2Corners%dx%sP%d_yShP%d", dim, basisNmMap[basisNm], polyOrder, yShPolyOrder)
   local funcSign = "(const double *xLimLoL, const double *xLimUpL, const double yLimLoL, const double yLimUpL, const double *xLimLoR, const double *xLimUpR, const double yLimLoR, const double yLimUpR, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const int matIdx)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select kernel that performs a mat-vec multiply to compute the contribution from a donor cell to a target cell.
function _M.selectTwistShiftMatVecMult(dim, basisNm, polyOrder)
   local funcType = "void"
   local funcNm   = string.format("twistShift_matVecMult%dx%sP%d", dim, basisNmMap[basisNm], polyOrder, yShPolyOrder)
   local funcSign = "(tsStruct *tsData, const int xIdx, const int matIdx, const double *fldDo, double *fldTar)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
