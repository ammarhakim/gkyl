-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for CartFieldSpectralTransform updater.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Updater.spectralTransformCalcData._CartFieldSpectralTransformCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select kernel to initialize Eigen objects. 
function _M.selectNewTransform(op, basisNm, DIM, polyOrder)
   local funcNm = string.format("new_spectralTransform")
   return ffi.C[funcNm]
end

-- Select kernel that assigns the elements of the maxx matrix.
function _M.selectLHSMatrixAssign(op, basisNm, DIM, polyOrder)
   local funcNm = string.format("assignLHSMatrix%s", basisNmMap[basisNm])
   return ffi.C[funcNm]
end

-- Select C++ function that computes the inverse of the mass matrix.
function _M.selectMatrixInverseCalc(op, basisNm, DIM, polyOrder)
   local funcNm = string.format("getLHSMatrixInverse")
   return ffi.C[funcNm]
end

-- Select kernel to assign the right-side vector elements (DG coefficients).
function _M.selectRHSassign(op, basisNm, CDIM, VDIM, polyOrder, dirNum)
   local funcNm = string.format("assignRHSMatrix%dx%dv%s_P%dOpDir%d",CDIM,VDIM,basisNmMap[basisNm],polyOrder,dirNum)
   return ffi.C[funcNm]
end

-- Select C++ function that solves the linear problem. 
function _M.selectTransformSolve(op, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("solveTransform")
   return ffi.C[funcNm]
end

-- Select kernel that redustributes solution from Eigen object to Gkeyll CartField data type.
function _M.selectSolutionReOrg(op, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("getSolution%dx%dv%s_P%d",CDIM,VDIM,basisNmMap[basisNm],polyOrder)
   return ffi.C[funcNm]
end

return _M
