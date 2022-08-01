-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for projecting a BiMaxwellian onto DG basis.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["tensor"] = "Tensor" }

local _M = {}

-- Select kernel functions to evaluate configuration space quantities at conf-space 
-- ordinates, and perform the phase space quadrature.
function _M.selectQuad(basisNm, cDim, vDim, polyOrder, quadType)
   local funcType = "void"

   local maxwellKernels = {}
   local funcNm = string.format("GkBiMaxwellianOnBasis%s%dx%dv%s_P%d_evAtConfOrd", quadType, cDim, vDim, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *den, const double *flowU, const double *vtparSq, const double *vtperpSq, const double *bmag, double *flowUOrd, double *vtparSqOrd, double *vtperpSqOrd, double *fMFacOrd, double *bmagOrd)"

   local CDefStr = ""
   CDefStr = CDefStr .. (funcType .. " " .. funcNm .. funcSign .. ";\n")
   ffi.cdef(CDefStr)

   local tmpKrn = ffi.C[funcNm]
   maxwellKernels[1] = tmpKrn

   local funcNm   = string.format("GkBiMaxwellianOnBasis%s%dx%dv%s_P%d_phaseQuad", quadType, cDim, vDim, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *flowUOrd, const double *vtparSqOrd, const double *vtperpSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dvx, double *fMOut)"
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   local tmpKrn = ffi.C[funcNm]
   maxwellKernels[2] = tmpKrn

   return maxwellKernels
end

return _M
