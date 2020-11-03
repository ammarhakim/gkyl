-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for projecting a Maxwellian onto DG basis.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Select kernel functions to evaluate configuration space quantities at conf-space 
-- ordinates, and perform the phase space quadrature.
function _M.selectQuad(basisNm, cDim, vDim, uDimIn, polyOrder, quadType, isGK)
   local funcType = "void"

   local maxwellKernels = {}
   local GKstr, strU
   if isGK then
      GKstr = "Gk"
      if uDimIn == 1 then
         strU = ""
      else
         strU = "Uz"
      end
   else
      GKstr = ""
      if vDim == uDimIn then
         strU = ""
      else
         strU = "Upar"
      end
   end
   local funcNm = string.format("%sMaxwellianOnBasis%s%dx%dv%s%s_P%d_evAtConfOrd", GKstr, quadType, cDim, vDim, basisNmMap[basisNm], strU, polyOrder)
   local funcSign = "(const double *den, const double *flowU, const double *vtSq, double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd)"
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   local tmpKrn = ffi.C[funcNm]
   maxwellKernels[1] = tmpKrn

   local funcNm   = string.format("%sMaxwellianOnBasis%s%dx%dv%s_P%d_phaseQuad", GKstr, quadType, cDim, vDim, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dvx, double *fMOut)"
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   local tmpKrn = ffi.C[funcNm]
   maxwellKernels[2] = tmpKrn

   return maxwellKernels
end

function _M.selectInnerLoop(isGK)
   local funcType = "void"
   local funcNm
   if isGK then
      funcNm = "GkMaxwellianInnerLoop"
   else
      funcNm = "MaxwellianInnerLoop"
   end
   local funcSign = "(double *n, double *u, double *vtSq, double *bmag, double m_, double *fItr, double *weights, double *dz, double *zc, double *ordinates, double *basisAtOrdinates, double *phaseToConfOrdMap, int numPhaseBasis, int numConfOrds, int numPhaseOrds, int numConfDims, int numPhaseDims)"
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
