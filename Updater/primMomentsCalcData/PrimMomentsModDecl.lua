-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing primitive moments u and vtSq
-- for cross collisions based on CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Select kernel function to compute the primitive moments for self-collision terms. 
function _M.selectSelfPrimMomentsCalc(kinSpecies, basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("%sSelfPrimMoments%dx%dv%s_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign
   if polyOrder==1 then
      funcSign = "(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq)"
   else
      funcSign = "(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq)"
   end

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select kernel function to compute the cross-collision primitive moments. 
function _M.selectCrossPrimMomentsCalc(operator, basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("%sCrossPrimMoments%dx%dv%s_P%d", operator, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign
   if (operator=="VmLBO") or (operator=="GkLBO") then
      if polyOrder==1 then
         funcSign = [[(binOpData_t *data, binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther)]]
      else
         funcSign = [[(binOpData_t *data, binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther)]]
      end
   elseif (operator=="VmBGK") or (operator=="GkBGK") then
      funcSign = [[(binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *m0Other, const double *uOther, const double *vtSqOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther)]]
   end

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
