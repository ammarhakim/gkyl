-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing the Lagrange
-- multiplier fix for a distribution function
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "ser", ["maximal-order"] = "max", ["tensor"] = "tensor" }

local _M = {}

-- Select kernel function.
function _M.selectLagrangeFixFunction(basisNm, CDIM, VDIM, polyOrder, speciesKind)
   local funcType = "void"
   local funcNm   = string.format("lagrangeFix_%s_%dx%dv_%s_p%d", speciesKind, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   if speciesKind == "vlasov" then
      funcSign = "(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f)"
   else
      funcSign = "(double *dm0, double *dm1, double *dm2, double *B, double mass, double *lo, double *L, double *Nv, double *vc, double *f)"
   end
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
