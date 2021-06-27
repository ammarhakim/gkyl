-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for projecting a Maxwellian onto DG basis.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "ser", ["maximal-order"] = "max", ["tensor"] = "tensor" }

local _M = {}

-- Select kernel functions to evaluate configuration space quantities at conf-space
-- ordinates, and perform the phase space quadrature.
function _M.selectQuad(basisNm, dim, polyOrder, quadType)
   local funcType = "void"

   local funcNm = string.format("sqrt_on_basis_%s_%dx_p%d_%s", quadType, dim, polyOrder, basisNmMap[basisNm])
   local funcSign = "(const double qExp, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
