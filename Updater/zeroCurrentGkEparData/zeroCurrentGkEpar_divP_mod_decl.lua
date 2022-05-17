-- Gkyl ------------------------------------------------------------------------
--
-- Select C++ kernels for computing the div{ P } contributions to the field
-- equation for computing the parallel electric field in 1D GK assuming
-- quasineutrality and vanishing currents.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "ser", ["tensor"] = "tensor" }

local _M = {}

-- Select kernel function to compute cell-wise constant nu.
function _M.selectDivPKer(basisNm, DIM, polyOrder)
   local funcType = "void"
   local funcNm   = string.format("zeroCurrentGkEpar_divP_%dx_p%d_%s", DIM, polyOrder, basisNmMap[basisNm])
   local funcSign = "(const double dz, const double charge, const double *delparFac, const double *dlnbmagdz, const double *m2par, const double *m2perp, double *divP)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
