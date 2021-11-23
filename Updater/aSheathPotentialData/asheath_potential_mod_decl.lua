-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernels for computing the potential assuming
-- ambipolar sheath flow and adiabatic electrons.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "ser", ["tensor"] = "tensor" }

local _M = {}

-- Select kernel function to compute expansion coefficients of nu. 
function _M.ModDecl.selectPhiSheathQuad(basisNm, CDIM, polyOrder, quadType)
   local funcType = "void"
   local funcSign = "(const double q_e, const double m_e, const double T_e, const double *Gamma_i, const double *m0Ion, double *m0IonS, double *phiS)"
   local funcNm = {}
   for i, b in ipairs({"lower","upper"}) do
      funcNm[i] = string.format("asheath_potential_gauss_phi_sheath_%s_%dx_p%d_%s", b, CDIM, polyOrder, basisNmMap[basisNm])
   end
   local CDefStr = ""
   for i = 1, #funcNm do CDefStr = CDefStr .. (funcType .. " " .. funcNm[i] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kerOut = {}
   for i, b in ipairs({"lower","upper"}) do
      local tmp = ffi.C[funcNm[i]]
      kerOut[b] = tmp
   end
   return kerOut 
end

-- Select kernel function to compute cell-wise constant nu. 
function _M.self._phiKer(basisNm, CDIM, polyOrder, quadType)
   local funcType = "void"
   local funcNm = string.format("asheath_potential_gauss_phi_%dx_p%d_%s", CDIM, polyOrder, basisNmMap[basisNm])
   local funcSign = "(const double q_e, const double T_e, const double *m0Ion, const double *m0IonS, const double *phiS, double *phi)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
