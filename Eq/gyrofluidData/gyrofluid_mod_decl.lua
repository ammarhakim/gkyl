-- Gkyl ------------------------------------------------------------------------
--
-- Selection of gyrofluid C/C++ kernels.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "ser", ["maximal-order"] = "max" , ["tensor"] = "tensor"}

local _M = {}

-- Select function to compute volume terms.
function _M.selectVol(basisNm, dim, polyOrder)
   local funcSign = "(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *w, const double *dx, const double *jac, const double *rBmag, const double *jacDbmag, const double *dBdz, const double *sMom, const double *phi, double *primMom, double *out)"

   local funcType = "double"
   local funcNm = string.format("gyrofluid_vol_%dx_p%d_%s", dim, polyOrder, basisNmMap[basisNm])

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

function _M.selectSurf(basisNm, dim, polyOrder)
   local funcType  = "double"
   local funcSign = "(const double q_, const double m_, const double kappaPar, const double kappaPerp, const double kperpSq, const double *wL1, const double *dxL1, const double *wR1, const double *dxR1, const double uMaxIn, const double *jacL, const double *rBmagL, const double *jacDbmagL, const double *sMomL1, const double *sMomR1, const double *phiL1, const double *phiR1, double *primMomL1, const double *primMomR1, double *outL, double *outR)"

   local funcNm = {}
   if dim == 1 then
      funcNm[1] = string.format("gyrofluid_surf_%dx_p%d_%s_x", dim, polyOrder, basisNmMap[basisNm])
   else
      assert(false, "Gyrofluid equation not implemented for this dimensionality!")
   end

   local CDefStr = ""
   for d = 1, #funcNm do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, #funcNm do
      local tmp = ffi.C[funcNm[d]]
      kernels[d] = tmp
   end
   return kernels
end

return _M
