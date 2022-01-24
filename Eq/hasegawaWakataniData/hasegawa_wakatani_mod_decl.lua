-- Gkyl ------------------------------------------------------------------------
--
-- Select the Hasegawa-Wakatani C++ kernel functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["tensor"] = "Tensor" }

local _M = {}

-- Select function to compute volume terms.
function _M.selectVol(basisNm, polyOrder)
   local funcSign = "(const double C_, const double kappa_, const double *xc, const double *dx, const double *phi, const double *fIn, double *out)"

   local funcType = "double"
   local funcNm = string.format("hasegawa_wakatani_vol_2x_p%d_%s", polyOrder, basisNmMap[basisNm])

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface terms (output is a table of functions).
function _M.selectSurf(basisNm, polyOrder)
   local funcType  = "double"
   local funcSign = "(const double C_, const double kappa_, const double cflL, const double cflR, const double *xc, const double *dx, const double amax_in, const double *phi, const double *fInL, const double *fInR, double *outL, double *outR)"

   local funcNm = {}
   funcNm[1] = string.format("hasegawa_wakatani_surf_2x_p%d_%s_x", polyOrder, basisNmMap[basisNm])
   funcNm[2] = string.format("hasegawa_wakatani_surf_2x_p%d_%s_y", polyOrder, basisNmMap[basisNm])

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
