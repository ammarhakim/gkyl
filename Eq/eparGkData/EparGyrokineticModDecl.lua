-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into EparGyrokinetic C++ kernel functions based on CDIM
-- basis functions and polyOrder.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["tensor"] = "Tensor"}

local _M = {}

-- Select function to compute volume terms.
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder)
   local funcSign = "(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *Epar, const double *f, double *out)"

   local funcType = "double"
   local funcNm = string.format("EparGyrokineticVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface terms (output is a table of functions).
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder)
   local funcType  = "double"
   local funcSign = "(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *Epar, const double *fL, const double *fR, double *outL, double *outR)"
   local funcNm = {}
   if CDIM == 1 and VDIM <= 2 then
      funcNm[1] = string.format("EparGyrokineticSurf%dx%dv%s_x_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      funcNm[2] = string.format("EparGyrokineticSurf%dx%dv%s_vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   elseif CDIM == 2 and VDIM == 2 then
      assert(false, "EparGyrokinetic equation not implemented for this dimensionality!")
   elseif CDIM == 3 and VDIM == 2 then
      assert(false, "EparGyrokinetic equation not implemented for this dimensionality!")
   elseif CDIM == 2 and VDIM == 0 then 
      assert(false, "EparGyrokinetic equation not implemented for this dimensionality!")
   else
      assert(false, "EparGyrokinetic equation not implemented for this dimensionality!")
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
