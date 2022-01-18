-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Vlasov Lenard-Bernstein operator C++ kernel functions based on CDIM, VDIM,
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local vvars = {"VX","VY","VZ"}

local _M = {}

-- Select function to compute volume terms.
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = string.format("VmLBOVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end
function _M.selectConstNuVol(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = string.format("VmLBOconstNuVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface terms (output is a table of functions).
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder, fluxType, isNonuniform)
   local fluxStr = fluxType=="penalty" and "" or "Upwind"
   local gridStr = isNonuniform and "Nonuniform" or ""

   local funcType = "double"
   local funcNm   = {}
   for d = 1, VDIM do
      funcNm[d] = string.format("VmLBO%sSurf%s%dx%dv%s_%s_P%d", fluxStr, gridStr, CDIM, VDIM, basisNmMap[basisNm], vvars[d], polyOrder)
   end
   local funcSign = "(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr)"

   local CDefStr = ""
   for d = 1, VDIM do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, VDIM do
      local tmp = ffi.C[funcNm[d]]
      kernels[d] = tmp
   end
   return kernels
end
function _M.selectConstNuSurf(basisNm, CDIM, VDIM, polyOrder, fluxType, isNonuniform)
   local fluxStr = fluxType=="penalty" and "" or "Upwind"
   local gridStr = isNonuniform and "Nonuniform" or ""

   local funcType = "double"
   local funcNm = {}
   for d = 1, VDIM do
      funcNm[d] = string.format("VmLBOconstNu%sSurf%s%dx%dv%s_%s_P%d", fluxStr, gridStr, CDIM, VDIM, basisNmMap[basisNm], vvars[d], polyOrder)
   end
   local funcSign = "(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr)"

   local CDefStr = ""
   for d = 1, VDIM do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, VDIM do
      local tmp = ffi.C[funcNm[d]]
      kernels[d] = tmp
   end
   return kernels
end

-- Select functions to compute boundary surface terms (output is a table of functions).
function _M.selectBoundarySurf(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = {}
   for d = 1, VDIM do
      funcNm[d] = string.format("VmLBOBoundarySurf%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], vvars[d], polyOrder)
   end
   local funcSign = "(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr)"

   local CDefStr = ""
   for d = 1, VDIM do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, VDIM do
      local tmp = ffi.C[funcNm[d]]
      kernels[d] = tmp
   end
   return kernels
end
function _M.selectConstNuBoundarySurf(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = {}
   for d = 1, VDIM do
      funcNm[d] = string.format("VmLBOconstNuBoundarySurf%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], vvars[d], polyOrder)
   end
   local funcSign = "(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const int edge, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr)"

   local CDefStr = ""
   for d = 1, VDIM do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, VDIM do
      local tmp = ffi.C[funcNm[d]]
      kernels[d] = tmp
   end
   return kernels
end


return _M
