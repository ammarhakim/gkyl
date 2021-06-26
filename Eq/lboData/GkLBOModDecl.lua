-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Gyrokinetic Lenard-Bernstein operator C++ kernel functions based on CDIM, VDIM,
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- "do nothing" function.
local function nullFunc(...) end

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" , ["tensor"] = "Tensor"}

local vvars = {"Vpar","Mu"}

local _M = {}

-- Select function to compute volume terms.
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = string.format("GkLBOVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end
function _M.selectConstNuVol(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = string.format("GkLBOconstNuVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface terms (output is a table of functions).
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder, isNonuniform, applyPos)
   local funcType = "double"
   local funcNm = {}
   for d = 1, VDIM do
      funcNm[d] = string.format("GkLBOSurf%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], vvars[d], polyOrder)
   end
   local funcSign = "(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr)"

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
function _M.selectConstNuSurf(basisNm, CDIM, VDIM, polyOrder, isNonuniform, applyPos)
   local gridStr = isNonuniform and "Nonuniform" or ""

   local posString = ""
   if applyPos then posString = "Positivity" end

   local funcType = "double"
   local funcNm = {}
   for d = 1, VDIM do
      funcNm[d] = string.format("GkLBOconstNuSurf%s%s%dx%dv%s_%s_P%d", gridStr, posString, CDIM, VDIM, basisNmMap[basisNm], vvars[d], polyOrder)
   end
   local funcSign = "(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr)"

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
function _M.selectBoundarySurf(basisNm, CDIM, VDIM, polyOrder, applyPos)
   local funcType = "double"
   local funcNm = {}
   for d = 1, VDIM do
      funcNm[d] = string.format("GkLBOBoundarySurf%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], vvars[d], polyOrder)
   end
   local funcSign = "(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr)"

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
function _M.selectConstNuBoundarySurf(basisNm, CDIM, VDIM, polyOrder, applyPos)
   local funcType = "double"
   local funcNm = {}
   for d = 1, VDIM do
      funcNm[d] = string.format("GkLBOconstNuBoundarySurf%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], vvars[d], polyOrder)
   end
   local funcSign = "(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr)"

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
