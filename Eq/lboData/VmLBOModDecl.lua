-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Vlasov Lenard-Bernstein operator C++ kernel functions based on CDIM, VDIM,
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Eq.lboData._VmLBOCdef"

-- "do nothing" function
local function nullFunc(...) end

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- select function to compute volume terms
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("VmLBOVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end
function _M.selectConstNuVol(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("VmLBOconstNuVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- select functions to compute surface terms (output is a table of functions)
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder)
   if VDIM == 1 then
      local funcNmX = string.format("VmLBOSurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], nullFunc, nullFunc }
   elseif VDIM == 2 then
      local funcNmX = string.format("VmLBOSurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("VmLBOSurf%dx%dv%s_VY_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], nullFunc }
   elseif VDIM == 3 then
      local funcNmX = string.format("VmLBOSurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("VmLBOSurf%dx%dv%s_VY_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("VmLBOSurf%dx%dv%s_VZ_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end
function _M.selectConstNuSurf(basisNm, CDIM, VDIM, polyOrder)
   if VDIM == 1 then
      local funcNmX = string.format("VmLBOconstNuSurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], nullFunc, nullFunc }
   elseif VDIM == 2 then
      local funcNmX = string.format("VmLBOconstNuSurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("VmLBOconstNuSurf%dx%dv%s_VY_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], nullFunc }
   elseif VDIM == 3 then
      local funcNmX = string.format("VmLBOconstNuSurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("VmLBOconstNuSurf%dx%dv%s_VY_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("VmLBOconstNuSurf%dx%dv%s_VZ_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end

-- select functions to compute boundary surface terms (output is a table of functions)
function _M.selectBoundarySurf(basisNm, CDIM, VDIM, polyOrder)
   if VDIM == 1 then
      local funcNmX = string.format("VmLBOBoundarySurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], nullFunc, nullFunc }
   elseif VDIM == 2 then
      local funcNmX = string.format("VmLBOBoundarySurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("VmLBOBoundarySurf%dx%dv%s_VY_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], nullFunc }
   elseif VDIM == 3 then
      local funcNmX = string.format("VmLBOBoundarySurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("VmLBOBoundarySurf%dx%dv%s_VY_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("VmLBOBoundarySurf%dx%dv%s_VZ_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end
function _M.selectConstNuBoundarySurf(basisNm, CDIM, VDIM, polyOrder)
   if VDIM == 1 then
      local funcNmX = string.format("VmLBOconstNuBoundarySurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], nullFunc, nullFunc }
   elseif VDIM == 2 then
      local funcNmX = string.format("VmLBOconstNuBoundarySurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("VmLBOconstNuBoundarySurf%dx%dv%s_VY_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], nullFunc }
   elseif VDIM == 3 then
      local funcNmX = string.format("VmLBOconstNuBoundarySurf%dx%dv%s_VX_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("VmLBOconstNuBoundarySurf%dx%dv%s_VY_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("VmLBOconstNuBoundarySurf%dx%dv%s_VZ_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end


return _M
