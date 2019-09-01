-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Gyrokinetic Lenard-Bernstein operator C++ kernel functions based on CDIM, VDIM,
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Eq.lboData._GkLBOCdef"

-- "do nothing" function.
local function nullFunc(...) end

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Select function to compute volume terms.
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("GkLBOVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end
function _M.selectConstNuVol(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("GkLBOconstNuVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select functions to compute surface terms (output is a table of functions).
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder)
   if VDIM == 1 then
      local funcNmX = string.format("GkLBOSurf%dx%dv%s_Vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], nullFunc, nullFunc }
   elseif VDIM == 2 then
      local funcNmX = string.format("GkLBOSurf%dx%dv%s_Vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("GkLBOSurf%dx%dv%s_Mu_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], nullFunc }
   end
end
function _M.selectConstNuSurf(basisNm, CDIM, VDIM, polyOrder, applyPos)
   local posString = ""
   if applyPos then posString = "Positivity" end
   if VDIM == 1 then
      local funcNmX = string.format("GkLBOconstNuSurf%s%dx%dv%s_Vpar_P%d", posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], nullFunc, nullFunc }
   elseif VDIM == 2 then
      local funcNmX = string.format("GkLBOconstNuSurf%s%dx%dv%s_Vpar_P%d", posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("GkLBOconstNuSurf%s%dx%dv%s_Mu_P%d", posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], nullFunc }
   end
end

-- Select functions to compute boundary surface terms (output is a table of functions).
function _M.selectBoundarySurf(basisNm, CDIM, VDIM, polyOrder, applyPos)
   if VDIM == 1 then
      local funcNmX = string.format("GkLBOBoundarySurf%dx%dv%s_Vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], nullFunc, nullFunc }
   elseif VDIM == 2 then
      local funcNmX = string.format("GkLBOBoundarySurf%dx%dv%s_Vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("GkLBOBoundarySurf%dx%dv%s_Mu_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], nullFunc }
   end
end
function _M.selectConstNuBoundarySurf(basisNm, CDIM, VDIM, polyOrder, applyPos)
   if VDIM == 1 then
      local funcNmX = string.format("GkLBOconstNuBoundarySurf%dx%dv%s_Vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], nullFunc, nullFunc }
   elseif VDIM == 2 then
      local funcNmX = string.format("GkLBOconstNuBoundarySurf%dx%dv%s_Vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("GkLBOconstNuBoundarySurf%dx%dv%s_Mu_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], nullFunc }
   end
end


return _M
