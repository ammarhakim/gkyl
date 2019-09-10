-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into ConstDiffusion C++ kernel functions based on CDIM,
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Eq.constDiffusionData._ConstDiffusionCdef"

-- "do nothing" function
local function nullFunc(...) end

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Select function to compute volume streaming terms.
function _M.selectVol(basisNm, CDIM, polyOrder)
   local funcNm = string.format("ConstDiffusionVol%dx%sP%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select functions to compute surface streaming terms (output is a table of functions).
function _M.selectSurf(basisNm, CDIM, polyOrder, applyPos)
   local posString = ""
   if applyPos then posString = "Positivity" end
   if CDIM == 1 then
      local funcNmX = string.format("ConstDiffusionSurf%s%dx%s_X_P%d", posString, CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX] }
   elseif CDIM == 2 then
      local funcNmX = string.format("ConstDiffusionSurf%s%dx%s_X_P%d", posString, CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("ConstDiffusionSurf%s%dx%s_Y_P%d", posString, CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY] }
   elseif CDIM == 3 then
      local funcNmX = string.format("ConstDiffusionSurf%s%dx%s_X_P%d", posString, CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("ConstDiffusionSurf%s%dx%s_Y_P%d", posString, CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("ConstDiffusionSurf%s%dx%s_Z_P%d", posString, CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end

-- Select functions to compute boundary surface streaming terms (output is a table of functions).
function _M.selectBoundarySurf(basisNm, CDIM, polyOrder, applyPos)
   if CDIM == 1 then
      local funcNmX = string.format("ConstDiffusionBoundarySurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX] }
   elseif CDIM == 2 then
      local funcNmX = string.format("ConstDiffusionBoundarySurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("ConstDiffusionBoundarySurf%dx%s_Y_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY] }
   elseif CDIM == 3 then
      local funcNmX = string.format("ConstDiffusionBoundarySurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("ConstDiffusionBoundarySurf%dx%s_Y_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("ConstDiffusionBoundarySurf%dx%s_Z_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end

return _M
