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

-- Select function to compute volume constDiffusion terms.
function _M.selectVol(basisNm, CDIM, polyOrder)
   local funcNm = string.format("ConstDiffusionVol%dx%sP%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select functions to compute surface constDiffusion terms (output is a table of functions).
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

-- Select functions to compute boundary surface constDiffusion terms (output is a table of functions).
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

-- Select kernels that implement BCs specific to constDiffusion term (output is a table of functions).
function _M.selectBCs(basisNm, CDIM, polyOrder, bcType)
   if CDIM == 1 then
      local funcNmXl = string.format("ConstDiffusionBC%dx%s%s_Xlower_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      local funcNmXu = string.format("ConstDiffusionBC%dx%s%s_Xupper_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      return { {ffi.C[funcNmXl],ffi.C[funcNmXu]} }                                                                     
   elseif CDIM == 2 then                                                                                               
      local funcNmXl = string.format("ConstDiffusionBC%dx%s%s_Xlower_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      local funcNmXu = string.format("ConstDiffusionBC%dx%s%s_Xupper_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      local funcNmYl = string.format("ConstDiffusionBC%dx%s%s_Ylower_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      local funcNmYu = string.format("ConstDiffusionBC%dx%s%s_Yupper_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      return { {ffi.C[funcNmXl],ffi.C[funcNmXu]}, {ffi.C[funcNmYl],ffi.C[funcNmYu]} }                                  
   elseif CDIM == 3 then                                                                                               
      local funcNmXl = string.format("ConstDiffusionBC%dx%s%s_Xlower_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      local funcNmXu = string.format("ConstDiffusionBC%dx%s%s_Xupper_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      local funcNmYl = string.format("ConstDiffusionBC%dx%s%s_Ylower_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      local funcNmYu = string.format("ConstDiffusionBC%dx%s%s_Yupper_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      local funcNmZl = string.format("ConstDiffusionBC%dx%s%s_Zlower_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      local funcNmZu = string.format("ConstDiffusionBC%dx%s%s_Zupper_P%d", CDIM, basisNmMap[basisNm], bcType, polyOrder)
      return { {ffi.C[funcNmXl],ffi.C[funcNmXu]}, {ffi.C[funcNmYl],ffi.C[funcNmYu]}, {ffi.C[funcNmZl],ffi.C[funcNmZu]} }
   end
end

return _M
