-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Gyrokinetic C++ kernel functions based on CDIM
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Eq.gkData._GyrokineticCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select function to compute volume  terms
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder, isElectromagnetic, positivity, Bvars)
   local emString = ""
   if isElectromagnetic then emString = "Em" end
   local posString = ""
   if positivity then posString = "Positivity" end
   local bvarString = "_Bvars"
   for k, v in ipairs(Bvars) do
      bvarString = bvarString .. "_" .. v
   end
   local funcNm = string.format("%sGyrokineticVol%s%dx%dv%s_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm..bvarString]
end

-- select functions to compute surface terms (output is a table of functions)
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder, isElectromagnetic, positivity, Bvars)
   local emString = ""
   if isElectromagnetic then emString = "Em" end
   local posString = ""
   if positivity then posString = "Positivity" end
   local bvarString = "_Bvars"
   for k, v in ipairs(Bvars) do
      bvarString = bvarString .. "_" .. v
   end
   if CDIM == 1 and VDIM <= 2 then
      local funcNmX = string.format("%sGyrokineticSurf%s%dx%dv%s_X_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVpar = string.format("%sGyrokineticSurf%s%dx%dv%s_Vpar_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX..bvarString], ffi.C[funcNmVpar..bvarString] }
   elseif CDIM == 2 and VDIM == 2 then
      local funcNmX = string.format("%sGyrokineticSurf%s%dx%dv%s_X_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("%sGyrokineticSurf%s%dx%dv%s_Y_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVpar = string.format("%sGyrokineticSurf%s%dx%dv%s_Vpar_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX..bvarString], ffi.C[funcNmY..bvarString], ffi.C[funcNmVpar..bvarString] }
   elseif CDIM == 3 and VDIM == 2 then
      local funcNmX = string.format("%sGyrokineticSurf%s%dx%dv%s_X_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("%sGyrokineticSurf%s%dx%dv%s_Y_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("%sGyrokineticSurf%s%dx%dv%s_Z_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVpar = string.format("%sGyrokineticSurf%s%dx%dv%s_Vpar_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX..bvarString], ffi.C[funcNmY..bvarString], ffi.C[funcNmZ..bvarString], ffi.C[funcNmVpar..bvarString] }
   elseif CDIM == 2 and VDIM == 0 then 
      local funcNmX = string.format("%sGyrokineticSurf%s%dx%dv%s_X_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("%sGyrokineticSurf%s%dx%dv%s_Y_P%d", emString, posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX..bvarString], ffi.C[funcNmY..bvarString]}
   else
      assert(false, "Gyrokinetic equation not implemented for this dimensionality!")
   end
end

function _M.selectStep2Vol(basisNm, CDIM, VDIM, polyOrder, positivity)
   local posString = ""
   if positivity then posString = "Positivity" end
   local funcNm = string.format("EmGyrokineticStep2Vol%s%dx%dv%s_P%d", posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

function _M.selectStep2Surf(basisNm, CDIM, VDIM, polyOrder, positivity, Bvars)
   local emString = "Em"
   local posString = ""
   if positivity then posString = "Positivity" end
   local bvarString = "_Bvars"
   for k, v in ipairs(Bvars) do
      bvarString = bvarString .. "_" .. v
   end
   local funcNm
   if polyOrder > 1 then
      funcNm = string.format("EmGyrokineticSurf%s%dx%dv%s_Vpar_P%d", posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   else
      funcNm = string.format("EmGyrokineticSurf%s%dx%dv%sStep2_Vpar_P%d", posString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   end
   return ffi.C[funcNm..bvarString]
end

function _M.selectSheathDeltaPhi(basisNm, CDIM, polyOrder)
   local funcNm = string.format("calcSheathDeltaPhi%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end
function _M.selectSheathPartialReflection(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("calcSheathPartialReflectionScaled%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end
function _M.selectSheathReflection(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("calcSheathReflection%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
