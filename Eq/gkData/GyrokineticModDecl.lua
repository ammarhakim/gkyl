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
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder, isElectromagnetic)
   if isElectromagnetic then emString = "Em" else emString = "" end
   local funcNm = string.format("%sGyrokineticVol%dx%dv%sP%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- select functions to compute surface terms (output is a table of functions)
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder, isElectromagnetic)
   if isElectromagnetic then emString = "Em" else emString = "" end
   if CDIM == 1 and VDIM == 1 then
      local funcNmX = string.format("%sGyrokineticSurf%dx%dv%s_X_P%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVpar = string.format("%sGyrokineticSurf%dx%dv%s_Vpar_P%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmVpar] }
   elseif CDIM == 2 and VDIM == 2 then
      local funcNmX = string.format("%sGyrokineticSurf%dx%dv%s_X_P%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("%sGyrokineticSurf%dx%dv%s_Y_P%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVpar = string.format("%sGyrokineticSurf%dx%dv%s_Vpar_P%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmVpar] }
   elseif CDIM == 3 and VDIM == 2 then
      local funcNmX = string.format("%sGyrokineticSurf%dx%dv%s_X_P%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("%sGyrokineticSurf%dx%dv%s_Y_P%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("%sGyrokineticSurf%dx%dv%s_Z_P%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVpar = string.format("%sGyrokineticSurf%dx%dv%s_Vpar_P%d", emString, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ], ffi.C[funcNmVpar] }
   else
      assert(false, "Gyrokinetic equation not implemented for this dimensionality!")
   end
end

function _M.select_dAdtVol(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("dAdtVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

function _M.select_dAdtSurf(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("dAdtSurf%dx%dv%s_Vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

function _M.selectSheathDeltaPhi(basisNm, CDIM, polyOrder)
   local funcNm = string.format("calcSheathDeltaPhi%dx%s_P%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end
function _M.selectSheathPartialReflection(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("calcSheathPartialReflection%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
