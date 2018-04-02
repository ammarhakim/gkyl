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
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("GyrokineticVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- select functions to compute surface terms (output is a table of functions)
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder)
   if CDIM == 1 and VDIM == 1 then
      local funcNmX = string.format("GyrokineticSurf%dx%dv%s_X_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVpar = string.format("GyrokineticSurf%dx%dv%s_Vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmVpar] }
   elseif CDIM == 2 and VDIM == 2 then
      local funcNmX = string.format("GyrokineticSurf%dx%dv%s_X_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("GyrokineticSurf%dx%dv%s_Y_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVpar = string.format("GyrokineticSurf%dx%dv%s_Vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmMu = string.format("GyrokineticSurf%dx%dv%s_Mu_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmVpar], ffi.C[funcNmMu] }
   elseif CDIM == 3 and VDIM == 2 then
      local funcNmX = string.format("GyrokineticSurf%dx%dv%s_X_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("GyrokineticSurf%dx%dv%s_Y_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("GyrokineticSurf%dx%dv%s_Z_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVpar = string.format("GyrokineticSurf%dx%dv%s_Vpar_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      local funcNmMu = string.format("GyrokineticSurf%dx%dv%s_Mu_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmVpar], ffi.C[funcNmMu] }
   else
      assert(false, "Gyrokinetic equation not implemented for this dimensionality!")
   end
end

return _M
