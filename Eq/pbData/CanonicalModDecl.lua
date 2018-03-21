-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Canonical C++ kernel functions based on CDIM
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Eq.pbData._CanonicalCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select function to compute volume  terms
function _M.selectVol(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("CanonicalVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- select functions to compute surface streaming terms (output is a table of functions)
function _M.selectSurf(basisNm, CDIM, VDIM, polyOrder)
   if CDIM == 1 and VDIM == 1 then
      local funcNmX = string.format("CanonicalSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVX = string.format("CanonicalSurf%dx%s_VX_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmVX] }
   elseif CDIM == 2 and VDIM == 2 then
      local funcNmX = string.format("CanonicalSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("CanonicalSurf%dx%s_Y_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVX = string.format("CanonicalSurf%dx%s_VX_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmVY = string.format("CanonicalSurf%dx%s_VY_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmVX], ffi.C[funcNmVY] }
   else
      assert(false, "Must have CDIM=VDIM<=2 for Hamiltonian equation with canonical poisson bracket!")
   end
end

return _M
