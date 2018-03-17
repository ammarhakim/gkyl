-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Canonical C++ kernel functions based on CDIM
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Eq.poissonBracketData.CanonicalCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select function to compute volume streaming terms
function _M.selectVol(basisNm, CDIM, polyOrder)
   local funcNm = string.format("CanonicalVol%dx%sP%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- select functions to compute surface streaming terms (output is a table of functions)
function _M.selectSurf(basisNm, CDIM, polyOrder)
   if CDIM == 1 then
      local funcNmX = string.format("CanonicalSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX] }
   elseif CDIM == 2 then
      local funcNmX = string.format("CanonicalSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("CanonicalSurf%dx%s_Y_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY] }
   elseif CDIM == 3 then
      local funcNmX = string.format("CanonicalSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("CanonicalSurf%dx%s_Y_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("CanonicalSurf%dx%s_Z_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end

return _M
