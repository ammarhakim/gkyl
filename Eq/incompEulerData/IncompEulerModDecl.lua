-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Incompressible Euler C++ kernel functions based on CDIM
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Eq.incompEulerData._IncompEulerCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Select function to compute volume terms.
function _M.selectVol(basisNm, CDIM, polyOrder)
   local funcNm = string.format("IncompEulerVol%dx%sP%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select functions to compute surface terms (output is a table of functions).
function _M.selectSurf(basisNm, CDIM, polyOrder, positivity)
   local posString = ""
   if positivity then posString = "Positivity" end
   if CDIM == 2 then 
      local funcNmX = string.format("IncompEulerSurf%s%dx%s_X_P%d", posString, CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("IncompEulerSurf%s%dx%s_Y_P%d", posString, CDIM, basisNmMap[basisNm], polyOrder)
      return { ffi.C[funcNmX], ffi.C[funcNmY]}
   else
      assert(false, "IncompEuler equation not implemented for this dimensionality!")
   end
end

return _M
