-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into PassiveAdvection C++ kernel functions based on CDIM
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Eq.passiveAdvectionData._PassiveAdvectionCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select function to compute volume  terms
function _M.selectVol(basisNm, CDIM, polyOrder)
   local funcNm = string.format("PassiveAdvectionVol%dx%sP%d", CDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- select functions to compute surface terms (output is a table of functions)
function _M.selectSurf(basisNm, CDIM, polyOrder, positivity)
   assert(CDIM<=5, "PassiveAdvection equation not implemented for this dimensionality!")

   local posString = ""
   if positivity then posString = "Positivity" end

   local funcTbl = {}

   for i = 1, CDIM do
      local funcNm = string.format("PassiveAdvectionSurf%s%dx%s_X%d_P%d", posString, CDIM, basisNmMap[basisNm], i, polyOrder)
      funcTbl[i] = ffi.C[funcNm]
   end

   return funcTbl
end

return _M
