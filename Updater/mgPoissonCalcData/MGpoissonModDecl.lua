-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into multigrid Poisson solver C++ kernel functions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Updater.mgPoissonCalcData._MGpoissonCdef"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select restriction operator kernel.
function _M.selectRestriction(basisNm, dim, polyOrder)
   local funcNm = string.format("MGpoissonRestrict%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select prolongation operator kernel.
function _M.selectProlongation(basisNm, dim, polyOrder)
   local funcNm = string.format("MGpoissonProlong%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select relaxation kernels.
function _M.selectRelaxation(basisNm, dim, polyOrder, kindOfRelax)
   local funcNm = string.format("MGpoisson%s%dx%s_P%d", kindOfRelax, dim, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- Select residue kernels.
function _M.selectResidueCalc(basisNm, dim, polyOrder)
   local funcNm = string.format("MGpoissonResidue%dx%s_P%d", dim, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
