-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing cell-ave Maxwellian.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- Kernel function to compute cell-average Maxwellian for Vlasov grid.
function _M.CellAvMax(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("MaxwellianCellAv%s%dx%dv_P%d", basisNmMap[basisNm], CDIM, VDIM, polyOrder)
   local funcSign = "(const double *w, const double *m0, const double *u, const double *vtSq, double *fMax)"
   
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Kernel function to compute cell-average Maxwellian for Gyrokinetic grid.
function _M.GkCellAvMax(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("GkMaxwellianCellAv%s%dx%dv_P%d", basisNmMap[basisNm], CDIM, VDIM, polyOrder)
   local funcSign = "(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax)"
     
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
