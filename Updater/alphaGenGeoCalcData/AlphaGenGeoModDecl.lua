-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Vlasov C++ kernel functions based on CDIM, VDIM,
-- basis functions and polyOrder.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- "do nothing" function.
local function nullFunc(...) end

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local cvars = {"X","Y","Z"}
local vvars = {"VX","VY","VZ"}

local _M = {}

-- Select function to compute gen geo alpha for streaming terms.
function _M.AlphaGenGeo(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("AlphaGenGeo%s%dx%dvP%d", basisNmMap[basisNm], CDIM, VDIM, polyOrder)
   local funcSign = "(const double *w, const double *dxv, const double *tvComp, double const *gxx, double const *gxy, double const *gxz, double const *gyy, double const *gyz, double const *gzz, double *alphaGeo)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
