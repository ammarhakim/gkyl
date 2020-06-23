-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Maxwell C++ kernel functions based on CDIM,
-- basis functions and polyOrder.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local cvars = {"X","Y","Z"}

local _M = {}

-- Select function to compute volume streaming terms.
function _M.selectVol(basisNm, CDIM, polyOrder)
   local funcType = "double"
   local funcSign = "(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out)"

   local funcNm = string.format("MaxwellVol%dx%sP%d", CDIM, basisNmMap[basisNm], polyOrder)

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface streaming terms using upwind fluxes.
function _M.selectSurf(basisNm, CDIM, polyOrder)
   local funcType = "double"
   local funcSign = "(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr)"

   local funcNm = {}
   for d = 1, CDIM do
      funcNm[d] = string.format("MaxwellSurf%dx%s_%s_P%d", CDIM, basisNmMap[basisNm], cvars[d], polyOrder)
   end

   local CDefStr = ""
   for d = 1, CDIM do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, CDIM do
      local tmp = ffi.C[funcNm[d]]
      kernels[d] = tmp 
   end
   return kernels
end

-- Select functions to compute surface streaming terms using central fluxes.
function _M.selectCentralSurf(basisNm, CDIM, polyOrder)
   local funcType = "double"
   local funcSign = "(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr)"

   local funcNm = {}
   for d = 1, CDIM do
      funcNm[d] = string.format("MaxwellCentralSurf%dx%s_%s_P%d", CDIM, basisNmMap[basisNm], cvars[d], polyOrder)
   end

   local CDefStr = ""
   for d = 1, CDIM do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, CDIM do
      local tmp = ffi.C[funcNm[d]]
      kernels[d] = tmp 
   end
   return kernels
end

return _M
