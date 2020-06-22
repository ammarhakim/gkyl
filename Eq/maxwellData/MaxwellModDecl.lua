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
   if CDIM == 1 then
      local funcNmX = string.format("MaxwellSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      ffi.cdef(funcType .. " " .. funcNmX .. funcSign .. ";\n")
      return { ffi.C[funcNmX] }
   elseif CDIM == 2 then
      local funcNmX = string.format("MaxwellSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("MaxwellSurf%dx%s_Y_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      ffi.cdef(funcType .. " " .. funcNmX .. funcSign .. ";\n" ..
               funcType .. " " .. funcNmY .. funcSign .. ";\n")
      return { ffi.C[funcNmX], ffi.C[funcNmY] }
   elseif CDIM == 3 then
      local funcNmX = string.format("MaxwellSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("MaxwellSurf%dx%s_Y_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("MaxwellSurf%dx%s_Z_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      ffi.cdef(funcType .. " " .. funcNmX .. funcSign .. ";\n" ..
               funcType .. " " .. funcNmY .. funcSign .. ";\n" ..
               funcType .. " " .. funcNmZ .. funcSign .. ";\n")
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end

-- Select functions to compute surface streaming terms using central fluxes.
function _M.selectCentralSurf(basisNm, CDIM, polyOrder)
   local funcType = "double"
   local funcSign = "(const MaxwellEq_t *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double *outl, double *outr)"
   if CDIM == 1 then
      local funcNmX = string.format("MaxwellCentralSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      ffi.cdef(funcType .. " " .. funcNmX .. funcSign .. ";\n")
      return { ffi.C[funcNmX] }
   elseif CDIM == 2 then
      local funcNmX = string.format("MaxwellCentralSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("MaxwellCentralSurf%dx%s_Y_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      ffi.cdef(funcType .. " " .. funcNmX .. funcSign .. ";\n" ..
               funcType .. " " .. funcNmY .. funcSign .. ";\n")
      return { ffi.C[funcNmX], ffi.C[funcNmY] }
   elseif CDIM == 3 then
      local funcNmX = string.format("MaxwellCentralSurf%dx%s_X_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmY = string.format("MaxwellCentralSurf%dx%s_Y_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      local funcNmZ = string.format("MaxwellCentralSurf%dx%s_Z_P%d", CDIM, basisNmMap[basisNm], polyOrder)
      ffi.cdef(funcType .. " " .. funcNmX .. funcSign .. ";\n" ..
               funcType .. " " .. funcNmY .. funcSign .. ";\n" ..
               funcType .. " " .. funcNmZ .. funcSign .. ";\n")
      return { ffi.C[funcNmX], ffi.C[funcNmY], ffi.C[funcNmZ] }
   end
end

return _M
