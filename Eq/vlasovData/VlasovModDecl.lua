-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into Vlasov C++ kernel functions based on CDIM, VDIM,
-- basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"

-- "do nothing" function
local function nullFunc(...) end

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local cvars = {"X","Y","Z"}
local vvars = {"VX","VY","VZ"}

local _M = {}

-- Select function to compute only volume streaming terms.
function _M.selectVolStream(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = string.format("VlasovVolStream%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *w, const double *dxv, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface streaming terms (output is a table of functions).
function _M.selectSurfStream(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = {}
   for d = 1, CDIM do
      funcNm[d] = string.format("VlasovSurfStream%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], cvars[d], polyOrder)
   end
   local funcSign = "(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr)"

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

-- Select function to compute total volume (streaming + EM field acceleration) term.
function _M.selectVolElcMag(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = string.format("VlasovVol%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *w, const double *dxv, const double *E, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select functions to compute surface EM field acceleration  terms (output is a table of functions).
function _M.selectSurfElcMag(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = {}
   for d = 1, VDIM do
      funcNm[d] = string.format("VlasovSurfElcMag%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], vvars[d], polyOrder)
   end
   local funcSign = "(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *E, const double *fl, const double *fr, double *outl, double *outr)"

   local CDefStr = ""
   for d = 1, VDIM do CDefStr = CDefStr .. (funcType .. " " .. funcNm[d] .. funcSign .. ";\n") end
   ffi.cdef(CDefStr)

   local kernels = {}
   for d = 1, VDIM do
      local tmp = ffi.C[funcNm[d]]
      kernels[d] = tmp
   end
   return kernels
end

-- Select function to compute only volume EM field acceleration terms.
function _M.selectVolForce(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "double"
   local funcNm = string.format("VlasovVolForce%dx%dv%sP%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *w, const double *dxv, const double *E, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end


return _M
