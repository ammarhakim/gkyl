-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into moment calculation C++ kernel functions based on
-- CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
if GKYL_HAVE_CUDA then
   _ = require "Updater.momentCalcData._DistFuncMomentCalcDeviceCdef"
end

local vvars = {"VX","VY","VZ"}

-- Map of basis function name -> function encoding.
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- Select function to compute specified moment.
function _M.selectMomCalc(mom, basisNm, CDIM, VDIM, polyOrder, calcOnDevice)
   local funcNm
   if not calcOnDevice then
      local funcType = "void"
      local funcSign = "(const double *w, const double *dxv, const double *f, double *out)"
      if mom=="FiveMoments" then
         funcSign = "(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2)"
      end
      funcNm = string.format("MomentCalc%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], mom, polyOrder)
      ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   else
      funcNm = string.format("cuda_MomentCalc%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], mom, polyOrder)
   end
   return ffi.C[funcNm]
end

-- Select function to compute specified gyrokinetic moment.
function _M.selectGkMomCalc(mom, basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("GkMomentCalc%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], string.sub(mom,3), polyOrder)
   local funcSign = "(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out)"
   if mom=="GkThreeMoments" then
      funcSign = "(const double *w, const double *dxv, const double m_, const double *Bmag, const double *f, double *out1, double *out2, double *out3)"
   end

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- Select function to compute integrated moments.
function _M.selectIntMomCalc(basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("IntMomentCalc%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   local funcSign = "(const double *w, const double *dxv, const double *f, double *out)"

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- With piecewise linear selfPrimMoments needs the star moments.
function _M.selectStarM0Calc(dir, kinSpecies, basisNm, CDIM, VDIM, applyPos)
   local funcType = "void"
   local posString = ""
   if applyPos then posString = "Positivity" end
   local funcNm = string.format("%sM0Star%s%dx%dv%s_%s", kinSpecies, posString, CDIM, VDIM, basisNmMap[basisNm], vvars[dir])
   local funcSign
   if kinSpecies=="Vm" then
      funcSign = "(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out)"
   elseif kinSpecies=="Gk" then
      funcSign = "(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out)"
   end

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end
function _M.selectStarM1iM2Calc(kinSpecies, basisNm, CDIM, VDIM)
   local funcType = "void"
   local funcNm = string.format("%sM1iM2Star%dx%dv%s", kinSpecies, CDIM, VDIM, basisNmMap[basisNm])
   local funcSign
   if kinSpecies=="Vm" then
      funcSign = "(const double *w, const double *dxv, const double *f, double *outM1i, double *outM2)"
   elseif kinSpecies=="Gk" then
      funcSign = "(const double *w, const double *dxv, const double intFac, const double m_, const double *Bmag, const double *f, double *outM1i, double *outM2)"
   end

   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

-- S1lfPrimMoments also need kernels to compute certain integrals along velocity boundaries.
function _M.selectBoundaryFintegral(dir, kinSpecies, basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("%sBoundaryIntegral%dx%dv%s_F_%s_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], vvars[dir], polyOrder)
   local funcSign
   if kinSpecies=="Vm" then
      funcSign = "(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out)"
   elseif kinSpecies=="Gk" then
      funcSign = "(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out)"
   end
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end
function _M.selectBoundaryVFintegral(dir, kinSpecies, basisNm, CDIM, VDIM, polyOrder)
   local funcType = "void"
   local funcNm = string.format("%sBoundaryIntegral%dx%dv%s_vF_%s_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], vvars[dir], polyOrder)
   local funcSign
   if kinSpecies=="Vm" then
      funcSign = "(const bool atLower, const double vBoundary, const double *dxv, const double *fIn, double *out)"
   elseif kinSpecies=="Gk" then
      funcSign = "(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out)"
   end
   ffi.cdef(funcType .. " " .. funcNm .. funcSign .. ";\n")
   return ffi.C[funcNm]
end

return _M
