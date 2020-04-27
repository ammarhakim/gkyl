-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into moment calculation C++ kernel functions based on
-- CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _   = require "Updater.momentCalcData._DistFuncCdef"
if GKYL_HAVE_CUDA then
   _ = require "Updater.momentCalcData._DistFuncCdefDevice"
end

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max", ["tensor"] = "Tensor" }

local _M = {}

-- select function to compute specified moment
function _M.selectMomCalc(mom, basisNm, CDIM, VDIM, polyOrder, calcOnDevice)
   local funcNm = nil
   if not calcOnDevice then
      funcNm = string.format("MomentCalc%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], mom, polyOrder)
   else
      funcNm = string.format("calcMom%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], mom, polyOrder)
   end
   return ffi.C[funcNm]
end

-- select function to compute specified gyrokinetic moment
function _M.selectGkMomCalc(mom, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("GkMomentCalc%dx%dv%s_%s_P%d", CDIM, VDIM, basisNmMap[basisNm], string.sub(mom,3), polyOrder)
   return ffi.C[funcNm]
end

-- select function to compute integrated moments
function _M.selectIntMomCalc(basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("IntMomentCalc%dx%dv%s_P%d", CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- With piecewise linear selfPrimMoments needs the star moments.
function _M.selectStarM0Calc(dir, kinSpecies, basisNm, CDIM, VDIM, applyPos)
   local posString = ""
   if applyPos then posString = "Positivity" end
   local funcNm = ""
   if dir == 1 then
      funcNm = string.format("%sM0Star%s%dx%dv%s_VX", kinSpecies, posString, CDIM, VDIM, basisNmMap[basisNm])
   elseif dir == 2 then                                                      
      funcNm = string.format("%sM0Star%s%dx%dv%s_VY", kinSpecies, posString, CDIM, VDIM, basisNmMap[basisNm])
   elseif dir == 3 then                                                      
      funcNm = string.format("%sM0Star%s%dx%dv%s_VZ", kinSpecies, posString, CDIM, VDIM, basisNmMap[basisNm])
   end
   return ffi.C[funcNm]
end
function _M.selectStarM1iM2Calc(kinSpecies, basisNm, CDIM, VDIM)
   local funcNm = string.format("%sM1iM2Star%dx%dv%s", kinSpecies, CDIM, VDIM, basisNmMap[basisNm])
   return ffi.C[funcNm]
end

-- selfPrimMoments also need kernels to compute certain integrals along velocity boundaries.
function _M.selectBoundaryFintegral(dir, kinSpecies, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = ""
   if dir == 1 then
      funcNm = string.format("%sBoundaryIntegral%dx%dv%s_F_VX_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   elseif dir == 2 then
      funcNm = string.format("%sBoundaryIntegral%dx%dv%s_F_VY_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   elseif dir == 3 then
      funcNm = string.format("%sBoundaryIntegral%dx%dv%s_F_VZ_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   end
   return ffi.C[funcNm]
end
function _M.selectBoundaryVFintegral(dir, kinSpecies, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = ""
   if dir == 1 then
      funcNm = string.format("%sBoundaryIntegral%dx%dv%s_vF_VX_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   elseif dir == 2 then
      funcNm = string.format("%sBoundaryIntegral%dx%dv%s_vF_VY_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   elseif dir == 3 then
      funcNm = string.format("%sBoundaryIntegral%dx%dv%s_vF_VZ_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   end
   return ffi.C[funcNm]
end

return _M
