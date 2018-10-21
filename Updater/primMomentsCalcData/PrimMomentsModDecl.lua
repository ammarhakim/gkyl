-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into C++ kernel functions for computing primitive moments u and vtSq
-- for cross collisions based on CDIM, VDIM, basis functions and polyOrder.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local _ = require "Updater.primMomentsCalcData._PrimMomentsCdef"

-- map of basis function name -> function encoding
local basisNmMap = { ["serendipity"] = "Ser", ["maximal-order"] = "Max" }

local _M = {}

-- select kernel function to compute the primitive moments for self-collision terms. 
function _M.selectSelfPrimMomentsCalc(kinSpecies, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("%sSelfPrimMoments%dx%dv%s_P%d", kinSpecies, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

-- With piecewise linear selfPrimMoments needs the star moments.
function _M.selectStarM0Calc(dir, kinSpecies, basisNm, CDIM, VDIM)
   local funcNm = ""
   if dir == 1 then 
      funcNm = string.format("%sM0Star%dx%dv%s_VX", kinSpecies, CDIM, VDIM, basisNmMap[basisNm])
   elseif dir == 2 then                                                          
      funcNm = string.format("%sM0Star%dx%dv%s_VY", kinSpecies, CDIM, VDIM, basisNmMap[basisNm])
   elseif dir == 3 then                                                          
      funcNm = string.format("%sM0Star%dx%dv%s_VZ", kinSpecies, CDIM, VDIM, basisNmMap[basisNm])
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

-- select kernel function to compute the cross primitive moments. 
function _M.selectCrossPrimMomentsCalc(collOp, collide, basisNm, CDIM, VDIM, polyOrder)
   local funcNm = string.format("CrossPrimMoments_%s%s_%dx%dv%s_P%d", collide, operator, CDIM, VDIM, basisNmMap[basisNm], polyOrder)
   return ffi.C[funcNm]
end

return _M
