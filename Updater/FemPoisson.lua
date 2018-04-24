-- Gkyl ------------------------------------------------------------------------
--
-- Wraps FemParPoisson and FemPerpPoisson for common cases in GK
-- To be generalized further...
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local FemParPoisson = require "Updater.FemParPoisson"
local FemPerpPoisson = require "Updater.FemPerpPoisson"
local DataStruct = require "DataStruct"
local ProjectOnBasis = require "Updater.ProjectOnBasis"

-- FEM Poisson solver updater object
local FemPoisson = Proto(UpdaterBase)

function FemPoisson:init(tbl)
   FemPoisson.super.init(self, tbl)

   self.ndim = tbl.onGrid:ndim()
   self.basis = tbl.basis
   local ndim = self.ndim
   self.laplacianWeight = tbl.laplacianWeight
   self.modifierConstant = tbl.modifierConstant
   self.zContinuous = tbl.zContinuous
  
   self.slvr = nil
   if ndim == 1 then 
      self.slvr = FemParPoisson {
        onGrid = tbl.onGrid,
        basis = self.basis,
        bcBack = tbl.bcBack,
        bcFront = tbl.bcFront,
        periodicDirs = tbl.periodicDirs,
        laplacianWeight = self.laplacianWeight,
        modifierConstant = self.modifierConstant,
      }
      -- set up constant dummy field
      self.unitWeight = DataStruct.Field {
           onGrid = tbl.onGrid,
           numComponents = self.basis:numBasis(),
           ghost = {1, 1},
      }
      local initUnit = ProjectOnBasis {
         onGrid = tbl.onGrid,
         basis = self.basis,
         evaluate = function (t,xn)
                       return 1.0
                    end
      }
      initUnit:advance(0.,0.,{},{self.unitWeight})
   elseif ndim == 2 then
      self.slvr = FemPerpPoisson {
        onGrid = tbl.onGrid,
        basis = self.basis,
        bcLeft = tbl.bcLeft,
        bcRight = tbl.bcRight,
        bcBottom = tbl.bcBottom,
        bcTop = tbl.bcTop,
        periodicDirs = tbl.periodicDirs,
        laplacianWeight = self.laplacianWeight, 
        modifierConstant = self.modifierConstant,
      }
   elseif ndim == 3 then
      self.slvr = FemPerpPoisson {
        onGrid = tbl.onGrid,
        basis = self.basis,
        bcLeft = tbl.bcLeft,
        bcRight = tbl.bcRight,
        bcBottom = tbl.bcBottom,
        bcTop = tbl.bcTop,
        periodicDirs = tbl.periodicDirs,
        laplacianWeight = self.laplacianWeight, 
        modifierConstant = self.modifierConstant,
        zContinuous = self.zContinuous,
      }
   else 
      assert(false, "Updater.FemPoisson: Requires ndim<=3")
   end   

   return self
end


---- advance method
function FemPoisson:_advance(tCurr, dt, inFld, outFld) 
   -- special case where solve is just algebraic
   if self.ndim == 1 and not self.zContinuous and self.laplacianWeight == 0.0 then
      assert(self.basis:polyOrder() == 1, "1D pseudo-solve only works for p=1 for now")
      local src = inFld[1]
      local sol = outFld[1]
      if inFld[2] == nil and inFld[3] == nil then 
         -- if no stiffWeight and no massWeight, just do scalar division
         if inFld[4] ~= nil then
            sol:combine(1.0/inFld[4], src)
         else 
            sol:combine(1.0/self.modifierConstant, src)
         end
         return true, GKYL_MAX_DOUBLE
      end
      local stiffWeight = inFld[2] or self.unitWeight
      local massWeight = inFld[3] or self.unitWeight
      stiffWeight:scale(self.laplacianWeight)
      massWeight:scale(self.modifierConstant)
      local localRange = src:localExtRange()
      local srcIndexer = src:genIndexer()
      local srcPtr = src:get(1)
      local massWeightIndexer = massWeight:genIndexer()
      local massWeightPtr = massWeight:get(1)
      local solIndexer = sol:genIndexer()
      local solPtr = sol:get(1)
 
      -- calculate sol := src / massWeight (via weak division)
      for idx in localRange:colMajorIter() do
         src:fill(srcIndexer(idx), srcPtr)
         massWeight:fill(massWeightIndexer(idx), massWeightPtr)
         sol:fill(solIndexer(idx), solPtr)

         solPtr:data()[1] = (math.sqrt(2)*srcPtr:data()[1]*massWeightPtr:data()[1]-math.sqrt(2)*srcPtr:data()[0]*massWeightPtr:data()[0])/(massWeightPtr:data()[1]^2-massWeightPtr:data()[0]^2)
         solPtr:data()[2] = (math.sqrt(2)*srcPtr:data()[0]*massWeightPtr:data()[1]-math.sqrt(2)*massWeightPtr:data()[0]*srcPtr:data()[1])/(massWeightPtr:data()[1]^2-massWeightPtr:data()[0]^2)
      end
      return true, GKYL_MAX_DOUBLE
   else
      return self.slvr:advance(tCurr, dt, inFld, outFld)
   end
end

return FemPoisson
