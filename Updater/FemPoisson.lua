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
local CartFieldBinOp = require "Updater.CartFieldBinOp"

-- FEM Poisson solver updater object
local FemPoisson = Proto(UpdaterBase)

function FemPoisson:init(tbl)
   FemPoisson.super.init(self, tbl)

   self.grid = tbl.onGrid
   self.ndim = self.grid:ndim()
   self.basis = tbl.basis
   local ndim = self.ndim
   self.laplacianWeight = tbl.laplacianWeight
   self.modifierConstant = tbl.modifierConstant
   self.zContinuous = tbl.zContinuous
  
   self.slvr = nil
   if ndim == 1 then 
      self.slvr = FemParPoisson {
        onGrid = self.grid,
        basis = self.basis,
        bcBack = tbl.bcBack,
        bcFront = tbl.bcFront,
        periodicDirs = tbl.periodicDirs,
        laplacianWeight = self.laplacianWeight,
        modifierConstant = self.modifierConstant,
      }
      -- set up constant dummy field
      self.unitWeight = DataStruct.Field {
           onGrid = self.grid,
           numComponents = self.basis:numBasis(),
           ghost = {1, 1},
      }
      local initUnit = ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = function (t,xn)
                       return 1.0
                    end
      }
      initUnit:advance(0.,0.,{},{self.unitWeight})
      -- set up weak division operator for special case when solve is algebraic
      self.weakDivide = CartFieldBinOp {
         onGrid = self.grid,
         weakBasis = self.basis,
         operation = "Divide",
         onGhosts = true,
      }
   elseif ndim == 2 then
      self.slvr = FemPerpPoisson {
        onGrid = self.grid,
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
        onGrid = self.grid,
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
      local massWeight = inFld[3] or self.unitWeight
      massWeight:scale(self.modifierConstant)

      self.weakDivide:advance(0, 0, {massWeight, src}, {sol})

      return true, GKYL_MAX_DOUBLE
   else
      return self.slvr:advance(tCurr, dt, inFld, outFld)
   end
end

return FemPoisson
