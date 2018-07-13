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

   local function contains(table, element)
     for _, value in pairs(table) do
       if value == element then
         return true
       end
     end
     return false
   end

   self.grid = tbl.onGrid
   self.ndim = self.grid:ndim()
   self.basis = tbl.basis
   local ndim = self.ndim
   self.laplacianWeight = tbl.laplacianWeight
   self.modifierConstant = tbl.modifierConstant
   self.zContinuous = tbl.zContinuous
   self.periodicDirs = tbl.periodicDirs   
   self.bcLeft = tbl.bcLeft
   self.bcRight = tbl.bcRight
   self.bcTop = tbl.bcTop
   self.bcBottom = tbl.bcBottom
   self.bcBack = tbl.bcBack
   self.bcFront = tbl.bcFront

   -- make sure BCs are specified consistently
   if contains(self.periodicDirs,1) == false and self.modifierConstant == 1.0 and self.laplacianWeight == 0.0 and self.bcLeft == nil and self.bcRight == nil then
     -- if not periodic, use neumann by default in this case (discont-to-cont projection)
     self.bcLeft = { T = "N", V = 0.0 }
     self.bcRight = { T = "N", V = 0.0 }
   elseif contains(self.periodicDirs,1) then
     assert(self.bcLeft == nil and self.bcRight == nil, "Cannot specify BCs if direction is periodic")
   end
   if contains(self.periodicDirs,2) == false and self.modifierConstant == 1.0 and self.laplacianWeight == 0.0 and self.bcTop == nil and self.bcBottom == nil then
     -- if not periodic, use neumann by default in this case (discont-to-cont projection)
     self.bcTop = { T = "N", V = 0.0 }
     self.bcBottom = { T = "N", V = 0.0 }
   elseif contains(self.periodicDirs,2) then
     assert(self.bcTop == nil and self.bcBottom == nil, "Cannot specify BCs if direction is periodic")
   end
   if contains(self.periodicDirs,3) == false and self.modifierConstant == 1.0 and self.laplacianWeight == 0.0 and self.bcBack == nil and self.bcFront == nil then
     -- if not periodic, use neumann by default in this case (discont-to-cont projection)
     self.bcBack = { T = "N", V = 0.0 }
     self.bcFront = { T = "N", V = 0.0 }
   elseif contains(self.periodicDirs,3) then
     assert(self.bcBack == nil and self.bcFront == nil, "Cannot specify BCs if direction is periodic")
   end
  
   self.slvr = nil
   if ndim == 1 then 
      self.slvr = FemParPoisson {
        onGrid = self.grid,
        basis = self.basis,
        bcBack = self.bcLeft,
        bcFront = self.bcRight,
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
        bcLeft = self.bcLeft,
        bcRight = self.bcRight,
        bcBottom = self.bcBottom,
        bcTop = self.bcTop,
        periodicDirs = self.periodicDirs,
        laplacianWeight = self.laplacianWeight, 
        modifierConstant = self.modifierConstant,
      }
   elseif ndim == 3 then
      self.slvr = FemPerpPoisson {
        onGrid = self.grid,
        basis = self.basis,
        bcLeft = self.bcLeft,
        bcRight = self.bcRight,
        bcBottom = self.bcBottom,
        bcTop = self.bcTop,
        periodicDirs = self.periodicDirs,
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
