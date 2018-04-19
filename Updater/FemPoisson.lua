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

-- FEM Poisson solver updater object
local FemPoisson = Proto(UpdaterBase)

function FemPoisson:init(tbl)
   FemPoisson.super.init(self, tbl)

   self.ndim = tbl.onGrid:ndim()
   local ndim = self.ndim
   self.laplacianWeight = tbl.laplacianWeight
   self.modifierConstant = tbl.modifierConstant
   self.zContinuous = tbl.zContinuous
  
   self.slvr = nil
   if ndim == 1 then 
      self.slvr = FemParPoisson {
        onGrid = tbl.onGrid,
        basis = tbl.basis,
        bcBack = tbl.bcBack,
        bcFront = tbl.bcFront,
        periodicDirs = tbl.periodicDirs,
        laplacianWeight = self.laplacianWeight,
        modifierConstant = self.modifierConstant,
      }
   elseif ndim == 2 then
      self.slvr = FemPerpPoisson {
        onGrid = tbl.onGrid,
        basis = tbl.basis,
        bcLeft = tbl.phiBcLeft,
        bcRight = tbl.phiBcRight,
        bcBottom = tbl.phiBcBottom,
        bcTop = tbl.phiBcTop,
        periodicDirs = tbl.periodicDirs,
        laplacianWeight = self.laplacianWeight, 
        modifierConstant = self.modifierConstant,
      }
   elseif ndim == 3 then
      self.slvr = FemPerpPoisson {
        onGrid = tbl.onGrid,
        basis = tbl.basis,
        bcLeft = tbl.phiBcLeft,
        bcRight = tbl.phiBcRight,
        bcBottom = tbl.phiBcBottom,
        bcTop = tbl.phiBcTop,
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
   if self.ndim == 1 and not self.zContinuous and self.laplacianWeight == 0.0 then
      local fin, fout = inFld[1], outFld[1]
      fout:combine(1/self.modifierConstant, fin)
      return true, GKYL_MAX_DOUBLE
   else
      return self.slvr:advance(tCurr, dt, inFld, outFld)
   end
end

return FemPoisson
