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

   local ndim = tbl.onGrid:ndim()
  
   self.slvr = nil
   if ndim == 1 then 
      self.slvr = FemParPoisson {
        onGrid = tbl.onGrid,
        basis = tbl.basis,
        bcBack = tbl.bcBack,
        bcFront = tbl.bcFront,
        periodicDirs = tbl.periodicDirs,
        laplacianWeight = tbl.laplacianWeight,
        modifierConstant = tbl.modifierConstant,
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
        laplacianWeight = tbl.laplacianWeight, 
        modifierConstant = tbl.modifierConstant,
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
        laplacianWeight = tbl.laplacianWeight, 
        modifierConstant = tbl.modifierConstant,
        zContinuous = tbl.zContinuous,
      }
   else 
      assert(false, "Updater.FemPoisson: Requires ndim<=3")
   end   

   return self
end


---- advance method
function FemPoisson:_advance(tCurr, dt, inFld, outFld) 
   return self.slvr:advance(tCurr, dt, inFld, outFld)
end

return FemPoisson
