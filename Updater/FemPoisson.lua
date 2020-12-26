-- Gkyl ------------------------------------------------------------------------
--
-- Wraps FemParPoisson and FemPerpPoisson for common cases in GK
-- To be generalized further...
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto          = require "Lib.Proto"
local UpdaterBase    = require "Updater.Base"
local FemParPoisson  = require "Updater.FemParPoisson"
local FemPerpPoisson = require "Updater.FemPerpPoisson"
local DataStruct     = require "DataStruct"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local CartFieldBinOp = require "Updater.CartFieldBinOp"
local xsys           = require "xsys"

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

   self.grid  = tbl.onGrid
   self.ndim  = self.grid:ndim()
   self.basis = tbl.basis
   local ndim = self.ndim
   self.zContinuous = xsys.pickBool(tbl.zContinuous, true)

   -- Boundary conditions.
   self.bcLower = tbl.bcLower
   self.bcUpper = tbl.bcUpper

   -- Metric coefficients.
   self.gxx = tbl.gxx
   self.gxy = tbl.gxy
   self.gyy = tbl.gyy

   self.smooth = xsys.pickBool(tbl.smooth, false)

   self.slvr = nil
   if ndim == 1 then 
      self.slvr = FemParPoisson {
         onGrid  = self.grid,
         basis   = self.basis,
         bcLower = self.bcLower,
         bcUpper = self.bcUpper,
         smooth  = self.smooth
      }
      -- Set up weak division operator for special case when solve is algebraic.
      self.weakDivide = CartFieldBinOp {
         onGrid    = self.grid,
         weakBasis = self.basis,
         operation = "Divide",
         onGhosts  = true,
      }
   elseif ndim == 2 then
      self.slvr = FemPerpPoisson {
        onGrid       = self.grid,
        basis        = self.basis,
        bcLower      = self.bcLower,
        bcUpper      = self.bcUpper,
        constStiff   = self.constStiff,
        gxx          = self.gxx,
        gxy          = self.gxy,
        gyy          = self.gyy,
        smooth       = self.smooth
      }
   elseif ndim == 3 then
      self.slvr = FemPerpPoisson {
        onGrid       = self.grid,
        basis        = self.basis,
        bcLower      = self.bcLower,
        bcUpper      = self.bcUpper,
        zContinuous  = self.zContinuous,
        constStiff   = self.constStiff,
        gxx          = self.gxx,
        gxy          = self.gxy,
        gyy          = self.gyy,
        smooth       = self.smooth
      }
   else 
      assert(false, "Updater.FemPoisson: Requires ndim<=3")
   end   

   return self
end


-- Advance method.
function FemPoisson:_advance(tCurr, inFld, outFld) 
   -- Special case where solve is just algebraic.
   if self.ndim == 1 and not self.zContinuous and self.slvr._hasLaplacian == false then
      local src = inFld[1]
      local sol = outFld[1]

      self.weakDivide:advance(0, {self.slvr:getModifierWeight(), src}, {sol})
   else
      self.slvr:advance(tCurr, inFld, outFld)
   end
end

function FemPoisson:setLaplacianWeight(weight)
   self.slvr:setLaplacianWeight(weight)
end
function FemPoisson:setModifierWeight(weight)
   self.slvr:setModifierWeight(weight)
end

function FemPoisson:getLaplacianWeight()
   return self.slvr:getLaplacianWeight()
end
function FemPoisson:getModifierWeight()
   return self.slvr:getModifierWeight()
end

return FemPoisson
