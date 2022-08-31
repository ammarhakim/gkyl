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
local GkFemPoisson = Proto(UpdaterBase)

function GkFemPoisson:init(tbl)
   GkFemPoisson.super.init(self, tbl)

   local function contains(table, element)
      for _, value in pairs(table) do
         if value == element then return true end
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
         onGrid = self.grid,   bcLower = self.bcLower,
         basis  = self.basis,  bcUpper = self.bcUpper,
         smooth = self.smooth,
         useG0 = tbl.useG0,
      }
      -- Set up weak division operator for special case when solve is algebraic.
      self.weakDivide = CartFieldBinOp {
         onRange   = self.slvr:getModifierWeight():localExtRange(),  onGhosts = true,
         weakBasis = self.basis,  operation = "Divide",
      }
   elseif ndim == 2 then
      self.slvr = FemPerpPoisson {
        onGrid     = self.grid,        gxx    = self.gxx,
        basis      = self.basis,       gxy    = self.gxy,
        bcLower    = self.bcLower,     gyy    = self.gyy,
        bcUpper    = self.bcUpper,     smooth = self.smooth,
        constStiff = self.constStiff,
        useG0 = tbl.useG0,
      }
   elseif ndim == 3 then
      self.slvr = FemPerpPoisson {
        onGrid      = self.grid,         constStiff = self.constStiff,
        basis       = self.basis,        gxx        = self.gxx,
        bcLower     = self.bcLower,      gxy        = self.gxy,
        bcUpper     = self.bcUpper,      gyy        = self.gyy,
        zContinuous = self.zContinuous,  smooth     = self.smooth,
        useG0 = tbl.useG0,
      }
   else 
      assert(false, "Updater.GkFemPoisson: Requires ndim<=3")
   end   

   -- Option to write development-related diagnostics.
   self.verbose = xsys.pickBool(tbl.verbose, false)

   return self
end

function GkFemPoisson:assemble(tCurr, inFld, outFld)
   -- Begin assembling the source vector and, if needed, the stiffness matrix.
   self.slvr:assemble(tCurr, inFld, outFld)
end

function GkFemPoisson:solve(tCurr, inFld, outFld) 
   -- Assuming the right-side vector (and if needed the stiffness matrix)
   -- has been assembled, this solves the linear problem.
   -- If the assembly initiated an MPI non-blocking reduce, this waits for it.
   self.slvr:solve(tCurr, inFld, outFld)
end

function GkFemPoisson:_advance(tCurr, inFld, outFld) 
   -- Advance method. This assembles the right-side source vector and, if needed,
   -- the stiffness matrix. Then it solves the linear problem.
   if self.ndim == 1 and not self.zContinuous and self.slvr._hasLaplacian == false then
      -- Special case where solve is just algebraic.
      local src = inFld[1]
      local sol = outFld[1]

      self.weakDivide:advance(0, {self.slvr:getModifierWeight(), src}, {sol})
   else
      self.slvr:advance(tCurr, inFld, outFld)
   end
end

function GkFemPoisson:_advanceOnDevice(tCurr, inFld, outFld)
   -- Advance method. This assembles the right-side source vector and, if needed,
   -- the stiffness matrix. Then it solves the linear problem.
   if self.ndim == 1 and not self.zContinuous and self.slvr._hasLaplacian == false then
      -- Special case where solve is just algebraic.
      local src = inFld[1]
      local sol = outFld[1]

      self.weakDivide:advance(0, {self.slvr:getModifierWeight(), src}, {sol})
   else
      self.slvr:advance(tCurr, inFld, outFld)
   end
end

function GkFemPoisson:setLaplacianWeight(weight)
   self.slvr:setLaplacianWeight(weight)
end
function GkFemPoisson:setModifierWeight(weight)
   self.slvr:setModifierWeight(weight)
end

function GkFemPoisson:getLaplacianWeight()
   return self.slvr:getLaplacianWeight()
end
function GkFemPoisson:getModifierWeight()
   return self.slvr:getModifierWeight()
end

function GkFemPoisson:printDevDiagnostics()
   if self.verbose then self.slvr:printDevDiagnostics() end
end

return GkFemPoisson
