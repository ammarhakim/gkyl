-- Gkyl ------------------------------------------------------------------------
--
-- Incompressible Euler equation in 2D
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local FluidSpecies = require "App.Species.FluidSpecies"
local IncompEulerEq = require "Eq.IncompEuler"
local Updater = require "Updater"
local DataStruct = require "DataStruct"

local IncompEulerSpecies = Proto(FluidSpecies)

function IncompEulerSpecies:fullInit(appTbl)
   IncompEulerSpecies.super.fullInit(self, appTbl)

   self.nMoments = 1
end

function IncompEulerSpecies:createSolver(hasE, hasB)
   -- create updater to advance solution by one time-step
   local eqn = IncompEulerEq {
      onGrid = self.grid,
      basis = self.basis,
   }

   self.solver = Updater.HyperDisCont {
      onGrid = self.grid,
      basis = self.basis,
      cfl = self.cfl,
      equation = eqn,
   }
end

-- nothing to calculate, just copy
function IncompEulerSpecies:calcCouplingMoments(tCurr, dt, rkIdx)
   local fIn = self:rkStepperFields()[rkIdx]
   self.couplingMoments:copy(fIn)
end

function IncompEulerSpecies:fluidMoments()
   return { self.couplingMoments }
end

-- for interfacing with GkField
function IncompEulerSpecies:getNumDensity(rkIdx)
   if rkIdx == nil then return self.couplingMoments 
   else 
      self.couplingMoments:copy(self:rkStepperFields()[rkIdx])
      return self.couplingMoments
   end
end

return IncompEulerSpecies
