-- Gkyl ------------------------------------------------------------------------
--
-- Species that evolves via passive advection (and optional diffusion)
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct    = require "DataStruct"
local FluidSpecies  = require "App.Species.FluidSpecies"
local PassiveAdvectionEq = require "Eq.PassiveAdvection"
local LinearDecomp  = require "Lib.LinearDecomp"
local Updater       = require "Updater"
local Mpi           = require "Comm.Mpi"
local Proto         = require "Lib.Proto"

local PassiveAdvectionSpecies = Proto(FluidSpecies)

function PassiveAdvectionSpecies:createGrid(...)
   PassiveAdvectionSpecies.super.createGrid(self, ...)

   self.nMoments = self.ndim + 1
end

function PassiveAdvectionSpecies:createSolver()
   -- Run the FluidSpecies 'createSolver()' to initialize the
   -- collisions (diffusion) solver.
   PassiveAdvectionSpecies.super.createSolver(self)

   -- Create updater to advance solution by one time-step.
   self.equation = PassiveAdvectionEq {
      onGrid     = self.grid,
      basis      = self.basis,
      nMoments   = self.nMoments,
   }

   self.solver = Updater.HyperDisCont {
      onGrid             = self.grid,
      basis              = self.basis,
      cfl                = self.cfl,
      equation           = self.equation,
      zeroFluxDirections = self.zeroFluxDirections,
      clearOut = false,
   }
end

-- Nothing to calculate, just copy.
function PassiveAdvectionSpecies:calcCouplingMomentsEvolve(tCurr, rkIdx)
   local fIn = self:rkStepperFields()[rkIdx]
   self.couplingMoments:copy(fIn)
end

function PassiveAdvectionSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self.tCurr = tCurr
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   -- clear RHS, because HyperDisCont was set up with clearOut = false
   fRhsOut:clear(0.0)

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      local em = emIn[1]:rkStepperFields()[inIdx]
      self.solver:advance(tCurr, {fIn, em}, {fRhsOut})
   end

   -- Perform the collision (diffusion) update.
   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
         c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
         c:advance(tCurr, fIn, species, fRhsOut)
      end
   end

   for _, src in pairs(self.sources) do src:advance(tCurr, fIn, species, fRhsOut) end
end

function PassiveAdvectionSpecies:fluidMoments()
   return { self.couplingMoments }
end

return PassiveAdvectionSpecies
