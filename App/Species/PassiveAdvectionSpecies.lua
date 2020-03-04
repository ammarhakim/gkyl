-- Gkyl ------------------------------------------------------------------------
--
-- Species that evolves via passive advection 
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

local SP_BC_ABSORB   = 1
local SP_BC_COPY     = 2
local SP_BC_ZEROFLUX = 3
PassiveAdvectionSpecies.bcAbsorb   = SP_BC_ABSORB   -- Absorb.
PassiveAdvectionSpecies.bcCopy     = SP_BC_COPY
PassiveAdvectionSpecies.bcZeroFlux = SP_BC_ZEROFLUX -- Zero flux.

function PassiveAdvectionSpecies:createGrid(...)
   PassiveAdvectionSpecies.super.createGrid(self, ...)

   self.nMoments = self.cdim + 1
end

function PassiveAdvectionSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- Need to wrap member functions so that self is passed.
   local function bcAbsorbFunc(...) return self:bcAbsorbFunc(...) end
   local function bcCopyFunc(...) return self:bcCopyFunc(...) end

   if bcType == SP_BC_ABSORB then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, { bcAbsorbFunc }, "pointwise"))
   elseif bcType == SP_BC_COPY then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, { bcCopyFunc }, "pointwise"))
   elseif bcType == SP_BC_ZEROFLUX then
      table.insert(self.zeroFluxDirections, dir)
   else
      assert(false, "PassiveAdvectionSpecies: Unsupported BC type!")
   end
end

function PassiveAdvectionSpecies:createSolver()
   -- Run the FluidSpecies 'createSolver()' to initialize the
   -- collisions (diffusion) solver.
   PassiveAdvectionSpecies.super.createSolver(self)

   -- Create updater to advance solution by one time-step.
   local eqn = PassiveAdvectionEq {
      onGrid     = self.grid,
      basis      = self.basis,
      positivity = self.positivity,
   }

   self.solver = Updater.HyperDisCont {
      onGrid             = self.grid,
      basis              = self.basis,
      cfl                = self.cfl,
      equation           = eqn,
      zeroFluxDirections = self.zeroFluxDirections,
   }
end

-- Nothing to calculate, just copy.
function PassiveAdvectionSpecies:calcCouplingMoments(tCurr, rkIdx)
   local fIn = self:rkStepperFields()[rkIdx]
   self.couplingMoments:copy(fIn)
end

function PassiveAdvectionSpecies:fluidMoments()
   return { self.couplingMoments }
end

-- For interfacing with GkField.
function PassiveAdvectionSpecies:getNumDensity(rkIdx)
end

function PassiveAdvectionSpecies:suggestDt(inIdx, outIdx)
   -- Loop over local region.
   local grid = self.grid
   self.dt[0] = GKYL_MAX_DOUBLE
   local tId              = grid:subGridSharedId() -- Local thread ID.
   local localRange       = self.cflRateByCell:localRange()
   local localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = localRange, numSplit = grid:numSharedProcs() }

   local fIn = self:rkStepperFields()[inIdx]
   local fRhsSurf = self:rkStepperFields()[outIdx]
   local fInPtr = fIn:get(1)
   local fRhsSurfPtr = fRhsSurf:get(1)
   local fIdxr = fIn:genIndexer()

   for idx in localRangeDecomp:rowMajorIter(tId) do
      -- Calculate local min dt from local cflRates.
      self.cflRateByCell:fill(self.cflRateIdxr(idx), self.cflRatePtr)
      self.dt[0] = math.min(self.dt[0], self.cfl/self.cflRatePtr:data()[0])
   end

   -- All reduce to get global min dt.
   Mpi.Allreduce(self.dt, self.dtGlobal, 1, Mpi.DOUBLE, Mpi.MIN, grid:commSet().comm)

   return math.min(self.dtGlobal[0], GKYL_MAX_DOUBLE)
end

function PassiveAdvectionSpecies:clearCFL()
   -- Clear cflRateByCell for next cfl calculation,
   self.cflRateByCell:clear(0.0)
end

return PassiveAdvectionSpecies
