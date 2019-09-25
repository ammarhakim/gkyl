-- Gkyl ------------------------------------------------------------------------
--
-- Incompressible Euler equation in 2D.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct    = require "DataStruct"
local FluidSpecies  = require "App.Species.FluidSpecies"
local IncompEulerEq = require "Eq.IncompEuler"
local LinearDecomp  = require "Lib.LinearDecomp"
local Updater       = require "Updater"
local Mpi           = require "Comm.Mpi"
local Proto         = require "Lib.Proto"

local IncompEulerSpecies = Proto(FluidSpecies)

local SP_BC_ABSORB   = 1
local SP_BC_COPY     = 2
local SP_BC_ZEROFLUX = 3
IncompEulerSpecies.bcAbsorb   = SP_BC_ABSORB   -- Absorb.
IncompEulerSpecies.bcCopy     = SP_BC_COPY
IncompEulerSpecies.bcZeroFlux = SP_BC_ZEROFLUX -- Zero flux.

function IncompEulerSpecies:fullInit(appTbl)
   IncompEulerSpecies.super.fullInit(self, appTbl)

   self.nMoments = 1
end

function IncompEulerSpecies:appendBoundaryConditions(dir, edge, bcType)
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
      assert(false, "IncompEulerSpecies: Unsupported BC type!")
   end
end

function IncompEulerSpecies:createSolver(hasE, hasB)
   -- Run the FluidSpecies 'createSolver()' to initialize the
   -- collisions (diffusion) solver.
   IncompEulerSpecies.super.createSolver(self)

   -- Create updater to advance solution by one time-step.
   local eqn = IncompEulerEq {
      onGrid     = self.grid,
      basis      = self.basis,
      charge     = self.charge,
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
function IncompEulerSpecies:calcCouplingMoments(tCurr, rkIdx)
   local fIn = self:rkStepperFields()[rkIdx]
   self.couplingMoments:copy(fIn)
end

function IncompEulerSpecies:fluidMoments()
   return { self.couplingMoments }
end

-- For interfacing with GkField.
function IncompEulerSpecies:getNumDensity(rkIdx)
   if rkIdx == nil then return self.couplingMoments 
   else 
      self.couplingMoments:copy(self:rkStepperFields()[rkIdx])
      return self.couplingMoments
   end
end

function IncompEulerSpecies:suggestDt()
   -- Loop over local region.
   local grid = self.grid
   self.dt[0] = GKYL_MAX_DOUBLE
   local tId              = grid:subGridSharedId() -- Local thread ID.
   local localRange       = self.cflRateByCell:localRange()
   local localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = localRange, numSplit = grid:numSharedProcs() }

   for idx in localRangeDecomp:rowMajorIter(tId) do
      -- Calculate local min dt from local cflRates.
      self.cflRateByCell:fill(self.cflRateIdxr(idx), self.cflRatePtr)
      self.dt[0] = math.min(self.dt[0], self.cfl/self.cflRatePtr:data()[0])
   end

   -- All reduce to get global min dt.
   Mpi.Allreduce(self.dt, self.dtGlobal, 1, Mpi.DOUBLE, Mpi.MIN, grid:commSet().comm)

   return math.min(self.dtGlobal[0], GKYL_MAX_DOUBLE)
end

function IncompEulerSpecies:clearCFL()
   -- Clear cflRateByCell for next cfl calculation,
   self.cflRateByCell:clear(0.0)
end

return IncompEulerSpecies
