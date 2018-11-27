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
local Mpi = require "Comm.Mpi"

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
      charge = self.charge,
      positivity = self.positivity,
   }

   self.solver = Updater.HyperDisCont {
      onGrid = self.grid,
      basis = self.basis,
      cfl = self.cfl,
      equation = eqn,
   }

   if self.positivity then 
      self.posRescaler = Updater.PositivityRescale {
         onGrid = self.grid,
         basis = self.basis,
      }
   end
end

-- nothing to calculate, just copy
function IncompEulerSpecies:calcCouplingMoments(tCurr, rkIdx)
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

function IncompEulerSpecies:suggestDt()
   -- loop over local region 
   local grid = self.grid
   self.dt[0] = GKYL_MAX_DOUBLE
   local tId = grid:subGridSharedId() -- local thread ID
   local localRange = self.cflRateByCell:localRange()
   local localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = localRange, numSplit = grid:numSharedProcs() }

   for idx in localRangeDecomp:colMajorIter(tId) do
      -- calculate local min dt from local cflRates
      self.cflRateByCell:fill(self.cflRateIdxr(idx), self.cflRatePtr)
      self.dt[0] = math.min(self.dt[0], self.cfl/self.cflRatePtr:data()[0])
   end

   -- all reduce to get global min dt
   Mpi.Allreduce(self.dt, self.dtGlobal, 1, Mpi.DOUBLE, Mpi.MIN, grid:commSet().comm)

   return math.min(self.dtGlobal[0], GKYL_MAX_DOUBLE)
end

function IncompEulerSpecies:clearCFL()
   -- clear cflRateByCell for next cfl calculation
   self.cflRateByCell:clear(0.0)
end

return IncompEulerSpecies
