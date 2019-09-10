-- Gkyl ------------------------------------------------------------------------
--
-- Incompressible Euler equation in 2D
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct = require "DataStruct"
local FluidSpecies = require "App.Species.FluidSpecies"
local IncompEulerEq = require "Eq.IncompEuler"
local LinearDecomp = require "Lib.LinearDecomp"
local Updater = require "Updater"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Lin = require "Lib.Linalg"

local IncompEulerSpecies = Proto(FluidSpecies)

local SP_BC_ABSORB = 1
local SP_BC_COPY = 2
local SP_BC_ZEROFLUX = 3
IncompEulerSpecies.bcAbsorb = SP_BC_ABSORB -- absorb 
IncompEulerSpecies.bcCopy = SP_BC_COPY
IncompEulerSpecies.bcZeroFlux = SP_BC_ZEROFLUX -- zero flux

function IncompEulerSpecies:fullInit(appTbl)
   IncompEulerSpecies.super.fullInit(self, appTbl)

   self.nMoments = 1
end

function IncompEulerSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- need to wrap member functions so that self is passed
   local function bcAbsorbFunc(...) return self:bcAbsorbFunc(...) end
   local function bcCopyFunc(...) return self:bcCopyFunc(...) end

   if bcType == SP_BC_ABSORB then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcAbsorbFunc }, "pointwise"))
   elseif bcType == SP_BC_COPY then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcCopyFunc }, "pointwise"))
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
      zeroFluxDirections = self.zeroFluxDirections,
   }

   if self.positivity then 
      self.posChecker = Updater.PositivityCheck {
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

function IncompEulerSpecies:suggestDt(inIdx, outIdx)
   -- loop over local region 
   local grid = self.grid
   self.dt[0] = GKYL_MAX_DOUBLE
   local tId = grid:subGridSharedId() -- local thread ID
   local localRange = self.cflRateByCell:localRange()
   local localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = localRange, numSplit = grid:numSharedProcs() }

   local fIn = self:rkStepperFields()[inIdx]
   local fRhs = self:rkStepperFields()[outIdx]
   local fInPtr = fIn:get(1)
   local fRhsPtr = fRhs:get(1)
   local fIdxr = fIn:genIndexer()
   local posDt = GKYL_MAX_DOUBLE
   local posIdx = Lin.IntVec(self.grid:ndim()) 

   for idx in localRangeDecomp:rowMajorIter(tId) do
     -- calculate local min dt from local cflRates
     self.cflRateByCell:fill(self.cflRateIdxr(idx), self.cflRatePtr)
     self.dt[0] = math.min(self.dt[0], self.cfl/self.cflRatePtr:data()[0])

     -- if using positivity, limit dt so that cell avg does not go negative, 
     -- but goes to exactly 0  
     if self.positivity then 
        fIn:fill(fIdxr(idx), fInPtr)
        fRhs:fill(fIdxr(idx), fRhsPtr)
        if fRhsPtr:data()[0]<0. and fInPtr:data()[0]>0. then
           -- calculate dt that will make cell avg go to exactly half its current value
           local dt = -fInPtr:data()[0]/fRhsPtr:data()[0]/2
           -- if this positivity-limited dt is smaller than the cfl-calculated dt, use it
           self.dt[0] = math.min(self.dt[0], dt)
        end
     end
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
