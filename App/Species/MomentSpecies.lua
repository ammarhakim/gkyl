-- Gkyl ------------------------------------------------------------------------
--
-- Species object constructed from moment equations
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct = require "DataStruct"
local Lin = require "Lib.Linalg"
local LinearTrigger = require "Lib.LinearTrigger"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local FluidSpecies = require "App.Species.FluidSpecies"
local Time = require "Lib.Time"
local Updater = require "Updater"
local xsys = require "xsys"

-- Species object treated as moment equations
local MomentSpecies = Proto(FluidSpecies)

-- add constants to object indicate various supported boundary conditions
local SP_BC_OPEN = 1
local SP_BC_COPY = SP_BC_OPEN
local SP_BC_WALL = 4

MomentSpecies.bcOpen = SP_BC_OPEN -- open BCs
MomentSpecies.bcCopy = SP_BC_COPY -- copy BCs
MomentSpecies.bcWall = SP_BC_WALL -- wall BCs

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function MomentSpecies:fullInit(appTbl)
   MomentSpecies.super.fullInit(self, appTbl)

   self.equation = self.tbl.equation -- equation system to evolve
   self.nMoments = self.tbl.equation:numEquations()
   self.nGhost = 2 -- we need two ghost-cells

   self.limiter = self.tbl.limiter and self.tbl.limiter or "monotonized-centered"
   self.hyperSlvr = {} -- list of solvers
end

function MomentSpecies:createSolver(hasE, hasB)
   local ndim = self.grid:ndim()
   for d = 1, ndim do
      self.hyperSlvr[d] = Updater.WavePropagation {
	 onGrid = self.grid,
	 equation = self.equation,
	 limiter = self.limiter,
	 cfl = self.cfl,
	 updateDirections = {d}
      }
   end
end

function MomentSpecies:forwardEuler(tCurr, dt, species, emIn, inIdx, outIdx)
   -- does nothing: perhaps when DG is supported this will need to be
   -- modified

   return true, GKYL_MAX_DOUBLE
end

function MomentSpecies:appendBoundaryConditions(dir, edge, bcType)
   local function bcCopyFunc(...) return self:bcCopyFunc(...) end

   if bcType == SP_BC_COPY then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, edge, { bcCopyFunc }))
   else
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, edge, bcType ))
   end
end

function MomentSpecies:updateInDirection(dir, tCurr, dt, fIn, fOut)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   if self.evolve then
      status, dtSuggested = self.hyperSlvr[dir]:advance(tCurr, dt, {fIn}, {fOut})
      self:applyBc(tCurr, dt, fOut)
   else
      fOut:copy(fIn)
   end
   return status, dtSuggested
end

function MomentSpecies:totalSolverTime()
   local tm = 0.0
   for d = 1, self.grid:ndim() do
      tm = tm + self.hyperSlvr[d].totalTime
   end
   return tm
end

return MomentSpecies
