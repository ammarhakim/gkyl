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
end

function MomentSpecies:updateInDirection(dir, tCurr, dt, fIn, fOut)
   if self.evolve then
      return self.hyperSlvr[dir]:advance(tCurr, dt, {fIn}, {fOut})
   else
      fOut:copy(fIn)
      return true, GKYL_MAX_DOUBLE
   end
end

return MomentSpecies
