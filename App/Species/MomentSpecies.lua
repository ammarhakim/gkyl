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
end

function MomentSpecies:initDist()
   local fld = self.moments[1]
   local grid = self.grid
   local meqn = self.equation:numEquations()
   local xn = Lin.Vec(fld:ndim())
   
   local fItr = fld:get(1)      
   local indexer = fld:genIndexer()
   for idx in fld:localExtRangeIter() do
      grid:setIndex(idx)
      grid:cellCenter(xn)
      -- evaluate supplied IC function
      local v = { self.initFunc(0.0, xn) } -- braces around function put return values in table
      -- set values in cell
      fld:fill(indexer(idx), fItr) -- pointer to data in cell
      for c = 1, meqn do fItr[c] = v[c] end
   end   
end

function MomentSpecies:forwardEuler(tCurr, dt, species, emIn, inIdx, outIdx)
   -- does nothing: perhaps when DG is supported this will need to be
   -- modified
end

function MomentSpecies:write(tm)
   if self.evolve then
      if self.diagIoTrigger(tm) then
	 self.moments[1]:write(string.format("%s_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame)
	 self.diagIoFrame = self.diagIoFrame + 1
      end
   else
      if self.diagIoFrame == 0 then
	 self.moments[1]:write(string.format("%s_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame)
	 self.diagIoFrame = self.diagIoFrame + 1
      end
   end
end

return MomentSpecies
