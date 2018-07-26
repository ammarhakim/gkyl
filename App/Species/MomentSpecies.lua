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

function MomentSpecies:forwardEuler(tCurr, dt, species, emIn, inIdx, outIdx)
   -- does nothing: perhaps when DG is supported this will need to be
   -- modified
end

return MomentSpecies
