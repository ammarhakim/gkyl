-- Gkyl ------------------------------------------------------------------------
--
-- Source term in a DG fluid model.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local DataStruct = require "DataStruct"
local Mpi        = require "Comm.Mpi"
local Projection = require "App.Projection.FluidProjection"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"

local FluidSource = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function FluidSource:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function FluidSource:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.srcIn = assert(tbl.source, "App.FluidSource: must specify source profile in 'source'.")

   self.timeDependence = tbl.timeDependence or function (t) return 1. end

   self.nMoments = 1

   self.timers = {accumulateSrc = 0.0}
end

function FluidSource:setName(nm) self.name = nm end
function FluidSource:setSpeciesName(nm) self.speciesName = nm end
function FluidSource:setConfBasis(basis) self.basis = basis end
function FluidSource:setConfGrid(grid) self.grid = grid end
function FluidSource:setCfl(cfl) self.cfl = cfl end

function FluidSource:createSolver(mySpecies, externalField)
   -- Source rate in each moment equation.
   self.momSource = mySpecies:allocVectorMoment(mySpecies.nMoments)

   local fluidProj
   if type(self.srcIn) == "function" then
       fluidProj = Projection.FunctionProjection {
         func = function(t, zn) return self.srcIn(t, zn) end
      }
   elseif type(self.srcIn) == "string" then
      fluidProj = Projection.ReadInput {
         inputFile = self.srcIn,
      }
   end
   fluidProj:fullInit(mySpecies)

   fluidProj:advance(0., {externalField}, {self.momSource})
   Mpi.Barrier(self.grid:commSet().sharedComm)

   if mySpecies.positivityRescale then
      mySpecies.posRescaler:advance(0.0, {self.momSource}, {self.momSource})
   end
end

function FluidSource:advance(tCurr, momIn, species, momRhsOut)
   local tm = Time.clock()

   Mpi.Barrier(self.grid:commSet().sharedComm)
   momRhsOut:accumulate(self.timeDependence(tCurr), self.momSource)

   self.timers.accumulateSrc = self.timers.accumulateSrc + Time.clock() - tm
end

function FluidSource:write(tm, frame, species)
   if tm == 0.0 then
      self.momSource:write(string.format("%s_source_0.bp", self.speciesName), tm, frame, true)
   end
end

function FluidSource:srcTime()
   return self.timers.accumulateSrc
end

return FluidSource
