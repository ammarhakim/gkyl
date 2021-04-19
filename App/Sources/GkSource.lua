-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Gyrokinetic source operator.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase     = require "App.Sources.SourceBase"
local DataStruct     = require "DataStruct"
local ffi            = require "ffi"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"
local Projection     = require "App.Projection"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"

local GkSource = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkSource:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkSource:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.
   self.density = assert(tbl.density, "App.GkSource: must specify density profile of source in 'density'.")
   self.temperature = assert(tbl.temperature, "App.GkSource: must specify temperature profile of source in 'density'.")
   self.power = tbl.power
   self.profile = Projection.GkProjection.MaxwellianProjection {
      density = self.density,
      temperature = self.temperature,
      power = self.power,
   }
   self.tmEvalSrc = 0.0
end

function GkSource:setName(nm)
   self.name = nm
end
function GkSource:setSpeciesName(nm)
   self.speciesName = nm
end
function GkSource:setConfBasis(basis)
   self.confBasis = basis
end
function GkSource:setConfGrid(grid)
   self.confGrid = grid
end

function GkSource:advance(tCurr, fIn, species, fRhsOut)
   local tm = Time.clock()
   self.timeDependence = species[self.speciesName].sourceTimeDependence
   Mpi.Barrier(self.confGrid:commSet().sharedComm)
   fRhsOut:accumulate(self.timeDependence(tCurr), species[self.speciesName].fSource)
   self.tmEvalSrc = self.tmEvalSrc + Time.clock() - tm
end

function GkSource:write(tm, frame, species)
   species.fSource:write(string.format("%s_fSource_0.bp", self.speciesName), tm, frame, true)
   if species.numDensitySrc then species.numDensitySrc:write(string.format("%s_srcM0_0.bp", self.speciesName), tm, frame) end
   if species.momDensitySrc then species.momDensitySrc:write(string.format("%s_srcM1_0.bp", self.speciesName), tm, frame) end
   if species.ptclEnergySrc then species.ptclEnergySrc:write(string.format("%s_srcM2_0.bp", self.speciesName), tm, frame) end
end

function GkSource:srcTime()
   return self.tmEvalSrc
end

return GkSource
