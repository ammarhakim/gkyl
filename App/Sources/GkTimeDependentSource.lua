local SourceBase     = require "App.Sources.SourceBase"
local DataStruct     = require "DataStruct"
local ffi            = require "ffi"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"
local Projection     = require "App.Projection"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"

local GkTimeDependentSource = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkTimeDependentSource:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkTimeDependentSource:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.
   self.density = assert(tbl.density, "App.GkTimeDependentSource: must specify density profile of source in 'density'.")
   self.temperature = assert(tbl.temperature, "App.GkTimeDependentSource: must specify temperature profile of source in 'density'.")
   self.power = assert(tbl.power, "App.GkTimeDependentSource: must specify source power in 'power'.")
   self.profile = Projection.GkProjection.MaxwellianProjection {
                   density = self.density,
                   temperature = self.temperature,
                   power = self.power,
                  }
   self.tmEvalSrc = 0.0
end

function GkTimeDependentSource:setName(nm)
   self.name = nm
end
function GkTimeDependentSource:setSpeciesName(nm)
   self.speciesName = nm
end
function GkTimeDependentSource:setConfBasis(basis)
   self.confBasis = basis
end
function GkTimeDependentSource:setConfGrid(grid)
   self.confGrid = grid
end

function GkTimeDependentSource:advance(tCurr, fIn, species, fRhsOut)
   self.timeDependence = species[self.speciesName].sourceTimeDependence
   Mpi.Barrier(self.confGrid:commSet().sharedComm)
   fRhsOut:accumulate(self.timeDependence(tCurr), species[self.speciesName].fSource)
end

function GkTimeDependentSource:srcTime()
   return self.tmEvalSource
end

return GkTimeDependentSource
