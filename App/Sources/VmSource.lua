local SourceBase     = require "App.Sources.SourceBase"
local DataStruct     = require "DataStruct"
local ffi            = require "ffi"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"
local Projection     = require "App.Projection"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"

local VmSource = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmTimeDependentSource:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmSource:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.
   self.density = assert(tbl.density, "App.VmSource: must specify density profile of source in 'density'.")
   self.temperature = assert(tbl.temperature, "App.VmSource: must specify temperature profile of source in 'density'.")
   self.power = tbl.power
   self.profile = Projection.VlasovProjection.MaxwellianProjection {
                   density = self.density,
                   temperature = self.temperature,
                   power = self.power,
                  }
   self.tmEvalSrc = 0.0
end

function VmSource:setName(nm)
   self.name = nm
end
function VmSource:setSpeciesName(nm)
   self.speciesName = nm
end
function VmSource:setConfBasis(basis)
   self.confBasis = basis
end
function VmSource:setConfGrid(grid)
   self.confGrid = grid
end

function VmSource:advance(tCurr, fIn, species, fRhsOut)
   self.timeDependence = species[self.speciesName].sourceTimeDependence
   Mpi.Barrier(self.confGrid:commSet().sharedComm)
   fRhsOut:accumulate(self.timeDependence(tCurr), species[self.speciesName].fSource)
end

function VmSource:srcTime()
   return self.tmEvalSource
end

return VmTimeDependentSource
