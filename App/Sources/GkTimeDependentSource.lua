local SourceBase     = require "App.Sources.SourceBase"
local DataStruct     = require "DataStruct"
local ffi            = require "ffi"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"
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
   print("made it! GK Init!")
   local tbl = self.tbl -- Previously stored table.
   self.source = Projection.GkProjection.MaxwellianProjection
   self.profile = function (t, zn) return tbl.source(t, zn, self) end
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
