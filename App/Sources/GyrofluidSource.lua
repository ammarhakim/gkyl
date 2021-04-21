-- Gkyl ------------------------------------------------------------------------
--
-- Source terms in the gyrofluid model. They correspond to
-- mass density, momentum density, total particle energy density and
-- perpendicular pressure (divided by magnetic field) source rates, all multiplied
-- by the Jacobian.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase     = require "App.Sources.SourceBase"
local DataStruct     = require "DataStruct"
local Mpi            = require "Comm.Mpi"
local Projection     = require "App.Projection.GyrofluidProjection"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"

local GyrofluidSource = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GyrofluidSource:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GyrofluidSource:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.denFunc   = assert(tbl.density, "App.GyrofluidSource: must specify density profile of source in 'density'.")
   self.uParFunc  = tbl.driftSpeed or function(t, xn) return 0. end
   self.TperpFunc = assert(tbl.perpendicularTemperature, "App.GyrofluidSource: must specify perpendicular temperature profile of source in 'perpendicularTemperature'.")
   self.TparFunc  = assert(tbl.parallelTemperature, "App.GyrofluidSource: must specify parallel temperature profile of source in 'parallelTemperature'.")

   self.timeDependence = tbl.timeDependence or function (t) return 1. end

   self.nMoments = 3+1

   self.gfProj = Projection.GyrofluidProjection {
      density = self.denFunc,
      driftSpeed = self.uParFunc, 
      parallelTemperature = self.TparFunc,
      perpendicularTemperature = self.TperpFunc,
   }
   self.gfProj:fullInit(speciesTbl)

   self.tmEvalSrc = 0.0
end

function GyrofluidSource:setName(nm) self.name = nm end
function GyrofluidSource:setSpeciesName(nm) self.speciesName = nm end
function GyrofluidSource:setConfBasis(basis) self.basis = basis end
function GyrofluidSource:setConfGrid(grid) self.grid = grid end
function GyrofluidSource:setCfl(cfl) self.cfl = cfl end

function GyrofluidSource:createSolver(externalField)
   -- Thermal speed squared times collisionality, summed over species.
   self.momSource = DataStruct.Field {
      onGrid        = self.grid,
      numComponents = self.basis:numBasis()*self.nMoments,
      ghost         = {1, 1},
   }

   self.gfProj:advance(0., {externalField}, {self.momSource})
   Mpi.Barrier(self.grid:commSet().sharedComm)

end

function GyrofluidSource:advance(tCurr, momIn, species, momRhsOut)
   local tm = Time.clock()
   Mpi.Barrier(self.grid:commSet().sharedComm)
   momRhsOut:accumulate(self.timeDependence(tCurr), self.momSource)
   self.tmEvalSrc = self.tmEvalSrc + Time.clock() - tm
end

function GyrofluidSource:write(tm, frame, species)
   self.momSource:write(string.format("%s_source_0.bp", self.speciesName), tm, frame, true)
end

function GyrofluidSource:srcTime()
   return self.tmEvalSrc
end

return GyrofluidSource
