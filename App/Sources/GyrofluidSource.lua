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

local SourceBase    = require "App.Sources.SourceBase"
local DataStruct    = require "DataStruct"
local Mpi           = require "Comm.Mpi"
local Projection    = require "App.Projection.GyrofluidProjection"
local DiagsApp      = require "App.Diagnostics.SpeciesDiagnostics"
local DiagsImplBase = require "App.Diagnostics.DiagnosticsImplBase"
local Updater       = require "Updater"
local Proto         = require "Lib.Proto"
local Time          = require "Lib.Time"

-- ............... IMPLEMENTATION OF DIAGNOSTICS ................. --
-- Diagnostics could be placed in a separate file if they balloon in
-- number. But if we only have one or two we can just place it here.

-- ~~~~ Source integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
local sourceDiagImpl = function()
   local _intSrc = Proto(DiagsImplBase)
   function _intSrc:fullInit(diagApp, mySpecies, fieldIn, srcIn)
      self.srcName  = string.gsub(srcIn.name, srcIn.speciesName.."_", "")
      self.field    = DataStruct.DynVector { numComponents = srcIn.nMoments }
      self.updaters = mySpecies.volIntegral.vector
      self.done     = false
   end
   function _intSrc:getType() return "integrated" end
   function _intSrc:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      self.updaters:advance(tm, {specIn.sources[self.srcName]:getSource()}, {self.field})
   end
   return {intSrc = _intSrc}
end

-- .................... END OF DIAGNOSTICS ...................... --

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

   self.timers = {accumulateSrc = 0.0}
end

function GyrofluidSource:setName(nm) self.name = self.speciesName.."_"..nm end
function GyrofluidSource:setSpeciesName(nm) self.speciesName = nm end
function GyrofluidSource:setConfBasis(basis) self.basis = basis end
function GyrofluidSource:setConfGrid(grid) self.grid = grid end
function GyrofluidSource:setCfl(cfl) self.cfl = cfl end

function GyrofluidSource:createSolver(mySpecies, externalField)
   -- Source rate in each moment equation.
   self.momSource = mySpecies:allocVectorMoment(mySpecies.nMoments)

   local gfProj = Projection.GyrofluidProjection {
      density                  = self.denFunc,
      driftSpeed               = self.uParFunc, 
      parallelTemperature      = self.TparFunc,
      perpendicularTemperature = self.TperpFunc,
   }
   gfProj:fullInit(mySpecies)

   gfProj:advance(0., {externalField}, {self.momSource})
   Mpi.Barrier(self.grid:commSet().sharedComm)

   self.momSource:write(string.format("%s_0.bp", self.name), 0.0, 0, true)
end

function GyrofluidSource:createDiagnostics(mySpecies, field)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = sourceDiagImpl()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

function GyrofluidSource:advance(tCurr, momIn, species, momRhsOut)
   local tm = Time.clock()

   Mpi.Barrier(self.grid:commSet().sharedComm)
   momRhsOut:accumulate(self.timeDependence(tCurr), self.momSource)

   self.timers.accumulateSrc = self.timers.accumulateSrc + Time.clock() - tm
end

function GyrofluidSource:getSource() return self.momSource end
function GyrofluidSource:srcTime() return self.timers.accumulateSrc end

return GyrofluidSource
