-- Gkyl ------------------------------------------------------------------------
--
-- Source term in a DG fluid model.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase    = require "App.Sources.SourceBase"
local DataStruct    = require "DataStruct"
local Mpi           = require "Comm.Mpi"
local Projection    = require "App.Projection.FluidProjection"
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
   function _intSrc:fullInit(diagApp, mySpecies, srcIn)
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

function FluidSource:setSpeciesName(nm) self.speciesName = nm end
function FluidSource:setName(nm) self.name = self.speciesName.."_"..nm end
function FluidSource:setConfBasis(basis) self.basis = basis end
function FluidSource:setConfGrid(grid) self.grid = grid end
function FluidSource:setCfl(cfl) self.cfl = cfl end

function FluidSource:createSolver(mySpecies, externalField)

   self.writeGhost = mySpecies.writeGhost

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

   self.momSource:write(string.format("%s_0.bp", self.name), 0.0, 0, self.writeGhost)
end

function FluidSource:advance(tCurr, momIn, species, momRhsOut)
   local tm = Time.clock()

   Mpi.Barrier(self.grid:commSet().sharedComm)
   momRhsOut:accumulate(self.timeDependence(tCurr), self.momSource)

   self.timers.accumulateSrc = self.timers.accumulateSrc + Time.clock() - tm
end

function FluidSource:getSource() return self.momSource end
function FluidSource:srcTime() return self.timers.accumulateSrc end

return FluidSource
