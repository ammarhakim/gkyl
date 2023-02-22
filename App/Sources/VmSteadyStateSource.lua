
-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Vlasov steady state source operator
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
local DiagsApp       = require "App.Diagnostics.SpeciesDiagnostics"
local VlasovDiags    = require "App.Diagnostics.VlasovDiagnostics"

local VmSteadyStateSource = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmSteadyStateSource:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmSteadyStateSource:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.
   self.sourceSpecies = assert(tbl.sourceSpecies, "App.VmSteadyStateSource: must specify names of species to in 'sourceSpecies'.")
   self.sourceLength = assert(tbl.sourceLength, "App.VmSteadyStateSource: must specify names of species to in 'sourceLength'.")
   assert(tbl.profile, "App.VmSteadyStateSource: must specify source profile in 'profile'")
   assert(type(tbl.profile) == "function", "App.VmSteadyStateSource: 'profile' must be a function")
   self.profile = Projection.KineticProjection.FunctionProjection { func = function(t, zn) return tbl.profile(t, zn, self) end, }

   self.tmEvalSrc = 0.0
end

function VmSteadyStateSource:setName(nm) self.name = self.speciesName.."_"..nm end
function VmSteadyStateSource:setSpeciesName(nm) self.speciesName = nm end
function VmSteadyStateSource:setConfBasis(basis) self.confBasis = basis end
function VmSteadyStateSource:setConfGrid(grid) self.confGrid = grid end

function VmSteadyStateSource:createSolver(mySpecies, extField)

   self.dgIntegratedMoms = Updater.DgMomentCalc {
      onGrid     = mySpecies.grid,   confBasis  = mySpecies.confBasis,
      phaseBasis = mySpecies.basis,  moment     = "M0",
      isIntegrated = true,
   }

   self.writeGhost = mySpecies.writeGhost
   
   self.profile:fullInit(mySpecies)
   self.fSource = mySpecies:allocDistf()
   self.profile:advance(0.0, {extField}, {self.fSource})

   local distf = mySpecies:getDistF()

   self.integMoms  = mySpecies:allocIntMoment(mySpecies.vdim+2)

   self.localEdgeFlux = ffi.new("double[1]")
   self.globalEdgeFlux = ffi.new("double[1]")
end

function VmSteadyStateSource:advance(tCurr, fIn, species, fRhsOut)
   local tm = Time.clock()

   self.dgIntegratedMoms:advance(tCurr, {fRhsOut}, {self.integMoms})
   
   self.localEdgeFlux[0] = 0.0
   for _, otherNm in ipairs(self.sourceSpecies) do
      
      local flux = self.integMoms

      if GKYL_USE_GPU then flux:copyDeviceToHost() end

      local fluxIndexer, fluxItr = flux:genIndexer(), flux:get(0)
      for idx in flux:localExtRangeIter() do
         if idx[1] == self.confGrid:numCells(1) + 1 then
            flux:fill(fluxIndexer(idx), fluxItr)
            self.localEdgeFlux[0] = self.localEdgeFlux[0] + fluxItr[1]
         elseif idx[1] == 0 then
            flux:fill(fluxIndexer(idx), fluxItr)
            self.localEdgeFlux[0] = self.localEdgeFlux[0] + fluxItr[1]
         end
      end
   end
   Mpi.Allreduce(self.localEdgeFlux, self.globalEdgeFlux, 1, Mpi.DOUBLE, Mpi.SUM, self.confGrid:commSet().comm)
   local densFactor = self.globalEdgeFlux[0]/self.sourceLength

   fRhsOut:accumulate(densFactor, self.fSource)
   
   self.tmEvalSrc = self.tmEvalSrc + Time.clock() - tm
end

function VmSteadyStateSource:createDiagnostics(thisSpecies, momTable)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = VlasovDiags()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end

   return self.diagnostics
end

-- These are needed to recycle the VlasovDiagnostics with VmSource.
function VmSteadyStateSource:rkStepperFields() return {self.fSource, self.fSource, self.fSource, self.fSource} end
function VmSteadyStateSource:getFlucF() return self.fSource end

function VmSteadyStateSource:write(tm, frame) end

function VmSteadyStateSource:srcTime() return self.tmEvalSrc end

return VmSteadyStateSource
