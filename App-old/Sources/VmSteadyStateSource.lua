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
local lume             = require "Lib.lume"

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
   self.profile = Projection.VlasovProjection.FunctionProjection { func = function(t, zn) return tbl.profile(t, zn, self) end, }

   self.timers = {advance = 0.}
end

function VmSteadyStateSource:setName(nm) self.name = self.speciesName.."_"..nm end
function VmSteadyStateSource:setSpeciesName(nm) self.speciesName = nm end
function VmSteadyStateSource:setConfBasis(basis) self.confBasis = basis end
function VmSteadyStateSource:setConfGrid(grid) self.confGrid = grid end

function VmSteadyStateSource:createSolver(mySpecies, extField)
   
   self.writeGhost = mySpecies.writeGhost
   
   self.profile:fullInit(mySpecies)
   self.fSource = mySpecies:allocDistf()
   self.profile:createSolver(mySpecies)
   self.profile:advance(0.0, {mySpecies, extField}, {self.fSource})

   if self.positivityRescale then
      mySpecies.posRescaler:advance(0.0, {self.fSource}, {self.fSource}, false)
   end

   if self.power then
      local calcInt = Updater.CartFieldIntegrate {
         onGrid = self.confGrid,  basis = self.confBasis,
      }
      local intKE = DataStruct.DynVector{numComponents = 1}
      mySpecies.ptclEnergyCalc:advance(0.0, {self.fSource}, {mySpecies.ptclEnergyAux})
      calcInt:advance(0.0, {mySpecies.ptclEnergyAux, mySpecies.mass/2}, {intKE})
      local _, intKE_data  = intKE:lastData()
      self.powerScalingFac = self.power/intKE_data[1]
      self.fSource:scale(self.powerScalingFac)
   end

   local momsSrc = mySpecies:allocVectorMoment(mySpecies.vdim+2)
   local fiveMomentsCalc = Updater.DistFuncMomentCalc {
      onGrid     = mySpecies.grid,   confBasis  = mySpecies.confBasis,
      phaseBasis = mySpecies.basis,  moment     = "FiveMoments",
   }
   fiveMomentsCalc:advance(0.0, {self.fSource}, {momsSrc})

   self.fSource:write(string.format("%s_0.bp", self.name), 0., 0, self.writeGhost)
   momsSrc:write(string.format("%s_Moms_0.bp", self.name), 0., 0)

   -- Need to define methods to allocate fields (used by diagnostics).
   self.allocMoment = function() return mySpecies:allocMoment() end
end

function VmSteadyStateSource:initCrossSpeciesCoupling(population)
   local species = population:getSpecies()

   self.srcBC = {}
   local hasNonPeriodic = false 
   for _, otherNm in ipairs(self.sourceSpecies) do
      for _, bc in lume.orderedIter(species[otherNm].nonPeriodicBCs) do
	 self.srcBC[bc:getEdge()] = bc
         self.srcBC[bc:getEdge()]:setSaveFlux(true)
         hasNonPeriodic = true
      end
   end
   assert(hasNonPeriodic, "App.Sources.VmSteadyStateSource: has to have a non-periodic BC.")
   lume.setOrder(self.srcBC)
end

function VmSteadyStateSource:createCouplingSolver(population, field, extField)
   self.localEdgeFlux = Lin.Vec(1)
   self.globalEdgeFlux = Lin.Vec(1)
   
   self.bcSrcFlux = {}
   for _, bc in lume.orderedIter(self.srcBC) do
      self.bcSrcFlux[bc:getEdge()] = bc:allocIntThreeMoments()
   end
end

function VmSteadyStateSource:advanceCrossSpeciesCoupling(tCurr, species, outIdx)
   local tmStart = Time.clock()

   local mySpecies = species[self.speciesName]

   local fRhsOut = mySpecies:rkStepperFields()[outIdx]
   self.localEdgeFlux:data()[0] = 0.0

   for _, bc in lume.orderedIter(self.srcBC) do
      local confRange = self.bcSrcFlux[bc:getEdge()]:localExtRange()
      local phaseRange = bc:getBoundaryFluxFields()[outIdx]:localExtRange()
      bc.integNumDensityCalc:advance(tCurr, {bc:getBoundaryFluxFields()[outIdx], confRange, phaseRange}, {self.bcSrcFlux[bc:getEdge()]})

      local flux = self.bcSrcFlux[bc:getEdge()]
      
      if GKYL_USE_GPU then flux:copyDeviceToHost() end

      local fluxIndexer, fluxItr = flux:genIndexer(), flux:get(0)
      for idx in flux:localExtRangeIter() do
         flux:fill(fluxIndexer(idx), fluxItr)
         self.localEdgeFlux:data()[0] = self.localEdgeFlux:data()[0] + math.abs(fluxItr[1])
      end
      Mpi.Allreduce(self.localEdgeFlux:data(), self.globalEdgeFlux:data(), 1, Mpi.DOUBLE, Mpi.SUM, self.confGrid:commSet().host)
   end
   local densFactor = self.globalEdgeFlux:data()[0]/self.sourceLength

   fRhsOut:accumulate(densFactor, self.fSource)
   
   self.timers.advance = self.timers.advance + Time.clock() - tmStart
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
