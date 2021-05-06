-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Gyrokinetic source operator.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase     = require "App.Sources.SourceBase"
local DataStruct     = require "DataStruct"
local lume           = require "Lib.lume"
local Mpi            = require "Comm.Mpi"
local Projection     = require "App.Projection"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"

local GkSource = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkSource:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkSource:fullInit(thisSpecies)
   local tbl = self.tbl -- Previously stored table.

   if tbl.timeDependence then
      self.timeDependence = tbl.timeDependence
   else
      self.timeDependence = function (t) return 1.0 end
   end

   self.power = tbl.power

   if tbl.profile then
      self.profile = tbl.profile
   elseif tbl.kind then
      self.density     = assert(tbl.density, "App.GkSource: must specify density profile of source in 'density'.")
      self.temperature = assert(tbl.temperature, "App.GkSource: must specify temperature profile of source in 'density'.")
      if tbl.kind == "Maxwellian" or tbl.kind == "maxwellian" then
         self.profile   = Projection.GkProjection.MaxwellianProjection {
            density     = self.density,
            temperature = self.temperature,
            power       = self.power,
         }
      else
         assert(false, "App.GkSource: Source kind not recognized.")
      end    
   else
      -- If user doesn't specify 'kind', default to Maxwellian.
      self.density     = assert(tbl.density, "App.GkSource: must specify density profile of source in 'density'.")
      self.temperature = assert(tbl.temperature, "App.GkSource: must specify temperature profile of source in 'density'.")
      self.profile     = Projection.GkProjection.MaxwellianProjection {
         density     = self.density,
         temperature = self.temperature,
         power       = self.power,
      }
   end
   self.tmEvalSrc = 0.0
end

function GkSource:setName(nm) self.name = nm end
function GkSource:setSpeciesName(nm) self.speciesName = nm end
function GkSource:setConfBasis(basis) self.confBasis = basis end
function GkSource:setConfGrid(grid) self.confGrid = grid end

function GkSource:createSolver(thisSpecies, extField)
   self.profile:fullInit(thisSpecies)
   self.profile:advance(0.0, {extField}, {thisSpecies.distf[2]})
   Mpi.Barrier(thisSpecies.grid:commSet().sharedComm)

   if not self.fSource then self.fSource = thisSpecies:allocDistf() end
   self.fSource:accumulate(1.0, thisSpecies.distf[2])

   if self.positivityRescale then
      thisSpecies.posRescaler:advance(0.0, {self.fSource}, {self.fSource}, false)
   end

   if self.power then
      local calcInt = Updater.CartFieldIntegratedQuantCalc {
         onGrid        = self.confGrid,
         basis         = self.confBasis,
         numComponents = 1,
         quantity      = "V",
      }
      local intKE = DataStruct.DynVector{numComponents = 1}
      thisSpecies.ptclEnergyCalc:advance(0.0, {self.fSource}, {thisSpecies.ptclEnergyAux})
      calcInt:advance(0.0, {thisSpecies.ptclEnergyAux, thisSpecies.mass/2}, {intKE})
      local _, intKE_data  = intKE:lastData()
      self.powerScalingFac = self.power/intKE_data[1]
      self.fSource:scale(self.powerScalingFac)
   end
   if thisSpecies.scaleInitWithSourcePower then thisSpecies.distf[1]:scale(self.powerScalingFac) end
end

function GkSource:advance(tCurr, fIn, species, fRhsOut)
   local tm = Time.clock()
   Mpi.Barrier(self.confGrid:commSet().sharedComm)
   fRhsOut:accumulate(self.timeDependence(tCurr), self.fSource)
   self.tmEvalSrc = self.tmEvalSrc + Time.clock() - tm
end

function GkSource:createDiagnostics(thisSpecies, momTable)
   local function contains(table, element) return lume.any(table, function(e) return e==element end) end

   self.diagnosticIntegratedMomentFields   = { }
   self.diagnosticIntegratedMomentUpdaters = { }
   self.diagnosticIntegratedMoments        = { }
   self.diagnosticIntegratedMomentTimeDep  = { }
   
   self.numDensitySrc = thisSpecies:allocMoment()
   self.momDensitySrc = thisSpecies:allocMoment()
   self.ptclEnergySrc = thisSpecies:allocMoment()
   thisSpecies.threeMomentsCalc:advance(0.0, {self.fSource}, {self.numDensitySrc, self.momDensitySrc, self.ptclEnergySrc})
   
   if contains(momTable, "intM0") or contains(momTable, "intSrcM0") then
      table.insert(self.diagnosticIntegratedMoments, "intSrcM0")
      self.diagnosticIntegratedMomentFields["intSrcM0"] = {}
      self.diagnosticIntegratedMomentFields["intSrcM0"].inFld = self.numDensitySrc
      self.diagnosticIntegratedMomentTimeDep["intSrcM0"] = function(t) return self.timeDependence(t) end
   end
   if contains(momTable, "intM1") or contains(momTable, "intSrcM1") then
      table.insert(self.diagnosticIntegratedMoments, "intSrcM1")
      self.diagnosticIntegratedMomentFields["intSrcM1"] = {}
      self.diagnosticIntegratedMomentFields["intSrcM1"].inFld = self.momDensitySrc
      self.diagnosticIntegratedMomentTimeDep["intSrcM1"] = function(t) return self.timeDependence(t) end
   end
   if contains(momTable, "intM2") or contains(momTable, "intSrcM2") then
      table.insert(self.diagnosticIntegratedMoments, "intSrcM2")
      self.diagnosticIntegratedMomentFields["intSrcM2"] = {}
      self.diagnosticIntegratedMomentFields["intSrcM2"].inFld = self.ptclEnergySrc
      self.diagnosticIntegratedMomentTimeDep["intSrcM2"] = function(t) return self.timeDependence(t) end
   end
   if contains(momTable, "intKE") or contains(momTable, "intSrcKE") then
      table.insert(self.diagnosticIntegratedMoments, "intSrcKE")
      self.diagnosticIntegratedMomentFields["intSrcKE"] = {}
      self.diagnosticIntegratedMomentFields["intSrcKE"].inFld = self.ptclEnergySrc
      self.diagnosticIntegratedMomentTimeDep["intSrcKE"] = function(t) return self.timeDependence(t)*thisSpecies.mass/2. end
   end
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
      self.diagnosticIntegratedMomentFields[mom].outFld = DataStruct.DynVector{numComponents = 1}
      self.diagnosticIntegratedMomentUpdaters[mom]      = Updater.CartFieldIntegratedQuantCalc {
         onGrid        = thisSpecies.confGrid,
         basis         = thisSpecies.confBasis,
         numComponents = 1,
         quantity      = "V",
         timeIntegrate = true,
      }
   end
end

function GkSource:calcDiagnosticIntegratedMoments(tm, thisSpecies)
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
      self.diagnosticIntegratedMomentUpdaters[mom]:advance(tm,
         {self.diagnosticIntegratedMomentFields[mom].inFld, self.diagnosticIntegratedMomentTimeDep[mom](tm)},
         {self.diagnosticIntegratedMomentFields[mom].outFld})
   end
end

function GkSource:writeDiagnosticIntegratedMoments(tm, frame)
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
      self.diagnosticIntegratedMomentFields[mom].outFld:write(string.format("%s_%s.bp", self.speciesName, mom), tm, frame)
   end
end

function GkSource:write(tm, frame)
   if tm == 0.0 then
      self.fSource:write(string.format("%s_fSource_0.bp", self.speciesName), tm, frame, true)
      self.numDensitySrc:write(string.format("%s_srcM0_0.bp", self.speciesName), tm, frame)
      self.momDensitySrc:write(string.format("%s_srcM1_0.bp", self.speciesName), tm, frame)
      self.ptclEnergySrc:write(string.format("%s_srcM2_0.bp", self.speciesName), tm, frame)
   end
end

function GkSource:srcTime() return self.tmEvalSrc end

return GkSource
