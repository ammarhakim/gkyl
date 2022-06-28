-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Gyrokinetic source operator.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local DataStruct = require "DataStruct"
local lume       = require "Lib.lume"
local Mpi        = require "Comm.Mpi"
local Projection = require "App.Projection"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"
local Updater    = require "Updater"
local DiagsApp   = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags    = require "App.Diagnostics.GkDiagnostics"

local GkSource = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkSource:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkSource:fullInit(thisSpecies)
   local tbl = self.tbl -- Previously stored table.

   self.timeDependence = tbl.timeDependence or function (t) return 1. end

   self.power = tbl.power

   if tbl.profile then
      if type(tbl.profile) == "function" then
	 self.profile = Projection.GkProjection.FunctionProjection {
	    func = function(t, zn) return tbl.profile(t, zn) end
	 }
      elseif type(tbl.profile) == "string" then
         self.profile = Projection.GkProjection.FunctionProjection{fromFile = tbl.profile,}
	 -- self.profile = Projection.ReadInput {
	 --    inputFile = tbl.profile,
	 -- }
      end
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

function GkSource:setName(nm) self.name = self.speciesName.."_"..nm end
function GkSource:setSpeciesName(nm) self.speciesName = nm end
function GkSource:setConfBasis(basis) self.confBasis = basis end
function GkSource:setConfGrid(grid) self.confGrid = grid end

function GkSource:createSolver(mySpecies, extField)

   self.writeGhost = mySpecies.writeGhost   

   self.profile:fullInit(mySpecies)

   self.fSource = mySpecies:allocDistf()

   self.profile:advance(0.0, {extField}, {self.fSource})
   Mpi.Barrier(mySpecies.grid:commSet().sharedComm)

   if self.positivityRescale then
      mySpecies.posRescaler:advance(0.0, {self.fSource}, {self.fSource}, false)
   end

   if self.power then
      local calcInt = Updater.CartFieldIntegratedQuantCalc {
         onGrid        = self.confGrid,
         basis         = self.confBasis,
         numComponents = 1,
         quantity      = "V",
      }
      local intKE = DataStruct.DynVector{numComponents = 1}
      mySpecies.ptclEnergyCalc:advance(0.0, {self.fSource}, {mySpecies.ptclEnergyAux})
      calcInt:advance(0.0, {mySpecies.ptclEnergyAux, mySpecies.mass/2}, {intKE})
      local _, intKE_data  = intKE:lastData()
      self.powerScalingFac = self.power/intKE_data[1]
      self.fSource:scale(self.powerScalingFac)
   end

   local momsSrc = mySpecies:allocVectorMoment(3)
   mySpecies.threeMomentsCalc:advance(0.0, {self.fSource}, {momsSrc})

   self.fSource:write(string.format("%s_0.bp", self.name), 0., 0, self.writeGhost)
   momsSrc:write(string.format("%s_Moms_0.bp", self.name), 0., 0)

   -- Need to define methods to allocate fields (used by diagnostics).
   self.allocMoment = function() return mySpecies:allocMoment() end
end

function GkSource:createDiagnostics(mySpecies, field)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = GkDiags()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

function GkSource:advance(tCurr, fIn, species, fRhsOut)
   local tm = Time.clock()
   Mpi.Barrier(self.confGrid:commSet().sharedComm)
   fRhsOut:accumulate(self.timeDependence(tCurr), self.fSource)
   self.tmEvalSrc = self.tmEvalSrc + Time.clock() - tm
end

-- These are needed to recycle the GkDiagnostics with GkSource.
function GkSource:rkStepperFields() return {self.fSource, self.fSource, self.fSource, self.fSource} end
function GkSource:getFlucF() return self.fSource end

function GkSource:write(tm, frame) end

function GkSource:srcTime() return self.tmEvalSrc end

return GkSource
