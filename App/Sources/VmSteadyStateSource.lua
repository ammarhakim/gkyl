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

   self.writeGhost = mySpecies.writeGhost   

   self.profile:fullInit(mySpecies)

   self.fSource = mySpecies:allocDistf()

   self.profile:advance(0.0, {extField}, {self.fSource})

   if self.positivityRescale then
      mySpecies.posRescaler:advance(0.0, {self.fSource}, {self.fSource}, false)
   end

   if self.power then
      local calcInt = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.confGrid,
         basis = self.confBasis,
         numComponents = 1,
         quantity = "V",
      }
      local intKE = DataStruct.DynVector{numComponents = 1}
      mySpecies.ptclEnergyCalc:advance(0.0, {self.fSource}, {mySpecies.ptclEnergyAux})
      calcInt:advance(0.0, {mySpecies.ptclEnergyAux, mySpecies.mass/2}, {intKE})
      local _, intKE_data = intKE:lastData()
      self.powerScalingFac = self.power/intKE_data[1]
      self.fSource:scale(self.powerScalingFac)
   end

   local numDensitySrc = mySpecies:allocMoment()
   local momDensitySrc = mySpecies:allocVectorMoment(mySpecies.vdim)
   local ptclEnergySrc = mySpecies:allocMoment()
   mySpecies.fiveMomentsCalc:advance(0.0, {self.fSource}, {numDensitySrc, momDensitySrc, ptclEnergySrc})

   self.fSource:write(string.format("%s_0.bp", self.name), 0., 0, self.writeGhost)
   numDensitySrc:write(string.format("%s_M0_0.bp", self.name), 0., 0)
   momDensitySrc:write(string.format("%s_M1_0.bp", self.name), 0., 0)
   ptclEnergySrc:write(string.format("%s_M2_0.bp", self.name), 0., 0)

   self.momDensityBuf = mySpecies:allocVectorMoment(mySpecies.vdim)

   local numConfDims  = self.confBasis:ndim()
   local numConfBasis = self.confBasis:numBasis()
   assert(numConfDims==1, "VlasovSpecies: The steady state source is available only for 1X.")

   local blowerV, bupperV = Lin.Vec(numConfDims), Lin.Vec(numConfDims)
   blowerV[1], bupperV[1] = -1.0, 1.0
   self.basisUpper, self.basisLower = Lin.Vec(numConfBasis),  Lin.Vec(numConfBasis)
   self.confBasis:evalBasis(bupperV:data(), self.basisUpper:data())
   self.confBasis:evalBasis(blowerV:data(), self.basisLower:data())

   self.localEdgeFlux = ffi.new("double[1]")
   self.globalEdgeFlux = ffi.new("double[1]")
end

function VmSteadyStateSource:advance(tCurr, fIn, species, fRhsOut)
   local tm = Time.clock()
   self.localEdgeFlux[0] = 0.0
   local numConfBasis = self.confBasis:numBasis()
   for _, otherNm in ipairs(self.sourceSpecies) do

      local moms = species[otherNm]:fluidMoments()
      local flux = self.momDensityBuf
      flux:combineOffset(1., moms, 1*numConfBasis)

      if GKYL_USE_GPU then flux:copyDeviceToHost() end

      local fluxIndexer, fluxItr = flux:genIndexer(), flux:get(1)
      for idx in flux:localRangeIter() do
         if idx[1] == self.confGrid:numCells(1) then
            flux:fill(fluxIndexer(idx), fluxItr)
            for k = 1, numConfBasis do
               self.localEdgeFlux[0] = self.localEdgeFlux[0] + fluxItr[k]*self.basisUpper[k]
            end
         elseif idx[1] == 1 then
            flux:fill(fluxIndexer(idx), fluxItr)
            for k = 1, numConfBasis do
               self.localEdgeFlux[0] = self.localEdgeFlux[0] - fluxItr[k]*self.basisLower[k]
            end
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
