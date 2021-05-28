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

local VmSteadyStateSource = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmSteadyStateSource:init(tbl)
   self.tbl = tbl
end

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

function VmSteadyStateSource:setName(nm)
   self.name = nm
end
function VmSteadyStateSource:setSpeciesName(nm)
   self.speciesName = nm
end
function VmSteadyStateSource:setConfBasis(basis)
   self.confBasis = basis
end
function VmSteadyStateSource:setConfGrid(grid)
   self.confGrid = grid
end

function VmSteadyStateSource:createSolver(mySpecies, extField)
   self.profile:fullInit(mySpecies)

   self.fSource = mySpecies:allocDistf()

   self.profile:advance(0.0, {extField}, {self.fSource})
   Mpi.Barrier(mySpecies.grid:commSet().sharedComm)

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
end

function VmSteadyStateSource:advance(tCurr, fIn, species, fRhsOut)
   local tm = Time.clock()
   local localEdgeFlux = ffi.new("double[3]")
   localEdgeFlux[0] = 0.0
   localEdgeFlux[1] = 0.0
   localEdgeFlux[2] = 0.0
   for sInd, otherNm in ipairs(self.sourceSpecies) do
      local numConfDims = self.confBasis:ndim()
      assert(numConfDims==1, "VlasovSpecies: The steady state source is available only for 1X.")
      local numConfBasis = self.confBasis:numBasis()
      local lower, upper = Lin.Vec(numConfDims), Lin.Vec(numConfDims)
      lower[1], upper[1] = -1.0, 1.0
      local basisUpper = Lin.Vec(numConfBasis)
      local basisLower = Lin.Vec(numConfBasis)

      self.confBasis:evalBasis(upper, basisUpper)
      self.confBasis:evalBasis(lower, basisLower)

      local flux = species[otherNm]:fluidMoments()[2]
      local fluxIndexer, fluxItr = flux:genIndexer(), flux:get(1)
      for idx in flux:localRangeIter() do
         if idx[1] == self.confGrid:numCells(1) then
            flux:fill(fluxIndexer(idx), fluxItr)
            for k = 1, numConfBasis do
               localEdgeFlux[0] = localEdgeFlux[0] + fluxItr[k]*basisUpper[k]
            end
         elseif idx[1] == 1 then
            flux:fill(fluxIndexer(idx), fluxItr)
            for k = 1, numConfBasis do
               localEdgeFlux[0] = localEdgeFlux[0] - fluxItr[k]*basisLower[k]
            end
         end
      end
   end
   local globalEdgeFlux = ffi.new("double[3]")
   Mpi.Allreduce(localEdgeFlux, globalEdgeFlux, 1, Mpi.DOUBLE, Mpi.MAX, self.confGrid:commSet().comm)
   local densFactor = globalEdgeFlux[0]/self.sourceLength
   fRhsOut:accumulate(densFactor, self.fSource)
   self.tmEvalSrc = self.tmEvalSrc + Time.clock() - tm
end

function VmSteadyStateSource:createDiagnostics(thisSpecies, momTable)
   self.diagnosticIntegratedMomentFields   = { }
   self.diagnosticIntegratedMomentUpdaters = { }
   self.diagnosticIntegratedMoments = { }
   
   self.numDensitySrc = thisSpecies:allocMoment()
   self.momDensitySrc = thisSpecies:allocVectorMoment(thisSpecies.vdim)
   self.ptclEnergySrc = thisSpecies:allocMoment()
   thisSpecies.fiveMomentsCalc:advance(0.0, {self.fSource}, {self.numDensitySrc, self.momDensitySrc, self.ptclEnergySrc})
end

function VmSteadyStateSource:writeDiagnosticIntegratedMoments(tm, frame)
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
       self.diagnosticIntegratedMomentFields[mom]:write(string.format("%s_%s.bp", self.speciesName, mom), tm, frame)
    end
end


function VmSteadyStateSource:write(tm, frame)
   if tm == 0.0 then
      self.fSource:write(string.format("%s_fSource_0.bp", self.speciesName), tm, frame, true)
      if self.numDensitySrc then self.numDensitySrc:write(string.format("%s_srcM0_0.bp", self.speciesName), tm, frame) end
      if self.momDensitySrc then self.momDensitySrc:write(string.format("%s_srcM1_0.bp", self.speciesName), tm, frame) end
      if self.ptclEnergySrc then self.ptclEnergySrc:write(string.format("%s_srcM2_0.bp", self.speciesName), tm, frame) end
   end
end

function VmSteadyStateSource:srcTime()
   return self.tmEvalSrc
end

return VmSteadyStateSource
