-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Vlasov source operator
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

local VmSource = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmSource:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmSource:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.
   if tbl.timeDependence then
      self.timeDependence = tbl.timeDependence
   else
      self.timeDependence = function (t) return 1.0 end
   end
   self.power = tbl.power
   if tbl.profile then
      self.profile = tbl.profile
   elseif tbl.type then
      self.density = assert(tbl.density, "App.VmSource: must specify density profile of source in 'density'.")
      self.temperature = assert(tbl.temperature, "App.VmSource: must specify temperature profile of source in 'density'.")
      if tbl.type == "Maxwellian" or tbl.type == "maxwellian" then
         self.profile = Projection.VlasovProjection.MaxwellianProjection {
            density = self.density,
            temperature = self.temperature,
            power = self.power,
         }
      else
         print("App.VmSource: Source type not recognized, defaulting to Maxwellian.")
	 self.profile = Projection.VlasovProjection.MaxwellianProjection {
            density = self.density,
            temperature = self.temperature,
            power = self.power,
         }
      end    
   else
      self.density = assert(tbl.density, "App.VmSource: must specify density profile of source in 'density'.")
      self.temperature = assert(tbl.temperature, "App.VmSource: must specify temperature profile of source in 'density'.")
      self.profile = Projection.VlasovProjection.MaxwellianProjection {
         density = self.density,
         temperature = self.temperature,
         power = self.power,
      }
   end
   self.tmEvalSrc = 0.0
end

function VmSource:setName(nm)
   self.name = nm
end
function VmSource:setSpeciesName(nm)
   self.speciesName = nm
end
function VmSource:setConfBasis(basis)
   self.confBasis = basis
end
function VmSource:setConfGrid(grid)
   self.confGrid = grid
end

function VmSource:advance(tCurr, fIn, species, fRhsOut)
   local tm = Time.clock()
   Mpi.Barrier(self.confGrid:commSet().sharedComm)
   fRhsOut:accumulate(self.timeDependence(tCurr), species[self.speciesName].fSource)
   self.tmEvalSrc = self.tmEvalSrc + Time.clock() - tm
end

function VmSource:createDiagnostics(thisSpecies, momTable)
   self.diagnosticIntegratedMomentFields   = { }
   self.diagnosticIntegratedMomentUpdaters = { }
   self.diagnosticIntegratedMoments = { }
   
   self.numDensitySrc = thisSpecies:allocMoment()
   self.momDensitySrc = thisSpecies:allocVectorMoment(thisSpecies.vdim)
   self.ptclEnergySrc = thisSpecies:allocMoment()
   thisSpecies.fiveMomentsCalc:advance(0.0, {thisSpecies.fSource}, {self.numDensitySrc, self.momDensitySrc, self.ptclEnergySrc})
end

function VmSource:writeDiagnosticIntegratedMoments(tm, frame)
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
       self.diagnosticIntegratedMomentFields[mom]:write(string.format("%s_%s.bp", self.speciesName, mom), tm, frame)
    end
end

function VmSource:write(tm, frame, thisSpecies)
   thisSpecies.fSource:write(string.format("%s_fSource_0.bp", self.speciesName), tm, frame, true)
   if self.numDensitySrc then self.numDensitySrc:write(string.format("%s_srcM0_0.bp", self.speciesName), tm, frame) end
   if self.momDensitySrc then self.momDensitySrc:write(string.format("%s_srcM1_0.bp", self.speciesName), tm, frame) end
   if self.ptclEnergySrc then self.ptclEnergySrc:write(string.format("%s_srcM2_0.bp", self.speciesName), tm, frame) end
end

function VmSource:srcTime()
   return self.tmEvalSrc
end

return VmSource
