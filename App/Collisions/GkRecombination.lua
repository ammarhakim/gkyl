-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Recombination operator ADAS data
--
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local Constants      = require "Lib.Constants"
local DataStruct     = require "DataStruct"
local DiagsImplBase  = require "App.Diagnostics.DiagnosticsImplBase"
local DiagsApp       = require "App.Diagnostics.SpeciesDiagnostics"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local Mpi            = require "Comm.Mpi"
local lume           = require "Lib.lume"
local xsys           = require "xsys"

-- GkRecombination -------------------------------------------------------------
--
-- Recombination operator.
--------------------------------------------------------------------------------

-- ............... IMPLEMENTATION OF DIAGNOSTICS ................. --
-- Diagnostics could be placed in a separate file if they balloon in
-- number. But if we only have one or two we can just place it here.

local gkRecombDiagImpl = function()
   local _M0 = Proto(DiagsImplBase)
   function _M0:fullInit(diagApp, mySpecies, fieldIn, owner)
      self.field    = mySpecies:allocMoment()
      self.updater  = mySpecies.numDensityCalc
      self.owner    = owner
      self.done     = false
   end
   function _M0:getType() return "grid" end
   function _M0:advance(tm, inFlds, outFlds)
      self.updater:advance(tm, {self.owner.recombSrc}, {self.field})
   end
   
   local _intM0 = Proto(DiagsImplBase)
   function _intM0:fullInit(diagApp, mySpecies, fieldIn, owner)
      self.fieldAux = mySpecies:allocMoment()
      self.updatersAux = mySpecies.numDensityCalc
      self.field    = DataStruct.DynVector { numComponents = 1 }
      self.updater  = mySpecies.volIntegral.scalar
      self.owner    = owner
      self.done     = false
   end
   function _intM0:getType() return "integrated" end
   function _intM0:advance(tm, inFlds, outFlds)
      self.updatersAux:advance(tm, {self.owner.recombSrc}, {self.fieldAux})
      self.updater:advance(tm, {self.fieldAux}, {self.field})
   end

   local _reactRate = Proto(DiagsImplBase)
   function _reactRate:fullInit(diagApp, mySpecies, fieldIn, owner)
      self.field = mySpecies:allocMoment()
      self.owner = owner
      self.done  = false
   end
   function _reactRate:getType() return "grid" end
   function _reactRate:advance(tm, inFlds, outFlds)
      if self.owner.reactRate then
	 self.field:copy(self.owner.reactRate)
      end
   end

   local _source = Proto(DiagsImplBase)
   function _source:fullInit(diagApp, mySpecies, fieldIn, owner)
      self.field = mySpecies:allocDistf()
      self.owner = owner
      self.done  = false
   end
   function _source:getType() return "grid" end
   function _source:advance(tm, inFlds, outFlds)
      self.field:copy(self.owner.recombSrc)
   end

   return {
      M0        = _M0,
      intM0     = _intM0,
      reactRate = _reactRate,
      source    = _source
   }
end

-- .................... END OF DIAGNOSTICS ...................... --

local GkRecombination = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkRecombination:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkRecombination:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously store table.

   self.cfl = 0.1
   self.collKind = "Recombination"

   self.collidingSpecies = assert(tbl.collideWith, "App.GkRecombination: Must specify names of species to collide with in 'collideWith'.")

   -- Set these values to be consistent with other collision apps
   self.selfCollisions  = false
   self.crossCollisions = true              
   self.varNu           = false
   self.timeDepNu       = false
   self.collFreqs       = {1}
   
   self.collideNm = tbl.collideWith[1]

   -- Intake for these values needs to be optimized
   self.elcNm  = assert(tbl.elcName, "App.GkRecombination: Must specify electron species name in 'elcName'.")
   self.ionNm = assert(tbl.ionName, "App.GkRecombination: Must specify donor species name in 'donorName'.")
   self.selfMass = tbl.selfMass
   self.plasma = tbl.plasma
   self.chargeState = tbl.chargeState
   self.selfSpecies = tbl.selfSpecies
   self.recvrGk = tbl.recvrGk
   
   self.timers = {nonSlvr = 0.}
end

function GkRecombination:createDiagnostics(mySpecies, field)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = gkRecombDiagImpl()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

function GkRecombination:setName(nm)
   self.name = self.speciesName.."_"..nm
   self.collNm = nm
end
function GkRecombination:setSpeciesName(nm) self.speciesName = nm end
function GkRecombination:setCfl(cfl) self.cfl = cfl end
function GkRecombination:setConfBasis(basis) self.confBasis = basis end
function GkRecombination:setConfGrid(grid) self.confGrid = grid end
function GkRecombination:setPhaseBasis(basis) self.phaseBasis = basis end
function GkRecombination:setPhaseGrid(grid) self.phaseGrid = grid end

function GkRecombination:createSolver(mySpecies, externalField)
   -- Geometry fields
   self.bmag = externalField.geo.bmag
   self.jacobTot = externalField.geo.jacobTot
   self.b_i = externalField.geo.b_i

   -- Reaction rate to be used for diagnostics
   self.reactRate = mySpecies:allocMoment()

   -- Recombination source term
   self.recombSrc = mySpecies:allocDistf()
   
   -- Get conf and phase range
   self.confRange = self.bmag:localRange()
   self.phaseRange = self.recombSrc:localRange()
   
   self.collisionSlvr = Updater.Recombination {
      onGrid = self.phaseGrid,        recvrGk = self.recvrGk,
      confBasis = self.confBasis,     phaseBasis = self.phaseBasis,
      confRange = self.confRange,     phaseRange = self.phaseRange,
      selfMass = self.selfMass,       plasma = self.plasma,
      chargeState = self.chargeState, selfSpecies = self.selfSpecies,
   }
end

function GkRecombination:advance(tCurr, fIn, population, out)
   local tmNonSlvrStart = Time.clock()
   local species = population:getSpecies()
   local momsElc = species[self.elcNm]:fluidMoments()
   local momsIon = species[self.ionNm]:fluidMoments()

   -- species[self.speciesName].distIo:write(momsDonor, string.format("%s_momsDonor_%d.bp", self.speciesName, species[self.speciesName].diagIoFrame), 0., species[self.speciesName].diagIoFrame, false) --true)
   local fRhsOut = out[1]
   local cflRateByCell = out[2]
   
   self.recombSrc:clear(0.0)
   self.collisionSlvr:advance(tCurr, {momsElc, momsIon, self.bmag, self.jacobTot, self.b_i, fIn},
   			      {self.recombSrc, cflRateByCell})

   fRhsOut:accumulate(1.0,self.recombSrc)

   self.timers.nonSlvr = self.timers.nonSlvr + Time.clock() - tmNonSlvrStart
end

function GkRecombination:write(tm, frame) end

function GkRecombination:setCfl(cfl)
   self.cfl = cfl
end

function GkRecombination:slvrTime()
   return 0.
end

function GkRecombination:nonSlvrTime()
   return self.timers.nonSlvr
end



return GkRecombination

