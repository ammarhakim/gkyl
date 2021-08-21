-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Charge exchange operator 
-- For details of model see:
-- H. L. Pauls, G. P. Zank, and L. L. Williams. Interaction of the solar wind
-- with the local interstellar medium. J. Geophys. Research,
-- 100:21,595 – 21,604, 1995
--
-- Meier, E. T. Modeling Plasmas with Strong Anisotropy, Neutral Fluid Effects,
-- and Open Boundaries. PhD thesis, U. Washington (2011).
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

-- VmChargeExchange  --------------------------------------------------------
--
-- Charge Exchange Operator
-----------------------------------------------------------------------------

-- ............... IMPLEMENTATION OF DIAGNOSTICS ................. --
-- Diagnostics could be placed in a separate file if they balloon in
-- number. But if we only have one or two we can just place it here.

local vmCxDiagImpl = function()
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
      self.field:copy(self.owner.sourceCX)
   end

   return {
      reactRate = _reactRate,
      source    = _source,
   }
end

-- .................... END OF DIAGNOSTICS ...................... --


local VmChargeExchange = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmChargeExchange:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmChargeExchange:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously store table.

   self.collKind = "CX"

   self.collidingSpecies = assert(tbl.collideWith, "App.VmChargeExchange: Must specify names of species to collide with in 'collideWith'.")

   self.collideNm = tbl.collideWith[1]

   self.plasma    = assert(tbl.plasma, "App.VmChargeExchange: Must specify plasma species in 'plasma' ('H', 'D', or 'Ne')")
   
   self.ionNm     = tbl.ions
   self.neutNm    = tbl.neutrals
   self.iMass     = tbl.ionMass
   self.nMass     = tbl.neutMass
   self.charge    = tbl.charge

   -- Set these values to be consistent with other collision apps
   self.selfCollisions  = false
   self.crossCollisions = true              
   self.varNu           = false
   self.timeDepNu       = false
   self.collFreqs       = {1}

   if self.plasma=='H' then
      self.a = 1.12e-18
      self.b = 7.15e-20
   elseif self.plasma=='D' then 
      self.a = 1.09e-18
      self.b = 7.15e-20
   elseif self.plasma =='Ne' then
      self.a = 7.95e-19
      self.b = 5.65e-20
   end

   self.timers = {nonSlvr = 0.}
end

function VmChargeExchange:createDiagnostics(mySpecies, field)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = vmCxDiagImpl()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

function VmChargeExchange:setName(nm)
   self.name = self.speciesName.."_"..nm
   self.collNm = nm
end
function VmChargeExchange:setSpeciesName(nm) self.speciesName = nm end
function VmChargeExchange:setCfl(cfl) self.cfl = cfl end
function VmChargeExchange:setConfBasis(basis) self.confBasis = basis end
function VmChargeExchange:setConfGrid(grid) self.confGrid = grid end
function VmChargeExchange:setPhaseBasis(basis) self.phaseBasis = basis end
function VmChargeExchange:setPhaseGrid(grid) self.phaseGrid = grid end

function VmChargeExchange:createSolver(funcField) --species)
   self.collisionSlvr = Updater.ChargeExchange {
      onGrid     = self.confGrid,
      confBasis  = self.confBasis,
      phaseBasis = self.phaseBasis,
      charge     = self.charge,
      a          = self.a,
      b          = self.b,
   }
   self.sourceCX = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData      = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.M0iDistFn = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData      = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.M0nDistFi = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData      = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.diffDistF = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData      = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   if (self.speciesName == self.ionNm) then
      self.reactRate =  DataStruct.Field {
	 onGrid        = self.confGrid,
	 numComponents = self.confBasis:numBasis(),
	 ghost         = {1, 1},
	 metaData = {
	    polyOrder = self.confBasis:polyOrder(),
	    basisType = self.confBasis:id()
	 },
      }
   end
end

function VmChargeExchange:advance(tCurr, fIn, species, fRhsOut)
   local tmNonSlvrStart = Time.clock()

   local neutM0    = species[self.neutNm]:fluidMoments()[1]
   local neutDistF = species[self.neutNm]:getDistF()
   local ionM0     = species[self.ionNm]:fluidMoments()[1]
   local ionDistF  = species[self.ionNm]:getDistF()
   
   species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {ionM0, neutDistF}, {self.M0iDistFn})
   species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {neutM0, ionDistF}, {self.M0nDistFi})
   self.diffDistF:combine(1.0, self.M0iDistFn, -1.0, self.M0nDistFi)
   species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {species[self.ionNm].collisions[self.collNm].reactRate, self.diffDistF}, {self.sourceCX})

   if (self.speciesName == self.ionNm) then
      fRhsOut:accumulate(1.0,self.sourceCX)
   else
      fRhsOut:accumulate(-self.iMass/self.nMass,self.sourceCX)
   end

   self.timers.nonSlvr = self.timers.nonSlvr + Time.clock() - tmNonSlvrStart
end

function VmChargeExchange:write(tm, frame) end

function VmChargeExchange:slvrTime()
   return 0.
end

function VmChargeExchange:nonSlvrTime()
   return self.timers.nonSlvr
end

function VmChargeExchange:projectMaxwellTime()
   return self.collisionSlvr:projectMaxwellTime()
end

return VmChargeExchange

