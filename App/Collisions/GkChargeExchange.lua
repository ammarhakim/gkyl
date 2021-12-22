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

-- GkChargeExchange  --------------------------------------------------------
--
-- Charge Exchange Operator
-----------------------------------------------------------------------------

-- ............... IMPLEMENTATION OF DIAGNOSTICS ................. --
-- Diagnostics could be placed in a separate file if they balloon in
-- number. But if we only have one or two we can just place it here.

local gkCxDiagImpl = function()
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

local GkChargeExchange = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkChargeExchange:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkChargeExchange:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously store table.

   self.collKind = "CX"

   self.collidingSpecies = assert(tbl.collideWith, "App.GkChargeExchange: Must specify names of species to collide with in 'collideWith'.")

   self.collideNm   = tbl.collideWith[1]

   self.plasma      = assert(tbl.plasma, "App.GkChargeExchange: Must specify plasma species in 'plasma' ('H', 'D', or 'Ne')")
   
   self.ionNm       = tbl.ions
   self.neutNm      = tbl.neutrals
   self.iMass       = tbl.ionMass
   self.nMass       = tbl.neutMass
   self.charge      = tbl.charge

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

function GkChargeExchange:createDiagnostics(mySpecies, field)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = gkCxDiagImpl()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

function GkChargeExchange:setName(nm)
   self.name = self.speciesName.."_"..nm
   self.collNm = nm
end
function GkChargeExchange:setSpeciesName(nm) self.speciesName = nm end
function GkChargeExchange:setCfl(cfl) self.cfl = cfl end
function GkChargeExchange:setConfBasis(basis) self.confBasis = basis end
function GkChargeExchange:setConfGrid(grid) self.confGrid = grid end
function GkChargeExchange:setPhaseBasis(basis) self.phaseBasis = basis end
function GkChargeExchange:setPhaseGrid(grid) self.phaseGrid = grid end

function GkChargeExchange:createSolver(funcField)
   self.sourceCX =  DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.M0iDistFn =  DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.M0nDistFi =  DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.diffDistF =  DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.collisionSlvr = Updater.ChargeExchange {
      onGrid         = self.confGrid,
      confBasis      = self.confBasis,
      phaseBasis     = self.phaseBasis,
      charge         = self.charge, 
      a              = self.a,
      b              = self.b,
   }
   
   if (self.speciesName == self.ionNm) then --ions
      self.fMaxNeut =  DataStruct.Field {
	 onGrid        = self.phaseGrid,
	 numComponents = self.phaseBasis:numBasis(),
	 ghost         = {1, 1},
	 metaData = {
	    polyOrder = self.phaseBasis:polyOrder(),
	    basisType = self.phaseBasis:id()
	 },
      }
      self.reactRate =  DataStruct.Field {
	 onGrid        = self.confGrid,
	 numComponents = self.confBasis:numBasis(),
	 ghost         = {1, 1},
	 metaData = {
	    polyOrder = self.confBasis:polyOrder(),
	    basisType = self.confBasis:id()
	 },
      }
   else --neutrals
      self.fMaxIon =  DataStruct.Field {
	 onGrid        = self.phaseGrid,
	 numComponents = self.phaseBasis:numBasis(),
	 ghost         = {1, 1},
	 metaData = {
	    polyOrder = self.phaseBasis:polyOrder(),
	    basisType = self.phaseBasis:id()
	 },
      }
      self.ionU =  DataStruct.Field {
	 onGrid        = self.confGrid,
	 numComponents = self.confBasis:numBasis()*self.confGrid:ndim(),
	 ghost         = {1, 1},
	 metaData = {
	    polyOrder = self.confBasis:polyOrder(),
	    basisType = self.confBasis:id()
	 },
      }
      self.bhat = species[self.neutNm].bHat
   end
end

function GkChargeExchange:advance(tCurr, fIn, species, fRhsOut)
   local tmNonSlvrStart = Time.clock()
   local reactRate = species[self.ionNm].collisions[self.collNm].reactRate

   -- Identify species and accumulate.
   if (self.speciesName == self.ionNm) then

      local neutM0   = species[self.neutNm]:fluidMoments()[1]
      local neutUpar = species[self.neutNm].uPar
      local neutVtSq = species[self.neutNm].vtSqGk
      local ionM0    = species[self.ionNm]:fluidMoments()[1]
      local ionDistF = species[self.ionNm]:getDistF()

      species[self.speciesName].calcMaxwell:advance(tCurr,
         {neutM0, neutUpar, neutVtSq, species[self.speciesName].bmag}, {self.fMaxNeut})

      species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {ionM0, self.fMaxNeut}, {self.M0iDistFn})
      species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {neutM0, ionDistF}, {self.M0nDistFi})
      self.diffDistF:combine(1.0, self.M0iDistFn, -1.0, self.M0nDistFi)
      species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {reactRate, self.diffDistF}, {self.sourceCX})

      fRhsOut:accumulate(1.0,self.sourceCX)

   elseif (self.speciesName == self.neutNm) then

      local ionM0     = species[self.ionNm]:fluidMoments()[1]
      local ionUpar   = species[self.ionNm]:selfPrimitiveMoments()[1] 
      local ionVtSq   = species[self.ionNm]:selfPrimitiveMoments()[2]
      local neutM0    = species[self.neutNm]:fluidMoments()[1]
      local neutDistF = species[self.neutNm]:getDistF()
      
      species[self.speciesName].confMult:advance(tCurr, {ionUpar,self.bhat}, {self.ionU})
      species[self.speciesName].calcMaxwell:advance(tCurr, {ionM0, self.ionU, ionVtSq}, {self.fMaxIon})

      species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {ionM0, neutDistF}, {self.M0iDistFn})
      species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {neutM0, self.fMaxIon}, {self.M0nDistFi})
      self.diffDistF:combine(1.0, self.M0iDistFn, -1.0, self.M0nDistFi)
      species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {reactRate, self.diffDistF}, {self.sourceCX})
      
      fRhsOut:accumulate(-self.iMass/self.nMass,self.sourceCX)

   end

   self.timers.nonSlvr = self.timers.nonSlvr + Time.clock() - tmNonSlvrStart
end

function GkChargeExchange:write(tm, frame) end

function GkChargeExchange:slvrTime()
   return 0
end

function GkChargeExchange:nonSlvrTime()
   return self.timers.nonSlvr
end

function GkChargeExchange:projectMaxwellTime()
   return self.collisionSlvr:projectMaxwellTime()
end

return GkChargeExchange

