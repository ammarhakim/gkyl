-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Ionization operator using Voronov
-- reaction rate.
-- See:
-- Voronov, G.S., 1997. A practical fit formula for ionization rate
-- coefficients of atoms and ions by electron impact: Z= 1–28. Atomic
-- Data and Nuclear Data Tables, 65(1), pp.1-35.
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

-- VmIonization -----------------------------------------------------------
--
-- Voronov ionization operator.
---------------------------------------------------------------------------

-- ............... IMPLEMENTATION OF DIAGNOSTICS ................. --
-- Diagnostics could be placed in a separate file if they balloon in
-- number. But if we only have one or two we can just place it here.

local vmIzDiagImpl = function()
   local _M0 = Proto(DiagsImplBase)
   function _M0:fullInit(diagApp, mySpecies, fieldIn, owner)
      self.field    = mySpecies:allocMoment()
      self.updater  = mySpecies.numDensityCalc
      self.owner    = owner
      self.done     = false
   end
   function _M0:getType() return "grid" end
   function _M0:advance(tm, inFlds, outFlds)
      self.updater:advance(tm, {self.owner.ionizSrc}, {self.field})
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
      self.updatersAux:advance(tm, {self.owner.ionizSrc}, {self.fieldAux})
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
      self.owner    = owner
      self.done     = false
   end
   function _source:getType() return "grid" end
   function _source:advance(tm, inFlds, outFlds)
      self.field:copy(self.owner.ionizSrc)
   end

   return {
      M0        = _M0,
      intM0     = _intM0,
      reactRate = _reactRate,
      source    = _source,
   }
end

-- .................... END OF DIAGNOSTICS ...................... --

local VmIonization = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmIonization:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmIonization:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously store table.

   self.cfl = 0.1
   self.collKind = "Ionization"

   self.collidingSpecies = assert(tbl.collideWith, "App.VmIonization: Must specify names of species to collide with in 'collideWith'.")

   -- Set these values to be consistent with other collision apps
   self.selfCollisions  = false
   self.crossCollisions = true              
   self.varNu           = false
   self.timeDepNu       = false
   self.collFreqs       = {1}
   
   self.collideNm = tbl.collideWith[1]
   
   self.elcNm  = assert(tbl.electrons, "App.VmIonization: Must specify electron species name in 'electrons'.")
   self.neutNm = assert(tbl.neutrals, "App.VmIonization: Must specify electron species name in 'neutrals'.")
   self.plasma = tbl.plasma
   self.mass   = tbl.elcMass
   self.charge = tbl.elemCharge

   self.timers = {nonSlvr = 0.}
end

function VmIonization:createDiagnostics(mySpecies, field)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = vmIzDiagImpl()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

function VmIonization:setName(nm)
   self.name = self.speciesName.."_"..nm
   self.collNm = nm
end
function VmIonization:setSpeciesName(nm) self.speciesName = nm end
function VmIonization:setCfl(cfl) self.cfl = cfl end
function VmIonization:setConfBasis(basis) self.confBasis = basis end
function VmIonization:setConfGrid(grid) self.confGrid = grid end
function VmIonization:setPhaseBasis(basis) self.phaseBasis = basis end
function VmIonization:setPhaseGrid(grid) self.phaseGrid = grid end

function VmIonization:createSolver(funcField)
   self.collisionSlvr = Updater.Ionization {
      confBasis  = self.confBasis,
      elcMass    = self.mass,
      elemCharge = self.charge,
      ionType    = self.plasma,
   }
   -- Reaction rate computed by each species for ease of species parallelization
   self.reactRate = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
      metaData      = { polyOrder = self.confBasis:polyOrder(),
                        basisType = self.confBasis:id() },
   }
   self.calcReactRate = function(tCurr, inFlds, outFlds)
      self.collisionSlvr:reactRateCoef(tCurr, inFlds, outFlds)
   end

   if (self.speciesName == self.elcNm) then
      -- Electrons need to calculate the ionization temperature
      self.vtSqIz =  DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
         metaData      = { polyOrder = self.confBasis:polyOrder(),
                           basisType = self.confBasis:id() },
      }
      self.calcIonizationTemp = function(tCurr, inFlds, outFlds)
         self.collisionSlvr:ionizationTemp(tCurr, inFlds, outFlds)
      end

      self.m0fMax =  DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
         metaData      = { polyOrder = self.confBasis:polyOrder(),
                           basisType = self.confBasis:id() },
      }
      self.m0mod =  DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
         metaData      = { polyOrder = self.confBasis:polyOrder(),
                           basisType = self.confBasis:id() },
      }
      self.fMaxElc = DataStruct.Field {
         onGrid        = self.phaseGrid,
         numComponents = self.phaseBasis:numBasis(),
         ghost         = {1, 1},
         metaData      = { polyOrder = self.confBasis:polyOrder(),
                           basisType = self.confBasis:id() },
      }
      self.sumDistF =  DataStruct.Field {
         onGrid        = self.phaseGrid,
         numComponents = self.phaseBasis:numBasis(),
         ghost         = {1, 1},
         metaData      = { polyOrder = self.confBasis:polyOrder(),
                           basisType = self.confBasis:id() },
      }
   end

   self.m0elc = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
      metaData      = { polyOrder = self.confBasis:polyOrder(),
                        basisType = self.confBasis:id() },
   }
   self.coefM0 = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
      metaData      = { polyOrder = self.confBasis:polyOrder(),
                        basisType = self.confBasis:id() },
   }
   self.ionizSrc = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData      = { polyOrder = self.phaseBasis:polyOrder(),
                        basisType = self.phaseBasis:id() },
   }
end

function VmIonization:advance(tCurr, fIn, species, fRhsOut)
   local tmNonSlvrStart = Time.clock()
   local reactRate = species[self.elcNm].collisions[self.collNm].reactRate
   local elcM0 = species[self.elcNm]:fluidMoments()[1]

   -- Check whether particle is electron, neutral or ion species.
   if (self.speciesName == self.elcNm) then
      -- Electrons.
      self.m0elc:copy(elcM0)
      local neutM0     = species[self.neutNm]:fluidMoments()[1]
      local neutU      = species[self.neutNm]:selfPrimitiveMoments()[1]
      local elcVtSq    = species[self.elcNm]:selfPrimitiveMoments()[2]
      local elcDistF   = species[self.elcNm]:getDistF()

      -- Calculate ioniz fMax
      self.calcIonizationTemp(tCurr, {elcVtSq}, {self.vtSqIz})
      species[self.speciesName].calcMaxwell:advance(tCurr, {self.m0elc, neutU, self.vtSqIz}, {self.fMaxElc})
      species[self.speciesName].numDensityCalc:advance(tCurr, {self.fMaxElc}, {self.m0fMax})
      species[self.speciesName].confWeakDivide:advance(tCurr, {self.m0fMax, self.m0elc}, {self.m0mod})
      species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {self.m0mod, self.fMaxElc}, {self.fMaxElc})
      
      self.sumDistF:combine(2.0,self.fMaxElc,-1.0,elcDistF)
 
      species[self.speciesName].confWeakMultiply:advance(tCurr, {reactRate, neutM0}, {self.coefM0})
      species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {self.coefM0, self.sumDistF}, {self.ionizSrc})
      -- Uncomment to test without fMaxwellian(Tiz).
      --self.confPhaseMult:advance(tCurr, {self.coefM0, elcDistF}, {self.ionizSrc})      

      fRhsOut:accumulate(1.0,self.ionizSrc)
   elseif (species[self.speciesName].charge == 0) then
      -- Neutrals.
      local neutDistF = species[self.neutNm]:getDistF()
      self.m0elc:copy(elcM0)
            
      species[self.speciesName].confWeakMultiply:advance(tCurr, {reactRate, self.m0elc}, {self.coefM0})
      species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {self.coefM0, neutDistF}, {self.ionizSrc})

      fRhsOut:accumulate(-1.0,self.ionizSrc)  
   else
      -- Ions.
      self.m0elc:copy(elcM0)
      local neutDistF = species[self.neutNm]:getDistF()

      species[self.speciesName].confWeakMultiply:advance(tCurr, {reactRate, self.m0elc}, {self.coefM0})
      species[self.speciesName].confPhaseWeakMultiply:advance(tCurr, {self.coefM0, neutDistF}, {self.ionizSrc})
      
      fRhsOut:accumulate(1.0,self.ionizSrc)
   end

   self.timers.nonSlvr = self.timers.nonSlvr + Time.clock() - tmNonSlvrStart
end
   
function VmIonization:write(tm, frame) end

function VmIonization:setCfl(cfl)
   self.cfl = cfl
end

function VmIonization:getIonizSrc()
   return self.ionizSrc
end

function VmIonization:slvrTime()
   local time = 0
   return time
end

function VmIonization:nonSlvrTime()
   return self.timers.nonSlvr
end

function VmIonization:projectMaxwellTime()
   local time = self.confMult.projectMaxwellTime()
   time = time + self.confPhaseMult.projectMaxwellTime
   return time
end

return VmIonization

