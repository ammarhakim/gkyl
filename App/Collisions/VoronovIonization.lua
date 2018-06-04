-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Voronov ionization operator
-- See:
-- Voronov, G.S., 1997. A practical fit formula for ionization rate
-- coefficients of atoms and ions by electron impact: Z= 1–28. Atomic
-- Data and Nuclear Data Tables, 65(1), pp.1-35.
--
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local Updater = require "Updater"
local CollisionsBase = require "App.Collisions.CollisionsBase"

-- VoronovIonization -----------------------------------------------------------
--
-- Voronov ionization operator
--------------------------------------------------------------------------------

local VoronovIonization = Proto(CollisionsBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function VoronovIonization:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function VoronovIonization:fullInit(collTbl)
   local tbl = self.tbl -- previously store table

   self.elcNm = tbl.electrons
   self.ionNm = tbl.ions
   self.neutOnElcNm = tbl.neutralOnElc
   self.neutOnIonNm = tbl.neutralOnIon
   self.plasma = tbl.plasma

   if self.plasma == "H" then
      self._E = 13.6
      self._P = 0
      self._A = 2.91e-14
      self._K = 0.39
      self._X = 0.232
   end
end

function VoronovIonization:setName(nm)
   self.name = nm
end

function VoronovIonization:setConfBasis(basis)
   self.confBasis = basis
end
function VoronovIonization:setConfGrid(cgrid)
   self.confGrid = cgrid
end

function VoronovIonization:setPhaseBasis(species)
   self.phaseBasis = {}
   self.phaseBasis['elc'] = species[self.elcNm].basis
   self.phaseBasis['ion'] = species[self.ionNm].basis
end

function VoronovIonization:setPhaseGrid(species)
   self.phaseGrid = {}
   self.phaseGrid['elc'] = species[self.elcNm].grid
   self.phaseGrid['ion'] = species[self.ionNm].grid
end

-- methods for Bgk collisions object

function VoronovIonization:createSolver(species)
   self.collisionSlvr = Updater.VoronovIonization {
      onGrid = species[self.elcNm].confGrid,
      confGrid = species[self.elcNm].confGrid,
      confBasis = species[self.elcNm].confBasis,
      phaseGrid = species[self.elcNm].grid,
      phaseBasis = species[self.elcNm].basis,
      elemCharge = math.abs(species[self.elcNm]:getCharge()),
      elcMass = species[self.elcNm]:getMass(),
      -- Voronov parameters
      A = self._A,
      E = self._E,
      K = self._K,
      P = self._P,
      X = self._X,
   }
end

function VoronovIonization:forwardEuler(tCurr, dt, idxIn, outIdx, species)
   local elcMomFields = species[self.elcNm]:fluidMoments()
   local spOutFields =  {}
   -- for nm, sp in pairs(species) do
   --    spOutFields[nm] = sp:rkStepperFields()[outIdx]
   -- end
   spOutFields['elc'] = species[self.elcNm]:rkStepperFields()[outIdx]
   spOutFields['ion'] = species[self.ionNm]:rkStepperFields()[outIdx]
   spOutFields['neutOnElc'] = species[self.neutOnElcNm]:rkStepperFields()[outIdx]
   spOutFields['neutOnIon'] = species[self.neutOnIonNm]:rkStepperFields()[outIdx]
   return self.collisionSlvr:advance(tCurr, dt, elcMomFields, spOutFields)
end

function VoronovIonization:totalSolverTime()
   return self.collisionSlvr.totalTime
end

function VoronovIonization:evalMomTime()
   return self.collisionSlvr:evalMomTime()
end

function VoronovIonization:projectMaxwellTime()
   return self.collisionSlvr:projectMaxwellTime()
end

return VoronovIonization
