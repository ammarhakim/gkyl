-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: BGK Collision operator
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local DataStruct = require "DataStruct"
local Proto = require "Lib.Proto"
local Updater = require "Updater"

-- BgkCollisions ---------------------------------------------------------------
--
-- Bhatnagar-Gross-Krook Collision operator
--------------------------------------------------------------------------------

local BgkCollisions = Proto(CollisionsBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function BgkCollisions:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function BgkCollisions:fullInit(collTbl)
   local tbl = self.tbl -- previously store table

   self.species = assert(tbl.species,
			 "Updater.BgkCollisions: Must specify species with 'species'")
   self.crossSpecies = tbl.crossSpecies
   self.collFreq = assert(tbl.collFreq,
			  "Updater.BgkCollisions: Must specify the collision frequency with 'collFreq'")
end

function BgkCollisions:setName(nm)
   self.name = nm
end

function BgkCollisions:setConfBasis(basis)
   self.confBasis = basis
end
function BgkCollisions:setConfGrid(cgrid)
   self.confGrid = cgrid
end

function BgkCollisions:setPhaseBasis(species)
   self.phaseBasis = species[self.species].basis
end

function BgkCollisions:setPhaseGrid(species)
   self.phaseGrid = species[self.species].grid
end

-- methods for Bgk collisions object

function BgkCollisions:createSolver()
   self.fMaxwell = DataStruct.Field {
      onGrid = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost = {1, 1},
   }
   self.maxwellian = Updater.MaxwellianOnBasis {
      onGrid = self.confGrid,
      confGrid = self.confGrid,
      confBasis = self.confBasis,
      phaseGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
   }
   self.collisionSlvr = Updater.BgkCollisions {
      onGrid = self.confGrid,
      confGrid = self.confGrid,
      confBasis = self.confBasis,
      phaseGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      collFreq = self.collFreq,
   }
end

function BgkCollisions:forwardEuler(tCurr, dt, idxIn, idxOut, species)
   local momFields = species[self.species]:fluidMoments()
   self.maxwellian:advance(tCurr, dt, {momFields[1], momFields[2], momFields[3]},
			   {self.fMaxwell})
   return self.collisionSlvr:advance(tCurr, dt,
				     {species[self.species]:rkStepperFields()[idxIn],
				      self.fMaxwell},
				     {species[self.species]:rkStepperFields()[idxOut]})
end

function BgkCollisions:totalSolverTime()
   return self.collisionSlvr.totalTime 
end

function BgkCollisions:evalMomTime()
   return self.collisionSlvr:evalMomTime()
end

function BgkCollisions:projectMaxwellTime()
   return self.collisionSlvr:projectMaxwellTime()
end

return BgkCollisions
