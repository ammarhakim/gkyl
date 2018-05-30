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
local xsys = require "xsys"

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
function BgkCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- previously store table

   self.selfCollisions = xsys.pickBool(tbl.selfCollisions, true) -- by default, self collisions are on
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
function BgkCollisions:setConfGrid(grid)
   self.confGrid = grid
end
function BgkCollisions:setPhaseBasis(basis)
   self.phaseBasis = basis
end
function BgkCollisions:setPhaseGrid(grid)
   self.phaseGrid = grid
end

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

function BgkCollisions:forwardEuler(tCurr, dt, fIn, momIn, fOut)
   self.maxwellian:advance(tCurr, dt,
			   {momIn[1], momIn[2], momIn[3]},
			   {self.fMaxwell})
   return self.collisionSlvr:advance(tCurr, dt,
				     {fIn, self.fMaxwell},
				     {fOut})
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
