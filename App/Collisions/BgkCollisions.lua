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

   self.cfl = 0.1 -- some default value
   self.selfCollisions = xsys.pickBool(tbl.selfCollisions, true) -- by default, self collisions are on
   self.crossSpecies = tbl.crossSpecies
   self.collFreq = assert(
      tbl.collFreq, "Updater.BgkCollisions: Must specify the collision frequency with 'collFreq'")
end

function BgkCollisions:setName(nm)
   self.name = nm
end
function BgkCollisions:setSpeciesName(nm)
   self.speciesName = nm
end

function BgkCollisions:setCfl(cfl)
   self.cfl = cfl
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
   -- temporary fields for the primitive moments
   self.numVelDims = self.phaseGrid:ndim() - self.confGrid:ndim()
   self.bulkVelocity = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.confBasis:numBasis()*self.numVelDims,
      ghost = {1, 1},
   }
   self.thermVelocity2 = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost = {1, 1},
   }
   self._kinEnergyDens = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost = {1, 1},
   }
   self._thermEnergyDens = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost = {1, 1},
   }
   -- Updaters for the primitive moments
   self.confDiv = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
   }
   self.confDotProduct = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "DotProduct",
   }

   -- Maxwellian field
   self.fMaxwell = DataStruct.Field {
      onGrid = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost = {1, 1},
   }
   -- Maxwellian solver
   self.maxwellian = Updater.MaxwellianOnBasis {
      onGrid = self.confGrid,
      confGrid = self.confGrid,
      confBasis = self.confBasis,
      phaseGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
   }

   -- BGK Collision solver itself
   self.collisionSlvr = Updater.BgkCollisions {
      onGrid = self.confGrid,
      confGrid = self.confGrid,
      confBasis = self.confBasis,
      phaseGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      collFreq = self.collFreq,
   }
end

function BgkCollisions:primMoments(M0, M1i, M2, u, vth2)
   self.confDiv:advance(0., 0., {M0, M1i}, {u})
   self.confDotProduct:advance(0., 0., {u, M1i},
			       {self._kinEnergyDens})
   self._thermEnergyDens:combine(1.0/self.numVelDims, M2,
				   -1.0/self.numVelDims, self._kinEnergyDens)
   self.confDiv:advance(0., 0., {M0, self._thermEnergyDens}, {vth2})
end

function BgkCollisions:forwardEuler(tCurr, dt, fIn, species, fOut)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   local selfMom = species[self.speciesName]:fluidMoments()

   self:primMoments(selfMom[1], selfMom[2], selfMom[3], 
		    self.bulkVelocity, self.thermVelocity2)
   if self.selfCollisions then
      self.maxwellian:advance(
	 tCurr, dt, {selfMom[1], self.bulkVelocity, self.thermVelocity2},
	 {self.fMaxwell})
      local tmpStatus, tmpDt = self.collisionSlvr:advance(
	 tCurr, dt, {fIn, self.fMaxwell}, {fOut})
      status = status and tmpStatus
      dtSuggested = math.min(dtSuggested, tmpDt)
   end
   if self.crossSpecies then
      -- Insert cross collisions here!
   end
   return status, dtSuggested
end

function BgkCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.maxwellian.totalTime
end

return BgkCollisions
