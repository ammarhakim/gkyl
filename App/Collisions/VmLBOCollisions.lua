-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Vlasov LB Collision operator
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local DataStruct = require "DataStruct"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local Updater = require "Updater"
local VmLBOconstNuEq = require "Eq.VmLBO"
local xsys = require "xsys"

-- VmLBOCollisions ---------------------------------------------------------------
--
-- Lenard-Bernstein Collision operator
-- Actually dates back to Lord Rayleigh, Philos. Mag. 32, 424 (1891).
--------------------------------------------------------------------------------

local VmLBOCollisions = Proto(CollisionsBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmLBOCollisions:init(tbl)
   self.tbl = tbl
   self._tmEvalMom = 0.0
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- previously stored table

   self.cfl = 0.1
   self.selfCollisions = xsys.pickBool(tbl.selfCollisions, true) -- by default, self collisions are on
   self.crossSpecies = tbl.crossSpecies
   self.collFreq = assert(tbl.collFreq,
			  "Updater.VmLBOCollisions: Must specify the collision frequency with 'collFreq'")
end

function VmLBOCollisions:setName(nm)
   self.name = nm
end
function VmLBOCollisions:setConfBasis(basis)
   self.confBasis = basis
end
function VmLBOCollisions:setConfGrid(grid)
   self.confGrid = grid
end
function VmLBOCollisions:setPhaseBasis(basis)
   self.phaseBasis = basis
end
function VmLBOCollisions:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function VmLBOCollisions:createSolver()
   self.vdim = self.phaseGrid:ndim() - self.confGrid:ndim()

   -- intemediate storage for output of collisions
   self.collOut = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Flow velocity in vdim directions.
   self.velocity = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis()*self.vdim,
      ghost         = {1, 1},
   }
   -- Thermal speed squared, vth=sqrt(T/m).
   self.vthSq = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Magnitude of kinetic energy density vector.
   self.kinEnergyDensM = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }
   self.thEnergyDens = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Zero-flux BCs in the velocity dimensions.
   local zfd = { }
   for d = 1, self.vdim do
      zfd[d] = self.confGrid:ndim() + d
   end

   -- Lenard-Bernestein equation.
   local vmLBOconstNuCalc = VmLBOconstNuEq {
      nu         = self.collFreq,
      phaseBasis = self.phaseBasis,
      confBasis  = self.confBasis,
   }
   self.collisionSlvr = Updater.HyperDisCont {
      onGrid             = self.phaseGrid,
      basis              = self.phaseBasis,
      cfl                = self.cfl,
      equation           = vmLBOconstNuCalc,
      onlyIncrement      = true,
      updateDirections   = zfd, -- only update velocity directions
      zeroFluxDirections = zfd,
   }

   -- Also need weak binary operations for primitive moments.
   -- Weak division of two configuration space fields.
   self.confDiv = Updater.CartFieldBinOp {
      onGrid     = self.confGrid,
      weakBasis  = self.confBasis,
      operation  = "Divide",
   }
   -- Dot product of two configuration space vector fields.
   self.confDotProduct = Updater.CartFieldBinOp {
      onGrid     = self.confGrid,
      weakBasis  = self.confBasis,
      operation  = "DotProduct",
   }
end

-- Computes primitive moments velocity and vth=sqrt(T/m) from zeroth,
-- first and second moments.
function VmLBOCollisions:primMoments(mom0, mom1, mom2)
   self.confDiv:advance(0, 0, {mom0, mom1}, {self.velocity})
   self.confDotProduct:advance(0, 0, {self.velocity, mom1}, {self.kinEnergyDensM})
   self.thEnergyDens:combine(1.0/self.vdim, mom2, -1.0/self.vdim, self.kinEnergyDensM)
   self.confDiv:advance(0, 0, {mom0, self.thEnergyDens}, {self.vthSq})
end

function VmLBOCollisions:forwardEuler(tCurr, dt, fIn, momIn, fOut)
   local tmEvalMomStart = Time.clock()
   self:primMoments(momIn[1], momIn[2], momIn[3]) -- compute primitive moments
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart

   -- compute increment from collisions
   local myStatus, myDt = self.collisionSlvr:advance(
      tCurr, dt, {fIn, self.velocity, self.vthSq}, {self.collOut})
   -- accumulate to output collisions
   fOut:accumulate(dt, self.collOut)

   return myStatus, myDt
end

function VmLBOCollisions:totalSolverTime()
   return self.collisionSlvr.totalTime + self._tmEvalMom
end

function VmLBOCollisions:evalMomTime()
   return self._tmEvalMom
end

return VmLBOCollisions
