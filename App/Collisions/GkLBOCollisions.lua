-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Gyrokinetic LB Collision operator
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local DataStruct     = require "DataStruct"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local GkLBOconstNuEq = require "Eq.GkLBO"
local xsys           = require "xsys"

-- GkLBOCollisions ---------------------------------------------------------------
--
-- Lenard-Bernstein Collision operator
-- Actually dates back to Lord Rayleigh, Philos. Mag. 32, 424 (1891).
--------------------------------------------------------------------------------

local GkLBOCollisions = Proto(CollisionsBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkLBOCollisions:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- previously stored table

   self.cfl = 0.0 -- will be replaced
   self.selfCollisions = xsys.pickBool(tbl.selfCollisions, true) -- by default, self collisions are on
   self.crossSpecies = tbl.crossSpecies
   self.collFreq = assert(
      tbl.collFreq, "Updater.GkLBOCollisions: Must specify the collision frequency with 'collFreq'")
   self.mass = speciesTbl.mass
   self.tmEvalMom = 0.0
end

function GkLBOCollisions:setName(nm)
   self.name = nm
end
function GkLBOCollisions:setSpeciesName(nm)
   self.speciesName = nm
end

function GkLBOCollisions:setCfl(cfl)
   self.cfl = cfl -- what should this be? - AHH
end
function GkLBOCollisions:setConfBasis(basis)
   self.confBasis = basis
end
function GkLBOCollisions:setConfGrid(grid)
   self.confGrid = grid
end
function GkLBOCollisions:setPhaseBasis(basis)
   self.phaseBasis = basis
end
function GkLBOCollisions:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function GkLBOCollisions:createSolver(funcField)
   self.vdim = self.phaseGrid:ndim() - self.confGrid:ndim()

   -- intemediate storage for output of collisions
   self.collOut = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Parallel flow velocity.
   self.uPar = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Thermal speed squared, vth=sqrt(T/m).
   self.vthSq = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Inverse of background magnetic field.
   self.bmag    = funcField.geo.bmag
   -- Inverse of background magnetic field.
   self.bmagInv = funcField.geo.bmagInv
      
   -- Zero-flux BCs in the velocity dimensions.
   local zfd = { }
   for d = 1, self.vdim do
      zfd[d] = self.confGrid:ndim() + d
   end

   -- Lenard-Bernestein equation.
   local gkLBOconstNuCalc = GkLBOconstNuEq {
      nu         = self.collFreq,
      phaseBasis = self.phaseBasis,
      confBasis  = self.confBasis,
      mass       = self.mass,
   }
   self.collisionSlvr = Updater.HyperDisCont {
      onGrid             = self.phaseGrid,
      basis              = self.phaseBasis,
      cfl                = self.cfl,
      equation           = gkLBOconstNuCalc,
      onlyIncrement      = true,
      updateDirections   = zfd, -- only update velocity directions
      zeroFluxDirections = zfd,
   }
   self.primMomSelf = Updater.SelfPrimMoments {
      onGrid     = self.confGrid,
      phaseGrid  = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confBasis  = self.confBasis,
      gkfacs     = {self.mass, self.bmag},
      operator   = "GkLBO",
   }
end

function GkLBOCollisions:forwardEuler(tCurr, dt, fIn, species, fOut)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   local selfMom = species[self.speciesName]:fluidMoments()

   if self.selfCollisions then
      local tmEvalMomStart = Time.clock()
      -- Compute primitive moments velocity and vthSq=T/m from zeroth,
      -- first and second moments, and distribution function.
      self.primMomSelf:advance(0.0, 0.0, {selfMom[1], selfMom[2], selfMom[3],fIn},
                                         {self.uPar,self.vthSq})
      self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

      -- LBO actually takes the product nu*velocity and nu*vthSq.
      -- For now just scale, but if nu is spatially dependent use
      -- binOp.Multiply.
      self.uPar:scale(self.collFreq)
      self.vthSq:scale(self.collFreq)
      
      -- compute increment from collisions and accumulate it into output
      local tmpStatus, tmpDt = self.collisionSlvr:advance(
	 tCurr, dt, {fIn, self.bmagInv, self.uPar, self.vthSq}, {self.collOut})
      status = status and tmpStatus
      dtSuggested = math.min(dtSuggested, tmpDt)

      fOut:accumulate(dt, self.collOut)
   end
   if self.crossSpecies then
      -- Insert cross collisions here!
   end
   return status, dtSuggested
end

function GkLBOCollisions:write(tm, frame)
   self.uPar:scale(1.0/self.collFreq)
   self.vthSq:scale(1.0/self.collFreq)
   self.uPar:write(string.format("%s_%s_%d.bp", self.speciesName, "uPar", frame), tm, frame)
   self.vthSq:write(string.format("%s_%s_%d.bp", self.speciesName, "vthSq", frame), tm, frame)
end

function GkLBOCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.tmEvalMom
end

return GkLBOCollisions
