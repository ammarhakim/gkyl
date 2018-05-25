-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Vlasov LB Collision operator
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto          = require "Lib.Proto"
local Updater        = require "Updater"
local CollisionsBase = require "App.Collisions.CollisionsBase"
local VmLBOconstNuEq = require "Eq.VmLBO"
local DataStruct     = require "DataStruct"
local Time           = require "Lib.Time"

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

   -- timings.
   self._tmEvalMom = 0.0
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmLBOCollisions:fullInit(collTbl)
   local tbl = self.tbl -- previously store table.

   self.cfl         = 0.1
   self.speciesList = tbl.species
   self.collFreq    = tbl.collFreq

   assert(#self.speciesList == #self.collFreq,
	  "'nu' must be defined for each 'species'")
end

function VmLBOCollisions:setName(nm)
   self.name = nm
end

function VmLBOCollisions:setConfBasis(basis)
   self.confBasis = basis
end
function VmLBOCollisions:setConfGrid(cgrid)
   self.confGrid = cgrid
end

function VmLBOCollisions:setPhaseBasis(species)
   self.phaseBasis = {}
   for _, nm in pairs(self.speciesList) do
      self.phaseBasis[nm] = species[nm].basis
   end
end

function VmLBOCollisions:setPhaseGrid(species)
   self.phaseGrid = {}
   for _, nm in pairs(self.speciesList) do
      self.phaseGrid[nm] = species[nm].grid
   end
end

-- methods for Bgk collisions object

function VmLBOCollisions:createSolver(species)
   local confBasis  = nil
   local confGrid   = nil
   local phaseBasis = nil
   local phaseGrid  = nil
   local zfd = { }
   for _, nm in pairs(self.speciesList) do
      confBasis  = species[nm].confBasis
      confGrid   = species[nm].confGrid
      phaseBasis = species[nm].basis
      phaseGrid  = species[nm].grid
   end

   self.vdim = phaseGrid:ndim()-confGrid:ndim()

   -- intemediate storage for output of collisions
   self.collOut = DataStruct.Field {
      onGrid        = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost         = {1, 1},
   }   
   
   -- Flow velocity in vdim directions.
   self.velocity = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis()*self.vdim,
      ghost         = {1, 1},
   }
   -- Thermal speed squared, vth=sqrt(T/m).
   self.vthSq = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Magnitude of kinetic energy density vector.
   self.kinEnergyDensM = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   self.thEnergyDens = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Lenard-Bernestein equation.
   local vmLBOconstNuCalc = VmLBOconstNuEq {
      nu         = self.collFreq[1],
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
   }
   -- Zero-flux BCs in the velocity dimensions.
   local zfd = { }
   for d = 1, (phaseGrid:ndim()-confGrid:ndim()) do
      zfd[d] = confGrid:ndim()+d
   end
   self.collisionSlvr = Updater.HyperDisCont {
      onGrid             = phaseGrid,
      basis              = phaseBasis,
      cfl                = self.cfl,
      equation           = vmLBOconstNuCalc,
      onlyIncrement      = true,
      zeroFluxDirections = zfd,
   }
   -- Also need weak binary operations for primitive moments.
   -- Weak division of two configuration space fields.
   self.confDiv = Updater.CartFieldBinOp {
      onGrid     = confGrid,
      weakBasis  = confBasis,
      operation  = "Divide",
   }
   -- Dot product of two configuration space vector fields.
   self.confDotProduct = Updater.CartFieldBinOp {
      onGrid     = confGrid,
      weakBasis  = confBasis,
      operation  = "DotProduct",
   }
end

-- This function computes the primitive moments velocity
-- and vth=sqrt(T/m) from the zeroth, first and second moments.
function primMoments(mom0,mom1,mom2)
   -- Compute the flow velocity using weak division of mom0 and mom1.
   self.confDiv:advance(0.,0.,{mom0,mom1},{self.velocity})
   -- Compute kinetic energy using weak multiplication of u and mom1.
   self.confDotProduct:advance(0.,0.,{self.velocity,mom1},{self.kinEnergyDensM})
   -- Thermal energy density =  mom2 - u*mom1.
   self.thEnergyDens:combine(1.0/self.vdim, mom2, -1.0/self.vdim, self.kinEnergyDensM)
   -- Compute thermal speed squared via weak division of thEnergy and mom0.
   self.confDiv:advance(0.,0.,{mom0,self.thEnergyDens},{self.vthSq})
end

function VmLBOCollisions:forwardEuler(tCurr, dt, idxIn, idxOut, species)
   -- Timings
   local tmEvalMomStart = 0.0
   local nm = self.speciesList[1] -- for now, only self-collisions
   local spMomFields = species[nm]:fluidMoments()
   local spInField = species[nm]:rkStepperFields()[idxIn]
   local spOutField = species[nm]:rkStepperFields()[idxOut]

   tmEvalMomStart = Time.clock()
   -- compute primitive moments
   self.confDiv:advance(0.,0.,{spMomFields[1], spMomFields[2]},{self.velocity})
   self.confDotProduct:advance(0.,0.,{self.velocity,spMomFields[2]},{self.kinEnergyDensM})
   self.thEnergyDens:combine(1.0/self.vdim, spMomFields[3], -1.0/self.vdim, self.kinEnergyDensM)
   self.confDiv:advance(0.,0.,{spMomFields[1],self.thEnergyDens},{self.vthSq})

   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
   local myStatus, myDt = self.collisionSlvr:advance(
      tCurr, dt, {spInField, self.velocity, self.vthSq}, {self.collOut})
   
   spOutField:accumulate(dt, self.collOut) -- accumulate output from collisions

   return myStatus, myDt
end

function VmLBOCollisions:totalSolverTime()
   return self.collisionSlvr.totalTime 
end

function VmLBOCollisions:evalMomTime()
--    return self.collisionSlvr:evalMomTime()
   return self._tmEvalMom
end

function VmLBOCollisions:projectMaxwellTime()
--   return self.collisionSlvr:projectMaxwellTime()
   return 0.0
end

return VmLBOCollisions
