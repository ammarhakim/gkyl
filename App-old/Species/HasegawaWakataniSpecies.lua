-- Gkyl ------------------------------------------------------------------------
--
-- Incompressible Euler equation in 2D.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct   = require "DataStruct"
local FluidSpecies = require "App.Species.FluidSpecies"
local HWEq         = require "Eq.HasegawaWakatani"
local LinearDecomp = require "Lib.LinearDecomp"
local Updater      = require "Updater"
local Mpi          = require "Comm.Mpi"
local Proto        = require "Lib.Proto"
local BasicBC      = require ("App.BCs.FluidBasic").FluidBasic
local lume         = require "Lib.lume"

local HWSpecies = Proto(FluidSpecies)

function HWSpecies:fullInit(appTbl)
   HWSpecies.super.fullInit(self, appTbl)

   self.nMoments = 2

   self.charge = -1.0   -- To interact with GkField.

   self.adiabaticity = assert(self.tbl.adiabaticity, "HasegawaWakataniSpecies: must specify the adiabaticity parameter with 'adiabaticity'.")
   self.gradientScale = assert(self.tbl.gradient, "HasegawaWakataniSpecies: must specify the normalized density gradient length scale (rho_s/L_n) parameter with 'gradient'.")
end

function HWSpecies:alloc(nRkDup)
   HWSpecies.super.alloc(self, nRkDup)   -- Call the FluidSpecies :alloc method.

   -- Allocate fields to store separate vorticity and density.
   self.vorticity = self:allocMoment()
   self.density   = self:allocMoment()

   -- These may be needed for interfacing with GkField.
   self.vorticityAux = self:allocMoment()
   self.densityAux   = self:allocMoment()

   -- Offsets needed to fetch specific moments from CartField
   -- containing the stepped moments (e.g. with :combineOffset).
   self.vorticityOff = 0
   self.densityOff   = 1*self.basis:numBasis()
   -- Package them into a single table for easier access.
   self.momOff = {self.vorticityOff,self.densityOff}
end

function HWSpecies:createSolver(field, externalField)
   -- Run the FluidSpecies 'createSolver()' to initialize the
   -- collisions (diffusion) solver.
   HWSpecies.super.createSolver(self, field, externalField)

   local hasE, hasB       = field:hasEB()
   local extHasE, extHasB = externalField:hasEB()

   -- Create updater to advance solution by one time-step.
   self.equation = HWEq {
      onGrid = self.grid,   adiabaticity  = self.adiabaticity,
      basis  = self.basis,  gradientScale = self.gradientScale,
   }

   self.solver = Updater.HyperDisCont {
      onGrid = self.grid,   cfl      = self.cfl,
      basis  = self.basis,  equation = self.equation,
   }

   self.advSolver = function(tCurr, momIn, em, momRhsOut)
      self.solver:advance(tCurr, {momIn, em}, {momRhsOut})
   end
end

function HWSpecies:getMomOff(momIdx)
   local offOut = momIdx and self.momOff[momIdx] or self.momOff
   return offOut
end

-- Nothing to calculate, just copy.
function HWSpecies:calcCouplingMomentsEvolve(tCurr, rkIdx)
   local momIn = self:rkStepperFields()[rkIdx]

   self.vorticity:combineOffset(1., momIn, self.vorticityOff)
   self.density:combineOffset(1., momIn, self.densityOff)
end

function HWSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self.tCurr = tCurr
   local momIn     = self:rkStepperFields()[inIdx]
   local momRhsOut = self:rkStepperFields()[outIdx]

   -- Clear RHS, because HyperDisCont set up with clearOut = false.
   momRhsOut:clear(0.0)

   -- Perform the collision (diffusion) update.
   for _, c in pairs(self.collisions) do
      c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      c:advance(tCurr, momIn, species, momRhsOut)
   end

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      local em = emIn[1]:rkStepperFields()[inIdx]
      self.advSolver(tCurr, momIn, em, momRhsOut)
   end

   for _, src in pairs(self.sources) do src:advance(tCurr, momIn, species, momRhsOut) end
end

function HWSpecies:fluidMoments()
   return { self.vorticity, self.density }
end

-- For interfacing with GkField, which we use to solve nabla^2_\perp(phi) = vorticity.
function HWSpecies:getNumDensity(rkIdx)
   -- If no rkIdx specified, assume numDensity has already been calculated.
   if rkIdx == nil then
      return self.vorticity
   end

   local momIn = self:rkStepperFields()[rkIdx]

   local tmStart = Time.clock()
   self.vorticityAux:combineOffset(1., momIn, self.vorticityOff)
   self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart

   return self.vorticityAux
end

function HWSpecies:getPolarizationWeight(linearized) return 1.0 end

function HWSpecies:suggestDt()
   return math.min(self.cfl/self.cflRateByCell:reduce('max')[1], GKYL_MAX_DOUBLE)
end

function HWSpecies:clearCFL()
   -- Clear cflRateByCell for next cfl calculation,
   self.cflRateByCell:clear(0.0)
end

return HWSpecies
