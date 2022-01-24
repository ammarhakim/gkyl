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

   self.adiabaticity = assert(self.tbl.adiabaticity, "HasegawaWakataniSpecies: must specify the adiabaticity parameter with 'adiabaticity'.")
   self.gradientScale = assert(self.tbl.gradient, "HasegawaWakataniSpecies: must specify the normalized density gradient length scale (rho_s/L_n) parameter with 'gradient'.")
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

   self.advSolver = function(tCurr, fIn, em, fRhsOut)
      self.solver:advance(tCurr, {fIn, em}, {fRhsOut})
   end
end

-- Nothing to calculate, just copy.
function HWSpecies:calcCouplingMomentsEvolve(tCurr, rkIdx)
   local fIn = self:rkStepperFields()[rkIdx]
   self.couplingMoments:copy(fIn)
end

function HWSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self.tCurr = tCurr
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   -- Clear RHS, because HyperDisCont set up with clearOut = false.
   fRhsOut:clear(0.0)

   -- Perform the collision (diffusion) update.
   for _, c in pairs(self.collisions) do
      c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      c:advance(tCurr, fIn, species, fRhsOut)
   end

   -- Complete the field solve.
   emIn[1]:phiSolve(tCurr, species, inIdx, outIdx)

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      local em = emIn[1]:rkStepperFields()[inIdx]
      self.advSolver(tCurr, fIn, em, fRhsOut)
   end

   for _, src in pairs(self.sources) do src:advance(tCurr, fIn, species, fRhsOut) end
end

function HWSpecies:fluidMoments() return { self.couplingMoments } end

-- For interfacing with GkField.
function HWSpecies:getNumDensity(rkIdx)
   if rkIdx == nil then return self.couplingMoments 
   else 
      self.couplingMoments:copy(self:rkStepperFields()[rkIdx])
      return self.couplingMoments
   end
end

function HWSpecies:suggestDt()
   return math.min(self.cfl/self.cflRateByCell:reduce('max')[1], GKYL_MAX_DOUBLE)
end

function HWSpecies:clearCFL()
   -- Clear cflRateByCell for next cfl calculation,
   self.cflRateByCell:clear(0.0)
end

return HWSpecies
