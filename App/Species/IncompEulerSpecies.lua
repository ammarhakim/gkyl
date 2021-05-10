-- Gkyl ------------------------------------------------------------------------
--
-- Incompressible Euler equation in 2D.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct    = require "DataStruct"
local FluidSpecies  = require "App.Species.FluidSpecies"
local IncompEulerEq = require "Eq.IncompEuler"
local LinearDecomp  = require "Lib.LinearDecomp"
local Updater       = require "Updater"
local Mpi           = require "Comm.Mpi"
local Proto         = require "Lib.Proto"
local BasicBC       = require "App.BCs.IncompEulerBasic"
local lume          = require "Lib.lume"

local IncompEulerSpecies = Proto(FluidSpecies)

-- ............. Backwards compatible treatment of BCs .....................--
-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB    = 1
local SP_BC_COPY      = 2
local SP_BC_ZEROFLUX  = 3
IncompEulerSpecies.bcAbsorb    = SP_BC_ABSORB       -- Absorb all particles.
IncompEulerSpecies.bcCopy      = SP_BC_COPY         -- Copy skin-cell values.
IncompEulerSpecies.bcZeroFlux  = SP_BC_ZEROFLUX     -- Zero flux.

function IncompEulerSpecies:makeBcApp(bcIn)
   local bcOut
   if bcIn == SP_BC_COPY then
      bcOut = BasicBC{kind="copy"}
   elseif bcIn == SP_BC_ABSORB then
      bcOut = BasicBC{kind="absorb"}
   elseif bcIn == SP_BC_ZEROFLUX then
      bcOut = BasicBC{kind="zeroFlux"}
   end
   return bcOut
end
-- ............. End of backwards compatibility for BCs .....................--

function IncompEulerSpecies:fullInit(appTbl)
   IncompEulerSpecies.super.fullInit(self, appTbl)

   self.nMoments = 1
end

function IncompEulerSpecies:createSolver(hasE, hasB)
   -- Run the FluidSpecies 'createSolver()' to initialize the
   -- collisions (diffusion) solver.
   IncompEulerSpecies.super.createSolver(self)

   -- Retrieve the zero-flux BCs (if any) from BC objects.
   local zeroFluxDirs = {}
   for _, bc in ipairs(self.nonPeriodicBCs) do 
      local bcDir = bc:getDir()
      if bc:getKind() == "zeroFlux" and not lume.any(zeroFluxDirs, function(e) return e==bcDir end) then
         table.insert(zeroFluxDirs, bcDir)
      end
   end

   if self.evolveCollisionless then
      -- Create updater to advance solution by one time-step.
      self.equation = IncompEulerEq {
         onGrid     = self.grid,
         basis      = self.basis,
         charge     = self.charge,
         positivity = self.positivity,
      }

      self.solver = Updater.HyperDisCont {
         onGrid             = self.grid,
         basis              = self.basis,
         cfl                = self.cfl,
         equation           = self.equation,
         zeroFluxDirections = zeroFluxDirs,
      }
   end
end

-- Nothing to calculate, just copy.
function IncompEulerSpecies:calcCouplingMoments(tCurr, rkIdx)
   local fIn = self:rkStepperFields()[rkIdx]
   self.couplingMoments:copy(fIn)
end

function IncompEulerSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self.tCurr = tCurr
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      local em = emIn[1]:rkStepperFields()[inIdx]
      if self.positivityRescale then
         self.posRescaler:advance(tCurr, {fIn}, {self.momPos}, false)
         self.solver:advance(tCurr, {self.momPos, em}, {fRhsOut})
      else
         self.solver:advance(tCurr, {fIn, em}, {fRhsOut})
      end
   else
      fRhsOut:clear(0.0) -- No RHS.
   end

   -- Perform the collision (diffusion) update.
   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
         c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
         c:advance(tCurr, fIn, species, fRhsOut)
      end
   end

   for _, src in pairs(self.sources) do src:advance(tCurr, fIn, species, fRhsOut) end
end

function IncompEulerSpecies:fluidMoments()
   return { self.couplingMoments }
end

-- For interfacing with GkField.
function IncompEulerSpecies:getNumDensity(rkIdx)
   if rkIdx == nil then return self.couplingMoments 
   else 
      self.couplingMoments:copy(self:rkStepperFields()[rkIdx])
      return self.couplingMoments
   end
end

function IncompEulerSpecies:suggestDt()
   return math.min(self.cfl/self.cflRateByCell:reduce('max')[1], GKYL_MAX_DOUBLE)
end

function IncompEulerSpecies:clearCFL()
   -- Clear cflRateByCell for next cfl calculation,
   self.cflRateByCell:clear(0.0)
end

return IncompEulerSpecies
