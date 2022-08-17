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
local BasicBC       = require ("App.BCs.FluidBasic").FluidBasic
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
      print("IncompEulerSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="copy"}
   elseif bcIn == SP_BC_ABSORB then
      print("IncompEulerSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="absorb"}
   elseif bcIn == SP_BC_ZEROFLUX or bcIn.tbl.kind=="zeroFlux" then
      bcOut = "zeroFlux"
      table.insert(self.zeroFluxDirections, dir)
   end
   return bcOut
end
-- ............. End of backwards compatibility for BCs .....................--

function IncompEulerSpecies:fullInit(appTbl)
   IncompEulerSpecies.super.fullInit(self, appTbl)

   self.nMoments = 1
end

function IncompEulerSpecies:createSolver(field, externalField)
   -- Run the FluidSpecies 'createSolver()' to initialize the
   -- collisions (diffusion) solver.
   IncompEulerSpecies.super.createSolver(self, field, externalField)

   local hasE, hasB       = field:hasEB()
   local extHasE, extHasB = externalField:hasEB()

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
         zeroFluxDirections = self.zeroFluxDirections,
      }

      if self.positivityRescale then
         self.advSolver = function(tCurr, fIn, em, fRhsOut)
            self.posRescaler:advance(tCurr, {fIn}, {self.momPos}, false)
            self.solver:advance(tCurr, {self.momPos, em}, {fRhsOut})
         end
      else
         self.advSolver = function(tCurr, fIn, em, fRhsOut)
            self.solver:advance(tCurr, {fIn, em}, {fRhsOut})
         end
      end
   end
end

-- Nothing to calculate, just copy.
function IncompEulerSpecies:calcCouplingMomentsEvolve(tCurr, rkIdx)
   local fIn = self:rkStepperFields()[rkIdx]
   self.couplingMoments:copy(fIn)
end

function IncompEulerSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
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

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      local em = emIn[1]:rkStepperFields()[inIdx]
      self.advSolver(tCurr, fIn, em, fRhsOut)
   end

   for _, src in pairs(self.sources) do src:advance(tCurr, fIn, species, fRhsOut) end
end

function IncompEulerSpecies:fluidMoments() return { self.couplingMoments } end

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
