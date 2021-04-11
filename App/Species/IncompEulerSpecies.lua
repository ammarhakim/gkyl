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
local ConstDiffusionModDecl = require "Eq.constDiffusionData.ConstDiffusionModDecl"

local IncompEulerSpecies = Proto(FluidSpecies)

local SP_BC_ABSORB    = 1
local SP_BC_COPY      = 2
local SP_BC_ZEROFLUX  = 3
--local SP_BC_WALL      = 4  -- Not yet available for IncompEulerSpecies.
local SP_BC_DIRICHLET = 5    -- Specify the value (currently only for diffusion term).
local SP_BC_NEUMANN   = 6    -- Specify the derivative (currently only for diffusion term).
IncompEulerSpecies.bcAbsorb    = SP_BC_ABSORB       -- Absorb.
IncompEulerSpecies.bcCopy      = SP_BC_COPY         -- Copy skin-cell values.
IncompEulerSpecies.bcZeroFlux  = SP_BC_ZEROFLUX     -- Zero flux.
IncompEulerSpecies.bcDirichlet = SP_BC_DIRICHLET    -- Specify the value (currently only for diffusion term).
IncompEulerSpecies.bcNeumann   = SP_BC_NEUMANN      -- Specify the derivative (currently only for diffusion term).

function IncompEulerSpecies:fullInit(appTbl)
   IncompEulerSpecies.super.fullInit(self, appTbl)

   self.nMoments = 1
   self.zeroFluxDirections = {}
end

function IncompEulerSpecies:bcDirichletFunc(dir, tm, idxIn, fIn, fOut)
   -- Impose f=fBC at the boundary.
   if (idxIn[dir] == 1) then
      self.constDiffDirichletBCs[dir][1](self.grid:dx(dir),fIn:data(),self.auxBCvalues[dir][1],fOut:data())
   else
      self.constDiffDirichletBCs[dir][2](self.grid:dx(dir),fIn:data(),self.auxBCvalues[dir][2],fOut:data())
   end
end

function IncompEulerSpecies:bcNeumannFunc(dir, tm, idxIn, fIn, fOut)
   -- Impose f'=fpBC at the boundary.
   if (idxIn[dir] == 1) then
      self.constDiffNeumannBCs[dir][1](self.grid:dx(dir),fIn:data(),self.auxBCvalues[dir][1],fOut:data())
   else
      self.constDiffNeumannBCs[dir][2](self.grid:dx(dir),fIn:data(),self.auxBCvalues[dir][2],fOut:data())
   end
end


function IncompEulerSpecies:appendBoundaryConditions(dir, edge, bcType)
   if bcType == SP_BC_ABSORB then
      local function bcAbsorbFunc(...)    return self:bcAbsorbFunc(...) end
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, { bcAbsorbFunc }, "pointwise"))
   elseif bcType == SP_BC_COPY then
      local function bcCopyFunc(...)      return self:bcCopyFunc(...) end
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, { bcCopyFunc }, "pointwise"))
   elseif bcType == SP_BC_ZEROFLUX then
      table.insert(self.zeroFluxDirections, dir)
   elseif bcType[1] == SP_BC_DIRICHLET then
      self.constDiffDirichletBCs = ConstDiffusionModDecl.selectBCs(self.basis:id(), self.basis:ndim(), self.basis:polyOrder(), 2, "Dirichlet")
      local function bcDirichletFunc(...) return self:bcDirichletFunc(...) end
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, { bcDirichletFunc }, "pointwise"))
   elseif bcType[1] == SP_BC_NEUMANN then
      self.constDiffNeumannBCs = ConstDiffusionModDecl.selectBCs(self.basis:id(), self.basis:ndim(), self.basis:polyOrder(), 2, "Neumann")
      local function bcNeumannFunc(...)   return self:bcNeumannFunc(...) end
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, { bcNeumannFunc }, "pointwise"))
   else
      assert(false, "IncompEulerSpecies: Unsupported BC type!")
   end
end

function IncompEulerSpecies:createSolver(hasE, hasB)
   -- Run the FluidSpecies 'createSolver()' to initialize the
   -- collisions (diffusion) solver.
   IncompEulerSpecies.super.createSolver(self)

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

   if self.mSource and self.evolveSources then
      -- Add source term to the RHS.
      -- Barrier over shared communicator before accumulate.
      Mpi.Barrier(self.grid:commSet().sharedComm)
      fRhsOut:accumulate(self.sourceTimeDependence(tCurr), self.mSource)
   end
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
