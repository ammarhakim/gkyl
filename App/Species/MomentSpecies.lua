-- Gkyl ------------------------------------------------------------------------
--
-- Species object constructed from moment equations.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct    = require "DataStruct"
local Lin           = require "Lib.Linalg"
local LinearTrigger = require "Lib.LinearTrigger"
local Mpi           = require "Comm.Mpi"
local Proto         = require "Lib.Proto"
local FluidSpecies  = require "App.Species.FluidSpecies"
local Time          = require "Lib.Time"
local Updater       = require "Updater"
local xsys          = require "xsys"
local ffi           = require "ffi"
local Euler         = require "Eq.Euler"
local TenMoment     = require "Eq.TenMoment"

-- Species object treated as moment equations.
local MomentSpecies = Proto(FluidSpecies)

-- Add constants to object indicate various supported boundary conditions.
local SP_BC_OPEN      = 1
local SP_BC_COPY      = SP_BC_OPEN
local SP_BC_WALL      = 4
-- The following two are not yet available for the MomentSpecies.
-- local SP_BC_DIRICHLET = 5    -- Specify the value (currently only for diffusion term).
-- local SP_BC_NEUMANN   = 6    -- Specify the derivative (currently only for diffusion term).

MomentSpecies.bcOpen = SP_BC_OPEN    -- Open BCs.
MomentSpecies.bcCopy = SP_BC_COPY    -- Copy BCs.
MomentSpecies.bcWall = SP_BC_WALL    -- Wall BCs.

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function MomentSpecies:fullInit(appTbl)
   MomentSpecies.super.fullInit(self, appTbl)

   self.equation = self.tbl.equation -- Equation system to evolve.
   self.nMoments = self.tbl.equation:numEquations()
   self.nGhost = 2     -- We need two ghost-cells.

   self.limiter   = self.tbl.limiter and self.tbl.limiter or "monotonized-centered"
   self.hyperSlvr = {} -- List of solvers.

   -- Invariant (positivity-preserving) equation system to evolve.
   self.equationInv  = self.tbl.equationInv
   self.hyperSlvrInv = {} -- List of solvers.
   self.limiterInv   = self.tbl.limiterInv and self.tbl.limiterInv or "zero"
   -- Always use invariant eqn.
   self.forceInv = self.tbl.forceInv and (self.equationInv ~= nil)
   -- Use invariant eqn. in next step; could change during run.
   self.tryInv = false

   self._myIsInv = ffi.new("int[2]")
   self._isInv   = ffi.new("int[2]")

   self._hasSsBnd  = xsys.pickBool(self.tbl.hasSsBnd, false)
   self._inOutFunc = self.tbl.inOutFunc
end

function MomentSpecies:createSolver(hasE, hasB)
   if self._hasSsBnd then
      self._inOut = DataStruct.Field {
         onGrid        = self.grid,
         numComponents = 1,
         ghost         = {2, 2}
      }
      local project = Updater.ProjectOnBasis {
         onGrid          = self.grid,
         basis           = self.basis,
         evaluate        = self._inOutFunc,
         projectOnGhosts = true,
      }
      project:advance(0.0, {}, {self._inOut})
      self.momIo:write(self._inOut, string.format("%s_inOut.bp", self.name), 0, 0)
   end

   local ndim = self.grid:ndim()
   for d = 1, ndim do
      self.hyperSlvr[d] = Updater.WavePropagation {
         onGrid           = self.grid,
         equation         = self.equation,
         limiter          = self.limiter,
         cfl              = self.cfl,
         updateDirections = {d},
         hasSsBnd         = self._hasSsBnd,
         inOut            = self._inOut,
      }
      if self.equationInv ~= nil then
         self.hyperSlvrInv[d] = Updater.WavePropagation {
            onGrid           = self.grid,
            equation         = self.equationInv,
            limiter          = self.limiterInv,
            cfl              = self.cfl,
            updateDirections = {d},
            hasSsBnd         = self._hasSsBnd,
            inOut            = self._inOut
         }
      end
   end

   -- This perhaps should go to createBCs but at that point _inOut is not
   -- created yet.
   if (self._hasSsBnd) then
      local function handleSsBc(dir, bcList)
         for _, bc in ipairs(bcList) do
            if bc then
               self:appendSsBoundaryConditions(dir, self._inOut, bc)
            end
         end
      end

      for d = 1, ndim do
         handleSsBc(d, self.ssBc)
      end
   end

end

function MomentSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   -- Does nothing: perhaps when DG is supported this will need to be
   -- modified.

   return true, GKYL_MAX_DOUBLE
end

function MomentSpecies:appendBoundaryConditions(dir, edge, bcType)
   local function bcCopyFunc(...) return self:bcCopyFunc(...) end

   if bcType == SP_BC_COPY then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, edge, { {bcCopyFunc} }))
   elseif bcType == SP_BC_WALL then
      -- FIXME better to define and use self.equation.bcWall.
      local bcWall
      if self.nMoments == 5 then
         bcWall = Euler.bcWall
      elseif self.nMoments == 10 then
         bcWall = TenMoment.bcWall
      else
         assert(false, "MomentSpecies: bcWall not provided by the equation!")
      end
      table.insert(self.boundaryConditions,
                   self:makeBcUpdater(dir, edge, bcWall))
   elseif type(bcType) == "table" then
      -- bcType can be literally a list of functions.
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, edge, bcType ))
   else
      assert(false, "MomentSpecies: Unsupported BC type!")
   end
end

-- TODO: merge into appendBoundaryConditions.
function MomentSpecies:appendSsBoundaryConditions(dir, inOut, bcType)
   local function bcCopyFunc(...) return self:bcCopyFunc(...) end

   if bcType == SP_BC_COPY then
      table.insert(self.ssBoundaryConditions,
		   self:makeSsBcUpdater(dir, inOut, { {bcCopyFunc} }))
   elseif bcType == SP_BC_WALL then
      -- FIXME better to define and use self.equation.bcWall.
      local bcWall
      if self.nMoments == 5 then
        bcWall = Euler.bcWall
      elseif self.nMoments == 10 then
        bcWall = TenMoment.bcWall
      else
        assert(false, "MomentSpecies: bcWall not provided by the equation!")
      end
      table.insert(self.ssBoundaryConditions,
		   self:makeSsBcUpdater(dir, inOut, bcWall))
   elseif type(bcType) == "table" then
      -- bcType can be literally a list of functions.
      table.insert(self.ssBoundaryConditions,
		   self:makeSsBcUpdater(dir, inOut, bcType ))
   else
      assert(false, "MomentSpecies: Unsupported BC type!")
   end
end

function MomentSpecies:updateInDirection(dir, tCurr, dt, fIn, fOut, tryInv)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   local tryInv_next         = false
   if self.evolve then
      self:applyBc(tCurr, fIn, dir)
      assert(self:checkInv(fIn))
      if self.forceInv or tryInv then
         self.hyperSlvrInv[dir]:setDtAndCflRate(dt, nil)
         status, dtSuggested = self.hyperSlvrInv[dir]:advance(tCurr, {fIn}, {fOut})
         -- If result is OK, do not try to use invariant eqn. in next step.
         tryInv_next = not status
         if status and (not self:checkInv(fOut)) then
            assert(false, "** Invalid output using Lax flux!")
         end
      else
         self.hyperSlvr[dir]:setDtAndCflRate(dt, nil)
         status, dtSuggested = self.hyperSlvr[dir]:advance(tCurr, {fIn}, {fOut})
         tryInv_next = status and not self:checkInv(fOut)
      end
   else
      fOut:copy(fIn)
   end
   return status, dtSuggested, tryInv_next
end

function MomentSpecies:totalSolverTime()
   local tm = 0.0
   for d = 1, self.grid:ndim() do
      tm = tm + self.hyperSlvr[d].totalTime
   end
   return tm
end

function MomentSpecies:checkInv(fIn)
   local fInIndexer = fIn:genIndexer()
   local fInPtr     = fIn:get(1)

   local isInv      = true
   local localRange = fIn:localRange()   
   for idx in localRange:rowMajorIter() do
      self.grid:setIndex(idx)

      fIn:fill(fInIndexer(idx), fInPtr)
      if not self.equation:isPositive(fInPtr) then
        isInv = false
        break
      end
  end

  self._myIsInv[0] = isInv and 1 or 0
  Mpi.Allreduce(self._myIsInv, self._isInv, 1, Mpi.INT, Mpi.LAND, self.grid:commSet().comm)
  return self._isInv[0] == 1 and true or false
end

return MomentSpecies
