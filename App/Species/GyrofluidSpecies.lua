-- Gkyl ------------------------------------------------------------------------
--
-- Gyrofluid species app.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto          = require "Lib.Proto"
local FluidSpecies   = require "App.Species.FluidSpecies"
local DiagsApp       = require "App.Species.Diagnostics.SpeciesDiagnostics"
local GyrofluidDiags = require "App.Species.Diagnostics.GyrofluidDiagnostics"
local Mpi            = require "Comm.Mpi"
local GyrofluidEq    = require "Eq.Gyrofluid"
local Updater        = require "Updater"
local DataStruct     = require "DataStruct"
local Time           = require "Lib.Time"
local Constants      = require "Lib.Constants"
local Lin            = require "Lib.Linalg"
local xsys           = require "xsys"
local lume           = require "Lib.lume"

local GyrofluidSpecies = Proto(FluidSpecies)

-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB = 1
local SP_BC_COPY   = 6
GyrofluidSpecies.bcAbsorb = SP_BC_ABSORB      -- Absorb all particles.
GyrofluidSpecies.bcCopy   = SP_BC_COPY        -- Copy stuff.

function GyrofluidSpecies:fullInit(appTbl)
   GyrofluidSpecies.super.fullInit(self, appTbl)

   self.n0 = self.tbl.n0 or n0
   assert(self.n0, "GyrofluidSpecies: must specify background density as global variable 'n0' in species table as 'n0 = ...'")
   self.kappaPar, self.kappaPerp = self.tbl.kappaPar, self.tbl.kappaPerp

   self.nMoments = 3+1
   self.zeroFluxDirections = {}
   self._firstMomentCalc = true  -- To avoid re-calculating moments when not evolving.
end

function GyrofluidSpecies:alloc(nRkDup)
   -- Allocate distribution function.
   GyrofluidSpecies.super.alloc(self, nRkDup)

   self.primMomSelf = self:allocVectorMoment(3)   -- upar, Tpar, Tperp of this species.

   -- Allocate fields to store separate stepped and primitive moments.
   self.mJacM0    = self:allocMoment()
   self.mJacM1    = self:allocMoment()
   self.mJacM2    = self:allocMoment()
   self.jacM2perp = self:allocMoment()
   self.uParSelf  = self:allocMoment()
   self.TparSelf  = self:allocMoment()
   self.TperpSelf = self:allocMoment()

   -- Intermediate fields needed in computing the primitive moments.
   self.pParJac    = self:allocMoment()
   self.pPerpJac   = self:allocMoment()
   self.mJacM2flow = self:allocMoment()

   -- For interfacing with gyrokinetics and GkField:
   self.jacM0 = self:allocMoment();  self.jacM0Aux = self:allocMoment()
   self.jacM1 = self:allocMoment();  self.jacM1Aux = self:allocMoment()
   self.jacM2 = self:allocMoment();  self.jacM2Aux = self:allocMoment()

   self.polarizationWeight = self:allocMoment()   -- Not used when using linearized poisson solve.

   self.first = true
end

function GyrofluidSpecies:createSolver(hasPhi, hasApar, externalField)
   -- Run the FluidSpecies 'createSolver()' to initialize the collisions solver.
   GyrofluidSpecies.super.createSolver(self,externalField)

   -- Set up Jacobian.
   if externalField then
      self.bmagFunc  = externalField.bmagFunc
      self.jacobFunc = externalField.jacobGeoFunc

      local xMid = {}
      for d = 1,self.ndim do xMid[d]=self.grid:mid(d) end
      self.B0 = externalField.bmagFunc(0.0, xMid)

      self.bmag     = assert(externalField.geo.bmag, "nil bmag")
      self.bmagInv  = externalField.geo.bmagInv
      self.jacob    = externalField.geo.jacobGeo
      self.jacobInv = externalField.geo.jacobGeoInv
   end

   -- Create updater to advance solution by one time-step.
   self.equation = GyrofluidEq {
      onGrid    = self.grid,
      basis     = self.basis,
      charge    = self.charge,
      mass      = self.mass,
      kappaPar  = self.kappaPar,
      kappaPerp = self.kappaPerp,
      bmagFunc  = self.bmagFunc,
   }

   self.solver = Updater.HyperDisCont {
      onGrid             = self.grid,
      basis              = self.basis,
      cfl                = self.cfl,
      equation           = self.equation,
      zeroFluxDirections = self.zeroFluxDirections,
      clearOut           = false,   -- Continue accumulating into output field.
      globalUpwind       = true,
   }

   -- Offsets needed to fetch specific moments from CartField
   -- containing the stepped moments (e.g. with :combineOffset).
   self.mJacM0Off    = 0 
   self.mJacM1Off    = 1*self.basis:numBasis() 
   self.mJacM2Off    = 2*self.basis:numBasis()  
   self.jacM2perpOff = 3*self.basis:numBasis() 
   -- Package them into a single table for easier access.
   self.momOff = {self.mJacM0Off,self.mJacM1Off,self.mJacM2Off,self.jacM2perpOff}

   self.timers = {couplingMom=0., weakMom=0., sources=0.}
end

function GyrofluidSpecies:getMomOff(momIdx)
   local offOut = momIdx and self.momOff[momIdx] or self.momOff
   return offOut
end

function GyrofluidSpecies:initCrossSpeciesCoupling(species)
   -- Nothing implemented yet.
end

function GyrofluidSpecies:uParCalc(tm, momIn, mJacM0, mJacM1, uParOut)
   -- Calculate the parallel flow speed.
   mJacM0:combineOffset(1, momIn, self.mJacM0Off)
   mJacM1:combineOffset(1, momIn, self.mJacM1Off)
   self.weakDivide:advance(tm, {mJacM0, mJacM1}, {uParOut})
end

function GyrofluidSpecies:pPerpJacCalc(tm, momIn, jacM2perp, pPerpJacOut)
   -- Compute the perpendicular pressure (times Jacobian).
   jacM2perp:combineOffset(1, momIn, self.jacM2perpOff)
   self.weakMultiply:advance(tm, {jacM2perp, self.bmag}, {pPerpJacOut})
end

function GyrofluidSpecies:pParJacCalc(tm, momIn, uPar, mJacM1, mJacM2flow, mJacM2, pPerpJac, pParJacOut)
   -- Compute the parallel pressure (times Jacobian).
   self.weakMultiply:advance(tm, {uPar, mJacM1}, {mJacM2flow})
   mJacM2:combineOffset(1, momIn, self.mJacM2Off)
   Mpi.Barrier(self.grid:commSet().sharedComm)   -- Barrier over sharedComm before combine.
   pParJacOut:combine(2., mJacM2, -2., pPerpJac, -1., mJacM2flow)
end

function GyrofluidSpecies:TparCalc(tm, mJacM0, pParJac, TparOut)
   -- Compute the parallel temperature.
   self.weakDivide:advance(tm, {mJacM0, pParJac}, {TparOut})
   TparOut:scale(self.mass)
end

function GyrofluidSpecies:TperpCalc(tm, mJacM0, pPerpJac, TperpOut)
   -- Compute the perpendicular temperature.
   self.weakDivide:advance(tm, {mJacM0, pPerpJac}, {TperpOut})
   TperpOut:scale(self.mass)
end

function GyrofluidSpecies:calcCouplingMoments(tCurr, rkIdx, species)
   local momIn = self:rkStepperFields()[rkIdx]
   -- Compute moments needed in coupling to fields and collisions.
   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()

      if self.deltaF then momIn:accumulate(-1.0, self.momBackground) end

      if not self.momentFlags[1] then -- No need to recompute if already computed.
         -- Calculate the parallel flow speed.
         self:uParCalc(tCurr, momIn, self.mJacM0, self.mJacM1, self.uParSelf)
         
         -- Get the perpendicular and parallel pressures (times Jacobian).
         self:pPerpJacCalc(tCurr, momIn, self.jacM2perp, self.pPerpJac)
         self:pParJacCalc(tCurr, momIn, self.uParSelf, self.mJacM1, self.mJacM2flow,
                                 self.mJacM2, self.pPerpJac, self.pParJac)

         -- Compute perpendicular and parallel temperatures.
         self:TparCalc(tCurr, self.mJacM0, self.pParJac, self.TparSelf)
         self:TperpCalc(tCurr, self.mJacM0, self.pPerpJac, self.TperpSelf)

         -- Package self primitive moments into a single field (expected by equation object).
         self.primMomSelf:combineOffset(1.,  self.uParSelf, 0*self.basis:numBasis(),
                                        1.,  self.TparSelf, 1*self.basis:numBasis(),
                                        1., self.TperpSelf, 2*self.basis:numBasis())

         -- Indicate that primitive moments have been computed.
         self.momentFlags[1] = true
      end

      if self.deltaF then momIn:accumulate(1.0, self.momBackground) end

      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   end
end

function GyrofluidSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self:setActiveRKidx(inIdx)
   self.tCurr = tCurr
   local momIn     = self:rkStepperFields()[inIdx]
   local momRhsOut = self:rkStepperFields()[outIdx]

   local em     = emIn[1]:rkStepperFields()[inIdx] -- Dynamic fields (e.g. phi, Apar)
   local emFunc = emIn[2]:rkStepperFields()[1]     -- Geometry/external field.

   if self.positivityRescale then
      self.posRescaler:advance(tCurr, {momIn}, {self.momPos}, false)
      momIn = self.momPos
   end

   -- Clear RHS, because HyperDisCont set up with clearOut = false.
   momRhsOut:clear(0.0)

   -- Perform the collision update. This includes terms that comes from
   -- collisions but also other objects in the Collisions App (e.g. diffusion).
   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
         c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
         c:advance(tCurr, momIn, species, momRhsOut)
      end
   end

   -- Complete the field solve.
   emIn[1]:phiSolve(tCurr, species, inIdx, outIdx)

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solver:advance(tCurr, {momIn, em, emFunc, self.primMomSelf}, {momRhsOut})
   else
      self.equation:setAuxFields({em, emFunc, self.primMomSelf})  -- Set auxFields in case they are needed by BCs/collisions.
   end

   -- Save boundary fluxes for diagnostics.
   if self.hasNonPeriodicBc and self.boundaryFluxDiagnostics then
      for _, bc in ipairs(self.boundaryConditions) do
         bc:storeBoundaryFlux(tCurr, outIdx, momRhsOut)
      end
   end

   if self.mSource and self.evolveSources then
      local tm = Time.clock()
      -- Add source term to the RHS.
      -- Barrier over shared communicator before accumulate.
      Mpi.Barrier(self.grid:commSet().sharedComm)
      momRhsOut:accumulate(self.sourceTimeDependence(tCurr), self.mSource)
      self.timers.sources = self.timers.sources + Time.clock() - tm
   end

end

function GyrofluidSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- Need to wrap member functions so that self is passed.
   local function bcAbsorbFunc(...)  return self:bcAbsorbFunc(...) end
   local function bcCopyFunc(...)    return self:bcCopyFunc(...) end

   if bcType == SP_BC_ABSORB then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, { bcAbsorbFunc }, "pointwise"))
   elseif bcType == SP_BC_COPY then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, { bcCopyFunc }, "pointwise"))
   else
      assert(false, "GyrofluidSpecies: Unsupported BC type!")
   end
end

function GyrofluidSpecies:createDiagnostics()  -- More sophisticated/extensive diagnostics go in Species/Diagnostics.
   -- Create this species' diagnostics.
   if self.tbl.diagnostics then
      self.diagnostics[self.name] = DiagsApp{}
      self.diagnostics[self.name]:fullInit(self, GyrofluidDiags)
   end

   -- Many diagnostics require dividing by the Jacobian (if present).
   -- Predefine the function that does that.
   self.calcNoJacMom = self.jacobInv
      and function(tm, rkIdx) self.weakMultiply:advance(tm, {self:getMoments(rkIdx), self.jacobInv}, {self.noJacMom}) end
      or function(tm, rkIdx) self.noJacMom:copy(self:getMoments(rkIdx)) end
end

function GyrofluidSpecies:fluidMoments()
   self.jacM0:combine(1./self.mass, self.mJacM0) 
   self.jacM1:combine(1./self.mass, self.mJacM1)
   self.jacM2:combine(1./self.mass, self.mJacM2)
   return { self.jacM0, self.jacM1, self.jacM2 }
end

function GyrofluidSpecies:selfPrimitiveMoments()
   return { self.uParSelf, self.TparSelf, self.TperpSelf }
end

function GyrofluidSpecies:crossPrimitiveMoments(otherSpeciesName)
end

function GyrofluidSpecies:getNumDensity(rkIdx)
   -- If no rkIdx specified, assume numDensity has already been calculated.
   if rkIdx == nil then
      self.jacM0:combine(1./self.mass, self.mJacM0) 
      return self.jacM0
   end

   local momIn = self:rkStepperFields()[rkIdx]

   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      self.jacM0Aux:combineOffset(1./self.mass, momIn, self.mJacM0Off)
      if self.deltaF then self.jacM0Aux:accumulateOffset(-1.0/self.mass, self.momBackground, self.mJacM0Off) end
      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end

   return self.jacM0Aux
end

function GyrofluidSpecies:getMomDensity(rkIdx)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then
      self.jacM1:combine(1./self.mass, self.mJacM1)
      return self.jacM1
   end

   local momIn = self:rkStepperFields()[rkIdx]

   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      self.jacM1Aux:combineOffset(1./self.mass, momIn, self.mJacM1Off)
      if self.deltaF then self.jacM1Aux:accumulate(-1.0/self.mass, self.momBackground, 1*self.basis:numBasis()) end
      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end

   return self.jacM1Aux
end

function GyrofluidSpecies:getPolarizationWeight(linearized)
   if linearized == false then
      self.weakMultiply:advance(0.0, {self.jacM0, self.bmagInv}, {self.polarizationWeight})
      self.weakMultiply:advance(0.0, {self.polarizationWeight, self.bmagInv}, {self.polarizationWeight})
      self.polarizationWeight:scale(self.mass)
      return self.polarizationWeight
   else
      return self.n0*self.mass/self.B0^2
   end
end

function GyrofluidSpecies:momCalcTime()
   local tm = self.timers.couplingMom
   for _, dOb in pairs(self.diagnostics) do
      tm = tm + dOb:getDiagTime()
   end
   return tm
end
function GyrofluidSpecies:solverVolTime()
   return self.equation.totalVolTime
end
function GyrofluidSpecies:solverSurfTime()
   return self.equation.totalSurfTime
end
function GyrofluidSpecies:totalSolverTime()
   local timer = self.solver.totalTime
   if self.solverStep2 then timer = timer + self.solverStep2.totalTime end
   if self.solverStep3 then timer = timer + self.solverStep3.totalTime end
   if self.posRescaler then timer = timer + self.posRescaler.totalTime end
   return timer
end

return GyrofluidSpecies
