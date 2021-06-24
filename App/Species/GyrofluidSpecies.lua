-- Gkyl ------------------------------------------------------------------------
--
-- Gyrofluid species app.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto          = require "Lib.Proto"
local FluidSpecies   = require "App.Species.FluidSpecies"
local DiagsApp       = require "App.Diagnostics.SpeciesDiagnostics"
local GyrofluidDiags = require "App.Diagnostics.GyrofluidDiagnostics"
local Mpi            = require "Comm.Mpi"
local GyrofluidEq    = require "Eq.Gyrofluid"
local Updater        = require "Updater"
local Time           = require "Lib.Time"
local Constants      = require "Lib.Constants"
local BasicBC        = require ("App.BCs.GyrofluidBasic").GyrofluidBasic
local Lin            = require "Lib.Linalg"
local xsys           = require "xsys"
local lume           = require "Lib.lume"

local GyrofluidSpecies = Proto(FluidSpecies)

-- ............. Backwards compatible treatment of BCs .....................-- 
-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB   = 1
local SP_BC_ZEROFLUX = 5
local SP_BC_COPY     = 6
GyrofluidSpecies.bcAbsorb   = SP_BC_ABSORB     -- Absorb all particles.
GyrofluidSpecies.bcCopy     = SP_BC_COPY       -- Copy stuff.
GyrofluidSpecies.bcZeroFlux = SP_BC_ZEROFLUX   -- Zero flux.

function GyrofluidSpecies:makeBcApp(bcIn)
   local bcOut
   if bcIn == SP_BC_COPY then
      print("GyrofluidSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="copy", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_ABSORB then
      print("GyrofluidSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="absorb", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_ZEROFLUX or bcIn.tbl.kind=="zeroFlux" then
      bcOut = "zeroFlux"
      table.insert(self.zeroFluxDirections, dir)
   end
   return bcOut
end

-- ............. End of backwards compatibility for BCs .....................-- 

function GyrofluidSpecies:fullInit(appTbl)
   GyrofluidSpecies.super.fullInit(self, appTbl)

   self.n0 = self.tbl.n0 or n0
   assert(self.n0, "GyrofluidSpecies: must specify background density as global variable 'n0' in species table as 'n0 = ...'")
   self.kappaPar, self.kappaPerp = self.tbl.kappaPar, self.tbl.kappaPerp

   self.nMoments = 3+1
end

function GyrofluidSpecies:alloc(nRkDup)
   GyrofluidSpecies.super.alloc(self, nRkDup)   -- Call the FluidSpecies :alloc method.

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

   -- Offsets needed to fetch specific moments from CartField
   -- containing the stepped moments (e.g. with :combineOffset).
   self.mJacM0Off    = 0 
   self.mJacM1Off    = 1*self.basis:numBasis() 
   self.mJacM2Off    = 2*self.basis:numBasis()  
   self.jacM2perpOff = 3*self.basis:numBasis() 
   -- Package them into a single table for easier access.
   self.momOff = {self.mJacM0Off,self.mJacM1Off,self.mJacM2Off,self.jacM2perpOff}

   self.first = true
end

function GyrofluidSpecies:createSolver(field, externalField)

   self.cSound = self:allocMoment()  -- Sound speed. Do this before calling FluidSpecies' :createSolver.

   -- Run the FluidSpecies 'createSolver()' to initialize the collisions solver.
   GyrofluidSpecies.super.createSolver(self, field, externalField)

   local hasE, hasB       = field:hasEB()
   local extHasE, extHasB = externalField:hasEB()

   local hasPhi  = hasE or extHasE
   local hasApar = hasB or extHasB

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
      onGrid    = self.grid,    kappaPar  = self.kappaPar,
      basis     = self.basis,   kappaPerp = self.kappaPerp,
      charge    = self.charge,  bmagFunc  = self.bmagFunc,
      mass      = self.mass,
   }

   self.solver = Updater.HyperDisCont {
      onGrid   = self.grid,      zeroFluxDirections = self.zeroFluxDirections,
      basis    = self.basis,     clearOut     = false,   -- Continue accumulating into output field.
      cfl      = self.cfl,       globalUpwind = true,
      equation = self.equation,
   }

   self.sqrtOnBasis = Updater.SqrtOnBasis{onGrid = self.grid,  basis = self.basis, onGhosts=true}

   if self.positivityRescale then
      self.returnPosRescaledMom = function(tCurr, momIn)
         self.posRescaler:advance(tCurr, {momIn}, {self.momPos}, false)
         return self.momPos
      end
   else
      self.returnPosRescaledMom = function(tCurr, momIn) return momIn end
   end

   if self.deltaF then
      self.minusBackgroundJacM0 = function(jacM0In)
         jacM0In:accumulateOffset(-1.0/self.mass, self.momBackground, self.mJacM0Off)
      end
      self.minusBackgroundJacM1 = function(jacM1In)
         jacM1In:accumulateOffset(-1.0/self.mass, self.momBackground, self.mJacM1Off)
      end
   else
      self.minusBackgroundJacM0 = function(jacM0In) end
      self.minusBackgroundJacM1 = function(jacM1In) end
   end

   self.timers = {couplingMom=0., weakMom=0., sources=0.}
end

function GyrofluidSpecies:getMomOff(momIdx)
   local offOut = momIdx and self.momOff[momIdx] or self.momOff
   return offOut
end

function GyrofluidSpecies:initCrossSpeciesCoupling(species)
   -- Save pointers electron and ion names. Will need them to compute sound speed.
   for nm, s in pairs(species) do
      if s.charge > 0. then
         self.ionNm = nm
      elseif s.charge < 0. then
         self.elcNm = nm
      end
   end
   self.ionMass = species[self.ionNm].mass
   -- Electrons may need to know the latest ion upar for sheath BCs.
   self.ionPrimMomSelf = species[self.ionNm].primMomSelf
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

function GyrofluidSpecies:pParJacCalc(tm, momIn, uPar, mJacM1, mJacM2flow, 
                                      mJacM2, pPerpJac, pParJacOut)
   -- Compute the parallel pressure (times Jacobian).
   self.weakMultiply:advance(tm, {uPar, mJacM1}, {mJacM2flow})
   mJacM2:combineOffset(1, momIn, self.mJacM2Off)
   Mpi.Barrier(self.grid:commSet().sharedComm)   -- Barrier over sharedComm before combine.
   pParJacOut:combine(2., mJacM2, -2., pPerpJac, -1., mJacM2flow)

   -- If the perpendicular pressure and/or the parallel flow energy density are sufficiently
   -- large, or if the kinetic energy density is sufficiently low, the parallel pressure
   -- could come out negative.
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

function GyrofluidSpecies:calcCouplingMomentsNoEvolve(tCurr, rkIdx, species) end
function GyrofluidSpecies:calcCouplingMomentsEvolve(tCurr, rkIdx, species)
   -- Compute moments needed in coupling to fields and collisions.
   if not self.momentFlags[1] then -- No need to recompute if already computed.
      local momIn   = self:rkStepperFields()[rkIdx]
      local tmStart = Time.clock()

      momIn = self.getMom_or_deltaMom(momIn)  -- Return full-F moments, or compute and return fluctuations.

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

      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   end
end

function GyrofluidSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self:setActiveRKidx(inIdx)
   self.tCurr = tCurr
   local momIn     = self:rkStepperFields()[inIdx]
   local momRhsOut = self:rkStepperFields()[outIdx]

   local em     = emIn[1]:rkStepperFields()[inIdx] -- Dynamic fields (e.g. phi, Apar)
   local extGeo = emIn[2]:rkStepperFields()[1]     -- Geometry/external field.

   momIn = self.returnPosRescaledMom(tCurr, momIn)

   -- Compute the sound speed sqrt((3*Ti_par+Te))/m_i.
   self.cSound:combineOffset(     3./self.ionMass, species[self.ionNm].primMomSelf, 1*self.basis:numBasis(),
                             1./(3.*self.ionMass), species[self.elcNm].primMomSelf, 1*self.basis:numBasis(),
                             2./(3.*self.ionMass), species[self.elcNm].primMomSelf, 2*self.basis:numBasis())
   self.sqrtOnBasis:advance(tCurr, {self.cSound}, {self.cSound})

   -- Clear RHS, because HyperDisCont set up with clearOut = false.
   momRhsOut:clear(0.0)

   -- Perform the collision update. This includes terms that comes from
   -- collisions but also other objects in the Collisions App (e.g. diffusion).
   for _, c in pairs(self.collisions) do
      c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      c:advance(tCurr, momIn, species, {em, extGeo}, momRhsOut)
   end

   -- Complete the field solve.
   emIn[1]:phiSolve(tCurr, species, inIdx, outIdx)

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solver:advance(tCurr, {momIn, em, extGeo, self.primMomSelf, self.cSound}, {momRhsOut})
   else
      self.equation:setAuxFields({em, extGeo, self.primMomSelf, self.cSound})  -- Set auxFields in case they are needed by BCs/collisions.
   end

   for _, bc in pairs(self.nonPeriodicBCs) do
      bc:storeBoundaryFlux(tCurr, outIdx, momRhsOut)   -- Save boundary fluxes.
   end

   for _, src in lume.orderedIter(self.sources) do src:advance(tCurr, momIn, species, momRhsOut) end
end

function GyrofluidSpecies:createDiagnostics(field)  -- More sophisticated/extensive diagnostics go in Species/Diagnostics.
   -- Create this species' diagnostics.
   if self.tbl.diagnostics then
      self.diagnostics[self.name] = DiagsApp{implementation = GyrofluidDiags()}
      self.diagnostics[self.name]:fullInit(self, field, self)
   end

   for srcNm, src in lume.orderedIter(self.sources) do
      self.diagnostics[self.name..srcNm] = src:createDiagnostics(self, field)
   end

   for bcNm, bc in lume.orderedIter(self.nonPeriodicBCs) do
      self.diagnostics[self.name..bcNm] = bc:createDiagnostics(self, field)
   end
   lume.setOrder(self.diagnostics)

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

function GyrofluidSpecies:selfPrimitiveMoments()
   return { self.uParSelf, self.TparSelf, self.TperpSelf }
end
function GyrofluidSpecies:getPressuresJac()
   return { self.pParJac, self.pPerpJac }
end

function GyrofluidSpecies:getNumDensity(rkIdx)
   -- If no rkIdx specified, assume numDensity has already been calculated.
   if rkIdx == nil then
      self.jacM0:combine(1./self.mass, self.mJacM0) 
      return self.jacM0
   end

   local momIn = self:rkStepperFields()[rkIdx]

   local tmStart = Time.clock()
   self.jacM0Aux:combineOffset(1./self.mass, momIn, self.mJacM0Off)
   self.minusBackgroundJacM0(self.jacM0Aux)
   self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart

   return self.jacM0Aux
end

function GyrofluidSpecies:getMomDensity(rkIdx)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then
      self.jacM1:combine(1./self.mass, self.mJacM1)
      return self.jacM1
   end

   local momIn = self:rkStepperFields()[rkIdx]

   local tmStart = Time.clock()
   self.jacM1Aux:combineOffset(1./self.mass, momIn, self.mJacM1Off)
   self.minusBackgroundJacM1(self.jacM1Aux)
   self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart

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
   for _, dOb in lume.orderedIter(self.diagnostics) do
      tm = tm + dOb:getDiagTime()
   end
   return tm
end
function GyrofluidSpecies:solverVolTime() return self.equation.totalVolTime end
function GyrofluidSpecies:solverSurfTime() return self.equation.totalSurfTime end
function GyrofluidSpecies:totalSolverTime()
   local timer = self.solver.totalTime
   if self.posRescaler then timer = timer + self.posRescaler.totalTime end
   return timer
end

return GyrofluidSpecies
