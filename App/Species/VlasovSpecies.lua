-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov species
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local KineticSpecies = require "App.Species.KineticSpecies"
local Mpi = require "Comm.Mpi"
local VlasovEq = require "Eq.Vlasov"
local Updater = require "Updater"
local DataStruct = require "DataStruct"
local Time = require "Lib.Time"

local VlasovSpecies = Proto(KineticSpecies)

-- add constants to object indicate various supported boundary conditions
local SP_BC_ABSORB = 1
local SP_BC_REFLECT = 3
local SP_BC_EXTERN = 4
local SP_BC_COPY = 5
-- AHH: This was 2 but seems that is unstable. So using plain copy
local SP_BC_OPEN = SP_BC_COPY
local SP_BC_ZEROFLUX = 6
local SP_BC_RESERVOIR = 7

VlasovSpecies.bcAbsorb = SP_BC_ABSORB -- absorb all particles
VlasovSpecies.bcOpen = SP_BC_OPEN -- zero gradient
VlasovSpecies.bcCopy = SP_BC_COPY -- copy stuff
VlasovSpecies.bcReflect = SP_BC_REFLECT -- specular reflection
VlasovSpecies.bcExternal = SP_BC_EXTERN -- load external BC file
VlasovSpecies.bcZeroFlux = SP_BC_ZEROFLUX
VlasovSpecies.bcReservoir = SP_BC_RESERVOIR

function VlasovSpecies:alloc(nRkDup)
   -- allocate distribution function
   VlasovSpecies.super.alloc(self, nRkDup)

   -- allocate fields to store coupling moments (for use in coupling
   -- to field and collisions)
   self.numDensity = self:allocMoment()
   self.momDensity = self:allocVectorMoment(self.vdim)
   self.ptclEnergy = self:allocMoment()

   -- allocate field to accumulate funcField if any
   self.totalEmField = self:allocVectorMoment(8) -- 8 components of EM field

   -- allocate field for external forces if any
   self.vExtForce = self:allocVectorMoment(self.vdim)

   -- allocate moment array for integrated moments (n, n*u_i, sum_i n*u_i^2, sum_i n*T_ii)
   self.flow = self:allocVectorMoment(self.vdim)
   self.kineticEnergyDensity = self:allocMoment()
   self.thermalEnergyDensity = self:allocMoment()

end

function VlasovSpecies:fullInit(appTbl)
   VlasovSpecies.super.fullInit(self, appTbl)

   local tbl = self.tbl
   -- if there is an external force, get the force function
   if tbl.vlasovExtForceFunc then
      self.vlasovExtForceFunc = tbl.vlasovExtForceFunc
   end

   local externalBC = tbl.externalBC
   if externalBC then
      self.wallFunction = require(externalBC)
   end
end

function VlasovSpecies:allocMomCouplingFields()
   return { currentDensity = self:allocVectorMoment(self.vdim) }
end


function VlasovSpecies:createSolver(hasE, hasB)
   -- run the KineticSpecies 'createSolver()' to initialize the
   -- collisions solver
   VlasovSpecies.super.createSolver(self)

   -- create updater to advance solution by one time-step
   local vlasovEqn = VlasovEq {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      charge = self.charge,
      mass = self.mass,
      hasElectricField = hasE,
      hasMagneticField = hasB,
   }

   -- must apply zero-flux BCs in velocity directions
   local zfd = { }
   for d = 1, self.vdim do 
     table.insert(self.zeroFluxDirections, self.cdim+d)
   end

   self.solver = Updater.HyperDisCont {
      onGrid = self.grid,
      basis = self.basis,
      cfl = self.cfl,
      equation = vlasovEqn,
      zeroFluxDirections = self.zeroFluxDirections,
   }

   -- create updaters to compute various moments
   self.numDensityCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "M0",
   }
   self.momDensityCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "M1i",
   }
   self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "M2",
   }
   -- create updater to compute M0, M1i, M2 moments sequentially
   -- this is used in calcCouplingMoments to reduce overhead and multiplications
   self.fiveMomentsCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "FiveMoments",
   }

   -- Updaters for the primitive moments
   -- These will be used to compute n*u^2 and n*T for computing integrated moments
   self.confDiv = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
   }
   self.confDotProduct = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "DotProduct",
   }

   if self.vlasovExtForceFunc then
      self.evalVlasovExtForce = Updater.ProjectOnBasis {
         onGrid = self.confGrid,
         basis = self.confBasis,
         evaluate = self.vlasovExtForceFunc,
         projectOnGhosts = false
      }
   end

   self.tmCouplingMom = 0.0 -- for timer 
end

function VlasovSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   local fIn = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   -- accumulate functional Maxwell fields (if needed)
   local emField = emIn[1]:rkStepperFields()[inIdx]
   local emFuncField = emIn[2]:rkStepperFields()[1]
   local totalEmField = self.totalEmField
   totalEmField:clear(0.0)

   local qbym = self.charge/self.mass

   if emField then totalEmField:accumulate(qbym, emField) end
   if emFuncField then totalEmField:accumulate(qbym, emFuncField) end
   
   -- if external force present (gravity, body force, etc.) accumulate it to electric field
   if self.vlasovExtForceFunc then
      local vExtForce = self.vExtForce
      self.evalVlasovExtForce:advance(tCurr, {}, {vExtForce})

      -- need to barrier over the shared communicator before accumulating force onto electric field
      Mpi.Barrier(self.confGrid:commSet().sharedComm)

      -- analogous to the current, the external force only gets accumulated onto the electric field
      local vItr, eItr = vExtForce:get(1), totalEmField:get(1)
      local vIdxr, eIdxr = vExtForce:genIndexer(), totalEmField:genIndexer()

      for idx in totalEmField:localRangeIter() do
         vExtForce:fill(vIdxr(idx), vItr)
         totalEmField:fill(eIdxr(idx), eItr)
         for i = 1, vExtForce:numComponents() do
            eItr[i] = eItr[i]+vItr[i]
         end
      end
   end

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solver:advance(tCurr, {fIn, totalEmField}, {fRhsOut})
   else
      fRhsOut:clear(0.0) -- no RHS
   end
   -- perform the collision update
   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
         c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
         c:advance(tCurr, fIn, species, fRhsOut)
         -- the full 'species' list is needed for the cross-species
         -- collisions
      end
   end

   if self.fSource and self.evolveSources then
      -- add source it to the RHS
      fRhsOut:accumulate(self.sourceTimeDependence(tCurr), self.fSource)
   end
end

function VlasovSpecies:createDiagnostics()
   -- create updater to compute volume-integrated moments
   -- function to check if integrated moment name is correct
   local function isIntegratedMomentNameGood(nm)
      if nm == "intM0" or nm == "intM1i" or nm == "intM2Flow" or nm == "intM2Thermal" or nm == "intL2" then
         return true
      end
      return false
   end

   local numCompInt = {}
   numCompInt["intM0"] = 1
   numCompInt["intM1i"] = self.vdim
   numCompInt["intM2Flow"] = 1
   numCompInt["intM2Thermal"] = 1
   numCompInt["intL2"] = 1

   self.diagnosticIntegratedMomentFields = { }
   self.diagnosticIntegratedMomentUpdaters = { } 
   -- allocate space to store moments and create moment updater
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
      if isIntegratedMomentNameGood(mom) then
         self.diagnosticIntegratedMomentFields[mom] = DataStruct.DynVector {
            numComponents = numCompInt[mom],
         }
         if mom == "intL2" then
            self.diagnosticIntegratedMomentUpdaters[mom] = Updater.CartFieldIntegratedQuantCalc {
               onGrid = self.grid,
               basis = self.basis,
               numComponents = numCompInt[mom],
               quantity = "V2"
            }
         else
            self.diagnosticIntegratedMomentUpdaters[mom] = Updater.CartFieldIntegratedQuantCalc {
               onGrid = self.confGrid,
               basis = self.confBasis,
               numComponents = numCompInt[mom],
               quantity = "V"
            }
         end
      else
         assert(false, string.format("Integrated Moment %s not valid", mom))
      end
   end

   -- function to check if moment name is correct
   local function isMomentNameGood(nm)
      return Updater.DistFuncMomentCalc:isMomentNameGood(nm)
   end

   local numComp = {}
   numComp["M0"] = 1
   numComp["M1i"] = self.vdim
   numComp["M2ij"] = self.vdim*(self.vdim+1)/2
   numComp["M2"] = 1
   numComp["M3i"] = self.vdim

   self.diagnosticMomentFields = { }
   self.diagnosticMomentUpdaters = { } 
   -- allocate space to store moments and create moment updater
   for i, mom in ipairs(self.diagnosticMoments) do
      if isMomentNameGood(mom) then
         self.diagnosticMomentFields[mom] = DataStruct.Field {
            onGrid = self.confGrid,
            numComponents = self.confBasis:numBasis()*numComp[mom],
            ghost = {1, 1}
         }
         self.diagnosticMomentUpdaters[mom] = Updater.DistFuncMomentCalc {
            onGrid = self.grid,
            phaseBasis = self.basis,
            confBasis = self.confBasis,
            moment = mom,
         }
      else
         assert(false, string.format("Moment %s not valid", mom))
      end
   end

   self.diagnosticWeakMoments = { }
   self.diagnosticAuxMoments = { }
end

-- BC functions
function VlasovSpecies:bcReflectFunc(dir, tm, idxIn, fIn, fOut)
   -- requires skinLoop = "flip"
   self.basis:flipSign(dir, fIn, fOut)
   self.basis:flipSign(dir+self.cdim, fOut, fOut)
end

function VlasovSpecies:bcExternFunc(dir, tm, idxIn, fIn, fOut)
   -- requires skinLoop = "flip"
   local velIdx = {}
   for d = 1, self.vdim do
      velIdx[d] = idxIn[self.cdim + d]
   end
   self.wallFunction[1](velIdx, fIn, fOut)
end

function VlasovSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- need to wrap member functions so that self is passed
   local function bcAbsorbFunc(...) return self:bcAbsorbFunc(...) end
   local function bcCopyFunc(...) return self:bcCopyFunc(...) end
   local function bcOpenFunc(...) return self:bcOpenFunc(...) end
   local function bcReflectFunc(...) return self:bcReflectFunc(...) end
   local function bcExternFunc(...) return self:bcExternFunc(...) end

   local vdir = dir + self.cdim

   if bcType == SP_BC_ABSORB then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcAbsorbFunc }, "pointwise", false))
   elseif bcType == SP_BC_OPEN then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcCopyFunc }, "pointwise", false))
   elseif bcType == SP_BC_COPY then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcCopyFunc }, "pointwise", false))
   elseif bcType == SP_BC_REFLECT then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcReflectFunc }, "flip", false))
   elseif bcType == SP_BC_EXTERN then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcExternFunc }, "flip", false))
   elseif bcType == SP_BC_ZEROFLUX then
      table.insert(self.zeroFluxDirections, dir)
   elseif bcType == SP_BC_RESERVOIR then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcCopyFunc }, "pointwise", true))
   else
      assert(false, "VlasovSpecies: Unsupported BC type!")
   end
end

function VlasovSpecies:calcCouplingMoments(tCurr, rkIdx)

   local tmStart = Time.clock()
   -- compute moments needed in coupling to fields and collisions
   local fIn = self:rkStepperFields()[rkIdx]
   if self.collisions then 
      self.fiveMomentsCalc:advance(tCurr, {fIn}, { self.numDensity, self.momDensity, self.ptclEnergy })
   else
      self.momDensityCalc:advance(tCurr, {fIn}, { self.momDensity })
   end
   self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart

end

-- function to compute n, u, nu^2, and nT for use in integrated moment routine
function VlasovSpecies:calcDiagnosticIntegratedMoments(tCurr)
   -- first compute M0, M1i, M2
   local fIn = self:rkStepperFields()[1]
   self.fiveMomentsCalc:advance(tCurr, {fIn}, { self.numDensity, self.momDensity, self.ptclEnergy })

   -- compute n*u^2 from n*u and n
   self.confDiv:advance(0., {self.numDensity, self.momDensity}, {self.flow})
   self.confDotProduct:advance(0., {self.flow, self.momDensity}, {self.kineticEnergyDensity})

   -- compute VDIM*n*T from M2 and kinetic energy density
   self.thermalEnergyDensity:combine(1.0, self.ptclEnergy, -1.0, self.kineticEnergyDensity)

   for i, mom in pairs(self.diagnosticIntegratedMoments) do
      if mom == "intM0" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.numDensity}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intM1i" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.momDensity}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intM2Flow" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.kineticEnergyDensity}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intM2Thermal" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.thermalEnergyDensity}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intL2" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.distf[1]}, {self.diagnosticIntegratedMomentFields[mom]})
      end
   end
end

function VlasovSpecies:fluidMoments()
   return { self.numDensity, self.momDensity, self.ptclEnergy }
end

function VlasovSpecies:getNumDensity(rkIdx)
   -- if no rkIdx specified, assume numDensity has already been calculated
   if rkIdx == nil then return self.numDensity end 

   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]
   self.numDensityCalc:advance(nil, {fIn}, { self.numDensity })

   self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart

   return self.numDensity
end

function VlasovSpecies:getMomDensity(rkIdx)
   -- if no rkIdx specified, assume momDensity has already been calculated
   if rkIdx == nil then return self.momDensity end 

   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]
   self.momDensityCalc:advance(nil, {fIn}, { self.momDensity })

   self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart

   return self.momDensity
end

function VlasovSpecies:momCalcTime()
   local tm = self.tmCouplingMom
   for i, mom in ipairs(self.diagnosticMoments) do
      tm = tm + self.diagnosticMomentUpdaters[mom].totalTime
   end
   return tm
end

-- please test this for higher than 1x1v... 
function VlasovSpecies:Maxwellian(xn, n0, T0, vdnIn)
   local vdn = vdnIn or {0, 0, 0}
   local vt2 = T0/self.mass
   local v2 = 0.0
   for d = self.cdim+1, self.cdim+self.vdim do
     v2 = v2 + (xn[d] - vdn[d-self.cdim])^2
   end
   return n0 / math.sqrt(2*math.pi*vt2)^self.vdim * math.exp(-v2/(2*vt2))
end

return VlasovSpecies
