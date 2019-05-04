-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov species
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto          = require "Lib.Proto"
local KineticSpecies = require "App.Species.KineticSpecies"
local Mpi = require "Comm.Mpi"
local VlasovEq = require "Eq.Vlasov"
local Updater = require "Updater"
local DataStruct = require "DataStruct"
local Time = require "Lib.Time"
local ffi = require "ffi"

local VlasovSpecies = Proto(KineticSpecies)

-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB  = 1
local SP_BC_REFLECT = 3
local SP_BC_EXTERN  = 4
local SP_BC_COPY    = 5
-- AHH: This was 2 but seems that is unstable. So using plain copy.
local SP_BC_OPEN      = SP_BC_COPY
local SP_BC_ZEROFLUX  = 6
local SP_BC_RESERVOIR = 7

VlasovSpecies.bcAbsorb    = SP_BC_ABSORB     -- Absorb all particles.
VlasovSpecies.bcOpen      = SP_BC_OPEN       -- Zero gradient.
VlasovSpecies.bcCopy      = SP_BC_COPY       -- Copy stuff.
VlasovSpecies.bcReflect   = SP_BC_REFLECT    -- Specular reflection.
VlasovSpecies.bcExternal  = SP_BC_EXTERN     -- Load external BC file.
VlasovSpecies.bcZeroFlux  = SP_BC_ZEROFLUX
VlasovSpecies.bcReservoir = SP_BC_RESERVOIR

function VlasovSpecies:alloc(nRkDup)
   -- Allocate distribution function.
   VlasovSpecies.super.alloc(self, nRkDup)

   -- Allocate fields to store coupling moments (for use in coupling
   -- to field and collisions).
   self.numDensity = self:allocMoment()
   self.momDensity = self:allocVectorMoment(self.vdim)
   self.ptclEnergy = self:allocMoment()

   -- Allocate field to accumulate funcField if any.
   self.totalEmField = self:allocVectorMoment(8)     -- 8 components of EM field.

   -- Allocate field for external forces if any.
   self.vExtForce = self:allocVectorMoment(self.vdim)

   -- Allocate moment array for integrated moments (n, n*u_i, sum_i n*u_i^2, sum_i n*T_ii).
   self.flow                 = self:allocVectorMoment(self.vdim)
   self.kineticEnergyDensity = self:allocMoment()
   self.thermalEnergyDensity = self:allocMoment()
end

function VlasovSpecies:fullInit(appTbl)
   VlasovSpecies.super.fullInit(self, appTbl)

   local tbl = self.tbl
   -- If there is an external force, get the force function.
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
   -- Run the KineticSpecies 'createSolver()' to initialize the
   -- collisions solver.
   VlasovSpecies.super.createSolver(self)

   -- Create updater to advance solution by one time-step.
   local vlasovEqn = VlasovEq {
      onGrid           = self.grid,
      phaseBasis       = self.basis,
      confBasis        = self.confBasis,
      charge           = self.charge,
      mass             = self.mass,
      hasElectricField = hasE,
      hasMagneticField = hasB,
   }

   -- Must apply zero-flux BCs in velocity directions.
   local zfd = { }
   for d = 1, self.vdim do 
     table.insert(self.zeroFluxDirections, self.cdim+d)
   end

   self.solver = Updater.HyperDisCont {
      onGrid             = self.grid,
      basis              = self.basis,
      cfl                = self.cfl,
      equation           = vlasovEqn,
      zeroFluxDirections = self.zeroFluxDirections,
   }

   -- Create updaters to compute various moments.
   self.numDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "M0",
   }
   self.momDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "M1i",
   }
   self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "M2",
   }
   -- Create updater to compute M0, M1i, M2 moments sequentially.
   -- If collisions are LBO, the following also computes boundary corrections and, if polyOrder=1, star moments.
   self.fiveMomentsCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "FiveMoments",
   }
   if self.needSelfPrimMom then
      -- This is used in calcCouplingMoments to reduce overhead and multiplications.
      -- If collisions are LBO, the following also computes boundary corrections and, if polyOrder=1, star moments.
      self.fiveMomentsLBOCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,
         phaseBasis = self.basis,
         confBasis  = self.confBasis,
         moment     = "FiveMomentsLBO",
      }
      if self.needCorrectedSelfPrimMom then
         self.primMomSelf = Updater.SelfPrimMoments {
            onGrid     = self.confGrid,
            phaseBasis = self.basis,
            confBasis  = self.confBasis,
            operator   = "VmLBO",
         }
      end
   end

   -- Updaters for the primitive moments.
   -- These will be used to compute n*u^2 and n*T for computing integrated moments.
   self.confDiv = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
   }
   self.confDotProduct = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "DotProduct",
   }

   if self.vlasovExtForceFunc then
      self.evalVlasovExtForce = Updater.ProjectOnBasis {
         onGrid          = self.confGrid,
         basis           = self.confBasis,
         evaluate        = self.vlasovExtForceFunc,
         projectOnGhosts = false
      }
   end

   self.tmCouplingMom = 0.0    -- For timer.
end

function VlasovSpecies:initCrossSpeciesCoupling(species)
   -- This method establishes the interaction between different
   -- species that is not mediated by the field (solver), like
   -- collisions.

   -- Function to find the index of an element in table.
   local function findInd(tblIn, el)
      for i, v in ipairs(tblIn) do
         if v == el then
            return i
         end
      end
      return #tblIn+1    -- If not found return a number larger than the length of the table.
   end

   -- Create a double nested table indicating who collides with whom.
   self.collPairs               = {}
   for sN, _ in pairs(species) do
      self.collPairs[sN]        = {}
      if next(species[sN].collisions) then -- Species[].collisions initialized as an empty table regardless. Need next. 
         -- This species collides with someone.
         local selfCollCurr, crossCollCurr, collSpecCurr
         -- Obtain the boolean indicating if self/cross collisions affect the sN species.
         for nm, _ in pairs(species[sN].collisions) do
            selfCollCurr  = species[sN].collisions[nm].selfCollisions
            crossCollCurr = species[sN].collisions[nm].crossCollisions
            collSpecsCurr = species[sN].collisions[nm].collidingSpecies
         end
         for sO, _ in pairs(species) do
            if sN == sO then
               self.collPairs[sN][sO] = selfCollCurr
            else
               if crossCollCurr then
                  local specInd = findInd(collSpecsCurr, sO)
                  if specInd < (#collSpecsCurr+1) then
                     self.collPairs[sN][sO] = true
                  else
                     self.collPairs[sN][sO] = false
                  end
               else
                  self.collPairs[sN][sO] = false
               end
            end
         end
      else
         -- This species does not collide with anyone.
         for sO, _ in pairs(species) do
            self.collPairs[sN][sO] = false
         end
      end
   end

   -- Determine if self primitive moments and boundary corrections are needed.
   -- If a pair of species only has cross-species collisions (no self-collisions)
   -- then the self-primitive moments may be computed without boundary corrections.
   self.needSelfPrimMom          = false
   self.needCorrectedSelfPrimMom = false
   if self.collPairs[self.name][self.name] then
      self.needSelfPrimMom          = true
      self.needCorrectedSelfPrimMom = true
   end
   for sO, _ in pairs(species) do
      if self.collPairs[self.name][sO] or self.collPairs[sO][self.name] then
         self.needSelfPrimMom = true
         if self.collPairs[sO][sO] then 
            self.needCorrectedSelfPrimMom = true
         end
      end
   end

   if self.needSelfPrimMom then
      -- Allocate fields to store self-species primitive moments.
      self.uSelf    = self:allocVectorMoment(self.vdim)
      self.vtSqSelf = self:allocMoment()

      -- Allocate fields for boundary corrections.
      self.m1Correction = self:allocVectorMoment(self.vdim)
      self.m2Correction = self:allocMoment()

      -- Allocate fields for star moments (only used with polyOrder=1).
      if (self.basis:polyOrder()==1) then
         self.m0Star = self:allocMoment()
         self.m1Star = self:allocVectorMoment(self.vdim)
         self.m2Star = self:allocMoment()
      end
   end

   -- Allocate fieds to store cross-species primitive moments.
   self.uCross    = {}
   self.vtSqCross = {}
   for sN, _ in pairs(species) do
      if sN ~= self.name then
         -- Flags for couplingMoments, boundary corrections, star moments,
         -- self primitive moments, cross primitive moments.
         self.momentFlags[5][sN] = false
      end

      for sO, _ in pairs(species) do
         -- Allocate space for this species' cross-primitive moments
         -- only if some other species collides with it.
         if (sN ~= sO) and (self.collPairs[sN][sO] or self.collPairs[sO][sN])
             and (self.uCross[string.gsub(sO .. sN, self.name, "")] == nil) then
            self.uCross[sO]    = self:allocVectorMoment(self.vdim)
            self.vtSqCross[sO] = self:allocMoment()
         end
      end
   end

end

function VlasovSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   -- Accumulate functional Maxwell fields (if needed).
   local emField      = emIn[1]:rkStepperFields()[inIdx]
   local emFuncField  = emIn[2]:rkStepperFields()[1]
   local totalEmField = self.totalEmField
   totalEmField:clear(0.0)

   local qbym = self.charge/self.mass

   if emField then totalEmField:accumulate(qbym, emField) end
   if emFuncField then totalEmField:accumulate(qbym, emFuncField) end
   
   -- If external force present (gravity, body force, etc.) accumulate it to electric field.
   if self.vlasovExtForceFunc then
      local vExtForce = self.vExtForce
      self.evalVlasovExtForce:advance(tCurr, {}, {vExtForce})

      -- Need to barrier over the shared communicator before accumulating force onto electric field.
      Mpi.Barrier(self.grid:commSet().sharedComm)

      -- Analogous to the current, the external force only gets accumulated onto the electric field.
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
      fRhsOut:clear(0.0)    -- No RHS.
   end
   -- Perform the collision update.
   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
         c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
         c:advance(tCurr, fIn, species, fRhsOut)
         -- The full 'species' list is needed for the cross-species
         -- collisions.
      end
   end
   if self.sourceSteadyState and self.evolveSources then
      local localEdgeFlux = ffi.new("double[3]")
      localEdgeFlux[0] = 0.0
      localEdgeFlux[1] = 0.0
      localEdgeFlux[2] = 0.0
      local flux = self:fluidMoments()[2]
      local fluxIndexer, fluxItr = flux:genIndexer(), flux:get(1)
      for idx in flux:localRangeIter() do
	 if idx[1] == self.grid:upper(1) then
	    flux:fill(fluxIndexer(idx), fluxItr)
	    localEdgeFlux[0] = fluxItr[1] / math.sqrt(2)
	 end
      end
      local globalEdgeFlux = ffi.new("double[3]")
      Mpi.Allreduce(localEdgeFlux, globalEdgeFlux, 1,
		    Mpi.DOUBLE, Mpi.MAX, self.grid:commSet().comm)

      local densFactor = globalEdgeFlux[0]/100.0--*self.dtGlobal[0] -- Hard-coded L !!
      fRhsOut:accumulate(densFactor, self.fSource)
   elseif self.fSource and self.evolveSources then
      -- add source it to the RHS
      fRhsOut:accumulate(self.sourceTimeDependence(tCurr), self.fSource)
   end
end

function VlasovSpecies:createDiagnostics()
   -- Create updater to compute volume-integrated moments
   -- function to check if integrated moment name is correct.
   local function isIntegratedMomentNameGood(nm)
      if nm == "intM0" or nm == "intM1i" or nm == "intM2Flow" or nm == "intM2Thermal" or nm == "intL2" then
         return true
      end
      return false
   end

   local numCompInt = {}
   numCompInt["intM0"]        = 1
   numCompInt["intM1i"]       = self.vdim
   numCompInt["intM2Flow"]    = 1
   numCompInt["intM2Thermal"] = 1
   numCompInt["intL2"]        = 1

   self.diagnosticIntegratedMomentFields   = { }
   self.diagnosticIntegratedMomentUpdaters = { } 
   -- Allocate space to store moments and create moment updater.
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
      if isIntegratedMomentNameGood(mom) then
         self.diagnosticIntegratedMomentFields[mom] = DataStruct.DynVector {
            numComponents = numCompInt[mom],
         }
         if mom == "intL2" then
            self.diagnosticIntegratedMomentUpdaters[mom] = Updater.CartFieldIntegratedQuantCalc {
               onGrid        = self.grid,
               basis         = self.basis,
               numComponents = numCompInt[mom],
               quantity      = "V2"
            }
         else
            self.diagnosticIntegratedMomentUpdaters[mom] = Updater.CartFieldIntegratedQuantCalc {
               onGrid        = self.confGrid,
               basis         = self.confBasis,
               numComponents = numCompInt[mom],
               quantity      = "V"
            }
         end
      else
         assert(false, string.format("Integrated Moment %s not valid", mom))
      end
   end

   -- Function to check if moment name is correct.
   local function isMomentNameGood(nm)
      return Updater.DistFuncMomentCalc:isMomentNameGood(nm)
   end
   -- Diagnostics computed with weak binary operations as diagnostic.
   -- Check if diagnostic name is correct.
   local function isWeakMomentNameGood(nm)
      return nm == "u" or nm == "vtSq"
   end
   local function isAuxMomentNameGood(nm)
      return nm == "uCross"
   end
   local function contains(table, element)
     for _, value in pairs(table) do
       if value == element then
         return true
       end
     end
     return false
   end

   local numComp   = {}
   numComp["M0"]   = 1
   numComp["M1i"]  = self.vdim
   numComp["M2ij"] = self.vdim*(self.vdim+1)/2
   numComp["M2"]   = 1
   numComp["M3i"]  = self.vdim
   numComp["u"]    = self.vdim
   numComp["vtSq"] = 1

   self.diagnosticMomentFields   = { }
   self.diagnosticMomentUpdaters = { } 
   self.diagnosticWeakMoments    = { }
   self.diagnosticAuxMoments     = { }
   self.weakMomentOpFields       = { }
   self.weakMomentScaleFac       = { }
   -- Create weak multiplication and division operations.
   self.weakDotProduct = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "DotProduct",
      onGhosts  = true,
   }
   self.weakDivision = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
      onGhosts  = true,
   }

   -- Sort moments into diagnosticWeakMoments and diagnosticAuxMoments.
   for i, mom in pairs(self.diagnosticMoments) do
      if isWeakMomentNameGood(mom) then
         -- Remove moment name from self.diagnosticMoments list, and add it to self.diagnosticWeakMoments list.
         self.diagnosticWeakMoments[mom] = true
         self.diagnosticMoments[i]       = nil
      elseif isAuxMomentNameGood(mom) then
         -- Remove moment name from self.diagnosticMoments list, and add it to self.diagnosticAuxMoments list.
         self.diagnosticAuxMoments[mom] = true
         self.diagnosticMoments[i]      = nil
      end
   end

   -- Make sure we have the updaters needed to calculate all the weak moments.
   for mom, _ in pairs(self.diagnosticWeakMoments) do
      -- All weak moments require M0 = density.
      if not contains(self.diagnosticMoments, "M0") then
         table.insert(self.diagnosticMoments, "M0")
      end

      if mom == "u" then
         if not contains(self.diagnosticMoments, "M1i") then
            table.insert(self.diagnosticMoments, "M1i")
         end
      end
      if mom == "vtSq" then
         if not contains(self.diagnosticMoments, "M2") then
            table.insert(self.diagnosticMoments, "M2")
         end
         if not self.diagnosticWeakMoments["u"] then
            self.diagnosticWeakMoments["u"] = true
         end
      end
   end

   -- Allocate space to store moments and create moment updater.
   for i, mom in ipairs(self.diagnosticMoments) do
      if isMomentNameGood(mom) then
         self.diagnosticMomentFields[mom] = DataStruct.Field {
            onGrid        = self.confGrid,
            numComponents = self.confBasis:numBasis()*numComp[mom],
            ghost         = {1, 1}
         }
         self.diagnosticMomentUpdaters[mom] = Updater.DistFuncMomentCalc {
            onGrid     = self.grid,
            phaseBasis = self.basis,
            confBasis  = self.confBasis,
            moment     = mom,
         }
      else
         assert(false, string.format("Moment %s not valid", mom))
      end
   end

   for mom, _ in pairs(self.diagnosticWeakMoments) do
      if isWeakMomentNameGood(mom) then
         self.diagnosticMomentFields[mom] = DataStruct.Field {
            onGrid        = self.confGrid,
            numComponents = self.confBasis:numBasis()*numComp[mom],
            ghost         = {1, 1}
         }
      else
         assert(false, string.format("Moment %s not valid", mom))
      end

      if mom == "u" then
         self.weakMomentOpFields["u"] = {self.diagnosticMomentFields["M0"], self.diagnosticMomentFields["M1i"]}
      elseif mom == "vtSq" then
         self.weakMomentOpFields["vtSq"] = {self.diagnosticMomentFields["M0"], self.diagnosticMomentFields["M2"]}
         self.weakMomentScaleFac["vtSq"] = 1.0/self.vdim
      end
   end

end

-- BC functions.
function VlasovSpecies:bcReflectFunc(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "flip".
   self.basis:flipSign(dir, fIn, fOut)
   self.basis:flipSign(dir+self.cdim, fOut, fOut)
end

function VlasovSpecies:bcExternFunc(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "flip".
   local velIdx = {}
   for d = 1, self.vdim do
      velIdx[d] = idxIn[self.cdim + d]
   end
   self.wallFunction[1](velIdx, fIn, fOut)
end

function VlasovSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- Need to wrap member functions so that self is passed.
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
   -- Compute moments needed in coupling to fields and collisions.
   local fIn = self:rkStepperFields()[rkIdx]
   if self.needSelfPrimMom then
      self.fiveMomentsLBOCalc:advance(tCurr, {fIn}, { self.numDensity, self.momDensity, self.ptclEnergy, 
                                                      self.m1Correction, self.m2Correction,
                                                      self.m0Star, self.m1Star, self.m2Star })
      if self.needCorrectedSelfPrimMom then
         -- Also compute self-primitive moments u and vtSq.
         self.primMomSelf:advance(tCurr, {self.numDensity, self.momDensity, self.ptclEnergy,
                                          self.m1Correction, self.m2Correction, 
                                          self.m0Star, self.m1Star, self.m2Star}, {self.uSelf, self.vtSqSelf})
      else
         -- Compute self-primitive moments with binOp updater.
         self.confDiv:advance(tCurr, {self.numDensity, self.momDensity}, {self.uSelf})
         self.confDotProduct:advance(tCurr, {self.uSelf, self.momDensity}, {self.kineticEnergyDensity})
         self.thermalEnergyDensity:combine( 1.0/self.vdim, self.ptclEnergy,
                                           -1.0/self.vdim, self.kineticEnergyDensity )
         self.confDiv:advance(tCurr, {self.numDensity, self.thermalEnergyDensity}, {self.vtSqSelf})
      end
      -- Indicate that moments, boundary corrections, star moments
      -- and self-primitive moments have been computed.
      for iF=1,4 do
         self.momentFlags[iF] = true
      end
   else
      self.momDensityCalc:advance(tCurr, {fIn}, { self.momDensity })
      -- Indicate that first moment has been computed.
      self.momentFlags[1] = true
   end
   self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart

end

-- Function to compute n, u, nu^2, and nT for use in integrated moment routine.
function VlasovSpecies:calcDiagnosticIntegratedMoments(tCurr)
   -- First compute M0, M1i, M2.
   local fIn = self:rkStepperFields()[1]
   self.fiveMomentsCalc:advance(tCurr, {fIn}, { self.numDensity, self.momDensity, self.ptclEnergy })

   -- Compute n*u^2 from n*u and n.
   self.confDiv:advance(0., {self.numDensity, self.momDensity}, {self.flow})
   self.confDotProduct:advance(0., {self.flow, self.momDensity}, {self.kineticEnergyDensity})

   -- Compute VDIM*n*T from M2 and kinetic energy density.
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

function VlasovSpecies:calcDiagnosticWeakMoments()
   VlasovSpecies.super.calcDiagnosticWeakMoments(self)
   if self.diagnosticWeakMoments["vtSq"] then
      -- Need to subtract (u^2)/vdim from vtSq (which at this point holds M2/(vdim*M0)).
      -- u is calculated in KineticEnergySpecies:calcDiagnositcWeakMoments().
      self.weakDotProduct:advance(0.0,
         {self.diagnosticMomentFields["u"], self.diagnosticMomentFields["u"]}, {self.kineticEnergyDensity})
      self.diagnosticMomentFields["vtSq"]:accumulate(-self.weakMomentScaleFac["vtSq"], self.kineticEnergyDensity)
   end
end

function VlasovSpecies:fluidMoments()
   return { self.numDensity, self.momDensity, self.ptclEnergy }
end

function VlasovSpecies:boundaryCorrections()
   return { self.m1Correction, self.m2Correction }
end

function VlasovSpecies:starMoments()
   return { self.m0Star, self.m1Star, self.m2Star }
end

function VlasovSpecies:selfPrimitiveMoments()
   return { self.uSelf, self.vtSqSelf }
end

function VlasovSpecies:crossPrimitiveMoments(otherSpeciesName)
   return { self.uCross[otherSpeciesName], self.vtSqCross[otherSpeciesName] }
end

function VlasovSpecies:getNumDensity(rkIdx)
   -- If no rkIdx specified, assume numDensity has already been calculated.
   if rkIdx == nil then return self.numDensity end 

   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]
   self.numDensityCalc:advance(nil, {fIn}, { self.numDensity })

   self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart

   return self.numDensity
end

function VlasovSpecies:getMomDensity(rkIdx)
   -- If no rkIdx specified, assume momDensity has already been calculated.
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
