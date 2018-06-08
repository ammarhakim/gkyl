local Proto = require "Lib.Proto"
local KineticSpecies = require "App.Species.KineticSpecies"
local VlasovEq = require "Eq.Vlasov"
local Updater = require "Updater"
local DataStruct = require "DataStruct"

local VlasovSpecies = Proto(KineticSpecies)

-- add constants to object indicate various supported boundary conditions
local SP_BC_ABSORB = 1
local SP_BC_REFLECT = 3
local SP_BC_REFLECT_QM = 4
local SP_BC_COPY = 5
-- AHH: This was 2 but seems that is unstable. So using plain copy
local SP_BC_OPEN = SP_BC_COPY

VlasovSpecies.bcAbsorb = SP_BC_ABSORB -- absorb all particles
VlasovSpecies.bcOpen = SP_BC_OPEN -- zero gradient
VlasovSpecies.bcCopy = SP_BC_COPY -- copy stuff
VlasovSpecies.bcReflect = SP_BC_REFLECT -- specular reflection
VlasovSpecies.bcReflectQM = SP_BC_REFLECT_QM -- specular reflection with < 1 probability

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
end

function VlasovSpecies:allocMomCouplingFields()
   -- only need currents for coupling to fields (returning a table
   -- with single entry, i.e. space to store currents)
   return {currentDensity = self:allocVectorMoment(self.vdim)}
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
   for d = 1, self.vdim do zfd[d] = self.cdim+d end

   self.solver = Updater.HyperDisCont {
      onGrid = self.grid,
      basis = self.basis,
      cfl = self.cfl,
      equation = vlasovEqn,
      zeroFluxDirections = zfd,
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

   -- create updater to evaluate source 
   if self.sourceFunc then 
      self.evalSource = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = self.sourceFunc,
         projectOnGhosts = true
      }
   end
end

function VlasovSpecies:forwardEuler(tCurr, dt, species, emIn, inIdx, outIdx)
   local fIn = self:rkStepperFields()[inIdx]
   local fOut = self:rkStepperFields()[outIdx]

   -- accumulate functional Maxwell fields (if needed)
   local emField = emIn[1]:rkStepperFields()[inIdx]
   local emFuncField = emIn[2]:rkStepperFields()[1]
   local totalEmField = nil
   if emFuncField then
      if emField then
	 self.totalEmField:combine(1.0, emFuncField, 1.0, emField)
	 totalEmField = self.totalEmField
      else
	 totalEmField = emFuncField
      end
   else
      totalEmField = emField
   end

   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   if self.evolveCollisionless then      
      status, dtSuggested = self.solver:advance(tCurr, dt,
						{fIn, totalEmField},
						{fOut})
      if self.sourceFunc then
        -- if there is a source, add it to the RHS
        local fSource = self.fSource
        self.evalSource:advance(tCurr, dt, {}, {fSource})
        fOut:accumulate(dt, fSource)
      end
   else
      fOut:copy(fIn) -- just copy stuff over
   end
   -- perform the collision update
   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
	 local collStatus, collDt = c:forwardEuler(tCurr, dt,
						   fIn, species,
						   fOut)
	 -- the full 'species' list is needed for the cross-species
	 -- collisions
	 status = status and collStatus
	 dtSuggested = math.min(dtSuggested, collDt)
      end
   end

   -- apply BCs
   self:applyBc(tCurr, dt, fOut)

   return status, dtSuggested
end

function VlasovSpecies:createDiagnostics()
   -- create updater to compute volume-integrated moments
   self.intMomentCalc = Updater.DistFuncIntegratedMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
   }

   -- function to check if moment name is correct
   local function isMomentNameGood(nm)
      if nm == "M0" or nm == "M1i" or nm == "M2ij" or nm == "M2" or nm == "M3i" then
         return true
      end
      return false
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
         self.diagnosticMomentFields[i] = DataStruct.Field {
            onGrid = self.confGrid,
            numComponents = self.confBasis:numBasis()*numComp[mom],
            ghost = {1, 1}
         }
         self.diagnosticMomentUpdaters[i] = Updater.DistFuncMomentCalc {
            onGrid = self.grid,
            phaseBasis = self.basis,
            confBasis = self.confBasis,
            moment = mom,
         }
      else
         assert(false, string.format("Moment %s not valid", mom))
      end
   end
end

-- BC functions
function VlasovSpecies:bcReflectFunc(dir, tm, idxIn, fIn, fOut)
   -- requires skinLoop = "flip"
   self.basis:flipSign(dir, fIn, fOut)
   self.basis:flipSign(dir+self.cdim, fOut, fOut)
end

function VlasovSpecies:bcReflectQMFunc(dir, tm, idxIn, fIn, fOut)
   -- requires skinLoop = "flip"
   self.basis:flipSign(dir, fIn, fOut)
   self.basis:flipSign(dir+self.cdim, fOut, fOut)

   local wallFunction = require "wall"
   local velIdx = {}
   velIdx[1] = idxIn[2]
   wallFunction[1](velIdx, fOut, fOut)
end

function VlasovSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- need to wrap member functions so that self is passed
   local function bcAbsorbFunc(...) return self:bcAbsorbFunc(...) end
   local function bcCopyFunc(...) return self:bcCopyFunc(...) end
   local function bcOpenFunc(...) return self:bcOpenFunc(...) end
   local function bcReflectFunc(...) return self:bcReflectFunc(...) end
   local function bcReflectQMFunc(...) return self:bcReflectQMFunc(...) end

   local vdir = dir + self.cdim

   if bcType == SP_BC_ABSORB then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcAbsorbFunc }, "pointwise"))
   elseif bcType == SP_BC_OPEN then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcCopyFunc }, "pointwise"))
   elseif bcType == SP_BC_COPY then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcCopyFunc }, "pointwise"))
   elseif bcType == SP_BC_REFLECT then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcReflectFunc }, "flip"))
   elseif bcType == SP_BC_REFLECT_QM then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcReflectQMFunc }, "flip"))
   else
      assert(false, "VlasovSpecies: Unsupported BC type!")
   end
end

function VlasovSpecies:calcCouplingMoments(tCurr, dt, rkIdx)
   -- compute moments needed in coupling to fields and collisions
   local fIn = self:rkStepperFields()[rkIdx]
   self.numDensityCalc:advance(tCurr, dt, {fIn}, { self.numDensity })
   self.momDensityCalc:advance(tCurr, dt, {fIn}, { self.momDensity })
   self.ptclEnergyCalc:advance(tCurr, dt, {fIn}, { self.ptclEnergy })
end

function VlasovSpecies:fluidMoments()
   return { self.numDensity, self.momDensity, self.ptclEnergy }
end

function VlasovSpecies:getNumDensity(rkIdx)
   -- if no rkIdx specified, assume numDensity has already been calculated
   if rkIdx == nil then return self.numDensity end 

   local fIn = self:rkStepperFields()[rkIdx]
   self.numDensityCalc:advance(tCurr, dt, {fIn}, { self.numDensity })
   return self.numDensity
end

function VlasovSpecies:getMomDensity(rkIdx)
   -- if no rkIdx specified, assume momDensity has already been calculated
   if rkIdx == nil then return self.momDensity end 

   local fIn = self:rkStepperFields()[rkIdx]
   self.momDensityCalc:advance(tCurr, dt, {fIn}, { self.momDensity })
   return self.momDensity
end

function VlasovSpecies:momCalcTime()
   local tm = self.momDensityCalc.totalTime + self.numDensityCalc.totalTime + self.ptclEnergyCalc.totalTime
   for i, mom in ipairs(self.diagnosticMoments) do
      tm = tm + self.diagnosticMomentUpdaters[i].totalTime
   end
   return tm
end

-- please test this for higher than 1x1v... 
function VlasovSpecies:Maxwellian(xn, n0, T0, vdn)
   local vdn = vdn or {0, 0, 0}
   local vt2 = T0/self.mass
   local v2 = 0.0
   for d = self.cdim+1, self.cdim+self.vdim do
     v2 = v2 + (xn[d] - vdn[d-self.cdim])^2
   end
   return n0/math.sqrt(2*math.pi*vt2)^self.cdim*math.exp(-v2/(2*vt2))
end

return VlasovSpecies
