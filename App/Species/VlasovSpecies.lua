local Proto = require "Lib.Proto"
local KineticSpecies = require "App.Species.KineticSpecies"
local VlasovEq = require "Eq.Vlasov"
local Updater = require "Updater"
local DataStruct = require "DataStruct"

local VlasovSpecies = Proto(KineticSpecies)

-- add constants to object indicate various supported boundary conditions
local SP_BC_ABSORB = 1
local SP_BC_OPEN = 2
local SP_BC_REFLECT = 3
local SP_BC_SEE = 4
VlasovSpecies.bcAbsorb = SP_BC_ABSORB -- absorb all particles
VlasovSpecies.bcOpen = SP_BC_OPEN -- zero gradient
VlasovSpecies.bcReflect = SP_BC_REFLECT -- specular reflection
VlasovSpecies.bcSEE = SP_BC_SEE -- specular reflection

function VlasovSpecies:alloc(nRkDup)
   -- allocate distribution function
   VlasovSpecies.super.alloc(self, nRkDup)

   -- allocate fields to store coupling moments (for use in coupling
   -- to field and collisions)
   self.numDensity = self:allocMoment()
   self.momDensity = self:allocVectorMoment(self.vdim)
   self.ptclEnergy = self:allocMoment()
end

function VlasovSpecies:allocMomCouplingFields()
   -- only need currents for coupling to fields (returning a table
   -- with single entry, i.e. space to store currents)
   return {currentDensity = self:allocVectorMoment(self.vdim)}
end


function VlasovSpecies:createSolver(hasE, hasB)
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
end

function VlasovSpecies:forwardEuler(tCurr, dt, fIn, emIn, fOut)
   -- accumulate functional Maxwell fields (if needed)
   local emFields = emIn[1]
   local emFuncFields = emIn[2]
   local totalEmField = nil
   if emFuncFields then
      if emFields then
         emFuncFields:accumulate(1.0, emFields)
      end
      totalEmField = emFuncFields
   else
      totalEmField = emFields
   end

   if self.evolve then
      return self.solver:advance(tCurr, dt, {fIn, totalEmField}, {fOut})
   else
      fOut:copy(fIn) -- just copy stuff over
      return true, GKYL_MAX_DOUBLE
   end
end

function VlasovSpecies:createDiagnostics()
   VlasovSpecies.super.createDiagnostics(self)

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
   -- get handle to function to compute basis functions at specified coordinates
   local bName = ""
   if self.basis:id() == "serendipity" then
      bName = "Serendip"
   else
      bName = "MaxOrder"
   end
   local fName = "Basis._data.Modal" .. bName .. "BasisReflect" .. self.cdim .. "x" .. self.vdim .. "v"
   local _m = require(fName)
   local polyOrder = self.basis:polyOrder()
   _m[polyOrder](dir, fIn, fOut) -- function to flip sign of both configuration and velocity component

   -- NRM: note that same behavior can be obtained from the following
   --self.basis:flipSign(dir, fIn, fOut)
   --self.basis:flipSign(dir+self.cdim, fOut, fOut)
end

function VlasovSpecies:bcSeeFunc(dir, tm, idxIn, fIn, fOut)
   -- requires skinLoop = "integrate"

end

function VlasovSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- need to wrap member functions so that self is passed
   local function bcAbsorbFunc(...) return self:bcAbsorbFunc(...) end
   local function bcOpenFunc(...) return self:bcOpenFunc(...) end
   local function bcReflectFunc(...) return self:bcReflectFunc(...) end
   local function bcSeeFunc(...) return self:bcSeeFunc(...) end

   local vdir = dir + self.cdim

   if bcType == SP_BC_ABSORB then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcAbsorbFunc }, "pointwise"))
   elseif bcType == SP_BC_OPEN then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcOpenFunc }, "pointwise"))
   elseif bcType == SP_BC_REFLECT then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcReflectFunc }, "flip"))
   elseif bcType == SP_BC_SEE then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcSeeFunc }, "integrate"))
   else
      assert(false, "VlasovSpecies: Unsupported BC type!")
   end
end

function VlasovSpecies:calcCouplingMoments(tCurr, dt, fIn)
   -- compute moments needed in coupling to fields and collisions
   self.numDensityCalc:advance(tCurr, dt, {fIn}, { self.numDensity })
   self.momDensityCalc:advance(tCurr, dt, {fIn}, { self.momDensity })
   self.ptclEnergyCalc:advance(tCurr, dt, {fIn}, { self.ptclEnergy })
end

function VlasovSpecies:fluidMoments()
   return { self.numDensity, self.momDensity, self.ptclEnergy } 
end

function VlasovSpecies:getNumDensity()
   return self.numDensity
end

function VlasovSpecies:getMomDensity()
   return self.momDensity
end

function VlasovSpecies:momCalcTime()
   local tm = self.momDensityCalc.totalTime
   for i, mom in ipairs(self.diagnosticMoments) do
      tm = tm + self.diagnosticMomentUpdaters[i].totalTime
   end
   return tm
end

return VlasovSpecies
