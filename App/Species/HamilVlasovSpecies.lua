local Proto = require "Lib.Proto"
local KineticSpecies = require "App.Species.KineticSpecies"
local HamilVlasovEq = require "Eq.HamilVlasov"
local Updater = require "Updater"
local DataStruct = require "DataStruct"

local HamilVlasovSpecies = Proto(KineticSpecies)

function HamilVlasovSpecies:fullInit(appTbl)
   HamilVlasovSpecies.super.fullInit(self, appTbl)

   self.hamilDisContDirs = self.tbl.hamilDisContDirs
end

function HamilVlasovSpecies:alloc(nRkDup)
   -- allocate distribution function
   HamilVlasovSpecies.super.alloc(self, nRkDup)

   -- allocate fields to store coupling moments (for use in coupling
   -- to field and collisions)
   self.numDensity = self:allocMoment()
   self.momDensity = self:allocVectorMoment(self.vdim)
   self.ptclEnergy = self:allocMoment()
end

function HamilVlasovSpecies:allocMomCouplingFields()
   assert(false, "HamilVlasovSpecies:allocMomCouplingFields should not be called. Field object should allocate its own coupling fields")
end

function HamilVlasovSpecies:createSolver(hasPhi, hasA)
   -- create updater to advance solution by one time-step
   self.vlasovEqn = HamilVlasovEq {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      charge = self.charge,
      mass = self.mass,
      hasPhi = hasPhi,
      hasA = hasA,
      hamilDisContDirs = self.hamilDisContDirs, 
   }

   -- must apply zero-flux BCs in velocity directions
   local zfd = { }
   for d = 1, self.vdim do zfd[d] = self.cdim+d end

   self.solver = Updater.HyperDisCont {
      onGrid = self.grid,
      basis = self.basis,
      cfl = self.cfl,
      equation = self.vlasovEqn,
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

function HamilVlasovSpecies:forwardEuler(tCurr, dt, fIn, emIn, fOut)
   if self.evolve then
      local em = emIn[1]
      local emFunc = emIn[2]
      return self.solver:advance(tCurr, dt, {fIn, em, emFunc}, {fOut})
   else
      fOut:copy(fIn) -- just copy stuff over
      return true, GKYL_MAX_DOUBLE
   end
end

function HamilVlasovSpecies:createDiagnostics()
   HamilVlasovSpecies.super.createDiagnostics(self)

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

function HamilVlasovSpecies:write(tm)
   HamilVlasovSpecies.super.write(self, tm)
   if self.distIoTrigger(tm) then
      self.vlasovEqn:writeHamiltonian(self.distIo, tm) 
   end
end

function HamilVlasovSpecies:calcCouplingMoments(tCurr, dt, fIn)
   -- compute moments needed in coupling to fields and collisions
   self.numDensityCalc:advance(tCurr, dt, {fIn}, { self.numDensity })
   self.momDensityCalc:advance(tCurr, dt, {fIn}, { self.momDensity })
   self.ptclEnergyCalc:advance(tCurr, dt, {fIn}, { self.ptclEnergy })
end

function HamilVlasovSpecies:fluidMoments()
   return { self.numDensity, self.momDensity, self.ptclEnergy } 
end

function HamilVlasovSpecies:getNumDensity()
   return self.numDensity
end

function HamilVlasovSpecies:getMomDensity()
   return self.momDensity
end

function HamilVlasovSpecies:getDens()
   return self.numDensity
end

function HamilVlasovSpecies:momCalcTime()
   local tm = self.numDensityCalc.totalTime
   tm = tm + self.momDensityCalc.totalTime
   for i, mom in ipairs(self.diagnosticMoments) do
      tm = tm + self.diagnosticMomentUpdaters[i].totalTime
   end
   return tm
end

return HamilVlasovSpecies
