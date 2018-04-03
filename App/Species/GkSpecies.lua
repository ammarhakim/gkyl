local Proto = require "Lib.Proto"
local KineticSpecies = require "App.Species.KineticSpecies"
local GkEq = require "Eq.Gyrokinetic"
local Updater = require "Updater"
local DataStruct = require "DataStruct"

local GkSpecies = Proto(KineticSpecies)

function GkSpecies:alloc(nRkDup)
   -- allocate distribution function
   GkSpecies.super.alloc(self, nRkDup)

   -- allocate fields to store coupling moments (for use in coupling
   -- to field and collisions)
   self.dens = self:allocMoment()
   self.upar = self:allocMoment()
   self.ppar = self:allocMoment()
   self.pperp = self:allocMoment()
end

function GkSpecies:allocMomCouplingFields()
   assert(false, "GkSpecies:allocMomCouplingFields should not be called. Field object should allocate its own coupling fields")
end

function GkSpecies:createSolver(hasPhi, hasApar)
   -- create updater to advance solution by one time-step
   self.gkEqn = GkEq {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      charge = self.charge,
      mass = self.mass,
      hasPhi = hasPhi,
      hasApar = hasApar,
   }

   -- no update in mu direction (last velocity direction if present)
   local upd = {}
   for d = 1, self.cdim + 1 do upd[d] = d end

   self.solver = Updater.HyperDisCont {
      onGrid = self.grid,
      basis = self.basis,
      cfl = self.cfl,
      equation = self.gkEqn,
      zeroFluxDirections = {self.cdim+1},
      updateDirections = upd,
   }
   
   -- create updaters to compute various moments
   self.calcDens = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkDens",
   }
   self.calcUpar = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkUpar",
   }
   self.calcPpar = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkPpar",
   }
   if self.vdim > 1 then
      self.calcPperp = Updater.DistFuncMomentCalc {
         onGrid = self.grid,
         phaseBasis = self.basis,
         confBasis = self.confBasis,
         moment = "GkPperp",
      }
   end
end

function GkSpecies:forwardEuler(tCurr, dt, fIn, emIn, fOut)
   if self.evolve then
      local em = emIn[1]
      local emFunc = emIn[2]
      return self.solver:advance(tCurr, dt, {fIn, em, emFunc}, {fOut})
   else
      fOut:copy(fIn) -- just copy stuff over
      return true, GKYL_MAX_DOUBLE
   end
end

function GkSpecies:createDiagnostics()
   GkSpecies.super.createDiagnostics(self)
   
   -- function to check if moment name is correct
   local function isMomentNameGood(nm)
      if nm == "GkDens" or nm == "GkUpar" or nm == "GkPpar" or nm == "GkPperp" then
         return true
      end
      return false
   end

   self.diagnosticMomentFields = { }
   self.diagnosticMomentUpdaters = { } 
   -- allocate space to store moments and create moment updater
   for i, mom in ipairs(self.diagnosticMoments) do
      if isMomentNameGood(mom) then
         self.diagnosticMomentFields[i] = DataStruct.Field {
            onGrid = self.confGrid,
            numComponents = self.confBasis:numBasis(),
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

function GkSpecies:write(tm)
   GkSpecies.super.write(self, tm)
   if self.distIoTrigger(tm) then
      self.gkEqn:writeHamiltonian(self.distIo, tm) 
   end
end

function GkSpecies:calcCouplingMoments(tCurr, dt, fIn)
   -- compute moments needed in coupling to fields and collisions
   self.calcDens:advance(tCurr, dt, {fIn}, { self.dens })
   self.calcUpar:advance(tCurr, dt, {fIn}, { self.upar })
end

function GkSpecies:fluidMoments()
   return { self.dens, self.upar, self.ppar, self.pperp } 
end

function GkSpecies:getDens()
   return self.dens
end

function GkSpecies:getUpar()
   return self.upar
end

function GkSpecies:momCalcTime()
   local tm = self.calcDens.totalTime
   tm = tm + self.calcUpar.totalTime
   for i, mom in ipairs(self.diagnosticMoments) do
      tm = tm + self.diagnosticMomentUpdaters[i].totalTime
   end
   return tm
end

return GkSpecies
