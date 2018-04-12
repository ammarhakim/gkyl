local Proto = require "Lib.Proto"
local KineticSpecies = require "App.Species.KineticSpecies"
local GkEq = require "Eq.Gyrokinetic"
local Updater = require "Updater"
local DataStruct = require "DataStruct"

local GkSpecies = Proto(KineticSpecies)

-- add constants to object indicate various supported boundary conditions
local SP_BC_ABSORB = 1
local SP_BC_OPEN = 2
local SP_BC_REFLECT = 3
local SP_BC_SHEATH = 4
GkSpecies.bcAbsorb = SP_BC_ABSORB -- absorb all particles
GkSpecies.bcOpen = SP_BC_OPEN -- zero gradient
GkSpecies.bcReflect = SP_BC_REFLECT -- specular reflection
GkSpecies.bcSheath = SP_BC_SHEATH -- specular reflection

function GkSpecies:alloc(nRkDup)
   -- allocate distribution function
   GkSpecies.super.alloc(self, nRkDup)

   -- allocate fields to store coupling moments (for use in coupling
   -- to field and collisions)
   self.dens = self:allocMoment()
   self.dens0 = self:allocMoment()
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

   -- account for factor of mass in definition of mu
   -- note that factor of 2*pi from gyrophase integration handled in GkMoment calculations
   self.momfac = 1/self.mass
   
   -- create updaters to compute various moments
   self.calcDens = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkDens",
      momfac = self.momfac,
   }
   self.calcUpar = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkUpar",
      momfac = self.momfac,
   }
   self.calcPpar = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkPpar",
      momfac = self.momfac,
   }
   if self.vdim > 1 then
      self.calcPperp = Updater.DistFuncMomentCalc {
         onGrid = self.grid,
         phaseBasis = self.basis,
         confBasis = self.confBasis,
         moment = "GkPperp",
         momfac = self.momfac,
      }
   end

   -- calculate background density averaged over simulation domain
   self.n0 = nil
   if self.f0 then
      self.calcDens:advance(0,0, {self.f0}, {self.dens0})
      local data
      local dynVec = DataStruct.DynVector { numComponents = 1 }
      -- integrate source
      local calcInt = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.confGrid,
         basis = self.confBasis,
         numComponents = 1,
      }
      calcInt:advance(0.0, 0.0, {self.dens0}, {dynVec})
      _, data = dynVec:lastData()
      self.n0 = data[1]/self.confGrid:gridVolume()
      print("Average density is " .. self.n0)
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
            momfac = self.momfac,
         }
      else
         assert(false, string.format("Moment %s not valid", mom))
      end
   end
end

-- BC functions
function GkSpecies:bcReflectFunc(dir, tm, xc, fIn, fOut)
   -- skinLoop should be "flip"
   self.basis:flipSign(dir, fIn, fOut)
   -- GK reflection only makes sense in z-vpar
   -- z is always last config space dir and vpar is always next
   if dir==self.cdim then 
      local vpardir=self.cdim+1 
      self.basis:flipSign(vpardir, fOut, fOut)
   end
end
function GkSpecies:bcPartialReflectFunc(dir, tm, xc, fIn, fOut, vcut)

end
function GkSpecies:bcSheathFunc(dir, tm, gridIn, fIn, fOut, phi)
   -- if phi not passed in, do nothing
   if phi == nil then return end

   -- calculate sheath cutoff velocity
   local vcut = nil

   local vpardir = self.cdim+1
   local vL = gridIn:cellLowerInDir(vpardir)
   local vR = gridIn:cellUpperInDir(vpardir)
   local vlower, vupper
   if abs(vR)>abs(vL) then
      vlower = abs(vL)
      vupper = abs(vR)
   else
      vlower = abs(vR)
      vupper = abs(vL)
   end
   if vcut < vlower then
      self:bcReflectFunc(dir, tm, xc, fIn, fOut)
   elseif vcut < vupper then
      self:bcPartialReflectFunc(dir, tm, xc, fIn, fOut, vcut)
   else
      self:bcAbsorbFunc(dir, tm, xc, fIn, fOut)
   end
end

function GkSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- need to wrap member functions so that self is passed
   local function bcAbsorbFunc(...) self:bcAbsorbFunc(...) end
   local function bcOpenFunc(...) self:bcOpenFunc(...) end
   local function bcReflectFunc(...) self:bcReflectFunc(...) end
   local function bcSheathFunc(...) self:bcSheathFunc(...) end
   
   local vdir = nil
   if dir==self.cdim then 
      vdir=self.cdim+1 
   end

   if bcType == SP_BC_ABSORB then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcAbsorb }, "pointwise"))
   elseif bcType == SP_BC_OPEN then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcOpen }, "pointwise"))
   -- note: reflection and sheath BCs only make sense in z direction,
   -- which is always last config space direction, i.e. dir = self.cdim
   elseif bcType == SP_BC_REFLECT and dir==self.cdim then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcReflect }, "flip"))
   elseif bcType == SP_BC_SHEATH and dir==self.cdim then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcSheath }, "flip"))
   else
      assert(false, "GkSpecies: Unsupported BC type!")
   end
end

-- sheath BCs need phi
function GkSpecies:getBcAux(potentials)
   return potentials.phi
end

function GkSpecies:calcCouplingMoments(tCurr, dt, fIn)
   -- compute moments needed in coupling to fields and collisions
   -- we only want perturbed moments so that we can calculate perturbed fields
   fIn:accumulate(-1, self.f0)
   self.calcDens:advance(tCurr, dt, {fIn}, { self.dens })
   self.calcUpar:advance(tCurr, dt, {fIn}, { self.upar })
   fIn:accumulate(1, self.f0)
end

function GkSpecies:calcDiagnosticMoments()
   self.distf[1]:accumulate(-1, self.f0)
   local numMoms = #self.diagnosticMoments
   for i = 1, numMoms do
      self.diagnosticMomentUpdaters[i]:advance(
	 0.0, 0.0, {self.distf[1]}, {self.diagnosticMomentFields[i]})
   end
   self.distf[1]:accumulate(1, self.f0)
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
