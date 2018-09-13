-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic species object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local KineticSpecies = require "App.Species.KineticSpecies"
local Gk = require "Eq.Gyrokinetic"
local Updater = require "Updater"
local DataStruct = require "DataStruct"
local Time = require "Lib.Time"

local GkSpecies = Proto(KineticSpecies)

-- add constants to object indicate various supported boundary conditions
local SP_BC_ABSORB = 1
local SP_BC_OPEN = 2
local SP_BC_REFLECT = 3
local SP_BC_SHEATH = 4
local SP_BC_ZEROFLUX = 5
GkSpecies.bcAbsorb = SP_BC_ABSORB -- absorb all particles
GkSpecies.bcOpen = SP_BC_OPEN -- zero gradient
GkSpecies.bcReflect = SP_BC_REFLECT -- specular reflection
GkSpecies.bcSheath = SP_BC_SHEATH -- specular reflection
GkSpecies.bcZeroFlux = SP_BC_ZEROFLUX -- zero flux

function GkSpecies:alloc(nRkDup)
   -- allocate distribution function
   GkSpecies.super.alloc(self, nRkDup)

   -- allocate fields to store coupling moments (for use in coupling
   -- to field and collisions)
   self.numDensity = self:allocMoment()
   self.momDensity = self:allocMoment()
   self.ptclEnergy = self:allocMoment()
   self.polarizationWeight = self:allocMoment()

   if self.gyavg then
      self.rho1 = self:allocDistf()
      self.rho2 = self:allocDistf()
      self.rho3 = self:allocDistf()
   end
end

function GkSpecies:allocMomCouplingFields()
   assert(false, "GkSpecies:allocMomCouplingFields should not be called. Field object should allocate its own coupling fields")
end

function GkSpecies:initDist()
   GkSpecies.super.initDist(self)

   -- calculate initial density averaged over simulation domain
   self.n0 = nil
   local dens0 = self:allocMoment()
   self.numDensityCalc:advance(0,0, {self.distf[1]}, {dens0})
   local data
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   -- integrate 
   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid = self.confGrid,
      basis = self.confBasis,
      numComponents = 1,
      quantity = "V"
   }
   calcInt:advance(0.0, 0.0, {dens0}, {dynVec})
   _, data = dynVec:lastData()
   self.n0 = data[1]/self.confGrid:gridVolume()
   --print("Average density is " .. self.n0)
end

function GkSpecies:createSolver(hasPhi, hasApar, funcField)
   -- run the KineticSpecies 'createSolver()' to initialize the
   -- collisions solver
   GkSpecies.super.createSolver(self,funcField)

   -- set up jacobian
   if funcField then
      -- save bmagFunc for later...
      self.bmagFunc = funcField.bmagFunc
      -- if vdim>1, get jacobian=bmag from geo, and multiply it by init functions for f0 and f
      self.jacobPhaseFunc = self.bmagFunc
      self.jacobGeoFunc = funcField.jacobGeoFunc
      if self.jacobPhaseFunc and self.vdim > 1 then
         local initFuncWithoutJacobian = self.initFunc
         self.initFunc = function (t, xn)
            local J = self.jacobPhaseFunc(t,xn)
            local f = initFuncWithoutJacobian(t,xn)
            return J*f
         end
         if self.initBackgroundFunc then
            local initBackgroundFuncWithoutJacobian = self.initBackgroundFunc
            self.initBackgroundFunc = function(t,xn)
               local J = self.jacobPhaseFunc(t,xn)
               local f0 = initBackgroundFuncWithoutJacobian(t,xn)
               return J*f0
            end
         end
      end
      if self.jacobGeoFunc then
         local initFuncWithoutJacobian = self.initFunc
         self.initFunc = function (t, xn)
            local J = self.jacobGeoFunc(t,xn)
            local f = initFuncWithoutJacobian(t,xn)
            return J*f
         end
         if self.initBackgroundFunc then
            local initBackgroundFuncWithoutJacobian = self.initBackgroundFunc
            self.initBackgroundFunc = function(t,xn)
               local J = self.jacobGeoFunc(t,xn)
               local f0 = initBackgroundFuncWithoutJacobian(t,xn)
               return J*f0
            end
         end
      end
      if self.cdim == 1 then 
         self.B0 = funcField.bmagFunc(0.0, {self.grid:mid(1)})
      elseif self.cdim == 2 then 
         self.B0 = funcField.bmagFunc(0.0, {self.grid:mid(1), self.grid:mid(1)})
      else
         self.B0 = funcField.bmagFunc(0.0, {self.grid:mid(1), self.grid:mid(1), self.grid:mid(2)})
      end
      self.bmag = assert(funcField.geo.bmag, "nil bmag")
   end

   if self.gyavg then
      -- set up geo fields needed for gyroaveraging
      local rho1Func = function (t, xn)
         local mu = xn[self.ndim]
         return math.sqrt(2*mu*self.mass*funcField.gxxFunc(t, xn)/(self.charge^2*funcField.bmagFunc(t, xn)))
      end
      local rho2Func = function (t, xn)
         local mu = xn[self.ndim]
         return funcField.gxyFunc(t,xn)*math.sqrt(2*mu*self.mass/(self.charge^2*funcField.gxxFunc(t, xn)*funcField.bmagFunc(t, xn)))
      end
      local rho3Func = function (t, xn)
         local mu = xn[self.ndim]
         return math.sqrt(2*mu*self.mass*(funcField.gxxFunc(t,xn)*funcField.gyyFunc(t,xn)-funcField.gxyFunc(t,xn)^2)/(self.charge^2*funcField.gxxFunc(t, xn)*funcField.bmagFunc(t, xn)))
      end
      local project1 = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = rho1Func,
         projectOnGhosts = true
      }
      project1:advance(0.0, 0.0, {}, {self.rho1})
      local project2 = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = rho2Func,
         projectOnGhosts = true
      }
      project2:advance(0.0, 0.0, {}, {self.rho2})
      local project3 = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = rho3Func,
         projectOnGhosts = true
      }
      project3:advance(0.0, 0.0, {}, {self.rho3})

      -- create solver for gyroaveraging potentials
      self.emGyavgSlvr = Updater.FemGyroaverage {
         onGrid = self.confGrid,
         confBasis = self.confBasis,
         phaseGrid = self.grid,
         phaseBasis = self.basis,
         rho1 = self.rho1,
         rho2 = self.rho2,
         rho3 = self.rho3,
         muOrder0 = true, -- cell-average in mu
      }

      -- create solver for gyroaveraging distribution function
      self.distfGyavgSlvr = Updater.FemGyroaverage {
         onGrid = self.confGrid,
         confBasis = self.confBasis,
         phaseGrid = self.grid,
         phaseBasis = self.basis,
         rho1 = self.rho1,
         rho2 = self.rho2,
         rho3 = self.rho3,
         integrate = true,
      }
   end

   -- create updater to advance solution by one time-step
   self.gkEqn = Gk.GkEq {
      onGrid = self.grid,
      confGrid = self.confGrid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      charge = self.charge,
      mass = self.mass,
      hasPhi = hasPhi,
      hasApar = hasApar,
      Bvars = funcField.bmagVars,
      hasSheathBcs = self.hasSheathBcs,
      positivity = self.positivity,
      gyavgSlvr = self.emGyavgSlvr,
   }

   -- no update in mu direction (last velocity direction if present)
   local upd = {}
   if hasApar then -- if electromagnetic only update conf dir surface terms on first step
      for d = 1, self.cdim do upd[d] = d end
   else
      for d = 1, self.cdim + 1 do upd[d] = d end
   end
   -- zero flux in vpar and mu
   table.insert(self.zeroFluxDirections, self.cdim+1)
   table.insert(self.zeroFluxDirections, self.cdim+2)

   self.solver = Updater.HyperDisCont {
      onGrid = self.grid,
      basis = self.basis,
      cfl = self.cfl,
      equation = self.gkEqn,
      zeroFluxDirections = self.zeroFluxDirections,
      updateDirections = upd,
      onlyIncrement = true, 
   }
   if hasApar then 
      -- set up solver that adds on volume term involving dApar/dt and the entire vpar surface term
      self.gkEqnStep2 = Gk.GkEqStep2 {
         onGrid = self.grid,
         phaseBasis = self.basis,
         confBasis = self.confBasis,
         charge = self.charge,
         mass = self.mass,
         Bvars = funcField.bmagVars,
         positivity = self.positivity,
      }
      -- note that the surface update for this term only involves the vpar direction
      self.solverStep2 = Updater.HyperDisCont {
         onGrid = self.grid,
         basis = self.basis,
         cfl = self.cfl,
         equation = self.gkEqnStep2,
         zeroFluxDirections = self.zeroFluxDirections,
         updateDirections = {self.cdim+1},
         clearOut = false,   -- continue accumulating into output field
         onlyIncrement = true, 
      }
   end
   
   -- create updaters to compute various moments
   self.numDensityCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkM0",
      gkfacs = {self.mass, self.bmag},
   }
   self.momDensityCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkM1",
      gkfacs = {self.mass, self.bmag},
   }
   self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkM2",
      gkfacs = {self.mass, self.bmag},
   }
   self.threeMomentsCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkThreeMoments",
      gkfacs = {self.mass, self.bmag},
   }
   
   self._firstMomentCalc = true  -- to avoid re-calculating moments when not evolving

   self.tmCouplingMom = 0.0 -- for timer 

   if self.positivity then 
      self.positivityRescale = Updater.PositivityRescale {
         onGrid = self.grid,
         basis = self.basis,
      }
   end
end

function GkSpecies:forwardEuler(tCurr, dt, species, emIn, inIdx, outIdx)
   local fIn = self:rkStepperFields()[inIdx]
   local fOut = self:rkStepperFields()[outIdx]

   local em = emIn[1]:rkStepperFields()[inIdx]
   local emFunc = emIn[2]:rkStepperFields()[1]

   local status, dtSuggested = true, GKYL_MAX_DOUBLE


   if self.evolveCollisionless then

      if self.positivity then 
         self.positivityRescale:advance(tCurr, dt, {fIn}, {self.fPos}) 
         self:applyBc(tCurr, dt, self.fPos)
         status, dtSuggested = self.solver:advance(tCurr, dt, {self.fPos, em, emFunc, emGy}, {fOut})
      else
         status, dtSuggested = self.solver:advance(tCurr, dt, {fIn, em, emFunc, emGy}, {fOut})
      end
      -- if step2, only compute RHS increment so that RHS can be used in step 2
      if not self.solverStep2 then fOut:scale(dt); fOut:accumulate(1.0, fIn) end -- fOut = fIn + dt*fOut
   else
      fOut:copy(fIn) -- just copy stuff over
      self.gkEqn:setAuxFields({em, emFunc})  -- set auxFields in case they are needed by BCs/collisions
   end

   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
	 local collStatus, collDt = c:forwardEuler(
	    tCurr, dt, fIn, species, fOut)
	 -- the full 'species' list is needed for the cross-species
	 -- collisions
	 status = status and collStatus
	 dtSuggested = math.min(dtSuggested, collDt)
      end
   end

   if self.sourceFunc and self.evolveSources then
     -- if there is a source, add it to the RHS
     fOut:accumulate(dt*self.sourceTimeDependence(tCurr), self.fSource)
   end

   -- apply BCs
   self:applyBc(tCurr, dt, fOut)
   
   return status, dtSuggested
end

function GkSpecies:forwardEulerStep2(tCurr, dt, species, emIn, inIdx, outIdx)
   local fIn = self:rkStepperFields()[inIdx]
   local fOut = self:rkStepperFields()[outIdx]

   if self.evolve then
      local em = emIn[1]:rkStepperFields()[inIdx]
      local emFunc = emIn[2]:rkStepperFields()[1]
      local status, dtSuggested
      status, dtSuggested = self.solverStep2:advance(tCurr, dt, {fIn, em, emFunc}, {fOut})
      if self.positivity then 
         self.positivityRescale:advance(tCurr, dt, {fIn}, {self.fPos}) 
         status, dtSuggested = self.solverStep2:advance(tCurr, dt, {self.fPos, em, emFunc}, {fOut})
      else
         status, dtSuggested = self.solverStep2:advance(tCurr, dt, {fIn, em, emFunc}, {fOut})
      end
      fOut:scale(dt); fOut:accumulate(1.0, fIn) -- fOut = fIn + dt*fOut

      -- apply BCs
      self:applyBc(tCurr, dt, fOut)

      return status, dtSuggested
   else
      fOut:copy(fIn) -- just copy stuff over
      return true, GKYL_MAX_DOUBLE
   end
end

function GkSpecies:createDiagnostics()
   -- create updater to compute volume-integrated moments -- NOT YET IMPLEMENTED FOR GK
   self.intMomentCalc = nil
   
   -- function to check if moment name is correct
   local function isMomentNameGood(nm)
      return Updater.DistFuncMomentCalc:isGkMomentNameGood(nm)
   end

   self.diagnosticMomentFields = { }
   self.diagnosticMomentUpdaters = { } 
   -- allocate space to store moments and create moment updater
   for i, mom in pairs(self.diagnosticMoments) do
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
            gkfacs = {self.mass, self.bmag},
         }
      else
         assert(false, string.format("Moment %s not valid", mom))
      end
   end
end

-- BC functions
function GkSpecies:bcReflectFunc(dir, tm, idxIn, fIn, fOut)
   -- skinLoop should be "flip"
   -- note that GK reflection only valid in z-vpar.
   -- this is checked when bc is created.

   self.basis:flipSign(dir, fIn, fOut)
   -- vpar is always first velocity dimension
   local vpardir=self.cdim+1 
   self.basis:flipSign(vpardir, fOut, fOut)
end
function GkSpecies:bcSheathFunc(dir, tm, idxIn, fIn, fOut)
   -- skinLoop should be "flip"
   -- note that GK reflection only valid in z-vpar.
   -- this is checked when bc is created.

   -- need to figure out if we are on lower or upper domain edge
   local edgeVal
   local globalRange = self.grid:globalRange()
   if idxIn[dir] == globalRange:lower(dir) then 
      -- this means we are at lower domain edge, 
      -- so we need to evaluate basis functions at z=-1
      edgeVal = -1 
   else 
      -- this means we are at upper domain edge
      -- so we need to evaluate basis functions at z=1
      edgeVal = 1 
   end
   local gkEqn = self.gkEqn
   -- calculate deltaPhi = phi - phiWall
   -- note: this gives surface-averaged scalar value of deltaPhi in this cell
   local deltaPhi = gkEqn:calcSheathDeltaPhi(idxIn, edgeVal)

   -- get vpar limits of cell
   local vpardir = self.cdim+1
   local gridIn = self.grid
   gridIn:setIndex(idxIn)
   local vL = gridIn:cellLowerInDir(vpardir)
   local vR = gridIn:cellUpperInDir(vpardir)
   local vlower, vupper
   -- this makes it so that we only need to deal with absolute values of vpar
   if math.abs(vR)>=math.abs(vL) then
      vlower = math.abs(vL)
      vupper = math.abs(vR)
   else
      vlower = math.abs(vR)
      vupper = math.abs(vL)
   end
   if -self.charge*deltaPhi > 0 then
      -- calculate cutoff velocity for reflection
      local vcut = math.sqrt(-2*self.charge*deltaPhi/self.mass)
      if vcut > vupper then
         -- reflect if vcut is above the velocities in this cell
         self:bcReflectFunc(dir, tm, nil, fIn, fOut)
      elseif vcut > vlower then
          -- partial reflect if vcut is in this velocity cell
          local fhat = self.fhatSheathPtr
          self.fhatSheath:fill(self.fhatSheathIdxr(idxIn), fhat)
          local w = gridIn:cellCenterInDir(vpardir)
          local dv = gridIn:dx(vpardir)
          -- calculate weak-equivalent distribution fhat
          gkEqn:calcSheathPartialReflection(w, dv, edgeVal, vcut, fIn, fhat)
          -- reflect fhat into skin cells
          self:bcReflectFunc(dir, tm, nil, fhat, fOut) 
      else
         -- absorb if vcut is below the velocities in this cell
         self:bcAbsorbFunc(dir, tm, nil, fIn, fOut)
      end
   else 
      -- entire species (usually ions) is lost
      self:bcAbsorbFunc(dir, tm, nil, fIn, fOut)
   end
end

function GkSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- need to wrap member functions so that self is passed
   local function bcAbsorbFunc(...) return self:bcAbsorbFunc(...) end
   local function bcOpenFunc(...) return  self:bcOpenFunc(...) end
   local function bcReflectFunc(...) return self:bcReflectFunc(...) end
   local function bcSheathFunc(...) return self:bcSheathFunc(...) end
   
   local vdir = nil
   if dir==self.cdim then 
      vdir=self.cdim+1 
   end

   if bcType == SP_BC_ABSORB then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcAbsorbFunc }, "pointwise"))
   elseif bcType == SP_BC_OPEN then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcOpenFunc }, "pointwise"))
   -- note: reflection and sheath BCs only make sense in z direction,
   -- which is always last config space direction, i.e. dir = self.cdim
   elseif bcType == SP_BC_REFLECT and dir==self.cdim then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcReflectFunc }, "flip"))
   elseif bcType == SP_BC_SHEATH and dir==self.cdim then
      self.fhatSheath = self:allocDistf()
      self.fhatSheathPtr = self.fhatSheath:get(1)
      self.fhatSheathIdxr = self.fhatSheath:genIndexer()
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcSheathFunc }, "flip"))
      self.hasSheathBcs = true
   elseif bcType == SP_BC_ZEROFLUX then
      table.insert(self.zeroFluxDirections, dir)
   else
      assert(false, "GkSpecies: Unsupported BC type!")
   end
end

function GkSpecies:calcCouplingMoments(tCurr, dt, rkIdx)
   local fIn = self:rkStepperFields()[rkIdx]

   -- compute moments needed in coupling to fields and collisions
   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      
      if self.collisions then 
         self.threeMomentsCalc:advance(tCurr, dt, {fIn}, { self.numDensity, self.momDensity, self.ptclEnergy })
      else
         self.numDensityCalc:advance(tCurr, dt, {fIn}, { self.numDensity })
      end

      self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
end

function GkSpecies:fluidMoments()
   return { self.numDensity, self.momDensity, self.ptclEnergy } 
end

function GkSpecies:getNumDensity(rkIdx)
   -- if no rkIdx specified, assume numDensity has already been calculated
   if rkIdx == nil then return self.numDensity end 

   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      self.numDensityCalc:advance(tCurr, dt, {self:rkStepperFields()[rkIdx]}, { self.numDensity })
      self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
   return self.numDensity
end

function GkSpecies:getBackgroundDens()
   return self.n0
end

function GkSpecies:getMomDensity(rkIdx)
   -- if no rkIdx specified, assume momDensity has already been calculated
   if rkIdx == nil then return self.momDensity end 
 
   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      self.momDensityCalc:advance(tCurr, dt, {self:rkStepperFields()[rkIdx]}, { self.momDensity })
      self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
   return self.momDensity
end

function GkSpecies:getPolarizationWeight(linearized)
   if linearized == false then 
     self.polarizationWeight:combine(self.mass/self.B0^2, self.numDensity)
     return self.polarizationWeight
   else 
     return self.n0*self.mass/self.B0^2
   end
end

function GkSpecies:momCalcTime()
   local tm = self.tmCouplingMom
   for i, mom in ipairs(self.diagnosticMoments) do
      tm = tm + self.diagnosticMomentUpdaters[i].totalTime
   end
   return tm
end

function GkSpecies:solverVolTime()
   return self.gkEqn.totalVolTime
end

function GkSpecies:solverSurfTime()
   return self.gkEqn.totalSurfTime
end
function GkSpecies:totalSolverTime()
   local timer = self.solver.totalTime
   if self.solverStep2 then timer = timer + self.solverStep2.totalTime end
   if self.positivityRescale then timer = timer + self.positivityRescale.totalTime end
   return timer
end

function GkSpecies:Maxwellian(xn, n0, T0, vd)
   local vd = vd or 0.0
   local vt2 = T0/self.mass
   local vpar = xn[self.cdim+1]
   local v2 = (vpar-vd)^2
   if self.vdim > 1 then 
     local mu = xn[self.cdim+2]
     v2 = v2 + 2*math.abs(mu)*self.bmagFunc(0,xn)/self.mass
     return n0*(2*math.pi*vt2)^(-3/2)*math.exp(-v2/(2*vt2))
   else
     return n0*(2*math.pi*vt2)^(-1/2)*math.exp(-v2/(2*vt2))
   end
end

return GkSpecies
