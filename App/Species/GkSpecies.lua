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
local Constants = require "Lib.Constants"

local GkSpecies = Proto(KineticSpecies)

-- add constants to object indicate various supported boundary conditions
local SP_BC_ABSORB = 1
local SP_BC_OPEN = 2
local SP_BC_REFLECT = 3
local SP_BC_SHEATH = 4
local SP_BC_ZEROFLUX = 5
local SP_BC_COPY = 6
GkSpecies.bcAbsorb = SP_BC_ABSORB -- absorb all particles
GkSpecies.bcOpen = SP_BC_OPEN -- zero gradient
GkSpecies.bcReflect = SP_BC_REFLECT -- specular reflection
GkSpecies.bcSheath = SP_BC_SHEATH -- specular reflection
GkSpecies.bcZeroFlux = SP_BC_ZEROFLUX -- zero flux
GkSpecies.bcCopy = SP_BC_COPY -- copy stuff

function GkSpecies:alloc(nRkDup)
   -- allocate distribution function
   GkSpecies.super.alloc(self, nRkDup)

   -- allocate fields to store coupling moments (for use in coupling
   -- to field and collisions)
   self.numDensity = self:allocMoment()
   self.numDensityAux = self:allocMoment()
   self.momDensity = self:allocMoment()
   self.momDensityAux = self:allocMoment()
   self.ptclEnergy = self:allocMoment()
   self.polarizationWeight = self:allocMoment() -- not used when using linearized poisson solve

   if self.gyavg then
      self.rho1 = self:allocDistf()
      self.rho2 = self:allocDistf()
      self.rho3 = self:allocDistf()
   end
end

function GkSpecies:allocMomCouplingFields()
   assert(false, "GkSpecies:allocMomCouplingFields should not be called. Field object should allocate its own coupling fields")
end

function GkSpecies:createSolver(hasPhi, hasApar, funcField)
   -- run the KineticSpecies 'createSolver()' to initialize the
   -- collisions solver
   GkSpecies.super.createSolver(self,funcField)

   -- set up jacobian
   if funcField then
      -- save bmagFunc for later...
      self.bmagFunc = funcField.bmagFunc
      -- if vdim>1, get jacobian=bmag from geo
      self.jacobPhaseFunc = self.bmagFunc
      self.jacobGeoFunc = funcField.jacobGeoFunc
      if self.cdim == 1 then 
         self.B0 = funcField.bmagFunc(0.0, {self.grid:mid(1)})
      elseif self.cdim == 2 then 
         self.B0 = funcField.bmagFunc(0.0, {self.grid:mid(1), self.grid:mid(1)})
      else
         self.B0 = funcField.bmagFunc(0.0, {self.grid:mid(1), self.grid:mid(1), self.grid:mid(2)})
      end
      self.bmag = assert(funcField.geo.bmag, "nil bmag")
      self.bmagInv = funcField.geo.bmagInv
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
      project1:advance(0.0, {}, {self.rho1})
      local project2 = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = rho2Func,
         projectOnGhosts = true
      }
      project2:advance(0.0, {}, {self.rho2})
      local project3 = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = rho3Func,
         projectOnGhosts = true
      }
      project3:advance(0.0, {}, {self.rho3})

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
   if hasApar and self.basis:polyOrder() > 1 then -- if electromagnetic only update conf dir surface terms on first step
      for d = 1, self.cdim do upd[d] = d end
   else
      for d = 1, self.cdim + 1 do upd[d] = d end
   end
   -- zero flux in vpar and mu
   table.insert(self.zeroFluxDirections, self.cdim+1)
   if self.vdim > 1 then table.insert(self.zeroFluxDirections, self.cdim+2) end

   self.solver = Updater.HyperDisCont {
      onGrid = self.grid,
      basis = self.basis,
      cfl = self.cfl,
      equation = self.gkEqn,
      zeroFluxDirections = self.zeroFluxDirections,
      updateDirections = upd,
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
      moment = "GkM1proj",
      gkfacs = {self.mass, self.bmag},
   }
   self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkM2",
      gkfacs = {self.mass, self.bmag},
   }
   self.M2parCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkM2par",
      gkfacs = {self.mass, self.bmag},
   }
   if self.vdim > 1 then
      self.M2perpCalc = Updater.DistFuncMomentCalc {
         onGrid = self.grid,
         phaseBasis = self.basis,
         confBasis = self.confBasis,
         moment = "GkM2perp",
         gkfacs = {self.mass, self.bmag},
      }
   end
   self.threeMomentsCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "GkThreeMoments",
      gkfacs = {self.mass, self.bmag},
   }
   
   self._firstMomentCalc = true  -- to avoid re-calculating moments when not evolving

   self.tmCouplingMom = 0.0 -- for timer 

   if self.positivityRescale or self.positivityDiffuse then 
      self.posRescaler = Updater.PositivityRescale {
         onGrid = self.grid,
         basis = self.basis,
      }
   end

   assert(self.n0, "Must specify background density as global variable 'n0' in species table as 'n0 = ...'")
end

function GkSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self.tCurr = tCurr
   local fIn = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em = emIn[1]:rkStepperFields()[inIdx]
   local dApardtPrev = emIn[1].dApardtPrev
   local emFunc = emIn[2]:rkStepperFields()[1]

   -- rescale slopes
   if self.positivityRescale then
      self.posRescaler:advance(tCurr, {fIn}, {self.fPos})
      fIn = self.fPos
   end

   -- do collisions first so that collisions contribution to cflRate is included in GK positivity
   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
         c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
         c:advance(tCurr, fIn, species, fRhsOut)
         -- the full 'species' list is needed for the cross-species
         -- collisions
      end
   end
   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solver:advance(tCurr, {fIn, em, emFunc, dApardtPrev}, {fRhsOut})
   else
      fRhsOut:clear(0.0) -- no RHS
      self.gkEqn:setAuxFields({em, emFunc, dApardtPrev})  -- set auxFields in case they are needed by BCs/collisions
   end

   if not self.solverStep2 then -- if step2, wait to do sources
      if self.fSource and self.evolveSources then
         -- add source it to the RHS
         fRhsOut:accumulate(self.sourceTimeDependence(tCurr), self.fSource)
      end
   end
end

function GkSpecies:advanceStep2(tCurr, species, emIn, inIdx, outIdx)
   local fIn = self:rkStepperFields()[inIdx]
   if self.positivityRescale then
      fIn = self.fPos
   end
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em = emIn[1]:rkStepperFields()[inIdx]
   local dApardtPrev = emIn[1].dApardtPrev
   local emFunc = emIn[2]:rkStepperFields()[1]

   if self.evolveCollisionless then
      self.solverStep2:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solverStep2:advance(tCurr, {fIn, em, emFunc, dApardtPrev}, {fRhsOut})
   else
      fRhsOut:clear(0.0)  -- no RHS
   end

   if self.fSource and self.evolveSources then
      -- add source it to the RHS
      fRhsOut:accumulate(self.sourceTimeDependence(tCurr), self.fSource)
   end
end

function GkSpecies:createDiagnostics()
   -- create updater to compute volume-integrated moments -- NOT YET IMPLEMENTED FOR GK
   self.intMomentCalc = nil
   
   -- function to check if moment name is correct
   local function isMomentNameGood(nm)
      return Updater.DistFuncMomentCalc:isGkMomentNameGood(nm)
   end
   local function isWeakMomentNameGood(nm)
      return nm == "GkUpar" or nm == "GkTpar" or nm == "GkTperp" or nm == "GkTemp"
   end
   local function isAuxMomentNameGood(nm)
      return nm == "GkBeta"
   end
   local function contains(table, element)
     for _, value in pairs(table) do
       if value == element then
         return true
       end
     end
     return false
   end

   self.diagnosticMomentFields = { }
   self.diagnosticMomentUpdaters = { } 
   self.diagnosticWeakMoments = { }
   self.diagnosticAuxMoments = { }
   self.weakMomentOpFields = { }
   self.weakMomentScaleFac = { }
   -- set up weak multiplication and division operators
   self.weakMultiplication = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Multiply",
      onGhosts = true,
   }
   self.weakDivision = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
      onGhosts = true,
   }
   -- sort moments into diagnosticMoments, diagnosticWeakMoments, and diagnosticAuxMoments
   for i, mom in pairs(self.diagnosticMoments) do
      if isWeakMomentNameGood(mom) then
         -- remove moment name from self.diagnosticMoments list, and add it to self.diagnosticWeakMoments list
         table.insert(self.diagnosticWeakMoments, mom)
         self.diagnosticMoments[i] = nil
      elseif isAuxMomentNameGood(mom) then
         -- remove moment name from self.diagnosticMoments list, and add it to self.diagnosticAuxMoments list
         table.insert(self.diagnosticAuxMoments, mom)
         self.diagnosticMoments[i] = nil
      end
   end

   -- make sure we have the updaters needed to calculate all the weak and aux moments
   for i, mom in pairs(self.diagnosticAuxMoments) do
      if mom == "GkBeta" then
         if not contains(self.diagnosticWeakMoments, "GkTemp") then 
            table.insert(self.diagnosticWeakMoments, "GkTemp")
         end
         if not contains(self.diagnosticMoments, "GkM0") then
            table.insert(self.diagnosticMoments, "GkM0")
         end
      end
   end

   for i, mom in pairs(self.diagnosticWeakMoments) do
      -- all GK weak moments require M0 = density
      if not contains(self.diagnosticMoments, "GkM0") then
         table.insert(self.diagnosticMoments, "GkM0")
      end

      if mom == "GkUpar" then
         if not contains(self.diagnosticMoments, "GkM1") then
            table.insert(self.diagnosticMoments, "GkM1")
         end
      end
      if mom == "GkTpar" then
         if not contains(self.diagnosticMoments, "GkM2par") then
            table.insert(self.diagnosticMoments, "GkM2par")
         end
         if not contains(self.diagnosticWeakMoments, "GkUpar") then
            table.insert(self.diagnosticWeakMoments, "GkUpar")
         end
      elseif mom == "GkTperp" then
         if not contains(self.diagnosticMoments, "GkM2perp") then
            table.insert(self.diagnosticMoments, "GkM2perp")
         end
      elseif mom == "GkTemp" then 
         if not contains(self.diagnosticMoments, "GkM2") then
            table.insert(self.diagnosticMoments, "GkM2")
         end      
         if not contains(self.diagnosticWeakMoments, "GkUpar") then
            table.insert(self.diagnosticWeakMoments, "GkUpar")
         end
      end
   end

   -- allocate space to store moments and create moment updaters
   for i, mom in pairs(self.diagnosticMoments) do
      if isMomentNameGood(mom) then
         self.diagnosticMomentFields[mom] = DataStruct.Field {
            onGrid = self.confGrid,
            numComponents = self.confBasis:numBasis(),
            ghost = {1, 1}
         }
         self.diagnosticMomentUpdaters[mom] = Updater.DistFuncMomentCalc {
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
   for i, mom in pairs(self.diagnosticWeakMoments) do
      if isWeakMomentNameGood(mom) then
         self.diagnosticMomentFields[mom] = DataStruct.Field {
            onGrid = self.confGrid,
            numComponents = self.confBasis:numBasis(),
            ghost = {1, 1}
         }
      else
         assert(false, string.format("Moment %s not valid", mom))
      end

      if mom == "GkUpar" then
         self.weakMomentOpFields["GkUpar"] = {self.diagnosticMomentFields["GkM0"], self.diagnosticMomentFields["GkM1"]}
      elseif mom == "GkTpar" then
         self.weakMomentOpFields["GkTpar"] = {self.diagnosticMomentFields["GkM0"], self.diagnosticMomentFields["GkM2par"]}
         self.weakMomentScaleFac["GkTpar"] = self.mass
      elseif mom == "GkTperp" then
         self.weakMomentOpFields["GkTperp"] = {self.diagnosticMomentFields["GkM0"], self.diagnosticMomentFields["GkM2perp"]}
         self.weakMomentScaleFac["GkTperp"] = self.mass
      elseif mom == "GkTemp" then 
         self.weakMomentOpFields["GkTemp"] = {self.diagnosticMomentFields["GkM0"], self.diagnosticMomentFields["GkM2"]}
         self.weakMomentScaleFac["GkTemp"] = self.mass/3
      end
   end
   for i, mom in pairs(self.diagnosticAuxMoments) do
      if isAuxMomentNameGood(mom) then
         self.diagnosticMomentFields[mom] = DataStruct.Field {
            onGrid = self.confGrid,
            numComponents = self.confBasis:numBasis(),
            ghost = {1, 1}
         }
      else
         assert(false, string.format("Moment %s not valid", mom))
      end
   end
end

function GkSpecies:calcDiagnosticWeakMoments()
   GkSpecies.super.calcDiagnosticWeakMoments(self)
   -- need to subtract m*Upar^2 from GkTemp and GkTpar
   if self.diagnosticWeakMoments["GkTemp"] or self.diagnosticWeakMoments["GkTpar"] then
      self.weakMultiplication:advance(0.0,
           {self.diagnosticMomentFields["GkUpar"], self.diagnosticMomentFields["GkUpar"]}, 
           {self.momDensityAux})
   end
   if self.diagnosticWeakMoments["GkTemp"] then
      self.diagnosticWeakMoments["GkTemp"]:accumulate(-self.mass/3, self.momDensityAux)
   end
   if self.diagnosticWeakMoments["GkTpar"] then
      self.diagnosticWeakMoments["GkTpar"]:accumulate(-self.mass, self.momDensityAux)
   end
end

function GkSpecies:calcDiagnosticAuxMoments()
   if self.diagnosticMomentFields["GkBeta"] then
      self.weakMultiplication:advance(0.0, 
           {self.diagnosticMomentFields["GkM0"], self.diagnosticMomentFields["GkTemp"]}, 
           {self.diagnosticMomentFields["GkBeta"]})
      self.weakMultiplication:advance(0.0, 
           {self.diagnosticMomentFields["GkBeta"], self.bmagInv}, 
           {self.diagnosticMomentFields["GkBeta"]})
      self.weakMultiplication:advance(0.0, 
           {self.diagnosticMomentFields["GkBeta"], self.bmagInv}, 
           {self.diagnosticMomentFields["GkBeta"]})
      self.diagnosticMomentFields["GkBeta"]:scale(2*Constants.MU0)
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
   local function bcCopyFunc(...) return self:bcCopyFunc(...) end
   
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
   elseif bcType == SP_BC_COPY then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcCopyFunc }, "pointwise"))
   else
      assert(false, "GkSpecies: Unsupported BC type!")
   end
end

function GkSpecies:calcCouplingMoments(tCurr, rkIdx)
   local fIn = self:rkStepperFields()[rkIdx]

   -- compute moments needed in coupling to fields and collisions
   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()

      if self.deltaF then
        fIn:accumulate(-1.0, self.f0)
      end
      
      if self.collisions then 
         self.threeMomentsCalc:advance(tCurr, {fIn}, { self.numDensity, self.momDensity, self.ptclEnergy })
      else
         self.numDensityCalc:advance(tCurr, {fIn}, { self.numDensity })
      end

      if self.deltaF then
        fIn:accumulate(1.0, self.f0)
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
   local fIn = self:rkStepperFields()[rkIdx]

   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      if self.deltaF then
        fIn:accumulate(-1.0, self.f0)
      end
      self.numDensityCalc:advance(nil, {fIn}, { self.numDensityAux })
      if self.deltaF then
        fIn:accumulate(1.0, self.f0)
      end
      self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
   return self.numDensityAux
end

function GkSpecies:getBackgroundDens()
   return self.n0
end

function GkSpecies:getMomDensity(rkIdx)
   -- if no rkIdx specified, assume momDensity has already been calculated
   if rkIdx == nil then return self.momDensity end 
   local fIn = self:rkStepperFields()[rkIdx]
 
   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      if self.deltaF then
        fIn:accumulate(-1.0, self.f0)
      end
      self.momDensityCalc:advance(nil, {fIn}, { self.momDensityAux })
      if self.deltaF then
        fIn:accumulate(1.0, self.f0)
      end
      self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
   return self.momDensityAux
end

function GkSpecies:getEmModifier(rkIdx)
   -- for p > 1, this is just numDensity
   if self.basis:polyOrder() > 1 then return self:getNumDensity(rkIdx) end

   local fIn = self.gkEqn.emMod

   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      if self.deltaF then
        fIn:accumulate(-1.0, self.f0)
      end
      self.momDensityCalc:advance(nil, {fIn}, { self.momDensityAux })
      if self.deltaF then
        fIn:accumulate(1.0, self.f0)
      end
      self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
   fIn:clear(0.0)
   return self.momDensityAux
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
      tm = tm + self.diagnosticMomentUpdaters[mom].totalTime
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
   if self.posRescaler then timer = timer + self.posRescaler.totalTime end
   return timer
end

function GkSpecies:Maxwellian(xn, n0, T0, vdIn)
   local vd = vdIn or 0.0
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
