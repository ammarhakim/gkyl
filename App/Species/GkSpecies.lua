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

function GkSpecies:initDist(geo)
   if geo then
      -- save bmagFunc for later...
      self.bmagFunc = geo.bmagFunc
      -- if vdim>1, get jacobian=bmag from geo, and multiply it by init functions for f0 and f
      self.jacobianFunc = self.bmagFunc
      if self.jacobianFunc and self.vdim > 1 then
         local initFuncWithoutJacobian = self.initFunc
         self.initFunc = function (t, xn)
            local J = self.jacobianFunc(t,xn)
            local f = initFuncWithoutJacobian(t,xn)
            return J*f
         end
         if self.initBackgroundFunc then
            local initBackgroundFuncWithoutJacobian = self.initBackgroundFunc
            self.initBackgroundFunc = function(t,xn)
               local J = self.jacobianFunc(t,xn)
               local f0 = initBackgroundFuncWithoutJacobian(t,xn)
               return J*f0
            end
         end
      end
   end

   GkSpecies.super.initDist(self)
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
   -- except in 1v there is no mu, so no factor of mass and also remove the factor of 2*pi... 
   if self.vdim == 1 then self.momfac = 1 / (2*math.pi) end
   
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
   self._firstMomentCalc = true  -- to avoid re-calculating moments when not evolving
 
   -- create updater to evaluate source 
   if self.sourceFunc then 
      self.evalSource = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = self.sourceFunc,
         projectOnGhosts = true
      }
   end

   -- calculate background density averaged over simulation domain
   self.n0 = nil
   if self.f0 then
      self.calcDens:advance(0,0, {self.f0}, {self.dens0})
      local data
      local dynVec = DataStruct.DynVector { numComponents = 1 }
      -- integrate 
      local calcInt = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.confGrid,
         basis = self.confBasis,
         numComponents = 1,
	 quantity = "V"
      }
      calcInt:advance(0.0, 0.0, {self.dens0}, {dynVec})
      _, data = dynVec:lastData()
      self.n0 = data[1]/self.confGrid:gridVolume()
      --print("Average density is " .. self.n0)
   end
end

function GkSpecies:forwardEuler(tCurr, dt, fIn, emIn, fOut)
   if self.evolve then
      local em = emIn[1]
      local emFunc = emIn[2]
      local status, dtSuggested
      status, dtSuggested = self.solver:advance(tCurr, dt, {fIn, em, emFunc}, {fOut})
      if self.sourceFunc then
        -- if there is a source, add it to the RHS
        local fSource = self.fSource
        self.evalSource:advance(tCurr, dt, {}, {fSource})
        fOut:accumulate(dt, fSource)
      end
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
   --local global = fOut:globalRange()
   local edgeVal
   if idxIn[dir] <= 2 then -- HACK! == global:lower(dir) then 
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
   else
      assert(false, "GkSpecies: Unsupported BC type!")
   end
end

function GkSpecies:calcCouplingMoments(tCurr, dt, fIn)
   -- compute moments needed in coupling to fields and collisions
   if self.evolve or self._firstMomentCalc then
      -- note: this is a hack until inaccuracy in initial conditions worked out
      if self.f0 then self.distf[1]:accumulate(-1, self.f0) end  --
      self.calcDens:advance(tCurr, dt, {fIn}, { self.dens })
      self.calcUpar:advance(tCurr, dt, {fIn}, { self.upar })
      if self.f0 then self.distf[1]:accumulate(1, self.f0) end  --
   end
   if not self.evolve then self._firstMomentCalc = false end
end

function GkSpecies:calcDiagnosticMoments()
   -- if there is a background distribution f0, we will output only
   -- perturbed moments by subtracting off background
   if self.f0 then self.distf[1]:accumulate(-1, self.f0) end
   local numMoms = #self.diagnosticMoments
   for i = 1, numMoms do
      self.diagnosticMomentUpdaters[i]:advance(
	 0.0, 0.0, {self.distf[1]}, {self.diagnosticMomentFields[i]})
   end
   if self.f0 then self.distf[1]:accumulate(1, self.f0) end
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
