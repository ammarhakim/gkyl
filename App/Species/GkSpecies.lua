-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic species object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto          = require "Lib.Proto"
local KineticSpecies = require "App.Species.KineticSpecies"
local Mpi            = require "Comm.Mpi"
local Gk             = require "Eq.Gyrokinetic"
local Updater        = require "Updater"
local DataStruct     = require "DataStruct"
local Time           = require "Lib.Time"
local Constants      = require "Lib.Constants"

local GkSpecies = Proto(KineticSpecies)

-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB   = 1
local SP_BC_OPEN     = 2
local SP_BC_REFLECT  = 3
local SP_BC_SHEATH   = 4
local SP_BC_ZEROFLUX = 5
local SP_BC_COPY     = 6
GkSpecies.bcAbsorb   = SP_BC_ABSORB      -- Absorb all particles.
GkSpecies.bcOpen     = SP_BC_OPEN        -- Zero gradient.
GkSpecies.bcReflect  = SP_BC_REFLECT     -- Specular reflection.
GkSpecies.bcSheath   = SP_BC_SHEATH      -- Specular reflection.
GkSpecies.bcZeroFlux = SP_BC_ZEROFLUX    -- Zero flux.
GkSpecies.bcCopy     = SP_BC_COPY        -- Copy stuff.

function GkSpecies:alloc(nRkDup)
   -- Allocate distribution function.
   GkSpecies.super.alloc(self, nRkDup)

   -- Allocate fields to store coupling moments (for use in coupling to field and collisions).
   self.numDensity         = self:allocMoment()
   self.momDensity         = self:allocMoment()
   self.ptclEnergy         = self:allocMoment()
   self.numDensityAux      = self:allocMoment()
   self.momDensityAux      = self:allocMoment()
   self.ptclEnergyAux      = self:allocMoment()
   if self.positivity then
      self.numDensityPos   = self:allocMoment()
      self.momDensityPos   = self:allocMoment()
      self.ptclEnergyPos   = self:allocMoment()
   end
   self.polarizationWeight = self:allocMoment() -- not used when using linearized poisson solve

   if self.gyavg then
      self.rho1 = self:allocDistf()
      self.rho2 = self:allocDistf()
      self.rho3 = self:allocDistf()
   end

   if self.vdim == 1 then
      self.vDegFreedom = 1.0
   else
      self.vDegFreedom = 3.0
   end

   self.first = true
end

function GkSpecies:allocMomCouplingFields()
   assert(false, "GkSpecies:allocMomCouplingFields should not be called. Field object should allocate its own coupling fields")
end

function GkSpecies:createSolver(hasPhi, hasApar, funcField)
   -- Run the KineticSpecies 'createSolver()' to initialize the
   -- collisions solver.
   GkSpecies.super.createSolver(self,funcField)

   -- Set up jacobian.
   if funcField then
      -- Save bmagFunc for later...
      self.bmagFunc = funcField.bmagFunc
      -- If vdim>1, get jacobian=bmag from geo.
      self.jacobPhaseFunc = self.bmagFunc
      self.jacobGeoFunc   = funcField.jacobGeoFunc
      if self.cdim == 1 then 
         self.B0 = funcField.bmagFunc(0.0, {self.grid:mid(1)})
      elseif self.cdim == 2 then 
         self.B0 = funcField.bmagFunc(0.0, {self.grid:mid(1), self.grid:mid(2)})
      else
         self.B0 = funcField.bmagFunc(0.0, {self.grid:mid(1), self.grid:mid(2), self.grid:mid(3)})
      end
      self.bmag    = assert(funcField.geo.bmag, "nil bmag")
      self.bmagInv = funcField.geo.bmagInv
   end

   if self.gyavg then
      -- Set up geo fields needed for gyroaveraging.
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
         onGrid          = self.grid,
         basis           = self.basis,
         evaluate        = rho1Func,
         projectOnGhosts = true
      }
      project1:advance(0.0, {}, {self.rho1})
      local project2 = Updater.ProjectOnBasis {
         onGrid          = self.grid,
         basis           = self.basis,
         evaluate        = rho2Func,
         projectOnGhosts = true
      }
      project2:advance(0.0, {}, {self.rho2})
      local project3 = Updater.ProjectOnBasis {
         onGrid          = self.grid,
         basis           = self.basis,
         evaluate        = rho3Func,
         projectOnGhosts = true
      }
      project3:advance(0.0, {}, {self.rho3})

      -- Create solver for gyroaveraging potentials.
      self.emGyavgSlvr = Updater.FemGyroaverage {
         onGrid     = self.confGrid,
         confBasis  = self.confBasis,
         phaseGrid  = self.grid,
         phaseBasis = self.basis,
         rho1       = self.rho1,
         rho2       = self.rho2,
         rho3       = self.rho3,
         muOrder0   = true, -- Cell-average in mu.
      }

      -- Create solver for gyroaveraging distribution function.
      self.distfGyavgSlvr = Updater.FemGyroaverage {
         onGrid     = self.confGrid,
         confBasis  = self.confBasis,
         phaseGrid  = self.grid,
         phaseBasis = self.basis,
         rho1       = self.rho1,
         rho2       = self.rho2,
         rho3       = self.rho3,
         integrate  = true,
      }
   end

   -- Create updater to advance solution by one time-step.
   self.gkEqn = Gk.GkEq {
      onGrid       = self.grid,
      confGrid     = self.confGrid,
      phaseBasis   = self.basis,
      confBasis    = self.confBasis,
      charge       = self.charge,
      mass         = self.mass,
      hasPhi       = hasPhi,
      hasApar      = hasApar,
      Bvars        = funcField.bmagVars,
      hasSheathBcs = self.hasSheathBcs,
      positivity   = self.positivity,
      gyavgSlvr    = self.emGyavgSlvr,
   }

   -- No update in mu direction (last velocity direction if present)
   local upd = {}
   if hasApar then    -- If electromagnetic only update conf dir surface terms on first step.
      for d = 1, self.cdim do upd[d] = d end
   else
      for d = 1, self.cdim + 1 do upd[d] = d end
   end
   -- Zero flux in vpar and mu.
   table.insert(self.zeroFluxDirections, self.cdim+1)
   if self.vdim > 1 then table.insert(self.zeroFluxDirections, self.cdim+2) end

   self.solver = Updater.HyperDisCont {
      onGrid             = self.grid,
      basis              = self.basis,
      cfl                = self.cfl,
      equation           = self.gkEqn,
      zeroFluxDirections = self.zeroFluxDirections,
      updateDirections   = upd,
      clearOut           = false,   -- Continue accumulating into output field.
   }
   if hasApar and self.basis:polyOrder()==1 then 
      -- This solver calculates vpar surface terms for Ohm's law. p=1 only!
      self.solverStep2 = Updater.HyperDisCont {
         onGrid             = self.grid,
         basis              = self.basis,
         cfl                = self.cfl,
         equation           = self.gkEqn,
         zeroFluxDirections = self.zeroFluxDirections,
         updateDirections   = {self.cdim+1},    -- Only vpar terms.
         updateVolumeTerm   = false,            -- No volume term.
         clearOut           = false,            -- Continue accumulating into output field.
      }
      -- Note that the surface update for this term only involves the vpar direction.
      self.solverStep3 = Updater.HyperDisCont {
         onGrid             = self.grid,
         basis              = self.basis,
         cfl                = self.cfl,
         equation           = self.gkEqn,
         eqnStep            = 2,  -- Use step2 functions from gkEqn.
         zeroFluxDirections = self.zeroFluxDirections,
         updateSurfaceTerm  = false,
         clearOut           = false,   -- Continue accumulating into output field.
      }
   elseif hasApar and self.basis:polyOrder()>1 then
      -- Set up solver that adds on volume term involving dApar/dt and the entire vpar surface term.
      -- Note that the surface update for this term only involves the vpar direction.
      self.solverStep2 = Updater.HyperDisCont {
         onGrid             = self.grid,
         basis              = self.basis,
         cfl                = self.cfl,
         equation           = self.gkEqn,
         eqnStep            = 2,  -- Use step2 functions from gkEqn.
         zeroFluxDirections = self.zeroFluxDirections,
         updateDirections   = {self.cdim+1},
         clearOut           = false,   -- Continue accumulating into output field.
      }
   end
   
   -- Create updaters to compute various moments.
   self.numDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "GkM0",
      gkfacs     = {self.mass, self.bmag},
   }
   self.momDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "GkM1",
      gkfacs     = {self.mass, self.bmag},
   }
   self.momProjDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "GkM1proj",
      gkfacs     = {self.mass, self.bmag},
   }
   self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "GkM2",
      gkfacs     = {self.mass, self.bmag},
   }
   self.M2parCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "GkM2par",
      gkfacs     = {self.mass, self.bmag},
   }
   if self.vdim > 1 then
      self.M2perpCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,
         phaseBasis = self.basis,
         confBasis  = self.confBasis,
         moment     = "GkM2perp",
         gkfacs     = {self.mass, self.bmag},
      }
   end
   self.threeMomentsCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "GkThreeMoments",
      gkfacs     = {self.mass, self.bmag},
   }
   if self.needSelfPrimMom then
      -- This is used in calcCouplingMoments to reduce overhead and multiplications.
      -- If collisions are LBO, the following also computes boundary corrections and, if polyOrder=1, star moments.
      self.threeMomentsLBOCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,
         phaseBasis = self.basis,
         confBasis  = self.confBasis,
         moment     = "GkThreeMomentsLBO",
         gkfacs     = {self.mass, self.bmag},
         positivity = self.positivity,
      }
      if self.needCorrectedSelfPrimMom then
         self.primMomSelf = Updater.SelfPrimMoments {
            onGrid     = self.confGrid,
            phaseBasis = self.basis,
            confBasis  = self.confBasis,
            operator   = "GkLBO",
            gkfacs     = {self.mass, self.bmag},
         }
      end
      -- Updaters for the primitive moments.
      self.confDiv = Updater.CartFieldBinOp {
         onGrid    = self.confGrid,
         weakBasis = self.confBasis,
         operation = "Divide",
      }
      self.confMul = Updater.CartFieldBinOp {
         onGrid    = self.confGrid,
         weakBasis = self.confBasis,
         operation = "Multiply",
      }
   end
   
   self._firstMomentCalc = true  -- To avoid re-calculating moments when not evolving.

   self.tmCouplingMom = 0.0      -- For timer.

   assert(self.n0, "Must specify background density as global variable 'n0' in species table as 'n0 = ...'")
end

function GkSpecies:initCrossSpeciesCoupling(species)
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

   -- Function to concatenate two tables.
   local function tableConcat(t1,t2)
      for i=1,#t2 do
         t1[#t1+1] = t2[i]
      end
      return t1
   end

   -- Create a double nested table of colliding species.
   -- In this table we will encode information about that collition such as:
   --   * does the collision take place?
   --   * Operator modeling the collision.
   --   * Is the collisionality constant in time?
   --   * Does it use spatially varying collisionality?
   --   * If using homogeneous collisionality, record its value/profile.
   -- Other features of a collision may be added in the future, such as
   -- velocity dependent collisionality, FLR effects, or some specific
   -- neutral/impurity effect.
   self.collPairs = {}
   for sN, _ in pairs(species) do
      self.collPairs[sN] = {}
      for sO, _ in pairs(species) do
         self.collPairs[sN][sO] = {}
         -- Need next below because species[].collisions is created as an empty table.
         if species[sN].collisions and next(species[sN].collisions) then
            -- This species collides with someone.
            local selfColl, crossColl, collSpecs = false, false, {}
            -- Obtain the boolean indicating if self/cross collisions affect the sN species.
            for collNm, _ in pairs(species[sN].collisions) do
               selfColl  = selfColl or species[sN].collisions[collNm].selfCollisions
               crossColl = crossColl or species[sN].collisions[collNm].crossCollisions
               collSpecs = tableConcat(collSpecs, species[sN].collisions[collNm].collidingSpecies)
            end

            -- Record if a specific binary collision is turned on.
            if sN == sO then
               self.collPairs[sN][sO].on = selfColl
            else
               if crossColl then
                  local specInd = findInd(collSpecs, sO)
                  if specInd < (#collSpecs+1) then
                     self.collPairs[sN][sO].on = true
                  else
                     self.collPairs[sN][sO].on = false
                  end
               else
                  self.collPairs[sN][sO].on = false
               end
            end

         else

            -- This species does not collide with anyone.
            self.collPairs[sN][sO].on = false

         end    -- end if next(species[sN].collisions) statement.
      end
   end

   -- Here we wish to record some properties of each collision in collPairs.
   for sN, _ in pairs(species) do
      -- Need next below because species[].collisions is created as an empty table.
      if species[sN].collisions and next(species[sN].collisions) then
         for sO, _ in pairs(species) do
            -- Find the kind of a specific collision, and the collision frequency it uses.
            for collNmN, _ in pairs(species[sN].collisions) do
               if self.collPairs[sN][sO].on then
                  local specInd = findInd(species[sN].collisions[collNmN].collidingSpecies, sO)
                  if specInd < (#species[sN].collisions[collNmN].collidingSpecies+1) then
                     -- Collision operator kind.
                     self.collPairs[sN][sO].kind      = species[sN].collisions[collNmN].collKind
                     -- Collision frequency time dependence (e.g. constant, time-varying).
                     self.collPairs[sN][sO].timeDepNu = species[sN].collisions[collNmN].timeDepNu
                     -- Collision frequency spatial dependence (e.g. homogeneous, spatially varying).
                     self.collPairs[sN][sO].varNu     = species[sN].collisions[collNmN].varNu
                     if (not self.collPairs[sN][sO].timeDepNu) then
                        -- Constant collisionality. Record it.
                        self.collPairs[sN][sO].nu = species[sN].collisions[collNmN].collFreqs[specInd]
                     else
                        -- Normalized collisionality to be scaled (e.g. by n_r/(v_{ts}^2+v_{tr}^2)^(3/2)).
                        if (species[sN].collisions[collNmN].userInputNormNu) then
                           self.collPairs[sN][sO].normNu = species[sN].collisions[collNmN].normNuIn[specInd]
                        else
                           self.collPairs[sN][sO].normNu = 0.0    -- Not used.
                        end
                     end
                  end
               elseif self.collPairs[sO][sN].on then
                  -- Species sN doesn't collide with sO, but sO collides with sN.
                  -- For computing cross-primitive moments, species sO may need the sN-sO
                  -- collision frequency. Set it such that m_sN*nu_{sN sO}=m_sO*nu_{sO sN}.
                  for collNmO, _ in pairs(species[sO].collisions) do
                     local specInd = findInd(species[sO].collisions[collNmO].collidingSpecies, sN)
                     if specInd < (#species[sO].collisions[collNmO].collidingSpecies+1) then
                        -- Collision operator kind.
                        self.collPairs[sO][sN].kind      = species[sO].collisions[collNmO].collKind
                        -- Collision frequency time dependence (e.g. constant, time-varying).
                        self.collPairs[sO][sN].timeDepNu = species[sO].collisions[collNmN].timeDepNu
                        self.collPairs[sN][sO].timeDepNu = species[sO].collisions[collNmN].timeDepNu
                        -- Collision frequency spatial dependence (e.g. homogeneous, spatially varying).
                        self.collPairs[sO][sN].varNu     = species[sO].collisions[collNmO].varNu
                        self.collPairs[sN][sO].varNu     = species[sO].collisions[collNmO].varNu
                        if (not self.collPairs[sN][sO].timeDepNu) then
                           -- Constant collisionality. Record it.
                           if (self.collPairs[sN][sO].varNu) then
                              -- We will need to first project the nu we do have, and later scale it by the mass ratio.
                              self.collPairs[sN][sO].nu = species[sO].collisions[collNmO].collFreqs[specInd]
                           else
                              self.collPairs[sN][sO].nu = (species[sO]:getMass()/species[sN]:getMass())*species[sO].collisions[collNmO].collFreqs[specInd]
                           end
                        else
                           -- Normalized collisionality to be scaled (e.g. by n_r/(v_{ts}^2+v_{tr}^2)^(3/2)).
                           if (species[sO].collisions[collNmO].userInputNormNu) then
                              self.collPairs[sN][sO].normNu = (species[sO]:getMass()/species[sN]:getMass())*species[sO].collisions[collNmO].normNuIn[specInd]
                           else
                              self.collPairs[sN][sO].normNu = 0.0    -- Not used.
                           end
                        end
                     end
                  end
               end
            end
         end    -- end if next(species[sN].collisions) statement.
      else
         -- This segment is needed when species sN doesn't have a collision object/table,
         -- but species sO collides with sN.
         -- For computing cross-primitive moments, species sO may need the sN-sO
         -- collision frequency. Set it such that m_sN*nu_{sN sO}=m_sO*nu_{sO sN}.
         for sO, _ in pairs(species) do
            if species[sO].collisions and next(species[sO].collisions) then
               for collNmO, _ in pairs(species[sO].collisions) do
                  if self.collPairs[sO][sN].on then
                     -- Species sO collides with sN. For computing cross-primitive moments,
                     -- species sO may need the sN-sO collision frequency. Set it such
                     -- that m_sN*nu_{sN sO}=m_sO*nu_{sO sN}.
                     local specInd = findInd(species[sO].collisions[collNmO].collidingSpecies, sN)
                     if specInd < (#species[sO].collisions[collNmO].collidingSpecies+1) then
                        self.collPairs[sO][sN].varNu = species[sO].collisions[collNmO].varNu
                        self.collPairs[sN][sO].varNu = species[sO].collisions[collNmO].varNu
                        if (not self.collPairs[sN][sO].timeDepNu) then
                           -- Constant collisionality. Record it.
                           if (self.collPairs[sN][sO].varNu) then
                              -- We will need to first project the nu we do have, and later scale it by the mass ratio.
                              self.collPairs[sN][sO].nu = species[sO].collisions[collNmO].collFreqs[specInd]
                           else
                              self.collPairs[sN][sO].nu = (species[sO]:getMass()/species[sN]:getMass())*species[sO].collisions[collNmO].collFreqs[specInd]
                           end
                        else
                           -- Normalized collisionality to be scaled (e.g. by n_r/(v_{ts}^2+v_{tr}^2)^(3/2)).
                           if (species[sO].collisions[collNmO].userInputNormNu) then
                              self.collPairs[sN][sO].normNu = (species[sO]:getMass()/species[sN]:getMass())*species[sO].collisions[collNmO].normNuIn[specInd]
                           else
                              self.collPairs[sN][sO].normNu = 0.0    -- Not used.
                           end
                        end
                     end
                  end
               end
            end
         end
      end
   end

   -- Determine if self primitive moments and boundary corrections are needed.
   -- If a pair of species only has cross-species collisions (no self-collisions)
   -- then the self-primitive moments may be computed without boundary corrections.
   -- Boundary corrections are only needed if there are LBO self-species collisions.
   self.needSelfPrimMom          = false
   self.needCorrectedSelfPrimMom = false
   -- Also check if spatially varying nu is needed, and if the user inputed a spatial
   -- profile for the collisionality (which needs to be projected).
   local needVarNu               = false
   local userInputNuProfile      = false
   if self.collPairs[self.name][self.name].on then
      self.needSelfPrimMom          = true
      if (self.collPairs[self.name][self.name].kind=="GkLBO") or
         (self.collPairs[self.name][self.name].kind=="VmLBO") then
         self.needCorrectedSelfPrimMom = true
      end
   end
   for sO, _ in pairs(species) do
      if self.collPairs[self.name][sO].on or self.collPairs[sO][self.name].on then
         self.needSelfPrimMom = true
         if ( self.collPairs[sO][sO].on and (self.collPairs[self.name][self.name].kind=="GkLBO" or
                                             self.collPairs[self.name][self.name].kind=="VmLBO") ) then
            self.needCorrectedSelfPrimMom = true
         end

         if self.collPairs[self.name][sO].varNu or self.collPairs[sO][self.name].varNu then
            needVarNu = true
            if (not self.collPairs[self.name][sO].timeDepNu) or (not self.collPairs[sO][self.name].timeDepNu) then
               userInputNuProfile = true
            end
         end
      end
   end

   if self.needSelfPrimMom then
      -- Allocate fields to store self-species primitive moments.
      self.uParSelf = self:allocMoment()
      self.vtSqSelf = self:allocMoment()

      -- Allocate fields for boundary corrections.
      self.m1Correction = self:allocMoment()
      self.m2Correction = self:allocMoment()

      -- Allocate fields for star moments (only used with polyOrder=1).
      if (self.basis:polyOrder()==1) then
         self.m0Star = self:allocMoment()
         self.m1Star = self:allocMoment()
         self.m2Star = self:allocMoment()
      end
   end

   -- Allocate fieds to store cross-species primitive moments.
   self.uParCross = {}
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
         if (sN ~= sO) and (self.collPairs[sN][sO].on or self.collPairs[sO][sN].on) then
            otherNm = string.gsub(sO .. sN, self.name, "")
            if (self.uParCross[otherNm] == nil) then
               self.uParCross[otherNm] = self:allocMoment()
            end
            if (self.vtSqCross[otherNm] == nil) then
               self.vtSqCross[otherNm] = self:allocMoment()
            end
         end
      end
   end

   if needVarNu then
      self.nuVarXCross = {}    -- Collisionality varying in configuration space.
      local projectNuX = nil
      if userInputNuProfile then
         projectNuX = Updater.ProjectOnBasis {
            onGrid          = self.confGrid,
            basis           = self.confBasis,
            evaluate        = function(t,xn) return 0.0 end, -- Function is set below.
            projectOnGhosts = false,
         }
      end
      for sN, _ in pairs(species) do
         if sN ~= self.name then
            -- Sixth moment flag is to indicate if spatially varying collisionality has been computed.
            self.momentFlags[6][sN] = false
         end

         for sO, _ in pairs(species) do
            -- Allocate space for this species' collision frequency
            -- only if some other species collides with it.
            if (sN ~= sO) and (self.collPairs[sN][sO].on or self.collPairs[sO][sN].on) then
               otherNm = string.gsub(sO .. sN, self.name, "")
               if self.nuVarXCross[otherNm] == nil then
                  self.nuVarXCross[otherNm] = self:allocMoment()
                  if (userInputNuProfile and (not self.collPairs[sN][sO].timeDepNu) or (not self.collPairs[sO][sN].timeDepNu)) then
                     projectNuX:setFunc(self.collPairs[self.name][otherNm].nu)
                     projectNuX:advance(0.0,{},{self.nuVarXCross[otherNm]})
                     if (not self.collPairs[self.name][otherNm].on) then
                        self.nuVarXCross[otherNm]:scale(species[self.name]:getMass()/species[otherNm]:getMass())
                     end
                     self.nuVarXCross[otherNm]:write(string.format("%s_nu-%s_%d.bp",self.name,otherNm,0),0.0,0,true)
                  end
               end
            end
         end
      end
   end

end

function GkSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self.tCurr = tCurr
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em          = emIn[1]:rkStepperFields()[inIdx]
   local dApardtProv = emIn[1].dApardtProv
   local emFunc      = emIn[2]:rkStepperFields()[1]

   -- Solvers specified with clearOut = false, so we need to zero out RHS here.
   fRhsOut:clear(0.0)

   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
         if (c.collKind == "GkLBO") then
            c.gkLBOconstNuCalcEq:clearRhsTerms()
         end
         c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
         c:advance(tCurr, fIn, species, fRhsOut)

         if (self.positivity) and (c.collKind == "GkLBO") then
            c.gkLBOconstNuCalcEq:getPositivityRhs(tCurr, self.dtGlobal[0], fIn, fRhsOut)

            -- Set the time step and CFL again.
            c.collisionSlvrStep2:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
            c:advanceStep2(tCurr, fIn, species, fRhsOut)
            c.gkLBOconstNuCalcEq:getPositivityRhsStep2(tCurr, self.dtGlobal[0], fIn, fRhsOut)
         end
      end
   end
   if self.evolveCollisionless then
      self.gkEqn:clearRhsTerms()
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solver:advance(tCurr, {fIn, em, emFunc, dApardtProv}, {fRhsOut})

      if self.positivity then
         self.gkEqn:getPositivityRhs(tCurr, self.dtGlobal[0], fIn, fRhsOut)
      end
   else
      self.gkEqn:setAuxFields({em, emFunc, dApardtProv})  -- Set auxFields in case they are needed by BCs/collisions.
   end

   -- combine positivity weights from collisionless and collisions (for use on next step)
   if self.positivity then
      self.gkEqn:setPositivityWeights(self.cflRateByCell)
      for _, c in pairs(self.collisions) do
         if (c.collKind == "GkLBO") then
            c.gkLBOconstNuCalcEq:setPositivityWeights(self.cflRateByCell)
         end
      end
   end

   if self.fSource and self.evolveSources then
      -- Add source it to the RHS.
      -- Barrier over shared communicator before accumulate.
      Mpi.Barrier(self.grid:commSet().sharedComm)
      fRhsOut:accumulate(self.sourceTimeDependence(tCurr), self.fSource)
   end
end

function GkSpecies:advanceStep2(tCurr, species, emIn, inIdx, outIdx)
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em          = emIn[1]:rkStepperFields()[inIdx]
   local dApardtProv = emIn[1].dApardtProv
   local emFunc      = emIn[2]:rkStepperFields()[1]

   if self.evolveCollisionless then
      self.solverStep2:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solverStep2:advance(tCurr, {fIn, em, emFunc, dApardtProv}, {fRhsOut})
   end
end

function GkSpecies:advanceStep3(tCurr, species, emIn, inIdx, outIdx)
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em          = emIn[1]:rkStepperFields()[inIdx]
   local dApardtProv = emIn[1].dApardtProv
   local emFunc      = emIn[2]:rkStepperFields()[1]

   if self.evolveCollisionless then
      self.gkEqn.ohmMod:clear(0.)
      self.gkEqn.fRhsSurfV:clear(0.)
 
      self.solverStep2:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solverStep2:advance(tCurr, {fIn, em, emFunc, dApardtProv}, {fRhsOut})
   end
end

function GkSpecies:advanceStep4(tCurr, species, emIn, inIdx, outIdx)
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em          = emIn[1]:rkStepperFields()[inIdx]
   local dApardtProv = emIn[1].dApardtProv
   local emFunc      = emIn[2]:rkStepperFields()[1]

   if self.evolveCollisionless then
      self.solverStep3:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solverStep3:advance(tCurr, {fIn, em, emFunc, dApardtProv}, {fRhsOut})

      if self.positivity then
         self.gkEqn:getPositivityRhsStep2(tCurr, self.dtGlobal[0], fIn, fRhsOut)
      end
   end
end

function GkSpecies:createDiagnostics()
   local function isIntegratedMomentNameGood(nm)
      if nm == "intM0" or nm == "intM1" or nm == "intM2" or nm == "intL1" or nm == "intL2"
      or nm == "intDelM0" or nm == "intDelM2" or nm == "intDelL2"
      or nm == "intDelPosM0" or nm == "intDelPosM2" or nm == "intDelPosL2" then
         return true
      end
      return false
   end
   self.diagnosticIntegratedMomentFields   = { }
   self.diagnosticIntegratedMomentUpdaters = { } 
   -- Allocate space to store integrated moments and create integrated moment updaters.
   for i, mom in pairs(self.diagnosticIntegratedMoments) do
      if isIntegratedMomentNameGood(mom) then
         self.diagnosticIntegratedMomentFields[mom] = DataStruct.DynVector {
            numComponents = 1,
         }
         if mom == "intL2" or "intDelL2" or "intDelPosL2" then
            self.diagnosticIntegratedMomentUpdaters[mom] = Updater.CartFieldIntegratedQuantCalc {
               onGrid        = self.grid,
               basis         = self.basis,
               numComponents = 1,
               quantity      = "RmsV"
            }
         elseif mom == "intL1" then
            self.diagnosticIntegratedMomentUpdaters[mom] = Updater.CartFieldIntegratedQuantCalc {
               onGrid        = self.grid,
               basis         = self.basis,
               numComponents = 1,
               quantity      = "AbsV"
            }
         else
            self.diagnosticIntegratedMomentUpdaters[mom] = Updater.CartFieldIntegratedQuantCalc {
               onGrid        = self.confGrid,
               basis         = self.confBasis,
               numComponents = 1,
               quantity      = "V"
            }
         end
      else
         assert(false, string.format("Integrated moment %s not valid", mom))
      end
   end
   
   -- Function to check if moment name is correct.
   local function isMomentNameGood(nm)
      return Updater.DistFuncMomentCalc:isGkMomentNameGood(nm)
   end
   local function isWeakMomentNameGood(nm)
      return nm == "GkUpar" or nm == "GkVtSq" or nm == "GkTpar" or nm == "GkTperp" or nm == "GkTemp"
   end
   local function isAuxMomentNameGood(nm)
      return nm == "GkBeta" or nm == "GkUparCross" or nm == "GkVtSqCross"
   end
   local function contains(table, element)
     for _, value in pairs(table) do
       if value == element then
         return true
       end
     end
     return false
   end

   self.diagnosticMomentFields   = { }
   self.diagnosticMomentUpdaters = { } 
   self.diagnosticWeakMoments    = { }
   self.diagnosticAuxMoments     = { }
   self.weakMomentOpFields       = { }
   self.weakMomentScaleFac       = { }
   -- Set up weak multiplication and division operators.
   self.weakMultiplication = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Multiply",
      onGhosts  = true,
   }
   self.weakDivision = Updater.CartFieldBinOp {
      onGrid     = self.confGrid,
      weakBasis  = self.confBasis,
      operation  = "Divide",
      onGhosts   = true,
      positivity = self.positivity,
   }
   self.weakDivisionPhase = Updater.CartFieldBinOp {
      onGrid     = self.grid,
      weakBasis  = self.basis,
      fieldBasis = self.confBasis,
      operation  = "Divide",
      onGhosts   = true,
   }

   -- Sort moments into diagnosticMoments, diagnosticWeakMoments, and diagnosticAuxMoments.
   for i, mom in pairs(self.diagnosticMoments) do
      if isWeakMomentNameGood(mom) then
         -- Remove moment name from self.diagnosticMoments list, and add it to self.diagnosticWeakMoments list.
         self.diagnosticWeakMoments[mom] = true
         self.diagnosticMoments[i]       = nil
      elseif isAuxMomentNameGood(mom) then
         -- Remove moment name from self.diagnosticMoments list, and add it to self.diagnosticAuxMoments list.
         if mom == "GkUparCross" then
            for nm, _ in pairs(self.uParCross) do
               -- Create one diagnostic for each cross velocity (used in collisions).
               self.diagnosticAuxMoments[mom .. "-" .. nm] = true
            end
         elseif mom == "GkVtSqCross" then
            for nm, _ in pairs(self.vtSqCross) do
               -- Create one diagnostic for each cross temperature (used in collisions).
               self.diagnosticAuxMoments[mom .. "-" .. nm] = true
            end
         else
            self.diagnosticAuxMoments[mom] = true
         end
         self.diagnosticMoments[i] = nil
      end
   end

   -- Make sure we have the updaters needed to calculate all the aux moments.
   for mom, _ in pairs(self.diagnosticAuxMoments) do
      if mom == "GkBeta" then
         if not self.diagnosticWeakMoments["GkTemp"] then 
            self.diagnosticWeakMoments["GkTemp"] = true
         end
         if not contains(self.diagnosticMoments, "GkM0") then
            table.insert(self.diagnosticMoments, "GkM0")
         end
      end
   end

   -- Make sure we have the updaters needed to calculate all the weak moments.
   for mom, _ in pairs(self.diagnosticWeakMoments) do
      -- All GK weak moments require M0 = density.
      if not contains(self.diagnosticMoments, "GkM0") then
         table.insert(self.diagnosticMoments, "GkM0")
      end

      if mom == "GkUpar" then
         if not contains(self.diagnosticMoments, "GkM1") then
            table.insert(self.diagnosticMoments, "GkM1")
         end
      end
      if mom == "GkVtSq" then
         if not contains(self.diagnosticMoments, "GkM2") then
            table.insert(self.diagnosticMoments, "GkM2")
         end
         if not self.diagnosticWeakMoments["GkUpar"] then
            self.diagnosticWeakMoments["GkUpar"] = true
         end
      elseif mom == "GkTpar" then
         if not contains(self.diagnosticMoments, "GkM2par") then
            table.insert(self.diagnosticMoments, "GkM2par")
         end
         if not self.diagnosticWeakMoments["GkUpar"] then
            self.diagnosticWeakMoments["GkUpar"] = true
         end
         if not contains(self.diagnosticMoments, "GkM1") then
            -- First moment is needed by GkUpar.
            table.insert(self.diagnosticMoments, "GkM1")
         end
      elseif mom == "GkTperp" then
         if not contains(self.diagnosticMoments, "GkM2perp") then
            table.insert(self.diagnosticMoments, "GkM2perp")
         end
      elseif mom == "GkTemp" then 
         if not contains(self.diagnosticMoments, "GkM2") then
            table.insert(self.diagnosticMoments, "GkM2")
         end      
         if not self.diagnosticWeakMoments["GkUpar"] then
            self.diagnosticWeakMoments["GkUpar"] = true
         end
         if not contains(self.diagnosticMoments, "GkM1") then
            -- First moment is needed by GkUpar.
            table.insert(self.diagnosticMoments, "GkM1")
         end
      end
   end

   -- Allocate space to store moments and create moment updaters.
   for i, mom in pairs(self.diagnosticMoments) do
      if isMomentNameGood(mom) then
         self.diagnosticMomentFields[mom] = DataStruct.Field {
            onGrid        = self.confGrid,
            numComponents = self.confBasis:numBasis(),
            ghost         = {1, 1},
	    metaData = {
	       polyOrder = self.basis:polyOrder(),
	       basisType = self.basis:id()
	    },	    
         }
         self.diagnosticMomentUpdaters[mom] = Updater.DistFuncMomentCalc {
            onGrid     = self.grid,
            phaseBasis = self.basis,
            confBasis  = self.confBasis,
            moment     = mom,
            gkfacs     = {self.mass, self.bmag},
         }
      else
         assert(false, string.format("Moment %s not valid", mom))
      end
   end
   for mom, _ in pairs(self.diagnosticWeakMoments) do
      if isWeakMomentNameGood(mom) then
         self.diagnosticMomentFields[mom] = DataStruct.Field {
            onGrid        = self.confGrid,
            numComponents = self.confBasis:numBasis(),
            ghost         = {1, 1},
	    metaData = {
	       polyOrder = self.basis:polyOrder(),
	       basisType = self.basis:id()
	    },	    
         }
      else
         assert(false, string.format("Moment %s not valid", mom))
      end

      if mom == "GkUpar" then
         self.weakMomentOpFields["GkUpar"] = {self.diagnosticMomentFields["GkM0"], self.diagnosticMomentFields["GkM1"]}
      elseif mom == "GkVtSq" then 
         self.weakMomentOpFields["GkVtSq"] = {self.diagnosticMomentFields["GkM0"], self.diagnosticMomentFields["GkM2"]}
         self.weakMomentScaleFac["GkVtSq"] = 1.0/self.vDegFreedom
      elseif mom == "GkTpar" then
         self.weakMomentOpFields["GkTpar"] = {self.diagnosticMomentFields["GkM0"], self.diagnosticMomentFields["GkM2par"]}
         self.weakMomentScaleFac["GkTpar"] = self.mass
      elseif mom == "GkTperp" then
         self.weakMomentOpFields["GkTperp"] = {self.diagnosticMomentFields["GkM0"], self.diagnosticMomentFields["GkM2perp"]}
         self.weakMomentScaleFac["GkTperp"] = self.mass
      elseif mom == "GkTemp" then 
         self.weakMomentOpFields["GkTemp"] = {self.diagnosticMomentFields["GkM0"], self.diagnosticMomentFields["GkM2"]}
         self.weakMomentScaleFac["GkTemp"] = self.mass/self.vDegFreedom
      end
   end

   for mom, _ in pairs(self.diagnosticAuxMoments) do
      self.diagnosticMomentFields[mom] = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
	 metaData = {
	       polyOrder = self.basis:polyOrder(),
	       basisType = self.basis:id()
	 },	    
      }
   end
end

function GkSpecies:calcDiagnosticIntegratedMoments(tCurr)
   -- First compute M0, M1, M2.
   local fIn = self:rkStepperFields()[1]
   self.threeMomentsCalc:advance(tCurr, {fIn}, { self.numDensity, self.momDensity, self.ptclEnergy })

   local fDel = self:rkStepperFields()[2]
   fDel:combine(1, fIn, -1, self.fPrev)
   if self.positivityDiffuse then
      fDel:accumulate(-1, self.fDelPos[1])
   end
   self.fPrev:copy(fIn)
   self.threeMomentsCalc:advance(tCurr, {fDel}, { self.numDensityAux, self.momDensityAux, self.ptclEnergyAux })

   if self.positivityDiffuse then
      self.threeMomentsCalc:advance(tCurr, {self.fDelPos[1]}, { self.numDensityPos, self.momDensityPos, self.ptclEnergyPos })
   end

   for i, mom in pairs(self.diagnosticIntegratedMoments) do
      if mom == "intM0" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.numDensity}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intM1" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.momDensity}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intM2" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.ptclEnergy}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intL1" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.distf[1]}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intL2" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.distf[1]}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intDelM0" then 
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.numDensityAux}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intDelM2" then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.ptclEnergyAux}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intDelL2" then 
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {fDel}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intDelPosL2" and self.positivityDiffuse then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.fDelPos[1]}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intDelPosM0" and self.positivityDiffuse then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.numDensityPos}, {self.diagnosticIntegratedMomentFields[mom]})
      elseif mom == "intDelPosM2" and self.positivityDiffuse then
         self.diagnosticIntegratedMomentUpdaters[mom]:advance(
            tCurr, {self.ptclEnergyPos}, {self.diagnosticIntegratedMomentFields[mom]})
      end
   end

end

function GkSpecies:calcDiagnosticWeakMoments()
   GkSpecies.super.calcDiagnosticWeakMoments(self)
   if self.diagnosticWeakMoments["GkVtSq"] then
      -- Need to subtract (uPar^2)/vdim from vtSq (which at this point holds M2/(vdim*M0)).
      -- uPar is calculated in KineticEnergySpecies:calcDiagnositcWeakMoments().
      self.weakMultiplication:advance(0.0,
           {self.diagnosticMomentFields["GkUpar"], self.diagnosticMomentFields["GkUpar"]}, 
           {self.momDensityAux})
      -- Barrier over shared communicator before accumulate.
      Mpi.Barrier(self.grid:commSet().sharedComm)
      self.diagnosticMomentFields["GkVtSq"]:accumulate(-self.weakMomentScaleFac["GkVtSq"], self.momDensityAux)
   end
   -- Need to subtract m*Upar^2 from GkTemp and GkTpar.
   if self.diagnosticWeakMoments["GkTemp"] or self.diagnosticWeakMoments["GkTpar"] then
      self.weakMultiplication:advance(0.0,
           {self.diagnosticMomentFields["GkUpar"], self.diagnosticMomentFields["GkUpar"]}, 
           {self.momDensityAux})
   end
   -- Barrier over shared communicator before accumulate.
   Mpi.Barrier(self.grid:commSet().sharedComm)
   if self.diagnosticWeakMoments["GkTemp"] then
      self.diagnosticMomentFields["GkTemp"]:accumulate(-self.weakMomentScaleFac["GkTemp"], self.momDensityAux)
   end
   if self.diagnosticWeakMoments["GkTpar"] then
      self.diagnosticMomentFields["GkTpar"]:accumulate(-self.mass, self.momDensityAux)
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
   for nm, _ in pairs(self.diagnosticAuxMoments) do
      if string.find(nm, "GkUparCross") then
         otherNm = string.gsub(nm, "GkUparCross%-", "")
         self.diagnosticMomentFields[nm]:copy(self.uParCross[otherNm])
      end
      if string.find(nm, "GkVtSqCross") then
         otherNm = string.gsub(nm, "GkVtSqCross%-", "")
         self.diagnosticMomentFields[nm]:copy(self.vtSqCross[otherNm])
      end
   end
end

-- BC functions.
function GkSpecies:bcReflectFunc(dir, tm, idxIn, fIn, fOut)
   -- skinLoop should be "flip"
   -- Note that GK reflection only valid in z-vpar.
   -- This is checked when bc is created.

   self.basis:flipSign(dir, fIn, fOut)
   -- vpar is always first velocity dimension
   local vpardir=self.cdim+1 
   self.basis:flipSign(vpardir, fOut, fOut)
end
function GkSpecies:bcSheathFunc(dir, tm, idxIn, fIn, fOut)
   -- skinLoop should be "flip"
   -- Note that GK reflection only valid in z-vpar.
   -- This is checked when bc is created.

   -- Need to figure out if we are on lower or upper domain edge
   local edgeVal
   local globalRange = self.grid:globalRange()
   if idxIn[dir] == globalRange:lower(dir) then 
      -- This means we are at lower domain edge, 
      -- so we need to evaluate basis functions at z=-1.
      edgeVal = -1 
   else 
      -- This means we are at upper domain edge
      -- so we need to evaluate basis functions at z=1.
      edgeVal = 1 
   end
   -- Get vpar limits of cell.
   local vpardir = self.cdim+1
   local gridIn  = self.grid
   gridIn:setIndex(idxIn)
   local vL = gridIn:cellLowerInDir(vpardir)
   local vR = gridIn:cellUpperInDir(vpardir)
   local vlowerSq, vupperSq
   -- This makes it so that we only need to deal with absolute values of vpar.
   if math.abs(vR)>=math.abs(vL) then
      vlowerSq = vL*vL
      vupperSq = vR*vR
   else
      vlowerSq = vR*vR
      vupperSq = vL*vL
   end
   local w = gridIn:cellCenterInDir(vpardir)
   local dv = gridIn:dx(vpardir)
   local fhat = self.fhatSheathPtr -- distribution function to be reflected
   self.fhatSheath:fill(self.fhatSheathIdxr(idxIn), fhat)
   -- calculate reflected distribution function fhat
   -- note: reflected distribution can be 
   -- 1) fhat=0 (no reflection, i.e. absorb), 
   -- 2) fhat=f (full reflection)
   -- 3) fhat=c*f (partial reflection)
   self.gkEqn:calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, self.charge, self.mass, idxIn, fIn, fhat)
   -- reflect fhat into skin cells
   self:bcReflectFunc(dir, tm, nil, fhat, fOut) 
end

function GkSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- Need to wrap member functions so that self is passed.
   local function bcAbsorbFunc(...)  return self:bcAbsorbFunc(...) end
   local function bcOpenFunc(...)    return  self:bcOpenFunc(...) end
   local function bcReflectFunc(...) return self:bcReflectFunc(...) end
   local function bcSheathFunc(...)  return self:bcSheathFunc(...) end
   local function bcCopyFunc(...)    return self:bcCopyFunc(...) end
   
   local vdir = nil
   if dir==self.cdim then 
      vdir = self.cdim+1 
   end

   if bcType == SP_BC_ABSORB then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcAbsorbFunc }, "pointwise"))
   elseif bcType == SP_BC_OPEN then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcOpenFunc }, "pointwise"))
   -- Note: reflection and sheath BCs only make sense in z direction,
   -- which is always last config space direction, i.e. dir = self.cdim.
   elseif bcType == SP_BC_REFLECT and dir==self.cdim then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, vdir, edge, { bcReflectFunc }, "flip"))
   elseif bcType == SP_BC_SHEATH and dir==self.cdim then
      self.fhatSheath     = self:allocDistf()
      self.fhatSheathPtr  = self.fhatSheath:get(1)
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

   -- Compute moments needed in coupling to fields and collisions.
   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()

      if self.deltaF then
        fIn:accumulate(-1.0, self.f0)
      end
      
      if self.needSelfPrimMom then
         self.threeMomentsLBOCalc:advance(tCurr, {fIn}, { self.numDensity, self.momDensity, self.ptclEnergy,
                                                          self.m1Correction, self.m2Correction,
                                                          self.m0Star, self.m1Star, self.m2Star })
         if self.needCorrectedSelfPrimMom then
            -- Also compute self-primitive moments uPar and vtSq.
            self.primMomSelf:advance(tCurr, {self.numDensity, self.momDensity, self.ptclEnergy,
                                             self.m1Correction, self.m2Correction,
                                             self.m0Star, self.m1Star, self.m2Star}, {self.uParSelf, self.vtSqSelf})
         else
            -- Compute self-primitive moments with binOp updaters.
            self.confDiv:advance(tCurr, {self.numDensity, self.momDensity}, {self.uParSelf})
            self.confMul:advance(tCurr, {self.uParSelf, self.momDensity}, {self.numDensityAux})
            -- Barrier over shared communicator before combine
            Mpi.Barrier(self.grid:commSet().sharedComm)
            self.momDensityAux:combine( 1.0/self.vDegFreedom, self.ptclEnergy,
                                       -1.0/self.vDegFreedom, self.numDensityAux )
            self.confDiv:advance(tCurr, {self.numDensity, self.momDensityAux}, {self.vtSqSelf})
         end
         -- Indicate that moments, boundary corrections, star moments
         -- and self-primitive moments have been computed.
         for iF=1,4 do
            self.momentFlags[iF] = true
         end
      else
         self.numDensityCalc:advance(tCurr, {fIn}, { self.numDensity })
         -- Indicate that first moment has been computed.
         self.momentFlags[1] = true
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

function GkSpecies:boundaryCorrections()
   return { self.m1Correction, self.m2Correction }
end

function GkSpecies:starMoments()
   return { self.m0Star, self.m1Star, self.m2Star }
end

function GkSpecies:selfPrimitiveMoments()
   return { self.uParSelf, self.vtSqSelf }
end

function GkSpecies:crossPrimitiveMoments(otherSpeciesName)
   return { self.uParCross[otherSpeciesName], self.vtSqCross[otherSpeciesName] }
end

function GkSpecies:getNumDensity(rkIdx)
   -- If no rkIdx specified, assume numDensity has already been calculated.
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

function GkSpecies:getMomDensity(rkIdx, addVolumeTerm)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then return self.momDensity end 
   local fIn = self:rkStepperFields()[rkIdx]
 
   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      if addVolumeTerm and self.gkEqn._isElectromagnetic then
        fIn:accumulate(1.0, self.gkEqn.fRhsVolV)
      elseif addVolumeTerm then
        fIn:accumulate(1.0, self.gkEqn.fRhsVol)
      end
      if self.deltaF then
        fIn:accumulate(-1.0, self.f0)
      end
      self.momDensityCalc:advance(nil, {fIn}, { self.momDensityAux })
      if addVolumeTerm and self.gkEqn._isElectromagnetic then
        fIn:accumulate(-1.0, self.gkEqn.fRhsVolV)
      elseif addVolumeTerm then
        fIn:accumulate(-1.0, self.gkEqn.fRhsVol)
      end
      if self.deltaF then
        fIn:accumulate(1.0, self.f0)
      end
      self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
   return self.momDensityAux
end

-- Like getMomDensity, but use GkM1proj instead of GkM1, which uses cell-average v_parallel in moment calculation.
function GkSpecies:getMomProjDensity(rkIdx, addVolumeTerm)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then return self.momDensity end 
   local fIn = self:rkStepperFields()[rkIdx]
 
   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      if addVolumeTerm and self.gkEqn._isElectromagnetic then
        fIn:accumulate(1.0, self.gkEqn.fRhsSurfV)
      elseif addVolumeTerm then
        fIn:accumulate(1.0, self.gkEqn.fRhsVol)
      end
      if self.deltaF then
        fIn:accumulate(-1.0, self.f0)
      end
      self.momProjDensityCalc:advance(nil, {fIn}, { self.momDensityAux })
      if addVolumeTerm and self.gkEqn._isElectromagnetic then
        fIn:accumulate(-1.0, self.gkEqn.fRhsSurfV)
      elseif addVolumeTerm then
        fIn:accumulate(-1.0, self.gkEqn.fRhsVol)
      end
      if self.deltaF then
        fIn:accumulate(1.0, self.f0)
      end
      self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
   return self.momDensityAux
end

function GkSpecies:getOhmModifier(rkIdx)
   -- For p > 1, this is just numDensity.
   if self.basis:polyOrder() > 1 then return self:getNumDensity(rkIdx) end

   local fIn = self.gkEqn.ohmMod
   --self.weakDivisionPhase:advance(0.0, {self.gkEqn.dApardtProv, fIn}, {fIn})

   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()
      if self.deltaF then
        fIn:accumulate(-1.0, self.f0)
      end
      self.momProjDensityCalc:advance(nil, {fIn}, { self.momDensityAux })
      if self.deltaF then
        fIn:accumulate(1.0, self.f0)
      end
      self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
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
   if self.solverStep3 then timer = timer + self.solverStep3.totalTime end
   if self.posRescaler then timer = timer + self.posRescaler.totalTime end
   return timer
end

function GkSpecies:Maxwellian(xn, n0, T0, vdIn)
   local vd   = vdIn or 0.0
   local vt2  = T0/self.mass
   local vpar = xn[self.cdim+1]
   local v2   = (vpar-vd)^2
   if self.vdim > 1 then 
     local mu = xn[self.cdim+2]
     v2 = v2 + 2*math.abs(mu)*self.bmagFunc(0,xn)/self.mass
     return n0*(2*math.pi*vt2)^(-3/2)*math.exp(-v2/(2*vt2))
   else
     return n0*(2*math.pi*vt2)^(-1/2)*math.exp(-v2/(2*vt2))
   end
end

return GkSpecies
