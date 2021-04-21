-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov species
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct     = require "DataStruct"
local Grid           = require "Grid"
local KineticSpecies = require "App.Species.KineticSpecies"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local VlasovEq       = require "Eq.Vlasov"
local ffi            = require "ffi"
local xsys           = require "xsys"
local lume           = require "Lib.lume"

local VlasovSpecies = Proto(KineticSpecies)

-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB  = 1
local SP_BC_REFLECT = 3
local SP_BC_EXTERN  = 4
local SP_BC_COPY    = 5
-- AHH: This was 2 but seems that is unstable. So using plain copy.
local SP_BC_OPEN      = SP_BC_COPY
local SP_BC_ZEROFLUX  = 6

VlasovSpecies.bcAbsorb    = SP_BC_ABSORB     -- Absorb all particles.
VlasovSpecies.bcOpen      = SP_BC_OPEN       -- Zero gradient.
VlasovSpecies.bcCopy      = SP_BC_COPY       -- Copy stuff.
VlasovSpecies.bcReflect   = SP_BC_REFLECT    -- Specular reflection.
VlasovSpecies.bcExternal  = SP_BC_EXTERN     -- Load external BC file.
VlasovSpecies.bcZeroFlux  = SP_BC_ZEROFLUX

function VlasovSpecies:alloc(nRkDup)
   -- Allocate distribution function.
   VlasovSpecies.super.alloc(self, nRkDup)

   -- Allocate fields to store coupling moments (for use in coupling
   -- to field and collisions).
   self.numDensity = self:allocMoment()
   self.momDensity = self:allocVectorMoment(self.vdim)
   self.ptclEnergy = self:allocMoment()

   -- Allocate field to accumulate externalField if any.
   --self.totalEmField = self:allocVectorMoment(8)     -- 8 components of EM field.

   -- Allocate field for external forces if any.
   if self.hasExtForce then 
      self.vExtForce = self:allocVectorMoment(self.vdim)
   end

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
      self.hasExtForce = true
   else
      self.hasExtForce = false
   end

   -- numVelFlux used for selecting which type of numerical flux function to use in velocity space
   -- defaults to "penalty" in Eq object, supported options: "penalty," "recovery"
   -- only used for DG Maxwell.
   self.numVelFlux = tbl.numVelFlux
end

function VlasovSpecies:allocMomCouplingFields()
   return { currentDensity = self:allocVectorMoment(self.vdim) }
end


function VlasovSpecies:createSolver(hasE, hasB, funcField, plasmaB)
   -- Run the KineticSpecies 'createSolver()' to initialize the
   -- collisions solver.
   VlasovSpecies.super.createSolver(self)

   -- External forces are accumulated to the electric field part of
   -- totalEmField
   if self.hasExtForce then
      hasE = true
      hasB = true
   end

   -- Allocate field to accumulate funcField if any.
   if hasB then
      self.totalEmField = self:allocVectorMoment(8)     -- 8 components of EM field.
   else
      self.totalEmField = self:allocVectorMoment(3)     -- Electric field only.
   end

   self.computePlasmaB = true and plasmaB   -- Differentiate plasma B from external B.

   -- Create updater to advance solution by one time-step.
   self.equation = VlasovEq {
      onGrid           = self.grid,
      phaseBasis       = self.basis,
      confBasis        = self.confBasis,
      charge           = self.charge,
      mass             = self.mass,
      hasElectricField = hasE,
      hasMagneticField = hasB,
      plasmaMagField   = plasmaB,
      numVelFlux       = self.numVelFlux,
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
      equation           = self.equation,
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
   self.calcMaxwell = Updater.MaxwellianOnBasis {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confGrid   = self.confGrid,
      confBasis  = self.confBasis,
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
         onGrid   = self.confGrid,
         basis    = self.confBasis,
         evaluate = self.vlasovExtForceFunc,
         onGhosts = false
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
   self.collPairs  = {}
   for sN, _ in lume.orderedIter(species) do
      self.collPairs[sN] = {}
      for sO, _ in lume.orderedIter(species) do
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
   for sN, _ in lume.orderedIter(species) do
      -- Need next below because species[].collisions is created as an empty table. 
      if species[sN].collisions and next(species[sN].collisions) then 
         for sO, _ in lume.orderedIter(species) do
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
         for sO, _ in lume.orderedIter(species) do
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
   for sO, _ in lume.orderedIter(species) do
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

   -- If ionization collision object exists, locate electrons
   local counterIz_elc = true
   local counterIz_neut = true
   for sN, _ in lume.orderedIter(species) do
      if species[sN].collisions and next(species[sN].collisions) then 
         for sO, _ in lume.orderedIter(species) do
	    if self.collPairs[sN][sO].on then
   	       if (self.collPairs[sN][sO].kind == 'Ionization') then
   		  for collNm, _ in pairs(species[sN].collisions) do
   		     if self.name==species[sN].collisions[collNm].elcNm and counterIz_elc then
   			self.neutNmIz         = species[sN].collisions[collNm].neutNm
   			self.needSelfPrimMom  = true
   			self.calcReactRate    = true
   			self.collNmIoniz      = collNm
			self.voronovReactRate = self:allocMoment()
   			self.vtSqIz           = self:allocMoment()
   			self.m0fMax           = self:allocMoment()
   			self.m0mod            = self:allocMoment()
   			self.fMaxwellIz       = self:allocDistf()
			self.intSrcIzM0       = DataStruct.DynVector{numComponents = 1}
			counterIz_elc         = false
		     elseif self.name==species[sN].collisions[collNm].neutNm and counterIz_neut then
			self.needSelfPrimMom = true
			counterIz_neut       = false
   		     end
   		  end
   	       end
   	    end
   	 end
      end
   end

   -- If Charge Exchange collision object exists, locate ions
   local counterCX = true
   for sN, _ in lume.orderedIter(species) do
      if species[sN].collisions and next(species[sN].collisions) then 
         for sO, _ in lume.orderedIter(species) do
   	    if self.collPairs[sN][sO].on then
   	       if (self.collPairs[sN][sO].kind == 'CX') then
   		  for collNm, _ in pairs(species[sN].collisions) do
   		     if self.name==species[sN].collisions[collNm].ionNm and counterCX then
  			self.calcCXSrc        = true			
			self.collNmCX         = collNm
   			self.neutNmCX         = species[sN].collisions[collNm].neutNm
   			self.needSelfPrimMom  = true
   			self.vSigmaCX         = self:allocMoment()
			species[self.neutNmCX].needSelfPrimMom = true
   			counterCX = false
    		     end
   		  end
   	       end
   	    end
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
   for sN, _ in lume.orderedIter(species) do
      if sN ~= self.name then
         -- Flags for couplingMoments, boundary corrections, star moments,
         -- self primitive moments, cross primitive moments.
         self.momentFlags[5][sN] = false
      end

      for sO, _ in lume.orderedIter(species) do
         -- Allocate space for this species' cross-primitive moments
         -- only if some other species collides with it.
         if (sN ~= sO) and (self.collPairs[sN][sO].on or self.collPairs[sO][sN].on) then
            otherNm = string.gsub(sO .. sN, self.name, "")
            if self.uCross[otherNm] == nil then
               self.uCross[otherNm] = self:allocVectorMoment(self.vdim)
            end
            if self.vtSqCross[otherNm] == nil then
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
            onGrid   = self.confGrid,
            basis    = self.confBasis,
            evaluate = function(t,xn) return 0.0 end, -- Function is set below.
            onGhosts = false,
         }
      end
      for sN, _ in lume.orderedIter(species) do
         if sN ~= self.name then
            -- Sixth moment flag is to indicate if spatially varying collisionality has been computed.
            self.momentFlags[6][sN] = false
         end

         for sO, _ in lume.orderedIter(species) do
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

function VlasovSpecies:setActiveRKidx(rkIdx)
   self.activeRKidx = rkIdx
end

function VlasovSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self:setActiveRKidx(inIdx)
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   -- Accumulate functional Maxwell fields (if needed).
   local emField         = emIn[1]:rkStepperFields()[inIdx]
   local emExternalField = emIn[2]:rkStepperFields()[1]
   local totalEmField    = self.totalEmField
   totalEmField:clear(0.0)

   local qbym = self.charge/self.mass

   if emField and self.computePlasmaB then totalEmField:accumulate(qbym, emField) end
   if emExternalField then totalEmField:accumulate(qbym, emExternalField) end

   -- If external force present (gravity, body force, etc.) accumulate it to electric field.
   if self.hasExtForce then
      local vExtForce = self.vExtForce
      self.evalVlasovExtForce:advance(tCurr, {}, {vExtForce})

      -- Need to barrier over the shared communicator before accumulating force onto electric field.
      Mpi.Barrier(self.grid:commSet().sharedComm)

      -- Analogous to the current, the external force only gets accumulated onto the electric field.
      local vItr, eItr   = vExtForce:get(1), totalEmField:get(1)
      local vIdxr, eIdxr = vExtForce:genIndexer(), totalEmField:genIndexer()

      for idx in totalEmField:localRangeIter() do
         vExtForce:fill(vIdxr(idx), vItr)
         totalEmField:fill(eIdxr(idx), eItr)
         for i = 1, vExtForce:numComponents() do eItr[i] = eItr[i]+vItr[i] end
      end
   end

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solver:advance(tCurr, {fIn, totalEmField, emField}, {fRhsOut})
   else
      fRhsOut:clear(0.0)    -- No RHS.
   end

   -- Perform the collision update.
   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
         c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
         c:advance(tCurr, fIn, species, fRhsOut)   -- 'species' needed for cross-species collisions.
      end
   end
   if self.evolveSources then
      for _, src in pairs(self.sources) do
         src:advance(tCurr, fIn, species, fRhsOut)
      end
   end

   if self.projSrc and self.evolveSources then
      -- add source it to the RHS
      -- Barrier over shared communicator before accumulate
      Mpi.Barrier(self.grid:commSet().sharedComm)
      fRhsOut:accumulate(self.sourceTimeDependence(tCurr), self.fSource)
   end

   -- Save boundary fluxes for diagnostics.
   if self.hasNonPeriodicBc and self.boundaryFluxDiagnostics then
      for _, bc in ipairs(self.boundaryConditions) do
         bc:storeBoundaryFlux(tCurr, outIdx, fRhsOut)
      end
   end
end

function VlasovSpecies:createDiagnostics()
   local function contains(table, element)
     for _, value in pairs(table) do
       if value == element then
         return true
       end
     end
     return false
   end

    if self.fSource then
      self.numDensitySrc = self:allocMoment()
      self.momDensitySrc = self:allocVectorMoment(self.vdim)
      self.ptclEnergySrc = self:allocMoment()
      self.fiveMomentsCalc:advance(0.0, {self.fSource}, {self.numDensitySrc, self.momDensitySrc, self.ptclEnergySrc})
    end
   -- Create updater to compute volume-integrated moments
   -- function to check if integrated moment name is correct.
   local function isIntegratedMomentNameGood(nm)
      if nm == "intM0" or nm == "intM1i" or nm == "intM2Flow" 
         or nm == "intM2Thermal" or nm == "intM2" or nm == "intL2" then
         return true
      end
      return false
   end

   local numCompInt = {}
   numCompInt["intM0"]        = 1
   numCompInt["intM1i"]       = self.vdim
   numCompInt["intM2Flow"]    = 1
   numCompInt["intM2Thermal"] = 1
   numCompInt["intM2"]        = 1
   numCompInt["intL2"]        = 1

   self.diagnosticIntegratedMomentFields   = { }
   self.diagnosticIntegratedMomentUpdaters = { } 
   -- Allocate space to store integrated moments and create integrated moment updaters.
   local function allocateDiagnosticIntegratedMoments(intMoments, bc, timeIntegrate)
      local label = ""
      local phaseGrid = self.grid
      local confGrid = self.confGrid
      if bc then
         label = bc:label()
         phaseGrid = bc:getBoundaryGrid()
         confGrid = bc:getConfBoundaryGrid()
      end
      local timeIntegrate = xsys.pickBool(timeIntegrate, false)
      for i, mom in ipairs(intMoments) do
         if isIntegratedMomentNameGood(mom) then
            self.diagnosticIntegratedMomentFields[mom..label] = DataStruct.DynVector {
               numComponents = numCompInt[mom],
            }
            local intCalc = Updater.CartFieldIntegratedQuantCalc {
                  onGrid        = self.confGrid,
                  basis         = self.confBasis,
                  numComponents = numCompInt[mom],
                  quantity      = "V",
                  timeIntegrate = timeIntegrate,
               }
            if mom == "intL2" then
               self.diagnosticIntegratedMomentUpdaters[mom..label] = Updater.CartFieldIntegratedQuantCalc {
                  onGrid        = self.grid,
                  basis         = self.basis,
                  numComponents = numCompInt[mom],
                  quantity      = "V2",
                  timeIntegrate = timeIntegrate,
               }
            else
               self.diagnosticIntegratedMomentUpdaters[mom..label] = intCalc
            end
         else
            assert(false, string.format("Error: integrated moment %s not valid", mom..label))
         end
      end
   end

   allocateDiagnosticIntegratedMoments(self.diagnosticIntegratedMoments)

   if self.hasNonPeriodicBc and self.boundaryFluxDiagnostics then
      for _, bc in ipairs(self.boundaryConditions) do
         bc:initBcDiagnostics(self.cdim)
         allocateDiagnosticIntegratedMoments(self.diagnosticIntegratedBoundaryFluxMoments, bc, true)
      end
   end

   -- Function to check if moment name is correct.
   local function isMomentNameGood(nm)
      return Updater.DistFuncMomentCalc:isMomentNameGood(nm)
   end
   -- weakMoments are diagnostics computed with weak binary operations.
   -- Check if diagnostic name is correct.
   local function isWeakMomentNameGood(nm)
      return nm == "u" or nm == "vtSq" or nm == "uCross" or nm == "vtSqCross" 
          or nm == "M2Flow" or nm == "M2Thermal" 
   end

   local numComp        = {}
   numComp["M0"]        = 1
   numComp["M1i"]       = self.vdim
   numComp["M2ij"]      = self.vdim*(self.vdim+1)/2
   numComp["M2"]        = 1
   numComp["M3i"]       = self.vdim
   numComp["u"]         = self.vdim
   numComp["vtSq"]      = 1
   numComp["uCross"]    = self.vdim
   numComp["vtSqCross"] = 1
   numComp["M2Flow"]    = 1
   numComp["M2Thermal"] = 1

   self.diagnosticMomentFields   = { }
   self.diagnosticMomentUpdaters = { } 
   self.diagnosticWeakMoments    = { }
   self.diagnosticWeakBoundaryFluxMoments = { }
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
   self.confPhaseMult = Updater.CartFieldBinOp {
      onGrid     = self.grid,
      weakBasis  = self.basis,
      fieldBasis = self.confBasis,
      operation  = "Multiply",
   }

   -- Sort moments into diagnosticWeakMoments.
   local function organizeDiagnosticMoments(moments, weakMoments, integratedMoments)
      -- At beginning, all moment names are in the 'moments' list.
      -- We want to remove the weak moments and put them in the 'weakMoments' list
      for i, mom in pairs(moments) do
         if isWeakMomentNameGood(mom) then
            -- Remove moment name from moments list, and add it to weakMoments list.
            if mom == "uCross" then
               for nm, _ in pairs(self.uCross) do
                  -- Create one diagnostic for each cross velocity (used in collisions).
                  table.insert(weakMoments, mom .. "-" .. nm)
               end
            elseif mom == "vtSqCross" then
               for nm, _ in pairs(self.vtSqCross) do
                  -- Create one diagnostic for each cross temperature (used in collisions).
                  table.insert(weakMoments, mom .. "-" .. nm)
               end
            else
               table.insert(weakMoments, mom)
            end
            moments[i] = nil
         end
      end

      -- Make sure we have moment updaters/fields needed to compute integrated moments.
      -- Note: this could result in extra moments being written out if they were not
      -- already requested.
      for i, mom in ipairs(integratedMoments) do
         -- integrated M0
         if mom == "intM0" then
            if not contains(moments, "M0") then
               table.insert(moments, "M0")
            end
         end
         -- integrated M1i
         if mom == "intM1i" then
            if not contains(moments, "M1i") then
               table.insert(moments, "M1i")
            end
         end
         -- integrated M2
         if mom == "intM2" then
            if not contains(moments, "M2") then
               table.insert(moments, "M2")
            end
         end
         -- integrated M2Flow
         if mom == "intM2Flow" then
            if not contains(weakMoments, "M2Flow") then
               table.insert(weakMoments, "M2Flow")
            end
         end
         -- integrated M2Thermal
         if mom == "intM2Thermal" then
            if not contains(weakMoments, "M2Thermal") then
               table.insert(weakMoments, "M2Thermal")
            end
         end
      end
   
      -- Make sure we have the updaters needed to calculate all the weak moments.
      -- Note: this could result in extra moments being written out if they were not
      -- already requested.
      for i, mom in ipairs(weakMoments) do
         -- All weak moments require M0 = density.
         if not contains(moments, "M0") then
            table.insert(moments, "M0")
         end
   
         if mom == "u" then
            if not contains(moments, "M1i") then
               table.insert(moments, "M1i")
            end
         end
         if mom == "vtSq" then
            if not contains(moments, "M2") then
               table.insert(moments, "M2")
            end
            if not contains(weakMoments, "M2Thermal") then
               table.insert(weakMoments, "M2Thermal")
            end
            if not contains(weakMoments, "M2Flow") then
               table.insert(weakMoments, "M2Flow")
            end
            -- M1i and u are needed for M2Flow
            if not contains(moments, "M1i") then
               table.insert(moments, "M1i")
            end
            if not contains(weakMoments, "u") then
               table.insert(weakMoments, "u")
            end
         end
         if mom == "M2Flow" then -- = n*u^2 = M1.u
            if not contains(moments, "M1i") then
               table.insert(moments, "M1i")
            end
            if not contains(weakMoments, "u") then
               table.insert(weakMoments, "u")
            end
         end
         if mom == "M2Thermal" then -- = VDIM*n*vtSq = M2 - M2Flow
            if not contains(moments, "M2") then
               table.insert(moments, "M2")
            end
            if not contains(weakMoments, "M2Flow") then
               table.insert(weakMoments, "M2Flow")
            end
            -- M1i and u are needed for M2Flow
            if not contains(moments, "M1i") then
               table.insert(moments, "M1i")
            end
            if not contains(weakMoments, "u") then
               table.insert(weakMoments, "u")
            end
         end
      end
   end

   -- Allocate space to store moments and create moment updater.
   local function allocateDiagnosticMoments(moments, weakMoments, bc)
      local label = ""
      local phaseGrid = self.grid
      local confGrid = self.confGrid
      if bc then
         label = bc:label()
         phaseGrid = bc:getBoundaryGrid()
         confGrid = bc:getConfBoundaryGrid()
      end

      for i, mom in pairs(moments) do
         if isMomentNameGood(mom) then
            self.diagnosticMomentFields[mom..label] = DataStruct.Field {
               onGrid        = confGrid,
               numComponents = self.confBasis:numBasis()*numComp[mom],
               ghost         = {1, 1},
               metaData = {
                  polyOrder = self.basis:polyOrder(),
                  basisType = self.basis:id(),
                  charge = self.charge,
                  mass = self.mass,
               },
            }
            self.diagnosticMomentUpdaters[mom..label] = Updater.DistFuncMomentCalc {
               onGrid     = phaseGrid,
               phaseBasis = self.basis,
               confBasis  = self.confBasis,
               moment     = mom,
            }
         else
            assert(false, string.format("Error: moment %s not valid", mom..label))
         end
      end

      for i, mom in ipairs(weakMoments) do
         self.diagnosticMomentFields[mom..label] = DataStruct.Field {
            onGrid        = confGrid,
            numComponents = self.confBasis:numBasis()*numComp[mom],
            ghost         = {1, 1},
            metaData = {
               polyOrder = self.basis:polyOrder(),
               basisType = self.basis:id(),
               charge = self.charge,
               mass = self.mass,
            },	    	    
         }

         self.diagnosticMomentUpdaters[mom..label] = {}
         -- Weak moments do not have their own MomentCalc updaters. 
         -- Instead they are computed via weak division of two other moments.
         -- Here we set up custom advance methods for each weak moment, 
         -- which call various weak ops.
         if mom == "u" then
            self.diagnosticMomentUpdaters["u"..label].advance = function (self, tm)
               if self.diagnosticMomentUpdaters["u"..label].tCurr == tm then return end -- return if already computed for this tm

               -- compute dependencies if not already computed: M0, M1i
               if self.diagnosticMomentUpdaters["M0"..label].tCurr ~= tm then
                  local fIn = self:rkStepperFields()[1]
                  if bc then
                     fIn = bc:getBoundaryFluxRate()
                  end
                  self.diagnosticMomentUpdaters["M0"..label]:advance(tm, {fIn}, {self.diagnosticMomentFields["M0"..label]})
               end
               if self.diagnosticMomentUpdaters["M1i"..label].tCurr ~= tm then
                  local fIn = self:rkStepperFields()[1]
                  if bc then
                     fIn = bc:getBoundaryFluxRate()
                  end
                  self.diagnosticMomentUpdaters["M1i"..label]:advance(tm, {fIn}, {self.diagnosticMomentFields["M1i"..label]})
               end

               -- do weak ops
               self.weakDivision:advance(tm, {self.diagnosticMomentFields["M0"..label], self.diagnosticMomentFields["M1i"..label]}, {self.diagnosticMomentFields["u"..label]})

               self.diagnosticMomentUpdaters["u"..label].tCurr = tm -- mark as complete for this tm
            end
         elseif mom == "M2Flow" then
            self.diagnosticMomentUpdaters["M2Flow"..label].advance = function (self, tm)
               if self.diagnosticMomentUpdaters["M2Flow"..label].tCurr == tm then return end -- return if already computed for this tm

               -- compute dependencies if not already computed: u (M1i will be computed via u if necessary)
               if self.diagnosticMomentUpdaters["u"..label].tCurr ~= tm then
                  self.diagnosticMomentUpdaters["u"..label].advance(self, tm)
               end

               -- do weak ops
               -- M2Flow = M1i.u
               self.weakDotProduct:advance(tm, {self.diagnosticMomentFields["M1i"..label], self.diagnosticMomentFields["u"..label]}, 
                                               {self.diagnosticMomentFields["M2Flow"..label]})

               self.diagnosticMomentUpdaters["M2Flow"..label].tCurr = tm -- mark as complete for this tm
            end
         elseif mom == "M2Thermal" then
            self.diagnosticMomentUpdaters["M2Thermal"..label].advance = function (self, tm)
               if self.diagnosticMomentUpdaters["M2Thermal"..label].tCurr == tm then return end -- return if already computed for this tm

               -- compute dependencies if not already computed: M2, M2Flow 
               if self.diagnosticMomentUpdaters["M2"..label].tCurr ~= tm then
                  local fIn = self:rkStepperFields()[1]
                  if bc then
                     fIn = bc:getBoundaryFluxRate()
                  end
                  self.diagnosticMomentUpdaters["M2"..label]:advance(tm, {fIn}, {self.diagnosticMomentFields["M2"..label]})
               end
               if self.diagnosticMomentUpdaters["M2Flow"..label].tCurr ~= tm then
                  self.diagnosticMomentUpdaters["M2Flow"..label].advance(self, tm)
               end

               -- do weak ops
               -- M2Thermal = VDIM*M0*vtSq = M2 - M2Flow               
               self.diagnosticMomentFields["M2Thermal"..label]:combine(1.0, self.diagnosticMomentFields["M2"..label], -1.0, self.diagnosticMomentFields["M2Flow"..label])

               self.diagnosticMomentUpdaters["M2Thermal"..label].tCurr = tm -- mark as complete for this tm
            end
         elseif mom == "vtSq" then 
            self.diagnosticMomentUpdaters["vtSq"..label].advance = function (self, tm)
               if self.diagnosticMomentUpdaters["vtSq"..label].tCurr == tm then return end -- return if already computed for this tm

               -- compute dependencies if not already computed: M0, M2Thermal
               if self.diagnosticMomentUpdaters["M0"..label].tCurr ~= tm then
                  local fIn = self:rkStepperFields()[1]
                  if bc then
                     fIn = bc:getBoundaryFluxRate()
                  end
                  self.diagnosticMomentUpdaters["M0"..label]:advance(tm, {fIn}, {self.diagnosticMomentFields["M0"..label]})
               end
               if self.diagnosticMomentUpdaters["M2Thermal"..label].tCurr ~= tm then
                  self.diagnosticMomentUpdaters["M2Thermal"..label].advance(self, tm)
               end

               -- do weak ops
               -- vtSq = M2Thermal/(M0*VDIM)
               self.weakDivision:advance(tm, {self.diagnosticMomentFields["M0"..label], self.diagnosticMomentFields["M2Thermal"..label]}, {self.diagnosticMomentFields["vtSq"..label]})
               self.diagnosticMomentFields["vtSq"..label]:scale(1.0/self.vdim)

               self.diagnosticMomentUpdaters["vtSq"..label].tCurr = tm -- mark as complete for this tm
            end
         elseif string.find(mom, "uCross") then
            self.diagnosticMomentUpdaters["uCross"..label].advance = function (self, tm)
               otherNm = string.gsub(mom, "uCross%-", "")
               local uCrossOther = self.uCross[otherNm]
               if bc then
                  uCrossOther = bc:evalOnConfBoundary(self.uCross[otherNm])
               end
               self.diagnosticMomentFields[mom..label]:copy(uCrossOther)
            end
         elseif string.find(mom, "vtSqCross") then
            self.diagnosticMomentUpdaters["vtSqCross"..label].advance = function (self, tm)
               otherNm = string.gsub(mom, "vtSqCross%-", "")
               local vtSqCrossOther = self.vtSqCross[otherNm]
               if bc then
                  vtSqCrossOther = bc:evalOnConfBoundary(self.vtSqCross[otherNm])
               end
               self.diagnosticMomentFields[mom..label]:copy(vtSqCrossOther)
            end
         end
      end
   end

   organizeDiagnosticMoments(self.diagnosticMoments, self.diagnosticWeakMoments, self.diagnosticIntegratedMoments)
   allocateDiagnosticMoments(self.diagnosticMoments, self.diagnosticWeakMoments)

   if self.hasNonPeriodicBc and self.boundaryFluxDiagnostics then
      organizeDiagnosticMoments(self.diagnosticBoundaryFluxMoments, self.diagnosticWeakBoundaryFluxMoments, self.diagnosticIntegratedBoundaryFluxMoments)
      for _, bc in ipairs(self.boundaryConditions) do
         allocateDiagnosticMoments(self.diagnosticBoundaryFluxMoments, self.diagnosticWeakBoundaryFluxMoments, bc)
      end
   end
end

-- Function to compute integrated moments.
function VlasovSpecies:calcDiagnosticIntegratedMoments(tm)
   local fIn = self:rkStepperFields()[1]

   local function computeIntegratedMoments(intMoments, fIn, label)
      local label = label or ""
      for i, mom in ipairs(intMoments) do
         if mom == "intM0" then
            self.diagnosticMomentUpdaters["M0"..label]:advance(
               tm, {fIn}, {self.diagnosticMomentFields["M0"..label]})
            self.diagnosticIntegratedMomentUpdaters[mom..label]:advance(
               tm, {self.diagnosticMomentFields["M0"..label]}, {self.diagnosticIntegratedMomentFields[mom..label]})
         elseif mom == "intM1i" then
            self.diagnosticMomentUpdaters["M1i"..label]:advance(
               tm, {fIn}, {self.diagnosticMomentFields["M1i"..label]})
            self.diagnosticIntegratedMomentUpdaters[mom..label]:advance(
               tm, {self.diagnosticMomentFields["M1i"..label]}, {self.diagnosticIntegratedMomentFields[mom..label]})
         elseif mom == "intM2" then
            self.diagnosticMomentUpdaters["M2"..label]:advance(
               tm, {fIn}, {self.diagnosticMomentFields["M2"..label]})
            self.diagnosticIntegratedMomentUpdaters[mom..label]:advance(
               tm, {self.diagnosticMomentFields["M2"..label]}, {self.diagnosticIntegratedMomentFields[mom..label]})
         elseif mom == "intM2Flow" then
            self.diagnosticMomentUpdaters["M2Flow"..label].advance(self, tm)
            self.diagnosticIntegratedMomentUpdaters[mom..label]:advance(
               tm, {self.diagnosticMomentFields["M2Flow"..label]}, {self.diagnosticIntegratedMomentFields[mom..label]})
         elseif mom == "intM2Thermal" then
            self.diagnosticMomentUpdaters["M2Thermal"..label].advance(self, tm)
            self.diagnosticIntegratedMomentUpdaters[mom..label]:advance(
               tm, {self.diagnosticMomentFields["M2Thermal"..label]}, {self.diagnosticIntegratedMomentFields[mom..label]})
         elseif mom == "intL2" then
            self.diagnosticIntegratedMomentUpdaters[mom]:advance(
               tm, {self.distf[1]}, {self.diagnosticIntegratedMomentFields[mom]})
         end
      end
   end

   computeIntegratedMoments(self.diagnosticIntegratedMoments, fIn)
   
   if self.hasNonPeriodicBc and self.boundaryFluxDiagnostics then
      for _, bc in ipairs(self.boundaryConditions) do
         computeIntegratedMoments(self.diagnosticIntegratedBoundaryFluxMoments, bc:getBoundaryFluxRate(), bc:label())
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
   local numBasis = self.basis:numBasis()
   local velIdx = Lin.IntVec(self.ndim)
   velIdx[1] = 1
   for d = 1, self.vdim do
      velIdx[d + 1] = idxIn[self.cdim + d]
   end
   local exIdxr = self.externalBCFunction:genIndexer()
   local externalBCFunction = self.externalBCFunction:get(exIdxr(velIdx))
   if velIdx[1] ~= 0 and velIdx[1] ~= self.grid:numCells(2) + 1 then
      for i = 1, numBasis do
	 fOut[i] = 0
         for j = 1, numBasis do
            fOut[i] = fOut[i] + fIn[j]*externalBCFunction[(i - 1)*numBasis + j]
         end
      end
   end
   return fOut
end

function VlasovSpecies:calcExternalBC()
   tbl = self.tbl
   local externalBC = assert(tbl.externalBC, "VlasovSpecies: Must define externalBC parameters")
   local lower = Lin.Vec(self.grid:ndim())
   local upper = Lin.Vec(self.grid:ndim())
   local cells = Lin.Vec(self.grid:ndim())
   local GridConstructor = Grid.RectCart
   local coordinateMap = {}
   if self.coordinateMap then
      for d = 1, self.cdim do
         table.insert(coordinateMap, function (z) return z end)
      end
      for d = 1, self.vdim do
         table.insert(coordinateMap, self.coordinateMap[d])
      end
      GridConstructor = Grid.NonUniformRectCart
   end
   for d = 1, self.cdim do
      lower[d] = 0
      upper[d] = 1
      cells[d] = 1
   end
   for d = 1, self.vdim do
      lower[d + self.cdim] = self.grid:lower(d + self.cdim)
      upper[d + self.cdim] = self.grid:upper(d + self.cdim)
      cells[d + self.cdim] = self.grid:numCells(d + self.cdim)
   end
   local grid = GridConstructor {
      lower = lower,
      upper = upper,
      cells = cells,
      periodicDirs = {},
      mappings = coordinateMap,
   }
   self.externalBCFunction = DataStruct.Field {
      onGrid = grid,
      numComponents = self.basis:numBasis()*self.basis:numBasis(),
      metaData = {
         polyOrder = self.basis:polyOrder(),
         basisType = self.basis:id(),
         charge = self.charge,
         mass = self.mass,
      },
   }
   -- Calculate Bronold and Fehske reflection function coefficients
   local evaluateBronold = Updater.EvaluateBronoldFehskeBC {
      onGrid = grid,
      basis = self.basis,
      cdim = self.cdim,
      vdim = self.vdim,
      electronAffinity = externalBC.electronAffinity,
      effectiveMass = externalBC.effectiveMass,
      elemCharge = externalBC.elemCharge,
      me = externalBC.electronMass,
   }
   evaluateBronold:advance(0.0, {}, {self.externalBCFunction})
end

function VlasovSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- Need to wrap member functions so that self is passed.
   local function bcAbsorbFunc(...)  return self:bcAbsorbFunc(...) end
   local function bcCopyFunc(...)    return self:bcCopyFunc(...) end
   local function bcOpenFunc(...)    return self:bcOpenFunc(...) end
   local function bcReflectFunc(...) return self:bcReflectFunc(...) end
   local function bcExternFunc(...)  return self:bcExternFunc(...) end

   local vdir = dir + self.cdim

   if type(bcType) == "function" then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcCopyFunc }, "pointwise", bcType))
   elseif bcType == SP_BC_ABSORB then
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
   elseif bcType == SP_BC_EXTERN then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, vdir, edge,
				      { bcExternFunc }, "flip"))
   elseif bcType == SP_BC_ZEROFLUX then
      table.insert(self.zeroFluxDirections, dir)
   else
      assert(false, "VlasovSpecies: Unsupported BC type!")
   end
end

function VlasovSpecies:calcCouplingMoments(tCurr, rkIdx, species)

   local tmStart = Time.clock()
   -- Compute moments needed in coupling to fields and collisions.
   local fIn = self:rkStepperFields()[rkIdx]
   if self.needSelfPrimMom and
      lume.any({unpack(self.momentFlags,1,4)},function(x) return x==false end) then -- No need to recompute if already computed.
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
         -- Barrier over shared communicator before combine
         Mpi.Barrier(self.grid:commSet().sharedComm)
         self.thermalEnergyDensity:combine( 1.0/self.vdim, self.ptclEnergy,
                                           -1.0/self.vdim, self.kineticEnergyDensity )
         self.confDiv:advance(tCurr, {self.numDensity, self.thermalEnergyDensity}, {self.vtSqSelf})
      end
      -- Indicate that moments, boundary corrections, star moments
      -- and self-primitive moments have been computed.
      for iF=1,4 do self.momentFlags[iF] = true end
   elseif self.momentFlags[1]==false then -- No need to recompute if already computed.
      if self.computePlasmaB then
         self.momDensityCalc:advance(tCurr, {fIn}, { self.momDensity })
      else
         self.numDensityCalc:advance(tCurr, {fIn}, { self.numDensity })
      end
      -- Indicate that the coupling moment has been computed.
      self.momentFlags[1] = true
   end

   -- For Ionization.
   if self.calcReactRate then
      local neuts = species[self.neutNmIz]
      if lume.any({unpack(neuts.momentFlags,1,4)},function(x) return x==false end) then
         -- Neutrals haven't been updated yet, so we need to compute their moments and primitive moments. 
         neuts:calcCouplingMoments(tCurr, rkIdx, species)
      end
      local neutM0   = neuts:fluidMoments()[1]
      local neutU    = neuts:selfPrimitiveMoments()[1]
      local neutVtSq = neuts:selfPrimitiveMoments()[2]
      
      if tCurr == 0.0 then
	 species[self.name].collisions[self.collNmIoniz].collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      end
      species[self.name].collisions[self.collNmIoniz].collisionSlvr:advance(tCurr, {neutM0, neutVtSq, self.vtSqSelf}, {self.voronovReactRate})
      species[self.name].collisions[self.collNmIoniz].calcIonizationTemp:advance(tCurr, {self.vtSqSelf}, {self.vtSqIz})

      self.calcMaxwell:advance(tCurr, {self.numDensity, neutU, self.vtSqIz}, {self.fMaxwellIz})
            
      self.numDensityCalc:advance(tCurr, {self.fMaxwellIz}, {self.m0fMax})
      self.confDiv:advance(tCurr, {self.m0fMax, self.numDensity}, {self.m0mod})
      self.confPhaseMult:advance(tCurr, {self.m0mod, self.fMaxwellIz}, {self.fMaxwellIz})
   end

   -- For charge exchange.
   if self.calcCXSrc then
      -- Calculate Vcx*SigmaCX.
      local neuts = species[self.neutNmCX]
      if lume.any({unpack(neuts.momentFlags,1,4)},function(x) return x==false end) then
         -- Neutrals haven't been updated yet, so we need to compute their moments and primitive moments. 
         neuts:calcCouplingMoments(tCurr, rkIdx, species)
      end
      local m0       = neuts:fluidMoments()[1]
      local neutU    = neuts:selfPrimitiveMoments()[1]
      local neutVtSq = neuts:selfPrimitiveMoments()[2]
      
      species[self.neutNmCX].collisions[self.collNmCX].collisionSlvr:advance(tCurr, {m0, self.uSelf, neutU, self.vtSqSelf, neutVtSq}, {self.vSigmaCX})
   end

   self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
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
   for i, mom in pairs(self.diagnosticMoments) do
      tm = tm + self.diagnosticMomentUpdaters[mom].totalTime
   end
   return tm
end

function VlasovSpecies:getVoronovReactRate()
   return self.voronovReactRate
end

function VlasovSpecies:getFMaxwellIz()
   return self.fMaxwellIz
end

function VlasovSpecies:getSrcCX()
   return self.srcCX
end

-- Please test this for higher than 1x1v... (MF: JJ?).
function VlasovSpecies:Maxwellian(xn, n0, T0, vdnIn)
   local vdn = vdnIn or {0, 0, 0}
   local vt2 = T0/self.mass
   local v2 = 0.0
   for d = self.cdim+1, self.cdim+self.vdim do
     v2 = v2 + (xn[d] - vdnIn[d-self.cdim])^2
   end
   return n0 / math.sqrt(2*math.pi*vt2)^self.vdim * math.exp(-v2/(2*vt2))
end

return VlasovSpecies
