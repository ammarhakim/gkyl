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
local Projection     = require "App.Projection"
local Proto          = require "Lib.Proto"
local Source         = require "App.Sources.VmSource"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local VlasovEq       = require "Eq.Vlasov"
local DiagsApp       = require "App.Diagnostics.SpeciesDiagnostics"
local VlasovDiags    = require "App.Diagnostics.VlasovDiagnostics"
local BasicBC        = require ("App.BCs.VlasovBasic").VlasovBasic
local BCsBase        = require "App.BCs.BCsBase"
local ffi            = require "ffi"
local xsys           = require "xsys"
local lume           = require "Lib.lume"

local VlasovSpecies = Proto(KineticSpecies)

local SP_BC_EXTERN  = 4
local SP_BC_RECYCLE = 7
VlasovSpecies.bcExternal = SP_BC_EXTERN     -- Load external BC file.
VlasovSpecies.bcRecycle  = SP_BC_RECYCLE

-- ............. Backwards compatible treatment of BCs .....................--
-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB    = 1
local SP_BC_REFLECT   = 3
local SP_BC_COPY      = 5
-- AHH: This was 2 but seems that is unstable. So using plain copy.
local SP_BC_OPEN      = SP_BC_COPY
local SP_BC_ZEROFLUX  = 6
VlasovSpecies.bcCopy      = SP_BC_COPY       -- Copy stuff.
VlasovSpecies.bcAbsorb    = SP_BC_ABSORB     -- Absorb all particles.
VlasovSpecies.bcOpen      = SP_BC_OPEN       -- Zero gradient.
VlasovSpecies.bcReflect   = SP_BC_REFLECT    -- Specular reflection.
VlasovSpecies.bcZeroFlux  = SP_BC_ZEROFLUX

function VlasovSpecies:makeBcApp(bcIn, dir, edge)
   local bcOut
   if type(bcIn) == "function" then
      bcOut = BasicBC{kind="function", bcFunction=bcIn, diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_COPY then
      print("VlasovSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="copy", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_ABSORB then
      print("VlasovSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="absorb", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_OPEN then
      print("VlasovSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      -- AHH: open seems unstable. So using plain copy.
      bcOut = BasicBC{kind="copy", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_REFLECT then
      print("VlasovSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="reflect", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_ZEROFLUX or bcIn.tbl.kind=="zeroFlux" then
      bcOut = "zeroFlux"
      table.insert(self.zeroFluxDirections, dir)
   end
   return bcOut
end

-- ............. End of backwards compatibility for BCs .....................--

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

   -- vFlux used for selecting which type of numerical flux function to use in velocity space
   -- defaults to "penalty" in Eq object, supported options: "penalty," "recovery," and "upwind"
   -- Note: "recovery" and "upwind" only supported by p=1 Serendipity and p=2 Tensor
   -- Note: only used for DG Vlasov-Maxwell.
   self.numVelFlux = tbl.vFlux
end

function VlasovSpecies:allocMomCouplingFields()
   return { currentDensity = self:allocVectorMoment(self.vdim) }
end


function VlasovSpecies:createSolver(field, externalField)
   -- Run the KineticSpecies 'createSolver()' to initialize the collisions solver
   VlasovSpecies.super.createSolver(self, field, externalField)

   if externalField.geo.bHat then -- for genGeo with neutrals
      self._isGenGeo = true
      self.bHat = externalField.geo.bHat
      self.gxx = externalField.geo.gxx
      self.gxy = externalField.geo.gxy
      self.gxz = externalField.geo.gxz
      self.gyy = externalField.geo.gyy
      self.gyz = externalField.geo.gyz
      self.gzz = externalField.geo.gzz
      self.jacobGeo = externalField.geo.jacobGeo
      self.tanVecComp = externalField.geo.tanVecComp

      self.calcAlphaGeo = Updater.AlphaGenGeoCalc {
	 confGrid = self.confGrid,
	 onGrid = self.grid,
	 confBasis = self.confBasis,
	 phaseBasis = self.basis,
      }
      local metadata = {polyOrder = self.basis:polyOrder(),
                     basisType = self.basis:id(),
                     charge    = self.charge,
                     mass      = self.mass}
      self.alphaGeo = self:allocCartField(self.grid,self.cdim*self.basis:numBasis(),{self.nGhost,self.nGhost},metaData)
      
      -- Calculate alphaGeo here. Use an Updater.
      self.calcAlphaGeo:advance(0.0, {self.tanVecComp, self.gxx, self.gxy, self.gxz, self.gyy, self.gyz, self.gzz, self.jacobGeo}, {self.alphaGeo})
   end
      
   local plasmaE, plasmaB = field:hasEB()
   local extHasE, extHasB = externalField:hasEB()

   local hasE = plasmaE or extHasE
   local hasB = plasmaB or extHasB
   
   -- External forces are accumulated to the electric field part of
   -- totalEmField
   if self.hasExtForce then
      hasE = true
      hasB = true
   end

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
   for d = 1, self.vdim do table.insert(self.zeroFluxDirections, self.cdim+d) end

   self.solver = Updater.HyperDisCont {
      onGrid             = self.grid,
      basis              = self.basis,
      cfl                = self.cfl,
      equation           = self.equation,
      zeroFluxDirections = self.zeroFluxDirections,
   }

   -- Create updaters to compute various moments.
   self.numDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  moment     = "M0",
   }
   self.momDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  moment     = "M1i",
   }
   self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  moment     = "M2",
   }
   self.M2ijCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis = self.confBasis,
      phaseBasis = self.basis,  moment    = "M2ij",
   }
   self.M3iCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis = self.confBasis,
      phaseBasis = self.basis,  moment    = "M3i",
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
      onGrid      = self.grid,
      phaseBasis  = self.basis,
      confGrid    = self.confGrid,
      confBasis   = self.confBasis,
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

      self.zIdx = Lin.IntVec(self.grid:ndim())

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

   -- Create an updater for volume integrals. Used by diagnostics.
   self.volIntegral = {
      scalar = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.confGrid,   numComponents = 1,
         basis  = self.confBasis,  quantity      = "V",
      },
      vector = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.confGrid,   numComponents = self.vdim,
         basis  = self.confBasis,  quantity      = "V",
      },
   }

   -- Create species source solvers.
   for _, src in lume.orderedIter(self.sources) do src:createSolver(self, externalField) end

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
			counterIz_elc         = false
		     elseif self.name==species[sN].collisions[collNm].neutNm and counterIz_neut then
			self.needSelfPrimMom = true
			if self.vdim == 3 then self.needGkMom = true end
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
			species[self.neutNmCX].needSelfPrimMom = true
			if self.vdim == 3 and (not self.needGkMom) then self.needGkMom = true end
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
      if self.needGkMom then
	 self.ptclEnergyNorm = self:allocMoment() -- for M2/M1
	 self.uPar = self:allocMoment()
	 self.uParSq = self:allocMoment()
	 self.vtSqGk = self:allocMoment()
	 self.confMult = Updater.CartFieldBinOp {
	    onGrid    = self.confGrid,
	    weakBasis = self.confBasis,
	    operation = "Multiply",
	 }
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

   -- Initialize the BC cross-coupling interactions.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:initCrossSpeciesCoupling(species) end

end

function VlasovSpecies:setActiveRKidx(rkIdx) self.activeRKidx = rkIdx end

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
      self.solver:advance(tCurr, {fIn, totalEmField, emField, self.alphaGeo}, {fRhsOut})
   else
      fRhsOut:clear(0.0)    -- No RHS.
   end

   -- Perform the collision update.
   for _, c in pairs(self.collisions) do
      c.collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      c:advance(tCurr, fIn, species, fRhsOut)   -- 'species' needed for cross-species collisions.
   end

   for _, src in lume.orderedIter(self.sources) do src:advance(tCurr, fIn, species, fRhsOut) end
   
   for _, bc in pairs(self.nonPeriodicBCs) do
      bc:storeBoundaryFlux(tCurr, outIdx, fRhsOut)   -- Save boundary fluxes.
   end
end

function VlasovSpecies:advanceCrossSpeciesCoupling(tCurr, species, emIn, inIdx, outIdx)
   -- Perform some operations after the updates have been computed, but before
   -- the combine RK (in PlasmaOnCartGrid) is called.

   for _, bc in pairs(self.nonPeriodicBCs) do bc:advanceCrossSpeciesCoupling(tCurr, species, outIdx) end
end

function VlasovSpecies:createDiagnostics(field)
   -- Run the KineticSpecies 'createDiagnostics()' (e.g. to create divideByJacobGeo()).
   VlasovSpecies.super.createDiagnostics(self, field)

   -- Create this species' diagnostics.
   if self.tbl.diagnostics then
      self.diagnostics[self.name] = DiagsApp{implementation = VlasovDiags()}
      self.diagnostics[self.name]:fullInit(self, field, self)
   end

   for srcNm, src in lume.orderedIter(self.sources) do
      self.diagnostics[self.name..srcNm] = src:createDiagnostics(self, field)
   end

   for bcNm, bc in lume.orderedIter(self.nonPeriodicBCs) do
      self.diagnostics[self.name..bcNm] = bc:createDiagnostics(self, field)
   end

   for collNm, coll in lume.orderedIter(self.collisions) do
      self.diagnostics[self.name..collNm] = coll:createDiagnostics(self, field)
   end
   lume.setOrder(self.diagnostics)
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
      if self.needGkMom then
	 -- When 3V Vlasov neutrals are coupled to GK species, calculate uPar for neutrals and
	 -- vtSq to use on GK grid for energy conservation.
	 self.confDotProduct:advance(tCurr, {self.uSelf,self.bHat}, {self.uPar})
	 self.confMult:advance(tCurr, {uPar, uPar}, {self.uParSq})
	 self.confDiv:advance(tCurr, {self.numDensity, self.ptclEnergy}, {self.ptclEnergyNorm})
	 -- Barrier over shared communicator before combine
         Mpi.Barrier(self.grid:commSet().sharedComm)
	 self.vtSqGk:combine(1/self.vdim, self.ptclEnergyNorm, -1/self.vdim, self.uParSq)
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
      local neutVtSq = neuts:selfPrimitiveMoments()[2]
      
      if tCurr == 0.0 then
	 species[self.name].collisions[self.collNmIoniz].collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      end
      species[self.name].collisions[self.collNmIoniz].collisionSlvr:advance(tCurr, {neutM0, neutVtSq, self.vtSqSelf}, {species[self.name].collisions[self.collNmIoniz].reactRate})
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
      
      species[self.neutNmCX].collisions[self.collNmCX].collisionSlvr:advance(tCurr, {m0, self.uSelf, neutU, self.vtSqSelf, neutVtSq}, {species[self.name].collisions[self.collNmCX].reactRate})
   end

   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:calcCouplingMoments(tCurr, rkIdx, species) end

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
   return tm
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

function VlasovSpecies:projToSource(proj)
   -- For backwards compatibility: in the case where the source is specified
   -- as a projection object in the input file, this function turns that
   -- projection object into a Source object.
   local tbl = proj.tbl
   local pow = tbl.power
   return Source { profile = proj, power = pow }
end

return VlasovSpecies
