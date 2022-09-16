-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov species
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Basis          = require "Basis"
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

-- Function to create basis functions.
local function createBasis(nm, cdim, vdim, polyOrder)
   local ndim = cdim+vdim
   if nm == "serendipity" then
      if polyOrder == 1 then
         return Basis.CartModalHybrid { cdim = cdim, vdim = vdim }
      else
         return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
      end
   elseif nm == "tensor" then
      return Basis.CartModalTensor { ndim = ndim, polyOrder = polyOrder }
   end
end

function VlasovSpecies:createBasis(nm, polyOrder)
   self.basis = createBasis(nm, self.cdim, self.vdim, polyOrder)
   for _, c in pairs(self.collisions) do c:setPhaseBasis(self.basis) end

   -- Output of grid file is placed here because as the file name is associated
   -- with a species, we wish to save the basisID and polyOrder in it. But these
   -- can only be extracted from self.basis after this is created.
   if self.grid:getMappings() then
      local metaData = {polyOrder = self.basis:polyOrder(),
                        basisType = self.basis:id(),
                        charge    = self.charge,
                        mass      = self.mass,
                        grid      = GKYL_OUT_PREFIX .. "_" .. self.name .. "_grid.bp"}
      self.grid:write(self.name .. "_grid.bp", 0.0, metaData)
   end
end

function VlasovSpecies:alloc(nRkDup)
   -- Allocate distribution function.
   VlasovSpecies.super.alloc(self, nRkDup)

   -- Allocate fields to store coupling moments (for use in coupling
   -- to field and collisions).
   self.numDensity = self:allocMoment()
   self.momDensity = self:allocVectorMoment(self.vdim)
   self.ptclEnergy = self:allocMoment()
   self.fiveMoments = self:allocVectorMoment(self.vdim+2)

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
   self.hasExtForce = false

   if tbl.vlasovExtForceFunc then
      self.vlasovExtForceFunc = tbl.vlasovExtForceFunc
      self.hasExtForce = true
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
   -- Run the KineticSpecies 'createSolver()' to initialize the collisions solver.
   VlasovSpecies.super.createSolver(self, field, externalField)

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
      --self.totalEmField = self:allocVectorMoment(3)     -- Electric field only.
      self.totalEmField = self:allocMoment()  -- Phi only (Vlasov-Poisson)
   end

   self.computePlasmaB = true and plasmaB   -- Differentiate plasma B from external B.

   ---- Create updater to advance solution by one time-step.
   self.solver = Updater.VlasovDG {
      onGrid     = self.grid,                       hasElectricField = hasE,
      confBasis  = self.confBasis,                  hasMagneticField = hasB,
      phaseBasis = self.basis,                      hasExtForce      = self.hasExtForce,
      confRange  = self.totalEmField:localRange(),  plasmaMagField   = plasmaB
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
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  moment     = "FiveMoments",
   }
   self.calcMaxwell = Updater.MaxwellianOnBasis {
      onGrid     = self.grid,   confGrid  = self.confGrid,
      phaseBasis = self.basis,  confBasis = self.confBasis,
   }
   if self.needFiveMoments then
      self.calcSelfCouplingMom = self.computePlasmaB
         and function(tCurr, fIn)
            -- Compute M0, M1i and M2.
            self.fiveMomentsCalc:advance(tCurr, {fIn}, {self.fiveMoments})
            -- Copy momentum density to its own field.
            self.momDensity:combineOffset(1., self.fiveMoments, 1*self.confBasis:numBasis())
         end
         or function(tCurr, fIn)
            -- Compute M0, M1i and M2.
            self.fiveMomentsCalc:advance(tCurr, {fIn}, {self.fiveMoments})
            -- Copy number density to its own field.
            self.numDensity:combineOffset(1., self.fiveMoments, 0)
         end
   else
      self.calcSelfCouplingMom = self.computePlasmaB
         and function(tCurr, fIn)
            self.momDensityCalc:advance(tCurr, {fIn}, { self.momDensity })
         end
         or function(tCurr, fIn)
            self.numDensityCalc:advance(tCurr, {fIn}, { self.numDensity })
         end
   end

   if self.vlasovExtForceFunc then
      self.evalVlasovExtForce = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = self.vlasovExtForceFunc,
         basis  = self.confBasis,  onGhosts = false
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

   -- Determine if M0, M1i and M2 are needed.
   self.needFiveMoments = false
   -- Create a double nested table of colliding species.
   -- In this table we will encode information about that collition such as:
   --   * does the collision take place?
   --   * Operator modeling the collision.
   -- Other features of a collision may be added in the future.
   self.collPairs = {}
   for sN, _ in lume.orderedIter(species) do
      self.collPairs[sN] = {}
      for sO, _ in lume.orderedIter(species) do
         self.collPairs[sN][sO] = {}
         -- Need next below because species[].collisions is created as an empty table. 
         if species[sN].collisions and next(species[sN].collisions) then 
            for collNm, _ in pairs(species[sN].collisions) do
               -- This species collides with someone.
               self.collPairs[sN][sO].on = lume.any(species[sN].collisions[collNm].collidingSpecies,
                                                    function(e) return e==sO end)
               if self.collPairs[sN][sO].on then
                  self.collPairs[sN][sO].kind = species[sO].collisions[collNm].collKind
                  self.needFiveMoments = true  -- MF 2022/09/16: at the moment all collision models need M0, M1, M2.
               end
            end
         else
            -- This species does not collide with anyone.
            self.collPairs[sN][sO].on = false
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
   			self.calcReactRate    = true
   			self.collNmIoniz      = collNm
			counterIz_elc         = false
		     elseif self.name==species[sN].collisions[collNm].neutNm and counterIz_neut then
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
			species[self.neutNmCX].needFiveMoments = true
   			counterCX = false
    		     end
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
   fRhsOut:clear(0.0)

   -- Accumulate functional Maxwell fields (if needed).
   local emField         = emIn[1]:rkStepperFields()[inIdx]
   local emExternalField = emIn[2]:rkStepperFields()[1]
   local totalEmField    = self.totalEmField
   totalEmField:clear(0.0)

   local qbym = self.charge/self.mass
   if emField then 
      totalEmField:accumulate(qbym, emField) 
   end
   if emExternalField then totalEmField:accumulate(qbym, emExternalField) end

   -- If external force present (gravity, body force, etc.) accumulate it to electric field.
   if self.hasExtForce then
      local vExtForce = self.vExtForce
      self.evalVlasovExtForce:advance(tCurr, {}, {vExtForce})

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
      self.solver:advance(tCurr, {fIn, totalEmField, emField}, {fRhsOut, self.cflRateByCell})
   end

   -- Perform the collision update.
   for _, c in pairs(self.collisions) do
      c:advance(tCurr, fIn, species, {fRhsOut, self.cflRateByCell})   -- 'species' needed for cross-species collisions.
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
   -- Compute moments needed in coupling to fields and collisions.
   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]

   -- Compute M0, M1i and/or M2 depending on what fields and collisions need.
   self.calcSelfCouplingMom(tCurr, fIn)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:calcCouplingMoments(tCurr, rkIdx, species)
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
      
      species[self.name].collisions[self.collNmIoniz].collisionSlvr:advance(tCurr, {neutM0, neutVtSq, self.vtSqSelf}, {species[self.name].collisions[self.collNmIoniz].reactRate, self.cflRateByCell})
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

function VlasovSpecies:fluidMoments() return self.fiveMoments end

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
