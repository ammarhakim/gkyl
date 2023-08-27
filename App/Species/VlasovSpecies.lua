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
   self.numDensity  = self:allocMoment()
   self.momDensity  = self:allocVectorMoment(self.vdim)
   self.fiveMoments = self:allocVectorMoment(self.vdim+2)

   self.ptclEnergyAux = self:allocMoment()

   -- Allocate field for external forces if any.
   if self.hasExtForce then 
      self.vExtForce = self:allocVectorMoment(self.vdim)
      self.vExtFptr, self.vExtFidxr = self.vExtForce:get(1), self.vExtForce:genIndexer()
   end
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
      hasE, hasB = true, true
   end

   if hasB then
      self.totalEmField = self:allocVectorMoment(8)     -- 8 components of EM field.
   else
      --self.totalEmField = self:allocVectorMoment(3)     -- Electric field only.
      self.totalEmField = self:allocMoment()  -- Phi only (Vlasov-Poisson)
   end
   if self.hasExtForce then
      self.totEmFptr, self.totEmFidxr = self.totalEmField:get(1), self.totalEmField:genIndexer()
   end

   self.computePlasmaB = true and plasmaB or extHasB   -- Differentiate plasma B from external B.

   -- Create table of pointers to fields needed in update
   self.fldPtrs = {self.totalEmField}
   -- Create updater to advance solution by one time-step.
   if self.evolveCollisionless then
      self.solver = Updater.VlasovDG {
         onGrid     = self.grid,                       hasElectricField = hasE,
         confBasis  = self.confBasis,                  hasMagneticField = hasB,
         phaseBasis = self.basis,                      hasExtForce      = self.hasExtForce,
         confRange  = self.totalEmField:localRange(),  phaseRange       = self.distf[1]:localRange(),  
         plasmaMagField = self.computePlasmaB,         fldPtrs          = self.fldPtrs, 
      }
      self.collisionlessAdvance = function(tCurr, inFlds, outFlds)
         self.solver:advance(tCurr, inFlds, outFlds)
      end

      -- Boundary flux updater.
      self.boundaryFluxSlvr = Updater.BoundaryFluxCalc {
         onGrid = self.grid,  equation    = self.solver:getEquation(), 
         cdim   = self.cdim,  equation_id = "vlasov",
      }
      self.collisionlessBoundaryAdvance = function(tCurr, inFlds, outFlds)
         self.boundaryFluxSlvr:advance(tCurr, inFlds, outFlds)
      end
   else
      self.solver = {totalTime = 0.}
      self.collisionlessAdvance = function(tCurr, inFlds, outFlds) end
      self.collisionlessBoundaryAdvance = function(tCurr, inFlds, outFlds) end
   end

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
   self.calcMaxwell = Updater.MaxwellianOnBasis {
      onGrid     = self.grid,   confGrid  = self.confGrid,
      phaseBasis = self.basis,  confBasis = self.confBasis,
   }
   if self.needFiveMoments then
      -- Create updater to compute M0, M1i, M2 moments.
      self.fiveMomentsCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis = self.confBasis,
         phaseBasis = self.basis,  moment    = "FiveMoments",
      }
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

   if self.hasExtForce then
      self.evalVlasovExtForce = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = self.vlasovExtForceFunc,
         basis  = self.confBasis,  onGhosts = false
      }

      self.accumulateExtForce = function(tCurr, totalEmField)
         local vExtForce  = self.vExtForce
         local vItr, eItr = self.vExtFptr, self.totEmFptr
         self.evalVlasovExtForce:advance(tCurr, {}, {vExtForce})

         -- Analogous to the current, the external force only gets accumulated onto the electric field.
         for idx in totalEmField:localRangeIter() do
            vExtForce:fill(self.vExtFidxr(idx), vItr)
            totalEmField:fill(self.totEmFidxr(idx), eItr)
            for i = 1, vExtForce:numComponents() do eItr[i] = eItr[i]+vItr[i] end
         end
      end
   else
      self.accumulateExtForce = function(tCurr, totalEmField) end
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

function VlasovSpecies:initCrossSpeciesCoupling(population)
   -- This method establishes the interaction between different
   -- species that is not mediated by the field (solver), like
   -- collisions.

   local species = population:getSpecies()

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
                  self.collPairs[sN][sO].kind = species[sN].collisions[collNm].collKind
                  self.needFiveMoments = true  -- MF 2022/09/16: currently all collision models need M0, M1, M2.
               end
            end
         else
            -- This species does not collide with anyone.
            self.collPairs[sN][sO].on = false
         end
      end
   end

   local isThisSpeciesMine = population:isSpeciesMine(self.name)
   -- Allocate fiveMoments if we collide with other species not in this rank.
   if self.needFiveMoments and (not isThisSpeciesMine) then
      self.fiveMoments = self:allocVectorMoment(self.vdim+2)
   end

   -- Create list of ranks we need to send/recv local fiveMoments to/from.
   self.fiveMomentsXfer = {}
   self.fiveMomentsXfer.destRank, self.fiveMomentsXfer.srcRank  = {}, {}
   self.fiveMomentsXfer.sendReqStat, self.fiveMomentsXfer.recvReqStat = nil, nil
   for sO, info in pairs(self.collPairs[self.name]) do
      local sOrank = population:getSpeciesOwner(sO)
      local selfRank = population:getSpeciesOwner(self.name)
      if sO~=self.name and info.on then
         if isThisSpeciesMine then
            -- Only species owned by this rank send fiveMoments to other ranks.
            if #self.fiveMomentsXfer.destRank == 0 and (not population:isSpeciesMine(sO)) then
               table.insert(self.fiveMomentsXfer.destRank, sOrank)
               self.fiveMomentsXfer.sendReqStat = Mpi.RequestStatus()
            end
         else
            -- Only species not owned by this rank receive fiveMoments from other ranks.
            if #self.fiveMomentsXfer.srcRank == 0 and (not population:isSpeciesMine(self.name)) then
               table.insert(self.fiveMomentsXfer.srcRank, selfRank)
               self.fiveMomentsXfer.recvReqStat = Mpi.RequestStatus()
            end
         end
      end
   end

   -- Initialize the BC cross-coupling interactions.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:initCrossSpeciesCoupling(species) end

   -- Initialize the source cross-coupling interactions.
   for _, src in lume.orderedIter(self.sources) do src:initCrossSpeciesCoupling(population) end
end

function VlasovSpecies:setActiveRKidx(rkIdx) self.activeRKidx = rkIdx end

function VlasovSpecies:advance(tCurr, population, emIn, inIdx, outIdx)
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
   if emField then totalEmField:accumulate(qbym, emField) end
   if emExternalField then totalEmField:accumulate(qbym, emExternalField) end

   -- If external force present (gravity, body force, etc.) accumulate it to electric field.
   self.accumulateExtForce(tCurr, totalEmField)

   self.collisionlessAdvance(tCurr, {fIn, totalEmField, emField}, {fRhsOut, self.cflRateByCell})

   -- Perform the collision update.
   for _, c in pairs(self.collisions) do
      c:advance(tCurr, fIn, population, {fRhsOut, self.cflRateByCell})   -- 'population' needed for cross-species collisions.
   end

   self.collisionlessBoundaryAdvance(tCurr, {fIn}, {fRhsOut})

   for _, src in lume.orderedIter(self.sources) do src:advance(tCurr, fIn, population:getSpecies(), fRhsOut) end
   
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do
      bc:storeBoundaryFlux(tCurr, outIdx, fRhsOut)   -- Save boundary fluxes.
   end

end

function VlasovSpecies:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx)
   -- Perform some operations after the updates have been computed, but before
   -- the combine RK (in PlasmaOnCartGrid) is called.

   local species = population:getSpecies()

   -- Wait to finish sending fiveMoments if needed.
   population:speciesXferField_waitSend(self.fiveMomentsXfer)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx)
   end

   for _, bc in pairs(self.nonPeriodicBCs) do bc:advanceCrossSpeciesCoupling(tCurr, species, outIdx) end

   for _, src in lume.orderedIter(self.sources) do src:advanceCrossSpeciesCoupling(tCurr, species, outIdx) end
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
      -- Neutrals haven't been updated yet, so we need to compute their moments and primitive moments. 
      neuts:calcCouplingMoments(tCurr, rkIdx, species)
      local neutM0   = neuts:fluidMoments()[1]
      local neutVtSq = neuts:selfPrimitiveMoments()[2]
      
      species[self.name].collisions[self.collNmIoniz].collisionSlvr:advance(tCurr, {neutM0, neutVtSq, self.vtSqSelf}, {species[self.name].collisions[self.collNmIoniz].reactRate, self.cflRateByCell})
   end

   -- For charge exchange.
   if self.calcCXSrc then
      -- Calculate Vcx*SigmaCX.
      local neuts = species[self.neutNmCX]
      -- Neutrals haven't been updated yet, so we need to compute their moments and primitive moments. 
      neuts:calcCouplingMoments(tCurr, rkIdx, species)
      local m0       = neuts:fluidMoments()[1]
      local neutU    = neuts:selfPrimitiveMoments()[1]
      local neutVtSq = neuts:selfPrimitiveMoments()[2]
      
      species[self.neutNmCX].collisions[self.collNmCX].collisionSlvr:advance(tCurr, {m0, self.uSelf, neutU, self.vtSqSelf, neutVtSq}, {species[self.name].collisions[self.collNmCX].reactRate})
   end

   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:calcCouplingMoments(tCurr, rkIdx, species) end

   self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
end

function VlasovSpecies:calcCrossCouplingMoments(tCurr, rkIdx, population)
   -- Perform cross-species calculation related to coupling moments that require the
   -- self-species coupling moments.

   -- Begin sending/receiving fiveMoments if needed.
   population:speciesXferField_begin(self.fiveMomentsXfer, self.fiveMoments, 22)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:calcCrossCouplingMoments(tCurr, rkIdx, population)
   end
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
function VlasovSpecies:Maxwellian(xn, n0, vdnIn, T0)
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
