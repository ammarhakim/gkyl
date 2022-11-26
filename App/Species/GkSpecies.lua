-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic species object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto          = require "Lib.Proto"
local KineticSpecies = require "App.Species.KineticSpecies"
local Basis          = require "Basis"
local Mpi            = require "Comm.Mpi"
local GyrokineticEq  = require "Eq.Gyrokinetic"
local Updater        = require "Updater"
local DataStruct     = require "DataStruct"
local Time           = require "Lib.Time"
local Constants      = require "Lib.Constants"
local Lin            = require "Lib.Linalg"
local xsys           = require "xsys"
local Source         = require "App.Sources.GkSource"
local DiagsApp       = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags        = require "App.Diagnostics.GkDiagnostics"
local BasicBC        = require("App.BCs.GkBasic").GkBasic
local BCsBase        = require "App.BCs.BCsBase"
local ffi            = require "ffi"
local lume           = require "Lib.lume"

local GkSpecies = Proto(KineticSpecies)

-- ............. Backwards compatible treatment of BCs .....................--
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
GkSpecies.bcSheath   = SP_BC_SHEATH      -- Sheath.
GkSpecies.bcZeroFlux = SP_BC_ZEROFLUX    -- Zero flux.
GkSpecies.bcCopy     = SP_BC_COPY        -- Copy stuff.

function GkSpecies:makeBcApp(bcIn, dir, edge)
   local bcOut
   if type(bcIn) == "function" then
      bcOut = BasicBC{kind="function", bcFunction=bcIn, diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_COPY then
      print("GkSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="copy", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_ABSORB then
      print("GkSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="absorb", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_OPEN then
      print("GkSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      -- AHH: open seems unstable. So using plain copy.
      bcOut = BasicBC{kind="copy", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_REFLECT then
      print("GkSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="reflect", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_SHEATH then
      print("GkSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="sheath", diagnostics={}, saveFlux=false}
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
         return Basis.CartModalGkHybrid { cdim = cdim, vdim = vdim }
      else
         return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
      end
   elseif nm == "tensor" then
      return Basis.CartModalTensor { ndim = ndim, polyOrder = polyOrder }
   end
end

function GkSpecies:createBasis(nm, polyOrder)
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

function GkSpecies:alloc(nRkDup)
   -- Allocate distribution function.
   GkSpecies.super.alloc(self, nRkDup)

   -- Allocate fields to store coupling moments (for use in coupling
   -- to field and collisions).
   self.numDensity    = self:allocMoment()
   self.numDensityAux = self:allocMoment()
   self.momDensity    = self:allocMoment()
   self.momDensityAux = self:allocMoment()
   self.ptclEnergy    = self:allocMoment()
   self.ptclEnergyAux = self:allocMoment()
   self.threeMoments  = self:allocVectorMoment(3)
   if self.positivity then
      self.numDensityPos = self:allocMoment()
      self.momDensityPos = self:allocMoment()
      self.ptclEnergyPos = self:allocMoment()
   end
			
   self.vDegFreedom = self.vdim == 1 and 1.0 or 3.0

   self.first = true
end

function GkSpecies:fullInit(appTbl)
   GkSpecies.super.fullInit(self, appTbl)
end

function GkSpecies:createSolver(field, externalField)
   -- Run the KineticSpecies 'createSolver()' to initialize the collisions solver.
   GkSpecies.super.createSolver(self, field, externalField)

   local hasE, hasB       = field:hasEB()
   local extHasE, extHasB = externalField:hasEB()

   local hasPhi  = hasE or extHasE
   local hasApar = hasB or extHasB

   -- Set up Jacobian.
   if externalField then
      self.bmagFunc = externalField.bmagFunc
      -- If vdim>1, get the phase-space Jacobian (=bmag) from geo.
      self.jacobPhaseFunc = self.bmagFunc
      self.jacobGeoFunc   = externalField.jacobGeoFunc

      self.bmag        = assert(externalField.geo.bmag, "nil bmag")
      self.bmagInv     = externalField.geo.bmagInv
      self.bmagInvSq   = externalField.geo.bmagInvSq
      self.jacobGeo    = externalField.geo.jacobGeo
      self.jacobGeoInv = externalField.geo.jacobGeoInv

      self.jacobGeoDbmagSq = self:allocMoment()
      self.confWeakMultiply:advance(0., {self.jacobGeo, self.bmagInvSq}, {self.jacobGeoDbmagSq})

      -- Compute the magnetic field in the center of the domain (e.g. for the Poisson equation).
      local xMid = {}
      for d = 1,self.cdim do xMid[d]=self.confGrid:mid(d) end
      self.bmagMid = self.bmagFunc(0.0, xMid)
   end

   self.hasSheathBCs = false
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do
      self.hasSheathBCs = self.hasSheathBCs or (bc.bcKind=="sheath" and true or false)
   end

   -- Create updater to advance solution by one time-step.
   if self.evolveCollisionless then
      self.solver = Updater.GyrokineticDG {
         onGrid     = self.grid,       confRange = self.bmag:localRange(), 
         confBasis  = self.confBasis,  charge    = self.charge,
         phaseBasis = self.basis,      mass      = self.mass,
      }
      self.collisionlessAdvance = function(tCurr, inFlds, outFlds)
         self.solver:advance(tCurr, inFlds, outFlds)
      end
      
--      if hasApar then
--         -- Set up solver that adds on volume term involving dApar/dt and the entire vpar surface term.
--         self.equationStep2 = GyrokineticEq.GkEqStep2 {
--            onGrid     = self.grid,       charge     = self.charge,
--            phaseBasis = self.basis,      mass       = self.mass,
--            confBasis  = self.confBasis,  positivity = self.positivity,
--         }
--
--         -- Note that the surface update for this term only involves the vpar direction.
--         self.solverStep2 = Updater.HyperDisCont {
--            onGrid   = self.grid,           zeroFluxDirections = self.zeroFluxDirections,
--            basis    = self.basis,          updateDirections   = {self.cdim+1},
--            cfl      = self.cfl,            clearOut           = false,   -- Continue accumulating into output field.
--            equation = self.equationStep2,  globalUpwind       = false,   -- Don't reduce max speed.
--         }
--      end
   else
      self.solver = {totalTime = 0.}
      self.collisionlessAdvance = function(tCurr, inFlds, outFlds) end
   end
   
   -- Create updaters to compute various moments.
   self.numDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM0", -- GkM0 = < f >
   }
   self.momDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM1", -- GkM1 = < v_parallel f > 
   }
   self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM2", -- GkM2 = < (v_parallel^2 + 2*mu*B/m) f >
   }
   self.M2parCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM2par", -- GkM2par = < v_parallel^2 f >
   }
   self.M3parCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM3par", -- GkM3par = < v_parallel^3 f >
   }
   if self.vdim > 1 then
      self.M2perpCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis  = self.confBasis,
         phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
         moment     = "GkM2perp", -- GkM2 = < (mu*B/m) f >
      }
      self.M3perpCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis  = self.confBasis,
         phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
         moment     = "GkM3perp", -- GkM3perp = < vpar*(mu*B/m) f >
      }
   end
   self.calcMaxwell = Updater.MaxwellianOnBasis {
      onGrid     = self.grid,   confGrid  = self.confGrid,
      phaseBasis = self.basis,  confBasis = self.confBasis,
      mass       = self.mass,
   }
   if self.needThreeMoments then
      -- Create updater to compute M0, M1, M2 moments.
      self.threeMomentsCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis = self.confBasis,
         phaseBasis = self.basis,  moment    = "GkThreeMoments",
         gkfacs     = {self.mass, self.bmag},
      }
      self.calcSelfCouplingMom = function(tCurr, fIn)
         -- Compute M0, M1i and M2.
         self.threeMomentsCalc:advance(tCurr, {fIn}, {self.threeMoments})
         -- Copy number density to its own field.
         self.numDensity:combineOffset(1., self.threeMoments, 0)
      end
   else
      self.calcSelfCouplingMom = function(tCurr, fIn)
         self.numDensityCalc:advance(tCurr, {fIn}, { self.numDensity })
      end
   end

   -- Create an updater for volume integrals. Used by diagnostics.
   -- Placed in a table with key 'scalar' to keep consistency with VlasovSpecies (makes diagnostics simpler).
   self.volIntegral = {
      scalar = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.confGrid,   numComponents = 1,
         basis  = self.confBasis,  quantity      = "V",
      }
   }

   -- Select the function that returns the mass density factor for the polarization (Poisson equation).
   -- Allow user to specify polarization density factor (n in polarization density).
   -- If the polarization is linearized, it should be specified in the input file (as a number, file, or function).
   self.polWeight = self:allocMoment() -- Polarization weight mass*jacobGeo*n0/B^2 (for Poisson equation).
   if field.linearizedPolarization then
      local evOnNodes = Updater.EvalOnNodes {
         onGrid = self.confGrid,   evaluate = function(t, xn) return 1. end,
         basis  = self.confBasis,  onGhosts = false, --true,
      }
      if self.tbl.polarizationDensityFactor == nil then
         print("*** App.Species.GkSpecies: WARNING... not specifying 'polarizationDensityFactor' and relying on n0 in the input file will be deprecated. Please change your input file to specify 'polarizationDensityFactor' (the density factor in the Poisson equation). ***")
         local den0 = self.tbl.n0 or n0
         evOnNodes:setFunc(function(t,xn) return den0*self.mass/(self.bmagMid^2) end)
         evOnNodes:advance(0., {}, {self.polWeight})
      else
         local polDenFacIn = self.tbl.polarizationDensityFactor
   
         if type(polDenFacIn) == "number" then
            evOnNodes:setFunc(function(t,xn) return polDenFacIn*self.mass/(self.bmagMid^2) end)
            evOnNodes:advance(0., {}, {self.polWeight})
         elseif type(polDenFacIn) == "string" or type(polDenFacIn) == "function" then
            if type(polDenFacIn) == "string" then
               self.distIo:read(self.polWeight, polDenFacIn) --, true)
            else
               evOnNodes:setFunc(polDenFacIn)
               evOnNodes:advance(0., {}, {self.polWeight})
            end
   
            -- Apply open BCs (although BCs here should not matter/be used).
            local function makeOpenBcUpdater(dir, edge)
               local bcOpen = function(dir, tm, idxIn, fIn, fOut)
                  self.confBasis:flipSign(dir, fIn, fOut)   -- Requires skinLoop = "pointwise".
               end
               return Updater.Bc {
                  onGrid   = self.confGrid,  edge = edge,
                  dir      = dir,            boundaryConditions = {bcOpen},
                  skinLoop = "pointwise",
               }
            end
            local openBCupdaters = {}
            for dir = 1, self.cdim do
               if not lume.any(self.confGrid:getPeriodicDirs(), function(t) return t==dir end) then
                  openBCupdaters["lower"] = makeOpenBcUpdater(dir, "lower")
                  openBCupdaters["upper"] = makeOpenBcUpdater(dir, "upper")
               end
            end
            for _, bc in pairs(openBCupdaters) do bc:advance(0., {}, {self.polWeight}) end
            self.polWeight:sync(true)
   
            self.distIo:write(self.polWeight, string.format("%s_polarizationDensityFactor_%d.bp", self.name, self.diagIoFrame), 0., self.diagIoFrame, false) --true)
   
            self.polWeight:scale(self.mass)
            self.confWeakMultiply:advance(0., {self.bmagInvSq, self.polWeight}, {self.polWeight})
         end
      end

      -- Include a factor of jacobGeo (included naturally when linearizedPolarization=false).
      self.confWeakMultiply:advance(0., {self.jacobGeo, self.polWeight}, {self.polWeight})

      self.getPolWeight = function() return self.polWeight end
   else
      self.getPolWeight = function()
         self.polWeight:combine(self.mass, self.numDensity)
         self.confWeakMultiply:advance(0., {self.jacobGeoDbmagSq, self.polWeight}, {self.polWeight})
         return self.polWeight
      end
   end

   -- Create species source solvers.
   for _, src in lume.orderedIter(self.sources) do src:createSolver(self, externalField) end

   self.timers = {couplingMom = 0., sources = 0.}
end

function GkSpecies:initCrossSpeciesCoupling(population)
   -- This method establishes the interaction between different
   -- species that is not mediated by the field (solver), like
   -- collisions.

   local species = population:getSpecies()

   -- Determine if M0, M1i and M2 are needed.
   self.needThreeMoments = false

   -- Create a double nested table of colliding species.
   -- In this table we will encode information about that collition such as:
   --   * does the collision take place?
   --   * Operator modeling the collision.
   -- Other features of a collision may be added in the future
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
                  self.needThreeMoments = true  -- MF 2022/09/16: currently all collision models need M0, M1, M2.
               end
            end
         else
            -- This species does not collide with anyone.
            self.collPairs[sN][sO].on = false
         end
      end
   end

   local isThisSpeciesMine = population:isSpeciesMine(self.name)
   -- Allocate threeMoments if we collide with other species not in this rank.
   if self.needThreeMoments and (not isThisSpeciesMine) then
      self.threeMoments = self:allocVectorMoment(3)
   end

   -- Create list of ranks we need to send/recv local threeMoments to/from.
   self.threeMomentsXfer = {}
   self.threeMomentsXfer.destRank, self.threeMomentsXfer.srcRank  = {}, {}
   self.threeMomentsXfer.sendReqStat, self.threeMomentsXfer.recvReqStat = nil, nil
   for sO, info in pairs(self.collPairs[self.name]) do
      local sOrank = population:getSpeciesOwner(sO)
      local selfRank = population:getSpeciesOwner(self.name)
      if sO~=self.name and info.on then
         if isThisSpeciesMine then
            -- Only species owned by this rank send threeMoments to other ranks.
            if #self.threeMomentsXfer.destRank == 0 and (not population:isSpeciesMine(sO)) then
               table.insert(self.threeMomentsXfer.destRank, sOrank)
               self.threeMomentsXfer.sendReqStat = Mpi.RequestStatus()
            end
         else
            -- Only species not owned by this rank receive threeMoments from other ranks.
            if #self.threeMomentsXfer.srcRank == 0 and (not population:isSpeciesMine(self.name)) then
               table.insert(self.threeMomentsXfer.srcRank, selfRank)
               self.threeMomentsXfer.recvReqStat = Mpi.RequestStatus()
            end
         end
       end
    end

   -- Initialize the BC cross-coupling interactions.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:initCrossSpeciesCoupling(species) end

end

function GkSpecies:advance(tCurr, population, emIn, inIdx, outIdx)
   self:setActiveRKidx(inIdx)
   self.tCurr = tCurr
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em     = emIn[1]:rkStepperFields()[inIdx] -- Dynamic fields (e.g. phi, Apar)
   local emFunc = emIn[2]:rkStepperFields()[1]     -- Geometry/external field.

   local dApardtProv = emIn[1].dApardtProv

   -- Rescale slopes.
   if self.positivityRescale then
      self.posRescaler:advance(tCurr, {fIn}, {self.fPos}, false)
      fIn = self.fPos
   end

   -- Clear RHS, because HyperDisCont set up with clearOut = false.
   fRhsOut:clear(0.0)

   -- Do collisions first so that collisions contribution to cflRate is included in GK positivity.
   for _, c in pairs(self.collisions) do
      c:advance(tCurr, fIn, population, {fRhsOut, self.cflRateByCell})
   end

   self.collisionlessAdvance(tCurr, {fIn, em, emFunc, dApardtProv}, {fRhsOut, self.cflRateByCell})

   for _, bc in pairs(self.nonPeriodicBCs) do
      bc:storeBoundaryFlux(tCurr, outIdx, fRhsOut)   -- Save boundary fluxes.
   end
   emIn[1]:useBoundaryFlux(tCurr, outIdx)  -- Some field objects need to use the boundary fluxes right away.

   for _, src in lume.orderedIter(self.sources) do src:advance(tCurr, fIn, population:getSpecies(), fRhsOut) end
end

function GkSpecies:advanceStep2(tCurr, population, emIn, inIdx, outIdx)
   local fIn = self:rkStepperFields()[inIdx]
   if self.positivityRescale then
      fIn = self.fPos
   end
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em = emIn[1]:rkStepperFields()[inIdx]
   local dApardtProv = emIn[1].dApardtProv
   local emFunc = emIn[2]:rkStepperFields()[1]

   if self.evolveCollisionless then
      self.solverStep2:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solverStep2:advance(tCurr, {fIn, em, emFunc, dApardtProv}, {fRhsOut})
   end
end

function GkSpecies:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx)
   -- Perform some operations after the updates have been computed, but before
   -- the combine RK (in PlasmaOnCartGrid) is called.

   local species = population:getSpecies()

   -- Wait to finish sending threeMoments if needed.
   population:speciesXferField_waitSend(self.threeMomentsXfer)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx)
   end

   for _, bc in pairs(self.nonPeriodicBCs) do bc:advanceCrossSpeciesCoupling(tCurr, species, outIdx) end
end

function GkSpecies:createDiagnostics(field)
   -- Run the KineticSpecies 'createDiagnostics()' (e.g. to create divideByJacobGeo()).
   GkSpecies.super.createDiagnostics(self, field)

   -- Create this species' diagnostics.
   if self.tbl.diagnostics then
      self.diagnostics[self.name] = DiagsApp{implementation = GkDiags()}
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

function GkSpecies:calcCouplingMoments(tCurr, rkIdx, species)
   -- Compute moments needed in coupling to fields and collisions.
   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]

   fIn = self.getF_or_deltaF(fIn)  -- Return full-F, or compute and return fluctuations.

   -- Compute M0, M1i and/or M2 depending on what fields and collisions need.
   self.calcSelfCouplingMom(tCurr, fIn)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:calcCouplingMoments(tCurr, rkIdx, species)
   end

   -- For ionization.
   if self.calcReactRate then
      local neuts = species[self.neutNmIz]
      -- Neutrals haven't been updated yet, so we need to compute their moments and primitive moments.
      neuts:calcCouplingMoments(tCurr, rkIdx, species)
      local neutM0   = neuts:fluidMoments()[1]
      local neutVtSq = neuts:selfPrimitiveMoments()[2]
         
      if tCurr == 0.0 then
         species[self.name].collisions[self.collNmIoniz].collisionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      end
     
      species[self.name].collisions[self.collNmIoniz].collisionSlvr:advance(tCurr, {neutM0, neutVtSq, self.vtSqSelf}, {species[self.name].collisions[self.collNmIoniz].reactRate})
   end

   if self.calcCXSrc then
      -- Calculate Vcx*SigmaCX.
      local neuts = species[self.neutNmCX]
      -- Neutrals haven't been updated yet, so we need to compute their moments and primitive moments.
      neuts:calcCouplingMoments(tCurr, rkIdx, species)
      local m0       = neuts:fluidMoments()[1]
      local neutU    = neuts:selfPrimitiveMoments()[1]
      local neutVtSq = neuts:selfPrimitiveMoments()[2]

      species[self.neutNmCX].collisions[self.collNmCX].collisionSlvr:advance(tCurr, {m0, self.uParSelf, neutU, self.vtSqSelf, neutVtSq}, {species[self.name].collisions[self.collNmCX].reactRate})
   end
   
   self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
end

function GkSpecies:calcCrossCouplingMoments(tCurr, rkIdx, population)
   -- Perform cross-species calculation related to coupling moments that require the
   -- self-species coupling moments.

   -- Begin sending/receiving threeMoments if needed.
   population:speciesXferField_begin(self.threeMomentsXfer, self.threeMoments, 22)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:calcCrossCouplingMoments(tCurr, rkIdx, population)
   end
end

function GkSpecies:fluidMoments() return self.threeMoments end

function GkSpecies:getNumDensity(rkIdx)
   -- If no rkIdx specified, assume numDensity has already been calculated.
   if rkIdx == nil then return self.numDensity end 
   local fIn = self:rkStepperFields()[rkIdx]

   local tmStart = Time.clock()

   fIn = self.getF_or_deltaF(fIn)
   self.numDensityCalc:advance(nil, {fIn}, { self.numDensityAux })

   self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   return self.numDensityAux
end

function GkSpecies:getMomDensity(rkIdx)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then return self.momDensity end 
   local fIn = self:rkStepperFields()[rkIdx]
 
   local tmStart = Time.clock()

   fIn = self.getF_or_deltaF(fIn)
   self.momDensityCalc:advance(nil, {fIn}, { self.momDensityAux })

   self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   return self.momDensityAux
end

-- Like getMomDensity, but use GkM1proj instead of GkM1, which uses cell-average v_parallel in moment calculation.
function GkSpecies:getMomProjDensity(rkIdx)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then return self.momDensity end 
   local fIn = self:rkStepperFields()[rkIdx]
 
   if self.evolve then
      local tmStart = Time.clock()

      fIn = self.getF_or_deltaF(fIn)
      self.momProjDensityCalc:advance(nil, {fIn}, { self.momDensityAux })

      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   end
   return self.momDensityAux
end

function GkSpecies:getEmModifier(rkIdx)
   -- For p > 1, this is just numDensity.
   if self.basis:polyOrder() > 1 then return self:getNumDensity(rkIdx) end

   local fIn = self.equation.emMod

   if self.evolve then
      local tmStart = Time.clock()

      fIn = self.getF_or_deltaF(fIn)
      self.momProjDensityCalc:advance(nil, {fIn}, { self.momDensityAux })

      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   end
   fIn:clear(0.0)
   return self.momDensityAux
end

function GkSpecies:getPolarizationWeight() return self.getPolWeight() end

function GkSpecies:getSrcCX() return self.srcCX end

function GkSpecies:getVSigmaCX() return self.vSigmaCX end

function GkSpecies:momCalcTime()
   local tm = self.timers.couplingMom
   return tm
end

function GkSpecies:solverVolTime() return self.equation.totalVolTime end

function GkSpecies:solverSurfTime() return self.equation.totalSurfTime end

function GkSpecies:totalSolverTime()
   local timer = self.solver and self.solver.totalTime or 0.
   if self.solverStep2 then timer = timer + self.solverStep2.totalTime end
   if self.posRescaler then timer = timer + self.posRescaler.totalTime end
   return timer
end

function GkSpecies:Maxwellian(xn, n0, vdIn, T0)
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

function GkSpecies:projToSource(proj)
   local tbl = proj.tbl
   local pow = tbl.power
   return Source { profile = proj, power = pow }
end

return GkSpecies
