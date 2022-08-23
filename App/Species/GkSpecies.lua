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
local VlasovEq       = require "Eq.Vlasov"
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

function GkSpecies:allocMomCouplingFields()
   assert(false, "GkSpecies:allocMomCouplingFields should not be called. Field object should allocate its own coupling fields")
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
   self.equation = GyrokineticEq.GkEq {
      onGrid     = self.grid,               mass         = self.mass,
      confGrid   = self.confGrid,           hasPhi       = hasPhi,
      phaseBasis = self.basis,              hasApar      = hasApar,
      confBasis  = self.confBasis,          hasSheathBCs = self.hasSheathBCs,
      confRange  = self.bmag:localRange(),  positivity   = self.positivity,
      charge     = self.charge,
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
      equation           = self.equation,
      zeroFluxDirections = self.zeroFluxDirections,
      updateDirections   = upd,
      clearOut           = false,   -- Continue accumulating into output field.
      globalUpwind       = not (self.basis:polyOrder()==1),   -- Don't reduce max speed.
   }
   
   if hasApar then
      -- Set up solver that adds on volume term involving dApar/dt and the entire vpar surface term.
      self.equationStep2 = GyrokineticEq.GkEqStep2 {
         onGrid     = self.grid,
         phaseBasis = self.basis,
         confBasis  = self.confBasis,
         charge     = self.charge,
         mass       = self.mass,
         positivity = self.positivity,
      }

      if self.basis:polyOrder()==1 then 
         -- This solver calculates vpar surface terms for Ohm's law. p=1 only!
         self.solverStep2 = Updater.HyperDisCont {
            onGrid             = self.grid,
            basis              = self.basis,
            cfl                = self.cfl,
            equation           = self.equation,
            zeroFluxDirections = self.zeroFluxDirections,
            updateDirections   = {self.cdim+1},    -- Only vpar terms.
            updateVolumeTerm   = false,            -- No volume term.
            clearOut           = false,            -- Continue accumulating into output field.
            globalUpwind       = not (self.basis:polyOrder()==1),   -- Don't reduce max speed.
         }
         -- Note that the surface update for this term only involves the vpar direction.
         self.solverStep3 = Updater.HyperDisCont {
            onGrid             = self.grid,
            basis              = self.basis,
            cfl                = self.cfl,
            equation           = self.equationStep2,
            zeroFluxDirections = self.zeroFluxDirections,
            updateDirections   = {self.cdim+1}, -- Only vpar terms.
            clearOut           = false,   -- Continue accumulating into output field.
            globalUpwind       = not (self.basis:polyOrder()==1),   -- Don't reduce max speed.
         }
      else
         -- Note that the surface update for this term only involves the vpar direction.
         self.solverStep2 = Updater.HyperDisCont {
            onGrid             = self.grid,
            basis              = self.basis,
            cfl                = self.cfl,
            equation           = self.equationStep2,
            zeroFluxDirections = self.zeroFluxDirections,
            updateDirections   = {self.cdim+1},
            clearOut           = false,   -- Continue accumulating into output field.
            globalUpwind       = false,   -- Don't reduce max speed.
         }
      end
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
   self.threeMomentsCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,
      phaseBasis = self.basis,
      confBasis  = self.confBasis,
      moment     = "GkThreeMoments",
      gkfacs     = {self.mass, self.bmag},
   }
   self.calcMaxwell = Updater.MaxwellianOnBasis {
      onGrid      = self.grid,
      phaseBasis  = self.basis,
      confGrid    = self.confGrid,
      confBasis   = self.confBasis,
      mass        = self.mass,
   }
   if self.needSelfPrimMom then
      -- This is used in calcCouplingMoments to reduce overhead and multiplications.
      -- If collisions are LBO, the following also computes boundary corrections and, if polyOrder=1, star moments.
--      self.threeMomentsLBOCalc = Updater.DistFuncMomentCalc {
--         onGrid     = self.grid,
--         phaseBasis = self.basis,
--         confBasis  = self.confBasis,
--         moment     = "GkThreeMomentsLBO",
--         gkfacs     = {self.mass, self.bmag},
--         positivity = self.positivity,
--      }
      local vbounds = ffi.new("double[4]")
      for i=0, self.vdim-1 do 
         vbounds[i]           = self.grid:lower(self.cdim+i+1)
         vbounds[i+self.vdim] = self.grid:upper(self.cdim+i+1)
      end
      self.primMomSelf = Updater.SelfPrimMoments {
         onGrid     = self.grid,
         phaseBasis = self.basis,
         confBasis  = self.confBasis,
         operator   = "GkLBO",
         mass       = self.mass,
         vbounds    = vbounds,
      }
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

   self._firstMomentCalc = true  -- To avoid re-calculating moments when not evolving.

   self.timers = {couplingMom = 0., sources = 0.}
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
   			self.neutNmIz = species[sN].collisions[collNm].neutNm
   			self.needSelfPrimMom  = true
			self.calcReactRate    = true
   			self.collNmIoniz      = collNm
			species[self.neutNmIz].calcIntSrcIz = true
			species[self.neutNmIz].collNmIoniz = collNm
   			counterIz_elc = false
		     elseif self.name==species[sN].collisions[collNm].neutNm and counterIz_neut then
			self.needSelfPrimMom = true
   		     end
   		  end
   	       end
   	    end
   	 end
      end
   end

   -- If Charge Exchange collision object exists, locate ions
   local counterCX_ion = true
   for sN, _ in lume.orderedIter(species) do
      if species[sN].collisions and next(species[sN].collisions) then 
         for sO, _ in lume.orderedIter(species) do
   	    if self.collPairs[sN][sO].on then
   	       if (self.collPairs[sN][sO].kind == 'CX') then
   		  for collNm, _ in pairs(species[sN].collisions) do
   		     if self.name==species[sN].collisions[collNm].ionNm and counterCX_ion then
   			self.calcCXSrc        = true			
   			self.collNmCX         = collNm
   			self.neutNmCX         = species[sN].collisions[collNm].neutNm
   			self.needSelfPrimMom  = true
			species[self.neutNmCX].needSelfPrimMom = true
			--self.vSigmaCX         = self:allocMoment()
   			species[self.neutNmCX].needSelfPrimMom = true
   			counterCX_ion = false
    		     end
   		  end
   	       end
   	    end
   	 end
      end
   end
   
   if self.needSelfPrimMom then
      -- Allocate fields to store self-species primitive moments.
      self.uParSelf = self:allocMoment()
      self.vtSqSelf = self:allocMoment()

      -- Allocate fields for boundary corrections.
      self.threeMomentsBoundaryCorrections = self:allocVectorMoment(3)

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
	       if species[sN].charge == 0 or species[sO].charge == 0 then
		  -- do nothing
               elseif self.nuVarXCross[otherNm] == nil then
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

function GkSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
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
   for nm, c in pairs(self.collisions) do
      c:advance(tCurr, fIn, species, {fRhsOut, self.cflRateByCell})
   end

   if self.evolveCollisionless then
      self.solver:advance(tCurr, {fIn, em, emFunc, dApardtProv}, {fRhsOut, self.cflRateByCell})
   else
      self.equation:setAuxFields({em, emFunc, dApardtProv})  -- Set auxFields in case they are needed by BCs/collisions.
   end

   for _, bc in pairs(self.nonPeriodicBCs) do
      bc:storeBoundaryFlux(tCurr, outIdx, fRhsOut)   -- Save boundary fluxes.
   end
   emIn[1]:useBoundaryFlux(tCurr, outIdx)  -- Some field objects need to use the boundary fluxes right away.

   for _, src in lume.orderedIter(self.sources) do src:advance(tCurr, fIn, species, fRhsOut) end
end

function GkSpecies:advanceStep2(tCurr, species, emIn, inIdx, outIdx)
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

function GkSpecies:advanceStep3(tCurr, species, emIn, inIdx, outIdx)
   local fIn = self:rkStepperFields()[inIdx]
   if self.positivityRescale then
      fIn = self.fPos
   end
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em = emIn[1]:rkStepperFields()[inIdx]
   local dApardtProv = emIn[1].dApardtProv
   local emFunc = emIn[2]:rkStepperFields()[1]

   if self.evolveCollisionless then
      self.solverStep3:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solverStep3:advance(tCurr, {fIn, em, emFunc, dApardtProv}, {fRhsOut})
   end
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
   if self.evolve or self._firstMomentCalc then
      local fIn     = self:rkStepperFields()[rkIdx]
      local tmStart = Time.clock()

      fIn = self.getF_or_deltaF(fIn)  -- Return full-F, or compute and return fluctuations.

      if self.needSelfPrimMom and
         lume.any({unpack(self.momentFlags,2,4)},function(x) return x==false end) then -- No need to recompute if already computed.

	 self.threeMomentsCalc:advance(tCurr, {fIn}, {self.threeMoments})

	 self.primMomSelf:advance(tCurr, {self.threeMoments, fIn, self.threeMomentsBoundaryCorrections}, {self.uParSelf, self.vtSqSelf})

         -- Indicate that moments, boundary corrections, star moments
         -- and self-primitive moments have been computed.
         for iF=2,4 do self.momentFlags[iF] = true end
      end
      if self.momentFlags[1]==false then -- No need to recompute if already computed.
         self.numDensityCalc:advance(tCurr, {fIn}, { self.numDensity })
         -- Indicate that first moment has been computed.
         self.momentFlags[1] = true
      end

      -- For ionization.
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

	 species[self.neutNmCX].collisions[self.collNmCX].collisionSlvr:advance(tCurr, {m0, self.uParSelf, neutU, self.vtSqSelf, neutVtSq}, {species[self.name].collisions[self.collNmCX].reactRate})
      end
      
      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
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

      fIn = self.getF_or_deltaF(fIn)
      self.numDensityCalc:advance(nil, {fIn}, { self.numDensityAux })

      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
   return self.numDensityAux
end

function GkSpecies:getMomDensity(rkIdx)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then return self.momDensity end 
   local fIn = self:rkStepperFields()[rkIdx]
 
   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()

      fIn = self.getF_or_deltaF(fIn)
      self.momDensityCalc:advance(nil, {fIn}, { self.momDensityAux })

      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
   return self.momDensityAux
end

-- Like getMomDensity, but use GkM1proj instead of GkM1, which uses cell-average v_parallel in moment calculation.
function GkSpecies:getMomProjDensity(rkIdx)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then return self.momDensity end 
   local fIn = self:rkStepperFields()[rkIdx]
 
   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()

      fIn = self.getF_or_deltaF(fIn)
      self.momProjDensityCalc:advance(nil, {fIn}, { self.momDensityAux })

      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
   return self.momDensityAux
end

function GkSpecies:getEmModifier(rkIdx)
   -- For p > 1, this is just numDensity.
   if self.basis:polyOrder() > 1 then return self:getNumDensity(rkIdx) end

   local fIn = self.equation.emMod

   if self.evolve or self._firstMomentCalc then
      local tmStart = Time.clock()

      fIn = self.getF_or_deltaF(fIn)
      self.momProjDensityCalc:advance(nil, {fIn}, { self.momDensityAux })

      self.timers.couplingMom = self.timers.couplingMom + Time.clock() - tmStart
   end
   if not self.evolve then self._firstMomentCalc = false end
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

function GkSpecies:projToSource(proj)
   local tbl = proj.tbl
   local pow = tbl.power
   return Source { profile = proj, power = pow }
end

return GkSpecies
