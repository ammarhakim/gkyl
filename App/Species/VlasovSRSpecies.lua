-- Gkyl ------------------------------------------------------------------------
--
-- Special Relativistic Vlasov species
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

local VlasovSRSpecies = Proto(KineticSpecies)

-- Function to create phase space basis functions.
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

-- Function to create phase space basis functions.
local function createVelBasis(nm, vdim, polyOrder)
   if nm == "serendipity" then
      if polyOrder == 1 then
         return Basis.CartModalSerendipity { ndim = vdim, polyOrder = 2 }
      else
         return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
      end
   elseif nm == "tensor" then
      return Basis.CartModalTensor { ndim = ndim, polyOrder = polyOrder }
   end
end

function VlasovSRSpecies:createBasis(nm, polyOrder)
   self.basis = createBasis(nm, self.cdim, self.vdim, polyOrder)
   self.velBasis = createBasis(nm, self.vdim, polyOrder)
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

function VlasovSRSpecies:alloc(nRkDup)
   -- Allocate distribution function.
   VlasovSRSpecies.super.alloc(self, nRkDup)

   -- Allocate fields to store coupling moments (for use in coupling
   -- to field and collisions).
   self.numDensity = self:allocMoment()
   self.momDensity = self:allocVectorMoment(self.vdim)
   -- Four-momentum density (GammaV*n, GammaV*n*Vx, GammaV*n*Vy, GammaV*n*Vz)
   self.Ni  = self:allocVectorMoment(self.vdim+1)

   -- Allocate fields to store velocity space arrays:
   -- p/gamma (relativistic velocity)
   -- gamma (particle Lorentz boost factor = sqrt(1 + p^2))
   -- gamma_inv = 1/gamma
   self.p_over_gamma  = self:allocVectorVelMoment(self.vdim)
   self.gamma  = self:allocVelMoment()
   self.gamma_inv  = self:allocVelMoment()

   -- Allocate additional auxiliary fields needed for certain relativistic moments
   -- V_drift (bulk velocity)
   -- GammaV2 (bulk velocity Lorentz boost factor squared = 1/(1 - V_drift^2))
   -- GammaV_inv = 1/GammaV = sqrt(1 - V_drift^2)
   self.V_drift = self:allocVectorMoment(self.vdim)
   self.GammaV2 = self:allocMoment()
   self.GammaV_inv = self:allocMoment()

   self.ptclEnergyAux = self:allocMoment()

   -- Allocate field for external forces if any.
   if self.hasExtForce then 
      self.vExtForce = self:allocVectorMoment(self.vdim)
      self.vExtFptr, self.vExtFidxr = self.vExtForce:get(1), self.vExtForce:genIndexer()
   end
end

function VlasovSRSpecies:fullInit(appTbl)
   VlasovSRSpecies.super.fullInit(self, appTbl)

   local tbl = self.tbl

   -- If there is an external force, get the force function.
   self.hasExtForce = false
   if tbl.vlasovExtForceFunc then
      self.vlasovExtForceFunc = tbl.vlasovExtForceFunc
      self.hasExtForce = true
   end
end

function VlasovSRSpecies:createSolver(field, externalField)
   -- Run the KineticSpecies 'createSolver()' to initialize the collisions solver.
   VlasovSRSpecies.super.createSolver(self, field, externalField)

   -- Initialize velocity-space arrays for relativistic Vlasov
   -- Only need to do this once, so we don't need to store the updater
   local initSRVarsCalc = Updater.CalcSRVars {
      velGrid = self.velGrid, 
      confBasis = self.confBasis, 
      velBasis = self.velBasis, 
      phaseBasis = self.basis, 
      op = "init", 
   }
   initSRVarsCalc:advance(0.0, {}, {self.p_over_gamma, self.gamma, self.gamma_inv})

   local plasmaE, plasmaB = field:hasEB()
   local extHasE, extHasB = externalField:hasEB()

   self.totalEmField = self:allocVectorMoment(8)     -- 8 components of EM field.

   if self.hasExtForce then
      self.totEmFptr, self.totEmFidxr = self.totalEmField:get(1), self.totalEmField:genIndexer()
   end

   -- Create updater to advance solution by one time-step.
   if self.evolveCollisionless then
      self.solver = Updater.VlasovDG {
         onGrid    = self.grid,                      confBasis = self.confBasis,                phaseBasis = self.basis, 
         confRange = self.totalEmField:localRange(), velRange = self.p_over_gamma:localRange(), phaseRange = self.distf[1]:localRange(), 
         model_id  = self.model_id,                  field_id = self.field_id, 
      }
      self.collisionlessAdvance = function(tCurr, inFlds, outFlds)
         self.solver:advance(tCurr, inFlds, outFlds)
      end
   else
      self.solver = {totalTime = 0.}
      self.collisionlessAdvance = function(tCurr, inFlds, outFlds) end
   end
   -- Create updaters to compute various moments.
   self.numDensityCalc = Updater.DistFuncMomentDG {
      onGrid    = self.grid,                      confBasis = self.confBasis,                phaseBasis = self.basis, 
      confRange = self.totalEmField:localRange(), velRange = self.p_over_gamma:localRange(), 
      model_id  = self.model_id,                  isIntegrated = false, 
      moment    = "M0", 
   }
   self.momDensityCalc = Updater.DistFuncMomentDG {
      onGrid    = self.grid,                      confBasis = self.confBasis,                phaseBasis = self.basis, 
      confRange = self.totalEmField:localRange(), velRange = self.p_over_gamma:localRange(), 
      model_id  = self.model_id,                  isIntegrated = false, 
      moment    = "M1i", 
   }
   self.ptclEnergyCalc = Updater.DistFuncMomentDG {
      onGrid    = self.grid,                      confBasis = self.confBasis,                phaseBasis = self.basis, 
      confRange = self.totalEmField:localRange(), velRange = self.p_over_gamma:localRange(), 
      model_id  = self.model_id,                  isIntegrated = false, 
      moment    = "Energy", 
   }
   self.NiCalc = Updater.DistFuncMomentDG {
      onGrid    = self.grid,                      confBasis = self.confBasis,                phaseBasis = self.basis, 
      confRange = self.totalEmField:localRange(), velRange = self.p_over_gamma:localRange(), 
      model_id  = self.model_id,                  isIntegrated = false, 
      moment    = "Ni", 
   }
   -- Create updater for the coupling moments 
   -- Computes four-momentum density, Ni, then sets M0 and M1i from the result
   self.calcSelfCouplingMom = function(tCurr, fIn)
      self.NiCalc:advance(tCurr, {fIn, self.p_over_gamma, self.gamma, self.gamma_inv, self.V_drift, self.GammaV2, self.GammaV_inv}, {self.Ni})
      -- Copy number density (GammaV*n) and momentum density (GammaV*n*Vi) to their own fields.
      self.momDensity:combineOffset(1., self.Ni, 1*self.confBasis:numBasis())
      self.numDensity:combineOffset(1., self.Ni, 0)
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

   self.tmCouplingMom = 0.0    -- For timer.
end

function VlasovSRSpecies:setActiveRKidx(rkIdx) self.activeRKidx = rkIdx end

function VlasovSRSpecies:advance(tCurr, population, emIn, inIdx, outIdx)
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

   self.collisionlessAdvance(tCurr, {fIn, totalEmField, self.p_over_gamma}, {fRhsOut, self.cflRateByCell})

end

function VlasovSRSpecies:calcCouplingMoments(tCurr, rkIdx, species)
   -- Compute moments needed in coupling to fields and collisions.
   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]

   -- Compute Ni.
   self.calcSelfCouplingMom(tCurr, fIn)

   self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart
end

function VlasovSRSpecies:fluidMoments() return self.fiveMoments end

function VlasovSRSpecies:getNumDensity(rkIdx)
   -- If no rkIdx specified, assume numDensity has already been calculated.
   if rkIdx == nil then return self.numDensity end 

   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]
   self.numDensityCalc:advance(nil, {fIn}, { self.numDensity })

   self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart

   return self.numDensity
end

function VlasovSRSpecies:getMomDensity(rkIdx)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then return self.momDensity end 

   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]
   self.momDensityCalc:advance(nil, {fIn}, { self.momDensity })

   self.tmCouplingMom = self.tmCouplingMom + Time.clock() - tmStart

   return self.momDensity
end

function VlasovSRSpecies:momCalcTime()
   local tm = self.tmCouplingMom
   return tm
end

-- ................... Classes meant as aliases to simplify input files ...................... --
-- Special Relativistic Vlasov-Maxwell (Cartesian geometry)
local VlasovSRMaxwellSpecies = Proto(VlasovSRSpecies)
function VlasovSRMaxwellSpecies:fullInit(mySpecies)
   self.model_id  = "GKYL_MODEL_SR"
   self.field_id  = "GKYL_FIELD_E_B"
   VlasovSRMaxwellSpecies.super.fullInit(self, mySpecies)
end

-- Neutral Special Relativistic Vlasov (Cartesian geometry)
local VlasovSRNeutralSpecies = Proto(VlasovSRSpecies)
function VlasovSRNeutralSpecies:fullInit(mySpecies)
   self.model_id  = "GKYL_MODEL_SR"
   self.field_id  = "GKYL_FIELD_NULL"
   VlasovSRNeutralSpecies.super.fullInit(self, mySpecies)
end

-- ................... End of VlasovSRSpecies alias classes .................... --

return {VlasovSRMaxwell = VlasovSRMaxwellSpecies,
        VlasovSRNeutral = VlasovSRNeutralSpecies}
