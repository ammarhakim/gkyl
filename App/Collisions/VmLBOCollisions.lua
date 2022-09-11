-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Vlasov LB Collision operator
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local Constants      = require "Lib.Constants"
local DataStruct     = require "DataStruct"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local xsys           = require "xsys"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"
local lume           = require "Lib.lume"
local DiagsApp       = require "App.Diagnostics.SpeciesDiagnostics"
local DiagsImplBase  = require "App.Diagnostics.DiagnosticsImplBase"

-- ............... IMPLEMENTATION OF DIAGNOSTICS ................. --
-- Diagnostics could be placed in a separate file if they balloon in
-- number. But if we only have one or two we can just place it here.

-- ~~~~ Source integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
local VmLBODiagsImpl = function()
   -- IMPORTANT: this diagnostic is only here for testing!!! do not use (MF).
   local _collOut = Proto(DiagsImplBase)
   function _collOut:fullInit(diagApp, mySpecies, fieldIn, owner)
      self.field = mySpecies:allocDistf()
      self.owner = owner 
      self.done  = false
   end
   function _collOut:getType() return "grid" end
   function _collOut:advance(tm, inFlds, outFlds)
      local specIn = inFlds[1]
      self.field:copy(self.owner.collOut)
   end

   return {collOut = _collOut}
end

-- .................... END OF DIAGNOSTICS ...................... --

-- VmLBOCollisions ---------------------------------------------------------------
--
-- Lenard-Bernstein Collision operator.
-- Actually dates back to Lord Rayleigh, Philos. Mag. 32, 424 (1891).
-- Really LBO=the Dougherty operator.
--------------------------------------------------------------------------------

local VmLBOCollisions = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmLBOCollisions:init(tbl)
   self.tbl = tbl
end

-- Function to find the index of an element in table.
local function findInd(tblIn, el)
   for i, v in ipairs(tblIn) do
      if v == el then
         return i
      end
   end
   return #tblIn+1    -- If not found return a number larger than the length of the table.
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.collKind = "VmLBO"    -- Type of collisions model. Useful at the species app level.

   self.collidingSpecies = assert(tbl.collideWith, "App.VmLBOCollisions: Must specify names of species to collide with in 'collideWith'.")

   self.selfCollisions = true -- MF: to be deleted.
   self.varNu          = true -- MF: to be deleted.   

   -- First determine if self-species and/or cross-species collisions take place,
   -- and (if cross-collisions=true) put the names of the other colliding species in a list.
   local selfSpecInd = findInd(self.collidingSpecies, self.speciesName)
   if selfSpecInd < (#self.collidingSpecies+1) then
      if #self.collidingSpecies > 1 then
         self.crossCollisions = true             -- Apply cross-species collisions.
         self.crossSpecies    = lume.clone(self.collidingSpecies)
         table.remove(self.crossSpecies, selfSpecInd)
      else
         self.crossCollisions = false            -- Don't apply cross-species collisions.
      end
   else
      assert(false, "App.VmLBOCollisions: must include self-collisions.")
   end

   self.collFreqs = tbl.frequencies -- List of collision frequencies, if using spatially constant nu.
   if self.collFreqs then
      -- Collisionality, provided by user, will remain constant in time.
      self.timeDepNu = false

      -- Ensure that collFreqs inputs are numbers or functions.
      for iC = 1,#self.collFreqs do
         local collFreqType = type(self.collFreqs[iC])
         assert(collFreqType=="number" or collFreqType=="function",
            "App.VmLBOCollisions: frequencies must either all be numbers, or all be functions")
         if (collFreqType == "number") then
            local val = self.collFreqs[iC]
            self.collFreqs[iC] = function(t, xn) return val end
         end
      end

      -- For now only cell-wise constant nu is implemented.
      self.cellConstNu  = true     -- Cell-wise constant nu?

      self.collFreqSelf = self.collFreqs[selfSpecInd]
      if self.crossCollisions then
         self.collFreqCross = lume.clone(self.collFreqs)
         table.remove(self.collFreqCross, selfSpecInd)
      end
   else
      -- Collisionality not provided by user. It will be calculated in time.
      self.timeDepNu = true

      self.mass         = speciesTbl.mass      -- Mass of this species.
      self.charge       = speciesTbl.charge    -- Charge of this species.
      -- For now only cell-wise constant nu is implemented.
      self.cellConstNu  = true     -- Cell-wise constant nu?
      -- If no time-constant collision frequencies provided ('frequencies'), user can specify
      -- 'normNu' list of collisionalities normalized by T_0^(3/2)/n_0 evaluated somewhere in the
      -- simulation (see Gkeyll website for exact normalization). Otherwise code compute Spitzer
      -- collisionality from scratch.
      self.normNuIn     = tbl.normNu
      -- normNuSelf, epsilon0 and elemCharge may not used, but are
      -- initialized to avoid if-statements in advance method.
      if self.normNuIn then
         self.userInputNormNu = true
         self.normNuSelf      = self.normNuIn[selfSpecInd]
         if self.crossCollisions then
            self.normNuCross = lume.clone(self.normNuIn)
            table.remove(self.normNuCross, selfSpecInd)
         end
      else
         self.userInputNormNu = false
         self.normNuSelf      = 0.0
         if self.crossCollisions then
            self.normNuCross = lume.clone(self.collidingSpecies)
            table.remove(self.normNuCross, selfSpecInd)
            for i, _ in ipairs(self.normNuCross) do self.normNuCross[i] = 0.0 end 
         end
      end
      -- Check for constants epsilon_0, elementary charge e, and Planck's constant/2pi. If not use default value.
      self.epsilon0   = tbl.epsilon0 and tbl.epsilon0 or Constants.EPSILON0
      self.elemCharge = tbl.elemCharge and tbl.elemCharge or Constants.ELEMENTARY_CHARGE
      self.hBar       = tbl.hBar and tbl.hBar or Constants.PLANCKS_CONSTANT_H/(2.0*Constants.PI)
   end

   if self.crossCollisions then
      self.mass       = speciesTbl.mass      -- Mass of this species.
      self.charge     = speciesTbl.charge    -- Charge of this species.
      -- Can specify 'betaGreene' free parameter in Grene cross-species collisions.
      self.betaGreene = tbl.betaGreene and tbl.betaGreene or 0.0
   end

   self.nuFrac = tbl.nuFrac and tbl.nuFrac or 1.0

   self.cfl = 0.0    -- Will be replaced.

   self.timers = {nonSlvr = 0.}
end

function VmLBOCollisions:setName(nm) self.name = self.speciesName.."_"..nm end
function VmLBOCollisions:setSpeciesName(nm) self.speciesName = nm end
function VmLBOCollisions:setCfl(cfl) self.cfl = cfl end
function VmLBOCollisions:setConfBasis(basis) self.confBasis = basis end
function VmLBOCollisions:setConfGrid(grid) self.confGrid = grid end
function VmLBOCollisions:setPhaseBasis(basis) self.phaseBasis = basis end
function VmLBOCollisions:setPhaseGrid(grid) self.phaseGrid = grid end

function VmLBOCollisions:createSolver(mySpecies, extField)
   self.vdim      = self.phaseGrid:ndim() - self.confGrid:ndim()

   self.cNumBasis = self.confBasis:numBasis()

   -- Maximum velocity of the velocity grid (and its square).
   self.vMax   = Lin.Vec(self.vdim)
   for vd = 1,self.vdim do
      self.vMax[vd]   = self.phaseGrid:upper(self.confGrid:ndim()+vd)
   end
   self.vMaxSq = self.vMax[1] 
   for vd = 1,self.vdim do
      if (self.vMaxSq < self.vMax[vd]) then
         self.vMaxSq = self.vMax[vd]
      end
   end
   self.vMaxSq = self.vMaxSq^2

   -- Collisionality, nu, summed over all species pairs.
   self.nuSum = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.cNumBasis,
      ghost         = {1, 1},
   }
   -- Sum of flow velocities in vdim directions multiplied by respective collisionalities.
   self.nuUSum = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.cNumBasis*self.vdim,
      ghost         = {1, 1},
   }
   -- Sum of squared thermal speeds, vthSq=T/m, multiplied by respective collisionalities.
   self.nuVtSqSum = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.cNumBasis,
      ghost         = {1, 1},
   }

   -- Zero-flux BCs in the velocity dimensions.
   local zfd = { }
   for d = 1, self.vdim do
      zfd[d] = self.confGrid:ndim() + d
   end

   -- Self-species collisionality, which varies in space.
   self.nuVarXSelf = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.cNumBasis,
      ghost         = {1, 1},
   }
   if self.timeDepNu then 
      -- Updater to compute spatially varying (Spitzer) nu.
      self.spitzerNu = Updater.SpitzerCollisionality {
         onGrid           = self.confGrid,         elemCharge = self.elemCharge,
         confBasis        = self.confBasis,        epsilon0   = self.epsilon0,
         useCellAverageNu = self.cellConstNu,      hBar       = self.hBar,
         willInputNormNu  = self.userInputNormNu,  nuFrac     = self.nuFrac,
      }
   else
      local projectUserNu = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = self.collFreqSelf,
         basis  = self.confBasis,  onGhosts = false
      }
      projectUserNu:advance(0.0, {}, {self.nuVarXSelf})
   end
   -- Weak multiplication to multiply nu(x) with u or vtSq.
   self.confMul = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "Multiply",
   }
   -- Lenard-Bernstein operator (LBO) collisions updater.
   self.collisionSlvr = Updater.VlasovLBO {
      onGrid     = self.phaseGrid,   confBasis = self.confBasis,
      phaseBasis = self.phaseBasis,  confRange = self.nuUSum:localRange(),
   }

   if self.crossCollisions then
      -- Temporary collisionality fields.
      self.nuCrossSelf = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis,
         ghost         = {1, 1},
      }
      self.nuCrossOther = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis,
         ghost         = {1, 1},
      }
      -- Cross-collision u and vtSq multiplied by collisionality.
      self.nuUCross = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis*self.vdim,
         ghost         = {1, 1},
      }
      self.nuVtSqCross = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Prefactor m_0s*delta_s in cross primitive moment calculation.
      self.m0s_deltas = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis,
         ghost         = {1, 1},
      }
      self.m0s_deltas_den = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis,
         ghost         = {1, 1},
      }
      -- Weak division to compute the pre-factor in cross collision primitive moments.
      self.confDiv = Updater.CartFieldBinOp {
         weakBasis = self.confBasis,           operation = "Divide",
         onRange   = self.nuSum:localRange(),  onGhosts  = false,
      }
      -- Updater to compute cross-species primitive moments.
      self.primMomCross = Updater.CrossPrimMoments {
         onGrid           = self.confGrid,    betaGreene       = self.betaGreene, 
         phaseBasis       = self.phaseBasis,  varyingNu        = self.varNu,
         confBasis        = self.confBasis,   useCellAverageNu = self.cellConstNu,
         operator         = "VmLBO",
      }
   end

   -- Number of cells in which number density was negative (somewhere).
   self.primMomLimitCrossings = DataStruct.DynVector {
      numComponents = 1,
   }
   self.primMomCrossLimitL = Lin.Vec(1)
   self.primMomCrossLimitG = Lin.Vec(1)
   -- Factor dividing zeroth-coefficient in configuration space cell average.
   self.cellAvFac          = 1.0/math.sqrt(2.0^self.confGrid:ndim())
end

function VmLBOCollisions:createDiagnostics(mySpecies, field)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = VmLBODiagsImpl()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

function VmLBOCollisions:advance(tCurr, fIn, species, out)

   local fRhsOut = out[1]
   local cflRateByCell = out[2]

   local tmNonSlvrStart = Time.clock()

   -- NOTE: The following code is commented out because Vm users don't seem
   -- to be as worried about limit crossings as Gk users, so counting them
   -- is disabled for now. See the 'write' method as well.
   ---- Determine whether primitive moments cross limits based on
   ---- parallel flow speed and thermal speed squared.
   --self.primMomCrossLimitL[1] = 0
   --self.primMomCrossLimitG[1] = 0
   --local confIndexer          = self.velocity:genIndexer()
   --local uItr                 = self.velocity:get(1)
   --local vthSqItr             = self.vthSq:get(1)
   --for idx in self.velocity:localRangeIter() do
   --   self.velocity:fill(confIndexer(idx), uItr)
   --   self.vthSq:fill(confIndexer(idx), vthSqItr)
   --   local primCrossingFound = false
   --   for vd = 1,self.vdim do
   --      if (math.abs(uItr[(vd-1)*self.cNumBasis+1]*self.cellAvFac)>self.vMax[vd]) then
   --         uCrossingFound = true
   --         break
   --      end
   --   end
   --   local vthSq0 = vthSqItr[1]*self.cellAvFac

   --   if (uCrossingFound or (vthSq0<0) or (vthSq0>self.vMaxSq)) then
   --      self.primMomCrossLimitL[1] = self.primMomCrossLimitL[1]+1
   --   end
   --end
   --Mpi.Allreduce(self.primMomCrossLimitL:data(), self.primMomCrossLimitG:data(), 1,
   --              Mpi.DOUBLE, Mpi.SUM, self.confGrid:commSet().comm)
   --self.primMomLimitCrossings:appendData(tCurr+dt, self.primMomCrossLimitG)

   -- Fetch coupling moments and primitive moments of this species.
   local momsSelf    = species[self.speciesName]:fluidMoments()
   local primMomSelf = species[self.speciesName]:selfPrimitiveMoments()

   if self.timeDepNu then
      -- Compute the Spitzer collisionality.
      self.spitzerNu:advance(tCurr, {self.charge, self.mass, momsSelf[1], primMomSelf[2],
                                     self.charge, self.mass, momsSelf[1], primMomSelf[2], self.normNuSelf}, {self.nuSum})
   else
      self.nuSum:copy(self.nuVarXSelf)
   end
   self.confMul:advance(tCurr, {self.nuSum, primMomSelf[1]}, {self.nuUSum})
   self.confMul:advance(tCurr, {self.nuSum, primMomSelf[2]}, {self.nuVtSqSum})

   if self.crossCollisions then

      local bCorrectionsSelf = species[self.speciesName]:boundaryCorrections()

      for sInd, otherNm in ipairs(self.crossSpecies) do

         local mOther            = species[otherNm]:getMass()
         local momsOther         = species[otherNm]:fluidMoments()
         local primMomOther      = species[otherNm]:selfPrimitiveMoments()
         local bCorrectionsOther = species[otherNm]:boundaryCorrections()

         if self.timeDepNu then
            -- Compute the Spitzer collisionality if another species hasn't already done so.
            local chargeOther = species[otherNm]:getCharge()
            if (not species[self.speciesName].momentFlags[6][otherNm]) then
               self.spitzerNu:advance(tCurr, {self.charge, self.mass, momsSelf[1], primMomSelf[2],
                                              chargeOther, mOther, momsOther[1], primMomOther[2], self.normNuCross[sInd]},
                                             {species[self.speciesName].nuVarXCross[otherNm]})
               species[self.speciesName].momentFlags[6][otherNm] = true
            end
            if (not species[otherNm].momentFlags[6][self.speciesName]) then
               self.spitzerNu:advance(tCurr, {chargeOther, mOther, momsOther[1], primMomOther[2],
                                              self.charge, self.mass, momsSelf[1], primMomSelf[2], species[otherNm].collPairs[otherNm][self.speciesName].normNu},
                                             {species[otherNm].nuVarXCross[self.speciesName]})
               species[otherNm].momentFlags[6][self.speciesName] = true
            end
         end
         self.nuCrossSelf:copy(species[self.speciesName].nuVarXCross[otherNm])
         self.nuCrossOther:copy(species[otherNm].nuVarXCross[self.speciesName])

         -- Compose the pre-factor:
         --   m0_s*delta_s = m0_s*(2*m_r*m0_r*nu_rs/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs))
         local deltas_num, deltas_den = self.m0s_deltas, self.m0s_deltas_den
         self.confMul:advance(tCurr, {momsSelf, self.nuCrossSelf, 1}, {deltas_den})
         self.confMul:advance(tCurr, {momsOther, self.nuCrossOther, 1}, {deltas_num})
         deltas_den:scale(self.mass)
         deltas_den:accumulate(mOther, deltas_num)
         deltas_num:scale(2.*mOther)
         self.confMul:advance(tCurr, {momsSelf, deltas_num, 1}, {deltas_num})
         self.confDiv:advance(tCurr, {deltas_den, deltas_num}, {self.m0s_deltas})

         self.primMomCross:advance(tCurr, {self.mass, self.nuCrossSelf, momsSelf, primMomSelf, bCorrectionsSelf,
                                           mOther, self.nuCrossOther, momsOther, primMomOther, bCorrectionsOther,
                                           self.m0s_deltas},
                                          {species[self.speciesName].uCross[otherNm], species[self.speciesName].vtSqCross[otherNm]})

         self.confMul:advance(tCurr, {self.nuCrossSelf, species[self.speciesName].uCross[otherNm]}, {self.nuUCross})
         self.confMul:advance(tCurr, {self.nuCrossSelf, species[self.speciesName].vtSqCross[otherNm]}, {self.nuVtSqCross})

         self.nuSum:accumulate(1.0, self.nuCrossSelf)
         self.nuUSum:accumulate(1.0, self.nuUCross)
         self.nuVtSqSum:accumulate(1.0, self.nuVtSqCross)

      end    -- end loop over other species that this species collides with.

   end    -- end if self.crossCollisions.
   self.timers.nonSlvr = self.timers.nonSlvr + Time.clock() - tmNonSlvrStart

   -- Compute increment from collisions and accumulate it into output.
   self.collisionSlvr:advance(
      tCurr, {fIn, self.nuUSum, self.nuVtSqSum, self.nuSum}, {fRhsOut, cflRateByCell})

end

function VmLBOCollisions:write(tm, frame)
-- Since this doesn't seem to be as big a problem in Vm as in Gk, we comment this out for now.
--   self.primMomLimitCrossings:write(string.format("%s_%s.bp", self.speciesName, "primMomLimitCrossings"), tm, frame)
end

function VmLBOCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.timers.nonSlvr
end
function VmLBOCollisions:slvrTime()
   return self.collisionSlvr.totalTime
end
function VmLBOCollisions:nonSlvrTime()
   return self.timers.nonSlvr
end

return VmLBOCollisions
