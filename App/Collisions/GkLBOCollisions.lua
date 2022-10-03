-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Gyrokinetic LB Collision operator.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local Constants      = require "Lib.Constants"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local xsys           = require "xsys"
local lume           = require "Lib.lume"
local ffi            = require "ffi"

-- GkLBOCollisions ---------------------------------------------------------------
--
-- Lenard-Bernstein Collision operator.
-- Actually dates back to Lord Rayleigh, Philos. Mag. 32, 424 (1891).
-- Really LBO=the Dougherty operator.
--------------------------------------------------------------------------------

local GkLBOCollisions = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkLBOCollisions:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- previously stored table

   self.selfCollisions = true -- MF: to be deleted.
   self.collKind = "GkLBO"    -- Type of collisions model. Useful at the species app level.

   -- For now only cell-wise constant nu is implemented.
   self.cellConstNu = true     -- Cell-wise constant nu?

   self.collidingSpecies = assert(tbl.collideWith, "App.GkLBOCollisions: Must specify names of species to collide with in 'collideWith'.")

   -- Determine if cross-species collisions take place,
   -- and put the names of the other colliding species in a list.
   local selfSpecInd = lume.find(self.collidingSpecies, self.speciesName)
   assert(selfSpecInd, "App.GkLBOCollisions: must include self-collisions.")

   if #self.collidingSpecies > 1 then
      self.crossCollisions = true             -- Apply cross-species collisions.
      self.crossSpecies    = lume.clone(self.collidingSpecies)
      table.remove(self.crossSpecies, selfSpecInd)
   else
      self.crossCollisions = false            -- Don't apply cross-species collisions.
   end
   
   self.collFreqs = tbl.frequencies -- List of collision frequencies.
   if self.collFreqs then
      -- Collisionality, provided by user, will remain constant in time.
      self.timeDepNu = false
      self.calcSelfNu = function(momsIn, nuOut) GkLBOCollisions['calcSelfNuTimeConst'](self,momsIn,nuOut) end
      self.calcCrossNu = self.crossCollisions
         and function(species, otherNm, mOther, momsOther,
                      vtSqOther, nuCrossSelf, nuCrossOther)
            GkLBOCollisions['calcCrossNuTimeConst'](self,species, otherNm,
              mOther, momsOther, vtSqOther, nuCrossSelf, nuCrossOther)
         end
         or nil

      -- Ensure that collFreqs inputs are numbers or functions.
      for iC = 1,#self.collFreqs do
         local collFreqType = type(self.collFreqs[iC])
         assert(collFreqType=="number" or collFreqType=="function",
            "App.GkLBOCollisions: frequencies must either all be numbers, or all be functions")
         if (collFreqType == "number") then
            local val = self.collFreqs[iC]
            self.collFreqs[iC] = function(t, xn) return val end
         end
      end

      self.collFreqSelf = self.collFreqs[selfSpecInd]
      if self.crossCollisions then
         self.collFreqCross = lume.clone(self.collFreqs)
         table.remove(self.collFreqCross, selfSpecInd)
      end
   else
      -- Collisionality not provided by user. It will be calculated in time.
      self.timeDepNu = true
      self.calcSelfNu = function(momsIn, nuOut) GkLBOCollisions['calcSelfNuTimeDep'](self,momsIn,nuOut) end
      self.calcCrossNu = self.crossCollisions
         and function(species, otherNm, mOther, momsOther,
                      vtSqOther, nuCrossSelf, nuCrossOther)
            GkLBOCollisions['calcCrossNuTimeDep'](self,species, otherNm,
              mOther, momsOther, vtSqOther, nuCrossSelf, nuCrossOther)
         end
         or nil

      self.charge = speciesTbl.charge    -- Charge of this species.
      -- If no time-constant collision frequencies provided ('frequencies'), user can specify
      -- 'normNu' list of collisionalities normalized by T_0^(3/2)/n_0 evaluated somewhere in the
      -- simulation (see Gkeyll website for exact normalization). Otherwise code compute Spitzer
      -- collisionality from scratch.
      self.normNuIn = tbl.normNu
      -- normNuSelf, epsilon0 and elemCharge may not used, but are
      -- initialized to avoid if-statements in advance method.
      if self.normNuIn then
         self.normNuSelf = self.normNuIn[selfSpecInd]
         if self.crossCollisions then
            local normNuCrossIn = lume.clone(self.normNuIn)
            table.remove(normNuCrossIn, selfSpecInd)
            self.normNuCross = {}  -- Need a name-value pairs table.
            for i, nm in ipairs(self.crossSpecies) do self.normNuCross[nm] = normNuCrossIn[i] end
         end
      end
      -- Check for constants epsilon_0, elementary charge e, and Planck's constant/2pi. If not use default value.
      self.epsilon0   = tbl.epsilon0 and tbl.epsilon0 or Constants.EPSILON0
      self.elemCharge = tbl.elemCharge and tbl.elemCharge or Constants.ELEMENTARY_CHARGE
      self.hBar       = tbl.hBar and tbl.hBar or Constants.PLANCKS_CONSTANT_H/(2.0*Constants.PI)
   end

   if self.crossCollisions then
      self.charge     = speciesTbl.charge    -- Charge of this species.
      -- Can specify 'betaGreene' free parameter in Grene cross-species collisions.
      self.betaGreene = tbl.betaGreene and tbl.betaGreene or 0.0
   end

   self.mass = speciesTbl.mass   -- Mass of this species.

   self.nuFrac = tbl.nuFrac and tbl.nuFrac or 1.0

   self.cfl = 0.0    -- Will be replaced.

   self.timers = {nonSlvr = 0.}
end

function GkLBOCollisions:setName(nm) self.name = self.speciesName.."_"..nm end
function GkLBOCollisions:setSpeciesName(nm) self.speciesName = nm end
function GkLBOCollisions:setCfl(cfl) self.cfl = cfl end
function GkLBOCollisions:setConfBasis(basis) self.confBasis = basis end
function GkLBOCollisions:setConfGrid(grid) self.confGrid = grid end
function GkLBOCollisions:setPhaseBasis(basis) self.phaseBasis = basis end
function GkLBOCollisions:setPhaseGrid(grid) self.phaseGrid = grid end

function GkLBOCollisions:createSolver(mySpecies, externalField)
   local vDim = self.phaseGrid:ndim() - self.confGrid:ndim()

   -- Inverse of background magnetic field.
   self.bmag    = externalField.geo.bmag
   -- Inverse of background magnetic field.
   self.bmagInv = externalField.geo.bmagInv
      
   -- Self-species collisionality, which varies in space.
   self.nuSelf = mySpecies:allocMoment()
   -- Allocate fields to store self-species primitive moments.
   self.uParSelf = mySpecies:allocMoment()
   self.vtSqSelf = mySpecies:allocMoment()
   -- Allocate fields for boundary corrections.
   self.boundCorrs = mySpecies:allocVectorMoment(2)

   local vbounds = ffi.new("double[4]")
   for i = 1, vDim do
      vbounds[i-1]      = self.phaseGrid:lower(self.confGrid:ndim()+i)
      vbounds[i-1+vDim] = self.phaseGrid:upper(self.confGrid:ndim()+i)
   end
   self.primMomSelf = Updater.SelfPrimMoments {
      onGrid     = self.phaseGrid,   operator = "GkLBO",
      phaseBasis = self.phaseBasis,  vbounds  = vbounds,
      confBasis  = self.confBasis,   mass     = self.mass,
   }

   local projectUserNu
   if self.timeDepNu then
      self.m0Self = mySpecies:allocMoment()  -- M0, to be extracted from fiveMoments.
      -- Updater to compute spatially varying (Spitzer) nu.
      self.spitzerNu = Updater.SpitzerCollisionality {
         onGrid           = self.confGrid,     elemCharge = self.elemCharge,
         confBasis        = self.confBasis,    epsilon0   = self.epsilon0,
         useCellAverageNu = self.cellConstNu,  hBar       = self.hBar,
         willInputNormNu  = self.normNuIn,     nuFrac     = self.nuFrac,
      }
   else
      projectUserNu = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = self.collFreqSelf,
         basis  = self.confBasis,  onGhosts = false
      }
      projectUserNu:advance(0.0, {}, {self.nuSelf})
   end

   -- Weak multiplication to multiply nu(x) with uPar or vtSq.
   self.confMul = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "Multiply",
   }
   self.collisionSlvr = Updater.GkLBO {
      onGrid     = self.phaseGrid,   confBasis = self.confBasis,
      phaseBasis = self.phaseBasis,  confRange = self.nuSelf:localRange(),
      mass       = self.mass,
   }

   if self.crossCollisions then
      -- Cross-collision u and vtSq multiplied by collisionality.
      self.nuUParCross = mySpecies:allocMoment()
      self.nuVtSqCross = mySpecies:allocMoment()
      -- Prefactor m_0s*delta_s in cross primitive moment calculation.
      self.m0s_deltas     = mySpecies:allocMoment()
      self.m0s_deltas_den = mySpecies:allocMoment()
      -- Weak division to compute the pre-factor in cross collision primitive moments.
      self.confDiv = Updater.CartFieldBinOp {
         weakBasis = self.confBasis,           operation = "Divide",
         onRange   = self.nuSelf:localRange(),  onGhosts  = false,
      }
      -- Updater to compute cross-species primitive moments.
      self.primMomCross = Updater.CrossPrimMoments {
         onGrid     = self.confGrid,    betaGreene       = self.betaGreene,
         phaseBasis = self.phaseBasis,  varyingNu        = true,
         confBasis  = self.confBasis,   useCellAverageNu = self.cellConstNu,
         operator   = "GkLBO",
      }

      -- Allocate (and assign if needed) cross-species collision frequencies,
      -- and cross primitive moments.
      self.nuCross, self.uParCross, self.vtSqCross = {}, {}, {}
      -- Flag to indicate if cross collision frequency has been computed.
      self.crossFlags = self.timeDepNu and {} or nil
      for ispec, otherNm in ipairs(self.crossSpecies) do
         self.nuCross[otherNm] = mySpecies:allocMoment()
         if self.timeDepNu then
            self.crossFlags[otherNm] = false
         else
            projectUserNu:setFunc(self.collFreqCross[ispec])
            projectUserNu:advance(0.0, {}, {self.nuCross[otherNm]})
            self.nuCross[otherNm]:write(string.format("%s_nu-%s_%d.bp",self.speciesName,otherNm,0),0.0,0)
         end
         self.uParCross[otherNm] = mySpecies:allocMoment()
         self.vtSqCross[otherNm] = mySpecies:allocMoment()
      end

      self.m0Other = self.timeDepNu and mySpecies:allocMoment() or nil  -- M0, to be extracted from fiveMoments.
   end

   -- Collisionality, nu, summed over all species pairs.
   self.nuSum = mySpecies:allocMoment()
   -- Sum of flow velocities multiplied by respective collisionalities.
   self.nuUParSum = mySpecies:allocMoment()
   -- Sum of squared thermal speeds, vthSq=T/m, multiplied by respective collisionalities.
   self.nuVtSqSum = mySpecies:allocMoment()

   self.m2Self = mySpecies:allocMoment() -- M2, to be extracted from threeMoments.
end

function GkLBOCollisions:createCouplingSolver(species, field, externalField)
   -- Store a pointer to the collision app in the other species, so we know
   -- where to find things stored in the collision app (e.g. primitive momemts, nu).
   if self.crossCollisions then
      self.collAppOther = {}
      for _, nm in ipairs(self.crossSpecies) do
         for _, app in pairs(species[nm].collisions) do
            if app.collKind == self.collKind then self.collAppOther[nm] = app end
         end
      end
   end
end

function GkLBOCollisions:boundaryCorrections() return self.boundCorrs end
function GkLBOCollisions:selfPrimitiveMoments() return self.uParSelf, self.vtSqSelf end
function GkLBOCollisions:crossFrequencies(speciesName) return self.nuCross[speciesName] end
function GkLBOCollisions:crossNormNu(speciesName) return self.normNuCross[speciesName] end
function GkLBOCollisions:crossFlags(speciesName) return self.crossFlags[speciesName] end

function GkLBOCollisions:createDiagnostics(mySpecies, field)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = GkLBODiagsImpl()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

function GkLBOCollisions:calcCouplingMoments(tCurr, rkIdx, species)
   -- Compute self-primitive moments u and vtSq.
   local fIn      = species[self.speciesName]:rkStepperFields()[rkIdx]
   local momsSelf = species[self.speciesName]:fluidMoments()

   self.primMomSelf:advance(tCurr, {momsSelf, fIn, self.boundCorrs},
                            {self.uParSelf, self.vtSqSelf})
end

function GkLBOCollisions:calcSelfNuTimeConst(momsSelf, nuOut) nuOut:copy(self.nuSelf) end

function GkLBOCollisions:calcSelfNuTimeDep(momsSelf, nuOut)
   -- Compute the Spitzer collisionality.
   self.m0Self:combineOffset(1., momsSelf, 0)
   self.spitzerNu:advance(tCurr, {self.charge, self.mass, self.m0Self, self.vtSqSelf,
                                  self.charge, self.mass, self.m0Self, self.vtSqSelf,
                                  self.normNuSelf, self.bmag}, {nuOut})
end

function GkLBOCollisions:calcCrossNuTimeConst(species, otherNm,
   mOther, momsOther, vtSqOther, nuCrossSelf, nuCrossOther) end

function GkLBOCollisions:calcCrossNuTimeDep(species, otherNm,
   mOther, momsOther, vtSqOther, nuCrossSelf, nuCrossOther)

   -- Compute the Spitzer collisionality if another species hasn't already done so.
   local chargeOther     = species[otherNm]:getCharge()
   local crossFlagsSelf  = self.crossFlags[otherNm]
   local crossFlagsOther = self.collAppOther[otherNm]:crossFlags(self.speciesName)
   self.m0Other:combineOffset(1., momsOther, 0)
   if not crossFlagsSelf then
      local crossNormNuSelf = self.normNuCross[otherNm]
      self.spitzerNu:advance(tCurr, {self.charge, self.mass, self.m0Self, self.vtSqSelf,
                                     chargeOther, mOther, self.m0Other, vtSqOther, crossNormNuSelf},
                                    {nuCrossSelf})
      crossFlagsSelf = true
   end
   if not crossFlagsOther then
      local crossNormNuOther = self.collAppOther[otherNm]:crossNormNu(self.speciesName)
      self.spitzerNu:advance(tCurr, {chargeOther, mOther, self.m0Other, vtSqOther,
                                     self.charge, self.mass, self.m0Self, self.vtSqSelf, crossNormNuOther},
                                    {nuCrossOther})
      crossFlagsOther = true
   end
end

function GkLBOCollisions:advance(tCurr, fIn, species, out)
   local tmNonSlvrStart = Time.clock()

   local fRhsOut = out[1]
   local cflRateByCell = out[2]

   -- Fetch coupling moments of this species.
   local momsSelf = species[self.speciesName]:fluidMoments()

   self.calcSelfNu(momsSelf, self.nuSum)
   self.confMul:advance(tCurr, {self.nuSum, self.uParSelf}, {self.nuUParSum})
   self.confMul:advance(tCurr, {self.nuSum, self.vtSqSelf}, {self.nuVtSqSum})

   if self.crossCollisions then

      local bCorrectionsSelf = self.boundCorrs

      for _, otherNm in ipairs(self.crossSpecies) do

         local mOther               = species[otherNm]:getMass()
         local momsOther            = species[otherNm]:fluidMoments()
         local uParOther, vtSqOther = self.collAppOther[otherNm]:selfPrimitiveMoments()

         local nuCrossSelf  = self.nuCross[otherNm]
         local nuCrossOther = self.collAppOther[otherNm]:crossFrequencies(self.speciesName)

         -- Calculate time-dependent collision frequency if needed.
         self.calcCrossNu(species, otherNm, mOther, momsOther, vtSqOther,
                          nuCrossSelf, nuCrossOther)


         -- Compose the pre-factor:
         --   m0_s*delta_s = m0_s*(2*m_r*m0_r*nu_rs/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs))
         local deltas_num, deltas_den = self.m0s_deltas, self.m0s_deltas_den
         self.confMul:advance(tCurr, {momsSelf, nuCrossSelf, 1}, {deltas_den})
         self.confMul:advance(tCurr, {momsOther, nuCrossOther, 1}, {deltas_num})
         deltas_den:scale(self.mass)
         deltas_den:accumulate(mOther, deltas_num)
         deltas_num:scale(2.*mOther)
         self.confMul:advance(tCurr, {momsSelf, deltas_num, 1}, {deltas_num})
         self.confDiv:advance(tCurr, {deltas_den, deltas_num}, {self.m0s_deltas})

         -- Cross-primitive moments for the collision of these two species has not been computed.
         self.primMomCross:advance(tCurr, {self.mass, nuCrossSelf, momsSelf, self.uParSelf, self.vtSqSelf, bCorrectionsSelf,
                                           mOther, nuCrossOther, momsOther, uParOther, vtSqOther, bCorrectionsOther,
                                           self.m0s_deltas},
                                          {self.uParCross[otherNm], self.vtSqCross[otherNm]})

         self.confMul:advance(tCurr, {nuCrossSelf, self.uParCross[otherNm]}, {self.nuUParCross})
         self.confMul:advance(tCurr, {nuCrossSelf, self.vtSqCross[otherNm]}, {self.nuVtSqCross})

         self.nuSum:accumulate(1.0, nuCrossSelf)
         self.nuUParSum:accumulate(1.0, self.nuUParCross)
         self.nuVtSqSum:accumulate(1.0, self.nuVtSqCross)

      end    -- end loop over other species that this species collides with.

   end    -- end if self.crossCollisions.
   self.timers.nonSlvr = self.timers.nonSlvr + Time.clock() - tmNonSlvrStart

   -- M2 self is needed as a precaution in GkLBO.
   self.m2Self:combineOffset(1., momsSelf, 2*self.confBasis:numBasis())

   -- Compute increment from collisions and accumulate it into output.
   self.collisionSlvr:advance(
      tCurr, {fIn, self.bmagInv, self.nuUParSum, self.nuVtSqSum, self.nuSum, self.m2Self}, {fRhsOut, cflRateByCell})
end

function GkLBOCollisions:write(tm, frame) end

function GkLBOCollisions:totalTime() return self.collisionSlvr.totalTime + self.timers.nonSlvr end

function GkLBOCollisions:slvrTime() return self.collisionSlvr.totalTime end

function GkLBOCollisions:nonSlvrTime() return self.timers.nonSlvr end

return GkLBOCollisions
