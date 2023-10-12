-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: BGK Collision operator for gyrokinetics.
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
local Mpi            = require "Comm.Mpi"
local lume           = require "Lib.lume"
local xsys           = require "xsys"


-- GkBGKCollisions ---------------------------------------------------------------
--
-- Bhatnagar-Gross-Krook Collision operator.
--------------------------------------------------------------------------------

local GkBGKCollisions = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkBGKCollisions:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function GkBGKCollisions:fullInit(speciesTbl)
   local tbl = self.tbl  -- Previously stored table.

   self.collKind = "GkBGK"  -- Type of collisions model. Useful at the species app level.

   -- (MF 2022/11/18: not entirely true) For now only cell-wise constant nu is implemented.
   self.cellConstNu = true     -- Cell-wise constant nu?

   self.collidingSpecies = assert(tbl.collideWith,
      "App.GkBGKCollisions: Must specify names of species to collide with in 'collideWith'.")

   -- Determine if cross-species collisions take place,
   -- and put the names of the other colliding species in a list.
   local selfSpecInd = lume.find(self.collidingSpecies, self.speciesName)
   assert(selfSpecInd, "App.GkBGKCollisions: must include self-collisions.")

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
      self.calcSelfNu = function(momsIn, nuOut) GkBGKCollisions['calcSelfNuTimeConst'](self,momsIn,nuOut) end
      self.calcCrossNu = self.crossCollisions
         and function(otherNm, qOther, mOther, momsOther,
                      primMomsOther, nuCrossSelf, nuCrossOther)
            GkBGKCollisions['calcCrossNuTimeConst'](self, otherNm, qOther,
              mOther, momsOther, primMomsOther, nuCrossSelf, nuCrossOther)
         end
         or nil

      -- Ensure that collFreqs inputs are numbers or functions.
      for iC = 1,#self.collFreqs do
         local collFreqType = type(self.collFreqs[iC])
         assert(collFreqType=="number" or collFreqType=="function",
            "App.GkBGKCollisions: frequencies must either all be numbers, or all be functions")
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
      self.calcSelfNu = function(momsIn, nuOut) GkBGKCollisions['calcSelfNuTimeDep'](self,momsIn,nuOut) end
      self.calcCrossNu = self.crossCollisions
         and function(otherNm, qOther, mOther, momsOther,
                      primMomsOther, nuCrossSelf, nuCrossOther)
            GkBGKCollisions['calcCrossNuTimeDep'](self, otherNm, qOther,
              mOther, momsOther, primMomsOther, nuCrossSelf, nuCrossOther)
         end
         or nil

      self.mass   = speciesTbl.mass      -- Mass of this species.
      self.charge = speciesTbl.charge    -- Charge of this species.
      -- If no time-constant collision frequencies provided ('frequencies'), user can specify
      -- 'normNu' list of collisionalities normalized by T_0^(3/2)/n_0 evaluated somewhere in the
      -- simulation (see Gkeyll website for exact normalization). Otherwise code compute Spitzer
      -- collisionality from scratch.
      self.normNuIn = tbl.normNu
      -- normNuSelf and epsilon0 may not used, but are
      -- initialized to avoid if-statements in advance method.
      self.normNuCross = self.crossCollisions and {} or nil  -- Need a name-value pairs table.
      if self.normNuIn then
         self.normNuSelf = self.normNuIn[selfSpecInd]
         if self.crossCollisions then
            local normNuCrossIn = lume.clone(self.normNuIn)
            table.remove(normNuCrossIn, selfSpecInd)
            for i, nm in ipairs(self.crossSpecies) do self.normNuCross[nm] = normNuCrossIn[i] end
         end
      end
      -- Check for constants epsilon_0, elementary charge e, and Planck's constant/2pi. If not use default value.
      self.epsilon0 = tbl.epsilon0 and tbl.epsilon0 or Constants.EPSILON0
      self.hBar     = tbl.hBar and tbl.hBar or Constants.PLANCKS_CONSTANT_H/(2.0*Constants.PI)
   end

   if self.crossCollisions then
      self.mass       = speciesTbl.mass      -- Mass of this species.
      self.charge     = speciesTbl.charge    -- Charge of this species.
      -- Can specify 'betaGreene' free parameter in Grene cross-species collisions.
      self.betaGreene = tbl.betaGreene and tbl.betaGreene or 0.0
   end

   self.nuFrac = tbl.nuFrac and tbl.nuFrac or 1.0

   self.mass = speciesTbl.mass   -- Mass of this species.

   -- DL 2023/09/26: Replace the original Lagrange fix by the iteration fix. 
   self.conservativeMaxwellian = xsys.pickBool(tbl.conservativeMaxwellian, true)

   self.timers = {mom = 0.,   momcross = 0.,   advance = 0.,}
end

function GkBGKCollisions:setName(nm) self.name = self.speciesName.."_"..nm end
function GkBGKCollisions:setSpeciesName(nm) self.speciesName = nm end
function GkBGKCollisions:setCfl(cfl) self.cfl = cfl end
function GkBGKCollisions:setConfBasis(basis) self.confBasis = basis end
function GkBGKCollisions:setConfGrid(grid) self.confGrid = grid end
function GkBGKCollisions:setPhaseBasis(basis) self.phaseBasis = basis end
function GkBGKCollisions:setPhaseGrid(grid) self.phaseGrid = grid end

function GkBGKCollisions:createSolver(mySpecies, externalField)
   -- Background magnetic field. Needed for spatially varying nu
   -- or to project Maxwellians with vdim>1.
   self.bmag     = externalField.geo.bmag
   self.jacobTot = externalField.geo.jacobTot

   if self.conservativeMaxwellian then
      -- Intermediate moments used in Lagrange fixing.
      self.dmoms = mySpecies:allocVectorMoment(3)
      -- Create updater to compute M0, M1, M2 moments sequentially.
      self.threeMomentsCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.phaseGrid,    confBasis = self.confBasis,
         phaseBasis = self.phaseBasis,   moment    = "GkThreeMoments",
         gkfacs     = {self.mass, self.bmag},
      }
      self.correctMaxwellian = Updater.CorrectMaxwellian {
         onGrid     = self.phaseGrid,   confGrid  = self.confGrid,
         phaseBasis = self.phaseBasis,  confBasis = self.confBasis,
         mass       = self.mass,
         bmag       = self.bmag,        jacobTot   = self.jacobTot,
         iter_max   = 100,              err_max    = 1e-14,
         useDevice  = false,
      }
   end

   -- Self-species collisionality, which varies in space.
   self.nuSelf = mySpecies:allocMoment()
   -- Allocate fields to store self-species primitive moments.
   self.primMomsSelf = mySpecies:allocVectorMoment(2)
   -- Allocate fields for boundary corrections.
   self.boundCorrs = mySpecies:allocVectorMoment(2)

   -- We'll use the LBO SelfPrimMom updater with zero boundary corrections.
   self.primMomsSelfCalc = Updater.SelfPrimMoments {
      onGrid     = self.phaseGrid,   operator  = "GkBGK",
      phaseBasis = self.phaseBasis,  confBasis = self.confBasis,
      confRange = self.nuSelf:localRange(),
   }

   -- Sum of Maxwellians multiplied by respective collisionalities.
   self.nufMaxwellSum = mySpecies:allocDistf()

   local projectUserNu
   if self.timeDepNu then
      self.m0Self    = mySpecies:allocMoment()  -- M0, to be extracted from threeMoments.
      self.vtSqSelf  = mySpecies:allocMoment()
      self.vtSqOther = self.crossCollisions and mySpecies:allocMoment() or nil
      -- Updater to compute spatially varying (Spitzer) nu.
      self.spitzerNu = Updater.SpitzerCollisionality {
         confBasis       = self.confBasis,  epsilon0 = self.epsilon0,
         hBar            = self.hBar,       nuFrac   = self.nuFrac,
         willInputNormNu = self.normNuIn,
      }
   else
      projectUserNu = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = self.collFreqSelf,
         basis  = self.confBasis,  onGhosts = false
      }
      projectUserNu:advance(0.0, {}, {self.nuSelf})
   end

   -- Weak multiplication in conf space.
   self.confMul = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "Multiply",
   }
   -- Weak multiplication to multiply nu(x) with fMaxwell.
   self.phaseMul = Updater.CartFieldBinOp {
      weakBasis  = self.phaseBasis,  operation = "Multiply",
      fieldBasis = self.confBasis,
   }

   -- Maxwellian projection.
   self.maxwellian = Updater.MaxwellianOnBasis {
      onGrid    = self.phaseGrid,  phaseBasis  = self.phaseBasis,
      confBasis = self.confBasis,  --usePrimMoms = true,
      mass      = self.mass, 
   }
   -- BGK Collision solver itself.
   self.collisionSlvr = Updater.BGKcollisions {
      onGrid   = self.phaseGrid,  confBasis  = self.confBasis,
      confGrid = self.confGrid,   phaseBasis = self.phaseBasis,
   }

   if self.crossCollisions then
      -- Cross-collision Maxwellian multiplied by collisionality.
      self.nufMaxwellCross = mySpecies:allocDistf()
      -- Prefactor m_0s*delta_s in cross primitive moment calculation.
      self.m0s_deltas = mySpecies:allocMoment()
      -- Weak division to compute the pre-factor in cross collision primitive moments.
      self.confDiv = Updater.CartFieldBinOp {
         weakBasis = self.confBasis,            operation = "Divide",
         onRange   = self.nuSelf:localRange(),  onGhosts  = false,
      }
      -- Updater to compute cross-species moments.
      self.MomCrossBGK = Updater.MomCrossBGK {
         phaseBasis = self.phaseBasis,   confBasis  = self.confBasis,
      }
      self.unitFld = mySpecies:allocMoment()
      local projectUnit = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = function(t,xn) return 1. end,
         basis  = self.confBasis,  onGhosts = false
      }
      projectUnit:advance(0.0, {}, {self.unitFld})

      -- Allocate (and assign if needed) cross-species collision frequencies,
      -- and cross moments.
      self.nuCross = {}
      for ispec, otherNm in ipairs(self.crossSpecies) do
         self.nuCross[otherNm] = mySpecies:allocMoment()
         if not self.timeDepNu then
            projectUserNu:setFunc(self.collFreqCross[ispec])
            projectUserNu:advance(0.0, {}, {self.nuCross[otherNm]})
--            self.nuCross[otherNm]:write(string.format("%s_nu-%s_%d.bp",self.speciesName,otherNm,0),0.0,0)
         end
      end
      self.momsCross = mySpecies:allocVectorMoment(6)   

      self.m0Other = self.timeDepNu and mySpecies:allocMoment() or nil  -- M0, to be extracted from threeMoments.

      if self.conservativeMaxwellian then
         self.momsCrossSelf = mySpecies:allocVectorMoment(3)   
      end
   end

   -- Collisionality, nu, summed over all species pairs.
   self.nuSum = mySpecies:allocMoment()
end

function GkBGKCollisions:createCouplingSolver(population, field, externalField)
   -- Store a pointer to the collision app in the other species, so we know
   -- where to find things stored in the collision app (e.g. primitive momemts, nu).
   local species = population:getSpecies()
   if self.crossCollisions then
      self.collAppOther = {}
      for _, nm in ipairs(self.crossSpecies) do
         for _, app in pairs(species[nm].collisions) do
            if app.collKind == self.collKind then self.collAppOther[nm] = app end
         end
      end

      local isThisSpeciesMine = population:isSpeciesMine(self.speciesName)
      if not isThisSpeciesMine then
         local mySpecies = population:getSpecies()[self.speciesName]
         -- Allocate primMoms if we collide with other species not in this rank.
         self.primMomsSelf = mySpecies:allocVectorMoment(2)

         -- Allocate (and assign if needed) cross-species collision frequencies
         -- of species not in this rank.
         local projectUserNu = Updater.ProjectOnBasis {
            onGrid = self.confGrid,   evaluate = function(t, xn) return 0. end,
            basis  = self.confBasis,  onGhosts = false
         }
         self.nuCross = {}
         for ispec, otherNm in ipairs(self.crossSpecies) do
            self.nuCross[otherNm] = mySpecies:allocMoment()
            if not self.timeDepNu then
               projectUserNu:setFunc(self.collFreqCross[ispec])
               projectUserNu:advance(0.0, {}, {self.nuCross[otherNm]})
            end
         end
      end

      local messenger = self.confGrid:getMessenger()
      -- Create list of ranks we need to send/recv local self primitive moments to/from.
      -- MF: We'll merge u and vtSq into a single CartField, and merge these two Xfer objects.
      self.primMomsSelfXfer = {}
      self.primMomsSelfXfer.destRank, self.primMomsSelfXfer.srcRank = {}, {}
      self.primMomsSelfXfer.sendReqStat, self.primMomsSelfXfer.recvReqStat = {}, {}
      for _, sO in ipairs(self.crossSpecies) do
         local sOrank = population:getSpeciesOwner(sO)
         local selfRank = population:getSpeciesOwner(self.speciesName)
         if isThisSpeciesMine then
            -- Only species owned by this rank send primMoms other ranks.
            if (not lume.any(self.primMomsSelfXfer.destRank, function(e) return e==sOrank end)) and
               (not population:isSpeciesMine(sO)) then
               table.insert(self.primMomsSelfXfer.destRank, sOrank)
               table.insert(self.primMomsSelfXfer.sendReqStat, messenger:newRequestStatus())
            end
         else
            -- Only species not owned by this rank receive primMoms from other ranks.
            if (not lume.any(self.primMomsSelfXfer.srcRank, function(e) return e==selfRank end)) and
               (not population:isSpeciesMine(self.speciesName)) then
               table.insert(self.primMomsSelfXfer.srcRank, selfRank)
               table.insert(self.primMomsSelfXfer.recvReqStat, messenger:newRequestStatus())
            end
         end
      end
   else
      self.primMomsSelfXfer = {}
      self.primMomsSelfXfer.destRank, self.primMomsSelfXfer.srcRank = {}, {}
   end
end

function GkBGKCollisions:selfPrimitiveMoments() return self.primMomsSelf end
function GkBGKCollisions:crossFrequencies(speciesName) return self.nuCross[speciesName] end
function GkBGKCollisions:crossNormNu(speciesName) return self.normNuCross[speciesName] end

function GkBGKCollisions:calcCouplingMoments(tCurr, rkIdx, species)
   -- Compute self-primitive moments u and vtSq.
   local tmStart = Time.clock()

   local fIn      = species[self.speciesName]:rkStepperFields()[rkIdx]
   local momsSelf = species[self.speciesName]:fluidMoments()

   self.primMomsSelfCalc:advance(tCurr, {momsSelf, fIn}, {self.boundCorrs, self.primMomsSelf})
   self.timers.mom = self.timers.mom + Time.clock() - tmStart
end

function GkBGKCollisions:calcCrossCouplingMoments(tCurr, rkIdx, population)
   -- Perform cross-species calculation related to coupling moments that require the
   -- self-species coupling moments.
   local tmStart = Time.clock()

   -- Begin sending/receiving drift velocity and thermal speed squared if needed.
   population:speciesXferField_begin(self.primMomsSelfXfer, self.primMomsSelf, 33)

   self.timers.momcross = self.timers.momcross + Time.clock() - tmStart
end

function GkBGKCollisions:calcSelfNuTimeConst(momsSelf, nuOut) nuOut:copy(self.nuSelf) end

function GkBGKCollisions:calcSelfNuTimeDep(momsSelf, nuOut)
   -- Compute the Spitzer collisionality.
   self.m0Self:combineOffset(1., momsSelf, 0)
   self.vtSqSelf:combineOffset(1., self.primMomsSelf, self.confBasis:numBasis())
   self.spitzerNu:advance(tCurr, {self.charge, self.mass, self.m0Self, self.vtSqSelf,
                                  self.charge, self.mass, self.m0Self, self.vtSqSelf,
                                  self.normNuSelf, self.bmag}, {nuOut})
end

function GkBGKCollisions:calcCrossNuTimeConst(otherNm, chargeOther,
   mOther, momsOther, primMomsOther, nuCrossSelf, nuCrossOther) end

function GkBGKCollisions:calcCrossNuTimeDep(otherNm, chargeOther,
   mOther, momsOther, primMomsOther, nuCrossSelf, nuCrossOther)
   -- Compute the Spitzer collisionality. There's some duplication here in which
   -- two species compute the same cross collision frequency, but that is simpler
   -- and better for species parallelization.
   self.m0Other:combineOffset(1., momsOther, 0)
   self.vtSqOther:combineOffset(1., primMomsOther, self.confBasis:numBasis())

   local crossNormNuSelf = self.normNuCross[otherNm]
   self.spitzerNu:advance(tCurr, {self.charge, self.mass, self.m0Self, self.vtSqSelf,
                                  chargeOther, mOther, self.m0Other, self.vtSqOther, crossNormNuSelf, self.bmag},
                                 {nuCrossSelf})

   local crossNormNuOther = self.collAppOther[otherNm]:crossNormNu(self.speciesName)
   self.spitzerNu:advance(tCurr, {chargeOther, mOther, self.m0Other, self.vtSqOther,
                                  self.charge, self.mass, self.m0Self, self.vtSqSelf, crossNormNuOther, self.bmag},
                                 {nuCrossOther})
end

function GkBGKCollisions:advance(tCurr, fIn, population, out)
   local tmStart = Time.clock()

   local fRhsOut, cflRateByCell = out[1], out[2]
   local species = population:getSpecies()
   
   -- Fetch coupling moments of this species
   local momsSelf = species[self.speciesName]:fluidMoments()

   self.calcSelfNu(momsSelf, self.nuSum) -- Compute the collision frequency (if needed).

   -- Compute the Maxwellian.
   --self.maxwellian:advance(tCurr, {momsSelf, self.primMomsSelf, self.bmag, self.jacobTot}, {self.nufMaxwellSum})
   self.maxwellian:advance(tCurr, {momsSelf, self.bmag, self.jacobTot}, {self.nufMaxwellSum})
   if self.conservativeMaxwellian then
      -- Barrier before manipulations to moments before passing them to Iteration Fix     updater
      --Mpi.Barrier(self.phaseGrid:commSet().sharedComm)
      -- Call the CorrectMaxwellian updater.
      self.correctMaxwellian:advance(tCurr, {self.nufMaxwellSum, momsSelf}, {self.nufMaxwellSum})
   end
   self.phaseMul:advance(tCurr, {self.nuSum, self.nufMaxwellSum}, {self.nufMaxwellSum})
   
   if self.crossSpecies then
      for sInd, otherNm in ipairs(self.crossSpecies) do

         -- If we don't own the other species wait for its moments to be transferred.
         if not population:isSpeciesMine(otherNm) then
            population:speciesXferField_waitRecv(population:getSpecies()[otherNm].threeMomentsXfer)
            population:speciesXferField_waitRecv(self.collAppOther[otherNm].primMomsSelfXfer)
         end

         local mOther        = species[otherNm]:getMass()
         local momsOther     = species[otherNm]:fluidMoments()
         local primMomsOther = self.collAppOther[otherNm]:selfPrimitiveMoments()

         local nuCrossSelf  = self.nuCross[otherNm]
         local nuCrossOther = self.collAppOther[otherNm]:crossFrequencies(self.speciesName)

         -- Calculate time-dependent collision frequency if needed.
         local chargeOther = species[otherNm]:getCharge()
         self.calcCrossNu(otherNm, chargeOther, mOther, momsOther, primMomsOther,
                          nuCrossSelf, nuCrossOther)

         -- Compute cross moments.
         self.MomCrossBGK:advance(tCurr, {self.betaGreene, self.mass, momsSelf, mOther, momsOther, nuCrossSelf, nuCrossOther}, {self.momsCross})

	 self.maxwellian:advance(tCurr, {momsSelf, self.bmag, self.jacobTot}, {self.nufMaxwellCross})
         if self.conservativeMaxwellian then
            self.momsCrossSelf:combineOffset(1., self.momsCross, 0, 1., self.momsCross, self.confBasis:numBasis(), 1., self.momsCross, 2*self.confBasis:numBasis())
            self.correctMaxwellian:advance(tCurr, {self.nufMaxwellCross, self.momsCrossSelf}, {self.nufMaxwellCross})
         end
         self.phaseMul:advance(tCurr, {nuCrossSelf, self.nufMaxwellCross}, {self.nufMaxwellCross})
        
         self.nuSum:accumulate(1.0, nuCrossSelf)
         self.nufMaxwellSum:accumulate(1.0, self.nufMaxwellCross)
      end    -- end loop over other species that this species collides with
   end    -- end if self.crossCollisions
   
   self.collisionSlvr:advance(tCurr, {self.nuSum, self.nufMaxwellSum, fIn}, {fRhsOut, cflRateByCell})
   self.timers.advance = self.timers.advance + Time.clock() - tmStart
   print("Stepped forward!")
end

function GkBGKCollisions:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx)
   -- Wait to finish sending selfPrimMoms if needed.
   print("flag 0")
   population:speciesXferField_waitSend(self.primMomsSelfXfer)
   print("flag 1")
end

function GkBGKCollisions:write(tm, frame) end

return GkBGKCollisions
