-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Vlasov LB Collision operator
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
local Mpi            = require "Comm.Mpi"
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
function VmLBOCollisions:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.collKind = "VmLBO"  -- Type of collisions model.

   -- (MF 2022/11/18: not entirely true) For now only cell-wise constant nu is implemented.
   self.cellConstNu = true     -- Cell-wise constant nu?

   self.collidingSpecies = assert(tbl.collideWith, "App.VmLBOCollisions: Must specify names of species to collide with in 'collideWith'.")

   -- Determine if cross-species collisions take place,
   -- and put the names of the other colliding species in a list.
   local selfSpecInd = lume.find(self.collidingSpecies, self.speciesName)
   assert(selfSpecInd, "App.VmLBOCollisions: must include self-collisions.")

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
      self.calcSelfNu = function(momsIn, nuOut) VmLBOCollisions['calcSelfNuTimeConst'](self,momsIn,nuOut) end
      self.calcCrossNu = self.crossCollisions
         and function(otherNm, qOther, mOther, momsOther,
                      primMomsOther, vtSqMinOther, nuCrossSelf, nuCrossOther)
            VmLBOCollisions['calcCrossNuTimeConst'](self, otherNm, qOther,
              mOther, momsOther, primMomsOther, vtSqMinOther, nuCrossSelf, nuCrossOther)
         end
         or nil

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

      self.collFreqSelf = self.collFreqs[selfSpecInd]
      if self.crossCollisions then
         self.collFreqCross = lume.clone(self.collFreqs)
         table.remove(self.collFreqCross, selfSpecInd)
      end
   else
      -- Collisionality not provided by user. It will be calculated in time.
      self.timeDepNu = true
      self.calcSelfNu = function(momsIn, nuOut) VmLBOCollisions['calcSelfNuTimeDep'](self,momsIn,nuOut) end
      self.calcCrossNu = self.crossCollisions
         and function(otherNm, qOther, mOther, momsOther,
                      primMomsOther, vtSqMinOther, nuCrossSelf, nuCrossOther)
            VmLBOCollisions['calcCrossNuTimeDep'](self, otherNm, qOther,
              mOther, momsOther, primMomsOther, vtSqMinOther, nuCrossSelf, nuCrossOther)
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

   self.timers = {mom = 0.,   momcross = 0.,   advance = 0.,}
end

function VmLBOCollisions:setSpeciesName(nm)   self.speciesName = nm end
function VmLBOCollisions:setName(nm)          self.name        = self.speciesName.."_"..nm end
function VmLBOCollisions:setConfBasis(basis)  self.confBasis   = basis end
function VmLBOCollisions:setConfGrid(grid)    self.confGrid    = grid end
function VmLBOCollisions:setPhaseBasis(basis) self.phaseBasis  = basis end
function VmLBOCollisions:setPhaseGrid(grid)   self.phaseGrid   = grid end

function VmLBOCollisions:createSolver(mySpecies, extField)
   self.vdim = self.phaseGrid:ndim() - self.confGrid:ndim()

   -- Collisionality, nu, summed over all species pairs.
   self.nuSum = mySpecies:allocMoment()
   -- Sum of flow velocities in vdim directions and squared
   -- thermal speeds multiplied by respective collisionalities.
   self.nuPrimMomsSum = mySpecies:allocVectorMoment(self.vdim+1)

   -- Self-species collisionality, which varies in space.
   self.nuSelf = mySpecies:allocMoment()
   -- Allocate fields to store self-species primitive moments.
   self.primMomsSelf = mySpecies:allocVectorMoment(self.vdim+1)
   -- Allocate fields for boundary corrections.
   self.boundCorrs = mySpecies:allocVectorMoment(self.vdim+1)

   local vbounds = ffi.new("double[6]")
   for i = 1, self.vdim do
      vbounds[i-1]           = self.phaseGrid:lower(self.confGrid:ndim()+i)
      vbounds[i-1+self.vdim] = self.phaseGrid:upper(self.confGrid:ndim()+i)
   end
   self.primMomsSelfCalc = Updater.SelfPrimMoments {
      onGrid     = self.phaseGrid,   operator  = "VmLBO",
      phaseBasis = self.phaseBasis,  vbounds   = vbounds,
      confBasis  = self.confBasis,   confRange = self.nuSelf:localRange(),
   }

   local projectUserNu
   if self.timeDepNu then 
      self.m0Self      = mySpecies:allocMoment()  -- M0, to be extracted from fiveMoments.
      self.vtSqSelf    = mySpecies:allocMoment()
      self.vtSqMinSelf = mySpecies:vtSqMin()
      self.vtSqOther   = self.crossCollisions and mySpecies:allocMoment() or nil
      -- Updater to compute spatially varying (Spitzer) nu.
      self.spitzerNu = Updater.SpitzerCollisionality {
         confBasis       = self.confBasis,  epsilon0 = self.epsilon0,
         hBar            = self.hBar,       nuFrac   = self.nuFrac,
         willInputNormNu = self.normNuIn,
      }
      -- SpitzerCollisionality expects a magnetic field (set to 0 for Vlasov).
      self.bmag = (not self.normNuIn) and mySpecies:allocMoment() or nil
   else
      projectUserNu = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = self.collFreqSelf,
         basis  = self.confBasis,  onGhosts = false
      }
      projectUserNu:advance(0.0, {}, {self.nuSelf})
   end

   -- Weak multiplication to multiply nu(x) with u or vtSq.
   self.confMul = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "Multiply",
   }
   -- Lenard-Bernstein operator (LBO) collisions updater.
   -- Create table of pointers to fields needed in update
   self.fldPtrs = {self.nuSum, self.nuPrimMomsSum}
   self.collisionSlvr = Updater.VlasovLBO {
      onGrid     = self.phaseGrid,   confBasis = self.confBasis,
      phaseBasis = self.phaseBasis,  confRange = self.nuSelf:localRange(),
      fldPtrs    = self.fldPtrs
   }

   if self.crossCollisions then
      -- Cross-collision u and vtSq multiplied by collisionality.
      self.nuPrimMomsCross = mySpecies:allocVectorMoment(self.vdim+1)
      -- Prefactor m_0s*delta_s in cross primitive moment calculation.
      self.m0s_deltas = mySpecies:allocMoment()
      -- Updater to compute cross-species primitive moments.
      self.primMomCross = Updater.CrossPrimMoments {
         onGrid     = self.confGrid,    betaGreene = self.betaGreene, 
         phaseBasis = self.phaseBasis,  operator   = "VmLBO",
         confBasis  = self.confBasis,   m0s_deltas = self.m0s_deltas,
      }

      -- Allocate (and assign if needed) cross-species collision frequencies,
      -- and cross primitive moments.
      self.nuCross = {}
      for ispec, otherNm in ipairs(self.crossSpecies) do
         self.nuCross[otherNm] = mySpecies:allocMoment()
         if not self.timeDepNu then 
            projectUserNu:setFunc(self.collFreqCross[ispec])
            projectUserNu:advance(0.0, {}, {self.nuCross[otherNm]})
--            self.nuCross[otherNm]:write(string.format("%s_nu-%s_%d.bp",self.speciesName,otherNm,0),0.0,0)
         end
      end
      self.primMomsCross = mySpecies:allocVectorMoment(self.vdim+1)

      self.m0Other = self.timeDepNu and mySpecies:allocMoment() or nil  -- M0, to be extracted from fiveMoments.
   end
end


function VmLBOCollisions:createCouplingSolver(population, field, externalField)
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
         self.primMomsSelf = mySpecies:allocVectorMoment(mySpecies:getVdim()+1)

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
            -- Only species owned by this rank send primMoms to other ranks.
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

function VmLBOCollisions:boundaryCorrections() return self.boundCorrs end
function VmLBOCollisions:selfPrimitiveMoments() return self.primMomsSelf end
function VmLBOCollisions:crossFrequencies(speciesName) return self.nuCross[speciesName] end
function VmLBOCollisions:crossNormNu(speciesName) return self.normNuCross[speciesName] end

function VmLBOCollisions:createDiagnostics(mySpecies, field)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = VmLBODiagsImpl()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

function VmLBOCollisions:calcCouplingMoments(tCurr, rkIdx, species)
   -- Compute self-primitive moments u and vtSq.
   local tmStart = Time.clock()

   local fIn      = species[self.speciesName]:rkStepperFields()[rkIdx]
   local momsSelf = species[self.speciesName]:fluidMoments()

   self.primMomsSelfCalc:advance(tCurr, {momsSelf, fIn}, {self.boundCorrs, self.primMomsSelf})
   self.timers.mom = self.timers.mom + Time.clock() - tmStart
end

function VmLBOCollisions:calcCrossCouplingMoments(tCurr, rkIdx, population)
   -- Perform cross-species calculation related to coupling moments that require the
   -- self-species coupling moments.
   local tmStart = Time.clock()

   -- Begin sending/receiving drift velocity and thermal speed squared if needed.
   population:speciesXferField_begin(self.primMomsSelfXfer, self.primMomsSelf, 33)

   self.timers.momcross = self.timers.momcross + Time.clock() - tmStart
end

function VmLBOCollisions:calcSelfNuTimeConst(momsSelf, nuOut) nuOut:copy(self.nuSelf) end

function VmLBOCollisions:calcSelfNuTimeDep(momsSelf, nuOut)
   -- Compute the Spitzer collisionality.
   self.m0Self:combineOffset(1., momsSelf, 0) 
   self.vtSqSelf:combineOffset(1., self.primMomsSelf, self.vdim*self.confBasis:numBasis()) 
   self.spitzerNu:advance(tCurr, {self.charge, self.mass, self.m0Self, self.vtSqSelf, self.vtSqMinSelf,
                                  self.charge, self.mass, self.m0Self, self.vtSqSelf, self.vtSqMinSelf,
				  self.normNuSelf, self.bmag}, {nuOut})
end

function VmLBOCollisions:calcCrossNuTimeConst(otherNm, chargeOther,
   mOther, momsOther, primMomsOther, vtSqMinOther, nuCrossSelf, nuCrossOther) end

function VmLBOCollisions:calcCrossNuTimeDep(otherNm, chargeOther,
   mOther, momsOther, primMomsOther, vtSqMinOther, nuCrossSelf, nuCrossOther)
   -- Compute the Spitzer collisionality. There's some duplication here in which
   -- two species compute the same cross collision frequency, but that is simpler
   -- and better for species parallelization.
   self.m0Other:combineOffset(1., momsOther, 0) 
   self.vtSqOther:combineOffset(1., primMomsOther, self.vdim*self.confBasis:numBasis()) 

   local crossNormNuSelf = self.normNuCross[otherNm]
   self.spitzerNu:advance(tCurr, {self.charge, self.mass, self.m0Self, self.vtSqSelf, self.vtSqMinSelf,
                                  chargeOther, mOther, self.m0Other, self.vtSqOther, vtSqMinOther, crossNormNuSelf, self.bmag},
                                 {nuCrossSelf})

   local crossNormNuOther = self.collAppOther[otherNm]:crossNormNu(self.speciesName)
   self.spitzerNu:advance(tCurr, {chargeOther, mOther, self.m0Other, self.vtSqOther, vtSqMinOther,
                                  self.charge, self.mass, self.m0Self, self.vtSqSelf, self.vtSqMinSelf, crossNormNuOther, self.bmag},
                                 {nuCrossOther})
end

function VmLBOCollisions:advance(tCurr, fIn, population, out)
   local tmStart = Time.clock()

   local fRhsOut, cflRateByCell = out[1], out[2]
   local species = population:getSpecies()

   -- Fetch coupling moments of this species.
   local momsSelf = species[self.speciesName]:fluidMoments()

   self.calcSelfNu(momsSelf, self.nuSum)
   self.confMul:advance(tCurr, {self.nuSum, self.primMomsSelf}, {self.nuPrimMomsSum})

   if self.crossCollisions then

      local bCorrectionsSelf = self.boundCorrs

      for _, otherNm in ipairs(self.crossSpecies) do

         -- If we don't own the other species wait for its moments to be transferred. 
         if not population:isSpeciesMine(otherNm) then
            population:speciesXferField_waitRecv(population:getSpecies()[otherNm].fiveMomentsXfer)
            population:speciesXferField_waitRecv(self.collAppOther[otherNm].primMomsSelfXfer)
         end

         local mOther        = species[otherNm]:getMass()
         local momsOther     = species[otherNm]:fluidMoments()
         local primMomsOther = self.collAppOther[otherNm]:selfPrimitiveMoments()
         local vtSqMinOther  = species[otherNm]:vtSqMin()

         local nuCrossSelf  = self.nuCross[otherNm]
         local nuCrossOther = self.collAppOther[otherNm]:crossFrequencies(self.speciesName)

         -- Calculate time-dependent collision frequency if needed.
         local chargeOther = species[otherNm]:getCharge()
         self.calcCrossNu(otherNm, chargeOther, mOther, momsOther, primMomsOther,
                          vtSqMinOther, nuCrossSelf, nuCrossOther)

         -- Calculate cross primitive moments.
         self.primMomCross:advance(tCurr, {self.mass, nuCrossSelf, momsSelf, self.primMomsSelf, bCorrectionsSelf,
                                           mOther, nuCrossOther, momsOther, primMomsOther, momsSelf},
                                          {self.primMomsCross})

         self.confMul:advance(tCurr, {nuCrossSelf, self.primMomsCross}, {self.nuPrimMomsCross})

         self.nuSum:accumulate(1.0, nuCrossSelf)
         self.nuPrimMomsSum:accumulate(1.0, self.nuPrimMomsCross)

      end  -- end loop over other species that this species collides with.

   end  -- end if self.crossCollisions.

   -- Compute increment from collisions and accumulate it into output.
   self.collisionSlvr:advance(tCurr, {fIn, self.nuPrimMomsSum, self.nuSum}, {fRhsOut, cflRateByCell})

   self.timers.advance = self.timers.advance + Time.clock() - tmStart
end

function VmLBOCollisions:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx)
   -- Wait to finish sending selfPrimMoms if needed.
   population:speciesXferField_waitSend(self.primMomsSelfXfer)
end

function VmLBOCollisions:write(tm, frame) end

return VmLBOCollisions
