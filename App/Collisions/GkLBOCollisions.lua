-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Gyrokinetic LB Collision operator.
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
local GkLBOconstNuEq = require "Eq.GkLBO"
local xsys           = require "xsys"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"
local lume           = require "Lib.lume"

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
function GkLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- previously stored table

   self.collKind = "GkLBO"    -- Type of collisions model. Useful at the species app level.

   self.collidingSpecies = assert(tbl.collideWith, "App.GkLBOCollisions: Must specify names of species to collide with in 'collideWith'.")

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
   
   -- Now establish if user wants constant or spatially varying collisionality.
   -- For constant nu, separate self and cross collision frequencies.
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
         self.normNuSelf  = self.normNuIn[selfSpecInd]
         if self.crossCollisions then
            self.normNuCross = lume.clone(self.normNuIn)
            table.remove(self.normNuCross, selfSpecInd)
         end
      else
         self.userInputNormNu = false
         self.normNuSelf       = 0.0
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
      self.charge     = speciesTbl.charge    -- Charge of this species.
      self.betaGreene = tbl.betaGreene and tbl.betaGreene or 0.0
   end

   self.nuFrac = tbl.nuFrac and tbl.nuFrac or 1.0

   self.mass = speciesTbl.mass   -- Mass of this species.

   self.cfl = 0.0    -- Will be replaced.

   self.timers = {nonSlvr = 0.}
end

function GkLBOCollisions:setName(nm) self.name = nm end
function GkLBOCollisions:setSpeciesName(nm) self.speciesName = nm end
function GkLBOCollisions:setCfl(cfl) self.cfl = cfl end
function GkLBOCollisions:setConfBasis(basis) self.confBasis = basis end
function GkLBOCollisions:setConfGrid(grid) self.confGrid = grid end
function GkLBOCollisions:setPhaseBasis(basis) self.phaseBasis = basis end
function GkLBOCollisions:setPhaseGrid(grid) self.phaseGrid = grid end

function GkLBOCollisions:createSolver(mySpecies, externalField)
   self.vDim = self.phaseGrid:ndim() - self.confGrid:ndim()

   -- Number of physical velocity dimensions.
   self.vDimPhys = 1.0
   if self.vDim == 2 then self.vDimPhys = 3.0 end

   -- Maximum velocity of the velocity grid (and its square).
   self.vParMax   = self.phaseGrid:upper(self.confGrid:ndim()+1)
   self.vParMaxSq = self.vParMax^2

   -- Collisionality, nu, summed over all species pairs.
   self.nuSum = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Parallel flow velocity times collisionality, summed over species.
   self.nuUParSum = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Thermal speed squared times collisionality, summed over species.
   self.nuVtSqSum = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Inverse of background magnetic field.
   self.bmag    = externalField.geo.bmag
   -- Inverse of background magnetic field.
   self.bmagInv = externalField.geo.bmagInv
      
   -- Zero-flux BCs in the velocity dimensions.
   local zfd = { }
   for d = 1, self.vDim do zfd[d] = self.confGrid:ndim() + d end

   -- Self-species collisionality, which varies in space.
   self.nuVarXSelf = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
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
   -- Weak multiplication to multiply nu(x) with uPar or vtSq.
   self.confMul = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "Multiply",
   }
   self.equation = {}
   self.collisionSlvr = Updater.GkLBO {
      onGrid     = self.phaseGrid,   confBasis = self.confBasis,
      phaseBasis = self.phaseBasis,  confRange = self.nuSum:localRange(),
      mass       = self.mass,
   }
   if self.crossCollisions then
      -- Temporary collisionality fields.
      self.nuCrossSelf = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      self.nuCrossOther = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Cross-species uPar and vtSq multiplied by collisionality.
      self.nuUParCross = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
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
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      self.m0s_deltas_den = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Weak division to compute the pre-factor in cross collision primitive moments.
      self.confDiv = Updater.CartFieldBinOp {
         weakBasis = self.confBasis,           operation = "Divide",
         onRange   = self.nuSum:localRange(),  onGhosts  = false,
      }
      -- Updater to compute cross-species primitive moments.
      self.primMomCross = Updater.CrossPrimMoments {
         onGrid     = self.confGrid,    betaGreene       = self.betaGreene,
         phaseBasis = self.phaseBasis,  varyingNu        = self.varNu,
         confBasis  = self.confBasis,   useCellAverageNu = self.cellConstNu,
         operator   = "GkLBO",
      }
   end

   -- Number of cells in which number density was negative (somewhere).
   self.primMomLimitCrossingsL = DataStruct.DynVector {
      numComponents = 1,
   }
   self.primMomLimitCrossingsG = DataStruct.DynVector {
      numComponents = 1,
   }
end

function GkLBOCollisions:advance(tCurr, fIn, species, out)

   local fRhsOut = out[1]
   local cflRateByCell = out[2]

   local tmNonSlvrStart = Time.clock()
   self.equation.primMomCrossLimit = 0.0

   -- Fetch coupling moments and primitive moments of this species.
   local momsSelf    = species[self.speciesName]:fluidMoments()
   local primMomSelf = species[self.speciesName]:selfPrimitiveMoments()

   if self.timeDepNu then
      -- Compute the Spitzer collisionality.
      self.spitzerNu:advance(tCurr, {self.charge, self.mass, momsSelf[1], primMomSelf[2],
                                     self.charge, self.mass, momsSelf[1], primMomSelf[2], 
                                     self.normNuSelf, self.bmag}, {self.nuSum})
   else
      self.nuSum:copy(self.nuVarXSelf)
   end
   self.confMul:advance(tCurr, {self.nuSum, primMomSelf[1]}, {self.nuUParSum})
   self.confMul:advance(tCurr, {self.nuSum, primMomSelf[2]}, {self.nuVtSqSum})

   if self.crossCollisions then

      local bCorrectionsSelf = species[self.speciesName]:boundaryCorrections()

      for sInd, otherNm in ipairs(self.crossSpecies) do

         local mOther            = species[otherNm]:getMass()
         local momsOther         = species[otherNm]:fluidMoments()
         local primMomOther      = species[otherNm]:selfPrimitiveMoments()
         local bCorrectionsOther = species[otherNm]:boundaryCorrections()

         if self.timeDepNu then
            -- Compute the collisionality if another species hasn't already done so.
            local chargeOther = species[otherNm]:getCharge()
            if (not species[self.speciesName].momentFlags[6][otherNm]) then
               self.spitzerNu:advance(tCurr, {self.charge, self.mass, momsSelf[1], primMomSelf[2],
                                              chargeOther, mOther, momsOther[1], primMomOther[2],
                                              self.normNuCross[sInd], self.bmag},
                                             {species[self.speciesName].nuVarXCross[otherNm]})
               species[self.speciesName].momentFlags[6][otherNm] = true
            end
            if (not species[otherNm].momentFlags[6][self.speciesName]) then
               self.spitzerNu:advance(tCurr, {chargeOther, mOther, momsOther[1], primMomOther[2],
                                              self.charge, self.mass, momsSelf[1], primMomSelf[2],
                                              species[otherNm].collPairs[otherNm][self.speciesName].normNu, self.bmag},
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

         -- Cross-primitive moments for the collision of these two species has not been computed.
         self.primMomCross:advance(tCurr, {self.mass, self.nuCrossSelf, momsSelf, primMomSelf, bCorrectionsSelf,
                                           mOther, self.nuCrossOther, momsOther, primMomOther, bCorrectionsOther,
                                           self.m0s_deltas},
                                          {species[self.speciesName].uParCross[otherNm], species[self.speciesName].vtSqCross[otherNm]})


         self.confMul:advance(tCurr, {self.nuCrossSelf, species[self.speciesName].uParCross[otherNm]}, {self.nuUParCross})
         self.confMul:advance(tCurr, {self.nuCrossSelf, species[self.speciesName].vtSqCross[otherNm]}, {self.nuVtSqCross})

         self.nuSum:accumulate(1.0, self.nuCrossSelf)
         self.nuUParSum:accumulate(1.0, self.nuUParCross)
         self.nuVtSqSum:accumulate(1.0, self.nuVtSqCross)

      end    -- end loop over other species that this species collides with.

   end    -- end if self.crossCollisions.
   self.timers.nonSlvr = self.timers.nonSlvr + Time.clock() - tmNonSlvrStart

   -- Compute increment from collisions and accumulate it into output.
   self.collisionSlvr:advance(
      tCurr, {fIn, self.bmagInv, self.nuUParSum, self.nuVtSqSum, self.nuSum}, {fRhsOut, cflRateByCell})
end

function GkLBOCollisions:write(tm, frame)
   Mpi.Allreduce(self.primMomLimitCrossingsL:data():data(), 
                 self.primMomLimitCrossingsG:data():data(), self.primMomLimitCrossingsG:size()*2,
                 Mpi.DOUBLE, Mpi.SUM, self.confGrid:commSet().comm)
   self.primMomLimitCrossingsG:write(string.format("%s_%s.bp", self.speciesName, "primMomLimitCrossings"), tm, frame)
   self.primMomLimitCrossingsL:clear(0.0)
   self.primMomLimitCrossingsG:clear(0.0)
end

function GkLBOCollisions:totalTime() return self.collisionSlvr.totalTime + self.timers.nonSlvr end

function GkLBOCollisions:slvrTime() return self.collisionSlvr.totalTime end

function GkLBOCollisions:nonSlvrTime() return self.timers.nonSlvr end

return GkLBOCollisions
