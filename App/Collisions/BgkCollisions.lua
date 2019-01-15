-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: BGK Collision operator.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local DataStruct     = require "DataStruct"
local Proto          = require "Lib.Proto"
local Updater        = require "Updater"
local xsys           = require "xsys"

-- BgkCollisions ---------------------------------------------------------------
--
-- Bhatnagar-Gross-Krook Collision operator.
--------------------------------------------------------------------------------

local BgkCollisions = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function BgkCollisions:init(tbl)
   self.tbl = tbl
end

-- Function to find the index of an element in table.
local function findInd(tbl, el)
   for i, v in ipairs(tbl) do
      if v == el then
         return i
      end
   end
   return #tbl+1    -- If not found return a number larger than the length of the table.
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function BgkCollisions:fullInit(speciesTbl)
   local tbl = self.tbl    -- Previously stored table.

   self.cfl = 0.1    -- Some default value.

   local collidingSpecies = assert(tbl.collideWith, "App.VmLBOCollisions: Must specify names of species to collide with in 'collideWith'.")

   -- First determine if self-species and/or cross-species collisions take place,
   -- and (if cross-collisions=true) put the names of the other colliding species in a list.
   local selfSpecInd = findInd(collidingSpecies, self.speciesName)
   if selfSpecInd < (#collidingSpecies+1) then
      self.selfCollisions = true                 -- Apply self-species collisions.
      if #collidingSpecies > 1 then
         self.crossCollisions = true             -- Apply cross-species collisions.
         self.crossSpecies    = collidingSpecies
         table.remove(self.crossSpecies, selfSpecInd)
      else
         self.crossCollisions = false            -- Don't apply cross-species collisions.
      end
   else
      self.selfCollisions  = false               -- Don't apply self-species collisions.
      self.crossCollisions = true                -- Apply cross-species collisions.
      self.crossSpecies    = collidingSpecies    -- All species in collidingSpecies must be cross-species.
   end

   -- Now establish if user wants constant or spatially varying collisionality.
   -- For constant nu, separate self and cross collision frequencies.
   local collFreqs          = tbl.frequencies -- List of collision frequencies, if using spatially constant nu.
   if collFreqs then
      self.varNu            = false    -- Not spatially varying nu.
      self.cellConstNu      = true     -- Cell-wise constant nu?
      if self.selfCollisions then
         self.collFreqSelf  = collFreqs[selfSpecInd]
      end
      if self.crossCollisions then
         self.collFreqCross = collFreqs
         table.remove(self.collFreqCross, selfSpecInd)
      end
   else
      self.varNu       = true    -- Spatially varying nu.
      self.mass        = speciesTbl.mass      -- Mass of this species.
      self.charge      = speciesTbl.charge    -- Charge of this species.
      -- For now only cell-wise constant nu is implemented.
      -- self.cellConstNu = assert(tbl.cellAvFrequencies, "App.GkLBOCollisions: Must specify 'useCellAverageNu=true/false' for using cellwise constant/expanded spatially varying collisionality.")
      self.cellConstNu = true
      -- If no constant collision frequencies provided ('frequencies'), user can specify 'normNu'
      -- list of collisionalities normalized by (T_0^(3/2)/n_0) evaluated somewhere in the
      -- simulation. Otherwise code compute Spitzer collisionality from scratch.
      local normNuIn   = tbl.normNu
      -- normNuSelf, epsilon0 and elemCharge may not used, but are
      -- initialized below to avoid if-statements in advance method.
      if normNuIn then
         self.userInputNormNu = true
         if self.selfCollisions then
            self.normNuSelf  = normNuIn[selfSpecInd]
         end
         if self.crossCollisions then
            self.normNuCross = normNuIn
            table.remove(self.normNuCross, selfSpecInd)
         end
         self.epsilon0   = 8.854187817620389850536563031710750260608e-12    -- Farad/meter.
         self.elemCharge = 1.602176487e-19    -- Coulomb.
      else
         self.userInputNormNu = false
         if self.selfCollisions then
            self.normNuSelf  = 0.0
         end
         if self.crossCollisions then
            self.normNuCross = collidingSpecies
            table.remove(self.normNuCross, selfSpecInd)
            for i, _ in ipairs(self.normNuCross) do self.normNuCross[i] = 0.0 end
         end
         self.epsilon0 = assert(
            tbl.epsilon0, "Updater.VmLBOCollisions: Must specify vacuum permittivity 'epsilon0' ('frequencies' and 'normNu' are not specified, so nu is calculated via Spitzer).")
         self.elemCharge = assert(
            tbl.elemCharge, "Updater.VmLBOCollisions: Must specify elementary charge with 'elemCharge' ('frequencies' and 'normNu' are not specified, so nu is calculated via Spitzer).")
--         self.coulombLog = assert(
--      	    tbl.coulombLog, "Updater.BgkCollisions: Must specify Coulomb logaritm with 'coulombLog' ('frequencies' and 'normNu' are not specified, so nu is calculated via Spitzer).")
      end
   end

   if self.crossCollisions then
      self.mass       = speciesTbl.mass      -- Mass of this species.
      self.charge     = speciesTbl.charge    -- Charge of this species.
      -- Currently only crossOpIn=Greene is available for BGK.
--      local crossOpIn = tbl.crossOption    -- Can specify 'crossOption' (Greene, GreeneSmallAngle, HeavyIons), formulas used to calculate cross-species primitive moments.
      crossOpIn = "Greene"
      if crossOpIn then
         self.crossMomOp  = crossOpIn
         if self.crossMomOp=="Greene" then
            local betaGreeneIn = tbl.betaGreene   -- Can specify 'betaGreene' free parameter in Grene cross-species collisions.
            if betaGreeneIn then
               self.betaGreene = betaGreeneIn
            else
               self.betaGreene = 1.0   -- Default value is the heavy ion, quasineutral limit.
            end
         else
            self.betaGreene = 0.0   -- Default value is the heavy ion, quasineutral limit.
         end
      else
         self.crossMomOp  = "Greene"    -- Default to Greene-type formulas.
         self.betaGreene  = 1.0         -- Default value is the heavy ion, quasineutral limit.
      end
   end

   if tbl.nuFrac then
      self.nuFrac = tbl.nuFrac
   else
      self.nuFrac = 1.0
   end

   self.exactLagFixM012 = xsys.pickBool(tbl.exactLagFixM012, true) 
end

function BgkCollisions:setName(nm)
   self.name = nm
end
function BgkCollisions:setSpeciesName(nm)
   self.speciesName = nm
end

function BgkCollisions:setCfl(cfl)
   self.cfl = cfl
end
function BgkCollisions:setConfBasis(basis)
   self.confBasis = basis
end
function BgkCollisions:setConfGrid(grid)
   self.confGrid = grid
end
function BgkCollisions:setPhaseBasis(basis)
   self.phaseBasis = basis
end
function BgkCollisions:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function BgkCollisions:createSolver()
   self.numVelDims = self.phaseGrid:ndim() - self.confGrid:ndim()

   local function createConfFieldCompV()
      return DataStruct.Field {
	 onGrid        = self.confGrid,
	 numComponents = self.confBasis:numBasis()*self.numVelDims,
	 ghost         = {1, 1},
      }
   end
   local function createConfFieldComp1()
      return DataStruct.Field {
	 onGrid        = self.confGrid,
	 numComponents = self.confBasis:numBasis(),
	 ghost         = {1, 1},
      }
   end

   self.uSelf     = createConfFieldCompV()
   self.dM1       = createConfFieldCompV()
   self.vthSqSelf = createConfFieldComp1()
   self.dM2       = createConfFieldComp1()
   self.dM0       = createConfFieldComp1()
   if self.crossCollisions then
      self.uOther     = createConfFieldCompV()
      self.uCross     = createConfFieldCompV()
      self.vthSqOther = createConfFieldComp1()
      self.vthSqCross = createConfFieldComp1()
   end

   if self.varNu then
      -- Collisionality, nu.
      self.collFreq = createConfFieldComp1()
      -- Updater to compute spatially varying (Spitzer) nu.
      self.spitzerNu = Updater.SpitzerCollisionality {
         onGrid           = self.confGrid,
         confBasis        = self.confBasis,
         useCellAverageNu = self.cellConstNu,
         willInputNormNu  = self.userInputNormNu,
         elemCharge       = self.elemCharge,
         epsilon0         = self.epsilon0,
      }
   else
      self.collFreq = 0.0    -- Assigned in advance method.
   end

   -- Dummy fields for the primitive moment calculator.
   self._kinEnergyDens = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }
   self._thermEnergyDens = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Updaters for the primitive moments.
   self.confDiv = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
   }
   self.confDotProduct = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "DotProduct",
   }

   -- Maxwellian field.
   self.fMaxwell = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Maxwellian solver.
   self.maxwellian = Updater.MaxwellianOnBasis {
      onGrid     = self.phaseGrid,
      confGrid   = self.confGrid,
      confBasis  = self.confBasis,
      phaseGrid  = self.phaseGrid,
      phaseBasis = self.phaseBasis,
   }

   -- BGK Collision solver itself.
   self.collisionSlvr = Updater.BgkCollisions {
      onGrid           = self.phaseGrid,
      confGrid         = self.confGrid,
      confBasis        = self.confBasis,
      phaseBasis       = self.phaseBasis,
      varyingNu        = self.varNu,
      useCellAverageNu = self.cellConstNu,
   }

   self.m0Calc = Updater.DistFuncMomentCalc {
      onGrid     = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confBasis  = self.confBasis,
      moment     = "M0",
   }
   self.m1Calc = Updater.DistFuncMomentCalc {
      onGrid     = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confBasis  = self.confBasis,
      moment     = "M1i",
   }
   self.m2Calc = Updater.DistFuncMomentCalc {
      onGrid     = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confBasis  = self.confBasis,
      moment     = "M2",
   }

   self.lagFix = Updater.LagrangeFix {
      onGrid     = self.phaseGrid,
      phaseGrid  = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confGrid   = self.confGrid,
      confBasis  = self.confBasis,
      mode       = 'Vlasov',
}
end

function BgkCollisions:primMoments(M0, M1i, M2, u, vthSq)
   self.confDiv:advance(0., {M0, M1i}, {u})
   self.confDotProduct:advance(0., {u, M1i},
			       {self._kinEnergyDens})
   self._thermEnergyDens:combine( 1.0/self.numVelDims, M2,
                                 -1.0/self.numVelDims, self._kinEnergyDens)
   self.confDiv:advance(0., {M0, self._thermEnergyDens}, {vthSq})
end

function BgkCollisions:advance(tCurr, fIn, species, fRhsOut)
   local selfMom = species[self.speciesName]:fluidMoments()

   self:primMoments(selfMom[1], selfMom[2], selfMom[3], 
		    self.uSelf, self.vthSqSelf)

   if self.selfCollisions then
      self.maxwellian:advance(
	 tCurr, {selfMom[1], self.uSelf, self.vthSqSelf}, {self.fMaxwell})
      if self.exactLagFixM012 then
	 self.m0Calc:advance(0.0, {self.fMaxwell}, {self.dM0})
	 self.dM0:scale(-1)
	 self.dM0:accumulate(1, selfMom[1])
	 self.m1Calc:advance(0.0, {self.fMaxwell}, {self.dM1})
	 self.dM1:scale(-1)
	 self.dM1:accumulate(1, selfMom[2])
	 self.m2Calc:advance(0.0, {self.fMaxwell}, {self.dM2})
	 self.dM2:scale(-1)
	 self.dM2:accumulate(1, selfMom[3])
	 self.lagFix:advance(0.0, {self.dM0, self.dM1, self.dM2}, {self.fMaxwell})
      end

      if self.varNu then
         -- Compute the collisionality.
         self.spitzerNu:advance(0.0, {self.mass, self.charge, selfMom[1], self.vthSqSelf, self.normNuSelf},{self.collFreq})
      else
         self.collFreq = self.collFreqSelf
      end

      self.collisionSlvr:advance(
	 tCurr, {fIn, self.fMaxwell, self.collFreq, self.nuFrac}, {fRhsOut})
   end

   if self.crossSpecies then
      for sInd, otherNm in ipairs(self.crossSpecies) do
	 local mOther   = species[otherNm]:getMass()
	 local otherMom = species[otherNm]:fluidMoments()
	 self:primMoments(otherMom[1], otherMom[2], otherMom[3], 
			  self.uOther, self.vthSqOther)
	 
	 -- Greene, Improved BGK model for electron-ion collisions, 1973.
	 self.uCross:combine(1.0, self.uSelf, -1.0, self.uOther)
	 self.confDotProduct:advance(0., {self.uCross, self.uCross}, {self._kinEnergyDens})
	 self.uCross:combine( 0.5, self.uSelf, 0.5, self.uOther,
                             -0.5*self.betaGreene, self.uSelf,
                              0.5*self.betaGreene, self.uOther)
	 self.vthSqCross:combine(mOther, self.vthSqSelf, mOther, self.vthSqOther,
                                -self.betaGreene*self.mass, self.vthSqSelf,
                                 self.betaGreene*mOther, self.vthSqOther,
                                 (1-self.betaGreene*self.betaGreene)/6*mOther, self._kinEnergyDens, 
                                 (1+self.betaGreene)*(1+self.betaGreene)/12*(mOther-self.mass), self._kinEnergyDens)
	 self.vthSqCross:scale(1.0/(mOther+self.mass))

         if self.varNu then
            -- Compute the collisionality.
            self.spitzerNu:advance(0., {self.mass, self.charge, otherMom[1], self.vthSqSelf, self.normNuCross[sInd]}, {self.collFreq})
         else
            self.collFreq = self.collFreqCross[sInd]
         end

	 self.maxwellian:advance(tCurr, {selfMom[1], self.uCross, self.vthSqCross}, {self.fMaxwell})
	 self.collisionSlvr:advance(tCurr, {fIn, self.fMaxwell, self.collFreq}, {fRhsOut})
      end
   end
end

function BgkCollisions:write(tm, frame)
end

function BgkCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.maxwellian.totalTime
end

return BgkCollisions
