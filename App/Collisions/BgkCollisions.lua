-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: BGK Collision operator
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local DataStruct = require "DataStruct"
local Proto = require "Lib.Proto"
local Updater = require "Updater"
local xsys = require "xsys"

-- BgkCollisions ---------------------------------------------------------------
--
-- Bhatnagar-Gross-Krook Collision operator
--------------------------------------------------------------------------------

local BgkCollisions = Proto(CollisionsBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function BgkCollisions:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function BgkCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- previously store table

   self.cfl = 0.1 -- some default value
   self.selfCollisions = xsys.pickBool(tbl.selfCollisions, true) -- by default, self collisions are on
   self.crossSpecies = tbl.crossSpecies
   if self.crossSpecies then
      self.beta = assert(
	 tbl.beta, "Updater.BgkCollisions: Must specify the beta from Green with 'beta' (cross-species are ON)")
   end

   self.collFreq = tbl.collFreq
   if tbl.nuFrac then
      self.nuFrac =  tbl.nuFrac
   else
      self.nuFrac = 1.0
   end
   if not self.collFreq then 
      self.mass = speciesTbl.mass
      self.charge = speciesTbl.charge
      self.epsilon0 = assert(
	 tbl.epsilon0, "Updater.BgkCollisions: Must specify vacuum permitivity with 'epsilon0' ('collFreq' is not specified, so classical nu is used instead)")
      self.elemCharge = assert(
	 tbl.elemCharge, "Updater.BgkCollisions: Must specify elementary charge with 'elemCharge' ('collFreq' is not specified, so classical nu is used instead)")
      -- self.coulombLog = assert(
      -- 	 tbl.coulombLog, "Updater.BgkCollisions: Must specify Coulomb logaritm with 'coulombLog' ('collFreq' is not specified, so classical nu is used instead)")
   end

   self.exactLagFixM012 = xsys.pickBool(tbl.exactLagFixM012, true) 

   self.tmEvalMom = 0.0
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
	 onGrid = self.confGrid,
	 numComponents = self.confBasis:numBasis()*self.numVelDims,
	 ghost = {1, 1},
      }
   end
   self.uSelf = createConfFieldCompV()
   self.dM1 = createConfFieldCompV()
   if self.selfCollisions then
      self.uOther = createConfFieldCompV()
      self.uCross = createConfFieldCompV()
   end

   local function createConfFieldComp1()
      return DataStruct.Field {
	 onGrid = self.confGrid,
	 numComponents = self.confBasis:numBasis(),
	 ghost = {1, 1},
      }
   end
   self.vth2Self = createConfFieldComp1()
   self.dM2 = createConfFieldComp1()
   self.dM0 = createConfFieldComp1()
   if self.selfCollisions then
      self.vth2Other = createConfFieldComp1()
      self.vth2Cross = createConfFieldComp1()
   end

   -- dummy fields for the primitive moment calculator
   self._kinEnergyDens = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost = {1, 1},
   }
   self._thermEnergyDens = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost = {1, 1},
   }
   -- Updaters for the primitive moments
   self.confDiv = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
   }
   self.confDotProduct = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "DotProduct",
   }

   -- Maxwellian field
   self.fMaxwell = DataStruct.Field {
      onGrid = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost = {1, 1},
   }
   -- Maxwellian solver
   self.maxwellian = Updater.MaxwellianOnBasis {
      onGrid = self.phaseGrid,
      confGrid = self.confGrid,
      confBasis = self.confBasis,
      phaseGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
   }

   -- BGK Collision solver itself
   self.collisionSlvr = Updater.BgkCollisions {
      onGrid = self.phaseGrid,
      confGrid = self.confGrid,
      confBasis = self.confBasis,
      phaseGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      collFreq = self.collFreq,
      mass = self.mass,
      charge = self.charge,
      elemCharge = self.elemCharge,
      epsilon0 = self.epsilon0,
      --coulombLog = self.coulombLog,
   }

   self.m0Calc = Updater.DistFuncMomentCalc {
      onGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confBasis = self.confBasis,
      moment = "M0",
   }
   self.m1Calc = Updater.DistFuncMomentCalc {
      onGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confBasis = self.confBasis,
      moment = "M1i",
   }
   self.m2Calc = Updater.DistFuncMomentCalc {
      onGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confBasis = self.confBasis,
      moment = "M2",
   }

   self.lagFix = Updater.LagrangeFix {
      onGrid = self.phaseGrid,
      phaseGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confGrid = self.confGrid,
      confBasis = self.confBasis,
      mode = 'Vlasov',
}
end

function BgkCollisions:primMoments(M0, M1i, M2, u, vth2)
   self.confDiv:advance(0., {M0, M1i}, {u})
   self.confDotProduct:advance(0., {u, M1i},
			       {self._kinEnergyDens})
   self._thermEnergyDens:combine(1.0/self.numVelDims, M2,
				   -1.0/self.numVelDims, self._kinEnergyDens)
   self.confDiv:advance(0., {M0, self._thermEnergyDens}, {vth2})
end

function BgkCollisions:advance(tCurr, fIn, species, fRhsOut)
   local selfMom = species[self.speciesName]:fluidMoments()

   local tmEvalMomStart = Time.clock()
   self:primMoments(selfMom[1], selfMom[2], selfMom[3], 
		    self.uSelf, self.vth2Self)
   self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

   if self.selfCollisions then
      self.maxwellian:advance(
	 tCurr, {selfMom[1], self.uSelf, self.vth2Self},
	 {self.fMaxwell})
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
	 self.lagFix:advance(0.0, {self.dM0, self.dM1, self.dM2},
			     {self.fMaxwell})
      end

      self.collisionSlvr:advance(
	 tCurr, {fIn, self.fMaxwell, self.nuFrac, selfMom[1], self.vth2Self}, {fRhsOut})
   end

   if self.crossSpecies then
      for _, otherNm in ipairs(self.crossSpecies) do
	 local mSelf = species[self.speciesName]:getMass()
	 local mOther = species[otherNm]:getMass()
	 local otherMom = species[otherNm]:fluidMoments()
         local tmEvalMomStart = Time.clock()
	 self:primMoments(otherMom[1], otherMom[2], otherMom[3], 
			  self.uOther, self.vth2Other)
	 
	 -- Greene, Improved BGK model for electron-ion collisions, 1973
	 self.uCross:combine(1.0, self.uSelf, -1.0, self.uOther)
	 self.confDotProduct:advance(0., {self.uCross, self.uCross},
				     {self._kinEnergyDens})
	 self.uCross:combine(0.5, self.uSelf, 0.5, self.uOther,
				-0.5*self.beta, self.uSelf,
			     0.5*self.beta, self.uOther)
	 self.vth2Cross:combine(mOther, self.vth2Self, mOther, self.vth2Other,
				   -self.beta*mSelf, self.vth2Self,
				self.beta*mOther, self.vth2Other,
				(1-self.beta*self.beta)/6*mOther,
				self._kinEnergyDens, 
				(1+self.beta)*(1+self.beta)/12*(mOther-mSelf),
				self._kinEnergyDens)
	 self.vth2Cross:scale(1.0/(mOther+mSelf))
         self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

	 self.maxwellian:advance(tCurr, {selfMom[1], self.uCross, self.vth2Cross}, {self.fMaxwell})
	 self.collisionSlvr:advance(tCurr, {fIn, self.fMaxwell}, {fRhsOut})
      end
   end
end

function BgkCollisions:write(tm, frame)
end

function BgkCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.maxwellian.totalTime + self.tmEvalMom
end

function BgkCollisions:slvrTime()
   return self.collisionSlvr.totalTime + self.maxwellian.totalTime
end

function BgkCollisions:momTime()
   return self.tmEvalMom
end

return BgkCollisions
