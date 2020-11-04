-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Charge exchange operator 
-- For details of model see:
-- H. L. Pauls, G. P. Zank, and L. L. Williams. Interaction of the solar wind
-- with the local interstellar medium. J. Geophys. Research,
-- 100:21,595 – 21,604, 1995
--
-- Meier, E. T. Modeling Plasmas with Strong Anisotropy, Neutral Fluid Effects,
-- and Open Boundaries. PhD thesis, U. Washington (2011).
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

-- GkChargeExchange  --------------------------------------------------------
--
-- Charge Exchange Operator
--------------------------------------------------------------------------------

local GkChargeExchange = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkChargeExchange:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkChargeExchange:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously store table.

   self.collKind = "CX"

   self.collidingSpecies = assert(tbl.collideWith, "App.GkChargeExchange: Must specify names of species to collide with in 'collideWith'.")

   self.collideNm   = tbl.collideWith[1]

   self.plasma      = assert(tbl.plasma, "App.GkChargeExchange: Must specify plasma species in 'plasma' ('H', 'D', or 'Ne')")
   
   self.ionNm       = tbl.ions
   self.neutNm      = tbl.neutrals
   self.iMass       = tbl.ionMass
   self.nMass       = tbl.neutMass
   self.charge      = tbl.charge

   -- Set these values to be consistent with other collision apps
   self.selfCollisions  = false
   self.crossCollisions = true              
   self.varNu           = false
   self.timeDepNu       = false
   self.collFreqs       = {1}

   if self.plasma=='H' then
      self.a = 1.12e-18
      self.b = 7.15e-20
   elseif self.plasma=='D' then 
      self.a = 1.09e-18
      self.b = 7.15e-20
   elseif self.plasma =='Ne' then
      self.a = 7.95e-19
      self.b = 5.65e-20
   end

   self._tmEvalMom = 0
end

function GkChargeExchange:setName(nm)
   self.name = nm
end

function GkChargeExchange:setSpeciesName(nm)
   self.speciesName = nm
end

function GkChargeExchange:setCfl(cfl)
   self.cfl = cfl
end

function GkChargeExchange:setConfBasis(basis)
   self.confBasis = basis
end

function GkChargeExchange:setConfGrid(grid)
   self.confGrid = grid
end

function GkChargeExchange:setPhaseBasis(basis)
   self.phaseBasis = basis
end

function GkChargeExchange:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function GkChargeExchange:createSolver(funcField)
   self.sourceCX =  DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.M0iDistFn =  DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.M0nDistFi =  DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.diffDistF =  DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.collisionSlvr = Updater.ChargeExchange {
      onGrid         = self.confGrid,
      confBasis      = self.confBasis,
      phaseBasis     = self.phaseBasis,
      charge         = self.charge, 
      a              = self.a,
      b              = self.b,
   }
   if (self.speciesName == self.ionNm) then --ions
      self.fMaxNeut =  DataStruct.Field {
	 onGrid        = self.phaseGrid,
	 numComponents = self.phaseBasis:numBasis(),
	 ghost         = {1, 1},
	 metaData = {
	    polyOrder = self.phaseBasis:polyOrder(),
	    basisType = self.phaseBasis:id()
	 },
      }
   else --neutrals
      self.fMaxIon =  DataStruct.Field {
	 onGrid        = self.phaseGrid,
	 numComponents = self.phaseBasis:numBasis(),
	 ghost         = {1, 1},
	 metaData = {
	    polyOrder = self.phaseBasis:polyOrder(),
	    basisType = self.phaseBasis:id()
	 },
      }
   end
end

function GkChargeExchange:advance(tCurr, fIn, species, fRhsOut)

   local writeOut = false
   -- Identify species and accumulate.
   if (self.speciesName == self.ionNm) then

      tmEvalMomStart = Time.clock()
      local neutM0   = species[self.neutNm]:fluidMoments()[1]
      local neutU    = species[self.neutNm]:selfPrimitiveMoments()[1] 
      local neutVtSq = species[self.neutNm]:selfPrimitiveMoments()[2]
      local ionM0    = species[self.ionNm]:fluidMoments()[1]
      local ionDistF = species[self.ionNm]:getDistF()

      species[self.speciesName].calcMaxwell:advance(tCurr,
         {neutM0, neutU, neutVtSq, species[self.speciesName].bmag}, {self.fMaxNeut})

      species[self.speciesName].confPhaseMult:advance(tCurr, {ionM0, self.fMaxNeut}, {self.M0iDistFn})
      species[self.speciesName].confPhaseMult:advance(tCurr, {neutM0, ionDistF}, {self.M0nDistFi})
      self.diffDistF:combine(1.0, self.M0iDistFn, -1.0, self.M0nDistFi)
      species[self.speciesName].confPhaseMult:advance(tCurr, {species[self.ionNm].vSigmaCX, self.diffDistF}, {self.sourceCX})

      if writeOut then
	 species[self.speciesName].distIo:write(self.fMaxNeut, string.format("%s_fMaxNeut_%d.bp",self.speciesName,tCurr*1e10),0,0,true)
	 species[self.speciesName].distIo:write(neutVtSq, string.format("%s_neutVtSq_%d.bp",self.speciesName,tCurr*1e10),0,0, true)
	 species[self.speciesName].distIo:write(self.sourceCX, string.format("%s_srcCX_%d.bp",self.speciesName,tCurr*1e10),0,0,true)
      end

      self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
      fRhsOut:accumulate(1.0,self.sourceCX)

   elseif (self.speciesName == self.neutNm) then

      tmEvalMomStart = Time.clock()      
      local ionM0     = species[self.ionNm]:fluidMoments()[1]
      local ionU      = species[self.ionNm]:selfPrimitiveMoments()[1] 
      local ionVtSq   = species[self.ionNm]:selfPrimitiveMoments()[2]
      local neutM0    = species[self.neutNm]:fluidMoments()[1]
      local neutDistF = species[self.neutNm]:getDistF()

      
      species[self.speciesName].calcMaxwell:advance(tCurr, {ionM0, ionU, ionVtSq}, {self.fMaxIon})

      species[self.speciesName].confPhaseMult:advance(tCurr, {ionM0, neutDistF}, {self.M0iDistFn})
      species[self.speciesName].confPhaseMult:advance(tCurr, {neutM0, self.fMaxIon}, {self.M0nDistFi})
      self.diffDistF:combine(1.0, self.M0iDistFn, -1.0, self.M0nDistFi)
      species[self.speciesName].confPhaseMult:advance(tCurr, {species[self.ionNm].vSigmaCX, self.diffDistF}, {self.sourceCX})
      
      if writeOut then
	 species[self.speciesName].distIo:write(neutDistF, string.format("%s_neutDistF_%d.bp",self.speciesName,tCurr*1e10),0,0)
	 species[self.speciesName].distIo:write(ionVtSq, string.format("%s_ionVtSq_%d.bp",self.speciesName,tCurr*1e10),0,0)
      end
      
      self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
      fRhsOut:accumulate(-self.iMass/self.nMass,self.sourceCX)

   end
   
end

function GkChargeExchange:write(tm, frame)
end

function GkChargeExchange:slvrTime()
   return 0
end

function GkChargeExchange:momTime()
   return self._tmEvalMom
end

function GkChargeExchange:projectMaxwellTime()
   return self.collisionSlvr:projectMaxwellTime()
end

return GkChargeExchange

