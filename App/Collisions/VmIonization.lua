-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Ionization operator using Voronov
-- reaction rate.
-- See:
-- Voronov, G.S., 1997. A practical fit formula for ionization rate
-- coefficients of atoms and ions by electron impact: Z= 1–28. Atomic
-- Data and Nuclear Data Tables, 65(1), pp.1-35.
--
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

-- VmIonization -----------------------------------------------------------
--
-- Voronov ionization operator.
---------------------------------------------------------------------------

local VmIonization = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmIonization:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmIonization:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously store table.

   self.cfl = 0.1
   self.collKind = "Ionization"

   self.collidingSpecies = assert(tbl.collideWith, "App.VmIonization: Must specify names of species to collide with in 'collideWith'.")

   -- Set these values to be consistent with other collision apps
   self.selfCollisions  = false
   self.crossCollisions = true              
   self.varNu           = false
   self.timeDepNu       = false
   self.collFreqs       = {1}
   
   self.collideNm   = tbl.collideWith[1]
   
   self.elcNm       = assert(tbl.electrons, "App.VmIonization: Must specify electron species name in 'electrons'.")
   self.neutNm      = assert(tbl.neutrals, "App.VmIonization: Must specify electron species name in 'neutrals'.")
   
   self.plasma      = tbl.plasma
   self.mass        = tbl.elcMass
   self.charge      = tbl.elemCharge

   if self.plasma == "H" then
      self._E = 13.6
      self._P = 0
      self._A = 0.291e-7
      self._K = 0.39
      self._X = 0.232
   end

   if self.plasma == "Ar" then
      self._E = 15.8
      self._P = 1
      self._A = 0.599e-7
      self._K = 0.26
      self._X = 0.136
   end

   self._tmEvalMom = 0.0 
end

function VmIonization:setName(nm)
   self.name = nm
end
function VmIonization:setSpeciesName(nm)
   self.speciesName = nm
end

function VmIonization:setCfl(cfl)
   self.cfl = cfl
end

function VmIonization:setConfBasis(basis)
   self.confBasis = basis
end
function VmIonization:setConfGrid(grid)
   self.confGrid = grid
end

function VmIonization:setPhaseBasis(basis)
   self.phaseBasis = basis
end

function VmIonization:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function VmIonization:createSolver(funcField)
   self.collisionSlvr = Updater.Ionization {
      onGrid     = self.confGrid,
      confBasis  = self.confBasis,
      phaseGrid     = self.phaseGrid,
      phaseBasis  = self.phaseBasis,
      elcMass    = self.mass,
      elemCharge = self.charge,
      reactRate  = true,
	 
      -- Voronov parameters
      A = self._A,
      E = self._E,
      K = self._K,
      P = self._P,
      X = self._X,
   }
   if (self.speciesName == self.elcNm) then
      self.calcIonizationTemp = Updater.Ionization {
	 onGrid     = self.confGrid,
	 confBasis  = self.confBasis,
	 phaseGrid  = self.phaseGrid,
	 phaseBasis = self.phaseBasis,
      	 elcMass    = self.mass,
      	 elemCharge = self.charge,
	 reactRate  = false, 
      	 E          = self._E,
      }
      self.sumDistF =  DataStruct.Field {
	 onGrid        = self.phaseGrid,
	 numComponents = self.phaseBasis:numBasis(),
	 ghost         = {1, 1},
	 metaData = {
	    polyOrder = self.phaseBasis:polyOrder(),
	    basisType = self.phaseBasis:id()
	 },
      }
   end
   self.confMult = Updater.CartFieldBinOp {
      onGrid     = self.confGrid,
      weakBasis  = self.confBasis,
      operation  = "Multiply",
   }
   self.confDiv = Updater.CartFieldBinOp {
      onGrid     = self.confGrid,
      weakBasis  = self.confBasis,
      operation  = "Divide",
   }
   self.confPhaseMult = Updater.CartFieldBinOp {
      onGrid     = self.phaseGrid,
      weakBasis  = self.phaseBasis,
      fieldBasis = self.confBasis,
      operation  = "Multiply",
   }
   self.m0elc = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.confBasis:polyOrder(),
	 basisType = self.confBasis:id()
      },
   }
   self.coefM0 = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.confBasis:polyOrder(),
	 basisType = self.confBasis:id()
      },
   }
   self.ionizSrc = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
end

function VmIonization:advance(tCurr, fIn, species, fRhsOut)
   local coefIz   = species[self.elcNm]:getVoronovReactRate()
   local elcM0    = species[self.elcNm]:fluidMoments()[1]
   local writeOut = false

   -- Check whether particle is electron, neutral or ion species.
   if (self.speciesName == self.elcNm) then
      -- Electrons.
      tmEvalMomStart   = Time.clock()
      local neutM0     = species[self.neutNm]:fluidMoments()[1]
      local elcDistF   = species[self.elcNm]:getDistF()
      local fMaxwellIz = species[self.elcNm]:getFMaxwellIz()
      
      self.sumDistF:combine(2.0,fMaxwellIz,-1.0,elcDistF)
 
      self.confMult:advance(tCurr, {coefIz, neutM0}, {self.coefM0})
      self.confPhaseMult:advance(tCurr, {self.coefM0, self.sumDistF}, {self.ionizSrc})
      -- Uncomment to test without fMaxwellian(Tiz).
      --self.confPhaseMult:advance(tCurr, {self.coefM0, elcDistF}, {self.ionizSrc})      
      if writeOut then
	 species[self.speciesName].distIo:write(self.ionizSrc, string.format("%s_izSrc_%d.bp",self.speciesName,tCurr*1e10),0,0,true)
      end

      self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
      fRhsOut:accumulate(1.0,self.ionizSrc)
   elseif (species[self.speciesName].charge == 0) then
      -- Neutrals.
      local neutDistF = species[self.neutNm]:getDistF()
      tmEvalMomStart  = Time.clock()
      self.m0elc:copy(elcM0)
            
      self.confMult:advance(tCurr, {coefIz, self.m0elc}, {self.coefM0})
      self.confPhaseMult:advance(tCurr, {self.coefM0, neutDistF}, {self.ionizSrc})
      if writeOut then
	 species[self.speciesName].distIo:write(self.ionizSrc, string.format("%s_izSrc_%d.bp",self.speciesName,tCurr*1e10),0,0,true)
      end

      self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
      fRhsOut:accumulate(-1.0,self.ionizSrc)  
   else
      -- Ions.
      tmEvalMomStart = Time.clock()
      self.m0elc:copy(elcM0)
      local neutDistF = species[self.neutNm]:getDistF()

      self.confMult:advance(tCurr, {coefIz, self.m0elc}, {self.coefM0})
      self.confPhaseMult:advance(tCurr, {self.coefM0, neutDistF}, {self.ionizSrc})

      if writeOut then
	 species[self.elcNm].distIo:write(neutVtSq, string.format("%s_neutVtSq_%d.bp",self.speciesName,tCurr*1e10),0,0)
	 species[self.speciesName].distIo:write(self.ionizSrc, string.format("%s_izSrc_%d.bp",self.speciesName,tCurr*1e10),0,0, true)
      end
      
      self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
      fRhsOut:accumulate(1.0,self.ionizSrc)
   end
end
   
function VmIonization:write(tm, frame)
end

function VmIonization:setCfl(cfl)
   self.cfl = cfl
end

function VmIonization:getIonizSrc()
   return self.ionizSrc
end

function VmIonization:slvrTime()
   local time = 0
   return time
end

function VmIonization:momTime()
   return self._tmEvalMom
end

function VmIonization:projectMaxwellTime()
   local time = self.confMult.projectMaxwellTime()
   time = time + self.confPhaseMult.projectMaxwellTime
   return time
end

return VmIonization

