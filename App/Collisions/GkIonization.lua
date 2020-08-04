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

-- GkIonization -----------------------------------------------------------
--
-- Voronov ionization operator.
--------------------------------------------------------------------------------

local GkIonization = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkIonization:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkIonization:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously store table.

   self.collKind = "Ionization"

   self.collidingSpecies = assert(tbl.collideWith, "App.GkIonization: Must specify names of species to collide with in 'collideWith'.")

   -- Set these values to be consistent with other collision apps
   self.selfCollisions  = false
   self.crossCollisions = true              
   self.varNu           = false
   self.timeDepNu       = false
   self.collFreqs       = {1}
   
   self.collideNm   = tbl.collideWith[1]
   
   self.elcNm       = assert(tbl.electrons, "App.GkIonization: Must specify electron species name in 'electrons'.")
   self.neutNm      = assert(tbl.neutrals, "App.GkIonization: Must specify electron species name in 'neutrals'.")
   
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

function GkIonization:setName(nm)
   self.name = nm
end
function GkIonization:setSpeciesName(nm)
   self.speciesName = nm
end

function GkIonization:setCfl(cfl)
   self.cfl = cfl
end

function GkIonization:setConfBasis(basis)
   self.confBasis = basis
end
function GkIonization:setConfGrid(grid)
   self.confGrid = grid
end

function GkIonization:setPhaseBasis(basis)
   self.phaseBasis = basis
end

function GkIonization:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function GkIonization:createSolver(funcField)
   if (self.speciesName == self.elcNm) then
      self.calcVoronovReactRate = Updater.Ionization {
	 onGrid     = self.confGrid,
	 confBasis  = self.confBasis,
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
      self.calcIonizationTemp = Updater.Ionization {
      	 onGrid     = self.confGrid,
      	 confBasis  = self.confBasis,
      	 elcMass    = self.mass,
      	 elemCharge = self.charge,
	 reactRate  = false, 
      	 E          = self._E,
      }
      self.sumDistF    = DataStruct.Field {
      	 onGrid        = self.phaseGrid,
      	 numComponents = self.phaseBasis:numBasis(),
      	 ghost         = {1, 1},
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
      operation = "Divide",
   }
   self.collisionSlvr = Updater.CartFieldBinOp {
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
   self.neutDistF = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.fNeutMax = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
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

function GkIonization:advance(tCurr, fIn, species, fRhsOut)
   local coefIz = species[self.elcNm]:getVoronovReactRate()
   local elcM0 = species[self.elcNm]:fluidMoments()[1]

   -- Check whether particle is electron, neutral or ion species
   if (self.speciesName == self.elcNm) then
      -- electrons
      tmEvalMomStart   = Time.clock()
      local neutM0     = species[self.neutNm]:fluidMoments()[1]
      local elcDistF   = species[self.elcNm]:getDistF()
      local fMaxwellIz = species[self.elcNm]:getFMaxwellIz()
      
      self.sumDistF:combine(2.0,fMaxwellIz,-1.0,elcDistF)
 
      self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
      
      self.confMult:advance(tCurr, {coefIz, neutM0}, {self.coefM0})
      self.collisionSlvr:advance(tCurr, {self.coefM0, self.sumDistF}, {self.ionizSrc})
      -- Uncomment to test without fMaxwellian(Tiz)
      --self.collisionSlvr:advance(tCurr, {self.coefM0, elcDistF}, {self.ionizSrc})      
      --species[self.elcNm].distIo:write(self.ionizSrc, string.format("electron_ionizSrc_%d.bp",species[self.elcNm].distIoFrame),0,0)
      
      fRhsOut:accumulate(1.0,self.ionizSrc)
   elseif (species[self.speciesName].charge == 0) then
      -- neutrals, on Vlasov grid
      local neutDistF = species[self.neutNm]:getDistF()
      tmEvalMomStart  = Time.clock()
      self.m0elc:copy(elcM0)
      
      self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
      
      self.confMult:advance(tCurr, {coefIz, self.m0elc}, {self.coefM0})
      self.collisionSlvr:advance(tCurr, {self.coefM0, neutDistF}, {self.ionizSrc})
      fRhsOut:accumulate(-1.0,self.ionizSrc)  
   else
      -- ions 
      tmEvalMomStart = Time.clock()
      self.m0elc:copy(elcM0)
      local neutM0   = species[self.neutNm]:fluidMoments()[1]
      local neutU    = species[self.neutNm]:selfPrimitiveMoments()[1] 
      local neutVtSq = species[self.neutNm]:selfPrimitiveMoments()[2]
      -- Include only z-component of neutU

      --species[self.elcNm].distIo:write(neutU, string.format("%s_neutU_%d.bp",self.speciesName,species[self.elcNm].distIoFrame),0,0)
      --species[self.elcNm].distIo:write(neutVtSq, string.format("%s_neutVtSq_%d.bp",self.speciesName,species[self.elcNm].distIoFrame),0,0)
      self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart

      self.confMult:advance(tCurr, {coefIz, self.m0elc}, {self.coefM0})
      species[self.speciesName].calcMaxwell:advance(tCurr, {neutM0, neutU, neutVtSq}, {self.fNeutMax})
      self.collisionSlvr:advance(tCurr, {self.coefM0, self.fNeutMax}, {self.ionizSrc})
      --species[self.elcNm].distIo:write(self.fNeutMax, string.format("%s_fNeutMax_%d.bp",self.speciesName,species[self.elcNm].distIoFrame),0,0)
      fRhsOut:accumulate(1.0,self.ionizSrc)
   end
end
   
function GkIonization:write(tm, frame)
end

function GkIonization:getIonizSrc()
   return self.ionizSrc
end

function GkIonization:slvrTime()
   local time = self.confMult.totalTime
   time = time + self.collisionSlvr.totalTime
   return time
end

function GkIonization:momTime()
    if (self.speciesName == self.elcNm) then 
       return self.calcVoronovReactRate:evalMomTime() + self._tmEvalMom
    else
       return self._tmEvalMom
    end
end

function GkIonization:projectMaxwellTime()
   local time = self.confMult.projectMaxwellTime()
   time = time + self.collisionSlvr.projectMaxwellTime
   return time
end

return GkIonization

