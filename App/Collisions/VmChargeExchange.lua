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

-- VmChargeExchange  --------------------------------------------------------
--
-- Charge Exchange Operator
--------------------------------------------------------------------------------

local VmChargeExchange = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmChargeExchange:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmChargeExchange:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously store table.

   self.collKind = "CX"

   self.collidingSpecies = assert(tbl.collideWith, "App.VmChargeExchange: Must specify names of species to collide with in 'collideWith'.")

   self.collideNm   = tbl.collideWith[1]

   self.plasma      = assert(tbl.plasma, "App.VmChargeExchange: Must specify plasma species in 'plasma' ('H', 'D', or 'Ne')")
   
   self.ionNm       = tbl.ions
   self.neutNm      = tbl.neutrals
   self.iMass       = tbl.ionMass
   self.nMass       = tbl.neutMass

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
   end

end

function VmChargeExchange:setName(nm)
   self.name = nm
end

function VmChargeExchange:setSpeciesName(nm)
   self.speciesName = nm
end

function VmChargeExchange:setCfl(cfl)
   self.cfl = cfl
end

function VmChargeExchange:setConfBasis(basis)
   self.confBasis = basis
end

function VmChargeExchange:setConfGrid(grid)
   self.confGrid = grid
end

function VmChargeExchange:setPhaseBasis(basis)
   self.phaseBasis = basis
end

function VmChargeExchange:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function VmChargeExchange:createSolver(funcField) --species)

   self.calcVrelProdCX = Updater.VrelProductCX {
	 onGrid         = self.phaseGrid,
	 confBasis      = self.confBasis,
	 phaseBasis     = self.phaseBasis,
	 kineticSpecies = 'Vm',
   }

   self.collisionSlvr = Updater.SigmaCX {
         onGrid         = self.confGrid,
         confBasis      = self.confBasis,
	 phaseBasis     = self.phaseBasis,
	 kineticSpecies = 'Vm',
	 a              = self.a,
	 b              = self.b,
   }
   self.sourceCX = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
end

function VmChargeExchange:advance(tCurr, fIn, species, fRhsOut)
   -- get CX source term from Vlasov species
   self.sourceCX = species[self.ionNm]:getSrcCX()
   
   -- identify species and accumulate
   if (self.speciesName == self.ionNm) then
      fRhsOut:accumulate(1.0,self.sourceCX)
   elseif (self.speciesName == self.neutNm) then
      fRhsOut:accumulate(-self.iMass/self.nMass,self.sourceCX)
   end
   
end
   
function VmChargeExchange:getSrcCX()
   return self.sourceCX
end

function VmChargeExchange:write(tm, frame)
end

function VmChargeExchange:slvrTime()
   return self.collisionSlvr.totalTime
end

function VmChargeExchange:momTime()
   return 0
end

function VmChargeExchange:projectMaxwellTime()
   return self.collisionSlvr:projectMaxwellTime()
end

return VmChargeExchange

