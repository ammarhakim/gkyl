-- Gkyl ------------------------------------------------------------------------
--
-- VlasovOnCartGrid support code: Various field objects
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis = require "Basis"
local BoundaryCondition = require "Updater.BoundaryCondition"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local LinearTrigger = require "LinearTrigger"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local PerfMaxwell = require "Eq.PerfMaxwell"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local Updater = require "Updater"
local date = require "Lib.date"
local xsys = require "xsys"

local CollisionsBase = Proto()

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
function BgkCollisions:fullInit(vlasovTbl)
   local tbl = self.tbl -- previously store table

   self.speciesList = tbl.species
   self.collFreq = tbl.collFreq

   assert(#self.speciesList == #self.collFreq, "'nu' must be defined for each 'species'")
end

function BgkCollisions:setName(nm)
   self.name = nm
end

function BgkCollisions:setConfBasis(basis)
   self.confBasis = basis
end
function BgkCollisions:setConfGrid(cgrid)
   self.confGrid = cgrid
end

function BgkCollisions:setPhaseBasis(species)
   self.phaseBasis = {}
   for _, nm in pairs(self.speciesList) do
      self.phaseBasis[nm] = species[nm].basis
   end
end

function BgkCollisions:setPhaseGrid(species)
   self.phaseGrid = {}
   for _, nm in pairs(self.speciesList) do
      self.phaseGrid[nm] = species[nm].grid
   end
end

-- methods for Bgk collisions object
function BgkCollisions:setCfl(cfl) self.cfl = cfl end

function BgkCollisions:createSolver(species)
   local confBasis = nil
   local confGrid = nil
   local phaseBasis = {}
   local phaseGrid = {}
   for _, nm in pairs(self.speciesList) do
      confBasis = species[nm].confBasis
      confGrid = species[nm].confGrid
      phaseBasis[nm] = species[nm].basis
      phaseGrid[nm] = species[nm].grid
   end
   self.collisionSlvr = Updater.BgkCollisions {
      confGrid = confGrid,
      confBasis = confBasis,
      phaseGrid = phaseGrid,
      phaseBasis = phaseBasis,
      speciesList = self.speciesList,
      collFreq = self.collFreq,
   }
end

function BgkCollisions:forwardEuler(tCurr, dt, idxIn, outIdx, speciesList)
   local spRkFields, spMomFields = {},  {}

   for i, nm in ipairs(species) do
      spOutFields[i] = speciesList[nm]:rkStepperFields()[outIdx]
      spMomFields[i] = speciesList[nm]:fluidMoments()
   end
   return self.collisionSlvr:advance(tCurr, dt, spMomFields, spOutFields)
end

   
function BgkCollisions:totalSolverTime()
   return self.collisionSlvr.totalTime + self.tmCurrentAccum
end

return {
   CollisionsBase = CollisionsBase,
   BgkCollisions = BgkCollisions,
}
