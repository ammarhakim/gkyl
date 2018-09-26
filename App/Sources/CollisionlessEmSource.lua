-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Collisionless EM sources for use in fluid sims
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local DataStruct = require "DataStruct"
local Proto = require "Lib.Proto"
local Updater = require "Updater"
local xsys = require "xsys"

-- CollisionlessEmSource ---------------------------------------------------------------
--
-- Coupled Lorentz and current sources
--------------------------------------------------------------------------------

local CollisionlessEmSource = Proto(SourceBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function CollisionlessEmSource:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function CollisionlessEmSource:fullInit(appTbl)
   local tbl = self.tbl -- previously store table

   -- set later by various set() methods
   self.cfl = nil
   self.grid = nil
   self.slvr = nil

   self.speciesList = tbl.species -- list of species to update

   -- (MAKE SURE NAME IS MATCHING ONE IN UPDATERS!!)
   self.timeStepper = tbl.timeStepper -- time-stepper to use
end

function CollisionlessEmSource:setName(nm)
   self.name = nm
end

function CollisionlessEmSource:setCfl(cfl)
   self.cfl = cfl
end
function CollisionlessEmSource:setConfGrid(grid)
   self.grid = grid
end

function CollisionlessEmSource:createSolver(species, field)
   local numSpecies = #self.speciesList
   local mass, charge = {}, {}

   for i, nm in ipairs(self.speciesList) do
      mass[i] = species[nm]:getMass()
      charge[i] = species[nm]:getCharge()
   end

   self.slvr = Updater.FiveMomentSrc {
      onGrid = self.grid,
      numFluids = numSpecies,
      charge = charge,
      mass = mass,
      epsilon0 = field:getEpsilon0(),
      elcErrorSpeedFactor = field:getElcErrorSpeedFactor(),
      mgnErrorSpeedFactor = field:getMgnErrorSpeedFactor(),
      scheme = self.timeStepper,
   }
end

function CollisionlessEmSource:forwardEuler(tCurr, dt, fIn, species, fOut)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   return status, dtSuggested
end

function CollisionlessEmSource:updateSource(tCurr, dt, speciesVar, fieldVar)
   local outVars = {}
   for i, nm in ipairs(self.speciesList) do
      outVars[i] = speciesVar[nm]
   end
   outVars[#self.speciesList+1] = fieldVar
   return self.slvr:advance(tCurr, dt, {}, outVars)
end

function CollisionlessEmSource:write(tm, frame)
end

function CollisionlessEmSource:totalTime()
   return 0.0
end

return CollisionlessEmSource
