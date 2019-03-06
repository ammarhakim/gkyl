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

-- TenMomentRelaxSource ---------------------------------------------------------------
--
-- Coupled Lorentz and current sources
--------------------------------------------------------------------------------

local TenMomentRelaxSource = Proto(SourceBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function TenMomentRelaxSource:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function TenMomentRelaxSource:fullInit(appTbl)
   local tbl = self.tbl -- previously store table

   -- set later by various set() methods
   self.cfl = nil
   self.grid = nil
   self.slvr = nil

   self.speciesList = tbl.species -- list of species to update

   -- (MAKE SURE NAME IS MATCHING ONE IN UPDATERS!!)
   self.timeStepper = tbl.timeStepper -- time-stepper to use

   self.k = tbl.k
   self.hasKField = tbl.hasKField
   self.hasEm = tbl.hasEm
   self.hasStatic = tbl.hasStatic
end

function TenMomentRelaxSource:setName(nm)
   self.name = nm
end

function TenMomentRelaxSource:setCfl(cfl)
   self.cfl = cfl
end
function TenMomentRelaxSource:setConfGrid(grid)
   self.grid = grid
end

function TenMomentRelaxSource:createSolver(species, field)
   local numSpecies = #self.speciesList
   local mass, charge = {}, {}

   local source_type
   for i, nm in ipairs(self.speciesList) do
      mass[i] = species[nm]:getMass()
      charge[i] = species[nm]:getCharge()
      if not source_type then
         source_type = species[nm].nMoments
      else
         -- FIXME Currently all species must have the same moments.
         assert(source_type == species[nm].nMoments)
      end
   end

   self.slvr = Updater.TenMomentRelax {
      onGrid = self.grid,
      numFluids = numSpecies,
      mass = mass,
      charge = charge,
      scheme = self.timeStepper,
      k = self.k,
      hasKField = self.hasKField,
      hasEm = self.hasEm,
      hasStatic = self.hasStatic,
   }
end

function TenMomentRelaxSource:forwardEuler(tCurr, dt, fIn, species, fOut)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   return status, dtSuggested
end

function TenMomentRelaxSource:updateSource(tCurr, dt, speciesVar, fieldVar)
   local outVars = {}
   for i, nm in ipairs(self.speciesList) do
      outVars[i] = speciesVar[nm]
   end
   outVars[#self.speciesList+1] = fieldVar
   self.slvr:setDtAndCflRate(dt, nil)
   return self.slvr:advance(tCurr, {}, outVars)
end

function TenMomentRelaxSource:write(tm, frame)
end

function TenMomentRelaxSource:totalTime()
   return self.slvr.totalTime
end

return TenMomentRelaxSource
