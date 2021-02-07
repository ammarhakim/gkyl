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

-- TenMomentGradSource ---------------------------------------------------------------
--
-- Coupled Lorentz and current sources
--------------------------------------------------------------------------------

local TenMomentGradSource = Proto(SourceBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function TenMomentGradSource:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function TenMomentGradSource:fullInit(appTbl)
   local tbl = self.tbl -- previously store table

   -- set later by various set() methods
   self.cfl = nil
   self.grid = nil
   self.slvr = nil

   self.speciesList = tbl.species -- list of species to update
   -- factor multiplying the thermal conductivity
   self.alpha = tbl.alpha
end

function TenMomentGradSource:setName(nm)
   self.name = nm
end

function TenMomentGradSource:setCfl(cfl)
   self.cfl = cfl
end
function TenMomentGradSource:setConfGrid(grid)
   self.grid = grid
end

function TenMomentGradSource:createSolver(species, field)
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

   self.slvr = Updater.TenMomentGrad {
      onGrid = self.grid,
      alpha = self.alpha
   }
   -- Array with ten components per cell to store symmetrized heat flux.
   self.q = DataStruct.Field {
      onGrid = self.grid,
      numComponents = 10,
      ghost = {2, 2},
   }
   self.q:clear(0.0)
end

function TenMomentGradSource:forwardEuler(tCurr, dt, fIn, species, fOut)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   return status, dtSuggested
end

function TenMomentGradSource:updateSource(tCurr, dt, speciesVar, fieldVar)
   local outVars = {}
   for i, nm in ipairs(self.speciesList) do
      outVars[i] = speciesVar[nm]
   end
   outVars[#self.speciesList+1] = fieldVar
   self.slvr:setDtAndCflRate(dt, nil)
   return self.slvr:advance(tCurr, {self.q}, outVars)
end

function TenMomentGradSource:write(tm, frame)
end

function TenMomentGradSource:totalTime()
   return self.slvr.totalTime
end

return TenMomentGradSource