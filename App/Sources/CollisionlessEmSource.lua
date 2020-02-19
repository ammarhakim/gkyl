-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Collisionless EM sources for use in fluid sims
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local DataStruct = require "DataStruct"
local Basis = require "Basis"
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
   self.evolve = tbl.evolve

   -- (MAKE SURE NAME IS MATCHING ONE IN UPDATERS!!)
   self.timeStepper = tbl.timeStepper -- time-stepper to use
   self.linSolType = tbl.linSolType

   self.hasStaticField = tbl.hasStaticField
   self.staticEmFunc = tbl.staticEmFunction
   self.staticEm = nil

   self.hasPressure = tbl.hasPressureField
   self.hasSigmaField = tbl.hasSigmaField
   self.hasAuxSourceFunction = tbl.hasAuxSourceFunction
   self.auxSourceFunction = tbl.auxSourceFunction
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
   local evolve = self.evolve

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
   if not evolve then
      evolve = {}
      for i, nm in ipairs(self.speciesList) do
         evolve[i] = species[nm]:getEvolve()
      end
   end

   if self.hasStaticField then
      local ndim = self.grid:ndim()
      local polyOrder = 0
      self.basis = Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
      self.staticEm = DataStruct.Field {
         onGrid = self.grid,
         numComponents = 6,
         ghost = {2, 2}
      }
      local project = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = self.staticEmFunc,
         projectOnGhosts = true,
      }
      project:advance(0.0, {}, {self.staticEm})
   end

   if source_type == 5 then
      self.slvr = Updater.FiveMomentSrc {
         onGrid = self.grid,
         numFluids = numSpecies,
         charge = charge,
         mass = mass,
         evolve = evolve,
         epsilon0 = field:getEpsilon0(),
         elcErrorSpeedFactor = field:getElcErrorSpeedFactor(),
         mgnErrorSpeedFactor = field:getMgnErrorSpeedFactor(),
         linSolType = self.linSolType,
         scheme = self.timeStepper,
         hasStaticField = self.hasStaticField,
         hasPressure = self.hasPressureField,
         hasSigmaField = self.hasSigmaField,
         sigmaField = self.sigmaField,
         hasAuxSourceFunction = self.hasAuxSourceFunction,
         auxSourceFunction = self.auxSourceFunction,
      }
   elseif source_type == 10 then
      self.slvr = Updater.TenMomentSrc {
         onGrid = self.grid,
         numFluids = numSpecies,
         charge = charge,
         mass = mass,
         evolve = evolve,
         epsilon0 = field:getEpsilon0(),
         elcErrorSpeedFactor = field:getElcErrorSpeedFactor(),
         mgnErrorSpeedFactor = field:getMgnErrorSpeedFactor(),
         linSolType = self.linSolType,
         scheme = self.timeStepper,
         hasStaticField = self.hasStaticField,
         hasPressure = self.hasPressureField,
         hasSigmaField = self.hasSigmaField,
         sigmaField = self.sigmaField,
         hasAuxSourceFunction = self.hasAuxSourceFunction,
         auxSourceFunction = self.auxSourceFunction,
      }
   else
      assert(false, string.format("source_type %s not supported.", source_type))
   end
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
   self.slvr:setDtAndCflRate(dt, nil)
   return self.slvr:advance(tCurr, {self.staticEm}, outVars)
end

function CollisionlessEmSource:write(tm, frame)
end

function CollisionlessEmSource:totalTime()
   return self.slvr.totalTime
end

return CollisionlessEmSource
