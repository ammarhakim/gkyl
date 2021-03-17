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

local AxisymmetricMomentsSource = Proto(SourceBase)

function AxisymmetricMomentsSource:init(tbl)
   self.tbl = tbl
end

function AxisymmetricMomentsSource:fullInit(appTbl)
   local tbl = self.tbl

   self.cfl = nil
   self.grid = nil
   self.slvr = nil

   self.speciesList = tbl.species
   self.evolve = tbl.evolve

   self.timeStepper = tbl.timeStepper
   self.gasGamma = tbl.gasGamma
   self.hasPressure = tbl.hasPressureField
end

function AxisymmetricMomentsSource:createSolver(species, field)
   local evolve = self.evolve

   local speciesLists = {}
   speciesLists[5] = {}
   speciesLists[10] = {}

   local evolveLists = {}
   evolveLists[5] = {}
   evolveLists[10] = {}

   for i, nm in ipairs(self.speciesList) do
      local nMoments = species[nm].nMoments
      assert(nMoments==5 or nMoments==10,
             string.format("nMoments %s not supported", nMoments))
      table.insert(speciesLists[nMoments], nm)
      if evolve then
         table.insert(evolveLists[nMoments], evolve[i])
      else
         table.insert(evolveLists[nMoments], species[nm]:getEvolve())
      end
   end
   self.speciesLists = speciesLists

   self.slvrs = {}
   if #speciesLists[5]>=0 then
      self.slvrs[5] = Updater.AxisymmetricFiveMomentSrc {
         onGrid = self.grid,
         scheme = self.timeStepper,
         numFluids = #speciesLists[5],
         evolve = evolveLists[5],
         gasGamma = self.gasGamma,
         hasPressure = self.hasPressureField,
      }
   elseif #speciesLists[10]>0 then
      assert(false, "10-moment not implemented")
   end
end

function AxisymmetricMomentsSource:updateSource(tCurr, dt, speciesVar, fieldVar)
   local status = true
   local dtSuggested = GKYL_MAX_DOUBLE

   for _,nMoments in ipairs({5, 10}) do
      local speciesList = self.speciesLists[nMoments]

      if #speciesList > 0 then
         local outVars = {}
         local slvr = self.slvrs[nMoments]

         for i, nm in ipairs(speciesList) do
            outVars[i] = speciesVar[nm]
         end
         slvr:setDtAndCflRate(dt, nil)

         local myStatus, myDtSuggested = slvr:advance(tCurr, {}, outVars)

         status = status and myStatus
         dtSuggested = math.min(dtSuggested, myDtSuggested)
      end
   end

   return status, dtSuggested
end

function AxisymmetricMomentsSource:setName(nm)
   self.name = nm
end

function AxisymmetricMomentsSource:totalTime()
   local totalTime = 0
   for _,slvr in ipairs(self.slvrs) do
      totalTime = totalTime + slvr.totalTime
  end
  return totalTime
end

function AxisymmetricMomentsSource:setConfGrid(grid)
   self.grid = grid
end

return AxisymmetricMomentsSource
