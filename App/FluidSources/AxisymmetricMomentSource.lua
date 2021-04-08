-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Geometric sources for axisymmetric 
-- multi-moment equations.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local FluidSourceBase = require "App.FluidSources.FluidSourceBase"
local DataStruct = require "DataStruct"
local Proto = require "Lib.Proto"
local Updater = require "Updater"

local AxisymmetricMomentsSource = Proto(FluidSourceBase)

function AxisymmetricMomentsSource:init(tbl)
   self.tbl = tbl
end

function AxisymmetricMomentsSource:fullInit(appTbl)
   self.cfl = nil
   self.grid = nil
   self.slvr = nil
end

function AxisymmetricMomentsSource:createSolver(species, field)
   local tbl = self.tbl

   -- Allow different species to have different number of moments.
   local speciesLists = {}
   speciesLists[5] = {}
   speciesLists[10] = {}

   local evolveLists = {}
   evolveLists[5] = {}
   evolveLists[10] = {}

   local evolve = tbl.evolve
   for i, nm in ipairs(tbl.species) do
      local nMoments = species[nm].nMoments
      assert(nMoments==5 or nMoments==10,
             string.format("%s-moment is not supported.", nMoments))
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
         scheme = tbl.timeStepper,
         numFluids = #speciesLists[5],
         evolve = evolveLists[5],
         gasGamma = tbl.gasGamma,
         hasPressure = tbl.hasPressureField,
      }
   elseif #speciesLists[10]>0 then
      assert(false, "10-moment is not supported yet.")
   end
end

function AxisymmetricMomentsSource:updateFluidSource(tCurr, dt, speciesVar, fieldVar)
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
