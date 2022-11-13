-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Inter-species friction for multifluid five-
-- and ten-moment models.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local FluidSourceBase = require "App.FluidSources.FluidSourceBase"
local Proto = require "Lib.Proto"
local Updater = require "Updater"
local DataStruct = require "DataStruct"
local Basis = require "Basis"

local MomentUserSource = Proto(FluidSourceBase)

function MomentUserSource:init(tbl)
   self.tbl = tbl
end

function MomentUserSource:fullInit(appTbl)
   self.cfl = nil
   self.grid = nil
   self.slvr = nil
end

function MomentUserSource:createSolver(species, field)
   local tbl = self.tbl

   -- Allow different species to have different number of moments.
   local speciesLists = {}
   speciesLists[5] = {}
   speciesLists[10] = {}

   local evolveLists = {}
   evolveLists[5] = {}
   evolveLists[10] = {}

   self.sourceLists = {}
   self.sourceLists[5] = {}
   self.sourceLists[10] = {}

   local evolve = tbl.evolve
   local mass = {}

   local basis = Basis.CartModalMaxOrder {
      ndim = self.grid:ndim(),
      polyOrder = 0
   }

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
      mass[i] = species[nm]:getMass()

      local source = DataStruct.Field {
         onGrid = self.grid,
         numComponents = 2,
         ghost = {0, 0}
      }
      local project = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = basis,
         evaluate = tbl.source_funcs[i],
      }
      project:advance(0.0, {}, {source})
      self.sourceLists[nMoments][nm] = source
   end
   self.speciesLists = speciesLists

   self.slvrs = {}
   if #speciesLists[5]>=0 then
      self.slvrs[5] = Updater.FiveMomentUserSrc {
         onGrid = self.grid,
         scheme = tbl.timeStepper,
         numFluids = #speciesLists[5],
         evolve = evolveLists[5],
         mass = mass,
         gasGamma = tbl.gasGamma,
         kBoltzmann = tbl.kBoltzmann,
      }
   end
   if #speciesLists[10]>0 then
      assert(false, "10-moment is not supported yet.")
   end
end

function MomentUserSource:updateFluidSource(tCurr, dt, speciesVar, fieldVar)
   local status = true
   local dtSuggested = GKYL_MAX_DOUBLE

   for _,nMoments in ipairs({5, 10}) do
      local speciesList = self.speciesLists[nMoments]

      if #speciesList > 0 then
         local outVars = {}
         local srcVars = {}
         local slvr = self.slvrs[nMoments]

         for i, nm in ipairs(speciesList) do
            outVars[i] = speciesVar[nm]
            srcVars[i] = self.sourceLists[nMoments][nm]
         end
         slvr:setDtAndCflRate(dt, nil)

         local myStatus, myDtSuggested = slvr:advance(tCurr, srcVars, outVars)

         status = status and myStatus
         dtSuggested = math.min(dtSuggested, myDtSuggested)
      end
   end

   return status, dtSuggested
end

function MomentUserSource:setName(nm)
   self.name = nm
end

function MomentUserSource:totalTime()
   local totalTime = 0
   for _,slvr in ipairs(self.slvrs) do
      totalTime = totalTime + slvr.totalTime
  end
  return totalTime
end

function MomentUserSource:setConfGrid(grid)
   self.grid = grid
end

return MomentUserSource
