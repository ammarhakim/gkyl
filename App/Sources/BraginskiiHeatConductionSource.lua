-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code:  Braginskii Heat Conduction Source, i.e., the
-- div(q) terms.
--
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local Proto = require "Lib.Proto"
local Updater = require "Updater"

local BraginskiiHeatConductionSource = Proto(SourceBase)

function BraginskiiHeatConductionSource:init(tbl)
   self.tbl = tbl
end

function BraginskiiHeatConductionSource:fullInit(appTbl)
   local tbl = self.tbl -- previously stored table

   self.cfl = nil
   self.grid = nil
   self.slvr = nil
end

function BraginskiiHeatConductionSource:setName(nm)
   self.name = nm
end

function BraginskiiHeatConductionSource:setCfl(cfl)
   self.cfl = cfl
end
function BraginskiiHeatConductionSource:setConfGrid(grid)
   self.grid = grid
end

function BraginskiiHeatConductionSource:createSolver(species, field)
   local tbl = self.tbl

   local mass, charge = {}, {}
   for i, nm in ipairs(tbl.species) do
      mass[i] = species[nm]:getMass()
      charge[i] = species[nm]:getCharge()
   end

   self.slvr = Updater.BraginskiiHeatConduction {
      onGrid           = self.grid,
      numFluids        = #tbl.species,
      mass             = mass,
      charge           = charge,
      gasGamma         = tbl.gasGamma,
      epsilon0         = tbl.epsilon0,
      tau              = tbl.tau,
      calcTau          = tbl.calcTau,
      coulombLogarithm = tbl.coulombLogarithm,
      coordinate       = tbl.coordinate,
   }
end

function BraginskiiHeatConductionSource:updateSource(
      tCurr, dt, speciesVar, fieldVar, speciesBuf, fieldBuf, allSpecies, field)
   local tbl = self.tbl

   local outVars = {}
   for i, nm in ipairs(tbl.species) do
      outVars[i] = speciesVar[nm]

      for d=1,self.grid:ndim() do
         allSpecies[nm]:applyBc(tCurr, speciesVar[nm], d)
      end
   end
   outVars[#tbl.species+1] = fieldVar
   field:applyBc(tCurr, fieldVar)

   local inVars = {}
   for i, nm in ipairs(tbl.species) do
      inVars[i] = speciesBuf[nm]
   end
   inVars[#tbl.species+1] = fieldBuf

   self.slvr:setDtAndCflRate(dt, nil)
   return self.slvr:advance(tCurr, inVars, outVars)
end

function BraginskiiHeatConductionSource:write(tm, frame)
end

function BraginskiiHeatConductionSource:totalTime()
   return self.slvr.totalTime
end

return BraginskiiHeatConductionSource
