-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Geometric sources for axisymmetric perfectly-
-- hyperbolic Maxwell's equations.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local SourceBase = require "App.Sources.SourceBase"
local Proto = require "Lib.Proto"
local Updater = require "Updater"

local AxisymmetricPhMaxwellsSource = Proto(SourceBase)

function AxisymmetricPhMaxwellsSource:init(tbl)
   self.tbl = tbl
end

function AxisymmetricPhMaxwellsSource:fullInit(appTbl)
   self.cfl = nil
   self.grid = nil
   self.slvr = nil
end

function AxisymmetricPhMaxwellsSource:createSolver(species, field)
   local tbl = self.tbl
   self.slvr = Updater.AxisymmetricPhMaxwellSrc {
      onGrid = self.grid,
      scheme = tbl.timeStepper,
      evolve = tbl.evolve,
      epsilon0 = tbl.epsilon0,
      mu0 = tbl.mu0,
      chi_e = tbl.chi_e,
      chi_m = tbl.chi_m
   }
end

function AxisymmetricPhMaxwellsSource:updateSource(
      tCurr, dt, speciesVar, fieldVar)
   local outVars = {fieldVar}
   self.slvr:setDtAndCflRate(dt, nil)
   local status, dtSuggested = self.slvr:advance(tCurr, {}, outVars)
   return status, dtSuggested
end

function AxisymmetricPhMaxwellsSource:setName(nm)
   self.name = nm
end

function AxisymmetricPhMaxwellsSource:totalTime()
  return self.slvr.totalTime
end

function AxisymmetricPhMaxwellsSource:setConfGrid(grid)
   self.grid = grid
end

return AxisymmetricPhMaxwellsSource
