-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to compute moments
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Basis = require "Basis"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_ser_1x1v()
   -- phase-space and config-space grids
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0, 6.0},
      cells = {1, 16},
   }
   local confGrid = Grid.RectCart {
      lower = { phaseGrid:lower(1) },
      upper = { phaseGrid:upper(1) },
      cells = { phaseGrid:numCells(1) },
   }
   -- basis functions
   local phaseBasis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   local confBasis = Basis.CartModalSerendipity { ndim = 1, polyOrder = 2 }
   -- fields
   local distf = DataStruct.Field {
      onGrid = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost = {0, 0},
   }
   local numDensity = DataStruct.Field {
      onGrid = confGrid,
      numComponents = confBasis:numBasis(),
      ghost = {0, 0},
   }

   -- updater to initialize distribution function
   local project = Updater.ProjectOnBasis {
      onGrid = phaseGrid,
      basis = phaseBasis,
      evaluate = function (t, xn)
	 return 1/(phaseGrid:upper(2)-phaseGrid:lower(2))
      end
   }
   project:advance(0.0, 0.0, {}, {distf})

   -- moment updater
   local calcNumDensity = Updater.DistFuncMomentCalc {
      onGrid = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis = confBasis,
      moment = "M0",
   }
   calcNumDensity:advance(0.0, 0.0, {distf}, {numDensity})

   local momIdxr = numDensity:genIndexer()
   local nItr = numDensity:get(momIdxr( {1} ))
   assert_equal(1, nItr[1]/math.sqrt(2), "Checking moment")
   
end

function test_max_1x1v()
   -- phase-space and config-space grids
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0, 6.0},
      cells = {1, 16},
   }
   local confGrid = Grid.RectCart {
      lower = { phaseGrid:lower(1) },
      upper = { phaseGrid:upper(1) },
      cells = { phaseGrid:numCells(1) },
   }
   -- basis functions
   local phaseBasis = Basis.CartModalMaxOrder { ndim = 2, polyOrder = 2 }
   local confBasis = Basis.CartModalMaxOrder { ndim = 1, polyOrder = 2 }
   -- fields
   local distf = DataStruct.Field {
      onGrid = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost = {0, 0},
   }
   local numDensity = DataStruct.Field {
      onGrid = confGrid,
      numComponents = confBasis:numBasis(),
      ghost = {0, 0},
   }

   -- moment updater
   local calcNumDensity = Updater.DistFuncMomentCalc {
      onGrid = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis = confBasis,
      moment = "M0",
   }
   
end

test_ser_1x1v()
test_max_1x1v()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
