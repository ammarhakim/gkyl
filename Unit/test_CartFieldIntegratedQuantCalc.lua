-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to compute integrated quantities.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

local function test_ser_V()
   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0,  6.0},
      cells = {1, 16},
   }
   local confGrid = Grid.RectCart {
      lower = { phaseGrid:lower(1) },
      upper = { phaseGrid:upper(1) },
      cells = { phaseGrid:numCells(1) },
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   local confBasis  = Basis.CartModalSerendipity { ndim = 1, polyOrder = 2 }
   -- Fields.
   local distf = DataStruct.Field {
      onGrid        = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost         = {0, 0},
   }
   local numDensity = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {0, 0},
   }

   -- Updater to initialize distribution function
   local project = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,
      basis    = phaseBasis,
      evaluate = function (t, xn)
         return 1/((phaseGrid:upper(1)-phaseGrid:lower(1))*(phaseGrid:upper(2)-phaseGrid:lower(2)))
      end
   }
   project:advance(0.0, {}, {distf})

   -- Moment updater.
   local calcNumDensity = Updater.DistFuncMomentCalc {
      advanceArgs = {{distf}, {numDensity}},
      onGrid      = phaseGrid,
      phaseBasis  = phaseBasis,
      confBasis   = confBasis,
      moment      = "M0",
   }
   calcNumDensity:advance(0.0, {distf}, {numDensity})

   -- Compute integrated f and integrated moment.
   local intF = DataStruct.DynVector {
      numComponents = 1,
   }
   local intM0 = DataStruct.DynVector {
      numComponents = 1,
   }
   local intQuantP = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = phaseGrid,
      basis         = phaseBasis,
      numComponents = 1,
      quantity      = "V",
   }
   local intQuantC = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = confGrid,
      basis         = confBasis,
      numComponents = 1,
      quantity      = "V",
   }

   intQuantP:advance(0.0, {distf}, {intF})
   intQuantC:advance(0.0, {numDensity}, {intM0})

   local _, NP = intF:lastData()
   local _, NC = intM0:lastData()

   assert_equal(1, NP[1], "Checking moment")
   assert_equal(1, NC[1], "Checking moment")
end

local function test_ser_V2()
   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0,  6.0},
      cells = {1, 16},
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   -- Fields.
   local distf = DataStruct.Field {
      onGrid        = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost         = {0, 0},
   }

   -- Updater to initialize distribution function
   local project = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,
      basis    = phaseBasis,
      evaluate = function (t, xn)
         return 1/math.sqrt((phaseGrid:upper(1)-phaseGrid:lower(1))*(phaseGrid:upper(2)-phaseGrid:lower(2)))
      end
   }
   project:advance(0.0, {}, {distf})

   -- Compute integrated f and integrated moment.
   local fL2 = DataStruct.DynVector {
      numComponents = 1,
   }
   local intQuantP = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = phaseGrid,
      basis         = phaseBasis,
      numComponents = 1,
      quantity      = "V2",
   }

   intQuantP:advance(0.0, {distf}, {fL2})

   local _, intL2 = fL2:lastData()

   assert_equal(1, intL2[1], "Checking moment")
end

local function test_ser_AbsV()
   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0,  6.0},
      cells = {1, 16},
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   -- Fields.
   local distf = DataStruct.Field {
      onGrid        = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost         = {0, 0},
   }

   -- Updater to initialize distribution function
   local project = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,
      basis    = phaseBasis,
      evaluate = function (t, xn)
         return -1/((phaseGrid:upper(1)-phaseGrid:lower(1))*(phaseGrid:upper(2)-phaseGrid:lower(2)))
      end
   }
   project:advance(0.0, {}, {distf})

   -- Compute integrated f and integrated moment.
   local absF = DataStruct.DynVector {
      numComponents = 1,
   }
   local intQuantP = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = phaseGrid,
      basis         = phaseBasis,
      numComponents = 1,
      quantity      = "AbsV",
   }

   intQuantP:advance(0.0, {distf}, {absF})

   local _, intAbsF = absF:lastData()

   assert_equal(1, intAbsF[1], "Checking moment")
end

test_ser_V()
test_ser_V2()
test_ser_AbsV()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
