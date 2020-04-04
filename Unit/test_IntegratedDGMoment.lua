-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to compute integrated moments.
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

local function test_one()
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
      onGrid     = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "M0",
   }
   calcNumDensity:advance(0.0, {distf}, {numDensity})

   -- Compute integrated f and integrated moment with CartFieldIntegratedQuantCalc.
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

   -- Compute integral of f and of zeroth with IntegratedDGMoment.
   local intFn = DataStruct.DynVector {
      numComponents = 1,
   }
   local intM0n = DataStruct.DynVector {
      numComponents = 1,
   }
   local intMomP = Updater.IntegratedDGMoment {
      onGrid = phaseGrid,
      basis  = phaseBasis,
      moment = "one",
   }
   local intMomC = Updater.IntegratedDGMoment {
      onGrid = confGrid,
      basis  = confBasis,
      moment = "one",
   }
   intMomP:advance(0.0, {distf}, {intFn})
   intMomC:advance(0.0, {numDensity}, {intM0n})
   local _, NPn = intFn:lastData()
   local _, NCn = intM0n:lastData()

   assert_equal(NP[1], NPn[1], "Checking moment")
   assert_equal(NP[1], NCn[1], "Checking moment")
end

local function test_vSq()
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
   local energyDensity = DataStruct.Field {
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
   local calcEnergyDensity = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "M2",
   }
   calcEnergyDensity:advance(0.0, {distf}, {energyDensity})

   -- Compute integrated f and integrated moment with CartFieldIntegratedQuantCalc.
   local intM2 = DataStruct.DynVector {
      numComponents = 1,
   }
   local intQuantC = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = confGrid,
      basis         = confBasis,
      numComponents = 1,
      quantity      = "V",
   }
   intQuantC:advance(0.0, {energyDensity}, {intM2})
   local _, NC = intM2:lastData()

   -- Compute integral of f and of zeroth with IntegratedDGMoment.
   local intFn = DataStruct.DynVector {
      numComponents = 1,
   }
   local intMomP = Updater.IntegratedDGMoment {
      onGrid    = phaseGrid,
      basis     = phaseBasis,
      confBasis = confBasis,
      moment    = "vSq",
   }
   intMomP:advance(0.0, {distf}, {intFn})
   local _, NPn = intFn:lastData()

   assert_equal(NC[1], NPn[1], "Checking moment")
end

local function test_v1()
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
   local momDensity = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
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

   -- Moment updater.
   local calcMomDensity = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "M1i",
   }
   calcMomDensity:advance(0.0, {distf}, {momDensity})

   -- Compute integrated M1_1 moment with CartFieldIntegratedQuantCalc.
   local intM1_1 = DataStruct.DynVector {
      numComponents = 1,
   }
   local intQuantC = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = phaseGrid,
      basis         = phaseBasis,
      numComponents = 1,
      quantity      = "V",
   }
   intQuantC:advance(0.0, {momDensity}, {intM1_1})
   local _, intM1 = intM1_1:lastData()

   -- Compute integral of v_1*f with IntegratedDGMoment.
   local intFn = DataStruct.DynVector {
      numComponents = 1,
   }
   local intMomP = Updater.IntegratedDGMoment {
      onGrid    = phaseGrid,
      basis     = phaseBasis,
      confBasis = confBasis,
      moment    = "v1",
   }
   intMomP:advance(0.0, {distf}, {intFn})
   local _, intM1n = intFn:lastData()

   local intFnA = DataStruct.DynVector {
      numComponents = 1,
   }
   local intMomPA = Updater.IntegratedDGMoment {
      onGrid    = phaseGrid,
      basis     = phaseBasis,
      confBasis = confBasis,
      moment    = "x2",
   }
   intMomPA:advance(0.0, {distf}, {intFnA})
   local _, intM1nA = intFnA:lastData()

   assert_equal(intM1[1], intM1n[1], "Checking moment")
   assert_equal(intM1n[1], intM1nA[1], "Checking moment")
end

local function test_v2()
   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0, -6},
      upper = {1.0,  6.0, -6},
      cells = {1, 16, 16},
   }
   local confGrid = Grid.RectCart {
      lower = { phaseGrid:lower(1) },
      upper = { phaseGrid:upper(1) },
      cells = { phaseGrid:numCells(1) },
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = 3, polyOrder = 2 }
   local confBasis  = Basis.CartModalSerendipity { ndim = 1, polyOrder = 2 }
   -- Fields.
   local distf = DataStruct.Field {
      onGrid        = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost         = {0, 0},
   }
   local momDensity = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
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

   -- Moment updater.
   local calcMomDensity = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "M1i",
   }
   calcMomDensity:advance(0.0, {distf}, {momDensity})

   -- Compute integrated M1_1 moment with CartFieldIntegratedQuantCalc.
   local intM1_1 = DataStruct.DynVector {
      numComponents = 1,
   }
   local intQuantC = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = phaseGrid,
      basis         = phaseBasis,
      numComponents = 1,
      quantity      = "V",
   }
   intQuantC:advance(0.0, {momDensity}, {intM1_1})
   local _, intM1 = intM1_1:lastData()

   -- Compute integral of v_2*f with IntegratedDGMoment.
   local intFn = DataStruct.DynVector {
      numComponents = 1,
   }
   local intMomP = Updater.IntegratedDGMoment {
      onGrid    = phaseGrid,
      basis     = phaseBasis,
      confBasis = confBasis,
      moment    = "v2",
   }
   intMomP:advance(0.0, {distf}, {intFn})
   local _, intM1n = intFn:lastData()

   local intFnA = DataStruct.DynVector {
      numComponents = 1,
   }
   local intMomPA = Updater.IntegratedDGMoment {
      onGrid    = phaseGrid,
      basis     = phaseBasis,
      confBasis = confBasis,
      moment    = "x3",
   }
   intMomPA:advance(0.0, {distf}, {intFnA})
   local _, intM1nA = intFnA:lastData()

   assert_equal(intM1[2], intM1n[1], "Checking moment")
   assert_equal(intM1n[1], intM1nA[1], "Checking moment")
end

test_one()
test_vSq()
test_v1()
test_v2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
