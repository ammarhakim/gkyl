-- Gkyl ------------------------------------------------------------------------
--
-- Test for diagnostic updater
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local CalcDiagnostic = require "Updater.CalcDiagnostic"
local ffi  = require "ffi"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"

local assert_equal = Unit.assert_equal
local stats = Unit.stats


function test_1()
   local decomp = DecompRegionCalc.CartProd { cuts = {1, 1} }
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {2, 2},
      decomposition = decomp
   }
   
   local q = DataStruct.Field {
      onGrid = grid,
      numComponents = 1,
      ghost = {1, 1}
   }
   q:clear(10.5)
   
   local diagCalc = CalcDiagnostic {
      onGrid = grid,
      diagnostic = function (t, v) return v[1] end
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }

   diagCalc:advance(0.0, 0.0, {q}, {dynVec})

   local tm, lv = dynVec:lastData()
   assert_equal(10.5*q:globalRange():volume(), lv[1], "Checking if diagCalc worked")
end

function test_2()
   local decomp = DecompRegionCalc.CartProd { cuts = {1, 1} }
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {2, 2},
      decomposition = decomp
   }
   
   local q = DataStruct.Field {
      onGrid = grid,
      numComponents = 1,
      ghost = {1, 1}
   }
   q:clear(10.5)
   local indexer = q:genIndexer()
   local qItr = q:get(indexer({1,1}))
   qItr[1] = 100.5

   local diagCalc = CalcDiagnostic {
      onGrid = grid,
      operator = "max",
      diagnostic = function (t, v) return v[1] end
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }

   diagCalc:advance(0.0, 0.0, {q}, {dynVec})

   local tm, lv = dynVec:lastData()
   assert_equal(100.5, lv[1], "Checking if diagCalc worked")
end

function test_3()
   local decomp = DecompRegionCalc.CartProd { cuts = {1, 1} }
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {1.0, 1.0},
      cells = {2, 2},
      decomposition = decomp
   }
   
   local q = DataStruct.Field {
      onGrid = grid,
      numComponents = 1,
      ghost = {1, 1}
   }
   q:clear(10.5)
   local indexer = q:genIndexer()
   local qItr = q:get(indexer({1,1}))
   qItr[1] = -200.5

   local diagCalc = CalcDiagnostic {
      onGrid = grid,
      operator = "min",
      diagnostic = function (t, v) return v[1] end
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }

   diagCalc:advance(0.0, 0.0, {q}, {dynVec})

   local tm, lv = dynVec:lastData()
   assert_equal(-200.5, lv[1], "Checking if diagCalc worked")
end

test_1()
test_2()
test_3()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
