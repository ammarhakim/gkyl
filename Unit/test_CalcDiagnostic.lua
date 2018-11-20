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
      upper = {2.0, 2.0},
      cells = {10, 10},
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
   assert_equal(4*10.5, lv[1], "Checking if diagCalc worked")
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
