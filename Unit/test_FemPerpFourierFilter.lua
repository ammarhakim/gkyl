-- Gkyl ------------------------------------------------------------------------
--
-- Tests for Fourier filter updater
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"
local Lin        = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats        = Unit.stats

function test_1(nx, ny)
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {2*math.pi, 2*math.pi},
      cells = {nx, ny},
   }
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 1 }
   local check = DataStruct.Field {
      onGrid = grid,
      numComponents = basis:numBasis(),
      ghost = {1, 1},
   }
   local src = DataStruct.Field {
      onGrid = grid,
      numComponents = basis:numBasis(),
      ghost = {1, 1},
   }
   local sol = DataStruct.Field {
      onGrid = grid,
      numComponents = basis:numBasis(),
      ghost = {1, 1},
   }
   local initSrc = Updater.EvalOnNodes {
      onGrid   = grid,
      basis    = basis,
      evaluate = function (t,xn)
         local x = xn[1]
         local y = xn[2]
         return 0.5*math.cos(y) + 1.1*math.cos(2*y) + .3*math.cos(3*y)
      end,
   }
   initSrc:advance(0.,{},{src})

   local initCheck = Updater.EvalOnNodes {
      onGrid   = grid,
      basis    = basis,
      evaluate = function (t,xn)
         local x = xn[1]
         local y = xn[2]
         return 0.5*math.cos(y) + .3*math.cos(3*y)
      end
   }
   initCheck:advance(0.,{},{check})

   local filter = Updater.FemPerpFourierFilter {
      onGrid = grid,
      basis = basis,
      kyfilter = {1, 3},
   }
   filter:advance(0., {src}, {sol})

   local err = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1}
   }

   --src:write("src.bp")
   --sol:write("sol.bp")

   err:combine(1.0, check, -1.0, sol)
   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = grid,
      basis         = basis,
      numComponents = 1,
      quantity      = 'V2',
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()

   local intErr = math.sqrt(lv[1])
   --io.write("Average RMS error = ", intErr, "\n")
   assert_close(intErr, 0., 1e-14, "checking integrated RMS error in filtered solution")
   return intErr
end   
   
test_1(32, 32)

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
