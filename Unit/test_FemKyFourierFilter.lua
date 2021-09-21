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

function test_2d(nx, ny)
   local grid = Grid.RectCart {
      lower = {0.0, 0.0},
      upper = {2*math.pi, 2},
      cells = {nx, ny},
   }
   local Ly = grid:upper(2) - grid:lower(2)
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
         local ky1 = 2*math.pi/Ly
         local ky2 = 8*math.pi/Ly
         local ky3 = 6*math.pi/Ly
         return 0.5*math.cos(ky1*y) + 1.1*math.cos(ky2*y)*x + .3*math.cos(ky3*y)*(x+2)
      end,
      onGhosts = true,
   }
   initSrc:advance(0.,{},{src})

   local initCheck = Updater.EvalOnNodes {
      onGrid   = grid,
      basis    = basis,
      evaluate = function (t,xn)
         local x = xn[1]
         local y = xn[2]
         local ky1 = 2*math.pi/Ly
         local ky3 = 6*math.pi/Ly
         return 0.5*math.cos(ky1*y) + .3*math.cos(ky3*y)*(x+2)
      end
   }
   initCheck:advance(0.,{},{check})

   local filter = Updater.FemKyFourierFilter {
      onGrid = grid,
      basis = basis,
      kyfilter = {2*math.pi/Ly, 6*math.pi/Ly},
   }
   filter:advance(0., {src}, {sol})

   local err = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1}
   }

   --src:write("src.bp")
   --sol:write("sol.bp")
   --check:write("check.bp")

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

function test_3d(nx, ny, nz)
   local grid = Grid.RectCart {
      lower = {0.0, 0.0, 0.0},
      upper = {2*math.pi, 2, 1.0},
      cells = {nx, ny, nz},
   }
   local Ly = grid:upper(2) - grid:lower(2)
   local basis = Basis.CartModalSerendipity { ndim = 3, polyOrder = 1 }
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
         local z = xn[3]
         local ky1 = 2*math.pi/Ly
         local ky2 = 8*math.pi/Ly
         local ky3 = 6*math.pi/Ly
         return 0.5*math.cos(ky1*y) + 1.1*math.cos(ky2*y) + 0.3*math.cos(ky3*y)*x*z + z
      end,
      onGhosts = true,
   }
   initSrc:advance(0.,{},{src})

   local initCheck = Updater.EvalOnNodes {
      onGrid   = grid,
      basis    = basis,
      evaluate = function (t,xn)
         local x = xn[1]
         local y = xn[2]
         local z = xn[3]
         local ky1 = 2*math.pi/Ly
         local ky3 = 6*math.pi/Ly
         return 0.5*math.cos(ky1*y) + 0.3*math.cos(ky3*y)*x*z
      end
   }
   initCheck:advance(0.,{},{check})

   local filter = Updater.FemKyFourierFilter {
      onGrid = grid,
      basis = basis,
      kyfilter = {2*math.pi/Ly, 6*math.pi/Ly},
   }
   filter:advance(0., {src}, {sol})

   local err = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {1, 1}
   }

   --src:write("src3d.bp")
   --sol:write("sol3d.bp")
   --check:write("check3d.bp")

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
   assert_close(intErr, 0., 1e-14, "checking integrated RMS error in filtered solution (3d)")
   return intErr
end   
   
test_2d(32, 32)
test_3d(32, 32, 32)

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
