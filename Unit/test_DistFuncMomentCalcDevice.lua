-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to compute moments on a device.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Don't do anything if we were not built with CUDA.
if GKYL_HAVE_CUDA == false then
   print("**** Can't run CUDA tests without CUDA enabled GPUs!")
   return 0
end

local Unit        = require "Unit"
local Grid        = require "Grid"
local DataStruct  = require "DataStruct"
local Basis       = require "Basis"
local Updater     = require "Updater"
local Lin         = require "Lib.Linalg"
local cudaRunTime = require "Cuda.RunTime"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

function test_ser_1x1v()
   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0, 6.0},
      cells = {8, 16},
   }
   local confGrid = Grid.RectCart {
      lower = { phaseGrid:lower(1) },
      upper = { phaseGrid:upper(1) },
      cells = { phaseGrid:numCells(1) },
   }
   -- Basis functions.
   local phaseBasis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 1 }
   local confBasis  = Basis.CartModalSerendipity { ndim = 1, polyOrder = 1 }
   -- Fields.
   local distF = DataStruct.Field {
      onGrid           = phaseGrid,
      numComponents    = phaseBasis:numBasis(),
      ghost            = {0, 0},
      createDeviceCopy = true,
   }
   local numDensity = DataStruct.Field {
      onGrid           = confGrid,
      numComponents    = confBasis:numBasis(),
      ghost            = {0, 0},
      createDeviceCopy = true,
   }

   -- Updater to initialize distribution function.
   local project = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,
      basis    = phaseBasis,
      evaluate = function (t, xn)
         x, v = xn[1], xn[2]
         n0   = 1.0
         u    = 0.0
         vt   = 0.8
         k    = 2.0
	 return (n0*math.cos(2.0*math.pi*k*x))/math.sqrt(2.0*math.pi*vt^2)
               *math.exp( -((v-u)^2)/(2.0*(vt^2)) )
      end
   }
   project:advance(0.0, {}, {distF})

   -- Copy distribution function to device.
   local err = distF:copyHostToDevice()

   -- Moment updater.
   local calcNumDensity = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "M0",
      onDevice   = true,
   }
   calcNumDensity:advance(0.0, {distF}, {numDensity})

   err = cudaRunTime.DeviceSynchronize()

   -- Copy moment from device to host.
   numDensity:copyDeviceToHost()

   numDensity:write("numDensity.bp",0.0)
--   local momIdxr = numDensity:genIndexer()
--   local nItr = numDensity:get(momIdxr( {1} ))
--   assert_equal(1, nItr[1]/math.sqrt(2), "Checking moment")

   err = cudaRunTime.DeviceReset()
   
end

test_ser_1x1v()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
