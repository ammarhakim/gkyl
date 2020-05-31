-- Gkyl ------------------------------------------------------------------------
--
-- Test for HyperDisCont updater
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local ffi  = require "ffi"
local Unit = require "Unit"
local Vlasov = require "Eq.Vlasov"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"
local cuda = nil
if GKYL_HAVE_CUDA then
  cuda = require "Cuda.RunTime"
end

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats = Unit.stats

function test_1()
   local confBasis = Basis.CartModalSerendipity { ndim = 1, polyOrder = 2 }
   local phaseBasis = Basis.CartModalSerendipity { ndim = 3, polyOrder = 2 }

   local grid = Grid.RectCart {
      lower = {0.0, -6.0, -6.0},
      upper = {1.0, 6.0, 6.0},
      cells = {16, 8, 8},
   }

   local confGrid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {16},
   }
   
   local vlasovEq = Vlasov {
      confBasis = confBasis,
      phaseBasis = phaseBasis,
      charge = -1.0,
      mass = 1.0,
   }

   local solver = Updater.HyperDisCont {
      onGrid = grid,
      basis = phaseBasis,
      cfl = 1.0,
      equation = vlasovEq,
      zeroFluxDirections = {2,3},
      clearOut = true,
      noPenaltyFlux = true, -- penalty flux not yet implemented on device
   }

   local distf = DataStruct.Field {
      onGrid = grid,
      numComponents = phaseBasis:numBasis(),
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   distf:clear(1)
   distf:copyHostToDevice()

   local fRhs = DataStruct.Field {
      onGrid = grid,
      numComponents = phaseBasis:numBasis(),
      ghost = {1, 1},
      createDeviceCopy = true,
   }

   local d_fRhs = DataStruct.Field {
      onGrid = grid,
      numComponents = phaseBasis:numBasis(),
      ghost = {1, 1},
      createDeviceCopy = true,
   }

   local cflRateByCell = DataStruct.Field {
      onGrid = grid,
      numComponents = phaseBasis:numBasis(),
      ghost = {1, 1},
      createDeviceCopy = true,
   }

   local emField = DataStruct.Field {
      onGrid = confGrid,
      numComponents = 8*confBasis:numBasis(),
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   emField:clear(2)
   emField:copyHostToDevice()

   solver:setDtAndCflRate(.1, cflRateByCell)

   local N = 1000
   local tmStart
   tmStart = Time.clock()
   for i = 1, N do
      solver:advance(0.0, {distf, emField}, {fRhs})
   end
   local totalCpuTime = (Time.clock()-tmStart)

   tmStart = Time.clock()
   for i = 1, N do
      solver:_advanceOnDevice(0.0, {distf, emField}, {d_fRhs})
   end
   -- Need to synchronize so that kernel actually runs!
   local err = cuda.DeviceSynchronize()
   local totalGpuTime = (Time.clock()-tmStart)
   assert_equal(0, err, "cuda error")

   d_fRhs:copyDeviceToHost()
   
   local indexer = fRhs:genIndexer()
   local d_indexer = d_fRhs:genIndexer()
   for idx in fRhs:localRangeIter() do
      local fitr = fRhs:get(indexer(idx))
      local d_fitr = d_fRhs:get(d_indexer(idx))
      for i = 0, fRhs:numComponents()-1 do
         assert_equal(fitr[i], d_fitr[i], string.format("index %d, component %d is incorrect", indexer(idx), i))
      end
   end

   print(string.format("Total CPU time for %d HyperDisCont calls = %f s", N, totalCpuTime))
   print(string.format("Total GPU time for %d HyperDisCont calls = %f s", N, totalGpuTime))
   print(string.format("GPU speed-up = %fx!", totalCpuTime/totalGpuTime))
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
