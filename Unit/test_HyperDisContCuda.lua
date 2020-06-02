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
   local nloop = 1 -- number of HyperDisCont calls to loop over
   local checkResult = false -- whether to check device result with host one, element-by-element. this can be expensive for large domains.
   local numThreads = 256 -- number of threads to use in HyperDisCont kernel configuration

   -- set up dimensionality and basis parameters. 
   -- these parameters needs to match what is hard-coded in kernel template at bottom of Eq/GkylVlasov.h for now.
   local cdim = 1 -- number of configuration space dimensions
   local vdim = 2 -- number of velocity space dimensions
   local polyOrder = 2 -- polynomial order of basis (currently 1, 2, or 3 is supported for Vlasov on GPU)

   local pdim = cdim + vdim -- total number of dimensions in phase space
   local confBasis = Basis.CartModalSerendipity { ndim = cdim, polyOrder = polyOrder }
   local phaseBasis = Basis.CartModalSerendipity { ndim = pdim, polyOrder = polyOrder }

   -- set up grids. adjust number of cells to increase domain size (more work for GPU).
   local nx = 64 -- number of configuration space dimensions in x
   local nvx = 128  -- number of velocity dimensions in vx
   local nvy = 128  -- number of velocity dimensions in vy 

   local grid = Grid.RectCart {
      cells = {nx, nvx, nvy},
   }
   local confGrid = Grid.RectCart {
      cells = {nx},
   }
   
   local vlasovEq = Vlasov {
      confBasis = confBasis,
      phaseBasis = phaseBasis,
      charge = -1.0,
      mass = 1.0,
   }

   local zfd = { }
   for d = 1, vdim do
     table.insert(zfd, cdim+d)
   end
   local solver = Updater.HyperDisCont {
      onGrid = grid,
      basis = phaseBasis,
      cfl = 1.0,
      equation = vlasovEq,
      zeroFluxDirections = zfd,
      clearOut = true,
      noPenaltyFlux = true, -- penalty flux not yet implemented on device
      numThreads = numThreads,
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

   local tmStart
   tmStart = Time.clock()
   for i = 1, nloop do
      solver:advance(0.0, {distf, emField}, {fRhs})
   end
   local totalCpuTime = (Time.clock()-tmStart)

   tmStart = Time.clock()
   for i = 1, nloop do
      solver:_advanceOnDevice(0.0, {distf, emField}, {d_fRhs})
   end
   -- Need to synchronize so that kernel actually runs!
   local err = cuda.DeviceSynchronize()
   local totalGpuTime = (Time.clock()-tmStart)
   assert_equal(0, err, "cuda error")

   d_fRhs:copyDeviceToHost()
   
   local indexer = fRhs:genIndexer()
   local d_indexer = d_fRhs:genIndexer()
   if checkResult then 
      for idx in fRhs:localRangeIter() do
         local fitr = fRhs:get(indexer(idx))
         local d_fitr = d_fRhs:get(d_indexer(idx))
         for i = 0, fRhs:numComponents()-1 do
            assert_close(fitr[i], d_fitr[i], 1e-11, string.format("index %d, component %d is incorrect", indexer(idx), i))
         end
      end
   end

   print(string.format("Total CPU time for %d HyperDisCont calls = %f s", nloop, totalCpuTime))
   print(string.format("Total GPU time for %d HyperDisCont calls = %f s", nloop, totalGpuTime))
   print(string.format("GPU speed-up = %fx!", totalCpuTime/totalGpuTime))
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
