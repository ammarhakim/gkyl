-- Gkyl ------------------------------------------------------------------------
--
-- Test for HyperDisCont updater
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Don't do anything if we were not built with CUDA.
if GKYL_HAVE_CUDA == false then
   print("**** Can't run CUDA tests without CUDA enabled GPUs!")
   return 0
end

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local ffi  = require "ffi"
local Unit = require "Unit"
local Maxwell = require "Eq.PerfMaxwell"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"
local xsys = require "xsys"
local Alloc      = require "Lib.Alloc"
local cuda = nil
if GKYL_HAVE_CUDA then
  cuda = require "Cuda.RunTime"
  cuAlloc = require "Cuda.Alloc"
end

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats = Unit.stats

function test_1()
   local nloop = NLOOP or 1 -- number of HyperDisCont calls to loop over
   local runCPU = xsys.pickBool(RUNCPU, true)
   local checkResult = runCPU and true -- whether to check device result with host one, element-by-element. this can be expensive for large domains.
   local numThreads = NTHREADS or 128 -- number of threads to use in HyperDisCont kernel configuration
   local useSharedMemory = xsys.pickBool(SHARED, false) -- whether to use device shared memory

   -- set up dimensionality and basis parameters. 
   -- these parameters needs to match what is hard-coded in kernel template at bottom of Eq/GkylVlasov.h for now.
   local cdim = 2 -- number of configuration space dimensions
   local polyOrder = 2 -- polynomial order of basis (currently 1, 2, or 3 is supported for Vlasov on GPU)

   local confBasis = Basis.CartModalSerendipity { ndim = cdim, polyOrder = polyOrder }

   -- set up grids. adjust number of cells to increase domain size (more work for GPU).
   local nx = 512 -- number of configuration space dimensions in x
   local ny = 512 -- number of configuration space dimensions in y

   local confGrid = Grid.RectCart {
      cells = {nx, ny},
   }
   
   local maxwellEq = Maxwell {
      basis = confBasis,
      lightSpeed = 1.0,
   }

   local solver = Updater.HyperDisCont {
      onGrid = confGrid,
      basis = confBasis,
      cfl = 1.0,
      equation = maxwellEq,
      clearOut = true,
      noPenaltyFlux = true, -- penalty flux not yet implemented on device
      numThreads = numThreads,
      useSharedDevice = useSharedMemory,
   }

   local emField = DataStruct.Field {
      onGrid = confGrid,
      numComponents = 8*confBasis:numBasis(),
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   emField:clear(2)   
   local indexer = emField:genIndexer()
   for idx in emField:localRangeIter() do
      local emitr = emField:get(indexer(idx))
      emitr[1] = idx[1]+2*idx[2]+1
      emitr[2] = idx[1]+2*idx[2]+2
      emitr[3] = idx[1]+2*idx[2]+3
   end
   emField:copyHostToDevice()

   local emFieldRhs = DataStruct.Field {
      onGrid = confGrid,
      numComponents = 8*confBasis:numBasis(),
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   
   local d_emFieldRhs = DataStruct.Field {
      onGrid = confGrid,
      numComponents = 8*confBasis:numBasis(),
      ghost = {1, 1},
      createDeviceCopy = true,
   }

   local cflRateByCell = DataStruct.Field {
      onGrid = confGrid,
      numComponents = 1,
      ghost = {1, 1},
      createDeviceCopy = true,
   }

   solver:setDtAndCflRate(.1, cflRateByCell)

   print("Running GPU kernel")
   tmStart = Time.clock()
   for i = 1, nloop do
      solver:_advanceOnDevice(0.0, {emField}, {d_emFieldRhs})
   end
   -- Need to synchronize so that kernel actually runs!
   local err = cuda.DeviceSynchronize()
   local totalGpuTime = (Time.clock()-tmStart)
   print("... done.", totalGpuTime)

   assert_equal(0, err, "cuda error")
   d_emFieldRhs:copyDeviceToHost()
   local d_cflRate = cuAlloc.Double(1)
   cflRateByCell:deviceReduce('max', d_cflRate)
   local cflRate_from_gpu = Alloc.Double(1)
   local err = d_cflRate:copyDeviceToHost(cflRate_from_gpu)

   local tmStart
   tmStart = Time.clock()
   if runCPU then
      for i = 1, nloop do
         solver:_advance(0.0, {emField}, {emFieldRhs})
      end
   end
   local totalCpuTime = (Time.clock()-tmStart)
   local cflRate = cflRateByCell:reduce('max')[1]

   local indexer = emFieldRhs:genIndexer()
   local d_indexer = d_emFieldRhs:genIndexer()
   if checkResult then 
      for idx in emFieldRhs:localRangeIter() do
         local fitr = emFieldRhs:get(indexer(idx))
         local d_fitr = d_emFieldRhs:get(d_indexer(idx))
         for i = 1, emFieldRhs:numComponents() do
            assert_close(fitr[i], d_fitr[i], 1e-6, string.format("index %d, component %d is incorrect", indexer(idx), i))
         end
      end
      assert_equal(cflRate, cflRate_from_gpu[1], "Checking max cflRate")
   end


   print(string.format("Total CPU time for %d HyperDisCont calls = %f s   (average = %f s)", nloop, totalCpuTime, totalCpuTime/nloop))
   print(string.format("Total GPU time for %d HyperDisCont calls = %f s   (average = %f s)", nloop, totalGpuTime, totalGpuTime/nloop))
   print(string.format("GPU speed-up = %fx!", totalCpuTime/totalGpuTime))
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
