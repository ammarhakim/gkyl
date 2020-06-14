-- Gkyl ------------------------------------------------------------------------
--
-- Test for WavePropagation updater
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
local Euler = require "Eq.Euler"
local Updater = require "Updater"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"
local xsys = require "xsys"
local Alloc      = require "Lib.Alloc"
local cuda = nil
local cuda
local cuAlloc
if GKYL_HAVE_CUDA then
  cuda = require "Cuda.RunTime"
  cuAlloc = require "Cuda.Alloc"
end

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats = Unit.stats

function test_1()
   local nloop = NLOOP or 1 -- number of WavePropagation calls to loop over
   local runCPU = xsys.pickBool(RUNCPU, true)
   local checkResult = runCPU and true -- whether to check device result with host one, element-by-element. this can be expensive for large domains.
   local numThreads = NTHREADS or 8 -- number of threads to use in WavePropagation kernel configuration
   local useSharedDevice = xsys.pickBool(SHARED, false) -- whether to use device shared memory

   -- local nx = 1024 -- number of configuration space dimensions in x
   -- local nx = 128*64 -- number of configuration space dimensions in x
   local nx = 16 -- number of configuration space dimensions in x

   local grid = Grid.RectCart {
      cells = {nx},
      periodicDirs = {1},
   }
   
   local eq = Euler {
      numericalFlux = "lax",
      gasGamma = 1.4,
   }
   local numEquations = 5

   local solver = Updater.WavePropagation {
      onGrid = grid,
      updateDirections = {1},
      cfl = 1.0,
      equation = eq,
      numThreads = numThreads,
      useSharedDevice = useSharedDevice,
      limiter = 'zero',
   }

   local qIn = DataStruct.Field {
      onGrid = grid,
      numComponents = numEquations,
      ghost = {2, 2},
   }
   qIn:clear(1)
   local indexer = qIn:genIndexer()
   local gg = eq._gasGamma
   local rho, u, v, w, pr = 1.0, 2.1, 3.1, 4.1, 0.1
   local er = pr/(gg-1) + 0.5*rho*(u*u+v*v+w*w)
   for idx in qIn:localRangeIter() do
      local fitr = qIn:get(indexer(idx))
      local x = idx[1]
      fitr[1] = rho * (1 + 0.1 * math.sin( x/nx * 2 * math.pi ))
      fitr[2] = rho*u * (1 + 0.11 * math.sin( x/nx * 2 * math.pi ))
      fitr[3] = rho*v * (1 + 0.12 * math.cos( x/nx * 2 * math.pi ))
      fitr[4] = rho*w * (1 + 0.13 * math.cos( x/nx * 2 * math.pi ))
      fitr[5] = pr/(gg-1) + 0.5 * (fitr[2]^2 + fitr[3]^2 + fitr[4]^2) / fitr[1]
   end
   qIn:sync()
   qIn:copyHostToDevice()

   local qOut = DataStruct.Field {
      onGrid = grid,
      numComponents = numEquations,
      ghost = {2, 2},
      createDeviceCopy = true,
   }

   local d_qOut = DataStruct.Field {
      onGrid = grid,
      numComponents = numEquations,
      ghost = {2, 2},
      createDeviceCopy = true,
   }

   solver:setDtAndCflRate(.01, nil)

   tmStart = Time.clock()
   for i = 1, nloop do
      solver:_advanceOnDevice(0.0, {qIn}, {d_qOut})
   end
   -- Need to synchronize so that kernel actually runs!
   local err = cuda.DeviceSynchronize()
   local totalGpuTime = (Time.clock()-tmStart)
   assert_equal(0, err, "cuda error")
   d_qOut:copyDeviceToHost()
   d_qOut:sync()
   -- local cflRate_from_gpu = Alloc.Double(1)
   -- local err = solver.dtPtr:copyDeviceToHost(cflRate_from_gpu)

   local tmStart
   tmStart = Time.clock()
   if runCPU then
      for i = 1, nloop do
         solver:_advance(0.0, {qIn}, {qOut})
      end
   end
   local totalCpuTime = (Time.clock()-tmStart)
   qOut:sync()

   local indexer = qOut:genIndexer()
   local d_indexer = d_qOut:genIndexer()
   if checkResult then 
      for idx in qOut:localRangeIter() do
         local fitr = qOut:get(indexer(idx))
         local d_fitr = d_qOut:get(d_indexer(idx))
         for i = 1, qOut:numComponents() do
            assert_close(fitr[i], d_fitr[i], 1e-10, string.format("index %d, component %d is incorrect", indexer(idx), i))
         end
      end
      --assert_equal(cflRate, cflRate_from_gpu[1], "Checking max cflRate")
   end


   print(string.format("Total CPU time for %d WavePropagation calls = %f s   (average = %f s)", nloop, totalCpuTime, totalCpuTime/nloop))
   print(string.format("Total GPU time for %d WavePropagation calls = %f s   (average = %f s)", nloop, totalGpuTime, totalGpuTime/nloop))
   print(string.format("GPU speed-up = %fx!", totalCpuTime/totalGpuTime))
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
