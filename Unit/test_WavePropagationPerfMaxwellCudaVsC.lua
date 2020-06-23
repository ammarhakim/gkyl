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
local PerfMaxwell = require "Eq.PerfMaxwell"
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
   local nloop = NLOOP or 1
   local runCPU = xsys.pickBool(RUNCPU, true)
   local check = xsys.pickBool(CHECK, runCPU and truex) -- check & report bad values
   local count = xsys.pickBool(COUNT, true) -- count bad/good values
   local useGlobalMemory = xsys.pickBool(GLOBALMEM, true)
   local numThreads = NTHREADS or 128
   local nx = NX or 128*32

   local grid = Grid.RectCart {
      cells = {nx},
      periodicDirs = {1},
   }
   
   local eq = PerfMaxwell {
      lightSpeed = 5.0,
      elcErrorSpeedFactor = 1.0,
      mgnErrorSpeedFactor = 3.0,
   }
   local numEquations = 8

   local solver = Updater.WavePropagation {
      onGrid = grid,
      updateDirections = {1},
      cfl = 1.0,
      equation = eq,
      limiter = 'min-mod',
      numThreads = numThreads,
      useGlobalMemory = useGlobalMemory,
   }

   local qIn = DataStruct.Field {
      onGrid = grid,
      numComponents = numEquations,
      ghost = {2, 2},
   }
   qIn:clear(1)
   local indexer = qIn:genIndexer()
   for idx in qIn:localRangeIter() do
      local fitr = qIn:get(indexer(idx))
      local x = idx[1]
      fitr[1] = 1.8 * (1 + 0.11 * math.sin( x/nx * 2 * math.pi ))
      fitr[2] = 2.7 * (1 + 0.12 * math.cos( x/nx * 2 * math.pi ))
      fitr[3] = 3.6 * (1 + 0.13 * math.sin( x/nx * 2 * math.pi ))
      fitr[4] = 4.5 * (1 + 0.14 * math.cos( x/nx * 2 * math.pi ))
      fitr[5] = 5.4 * (1 + 0.15 * math.sin( x/nx * 2 * math.pi ))
      fitr[6] = 6.3 * (1 + 0.16 * math.cos( x/nx * 2 * math.pi ))
      fitr[7] = 7.2 * (1 + 0.17 * math.sin( x/nx * 2 * math.pi ))
      fitr[8] = 8.1 * (1 + 0.18 * math.cos( x/nx * 2 * math.pi ))
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

   solver:setDtAndCflRate(1e-5, nil)

   tmStart = Time.clock()
   for i = 1, nloop do
      solver:_advanceOnDevice(0.0, {qIn}, {d_qOut})
   end
   -- Need to synchronize so that kernel actually runs!
   local err = cuda.DeviceSynchronize()
   local totalGpuTime = (Time.clock()-tmStart)
   if err~=0 then
      print(string.format("_advanceOnDevice returned error code %d; not checking results", err))
      assert_equal(0, err, "cuda error")
      return
   end
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
   if check or count then 
      local good =0
      local bad = 0
      for idx in qOut:localRangeIter() do
         local fitr = qOut:get(indexer(idx))
         local d_fitr = d_qOut:get(d_indexer(idx))
         for i = 1, qOut:numComponents() do
            if count then
               if math.abs(fitr[i]-d_fitr[i])>1e-10 then
                  bad = bad+1
               else
                  good = good+1
               end
            end
            if check then
               assert_close(
                  fitr[i], d_fitr[i], 1e-10, 
                  string.format(
                  "index %d, component %d is incorrect", indexer(idx), i))
            end

         end
      end
      --assert_equal(cflRate, cflRate_from_gpu[1], "Checking max cflRate")
      if count then
         print("# good values:", good, "; # bad values", bad)
      end
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
