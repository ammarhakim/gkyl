-- Gkyl ------------------------------------------------------------------------
--
-- Test for CUDA runtime API wrappers
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- don't do anything if we were not built with CUDA
if GKYL_HAVE_CUDA == false then
   print("**** Can't run CUDA tests without CUDA enabled GPUs!")
   return 0
end

local Unit = require "Unit"
local Alloc = require "Lib.Alloc"
local ffi = require "ffi"
local cuda = require "Cuda.RunTime"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

ffi.cdef [[
  void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y);
]]

function test_1()
    assert_equal(GKYL_CUDA_DRIVER_VERSION, cuda.DriverGetVersion(), "Checking CUDA driver version")
end

function test_2()
   assert_equal(0, cuda.MemcpyHostToHost, "Checking MemcpyKind")
   assert_equal(1, cuda.MemcpyHostToDevice, "Checking MemcpyKind")
   assert_equal(2, cuda.MemcpyDeviceToHost, "Checking MemcpyKind")
   assert_equal(3, cuda.MemcpyDeviceToDevice, "Checking MemcpyKind")
   assert_equal(4, cuda.MemcpyDefault, "Checking MemcpyKind")

   local len = 1e6

   -- allocate memory on host
   local h_x, h_y = Alloc.Double(len), Alloc.Double(len)
   for i = 1, len do
       h_x[i], h_y[i] = 10.5*i, 0.0
   end

   -- check if kernel worked (this is a null test to ensure that the kernel
   -- call below is not "passing" erroneously
   local pass = true
   for i = 1, h_y:size() do
      if h_y[i] ~= 2.5*h_x[i] then
        pass = false
        break
      end
   end
   assert_equal(false, pass, "Checking if sumArray kernel worked")

   -- allocate memory on device
   local elmSz = ffi.sizeof("double")
   local d_x, d_y = cuda.Malloc(elmSz*len), cuda.Malloc(elmSz*len)

   -- Copy host -> device memory
   local err = cuda.Memcpy(d_x, h_x:data(), h_x:size()*elmSz, cuda.MemcpyHostToDevice)
   assert_equal(cuda.Success, err, "Checking if Memcpy worked")
   local err = cuda.Memcpy(d_y, h_y:data(), h_y:size()*elmSz, cuda.MemcpyHostToDevice)
   assert_equal(cuda.Success, err, "Checking if Memcpy worked")

   -- call kernel to do sum
   local numThread = 256
   local numBlock = math.floor(len/numThread)+1
   ffi.C.unit_sumArray(numBlock, numThread, len, 2.5, d_x, d_y)

   -- Copy device -> host memory
   local err = cuda.Memcpy(h_x:data(), d_x, h_x:size()*elmSz, cuda.MemcpyDeviceToHost)
   assert_equal(cuda.Success, err, "Checking if Memcpy worked")
   local err = cuda.Memcpy(h_y:data(), d_y, h_y:size()*elmSz, cuda.MemcpyDeviceToHost)
   assert_equal(cuda.Success, err, "Checking if Memcpy worked")

   -- check if kernel worked
   local pass = true
   for i = 1, h_y:size() do
      if h_y[i] ~= 2.5*h_x[i] then
        pass = false
        break
      end
   end
   assert_equal(true, pass, "Checking if sumArray kernel worked")

   cuda.Free(d_x); cuda.Free(d_y)
end

-- Run tests
test_1()
test_2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
