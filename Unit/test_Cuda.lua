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

function test_1()
    assert_equal(GKYL_CUDA_DRIVER_VERSION, cuda.DriverGetVersion(), "Checking CUDA driver version")
end

function test_2()
   assert_equal(0, cuda.MemcpyHostToHost, "Checking MemcpyKind")
   assert_equal(1, cuda.MemcpyHostToDevice, "Checking MemcpyKind")
   assert_equal(2, cuda.MemcpyDeviceToHost, "Checking MemcpyKind")
   assert_equal(3, cuda.MemcpyDeviceToDevice, "Checking MemcpyKind")
   assert_equal(4, cuda.MemcpyDefault, "Checking MemcpyKind")

   local len = 100

   local hostMem = Alloc.Double(len)
   for i = 1, hostMem:size() do
       hostMem[i] = 10.5*i
   end

   -- Allocate memory
   local elmSz = ffi.sizeof("double")
   local devMem = cuda.Malloc(elmSz*len)

   -- Copy host -> device memory
   local err = cuda.Memcpy(devMem, hostMem:data(), hostMem:size()*elmSz, cuda.MemcpyHostToDevice);
   assert_equal(cuda.Success, err, "Checking if Memcpy worked")

   cuda.Free(d)
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
