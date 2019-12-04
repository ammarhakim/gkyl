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
local cuAlloc = require "Cuda.Alloc"
local Range = require "Lib.Range"
local Grid = require "Grid"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

ffi.cdef [[
  void unit_sumArray(int numBlocks, int numThreads, int n, double a, double *x, double *y);
  void unit_sayHello();
  void unit_showRange(Range_t *range);

    typedef struct {
        int32_t ndim;
        int32_t cells[6];
        double lower[6], upper[6];
        double vol, dx[6];
    } GkylRectCart;

  void unit_showGrid(GkylRectCart *grid);
]]

-- basis tests
function test_1()
   assert_equal(GKYL_CUDA_DRIVER_VERSION, cuda.DriverGetVersion(), "Checking CUDA driver version")
   local devNum = cuda.GetDevice()
   assert_equal(0, devNum, "Checking device number")

   local prop, err = cuda.GetDeviceProperties(devNum)
   assert_equal(cuda.Success, err, "Checking of GetDeviceProperties worked")

   ffi.C.unit_sayHello()
   
   local range = Range.Range({0, 0, 0}, {1, 5, 20})
   local cuRange = Range.copyHostToDevice(range)
   ffi.C.unit_showRange(cuRange)

   cuda.Free(cuRange)

   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20}
   }
   local cuGrid = grid:copyHostToDevice()
   ffi.C.unit_showGrid(cuGrid)
end

-- raw mem management and kernel launches
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

-- managed mem management and kernel launches
function test_3()
   local devNum = cuda.GetDevice()
   local prop, err = cuda.GetDeviceProperties(devNum)
   -- return if this device does not support managed memory
   if prop.managedMemory ~= 1 then
      return 
   end

   local dtype = ffi.typeof("struct { double *_data; }")
   local m_x, m_y = ffi.new(dtype), ffi.new(dtype)

   local len = 1e6 
   -- allocate managed memory
   local elmSz = ffi.sizeof("double")
   local err
   m_x._data, err = cuda.MallocManaged(elmSz*len)
   assert_equal(cuda.Success, err, "Checking if MallocManaged worked")
   m_y._data, err = cuda.MallocManaged(elmSz*len)
   assert_equal(cuda.Success, err, "Checking if MallocManaged worked")

   for i = 1, len do
      m_x._data[i-1], m_y._data[i-1] = 10.5*i, 0.0
   end

   -- call kernel to do sum
   local numThread = 256
   local numBlock = math.floor(len/numThread)+1
   ffi.C.unit_sumArray(numBlock, numThread, len, 2.5, m_x._data, m_y._data)

   err = cuda.DeviceSynchronize()

   -- check if kernel worked
   local pass = true
   for i = 1, len do
      if m_y._data[i-1] ~= 2.5*m_x._data[i-1] then
	 pass = false
	 break
      end
   end
   assert_equal(true, pass, "Checking if sumArray kernel worked")

   cuda.Free(m_x._data); cuda.Free(m_y._data)
end

-- mem management via Cuda.Alloc and kernel launches
function test_4()
   local len = 1e6

   -- allocate memory on host
   local h_x, h_y = Alloc.Double(len), Alloc.Double(len)
   for i = 1, len do
      h_x[i], h_y[i] = 10.5*i, 0.0
   end

   -- check if kernel worked (this is a null test to ensure that the kernel
   -- call below is not "passing" erroneously
   local pass = true
   for i = 1, len do
      if h_y[i] ~= 2.5*h_x[i] then
	 pass = false
	 break
      end
   end
   assert_equal(false, pass, "Checking if sumArray kernel worked")

   -- allocate memory on device
   local d_x, d_y = cuAlloc.Double(len), cuAlloc.Double(len)
   assert_equal(len, d_x:size(), "Checking size of cuAlloc")
   assert_equal(len, d_y:size(), "Checking size of cuAlloc")

   -- Copy from host memory
   local err = d_x:copyHostToDevice(h_x)
   assert_equal(cuda.Success, err, "Checking if Memcpy worked")
   local err = d_y:copyHostToDevice(h_y)
   assert_equal(cuda.Success, err, "Checking if Memcpy worked")

   -- call kernel to do sum
   local numThread = 256
   local numBlock = math.floor(len/numThread)+1
   ffi.C.unit_sumArray(numBlock, numThread, len, 2.5, d_x:data(), d_y:data())

   -- Copy to host memory
   local err = d_x:copyDeviceToHost(h_x)
   assert_equal(cuda.Success, err, "Checking if Memcpy worked")
   local err = d_y:copyDeviceToHost(h_y)
   assert_equal(cuda.Success, err, "Checking if Memcpy worked")

   -- check if kernel worked
   local pass = true
   for i = 1, len do
      if h_y[i] ~= 2.5*h_x[i] then
   	 pass = false
   	 break
      end
   end
   assert_equal(true, pass, "Checking if sumArray kernel worked")
end

-- managed mem management via Cuda.Alloc and kernel launches
function test_5()
   local len = 1e6

   -- allocate managed memory on host
   local d_x, d_y = cuAlloc.ManagedDouble(len, true), cuAlloc.ManagedDouble(len, true)
   for i = 1, len do
      d_x[i], d_y[i] = 10.5*i, 0.0
   end

   -- call kernel to do sum
   local numThread = 256
   local numBlock = math.floor(len/numThread)+1
   ffi.C.unit_sumArray(numBlock, numThread, len, 2.5, d_x:data(), d_y:data())

   err = cuda.DeviceSynchronize()

   -- check if kernel worked
   local pass = true
   for i = 1, len do
      if d_y[i] ~= 2.5*d_x[i] then
   	 pass = false
   	 break
      end
   end
   assert_equal(true, pass, "Checking if sumArray kernel worked")
end

-- Run tests
test_1()
test_2()
test_3()
test_4()
test_5()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
