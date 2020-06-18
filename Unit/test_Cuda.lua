-- Gkyl ------------------------------------------------------------------------
--
-- Test for CUDA runtime API wrappers.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Don't do anything if we were not built with CUDA.
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
  void unit_showRange(GkylRange_t *range);

  void unit_test_BasisTypes_1xp1_ser();
  void unit_test_BasisTypes_2xp2_ser();
  void unit_test_BasisTypes_5xp2_ser();

  void unit_showGrid(GkylRectCart_t *grid);
  void unit_getCellCenter(GkylRectCart_t *grid, int *idx, double *xc);

  typedef double (*sumFunc_t)(void *self, double a, double b);
  typedef struct {
      void *child;
      sumFunc_t sumFunc;
  } SimpleEquation_t;
  typedef struct {
      double gamma;
  } SimpleEulerEquation_t;

  double unit_test_SimpleEquation(SimpleEquation_t *eqn, double* ab);
  sumFunc_t getEulerSumFuncOnDevice();
]]

function test_0()
  assert_equal(0, cuda.Success, "Checking Success code")
  assert_equal(1, cuda.ErrorInvalidValue, "Checking ErrorInvalidValue code")
  assert_equal(2, cuda.ErrorMemoryAllocation, "Checking ErrorMemoryAllocation code")
  assert_equal(13, cuda.ErrorInvalidSymbol, "Checking ErrorInvalidSymbol code")
  assert_equal(21, cuda.ErrorInvalidMemcpyDirection, "Checking ErrorInvalidMemcpyDirection code")
end

-- basic tests
function test_1()
   assert_equal(GKYL_CUDA_DRIVER_VERSION, cuda.DriverGetVersion(), "Checking CUDA driver version")
   local devNum, err = cuda.GetDevice()
   assert_equal(0, devNum, "Checking device number")

   local err = cuda.SetDevice(devNum)
   assert_equal(cuda.Success, err, "Setting device")

   local prop, err = cuda.GetDeviceProperties(devNum)
   assert_equal(cuda.Success, err, "Checking of GetDeviceProperties worked")

   ffi.C.unit_sayHello()
   
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
   local devNum, _ = cuda.GetDevice()
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
   local d_x, d_y = cuAlloc.ManagedDouble(len), cuAlloc.ManagedDouble(len)
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

-- range tests
function test_6()
   local range = Range.Range({5}, {20})
   local cuRange = Range.copyHostToDevice(range)
   --ffi.C.unit_showRange(cuRange)

   cuda.Free(cuRange)
end

-- grid tests
function test_7()
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20}
   }
   local cuGrid = grid:copyHostToDevice()
   --ffi.C.unit_showGrid(cuGrid)

   local h_idx = Alloc.Int(2)
   h_idx[1] = 1
   h_idx[2] = 2
   local d_idx = cuAlloc.Int(2)
   d_idx:copyHostToDevice(h_idx)
 
   local d_xc = cuAlloc.Double(2)
   ffi.C.unit_getCellCenter(cuGrid, d_idx:data(), d_xc:data())

   local h_xc = Alloc.Double(2)
   d_xc:copyDeviceToHost(h_xc)

   assert_equal(h_xc[1]==0.1 and h_xc[2]==1.3, true, "Checking cell centers")
end

function test_8()
   ffi.C.unit_test_BasisTypes_1xp1_ser()
   ffi.C.unit_test_BasisTypes_2xp2_ser()
   ffi.C.unit_test_BasisTypes_5xp2_ser()
end

function test_9()
   -- Create Euler equation object on host
   local eulerEqn = ffi.new("SimpleEulerEquation_t")
   eulerEqn.gamma = 1.4
   -- copy EulerEqn to device
   local sz = ffi.sizeof("SimpleEulerEquation_t")
   local d_eulerEqn = cuda.Malloc(sz)
   local err = cuda.Memcpy(d_eulerEqn, eulerEqn, sz, cuda.MemcpyHostToDevice)

   -- Create Equation object on host
   local eqn = ffi.new("SimpleEquation_t")
   eqn.child = d_eulerEqn
   eqn.sumFunc = ffi.C.getEulerSumFuncOnDevice()

   -- copy equation object to device
   sz = ffi.sizeof("SimpleEquation_t")
   local d_eqn = cuda.Malloc(sz)
   err = cuda.Memcpy(d_eqn, eqn, sz, cuda.MemcpyHostToDevice)

   -- allocate and set a,b values
   local h_ab = Alloc.Double(2)
   h_ab[1], h_ab[2] = 1.0, 2.0

   -- allocate and copy to host
   local d_ab = cuAlloc.Double(2)
   local err = d_ab:copyHostToDevice(h_ab)
   
   local res = ffi.C.unit_test_SimpleEquation(d_eqn, d_ab:data())
   err = cuda.DeviceSynchronize()

   assert_equal(eulerEqn.gamma*(h_ab[1]+h_ab[2]), res, "Checking function pointer abstraction")
end

-- Run tests
test_0()
test_1()
test_2()
test_3()
test_4()
test_5()
test_6()
test_7()
test_8()
test_9()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
