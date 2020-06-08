-- Gkyl ------------------------------------------------------------------------
--
-- Test for fields on cartesian grids on CUDA device
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
--
-- Don't do anything if we were not built with CUDA.
if GKYL_HAVE_CUDA == false then
   print("**** Can't run CUDA tests without CUDA enabled GPUs!")
   return 0
end

local ffi = require "ffi"
local Unit = require "Unit"
local Grid = require "Grid"
local DataStruct = require "DataStruct"
local cuda = require "Cuda.RunTime"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

ffi.cdef [[
  void unit_showFieldRange(GkylCartField_t *f, double *g);
  void unit_showFieldGrid(GkylCartField_t *f);
  void unit_readAndWrite(int numBlocks, int numThreads, GkylCartField_t *f, GkylCartField_t *res);
  void unit_readAndWrite_shared(int numBlocks, int numThreads, int sharedSize, GkylCartField_t *f, GkylCartField_t *res);
]]

function test_1()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {1000},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   local field2 = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 1},
      createDeviceCopy = true,
   }

   local localRange = field:localRange()
   local indexer = field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      fitr[1] = i+1
      fitr[2] = i+2
      fitr[3] = i+3
   end

   -- copy stuff to device
   local err = field:copyHostToDevice()
   local err2 = field2:copyHostToDevice()
   assert_equal(0, err, "Checking if copy to device worked")
   assert_equal(0, err2, "Checking if copy to device worked")

   field:deviceScale(-0.5)
   field:deviceAbs()
   field2:deviceClear(1.0)
   field:deviceAccumulate(2.0, field2, 0.5, field2)
   field:copyDeviceToHost()

   if GKYL_HAVE_CUDA then
      for i = localRange:lower(1), localRange:upper(1) do
	 local fitr = field:get(indexer(i))
	 assert_equal(0.5*(i+1)+2.5, fitr[1], "Checking deviceScale")
	 assert_equal(0.5*(i+2)+2.5, fitr[2], "Checking deviceScale")
	 assert_equal(0.5*(i+3)+2.5, fitr[3], "Checking deviceScale")
      end
   end
end

function test_2()
   local grid = Grid.RectCart {
      lower = {0.0},
      upper = {1.0},
      cells = {20},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   field:clear(1.0)
   field:copyHostToDevice()
   ffi.C.unit_showFieldRange(field._onDevice, field._devAllocData:data())
   ffi.C.unit_showFieldGrid(field._onDevice)
end

-- read a CartField, and then write it to a result CartField (basically a copy via a kernel)
function test_3()
   local grid = Grid.RectCart {
      cells = {4, 32},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 32,
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   field:clear(-1)
   local indexer = field:genIndexer()
   for idx in field:localExtRangeIter() do
      local fitr = field:get(indexer(idx))
      fitr[1] = 5*idx[1] + 2*idx[2]+1
      fitr[2] = 5*idx[1] + 2*idx[2]+2
      fitr[3] = 5*idx[1] + 2*idx[2]+3
   end
   field:copyHostToDevice()
   local result = DataStruct.Field {
      onGrid = grid,
      numComponents = 32,
      ghost = {1, 1},
      createDeviceCopy = true,
   }

   local numCellsLocal = field:localRange():volume()
   local numThreads = 32
   local numBlocks  = math.ceil(numCellsLocal/numThreads)
   ffi.C.unit_readAndWrite(numBlocks, numThreads, field._onDevice, result._onDevice)
   result:copyDeviceToHost()

   for idx in field:localRangeIter() do
      local fitr = field:get(indexer(idx))
      local ritr = result:get(indexer(idx))
      for i = 1, field:numComponents() do
         assert_equal(fitr[i], ritr[i], string.format("readAndWrite test: incorrect element at index %d, component %d", indexer(idx), i))
      end
   end
end

function test_4()
   local grid = Grid.RectCart {
      cells = {8, 8},
   }
   local field = DataStruct.Field {
      onGrid = grid,
      numComponents = 32,
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   field:clear(-1)
   local indexer = field:genIndexer()
   for idx in field:localExtRangeIter() do
      local fitr = field:get(indexer(idx))
      fitr[1] = 5*idx[1] + 2*idx[2]+1
      fitr[2] = 5*idx[1] + 2*idx[2]+2
      fitr[3] = 5*idx[1] + 2*idx[2]+3
   end
   field:copyHostToDevice()
   local result = DataStruct.Field {
      onGrid = grid,
      numComponents = 32,
      ghost = {1, 1},
      createDeviceCopy = true,
   }

   local numCellsLocal = field:localRange():volume()
   local numThreads = 32
   local numBlocks  = math.ceil(numCellsLocal/numThreads)
   local sharedSize = numThreads*field:numComponents()
   ffi.C.unit_readAndWrite_shared(numBlocks, numThreads, sharedSize, field._onDevice, result._onDevice)
   local err = cuda.DeviceSynchronize()
   assert_equal(0, err, "cuda error")
   result:copyDeviceToHost()

   for idx in field:localRangeIter() do
      local fitr = field:get(indexer(idx))
      local ritr = result:get(indexer(idx))
      for i = 1, field:numComponents() do
         assert_equal(fitr[i], ritr[i], string.format("readAndWrite_shared test: incorrect element at index %d, component %d", indexer(idx)-1, i))
      end
   end
end

test_1()
test_2()
test_3()
test_4()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
