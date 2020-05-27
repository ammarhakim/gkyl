-- Gkyl ------------------------------------------------------------------------
--
-- Test of a simple updater that loops through the grid and does some operation.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Don't do anything if we were not built with CUDA.
if GKYL_HAVE_CUDA == false then
   print("**** Can't run CUDA tests without CUDA enabled GPUs!")
   return 0
end

local Unit        = require "Unit"
local Grid        = require "Grid"
local DataStruct  = require "DataStruct"
local Basis       = require "Basis"
local Updater     = require "Updater"
local Lin         = require "Lib.Linalg"
local cudaRunTime = require "Cuda.RunTime"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

function test_1()
   local grid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0, 6.0},
      cells = {8, 16},
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
   local result = DataStruct.Field {
      onGrid = grid,
      numComponents = 3,
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   
   -- initial fields with some data (on host)
   local localRange = field:localRange()
   local indexer = field:indexer()
   for i = localRange:lower(1), localRange:upper(1) do
      local fitr = field:get(indexer(i))
      local f2itr = field2:get(indexer(i))
      fitr[1] = i+1
      fitr[2] = i+2
      fitr[3] = i+3
   
      f2itr[1] = 4*i+1
      f2itr[2] = 5*i+2
      f2itr[3] = 6*i+3
   end
   
   -- copy stuff to device
   local err = field:copyHostToDevice()
   local err2 = field2:copyHostToDevice()

   local range = field:localRange()
   
   local grid_d = grid:copyHostToDevice()
   local range_d = Range.copyHostToDevice(range)

   local numCellsLocal = range:volume()

   local numThreads = GKYL_DEFAULT_NUM_THREADS
   local numBlocks = math.floor(numCellsLocal/numThreads)+1

   ffi.C.cuda_looper(grid_d, range_d, numBlocks, numThreads, field:deviceDataPointer(), field2:deviceDataPointer(), result:deviceDataPointer())

   cudaRunTime.Free(range_d)
end
